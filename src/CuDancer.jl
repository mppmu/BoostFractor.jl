
#
#   BoostFractor Beam Propagation Core + PaD addition (work-in-progress)
#
#   V: 2019-04-05 (core) + 2019-06-09 (add-on)
#
#   Stefan Knirck
#		add-on: Dominik Bergermann

export SetupBoundaries, SeedSetupBoundaries, dancer, dance_intro, propagator, propagator1D,init_coords



"""
Propagates the fields through the system

* `amin`:           Mimum (local) amplitude of a field, in order to be propagated
* `nmax`:           Maximum number of beam iteration steps, directly equivalent to how many boundaries a beam 'sees' at most
* `bdry`:           SetupBoundaries-Objekt, containing all relevant geometrical information (disk positions, epsilon, etc)
* `f`:              Frequency
* `prop`:           Propagator Function to use. Standard is propagator()
* `reflect`:        If nothing (standar value), the axion-induced signal is computed.
                    If set, this field defines a beam, for which the reflected beam will be calculated
* `Xset`, `Yset`:   Explicitly set the coordinate system for the fields
* `returnsum`:      If false, the out-propagating contributions after each iteration will be returned, without summing.
* `immediatesum`:   If false, the out-propagating contributions will be saved and summed up at the end.
"""
function cudancer(amin, nmax, bdry::SetupBoundaries; f=10.0e9, prop=propagator, emit=nothing, reflect=nothing, Xset=X, Yset=Y, diskR=0.1, returnsum=true, immediatesum=true)
    init_coords(Xset, Yset);

    rightmoving = 1
    leftmoving = 2

    fields = Array{Complex{Float64}}(zeros(length(bdry.distance), 2, length(X), length(Y)))
    #fields = SharedArray(Complex{Float64},  length(bdry.distance), 2, length(X), length(Y))
    #                       number of regions -----^                ^              ^
    #                       number of propagation directions -------^              ^
    #                       dimensions of the fields at each position -------------^

	#Pre-allocate memory on GPU
	if immediatesum
		Eout = CuArray{ComplexF64, 3}(undef, length(X), length(Y), 1)
	else
		Eout = CuArray{ComplexF64, 3}(undef, length(X), length(Y), nmax+1)
	end

    # TODO: propagation through and emission from last bdry to the right
    if reflect == nothing
        if emit == nothing
            fields = dance_intro(bdry,X,Y)
        else
            fields = emit
        end
        # Eout =
    else
        fields[length(bdry.distance), leftmoving, :, :] = reflect
        # Eout =
    end

    c = 299792458.
    lambda = c/f

    n = 0   # Iteration Count
    a = 1.0 # Maximum Field Amplitude in the iteration

	#Copy data to GPU
	fields = CuArray(fields)
	FAL = CuArray{ComplexF64, 3}(undef, length(bdry.distance) - 1, length(X), length(Y))
	tilt_x = CuArray(bdry.relative_tilt_x)
	tilt_y = CuArray(bdry.relative_tilt_y)
	dist = CuArray(bdry.distance)
	ref = CuArray(bdry.r)
	surfaces = CuArray(bdry.relative_surfaces)

	k0 = CuArray(2*pi/lambda*sqrt.(bdry.eps))

	#allocate memory for the phase shift arrays
	surface_prop = CuArray{ComplexF32, 3}(undef, size(bdry.relative_surfaces))
	surface_prop_out = CuArray{ComplexF64, 2}(undef, length(X), length(Y))
	k0_prop = CuArray{ComplexF32, 3}(undef, size(bdry.relative_surfaces))

	#calculate block configuration (Grid dimensions)
	block_xy = min(64, length(X)÷16)
	block_z = min(length(bdry.distance), 4096 ÷ block_xy^2)

	#generate surface shift arrays and perform an initail cutoff
	@cuda threads=(16,16) blocks=(block_xy, block_xy, block_z) gen_prop_cut(
				fields, k0_prop, surface_prop, surface_prop_out,
				k0, dist, diskR,
				tilt_x, tilt_y, surfaces,
				X, Y, coordsKx, coordsKy)

	#these were just used to generate the surface shifts and are not needed anymore
	CUDA.unsafe_free!(dist)
	CUDA.unsafe_free!(tilt_x)
	CUDA.unsafe_free!(tilt_y)
	CUDA.unsafe_free!(k0)
	CUDA.unsafe_free!(surfaces)

	#pre-plan FFT and inverse
	plan = plan_fft!(fields, (3,4))
	inv_plan = CUDA.CUFFT.plan_inv(plan)

    # The measure of music (Takt) .... da capo al fine
	#CUDA.@time begin
    while n <= nmax && a >= amin

		plan * fields #Execute FFT

		#apply momentum space phase shifts
		fields[:, 1, :, :] .*= k0_prop
		fields[:, 2, :, :] .*= k0_prop

		#iFFT
		inv_plan * fields #Execute iFFT

		#return Eout before surface and tilt phase shifts for the last region are applied
		if immediatesum
			Eout[:,:,1]    .+= (1+ref[length(bdry.r)]).*fields[length(bdry.r)-1, rightmoving,:,:]
		else
			Eout[:,:,n+1]    = (1+ref[length(bdry.r)]).*fields[length(bdry.r)-1, rightmoving,:,:]
		end

		#surface and tilt phase shifts, cutoff at the disk radius
		fields[:, 1, :, :] .*= surface_prop
		fields[:, 2, :, :] .*= surface_prop

		#Reflection and Transmission

		# Always the field on the left gets a copy
		FAL = fields[1:length(bdry.r)-2, rightmoving,:,:]

        # Then we in-place update it
        fields[1:length(bdry.r)-1, rightmoving,:,:] .*= ref[2:length(bdry.r)] #combined with reflection rightmost
        fields[1:length(bdry.r)-2, rightmoving,:,:] .+= (1 .- ref[2:length(bdry.r)-1]).*fields[2:length(bdry.r)-1, leftmoving,:,:]

        # Then we in-place update the other one with the previously made backup
        fields[1:length(bdry.r)-1, leftmoving,:,:] .*= -ref[1:length(bdry.r)-1] #combined with reflection leftmost
        fields[2:length(bdry.r)-1, leftmoving,:,:] .+= (1 .+ ref[2:length(bdry.r)-1]).*FAL

        # Now the directions have changed
        lmv = copy(leftmoving)
        leftmoving = rightmoving
        rightmoving = lmv

        # a is the maximum field that still occurs in the system
        if amin != 0
            a = maximum(abs.(fields))
        end
        # n is the iteration count
        n += 1
    end
	#end

	#sufrace and tilt phase shift without cutoff for the last region
	Eout .*= surface_prop_out

	#destroy plans
	finalize(plan)
	finalize(inv_plan)

	#free GPU Memory
	CUDA.unsafe_free!(fields)
	CUDA.unsafe_free!(surface_prop)
	CUDA.unsafe_free!(surface_prop_out)
	CUDA.unsafe_free!(k0_prop)
	CUDA.unsafe_free!(ref)
	CUDA.unsafe_free!(FAL)

	#move Eout back to CPU
	out = Array(Eout)
	CUDA.unsafe_free!(Eout)

	GC.gc(false)

    # Dance Outro
    if returnsum
        return sum(reverse(out, dims = 3), dims = 3)[:,:,1]
    else
        return out
    end
end

