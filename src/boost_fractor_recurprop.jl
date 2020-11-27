#
# Recursive Fourier Propagation for BoostFractor
#
# V: 2019-08-06
#
# Stefan Knirck
#

export dancer, dance_intro, cheerleader

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
function dancer(amin, nmax, sbdry::SetupBoundaries, coords::CoordinateSystem; f=10.0e9, prop=propagator, emit=nothing, reflect=nothing, diskR=0.1, returnsum=true, immediatesum=true)

    # Make dancer swallow the same SetupBoundaries object as cheerleader and transformer
    bdry = deepcopy(sbdry)
    bdry.eps = bdry.eps[2:end]
    bdry.distance = bdry.distance[2:end]
    bdry.relative_tilt_x = bdry.relative_tilt_x[2:end]
    bdry.relative_tilt_y = bdry.relative_tilt_y[2:end]
    bdry.relative_surfaces = bdry.relative_surfaces[2:end,:,:]
    append!(bdry.r, 0.0)

    rightmoving = 1
    leftmoving = 2

    fields = Array{Complex{Float64}}(zeros(length(coords.X), length(coords.Y), length(bdry.distance), 2))
    #fields = SharedArray(Complex{Float64}, length(coords.X), length(coords.Y), length(bdry.distance), 2)
    #           dimensions of the fields at each position ----^
    #                       number of regions ----------------------------------------^
    #                                  number of propagation directions ------------------------------^

    fields_arrived_left = Array{ComplexF64, 3}(undef, length(coords.X), length(coords.Y), length(bdry.distance) - 1)

    # Pre-allocate memory
    if immediatesum
        Eout = Array{Complex{Float64}}(zeros(length(coords.X), length(coords.Y),1))
    else
        Eout = Array{Complex{Float64}}(zeros(length(coords.X), length(coords.Y),nmax+1))
    end

    #create fft plans
    plan = plan_fft!(fields[:,:,1,rightmoving]) #plan_fft!(fields, (1, 2)) when transforming all boundaries at once
    i_plan = plan_ifft!(fields[:,:,1,rightmoving])

    # TODO: propagation through and emission from last bdry to the right
    if reflect === nothing
        if emit === nothing
           fields = dance_intro(bdry,coords;diskR=diskR)
        else
            fields = emit
        end
        # Eout =
    else
        fields[:, :, length(bdry.distance), leftmoving] = reflect
        # Eout =
    end

    lambda = wavelength(f)

    #generate phase shift arrays
    phase = init_prop(bdry, coords, lambda, diskR, output_surfaces = [length(bdry.distance)])

    n = 0   # Iteration Count
    a = 1.0 # Maximum Field Amplitude in the iteration

    # The measure of music (Takt) .... da capo al fine
    while n <= nmax && a >= amin
""" alternate version without propagation functions for a single boundary
        plan * fields

        fields[:, :, :, rightmoving] .*= phase.k_prop
        fields[:, :, :, leftmoving] .*= phase.k_prop

        i_plan * fields

        if immediatesum
			Eout[:,:,1]    .+= (1+bdry.r[end]) .* fields[:, :, end, rightmoving]
		else
			Eout[:,:,n+1]    = (1+bdry.r[end]) .* fields[:, :, end, rightmoving]
		end

        fields[phase.x_sub, phase.y_sub, :, rightmoving] .*= phase.surface
        fields[phase.x_sub, phase.y_sub, :, leftmoving] .*= phase.surface
"""

		for i in 1:1:length(bdry.distance)
			prop(view(fields, :, :, i, rightmoving), i, coords, phase, plan, i_plan)
			prop(view(fields, :, :, i, leftmoving ), i, coords, phase, plan, i_plan)
		end

		if immediatesum
			Eout[:,:,1]    .+= (1+bdry.r[end]) .* fields[:, :, end, rightmoving]
		else
			Eout[:,:,n+1]    = (1+bdry.r[end]) .* fields[:, :, end, rightmoving]
		end

        fields_arrived_left = fields[:, :, 1:end-1, rightmoving]

        for n in 1:1:(length(bdry.distance) - 1)
            fields[:, :, n, rightmoving] .*= bdry.r[n + 1]
            fields[:, :, n, rightmoving] .+= (1 - bdry.r[n + 1]) .* fields[:, :, n + 1, leftmoving]

            fields[:, :, n + 1, leftmoving] .*= -bdry.r[n + 1]
            fields[:, :, n + 1, leftmoving] .+= (1 + bdry.r[n + 1]) .* fields_arrived_left[:, :, n]
        end

        fields[:, :, end, rightmoving] .*= bdry.r[end]
        fields[:, :, 1, leftmoving] .*= -bdry.r[1]

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

	Eout .*= phase.surface_out[:,:,1] #surface misalignments for the last boundary outside the disk

    # Dance Outro
    # The summation at the end is to give the opportunity to intelligently sum
    # since different rays might have very different orders of magnitude
    # The other reason is, that we are able to return the different components
    # such that we can get immediate access to less accurate results
    # TODO: See if Julia is really using the most intelligent summing strategy
    if returnsum
        return sum(reverse(Eout, dims = 3), dims = 3)[:,:,1]
    else
        return Eout
    end
end

"""
    dance_intro(bdry::SetupBoundaries, X, Y; bfield=nothing, velocity_x=0, f=10e9,diskR=0.1)

Initialize EM fields. Can include velocity effects.
"""
function dance_intro(bdry::SetupBoundaries, coords::CoordinateSystem; bfield=nothing, velocity_x=0, f=10e9,diskR=0.1)
    # Initialize the variable we want to return
    fields_initial = Array{Complex{Float64}}(zeros(length(coords.X), length(coords.Y), length(bdry.distance), 2))

    # Inaccuracies of the emitted fields, BField and Velocity Effects ###################
    if bfield === nothing
        # Make sure that there is only emission on the disk surfaces
        bfield = [x^2 + y^2 > diskR^2 ? 0.0 : 1.0 for x in coords.X, y in coords.Y, z in ones(length(bdry.distance)+1)]
    end
    if velocity_x != 0
        lambda = wavelength(f)
        Ma_PerMeter = 2pi/lambda # k = 2pi/lambda (c/f = lambda)
        bfield .*= [exp(-1im*Ma_PerMeter*(-velocity_x)*x) for z in ones(x in coords.X, y in coords.Y, length(bdry.distance)+1)]
    end
    ####################################################################################

    # Iterate over the gaps and initialize the emissions from them #####################
    # This implements Theoretical Foundations (arxiv: 1612.07057) (3.3)
    for n in 1:length(bdry.distance)
        ax_rightmoving = 0
        if n == 1
            if bdry.r[1] == 1
                ax_rightmoving = -1.0
            else #TODO
                ax_rightmoving = 0
            end
        else
            # Right-Moving in that gap
            eps_i = bdry.eps[n-1]
            eps_m = bdry.eps[n]
            denominator = eps_i * sqrt(eps_m) + eps_m * sqrt(eps_i)
            #ax_i = sqrt(eps_i) * (1 - eps_m/eps_i) / denominator + 0im
            ax_rightmoving = sqrt(eps_m) * (1 - eps_i/eps_m) / denominator +0.0im
        end

        # Left-Moving in that gap
        eps_i = bdry.eps[n]
        eps_m = (n == length(bdry.distance)) ? 1 : bdry.eps[n+1] # Rightmost epsilon is 1.
        denominator = eps_i * sqrt(eps_m) + eps_m * sqrt(eps_i)
        ax_leftmoving = sqrt(eps_i) * (1 - eps_m/eps_i) / denominator +0.0im
        #ax_m = sqrt(eps_m) * (1 - eps_i/eps_m) / denominator +0.0im

        # Fill the initial array
        fields_initial[:,:,n,1] = [1.0*ax_rightmoving + 0.0im for x in coords.X, y in coords.Y] .* bfield[:,:,n]
        fields_initial[:,:,n,2] = [1.0*ax_leftmoving + 0.0im for x in coords.X, y in coords.Y] .* bfield[:,:,n+1]
    end

    return fields_initial
end

"""
    cheerleader(amin, nmax, bdry::SetupBoundaries; f=10.0e9, prop=propagator, emit=nothing, reflect=nothing, Xset=X, Yset=Y, diskR=0.1, returnboth=false)

New Recursive Fourier Propagation implementation.

# Arguments:
- `amin`: Mimum (local) amplitude of a field, in order to be propagated
- `nmax`: Maximum number of beam iteration steps, directly equivalent to how many boundaries a beam 'sees' at most
- `bdry::SetupBoundaries`: Properties of dielectric boundaries
- `f::Float64` ```> 0```: Frequency of EM radiation
- `prop`: Propagator Function to use. Standard is propagator().
- `emit`:  Explicitly set the axion-induced fields emitted from each boundary (to the left and to the right).
               If ``nothing`` fields are initialized according to uniform,
               homogeneous external B-field with zero axion velocity.
- `reflect`: If `nothing` (standar value), the axion-induced signal is computed.
             If set, this field defines a beam, for which the reflected beam will be calculated
- `Xset` and `Yset`: Explicitly set the coordinate system for the fields
- `diskR`: Radius of dielectric disk
- `returnboth::Bool`: If `true` cheerleader returns fields leaving on left and right.
                      If `false` only returns fields leaving on right.

See [`dancer`](@ref) for old version.
"""
function cheerleader(amin, nmax, bdry::SetupBoundaries, coords::CoordinateSystem; f=10.0e9, prop=propagator, emit=nothing, reflect=nothing, diskR=0.1, returnboth=false)

    # Before speed of light was 3e8 here, but this means an error at the permil level, i.e. order ~20MHz at 20GHz,
    # if fixing lambda to 1.5 cm, one gets a shift of roughly 10MHz
    lambda = wavelength(f)

    # Pre-allocate memory
    # Note that fields[0,:,:] contains the fields leaving the system on the left
    # and fields[length(bdry.distance)+1,:,:] the fields leaving on the right
    #fields = OffsetArray(::Array{Complex{Float64}}}, 0:length(bdry.distance)+1, 1:length(X), 1:length(Y))
    fields = Array{Complex{Float64}}(zeros(length(coords.X), length(coords.Y), length(bdry.distance)+2))
    # dimensions of the fields at each position --^
    # number of regions + 2 outporpagating ------------------------------------^

    # In a next step this could/should be generalized in the SetupBoundaries structure..
    reflectivities_leftmoving = -bdry.r
    reflectivities_rightmoving = bdry.r
    transmissivities_leftmoving = 1 .- bdry.r
    transmissivities_rightmoving = 1 .+ bdry.r

	#initialize the arrays with the propagator phases
	phase = init_prop(bdry, coords, lambda, diskR, output_surfaces = [1, length(bdry.distance)])

	#create fft plans for single-slice transforms
	plan = plan_fft!(fields[:,:,1])
	inverse_plan = plan_ifft!(fields[:,:,1])

    ### Indexing of arrays, note the specific indexing of the regions, different from dancer()!.
    #Boundaries:     [ 1,      2,      3,      4 ]
    #          ####### |       |       |       |   #######
    #          # <-- # |       |       |       |   # --> #
    #          ####### |       |       |       |   #######
    #Regions:    [  1  ,   2   ,   3   ,   4   ,   5  ]


    # emit in this function expects for each region left- and right-propagating fields.
    # Note that this is different from dancer() since here we do not take the mirror
    # explicit, such that also a transmissivity can be calculated


    # TODO: propagation through and emission from last bdry to the right
    if reflect === nothing
        if emit === nothing
            emit = Array{Complex{Float64}}(zeros(length(coords.X), length(coords.Y), length(bdry.distance), 2))
            # we may reuse the dance_intro-function in the standard case
            bdrycpy = deepcopy(bdry)
            bdrycpy.distance = bdrycpy.distance[2:end]
            bdrycpy.eps = bdrycpy.eps[2:end]
            emit[:,:,2:end,:] = dance_intro(bdrycpy,coords,diskR=diskR)
        end
    else
        emit = Array{Complex{Float64}}(zeros(length(coords.X), length(coords.Y), length(bdry.distance), 2))
        # This should be right like this, it is explicitly taken care of at the bottom...
        emit[:, :, length(bdry.distance), 2] = reflect
    end

    # Initialize fields ........................................

    # Add everything which is emitted to the right
    fields = copy(emit[:,:,:,1])

    # Push forward the fields rightmoving...
    for i in 2:1:length(bdry.distance)-1
        # Propagate to the right
        prop(view(fields,:,:,i), i, coords, phase, plan, inverse_plan)
        # Reflect and Transmit
        fields[:,:,i+1] .+= transmissivities_rightmoving[i].*fields[:,:,i]
        fields[:,:,i]   .*= reflectivities_rightmoving[i]
    end

    # Add everything which is emitted to the left (except last gap)
    fields[:,:,1:length(bdry.distance)-1] .+= emit[:,:,1:length(bdry.distance)-1,2]
    #fields .+= emit[:,:,:,2]

    # The last gap is not supposed to contain left-propagating stuff,
    # so we have to treat it explicitly (not elegant, but ok)
    let i = length(bdry.distance)-1
        # Transmit
        fields[:,:,i]   .+= transmissivities_leftmoving[i].*emit[:,:,i+1,2]
        # Reflect
        fields[:,:,i+1] .+= reflectivities_leftmoving[i].*emit[:,:,i+1,2]
    end

    # Main iterations ...............................................
    n = 0   # Iteration Count
    a = 1.0 # Maximum Field Amplitude in the iteration

    while n <= nmax && a >= amin
        # Iterate over regions

        # The first region explicitly
        let i = 1
            # Do NOT Propagate to the fields in the current gap to the right
            # Region 1 always contains left-moving fields.
            # Propagate to the fields in the next gap to the left
            prop(view(fields,:,:,i+1), i+1, coords, phase, plan, inverse_plan)
            # Reflect and Transmit
            # No backup copy needed since nothing is coming from the i = 1 region
            # Transmit
            fields[:,:,i]   .+= transmissivities_leftmoving[i].*fields[:,:,i+1]
            # Reflect
            fields[:,:,i+1] .*= reflectivities_leftmoving[i]
        end

        # Push forward...
        for i in 2:1:length(bdry.distance)-2
            # Propagate to the fields in the current gap to the right
            prop(view(fields,:,:,i), i, coords, phase, plan, inverse_plan)
            # Propagate to the fields in the next gap to the left
            prop(view(fields,:,:,i+1), i+1, coords, phase, plan, inverse_plan)
            # Reflect and Transmit
            fields_next_gap_copy = copy(fields[:,:,i+1])
            fields[:,:,i+1] .*= reflectivities_leftmoving[i]
            fields[:,:,i+1] .+= transmissivities_rightmoving[i].*fields[:,:,i]
            fields[:,:,i]   .*= reflectivities_rightmoving[i]
            fields[:,:,i]   .+= transmissivities_leftmoving[i].*fields_next_gap_copy
        end

        let i = length(bdry.distance)-1
            # Propagate to the fields in the current gap to the right
            prop(view(fields,:,:,i), i, coords, phase, plan, inverse_plan)
            # DO NOT Propagate to the fields in the next gap to the left
            # Since it only contains fields propagating to the right
            # Reflect and Transmit
            fields[:,:,i+1] .+= transmissivities_rightmoving[i].*fields[:,:,i]
            fields[:,:,i]   .*= reflectivities_rightmoving[i]
        end

        # The last region always contains right-moving fields and is automatically added up correctly

        # Check convergence
        if amin != 0
        	a = maximum(abs.(fields))
		end
        # n is the iteration count
        n += 1
    end

	fields[:,:,[1, length(bdry.distance)]] .*= phase.surface_out

    if returnboth
        return fields[:,:,1], fields[:,:,length(bdry.distance)]
    else
        return fields[:,:,length(bdry.distance)]
    end
end
