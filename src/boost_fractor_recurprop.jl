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
function dancer(amin, nmax, bdry::SetupBoundaries; f=10.0e9, prop=propagator, emit=nothing, reflect=nothing, Xset=X, Yset=Y, diskR=0.1, returnsum=true, immediatesum=true)
    init_coords(Xset, Yset);

    rightmoving = 1
    leftmoving = 2

    fields = Array{Complex{Float64}}(zeros(length(bdry.distance), 2, length(X), length(Y)))
    #fields = SharedArray(Complex{Float64},  length(bdry.distance), 2, length(X), length(Y))
    #                       number of regions -----^                ^              ^
    #                       number of propagation directions -------^              ^
    #                       dimensions of the fields at each position -------------^

    # Pre-allocate memory
    if immediatesum
        Eout = Array{Complex{Float64}}(zeros(length(X), length(Y),1))
    else
        Eout = Array{Complex{Float64}}(zeros(length(X), length(Y),nmax+1))
    end

    # TODO: propagation through and emission from last bdry to the right
    if reflect == nothing
        if emit == nothing
            fields = dance_intro(bdry,X,Y;diskR=diskR)
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

    # The measure of music (Takt) .... da capo al fine
    while n <= nmax && a >= amin
        # The "legs" of the dancer dancing in parallel
        # @sync is important to make sure that the outer loop waits for the inner one to finish...
        #@sync @parallel
        for i in 1:1:length(bdry.r)
            # The beats

            # Move left leg (rightmoving)
            if i > 1
                # In place propagate left leg:
                fields[i-1,rightmoving,:,:] = prop(fields[i-1,rightmoving,:,:], bdry.distance[i-1], diskR, bdry.eps[i-1], bdry.relative_tilt_x[i-1], bdry.relative_tilt_y[i-1], bdry.relative_surfaces[i-1,:,:], lambda)
            end

            # Move right leg (leftmoving)
            if i < length(bdry.r)
                fields[i,leftmoving,:,:] = prop(fields[i,leftmoving,:,:], bdry.distance[i], diskR, bdry.eps[i], bdry.relative_tilt_x[i], bdry.relative_tilt_y[i], bdry.relative_surfaces[i,:,:], lambda)
            end

            # Reflection & Transmission
            if i == 1 # Leftmost case
                fields[i, leftmoving,:,:] .*= -bdry.r[i]
                # The right-moving does not exist and about the transmitted one we dont care (yet)
            elseif i == length(bdry.r) # Rightmost case
                # The rightmoving is transmitted to Eout (yeah! that'll be our result!)
                # old version without pre-allocation
                #Eout                         = cat(3,Eout, (1+bdry.r[i])*fields[i-1, rightmoving,:,:])
                # with pre-allocated array:
                if immediatesum
                    Eout[:,:,1]              .+= (1+bdry.r[i]).*fields[i-1, rightmoving,:,:]
                else
                    Eout[:,:,n+1]              = (1+bdry.r[i]).*fields[i-1, rightmoving,:,:]
                end
                # And reflected as well.
                fields[i-1, rightmoving,:,:] .*= bdry.r[i]
            else # Standard Case

                # field(i-1) = transmit(field(i)) + reflect(field(i-1))
                # field(i) = transmit(field(i-1)) + reflect(field(i))
                # Basically the above, but as much in-place as possible

                # Always the field on the left gets a copy
                FieldArrivedLeft = copy(fields[i-1, rightmoving,:,:])
                # Then we in-place update it
                fields[i-1, rightmoving,:,:] .*= bdry.r[i]
                fields[i-1, rightmoving,:,:] .+= (1-bdry.r[i]).*fields[i, leftmoving,:,:]
                # Then we in-place update the other one with the previously made backup
                fields[i, leftmoving,:,:] .*= -bdry.r[i]
                fields[i, leftmoving,:,:] .+= (1+bdry.r[i]).*FieldArrivedLeft
            end
        end

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

function dance_intro(bdry::SetupBoundaries, X, Y; bfield=nothing, velocity_x=0, f=10e9,diskR=0.1)
    # Initialize the variable we want to return
    fields_initial = Array{Complex{Float64}}(zeros(length(bdry.distance), 2, length(X), length(Y)))

    # Inaccuracies of the emitted fields, BField and Velocity Effects ###################
    if bfield == nothing
        # Make sure that there is only emission on the disk surfaces
        bfield = [x^2 + y^2 > diskR^2 ? 0.0 : 1.0 for z in ones(length(bdry.distance)+1), x in X, y in Y] 
    end
    if velocity_x != 0
        c = 299792458.
        Ma_PerMeter = 2pi*f/c # k = 2pi/lambda (c/f = lambda)
        bfield .*= [exp(-1im*Ma_PerMeter*(-velocity_x)*x) for z in ones(length(bdry.distance)+1), x in X, y in Y]
    end
    ####################################################################################

    # Iterate over the gaps and initialize the emissions from them #####################
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

        # Right-Moving in that gap
        eps_i = bdry.eps[n]
        eps_m = (n == length(bdry.distance)) ? 1 : bdry.eps[n+1] # Rightmost epsilon is 1.
        denominator = eps_i * sqrt(eps_m) + eps_m * sqrt(eps_i)
        ax_leftmoving = sqrt(eps_i) * (1 - eps_m/eps_i) / denominator +0.0im
        #ax_m = sqrt(eps_m) * (1 - eps_i/eps_m) / denominator +0.0im

        # Fill the initial array
        fields_initial[n,1,:,:] = [1.0*ax_rightmoving + 0.0im for x in X, y in Y] .* bfield[n,:,:]
        fields_initial[n,2,:,:] = [1.0*ax_leftmoving + 0.0im for x in X, y in Y] .* bfield[n+1,:,:]
    end

    return fields_initial
end

function cheerleader(amin, nmax, bdry::SetupBoundaries; f=10.0e9, prop=propagator, emit=nothing, reflect=nothing, Xset=X, Yset=Y, diskR=0.1, returnboth=false)
    # Before speed of light was 3e8 here, but this means an error at the permil level, i.e. order ~20MHz at 20GHz,
    # if fixing lambda to 1.5 cm, one gets a shift of roughly 10MHz
    c = 299792458.
    lambda = c/f

    init_coords(Xset, Yset);

    # Pre-allocate memory
    # Note that fields[0,:,:] contains the fields leaving the system on the left
    # and fields[length(bdry.distance)+1,:,:] the fields leaving on the right
    #fields = OffsetArray(::Array{Complex{Float64}}}, 0:length(bdry.distance)+1, 1:length(X), 1:length(Y))
    fields = Array{Complex{Float64}}(zeros(     length(bdry.distance)+2, length(X), length(Y)))
    # number of regions + 2 outporpagating --^                 ^
    # dimensions of the fields at each position ---------------^

    # In a next step this could/should be generalized in the SetupBoundaries structure..
    reflectivities_leftmoving =  -bdry.r
    reflectivities_rightmoving = bdry.r
    transmissivities_leftmoving = 1 .- bdry.r
    transmissivities_rightmoving = 1 .+ bdry.r

    ### Indexing of arrays, note the specific indexing of the regions, different from dancer()!.
    #Boundaries:     [ 1,      2,      3,      4 ]
    #          ####### |       |       |       |   #######
    #          # <-- # |       |       |       |   # --> #
    #          ####### |       |       |       |   #######
    #Regions:    [  1  ,   2   ,   3   ,   4   ,   5  ]


    # emit in this function expects for each region left- and right-propagating fields.
    # Note that this is different from dancer() since here we do not take the mirror
    # explicit, such that also a transmissivity can be calculated

    #emit = zeros(length(bdry.distance), 2, length(X), length(Y))
    # TODO: propagation through and emission from last bdry to the right
    if reflect == nothing
        if emit == nothing
            emit = Array{Complex{Float64}}(zeros(length(bdry.distance), 2, length(X), length(Y)))
            # we may reuse the dance_intro-function in the standard case
            bdrycpy = deepcopy(bdry)
            bdrycpy.distance = bdrycpy.distance[2:end]
            bdrycpy.eps = bdrycpy.eps[2:end]
            emit[2:end,:,:,:] = dance_intro(bdrycpy,X,Y,diskR=diskR)
        end
    else
        emit = Array{Complex{Float64}}(zeros(length(bdry.distance), 2, length(X), length(Y)))
        # This should be right like this, it is explicitly taken care of at the bottom...
        emit[length(bdry.distance), 2, :, :] = reflect
    end

    # Initialize fields ........................................

    # Add everything which is emitted to the right
    fields = copy(emit[:,1,:,:])

    # Push forward the fields rightmoving...
    for i in 2:1:length(bdry.distance)-1
        # Propagate to the right
        fields[i,:,:] = prop(fields[i,:,:], bdry.distance[i], diskR, bdry.eps[i], bdry.relative_tilt_x[i], bdry.relative_tilt_y[i], bdry.relative_surfaces[i,:,:], lambda)
        # Reflect and Transmit
        fields[i+1,:,:] .+= transmissivities_rightmoving[i].*fields[i,:,:]
        fields[i,:,:]   .*= reflectivities_rightmoving[i]
    end

    # Add everything which is emitted to the left (except last gap)
    fields[1:length(bdry.distance)-1,:,:] .+= emit[1:length(bdry.distance)-1,2,:,:]
    #fields .+= emit[:,2,:,:]

    # The last gap is not supposed to contain left-propagating stuff,
    # so we have to treat it explicitly (not elegant, but ok)
    let i = length(bdry.distance)-1
        # Transmit
        fields[i,:,:]   .+= transmissivities_leftmoving[i].*emit[i+1,2,:,:]
        # Reflect
        fields[i+1,:,:] .+= reflectivities_leftmoving[i].*emit[i+1,2,:,:]
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
            fields[i+1,:,:] = prop(fields[i+1,:,:], bdry.distance[i+1], diskR, bdry.eps[i+1], bdry.relative_tilt_x[i+1], bdry.relative_tilt_y[i+1], bdry.relative_surfaces[i+1,:,:], lambda)
            # Reflect and Transmit
            # No backup copy needed since nothing is coming from the i = 1 region
            # Transmit
            fields[i,:,:]   .+= transmissivities_leftmoving[i].*fields[i+1,:,:]
            # Reflect
            fields[i+1,:,:] .*= reflectivities_leftmoving[i]
        end

        # Push forward...
        for i in 2:1:length(bdry.distance)-2
            # Propagate to the fields in the current gap to the right
            fields[i,:,:] = prop(fields[i,:,:], bdry.distance[i], diskR, bdry.eps[i], bdry.relative_tilt_x[i], bdry.relative_tilt_y[i], bdry.relative_surfaces[i,:,:], lambda)
            # Propagate to the fields in the next gap to the left
            fields[i+1,:,:] = prop(fields[i+1,:,:], bdry.distance[i+1], diskR, bdry.eps[i+1], bdry.relative_tilt_x[i+1], bdry.relative_tilt_y[i+1], bdry.relative_surfaces[i+1,:,:], lambda)
            # Reflect and Transmit
            fields_next_gap_copy = copy(fields[i+1,:,:])
            fields[i+1,:,:] .*= reflectivities_leftmoving[i]
            fields[i+1,:,:] .+= transmissivities_rightmoving[i].*fields[i,:,:]
            fields[i,:,:]   .*= reflectivities_rightmoving[i]
            fields[i,:,:]   .+= transmissivities_leftmoving[i].*fields_next_gap_copy
        end

        let i = length(bdry.distance)-1
            # Propagate to the fields in the current gap to the right
            fields[i,:,:] = prop(fields[i,:,:], bdry.distance[i], diskR, bdry.eps[i], bdry.relative_tilt_x[i], bdry.relative_tilt_y[i], bdry.relative_surfaces[i,:,:], lambda)
            # DO NOT Propagate to the fields in the next gap to the left
            # Since it only contains fields propagating to the right
            # Reflect and Transmit
            fields[i+1,:,:] .+= transmissivities_rightmoving[i].*fields[i,:,:]
            fields[i,:,:]   .*= reflectivities_rightmoving[i]
        end

        # The last region always contains right-moving fields and is automatically added up correctly

        # Check convergence
        if amin != 0
        	a = maximum(abs.(fields))
		end
        # n is the iteration count
        n += 1
    end


    if returnboth
        return fields[1,:,:], fields[length(bdry.distance),:,:]
    else
        return fields[length(bdry.distance),:,:]
    end
end
