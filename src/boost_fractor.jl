
#
#   BoostFractor Beam Propagation Core + PaD addition (work-in-progress)
#
#   V: 2019-04-05 (core) + 2019-06-09 (add-on)
#
#   Stefan Knirck
#		add-on: Dominik Bergermann

export SetupBoundaries, SeedSetupBoundaries, dancer, dancer_pad, dance_intro, propagator, propagator1D,init_coords


mutable struct SetupBoundaries
    distance::Array{Float64,1} # = [15e-3, 5e-3,0]
    # Boundary reflection coefficient for right-propagating wave
    r::Array{Complex{Float64},1}   # = [1,-0.5,0.5,0]
    # Epsilon to the right of each boundary 
    eps::Array{Complex{Float64},1}   # = [1,9,1]
    # Tilts
    relative_tilt_x # = [0,0]
    relative_tilt_y # = [0,0]
    # Surface Roughness ...
    relative_surfaces::Array{Complex{Float64},3} # = [z, x,y ]
    # etc.
end
SetupBoundaries(distance::Array{Float64,1}, r::Array{Complex{Float64},1}, eps::Array{Complex{Float64},1}) = SetupBoundaries(distance,r,eps, zeros(length(distance)), zeros(length(distance)), zeros(length(distance), length(X), length(Y) ));
SetupBoundaries(distance::Array{Float64,1}) = SetupBoundaries(distance,[1,-0.5,0.5,0],[1,9,1], [0.0,0.0,0.0], [0.0,0.0,0.0]);
SetupBoundaries() = SetupBoundaries([15e-3, 5e-3,0]);

mutable struct DiskDefiniton
    thickness::Float64 # = 1e-3
    # Boundary reflection coefficient for right-propagating wave
    eps::Complex{Float64} # = 9
end
DiskDefiniton() = DiskDefiniton(1e-3, 9)

## Convenient tools ################################################################################
function SeedSetupBoundaries(diskno=3)
    
    distances = [ x % 2 == 1 ? 8e-3 : 1e-3 for x in 1:2*(diskno) ]
    append!(distances, 0e-3) #8e-3)
    
    reflectivities = [1.0]
    append!(reflectivities, [ x % 2 == 1 ? -0.5 : 0.5 for x in 1:2*(diskno) ])
    append!(reflectivities, 0)
    
    epsilon = [ x % 2 == 1.0 ? 1.0 : 9.0 for x in 1:2*(diskno) ]
    append!(epsilon, 1.0)
    
    return SetupBoundaries(distances, Array{Complex{Float64}}(reflectivities), Array{Complex{Float64}}(epsilon))
end


## The heart of it #################################################################################
global X = -0.5:0.01:0.5
global Y = -0.5:0.01:0.5

global minimum_Kx = 2*pi/(maximum(X)*2)
global maximum_Kx = minimum_Kx * (length(X)-1)/2
global coordsKx = -maximum_Kx:minimum_Kx:maximum_Kx
global minimum_Ky = 2*pi/(maximum(Y)*2)
global maximum_Ky = minimum_Ky * (length(Y)-1)/2
global coordsKy = -maximum_Ky:minimum_Ky:maximum_Ky

function init_coords(Xset,Yset)	
    global X=Xset
    global Y=Yset
    global minimum_Kx = 2*pi/(maximum(X)*2)
    global maximum_Kx = minimum_Kx * (length(X)-1)/2
    global coordsKx = -maximum_Kx:minimum_Kx:maximum_Kx
    global minimum_Ky = 2*pi/(maximum(Y)*2)
    global maximum_Ky = minimum_Ky * (length(Y)-1)/2
    global coordsKy = -maximum_Ky:minimum_Ky:maximum_Ky
end

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
		#print("\r Subprogress $(round(n*100/nmax,digits=1))% \r") #neat for unparalleled calcs
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

#Implementation of the propagation around the disc. To be used in the same way as the dancer-function.
#NOTE: Only supports untilted discs with no surface roughness as of now.

function dancer_pad(amin, nmax, bdry::SetupBoundaries; f=10.0e9, prop=propagator, emit=nothing, reflect=nothing, Xset=X, Yset=Y, diskR=0.1, returnsum=true, immediatesum=true)
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

    out = [abs(x^2 + y^2) >= diskR^2 for x in X, y in Y] #required for cutting
    ins = [abs(x^2 + y^2) <  diskR^2 for x in X, y in Y]
    for j in 1:length(bdry.distance), k in 1:2
        fields[j,k,:,:] .*= ins #axion induced fields are generated across whole grid. cut the outer part
    end
    # The measure of music (Takt) .... da capo al fine
    while n <= nmax && a >= amin
        backupfields = copy(fields) #part that gets propagated outside, only outer part required (cut when needed, faster. optimisable?)
        for i in 1:length(bdry.r)
            if i > 1
                if abs(bdry.eps[i-1]) > 1 #cut accordingly if inside of of disc
                    fields[i-1,rightmoving,:,:] .*= ins
                    fields[i-1,rightmoving,:,:] = pfad(fields[i-1,rightmoving,:,:], bdry.distance[i-1], diskR, bdry.eps[i-1], bdry.relative_tilt_x[i-1], bdry.relative_tilt_y[i-1], bdry.relative_surfaces[i-1,:,:], lambda) .* ins
                else
                    fields[i-1,rightmoving,:,:] = pfad(fields[i-1,rightmoving,:,:], bdry.distance[i-1], diskR, bdry.eps[i-1], bdry.relative_tilt_x[i-1], bdry.relative_tilt_y[i-1], bdry.relative_surfaces[i-1,:,:], lambda)
                end
            end
            if i < length(bdry.r)
                if abs(bdry.eps[i]) > 1 #needs to be adjusted if planning on having other medium than vacuum around discs
                    fields[i,leftmoving,:,:] .*= ins
                    fields[i,leftmoving,:,:] = pfad(fields[i,leftmoving,:,:], bdry.distance[i], diskR, bdry.eps[i], bdry.relative_tilt_x[i], bdry.relative_tilt_y[i], bdry.relative_surfaces[i,:,:], lambda) .* ins
                else
                    fields[i,leftmoving,:,:] = pfad(fields[i,leftmoving,:,:], bdry.distance[i], diskR, bdry.eps[i], bdry.relative_tilt_x[i], bdry.relative_tilt_y[i], bdry.relative_surfaces[i,:,:], lambda)
                end
            end
            # Reflection & Transmission
            if i == 1 #leftmost
                fields[i,leftmoving,:,:] .*= -bdry.r[i]
            elseif i == length(bdry.r) #rightmost
                if immediatesum             
                    Eout[:,:,1] .+= (1+bdry.r[i]) .* fields[i-1, rightmoving,:,:]
                else            
                    Eout[:,:,n+1] = (1+bdry.r[i]) .* fields[i-1, rightmoving,:,:]
                end
                fields[i-1,rightmoving,:,: ] .*= bdry.r[i]
            else #standard
                FieldArrivedLeft = copy(fields[i-1, rightmoving,:,:])
                if abs(bdry.eps[i]) > 1
                    fields[i-1, rightmoving,:,:] .*= ins #only the inner part
                    fields[i-1, rightmoving,:,:] .*= bdry.r[i] #gets reflected
                    fields[i-1, rightmoving,:,:] .+= (1-bdry.r[i]) .* fields[i, leftmoving,:,:] 
                    fields[i-1, rightmoving,:,:] .+= backupfields[i,leftmoving,:,:] .* exp(-1im*conj(sqrt((2pi/lambda)^2))*bdry.distance[i]*sqrt(1)) .* out

                    fields[i,leftmoving,:,:] .*= -bdry.r[i] #as this is already only the inner part (see above), no need to seperate
                    fields[i,leftmoving,:,:] .+= (1+bdry.r[i]) .* FieldArrivedLeft .* ins
                    fields[i,leftmoving,:,:] .+= FieldArrivedLeft .* out
                else
                    fields[i-1,rightmoving,:,:] .*= bdry.r[i]
                    fields[i-1,rightmoving,:,:] .+= (1-bdry.r[i]) .* fields[i,leftmoving,:,:] .* ins
                    fields[i-1,rightmoving,:,:] .+= fields[i,leftmoving,:,:] .* out

                    fields[i,leftmoving,:,:] .*= ins
                    fields[i,leftmoving,:,:] .*= -bdry.r[i]
                    fields[i,leftmoving,:,:] .+= (1+bdry.r[i]) .* FieldArrivedLeft .* ins
                    fields[i,leftmoving,:,:] .+= backupfields[i-1,rightmoving,:,:] .* exp(-1im*conj(sqrt((2pi/lambda)^2))*bdry.distance[i-1]*sqrt(1)) .* out
                end
            end
        end
        lmv = copy(leftmoving) #direction change
        leftmoving = rightmoving
        rightmoving = lmv
        if amin != 0
            a = maximum(abs.(fields))
        end
        n += 1
        #print("\r Subprogress $(round(n*100/nmax,digits=1))% \r") #neat for unparalleled calcs
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

function dance_intro(bdry::SetupBoundaries, X, Y; bfield=nothing, velocity_x=0, f=10e9)
    # Initialize the variable we want to return 
    fields_initial = Array{Complex{Float64}}(zeros(length(bdry.distance), 2, length(X), length(Y)))

    # Inaccuracies of the emitted fields, BField and Velocity Effects ###################
    if bfield == nothing
        bfield = Array{Complex{Float64}}(ones(length(bdry.distance)+1, length(X), length(Y)))
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

function initialize_reflection_transmission(freq::Float64, bdry::SetupBoundaries, disk::DiskDefiniton)
    if disk == nothing
        # Initilize reflection coefficients according to epsilon
        r_left = ones(length(bdry.eps))
        r_left[1] = -1
        for i in 2:length(bdry.eps)
            # The reflectivity at this point
            r_left = (sqrt(bdry.eps[i-1])-sqrt(bdry.eps[i]))/(sqrt(bdry.eps[i])+sqrt(bdry.eps[i-1]))
        end
        r_right = -r_left
        t_left = 1. + r_left
        t_right = 1 .+ r_right
    else
        # Initialize reflection coefficients according to disk model
        ref, trans = reflectivity_transmissivity_1d(freq, disk.thickness)
        r_left = ones(length(bdry.eps),length(coordsKx),length(coordsKy)).*ref
        r_right = r_left
        t_left = ones(length(bdry.eps),length(coordsKx),length(coordsKy)).*trans
        t_right = t_left
    end
    return r_left, r_right, t_left, t_right
end



## Propagators ########################################################################################

"""
Does the FFT of E0 on a disk, propagates the beam a given distance and does the iFFT
Note that this method is in-place. If it should be called more than one time on the 
same fields, use propagator(copy(E0), ...).
"""
function propagator(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda)    
    k0 = 2*pi/lambda*sqrt(eps)
    # Call the propagator and add a phase imposed by the tilt
    # Assumptions: Tilt is small, the additional phase is obtained by propagating
    #              all fields just with k0 to the tilted surface (only valid if diffraction effects are small)
    E0 = propagatorNoTilts(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda)
    # Tilts:
    E0 .*= [exp(-1im*k0*tilt_x*x) * exp(-1im*k0*tilt_y*y) for x in X, y in Y]
    # More general: Any surface misalignments:
    E0 .*= exp.(-1im*k0*surface) #(the element wise (exp.) is important, otherwise "surface" is treated as a matrix!)
    return E0
end

function propagatorNoTilts(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda)
    # Diffract at the Disk. Only the disk is diffracting.
    E0 .*= [abs(x^2 + y^2) < diskR^2 for x in X, y in Y]
    # FFT the E-Field to spatial frequencies
    #print(E0)
    FFTW.fft!(E0)
    #print(E0)
    E0 = FFTW.fftshift(E0)

    # This should now work with the global variable 
    #minimum_Kx = 2*pi/(maximum(X)*2)
    #maximum_Kx = minimum_Kx * (length(X)-1)/2
    #coordsKx = -maximum_Kx:minimum_Kx:maximum_Kx
    #minimum_Ky = 2*pi/(maximum(Y)*2)
    #maximum_Ky = minimum_Ky * (length(Y)-1)/2
    #coordsKy = -maximum_Ky:minimum_Ky:maximum_Ky
    
    # TODO: If maximum k is higher than k0, then it is not defined
    #       what happens with this mode
    #       We should give a warning and handle this here
    #       At the moment the script will just also propagate with a loss for those components
    # Propagate through space
    k0 = 2*pi/lambda*sqrt(eps)
    k_prop = [conj(sqrt( Complex{Float64}(k0^2 - Kx^2 - Ky^2) )) for Kx in coordsKx, Ky in coordsKy]
    E0 = E0 .* exp.(-1im*k_prop*dz)
    # Backtransform
    E0 = FFTW.ifftshift(E0)
    FFTW.ifft!(E0)
    return E0
end

"""
Propagator that assumes E0 is already in momentum space.
"""
function propagatorMomentumSpace(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda)
    # Propagate through space
    k0 = 2*pi/lambda*sqrt(eps)
    k_prop = [conj(sqrt( Complex{Float64}(k0^2 - Kx^2 - Ky^2) )) for Kx in coordsKx, Ky in coordsKy]
    E0 = E0 .* exp.(-1im*k_prop*dz)

    # Transform to position space
    E0 = FFTW.ifftshift(E0)
    FFTW.ifft!(E0)

    # Diffract at the Disk. Only the disk is diffracting.   
    E0 .*= [abs(x^2 + y^2) < diskR^2 for x in X, y in Y]

    # Kick (tilt)
    if tilt_x != 0 || tilt_y != 0
        E0 .*= [exp(-1im*k0*tilt_x*x) * exp(-1im*k0*tilt_y*y) for x in X, y in Y]
    end
    
    # FFT the E-Field to spatial frequencies / momentum space
    FFTW.fft!(E0)
    E0 = FFTW.fftshift(E0)

    return E0
end

"""
This propagator just does the phase propagation.
"""
function propagator1D(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda)
    # Version of the propagator without the fft
    # should be faster and easy to check consistency with 1D calc
    
    # Propagate through space
    k0 = 2*pi/lambda
    k_prop = conj(sqrt.(k0^2))
    e1 = E0.*exp(-1im*k_prop*dz*sqrt(eps))
    return e1
    
end

"""
This propagator is used for propagation around the disc.
"""
function pfad(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda) #prop for around disk
    k0 = 2*pi/lambda*sqrt(eps)                                     #performs propagator and propagatorNoTilts w/o cutting
    FFTW.fft!(E0)                                                  #necessary cutting is performed in the dancer_pad
    E0 = FFTW.fftshift(E0)
    k_prop = [conj(sqrt( Complex{Float64}(k0^2 - Kx^2 - Ky^2) )) for Kx in coordsKx, Ky in coordsKy]
    E0 = E0 .* exp.(-1im*k_prop*dz)
    E0 = FFTW.ifftshift(E0)
    FFTW.ifft!(E0)
    E0 .*= [exp(-1im*k0*tilt_x*x) * exp(-1im*k0*tilt_y*y) for x in X, y in Y]
    E0 .*= exp.(-1im*k0*surface) #(the element wise (exp.) is important, otherwise "surface" is treated as a matrix!)
    return E0
end
