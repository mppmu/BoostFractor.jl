#
#   BoostFractor
#
#   V: 2019-11-29
#
#   Stefan Knirck
#
#
#   This file defines common interfaces for all algorithms, such as (x,y)
#   coordinate system, the SetupBoundaries structure, Propagators, ...
#


export SetupBoundaries, SeedSetupBoundaries, propagator, propagator1D, CoordinateSystem, SeedCoordinateSystem, wavelength

function wavelength(frequency::Float64)
    speed_of_light = 299792458. # [m]
    return speed_of_light / frequency
end

struct CoordinateSystem
    # Real Space Coordinates
    X::Array{Float64,1}
    Y::Array{Float64,1}
    # K space Coordinates
    kX::Array{Float64,1}
    kY::Array{Float64,1}
end

function SeedCoordinateSystem(;X = -0.5:0.01:0.5, Y = -0.5:0.01:0.5) # units [m]
    kX = get_kspace_coords(X)
    kY = get_kspace_coords(Y)
    return CoordinateSystem(X, Y, kX, kY)
end


"""
    get_kspace_coords(RealSpaceCoords)

Calculate k space coordinate system from real space one. Helper function for SeedCoordinateSystem.
"""
function get_kspace_coords(RealSpaceCoords)
    minimum_Kspace = 2*pi/(maximum(RealSpaceCoords)*2)
    maximum_Kspace = minimum_Kspace * (length(RealSpaceCoords)-1)/2
    coordsKspace = -maximum_Kspace:minimum_Kspace:maximum_Kspace
    return coordsKspace
end


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


## Convenient tools ################################################################################

function SeedSetupBoundaries(coords::CoordinateSystem; diskno=3, distance=nothing, reflectivities=nothing, epsilon=nothing, relative_tilt_x=zeros(2*diskno+2), relative_tilt_y=zeros(2*diskno+2), relative_surfaces=zeros(2*diskno+2 , length(coords.X), length(coords.Y)))

    # Initialize SetupBoundaries entries with default values given diskno. Rest was already defined in function definition.
    if distance == nothing # [m]
        distance = [0.0]
        append!(distance, [ x % 2 == 1 ? 8e-3 : 1e-3 for x in 1:2*(diskno) ])
        append!(distance, 0e-3) #8e-3)
    end

    if reflectivities == nothing
        reflectivities = [1.0]
        append!(reflectivities, [ x % 2 == 1 ? -0.5 : 0.5 for x in 1:2*(diskno) ])
    end

    if epsilon == nothing
        epsilon = [NaN]
        append!(epsilon, [ x % 2 == 1.0 ? 1.0 : 9.0 for x in 1:2*(diskno) ])
        append!(epsilon, 1.0)
    end

    # Check if initialization was self-consistent
    length(distance) == length(reflectivities)+1 == length(epsilon) == length(relative_tilt_x) == length(relative_tilt_y) == size(relative_surfaces, 1) || throw(DimensionMismatch("the arrays in your SetupBoundaries objects don't fit together!"))

    return SetupBoundaries(distance, Array{Complex{Float64}}(reflectivities), Array{Complex{Float64}}(epsilon), relative_tilt_x, relative_tilt_y, relative_surfaces)
end


## The heart of it #################################################################################


#TODO: This function is no longer supported. reflectivity_transmissivity_1d does not exist; DiskDefinition has merged into SetuoBoundaries.
"""
OUTDATED! Does not work, do not use!
"""
function initialize_reflection_transmission(freq::Float64, bdry::SetupBoundaries, coords::CoordinateSystem, disk)#::DiskDefiniton)
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
        # Initilize reflection coefficients according to disk model
        ref, trans = reflectivity_transmissivity_1d(freq, disk.thickness)
        r_left = ones(length(bdry.eps),length(coords.kX),length(coords.kY)).*ref
        r_right = r_left
        t_left = ones(length(bdry.eps),length(coords.kX),length(coords.kY)).*trans
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
function propagator(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda, coords::CoordinateSystem)
    k0 = 2*pi/lambda*sqrt(eps)
    # Call the propagator and add a phase imposed by the tilt
    # Assumptions: Tilt is small, the additional phase is obtained by propagating
    #              all fields just with k0 to the tilted surface (only valid if diffraction effects are small)
    E0 = propagatorNoTilts(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda, coords)
    # Tilts:
    E0 .*= [exp(-1im*k0*tilt_x*x) * exp(-1im*k0*tilt_y*y) for x in coords.X, y in coords.Y]
    # More general: Any surface misalignments:
    E0 .*= exp.(-1im*k0*surface) #(the element wise (exp.) is important, otherwise "surface" is treated as a matrix!)
    return E0
end

function propagatorNoTilts(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda, coords::CoordinateSystem)
    # Diffract at the Disk. Only the disk is diffracting.
    E0 .*= [abs(x^2 + y^2) < diskR^2 for x in coords.X, y in coords.Y]
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
    k_prop = [conj(sqrt( Complex{Float64}(k0^2 - Kx^2 - Ky^2) )) for Kx in coords.kX, Ky in coords.kY]
    E0 = E0 .* exp.(-1im*k_prop*dz)
    # Backtransform
    E0 = FFTW.ifftshift(E0)
    FFTW.ifft!(E0)
    return E0
end

"""
Propagator that assumes E0 is already in momentum space.
"""
function propagatorMomentumSpace(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda, coords::CoordinateSystem)
    # Propagate through space
    k0 = 2*pi/lambda*sqrt(eps)
    k_prop = [conj(sqrt( Complex{Float64}(k0^2 - Kx^2 - Ky^2) )) for Kx in coords.kX, Ky in coords.kY]
    E0 = E0 .* exp.(-1im*k_prop*dz)

    # Transform to position space
    E0 = FFTW.ifftshift(E0)
    FFTW.ifft!(E0)

    # Diffract at the Disk. Only the disk is diffracting.
    E0 .*= [abs(x^2 + y^2) < diskR^2 for x in coords.X, y in coords.Y]

    # Kick (tilt)
    if tilt_x != 0 || tilt_y != 0
        E0 .*= [exp(-1im*k0*tilt_x*x) * exp(-1im*k0*tilt_y*y) for x in coords.X, y in coords.Y]
    end

    # FFT the E-Field to spatial frequencies / momentum space
    FFTW.fft!(E0)
    E0 = FFTW.fftshift(E0)

    return E0
end

"""
This propagator just does the phase propagation.
"""
function propagator1D(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda, coords::CoordinateSystem)
    # Version of the propagator without the fft
    # should be faster and easy to check consistency with 1D calc

    # Propagate through space
    k0 = 2*pi/lambda
    k_prop = conj(sqrt.(k0^2))
    e1 = E0.*exp(-1im*k_prop*dz*sqrt(eps))
    return e1

end
