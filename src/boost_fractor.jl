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


export SetupBoundaries, SeedSetupBoundaries, propagator, propagator1D,init_coords

@doc raw"""
# Summary
    mutable struct SetupBoundaries <: Any

Define properties of dielectric boundaries. Coordinate system?

# Fields:
- `distance::Array{Float64,1}` ```> 0```: Distance in z direction to boundary
- `r::Array{Complex{Float64},1}` ```[0, 1]```: Boundary reflection coefficient for right-propagating wave
- `eps::Array{Complex{Float64},1}` ```≥ 1```: Dielectric permittivity to the right of each boundary"
- `relative_tilt_x` ```> 0```: Tilt in x direction [rad?]
- `relative_tilt_y` ```> 0```: Tilt in y direction [rad?]
- `relative_surfaces::Array{Complex{Float64},3}` ```?```: Surface roughness. z offset (1st dim) at x,y (2nd, 3rd dims)
"""
mutable struct SetupBoundaries
    distance::Array{Float64,1} # = [15e-3, 5e-3,0]
    r::Array{Complex{Float64},1}   # = [1,-0.5,0.5,0]
    eps::Array{Complex{Float64},1}   # = [1,9,1]
    relative_tilt_x # = [0,0]
    relative_tilt_y # = [0,0]
    relative_surfaces::Array{Complex{Float64},3} # = [z, x,y ]
    # etc.
end
SetupBoundaries(distance::Array{Float64,1}, r::Array{Complex{Float64},1}, eps::Array{Complex{Float64},1}) = SetupBoundaries(distance,r,eps, zeros(length(distance)), zeros(length(distance)), zeros(length(distance), length(X), length(Y) ));
SetupBoundaries(distance::Array{Float64,1}) = SetupBoundaries(distance,[1,-0.5,0.5,0],[1,9,1], [0.0,0.0,0.0], [0.0,0.0,0.0]);
SetupBoundaries() = SetupBoundaries([15e-3, 5e-3,0]);


@doc raw"""
# Summary
    mutable struct DiskDefinition <: Any

Define properties of dielectric discs (same for all discs).

# Fields:
- `thickness::Float64` ```> 0```: Thickness of discs
- `eps::Complex{Float64}` ```> 1```: Dielectric permittivity
"""
mutable struct DiskDefiniton
    thickness::Float64 # = 1e-3
    # ??? Boundary reflection coefficient for right-propagating wave
    eps::Complex{Float64} # = 9
end
DiskDefiniton() = DiskDefiniton(1e-3, 9)

## Convenient tools ################################################################################
@doc raw"""
    SeedSetupBoundaries(diskno=3)

Initialize `mutable struct SetupBoundaries` with sensible values.

# Arguments
- `diskno::Int` ```> 0```: Number of dielectric discs
"""
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

"""
    init_coords(Xset, Yset)

Initialize coordinate system in real and fourier space.

# Arguments
- `Xset::AbstractRange{Float}`: x coordinates
- `Yset::AbstractRange{Float}`: y coordinates
"""
function init_coords(Xset,Yset)
    #TODO: This function is definitely bad style!
    global X=Xset
    global Y=Yset
    global minimum_Kx = 2*pi/(maximum(X)*2)
    global maximum_Kx = minimum_Kx * (length(X)-1)/2
    global coordsKx = -maximum_Kx:minimum_Kx:maximum_Kx
    global minimum_Ky = 2*pi/(maximum(Y)*2)
    global maximum_Ky = minimum_Ky * (length(Y)-1)/2
    global coordsKy = -maximum_Ky:minimum_Ky:maximum_Ky
end

@doc raw"""
    initialize_reflection_transmission(freq::Float64, bdry::SetupBoundaries, disk::DiskDefiniton)

Calculate reflection and transmission coefficients.

# Arguments
- `freq::Float64` ```> 0```: Frequency of EM radiation
- `bdry::SetupBoundaries`: Properties of dielectric boundaries
- `disk::DiskDefiniton`: Properties of dielectric discs
"""
function initialize_reflection_transmission(freq::Float64, bdry::SetupBoundaries, disk::DiskDefiniton)
    if disk == nothing
        # Iniatilize reflection coefficients according to epsilon
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
        # Initailize reflection coefficients according to disk model
        ref, trans = reflectivity_transmissivity_1d(freq, disk.thickness)
        r_left = ones(length(bdry.eps),length(coordsKx),length(coordsKy)).*ref
        r_right = r_left
        t_left = ones(length(bdry.eps),length(coordsKx),length(coordsKy)).*trans
        t_right = t_left
    end
    return r_left, r_right, t_left, t_right
end

## Propagators ########################################################################################
#TODO: propagator and propagator NoTilts are in-place, julia convention is name! instead of name
"""
    propagator(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda)

Do the FFT of E0 on a disk, propagate the beam a given distance and do the iFFT.
Note that this method is in-place. If it should be called more than one time on the
same fields, use propagator(copy(E0), ...).

Assume: Tilt is small, additional phase is obtained by propagating all fields just
with k0 to the tilted surface (only valid if diffraction effects are small).

# Arguments
- `E0::Array{Float64,1}`: Electric field before propagation
- `dz`: Distance propagated in z direction
- `diskR`: Radius of discs
- `eps`: Dielectric permittivity
- `tilt_x`: Disc tilt in x direction
- `tilt_y`: Disc tilt in y direction
- `surface`: Surface roughness of disc (only at end of propagation)
- `lambda`: Wavelength of electric field

"""
function propagator(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda)
    k0 = 2*pi/lambda*sqrt(eps)
    # Call the propagator and add a phase imposed by the tilt
    E0 = propagatorNoTilts(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda)
    # Tilts:
    E0 .*= [exp(-1im*k0*tilt_x*x) * exp(-1im*k0*tilt_y*y) for x in X, y in Y]
    # More general: Any surface misalignments:
    E0 .*= exp.(-1im*k0*surface) #(the element wise (exp.) is important, otherwise "surface" is treated as a matrix!)
    return E0
end

"""
    propagatorNoTilts(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda)

"""
function propagatorNoTilts(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda)
    # Diffract at the Disk. Only the disk is diffracting.
    E0 .*= [abs(x^2 + y^2) < diskR^2 for x in X, y in Y]
    # FFT the E-Field to spatial frequencies
    #print(E0)
    FFTW.fft!(E0)
    #print(E0)
    E0 = FFTW.fftshift(E0)

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
