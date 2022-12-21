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

@doc raw"""
# Summary
    mutable struct SetupBoundaries <: Any

Define properties of dielectric boundaries. Coordinate system?

# Fields:
- `distance::Array{Float64,1}` ```> 0```: Distance in z direction to boundary
- `r::Array{Complex{Float64},1}` ```[0, 1]```: Boundary reflection coefficient for right-propagating wave
- `eps::Array{Complex{Float64},1}`: Dielectric permittivity to the right of each boundary"
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


## Convenient tools ################################################################################
@doc raw"""
    SeedSetupBoundaries(diskno=3)

Initialize `mutable struct SetupBoundaries` with sensible values.

# Arguments
- `diskno::Int` ```> 0```: Number of dielectric discs
"""
function SeedSetupBoundaries(coords::CoordinateSystem, distance, epsilon; relative_tilt_x=zeros(length(distance)), relative_tilt_y=zeros(length(distance)), relative_surfaces=zeros(length(distance) , length(coords.X), length(coords.Y)))
    epsilon = Array{Complex{Float64}}(epsilon)
    epsilon[imag.(epsilon) .== 0.0] .= real.(epsilon[imag.(epsilon) .== 0.0]) .- 0.0im#Make sure that a real epsilon is casted to Re(eps)-0.0im.
    reflectivities = complex([1.0])
    R = [(sqrt(epsilon[i-1]) - sqrt(epsilon[i])) / (sqrt(epsilon[i-1]) + sqrt(epsilon[i])) for i in 3:length(epsilon)]
    append!(reflectivities, R)
    # Check if initialization was self-consistent
    length(distance) == length(reflectivities)+1 == length(epsilon) == length(relative_tilt_x) == length(relative_tilt_y) == size(relative_surfaces, 1) || throw(DimensionMismatch("the arrays in your SetupBoundaries objects don't fit together!"))
    return SetupBoundaries(distance, Array{Complex{Float64}}(reflectivities), epsilon, relative_tilt_x, relative_tilt_y, relative_surfaces)
end


## The heart of it #################################################################################

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
- `E0::Array{Float64,2}`: Electric field before propagation
- `dz`: Distance propagated in z direction
- `diskR`: Radius of discs
- `eps`: Dielectric permittivity
- `tilt_x`: Disc tilt in x direction
- `tilt_y`: Disc tilt in y direction
- `surface`: Surface roughness of disc (only at end of propagation)
- `lambda`: Wavelength of electric field

See also: [`propagatorMomentumSpace`](@ref)
"""
function propagator(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda, coords::CoordinateSystem)
    #For real eps we need to make sure it's represented as real(eps)-0.0im so that sqrt(k0^2- Kx^2 - Ky^2)[propagatorNoTilts] has negative imaginary part for
    #Kx^2 + Ky^2>k0 and these modes propagate with a loss:
    imag(eps) == 0.0 ? eps = real(eps) - 0.0im : nothing
    k0 = 2*pi/lambda*sqrt(eps)
    # Call the propagator and add a phase imposed by the tilt
    E0 = propagatorNoTilts(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda, coords)
    # Tilts:
    E0 .*= [exp(-1im*k0*tilt_x*x) * exp(-1im*k0*tilt_y*y) for x in coords.X, y in coords.Y]
    # More general: Any surface misalignments:
    E0 .*= exp.(-1im*k0*surface) #(the element wise (exp.) is important, otherwise "surface" is treated as a matrix!)
    return E0
end

"""
    propagatorNoTilts(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda)

Wrapped by [`propagator`](@ref). Go there for documentation. Tilt arguments to be
compatible with other propagators.
"""
function propagatorNoTilts(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda, coords::CoordinateSystem)
    # Diffract at the Disk. Only the disk is diffracting.
    E0 .*= [abs(x^2 + y^2) < diskR^2 for x in coords.X, y in coords.Y]
    # FFT the E-Field to spatial frequencies
    # fft! and ifft! in the current release (1.2.2) only work with type ComplexF32 and ComplexF64
    # fft and ifft seem more stable
    FFTW.fft!(E0)
    E0 = FFTW.fftshift(E0)

    # TODO: If maximum k is higher than k0, then it is not defined
    #       what happens with this mode
    #       We should give a warning and handle this here
    #       At the moment the script will just also propagate with a loss for those components
    # Propagate through space
    k0 = 2*pi/lambda*sqrt(eps)
    k_prop = [sqrt(k0^2 - Kx^2 - Ky^2) for Kx in coords.kX, Ky in coords.kY]

    E0 = E0 .* exp.(-1im*k_prop*dz)
    # Backtransform
    E0 = FFTW.ifftshift(E0)
    FFTW.ifft!(E0)
    return E0
end

"""
    propagatorMomentumSpace(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda)

Propagator that assumes E0 is already in momentum space. Mix between [`propagator`](@ref)
and [`propagatorNoTilts`](@ref). Go to [`propagator`](@ref) for documentation.
"""
function propagatorMomentumSpace(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda, coords::CoordinateSystem)
    # Propagate through space
    k0 = 2*pi/lambda*sqrt(eps)
    k_prop = [sqrt(k0^2 - Kx^2 - Ky^2) for Kx in coords.kX, Ky in coords.kY]
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
    propagator1D(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda)

This propagator just does the phase propagation. Go to [`propagator`](@ref)
for documentation. 3D arguments to be compatible with other propagators.
"""
function propagator1D(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda, coords::CoordinateSystem)
    # Version of the propagator without the fft
    # should be faster and easy to check consistency with 1D calc

    # Propagate through space
    k0 = 2*pi/lambda*sqrt(eps)
    e1 = E0.*exp(-1im*k0*dz)
    return e1

end
