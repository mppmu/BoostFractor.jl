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


export SetupBoundaries, SeedSetupBoundaries, propagator, propagator1D, propagatorNoTilts, CoordinateSystem, SeedCoordinateSystem, wavelength, propagators, init_prop



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
    kX = ifftshift(get_kspace_coords(X))
    kY = ifftshift(get_kspace_coords(Y)) #transform k-space coordinates to fft-ordering
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
    relative_surfaces::Array{Complex{Float64},3} # = [z, x, y]
    # etc.
end

struct Phase_shifts
    #pre-calculate k0 for each boundary.
    k0::Array{ComplexF64, 1}

    k_prop::Array{ComplexF64, 3}
    #combine surface misalignments and tilts into one array of phase-shifts
    surface::Array{ComplexF64, 3}
    #surface misalignments can usually be applied for only the area within the diskradius
    #for dancer and cheerleader: need to account for the first/last boundaries.
    #single 2D slice for dancer, two for cheerleader, unused for transformer
    surface_out::Array{ComplexF64,3}
    #two more values that are needed
    diskR::Float64
    dz::Array{Float64,1}
    #convenient arrays for slicing to get the subarray enclosing the disk area
    x_sub::Array{Int64, 1}
    y_sub::Array{Int64, 1}
end

#calcualte the phase shifts for every boundary
function init_prop(bdry::SetupBoundaries, coords::CoordinateSystem, lambda, radius; output_surfaces = [])
    k0 = 2 * pi * sqrt.(bdry.eps) / lambda

    #t_x = exp.(-1im * coords.X * transpose(k0 .* bdry.relative_tilt_x))
    #t_y = exp.(-1im * coords.Y * transpose(k0 .* bdry.relative_tilt_x))

    #surface_prop = [exp(-1im * k0[z] * (bdry.relative_tilt_y[z] * coords.Y[y] + bdry.relative_tilt_x[z] * coords.X[x] + bdry.relative_surfaces[z, x, y])) for x in 1:length(coords.X), y in 1:length(coords.Y), z in 1:length(bdry.distance)]

    #In most cases the field outside the disk is cut off right after applying surface misalignments,
    #because of this it is more efficient to apply surface misalignments only for the subarry containing the disk
    x_sub = coords.X[abs.(coords.X) .<= radius] #X / Y coordinates for the subarray
    y_sub = coords.Y[abs.(coords.Y) .<= radius]

    x_sub_min = (length(coords.X) - length(x_sub)) ÷ 2 + 1 #maximum and minimum indices of the subarry
    x_sub_max = x_sub_min + length(x_sub) - 1

    y_sub_min = (length(coords.Y) - length(y_sub)) ÷ 2 + 1
    y_sub_max = y_sub_min + length(y_sub) - 1

    x_sub_index = x_sub_min:1:x_sub_max # for slicing, field[x_sub_index, y_sub_index] will return the field
    y_sub_index = y_sub_min:1:y_sub_max # in the square subarray exactly enclosing the disk-area

    t_x = [exp(-1im * k0[z] * bdry.relative_tilt_x[z] * x) for x in x_sub, z in 1:length(k0)]
    t_y = [exp(-1im * k0[z] * bdry.relative_tilt_y[z] * y) for y in y_sub, z in 1:length(k0)]
    surface_phase = [exp(-1im * k0[z] * bdry.relative_surfaces[z, x, y]) for x in x_sub_index, y in y_sub_index, z in 1:length(k0)]
    surface_phase = [t_x[x,z] * t_y[y,z] * surface_phase[x,y,z] for x in 1:length(x_sub), y in 1:length(y_sub), z in 1:length(k0)]

    k_prop = [exp(-1im * bdry.distance[z] * conj(sqrt( ComplexF64(k0[z]^2 - kx^2 - ky^2)))) for kx in coords.kX, ky in coords.kY, z in 1:length(bdry.distance)]

    #for dancer and cheerleader we need the surface misalignments outside the disk area for the rightmost/leftmost boundaries
    surface_phase_out = Array{ComplexF64, 3}(undef, length(coords.X), length(coords.Y), length(output_surfaces))
    for i in 1:length(output_surfaces)
        z = output_surfaces[i]
        t_x = [exp(-1im * k0[z] * bdry.relative_tilt_x[z] * x) for x in coords.X]
        t_y = [exp(-1im * k0[z] * bdry.relative_tilt_y[z] * y) for y in coords.Y]
        surface_phase_out[:,:,i] = [exp(-1im * k0[z] * bdry.relative_surfaces[z, x, y]) * t_x[x] * t_y[y] for x in 1:length(coords.X), y in 1:length(coords.Y)]
        surface_phase_out[x_sub_index, y_sub_index, i] .= 1+0im # for the sub array the misalignments were already applied
    end
    return Phase_shifts(k0, k_prop, surface_phase, surface_phase_out, radius, bdry.distance, x_sub_index, y_sub_index)#, surface_prop, surface_out_prop, radius)
end


## Convenient tools ################################################################################
@doc raw"""
    SeedSetupBoundaries(diskno=3)

Initialize `mutable struct SetupBoundaries` with sensible values.

# Arguments
- `diskno::Int` ```> 0```: Number of dielectric discs
"""
function SeedSetupBoundaries(coords::CoordinateSystem; diskno=3, distance=nothing, epsilon=nothing, relative_tilt_x=zeros(2*diskno+2), relative_tilt_y=zeros(2*diskno+2), relative_surfaces=zeros(2*diskno+2 , length(coords.X), length(coords.Y)))

    # Initialize SetupBoundaries entries with default values given diskno. Rest was already defined in function definition.
    if distance === nothing # [m]
        distance = [0.0]
        append!(distance, [ x % 2 == 1 ? 8e-3 : 1e-3 for x in 1:2*(diskno) ])
        append!(distance, 0e-3) #8e-3)
    end

    if epsilon === nothing
        epsilon = [NaN]
        append!(epsilon, [ x % 2 == 1.0 ? 1.0 : 9.0 for x in 1:2*(diskno) ])
        append!(epsilon, 1.0)
    end

    reflectivities = complex([1.0])
    R = [(sqrt(epsilon[i-1]) - sqrt(epsilon[i])) / (sqrt(epsilon[i-1]) + sqrt(epsilon[i])) for i in 3:length(epsilon)]
    append!(reflectivities, R)

    # Check if initialization was self-consistent
    length(distance) == length(reflectivities)+1 == length(epsilon) == length(relative_tilt_x) == length(relative_tilt_y) == size(relative_surfaces, 1) || throw(DimensionMismatch("the arrays in your SetupBoundaries objects don't fit together!"))

    return SetupBoundaries(distance, Array{Complex{Float64}}(reflectivities), Array{Complex{Float64}}(epsilon), relative_tilt_x, relative_tilt_y, relative_surfaces)
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
function propagator(E0, i, coords::CoordinateSystem, prop::Phase_shifts, plan, i_plan)

    # Call the propagator and add a phase imposed by the tilt
    E0 = propagatorNoTilts(E0, i, coords, prop, plan, i_plan)

    # More general: Any surface misalignments:
    E0[prop.x_sub, prop.y_sub] .*= prop.surface[:,:,i] #only applied on a subarry enclosing the disk area
    return E0
end

"""
    propagatorNoTilts(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda)

Wrapped by [`propagator`](@ref). Go there for documentation. Tilt arguments to be
compatible with other propagators.
"""
function propagatorNoTilts(E0, i, coords::CoordinateSystem, prop::Phase_shifts, plan, i_plan)
    # Diffract at the Disk. Only the disk is diffracting.
    E0 .*= [abs(x^2 + y^2) < prop.diskR^2 for x in coords.X, y in coords.Y]
    # FFT the E-Field to spatial frequencies
    # fft! and ifft! in the current release (1.2.2) only work with type ComplexF32 and ComplexF64
    # fft and ifft seem more stable

    plan * E0

    # TODO: If maximum k is higher than k0, then it is not defined
    #       what happens with this mode
    #       We should give a warning and handle this here
    #       At the moment the script will just also propagate with a loss for those components
    # Propagate through space
    #k = [conj(sqrt( ComplexF64(prop.k0[i]^2 - kx^2 - ky^2))) for kx in coords.kX, ky in coords.kY]
    #E0 .*= exp.(-1im * prop.dz[i] * k)
    E0 .*= prop.k_prop[:,:,i]
    # Backtransform
    #E0 = FFTW.ifftshift(E0)
    i_plan * E0
    return E0
end

"""
    propagatorMomentumSpace(E0, i, coords::CoordinateSystem, prop::Phase_shifts, plan, i_plan)

Propagator that assumes E0 is already in momentum space. Mix between [`propagator`](@ref)
and [`propagatorNoTilts`](@ref). Go to [`propagator`](@ref) for documentation.
"""
function propagatorMomentumSpace(E0, i, coords::CoordinateSystem, prop::Phase_shifts, plan, i_plan)
    # Propagate through space
    E0 .* exp.(-1im*prop.k_prop[i]*prop.dz[i])

    # Transform to position space
    i_plan * E0

    # Diffract at the Disk. Only the disk is diffracting.
    E0 .*= [abs(x^2 + y^2) < prop.diskR^2 for x in coords.X, y in coords.Y]

    #tilts and surface misalignments, for now without check if neccessary
    view(E0, prop.x_sub, prop.y_sub) .*= prop.surface[:,:,i]

    # FFT the E-Field to spatial frequencies / momentum space
    plan * E0

    return E0
end

"""
    propagator1D(E0, dz, diskR, eps, tilt_x, tilt_y, surface, lambda)

This propagator just does the phase propagation. Go to [`propagator`](@ref)
for documentation. 3D arguments to be compatible with other propagators.
"""
function propagator1D(E0, i, coords::CoordinateSystem, prop::Phase_shifts, plan, i_plan)
    # Version of the propagator without the fft
    # should be faster and easy to check consistency with 1D calc

    # Propagate through space
    k_prop = conj(sqrt.(prop.k0[i]^2))
    E0 .*= exp(-1im*k_prop*dz)

end


#TODO: This function is no longer supported. reflectivity_transmissivity_1d does not exist; DiskDefinition has merged into SetuoBoundaries.
@doc raw"""
    initialize_reflection_transmission(freq::Float64, bdry::SetupBoundaries, disk::DiskDefiniton)

OUTDATED! Does not work, do not use!
Calculate reflection and transmission coefficients.

# Arguments
- `freq::Float64` ```> 0```: Frequency of EM radiation
- `bdry::SetupBoundaries`: Properties of dielectric boundaries
- `disk::DiskDefiniton`: Properties of dielectric discs
"""
# function initialize_reflection_transmission(freq::Float64, bdry::SetupBoundaries, coords::CoordinateSystem)#, disk::DiskDefiniton)
#     if disk === nothing
#         # Iniatilize reflection coefficients according to epsilon
#         r_left = ones(length(bdry.eps))
#         r_left[1] = -1
#         for i in 2:length(bdry.eps)
#             # The reflectivity at this point
#             r_left = (sqrt(bdry.eps[i-1])-sqrt(bdry.eps[i]))/(sqrt(bdry.eps[i])+sqrt(bdry.eps[i-1]))
#         end
#         r_right = -r_left
#         t_left = 1. + r_left
#         t_right = 1 .+ r_right
#     else
#         # Initailize reflection coefficients according to disk model
#         ref, trans = reflectivity_transmissivity_1d(freq, disk.thickness)
#         r_left = ones(length(bdry.eps),length(coords.kX),length(coords.kY)).*ref
#         r_right = r_left
#         t_left = ones(length(bdry.eps),length(coords.kX),length(coords.kY)).*trans
#         t_right = t_left
#     end
#     return r_left, r_right, t_left, t_right
# end
