#
# Transfomer Alogrithm for BoostFractor
#
# V: 2019-11-30
#
# Stefan Knirck
#
export transformer,calc_propagation_matrices,field2modes,modes2field, Modes, SeedModes

# Transformation Matrices
using LinearAlgebra
# Eigenmodes (Bessel Functions):
using SpecialFunctions, FunctionZeros

#######################################################################################################
## Modes ##############################################################################################

# ---------------------- Initializations ------------------------------------------------

"""
    Gives out E-field dimensions in a more readable way. Output should be either 1 or 3.
"""
e_field_dimensions(modes) = size(modes.mode_patterns)[5]

# Pre-calculate the mode patterns to speed up the matching calculations
"""
    Modes

Doc: TODO
"""
mutable struct Modes
    M::Int64
    L::Int64
    mode_patterns::Array{Complex{Float64}, 5}
    mode_kt::Array{Complex{Float64}, 2}
    id::Array{Complex{Float64}, 2}
    zeromatrix::Array{Complex{Float64}, 2}
end

"""
    SeedModes(coords::CoordinateSystem;ThreeDim=false, Mmax=1, Lmax=0, diskR=0.15)

Doc: TODO
"""

function SeedModes(coords::CoordinateSystem;ThreeDim=false, Mmax=1, Lmax=0, diskR=0.15,pattern_input=nothing)
    if ThreeDim
        M, L = Mmax,Lmax
        if pattern_input === nothing # If no specific pattern is put in, assume that E field is treated as scalar. Will change in some future version.
            mode_patterns = Array{Complex{Float64}}(zeros(M,2*L+1, length(coords.X), length(coords.Y),1))
        else
            mode_patterns = Array{Complex{Float64}}(zeros(M,2*L+1, length(coords.X), length(coords.Y),size(pattern_input)[end]))
        end
        mode_kt = Array{Complex{Float64}}(zeros(M,2*L+1))
        for m in 1:M, l in -L:L
            mode_kt[m,l+L+1], mode_patterns[m,l+L+1,:,:,:] = mode(m,l,L,coords;diskR=diskR,pattern_input=pattern_input)
        end
    else
        M = 1
        L = 0
        mode_patterns = Array{Complex{Float64}}(ones(M,2*L+1, length(coords.X), length(coords.Y),1))
        mode_kt = Array{Complex{Float64}}(zeros(M,2*L+1))
    end

    id = Matrix{Float64}(I, (M*(2L+1),M*(2L+1))) # I: Identity matrix
    zeromatrix = zeros(M*(2L+1),M*(2L+1))

    return Modes(M,L,mode_patterns,mode_kt,id,zeromatrix)
end

"""
    index(modes::Modes, k::Int64)

Indexing function to get sub-matrices. (Compare Knirck 6.18.)
"""
function index(modes::Modes, k::Int64)
    return ((k-1)*modes.M*(2modes.L+1)+1):(k*modes.M*(2modes.L+1))
end

# ---------------------- Functionality ------------------------------------------------

"""
    mode(m,l, coords::CoordinateSystem; diskR=0.15, k0=2pi/0.03)

Calculate the transverse k and field distribution for a mode.
Implements Knirck 6.12.
"""
function mode(m,l,L, coords::CoordinateSystem; diskR=0.15, k0=2pi/0.03,pattern_input=nothing)
    RR = [sqrt(x^2 + y^2) for x in coords.X, y in coords.Y]
    Phi = [atan(y,x) for x in coords.X, y in coords.Y]
    kr = besselj_zero(abs.(l),m)/diskR
    ktransversal = kr
    if pattern_input === nothing
        pattern = besselj.(l,kr.*RR).*exp.(-1im*l*Phi)
        pattern[RR .> diskR] .*= 0 # Cutoff
        pattern ./= sqrt(sum(abs2.(pattern))) #*dx*dy)
        pattern = reshape(pattern, (size(pattern)...,1)) # Expand dims without putting additional values in
    else
        pattern = pattern_input[m,l+L+1,:,:,:]
    end
    return ktransversal, pattern
end

"""
Calculate the vector of mode-coefficients that describe the uniform axion-induced
field and that is normalized to power 1.

Implements Knirck 6.17.
"""
function axion_induced_modes(coords::CoordinateSystem, modes::Modes;B=nothing, velocity_x=0, diskR=0.15,f=20e9)

    if B === nothing
        if e_field_dimensions(modes) == 1
            B = ones(length(coords.X), length(coords.Y), e_field_dimensions(modes))
        elseif e_field_dimensions(modes) == 3
            B = zeros(length(coords.X), length(coords.Y), e_field_dimensions(modes))
            B[:,:,2] = ones(length(coords.X), length(coords.Y))
        end
    end

    # Inaccuracies of the emitted fields: BField and Velocity Effects ###################
    if velocity_x != 0
        B = Array{Complex{Float64}}(B)
        lambda = wavelength(f)
        Ma_PerMeter = 2pi/lambda # k = 2pi/lambda (c/f = lambda)
        B .*= [exp(-1im*Ma_PerMeter*(-velocity_x)*x) for x in coords.X, y in coords.Y]
    end

    # Only the axion-induced field on the disks matters:
    B .*= [sqrt(x^2 + y^2) <= diskR for x in coords.X, y in coords.Y]

    # Note: in all the normalizations the dx factors are dropped, since they should drop out in the final
    #       integrals
    B ./= sqrt.( sum(abs2.(B)) )


    modes_initial = Array{Complex{Float64}}(zeros(modes.M*(2modes.L+1)))
    for m in 1:modes.M, l in -modes.L:modes.L
        # (m-1)*(2modes.L+1)+l+modes.L+1 walks over all possible m,l combinations
        modes_initial[(m-1)*(2modes.L+1)+l+modes.L+1] =
                sum( conj.(modes.mode_patterns[m,l+modes.L+1,:,:,:]) .* B )
    end


    return modes_initial
end

"""
Get the mode coefficients for a given field distribution E(x,y)
"""
field2modes(pattern, coords::CoordinateSystem, modes::Modes;diskR=0.15) = axion_induced_modes(coords, modes;B=pattern,diskR=diskR)

"""
Get the field distribution E(x,y) for a given vector of mode coefficients
"""
function modes2field(mode_coeffs, coords::CoordinateSystem, modes::Modes)
    result = Array{Complex{Float64}}(zeros(length(coords.X), length(coords.Y), e_field_dimensions(modes)))
    for m in 1:modes.M, l in -modes.L:modes.L
            result[:,:,:] .+= mode_coeffs[(m-1)*(2modes.L+1)+l+modes.L+1].*modes.mode_patterns[m,l+modes.L+1,:,:,:]
    end
    return result
end


# Mode Mixing Matrix for the Propagation in a gap
#TODO: tilts, eps and surface should come from SetupBoundaries object?
"""
"""
function propagation_matrix(dz, diskR, eps, tilt_x, tilt_y, surface, lambda, coords::CoordinateSystem, modes::Modes; is_air=(eps==1), onlydiagonal=false, prop=propagator)
    matching_matrix = Array{Complex{Float64}}(zeros(modes.M*(2modes.L+1),modes.M*(2modes.L+1)))

    k0 = 2pi/lambda*sqrt(eps)

    # Define the propagation function
    propfunc = nothing # initialize
    if is_air
        # In air use the propagator we get
        function propagate(x)
            return prop(copy(x), dz, diskR, eps, tilt_x, tilt_y, surface, lambda, coords)
        end
        propfunc = propagate
    else
        # In the disk the modes are eigenmodes, so we only have to apply the
        # inaccuracies and can apply the propagation later separately
        propfunc(efields) = efields.*[exp(-1im*k0*tilt_x*x) * exp(-1im*k0*tilt_y*y) for x in coords.X, y in coords.Y].*exp.(-1im*k0*surface)
        # Applying exp element-wise to the surface is very important otherwise it is e^M with M matrix
    end

    # Calculate the mixing matrix
    for m_prime in 1:modes.M, l_prime in -modes.L:modes.L
        # Propagated fields of mode (m_prime, l_prime)
        for i in 1:e_field_dimensions(modes)
            # If 3D E-field, fields propagate separately. For interaction need to implement in propagator.
            propagated = propfunc(modes.mode_patterns[m_prime,l_prime+modes.L+1,:,:,i])

            for m in (onlydiagonal ? [m_prime] : 1:modes.M), l in (onlydiagonal ? [l_prime] : -modes.L:modes.L)
                # P_ml^m'l' = <E_ml | E^p_m'l' > = ∫ dA \bm{E}_ml^* ⋅ \bm{E}^p_m'l' = ∑_{j = x,y,z} ∫ dA ({E}_j)_ml^* ⋅ ({E}_j)^p_m'l'
                # 6.15 in Knirck
                matching_matrix[(m-1)*(2modes.L+1)+l+modes.L+1, (m_prime-1)*(2modes.L+1)+l_prime+modes.L+1] +=
                        sum( conj.(modes.mode_patterns[m,l+modes.L+1,:,:,i]) .* propagated ) #*dx*dy

                #v = 1-abs2.(matching_matrix[(m-1)*(2L+1)+l+L+1, (m_prime-1)*(2L+1)+l_prime+L+1])
            end
        end
    end

    if !is_air
        propagation_matrix = Array{Complex{Float64}}(zeros(modes.M*(2modes.L+1),modes.M*(2modes.L+1)))
        #The propagation within the disk is still missing
        for m in 1:modes.M, l in -modes.L:modes.L
            kz = sqrt(k0^2 - modes.mode_kt[m,l+modes.L+1])
            propagation_matrix[(m-1)*(2modes.L+1)+l+modes.L+1, (m-1)*(2modes.L+1)+l+modes.L+1] = exp(-1im*kz*dz)
        end
        # It is important to note the multiplication from the left
        matching_matrix = propagation_matrix*matching_matrix
    end

    return matching_matrix

    # TODO: This only takes surface roughness at the end of the propagation into account, not at its
    # start. General problem in all the codes so far.
end


#########################################################################################################
## The heart of the transformation code #################################################################

# Comments refer to the formuli in the theoretical foundations paper (arXiv:1612.07057v2)
# but generalized to 3D making L and R vectors.
"""
Calculate Transfer matrix like in 4.9.
"""
function get_boundary_matrix(n_left, n_right, diffprop, modes::Modes)
    # We calculate G_r P_r analogous to eqs. 4.7
    # where n_right = n_{r+1}, n_left = n_r

    # The matrix encoding reflection and transmission (Knirck 6.18)
    G = (( (1. /(2*n_right)).*[(n_right+n_left)*modes.id (n_right-n_left)*modes.id ; (n_right-n_left)*modes.id (n_right+n_left)*modes.id] ))


    # The product, i.e. transfer matrix
    return G * [diffprop modes.zeromatrix; modes.zeromatrix inv(Array{Complex{Float64}}(diffprop))]

    # Note: we build up the system from the end (Lm) downwards until L0
    # so this makes a transfer matrix from interface n -> m to a function that goes from interface n-1 ->m
    # This is convenient, because using this iteratively one arrives at exactly the T_n^m matrix from
    # the theoretical foundations paper
end

"""
Calculates one summand of the term (M[2,1]+M[1,1]) E_0 = sum{s=1...m} (T_s^m[2,1]+T_s^m[1,1]) E_0
as in equation 4.14a
"""
function axion_contrib(T,n1,n0, initial, modes::Modes)
    axion_beam = (1. /n1^2 - 1. /n0^2)/2 .* (T[index(modes,2),index(modes,1)]*(copy(initial)) + T[index(modes,2),index(modes,2)]*(copy(initial)))
    return axion_beam
end

"""
Pre-calculates all the propagation matrices.
Useful, if they should be altered later (e.g. take out modes, add some additional mixing, etc.)
"""
function calc_propagation_matrices(bdry::SetupBoundaries, coords::CoordinateSystem, modes::Modes; f=10.0e9, prop=propagator, diskR=0.15)
    Nregions = length(bdry.eps)
    lambda = wavelength(f)
    return [ propagation_matrix(bdry.distance[i], diskR, bdry.eps[i],
            bdry.relative_tilt_x[i], bdry.relative_tilt_y[i], bdry.relative_surfaces[i,:,:], lambda, coords, modes;
            prop=prop) for i in 1:(Nregions) ]
end

"""
Transformer Algorithm using Transfer Matrices and Modes to do the 3D Calculation.
"""
function transformer(bdry::SetupBoundaries, coords::CoordinateSystem, modes::Modes; f=10.0e9, velocity_x=0, prop=propagator, propagation_matrices=nothing, diskR=0.15, emit=axion_induced_modes(coords,modes;B=nothing,velocity_x=velocity_x,diskR=diskR), reflect=nothing)
    # For the transformer the region of the mirror must contain a high dielectric constant,
    # as the mirror is not explicitly taken into account
    # To have same SetupBoundaries object for all codes and cheerleader assumes NaN, just define a high constant
    bdry.eps[isnan.(bdry.eps)] .= 1e30

    #Definitions
    transmissionfunction_complete = [modes.id modes.zeromatrix ; modes.zeromatrix modes.id ]
    lambda = wavelength(f)

    initial = emit
    #println(initial)
    axion_beam = Array{Complex{Float64}}(zeros((modes.M)*(2modes.L+1)))
    #println(axion_beam)

    #=
        To have the different algorithms consistent with each other,
        we always assume that the mirror is at the left and we want to calculated
        the boost factor / reflectivity /  ... on the right.
        For the transfer matrices it is however more convenient to calcluate them
        boost factor on the left (cf eq. 4.14a). Therefore, we use the following
        little dummy function to "reindex" the SetupBoundaries arrays such that
        we see it ordered the other way round.
    =#
    Nregions = length(bdry.eps)
    idx_reg(s) = Nregions-s+1

    #=
        Essentially we want to calculate 4.14a
        We iteratively calculate the T matrices, in each step directly computing its axion contribution.
        We calculate T_{m-1}^m first, and then expand iteratively to get T_{n}^m.
        Notice from eq. (4.9) that going one region down is a multiplication from the right, not the left, e.g.
        T_3^5 = T_4^5 G_3 P_3.

        I agree, iterating in reverse order and changing before the indices to
        reverse, is not neccessary. Though, like this "s" follows the same Indexing
        than in the theoretical foundations paper and we can call the function with
        a bdry structure having the mirror at the lowest index.
    =#
    for s in (Nregions-1):-1:1
        # Add up the summands of (M[2,1]+M[1,1]) E_0
        # (M is a sum over T_{s+1}^m S_s from s=1 to m) and we have just calculated
        #  T_{s+1}^m in the previous iteration)
        axion_beam .+= axion_contrib(transmissionfunction_complete, sqrt(bdry.eps[idx_reg(s+1)]), sqrt(bdry.eps[idx_reg(s)]), initial, modes)

        # calculate T_s^m ---------------------------

        # Propagation matrix (later become the subblocks of P)
        diffprop = (propagation_matrices === nothing ?
                        propagation_matrix(bdry.distance[idx_reg(s)], diskR, bdry.eps[idx_reg(s)], bdry.relative_tilt_x[idx_reg(s)], bdry.relative_tilt_y[idx_reg(s)], bdry.relative_surfaces[idx_reg(s),:,:], lambda, coords, modes; prop=prop) :
                        propagation_matrices[idx_reg(s)])

        # T_s^m = T_{s+1}^m G_s P_s
        transmissionfunction_complete *= get_boundary_matrix(sqrt(bdry.eps[idx_reg(s)]), sqrt(bdry.eps[idx_reg(s+1)]), diffprop, modes)
    end

    # The rest of 4.14a
    boost =  - (transmissionfunction_complete[index(modes,2),index(modes,2)]) \ (axion_beam)
    # The backslash operator A\b solves the linear system Ax = b for x
    # Alternative ways are e.g.
    #rtol = sqrt(eps(real(float(one(eltype(transmissionfunction_complete[index(modes,2),index(modes,2)]))))))
    #boost = - pinv(transmissionfunction_complete[index(modes,2),index(modes,2)], rtol=rtol) * (axion_beam)

    # If no reflectivity is ought to be calculated, we only return the axion field
    if reflect === nothing
        return boost
    end

    refl = transmissionfunction_complete[index(modes,2),index(modes,2)] \
           ((transmissionfunction_complete[index(modes,2),index(modes,1)]) * (reflect))
    return boost, refl
end


"""
Traces back the field of a solution for the reflectivity in the system.
Same arguments as transformer(), but first two arguments are the reflected beam
solution and the input beam. Returns an array with the left- and rightgoing field
amplitudes in each region.
"""
function transformer_trace_back(refleced_beam, input_beam,
    bdry::SetupBoundaries, coords::CoordinateSystem, modes::Modes;
    f=10.0e9, velocity_x=0, prop=propagator, propagation_matrices=nothing, diskR=0.15,
    emit=axion_induced_modes(coords,modes;B=ones(length(coords.X),length(coords.Y)),diskR=diskR))
    # Note: So far this is only implemented for the reflectivity. Axion-contributions are neglected!

    #Definitions
    transmissionfunction_complete = [modes.id modes.zeromatrix ; modes.zeromatrix modes.id ]

    c = 299792458.
    lambda = c / f


    Nregions = length(bdry.eps)
    idx_reg(s) = Nregions-s+1

    # The thing which we want to give back
    fields_regions = Array{Complex{Float64}}(zeros(Nregions,2,(modes.M)*(2modes.L+1)))

    #=
        We start with the leftmost region where the solution is known and successively transform through the system:
    =#
    solution_current = permutedims(hcat(input_beam, refleced_beam), [2,1])
    fields_regions[idx_reg(1),:,:] = copy(solution_current)
    for s in 1:Nregions-1 # (Nregions-1):-1:1
        # Propagation matrix (later become the subblocks of P)
        diffprop = (propagation_matrices === nothing ?
                        propagation_matrix(bdry.distance[idx_reg(s)], diskR, bdry.eps[idx_reg(s)], bdry.relative_tilt_x[idx_reg(s)], bdry.relative_tilt_y[idx_reg(s)], bdry.relative_surfaces[idx_reg(s),:,:], lambda, coords, modes; prop=prop) :
                        propagation_matrices[idx_reg(s)])

        # G_s P_s
        transmissionfunction_bdry = get_boundary_matrix(sqrt(bdry.eps[idx_reg(s)]), sqrt(bdry.eps[idx_reg(s+1)]), diffprop, modes)

        # Apply the transmission function to the solution
        # (R_s+1, L_s+1) = G_s P_s (R_s, L_s)
        solution_current = transmissionfunction_bdry * solution_current

        fields_regions[idx_reg(s+1),:,:] = copy(solution_current)
    end

    return fields_regions
end
