#
# Transfomer Alogrithm for BoostFractor
#
# V: 2019-08-06
#
# Stefan Knirck
#
export transformer,init_waveguidemodes_3d,init_modes_1d,calc_propagation_matrices,field2modes,modes2field

# Transformation Matrices
using LinearAlgebra
# Eigenmodes (Bessel Functions):
using SpecialFunctions, FunctionZeros

#######################################################################################################
## Modes ##############################################################################################

# ---------------------- Initializations ------------------------------------------------

# Pre-calculate the mode patterns to speed up the matching calculations

# Number of modes to take into account
global M, L = 6,0
global mode_patterns = Array{Complex{Float64}}(zeros(M,2*L+1, length(X), length(Y)))
global mode_kt = Array{Complex{Float64}}(zeros(M,2*L+1))

global id = Matrix{Float64}(I, (M*(2L+1),M*(2L+1)))
global zeromatrix = zeros(((M*(2L+1),M*(2L+1))))
# Indexing function to get sub-matrices
global i(k) = ((k-1)*M*(2L+1)+1):((k)*M*(2L+1))

"""
Initialize waveguide-modes for a 3D calculation
"""
function init_waveguidemodes_3d(Mmax,Lmax;diskR=0.15)
    global M, L = Mmax,Lmax
    global mode_patterns = Array{Complex{Float64}}(zeros(M,2*L+1, length(X), length(Y)))
    global mode_kt = Array{Complex{Float64}}(zeros(M,2*L+1))
    #mode_patterns = Array{Complex{Float64}}(zeros(M,2*L+1, length(X), length(Y)))
    #mode_kt = Array{Complex{Float64}}(zeros(M,2*L+1))
    for m in 1:M, l in -L:L
        mode_kt[m,l+L+1], mode_patterns[m,l+L+1,:,:] = waveguidemode(m,l; X=X,Y=Y,diskR=diskR)
    end

    global id = Matrix{Float64}(I, (M*(2L+1),M*(2L+1)))
    global zeromatrix = zeros(((M*(2L+1),M*(2L+1))))

    # Indexing function to get sub-matrices
    #global i(k) = ((k-1)*M*(2L+1)+1):((k)*M*(2L+1))
end

"""
Initialize 1D modes for a 1D calculation.
"""
function init_modes_1d()
    global M = 1
    global L = 0
    global mode_patterns = Array{Complex{Float64}}(ones(1,1, length(X), length(Y)))
    global mode_kt = Array{Complex{Float64}}(zeros(1,1))


    global id = Matrix{Float64}(I, (M*(2L+1),M*(2L+1)))
    global zeromatrix = zeros(((M*(2L+1),M*(2L+1))))

    # Indexing function to get sub-matrices
    #global i(k) = ((k-1)*M*(2L+1)+1):((k)*M*(2L+1))
end


# ---------------------- Functionality ------------------------------------------------

"""
Calculate the transverse k and field distribution for a dielectric waveguidemode
"""
function waveguidemode(m,l; X=-0.5:0.007:0.5,  Y=-0.5:0.007:0.5,dx=0.007,dy=dx, diskR=0.15, k0=2pi/0.03)
    RR = [sqrt(x^2 + y^2) for x in X, y in Y]
    Phi = [atan(y,x) for x in X, y in Y]
    kr = besselj_zero(abs.(l),m)/diskR
    pattern = besselj.(l,kr.*RR).*exp.(-1im*l*Phi)
    pattern[RR .> diskR] .*= 0 # Cutoff
    ktransversal = kr

    pattern ./= sqrt(sum(abs2.(pattern))) #*dx*dy)
    return ktransversal, pattern
end

"""
Calculate the vector of mode-coefficients that describe the uniform axion-induced
field and that is normalized to power 1.
"""
function axion_induced_modes(;B=ones(length(X), length(Y)), velocity_x=0, diskR=0.15,f=20e9)

    # Inaccuracies of the emitted fields: BField and Velocity Effects ###################
    if velocity_x != 0
        B = Array{Complex{Float64}}(B)
        c = 299792458.
        Ma_PerMeter = 2pi*f/c # k = 2pi/lambda (c/f = lambda)
        B .*= [exp(-1im*Ma_PerMeter*(-velocity_x)*x) for x in X, y in Y]
    end

    # Only the axion-induced field on the disks matters:
    B .*= [sqrt(x^2 + y^2) <= diskR for x in X, y in Y]

    # Note: in all the normalizations the dx factors are dropped, since they should drop out in the final
    #       integrals
    B ./= sqrt.( sum(abs2.(B)) )

    modes_intital = Array{Complex{Float64}}(zeros(M*(2L+1)))
    for m in 1:M, l in -L:L
        modes_intital[(m-1)*(2L+1)+l+L+1] =
                sum( conj.(mode_patterns[m,l+L+1,:,:]) .* B )
    end

    return modes_intital
end

"""
Get the mode coefficients for a given field distribution E(x,y)
"""
field2modes(pattern;diskR=0.15) = axion_induced_modes(;B=pattern,diskR=diskR)

"""
Get the field distribution E(x,y) for a given vector of mode coefficients
"""
function modes2field(modes)
    result = Array{Complex{Float64}}(zeros(length(X), length(Y)))
    for m in 1:M, l in -L:L
        result .+= modes[(m-1)*(2L+1)+l+L+1].*mode_patterns[m,l+L+1,:,:]
    end
    return result
end


# Mode Mixing Matrix for the Propagation in a gap
function propagation_matrix(dz, diskR, eps, tilt_x, tilt_y, surface, lambda; is_air=(eps==1), onlydiagonal=false, prop=propagator)
    matching_matrix = Array{Complex{Float64}}(zeros(M*(2L+1),M*(2L+1)))

    k0 = 2pi/lambda*sqrt(eps)

    # Define the propagation function
    propfunc = nothing # initialize
    if is_air
        # In air use the propagator we get
        function propagate(x)
            return prop(copy(x), dz, diskR, eps, tilt_x, tilt_y, surface, lambda)
        end
        propfunc = propagate
    else
        # In the disk the modes are eigenmodes, so we only have to apply the
        # inaccuracies and can apply the propagation later seperately
        propfunc(efields) = efields.*[exp(-1im*k0*tilt_x*x) * exp(-1im*k0*tilt_y*y) for x in X, y in Y].*exp(-1im*surface)
    end

    # Calculate the mixing matrix
    for m_prime in 1:M, l_prime in -L:L
        # Propagated fields of mode (m_prime, l_prime)
        propagated = propfunc(mode_patterns[m_prime,l_prime+L+1,:,:])

        for m in (onlydiagonal ? [m_prime] : 1:M), l in (onlydiagonal ? [l_prime] : -L:L)

            # P_ml^m'l' = int dA E_ml* propagated(E_m'l')
            matching_matrix[(m-1)*(2L+1)+l+L+1, (m_prime-1)*(2L+1)+l_prime+L+1] =
                    sum( conj.(mode_patterns[m,l+L+1,:,:]) .* propagated ) #*dx*dy

            #v = 1-abs2.(matching_matrix[(m-1)*(2L+1)+l+L+1, (m_prime-1)*(2L+1)+l_prime+L+1])
        end
    end

    if !is_air
        propagation_matrix = Array{Complex{Float64}}(zeros(M*(2L+1),M*(2L+1)))
        #The propagation within the disk is still missing
        for m in 1:M, l in -L:L
            kz = sqrt(k0^2 - mode_kt[m,l+L+1])
            propagation_matrix[(m-1)*(2L+1)+l+L+1, (m-1)*(2L+1)+l+L+1] = exp(-1im*kz*dz)
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

function add_boundary(transm, n_left, n_right, diffprop)
    # We calculate G_r P_r analogous to eqs. 4.7
    # where n_right = n_{r+1}, n_left = n_r

    # The matrix encoding reflection and transmission
    G = (( (1. /(2*n_right)).*[(n_right+n_left)*id (n_right-n_left)*id ; (n_right-n_left)*id (n_right+n_left)*id] ))

    # The product, i.e. transfer matrix
    transm *= G * [diffprop zeromatrix; zeromatrix inv(Array{Complex{Float64}}(diffprop))]

    # Note: we build up the system from the end (Lm) downwards until L0
    # so this makes a transfer matrix from interface n -> m to a function that goes from interface n-1 ->m
    # This is convenient, because using this iteratively one arrives at exactly the T_n^m matrix from
    # the theoretical foundations paper

    return transm
end

function axion_contrib(T,n1,n0, initial)
    # Calculas one summand of the term (M[2,1]+M[1,1]) E_0 = \sum{s=1...m} (T_s^m[2,1]+T_s^m[1,1]) E_0
    # as in equation 4.14a
    return (1. /n1^2 - 1. /n0^2)/2 .* (T[i(2),i(1)]*(copy(initial)) + T[i(2),i(2)]*(copy(initial)))
end

"""
Pre-calculates all the propagation matrices.
Useful, if they should be altered later (e.g. take out modes, add some additional mixing, etc.)
"""
function calc_propagation_matrices(bdry::SetupBoundaries; f=10.0e9, prop=propagator, diskR=0.15)
    Nregions = length(bdry.eps)
    c = 299792458.
    wavelength = c / f
    return [ propagation_matrix(bdry.distance[i], diskR, bdry.eps[i],
            bdry.relative_tilt_x[i], bdry.relative_tilt_y[i], bdry.relative_surfaces[i,:,:], wavelength;
            prop=prop) for i in 1:(Nregions) ]
end

"""
Transformer Algorithm using Transfer Matrices and Modes to do the 3D Calculation.
"""
function transformer(bdry::SetupBoundaries; f=10.0e9, velocity_x=0, prop=propagator, propagation_matrices=nothing, Xset=X, Yset=Y, diskR=0.15, emit=axion_induced_modes(;B=ones(length(X),length(Y)),diskR=diskR), reflect=nothing)
    # For the transformer the region of the mirror must contain a high dielectric constant,
    # as the mirror is not explicitly taken into account

    #Definitions
    transmissionfunction_complete = [id zeromatrix ; zeromatrix id ]

    c = 299792458.
    wavelength = c / f

    initial = emit
    #println(initial)
    axion_beam = Array{Complex{Float64}}(zeros((M)*(2L+1)))
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
        axion_beam .+= axion_contrib(transmissionfunction_complete, sqrt(bdry.eps[idx_reg(s+1)]), sqrt(bdry.eps[idx_reg(s)]), initial)

        # calculate T_s^m ---------------------------

        # Propagation matrix (later become the subblocks of P)
        diffprop = (propagation_matrices == nothing ?
                        propagation_matrix(bdry.distance[idx_reg(s)], diskR, bdry.eps[idx_reg(s)], bdry.relative_tilt_x[idx_reg(s)], bdry.relative_tilt_y[idx_reg(s)], bdry.relative_surfaces[idx_reg(s),:,:], wavelength; prop=prop) :
                        propagation_matrices[idx_reg(s)])

        # T_s^m = T_{s+1}^m G_s P_s
        transmissionfunction_complete = add_boundary(transmissionfunction_complete,
                         sqrt(bdry.eps[idx_reg(s)]), sqrt(bdry.eps[idx_reg(s+1)]), diffprop)
    end

    # The rest of 4.14a
    boost =  - (transmissionfunction_complete[i(2),i(2)]) \ (axion_beam)
    # The backslash operator A\b solves the linear system Ax = b for x
    # Alaternative ways are e.g.
    #rtol = sqrt(eps(real(float(one(eltype(transmissionfunction_complete[i(2),i(2)]))))))
    #boost = - pinv(transmissionfunction_complete[i(2),i(2)], rtol=rtol) * (axion_beam)

    # If no reflectivity is ought to be calculated, we only return the axion field
    if reflect == nothing
        return boost
    end

    refl = transmissionfunction_complete[i(2),i(2)] \
           ((transmissionfunction_complete[i(2),i(1)]) * (reflect))
    return boost, refl
end
