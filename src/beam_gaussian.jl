#
#   Gaussian Beam
#
#   V: 2019-04-05
#
#   Stefan Knirck
#


export gauss_profile

"""
    gauss_profile(X = -0.5:0.01:0.5, Y = -0.5:0.01:0.5, z = 0.000000001, omega0 = 0.1, f = 10e9)

Compute normalized gaussian beam profile.

# Arguments
- `X::AbstractRange{Float}`: x coordinates
- `Y::AbstractRange{Float}`: y coordinates
- `z::Float`: Distance from beam waist
- `omega0::Float`: Beam waist radius
- `f::Float`: Frequency of monochromatic EM radiation
"""
function gauss_profile(;X = -0.5:0.01:0.5, Y = -0.5:0.01:0.5, z = 0.000000001, omega0 = 0.1, f = 10e9)
#TODO: Documentation in Jupyter notenbook doesnt match implementation (which is equal to wikipedia).
	radius = [sqrt(x^2 + y^2) for x in X, y in Y]

    c = 299792458.
    lambda = c/f
    #Gaussian Beam profile:
    k = 2*pi/lambda
    zR = pi*omega0^2/lambda
    omega(z) = omega0*sqrt(1+(z/zR)^2)
    R(z)     = z*(1+(zR/z)^2)
    Psi(z)   = atan(z/zR)

   return omega0/omega(z) .* exp.(-radius.^2/omega(z)^2) .* exp.(-1im .* (k*z .+ k*radius.^2. /(2*R.(z)) .- Psi.(z)))
end
