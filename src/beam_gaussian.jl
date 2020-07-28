#
#   Gaussian Beam
#
#   V: 2019-04-05
#
#   Stefan Knirck
#


export gauss_profile

@doc raw"""
    gauss_profile(X = -0.5:0.01:0.5, Y = -0.5:0.01:0.5, z = 0.000000001, omega0 = 0.1, f = 10e9)

Compute normalized 2D gaussian beam profile at given z position (waist at ```z = 0```).

# Arguments
- `X::AbstractRange{Float}`: x coordinates
- `Y::AbstractRange{Float}`: y coordinates
- `z`: Distance from beam waist
- `omega0` ```> 0```: Beam waist radius
- `f` ```> 0```: Frequency of monochromatic EM radiation
"""
function gauss_profile(;X = -0.5:0.01:0.5, Y = -0.5:0.01:0.5, z = 0.000000001, omega0 = 0.1, f = 10e9)
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
