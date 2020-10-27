# 3D Algorithms
This file describes how to use the different implemented 3D algorithms.

## Cheerleader and Dancer
Iteratively propagate fields between boundaries for `nmax` iterations or until output amplitude is smaller than `amin`.

Arguments:
* `amin`:           Mimum (local) amplitude of a field, in order to be propagated
* `nmax`:           Maximum number of beam iteration steps, directly equivalent to how many boundaries a beam 'sees' at most
* `bdry`:           SetupBoundaries object, containing all relevant geometrical information (disk positions, epsilon, etc).
* `coords`:         CoordinateSystem object containing X, Y coordinates in real and K space
* `f`:              Frequency in Hz
* `prop`:           Propagator function to use. Standard is `propagator()`
* `reflect`:        If nothing (standar value), the axion-induced signal is computed.
                    If set, this field defines a beam, for which the reflected beam will be calculated
* `emit`:           If nothing (standar value), the axion-induced signal is computed.
                    If set, this array contains all emitted fields
* `returnsum`:      (only `dancer`) If false, the out-propagating contributions after each iteration will be returned, without summing.
* `immediatesum`:   (only `dancer`) If false, the out-propagating contributions will be saved and summed up at the end.
* `returnboth`:     (only `cheerleader`) If true, the out-propagating waves on the right and on the left are returned. Otherwise only the out-propagating wave on the right. Can be used to calculate transmission.

Returns:

A complex array `emitted_fields` of dimension `(length(X), length(Y))` containing the emitted beam shape.
Normalization: For the boost factor, the fields are normalized to the axion-induced field `E0`, i.e. the total power contained in the beam is
```julia
  sum(abs2.(emitted_fields * E0)*dx*dx)
```
A perfect mirror would emit
```julia
  abs2.(E0)*pi*diskR^2
```
Therefore, the boost factor is
```julia
  sum(abs2.(emitted_fields * E0)*dx*dx) / (pi*diskR^2)
```

### Initialization
#### `cheerleader()`
`cheerleader` can calculate the fields that leave the system on both sides, and can also calculate transmission.

Example for `cheerleader()`:
```julia
# Coordinate System
### Test CHEERLEADER ###
@everywhere begin
    # Initialize Coordinate System
    dx = 0.01
    coords = SeedCoordinateSystem(X = -0.5:dx:0.5, Y = -0.5:dx:0.5)
    
    diskR = 0.15
    
    # SetupBoundaries
    # No Tilt or surface roughness here.
    epsilon = 9
    eps = Array{Complex{Float64}}([NaN,1,epsilon,1])
    R = -(1 - sqrt(epsilon)) / (1 + sqrt(epsilon))
    r = Array{Complex{Float64}}([1.0, -R, R]);
    distance = [0.0, 0.15, 0.15/sqrt.(epsilon), 0.0]*1e-2
    
    sbdry = SeedSetupBoundaries(coords, diskno=1, distance=distance, epsilon=eps)
    
end
```

#### `dancer()`
**This algorithm is the one from the previous version. `cheerleader` needs only half as much iterations.**

`dancer` only calculates fields that leave the system to the right, and can only calculate reflectivities. Therefore, the mirror is not a region and not initialized as such.
`dancer` now calculates reflectivities from permittivities.

Example for `dancer()`:
```julia
# Coordinate System
@everywhere begin
    # Initialize Coordinate System
    dx = 0.01
    coords = SeedCoordinateSystem(X = -0.5:dx:0.5, Y = -0.5:dx:0.5)
    
    diskR = 0.15
    
    # SetupBoundaries
    # No Tilt or surface roughness here.
    epsilon = 9
    eps = Array{Complex{Float64}}([NaN,1,epsilon,1])
    R = -(1 - sqrt(epsilon)) / (1 + sqrt(epsilon))
    r = Array{Complex{Float64}}([1.0, -R, R]);
    distance = [0.0, 0.15, 0.15/sqrt.(epsilon), 0.0]*1e-2
    
    sbdry = SeedSetupBoundaries(coords, diskno=1, distance=distance, epsilon=eps)
    
end
```

### Convergence Remarks
Until now I always used `amin = 0`, but a different one might be also useful, but increase the runtime a bit although I did not study this due to time. For one disk + mirror I recommend ``nmax = 100`` (foŕ `dancer()` double those numbers, but don't use dancer anymore...).
Empirically for more disks: for 5 Sapphire disks ``nmax=900`` seems still good enough. For higher number it should roughly scale quadratically with the number of disks, so ``nmax=12800`` is good for 20 disks.

The X and Y grid should be set in such a way that the resolution is at least half a wavelength, i.e. for this example it should be sufficient to just use ``X=-0.5:0.01:0.5,Y=-0.5:0.01:0.5``, since the wavelength at 10GHz is roughly 3cm.



## Transformer
Apply transfer matrices to match modes between different regions.
The transfer matrices for reflection/transmission are directly calculated using the permittivities. Therefore, there is no need to set the reflection coefficients anymore. Note: epsilon = NaN in initialization to match behavior of cheerleader. However epsilon is still set in the code, so if you want to consider something other than a mirror, you can.

Arguments:
* `bdry`:           SetupBoundaries object, containing all relevant geometrical information (disk positions, epsilon, etc).
* `coords`:         CoordinateSystem object containing X, Y coordinates in real and K space
* `modes`:       Modes object containing important information for modes
* `f`:              Frequency in Hz
* `prop`:           Propagator function to use. Standard is `propagator()`
* `reflect`:        If nothing (standar value), the axion-induced signal is computed.
                    If set, this field defines a beam, for which the reflected beam will be calculated.
                    A vector with the mode coefficients of the modes. To compute it from an arbitrary field distribution use `field2modes()`.
* `emit`:           If nothing (standar value), the axion-induced signal is computed.
                    If set, this array contains all emitted fields.
                    A vector with the mode coefficients of the modes. To compute it from an arbitrary field distribution use `field2modes()`.

Returns:

Vector with the mode coefficients of the boosted fields. To compute the field distribution from this use `field2modes()`.
Normalization: For the boost factor the vector is normalized such that its absolute `abs2.(vec) == 1` for a mirror, i.e. its absolute directly gives the boost factor. Note that this is different from `cheerleader` which returns fields normalized to the axion-induced field `E0` (see above).


### Initialization
```julia
@everywhere begin
    # Initialize Coordinate System
    dx = 0.01
    coords = SeedCoordinateSystem(X = -0.5:dx:0.5, Y = -0.5:dx:0.5)
    
    diskR = 0.15
    
    # SetupBoundaries
    # No Tilt or surface roughness here.
    epsilon = 9
    eps = Array{Complex{Float64}}([NaN,1,epsilon,1])
    R = -(1 - sqrt(epsilon)) / (1 + sqrt(epsilon))
    r = Array{Complex{Float64}}([1.0, -R, R]);
    distance = [0.0, 0.15, 0.15/sqrt.(epsilon), 0.0]*1e-2
    
    sbdry = SeedSetupBoundaries(coords, diskno=1, distance=distance, epsilon=eps)
    
    # Initialize modes
    Mmax = 10
    Lmax = 0 # l-Modes are irrelevant for the azimuthally symmetric haloscope
    # For a 1D calculation:
    #modes = SeedModes(coords, ThreeDim=false, Mmax=Mmax, Lmax=Lmax, diskR=diskR)
    # For 3D:
    modes = SeedModes(coords, ThreeDim=true, Mmax=Mmax, Lmax=Lmax, diskR=diskR)
    
end
```

### Convergence Remarks
The transfer matrix system is solved numerically using simply Julia's backslash `\` operator. For the low number of modes considered here (<~ 100) this should be fairly accurate and give smaller numerical errors than the systematic error coming from neglecting higher modes.

There might be a numerical error in calculating the mixing matrices between the modes. It could be estimated by for example varying the grid size of the (X,Y) grid.

To estimate the error from neglecting higher modes, simply include more modes. For the emitted power of the disks one can say that for higher `m` the coupling of the axion-induced field and the mode gets smaller. Therefore, for the result to be valid is a neccessary condition, that including higher modes does not contribute significant power. This should be checked at least for one more higher mode for each calculation.
