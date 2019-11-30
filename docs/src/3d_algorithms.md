# 3D Algorithms
This little docfile describes the different algorithms and (in particular the differences) how to initilize them with a minimal example.

## Dancer and Cheerleader
Iteratively propagate fields between boundaries for `nmax` iterations or until output amplitude is smaller than `amin`.

Arguments
* `amin`:           Mimum (local) amplitude of a field, in order to be propagated
* `nmax`:           Maximum number of beam iteration steps, directly equivalent to how many boundaries a beam 'sees' at most
* `bdry`:           SetupBoundaries object, containing all relevant geometrical information (disk positions, epsilon, etc). **Each algorithm expects a slightly differently initialized SetupBoundaries object. See below for details.**
* `f`:              Frequency in Hz
* `prop`:           Propagator function to use. Standard is `propagator()`
* `reflect`:        If nothing (standar value), the axion-induced signal is computed.
                    If set, this field defines a beam, for which the reflected beam will be calculated
* `emit`:           If nothing (standar value), the axion-induced signal is computed.
                    If set, this array contains all emitted fields
* `Xset`, `Yset`:   Explicitly set the coordinate system for the fields
* `returnsum`:      (only `dancer`) If false, the out-propagating contributions after each iteration will be returned, without summing.
* `immediatesum`:   (only `dancer`) If false, the out-propagating contributions will be saved and summed up at the end.
* `returnboth`:     (only `cheerleader`) If true, the out-propagating waves on the right and on the left are returned. Otherwise only the out-propagating wave on the right.

### Initilization for `dancer()`
**This algorithm is the one provided in the previous version. `cheerleader` is by a factor 2 faster and can calculate transmission.**

`dancer` only calculates fields that leave the system to the right, and can only calculate reflectivities. Therefore, the mirror is not a region and not initialized as such.
`dancer` needs to know the reflectivities between the regions and does not calculte them from the permittivities. This might be useful, if one wants to include an arbitrary surface loss which is independent of `epsilon`.

Example for `dancer()`:
```julia
# Coordinate System
X = [0.001]
Y = [0.001]
@everywhere sbdry = SeedSetupBoundaries(1)

# Permittivity
epsilon = 9
sbdry.eps = Array{Complex{Float64}}([1,epsilon,1]) #<-- The mirror is not included as a region

# Refelectivity
R = -(1 - sqrt(epsilon)) / (1 + sqrt(epsilon))
sbdry.r = Array{Complex{Float64}}([1.0, -R, R, 0]);
# The order is: Ref. Coeff as seen from [inside the mirror to vacuum, vacuum to disk, disk to vacuum, vacuum to receiver (should be 0)]

# We can also set relative tilts, etc.
sbdry.relative_tilt_x = zeros(3);
sbdry.relative_tilt_y = zeros(3);


# Surface Roughness
sbdry.relative_surfaces = zeros(3,length(X),length(Y));


# Initialize Coordiates
init_coords(X,Y)
```

### Initilization for `cheerleader()`
`cheerleader` can calculate the fields that leave the system on both sides, and can also calculate transmission.

Example for `cheerleader()`:
```julia
# Coordinate System
X = [0.001]
Y = [0.001]


@everywhere bdry = SeedSetupBoundaries(1)


# Permittivity
epsilon = 9
sbdry.eps = [NaN,1.0,epsilon, 1.0] #<-- Note that the mirror is included with a NaN, but we have 4 regions
# Reflectivity
R = -(1 - sqrt(epsilon)) / (1 + sqrt(epsilon))
sbdry.r = Array{Complex{Float64}}([1.0, -R, R]); #<-- There are only 3 boundaries between regions, so 3 reflectivities

# Tilts
sbdry.relative_tilt_x = zeros(4);
sbdry.relative_tilt_y = zeros(4);

# Surface Roughness
sbdry.relative_surfaces = zeros(4,length(X),length(Y));

# Initialize Coordiates
init_coords(X,Y)
```


## Transformer
Apply transfer matrices to match modes between different regions.
The transfer matrices for reflection/transmission are dierectly calculated using the permittivities. Therefore, there is no need to set the reflection coefficients anymore.

Arguments:
* `bdry`:           SetupBoundaries object, containing all relevant geometrical information (disk positions, epsilon, etc). **Each algorithm expects a slightly differently initialized SetupBoundaries object. See below for details.**
* `f`:              Frequency in Hz
* `prop`:           Propagator function to use. Standard is `propagator()`
* `reflect`:        If nothing (standar value), the axion-induced signal is computed.
                    If set, this field defines a beam, for which the reflected beam will be calculated.
                    A vector with the mode coefficients of the modes. To compute it from an arbitrary field distribution use `field2modes()`.
* `emit`:           If nothing (standar value), the axion-induced signal is computed.
                    If set, this array contains all emitted fields.
                    A vector with the mode coefficients of the modes. To compute it from an arbitrary field distribution use `field2modes()`.
* `Xset`, `Yset`:   Explicitly set the coordinate system for the fields
* `returnsum`:      (only `dancer`) If false, the out-propagating contributions after each iteration will be returned, without summing.
* `immediatesum`:   (only `dancer`) If false, the out-propagating contributions will be saved and summed up at the end.
* `returnboth`:     (only `cheerleader`) If true, the out-propagating waves on the right and on the left are returned. Otherwise only the out-propagating wave on the right.

Returns:
Vector with the mode coefficients of the boosted fields. To compute the field distribution from this use `field2modes()`.

### Initilization
```julia
# Coordinate System
X = [0.001]
Y = [0.001]

@everywhere sbdry = SeedSetupBoundaries(1)

#Permittivity
epsiolon = 9
sbdry.eps = [1e30,1.0,epsilon, 1.0] # <-- The mirror is included as a region with high permittivity
#Tilts
sbdry.relative_tilt_x = zeros(4);
sbdry.relative_tilt_y = zeros(4);
#Surface Roughness
sbdry.relative_surfaces = zeros(4,length(X),length(Y));

# Initialize Coordiates and Modes
init_coords(X,Y)
init_modes_1d()
```
