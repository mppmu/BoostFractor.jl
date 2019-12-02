# BoostFractor.jl

<img src="img/boost_fractor_logo.png" alt="BoostFractor" width=300> <!--BoostFractor.jl-->

MadMax electrodynamics simulation in Julia.


## 1D Calculations
[1D Code](1d_model.md) based on impedance transformation / 1D ray tracing.

All 3D Codes below generalize also to 1D, but may be less performant.

## 3D Calculations
#### Recursive (Fourier) Propagation
Iteratively propagate electromagnetic fields between dielectric interfaces.
Implementations:
 * [Dancer](3d_algorithms.md#dancer_and_cheerleader)
 <!-- --->
 * [Cheerleader](3d_algorithms.md#dancer_and_cheerleader): somewhat more elegant implementation, needs only half as much iterations

More details on how the algorithms work can be found in [literature](resources.md).
#### Mode Matching
Solve boundary conditions on disk interfaces by matching (eigen)modes between the different regions. Since only a few modes need to be considered for a large boost factor, transfer matrices can be used.
 * [Transformer](3d_algorithms.md#transformer): Generalized Transfer-Matrix implementation in 3D.
 
## Example Code
The `examples` folder contains Jupyter notebooks with examples on how to use `BoostFractor`.
A summary of the contained examples can be found [here](examples.md)

## Disk Position Optimization
[1D Code](1d_model.md) also contains some examples how to use julia to do optimizations (of the disk spacings for example).
