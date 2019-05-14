# 1D Model

The `src` directory contains a file `Analytical1D.py` containing a 1D calculation for various disk spacings in Python.

The easiest way to use it in julia is to copy it to your project directory and call it via
```julia
using PyCall
A1DM = pyimport("Analytical1D")
frequencies = 10:0.0001:30
spacings = [2,3,4,5,6,2,3,4,5,6,7,0]*1e-3 # First value is distance between mirror and first disk
ref, ax = A1DM.disk_system(Array{Float64}(frequencies*1e9),
        num_disk=11, disk_epsilon=25, mirror=true,
        spacings=spacings)
```


The last spacing is the distance to the antenna and only changes the phase of the result.
The disk thickness is 1mm, but can be changed (see code...).,



`ref` contains an array with the _amplitude reflecitivtiy_ of the system depending on frequency.

`ax` contains an array with the _amplitude boost factor_ of the system depending on frequency.

To get the reflected power / power boost, the absolute value of these quantities need to be squared.
To get the group delay, you have to take numerically the derivative of the phase of the reflecticity over the angular frequency.


### ToDo
It remains to clean up this file, translate it to julia such that maybe a automatic differentiation could be applied. This could be very, very useful to do optimizations of the disks spacings.
