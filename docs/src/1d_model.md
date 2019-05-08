# 1D Model

The `src` directory contains a file containing a 1D calculation for various disk spacings.

You can use it via
'''julia
using PyCall
A1DM = pyimport("Analytical1DMod")
frequencies = 10:0.0001:30
spacings = [2,3,4,5,6,2,3,4,5,6,7,0]*1e-3 # First value is distance between mirror and first disk
ref, ax = A1DMM.disk_system(Array{Float64}(frequencies*1e9),
        num_disk=11, disk_epsilon=25, mirror=true,
        spacings=ones(12)*1e-2)

'''
The last spacing is the distance to the antenna and only changes the phase of the result.
The disk thickness is 1mm, but can be changed (see code...).,

`ref` contains an array with the _amplitude reflecitivtiy_ of the system depending on frequency.
`ax` contains an array with the _amplitude boost factor_ of the system depending on frequency.
To get the reflected power / power boost, the absolute value of these quantities need to be squared.
To get the group delay, you have to take numerically the derivative of the phase of the reflecticity over the angular frequency.

