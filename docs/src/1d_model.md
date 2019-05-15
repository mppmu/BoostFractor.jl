# 1D Model

## Usage
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

## Optimization
In order to compare our simulations, I think it is good to use always the same disk spacing configurations. A list of spacings can be found on Confluence, see [MADMAX / Working Groups / Theory, simulations and analysis / Disk Positions for various Boost Factors](https://confluence.desy.de/pages/viewpage.action?pageId=118280190).

Work on the optimizations should focus on improving the optimization procedure itself together with the proof of principle setup or learning about the 1D model.

You can use any Python or Julia optimization package to optimize boost factors.

For example you can run

```julia
# Optimization package for Julia
using Optim
# 1D Model
A1DM = pyimport("Analytical1D")


# Frequencies to look at (e.g. for a plot)
frequencies = 21.98:0.002:22.12
# Frequencies to maximize the boost within
freqmin = 22.0
freqmax = 22.045
cond = (frequencies .> freqmin).*(frequencies .< freqmax)

# A good initial guess of course may speed up the optimization drastically.
init = [1.00334, 6.94754, 7.1766, 7.22788, 7.19717,
        7.23776, 7.07746, 7.57173, 7.08019, 7.24657,
        7.21708, 7.18317, 7.13025, 7.2198, 7.45585,
        7.39873, 7.15403, 7.14252, 6.83105, 7.42282]*1e-3 

# The function to optimize
function minimizeme(arr)
    # You do not want to optimize the last spacing to the antenna,
    # because it does not change the boost factor
    # So arr only contains e.g. 20 spacings
    # but the 1D code needs 21, so we add one:
    
    distances = copy(arr)
    push!(distances,0)
    
    # We run the 1D code for the frequencies where we want to maximize the boost 
    ref, ax = A1DM.disk_system(Array{Float64}(frequencies*1e9)[cond], num_disk=20, disk_epsilon=24, spacings=distances)
    
    # Take the minimum value of the boost in that region and
    # return it after multiplying with (-1).
    # This value shall be minimized.
    return -minimum(abs2.(ax))
end

# Do the optimization
optim_result = Optim.optimize(minimizeme, init, NelderMead(), Optim.Options(x_tol=100e-9, f_tol=1e-2,g_tol=1e-5,show_trace=true,show_every=100,iterations=Int(5e4) ) )
# Get the resulting spacings
optimized_spacings = Optim.minimizer(optim_result)

# ... 
```

Of course, you can do this more sophiscticated, e.g. if you get negative spacings or multiple wavelength spacings, you may want to constrain the optimization range. You can also set a time limit for the optimization etc.
```julia
lower = zeros(20) # negative spacings would be a bit, well, impractical..
upper = 0.036*ones(20) # more than 2 wavelengths really isn't neccessary...
optim_result = Optim.optimize(minimizeme, lower, upper, init,  Fminbox(NelderMead()), Optim.Options(x_tol=100e-9, f_tol=1e-3,g_tol=1e-5,show_trace=true,show_every=100,time_limit=43200. ) )
```
For more details see [the `Optim` documentation](http://julianlsolvers.github.io/Optim.jl/stable/).

## ToDo
It remains to clean up the 1D model file, translate it to julia such that maybe a automatic differentiation could be applied. This could be very, very useful to do optimizations of the disks spacings.
