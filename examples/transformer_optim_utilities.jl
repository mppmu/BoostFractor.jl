
using BoostFractor
using Interpolations
using DSP


"""
Calculates propagation matrices for all regions and frequencies at given relative spacings (optionally relative tilts aswell)
The returned matrix is (n_regions x n_freq x n_spacing x n_tilt_x x n_tilt_y x n_mode x n_mode)      
"""
function calc_propagation_matrices_grid(sbdry::SetupBoundaries,coords::CoordinateSystem, modes::Modes,spacing_grid,frequencies;tilt_x_grid=0,tilt_y_grid=0, prop=propagator, diskR=0.15)
    n_region = length(sbdry.distance)
    n_disk = (n_region-2)÷2
    n_freq = length(frequencies)
    n_spacing = length(spacing_grid)
    n_tilt_x = length(tilt_x_grid)
    n_tilt_y = length(tilt_y_grid)
    n_mode = modes.M*(2*modes.L+1)
    #Split booster into air gaps and solids(disks+mirror) since the latter stay constant and need just one calculation
    #for each frequency
    ind_gap = 2:2:n_region
    ind_solid = 1:2:n_region-1

    sbdry_gap = SeedSetupBoundaries(coords, diskno=n_disk, distance=sbdry.distance[ind_gap],
                                        epsilon=sbdry.eps[ind_gap],relative_tilt_x=sbdry.relative_tilt_x[ind_gap],
                                        relative_tilt_y=sbdry.relative_tilt_y[ind_gap],
                                        relative_surfaces=sbdry.relative_surfaces[ind_gap,:,:])

    sbdry_solid = SeedSetupBoundaries(coords, diskno=n_disk, distance=sbdry.distance[ind_solid],
                                        epsilon=sbdry.eps[ind_solid],relative_tilt_x=sbdry.relative_tilt_x[ind_solid],
                                        relative_tilt_y=sbdry.relative_tilt_y[ind_solid],
                                        relative_surfaces=sbdry.relative_surfaces[ind_solid,:,:])

    distance_0 = copy(sbdry_gap.distance)
    tilt_x_0 = copy(sbdry_gap.relative_tilt_x)
    tilt_y_0 = copy(sbdry_gap.relative_tilt_y)
    prop_matrix_grid = Array{Complex{Float64},7}(undef,n_region,n_freq,n_spacing,n_tilt_x,n_tilt_y,n_mode,n_mode)


    Threads.@threads for f in 1:n_freq
        prop_matrix_solid = calc_propagation_matrices(sbdry_solid,coords,modes;f=frequencies[f],prop=prop,diskR=diskR)
        for s in 1:n_spacing, tx in 1:n_tilt_x, ty in 1:n_tilt_y
            sbdry_i = copy_setup_boundaries(sbdry_gap,coords) #For thread safety
            sbdry_i.distance = distance_0 .+ spacing_grid[s]
            sbdry_i.relative_tilt_x = tilt_x_0 .+ tilt_x_grid[tx]
            sbdry_i.relative_tilt_y = tilt_y_0 .+ tilt_y_grid[ty]

            prop_matrix_gap = calc_propagation_matrices(sbdry_i,coords,modes;f=frequencies[f],prop=prop,diskR=diskR)
            prop_matrix_grid[ind_gap,f,s,tx,ty,:,:] = [(prop_matrix_gap[r][k,j]) for r in 1:n_disk+1, k in 1:n_mode, j in 1:n_mode]
            prop_matrix_grid[ind_solid,f,s,tx,ty,:,:] = [(prop_matrix_solid[r][k,j]) for r in 1:n_disk+1, k in 1:n_mode, j in 1:n_mode]
        end;
    end

    return prop_matrix_grid
end;


"""
Constructs the interpolation object from Interpolations without tilts
"""
function construct_prop_matrix_interpolation(prop_matrix_grid::Array{Complex{Float64},7}, spacing_grid)
    n_region = size(prop_matrix_grid,1)
    n_freq = size(prop_matrix_grid,2)
    n_spacing = size(prop_matrix_grid,3)
    n_tilt_x = size(prop_matrix_grid,4)
    n_tilt_y = size(prop_matrix_grid,5)
    n_mode = size(prop_matrix_grid,6)    
    #Construct the interpolation object
    itp_fun = Interpolations.BSpline(Cubic(Natural(OnCell())))
    itp = Interpolations.interpolate(prop_matrix_grid, (NoInterp(),NoInterp(),
                                        itp_fun,NoInterp(),
                                        NoInterp(),NoInterp(),NoInterp()))
    itp = Interpolations.scale(itp,1:n_region,1:n_freq,spacing_grid,1:n_tilt_x,1:n_tilt_y,1:n_mode,1:n_mode)  
    return itp
end;

"""
Calculate interpolated propagation matrices set without tilts
"""
function interpolate_prop_matrix(itp,dist_shift::Array{T,1}) where T<:Real
    n_region = size(itp,1)
    n_freq = size(itp,2)
    n_mode = size(itp,6)
    #Disk thickness stays constant and last air gap stay constant 
    dist_shift_all = [(r+1)%2==0 ? 0.0 : r==n_region ? 0 : dist_shift[r÷2] for r in 1:n_region]
    prop_matrix_set_interp = Array{Array{Complex{T},2}}(undef,n_region,n_freq)
    for f=1:n_freq
        prop_matrix_set_interp[:,f] = [itp(r,f,dist_shift_all[r],1,1,1:n_mode,1:n_mode) for r in 1:n_region]
    end
    return prop_matrix_set_interp
end;



"""
Calculates boostfactor for given frequencies and booster. Note that the distances in sbdry are meaningless since
propagation_matrices_set already contains the effect of spacings.
"""        
function calc_boostfactor_modes(sbdry,coords,modes, frequencies, prop_matrices_set::Array{Array{Complex{T},2},2}; prop=propagator, diskR=0.15) where T<:Real
    n_freq = length(frequencies)
    n_modes = size(prop_matrices_set[1,1])[1]
    EoutModes0 = Array{Complex{T},3}(undef,1,n_modes,n_freq)
    # Sweep over frequency
    Threads.@threads for f in 1:n_freq
        boost = transformer(sbdry,coords,modes; prop=prop,diskR=diskR,f=frequencies[f],propagation_matrices=prop_matrices_set[:,f],reflect=nothing)
        EoutModes0[1,:,f] =  boost
    end
    return EoutModes0
end;


"""
Calculates boostfactor and reflectivity
"""
function calc_modes(sbdry,coords,modes, frequencies, prop_matrices_set::Array{Array{Complex{T},2},2},reflect;prop=propagator, diskR=0.15) where T<:Real
    n_freq = length(frequencies)
    n_modes = size(prop_matrices_set[1,1])[1]
    EoutModes0 = Array{Complex{T},3}(undef,2,n_modes,n_freq)
    # Sweep over frequency
    Threads.@threads for f in 1:n_freq
        boost, refl = transformer(sbdry,coords,modes; prop=prop,diskR=diskR,f=frequencies[f],propagation_matrices=prop_matrices_set[:,f],reflect=reflect)
        EoutModes0[1,:,f] =  boost
        EoutModes0[2,:,f] =  refl
    end 
    return EoutModes0
end;



function calc_boostfactor_cost(dist_shift::Array{T,1},itp,frequencies,sbdry::SetupBoundaries,coords::CoordinateSystem,modes::Modes,m_reflect;prop=propagator, diskR=0.15) where T<:Real
    dist_bound_hard = Interpolations.bounds(itp)[3]
    #Return hard penalty when exceeding interpolation bounds
    if any(.!(dist_bound_hard[1] .< dist_shift .< dist_bound_hard[2])) 
        return 1000.0
    end
    #Add soft penalty when approaching interpolation bounds
    penalty = soft_box_penalty(dist_shift,dist_bound_hard)

    prop_matrices_set_interp = interpolate_prop_matrix(itp,dist_shift);
    Eout = calc_boostfactor_modes(sbdry,coords,modes,frequencies,prop_matrices_set_interp,prop=prop, diskR=diskR)
    cpld_pwr = abs2.(sum(conj.(Eout[1,:,:]).*m_reflect, dims=1)[1,:])
    cost =  -p_norm(cpld_pwr,-20)*penalty
    return cost
end;

"""
This adds a soft barrier penalty when dist_shift is approaching the maximum shifts allowed.
Hard box contrains usually confuse optimizers
""" 
function soft_box_penalty(shift::Array{T,1},bound_hard) where T<:Real
    soft_bound_depth = (bound_hard[2]-bound_hard[1])*0.05
    l_soft_bound = bound_hard[1] + soft_bound_depth
    u_soft_bound = bound_hard[2] - soft_bound_depth

    excess_pos = maximum([shift;sum(shift)]) - u_soft_bound
    excess_neg = -(minimum([shift;sum(shift)]) - l_soft_bound)
    excess = maximum([excess_pos,excess_neg,0])
    penalty = 1 - (excess/soft_bound_depth)^6
    return penalty
end





function surface_roughness(X,Y,ngaps; mag=1e-4,trunc=1e-2, xi=nothing, diskR=0.15)
    # TODO: Some kind of correlation length etc to have better control
    #over the power spectrum
    init = mag*[randn() for i in 1:ngaps, x in X, y in Y];
    # Trunctrate standard deviation such that huge outliers can't happen
    # TODO: This makes delta peaks at +- trunc with magnitude of what
    #we have left.numerical stable softmax
    # Better: redistribute this or directly use sampling from
    #trunctrated distribution
    if xi !== nothing
        # Now convolute with a gaussian
        g = [exp(-(x^2 + y^2)/(2*xi^2)) for x in X, y in Y]

        for i in 1:ngaps
        init[i,:,:] = conv2(init[i,:,:],
            g)[Int(floor(length(X)/2)):Int(floor(length(X)*1.5))-1,Int(floor(length(X)/2)):Int(floor(length(X)*1.5))-1]
        end

        # Now renormalize to have the rms mag
        init ./= sqrt.(sum(abs2.(init),
        dims=(2,3))./(length(X)*length(Y)))
        init .*= mag
    end

    init[init .> trunc] .= trunc
    init[init .< -trunc] .= -trunc

    init .*= [(x^2 + y^2) <= diskR^2 for i in 1:ngaps, x in X, y in Y]
    return init
end;

"""
Helper function to quickly copy SetupBoundaries
"""
function copy_setup_boundaries(sbdry::SetupBoundaries,coords::CoordinateSystem)
    n_disk = (length(sbdry.distance)-2)/2
    sbdry_new = SeedSetupBoundaries(coords, diskno=n_disk, distance=copy(sbdry.distance), epsilon=copy(sbdry.eps),relative_tilt_x=copy(sbdry.relative_tilt_x), relative_tilt_y=copy(sbdry.relative_tilt_y), relative_surfaces=copy(sbdry.relative_surfaces))
    return sbdry_new
end;

"""
p norm (https://en.wikipedia.org/wiki/Lp_space#The_p-norm_in_finite_dimensions)
For large positive/negative p this is a differentiable approximatation of max(X)/min(X) 
Beware of numerical instability for large p. p~20 seems fine.
"""
function p_norm(X::Array{T,1},p) where T<:Real        
    magnitude = maximum(X)
    X_norm = X ./ magnitude
    return magnitude*(1/length(X) * sum(X_norm.^p))^(1/p)
end

 ############################################################################

