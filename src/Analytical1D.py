import numpy as np

c = 299792458.

def transform_surface(freq, gamma, axion_out=None, eps_i=1., eps_m=1., l=0, surfaceloss=0, USE_RAY=True):
    if USE_RAY:
        # Ray Tracing
        
        r_surface = (np.sqrt(eps_m)-np.sqrt(eps_i))/(np.sqrt(eps_m)+np.sqrt(eps_i))
        r_internal = gamma
        
        # Any transmission or reflection factor simply gets multiplied with this loss factor
        g = (1- surfaceloss)
        
        # The full relflectivity in front of the surface
        gamma = g*r_surface + g**2 * (1+r_surface)*(1-r_surface) * r_internal / (1 + r_internal*r_surface*g)
        
        # Propagate the reflecitivity is the same as for the impedance trfo...
    else:
        # Impedance Transformation
        
        # Calculate impedance from new reflection coefficient
        z = (1+gamma)/(1-gamma)
        
        Z_i = 1. / np.sqrt(eps_i) 
        Z_m = 1. / np.sqrt(eps_m)

        # Change basis to new medium
        z *= Z_i / Z_m

        # Calculate reflection coefficient
        # According to ray calculation that is not fully correct, only if gamma = 1
        gamma = ((z-1)/(z+1))*(1-surfaceloss)

    # Propagation (bidirectional)
    wavel = c / (freq * np.sqrt(eps_m))
    gamma *= np.exp(-1j*4*np.pi*l/wavel)

    # Calculate axion contribution
    if axion_out is not None:
        if not USE_RAY:
            error('Need to use ray approx.')
            
        denominator = eps_i * np.sqrt(eps_m) + eps_m * np.sqrt(eps_i)
        ax_i = np.sqrt(eps_i) * (1 - eps_m/eps_i) / denominator
        ax_m = np.sqrt(eps_m) * (1 - eps_i/eps_m) / denominator
        
        #Propagation of the internal new axion field inside and through the surface
        ax_i *= g*r_internal*(1-r_surface)/(1+g*r_internal*r_surface)
        #ax_i *= -g*(1-r_surface)*(-g*r_surface*r_internal)/((-g*r_surface*r_internal)-1)
        
        # Propagation of the present axion_out through the surface
        axion_out *= g*(1-r_surface)/(1 + r_internal*r_surface*g)
        
        # Total outgoing axion amplitude
        axion_out += ax_i + ax_m
        
        # Propagate the axion wave
        # Propagates in the direction of the reflected wave therefore the leading sign in the exponent is minus
        axion_out *= np.exp(-1j*2*np.pi*l/wavel)
        #print(axion_out)

    return gamma

def lossy_eps(freq, eps, tand, SvenssonDjordjevic=False):
    if not SvenssonDjordjevic:
        return (eps  -1j*eps*tand) #/np.sqrt(1+tand**2) <-- ADS is not deviding, apparently...
    else:
        FreqForEpsrTanD = 1e9 #Hz
        HighFreqForTanD = 1e12 #Hz
        LowFreqForTanD = 1e3 #Hz
        # Svensson/Djordjevic Model
        # http://edadocs.software.keysight.com/display/ads2009/About+Dielectric+Loss+Models
        L = np.log((HighFreqForTanD + 1j*FreqForEpsrTanD)/(LowFreqForTanD + 1j*FreqForEpsrTanD))
        a = - (eps*tand)/np.imag(L)
        Einf = eps - a*np.real(L)
        return (Einf  + a * np.log((HighFreqForTanD + 1j*freq)/(LowFreqForTanD + 1j*freq)) )

def disk_system(freq, tand=0, num_disk=5, non_uniform_surfaceloss=None, spacings=None, mirror=True, disk_thickness=0.001,disk_epsilon=9,**kwargs):
    if spacings is None:
        spacings = 8e-3*np.ones(num_disk+1)
    
    Z_0 = 0.001
    eps_1 = lossy_eps(freq, disk_epsilon, 0, SvenssonDjordjevic=False) #9.
    l_1 = 0.001
    # Traditional Loss Model
    eps_2 = 1.
    eps_2_tr = lossy_eps(freq, 1.0, tand, SvenssonDjordjevic=False)
    #eps_2_sd = lossy_eps(freq, 1.0, tand, SvenssonDjordjevic=True)
    l_2 = 0.
    
    # Air 0
    #gamma = transform_surface(freq, 0., eps_i=1., eps_m=1., l=0, **kwargs)
    
    #if mirror:
    #    axion_out=np.ones(freq.shape) + 0*1j
    #else:
    axion_out=np.zeros(freq.shape) + 0*1j
    
    # Metal Disk
    if mirror:
        gamma = transform_surface(freq, 0.0, axion_out=axion_out, eps_i=1., eps_m=1e20, l=100e-3, **kwargs)
    else:
        gamma = 0
        
    # print(axion_out)
        
    if non_uniform_surfaceloss is None:
        
        # Air 1
        # 1. = Z_0 / Z_0 (at the beginning the normalized impedance is always 1)
        gamma = transform_surface(freq, gamma, axion_out=axion_out, eps_i=(1e20 if mirror else 1.), eps_m=eps_2_tr, l=spacings[0], **kwargs)
        #print (gamma)
        for i in np.arange(0,num_disk):
            # Disk i+1
            gamma = transform_surface(freq, gamma, axion_out=axion_out, eps_i=eps_2_tr, eps_m=eps_1, l=disk_thickness, **kwargs)
            # Air i+2
            gamma = transform_surface(freq, gamma, axion_out=axion_out, eps_i=eps_1, eps_m=eps_2_tr, l=spacings[i+1], **kwargs)
    else:        
        # Air 1
        # 1. = Z_0 / Z_0 (at the beginning the normalized impedance is always 1)
        gamma = transform_surface(freq, gamma, axion_out=axion_out, eps_i=1e20, eps_m=eps_2_tr, l=spacings[0], surfaceloss=non_uniform_surfaceloss[0], **kwargs)
    
        for i in np.arange(0,num_disk):
            # Disk i+1
            gamma = transform_surface(freq, gamma, axion_out=axion_out, eps_i=eps_2_tr, eps_m=eps_1, l=disk_thickness,  surfaceloss=non_uniform_surfaceloss[2*i+1],**kwargs)
            # Air i+2
            gamma = transform_surface(freq, gamma, axion_out=axion_out, eps_i=eps_1, eps_m=eps_2_tr, l=spacings[i+1],  surfaceloss=non_uniform_surfaceloss[2*i+2],**kwargs)
        
        
    return gamma, axion_out

def disk_system_phase_depths(d_air,d_disk=np.pi/2, tand=0, num_disk=5, non_uniform_surfaceloss=None, **kwargs):

    # The result should be independent of frequency now...
    freq = 20e9
    wavel = c / (freq)
    l_air  = d_air*wavel/(2*np.pi)
    l_disk = d_disk*wavel/(2*np.pi*np.sqrt(9.4))
    
    
    Z_0 = 0.001
    eps_1 = lossy_eps(freq, 9.4, 0, SvenssonDjordjevic=False) #9.
    l_1 = 0.001
    # Traditional Loss Model
    eps_2 = 1.
    eps_2_tr = lossy_eps(freq, 1.0, tand, SvenssonDjordjevic=False)
    #eps_2_sd = lossy_eps(freq, 1.0, tand, SvenssonDjordjevic=True)
    l_2 = 0.
    
    
    # Air 0
    #gamma = transform_surface(freq, 0., eps_i=1., eps_m=1., l=0, **kwargs)
    
    axion_out=np.zeros(d_air.shape) + 0*1j
    # Metal Disk
    gamma = transform_surface(freq, 0.0, axion_out=axion_out, eps_i=1., eps_m=1e20, l=100e-3, **kwargs)
    
    if non_uniform_surfaceloss is None:
        
        # Air 1
        # 1. = Z_0 / Z_0 (at the beginning the normalized impedance is always 1)
        gamma = transform_surface(freq, gamma, axion_out=axion_out, eps_i=1e20, eps_m=eps_2_tr, l=l_air, **kwargs)
    
        for i in np.arange(0,num_disk):
            # Disk i+1
            gamma = transform_surface(freq, gamma, axion_out=axion_out, eps_i=eps_2_tr, eps_m=eps_1, l=l_disk, **kwargs)
            # Air i+2
            gamma = transform_surface(freq, gamma, axion_out=axion_out, eps_i=eps_1, eps_m=eps_2_tr, l=l_air, **kwargs)
    else:        
        # Air 1
        # 1. = Z_0 / Z_0 (at the beginning the normalized impedance is always 1)
        gamma = transform_surface(freq, gamma, axion_out=axion_out, eps_i=1e20, eps_m=eps_2_tr, l=8e-3, surfaceloss=non_uniform_surfaceloss[0], **kwargs)
    
        for i in np.arange(0,num_disk):
            # Disk i+1
            gamma = transform_surface(freq, gamma, axion_out=axion_out, eps_i=eps_2_tr, eps_m=eps_1, l=l_disk,  surfaceloss=non_uniform_surfaceloss[2*i+1],**kwargs)
            # Air i+2
            gamma = transform_surface(freq, gamma, axion_out=axion_out, eps_i=eps_1, eps_m=eps_2_tr, l=l_air,  surfaceloss=non_uniform_surfaceloss[2*i+2],**kwargs)
        
        
    return gamma, axion_out

