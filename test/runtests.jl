# This file is a part of BoostFractor.jl, licensed under the MIT License (MIT).

# Note: Travis doesn't seem to like getting no output for 10 mins.
#       Re-added status updates for CHEERLEADER and DANCER for this reason.

using Test
using Distributed
using FFTW

@testset "Package BoostFractor" begin

@everywhere begin
    PATH = ".."
    #PATH = "BoostFractor.jl"
    # General Interface
    include(PATH*"/src/boost_fractor.jl")

    # Differnt algorithms
    include(PATH*"/src/boost_fractor_recurprop.jl")
    include(PATH*"/src/boost_fractor_transformer.jl")

    # Some convenient tools
    include(PATH*"/src/beam_gaussian.jl")
    include(PATH*"/src/beam_coupling.jl")
end

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
    distance = [0.0, 0.15, 0.15/sqrt.(epsilon), 0.0]*1e-2
    
    sbdry = SeedSetupBoundaries(coords, diskno=1, distance=distance, epsilon=eps)
    
end

# Needs to be behind cheerleader init to define coords!
gauss_shape = gauss_profile(coords; z = 1e-9, omega0 = 0.097, f = 22e9)
@test gauss_shape[1] ≈ 8.343077704841396e-24 - 4.0474183340928325e-30im

# Let us sweep over different disk phase depths for a disk at distance lambda/2 infront of the mirror
# Of course the same way you can create a parallelized frequency sweep, etc.
@everywhere Deps = 0:0.1:3*pi
@everywhere zcat(args...) = cat(dims = 3, args...)
@time result = @sync @distributed (zcat) for i in 1:length(Deps)
    deps = Deps[i]
    dnu = pi # One could optimize this for maximal boost in 1D before; its left out here for simplicitly and just set to pi
    print("Running dEps = $deps (dNu = $dnu)... \n")

    # Set the distances (the wavelength in free space is 30mm, in the medium 10mm (free space devided by sqrt(epsilon))
    sbdry.distance = [0.0,(dnu)*3e-2/(2*pi), (deps)*1e-2/(2*pi),0.e-3];

    nmax = 100 # Maximum Number of Iterations
    #1D (if you calculate in 1D, use (vcat) above instead of (zcat) )
    #cheerleader(0,nmax,sbdry; prop=propagator1D,f=10e9,Xset=X,Yset=Y)
    #3D    
    cheerleader(0,nmax,sbdry,coords; prop=propagator,f=10e9,diskR=diskR)
end;

boost_factor_cheerleader = sum(abs2.(result), dims=(1,2))[1,1,:]*dx*dx/(pi*diskR^2);

println("Finished testing for CHEERLEADER")

#--------------------------------------------------------------------------------------------
### Test DANCER ###
@everywhere begin
    # Initialize Coordinate System
    dx = 0.01
    coords = SeedCoordinateSystem(X = -0.5:dx:0.5, Y = -0.5:dx:0.5)
    
    diskR = 0.15
    
    # SetupBoundaries
    # No Tilt or surface roughness here.
    epsilon = 9
    eps = Array{Complex{Float64}}([NaN,1,epsilon,1])
    distance = [0.0, 0.15, 0.15/sqrt.(epsilon), 0.0]*1e-2
    
    sbdry = SeedSetupBoundaries(coords, diskno=1, distance=distance, epsilon=eps)
    
end

# Let us sweep over different disk phase depths for a disk at distance lambda/2 infront of the mirror
# Of course the same way you can create a parallelized frequency sweep, etc.
@everywhere Deps = 0:0.1:3*pi
@everywhere zcat(args...) = cat(dims = 3, args...)
@time result_dancer = @sync @distributed (zcat) for i in 1:length(Deps)
    deps = Deps[i]
    dnu = pi # One could optimize this for maximal boost in 1D before; its left out here for simplicitly and just set to pi
    print("Running dEps = $deps (dNu = $dnu)... \n")

    # Set the distances (the wavelength in free space is 30mm, in the medium 10mm (free space devided by sqrt(epsilon))
    sbdry.distance = [0.0, (dnu)*3e-2/(2*pi), Deps[i]*1e-2/(2*pi),0];

    nmax = 200 # Maximum Number of Iterations
    #1D (if you calculate in 1D, use (vcat) above instead of (zcat) )
    #dancer(0,nmax,sbdry; prop=propagator1D,f=10e9,Xset=[0.0000001],Yset=[0.0000001])
    #3D    
    dancer(0,nmax,sbdry,coords; prop=propagator,f=10e9,diskR=diskR)
end;

boost_factor_dancer = sum(abs2.(result_dancer), dims=(1,2))[1,1,:]*dx*dx/(pi*diskR^2);

println("Finished testing for DANCER")

#--------------------------------------------------------------------------------------------
### Test TRANSFORMER ###
@everywhere begin
    # Initialize Coordinate System
    dx = 0.01
    coords = SeedCoordinateSystem(X = -0.5:dx:0.5, Y = -0.5:dx:0.5)
    
    diskR = 0.15
    
    # SetupBoundaries
    # No Tilt or surface roughness here.
    epsilon = 9
    eps = Array{Complex{Float64}}([NaN,1,epsilon,1])
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

# Let us sweep over different disk phase depths for a disk at distance lambda/2 infront of the mirror
# Of course the same way you can create a parallelized frequency sweep, etc.
@everywhere Deps = 0:0.1:3*pi
@time result_transformer = @sync @distributed (hcat) for i in 1:length(Deps)
    deps = Deps[i]
    dnu = pi # One could optimize this for maximal boost in 1D before; its left out here for simplicitly and just set to pi
    print("Running dEps = $deps (dNu = $dnu)... \n")

    # Set the distances (the wavelength in free space is 30mm, in the medium 10mm (free space devided by sqrt(epsilon))
    sbdry.distance = [0.0, (dnu)*3e-2/(2*pi), Deps[i]*1e-2/(2*pi),0];

    #1D (if you calculate in 1D, use (vcat) above instead of (hcat) )
    #transformer(sbdry; prop=propagator1D,f=10e9,Xset=[1e-7],Yset=[1e-7])
    #3D    
    transformer(sbdry,coords,modes; prop=propagator,f=10e9,diskR=diskR)
end;

boost_factor_transformer = sum(abs2.(result_transformer), dims=1)[1,:];

println("Finished testing for TRANSFORMER")
#--------------------------------------------------------------------------------------------

@test boost_factor_cheerleader ≈ [0.9541048836881387, 0.9428833993720742, 0.968920462905896, 1.0341129917881549, 1.144486522206282, 1.311281404330431, 1.5530311134049792, 1.8991893743757546, 2.3962627706758717, 3.1180836339325797, 4.182537182036321, 5.77631829548155, 8.180541344734435, 11.750104248034686, 16.683053458978392, 22.328276414700372, 26.584187537779602, 27.43539899358354, 25.287916750692773, 21.957711497592978, 18.70293960502623, 15.968880656976728, 13.809281087444926, 12.13342947587664, 10.836197159631562, 9.835302811952268, 9.062953656331851, 8.469762044594315, 8.021400572162888, 7.692579478167976, 7.465918522414759, 7.329994000109685, 7.2782449520717565, 7.308492142283055, 7.422848236181687, 7.6280184632230785, 7.936389069538961, 8.36766853387498, 8.951220499791912, 9.729781130188297, 10.764938064402624, 12.144089632575064, 13.986113592200057, 16.433651816395216, 19.591343308612494, 23.3074109612804, 26.695661251021864, 27.883116341877454, 25.476205860505605, 20.569262201812844, 15.415614620919474, 11.236389682727756, 8.18288180649727, 6.037430281779245, 4.5362561494977, 3.473161524397452, 2.713994325825552, 2.16586121911677, 1.7663720636616334, 1.4758050703937462, 1.2667394643264536, 1.1207427057453578, 1.025753566948445, 0.9741434315936932, 0.961917112768637, 0.9883369519635989, 1.0560617558520191, 1.1717821857313868, 1.3476155646302996, 1.6035770542564973, 1.9715312870420645, 2.501814156824845, 3.2745053970180558, 4.417660301811236, 6.13329514001611, 8.719975781044885, 12.530192745139441, 17.668994151816275, 23.226061744057777, 26.92556913789544, 27.119751483329694, 24.648270973331165, 21.3056595232044, 18.163653151444, 15.554966421842678, 13.502429127893384, 11.911407897285025, 10.675141498985058, 9.71876066215671, 8.981342488510444, 8.413847982210477, 7.984858529840069, 7.671899347016415, 7.458149018169155, 7.33352793792804]

@test boost_factor_dancer ≈ [0.9465597604600772, 0.9353352136853974, 0.9613500833206581, 1.0265038271281812, 1.1368236550259119, 1.3035506496401894, 1.5452177859585288, 1.8912763453371952, 2.3882277762141007, 3.1098957209819647, 4.174151879779258, 5.767671034944011, 8.171539715302092, 11.740624274034273, 16.6729565465068, 22.31747736268758, 26.572764341320273, 27.423602923437432, 25.276015311731925, 21.94587203360396, 18.69124362921158, 15.957338828020001, 13.797876691456961, 12.122159028174904, 10.825042594928846, 9.824242192983457, 9.05197733385444, 8.458858565994458, 8.010564389351762, 7.68180970131266, 7.4552156931162274, 7.3193608505429895, 7.267685228467108, 7.298011723964006, 7.412455350301969, 7.6177193601016056, 7.926182866326155, 8.357549109013213, 8.941179575288242, 9.719807990020836, 10.755018242265454, 12.134204765439094, 13.97624329007452, 16.423779229700262, 19.58146822920849, 23.29757595456879, 26.685987721853778, 27.873800005446203, 25.4673997532372, 20.56095299917797, 15.407657470121954, 11.228632064574768, 8.175203196725251, 6.029722950269289, 4.52845052200948, 3.465206267228134, 2.7058182554234054, 2.157406992392685, 1.7575956843424423, 1.4666424302967727, 1.2571295318833156, 1.1106325891076447, 1.015090899155427, 0.9628837252091647, 0.9500194106597534, 0.9757574022771622, 1.0427511871062949, 1.1576808370114224, 1.3326197879771289, 1.5875240918437874, 1.954222605311172, 2.48302081927242, 3.2539303359739677, 4.394881308723679, 6.107709566775495, 8.690766655560651, 12.496389375253083, 17.629819645108952, 23.18174731282909, 26.878049544460755, 27.071896250013232, 24.602197593707817, 21.262213904229338, 18.123030553046704, 15.516867175830251, 13.466334807311917, 11.877078481350779, 10.642339693624438, 9.687101909202692, 8.950609350269554, 8.383897489809803, 7.955480021177518, 7.6429588384112215, 7.4295446543075006, 7.30513529006866]

@test boost_factor_transformer ≈ [0.9219455134201293, 0.9088200235288822, 0.932472788773596, 0.9944207593641872, 1.1003581017820354, 1.2611863926386724, 1.4949283227109704, 1.830218075631676, 2.312297724310878, 3.013071038825134, 4.04748065766292, 5.597811746319159, 7.939063571323361, 11.419642965451624, 16.237339442286526, 21.762153491791135, 25.94078520367847, 26.792095591972686, 24.702468292477356, 21.447386690530784, 18.261152610706343, 15.58295201463342, 13.466455659627735, 11.823057780821372, 10.551035900830932, 9.569460134669937, 8.811943274009819, 8.230495827772176, 7.791034891609431, 7.468848293066695, 7.246839258650652, 7.113709251557361, 7.063016802684704, 7.092586196139127, 7.2042465260699515, 7.404537357496414, 7.705836236739639, 8.127673560768761, 8.698972696379933, 9.461765129881503, 10.47662887719731, 11.829593496954038, 13.637801537929464, 16.041996446339745, 19.145637191639622, 22.79913203746825, 26.127752327000355, 27.285415126839293, 24.90699832909787, 20.08049151233705, 15.022228781938072, 10.924340523719358, 7.934461008467611, 5.835365830636519, 4.365617413305234, 3.326589761075334, 2.585469801238749, 2.0507559250463068, 1.6628485969509776, 1.3823877426172526, 1.1827220364110236, 1.0458874408021794, 0.9597096314284665, 0.9162076426336727, 0.9108154369913783, 0.942226086594848, 1.012766421097211, 1.1291151754297613, 1.3034670791757523, 1.5557165291069917, 1.9173814408054293, 2.4382970793440837, 3.1977944804460456, 4.322700464678026, 6.013034962036725, 8.564165028005242, 12.323314613788607, 17.38715683051501, 22.83894840807172, 26.422258677978554, 26.539784031486974, 24.058135524654077, 20.74801236360884, 17.65094597322297, 15.093879494161175, 13.086093286518972, 11.528435507644087, 10.323654342883119, 9.394058173303977, 8.67688545686285, 8.127709799564812, 7.714600528320907, 7.414617399406159, 7.2120011630079555, 7.096533830554453]

end # testset
