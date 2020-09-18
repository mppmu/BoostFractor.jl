#include <iostream>
#include <cstdio>
#include <cmath>
#include <math.h> 
/* try using FP_FAST_FMA -fused multipy add- for faster and more accurate multipy add operations */
#include <cufftdx.hpp>
#include <complex.h>
#include <chrono>

using namespace std;
using namespace cufftdx;

const unsigned int Grid_Size = 64;
const unsigned int Elements_Per_Thread = 4;

typedef std::complex<float> c_type;
typedef float f_type;
typedef float2 v_type;
typedef unsigned int uint;

const int disknumber = 10;
const unsigned int iterations = 1000;


//some constant memory on the GPU
__constant__ float r[disknumber * 2 + 2];
/* 
   C++ / CUDA implementation of the dancer algorithm.  

*/

/*
Some usefull strucs for the experimental setup (boundarys), coordiante system (coords), 
fields on the plates (fields) and the pre-generated phase shifts (phase_shifts)
*/
struct coords 
{ /* Coordinate grids in real and momentum space. square and symmetric around 0 */
    unsigned int size;

    float X_step;
    float X_max;

    float kX_step;
    float kX_max;

    coords (float max, float step) /* Constructor equivalent to init_coords */
    {
        X_step = step;
        X_max = max;

        size = 2 * X_max / X_step;

        kX_step = M_PI / max;
        kX_max = kX_step * (size - 1) / 2;
    }

    coords (float max, int len) /* alternative constructor using the number of elements along one dimension */
    {
        X_step = 2 * max / (len - 1);
        X_max = max;

        size = len;

        kX_step = M_PI / max;
        kX_max = kX_step * (size - 1) / 2;
    }

    /* functions for acessing values like in arrays */
    f_type x(int i)
    {
        return i * X_step - X_max;
    }

    f_type k(int i)
    {
        return i * kX_step - kX_max;
    }

    /* no dynamic memory so no destrucor needed */
};

struct boundarys 
{ /* boundary object can be created at compile time */ 

    int disknumber;
    float * distance; /* dynamic arrays because size not known at compile time */
    float radius;

    c_type * eps;
    float * r;

    float * tilt_x;
    float * tilt_y;
    c_type * surfaces; /* 3D array, linear in memory */

    boundarys (int disk_n, float rad, coords &cord) /* basic constructor for n disks. Maybe calculate r from eps */
    {
        disknumber = disk_n;
        radius = rad;

        /* dynamically allocate space for 1D arrays*/
        distance = new float[2 * disk_n + 1]; 
        eps = new c_type[2 * disk_n + 1];
        r = new float[2 * disk_n + 2];

        /* 
        initialize titlts to zero
        The standart zero initialisation new double [...]() seems to be causing issues */
        tilt_x = new float[2 * disk_n + 1]();
        tilt_y = new float[2 * disk_n + 1]();


        surfaces = new c_type [(2 * disk_n + 1) * cord.size * cord.size]();

        /* set values for distances and dielectric constant */
        for(int i = 0; i < disk_n; i++) 
        {
            distance[2*i] = 8e-3;
            distance[2*i+1] = 1e-3;

            eps[2*i] = 1.;
            eps[2*i+1] = 9;

            r[2*i+1] = -0.5;
            r[2*i+2] = 0.5;
        }
        distance[2*disk_n] = 0.;
        eps[2*disk_n] = 1.;
        r[0] = 1.;
        r[2*disk_n+1] = 0.;
    }

    /* write destructor for deleting dynamic memory */
    ~boundarys ()
    {
        delete[] distance;
        delete[] r;
        delete[] eps;
        delete[] tilt_x;
        delete[] tilt_y;
        delete[] surfaces;
    }
};

struct phase_shifts
{
  cudaTextureObject_t k_tex = 0;
  cudaTextureObject_t surface_tex = 0;

  cudaArray * k_arr;
  cudaArray * surface_arr;
  v_type * surface_out;

  phase_shifts (boundarys &bdry, coords &cord, float lambda)
  {
    //First Calcualte Phase shifts on CPU
    int n_regions = 2 * bdry.disknumber + 1;
    const unsigned int size2 = cord.size * cord.size;

    float2 * k;
    float2 * surface;
    
    k = new float2 [n_regions * size2];
    surface = new float2 [n_regions * size2];
    surface_out = new v_type [size2];

    int x_shift;
    int y_shift;
    int lin_index;

    c_type k0;
    c_type t_y;
    c_type t_x;
    c_type s;
    c_type a;

    float disk_cutoff;

    for(int z = 0; z < n_regions; z++)
      {
	k0 = 2 * float(M_PI) * sqrt(bdry.eps[z]) / lambda;

	for(int y = 0; y < cord.size; y++)
	  {
	    t_y = exp(c_type(0,-1) * k0 * bdry.tilt_y[z] * cord.x(y));

	    for(int x = 0; x < cord.size; x++)
	      {
		lin_index = x + cord.size * (y + cord.size * z);
		y_shift = (y + cord.size / 2) % cord.size;
		x_shift = (x + cord.size / 2) % cord.size;


		a = conj(sqrt( pow(k0,2) - pow(cord.k(x_shift), 2) - pow(cord.k(y_shift), 2)));

		k[lin_index].x = real(exp(c_type(0,-1) * bdry.distance[z] * a));
		k[lin_index].y = imag(exp(c_type(0,-1) * bdry.distance[z] * a));


		disk_cutoff = (pow(bdry.radius,2) > pow(cord.x(x),2) + pow(cord.x(y),2));

		t_x = exp(c_type(0,-1) * k0 * bdry.tilt_x[z] * cord.x(x));
		s = exp(c_type(0,-1) * k0 * bdry.surfaces[lin_index]);

		//Important: use 1.f in division for single pres. floating point results
		surface[lin_index].x = (1.f / size2) * disk_cutoff * real(t_x * t_y * s);
		surface[lin_index].y = (1.f / size2) * disk_cutoff * imag(t_x * t_y * s);

		if(z == n_regions - 1)
		  {
		    surface_out[x + y * cord.size].x = (1.f / size2) * real(t_x * t_y * s);
		    surface_out[x + y * cord.size].y = (1.f / size2) * imag(t_x * t_y * s);
		  }
	      }
	  }
      }

    //Allocate CUDA Array
    cudaChannelFormatDesc channel =
      cudaCreateChannelDesc<float2>();//(32, 32, 0, 0, cudaChannelFormatKindFloat);

    cudaExtent ext =
      make_cudaExtent(cord.size, cord.size, n_regions);

    cudaMalloc3DArray(&k_arr, &channel, ext, cudaArrayLayered);
    cudaMalloc3DArray(&surface_arr, &channel, ext, cudaArrayLayered);

    //Copy Host Data to Device

    //cudaMemcpyToArray(k_arr, 0, 0, k, size2 * n_regions * sizeof(float2),cudaMemcpyHostToDevice);

    //cudaMemcpyToArray(surface_arr, 0, 0, surface, size2 * n_regions * sizeof(float2), cudaMemcpyHostToDevice);

    cudaMemcpy3DParms copy_k = {0};
    cudaMemcpy3DParms copy_s = {0};

    copy_k.kind = copy_s.kind = cudaMemcpyHostToDevice;
    copy_k.extent = copy_s.extent = ext;

    copy_k.dstArray = k_arr;
    copy_s.dstArray = surface_arr;

    copy_k.srcPtr = make_cudaPitchedPtr((void*)k, ext.width * sizeof(float2), ext.width, ext.height);
    
    copy_s.srcPtr = make_cudaPitchedPtr((void*)surface, ext.width * sizeof(float2), ext.width, ext.height);

    cudaMemcpy3D(&copy_k);
    cudaMemcpy3D(&copy_s);

    // Specify Texture
    struct cudaResourceDesc resDesc_k;
    memset(&resDesc_k, 0, sizeof(resDesc_k));
    resDesc_k.resType = cudaResourceTypeArray;
    resDesc_k.res.array.array = k_arr;

    struct cudaResourceDesc resDesc_surf;
    memset(&resDesc_surf, 0, sizeof(resDesc_surf));
    resDesc_surf.resType = cudaResourceTypeArray;
    resDesc_surf.res.array.array = surface_arr;

    struct cudaTextureDesc texDesc;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.addressMode[0] = cudaAddressModeBorder;
    texDesc.addressMode[1] = cudaAddressModeBorder;
    texDesc.filterMode = cudaFilterModePoint;
    texDesc.readMode = cudaReadModeElementType;
    texDesc.normalizedCoords = false;

    cudaCreateTextureObject(& k_tex, & resDesc_k, & texDesc, NULL);
    cudaCreateTextureObject(& surface_tex, & resDesc_surf, & texDesc, NULL);
  }
  /* destructor */
  ~phase_shifts() 
  {
    cudaDestroyTextureObject(k_tex);
    cudaFreeArray(k_arr);
    
    cudaDestroyTextureObject(surface_tex);
    cudaFreeArray(surface_arr);
    
    delete[] surface_out;
  }
};

/* note that for the phase shifts, these are mutiplied each by 1/size^2 as fftw resusts in unnormalized ffts. */

/* simple function used for printing field slices to console for use in testing / debugging */
void print_field (float2 * field, coords & cord)
{
    int n_x = cord.size;

    for(int i = 0; i < n_x; i++)
    {
        for(int j = 0; j < n_x; j++)
        {
            printf( "%5.1f %5.1fi   ", field[n_x * i + j].x, field[n_x * i + j].y);
        }
        cout << endl;
    }
    cout << endl;
}

void print_prop (float2 * field, coords & cord, int zmax)
{
    int n_x = cord.size;
    for(int z = 0; z < zmax; z++)
    {
        for(int i = 0; i < n_x; i++)
        {
            for(int j = 0; j < n_x; j++)
            {
                printf( "%7.3f %7.3fi   ", field[n_x * (n_x * z + i) + j].x, field[n_x * (n_x * z + i) + j].y );
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << endl << endl;
}

/* 
Initialise the emission from the boundarys.
Combine with initial cutoff at the disk radius 
*/

void set_field (float2 * field, coords & cord, boundarys & bdry, c_type val)
{

    for(int i = 0; i < (cord.size * cord.size); i++)
    {
        /* Using boolean arithmetic may be faster than using if statments (see branchless programming) */
        bool disk_cutoff = pow(bdry.radius, 2) > (pow(cord.x(i % cord.size), 2) + pow(cord.x(i / cord.size), 2));
        field[i].x = disk_cutoff * real(val);
        field[i].y = disk_cutoff * imag(val);
    }
}
    
c_type ax (c_type eps_i, c_type eps_m)
{
    return sqrt(eps_m) * (1.f - eps_i / eps_m) / (eps_i * sqrt(eps_m) + eps_m * sqrt(eps_i));
}

void setup_fields (float2 * fields, boundarys & bdry, coords & cord) 
{
    const unsigned int n_regions = 2 * bdry.disknumber + 1;
    const unsigned int disk_offset = cord.size * cord.size;
    const unsigned int dir_offset = n_regions * disk_offset;

    set_field(&fields[0], cord, bdry, - 1.f);
    set_field(&fields[dir_offset], cord, bdry, ax(bdry.eps[1], bdry.eps[0]));

    /* standart cases */
    for(int n = 1; n < n_regions - 1; n++)
    {
        set_field(&fields[disk_offset * n], cord, bdry, ax(bdry.eps[n-1], bdry.eps[n]));
        set_field(&fields[dir_offset + disk_offset * n], cord, bdry, ax(bdry.eps[n+1], bdry.eps[n]));
    }

    set_field(&fields[disk_offset * (n_regions-1)], cord, bdry, ax(bdry.eps[n_regions-2], bdry.eps[n_regions-1]));
    set_field(&fields[dir_offset + disk_offset * (n_regions-1)], cord, bdry, ax(1., bdry.eps[n_regions-1]));
}

/*
The heart of it
*/
template<typename complex_type>
__device__ void phase_shift (complex_type* rightmoving, complex_type* leftmoving, cudaTextureObject_t tex, uint elements_per_thread, uint stride)
{
  const unsigned int index_left = (blockIdx.x + 1) % gridDim.x;

  for (int i = 0; i < elements_per_thread; ++i){

    float buffer = rightmoving[i].x;
    //Read the transpose texture
    float2 texture = tex2DLayered<float2>(tex, threadIdx.y, threadIdx.x + i * stride, blockIdx.x);
    
    rightmoving[i].x = rightmoving[i].x * texture.x - rightmoving[i].y * texture.y;
    rightmoving[i].y = buffer * texture.y + rightmoving[i].y * texture.x;

    buffer = leftmoving[i].x;
    texture = tex2DLayered<float2>(tex, threadIdx.y, threadIdx.x + i * stride, index_left);

    leftmoving[i].x = leftmoving[i].x * texture.x - leftmoving[i].y * texture.y;
    leftmoving[i].y = buffer * texture.y + leftmoving[i].y * texture.x;    
  }
}
   
template<typename complex_type>
__device__ void transpose (complex_type* rightmoving, complex_type* leftmoving, uint elements_per_thread, uint size, uint stride, complex_type* shared_mem)
{
  for (int n = 0; n < elements_per_thread; ++n){

    unsigned int i = (threadIdx.x + threadIdx.y) % stride;
    unsigned int x = ((i + n * stride) % size - threadIdx.y + size) % size;   
    unsigned int pos = (x - threadIdx.x + size) % size / stride;
      
    shared_mem[threadIdx.y + size * i].x = rightmoving[pos].x;
    shared_mem[threadIdx.y + size * i].y = rightmoving[pos].y;

    __syncthreads();

    rightmoving[pos].x = shared_mem[x + size * i].x;
    rightmoving[pos].y = shared_mem[x + size * i].y;

  }
    //again for leftmoving
  for (int n = 0; n < elements_per_thread; ++n){

    unsigned int i = (threadIdx.x + threadIdx.y) % stride;
    unsigned int x = ((i + n * stride) % size - threadIdx.y + size) % size;   
    unsigned int pos = (x - threadIdx.x + size) % size / stride;    

    shared_mem[threadIdx.y + size * i].x = leftmoving[pos].x;
    shared_mem[threadIdx.y + size * i].y = leftmoving[pos].y;

    __syncthreads();

    leftmoving[pos].x = shared_mem[x + size * i].x;
    leftmoving[pos].y = shared_mem[x + size * i].y;    
  }
}


template<class FFT, class iFFT>
__launch_bounds__(FFT::max_threads_per_block)
__global__ void propagator_kernel(typename FFT::value_type * fields, float2 * out, cudaTextureObject_t k, cudaTextureObject_t surface)
{
  using complex_type = typename FFT::value_type;

  complex_type thread_rightmoving[FFT::storage_size];
  complex_type thread_leftmoving[FFT::storage_size];  

  //some convinient variables
  const unsigned int stride = size_of<FFT>::value / FFT::elements_per_thread;
  const unsigned int size = size_of<FFT>::value;
  const unsigned int disk_offset = size_of<FFT>::value * size_of<FFT>::value;
  const unsigned int dir_offset = disk_offset * gridDim.x;
  const unsigned int index_left = (blockIdx.x + 1) % gridDim.x;
  
  for (int i = 0; i < FFT::elements_per_thread; ++i){
    //rightmoving components
    thread_rightmoving[i].x = fields[threadIdx.y * size + threadIdx.x + i * stride + disk_offset * blockIdx.x].x;
    thread_rightmoving[i].y = fields[threadIdx.y * size + threadIdx.x + i * stride + disk_offset * blockIdx.x].y;

    //leftmoving components

    thread_leftmoving[i].x = fields[threadIdx.y * size + threadIdx.x + i * stride + disk_offset * index_left + dir_offset].x;
    thread_leftmoving[i].y = fields[threadIdx.y * size + threadIdx.x + i * stride + disk_offset * index_left + dir_offset].y;    
  };

  extern __shared__ complex_type shared_mem[];

  //execute first round (dimension) of Fourier transforms
  FFT().execute(thread_rightmoving, shared_mem);
  FFT().execute(thread_leftmoving, shared_mem);

  //Transpose thread_data matrix

  transpose<complex_type>(thread_rightmoving, thread_leftmoving, FFT::elements_per_thread, size, stride, shared_mem);
  

  FFT().execute(thread_rightmoving, shared_mem);
  FFT().execute(thread_leftmoving, shared_mem);
  
  //elementwise product with the texture
  phase_shift<complex_type>(thread_rightmoving, thread_leftmoving, k, FFT::elements_per_thread, stride);

  //inverse Transform
  iFFT().execute(thread_rightmoving, shared_mem);
  iFFT().execute(thread_leftmoving, shared_mem);

  //transpose
  transpose<complex_type>(thread_rightmoving, thread_leftmoving, FFT::elements_per_thread, size, stride, shared_mem);

  iFFT().execute(thread_rightmoving, shared_mem);
  iFFT().execute(thread_leftmoving, shared_mem);

  
  //signal to Eout
  if (blockIdx.x == gridDim.x - 1){
    for (int i = 0; i < FFT::elements_per_thread; ++i){
      //rightmoving components
      out[threadIdx.y * size + threadIdx.x + i * stride].x += (1 + r[gridDim.x]) * thread_rightmoving[i].x;
      out[threadIdx.y * size + threadIdx.x + i * stride].y += (1 + r[gridDim.x]) * thread_rightmoving[i].y;
    }
  }

  //surface phase shift
  phase_shift<complex_type>(thread_rightmoving, thread_leftmoving, surface, FFT::elements_per_thread, stride);

  //reflection and transmission and write back to global memory
  if (blockIdx.x == gridDim.x - 1){
    //special cases for first / last boundary
    for (int i = 0; i < FFT::elements_per_thread; ++i){
      fields[threadIdx.y * size + threadIdx.x + i * stride + disk_offset * index_left].x
	= thread_leftmoving[i].x * -1 * r[index_left];
      fields[threadIdx.y * size + threadIdx.x + i * stride + disk_offset * index_left].y
	= thread_leftmoving[i].y * -1 * r[index_left];

      //leftmoving components reflect on the last boundary
      fields[threadIdx.y * size + threadIdx.x + i * stride + disk_offset * blockIdx.x + dir_offset].x
	= thread_rightmoving[i].x * r[gridDim.x];
      fields[threadIdx.y * size + threadIdx.x + i * stride + disk_offset * blockIdx.x + dir_offset].y
	= thread_rightmoving[i].y * r[gridDim.x];
    }
  }
  else {
    for (int i = 0; i < FFT::elements_per_thread; ++i){
      //rightmoving components
      fields[threadIdx.y * size + threadIdx.x + i * stride + disk_offset * index_left].x
	= thread_leftmoving[i].x * -1 * r[index_left] + thread_rightmoving[i].x * (1 + r[index_left]);
      fields[threadIdx.y * size + threadIdx.x + i * stride + disk_offset * index_left].y
	= thread_leftmoving[i].y * -1 * r[index_left] + thread_rightmoving[i].y * (1 + r[index_left]);

      //leftmoving components
      fields[threadIdx.y * size + threadIdx.x + i * stride + disk_offset * blockIdx.x + dir_offset].x
	= thread_rightmoving[i].x * r[index_left] + thread_leftmoving[i].x * (1 - r[index_left]);
      fields[threadIdx.y * size + threadIdx.x + i * stride + disk_offset * blockIdx.x + dir_offset].y
	= thread_rightmoving[i].y * r[index_left] + thread_leftmoving[i].y * (1 - r[index_left]);
    };
  }
}


float2 * dancer (int nmax, boundarys & bdry, coords & cord, float lambda)
{
    const int size2 = cord.size * cord.size;
    const int n_regions = 2 * bdry.disknumber + 1;

    //Managed memory can be accessed by both the host and device
    float2 * field;
    cudaMallocManaged(&field, size2 * n_regions * 2 * sizeof(float2));

    // maybe just use regular device memory for Eout
    float2 * Eout;
    cudaMallocManaged(&Eout, size2 * sizeof(float2));

    /* This holds the pre-generated phase shifts (generation is handeld by the constructor)*/
    phase_shifts phases(bdry, cord, lambda);

    /* initialize axion emissions from disks */
    setup_fields(field, bdry, cord);
    
    using FFT = decltype(Size<Grid_Size>() + Precision<float>() + Type<fft_type::c2c>()
			 + Direction<fft_direction::forward>() + SM<700>()
			 + Block() + ElementsPerThread<Elements_Per_Thread>() + FFTsPerBlock<Grid_Size>());

    using iFFT = decltype(Size<Grid_Size>() + Precision<float>() + Type<fft_type::c2c>()
			 + Direction<fft_direction::inverse>() + SM<700>()
			 + Block() + ElementsPerThread<Elements_Per_Thread>() + FFTsPerBlock<Grid_Size>());

    using complex_type = typename FFT::value_type;

    dim3 block_dim(size_of<FFT>::value/FFT::elements_per_thread, FFT::ffts_per_block, 1);
    
    dim3 grid_dim(n_regions, 1, 1);
    
    const unsigned int matr_transpose_mem = sizeof(complex_type) * size2 / FFT::elements_per_thread;
    
    const unsigned int shared_mem_size = max(FFT::shared_memory_size, matr_transpose_mem);
    
    //printf("%i FFT, %i transpose \n\n", FFT::shared_memory_size, matr_transpose_mem);

    //copy reflectivitys to constant memory
    cudaMemcpyToSymbol(r, bdry.r, (n_regions + 1) * sizeof(float));

    /* main loop */ 
    for(int n = 0; n <= nmax; n++)
    {
      propagator_kernel<FFT, iFFT><<<grid_dim, block_dim, shared_mem_size>>>((complex_type *) field, Eout, phases.k_tex, phases.surface_tex);

      //probably not necessary
      cudaDeviceSynchronize();
    }
    
    cudaDeviceSynchronize();
    
    //print_field(&field[0], cord);
    
    //copy Eout back to host;
    float2 * output;
    output = new float2 [size2];

    //apply surface phase shift for the last region on the host
    for (int i = 0; i < size2; ++i){
      output[i].x = Eout[i].x * phases.surface_out[i].x - Eout[i].y * phases.surface_out[i].y;

      output[i].y = Eout[i].x * phases.surface_out[i].y + Eout[i].y * phases.surface_out[i].x;
    }
    
    cudaFree(field);
    cudaFree(Eout);
    
    return output;
}

int main ()
{
  //initialise all constants and disk settings
  float freq = 10e9;
  float c = 299792458;

  float lambda = c / freq;

  
  /* generates coordinate grid in real / frequency space with arg number of elements and values between -0.5 and 0.5 (position space) */
  coords coordinates(0.5, 64);

  /* Like in julia this holds all relevant data of the experimental setup including the disk radius */
  boundarys bdry(disknumber, 0.4, coordinates);

  /* perform some tests to ensure all array acssesses are in bounds and results are as expected */
  float2 * Eout;

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();


  Eout = dancer(iterations, bdry, coordinates, lambda);


  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  //print_field(Eout, coordinates);

  std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[µs]" << std::endl;


}
