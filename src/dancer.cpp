#include <iostream>
#include <cstdio>
#include <cmath>
#include <math.h> 
/* try using FP_FAST_FMA -fused multipy add- for faster and more accurate multipy add operations */
#include <complex.h>
#include <fftw3.h>
#include <chrono>

using namespace std;
/* 
C++ implementation of the dancer algorithm. 

    TODO simplify by using some gobal constants (these should not impact perfomance as long as they are constant)
    TODO implement some of the special cases from the dancer function

*/

/*
Some usefull strucs for the experimental setup (boundarys), coordiante system (coords), 
fields on the plates (fields) and the pre-generated phase shifts (phase_shifts)
*/
struct coords 
{ /* Coordinate grids in real and momentum space. square and symmetric around 0 */
    int size;

    double X_step;
    double X_max;

    double kX_step;
    double kX_max;

    coords (double max, double step) /* Constructor equivalent to init_coords */
    {
        X_step = step;
        X_max = max;

        size = 2 * X_max / X_step;

        kX_step = M_PI / max;
        kX_max = kX_step * (size - 1) / 2;
    }

    coords (double max, int len) /* alternative constructor using the number of elements along one dimension */
    {
        X_step = 2 * max / (len - 1);
        X_max = max;

        size = len;

        kX_step = M_PI / max;
        kX_max = kX_step * (size - 1) / 2;
    }

    /* functions for acessing values like in arrays */
    double x(int i)
    {
        return i * X_step - X_max;
    }

    double k(int i)
    {
        return i * kX_step - kX_max;
    }

    /* no dynamic memory so no destrucor needed */
};

struct boundarys 
{ /* boundary object can be created at compile time */ 

    int disknumber;
    double * distance; /* dynamic arrays because size not known at compile time */
    double radius;

    complex<double> * eps;
    double * r;

    double * tilt_x;
    double * tilt_y;
    complex<double> * surfaces; /* 3D array, linear in memory */

    boundarys (int disk_n, double rad, coords &cord) /* basic constructor for n disks. Maybe calculate r from eps */
    {
        disknumber = disk_n;
        radius = rad;

        /* dynamically allocate space for 1D arrays*/
        distance = new double[2 * disk_n + 1]; 
        eps = new complex<double>[2 * disk_n + 1];
        r = new double[2 * disk_n + 2];

        /* 
        initialize titlts to zero
        The standart zero initialisation new double [...]() seems to be causing issues */
        tilt_x = new double[2 * disk_n + 1]();
        tilt_y = new double[2 * disk_n + 1]();


        surfaces = new complex<double> [(2 * disk_n + 1) * cord.size * cord.size]();

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

struct fields
{
    /* Fftw needs data as linear array. Maybe convert all dynamic arrays to linear memory for speed */ 
    fftw_complex * lin;

    /* Add pointers to the field on every disk for convinient acssess */
    fftw_complex *** ptr;

    fields(boundarys & bdry, coords & cord)
    {
        int n_regions = 2 * bdry.disknumber + 1;
 
        lin = (fftw_complex*) fftw_malloc(2 * n_regions * cord.size * cord.size * sizeof(fftw_complex));


        ptr = new fftw_complex ** [2];

        ptr[0] = new fftw_complex * [n_regions];
        ptr[1] = new fftw_complex * [n_regions];

        for(int i = 0; i < n_regions; i++)
        {
            ptr[0][i] = &lin[i * cord.size * cord.size];
            ptr[1][i] = ptr[0][i] + (n_regions * cord.size * cord.size);
        }
    }

    ~fields()
    {
        fftw_free(lin);
        delete[] ptr[0];
        delete[] ptr[1];
        delete[] ptr;
        /* ptr[i][j] already points to memory that was de-allocated by fftw_free(lin) (I hope) */
    }
};

struct phase_shifts
{
    complex<double> * k;
    complex<double> * surface;
    complex<double> * surface_out;

    phase_shifts (boundarys &bdry, coords &cord, double lambda)
    {
        int n_regions = 2 * bdry.disknumber + 1;
        int size2 = cord.size * cord.size;

        k = new complex<double> [n_regions * size2];
        surface = new complex<double> [n_regions * size2];
        surface_out = new complex<double> [size2];

        int x_shift;
        int y_shift;
        int lin_index;

        complex<double> k0;
        complex<double> t_y;
        complex<double> t_x;
        complex<double> s;
        complex<double> a;

        double disk_cutoff;

        for(int z = 0; z < n_regions; z++)
        {
            k0 = 2 * M_PI * sqrt(bdry.eps[z]) / lambda;

            for(int y = 0; y < cord.size; y++)
            {
                t_y = exp(-1i * k0 * bdry.tilt_y[z] * cord.x(y));

                for(int x = 0; x < cord.size; x++)
                {
                    lin_index = x + cord.size * (y + cord.size * z);
                    y_shift = (y + cord.size / 2) % cord.size;
                    x_shift = (x + cord.size / 2) % cord.size;


                    a = conj(sqrt( pow(k0,2) - pow(cord.k(x_shift), 2) - pow(cord.k(y_shift), 2)));

                    k[lin_index] = exp(-1i * bdry.distance[z] * a);


                    disk_cutoff = (pow(bdry.radius,2) > pow(cord.x(x),2) + pow(cord.x(y),2));

                    t_x = exp(-1i * k0 * bdry.tilt_x[z] * cord.x(x));
                    s = exp(-1i * k0 * bdry.surfaces[lin_index]);

                    surface[lin_index] = (1. / size2) * disk_cutoff * t_x * t_y * s;

                    if(z == n_regions - 1)
                    {
                        surface_out[x + y * cord.size] = (1. / size2) * t_x * t_y * s;
                    }
                }
            }
        }
    }
    /* destructor */
    ~phase_shifts() 
    {
        delete[] k;
        delete[] surface;
        delete[] surface_out;
    }
};

/* note that for the phase shifts, these are mutiplied each by 1/size^2 as fftw resusts in unnormalized ffts. */

/* simple function used for printing field slices to console for use in testing / debugging */
void print_field (fftw_complex * field, coords & cord)
{
    int n_x = cord.size;

    for(int i = 0; i < n_x; i++)
    {
        for(int j = 0; j < n_x; j++)
        {
            printf( "%7.3f %7.3fi   ", field[n_x * i + j][0], field[n_x * i + j][1] );
        }
        cout << endl;
    }
    cout << endl;
}

void print_prop (complex<double> * field, coords & cord, int zmax)
{
    int n_x = cord.size;
    for(int z = 0; z < zmax; z++)
    {
        for(int i = 0; i < n_x; i++)
        {
            for(int j = 0; j < n_x; j++)
            {
                printf( "%7.3f %7.3fi   ", real(field[n_x * (n_x * z + i) + j]), imag(field[n_x * (n_x * z + i) + j]) );
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

void set_field (fftw_complex * field, coords & cord, boundarys & bdry, complex<double> val)
{

    for(int i = 0; i < (cord.size * cord.size); i++)
    {
        /* Using boolean arithmetic may be faster than using if statments (see branchless programming) */
        bool disk_cutoff = pow(bdry.radius, 2) > (pow(cord.x(i % cord.size), 2) + pow(cord.x(i / cord.size), 2));
        field[i][0] = disk_cutoff * real(val);
        field[i][1] = disk_cutoff * imag(val);
    }
}
    
complex<double> ax (complex<double> eps_i, complex<double> eps_m)
{
    return sqrt(eps_m) * (1. - eps_i / eps_m) / (eps_i * sqrt(eps_m) + eps_m * sqrt(eps_i));
}

void setup_fields (fftw_complex *** fields, boundarys & bdry, coords & cord) 
{
    int n_regions = 2 * bdry.disknumber + 1;

    set_field(fields[0][0], cord, bdry, -1.);
    set_field(fields[1][0], cord, bdry, ax(bdry.eps[1], bdry.eps[0]));

    /* standart cases */
    for(int n = 1; n < n_regions - 1; n++)
    {
        set_field(fields[0][n], cord, bdry, ax(bdry.eps[n-1], bdry.eps[n]));
        set_field(fields[1][n], cord, bdry, ax(bdry.eps[n+1], bdry.eps[n]));
    }

    set_field(fields[0][n_regions-1], cord, bdry, ax(bdry.eps[n_regions-2], bdry.eps[n_regions-1]));
    set_field(fields[1][n_regions-1], cord, bdry, ax(1., bdry.eps[n_regions-1]));
}

/*
Some usefull functions for multipling the phase-shift on to the fields
*/
void mult_phase(fftw_complex * in, fftw_complex * out, complex<double> * phase, int len)
{
    for(int i = 0; i < len; i++)
    {
        double a = real(phase[i]);
        double b = imag(phase[i]);
        double c = in[i][0]; /* neccessary if in = out because the first operation overrides the real part needed to calcualte the imagnary part */ 

        out[i][0] = a*c - b*in[i][1];
        out[i][1] = b*c + a*in[i][1];

        c = in[i + len][0];//repeat on the other direction to apply phase shift to the entire array.

        out[i + len][0] = a*c - b*in[i + len][1];
        out[i + len][1] = b*c + a*in[i + len][1];
    } 
}

void mult(fftw_complex * in, fftw_complex * out, double x, int len, bool dbl)
{
    if(dbl)
    {
        for(int i = 0; i < len; i++)
        {
            out[i][0] *= x;
            out[i][1] *= x;
            out[i + len][0] *= x;
            out[i + len][1] *= x;
        }
    }
    else
    {
        for(int i = 0; i < len; i++)
        {
            out[i][0] *= x;
            out[i][1] *= x;
            out[i + len][0] *= x;
            out[i + len][1] *= x;
        }
    }  
}
/*
The heart of it
*/
complex<double> * dancer (int nmax, boundarys & bdry, coords & cord, double lambda)
{

    int disk_index = 2 * bdry.disknumber;
    /* 
    field.lin saves the fields as fftw_complex in a single, linear array. fftw_malloc ensures good alingment for SIMD operations.
    Individual fields for a single region can be accessed by field.ptr[direction][region]
    The order of these is reversed compared to the julia version to optimise memory access along a single direction.
    */
    fields field(bdry, cord);

    fftw_complex * field_arrived_left; /* pre allocated space for the intermedeary array during reflection and transmission */
    field_arrived_left = (fftw_complex*) fftw_malloc(disk_index * pow(cord.size, 2) * sizeof(fftw_complex));

    complex<double> * Eout;
    Eout = (complex<double>*) fftw_malloc(pow(cord.size, 2) * sizeof(complex<double>));

    /* This holds the pre-generated phase shifts (generation is handeld by the constructor)*/
    phase_shifts phases(bdry, cord, lambda);

    /* 
    Create fft plans. 
    fft_plan_many_dft is the fftw equivalent of cuda batched fft with many of the same arguments.
    See http://fftw.org/fftw3_doc/Advanced-Complex-DFTs.html#Advanced-Complex-DFTs for an explanation of the arguments.
    Plan creation must happen before initialisation of data becaudse FFTW_MEASURE may overwrite data.
    Using FFTW_MEASURE is much slower in plan creation but can be faster in execution than FFTW_ESTIMATE. 
    More flags can be found in the FFTW documentation (such as FFTW_PATIENT for even faster execution).
    Alternativley one may use FFTW_EXAUSTIVE once to determine a near optimal plan for a given size and save it to an exteranl file
    (see FFTW wisdom) for future use.
    Note that these transforms are un-normalised.
    */
    const int n[2] = {cord.size, cord.size};

    fftw_plan plan = fftw_plan_many_dft(2, n, 2 * (2 * bdry.disknumber + 1), 
        field.lin, NULL, 1, pow(cord.size, 2),
        field.lin, NULL, 1, pow(cord.size, 2),
        FFTW_FORWARD, FFTW_MEASURE);

    fftw_plan inv_plan = fftw_plan_many_dft(2, n, 2 * (2 * bdry.disknumber + 1), 
        field.lin, NULL, 1, pow(cord.size, 2),
        field.lin, NULL, 1, pow(cord.size, 2),
        FFTW_BACKWARD, FFTW_MEASURE);

    /* initialize axion emissions from disks */
    setup_fields(field.ptr, bdry, cord);

    bool rightmoving = false; /* rightmoving starts out as 0 */

    /* main loop */ 
    int size2 = cord.size * cord.size;
    int length = (2 * bdry.disknumber + 1) * size2;

    for(int n = 0; n <= nmax; n++)
    {
        /*cout << "initialization" << endl;
        print_field(field.ptr[rightmoving][1], cord);*/

        /* propagator(field.lin, Eout, phases, length, plan, inv_plan); */
        fftw_execute(plan);

        /*cout << "fft" << endl;
        print_field(field.ptr[rightmoving][1], cord); */       

        mult_phase(field.lin, field.lin, phases.k, length);

        /*cout << "prop" << endl;
        print_field(field.ptr[rightmoving][1], cord);*/

        fftw_execute(inv_plan);

        /*cout << "ifft" << endl;
        print_field(field.ptr[rightmoving][1], cord);*/    


        for(int i = 0; i < size2; i++)
        {
            Eout[i] += (1 + bdry.r[disk_index + 1]) * phases.surface_out[i] *
                (field.ptr[rightmoving][disk_index][i][0] + 1i * field.ptr[rightmoving][disk_index][i][1]) ;
        }

        mult_phase(field.lin, field.lin, phases.surface, length);

        /*cout << "surface" << endl;
        print_field(field.ptr[rightmoving][1], cord);*/

        /* 
            Potentially oveload operators to make calculations involving fftw_complex less tedious.
            Copy to field_arrived_left reflect rightmoving components (exept for one) and ass transmissions
        */

        for(int i = 0; i <= length-size2; i++)/* copy to fields arrived left */
        {
            field_arrived_left[i][0] = field.ptr[rightmoving][0][i][0];
            field_arrived_left[i][1] = field.ptr[rightmoving][0][i][1];
        }
        for(int i = 0; i <= length; i++)/* reflect rightmoving components */
        {
            field.ptr[rightmoving][0][i][0] *= bdry.r[1 + i / size2];
            field.ptr[rightmoving][0][i][1] *= bdry.r[1 + i / size2];
        }
        for(int i = 0; i <= length-size2; i++)/* combine transmitted leftmoving with rightmoving */
        {
            field.ptr[rightmoving][0][i][0] += (1 - bdry.r[1 + i / size2]) * field.ptr[!rightmoving][1][i][0];
            field.ptr[rightmoving][0][i][1] += (1 - bdry.r[1 + i / size2]) * field.ptr[!rightmoving][1][i][1];
        }
        for(int i = 0; i <= length; i++)/* reflect leftmoving */
        {
            field.ptr[!rightmoving][0][i][0] *= -bdry.r[i / size2];
            field.ptr[!rightmoving][0][i][1] *= -bdry.r[i / size2];
        }
        for(int i = 0; i <= length-size2; i++)/* combine transmitted rightmoving with leftmoving */
        {
            field.ptr[!rightmoving][1][i][0] += (1 + bdry.r[1 + i / size2]) * field_arrived_left[i][0];
            field.ptr[!rightmoving][1][i][1] += (1 + bdry.r[1 + i / size2]) * field_arrived_left[i][1];
        }
        /* switch directions */
        rightmoving = !rightmoving;
    } 

    /* Destroy plans */
    fftw_destroy_plan(plan);
    fftw_destroy_plan(inv_plan);

    fftw_free(field_arrived_left); 

    return Eout;
}

int main ()
{
    /* initialise all constants and disk settings */
    double freq = 10e9;
    double c = 299792458;

    int disk_number = 10;

    double lambda = c / freq;

    /* generates coordinate grid in real / frequency space with arg number of elements and values between -0.5 and 0.5 (position space) */
    coords coordinates(0.5, 128);

    /* Like in julia this holds all relevant data of the experimental setup including the disk radius */
    boundarys bdry(disk_number, 0.4, coordinates);

    /* perform some tests to ensure all array acssesses are in bounds and results are as expected */
    complex<double> * Eout; 

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    Eout = dancer(1000, bdry, coordinates, lambda);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[µs]" << std::endl;
}
