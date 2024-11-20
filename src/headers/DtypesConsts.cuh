#ifndef _DTYPESCONSTSCUH
#define _DTYPESCONSTSCUH

 
//  * Complex data type: a set of datatypes are
//  * defined to make the code more readable.
//  *
//  * Definitions for numbers
//  * real_t    : datatype for real numbers
//  * complex_t : datatype for complex numbers
//  * 
//  * Definitions for vectors:
//  * 
//  * rVech_t   : real vector host
//  * rVecd_t   : real vector device
//  * cVech_t   : complex vector host
//  * cVecd_t   : complex vector device

using real_t = float;
using complex_t = cufftComplex;
using rVech_t = thrust::host_vector<real_t>;
using rVecd_t = thrust::device_vector<real_t>;
using cVech_t = thrust::host_vector<complex_t>;
using cVecd_t = thrust::device_vector<complex_t>;


// Define grid size

#ifdef DIFFRACTION
const uint32_t NX     = 1 << 7;		// size in X direction
const uint32_t NY     = NX; 		// size in Y direction
#ifdef DISPERSION
const uint32_t NT     = 1 << 11;	// vector size
#else
const uint32_t NT     = 1;	// vector size
#endif
const uint32_t SIZE   = NX*NY*NT;	// size of A(x,y,t) for fixed z
#else
const uint32_t NT     = 1 << 8;	// vector size
const uint32_t NX     = 1;		    // size in X direction
const uint32_t NY     = 1;		    // size in Y direction
const uint32_t SIZE   = NT;	        // size of A(x,y,t) for fixed z
#endif
const uint32_t NZ     = 100;		// size discretization


// Define block size for CUDA kernels
const uint32_t BLKT   = 1 << 5;		// block dimensions for kernels
const uint32_t BLKX   = 1 << 4;		// block dimensions for kernels
const uint32_t BLKY   = BLKX;		// block dimensions for kernels


// Define global constants
const uint32_t NRT    = 10000;                       // number of round trips    
const real_t PI       = 3.14159265358979323846;		// pi
const real_t C        = 299792458*1E6/1E12;			// speed of ligth in vacuum [um/ps]
const real_t EPS0     = 8.8541878128E-12*1E12/1E6;	// vacuum pertivity [W.ps/V²μm] 


// Memory size for vectors and matrices
const size_t nBytes1Dr = sizeof(real_t) * NT;		// real 1D
const size_t nBytes1Dc = sizeof(complex_t) * NT;		// complex 1D
#ifdef DIFFRACTION
const size_t nBytes2Dr = sizeof(real_t) * NX * NY;		// real 2D 
const size_t nBytes2Dc = sizeof(complex_t) * NX * NY;	// complex 2D
const size_t nBytes3Dr = sizeof(real_t) * SIZE;	        // real 3D
const size_t nBytes3Dc = sizeof(complex_t) * SIZE;	    // complex 3D
#endif


#endif // -> #ifdef _DTYPESCONSTS
