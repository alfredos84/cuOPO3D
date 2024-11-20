#include <sys/time.h>

#ifndef _COMMON_H
#define _COMMON_H

#define CHECK(call)                                                            \
{                                                                              \
    const cudaError_t error = call;                                            \
    if (error != cudaSuccess)                                                  \
    {                                                                          \
        fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__);                 \
        fprintf(stderr, "code: %d, reason: %s\n", error,                       \
                cudaGetErrorString(error));                                    \
        exit(1);                                                               \
    }                                                                          \
}

#define CHECK_CUBLAS(call)                                                     \
{                                                                              \
    cublasStatus_t err;                                                        \
    if ((err = (call)) != CUBLAS_STATUS_SUCCESS)                               \
    {                                                                          \
        fprintf(stderr, "Got CUBLAS error %d at %s:%d\n", err, __FILE__,       \
                __LINE__);                                                     \
        exit(1);                                                               \
    }                                                                          \
}

#define CHECK_CURAND(call)                                                     \
{                                                                              \
    curandStatus_t err;                                                        \
    if ((err = (call)) != CURAND_STATUS_SUCCESS)                               \
    {                                                                          \
        fprintf(stderr, "Got CURAND error %d at %s:%d\n", err, __FILE__,       \
                __LINE__);                                                     \
        exit(1);                                                               \
    }                                                                          \
}

#define CHECK_CUFFT(call)                                                      \
{                                                                              \
    cufftResult err;                                                           \
    if ( (err = (call)) != CUFFT_SUCCESS)                                      \
    {                                                                          \
        fprintf(stderr, "Got CUFFT error %d at %s:%d\n", err, __FILE__,        \
                __LINE__);                                                     \
        exit(1);                                                               \
    }                                                                          \
}

#define CHECK_CUSPARSE(call)                                                   \
{                                                                              \
    cusparseStatus_t err;                                                      \
    if ((err = (call)) != CUSPARSE_STATUS_SUCCESS)                             \
    {                                                                          \
        fprintf(stderr, "Got error %d at %s:%d\n", err, __FILE__, __LINE__);   \
        cudaError_t cuda_err = cudaGetLastError();                             \
        if (cuda_err != cudaSuccess)                                           \
        {                                                                      \
            fprintf(stderr, "  CUDA error \"%s\" also detected\n",             \
                    cudaGetErrorString(cuda_err));                             \
        }                                                                      \
        exit(1);                                                               \
    }                                                                          \
}


void printLineOnScreen(){
    std::cout << "*-----------------------------------------------------------------*" << std::endl;
   return ;
}


double seconds()
{   // Time meter
    static auto start = std::chrono::high_resolution_clock::now();
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = now - start;
    return elapsed.count();
}


void TimingCode( double iElaps)
{   // Measure the simulation runtime
	if ( iElaps < 60. ){std::cout << "\n\nTime elapsed " <<  iElaps << " seconds\n\n " << std::endl;}
	else if ( iElaps >= 60. and iElaps < 3600. ){std::cout << "\n\nTime elapsed " <<  iElaps/60 << " minutes\n\n " << std::endl;}
	else{std::cout << "\n\nTime elapsed " <<  iElaps/3600 << " hours\n\n " << std::endl;}
	
	return ;
}


real_t measuredTime(real_t tic) { return (seconds() - tic); }


void RemainingTime(uint rt, real_t toc)
{   // Measure remaining time for `rt` round trips and `toc` elapsed time from current simulation
    real_t iElaps = toc*(NRT-rt);
	if ( iElaps < 60. ){std::cout << "\n\nEstimated remaining time " <<  iElaps << " seconds\n\n " << std::endl;}
	else if ( iElaps >= 60. and iElaps < 3600. ){std::cout << "\n\nEstimated remaining time " <<  iElaps/60 << " minutes\n\n " << std::endl;}
	else{std::cout << "\n\nEstimated remaining time " <<  iElaps/3600 << " hours\n\n " << std::endl;}
	printLineOnScreen();

	return ;
}


// Linear spacing for time vectors
template<typename T>
void linspace( T& Vec, real_t xmin, real_t xmax)
{   // Function for filling a vector with linear-spaced values
    uint size = Vec.size();
	for (uint32_t i = 0; i < Vec.size(); i++)
		Vec[i] = xmin + i * (xmax - xmin)/(size-1);
	
	return ;
}


template<typename T>
T linspace( real_t xmin, real_t xmax, uint size)
{   // Function for generating linear-spaced values
    T Vec (size);
	for (uint32_t i = 0; i < Vec.size(); i++)
		Vec[i] = xmin + i * (xmax - xmin)/(size-1);
	
	return Vec ;
}


struct PrintComplex
{
    __host__ __device__
    void operator()(const complex_t& x) const
    {
        printf("(%f, %f)\n", x.x, x.y);
    }
};


struct PrintReal
{
    __host__ __device__
    void operator()(const real_t& x) const
    {
        printf("%f\n", x);
    }
};


void runOPO_status( uint r, uint print_each )
{   // Show status during simulations in cavity
    if( ((r%print_each == 0) or (r == NRT-1)) and (r != 0) )
        std::cout << "\n\n# Round trip: " << r << " - Completed " << r*100/NRT << "%" << "\t\r" << std::flush;
    return ;
}


template<typename T>
void printVarOnScreen( std::string text,  T var)
{   // Print value of variable on screen
    std::cout << text << var << std::endl;
    return ;
}


__host__ __device__  inline uint32_t IDX( uint x,  uint y, uint z)
{   // Function for indexing matrices
	return ((z*(NX*NY))+(y*NX)+x);
}


__host__ __device__  inline uint32_t IDX( uint x,  uint y )
{   // Function for indexing matrices in 2D <-> 3D conversion
	return ( ( y * NT ) + x );
}


__host__ __device__  inline uint32_t IDX( uint t)
{   // Function for indexing vectors
	return t ;
}


template<typename T>
void printVectorOnScreen( T Vec )
{   // Print real vector on screen
    rVech_t Vh = Vec;
    for (uint32_t i = 0; i < Vh.size(); i++)
            std::cout << std::setprecision(5) << Vh[i] << "\n";

    return ;
}


template<typename T>
void printComplexVectorOnScreen( T Vec )
{   // Print complex vector on screen   
    cVech_t Vh = Vec;
    for (uint32_t i = 0; i < Vh.size(); i++)
            std::cout << std::setprecision(5) << Vh[i].x << std::showpos << " " << Vh[i].y << " j" << std::endl;

    return ;
}


template<typename T>
void printRealPartMatrixOnScreen( T Vec, uint slice )
{   // Print real matrix on screen
    cVech_t Vh = Vec;
    for (uint32_t i = 0; i < NX; i++) {
        for (uint32_t j = 0; j < NY; j++) {
            std::cout << std::setprecision(20) << Vh[IDX(i,j,slice)].x << " ";
        }
        std::cout << std::endl; 
    }
    return ;
}


void is_gpu_available()
{   // Function to verify the availability of a GPU
    int dev = cudaGetDevice(&dev);	// Set up device (GPU)
    // if(dev>0)
    {cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    printLineOnScreen();
    std::cout << "\n\nUsing Device " << dev << ": GPU " << deviceProp.name << "\n\n" << std::endl;
    
    cudaSetDevice(dev);}
    return ;
}


inline void printGrid()
{   // Print set grid on screen
    printLineOnScreen();
    std::cout << "\nSpatio-temporal and block grids:" << std::endl;
    std::cout << "        ---> (NT, NX, NY, NZ)   = (" << NT << ", "<< NX << ", " << NY << ", " << NZ << ")" << std::endl;
    std::cout << "        ---> (BLKT, BLKX, BLKY) = (" << BLKT << ", "<< BLKX << ", " << BLKY << ")" << std::endl;
    std::cout << "        ---> SIZE = NT*NX*NY    = "  << NT * NX * NY << " = 2^(" << log2(SIZE) << ")\n" << std::endl;
}

#endif // _COMMON_H