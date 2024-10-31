/*---------------------------------------------------------------------------*/
// * This file contains functions to -----
/*---------------------------------------------------------------------------*/



#ifndef _CAVITYCUH
#define _CAVITYCUH


/** Add phase (phase) and mirror losses (R) after a single-pass */
__global__ void kernelRoundTripChange(complex_t *A, complex_t *aux, real_t R, real_t phase)
{
	uint idt = threadIdx.x + blockDim.x*blockIdx.x;
	
	if( idt < NT ){aux[idt] = (sqrtf(R) * CpxExp(phase)) * A[idt];}
	if( idt < NT ){A[idt] = aux[idt];}

	return ;
}


template<typename Crystal>
class Cavity
{	// Difine the class Cavity.

public:	
	real_t R, Lcav, trt, fsr; // cavity 
	real_t gamma;             // GDD comp
	real_t beta, fpm;		  // EOM	
	bool gdd, satabs, eom;    // presence of intracavity elements
	Crystal *Cr;

	Cavity( Crystal *_Cr, real_t _Lcav, real_t _R ) : Cr(_Cr), Lcav(_Lcav), R(_R)
	{	// Constructor
		trt = ( Lcav + (Cr->Lcr)*((Cr->ns) - 1) )/C;
		fsr = 1./trt;
		printLineOnScreen();
		printf("\nInstance of the class Cavity.\n");
	}

	~Cavity()
	{	// Destructor
		printf("Cavity Destructor.\n");
	}

	// Methods definition
	void applyReflectivity(cVecd_t &A, real_t phase);
	void getCavityProp();
	void setEOM(real_t _beta, real_t _fpm);
	void setGDD (real_t _gamma);
};


// Methods declaration
template<typename Crystal>
void Cavity<Crystal>::applyReflectivity(cVecd_t &A, real_t phase)
{
	dim3 block1D(BLKT);	dim3 grid1D((NT+BLKT-1)/BLKT);
	
	complex_t * aux_ptr;
	CHECK(cudaMalloc((void **)&aux_ptr, nBytes1Dc ));
	complex_t * A_ptr = thrust::raw_pointer_cast(A.data());
	kernelRoundTripChange<<<grid1D,block1D>>>(A_ptr, aux_ptr, this->R, phase);
	CHECK(cudaDeviceSynchronize());
	CHECK(cudaFree(aux_ptr));
	
	return ;
}


template<typename Crystal>
void Cavity<Crystal>::setEOM(real_t _beta, real_t _fpm)
{	
	this->eom = true; 
	this->fpm = _fpm; this->beta = _beta; 
	std::cout << "Set intracavity EOM:" << std::endl;
	std::cout << "        ---> \u03B2 = " << this->beta << std::endl;
	std::cout << "        ---> fpm = " << (this->fsr)*1e3 << " GHz" << std::endl;
	return ;
}


template<typename Crystal>
void Cavity<Crystal>::setGDD (real_t _gamma)
{
	// set gamma between 0 and 1
	this->gdd = true;
	this->gamma = _gamma;
	std::cout << "\nSet " << static_cast<int>(100*this->gamma)<< "% intracavity GDD compensation.\n" << std::endl;
	return ;
}


template<typename Crystal>
void Cavity<Crystal>::getCavityProp ()
{
	std::cout << "Set OPO cavity:" << std::endl;
	std::cout << "        ---> Cavity length         = " << (this->Lcav)*1e-6 << " m" << std::endl;
	std::cout << "        ---> Round trip time       = " << this->trt << " ps" << std::endl;
	std::cout << "        ---> Free spectral range   = " << (this->fsr)*1e3 << " GHz" << std::endl;
	std::cout << "        ---> Reflectivity          = " << (this->R)*100 << " %\n" << std::endl;

	return ;
}


#endif // -> #ifdef _CAVITYCUH