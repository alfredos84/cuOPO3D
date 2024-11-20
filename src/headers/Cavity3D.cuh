/*---------------------------------------------------------------------------*/
// * This file contains functions to -----
/*---------------------------------------------------------------------------*/



#ifndef _CAVITYCUH
#define _CAVITYCUH


/** Add phase (phase) and mirror losses (R) after a single-pass */
__global__ void roundTripChange(complex_t *A, complex_t *aux, real_t R, real_t phase)
{
	
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint32_t idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY){
			aux[IDX(idx,idy,idt)] = (sqrtf(R) * CpxExp(phase)) * A[IDX(idx,idy,idt)];
		}
		if( idx < NX and idy < NY){
			A[IDX(idx,idy,idt)] = aux[IDX(idx,idy,idt)];
		}
	}
	
	return ;
}


template<typename Crystal>
class Cavity
{	// Difine the class Cavity....

public:	
	real_t Rs, Ri, Lcav, Lfs, trt, fsr, dt; 	// cavity 
	real_t gamma;      		       	// GDD comp
	real_t beta, fpm;			  	// EOM	
	bool gdd, satabs, eom;    		// presence of intracavity elements
	Crystal *Cr;

	Cavity( Crystal *_Cr, real_t _Lcav, real_t _Rs, real_t _Ri) : Cr(_Cr), Lcav(_Lcav), Rs(_Rs), Ri(_Ri)
	{	// Constructor
		trt = ( Lcav + (Cr->Lcr)*((Cr->ns) - 1) )/C;
		dt	= trt/(NT-1);
		fsr = 1./trt;
		Lfs = Lcav-(Cr->Lcr); 
		printLineOnScreen();
		printf("\nInstance of the class Cavity.\n");
	}

	~Cavity()
	{	// Destructor
		printf("Cavity Destructor.\n");
	}

	// Methods definition
	void applyReflectivity(cVecd_t &A, real_t Reflectivity, real_t phase);
	void getCavityProp();
	// void setEOM(real_t _beta, real_t _fpm);
	void setGDDCompensation (real_t _gamma);
};


// Methods declaration

template<typename Crystal>
void Cavity<Crystal>::applyReflectivity(cVecd_t &A, real_t R, real_t phase)
{
	dim3 block2D(BLKX, BLKY);	dim3 grid2D((NX+BLKX-1)/BLKX, (NY+BLKY-1)/BLKY);
	
	complex_t * aux_ptr;
	CHECK(cudaMalloc((void **)&aux_ptr, A.size()*sizeof(complex_t) ));
	complex_t * A_ptr = thrust::raw_pointer_cast(A.data());
	roundTripChange<<<grid2D,block2D>>>(A_ptr, aux_ptr, R, phase);
	CHECK(cudaDeviceSynchronize());
	CHECK(cudaFree(aux_ptr));
	
	return ;
}


template<typename Crystal>
void Cavity<Crystal>::setGDDCompensation (real_t _gamma)
{
	// set gamma between 0 and 1
	this->gdd = true;
	this->gamma = _gamma;
	return ;
}


// template<typename Crystal>
// void Cavity<Crystal>::setEOM(real_t _beta, real_t _fpm)
// {	
// 	this->eom = true; 
// 	this->fpm = _fpm; this->beta = _beta; 
// 	std::cout << "Set intracavity EOM with \u03B2 = " << this->beta << "- fpm = " << this->fpm << std::endl;
// 	return ;
// }


template<typename Crystal>
void Cavity<Crystal>::getCavityProp ()
{
	std::cout << "Set OPO cavity:" << std::endl;
	std::cout << "        ---> Cavity length         = " << (this->Lcav)*1e-6 << " m" << std::endl;
	std::cout << "        ---> Free-space length     = " << (this->Lfs)*1e-6 << " m" << std::endl;
	std::cout << "        ---> Round trip time       = " << this->trt << " ps" << std::endl;
	std::cout << "        ---> Free spectral range   = " << (this->fsr)*1e3 << " GHz" << std::endl;
	std::cout << "        ---> Signal reflectivity   = " << (this->Rs)*100 << " %" << std::endl;
	std::cout << "        ---> Idler reflectivity    = " << (this->Ri)*100 << " %" << std::endl;
	if(this->gdd){
		std::cout << "        ---> Set " << static_cast<int>(100*this->gamma)<< "% Intracavity GDD compensation." << std::endl;}
	std::cout << "\n" << std::endl;
	return ;
}


#endif // -> #ifdef _CAVITYCUH