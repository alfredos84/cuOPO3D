/*---------------------------------------------------------------------------*/
// * This file contains functions to solve the Split-Step Fourier method (SSMF)
// * needed to calculate the electric fields evolution through the nonlinear Cr.
// * 
// * In particular, this file should be used when only three equation describes the 
// * problem, i.e., sum or difference frequency generation (SFG or DFG).
// *
// * For any specific process, please check the form of the Couple wave equations
// * in the first function called dAdz(). 
/*---------------------------------------------------------------------------*/



#ifndef _SOLVERCUH
#define _SOLVERCUH

// Scales a vector after Fourier transforms in the time-frequency domaing
// (CUFFT_INVERSE mode)
__global__ void cuFFT1DScale(complex_t *Ax)
{		
	
	real_t size = static_cast<real_t>(NT);
	
	uint idt = threadIdx.x + blockDim.x*blockIdx.x;
			
	if(idt < NT){Ax[idt] = Ax[idt] / size;}
	
	return ;
	
}


// Scales a vector after Fourier transforms in the time-frequency domaing
// (CUFFT_INVERSE mode)
__global__ void cuFFT1DScaleAllFields(complex_t *Ap, complex_t *As, complex_t *Ai)
{		
	
	real_t size = static_cast<real_t>(NT);
	
	uint idt = threadIdx.x + blockDim.x*blockIdx.x;
			
	if(idt < NT){Ap[idt] = Ap[idt] / size; As[idt] = As[idt] / size; Ai[idt] = Ai[idt] / size;}
	
	return ;
	
}


/** This function compensates the GVD after a single-pass */
__global__ void kernelAddGDD(complex_t *A, complex_t *aux, real_t *w, real_t GDD)
{
	
	uint idt = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (idt < NT){aux[idt] = A[idt] * CpxExp( 0.5*w[idt]*w[idt]*GDD );}
	if (idt < NT){A[idt] = aux[idt];}
	
	return ;
}


/** Applies an electro optical modulator to an electric field after one round trip. */
__global__ void KernelEOM(complex_t *A, complex_t *aux, real_t beta, real_t fpm, real_t *t)
{
	
	uint idt = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (idt < NT){aux[idt] = A[idt] * CpxExp(beta*sinf(2*PI*fpm*t[idt]));}
	if (idt < NT){A[idt] = aux[idt];}

	return ;
}


// Product of complex numbers in GPU for Dispersion propagator
__global__ void kernelDispersionPropagatorProduct ( complex_t *Awp_propagated, complex_t *eiLz_p, complex_t *Awp, 
													complex_t *Aws_propagated, complex_t *eiLz_s, complex_t *Aws, 
													complex_t *Awi_propagated, complex_t *eiLz_i, complex_t *Awi )
{	
	
	uint idt = threadIdx.x + blockDim.x*blockIdx.x;

	if( idt < NT ){
		Awp_propagated[idt] = eiLz_p[idt] * Awp[idt];
		Aws_propagated[idt] = eiLz_s[idt] * Aws[idt];
		Awi_propagated[idt] = eiLz_i[idt] * Awi[idt];
	}
	
	return ;
	
}

/** Computes the nonlinear part: dA/dz=i.κ.Ax.Ay.exp(i.Δk.L) and saves the result in dAx (x represents the different fields) */
__global__ void dAdz( complex_t *dAp, complex_t *dAs,  complex_t *dAi, complex_t *Ap, complex_t *As, complex_t *Ai, 
	real_t kappa_p, real_t kappa_s, real_t kappa_i, real_t dk, real_t z )
{
	complex_t Im; Im.x = 0; Im.y = 1;
	
	uint idt = threadIdx.x + blockDim.x*blockIdx.x;

	if( idt < SIZE ){	
		dAp[idt]  = Im * kappa_p * As[idt] * Ai[idt] * CpxExp(-0.*dk*z) ;
		dAs[idt]  = Im * kappa_s * Ap[idt] * CpxConj(Ai[idt]) * CpxExp(+dk*0.*z);
		dAi[idt]  = Im * kappa_i * Ap[idt] * CpxConj(As[idt]) * CpxExp(+dk*0.*z);
	}

	
	return ;
}


/** Computes a linear combination Ax + s.kx and saves the result in aux_x */
__global__ void LinealCombination( complex_t *auxp, complex_t *auxs, complex_t *auxi,
								   complex_t *Ap, complex_t *As, complex_t *Ai, 
								   complex_t *kp, complex_t *ks, complex_t *ki, real_t s )
{
	
	uint idt = threadIdx.x + blockDim.x*blockIdx.x;
	
	if( idt < SIZE ){		
		auxp[idt] = Ap[idt] + kp[idt] * s;
		auxs[idt] = As[idt] + ks[idt] * s;
		auxi[idt] = Ai[idt] + ki[idt] * s;
	}


	return ;
}


/** This kernel computes the final sum after appling the Rounge-Kutta algorithm */
__global__ void rk4(complex_t *Ap, complex_t *As,  complex_t *Ai, 
					complex_t *k1p, complex_t *k1s, complex_t *k1i, 
					complex_t *k2p, complex_t *k2s, complex_t *k2i, 
					complex_t *k3p, complex_t *k3s, complex_t *k3i,
					complex_t *k4p, complex_t *k4s, complex_t *k4i,
					real_t dz)
{
	
	uint idt = threadIdx.x + blockDim.x*blockIdx.x;
	
	if( idt < SIZE ){	
		Ap[idt] = Ap[idt] + (k1p[idt] + 2*k2p[idt] + 2*k3p[idt] + k4p[idt]) * dz / 6.0f;
		As[idt] = As[idt] + (k1s[idt] + 2*k2s[idt] + 2*k3s[idt] + k4s[idt]) * dz / 6.0f;
		Ai[idt] = Ai[idt] + (k1i[idt] + 2*k2i[idt] + 2*k3i[idt] + k4i[idt]) * dz / 6.0f;
	}
	
	return ;
}



// Difine the class Solver
template<typename Crystal>
class Solver
{	
public:	

	cVecd_t k1p, k2p, k3p, k4p;
	cVecd_t k1s, k2s, k3s, k4s;
	cVecd_t k1i, k2i, k3i, k4i;
	cVecd_t auxp, auxs, auxi;

	Crystal *Cr;
	EFields<Crystal> *A;
	Cavity<Crystal> *Cav;

	Solver(Crystal *_Cr, EFields<Crystal> *_A) : Cr(_Cr), A(_A)
	{	// Constructor

		k1p.resize(SIZE); k2p.resize(SIZE); k3p.resize(SIZE); k4p.resize(SIZE);
		k1s.resize(SIZE); k2s.resize(SIZE); k3s.resize(SIZE); k4s.resize(SIZE);
		k1i.resize(SIZE); k2i.resize(SIZE); k3i.resize(SIZE); k4i.resize(SIZE);
		auxp.resize(SIZE); auxs.resize(SIZE); auxi.resize(SIZE);

		printLineOnScreen();
		printf("\nInstance of the class Solver.\n");
		std::cout << "        ---> Number of crystal slices: " << NZ << std::endl;
	}


	Solver(Crystal *_Cr, Cavity<Crystal> *_Cav, EFields<Crystal> *_A) : Cr(_Cr), Cav(_Cav), A(_A)
	{	// Constructor

		k1p.resize(SIZE); k2p.resize(SIZE); k3p.resize(SIZE); k4p.resize(SIZE);
		k1s.resize(SIZE); k2s.resize(SIZE); k3s.resize(SIZE); k4s.resize(SIZE);
		k1i.resize(SIZE); k2i.resize(SIZE); k3i.resize(SIZE); k4i.resize(SIZE);
		auxp.resize(SIZE); auxs.resize(SIZE); auxi.resize(SIZE);
		
		printLineOnScreen();
		printf("\nInstance of the class Solver.\n");
		std::cout << "        ---> Number of roundtrips: " << NRT << std::endl;
		std::cout << "        ---> Number of crystal slices: " << NZ << std::endl;
	}
	
	~Solver(){ printf("Solver Destructor.\n"); }

	// Methods definition
	void EOM ( cVecd_t &Ax );
	void GDD( cVecd_t &Ax );
	bool checkNaN( uint r );
	void solverRK4( real_t z );	
	void dispersion();
	void SSFM ();
	void runOPO(const std::vector<int>& save_roundtrips);
	void runOPO();
	void runSinglePass();

};


// Methods declaration

// Check NaN subrutine: 
// 				----> r: 	roundtrip
// 				----> each: How many steps do I want to check?
template<typename Crystal>
inline bool Solver<Crystal>::checkNaN ( uint r )
{
	bool check = false;
	cVech_t test_vec = A->As;	real_t test_nan =  test_vec[0].x;
	if(isnan(test_nan)){
		std::cout << "Simulation failed! \nError in NaN check at round-trip = " << r << std::endl;
		check = true;
	}	
	
	return check;
}


template<typename Crystal>
void Solver<Crystal>::EOM ( cVecd_t &Ax )
{
	// Parameters for kernels 1D	
	dim3 block1D(BLKT);	dim3 grid1D((NT+BLKT-1)/BLKT);

	complex_t *aux_ptr;
	CHECK( cudaMalloc((void **)&aux_ptr, Ax.size()*sizeof(complex_t)) );
	real_t* t_ptr = thrust::raw_pointer_cast( A->t.data() );
	complex_t* A_ptr = thrust::raw_pointer_cast( Ax.data() );
			
	KernelEOM<<<grid1D,block1D>>>( A_ptr, aux_ptr, Cav->beta, Cav->fpm, t_ptr );
	CHECK(cudaDeviceSynchronize());
	CHECK(cudaFree(aux_ptr));

	return ;
}


template<typename Crystal>
void Solver<Crystal>::GDD( cVecd_t &Ax )
{
	// Parameters for kernels 1D	
	dim3 block1D(BLKT);	dim3 grid1D((NT+BLKT-1)/BLKT);

	// Set plan for cuFFT 1D//
	cufftHandle plan1D; cufftPlan1d(&plan1D, NT, CUFFT_C2C, 1);

	complex_t *A_ptr = thrust::raw_pointer_cast(Ax.data());
	real_t* w_ptr = thrust::raw_pointer_cast(A->w.data());
	
	complex_t *Aw_ptr , *aux_ptr;
	CHECK(cudaMalloc((void **)&Aw_ptr, Ax.size()*sizeof(complex_t) ));
	CHECK(cudaMalloc((void **)&aux_ptr, Ax.size()*sizeof(complex_t) ));

	cufftExecC2C(plan1D, (complex_t *)A_ptr, (complex_t *)Aw_ptr, CUFFT_INVERSE);
	CHECK(cudaDeviceSynchronize());

	cuFFT1DScale<<<grid1D,block1D>>>(Aw_ptr);
	CHECK(cudaDeviceSynchronize());
	
	real_t GDD = -(Cav->gamma)*(Cr->b2s)*(Cr->Lcr);
	kernelAddGDD<<<grid1D,block1D>>>(Aw_ptr, aux_ptr, w_ptr, GDD);
	CHECK(cudaDeviceSynchronize());
	
	cufftExecC2C(plan1D, (complex_t *)Aw_ptr, (complex_t *)A_ptr, CUFFT_FORWARD);
	CHECK(cudaDeviceSynchronize());
	
	cufftDestroy(plan1D);
	CHECK(cudaFree(Aw_ptr)); CHECK(cudaFree(aux_ptr));

	return ;
}


template<typename Crystal>
void Solver<Crystal>::solverRK4( real_t z )
{	
	// A->Applies the Fourh-order Runge-Kutta (RK4) Method
	// with fixed step size dz.
	
	// Parameters for kernels
	dim3 block1D(BLKT);	dim3 grid1D((NT+BLKT-1)/BLKT);

	// Define pointers to use them in kernels
	complex_t * Ap_ptr  = thrust::raw_pointer_cast(A->Ap.data());
	complex_t * As_ptr  = thrust::raw_pointer_cast(A->As.data());
	complex_t * Ai_ptr  = thrust::raw_pointer_cast(A->Ai.data());

	complex_t * k1p_ptr = thrust::raw_pointer_cast(this->k1p.data());
	complex_t * k2p_ptr = thrust::raw_pointer_cast(this->k2p.data());
	complex_t * k3p_ptr = thrust::raw_pointer_cast(this->k3p.data());
	complex_t * k4p_ptr = thrust::raw_pointer_cast(this->k4p.data());
	
	complex_t * k1s_ptr = thrust::raw_pointer_cast(this->k1s.data());
	complex_t * k2s_ptr = thrust::raw_pointer_cast(this->k2s.data());
	complex_t * k3s_ptr = thrust::raw_pointer_cast(this->k3s.data());
	complex_t * k4s_ptr = thrust::raw_pointer_cast(this->k4s.data());

	complex_t * k1i_ptr = thrust::raw_pointer_cast(this->k1i.data());
	complex_t * k2i_ptr = thrust::raw_pointer_cast(this->k2i.data());
	complex_t * k3i_ptr = thrust::raw_pointer_cast(this->k3i.data());
	complex_t * k4i_ptr = thrust::raw_pointer_cast(this->k4i.data());
	
	complex_t * auxp_ptr = thrust::raw_pointer_cast(this->auxp.data());
	complex_t * auxs_ptr = thrust::raw_pointer_cast(this->auxs.data());
	complex_t * auxi_ptr = thrust::raw_pointer_cast(this->auxi.data());
	
	real_t dz = Cr->dz;	real_t dk = A->dk; 
	real_t kappa_p = A->kappa_p, kappa_s = A->kappa_s, kappa_i = A->kappa_i;

	//k1 = dAdz(kappas,dk,z,A)
	dAdz<<<grid1D,block1D>>>(	k1p_ptr, k1s_ptr, k1i_ptr,
								Ap_ptr, As_ptr, Ai_ptr,
								kappa_p, kappa_s, kappa_i, dk,
								z ); CHECK(cudaDeviceSynchronize()); 

	//k2 = dAdz(kappas,dk,z+dz/2,A+k1/2) -> aux = A+k1/2
	LinealCombination<<<grid1D,block1D>>>( 	auxp_ptr, auxs_ptr, auxi_ptr,
											Ap_ptr, As_ptr, Ai_ptr,
											k1p_ptr, k1s_ptr, k1i_ptr,
											0.5f );	CHECK(cudaDeviceSynchronize());   
											
	dAdz<<<grid1D,block1D>>>(	k2p_ptr, k2s_ptr, k2i_ptr,
								auxp_ptr, auxs_ptr, auxi_ptr,
								kappa_p, kappa_s, kappa_i, dk,
								z+dz/4 ); CHECK(cudaDeviceSynchronize());

	// k3 = dAdz(kappas,dk,z+dz/2,A+k2/2)
	LinealCombination<<<grid1D,block1D>>>(	auxp_ptr, auxs_ptr, auxi_ptr,
											Ap_ptr, As_ptr, Ai_ptr, 
											k2p_ptr, k2s_ptr, k2i_ptr, 
											0.5f ); CHECK(cudaDeviceSynchronize());   
	
	dAdz<<<grid1D,block1D>>>(	k3p_ptr, k3s_ptr, k3i_ptr,
								auxp_ptr, auxs_ptr, auxi_ptr,
								kappa_p, kappa_s, kappa_i, dk,
								z+dz/4 ); CHECK(cudaDeviceSynchronize());

	// k4 = dAdz(kappas,dk,z+dz,A+k3)
	LinealCombination<<<grid1D,block1D>>>(	auxp_ptr, auxs_ptr, auxi_ptr, 
											Ap_ptr, As_ptr, Ai_ptr, 
											k3p_ptr, k3s_ptr, k3i_ptr, 
											1.0f );	CHECK(cudaDeviceSynchronize());

	dAdz<<<grid1D,block1D>>>(	k4p_ptr, k4s_ptr, k4i_ptr,
								auxp_ptr, auxs_ptr, auxi_ptr,
								kappa_p, kappa_s, kappa_i, dk,
								z+dz/2 ); CHECK(cudaDeviceSynchronize());

	// A = A + (k1+k2+k3+k4)/6
	rk4<<<grid1D,block1D>>>(	Ap_ptr, As_ptr, Ai_ptr, 
								k1p_ptr, k1s_ptr, k1i_ptr, 
								k2p_ptr, k2s_ptr, k2i_ptr, 
								k3p_ptr, k3s_ptr, k3i_ptr,
								k4p_ptr, k4s_ptr, k4i_ptr, 
								dz/2 );	CHECK(cudaDeviceSynchronize());

	return ;
	
}


// Applies the dispersion term to the electric fields
template<typename Crystal>
void Solver<Crystal>::dispersion ()
{
	// Parameters for kernels 1D	
	dim3 block1D(BLKT);	dim3 grid1D((NT+BLKT-1)/BLKT);

	// Set plan for cuFFT 1D//
	cufftHandle plan1D; cufftPlan1d(&plan1D, NT, CUFFT_C2C, 1);

	complex_t *Aux_p; CHECK(cudaMalloc((void **)&Aux_p, nBytes1Dc));
	complex_t *Aux_s; CHECK(cudaMalloc((void **)&Aux_s, nBytes1Dc));
	complex_t *Aux_i; CHECK(cudaMalloc((void **)&Aux_i, nBytes1Dc));
	
	complex_t *Ap_ptr = thrust::raw_pointer_cast(A->Ap.data());
	complex_t *Awp_ptr = thrust::raw_pointer_cast(A->Awp.data());
	complex_t *eiLz_p_ptr = thrust::raw_pointer_cast(A->eiLz_p.data());

	complex_t *As_ptr = thrust::raw_pointer_cast(A->As.data());
	complex_t *Aws_ptr = thrust::raw_pointer_cast(A->Aws.data());
	complex_t *eiLz_s_ptr = thrust::raw_pointer_cast(A->eiLz_s.data());

	complex_t *Ai_ptr = thrust::raw_pointer_cast(A->Ai.data());
	complex_t *Awi_ptr = thrust::raw_pointer_cast(A->Awi.data());
	complex_t *eiLz_i_ptr = thrust::raw_pointer_cast(A->eiLz_i.data());


	cufftExecC2C(plan1D, (complex_t *)Ap_ptr, (complex_t *)Awp_ptr, CUFFT_INVERSE);
	CHECK(cudaDeviceSynchronize());
	cufftExecC2C(plan1D, (complex_t *)As_ptr, (complex_t *)Aws_ptr, CUFFT_INVERSE);
	CHECK(cudaDeviceSynchronize());	
	cufftExecC2C(plan1D, (complex_t *)Ai_ptr, (complex_t *)Awi_ptr, CUFFT_INVERSE);
	CHECK(cudaDeviceSynchronize());

	cuFFT1DScaleAllFields<<<grid1D, block1D>>>(Awp_ptr, Aws_ptr, Awi_ptr);
	CHECK(cudaDeviceSynchronize());	

	kernelDispersionPropagatorProduct<<<grid1D, block1D>>> (Aux_p, eiLz_p_ptr, Awp_ptr, 
															Aux_s, eiLz_s_ptr, Aws_ptr,
															Aux_i, eiLz_i_ptr, Awi_ptr);
	CHECK(cudaDeviceSynchronize());

	CHECK(cudaMemcpy(Awp_ptr, Aux_p, nBytes1Dc, cudaMemcpyDeviceToDevice));	
	CHECK(cudaMemcpy(Aws_ptr, Aux_s, nBytes1Dc, cudaMemcpyDeviceToDevice));	
	CHECK(cudaMemcpy(Awi_ptr, Aux_i, nBytes1Dc, cudaMemcpyDeviceToDevice));	

	cufftExecC2C(plan1D, (complex_t *)Awp_ptr, (complex_t *)Ap_ptr, CUFFT_FORWARD);
	CHECK(cudaDeviceSynchronize());
	cufftExecC2C(plan1D, (complex_t *)Aws_ptr, (complex_t *)As_ptr, CUFFT_FORWARD);
	CHECK(cudaDeviceSynchronize());
	cufftExecC2C(plan1D, (complex_t *)Awi_ptr, (complex_t *)Ai_ptr, CUFFT_FORWARD);
	CHECK(cudaDeviceSynchronize());


	CHECK(cudaFree(Aux_p));	CHECK(cudaFree(Aux_s));	CHECK(cudaFree(Aux_i));
	
	cufftDestroy(plan1D);
		
	return ;
}


template<typename Crystal>
void Solver<Crystal>::SSFM ( )
{
	for (real_t z = 0; z <= Cr->Lcr; z+=(Cr->dz))
	{	
		solverRK4(z); // RK4 in dz/2
		dispersion(); // Dispersion in dz
		solverRK4(z); // RK4 in dz/2
	}

	return ;
}


template<typename Crystal>
void Solver<Crystal>::runOPO(const std::vector<int>& save_roundtrips)
{
	complex_t * Ap_ptr  = thrust::raw_pointer_cast(A->Ap.data());
	complex_t * Api_ptr  = thrust::raw_pointer_cast(A->Api.data());

	std::unordered_set<int> save_set(save_roundtrips.begin(), save_roundtrips.end());

	real_t iStart = seconds();	// Timing code
	
	for( uint r = 0; r < NRT; r++ )
	{
		
		runOPO_status (r, NRT/10); // print simulation status on screen
		SSFM(); // compute couple-wave equations


		if(Cav->eom){EOM(A->As); EOM(A->Ai);}			// intracavity elements: EOM
		if(Cav->gdd){GDD(A->As); GDD(A->Ai);}			// GDD
				
		// Apply boundary conditions to resonant fields	
		if(std::get<1>(A->resFields)){Cav->applyReflectivity(A->As, 0.0f);}
		if(std::get<2>(A->resFields)){Cav->applyReflectivity(A->Ai, 0.0f);}
		A->Ap = A->Api;
		
		if ( r%NRT/10 == 0 ){if(checkNaN(r) ){r = NRT;}}
		if (r == 0){RemainingTime( timingOneRoundtrip(iStart));}

	}

	real_t iElaps = seconds() - iStart;	// finish timing
	TimingCode( iElaps); // print time
	
	return ;

}


// Overloaded runOPO() without saving data
template<typename Crystal>
void Solver<Crystal>::runOPO() {runOPO(std::vector<int>());}


template<typename Crystal>
void Solver<Crystal>::runSinglePass()
{
	double iStart = seconds();	// Timing code

	SSFM();

	double iElaps = seconds() - iStart;	// finish timing
	TimingCode( iElaps); // print time
	
	return ;
}



#endif // -> #ifdef _SOLVERCUH