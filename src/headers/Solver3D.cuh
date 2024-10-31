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
__global__ void cuFFT1DScaleAllFields(complex_t *Ap, complex_t *As, complex_t *Ai)
{		
	
	real_t size = static_cast<real_t>(NT);
	
	uint idt = threadIdx.x + blockDim.x*blockIdx.x;
			
	if(idt < NT){Ap[idt] /= size; As[idt] /= size; Ai[idt] /= size;}
	
	return ;
	
}


// Scales a vector after Fourier transforms in the space-reciprocal space
// (CUFFT_INVERSE mode)	
__global__ void cuFFT2DScaleAllFields( complex_t *Ap, complex_t *As, complex_t *Ai )
{	
	
	real_t size = static_cast<real_t>(NX*NY);
	
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;
		
	for (uint idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY){
			Ap[IDX(idx,idy,idt)] /= size;
			As[IDX(idx,idy,idt)] /= size;
			Ai[IDX(idx,idy,idt)] /= size;
		}
	}
	
	return ;
	
}


// Scales a vector after Fourier transforms in the space-reciprocal space
// (CUFFT_INVERSE mode)	
__global__ void cuFFT2DScaleSignalIdler( complex_t *As, complex_t *Ai )
{	
	
	real_t size = static_cast<real_t>(NX*NY);
	
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;
		
	for (uint idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY){
			As[IDX(idx,idy,idt)] /= size;
			Ai[IDX(idx,idy,idt)] /= size;
		}
	}
	
	return ;
	
}


#ifdef DISPERSION
// This kernel set the linear propagator operators useful in dispersion calculations
__global__ void dispersionPropagators ( complex_t *eiLz_p, complex_t *eiLz_s, complex_t *eiLz_i, 
										real_t vp, real_t vs, real_t vi, 
										real_t b2p, real_t b2s, real_t b2i,
										real_t alpha_crp, real_t alpha_crs, real_t alpha_cri,
										real_t dz, real_t *w )
{	
			
	uint idw = threadIdx.x + blockDim.x*blockIdx.x;

	if( idw < NT ){
		eiLz_p[idw] = CpxExp(dz*w[idw]*((1/vp-1/vi)+0.5*w[idw]*b2p )) * expf(-0.5*alpha_crp*dz);
		eiLz_s[idw] = CpxExp(dz*w[idw]*((1/vs-1/vs)+0.5*w[idw]*b2s )) * expf(-0.5*alpha_crs*dz);
		eiLz_i[idw] = CpxExp(dz*w[idw]*((1/vi-1/vs)+0.5*w[idw]*b2i )) * expf(-0.5*alpha_cri*dz); 
	}

	return ;	
}


// Product of complex numbers in GPU for Dispersion propagator
__global__ void kernelDispersionPropagatorProduct ( complex_t *Awp_propagated, complex_t *eiLz_p, complex_t *Awp,
													complex_t *Aws_propagated, complex_t *eiLz_s, complex_t *Aws,
													complex_t *Awi_propagated, complex_t *eiLz_i, complex_t *Awi  )
{	
	
	uint idt = threadIdx.x + blockDim.x*blockIdx.x;

	if( idt < NT ){
		Awp_propagated[idt] = eiLz_p[idt] * Awp[idt];
		Aws_propagated[idt] = eiLz_s[idt] * Aws[idt];
		Awi_propagated[idt] = eiLz_i[idt] * Awi[idt];
	}
	
	return ;
	
}


// Copy a slice from a 3D tensor to a 1D vector of complex numbers in GPU
__global__ void kernelGetVectorT (	complex_t *Aux_Ap, complex_t *Ap_ptr, 
									complex_t *Aux_As, complex_t *As_ptr, 
									complex_t *Aux_Ai, complex_t *Ai_ptr, 
									uint x, uint y )
{	
	uint idt = threadIdx.x + blockDim.x*blockIdx.x;

	if( idt < NT ){
		Aux_Ap[idt] = Ap_ptr[IDX(x,y,idt)] ;
		Aux_As[idt] = As_ptr[IDX(x,y,idt)] ;
		Aux_Ai[idt] = Ai_ptr[IDX(x,y,idt)] ;
	}
	
	return ;
}


// Copy from a 1D vector to a specific (x,y) of a 3D tensor of complex numbers in GPU
__global__ void kernelSetVectorT (	complex_t *Ap_ptr, complex_t *Aux_Ap, 
									complex_t *As_ptr, complex_t *Aux_As,
									complex_t *Ai_ptr, complex_t *Aux_Ai, 
									uint x, uint y )
{		
	uint idt = threadIdx.x + blockDim.x*blockIdx.x;

	if( idt < NT ){
		Ap_ptr[IDX(x,y,idt)] = Aux_Ap[idt];
		As_ptr[IDX(x,y,idt)] = Aux_As[idt];
		Ai_ptr[IDX(x,y,idt)] = Aux_Ai[idt];
	}
	
	return ;
}
#endif


// This kernel set the beam propagator operators useful in diffraction calculations
__global__ void kernelSetDiffractionPropagator ( complex_t *eiQz_p, complex_t *eiQz_s, complex_t *eiQz_i, 
										real_t uX, real_t uY,
										real_t dx, real_t dy, real_t dz,
										real_t kp, real_t ks, real_t ki,
										real_t ap, real_t as, real_t ai )
{	
	real_t dfX  = 1/dx/NX;	real_t dfY  = 1/dy/NY;
		
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;

	if( idx < NX and idy < NY){
		eiQz_p[IDX(idx,idy,0)] = CpxExp(-dz*(2*powf(PI,2)/kp * ( dfX*dfX*powf(idx - uX,2) + dfY*dfY*powf(idy - uY,2))))*expf(-0.5f*ap*dz); 
		eiQz_s[IDX(idx,idy,0)] = CpxExp(-dz*(2*powf(PI,2)/ks * ( dfX*dfX*powf(idx - uX,2) + dfY*dfY*powf(idy - uY,2))))*expf(-0.5f*as*dz); 
		eiQz_i[IDX(idx,idy,0)] = CpxExp(-dz*(2*powf(PI,2)/ki * ( dfX*dfX*powf(idx - uX,2) + dfY*dfY*powf(idy - uY,2))))*expf(-0.5f*ai*dz); 
	}
	
	
	return ;	
}


// This kernel set the beam propagator operators useful in diffraction calculations
__global__ void kernelSetDiffractionPropagatorInFreeSpace ( complex_t *eiQz_s, complex_t *eiQz_i, 
															real_t uX, real_t uY,
															real_t dx, real_t dy,
															real_t ks, real_t ki,
															real_t length )
{	
	real_t dfX  = 1/dx/NX;	real_t dfY  = 1/dy/NY;
		
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;

	if( idx < NX and idy < NY){ 
		eiQz_s[IDX(idx,idy,0)] = CpxExp(length*(2*powf(PI,2)/(ks) * ( dfX*dfX*powf(idx - uX,2) + dfY*dfY*powf(idy - uY,2)))); 
		eiQz_i[IDX(idx,idy,0)] = CpxExp(length*(2*powf(PI,2)/(ki) * ( dfX*dfX*powf(idx - uX,2) + dfY*dfY*powf(idy - uY,2)))); 
	}
	
	
	return ;	
}


// Product of complex numbers in GPU for Diffraction propagator
__global__ void kernelDiffractionPropagatorProduct (complex_t *AQp_propagated, complex_t *eiQz_p, complex_t *AQp,
													complex_t *AQs_propagated, complex_t *eiQz_s, complex_t *AQs,
													complex_t *AQi_propagated, complex_t *eiQz_i, complex_t *AQi )
{	
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){
			AQp_propagated[IDX(idx,idy,idt)] = eiQz_p[IDX(idx,idy,0)] * AQp[IDX(idx,idy,idt)] ;
			AQs_propagated[IDX(idx,idy,idt)] = eiQz_s[IDX(idx,idy,0)] * AQs[IDX(idx,idy,idt)] ;
			AQi_propagated[IDX(idx,idy,idt)] = eiQz_i[IDX(idx,idy,0)] * AQi[IDX(idx,idy,idt)] ;
		}
	}

	return ;
}


// Product of complex numbers in GPU for Diffraction propagator
__global__ void kernelDiffractionPropagatorProductInFreeSpace (	complex_t *AQs_propagated, complex_t *eiQz_s, complex_t *AQs,
																complex_t *AQi_propagated, complex_t *eiQz_i, complex_t *AQi )
{	
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){
			AQs_propagated[IDX(idx,idy,idt)] = eiQz_s[IDX(idx,idy,0)] * AQs[IDX(idx,idy,idt)] ;
			AQi_propagated[IDX(idx,idy,idt)] = eiQz_i[IDX(idx,idy,0)] * AQi[IDX(idx,idy,idt)] ;
		}
	}

	return ;
}


/** Computes the nonlinear part: dA/dz=i.κ.Ax.Ay.exp(i.Δk.L) and saves the result in dAx (x represents the different fields) */
__global__ void dAdz( complex_t *dAp, complex_t *dAs,  complex_t *dAi, complex_t *Ap, complex_t *As, complex_t *Ai, 
	real_t kappa_p, real_t kappa_s, real_t kappa_i, real_t dk, real_t z )
{
	complex_t Im; Im.x = 0; Im.y = 1;
	
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){	
			dAp[IDX(idx,idy,idt)]  = Im * kappa_p * As[IDX(idx,idy,idt)] * Ai[IDX(idx,idy,idt)] * CpxExp(-0.*dk*z) ;
			dAs[IDX(idx,idy,idt)]  = Im * kappa_s * Ap[IDX(idx,idy,idt)] * CpxConj(Ai[IDX(idx,idy,idt)]) * CpxExp(+dk*0.*z);
			dAi[IDX(idx,idy,idt)]  = Im * kappa_i * Ap[IDX(idx,idy,idt)] * CpxConj(As[IDX(idx,idy,idt)]) * CpxExp(+dk*0.*z);
		}
	}
	
	return ;
}


/** Computes a linear combination Ax + s.kx and saves the result in aux_x */
__global__ void kernelLinealCombination(complex_t *auxp, complex_t *auxs, complex_t *auxi,
								   		complex_t *Ap, complex_t *As, complex_t *Ai, 
								   		complex_t *kp, complex_t *ks, complex_t *ki, real_t s )
{
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){	
			auxp[IDX(idx,idy,idt)] = Ap[IDX(idx,idy,idt)] + kp[IDX(idx,idy,idt)] * s;
			auxs[IDX(idx,idy,idt)] = As[IDX(idx,idy,idt)] + ks[IDX(idx,idy,idt)] * s;
			auxi[IDX(idx,idy,idt)] = Ai[IDX(idx,idy,idt)] + ki[IDX(idx,idy,idt)] * s;
		}
	}

	return ;
}


/** This kernel computes the final sum after appling the Rounge-Kutta algorithm */
__global__ void rk4(complex_t *Ap, complex_t *As,  complex_t *Ai, 
					complex_t *k1p, complex_t *k1s, complex_t *k1i, 
					complex_t *k2p, complex_t *k2s, complex_t *k2i, 
					complex_t *k3p, complex_t *k3s, complex_t *k3i,
					complex_t *k4p, complex_t *k4s, complex_t *k4i, real_t dz )
{
	
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){
			Ap[IDX(idx,idy,idt)] += (k1p[IDX(idx,idy,idt)] + 2.0f*k2p[IDX(idx,idy,idt)] + 2.0f*k3p[IDX(idx,idy,idt)] + k4p[IDX(idx,idy,idt)]) * dz / 6.0f;
			As[IDX(idx,idy,idt)] += (k1s[IDX(idx,idy,idt)] + 2.0f*k2s[IDX(idx,idy,idt)] + 2.0f*k3s[IDX(idx,idy,idt)] + k4s[IDX(idx,idy,idt)]) * dz / 6.0f;
			Ai[IDX(idx,idy,idt)] += (k1i[IDX(idx,idy,idt)] + 2.0f*k2i[IDX(idx,idy,idt)] + 2.0f*k3i[IDX(idx,idy,idt)] + k4i[IDX(idx,idy,idt)]) * dz / 6.0f;
		}
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

	Crystal *Cr;	EFields<Crystal> *A;	Cavity<Crystal> *Cav;

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
	// void EOM ( cVecd_t &Ax );
	// void GDD( cVecd_t &Ax );
	#ifdef DISPERSION
	void setDispersionPropagators();
	void dispersion();
	#endif
	void setDiffractionPropagators();
	void setDiffractionPropagatorsInFreeSpace(real_t length);
	void diffractionInCrystal();
	void diffractionInFreeSpace();
	bool checkCourantStability(real_t dz , real_t kpdxdy2);
	void solverRK4( real_t z );
	void SSFM ();
	void runOPO(const std::vector<int>& save_roundtrips);
	void runOPO();
	void runSinglePass();

};


// Methods declaration
template<typename Crystal>
inline bool Solver<Crystal>::checkCourantStability(real_t dz , real_t kpdxdy)
{
	bool condition;
	
	if(dz < kpdxdy){condition = true;}
	else{condition = false;}

	return condition;
}


#ifdef DISPERSION
// Set vectors for dispersion propagators
template<typename Crystal>
void Solver<Crystal>::setDispersionPropagators()
{	
	// Parameters for kernels 1D
	dim3 block1D(BLKT);	dim3 grid1D((NT+BLKT-1)/BLKT);
	
	real_t *w_ptr = thrust::raw_pointer_cast(A->w.data());
	complex_t *eiLz_p_ptr = thrust::raw_pointer_cast(A->eiLz_p.data());
	complex_t *eiLz_s_ptr = thrust::raw_pointer_cast(A->eiLz_s.data());
	complex_t *eiLz_i_ptr = thrust::raw_pointer_cast(A->eiLz_i.data());

	dispersionPropagators<<<grid1D, block1D>>> (eiLz_p_ptr, eiLz_s_ptr, eiLz_i_ptr, 
												A->vp, A->vs, A->vi, A->b2p, A->b2s, A->b2i,
												A->alpha_crp, A->alpha_crs, A->alpha_cri,
												A->dz, w_ptr );
	CHECK(cudaDeviceSynchronize());
	
	return ;
}


// Applies the dispersion term to the electric fields
template<typename Crystal>
void Solver<Crystal>::dispersion ( )
{
	// Parameters for kernels 1D	
	dim3 block1D(BLKT);	dim3 grid1D((NT+BLKT-1)/BLKT);

	// Set plan for cuFFT 1D//
	cufftHandle plan1D; cufftPlan1d(&plan1D, NT, CUFFT_C2C, 1);

	complex_t *Aux_p; CHECK(cudaMalloc((void **)&Aux_p, nBytes1Dc));
	complex_t *Aux_Ap; CHECK(cudaMalloc((void **)&Aux_Ap, nBytes1Dc));
	complex_t *Aux_Awp; CHECK(cudaMalloc((void **)&Aux_Awp, nBytes1Dc));
	complex_t *Ap_ptr = thrust::raw_pointer_cast(A->Ap.data());
	complex_t *Awp_ptr = thrust::raw_pointer_cast(A->Awp.data());
	complex_t *eiLz_p_ptr = thrust::raw_pointer_cast(A->eiLz_p.data());

	complex_t *Aux_s; CHECK(cudaMalloc((void **)&Aux_s, nBytes1Dc));
	complex_t *Aux_As; CHECK(cudaMalloc((void **)&Aux_As, nBytes1Dc));
	complex_t *Aux_Aws; CHECK(cudaMalloc((void **)&Aux_Aws, nBytes1Dc));
	complex_t *As_ptr = thrust::raw_pointer_cast(A->As.data());
	complex_t *Aws_ptr = thrust::raw_pointer_cast(A->Aws.data());
	complex_t *eiLz_s_ptr = thrust::raw_pointer_cast(A->eiLz_s.data());

	complex_t *Aux_i; CHECK(cudaMalloc((void **)&Aux_i, nBytes1Dc));
	complex_t *Aux_Ai; CHECK(cudaMalloc((void **)&Aux_Ai, nBytes1Dc));
	complex_t *Aux_Awi; CHECK(cudaMalloc((void **)&Aux_Awi, nBytes1Dc));
	complex_t *Ai_ptr = thrust::raw_pointer_cast(A->Ai.data());
	complex_t *Awi_ptr = thrust::raw_pointer_cast(A->Awi.data());
	complex_t *eiLz_i_ptr = thrust::raw_pointer_cast(A->eiLz_i.data());

	for (uint y = 0; y < NY; y++){
		for (uint x = 0; x < NX; x++){

			kernelGetVectorT<<<grid1D, block1D>>> (Aux_Ap, Ap_ptr, Aux_As, As_ptr, Aux_Ai, Ai_ptr, x, y);
			CHECK(cudaDeviceSynchronize());

			cufftExecC2C(plan1D, (complex_t *)Aux_Ap, (complex_t *)Aux_Awp, CUFFT_INVERSE);
			CHECK(cudaDeviceSynchronize());
			cufftExecC2C(plan1D, (complex_t *)Aux_As, (complex_t *)Aux_Aws, CUFFT_INVERSE);
			CHECK(cudaDeviceSynchronize());
			cufftExecC2C(plan1D, (complex_t *)Aux_Ai, (complex_t *)Aux_Awi, CUFFT_INVERSE);
			CHECK(cudaDeviceSynchronize());

			cuFFT1DScaleAllFields<<<grid1D, block1D>>>( Aux_Awp, Aux_Aws, Aux_Awi );
			CHECK(cudaDeviceSynchronize());	

			kernelDispersionPropagatorProduct<<<grid1D, block1D>>> (Aux_p, eiLz_p_ptr, Aux_Awp,
																	Aux_s, eiLz_s_ptr, Aux_Aws,
																	Aux_i, eiLz_i_ptr, Aux_Awi);
			CHECK(cudaDeviceSynchronize());

			CHECK(cudaMemcpy(Aux_Awp, Aux_p, nBytes1Dc, cudaMemcpyDeviceToDevice));	
			CHECK(cudaMemcpy(Aux_Aws, Aux_s, nBytes1Dc, cudaMemcpyDeviceToDevice));	
			CHECK(cudaMemcpy(Aux_Awi, Aux_i, nBytes1Dc, cudaMemcpyDeviceToDevice));	

			cufftExecC2C(plan1D, (complex_t *)Aux_Awp, (complex_t *)Aux_Ap, CUFFT_FORWARD);
			CHECK(cudaDeviceSynchronize());
			cufftExecC2C(plan1D, (complex_t *)Aux_Aws, (complex_t *)Aux_As, CUFFT_FORWARD);
			CHECK(cudaDeviceSynchronize());
			cufftExecC2C(plan1D, (complex_t *)Aux_Awi, (complex_t *)Aux_Ai, CUFFT_FORWARD);
			CHECK(cudaDeviceSynchronize());

			kernelSetVectorT<<<grid1D, block1D>>> (Ap_ptr, Aux_Ap, As_ptr, Aux_As, Ai_ptr, Aux_Ai, x, y);
			CHECK(cudaDeviceSynchronize());

		}
	}

	CHECK(cudaFree(Aux_p)); CHECK(cudaFree(Aux_Ap)); CHECK(cudaFree(Aux_Awp));
	CHECK(cudaFree(Aux_s)); CHECK(cudaFree(Aux_As)); CHECK(cudaFree(Aux_Aws));
	CHECK(cudaFree(Aux_i)); CHECK(cudaFree(Aux_Ai)); CHECK(cudaFree(Aux_Awi));
	
	cufftDestroy(plan1D);
		
	return ;
}
#endif


// Set vectors for diffraction propagators
template<typename Crystal>
void Solver<Crystal>::setDiffractionPropagators()
{	
	
	// Parameters for kernels 2D
	dim3 block2D(BLKX, BLKY);	dim3 grid2D((NX+BLKX-1)/BLKX, (NY+BLKY-1)/BLKY);
	// 		* Diffraction operation ∇²_{xy}: Beam propagator for Pump and Signal
	
	complex_t *eiQz_p_ptr = thrust::raw_pointer_cast(A->eiQz_p.data());
	complex_t *eiQz_s_ptr = thrust::raw_pointer_cast(A->eiQz_s.data());
	complex_t *eiQz_i_ptr = thrust::raw_pointer_cast(A->eiQz_i.data());

	kernelSetDiffractionPropagator<<<grid2D, block2D>>> (eiQz_p_ptr, eiQz_s_ptr, eiQz_i_ptr, 
												0.5f*real_t(NX),  0.5f*real_t(NY), 
												A->dx, A->dy, A->dz,
												A->kp, A->ks, A->ki,
												Cr->alpha_crp, Cr->alpha_crs, Cr->alpha_cri);
	CHECK(cudaDeviceSynchronize());
	
	A->fftShift2D ( A->eiQz_p );	A->fftShift2D ( A->eiQz_s ); A->fftShift2D ( A->eiQz_i );

	return ;
}


// Applies the diffraction term to the electric fields
template<typename Crystal>
void Solver<Crystal>::diffractionInCrystal ()
{	
	// Parameters for kernels 2D
	dim3 block2D(BLKX, BLKY);	dim3 grid2D((NX+BLKX-1)/BLKX, (NY+BLKY-1)/BLKY);
	
	// Set plan for cuFFT 2D//
	cufftHandle plan2D; 
	int batch = NT;	int rank = 2;
	int nRows = NY;	int nCols = NX; int n[2] = {NY, NX};
	int idist = NX * NY; int odist = NX * NY;
	int inembed[] = {NY , NX}; 	int onembed[] = {NY, NX};
	int istride = 1; int ostride = 1;

	cufftPlanMany(&plan2D,  rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, batch);
	
	complex_t *AuxQp; CHECK(cudaMalloc((void **)&AuxQp, nBytes3Dc));
	complex_t *Ap_ptr = thrust::raw_pointer_cast(A->Ap.data());
	complex_t *AQp_ptr = thrust::raw_pointer_cast(A->AQp.data());
	complex_t *eiQz_p_ptr = thrust::raw_pointer_cast(A->eiQz_p.data());

	complex_t *AuxQs; CHECK(cudaMalloc((void **)&AuxQs, nBytes3Dc));
	complex_t *As_ptr = thrust::raw_pointer_cast(A->As.data());
	complex_t *AQs_ptr = thrust::raw_pointer_cast(A->AQs.data());
	complex_t *eiQz_s_ptr = thrust::raw_pointer_cast(A->eiQz_s.data());

	complex_t *AuxQi; CHECK(cudaMalloc((void **)&AuxQi, nBytes3Dc));
	complex_t *Ai_ptr = thrust::raw_pointer_cast(A->Ai.data());
	complex_t *AQi_ptr = thrust::raw_pointer_cast(A->AQi.data());
	complex_t *eiQz_i_ptr = thrust::raw_pointer_cast(A->eiQz_i.data());


	cufftExecC2C(plan2D, (complex_t *)Ap_ptr, (complex_t *)AQp_ptr, CUFFT_FORWARD);
	CHECK(cudaDeviceSynchronize());
	cufftExecC2C(plan2D, (complex_t *)As_ptr, (complex_t *)AQs_ptr, CUFFT_FORWARD);
	CHECK(cudaDeviceSynchronize());
	cufftExecC2C(plan2D, (complex_t *)Ai_ptr, (complex_t *)AQi_ptr, CUFFT_FORWARD);
	CHECK(cudaDeviceSynchronize());
	
	kernelDiffractionPropagatorProduct<<<grid2D, block2D>>> (	AuxQp, eiQz_p_ptr, AQp_ptr,
																AuxQs, eiQz_s_ptr, AQs_ptr,	
																AuxQi, eiQz_i_ptr, AQi_ptr	);
	CHECK(cudaDeviceSynchronize());
	
	CHECK(cudaMemcpy(AQp_ptr, AuxQp, nBytes3Dc, cudaMemcpyDeviceToDevice));	
	CHECK(cudaMemcpy(AQs_ptr, AuxQs, nBytes3Dc, cudaMemcpyDeviceToDevice));
	CHECK(cudaMemcpy(AQi_ptr, AuxQi, nBytes3Dc, cudaMemcpyDeviceToDevice));
	
	cufftExecC2C(plan2D, (complex_t *)AQp_ptr, (complex_t *)Ap_ptr, CUFFT_INVERSE);
	CHECK(cudaDeviceSynchronize());
	cufftExecC2C(plan2D, (complex_t *)AQs_ptr, (complex_t *)As_ptr, CUFFT_INVERSE);
	CHECK(cudaDeviceSynchronize());
	cufftExecC2C(plan2D, (complex_t *)AQi_ptr, (complex_t *)Ai_ptr, CUFFT_INVERSE);
	CHECK(cudaDeviceSynchronize());

	cuFFT2DScaleAllFields<<<grid2D, block2D>>>(Ap_ptr, As_ptr, Ai_ptr);
	CHECK(cudaDeviceSynchronize());		
	
	CHECK(cudaFree(AuxQp)); CHECK(cudaFree(AuxQs)); CHECK(cudaFree(AuxQi));
	
	cufftDestroy(plan2D);
		
	return ;
}


template<typename Crystal>
void Solver<Crystal>::setDiffractionPropagatorsInFreeSpace(real_t length)
{	
	
	// Parameters for kernels 2D
	dim3 block2D(BLKX, BLKY);	dim3 grid2D((NX+BLKX-1)/BLKX, (NY+BLKY-1)/BLKY);
	// 		* Diffraction operation ∇²_{xy}: Beam propagator for Pump and Signal
	
	complex_t *eiQz_p_ptr = thrust::raw_pointer_cast(A->eiQz_p.data());
	complex_t *eiQz_s_ptr = thrust::raw_pointer_cast(A->eiQz_s.data());
	complex_t *eiQz_i_ptr = thrust::raw_pointer_cast(A->eiQz_i.data());

	kernelSetDiffractionPropagatorInFreeSpace<<<grid2D, block2D>>> (eiQz_s_ptr, eiQz_i_ptr, 
																	0.5f*real_t(NX),  0.5f*real_t(NY), 
																	A->dx, A->dy,
																	A->ks, A->ki,
																	length );
	CHECK(cudaDeviceSynchronize());
	
	A->fftShift2D ( A->eiQz_s ); A->fftShift2D ( A->eiQz_i );

	return ;
}


// Applies the diffraction term to the electric fields
template<typename Crystal>
void Solver<Crystal>::diffractionInFreeSpace ()
{	
	// Parameters for kernels 2D
	dim3 block2D(BLKX, BLKY);	dim3 grid2D((NX+BLKX-1)/BLKX, (NY+BLKY-1)/BLKY);
	
	// Set plan for cuFFT 2D//
	cufftHandle plan2D; 
	int batch = NT;	int rank = 2;
	int nRows = NY;	int nCols = NX; int n[2] = {NY, NX};
	int idist = NX * NY; int odist = NX * NY;
	int inembed[] = {NY , NX}; 	int onembed[] = {NY, NX};
	int istride = 1; int ostride = 1;

	cufftPlanMany(&plan2D,  rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, batch);
	// cufftPlan2d(&plan2D, NX, NY, CUFFT_C2C );

	complex_t *AuxQs; CHECK(cudaMalloc((void **)&AuxQs, nBytes3Dc));
	complex_t *As_ptr = thrust::raw_pointer_cast(A->As.data());
	complex_t *AQs_ptr = thrust::raw_pointer_cast(A->AQs.data());
	complex_t *eiQz_s_ptr = thrust::raw_pointer_cast(A->eiQz_s.data());

	complex_t *AuxQi; CHECK(cudaMalloc((void **)&AuxQi, nBytes3Dc));
	complex_t *Ai_ptr = thrust::raw_pointer_cast(A->Ai.data());
	complex_t *AQi_ptr = thrust::raw_pointer_cast(A->AQi.data());
	complex_t *eiQz_i_ptr = thrust::raw_pointer_cast(A->eiQz_i.data());

	cufftExecC2C(plan2D, (complex_t *)As_ptr, (complex_t *)AQs_ptr, CUFFT_FORWARD);
	CHECK(cudaDeviceSynchronize());
	cufftExecC2C(plan2D, (complex_t *)Ai_ptr, (complex_t *)AQi_ptr, CUFFT_FORWARD);
	CHECK(cudaDeviceSynchronize());
	
	kernelDiffractionPropagatorProductInFreeSpace<<<grid2D, block2D>>> (AuxQs, eiQz_s_ptr, AQs_ptr,	
																		AuxQi, eiQz_i_ptr, AQi_ptr);
	CHECK(cudaDeviceSynchronize());
	
	CHECK(cudaMemcpy(AQs_ptr, AuxQs, nBytes3Dc, cudaMemcpyDeviceToDevice));
	CHECK(cudaMemcpy(AQi_ptr, AuxQi, nBytes3Dc, cudaMemcpyDeviceToDevice));
	
	cufftExecC2C(plan2D, (complex_t *)AQs_ptr, (complex_t *)As_ptr, CUFFT_INVERSE);
	CHECK(cudaDeviceSynchronize());
	cufftExecC2C(plan2D, (complex_t *)AQi_ptr, (complex_t *)Ai_ptr, CUFFT_INVERSE);
	CHECK(cudaDeviceSynchronize());

	cuFFT2DScaleSignalIdler<<<grid2D, block2D>>>(As_ptr, Ai_ptr);
	CHECK(cudaDeviceSynchronize());		
	
	CHECK(cudaFree(AuxQs)); CHECK(cudaFree(AuxQi));
	
	cufftDestroy(plan2D);
		
	return ;
}


template<typename Crystal>
void Solver<Crystal>::solverRK4( real_t z )
{	
	// A->Applies the Fourh-order Runge-Kutta Method with fixed step size dz
	// This function apply the fourth-order Runge-Kutta method	
	
	// Parameters for kernels
	dim3 block2D(BLKX, BLKY);	dim3 grid2D((NX+BLKX-1)/BLKX, (NY+BLKY-1)/BLKY);


	// Define pointers to use them in kernels
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
	
	complex_t * Ap_ptr  = thrust::raw_pointer_cast(A->Ap.data());
	complex_t * As_ptr  = thrust::raw_pointer_cast(A->As.data());
	complex_t * Ai_ptr  = thrust::raw_pointer_cast(A->Ai.data());

	complex_t * auxp_ptr = thrust::raw_pointer_cast(this->auxp.data());
	complex_t * auxs_ptr = thrust::raw_pointer_cast(this->auxs.data());
	complex_t * auxi_ptr = thrust::raw_pointer_cast(this->auxi.data());
	
	real_t dz = Cr->dz;
	real_t dk = A->dk; 
	real_t kp = A->kappa_p, ks = A->kappa_s, ki = A->kappa_i;

	//k1 = dAdz(kappas,dk,z,A)
	dAdz<<<grid2D,block2D>>>( k1p_ptr, k1s_ptr, k1i_ptr, Ap_ptr, As_ptr, Ai_ptr, kp, ks, ki, dk, z );
	CHECK(cudaDeviceSynchronize()); 

	//k2 = dAdz(kappas,dk,z+dz/2,A+k1/2) -> aux = A+k1/2
	kernelLinealCombination<<<grid2D,block2D>>>( auxp_ptr, auxs_ptr, auxi_ptr, Ap_ptr, As_ptr, Ai_ptr, k1p_ptr, k1s_ptr, k1i_ptr, 0.5f );
	CHECK(cudaDeviceSynchronize());   
	dAdz<<<grid2D,block2D>>>( k2p_ptr, k2s_ptr, k2i_ptr, auxp_ptr, auxs_ptr, auxi_ptr, kp, ks, ki, dk, z+dz/4.0f );
	CHECK(cudaDeviceSynchronize());

	// k3 = dAdz(kappas,dk,z+dz/2,A+k2/2)
	kernelLinealCombination<<<grid2D,block2D>>>( auxp_ptr, auxs_ptr, auxi_ptr, Ap_ptr, As_ptr, Ai_ptr, k2p_ptr, k2s_ptr, k2i_ptr, 0.5f );
	CHECK(cudaDeviceSynchronize());   
	dAdz<<<grid2D,block2D>>>( k3p_ptr, k3s_ptr, k3i_ptr, auxp_ptr, auxs_ptr, auxi_ptr, kp, ks, ki, dk, z+dz/4.0f );
	CHECK(cudaDeviceSynchronize());

	// k4 = dAdz(kappas,dk,z+dz,A+k3)
	kernelLinealCombination<<<grid2D,block2D>>>( auxp_ptr, auxs_ptr, auxi_ptr, Ap_ptr, As_ptr, Ai_ptr, k3p_ptr, k3s_ptr, k3i_ptr, 1.0f );
	CHECK(cudaDeviceSynchronize());   
	dAdz<<<grid2D,block2D>>>( k4p_ptr, k4s_ptr, k4i_ptr, auxp_ptr, auxs_ptr, auxi_ptr, kp, ks, ki, dk, z+dz/2.0f );
	CHECK(cudaDeviceSynchronize());

	// A = A + (k1 + 2*k2 + 2*k3 + k4)/6
	rk4<<<grid2D,block2D>>>(Ap_ptr, As_ptr, Ai_ptr, 
							k1p_ptr, k1s_ptr, k1i_ptr,
							k2p_ptr, k2s_ptr, k2i_ptr, 
							k3p_ptr, k3s_ptr, k3i_ptr,
							k4p_ptr, k4s_ptr, k4i_ptr, 
							dz/2 );
	CHECK(cudaDeviceSynchronize());

	return ;
	
}


template<typename Crystal>
void Solver<Crystal>::SSFM ( )
{
	for (real_t z = 0; z <= Cr->Lcr; z+=(Cr->dz))
	{	
		solverRK4(z); // RK4 in dz/2
		#ifdef DISPERSION
		dispersion(); // Dispersion in dz
		#endif
		diffractionInCrystal(); // Diffraction in dz
		solverRK4(z); // RK4 in dz/2
	}

	return ;
}


template<typename Crystal>
void Solver<Crystal>::runOPO(const std::vector<int>& save_roundtrips)
{
	
	std::unordered_set<int> save_set(save_roundtrips.begin(), save_roundtrips.end());
	
	#ifdef DISPERSION
	setDispersionPropagators();
	#endif

	if (checkCourantStability(Cr->dz , 0.2f*A->kp*Cr->dx*Cr->dy))
	{
		std::cout << "        ---> Courant convergence criteria: passed" << std::endl;
	
		real_t iStart = seconds();	// Timing code
		
		for( uint rt = 0; rt < NRT; rt++ ) // rt = round trip iteration
		{
			runOPO_status (rt, NRT/10); // print simulation status on screen
			
			A->Ap = A->Api; 
			// setDiffractionPropagators();			

			// Evolution in nonlinear crystal
			SSFM(); // compute couple-wave equations
			
			if(rt != (NRT-1)) // Evolution in free space for resonant fields
			{
				setDiffractionPropagatorsInFreeSpace(Cr->Lcr);
				diffractionInFreeSpace (); // From crystal output to concave mirror
			}

		 	// if(Cav->eom){EOM(A->As); EOM(A->Ai);}			// intracavity elements: EOM
		 	// if(Cav->gdd){GDD(A->As); GDD(A->Ai);}			// GDD
					
			// Apply boundary conditions to resonant fields	
			if(std::get<1>(A->resFields)){Cav->applyReflectivity(A->As, Cav->Rs, 0.0f);}
			if(std::get<2>(A->resFields)){Cav->applyReflectivity(A->Ai, Cav->Ri, 0.0f);}

			if ( rt%NRT/100 == 0 ){ // Check NaN subrutine
				cVech_t t_vec = A->As;	real_t test_nan =  t_vec[0].x;
				if(isnan(test_nan)){
					std::cout << "Simulation failed. Error in NaN check at r = " << rt << std::endl;
					rt = NRT;
				}
			}
			
			if (rt == 0){RemainingTime( rt, measuredTime(iStart) );}
		}

		real_t iElaps = seconds() - iStart;	// finish timing
		TimingCode( iElaps); // print time
	}
	else{
		std::cout << "        ---> Courant convergence criteria: failed." << std::endl;
		std::cout << "        ---> Please change step sizes. Courant factor = " << (0.2f*A->kp*Cr->dx*Cr->dy)/Cr->dz << " < 1" << std::endl;
	}
	

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



// /** This function compensates the GVD after a single-pass */
// __global__ void GDDKernel(complex_t *Aw, complex_t *aux, real_t *w, real_t GDD)
// {
	
// 	uint idt = threadIdx.x + blockDim.x*blockIdx.x;

// 	for (uint idy = 0; idy < NY; idy++){
// 		for (uint idx = 0; idx < NX; idx++){
// 			if( idt < NT ){
// 				aux[IDX(idx,idy,idt)] = CpxExp( 0.5*w[idt]*w[idt]*GDD ) * Aw[IDX(idx,idy,idt)];
// 			}
// 			if( idt < NT ){
// 				Aw[IDX(idx,idy,idt)] = aux[IDX(idx,idy,idt)];
// 			}
// 		}
// 	}

// 	return ;

// }

// /** Applies an electro optical modulator to an electric field after one round trip. */
// __global__ void EOMKernel(complex_t *A, complex_t *aux, real_t beta, real_t fpm, real_t *t)
// {
// 	uint idt = threadIdx.x + blockDim.x*blockIdx.x;
	
// 	for (uint idy = 0; idy < NY; idy++){
// 		for (uint idx = 0; idx < NX; idx++){
// 			if( idt < NT ){
// 				aux[IDX(idx,idy,idt)] = CpxExp(beta*sinf(2*PI*fpm*t[idt])) * A[IDX(idx,idy,idt)];
// 			}
// 			if( idt < NT ){
// 				A[IDX(idx,idy,idt)] = aux[IDX(idx,idy,idt)];
// 			}
// 		}
// 	}	
		
// 	return ;
// }


// template<typename Crystal>
// void Solver<Crystal>::EOM ( cVecd_t &Ax )
// {
// 	if(Cav->eom){
// 		dim3 block(BLKX); dim3 grid((SIZE+BLKX-1)/BLKX);

// 		complex_t *aux_ptr;
// 		CHECK( cudaMalloc((void **)&aux_ptr, Ax.size()*sizeof(complex_t)) );
// 		real_t* t_ptr = thrust::raw_pointer_cast( A->t.data() );
// 		complex_t* A_ptr = thrust::raw_pointer_cast( Ax.data() );
				
// 		EOMKernel<<<grid,block>>>( A_ptr, aux_ptr, Cav->beta, Cav->fpm, t_ptr );
// 		CHECK(cudaDeviceSynchronize());
// 		CHECK(cudaFree(aux_ptr));
// 	}
// 	return ;
// }


// template<typename Crystal>
// void Solver<Crystal>::GDD( cVecd_t &Ax )
// {
// 	if(Cav->gdd){
// 		// Parameters for kernels 1D
// 		dim3 block1D(BLKT);	dim3 grid1D((NT+BLKT-1)/BLKT);
// 		cufftHandle plan;	cufftPlan1d(&plan, SIZE, CUFFT_C2C, 1);

// 		complex_t *A_ptr = thrust::raw_pointer_cast(Ax.data());
// 		real_t* w_ptr = thrust::raw_pointer_cast(A->w.data());
		
// 		complex_t *Aw_ptr , *aux_ptr;
// 		CHECK(cudaMalloc((void **)&Aw_ptr, nBytes3Dc ));
// 		CHECK(cudaMalloc((void **)&aux_ptr, nBytes3Dc ));

// 		cufftExecC2C(plan, (complex_t *)A_ptr, (complex_t *)Aw_ptr, CUFFT_INVERSE);
// 		CHECK(cudaDeviceSynchronize());

// 		realScaling<<<grid1D,block1D>>>( Aw_ptr );
// 		CHECK(cudaDeviceSynchronize());
		
// 		real_t GDD = -(Cav->gamma)*(Cr->b2s)*(Cr->Lcr);
// 		GDDKernel<<<grid1D,block1D>>>(Aw_ptr, aux_ptr, w_ptr, GDD);
// 		CHECK(cudaDeviceSynchronize());
		
// 		cufftExecC2C(plan, (complex_t *)Aw_ptr, (complex_t *)A_ptr, CUFFT_FORWARD);
// 		CHECK(cudaDeviceSynchronize());
		
// 		cufftDestroy(plan);
// 		CHECK(cudaFree(Aw_ptr)); CHECK(cudaFree(aux_ptr));
// 	}

// 	return ;
// }
