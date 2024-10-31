/*---------------------------------------------------------------------------*/
// * This file contains the class Twm which models the three-wave mixing procces
// * in a nonlinear Cr
/*---------------------------------------------------------------------------*/


#ifndef _EFIELDSCUH
#define _EFIELDSCUH

#pragma once

///////////////////////////////////////////////////////////////////////////////////////////////////
// Kernel A² from A.
__global__ void kernelVectorPower2 ( real_t *A2, complex_t *A )
{	
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;
		
	for (uint idt = 0; idt < NT; idt++){	
		if( idx < NX and idy < NY ){
			A2[IDX(idx,idy,idt)] = CpxAbs2(A[IDX(idx,idy,idt)]) ;
		}
	}

	return ;
}


// Scales a vector after Fourier transforms in the time-frequency domaing
// (CUFFT_INVERSE mode)
__global__ void cuFFT1D_Scale(complex_t *A)
{		
	
	real_t size = static_cast<real_t>(NT);
	
	uint idt = threadIdx.x + blockDim.x*blockIdx.x;
			
	if( idt < SIZE ) {A[idt] = A[idt] / size;}
	
	return ;
	
}


// Copy a slice from a 3D tensor to a 2D matrix of complex numbers in GPU
__global__ void kernelGetSlice ( complex_t *Aux_Ax, complex_t *Ax_ptr, int slice )
{	
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;
		
	if( idx < NX and idy < NY ){
		Aux_Ax[IDX(idx,idy,0)] = Ax_ptr[IDX(idx,idy,slice)] ;
	}
	
	return ;
}


// Copy from a 2D matrix to a specific slice of a 3D tensor of complex numbers in GPU
__global__ void kernelSetSlice ( complex_t *Ax_ptr, complex_t *Aux_Ax, int slice )
{		
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;
		
	if( idx < NX and idy < NY ){
		Ax_ptr[IDX(idx,idy,slice)] = Aux_Ax[IDX(idx,idy,0)] ;
	}
	
	return ;
}


// Flips a vector for Fourier transforms
template<typename T>
void fftshift( T& V_flip, T V )
{
	int i, c = V.size()/2;
	for ( i = 0; i < V.size()/2; i++ ){
		V_flip[i+c] = V[i];
		V_flip[i]   = V[i+c];
	}
	
	return ;
}


// Set initial pump as a Gaussian beam in time and plane in XY cordinates
__global__ void setWavePlaneCW( complex_t *Ap_ptr, real_t Ap0, real_t waist, real_t uX, real_t uY, real_t dx, real_t dy )
{		
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){			
			Ap_ptr[IDX(idx,idy,idt)].x = (Ap0) * expf(((-powf((idx-uX)*dx,2)-powf((idy-uY)*dy,2))/(waist*waist)));
			Ap_ptr[IDX(idx,idy,idt)].y = 0.0f;
		}
	}
	
	return;
}


// Set initial pump as a Gaussian beam in time and plane in XY cordinates
__global__ void setWavePlanePulsed( complex_t *Ap_ptr, real_t *t_ptr, real_t Ap0, real_t waist, real_t tau, real_t uX, real_t uY, real_t dx, real_t dy )
{		
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){			
			Ap_ptr[IDX(idx,idy,idt)].x = (Ap0) * expf( -powf(t_ptr[idt]/tau, 2) - (powf((idx-uX)*dx,2)+powf((idy-uY)*dy,2))/(waist*waist) );
			Ap_ptr[IDX(idx,idy,idt)].y = 0.0f;
		}
		
	}
	
	return;
}


// Set initial pump as a CW beam in time and XY cordinates
__global__ void setFocusedCW( complex_t *Ap_ptr, real_t Ap0, complex_t MX, real_t waist, real_t uX, real_t uY, real_t dx, real_t dy )
{	
	
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){			
			Ap_ptr[IDX(idx,idy,idt)] = (Ap0/MX) * CpxExp(((-powf((idx-uX)*dx,2)-powf((idy-uY)*dy,2))/(waist*waist*MX)));
		}
	}
	
	return;
}


// Set initial pump as a Gaussian beam in time and XY cordinates
__global__ void setFocusedPulsed( complex_t *Ap_ptr, real_t *t_ptr, real_t Ap0, complex_t MX, real_t waist, real_t tau, real_t uX, real_t uY, real_t dx, real_t dy )
{	
	
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;
	
	for (uint idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){			
			Ap_ptr[IDX(idx,idy,idt)] = (Ap0/MX) * expf(-powf(t_ptr[idt]/tau, 2)) * CpxExp(((-powf((idx-uX)*dx,2)-powf((idy-uY)*dy,2))/(waist*waist*MX)));
		}
	}
	
	return;
}


// Swap horizontally the values un a matrix
__global__ void fftShift2DH( complex_t *Field, complex_t *aux)
{	 
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;

	uint c = (int) floor((real_t)NX/2);
	uint idt = 0;
	if (idx < c and idy < NY){
		Field[IDX(idx+c,idy,idt)]  =  aux[IDX(idx,idy,idt)];
		Field[IDX(idx,idy,idt)]    =  aux[IDX(idx+c,idy,idt)];
	}

	return ;
}


// Swap vertically the values un a matrix
__global__ void fftShift2DV( complex_t *Field, complex_t *aux)
{	// Swap vertically the values un a matrix
	uint idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint idy = threadIdx.y + blockDim.y*blockIdx.y;
	
	uint r = (int) floor((real_t)NY/2);
	uint idt = 0;
	if (idy < r and idx < NX){
		Field[IDX(idx,idy+r,idt)]  =  aux[IDX(idx,idy,idt)];
		Field[IDX(idx,idy,idt)]  =  aux[IDX(idx,idy+r,idt)];
	}  
	
	return ;
}


// Class Efields
template<typename Crystal>
class EFields
{
public:
	cVecd_t Api;
    cVecd_t Ap, As, Ai;
	cVecd_t Awp, Aws, Awi;
	cVecd_t AQp, AQs, AQi;
	cVecd_t eiQz_p, eiQz_s, eiQz_i;
	cVecd_t eiLz_p, eiLz_s, eiLz_i;
    rVecd_t t, F, w;

	real_t lp, ls, li, waist, Power;
	real_t np, ns, ni;
	real_t vp, vs, vi;
	real_t b2p, b2s, b2i;
	real_t alpha_crp, alpha_crs, alpha_cri;
	real_t kp, ks, ki;
	real_t kappa_p, kappa_s, kappa_i;
	real_t dx, dy, dz;
	real_t dk, dkp;	// mismatch and group-velocity mismatch
	std::tuple <bool, bool, bool> resFields;  // is resonant? <pump,signal,idler>

    // Constructor
    EFields(real_t _lp, real_t _ls, real_t _li, real_t _Power, real_t _waist, Crystal *Cr) :
			lp(_lp), ls(_ls), li(_li), Power(_Power), waist(_waist)
    {
        // Initialization or other constructor logic if needed
		this->Api.resize(SIZE); // initial pump efield
		this->Ap.resize(SIZE); this->As.resize(SIZE); this->Ai.resize(SIZE);
		this->Awp.resize(SIZE); this->Aws.resize(SIZE); this->Awi.resize(SIZE);
		this->AQp.resize(SIZE); this->AQs.resize(SIZE); this->AQi.resize(SIZE);
		this->eiQz_p.resize(NX*NY); this->eiQz_s.resize(NX*NY); this->eiQz_i.resize(NX*NY);
		this->eiLz_p.resize(NT); this->eiLz_s.resize(NT); this->eiLz_i.resize(NT);
		
		this->t.resize(NT); this->F.resize(NT); this->w.resize(NT);
		
		printLineOnScreen();
		printf("\nInstance of the class Efields.\n");
		np  = Cr->np; ns = Cr->ns; ni = Cr->ni;
		vp  = (Cr->GV(lp, Cr->T)); vs = (Cr->GV(ls, Cr->T)); vi = (Cr->GV(li, Cr->T));
		b2p  = (Cr->GVD(lp, Cr->T)); b2s = (Cr->GVD(ls, Cr->T)); b2i = (Cr->GVD(li, Cr->T));
		alpha_crp = Cr->alpha_crp; alpha_crs = Cr->alpha_crs; alpha_cri = Cr->alpha_cri;
		dx  = Cr->dx; dy = Cr->dy; dz = Cr->dz;
		kappa_p  = 2*PI*Cr->dQ/(np*lp);   // pump   kappa [1/V]
		kappa_s  = 2*PI*Cr->dQ/(ns*ls);   // signal kappa [1/V]
		kappa_i  = 2*PI*Cr->dQ/(ni*li);   // idler  kappa [1/V]
		kp  = 2.0f*PI*ns/lp;			  // pump   wavevector
		ks  = 2.0f*PI*ns/ls; 			  // signal wavevector
		ki  = 2.0f*PI*ni/ls;			  // idler  wavevector
		dk  = kp - ks - ki - 2.0f*PI/(Cr->Lambda); // mismatch factor
		dkp = 1/vp-1/vs;
    }

    // Destructor
	~EFields(){	printf("Efields Destructor.\n"); }

	// Methods definition
	void setPumpField( real_t Power, real_t FWHM, real_t waist, real_t focalpoint, std::string mode );
	void noiseGenerator ( cVecd_t& Vec );
	void setTimeFreqVectors(real_t trt);
	void fftShift2D ( cVecd_t& propagator );
	real_t averagePower(cVecd_t A);

};


// Methods declaration
template<typename Crystal>
void EFields<Crystal>::setPumpField( real_t Power, real_t FWHM, real_t waist, real_t focalpoint, std::string mode )
{
	complex_t Im; Im.x = 0; Im.y = 1;
	real_t tau    = FWHM*sqrtf(2)/(2*sqrtf(2*logf(2)));
	real_t Inten  = Power/(PI*waist*waist);
	real_t w02    = waist*waist;
	real_t zR     = PI*(this->np)*w02/lp;
	real_t eta    = focalpoint/zR;
	real_t Ap0    = sqrtf(4*Power/(EPS0*C*PI*(this->np)*w02));
	complex_t MX  = (1-Im*eta);

	complex_t *Ap_ptr = thrust::raw_pointer_cast(this->Api.data());
	real_t *t_ptr = thrust::raw_pointer_cast(this->t.data());
	dim3 block2D(BLKX, BLKY);	dim3 grid2D((NX+BLKX-1)/BLKX, (NY+BLKY-1)/BLKY);

	std::string mode1 = "waveplane-cw"; std::string mode2 = "waveplane-pulsed";
	std::string mode3 = "focused-cw"; std::string mode4 = "focused-pulsed";

	if (mode.compare(mode1) == 0){
		setWavePlaneCW<<<grid2D,block2D>>>( Ap_ptr, Ap0, waist, 0.5*NX, 0.5*NY, this->dx, this->dy );
		CHECK(cudaDeviceSynchronize()); 
		std::cout << "        ---> Pump field mode: " + mode << std::endl;
		std::cout << "        ---> Pump Power: " << Power << " W" << std::endl;
		std::cout << "        ---> Beam waist: " << waist << " \u03BCm\n" << std::endl;
	}
	else if (mode.compare(mode2) == 0){
		setWavePlanePulsed<<<grid2D,block2D>>>( Ap_ptr, t_ptr, Ap0, waist, tau, 0.5*NX, 0.5*NY, this->dx, this->dy );
		CHECK(cudaDeviceSynchronize()); 
		std::cout << "        ---> Pump field mode: " + mode << std::endl;
		std::cout << "        ---> Pump Power: " << Power << " W" << std::endl;
		std::cout << "        ---> Beam waist: " << waist << " \u03BCm\n" << std::endl;
	}
	else if (mode.compare(mode3) == 0){
		setFocusedCW<<<grid2D,block2D>>>( Ap_ptr, Ap0, MX, waist, 0.5*NX, 0.5*NY, this->dx, this->dy );
		CHECK(cudaDeviceSynchronize());
		std::cout << "        ---> Pump field mode: " + mode << std::endl;
		std::cout << "        ---> Setting Pump e-field with \u03BE = " << eta << std::endl;
		std::cout << "        ---> Pump Power: " << Power << " W" << std::endl;
		std::cout << "        ---> Beam waist: " << waist << " \u03BCm\n" << std::endl;
	}
	else if (mode.compare(mode4) == 0){
		setFocusedPulsed<<<grid2D,block2D>>>( Ap_ptr, t_ptr, Ap0, MX, waist, tau, 0.5*NX, 0.5*NY, this->dx, this->dy );
		CHECK(cudaDeviceSynchronize()); 
		std::cout << "        ---> Pump field mode: " + mode << std::endl;
		std::cout << "        ---> Setting Pump e-field with \u03BE = " << eta << std::endl;
		std::cout << "        ---> Pump Power: " << Power << " W" << std::endl;
		std::cout << "        ---> Beam waist: " << waist << " \u03BCm\n" << std::endl;
		
	}
	else{std::cout << "Invalid option for pump electric field!!\n";}

	return ;
}


template<typename Crystal>
void EFields<Crystal>::noiseGenerator ( cVecd_t& Vec )
{	// Noise generator for initial signal/idler vectors 
	uint seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::normal_distribution<real_t> distribution(0.0,1.0e-15);
	
	cVech_t A(Vec.size());
	real_t nsx, nsy;    
	for (int i=0; i<Vec.size(); ++i) {
		nsx = distribution(generator); A[i].x = static_cast<real_t>(nsx);
		nsy = distribution(generator); A[i].y = static_cast<real_t>(nsy);
	}
	thrust::copy(A.begin(), A.end(), Vec.begin());

	return ;	
}


template<typename Crystal>
void EFields<Crystal>::setTimeFreqVectors(real_t trt)
{
	std::cout << "\t---> Setting time and frequency vectors...";
	this->t = linspace<decltype(this->t)>( -trt*0.5, trt*0.5, this->t.size() );
	this->F = linspace<decltype(this->F)>( -0.5*this->F.size()/trt, +0.5*this->F.size()/trt, this->F.size() );
	fftshift<decltype(this->F)>(this->w, this->F) ;
	this->w *= (2.0*PI);
	std::cout << "done.\n" << std::flush;
	return ;
}


template<typename Crystal>
void EFields<Crystal>::fftShift2D ( cVecd_t& propagator )
{	// Standard fftshift in 2D

	complex_t *propagator_ptr = thrust::raw_pointer_cast(propagator.data());
	complex_t *aux;	CHECK(cudaMalloc((void **)&aux, nBytes2Dc));
	
	dim3 block2D(BLKX, BLKY);	dim3 grid2D((NX+BLKX-1)/BLKX, (NY+BLKY-1)/BLKY);
	
	CHECK(cudaMemcpy(aux, propagator_ptr, nBytes2Dc, cudaMemcpyDeviceToDevice));
	fftShift2DV<<<grid2D, block2D>>>(propagator_ptr, aux);
	cudaDeviceSynchronize();
	
	CHECK(cudaMemcpy(aux, propagator_ptr, nBytes2Dc, cudaMemcpyDeviceToDevice));
	fftShift2DH<<<grid2D, block2D>>>(propagator_ptr, aux);
	cudaDeviceSynchronize();
	
	CHECK(cudaFree(aux));
	
	return ;	
}


// Compute average power
template<typename Crystal>
real_t EFields<Crystal>::averagePower(cVecd_t A)
{
	// Parameters for kernels 2D
	dim3 block2D(BLKX, BLKY);	dim3 grid2D((NX+BLKX-1)/BLKX, (NY+BLKY-1)/BLKY);

	complex_t *A_ptr = thrust::raw_pointer_cast(A.data());
	real_t *A2; CHECK(cudaMalloc((void **)&A2, nBytes3Dr));
	thrust::device_ptr<real_t> A2_ptr(A2);

	kernelVectorPower2<<<grid2D, block2D>>> ( A2, A_ptr );
	CHECK(cudaDeviceSynchronize());

	real_t Power_avg = 0.5f*EPS0*C*(this->ns)*(this->dx)*(this->dy)*(thrust::reduce(A2_ptr, A2_ptr + SIZE));

	return Power_avg;
}


#endif // -> #ifdef _EFIELDSCUH