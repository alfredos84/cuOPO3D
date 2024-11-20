/*---------------------------------------------------------------------------*/
// * This file contains the class Twm which models the three-wave mixing procces
// * in a nonlinear Cr
/*---------------------------------------------------------------------------*/


#ifndef _EFIELDSCUH
#define _EFIELDSCUH

#pragma once

///////////////////////////////////////////////////////////////////////////////////////////////////

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
__global__ void setWavePlaneCW( complex_t *Ap_ptr, real_t Ap0  )
{		
	uint idt = threadIdx.x + blockDim.x*blockIdx.x;

	if( idt < SIZE ){			
		Ap_ptr[idt].x = Ap0 ;
		Ap_ptr[idt].y = 0.0f;
	}

	return;
}


// Set initial pump as a Gaussian beam in time and plane in XY cordinates
__global__ void setWavePlanePulsed( complex_t *Ap_ptr, real_t *t_ptr, real_t Ap0, real_t tau )
{		
	uint idt = threadIdx.x + blockDim.x*blockIdx.x;
	
	if( idt < SIZE ){
		Ap_ptr[idt].x = (Ap0) * expf( -powf(t_ptr[idt]/tau, 2) ) ;
		Ap_ptr[idt].y = 0.0f;
	}

	return;
}


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



// Class Efields
template<typename Crystal>
class EFields
{
public:
	cVecd_t Api;
    cVecd_t Ap, As, Ai;
	cVecd_t Awp, Aws, Awi;
	cVecd_t eiLz_p, eiLz_s, eiLz_i;
    rVecd_t t, F, w;

	real_t lp, ls, li, waist, Power;
	real_t np, ns, ni;
	real_t vp, vs, vi;
	real_t b2p, b2s, b2i;
	real_t alpha_crp, alpha_crs, alpha_cri;
	real_t kp, ks, ki;
	real_t kappa_p, kappa_s, kappa_i;
	real_t dz;
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
		this->eiLz_p.resize(NT); this->eiLz_s.resize(NT); this->eiLz_i.resize(NT);
		
		this->t.resize(NT); this->F.resize(NT); this->w.resize(NT);
		
		printLineOnScreen();
		printf("\nInstance of the class Efields.\n");
		np  = Cr->np; ns = Cr->ns; ni = Cr->ni;
		vp  = (Cr->GV(lp, Cr->T)); vs = (Cr->GV(ls, Cr->T)); vi = (Cr->GV(li, Cr->T));
		b2p  = (Cr->GVD(lp, Cr->T)); b2s = (Cr->GVD(ls, Cr->T)); b2i = (Cr->GVD(li, Cr->T));
		alpha_crp = Cr->alpha_crp; alpha_crs = Cr->alpha_crs; alpha_cri = Cr->alpha_cri;
		dz = Cr->dz;
		kappa_p  = 2*PI*Cr->dQ/(np*lp);   // pump   kappa [1/V]
		kappa_s  = 2*PI*Cr->dQ/(ns*ls);   // signal kappa [1/V]
		kappa_i  = 2*PI*Cr->dQ/(ni*li);   // idler  kappa [1/V]
		kp  = 2.0f*PI*ns/lp;			  // pump   wavevector
		ks  = 2.0f*PI*ns/ls; 			  // signal wavevector
		ki  = 2.0f*PI*ns/ls;			  // idler  wavevector
		dk  = kp - ks - ki - 2.0f*PI/(Cr->Lambda); // mismatch factor
		dkp = 1/vp-1/vs;
    }

    // Destructor
	~EFields(){	printf("Efields Destructor.\n"); }

	// Methods definition
	void setPumpField( real_t Power, real_t FWHM, real_t waist, std::string mode );
	void noiseGenerator ( cVecd_t& Vec );
	void setTimeFreqVectors(real_t trt);
	void setDispersionPropagators();
	// void dispersion();

};


// Methods declaration
template<typename Crystal>
void EFields<Crystal>::setPumpField( real_t Power, real_t FWHM, real_t waist, std::string mode )
{
	real_t tau    = FWHM*sqrtf(2)/(2*sqrtf(2*logf(2)));
	real_t Inten  = Power/(PI*waist*waist);
	real_t w02    = waist*waist;
	real_t Ap0    = sqrtf(2*Inten/(EPS0*C*(this->np)));


	complex_t *Ap_ptr = thrust::raw_pointer_cast(this->Api.data());
	real_t *t_ptr = thrust::raw_pointer_cast(this->t.data());
	
	dim3 block1D(BLKT);	dim3 grid1D((NT+BLKT-1)/BLKT);

	std::string mode1 = "waveplane-cw"; std::string mode2 = "waveplane-pulsed";

	std::cout << "        ---> Setting Pump e-field" << std::endl;
	if (mode.compare(mode1) == 0){
		setWavePlaneCW<<<grid1D,block1D>>>( Ap_ptr, Ap0 );
		CHECK(cudaDeviceSynchronize()); 
		std::cout << "        ---> Pump field mode: " + mode << std::endl;
		std::cout << "        ---> Pump Power: " << Power << " W\n" << std::endl;
	}
	else if (mode.compare(mode2) == 0){
		setWavePlanePulsed<<<grid1D,block1D>>>( Ap_ptr, t_ptr, Ap0, tau );
		CHECK(cudaDeviceSynchronize()); 
		std::cout << "        ---> Pump field mode: " + mode << std::endl;
		std::cout << "        ---> Pump Power: " << Power << " W\n" << std::endl;
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


// Set vectors for dispersion propagators
template<typename Crystal>
void EFields<Crystal>::setDispersionPropagators()
{	
	// Parameters for kernels 2D
	dim3 block1D(BLKT);	dim3 grid1D((NT+BLKT-1)/BLKT);
	
	real_t *w_ptr = thrust::raw_pointer_cast(this->w.data());
	complex_t *eiLz_p_ptr = thrust::raw_pointer_cast(this->eiLz_p.data());
	complex_t *eiLz_s_ptr = thrust::raw_pointer_cast(this->eiLz_s.data());
	complex_t *eiLz_i_ptr = thrust::raw_pointer_cast(this->eiLz_i.data());

	dispersionPropagators<<<grid1D, block1D>>> (eiLz_p_ptr, eiLz_s_ptr, eiLz_i_ptr, 
				this->vp, this->vs, this->vi, this->b2p, this->b2s, this->b2i,
				this->alpha_crp, this->alpha_crs, this->alpha_cri,
				this->dz, w_ptr );
	CHECK(cudaDeviceSynchronize());
	
	return ;
}



#endif // -> #ifdef _EFIELDSCUH

