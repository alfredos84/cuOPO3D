/*---------------------------------------------------------------------------*/
// * This file contains a set of functions based on the 
// * Sellmeier equations for the PPLN nonlinear crystal and other 
// * properties of the χ⁽²⁾ material. Sellmeier equations from reference 
// * Dieter H. Jundt: Temperature-dependent Sellmeier equation for the index
// * of refraction, ne, in congruent lithium niobate
/*---------------------------------------------------------------------------*/

// All the functions have two input arguments:
// *     L: wavelenght in um
// *     T: temperature in degrees


#ifndef _PPLN 
#define _PPLN 

#pragma once


class PPLN
{
	public:
	real_t Lcr;			// crystal length [μm]
	real_t LX;          // Crystal width [μm]
	real_t LY;          // Crystal heigth [μm]
	real_t T; 			// crystal temperature [ºC]
	real_t Lambda;	   	// grating period for QPM [um]
	real_t lp, ls, li;  // wavelenghts

	// Define constants
	real_t d33;				// Eff. second-order susceptibility (d33) [um/V] [Ref]
	real_t dQ;				// Eff. second-order susceptibility for QPM [um/V]
	real_t alpha_crp; 		// pump linear absorption [1/μm]
	real_t alpha_crs;  		// signal linear absorption [1/μm]
	real_t alpha_cri;  		// idler linear absorption [1/μm]
	real_t beta_crs;		// signal 2-photons absorption [μm/W]
	real_t chi3p;			// χ⁽³⁾ in [um²/V²]
	real_t chi3s;			// χ⁽³⁾ in [um²/V²]
	real_t chi3i;			// χ⁽³⁾ in [um²/V²]  
	real_t np, ns, ni;  	// relevant refractive indexes
	real_t vp, vs, vi;  	// relevant group-velocities
	real_t b2p, b2s, b2i;  	// relevant group-velocities
	real_t b3p, b3s, b3i;  	// relevant group-velocities
	real_t dx;				// x step [μm]
    real_t dy;				// y step [μm]
	real_t dz; 			   	// step size
	real_t kp, ks, ki, dk; 	// mismatch factor
	std::string name;		// crystal name
	bool QPM;

	PPLN(real_t _LX, real_t _LY, real_t _Lcr, real_t _T, real_t _Lambda, real_t _lp, real_t _ls, real_t _li) :
			LX(_LX), LY(_LY), Lcr(_Lcr), T(_T), Lambda(_Lambda), lp(_lp), ls(_ls), li(_li)
	{	// Constructor
		printLineOnScreen();
		printf("\nInstance of the class PPLN.\n");

		name = "PPLN";
		QPM = true;
		d33 = 25.20e-6;	 dQ = 2.0*d33/PI;
		alpha_crp = 0.025e-4; alpha_crs = 0.002e-4; alpha_cri = 0.002e-4;
		#ifdef DIFFRACTION
		dx = static_cast<real_t> (LX/(NX-1.0f));
		dy = static_cast<real_t> (LY/(NY-1.0f));
		#endif
		dz = static_cast<real_t> (Lcr/(NZ-1.0f));
		this->LX = LX; this->LY = LY;
		this->Lcr = Lcr; this->T = T; this->Lambda = Lambda;
		this->lp = lp; this->ls = ls; this->li = li;
		
		np = this->n(lp, T); ns = this->n(ls, T); ni = this->n(li, T);
		vp = this->GV(lp, T); vs = this->GV(ls, T); vi = this->GV(li, T);
		b2p = this->GVD(lp, T); b2s = this->GVD(ls, T); b2i = this->GVD(li, T);
		b3p = this->TOD(lp, T); b3s = this->TOD(ls, T); b3i = this->TOD(li, T);b3p = this->TOD(lp, T); b3s = this->TOD(ls, T); b3i = this->TOD(li, T);
		kp  = 2.0f*PI*ns/lp;			  // pump   wavevector
		ks  = 2.0f*PI*ns/ls; 			  // signal wavevector
		ki  = 2.0f*PI*ni/ls;			  // idler  wavevector
		dk  = kp - ks - ki - 2.0f*PI/this->Lambda;
	}
	

	~PPLN(){printf("Destructor PPLN\n");}

	real_t n(real_t L,real_t T);
	real_t resonances(real_t L,real_t T, int p);
	real_t dndl(real_t L,real_t T);
	real_t d2ndl2(real_t L,real_t T);
	real_t d3ndl3(real_t L,real_t T);
	real_t GV(real_t L,real_t T);
	real_t GVD(real_t L,real_t T);
	real_t TOD(real_t L,real_t T);
	void getCrystalProp ();

};	

/** This function returns the PPLN extraordinary refractive index */
//__host__ __device__ 
real_t PPLN::n(real_t L,real_t T)
{
	
	real_t f = (T - 24.5) * (T + 570.82);
	real_t a1 = 5.35583;
	real_t a2 = 0.100473;
	real_t a3 = 0.20692;
	real_t a4 = 100.0;
	real_t a5 = 11.34927;
	real_t a6 =  1.5334e-2;
	real_t b1 =  4.629e-7;
	real_t b2 =  3.862e-8;
	real_t b3 =  -0.89e-8;
	real_t b4 =  2.657e-5;
	real_t G1 = a1 + b1*f;
	real_t G2 = a2 + b2*f;
	real_t G3 = a3 + b3*f;
	real_t G4 = a4 + b4*f;
	return sqrtf(G1+G2/(powf(L,2) - powf(G3,2))+G4/(powf(L,2) - powf(a5,2))-a6*L*L);
	
}


/** This function is an auxiliary function related with the resonances */
//__host__ __device__ 
real_t PPLN::resonances(real_t L,real_t T, int p)
{
	
	real_t f = (T - 24.5) * (T + 570.82);
	real_t a1 = 5.35583;
	real_t a2 = 0.100473;
	real_t a3 = 0.20692;
	real_t a4 = 100.0;
	real_t a5 = 11.34927;
	real_t a6 =  1.5334e-2;
	real_t b1 =  4.629e-7;
	real_t b2 =  3.862e-8;
	real_t b3 =  -0.89e-8;
	real_t b4 =  2.657e-5;
	real_t G1 = a1 + b1*f;
	real_t G2 = a2 + b2*f;
	real_t G3 = a3 + b3*f;
	real_t G4 = a4 + b4*f;

	return G2/powf((powf(L,2) - powf(G3,2)), p) + G4/powf((powf(L,2) - powf(a5,2)), p);
	
}


/** Returns the first-order derivative of the 
 * refractive index respect to the wavelength dn/dλ. */
//__host__ __device__ 
real_t PPLN::dndl(real_t L,real_t T)
{
	
	real_t f = (T - 24.5) * (T + 570.82);
	real_t a1 = 5.35583;
	real_t a2 = 0.100473;
	real_t a3 = 0.20692;
	real_t a4 = 100.0;
	real_t a5 = 11.34927;
	real_t a6 =  1.5334e-2;
	real_t b1 =  4.629e-7;
	real_t b2 =  3.862e-8;
	real_t b3 =  -0.89e-8;
	real_t b4 =  2.657e-5;
	real_t G1 = a1 + b1*f;
	real_t G2 = a2 + b2*f;
	real_t G3 = a3 + b3*f;
	real_t G4 = a4 + b4*f;

	return -L*(resonances(L, T, 2) + a6)/n(L, T);
	
}


/** Returns the second-order derivative of the
 * refractive index respect to the wavelength d²n/dλ². */
// __host__ __device__
real_t PPLN::d2ndl2(real_t L,real_t T)
{
	
	real_t f = (T - 24.5) * (T + 570.82);
	real_t a1 = 5.35583;
	real_t a2 = 0.100473;
	real_t a3 = 0.20692;
	real_t a4 = 100.0;
	real_t a5 = 11.34927;
	real_t a6 =  1.5334e-2;
	real_t b1 =  4.629e-7;
	real_t b2 =  3.862e-8;
	real_t b3 =  -0.89e-8;
	real_t b4 =  2.657e-5;
	real_t G1 = a1 + b1*f;
	real_t G2 = a2 + b2*f;
	real_t G3 = a3 + b3*f;
	real_t G4 = a4 + b4*f;


	real_t A  = (L*dndl(L,T)/powf(n(L,T),2)-1/n(L,T))*(resonances(L, T, 2) + a6);
	real_t B  = 4*L*L/n(L,T) * resonances(L, T, 3);
	
	return A+B;
	
}


/** Returns the third-order derivative of the
 * refractive index respect to the wavelength d³n/dλ³. */
// __host__ __device__ 
real_t PPLN::d3ndl3(real_t L,real_t T)
{
	
	real_t f = (T - 24.5) * (T + 570.82);
	real_t a1 = 5.35583;
	real_t a2 = 0.100473;
	real_t a3 = 0.20692;
	real_t a4 = 100.0;
	real_t a5 = 11.34927;
	real_t a6 =  1.5334e-2;
	real_t b1 =  4.629e-7;
	real_t b2 =  3.862e-8;
	real_t b3 =  -0.89e-8;
	real_t b4 =  2.657e-5;
	real_t G1 = a1 + b1*f;
	real_t G2 = a2 + b2*f;
	real_t G3 = a3 + b3*f;
	real_t G4 = a4 + b4*f;

	real_t A1 = (2*dndl(L,T)+L*d2ndl2(L,T))/powf(n(L,T),2);
	real_t A2 = -2*L*powf(dndl(L,T),2)/powf(n(L,T),3);
	real_t AA = (A1 + A2)*(resonances(L,T,2)+a6);
	real_t B1 = 12*L/n(L,T);
	real_t B2 = - 4*L*L/n(L,T)*d2ndl2(L,T)*(1-1/n(L,T));
	real_t BB = (B1+B2)*resonances(L,T,3);
	real_t CC = -24*L*L*L/n(L,T)*resonances(L,T,4);

	return AA + BB + CC;	
}


/** Returns the group-velocity vg(λ) = c/(n(λ)-λdn/dλ). */
// __host__ __device__ 
real_t PPLN::GV(real_t L,real_t T)
{		
	return C/(n(L,T)-L*dndl(L,T));
}


/** Returns the group-velocity β2(λ)=λ^3/(2πc²)(d²n/dλ²). */
// __host__ __device__ 
real_t PPLN::GVD(real_t L,real_t T)
{
	return powf(L,3)*d2ndl2(L, T)/(2*PI*C*C);
}


/** Returns the TOD β3(λ)=-λ^4/(4π²c³)[3.d²n/dλ² + λ.d³n/dλ³]. */
// __host__ __device__ 
real_t PPLN::TOD(real_t L,real_t T)
{
	return -powf(L,4)/(4*PI*PI*C*C*C)*(3*d2ndl2(L, T)+L*d3ndl3(L, T));
}


void PPLN::getCrystalProp () 
{
	std::cout << "Crystal name = " << this->name << std::endl;
	std::cout << "        ---> Temp            = " << T << " ºC" << std::endl;
	std::cout << "        ---> lp              = " << this->lp << " um" << std::endl;
	std::cout << "        ---> ls              = " << this->ls << " um" << std::endl;
	std::cout << "        ---> li              = " << this->li << " um" << std::endl;
	std::cout << "        ---> np              = " << this->np << std::endl;
	std::cout << "        ---> ns              = " << this->ns << std::endl;
	std::cout << "        ---> ni              = " << this->ni << std::endl;
	std::cout << "        ---> vgp             = " << this->vp << " um/ps" << std::endl;
	std::cout << "        ---> vgs             = " << this->vs << " um/ps" << std::endl;
	std::cout << "        ---> vgi             = " << this->vi << " um/ps" << std::endl;
	std::cout << "        ---> b2p             = " << this->b2p << " ps²/um" << std::endl;
	std::cout << "        ---> b2s             = " << this->b2s << " ps²/um" << std::endl;
	std::cout << "        ---> b2i             = " << this->b2i << " ps²/um" << std::endl;
	std::cout << "        ---> dQ              = " << dQ*1e6 << " pm/V"  << std::endl;
	std::cout << "        ---> \u039B               = " << Lambda << " \u03BCm"  << std::endl;
	std::cout << "        ---> \u03B1cp             = " << alpha_crp << " \u03BCm⁻¹"  << std::endl;
	std::cout << "        ---> \u03B1cs             = " << alpha_crs << " \u03BCm⁻¹" << std::endl;
	std::cout << "        ---> \u03B1ci             = " << alpha_cri << " \u03BCm⁻¹" << std::endl;
	std::cout << "        ---> dx              = " << dx << " \u03BCm"  << std::endl;
	std::cout << "        ---> dy              = " << dy << " \u03BCm"  << std::endl;
	std::cout << "        ---> dz              = " << dz << " \u03BCm"  << std::endl;
	std::cout << "        ---> Crystal length  = " << Lcr*1e-3 << " mm\n"  << std::endl;
	return ;
	
}

#endif // -> #ifdef _PPLN