/*---------------------------------------------------------------------------*/
// * This file contains a set of functions based on the 
// * Sellmeier equations for the OP-GaP nonlinear crystal and other 
// * properties of the χ⁽²⁾ material. Sellmeier equations from reference 
// * J. Wei: Temperature dependent Sellmeier equation for the refractive
// * index of GaP.
/*---------------------------------------------------------------------------*/

// All the functions have two input arguments:
// *     L: wavelenght in um
// *     T: temperature in degrees


#ifndef _OPGAP
#define _OPGAP

#pragma once


class OPGaP{
public:
	
	real_t Lcr;			// crystal length [μm]
	real_t LX;          // Crystal width [μm]
	real_t LY;          // Crystal heigth [μm]
	real_t lp, ls, li;  // wavelenghts
	std::tuple<char, char, char> pol; // type process, e.g. e->eo is <'e','e','o'>
	
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
	std::string name;		// crystal name



	OPGaP(real_t _LX, real_t _LY, real_t _Lcr, std::tuple<char,char,char> _pol, real_t _lp, real_t _ls, real_t _li) :
			LX(_LX), LY(_LY), Lcr(_Lcr), pol(_pol), lp(_lp), ls(_ls), li(_li)
	{	// Constructor	
		printLineOnScreen();
		printf("\nInstance of the class OPGaP.\n");

		name = "OPGaP";
		dQ = 70.6e-6;
		alpha_crp = 0.025e-4; alpha_crs = 0.002e-4; alpha_cri = 0.002e-4;
		#ifdef DIFFRACTION
		dx = static_cast<real_t> (LX/(NX-1.0f));
		dy = static_cast<real_t> (LY/(NY-1.0f));
		#endif
		dz = static_cast<real_t> (Lcr/NZ);
		this->LX = LX; this->LY = LY;
		char pol_p, pol_s, pol_i;
		this->Lcr = Lcr; tie(pol_p, pol_s, pol_i) = pol;
		this->lp = lp; this->ls = ls; this->li = li;
		
		np = this->n(lp, pol_p); ns = this->n(ls, pol_s); ni = this->n(li, pol_i);
		vp = this->GV(lp, pol_p); vs = this->GV(ls, pol_s); vi = this->GV(li, pol_i);
		b2p = this->GVD(lp, pol_p); b2s = this->GVD(ls, pol_s); b2i = this->GVD(li, pol_i);
		b3p = this->TOD(lp, pol_p); b3s = this->TOD(ls, pol_s); b3i = this->TOD(li, pol_i);
	}


	~OPGaP(){printf("Destructor OPGaP\n");}

	void setPolarizations(char p, char s, char i);
	real_t n( real_t L, char pol );
	real_t resonances(real_t L, char pol, int p);
	real_t dndl(real_t L, char pol);
	real_t d2ndl2(real_t L, char pol);
	real_t d3ndl3(real_t L, char pol);
	real_t GV(real_t L, char pol);
	real_t GVD(real_t L, char pol);
	real_t TOD(real_t L, char pol);
	void getCrystalProp ();
	
};



/** This function returns the MgO:PPLN extraordinary refractive index */
// __host__ __device__ 
real_t OPGaP::n(real_t L,char pol)
{
	
	T += 273;
	real_t A = 10.926 + (7.0787e-4*T)+ (1.8594e-7*T*T);
	real_t B = 0.53718 + (5.8035e-5*T)+ (1.9819e-7*T*T);
	real_t C = 0.0911014;
	real_t D = 1504 + (0.25935*T) - (0.00023326*T*T);
	real_t E = 758.048;
	
	return sqrtf( A + (B/((L*L)-C)) + (D/((L*L)-E)) );

}


/** This function is an auxiliary function related with the resonances */
// __host__ __device__ 
real_t OPGaP::resonances(real_t L,char pol, int p)
{
	
	T += 273;
	real_t B = 0.53718 + (5.8035e-5*T)+ (1.9819e-7*T*T);
	real_t C = 0.0911014;
	real_t D = 1504+(0.25935*T) - (0.00023326*T*T);
	real_t E = 758.048;

	return powf(B/((L*L)-C), p) + powf(D/((L*L)-E), p);
	
}


/** Returns the first-order derivative of the 
 * refractive index respect to the wavelength dn/dλ. */
// __host__ __device__ 
real_t OPGaP::dndl(real_t L,char pol)
{
	
	T += 273;
	return -L*(resonances(L, T, 2))/n(L, T);
	
}


/** Returns the second-order derivative of the
 * refractive index respect to the wavelength d²n/dλ². */
// __host__ __device__ 
real_t OPGaP::d2ndl2(real_t L,char pol)
{
	
	T += 273;
	real_t A  = (L*dndl(L,T)/powf(n(L,T),2)-1/n(L,T))*resonances(L, T, 2);
	real_t B  = 4*L*L/n(L,T) * resonances(L, T, 3);
	
	return A+B;
	
}


/** Returns the third-order derivative of the
 * refractive index respect to the wavelength d³n/dλ³. */
// __host__ __device__ 
real_t OPGaP::d3ndl3(real_t L,char pol)
{
	T += 273;
	real_t A1 = (2*dndl(L,T)+L*d2ndl2(L,T))/powf(n(L,T),2);
	real_t A2 = -2*L*powf(dndl(L,T),2)/powf(n(L,T),3);
	real_t AA = (A1 + A2)*(resonances(L,T,2));
	real_t B1 = 12*L/n(L,T);
	real_t B2 = - 2*powf(2*L/n(L,T),2)*dndl(L,T);
	real_t BB = (B1+B2)*resonances(L,T,3);
	real_t CC = -24*L*L*L/n(L,T)*resonances(L,T,4);

	return AA + BB + CC;	
}


/** Returns the group-velocity vg(λ) = c/(n(λ)-λdn/dλ). */
// __host__ __device__ 
real_t OPGaP::GV(real_t L,char pol)
{
	
	return C/(n(L,T)-L*dndl(L,T));
}


/** Returns the group-velocity β2(λ)=λ^3/(2πc²)(d²n/dλ²). */
// __host__ __device__ 
real_t OPGaP::GVD(real_t L,char pol)
{
	return powf(L,3)*d2ndl2(L, T)/(2*PI*C*C);
}


/** Returns the TOD β3(λ)=-λ^4/(4π²c³)[3.d²n/dλ² + λ.d³n/dλ³]. */
// __host__ __device__ 
real_t OPGaP::TOD(real_t L,char pol)
{
	return -powf(L,4)/(4*PI*PI*C*C*C)*(3*d2ndl2(L, T)+L*d3ndl3(L, T));
}


void OPGaP::getCrystalProp ()
{
	std::cout << "Crystal name = " << this->name << std::endl;
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
	std::cout << "        ---> \u03B1cp             = " << alpha_crp << " \u03BCm⁻¹"  << std::endl;
	std::cout << "        ---> \u03B1cs             = " << alpha_crs << " \u03BCm⁻¹" << std::endl;
	std::cout << "        ---> \u03B1ci             = " << alpha_cri << " \u03BCm⁻¹" << std::endl;
	std::cout << "        ---> dx              = " << dx << " \u03BCm"  << std::endl;
	std::cout << "        ---> dy              = " << dy << " \u03BCm"  << std::endl;
	std::cout << "        ---> dz              = " << dz << " \u03BCm"  << std::endl;
	std::cout << "        ---> Crystal length  = " << Lcr*1e-3 << " mm\n"  << std::endl;
	return ;	
}


/**
 * Letter   Description  Escape-Sequence
 * -------------------------------------
 * A        Alpha        \u0391
 * B        Beta         \u0392
 * Γ        Gamma        \u0393
 * Δ        Delta        \u0394
 * Ε        Epsilon      \u0395
 * Ζ        Zeta         \u0396
 * Η        Eta          \u0397
 * Θ        Theta        \u0398
 * Ι        Iota         \u0399
 * Κ        Kappa        \u039A
 * TOD idlerΛ        Lambda       \u039B
 * Μ        Mu           \u039C
 * Ν        Nu           \u039D
 * Ξ        Xi           \u039E
 * Ο        Omicron      \u039F
 * Π        Pi           \u03A0
 * Ρ        Rho          \u03A1
 * Σ        Sigma        \u03A3
 * Τ        Tau          \u03A4
 * Υ        Upsilon      \u03A5
 * Φ        Phi          \u03A6
 * Χ        Chi          \u03A7
 * Ψ        Psi          \u03A8
 * Ω        Omega        \u03A9 
 * -------------------------------------
 * Letter   Description  Escape-Sequence
 * -------------------------------------
 * α        Alpha        \u03B1
 * β        Beta         \u03B2
 * γ        Gamma        \u03B3
 * δ        Delta        \u03B4
 * ε        Epsilon      \u03B5
 * ζ        Zeta         \u03B6
 * η        Eta          \u03B7
 * θ        Theta        \u03B8
 * ι        Iota         \u03B9
 * κ        Kappa        \u03BA
 * λ        Lambda       \u03BB
 * μ        Mu           \u03BC
 * ν        Nu           \u03BD
 * ξ        Xi           \u03BE
 * ο        Omicron      \u03BF
 * π        Pi           \u03C0
 * ρ        Rho          \u03C1
 * σ        Sigma        \u03C3
 * τ        Tau          \u03C4
 * υ        Upsilon      \u03C5
 * φ        Phi          \u03C6
 * χ        Chi          \u03C7
 * ψ        Psi          \u03C8
 * ω        Omega        \u03C9
 * -------------------------------------
 */


#endif // -> #ifdef _OPGAP