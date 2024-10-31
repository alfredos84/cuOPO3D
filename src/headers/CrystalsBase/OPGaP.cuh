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

// z discretization, time and frequency discretization
__constant__ const real_t Twin		= 10.0; 		// [ps]
// __constant__ const real_t Lcr		= 10e3;		// crystal length [um]
// __constant__ const real_t dz		= 33Lcr/(NZ-1);	// number of z-steps in the crystal
__constant__ const real_t dT		= Twin/SIZE;	// time step in [ps]
__constant__ const real_t dF		= 1/Twin;		// frequency step in [THz]


// Define global constants
__constant__ const real_t d14		= 70.6e-6;			// Eff. χ⁽²⁾ (d14) [um/V] [Ref2]
__constant__ const real_t dQ		= 2.0*d14/PI;		// Eff. χ⁽²⁾ for Orientation patterned Gap [um/V]
__constant__ const real_t k			= 114e-6;    		// thermal conductivity [W/um K]
__constant__ const real_t alpha_crp	= 0.025e-4; 		// pump linear absorption [1/μm]
__constant__ const real_t alpha_crs	= 0.002e-4;  		// signal linear absorption [1/μm]
__constant__ const real_t alpha_cri	= 0.002e-4;  		// idler linear absorption [1/μm]
__constant__ const real_t chi3p		= 3.0*1e-9;			// χ⁽³⁾ in [um²/V²]
__constant__ const real_t chi3s		= 1.5*1e-9;			// χ⁽³⁾ in [um²/V²]
__constant__ const real_t chi3i		= 1.0*1e-9;			// χ⁽³⁾ in [um²/V²]
// __constant__ const real_t Lambda	= 30.0;	   		// grating period for QPM [um]  

__constant__ const real_t beta_crs	= 0;			// signal 2-photons absorption [1/μm]
__constant__ const real_t rho		= 0;			// walk-off angle [rad] 


/** This function returns the MgO:PPLN extraordinary refractive index */
__host__ __device__ real_t n(real_t L,real_t T)
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
__host__ __device__ real_t resonances(real_t L,real_t T, int p)
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
__host__ __device__ real_t dndl(real_t L,real_t T)
{
	
	T += 273;
	return -L*(resonances(L, T, 2))/n(L, T);
	
}


/** Returns the second-order derivative of the
 * refractive index respect to the wavelength d²n/dλ². */
__host__ __device__ real_t d2ndl2(real_t L,real_t T)
{
	
	T += 273;
	real_t A  = (L*dndl(L,T)/powf(n(L,T),2)-1/n(L,T))*resonances(L, T, 2);
	real_t B  = 4*L*L/n(L,T) * resonances(L, T, 3);
	
	return A+B;
	
}


/** Returns the third-order derivative of the
 * refractive index respect to the wavelength d³n/dλ³. */
__host__ __device__ real_t d3ndl3(real_t L,real_t T)
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
__host__ __device__ real_t GV(real_t L,real_t T)
{
	
	return C/(n(L,T)-L*dndl(L,T));
}


/** Returns the group-velocity β2(λ)=λ^3/(2πc²)(d²n/dλ²). */
__host__ __device__ real_t GVD(real_t L,real_t T)
{
	return powf(L,3)*d2ndl2(L, T)/(2*PI*C*C);
}


/** Returns the TOD β3(λ)=-λ^4/(4π²c³)[3.d²n/dλ² + λ.d³n/dλ³]. */
__host__ __device__ real_t TOD(real_t L,real_t T)
{
	return -powf(L,4)/(4*PI*PI*C*C*C)*(3*d2ndl2(L, T)+L*d3ndl3(L, T));
}


void getCrystalProp ( real_t lp, real_t ls, real_t li, real_t Temp, real_t Lambda )
{
	std::cout << "\n\nUsing a OP_GaP nonlinear crystal\n\n " << std::endl;
	std::cout << "Using N                 = " << SIZE << " points" << std::endl;
	std::cout << "Pump wavelength         = " << lp*1e3 << " nm" << std::endl;
	std::cout << "Signal wavelength       = " << ls*1e3 << " nm" << std::endl;
	std::cout << "Idler wavelength        = " << li*1e3 << " nm" << std::endl;
	std::cout << "Temp                    = " << Temp << " ºC" << std::endl;
	std::cout << "np                      = " << n(lp, Temp) << std::endl;
	std::cout << "ns                      = " << n(ls, Temp) << std::endl;
	std::cout << "ni                      = " << n(li, Temp) << std::endl;
	std::cout << "\u03BD⁻¹ pump                = " << 1.0/GV(lp, Temp) << " ps/\u03BCm" << std::endl;
	std::cout << "\u03BD⁻¹ signal              = " << 1.0/GV(ls, Temp) << " ps/\u03BCm" << std::endl;
	std::cout << "\u03BD⁻¹ idler               = " << 1.0/GV(li, Temp) << " ps/\u03BCm" << std::endl;		
	#ifndef TWOEQS
	std::cout << "\u0394k                      = " << 2*PI*( abs(n(lp, Temp)/lp-n(ls, Temp)/ls-n(li, Temp)/li)-1/Lambda ) << " \u03BCm⁻¹" << std::endl;
	#else
	std::cout << "\u0394k                      = " << 2*PI*(abs(n(ls, Temp)/ls-2*n(lp, Temp)/lp) - 1/Lambda) << " \u03BCm⁻¹" << std::endl;
	#endif
	std::cout << "\u0394k'                     = " << 1/GV(lp, Temp)-1/GV(li, Temp) << " ps/\u03BCm" << std::endl;	
	std::cout << "GVD pump                = " << GVD(lp, Temp) << " ps²/\u03BCm" << std::endl;
	std::cout << "GVD signal              = " << GVD(ls, Temp) << " ps²/\u03BCm" << std::endl;
	std::cout << "GVD idler               = " << GVD(li, Temp) << " ps²/\u03BCm" << std::endl;		
	std::cout << "TOD pump                = " << TOD(lp, Temp) << " ps²/\u03BCm" << std::endl;
	std::cout << "TOD signal              = " << TOD(ls, Temp) << " ps²/\u03BCm" << std::endl;
	std::cout << "TOD idler               = " << TOD(li, Temp) << " ps²/\u03BCm" << std::endl;
	std::cout << "\u03A7p⁽³⁾                   = " << chi3p << " [\u03BCm²/V²]" << std::endl;
	std::cout << "dQ			= " << dQ*1e6 << " pm/V"  << std::endl;
	std::cout << "\u039B                       = " << Lambda << " \u03BCm"  << std::endl;
	std::cout << "\u03B1cp                     = " << alpha_crp << " \u03BCm⁻¹"  << std::endl;
	std::cout << "\u03B1cs                     = " << alpha_crs << " \u03BCm⁻¹" << std::endl;
	std::cout << "\u03B1ci                     = " << alpha_cri << " \u03BCm⁻¹" << std::endl;
	std::cout << "Crystal length, Lcr     = " << Lcr*1e-3 << " mm"  << std::endl;
	std::cout << "\u0394z                      = " << dz << " \u03BCm"  << std::endl;
	std::cout << "dT                      = " << dT << " ps" << std::endl;
	
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