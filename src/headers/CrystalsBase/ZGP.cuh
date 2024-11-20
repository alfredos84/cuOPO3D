/*---------------------------------------------------------------------------*/
// * This file contains a set of functions based on the 
// * Sellmeier equations for the ZGP nonlinear crystal and other 
// * properties of the χ⁽²⁾ material. Sellmeier equations from reference 
// * D. Zelmon et al: Refractive-index measurements and Sellmeier coefficients
// * for zinc germanium phosphide from 2 to 9 mm with implications for phase 
// * matching in optical frequency-conversion devices.
/*---------------------------------------------------------------------------*/

// All the functions have two input arguments:
// *     L: wavelenght in um


#ifndef _ZGP
#define _ZGP

class ZGP{
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
	real_t kp, ks, ki, dk; 	// mismatch factor
	std::string name;		// crystal name
	bool QPM;


	ZGP(real_t _LX, real_t _LY, real_t _Lcr, std::tuple<char,char,char> _pol, real_t _lp, real_t _ls, real_t _li) :
			LX(_LX), LY(_LY), Lcr(_Lcr), pol(_pol), lp(_lp), ls(_ls), li(_li)

	{	// Constructor	
		printLineOnScreen();
		printf("\nInstance of the class MgOsPPLT.\n");

		name = "ZPG";
		QPM	= false;
		d33 = 75.0e-6; dQ = 2.0*d33/PI;
		alpha_crp = 1.57e-4; alpha_crs = 0.17e-6; alpha_cri = 0.17e-6;
		#ifdef DIFFRACTION
		dx = static_cast<real_t> (LX/(NX-1.0f));
		dy = static_cast<real_t> (LY/(NY-1.0f));
		#endif
		dz = static_cast<real_t> (Lcr/(NZ-1.0f));
		this->LX = LX; this->LY = LY; this->Lcr = Lcr; 
		char pol_p, pol_s, pol_i; std::tie(pol_p, pol_s, pol_i) = pol;
		this->lp = lp; this->ls = ls; this->li = li;
		
		np = this->n(lp, pol_p); ns = this->n(ls, pol_s); ni = this->n(li, pol_i);
		vp = this->GV(lp, pol_p); vs = this->GV(ls, pol_s); vi = this->GV(li, pol_i);
		b2p = this->GVD(lp, pol_p); b2s = this->GVD(ls, pol_s); b2i = this->GVD(li, pol_i);
		kp  = 2.0f*PI*ns/lp;			  // pump   wavevector
		ks  = 2.0f*PI*ns/ls; 			  // signal wavevector
		ki  = 2.0f*PI*ni/ls;			  // idler  wavevector
		dk  = kp - ks - ki;
	}


	~ZGP(){printf("Destructor ZPG\n");}

	void setPolarizations(char p, char s, char i);
	real_t n( real_t L, char pol );
	real_t fb(real_t L, char pol, int p);
	real_t fd(real_t L, char pol, int p);
	real_t dndl(real_t L, char pol);
	real_t d2ndl2(real_t L, char pol);
	real_t GV(real_t L, char pol);
	real_t GVD(real_t L, char pol);
	void getCrystalProp ( );
	void getRefIndexRange ( real_t lmin, real_t lmax, char pol, uint nelem, std::string Filename );
	void getGVRange ( real_t lmin, real_t lmax, char pol, uint nelem, std::string Filename );
	
};


void ZGP::setPolarizations(char p, char s, char i)
{
	this->pol = std::make_tuple(p, s, i);
	return ;
}


/** This function returns the ZGP extraordinary refractive index */
//__host__ __device__ 
real_t ZGP::n( real_t L, char pol )
{
	real_t result;
	if ( pol == 'o' )
	{
		real_t A = 8.0409;
		real_t B = 1.68625;
		real_t CC = 0.40824;
		real_t D = 1.2880;
		real_t E = 611.05;
		result = sqrtf( A + (B*L*L/((L*L)-CC)) + (D*L*L/((L*L)-E)) );
	}
	if ( pol == 'e' )
	{
		real_t A = 8.0929;
		real_t B = 1.8649;
		real_t CC = 0.41468;
		real_t D = 0.84052;
		real_t E = 452.05;
		result = sqrtf( A + (B*L*L/((L*L)-CC)) + (D*L*L/((L*L)-E)) );
	}

	return result;
}


/** This function is an auxiliary function related with the resonances */
//__host__ __device__ 
real_t ZGP::fb(real_t L, char pol, int p)
{
	real_t result;

	if ( pol == 'o' )
	{
		real_t B = 1.68625;
		real_t CC = 0.40824;
		
		result = B*L*L/powf((L*L)-CC, p);
	}
	if ( pol == 'e' )
	{
		real_t B = 1.8649;
		real_t CC = 0.41468;
		
		result = B*L*L/powf((L*L)-CC, p);
	}
	
	return result;
}


//__host__ __device__ 
real_t ZGP::fd(real_t L, char pol, int p)
{
	real_t result;

	if ( pol == 'o' )
	{
		real_t D = 1.2880;
		real_t E = 611.05;
		
		result = D*L*L/powf((L*L)-E, p);
	}
	if ( pol == 'e' )
	{
		real_t D = 0.84052;
		real_t E = 452.05;
		
		result = D*L*L/powf((L*L)-E, p);
	}
	
	return result;
}


/** Returns the first-order derivative of the 
 * refractive index respect to the wavelength dn/dλ. */
//__host__ __device__ 
real_t ZGP::dndl(real_t L, char pol)
{
	real_t result;

	if ( pol == 'o' )
	{
		real_t CC = 0.40824;
		real_t E = 611.05;

		result = -( CC*fb(L, pol, 2) + E*fd(L, pol, 2) ) / ( n(L, pol)* L);
	}
	if ( pol == 'e' )
	{
		real_t CC = 0.41468;
		real_t E = 452.05;
	
		result = -( CC*fb(L, pol, 2) + E*fd(L, pol, 2) ) / ( n(L, pol)* L);
	}
	
	return result;
}


/** Returns the second-order derivative of the
 * refractive index respect to the wavelength d²n/dλ². */
//__host__ __device__ 
real_t ZGP::d2ndl2(real_t L, char pol)
{
	real_t result;

	if ( pol == 'o' )
	{
		real_t A = 8.0409;
		real_t B = 1.68625;
		real_t CC = 0.40824;
		real_t D = 1.2880;
		real_t E = 611.05;

		result = ( CC*fb(L, pol, 3)*(CC/(L*L)+3) + E*fd(L, pol, 3)*(E/(L*L)+3) - powf( dndl(L, pol), 2 ) ) / n(L, pol);
	}
	if ( pol == 'e' )
	{
		real_t A = 8.0929;
		real_t B = 1.8649;
		real_t CC = 0.41468;
		real_t D = 0.84052;
		real_t E = 452.05;
		result = ( CC*fb(L, pol, 3)*(CC/(L*L)+3) + E*fd(L, pol, 3)*(E/(L*L)+3) - powf( dndl(L, pol), 2 ) ) / n(L, pol);
	}
	
	return result;
}


/** Returns the group-velocity vg(λ) = c/(n(λ)-λdn/dλ). */
//__host__ __device__ 
real_t ZGP::GV(real_t L, char pol)
{
	
	return C/(n(L,pol)-L*dndl(L,pol));
}


/** Returns the group-velocity β2(λ)=λ^3/(2πc²)(d²n/dλ²). */
//__host__ __device__ 
real_t ZGP::GVD(real_t L, char pol)
{
	return powf(L,3)*d2ndl2(L, pol)/(2*PI*C*C);
}


void ZGP::getCrystalProp ()
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


#endif // -> #ifdef _ZGP