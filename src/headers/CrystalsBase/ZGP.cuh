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
	std::tuple<char, char, char> pol; // type process, e.g. e->eo is <'e','e','o'>
	real_t lp, ls, li;  // wavelenghts

	// Define constants
	real_t dQ;				// Eff. second-order susceptibility for QPM [um/V]
	real_t alpha_crp; 		// pump linear absorption [1/μm]
	real_t alpha_crs;  		// signal linear absorption [1/μm]
	real_t alpha_cri;  		// idler linear absorption [1/μm]
	real_t beta_crs;		// signal 2-photons absorption [μm/W]
	real_t chi3p;			// χ⁽³⁾ in [um²/V²]
	real_t chi3s;			// χ⁽³⁾ in [um²/V²]
	real_t chi3i;			// χ⁽³⁾ in [um²/V²]  
	char pol_p, pol_s, pol_i;  	// field polarization
	real_t np, ns, ni;  	// relevant refractive indexes
	real_t vp, vs, vi;  	// relevant group-velocities
	real_t b2p, b2s, b2i;  	// relevant group-velocities
	real_t b3p, b3s, b3i;  	// relevant group-velocities
	real_t dz; 			   	// step size
	std::string name;		// crystal name



	ZGP(real_t _Lcr, std::tuple<char,char,char> _pol, real_t _lp, real_t _ls, real_t _li) :
		Lcr(_Lcr), pol(_pol), lp(_lp), ls(_ls), li(_li)
	{	// Constructor
		printf("\nInstance of the class ZPG.\n\n");

		name = "ZPG";
		dQ = 75.00e-6;
		alpha_crp = 0.025e-4; alpha_crs = 0.002e-4; alpha_cri = 0.002e-4;
		dz = static_cast<real_t> (Lcr/NZ);
		this->Lcr = Lcr; tie(pol_p, pol_s, pol_i) = pol;
		this->lp = lp; this->ls = ls; this->li = li;

		np = this->n(lp, pol_p); ns = this->n(ls, pol_s); ni = this->n(li, pol_i);
		vp = this->GV(lp, pol_p); vs = this->GV(ls, pol_s); vi = this->GV(li, pol_i);
		b2p = this->GVD(lp, pol_p); b2s = this->GVD(ls, pol_s); b2i = this->GVD(li, pol_i);
		b3p = this->TOD(lp, pol_p); b3s = this->TOD(ls, pol_s); b3i = this->TOD(li, pol_i);
	}


	~ZPG(){printf("Destructor ZPG\n");}

	void setPolarizations(char p, char s, char i);
	real_t n( real_t L, chat pol );
	real_t fb(real_t L, chat pol, int p);
	real_t fd(real_t L, chat pol, int p);
	real_t dndl(real_t L, chat pol);
	real_t d2ndl2(real_t L, chat pol);
	real_t GV(real_t L, chat pol);
	real_t GVD(real_t L, chat pol);
	void getCrystalProp ( real_t lp, real_t ls, real_t li, chat pol, real_t Lambda );
	void getRefIndexRange ( real_t lmin, real_t lmax, chat pol, uint nelem, std::string Filename );
	void getGVRange ( real_t lmin, real_t lmax, chat pol, uint nelem, std::string Filename );
	
};


void ZGP::setPolarizations(char p, char s, char i)
{
	this->pol = make_tuple(p, s, i);
	return ;
}

/** This function returns the ZGP extraordinary refractive index */
__host__ __device__ 
real_t ZGP::n( real_t L, chat pol )
{
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
__host__ __device__ 
real_t ZGP::fb(real_t L, chat pol, int p)
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


__host__ __device__ 
real_t ZGP::fd(real_t L, chat pol, int p)
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
__host__ __device__ 
real_t ZGP::dndl(real_t L, chat pol)
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
__host__ __device__ 
real_t ZGP::d2ndl2(real_t L, chat pol)
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
__host__ __device__ 
real_t ZGP::GV(real_t L, chat pol)
{
	
	return C/(n(L,pol)-L*dndl(L,pol));
}


/** Returns the group-velocity β2(λ)=λ^3/(2πc²)(d²n/dλ²). */
__host__ __device__ 
real_t ZGP::GVD(real_t L, chat pol)
{
	return powf(L,3)*d2ndl2(L, pol)/(2*PI*C*C);
}


void ZGP::getCrystalProp ( real_t lp, real_t ls, real_t li, chat pol, real_t Lambda )
{
	std::cout << "\n\nUsing a OP_GaP nonlinear crystal\n\n " << std::endl;
	std::cout << "Using N                 = " << SIZE << " points" << std::endl;
	std::cout << "Pump wavelength         = " << lp*1e3 << " nm" << std::endl;
	std::cout << "Signal wavelength       = " << ls*1e3 << " nm" << std::endl;
	std::cout << "Idler wavelength        = " << li*1e3 << " nm" << std::endl;
	std::cout << "np                      = " << n(lp, 0) << std::endl;
	std::cout << "ns                      = " << n(ls, 1) << std::endl;
	std::cout << "ni                      = " << n(li, 0) << std::endl;
	std::cout << "\u03BD⁻¹ pump                = " << 1.0/GV(lp, 0) << " ps/\u03BCm" << std::endl;
	std::cout << "\u03BD⁻¹ signal              = " << 1.0/GV(ls, 1) << " ps/\u03BCm" << std::endl;
	std::cout << "\u03BD⁻¹ idler               = " << 1.0/GV(li, 0) << " ps/\u03BCm" << std::endl;		
	#ifndef TWOEQS
	std::cout << "\u0394k                      = " << 2*PI*( abs(n(lp, 0)/lp-n(ls, 1)/ls-n(li, 0)/li) ) << " \u03BCm⁻¹" << std::endl;
	#else
	std::cout << "\u0394k                      = " << 2*PI*(abs(n(ls, 1)/ls-2*n(lp, 0)/lp)) << " \u03BCm⁻¹" << std::endl;
	#endif
	std::cout << "\u0394k'                     = " << 1/GV(lp, 0)-1/GV(li, 0) << " ps/\u03BCm" << std::endl;	
	std::cout << "GVD pump                = " << GVD(lp, 0) << " ps²/\u03BCm" << std::endl;
	std::cout << "GVD signal              = " << GVD(ls, 0) << " ps²/\u03BCm" << std::endl;
	std::cout << "GVD idler               = " << GVD(li, 0) << " ps²/\u03BCm" << std::endl;
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


void ZGP::getRefIndexRange ( real_t lmin, real_t lmax, chat pol, uint nelem, std::string Filename )
{
	real_t *lambda = (real_t*)malloc(nBytesr);
	real_t *refind = (real_t*)malloc(nBytesr);
	linspace( lambda, nelem, lmin, lmax);
	std::ofstream myfile;	myfile.open(Filename);
	for (int i = 0; i < nelem; i++)
	{
		refind[i] = n(lambda[i], pol);
		myfile << std::setprecision(10) << lambda[i] << "\t" << refind[i] << "\n";
	}
	myfile.close();	free (lambda); free (refind);
	return ;
}


void ZGP::getGVRange ( real_t lmin, real_t lmax, chat pol, uint nelem, std::string Filename )
{
	real_t *lambda = (real_t*)malloc(nBytesr);
	real_t *group_v = (real_t*)malloc(nBytesr);
	linspace( lambda, nelem, lmin, lmax);
	std::ofstream myfile;	myfile.open(Filename);
	for (int i = 0; i < nelem; i++)
	{
		group_v[i] = GV(lambda[i], pol);
		myfile << std::setprecision(10) << lambda[i] << "\t" << group_v[i] << "\n";
	}
	myfile.close();	free (lambda); free (group_v);
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