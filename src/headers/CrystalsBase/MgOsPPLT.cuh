/*---------------------------------------------------------------------------*/
// * This file contains a set of functions based on the 
// * Sellmeier equations for the MgO:sPPLT nonlinear crystal and other 
// * properties of the χ⁽²⁾ material. Sellmeier equations from reference 
// * Bruner et. al.: Temperature-dependent Sellmeier equation for the 
// * refractive index of stoichiometric lithium tantalate.
/*---------------------------------------------------------------------------*/

// All the functions have two input arguments:
// *     L: wavelenght in um
// *     T: temperature in degrees


#ifndef _MGOSPPLT 
#define _MGOSPPLT 

#pragma once


class MgOsPPLT
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

	MgOsPPLT(real_t _LX, real_t _LY, real_t _Lcr, real_t _T, real_t _Lambda, real_t _lp, real_t _ls, real_t _li) :
			LX(_LX), LY(_LY), Lcr(_Lcr), T(_T), Lambda(_Lambda), lp(_lp), ls(_ls), li(_li)
	{	// Constructor
		printLineOnScreen();
		printf("\nInstance of the class MgOsPPLT.\n");

		name = "MgOsPPLT";
		QPM = true;
		d33 = 11.0e-6; dQ = 2.0*d33/PI;
		alpha_crp = 1.57e-4; alpha_crs = 0.17e-6; alpha_cri = 0.17e-6;
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
		b3p = this->TOD(lp, T); b3s = this->TOD(ls, T); b3i = this->TOD(li, T);
		kp  = 2.0f*PI*ns/lp;			  // pump   wavevector
		ks  = 2.0f*PI*ns/ls; 			  // signal wavevector
		ki  = 2.0f*PI*ni/ls;			  // idler  wavevector
		dk  = kp - ks - ki - 2.0f*PI/this->Lambda;
	}
	

	~MgOsPPLT(){printf("MgOsPPLT Destructor\n");}
	
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

//__host__ __device__
real_t MgOsPPLT::n(real_t L,real_t T)
{

real_t A =  4.502483;
real_t B =  0.007294;
real_t C =  0.185087;
real_t D =  -0.02357;
real_t E =  0.073423;
real_t F =  0.199595;
real_t G =  0.001;
real_t H =  7.99724;
real_t b =  3.483933e-8 * pow(T + 273.15,2);
real_t c =  1.607839e-8 * pow(T + 273.15,2);

return sqrt( A + (B+b)/(pow(L,2)-pow((C+c),2)) + E/(pow(L,2)-pow(F,2)) + G/(pow(L,2)-pow(H,2))+ D*pow(L,2));

}


/** This function is an auxiliary function related with the resonances */
//__host__ __device__
real_t MgOsPPLT::resonances(real_t L,real_t T, int p)
{

real_t A =  4.502483;
real_t B =  0.007294;
real_t C =  0.185087;
real_t D =  -0.02357;
real_t E =  0.073423;
real_t F =  0.199595;
real_t G =  0.001;
real_t H =  7.99724;
real_t b =  3.483933e-8 * pow(T + 273.15,2);
real_t c =  1.607839e-8 * pow(T + 273.15,2);

return (B+b)/powf((powf(L,2) - powf((C+c),2)), p) + E/powf((powf(L,2) - powf(F,2)), p) + G/powf((powf(L,2) - powf(H,2)), p);
}


/** Returns the first-order derivative of the
 * refractive index respect to the wavelength dn/dλ. */
//__host__ __device__
real_t MgOsPPLT::dndl(real_t L,real_t T)
{
	
	real_t B =  0.007294;
	real_t C =  0.185087;
	real_t D =  -0.02357;
	real_t E =  0.073423;
	real_t F =  0.199595;
	real_t G =  0.001;
	real_t H =  7.99724;
	real_t b =  3.483933e-8 * pow(T + 273.15,2);
	real_t c =  1.607839e-8 * pow(T + 273.15,2);
	
	return -L/n(L,T)*( resonances(L,T,2) - D );
}


/** Returns the second-order derivative of the
 * refractive index respect to the wavelength d²n/dλ². */
//__host__ __device__
real_t MgOsPPLT::d2ndl2(real_t L,real_t T)
{
	
	real_t A =  4.502483;
	real_t B =  0.007294;
	real_t C =  0.185087;
	real_t D =  -0.02357;
	real_t E =  0.073423;
	real_t F =  0.199595;
	real_t G =  0.001;
	real_t H =  7.99724;
	real_t b =  3.483933e-8 * pow(T + 273.15,2);
	real_t c =  1.607839e-8 * pow(T + 273.15,2);


	return (L/powf(n(L,T),2)*dndl(L,T)-1/n(L,T))*(resonances(L,T,2)-D) + 4*L*L/n(L,T)*resonances(L,T,3);	
}


/** Returns the third-order derivative of the
 * refractive index respect to the wavelength d³n/dλ³. */
//__host__ __device__
real_t MgOsPPLT::d3ndl3(real_t L,real_t T)
{
	
	real_t A =  4.502483;
	real_t B =  0.007294;
	real_t C =  0.185087;
	real_t D =  -0.02357;
	real_t E =  0.073423;
	real_t F =  0.199595;
	real_t G =  0.001;
	real_t H =  7.99724;
	real_t b =  3.483933e-8 * powf(T + 273.15,2);
	real_t c =  1.607839e-8 * powf(T + 273.15,2);

	real_t A1 = (2*dndl(L,T)+L*d2ndl2(L,T))/powf(n(L,T),2);
	real_t A2 = -2*L*powf(dndl(L,T),2)/powf(n(L,T),3);
	real_t AA = (A1 + A2)*(resonances(L,T,2)-D);
	real_t B1 = 12*L/n(L,T);
	real_t B2 = -4*L*L/n(L,T)*d2ndl2(L,T)*(1-1/n(L,T));
	real_t BB = (B1+B2)*resonances(L,T,3);
	real_t CC = -24*L*L*L/n(L,T)*resonances(L,T,4);

	return AA + BB + CC;	
}

/** Returns the group-velocity vg(λ) = c/(n(λ)-λdn/dλ). */
//__host__ __device__
real_t MgOsPPLT::GV(real_t L,real_t T)
{
	return C/(n(L,T)-L*dndl(L,T));
}


/** Returns the group-velocity β(λ)=λ^3/(2πc²)(d²n/dλ²). */
//__host__ __device__
real_t MgOsPPLT::GVD(real_t L,real_t T)
{
	return pow(L,3)*d2ndl2(L, T)/(2*PI*C*C);
}


/** Returns the TOD β3(λ)=-λ^4/(4π²c³)[3.d²n/dλ² + λ.d³n/dλ³]. */
//__host__ __device__
real_t MgOsPPLT::TOD(real_t L,real_t T)
{
	return -powf(L,4)/(4*PI*PI*C*C*C)*(3*d2ndl2(L, T)+L*d3ndl3(L, T));
}


void MgOsPPLT::getCrystalProp () 
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


#endif // -> #ifdef _MGOSPPLT