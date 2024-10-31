// * Author name: Alfredo Daniel Sanchez
// * email:       alfredo.daniel.sanchez@gmail.com

#include "headers/Libraries.cuh"		// Required libraries
#include "headers/PackageLibraries.cuh"	// Required package libraries

int main(int argc, char *argv[]){

	is_gpu_available();

	///////////////////////////////////////////////////////////////////////////////////
	// 1. Build the OPO

	// a. set crystal
	real_t LX = 1.e3f, LY = 1.e3f, Lcr = 29.e3f; // crystal dimensions
	printGrid();
	real_t T = 27.0f, Lambda = 6.97f;	// temperature and grating period for QPM
	real_t lp = 0.532, ls = 2*lp, li = ls*lp/(ls-lp); // wavelengths

	MgOsPPLT * cr1 = new MgOsPPLT(LX, LY, Lcr, T, Lambda, lp, ls, li);
	cr1->getCrystalProp();


	// b. set cavity
	real_t Lcav = 4*(cr1->Lcr), Rs = 0.98f, Ri = Rs;
	Cavity<MgOsPPLT> * cav1 = new Cavity<MgOsPPLT>( cr1, Lcav, Rs, Ri );
	cav1->getCavityProp ();
	// cav1->setEOM( 0.8*PI, cav1->fsr ); 
	// cav1->setGDD( 1.0 );


	// c. set electric fields	
	real_t alphas = 0.5*( (1-cav1->Rs)+cr1->alpha_crs*cr1->Lcr ), alphai = alphas; 
	real_t Ith   = EPS0*C*(cr1->np)*(cr1->ns)*(cr1->ni)*(cr1->ls)*(cr1->li)*alphas*alphai/(8*powf((cr1->dQ*cr1->Lcr*PI),2));

	real_t Power = atof(argv[1]), FWHM = 0.20f, focalpoint = Lcr/2.f, waist=atof(argv[2]);
	real_t Power_th = Ith*(waist*waist*PI); printVarOnScreen("Power threshold = ", Power_th);


	EFields<MgOsPPLT> * A = new EFields<MgOsPPLT>(lp, ls, li, Power, waist, cr1);
	A->resFields = std::make_tuple(false, true, true); //are (pump,signal,idler) resonant?
	A->setTimeFreqVectors( cav1->trt );

	// Set pump mode: (1) "waveplane-cw"; (2) "waveplane-pulsed"; (3) "focused-cw"; (4) "focused-pulsed"
	std::string mode = "focused-cw";
	#ifdef DIFFRACTION
	A->setPumpField( Power, FWHM, waist, focalpoint, mode );
	#else
	A->setPumpField( Power, FWHM, waist, mode ); A->Ap = A->Api;
	#endif

	A->noiseGenerator( A->As ); A->Ai = A->As;

	///////////////////////////////////////////////////////////////////////////////////
	// 2. Run OPO
	
	Solver<MgOsPPLT> * solv1 = new Solver<MgOsPPLT>(cr1, cav1, A);
	solv1->runOPO();

	// std::vector<int> save_roundtrips(128); // set of roundtrips to save
	// for (uint i = 0; i < save_roundtrips.size(); i++)
	// 	save_roundtrips[i] = NRT - save_roundtrips.size() + i;
	// solv1->runOPO(save_roundtrips);

	


	///////////////////////////////////////////////////////////////////////////////////
	// 3. Save output data

	// Save output signal electric field
	for (uint slice = 0; slice <= NT; slice += 16)
		saveMatrixComplex_TimeSlice ( A->As, NT/2, "Output_signal_Pump_power_" + std::to_string(atof(argv[1])) + "_W" );

	for(uint x = 0; x < NX; x++){
		for(uint y = 0; y < NY; y++){
			SaveElectricField_time3D_GPU  (A->As, x, y, "As_temporal_x_" + std::to_string(x) + "_y_" + std::to_string(y));
		}
	}
	

	// saveAppendedValues2Columns("Average_power_Waist_" + std::to_string(atoi(argv[2])), Power, (A->averagePower(A->As))*(1-cav1->Rs) );


	///////////////////////////////////////////////////////////////////////////////////
	// 4. Delete objects

	delete solv1; delete A, delete cav1; delete cr1;

	
	
	///////////////////////////////////////////////////////////////////////////////////
	
	return 0;
	
}