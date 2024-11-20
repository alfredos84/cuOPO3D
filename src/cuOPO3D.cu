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
	real_t lp = 0.532, ls = 2*lp, li = ls*lp/(ls-lp); // wavelengths

	real_t T = 27.0f, Lambda = 6.97f;	// temperature and grating period for QPM
	MgOsPPLT * cr1 = new MgOsPPLT(LX, LY, Lcr, T, Lambda, lp, ls, li);

	// std::tuple<char, char, char> pol = std::make_tuple('e','e','o'); // beam polarizations
	// ZGP * cr1 = new ZGP(LX, LY, Lcr, pol, lp, ls, li);

	cr1->getCrystalProp();


	// b. set cavity
	real_t Lcav = 4*(cr1->Lcr), Rs = 0.98f, Ri = Rs;
	Cavity<MgOsPPLT> * cav1 = new Cavity<MgOsPPLT>( cr1, Lcav, Rs, Ri );
	cav1->setGDDCompensation( 1.0f );
	cav1->getCavityProp ();

	// c. set electric fields	
	real_t alphas = 0.5*( (1-cav1->Rs)+cr1->alpha_crs*cr1->Lcr ), alphai = alphas; 
	real_t Power = atof(argv[1])*1e-3f, FWHM = 0.20f, focalpoint = Lcr/2.f, waist=atof(argv[2]);


	EFields<MgOsPPLT> * A = new EFields<MgOsPPLT>(lp, ls, li, Power, waist, cr1);
	A->resFields = std::make_tuple(false, true, true); //are (pump,signal,idler) resonant?
	A->setTimeFreqVectors( cav1->trt );

	// Set pump mode: (1) "waveplane-cw"; (2) "waveplane-pulsed"; (3) "focused-cw"; (4) "focused-pulsed"
	std::string mode = "focused-cw";
	#ifdef DIFFRACTION
	A->setPumpField( Power, FWHM, waist, focalpoint, mode ); A->Ap = A->Api;
	#else
	A->setPumpField( Power, FWHM, waist, mode ); A->Ap = A->Api;
	#endif
	// saveMatrixComplex_TimeSlice ( A->Ap, 0, "Input_pump_power" );
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

	// Save output signal electric field in x-y plane.
	// Temporal `slice` takes an instant in the time vector and save the beam profile.
	// By default is set `slice=NT/2` here. Users can uncomment the `for()` and replace 
	// `NT/2` by `slice` variable in case of need all the field evolution.	
	
	// for (uint slice = 0; slice < NT; slice ++)
	saveMatrixComplex_TimeSlice ( A->As, NT/2, "Output_signal_power_" + std::to_string(atoi(argv[1])) + "_mW" );


	// Save output signal electric field in the temporal domain.
	// This subroutine take each `x,y` pair and save the corresponding the time vector.
	for(uint x = 0; x < NX; x++){
		for(uint y = 0; y < NY; y++){
			SaveElectricField_time3D_GPU  (A->As, x, y, "As_temporal_x_" + std::to_string(x) + "_y_" + std::to_string(y));
		}
	}
	
	// saveAppendedValues2Columns("Average_power_Waist_" + std::to_string(atoi(argv[2])), Power, (solv1->averagePowerCutXY(A->As))*(1-cav1->Rs) );


	///////////////////////////////////////////////////////////////////////////////////
	// 4. Delete objects

	delete solv1; delete A, delete cav1; delete cr1;
	
	///////////////////////////////////////////////////////////////////////////////////
	
	return 0;
	
}