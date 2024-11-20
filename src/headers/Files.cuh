/*---------------------------------------------------------------------------*/
// * This file contains four functions that save files in .dat extension

// Inputs:
// - Vector   : vector to save (stored on CPU, or GPU if sufix appears)
// - Filename : name of the saved file
/*---------------------------------------------------------------------------*/


#ifndef _FILESCUH
#define _FILESCUH

#pragma once


// Save real vector. 
void SaveVectorReal (rVech_t V, std::string Filename)
{
	std::ofstream myfile;
	std::string extension = ".dat";
	myfile.open(Filename+extension);
	for (uint32_t i = 0; i < V.size(); i++)
		myfile << std::setprecision(20) << V[i] << "\n";
	myfile.close();
	
	return;
}


// Save real vector on GPU.
void SaveVectorRealGPU (rVecd_t V, std::string Filename)
{
	rVech_t Vh = V;
	SaveVectorReal ( Vh, Filename );
	
	return;
}


// Save complex vector.
void SaveVectorComplex (cVech_t V, std::string Filename)
{
	std::ofstream myfile;
	std::string extension_r = "_r.dat", extension_i = "_i.dat";
	myfile.open(Filename+extension_r);
	for (uint32_t iy = 0; iy < V.size(); iy++)
		myfile << std::setprecision(20) << V[iy].x << "\n";
	myfile.close();
	myfile.open(Filename+extension_i);
	for (uint32_t iy = 0; iy < V.size(); iy++)
		myfile << std::setprecision(20) << V[iy].y << "\n";
	myfile.close();
	
	return;
	
}


// Save real vector on GPU.
void SaveVectorComplexGPU (cVecd_t V, std::string Filename)
{

	cVech_t Vh = V;
	SaveVectorComplex ( Vh, Filename );
	
	return;
}


// Save complex matrix ((x-y) beam profile)).
// For (3+1)D simulations, where time dependency is considered, `slice` is the 
// temporal slice to be saved.
void saveMatrixComplex_TimeSlice (cVech_t V, uint slice, std::string Filename)
{
	std::ofstream myfile;
	std::string filenamer = "_r.dat", filenamei = "_i.dat";
	myfile.open(Filename+filenamer);
	for (uint32_t iy = 0; iy < NY; iy++){
		for (uint32_t ix = 0; ix < NX; ix++)
			myfile << std::setprecision(20) << V[IDX(ix,iy,slice)].x << "\t";
		myfile << "\n"; 
	}
	myfile.close();
	myfile.open(Filename+filenamei);
	for (uint32_t iy = 0; iy < NY; iy++){
		for (uint32_t ix = 0; ix < NX; ix++)
			myfile << std::setprecision(20) << V[IDX(ix,iy,slice)].y << "\t";
		myfile << "\n"; 
	}
	myfile.close();
	
	return;
	
}


// Save complex matrix ((x-y) beam profile)) on GPU.
// For (3+1)D simulations, where time dependency is considered, `slice` is the 
// temporal slice to be saved.
void saveMatrixComplex_TimeSlice_GPU (cVecd_t V, uint slice, std::string Filename)
{
	cVech_t Vh = V;
	saveMatrixComplex_TimeSlice ( Vh, slice, Filename );
	
	return;
}


// Save complex matrix ((x-y) beam profile)).
void saveMatrixComplex (cVech_t Matrix, std::string Filename)
{
	std::cout << "Saving " + Filename << std::endl;
	std::ofstream myfile;
	std::string filenamer = "_r.dat", filenamei = "_i.dat";
	myfile.open(Filename+filenamer);
	for (uint32_t iy = 0; iy < NY; iy++){
		for (uint32_t ix = 0; ix < NX; ix++)
			myfile << std::setprecision(20) << Matrix[IDX(ix,iy,0)].x << "\t";
		myfile << "\n"; 
	}
	myfile.close();
	myfile.open(Filename+filenamei);
	for (uint32_t iy = 0; iy < NY; iy++){
		for (uint32_t ix = 0; ix < NX; ix++)
			myfile << std::setprecision(20) << Matrix[IDX(ix,iy,0)].y << "\t";
		myfile << "\n"; 
	}
	myfile.close();
	
	return;
	
}


// Save complex matrix ((x-y) beam profile)) on GPU.
void saveMatrixComplexGPU (cVecd_t Matrix_gpu, std::string Filename)
{
	cVech_t Matrix = Matrix_gpu;
	saveMatrixComplexGPU ( Matrix, Filename );
	
	return;
	
}


// Save complex electric field vector (similar to `SaveVectorReal()`).
void SaveElectricField1D_time (cVech_t V, std::string Filename)
{
	std::ofstream myfile;
	std::string extension_r = "_r.dat", extension_i = "_i.dat";
	myfile.open(Filename+extension_r);
	for (uint32_t idt = 0; idt < NT; idt++)
		myfile << std::setprecision(20) << V[idt].x << "\n";
	myfile.close();
	myfile.open(Filename+extension_i);
	for (uint32_t idt = 0; idt < NT; idt++)
		myfile << std::setprecision(20) << V[idt].y << "\n";
	myfile.close();
	
	return;
	
}


// Save complex electric field vector on GPU (similar to `SaveVectorRealGPU()`).
void SaveElectricField1D_time_GPU  (cVecd_t V, std::string Filename)
{
	cVech_t Vh = V;
	SaveElectricField1D_time ( Vh, Filename);

	return;
}


// Save complex electric field vector in (3+1)D simulations.
// `uint x` and `uint y` variables are the (x_i,y_j) cordinates of the full
// electric field. The saved vector is that A(x_i,y_j,z=output,t).
// For saving the full electric field in the time domain, use this function
// in a for-loop spanning all the posible x and y values (from 0 to NX-1,NY-1)
void SaveElectricField3D_time (cVech_t V, uint x, uint y, std::string Filename)
{
	std::ofstream myfile;
	std::string extension_r = "_r.dat", extension_i = "_i.dat";
	myfile.open(Filename+extension_r);
	for (uint32_t idt = 0; idt < NT; idt++)
		myfile << std::setprecision(20) << V[IDX(x,y,idt)].x << "\n";
	myfile.close();
	myfile.open(Filename+extension_i);
	for (uint32_t idt = 0; idt < NT; idt++)
		myfile << std::setprecision(20) << V[IDX(x,y,idt)].y << "\n";
	myfile.close();
	
	return;
	
}


// Save complex electric field vector in (3+1)D simulations in GPU.
// `uint x` and `uint y` variables are the (x_i,y_j) cordinates of the full
// electric field. The saved vector is that A(x_i,y_j,z=output,t).
// For saving the full electric field in the time domain, use this function
// in a for-loop spanning all the posible x and y values (from 0 to NX-1,NY-1)
void SaveElectricField_time3D_GPU  (cVecd_t V, uint x, uint y, std::string Filename)
{
	cVech_t Vh = V;
	SaveElectricField3D_time ( Vh, x, y, Filename);

	return;
}


// This function is useful to save and append 2 values in a file. 
// For instance, this function can be used, in combination with `Solver->averagePower()`
// function in massive simulations to compute the power threshold.
void saveAppendedValues2Columns(std::string Filename, real_t value_column1, real_t value_column2 )
{

	// Open a file in append mode
	std::ofstream outputFile(Filename + ".txt", std::ios::app);
	if (outputFile.is_open()) {
		// Write data to the end of the file
		outputFile << value_column1 << "\t" << value_column2 << std::endl;
		// Close the file
		outputFile.close();
		
		std::cout << "Data successfully written to the end of the file." << std::endl;
	} 
	else {
		std::cerr << "Unable to open the file for writing." << std::endl;
	}
	return;
}

#endif // -> #ifdef _FILESCUH