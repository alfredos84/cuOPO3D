/*---------------------------------------------------------------------------*/
// * This file contains four functions that save files in .dat extension
// * 1 - SaveFileVectorReal()       : save CPU real vectors 
// * 2 - SaveFileVectorRealGPU()    : save GPU real vectors
// * 3 - SaveFileVectorComplex()    : save CPU complex vectors
// * 4 - SaveFileVectorComplexGPU() : save GPU complex vectors

// Inputs:
// - Vector   : vector to save (stored on CPU or GPU)
// - N        : vector size
// - Filename : name of the saved file
/*---------------------------------------------------------------------------*/


#ifndef _FILESCUH
#define _FILESCUH

#pragma once


void SaveVectorReal (rVech_t V, std::string Filename)
{
	std::ofstream myfile;
	std::string extension = ".dat";
	myfile.open(Filename+extension);
	for (int i = 0; i < V.size(); i++)
		myfile << std::setprecision(20) << V[i] << "\n";
	myfile.close();
	
	return;
}


void SaveVectorRealGPU (rVecd_t V, std::string Filename)
{
	rVech_t Vh = V;
	SaveVectorReal ( Vh, Filename );
	
	return;
}


void SaveVectorComplex (cVech_t V, std::string Filename)
{
	std::ofstream myfile;
	std::string extension_r = "_r.dat", extension_i = "_i.dat";
	myfile.open(Filename+extension_r);
	for (int iy = 0; iy < V.size(); iy++)
		myfile << std::setprecision(20) << V[iy].x << "\n";
	myfile.close();
	myfile.open(Filename+extension_i);
	for (int iy = 0; iy < V.size(); iy++)
		myfile << std::setprecision(20) << V[iy].y << "\n";
	myfile.close();
	
	return;
	
}


void SaveVectorComplexGPU (cVecd_t V, std::string Filename)
{

	cVech_t Vh = V;
	SaveVectorComplex ( Vh, Filename );
	
	return;
}



void saveMatrixComplex_TimeSlice (cVech_t V, uint slice, std::string Filename)
{
	std::ofstream myfile;
	std::string filenamer = "_r.dat", filenamei = "_i.dat";
	myfile.open(Filename+filenamer);
	for (int iy = 0; iy < NY; iy++){
		for (int ix = 0; ix < NX; ix++)
			myfile << std::setprecision(20) << V[IDX(ix,iy,slice)].x << "\t";
		myfile << "\n"; 
	}
	myfile.close();
	myfile.open(Filename+filenamei);
	for (int iy = 0; iy < NY; iy++){
		for (int ix = 0; ix < NX; ix++)
			myfile << std::setprecision(20) << V[IDX(ix,iy,slice)].y << "\t";
		myfile << "\n"; 
	}
	myfile.close();
	
	return;
	
}


void saveMatrixComplex_TimeSlice_GPU (cVecd_t V, uint slice, std::string Filename)
{
	cVech_t Vh = V;
	saveMatrixComplex_TimeSlice ( Vh, slice, Filename );
	
	return;
	
}


void saveMatrixComplex (cVech_t Matrix, std::string Filename)
{
	std::cout << "Saving " + Filename << std::endl;
	std::ofstream myfile;
	std::string filenamer = "_r.dat", filenamei = "_i.dat";
	myfile.open(Filename+filenamer);
	for (int iy = 0; iy < NY; iy++){
		for (int ix = 0; ix < NX; ix++)
			myfile << std::setprecision(20) << Matrix[IDX(ix,iy,0)].x << "\t";
		myfile << "\n"; 
	}
	myfile.close();
	myfile.open(Filename+filenamei);
	for (int iy = 0; iy < NY; iy++){
		for (int ix = 0; ix < NX; ix++)
			myfile << std::setprecision(20) << Matrix[IDX(ix,iy,0)].y << "\t";
		myfile << "\n"; 
	}
	myfile.close();
	
	return;
	
}


void saveMatrixComplexGPU (cVecd_t Matrix_gpu, std::string Filename)
{
	// 	uint nBytes2D = NX * NY * sizeof(complex_t);
	cVech_t Matrix = Matrix_gpu;
	saveMatrixComplexGPU ( Matrix, Filename );
	
	return;
	
}


void SaveElectricField1D_time (cVech_t V, std::string Filename)
{
	std::ofstream myfile;
	std::string extension_r = "_r.dat", extension_i = "_i.dat";
	myfile.open(Filename+extension_r);
	for (int idt = 0; idt < NT; idt++)
		myfile << std::setprecision(20) << V[idt].x << "\n";
	myfile.close();
	myfile.open(Filename+extension_i);
	for (int idt = 0; idt < NT; idt++)
		myfile << std::setprecision(20) << V[idt].y << "\n";
	myfile.close();
	
	return;
	
}


void SaveElectricField1D_time_GPU  (cVecd_t V, std::string Filename)
{
	cVech_t Vh = V;
	SaveElectricField1D_time ( Vh, Filename);

	return;
}


void SaveElectricField3D_time (cVech_t V, uint x, uint y, std::string Filename)
{
	std::ofstream myfile;
	std::string extension_r = "_r.dat", extension_i = "_i.dat";
	myfile.open(Filename+extension_r);
	for (int idt = 0; idt < NT; idt++)
		myfile << std::setprecision(20) << V[IDX(x,y,idt)].x << "\n";
	myfile.close();
	myfile.open(Filename+extension_i);
	for (int idt = 0; idt < NT; idt++)
		myfile << std::setprecision(20) << V[IDX(x,y,idt)].y << "\n";
	myfile.close();
	
	return;
	
}


void SaveElectricField_time3D_GPU  (cVecd_t V, uint x, uint y, std::string Filename)
{
	cVech_t Vh = V;
	SaveElectricField3D_time ( Vh, x, y, Filename);

	return;
}


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