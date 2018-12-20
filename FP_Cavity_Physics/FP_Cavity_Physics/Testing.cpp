#ifndef ATTACH_H
#include "Attach.h"
#endif

void testing::matrix_file_IO()
{
	// check that file IO for reading into 2D matrix is working correctly
	// R. Sheehan 18 - 12 - 2018

	std::string filename;

	filename = "RI_Data\\Si.txt";
	//filename = "Silicon_Refractive_Index_Data.csv";

	int nc, nr;
	std::vector<std::vector<double>> data;

	vecut::read_into_matrix(filename, data, nr, nc, true);

	std::vector<double> xa, ya; 

	xa = vecut::get_col(data, 0); ya = vecut::get_col(data, 1); 

	double x, y, dy; 

	x = 1567; 

	interpolation::polint(xa, ya, x, y, dy); 

	std::cout << "x: " << x << ", y: " << y << ", dy: " << dy << "\n"; 

	data.clear();
}

void testing::material_values() 
{
	// How do you access the material classes in a generic way without having to have a switch statement of some kind? 
	// R. Sheehan 20 - 12 - 2018
	
	material *mat1; 
	material *mat2; 

	Si smpl1; 
	SiO2 smpl2; 

	mat1 = &smpl1; 
	mat2 = &smpl2; 

	double wavelength = 1.55; 

	mat1->set_wavelength(wavelength);
	mat2->set_wavelength(wavelength); 

	std::cout<<mat1->refractive_index()<<"\n"; // this should return the RI of Si

	smpl1.set_wavelength(wavelength);
	std::cout << smpl1.refractive_index() << "\n";

}
