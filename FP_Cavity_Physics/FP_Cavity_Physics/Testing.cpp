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

	InGaAs smpl3;

	mat1 = &smpl1; 
	mat2 = &smpl3; 

	double wavelength = 1.55; 

	mat1->set_wavelength(wavelength);
	mat2->set_wavelength(wavelength); 

	std::cout<<mat1->refractive_index()<<"\n"; // this should return the RI of Si

	smpl1.set_wavelength(wavelength);
	std::cout << smpl1.refractive_index() << "\n\n";

	std::cout << mat2->refractive_index() << "\n"; // this should return the RI of SiO2
	std::cout << mat2->refractive_index(0.5) << "\n"; 

	smpl2.set_wavelength(wavelength);
	std::cout << smpl2.refractive_index() << "\n";

}

void testing::fresnel_values()
{
	// run some tests on the Fresnel equations

	double wavelength = 1.55;

	// Declarate some material objects

	Air smpl1;
	SiO2 smpl2;

	/*InP smpl2; 
	AlN smpl1;*/ 

	/*material *mat1;
	material *mat2;

	mat1 = &smpl1; mat2 = &smpl2;

	mat1->set_wavelength(wavelength);
	mat2->set_wavelength(wavelength);*/

	smpl1.set_wavelength(wavelength);  smpl2.set_wavelength(wavelength);

	fresnel calc; 

	calc.set_params(wavelength, smpl1.refractive_index(), smpl2.refractive_index()); 

	std::cout << "Power Reflection using Fresnel Equations\nwavelength: " << calc.get_wl() << " um\n"; 
	std::cout << "n1: " << calc.get_n1() << "\nn2: " << calc.get_n2() << "\n"; 
	std::cout << "RI ratio: " << calc.get_n_ratio() << "\n";
	std::cout << "Critical angle: " << calc.get_critical_angle()<<" rad, "<< calc.get_critical_angle()*RAD_TO_DEG << " deg\n";
	std::cout << "Brewster angle: " << calc.get_brewster_angle() << " rad, " << calc.get_brewster_angle()*RAD_TO_DEG << " deg\n\n";

	double angle = 0.0*DEG_TO_RAD; 

	std::cout << "Angle of incidence: " << angle << " rad, "<< angle * RAD_TO_DEG<<" deg\n";
	std::cout << "Angle of transmission: " << calc.transmission_angle(angle) << " rad, " << calc.transmission_angle(angle)*RAD_TO_DEG << " deg\n"; 
	std::cout << "TE reflectance: " << calc.reflectance(angle, TE) << ", TE transmittance: " << calc.transmittance(angle, TE) <<"\n";
	std::cout << "TM reflectance: " << calc.reflectance(angle, TM) << ", TM transmittance: "<< calc.transmittance(angle, TM) << "\n\n";
}

void testing::fp_test()
{
	// start testing the FP cavity implementation

	double incident_angle, cav_length, loss_fac, wl_start, wl_end;

	incident_angle = 0.0; cav_length = 17; loss_fac = 0.0; wl_start = 1.4; wl_end = 1.6; 

	/*InP smpl2;
	AlN smpl1;*/

	Air smpl1;
	Si smpl2;

	material *mat1;
	material *mat2;

	mat1 = &smpl1; mat2 = &smpl2;

	fp_cavity etalon(incident_angle, cav_length, loss_fac, wl_start, wl_end, mat1, mat2);

	etalon.compute_spectrum(true); 
}
