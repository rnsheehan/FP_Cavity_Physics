#ifndef ATTACH_H
#include "Attach.h"
#endif

void testing::matrix_file_IO()
{
	// check that file IO for reading into 2D matrix is working correctly
	// R. Sheehan 18 - 12 - 2018

	std::string filename;

	filename = "Ag.txt";
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