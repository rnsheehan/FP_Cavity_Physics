#ifndef ATTACH_H
#include "Attach.h"
#endif // !ATTACH_H

int main()
{
	std::string filename; 

	//filename = "Ag.txt";
	filename = "Silicon_Refractive_Index_Data.csv";

	int nc, nr; 
	std::vector<std::vector<double>> data; 

	vecut::read_into_matrix(filename, data, nr, nc, true); 

	std::cout << "Press enter to close.\n";
	std::cin.get(); 

	return 0; 
}
