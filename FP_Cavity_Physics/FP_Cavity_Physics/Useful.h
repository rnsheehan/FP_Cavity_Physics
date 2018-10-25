#ifndef USEFUL_H
#define USEFUL_H

// Library of functions that are very useful
// R. Sheehan 4 - 7 - 2011

namespace useful_funcs{
	
	std::string TheTime();

	void exit_failure_output(std::string reason);

	void read_into_vector(std::string &filename, std::vector<double> &data, int &n_pts, bool loud = false); 

	//double test_func( double (*f)(int, int), int a, int b); 

}

#endif