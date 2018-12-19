#ifndef NEVILLE_H
#define NEVILLE_H

// File contains the definition of a namespace for doing interpolation based on Neville's algorithm
// Implementation is based on that given in NRinC
// R. Sheehan 19 - 12 - 2018

namespace interpolation {
	void polint(std::vector<double> &xa, std::vector<double> &ya, double x, double &y, double &dy); 
}

#endif

