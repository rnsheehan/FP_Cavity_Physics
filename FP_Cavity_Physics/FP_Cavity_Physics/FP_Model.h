#ifndef FP_MODEL_H
#define FP_MODEL_H

// The plan is to implement a class that describes an asymmetric Fabry-Perot etalon according to the model
// given in "Photonics" by Yariv and Yeh. The model includes mirrors with different reflectivities and can include 
// a loss / gain section. The model should return the FP spectrum based on the input parameters, as well as 
// characteristics of the FP cavity such as Frequency / Wavelength Free Spectral Range, Resonance Wavelengths, Resonance Widths
// R. Sheehan 25 - 10 - 2018

class fp_cavity {
public:
	fp_cavity();

private:
	int inside_mat; // this will describe the spectrally sensitive refractive index inside the cavity
	int outside_mat; // this will describe the spectrally sensitive refractive index outside the cavity

	double length; // total length of FP cavity
	double R1; // reflectivity mirror 1
	double R2; // reflectivity mirror 2
	double phase; // optical phase inside the cavity
	double alpha; // intensity loss coefficient (gain if alpha < 0)


};

#endif
