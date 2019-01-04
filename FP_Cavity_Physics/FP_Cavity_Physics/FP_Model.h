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
	fp_cavity(double incident_angle, double cav_length, double loss_fac, double wl_start, double wl_end, material *m1, material *m2);

	void set_params(double incident_angle, double cav_length, double loss_fac, double wl_start, double wl_end, material *m1, material *m2);

	void compute_spectrum(bool loud = false); 

private:
	void dispersion(double wavelength); // compute change in RI and R with wavelength
	double finesse(double wavelength); // compute the FP cavity finesse
	double airy(double F); // compute the FP cavity Airy function	
	void reflection(double wavelength, double &F, double &A, double &RFP, bool loud = false); // compute the FP cavity reflection

private:
	bool params_defined; 
	int nwl; // number of wavelength subdivisions
	double angle; // incident angle for light propagating into FP cavity
	double length; // total length of FP cavity
	double wl1; // initial wavelength for computing the FP spectrum
	double wl2; // final wavelength for computing the FP spectrum
	double dwl; // wavelength spacing
	double n1; // RI of external material
	double n2; // RI of internal material
	double r; // reflectivity n1 -> n2
	double rpr; // reflectivity n2 -> n1
	double R; // cavity mirror reflectance
	double phase; // optical phase inside the cavity
	double alpha; // intensity loss coefficient (gain if alpha < 0)

	material *mat1; // external material RI = n1
	material *mat2; // internal material RI = n2
};

#endif
