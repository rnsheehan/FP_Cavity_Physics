#ifndef FRESNEL_EQN_H
#define FRESNEL_EQN_H

// Implementation of the Fresnel equations for power reflection and transmission
// R. Sheehan 20 - 12 - 2018

class fresnel {
public:
	fresnel();
	fresnel(double lambda, material *m1, material *m2);

	void set_params(double lambda, material *m1, material *m2);

	double reflectance(double angle, bool polarisation); 

private:
	double theta_in; // angle of incidence = angle of reflection
	double theta_t; // angle of transmission
	double n1; // RI of external material
	double n2; // RI of internal material
	double nrat; // ratio n1/n2; 
	double wavelength; // operating wavelength, should be specified in units of um

	material *mat1; // external material RI = n1
	material *mat2; // internal material RI = n2
};

#endif
