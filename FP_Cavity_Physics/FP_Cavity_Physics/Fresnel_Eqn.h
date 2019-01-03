#ifndef FRESNEL_EQN_H
#define FRESNEL_EQN_H

// Implementation of the Fresnel equations for power reflection and transmission
// R. Sheehan 20 - 12 - 2018

class fresnel {
public:
	fresnel();	
	fresnel(double lambda, double indx_1, double indx_2); 

	//fresnel(double lambda, material *m1, material *m2);
	//void set_params(double lambda, material *m1, material *m2);
	//void update_wavelength(double lambda);
	
	void set_params(double lambda, double indx_1, double indx_2);
	
	double transmission_angle(double angle); 
	double reflectance(double angle, bool polarisation); 
	double transmittance(double angle, bool polarisation);

	inline double get_wl() { return wavelength;  }
	inline double get_n1() { return n1;  }
	inline double get_n2() { return n2; }
	inline double get_n_ratio() { return nrat;  }
	inline double get_critical_angle() { return theta_critical;}
	inline double get_brewster_angle() { return theta_brewster;}

private:
	double n1; // RI of external material
	double n2; // RI of internal material
	double wavelength; // operating wavelength, should be specified in units of um
	// wavelength not necessary for the calculation but is there to remind the user that RI is wavelength dependent

	double nrat; // ratio n2/n1, nrat > 1 => external reflection, nrat < 1 => internal reflection 
	double nrat_sqr; // square of RI ratio	

	// critical angle only relevant in the case of internal reflection n2 < n1, nrat < 1
	// internal reflection is real valued for angle < critical_angle, otherwise reflection is complex valued
	// nrat > 1 => external reflection, nrat < 1 => internal reflection
	double theta_critical; // critical angle
	
	// compute the brewster angle based on the input refractive index values
	// TM reflection = 0 for light incident at Brewster angle
	// Unpolarised light incident at Brewster angle is reflected in linearly polarised TE state
	// however, TE reflectivity is not great, generation of linearly polarised light in this manner is not efficient
	double theta_brewster; // brewster angle

	double theta_t; // angle of transmission computed using Snell's law

	// It is assumed that the light is propagating from material with RI = n1 into material with RI = n2
	// I'm going to change the class so that RI values can be input as parameters
	// This removes the need to have any material classes as members
	// The reason for this is that Ternary and Quaternary materials need extra input parameters that will make the code
	// awkward to implement if I was to include them. 
	//material *mat1; // external material RI = n1
	//material *mat2; // internal material RI = n2
};

#endif
