#ifndef MATERIALS_H
#define MATERIALS_H

// Code that implements a materials class new material classes can be added dynamically
// R. Sheehan 14 - 6 - 2016


// class that contains material interpolation coefficients that can be passed to interpolation function
class coefficients{
public:
	coefficients(); 

	void set_coeffs(double v1, double v2, double v3, double v4); 

	double planar_interpolation(double x, double y); 

	double p1, p2, p3, p4; 
}; 

// Declaration of base class material

class material{
public:

	material(); // Constructor

	// Declare pure virtual functions that will enable declaration of methods for RI and Eg
	// every class that inherits material must have an implementation of each of these functions
	//virtual double refractive_index(void)=0; 
	//virtual double bandgap_energy(void)=0;
	//virtual double bandgap_energy_with_strain(void)=0; 

	inline double get_wavelength(){return wavelength; } 
	inline double get_energy(){return energy; } 
	inline double get_temperature(){return temperature; } 

	inline void set_wavelength(double &val){wavelength = val; } // setter for wavelength
	inline void set_energy(double &val){energy = val; } // setter for lambda
	inline void set_temperature(double &val){temperature = val; } // setter for temperature

	//inline void convert_energy_to_wavelength(){wavelength = 1240.0 / energy;} // 

	void convert_energy_to_wavelength(); // convert a bandgap energy expressed in eV to a wavelength expressed in nm
	void convert_wavelength_to_energy(); // convert a wavelength expressed in nm to a bandgap energy expressed in eV
	void convert_C_K(double &Cvalue); // convert temperature in Celcius to Kelvin scale
	void convert_K_C(double &Kvalue); // convert temperature in Kelvin to Celcius scale

	double data_based_ri(std::string &filename); // ri based on data set

	// Methods needed to compute Strain Compensation in GaAsN and InGaAsN
	double hydrostatic_strain(double x, double y); 
	double shear_strain(double x, double y); 
	double stress_tensor(double x, double y); 
	double strain_correction(double x, double y); 
	double hh_lh_splitting(double x, double y);		

	// derived classes need access to these parameters
protected:
	double wavelength; // wavelength [ nm ] at which material properties will be computed
	double energy; // energy [ eV ] at which material properties will be computed
	double temperature; // temperature [ K ] at which material properties will be computed

	double delta; // spin split off energy 
	double LL; // stress tensor parameter

	coefficients C11; // elastic constants associated with stress tensor
	coefficients C12; // elastic constants associated with stress tensor
	coefficients alc; // material dependent lattice constant
	coefficients a; // deformation potential due to hydrostatic strain
	coefficients b; // deformation potential due to shear deformation

	// these parameters are only accessed by data_base_ri
private: 
	int n_data_points; 
	int n_cols; 
	bool data_in_memory; 
	vector<vector<double>> ri_data;
};

/****************************************************************************************************************************/
/*                                         Binaries                                                                         */
/****************************************************************************************************************************/

// Declaration of derived class GaAs

class GaAs : public material{
public:
	GaAs();
	GaAs(double wavelength); 

	double refractive_index(); 

	double bandgap_energy(); 

	//double bandgap_energy_with_strain(); 

private:
	
};

// Declaration of derived class AlAs

class AlAs : public material{
public:
	AlAs();
	AlAs(double wavelength); 

	double refractive_index(); 

	double bandgap_energy(); 

private:
	
};

// Declaration of derived class SiO2

class SiO2 : public material{
public:
	SiO2();
	SiO2(double wavelength); 

	double refractive_index(); 

	double bandgap_energy(); 

private:
	
};

// Declaration of derived class AlN

class AlN : public material{
public:
	AlN();
	AlN(double wavelength); 

	double refractive_index(); 

	double bandgap_energy(); 

private:
	
};

// Declaration of derived class AlSb

class AlSb : public material{
public:
	AlSb();
	AlSb(double wavelength); 

	double refractive_index(); 

	double bandgap_energy(); 

private:
	
};

// Declaration of derived class GaN

class GaN : public material{
public:
	GaN();
	GaN(double wavelength); 

	double refractive_index(); 

	double bandgap_energy(); 

private:
	
};

// Declaration of derived class GaSb

class GaSb : public material{
public:
	GaSb();
	GaSb(double wavelength); 

	double refractive_index(); 

	double bandgap_energy(); 

private:
	
};

// Declaration of derived class InP

class InP : public material{
public:
	InP();
	InP(double wavelength); 

	double refractive_index(); 

	double bandgap_energy(); 

private:
	
};

// Declaration of derived class InAs

class InAs : public material{
public:
	InAs();
	InAs(double wavelength); 

	double refractive_index(); 

	double bandgap_energy(); 

private:
	
};

/****************************************************************************************************************************/
/*                                         Ternaries                                                                        */
/****************************************************************************************************************************/

// Declaration of derived class AlGaAs

class AlGaAs : public material{
public:
	AlGaAs();
	AlGaAs(double wavelength); 

	double refractive_index(double alfrac); 

	double bandgap_energy(double alfrac); 

private:
	
};

// Declaration of derived class AlInAs

class AlInAs : public material{
public:
	AlInAs();
	AlInAs(double wavelength); 

	double refractive_index(double alfrac); 

	double bandgap_energy(double alfrac); 

private:
	
};

// Declaration of derived class InGaAs

class InGaAs : public material{
public:
	InGaAs();
	InGaAs(double wavelength); 

	double refractive_index(double infrac); 

	double bandgap_energy(double infrac); 

private:
	
};

// Declaration of derived class GaAsN

class GaAsN : public GaAs{
public:
	GaAsN();
	GaAsN(double wavelength); 

	double refractive_index(double nfrac); 

	double bandgap_energy(double nfrac); 

private:
	 
};

// Declaration of derived class InAsN

class InAsN : public InAs{
public:
	InAsN();
	InAsN(double wavelength); 

	double refractive_index(double nfrac); 

	double bandgap_energy(double nfrac); 

private:
	
};

/****************************************************************************************************************************/
/*                                         Quaternaries                                                                     */
/****************************************************************************************************************************/

// Declaration of derived class InGaAsP
class InGaAsP : public material{
public:
	InGaAsP();
	InGaAsP(double wavelength); 

	double refractive_index(double asfrac); 

	double bandgap_energy(double asfrac); 

private:
	
};

// Declaration of derived class AlInGaAs
class AlInGaAs : public material{
public:
	AlInGaAs();
	AlInGaAs(double wavelength); 

	double refractive_index(double asfrac); 

	double bandgap_energy(double asfrac); 

private:
	
};

// Declaration of derived class InGaAsN
class InGaAsN : public GaAsN, public InAsN{
public:
	InGaAsN();
	InGaAsN(double wavelength); 

	double refractive_index(double infrac, double nfrac); 

	double bandgap_energy(double infrac, double nfrac); 

	//double bandgap_energy_with_strain(double infrac, double nfrac); 

private:
	double bowing; // parameter associated with interpolation of InAsN and GaAsN
};

#endif