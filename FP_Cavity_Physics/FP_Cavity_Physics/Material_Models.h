#ifndef MATERIAL_MODELS_H
#define MATERIAL_MODELS_H

// Code that implements a materials class new material classes can be added dynamically
// R. Sheehan 14 - 6 - 2016

// Declaration of base class material

class material {
public:
	material(); // Constructor

	// Declare pure virtual functions that will enable declaration of methods for RI and Eg
	// every class that inherits material must have an implementation of each of these functions
	//virtual double refractive_index(void)=0; 
	//virtual double bandgap_energy(void)=0;
	//virtual double bandgap_energy_with_strain(void)=0; 

	inline double get_wavelength() { return wavelength; }
	inline double get_energy() { return energy; }
	inline double get_temperature() { return temperature; }

	inline void set_wavelength(double &val) { wavelength = val; } // setter for wavelength
	inline void set_energy(double &val) { energy = val; } // setter for lambda
	inline void set_temperature(double &val) { temperature = val; } // setter for temperature

	//inline void convert_energy_to_wavelength(){wavelength = 1240.0 / energy;} // 

	void convert_energy_to_wavelength(); // convert a bandgap energy expressed in eV to a wavelength expressed in nm
	void convert_wavelength_to_energy(); // convert a wavelength expressed in nm to a bandgap energy expressed in eV
	void convert_C_K(double &Cvalue); // convert temperature in Celcius to Kelvin scale
	void convert_K_C(double &Kvalue); // convert temperature in Kelvin to Celcius scale

	double data_based_ri(std::string &filename); // ri based on data set

	// derived classes need access to these parameters
protected:
	double wavelength; // wavelength [ nm ] at which material properties will be computed
	double energy; // energy [ eV ] at which material properties will be computed
	double temperature; // temperature [ K ] at which material properties will be computed
	
	// these parameters are only accessed by data_base_ri
private:
	int n_data_points;
	int n_cols;
	bool data_in_memory;
	std::vector< std::vector< double > > ri_data;
};

#endif
