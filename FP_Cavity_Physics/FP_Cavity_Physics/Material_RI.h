#ifndef MATERIAL_MODELS_H
#define MATERIAL_MODELS_H

// Implementation of known Refractive Index models for various different materials
// Based on measured data or known formulae where available
// R. Sheehan 27 - 10 - 2015

// This could be implemented in a dynamic binding fashion
// R. Sheehan 24 - 2 - 2016

class ri_model{
public:
	ri_model(); 
	~ri_model(); 

	double ri_value(int mat_type, double wavelength, double xfr = 0.0, double yfr = 0.0); // user interacts with object through this method only

private:
	// ri based on a formula

	// Binaries
	double GaAs(double wavelength); 
	double AlAs(double wavelength); 
	double SiO2(double wavelength);
	double SiGeO2(double wavelength);

	// Ternaries
	double AlGaAs(double alfrac, double wavelength); 
	double InGaAs(double infrac, double wavelength);
	double GaAsN_1300(double nfrac);

	// Quaternaries
	double InGaAsP(double asfrac, double wavelength); 
	double InGaAsN(double infrac, double nfrac, double wavelength); 
	double AlInGaAs(double alfrac, double wavelength); 

	// ri based on data set
	double data_based_ri(int mat_type, std::string &filename, double wavelength); 

	void fill_material_list(); 
	
private:
	list<int> materials; // need a list of materials for which the ri data is known

	// If you want to sweep over a set of wavelength values you need to keep recorded data in memory
	// No sense in reading in data everytime a wavelength changes by 0.5 um if the material does not change
	int material_in_memory; 
	int n_data_points; 
	int n_cols; 
	bool data_in_memory; 
	vector<vector<double>> ri_data;
};
#endif