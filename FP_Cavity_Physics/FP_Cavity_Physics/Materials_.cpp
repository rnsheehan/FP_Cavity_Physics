#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of coefficients class
coefficients::coefficients()
{
	p1 = p2 = p3 = p4 = 0.0; 
}

void coefficients::set_coeffs(double v1, double v2, double v3, double v4)
{
	p1 = v1; p2 = v2; p3 = v3; p4 = v4; 
}

double coefficients::planar_interpolation(double x, double y)
{
	// Compute the value of a planar interpolation based on the values of the members in coefficients

	double t1 = p1 * x * y; 
	double t2 = p2 * y * ( 1 - x );	
	double t3 = p3 * x * ( 1 - y );
	double t4 = p4 * ( 1- x ) * ( 1 - y );

	return ( t1 + t2 + t3 + t4 ); 
}

// Definition of material class
material::material()
{
	// Default constructor

	data_in_memory = false; 

	n_data_points = n_cols = 0; 

	temperature = 298.15; // Assume that temperature is 25 ( C ) unless otherwise stated

	energy = wavelength = 0.0; 

	//delta = LL = 0.0; // material dependent parameters

	LL = 5.6533; // lattice constant parameter
	delta = 17e-3; // spin-split off energy

	// define interpolation coefficients
	C11.set_coeffs(18.7, 8.33, 29.3, 11.9); 
	C12.set_coeffs(12.5, 4.53, 15.9, 5.38); 
	alc.set_coeffs(4.98, 6.0583, 4.5, 5.6533);
	a.set_coeffs(-0.35, -6.08, 3.0, -8.33);
	b.set_coeffs(-1.2, -1.8, -2.2, -2.0); 
}

void material::convert_energy_to_wavelength()
{
	// convert a bandgap energy expressed in eV to a wavelength expressed in nm
	// R. Sheehan 14 - 6 - 2016

	try{		
		if(energy > 0.0){
			wavelength = 1240.0 / energy; 
		}
		else{
			std::string reason = "Value stored in energy is not correct\n"; 
			reason += "energy = " + template_funcs::toString(energy, 2) + "\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void material::convert_wavelength_to_energy()
{
	// convert a wavelength expressed in nm to a bandgap energy expressed in eV
	// R. Sheehan 14 - 6 - 2016

	try{		
		if(wavelength > 0.0){
			energy = 1240.0 / wavelength; 
		}
		else{
			std::string reason = "Value stored in energy is not correct\n"; 
			reason += "energy = " + template_funcs::toString(energy, 2) + "\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void material::convert_C_K(double &Cvalue)
{
	// convert temperature in Celcius to Kelvin scale
	// R. Sheehan 14 - 6 - 2016

	try{		
		if(Cvalue > 0.0){
			temperature = Cvalue + 273.15; // temperature is expressed in units of K
		}
		else{
			std::string reason = "Value stored in Cvalue is not correct\n"; 
			reason += "Cvalue = " + template_funcs::toString(Cvalue, 2) + "\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void material::convert_K_C(double &Kvalue)
{
	// convert temperature in Kelvin to Celcius scale
	// R. Sheehan 14 - 6 - 2016

	try{		
		if(Kvalue > 0.0){
			temperature = Kvalue - 273.15; // temperature is expressed in units of C
		}
		else{
			std::string reason = "Value stored in Cvalue is not correct\n"; 
			reason += "Cvalue = " + template_funcs::toString(Kvalue, 2) + "\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double material::data_based_ri(string &filename)
{
	// this method reads in the data stored in file filename and performs an interpolation over the data set at position wavelength
	// source for data is one of the following
	// http://refractiveindex.info/ // excel csv files
	// http://www.filmetrics.com/refractive-index-database/ // notepad csv files
	// I'm going to work with the filmetrics data sets as they are in a more convenient format, data is the same in both
	
	try{
		wavelength *= 1000.0; // convert to nm because wavelength data is stored in nanometres

		int wavepos; 

		if(!data_in_memory){
			// read in new data only if no data present or material type has changed

			//cout<<"\nReading data in from memory\n"; 
			
			ri_data.clear(); 

			n_data_points = n_cols = 0; 

			// read the data from the file
			array_funcs::vvd_read_matrix_from_file(filename, ri_data, n_data_points, n_cols); 

			// first line contains no text so actual n_data_points is less than measured from file
			n_data_points -= 1;	

			data_in_memory = true; 
		} 

		// Program will close before this point if file cannot be found

		// store the wavelength values in the array vals for searching
		double *waves = new(double [n_data_points+1]); 

		for(int i=1; i<=n_data_points; i++){
			*(waves+i) = ri_data[i+1][1]; 
		}

		// find the position in vals that is closest to wavelength
		wavepos = useful_funcs::binary_search(waves, n_data_points, wavelength); 
	
		// perform interpolation on ri_data
		if(wavepos != -1){

			double ri_value = 0.0; 
			double delta_ri_value = 0.0; 
		
			// store the wavelength values in the array vals for searching
			double *vals = new(double [n_data_points+1]); 

			for(int i=1; i<=n_data_points; i++){
				*(vals+i) = ri_data[i+1][2]; // store refractive index values
			}

			useful_funcs::polint(waves, vals, n_data_points, wavelength, ri_value, delta_ri_value); 

			delete[] vals; 

			return ri_value; 
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in material::data_based_ri(string &filename)";
			reason += "wavelength = " + template_funcs::toString(wavelength, 3) + " (um) is not found in  file = " + filename + " \n"; 
			throw std::range_error(reason); 
		}

		delete[] waves; 

		// Under what circumstances do you de-allocate ri_data? 
		//delete[] ri_data
		//data_in_memory = false; 
		//material_in_memory = -1; 
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double material::hydrostatic_strain(double x, double y)
{
	// compute the material dependent hydrostatic-strain H(x, y)

	try{
		
		double c11_val = C11.planar_interpolation(x, y); 

		if( fabs(c11_val) > EPS ){

			double t = ( 2.0*( c11_val - C12.planar_interpolation(x, y) ) ) / c11_val; 

			double H = a.planar_interpolation(x, y) * t * stress_tensor(x, y); 

			return H; 
		}
		else{
			std::string reason;
			reason = "C_{11} = " + template_funcs::toString(c11_val) + 
				" is equal to zero\nmaterial::hydrostatic_strain(double x, double y)\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	} 
}

double material::shear_strain(double x, double y)
{
	// compute the material dependent shear-strain S(x, y)

	try{
		
		double c11_val = C11.planar_interpolation(x, y); 

		if( fabs(c11_val) > EPS ){

			double t = ( c11_val + 2.0 * C12.planar_interpolation(x, y) ) / c11_val; 

			double S = b.planar_interpolation(x, y) * t * stress_tensor(x, y); 

			return S; 
		}
		else{
			std::string reason;
			reason = "C_{11} = " + template_funcs::toString(c11_val) + 
				" is equal to zero\nmaterial::shear_strain(double x, double y)\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	} 
}

double material::stress_tensor(double x, double y)
{
	// compute the material dependent stress tensor \sigma(x, y)

	try{
		
		double lattice_constant = alc.planar_interpolation(x, y); 

		if( lattice_constant > 0.0 ){

			return (-1.0 + ( LL / lattice_constant)); 
		}
		else{
			std::string reason;
			reason = "lattice constant = " + template_funcs::toString(lattice_constant) + 
				" is not positive in\nmaterial::stress_tensor(double x, double y)\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double material::strain_correction(double x, double y)
{
	// total correction due to strain compensation H - S

	return ( hydrostatic_strain(x, y) - shear_strain(x, y) ); 
}

double material::hh_lh_splitting(double x, double y)
{
	// correction due to hh-lf splitting E_{hh-lh}

	double s_sqr = template_funcs::DSQR( shear_strain(x, y) ); 

	return ( 2.0 * s_sqr - ( s_sqr / delta ) ); 
}

/****************************************************************************************************************************/
/*                                         Binaries                                                                         */
/****************************************************************************************************************************/

// Definition of class GaAs
GaAs::GaAs()
{
	// Default Constructor
}

GaAs::GaAs(double wavelength)
{
	set_wavelength(wavelength); 
}

double GaAs::refractive_index()
{
	// source: http://www.batop.com/information/n_GaAs.html
    // type: fit of experimental data
    // model valid for at least lambda (um) in [0.65, 1.8]

	try{
		// Model is only valid on certain wavelength range

		if( wavelength > 0.64 && wavelength < 1.9){

			double A = 8.95; 
			double B = 2.054;
			double Csqr = 0.39; 

			double t1 = ( 1.0 - ( Csqr / template_funcs::DSQR( wavelength ) ) ); 
			double t2 = ( A + (B/t1) ); 

			return sqrt(t2); 
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in GaAs::refractive_index()\n";
			reason += "wavelength = " + template_funcs::toString( wavelength , 3) + " (um)\n"; 
			reason += "Allowed range for GaAs model is [ 0.65, 1.8 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double GaAs::bandgap_energy()
{
	// model for the band-gap energy of the material GaAs
    // source: http://www.batop.com/information/Eg_GaAs.html
    // energy is expressed in eV, model valid for T in [0, 700] Kelvin

	try{
		// Model is only valid on certain temperature range

		if(get_temperature() > 0.0 && get_temperature() < 700.0){

			double c1 = 1.519;
			double c2 = 5.408e-4;
			double c3 = 204.0;
			double t1 = c2 * template_funcs::DSQR( get_temperature() );
			double t2 = ( get_temperature() + c3 );

			return ( c1 - (t1/t2) );
		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in GaAs::bandgap_energy()\n";
			reason += "temperature = " + template_funcs::toString(get_temperature() , 3) + " (um)\n"; 
			reason += "Allowed range for GaAs model is [ 0, 700 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Definition of class AlAs
AlAs::AlAs()
{
	// Default Constructor
}

AlAs::AlAs(double wavelength)
{
	set_wavelength(wavelength); 
}

double AlAs::refractive_index()
{
	// source: http://www.batop.com/information/n_AlAs.html
    // type: theoretical fit
    // model valid for at least lambda (um) in [0.45, 2.0]

	try{
		// Model is only valid on certain wavelength range

		if( wavelength > 0.44 && wavelength < 2.1){

			// PE = HC/wavelength # photon energy in Joules
			double A0 = 25.3; // constant, determined by fitting with experimental data
			double B0 = -0.8; // constant, determined by fitting with experimental data
			double E0 = 2.95; // fundamental band gap at G-point expressed in eV
			double S0 = 3.25; // fundamental band gap at G-point plus D0 spin-orbit splitting energy expressed in eV
			//double ER = math.pow((E0/S0),1.5) // (E0 / E0+D0)^{3/2}
			double ER = pow((E0/S0),1.5); 

			double chi = HC/(wavelength*E0); 
			double chisqr = template_funcs::DSQR(chi); 

			double t1 = sqrt(1.0+chi); 
			double t2 = sqrt(1.0-chi);
			double fchi = (2.0-t1-t2)/chisqr;

			double chiS0 = HC/(wavelength*S0);
			double chisqrS0 = template_funcs::DSQR(chiS0);

			double t1S0 = sqrt(1.0+chiS0);
			double t2S0 = sqrt(1.0-chiS0);
			double fchiS0 = (2.0-t1S0-t2S0)/chisqrS0;

			double t3 = ( fchi + 0.5*fchiS0*ER );

			return sqrt( A0*t3 + B0 ); 	
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in AlAs::refractive_index()\n";
			reason += "wavelength = " + template_funcs::toString( wavelength , 3) + " (um)\n"; 
			reason += "Allowed range for AlAs model is [ 0.45, 2.0 ]\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

double AlAs::bandgap_energy()
{
	// model for the band-gap energy of the material AlAs
    // source: http://www.batop.com/information/Eg_AlGaAs.html

	try{
		// Model is only valid on certain temperature range

		if(get_temperature() > 0.0 && get_temperature() < 700.0){

			return 2.16; // this is the value of the bandgap at T = 300 (K)

		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in AlAs::bandgap_energy()\n";
			reason += "temperature = " + template_funcs::toString(get_temperature() , 3) + " (um)\n"; 
			reason += "Allowed range for AlAs model is [ 0, 700 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Definition of class SiO2
SiO2::SiO2()
{
	// Default Constructor
}

SiO2::SiO2(double wavelength)
{
	set_wavelength(wavelength); 
}

double SiO2::refractive_index()
{
	// Refractive Index of Pure Silica
	// Model is based on a Sellmeier polynomial with experimentally derived coefficients
	// Taken from K. Okamoto, "Fundamanetals of Optical Waveguides", 2006
	// wavelength is input in units of microns
	// R. Sheehan 12 - 4 - 2012

	try{
		// Model is only valid on certain wavelength range

		if( wavelength > 0.4 && wavelength < 2.5){

			int ncoeffs = 3;

			double a[4]; 
			double b[4]; 
	
			// Taken from K. Okamoto, "Fundamanetals of Optical Waveguides", 2006
			// Valid wavelength range 0.4 < \lambda < 2.5
			a[1]=0.6965325; a[2]=0.4083099; a[3]=0.8968766; 
			b[1]=4.368309e-3; b[2]=1.394999e-2; b[3]=97.93399;

			// Taken from refractiveindex.info/?group=GLASSES&material=F_SILICA
			// Valid wavelength range 0.21 < \lambda < 3.71
			/*a[1]=0.6961663; a[2]=0.4079426; a[3]=0.8974794; 
			b[1]=DSQR(0.0684043); b[2]=DSQR(0.1162414); b[3]=DSQR(9.896161);*/
	
			double nsqr=1.0;
			double lsqr = template_funcs::DSQR( wavelength );
	
			for(int i=1; i<=ncoeffs; i++){
				nsqr +=  ( a[i]*lsqr ) / ( lsqr - b[i] ) ;
			}
	
			return sqrt(nsqr);	
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in GaAs::refractive_index()\n";
			reason += "wavelength = " + template_funcs::toString( wavelength , 3) + " (um)\n"; 
			reason += "Allowed range for SiO2 model is [ 0.65, 1.8 ]\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

double SiO2::bandgap_energy()
{
	// model for the band-gap energy of the material SiO2
    // source: https://en.wikipedia.org/wiki/Band_gap, which references http://journals.aps.org/prb/abstract/10.1103/PhysRevB.83.174201

	try{

		if(get_temperature() > 0.0 && get_temperature() < 700.0){
			return 9.0; 
		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in SiO2::bandgap_energy()\n";
			reason += "temperature = " + template_funcs::toString(get_temperature() , 3) + " (um)\n"; 
			reason += "Allowed range for SiO2 model is [ 0, 700 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Definition of class AlN
AlN::AlN()
{
	// Default Constructor
}

AlN::AlN(double wavelength)
{
	set_wavelength(wavelength); 
}

double AlN::refractive_index()
{
	// Refractive Index of AlN
	// source: http://www.filmetrics.com/refractive-index-database/

	try{
		// Model is only valid on certain wavelength range

		if( wavelength > 0.0 && wavelength < 3.0){

			std::string ri_file = "AlN.txt"; 

			double ret_val = data_based_ri(ri_file); 
	
			return ret_val;	
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in GaAs::refractive_index()\n";
			reason += "wavelength = " + template_funcs::toString( wavelength , 3) + " (um)\n"; 
			reason += "Allowed range for GaAs model is [ 0.65, 1.8 ]\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

double AlN::bandgap_energy()
{
	// model for the band-gap energy of the material SiO2
    // source: http://www.ioffe.ru/SVA/NSM/Semicond/AlN/index.html

	try{

		if(get_temperature() > 0.0 && get_temperature() < 700.0){
			return 6.2; 
		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in SiO2::bandgap_energy()\n";
			reason += "temperature = " + template_funcs::toString(get_temperature() , 3) + " (um)\n"; 
			reason += "Allowed range for AlAs model is [ 0, 700 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Definition of class AlSb
AlSb::AlSb()
{
	// Default Constructor
}

AlSb::AlSb(double wavelength)
{
	set_wavelength(wavelength); 
}

double AlSb::refractive_index()
{
	// Refractive Index of AlSb
	// source: http://www.filmetrics.com/refractive-index-database/

	try{
		// Model is only valid on certain wavelength range

		if( wavelength > 0.0 && wavelength < 3.0){

			std::string ri_file = "AlSb.txt"; 

			double ret_val = data_based_ri(ri_file); 
	
			return ret_val;	
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in AlSb::refractive_index()\n";
			reason += "wavelength = " + template_funcs::toString( wavelength , 3) + " (um)\n"; 
			reason += "Allowed range for AlSb model is [ 0.65, 1.8 ]\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

double AlSb::bandgap_energy()
{
	// model for the band-gap energy of the material AlSb
    // source: http://www.semiconductors.co.uk/propiiiv5653.htm

	try{

		if(get_temperature() > 0.0 && get_temperature() < 700.0){
			return 1.615; 
		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in AlSb::bandgap_energy()\n";
			reason += "temperature = " + template_funcs::toString(get_temperature() , 3) + " (um)\n"; 
			reason += "Allowed range for AlSb model is [ 0, 700 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Definition of class GaN
GaN::GaN()
{
	// Default Constructor
}

GaN::GaN(double wavelength)
{
	set_wavelength(wavelength); 
}

double GaN::refractive_index()
{
	// Refractive Index of GaN
	// source: http://www.filmetrics.com/refractive-index-database/

	try{
		// Model is only valid on certain wavelength range

		if( wavelength > 0.0 && wavelength < 3.0){

			std::string ri_file = "GaN.txt"; 

			double ret_val = data_based_ri(ri_file); 
	
			return ret_val;	
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in GaN::refractive_index()\n";
			reason += "wavelength = " + template_funcs::toString( wavelength , 3) + " (um)\n"; 
			reason += "Allowed range for GaN model is [ 0.65, 1.8 ]\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

double GaN::bandgap_energy()
{
	// model for the band-gap energy of the material GaN
    // source: http://www.ioffe.ru/SVA/NSM/Semicond/GaN/index.html
	try{

		if(get_temperature() > 0.0 && get_temperature() < 700.0){
			return 3.4; 
		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in GaN::bandgap_energy()\n";
			reason += "temperature = " + template_funcs::toString(get_temperature() , 3) + " (um)\n"; 
			reason += "Allowed range for GaN model is [ 0, 700 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Definition of class GaSb
GaSb::GaSb()
{
	// Default Constructor
}

GaSb::GaSb(double wavelength)
{
	set_wavelength(wavelength); 
}

double GaSb::refractive_index()
{
	// Refractive Index of GaSb
	// source: http://www.filmetrics.com/refractive-index-database/

	try{
		// Model is only valid on certain wavelength range

		if( wavelength > 0.0 && wavelength < 3.0){

			std::string ri_file = "GaSb.txt"; 

			double ret_val = data_based_ri(ri_file); 
	
			return ret_val;	
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in GaSb::refractive_index()\n";
			reason += "wavelength = " + template_funcs::toString( wavelength , 3) + " (um)\n"; 
			reason += "Allowed range for GaSb model is [ 0.65, 1.8 ]\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

double GaSb::bandgap_energy()
{
	// model for the band-gap energy of the material GaSb
    // source: http://www.ioffe.ru/SVA/NSM/Semicond/GaSb/index.html
	try{

		if(get_temperature() > 0.0 && get_temperature() < 700.0){
			return 0.726; 
		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in GaSb::bandgap_energy()\n";
			reason += "temperature = " + template_funcs::toString(get_temperature() , 3) + " (um)\n"; 
			reason += "Allowed range for GaSb model is [ 0, 700 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Definition of class InP
InP::InP()
{
	// Default Constructor
}

InP::InP(double wavelength)
{
	set_wavelength(wavelength); 
}

double InP::refractive_index()
{
	// Refractive Index of InP
	// source: http://www.filmetrics.com/refractive-index-database/

	try{
		// Model is only valid on certain wavelength range

		if( wavelength > 0.0 && wavelength < 3.0){

			std::string ri_file = "InP.txt"; 

			double ret_val = data_based_ri(ri_file); 
	
			return ret_val;	
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in InP::refractive_index()\n";
			reason += "wavelength = " + template_funcs::toString( wavelength , 3) + " (um)\n"; 
			reason += "Allowed range for InP model is [ 0.65, 1.8 ]\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

double InP::bandgap_energy()
{
	// model for the band-gap energy of the material GaSb
    // source: http://www.ioffe.ru/SVA/NSM/Semicond/InP/index.html
	try{

		if(get_temperature() > 0.0 && get_temperature() < 700.0){
			return 1.344; 
		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in InP::bandgap_energy()\n";
			reason += "temperature = " + template_funcs::toString(get_temperature() , 3) + " (um)\n"; 
			reason += "Allowed range for InP model is [ 0, 700 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Definition of class InAs
InAs::InAs()
{
	// Default Constructor
}

InAs::InAs(double wavelength)
{
	set_wavelength(wavelength); 
}

double InAs::refractive_index()
{
	// Refractive Index of InP
	// source: http://www.filmetrics.com/refractive-index-database/

	try{
		// Model is only valid on certain wavelength range

		if( wavelength > 0.0 && wavelength < 3.0){

			std::string ri_file = "InAs.txt"; 

			double ret_val = data_based_ri(ri_file); 
	
			return ret_val;	
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in InAs::refractive_index()\n";
			reason += "wavelength = " + template_funcs::toString( wavelength , 3) + " (um)\n"; 
			reason += "Allowed range for InAs model is [ 0.65, 1.8 ]\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

double InAs::bandgap_energy()
{
	// model for the band-gap energy of the material InAs
    // source: http://www.ioffe.ru/SVA/NSM/Semicond/InAs/index.html
	try{

		if(get_temperature() > 0.0 && get_temperature() < 700.0){
			return 0.355; 
		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in InAs::bandgap_energy()\n";
			reason += "temperature = " + template_funcs::toString(get_temperature() , 3) + " (um)\n"; 
			reason += "Allowed range for InAs model is [ 0, 700 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

/****************************************************************************************************************************/
/*                                         Ternaries                                                                        */
/****************************************************************************************************************************/

// Definition of class AlGaAs
AlGaAs::AlGaAs()
{
	// Default Constructor
}

AlGaAs::AlGaAs(double wavelength)
{
	set_wavelength(wavelength); 
}

double AlGaAs::refractive_index(double alfrac)
{
	// source: http://www.batop.com/information/n_AlGaAs.html
    // type: theoretical fit
    // model valid for at least lambda (um) in [0.65, 2.0]
    // model describes refractive index of Al_{x}Ga_{1-x}As
    // As alfrac -> 0 Al_{x}Ga_{1-x}As -> GaAs
    // As alfrac -> 1 Al_{x}Ga_{1-x}As -> AlAs

	try{
		// Model is only valid on certain wavelength range

		if( wavelength > 0.64 && wavelength < 2.1 && alfrac >= 0.0 && alfrac <= 1.0){

			// PE = HC/wavelength # photon energy in Joules
			double A0 = 6.0 + (19.0*alfrac); // fitting parameter, value dependent on Alfrac
			double B0 = 9.4 - (10.2*alfrac); // fitting parameter, value dependent on Alfrac
			double alfracsqr = template_funcs::DSQR(alfrac);
			double E0 = 1.425 + (1.155*alfrac) + (0.37*alfracsqr); // fundamental band gap at G-point expressed in eV
			double S0 = 1.765 + (1.115*alfrac) + (0.37*alfracsqr); // fundamental band gap at G-point plus D0 spin-orbit splitting energy expressed in eV
			double ER = pow((E0/S0),1.5); // (E0 / E0+D0)^{3/2}

			double chi = HC/(wavelength*E0);
			double chisqr = template_funcs::DSQR(chi);

			double t1 = sqrt(1.0+chi);
			double t2 = sqrt(1.0-chi);
			double fchi = (2.0-t1-t2)/chisqr;

			double chiS0 = HC/(wavelength*S0);
			double chisqrS0 = template_funcs::DSQR(chiS0);

			double t1S0 = sqrt(1.0+chiS0);
			double t2S0 = sqrt(1.0-chiS0);
			double fchiS0 = (2.0-t1S0-t2S0)/chisqrS0;

			double t3 = ( fchi + 0.5*fchiS0*ER );

			return sqrt(A0*t3+B0); 
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in AlGaAs::refractive_index()\n";
			reason += "wavelength = " + template_funcs::toString( wavelength , 3) + " (um)\n"; 
			reason += "Allowed range for GaAs model is [ 0.65, 2.0 ]\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

double AlGaAs::bandgap_energy(double alfrac)
{
	// model for the band-gap energy of the material Al_{x}Ga_{1-x}As
    // source: http://www.batop.com/information/Eg_AlGaAs.html
    // energy is expressed in eV, model valid for x in [0, 0.45] which is a direct BG material, for [0.45, 1] AlGaAs has an indirect BG
    // at a temperature of T = 300 K
    // As alfrac -> 0 Al_{x}Ga_{1-x}As -> GaAs
    // As alfrac -> 1 Al_{x}Ga_{1-x}As -> AlAs

	try{
		// Model is only valid on certain alfrac range

		if(alfrac >= 0.0 && alfrac <= 0.45){

			return (1.422 + 1.2475*alfrac);
		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in AlGaAs::bandgap_energy()\n";
			reason += "alfrac = " + template_funcs::toString( alfrac , 3) + " (um)\n"; 
			reason += "Allowed range for AlGaAs model is [ 0, 0.45 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Definition of class AlInAs
AlInAs::AlInAs()
{
	// Default Constructor
}

AlInAs::AlInAs(double wavelength)
{
	set_wavelength(wavelength); 
}

double AlInAs::refractive_index(double alfrac)
{
	// no model available yet

	try{
		// Model is only valid on certain wavelength range

		if( wavelength > 0.64 && wavelength < 2.1 && alfrac >= 0.0 && alfrac <= 1.0){
			return 0.0; 
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in AlInAs::refractive_index()\n";
			reason += "wavelength = " + template_funcs::toString( wavelength , 3) + " (um)\n"; 
			reason += "Allowed range for AlInAs model is [ 0.65, 2.0 ]\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

double AlInAs::bandgap_energy(double alfrac)
{
	// model for the band-gap energy of the material Al_{x}In_{1-x}As
    // source: Chuang
    // energy is expressed in eV, model valid for x in [0, 1.0]
    // at a temperature of T = 300 K
    // As alfrac -> 0 Al_{x}In_{1-x}As -> InAs
    // As alfrac -> 1 Al_{x}In_{1-x}As -> AlAs

	try{
		// Model is only valid on certain alfrac range

		if(alfrac >= 0.0 && alfrac <= 1.0){

			return (0.36 + 2.35*alfrac + 0.24*template_funcs::DSQR(alfrac) );

		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in AlInAs::bandgap_energy()\n";
			reason += "alfrac = " + template_funcs::toString( alfrac , 3) + " (um)\n"; 
			reason += "Allowed range for AlInAs model is [ 0, 1.0 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Definition of class InGaAs
InGaAs::InGaAs()
{
	// Default Constructor
}

InGaAs::InGaAs(double wavelength)
{
	set_wavelength(wavelength); 
}

double InGaAs::refractive_index(double infrac)
{
	// source: http://www.batop.com/information/n_InGaAs.html
    // type: theoretical fit
    // model valid for at least lambda (um) in [1.15, 2.0]
    // model describes refractive index of In_{x}Ga_{1-x}As
    // As infrac -> 0 In_{x}Ga_{1-x}As -> GaAs
    // As infrac -> 1 In_{x}Ga_{1-x}As -> InAs

	try{
		// Model is only valid on certain wavelength range

		if( wavelength > 1.14 && wavelength < 2.1 && infrac >= 0.0 && infrac <= 1.0){

			double A = 8.950; // empirical coefficient
			double B = 2.054; // empirical coefficient
			double C = 0.6245; // empirical coefficient
			double EgGaAs = 1.424; // fundamental band gap of GaAs at room temperature (300 K)

			double infracsqr = template_funcs::DSQR(infrac);
			double Eg = EgGaAs - (1.501*infrac) + (0.436*infracsqr);

			double t1 = ( C * EgGaAs)/( wavelength * Eg );
			double t2 = 1.0 - template_funcs::DSQR(t1);
			double t3 = B/t2;
			double t4 = A+t3;

			return sqrt(t4); 
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in InGaAs::refractive_index()\n";
			reason += "wavelength = " + template_funcs::toString( wavelength , 3) + " (um)\n"; 
			reason += "Allowed range for GaAs model is [ 1.15, 2.0 ]\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

double InGaAs::bandgap_energy(double infrac)
{
	// model for the band-gap energy of the material In_{x}Ga_{1-x}As
    // source: http://www.batop.com/information/Eg_InGaAs.html
    // energy is expressed in eV, model valid for x in [0, 1.0]
    // at a temperature of T = 300 K
    // As infrac -> 0 In_{x}Ga_{1-x}As -> GaAs
    // As infrac -> 1 In_{x}Ga_{1-x}As -> InAs

	try{
		// Model is only valid on certain infrac range

		if(infrac >= 0.0 && infrac <= 1.0){

			return ( 1.425 - 1.501*infrac + 0.436*template_funcs::DSQR(infrac)  );
		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in InGaAs::bandgap_energy()\n";
			reason += "infrac = " + template_funcs::toString( infrac , 3) + " (um)\n"; 
			reason += "Allowed range for InGaAs model is [ 0, 1 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Definition of class GaAsN
GaAsN::GaAsN()
{
	// Default Constructor

	//GaAsN::LL = 5.6533; // lattice constant parameter
	//GaAsN::delta = 17e-3; // spin-split off energy

	//// define interpolation coefficients
	//C11.set_coeffs(18.7, 8.33, 29.3, 11.9); 
	//C12.set_coeffs(12.5, 4.53, 15.9, 5.38); 
	//alc.set_coeffs(4.98, 6.0583, 4.5, 5.6533);
	//a.set_coeffs(-0.35, -6.08, 3.0, -8.33);
	//b.set_coeffs(-1.2, -1.8, -2.2, -2.0); 
}

GaAsN::GaAsN(double wavelength)
{
	set_wavelength(wavelength); 
}

double GaAsN::refractive_index(double nfrac)
{
	// source: data stripped from various papers
	// valid between 0% and 5% N at a wavelength of 1300 nm
	// see document EAM_Results_21_10_2014_NB.pdf

	try{
		// Model is only valid on certain wavelength range

		if( wavelength > 1.0 && wavelength < 1.5 && nfrac >= 0.0 && nfrac < 0.06){

			return (3.41326 + 1.94069*nfrac); 
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in GaAsN::refractive_index()\n";
			reason += "wavelength = " + template_funcs::toString( wavelength , 3) + " (um)\n"; 
			reason += "Allowed range for GaAsN model is [ 1.1, 1.5 ]\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

double GaAsN::bandgap_energy(double nfrac)
{
	// model for the band-gap energy of the material GaAs_{1-x}N_{x}
    // source: R. Kudrawiec, J. Appl. Phys., 101 (023522), 2007
    // energy is expressed in eV, model valid for x in [0, 5/100]
    // at a temperature of T = 300 K
    // As nfrac -> 0 GaAs_{1-x}N_{x} -> GaAs

	try{
		// Model is only valid on certain infrac range

		if(nfrac >= 0.0 && nfrac < 0.06){

			double E_N_GaAs = 1.65; // contribution of N relative to GaAs
			double C_MN_GaAs = 2.7; // contribution due to bowing caused by N

			double e1 = E_N_GaAs + GaAs::bandgap_energy(); 
			double e2 = std::sqrt( template_funcs::DSQR( E_N_GaAs - GaAs::bandgap_energy() ) + 4.0 * template_funcs::DSQR( C_MN_GaAs )*nfrac ); 

			return 0.5*(e1 - e2);
		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in GaAsN::bandgap_energy()\n";
			reason += "infrac = " + template_funcs::toString( nfrac , 3) + " (um)\n"; 
			reason += "Allowed range for InGaAs model is [ 0, 1 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Definition of class InAsN
InAsN::InAsN()
{
	// Default Constructor
	InAsN::LL = 5.6533; // lattice constant parameter
	InAsN::delta = 17e-3; // spin-split off energy
}

InAsN::InAsN(double wavelength)
{
	set_wavelength(wavelength); 
}

double InAsN::refractive_index(double nfrac)
{
	// no model available for InAsN refractive index
	// use the data for InAs and assume it's a good approximation
	// R. Sheehan 14 - 6 - 2016

	try{
		// Model is only valid on certain wavelength range

		if( wavelength > 1.0 && wavelength < 1.5 && nfrac >= 0.0 && nfrac < 0.06){

			return InAs::refractive_index(); 
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in GaAsN::refractive_index()\n";
			reason += "wavelength = " + template_funcs::toString( wavelength , 3) + " (um)\n"; 
			reason += "Allowed range for GaAsN model is [ 1.1, 1.5 ]\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

double InAsN::bandgap_energy(double nfrac)
{
	// model for the band-gap energy of the material InAs_{1-x}N_{x}
    // source: R. Kudrawiec, J. Appl. Phys., 101 (023522), 2007
    // energy is expressed in eV, model valid for x in [0, 5/100]
    // at a temperature of T = 300 K
    // As nfrac -> 0 InAs_{1-x}N_{x} -> InAs

	try{
		// Model is only valid on certain infrac range

		if(nfrac >= 0.0 && nfrac < 0.06){

			double E_N_InAs = 1.44; // contribution of N relative to GaAs
			double C_MN_InAs = 2.0; // contribution due to bowing caused by N

			double e1 = E_N_InAs + InAs::bandgap_energy(); 
			double e2 = std::sqrt( template_funcs::DSQR( E_N_InAs - InAs::bandgap_energy() ) + 4.0 * template_funcs::DSQR( C_MN_InAs )*nfrac ); // BAC model

			return 0.5*(e1 - e2);
		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in InGaAs::bandgap_energy()\n";
			reason += "infrac = " + template_funcs::toString( nfrac , 3) + " (um)\n"; 
			reason += "Allowed range for InGaAs model is [ 0, 1 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

/****************************************************************************************************************************/
/*                                         Quaternaries                                                                     */
/****************************************************************************************************************************/

// Definition of class InGaAsP
InGaAsP::InGaAsP()
{
	// Default Constructor
}

InGaAsP::InGaAsP(double wavelength)
{
	set_wavelength(wavelength); 
}

double InGaAsP::refractive_index(double asfrac)
{
	//Model for the refactive index of In_{1-x}Ga_{x}As_{y}P_{1-y}
	//Taken from Broberg and Lindgren, J. Appl. Phys., 55(9), 1984
	//The model is valid for 0 <= y < 1 and for wavelengths above the band-gap wavelength l > lg
	//R. Sheehan 24 - 3 - 2010

	try{
		// Model is only valid on certain wavelength range

		double y = asfrac; // As fraction
		double Eg=( 1.35 - 0.72*y + 0.12*template_funcs::DSQR(y) ); // Band gap energy in units of eV
		double lg = (1.24/Eg); // Band gap wavelength

		if( wavelength > lg && asfrac >= 0.0 && asfrac <= 1.0){

			double Egsqr=template_funcs::DSQR(Eg);
			double x=((0.1894*y)/(0.4184-0.013*y)); // Mole fraction x as a function of y
			double En=(1.24/wavelength); // Wavelength in units of energy
			double Ensqr=template_funcs::DSQR(En);
			double Enfourth=template_funcs::DSQR(Ensqr);
			
			double Ed=((12.36*x-12.71)*y+7.54*x+28.91); // Some other parameter
			double Eo=(0.595*template_funcs::DSQR(x)*(1-y)+1.626*x*y-1.891*y+0.524*x+3.391); // Some other parameter
			
			double Eosqr=template_funcs::DSQR(Eo);
			double Eocube=Eo*Eosqr;
			double Et=((PI*Ed)/(2*Eocube*(Eosqr-Egsqr)));

			double v1=(Ed/Eo);
			double v2=((Ed*Ensqr)/Eocube);
			double v3=((Et*Enfourth)/PI);
			double v4=log((2.0*Eosqr-Egsqr-Ensqr)/(Egsqr-Ensqr));

			return sqrt(1.0+v1+v2+v3*v4); // This is the refractive index
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in InGaAsP::refractive_index(double asfrac)\n";
			reason += "wavelength = " + template_funcs::toString(wavelength, 3) + " (um)\n"; 
			reason += "Allowed range for InGaAsP model is l > lg(asfrac)\n"; 
			reason += "Or As fraction is too large: asfrac = " + template_funcs::toString(asfrac, 3) + "\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

double InGaAsP::bandgap_energy(double asfrac)
{
	// model for the band-gap energy of the material In_{x}Ga_{1-x}As_{y}P_{1-y}
    // source: http://www.ioffe.ru/SVA/NSM/Semicond/GaInAsP/index.html, or Chuang
	// several models available depending on what substrate is required
	// for telecomms purposes you generally want a the model to be lattice matched to InP
	// R. Sheehan 14 - 6 - 2016
    
	try{
		// Model is only valid on certain infrac range

		if(asfrac >= 0.0 && asfrac <= 1.0){

			return 0.0;
		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in InGaAs::bandgap_energy()\n";
			reason += "infrac = " + template_funcs::toString( asfrac , 3) + " (um)\n"; 
			reason += "Allowed range for InGaAs model is [ 0, 1 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Definition of class AlInGaAs
AlInGaAs::AlInGaAs()
{
	// Default Constructor
}

AlInGaAs::AlInGaAs(double wavelength)
{
	set_wavelength(wavelength); 
}

double AlInGaAs::refractive_index(double alfrac)
{
	//Model for the refactive index of Al_{1-x}In_{x}Ga_{y}As_{1-y}
	//Taken from "Empirical modeling of the refractive index for (AlGaIn)As lattice matched to InP", Semicond. Sci. Technol. 25 (2010)
	//The model is valid for 0 <= alfrac < 1 
	//R. Sheehan 22 - 9 - 2011

	try{
		// Model is only valid on certain wavelength range

		double a, b, csqr, lsqr, numer, denom; 

		// add in correct polynomial calculation method to Useful
		double A[4] = {10.761, -2.2617, -5.3678, 5.4943};
		double B[4] = {0.9485, -0.2369, 5.2374, -5.1947};
		double C[4] = {1.5986, -3.0266, 2.3326, -0.3629};

		// evaluate the polynomials
		a = A[3]; b = B[3]; csqr = C[3]; 
		for(int j=2; j>=0; j--){
			a = a*alfrac + A[j]; 
			b = b*alfrac + B[j];
			csqr = csqr*alfrac + C[j]; 
		}

		csqr = template_funcs::DSQR(csqr); 

		lsqr = template_funcs::DSQR(wavelength); 

		numer = (b*lsqr); 
		denom = (lsqr - csqr);

		if( fabs(denom) > EPS && alfrac >= 0.0 && alfrac <= 1.0 ){

			return sqrt(a + numer / denom);
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in AlInGaAs::refractive_index(double asfrac)\n";
			reason += "wavelength = " + template_funcs::toString(wavelength, 3) + " (um)\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

double AlInGaAs::bandgap_energy(double asfrac)
{
	// model for the band-gap energy of the material AlInGaAs
    // source chuang
    
	try{
		// Model is only valid on certain infrac range

		if(asfrac >= 0.0 && asfrac <= 1.0){

			return 0.0;
		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in InGaAs::bandgap_energy()\n";
			reason += "infrac = " + template_funcs::toString( asfrac , 3) + " (um)\n"; 
			reason += "Allowed range for InGaAs model is [ 0, 1 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Definition of class InGaAsN
InGaAsN::InGaAsN()
{
	// Default Constructor
	bowing = 0.477; // parameter associated with interpolation of InAsN and GaAsN
}

InGaAsN::InGaAsN(double wavelength)
{
	GaAsN::set_wavelength(wavelength); 
}

double InGaAsN::refractive_index(double infrac, double nfrac)
{
	// Model is from an unknown source

	try{

		if( infrac >= 0.0 && infrac <= 1.0 ){

			return 0.0;
		}
		else{
			std::string reason; 
			reason = "Attempting to compute RI outside of allowed range in AlInGaAs::refractive_index(double asfrac)\n";
			reason += "wavelength = " + template_funcs::toString(GaAsN::wavelength, 3) + " (um)\n"; 
			throw std::range_error(reason); 
		}

	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}

}

double InGaAsN::bandgap_energy(double infrac, double nfrac)
{
	// model for the band-gap energy of the material In_{y}Ga_{1-y}As_{1-x}N_{x}
    // source: R. Kudrawiec, J. Appl. Phys., 101 (023522), 2007
    // energy is expressed in eV, model valid for x in [0, 5/100]
    // at a temperature of T = 300 K	
    
	try{
		// Model is only valid on certain infrac range

		if(infrac >= 0.0 && infrac <= 1.0 && nfrac >= 0.0 && nfrac < 0.06){

			double t1 = (1.0 - infrac)*GaAsN::bandgap_energy(nfrac); 
			double t2 = infrac * InAsN::bandgap_energy(nfrac); 
			double t3 = infrac * (1.0 - infrac) * bowing; 
			double t4 = GaAsN::strain_correction(nfrac, infrac); 

			return ( t1 + t2 - t3 + t4 ); 
		}
		else{
			std::string reason; 
			reason = "Attempting to compute Eg outside of allowed range in InGaAs::bandgap_energy()\n";
			reason += "infrac = " + template_funcs::toString( infrac , 3) + " (um)\n"; 
			reason += "Allowed range for InGaAs model is [ 0, 1 ]\n"; 
			throw std::range_error(reason); 
		}
	}
	catch(std::range_error &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

//double InGaAsN::bandgap_energy_with_strain(double infrac, double nfrac)
//{
//	// model for the band-gap energy of the material InAs_{1-x}N_{x}
//    // source: R. Kudrawiec, J. Appl. Phys., 101 (023522), 2007
//    // energy is expressed in eV, model valid for x in [0, 5/100]
//    // at a temperature of T = 300 K
//	// This method includes strain compensation
//    
//	try{
//		// Model is only valid on certain infrac range
//
//		if(infrac >= 0.0 && infrac <= 1.0 && nfrac >= 0.0 && nfrac < 0.06){
//
//			double t1 = InGaAsN::bandgap_energy(infrac, nfrac); 
//			double t2 = GaAsN::strain_correction(nfrac, infrac); 
//
//			return (t1 + t2);
//		}
//		else{
//			std::string reason; 
//			reason = "Attempting to compute Eg outside of allowed range in InGaAs::bandgap_energy()\n";
//			reason += "infrac = " + template_funcs::toString( infrac , 3) + " (um)\n"; 
//			reason += "Allowed range for InGaAs model is [ 0, 1 ]\n"; 
//			throw std::range_error(reason); 
//		}
//	}
//	catch(std::range_error &e){
//		useful_funcs::exit_failure_output(e.what()); 
//		exit(EXIT_FAILURE); 
//	}
//}