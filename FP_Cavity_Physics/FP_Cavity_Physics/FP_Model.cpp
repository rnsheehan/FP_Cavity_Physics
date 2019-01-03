#ifndef ATTACH_H
#include "Attach.h"
#endif

fp_cavity::fp_cavity()
{
	// Default Constructor
	nwl = 0; 
	wl1 = wl2 = dwl = n1 = n2 = length = alpha = R = r = rpr = phase = 0.0;
}

fp_cavity::fp_cavity(double incident_angle, double cav_length, double loss_fac, double wl_start, double wl_end, material *m1, material *m2)
{
	set_params(incident_angle, cav_length, loss_fac, wl_start, wl_end, m1, m2);
}

void fp_cavity::set_params(double incident_angle, double cav_length, double loss_fac, double wl_start, double wl_end, material *m1, material *m2)
{
	try {		
		bool c1 = wl_start > 0.0 ? true : false; 
		bool c2 = wl_end > wl_start ? true : false; 
		bool c3 = incident_angle >= 0.0 ? true : false; 
		bool c4 = cav_length > 0.0 ? true : false; 
		bool c9 = c1 && c2 && c3 && c4; 

		if (c9) {
			angle = incident_angle; 
			length = cav_length; 
			alpha = loss_fac; 
			wl1 = wl_start; wl2 = wl_end; // define endpoints of FP spectrum
			nwl = 301; // specify num. subdivisions
			dwl = (wl2 - wl1) / ((double)(nwl - 1)); // specify wl spacing 
			mat1 = m1; 
			mat2 = m2; 
		}
		else {
			std::string reason;
			reason = "Error: void fp_cavity::set_params(double incident_angle, double cav_length, double loss_fac, double wl_start, double wl_end, material *m1, material *m2)\n";
			if (!c2) reason += "Wavelength not correctly defined\n";
			if (!c3) reason += "Incident angle not correctly defined\n";
			if (!c4) reason += "Cavity length not correctly defined\n";
			throw std::invalid_argument(reason);
		}		
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void fp_cavity::dispersion(double wavelength)
{
	try {
		bool c1 = wavelength >= wl1 && wavelength <= wl2 ? true : false; 

		if (c1) {
			// update the material wavelength
			mat1->set_wavelength(wavelength); mat2->set_wavelength(wavelength); 

			// update the local RI values
			n1 = mat1->refractive_index(); n2 = mat2->refractive_index(); 

			phase = (Two_PI * n2 * length) / wavelength; 
			
			// Compute the wavelength dependent mirror reflectivity
			fresnel calc; 

			calc.set_params(wavelength, n1, n2); 

			r = calc.ref_coeff(angle, TE);

			calc.set_params(wavelength, n2, n1);

			rpr = calc.ref_coeff(angle, TE);

			R = r*rpr; // Cavity Reflectance

			std::cout << "n1: " << n1 << " , n2: " << n2 << "\n"; 

			std::cout << wavelength << " , " << r << " , " << rpr << " , " << R << "\n"; 
		}
		else {
			std::string reason;
			reason = "Error: double fp_cavity::dispersion(double wavelength)\n";
			reason += "Attempting to dispersion beyond allowed range\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}