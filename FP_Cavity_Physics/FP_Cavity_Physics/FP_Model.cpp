#ifndef ATTACH_H
#include "Attach.h"
#endif

fp_cavity::fp_cavity()
{
	// Default Constructor
	nwl = 0; 
	params_defined = false; 
	wl1 = wl2 = dwl = wlc = n1 = n2 = length = alpha = R = r = rpr = phase = 0.0;
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
			wlc = 0.5*(wl1 + wl2); // central wavelength is the average of the spectral endpoints
			nwl = 201; // specify num. subdivisions
			dwl = (wl2 - wl1) / ((double)(nwl - 1)); // specify wl spacing 
			mat1 = m1; 
			mat2 = m2;
			params_defined = true; 
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

void fp_cavity::dispersion(bool ignore_dispersion, double wavelength)
{
	// ignore dispersion == true means that dispersion is ignored when computing the material RI
	// wavelength is only accounted for when the phase term is computed

	try {
		bool c1 = useful_funcs::test_gleq(GEQ, wavelength, wl1) && useful_funcs::test_gleq(LEQ, wavelength, wl2) ? true : false;

		if (c1) {
			if (ignore_dispersion) {
				// update the material wavelength at the central wavelength only
				mat1->set_wavelength(wlc); mat2->set_wavelength(wlc);
			}
			else {
				// update the material wavelength
				mat1->set_wavelength(wavelength); mat2->set_wavelength(wavelength);
			}		

			// update the local RI values
			n1 = mat1->refractive_index(); n2 = mat2->refractive_index(); 

			// the phase value is actually the argument of the sin function, however, taking the sin^{2} value for the sake of computational efficiency	
			phase = template_funcs::DSQR( sin( (Two_PI * n2 * length) / wavelength ) );

			// Compute the wavelength dependent mirror reflectivity
			fresnel calc; 

			calc.set_params(wavelength, n1, n2); 

			r = calc.ref_coeff(angle, TE);

			calc.set_params(wavelength, n2, n1);

			rpr = calc.ref_coeff(angle, TE);

			R = sqrt(template_funcs::DSQR(r*rpr)); // Cavity Reflectance

			//std::cout << "n1: " << n1 << " , n2: " << n2 << "\n"; 
			//std::cout << wavelength << " , " << r << " , " << rpr << " , " << R << "\n"; 
		}
		else {
			std::string reason;
			reason = "Error: double fp_cavity::dispersion(double wavelength)\n";
			reason += "Attempting to compute dispersion beyond allowed range\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double fp_cavity::finesse(double wavelength)
{
	// compute the FP cavity finesse

	try {
		bool c1 = wavelength >= wl1 && wavelength <= wl2 ? true : false;
		bool c2 = R > 0.0 && R < 1 ? true : false; 

		if (c1 && c2) {
			return  (4.0 * R) / (template_funcs::DSQR( 1.0 - R )); 
		}
		else {
			return 0.0; 
			std::string reason;
			reason = "Error: double fp_cavity::finesse(double wavelength)\n";
			if(!c1) reason += "Attempting to compute finesse beyond allowed range\n";
			if (!c2) reason += "R is not yet defined\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double fp_cavity::airy(double F)
{
	// compute the FP cavity Airy function
	// F is the cavity finesse
	try {
		bool c1 = F > 0.0 ? true : false;

		if (c1) {
			double denom = 1.0 + F * phase; 
			return fabs(denom) > 0.0 ? 1.0 / denom : 0.0; 
		}
		else {
			return 0.0; 
			std::string reason;
			reason = "Error: double fp_cavity::airy(double wavelength)\n";
			reason += "Finesse F not correctly defined\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void fp_cavity::reflection(double wavelength, double &F, double &A, double &RFP, bool ignore_dispersion, bool loud)
{
	// compute the FP cavity reflection
	try {
		bool c1 = useful_funcs::test_gleq(GEQ, wavelength, wl1) && useful_funcs::test_gleq(LEQ, wavelength, wl2) ? true : false;

		if (c1) {
			dispersion(ignore_dispersion, wavelength); // define all parameters needed for the input wavelength
			F = finesse(wavelength); // compute the finesse
			A = airy(F); // compute the value of the Airy function
			RFP = F * A * phase; 

			if (loud) {
				std::cout << "wavelength: " << wavelength << "\n";
				std::cout << "n1: " << n1 << " , n2: " << n2 << " , phase: " << phase << "\n";
				std::cout << "r: " << r << " , r': " << rpr << " , R: " << R << "\n";
				std::cout << "F: " << F << " , A(F): " << A << "\n"; 
				std::cout << "Cavity Reflectance: " << RFP << "\n\n"; 
			}
		}
		else {
			std::string reason;
			reason = "Error: void fp_cavity::reflection(double wavelength, double &F, double &A, double &RFP, bool loud)\n";
			reason += "Attempting to compute reflection beyond allowed range\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void fp_cavity::compute_spectrum(bool ignore_dispersion, bool loud)
{
	// compute the FP cavity reflection spectrum
	// ignore dispersion == true means that the refractive index variation due to wavelength is not included in the calculation

	try {
		if (params_defined) {
			std::vector<double> lambda(nwl, 0.0); 
			std::vector<double> refl(nwl, 0.0);
			std::vector<double> airyvals(nwl, 0.0);

			double wl = wl1;
			double airy, fine, rfp;
			for (int i = 0; i < nwl; i++) {
				
				reflection(wl, fine, airy, rfp); 
				
				lambda[i] = wl; airyvals[i] = airy; refl[i] = rfp;
				
				if(loud && i%5 == 0) std::cout << wl << " , " << fine<< " , " << airy << " , " << rfp << "\n"; 
				
				wl += dwl; 
			}

			std::string data_file; 

			data_file = "wavelength.txt"; 
			vecut::write_into_file(data_file, lambda);

			data_file = "fp_airy.txt";
			vecut::write_into_file(data_file, airyvals);

			data_file = "fp_reflectivity.txt";
			vecut::write_into_file(data_file, refl);
		}
		else {
			std::string reason;
			reason = "Error: void fp_cavity::compute_spectrum()\n";
			reason += "Device parameters not defined\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}