#ifndef ATTACH_H
#include "Attach.h"
#endif

fresnel::fresnel()
{
	// Default constructor

	theta_in = theta_t = wavelength = n1 = n2 = nrat = wavelength; 
}

fresnel::fresnel(double lambda, material *m1, material *m2)
{

	set_params(lambda, m1, m2); 
}

void fresnel::set_params(double lambda, material *m1, material *m2)
{
	// assign the parameters for the class
	// wavelength should be in units of um

	try {
		if (lambda > 0.0) {
			wavelength = lambda;
			mat1 = m1;
			mat2 = m2;
			mat1->set_wavelength(wavelength);
			mat1->set_wavelength(wavelength);
			n1 = mat1->refractive_index();
			n2 = mat2->refractive_index();
			if (n2 > 0.0) {
				nrat = n1 / n2;
			}
			else {
				std::string reason;
				reason = "Error: fresnel::set_params(double lambda, material *m1, material *m2)\n";
				reason += "n2: " + template_funcs::toString(n2) + " is not positive\n";
				throw std::runtime_error(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: fresnel::set_params(double lambda, material *m1, material *m2)\n";
			reason += "lambda: " + template_funcs::toString(lambda) + " is not positive\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what(); 
	}
}

double fresnel::reflectance(double angle, bool polarisation)
{
	// compute the power reflection for a given polarisation and angle of incidence

	try {
		bool c1 = angle >= 0.0 && angle <= PI_2 ? true : false;
		bool c2 = wavelength > 0.0 ? true : false;
		bool c3 = nrat > 0.0 ? true : false;
		bool c9 = c1 && c2 && c3 ? true : false; 

		if (c9) {
			double cangle = cos(angle); 
			double sangle = sin(angle); 
			double t1 = polarisation == TE ? n1 * cangle : n2*cangle; 
			double t2 = sqrt( 1.0 - template_funcs::DSQR(nrat * sangle) );
			double t3 = polarisation == TE ? n2 * t1 : n1 * t1; 
			double numer = polarisation == TE ? t1 - t3 : t3 - t1; 
			double denom = polarisation == TM ? t1 + t3 : t3 + t1; 
			if (fabs(denom) > 0.0) {
				return numer / denom; 
			}
			else {
				return 0.0; 
				std::string reason = "Error: double fresnel::reflectance(double angle, bool polarisation)\n"; 
				reason += "Attempt to divide by zero\n"; 
				throw std::runtime_error(reason); 
			}
		}
		else {
			return 0.0; 
			std::string reason;
			reason = "Error: double fresnel::reflectance(double angle, bool polarisation)\n";
			if (!c1) reason += "Angle of incidence not defined correctly\n";
			if (!c2) reason += "Operating wavelength is not defined\n";
			if (!c3) reason += "RI ratio is not defined\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what(); 
	}
}