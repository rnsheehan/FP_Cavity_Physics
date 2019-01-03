#ifndef ATTACH_H
#include "Attach.h"
#endif

fresnel::fresnel()
{
	// Default constructor

	nrat_sqr = theta_critical = theta_brewster = theta_t = wavelength = n1 = n2 = nrat = wavelength = 0.0;
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
			mat1 = m1;
			mat2 = m2;
			update_wavelength(lambda); // recompute the refractive index values each time the wavelength changes
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
}

void fresnel::update_wavelength(double lambda)
{
	// recompute the refractive index values each time the wavelength changes without changing the materials

	try {
		if (lambda > 0.0) {
			wavelength = lambda;
			mat1->set_wavelength(wavelength);
			mat1->set_wavelength(wavelength);
			// It is assumed that the light is propagating from material with RI = n1 into material with RI = n2
			n1 = mat1->refractive_index();
			n2 = mat2->refractive_index();
			if (n1 > 0.0) {
				nrat = n2 / n1; // nrat > 1 => external reflection, nrat < 1 => internal reflection
				nrat_sqr = template_funcs::DSQR(nrat);
				theta_critical = nrat < 1.0 ? asin(nrat) : PI_2;
				theta_brewster = atan(nrat);
			}
			else {
				std::string reason;
				reason = "Error: void fresnel::update_wavelength(double lambda)\n";
				reason += "n1: " + template_funcs::toString(n1) + " is not positive\n";
				throw std::invalid_argument(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: void fresnel::update_wavelength(double lambda)\n";
			reason += "lambda: " + template_funcs::toString(lambda) + " is not positive\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double fresnel::transmission_angle(double angle)
{
	// compute the transmission angle after reflection
	try {
		bool c1 = angle >= 0.0 ? true : false;
		bool c3 = nrat > 0.0 ? true : false;
		bool c4 = angle < theta_critical ? true : false;
		bool c9 = c1 && c3 && c4 ? true : false;

		if (c9) {
			return asin(sin(angle)/nrat); 
		}
		else {
			if (!c4) {
				return 0.0; 
			}
			else {
				return 0.0; 
				std::string reason;
				reason = "Error: double fresnel::transmission_angle(double angle)\n";
				if (!c1) reason += "Angle of incidence not defined correctly\n";
				if (!c3) reason += "RI ratio is not defined\n";
				throw std::invalid_argument(reason);
			}
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double fresnel::reflectance(double angle, bool polarisation)
{
	// compute the power reflection for a given polarisation and angle of incidence
	// angle should be input in units of radians
	// See Fowles, section 2.7, eqns 2.58 and 2.59, page 44

	try {
		bool c1 = angle >= 0.0 ? true : false;
		bool c2 = wavelength > 0.0 ? true : false;
		bool c3 = nrat > 0.0 ? true : false;
		bool c4 = angle < theta_critical ? true : false; 
		bool c9 = c1 && c2 && c3 && c4 ? true : false; 

		if (c9) {
			double cangle = cos(angle); 
			double sangle = sin(angle); 
			double t1 = polarisation == TE ? cangle : nrat_sqr*cangle;
			double t2 = sqrt( nrat_sqr - template_funcs::DSQR(sangle) );
			double numer = polarisation == TE ? t1 - t2 : t2 - t1; 
			double denom = t1 + t2; 
			if (fabs(denom) > 0.0) {
				return template_funcs::DSQR(numer / denom); 
			}
			else {
				return 0.0; // nrat > 1 => external reflection, nrat < 1 => internal reflection 
				std::string reason = "Error: double fresnel::reflectance(double angle, bool polarisation)\n"; 
				reason += "Attempt to divide by zero\n"; 
				throw std::runtime_error(reason); 
			}
		}
		else {
			if (!c4) {
				// angle > critical angle => no transmission
				return 1.0; // nrat > 1 => external reflection, nrat < 1 => internal reflection 
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
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what(); 
	}
}

double fresnel::transmittance(double angle, bool polarisation)
{
	// compute the power reflection for a given polarisation and angle of incidence
	// T = 1 - R by conservation of energy

	return 1.0 - reflectance(angle, polarisation); 
}