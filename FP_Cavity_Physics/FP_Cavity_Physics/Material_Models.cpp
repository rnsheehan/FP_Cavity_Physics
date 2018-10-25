#ifndef ATTACH_H
#include "Attach.h"
#endif

material::material()
{
	// Default constructor

	data_in_memory = false;

	n_data_points = n_cols = 0;

	temperature = 298.15; // Assume that temperature is 25 ( C ) unless otherwise stated

	energy = wavelength = 0.0;
}

void material::convert_energy_to_wavelength()
{
	// convert a bandgap energy expressed in eV to a wavelength expressed in nm
	// R. Sheehan 14 - 6 - 2016

	try {
		if (energy > 0.0) {
			wavelength = 1240.0 / energy;
		}
		else {
			std::string reason = "Value stored in energy is not correct\n";
			reason += "energy = " + template_funcs::toString(energy, 2) + "\n";
			throw std::range_error(reason);
		}

	}
	catch (std::range_error &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void material::convert_wavelength_to_energy()
{
	// convert a wavelength expressed in nm to a bandgap energy expressed in eV
	// R. Sheehan 14 - 6 - 2016

	try {
		if (wavelength > 0.0) {
			energy = 1240.0 / wavelength;
		}
		else {
			std::string reason = "Value stored in energy is not correct\n";
			reason += "energy = " + template_funcs::toString(energy, 2) + "\n";
			throw std::range_error(reason);
		}

	}
	catch (std::range_error &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void material::convert_C_K(double &Cvalue)
{
	// convert temperature in Celcius to Kelvin scale
	// R. Sheehan 14 - 6 - 2016

	try {
		if (Cvalue > 0.0) {
			temperature = Cvalue + 273.15; // temperature is expressed in units of K
		}
		else {
			std::string reason = "Value stored in Cvalue is not correct\n";
			reason += "Cvalue = " + template_funcs::toString(Cvalue, 2) + "\n";
			throw std::range_error(reason);
		}

	}
	catch (std::range_error &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void material::convert_K_C(double &Kvalue)
{
	// convert temperature in Kelvin to Celcius scale
	// R. Sheehan 14 - 6 - 2016

	try {
		if (Kvalue > 0.0) {
			temperature = Kvalue - 273.15; // temperature is expressed in units of C
		}
		else {
			std::string reason = "Value stored in Cvalue is not correct\n";
			reason += "Cvalue = " + template_funcs::toString(Kvalue, 2) + "\n";
			throw std::range_error(reason);
		}

	}
	catch (std::range_error &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}