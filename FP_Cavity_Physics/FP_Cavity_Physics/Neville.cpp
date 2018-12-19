#ifndef ATTACH_H
#include "Attach.h"
#endif

void interpolation::polint(std::vector<double> &xa, std::vector<double> &ya, double x, double &y, double &dy)
{
	// Given vectors xa, ya compute an approximation to y = f(x) with error dy at the point x 
	// Scheme is based on Neville's method. 
	// Technically, Neville's algorithm will extrapolate to points beyond the data set, however, extrapolation is not allowed here

	try {
		bool c1 = xa.size() == ya.size() ? true : false; 
		bool c2 = x >= xa.front() ? true : false; 
		bool c3 = x <= xa.back() ? true : false; 
		bool c4 = !xa.empty() ? true : false; 
		bool c5 = !ya.empty() ? true : false;
		bool c9 = c1 && c2 && c3 && c4 && c5 ? true : false;

		if (c9) {
			int i, m, ns = 1, n = static_cast<int>(xa.size());
			double den, dif, dift, ho, hp, w;
			double dyold = 10, yold; // Code to prevent underflow
			std::vector<double> c(n, 0.0);
			std::vector<double> d(n, 0.0);

			dif = fabs(x - xa[0]);

			for (i = 0; i < n; i++) {
				if ((dift = fabs(x - xa[i]))<dif) {
					ns = i;
					dif = dift;
				}
				c[i] = ya[i];
				d[i] = ya[i];
			}
			y = ya[ns--];
			for (m = 1; m < n; m++) {
				for (i = 0; i < n - m; i++) {
					ho = xa[i] - x;
					hp = xa[i + m] - x;
					w = c[i + 1] - d[i];
					if ((den = ho - hp) == 0.0) {
						std::string reason; 
						reason = "Error: void interpolation::polint(std::vector<double> &xa, std::vector<double> &ya, double x, double &y, double &dy)\n"; 
						reason += "den == 0.0"; 
						throw std::runtime_error(reason); 
					}
					den = w / den;
					d[i] = hp * den;
					c[i] = ho * den;
				}
				y += (dy = (2 * (ns+1)<(n - m) ? c[ns + 1] : d[ns--]));
				//Code to prevent underflow
				if (fabs(dy)<1.0e-8) {
					break;
				}
				else if (fabs(dy)>fabs(dyold)) {
					y = yold;
					break;
				}
				else {
					yold = y;
					dyold = dy;
				}
			}
			c.clear(); d.clear(); 
		}
		else {
			std::string reason;
			reason = "Error: void interpolation::polint(std::vector<double> &xa, std::vector<double> &ya, double x, double &y, double &dy)\n";
			if (!c1) reason += "input arrays are not the same dimension\n"; 
			if (!c2 || !c3) reason += "Position: " + template_funcs::toString(x, 3) + " is not in the data range\n"; 
			if (!c4 || !c5) reason += "Input arrays are empty\n"; 
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