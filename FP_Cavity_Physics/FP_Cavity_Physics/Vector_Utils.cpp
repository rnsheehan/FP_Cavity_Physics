#ifndef ATTACH_H
#include "Attach.h"
#endif

std::vector<double> vecut::get_row(std::vector<std::vector<double>> &data, int r_num)
{
	// extract a row of data from an existing 2D array

	try {
		bool c1 = !data.empty() ? true : false; 
		bool c2 = r_num > -1 && r_num < (int)(data.size()) ? true : false; 
		bool c3 = c1 && c2 ? true : false; 

		if (c3) {
			return data[r_num]; 
		}
		else {
			return std::vector<double>(); 
			std::string reason; 
			reason = "Error: std::vector<double> get_row(std::vector<std::vector<double>> &data, int r_num)\n"; 
			if (!c1) reason += "input array not defined\n"; 
			if (!c2) reason += "attempting to access beyond bounds of array\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

std::vector<double> vecut::get_col(std::vector<std::vector<double>> &data, int c_num)
{
	// extract a column of data from an existing 2D array

	try {
		bool c1 = !data.empty() ? true : false;
		bool c2 = c_num > -1 && c_num < (int)(data[0].size()) ? true : false;
		bool c3 = c1 && c2 ? true : false;

		if (c3) {
			//int n_rows = (int)(data.size());
			std::vector<double> data_col(data.size());

			for (size_t i = 0; i < data.size(); i++) {
				data_col[i] = data[i][c_num]; 
			}

			return data_col; 
		}
		else {
			return std::vector<double>();
			std::string reason;
			reason = "Error: std::vector<double> get_col(std::vector<std::vector<double>> &data, int c_num)\n";
			if (!c1) reason += "input array not defined\n";
			if (!c2) reason += "attempting to access beyond bounds of array\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void vecut::read_into_vector(std::string &filename, std::vector<double> &data, int &n_pts, bool loud)
{
	// read a single column of data from a file into a vector
	// It is assumed that data is empty when the values are being read in
	// R. Sheehan 11 - 9 - 2017

	try {
		std::ifstream the_file;
		the_file.open(filename, std::ios_base::in);

		if (the_file.is_open()) {

			if (loud) std::cout << filename << " opened for reading\n";

			if (!data.empty()) data.clear(); // make sure that data is empty when values are being read in

			double value;
			n_pts = 0;
			while (the_file >> value) {
				data.push_back(value);
				n_pts++;
			}

			if (loud) std::cout << template_funcs::toString(n_pts) << " data were read from " << filename << "\n";

			the_file.close();
		}
		else {
			std::string reason;
			reason = "Error: void vecut::read_into_vector(std::string &filename, std::vector<double> &data, int &n_pts, bool loud)\n";
			reason += "Cannot open: " + filename + "\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void vecut::read_into_matrix(std::string &filename, std::vector<std::vector<double>> &data, int &n_rows, int &n_cols, bool loud)
{
	// read an array of data from a file
	// store the data in a matrix of size n_rows * n_cols
	// R. Sheehan 18 - 12 - 2018

	try {
		std::ifstream the_file;
		the_file.open(filename, std::ios_base::in);

		if (the_file.is_open()) {
			
			if (loud) std::cout << filename << " opened for reading\n";

			if (!data.empty()) data.clear(); // make sure that data is empty when values are being read in

			bool has_header = false; 
			std::string line, item;

			char endline = '\n';			
			char tab_token = '\t'; 
			char comma_token = ','; 
			char sep_token; 

			// Initialise the values to zero
			n_rows = 0; n_cols = 0;

			// Determine which token separates data in the file
			the_file.seekg(1, std::ios::beg); // move to the start of the file
			std::getline(the_file, line, endline); 
			if (line.find(tab_token) != std::string::npos) sep_token = tab_token; 
			if (line.find(comma_token) != std::string::npos) sep_token = comma_token;

			if (loud) std::cout << filename << " uses " << sep_token << " as a token\n";

			// Determine if first line is a file header
			if (isalpha(line[2])) has_header = true;

			// Count the number of rows and columns
			// This only seems to work when the data are separated by ',' also works for '\t' and ' '
			// http://www.cplusplus.com/reference/string/string/getline/
			// getline (istream& is, string& str, char delim)
			// Extracts characters from is and stores them into str until the delimitation character delim is found
			the_file.seekg( has_header ? 1 : 0, std::ios::beg); // move to the start of the file
			while (std::getline(the_file, line, endline)) {
				n_rows++;
				std::istringstream linestream(line);
				if (n_rows == 1) {
					while (std::getline(linestream, item, sep_token)) {
						n_cols++;
					}
				}
			}

			if (loud) std::cout << filename << " contains " << n_rows << " rows and " << n_cols << " columns\n"; 

			if (n_rows > 1 && n_cols > 0) {
				// Allocate the memory required to store the data
				data.resize(n_rows);
				for (size_t k = 0; k < data.size(); k++) {
					data[k].resize(n_cols, 0.0); 
				}

				the_file.clear(); // empty a buffer, needed to ensure data can be read from the file
				the_file.seekg(has_header ? 1 : 0, std::ios::beg); // move to the start of the file

				int i, j;
				
				i = 0;
				while (std::getline(the_file, line, endline)) {
					std::istringstream linestream(line);
					j = 0;
					while (std::getline(linestream, item, sep_token)) {
						data[i][j] = atof(item.c_str());
						j++;
					}
					i++;
				}

				the_file.close(); 
			}
			else {
				std::string reason;
				reason = "Error: void read_into_matrix(std::string &filename, std::vector<std::vector<double>> &data, int &n_rows, int &n_cols, bool loud = false)\n";
				reason = filename + " contains no data\n"; 
				reason += "n_rows: " + template_funcs::toString(n_rows) + ", n_cols: " + template_funcs::toString(n_cols) + "\n"; 
				throw std::runtime_error(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: void read_into_matrix(std::string &filename, std::vector<std::vector<double>> &data, int &n_rows, int &n_cols, bool loud = false)\n";
			reason += "Cannot open: " + filename + "\n";
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

void vecut::write_into_file(std::string &filename, std::vector<double> &data, bool loud)
{
	// write a single column of data to a new file

	try {

		if ( !data.empty() && filename != empty_str && useful_funcs::valid_filename_length(filename) ) {

			std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);

			if (write.is_open()) {

				for (size_t k = 0; k < data.size(); k++) {
					write << std::setprecision(10) << data[k] << "\n"; 
				}
			
				write.close(); 
			}
			else {
				std::string reason; 
				reason = "Error: void vecut::write_into_file(std::string &filename, std::vector<double> &data, bool loud)\n";
				reason += "Could not open file: " + filename + "\n"; 
				throw std::runtime_error(reason); 
			}
		}
		else {
			std::string reason;
			reason = "Error: void vecut::write_into_file(std::string &filename, std::vector<double> &data, bool loud)\n";
			reason += "Filename: " + filename + " is not valid\n";
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