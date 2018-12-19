#ifndef ATTACH_H
#include "Attach.h"
#endif

// Source definitions for the material model object
// R. Sheehan 27 - 10 - 2015

ri_model::ri_model()
{
	// Default constructor

	data_in_memory = false; 
	
	n_data_points = n_cols = material_in_memory = 0; 

	fill_material_list(); 
}

ri_model::~ri_model()
{
	// Deconstructor
	if(data_in_memory){

		//cout<<"\nClearing data from memory\n"; 
		
		ri_data.clear(); 

		data_in_memory = false; 

		n_data_points = n_cols = material_in_memory = 0; 
	}

	materials.clear(); 
}

double ri_model::ri_value(int mat_type, double wavelength, double xfr, double yfr)
{
	// interface method for the ri_model object

	try{
		// Check if the material model is available
		bool found = false; 
		found = ( std::find( materials.begin(), materials.end(), mat_type) != materials.end() ); 

		if(found){

			string the_file = null_string; 

			// Switch the output to the material model you want
			switch(mat_type){
				case RI_GaAs:
					return GaAs(wavelength); 
					break; 
				case RI_AlAs:
					return AlAs(wavelength); 
					break; 
				case RI_SiO2:
					return SiO2(wavelength); 
					break; 
				case RI_SiGeO2:
					return SiGeO2(wavelength); 
					break;
				case RI_AlGaAs:
					return AlGaAs(xfr, wavelength); 
					break; 
				case RI_InGaAs:
					return InGaAs(xfr, wavelength); 
					break;
				case RI_GaAsN1300:
					return GaAsN_1300(xfr); 
					break;
				case RI_InGaAsP:
					return InGaAsP(xfr, wavelength); 
					break; 
				case RI_InGaAsN:
					return InGaAsN(xfr, yfr, wavelength); 
					break;
				case RI_AlInGaAs:
					return AlInGaAs(xfr, wavelength); 
					break;
				case RI_InP:
					the_file = "RI_Data\\InP.txt"; 
					return data_based_ri(RI_InP, the_file, wavelength); 
					break; 
				case RI_GaP:
					the_file = "RI_Data\\GaP.txt"; 
					return data_based_ri(RI_GaP, the_file, wavelength); 
					break; 
				case RI_AlN:
					the_file = "RI_Data\\AlN.txt"; 
					return data_based_ri(RI_AlN, the_file, wavelength); 
					break; 
				case RI_Si:
					the_file = "RI_Data\\Si.txt"; 
					return data_based_ri(RI_Si, the_file, wavelength); 
					break; 
				default:
					cout<<"You should not be able to see this statement\n"; 
					return 0; 
			}
			
		}
		else{
			throw assignment_error(); 
		}
	}
	catch(assignment_error){
		string reason;
		reason = "Material model / data not available in ri_model::ri_value(int mat_type, double wavelength, double xfr, double yfr)\n"; 
		reason += "Material choice = " + toString(mat_type) + "\n"; 
		exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

// ri based on a formula
double ri_model::GaAs(double wavelength)
{
	// source: http://www.batop.com/information/n_GaAs.html
    // type: fit of experimental data
    // model valid for at least lambda (um) in [0.65, 1.8]

	try{
		// Model is only valid on certain wavelength range

		if(wavelength > 0.64 && wavelength < 1.9){

			double A = 8.95; 
			double B = 2.054;
			double Csqr = 0.39; 

			double t1 = ( 1.0 - ( Csqr / DSQR(wavelength) ) ); 
			double t2 = ( A + (B/t1) ); 

			return sqrt(t2); 
		}
		else{
			throw assignment_error(); 
		}

	}
	catch(assignment_error){
		string reason; 
		reason = "Attempting to compute RI outside of allowed range in ri_model::GaAs(double wavelength)\n";
		reason += "wavelength = " + toString(wavelength, 3) + " (um)\n"; 
		reason += "Allowed range for GaAs model is [ 0.65, 1.8 ]\n"; 
		exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}

}

double ri_model::AlAs(double wavelength)
{
	// source: http://www.batop.com/information/n_AlAs.html
    // type: theoretical fit
    // model valid for at least lambda (um) in [0.45, 2.0]

	try{
		// Model is only valid on certain wavelength range

		if(wavelength > 0.44 && wavelength < 2.1){
			// PE = HC/wavelength # photon energy in Joules
			double A0 = 25.3; // constant, determined by fitting with experimental data
			double B0 = -0.8; // constant, determined by fitting with experimental data
			double E0 = 2.95; // fundamental band gap at G-point expressed in eV
			double S0 = 3.25; // fundamental band gap at G-point plus D0 spin-orbit splitting energy expressed in eV
			//double ER = math.pow((E0/S0),1.5) // (E0 / E0+D0)^{3/2}
			double ER = pow((E0/S0),1.5); 

			double chi = HC/(wavelength*E0); 
			double chisqr = DSQR(chi); 

			double t1 = sqrt(1.0+chi); 
			double t2 = sqrt(1.0-chi);
			double fchi = (2.0-t1-t2)/chisqr;

			double chiS0 = HC/(wavelength*S0);
			double chisqrS0 = DSQR(chiS0);

			double t1S0 = sqrt(1.0+chiS0);
			double t2S0 = sqrt(1.0-chiS0);
			double fchiS0 = (2.0-t1S0-t2S0)/chisqrS0;

			double t3 = ( fchi + 0.5*fchiS0*ER );

			return sqrt( A0*t3 + B0 ); 			
		}
		else{
			throw assignment_error(); 
		}

	}
	catch(assignment_error){
		string reason; 
		reason = "Attempting to compute RI outside of allowed range in ri_model::AlAs(double wavelength)\n";
		reason += "wavelength = " + toString(wavelength, 3) + " (um)\n"; 
		reason += "Allowed range for AlAs model is [ 0.45, 2.0 ]\n"; 
		exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

double ri_model::SiO2(double wavelength)
{
	// Refractive Index of Pure Silica
	// Model is based on a Sellmeier polynomial with experimentally derived coefficients
	// Taken from K. Okamoto, "Fundamanetals of Optical Waveguides", 2006
	// wavelength is input in units of microns
	// R. Sheehan 12 - 4 - 2012

	try{

		if(wavelength > 0.4 && wavelength < 2.5){
	
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
			double lsqr = DSQR(wavelength);
	
			for(int i=1; i<=ncoeffs; i++){
				nsqr +=  ( a[i]*lsqr ) / ( lsqr - b[i] ) ;
			}
	
			return sqrt(nsqr);

		}
		else{
			throw assignment_error(); 
		}

	}
	catch(assignment_error){
		string reason; 
		reason = "Attempting to compute RI outside of allowed range in ri_model::SiO2(double wavelength)\n";
		reason += "wavelength = " + toString(wavelength, 3) + " (um)\n"; 
		reason += "Allowed range for AlAs model is [ 0.4, 2.5 ]\n"; 
		exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

double ri_model::SiGeO2(double wavelength)
{
	// Refractive Index of Ge02 Doped Silica
	// Model is based on a Sellmeier polynomial with experimentally derived coefficients
	// Taken from K. Okamoto, "Fundamanetals of Optical Waveguides", 2006
	// wavelength is input in units of microns
	// R. Sheehan 12 - 4 - 2012
	
	try{

		if(wavelength > 0.4 && wavelength < 2.5){

			int ncoeffs = 3; 
	
			double a[4]; 
			double b[4]; 
	
			a[1]=0.7083952; a[2]=0.4203993; a[3]=0.8663412; 
			b[1]=7.290464e-3; b[2]=1.050294e-2; b[3]=97.93428; 
	
			double nsqr=1.0;
			double lsqr = DSQR(wavelength);
	
			for(int i=1; i<=ncoeffs; i++){
				nsqr += ( a[i]*lsqr ) / ( lsqr - b[i] );
			}
	
			return sqrt(nsqr);

		}
		else{
			throw assignment_error(); 
		}

	}
	catch(assignment_error){
		string reason; 
		reason = "Attempting to compute RI outside of allowed range in ri_model::SiO2(double wavelength)\n";
		reason += "wavelength = " + toString(wavelength, 3) + " (um)\n"; 
		reason += "Allowed range for AlAs model is [ 0.4, 2.5 ]\n"; 
		exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

double ri_model::AlGaAs(double alfrac, double wavelength)
{
	// source: http://www.batop.com/information/n_AlGaAs.html
    // type: theoretical fit
    // model valid for at least lambda (um) in [0.65, 2.0]
    // model describes refractive index of Al_{x}Ga_{1-x}As
    // As alfrac -> 0 Al_{x}Ga_{1-x}As -> GaAs
    // As alfrac -> 1 Al_{x}Ga_{1-x}As -> AlAs

	try{
		// Model is only valid on certain wavelength range

		if(wavelength > 0.64 && wavelength < 2.1 && alfrac >= 0.0 && alfrac <= 1.0){

			// PE = HC/wavelength # photon energy in Joules
			double A0 = 6.0 + (19.0*alfrac); // fitting parameter, value dependent on Alfrac
			double B0 = 9.4 - (10.2*alfrac); // fitting parameter, value dependent on Alfrac
			double alfracsqr = DSQR(alfrac);
			double E0 = 1.425 + (1.155*alfrac) + (0.37*alfracsqr); // fundamental band gap at G-point expressed in eV
			double S0 = 1.765 + (1.115*alfrac) + (0.37*alfracsqr); // fundamental band gap at G-point plus D0 spin-orbit splitting energy expressed in eV
			double ER = pow((E0/S0),1.5); // (E0 / E0+D0)^{3/2}

			double chi = HC/(wavelength*E0);
			double chisqr = DSQR(chi);

			double t1 = sqrt(1.0+chi);
			double t2 = sqrt(1.0-chi);
			double fchi = (2.0-t1-t2)/chisqr;

			double chiS0 = HC/(wavelength*S0);
			double chisqrS0 = DSQR(chiS0);

			double t1S0 = sqrt(1.0+chiS0);
			double t2S0 = sqrt(1.0-chiS0);
			double fchiS0 = (2.0-t1S0-t2S0)/chisqrS0;

			double t3 = ( fchi + 0.5*fchiS0*ER );

			return sqrt(A0*t3+B0); 
		}
		else{
			throw assignment_error(); 
		}

	}
	catch(assignment_error){
		string reason; 
		reason = "Attempting to compute RI outside of allowed range in ri_model::AlGaAs(double alfrac, double wavelength)\n";
		reason += "wavelength = " + toString(wavelength, 3) + " (um)\n"; 
		reason += "Allowed range for AlGaAs model is [ 0.65, 2.0 ]\n"; 
		reason += "Or Al fraction is too large: alfrac = " + toString(alfrac, 3) + "\n";
		exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

double ri_model::InGaAs(double infrac, double wavelength)
{
	// source: http://www.batop.com/information/n_InGaAs.html
    // type: theoretical fit
    // model valid for at least lambda (um) in [1.15, 2.0]
    // model describes refractive index of In_{x}Ga_{1-x}As
    // As infrac -> 0 In_{x}Ga_{1-x}As -> GaAs
    // As infrac -> 1 In_{x}Ga_{1-x}As -> InAs

	try{
		// Model is only valid on certain wavelength range

		if(wavelength > 1.14 && wavelength < 2.1 && infrac >= 0.0 && infrac <= 1.0){

			double A = 8.950; // empirical coefficient
			double B = 2.054; // empirical coefficient
			double C = 0.6245; // empirical coefficient
			double EgGaAs = 1.424; // fundamental band gap of GaAs at room temperature (300 K)

			double infracsqr = DSQR(infrac);
			double Eg = EgGaAs - (1.501*infrac) + (0.436*infracsqr);

			double t1 = (C*EgGaAs)/(wavelength*Eg);
			double t2 = 1.0 - DSQR(t1);
			double t3 = B/t2;
			double t4 = A+t3;

			return sqrt(t4); 
		}
		else{
			throw assignment_error(); 
		}

	}
	catch(assignment_error){
		string reason; 
		reason = "Attempting to compute RI outside of allowed range in ri_model::InGaAs(double infrac, double wavelength)\n";
		reason += "wavelength = " + toString(wavelength, 3) + " (um)\n"; 
		reason += "Allowed range for InGaAs model is [ 1.15, 2.0 ]\n"; 
		reason += "Or In fraction is too large: infrac = " + toString(infrac, 3) + "\n";
		exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

double ri_model::GaAsN_1300(double nfrac)
{
	// source: data stripped from various papers
	// valid between 0% and 5% N at a wavelength of 1300 nm
	// see document EAM_Results_21_10_2014_NB.pdf


	try{
		// Model is only valid on certain wavelength range

		if(nfrac >= 0.0 && nfrac < 0.06){

			return (3.41326 + 1.94069*nfrac); 
		}
		else{
			throw assignment_error(); 
		}

	}
	catch(assignment_error){
		string reason; 
		reason = "Attempting to compute RI outside of allowed range in ri_model::InGaAs(double nfrac)\n";
		reason += "wavelength = 1300 (um)\n"; 
		reason += "N fraction outside allowed range: nfrac = " + toString(100*nfrac, 3) + "\n";
		reason += "allowed range 0 < %N < 5\n";
		exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

double ri_model::InGaAsP(double asfrac, double wavelength)
{
	//Model for the refactive index of In_{1-x}Ga_{x}As_{y}P_{1-y}
	//Taken from Broberg and Lindgren, J. Appl. Phys., 55(9), 1984
	//The model is valid for 0 <= y < 1 and for wavelengths above the band-gap wavelength l > lg
	//R. Sheehan 24 - 3 - 2010

	try{
		// Model is only valid on certain wavelength range

		double y = asfrac; // As fraction
		double Eg=( 1.35 - 0.72*y + 0.12*DSQR(y) ); // Band gap energy in units of eV
		double lg = (1.24/Eg); // Band gap wavelength

		if(wavelength > lg && asfrac >= 0.0 && asfrac <= 1.0){
			
			double Egsqr=DSQR(Eg);
			double x=((0.1894*y)/(0.4184-0.013*y)); // Mole fraction x as a function of y
			double En=(1.24/wavelength); // Wavelength in units of energy
			double Ensqr=DSQR(En);
			double Enfourth=DSQR(Ensqr);
			
			double Ed=((12.36*x-12.71)*y+7.54*x+28.91); // Some other parameter
			double Eo=(0.595*DSQR(x)*(1-y)+1.626*x*y-1.891*y+0.524*x+3.391); // Some other parameter
			
			double Eosqr=DSQR(Eo);
			double Eocube=Eo*Eosqr;
			double Et=((PI*Ed)/(2*Eocube*(Eosqr-Egsqr)));

			double v1=(Ed/Eo);
			double v2=((Ed*Ensqr)/Eocube);
			double v3=((Et*Enfourth)/PI);
			double v4=log((2.0*Eosqr-Egsqr-Ensqr)/(Egsqr-Ensqr));

			return sqrt(1.0+v1+v2+v3*v4); // This is the refractive index
		}
		else{
			throw assignment_error(); 
		}

	}
	catch(assignment_error){
		string reason; 
		reason = "Attempting to compute RI outside of allowed range in ri_model::InGaAsP(double asfrac, double wavelength)\n";
		reason += "wavelength = " + toString(wavelength, 3) + " (um)\n"; 
		reason += "Allowed range for InGaAsP model is l > lg(asfrac)\n"; 
		reason += "Or As fraction is too large: asfrac = " + toString(asfrac, 3) + "\n"; 
		exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

double ri_model::InGaAsN(double infrac, double nfrac, double wavelength)
{
	// y is the N fraction
    // 1-y is the As fraction
    // x is the In fraction
    // 1-x is the Ga fraction
    // w is the wavelength

    // this model is not accurate for GaAsN
    // see comparison of measured data scraped from literature
    // model does not reproduce correct RI of GaAs
    // R. Sheehan 15 - 10 - 2014

	try{
		// Model is only valid on certain wavelength range

		if(wavelength > 1.0 && infrac > 0.0 && infrac < 1.0 && nfrac > 0.0 && nfrac < 0.07){
			
			double E0 = (3.9 - 3.75*infrac - 15.9*nfrac );
			double Ed = (36.1 - 19.9*infrac );
			//double E = (1.23985 / (wavelength*1000.0) ); // assumes that wavelength is in nm
			double E = (1.23985 / (wavelength) ); // assumes that wavelength is in um
			double X; 

			if( fabs(E0 - E) > EPS ){
				// E0 != E
				X = (E0*Ed)/( DSQR(E0) - DSQR(E) ); 
			}
			else{
				// E0 = E
				X = (E0*Ed)/(1.0E-30); 
			}

			return sqrt(1.0 + X); 
		}
		else{
			throw assignment_error(); 
		}

	}
	catch(assignment_error){
		string reason; 
		reason = "Attempting to compute RI outside of allowed range in ri_model::InGaAsN(double asfrac, double wavelength)\n";
		reason += "wavelength = " + toString(wavelength, 3) + " (um)\n"; 
		reason += "Allowed range for InGaAsP model is l > lg(asfrac)\n"; 
		reason += "Or In fraction is too large: infrac = " + toString(infrac, 3) + "\n"; 
		reason += "Or N fraction is too large: nfrac = " + toString(nfrac, 3) + "\n"; 
		exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

double ri_model::AlInGaAs(double alfrac, double wavelength)
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
		
		/*a = 10.761 - 2.2617*alfrac - 5.3678*DSQR(alfrac) + 5.4943*alfrac*DSQR(alfrac); 
		b = 0.9485 - 0.2369*alfrac + 5.2374*DSQR(alfrac) - 5.1947*alfrac*DSQR(alfrac); 
		csqr = DSQR( 1.5986 - 3.0266*alfrac + 2.3326*DSQR(alfrac) - 0.3629*alfrac*DSQR(alfrac) );*/ 

		// evaluate the polynomials
		a = A[3]; b = B[3]; csqr = C[3]; 
		for(int j=2; j>=0; j--){
			a = a*alfrac + A[j]; 
			b = b*alfrac + B[j];
			csqr = csqr*alfrac + C[j]; 
		}
		csqr = DSQR(csqr); 

		lsqr = DSQR(wavelength); 

		numer = (b*lsqr); 
		denom = (lsqr - csqr); 
		
		if( fabs(denom) > EPS && alfrac >= 0.0 && alfrac <= 1.0){	
			return sqrt(a + numer / denom);
		}
		else{
			throw assignment_error(); 
		}

	}
	catch(assignment_error){
		string reason; 
		reason = "Attempting to compute RI outside of allowed range in ri_model::AlInGaAs(double alfrac, double wavelength)\n";
		reason += "wavelength = " + toString(wavelength, 3) + " (um)\n"; 
		reason += "Allowed range for AlInGaAs model is l > lg(alfrac)\n"; 
		reason += "Or Al fraction is too large: alfrac = " + toString(alfrac, 3) + "\n"; 
		exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

double ri_model::data_based_ri(int mat_type, string &filename, double wavelength)
{
	// this method reads in the data stored in file filename and performs an interpolation over the data set at position wavelength
	// source for data is one of the following
	// http://refractiveindex.info/ // excel csv files
	// http://www.filmetrics.com/refractive-index-database/ // notepad csv files
	// I'm going to work with the filmetrics data sets as they are in a more convenient format, data is the same in both

	
	try{
		wavelength *= 1000.0; // convert to nm because wavelength data is stored in nanometres

		int wavepos; 

		if(!data_in_memory || material_in_memory != mat_type){
			// read in new data only if no data present or material type has changed

			//cout<<"\nReading data in from memory\n"; 
			
			ri_data.clear(); 

			n_data_points = n_cols = 0; 

			// read the data from the file
			vvd_read_matrix_from_file(filename, ri_data, n_data_points, n_cols); 

			// first line contains no text so actual n_data_points is less than measured from file
			n_data_points -= 1;	

			data_in_memory = true; 

			material_in_memory = mat_type; 
		} 

		// store the wavelength values in the array vals for searching
		double *waves = new(double [n_data_points+1]); 

		for(int i=1; i<=n_data_points; i++){
			*(waves+i) = ri_data[i+1][1]; 
		}

		// find the position in vals that is closest to wavelength
		wavepos = binary_search(waves, n_data_points, wavelength); 
	
		// perform interpolation on ri_data
		if(wavepos != -1){

			double ri_value = 0.0; 
			double delta_ri_value = 0.0; 
		
			// store the wavelength values in the array vals for searching
			double *vals = new(double [n_data_points+1]); 

			for(int i=1; i<=n_data_points; i++){
				*(vals+i) = ri_data[i+1][2]; // store refractive index values
			}

			polint(waves, vals, n_data_points, wavelength, ri_value, delta_ri_value); 

			delete[] vals; 

			return ri_value; 
		}
		else{
			throw assignment_error(); 
		}

		delete[] waves; 

		// Under what circumstances do you de-allocate ri_data? 
		//delete[] ri_data
		//data_in_memory = false; 
		//material_in_memory = -1; 
	}
	catch(assignment_error){
		string reason; 
		reason = "Attempting to compute RI outside of allowed range in ri_model::data_based_ri(string &filename, double wavelength)";
		reason += "wavelength = " + toString(wavelength, 3) + " (um) is not found in  file = " + filename + " \n"; 
		exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}

}

void ri_model::fill_material_list()
{
	// list of materials for which the material model / data is known
	materials.push_back(RI_GaAs); 
	materials.push_back(RI_AlAs);
	materials.push_back(RI_SiO2);
	materials.push_back(RI_SiGeO2);
	materials.push_back(RI_AlGaAs); 
	materials.push_back(RI_InGaAs); 
	materials.push_back(RI_InGaAsP);
	materials.push_back(RI_InGaAsN); 
	materials.push_back(RI_GaAsN1300); 
	materials.push_back(RI_InP); 
	materials.push_back(RI_GaP); 
	materials.push_back(RI_AlN); 
	materials.push_back(RI_Si); 
	materials.push_back(RI_AlInGaAs); 
}