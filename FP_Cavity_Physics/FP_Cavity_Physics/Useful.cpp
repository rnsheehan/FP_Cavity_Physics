#ifndef ATTACH_H
#include "Attach.h"
#endif

std::string useful_funcs::TheTime()
{
	// Implementation of a function returning the current time as a string
	// This is just a neater way of ensuring that the time can be correctly and easily accessed
	// without being hassled about whether or not you've remembered to use non-deprecated versions 
	// of certain functions
	// R. Sheehan 4 - 7 - 2011
	
	const int N=30;	
	char time_str[N];	
	size_t bytes=( N*sizeof(char) );
	
	time_t rawtime;
	
	struct tm timeinfo;
	struct tm *timeinfo_ptr;
	
	timeinfo_ptr=&timeinfo;
	
	// Get current time information
	time(&rawtime);
	
	localtime_s(timeinfo_ptr,&rawtime);
	
	asctime_s(time_str,bytes,timeinfo_ptr);
	
	// Deprecated calls
	//timeinfo=localtime(&rawtime);
	//asctime(timeinfo);
	
	std::string the_time;
	the_time.append(time_str);
	
	return the_time;
}

void useful_funcs::exit_failure_output(std::string reason)
{
	// Code that creates a file and writes a reason in it why the program crashed
	// If it is called of course
	// Call before using the exit(EXIT_FAILURE) command

	// This function outputs to a file an explanation of why the program exited with an EXIT_FAILURE
	// R. Sheehan 17 - 5 - 2011
	
	// Get current time information
	std::string time = TheTime();

	std::ofstream write; // open file for writing
	
	write.open("Exit_Failure_Explanation.txt",std::ios_base::out|std::ios_base::trunc);
	
	//if(!write){
	//	std::cout<<"You're not going to see this statement\n";
	//	std::cout<<"\n";
	//}
	//else{
	//	//printf ( "Current local time and date: %s", asctime (timeinfo) );
	//	write<<"Program Exit Explanation\n\n";
	//	write<<"Error occurred "<<time<<"\n";
	//	write<<reason<<"\n";
	//	write.close();
	//}

	if( write.is_open() ){
		
		write<<"Program Exit Explanation\n\n";
		write<<"Error occurred: "<<time<<"\n";
		write<<reason<<"\n";

		write.close();
	}
}

void useful_funcs::read_into_vector(std::string &filename, std::vector<double> &data, int &n_pts, bool loud)
{
	// read data from a file into a vector
	// R. Sheehan 11 - 9 - 2017

	try{
		std::ifstream the_file; 
		the_file.open(filename, std::ios_base::in);

		if(the_file.is_open()){

			if(loud) std::cout<<filename<<" opened for reading\n"; 

			double value; 
			n_pts = 0;
			while(the_file >> value){
				data.push_back(value);
				n_pts++; 
			}

			if(loud) std::cout<<template_funcs::toString(n_pts)<<" data were read from "<<filename<<"\n"; 

			the_file.close(); 		
		}
		else{
			std::string reason; 
			reason = "Error: void read_into_vector(std::string &filename, std::vector<double> &data)\n"; 
			reason += "Cannot open: " + filename + "\n"; 
			throw std::invalid_argument(reason); 
		}
		
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE); 
	}
}

//double useful_funcs::test_func(double (*f)(int, int), int a, int b)
//{
//	return (*f)(a, b); 
//}