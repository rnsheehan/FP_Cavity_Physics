#ifndef TEMPLATES_H
#define TEMPLATES_H

// Template fucntions for commonly performed mathematical operations
// R. Sheehan 25 - 3 - 2013

namespace template_funcs{
	
	template <class T> bool LimitSafe(T a)
	{
		// check whether a given value is within acceptable numeric limits
		return (a < std::numeric_limits<T>::max() && a > std::numeric_limits<T>::min() ) ? true : false; 
	}

	template <class T> T Signum(T a)
	{
		// The sign operator
		T darg;
		//return ((darg=(a))==0.0?0.0:(darg=(a))>=0.0?1.0:-1.0);
		return ( (darg=(a)) >= (T)(0) ? (T)(1) : -(T)(1) ); // Setting the Sign of zero to be 1
	}

	template <class T> T DSQR(T a)
	{
		// Efficient squaring operator
		// Write injuries in dust, benefits in marble
		T darg;
		return ( (darg=(a)) == (T)(0) ? (T)(0) : darg*darg );
	}

	template <class T> T SIGN(T a,T b)
	{
		return ((b)>(T)(0)?fabs(a):-fabs(a));
	}

	template <class T> std::string toString(const T & t)
	{
		// This is just too convenient not to use
		// Is there a version that can include something similar to %0.5d ? 
		// There has to be, look into the setw method for strings
		// Requires the string-stream (sstream) standard library
		// R. Sheehan 16 - 5 - 2011
    
		std::ostringstream oss; // create a stream
		oss << t;				// insert value to stream
		return oss.str();		// return as a string
	}

	template <class T> std::string toString(const T &t,int places)
	{
		// toString function that allows for the
		// number of decimal places to be specified 
		// far too convenient
		// R. Sheehan 17 - 5 - 2011

		std::ostringstream oss; // create a stream

		oss<<std::fixed<<std::setprecision(places)<<t; // insert value to stream

		return oss.str(); // return as a string
	}

	template <typename T> void sort2(std::vector<T> &vec1, std::vector<T> &vec2)
	{
		// Sort vec1 into ascending order using while making the corresponding rearrangement of vec2
		// This is an implementation of the NRinC::sort2 using std::sort
		// vec1 and vec2 must have the same number of elements
		// R. Sheehan 10 - 1 - 2018

		// this template function can take vector of arbitrary type as argument, but since it has no return value
		// it must be declared as type void to ensure it compiles correctly
		// see https://stackoverflow.com/questions/1640758/passing-stdvector-for-any-type-to-a-function

		try{
			bool c1 = vec1.size() == vec2.size() ? true : false; // vec1 and vec2 must have the same dimensions!
			bool c2 = vec1.size() > 0 ? true : false; 
			bool c3 = vec2.size() > 0 ? true : false; 
			bool c4 = (c1 && c2 && c3) ? true : false;

			if(c4){
				// Join the data sets together to form a single container
				// use the pair data type rather than std::vector<std::vector<double>> Data; 
				std::vector< std::pair< T, T > > Data; 

				// Store the data from vec1 and vec2 in Data
				// vec1 and vec2 must have the same dimensions!
				// vec1 corresponds to the first column, vec2 to the second				
				//std::vector<T>::iterator it1 = vec1.begin(); 
				//std::vector<T>::iterator it2 = vec2.begin();
				std::iterator it1 = vec1.begin(); 
				std::iterator it2 = vec2.begin();
				for(it1; it1 != vec1.end(); it1++){
					Data.push_back( std::make_pair(*it1, *it2) ); 
					it2++; 
				}

				// sort the data using std::sort
				// this will sort according to the data in the first element of the pair
				// the same rearrangement will be applied to the data in the second element
				std::sort(Data.begin(), Data.end()); 

				// place the sorted data back into the original containers		
				//std::vector< std::pair<T,T> >::iterator pit = Data.begin(); 
				std::iterator pit = Data.begin();
				int i=0; 
				for(pit; pit != Data.end(); pit++){
					vec1[i] = pit->first; vec2[i] = pit->second; 
					i++; 
				}

				// no longer need Data
				Data.clear(); 
			}
			else{
				std::string reason; 
				reason = "Error: template_funcs::Sort2\n"; 
				if(!c1) reason+="vec1 and vec2 have different sizes\n"; 
				if(!c2) reason+="vec1 has no elements\n";
				if(!c3) reason+="vec1 has no elements\n";
				throw std::invalid_argument(reason); 
			}
		}
		catch(std::invalid_argument &e){
			std::cerr<<e.what(); 
		}
	}
}

#endif