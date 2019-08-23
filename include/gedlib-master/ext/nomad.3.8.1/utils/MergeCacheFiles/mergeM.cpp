/*-----------------------------------------------------*/
/*  Merge cache files from NOMAD                       */
/*-----------------------------------------------------*/
#include "mex.h"
#include "nomad.hpp"

using namespace std;
// using namespace NOMAD; avoids putting NOMAD:: everywhere


struct printfbuf : std::streambuf {
public:
     //Constructor
	 printfbuf() 
	 { 
		 setp(m_buffer, m_buffer + s_size - 2); 
	 }	 
private:
    enum { s_size = 1024 }; //not sure on this size
	char m_buffer[s_size];
	int_type overflow(int_type c) 
	{
		if (!traits_type::eq_int_type(c, traits_type::eof())) {
			*pptr() = traits_type::to_char_type(c);
			pbump(1);
		}
		return sync() != -1 ? traits_type::not_eof(c) : traits_type::eof();
	}

	int sync() {
		*pptr() = 0;
		mexPrintf(pbase());
        mexEvalString("drawnow;");
		setp(m_buffer, m_buffer + s_size - 2);
		return 0;
	} 
 };

// Main Entry Function
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // display:

    
    try
    {
        
        if(nrhs != 3 || !mxIsChar(prhs[0]) || !mxIsChar(prhs[1]) || !mxIsChar(prhs[2]) )
            mexErrMsgTxt("You must supply the name of 3 cache files for merging. 2 for input and 1 for output \n");
        
        
        char *cache_file_name_1 = NULL;
        char *cache_file_name_2 = NULL;
        char *cache_file_name_3 = NULL;
        cache_file_name_1 = mxArrayToString(prhs[0]);
        cache_file_name_2 = mxArrayToString(prhs[1]);
        cache_file_name_3 = mxArrayToString(prhs[2]);
        
        printfbuf buf;
        std::streambuf *cout_sbuf = std::cout.rdbuf(); //keep existing buffer
        std::cout.rdbuf(&buf); //redirect buffer
        
        NOMAD::Display out ( std::cout );
        out.precision ( NOMAD::DISPLAY_PRECISION_STD );
        
        // parameters creation:
        NOMAD::Cache cache1 ( out );
        NOMAD::Cache cache2 ( out );
        NOMAD::Cache cache3 ( out );
        
        cache1.load(cache_file_name_1);
        cache2.load(cache_file_name_2);
        cache3.load(cache_file_name_3);
        
        
        cache3.insert(cache1);
        cache3.insert(cache2);
        cache3.save(true);
        
        out << "After merging \n " << cache_file_name_3 << " contains " << cache3.size() << " points" << endl;
        
        
    }
    catch ( exception & e )
    {
        cerr << "\n Merge triggered an exception (" << e.what() << ")\n\n";
    }
  
}
