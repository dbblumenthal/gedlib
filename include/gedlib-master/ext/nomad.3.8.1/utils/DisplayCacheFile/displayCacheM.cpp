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
        
        if( nrhs != 1 || !mxIsChar(prhs[0]) )
            mexErrMsgTxt("You must supply the name of a cache file for display \n");
        
        
        char *cache_file_name_1 = NULL;
        cache_file_name_1 = mxArrayToString(prhs[0]);
        
        printfbuf buf;
        std::streambuf *cout_sbuf = std::cout.rdbuf(); //keep existing buffer
        std::cout.rdbuf(&buf); //redirect buffer
        
        NOMAD::Display out ( std::cout );
        out.precision ( NOMAD::DISPLAY_PRECISION_STD );
        
        // parameters creation:
        NOMAD::Cache cache1 ( out );
       
        cache1.load(cache_file_name_1);
        
        const NOMAD::Eval_Point * cur = cache1.begin();
        int  nb = cache1.size();
        int cnt = 0;
        while ( cur )
        {
            out << "point ";
            out.display_int_w ( ++cnt , nb );
            out << "/" << nb << ": ";
            cur->display_eval ( out , false );
            out << std::endl;
            cur = cache1.next();
        }
        out.close_block();

        
    }
    catch ( exception & e )
    {
        cerr << "\n Cache display triggered an exception (" << e.what() << ")\n\n";
    }
  
}
