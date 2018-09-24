/*-----------------------------------------------------*/
/*  how to use the NOMAD library with a user function  */
/*-----------------------------------------------------*/
#include "nomad.hpp"
using namespace std;
// using namespace NOMAD; avoids putting NOMAD:: everywhere


/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv )
{
    
    // display:
    NOMAD::Display out ( std::cout );
    out.precision ( NOMAD::DISPLAY_PRECISION_STD );
    
    try
    {
        
        // NOMAD initializations:
        NOMAD::begin ( argc , argv );
        
        if ( argc != 2 )
        {
            out << "\n Provide a cache file name on the commande line" <<endl;
            return EXIT_FAILURE;
        }
        
        // parameters file:
        std::string cache_file_name = argv[1];
        
        NOMAD::Display out ( std::cout );
        out.precision ( NOMAD::DISPLAY_PRECISION_STD );
        
        // parameters creation:
        NOMAD::Cache cache1 ( out );
        
        cache1.load(cache_file_name);
        
        const NOMAD::Eval_Point * cur = cache1.begin();
        int  nb = cache1.size();
        int cnt = 0;
        while ( cur )
        {
            out << "point ";
            out.display_int_w ( ++cnt , nb );
            out << "/" << nb << ": ";
            cur->display_eval ( out , false );
            out << " obj=" << cur->get_f() << std::endl;
            
            cur = cache1.next();
        }
        out.close_block();
        
        
    }
    catch ( exception & e )
    {
        cerr << "\n Cache display triggered an exception (" << e.what() << ")\n\n";
    }
    
        
     return EXIT_SUCCESS;
}
