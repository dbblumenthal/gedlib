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
        
        // parameters creation:
        NOMAD::Cache cache1 ( out );
        NOMAD::Cache cache2 ( out );
        NOMAD::Cache cache3 ( out );
        
        cache1.load("cache1.bin");
        out << "Before merging \n cache1.bin contains " << cache1.size() << " points" << endl;
        
        cache2.load("cache2.bin");
        out << " cache2.bin contains " << cache2.size() << " points" << endl;
        
        cache3.load("cache3.bin");
        out << " cache3.bin contains " << cache3.size() << " points" << endl;

        
        
        cache3.insert(cache1);
        cache3.insert(cache2);
        
        cache3.save(true);
        
        out << "After merging \n cache3.bin contains " << cache3.size() << " points" << endl;
        
        
    }
    catch ( exception & e )
    {
        cerr << "\n Merge triggered an exception (" << e.what() << ")\n\n";
    }
    
    return EXIT_SUCCESS;
}
