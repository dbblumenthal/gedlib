/***************************************************************************/
/* CONFIGURATION OF THE OPERATING SYSTEM                                   */
/*                                                                         */
/* If Windows platform, put a "0"                                          */
/* If Linux, Unix or Mac platform, put a "1"                               */
/*                                                                         */
#define PLATFORM 1                                /*                       */
/***************************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>

// #include <iomanip>
// #include <cstring>
// #include <cctype>
#include <cstdlib>
// #include <ctime>
using namespace std;


double arrondi ( double x , int n );


#define ARRONDI 1 // arrondir ou pas

#define DEBUG 0		//choose the debug (1) or normal (0) mode
#define MUTE 0      //choose 1 to avoid all screen output, 0 to have the normal display
#define EPS 1e-5           //the software "0"

// Platform-dependant declarations
#define WIN 0           /* all versions                                    */
#define OTHER 1         /* Unix, Linux, SunOS, MacIntosh                   */
#define MESSAGES "runtime/messages.r"
#define RUNTIME "runtime/"
#define DATA "data/"

//useful constants
#define R 0.0821            //ideal gas constant in atm.l/mol.K
#define MAX_TEMP 3000.0     //Maximal temperature in the process
#define pi 3.14159265358979323846 //the pi number
#define MAX_DIM 30
#define MAX_STREAM 30
#define MAX_UNIT 30

//For the chemical class
#define MAX_ERROR 0
#define MAX_WARNING 10

//For the secant solver
#define TOL_SECANT 1e-3
#define MAX_ITER_SECANT 40

//For the bissection solver
#define TOL_BISSECTION 1e-3
#define MAX_ITER_BISSECTION 40

//For the Newton solver
#define TOL_NEWTON 1e-3
#define MAX_ITER_NEWTON 40
#define STEP_NEWTON 1e-3

//For the Runge-Kutta solver
#define N_INTER 100
#define MAX_ITER_RK N_INTER+1

//For the Wegstein solver
#define MIN_THETA -3.0
#define MAX_THETA 1.0
#define TOL_WEGSTEIN 1e-3
#define MAX_ITER_WEGSTEIN 50

//For the stream flashing
#define TOL_BP 1e-3
#define TOL_DP 1e-3

//For the burner
#define TOL_BURN 1e-4

//For the column
#define MAX_PLATES 500
#define MIN_PLATES 1

//For the flash
#define TOL_FLASH 1e-2

//For the cost estimtiors
#define MS_2001 1094.0
#define MS_YEAR 1139.0
