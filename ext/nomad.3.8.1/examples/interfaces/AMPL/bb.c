// black-box code using the AMPL interface library
// (written by D. Orban and lightly modified for NOMAD by S. Le Digabel)

// before compilation: define MODEL_NAME

// debug by defining DISPLAY


#include <stdio.h>
#include <stdlib.h>
#include "asl.h"
#include "nlp.h"
#include "getstub.h"

// Defines
#define CHR (char*)

#define MODEL_NAME "./model/optrisk"

// #define DISPLAY  // to see more outputs

// Global variables
ASL *asl;
fint showgrad = (fint)0;  // A command-line option (integer)
fint showname = (fint)0;  // A command-line option (integer)

keyword keywds[] = {  // MUST appear in alphabetical order!
  KW(CHR"showgrad", L_val, &showgrad, CHR"Evaluate gradient"),
  KW(CHR"showname", L_val, &showname, CHR"Display objective name")
};

Option_Info Oinfo = {
  CHR"miniampl", CHR"Mini AMPL Example",
  CHR"miniampl_options", keywds, nkeywds, 0,
  CHR"0.1", 0, 0, 0, 0, 0, 20091021
};


int main ( int argc, char **argv ) {

  FILE * nl;
  char * stub;
  FILE * point_file;
  char * point_file_name;
  int    point_file_name_size;
  fint   nerror = (fint)0;
  int    n_badvals = 0;
  int    n_con_tmp = 0;
  int    i;
  real   f;
  real * R;

  if( argc < 2 ) {
    fprintf ( stderr , "Usage: %s x.txt\n" , argv[0] );
    return 1;
  }

  // get the point file name:
  point_file_name_size = strlen(argv[1]) + 1;
  point_file_name      = (char*)Malloc(point_file_name_size * sizeof(char));
  strcpy ( point_file_name , argv[1] );
  strcpy ( argv[1] , MODEL_NAME );

  // Read objectives and first derivative information.
  if( !(asl = ASL_alloc(ASL_read_fg)) ) exit(1);
  stub = getstub(&argv, &Oinfo);
  nl   = jac0dim(stub, (fint)strlen(stub));

  // Get command-line options.
  if (getopts(argv, &Oinfo)) exit(1);

  // Check command-line options.
  if( showgrad < 0 || showgrad > 1 ) {
    Printf("Invalid value for showgrad: %d\n", showgrad);
    n_badvals++;
  }
  if( showname < 0 || showname > 1 ) {
    Printf("Invalid value for showgrad: %d\n", showgrad);
    n_badvals++;
  }

  if(n_badvals) {
    Printf("Found %d errors in command-line options.\n", n_badvals);
    exit(1);
  }

  // Allocate memory for problem data.
  // The variables below must have precisely THESE names.
  X0    = (real*)Malloc(n_var * sizeof(real));  // Initial guess
  pi0   = (real*)Malloc(n_con * sizeof(real));  // Initial multipliers
  LUv   = (real*)Malloc(n_var * sizeof(real));  // Lower bounds on variables
  Uvx   = (real*)Malloc(n_var * sizeof(real));  // Upper bounds on variables
  LUrhs = (real*)Malloc(n_con * sizeof(real));  // Lower bounds on constraints
  Urhsx = (real*)Malloc(n_con * sizeof(real));  // Upper bounds on constraints
  R     = (real*)Malloc(n_con * sizeof(real));  // constraints

  want_xpi0 = 3;

  // Read in ASL structure - trap read errors
  if( fg_read(nl, 0) ) {
    fprintf(stderr, "Error fg-reading nl file\n");
    goto bailout;
  }

#ifdef DISPLAY

  n_con_tmp = 0;
  for ( i = 0 ; i < n_con ; ++i ) {
    if ( LUrhs[i] > -Infinity )
      ++n_con_tmp;
    if ( Urhsx[i] < Infinity )
      ++n_con_tmp;
  }

  printf ( "n_obj=%i\nn_var=%i\nn_con=%i\nx0=[" , n_obj , n_var , n_con_tmp );
  for ( i = 0 ; i < n_var ; ++i )
    printf ( "%g " , X0[i] );
  printf ( "]\n" );
#endif

  // read x:
  if ((point_file = fopen(point_file_name,"r")) == NULL) {
    fprintf(stderr, "Cannot open file %s.\n",point_file_name);
    goto bailout;
  }

  for ( i = 0 ; i < n_var ; ++i )
    fscanf ( point_file , "%lf" , &X0[i] );

  fclose(point_file);
  free ( point_file_name );


#ifdef DISPLAY
  printf ( "x =[" );
  for ( i = 0 ; i < n_var ; ++i )
    printf ( "%g " , X0[i] );
  printf ( "]\n" );
#endif

  // objective functions:
  for ( i = 0 ; i < n_obj ; ++i ) {
    f = objval ( i , X0 , &nerror ); 

    if ( nerror ) {
      fprintf(stderr, "Error while evaluating objective.\n");
      goto bailout;
    }

#ifdef DISPLAY
    Printf("f%i(x) = %21.15e\n", i , f );
#else
    Printf("%21.15e\n", f );
#endif
  }

  // constraints:
  conval ( X0 , R , &nerror );

  for ( i = 0 ; i < n_con ; ++i ) {

#ifdef DISPLAY
    printf ("%g <= %g <= %g\n" ,  LUrhs[i] , R[i] , Urhsx[i] );
#else
    if ( LUrhs[i] > -Infinity )
      Printf("%21.15e\n", LUrhs[i]-R[i] );
    if ( Urhsx[i] < Infinity )
      Printf("%21.15e\n", R[i]-Urhsx[i] );
#endif
  }

 bailout:
  // Free data structure. DO NOT use free() on X0, pi0, etc.
  ASL_free((ASL**)(&asl));

  return 0;
}
