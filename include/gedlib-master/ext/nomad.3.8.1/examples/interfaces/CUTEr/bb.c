
/* ====================================================
 * CUTEr interface for generic package     Feb. 3, 2003
 *
 * D. Orban
 *
 * Take a look at $CUTER/common/include/cuter.h     and
 * $CUTER/common/src/tools/loqoma.c  for more examples.
 * ====================================================
 */

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

#include "cuter.h"

  integer CUTEr_nvar;        /* number of variables */
  integer CUTEr_ncon;        /* number of constraints */

  int MAINENTRY( int argc , char ** argv ) {

    FILE * f;

    char   *fname = "OUTSDIF.d"; /* CUTEr data file */
    integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
    integer iout  = 6;          /* FORTRAN unit number for error output */
    integer ierr;              /* Exit flag from OPEN and CLOSE */


    integer ncon_dummy;
    doublereal *x, *bl, *bu;
    doublereal *v = NULL, *cl = NULL, *cu = NULL;
    logical *equatn = NULL, *linear = NULL;
    logical efirst = FALSE_, lfirst = FALSE_, nvfrst = FALSE_, grad;
    logical constrained = FALSE_;

    real calls[7], cpu[2];
    integer nlin = 0, nbnds = 0, neq = 0;
    doublereal dummy;
    integer ExitCode;
    int i;

    /* Open problem description file OUTSDIF.d */
    ierr = 0;
    FORTRAN_OPEN( &funit, fname, &ierr );
    if( ierr ) {
      /* printf("Error opening file OUTSDIF.d.\nAborting.\n"); */
      printf ( "1e+20 " );
      for ( i = 0 ; i < CUTEr_ncon ; i++ )
	printf ( "1e+20 " );
      printf ("\n");
      exit(1); 
    }

    /* Determine problem size */
    CDIMEN( &funit, &CUTEr_nvar, &CUTEr_ncon );

    /* Determine whether to call constrained or unconstrained tools */
    if( CUTEr_ncon ) constrained = TRUE_;

    /* Reserve memory for variables, bounds, and multipliers */
    /* and call appropriate initialization routine for CUTEr */
    MALLOC( x,      CUTEr_nvar, doublereal );
    MALLOC( bl,     CUTEr_nvar, doublereal );
    MALLOC( bu,     CUTEr_nvar, doublereal );
    
    if( constrained ) {
      MALLOC( equatn, CUTEr_ncon+1, logical    );
      MALLOC( linear, CUTEr_ncon+1, logical    );
      MALLOC( v,      CUTEr_ncon+1, doublereal );
      MALLOC( cl,     CUTEr_ncon+1, doublereal );
      MALLOC( cu,     CUTEr_ncon+1, doublereal );
      ncon_dummy = CUTEr_ncon + 1;
      CSETUP( &funit, &iout, &CUTEr_nvar, &CUTEr_ncon, x, bl, bu,
	      &CUTEr_nvar, equatn, linear, v, cl, cu, &ncon_dummy,
	      &efirst, &lfirst, &nvfrst );
    }
    else {
      MALLOC( equatn, 1, logical    );
      MALLOC( linear, 1, logical    );
      MALLOC( cl, 1, doublereal );
      MALLOC( cu, 1, doublereal );
      USETUP( &funit, &iout, &CUTEr_nvar, x, bl, bu, &CUTEr_nvar );
    }

    ierr = 0;
    FORTRAN_CLOSE( &funit, &ierr );

    /*----------------------------------------------------------*/
/*     printf(" # variables             = %-10d\n", CUTEr_nvar); */
/*     printf(" # constraints           = %-10d\n", CUTEr_ncon); */
/*     for ( i = 0 ; i < CUTEr_nvar ; i++ ) */
/*       printf( "var %d : lb=%g, x0=%g, ub=%g\n",i,bl[i],x[i],bu[i]); */
/*     printf("\n"); */
/*     for ( i = 0 ; i < CUTEr_ncon ; i++ ) */
/*       printf( "cstr %d : equatn=%d, cl=%g, cu=%g\n", */
/* 	      i,equatn[i],cl[i],cu[i]); */
/*     printf("\n"); */
    /*----------------------------------------------------------*/

    if( argc != 2 ) {
      printf ( "1e+20 " );
      for ( i = 0 ; i < CUTEr_ncon ; i++ )
	printf ( "1e+20 " );
      printf ("\n");
      exit(1); 
    }

    f = fopen ( argv[1] , "r" );
    
    for ( i = 0 ; i < CUTEr_nvar ; i++ )
      fscanf ( f , "%lf" , &x[i] );

/*     for ( i = 0 ; i < CUTEr_nvar ; i++ ) */
/*       printf("x[%d]=%g\n",i,x[i]); */

    fclose(f);

    if( constrained ) {
      CFN ( &CUTEr_nvar , &CUTEr_ncon , x , &dummy , &CUTEr_ncon , v );
      printf("%-0.8lf ",dummy);
      for ( i = 0 ; i < CUTEr_ncon ; i++ )
	printf("%-0.8lf ",v[i]);
      printf("\n");
    }
    else {
      UFN ( &CUTEr_nvar , x , &dummy );
      printf("%-15.15g\n",dummy);
    }


    /* Free workspace */

    FREE( x );
    FREE( bl );
    FREE( bu );
    FREE( v );
    FREE( cl );
    FREE( cu );
    FREE( equatn );
    FREE( linear );

    return 0;
    
  }

#ifdef __cplusplus
}    /* Closing brace for  extern "C"  block */
#endif
