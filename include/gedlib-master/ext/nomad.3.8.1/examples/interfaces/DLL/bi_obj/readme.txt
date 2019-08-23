
This example illustrates how to use NOMAD with a black-box
problem that is coded inside a DLL file.

The problem has two objective values to minimize.

The example runs on windows with visual C++ and minGW.


DLL characteristics:

/*----------*/

/*  bb.dll  */

/*----------*/



// init/clear:

// -----------

void init_dll  ( void );

void clear_dll ( void );

// problem characteristics:

// ------------------------


int get_n_dll ( void ); // number of variables

int get_m _dll( void ); // number of outputs: output #0 and #1: objective values to minimize

                        //                    output #2..m-1: PB constraints with format g(x) <= 0



// evaluation: f(x):

// -----------------


void eval_x_dll ( double * x , double * fx , bool & count_eval );



//          x: input , double array of size n

//         fx: output, double array of size m

// count_eval: output, set to true if the evaluation has to be counted, false else
