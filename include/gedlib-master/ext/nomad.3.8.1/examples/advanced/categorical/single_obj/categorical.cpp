/*-------------------------------------------------------------------*/
/*            Example of a problem with categorical variables        */
/*-------------------------------------------------------------------*/
/*                                                                   */
/*  . portfolio problem with 3 assets                                */
/*                                                                   */
/*  . NOMAD is used in library mode                                  */
/*                                                                   */
/*  . the number of variables can be 3,5, or 7, depending on the     */
/*    number of assets considered in the portfolio                   */
/*                                                                   */
/*  . variables are of the form (n t0 v0 t1 v1 t2 v2) where n is the */
/*    number of assets, ti is the type of an asset, and vi is the    */
/*    money invested in this asset                                   */
/*                                                                   */
/*  . categorical variables are n and the ti's                       */
/*                                                                   */
/*  . with a $10,000 budget, the problem consists in minimizing      */
/*    some measure of the risk and of the revenue                    */
/*                                                                   */
/*  . two classes are defined:                                       */
/*                                                                   */
/*    1. My_Evaluator, wich inherits from the NOMAD class Evaluator, */
/*       in order to define the problem via the virtual function     */
/*       eval_x()                                                    */
/*                                                                   */
/*    2. My_Extended_Poll, which inherits from the NOMAD class       */
/*       Extended_Poll, in order to define the categorical           */
/*       variables neighborhoods, via the virtual function           */
/*       construct_extended_points()                                 */
/*                                                                   */
/*  . My_Extended_Poll also defines 3 signatures, for the solutions  */
/*    with 3, 5, and 7 variables                                     */
/*-------------------------------------------------------------------*/
/*  . compile the scalar version with 'make'                         */
/*  . compile the parallel version with 'make mpi'                   */
/*-------------------------------------------------------------------*/
/*  . execute the scalar version with './categorical'                */
/*  . execute the parallel version with 'mpirun -np p ./categorical' */
/*    with p > 1                                                     */
/*-------------------------------------------------------------------*/
#include "nomad.hpp"

using namespace std;
using namespace NOMAD;

#define USE_SURROGATE false

/*----------------------------------------*/
/*               the problem              */
/*----------------------------------------*/
class My_Evaluator : public Evaluator
{
public:
	My_Evaluator ( const Parameters & p ) :
    Evaluator ( p ) {}
	
	~My_Evaluator ( void ) {}
	
	bool eval_x ( Eval_Point   & x          ,
				 const Double & h_max      ,
				 bool         & count_eval   ) const;
};

/*--------------------------------------------------*/
/*  user class to define categorical neighborhoods  */
/*--------------------------------------------------*/
class My_Extended_Poll : public Extended_Poll
{
	
private:
	
	// signatures for 1, 2, and 3 assets:
	Signature * _s1 , * _s2 , * _s3;
	
public:
	
	// constructor:
	My_Extended_Poll ( Parameters & );
	
	// destructor:
	virtual ~My_Extended_Poll ( void ) { delete _s1; delete _s2; delete _s3; }
	
	// construct the extended poll points:
	virtual void construct_extended_points ( const Eval_Point & );
	
};

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv )
{
	
	// NOMAD initializations:
	begin ( argc , argv );
	
	// display:
	Display out ( cout );
	out.precision ( DISPLAY_PRECISION_STD );
	
	// check the number of processess:
#ifdef USE_MPI
	if ( Slave::get_nb_processes() < 2 )
	{
		if ( Slave::is_master() )
			cerr << "usage: \'mpirun -np p ./categorical\' with p>1"
			<< endl;
		end();
		return EXIT_FAILURE;
	}
#endif
	
	try
	{
		
		// parameters creation:
		Parameters p ( out );
		
		if ( USE_SURROGATE )
			p.set_HAS_SGTE ( true );
		
        // p.set_DISPLAY_DEGREE ( FULL_DISPLAY );
		
		p.set_DIMENSION (3);
		
		vector<bb_output_type> bbot (3);
		bbot[0] = EB;  // budget constraint
		bbot[1] = EB;  // total value >= 1$
		bbot[2] = OBJ; // objective
		p.set_BB_OUTPUT_TYPE ( bbot );
		
		
		// categorical variables:
		p.set_BB_INPUT_TYPE ( 0 , CATEGORICAL );
		p.set_BB_INPUT_TYPE ( 1 , CATEGORICAL );
		
		// initial point
		Point x0 ( 3 , 0 );
		x0[0] = 1;     // 1 asset
		x0[1] = 0;     // of type 0
		x0[2] = 100;   // 100$
		p.set_X0 ( x0 );
		
		Point lb ( 3 );
		Point ub ( 3 );
		// Categorical variables 0 and 1 don't need bounds
		lb[2] = 0; ub[2] = 10000;
		
		p.set_LOWER_BOUND ( lb );
		p.set_UPPER_BOUND ( ub );
		
		p.set_DISPLAY_STATS ( "bbe ( sol ) obj" );
		
		// extended poll trigger:
		p.set_EXTENDED_POLL_TRIGGER ( 10 , false );
		
		// parameters validation:
		p.check();
		
		// custom evaluator creation:
		My_Evaluator ev ( p );
		
		// extended poll:
		My_Extended_Poll ep ( p );
		
		// algorithm creation and execution:
		Mads mads ( p , &ev , &ep , NULL , NULL );
		mads.run();
	}
	catch ( exception & e ) {
		string error = string ( "NOMAD has been interrupted: " ) + e.what();
		if ( Slave::is_master() )
			cerr << endl << error << endl << endl;
	}
	
	
	Slave::stop_slaves ( out );
	end();
	
	return EXIT_SUCCESS;
}

/*----------------------------------------------------*/
/*                         eval_x                     */
/*----------------------------------------------------*/
bool My_Evaluator::eval_x ( Eval_Point   & x          ,
						   const Double & h_max      ,
						   bool         & count_eval   ) const
{
	
	// number of assets:
	int n = static_cast<int> ( x[0].value() );
	
	count_eval=false;
	
	// get the asset types and values:
	Point v ( 3 , 0.0 );
	Double vmin = 10000 , tmp;
	for ( int i = 0 ; i < n ; ++i )
	{
		tmp = v [ static_cast<int> ( x[2*i+1].value() ) ] = x[2*i+2];
		if ( tmp < vmin )
			vmin = tmp;
	}
	
	// constraints (budget and each asset is considered with at least 1$):
	Double vt = v[0] + v[1] + v[2];
	Double h  = vt - 10000;
	x.set_bb_output  ( 0 , h );
	x.set_bb_output  ( 1 , 1-vmin );
	
	if ( h <= 0 && vmin >= 1 )
	{
		
		// compute the risk and revenue:
		Double vt2  = vt.pow2();
		Double rev  = v[0] * 0.0891 + v[1] * 0.2137 + v[2] * 0.2346;
		Double risk = 0.01 * (v[0]/vt).pow2() +
		0.05 * (v[1]/vt).pow2() +
		0.09 * (v[2]/vt).pow2() +
		0.02 * (v[0]*v[1]/vt2)  +
		0.02 * (v[0]*v[2]/vt2)  +
		0.10 * (v[1]*v[2]/vt2);
		
		// the objective is taken as a scaled distance
		// between (risk,revenue) and (risk_best,rev_best):
		Double a = ( risk - 0.01 ) * 100 / 0.08;
		Double b = ( rev  - 891  ) * 100 / 1455;
		
		x.set_bb_output  ( 2 , ( a.pow2() + (100-b).pow2() ).sqrt() );
		
		count_eval = true;
	}
	else
		x.set_bb_output  ( 2 , 145 );
	
	// simulation of a surrogate:
	if ( USE_SURROGATE && x.get_eval_type() == SGTE )
	{
		Double f = x.get_bb_outputs()[2];
		f.round();
		f += ( (rand()%2) ? -1.0 : 1.0 ) * ( (rand()%1000) / 1000.0 );
		x.set_bb_output  ( 2 , f );
	}
	
	return true;
}

/*-----------------------------------------*/
/*  constructor: creates the 3 signatures  */
/*-----------------------------------------*/
My_Extended_Poll::My_Extended_Poll ( Parameters & p )
: Extended_Poll ( p    ) ,
_s1           ( NULL ) ,
_s2           ( NULL ) ,
_s3           ( NULL )
{
	
	// signature for 1 asset:
	// ----------------------
	vector<bb_input_type> bbit_1 (3);
	bbit_1[0] = bbit_1[1] = CATEGORICAL;
	bbit_1[2] = CONTINUOUS;
	
	const Point & d0_1 = p.get_initial_poll_size();
	const Point & lb_1 = p.get_lb();
	const Point & ub_1 = p.get_ub();
	
	_s1 = new Signature ( 3                          ,
						 bbit_1                     ,
						 d0_1                       ,
						 lb_1                       ,
						 ub_1                       ,
						 p.get_direction_types   () ,
						 p.get_sec_poll_dir_types() ,
						 _p.out()                     );
	
	// signature for 2 assets:
	// -----------------------
	{
		vector<bb_input_type> bbit_2 (5);
		Point d0_2 (5);
		Point lb_2 (5);
		Point ub_2 (5);
		
		// Categorical variables don't need bounds
		for ( int i = 0 ; i < 5 ; ++i )
		{
			bbit_2[i] = bbit_1[1+(i-1)%2];
			d0_2  [i] = d0_1  [1+(i-1)%2];
			lb_2  [i] = lb_1  [1+(i-1)%2];
			ub_2  [i] = ub_1  [1+(i-1)%2];
		}
		
		_s2 = new Signature ( 5                          ,
							 bbit_2                     ,
							 d0_2                       ,
							 lb_2                       ,
							 ub_2                       ,
							 p.get_direction_types   () ,
							 p.get_sec_poll_dir_types() ,
							 _p.out()                     );
	}
	
	// signature for 3 assets:
	// -----------------------
	{
		vector<bb_input_type> bbit_3 (7);
		Point d0_3 (7);
		Point lb_3 (7);
		Point ub_3 (7);
		
		// Categorical variables don't need bounds
		for ( int i = 0 ; i < 7 ; ++i )
		{
			bbit_3[i] = bbit_1[1+(i-1)%2];
			d0_3  [i] = d0_1  [1+(i-1)%2];
			lb_3  [i] = lb_1  [1+(i-1)%2];
			ub_3  [i] = ub_1  [1+(i-1)%2];
		}
		
		_s3 = new Signature ( 7                          ,
					 bbit_3                     ,
							 d0_3                       ,
							 lb_3                       ,
							 ub_3                       ,
							 p.get_direction_types   () ,
							 p.get_sec_poll_dir_types() ,
							 _p.out()                     );
	}
}

/*--------------------------------------*/
/*  construct the extended poll points  */
/*      (categorical neighborhoods)     */
/*--------------------------------------*/
void My_Extended_Poll::construct_extended_points ( const Eval_Point & x )
{
	
	// number of assets:
	int n = static_cast<int> ( x[0].value() );
	
	// 1 asset:
	// --------
	if ( n==1 )
	{
		
		int cur_type = static_cast<int> ( x[1].value() );
		
		// this vector contains the types of the other assets:
		vector<int> other_types;
		switch ( cur_type )
		{
			case 0:
			case 2:
				other_types.push_back(1);
				break;
			default:
				other_types.push_back(0);
				other_types.push_back(2);
		}
		
		// add 1 asset (1 or 2 neighbors):
		for ( size_t k = 0 ; k < other_types.size() ; ++k )
		{
			Point y (5);
			y[0] = 2;
			y[1] = cur_type;
			y[3] = other_types[k];
			y[2] = y[4] = x[2]/2.0;
			
			add_extended_poll_point ( y , *_s2 );
		}
		
		// change the type of the asset to the other types (1 or 2 neighbors):
		for ( size_t k = 0 ; k < other_types.size() ; ++k )
		{
			Point y = x ;
			y[1] = other_types[k];
			
			add_extended_poll_point ( y , *_s1 );
		}
	}
	
	// 2 assets:
	// ---------
	else if ( n == 2 )
	{
		
		int other_type = static_cast<int> ( (3 - x[1] - x[3]).value() );
		
		// change the type of one asset (2 neighbors):
		{
			Point y1 = x;
			Point y2 = x;
			
			y1[1] = y2[3] = other_type;
			
			add_extended_poll_point ( y1 , *_s2 );
			add_extended_poll_point ( y2 , *_s2 );
		}
		
		// remove 1 asset (2 neighbors):
		{
			Point y1(3) , y2(3);
			y1[0] = y2[0] = 1;
			y1[1] = x[1];
			y2[1] = x[3];
			y1[2] = y2[2] = x[2]+x[4];
			
			add_extended_poll_point ( y1 , *_s1 );
			add_extended_poll_point ( y2 , *_s1 );
		}
		
		// add one asset (1 neighbor):
		{
			Point y(7);
			y[0] = 3;
			
			y[1] = x[1];
			y[3] = x[3];
			y[5] = other_type;
			
			y[2] = y[4] = y[6] = (x[2]+x[4])/3.0;
			
			add_extended_poll_point ( y , *_s3 );
		}
		
	}
	
	// 3 assets:
	// ---------
	else {
		
		// remove one asset (3 neighbors):
		Point y1(5);
		Point y2(5);
		Point y3(5);
		
		y1[0] = y2[0] = y3[0] = 2;
		
		y1[1] = x[1];
		y1[3] = x[3];
		
		y2[1] = x[1];
		y2[3] = x[5];
		
		y3[1] = x[3];
		y3[3] = x[5];
		
		y1[2] = y1[4] = y2[2] = y2[4] = y3[2] = y3[4] = (x[2]+x[4]+x[6]) / 2.0;
		
		add_extended_poll_point ( y1 , *_s2 );
		add_extended_poll_point ( y2 , *_s2 );
		add_extended_poll_point ( y3 , *_s2 );
	}
}
