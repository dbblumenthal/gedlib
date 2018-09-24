/*--------------------------------------------------------------------------*/
/*  example of a program that makes NOMAD restarts after failed iterations  */
/*--------------------------------------------------------------------------*/
#include "nomad.hpp"
using namespace std;
using namespace NOMAD;

/*----------------------------------------*/
/*               the problem              */
/*----------------------------------------*/
class My_Evaluator : public Evaluator
{
private:
    
    
public:
    
    My_Evaluator  ( const Parameters & p ) : Evaluator ( p ) { }
    
    
    
    ~My_Evaluator ( void ) { }
    
    
    virtual bool eval_x ( Eval_Point   & x          ,
                         const Double & h_max      ,
                         bool         & count_eval   ) const;
    
    virtual void update_iteration ( success_type              success      ,
                                   const Stats             & stats        ,
                                   const Evaluator_Control & ev_control   ,
                                   const Barrier           & true_barrier ,
                                   const Barrier           & sgte_barrier ,
                                   const Pareto_Front      & pareto_front ,
                                   bool                    & stop           );
};

/*----------------------------------------*/
/*           user-defined eval_x          */
/*----------------------------------------*/
bool My_Evaluator::eval_x ( Eval_Point   & x          ,
                           const Double & h_max      ,
                           bool         & count_eval   ) const
{
    Double c1 = 0.0 , c2 = 0.0;
    for ( int i = 0 ; i < 5 ; i++ )
    {
        
        c1 += (x[i]-1).pow2();
        c2 += (x[i]+1).pow2();
    }
    x.set_bb_output  ( 0 , x[4]  ); // objective value
    x.set_bb_output  ( 1 , c1-25 ); // constraint 1
    x.set_bb_output  ( 2 , 25-c2 ); // constraint 2
    
    count_eval = true; // count a black-box evaluation
    
    return true;       // the evaluation succeeded
}

/*----------------------------------------*/
/*       updates at each iteration        */
/*----------------------------------------*/
void My_Evaluator::update_iteration ( success_type              success      ,
                                     const Stats             & stats        ,
                                     const Evaluator_Control & ev_control   ,
                                     const Barrier           & true_barrier ,
                                     const Barrier           & sgte_barrier ,
                                     const Pareto_Front      & pareto_front ,
                                     bool                    & stop           )
{
    
    if ( success == UNSUCCESSFUL  && stats.get_bb_eval() > 10 )
        stop = true;
}

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv )
{
    // display:
    Display out ( std::cout );
    out.precision ( DISPLAY_PRECISION_STD );
    
    try {
        
        // NOMAD initializations:
        begin ( argc , argv );
        
        // parameters creation:
        Parameters p ( out );
        
        p.set_DIMENSION (5);             // number of variables
        
        vector<bb_output_type> bbot (3); // definition of
        bbot[0] = OBJ;                   // output types
        bbot[1] = EB;
        bbot[2] = EB;
        p.set_BB_OUTPUT_TYPE ( bbot );
        
        p.set_DISPLAY_DEGREE ( 2 );
        
        p.set_DISPLAY_STATS ( "bbe ( sol ) obj" );
        
        p.set_X0 ( Point ( 5 , 0.0 ) );  // starting point
        
        p.set_LOWER_BOUND ( Point ( 5 , -6.0 ) ); // all var. >= -6
        Point ub ( 5 );                           // x_4 and x_5 have no bounds
        ub[0] = 5.0;                              // x_1 <= 5
        ub[1] = 6.0;                              // x_2 <= 6
        ub[2] = 7.0;                              // x_3 <= 7
        p.set_UPPER_BOUND ( ub );
       

 
        p.set_MAX_BB_EVAL (100); // the algorithm terminates after
        // 100 black-box evaluations
        
        // p.set_DISABLE_MODELS();
        
        // parameters validation:
        p.check();
        
        OrthogonalMesh * oMesh=p.get_signature()->get_mesh();
        
        // custom evaluator creation:
        My_Evaluator ev ( p );
        
        // best solutions:
        const Point * bf = NULL , * bi = NULL;
        
        // algorithm creation:
        Mads mads ( p , &ev );
        
        
        // successive runs:
        for ( int i = 0 ; i < 6 ; ++i )
        {
            
            out << endl << open_block ( "MADS run #" + NOMAD::itos(i) );
            
            // not for the first run:
            if ( i > 0 )
            {
                
                // new starting points:
                p.reset_X0();
                if ( bf )
                    p.set_X0 ( *bf );
                if ( bi )
                    p.set_X0 ( *bi );
                
                if (!bf && !bi)
                    p.set_LH_SEARCH(1,0);  // at least one evaluation is conducted if point is bf and bi are null
                else
                    p.set_LH_SEARCH(0,0);
                
                // Update the mesh for an unsuccessful iteration and put current
                // mesh and poll sizes as initial mesh and poll sizes for the next start.
                oMesh->update(UNSUCCESSFUL);
                Point delta_0, Delta_0;
                oMesh->get_delta(delta_0);
                oMesh->set_delta_0(delta_0);
                oMesh->get_Delta(Delta_0);
                oMesh->set_Delta_0(Delta_0);
                
                // parameters validation:
                p.check();
                
                // reset the Mads object:
                mads.reset ( );
                
            }
            
            // the run:
            mads.run();
            bf = mads.get_best_feasible();
            bi = mads.get_best_infeasible();
            
            out.close_block();
        }
    }
    catch ( exception & e ) 
    {
        
        cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
    }
    
    Slave::stop_slaves ( out );
    end();
    
    return EXIT_SUCCESS;
}
