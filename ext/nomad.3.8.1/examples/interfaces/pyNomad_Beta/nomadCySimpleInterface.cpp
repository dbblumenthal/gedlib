/*------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search -             */
/*          version 3.8.1                                                       */
/*                                                                              */
/*  NOMAD - version 3.8.1 has been created by                                   */
/*                 Charles Audet        - Ecole Polytechnique de Montreal       */
/*                 Sebastien Le Digabel - Ecole Polytechnique de Montreal       */
/*                 Christophe Tribes    - Ecole Polytechnique de Montreal       */
/*                                                                              */
/*  The copyright of NOMAD - version 3.8.1 is owned by                          */
/*                 Sebastien Le Digabel - Ecole Polytechnique de Montreal       */
/*                 Christophe Tribes    - Ecole Polytechnique de Montreal       */
/*                                                                              */
/*  NOMAD v3 has been funded by AFOSR, Exxon Mobil, Hydro Québec, Rio Tinto     */
/*  and IVADO.                                                                  */
/*                                                                              */
/*  NOMAD v3 is a new version of NOMAD v1 and v2. NOMAD v1 and v2 were created  */
/*  and developed by Mark Abramson, Charles Audet, Gilles Couture, and John E.  */
/*  Dennis Jr., and were funded by AFOSR and Exxon Mobil.                       */
/*                                                                              */
/*  Contact information:                                                        */
/*    Ecole Polytechnique de Montreal - GERAD                                   */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada           */
/*    e-mail: nomad@gerad.ca                                                    */
/*    phone : 1-514-340-6053 #6928                                              */
/*    fax   : 1-514-340-5665                                                    */
/*                                                                              */
/*  This program is free software: you can redistribute it and/or modify it     */
/*  under the terms of the GNU Lesser General Public License as published by    */
/*  the Free Software Foundation, either version 3 of the License, or (at your  */
/*  option) any later version.                                                  */
/*                                                                              */
/*  This program is distributed in the hope that it will be useful, but WITHOUT */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License */
/*  for more details.                                                           */
/*                                                                              */
/*  You should have received a copy of the GNU Lesser General Public License    */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.        */
/*                                                                              */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad        */
/*------------------------------------------------------------------------------*/

#define NOMAD_PYTHON_VERSION "1.0 (beta)  [Feb, 3rd, 2017]"

#include "nomad.hpp"
#include "defines.hpp"
#include <stdio.h>
#include <string.h>

using namespace std;

typedef int (*Callback)(void * apply, NOMAD::Eval_Point &sv);
typedef int (*CallbackL)(void * apply, std::list<NOMAD::Eval_Point *> &sv);

//Print Nomad general Information
//Print Nomad general Information
static void printPyNomadVersion()
{
    printf("-----------------------------------------------------------\n");
    printf("\n Python Interface version %s: C. Tribes \n",NOMAD_PYTHON_VERSION);
    printf("-----------------------------------------------------------\n");
}


static void printPyNomadInfo()
{
    printf("\n-----------------------------------------------------------\n");
    printf(" NOMAD: Nonlinear Optimization using the MADS Algorithm [v%s]\n",NOMAD::VERSION.c_str());
    printf("  - Released under the GNU Lesser General Public License: http://www.gnu.org/copyleft/lesser.html\n");
    printf("  - Source available from: https://www.gerad.ca/nomad/\n");
    
    printPyNomadVersion();
    
    printf(" % NOMAD solves a Global MINLP/NLP in the form               \n");
    printf("          min f(x)                                           \n");
    printf("          subject to:                                        \n");
    printf("                   nlcon(x) <= 0                             \n");
    printf("                   lb <= x <= ub                             \n");
    printf("                   x in R                                    \n");
    printf("\n-----------------------------------------------------------\n");
    
}

//Print Nomad usage
static void printPyNomadUsage()
{
	printf("\n----------------------------------------------------------------------------------------------------------------------\n");
    printf(" PyNomad interface usage:\n");
    printf("------------------------------------------------------------------------------------------------------------------------\n");
    printf(" Info        : PyNomad.info()\n");
    printf(" Version     : PyNomad.version()\n");
    printf(" Help        : PyNomad.help( 'keywords' or 'all' or empty ) \n");
    printf(" Usage       : PyNomad.usage() \n");
    printf("------------------------------------------------------------------------------------------------------------------------\n");
    printf(" Run NOMAD   : [x_best,f_best,h_best,nb_evals,nb_iters,exit_status] = PyNomad.optimize(bb,x0,lb,ub,param)               \n\n");
    printf("               Input parameters:                                                                                        \n");
    printf("                  bb      ---> name of the blackbox function (see below for details) {not optional}                     \n");
    printf("                  x0      ---> list of values for initial point (e.g. x0=[0,0,0,0,0]) {[]}                              \n");
    printf("                  lb      ---> list of values for initial point (e.g. lb=[-1,-2,-3,-4,-5]) {[]}                         \n");
    printf("                  ub      ---> list of values for initial point (e.g. ub=[1,2,3,4,5]) {[]}                              \n");
    printf("                  params  ---> list of Nomad parameters (e.g. param =['DIMENSION 5','BB_OUTPUT_TYPE OBJ']) {not optional}\n");
    printf("                              Consult Nomad documentation or type PyNomad.help() to see all parameters                 \n\n");
    printf("               Output parameters:                                                                                        \n");
    printf("                  x_best   ---> list of values for the best feasible or infeasible point at the end of Nomad optimization\n");
    printf("                  f_best   ---> Objective function value for the best point obtained                                     \n");
    printf("                  h_best   ---> Infeasibility measure value for best point obtained ( 0 if feasible)                     \n");
    printf("                  nb_evals ---> Number of blackbox evaluations                                                           \n");
    printf("                  nb_iters ---> Number of iterations of the Mads algorithm                                               \n");
    printf("                  exit_status  ---> Exit status for Nomad termination criterion (see stop_type in $NOMAD_HOME/src/defines.hpp)\n\n");
    printf("               Blackbox evaluation (example for form 1): a single point is passed to blackbox function at a time         \n");
    printf("                  def bb(x):                                                                                             \n");
    printf("                     f = sum([x.get_coord(i)**2 for i in {0,1,2}])    # the coordinate of the point are accessed         \n");
    printf("                                                                      # with the get_coord function                      \n");
    printf("                     x.set_bb_output(0, f )    # The first output is set. The outputs are given in same order            \n");
    printf("                                               # as defined with the BB_OUTPUT_TYPE parameters                           \n");
    printf("                     return 1    # 1 ->success, 0 -> failed                                                            \n\n");
    printf("               Blackbox evaluation for multiprocess (example for form 2): several bb can be called in parallel           \n");
    printf("                  def bb(x,bb_out):        # blackbox function requires 2 arguments                                      \n");
    printf("                                           # First argument is similar to form 1                                         \n");
    printf("                                           # Second argument is a list of the blackbox outputs                           \n"); 
    printf("                     bb_out.put(sum([x.get_coord(i)**2 for i in {0,1,2}]))   # bb_out is a queue to work in multiprocess \n");    
    printf("                                                                             # The put function is thread safe           \n");
    printf("------------------------------------------------------------------------------------------------------------------------ \n"); 
}

//Print Nomad usage
static void printNomadHelp( string about )
{
	NOMAD::Parameters p ( cout );
	p.help ( about );
}


//Python Evaluator Class
class pyEval : public NOMAD::Evaluator
{
private:
	Callback cb;
	CallbackL cbL;
	void * apply;

    
public:
    //Constructor
    pyEval(const NOMAD::Parameters &p, Callback _cb, CallbackL _cbL, void * _apply) : cb(_cb),cbL(_cbL),apply(_apply),NOMAD::Evaluator(p){}
      
    //Destructor
    ~pyEval(void) {}    
    
	bool eval_x(NOMAD::Eval_Point &x, const NOMAD::Double &h_max, bool &count_eval)
	{
		int success=0;

		//Call Python blackbox function on a single Eval_Point
		try
		{
			count_eval = true;           
			success = cb(apply,x);
            if ( success == -1 )
            {
                printf("Unrecoverable Error from Objective / Blackbox Callback, Exiting NOMAD...\n\n");
                //Force exit
                raise(SIGINT);
                return false;
            }
        }
		//Note if these errors occur it is due to errors in python code
		catch(...)
		{
			printf("Unrecoverable Error from Objective / Blackbox Callback, Exiting NOMAD...\n\n");
			//Force exit
			raise(SIGINT);
			return false;
		}
		if ( success == 1 )
			return true;	
		else 
			return false;
	}
	
	bool eval_x ( std::list<NOMAD::Eval_Point *>	&list_x  ,
				 const NOMAD::Double				& h_max,
				 std::list<bool>					& list_count_eval ) const
	{

		int success=0;
	
		//Call Python blackbox function on a list of Eval_Points
		try
		{
			success = cbL(apply,list_x);
            if ( success == -1 )
            {
                printf("Unrecoverable Error from Objective / Blackbox Callback, Exiting NOMAD...\n\n");
                //Force exit
                raise(SIGINT);
                return false;
            }
		
			std::list<bool>::iterator it_count;
			for (it_count=list_count_eval.begin(); it_count!=list_count_eval.end(); ++it_count)
				(*it_count)=true;

		}
		//Note if these errors occur it is due to errors in python code
		catch(...)
		{
			printf("Unrecoverable Error from Objective / Blackbox Callback, Exiting NOMAD...\n\n");
			//Force exit
			raise(SIGINT);
			return false;
		}
		if ( success == 1 )
			return true;	
		else 
			return false;
	}

  
  
 
};


static int runNomad(Callback cb, CallbackL cbL, void * apply, std::vector<double> X0 , std::vector<double> LB, std::vector<double> UB , const std::vector<std::string> & param , NOMAD::Eval_Point *& best_feas_sol , NOMAD::Eval_Point *& best_infeas_sol , int & nb_evals, int & nb_iters )
{
    
    NOMAD::Display out(std::cout);
    NOMAD::Parameters p(out);
      
    size_t nobj; 
    int dimension;
    bool dimension_is_defined = false;
       
    try
    {    
		size_t ndec = X0.size();	
		if ( ndec != 0 )
		{
			NOMAD::Point px0(ndec);
		
     	    dimension = (int)ndec;
     	    dimension_is_defined = true;
		
			for(int i=0 ; i < ndec ; i++)
				px0[i] = X0[i];
		 
			p.set_X0(px0); 
		}  
	
		size_t nlb = LB.size();	
		if ( nlb != 0 )
		{    	        
			NOMAD::Point plb( (int)nlb );
			
			if ( ! dimension_is_defined )
			{
				dimension = (int)nlb;
     	    	dimension_is_defined = true;
     	    }
			
			if ( ndec != 0 && nlb !=ndec )
				throw NOMAD::Exception("",0,"The lower bound size is inconsistent with X0 size");    			
				
			for(int i=0 ; i < nlb ; i++)
				plb[i] = LB[i];
			
			p.set_LOWER_BOUND( plb );	
		}
		
		size_t nub = UB.size();	
		if ( nub != 0 )
		{         	
			NOMAD::Point pub( (int)nub );
			
			if ( ! dimension_is_defined )
			{
				dimension = (int)nub;
     	    	dimension_is_defined = true;
     	    }
  		
  		    if ( nlb != 0 && nub != nlb )
				throw NOMAD::Exception("",0,"The upper bound size is inconsistent with lower bound size");    			

			if ( ndec != 0 && nub != ndec )
				throw NOMAD::Exception("",0,"The upper bound size is inconsistent with X0 size");    			
  
			for(int i=0 ; i < nub ; i++)
				pub[i] = UB[i];
			
			p.set_UPPER_BOUND( pub );
		
		}
		
		if ( dimension_is_defined )
			p.set_DIMENSION( dimension );
		
		// Make sure that evaluation numbers are reset
		NOMAD::Eval_Point::reset_tags_and_bbes();
	
		// The seed will always be to its default value
		NOMAD::RNG::reset_private_seed_to_default();
		 
		NOMAD::Parameter_Entries entries;
		NOMAD::Parameter_Entry * pe;
		
		// Interpret parameters given in vector
		if ( param.size() > 0 )
		{
			for (int i=0; i < param.size() ; i++ )
			{
				pe = new NOMAD::Parameter_Entry ( param[i] );
				entries.insert( pe );
			}
			p.read( entries );
		}
		
		p.check();
	    	
		nobj = p.get_nb_obj();
		
		if ( nobj != 1 )
			throw NOMAD::Exception ( "" , 0 ,"The python interface for BiMads (biobjective) is not yet implemented");
			
		int noutput = p.get_bb_output_type().size();
			
    	// Create the best solutions (delete in Cython interface    	
    	best_feas_sol =new NOMAD::Eval_Point(p.get_dimension(),static_cast<int>(noutput));
    	best_infeas_sol =new NOMAD::Eval_Point(p.get_dimension(),static_cast<int>(noutput));

			
	}
    catch(exception &e)
    {
        printf("NOMAD Parameter Error:\n%s\n",e.what());
        return -1;        
    } 
    
    
    try
    {
    	pyEval *mSEval = new pyEval(p,cb,cbL,apply);
    
    	NOMAD::Mads mads (p, mSEval);
    
    	NOMAD::stop_type stopflag = mads.run();
    	
    	nb_evals = mads.get_stats().get_bb_eval();
    	nb_iters = mads.get_stats().get_iterations();
        	
    	// Set the best feasible solution or the best infeasible solution.
    	// One of the pointer is set to null to identify the type of solution
  	    const NOMAD::Eval_Point * bf = mads.get_best_feasible();
    	const NOMAD::Eval_Point *bimv = mads.get_best_infeasible_min_viol();
    	if ( bimv )
    	{
    	   	best_infeas_sol->set_f( bimv->get_f());
    	   	best_infeas_sol->set_h( bimv->get_h());
    	   	best_infeas_sol->set( *bimv ,static_cast<int>(nobj)); 
	        if ( !bf )
	        {
	        	delete best_feas_sol;
	        	best_feas_sol = NULL; 
	        }
    	}
    	if ( bf )
    	{
    		best_feas_sol->set_f(bf->get_f());
    	    best_feas_sol->set( *bf ,static_cast<int>(nobj));  	
    	    delete best_infeas_sol;
    	   	best_infeas_sol = NULL;
		}
    	return stopflag;
    	
    }
    catch(exception &e)
    {
        printf("NOMAD exception (report to developper):\n%s\n",e.what()); 
    }
    return -1; 
        
}
