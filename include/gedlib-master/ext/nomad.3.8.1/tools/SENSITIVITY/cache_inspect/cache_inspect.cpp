/*------------------------------------------------------------------*/
/*                         Cache inspect                            */
/*------------------------------------------------------------------*/
/*  This program inspects a NOMAD cache and displays (cj,f) points. */
/*  Output is formatted as cj f 0/1 (1 for a dominating point)      */
/*  Usage:                                                          */
/*    cache_inspect cache_file                                      */
/*                  obj_index                                       */
/*                  cstr_index      (j)                             */
/*                  only_dominating (optional, 0 or 1)              */
/*------------------------------------------------------------------*/
#include "nomad.hpp"
#include <cstring>

using std::endl;

// prototypes:
bool process_cache ( const NOMAD::Display & ,
		     const NOMAD::Cache   & ,
		     int                    ,
		     int                    ,
		     const NOMAD::Point   & ,
		     const NOMAD::Point   & ,
		     bool                     );

bool check_args ( int                    ,
		  char                ** ,
		  const NOMAD::Display & ,
		  std::string          & ,
		  int                  & ,
		  int                  & ,
		  NOMAD::Point         & ,
		  NOMAD::Point         & ,
		  bool                 &   );

bool dominating ( const NOMAD::Eval_Point                  * ,
		  const std::multiset<NOMAD::Filter_Point> &   );

void display_usage ( const NOMAD::Display & , char ** );

/*-------------------------------------*/
/*            main function            */
/*-------------------------------------*/
int main ( int argc , char ** argv ) {

  // display:
  NOMAD::Display out ( std::cout );
  out.precision ( NOMAD::DISPLAY_PRECISION_STD );

  // check and interpret the arguments:
  std::string  cache_file;
  int          obj_index;
  int          cstr_index;
  NOMAD::Point lb , ub;
  bool         only_dominating;
  if ( !check_args ( argc            ,
		     argv            ,
		     out             ,
		     cache_file      ,
		     obj_index       ,
		     cstr_index      ,
		     lb              ,
		     ub              ,
		     only_dominating   ) )
    return EXIT_FAILURE;

  // read the cache file:
  NOMAD::Cache cache ( out );
  if ( !cache.load ( cache_file , NULL , false ) )
    return EXIT_FAILURE;

  // process the cache:
  if ( !process_cache ( out             ,
			cache           ,
			obj_index       ,
			cstr_index      ,
			lb              ,
			ub              ,
			only_dominating   ) )
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}

/*-------------------------------------*/
/*           process the cache         */
/*-------------------------------------*/
bool process_cache ( const NOMAD::Display & out             ,
		     const NOMAD::Cache   & cache           ,
		     int                    obj_index       ,
		     int                    cstr_index      ,
		     const NOMAD::Point   & lb              ,
		     const NOMAD::Point   & ub              ,
		     bool                   only_dominating   ) {

  // browse the cache:
  std::string                        error;
  std::multiset<NOMAD::Filter_Point> filter;
  bool                               feas;
  int                                i;
  int                                m   = -1;
  int                                n   = lb.size();
  const NOMAD::Eval_Point          * cur = cache.begin();

  if ( ub.size() != n )
    error = "problem with the bounds";

  // n may be equal to zero if bounds are not defined.

  while ( error.empty() && cur ) {

    // current cache point:
    NOMAD::Point bbo = cur->get_bb_outputs();

    // number of outputs:
    if ( m == -1 ) {
      m = bbo.size();
      if ( obj_index < 0 || obj_index >= m ) {
	error = "bad objective index";
	break;
      }

      if ( cstr_index >= m ) {

	if ( n == 0 ) {
	  error = "bounds not defined";
	  break;
	}

	int var_index = (cstr_index-m)/2;
	if ( var_index >= n ) {
	  error = "bad constraint index";
	  break;
	}

	if ( (cstr_index - m)%2 == 0 ) {
	  if ( !lb[var_index].is_defined() ) {
	    error = "lower bound undefined for constraint index";
	    break;
	  }
	}
	else if ( !ub[var_index].is_defined() ) {
	  error = "upper bound undefined for constraint index";
	  break;
	}
	
      }

      if ( cstr_index < 0 || cstr_index >= m+2*n ) {
	error = "bad constraint index";
	break;
      }
    }
    else if ( m != bbo.size() ) {
      error = "some entries in the cache have not the same number of outputs";
      break;
    }

    if ( n == 0 || cur->size() == n ) {

      // process bbo to include bound constraints:
      if ( n > 0 ) {
	bbo.resize ( m+2*n );
	for ( i = 0 ; i < n ; ++i ) {
	  if ( (*cur)[i].is_defined() ) {
	    if ( lb[i].is_defined() )
	      bbo[m+2*i] = lb[i] - (*cur)[i];
	    else
	      bbo[m+2*i] = 0.0;
	    if ( ub[i].is_defined() )
	      bbo[m+2*i+1] = (*cur)[i] - ub[i];
	    else
	      bbo[m+2*i+1] = 0.0;
	  }
	}
      }

      // consider only points satisfying all the
      // constraints except cj, with no undefined outputs:
      feas = ( cur->get_eval_status() == NOMAD::EVAL_OK );

      if ( feas )
	for ( i = 0 ; i < m+2*n ; ++i )
	  if ( !bbo[i].is_defined() || bbo[i] >= 1e+20 ||
	       ( i != obj_index && i != cstr_index && bbo[i] > 0.0 ) ) {
	    feas = false;
	    break;
	  }

      if ( feas ) {
	
	// set f and h:
	{
	  NOMAD::Eval_Point * x = const_cast<NOMAD::Eval_Point *> ( cur );
	  x->set_f ( bbo[obj_index ] );
	  x->set_h ( bbo[cstr_index] );
	}
	
	// insert in filter:
	{
	  NOMAD::Filter_Point tmp ( cur );
	  filter.insert ( tmp );
	}
      }
    }
    cur = cache.next();
  }
  
  if ( !error.empty() ) {
    std::cerr << endl << "ERROR: "
	      << error << endl << endl;
    return false;
  }

  // display the points:
  std::multiset<NOMAD::Filter_Point>::const_iterator
    it  = filter.begin() ,
    end = filter.end();
 
  while ( it != end ) {
    cur = (*it).get_point();

    if ( only_dominating ) {
      if ( dominating ( cur , filter ) )
	out << std::setw(17) << cur->get_h() << "\t"
	    << std::setw(17) << cur->get_f() << "\t"
	    << true << endl;
    }
    else
      out << std::setw(17) << cur->get_h() << "\t"
	  << std::setw(17) << cur->get_f() << "\t"
	  << dominating ( cur , filter )
	  << endl;

    ++it;
  }

  return true;
}

/*-------------------------------------*/
/*  check and interpret the arguments  */
/*-------------------------------------*/
bool check_args ( int                    argc            ,
		  char                ** argv            ,
		  const NOMAD::Display & out             ,
		  std::string          & cache_file      ,
		  int                  & obj_index       ,
		  int                  & cstr_index      ,
		  NOMAD::Point         & lb              ,
		  NOMAD::Point         & ub              ,
		  bool                 & only_dominating   ) {

  if ( argc < 4 || argc > 6 ) {
    display_usage ( out , argv );
    return false;
  }

  cache_file = std::string ( argv[1] );

  int i , n = strlen ( argv[2] );
  for ( i = 0 ; i < n ; ++i )
    if ( !isdigit ( argv[2][i] ) ) {
      display_usage ( out , argv );
      out << "(error with obj_index)" << endl << endl;
      return false;
    }
  obj_index = atoi ( argv[2] );
  if ( obj_index < 0 ) {
    display_usage ( out , argv );
    out << "(error with obj_index: < 0)" << endl << endl;
    return false;
  }

  n = strlen ( argv[3] );
  for ( i = 0 ; i < n ; ++i )
    if ( !isdigit ( argv[3][i] ) ) {
      display_usage ( out , argv );
      out << "(error with cstr_index)" << endl << endl;
      return false;
    }
  cstr_index = atoi ( argv[3] );
  if ( cstr_index < 0 ) {
    display_usage ( out , argv );
    out << "(error with cstr_index: < 0)" << endl << endl;
    return false;
  }

  // read the bounds:
  if ( argc > 4 ) {

    std::ifstream in ( argv[4] );

    if ( in.fail() ) {
      display_usage ( out , argv );
      out << "(error with bnds_file)" << endl << endl;
      return false;
    }

    std::string              s1 , s2;
    std::vector<std::string> ls;
    while ( !in.eof() ) {
      in >> s1 >> std::ws;
      if ( !in.eof() ) {
	in >> s2 >> std::ws;
	ls.push_back ( s1 );
	ls.push_back ( s2 );
      }
    }

    n = ls.size() / 2;
    lb.reset ( n );
    ub.reset ( n );
    for ( i = 0 ; i < n ; ++i ) {
      if ( !lb[i].atof ( ls[2*i  ] ) ||
	   !ub[i].atof ( ls[2*i+1] )    ) {
	display_usage ( out , argv );
	out << "(error with bnds_file)" << endl << endl;
	return false;
      }
    }

    in.close();
  }

  only_dominating = false;
  if ( argc > 5 )
    only_dominating = (atoi(argv[5]) == 1);

  return true;
}

/*----------------------------------*/
/*           display usage          */
/*----------------------------------*/
void display_usage ( const NOMAD::Display & out , char ** argv ) {
  out << endl << "usage: " << argv[0]
      << " cache_file obj_index cstr_index bounds_file only_dominating"
      << endl
      << "       cache_file     : binary cache file" << endl
      << "       obj_index      : objective index in [0,m-1]; m is the"
      << " number of outputs" << endl
      << "       cstr_index     : index of output j in [0,m-1]" << endl
      << "                        j=m+2i corresponds to the lower"
      << " bound of variable i" << endl
      << "                        j=m+2i+1 corresponds to the upper"
      << " bound of variable i" << endl
      << "       bounds_file    : optional: one line in this file"
      << " corresponds to a lower and to an upper bound" << endl
      << "       only_dominating: optional: set to 1 to display"
      << " only dominating points, 0 otherwise" << endl
      << endl;
}

/*----------------------------------*/
/*  check if a point is dominating  */
/*----------------------------------*/
bool dominating ( const NOMAD::Eval_Point                  * x      ,
		  const std::multiset<NOMAD::Filter_Point> & filter   ) {
  const NOMAD::Eval_Point * y;
  std::set<NOMAD::Filter_Point>::const_iterator
    it  = filter.begin() ,
    end = filter.end();
  while ( it != end ) {
    y = (*it).get_point();
    if ( (*y) < (*x) )
      return false; 
    ++it;
  }
  return true;
}
