#include "Master_Slaves.hpp"
using namespace std;
using namespace NOMAD;

/*--------------------------*/
/*  tags and signal values  */
/*--------------------------*/
const int  Master_Slaves::TAG_SIGNAL       = 100;
const int  Master_Slaves::TAG_I1           = 101;
const int  Master_Slaves::TAG_R1           = 102;
const int  Master_Slaves::TAG_CSTOP        = 103;
const int  Master_Slaves::TAG_D1           = 105;
const int  Master_Slaves::TAG_I2           = 106;
char       Master_Slaves::STOP_SIGNAL      = 'S';
char       Master_Slaves::OPTI_RES_SIGNAL  = 'O';
char       Master_Slaves::OPTI_DATA_SIGNAL = 'D';

/*--------------------------------*/
/*  start the master (process 0)  */
/*--------------------------------*/
void Master_Slaves::start ( void ) const {
    
    if ( _rank != 0 )
        return;
    
    MPI_Status status;
    int        source;
    char       signal;
    int        n                   = _p.get_dimension();
    int        nb_stops            = 0;
    int        pollster_mesh_index = 0;
    bool       stop_algo           = false;
    double   * best_feasible       = NULL;
    
    // the first best infeasible point is initialized with x0:
    double      * best_infeasible = new double [n+2];
    const Point * x0              = (_p.get_x0s())[0];
    for ( int i = 0 ; i < n ; ++i )
        best_infeasible[i] = (*x0)[i].value();
    best_infeasible[n  ] = INF;
    best_infeasible[n+1] = INF;
    
    /*-------------*/
    /*  main loop  */
    /*-------------*/
    while ( nb_stops != _np-2 ) {
        
        MPI_Recv ( &signal , 1 , MPI_CHAR , MPI_ANY_SOURCE ,
                  Master_Slaves::TAG_SIGNAL , MPI_COMM_WORLD , &status );
        
        source = status.MPI_SOURCE;
        
        // stop signal:
        // ------------
        if ( signal == Master_Slaves::STOP_SIGNAL ) {
            if ( _debug ) {
                _p.out() << "MASTER: STOP SIGNAL FROM RANK " << source;
                if ( source == 1 )
                    _p.out() << " (POLLSTER)";
                _p.out() << endl;
            }
            ++nb_stops;
        }
        
        // optimization data signal:
        // -------------------------
        else if ( signal == Master_Slaves::OPTI_DATA_SIGNAL ) {
            if ( _debug ) {
                _p.out() << "MASTER: receive optimization data request from slave "
                << source;
                if ( source == 1 )
                    _p.out() << " (POLLSTER)";
                _p.out() << endl;
            }
            
            send_optimization_data ( pollster_mesh_index ,
                                    stop_algo           ,
                                    best_feasible       ,
                                    best_infeasible     ,
                                    source                );
        }
        
        // optimization result signal:
        // ---------------------------
        else if ( signal == Master_Slaves::OPTI_RES_SIGNAL ) {
            if ( _debug ) {
                _p.out() << "MASTER: receive optimization result from slave "
                << source;
                if ( source == 1 )
                    _p.out() << " (POLLSTER)";
                _p.out() << endl;
            }
            
            // pollster:
            if ( source == 1 )
                receive_optimization_result ( pollster_mesh_index ,
                                             stop_algo           ,
                                             best_feasible       ,
                                             best_infeasible     ,
                                             source                );
            
            // other slaves:
            else {
                int  tmp1;
                bool tmp2;
                receive_optimization_result ( tmp1            ,
                                             tmp2            ,
                                             best_feasible   ,
                                             best_infeasible ,
                                             source            );
            }
        }
    }
    
    if ( best_feasible )
        delete [] best_feasible;
    if ( best_infeasible )
        delete [] best_infeasible;
}

/*------------------------------------------*/
/*  stop the master (processes 1 to _np-2)  */
/*------------------------------------------*/
void Master_Slaves::stop ( void ) const {
    if ( _rank == 0 || _rank == _np-1 )
        return;
    
    MPI_Send ( &Master_Slaves::STOP_SIGNAL , 1 , MPI_CHAR ,
              0 , Master_Slaves::TAG_SIGNAL , MPI_COMM_WORLD );
}

/*----------------------------------------*/
/*                MADS run                */
/*----------------------------------------*/
void Master_Slaves::mads_run ( Cache & cache ) {
    
    const Eval_Point * best_feasible   = NULL;
    const Eval_Point * best_infeasible = NULL;
    Double             old_f           = INF;
    bool               stop_algo       = false;
    int                run_index       = 0;
    int                mesh_index      = 0;
    double             default_eps     = Double::get_epsilon();
    int                n               = _p.get_dimension();
    int              * free_vars       = new int [_ns];
    Point              x0 ( n ) , delta_0 , delta_min;
    
    /*---------------------*/
    /*      main loop      */
    /*---------------------*/
    while ( !stop_algo ) {
        
        best_feasible   = NULL;
        best_infeasible = NULL;
        
        // Seed:
        _p.set_SEED ( 99 * _rank + 20 * run_index );
        
        // first run:
        if ( run_index == 0 ) {
            
            // max number of evaluations for regular slaves:
            if ( _rank != 1 )
                _p.set_MAX_BB_EVAL ( _bbe );
            
            // display:
            // p.set_DISPLAY_STATS ( "process #" + itos(_rank) + " BBE OBJ" );
            // p.set_DISPLAY_DEGREE ( FULL_DISPLAY );
            _p.set_DISPLAY_DEGREE ( NO_DISPLAY );
        }
        
        // force the parameters check:
        _p.force_check_flag();
        
        /*------------------*/
        /*  pollster slave  */
        /*------------------*/
        if ( _rank == 1 ) {
            
            stop_type stop_reason = UNKNOWN_STOP_REASON;
            
            _p.check();   // must do check to get a valid signature
            _p.get_signature()->get_mesh()->set_mesh_indices( NOMAD::Point( n,mesh_index ) );
            delta_0=_p.get_signature()->get_mesh()->get_delta ( );
            
            
            
            // we set a very small epsilon in order to accept
            // very small initial mesh sizes:
            Double::set_epsilon ( 1e-16 );
            
            if ( !check_delta ( delta_0 ) )
                stop_algo = true;
            
            else {
                
                // first run:
                if ( run_index == 0 ) {
                    
                    // directions:
                    {
                        bool use_orthomads = _p.has_orthomads_directions();
                        _p.reset_directions   ( );
                        _p.set_DIRECTION_TYPE ( (use_orthomads) ? ORTHO_1 : LT_1 );
                    }
                    
                    // cache search:
                    _p.set_CACHE_SEARCH               ( true  );
                    _p.set_OPPORTUNISTIC_CACHE_SEARCH ( false );
                }
                
                // other runs:
                else {
                    
                    // stop_algo may be set to 'false' here:
                    receive_optimization_data ( stop_algo , x0 , old_f );
                    
                    // starting point:
                    _p.reset_X0();
                    _p.set_X0 ( x0 );
                }
                
                if ( !stop_algo ) {
                    
                    // check the parameters:
                    _p.check();
                    
                    _p.get_signature()->get_mesh()->set_min_mesh_sizes( delta_0 );
                    
                    if ( mesh_index <= 0 )
                        _p.get_signature()->get_mesh()->set_limit_mesh_index ( mesh_index );
                    else
                        _p.get_signature()->get_mesh()->set_limit_mesh_index ( 0 );
                    
                    
                    Double::set_epsilon ( default_eps );
                    
                    Mads mads ( _p , &_ev , NULL , &cache , NULL );
                    stop_reason     = mads.run();
                    best_feasible   = mads.get_best_feasible();
                    best_infeasible = mads.get_best_infeasible();
                    
                    bool success = false;
                    
                    if ( best_feasible ) {
                        
                        success = (best_feasible->get_f() < old_f);
                        
                        if ( _debug )
                            _p.out() << "POLLSTER: ELL="
                            << mesh_index << " BBE="
                            << mads.get_stats().get_bb_eval()
                            << " OLD_F=" << old_f << " NEW_F="
                            << best_feasible->get_f()
                            << " SUCCESS=" << success << endl;
                    }
                    
                    // pollster mesh update:
                    if ( success )
                        ++mesh_index;
                    else
                        --mesh_index;
                    
                }
            }
            
            send_optimization_result ( mesh_index      ,
                                      stop_algo       ,
                                      best_feasible   ,
                                      best_infeasible ,
                                      stop_reason       );
        }
        
        /*------------------*/
        /*  regular slaves  */
        /*------------------*/
        else {
            
            int i , j , pollster_mesh_index;
            
            receive_optimization_data ( stop_algo           ,
                                       x0                  ,
                                       old_f               ,
                                       pollster_mesh_index ,
                                       free_vars             );
            
            if ( _debug ) {
                _p.out() << "SLAVE #" << _rank
                << ": OPTIM. DATA: [STOP=" << stop_algo
                << "] [POLLSTER_MESH_INDEX=" << pollster_mesh_index
                << "] [X0=" << x0 << " ] [f(X0)="
                << old_f << "] [FREE VARS= ";
                for ( i = 0 ; i < _ns ; ++i )
                    _p.out() << free_vars[i] << " ";
                _p.out() << " ]" << endl;
            }
            
            if ( !stop_algo ) {
                
                // starting point:
                _p.reset_X0();
                _p.set_X0 ( x0 );
                
                // mesh of the regular slave:
                int ell_0 = 0;
                if ( pollster_mesh_index > ell_0 )
                    ell_0 = pollster_mesh_index;
                
                _p.check();  // Must do check to access signature
                _p.get_signature()->get_mesh()->set_mesh_indices( NOMAD::Point( n,ell_0 ) );
                delta_0=_p.get_signature()->get_mesh()->get_delta();
                
                _p.get_signature()->get_mesh()->set_mesh_indices( NOMAD::Point( n ,pollster_mesh_index ) );
                delta_min=_p.get_signature()->get_mesh()->get_delta();
                
                
                
                Double::set_epsilon ( 1e-16 );
                if ( !check_delta ( delta_0 ) || !check_delta ( delta_min ) )
                    stop_algo = true;
                
                else {
                    
                    // free variables:
                    {
                        _p.reset_fixed_variables();
                        bool fix_var;
                        for ( i = 0 ; i < n ; ++i ) {
                            fix_var = true;
                            for ( j = 0 ; j < _ns ; ++j )
                                if ( free_vars[j] == i ) {
                                    fix_var = false;
                                    break;
                                }
                            if ( fix_var )
                                _p.set_FIXED_VARIABLE ( i );
                        }
                    }
                    
                    // check the parameters:
                    _p.check();
                    
                    // modify mesh termination criterions
                    _p.get_signature()->get_mesh()->set_mesh_indices( NOMAD::Point( n,ell_0 ) );
                    
                    if ( pollster_mesh_index <=0 )
                        _p.get_signature()->get_mesh()->set_limit_mesh_index( pollster_mesh_index );
                    else
                        _p.get_signature()->get_mesh()->set_limit_mesh_index( 0 );
                    
                    _p.get_signature()->get_mesh()->set_min_mesh_sizes( delta_min );
                    _p.get_signature()->get_mesh()->set_delta_0 ( delta_0 );
                    
                    
                    
                    Double::set_epsilon ( default_eps );
                    
                    // MADS run:
                    Mads mads ( _p , &_ev , NULL , &cache , NULL );
                    mads.run();
                    best_feasible   = mads.get_best_feasible();
                    best_infeasible = mads.get_best_infeasible();
                    
                    if ( _debug && best_feasible ) {
                        _p.out() << "RANK #" << _rank << ": POLLSTER_ELL="
                        << pollster_mesh_index << " VARS = [";
                        for ( i = 0 ; i < _ns ; ++i )
                            _p.out() << free_vars[i] << " ";
                        _p.out() << " ] BBE=" << mads.get_stats().get_bb_eval()
                        << " OLD_F=" << old_f << " NEW_F="
                        << best_feasible->get_f()
                        << " SUCCESS="
                        << (best_feasible->get_f() < old_f)
                        << endl;
                    }
                }
            }
            
            {
                int       tmp1 = -1;
                bool      tmp2 = false;
                stop_type tmp3 = UNKNOWN_STOP_REASON;
                send_optimization_result ( tmp1            ,
                                          tmp2            ,
                                          best_feasible   ,
                                          best_infeasible ,
                                          tmp3              );
            }
        }
        
        // loop increment:
        ++run_index;
    }
    
    delete [] free_vars;
}

/*----------------------------------------------------*/
/*  receive an optimization result from the pollster  */
/*  POLLSTER --> MASTER                               */
/*----------------------------------------------------*/
void Master_Slaves::receive_optimization_result
( int     & pollster_mesh_index ,
 bool    & stop_algo           ,
 double *& best_feasible       ,
 double *& best_infeasible     ,
 int       source                ) const {
    
    int        itab[5];
    MPI_Status status;
    
    MPI_Recv ( itab , 5 , MPI_INT , source ,
              Master_Slaves::TAG_I1 , MPI_COMM_WORLD , &status );
    
    pollster_mesh_index = itab[0];
    
    // stop the algorithm ?
    stop_algo = ( itab[4] == 1 );
    
    if ( !stop_algo ) {
        
        stop_type stop_reason = static_cast<stop_type>(itab[1]);
        
        switch ( stop_reason ) {
            case ERROR:
            case UNKNOWN_STOP_REASON:
            case CTRL_C:
            case MESH_PREC_REACHED:
            case X0_FAIL:
            case P1_FAIL:
            case L_MAX_REACHED:
            case L_LIMITS_REACHED:
            case MAX_TIME_REACHED:
            case MAX_BB_EVAL_REACHED:
            case MAX_SGTE_EVAL_REACHED:
            case MAX_EVAL_REACHED:
            case MAX_SIM_BB_EVAL_REACHED:
            case MAX_ITER_REACHED:
            case FEAS_REACHED:
            case F_TARGET_REACHED:
            case STAT_SUM_TARGET_REACHED:
            case L_CURVE_TARGET_REACHED:
            case MULTI_MAX_BB_REACHED:
            case MULTI_NB_MADS_RUNS_REACHED:
            case MULTI_STAGNATION:
            case MULTI_NO_PARETO_PTS:
            case MAX_CACHE_MEMORY_REACHED:
                stop_algo = true;
            default:
                stop_algo = false;
        }
    }
    
    int i , nb_pts = 0;
    if ( itab[2] == 1 )
        ++nb_pts;
    if ( itab[3] == 1 )
        ++nb_pts;
    
    // itab[2] == 1 --> bf != NULL
    // itab[3] == 1 --> bi != NULL
    
    if ( nb_pts > 0 ) {
        
        int      n    = _p.get_dimension();
        double * rtab = new double [(n+2)*nb_pts];
        
        MPI_Recv ( rtab , (n+2)*nb_pts , MPI_DOUBLE , source ,
                  Master_Slaves::TAG_R1 , MPI_COMM_WORLD , &status );
        
        if ( nb_pts == 2 ) {
            
            // best feasible and infeasible updates:
            bool update = false;
            
            if ( best_feasible ) {
                Double old_f = best_feasible[n+1];
                Double new_f = rtab         [n+1];
                if ( new_f < old_f )
                    update = true;
            }
            else {
                best_feasible = new double[n+2];
                update = true;
            }
            
            if ( update ) {
                for ( i = 0 ; i < n ; ++i )
                    best_feasible[i] = rtab[i];
                best_feasible[n  ] = rtab[n];
                best_feasible[n+1] = rtab[n+1];
            }
            
            update = false;
            
            if ( best_infeasible ) {
                
                Double old_h = best_infeasible[n];
                Double new_h = rtab           [2*n+2];
                if ( new_h < old_h )
                    update = true;
            }
            else {
                best_infeasible = new double[n+2];
                update = true;
            }
            
            if ( update ) {
                int cur = n+2;
                for ( i = 0 ; i < n ; ++i )
                    best_infeasible[i] = rtab[cur++];
                best_infeasible[n  ] = rtab[cur++];
                best_infeasible[n+1] = rtab[cur];
            }
            delete [] rtab;
        }
        else {
            
            // best feasible update:
            if ( itab[2] == 1 ) {
                if ( best_feasible ) {
                    Double old_f = best_feasible[n+1];
                    Double new_f = rtab         [n+1];
                    if ( new_f < old_f ) {
                        delete [] best_feasible;
                        best_feasible = rtab;
                    }
                    else
                        delete [] rtab;
                }
                else
                    best_feasible = rtab;
            }
            
            // best infeasible update:
            else {
                if ( best_infeasible ) {
                    Double old_h = best_infeasible[n];
                    Double new_h = rtab           [n];
                    if ( new_h < old_h ) {
                        delete [] best_infeasible;
                        best_infeasible = rtab;
                    }
                    else
                        delete [] rtab;
                }
                else
                    best_infeasible = rtab;
            }
        }
    }
}

/*---------------------------------------------*/
/*  send an optimization result to the master  */
/*  POLLSTER --> MASTER                        */
/*---------------------------------------------*/
void Master_Slaves::send_optimization_result
( int                pollster_mesh_index ,
 bool               stop_algo           ,
 const Eval_Point * bf                  ,
 const Eval_Point * bi                  ,
 stop_type          st                    ) const {
    
    // send a signal to the master:
    MPI_Send ( &Master_Slaves::OPTI_RES_SIGNAL , 1 , MPI_CHAR ,
              0 , Master_Slaves::TAG_SIGNAL , MPI_COMM_WORLD );
    
    // send the data:
    int itab[5];
    
    itab[0] = pollster_mesh_index;
    itab[1] = static_cast<int>(st);
    
    int nb_pts = 0;
    if ( bf ) {
        ++nb_pts;
        itab[2] = 1;
    }
    else
        itab[2] = 0;
    
    if ( bi ) {
        ++nb_pts;
        itab[3] = 1;
    }
    else
        itab[3] = 0;
    
    itab[4] = (stop_algo) ? 1 : 0;
    
    MPI_Send ( itab , 5 , MPI_INT , 0 ,
              Master_Slaves::TAG_I1 , MPI_COMM_WORLD );
    
    if ( nb_pts > 0 ) {
        
        int      n    = _p.get_dimension();
        double * rtab = new double [(n+2)*nb_pts];
        
        int i , cur = 0;
        if ( bf ) {
            for ( i = 0 ; i < n ; ++i )
                rtab[cur++] = (*bf)[i].value();
            rtab[cur++] = bf->get_h().value();
            rtab[cur++] = bf->get_f().value();
        }
        
        if ( bi ) {
            for ( i = 0 ; i < n ; ++i )
                rtab[cur++] = (*bi)[i].value();
            rtab[cur++] = bi->get_h().value();
            rtab[cur  ] = bi->get_f().value();
        }
        
        MPI_Send ( rtab , cur , MPI_DOUBLE , 0 ,
                  Master_Slaves::TAG_R1 , MPI_COMM_WORLD );
        
        delete [] rtab;
    }
}

/*---------------------------------------------*/
/*  receive optimization data from the master  */
/*  MASTER --> REGULAR SLAVE                   */
/*---------------------------------------------*/
void Master_Slaves::receive_optimization_data
( bool   & stop_algo           ,
 Point  & x0                  ,
 Double & fx0                 ,
 int    & pollster_mesh_index ,
 int    * free_vars             ) const {
    
    // step 1/2: receive common pollster data:
    receive_optimization_data ( stop_algo , x0 , fx0 );
    
    int i;
    
    // step 2/2: receive additional data for regular slaves:
    if ( !stop_algo ) {
        
        int      * itab = new int [_ns+1];
        MPI_Status status;
        
        MPI_Recv ( itab , _ns+1 , MPI_INT , 0 ,
                  Master_Slaves::TAG_I2 , MPI_COMM_WORLD , &status );
        
        for ( i = 0 ; i < _ns ; ++i )
            free_vars[i] = itab[i];
        
        pollster_mesh_index = itab[_ns];
        
        delete [] itab;
    }
    else {
        pollster_mesh_index = 1;
        for ( i = 0 ; i < _ns ; ++i )
            free_vars[i] = -1;
    }
}

/*---------------------------------------------*/
/*  receive optimization data from the master  */
/*  MASTER --> POLLSTER                        */
/*---------------------------------------------*/
void Master_Slaves::receive_optimization_data ( bool   & stop_algo ,
                                               Point  & x0        ,
                                               Double & fx0         ) const {
    char        c_stop;
    MPI_Status  status;
    MPI_Request req    = MPI_REQUEST_NULL;
    int         i , n  = _p.get_dimension();
    
    // send a request for c_stop:
    MPI_Irecv ( &c_stop , 1 , MPI_CHAR , 0 ,
               Master_Slaves::TAG_CSTOP , MPI_COMM_WORLD , &req );
    
    // send a signal to the master:
    MPI_Send ( &Master_Slaves::OPTI_DATA_SIGNAL , 1 , MPI_CHAR ,
              0 , Master_Slaves::TAG_SIGNAL , MPI_COMM_WORLD );
    
    MPI_Wait ( &req , &status );
    
    // stop:
    if ( c_stop == '1' ) {
        stop_algo = true;
        x0.reset();
    }
    
    // continue:
    else {
        
        stop_algo     = false;
        double * rtab = new double [n+2];
        
        MPI_Recv ( rtab , n+2 , MPI_DOUBLE , 0 ,
                  Master_Slaves::TAG_D1 , MPI_COMM_WORLD , &status );
        
        for ( i = 0 ; i < n ; ++i )
            x0[i] = rtab[i];
        
        Double h = rtab[n];
        fx0 = ( h <= _p.get_h_min() ) ? rtab[n+1] : INF;
        
        // get the best feasible point from the cache server:
        {
            char pt_flag;
            int  npm1 = _np-1;
            MPI_Irecv ( &pt_flag , 1 , MPI_CHAR , npm1 ,
                       Cache_Server::TAG_BF , MPI_COMM_WORLD , &req );
            
            MPI_Send ( &Cache_Server::BF_SIGNAL , 1 , MPI_CHAR ,
                      npm1 , Cache_Server::TAG_SIGNAL , MPI_COMM_WORLD );
            
            MPI_Wait ( &req , &status );
            
            if ( pt_flag == '1' ) {
                
                MPI_Recv ( rtab , n+2 , MPI_DOUBLE , npm1 ,
                          Cache_Server::TAG_X7 , MPI_COMM_WORLD , &status );
                
                for ( i = 0 ; i < n ; ++i )
                    x0[i] = rtab[i];
                
                Double h = rtab[n];
                fx0 = ( h <= _p.get_h_min() ) ? rtab[n+1] : INF;
                
            }
        }
        delete [] rtab;
        
        // check the bounds:
        const Point & lb = _p.get_lb();
        const Point & ub = _p.get_ub();
        
        for ( i = 0 ; i < n ; ++i ) {
            if ( lb[i].is_defined() && x0[i].value() < lb[i].value() )
                x0[i] = lb[i];
            if ( ub[i].is_defined() && x0[i].value() > ub[i].value() )
                x0[i] = ub[i];
        }
    }
}

/*-----------------------------------------------------*/
/*  send optimization data from the master to a slave  */
/*  MASTER --> POLLSTER or MASTER --> SLAVE            */
/*-----------------------------------------------------*/
void Master_Slaves::send_optimization_data
( int            pollster_mesh_index ,
 bool           stop_algo           ,
 const double * best_feasible       ,
 const double * best_infeasible     ,
 int            source                ) const {
    
    char c_stop = (stop_algo) ? '1' : '0';
    
    MPI_Rsend ( &c_stop , 1 , MPI_CHAR , source ,
               Master_Slaves::TAG_CSTOP , MPI_COMM_WORLD );
    
    // continue:
    if ( !stop_algo ) {
        
        int n = _p.get_dimension();
        
        // data for pollster and regular slaves:
        if ( best_feasible )
            MPI_Send ( const_cast<double*>(best_feasible) , n+2 , MPI_DOUBLE ,
                      source , Master_Slaves::TAG_D1 , MPI_COMM_WORLD );
        else
            MPI_Send ( const_cast<double*>(best_infeasible) , n+2 , MPI_DOUBLE ,
                      source , Master_Slaves::TAG_D1 , MPI_COMM_WORLD );
        
        // additional data for regular slaves:
        if ( source != 1 ) {
            
            int * itab = new int [_ns+1];
            
            // choose the free varables:
            {
                Random_Pickup rp ( n );
                
                for ( int i = 0 ; i < _ns ; ++i )
                    itab[i] = rp.pickup(); // index of the ith free variable
            }
            
            itab[_ns] = pollster_mesh_index;
            
            MPI_Send ( itab , _ns+1 , MPI_INT , source ,
                      Master_Slaves::TAG_I2 , MPI_COMM_WORLD );
            
            delete [] itab;
        }
    }
}

/*-----------------------------------------------*/
/*        check the initial mesh size values     */
/*-----------------------------------------------*/
bool Master_Slaves::check_delta ( const Point & delta ) {
    int n = delta.size();
    for ( int i = 0 ; i < n ; ++i )
        if ( delta[i].value() < Double::get_epsilon() ||
            delta[i].value() <= 0.0 )
            return false;
    return true;
}

