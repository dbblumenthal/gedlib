#ifndef __MASTER_SLAVES__
#define __MASTER_SLAVES__

#include "Cache_Server.hpp"

using namespace NOMAD;
using namespace std;

// Cache server:
class Master_Slaves {
    
private:
    
    int                  _rank; // process rank
    int                    _np; // number of processes
    
    int                   _bbe; // max number of evaluations for each process
    int                    _ns; // number of free variables for each process
    
    Parameters            & _p; // parameters
    bool                _debug; // debug display flag
    
    
    static const int  TAG_SIGNAL;
    static const int  TAG_I1;
    static const int  TAG_I2;
    static const int  TAG_R1;
    static const int  TAG_D1;
    static const int  TAG_CSTOP;
    static       char STOP_SIGNAL;
    static       char OPTI_RES_SIGNAL;
    static       char OPTI_DATA_SIGNAL;
    
    // Receive an optimization result from the pollster:
    void receive_optimization_result ( int       & pollster_mesh_index ,
                                      bool      & stop_algo           ,
                                      double   *& best_feasible       ,
                                      double   *& best_infeasible     ,
                                      int         source ) const;
    
    // Send an optimization result to the master:
    void send_optimization_result ( int                pollster_mesh_index ,
                                   bool               stop_algo           ,
                                   const Eval_Point * bf                  ,
                                   const Eval_Point * bi                  ,
                                   stop_type          st ) const;
    
    // Send optimization data from the master to a slave:
    void send_optimization_data ( int            pollster_mesh_index ,
                                 bool           stop_algo           ,
                                 const double * best_feasible       ,
                                 const double * best_infeasible     ,
                                 int            source ) const;
    
    // Receive optimization data from the master:
    void receive_optimization_data ( bool   & stop_algo ,
                                    Point  & x0        ,
                                    Double & fx0 ) const;
    
    void receive_optimization_data ( bool   & stop_algo           ,
                                    Point  & x0                  ,
                                    Double & fx0                 ,
                                    int    & pollster_mesh_index ,
                                    int    * free_vars ) const;
    
    // Check the initial mesh size values:
    static bool check_delta ( const Point & delta );
    
public:
    
    // Constructor:
    Master_Slaves ( int       rank ,
                   int         np ,
                   int        bbe ,
                   int         ns ,
                   Parameters & p ,
                   bool     debug   )
    : _rank               ( rank                              ) ,
    _np                 ( np                                ) ,
    _bbe                ( bbe                               ) ,
    _ns                 ( ns                                ) ,
    _p                  ( p                                 ) ,
    _debug              ( debug                             ) {}
    
    // Destructor:
    virtual ~Master_Slaves ( void ) {}
    
    
    // Start the master:
    void start ( void ) const;
    
    // Stop the master:
    void stop ( void ) const;
    
    // MADS run:
    void mads_run ( Cache & cache , Cache & init_cache );
    
};


#endif
