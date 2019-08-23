#ifndef __CACHE_SERVER__
#define __CACHE_SERVER__

#include "mpi.h"
#include "nomad.hpp"

using namespace NOMAD;
using namespace std;

// Cache server:
class Cache_Server : public Cache
{
    
private:
    
    int                       _rank;             // process rank
    int                       _np;               // number of processes
    
    bool                      _debug;            // debug display flag
    
    Clock                     _clock;            // clock
    
    Double                    _h_min;            // h_min (min feasibility)
    int                       _max_bbe;          // max number of bb evaluations
    
    const Eval_Point        * _bf;               // best points
    const Eval_Point        * _bi1;
    const Eval_Point        * _bi2;
    
    Eval_Point              * _stop_point;       // stopping point
    
    mutable int               _multiple_evals;   // number of multiple evaluations
    mutable int               _cache_hits;       // number of cache hits
    mutable int               _cache_search_pts; // number of cache search points
    
    Point                  ** _waited_pts;      // list of points beeing evaluated
    list<const Eval_Point*> * _clients_ext_pts; // replaces _extern_pts
    
    string                    _history_file;
    
    // process the best feasible point signal:
    void process_bf_signal ( int source ) const;
    
    // process the extern point signal:
    void process_ep_signal ( int source ) const;
    
    // process the find signal:
    void process_find_signal ( int source ) const;
    
    // process the insertion signal:
    void process_insert_signal ( int source );
    
    // update and display the best points:
    void update_best_points ( const Eval_Point & x , int source );
    
public:
    
    static const int  TAG_SIGNAL;
    static const int  TAG_CACHE_HIT;
    static const int  TAG_X1;
    static const int  TAG_X2;
    static const int  TAG_X3;
    static const int  TAG_X4;
    static const int  TAG_X5;
    static const int  TAG_X6;
    static const int  TAG_X7;
    static const int  TAG_BBOR;
    static const int  TAG_BBOC;
    static const int  TAG_NB_EP;
    static const int  TAG_EP;
    static const int  TAG_BF;
    static       char STOP_SIGNAL;
    static       char FIND_SIGNAL;
    static       char INSERT_SIGNAL;
    static       char NB_EP_SIGNAL;
    static       char EP_SIGNAL;
    static       char BF_SIGNAL;
    
    // Constructor:
    Cache_Server ( const Display & out                  ,
                  int             rank                 ,
                  int             np                   ,
                  const Double  & h_min                ,
                  int             max_bbe              ,
                  const string  & history_file         ,
                  bool            allow_multiple_evals ,
                  bool            debug                  );
    
    // Destructor:
    virtual ~Cache_Server ( void );
    
    // Start the server:
    void start ( void );
    
    // Stop the server:
    void stop ( void ) const;
    
    // Display the clients extern points:
    void display_extern_pts ( void ) const { display_extern_pts(_out); }
    void display_extern_pts ( const Display & out ) const;
    
    // Display the best points:
    void display_best_points ( void ) const { display_best_points(_out); }
    void display_best_points ( const Display & out ) const;
    
    // Write the history file
    void write_his_file ( const NOMAD::Eval_Point & x ) const;
    
    // Display the current best solution:
    void display_current_solution ( void ) const;
    
    // Find a point:
    virtual const Eval_Point * find ( const Eval_Point & x ) const;
    
    // Insert a point:
    virtual void insert ( const NOMAD::Eval_Point & x );
    
    // Get the number of extern points:
    virtual int get_nb_extern_points ( void ) const;
    
    // Get and remove an extern point:
    virtual const Eval_Point * get_and_remove_extern_point ( void ) const;
    
};


#endif
