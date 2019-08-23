#include "Cache_Server.hpp"
using namespace std;
using namespace NOMAD;

/*--------------------------*/
/*  tags and signal values  */
/*--------------------------*/
const int  Cache_Server::TAG_SIGNAL    = 0;
const int  Cache_Server::TAG_CACHE_HIT = 1;
const int  Cache_Server::TAG_X1        = 2;
const int  Cache_Server::TAG_X2        = 3;
const int  Cache_Server::TAG_X3        = 4;
const int  Cache_Server::TAG_X4        = 5;
const int  Cache_Server::TAG_X5        = 6;
const int  Cache_Server::TAG_X6        = 7;
const int  Cache_Server::TAG_X7        = 8;
const int  Cache_Server::TAG_BBOR      = 9;
const int  Cache_Server::TAG_BBOC      = 10;
const int  Cache_Server::TAG_NB_EP     = 11;
const int  Cache_Server::TAG_EP        = 12;
const int  Cache_Server::TAG_BF        = 13;
char       Cache_Server::STOP_SIGNAL   = 'S';
char       Cache_Server::FIND_SIGNAL   = 'F';
char       Cache_Server::INSERT_SIGNAL = 'I';
char       Cache_Server::NB_EP_SIGNAL  = 'N';
char       Cache_Server::EP_SIGNAL     = 'E';
char       Cache_Server::BF_SIGNAL     = 'B';

/*-----------------------------------*/
/*            constructor            */
/*-----------------------------------*/
Cache_Server::Cache_Server ( const Display & out                  ,
                            int             rank                 ,
                            int             np                   ,
                            const Double  & h_min                ,
                            int             max_bbe              ,
                            const string  & history_file         ,
                            bool            allow_multiple_evals ,
                            bool            debug                  )
: Cache             ( out , TRUTH ) ,
_rank             ( rank        ) ,
_np               ( np          ) ,
_debug            ( debug       ) ,
_h_min            ( h_min       ) ,
_max_bbe          ( max_bbe     ) ,
_history_file     ( history_file) ,
_bf               ( NULL        ) ,
_bi1              ( NULL        ) ,
_bi2              ( NULL        ) ,
_multiple_evals   ( 0           ) ,
_cache_hits       ( 0           ) ,
_cache_search_pts ( 0           ) ,
_waited_pts       ( NULL        ) ,
_clients_ext_pts  ( NULL        )
{
    
    // cache server:
    if ( _rank == _np - 1 )
    {
        
        _clients_ext_pts = new list<const Eval_Point*> [_np];
        
        if ( !allow_multiple_evals )
        {
            _waited_pts = new Point * [_np];
            for ( int i = 0 ; i < _np ; ++i )
                _waited_pts[i] = NULL;
        }
        
    }
}

/*-----------------------------------*/
/*             destructor            */
/*-----------------------------------*/
Cache_Server::~Cache_Server ( void )
{
    if ( _waited_pts )
    {
        for ( int i = 0 ; i < _np ; ++i )
            if ( _waited_pts[i] )
                delete _waited_pts;
        delete [] _waited_pts;
    }
    
    if ( _clients_ext_pts )
        delete [] _clients_ext_pts;
}

/*-----------------------------------*/
/*  start the server (process np-1)  */
/*-----------------------------------*/
void Cache_Server::start ( void )
{
    
    int npm1 = _np-1;
    
    if ( _rank != npm1 )
        return;
    
    MPI_Status status;
    int        nb_stops = 0;
    int        source;
    char       signal;
    
    /*-------------*/
    /*  main loop  */
    /*-------------*/
    while ( nb_stops != npm1 )
    {
        
        MPI_Recv ( &signal , 1 , MPI_CHAR , MPI_ANY_SOURCE ,
                  Cache_Server::TAG_SIGNAL , MPI_COMM_WORLD , &status );
        
        source = status.MPI_SOURCE;
        
        // stop signal:
        // ------------
        if ( signal == Cache_Server::STOP_SIGNAL )
        {
            if ( _debug )
                _out << "CACHE SERVER: STOP SIGNAL FROM RANK " << source << endl;
            ++nb_stops;
        }
        
        // find signal:
        // ------------
        else if ( signal == Cache_Server::FIND_SIGNAL )
        {
            if ( _debug )
                _out << "CACHE SERVER: FIND SIGNAL FROM RANK " << source << endl;
            process_find_signal ( source );
        }
        
        // insert signal:
        // --------------
        else if ( signal == Cache_Server::INSERT_SIGNAL )
        {
            if ( _debug )
            {
                _out << "CACHE SERVER: INSERT SIGNAL FROM RANK "
                << source;
                if ( source == 1 )
                    _out << " (POLLSTER)";
                _out << endl;
            }
            process_insert_signal ( source );
        }
        
        // number of extern points signal:
        // -------------------------------
        else if ( signal == Cache_Server::NB_EP_SIGNAL )
        {
            if ( _debug )
            {
                _out << "CACHE SERVER: NB EXTERN POINTS SIGNAL FROM RANK "
                << source;
                if ( source == 1 )
                    _out << " (POLLSTER)";
                _out << endl;
            }
            int nb_client_extern_pts = _clients_ext_pts[source].size();
            MPI_Rsend ( &nb_client_extern_pts , 1 , MPI_INT , source ,
                       Cache_Server::TAG_NB_EP , MPI_COMM_WORLD );
        }
        
        // extern point signal:
        // --------------------
        else if ( signal == Cache_Server::EP_SIGNAL )
        {
            if ( _debug )
            {
                _out << "CACHE SERVER: EXTERN POINT SIGNAL FROM RANK "
                << source;
                if ( source == 1 )
                    _out << " (POLLSTER)";
                _out << endl;
            }
            process_ep_signal ( source );
        }
        
        // best feasible point signal:
        else if ( signal == Cache_Server::BF_SIGNAL )
        {
            if ( _debug )
            {
                _out << "CACHE SERVER: BEST FEASIBLE POINT SIGNAL FROM RANK "
                << source;
                if ( source == 1 )
                    _out << " (POLLSTER)";
            }
            process_bf_signal ( source );
        }
    }
}

/*---------------------------------*/
/*    stop the server (clients)    */
/*---------------------------------*/
void Cache_Server::stop ( void ) const
{
    
    int npm1 = _np-1;
    
    if ( _rank == npm1 )
        return;
    
    MPI_Send ( &Cache_Server::STOP_SIGNAL , 1 , MPI_CHAR ,
              npm1 , Cache_Server::TAG_SIGNAL , MPI_COMM_WORLD );
}

/*----------------------------------------------*/
/*    process the best feasible point signal    */
/*----------------------------------------------*/
void Cache_Server::process_bf_signal ( int source ) const
{
    
    char pt_flag = (_bf) ? '1' : '0';
    
    MPI_Rsend ( &pt_flag , 1 , MPI_CHAR , source ,
               Cache_Server::TAG_BF , MPI_COMM_WORLD );
    
    if ( _bf )
    {
        int n = _bf->size();
        double * rtab = new double[n+2];
        for ( int i = 0 ; i < n ; ++i )
            rtab[i] = (*_bf)[i].value();
        rtab[n] = _bf->get_h().value();
        rtab[n+1] = _bf->get_f().value();
        
        MPI_Send ( rtab , n+2 , MPI_DOUBLE , source ,
                  Cache_Server::TAG_X7 , MPI_COMM_WORLD );
        
        delete [] rtab;
    }
    
}

/*---------------------------------------*/
/*    process the extern point signal    */
/*---------------------------------------*/
void Cache_Server::process_ep_signal ( int source ) const
{
    
    int nb_pt = ( _clients_ext_pts[source].size() > 0 ) ? 1 : 0;
    
    MPI_Rsend ( &nb_pt , 1 , MPI_INT , source ,
               Cache_Server::TAG_EP , MPI_COMM_WORLD );
    
    // send and remove the extern point:
    if ( nb_pt > 0 )
    {
        
        const Eval_Point * x = *(_clients_ext_pts[source].begin());
        
        ++_cache_search_pts;
        
        // send the point :
        int i , n = x->size() , m = x->get_m();
        int itab[5];
        itab[0] = n;
        itab[1] = m;
        itab[2] = ( x->is_eval_ok() ) ? 1 : 0;
        
        if ( x->get_signature() && (x->get_signature()->get_mesh()->get_mesh_indices())[0].is_defined() )
        {
            itab[3]=1;
            itab[4]=static_cast<int>((x->get_signature()->get_mesh()->get_mesh_indices())[0].value());
        }
        else
            itab[3] = itab[4] = 0;
        
        double * rtab = new double[n+2*m];
        for ( i = 0 ; i < n ; ++i )
            rtab[i] = (*x)[i].value();
        
        const Point & bbo = x->get_bb_outputs();
        
        for ( i = 0 ; i < m ; ++i )
            if ( bbo[i].is_defined() )
            {
                rtab[2*i+n  ] = 1.0;
                rtab[2*i+n+1] = bbo[i].value();
            }
            else
                rtab[2*i+n] = rtab[2*i+n+1] = -1.0;
        
        MPI_Send ( itab , 5 , MPI_INT , source ,
                  Cache_Server::TAG_X5 , MPI_COMM_WORLD );
        
        MPI_Send ( rtab , n+2*m , MPI_DOUBLE , source ,
                  Cache_Server::TAG_X6 , MPI_COMM_WORLD );
        
        // remove the point:
        _clients_ext_pts[source].pop_front();
    }
}

/*-------------------------------*/
/*    process the find signal    */
/*-------------------------------*/
void Cache_Server::process_find_signal ( int source ) const
{
    
    MPI_Status status;
    int i;
    
    // receive the point coordinates:
    int itab[2];
    MPI_Recv ( itab , 2 , MPI_INT , source ,
              Cache_Server::TAG_X1 , MPI_COMM_WORLD , &status );
    
    int n = itab[0];
    int m = itab[1];
    
    double * rtab = new double[n];
    
    MPI_Recv ( rtab , n , MPI_DOUBLE , source ,
              Cache_Server::TAG_X2 , MPI_COMM_WORLD , &status );
    
    // create the Eval_Point to search:
    Eval_Point * x = new Eval_Point ( n , m );
    
    for ( i = 0 ; i < n ; ++i )
        (*x)[i] = rtab[i];
    
    delete [] rtab;
    
    // search in cache, or stop the algorithm:
    const Eval_Point * cache_x;
    
    if ( _max_bbe > 0 && size() >= _max_bbe )
    {
        Eval_Point * stop_point = new Eval_Point ( n , m );
        for ( i = 0 ; i < n ; ++i )
            (*stop_point)[i] = (*x)[i];
        stop_point->set_eval_status ( EVAL_FAIL );
        cache_x = stop_point;
    }
    else
        cache_x = Cache::find ( *x );
    
    // cache hit signal :
    int cache_hit;
    
    // point in cache :
    if ( cache_x )
    {
        
        delete x;
        
        cache_hit = 1;
        
        ++_cache_hits;
        
        MPI_Rsend ( &cache_hit , 1 , MPI_INT , source ,
                   Cache_Server::TAG_CACHE_HIT , MPI_COMM_WORLD );
        
        // bb output values, defined values and eval_ok flag:
        rtab               = new double[m];
        char        * ctab = new char  [m+1];
        const Point & bbo  = cache_x->get_bb_outputs();
        
        for ( i = 0 ; i < m ; ++i )
        {
            if ( bbo[i].is_defined() )
            {
                rtab[i] = bbo[i].value();
                ctab[i] = '1';
            }
            else
            {
                rtab[i] = INF;
                ctab[i] = '0';
            }
        }
        
        ctab[m] = ( cache_x->is_eval_ok() ) ? '1' : '0';
        
        MPI_Send ( rtab , m , MPI_DOUBLE , source ,
                  Cache_Server::TAG_BBOR , MPI_COMM_WORLD );
        
        MPI_Send ( ctab , m+1 , MPI_CHAR , source ,
                  Cache_Server::TAG_BBOC , MPI_COMM_WORLD );
        
        delete [] rtab;
        delete [] ctab;
        
        // remove this point from _clients_ext_pts:
        {
            list<const Eval_Point *>::iterator
            it  = _clients_ext_pts[source].begin() ,
            end = _clients_ext_pts[source].end  ();
            while ( it != end )
            {
                if ( *it == cache_x )
                {
                    _clients_ext_pts[source].erase(it);
                    break;
                }
                ++it;
            }
        }
    }
    
    // point not in cache :
    else
    {
        
        cache_hit = 0;
        
        // evaluation in progress ?
        if ( _waited_pts )
        {
            
            for ( i = 0 ; i < _np ; ++i )
                if ( _waited_pts[i] && *_waited_pts[i] == *x )
                {
                    cache_hit = -1;
                    break;
                }
            
            if ( cache_hit == 0 )
                _waited_pts[source] = x;
            else
                delete x;
        }
        else
            delete x;
        
        MPI_Rsend ( &cache_hit , 1 , MPI_INT , source ,
                   Cache_Server::TAG_CACHE_HIT , MPI_COMM_WORLD );
    }
}

/*--------------------*/
/*    find a point    */
/*--------------------*/
const Eval_Point * Cache_Server::find ( const Eval_Point & x ) const
{
    
    int npm1 = _np-1;
    
    // server:
    if ( _rank == npm1 )
        return Cache::find ( x );
    
    // A. search in local cache:
    const Eval_Point * cache_x = Cache::find ( x );
    if ( cache_x )
        return cache_x;
    
    // B. ask the server.
    int i , n = x.size() , m = x.get_m();
    int itab[2];
    itab[0] = n;
    itab[1] = m;
    double * rtab = new double[n];
    for ( i = 0 ; i < n ; ++i )
        rtab[i] = x[i].value();
    
    int      cache_hit = -1;
    MPI_Request    req = MPI_REQUEST_NULL;
    MPI_Status  status;
    
    while ( cache_hit < 0 )
    {
        
        // B1. send a request for cache_hit:
        MPI_Irecv ( &cache_hit , 1 , MPI_INT , npm1 ,
                   Cache_Server::TAG_CACHE_HIT , MPI_COMM_WORLD , &req );
        
        // B2. send the find signal:
        MPI_Send ( &Cache_Server::FIND_SIGNAL , 1 , MPI_CHAR ,
                  npm1 , Cache_Server::TAG_SIGNAL , MPI_COMM_WORLD );
        
        // B3. send the point coordinates:
        MPI_Send ( itab , 2 , MPI_INT , npm1 ,
                  Cache_Server::TAG_X1 , MPI_COMM_WORLD );
        MPI_Send ( rtab , n , MPI_DOUBLE , npm1 ,
                  Cache_Server::TAG_X2 , MPI_COMM_WORLD );
        
        // B4. wait for the cache_hit request:
        MPI_Wait ( &req , &status );
        
        // cache hit possible values:
        //
        //  -1: point is being evaluated by another process (wait)
        //   0: point not in cache server
        //   1: point in cache server
    }
    
    delete [] rtab;
    
    // C. cache hit: receive the point outputs:
    if ( cache_hit == 1 )
    {
        
        // C.1. bb output values and eval status:
        rtab = new double[m];
        MPI_Recv ( rtab , m , MPI_DOUBLE , npm1 ,
                  Cache_Server::TAG_BBOR , MPI_COMM_WORLD , &status );
        
        char * ctab = new char[m+1];
        MPI_Recv ( ctab , m+1 , MPI_CHAR , npm1 ,
                  Cache_Server::TAG_BBOC , MPI_COMM_WORLD , &status );
        
        Point bbo(m);
        for ( i = 0 ; i < m ; ++i )
            if ( ctab[i]=='1' )
                bbo[i] = rtab[i];
        
        delete [] rtab;
        
        //  C.2. eval point construction:
        Eval_Point * y = new Eval_Point ( n , m );
        y->set_bb_output ( bbo );
        for ( i = 0 ; i < n ; ++i )
            (*y)[i] = x[i];
        
        y->set_eval_status ( (ctab[m]=='1') ? EVAL_OK : EVAL_FAIL );
        
        y->set_current_run ( x.get_current_run() );
        
        delete [] ctab;
        
        cache_x = y;
        
        
        // C.3. insertion in local cache:
        const_cast<Cache_Server*>(this)->Cache::insert ( *cache_x );
    }
    
    return cache_x;
}

/*------------------------------------*/
/*    process the insertion signal    */
/*------------------------------------*/
void Cache_Server::process_insert_signal ( int source )
{
    
    if ( _waited_pts && _waited_pts[source] )
    {
        delete _waited_pts[source];
        _waited_pts[source] = NULL;
    }
    
    MPI_Status status;
    
    // receive the evaluation point:
    int itab[7];
    MPI_Recv ( itab , 7 , MPI_INT , source ,
              Cache_Server::TAG_X3 , MPI_COMM_WORLD , &status );
    
    int n = itab[0];
    int m = itab[1];
    
    double * rtab = new double[n+2*m+2];
    
    MPI_Recv ( rtab , n+2*m+2 , MPI_DOUBLE , source ,
              Cache_Server::TAG_X4 , MPI_COMM_WORLD , &status );
    
    // create the Eval_Point to insert:
    Eval_Point * x = new Eval_Point ( n , m );
    
    int i;
    for ( i = 0 ; i < n ; ++i )
        (*x)[i] = rtab[i];
    
    for ( i = 0 ; i < m ; ++i )
        if ( rtab[2*i+n] > 0 )
            x->set_bb_output ( i , rtab[2*i+n+1] );
    
    if ( itab[5] == 1 )
        x->set_f ( rtab[n+2*m  ] );
    
    if ( itab[6] == 1 )
        x->set_h ( rtab[n+2*m+1] );
    
    delete [] rtab;
    
    x->set_eval_status ( ( itab[2] == 1 ) ? EVAL_OK : EVAL_FAIL );
    
    // Eval_Point insertion in cache and multiple_evals detection:
    const Eval_Point * cache_x = Cache::find ( *x );
    if ( cache_x )
    {
        ++_multiple_evals;
        delete x;
        x = const_cast<Eval_Point *>(cache_x);
    }
    else
    {
        Cache::insert ( *x );
        write_his_file ( *x );
    }
    
    // update the best points:
    update_best_points ( *x , source );
}

/*--------------------*/
/*   insert a point   */
/*--------------------*/
void Cache_Server::insert ( const NOMAD::Eval_Point & x )
{
    
    // insertion in local cache :
    Cache::insert ( x );
    
    int npm1 = _np-1;
    if ( _rank == npm1 )
        return;
    
    // insert signal :
    MPI_Send ( &Cache_Server::INSERT_SIGNAL , 1 , MPI_CHAR ,
              npm1 , Cache_Server::TAG_SIGNAL , MPI_COMM_WORLD );
    
    // send the point :
    int i , n = x.size() , m = x.get_m();
    int itab[7];
    itab[0] = n;
    itab[1] = m;
    itab[2] = ( x.is_eval_ok() ) ? 1 : 0;
    
    if ( x.get_signature() && (x.get_signature()->get_mesh()->get_mesh_indices())[0].is_defined() )
    {
        itab[3]=1;
        itab[4]=static_cast<int>((x.get_signature()->get_mesh()->get_mesh_indices())[0].value());
    }
    else
        itab[3] = itab[4] = 0;
    
    double * rtab = new double[n+2*m+2];
    for ( i = 0 ; i < n ; ++i )
        rtab[i] = x[i].value();
    
    const Point & bbo = x.get_bb_outputs();
    
    for ( i = 0 ; i < m ; ++i )
        if ( bbo[i].is_defined() )
        {
            rtab[2*i+n  ] = 1.0;
            rtab[2*i+n+1] = bbo[i].value();
        }
        else
            rtab[2*i+n] = rtab[2*i+n+1] = -1.0;
    
    // f and h values:
    if ( x.get_f().is_defined() )
    {
        itab[5    ] = 1;
        rtab[n+2*m] = x.get_f().value();
    }
    else
    {
        itab[5    ] = 0;
        rtab[n+2*m] = INF;
    }
    if ( x.get_h().is_defined() )
    {
        itab[6      ] = 1;
        rtab[n+2*m+1] = x.get_h().value();
    }
    else
    {
        itab[6      ] = 0;
        rtab[n+2*m+1] = INF;
    }
    
    MPI_Send ( itab , 7 , MPI_INT , npm1 ,
              Cache_Server::TAG_X3 , MPI_COMM_WORLD );
    
    MPI_Send ( rtab , n+2*m+2 , MPI_DOUBLE , npm1 ,
              Cache_Server::TAG_X4 , MPI_COMM_WORLD );
    
    delete [] rtab;
}

/*--------------------------------------*/
/*    get and remove an extern point    */
/*--------------------------------------*/
const Eval_Point * Cache_Server::get_and_remove_extern_point ( void ) const
{
    
    int npm1 = _np-1;
    
    if ( _rank == npm1 )
        return NULL;
    
    // extern point from the client:
    // -----------------------------
    if ( Cache::get_nb_extern_points() > 0 )
        return Cache::get_and_remove_extern_point();
    
    // extern point from the server:
    // -----------------------------
    
    int nb_pt;
    
    // send a request for an extern point:
    MPI_Request req;
    MPI_Irecv ( &nb_pt , 1 , MPI_INT , npm1 ,
               Cache_Server::TAG_EP , MPI_COMM_WORLD , &req );
    
    // extern points signal :
    MPI_Send ( &Cache_Server::EP_SIGNAL , 1 , MPI_CHAR ,
              npm1 , Cache_Server::TAG_SIGNAL , MPI_COMM_WORLD );
    
    // wait for the request:
    MPI_Status status;
    MPI_Wait ( &req , &status );
    
    if ( nb_pt == 0 )
        return NULL;
    
    // receive the extern point:
    int itab[5];
    MPI_Recv ( itab , 5 , MPI_INT , npm1 ,
              Cache_Server::TAG_X5 , MPI_COMM_WORLD , &status );
    
    int n = itab[0];
    int m = itab[1];
    
    double * rtab = new double[n+2*m];
    
    MPI_Recv ( rtab , n+2*m , MPI_DOUBLE , npm1 ,
              Cache_Server::TAG_X6 , MPI_COMM_WORLD , &status );
    
    // create the Eval_Point:
    Eval_Point * x = new Eval_Point ( n , m );
    
    int i;
    for ( i = 0 ; i < n ; ++i )
        (*x)[i] = rtab[i];
    
    for ( i = 0 ; i < m ; ++i )
        if ( rtab[2*i+n] > 0 )
            x->set_bb_output ( i , rtab[2*i+n+1] );
    
    delete [] rtab;
    
    x->set_eval_status ( ( itab[2] == 1 ) ? EVAL_OK : EVAL_FAIL );
    
    // insert the point in local cache:
    const Eval_Point * cache_x = Cache::find ( *x );
    if ( cache_x )
    {
        delete x;
        return cache_x;
    }
    
    x->set_current_run ( true );
    const_cast<Cache_Server*>(this)->Cache::insert ( *x );
    x->set_current_run ( false );
    
    return x;
}

/*---------------------------------------*/
/*    get the number of extern points    */
/*---------------------------------------*/
int Cache_Server::get_nb_extern_points ( void ) const
{
    
    int nb_client_extern_pts = Cache::get_nb_extern_points();
    int nb_server_extern_pts = 0;
    int npm1                 = _np-1;
    
    if ( _rank != npm1 )
    {
        
        // send a request for the number of extern points:
        MPI_Request req;
        MPI_Irecv ( &nb_server_extern_pts , 1 , MPI_INT , npm1 ,
                   Cache_Server::TAG_NB_EP , MPI_COMM_WORLD , &req );
        
        // number of extern points signal :
        MPI_Send ( &Cache_Server::NB_EP_SIGNAL , 1 , MPI_CHAR ,
                  npm1 , Cache_Server::TAG_SIGNAL , MPI_COMM_WORLD );
        
        // wait for the request:
        MPI_Status status;
        MPI_Wait ( &req , &status );
    }
    
    return nb_client_extern_pts + nb_server_extern_pts;
}

/*---------------------------------*/
/*    display the extern points    */
/*---------------------------------*/
void Cache_Server::display_extern_pts ( const Display & out ) const
{
    
    int npm1 = _np-1;
    
    // server:
    // -------
    if ( _rank == npm1 )
    {
        
        list<const Eval_Point*>::const_iterator it;
        out << endl << open_block ("Clients extern points");
        
        for ( int i = 0 ; i < npm1 ; ++i )
        {
            out.open_block ( "client #"+itos(i) );
            for ( it  = _clients_ext_pts[i].begin() ;
                 it != _clients_ext_pts[i].end  () ;
                 ++it )
            {
                out << "#" << (*it)->get_tag() << " ( ";
                (*it)->Point::display ( out );
                out << " ) " << " ["
                << (*it)->get_bb_outputs() << " ] f="
                << (*it)->get_f() << " h="
                << (*it)->get_h() << endl;
            }
            out.close_block();
        }
    }
    
    // clients:
    else
    {
        
        out << endl
        << open_block ( "Process #" + itos(_rank) + ": extern points" );
        
        out << "number of points = "
        << get_nb_extern_points() << endl;
        
        const Eval_Point * extern_pt = get_and_remove_extern_point();
        
        while ( extern_pt )
        {
            
            out << "#" << extern_pt->get_tag() << " ( ";
            extern_pt->Point::display ( out );
            out << " ) " << " ["
            << extern_pt->get_bb_outputs() << " ] f="
            << extern_pt->get_f() << " h="
            << extern_pt->get_h() << endl;
            
            extern_pt = get_and_remove_extern_point();
        }
    }
    out << close_block() << endl;
}

/*--------------------------------------*/
/*  update and display the best points  */
/*--------------------------------------*/
void Cache_Server::update_best_points ( const Eval_Point & x      ,
                                       int                 source   )
{
    const Double & f = x.get_f();
    const Double & h = x.get_h();
    
    if ( !f.is_defined() || !h.is_defined() )
        return;
    
    int  i;
    bool add_x = false;
    
    // feasible:
    if ( h <= _h_min )
    {
        
        // new best feasible point:
        if ( !_bf || f < _bf->get_f() )
        {
            _bf   = &x;
            add_x = true;
            display_current_solution();
        }
    }
    
    // infeasible:
    else
    {
        if ( !_bi1 || h < _bi1->get_h() )
        {
            _bi1  = &x;
            add_x = true;
        }
        if ( !_bi2 || f < _bi2->get_f() )
        {
            _bi2  = &x;
            add_x = true;
        }
    }
    
    if ( add_x )
        for ( i = 0 ; i < _np-1 ; ++i )
            if ( i != source )
                _clients_ext_pts[i].push_front ( &x );
}

/*-----------------------------------------*/
/*    display the current best solution    */
/*-----------------------------------------*/
void Cache_Server::display_current_solution ( void ) const
{
    if ( _rank != _np-1 || !_bf )
        return;
    if ( _debug )
        _out << "CACHE SERVER: CURRENT SOLUTION: \t";
    _out << _clock.get_real_time() << "\t"
    << size() << "\t" << _bf->get_f() << endl;
}

/*-------------------------------*/
/*    display the best points    */
/*-------------------------------*/
void Cache_Server::display_best_points ( const Display & out ) const
{
    if ( _rank != _np-1 )
        return;
    
    // display the last solution:
    display_current_solution();
    
    // stats:
    out << "evaluations: " << size()      << endl
    << "cache hits : " << _cache_hits << endl;
    
    // best feasible solution:
    out << "best feasible solution: ";
    if ( _bf )
    {
        out << "x=( ";
        _bf->Point::display(out);
        out << " )"
        << " F(x)=[ " << _bf->get_bb_outputs() << " ] h="
        << _bf->get_h() << " f=" << _bf->get_f() << endl;
    }
    else
    {
        
        out << "NULL" << endl;
        
        // best infeasible solutions:
        if ( _bi1 )
        {
            out << "best infeas. sol. #1  : x=( ";
            _bi1->Point::display(out);
            out << " )"
            << " F(x)=[ " << _bi1->get_bb_outputs() << " ] h="
            << _bi1->get_h() << " f=" << _bi1->get_f() << endl;
        }
        
        if ( _bi2 && _bi2 != _bi1 )
        {
            out << "best infeas. sol. #2  : x=( ";
            _bi2->Point::display(out);
            out << " )"
            << " F(x)=[ " << _bi2->get_bb_outputs() << " ] h="
            << _bi2->get_h() << " f=" << _bi2->get_f() << endl;
        }
    }
}


void Cache_Server::write_his_file ( const NOMAD::Eval_Point & x ) const
{
    
    if ( ! _history_file.empty()  )
    {
        std::ofstream fout;
        
        fout.open ( _history_file.c_str() , std::ios::app );
        
        cout << "about to write in history file ..." ;
        
        if ( !fout.fail() )
        {
            fout.setf      ( std::ios::fixed             );
            fout.precision ( NOMAD::DISPLAY_PRECISION_BB );
            
            x.Point::display ( fout , " " , -1 , -1 );
            fout << " ";
            x.get_bb_outputs().Point::display ( fout , " " , -1 , -1 );
            fout << std::endl;
            
            cout << " done" <<endl;
            fout.close();
        }
        else
        {
            cout << "Warning ( Cache_Server.cpp, line " << __LINE__
            << "): could not update the history"
            << " in \'"
            << _history_file << "\'" << std::endl << std::endl;
        }
        
    }
    
}
