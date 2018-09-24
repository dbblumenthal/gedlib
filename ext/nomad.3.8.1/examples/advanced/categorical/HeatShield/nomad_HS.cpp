/*------------------------------------------*/
/*  The Heatshield problem in library mode  */
/*------------------------------------------*/
#include "nomad.hpp"
#include "HS_Param.hpp"
using namespace std;
using namespace NOMAD;

const int LOAD_FLAG = 2;

const double L      = 100;
const double T_COLD = 4.2;
const double T_HOT  = 300;

const double EPS = 1e-13;

const int N_MAX  = 10;
const int NB_MAT = 7;

// BB_OUTPUT_TYPE depends on Load_Flag values:
// -------------------------------------------
//   0 : OBJ EB EB
//   1 : OBJ EB EB
//   2 : OBJ EB EB PB PB
//   3 : OBJ EB EB PB PB PB

// Variables indexes:
// ------------------
// 0      n
// ------------
// 1      M[0]
// ...
// n+1    M[n]     (n+1 materials)
// ------------
// n+2    dx[0]
// ...
// 2n+1   dx[n-1]  (dx's)
// ------------
// 2n+2   T[1]     (T's)
// ...
// 3n+1   T[n]
// ------------

// total: 3n+2 variables

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
class My_Evaluator : public Evaluator {
    
private:
    
    HS_Param _param;
    
public:
    
    My_Evaluator  ( const Parameters & p ) :
    Evaluator ( p ) ,
    _param    ( LOAD_FLAG ,
               250000    ,       // Load
               10        ,       // Max_Weight
               0.05      ,       // TE_Limit
               T_COLD    ,       // T_cold
               T_HOT     ,       // T_hot
               L           )
    {
    }  // L
    
    ~My_Evaluator ( void )
    {
    }
    
    bool eval_x ( Eval_Point          & x          ,
                 const NOMAD::Double & h_max      ,
                 bool                & count_eval   ) const;
};

/*--------------------------------------------------*/
/*  user class to define categorical neighborhoods  */
/*--------------------------------------------------*/
class My_Extended_Poll : public Extended_Poll {
    
private:
    
    Signature ** _signatures;
    int          _km;
    int          _km_tag;
    
    // create 3 signatures from the current one (n , n-1, n+1):
    void update_signatures ( const Signature & cur_signature , int n );
    
    // create _signatures[n] :
    void create_signature ( int n );
    
public:
    
    // constructor:
    My_Extended_Poll ( Parameters & p );
    
    // destructor:
    virtual ~My_Extended_Poll ( void );
    
    // construct the extended poll points:
    virtual void construct_extended_points ( const Eval_Point & xk );
    
};

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv )
{
    
    
    // NOMAD initializations:
    NOMAD::begin ( argc , argv );
    
    // display:
    NOMAD::Display out ( cout );
    out.precision ( DISPLAY_PRECISION_STD );
    
    try {
        
        // parameters creation:
        Parameters p ( out );
        
        // (initial) number of variables:
        int i , n = 1 , dim = 3*n+2;
        p.set_DIMENSION ( dim );
        
        // definition of output types:
        {
            vector<bb_output_type> bbot;
            bbot.push_back ( OBJ );
            bbot.push_back ( EB  );
            bbot.push_back ( EB  );
            if ( LOAD_FLAG > 1 )
            {
                
                bbot.push_back ( PB );
                bbot.push_back ( PB );
                if ( LOAD_FLAG == 3 )
                    bbot.push_back ( PB );
            }
            p.set_BB_OUTPUT_TYPE ( bbot );
        }
        
        // starting point:
        Point x0 ( dim );
        x0[0] = n;
        x0[1] = 0;
        x0[2] = 1;
        x0[3] = 50;
        x0[4] = 150;
        
        p.set_X0 ( x0 );
        
        // categorical variables (n and materials):
        for ( i = 0 ; i <= n+1 ; ++i )
            p.set_BB_INPUT_TYPE ( i , CATEGORICAL );
        
        // display degree:
        p.set_DISPLAY_DEGREE ( NORMAL_DISPLAY );
        
        // display stats:
        p.set_DISPLAY_STATS ( "bbe ( sol ) obj" );
        
        // bounds:
        Point lb ( dim , 0 );
        Point ub ( dim );
        for ( i = 2*n+2 ; i <= 3*n+1 ; ++i )
        {
            
            lb[i] = T_COLD;
            ub[i] = T_HOT;
        }
        ub[0] = N_MAX;
        for ( i = 1 ; i <= n+1 ; ++i )
            ub[i] = NB_MAT-1;
        for ( i = n+2 ; i <= 2*n+1 ; ++i )
            ub[i] = L;
        
        p.set_LOWER_BOUND ( lb );
        p.set_UPPER_BOUND ( ub );
        
       
        // max number of evaluations:
        p.set_MAX_BB_EVAL ( 30000 );
        
        // extended poll trigger:
        p.set_EXTENDED_POLL_TRIGGER ( 1.1 , true );
        
        p.set_OPPORTUNISTIC_EVAL(false);
       p.set_DIRECTION_TYPE(ORTHO_2N);
	p.set_DISABLE_EVAL_SORT();
//	p.set_DIRECTION_TYPE(LT_2N);
//	p.set_DIRECTION_TYPE(ORTHO_NP1_QUAD);        
        // parameters validation:
        p.check();
        
        // custom evaluator creation:
        My_Evaluator ev   ( p );
        
        // extended poll:
        My_Extended_Poll ep ( p );
        
        // algorithm creation and execution:
        Mads mads ( p , &ev , &ep , NULL , NULL );
        mads.run();
        
        // display the solution and the matlab format:
        const Eval_Point * xe = mads.get_best_feasible();
        if ( xe )
        {
            
            
            // the solution:
            out << endl << "solution:" << endl
            << "\tf=" << xe->get_f() << endl << endl;
            n = static_cast<int> ( (*xe)[0].value() );
            out << "\tn=" << n << endl << endl;
            for ( i = 1 ; i <= n+1 ; ++i )
                out << "\tmat[" << setw(2) << i-1 << "]=" << (*xe)[i] << " ("
                << static_cast<material_type> ( (*xe)[i].round() )
                << ")" << endl;
            out << endl;
            for ( i = n+2 ; i <= 2*n+1 ; ++i )
                out << "\tdx [" << setw(2) << i-n-2 << "]=" << (*xe)[i] << endl;
            out << endl;
            for ( i = 2*n+2 ; i <= 3*n+1 ; ++i )
                out << "\tT  [" << setw(2) << i-2*n-1 << "]=" << (*xe)[i] << endl;
            
            // the solution in matlab format:
            cout << endl << "Matlab:" << endl << "x = [ ";
            for ( i = n+2 ; i <= 2*n+1 ; ++i )
                cout << (*xe)[i] << " ";
            for ( i = 2*n+2 ; i <= 3*n+1 ; ++i )
                cout << (*xe)[i] << " ";
            cout << "]'" << endl;
            
            cout << "p = { [" << n << "] ";
            for ( i = 1 ; i <= n+1 ; ++i )
                cout << ", \'" << static_cast<material_type> ( (*xe)[i].round() )
                << "\' ";
            cout << "}" << endl <<endl;
        }
    }
    catch ( exception & e )
    {
        
        cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
        return EXIT_FAILURE;
    }
    
    NOMAD::end();
    
    return EXIT_SUCCESS;
}

/*------------------*/
/*  the evaluation  */
/*------------------*/
bool My_Evaluator::eval_x ( Eval_Point          & x          ,
                           const NOMAD::Double & h_max      ,
                           bool                & count_eval   ) const
{
    
    // 1. number of heat intercepts:
    int i , cur = 1 , n = int(x[0].value());
    
    // 2. materials:
    int * material = new int[n+1];
    for ( i = 0 ; i <= n ; ++i )
        material[i] = int(x[cur++].value());;
    
    // 3. dx:
    double * dx     = new double[n+1];
    double   sumx   = 0.0;
    double   min_dx = 999;
    for ( i = 0 ; i < n ; ++i )
    {
        
        dx[i] = x[cur++].value();
        if ( dx[i] < min_dx )
            min_dx = dx[i];
        sumx += dx[i];
    }
    dx[n] = _param.get_L() - sumx;
    if ( dx[n] < min_dx )
        min_dx = dx[n];
    
    // 4. T:
    double * T       = new double[n+2];
    double * objterm = new double[n+2];
    double   cycle , v , max_T_viol = 0.0;
    
    for ( i = 0 ; i <= n+1 ; ++i )
    {
        
        
        if ( i==0 )
            T[0  ] = _param.get_T_cold();
        else
        {
            
            if ( i==n+1 )
                T[n+1] = _param.get_T_hot();
            else
                T[i] = x[cur++].value();
            
            v = T[i-1] - T[i];
            if ( v > max_T_viol )
                max_T_viol = v;
        }
        
        if ( T[i] <= 4.2 )
            cycle = 5;
        else cycle = ( T[i] >= 71 ) ? 2.5 : 4;
        
        objterm[i] = cycle * ( _param.get_T_hot() / T[i] - 1 );
    }
    
    /*--------------*/
    /*  evaluation  */
    /*--------------*/
    
    double        Ti , Tip1;
    int           indfirst , indlast;
    HS_Material * mat;
    
    double        TC_Integral , TE_Integral , q , a;
    double        power       = 0;
    double        cx1         = 0;
    double        cx2         = 0;
    double        cx3         = 0;
    
    if ( dx[n] > -EPS && max_T_viol < EPS && min_dx > EPS )
    {
        
        
        for ( i = 0 ; i <= n ; ++i )
        {
            
            
            Ti   = T[i  ];
            Tip1 = T[i+1];
            
            double  * x;
            int      nx;
            
            mat = _param.get_material ( material[i] );
            
            if ( mat->get_T_indexes ( Ti          ,
                                     Tip1        ,
                                     indfirst    ,
                                     indlast       ) )
            {
                
                nx   = 4;
                x    = new double [4];
                x[0] = Ti;
                x[1] = ( Ti   + mat->get_T ( indfirst ) ) / 2.0;
                x[2] = ( Tip1 + mat->get_T ( indlast  ) ) / 2.0;
                x[3] = Tip1;
            }
            else
            {
                
                nx   = 3;
                x    = new double[3];
                x[0] = Ti;
                x[1] = (Ti + Tip1) / 2.0;
                x[2] = Tip1;
            }
            
            double * y = new double [nx];
            
            mat->thermal_conductivity ( nx , x , y );
            
            double sum_TCI_Ends = ( nx==3 ) ?
            ( y[0] + 4*y[1] + y[2] ) * (Tip1-Ti) / 6 :
            ( y[0] + 4*y[1] + mat->get_k(indfirst) ) *
            ( mat->get_T(indfirst) - Ti  ) / 6
            + ( mat->get_k(indlast) + 4*y[2] + y[3]  ) *
            ( Tip1 - mat->get_T(indlast) ) / 6;
            
            double sum_TEI_Ends = 0;
            
            if ( _param.get_load_flag() >= 2 )
            {
                
                
                mat->thermal_expansion ( nx , x , y );
                
                sum_TEI_Ends = ( nx==3 ) ?
                ( y[0] + 4*y[1] + y[2] ) * (Tip1-Ti) / 6 :
                ( y[0] + 4*y[1] + mat->get_k(indfirst)*mat->get_e(indfirst) ) *
                ( mat->get_T(indfirst) - Ti  ) / 6
                + ( mat->get_k(indlast)*mat->get_e(indlast) + 4*y[2] + y[3]  ) *
                ( Tip1 - mat->get_T(indlast) ) / 6;
            }
            
            delete [] x;
            delete [] y;
            
            TC_Integral = sum_TCI_Ends;
            TE_Integral = sum_TEI_Ends;
            for ( int k = indfirst ; k < indlast ; ++k )
            {
                
                TC_Integral += mat->get_TCI(k);
                TE_Integral += mat->get_TEI(k);
            }
            
            double stress;
            
            if ( TC_Integral < EPS )
            {
                
                stress = -1e20;
                cx2 = cx3 = 1e20;
                a = q = 0;
            }
            
            else
            {
                
                
                stress = mat->thermal_stress ( Tip1 );
                
                if ( _param.get_load_flag() == 3 )
                {
                    
                    
                    double xtmp[2];
                    double ytmp[2];
                    xtmp[0] = Ti;
                    xtmp[1] = Tip1;
                    mat->thermal_YM ( 2 , xtmp , ytmp );
                    
                    double young_modulus = ( ytmp[0] > ytmp[1] ) ? ytmp[0] : ytmp[1];
                    
                    stress -= young_modulus * TE_Integral / TC_Integral;
                    
                    if ( stress < 0 )
                        cx3 += stress*stress;
                }
                
                if ( _param.get_load_flag() == 0 )
                    q = TC_Integral / dx[i];
                else
                {
                    
                    a = _param.get_load() / stress;
                    q = a * TC_Integral / dx[i];
                    
                    if ( _param.get_load_flag() > 1 )
                    {
                        
                        cx1 += mat->get_density() * dx[i] * a;
                        cx2 += (dx[i]*TE_Integral/TC_Integral);
                    }
                }
            }
            power += q * ( objterm[i] - objterm[i+1] );
        }
        
        count_eval = true;
        
        dx[n] = -0.0;
        
    }
    else
    {
        
        power      = 1e20 / L;
        count_eval = false;
    }
    
    // constraints:
    x.set_bb_output  ( 1 , -dx[n] );
    x.set_bb_output  ( 2 , max_T_viol*max_T_viol );
    
    if ( _param.get_load_flag() > 1 )
    {
        
        if ( _param.get_load_flag() == 3 )
            x.set_bb_output  ( 5 , cx3 );
        x.set_bb_output  ( 3 , cx1 / _param.get_max_weight() - 1                );
        x.set_bb_output  ( 4 , cx2 / (_param.get_TE_limit()*_param.get_L()) - 1 );
    }
    
    // f(x):
    x.set_bb_output  ( 0 , power * L );
    
    // clean memory:
    delete [] material;
    delete [] dx;
    delete [] T;
    delete [] objterm;
    
    return true; // the evaluation succeeded
}

/*--------------------------------*/
/*  My_Extended_Poll constructor  */
/*--------------------------------*/
My_Extended_Poll::My_Extended_Poll ( Parameters & p ) :
Extended_Poll ( p                         ) ,
_signatures   ( new Signature * [N_MAX+1] ) ,
_km           (  1                        ) ,
_km_tag       ( -1                        )
{
    for ( int i = 0 ; i <= N_MAX ; ++i )
        _signatures[i] = NULL;
}

/*-------------------------------*/
/*  My_Extended_Poll destructor  */
/*-------------------------------*/
My_Extended_Poll::~My_Extended_Poll ( void )
{
    
    for ( int i = 0 ; i <= N_MAX ; ++i )
        delete _signatures[i];
    delete [] _signatures;
}

/*-------------------------*/
/*  create _signatures[n]  */
/*-------------------------*/
void My_Extended_Poll::create_signature ( int n )
{
    
    // C.Tribes nov 25, 2015 --- Several modifs for the new versions
    
    
    NOMAD::Double dm0_dx , dm0_T , dmmin_dx , dmmin_T , dpmin_dx , dpmin_T;
    
    {
        int nref = (_p.get_dimension() - 2) / 3;
        
        const Point & delta_0   = _p.get_initial_mesh_size();
        const Point & delta_min = _p.get_min_mesh_size();
        const Point & Delta_min = _p.get_min_poll_size();
        
        bool delta_min_def = delta_min.is_defined();
        bool Delta_min_def = Delta_min.is_defined();
        
        dm0_dx   = delta_0 [   nref+2 ];
        dm0_T    = delta_0 [ 2*nref+2 ];
        dmmin_dx = (delta_min_def) ? delta_min [   nref+2 ] : NOMAD::Double();
        dmmin_T  = (delta_min_def) ? delta_min [ 2*nref+2 ] : NOMAD::Double();
        dpmin_dx = (Delta_min_def) ? Delta_min [   nref+2 ] : NOMAD::Double();
        dpmin_T  = (Delta_min_def) ? Delta_min [ 2*nref+2 ] : NOMAD::Double();
    }
    
    int                   i , dim = 3*n+2;
    vector<bb_input_type> bbit               ( dim );
    Point                 initial_poll_size  ( dim , 1.0 );
    Point                 min_mesh_size      ( dim );
    Point                 min_poll_size      ( dim );
    Point                 lb                 ( dim , 0.0   );
    Point                 ub                 ( dim , N_MAX );
    Point                 fixed_variables    ( dim );
    vector<bool>          periodic_variables ( dim );
    
//    set<int>              var_ind_1;
//    set<int>              var_ind_2;
    
    // n and materials:
    for ( i = 0 ; i <= n+1 ; ++i )
    {
        
        bbit[i] = CATEGORICAL;
//       var_ind_1.insert ( i );
    }
    
    // dx's:
    for ( i = n+2 ; i <= 2*n+1 ; ++i )
    {
        
        bbit             [i] = CONTINUOUS;
        initial_poll_size[i] = dm0_dx;
        min_mesh_size    [i] = dmmin_dx;
        min_poll_size    [i] = dpmin_dx;
        ub               [i] = L;
//        var_ind_2.insert ( i );
    }
    
    // T's:
    for ( i = 2*n+2 ; i < dim ; ++i )
    {
        
        bbit             [i] = CONTINUOUS;
        initial_poll_size[i] = dm0_T;
        min_mesh_size    [i] = dmmin_T;
        min_poll_size    [i] = dpmin_T;
        lb               [i] = T_COLD;
        ub               [i] = T_HOT;
//        var_ind_2.insert ( i );
    }

// C.Tribes nov 25, 2015 --- no halton seed in new version
//    // construct a new Halton seed:
//    int * primes  = new int [dim];
//    construct_primes ( dim , primes );
//    int halton_seed_ref = primes[dim-1];
//    delete [] primes;
    
//    Variable_Group * vg1 = new Variable_Group ( var_ind_1                   ,
//                                               _p.get_direction_types   () ,
//                                               _p.get_sec_poll_dir_types() ,
//                                               _p.out()                        );
//    
//    Variable_Group * vg2 = new Variable_Group ( var_ind_2                            ,
//                                               _p.get_direction_types()             ,
//                                               _p.get_sec_poll_dir_types()          ,
//                                               _p.out()                            );
//    
//    set<Variable_Group*,VG_Comp> var_groups;
//    var_groups.insert ( vg1 );
//    var_groups.insert ( vg2 );
    
   
    NOMAD::Double mesh_update_basis = 1.0;
    NOMAD::Double poll_update_basis = 1.0;
    int mesh_coarsening_exponent = 1;
    int mesh_refining_exponent = 1;

    
//// create the signature:
//    _signatures[n] = new Signature ( dim                ,
//                                    bbit               ,
//                                    lb                 ,
//                                    ub                 ,
//                                    NOMAD::XMESH       ,
//                                    false              ,
//                                    NOMAD::Point()     ,
//                                    initial_poll_size  ,
//                                    min_poll_size      ,
//                                    min_mesh_size      ,
//                                    mesh_update_basis  ,
//                                    poll_update_basis  ,
//                                    mesh_coarsening_exponent,
//                                    mesh_refining_exponent,
//                                    0                  ,
//                                    NOMAD::Point()     ,
//                                    fixed_variables    ,
//                                    periodic_variables ,
//                                    var_groups            );
//    
//    
//    delete vg1;
//    delete vg2;
    
    
    // create the signature:
    _signatures[n] = new Signature ( dim               ,
                                    bbit               ,
                                    initial_poll_size  ,
                                    lb                 ,
                                    ub                 ,
                                    _p.get_direction_types   () ,
                                    _p.get_sec_poll_dir_types() ,
                                    _p.out() );

    
}

/*-----------------------------------------------------------*/
/*  create 3 signatures from the current one (n , n-1, n+1)  */
/*-----------------------------------------------------------*/
void My_Extended_Poll::update_signatures ( const Signature & cur_signature , int n )
{
    
    
    // n :
    // ---
    if ( !_signatures[n] )
        _signatures[n] = new Signature ( cur_signature );
    
    // n-1 :
    // -----
    if ( n > 1 && !_signatures[n-1] )
        create_signature ( n-1 );    
    
    // n+1 :
    // -----
    if ( n < N_MAX && !_signatures[n+1] )
        create_signature ( n+1 );
}

/*--------------------------------------*/
/*  construct the extended poll points  */
/*--------------------------------------*/
void My_Extended_Poll::construct_extended_points ( const Eval_Point & xk ) 
{
    
    
    int tag = xk.get_tag();
    int i , k , mat , n = xk[0].round();
    
    // create/update the signatures:
    update_signatures ( *xk.get_signature() , n );
    
    
    // cout << endl << "NEIGHBORS OF " << xk << " :\n";
    
    // 1. include neighbors in which 1 insulation material is changed:
    {
        int incr_m;
        if ( tag == _km_tag ) 
        {
            
            incr_m = ++_km;
        }
        else 
        {
            
            incr_m = _km = 1;
        }
        if ( incr_m > NB_MAT ) 
        {
            
            incr_m = _km = 1;
        }
        for ( i = 1 ; i <= n+1 ; ++i ) 
        {
            
            mat = static_cast<int>(xk[i].value()+incr_m)%NB_MAT;
            Point ep1 ( xk );
            ep1[i] = mat;
            // cout << "\tN1(" << i << "/" << n+1 << ") = " << ep1 << endl;
            add_extended_poll_point ( ep1 , *_signatures[n] );
        }
    }
    
    // 2. add an intercept:
    if ( n < N_MAX ) 
    {
        
        
        NOMAD::Double xn = L;
        for ( k = n+2 ; k <= 2*n+1 ; ++k )
            xn -= xk[k];
        
        for ( k = 0 ; k <= n ; ++k ) 
        {
            
            
            Point ep2 ( 3*n + 5 );
            ep2[0] = n+1;
            
            // materials:
            for ( i = 1 ; i <= k+1 ; ++i )
                ep2[i] = xk[i];
            for ( i = k+2 ; i <= n+2 ; ++i )
                ep2[i] = xk[i-1];
            
            // dx's:
            for ( i = n+3 ; i < k+n+3 ; ++i )
                ep2[i] = xk[i-1];
            ep2[k+n+3] = ( k+n+2 <= 2*n+1 ) ? xk[k+n+2]/2 : xn/2;
            if ( k+n+4 <= 2*n+3 )
                ep2[k+n+4] = xk[k+n+2]/2;
            for ( i = k+n+5 ; i <= 2*n+3 ; ++i )
                ep2[i] = xk[i-2];
            
            // T's:
            if ( k == 0 )
                ep2[2*n+4] = 0.5 * ( T_COLD + xk[2*n+2] );
            else if ( k == n )
                ep2[3*n+4] = 0.5 * ( T_HOT  + xk[3*n+1] );
            else 
            {
                
                ep2[2*n+4+k] = 0.5 * ( xk[2*n+1+k] + xk[2*n+2+k] );
            }
            for ( i = 2*n+4 ; i < 2*n+4+k ; ++i )
                ep2[i] = xk[i-2];
            for ( i = 2*n+5+k ; i <= 3*n+4 ; ++i )
                ep2[i] = xk[i-3];
            
            // cout << "\tN2(" << k+1 << "/" << n+1 << ") = " << ep2 << endl;
            add_extended_poll_point ( ep2 , *_signatures[n+1] );
            
        }
    }
    
    // 3. remove an intercept:
    if ( n > 1 ) 
    {
        
        
        for ( k = 0 ; k < n ; ++k ) 
        {
            
            
            Point ep3 ( 3*n - 1 );
            ep3[0] = n-1;
            
            // materials:
            for ( i = 1 ; i <= k ; ++i )
                ep3[i] = xk[i];
            for ( i = k+1 ; i <= n ; ++i )
                ep3[i] = xk[i+1];
            
            // dx's:
            NOMAD::Double d = xk[n+2+k] / (n+1);
            for ( i = n+1+k ; i <= 2*n-1 ; ++i )
                ep3[i] = xk[i+2] + d;
            for ( i = n+1 ; i <= n+k ; ++i )
                ep3[i] = xk[i+1] + d;
            
            // T's:
            for ( i = 2*n+k ; i <= 3*n-2 ; ++i )
                ep3[i] = xk[i+3];
            for ( i = 2*n ; i < 2*n+k ; ++i )
                ep3[i] = xk[i+2];
            
            // cout << "\tN3(" << k+1 << "/" << n << ") = " << ep3 << endl;
            add_extended_poll_point ( ep3 , *_signatures[n-1] );
        }
    }
    
    _km_tag = tag;
    
    // cout << endl;
}
