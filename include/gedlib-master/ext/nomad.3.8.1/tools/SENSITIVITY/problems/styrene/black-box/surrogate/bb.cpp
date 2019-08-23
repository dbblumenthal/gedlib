#include "servor.hpp"
#include "profitability.hpp"
using namespace std;

/*------------------------------------------------------------------------*/
/*                              fonction principale                       */
/*------------------------------------------------------------------------*/
int main ( int argc , char ** argv ) {

  double g0  = 1e+20;
  double g1  = 1e+20;
  double g2  = 1e+20;
  double g3  = 1e+20;
  double g4  = 1e+20;
  double g5  = 1e+20;
  double g6  = 1e+20;
  double g7  = 1e+20;
  double g8  = 1e+20;
  double g9  = 1e+20;
  double g10 = 1e+20;
  double f   = 1e+20;

  bool        OK;
  ifstream    in;
  double      d;
  int         i;
  int         i_stream , i_chem;
  double      x  [8];
  long double tmp[8];
  double      raw_cost;
  double      util_cost;
  double      Itot;
  double      Coper;
  double      Rtot;
  double      max;
  double      mtot;
  double      m;
  double      purity;
  int         nb_chem  = 7;
  int         nb_s     = 15;
  int         nb_u     = 11;
  int         n        = 8;
  double      l    [8] = { 600  , 2  , 0.0001 , 0.0001 , 2  , 0.01 , 0.1 , 300 };
  double      u    [8] = { 1100 , 20 , 0.1    , 0.1    , 20 , 0.5  , 5   , 500 };
  double      list [7] = { 0.5 , 0 , 0 , 0 , 0 , 0 , 0 };
  double      price[7] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 };
  profitability * P        = NULL;
  cashflow      * F        = NULL;
  chemical     ** chem     = NULL;
  stream       ** s        = NULL;
  servor        * units    = NULL;
  string        * safe     = NULL;

  // verif. du nombre d'arguments :
  if (argc!=2) {
    // cout << "\nargc != 2\n\n";
    goto TERMINATE;
  }

  // lecture et scaling de x :
  // -------------------------
  in.open (argv[1]);

  for ( i = 0 ; i < n ; i++ )
    in >> tmp[i];

  in.close();

  if (in.fail()) {
    // cout << "\nin.fail (1)\n\n";
    goto TERMINATE;
  }

  //for ( i = 0 ; i < n ; i++ )
  //  cout << "tmp[" << i << "] = " << tmp[i] << endl;
  //cout << endl;

  for ( i = 0 ; i < n ; i++ )
    x[i] = (u[i]-l[i])*((double)tmp[i]/100.0) + l[i];
  
  //for ( i = 0 ; i < n ; i++ )
  //  cout << "x[" << i << "] = " << x[i] << endl;

  // verifs de x (j'ai pris ca dans checkup et c'est tout ce dont
  // ca a besoin de checker ici) :
  if ( x[6] < EPS || x[2] < 0 || x[2] > 1 || x[3] < 0 || x[3] > 1 )
    goto TERMINATE;

  // chemicals :
  // -----------
  chem    = new chemical * [nb_chem];
  chem[0] = new chemical ("100-41-4" );
  chem[1] = new chemical ("1333-74-0");
  chem[2] = new chemical ("108-88-3" );
  chem[3] = new chemical ("74-82-8"  );
  chem[4] = new chemical ("71-43-2"  );
  chem[5] = new chemical ("74-85-1"  );
  chem[6] = new chemical ("100-42-5" );

  price[6] = 1.39;
  price[2] = 0.64;
  price[0] = 0.11;
  price[4] = 1.19;

  // streams :
  // ---------
  s     = new stream * [nb_s];
  s[ 0] = new stream ( "feed"    , nb_chem , chem );
  s[ 1] = new stream ( "2"       , nb_chem , chem );
  s[ 2] = new stream ( "3"       , nb_chem , chem );
  s[ 3] = new stream ( "4"       , nb_chem , chem );
  s[ 4] = new stream ( "5"       , nb_chem , chem );
  s[ 5] = new stream ( "6"       , nb_chem , chem );
  s[ 6] = new stream ( "7"       , nb_chem , chem );
  s[ 7] = new stream ( "8"       , nb_chem , chem );
  s[ 8] = new stream ( "9"       , nb_chem , chem );
  s[ 9] = new stream ( "10"      , nb_chem , chem );
  s[10] = new stream ( "back"    , nb_chem , chem );
  s[11] = new stream ( "12"      , nb_chem , chem );
  s[12] = new stream ( "stack"   , nb_chem , chem );
  s[13] = new stream ( "out-bz"  , nb_chem , chem );
  s[14] = new stream ( "out-sty" , nb_chem , chem );

  // initial conditions on streams :
  // -------------------------------
  s[0]->P = 1.0;
  s[0]->T = 298;
  s[0]->set(list);

  // units settings and calculation sequence :
  // -----------------------------------------
  units = new servor ( nb_u , nb_s , s );
  safe  = new string[nb_u];

  units->type[ 0] = "mix";
  units->name[ 0] = safe[ 0] = "mixfeed";
  units->type[ 1] = "pump";
  units->name[ 1] = safe[ 1] = "pump";
  units->type[ 2] = "heatx";
  units->name[ 2] = safe[ 2] = "heater";
  units->type[ 3] = "reactor";
  units->name[ 3] = safe[ 3] = "pfr";
  units->type[ 4] = "heatx";
  units->name[ 4] = safe[ 4] = "cooler";
  units->type[ 5] = "flash";
  units->name[ 5] = safe[ 5] = "degasor";
  units->type[ 6] = "column";
  units->name[ 6] = safe[ 6] = "sep-sty";
  units->type[ 7] = "split";
  units->name[ 7] = safe[ 7] = "spliter";
  units->type[ 8] = "column";
  units->name[ 8] = safe[ 8] = "sep-bz";
  units->type[ 9] = "loop";
  units->name[ 9] = safe[ 9] = "looping";
  units->type[10] = "burner";
  units->name[10] = safe[10] = "fire";

  // executing the calculation sequence :
  //-------------------------------------
  double y[14];
  y[7] = y[8] = 1000.0;

  if (!units->solve_process(x,y))
    goto TERMINATE;

  // on recupere les resultats :
  // ---------------------------
  y[6] = s[0]->m;
  y[3] = s[13]->m;
  y[2] = s[13]->chem[4]->m;
  y[1] = s[14]->m;
  y[0] = s[14]->chem[6]->m;

  y[9 ] = 1e20;
  y[10] = 1e20;
  y[11] = units->get_costs_sum() * 6.192;
  y[12] = 1e20;
  y[13] = 1e20;


  // arrondis :
  if (ARRONDI) {
    y[ 6] = arrondi ( y[ 6] , 4 );
    y[ 3] = arrondi ( y[ 3] , 4 );
    y[ 2] = arrondi ( y[ 2] , 4 );
    y[ 1] = arrondi ( y[ 1] , 4 );
    y[ 0] = arrondi ( y[ 0] , 4 );
    y[11] = arrondi ( y[11] , 6 );
  }

  g0 = ( y[0] > 0 && y[0] < 1e+20 ) ? 0.0 : 1e20;
  g1 = ( y[4] <= 80 ) ? 0.0 : 1.0;
  g2 = ( y[5] <= 80 ) ? 0.0 : 1.0;
  g3 = ( y[7] <= 200 && y[8] <= 8 ) ? 0.0 : 1.0;
  g4 = ( y[1] > 0 && y[1] < 1e20 ) ? (0.99-y[0]/y[1])/0.99 : 1e20;
  g5 = ( y[3] > 0 && y[3] < 1e20 ) ? (0.99-y[2]/y[3])/0.99 : 1e20;
  g6 = ( y[6] > 0 && y[6] < 1e20 ) ? (0.6-y[0]/y[6])/0.6 : 1e20;

  // bloc econo :
  // ------------

  // raw_cost :
  raw_cost = 0.0;
  if ( s[0]->m > EPS ) {
    if ( s[0]->chem[6]->m > EPS ) {
      d = (ARRONDI) ? arrondi ( s[0]->chem[6]->m , 4 ) : s[0]->chem[6]->m;
      raw_cost += d * 1.39;
    }
    if ( s[0]->chem[2]->m > EPS ) {
      d = (ARRONDI) ? arrondi ( s[0]->chem[2]->m , 4 ) : s[0]->chem[2]->m;
      raw_cost += d * 0.64;
    }
    if ( s[0]->chem[0]->m > EPS ) {
      d = (ARRONDI) ? arrondi ( s[0]->chem[0]->m , 4 ) : s[0]->chem[0]->m;
      raw_cost += d * 0.11;
    }
    if ( s[0]->chem[4]->m > EPS ) {
      d = (ARRONDI) ? arrondi ( s[0]->chem[4]->m , 4 ) : s[0]->chem[4]->m;
      raw_cost += d * 1.19;
    }
  }
  // raw_cost = raw_cost*3600.0*nb_h*nb_d :
  raw_cost *= 25920000;

  // util_cost :
  util_cost = (units->get_power_sum() * 0.000125 + units->get_water_sum() * 0.00008) * 25920000;

  // Coper :
  Itot  = y[11];
  Coper = 0.16 * Itot + 2.26 * raw_cost + util_cost;


  // Rtot :
  Rtot = 0.0;
  for ( i_stream = 13 ; i_stream < 15 ; i_stream++ ) {
    i_chem = 0;
    max    = 0.0;
    mtot   = 0.0;
    for ( i = 0 ; i < nb_chem ; i++ ) {
      m = (ARRONDI) ? arrondi ( s[i_stream]->chem[i]->m , 4 ) : s[i_stream]->chem[i]->m;
	
      if ( m > EPS ) {
	mtot += m;
	if ( m > max ) {
	  max    = m;
	  i_chem = i;
	}
      }
    }
    if (mtot > EPS ) {
      purity = max/mtot;
      d = price[i_chem] * max;
      if ( purity < 0.5 )
	d = 0.0;
      else if ( purity < 0.6 )
	d *= 0.05;
      else if ( purity < 0.7 )
	d *= 0.1;
      else if ( purity < 0.8 )
	d *= 0.15;
      else if ( purity < 0.9 )
	d *= 0.3;
      else if ( purity < 0.95 )
	d *= 0.5;
      else if ( purity < 0.98 )
	d *= 0.6;
      else if ( purity < 0.99 )
	d *= 0.8;
      else if ( purity > 0.995 )
	d *= 1.1;
      Rtot += d*25920000;
    }
  }
  
  // cash flows :
  F = new cashflow(15);
  F->set_rates(0.1,0.4);
  F->set_basics(Itot, Coper, Rtot);
  OK = F->run();

  // calculating profitability indicators :
  if(OK) {
    P = new profitability(F);
    OK = P->run(y);
    delete P;
  }

  delete F;

  // arrondis :
  if (ARRONDI) {
    y[12] = round(y[12]);
    y[13] = round(y[13]);
    y[ 9] = round(y[ 9]*10)/10.0;
    y[10] = round(y[10]*10)/10.0;
  }


  g7  = ( y[ 9] < 1e20 ) ? ( y[9] - 4.0 ) / 4.0 : 1e20;
  g8  = ( y[10] < 1e20 ) ? (0.2-y[10])/0.2 : 1e20;
  g9  = ( y[11] < 1e20 ) ? (y[11]-10e6)/10e6 : 1e20;
  g10 = ( y[12] > 0 && y[12] < 1e20 ) ? (y[12]-15e6) / 15e6 : 1e20;

  f   = ( y[13] > 0 && y[13] < 1e20 ) ? -y[13] : 1e20;

//   cout << setprecision(10);
//   cout << "\n\n";
//   for ( int i = 0 ; i < 14 ; i++ )
//     cout << "y[" << i << "] = " << y[i] << endl;


  // menage et affichage du resultat de la boite noire :
 TERMINATE:
  
  //   cout << "\n\n";
  //   cout << " g0 = " << g0  << endl
  //        << " g1 = " << g1  << endl
  //        << " g2 = " << g2  << endl
  //        << " g3 = " << g3  << endl
  //        << " g4 = " << g4  << endl
  //        << " g5 = " << g5  << endl
  //        << " g6 = " << g6  << endl
  //        << " g7 = " << g7  << endl
  //        << " g8 = " << g8  << endl
  //        << " g9 = " << g9  << endl
  //        << "g10 = " << g10 << endl
  //        << "  f = " << f   << endl;


  cout << g0  << " "
       << g1  << " "
       << g2  << " "
       << g3  << " "
       << g4  << " "
       << g5  << " "
       << g6  << " "
       << g7  << " "
       << g8  << " "
       << g9  << " "
       << g10 << " "
       << f   << endl;

  if (units)
    delete units;
  if (safe)
    delete [] safe;

  if (chem) {
    for ( i = 0 ; i < nb_chem ; i++ )
      delete chem[i];
    delete [] chem;
  }

  if (s) {
    for ( i = 0 ; i < nb_s ; i++ )
      delete s[i];
    delete [] s;
  }

  return 0;
}
