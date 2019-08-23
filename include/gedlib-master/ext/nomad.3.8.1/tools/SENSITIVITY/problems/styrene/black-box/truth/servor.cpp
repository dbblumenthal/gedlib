#include "servor.hpp"
#include "reactor.cpp"
using namespace std;

/*---------------------------------------------------------------------*/
servor::servor ( int nb_u , int n2 , stream ** streams ) {
  nb   = nb_u;
  nb_s = n2;
  type = new string[nb];
  name = new string[nb];
  s    = streams;


//   for (int i=0; i<nb; i++)
//     {
//       type[i]=new char[31];
//       name[i]=new char[31];
//     }
//   cursor=0;
//   s = s_list;
//   // end = new terminator("\0");
//   mix1=NULL;
//   split1=NULL;
//   flash1=NULL;
//   pump1=NULL;
//   col=NULL;
//   react_pfr=NULL;
//   react_cstr=NULL;
//   heat=NULL;

  burn = new burner ( s[0]->nb , s[0]->chem );
}

/*---------------------------------------------------------------------*/
servor::~servor() {
  delete [] type;
  delete [] name;
  delete burn;
}

/*---------------------------------------------------------------------*/
bool servor::solve_process ( const double * x , double * y ) {


  for ( i = 0 ; i < 8 ; i++ )
    costs[i] = 0.0;

  k    = 0;
  norm = 1.0 / TOL_WEGSTEIN;

  for ( cursor = 0 ; cursor < nb ; cursor++ ) {

    if (type[cursor] == "mix" ) {
      do_mix_process(x);
    }
    else if ( type[cursor] == "split" ) {
      do_split_process(x);
    }
    else if ( type[cursor] == "flash" ) {
      do_flash_process(x);
    }
    else if ( type[cursor] == "pump" ) {
      do_pump_process(x);
    }
    else if ( type[cursor] == "heatx") {
      do_heatx_process(x);
    }
    else if ( type[cursor] == "burner" ) {
      do_burner_process(x,y);
    }
    else if ( type[cursor] == "column" ) {
      do_column_process(x,y);
    }
    else if ( type[cursor] == "reactor" ) {
      do_reactor_process(x);
    }
    else if ( type[cursor] == "loop" ) {

      recycle     = 10;
      end_recycle = 0;
      do_loop_process(x);
    }
    else {
      cout << "ERROR 18\n\n";
      exit(0);
    }

  }

  return true;
}


/*---------------------------------------------------------------------*/
void servor::do_loop_process ( const double * x ) {

  // structure in the input file :
  // loop
  // name                              : "looping"
  // index of recycle stream           : 11
  // index of stream's beginning block : 7
  // index of stream's end block       : 1

  // TOTO
//   if (k==0)
//     cout << endl << "       -> Wegstein iterations "; 
//   else if ( k <= 2 )
//     cout << endl << "       -> loop " << setw(3) << k;
//   else
//     cout << endl << "       -> loop " << setw(3) << k << "  > error " << norm;

  // Get the two fisrst points
  if ( k==0 ) {
    for ( i = 0 ; i < s[0]->nb ; i++ )
      x_last[i] = s[recycle]->chem[i]->m; 
    end_loop = false;
  }

  if ( k==1 ) {
    for ( i = 0 ; i < s[0]->nb ; i++ ) {
      x_now [i] = s[recycle]->chem[i]->m;
      g_last[i] = s[recycle]->chem[i]->m;
    }

    end_loop = false;
  }

  if ( k == 2 ) {
    for ( i = 0 ; i < s[0]->nb ; i++ )
      g_now[i] = s[recycle]->chem[i]->m;
    end_loop=false;
  }
  k++;


  // run the Wegstein algorithm
  if ( k > 2 ) {
    for ( i = 0 ; i < s[0]->nb ; i++ ) {
      g_now[i] = s[recycle]->chem[i]->m;

      if ( fabs (x_now[i]-x_last[i]) > EPS )
	slope[i] = ( g_now[i] - g_last[i] ) / ( x_now[i] - x_last[i] );
      else
	slope[i] = 0;

      theta[i] = 1.0 / (1.0-slope[i]);
      if ( theta[i] < MIN_THETA )
	theta[i] = MIN_THETA;
      if ( theta[i] > MAX_THETA )
	theta[i] = MAX_THETA;
      x_next[i] = (1.0-theta[i])*x_now[i] + theta[i]*g_now[i];
    }
    norm = 0.0;
    for ( i = 0 ; i < s[0]->nb ; i++ )
      if ( fabs(x_now[i]) > EPS )
	norm += fabs (x_next[i]-x_now[i]) / fabs(x_now[i]);

    if ( norm > TOL_WEGSTEIN && k < MAX_ITER_WEGSTEIN ) {
      s[recycle]->m = 0.0;
      for ( i = 0 ; i < s[0]->nb ; i++ ) {
	s[recycle]->chem[i]->m = x_next[i];
	s[recycle]->m += x_next[i];
	x_last[i] = x_now[i];
	g_last[i] = g_now[i];
	x_now[i] = x_next[i];
      }
      end_loop=false;
    }
    if ( norm <= TOL_WEGSTEIN && k < MAX_ITER_WEGSTEIN )
      end_loop = true;

  }

 
  if ( end_loop ) {
    if ( k < MAX_ITER_WEGSTEIN && k > 3 )  {
      // cout<<" OK"; // TOTO
      // s[recycle]->write(); // WRITE TOTO

//       // WRITE TOTO :
//       cout << "WRITE FILE " << RUNTIME << name[cursor] << ".unit" << " :\n\tBEGIN\n";
//       cout << "\t>>         " << name[cursor];
//       cout << endl << "\t>>           from block " << cursor+1 << " to block " << end_recycle+1;
//       cout << endl << "\t>>           Wegstein converged in "
// 	   << k << " iterations (rel. err. " << norm << ").";
//       cout << "\n\tEND\n\n";

      k = 0;
      norm = 1.0/TOL_WEGSTEIN;
    }
//    else {
//       log.open(MESSAGES, ios::app);
//       log<<"   ==> Error <== Wegstein algorithm did not converge.\n";
//       log.close();      
//    }      
  }


  if ( !end_loop && k < MAX_ITER_WEGSTEIN )
    cursor = end_recycle-1;
  if ( !end_loop && k==MAX_ITER_WEGSTEIN ) {
//       log.open(MESSAGES, ios::app);
//       log<<"   ==> Error <== Wegstein algorithm did not converge.\n";
//       log.close();
    k=0;
    norm = 1.0/TOL_WEGSTEIN;
    cursor=nb;
  }

}

/*---------------------------------------------------------------------*/
void servor::do_burner_process ( const double * x , double * y ) {

  // cout << endl << "         -- " << name[cursor] << "... "; // TOTO

  // read parameters :
  int    i1 = 8;
  int    i2 = 13;
  double f  = x[6];

  burn->set ( s[i1-1] , s[i2-1] );
  burn->set(f);
  burn->set(name[cursor]);

  if ( burn->solve(y) ) {
    // cout << "OK"; // TOTO
    // burn->write(); // WRITE TOTO
    costs[7] = burn->get_cost();
  }
  else {
    cout << "ERROR 20\n\n";
    exit(0);
  }
}

/*---------------------------------------------------------------------*/
void servor::do_split_process ( const double * x )
{
  // cout<<endl<<"         -- "<<name[cursor]<<"... ";  // TOTO

  //Read parameters
  int i1 = 9;
  int i2 = 2;
  int i_tab[2] = { 10 , 11 };
  double f_tab[2];
  f_tab[0] = x[5];
  f_tab[1] = 1-x[5];

  stream * list1[2];
  list1[0] = s[i_tab[0]-1];
  list1[1] = s[i_tab[1]-1];

  split * split1 = new split ( i2 , s[i1-1] , list1 );

  split1->set(f_tab);
  split1->set(name[cursor]);
  if ( split1->solve() ) {
    // cout<<"OK"; // TOTO
    // split1->write(); // WRITE TOTO
  }
  else {
    cout << "ERROR 19\n\n";
    exit(0);
  }
  delete split1;
}

/*---------------------------------------------------------------------*/
void servor::do_column_process ( const double * x , double * y ) {
  // cout<<endl<<"         -- "<<name[cursor]<<"... "; // TOTO
  
  //Read parameters
  double f1 , f2;
  int    i , i1 , i2;
  int    i_tab[2];
  double f = 1.0;

  if (name[cursor]=="sep-sty") {
    i        =  7;
    i1       = 15;
    i2       =  9;
    i_tab[0] =  1;
    i_tab[1] =  7;
    f1 = f2 = x[2];
  }
  else if (name[cursor]=="sep-bz") {
    i        = 10;
    i1       = 12;
    i2       = 14;
    i_tab[0] =  5;
    i_tab[1] =  1;
    f1 = f2 = x[3];
  }
  else {
    cout << "ERROR 17\n\n";
    exit(0);
  }

  column * col = new column ( s[i-1] , s[i1-1] , s[i2-1] );
  col->set ( f , i_tab[0] , f1 , i_tab[1] , f2 );
  col->set(name[cursor]);
  if ( col->solve() ) {
    //cout<<"OK"; // TOTO
    //col->write(); // WRITE TOTO

    if (name[cursor]=="sep-sty") {
      y[4] = col->get_N();
      costs[5] = col->get_cost();
      power[5] = col->get_power();
      water[5] = col->get_water();
    }
    else {
      y[5] = col->get_N();
      costs[6] = col->get_cost();
      power[4] = col->get_power();
      water[4] = col->get_water();
    }
  }
  else {
    cout << "ERROR 15\n\n";
    exit(0);
  }
  delete col;
}

/*---------------------------------------------------------------------*/
void servor::do_flash_process ( const double * x ) {
  // cout<<endl<<"         -- "<<name[cursor]<<"... ";  // TOTO

  //Read parameters
  double f1        = 1.0;
  double f2        = x[7];
  // int    i_tab [3] = { 6 , 7 , 8 };
  // flash * flash1 = new flash ( s[i_tab[0]-1] , s[i_tab[1]-1] , s[i_tab[2]-1] );
  flash * flash1 = new flash ( s[5] , s[6] , s[7] );

  flash1->set(f1, f2);
  flash1->set(name[cursor]);
  if ( flash1->solve() ) {
    // cout<<"OK";  // TOTO
    // flash1->write(); // WRITE TOTO
    costs[4] = flash1->get_cost();
    power[2] = flash1->get_power();
    water[1] = flash1->get_water();
  }
  else {
    cout << "ERROR 14\n\n";
    exit(0);
  }
  delete flash1;
}


/*---------------------------------------------------------------------*/
void servor::do_mix_process ( const double * x ) {

  // cout << endl << "         -- " << name[cursor] << "... "; // TOTO

  // Read parameters (hardcode avec eb2sty.process) :
  double f1       = 1.0;
  int    i1       = 3;
  int    i_tab[3] = { 1 , 12 , 11 };
  int    i2       = 2;

  // We can solve the unit
  stream * list2 = s[i2-1] , ** list1 = new stream * [i1];
  for ( int i = 0 ; i < i1 ; i++ )
    list1[i] = s[i_tab[i]-1];
  mix * mix1 = new mix ( i1 , list1 , list2 );
  mix1->set(f1);
  mix1->set(name[cursor]);
  if (mix1->solve()) {
    // cout<<"OK"; // TOTO
    // mix1->write(); // WRITE TOTO
  }
  else {
    cout << "ERROR 6\n\n";
    exit(0);
  }
  delete mix1;
  delete [] list1;
}


/*---------------------------------------------------------------------*/
void servor::do_pump_process ( const double * x ) {

  // cout<<endl<<"         -- "<<name[cursor]<<"... ";  // TOTO

  //Read parameters (hardcode avec eb2sty.process) :
  int    i1 = 2;
  int    i2 = 3;
  double f1 = x[4];
  double f2 = 0.75;

  pump * pump1 = new pump ( s[i1-1] , s[i2-1] );
  pump1->set(f1,f2);
  pump1->set(name[cursor]);



  // solve :
  if ( pump1->solve() ) {
    // cout<<"OK";  // TOTO
    // pump1->write(); // WRITE TOTO

    power[0] = pump1->get_power();
    costs[0] = pump1->get_cost();
  }
  else {
    cout << "ERROR 7\n\n";
    exit(0);
  }

  delete pump1;
}

/*---------------------------------------------------------------------*/
void servor::do_heatx_process ( const double * x ) {

  // cout << endl << "         -- " << name[cursor] << "... "; // TOTO

  //Read parameters (idem) :

  bool   b = false;

  // heater :
  int    i1;
  int    i2;
  double f1;

  // heater :
  if (name[cursor]=="heater") {
    i1 = 3;
    i2 = 4;
    f1 = x[0];
  }

  // cooler :
  else if (name[cursor]=="cooler") {
    i1 = 5;
    i2 = 6;
    f1 = x[7];
  }
  else {
    cout << "ERROR 16\n\n";
    exit(0);
  }

  double f2 = 0.85;

  // solve :
  heatx * heat = new heatx ( b , s[i1-1] , s[i2-1] );
  heat->set(f1,f2);
  heat->set(name[cursor]);
  if (heat->solve()) {
    // cout<<"OK"; // TOTO
    // heat->write();  // WRITE TOTO
    if (name[cursor]=="heater") {
      costs[1] = heat->get_cost();
      power[3] = heat->get_power();
      water[2] = heat->get_water();
    }
    else {
      costs[3] = heat->get_cost();
      power[1] = heat->get_power();
      water[0] = heat->get_water();
    }
  }
  else {
    cout << "ERROR 8\n\n";
    exit(0);
  }
  delete heat;
}

/*---------------------------------------------------------------------*/
void servor::do_reactor_process ( const double * x ) {

  // cout<<endl<<"         -- "<<name[cursor]<<"... "; // TOTO

  // Read parameters (idem) :
  int    i1   = 4;
  int    i2   = 5;

  reactor<pfr> * react_pfr = new reactor<pfr> ( s[i1-1] , s[i2-1] );

  react_pfr->set(name[cursor]);

  double f1 = x[1];
  double f2 = 0.5;

  string list[5] = { "eb2sty" , "sty2eb" , "eb2bz" , "eb2tol" , "tol2bz" };

  i1 = 5;

  react_pfr->set(f1,f2,i1,list);

  f1 = 0.0;
  f2 = 300.0;

  react_pfr->set(f1,f2);

  if ( react_pfr->solve() ) {
    // cout<<"OK"; // TOTO
    // react_pfr->write(); // WRITE TOTO
    costs[2] = react_pfr->get_cost();
    water[3] = react_pfr->get_water();
  }
  else {
    cout << "ERROR 9\n\n";
    exit(0);
  }
  delete react_pfr;
}


double  servor::get_costs_sum ( void ) const {
  double sum = 0.0;
  for ( int i = 0 ; i < 8 ; i++ )
    sum += ( (ARRONDI) ? arrondi(costs[i],6) : costs[i] );
  return sum;
}

double  servor::get_power_sum ( void ) const {
  double sum = 0.0;
  for ( int i = 0 ; i < 6 ; i++ )
    sum += ( (ARRONDI) ? arrondi(power[i],6) : power[i] );
  return sum;
}

double servor::get_water_sum ( void ) const {
  double sum = 0.0;
  for ( int i = 0 ; i < 6 ; i++ )
    sum += ( (ARRONDI) ? arrondi(water[i],6) : water[i] );
  return sum;
}
