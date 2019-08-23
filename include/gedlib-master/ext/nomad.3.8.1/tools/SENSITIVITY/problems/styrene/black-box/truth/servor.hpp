#ifndef SERVOR_H
#define SERVOR_H

//Possible units
#include "mix.hpp"
#include "split.hpp"
#include "pump.hpp"
#include "column.hpp"
#include "reactor.hpp"
#include "heatx.hpp"
#include "burner.hpp"

using namespace std;

class servor {

public:

  double    norm;
  int       k;
  int       nb;
  int       nb_s;
  string  * type;
  string  * name;
  int       cursor;
  stream ** s;

  int    recycle;
  int    end_recycle;
  int    i;
  bool   end_loop;

  double x_last[MAX_DIM];
  double x_now [MAX_DIM];
  double x_next[MAX_DIM];
  double g_last[MAX_DIM];
  double g_now [MAX_DIM];
  double slope [MAX_DIM];
  double theta [MAX_DIM];

  double costs[8];
  double power[6];
  double water[6];


  // clock_t beg_time, end_time;
  //   bool OK, b, solve_OK;
  //   char filename[31], kind[10], **list, process_name[31];
  //   double f, f1, f2, f_tab[10], t;
  //   int i, i1, i2, i_tab[10], time;
  //   stream *list2;

//   double theta[MAX_DIM], slope[MAX_DIM], norm;
//   bool end_loop;
//   int k, recycle, end_recycle;
//   char loop_name[31];
 
//   mix *mix1; void do_mix();
//   split *split1;
//   flash *flash1; void do_flash();
//   pump *pump1; void do_pump();
//   column *col; void do_column();
//   reactor<pfr> *react_pfr; void do_reactor();
//   reactor<cstr> *react_cstr;
//   heatx *heat; void do_heatx();
  burner * burn;
	    
//   stream *s; //list of streams
//   friend void out_of_memory(void){exit(0);}


  // constructeur :
  servor ( int , int , stream ** );

  // destructeur :
  ~servor();


//   void set(int t) {time=t;}
//   void set(char n[31]) {strcpy(process_name, n);}
  void do_split_process   ( const double * x );
  void do_column_process  ( const double * x , double * y );
  void do_flash_process   ( const double * x );
  void do_mix_process     ( const double * x );
  void do_pump_process    ( const double * x );
  void do_heatx_process   ( const double * x );
  void do_reactor_process ( const double * x );
  void do_burner_process  ( const double * x , double * y );
  void do_loop_process    ( const double * x );

  bool solve_process      ( const double * x , double * y );		//main solver of the software.	  

  double get_costs_sum ( void ) const;

  double get_power_sum ( void ) const;

  double get_water_sum ( void ) const;

  
};
#endif
