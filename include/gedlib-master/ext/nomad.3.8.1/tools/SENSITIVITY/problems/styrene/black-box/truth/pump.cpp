#include "pump.hpp"
using namespace std;

bool pump::solve() {

  // out->chem = in->chem;
  out->set ( in->nb , in->chem );
  out->m = in->m;


  in->set(in->P,in->T);
  for ( i = 0 ; i < in->nb ; i++ )
    if(in->chem[i]->m>EPS) {
      in->chem[i]->find_v();
      if(in->chem[i]->state==1) {
	W+=in->chem[i]->gamma()*in->T*0.0083144*in->chem[i]->n()/
	  (in->chem[i]->gamma()-1.0)*(pow(P/in->P, 1.0-1.0/in->chem[i]->gamma())-1.0);
	tmp += in->chem[i]->gamma();
	n++;
      }
      if(in->chem[i]->state==0)
	W+=in->chem[i]->v*(P-in->P)*101.325;
    }
  if (fabs(state-1)<EPS) //compressing gases
    out->T = in->T*pow(P/in->P, 1.0-1.0/(tmp/n));
  else //compressing liquids
    out->T=in->T;
  out->set(P, out->T);
  if(eta>EPS)
    W /= eta;
  else
    success=false;
  // out->write(); // WRITE TOTO
  return success;
}

void pump::write() {

  cout << setprecision(6);

  string file_name = RUNTIME + name + ".unit";
  cout << "WRITE FILE " << file_name << " :\n\tBEGIN\n";
  cout <<"\t>>         " << name;
  cout << endl << "\t>>           stream in: "<<in->name<<"   out: "<<out->name;
  cout << endl << "\t>>           P(in) = "<<in->P<<"   P(out) = "<<out->P<<"  atm";
  cout << endl << "\t>>           T(in) = "<<in->T<<"   T(out) = "<<out->T<<"  K";
  cout << endl << "\t>>           Shaft work = "<<W;
  if (success)
    cout <<" kW (converge normally)";
  cout << "\n\tEND\n\n";

  power();
  cost();
}


double pump::get_cost ( void ) {

  if ( fabs(state-1) < EPS ) {
    if(W<450) W=450; if(W>3000)W=3000;
    tmp=2.2891+1.3604*log10(W)-0.1027*pow(log10(W),2);
    tmp=3.2*pow(10.0, tmp);
    tmp1=2.4604+1.4191*log10(W)-0.1798*pow(log10(W),2);
    tmp1=1.5*pow(10.0, tmp1);
    tmp+=tmp1;
  }
  else {
    if(W<1) W=1; if(W>300)W=300;
    tmp=3.3892+0.0536*log10(W)+0.1538*pow(log10(W),2);
    tmp=pow(10.0, tmp);
    P=(P-1.0)*101.325/100.0;
    if (P<EPS) P=1; if(P>100) P=100;
    W = -0.3925+0.3957*log10(P)-0.00226*pow(log10(P),2);
    W=pow(10.0, W); if(W<1) W=1;
    tmp*=(1.89+1.35*W*1.8);
  }
  tmp = tmp*MS_YEAR/MS_2001;
  return tmp;
}


void pump::cost() {
  string file_name = RUNTIME + name + ".cost";
  cout << "WRITE FILE " << file_name << " :\n\tBEGIN\n";
  cout << "\t>>" << get_cost();
  cout << "\n\tEND\n\n";
}

void pump::power() {
  string file_name = RUNTIME + name + ".power";
  cout << "WRITE FILE " << file_name << " :\n\tBEGIN\n";
  cout << "\t>>" << W;
  cout << "\n\tEND\n\n";
}
