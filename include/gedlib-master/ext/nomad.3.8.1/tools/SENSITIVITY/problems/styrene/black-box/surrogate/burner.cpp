#include "burner.hpp"
using namespace std;

/*----------------------------------------------------------*/
/*  arrondi pour ne conserver que n chiffres significatifs  */
/*----------------------------------------------------------*/
double arrondi ( double x , int n ) {
  if (fabs(x) < EPS)
    return 0.0;
  double m = pow ( 10 , ceil(-log10(x)) + n - 1 );
  return round(m*x)/m;
}

/*---------------------------------------------------------------*/
burner::burner ( int nb , chemical ** chem ) {

  rem_nb = nb;

  NO=NO2=N2O=CO=0.0;
  m = new double[nb];
  can_burn = new bool[nb];

  // combustion.prop :
  // -----------------
  // 64-17-5 3 2 3
  // 74-82-8 2 1 2
  // 1333-74-0 0.5 0 1
  // 100-42-5 10 8 4
  // 74-85-1 3 2 2
  // 108-88-3 9 7 4
  // 100-41-4 10.5 8 5
  // 71-43-2 7.5 6 3

  for ( i = 0 ; i < nb ; i++ ) {

    can_burn[i] = false;

    if ( chem[i]->CAS == "64-17-5"   ||
	 chem[i]->CAS == "74-82-8"   ||
	 chem[i]->CAS == "1333-74-0" ||
	 chem[i]->CAS == "100-42-5"  ||
	 chem[i]->CAS == "74-85-1"   ||
	 chem[i]->CAS == "108-88-3"  ||
	 chem[i]->CAS == "100-41-4"  ||
	 chem[i]->CAS == "71-43-2" )
      can_burn[i] = true;
  }

  O2  = new chemical ("7782-44-7");
  N2  = new chemical ("7727-37-9");
  CO2 = new chemical ("124-38-9");
  H2O = new chemical ("7732-18-5");

  // Construct the rx array;
  rx = new combrx * [nb];
  for ( i = 0 ; i < nb ; i++ ) {
    if ( can_burn[i] )
      rx[i] = new combrx ( chem[i]->CAS );
    else
      rx[i] = NULL;
  }
}

/*---------------------------------------------------------------*/
burner::~burner ( void ) {

  delete [] m;
  delete [] can_burn;

  for ( i = 0 ; i < rem_nb ; i++ )
    if (rx[i])
      delete rx[i];
  delete [] rx;
  
  delete O2;
  delete N2;
  delete CO2;
  delete H2O;
}

/*---------------------------------------------------------------*/
bool burner::solve(double * y)
{
   OK=true;
   //perform mass balance (neglect pollutants flows)
   out->m = 0.0;
   for(i=0;i<in->nb;i++)
   {
     if (!can_burn[i]) {
       out->chem[i]->m = in->chem[i]->m;
       out->m+=out->chem[i]->m;
     }
     else {
       out->chem[i]->m=0.0;
       O2->m+=rx[i]->O2_flow()*in->chem[i]->n();
       N2->m+=rx[i]->N2_flow()*in->chem[i]->n();
       CO2->m+=rx[i]->CO2_flow()*in->chem[i]->n();
       H2O->m+=rx[i]->H2O_flow()*in->chem[i]->n();
     }
   }
   N2->m*=(1.0+eta);
   O2->m*=(1.0+eta);
   //perform energy balance to find Tout
   T = in->T;

   step=10;
   Q=1;
   // find temperature
   while (fabs(step)>TOL_BURN && fabs(Q)>TOL_BURN && T<MAX_TEMP)
   {
      T+=step;

      if(T>MAX_TEMP)
	T=MAX_TEMP;
      Q = 0.0;
      for ( i = 0 ; i < in->nb ; i++ )
	Q += in->chem[i]->dH ( in->T , T , in->P ) * in->chem[i]->n();

      for ( i = 0 ; i < in->nb ; i++ )
	if ( can_burn[i] )
	  Q += rx[i]->Hcomb(T) * in->chem[i]->n();

      Q += O2->dH ( 293 , T , in->P ) * O2->n();
      Q += N2->dH ( 293 , T , in->P ) * N2->n();


      if (step/fabs(step)*Q >0)
	step*= -0.1;
      else if (fabs(Q)<10)
	step*=0.25;

   }



   out->set_thermo(in->thermo);
   // out->thermo = in->thermo;

   out->set(in->P, T);
   O2->P=in->P; O2->T=T; O2->state=1; O2->find_v();
   N2->P=in->P; N2->T=T; N2->state=1; N2->find_v();
   CO2->P=in->P; CO2->T=T; CO2->state=1; CO2->find_v();
   H2O->P=in->P; H2O->T=T; H2O->state=1; H2O->find_v();
   //check if mixture can burn
   m_can_burn = 0.0;
   for(i=0;i<in->nb;i++) if(can_burn[i]) m_can_burn+=in->chem[i]->n();
   LFLmix=0.0;
   for(i=0;i<in->nb;i++) if(can_burn[i]) LFLmix+=in->chem[i]->n()/m_can_burn*rx[i]->LFL(in->P,T);
   UFLmix=0.0;
   for(i=0;i<in->nb;i++) if(can_burn[i]) UFLmix+=in->chem[i]->n()/m_can_burn*rx[i]->UFL(in->P,T);
   num = 0.0;
   buff=in->T; in->set(in->P, T);
   for(i=0;i<in->nb;i++) if(can_burn[i]) num+=in->chem[i]->n()/in->n()*in->v;
   in->set(in->P, buff);
   den = O2->v+N2->v+out->v;
   composition = num/den;
   if(!(LFLmix<=composition && composition<=UFLmix) || T==MAX_TEMP)
   {
//       logf.open(MESSAGES,ios::app);
//       logf<<"   --> Warning <--  Mixture in "<<name<<" can't burn (LFL="<<LFLmix
//           <<" UFL="<<UFLmix<<" x="<<composition<<").\n";
//       logf.close();
      T=in->T;
      filename = out->name;
      *out=*in;
      out->set(filename);
      // out->write();  // WRITE TOTO
      OK=false;
   }
   else
   {
      O2->P=in->P; O2->T=T; O2->state=1; O2->find_v();
      N2->P=in->P; N2->T=T; N2->state=1; N2->find_v();
      // out->write(); // WRITE TOTO
   }
   if(OK) //compute the pollutants production
   {
      fill_K_array();
      NO = 1e6*sqrt(K[0]*(N2->n()/den)*(O2->n()/den))*den*0.03/(O2->m+N2->m+out->m+H2O->m+CO2->m);
      N2O = 1e6*K[1]*(N2->n()/den)*sqrt(O2->n()/den)*den*0.044/(O2->m+N2->m+out->m+H2O->m+CO2->m);
      NO2 = 1e6*K[2]*sqrt(N2->n()/den)*(O2->n()/den)*den*0.046/(O2->m+N2->m+out->m+H2O->m+CO2->m);
      CO = 1e6*K[3]*(CO2->n()/den)*den/sqrt(O2->n()/den)*0.028/(O2->m+N2->m+out->m+H2O->m+CO2->m);
   }
//    logf.open(MESSAGES,ios::app);
   if (NO>EPS && NO2>EPS && N2O>EPS) {
     // logf<<"   --> Warning <-- Presence of NOx: "<<(NO+NO2+N2O)<<" ppm in "<<name<<".\n";
     y[7] = NO+NO2+N2O;

     if (ARRONDI)
       y[7] = arrondi ( y[7] , 6 );

   }
   if (CO>EPS) {
     // logf<<"   --> Warning <-- Presence of CO: "<<CO<<" ppm in "<<name<<".\n";
     y[8] = CO;
     if (ARRONDI)
       y[8] = arrondi ( y[8] , 6 );

   }
//    logf.close();
   return OK;
}

void burner::fill_K_array()
{
  a[0]=1.0; a[1]=1.0; a[2]=0.5; a[3]=1.0;
  b[0]=1.0; b[1]=0.5; b[2]=1.0; b[3]=-0.5;
  c[0]=2.0; c[1]=1.0; c[2]=1.0; c[3]=1.0;
  K[0] = exp(-120.27*(173.38-T*0.012)/T);
  K[1] = exp(-120.27*(103.64+T*0.074)/T);
  K[2] = exp(-120.27*(51.96+T*0.061)/T);
  K[3] = exp(-120.27*(283.84-T*0.087)/T);
  for(i=0;i<4;i++)
    K[i]*=pow(1000, c[i]-a[i]-b[i]);
}

void burner::write() {

  cout << setprecision(6);

  cout << "WRITE FILE " << RUNTIME << name << ".unit" << " :\n\tBEGIN\n";
  cout << "\t>>         " << name;
  cout << endl << "\t>>           stream in : " << in->name;
  cout << endl << "\t>>           streams out : " << out->name;
  cout << endl << "\t>>           P = " << in->P << " atm,  T(in) = " << in->T
       << "   T(out) = " << T << " K";
  O2->P = 1;
  O2->T = 293;
  O2->state = 1;
  O2->find_v();
  N2->P=1;
  N2->T=293;
  N2->state=1;
  N2->find_v();
  cout << endl << "\t>>           Required air flow = "
       << (O2->m+N2->m) << " kg/s  (" << (O2->v+N2->v) << " m3/s)";
  O2->P=in->P;
  O2->T=T;
  O2->state=1;
  O2->find_v();
  N2->P=in->P;
  N2->T=T;
  N2->state=1;
  N2->find_v();
  step=(eta*O2->v/(1+eta)+N2->v+H2O->v+CO2->v+out->v);
  cout << endl << "\t>>           Total flue gases  = "
       << (out->m+CO2->m+H2O->m+N2->m+eta*O2->m/(1+eta))
       <<" kg/s  (" << step << " m3/s)";
  cout << "\n\tEND\n\n";
  cost();
}


double burner::get_cost ( void ) {


  O2->P = 1;
  O2->T = 293;
  O2->state = 1;
  O2->find_v();
  N2->P=1;
  N2->T=293;
  N2->state=1;
  N2->find_v();
  O2->P=in->P;
  O2->T=T;
  O2->state=1;
  O2->find_v();
  N2->P=in->P;
  N2->T=T;
  N2->state=1;
  N2->find_v();
  step=(eta*O2->v/(1+eta)+N2->v+H2O->v+CO2->v+out->v);

  buff = 3.1761-0.1373*log10(step) + 0.3414*pow(log10(step),2);
  buff = 2.7*pow(10, buff);
  buff = buff*MS_YEAR/MS_2001;

  return buff;
}


void burner::cost ( void ) {
  cout << "WRITE FILE " << RUNTIME << name << ".cost" << " :\n\tBEGIN\n";
  cout << "\t>>" << get_cost();
  cout << "\n\tEND\n\n";
}
