#include "chemical.hpp"
using namespace std;

void chemical::check_error ( void ) {
  if (error>MAX_ERROR) {
    cout << "ERROR 2\n\n";
    exit(0);
  }  
  if (warning>MAX_WARNING) {
    cout << "ERROR 3\n\n";
    exit(0);
  }
}

// copy-constr. :
chemical::chemical ( const chemical & chem ) {

  CAS           = chem.CAS;
  name          = chem.name;
  M             = chem.M;

  state         = chem.state;
  Tm            = chem.Tm;
  Tb            = chem.Tb;
  Tc            = chem.Tc;
  Pc            = chem.Pc;
  Ho            = chem.Ho;
  rho_liq       = chem.rho_liq;
  dHvap         = chem.dHvap;

  mu_param[0]   = chem.mu_param[0];
  mu_param[1]   = chem.mu_param[1];
  Cp_param[0]   = chem.Cp_param[0];
  Cp_param[1]   = chem.Cp_param[1];
  Cp_param[2]   = chem.Cp_param[2];
  Cp_param[3]   = chem.Cp_param[3];
  Cp_liq        = chem.Cp_liq;
  Psat_param[0] = chem.Psat_param[0];
  Psat_param[1] = chem.Psat_param[1];
  Psat_param[2] = chem.Psat_param[2];
    
  thermo = new thermolib();
  thermo->send(Pc,Tc, omega());

  P = chem.P;
  T = chem.T;
  m = chem.m;
  v = chem.v;

  warning = chem.warning;
  error   = chem.error;
  tmp     = chem.tmp;
}

chemical::chemical ( const string & chem_name ) {

  CAS = chem_name;
	
	// C. Tribes add initialization for more robustness (variables may be initialized uncorrectly dependent on the execution)
	P=T=m=v=0.0;	

  // 1/12 :
  if (CAS=="100-41-4") {
    name = "ethylbenzene";
    M = 106.17;
    state = 0;
    Tm = 178.2;
    Tb = 409.3;
    Tc = 617.1;
    Pc = 35.6;
    Ho = 29.79;
    rho_liq = 867.0;
    dHvap = 35.56;
    mu_param[0] = 472.82;
    mu_param[1] = 264.22;
    Cp_param[0] = -43.069;
    Cp_param[1] = 7.067e-01;
    Cp_param[2] = -4.807e-04;
    Cp_param[3] = 1.30e-07;
    Cp_liq = 190.23;
    Psat_param[0] = 16.0195;
    Psat_param[1] = 3279.47;
    Psat_param[2] = -59.95;
  }

  // 2/12 :
  else if (CAS=="1333-74-0") {
    name = "hydrogen";
    M = 2.02;
    state = 1;
    Tm = 14.0;
    Tb = 20.4;
    Tc = 33.2;
    Pc = 12.8;
    Ho = 0.0;
    rho_liq = 71.0;
    dHvap = 0.9;
    mu_param[0] = 13.82;
    mu_param[1] = 5.39;
    Cp_param[0] = 27.124;
    Cp_param[1] = 9.267e-03;
    Cp_param[2] = -1.380e-05;
    Cp_param[3] = 7.64e-09;
    Cp_liq = 0.0;
    Psat_param[0] = 13.6333;
    Psat_param[1] = 164.90;
    Psat_param[2] = 3.19;
  }

  // 3/12 :
  else if (CAS=="108-88-3") {
    name = "toluene";
    M =92.14;
    state = 0;
    Tm = 178.0;
    Tb = 383.8;
    Tc = 591.7;
    Pc = 40.6;
    Ho = 50.0;
    rho_liq = 867;
    dHvap = 33.18;
    mu_param[0] = 467.33;
    mu_param[1] = 255.24;
    Cp_param[0] = -24.338;
    Cp_param[1] = 5.121e-1;
    Cp_param[2] = -2.763e-4;
    Cp_param[3] = 4.91e-8;
    Cp_liq      = 159.85;
    Psat_param[0] = 16.0137;
    Psat_param[1] = 3096.52;
    Psat_param[2] = -53.67;
  }

  // 4/12 :
  else if (CAS=="74-82-8") {
    name = "methane";
    M =16.04;
    state =1;
    Tm = 90.7;
    Tb = 111.7;
    Tc = 190.6;
    Pc = 45.4;
    Ho = -74.85;
    rho_liq = 425;
    dHvap = 8.18;
    mu_param[0] = 114.14;
    mu_param[1] = 57.60;
    Cp_param[0] = 19.238;
    Cp_param[1] = 5.209e-02;
    Cp_param[2] = 1.197e-05;
    Cp_param[3] = -1.13e-08;
    Cp_liq      = 0.0;
    Psat_param[0] = 15.2243;
    Psat_param[1] = 897.84;
    Psat_param[2] = -7.16;
  }

  // 5/12 :
  else if (CAS=="71-43-2") {
    name = "benzene";
    M = 78.11;
    state = 0;
    Tm = 278.7;
    Tb = 353.3;
    Tc = 562.1;
    Pc = 48.3;
    Ho = 82.93;
    rho_liq = 885;
    dHvap = 30.76;
    mu_param[0] = 545.64;
    mu_param[1] = 265.24;
    Cp_param[0] = 33.894;
    Cp_param[1] = 4.74e-1;
    Cp_param[2] = -3.015e-4;
    Cp_param[3] = 7.13e-8;
    Cp_liq      = 116.03;
    Psat_param[0] = 15.9008;
    Psat_param[1] = 2788.51;
    Psat_param[2] = -52.36;
  }

  // 6/12 :
  else if (CAS=="74-85-1") {
    name = "ethylene";
    M = 28.05;
    state =1;
    Tm = 104.0;
    Tb = 169.4;
    Tc = 282.4;
    Pc = 49.7;
    Ho =52.3;
    rho_liq = 577.0;
    dHvap = 13.54;
    mu_param[0] = 168.98;
    mu_param[1] = 93.94;
    Cp_param[0] = 3.803;
    Cp_param[1] = 1.565e-01;
    Cp_param[2] = -8.343e-05;
    Cp_param[3] = 1.75e-08;
    Cp_liq      = 0.0;
    Psat_param[0] =15.5368;
    Psat_param[1] = 1347.01;
    Psat_param[2] = -18.15;
  }

  // 7/12 :
  else if (CAS=="100-42-5") {
    name = "styrene";
    M =104.15;
    state = 0;
    Tm = 242.5;
    Tb = 418.3;
    Tc =647.0;
    Pc =39.4;
    Ho =147.36;
    rho_liq =906.0;
    dHvap = 36.82;
    mu_param[0] = 528.64;
    mu_param[1] = 276.71;
    Cp_param[0] =-28.229;
    Cp_param[1] =6.155e-01;
    Cp_param[2] = -4.020e-04;
    Cp_param[3] = 9.93e-08;
    Cp_liq      =166.13;
    Psat_param[0] = 16.0193;
    Psat_param[1] = 3328.57;
    Psat_param[2] =-63.72;
  }

  // 8/12 :
  else if (CAS=="7782-44-7") {
    name = "oxygen";
    M = 32.00;
    state = 1;
    Tm = 54.4;
    Tb = 90.2;
    Tc = 154.6;
    Pc = 49.8;
    Ho =0.0 ;
    rho_liq =1149.1 ;
    dHvap =6.82 ;
    mu_param[0] = 85.68;
    mu_param[1] =  51.50;
    Cp_param[0] = 28.087;
    Cp_param[1] = -3.678e-06 ;
    Cp_param[2] = 1.745e-05;
    Cp_param[3] = -1.06e-08;
    Cp_liq      =0.0 ;
    Psat_param[0] = 15.4075;
    Psat_param[1] =  734.55 ;
    Psat_param[2] =-6.45 ;
  }

  // 9/12 :
  else if (CAS=="7727-37-9") {
    name = "nitrogen";
    M = 28.01;
    state = 1;
    Tm = 63.3;
    Tb = 77.4;
    Tc =  126.2;
    Pc =  33.5;
    Ho = 0.0;
    rho_liq = 804.0;
    dHvap = 5.58;
    mu_param[0] = 90.30;
    mu_param[1] = 46.41;
    Cp_param[0] = 31.128;
    Cp_param[1] = -1.356e-02 ;
    Cp_param[2] = 2.678e-05;
    Cp_param[3] =-1.17e-08 ;
    Cp_liq      = 0.0;
    Psat_param[0] = 14.9342;
    Psat_param[1] =  588.72;
    Psat_param[2] = -6.60;
  }

  // 10/12 :
  else if (CAS=="124-38-9") {
    name = "carbon-dioxide";
    M =44.01;
    state = 1;
    Tm = 216.6;
    Tb = 194.4;
    Tc = 304.2;
    Pc =  72.8;
    Ho =  -393.41;
    rho_liq = 777.0;
    dHvap = 17.15;
    mu_param[0] = 578.08;
    mu_param[1] = 185.24 ;
    Cp_param[0] = 19.782;
    Cp_param[1] = 7.339e-02;
    Cp_param[2] = -5.598e-05;
    Cp_param[3] = 1.71e-08;
    Cp_liq      = 0.0;
    Psat_param[0] = 22.5898;
    Psat_param[1] =3103.39 ;
    Psat_param[2] =  -0.16;
  }

  // 11/12 :
  else if (CAS=="7732-18-5") {
    name = "water";
    M =18.02;
    state = 0;
    Tm = 273.15;
    Tb = 373.15;
    Tc = 647.4;
    Pc = 217.6;
    Ho = -241.83;
    rho_liq = 998 ;
    dHvap = 40.66;
    mu_param[0] = 658.25;
    mu_param[1] = 283.16;
    Cp_param[0] = 32.220;
    Cp_param[1] = 1.923e-03 ;
    Cp_param[2] = 1.055e-05;
    Cp_param[3] = -3.59e-09;
    Cp_liq      = 75.24;
    Psat_param[0] = 18.3036;
    Psat_param[1] = 3816.44;
    Psat_param[2] = -46.13;
  }

  // 12/12 :
  else if (CAS=="64-17-5") {
    name = "ethanol";
    M =46.07;
    state =0 ;
    Tm =159.1 ;
    Tb = 351.5;
    Tc =  516.2;
    Pc =63.0 ;
    Ho =  -234.8;
    rho_liq = 789.0;
    dHvap =38.74 ;
    mu_param[0] = 686.64;
    mu_param[1] = 300.88;
    Cp_param[0] = 9.008;
    Cp_param[1] =  2.139e-01;
    Cp_param[2] = -8.385e-05 ;
    Cp_param[3] =  1.37e-09;
    Cp_liq      = 2.22;
    Psat_param[0] = 18.9119;
    Psat_param[1] =  3803.98;
    Psat_param[2] = -41.68;
  }

  else {
    cout << "ERROR 1\n\n";
    exit(0);
  }

  thermo = new thermolib();
  thermo->send(Pc,Tc, omega());


}

double chemical::K()
{
      thermo->set(P,T,v,n());
      return thermo->K();
}

double chemical::mu()
{
	// Returns the fluid's viscosity in Pa.s
   if (Tm<=T && T<=Tboil(P))
      return pow(10,(mu_param[0]*(1.0/T-1.0/mu_param[1])-3));
   else
   {
	   ofstream logf;
      logf.open(MESSAGES, ios::app);
      logf<<"   --> Warning <--  Cannot compute viscosity of "<<name<<".\n";
      logf.close();
      warning++;
      check_error();
      return 1.0;
   }
}

double chemical::rho()
{
	// Returns the fluid's density in kg/m3, wether it's liquid or gas
   if(state==0) tmp= rho_liq;
   if(state==1)
   {
      find_v();
      if (v>EPS) tmp= m/v;
      else tmp= 0.0;
   }
   return  tmp;
}

double chemical::Cp() {
  // cout<<endl<<"Cp de "<<name<<" a "<<T;
  // Returns the fluid's Cp in J/mol.K
  if(state==0) {

    // tmp = Cp_liq;  // BUG : boucle infinie !!!
    return Cp_liq; // SEB
    
  }
  if(state==1 || T>Tboil(P)) {
    tmp=0;
    for (int i=0;i<4;i++)
      tmp+=Cp_param[i]*pow(T,i);
  }
  else {
    T=Tb;

    tmp = Cp(); // ici boucle infinie si state==0 !!!

  }
  return tmp;
}

double chemical::Cp(bool q)
{
	// Returns the fluid's Cp in J/mol.K
   if(q==0) tmp=Cp_liq;
   if(q==1)
   {
   	   tmp=0;
      for (int i=0;i<4;i++) tmp+=Cp_param[i]*pow(T,i);
   }
   return tmp;
}

double chemical::Psat()
{
   // Returns the fluid's vapor pressure in atm, using Antoine's equation
   if(Tm<=T && T<=Tc)
      return (exp(Psat_param[0]-Psat_param[1]/(T+Psat_param[2]))/760.01);
   else
   {

      return Psat(Tb);
   }
}
double chemical::Psat(double t)
{
   // Returns the fluid's vapor pressure in atm, using Antoine's equation
      return (exp(Psat_param[0]-Psat_param[1]/(t+Psat_param[2]))/760.01);
}

double chemical::dH(double T1,double T2, double pres)
{
   //Enthalpy variation in kJ/mol. Does not affect any attributes of current object.
   double energy=0, TT=Tboil(pres), vap=Hvap(TT);
   int sign=1, i;
   if (T2<T1) {sign = -1; energy=T1; T1=T2; T2=energy; energy=0;}
   if (T1==T2) energy = 0.0;
   if (T2<TT) energy = Cp_liq*(T2-T1)/1000;
   if (TT<T1) for (i=1;i<=4;i++) energy+=Cp_param[i-1]*(pow(T2,i)-pow(T1,i))/i/1000;
   if(T1<=TT && TT<=T2)
   {
      energy=Cp_liq*(TT-T1)/1000;
      energy+=vap;
      for (i=1;i<=4;i++) energy+=Cp_param[i-1]*(pow(T2,i)-pow(TT,i))/i/1000;
   }
   return energy*sign;
}

void chemical::find_T()
{
   if(n()>EPS && P>EPS)
   {
      thermo->set(P,T,v,n());
      T=thermo->T();
   }
   else
   {
	   ofstream logf;
      logf.open(MESSAGES, ios::app);
      logf<<"   --> Warning <--  Cannot find T of "<<name<<".\n";
      logf.close();
      warning++;
   }
   check_error();
}

void chemical::find_P()
{
   if(n()>EPS && T>EPS)
   {
      thermo->set(P,T,v,n());
      P=thermo->P();
   }
   else
   {
	   ofstream logf;
      logf.open(MESSAGES, ios::app);
      logf<<"   --> Warning <--  Cannot find P of "<<name<<".\n";
      logf.close();
      warning++;
   }
   check_error();
}

void chemical::find_v()
{

   if(state==0) v=m/rho_liq;
   if(state==1 && P>EPS && T>EPS && m>EPS)
   {
      thermo->set(P,T,v,n());
      v=thermo->v();
   }
}

void chemical::find_state()
{
   ofstream logf;
   if (T>Tc || P>Pc) state = 1;      //T or P is bigger than Tc or Pc
   if (T<=Tm)                         //T is smaller than melting point
   {
	   ofstream logf;
      logf.open(MESSAGES, ios::app);
      logf<<"   --> Warning <--  The chemical "<<name<<" is solid.\n";
      logf.close();
      warning++;
   }
   check_error();
   if (T<Tboil(P)) state=0;
   else state=1;
}
