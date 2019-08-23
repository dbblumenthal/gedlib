#include "thermolib.hpp"
#include "defines.hpp"
#include "secant.cpp"
using namespace std;




// affectation :
thermolib & thermolib::operator = ( const thermolib & t ) {

  if (dim!=t.dim)
    reset(t.dim);
  
  for ( i = 0 ; i < dim ; i++ ) {
    Pc[i] = t.Pc[i];
    Tc[i] = t.Tc[i];
    omega[i] = t.omega[i];
    molefrac[i] = t.molefrac[i];
  }

  return *this;
}

thermolib::~thermolib() {
  delete solver;
  delete [] Pc;
  delete [] Tc;
  delete [] omega;
  delete [] molefrac;
}

void thermolib::construct ( void ) {
  molefrac = new double[dim];
  Pc = new double[dim];
  Tc = new double[dim];
  omega = new double[dim];
	
	// C. Tribes add this for more robustness (variables may be initialized uncorrectly dependent on the execution)
	for ( i = 0 ; i < dim ; i++ ) 
	{
		Pc[i] = 0.0;
		Tc[i] = 0.0;
		omega[i] = 0.0;
		molefrac[i] = 0.0;
  }
	
  solver = new secant<thermolib>();
}

void thermolib::reset(int b)
{
   delete [] molefrac;
   delete [] Pc;
   delete [] Tc;
   delete [] omega;
   delete solver;
   dim = b;
   construct();
}

double thermolib::a_mix()
{
   if (dim>1)
   {
      tmp=0;
      for (i=0;i<dim;i++)
         for (j=0;j<dim;j++)
            tmp += molefrac[i]*molefrac[j]*sqrt(a(i)*a(j));
      return tmp;
   }
   else return a(0);
}

double thermolib::b_mix()
{
   if (dim>1)
   {
      tmp=0;
      for (i=0;i<dim;i++)
         tmp += molefrac[i]*b(i);
      return tmp;
   }
   else return b(0);
}

void thermolib::send(double* pc, double* tc, double* w, double* y)
{
   for (i=0;i<dim;i++)
   {
      Pc[i] = pc[i]*101.325;
      Tc[i] = tc[i];
      omega[i] = w[i];
      molefrac[i] = y[i];
   }
}

double thermolib::P()
{
   task=0;
   pressure = 8.3144*temperature/molevolume;
   solver->set(this, pressure, 1.001*pressure);
   success=solver->run();
   return pressure/101.325;
}

double thermolib::T()
{
   task=1;
   temperature = pressure*molevolume/8.144;
   solver->set(this, temperature, 1.001*temperature);
   success=solver->run();
   return temperature;
}

double thermolib::v()
{
   if (mole>EPS)
   {
      task=2;
      molevolume = 8.3144*temperature/pressure;
      solver->set(this, molevolume, 1.001*molevolume);
      success=solver->run();
      return 0.001*mole*molevolume;
   }
   else return 0.0;
}

double thermolib::Zv()
{
   task=4;
   solver->set(this, 1.0, 0.99);
   success=solver->run();
   return Z;
}

double thermolib::phiV(int i)
{
   return exp((Z-1)*B(i)/B() - log(Z-B()) - A()/B()*(2*sqrt(A(i)/A())-B(i)/B())*log(1+B()/Z));
}

double thermolib::phiL(int i)
{
   Pr = pressure/Pc[i];
   Tr = temperature/Tc[i];
   tmp = 2.05135 - 2.10899/Tr - 0.19396*pow(Tr,2) + 0.02282*pow(Tr,3) + (0.08852 - 0.00872*pow(Tr,2))*Pr + (-0.00353 - 0.00203*Tr)*pow(Pr,2) - log10(Pr);
   tmp += omega[i]*(-4.23893 + 8.65808*Tr - 1.2206/Tr - 3.15224*pow(Tr,3) - 0.025*(Pr-0.6));
   return pow(10, tmp);
}

double thermolib::f(double x)
{
   if (task==0)
   {
      pressure=x;
      x= 8.3144*temperature/(molevolume-b_mix()) - a_mix()/(pow(molevolume,2)+b_mix()*molevolume) - x;
   }
   if (task==1)
   {
      temperature=x;
      x= 8.3144*x/(molevolume-b_mix()) - a_mix()/(pow(molevolume,2)+b_mix()*molevolume) - pressure;
   }
   if (task==2)
   {
      molevolume=x;
      x= 8.3144*temperature/(x-b_mix()) - a_mix()/(pow(x,2)+b_mix()*x) - pressure;
   }
   if(task==4)
   {
      Z=x;
      x= (pow(x,3)-pow(x,2)+(A()-B()-pow(B(),2))*x-A()*B());
   }
   return x;
}
