#include "profitability.hpp"
#include "secant.cpp"
using namespace std;

bool profitability::run ( double * y )
{
   OK=true;

   // cout<<setiosflags(ios::fixed);
   // cout<<endl<<endl<<"   PROFITABILITY ANALYSIS"<<endl;
   // cout<<endl<<setprecision(1)<<"      Return on investment   (%)= "<<ROI()*100.0;
   // cout<<endl<<"      Rate of return         (%)= "<<RR()*100.0; // y[10]

   y[10] = RR()*100.0;  // y[10]

   ROI();
   RR();

   // cout<<endl<<"      Discounted flow rate   (%)= "<<DFR()*100.0;

   DFR();

   // cout<<endl<<endl<<"      Payout time            (y)= "<<PT();              // y[ 9]
   y[9] = PT();

   // cout<<endl<<setprecision(0)<<"      Annual equivalent cost ($)= "<<AEC();  // y[12]
   y[12] = AEC();

   // cout<<endl<<"      Net present value      ($)= "<<NPV()<<endl;             // y[13]
   y[13] = NPV();

   return OK;
}

double profitability::ROI()
{
  // if(!MUTE)cout<<endl<<"        return on investment...";
   num=den=0.0;
   for(i=0;i<C->N;i++)
   {
      if(C->Inv[i]>EPS) den+=C->Inv[i];
      num+=(C->Rev[i]-C->Coper[i]-C->Amort[i]);
   }
   if (num>EPS && den>EPS && C->N>0) {
     // if(!MUTE)cout<<" OK";
     return num/C->N/den;
   }
   else return 0.0;
}

double profitability::RR()
{
   // if(!MUTE)cout<<endl<<"        rate of return...";
   num=den=0.0;
   for(i=0;i<C->N;i++)
   {
      num+=(C->Rev[i]-C->Coper[i])/pow(1.0+C->i_rate, i);
      den+=C->Inv[i]/pow(1.0+C->i_rate, i);
   }
   if(num>EPS && den>EPS) {
     // if(!MUTE)cout<<" OK";
     return num/den;
   }
   else return 0.0;
}

double profitability::DFR()
{
  //if(!MUTE)cout<<endl<<"        discounted cash flow rate...";
   solver = new secant<profitability>();
   solver->set(this, 0.0, 0.01);
   OK = solver->run();

   if ( OK && num>EPS && num < 1e20 ) {
     // if(!MUTE)cout<<" OK";
     return num;
   }
   else return 0.0;
}

double profitability::f(double x)
{
   num=x;
   sum=0.0;
   for(i=0;i<C->N;i++)
     sum += C->Flow[i]/pow(1.0+x, i);
   return sum;
}

double profitability::PT()
{
  // if(!MUTE)cout<<endl<<"        payout time...";
   sum=0.0;
   for(i=0;i<C->N;i++)
   {
      if((sum+C->Flow[i])>0.0)
      {
         den=0.0;
         while(sum+den*C->Flow[i]<=0.0) den+=0.001;
         den+=double(i-1);
         i=C->N;
      }
      else sum+=C->Flow[i];
   }

   if(den>EPS) {
     // if(!MUTE)cout<<" OK";
     return den;
   }
   else return 0.0;
}

double profitability::AEC()
{
   //if(!MUTE)cout<<endl<<"        annual equivalent cost...";
   sum=0.0;
   for(i=0;i<C->N;i++) sum+=(C->Coper[i]+C->Inv[i])/pow(1.0+C->i_rate, i);
   if (sum>EPS) {
//      if(!MUTE)
//        cout<<" OK";
     return sum*(C->i_rate*pow(1.0+C->i_rate,C->N))/(pow(1.0+C->i_rate,C->N)-1.0);
   }
   else return 0.0;
}

double profitability::NPV()
{
  // if(!MUTE)cout<<endl<<"        net present value...";
   sum=0.0;
   for ( i = 0 ; i < C->N ; i++ )
     sum += C->Flowact[i];
   if ( sum > EPS ) {
//      if(!MUTE)
//        cout<<" OK";
     return sum;
   }
   return 0.0;
}
