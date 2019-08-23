#include "cashflow.hpp"
using namespace std;

cashflow::cashflow(int n)
{
   N=n+1;
   Inv = new double[N];
   Coper = new double[N];
   Amort = new double[N];
   Rev = new double[N];
   Flow = new double[N];
   Flowact = new double[N];
   Itot=Ctot=Rtot=i_rate=a_rate=0.0;


  yield_tab[ 0] = 0.515;
  yield_tab[ 1] = 0.778;
  yield_tab[ 2] = 0.812;
  yield_tab[ 3] = 0.893;
  yield_tab[ 4] = 0.985;
  yield_tab[ 5] = 0.837;
  yield_tab[ 6] = 0.849;
  yield_tab[ 7] = 0.746;
  yield_tab[ 8] = 0.812;
  yield_tab[ 9] = 0.954;
  yield_tab[10] = 0.999;
  yield_tab[11] = 0.961;
  yield_tab[12] = 0.815;
  yield_tab[13] = 0.886;
  yield_tab[14] = 0.922;

}

cashflow::~cashflow ( void ) {
  delete [] Inv;
  delete [] Coper;
  delete [] Amort;
  delete [] Rev;
  delete [] Flow;
  delete [] Flowact;
}

bool cashflow::run()
{
   if(Itot<EPS || Ctot<EPS || Rtot<EPS || i_rate<EPS || a_rate<EPS)
      OK=false;
   else
   {
      //if(!MUTE)cout<<endl<<"       investments flow... OK";
      set_Inv();
      //if(!MUTE)cout<<endl<<"       depreciation flow... OK";
      set_Amort();
      set_C_R();
      //if(!MUTE)cout<<endl<<"       costs flow... OK";
      //if(!MUTE)cout<<endl<<"       revenus flow... OK";
      for(i=0;i<N;i++)
      {
         Flow[i] = (Rev[i]-Coper[i])*(1.0-a_rate)-(Inv[i]-a_rate*Amort[i]);
         Flowact[i] = Flow[i]/pow(1.0+i_rate, i);
      }
      //if(!MUTE)cout<<endl<<"       cash flow... OK";
      //if(!MUTE)cout<<endl<<"       actualizing cash flow... OK";
      OK=true;
      
//       cout<<endl<<endl<<"   CASH FLOW DETAILS"<<endl;
//       cout<<endl<<"      "<<setfill('-')<<setw(76)<<" ";
//       cout<<endl<<"      "<<" i "<<" Investment "<<"    Deprec. "
// 	  <<"   Expenses "<<"    Revenus "<<"  Cash flow "<<"  Act. flow ";
//       cout<<endl<<"      "<<setfill('-')<<setw(76)<<" ";
//       cout<<setfill(' ')<<setiosflags(ios::fixed)<<setprecision(0);
//       for(i=0;i<N;i++)
//          cout<<endl
// 	     <<"      "<<setw(2)<<i<<" "<<setw(11)
// 	     <<Inv[i]<<" "<<setw(11)<<Amort[i]<<" "
// 	     <<setw(11)<<Coper[i]<<" "<<setw(11)<<Rev[i]<<" "<<setw(11)<<Flow[i]<<" "<<setw(11)<<Flowact[i];
//       cout<<endl<<"      "<<setfill('-')<<setw(76)<<" ";
   }
   return OK;
}

void cashflow::set_Amort()
{
   Amort[0] = 0.0;
   temp=Itot;
   for(i=1;i<N-1;i++)
   {
      temp+=Inv[i];
      Amort[i] = temp/double(N-i);
      temp-=Amort[i];
   }
   Amort[N-1]=Amort[N-2];
}

void cashflow::set_Inv()
{
   Inv[0] = Itot;
   for(i=1;i<N-1;i++)
   {
      if((i)%5==0) Inv[i]=0.1*Itot;
      else Inv[i]=0.0;
   }
   Inv[N-1]=0.0;
   for(i=0;i<N-1;i++) Inv[N-1]-=0.1*Inv[i];
}

void cashflow::set_C_R()
{
   Coper[0] = Rev[0] = 0.0;
   for(i=1;i<N;i++)
   {
      Coper[i] = yield(i)*Ctot;
      Rev[i] = yield(i)*Rtot;
   }
}
