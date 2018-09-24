#include "split.hpp"
using namespace std;

split::split(int n, stream* s1, stream** s2)
{
  nb_out=n;
  in = s1;
  out= s2;
  success = true;
}

bool split::solve()
{
   tmp=0;
   for (i=0;i<nb_out; i++) tmp+=frac[i];
   if(fabs(1-tmp)<=EPS)
   {
      success = true;
      for (i=0; i<nb_out;i++)
      {
         out[i]->m=0;
         for (j=0;j<in->nb;j++)
         {

	   out[i]->chem[j]->m = frac[i]*in->chem[j]->m;
	   out[i]->m += out[i]->chem[j]->m;
         }
         out[i]->set(in->P, in->T);
         // out[i]->write(); // TOTO
      }
   }
   tmp=0; for(i=0;i<nb_out;i++) tmp+=out[i]->m;
   if(fabs(tmp-in->m)>EPS)
   {
//       logf.open(MESSAGES, ios::app);
//       logf<<"   --> Warning <--  Block "<<name<<" is not in mass balance ("<<fabs(tmp-in->m)/tmp<<").\n";
//       logf.close();
      success = false;
   }
   else success = true;
   return success;
}

void split::write()
{
  cout << "WRITE FILE " << RUNTIME << name << ".unit" << " :\n\tBEGIN\n";
  cout <<"\t>>         "<<name;
  cout << endl<<"\t>>           stream in: "<<in->name;;
  cout<<endl<<"\t>>           streams out: "<<setprecision(3);
  for ( i = 0 ; i < nb_out ; i++ )
    cout << out[i]->name<<" ("<<frac[i]<<")  ";
  cout << "\n\tEND\n\n";
}
