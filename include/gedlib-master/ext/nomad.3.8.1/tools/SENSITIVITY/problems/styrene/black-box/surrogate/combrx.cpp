#include "combrx.hpp"
using namespace std;

combrx::combrx ( const string & cas ) {

  // combustion.prop :
  // -----------------
  CAS = cas;

  if (CAS=="64-17-5") {
    nO2  = 3;
    nCO2 = 2;
    nH2O = 3;
  }
  else if (CAS=="74-82-8") {
    nO2  = 2;
    nCO2 = 1;
    nH2O = 2;
  }
  else if (CAS=="1333-74-0") {
    nO2  = 0.5;
    nCO2 = 0;
    nH2O = 1;
  }
  else if (CAS=="100-42-5") {
    nO2  = 10;
    nCO2 = 8;
    nH2O = 4;
  }
  else if (CAS=="74-85-1") {
    nO2  = 3;
    nCO2 = 2;
    nH2O = 2;
  }
  else if (CAS=="108-88-3") {
    nO2  = 9;
    nCO2 = 7;
    nH2O = 4;
  }
  else if (CAS=="100-41-4") {
    nO2  = 10.5;
    nCO2 = 8;
    nH2O = 5;
  }
  else if (CAS=="71-43-2") {
    nO2  = 7.5;
    nCO2 = 6;
    nH2O = 3;
  }
  else {
    cout << "ERROR 21" << endl;
    exit(0);
  }

  COMB = new chemical(CAS);
  O2   = new chemical("7782-44-7");
  N2   = new chemical("7727-37-9");
  CO2  = new chemical("124-38-9");
  H2O  = new chemical("7732-18-5");
  Hro  = CO2->Ho*nCO2 + nH2O*(H2O->Ho - H2O->dHvap) - COMB->Ho;
  LFLo = -3420.0/Hro + 0.569e-3*Hro + 0.0538e-6*pow(Hro,2) + 1.8;
  LFLo = LFLo/100.0;
  UFLo = 0.0063*Hro + 0.567e-6*pow(Hro, 2) + 23.5;
  UFLo = UFLo/100.0;

}
