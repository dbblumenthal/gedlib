
Set t time periods
/
2005, 2015, 2025, 2035, 2045, 2055, 2065, 2075, 2085, 2095,
2105, 2115, 2125, 2135, 2145, 2155, 2165, 2175, 2185, 2195,
2205, 2215, 2225, 2235, 2245, 2255, 2265, 2275, 2285, 2295,
2305, 2315, 2325, 2335, 2345, 2355, 2365, 2375, 2385, 2395,
2405
/

subtGDP(t) indice sommation GDP
/
2005, 2015, 2025, 2035, 2045, 2055, 2065, 2075, 2085, 2095
/

subtem(t) indice sommation emissions
/
2005, 2015, 2025, 2035, 2045, 2055, 2065, 2075, 2085
/
;

Set ts(t) time periods
/
2005, 2015, 2025, 2035, 2045, 2055, 2065, 2075, 2085, 2095,
2105, 2115, 2125, 2135, 2145, 2155, 2165, 2175, 2185, 2195,
2205, 2215, 2225, 2235, 2245, 2255, 2265, 2275, 2285, 2295,
2305, 2315, 2325, 2335, 2345, 2355, 2365, 2375, 2385, 2395,
2405
/;

Set tfirst(t) first time period;
tfirst(t) = YES$(ord(t) eq 1);

Set tlast(t) last time period;
tlast(t) = YES$(ord(t) eq card(t));

Parameter nbp(t) number of periods since 2005
/
  2005 0, 2015 1, 2025 2, 2035 3, 2045 4, 2055 5, 2065 6, 2075 7, 2085 8,
  2095 9, 2105 10, 2115 11, 2125 12, 2135 13, 2145 14, 2155 15, 2165 16,
  2175 17, 2185 18, 2195 19, 2205 20, 2215 21, 2225 22, 2235 23, 2245 24,
  2255 25, 2265 26, 2275 27, 2285 28, 2295 29, 2305 30, 2315 31, 2325 32,
  2335 33, 2345 34, 2355 35, 2365 36, 2375 37, 2385 38, 2395 39, 2405 40
/;

Parameter nby(t) number of years since 2005
/
 2005 0, 2015 10, 2025 20, 2035 30, 2045 40, 2055 50, 2065 60, 2075 70, 2085 80,
 2095 90, 2105 100, 2115 110, 2125 120, 2135 130, 2145 140, 2155 150, 2165 160,
 2175 170, 2185 180, 2195 190, 2205 200, 2215 210, 2225 220, 2235 230, 2245 240,
 2255 250, 2265 260, 2275 270, 2285 280, 2295 290, 2305 300, 2315 310, 2325 320,
 2335 330, 2345 340, 2355 350, 2365 360, 2375 370, 2385 380, 2395 390, 2405 400
/;

Parameter xi(t) indicator for the clean economy
/
  2005 1, 2015 1, 2025 1, 2035 1, 2045 1, 2055 1, 2065 1, 2075 1, 2085 1,
  2095 1, 2105 1, 2115 1, 2125 1, 2135 1, 2145 1, 2155 1, 2165 1, 2175 1,
  2185 1, 2195 1, 2205 1, 2215 1, 2225 1, 2235 1, 2245 1, 2255 1, 2265 1,
  2275 1, 2285 1, 2295 1, 2305 1, 2315 1, 2325 1, 2335 1, 2345 1, 2355 1,
  2365 1, 2375 1, 2385 1, 2395 1, 2405 1
/;
$ontext
/
  2005 0, 2015 0, 2025 0, 2035 0, 2045 0, 2055 0, 2065 0, 2075 0, 2085 0,
  2095 0, 2105 0, 2115 0, 2125 0, 2135 0, 2145 0, 2155 0, 2165 0, 2175 0,
  2185 0, 2195 0, 2205 0, 2215 0, 2225 0, 2235 0, 2245 0, 2255 0, 2265 0,
  2275 0, 2285 0, 2295 0, 2305 0, 2315 0, 2325 0, 2335 0, 2345 0, 2355 0,
  2365 0, 2375 0, 2385 0, 2395 0, 2405 0
/;
$offtext

Parameter GDP(t) GDP donne par le scenario B2MESSAGE par moyenne simple
/
  2005 49.84, 2015 66.53, 2025 86.94, 2035 112.9, 2045 145.3, 2055 182.0, 2065 220.7, 2075 259.1, 2085 295.5,
  2095 331.7, 2105 0, 2115 0, 2125 0, 2135 0, 2145 0, 2155 0, 2165 0, 2175 0,
  2185 0, 2195 0, 2205 0, 2215 0, 2225 0, 2235 0, 2245 0, 2255 0, 2265 0,
  2275 0, 2285 0, 2295 0, 2305 0, 2315 0, 2325 0, 2335 0, 2345 0, 2355 0,
  2365 0, 2375 0, 2385 0, 2395 0, 2405 0
/;

Parameter emissionB2(t) emission CO2 par decennie donne par le scenario B2MESSAGE par moyenne simple
/
  2005 8.65, 2015 9.19, 2025 9.88, 2035 10.57, 2045 11.05, 2055 11.40, 2065 11.72, 2075 12.18, 2085 12.76,
  2095 0, 2105 0, 2115 0, 2125 0, 2135 0, 2145 0, 2155 0, 2165 0, 2175 0,
  2185 0, 2195 0, 2205 0, 2215 0, 2225 0, 2235 0, 2245 0, 2255 0, 2265 0,
  2275 0, 2285 0, 2295 0, 2305 0, 2315 0, 2325 0, 2335 0, 2345 0, 2355 0,
  2365 0, 2375 0, 2385 0, 2395 0, 2405 0
/;

Scalars
    borneinf borne inferieure de DICE                      /0.9/
    bornesup borne superieure de DICE                      /1.1/
    deltaalpha ecart entre modeles DICE et GERAD  /0/
    deltaGDP ecart DICE GERAD pour le GDP /0/
    deltaem idem                         /0/
    moyenneGDP moyenne GDP B2-message                /175.05/
    moyenneem idem                             /10.82/

    vb       slope w.r.t K2(t) of the probability rate of discovery /0.0019/
    wb       initial probability rate of discovery                  /0.05/
    beta     marginal atmospheric retention rate                    /0.64/
    catM     catastrophic carbon concentration                      /2059.0/
    delta_M  carbon removal rate per decade                         /0.036/
*    mtargM   target for carbon concentration                        /1025.00/
    mtargM   target for carbon concentration                        /5000.00/
    piE1     initial energy price in dirty economy                          /0.35/
    piE2     initial energy price in clean economy                          /0.60/
    rhop  discount rate /0.025/

    prstp    initial rate of social time preference per year        /0.015/
    dr       decline rate of social time preference per year        /0.000001/
    rho      discount rate                                          /0.03/

    L0       2005 world population (millions)                       /6409/
    gL0      initial value for growth rate of population            /0.08/
    dgL      rate of decrease for growth rate of population         /0.3/

    scale1 scaling coefficient in objective function                /300/
    scale2 scaling coefficient in objective function                /1.00E+07/

*Scalaires calibres
$include "input.txt";

;

*Initialisation des parametres
Parameter L(t) labor - world population (millions);
L(t) =  L0 * exp( (gl0/dgL)*(1-exp( -dgL*nbp(t) )) );

Parameter r(t) instantaeous rate of social time preference;
r(t)=prstp*EXP(-dr*10*(ord(t)-1));

Parameter rr(t) average utility social discount rate;
rr(tfirst)=1;
loop(t, rr(t+1) = rr(t) / ((1+r(t))**10) );

Parameter A(t) total factor productivity;
A(tfirst) = 1.38734*A0;
*loop(t, A(t+1) = A(t) / (1-gA0*exp(-dgA*10*nbp(t))) );
A(t) = 1.38734*A0 * exp( (gA0/dgA)*(1-exp( -dgA*nbp(t) )) );

Parameter phid(t) energy efficiency in the dirty economy;
phid(t) = (phid0 * exp( (-gphid0/dgphid)*exp(-dgphid*(nbp(t)+1)) ) ) /
          exp(-gphid0/dgphid);

Parameter phic(t) energy efficiency in the clean economy;
phic(t) = (phic0 * exp( (-gphic0/dgphic)*exp(-dgphic*(nbp(t)+1)) ) ) /
          exp(-gphic0/dgphic);

Parameter thetad(t) dirty energy elasticity in production function;
thetad(t) = (thetad0 * exp( (gthetad0/dgthetad)*exp(-dgthetad*(nbp(t)+1)) ) ) /
            exp(gthetad0/dgthetad);

Parameter thetac(t) clean energy elasticity in production function;
thetac(t) = (thetac0 * exp( (gthetac0/dgthetac)*exp(-dgthetac*(nbp(t)+1)) ) ) /
            exp(gthetac0/dgthetac);

Parameter EM(t) total emissions (over 10 years);
EM(t) = 0;

Parameter TK(t) total capital;
TK(t) = 0;

Parameter deltaalphaGDP(t) ecart courant entre GDP de DICE et GERAD;
deltaalphaGDP(t) = 0;

*Definition du modele
Positive variables
    C(t)  total consumption (trillion dollars)
    ELF(t) economic loss factor (%)
    E1(t) carbon emissions from the dirty economy (GtC)
    E2(t) carbon emissions from the clean economy (GtC)
    I1(t) investment in dirty capital K1 (trillion dollars)
    I2(t) investment in clean capital K2 (trillion dollars)
    K1(t) physical stock of dirty capital (trillion dollars)
    K2(t) physical stock of clean capital (trillion dollars)
    L1(t) labor for the dirty economy (millions)
    L2(t) labor for the clean economy (millions)
    M(t)  atmospheric carbon concentration
    Y(t)  economic output (trillion dollars);

Variables
    PW(t)   welfare of current period
    W2      total discounted welfare after the last jump (V_2-1-l);

Equations
    ALOCLABOR(t) compute allocation of total labor
    CARBDYNAM(t) compute evolution of carbon concentrations
    E1FIRST(t)    initial value for carbon emissions from dirty capital
    E2FIRST(t)    initial value for carbon emissions from clean capital
    ELFCOMP(t)   compute economic loss factor
    I2MAX(t)      bound for investment on clean capital
    K1ACCUMUL(t) compute dirty capital accumulation
    K1FIRST(t)   initial value for dirty capital
    K1LAST(t)    final value for dirty capital
    K2ACCUMUL(t) compute clean capital accumulation
    K2FIRST(t)   initial value for clean capital
    K2LAST(t)    final value for clean capital
    MFIRST(t)    initial value for carbon concentration
    TARGETM(t)   impose cap on atmospheric carbon concentrations
    TOTALCONS(t) compute total consumption
    TOTALPROD(t) compute total production
    PWEQ(t)      compute welfare of current period
    WELFARE2     compute total welfare after the last jump;

ALOCLABOR(t)..
  L1(t) + L2(t) =l= L(t);

CARBDYNAM(t+1)..
  M(t+1) =e= beta*10*(E1(t)+E2(t)) + M(t)*(1-delta_M) + delta_M*590;

MFIRST(tfirst)..
* 2005 value around 385 ppmv = 808.8 GtC
  M(tfirst) =e= 808.8;

ELFCOMP(t)..
  ELF(t) =e= 1;

E1FIRST(tfirst)..
* DICE-2007 value
  E1(tfirst) =e= 6.7;

E2FIRST(tfirst)..
* 0.0001 of DICE-2007 value
  E2(tfirst) =e= 0.00067;

I2MAX(t)..
  I2(t) =l= 0.97*Y(t);

K1FIRST(tfirst)..
  K1(tfirst) =e= 137.0;

K1ACCUMUL(t)..
  K1(t+1) =l= K1(t)*(1-dk)**10 + 10*I1(t);

K1LAST(tlast)..
  0.02*K1(tlast) =l= I1(tlast);

K2FIRST(tfirst)..
  K2(tfirst) =e= 0.0137;

K2ACCUMUL(t)..
  K2(t+1) =l= K2(t)*(1-dk)**10 + 10*I2(t);

K2LAST(tlast)..
  0.02*K2(tlast) =l= I2(tlast);

TARGETM(t)..
* Constraint is always active
  M(t) =l= mtargM;

TOTALPROD(t)..
  Y(t) =e= A(t) * (
  K1(t)**alphad * (phid(t)*E1(t))**thetad(t) * L1(t)**(1-alphad-thetad(t)) +
  xi(t)*K2(t)**alphac * (phic(t)*E2(t))**thetac(t) * L2(t)**(1-alphac-thetac(t))
                  );

TOTALCONS(t)..
  C(t) =e= Y(t) - I1(t) - I2(t) - piE1*phid(t)*E1(t) - piE2*phic(t)*E2(t);

PWEQ(t)..
  PW(t) =e= ( (ELF(t)*C(t)/L(t))**(1-elasmu) -1) / (1-elasmu);

WELFARE2..
  W2 =e= sum(t, 10*exp(-rho*nby(t))*L(t)*(PW(t))/scale1) + scale2 ;

Model ECM /all/ ;
ECM.OPTFILE =1;

C.lo(t)    = 1;
C.up(t)    = 1000;
E1.lo(t)   = 0.0000001;
E1.up(t)   = 500;
E2.lo(t)   = 0.0000001;
E2.up(t)   = 100;
ELF.lo(t)  = 0.5;
ELF.up(t)  = 1.0;
I1.lo(t)   = 0.000000001;
I1.up(t)   = 1000;
I2.lo(t)   = 0.000000001;
I2.up(t)   = 1000;
K1.lo(t)   = 0.01;
K1.up(t)   = 3000;
K2.lo(t)   = 0.000000001;
K2.up(t)   = 3000;
L1.lo(t)   = 0.1;
L1.up(t)   = 10000;
L2.lo(t)   = 0.000000001;
L2.up(t)   = 10000;
M.lo(t)    = 300;
M.up(t)    = 4000;
Y.lo(t)    = 50;
Y.up(t)    = 1000;

option iterlim = 99900;
option reslim = 99999;
option solprint = off;
option limrow = 0;
option limcol = 0;

Solve ECM using nlp maximising W2;

EM(t) = E1.L(t)+E2.L(t);
TK(t) = K1.L(t)+K2.L(t);

deltaalphaGDP(t) = Y.L(t)-GDP(t);
deltaGDP = sum ( subtGDP(t) , sqr(deltaalphaGDP(t)) );
deltaem = sum ( subtem(t) , sqr(EM(t)-emissionB2(t)) );
deltaalpha = deltaGDP / moyenneGDP  + deltaem / moyenneem ;

file FOUT / output.txt /;
put FOUT;
put deltaalpha:13:6;
putclose;






