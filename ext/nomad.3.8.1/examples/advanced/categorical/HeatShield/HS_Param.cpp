#include "HS_Param.hpp"
using namespace std;

/*---------------*/
/*  constructor  */
/*---------------*/
HS_Param::HS_Param ( int    load_flag  ,
		     double load       ,
		     double max_weight ,
		     double TE_limit   ,
		     double T_cold     ,
		     double T_hot      ,
		     double L            ) : _load_flag  ( load_flag  ) ,
					     _load       ( load       ) ,
					     _max_weight ( max_weight ) ,
					     _TE_limit   ( TE_limit   ) ,
					     _T_cold     ( T_cold     ) ,
					     _T_hot      ( T_hot      ) ,
					     _L          ( L          )   {
  // materials:
  // ----------
  _material = new HS_Material * [_CARBON_STEEL_ + 1];

  double T1[18] = { 0 , 7.2 , 40 , 80 , 100 , 140 , 180 , 200 , 240 ,
		    280 , 300 , 340 , 380 , 400 , 440 , 480 , 500 , 540 };
  double T2[ 5] = { 0 , 20 , 77 , 195 , 297 };

  double k1[18] = { 0 , 0.0072 , 0.064 , 0.125 , 0.148 , 0.17 , 0.185 ,
		    0.188 , 0.194 , 0.2 , 0.202 , 0.202 , 0.202 , 0.202,
		    0.202 , 0.202 , 0.202 , 0.202 };
  double k2[18] = { 0 , 0.0266 , 0.0855 , 0.117 , 0.125 , 0.135 , 0.142 ,
		    0.142 , 0.143 , 0.146 , 0.148 ,
		    0.15 , 0.15 , 0.15 , 0.15 , 0.15 , 0.15 , 0.15 };
  double k3[ 5] = { 0 , 0.0015 , 0.002 , 0.003 , 0.0035 };
  double k4[ 5] = { 0 , 0.002 , 0.0026 , 0.0044 , 0.005 };
  double k5[18] = { 0 , 20 , 105 , 158 , 159 , 136 , 121 , 119 , 117 ,
		    116 , 116 , 116 , 116 , 116 , 116  ,116 , 116 , 116 };
  double k6[18] = { 0 , 0.14 , 1.26 , 3 , 3.7 , 4.68 , 5.49 , 5.72 , 5.24 ,
		    6.53 , 7 , 7.22 , 7.64 , 7.75 , 8.02 , 8.26 , 8.43 , 8.67 };
  double k7[18] = { 0 , 1.7 , 15 , 25.1 , 29.5 , 34.4 , 36.4 , 37 ,
		    37.6 , 37.6 , 37.6 , 37.6 , 37.6 ,
		    37.6 , 37.6 , 37.6 , 37.6 , 37.6 };

  double s1[18] = { 29.9 , 29.8 , 29.6 , 29.1 , 28.9 , 28 , 27 , 26.5 , 25 ,
		    23.9 , 23.2 , 21 , 19.3 , 18.5 , 15.8 , 13.7 , 12 , 9  };
  double s2[18] = { 18.8 , 18.6 , 17.8 , 16.8 , 16.2 , 14.8 , 13.9 , 13 , 11.5 ,
		    10.6 , 9.6 , 8.2 , 7.6 , 6.4 , 4.8 , 3.9 , 3.2 , 1.9 };
  double s3[4]  = { 117 , 125 , 90 , 81 };
  double s5[18] = { 51 , 50 , 43 , 36 , 35 , 32 , 29.5 , 29 ,  28 ,  28 , 28 ,
		    28 , 28 , 28 , 28 , 28 , 28 , 28 };
  double s6[18] = { 195.5 , 195 , 194 , 193 , 192 , 190 , 189 , 188 , 187 ,
		    186 , 186 , 182 , 181 , 180 , 178 , 176 , 175 , 173 };
  double s7[18] = { 102 , 100 , 92 , 87 , 82 , 78 , 72 , 70 , 66 , 63 , 60 ,
		    58 , 55 , 53 , 51 , 49 , 48 , 47 };
  double e1[17] = { 0.0145 , 0.01437, 0.01404 , 0.01376 , 0.01314 ,
		    0.01233 , 0.01186 , 0.01078 , 0.00979 , 0.00918 ,
		    0.00798 , 0.00664 , 0.00595 , 0.00444 ,
		    0.00273 , 0.00184 , 0 };
  double e2[17] = { 0.02695 , 0.02643 , 0.02563 , 0.02514 , 0.02406 ,
		    0.02285 , 0.02221 , 0.02087 , 0.0194 , 0.01858 ,
		    0.01671 , 0.01432 , 0.013 , 0.01024 , 0.00735 ,
		    0.00591 , 0 };
  double e3[5] = { 0.00249 , 0.002463 , 0.00214 , 0.001123 , 0 };
  double e4[5] = { 0.00249 , 0.002463 , 0.00214 , 0.001123 , 0 };

  double e5[17] = { 0.00431 , 0.00431 , 0.00428 , 0.00423 , 0.004090 ,
		    0.00385 , 0.00371 , 0.00338 , 0.00301 , 0.00281 ,
		    0.00239 , 0.00195 , 0.00172 , 0.00125 , 0.000750 ,
		    0.0005 , 0 };
  double e6[17] = { 0.00304 , 0.003051 , 0.003051 , 0.003021 , 0.00291 ,
		    0.00274 , 0.00263 , 0.00239 , 0.00212 , 0.00198 ,
		    0.00167 , 0.00137 , 0.00121 , 0.00088 , 0.00053 ,
		    0.00035 , 0 };
  double e7[17] = { 0.0021 , 0.0021 , 0.00209 , 0.00207 , 0.00201 ,
		    0.0019 , 0.00183 , 0.00168 , 0.00151 , 0.00141 ,
		    0.0012 , 0.00098 , 0.00087 , 0.00064 , 0.00036 ,
		    0.00025 , 0 };
  double y1[18] = { 1.41 , 1.40 , 1.35 , 1.26 , 1.2 , 1.1 , 1.02 ,
		    .96 , .83 , .76 , .71 , .58 , .52 , .50 , .46 ,
		    .435 , .425 , .42 };  
  double y2[18] = { .84 , .83 , .81 , .78 , .76 , .73 , .71 , .68 ,
		    .60 , .55 , .48 , .32 , .235 , .19 , .12 ,
		    .09 , .07 , .05 };
  double y3[4] = { 4.27 , 3.98 , 3.50 , 3.29 };
  double y5[18] = { 11 , 11 , 11 , 11 , 11 , 11 , 11 , 11 , 10.8 ,
		    10.7 , 10.6 , 10.5 , 10.4 , 10.3 , 10.2 , 10.1 ,
		    10 , 10 };
  double y6[18] = { 29 , 28.9 , 28.6 , 28.2 , 28 , 27.6 , 27.2 , 27 ,
		    26.6 , 26.2 , 26 , 25.6 , 25.2 , 25 , 24.6 , 24.2 ,
		    24 , 23.6 };
  double y7[18] = { 31.58 , 31.56 , 31.46 , 31.31 , 31.24 , 31.1 , 30.96 ,
		    30.89 , 30.76 , 30.62 , 30.55 , 30.41 , 30.28 , 30.21 ,
		    30.07 , 29.93 , 29.86 , 29.72 };

  double tmp1 , tmp2;

  // nylon:
  {
    _material[_NYLON_] = new HS_Material ( _NYLON_ , 0.0379 );
    _material[_NYLON_]->set_T ( 18 , T1 );
    _material[_NYLON_]->set_k ( 18 , k1 );
    _material[_NYLON_]->set_s ( 18 , s1 );
    _material[_NYLON_]->set_kABC();
    _material[_NYLON_]->set_sABC();
    _material[_NYLON_]->set_TCI();
    _material[_NYLON_]->set_e ( 17 , e1 );
    _material[_NYLON_]->set_ekABC1();
    tmp1 = 4;
    _material[_NYLON_]->thermal_expansion ( 1 , &tmp1 , &tmp2 );
    double etmp[18] = { 0.0145 , -1 , 0.01437, 0.01404 , 0.01376 , 0.01314 ,
			0.01233 , 0.01186 , 0.01078 , 0.00979 , 0.00918 ,
			0.00798 , 0.00664 , 0.00595 , 0.00444 ,
			0.00273 , 0.00184 , 0 };
    etmp[1] = tmp2;
    _material[_NYLON_]->set_e ( 18 , etmp );
    _material[_NYLON_]->set_ekABC2();

    _material[_NYLON_]->set_ksABC();
    _material[_NYLON_]->set_TEI();

    _material[_NYLON_]->set_y ( 18 , y1 );
    _material[_NYLON_]->set_yABC();
  }

  // teflon:
  {
    _material[_TEFLON_] = new HS_Material ( _TEFLON_ , 0.054 );
    _material[_TEFLON_]->set_T ( 18 , T1 );
    _material[_TEFLON_]->set_k ( 18 , k2 );
    _material[_TEFLON_]->set_s ( 18 , s2 );
    _material[_TEFLON_]->set_kABC();
    _material[_TEFLON_]->set_sABC();
    _material[_TEFLON_]->set_TCI();
    _material[_TEFLON_]->set_e ( 17 , e2 );
    _material[_TEFLON_]->set_ekABC1();

    tmp1 = 4;
    _material[_TEFLON_]->thermal_expansion ( 1 , &tmp1 , &tmp2 );
    double etmp[18] = { 0.02695 , -1 , 0.02643 , 0.02563 , 0.02514 , 0.02406 ,
			0.02285 , 0.02221 , 0.02087 , 0.0194 , 0.01858 ,
			0.01671 , 0.01432 , 0.013 , 0.01024 , 0.00735 ,
			0.00591 , 0 };
    etmp[1] = tmp2;
    _material[_TEFLON_]->set_e ( 18 , etmp );
    _material[_TEFLON_]->set_ekABC2();

    _material[_TEFLON_]->set_ksABC();
    _material[_TEFLON_]->set_TEI();

    _material[_TEFLON_]->set_y ( 18 , y2 );
    _material[_TEFLON_]->set_yABC();
  }

  // epoxy normal:
  {
    _material[_EPOXY_NORMAL_] = new HS_Material ( _EPOXY_NORMAL_ , 0.064 );
    _material[_EPOXY_NORMAL_]->set_T ( 5 , T2 , false );
    _material[_EPOXY_NORMAL_]->set_k ( 5 , k3 , false );
    _material[_EPOXY_NORMAL_]->set_s ( 4 , s3 );
    _material[_EPOXY_NORMAL_]->set_kABC();
    _material[_EPOXY_NORMAL_]->set_sABC();

    _material[_EPOXY_NORMAL_]->set_e ( 5 , e3 );
    _material[_EPOXY_NORMAL_]->set_ekABC2();
  
    double s32[5];
    s32[0]  = _material[_EPOXY_NORMAL_]->thermal_stress(0);
    s32[0] /= HS_Material::S_CC; 
    for ( int i = 1 ; i < 5 ; ++i )
      s32[i] = s3[i-1];
    _material[_EPOXY_NORMAL_]->set_s ( 5 , s32 );
    _material[_EPOXY_NORMAL_]->set_sABC();

    _material[_EPOXY_NORMAL_]->set_TCI();

    _material[_EPOXY_NORMAL_]->set_ksABC();
    _material[_EPOXY_NORMAL_]->set_TEI();

    _material[_EPOXY_NORMAL_]->set_y ( 4 , y3 );
    _material[_EPOXY_NORMAL_]->set_yABC();

    tmp1 = 0;
    _material[_EPOXY_NORMAL_]->thermal_YM ( 1 , &tmp1 , &tmp2 );
    tmp2 /= HS_Material::YM_CC;
    double ytmp[5] = { -1 , 4.27 , 3.98 , 3.50 , 3.29 };
    ytmp[0] = tmp2;
    _material[_EPOXY_NORMAL_]->set_y ( 5 , ytmp );
    _material[_EPOXY_NORMAL_]->set_yABC();
  }

  // epoxy plane:
  {
    _material[_EPOXY_PLANE_] = new HS_Material ( _EPOXY_PLANE_ , 0.064 );
    _material[_EPOXY_PLANE_]->set_T ( 5 , T2 , false );
    _material[_EPOXY_PLANE_]->set_k ( 5 , k4 , false );
    _material[_EPOXY_PLANE_]->set_s ( 4 , s3 );
    _material[_EPOXY_PLANE_]->set_kABC();
    _material[_EPOXY_PLANE_]->set_sABC();

    _material[_EPOXY_PLANE_]->set_e ( 5 , e4 );
    _material[_EPOXY_PLANE_]->set_ekABC2();

    double s32[5];
    s32[0] = _material[_EPOXY_PLANE_]->thermal_stress(0);
    s32[0] /= HS_Material::S_CC; 
    for ( int i = 1 ; i < 5 ; ++i )
      s32[i] = s3[i-1];
    _material[_EPOXY_PLANE_]->set_s ( 5 , s32 );
    _material[_EPOXY_PLANE_]->set_sABC();

    _material[_EPOXY_PLANE_]->set_TCI();

    _material[_EPOXY_PLANE_]->set_ksABC();
    _material[_EPOXY_PLANE_]->set_TEI();

    _material[_EPOXY_PLANE_]->set_y ( 4 , y3 );
    _material[_EPOXY_PLANE_]->set_yABC();

    tmp1 = 0;
    _material[_EPOXY_PLANE_]->thermal_YM ( 1 , &tmp1 , &tmp2 );
    tmp2 /= HS_Material::YM_CC;
    double ytmp[5] = { -1 , 4.27 , 3.98 , 3.50 , 3.29 };
    ytmp[0] = tmp2;
    _material[_EPOXY_PLANE_]->set_y ( 5 , ytmp );
    _material[_EPOXY_PLANE_]->set_yABC();
  }

  // aluminium:
  {
    _material[_ALUMINIUM_] = new HS_Material ( _ALUMINIUM_ , 0.098 );
    _material[_ALUMINIUM_]->set_T ( 18 , T1 );
    _material[_ALUMINIUM_]->set_k ( 18 , k5 );
    _material[_ALUMINIUM_]->set_s ( 18 , s5 );
    _material[_ALUMINIUM_]->set_kABC();
    _material[_ALUMINIUM_]->set_sABC();
    _material[_ALUMINIUM_]->set_TCI();

    _material[_ALUMINIUM_]->set_e ( 17 , e5 );
    _material[_ALUMINIUM_]->set_ekABC1();
    tmp1 = 4;
    _material[_ALUMINIUM_]->thermal_expansion ( 1 , &tmp1 , &tmp2 );
    double etmp[18] = { 0.00431 , -1 , 0.00431 , 0.00428 , 0.00423 , 0.004090 ,
			0.00385 , 0.00371 , 0.00338 , 0.00301 , 0.00281 ,
			0.00239 , 0.00195 , 0.00172 , 0.00125 , 0.000750 ,
			0.0005 , 0 };
    etmp[1] = tmp2;
    _material[_ALUMINIUM_]->set_e ( 18 , etmp );
    _material[_ALUMINIUM_]->set_ekABC2();

    _material[_ALUMINIUM_]->set_ksABC();
    _material[_ALUMINIUM_]->set_TEI();

    _material[_ALUMINIUM_]->set_y ( 18 , y5 );
    _material[_ALUMINIUM_]->set_yABC();
  }

  // steel:
  {
    _material[_STEEL_] = new HS_Material ( _STEEL_ , 0.28 );
    _material[_STEEL_]->set_T ( 18 , T1 );
    _material[_STEEL_]->set_k ( 18 , k6 );
    _material[_STEEL_]->set_s ( 18 , s6 );
    _material[_STEEL_]->set_kABC();
    _material[_STEEL_]->set_sABC();
    _material[_STEEL_]->set_TCI();

    _material[_STEEL_]->set_e ( 17 , e6 );
    _material[_STEEL_]->set_ekABC1();
    tmp1 = 4;
    _material[_STEEL_]->thermal_expansion ( 1 , &tmp1 , &tmp2 );
    double etmp[18] = { 0.00304 , -1 , 0.003051 , 0.003051 , 0.003021 ,
			0.00291 , 0.00274 , 0.00263 , 0.00239 , 0.00212 ,
			0.00198 , 0.00167 , 0.00137 , 0.00121 , 0.00088 ,
			0.00053 , 0.00035 , 0 };
    etmp[1] = tmp2;
    _material[_STEEL_]->set_e ( 18 , etmp );
    _material[_STEEL_]->set_ekABC2();

    _material[_STEEL_]->set_ksABC();
    _material[_STEEL_]->set_TEI();

    _material[_STEEL_]->set_y ( 18 , y6 );
    _material[_STEEL_]->set_yABC();
  }

  // carbon steel:
  {
    _material[_CARBON_STEEL_] = new HS_Material ( _CARBON_STEEL_ , 0.282 );
    _material[_CARBON_STEEL_]->set_T ( 18 , T1 );
    _material[_CARBON_STEEL_]->set_k ( 18 , k7 );
    _material[_CARBON_STEEL_]->set_s ( 18 , s7 );
    _material[_CARBON_STEEL_]->set_kABC();
    _material[_CARBON_STEEL_]->set_sABC();
    _material[_CARBON_STEEL_]->set_TCI();

    _material[_CARBON_STEEL_]->set_e ( 17 , e7 );
    _material[_CARBON_STEEL_]->set_ekABC1();
    tmp1 = 4;
    _material[_CARBON_STEEL_]->thermal_expansion ( 1 , &tmp1 , &tmp2 );
    double etmp[18] = { 0.0021 , -1 , 0.0021 , 0.00209 , 0.00207 , 0.00201 ,
			0.0019 , 0.00183 , 0.00168 , 0.00151 , 0.00141 ,
			0.0012 , 0.00098 , 0.00087 , 0.00064 , 0.00036 ,
			0.00025 , 0 };

    etmp[1] = tmp2;
    _material[_CARBON_STEEL_]->set_e ( 18 , etmp );
    _material[_CARBON_STEEL_]->set_ekABC2();

    _material[_CARBON_STEEL_]->set_ksABC();
    _material[_CARBON_STEEL_]->set_TEI();

    _material[_CARBON_STEEL_]->set_y ( 18 , y7 );
    _material[_CARBON_STEEL_]->set_yABC();
  }
}

/*--------------*/
/*  destructor  */
/*--------------*/
HS_Param::~HS_Param ( void ) 
 {
  for ( int t = _NYLON_ ; t <= _CARBON_STEEL_ ; ++t )
    delete _material[t];
  delete [] _material;
}
