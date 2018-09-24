#ifndef __HS_PARAM__
#define __HS_PARAM__

#include "HS_Material.hpp"

class HS_Param {

private:

  int            _load_flag;
  double         _load;
  double         _max_weight;
  double         _TE_limit;
  double         _T_cold;
  double         _T_hot;
  double         _L;
  HS_Material ** _material;

public:

  // constructor:
  HS_Param ( int    load_flag  ,
	     double load       ,
	     double max_weight ,
	     double TE_limit   ,
	     double T_cold     ,
	     double T_hot      ,
	     double L            );

  // destructor:
  ~HS_Param ( void );

  // GET methods:
  int           get_load_flag  ( void  ) const { return _load_flag;   }
  double        get_load       ( void  ) const { return _load;        }
  double        get_max_weight ( void  ) const { return _max_weight;  }
  double        get_TE_limit   ( void  ) const { return _TE_limit;    }
  double        get_T_cold     ( void  ) const { return _T_cold;      }
  double        get_T_hot      ( void  ) const { return _T_hot;       }
  double        get_L          ( void  ) const { return _L;           }
  HS_Material * get_material   ( int m ) const { return _material[m]; }
};

#endif
