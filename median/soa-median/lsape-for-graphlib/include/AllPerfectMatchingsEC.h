/**
 * @file AllPerfectMatchings.h
 * @author Évariste <<21511575@etu.unicaen.fr>>
 */

#ifndef __ALLPERFECTMATCHINGSEC_H__
#define __ALLPERFECTMATCHINGSEC_H__


#include "AllPerfectMatchings.h"

/**
 * @class AllPerfectMatchingsEC
 * @autor evariste
 * @brief Error correcting version of AllPerfectMatchings
 *
 *   This class allows to avoid all epsilon-to-epsilon matching enumeration
 *   in an error correcting version of the LSAP. 
 */
template <typename IT>
class AllPerfectMatchingsEC : public AllPerfectMatchings<IT> 
{

protected:

  IT _n; //< Number of elements in G1
  IT _m; //< Number of elements in G2

protected:
  
  virtual bool
  chooseEdge( const cDigraph<IT>& gm, IT* x, IT* y );
  
  
public :

  AllPerfectMatchingsEC(cDigraph<IT>& gm) : 
    AllPerfectMatchings<IT>(gm),
    _n(gm.rows()/2),
    _m(gm.cols()/2)
  {}
  
  
  AllPerfectMatchingsEC(cDigraph<IT>& gm, const IT& n, const IT& m ) : 
    AllPerfectMatchings<IT>(gm),
    _n(n), _m(m)
  {}
  
  
};



/*********** Implementation ***************/


template <typename IT>
bool
AllPerfectMatchingsEC<IT>::chooseEdge( const cDigraph<IT>& gm, IT* x, IT* y ) 
{
  for (IT i=0; i<gm.rows(); i++){
    for (IT j=0; j<gm.cols(); j++){
      if((i<_n || j<_m) && gm(i,j) == 1){
        *x = i;
        *y = j;
        return true;
      }
    }
  }
  return false;
}



#endif
