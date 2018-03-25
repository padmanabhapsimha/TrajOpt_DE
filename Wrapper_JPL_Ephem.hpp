/**    ntarg = integer number of 'target' point.                             **
**                                                                          **
**    ncent = integer number of center point.                               **
**                                                                          **
**    The numbering convention for 'ntarg' and 'ncent' is:                  **
**                                                                          **
**            1 = mercury           8 = neptune                             **
**            2 = venus             9 = pluto                               **
**            3 = earth            10 = moon                                **
**            4 = mars             11 = sun                                 **
**            5 = jupiter          12 = solar-system barycenter             **
**            6 = saturn           13 = earth-moon barycenter               **
**            7 = uranus           14 = nutations (longitude and obliq)     **
**                                 15 = librations, if on eph. file         **
**                                 16 = lunar mantle omega_x,omega_y,omega_z**
**                                 17 = TT-TDB, if on eph. file             **
**                                                                          **
**            (If nutations are wanted, set ntarg = 14.                     **
**             For librations, set ntarg = 15. set ncent= 0.                **
**             For TT-TDB,  set ntarg = 17.  I've not actually              **
**             seen an ntarg = 16 case yet.)                                **
**                                                                          **
**/

#ifndef __WRAPPER_JPL__
#define __WRAPPER_JPL__

#include<array>
#include<string>
#include<cmath>
#include<iostream>
#include"jpleph.h"

class wrappy_jpl{
    void* JPLeph;
public:
    wrappy_jpl(){///default constructor
        JPLeph=nullptr;
    }
    void wrappy_jpl_default(){
        JPLeph=jpl_init_ephemeris("linux_p1550p2650.430",nullptr,nullptr);
    }
    wrappy_jpl(const std::string ephemfile){///init constructor
        JPLeph=jpl_init_ephemeris(ephemfile.c_str(),nullptr,nullptr);
    }
    wrappy_jpl(const wrappy_jpl &inputobj){
        jpl_close_ephemeris(JPLeph);
        JPLeph=inputobj.JPLeph;
    }///copy constructor
    void JPL_ecliptic_m_get(const double &epoch,const int &ntarg,const int &ncent,std::array<double,6> &state,
                            const int &vCompute);///STATE VECTOR IN ECLIPTIC FRAME
    void JPL_equatorial_m_get(const double &epoch,const int &ntarg,const int &ncent,std::array<double,6> &state,
                            const int &vCompute);///STATE VECTOR IN EQUATORIAL FRAME
    ~wrappy_jpl(){
        jpl_close_ephemeris(JPLeph);
    }///destructor
};


#endif // __WRAPPER_JPL__

