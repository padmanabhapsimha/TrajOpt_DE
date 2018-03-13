#include<vector>
#include<string>
#include<fstream>
#include<cstdlib>
#include<iostream>
#include"Wrapper_JPL_Ephem.hpp"

using namespace std;

void wrappy_jpl::JPL_ecliptic_m_get(const double &epoch,const int &ntarg,const int &ncent,std::array<double,6> &state,
                                    const int &vCompute)
{
        jpl_pleph(JPLeph,epoch,ntarg,ncent,state.begin(),vCompute);
    //      constexpr double e=23.439291111111111111111111111111*3.1415926535897932384626433832795028/180.0;
        constexpr double e=0.40909280422232893747235670152584;
        auto posnew=state;
        posnew.at(1)=cos(e)*state.at(1)+sin(e)*state.at(2);
        posnew.at(2)=-sin(e)*state.at(1)+cos(e)*state.at(2);
        posnew.at(4)=cos(e)*state.at(4)+sin(e)*state.at(5);
        posnew.at(5)=-sin(e)*state.at(4)+cos(e)*state.at(5);
        state=posnew;
        state.at(0)*=149597870700.0;
        state.at(1)*=149597870700.0;
        state.at(2)*=149597870700.0;
        state.at(3)*=1731456.8368055555555555555555556;
        state.at(4)*=1731456.8368055555555555555555556;
        state.at(5)*=1731456.8368055555555555555555556;
}
