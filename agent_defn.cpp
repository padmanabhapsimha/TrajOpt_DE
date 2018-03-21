

#include"agent_defn.hpp"
#include"dataTYPE.hpp"
#include<vector>
#include<string>
#include<fstream>
#include<cstdlib>
#include<iostream>
#include"Wrapper_JPL_Ephem.hpp"
template<typename Agent_data>
inline Agent<Agent_data>::Agent(){;}///DEFAULT CONSTRUCTOR

///COPY CONSTRUCTOR - C++ MAKES IT'S OWN - NO NEED TO WRITE THIS
template<typename Agent_data>
inline Agent<Agent_data>::Agent(const Agent &input_agent){
    vals=input_agent.vals;cost_val=input_agent.cost_val;
    params=input_agent.params;
}

///SOME OTHER CONSTRUCTOR
template<typename Agent_data>
inline Agent<Agent_data>::Agent(const std::vector<Agent_data> &input_vec){
    if(!input_vec.empty()){
        for(auto itr=input_vec.begin();itr!=input_vec.end();++itr){
            vals.push_back(*itr);
        }
    }
}

///INITIALIZE AGENTS WITH THEIR COSTS
template<typename Agent_data>
inline void Agent<Agent_data>::cost_init(const int &cost_select){
    cost_val=cost_fn(cost_select);
}


using namespace std;
///INITIALIZE SYSTEM PARAMETERS WITH THEIR VALUES
template<typename Agent_data>
inline void Agent<Agent_data>::params_init(const string inpfile){
    fstream fin(inpfile,ios::in);
    if(fin.good()){
        fin>>params.mu>>params.mi>>params.Isp>>params.g0>>params.kfactor>>params.rpfac>>params.ecc>>params.nu
        >>params.rinit>>params.rf_fac>>params.NEP>>params.eccf>>params.omegaf>>params.maxtime
        >>params.a_tol>>params.integ_tol>>params.maxstep>>params.termination_type>>params.inclf>>params.raanf
        >>params.raani>>params.omegai>>params.incli>>params.mu_e>>params.startbody>>params.ephem_file
        >>params.a_second>>params.julian_start>>params.rfinal>>params.endbody>>params.printframe>>params.maxstepcount
        >>params.mu_m>>params.a_third;
        params.raanf*=3.1415926535897932384626433832795028/180.0;
        params.inclf*=3.1415926535897932384626433832795028/180.0;
        params.omegaf*=3.1415926535897932384626433832795028/180.0;
        params.raani*=3.1415926535897932384626433832795028/180.0;
        params.omegai*=3.1415926535897932384626433832795028/180.0;
        params.incli*=3.1415926535897932384626433832795028/180.0;
    }
    else{
        cout<<endl<<"Unable to open params file";
        exit(-1);
    }
    fin.close();
}

