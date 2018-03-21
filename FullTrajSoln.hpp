#ifndef _FullTrajSoln_hpp
#define _FullTrajSoln_hpp

#include"agent_defn.hpp"
#include"dataTYPE.hpp"

#include<cmath>
#include<vector>
#include<boost/numeric/odeint.hpp>

Agent_datatype _fullSoln_cost(sys_pars<Agent_datatype> &params,const std::vector<Agent_datatype> &vals,const bool printval,
                           Agent_datatype &massfrac);


#endif // _FullTrajSoln_hpp

