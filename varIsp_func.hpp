#ifndef _varIsp_func_hpp
#define _varIsp_func_hpp

#include"agent_defn.hpp"
#include"dataTYPE.hpp"

#include<cmath>
#include<vector>
#include<boost/numeric/odeint.hpp>

Agent_datatype _varIsp_cost(sys_pars<Agent_datatype> &params,const std::vector<Agent_datatype> &vals,const bool printval,
                           Agent_datatype &massfrac);


#endif // _varIsp_func_hpp
