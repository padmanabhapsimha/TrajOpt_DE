#ifndef _3bodyfunction_hpp
#define _3bodyfunction_hpp

#include"agent_defn.hpp"
#include"dataTYPE.hpp"

#include<cmath>
#include<vector>
#include<boost/numeric/odeint.hpp>

Agent_datatype _3body_cost(sys_pars<Agent_datatype> &params,const std::vector<Agent_datatype> &vals,const bool printval,
                           Agent_datatype &massfrac);


#endif // _3bodyfunction_hpp
