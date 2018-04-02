#ifndef _PARAMS_READER_HPP
#define _PARAMS_READER_HPP
#include<fstream>
/**< ******************************************** */
/**< CODE IS SELF EXPLANATORY */
template<typename Agent_data>
void params_reader(unsigned int &generations,unsigned int &dimensions,unsigned int &agents,Agent_data &CR,Agent_data &F,
                   std::vector<Agent_data> &L,std::vector<Agent_data> &H,int &cost_select,unsigned int &threads_no,
                   std::string filename,int &seedval,int &agtinit,std::string &agtinputfile,std::string &bestagtfile,
                   int &rebound,std::vector<Agent_data> &newboundlimit);
/**< ******************************************** */
#endif // _PARAMS_READER_HPP
