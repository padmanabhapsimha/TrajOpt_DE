#ifndef _DE_SOLVER_HPP
#define _DE_SOLVER_HPP
#include<vector>
#include<mutex>
/**< ******************************************** */
template<typename Agent_data>
void de_looper(const unsigned int agents_no,const unsigned int gens,const unsigned int dim,const Agent_data CR,
               const Agent_data F,const std::vector<Agent_data> L,const std::vector<Agent_data> H,
               std::vector<Agent<Agent_data>> &allAgents,const int cost_select,const int thread_no,std::mutex &myMutex,
               const unsigned int low,const unsigned int high,const unsigned int threadid,const unsigned int seedval);
/**< ******************************************** */
template<typename Agent_data>
void de_threader(const unsigned int total_threads,const unsigned int agents_no,const unsigned int gens,const unsigned int dim,
                 const Agent_data CR,const Agent_data F,const std::vector<Agent_data> &L,const std::vector<Agent_data> &H,
                 std::vector<Agent<Agent_data>> &allAgents,const int cost_select,std::mutex &myMutex,const unsigned int &seedval);
/**< ******************************************** */
#endif // _DE_SOLVER_HPP
