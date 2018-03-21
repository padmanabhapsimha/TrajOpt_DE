#include"agent_defn.hpp"
#include"dataTYPE.hpp"

#include<vector>
#include<algorithm>
#include<iterator>
#include<ctime>
#include<thread>
#include<mutex>
#include<cstdlib>
#include<cmath>
#include<string>
#include<fstream>

#include<boost/random.hpp>
#include<boost/generator_iterator.hpp>
/**< ******************************************** */
template<typename Agent_data>
void de_looper(const unsigned int agents_no,const unsigned int gens,const unsigned int dim,const Agent_data CR,
               const Agent_data F,const std::vector<Agent_data> L,const std::vector<Agent_data> H,
               std::vector<Agent<Agent_data>> &allAgents,const int cost_select,const int thread_no,std::mutex &myMutex,
               const unsigned int low,const unsigned int high,const unsigned int threadid,const unsigned int seedval)
{
    std::string prefix="Populations\\";
    std::string suffix=".txt";
    std::fstream fout(prefix+std::to_string(threadid)+suffix,std::ios::out);
    fout.close();
    /**< RANDOM NUMBER GENERATORS */
    typedef boost::mt19937 RNGType;
    /// /////////////////////////////
    /**< WHOLE BUNCH OF RANDOM NUMBER GENERATORS INDEPENDENTLY SEEDED */
    RNGType RNG(thread_no+seedval);
    boost::uniform_real<> uni_dist(0,1);
    boost::variate_generator<RNGType,boost::uniform_real<>>randgen_real(RNG,uni_dist);
    RNGType RNG2(5*thread_no+17+seedval);
    boost::uniform_int<> uni_distint2(0,dim-1);
    boost::variate_generator<RNGType,boost::uniform_int<>>randgen_int2(RNG2,uni_distint2);
    RNGType RNG3(10*thread_no+34+seedval);
    boost::uniform_int<> uni_distint3(0,agents_no-2);
    boost::variate_generator<RNGType,boost::uniform_int<>>randgen_int3(RNG3,uni_distint3);
    RNGType RNG4(11*thread_no+341+seedval);
    boost::uniform_int<> uni_distint4(0,agents_no-3);
    boost::variate_generator<RNGType,boost::uniform_int<>>randgen_int4(RNG4,uni_distint4);
    RNGType RNG5(12*thread_no+1235+seedval);
    boost::uniform_int<> uni_distint5(0,agents_no-4);
    boost::variate_generator<RNGType,boost::uniform_int<>>randgen_int5(RNG5,uni_distint5);
    /// ////////////////////////////
    /**< PREPARATIONS FOR FISCHER-YATES TYPE SHUFFLING */
    Agent<Agent_data> Agent_x,Agent_y,Agent_a,Agent_b,Agent_c;
    Agent_x=allAgents.at(0);
    Agent_y=Agent_x;
    unsigned int ind_a,ind_b,ind_c;
    /**< INDEX VECTOR FOR THE FISCHER-YATES TYPE SHUFFLE */
    std::vector<unsigned int> indices;
    for(unsigned int i=0;i<agents_no;i++){
        indices.push_back(i);
    }
/**< LAMBDA EXPRESSION FOR EASE OF CODING - INSTEAD OF SEPARATE NATIVE FUNCTION*/
auto agent_selector=[&](const unsigned int &agntno){
    /**< ***FISCHER-YATES TYPE SHUFFLE*** */
    /**< ***MARGINALLY BETTER THAN NAIVE METHOD FOR THIS SMALL PROBLEM*** */
    /**< IN CASE MORE DISTINCT RANDOM AGEENTS ARE NEEDED, THIS WILL OUTPERFORM THE NAIVE TECHNIQUE */
    std::swap(indices.at(agntno),(indices.at(agents_no-1)));
    ind_a=randgen_int3();///rand numbers from 0 to np-1
    std::swap(indices.at(ind_a),(indices.at(agents_no-2)));
    ind_b=randgen_int4();///rand numbers from 0 to np-2
    std::swap(indices.at(ind_b),(indices.at(agents_no-3)));
    ind_c=randgen_int5();///rand numbers from 0 to np-3
    ///UNSWAP INDICES VECTOR IN REVERSE ORDER FOR NEXT USE
    std::swap(indices.at(ind_b),(indices.at(agents_no-3)));
    std::swap(indices.at(ind_a),(indices.at(agents_no-2)));
    std::swap(indices.at(agntno),(indices.at(agents_no-1)));
    /**< BOTH ARE ALMOST EQUAL IN PERFORMANCE BUT FISCHER-YATES IS EASILY EXTENSIBLE
         FISHER-YATES IS ALSO OF CONSTANT COMPLEXITY COMPARED TO THE NAIVE TECHNIQUE
     */
/// ////////////////////////////
    Agent_a=std::ref(allAgents.at(ind_a));
    Agent_b=std::ref(allAgents.at(ind_b));
    Agent_c=std::ref(allAgents.at(ind_c));
};
/**< END OF LAMBDA EXPRESSION */
    /**< START OF ACTUAL DIFFERENTIAL EVOLUTION LOOPS */
    /**< ALL THE COST FUNCTION EVALS ALONG WITH INITIALIZATIONS IS DONE INDEPENDENTLY PER THREAD */
    bool repick_flag=false;
    for(unsigned int j=low;j<high;j++){///INITIALIZE ALL THE AGENTS WITH THEIR COST VALUES
        allAgents.at(j).cost_init(cost_select);
    }
    for(unsigned int i=0;i<gens;i++){///LOOP OVER GENERATIONS
        for(unsigned int j=low;j<high;j++){///LOOP OVER GIVEN AGENTS
///GOTO MAY PREVENT LOOP VECTORIZATION
///CURRENT IMPLEMENTATION IS JUST AS FAST AS THE GOTO STUFF
///re_pick:///TO REDO THE AGENT SELECTION IF AGENT Y GOES OUT OF BOUNDING HYPERBOX
            Agent_x=allAgents.at(j);repick_flag=false;
            agent_selector(j);///SELECT 3 DISTINCT AGENTS FROM Agent_x
            for(unsigned int k=0;k<dim;k++){///LOOP OVER DIMENSIONS - IE COMPONENTS OF THE AGENT
                if((randgen_real()<=CR)||(static_cast<unsigned int>(randgen_int2())==k)){///MUTATE AND CROSSOVER
                    Agent_y.vals.at(k)=Agent_a.vals.at(k)+F*(Agent_b.vals.at(k)-Agent_c.vals.at(k));
                    if(L.at(k)>Agent_y.vals.at(k)||H.at(k)<Agent_y.vals.at(k)){
                        repick_flag=true;break;///IF ANY COMPONENT GOES OUT OF HYPERBOX, PICK AGENTS AGAIN
                        ///goto re_pick;///DON'T USE THIS
                    }
                }
                else{
                    Agent_y.vals.at(k)=Agent_x.vals.at(k);
                }
//                if(repick_flag){break;}
            }
            if(repick_flag){
                j?(--j):j;///TO PREVENT INDEX FROM GOING BELOW ZERO - MAYBE NOT NEEDED - DIDN'T THINK ABOUT THIS TOO MUCH
                continue;///REDO THIS LOOP ITERATION -> LOOP INDEX IS DECREMENTED
                ///IS ALL THIS ACTUALLY BETTER THAN A SIMPLE GOTO???? MAYBE NOT
            }
            ///SELECTION STEP - DE HAS INHERENT ELITISM TYPE OF SELECTION
            ///1 COST FUNCTION EVALUATION PER AGENT PER GENERATION
            if(Agent_y.cost_fn(cost_select)<=Agent_x.cost_val)
//            if(Agent_y.cost_fn(cost_select)<=Agent_x.cost_fn(cost_select))
            {
                std::lock_guard<std::mutex> lock(myMutex);///PREVENT DATA RACES AND CONFLICTS
                allAgents.at(j)=Agent_y;///OVERWRITE THE OLD AGENT WITH NEW AGENT
            }
        }
        /// ///////////////////////////////////////
        ///PRINT AGENTS IN RESPECTIVE FILES
        /**
        if(i%10==0){
            std::lock_guard<std::mutex> lock(myMutex);///PREVENT DATA RACES AND CONFLICTS
            std::fstream fout(prefix+std::to_string(threadid)+suffix,std::ios::app);
            fout.precision(17);
            fout<<"Generation#\t"<<i+1<<std::endl;
            for(unsigned int i=0;i<agents_no;i++){
                fout<<i<<"\t";
                for(unsigned int j=0;j<dim;j++){
                    fout<<(allAgents.at(i).vals).at(j)<<"\t";
                }
                fout<<allAgents.at(i).cost_val<<"\n";
            }
            std::vector<Agent_data> costvec;
            for(unsigned int i=0;i<agents_no;i++){
                costvec.push_back(allAgents.at(i).cost_fn(cost_select));
            }
            auto itrmin=std::min_element(costvec.begin(),costvec.end());
            auto itrmax=std::max_element(costvec.begin(),costvec.end());
            unsigned int minagtno=static_cast<int>(itrmin-costvec.begin());
            unsigned int maxagtno=static_cast<int>(itrmax-costvec.begin());
            fout<<"Min member\t"<<minagtno<<"\t"<<allAgents.at(minagtno).cost_fn(cost_select)<<std::endl;
            fout<<"Max member\t"<<maxagtno<<"\t"<<allAgents.at(maxagtno).cost_fn(cost_select)<<"\t"<<std::endl<<std::endl<<std::endl;
            fout.close();
        }
        **/
        /// ///////////////////////////////////////
    }
}
/**< ******************************************** */
/**< ******************************************** */
template<typename Agent_data>
void de_threader(const unsigned int total_threads,const unsigned int agents_no,const unsigned int gens,const unsigned int dim,
                 const Agent_data CR,const Agent_data F,const std::vector<Agent_data> &L,const std::vector<Agent_data> &H,
                 std::vector<Agent<Agent_data>> &allAgents,const int cost_select,std::mutex &myMutex,const unsigned int &seedval)
{
    std::vector<std::thread> Tvector;
    std::vector<unsigned int> l;
    std::vector<unsigned int> h;
    for(unsigned int i=0;i<total_threads;i++){
        l.push_back(static_cast<unsigned int>(i*agents_no/total_threads));
        h.push_back(static_cast<unsigned int>((i+1)*agents_no/total_threads));
    }
    l.at(0)=0;
    h.at(total_threads-1)=agents_no;///TO ADJUST FOR NON DIVIDING CASES
    /**< PUSH BACK THREADS AND WAIT FOR COMPLETION */
    for(unsigned int i=0;i<total_threads;i++){
        Tvector.push_back(std::thread(&de_looper<Agent_data>,agents_no,gens,dim,CR,F,L,H,
                                      std::ref(allAgents),cost_select,1,std::ref(myMutex),l.at(i),h.at(i),i+1,seedval));
    }
    /**< JOIN THE THREADS */
    for(unsigned int i=0;i<Tvector.size();i++){
        if(Tvector.at(i).joinable()){Tvector.at(i).join();}
    }
}
/**< ******************************************** */
///TEMPLATE FUNCTION INSTANTIATIONS
template void de_threader(const unsigned int total_threads,const unsigned int agents_no,const unsigned int gens,
                          const unsigned int dim,const Agent_datatype CR,const Agent_datatype F,
                          const std::vector<Agent_datatype> &L,const std::vector<Agent_datatype> &H,
                          std::vector<Agent<Agent_datatype>> &allAgents,const int cost_select,std::mutex &myMutex,
                          const unsigned int &seedval);

template void de_looper(const unsigned int agents_no,const unsigned int gens,const unsigned int dim,const Agent_datatype CR,
                        const Agent_datatype F,const std::vector<Agent_datatype> L,const std::vector<Agent_datatype> H,
                        std::vector<Agent<Agent_datatype>> &allAgents,const int cost_select,const int thread_no,
                        std::mutex &myMutex,const unsigned int low,const unsigned int high,const unsigned int threadid,
                        const unsigned int seedval);
/**< ******************************************** */

