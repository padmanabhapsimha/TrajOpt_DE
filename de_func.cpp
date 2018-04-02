/**< ******************************************************** */
/**< ***DIFFERENTIAL EVOLUTION MINIMIZATION CODE ************ */
/**< ******************************************************** */
/**< ***DEVELOPED BY PADMANABHA PRASANNA SIMHA (SC14B034)**** */
/**< ***UNDER THE GUIDANCE OF DR. R. V. RAMANAN************** */
/**< ***FINAL SEMESTER BTECH PROJECT FOR LOW THRUST TRANSFERS */
/**< ***INDIAN INSTITUTE OF SPACE SCIENCE AND TECHNOLOGY***** */
/**< ***C++ 14 - STL AND BOOST LIBRARIES ARE USED************ */
/**< ***USE WITH 64BIT GCC COMPILER & LINK BOOST LIBRARIES*** */
/**< ***AS OF NOW, CODE IS THREAD SAFE, TESTED & VERIFIED**** */
/**< ***READ STL AND BOOST DOCUMENTATION FOR FURTHER INFO**** */
/**< ***RUN IN MULTI THREADED MODE FOR LARGE AGENT COUNTS**** */
/**< ***ELSE SINGLE THREADED MODE IS BETTER****************** */
/**< ***DUE TO BOOST LIBRARIES, LARGE COMPILATION TIME******* */
/**< ******************************************************** */
#include"de_headers.hpp"
using namespace std;
void de_function(string inputfile,string agentoutputfile,string bestagentfile,string sysparamsfile)
{
    ///RANDOM NUMBER GENERATOR TYPE
    typedef boost::mt19937 RNGType;
    ///DE PARAMS
    unsigned int gens,dim,agents_no=0;
    unsigned int threads_no;
    int cost_select,seedval,agtinit,rebound;
    string agtinputfile,bestagtfile;
    Agent_datatype CR,F;
    vector<Agent_datatype> L,H,newboundlimit;
    ///READ ALL THESE FROM FILE
    params_reader(gens,dim,agents_no,CR,F,L,H,cost_select,threads_no,inputfile,seedval,agtinit,agtinputfile,
                  bestagtfile,rebound,newboundlimit);
    if(seedval<0){
        seedval=static_cast<unsigned int>(clock());
    }
    else{
        seedval=static_cast<unsigned int>(seedval);
    }
    /// /////////////////////////////////////////////////////
    ///RANDOM NUMBERS GENERATORS
    RNGType RNG(seedval);
    boost::uniform_real<> uni_dist(0,1);
    boost::variate_generator<RNGType,boost::uniform_real<>>randgen_real(RNG,uni_dist);
    RNGType RNG1(1+seedval);
    boost::uniform_int<> uni_distint(0,agents_no-1);
    boost::variate_generator<RNGType,boost::uniform_int<>>randgen_int(RNG1,uni_distint);
    RNGType RNG2(79+seedval);
    boost::uniform_int<> uni_distint2(0,dim-1);
    boost::variate_generator<RNGType,boost::uniform_int<>>randgen_int2(RNG2,uni_distint2);
    /// /////////////////////////////////////////////////////
    /**< FILL UP AGENTS WITH RANDOM VALUES IN GIVEN HYPERBOX */
    typedef Agent<Agent_datatype> MyAgent;
    vector<MyAgent> allAgents;
    vector<Agent_datatype> holder;
    for(unsigned int i=0;i<agents_no;i++){
        holder.clear();
        for(unsigned int j=0;j<dim;j++){
            holder.push_back(L.at(j)+(H.at(j)-L.at(j))*randgen_real());
        }
        allAgents.push_back(MyAgent(holder));
    }
    ///REWRITE AGENTS WITH FILE INPUT IF DEEMED SO
    double inputval;
    int agtnumber;
    fstream fileinput(agtinputfile,ios::in);
    if(!fileinput.good()){
        cout<<endl<<"Something bad has happened. Unable to open "<<agtinputfile<<endl;
        exit(-1);
    }
    if(agtinit){
        for(unsigned int i=0;i<agents_no;i++){
            fileinput>>agtnumber;///junk data -- agent number
            for(unsigned int j=0;j<dim;j++){
                fileinput>>inputval;///useful data
                allAgents.at(i).vals.at(j)=inputval;
            }
            fileinput>>inputval;///junk data -- agent cost value
            allAgents.at(i).cost_val=inputval;
        }
    }
    ///RECENTER BOUNDS IF NEEDED
    fileinput.close();
    fileinput.open(bestagtfile,ios::in);
    if(!fileinput.good()){
        cout<<endl<<"Something bad has happened. Unable to open "<<agtinputfile<<endl;
        exit(-1);
    }
    if(agtinit){
        for(unsigned int j=0;j<dim;j++){
            fileinput>>inputval;///junk
        }
        for(unsigned int j=0;j<dim;j++){
            fileinput>>inputval;///lower bounds
            L.at(j)=inputval;
        }
        for(unsigned int j=0;j<dim;j++){
            fileinput>>inputval;///upper bounds
            H.at(j)=inputval;
        }
    }
    fileinput.close();
    fileinput.open(bestagtfile,ios::in);
    if(rebound){
        for(unsigned int i=0;i<dim;i++){
            fileinput>>inputval;///useful data
            L.at(i)=inputval-newboundlimit.at(i)*abs(inputval);
            H.at(i)=inputval+newboundlimit.at(i)*abs(inputval);
        }
        ///reseed agents within bounds
        for(unsigned int i=0;i<agents_no;i++){
            holder.clear();
            for(unsigned int j=0;j<dim;j++){
                holder.push_back(L.at(j)+(H.at(j)-L.at(j))*randgen_real());
            }
            allAgents.at(i)=holder;
        }
    }
    ///INITIALIZE AGENTS WITH SYSTEM PARAMETERS
    for(unsigned int i=0;i<agents_no;i++){
        allAgents.at(i).params_init(sysparamsfile);
        allAgents.at(i).cost_fn(cost_select);
    }
    /**< START OF DIFFERENTIAL EVOLUTION */
    mutex myMutex;
    chrono::high_resolution_clock::time_point t1=chrono::high_resolution_clock::now();
    ///SINGLE THREADED VERSION
    if(threads_no==1){
        de_looper(agents_no,gens,dim,CR,F,L,H,allAgents,cost_select,1,myMutex,0,agents_no,1,seedval);
    }
    ///MULTI THREADED VERSION - SHARED MEMORY PARALLEL
    ///OPENMP IS TOO MUCH OF A BLACK BOX FOR ME TO BE COMFORTABLE WITH
    ///C++ 14 STANDARD THREADS ARE USED
    else{
        de_threader(threads_no,agents_no,gens,dim,CR,F,L,H,allAgents,cost_select,myMutex,seedval);
    }
    chrono::high_resolution_clock::time_point t2=chrono::high_resolution_clock::now();
    chrono::duration<double> time_span=chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    /// /////////////////////////////////////////////////////
    /**< WRITE ALL AGENTS TO FILE */
    fstream fout(agentoutputfile,ios::out);
    fout.precision(18);
    for(unsigned int i=0;i<agents_no;i++){
        fout<<i<<"\t";
        for(unsigned int j=0;j<dim;j++){
            fout<<(allAgents.at(i).vals).at(j)<<"\t";
        }
        fout<<allAgents.at(i).cost_val<<"\n";
    }
    fout.close();
    /// /////////////////////////////////////////////////////
    /**< FIND OUT BEST AGENT AND DISPLAY DATA TO SCREEN */
    vector<Agent_datatype> costvec;
    for(unsigned int i=0;i<agents_no;i++){
        costvec.push_back(allAgents.at(i).cost_val);
    }
    auto itr=min_element(costvec.begin(),costvec.end());
    unsigned int optagtno=static_cast<int>(itr-costvec.begin());//distance(costvec.begin(),itr);
    cout<<endl<<"# Agents = "<<agents_no<<endl<<"Dimension = "<<dim;
    cout<<endl<<"CR = "<<CR<<endl<<"F = "<<F<<endl<<"Generations = "<<gens;
    cout<<endl<<"Number of threads used = "<<threads_no<<endl;
    cout<<endl<<"Problem number = "<<cost_select<<endl;
    cout.precision(numeric_limits<Agent_datatype>::max_digits10);
    cout<<endl<<endl<<"Agent #"<<optagtno<<"; cost = "<<allAgents.at(optagtno).cost_val<<endl;
    cout<<endl<<allAgents.at(optagtno).massfrac<<endl;
    cout.precision(10);
    fstream optagtout(bestagentfile,ios::out);
    optagtout.precision(numeric_limits<Agent_datatype>::max_digits10);
    for(unsigned int i=0;i<dim;i++){
        cout<<allAgents.at(optagtno).vals.at(i)<<"\t";
        optagtout<<allAgents.at(optagtno).vals.at(i)<<"\t";
    }
    optagtout<<endl;
    for(unsigned int i=0;i<dim;i++){
        optagtout<<L.at(i)<<"\t";
    }
    optagtout<<endl;
    for(unsigned int i=0;i<dim;i++){
        optagtout<<H.at(i)<<"\t";
    }
    allAgents.at(optagtno).cost_fn_print(cost_select);///PRINT OUTPUT IF ANY TO FILE
    /// /////////////////////////////////////////////////////
    cout<<endl<<endl<<endl<<"Time elapsed for solution = "<<time_span.count()<<" s\n\n\n";
}
