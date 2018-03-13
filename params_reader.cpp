#include"agent_defn.hpp"
#include"dataTYPE.hpp"

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>

template<typename Agent_data>
void params_reader(unsigned int &generations,unsigned int &dimensions,unsigned int &agents,Agent_data &CR,Agent_data &F,
                   std::vector<Agent_data> &L,std::vector<Agent_data> &H,int &cost_select,unsigned int &threads_no,
                   std::string filename,int &seedval)
{
    std::fstream fin(filename,std::ios::in);
    Agent_data Ltemp,Htemp;
    if(fin.good())
    {
        fin>>generations>>dimensions>>agents>>CR>>F;
        std::string selector;
        fin>>selector;
        if(selector=="carry_toAll"){
            fin>>Ltemp>>Htemp;
            for(unsigned int i=0;i<dimensions;i++){
                L.push_back(Ltemp);H.push_back(Htemp);
            }
        }
        else{
            for(unsigned int i=0;i<dimensions;i++){
                fin>>Ltemp;L.push_back(Ltemp);
            }
            for(unsigned int i=0;i<dimensions;i++){
                fin>>Htemp;H.push_back(Htemp);
            }
        }
        fin>>cost_select>>threads_no>>seedval;
    }
    else
    {
        std::cout<<std::endl<<"Error opening input file \""<<filename<<"\""<<std::endl;
        exit(-101);
    }
    fin.close();
}


///TEMPLATE FUNCTION INSTANTIATION
template void params_reader(unsigned int &generations,unsigned int &dimensions,unsigned int &agents,Agent_datatype &CR,
                            Agent_datatype &F,std::vector<Agent_datatype> &L,std::vector<Agent_datatype> &H,
                            int &cost_select,unsigned int &threads_no,std::string filename,int &seedval);
