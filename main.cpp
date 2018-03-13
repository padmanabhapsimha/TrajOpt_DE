#include"de_func.hpp"
#include"Wrapper_JPL_Ephem.hpp"
#include<iostream>
#include<array>
#include<cmath>
#include"jpleph.h"


using namespace std;
/**< ******************************************** */
int main()
{
    de_function("de_input.txt","Agent_output.txt","Best_agent.txt","System_params.txt");


    return 0;
}
