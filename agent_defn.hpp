#ifndef _AGENT_DEFN_HPP
#define _AGENT_DEFN_HPP
/**< ******************************************** */
#include<vector>
#include<string>
#include<thread>
#include"Wrapper_JPL_Ephem.hpp"
///THIS STRUCT HOLDS DATA OF SYSTEM PARAMETERS FOR THE AGENT
template<typename Agent_data>
struct sys_pars{
public:
    Agent_data mu;///gravitational parameter of first major body
    Agent_data mi;///initial spacecraft mass
    Agent_data Isp;///isp of propulsion system
    Agent_data g0;///acceleration due to gravity - required for propulsion system performance
    Agent_data kfactor;///or max thrust
    Agent_data rpfac;///initial orbit periapsis multiplication factor for rinit
    Agent_data ecc;///initial orbit eccentricity
    Agent_data nu;///starting true anomaly
    Agent_data rinit;///periapsis radius of initial orbit
    Agent_data rf_fac;///final orbit periapsis radius multiplying factor with rinit
    int NEP;
    Agent_data eccf;///final orbit eccentricity
    Agent_data omegaf;///final orbit argument of periapsis
    Agent_data maxtime;///maximum integration time for some integration terminations modes
    Agent_data a_tol;///tolerance for integration termination in semi major axis
    Agent_data integ_tol;///tolerance for integration step size control
    Agent_data maxstep;///maximum step size for integration
    int termination_type;///termination type - whether radius or de_input time
    Agent_data inclf;///final orbit inclination
    Agent_data raanf;///final orbit right ascension of ascending node
    Agent_data raani;///initial orbit right ascension of ascending node
    Agent_data omegai;///initial orbit argument of periapsis
    Agent_data incli;///initial orbit inclination
    Agent_data mu_e;///gravitational parameter of second major body
    int startbody;///starting central body
    std::string ephem_file;///ephemeris filename
    Agent_data a_second;///second body semi-major axis from overall central
    Agent_data julian_start;///starting Julian date
    Agent_data rfinal;///final periapsis size
    int endbody;///final central body
    int printframe;///frame to print in
    int maxstepcount;///maximum number of integration steps
    Agent_data mu_m;///gravitational parameter of third major body
    Agent_data a_third;///third body semi-major axis from overall central
    int parking_transfer;///rendezvous with planet orbit or transfer to parking orbit
    int finaliz_rendez;///perform final orbit corrections to planetary rendezvous
    std::string outputfile;///timestep output file
    int initwrite;///write finalize rendezvous init conditions to file
    Agent_data Pmax;///variable isp max available power
    Agent_data effic;///variable isp nominal efficiency
    Agent_data IspMax;///variable isp max available isp
    Agent_data IspMin;///variable isp min available isp
};
/**< THIS CLASS IS USED FOR CONVENIENT HANDLING OF AGENTS */
/**< CLASS IS TEMPLATED SO THAT EITHER FLOAT, DOUBLE, LONG DOUBLE ETC CAN BE CHOSEN */
/**< FOR THE DATATYPE THAT IS BEING STORED */
template<typename Agent_data>
class Agent     /**< CAN BE REPLACED WITH STRUCT */
{
    sys_pars<Agent_data> params;
public:
    ///DATA
    Agent_data cost_val;///COST OF AGENT
    Agent_data massfrac;
    std::vector<Agent_data> vals;///COMPONENTS OF THE AGENT
    ///METHODS
    Agent();///DEFAULT CONSTRUCTOR
    Agent(const Agent&);///COPY CONSTRUCTOR - C++ MAKES IT'S OWN - NO NEED TO WRITE THIS
    Agent(const std::vector<Agent_data>&);///SOME OTHER CONSTRUCTOR
    Agent_data cost_fn(const int&);///COST FUNCTION DECLARATION
    Agent_data cost_fn_print(const int&);///COST FUNCTION PRINTER DECLARATION
    void cost_init(const int&);///INITIALIZE ALL AGENTS WITH THEIR COSTS
    void params_init(const std::string);///INITIALIZE SYSTEM PARAMETERS
};
/**< ******************************************** */
#endif // _AGENT_DEFN_HPP
