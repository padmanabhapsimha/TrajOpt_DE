#ifndef _DATATYPE_HPP
#define _DATATYPE_HPP
#include<boost/multiprecision/float128.hpp>
#include<boost/multiprecision/cpp_dec_float.hpp>
/**< ******************************************** */
///VARIABLE TYPE FOR AGENT
///ONLY CHANGE THE TYPE SPECIFIER
///KEEP ALIAS INTACT OR THINGS WILL BLOW UP
///USE WITH FLOAT, DOUBLE OR LONG DOUBLE
typedef double Agent_datatype;
///using namespace boost::multiprecision;
///typedef float128 Agent_datatype;///QUAD PRECISION
/**< ******************************************** */
///TEMPLATE CLASS EXPLICIT  INSTANTIATION
///DO NOT MESS WITH THIS
template class sys_pars<Agent_datatype>;
template class Agent<Agent_datatype>;
/**< ******************************************** */
#endif // _DATATYPE_HPP
