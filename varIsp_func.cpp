///variable isp fuel optimal formulation
///right now single central body
#include"varIsp_func.hpp"
#include<fstream>
#include<cstdlib>
#include<thread>
#include<array>
#include<string>
#include<cmath>
#include<iostream>
#include<boost/numeric/odeint.hpp>
/// ///////////////////////////////////////////////////////////////////////////////
    std::fstream fileobjvar("Output.txt",std::ios::out);///file handle auto closed due to RAII
    auto write_cout=[&](const std::vector<Agent_datatype> &x, const Agent_datatype t)
    {
        fileobjvar.precision(17);
        for(unsigned int i=0;i<x.size();i++){
            fileobjvar<<x.at(i)<<"\t";
        }
        fileobjvar<<t<<std::endl;
    };
/// ///////////////////////////////////////////////////////////////////////////////
auto sqr=[](const auto &x){return x*x;};
auto cube=[](const auto &x){return x*x*x;};
/// ///////////////////////////////////////////////////////////////////////////////
typedef Agent_datatype state_type;
///RUNGE KUTTA SETUP
using namespace boost::numeric::odeint;
//typedef runge_kutta_fehlberg78<std::vector<state_type>> switch_stepper;
typedef runge_kutta_fehlberg78<std::vector<state_type>> stepper_type;
/// ///////////////////////////////////////////////////////////////////////////////
const double pi=3.1415926535897932384626433832795028;
/// ///////////////////////////////////////////////////////////////////////////////
Agent_datatype _varIsp_cost(sys_pars<Agent_datatype> &params,const std::vector<Agent_datatype> &vals,const bool printval,
                           Agent_datatype &massfrac)
{
    using namespace boost::numeric::odeint;
    /// ///////////////////////////////////////////////////////////////////////////////
    auto controlled_stepper=make_controlled(params.integ_tol,params.integ_tol,params.maxstep*86400.0,stepper_type());
    /// ///////////////////////////////////////////////////////////////////////////////
    const double pi=acos(-1.0);
    const double mu=params.mu;
    double g0Isp=params.g0*params.Isp;
    const double mi=params.mi;
    double T=params.kfactor;
    double Pmax;
    const double ri=params.rinit;

    double r,m,r2,r3,r4,msq;
    double l1,k,Isp,kmax,kmin,l3;
    double costsqrt,rx,ry,rz;
    double ax,ay,az,accl;
    ///RHS FOR OF ODE SYSTEM
    auto rhs=[&](const auto x,auto &dxdt, const state_type t){
        r=sqrt(sqr(x.at(0))+sqr(x.at(1))+sqr(x.at(2)));
        r2=sqr(r);r3=cube(r);r4=sqr(r2);
        m=x.at(6);msq=sqr(m);
        costsqrt=sqrt(sqr(x.at(10))+sqr(x.at(11))+sqr(x.at(12)));
        ///MAX AVAILABLE POWER DETERMINATION
        if(params.NEP==1){
            ///NEP
            Pmax=params.Pmax;
        }
        else{
            ///SEP INVERSE SQR LAW
            Pmax=params.Pmax*sqr(ri/r);
        }
        ///CONTROL LAW
        kmax=-((2.0*Pmax*params.effic)/(params.g0*params.IspMax))/(m*costsqrt);
        kmin=-((2.0*Pmax*params.effic)/(params.g0*params.IspMin))/(m*costsqrt);
        l3=-kmax*m*costsqrt*((2.0*(1.0-x.at(13)/(params.g0*params.IspMax)))-(costsqrt/(m)))/params.IspMax;
        if(l3<0){
            Isp=params.IspMin;
            k=-((2.0*Pmax*params.effic)/(params.g0*Isp))/(m*costsqrt);
            l1=((costsqrt/(m))-(1.0-x.at(13))/(params.g0*Isp))/Isp;
            if(l1>=0){
                ax=kmin*x.at(10);ay=kmin*x.at(11);az=kmin*x.at(12);
                g0Isp=params.g0*params.IspMin;
            }
            else{
                ax=0.0;ay=0.0;az=0.0;g0Isp=params.g0*params.IspMin;
            }
        }
        else{
            Isp=params.IspMax;
            k=-((2.0*Pmax*params.effic)/(params.g0*Isp))/(m*costsqrt);
            l1=((costsqrt/(m))-(1.0-x.at(13))/(params.g0*Isp))/Isp;
            if(l1>=0){
                ax=kmax*x.at(10);ay=kmax*x.at(11);az=kmax*x.at(12);
                g0Isp=params.g0*params.IspMax;
            }
            else{
                ax=0.0;ay=0.0;az=0.0;g0Isp=params.g0*params.IspMax;
            }
        }
        accl=sqrt(sqr(ax)+sqr(ay)+sqr(az));
        ///MINIMUM FINAL MASS
        if(m<0.025*mi){
            accl=0.0;ax=0.0;ay=0.0;az=0.0;Isp=0.5*(params.IspMax+params.IspMin);
        }
        ///STATE EQUATIONS
        dxdt.at(0)=x.at(3);
        dxdt.at(1)=x.at(4);
        dxdt.at(2)=x.at(5);
        dxdt.at(3)=(-mu*x.at(0)/r3)+(ax);
        dxdt.at(4)=(-mu*x.at(1)/r3)+(ay);
        dxdt.at(5)=(-mu*x.at(2)/r3)+(az);
        dxdt.at(6)=(-m*accl/g0Isp);
        ///COSTATE EQUATIONS
        rx=x.at(0)/r;ry=x.at(1)/r;rz=x.at(2)/r;
dxdt.at(7)=-x.at(10)*((-mu/r3)+(3.0*mu*x.at(0)*rx/r4))-x.at(11)*(3.0*mu*x.at(1)*rx/r4)-x.at(12)*(3.0*mu*x.at(2)*rx/r4);
dxdt.at(8)=-x.at(10)*(3.0*mu*x.at(0)*ry/r4)-x.at(11)*((-mu/r3)+(3.0*mu*x.at(1)*ry/r4))-x.at(12)*(3.0*mu*x.at(2)*ry/r4);
dxdt.at(9)=-x.at(10)*(3.0*mu*x.at(0)*rz/r4)-x.at(11)*(3.0*mu*x.at(1)*rz/r4)-x.at(12)*((-mu/r3)+(3.0*mu*x.at(2)*rz/r4));
        dxdt.at(10)=-x.at(7);
        dxdt.at(11)=-x.at(8);
        dxdt.at(12)=-x.at(9);
        dxdt.at(13)=-(1.0-x.at(13))*accl/g0Isp;
    };
    ///INITIAL CONDITIONS
    double x0,y0,z0,vx0,vy0,vz0,rp0,H0,aval;
    rp0=params.rpfac*params.rinit;
    aval=rp0/(1.0-params.ecc);
    H0=sqrt(mu*aval*(1.0-sqr(params.ecc)));
    double pval=sqr(H0)/mu;
    double nui=params.nu*3.1415926535897932384626433832795028/180.0;
    double rval=pval/(1.0+params.ecc*cos(nui));
    x0=rval*(cos(params.raani)*cos(params.omegai+nui)-sin(params.raani)*sin(params.omegai+nui)*cos(params.incli));
    y0=rval*(sin(params.raani)*cos(params.omegai+nui)+cos(params.raani)*sin(params.omegai+nui)*cos(params.incli));
    z0=rval*(sin(params.incli)*sin(params.omegai+nui));
    vx0=(x0*H0*params.ecc*sin(nui)/(rval*pval))-(H0/rval)*(cos(params.raani)*sin(params.omegai+nui)+sin(params.raani)*cos(params.omegai+nui)*cos(params.incli));
    vy0=(y0*H0*params.ecc*sin(nui)/(rval*pval))-(H0/rval)*(sin(params.raani)*sin(params.omegai+nui)-cos(params.raani)*cos(params.omegai+nui)*cos(params.incli));
    vz0=(z0*H0*params.ecc*sin(nui)/(rval*pval))+(H0/rval)*(sin(params.incli)*cos(params.omegai+nui));
    std::vector<double> x;
    x.clear();
    for(unsigned int i=0;i<14;i++){
        x.push_back(0.0);
    }
    x.at(0)=(x0);///X COORDINATE
    x.at(1)=(y0);///Y COORDINATE
    x.at(2)=(z0);///Z COORDINATE
    x.at(3)=(vx0);///X VELOCITY
    x.at(4)=(vy0);///Y VELOCITY
    x.at(5)=(vz0);///Z VELOCITY
    x.at(6)=(mi);///INIT MASS
    for(unsigned int i=0;i<7;i++){
        x.at(i+7)=(vals.at(i));
    }
    /// ///////////////////////////////////////////////////////////////////////////////
    double rp=params.rinit*params.rf_fac;
    double vp=sqrt(mu*(1.0+params.eccf)/rp);
    double rfinal=(rp/(1.0-params.eccf));
    double endtime=86400.0*vals.at(7);
    if(printval){
        integrate_adaptive(controlled_stepper,rhs,x,0.0,endtime,0.001,write_cout);
    }
    else{
        integrate_adaptive(controlled_stepper,rhs,x,0.0,endtime,0.001);
    }
    /// ///////////////////////////////////////////////////////////////////////////////
    ///COST FUNCTION EVALUATION
    double xf=x.at(0),yf=x.at(1),zf=x.at(2);
    double vxf=x.at(3),vyf=x.at(4),vzf=x.at(5);
    double rfc=sqrt(sqr(xf)+sqr(yf)+sqr(zf));
    double vfc=sqrt(sqr(vxf)+sqr(vyf)+sqr(vzf));
    double rpf=params.rinit*params.rf_fac;
    double a0=rp/(1.0-params.eccf);
    double af=1.0/((2.0/rfc)-(sqr(vfc)/mu));
    double Hf=sqrt(mu*rpf*(1.0+params.eccf));
    double xi,yi,zi,vxi,vyi,vzi,Hxi,Hyi,Hzi,exi,eyi,ezi;
    double Hx,Hy,Hz,ex,ey,ez;
    xi=rpf*(cos(params.raanf)*cos(params.omegaf)-sin(params.raanf)*sin(params.omegaf)*cos(params.inclf));
    yi=rpf*(sin(params.raanf)*cos(params.omegaf)+cos(params.raanf)*sin(params.omegaf)*cos(params.inclf));
    zi=rpf*(sin(params.inclf)*sin(params.omegaf));
    vxi=-(Hf/rpf)*(cos(params.raanf)*sin(params.omegaf)+sin(params.raanf)*cos(params.omegaf)*cos(params.inclf));
    vyi=-(Hf/rpf)*(sin(params.raanf)*sin(params.omegaf)-cos(params.raanf)*cos(params.omegaf)*cos(params.inclf));
    vzi=(Hf/rpf)*(sin(params.inclf)*cos(params.omegaf));
    Hxi=(yi*vzi-zi*vyi);Hyi=(zi*vxi-xi*vzi);Hzi=(xi*vyi-yi*vxi);
    exi=((vyi*Hzi-vzi*Hyi)/mu)-(xi/rpf);
    eyi=((vzi*Hxi-vxi*Hzi)/mu)-(yi/rpf);
    ezi=((vxi*Hyi-vyi*Hxi)/mu)-(zi/rpf);
    Hx=(yf*vzf-zf*vyf);Hy=(zf*vxf-xf*vzf);Hz=(xf*vyf-yf*vxf);
    double Hfc=sqrt(sqr(Hx)+sqr(Hy)+sqr(Hz));
    Hf=sqrt(sqr(Hxi)+sqr(Hyi)+sqr(Hzi));
    ex=((vyf*Hz-vzf*Hy)/mu)-(xf/rfc);
    ey=((vzf*Hx-vxf*Hz)/mu)-(yf/rfc);
    ez=((vxf*Hy-vyf*Hx)/mu)-(zf/rfc);
    double incl=acos(Hz/Hfc);
    double cost_val=sqr(1.0-af/a0);
    cost_val+=sqr(Hxi/Hf-Hx/Hfc);
    cost_val+=sqr(Hyi/Hf-Hy/Hfc);
    cost_val+=sqr(Hzi/Hf-Hz/Hfc);
    cost_val+=sqr(incl-params.inclf)/pi;
    cost_val+=sqr(ex-exi);
    cost_val+=sqr(ey-eyi);
    cost_val+=sqr(ez-ezi);
    massfrac=(mi-x.at(6))/mi;
//    if(massfrac<0.001){
//        cost_val+=25.0;
//    }
    return sqrt(cost_val);
}
/// ///////////////////////////////////////////////////////////////////////////////

