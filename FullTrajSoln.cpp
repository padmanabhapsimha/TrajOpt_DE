/// ///////////////////////////////////////////////////////////////////////////////
/**<
 *  THREE MAIN BODIES, RIGHT NOW EARTH-SUN-MARS
 *  MEANT FOR PARKING ORBIT-PARKING ORBIT TYPE TRANSFERS
 *  FRAME SWITCHING BASED ON SOI
 *  FORCE MODEL CAN INCLUDE HELIOCENTRIC PERTURBATION IF NECESSARY
 *  JPL EPHEMERIS DE430 USED
 *  CODE DEVELOPED BY PADMANABHA PRASANNA SIMHA
 *  FRAMEID 0=SUN
 *  FRAMEID 1=EARTH
 *  FRAMEID 2=MARS
 */
/// ///////////////////////////////////////////////////////////////////////////////
#include"FullTrajSoln.hpp"
#include<fstream>
#include<cstdlib>
#include<thread>
#include<array>
#include<string>
#include<cmath>
#include<iostream>
#include<boost/numeric/odeint.hpp>
/// ///////////////////////////////////////////////////////////////////////////////
namespace fileLocal_FullTrajSoln{
    auto write_cout=[&](const auto &xfr1,const auto &xfr2,const auto &xfr3,const auto &xE,const auto &xM,const double &t,
                        auto &fileout2,const auto &julianstart,const bool &initCondprint,const auto&initwrite)
    {
        auto printer_lambda=[&](auto &fileoutobj){
            fileoutobj.precision(17);
            for(unsigned int i=0;i<xfr1.size();i++){
                fileoutobj<<xfr1.at(i)<<"\t";
            }
            for(unsigned int i=0;i<xfr2.size();i++){
                fileoutobj<<xfr2.at(i)<<"\t";
            }
            for(unsigned int i=0;i<xfr3.size();i++){
                fileoutobj<<xfr3.at(i)<<"\t";
            }
            for(unsigned int i=0;i<xE.size();i++){
                fileoutobj<<xE.at(i)<<"\t";
            }
            for(unsigned int i=0;i<xM.size();i++){
                fileoutobj<<xM.at(i)<<"\t";
            }
            fileoutobj<<t<<std::endl;
        };
        printer_lambda(fileout2);
        if(initCondprint&&initwrite){///for reinitialization purposes
            std::fstream fileout3("InitConds.txt",std::ios::out);
            printer_lambda(fileout3);
            fileout3<<"\t"<<julianstart;
        }
    };
    auto read_initConds=[&](auto &xfr1,auto &xfr2,auto &xfr3,auto &julianstart){
        std::fstream filein("InitConds.txt",std::ios::in);
        if(!filein.good()){
            exit(-1);
        }
        for(unsigned int i=0;i<14;i++){
            filein>>xfr1.at(i);
        }
        for(unsigned int i=0;i<14;i++){
            filein>>xfr2.at(i);
        }
        for(unsigned int i=0;i<14;i++){
            filein>>xfr3.at(i);
        }
        double dummy;
        for(unsigned int i=0;i<2*6;i++){
            filein>>dummy;
        }
        filein>>dummy>>julianstart;
    };
}
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
Agent_datatype _fullSoln_cost(sys_pars<Agent_datatype> &params,const std::vector<Agent_datatype> &vals,const bool printval,
                           Agent_datatype &massfrac)
{
    using namespace fileLocal_FullTrajSoln;
    std::fstream fileobj(params.outputfile,std::ios::out);///file handle auto closed due to RAII
    if(printval){
        fileobj.close();
        fileobj.open(params.outputfile,std::ios::app);
    }
    /// ///////////////////////////////////////////////////////////////////////////////
    /// CONTROLLED STEPPERS
//    auto controlled_switch_stepper=make_controlled(params.integ_tol,params.integ_tol,0.5*params.maxstep*86400.0,switch_stepper());
    auto controlled_stepper=make_controlled(params.integ_tol,params.integ_tol,params.maxstep*86400.0,stepper_type());
    /// ///////////////////////////////////////////////////////////////////////////////
    wrappy_jpl Ephem_obj(params.ephem_file);///ephemeris object
    std::array<double,6> state_e,state_m,state;///of moving frame
    /// PARAMETERS AND VARIABLES
    const double rsol=params.a_second;
    const double mi=params.mi;
    double julianstart=vals.at(8);
    double mu_s=params.mu;
    double mu_p=params.mu_e;///or params.mu_m;
    double g0Isp=params.g0*params.Isp;
    double x,y,z,vx,vy,vz,m,ax,ay,az,accl,lx,ly,lz,lvx,lvy,lvz,lm;
    double r,r3,r4;
    double xp_e,yp_e,zp_e,vxp_e,vyp_e,vzp_e,rp,rp3,rs,rs3,rrp,rrp3,rrp4,rrs,rrs3,rrs4;
    double xp_m,yp_m,zp_m,vxp_m,vyp_m,vzp_m;
    double xs,ys,zs;
    double xxp_e,yyp_e,zzp_e,xxs,yys,zzs;
    double rx,ry,rz;
    double valx,valy,valz,muval;
    double l,k,costsqrt,T;
    /// ///////////////////////////////////////////////////////////////////////////////
    int frameID=params.startbody;///this needs to be looked up by the planet's dynamics
    auto rhs_helio=[&](const auto &xstate,auto &dxdt,const auto &t){
        /// VARIABLES SETUP
        x=xstate.at(0);y=xstate.at(1);z=xstate.at(2);vx=xstate.at(3);vy=xstate.at(4);vz=xstate.at(5);m=xstate.at(6);
        lx=xstate.at(7);ly=xstate.at(8);lz=xstate.at(9);lvx=xstate.at(10);lvy=xstate.at(11);lvz=xstate.at(12);lm=xstate.at(13);
        r=sqrt(sqr(x)+sqr(y)+sqr(z));r3=cube(r);r4=r3*r;
        ///PLANETARY PERTURBATION. DON'T INCLUDE YET
//        Ephem_obj.JPL_ecliptic_m_get(julianstart+t/86400.0,closeplanetID,12,state,0);
//        xp_e=state.at(0);yp_e=state.at(1);zp_e=state.at(2);
//        rp=sqrt(sqr(xp_e)+sqr(yp_e)+sqr(zp_e));rp3=cube(rp);
//        xxp_e=x-xp_e;yyp_e=y-yp_e;zzp_e=z-zp_e;
//        rrp=sqrt(sqr(xxp_e)+sqr(yyp_e)+sqr(zzp_e));rrp3=cube(rrp);rrp4=rrp3*rrp;
        /// MAX THRUST DETERMINATION
        if(params.NEP==1){
            T=params.kfactor;
        }
        else{
            T=params.kfactor*sqr(rsol/r);
        }
        /// CONTROL LAW
        costsqrt=sqrt(sqr(lvx)+sqr(lvy)+sqr(lvz));
        l=(costsqrt/m)-(1.0-lm)/g0Isp;
        k=-(T/m)/costsqrt;
        if(l>=0.0){
            ax=k*lvx;ay=k*lvy;az=k*lvz;
        }
        else{
            ax=0.0;ay=0.0;az=0.0;
        }
        if(m<0.025*mi){
            ax=0.0;ay=0.0;az=0.0;
        }
        accl=sqrt(sqr(ax)+sqr(ay)+sqr(az));
        /// STATE EQUATIONS
        dxdt.at(0)=vx;dxdt.at(1)=vy;dxdt.at(2)=vz;
        dxdt.at(3)=-(mu_s*x/r3)/**-(mu_p*(xxp_e)/rrp3)**/+ax;
        dxdt.at(4)=-(mu_s*y/r3)/**-(mu_p*(yyp_e)/rrp3)**/+ay;
        dxdt.at(5)=-(mu_s*z/r3)/**-(mu_p*(zzp_e)/rrp3)**/+az;
        dxdt.at(6)=-m*accl/g0Isp;
        /// COSTATE EQUATIONS
        rx=x/r;ry=y/r;rz=z/r;
        valx=((mu_s*x/r4)/**+(mu_p*xxp_e/rrp4)**/);
        valy=((mu_s*y/r4)/**+(mu_p*yyp_e/rrp4)**/);
        valz=((mu_s*z/r4)/**+(mu_p*zzp_e/rrp4)**/);
        muval=-(mu_s/r3)/**-(mu_p/rrp3)**/;
        dxdt.at(7)=-lvx*(muval+3.0*rx*valx)-lvy*(3.0*rx*valy)-lvz*(3.0*rx*valz);
        dxdt.at(8)=-lvx*(3.0*ry*valx)-lvy*(muval+3.0*ry*valy)-lvz*(3.0*ry*valz);
        dxdt.at(9)=-lvx*(3.0*rz*valx)-lvy*(3.0*rz*valy)-lvz*(muval+3.0*rz*valz);
        dxdt.at(10)=-lx;dxdt.at(11)=-ly;dxdt.at(12)=-lz;
        dxdt.at(13)=-(1.0-lm)*accl/g0Isp;
    };
    /// ///////////////////////////////////////////////////////////////////////////////
    auto rhs_planeto=[&](const auto &xstate,auto &dxdt,const auto &t){
        /// VARIABLES SETUP
        x=xstate.at(0);y=xstate.at(1);z=xstate.at(2);vx=xstate.at(3);vy=xstate.at(4);vz=xstate.at(5);m=xstate.at(6);
        lx=xstate.at(7);ly=xstate.at(8);lz=xstate.at(9);lvx=xstate.at(10);lvy=xstate.at(11);lvz=xstate.at(12);lm=xstate.at(13);
        r=sqrt(sqr(x)+sqr(y)+sqr(z));r3=cube(r);r4=r3*r;
        ///SUN'S PERTURBATION EFFECT
        int planetval;
        if(frameID==1){
            planetval=3;
        }
        else if(frameID==2){
            planetval=4;
        }
        Ephem_obj.JPL_ecliptic_m_get(julianstart+t/86400.0,11,planetval,state,0);
        xs=state.at(0);ys=state.at(1);zs=state.at(2);
        rs=sqrt(sqr(xs)+sqr(ys)+sqr(zs));rs3=cube(rs);
        xxs=x-xs;yys=y-ys;zzs=z-zs;
        rrs=sqrt(sqr(xxs)+sqr(yys)+sqr(zzs));rrs3=cube(rrs);rrs4=rrs3*rrs;
        /// MAX THRUST DETERMINATION
        if(params.NEP==1){
            T=params.kfactor;
        }
        else{
            if(mu_p==params.mu_e){
                T=params.kfactor*sqr(rsol)/(sqr(x+xp_e)+sqr(y+yp_e)+sqr(z+zp_e));
            }
            else if(mu_p==params.mu_m){
                T=params.kfactor*sqr(rsol)/(sqr(x+xp_m)+sqr(y+yp_m)+sqr(z+zp_m));
            }
        }
        /// CONTROL LAW
        costsqrt=sqrt(sqr(lvx)+sqr(lvy)+sqr(lvz));
        l=(costsqrt/m)-(1.0-lm)/g0Isp;
        k=-(T/m)/costsqrt;
        if(l>=0.0){
            ax=k*lvx;ay=k*lvy;az=k*lvz;
        }
        else{
            ax=0.0;ay=0.0;az=0.0;
        }
        if(m<0.025*mi){
            ax=0.0;ay=0.0;az=0.0;
        }
        accl=sqrt(sqr(ax)+sqr(ay)+sqr(az));
        /// STATE EQUATIONS
        dxdt.at(0)=vx;dxdt.at(1)=vy;dxdt.at(2)=vz;
        dxdt.at(3)=-(mu_p*x/r3)+mu_s*(-(xs/rs3)-(xxs/rrs3))+ax;
        dxdt.at(4)=-(mu_p*y/r3)+mu_s*(-(ys/rs3)-(yys/rrs3))+ay;
        dxdt.at(5)=-(mu_p*z/r3)+mu_s*(-(zs/rs3)-(zzs/rrs3))+az;
        dxdt.at(6)=-m*accl/g0Isp;
        /// COSTATE EQUATIONS
        rx=x/r;ry=y/r;rz=z/r;
        valx=((mu_p*x/r4)/**+(mu_s*xxs/rrs4)**/);
        valy=((mu_p*y/r4)/**+(mu_s*yys/rrs4)**/);
        valz=((mu_p*z/r4)/**+(mu_s*zzs/rrs4)**/);
        muval=-(mu_p/r3)/**-(mu_s/rrs3)**/;
        dxdt.at(7)=-lvx*(muval+3.0*rx*valx)-lvy*(3.0*rx*valy)-lvz*(3.0*rx*valz);
        dxdt.at(8)=-lvx*(3.0*ry*valx)-lvy*(muval+3.0*ry*valy)-lvz*(3.0*ry*valz);
        dxdt.at(9)=-lvx*(3.0*rz*valx)-lvy*(3.0*rz*valy)-lvz*(muval+3.0*rz*valz);
        dxdt.at(10)=-lx;dxdt.at(11)=-ly;dxdt.at(12)=-lz;
        dxdt.at(13)=-(1.0-lm)*accl/g0Isp;
    };
    /// ///////////////////////////////////////////////////////////////////////////////
    std::array<double,14> xstate,xfr1,xfr2,xfr3;///xstate-spacecraft local state,xfr1-helio,xfr2-earth,xfr3-areo
    /// ///////////////////////////////////////////////////////////////////////////////
    /// INITIALIZE VALUES HERE
    double x0,y0,z0,vx0,vy0,vz0,rp0,H0,aval,pval;
    rp0=params.rpfac*params.rinit;
    aval=rp0/(1.0-params.ecc);
    if(frameID==0){
        H0=sqrt(mu_s*aval*(1.0-sqr(params.ecc)));
        pval=sqr(H0)/mu_s;
    }
    else if(frameID==1){
        mu_p=params.mu_e;
        H0=sqrt(mu_p*aval*(1.0-sqr(params.ecc)));
        pval=sqr(H0)/mu_p;
    }
    else if(frameID==2){
        mu_p=params.mu_m;
        H0=sqrt(mu_p*aval*(1.0-sqr(params.ecc)));
        pval=sqr(H0)/mu_p;
    }
    else{
        std::cout<<std::endl<<"Not supposed to happen"<<std::endl;
    }
    double nui=params.nu*pi/180.0;
    double rval=pval/(1.0+params.ecc*cos(nui));
    x0=rval*(cos(params.raani)*cos(params.omegai+nui)-sin(params.raani)*sin(params.omegai+nui)*cos(params.incli));
    y0=rval*(sin(params.raani)*cos(params.omegai+nui)+cos(params.raani)*sin(params.omegai+nui)*cos(params.incli));
    z0=rval*(sin(params.incli)*sin(params.omegai+nui));
    vx0=(x0*H0*params.ecc*sin(nui)/(rval*pval))-(H0/rval)*(cos(params.raani)*sin(params.omegai+nui)+sin(params.raani)*cos(params.omegai+nui)*cos(params.incli));
    vy0=(y0*H0*params.ecc*sin(nui)/(rval*pval))-(H0/rval)*(sin(params.raani)*sin(params.omegai+nui)-cos(params.raani)*cos(params.omegai+nui)*cos(params.incli));
    vz0=(z0*H0*params.ecc*sin(nui)/(rval*pval))+(H0/rval)*(sin(params.incli)*cos(params.omegai+nui));
    xstate.at(0)=(x0);///X COORDINATE
    xstate.at(1)=(y0);///Y COORDINATE
    xstate.at(2)=(z0);///Z COORDINATE
    xstate.at(3)=(vx0);///X VELOCITY
    xstate.at(4)=(vy0);///Y VELOCITY
    xstate.at(5)=(vz0);///Z VELOCITY
    xstate.at(6)=(mi);///INIT MASS
    for(unsigned int i=0;i<7;i++){
        xstate.at(i+7)=(vals.at(i));///INIT COSTATES
    }
    /// ///////////////////////////////////////////////////////////////////////////////
    double endtime=vals.at(7)*86400.0;
    const double r_soi_sec=params.a_second*pow((params.mu_e/mu_s),0.4);
    const double r_soi_thr=params.a_third*pow((params.mu_m/mu_s),0.4);
    double r2p_e,r2p_m,r2p_e_min,r2p_m_min;
    double t=0.0,dt=1e-3;
    /// ///////////////////////////////////////////////////////////////////////////////
    Ephem_obj.JPL_ecliptic_m_get(julianstart,3,12,state_e,1);
    xp_e=state_e.at(0);yp_e=state_e.at(1);zp_e=state_e.at(2);vxp_e=state_e.at(3);vyp_e=state_e.at(4);vzp_e=state_e.at(5);
    Ephem_obj.JPL_ecliptic_m_get(julianstart,4,12,state_m,1);
    xp_m=state_m.at(0);yp_m=state_m.at(1);zp_m=state_m.at(2);vxp_m=state_m.at(3);vyp_m=state_m.at(4);vzp_m=state_m.at(5);
    /// ///////////////////////////////////////////////////////////////////////////////
    /// IF SPECIFIED INITIAL CONDITIONS ARE AVAILABLE
    if(params.finaliz_rendez){
        read_initConds(xfr1,xfr2,xfr3,julianstart);
        r2p_e=sqrt(sqr(xfr2.at(0))+sqr(xfr2.at(1))+sqr(xfr2.at(2)));
        r2p_m=sqrt(sqr(xfr3.at(0))+sqr(xfr3.at(1))+sqr(xfr3.at(2)));
        if(r2p_m<r_soi_thr){
            frameID=2;
            xstate=xfr3;
        }
        else if(r2p_m<r_soi_sec){
            frameID=1;
            xstate=xfr2;
        }
        else{
            frameID=0;
            xstate=xfr1;
        }
        Ephem_obj.JPL_ecliptic_m_get(julianstart,3,12,state_e,1);
        xp_e=state_e.at(0);yp_e=state_e.at(1);zp_e=state_e.at(2);vxp_e=state_e.at(3);vyp_e=state_e.at(4);vzp_e=state_e.at(5);
        Ephem_obj.JPL_ecliptic_m_get(julianstart,4,12,state_m,1);
        xp_m=state_m.at(0);yp_m=state_m.at(1);zp_m=state_m.at(2);vxp_m=state_m.at(3);vyp_m=state_m.at(4);vzp_m=state_m.at(5);
        for(unsigned int i=0;i<7;i++){
            xstate.at(i+7)=(vals.at(i));///INIT COSTATES
        }
    }
    /// ///////////////////////////////////////////////////////////////////////////////
    if(frameID==0){
        xfr1=xstate;
        xfr2=xfr1;
        xfr2.at(0)-=xp_e;xfr2.at(1)-=yp_e;xfr2.at(2)-=zp_e;xfr2.at(3)-=vxp_e;xfr2.at(4)-=vyp_e;xfr2.at(5)-=vzp_e;
        xfr3=xfr1;
        xfr3.at(0)-=xp_m;xfr3.at(1)-=yp_m;xfr3.at(2)-=zp_m;xfr3.at(3)-=vxp_m;xfr3.at(4)-=vyp_m;xfr3.at(5)-=vzp_m;
    }
    else if(frameID==1){
        xfr2=xstate;
        xfr1=xfr2;
        xfr1.at(0)+=xp_e;xfr1.at(1)+=yp_e;xfr1.at(2)+=zp_e;xfr1.at(3)+=vxp_e;xfr1.at(4)+=vyp_e;xfr1.at(5)+=vzp_e;
        xfr3=xfr1;
        xfr3.at(0)-=xp_m;xfr3.at(1)-=yp_m;xfr3.at(2)-=zp_m;xfr3.at(3)-=vxp_m;xfr3.at(4)-=vyp_m;xfr3.at(5)-=vzp_m;
    }
    else if(frameID==2){
        xfr3=xstate;
        xfr1=xfr3;
        xfr1.at(0)+=xp_m;xfr1.at(1)+=yp_m;xfr1.at(2)+=zp_m;xfr1.at(3)+=vxp_m;xfr1.at(4)+=vyp_m;xfr1.at(5)+=vzp_m;
        xfr2=xfr1;
        xfr2.at(0)-=xp_e;xfr2.at(1)-=yp_e;xfr2.at(2)-=zp_e;xfr2.at(3)-=vxp_e;xfr2.at(4)-=vyp_e;xfr2.at(5)-=vzp_e;
    }
    if(printval){
        write_cout(xfr1,xfr2,xfr3,state_e,state_m,t,fileobj,julianstart+t/86400.0,false,params.initwrite);
    }
    /// ///////////////////////////////////////////////////////////////////////////////
    /// INTEGRATION ROUTINE
    int counter=0;
    ///PLANET CENTRIC DISTANCES
    r2p_e=sqrt(sqr(xfr2.at(0))+sqr(xfr2.at(1))+sqr(xfr2.at(2)));
    r2p_m=sqrt(sqr(xfr3.at(0))+sqr(xfr3.at(1))+sqr(xfr3.at(2)));
    r2p_e_min=r2p_e;r2p_m_min=r2p_m;
    while(t+dt<endtime){
        counter++;
        if(counter>params.maxstepcount){
            break;
        }
        ///DYNAMICS INTEGRATION
        if(r2p_e>=r_soi_sec&&r2p_m>=r_soi_thr){
            ///HELIOCENTRIC PHASE --- largest body
            ///UNIFIED FRAME SWITCHING ROUTINE
            xstate=xfr1;
            frameID=0;
            controlled_stepper.try_step(rhs_helio,xstate,t,dt);
            xfr1=xstate;
            Ephem_obj.JPL_ecliptic_m_get(julianstart+t/86400.0,3,12,state_e,1);
            Ephem_obj.JPL_ecliptic_m_get(julianstart+t/86400.0,4,12,state_m,1);
            xp_e=state_e.at(0);yp_e=state_e.at(1);zp_e=state_e.at(2);vxp_e=state_e.at(3);vyp_e=state_e.at(4);vzp_e=state_e.at(5);
            xp_m=state_m.at(0);yp_m=state_m.at(1);zp_m=state_m.at(2);vxp_m=state_m.at(3);vyp_m=state_m.at(4);vzp_m=state_m.at(5);
            xfr2=xfr1;
            xfr2.at(0)-=xp_e;xfr2.at(1)-=yp_e;xfr2.at(2)-=zp_e;xfr2.at(3)-=vxp_e;xfr2.at(4)-=vyp_e;xfr2.at(5)-=vzp_e;
            xfr3=xfr1;
            xfr3.at(0)-=xp_m;xfr3.at(1)-=yp_m;xfr3.at(2)-=zp_m;xfr3.at(3)-=vxp_m;xfr3.at(4)-=vyp_m;xfr3.at(5)-=vzp_m;
        }
        else if(r2p_e<r_soi_sec&&r2p_m>=r_soi_thr){
            ///GEOCENTRIC PHASE --- second largest
            mu_p=params.mu_e;
            ///UNIFIED FRAME SWITCHING ROUTINE
            xstate=xfr2;
            frameID=1;
            controlled_stepper.try_step(rhs_planeto,xstate,t,dt);
            xfr2=xstate;
            Ephem_obj.JPL_ecliptic_m_get(julianstart+t/86400.0,3,12,state_e,1);
            Ephem_obj.JPL_ecliptic_m_get(julianstart+t/86400.0,4,12,state_m,1);
            xp_e=state_e.at(0);yp_e=state_e.at(1);zp_e=state_e.at(2);vxp_e=state_e.at(3);vyp_e=state_e.at(4);vzp_e=state_e.at(5);
            xp_m=state_m.at(0);yp_m=state_m.at(1);zp_m=state_m.at(2);vxp_m=state_m.at(3);vyp_m=state_m.at(4);vzp_m=state_m.at(5);
            xfr1=xfr2;
            xfr1.at(0)+=xp_e;xfr1.at(1)+=yp_e;xfr1.at(2)+=zp_e;xfr1.at(3)+=vxp_e;xfr1.at(4)+=vyp_e;xfr1.at(5)+=vzp_e;
            xfr3=xfr1;
            xfr3.at(0)-=xp_m;xfr3.at(1)-=yp_m;xfr3.at(2)-=zp_m;xfr3.at(3)-=vxp_m;xfr3.at(4)-=vyp_m;xfr3.at(5)-=vzp_m;
        }
        else if(r2p_m<r_soi_thr){
            ///AREOCENTRIC PHASE --- smallest body
            mu_p=params.mu_m;
            ///UNIFIED FRAME SWITCHING ROUTINE
            xstate=xfr3;
            frameID=2;
            controlled_stepper.try_step(rhs_planeto,xstate,t,dt);
            xfr3=xstate;
            Ephem_obj.JPL_ecliptic_m_get(julianstart+t/86400.0,3,12,state_e,1);
            Ephem_obj.JPL_ecliptic_m_get(julianstart+t/86400.0,4,12,state_m,1);
            xp_e=state_e.at(0);yp_e=state_e.at(1);zp_e=state_e.at(2);vxp_e=state_e.at(3);vyp_e=state_e.at(4);vzp_e=state_e.at(5);
            xp_m=state_m.at(0);yp_m=state_m.at(1);zp_m=state_m.at(2);vxp_m=state_m.at(3);vyp_m=state_m.at(4);vzp_m=state_m.at(5);
            xfr1=xfr3;
            xfr1.at(0)+=xp_m;xfr1.at(1)+=yp_m;xfr1.at(2)+=zp_m;xfr1.at(3)+=vxp_m;xfr1.at(4)+=vyp_m;xfr1.at(5)+=vzp_m;
            xfr2=xfr1;
            xfr2.at(0)-=xp_e;xfr2.at(1)-=yp_e;xfr2.at(2)-=zp_e;xfr2.at(3)-=vxp_e;xfr2.at(4)-=vyp_e;xfr2.at(5)-=vzp_e;
        }
        else{
            ///SHOULD NOT BE POSSIBLE
            ///BUT I HAVE SEEN IT HAPPEN WITH R2P VALUE GOING TO NAN
            ///FOR SOME WEIRD INPUT CASE THAT I DON'T UNDERSTAND SQUAT ABOUT
            std::cout<<std::endl<<"This should not happen during integration "<<r2p_e<<"\t"<<r2p_m<<std::endl;
        }
        if(printval){///PRINT TO FILE IF SPECIFIED
            write_cout(xfr1,xfr2,xfr3,state_e,state_m,t,fileobj,julianstart+t/86400.0,false,params.initwrite);
        }
        ///PLANET CENTRIC DISTANCES
        r2p_e=sqrt(sqr(xfr2.at(0))+sqr(xfr2.at(1))+sqr(xfr2.at(2)));
        r2p_m=sqrt(sqr(xfr3.at(0))+sqr(xfr3.at(1))+sqr(xfr3.at(2)));
        r2p_e_min=(r2p_e<r2p_e_min)?r2p_e:r2p_e_min;
        r2p_m_min=(r2p_m<r2p_m_min)?r2p_m:r2p_m_min;
        dt=(dt>86400.0*params.maxstep)?86400.0*params.maxstep:dt;
        dt=((r2p_e<=1.25*r_soi_sec||r2p_m<7.5*r_soi_thr)&&dt>8640.0*params.maxstep)?8640.0*params.maxstep:dt;
    }
    dt=endtime-t;
    double dtlast=dt;
    if(frameID==0){
        xstate=xfr1;
        controlled_stepper.try_step(rhs_helio,xstate,t,dt);
        xfr1=xstate;
        Ephem_obj.JPL_ecliptic_m_get(julianstart+endtime/86400.0,3,12,state_e,1);
        Ephem_obj.JPL_ecliptic_m_get(julianstart+endtime/86400.0,4,12,state_m,1);
        xp_e=state_e.at(0);yp_e=state_e.at(1);zp_e=state_e.at(2);vxp_e=state_e.at(3);vyp_e=state_e.at(4);vzp_e=state_e.at(5);
        xp_m=state_m.at(0);yp_m=state_m.at(1);zp_m=state_m.at(2);vxp_m=state_m.at(3);vyp_m=state_m.at(4);vzp_m=state_m.at(5);
        xfr2=xfr1;
        xfr2.at(0)-=xp_e;xfr2.at(1)-=yp_e;xfr2.at(2)-=zp_e;xfr2.at(3)-=vxp_e;xfr2.at(4)-=vyp_e;xfr2.at(5)-=vzp_e;
        xfr3=xfr1;
        xfr3.at(0)-=xp_m;xfr3.at(1)-=yp_m;xfr3.at(2)-=zp_m;xfr3.at(3)-=vxp_m;xfr3.at(4)-=vyp_m;xfr3.at(5)-=vzp_m;
    }
    else if(frameID==1){
        xstate=xfr2;
        mu_p=params.mu_e;
        controlled_stepper.try_step(rhs_planeto,xstate,t,dt);
        xfr2=xstate;
        Ephem_obj.JPL_ecliptic_m_get(julianstart+endtime/86400.0,3,12,state_e,1);
        Ephem_obj.JPL_ecliptic_m_get(julianstart+endtime/86400.0,4,12,state_m,1);
        xp_e=state_e.at(0);yp_e=state_e.at(1);zp_e=state_e.at(2);vxp_e=state_e.at(3);vyp_e=state_e.at(4);vzp_e=state_e.at(5);
        xp_m=state_m.at(0);yp_m=state_m.at(1);zp_m=state_m.at(2);vxp_m=state_m.at(3);vyp_m=state_m.at(4);vzp_m=state_m.at(5);
        xfr1=xfr2;
        xfr1.at(0)+=xp_e;xfr1.at(1)+=yp_e;xfr1.at(2)+=zp_e;xfr1.at(3)+=vxp_e;xfr1.at(4)+=vyp_e;xfr1.at(5)+=vzp_e;
        xfr3=xfr1;
        xfr3.at(0)-=xp_m;xfr3.at(1)-=yp_m;xfr3.at(2)-=zp_m;xfr3.at(3)-=vxp_m;xfr3.at(4)-=vyp_m;xfr3.at(5)-=vzp_m;
    }
    else if(frameID==2){
        xstate=xfr3;
        mu_p=params.mu_m;
        controlled_stepper.try_step(rhs_planeto,xstate,t,dt);
        xfr3=xstate;
        Ephem_obj.JPL_ecliptic_m_get(julianstart+endtime/86400.0,3,12,state_e,1);
        Ephem_obj.JPL_ecliptic_m_get(julianstart+endtime/86400.0,4,12,state_m,1);
        xp_e=state_e.at(0);yp_e=state_e.at(1);zp_e=state_e.at(2);vxp_e=state_e.at(3);vyp_e=state_e.at(4);vzp_e=state_e.at(5);
        xp_m=state_m.at(0);yp_m=state_m.at(1);zp_m=state_m.at(2);vxp_m=state_m.at(3);vyp_m=state_m.at(4);vzp_m=state_m.at(5);
        xfr1=xfr3;
        xfr1.at(0)+=xp_m;xfr1.at(1)+=yp_m;xfr1.at(2)+=zp_m;xfr1.at(3)+=vxp_m;xfr1.at(4)+=vyp_m;xfr1.at(5)+=vzp_m;
        xfr2=xfr1;
        xfr2.at(0)-=xp_e;xfr2.at(1)-=yp_e;xfr2.at(2)-=zp_e;xfr2.at(3)-=vxp_e;xfr2.at(4)-=vyp_e;xfr2.at(5)-=vzp_e;
    }
    r2p_e=sqrt(sqr(xfr2.at(0))+sqr(xfr2.at(1))+sqr(xfr2.at(2)));
    r2p_m=sqrt(sqr(xfr3.at(0))+sqr(xfr3.at(1))+sqr(xfr3.at(2)));
    r2p_e_min=(r2p_e<r2p_e_min)?r2p_e:r2p_e_min;
    r2p_m_min=(r2p_m<r2p_m_min)?r2p_m:r2p_m_min;
    if(printval){
        write_cout(xfr1,xfr2,xfr3,state_e,state_m,endtime,fileobj,julianstart+t/86400.0,true,params.initwrite);
    }
    /// ///////////////////////////////////////////////////////////////////////////////
    /// COST FUNCTION EVALUTATION
    double mu;
    if(params.endbody==0){
        mu=mu_s;
        xstate=xfr1;
    }
    else if(params.endbody==1){
        xstate=xfr2;
        mu=params.mu_e;
    }
    else if(params.endbody==2){
        xstate=xfr3;
        mu=params.mu_m;
    }
    else{
        std::cout<<std::endl<<"Not supposed to happen"<<std::endl;
    }
    double rpfin=params.rfinal*params.rf_fac;
    double vpfin=sqrt(mu*(1.0+params.eccf)/rpfin);
    double rfinal=(rpfin/(1.0-params.eccf));
    double xf=xstate.at(0),yf=xstate.at(1),zf=xstate.at(2);
    double vxf=xstate.at(3),vyf=xstate.at(4),vzf=xstate.at(5);
    double rfc=sqrt(sqr(xf)+sqr(yf)+sqr(zf));
    double vfc=sqrt(sqr(vxf)+sqr(vyf)+sqr(vzf));
    double rpf=params.rinit*params.rf_fac;
    double a0=rpfin/(1.0-params.eccf);
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
    double cost_val=(sqr(1.0-af/a0));
    cost_val+=(sqr(Hxi/Hf-Hx/Hfc));
    cost_val+=(sqr(Hyi/Hf-Hy/Hfc));
    cost_val+=(sqr(Hzi/Hf-Hz/Hfc));
    cost_val+=(sqr(ex-exi));
    cost_val+=(sqr(ey-eyi));
    cost_val+=(sqr(ez-ezi));
    cost_val+=sqr(incl-params.inclf)/pi;
    massfrac=1.0-xstate.at(6)/mi;
    double r2p=sqrt(sqr(xstate.at(0))+sqr(xstate.at(1))+sqr(xstate.at(2)));
    double v2p=sqrt(sqr(xstate.at(3))+sqr(xstate.at(4))+sqr(xstate.at(5)));
    double energy=0.5*sqr(v2p)-mu/r2p;
    if(params.endbody==1){
        Ephem_obj.JPL_ecliptic_m_get(julianstart+endtime/86400.0,3,12,state_e,1);
        xp_m=state_e.at(0);yp_m=state_e.at(1);zp_m=state_e.at(2);vxp_m=state_e.at(3);vyp_m=state_e.at(4);vzp_m=state_e.at(5);
    }
    else if(params.endbody==2){
        Ephem_obj.JPL_ecliptic_m_get(julianstart+endtime/86400.0,4,12,state_m,1);
        xp_m=state_m.at(0);yp_m=state_m.at(1);zp_m=state_m.at(2);vxp_m=state_m.at(3);vyp_m=state_m.at(4);vzp_m=state_m.at(5);
    }
//    r2p/=sqrt(sqr(xp_m)+sqr(yp_m)+sqr(zp_m));
    double cost_val2=abs(sqr(xfr1.at(0)-xp_m)+sqr(xfr1.at(1)-yp_m)+sqr(xfr1.at(2)-zp_m))/(sqr(xp_m)+sqr(yp_m)+sqr(zp_m));
    cost_val2+=(sqr(xfr1.at(3)-vxp_m)+sqr(xfr1.at(4)-vyp_m)+sqr(xfr1.at(5)-vzp_m))/(sqr(vxp_m)+sqr(vyp_m)+sqr(vzp_m));
    double retval;
    if(frameID!=params.endbody){
        retval=(1.0+cost_val2)*(1e+8);
    }
    else{
        retval=cost_val;
    }
//    if(frameID==params.startbody&&params.startbody!=params.endbody&&params.startbody!=0){
//        retval=1e8*(1.0+abs(energy));
//    }
    if(params.parking_transfer==0){
        retval=cost_val2;
    }
    else if(params.parking_transfer==2){
        if(params.startbody==1){
            retval=-(0.5*(sqr(xfr2.at(3))+sqr(xfr2.at(4))+sqr(xfr2.at(5)))-params.mu_e/sqrt(sqr(xfr2.at(0))+sqr(xfr2.at(1))+sqr(xfr2.at(2))));
        }
        else if(params.startbody==2){
            retval=-(0.5*(sqr(xfr3.at(3))+sqr(xfr3.at(4))+sqr(xfr3.at(5)))-params.mu_m/sqrt(sqr(xfr3.at(0))+sqr(xfr3.at(1))+sqr(xfr3.at(2))));
        }
    }
    if(abs(dtlast)>params.maxstep*86400.0){
        retval+=(sqr(dtlast));
    }
    return retval;
}
/// ///////////////////////////////////////////////////////////////////////////////
