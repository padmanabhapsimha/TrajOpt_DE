/// ///////////////////////////////////////////////////////////////////////////////
#include"3bodyfunction.hpp"
#include<fstream>
#include<thread>
#include<array>
#include<string>
#include<cmath>
#include<iostream>
#include<boost/numeric/odeint.hpp>
/// ///////////////////////////////////////////////////////////////////////////////
std::fstream fileout2("Output.txt",std::ios::out);
auto write_cout=[&](const auto &xfr1,const auto &xfr2,const double &t)
{
    fileout2.precision(17);
    for(unsigned int i=0;i<xfr1.size();i++){
        fileout2<<xfr1.at(i)<<"\t";
    }
    for(unsigned int i=0;i<xfr2.size();i++){
        fileout2<<xfr2.at(i)<<"\t";
    }
    fileout2<<t<<std::endl;
};
/// ///////////////////////////////////////////////////////////////////////////////
auto sqr=[](const auto &x){return x*x;};
auto cube=[](const auto &x){return x*x*x;};
/// ///////////////////////////////////////////////////////////////////////////////
typedef Agent_datatype state_type;
///RUNGE KUTTA SETUP
using namespace boost::numeric::odeint;
typedef runge_kutta_fehlberg78<std::vector<state_type>> switch_stepper;
typedef runge_kutta_fehlberg78<std::vector<state_type>> stepper_type;
/// ///////////////////////////////////////////////////////////////////////////////
const double pi=3.1415926535897932384626433832795028;
/// ///////////////////////////////////////////////////////////////////////////////
Agent_datatype _3body_cost(sys_pars<Agent_datatype> &params,const std::vector<Agent_datatype> &vals,const bool printval,
                           Agent_datatype &massfrac)
{
    /// ///////////////////////////////////////////////////////////////////////////////
    /// CONTROLLED STEPPERS
    auto controlled_switch_stepper=make_controlled(params.integ_tol,params.integ_tol,params.maxstep*86400.0,switch_stepper());
    auto controlled_stepper=make_controlled(params.integ_tol,params.integ_tol,params.maxstep*86400.0,stepper_type());
    /// ///////////////////////////////////////////////////////////////////////////////
    wrappy_jpl Ephem_obj(params.ephem_file);///ephemeris object
    std::array<double,6> state;///of moving frame
    /// PARAMETERS AND VARIABLES
    const double rsol=(params.startbody==0)?params.rinit:params.rfinal;
    const double mi=params.mi;
    const double julianstart=vals.at(8);
    double mu_s=params.mu;
    double mu_p=params.mu_e;
    double g0Isp=params.g0*params.Isp;
    double x,y,z,vx,vy,vz,m,ax,ay,az,accl,lx,ly,lz,lvx,lvy,lvz,lm;
    double r,r3,r4;
    double xp,yp,zp,vxp,vyp,vzp,rp,rp3,rs,rs3,rrp,rrp3,rrp4,rrs,rrs3,rrs4;
    double xs,ys,zs;
    double xxp,yyp,zzp,xxs,yys,zzs;
    double rx,ry,rz;
    double valx,valy,valz,muval;
    double l,k,costsqrt,T;
    /// ///////////////////////////////////////////////////////////////////////////////
    auto rhs_helio=[&](const auto &xstate,auto &dxdt,const auto &t){
        /// VARIABLES SETUP
        x=xstate.at(0);y=xstate.at(1);z=xstate.at(2);vx=xstate.at(3);vy=xstate.at(4);vz=xstate.at(5);m=xstate.at(6);
        lx=xstate.at(7);ly=xstate.at(8);lz=xstate.at(9);lvx=xstate.at(10);lvy=xstate.at(11);lvz=xstate.at(12);lm=xstate.at(13);
        r=sqrt(sqr(x)+sqr(y)+sqr(z));r3=cube(r);r4=r3*r;
        Ephem_obj.JPL_ecliptic_m_get(julianstart+t/86400.0,3,12,state,0);
        xp=state.at(0);yp=state.at(1);zp=state.at(2);
        rp=sqrt(sqr(xp)+sqr(yp)+sqr(zp));rp3=cube(rp);
        xxp=x-xp;yyp=y-yp;zzp=z-zp;
        rrp=sqrt(sqr(xxp)+sqr(yyp)+sqr(zzp));rrp3=cube(rrp);rrp4=rrp3*rrp;
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
        dxdt.at(3)=-(mu_s*x/r3)-(mu_p*(xxp)/rrp3)+ax;
        dxdt.at(4)=-(mu_s*y/r3)-(mu_p*(yyp)/rrp3)+ay;
        dxdt.at(5)=-(mu_s*z/r3)-(mu_p*(zzp)/rrp3)+az;
        dxdt.at(6)=-m*accl/g0Isp;
        /// COSTATE EQUATIONS
        rx=x/r;ry=y/r;rz=z/r;
        valx=((mu_s*x/r4)+(mu_p*xxp/rrp4));
        valy=((mu_s*y/r4)+(mu_p*yyp/rrp4));
        valz=((mu_s*z/r4)+(mu_p*zzp/rrp4));
        muval=-(mu_s/r3)-(mu_p/rrp3);
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
//        Ephem_obj.JPL_ecliptic_m_get(julianstart+t/86400.0,11,3,state,0);
//        xs=state.at(0);ys=state.at(1);zs=state.at(2);
//        rs=sqrt(sqr(xs)+sqr(ys)+sqr(zs));rs3=cube(rs);
//        xxs=x-xs;yys=y-ys;zzs=z-zs;
//        rrs=sqrt(sqr(xxs)+sqr(yys)+sqr(zzs));rrs3=cube(rrs);rrs4=rrs3*rrs;
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
        if(m<0.025*mi||r<6578000.0){
            ax=0.0;ay=0.0;az=0.0;
        }
        accl=sqrt(sqr(ax)+sqr(ay)+sqr(az));
        /// STATE EQUATIONS
        dxdt.at(0)=vx;dxdt.at(1)=vy;dxdt.at(2)=vz;
        dxdt.at(3)=-(mu_p*x/r3)/**+mu_s*(-(xs/rs3)-(xxs/rrs3))**/+ax;
        dxdt.at(4)=-(mu_p*y/r3)/**+mu_s*(-(ys/rs3)-(yys/rrs3))**/+ay;
        dxdt.at(5)=-(mu_p*z/r3)/**+mu_s*(-(zs/rs3)-(zzs/rrs3))**/+az;
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
    std::array<double,14> xstate,xfr1,xfr2;
    /// ///////////////////////////////////////////////////////////////////////////////
    /// INITIALIZE VALUES HERE
    int frameID=params.startbody;
    double x0,y0,z0,vx0,vy0,vz0,rp0,H0,aval,pval;
    rp0=params.rpfac*params.rinit;
    aval=rp0/(1.0-params.ecc);
    if(frameID==0){
        H0=sqrt(mu_s*aval*(1.0-sqr(params.ecc)));
        pval=sqr(H0)/mu_s;
    }
    else if(frameID==1){
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
    double t=0.0,dt=1e-3;
    /// ///////////////////////////////////////////////////////////////////////////////
    Ephem_obj.JPL_ecliptic_m_get(julianstart,3,12,state,1);
    xp=state.at(0);yp=state.at(1);zp=state.at(2);vxp=state.at(3);vyp=state.at(4);vzp=state.at(5);
    if(frameID==0){
        xfr1=xstate;
        xfr2=xfr1;
        xfr2.at(0)-=xp;xfr2.at(1)-=yp;xfr2.at(2)-=zp;xfr2.at(3)-=vxp;xfr2.at(4)-=vyp;xfr2.at(5)-=vzp;
    }
    else if(frameID==1){
        xfr2=xstate;
        xfr1=xfr2;
        xfr1.at(0)+=xp;xfr1.at(1)+=yp;xfr1.at(2)+=zp;xfr1.at(3)+=vxp;xfr1.at(4)+=vyp;xfr1.at(5)+=vzp;
    }
    if(printval){
        write_cout(xfr1,xfr2,t);
    }
    /// ///////////////////////////////////////////////////////////////////////////////
    /// INTEGRATION ROUTINE
    const double r_soi_sec=params.a_second*pow((mu_p/mu_s),0.4);
    const double switchpoint=r_soi_sec;
    const double deltaswitch=0.325*r_soi_sec;
    double r2p;
    int counter=0;
    double radfinal=params.rf_fac*params.rfinal;
    while(t+dt<endtime){
        r=sqrt(sqr(xfr1.at(0))+sqr(xfr1.at(1))+sqr(xfr1.at(2)));
        counter++;
        if(counter>params.maxstepcount){
            break;
        }
        Ephem_obj.JPL_ecliptic_m_get(julianstart+t/86400.0,3,12,state,1);
        xp=state.at(0);yp=state.at(1);zp=state.at(2);vxp=state.at(3);vyp=state.at(4);vzp=state.at(5);
        r2p=sqrt(sqr(xfr2.at(0))+sqr(xfr2.at(1))+sqr(xfr2.at(2)));
        if(r2p>=switchpoint){
            if(frameID==1){
                xstate.at(0)+=xp;xstate.at(1)+=yp;xstate.at(2)+=zp;xstate.at(3)+=vxp;xstate.at(4)+=vyp;xstate.at(5)+=vzp;
                frameID=0;
            }
            if(abs(r2p-switchpoint)>deltaswitch){
                controlled_stepper.try_step(rhs_helio,xstate,t,dt);
            }
            else{
                controlled_switch_stepper.try_step(rhs_helio,xstate,t,dt);
            }
            xfr1=xstate;xfr2=xfr1;
            Ephem_obj.JPL_ecliptic_m_get(julianstart+t/86400.0,3,12,state,1);
            xp=state.at(0);yp=state.at(1);zp=state.at(2);vxp=state.at(3);vyp=state.at(4);vzp=state.at(5);
            xfr2.at(0)-=xp;xfr2.at(1)-=yp;xfr2.at(2)-=zp;xfr2.at(3)-=vxp;xfr2.at(4)-=vyp;xfr2.at(5)-=vzp;
            frameID=0;
        }
        else if(r2p<switchpoint){
            if(frameID==0){
                xstate.at(0)-=xp;xstate.at(1)-=yp;xstate.at(2)-=zp;xstate.at(3)-=vxp;xstate.at(4)-=vyp;xstate.at(5)-=vzp;
                frameID=1;
            }
            if(abs(r2p-switchpoint)>deltaswitch){
                controlled_stepper.try_step(rhs_planeto,xstate,t,dt);
            }
            else{
                controlled_switch_stepper.try_step(rhs_planeto,xstate,t,dt);
            }
            xfr2=xstate;xfr1=xfr2;
            Ephem_obj.JPL_ecliptic_m_get(julianstart+t/86400.0,3,12,state,1);
            xp=state.at(0);yp=state.at(1);zp=state.at(2);vxp=state.at(3);vyp=state.at(4);vzp=state.at(5);
            xfr1.at(0)+=xp;xfr1.at(1)+=yp;xfr1.at(2)+=zp;xfr1.at(3)+=vxp;xfr1.at(4)+=vyp;xfr1.at(5)+=vzp;
            frameID=1;
        }
        if(printval){
            write_cout(xfr1,xfr2,t);
        }
        dt=(dt>86400.0*params.maxstep)?86400.0*params.maxstep:dt;
    }
    dt=endtime-t;
    if(frameID==0){
        controlled_stepper.try_step(rhs_helio,xstate,t,dt);
        xfr1=xstate;xfr2=xfr1;
        Ephem_obj.JPL_ecliptic_m_get(julianstart+endtime/86400.0,3,12,state,1);
        xp=state.at(0);yp=state.at(1);zp=state.at(2);vxp=state.at(3);vyp=state.at(4);vzp=state.at(5);
        xfr2.at(0)-=xp;xfr2.at(1)-=yp;xfr2.at(2)-=zp;xfr2.at(3)-=vxp;xfr2.at(4)-=vyp;xfr2.at(5)-=vzp;
    }
    else if(frameID==1){
        controlled_stepper.try_step(rhs_planeto,xstate,t,dt);
        xfr2=xstate;xfr1=xfr2;
        Ephem_obj.JPL_ecliptic_m_get(julianstart+endtime/86400.0,3,12,state,1);
        xp=state.at(0);yp=state.at(1);zp=state.at(2);vxp=state.at(3);vyp=state.at(4);vzp=state.at(5);
        xfr1.at(0)+=xp;xfr1.at(1)+=yp;xfr1.at(2)+=zp;xfr1.at(3)+=vxp;xfr1.at(4)+=vyp;xfr1.at(5)+=vzp;
    }
    if(printval){
        write_cout(xfr1,xfr2,t);
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
        mu=mu_p;
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
    double cost_val=(sqr(1.0-af/a0));
    cost_val+=(sqr(Hxi/Hf-Hx/Hfc));
    cost_val+=(sqr(Hyi/Hf-Hy/Hfc));
    cost_val+=(sqr(Hzi/Hf-Hz/Hfc));
    cost_val+=(sqr(ex-exi));
    cost_val+=(sqr(ey-eyi));
    cost_val+=(sqr(ez-ezi));
    massfrac=1.0-xstate.at(6)/mi;
    cost_val+=(sqr(t-endtime));
//    cost_val+=abs(frameID-params.endbody);
    return cost_val;
}
/// ///////////////////////////////////////////////////////////////////////////////
