#include"agent_defn.hpp"
#include"dataTYPE.hpp"
#include"3bodyfunction.hpp"
#include"FullTrajSoln.hpp"
#include"varIsp_func.hpp"
#include<cmath>
#include<vector>
#include<boost/numeric/odeint.hpp>
/**< ******************************************** */
/**< GLOBAL LAMBDA EXPRESSIONS */
auto sqr=[](const auto &x){return x*x;};
auto cube=[](const auto &x){return x*x*x;};
/**< BIG FUNCTION FOR DIFFERENT COST FUNCTION DEFINITIONS */
template<typename Agent_data>
Agent_data Agent<Agent_data>::cost_fn(const int &selector)
{
    cost_val=0;
///    for(double i=0;i<1e4;i++){double j=0;j++;}///MIMICKS EXPENSIVE COST FUNCTION EVALUATION
    switch(selector)
    {
    case 1:{
            /**< SPHERE */
            for(unsigned int i=0;i<vals.size();i++){
                cost_val+=sqr(vals.at(i));
            }
            break;
        }
    case 2:{
            /**< ACKERLEY */
            const Agent_data pi=acos(-1.0L);
            Agent_data x1=vals.at(0);
            Agent_data x2=vals.at(1);
            cost_val=20.0+exp(1.0)-20.0*exp(-0.2*sqrt(0.5*(sqr(x1)+sqr(x2))))-exp(0.5*(cos(2*pi*x1)+cos(2*pi*x2)));///ACKERLEY
            break;
        }
    case 3:{
            /**< ROSENBROCK */
            for(unsigned int i=0;i<vals.size()-1;i++){
                cost_val+=(100*sqr(vals.at(i+1)-sqr(vals.at(i)))+sqr(vals.at(i)-1));
            }
            break;
        }
    case 4:{
            /**< EASOM */
            const Agent_data pi=acos(-1.0L);
            cost_val=-cos(vals.at(0))*cos(vals.at(1))*exp(-(sqr(vals.at(0)-pi)+sqr(vals.at(1)-pi)));
            break;
        }
    case 5:{
            /**< SCHAFFER FUNCTION N. 2 */
            cost_val=0.50+(sqr(sin(sqr(vals.at(0))-sqr(vals.at(1))))-0.50)/sqr(1.0+0.001*(sqr(vals.at(0))+sqr(vals.at(1))));
            break;
        }
    case 6:{
            /**< STYBLINSKI-TANG */
            for(unsigned int i=0;i<vals.size();i++){
                cost_val+=(sqr(sqr(vals.at(i)))-16.0*sqr(vals.at(i))+5.0*vals.at(i));
            }
            cost_val/=2.0;
            break;
        }
    case 7:{
            /**< EGGHOLDER FUNCTION */
            cost_val=-(vals.at(1)+47.0)*sin(sqrt(abs(0.5*vals.at(0)+vals.at(1)+47.0)))-vals.at(0)*sin(sqrt(abs(vals.at(0)-vals.at(1)-47.0)));
            break;
        }
    case 8:{
            /**< HIMMELBLAU */
            cost_val=sqr(sqr(vals.at(0))+vals.at(1)-11.0)+sqr(vals.at(0)+sqr(vals.at(1))-7);
            break;
        }
    case 9:{
            /**< RASTRIGIN */
            const Agent_data pi=acos(-1.0L);
            for(unsigned int i=0;i<vals.size();i++){
                cost_val+=sqr(vals.at(i))-10.0*cos(2*pi*vals.at(i));
            }
            cost_val+=10.0*vals.size();
            break;
        }
    case 10:{
            /**< SCHWEFEL - difficult */
            for(unsigned int i=0;i<vals.size();i++){
                cost_val-=vals.at(i)*sin(sqrt(abs(vals.at(i))));
            }
            cost_val+=418.9829*vals.size();
            ///min is at all indept vars = 420.97
            break;
        }
    case 11:{
            /**< CROSS IN TRAY */
            const Agent_data pi=acos(-1.0L);
            Agent_data x=vals.at(0);Agent_data y=vals.at(1);
            cost_val=-0.0001*pow(abs(sin(x)*sin(y)*exp(abs(100.0-sqrt(sqr(x)+sqr(y))/pi)))+1.0,0.1);
            break;
        }
    case 12:{
            /**< CHICHINADZE */
            const Agent_data pi=acos(-1.0L);
            Agent_data x=vals.at(0);Agent_data y=vals.at(1);
            cost_val=sqr(x)-12.0*x+8.0*sin(2.5*pi*x)+10*cos(0.5*pi*x)+11.0-0.2*sqrt(5)/(exp(0.5*sqr(y-0.5)));
            break;
        }
    case 13:{
            /**< CUBE - similar to 2D Rosenbrock */
            cost_val=100.0*sqr(vals.at(1)-cube(vals.at(0)))+sqr(1-vals.at(0));
            break;
        }
    case 14:{
            /**< DAMAVANDI - very difficult to arrive at global optimum */
            const Agent_data pi=acos(-1.0L);
            Agent_data x=vals.at(0);Agent_data y=vals.at(1);
            cost_val=(1.0-abs(pow((sin(pi*(x-2.0))*sin(pi*(y-2.0)))/(sqr(pi)*(x-2.0)*(y-2.0)),5.0)))*(2.0+sqr(x-7.0)+2*sqr(y-7.0));
            break;
        }
    case 15:{
            /**< HOLDER TABLE FUNCTION */
            const Agent_data pi=acos(-1.0L);
            Agent_data x=vals.at(0);Agent_data y=vals.at(1);
            cost_val=-abs(sin(x)*cos(y)*exp(abs(1.0-sqrt(sqr(x)+sqr(y))/pi)));
            break;
        }
    case 16:{
            /**< OPTIMAL TIME CAR ACCELERATION */
            /**< 3D PROBLEM, 2 COSTATES 1 TIME VARIABLE */
            /**< 2s IS OPTIMUM TIME */
            using namespace boost::numeric::odeint;
            typedef Agent_data state_type;
            typedef runge_kutta_cash_karp54<std::vector<state_type>> stepper_type;
            typedef controlled_runge_kutta<stepper_type> controlled_stepper_type;
            auto rhs=[&](const std::vector<state_type> x, std::vector<state_type> &dxdt, const state_type){
                dxdt.at(0)=x.at(1);
                dxdt.at(1)=-1.0*copysign(1.0,static_cast<state_type>(x.at(2)));
                dxdt.at(2)=-vals.at(0);
            };
            state_type abs_err = 1.0e-16,rel_err = 1.0e-16,a_x = 1.0,a_dxdt = 1.0;
            controlled_stepper_type controlled_stepper(default_error_checker<state_type,
                                                       range_algebra,default_operations>(abs_err,rel_err,a_x,a_dxdt));
            ///SHIFT STATE VECTOR AS AGENT CLASS PRIVATE VARIABLE
            ///THAT'LL PROBABLY INCREASE SPEED
            std::vector<state_type> x;
            x.push_back(0.0);
            x.push_back(0.0);
            x.push_back(vals.at(1));
            integrate_adaptive(controlled_stepper,rhs,x,0.0,vals.at(2),1e-3);

            cost_val=abs(x.at(0)-1.0)+abs(x.at(1));
            ///cost_val-=vals.at(2)*(cost_val<=1e-9)?1.0:0.0;
            break;
        }
    case 17:{
            /**< 2D SOLAR ELECTRIC LOW THRUST MIN TIME */
                using namespace boost::numeric::odeint;
                typedef Agent_datatype state_type;
                typedef runge_kutta_fehlberg78<std::vector<state_type>> stepper_type;
                typedef controlled_runge_kutta<stepper_type> controlled_stepper_type;
                static const double mu=1.32712440018e+20;
                static const double g0=9.806;
                static const double Isp=1500.0;
                static const double g0Isp=g0*Isp;
                static double mi=1000.0;
                double re=149.6e+9;
                double ve=sqrt(mu/re);
                double rm=1.524*re;
                double vm=sqrt(mu/rm);
                static const double k=0.1*sqr(re);
                static const double pi=acos(-1);
                double r;double alph;double T;
                double a1,b1,c1,a2,b2,c2;
                double r4,r3,msq;
                state_type abs_err=1.0e-16,rel_err=1.0e-16,a_x=1.0,a_dxdt=1.0;
                controlled_stepper_type controlled_stepper(default_error_checker<state_type,
                                                           range_algebra,default_operations>(abs_err,rel_err,a_x,a_dxdt));

            auto rhs=[&](const std::vector<state_type> x, std::vector<state_type> &dxdt, const state_type){
                r=sqrt(sqr(x.at(0))+sqr(x.at(1)));///RADIUS
                alph=pi+atan2(x.at(8),x.at(7));///OPTIMAL CONTROL ANGLE
                T=k/sqr(r);
                ///STATE EQUATIONS
                dxdt.at(0)=x.at(2);
                dxdt.at(1)=x.at(3);
                dxdt.at(2)=(-mu*x.at(0)/cube(r))+(T*cos(alph)/x.at(4));
                dxdt.at(3)=(-mu*x.at(1)/cube(r))+(T*sin(alph)/x.at(4));
                dxdt.at(4)=-T/(g0Isp);
                ///COSTATE EQUATIONS
                r4=sqr(sqr(r));r3=cube(r);msq=sqr(x.at(4));
                a1=(-mu/r3)+((3.0*mu*x.at(0)/(r4)*(x.at(0)/r))+(cos(alph)*-2.0*k*x.at(0)/(x.at(4)*r4)));
                a2=(3.0*mu*x.at(0)/(r4)*(x.at(1)/r))+(cos(alph)*-2.0*k*x.at(1)/(x.at(4)*r4));
                b1=(3.0*mu*x.at(1)/(r4)*(x.at(0)/r))+(sin(alph)*-2.0*k*x.at(0)/(x.at(4)*r4));
                b2=(-mu/r3)+((3.0*mu*x.at(1)/(r4)*(x.at(1)/r))+(sin(alph)*-2.0*k*x.at(1)/(x.at(4)*r4)));
                c1=-2.0*k*x.at(0)/(r4*g0Isp);
                c2=-2.0*k*x.at(1)/(r4*g0Isp);
                dxdt.at(5)=-x.at(7)*a1-x.at(8)*b1-x.at(9)*c1;
                dxdt.at(6)=-x.at(7)*a2-x.at(8)*b2-x.at(9)*c2;
                dxdt.at(7)=-x.at(5);
                dxdt.at(8)=-x.at(6);
                dxdt.at(9)=-x.at(7)*(-T*cos(alph)/msq)-x.at(8)*(-T*sin(alph)/msq);
            };
            std::vector<state_type> x;

            x.push_back(re);///X COORDINATE
            x.push_back(0.0);///Y COORDINATE
            x.push_back(0.0);///X VELOCITY
            x.push_back(ve);///Y VELOCITY
            x.push_back(mi);///INIT MASS
            for(unsigned int i=0;i<5;i++){
                x.push_back(vals.at(i));
            }
            double endtime=86400.0*vals.at(5);
            integrate_adaptive(controlled_stepper,rhs,x,0.0,endtime,endtime/1000.0);

            cost_val=abs(-1.0+(sqr(x.at(0))+sqr(x.at(1)))/sqr(rm));///achieve position
            cost_val+=abs(-1.0+(sqr(x.at(2))+sqr(x.at(3)))/sqr(vm));///achieve velocity
            cost_val+=abs(-1.0+sqr(abs(x.at(0)*x.at(3)-x.at(1)*x.at(2))/(rm*vm)));///achieve angular momentum
            break;
        }
    case 18:{
            /**< 2D NUCLEAR ELECTRIC LOW THRUST MIN TIME */
                using namespace boost::numeric::odeint;
                typedef Agent_datatype state_type;
                typedef runge_kutta_fehlberg78<std::vector<state_type>> stepper_type;
//                typedef controlled_runge_kutta<stepper_type> controlled_stepper_type;
                static const double mu=params.mu;//12440018e+20;
                static const double g0=params.g0;
                static const double Isp=params.Isp;
                static const double g0Isp=g0*Isp;
                double re=params.rinit;
//                double ve=sqrt(mu/re);
                double rm=params.rf_fac*params.rinit;
                double vm=sqrt(mu/rm);
                static double mi=params.mi;
                static const double pi=acos(-1.0);
                double r;double alph;double T;
                double a1,b1,c1,a2,b2,c2;
                double r4,r3,msq,cth,sth;
                ///state_type abs_err=1.0e-16,rel_err=1.0e-16,a_x=1.0,a_dxdt=1.0;
                ///controlled_stepper_type controlled_stepper(default_error_checker<state_type,
                ///                                           range_algebra,default_operations>(abs_err,rel_err,a_x,a_dxdt));

                auto controlled_stepper=make_controlled(1.0e-16,1.0e-16,86400.0*params.maxstep,stepper_type());

            auto rhs=[&](const std::vector<state_type> x, std::vector<state_type> &dxdt, const state_type){
                r=sqrt(sqr(x.at(0))+sqr(x.at(1)));///RADIUS
                alph=pi+atan2(x.at(8),x.at(7));///OPTIMAL CONTROL ANGLE
                T=params.kfactor;//k;
                cth=cos(alph);sth=sin(alph);
                ///STATE EQUATIONS
                dxdt.at(0)=x.at(2);
                dxdt.at(1)=x.at(3);
                dxdt.at(2)=(-mu*x.at(0)/cube(r))+(T*cth/x.at(4));
                dxdt.at(3)=(-mu*x.at(1)/cube(r))+(T*sth/x.at(4));
                dxdt.at(4)=-T/(g0Isp);
                ///COSTATE EQUATIONS
                r4=sqr(sqr(r));r3=cube(r);msq=sqr(x.at(4));
                a1=(-mu/r3)+((3.0*mu*x.at(0)/(r4)*(x.at(0)/r)));
                a2=(3.0*mu*x.at(0)/(r4)*(x.at(1)/r));
                b1=(3.0*mu*x.at(1)/(r4)*(x.at(0)/r));
                b2=(-mu/r3)+((3.0*mu*x.at(1)/(r4)*(x.at(1)/r)));
                c1=0.0;
                c2=0.0;
                dxdt.at(5)=-x.at(7)*a1-x.at(8)*b1-x.at(9)*c1;
                dxdt.at(6)=-x.at(7)*a2-x.at(8)*b2-x.at(9)*c2;
                dxdt.at(7)=-x.at(5);
                dxdt.at(8)=-x.at(6);
                dxdt.at(9)=(-x.at(7)*(-T*cth/msq)-x.at(8)*(-T*sth/msq));
            };
            std::vector<state_type> x;
            double a,rad,v,x0,y0,vx0,vy0;
            a=re*params.rpfac/(1.0-params.ecc);
            rad=a*(1.0-sqr(params.ecc))/(1.0+params.ecc*cos(params.nu*pi/180));
            x0=rad*cos(params.nu*pi/180);
            y0=rad*sin(params.nu*pi/180);
            v=sqrt(mu*((2.0/rad)-(1.0/a)));
            double gamma=atan2(params.ecc*sin(params.nu*pi/180.0),1.0+params.ecc*cos(params.nu*pi/180.0));
            vx0=v*cos(0.5*pi-gamma+params.nu*pi/180.0);vy0=v*sin(0.5*pi-gamma+params.nu*pi/180.0);
            x.push_back(x0);///X COORDINATE
            x.push_back(y0);///Y COORDINATE
            x.push_back(vx0);///X VELOCITY
            x.push_back(vy0);///Y VELOCITY
            x.push_back(mi);///INIT MASS
            for(unsigned int i=0;i<5;i++){
                x.push_back(vals.at(i));
            }
            double endtime=86400.0*vals.at(5);
            integrate_adaptive(controlled_stepper,rhs,x,0.0,endtime,0.001);
            cost_val=sqr(-1.0+((sqr(x.at(0))+sqr(x.at(1)))/sqr(rm)));///achieve position
            cost_val+=sqr(-1.0+((sqr(x.at(2))+sqr(x.at(3)))/sqr(vm)));///achieve velocity
            cost_val+=sqr(-1.0+sqr(abs(x.at(0)*x.at(3)-x.at(1)*x.at(2))/(rm*vm)));///achieve angular momentum
            massfrac=(mi-x.at(4))/mi;


            double rp=params.rf_fac*params.rinit;
            double vp=sqrt(mu*(1.0+params.eccf)/rp);
            double nu_c=atan2(x.at(1),x.at(0))-params.omegaf;
            double ra=rp*(1.0+params.eccf)/(1.0+params.eccf*cos(nu_c));
            double va=sqrt(mu*((2.0/ra)-((1.0-params.eccf)/rp)));
            double Ha=rp*vp;
            double rc=sqrt(sqr(x.at(0))+sqr(x.at(1)));
            double vc=sqrt(sqr(x.at(2))+sqr(x.at(3)));
            double Hc=((x.at(0)*x.at(3)-x.at(1)*x.at(2)));
            double ex=(Hc*x.at(3)/mu)-x.at(0)/rc;
            double ey=(-Hc*x.at(2)/mu)-x.at(1)/rc;
            double af=rp/(1.0-params.eccf);
            double ac=1.0/((2.0/rc)-sqr(vc)/mu);
            double eccf=sqrt(sqr(ex)+sqr(ey));

            cost_val=(abs(-1.0+sqr(rc/ra)));///achieve position
            cost_val+=(abs(-1.0+sqr(vc/va)));///achieve velocity
            cost_val+=(sqr(-1.0+(Hc/Ha)));///achieve angular momentum
            cost_val+=sqr(ex-params.eccf*cos(params.omegaf));
            cost_val+=sqr(ey-params.eccf*sin(params.omegaf));


//            cost_val=abs((ra-rc)/1000.0);
//            double velx,vely,rex,rey;
//            rex=x.at(0)/rc;rey=x.at(1)/rc;
//            velx=-1.0*rex;vely=-1.0*rey;
//            cost_val+=abs(x.at(2)-velx)+abs(x.at(3)-vely);


            break;
        }
    case 19:{
            /**< 2D CONTINUOUS THRUST MIN EFFORT */
            using namespace boost::numeric::odeint;
            typedef Agent_datatype state_type;
            typedef runge_kutta_fehlberg78<std::vector<state_type>> stepper_type;
            static const double mu=params.mu;
            static const double g0=params.g0;
            static const double Isp=params.Isp;
            static const double g0Isp=g0*Isp;
            static double mi=params.mi;
            double re=params.rinit;
            double rm=params.rf_fac*params.rinit;
            double vm=sqrt(mu/rm);
            static const double pi=acos(-1.0);
            double r,alph,T;
            double a1,b1,a2,b2,c1,c2;
            double r2,r4,r3,msq,m;
            double costsqrt,p,l1,l2,calp,salp;
            double dTx,dTy;
//            int NEP=params.NEP;
            auto controlled_stepper=make_controlled(1.0e-16,1.0e-16,1.0*86400.0,stepper_type());

            auto rhs=[&](const std::vector<state_type> x, std::vector<state_type> &dxdt, const state_type){
                r=sqrt(sqr(x.at(0))+sqr(x.at(1)));///RADIUS
                r2=sqr(r);r4=(sqr(r2));r3=cube(r);
                m=x.at(4);msq=sqr(x.at(4));
                if(params.NEP==1){
                    ///NEP
                    T=params.kfactor;
                    dTx=0.0;dTy=0.0;
                }
                else{
                    ///SEP INVERSE SQR LAW
                    T=params.kfactor*sqr(re/r);
                    if(r<0.25*re){///SATURATE THRUST
                        T=params.kfactor*sqr(1.0/0.25);
                    }
                    dTx=-2.0*T*x.at(0)/r2;dTy=-2.0*T*x.at(1)/r2;
                }
                ///condition to prevent bad stuff
                if(m<0.1*mi){
                    T=0.0;
                }
                costsqrt=sqrt(sqr(x.at(7))+sqr(x.at(8)));
                alph=pi+atan2(x.at(8),x.at(7));
                calp=cos(alph);salp=sin(alph);
                p=1.0;
                l1=T*((costsqrt/m)+(x.at(9)/g0Isp)-T);
                ///l2=-T*((costsqrt/m)+(x.at(9)/g0Isp));
                l2=-l1-sqr(T);
                if(l2>=0.0){
                    p=0.0;
                }
                else if(l1<0.0&&l2<0.0){
                    p=-l2/sqr(T);
                }
                ///STATE EQUATIONS
                dxdt.at(0)=x.at(2);
                dxdt.at(1)=x.at(3);
                dxdt.at(2)=(-mu*x.at(0)/r3)+(p*T*calp/m);
                dxdt.at(3)=(-mu*x.at(1)/r3)+(p*T*salp/m);
                dxdt.at(4)=-p*T/(g0Isp);
                ///COSTATE EQUATIONS
                a1=(-mu/r3)+((3.0*mu*x.at(0)/(r4)*(x.at(0)/r)))+(p*calp*dTx/m);
                a2=(3.0*mu*x.at(0)/(r4)*(x.at(1)/r))+(p*calp*dTy/m);
                b1=(3.0*mu*x.at(1)/(r4)*(x.at(0)/r))+(p*salp*dTx/m);
                b2=(-mu/r3)+((3.0*mu*x.at(1)/(r4)*(x.at(1)/r)))+(p*salp*dTy/m);
                c1=(-p*dTx/g0Isp);c2=(-p*dTy/g0Isp);
                dxdt.at(5)=-x.at(7)*a1-x.at(8)*b1-x.at(9)*c1;
                dxdt.at(6)=-x.at(7)*a2-x.at(8)*b2-x.at(9)*c1;
                dxdt.at(7)=-x.at(5);
                dxdt.at(8)=-x.at(6);
                dxdt.at(9)=-p*T*costsqrt/msq;
            };
            std::vector<state_type> x;
            double a,rad,v,x0,y0,vx0,vy0;
            a=re*params.rpfac/(1.0-params.ecc);
            rad=a*(1.0-sqr(params.ecc))/(1.0+params.ecc*cos(params.nu*pi/180));
            x0=rad*cos(params.nu*pi/180);
            y0=rad*sin(params.nu*pi/180);
            v=sqrt(mu*((2.0/rad)-(1.0/a)));
            double gamma=atan2(params.ecc*sin(params.nu*pi/180.0),1.0+params.ecc*cos(params.nu*pi/180.0));
            vx0=v*cos(0.5*pi-gamma+params.nu*pi/180.0);vy0=v*sin(0.5*pi-gamma+params.nu*pi/180.0);
            x.push_back(x0);///X COORDINATE
            x.push_back(y0);///Y COORDINATE
            x.push_back(vx0);///X VELOCITY
            x.push_back(vy0);///Y VELOCITY
            x.push_back(mi);///INIT MASS
            for(unsigned int i=0;i<5;i++){
                x.push_back(vals.at(i));
            }
            double endtime=86400.0*vals.at(5);
            integrate_adaptive(controlled_stepper,rhs,x,0.0,endtime,0.1);
            /// /////////////////////////////////////////////////////////////
            double rp=params.rf_fac*params.rinit;
            double vp=sqrt(mu*(1.0+params.eccf)/rp);
            double nu_c=atan2(x.at(1),x.at(0))-params.omegaf*pi/180.0;
            double ra=rp*(1.0+params.eccf)/(1.0+params.eccf*cos(nu_c));
            double va=sqrt(mu*((2.0/ra)-((1.0-params.eccf)/rp)));
            double Ha=rp*vp;
            double rc=sqrt(sqr(x.at(0))+sqr(x.at(1)));
            double vc=sqrt(sqr(x.at(2))+sqr(x.at(3)));
            double Hc=((x.at(0)*x.at(3)-x.at(1)*x.at(2)));
            double ex=(Hc*x.at(3)/mu)-x.at(0)/rc;
            double ey=(-Hc*x.at(2)/mu)-x.at(1)/rc;
            double af=rp/(1.0-params.eccf);
            double ac=1.0/((2.0/rc)-sqr(vc)/mu);
            double eccf=sqrt(sqr(ex)+sqr(ey));

            cost_val=(abs(-1.0+sqr(rc/ra)));///achieve position
            cost_val+=(abs(-1.0+sqr(vc/va)));///achieve velocity
            cost_val+=(sqr(-1.0+(Hc/Ha)));///achieve angular momentum
            cost_val+=sqr(ex-params.eccf*cos(params.omegaf*pi/180.0));
            cost_val+=sqr(ey-params.eccf*sin(params.omegaf*pi/180.0));

            massfrac=(mi-x.at(4))/mi;
            cost_val*=(cost_val<1e-2)?massfrac:1.0;
            break;
        }
    case 20:{
            /**< 2D CONTINUOUS THRUST TANGENTIAL */
                using namespace boost::numeric::odeint;
                typedef Agent_datatype state_type;
                typedef runge_kutta_fehlberg78<std::vector<state_type>> stepper_type;
                typedef controlled_runge_kutta<stepper_type> controlled_stepper_type;
                static const double mu=params.mu;
                static const double g0=params.g0;
                static const double Isp=params.Isp;
                static const double g0Isp=g0*Isp;
                static double mi=params.mi;
                double r,alph,T;
                const double pi=acos(-1.0);
                double m;
                double p,calp,salp;
                state_type abs_err=1.0e-16,rel_err=1.0e-16,a_x=1.0,a_dxdt=1.0;
                controlled_stepper_type controlled_stepper(default_error_checker<state_type,
                                                           range_algebra,default_operations>(abs_err,rel_err,a_x,a_dxdt));

            auto rhs=[&](const std::vector<state_type> x, std::vector<state_type> &dxdt, const state_type){
                r=sqrt(sqr(x.at(0))+sqr(x.at(1)));///RADIUS
                T=params.kfactor*sqr(params.rinit/r);
                m=x.at(4);
                ///condition to prevent bad stuff
                if(m<0.1*mi){
                    T=0.0;
                }
                alph=atan2(x.at(3),x.at(2));
                calp=cos(alph);salp=sin(alph);
                p=1.0;
                ///STATE EQUATIONS
                dxdt.at(0)=x.at(2);
                dxdt.at(1)=x.at(3);
                dxdt.at(2)=(-mu*x.at(0)/cube(r))-(p*T*calp/m);
                dxdt.at(3)=(-mu*x.at(1)/cube(r))-(p*T*salp/m);
                dxdt.at(4)=-p*T/(g0Isp);
            };
            const double re=params.rinit;
            double a,rad,v,x0,y0,vx0,vy0;
            a=re*params.rpfac/(1.0-params.ecc);
            rad=a*(1.0-sqr(params.ecc))/(1.0+params.ecc*cos(params.nu*pi/180));
            x0=rad*cos(params.nu*pi/180);
            y0=rad*sin(params.nu*pi/180);
            v=sqrt(mu*((2.0/rad)-(1.0/a)));
            double gamma=atan2(params.ecc*sin(params.nu*pi/180.0),1.0+params.ecc*cos(params.nu*pi/180.0));
            vx0=v*cos(0.5*pi-gamma+params.nu*pi/180.0);vy0=v*sin(0.5*pi-gamma+params.nu*pi/180.0);
            std::vector<Agent_data> x;
            x.push_back(x0);///X COORDINATE
            x.push_back(y0);///Y COORDINATE
            x.push_back(vx0);///X VELOCITY
            x.push_back(vy0);///Y VELOCITY
            x.push_back(mi);///INIT MASS
            for(unsigned int i=0;i<5;i++){
                x.push_back(vals.at(i));
            }
            double endtime=86400.0*vals.at(5);
            integrate_adaptive(controlled_stepper,rhs,x,0.0,endtime,endtime/1000.0);
            double rm=params.rf_fac*params.rinit;
            cost_val=(abs(-1.0+((sqr(x.at(0))+sqr(x.at(1)))/sqr(rm))));///achieve position
            if(x.at(4)<0.1*mi){
                cost_val+=1e5;
            }
            massfrac=(mi-x.at(4))/mi;
            break;
        }
    case 21:{
            /**< 2D CONTINUOUS THRUST MIN effort */
            /**< ax, ay FORMULATION WITH THRUST LIMIT */
            using namespace boost::numeric::odeint;
            typedef Agent_datatype state_type;
            typedef runge_kutta_fehlberg78<std::vector<state_type>> stepper_type;
//            typedef controlled_runge_kutta<stepper_type> controlled_stepper_type;
            static const double mu=params.mu;
            static const double g0=params.g0;
            static const double Isp=params.Isp;
            static const double g0Isp=g0*Isp;
            static double mi=params.mi;
            double ri=params.rinit;
//            double ve=sqrt(mu/re);
            double rm=params.rf_fac*params.rinit;
            double vm=sqrt(mu/rm);
            static const double pi=acos(-1.0);
            double r,T,accl;
            double r4,r3,m;
            double l,costsqrt;
            double ax,ay;
            double a1,a2,b1,b2;

            auto controlled_stepper=make_controlled(params.integ_tol,params.integ_tol,params.maxstep*86400.0,stepper_type());

            auto rhs=[&](const std::vector<state_type> x, std::vector<state_type> &dxdt, const state_type){
                r=sqrt(sqr(x.at(0))+sqr(x.at(1)));///RADIUS
                r4=sqr(sqr(r));r3=cube(r);
                m=x.at(4);
                if(params.NEP==1){
                    ///NEP
                    T=params.kfactor;
                }
                else{
                    ///SEP INVERSE SQR LAW
                    T=params.kfactor*sqr(ri/r);
                    ///COVERSTONE CARROLL
//                    T=(params.kfactor)*((0.7119713 + 0.4089753/(r/ri) + -0.0783905/(sqr(r/ri)))/(1.0+-0.0201896*(r/ri)+0.0623406*sqr(r/ri)))/(sqr(r/ri));
                }
                costsqrt=sqrt(sqr(x.at(7))+sqr(x.at(8)));
                l=(costsqrt/m)-(x.at(9)/(g0Isp));
                if(l>=0.0){
                    ax=-x.at(7)*T/(costsqrt*m);
                    ay=-x.at(8)*T/(costsqrt*m);
                }
                else{
                    ax=0.0;ay=0.0;
                }
                accl=sqrt(sqr(ax)+sqr(ay));
                if(m<0.025*mi){
                    ax=0.0;ay=0.0;accl=0.0;
                }
                ///STATE EQUATIONS
                dxdt.at(0)=x.at(2);
                dxdt.at(1)=x.at(3);
                dxdt.at(2)=(-mu*x.at(0)/cube(r))+(ax);
                dxdt.at(3)=(-mu*x.at(1)/cube(r))+(ay);
                dxdt.at(4)=-x.at(4)*accl/(g0Isp);
                ///COSTATE EQUATIONS
                a1=(-mu/r3)+((3.0*mu*x.at(0)/(r4)*(x.at(0)/r)));
                a2=(3.0*mu*x.at(0)/(r4)*(x.at(1)/r));
                b1=(3.0*mu*x.at(1)/(r4)*(x.at(0)/r));
                b2=(-mu/r3)+((3.0*mu*x.at(1)/(r4)*(x.at(1)/r)));
                dxdt.at(5)=-x.at(7)*a1-x.at(8)*b1;
                dxdt.at(6)=-x.at(7)*a2-x.at(8)*b2;
                dxdt.at(7)=-x.at(5);
                dxdt.at(8)=-x.at(6);
                dxdt.at(9)=+x.at(9)*accl/g0Isp;
            };
            std::vector<state_type> x;
            double a,rad,v,x0,y0,vx0,vy0;
            a=ri*params.rpfac/(1.0-params.ecc);
            rad=a*(1.0-sqr(params.ecc))/(1.0+params.ecc*cos(params.nu*pi/180));
            x0=rad*cos(params.nu*pi/180);
            y0=rad*sin(params.nu*pi/180);
            v=sqrt(mu*((2.0/rad)-(1.0/a)));
            double gamma=atan2(params.ecc*sin(params.nu*pi/180.0),1.0+params.ecc*cos(params.nu*pi/180.0));
            vx0=v*cos(0.5*pi-gamma+params.nu*pi/180.0);vy0=v*sin(0.5*pi-gamma+params.nu*pi/180.0);
            x.push_back(x0);///X COORDINATE
            x.push_back(y0);///Y COORDINATE
            x.push_back(vx0);///X VELOCITY
            x.push_back(vy0);///Y VELOCITY
            x.push_back(mi);///INIT MASS
            for(unsigned int i=0;i<5;i++){
                x.push_back(vals.at(i));
            }
            double endtime=86400.0*vals.at(5);
            integrate_adaptive(controlled_stepper,rhs,x,0.0,endtime,endtime/1000.0);

            double rp=params.rf_fac*params.rinit;
            double vp=sqrt(mu*(1.0+params.eccf)/rp);
            double nu_c=atan2(x.at(1),x.at(0))-params.omegaf;
            double ra=rp*(1.0+params.eccf)/(1.0+params.eccf*cos(nu_c));
            double va=sqrt(mu*((2.0/ra)-((1.0-params.eccf)/rp)));
            double Ha=rp*vp;
            double rc=sqrt(sqr(x.at(0))+sqr(x.at(1)));
            double vc=sqrt(sqr(x.at(2))+sqr(x.at(3)));
            double Hc=((x.at(0)*x.at(3)-x.at(1)*x.at(2)));
            double ex=(Hc*x.at(3)/mu)-x.at(0)/rc;
            double ey=(-Hc*x.at(2)/mu)-x.at(1)/rc;
            double af=rp/(1.0-params.eccf);
            double ac=1.0/((2.0/rc)-sqr(vc)/mu);
            double eccf=sqrt(sqr(ex)+sqr(ey));

            cost_val=(abs(-1.0+sqr(rc/ra)));///achieve position
            cost_val+=(abs(-1.0+sqr(vc/va)));///achieve velocity
            cost_val+=(sqr(-1.0+(Hc/Ha)));///achieve angular momentum
            cost_val+=sqr(ex-params.eccf*cos(params.omegaf));
            cost_val+=sqr(ey-params.eccf*sin(params.omegaf));

            massfrac=(mi-x.at(4))/mi;
            cost_val*=(cost_val<5e-2)?massfrac:1.0;

///optimal descent formulation
//            cost_val=abs((ra-rc)/1000.0);
//            double velx,vely,rex,rey;
//            rex=x.at(0)/rc;rey=x.at(1)/rc;
//            velx=-25.00*rex;vely=-25.00*rey;
//            cost_val+=abs(x.at(2)-velx)+abs(x.at(3)-vely);

            break;
        }
    case 22:{
            /**< 3D MIN CONTROL EFFORT OPTIMAL CONTROL FOR SPACECRAFT */
            using namespace boost::numeric::odeint;
            typedef Agent_datatype state_type;
            typedef runge_kutta_fehlberg78<std::vector<state_type>> stepper_type;
//            typedef controlled_runge_kutta<stepper_type> controlled_stepper_type;
            /// ///////////////////////////////////////////////////////////////////////////////
            auto controlled_stepper=make_controlled(1.0e-16,1.0e-16,10*86400.0,stepper_type());
            /// ///////////////////////////////////////////////////////////////////////////////
            const double pi=acos(-1.0);
            const double mu=params.mu;
            const double g0Isp=params.g0*params.Isp;
            const double mi=params.mi;
            double T=params.kfactor;
            const double ri=params.rinit;
//            const double rf=ri*params.rf_fac;
//            const double ecc=params.ecc;
//            const double nu=params.nu;

            double r,alpha,beta,m,r2,r3,r4,msq,p,l1,l2;
            double costsqrt,dTx,dTy,dTz,rx,ry,rz,ca,sa,cb,sb;
            double a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3;
            ///RHS FOR OF ODE SYSTEM
            auto rhs=[&](const std::vector<state_type> x, std::vector<state_type> &dxdt, const state_type t){
                r=sqrt(sqr(x.at(0))+sqr(x.at(1))+sqr(x.at(2)));
                r2=sqr(r);r3=cube(r);r4=sqr(r2);
                m=x.at(6);msq=sqr(m);
                costsqrt=sqrt(sqr(x.at(10))+sqr(x.at(11))+sqr(x.at(12)));
                ///MAX AVAILABLE THRUST DETERMINATION
                if(params.NEP==1){
                    ///NEP
                    T=params.kfactor;
                    dTx=0.0;dTy=0.0;dTz=0.0;
                }
                else{
                    ///SEP INVERSE SQR LAW
                    T=params.kfactor*sqr(ri/r);
                    if(r<0.25*ri){///SATURATE THRUST
                        T=params.kfactor*sqr(1.0/0.25);
                    }
                    dTx=-2.0*T*x.at(0)/r2;dTy=-2.0*T*x.at(1)/r2;dTz=-2.0*T*x.at(2)/r2;
                }
                ///CONTROL LAW
                l1=-T*((costsqrt/m)+(x.at(13)/g0Isp));
                l2=-l1-sqr(T);
                alpha=pi+atan2(x.at(11),x.at(10));
                beta=pi+0.0*atan2(x.at(12),-sqrt(sqr(x.at(10))+sqr(x.at(11))));
                p=1.0;
                if(l1>=0){
                    p=0.0;
                }
                if(l2>=0){
                    p=1.0;
                }
                if(l1<0&&l2<0){
                    p=-l1/sqr(T);
                }
                ///MINIMUM FINAL MASS
                if(m<0.10*mi){
                    p=0.0;
                }
                ca=cos(alpha);sa=sin(alpha);cb=cos(beta);sb=sin(beta);
                ///STATE EQUATIONS
                dxdt.at(0)=x.at(3);
                dxdt.at(1)=x.at(4);
                dxdt.at(2)=x.at(5);
                dxdt.at(3)=(-mu*x.at(0)/r3)+(p*T*ca*cb/m);
                dxdt.at(4)=(-mu*x.at(1)/r3)+(p*T*sa*cb/m);
                dxdt.at(5)=(-mu*x.at(2)/r3)+(p*T*sb/m);
                dxdt.at(6)=(-p*T/g0Isp);
                ///COSTATE EQUATIONS
                rx=x.at(0)/r;ry=x.at(1)/r;rz=x.at(2)/r;
                a1=(-mu/r3)+(3.0*mu*x.at(0)*rx/r4)+(p*ca*cb*dTx/m);
                a2=(3.0*mu*x.at(0)*ry/r4)+(p*ca*cb*dTy/m);
                a3=(3.0*mu*x.at(0)*rz/r4)+(p*ca*cb*dTz/m);
                b1=(3.0*mu*x.at(1)*rx/r4)+(p*sa*cb*dTx/m);
                b2=(-mu/r3)+(3.0*mu*x.at(1)*ry/r4)+(p*sa*cb*dTy/m);
                b3=(3.0*mu*x.at(1)*rz/r4)+(p*sa*cb*dTz/m);
                c1=(3.0*mu*x.at(1)*rx/r4)+(p*sb*dTx/m);
                c2=(3.0*mu*x.at(1)*ry/r4)+(p*sb*dTy/m);
                c3=(-mu/r3)+(3.0*mu*x.at(1)*rz/r4)+(p*sb*dTz/m);
                d1=(-p*dTx/g0Isp);d2=(-p*dTy/g0Isp);d3=(-p*dTz/g0Isp);
                dxdt.at(7)=-x.at(10)*a1-x.at(11)*b1-x.at(12)*c1-x.at(13)*d1;
                dxdt.at(8)=-x.at(10)*a2-x.at(11)*b2-x.at(12)*c2-x.at(13)*d2;
                dxdt.at(9)=-x.at(10)*a3-x.at(11)*b3-x.at(12)*c3-x.at(13)*d3;
                dxdt.at(10)=-x.at(7);
                dxdt.at(11)=-x.at(8);
                dxdt.at(12)=-x.at(9);
                dxdt.at(13)=-p*(costsqrt)/msq;
            };
            ///INITIAL CONDITIONS
            const double re=params.rinit;
//            const double ve=sqrt(mu/re);
            const double rm=params.rf_fac*params.rinit;
            const double vm=sqrt(mu/rm);
            double a,rad,v,x0,y0,vx0,vy0;
            a=re*params.rpfac/(1.0-params.ecc);
            rad=a*(1.0-sqr(params.ecc))/(1.0+params.ecc*cos(params.nu*pi/180));
            x0=rad*cos(params.nu*pi/180);
            y0=rad*sin(params.nu*pi/180);
            v=sqrt(mu*((2.0/rad)-(1.0/a)));
            double gamma=atan2(params.ecc*sin(params.nu*pi/180.0),1.0+params.ecc*cos(params.nu*pi/180.0));
            vx0=v*cos(0.5*pi-gamma+params.nu*pi/180.0);vy0=v*sin(0.5*pi-gamma+params.nu*pi/180.0);
            std::vector<Agent_data> x;
            x.push_back(x0);///X COORDINATE
            x.push_back(y0);///Y COORDINATE
            x.push_back(0.0);///Z COORDINATE
            x.push_back(vx0);///X VELOCITY
            x.push_back(vy0);///Y VELOCITY
            x.push_back(0.0);///Z VELOCITY
            x.push_back(mi);///INIT MASS
            for(unsigned int i=0;i<7;i++){
                x.push_back(vals.at(i));
            }
            double endtime=86400.0*vals.at(7);
            /// ///////////////////////////////////////////////////////////////////////////////
            integrate_adaptive(controlled_stepper,rhs,x,0.0,endtime,1.0);
            /// ///////////////////////////////////////////////////////////////////////////////
            cost_val=1.5*(abs(-1.0+((sqr(x.at(0))+sqr(x.at(1))+sqr(x.at(2)))/sqr(rm))));///achieve position
            cost_val+=1.0*(abs(-1.0+((sqr(x.at(3))+sqr(x.at(4))+sqr(x.at(5)))/sqr(vm))));///achieve velocity
            cost_val+=0.25*tanh(abs(x.at(2)/rm)+abs(x.at(5)/vm));///out of plane motion
            cost_val+=1.0*(abs(-1.0+sqr(abs(x.at(0)*x.at(4)-x.at(1)*x.at(3))/(rm*vm))));///achieve angular momentum
            cost_val*=(cost_val<=1e-2)?(0.15*(mi-x.at(6))/mi):1.0;
//            cost_val+=(abs(x.at(1)*x.at(5)-x.at(2)*x.at(4))/(rm*vm));///ang momentum
//            cost_val+=(abs(x.at(2)*x.at(3)-x.at(0)*x.at(5))/(rm*vm));///ang momentum
            break;
        }
    case 23:{
            /**< 2D MIN FUEL OPTIMAL CONTROL FOR SPACECRAFT */
            /**< NO INTERMEDIATE STATE */
            using namespace boost::numeric::odeint;
            typedef Agent_datatype state_type;
            typedef runge_kutta_fehlberg78<std::vector<state_type>> stepper_type;
            /// ///////////////////////////////////////////////////////////////////////////////
            auto controlled_stepper=make_controlled(1.0e-16,1.0e-16,1.0*86400.0,stepper_type());
            /// ///////////////////////////////////////////////////////////////////////////////
            const double pi=acos(-1.0);
            const double mu=params.mu;
            const double g0Isp=params.g0*params.Isp;
            const double mi=params.mi;
            double T=params.kfactor;
            const double ri=params.rinit;

            double r,alpha,m,r2,r3,r4,msq,p,l1;
            double costsqrt,dTx,dTy,rx,ry,ca,sa;
            double a1,a2,b1,b2,c1,c2;
            ///RHS FOR OF ODE SYSTEM
            auto rhs=[&](const std::vector<state_type> x, std::vector<state_type> &dxdt, const state_type t){
                r=sqrt(sqr(x.at(0))+sqr(x.at(1)));
                r2=sqr(r);r3=cube(r);r4=sqr(r2);
                m=x.at(4);msq=sqr(m);
                costsqrt=sqrt(sqr(x.at(7))+sqr(x.at(8)));
                ///MAX AVAILABLE THRUST DETERMINATION
                if(params.NEP==1){
                    ///NEP
                    T=params.kfactor;
                    dTx=0.0;dTy=0.0;
                }
                else{
                    ///SEP INVERSE SQR LAW
                    T=params.kfactor*sqr(ri/r);
                    if(r<0.25*ri){///SATURATE THRUST
                        T=params.kfactor*sqr(1.0/0.25);
                    }
                    dTx=-2.0*T*x.at(0)/r2;dTy=-2.0*T*x.at(1)/r2;
                }
                ///CONTROL LAW
                l1=T*((costsqrt/m)+((x.at(9)-1.0)/g0Isp));
                alpha=pi+atan2(x.at(8),x.at(7));
                if(l1>=0){
                    p=1.0;
                }
                else{
                    p=0.0;
                }
                ///MINIMUM FINAL MASS
                if(m<0.025*mi){
                    p=0.0;
                }
                ca=cos(alpha);sa=sin(alpha);
                ///STATE EQUATIONS
                dxdt.at(0)=x.at(2);
                dxdt.at(1)=x.at(3);
                dxdt.at(2)=(-mu*x.at(0)/r3)+(p*T*ca/m);
                dxdt.at(3)=(-mu*x.at(1)/r3)+(p*T*sa/m);
                dxdt.at(4)=(-p*T/g0Isp);
                ///COSTATE EQUATIONS
                rx=x.at(0)/r;ry=x.at(1)/r;
                a1=(-mu/r3)+(3.0*mu*x.at(0)*rx/r4)+(p*ca*dTx/m);
                a2=(3.0*mu*x.at(0)*ry/r4)+(p*ca*dTy/m);
                b1=(3.0*mu*x.at(1)*rx/r4)+(p*sa*dTx/m);
                b2=(-mu/r3)+(3.0*mu*x.at(1)*ry/r4)+(p*sa*dTy/m);
                c1=(-p*dTx/g0Isp);c2=(-p*dTy/g0Isp);
                dxdt.at(5)=-x.at(7)*a1-x.at(8)*b1-x.at(9)*c1;
                dxdt.at(6)=-x.at(7)*a2-x.at(8)*b2-x.at(9)*c2;
                dxdt.at(7)=-x.at(5);
                dxdt.at(8)=-x.at(6);
                dxdt.at(9)=-p*T*(costsqrt)/msq;
            };
            ///INITIAL CONDITIONS
            const double re=params.rinit;
            double a,rad,v,x0,y0,vx0,vy0;
            a=re*params.rpfac/(1.0-params.ecc);
            rad=a*(1.0-sqr(params.ecc))/(1.0+params.ecc*cos(params.nu*pi/180));
            x0=rad*cos(params.nu*pi/180);
            y0=rad*sin(params.nu*pi/180);
            v=sqrt(mu*((2.0/rad)-(1.0/a)));
            double gamma=atan2(params.ecc*sin(params.nu*pi/180.0),1.0+params.ecc*cos(params.nu*pi/180.0));
            vx0=v*cos(0.5*pi-gamma+params.nu*pi/180.0);vy0=v*sin(0.5*pi-gamma+params.nu*pi/180.0);
            std::vector<Agent_data> x;
            x.push_back(x0);///X COORDINATE
            x.push_back(y0);///Y COORDINATE
            x.push_back(vx0);///X VELOCITY
            x.push_back(vy0);///Y VELOCITY
            x.push_back(mi);///INIT MASS
            for(unsigned int i=0;i<5;i++){
                x.push_back(vals.at(i));
            }
            double endtime=86400.0*vals.at(5);
            /// ///////////////////////////////////////////////////////////////////////////////
            integrate_adaptive(controlled_stepper,rhs,x,0.0,endtime,0.001);
            /// ///////////////////////////////////////////////////////////////////////////////
            double rp=params.rf_fac*params.rinit;
            double vp=sqrt(mu*(1.0+params.eccf)/rp);
            double nu_c=atan2(x.at(1),x.at(0))-params.omegaf*pi/180.0;
            double ra=rp*(1.0+params.eccf)/(1.0+params.eccf*cos(nu_c));
            double va=sqrt(mu*((2.0/ra)-((1.0-params.eccf)/rp)));
            double Ha=rp*vp;
            double rc=sqrt(sqr(x.at(0))+sqr(x.at(1)));
            double vc=sqrt(sqr(x.at(2))+sqr(x.at(3)));
            double Hc=((x.at(0)*x.at(3)-x.at(1)*x.at(2)));
            double ex=(Hc*x.at(3)/mu)-x.at(0)/rc;
            double ey=(-Hc*x.at(2)/mu)-x.at(1)/rc;
            double af=rp/(1.0-params.eccf);
            double ac=1.0/((2.0/rc)-sqr(vc)/mu);
            double eccf=sqrt(sqr(ex)+sqr(ey));

            cost_val=(abs(-1.0+sqr(rc/ra)));///achieve position
            cost_val+=(abs(-1.0+sqr(vc/va)));///achieve velocity
            cost_val+=(sqr(-1.0+(Hc/Ha)));///achieve angular momentum
            cost_val+=sqr(ex-params.eccf*cos(params.omegaf*pi/180.0));
            cost_val+=sqr(ey-params.eccf*sin(params.omegaf*pi/180.0));

            massfrac=(mi-x.at(4))/mi;
            cost_val*=(cost_val<1e-2)?massfrac:1.0;

///opt landing
//            cost_val=abs((ra-rc)/1000.0);
//            double velx,vely,rex,rey;
//            rex=x.at(0)/rc;rey=x.at(1)/rc;
//            velx=-0.01*rex;vely=-0.01*rey;
//            cost_val+=abs(x.at(2)-velx)+abs(x.at(3)-vely);


            break;
        }
    case 24:{
            /**< 2D MIN FUEL OPTIMAL CONTROL FOR SPACECRAFT */
            /**< INTEGRATION TERMINATION BASED ON RADIUS */
            /**< OUTWARD PROPAGATING */
            /**< NO INTERMEDIATE STATE */
            using namespace boost::numeric::odeint;
            typedef Agent_datatype state_type;
            typedef runge_kutta_fehlberg78<std::vector<state_type>> stepper_type;
            /// ///////////////////////////////////////////////////////////////////////////////
            auto controlled_stepper=make_controlled(params.integ_tol,params.integ_tol,params.maxstep*86400.0,stepper_type());
            /// ///////////////////////////////////////////////////////////////////////////////
            const double pi=acos(-1.0);
            const double mu=params.mu;
            const double g0Isp=params.g0*params.Isp;
            const double mi=params.mi;
            double T=params.kfactor;
            const double ri=params.rinit;

            double r,alpha,m,r2,r3,r4,msq,p,l1;
            double costsqrt,dTx,dTy,rx,ry,ca,sa;
            double a1,a2,b1,b2,c1,c2;
            ///RHS FOR OF ODE SYSTEM
            auto rhs=[&](const std::vector<state_type> x, std::vector<state_type> &dxdt, const state_type t){
                r=sqrt(sqr(x.at(0))+sqr(x.at(1)));
                r2=sqr(r);r3=cube(r);r4=sqr(r2);
                m=x.at(4);msq=sqr(m);
                costsqrt=sqrt(sqr(x.at(7))+sqr(x.at(8)));
                ///MAX AVAILABLE THRUST DETERMINATION
                if(params.NEP==1){
                    ///NEP
                    T=params.kfactor;
                    dTx=0.0;dTy=0.0;
                }
                else{
                    ///SEP INVERSE SQR LAW
                    T=params.kfactor*sqr(ri/r);
                    if(r<0.25*ri){///SATURATE THRUST
                        T=params.kfactor*sqr(1.0/0.25);
                    }
                    dTx=-2.0*T*x.at(0)/r2;dTy=-2.0*T*x.at(1)/r2;
                }
                ///CONTROL LAW
                l1=T*((costsqrt/m)+((x.at(9)-1.0)/g0Isp));
                alpha=pi+atan2(x.at(8),x.at(7));
                if(l1>=0){
                    p=1.0;
                }
                else{
                    p=0.0;
                }
                ///MINIMUM FINAL MASS
                if(m<0.025*mi){
                    p=0.0;
                }
                ca=cos(alpha);sa=sin(alpha);
                ///STATE EQUATIONS
                dxdt.at(0)=x.at(2);
                dxdt.at(1)=x.at(3);
                dxdt.at(2)=(-mu*x.at(0)/r3)+(p*T*ca/m);
                dxdt.at(3)=(-mu*x.at(1)/r3)+(p*T*sa/m);
                dxdt.at(4)=(-p*T/g0Isp);
                ///COSTATE EQUATIONS
                rx=x.at(0)/r;ry=x.at(1)/r;
                a1=(-mu/r3)+(3.0*mu*x.at(0)*rx/r4)+(p*ca*dTx/m);
                a2=(3.0*mu*x.at(0)*ry/r4)+(p*ca*dTy/m);
                b1=(3.0*mu*x.at(1)*rx/r4)+(p*sa*dTx/m);
                b2=(-mu/r3)+(3.0*mu*x.at(1)*ry/r4)+(p*sa*dTy/m);
                c1=(-p*dTx/g0Isp);c2=(-p*dTy/g0Isp);
                dxdt.at(5)=-x.at(7)*a1-x.at(8)*b1-x.at(9)*c1;
                dxdt.at(6)=-x.at(7)*a2-x.at(8)*b2-x.at(9)*c2;
                dxdt.at(7)=-x.at(5);
                dxdt.at(8)=-x.at(6);
                dxdt.at(9)=-p*T*(costsqrt)/msq;
            };
            ///INITIAL CONDITIONS
            const double re=params.rinit;
            double a,rad,v,x0,y0,vx0,vy0;
            a=re*params.rpfac/(1.0-params.ecc);
            rad=a*(1.0-sqr(params.ecc))/(1.0+params.ecc*cos(params.nu*pi/180));
            x0=rad*cos(params.nu*pi/180);
            y0=rad*sin(params.nu*pi/180);
            v=sqrt(mu*((2.0/rad)-(1.0/a)));
            double gamma=atan2(params.ecc*sin(params.nu*pi/180.0),1.0+params.ecc*cos(params.nu*pi/180.0));
            vx0=v*cos(0.5*pi-gamma+params.nu*pi/180.0);vy0=v*sin(0.5*pi-gamma+params.nu*pi/180.0);
            std::vector<Agent_data> x;
            x.push_back(x0);///X COORDINATE
            x.push_back(y0);///Y COORDINATE
            x.push_back(vx0);///X VELOCITY
            x.push_back(vy0);///Y VELOCITY
            x.push_back(mi);///INIT MASS
            for(unsigned int i=0;i<5;i++){
                x.push_back(vals.at(i));
            }
            /// ///////////////////////////////////////////////////////////////////////////////
            double rac,rprev;
            double t,dt,tprev,dtprev;
            std::vector<double> xprev;
            t=0.0;dt=1e-3;
            ///r=1.0/((2.0/(sqrt((x.at(0))+(x.at(1)))))-(sqr(x.at(2))+sqr(x.at(3)))/mu);
            rac=sqrt(sqr(x.at(0))+sqr(x.at(1)));
            /// ///////////////////////////////////////
            double rp=params.rf_fac*params.rinit;
            double vp=sqrt(mu*(1.0+params.eccf)/rp);
            double rfinal=(rp/(1.0-params.eccf));
            /// ///////////////////////////////////////
            while(1){
                tprev=t;dtprev=dt;xprev=x;
                rprev=rac;
                controlled_stepper.try_step(rhs,x,t,dt);
//                r=1.0/((2.0/(sqrt((x.at(0))+(x.at(1)))))-(sqr(x.at(2))+sqr(x.at(3)))/mu);
                rac=sqrt(sqr(x.at(0))+sqr(x.at(1)));
                dt=(dt>1.0*86400)?(1.0*86400):dt;
                if(rac>rfinal){
                    t=tprev;
                    x=xprev;
                    rac=rprev;
                    dt=0.5*(dtprev);
                }
                if(abs(rac-rfinal)/abs(rfinal)<params.a_tol||t>86400.0*params.maxtime){
                    break;
                }
            }
            /// ///////////////////////////////////////////////////////////////////////////////

            double nu_c=atan2(x.at(1),x.at(0))-params.omegaf*pi/180.0;
            double ra=rp*(1.0+params.eccf)/(1.0+params.eccf*cos(nu_c));
            double va=sqrt(mu*((2.0/ra)-((1.0-params.eccf)/rp)));
            double Ha=rp*vp;
            double rc=sqrt(sqr(x.at(0))+sqr(x.at(1)));
            double vc=sqrt(sqr(x.at(2))+sqr(x.at(3)));
            double Hc=((x.at(0)*x.at(3)-x.at(1)*x.at(2)));
            double ex=(Hc*x.at(3)/mu)-x.at(0)/rc;
            double ey=(-Hc*x.at(2)/mu)-x.at(1)/rc;

            double ac=1.0/((2.0/rc)-sqr(vc)/mu);
            double eccf=sqrt(sqr(ex)+sqr(ey));

            cost_val=(abs(-1.0+sqr(rc/ra)));///achieve position
            cost_val+=(abs(-1.0+sqr(vc/va)));///achieve velocity
            cost_val+=(sqr(-1.0+(Hc/Ha)));///achieve angular momentum
            cost_val+=sqr(ex-params.eccf*cos(params.omegaf*pi/180.0));
            cost_val+=sqr(ey-params.eccf*sin(params.omegaf*pi/180.0));

            massfrac=(mi-x.at(4))/mi;
            cost_val*=(cost_val<1e-2)?massfrac:1.0;
            break;
        }
    case 25:{
            /**< 3D MIN FUEL OPTIMAL CONTROL FOR SPACECRAFT */
            /**< NO INTERMEDIATE THROTTLE */
            using namespace boost::numeric::odeint;
            typedef Agent_datatype state_type;
            typedef runge_kutta_fehlberg78<std::vector<state_type>> stepper_type;
            /// ///////////////////////////////////////////////////////////////////////////////
            auto controlled_stepper=make_controlled(params.integ_tol,params.integ_tol,params.maxstep*86400.0,stepper_type());
            /// ///////////////////////////////////////////////////////////////////////////////
            const double pi=acos(-1.0);
            const double mu=params.mu;
            double g0Isp=params.g0*params.Isp;
            const double mi=params.mi;
            double T=params.kfactor;
            const double ri=params.rinit;

            double r,m,r2,r3,r4,msq,l,rf;
            double costsqrt,rx,ry,rz;
            double ax,ay,az,accl;
            ///RHS FOR OF ODE SYSTEM
            auto rhs=[&](const auto x,auto &dxdt, const state_type t){
                r=sqrt(sqr(x.at(0))+sqr(x.at(1))+sqr(x.at(2)));
                r2=sqr(r);r3=cube(r);r4=sqr(r2);
                m=x.at(6);msq=sqr(m);
                costsqrt=sqrt(sqr(x.at(10))+sqr(x.at(11))+sqr(x.at(12)));
                ///MAX AVAILABLE THRUST DETERMINATION
                if(params.NEP==1){
                    ///NEP
                    T=params.kfactor;
                }
                else{
                    ///SEP INVERSE SQR LAW
                    T=params.kfactor*sqr(ri/r);
//                    T=(T>0.236)?0.236:T;
//                    T=(T<0.040)?0.040:T;
//                    g0Isp=1900.0+(T-0.017)*10456.0;g0Isp*=9.806;
//                    if(r<0.25*ri){///SATURATE THRUST
//                        T=params.kfactor*sqr(1.0/0.25);
//                    }
                }
                ///CONTROL LAW
                l=(costsqrt/m)-(x.at(13)/g0Isp);
                rf=-(T/m)/costsqrt;
                if(l>=0.0){
                    ax=rf*x.at(10);ay=rf*x.at(11);az=rf*x.at(12);
                }
                else{
                    ax=0.0;ay=0.0;az=0.0;
                }
                accl=sqrt(sqr(ax)+sqr(ay)+sqr(az));
                ///MINIMUM FINAL MASS
                if(m<0.025*mi){
                    accl=0.0;ax=0.0;ay=0.0;az=0.0;
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
                dxdt.at(13)=+x.at(13)*accl/g0Isp;
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
            thread_local std::array<Agent_data,14> x;
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
            if(params.termination_type==2){
                double r,rprev;
                double t,dt,tprev,dtprev;
                thread_local std::array<double,14> xprev;
                t=0.0;dt=1e-3;
                r=sqrt(sqr(x.at(0))+sqr(x.at(1))+sqr(x.at(2)));
                /// ///////////////////////////////////////
            /// ///////////////////////////////////////
                while(1){
                    tprev=t;dtprev=dt;xprev=x;
                    rprev=r;
                    controlled_stepper.try_step(rhs,x,t,dt);
                    r=sqrt(sqr(x.at(0))+sqr(x.at(1))+sqr(x.at(2)));
                    dt=(dt>864000.0)?(864000):dt;
                    if(r>rfinal){
                        t=tprev;
                        x=xprev;
                        r=rprev;
                        dt=0.5*(dtprev);
                    }
                    if(abs(r-rfinal)/abs(rfinal)<params.a_tol||t>86400.0*params.maxtime){
                        break;
                    }
                }
            }
            else{
                double endtime=86400.0*vals.at(7);
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
            cost_val=sqr(1.0-af/a0);
//            cost_val+=sqr(1.0-Hfc/Hf);
            cost_val+=sqr(Hxi/Hf-Hx/Hfc);
            cost_val+=sqr(Hyi/Hf-Hy/Hfc);
            cost_val+=sqr(Hzi/Hf-Hz/Hfc);
            cost_val+=sqr(incl-params.inclf)/pi;
            cost_val+=sqr(ex-exi);
            cost_val+=sqr(ey-eyi);
            cost_val+=sqr(ez-ezi);
            massfrac=(mi-x.at(6))/mi;
//            cost_val*=(cost_val<1e-2)?massfrac:1.0;
            ///UNCOMMENT BELOW LINE FOR PLANETARY ESCAPE
//            cost_val=-(0.5*sqr(vfc)-mu/rfc)/1e8;

///optimal landing
//            cost_val=abs((rpf-rfc)/1000.0);
//            double velx,vely,rex,rey;
//            rex=x.at(0)/rfc;rey=x.at(1)/rfc;
//            velx=-25.00*rex;vely=-25.0*rey;
//            cost_val+=abs(x.at(3)-velx)+abs(x.at(4)-vely);



            ///UNCOMMENT BELOW LINE FOR ACHIEVING LONGITUDE IN GTO-GSO TRANSFER
            ///ASSUMING GSO is in X-Y plane
//            double theta=2.0*pi*vals.at(7);
//            double perigee_offset=180.0*(pi/180.0);
//            double xfin=rpf*cos(theta+perigee_offset);
//            double yfin=rpf*sin(theta+perigee_offset);
//            double zfin=0.0;
//            double vtot=sqrt(params.mu/rpf);
//            double vxfin=vtot*cos(theta+perigee_offset+pi/2.0);
//            double vyfin=vtot*sin(theta+perigee_offset+pi/2.0);
//            double vzfin=0.0;
//            double cost_val2=(sqr(xfin-x.at(0))+sqr(yfin-x.at(1))+sqr(zfin-x.at(2)))/sqr(rpf);
//            cost_val2+=(sqr(vxfin-x.at(3))+sqr(vyfin-x.at(4))+sqr(vzfin-x.at(5)))/sqr(vtot);
//            cost_val+=sqrt(cost_val2);

            break;
        }
    case 26:{
            ///3 body fuel optimal soln
            cost_val=_3body_cost(params,vals,false,massfrac);
            break;
        }
    case 27:{
            ///complete parking orbit to parking orbit transfer
            ///ephemeris model type
            cost_val=_fullSoln_cost(params,vals,false,massfrac);
            break;
        }
    case 28:{
            ///variable specific impulse fuel optimal formulation
            cost_val=_varIsp_cost(params,vals,false,massfrac);
            break;
        }
    default:{
            /**< INSERT CUSTOM FUNCTION CALL */
            cost_val=0.0;
        }
    }
    return cost_val;
}
/**< ******************************************** */
