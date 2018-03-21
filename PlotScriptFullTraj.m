close all;clear;clc
data=importdata('Output.txt');
sys=importdata('System_params.txt','\n');
close all;
%%
mu_s=1.32712440018e+20;
mu_e=398600e+9;
mu_m=4.282837e+13;
g0=sys(4);Isp=sys(3);g0Isp=g0*Isp;
kfac=sys(5);
a_earthsun=149.597870700e+9;
%%
X=data(:,1);Y=data(:,2);Z=data(:,3);VX=data(:,4);VY=data(:,5);VZ=data(:,6);M=data(:,7);
LX=data(:,8);LY=data(:,9);LZ=data(:,10);LVX=data(:,11);LVY=data(:,12);LVZ=data(:,13);LM=data(:,14);
x=data(:,15);y=data(:,16);z=data(:,17);vx=data(:,18);vy=data(:,19);vz=data(:,20);m=data(:,21);
lx=data(:,22);ly=data(:,23);lz=data(:,24);lvx=data(:,25);lvy=data(:,26);lvz=data(:,27);lm=data(:,28);
xm=data(:,29);ym=data(:,30);zm=data(:,31);vxm=data(:,32);vym=data(:,33);vzm=data(:,34);mm=data(:,35);
lxm=data(:,36);lym=data(:,37);lzm=data(:,38);lvxm=data(:,39);lvym=data(:,40);lvzm=data(:,41);lmm=data(:,42);
xE=data(:,43);yE=data(:,44);zE=data(:,45);vxE=data(:,46);vyE=data(:,47);vzE=data(:,48);
xM=data(:,49);yM=data(:,50);zM=data(:,51);vxM=data(:,52);vyM=data(:,53);vzM=data(:,54);
t=data(:,55);
%% geocentric traj
figure(1)
subplot(1,3,1)
plot3(x,y,z,'k.-')
axis equal;
xlim([-max(abs([x;y])) max(abs([x;y]))])
ylim([-max(abs([x;y])) max(abs([x;y]))])
grid on;grid minor;
view([0 0 90]);
xlabel('X axis');ylabel('Y axis');title('Geocentric trajectory');
%% areocentric traj
subplot(1,3,2)
plot3(xm,ym,zm,'k.-')
axis equal;
xlim([-max(abs([xm;ym])) max(abs([xm;ym]))])
ylim([-max(abs([xm;ym])) max(abs([xm;ym]))])
grid on;grid minor;
view([0 0 90]);
xlabel('X axis');ylabel('Y axis');title('Mars-centric trajectory');
%% heliocentric traj
figure(1)
subplot(1,3,3)
plot3(X,Y,Z,'k.-',xE,yE,zE,'b.-',xM,yM,zM,'r.-')
axis equal;
xlim([-max(abs([X;Y;xM;yM;xE;yE])) max(abs([X;Y;xM;yM;xE;yE]))])
ylim([-max(abs([X;Y;xM;yM;xE;yE])) max(abs([X;Y;xM;yM;xE;yE]))])
grid on;grid minor;
view([0 0 90]);
xlabel('X axis');ylabel('Y axis');title('Heliocentric trajectory');
%% distances
figure(2)
R=sqrt(X.^2+Y.^2+Z.^2);
V=sqrt(VX.^2+VY.^2+VZ.^2);
r=sqrt(x.^2+y.^2+z.^2);
v=sqrt(vx.^2+vy.^2+vz.^2);
rm=sqrt(xm.^2+ym.^2+zm.^2);
vm=sqrt(vxm.^2+vym.^2+vzm.^2);
subplot(3,1,1)
plot(t/86400,R/1000,'k-');
grid on;grid minor;
xlabel('Time (days)');ylabel('Radial distance (km)');title('Heliocentric');
subplot(3,1,2)
plot(t/86400,r/1000,'k-');
grid on;grid minor;
xlabel('Time (days)');ylabel('Radial distance (km)');title('Geocentric');
subplot(3,1,3)
plot(t/86400,rm/1000,'k-');
grid on;grid minor;
xlabel('Time (days)');ylabel('Radial distance (km)');title('Mars-centric');
%% velocity
figure(3)
subplot(3,1,1)
plot(t/86400,V/1000,'k-');
grid on;grid minor;
xlabel('Time (days)');ylabel('Velocity (km/s)');title('Heliocentric');
subplot(3,1,2)
plot(t/86400,v/1000,'k-');
grid on;grid minor;
xlabel('Time (days)');ylabel('Velocity (km/s)');title('Geocentric');
subplot(3,1,3)
plot(t/86400,vm/1000,'k-');
grid on;grid minor;
xlabel('Time (days)');ylabel('Velocity (km/s)');title('Mars-centric');
%% eccentricity
A=1./((2./R)-(V.^2/mu_s));
a=1./((2./r)-(v.^2/mu_e));
am=1./((2./rm)-(vm.^2/mu_m));
figure(4)
subplot(3,1,1)
plot(t/86400,A/1000,'k-');
grid on;grid minor;
xlabel('Time (days)');ylabel('Semi-major axis (km)');title('Heliocentric');
subplot(3,1,2)
plot(t/86400,a/1000,'k-');
grid on;grid minor;
xlabel('Time (days)');ylabel('Semi-major axis (km)');title('Geocentric');
subplot(3,1,3)
plot(t/86400,am/1000,'k-');
grid on;grid minor;
xlabel('Time (days)');ylabel('Semi-major axis (km)');title('Mars-centric');
%% hamiltonian
COSTSQRT=sqrt(LVX.^2+LVY.^2+LVZ.^2);
costsqrt=sqrt(lvx.^2+lvy.^2+lvz.^2);
L=(COSTSQRT./M)-(1-LM)/g0Isp;
l=(costsqrt./m)-(1-lm)/g0Isp;
T=kfac;
K=-(T./M)./COSTSQRT;
AX=K.*LVX.*heaviside(L);AY=K.*LVY.*heaviside(L);AZ=K.*LVZ.*heaviside(L);
k=-(T./m)./costsqrt;
ax=k.*lvx.*heaviside(l);ay=k.*lvy.*heaviside(l);az=k.*lvz.*heaviside(l);
figure(1)
subplot(1,3,1)
hold on
quiver3(x,y,z,ax,ay,az,0.0025);
quiver3(x,y,z,vx,vy,vz,0.0025);
subplot(1,3,3)
hold on
quiver3(X,Y,Z,AX,AY,AZ,0.0125);
quiver3(X,Y,Z,VX,VY,VZ,0.0125);
%% Hamiltonian
r_soi=a_earthsun*(mu_e/mu_s)^0.4;
ACC=sqrt(AX.^2+AY.^2+AZ.^2);
acc=sqrt(ax.^2+ay.^2+az.^2);
H=(1-LM).*M.*ACC/g0Isp+LX.*VX+LY.*VY+LZ.*VZ+LVX.*(-mu_s.*X./R.^3+AX)+LVY.*(-mu_s.*Y./R.^3+AY)+LVZ.*(-mu_s.*Z./R.^3+AZ);
h=(1-lm).*m.*acc/g0Isp+lx.*vx+ly.*vy+lz.*vz+lvx.*(-mu_e.*x./r.^3+ax)+lvy.*(-mu_e.*y./r.^3+ay)+lvz.*(-mu_e.*z./r.^3+az);
Hinteg=H.*heaviside(r-r_soi)+h.*heaviside(r_soi-r);
figure(5)
subplot(3,1,1)
plot(t/86400,H,'k');
grid on;grid minor;
xlabel('Time');ylabel('Hamiltonian');title('Heliocentric');
subplot(3,1,2)
plot(t/86400,h,'k');
grid on;grid minor;
xlabel('Time');ylabel('Hamiltonian');title('Geocentric');
subplot(3,1,3)
plot(t/86400,Hinteg,'k');
grid on;grid minor;
xlabel('Time');ylabel('Hamiltonian');title('Integrated Helio-Geo');
%% costates
figure(6)
subplot(3,2,1)
plot(t/86400,LX,'k');
grid on;grid minor;title('\lambda_x');
subplot(3,2,2)
plot(t/86400,LY,'k');
grid on;grid minor;title('\lambda_y');
subplot(3,2,3)
plot(t/86400,LZ,'k');
grid on;grid minor;title('\lambda_z');
subplot(3,2,4)
plot(t/86400,LVX,'k');
grid on;grid minor;title('\lambda_{v_x}');
subplot(3,2,5)
plot(t/86400,LVY,'k');
grid on;grid minor;title('\lambda_{v_y}');
subplot(3,2,6)
plot(t/86400,LVZ,'k');
grid on;grid minor;title('\lambda_{v_z}');
%% thrust
figure(7)
plot(t/86400,M.*ACC,'k')
grid on;grid minor;
xlabel('Time (days)');ylabel('Thrust (N)');
%% conserved qty
figure(8)
plot(t/86400,M.*(1-LM),'k')
grid on;grid minor;
xlabel('Time (days)');ylabel('Conserved qty');title('m(1-\lambda_m)');