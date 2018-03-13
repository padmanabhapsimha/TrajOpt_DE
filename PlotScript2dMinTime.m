close all;clear;
data=importdata('Output.txt');
sys=importdata('System_params.txt');
%%
close all;
figure(1);
plot(data(:,1),data(:,2),'.-.k');
view([0 0 90]);
hold on
T=sys(5);
mu=sys(1);
g0Isp=sys(3)*sys(4);
re=sys(9)*sys(6);
ecci=sys(7);
eccf=sys(12);
omgaf=sys(13)*pi/180;
rm=sys(9)*sys(10);
nu=0:0.001:2*pi;
ai=re/(1-ecci);
af=rm/(1-eccf);
ra=ai*(1-ecci^2)./(1+ecci*cos(nu));
rf=af*(1-eccf^2)./(1+eccf*cos(nu));
plot(ra.*cos(nu),ra.*sin(nu),'Color',[0 0.4039 0.9333],'LineWidth',1.5)
plot(rf.*cos(nu+omgaf),rf.*sin(nu+omgaf),'r','LineWidth',1.5)
plot(0,0,'rp','MarkerSize',20,'MarkerFaceColor',[0.9333 0.9333 0])
xlabel('x axis (m)');ylabel('y axis (m)');title('Trajectory')
axis equal;grid on;grid minor;
%%
figure(2)
plot(data(:,end)/86400,data(:,5),'k.-')
xlabel('Days');ylabel('Spacecraft mass (kg)')
title('Spacecraft mass variation');
grid on;grid minor;
disp(data(1,5)-data(end,5));
r=sqrt(data(:,1).^2+data(:,2).^2);
v=sqrt(data(:,3).^2+data(:,4).^2);
%%
figure(3)
subplot(2,1,1)
plot(data(:,end)/86400,r/re,'k.-')
grid on;grid minor;
xlabel('Days');ylabel('AU from Center');
title('Distance vs time');
subplot(2,1,2)
plot(data(:,end)/86400,v/1000,'k.-')
grid on;grid minor;
xlabel('Days');ylabel('Speed wrt Center');
title('Speed vs time (km/s vs days)');
%%
figure(4)
subplot(2,3,1)
plot(data(:,end)/86400,data(:,6),'k.-')
grid on;grid minor;
xlabel('Days');ylabel('\lambda_x');
title('\lambda_x vs time')
subplot(2,3,2)
plot(data(:,end)/86400,data(:,7),'k.-')
grid on;grid minor;
xlabel('Days');ylabel('\lambda_y');
title('\lambda_y vs time')
subplot(2,3,3)
plot(data(:,end)/86400,data(:,8),'k.-')
grid on;grid minor;
xlabel('Days');ylabel('\lambda_v_x');
title('\lambda_v_x vs time')
subplot(2,3,4)
plot(data(:,end)/86400,data(:,9),'k.-')
grid on;grid minor;
xlabel('Days');ylabel('\lambda_v_y');
title('\lambda_v_y vs time')
subplot(2,3,5)
plot(data(:,end)/86400,data(:,10),'k.-')
grid on;grid minor;
xlabel('Days');ylabel('\lambda_m');
title('\lambda_m vs time')
subplot(2,3,6)

lambdamdot=-T.*sqrt(data(:,8).^2+data(:,9).^2)./data(:,5).^2;
plot(data(:,end)/86400,lambdamdot,'k.-')
grid on;grid minor;
xlabel('Days');ylabel('\lambda_m rate');
title('Rate of change of \lambda_m vs time')
%%
figure(5)
subplot(2,1,1)
plot(data(:,end)/86400,T./T,'k.-')
grid on;grid minor;
xlabel('Days');ylabel('Throttle factor');
title('Throttle factor vs time');
subplot(2,1,2)
alpha=pi+atan2(data(:,9),data(:,8));
plot(data(:,end)/86400,alpha*180/pi,'k.-')
grid on;grid minor;
xlabel('Days');ylabel('Thrust angle absolute');
title('Absolute thrust angle vs time (degrees)');
figure(1)
hold on
quiver(data(:,1),data(:,2),T.*cos(alpha),T.*sin(alpha));
figure(5)
%%
H=1 + data(:,6).*data(:,3) + data(:,7).*data(:,4) + data(:,8).*(-mu.*data(:,1)./r.^3 + T.*cos(alpha)./data(:,5));
H=H+data(:,9).*(-mu.*data(:,2)./r.^3 + T.*sin(alpha)./data(:,5))+data(:,10).*(-T./g0Isp);
figure(6)
% plot(data(:,end)/86400,data(:,6).*data(:,3),'k.-')
% hold on
% plot(data(:,end)/86400,data(:,7).*data(:,4),'r.-')
% plot(data(:,end)/86400,data(:,8).*(-mu.*data(:,1)./r.^3 + T.*cos(alpha)./data(:,5)),'b.-')
% plot(data(:,end)/86400,data(:,9).*(-mu.*data(:,2)./r.^3 + T.*sin(alpha)./data(:,5)),'m.-')
% plot(data(:,end)/86400,data(:,10).*(-T./g0Isp),'c.-')
plot(data(:,end)/86400,H,'ko-')
grid on;grid minor;
xlabel('Days');ylabel('Hamiltonian');
title('Hamiltonian vs time');
%%
% figure(9)
% angt=alpha*180/pi;
% angv=atan2d(data(:,4),data(:,3));
% angr=atan2d(data(:,2),data(:,1));
% plot(data(:,end)/86400,(angt-angr)-360*heaviside(angt-angr-180),'k.-')
% ylim([-180 180])
% dat1=importdata('H:\Sem _8\Work\Week2\Validation\Img_dat.txt');
% dat2=importdata('H:\Sem _8\Work\Week2\Validation\Img_dat_2.txt');
% dat=[dat1;dat2];
% hold on
% plot(dat(:,1)*180/pi,dat(:,2)*360/pi,'r-');
%%
figure
plot(data(:,end)/86400,r-2410000,'k.-')
disp(max((r-2410000)/1000));
grid on;grid minor
figure(1)
hold on
quiver(data(end,1),data(end,2),data(end,3),data(end,4))
%%
disp(atan2d(data(end,2),data(end,1)))