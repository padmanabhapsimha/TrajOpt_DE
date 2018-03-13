close all;clear;
data=importdata('Output.txt');
sys=importdata('System_params.txt');
%%
% close all;
figure(1);
plot(data(:,1),data(:,2),'k.-','MarkerSize',4.5);
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
plot(data(end,1),data(end,2),'bo','MarkerSize',5,'MarkerFaceColor',[0 0 0])
xlabel('x axis (m)');ylabel('y axis (m)');title('Trajectory')
axis equal;grid on;grid minor;
%%
figure(2)
plot(data(:,end)/86400,data(:,5),'k.-')
xlabel('Days');ylabel('Spacecraft mass (kg)')
title('Spacecraft mass variation');
grid on;grid minor;
disp(100*(data(1,5)-data(end,5))/data(1,5));
r=sqrt(data(:,1).^2+data(:,2).^2);
v=sqrt(data(:,3).^2+data(:,4).^2);
%%
figure(3)
subplot(2,1,1)
plot(data(:,end)/86400,r/sys(9),'k.-')
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
xlabel('Days');ylabel('\lambda_m^`');
title('\lambda_m^` vs time')
subplot(2,3,6)

l1=((sqrt((data(:,8).^2+data(:,9).^2))./data(:,5))-(data(:,10)/g0Isp));
p=ones(length(l1),1);
alpha=pi+atan2(data(:,9),data(:,8));

for i=1:1:length(l1)
    if l1(i)>=0
        p(i)=1.0;
    else
        p(i)=0.0;
    end
end
lambdamdot=p.*data(:,10).*sqrt(data(:,8).^2+data(:,9).^2)./(g0Isp*data(:,5));
plot(data(:,end)/86400,lambdamdot,'k.-')
grid on;grid minor;
xlabel('Days');ylabel('\lambda_m^` rate');
title('Rate of change of \lambda_m^` vs time')
%%
figure(5)
subplot(2,1,1)
plot(data(:,end)/86400,p,'k.-')
grid on;grid minor;
xlabel('Days');ylabel('Throttle factor');
title('Throttle factor vs time');
subplot(2,1,2)
plot(data(:,end)/86400,alpha*180/pi,'k.-')
grid on;grid minor;
xlabel('Days');ylabel('Thrust angle absolute');
title('Thrust angle vs time (degrees)');
figure(1)
hold on
quiver(data(:,1),data(:,2),p.*cos(alpha),p.*sin(alpha),0.75,'Color',[0.04314 0.6471 0]);
figure(5)
%%
a=1./((2./r)-(v.^2/mu));
figure(7)
plot(data(:,end)/86400,a/1000,'k.-')
grid on;grid minor;
xlabel('Days');ylabel('Semimajor axis');
title('Semimajor axis vs time');
%%
figure(8)
conqty=data(:,5).*data(:,10);
plot(data(:,end)/86400,conqty,'k.-')
grid on;grid minor;
xlabel('Days');ylabel('Conserved qty');
title('Variation of conserved qty')
% ylim([0.999975*min(conqty) 1.000025*max(conqty)+1])
