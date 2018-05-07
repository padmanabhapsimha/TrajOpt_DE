close all;clear;clc
data=importdata('Output_dispOnly.txt');
sys=importdata('System_params.txt');
%%
% data(:,end)=data(:,end)*86400;
close all;
figure(1);
plot3(data(:,1),data(:,2),data(:,3),'k.-','MarkerSize',4.5);
view([0 0 90]);
hold on
T=sys(5);
mu=sys(1);
% mu=398600e+9;
% mu=4.282837e+13;
g0Isp=sys(3)*sys(4);
re=sys(9)*sys(6);
ecci=sys(7);
eccf=sys(12);
omgaf=sys(13)*pi/180;
inclf=sys(19)*pi/180;
raanf=sys(20)*pi/180;
raani=sys(21)*pi/180;
omgai=sys(22)*pi/180;
incli=sys(23)*pi/180;
rm=sys(9)*sys(10);
nu=(0:0.001:2*pi);
ai=re/(1-ecci);
af=rm/(1-eccf);
ri=ai*(1-ecci^2)./(1+ecci*cos(nu));
rf=af*(1-eccf^2)./(1+eccf*cos(nu));
rfa=af*(1-eccf);rfd=af*(1+eccf);
xf=rf.*(cos(raanf)*cos(omgaf+nu)-sin(raanf)*sin(omgaf+nu)*cos(inclf));
yf=rf.*(sin(raanf)*cos(omgaf+nu)+cos(raanf)*sin(omgaf+nu)*cos(inclf));
zf=rf.*(sin(omgaf+nu)*sin(inclf));
xi=ri.*(cos(raani)*cos(omgai+nu)-sin(raani)*sin(omgai+nu)*cos(incli));
yi=ri.*(sin(raani)*cos(omgai+nu)+cos(raani)*sin(omgai+nu)*cos(incli));
zi=ri.*(sin(omgai+nu)*sin(incli));
plot3(xi,yi,zi,'b','LineWidth',1.5)
plot3(xf,yf,zf,'r','LineWidth',1.5)
plot3(0,0,0,'rp','MarkerSize',20,'MarkerFaceColor',[0.9333 0.9333 0])
plot3(data(end,1),data(end,2),data(end,3),'bo','MarkerSize',5,'MarkerFaceColor',[0 0 0])
xlabel('x axis (m)');ylabel('y axis (m)');zlabel('z axis (m)');title('Trajectory')
% % % 3d angle visualization help
% xpt=linspace(-max(xf),max(xf),250);ypt=linspace(-max(yf),max(yf),250);
% [X,Y]=meshgrid(xpt,ypt);
% xa=rfa.*(cos(raanf)*cos(omgaf+0)-sin(raanf)*sin(omgaf+0)*cos(inclf));
% ya=rfa.*(sin(raanf)*cos(omgaf+0)+cos(raanf)*sin(omgaf+0)*cos(inclf));
% za=rfa.*(sin(omgaf+0)*sin(inclf));
% xd=rfd.*(cos(raanf)*cos(omgaf+pi)-sin(raanf)*sin(omgaf+pi)*cos(inclf));
% yd=rfd.*(sin(raanf)*cos(omgaf+pi)+cos(raanf)*sin(omgaf+pi)*cos(inclf));
% zd=rfd.*(sin(omgaf+pi)*sin(inclf));
% plot3(xa,ya,za,'bo','MarkerSize',5,'MarkerFaceColor',[1 0 0])
% plot3(xd,yd,zd,'bo','MarkerSize',5,'MarkerFaceColor',[0 1 0])
% plot3(X,Y,zeros(size(X)),'k.')
% plot3(X,zeros(size(X)),Y,'k.')
% plot3(zeros(size(X)),X,Y,'k.')
% % % more plot options
axis equal;grid on;grid minor;
%%
figure(2)
plot(data(:,end)/86400,data(:,7),'k.-')
xlabel('Days');ylabel('Spacecraft mass (kg)')
title('Spacecraft mass variation');
grid on;grid minor;
disp(data(1,7)-data(end,7));
%%
r=sqrt(data(:,1).^2+data(:,2).^2+data(:,3).^2);
v=sqrt(data(:,4).^2+data(:,5).^2+data(:,6).^2);
figure(3)
subplot(2,1,1)
plot(data(:,end)/86400,r/1000,'k-')
grid on;grid minor;
xlabel('Days');ylabel('Distance (km)');
title('Distance vs time');
subplot(2,1,2)
plot(data(:,end)/86400,v/1000,'k-')
grid on;grid minor;
xlabel('Days');ylabel('Velocity');
title('Velocity vs time (km/s vs days)');
%%
figure(4)
subplot(2,1,1)
plot(data(:,end)/86400,data(:,3)/1000,'k-')
grid on;grid minor;
xlabel('Days');ylabel('Out of plane distance (km)');
title('Out of plane movement');
subplot(2,1,2)
plot(data(:,end)/86400,data(:,6),'k-')
grid on;grid minor;
xlabel('Days');ylabel('Out of plane speed (m/s)');
title('Out of plane velocity');
%%
figure(5)
subplot(2,3,1)
plot(data(:,end)/86400,data(:,8),'k-')
grid on;grid minor;
xlabel('Days');ylabel('\lambda_x');
title('\lambda_x vs time')
subplot(2,3,2)
plot(data(:,end)/86400,data(:,9),'k-')
grid on;grid minor;
xlabel('Days');ylabel('\lambda_y');
title('\lambda_y vs time')
subplot(2,3,3)
plot(data(:,end)/86400,data(:,11),'k-')
grid on;grid minor;
xlabel('Days');ylabel('\lambda_v_x');
title('\lambda_v_x vs time')
subplot(2,3,4)
plot(data(:,end)/86400,data(:,12),'k-')
grid on;grid minor;
xlabel('Days');ylabel('\lambda_v_y');
title('\lambda_v_y vs time')
subplot(2,3,5)
plot(data(:,end)/86400,data(:,14),'k-')
grid on;grid minor;
xlabel('Days');ylabel('\lambda_m');
title('\lambda_m vs time')

subplot(2,3,6)
line(data(:,end)/86400,data(:,13),'Color','r')
ax1 = gca; % current axes
ax1.XColor = 'r';
ax1.YColor = 'r';
ax1_pos = ax1.Position; % position of first axes
xlabel('Days');ylabel('\lambda_v_x');
ax2 = axes('Position',ax1_pos,...
'XAxisLocation','top',...
'YAxisLocation','right',...
'Color','none');
hold on;
line(data(:,end)/86400,data(:,10),'Parent',ax2,'Color','k')
ylabel('\lambda_x');
grid on;grid minor


title('\lambda_z and \lambda_v_z vs time')
l=(sqrt(data(:,11).^2+data(:,12).^2+data(:,13).^2)./data(:,7))-(1-(data(:,14))/g0Isp);
% solar electric propulsion
% T=sys(5).*(sys(9)./r).^2;
% for i=1:1:length(T)
% if T(i)>0.236
%     T(i)=0.236;
% end
% if T(i)<0.040
%     T(i)=0.040;
% end
% end
%
rf=-(T./data(:,7))./sqrt(data(:,11).^2+data(:,12).^2+data(:,13).^2);
ax=zeros(length(rf),1);ay=ax;az=ax;
for i=1:1:length(l)
if l(i)>=0.0
	ax(i)=rf(i).*data(i,11);ay(i)=rf(i).*data(i,12);az(i)=rf(i).*data(i,13);
else
	ax(i)=0.0;ay(i)=0.0;az(i)=0.0;
end
end
accl=sqrt((ax).^2+(ay).^2+(az).^2);
figure(1)
% plot3(xi,yi,zi,'b','LineWidth',1.5)
% hold on
% plot3(xf,yf,zf,'r','LineWidth',1.5)
% quiver3(data(:,1),data(:,2),data(:,3),ax,ay,az,0.75,'Color',[0.04314 0.6471 0]);
% grid on;grid minor;view([0 0 90]);
% legend('Init orbit','Final orbit','Thrust vector');axis equal;
%%
figure(6)
hx=data(:,2).*data(:,6)-data(:,3).*data(:,5);
hy=data(:,3).*data(:,4)-data(:,1).*data(:,6);
hz=data(:,1).*data(:,5)-data(:,2).*data(:,4);
h=sqrt(hx.^2+hy.^2+hz.^2);
plot(data(:,end)/86400,acosd(hz./h),'k-')
grid on;grid minor;
xlabel('Days');ylabel('Inclination (degrees)');
title('Inclination vs time');
%%
a=1./((2./r)-(v.^2/mu));
figure(7)
plot(data(:,end)/86400,a/1000,'k-')
grid on;grid minor;
xlabel('Days');ylabel('Semimajor axis (km)');
title('Semimajor axis vs time');

%%
% figure(8)
% conqty=data(:,7).*data(:,14);
% plot(data(:,end)/86400,conqty,'k-')
% grid on;grid minor;
% xlabel('Days');ylabel('Conserved qty');
% title('Variation of conserved qty')
% ylim([0.999975*min(conqty) 1.000025*max(conqty)+1])
%%
figure(9)
F=data(:,7).*accl;
plot(data(:,end)/86400,F,'k-')
grid on;grid minor;
xlabel('Days'),ylabel('Thrust (mN)')
title('Thrust profile')
%%
figure(10)
ex=((data(:,5).*hz-data(:,6).*hy)/mu)-data(:,1)./r;
ey=((data(:,6).*hx-data(:,4).*hz)/mu)-data(:,2)./r;
ez=((data(:,4).*hy-data(:,5).*hx)/mu)-data(:,3)./r;
eccen=sqrt(ex.^2+ey.^2+ez.^2);
plot(data(:,end)/86400,eccen,'k-')
grid on;grid minor;
xlabel('Days'),ylabel('Eccentricity')
title('Eccentricity profile')
%%
% figure(11)
% plot(data(:,end)/86400,atan2d(az,(sqrt(ax.^2+ay.^2+az.^2))),'k-')
% grid on;grid minor
% figure(12)
% plot(data(:,end)/86400,atan2d(ay,ax),'k-')
% grid on;grid minor