%% Spacecraft Guidance and Navigation - Uncertainty propagation (2022/2023)
% Assignment:     2
% Exercise:       Uncertainty propagation
% Author:         Giovanni Fereoli
% ID number:      990939

%% Ex.1
clear; clc; close;

%Data and options
r0=[6054.30795817484; -3072.03883303992; -133.115352431876];
v0=[4.64750094824087; 9.18608475681236; -0.62056520749034];
x0=[r0; v0];
P0=[5.6*1e-3, 3.5*1e-3, -7.1*1e-4, 0, 0, 0;...
    3.5*1e-3, 9.7*1e-3, 7.6*1e-4, 0, 0, 0;...
    -7.1*1e-4, 7.6*1e-4, 8.1*1e-4, 0, 0, 0;...
    0, 0, 0, 2.8*1e-7, 0, 0;...
    0, 0, 0, 0,2.7*1e-7, 0;...
    0, 0, 0, 0, 0, 9.6*1e-8];
mu=398600;
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);

%Period orbit
a=1/(2/norm(r0)-norm(v0)^2/mu);
T=2*pi*sqrt(a^3/mu);

%Propagation orbit
tof=4*T;
tvec_orb=linspace(0,tof,2000);
[ttM, xxM] = ode113(@(t,xM) STMexact_RTBP(t,xM,mu),...
tvec_orb,[x0; reshape(eye(6),[],1)], options);
xx=xxM(:,1:6);    
STMend=(reshape(xxM(end,7:42),6,6))';

%Plot orbit
figure(1);
[X, Y, Z] = sphere;
EarthRadius = 6371;
hSurface = surf(X*EarthRadius(1), Y*EarthRadius(1), Z*EarthRadius(1));
set(hSurface,'FaceColor',[0.5 0.5 0.5]);
hold on;
plot3(xx(:,1), xx(:,2),xx(:,3),'b','LineWidth',1);
xlabel('x [Km]');
ylabel('y [Km]');
zlabel('z [Km]');
grid on; grid minor;
axis equal;

%UP with LinCov
tvec=[0,tof/8,tof/4,3*tof/8,tof/2,5*tof/8,6*tof/8,7*tof/8,tof];
[XmeanLC, PLC]=LinCov(x0, P0, tvec);
r_dispLC=zeros(length(XmeanLC),1);
v_dispLC=zeros(length(XmeanLC),1);
for i=1:length(XmeanLC)
    r_dispLC(i)=sqrt(trace(PLC(1:3,1:3,i)));
    v_dispLC(i)=sqrt(trace(PLC(4:6,4:6,i)));
end

%UP with UT
[XmeanUT, PUT]=UT(x0, P0, tvec);
r_dispUT=zeros(length(XmeanUT),1);
v_dispUT=zeros(length(XmeanUT),1);
for i=1:length(XmeanLC)
    r_dispUT(i)=sqrt(trace(PUT(1:3,1:3,i)));   
    v_dispUT(i)=sqrt(trace(PUT(4:6,4:6,i)));
end

%% Ex.2

%UP with MCM
[XmeanMCM, PMCM, Ysamples]=MCM(x0, P0, tvec);
r_dispMCM=zeros(length(XmeanMCM),1);
v_dispMCM=zeros(length(XmeanMCM),1);
for i=1:length(XmeanMCM)
    r_dispMCM(i)=sqrt(trace(PMCM(1:3,1:3,i)));
    v_dispMCM(i)=sqrt(trace(PMCM(4:6,4:6,i)));
end

%Plots covariance ellipsoid at each perigee/apogee in ECI
for j=2:length(tvec)
    %Clouds sampled points with 1-sigma
    disc=500;
    X_cloudLC=mvnrnd(XmeanLC(j,1:3),PLC(1:3,1:3,j),disc);
    X_cloudUT=mvnrnd(XmeanUT(j,1:3),...
        0.5*(PUT(1:3,1:3,j)+PUT(1:3,1:3,j)'),disc);
    X_cloudMCM=mvnrnd(XmeanMCM(j,1:3),PMCM(1:3,1:3,j),disc);
    
    %Covariance ellipsoids plots with MCM samples
    figure(j+1);
    subplot(1,3,1);
    scatter3(Ysamples(j,1,:),Ysamples(j,2,:),Ysamples(j,3,:),3,'k','*');
    hold on;
    scatter3(X_cloudLC(:,1),X_cloudLC(:,2),X_cloudLC(:,3),5,'b','*');
    plot3(XmeanLC(j,1),XmeanLC(j,2),XmeanLC(j,3),'r.',...
        'MarkerSize',25);
    box on;
    hold on;
    xlabel('x [km]','Interpreter','latex');
    ylabel('y [km]','Interpreter','latex');
    zlabel('z [km]','Interpreter','latex');
    legend('MC Points', 'Sampled Cov', 'Mean','Interpreter','latex',...
        'Location','southeast');
    grid on; grid minor;
    view(2);
    title('LinCov','Interpreter','latex');
    subplot(1,3,2)
    scatter3(Ysamples(j,1,:),Ysamples(j,2,:),Ysamples(j,3,:),3,'k','*');
    hold on;
    scatter3(X_cloudUT(:,1),X_cloudUT(:,2),X_cloudUT(:,3),5,'b','*');
    plot3(XmeanUT(j,1),XmeanUT(j,2),XmeanUT(j,3),'r.',...
        'MarkerSize',25);
    box on;
    xlabel('x [km]','Interpreter','latex');
    ylabel('y [km]','Interpreter','latex');
    zlabel('z [km]','Interpreter','latex');
    legend('MC Points', 'Sampled Cov', 'Mean','Interpreter','latex',...
        'Location','southeast');
    grid on; grid minor;
    view(2);
    title('UT','Interpreter','latex');
    subplot(1,3,3)
    scatter3(Ysamples(j,1,:),Ysamples(j,2,:),Ysamples(j,3,:),3,'k','*');
    hold on;
    scatter3(X_cloudMCM(:,1),X_cloudMCM(:,2),X_cloudMCM(:,3),5,'b','*');
    plot3(XmeanMCM(j,1),XmeanMCM(j,2),XmeanMCM(j,3),'r.',...
        'MarkerSize',25);
    box on;
    xlabel('x [km]','Interpreter','latex');
    ylabel('y [km]','Interpreter','latex');
    zlabel('z [km]','Interpreter','latex');
    legend('MC Points', 'Sampled Cov', 'Mean','Interpreter','latex',...
        'Location','southeast');
    grid on; grid minor;
    view(2);
    title('MCM','Interpreter','latex');
end

%% Functions

function dxdt=RTBP(~,xx,mu)
%Initialize
dxdt=zeros(6,1);

%Unpacking
x=xx(1);
y=xx(2);
z=xx(3);
xdot=xx(4);
ydot=xx(5);
zdot=xx(6);

% RTBP dynamics
r_norm=sqrt(x^2+y^2+z^2);
dxdt(1:3)=[xdot;ydot;zdot];
dxdt(4:6)=-mu/r_norm^3*[x;y;z];
end

function dxMdt=STMexact_RTBP(~,xM,mu)
%Initialize
dxMdt=zeros(42,1);

%Unpacking
x=xM(1);
y=xM(2);
z=xM(3);
xdot=xM(4);
ydot=xM(5);
zdot=xM(6);
M=(reshape(xM(7:42),6,6))';   

%RTBP dynamics
r_norm=sqrt(x^2+y^2+z^2);
dxMdt(1:3)=[xdot;ydot;zdot];
dxMdt(4:6)=-mu/r_norm^3*[x;y;z];

%Variational equations, STM dynamics
A1=zeros(3);
A2=eye(3);
A3=(3*mu/r_norm^5)*[x^2, x*y, x*z;...
    x*y, y^2, y*z;...
    x*z, y*z, z^2]-(mu/r_norm^3)*eye(3);
A4=zeros(3);
A=[A1, A2;...
    A3, A4];

dMdt=A*M;
dxMdt(7:12)=dMdt(1,1:6)';
dxMdt(13:18)=dMdt(2,1:6)';
dxMdt(19:24)=dMdt(3,1:6)';
dxMdt(25:30)=dMdt(4,1:6)';
dxMdt(31:36)=dMdt(5,1:6)';
dxMdt(37:42)=dMdt(6,1:6)';
end

function [Xmean, P]=LinCov(X0, P0, tvec)
%Propagation
mu=398600;
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~, xxM] = ode113(@(t,xM) STMexact_RTBP(t,xM,mu),...
tvec,[X0; reshape(eye(6),[],1)], options);
X=xxM(:,1:6);

%Mean of LinCov method
Xmean=X;

%Covarianze of LinCov method
P=zeros(6,6,length(tvec));
for i=1:length(tvec)
   STM=(reshape(xxM(i,7:42),6,6))'; 
   P(:,:,i)=STM*P0*STM';
end
end

function [Xmean, P]=UT(X0, P0, tvec)
%Initialization
n=6;
alpha=1e-3;
k=0;
lambda=alpha^2*(n+k)-n;
beta=2;
mu=398600;

%Sigma initial points
UTmatrix=sqrtm((n+lambda)*P0);
chi0=X0;
chi1=X0+UTmatrix(:,1);
chi2=X0+UTmatrix(:,2);
chi3=X0+UTmatrix(:,3);
chi4=X0+UTmatrix(:,4);
chi5=X0+UTmatrix(:,5);
chi6=X0+UTmatrix(:,6);
chi7=X0-UTmatrix(:,1);
chi8=X0-UTmatrix(:,2);
chi9=X0-UTmatrix(:,3);
chi10=X0-UTmatrix(:,4);
chi11=X0-UTmatrix(:,5);
chi12=X0-UTmatrix(:,6);

%Weigths
w0m=lambda/(n+lambda);
w0p=lambda/(n+lambda)+(1-alpha^2+beta);
wim=1/(2*(n+lambda));
wip=1/(2*(n+lambda));

%Propagation
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
%chi0
[~, XX] = ode113(@(t,xx) RTBP(t,xx,mu),tvec,chi0, options);
jota0=XX;  
%chi1
[~, XX] = ode113(@(t,xx) RTBP(t,xx,mu),tvec,chi1, options);
jota1=XX;   
%chi2
[~, XX] = ode113(@(t,xx) RTBP(t,xx,mu),tvec,chi2, options);
jota2=XX; 
%chi3
[~, XX] = ode113(@(t,xx) RTBP(t,xx,mu),tvec,chi3, options);
jota3=XX; 
%chi4
[~, XX] = ode113(@(t,xx) RTBP(t,xx,mu),tvec,chi4, options);
jota4=XX;  
%chi5
[~, XX] = ode113(@(t,xx) RTBP(t,xx,mu),tvec,chi5, options);
jota5=XX; 
%chi6
[~, XX] = ode113(@(t,xx) RTBP(t,xx,mu),tvec,chi6, options);
jota6=XX; 
%chi6
[~, XX] = ode113(@(t,xx) RTBP(t,xx,mu),tvec,chi7, options);
jota7=XX; 
%chi6
[~, XX] = ode113(@(t,xx) RTBP(t,xx,mu),tvec,chi8, options);
jota8=XX; 
%chi6
[~, XX] = ode113(@(t,xx) RTBP(t,xx,mu),tvec,chi9, options);
jota9=XX; 
%chi6
[~, XX] = ode113(@(t,xx) RTBP(t,xx,mu),tvec,chi10, options);
jota10=XX; 
%chi6
[~, XX] = ode113(@(t,xx) RTBP(t,xx,mu),tvec,chi11, options);
jota11=XX; 
%chi6
[~, XX] = ode113(@(t,xx) RTBP(t,xx,mu),tvec,chi12, options);
jota12=XX; 

%Sample Mean and Sample Covariance
Xmean=zeros(length(tvec),6);
P=zeros(6,6,length(tvec));
for i=1:length(tvec)
Xmean(i,:)=w0m*jota0(i,:)+wim*(jota1(i,:)+jota2(i,:)+jota3(i,:)+...
    jota4(i,:)+jota5(i,:)+jota6(i,:)+jota7(i,:)+jota8(i,:)+jota9(i,:)+...
    jota10(i,:)+jota11(i,:)+jota12(i,:));
P(:,:,i)=w0p*(jota0(i,:)-Xmean(i,:))'*(jota0(i,:)-Xmean(i,:))+...
    wip*((jota1(i,:)-Xmean(i,:))'*(jota1(i,:)-Xmean(i,:))+...
    (jota2(i,:)-Xmean(i,:))'*(jota2(i,:)-Xmean(i,:))+...
    (jota3(i,:)-Xmean(i,:))'*(jota3(i,:)-Xmean(i,:))+...
    (jota4(i,:)-Xmean(i,:))'*(jota4(i,:)-Xmean(i,:))+...
    (jota5(i,:)-Xmean(i,:))'*(jota5(i,:)-Xmean(i,:))+...
    (jota6(i,:)-Xmean(i,:))'*(jota6(i,:)-Xmean(i,:))+...
    (jota7(i,:)-Xmean(i,:))'*(jota7(i,:)-Xmean(i,:))+...
    (jota8(i,:)-Xmean(i,:))'*(jota8(i,:)-Xmean(i,:))+...
    (jota9(i,:)-Xmean(i,:))'*(jota9(i,:)-Xmean(i,:))+...
    (jota10(i,:)-Xmean(i,:))'*(jota10(i,:)-Xmean(i,:))+...
    (jota11(i,:)-Xmean(i,:))'*(jota11(i,:)-Xmean(i,:))+...
    (jota12(i,:)-Xmean(i,:))'*(jota12(i,:)-Xmean(i,:)));
end
end

function [Xmean, P, Ysamples]=MCM(X0, P0, tvec)
%Initialization
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
mu=398600;

%Sample generation
n = 1000;
Xsamples=mvnrnd(X0, P0, n); 

%Sample propagation
Ysamples=zeros(length(tvec),6,n);
for i=1:n
[~, XX] = ode113(@(t,xx) RTBP(t,xx,mu),tvec,Xsamples(i,:)', options);
Ysamples(:,:,i)=XX;
end

%Sample Mean and Sample Covariance
Xmean=zeros(length(tvec),6);
P=zeros(6,6,length(tvec));
for i=1:length(tvec)
    Ysamples_it=(squeeze(Ysamples(i,:,:)))';
    Xmean(i,:)=mean(Ysamples_it);                                           
    P(:,:,i)=cov(Ysamples_it);
end
end
