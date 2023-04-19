%% Spacecraft Guidance and Navigation - Sequential filtering (2022/2023)
% Assignment:     2
% Exercise:       Sequential filtering
% Author:         Giovanni Fereoli
% ID number:      990939

%% Ex.1
clear; close all; clc; cspice_kclear;

%Add kernels and paths
cspice_furnsh('assignment02.tm');
addpath('sgp4');
addpath('simulator_ex3');

%Data GEO
mu=398600;
R=42241.08;
n=sqrt(mu/R^3);

%Reference Time
t0_str='UTC 2023/04/01-14:55:12.023';
tf_str='UTC 2023/04/02-14:55:12.023';
t0=cspice_str2et(t0_str);
tf=cspice_str2et(tf_str);

%Initial attitude dynamics
omega0=[-0.001262427155865; 0.001204540074343; -0.000039180139156]; 
q0=[0.674156352338764; 0.223585877389611; 0.465489474399161; ...
    0.528055032413102]; 
x0_att=[omega0; q0];

%Initial relative mean and covariance
r0_mean=[15.792658268071492;-59.044939772661586;3.227106250277039];
v0_mean=[-0.053960274403210;-0.053969644762889;-0.089140748762173];
x0_mean=[r0_mean;v0_mean];
P0=diag([10,10,10,0.1,0.1,0.1]);

%Initial relative reference state
r0_ref= [12; -60; 0];
v0_ref=[1e-4; -2*n*r0_ref(1); -1.2e-3];
x0_ref=[r0_ref; v0_ref];

%Camera data
Cam.Cframe=[1, 0, 0;0, 0, -1;0, 1, 0];
Cam.f=30; %[mm]
Cam.d=54; %[pix/mm]
Cam.b=1; %[m]
Cam.p0=[960;600];
Sensor_size=[1920;1200]; %[pix]
Cam.R = diag([10,10,10]); %[pix2]

%Propagation SGN attitude dynamics
tvec=t0:1:tf;
disc=length(tvec);
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~, XX_att] = ode113(@(t,X) dynamicsSGN(t,X), tvec, x0_att, options);

%Unpacking SGN attitude dynamics
omega_vec=XX_att(:,1:3);
q_vec=XX_att(:,4:7);
Adcm=zeros(3,3,disc);
for i=1:disc
    Adcm(:,:,i)=quat2dcm(q_vec(i,:));
end

%Plot SGN attitude dynamics
figure(1);
plot((tvec-t0)/60^2, omega_vec(:,1),'b','LineWidth',1);
hold on;
plot((tvec-t0)/60^2, omega_vec(:,2),'r','LineWidth',1);
plot((tvec-t0)/60^2, omega_vec(:,3),'k','LineWidth',1);
xlabel('Time [h]','Interpreter','latex');
ylabel('$\mathbf{\omega}$ [rad/s]','Interpreter','latex');
grid on; grid minor;
figure(2);
plot((tvec-t0)/60^2, q_vec(:,1),'b','LineWidth',1);
hold on;
plot((tvec-t0)/60^2, q_vec(:,2),'r','LineWidth',1);
plot((tvec-t0)/60^2, q_vec(:,3),'k','LineWidth',1);
plot((tvec-t0)/60^2, q_vec(:,4),'g','LineWidth',1);
xlabel('Time [h]','Interpreter','latex');
ylabel('$\mathbf{q}$ [-]','Interpreter','latex');
grid on; grid minor;

%Propagation SGM CoM dynamics
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~, XX_ref] = ode113(@(t,X) CW(t,X), tvec, x0_ref, options);

%Plot SGN CoM dynamics
figure(3);
plot3(XX_ref(:,1), XX_ref(:,2), XX_ref(:,3),'b','LineWidth',2);
hold on;
plot3(0,0,0,'r.','MarkerSize',15);
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
zlabel('z [m]','Interpreter','latex');
legend('Chaser OSR','Target SGN','Interpreter','latex')
grid on; grid minor;

%Camera measuremnts during acquisition period
meas_vec=zeros(3,8,disc);
visible_vec=zeros(disc,8);
vis_corners=zeros(disc,1);
for i=1:disc
    meas=meas_sim(n,XX_ref(i,1:3)',q_vec(i,:)',tvec(i)-t0,t0,Cam);
    meas_vec(:,1:length(meas.visible),i)=meas.y;
    visible_vec(i,1:length(meas.visible))=meas.visible;
    vis_corners(i)=length(meas.visible);
end

%% Ex.2

%Corners Outside FOV count
Count=0;
for i=1:length(tvec)
    for j=1:8
        if meas_vec(1,j,i)<0 || meas_vec(2,j,i)<0
        Count=Count+1;
        end
        if meas_vec(1,j,i)>Sensor_size(1) || meas_vec(2,j,i)>Sensor_size(2)
        Count=Count+1;
        end
    end
end

%Plot: corners inside FoV during acquisition period
figure(1);
for i=1:8
    scatter(squeeze(meas_vec(1,i,:)), squeeze(meas_vec(2,i,:)),5,'b','.');
    hold on;
end
grid on;
grid minor;
box on;
xlim([0,1920]);
ylim([0,1200]);
xlabel('$q_u$ [pix]','Interpreter','latex');
ylabel('$q_v$ [pix]','Interpreter','latex');

%Plot: number of visible corners during acquisition period
figure(2);
plot((tvec-t0)/60^2,vis_corners,'r.','LineWidth',2);
xlabel('Time [s]','Interpreter','latex');
ylabel('$n_{corners}$ [-]','Interpreter','latex');
grid on; grid minor;

%% Ex.3

%Sequential filtering: EKF
tic;
[XX_ekf, Pekf]=EKF(x0_mean, P0, q_vec, meas_vec, visible_vec, tvec, Cam);
dT_EKF=toc;
XX_ekf=XX_ekf';

%Sequential filtering: UKF
tic;
[XX_ukf, Pukf]=UKF(x0_mean, P0, q_vec, meas_vec, visible_vec, tvec, Cam);
dT_UKF=toc;
XX_ukf=XX_ukf';

%% Ex.4

%Plot results EKF
figure(1);
plot((tvec-t0)/60^2,XX_ekf(:,1),'b','LineWidth',1);
hold on;
plot((tvec-t0)/60^2,XX_ekf(:,2),'r','LineWidth',1);
plot((tvec-t0)/60^2,XX_ekf(:,3),'k','LineWidth',1);
xlabel('Time [s]','Interpreter','latex');
ylabel('$\mathbf{r}$ [m]','Interpreter','latex');
legend('$x$','$y$','$z$','Interpreter','latex','Location','Southeast');
grid on; grid minor;
figure(2);
plot((tvec-t0)/60^2,XX_ekf(:,4),'b','LineWidth',1);
hold on;
plot((tvec-t0)/60^2,XX_ekf(:,5),'r','LineWidth',1);
plot((tvec-t0)/60^2,XX_ekf(:,6),'k','LineWidth',1);
xlabel('Time [s]','Interpreter','latex');
ylabel('$\mathbf{v}$ [m/s]','Interpreter','latex');
legend('$v_x$','$v_y$','$v_z$','Interpreter','latex',...
    'Location','Southeast');
grid on; grid minor;

%Plot results UKF
figure(3);
plot((tvec-t0)/60^2,XX_ukf(:,1),'b','LineWidth',1);
hold on;
plot((tvec-t0)/60^2,XX_ukf(:,2),'r','LineWidth',1);
plot((tvec-t0)/60^2,XX_ukf(:,3),'k','LineWidth',1);
xlabel('Time [s]','Interpreter','latex');
ylabel('$\mathbf{r}$ [m]','Interpreter','latex');
legend('$x$','$y$','$z$','Interpreter','latex','Location','Southeast');
grid on; grid minor;
figure(4);
plot((tvec-t0)/60^2,XX_ukf(:,4),'b','LineWidth',1);
hold on;
plot((tvec-t0)/60^2,XX_ukf(:,5),'r','LineWidth',1);
plot((tvec-t0)/60^2,XX_ukf(:,6),'k','LineWidth',1);
xlabel('Time [s]','Interpreter','latex');
ylabel('$\mathbf{v}$ [m/s]','Interpreter','latex');
legend('$v_x$','$v_y$','$v_z$','Interpreter','latex',...
    'Location','Southeast');
grid on; grid minor;

%Analysis error EKF
err_ekf_r=zeros(disc,1);
sigma_ekf_r=zeros(disc,1);
err_ekf_v=zeros(disc,1);
sigma_ekf_v=zeros(disc,1);
for i=1:disc
    err_ekf_r(i)=norm(XX_ekf(i,1:3)-XX_ref(i,1:3));
    sigma_ekf_r(i)=sqrt(trace(Pekf(1:3,1:3,i)));
    err_ekf_v(i)=norm(XX_ekf(i,4:6)-XX_ref(i,4:6));
    sigma_ekf_v(i)=sqrt(trace(Pekf(4:6,4:6,i)));
end
figure(5);
semilogy((tvec-t0)/60^2,err_ekf_r,'b','LineWidth',1);
hold on;
semilogy((tvec-t0)/60^2,3*sigma_ekf_r,'k','LineWidth',1);
xlim([0,24]);
xlabel('Time [s]','Interpreter','latex');
ylabel('$\varepsilon_{r}$ [-]','Interpreter','latex');
legend('$\varepsilon_{r}$','$3\sigma_{r}$','Interpreter','latex');
grid on; grid minor;
figure(6);
semilogy((tvec-t0)/60^2,err_ekf_v,'b','LineWidth',1);
hold on;
semilogy((tvec-t0)/60^2,3*sigma_ekf_v,'k','LineWidth',1);
xlim([0,24]);
xlabel('Time [s]','Interpreter','latex');
ylabel('$\varepsilon_{v}$ [-]','Interpreter','latex');
legend('$\varepsilon_{v}$','$3\sigma_{v}$','Interpreter','latex');
grid on; grid minor;

%Analysis error UKF
err_ukf_r=zeros(disc,1);
sigma_ukf_r=zeros(disc,1);
err_ukf_v=zeros(disc,1);
sigma_ukf_v=zeros(disc,1);
for i=1:disc
    err_ukf_r(i)=norm(XX_ukf(i,1:3)-XX_ref(i,1:3));
    sigma_ukf_r(i)=sqrt(trace(Pukf(1:3,1:3,i)));
    err_ukf_v(i)=norm(XX_ukf(i,4:6)-XX_ref(i,4:6));
    sigma_ukf_v(i)=sqrt(trace(Pukf(4:6,4:6,i)));
end
figure(7);
semilogy((tvec-t0)/60^2,err_ukf_r,'r','LineWidth',1);
hold on;
semilogy((tvec-t0)/60^2,3*sigma_ukf_r,'k','LineWidth',1);
xlim([0,24]);
xlabel('Time [s]','Interpreter','latex');
ylabel('$\varepsilon_{r}$ [-]','Interpreter','latex');
legend('$\varepsilon_{r}$','$3\sigma_{r}$','Interpreter','latex');
grid on; grid minor;
figure(8);
semilogy((tvec-t0)/60^2,err_ukf_v,'r','LineWidth',1);
hold on;
semilogy((tvec-t0)/60^2,3*sigma_ukf_v,'k','LineWidth',1);
xlim([0,24]);
xlabel('Time [s]','Interpreter','latex');
ylabel('$\varepsilon_{v}$ [-]','Interpreter','latex');
legend('$\varepsilon_{v}$','$3\sigma_{v}$','Interpreter','latex');
grid on; grid minor;

%% Functions

function dXdt=dynamicsSGN(~,X)
%Initialization
dXdt=zeros(7,1);

%Data SGN
l=10;
h=5;
d=3;
m=213000;
J=(m/12)*diag([d^2+h^2,l^2+h^2,l^2+d^2]);

%Unpacking state
wx=X(1);
wy=X(2);
wz=X(3);
q1=X(4);
q2=X(5);
q3=X(6);
q4=X(7);

%Composing dXdt
K=1;
eps=1-(q1^2+q2^2+q3^2+q4^2);
Omega=[0, -wx, -wy, -wz;...
    wx, 0, wz, -wy;...
    wy, -wz, 0, wx;...
    wz, wy, -wx, 0];
dXdt(1:3)=-J\(cross([wx;wy;wz],J*[wx;wy;wz]));
dXdt(4:7)=0.5*Omega*[q1; q2; q3; q4]+K*eps*[q1; q2; q3; q4];
end

function dXdt=CW(~,X)
%Data GEO orbit
mu=398600;
R=42241.08;
n=sqrt(mu/R^3);

%Initialization
dXdt=zeros(6,1);

%Unpacking state
x=X(1);
% y=X(2);
z=X(3);
xdot=X(4);
ydot=X(5);
zdot=X(6);

%Composing dXdt
dXdt(1:3)=[xdot; ydot; zdot];
dXdt(4:6)=[3*n^2*x+2*n*ydot;...
    -2*n*xdot;...
    -n^2*z];
end

function dxMdt=STMexact_CW(~,xM)
%Data GEO orbit
mu=398600;
R=42241.08;
n=sqrt(mu/R^3);

%Initialize
dxMdt=zeros(42,1);

%Unpacking state
x=xM(1);
% y=xM(2);
z=xM(3);
xdot=xM(4);
ydot=xM(5);
zdot=xM(6);
M=(reshape(xM(7:42),6,6))';   %From equations to STM

%CW dynamics
dxMdt(1:3)=[xdot; ydot; zdot];
dxMdt(4:6)=[3*n^2*x+2*n*ydot;...
    -2*n*xdot;...
    -n^2*z];

% Variational equations: STM dynamics
A1=zeros(3);
A2=eye(3);
A3=[3*n^2, 0, 0;...
    0, 0, 0;...
    0, 0, -n^2];
A4=[0, 2*n, 0;...
    -2*n, 0, 0;...
    0, 0, 0];
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

function [meas_vert,H]=meas_est(n,X_rel, q,t,~,Cam)
%Camera data
C=Cam.Cframe;
f=Cam.f;
d=Cam.d;
b=Cam.b;
p0=Cam.p0;

%Vertices position
l_ver=10;
h_ver=5;
d_ver=3;
Vert=0.5*[l_ver, l_ver, l_ver, l_ver, -l_ver, -l_ver, -l_ver, -l_ver;...
    -d_ver, d_ver, d_ver, -d_ver, -d_ver, d_ver, d_ver, -d_ver;...
    -h_ver, -h_ver, h_ver, h_ver, -h_ver, -h_ver, h_ver, h_ver];

%Aln
Aln=[cos(n*t), sin(n*t), 0;...
    -sin(n*t), cos(n*t),0;...
    0, 0, 1];

%Abn
Abn=quat2dcm(q');

%Abl
Abl=Abn/Aln;

%Vertices in Cam.frame wrt Chaser
X_vert=zeros(3,8);
for j=1:8
    X_vert(:,j)=-C*X_rel(1:3)+C*(Abl\Vert(:,j));
end

%Measurements composition
meas_eq=@(x,y,z) [p0(1)-d*f*y/z; p0(2)+d*f*x/z; b*d*f/z];
meas_vert=zeros(3,8);
for j=1:8
    meas_vert(:,j)=meas_eq(X_vert(1,j),X_vert(2,j),X_vert(3,j));
end

%Jacobian composition
if nargout>1
    H=zeros(3,6,8);
    H1=@(x,y,z) [0, -d*f/z, d*f*y/z^2;...
        d*f/z, 0, -d*f*x/z^2;...
        0, 0, -b*d*f/z^2];
    for j=1:8
        H(:,:,j)=[-H1(X_vert(1,j),X_vert(2,j),X_vert(3,j))*C, zeros(3,3)];
    end
end
end

function [X, P]=EKF(X0, P0, q, meas, vis, tvec, Cam)   
%Initialization
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
disc=length(tvec);
X=zeros(6,disc);
P=zeros(6,6,disc);
X(:,1)=X0;
P(:,:,1)=P0;
t0_str='UTC 2023/04/01-14:55:12.023';
t0=cspice_str2et(t0_str);

%Camera measurements noise
Rsing=Cam.R;

%Data GEO orbit
mu=398600;
Re=42241.08;
n=sqrt(mu/Re^3);

for i=2:disc
    %Old estimated states
    X_old=X(:,i-1);
    P_old=squeeze(P(:,:,i-1));
    
    %Prediction
    [~, xxM] = ode113(@(t,xM) STMexact_CW(t,xM),...
    [tvec(i-1), tvec(i)],[X_old; reshape(eye(6),[],1)], options);
    
    X_min=xxM(end,1:6);
    STM=(reshape(xxM(end,7:42),6,6))';
    P_min=STM*P_old*STM';
    
    %Real measurements
    Vis_corners=nonzeros(vis(i,:));
    Yreal=nonzeros(meas(:,:,i));
    
    %Measurement equation
    [meas_vert,H_full]=meas_est(n, X_min(1:3)', q(i,:)',...
        tvec(i)-t0, tvec(1), Cam);
    Ymeas=nonzeros(meas_vert(:,Vis_corners));
    H=zeros(3*length(Vis_corners),6);
    for k=1:length(Vis_corners)
        H((3*k-2):3*k,:)=H_full(:,:,Vis_corners(k));
    end
    
    %Kalman gain
    Rcell=repmat({Rsing},1,length(Vis_corners));
    R=blkdiag(Rcell{:});
    K=P_min*H'/(H*P_min*H'+R);
    
    %Correction
    X_plus=X_min'+K*(Yreal-Ymeas);
    P_plus=(eye(6)-K*H)*P_min;
    
    %Estimated state and covariance
    X(:,i)=X_plus;
    P(:,:,i)=P_plus;
end
end

function [X, P]=UKF(X0, P0, q, meas, vis, tvec, Cam)   
%Initialization
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
disc=length(tvec);
X=zeros(6,disc);
P=zeros(6,6,disc);
X(:,1)=X0;
P(:,:,1)=P0;
t0_str='UTC 2023/04/01-14:55:12.023';
t0=cspice_str2et(t0_str);

%Initialization UT algorithm
nUT=6;
alpha=1e-3;
k=0;
lambda=alpha^2*(nUT+k)-nUT;
beta=2;

%UT weigths
w0m=lambda/(nUT+lambda);
w0p=lambda/(nUT+lambda)+(1-alpha^2+beta);
wim=1/(2*(nUT+lambda));
wip=1/(2*(nUT+lambda));

%Camera measurement noise
Rsing=Cam.R;

%Data GEO orbit
mu=398600;
Re=42241.08;
n=sqrt(mu/Re^3);

for i=2:disc
    %Old estimated states
    X_old=X(:,i-1);
    P_old=squeeze(P(:,:,i-1));

    %Sigma Initial Points
    UTmatrix=sqrtm((nUT+lambda)*P_old);
    chi0=X_old;
    chi1=X_old+UTmatrix(:,1);
    chi2=X_old+UTmatrix(:,2);
    chi3=X_old+UTmatrix(:,3);
    chi4=X_old+UTmatrix(:,4);
    chi5=X_old+UTmatrix(:,5);
    chi6=X_old+UTmatrix(:,6);
    chi7=X_old-UTmatrix(:,1);
    chi8=X_old-UTmatrix(:,2);
    chi9=X_old-UTmatrix(:,3);
    chi10=X_old-UTmatrix(:,4);
    chi11=X_old-UTmatrix(:,5);
    chi12=X_old-UTmatrix(:,6);

    %Propagation Sigma Points
    %chi0
    [~, XX] = ode113(@(t,xx) CW(t,xx),[tvec(i-1), tvec(i)],chi0, options);
    jota0=XX(end,:);  
    %chi1
    [~, XX] = ode113(@(t,xx) CW(t,xx),[tvec(i-1), tvec(i)],chi1, options);
    jota1=XX(end,:);   
    %chi2
    [~, XX] = ode113(@(t,xx) CW(t,xx),[tvec(i-1), tvec(i)],chi2, options);
    jota2=XX(end,:); 
    %chi3
    [~, XX] = ode113(@(t,xx) CW(t,xx),[tvec(i-1), tvec(i)],chi3, options);
    jota3=XX(end,:); 
    %chi4
    [~, XX] = ode113(@(t,xx) CW(t,xx),[tvec(i-1), tvec(i)],chi4, options);
    jota4=XX(end,:);  
    %chi5
    [~, XX] = ode113(@(t,xx) CW(t,xx),[tvec(i-1), tvec(i)],chi5, options);
    jota5=XX(end,:); 
    %chi6
    [~, XX] = ode113(@(t,xx) CW(t,xx),[tvec(i-1), tvec(i)],chi6, options);
    jota6=XX(end,:); 
    %chi6
    [~, XX] = ode113(@(t,xx) CW(t,xx),[tvec(i-1), tvec(i)],chi7, options);
    jota7=XX(end,:); 
    %chi6
    [~, XX] = ode113(@(t,xx) CW(t,xx),[tvec(i-1), tvec(i)],chi8, options);
    jota8=XX(end,:); 
    %chi6
    [~, XX] = ode113(@(t,xx) CW(t,xx),[tvec(i-1), tvec(i)],chi9, options);
    jota9=XX(end,:); 
    %chi6
    [~, XX] = ode113(@(t,xx) CW(t,xx),[tvec(i-1), tvec(i)],chi10, options);
    jota10=XX(end,:); 
    %chi6
    [~, XX] = ode113(@(t,xx) CW(t,xx),[tvec(i-1), tvec(i)],chi11, options);
    jota11=XX(end,:); 
    %chi6
    [~, XX] = ode113(@(t,xx) CW(t,xx),[tvec(i-1), tvec(i)],chi12, options);
    jota12=XX(end,:); 

    %Real measurements
    Vis_corners=nonzeros(vis(i,:));
    Yreal=nonzeros(meas(:,:,i));
    
    %Measurement equations ON Sigma Points
    %chi0
    meas_vert=meas_est(n, jota0(1:3)', q(i,:)', tvec(i)-t0, tvec(1), Cam);
    Ymeas0=nonzeros(meas_vert(:,Vis_corners));
    %chi1
    meas_vert=meas_est(n, jota1(1:3)', q(i,:)', tvec(i)-t0, tvec(1), Cam);
    Ymeas1=nonzeros(meas_vert(:,Vis_corners));
    %chi2
    meas_vert=meas_est(n, jota2(1:3)', q(i,:)', tvec(i)-t0, tvec(1), Cam);
    Ymeas2=nonzeros(meas_vert(:,Vis_corners));
    %chi3
    meas_vert=meas_est(n, jota3(1:3)', q(i,:)', tvec(i)-t0, tvec(1), Cam);
    Ymeas3=nonzeros(meas_vert(:,Vis_corners));
    %chi4
    meas_vert=meas_est(n, jota4(1:3)', q(i,:)', tvec(i)-t0, tvec(1), Cam);
    Ymeas4=nonzeros(meas_vert(:,Vis_corners));
    %chi5
    meas_vert=meas_est(n, jota5(1:3)', q(i,:)', tvec(i)-t0, tvec(1), Cam);
    Ymeas5=nonzeros(meas_vert(:,Vis_corners));
    %chi6
    meas_vert=meas_est(n, jota6(1:3)', q(i,:)', tvec(i)-t0, tvec(1), Cam);
    Ymeas6=nonzeros(meas_vert(:,Vis_corners));
    %chi7
    meas_vert=meas_est(n, jota7(1:3)', q(i,:)', tvec(i)-t0, tvec(1), Cam);
    Ymeas7=nonzeros(meas_vert(:,Vis_corners));
    %chi8
    meas_vert=meas_est(n, jota8(1:3)', q(i,:)', tvec(i)-t0, tvec(1), Cam);
    Ymeas8=nonzeros(meas_vert(:,Vis_corners));
    %chi9
    meas_vert=meas_est(n, jota9(1:3)', q(i,:)', tvec(i)-t0, tvec(1), Cam);
    Ymeas9=nonzeros(meas_vert(:,Vis_corners));
    %chi10
    meas_vert=meas_est(n, jota10(1:3)', q(i,:)', tvec(i)-t0, tvec(1), Cam);
    Ymeas10=nonzeros(meas_vert(:,Vis_corners));
    %chi11
    meas_vert=meas_est(n, jota11(1:3)', q(i,:)', tvec(i)-t0, tvec(1), Cam);
    Ymeas11=nonzeros(meas_vert(:,Vis_corners));
    %chi12
    meas_vert=meas_est(n, jota12(1:3)', q(i,:)', tvec(i)-t0, tvec(1), Cam);
    Ymeas12=nonzeros(meas_vert(:,Vis_corners));

    %Mean and Covariances MINUS
    %Means
    X_min=w0m*jota0+wim*(jota1+jota2+jota3+...
    jota4+jota5+jota6+jota7+jota8+jota9+...
    jota10+jota11+jota12);

    Y_min=w0m*Ymeas0+wim*(Ymeas1+Ymeas2+Ymeas3+...
    Ymeas4+Ymeas5+Ymeas6+Ymeas7+Ymeas8+Ymeas9+...
    Ymeas10+Ymeas11+Ymeas12);

    %Covariances
    Rcell=repmat({Rsing},1,length(Vis_corners));
    R=blkdiag(Rcell{:});

    P_min=w0p*(jota0-X_min)'*(jota0-X_min)+...
    wip*((jota1-X_min)'*(jota1-X_min)+...
    (jota2-X_min)'*(jota2-X_min)+...
    (jota3-X_min)'*(jota3-X_min)+...
    (jota4-X_min)'*(jota4-X_min)+...
    (jota5-X_min)'*(jota5-X_min)+...
    (jota6-X_min)'*(jota6-X_min)+...
    (jota7-X_min)'*(jota7-X_min)+...
    (jota8-X_min)'*(jota8-X_min)+...
    (jota9-X_min)'*(jota9-X_min)+...
    (jota10-X_min)'*(jota10-X_min)+...
    (jota11-X_min)'*(jota11-X_min)+...
    (jota12-X_min)'*(jota12-X_min));

    P_ee=w0p*(Ymeas0-Y_min)*(Ymeas0-Y_min)'+...
    wip*((Ymeas1-Y_min)*(Ymeas1-Y_min)'+...
    (Ymeas2-Y_min)*(Ymeas2-Y_min)'+...
    (Ymeas3-Y_min)*(Ymeas3-Y_min)'+...
    (Ymeas4-Y_min)*(Ymeas4-Y_min)'+...
    (Ymeas5-Y_min)*(Ymeas5-Y_min)'+...
    (Ymeas6-Y_min)*(Ymeas6-Y_min)'+...
    (Ymeas7-Y_min)*(Ymeas7-Y_min)'+...
    (Ymeas8-Y_min)*(Ymeas8-Y_min)'+...
    (Ymeas9-Y_min)*(Ymeas9-Y_min)'+...
    (Ymeas10-Y_min)*(Ymeas10-Y_min)'+...
    (Ymeas11-Y_min)*(Ymeas11-Y_min)'+...
    (Ymeas12-Y_min)*(Ymeas12-Y_min)')+R;

    P_xy=w0p*(jota0-X_min)'*(Ymeas0-Y_min)'+...
    wip*((jota1-X_min)'*(Ymeas1-Y_min)'+...
    (jota2-X_min)'*(Ymeas2-Y_min)'+...
    (jota3-X_min)'*(Ymeas3-Y_min)'+...
    (jota4-X_min)'*(Ymeas4-Y_min)'+...
    (jota5-X_min)'*(Ymeas5-Y_min)'+...
    (jota6-X_min)'*(Ymeas6-Y_min)'+...
    (jota7-X_min)'*(Ymeas7-Y_min)'+...
    (jota8-X_min)'*(Ymeas8-Y_min)'+...
    (jota9-X_min)'*(Ymeas9-Y_min)'+...
    (jota10-X_min)'*(Ymeas10-Y_min)'+...
    (jota11-X_min)'*(Ymeas11-Y_min)'+...
    (jota12-X_min)'*(Ymeas12-Y_min)');

    %Kalman gain
    K=P_xy/P_ee;
    
    %Correction
    X_plus=X_min'+K*(Yreal-Y_min);
    P_plus=P_min-K*P_ee*K';

    %Estimated state and covariance
    X(:,i)=X_plus;
    P(:,:,i)=P_plus;
end
end



