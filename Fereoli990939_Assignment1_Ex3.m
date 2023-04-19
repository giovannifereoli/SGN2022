%% Spacecraft Guidance and Navigation - Continuos Guidance (2022/2023)
% Assignment:     1
% Exercise:       Continuos Guidance
% Author:         Giovanni Fereoli
% ID number:      990939

%NB: Initial costate and final time guesses generation algorithm 
%    is commented. This process can require ~1/5 min (it's a matter of
%    luck), therefore for code usability reasons a solution is already
%    reported. To verify the procedure remove '%' and run.

%% Ex.2
clear; clc; close;

%Load SPICE Kernels
cspice_furnsh('sgn_Assignment1_Ex3.tm');

%Dimension units
DU=149597870700*1e-3; 
TU=5.0229*1e6;  
MU=1;  

%Data
ti = cspice_str2et('2022-08-03-12:45:20.000 UTC');
Xi_planet1=cspice_spkezr('Earth', ti, 'ECLIPJ2000', 'NONE', 'SSB');
Xi_planet1(1:3)=Xi_planet1(1:3)/DU;
Xi_planet1(4:6)=Xi_planet1(4:6)/DU*TU;

Data.mu=cspice_bodvrd('Sun','GM',1)/(DU^3/TU^2);
Data.m0=1500/MU;
Data.T=0.15*10^-3/(MU*DU/TU^2);
Data.N=4;
Data.Isp=3000/TU;
Data.ti=ti/TU;
Data.Earth=Xi_planet1;
Data.g0=(9.81*1e-3)/(DU/TU^2);

% % Initial Costate Generation and Shooting Function Resolution
% err=10;
% while err>1e-6
%     %Guess generation
%     tf_guess=Data.ti+((90+100*rand(1))*86400)/TU;
%     lambda0_guess=zeros(7,1);
%     lambda0_guess(1)=10^(6*(2*rand(1)-1))*(2*rand(1)-1);
%     lambda0_guess(2)=10^(6*(2*rand(1)-1))*(2*rand(1)-1);
%     lambda0_guess(3)=10^(6*(2*rand(1)-1))*(2*rand(1)-1);
%     lambda0_guess(4)=10^(6*(2*rand(1)-1))*(2*rand(1)-1);
%     lambda0_guess(5)=10^(6*(2*rand(1)-1))*(2*rand(1)-1);
%     lambda0_guess(6)=10^(6*(2*rand(1)-1))*(2*rand(1)-1);
%     lambda0_guess(7)=10^(6*(2*rand(1)-1));
%     
%     %Shooting Function Roots Search
%     options_fsolve=optimset('Algorithm','levenberg-marquardt','Display',...
%         'iter','MaxIter', 100,'MaxFunEvals',1000); 
%     [Y_opt, fval, exitflag, output]=fsolve(@(Y) const_OCP(Y, Data),...
%         [lambda0_guess; tf_guess], options_fsolve);
%     lambda0_opt=Y_opt(1:7);
%     tf_opt=Y_opt(8);
%     
%     %Exit Criteria
%     Res=const_OCP(Y_opt, Data);
%     err=norm(Res);
% end

%Shooting Function solutions (already provided)
lambda0_opt=[-1.542836277927404e+02;1.093792387924383e+02;...
    4.836122581483903;-1.903491311643344e+02;-76.375740676180540;...
    -0.147450732027128;0.014904763393519];
tf_opt=1.472024733839468e+02;

%Error analysis
Res=const_OCP([lambda0_opt; tf_opt], Data);
errPos1=norm(Res(1:3))*DU;
errVel1=norm(Res(4:6))*DU/TU*1e3;

%Optimized trajectory integration
X0=[Xi_planet1; Data.m0; lambda0_opt];
options_ode=odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[tt, YY]=ode113(@(t, Y) RTBP_OCP(t, Y, Data), ...
    linspace(Data.ti, tf_opt, 200), X0, options_ode);

%Earth/Mars positions and trajectories
X0_planet1=cspice_spkezr('Earth', Data.ti*TU, 'ECLIPJ2000', 'NONE', 'SSB');
Xf_planet2=cspice_spkezr('Mars', tf_opt*TU, 'ECLIPJ2000', 'NONE', 'SSB');
rt_planet1=cspice_spkpos('Earth', linspace(Data.ti*TU, tf_opt*TU, 1000),...
    'ECLIPJ2000', 'NONE', 'SSB'); 
rt_planet2=cspice_spkpos('Mars', linspace(Data.ti*TU, tf_opt*TU, 1000),...
    'ECLIPJ2000', 'NONE', 'SSB'); 

%Post-analysis
Isp=3000;
g0=9.81;
m0=1500;
ToF1=(tf_opt-Data.ti)*TU/86400;        %Days
tf_data1=cspice_et2utc(tf_opt*TU,'C', 0);
mf=YY(end,7)*MU;
dV1=(Isp*g0*log(m0/mf))/1000;          %km/s
Ht=zeros(length(YY),1);
primevec=zeros(length(YY),3);
alpha=zeros(length(YY),1);
beta=zeros(length(YY),1);
rho=zeros(length(YY),1);
for i=1:length(YY)
 primevec_i=-YY(i,11:13)./norm(YY(i,11:13));
 primevec(i,:)=primevec_i;
 [alpha(i), beta(i), rho(i)]=cart2sph(primevec_i(1), primevec_i(2),...
     primevec_i(3));
 alpha(i)=rad2deg(alpha(i));
 beta(i)=rad2deg(beta(i));
Ht(i)=Hamiltonian(YY(i,:), Data);
end

%Plot trajectory
figure(1)
plot3(YY(:,1), YY(:,2), YY(:,3),'b',...
    'LineWidth',2);
xlabel('x [AU]');
ylabel('y [AU]');
zlabel('y [AU]');
grid on;
real axis;
hold on;
plot3(rt_planet1(1,:)/DU, rt_planet1(2,:)/DU, rt_planet1(3,:)/DU,'--',...
    'LineWidth',2);
plot3(rt_planet2(1,:)/DU, rt_planet2(2,:)/DU, rt_planet2(3,:)/DU,'--',...
    'LineWidth',2);
quiver3(YY(:,1), YY(:,2), YY(:,3), primevec(:,1), primevec(:,2),...
    primevec(:,3));
plot3(X0_planet1(1)/DU, X0_planet1(2)/DU, X0_planet1(3)/DU,...
    'k.','MarkerSize',10);
plot3(Xf_planet2(1)/DU, Xf_planet2(2)/DU, Xf_planet2(3)/DU,...
    'k.','MarkerSize',10);
plot3(0, 0, 0,'k.','MarkerSize',10);
text(X0_planet1(1)/DU, X0_planet1(2)/DU, X0_planet1(3)/DU,'EARTH',...
    'HorizontalAlignment', 'center','VerticalAlignment','top');
text(Xf_planet2(1)/DU, Xf_planet2(2)/DU, Xf_planet2(3)/DU,'MARS',...
    'HorizontalAlignment', 'center','VerticalAlignment','top');
text(0, 0, 0,'SUN','HorizontalAlignment', 'center','VerticalAlignment',...
    'top');
axis equal;
view(2);

%Plot Thrusting angles
figure(2);
plot((tt-Data.ti)*TU/86400, alpha,'r','LineWidth',2);
hold on;
plot((tt-Data.ti)*TU/86400, beta,'b','LineWidth',2);
grid on;
legend('In-Plane','Out-of-Plane', 'Location', 'southwest');
xlabel('Time [Days]');
ylabel('Angle [Deg]');

%Plot Hamiltonian
figure(3);
plot((tt-Data.ti)*TU/86400,Ht,'LineWidth',2);
xlabel('Time [Days]');
ylabel('H [-]');
axis equal;
grid on;

%Clear Kernels
cspice_kclear();

%% Ex.3
clear; clc; close;

%Load SPICE Kernels
cspice_furnsh('sgn_Assignment1_Ex3.tm');

%Dimension units
DU=149597870700*1e-3; 
TU=5.0229*1e6;  
MU=1;  

%Data
ti = cspice_str2et('2022-08-03-12:45:20.000 UTC');
Xi_planet1=cspice_spkezr('Earth', ti, 'ECLIPJ2000', 'NONE', 'SSB');
Xi_planet1(1:3)=Xi_planet1(1:3)/DU;
Xi_planet1(4:6)=Xi_planet1(4:6)/DU*TU;

Data.mu=cspice_bodvrd('Sun','GM',1)/(DU^3/TU^2);
Data.m0=1500/MU;
Data.T=0.15*10^-3/(MU*DU/TU^2);
Data.N=3;
Data.Isp=3000/TU;
Data.ti=ti/TU;
Data.Earth=Xi_planet1;
Data.g0=(9.81*1e-3)/(DU/TU^2);

% % Initial Costate Generation and Shooting Function resolution
% err=10;
% while err>1e-6
%     %Guess generation
%     tf_guess=Data.ti+((300+300*rand(1))*86400)/TU;
%     lambda0_guess=zeros(7,1);
%     lambda0_guess(1)=10^(6*(2*rand(1)-1))*(2*rand(1)-1);
%     lambda0_guess(2)=10^(6*(2*rand(1)-1))*(2*rand(1)-1);
%     lambda0_guess(3)=10^(6*(2*rand(1)-1))*(2*rand(1)-1);
%     lambda0_guess(4)=10^(6*(2*rand(1)-1))*(2*rand(1)-1);
%     lambda0_guess(5)=10^(6*(2*rand(1)-1))*(2*rand(1)-1);
%     lambda0_guess(6)=10^(6*(2*rand(1)-1))*(2*rand(1)-1);
%     lambda0_guess(7)=10^(6*(2*rand(1)-1));
%     
%     %Shooting Function Roots Search
%     options_fsolve=optimset('Algorithm','levenberg-marquardt','Display',...
%         'iter','MaxIter',10000,'MaxFunEvals',10000); 
%     [Y_opt, fval, exitflag, output]=fsolve(@(Y) const_OCP(Y, Data),...
%         [lambda0_guess; tf_guess], options_fsolve);
%     lambda0_opt=Y_opt(1:7);
%     tf_opt=Y_opt(8);
%     
%     %Exit Criteria
%     Res=const_OCP(Y_opt, Data);
%     err=norm(Res);
% end

%Shooting Function solutions (already provided)
lambda0_opt=[-5.856827654680059;5.229510677613138;0.052751724221299;...
    -7.711925502000857;-4.394719988347892;0.057426491200811;...
    6.963315927870564e-04];
tf_opt=1.501023792996373e+02;

%Error analysis
Res=const_OCP([lambda0_opt; tf_opt], Data);
errPos2=norm(Res(1:3))*DU;
errVel2=norm(Res(4:6))*DU/TU*1e3;

%Optimized trajectory integration
X0=[Xi_planet1; Data.m0; lambda0_opt];
options_ode=odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[tt, YY]=ode113(@(t, Y) RTBP_OCP(t, Y, Data), ...
    linspace(Data.ti, tf_opt, 200), X0, options_ode);

%Earth/Mars positions and trajectories
X0_planet1=cspice_spkezr('Earth', Data.ti*TU, 'ECLIPJ2000', 'NONE', 'SSB');
Xf_planet2=cspice_spkezr('Mars', tf_opt*TU, 'ECLIPJ2000', 'NONE', 'SSB');
rt_planet1=cspice_spkpos('Earth', linspace(Data.ti*TU, tf_opt*TU, 1000),...
    'ECLIPJ2000', 'NONE', 'SSB'); 
rt_planet2=cspice_spkpos('Mars', linspace(Data.ti*TU, tf_opt*TU, 1000),...
    'ECLIPJ2000', 'NONE', 'SSB'); 

%Post-analysis
Isp=3000;
g0=9.81;
m0=1500;
ToF2=(tf_opt-Data.ti)*TU/86400;        %Days
tf_data2=cspice_et2utc(tf_opt*TU,'C', 0);
mf=YY(end,7)*MU;
dV2=(Isp*g0*log(m0/mf))/1000;          %km/s
Ht=zeros(length(YY),1);
primevec=zeros(length(YY),3);
alpha=zeros(length(YY),1);
beta=zeros(length(YY),1);
rho=zeros(length(YY),1);
for i=1:length(YY)
 primevec_i=-YY(i,11:13)./norm(YY(i,11:13));
 primevec(i,:)=primevec_i;
 [alpha(i), beta(i), rho(i)]=cart2sph(primevec_i(1), primevec_i(2),...
     primevec_i(3));
 alpha(i)=rad2deg(alpha(i));
 beta(i)=rad2deg(beta(i));
Ht(i)=Hamiltonian(YY(i,:), Data);
end

%Plot trajectory
figure(1)
plot3(YY(:,1), YY(:,2), YY(:,3),'b',...
    'LineWidth',2);
xlabel('x [AU]');
ylabel('y [AU]');
zlabel('y [AU]');
grid on;
real axis;
hold on;
plot3(rt_planet1(1,:)/DU, rt_planet1(2,:)/DU, rt_planet1(3,:)/DU,'--',...
    'LineWidth',2);
plot3(rt_planet2(1,:)/DU, rt_planet2(2,:)/DU, rt_planet2(3,:)/DU,'--',...
    'LineWidth',2);
quiver3(YY(:,1), YY(:,2), YY(:,3), primevec(:,1), primevec(:,2),...
    primevec(:,3));
plot3(X0_planet1(1)/DU, X0_planet1(2)/DU, X0_planet1(3)/DU,...
    'k.','MarkerSize',10);
plot3(Xf_planet2(1)/DU, Xf_planet2(2)/DU, Xf_planet2(3)/DU,...
    'k.','MarkerSize',10);
plot3(0, 0, 0,'k.','MarkerSize',10);
text(X0_planet1(1)/DU, X0_planet1(2)/DU, X0_planet1(3)/DU,'EARTH',...
    'HorizontalAlignment', 'center','VerticalAlignment','top');
text(Xf_planet2(1)/DU, Xf_planet2(2)/DU, Xf_planet2(3)/DU,'MARS',...
    'HorizontalAlignment', 'center','VerticalAlignment','top');
text(0, 0, 0,'SUN','HorizontalAlignment', 'center','VerticalAlignment',...
    'top');
axis equal;
view(2);

%Plot Thrusting angles
figure(2);
plot((tt-Data.ti)*TU/86400, alpha,'r','LineWidth',2);
hold on;
plot((tt-Data.ti)*TU/86400, beta,'b','LineWidth',2);
grid on;
legend('In-Plane','Out-of-Plane', 'Location', 'southwest');
xlabel('Time [Days]');
ylabel('Angle [Deg]');

%Plot Hamiltonian
figure(3);
plot((tt-Data.ti)*TU/86400,Ht,'LineWidth',2);
xlabel('Time [Days]');
ylabel('H [-]');
axis equal;
grid on;

%Clear Kernels
cspice_kclear();

%% Functions 

function H=Hamiltonian(Y, Data)
%Data
Tmax=Data.T*Data.N;
Isp=Data.Isp;
g0=Data.g0;
mu=Data.mu;
u_star=1;

%Unpacking
r=Y(1:3); r_norm=norm(r);
v=Y(4:6);
m=Y(7);
lambda_r=Y(8:10);
lambda_v=Y(11:13); lambda_vnorm=norm(lambda_v);
lambda_m=Y(14);

%Composing
H=1+dot(lambda_r,v)-(mu/r_norm^3)*dot(lambda_v,r)+(Tmax/(Isp*g0))*...
    u_star*(-lambda_vnorm*Isp*g0/m-lambda_m);
end

function C=const_OCP(Y, Data) 
%Data
mu=Data.mu;
m0=Data.m0;
Isp=Data.Isp;
ti=Data.ti;
X0=Data.Earth;
g0=Data.g0;
Tmax=Data.T*Data.N;
DU=149597870700*1e-3;
TU=5.0229*1e6; 
u_star=1;

%Unpacking
lambda0=Y(1:7);
tf=Y(8);

%Propagation
Y0=[X0; m0; lambda0];
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~, YY] = ode113(@(t,Y) RTBP_OCP(t, Y, Data), [ti tf], Y0, options);
YYend=YY(end,:)';

%Target State
Xf_planet2=cspice_spkezr('Mars', tf*TU, 'ECLIPJ2000', 'NONE', 'SSB');
Xf_planet2(1:3)=Xf_planet2(1:3)/DU;     rf_planet2=Xf_planet2(1:3);
Xf_planet2(4:6)=Xf_planet2(4:6)/DU*TU;  vf_planet2=Xf_planet2(4:6);
af_planet2=-(mu/norm(rf_planet2)^3)*rf_planet2;

%Useful quantities for C
rf=YYend(1:3); vf=YYend(4:6); mf=YYend(7);
lambda_rf=YYend(8:10); lambda_vf=YYend(11:13); lambda_mf=YYend(14);
Hf=1+dot(lambda_rf,vf)--mu/norm(rf)^3*dot(lambda_vf,rf)...
    -u_star*Tmax*norm(lambda_vf)/mf-lambda_mf*u_star*Tmax/(Isp*g0);

%Composition of C
C=zeros(8,1);
C(1:3)=rf-rf_planet2;
C(4:6)=vf-vf_planet2;
C(7)=lambda_mf;
C(8)=Hf-dot(lambda_rf, vf_planet2)-dot(lambda_vf, af_planet2);
end

function dYdt=RTBP_OCP(~,Y, Data)
%Data
mu=Data.mu;
Tmax=Data.T*Data.N;
Isp=Data.Isp;         
g0=Data.g0;

%Initialize
dYdt=zeros(14,1);

%Unpacking
x=Y(1);
y=Y(2);
z=Y(3);
xdot=Y(4);
ydot=Y(5);
zdot=Y(6);
m=Y(7);
lambdax=Y(8);
lambday=Y(9);
lambdaz=Y(10);
lambdaxdot=Y(11);
lambdaydot=Y(12);
lambdazdot=Y(13);
%lambdam=Y(14);

%Useful quantities
u_star=1;
lambdav_norm=sqrt(lambdaxdot^2+lambdaydot^2+lambdazdot^2);
prime_vec=-[lambdaxdot; lambdaydot; lambdazdot]/lambdav_norm;                
r_norm=sqrt(x^2+y^2+z^2);

%RTBP and OCP dynamics (u_star, prime_vec already embedded)
dYdt(1:3)=[xdot; ydot; zdot];
dYdt(4:6)=-(mu/r_norm^3)*[x; y; z]+(u_star*Tmax/m)*prime_vec;
dYdt(7)=-u_star*Tmax/(Isp*g0);
dYdt(8:10)=-(3*mu/r_norm^5)*dot([lambdaxdot; lambdaydot; lambdazdot],...
    [x; y; z])*[x; y; z]+(mu/r_norm^3)*[lambdaxdot; lambdaydot;...
    lambdazdot];
dYdt(11:13)=-[lambdax; lambday; lambdaz];
dYdt(14)=-u_star*lambdav_norm*Tmax/m^2;
end
