%% Spacecraft Guidance and Navigation - Batch filters (2022/2023)
% Assignment:     2
% Exercise:       Batch filters
% Author:         Giovanni Fereoli
% ID number:      990939

%% Ex.1
clear; clc; close all;

%Initialization
addpath('sgp4');
arcsec2rad = pi / (180*3600);
typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)
cspice_furnsh('assignment02.tm');
mu= cspice_bodvrd('EARTH','GM',3);
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);

%Initial state of EX.1
r0=[6054.30795817484; -3072.03883303992; -133.115352431876];
v0=[4.64750094824087; 9.18608475681236; -0.62056520749034];
x0_real=[r0; v0];
t0_real=cspice_str2et('2022-11-11-19:08:49.824 UTC');

%Initial state for EX.2
t0=cspice_str2et('2022-11-12-04:30:00.000 UTC');
th=cspice_str2et('2022-11-12-04:31:00.000 UTC');
tf=cspice_str2et('2022-11-14-16:30:00.000 UTC');
[~, xx] = ode113(@(t,x) RTBP(t,x,mu),[t0_real,t0],x0_real, options);
x0=xx(end,:);

%Integration of spacecraft motion in ECI
h=th-t0;
tvec = t0:h:tf;
npoints=length(tvec);
[~, xx] = ode113(@(t,x) RTBP(t,x,mu),tvec,x0, options);

%Station definition
stat1='KOUROU';
stat2='PERTH';

%Antenna pointings
[Az1, El1, ~, ~]=antenna(stat1, tvec, xx);
[Az2, El2, ~, ~]=antenna(stat2, tvec, xx);

%Plot Orbit
figure(1);
[X, Y, Z] = sphere;
EarthRadius = cspice_bodvrd('EARTH','RADII',3);
hSurface = surf(X*EarthRadius(1), Y*EarthRadius(1), Z*EarthRadius(1));
set(hSurface,'FaceColor',[0.5 0.5 0.5]);
hold on; 
grid on; grid minor;
axis equal;
plot3(xx(:, 1), xx(:, 2), xx(:, 3),'r','LineWidth',1);
xlabel('x [km]','Interpreter','latex');
ylabel('y [km]','Interpreter','latex');
zlabel('z [km]','Interpreter','latex');

%Visbilities windows of Keplerian Motion                                     
i_visibility1 = El1 > deg2rad(10); 
i_visibility2 = El2 > deg2rad(5);
on_tvec1=(tvec(strfind(double(i_visibility1)',[0,1]))-t0)/60^2;
off_tvec1=(tvec(strfind(double(i_visibility1)',[1,0]))-t0)/60^2;
on_tvec2=(tvec(strfind(double(i_visibility2)',[0,1]))-t0)/60^2;
off_tvec2=(tvec(strfind(double(i_visibility2)',[1,0]))-t0)/60^2;

%Plots
figure(2);
plot(Az1(i_visibility1)*cspice_dpr(),...
    El1(i_visibility1)*cspice_dpr(),'b*');
axis([-180,180,0, 90]);
xlabel('Azimuth [deg]','Interpreter','latex');
ylabel('Elevation [deg]','Interpreter','latex');
grid on; grid minor;

figure(3);
plot(Az2(i_visibility2)*cspice_dpr(),...
    El2(i_visibility2)*cspice_dpr(),'r*');
axis([-180,180,0, 90]);
xlabel('Azimuth [deg]','Interpreter','latex');
ylabel('Elevation [deg]','Interpreter','latex');
grid on; grid minor;

%% Ex.2a                                                                                                                      
%Initialization TLE
longstr1 = '1 87654U 22110B   22316.00967942  .00000002  00000-0  32024-3 0  9990';
longstr2 = '2 87654   3.6309 137.1541 8138191 196.1632  96.6141  1.26411866   834';

%Read TLE
satrec = twoline2rv(longstr1, longstr2, typerun,'e', opsmode, whichconst);

%Epoch TLE
[year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f',...
    [year,mon,day,hr,min,sec]);
sat_epoch_et = cspice_str2et(sat_epoch_str);

%Initialization Nutation Correction parameters
ddpsi = -0.113613*arcsec2rad; %  [rad]
ddeps = -0.007056*arcsec2rad; %  [rad]

%SPG4 propagation
reci= zeros(3, npoints);
veci = zeros(3, npoints);

for i = 1:npoints
    %Propagation
    tsince = (tvec(i) - sat_epoch_et)/60.0; % minutes from TLE epoch
    [~, rteme, vteme] = sgp4(satrec,  tsince);
    %Centuries from 2000-01-01T00:00:00.00 TDT
    ttt = cspice_unitim(tvec(i), 'ET', 'TDT')/cspice_jyear()/100;
    %TEME to ECI
    [reci(:,i), veci(:,i)] = ...
        teme2eci(rteme, vteme, [0.0;0.0;0.0],  ttt, ddpsi, ddeps);
end

%Antenna pointings
[Az1, El1, rho1, ~]=antenna(stat1, tvec, [reci',veci']);
[Az2, El2, rho2, ~]=antenna(stat2, tvec, [reci',veci']);

%Visibility of SGP4 with respect to Keplerian Motion
i_visibility1new = El1 > deg2rad(10);   
i_visibility1=logical(i_visibility1new.*i_visibility1);
tvec1=tvec(i_visibility1);
i_visibility2new = El2 > deg2rad(5);
i_visibility2=logical(i_visibility2new.*i_visibility2);
tvec2=tvec(i_visibility2);

%Measurements visible
rho1_vis=rho1(i_visibility1);
rho2_vis=rho2(i_visibility2);
Az1_vis=Az1(i_visibility1);
Az2_vis=Az2(i_visibility2);
El1_vis=El1(i_visibility1);
El2_vis=El2(i_visibility2);
npoints_vis1=length(rho1_vis);
npoints_vis2=length(rho2_vis);

%Plots angular measurements
figure(4);
plot(Az1_vis*cspice_dpr(),...
    El1_vis*cspice_dpr(),'b*');
axis([-180,180,0, 90]);
xlabel('Azimuth [deg]','Interpreter','latex');
ylabel('Elevation [deg]','Interpreter','latex');
grid on; grid minor;

figure(5);
plot(Az2_vis*cspice_dpr(),...
    El2_vis*cspice_dpr(),'r*');
axis([-180,180,0, 90]);
xlabel('Azimuth [deg]','Interpreter','latex');
ylabel('Elevation [deg]','Interpreter','latex');
grid on; grid minor;

%Plots range measurements
figure(6);
plot((tvec1-t0)/60^2, rho1_vis,'b*');
xlabel('Time [h]','Interpreter','latex');
ylabel('Range [km]','Interpreter','latex');
grid on; grid minor;

figure(7);
plot((tvec2-t0)/60^2, rho2_vis,'r*');
xlabel('Time [h]','Interpreter','latex');
ylabel('Range [km]','Interpreter','latex');
grid on; grid minor;

%% Ex.2b          

%Covariance matrix
sig_rho=0.01;
sig_ang=100*1e-3;
sigma_meas=[sig_rho, sig_ang, sig_ang];

%Componing measurements states
meas1=[rho1_vis, rad2deg(Az1_vis), rad2deg(El1_vis)];
meas2=[rho2_vis, rad2deg(Az2_vis), rad2deg(El2_vis)];

%Adding noise
meas1_noise=zeros(npoints_vis1,3);
meas2_noise=zeros(npoints_vis2,3);
for i=1:npoints_vis1
    meas1_noise(i,:)=mvnrnd(meas1(i,:), sigma_meas);
end
for i=1:npoints_vis2
    meas2_noise(i,:)=mvnrnd(meas2(i,:), sigma_meas);
end

%Plots angular measurements with noise
figure(8);
plot(meas1_noise(:,2),...
    meas1_noise(:,3),'b*');
axis([-180,180,0, 90]);
xlabel('Azimuth [deg]','Interpreter','latex');
ylabel('Elevation [deg]','Interpreter','latex');
grid on; grid minor;

figure(9);
plot(meas2_noise(:,2),...
    meas2_noise(:,3),'r*');
axis([-180,180,0, 90]);
xlabel('Azimuth [deg]','Interpreter','latex');
ylabel('Elevation [deg]','Interpreter','latex');
grid on; grid minor;

%Plots range measurements with noise
figure(10);
plot((tvec1-t0)/60^2, meas1_noise(:,1),'b*');
xlabel('Time [h]','Interpreter','latex');
ylabel('Range [km]','Interpreter','latex');
grid on; grid minor;

figure(11);
plot((tvec2-t0)/60^2, meas2_noise(:,1),'r*');
xlabel('Time [h]','Interpreter','latex');
ylabel('Range [km]','Interpreter','latex');
grid on; grid minor;

%% Ex.3a      

%Initalization
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', ...
    'Display', 'iter','TolFun', 1e-12);

%Reference state
x0_tle=[reci(:,1);veci(:,1)];

%Applyin Batch Filter only with PERTH measurements
[xA,resnormA,resA,exitflagA,~,~,jacA]=lsqnonlin(@(x) Jcost1(x,...
    stat2, meas2_noise, tvec, i_visibility2, sigma_meas), x0_tle,...
    [], [], options);

%Resulting covariance
JacA = full(jacA);
P_lsA = resnormA / (length(resA)-length(x0_tle)) .* inv(JacA'*JacA);
sig_rA=sqrt(trace(P_lsA(1:3,1:3)));
sig_vA=sqrt(trace(P_lsA(4:6,4:6)));

%% Ex.3b

%Initalization
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', ...
    'Display', 'iter','TolFun',1e-12);

%Reference state
x0_tle=[reci(:,1);veci(:,1)];

%Applyin Batch Filter only with PERTH and KOUROU measurements
[xB,resnormB,resB,exitflagB,~,~,jacB]=lsqnonlin(@(x) Jcost2(x, stat1,...
    stat2, meas1_noise, meas2_noise, tvec, i_visibility1, i_visibility2,...
    sigma_meas), x0_tle, [], [], options);

%Resulting covariance
JacB = full(jacB);
P_lsB = resnormB / (length(resB)-length(x0_tle)) .* inv(JacB'*JacB);
sig_rB=sqrt(trace(P_lsB(1:3,1:3)));
sig_vB=sqrt(trace(P_lsB(4:6,4:6)));

%% Ex.3c

%Initalization
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', ...
    'Display', 'iter','TolFun',1e-12);

%Reference state
x0_tle=[reci(:,1);veci(:,1)];

%Applyin Batch Filter only with PERTH and KOUROU measurements
[xC,resnormC,resC,exitflagC,~,~,jacC]=lsqnonlin(@(x) Jcost3(x, stat1,...
    stat2, meas1_noise, meas2_noise, tvec, i_visibility1, i_visibility2,...
    sigma_meas), x0_tle, [], [], options);

%Resulting covariance
JacC = full(jacC);
P_lsC = resnormC / (length(resC)-length(x0_tle)) .* inv(JacC'*JacC);
sig_rC=sqrt(trace(P_lsC(1:3,1:3)));
sig_vC=sqrt(trace(P_lsC(4:6,4:6)));

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

%RTBP dynamics
r_norm=sqrt(x^2+y^2+z^2);
dxdt(1:3)=[xdot;ydot;zdot];
dxdt(4:6)=-mu/r_norm^3*[x;y;z];
end

function dxdt = RTBP_J2(t,xx,mu)
%Initialization
dxdt=zeros(6,1);

%Data
EarthRadius=cspice_bodvrd('EARTH','RADII',3);
Re=EarthRadius(1);
Jtwo=1.082616e-3;

%Unpack state
x=xx(1);
y=xx(2);
z=xx(3);
xdot=xx(4);
ydot=xx(5);
zdot=xx(6);
r_norm = sqrt(x^2+y^2+z^2);

%Jtwo Effect
rotm=cspice_pxform('J2000','ITRF93',t);
pos_ECEF=rotm*[x;y;z];
pos_ECEF_norm=norm(pos_ECEF);
a_Jtwo=3/2*mu*Jtwo*pos_ECEF/pos_ECEF_norm^3*(Re/pos_ECEF_norm)^2.*...
    (5*(pos_ECEF(3)/pos_ECEF_norm)^2-[1;1;3]);
a_Jtwo=rotm\a_Jtwo;

%Composition of dxdt
dxdt(1:3)=[xdot; ydot; zdot];
dxdt(4:6)=-mu/r_norm^3*[x;y;z]+a_Jtwo;
end

function [Az,El,rho,rhodot]=antenna(stationName, et_vec, rv_sat_eci)
%Initialization
rv_sat_eci=rv_sat_eci';
topoFrame = [stationName, '_TOPO'];
disc=length(et_vec);
Az=zeros(disc,1);
El=zeros(disc,1);
rho=zeros(disc,1);
rhodot=zeros(disc,1);

for i=1:disc
%Rotation matrix: ECI to TOPO
ROT_ECI2TOPO = cspice_sxform('J2000', topoFrame, et_vec(i));

%Station position in ECI
rv_station_eci = cspice_spkezr(stationName, et_vec(i), 'J2000',...
    'NONE', 'EARTH');

%Satellite position with respect to station in ECI
rv_station_sat_eci = rv_sat_eci(:,i) - rv_station_eci;

%Satellite position with respect to station in TOPO
rv_station_sat_topo=ROT_ECI2TOPO*rv_station_sat_eci;

% Compute range, azimuth and elevation
rll_station_sat = cspice_xfmsta(rv_station_sat_topo,...
    'RECTANGULAR','LATITUDINAL','EARTH');

%Outputs
rho(i)=rll_station_sat(1);
Az(i)=rll_station_sat(2);
El(i)=rll_station_sat(3);
rhodot(i)=rll_station_sat(4);
end
end

function residual = Jcost1(x, stat, meas_real, tvec, i_visibility, sigma_meas)
%Propagate to measurement epochs (without J2)
mu=398600;
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~,xx_prop] = ode113(@(t,xx) RTBP(t,xx,mu),tvec,x,options);

%Compute predicted measurements (with only Perth)
[Az,El,rho] = antenna(stat, tvec, xx_prop);
meas_pred(:,1)=rho(i_visibility);
meas_pred(:,2)=Az(i_visibility)*cspice_dpr();
meas_pred(:,3)=El(i_visibility)*cspice_dpr();

%Compute the residual
W_m = diag(1./sigma_meas);
weighted_res = W_m*(meas_pred-meas_real)';
residual = weighted_res(:);
end

function residual = Jcost2(xx0, stat1, stat2, meas_real1, meas_real2, tvec, i_visibility1, i_visibility2, sigma_meas)  
%Propagate to measurement epochs (without J2)
mu=398600;
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~,xx_prop] = ode113(@(t,xx) RTBP(t,xx,mu),tvec,xx0,options);

%Compute predicted measurements (with both Perth and Kourou)
[Az1,El1,rho1] = antenna(stat1,tvec,xx_prop);
[Az2,El2,rho2] = antenna(stat2,tvec,xx_prop);
meas_pred1(:,1)=rho1(i_visibility1);
meas_pred1(:,2)=Az1(i_visibility1)*cspice_dpr();
meas_pred1(:,3)=El1(i_visibility1)*cspice_dpr();
meas_pred2(:,1)=rho2(i_visibility2);
meas_pred2(:,2)=Az2(i_visibility2)*cspice_dpr();
meas_pred2(:,3)=El2(i_visibility2)*cspice_dpr();

%Compute the residual
W_m = diag(1./sigma_meas);
residual2 = W_m*(meas_pred2-meas_real2)';
residual1 = W_m*(meas_pred1-meas_real1)';
residual = [residual2(:); residual1(:)];
end

function residual = Jcost3(xx0, stat1, stat2, meas_real1, meas_real2, tvec, i_visibility1, i_visibility2, sigma_meas)  
%Propagate to measurement epochs (with J2)
mu=398600;
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~,xx_prop] = ode113(@(t,xx)  RTBP_J2(t,xx,mu),tvec,xx0,options);

%Compute predicted measurements (with both Perth and Kourou)
[Az1,El1,rho1] = antenna(stat1,tvec,xx_prop);
[Az2,El2,rho2] = antenna(stat2,tvec,xx_prop);
meas_pred1(:,1)=rho1(i_visibility1);
meas_pred1(:,2)=Az1(i_visibility1)*cspice_dpr();
meas_pred1(:,3)=El1(i_visibility1)*cspice_dpr();
meas_pred2(:,1)=rho2(i_visibility2);
meas_pred2(:,2)=Az2(i_visibility2)*cspice_dpr();
meas_pred2(:,3)=El2(i_visibility2)*cspice_dpr();

%Compute the residual
W_m = diag(1./sigma_meas);
residual2 = W_m*(meas_pred2-meas_real2)';
residual1 = W_m*(meas_pred1-meas_real1)';
residual = [residual2(:); residual1(:)];
end