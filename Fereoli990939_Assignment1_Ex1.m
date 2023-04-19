%% Spacecraft Guidance and Navigation - Periodic Orbit (2022/2023)
% Assignment:     1
% Exercise:       Periodic Orbit
% Author:         Giovanni Fereoli
% ID number:      990939

%% Ex.1 - L2 Position
clear; clc; close;

%Initialization
mu=3.0359*1e-6; %Sun-Earth CR3BP

%L2 position
xL2=L2position(mu);
   
%Plots
figure(1)
plot3(xL2, 0, 0,'k.','MarkerSize',10);
xlabel('x [DU]');
ylabel('y [DU]');
zlabel('z [DU]');
grid on;
real axis;
view(2);
hold on;
plot3(1-mu,0,0,'k.','MarkerSize',10);
text(L2position(mu), 0, 0,'L2','HorizontalAlignment',...
    'center','VerticalAlignment','top');
text(1-mu,0,0,'EARTH','HorizontalAlignment',...
    'center','VerticalAlignment','top')

%% Ex.2 - Single Halo Orbit
clear; clc; close;

%Initialization
mu=3.0359*1e-6; %Sun-Earth 3BP
options1 = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
options2 = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14,...
        'Events',@ev_orb2D);

%Initial condition
x0 = [1.008296144180133; 0;0.001214294450297; 0; 0.010020975499502; 0];

%Check
check_value = STMexact_CRTBP(0, [x0; reshape(eye(6),[],1)], mu);
    
%Differential correction scheme
xxM=ones(1,6);
dX0=zeros(6,1);

while abs(xxM(end,4))>1e-12 || abs(xxM(end,6))>1e-12
    %Integration
    [ttM, xxM] = ode113(@(t,xM) STMexact_CRTBP(t,xM, mu), [0 2*pi],...
    [x0+dX0; reshape(eye(6),[],1)], options2);
    
    %Reshape
    STMf=(reshape(xxM(end,7:42),6,6))';
    
    %Correction
    Corr=[STMf(4,1), STMf(4,5); STMf(6,1), STMf(6,5)]\...
    [-xxM(end,4);-xxM(end,6)];
    dX0=dX0+[Corr(1); 0; 0; 0; Corr(2); 0];
end

dX0=zeros(6,1);
%Correction initial guess
X0_corr=x0+dX0;

%Final integration
[ttMcorr_full, xxMcorr_full] = ode113(@(t,xM) STMexact_CRTBP(t,xM,mu),...
[0 pi],[x0+dX0; reshape(eye(6),[],1)], options1);
xxcorr_full=xxMcorr_full(:,1:6);    

%Plots
figure(1)
plot3(xxcorr_full(:,1), xxcorr_full(:,2),xxcorr_full(:,3),'b',...
'LineWidth',1);
xlabel('x [DU]');
ylabel('y [DU]');
zlabel('z [DU]');
grid on;
real axis;
hold on;
plot3(L2position(mu), 0, 0,'k.','MarkerSize',10);
plot3(1-mu,0,0,'k.','MarkerSize',10);
text(L2position(mu), 0, 0,'L2','HorizontalAlignment',...
'center','VerticalAlignment','top');
text(1-mu,0,0,'EARTH','HorizontalAlignment',...
'center','VerticalAlignment','top');
view(2);
figure(2);
plot3(xxcorr_full(:,1), xxcorr_full(:,2),xxcorr_full(:,3),'b',...
'LineWidth',1);
xlabel('x [DU]');
ylabel('y [DU]');
zlabel('z [DU]');
grid on;
real axis;
hold on;
plot3(L2position(mu), 0, 0,'k.','MarkerSize',10);
plot3(1-mu,0,0,'k.','MarkerSize',10);
text(L2position(mu), 0, 0,'L2','HorizontalAlignment',...
'center','VerticalAlignment','top');
text(1-mu,0,0,'EARTH','HorizontalAlignment',...
'center','VerticalAlignment','top');
view(+37.5, 30);

%% Ex.3 - Set of Halo Orbits
clear; clc; close;

%Initialization
mu=3.0359*1e-6; %Sun-Earth 3BP
dX0=zeros(6,1);
options1 = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
options2 = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14,...
        'Events',@ev_orb2D);
x0new=zeros(6,1);
err=zeros(17,1);
i=0;

%Initial condition
x0 = [1.008296144180133; 0;0.001214294450297; 0; 0.010020975499502; 0];

%Halo orbits family search (Numerical Continuation)
for z0=0.001214294450297:0.0002:0.0046  
    %Initialization
    i=i+1;
    x0(3)=z0;

    %Differential correction
    xxM=ones(1,6);
    dX0=zeros(6,1);

    while abs(xxM(end,4))>1e-12 || abs(xxM(end,6))>1e-12
        %Integration
        [ttM, xxM] = ode113(@(t,xM) STMexact_CRTBP(t,xM, mu), [0 2*pi],...
        [x0+dX0; reshape(eye(6),[],1)], options2);

        %Reshape
        STMf=(reshape(xxM(end,7:42),6,6))';

        %Correction
        Corr=[STMf(4,1), STMf(4,5); STMf(6,1), STMf(6,5)]\...
            [-xxM(end,4);-xxM(end,6)];
        dX0=dX0+[Corr(1); 0; 0; 0; Corr(2); 0];
    end
   
    %Final integration
    [ttMcorr_full, xxMcorr_full] = ode113(@(t,xM) STMexact_CRTBP(t,xM,mu),...
        [0 pi],[x0+dX0; reshape(eye(6),[],1)], options1);
    xxcorr_full=xxMcorr_full(:,1:6);  
    
    %Plots
    figure(1)
    plot3(xxcorr_full(:,1), xxcorr_full(:,2),xxcorr_full(:,3),'b',...
        'LineWidth',1);
    xlabel('x [DU]');
    ylabel('y [DU]');
    zlabel('z [DU]');
    grid on;
    real axis;
    hold on;
    plot3(L2position(mu), 0, 0,'k.','MarkerSize',10);
    plot3(1-mu,0,0,'k.','MarkerSize',10);
    text(L2position(mu), 0, 0,'L2','HorizontalAlignment',...
        'center','VerticalAlignment','top');
    text(1-mu,0,0,'EARTH','HorizontalAlignment',...
        'center','VerticalAlignment','top');
    view(0,0);
    figure(2);
    plot3(xxcorr_full(:,1), xxcorr_full(:,2),xxcorr_full(:,3),'b',...
        'LineWidth',1);
    xlabel('x [DU]');
    ylabel('y [DU]');
    zlabel('z [DU]');
    grid on;
    real axis;
    hold on;
    plot3(L2position(mu), 0, 0,'k.','MarkerSize',10);
    plot3(1-mu,0,0,'k.','MarkerSize',10);
    text(L2position(mu), 0, 0,'L2','HorizontalAlignment',...
        'center','VerticalAlignment','top');
    text(1-mu,0,0,'EARTH','HorizontalAlignment',...
        'center','VerticalAlignment','top');
   view(+37.5, 30);

   %Numerical continuation
    x0=x0+dX0;
end


%% Functions

function xL2=L2position(mu)

%Gradient of U in x
f=@(x) x-(1-mu)*(x+mu)/(x+mu)^3-mu*(x+mu-1)/(x+mu-1)^3;

%Inital guesses for L2
xL20=1.5;

%Zeros computation
options = optimset('TolX',1e-22);
xL2=fzero(f, xL20, options);
end

function [value,isterminal,direction]=ev_orb2D(~,x_vec,~)

value=x_vec(2);   %When happens that y=0...
isterminal=1;
direction=-1;

end

function dxMdt=STMexact_CRTBP(~,xM,mu)

% -------------------------------------------------------------------------
% STMexact_CRTBP.m
%
% PROTOTYPE:
%   [dxMdt]=STMexact_CRTBP(~,xM,mu)
%
% DESCRIPTION:
%  State Transition Matrix computed through a Variational Approach needs
%  an integration of a [42x1] right-hand side. This function gives as 
%  output its vectorial field.
%
% INPUT:
%   t      : [1,1] Ephemeris time
%   xM     : [42,1] Cartesian state and STM vectorially re-arranged
%   mu     : [1]   Gravitational parameter
%
% OUTPUT:
%   dxMdt   : [6,1] RHS for STM variational approach
%
% CALLED FUNCTIONS:
%   None
%
% REFERENCES:
%   - Battin R., "An Introduction to the Mathematics and Methods of
%       Astrodynamics, Revised Edition", 1999.
%   - Curtis H., "Orbital Mechanics for Engineering Students, 4th Edition",
%       2019.
%
% FUTURE DEVELOPMENT:
%   None
%
% ORIGINAL VERSION:
%   Giovanni Fereoli, 30/09/2022, MATLAB, STMexact_CRTBP.m
%
% AUTHOR:
%   Giovanni Fereoli, 30/09/2022, MATLAB, STMexact_CRTBP.m
%
% CHANGELOG:
%  None
%
% Note: Please if you have got any changes that you would like to be done,
%   do not change the function, please contact the author.
%
% -------------------------------------------------------------------------

%Initialize
dxMdt=zeros(42,1);

x=xM(1);
y=xM(2);
z=xM(3);
xdot=xM(4);
ydot=xM(5);
zdot=xM(6);

M=(reshape(xM(7:42),6,6))';   %From equations to STM

% CRTBP dynamics
r1_norm=sqrt((x+mu)^2+y^2+z^2);
r2_norm=sqrt((x+mu-1)^2+y^2+z^2);

dxMdt(1:3)=[xdot;ydot;zdot];
dxMdt(4:6)=[2*ydot+x-(1-mu)*(x+mu)/r1_norm^3-mu*(x+mu-1)/r2_norm^3;...
    -2*xdot+y-(1-mu)*y/r1_norm^3-mu*y/r2_norm^3;...
    -(1-mu)*z/r1_norm^3-mu*z/r2_norm^3];

% Variational equations
df4dx=1-(1-mu)/r1_norm^3+3*(1-mu)*(x+mu)^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*(x+mu-1)^2/r2_norm^5;
df4dy=3*(1-mu)*(x+mu)*y/r1_norm^5+3*mu*(x+mu-1)*y/r2_norm^5;
df4dz=3*(1-mu)*(x+mu)*z/r1_norm^5+3*mu*(x+mu-1)*z/r2_norm^5;
df5dy=1-(1-mu)/r1_norm^3+3*(1-mu)*y^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*y^2/r2_norm^5;
df5dz=3*(1-mu)*y*z/r1_norm^5+3*mu*y*z/r2_norm^5;
df6dz=-(1-mu)/r1_norm^3+3*(1-mu)*z^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*z^2/r2_norm^5;

A=[0, 0, 0, 1, 0, 0;...
    0, 0, 0, 0, 1, 0;...
    0, 0, 0, 0, 0, 1;...
    df4dx, df4dy, df4dz, 0, 2, 0;...
    df4dy, df5dy, df5dz, -2, 0, 0;...
    df4dz, df5dz, df6dz, 0, 0, 0];

dMdt=A*M;
dxMdt(7:12)=dMdt(1,1:6)';
dxMdt(13:18)=dMdt(2,1:6)';
dxMdt(19:24)=dMdt(3,1:6)';
dxMdt(25:30)=dMdt(4,1:6)';
dxMdt(31:36)=dMdt(5,1:6)';
dxMdt(37:42)=dMdt(6,1:6)';

end




