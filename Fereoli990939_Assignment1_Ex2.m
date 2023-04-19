%% Spacecraft Guidance and Navigation - Impulsive Guidance (2022/2023)
% Assignment:     1
% Exercise:       Impulsive Guidance
% Author:         Giovanni Fereoli
% ID number:      990939

%% Ex.1 - Guess Generation
clear; clc; close;

%Data
mu=1.21506683*1e-2;
alpha=1.5*pi;
beta=1.41;
delta=7;
ti=0;

%Guess generation in the rotating frame
[x0, Ti, Tf]=firstguess_gen(alpha, beta, delta, ti, mu);

%Change of r.f.: guess in the inertial frame Earth-centered
[X0]=change_frame(x0, Ti, mu);

%Integration
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[tt, xx] = ode113(@(t,x) PBRFBP_rhs(t, x, mu), [Ti Tf], x0, options);

%Plot in the rotating frame
figure(1)
plot(xx(:,1), xx(:,2),'b',...
    'LineWidth',2);
xlabel('x [DU]');
ylabel('y [DU]');
grid on;
real axis;
hold on;
plot(-mu, 0,'k.','MarkerSize',10);
plot(1-mu,0,'k.','MarkerSize',10);
text(-mu, 0,'EARTH','HorizontalAlignment',...
    'center','VerticalAlignment','top');
text(1-mu,0,'MOON','HorizontalAlignment',...
    'center','VerticalAlignment','top');

%Plot in the inertial r.f. earth-centered
disc=length(xx);
XX=zeros(disc,4);
for i=1:length(XX)
    XX(i,:)=change_frame(xx(i,:), tt(i), mu);
end
figure(2)
plot(XX(:,1), XX(:,2),'b', 'LineWidth',2);
xlabel('x [DU]');
ylabel('y [DU]');
grid on;
real axis;
hold on;
plot(0, 0,'k.','MarkerSize',10);
text(0, 0,'EARTH','HorizontalAlignment',...
    'center','VerticalAlignment','top');


%% Ex.2A - Single Shooting without Derivatives
clc; close all;

%Solving NLP
Xopt0=[x0(1:2); x0(3:4); Ti; Tf];
options=optimset('Algorithm','active-set','LargeScale','on',...
    'Display','iter','TolCon',1e-10, 'MaxIter',100,'MaxFunEvals',2000);
[Xopt, ~]=fmincon(@obj_ss, Xopt0, [],[],[],[],[],[], @nonlcons_ss,...
    options);

%Integration
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~, xx1] = ode113(@(t,x) PBRFBP_rhs(t, x, mu), [Xopt(5) Xopt(6)],...
    [Xopt(1); Xopt(2); Xopt(3); Xopt(4)], options);

%Plot in the rotating frame
figure(1);
plot(xx1(:,1), xx1(:,2),'b', 'LineWidth',2);
xlabel('x [DU]');
ylabel('y [DU]');
grid on;
real axis;
hold on;
plot(-mu, 0, 'k.','MarkerSize',10);
plot(1-mu,0,'k.','MarkerSize',10);
text(-mu, 0, 'EARTH','HorizontalAlignment',...
    'center','VerticalAlignment','top');
text(1-mu,0,'MOON','HorizontalAlignment',...
    'center','VerticalAlignment','top');

%% Ex.2B - Single Shooting with Derivatives
clc; close all;

%Solving NLP
Xopt0=[x0(1:2); x0(3:4); Ti; Tf];
options=optimset('Algorithm','active-set','LargeScale','on',...
    'Display','iter','TolCon',1e-10,'MaxIter',100,'MaxFunEvals',2000,...
    'GradObj','on','GradConstr','on');
[Xopt, ~]=fmincon(@obj_ss2, Xopt0, [],[],[],[],[],[],...
    @nonlcons_ss2, options);

%Integration
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~, xx2] = ode113(@(t,x) PBRFBP_rhs(t, x, mu), [Xopt(5) Xopt(6)],...
    [Xopt(1); Xopt(2); Xopt(3); Xopt(4)], options);

%Plot in the rotating frame of both Single Shooting solutions
figure(1);
plot(xx2(:,1), xx2(:,2),'r', 'LineWidth',2);
xlabel('x [DU]');
ylabel('y [DU]');
ylim([-2.5,3]);
grid on;
hold on;
real axis;
plot(-mu, 0,'k.','MarkerSize',10);
plot(1-mu,0,'k.','MarkerSize',10);
text(-mu, 0, 'EARTH','HorizontalAlignment',...
    'center','VerticalAlignment','top');
text(1-mu,0,'MOON','HorizontalAlignment',...
    'center','VerticalAlignment','top');

%% Ex.3 Multiple Shooting with Deriatives
clc; close all;

%Multiple shooting guess generation
Xopt0=[x0(1:2); x0(3:4); Ti; Tf];
[Xopt0_ms]=firstguess_gen_ms(Xopt0);

%Solving MLP
options=optimset('Algorithm','active-set','LargeScale','on',...
    'Display','iter','TolCon',1e-10,'MaxIter',1000,'MaxFunEvals',2000,...
    'GradObj','on','GradConstr','on');
[Xopt_ms, fval]=fmincon(@obj_ms, Xopt0_ms, [],[],[],[],[],[],...
    @nonlcons_ms, options);

%Integration
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~, xx3] = ode113(@(t,x) PBRFBP_rhs(t, x, mu), [Xopt_ms(17) Xopt_ms(18)],...
    [Xopt_ms(1); Xopt_ms(2); Xopt_ms(3); Xopt_ms(4)], options);

%Plot in the rotating frame
figure(1)
plot(xx3(:,1), xx3(:,2),'b', 'LineWidth',2);
hold on;
plot(Xopt_ms(1), Xopt_ms(2), 'r.', 'MarkerSize', 20);
plot(Xopt_ms(5), Xopt_ms(6), 'r.', 'MarkerSize', 20);
plot(Xopt_ms(9), Xopt_ms(10), 'r.', 'MarkerSize', 20);
plot(Xopt_ms(13), Xopt_ms(14), 'r.', 'MarkerSize', 20);
xlabel('x [DU]');
ylabel('y [DU]');
grid on;
real axis;
plot(-mu, 0, 'k.','MarkerSize',10);
plot(1-mu,0,'k.','MarkerSize',10);
text(-mu, 0, 'EARTH','HorizontalAlignment',...
    'center','VerticalAlignment','top');
text(1-mu,0,'MOON','HorizontalAlignment',...
    'center','VerticalAlignment','top');

%% Functions

%Ex.1
function [x0,Ti,ToF]=firstguess_gen(alpha, beta, delta, ti, mu)
%Data
Re=6378;            %Earth radius
hi=167;             %Earth Parking orbit
DU=3.84405000*1e5;  %Distance unit
r0=(Re+hi)/DU;

%Initialization
x0=zeros(4,1); 
v0=beta*sqrt((1-mu)/r0);             %Velocity magnitude (already with DU)

%Guess generation
x0(1)=r0*cos(alpha)-mu;              %State vector
x0(2)=r0*sin(alpha);
x0(3)=-(v0-r0)*sin(alpha);
x0(4)=(v0-r0)*cos(alpha);
Ti=ti;                               %Initial/Final time
ToF=delta;
end

function [X0]=change_frame(x0, Ti, mu)
%Initialization
X0=zeros(4,1);

%Unpack
x=x0(1);
y=x0(2);
xdot=x0(3);
ydot=x0(4);

%Change of r.f.
X0(1)=(x+mu)*cos(Ti)-y*sin(Ti);
X0(2)=(x+mu)*sin(Ti)+y*cos(Ti);
X0(3)=(xdot-y)*cos(Ti)-(ydot+x+mu)*sin(Ti);
X0(4)=(xdot-y)*sin(Ti)+(ydot+x+mu)*cos(Ti);
end

function [dxdt] = PBRFBP_rhs(t, x_vec, mu)
% Initialize
dxdt=zeros(4,1);  %I always initialize the output

x=x_vec(1);
y=x_vec(2);
xdot=x_vec(3);
ydot=x_vec(4);

r1_norm=sqrt((x+mu)^2+y^2);
r2_norm=sqrt((x+mu-1)^2+y^2);

%Additional data for PBRFBP
ms=3.28900541*1e5;
rho=3.88811143*1e2;
ws=-9.25195985*1e-1;

%Composition of a PCRTBP
dxdt(1:2)=[xdot;ydot];
dxdt(3:4)=[2*ydot+x-(1-mu)*(x+mu)/r1_norm^3-mu*(x+mu-1)/r2_norm^3;...
    -2*xdot+y-(1-mu)*y/r1_norm^3-mu*y/r2_norm^3];

%BRFBP through PCRTBP perturbation
r3=sqrt((x-rho*cos(ws*t))^2 + (y-rho*sin(ws*t))^2);
dxdt(3)=dxdt(3)-ms*(x-rho*cos(ws*t))/r3^3-ms*cos(ws*t)/rho^2;
dxdt(4)=dxdt(4)-ms*(y-rho*sin(ws*t))/r3^3-ms*sin(ws*t)/rho^2;

end

function [dxMdt]=STMexact_PBRFBP(t,xM,mu)
%Initialize
dxMdt=zeros(20,1);

x=xM(1);
y=xM(2);
xdot=xM(3);
ydot=xM(4);

M=(reshape(xM(5:20),4,4))';   %From equations to STM

%Additional data for PBRFBP
ms=3.28900541*1e5;
rho=3.88811143*1e2;
ws=-9.25195985*1e-1;

%PCRTBP dynamics
r1_norm=sqrt((x+mu)^2+y^2);
r2_norm=sqrt((x+mu-1)^2+y^2);

dxMdt(1:2)=[xdot;ydot];
dxMdt(3:4)=[2*ydot+x-(1-mu)*(x+mu)/r1_norm^3-mu*(x+mu-1)/r2_norm^3;...
    -2*xdot+y-(1-mu)*y/r1_norm^3-mu*y/r2_norm^3];

%PBRFBP through PCRTBP perturbation
r3=sqrt((x-rho*cos(ws*t))^2 + (y-rho*sin(ws*t))^2);
dxMdt(3)=dxMdt(3)-ms*(x-rho*cos(ws*t))/r3^3-ms*cos(ws*t)/rho^2;
dxMdt(4)=dxMdt(4)-ms*(y-rho*sin(ws*t))/r3^3-ms*sin(ws*t)/rho^2;

%Variational equations of PCRTBP
df1dx=1-(1-mu)/r1_norm^3+3*(1-mu)*(x+mu)^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*(x+mu-1)^2/r2_norm^5;
df1dy=3*(1-mu)*(x+mu)*y/r1_norm^5+3*mu*(x+mu-1)*y/r2_norm^5;
df2dy=1-(1-mu)/r1_norm^3+3*(1-mu)*y^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*y^2/r2_norm^5;

%Variational equations of BRFBP
df1dx=df1dx+ms*(2*(x-rho*cos(ws*t))^2-(y-rho*sin(ws*t))^2)/r3^5;
df1dy=df1dy+3*ms*(x-rho*cos(ws*t))*(y-rho*sin(ws*t))/r3^5;
df2dy=df2dy+ms*(2*(y-rho*sin(ws*t))^2-(x-rho*cos(ws*t))^2)/r3^5;

A=[0, 0, 1, 0;...
    0, 0, 0, 1;...
    df1dx, df1dy, 0, 2;...
    df1dy, df2dy, -2, 0];

dMdt=A*M;
dxMdt(5:8)=dMdt(1,1:4)';
dxMdt(9:12)=dMdt(2,1:4)';
dxMdt(13:16)=dMdt(3,1:4)';
dxMdt(17:20)=dMdt(4,1:4)';
end


%Ex.2A - Without Derivatives (ODEs without STM to improve speed here)

function [f] = obj_ss(Xopt)
%Data
mu=1.21506683*1e-2;

%Upacking state
xi=Xopt(1);
yi=Xopt(2);
xidot=Xopt(3);
yidot=Xopt(4);
Ti=Xopt(5);
Tf=Xopt(6);

%Earth orbit
Re=6378;            %Earth radius
hi=167;             %Earth Parking orbit
DU=3.84405000*1e5;  %Distance unit
ri=(Re+hi)/DU;

%Integration
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~, Xf] = ode113(@(t,x) PBRFBP_rhs(t, x, mu), [Ti Tf], [xi; yi; xidot;...
    yidot], options);
Xf=Xf(end,:);

%Final state
xf=Xf(1);
yf=Xf(2);
xfdot=Xf(3);
yfdot=Xf(4);

%Moon orbit
Rm=1738;            %Earth radius
hf=100;             %Earth Parking orbit
rf=(Rm+hf)/DU;

%OF
dV1=abs(sqrt((xidot-yi)^2+(yidot+xi+mu)^2)-sqrt((1-mu)/ri));
dV2=abs(sqrt((xfdot-yf)^2+(yfdot+xf+mu-1)^2)-sqrt(mu/rf));
f=dV1+dV2;
end

function [c, ceq] = nonlcons_ss(Xopt)
%Data
mu=1.21506683*1e-2;

%Upacking state
xi=Xopt(1);
yi=Xopt(2);
xidot=Xopt(3);
yidot=Xopt(4);
Ti=Xopt(5);
Tf=Xopt(6);

%Earth orbit
Re=6378;            %Earth radius
hi=167;             %Earth Parking orbit
DU=3.84405000*1e5;  %Distance unit
ri=(Re+hi)/DU;

%Integration
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~, Xf] = ode113(@(t,x) PBRFBP_rhs(t, x, mu), [Ti Tf], [xi; yi; xidot;...
    yidot], options);
Xf=Xf(end, :);

%Final state
xf=Xf(1);
yf=Xf(2);
xfdot=Xf(3);
yfdot=Xf(4);

%Moon orbit
Rm=1738;            %Earth radius
hf=100;             %Earth Parking orbit
rf=(Rm+hf)/DU;

%C
c=[];

%Ceq
ceq=[(xi+mu)^2+yi^2-ri^2;...
    (xi+mu)*(xidot-yi)+yi*(yidot+xi+mu);...
    (xf+mu-1)^2+yf^2-rf^2;...
    (xf+mu-1)*(xfdot-yf)+yf*(yfdot+xf+mu-1)];
end

%Ex.2B - With Derivatives

function [f, df] = obj_ss2(Xopt)
%Data
mu=1.21506683*1e-2;

%Upacking state
xi=Xopt(1);
yi=Xopt(2);
xidot=Xopt(3);
yidot=Xopt(4);
Ti=Xopt(5);
Tf=Xopt(6);

%Earth orbit
Re=6378;            %Earth radius
hi=167;             %Earth Parking orbit
DU=3.84405000*1e5;  %Distance unit
ri=(Re+hi)/DU;

%Integration
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~, xxM] = ode113(@(t,xM) STMexact_PBRFBP(t,xM, mu), [Ti Tf],...
        [[xi; yi; xidot; yidot]; reshape(eye(4),[],1)], options);
Xf=xxM(end,1:4);
STMf=(reshape(xxM(end,5:20),4,4))';   %From equations to STM

%Final state
xf=Xf(1);
yf=Xf(2);
xfdot=Xf(3);
yfdot=Xf(4);

%Moon orbit
Rm=1738;            %Earth radius
hf=100;             %Earth Parking orbit
rf=(Rm+hf)/DU;

%OF
dV1=abs(sqrt((xidot-yi)^2+(yidot+xi+mu)^2)-sqrt((1-mu)/ri));
dV2=abs(sqrt((xfdot-yf)^2+(yfdot+xf+mu-1)^2)-sqrt(mu/rf));
f=dV1+dV2;

%DOF
if nargout>1
    
    %Velocity at Ti
    vi=sqrt((xidot-yi)^2+(yidot+xi+mu)^2);
    vxi=xidot-yi; 
    vyi=yidot+xi+mu;

    %Velocity at Tf
    vf=sqrt((xfdot-yf)^2+(yfdot+xf+mu-1)^2);
    vxf=xfdot-yf; 
    vyf=yfdot+xf+mu-1;
    
    %Accelerations at Ti, Tf in 2D
    rhs1=PBRFBP_rhs(Ti,[xi; yi; xidot; yidot], mu);
    rhs2=PBRFBP_rhs(Tf,Xf, mu);

    %Components df
    dfdX=[vyi/vi; -vxi/vi; vxi/vi; vyi/vi]+STMf'*...
        [vyf/vf; -vxf/vf; vxf/vf; vyf/vf];
    dfdti=-(STMf*rhs1)'*[vyf/vf; -vxf/vf; vxf/vf; vyf/vf];
    dfdtf=rhs2'*[vyf/vf; -vxf/vf; vxf/vf; vyf/vf];

    df=[dfdX; dfdti; dfdtf];
end
end

function [c, ceq, gradc, gradceq] = nonlcons_ss2(Xopt)
%Data
mu=1.21506683*1e-2;

%Upacking state
xi=Xopt(1);
yi=Xopt(2);
xidot=Xopt(3);
yidot=Xopt(4);
Ti=Xopt(5);
Tf=Xopt(6);

%Earth orbit
Re=6378;            %Earth radius
hi=167;             %Earth Parking orbit
DU=3.84405000*1e5;  %Distance unit
ri=(Re+hi)/DU;

%Integration
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~, xxM] = ode113(@(t,xM) STMexact_PBRFBP(t,xM, mu), [Ti Tf],...
        [[xi; yi; xidot; yidot]; reshape(eye(4),[],1)], options);
Xf=xxM(end,1:6);
STMf=(reshape(xxM(end,5:20),4,4))';   %From equations to STM

%Final state
xf=Xf(1);
yf=Xf(2);
xfdot=Xf(3);
yfdot=Xf(4);

%Moon orbit
Rm=1738;            %Earth radius
hf=100;             %Earth Parking orbit
rf=(Rm+hf)/DU;

%C
c=[];

%Ceq
ceq=[(xi+mu)^2+yi^2-ri^2;...
    (xi+mu)*(xidot-yi)+yi*(yidot+xi+mu);...
    (xf+mu-1)^2+yf^2-rf^2;...
    (xf+mu-1)*(xfdot-yf)+yf*(yfdot+xf+mu-1)];

if nargout>2
    %Gradient inequality constraints
    gradc=[];

    %Accelerations at Ti, Tf in 2D
    rhs1=PBRFBP_rhs(Ti,[xi; yi; xidot; yidot], mu);
    rhs2=PBRFBP_rhs(Tf,Xf, mu);

    %Gradient equality constraints
    gradceq_12=[2*(xi+mu), 2*yi, 0, 0, 0, 0;...
        xidot, yidot, xi+mu, yi, 0, 0];

    gradceq_3_state=2*[xf+mu-1, yf, 0, 0]*STMf;
    gradceq_3_time=2*[xf+mu-1, yf, 0, 0]*[-STMf*rhs1, rhs2];
    gradceq_3=[gradceq_3_state, gradceq_3_time];

    gradceq_4_state=[xfdot, yfdot, xf+mu-1, yf]*STMf;
    gradceq_4_time=[xfdot, yfdot, xf+mu-1, yf]*[-STMf*rhs1, rhs2];
    gradceq_4=[gradceq_4_state, gradceq_4_time];

    gradceq=[gradceq_12; gradceq_3; gradceq_4];
    gradceq=gradceq';
end

end

%Ex.3

function [Xopt0_ms]=firstguess_gen_ms(Xopt0)
%Data
mu=1.21506683*1e-2;

%Unpacking state
X10=[Xopt0(1); Xopt0(2); Xopt0(3); Xopt0(4)];
T1=Xopt0(5);
T4=Xopt0(6);

%Definition of T2, T3
T2=T1+(T4-T1)/3;
T3=T1+2*(T4-T1)/3;

%Integrations
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[~, xxM1] = ode113(@(t,xM) STMexact_PBRFBP(t,xM, mu), [T1, T2, T3, T4],...
        [X10; reshape(eye(4),[],1)], options);
X10=xxM1(1,1:4);
X20=xxM1(2,1:4);
X30=xxM1(3,1:4);
X40=xxM1(4,1:4);

%Composition
Xopt0_ms=zeros(18,1);
Xopt0_ms(1:4,1)=X10;
Xopt0_ms(5:8,1)=X20;
Xopt0_ms(9:12,1)=X30;
Xopt0_ms(13:16,1)=X40;
Xopt0_ms(17,1)=T1;
Xopt0_ms(18,1)=T4;
end

function [f, df] = obj_ms(Xopt)
%Data
mu=1.21506683*1e-2;

%Upacking state
x1=Xopt(1);
y1=Xopt(2);
x1dot=Xopt(3);
y1dot=Xopt(4);

%Earth orbit
Re=6378;            %Earth radius
hi=167;             %Earth Parking orbit
DU=3.84405000*1e5;  %Distance unit
ri=(Re+hi)/DU;

%Final state
x4=Xopt(13);
y4=Xopt(14);
x4dot=Xopt(15);
y4dot=Xopt(16);

%Moon orbit
Rm=1738;            %Earth radius
hf=100;             %Earth Parking orbit
rf=(Rm+hf)/DU;

%OF
dV1=abs(sqrt((x1dot-y1)^2+(y1dot+x1+mu)^2)-sqrt((1-mu)/ri));
dV2=abs(sqrt((x4dot-y4)^2+(y4dot+x4+mu-1)^2)-sqrt(mu/rf));
f=dV1+dV2;

%DOF
if nargout>1
    P1=[y1dot+x1+mu; y1-x1dot; x1dot-y1; y1dot+x1+mu]/...
        sqrt((x1dot-y1)^2+(y1dot+x1+mu)^2);
    P4=[y4dot+x4+mu-1; y4-x4dot; x4dot-y4; y4dot+x4+mu-1]/...
        sqrt((x4dot-y4)^2+(y4dot+x4+mu-1)^2);

    df=[P1; zeros(8,1); P4; 0; 0];
end
end

function [c, ceq, gradc, gradceq] = nonlcons_ms(Xopt)
%Data
mu=1.21506683*1e-2;

%Upacking state
X1=Xopt(1:4);
X2=Xopt(5:8);
X3=Xopt(9:12);
X4=Xopt(13:16);
Ti=Xopt(17);
Tf=Xopt(18);

%Time instant
T1=Ti;
T2=T1+(Tf-Ti)/3;
T3=T1+2*(Tf-Ti)/3;
T4=Tf;

%Earth orbit
Re=6378;            %Earth radius
hi=167;             %Earth Parking orbit
DU=3.84405000*1e5;  %Distance unit
ri=(Re+hi)/DU;

%Moon orbit
Rm=1738;            %Earth radius
hf=100;             %Earth Parking orbit
rf=(Rm+hf)/DU;

%Integrations
options = odeset('RelTol', 2.24*1e-12, 'AbsTol', 2.24*1e-12);
%Flow 1
[~, xxM1] = ode113(@(t,xM) STMexact_PBRFBP(t,xM, mu), [T1 T2],...
        [[X1(1); X1(2); X1(3); X1(4)]; reshape(eye(4),[],1)],...
        options);
Xf1=(xxM1(end,1:4))';
STMf1=(reshape(xxM1(end,5:20),4,4))';   %From equations to STM

%Flow 2
[~, xxM2] = ode113(@(t,xM) STMexact_PBRFBP(t,xM, mu), [T2 T3],...
        [[X2(1); X2(2); X2(3); X2(4)]; reshape(eye(4),[],1)],...
        options);
Xf2=(xxM2(end,1:4))';
STMf2=(reshape(xxM2(end,5:20),4,4))';   %From equations to STM

%Flow 3
[~, xxM3] = ode113(@(t,xM) STMexact_PBRFBP(t,xM, mu), [T3 T4],...
        [[X3(1); X3(2); X3(3); X3(4)]; reshape(eye(4),[],1)],...
        options);
Xf3=(xxM3(end,1:4))';
STMf3=(reshape(xxM3(end,5:20),4,4))';   %From equations to STM

%C
c=[];

%Ceq
ceq=[Xf1-X2;...
    Xf2-X3;...
    Xf3-X4;...
    (X1(1)+mu)^2+X1(2)^2-ri^2;...
    (X1(1)+mu)*(X1(3)-X1(2))+X1(2)*(X1(4)+X1(1)+mu);...
    (X4(1)+mu-1)^2+X4(2)^2-rf^2;...
    (X4(1)+mu-1)*(X4(3)-X4(2))+X4(2)*(X4(4)+X4(1)+mu-1)];

if nargout>2
    %Gradient inequality costraint
    gradc=[];

    %Gradient equality constraint
    %Phi_i and Phi_f wrt X
    R1=[2*(X1(1)+mu), 2*X1(2), 0, 0;...
        X1(3), X1(4), X1(1)+mu, X1(2)];
    R4=[2*(X4(1)+mu-1), 2*X4(2), 0, 0;...
        X4(3), X4(4), X4(1)+mu-1, X4(2)];

    %Eps_1,2,3 wrt Ti and Tf
    rhs1=PBRFBP_rhs(T1, [X1(1); X1(2); X1(3); X1(4)], mu);
    rhs2=PBRFBP_rhs(T2, [X2(1); X2(2); X2(3); X2(4)], mu);
    rhs3=PBRFBP_rhs(T3, [X3(1); X3(2); X3(3); X3(4)], mu);

    rhs_phi1=PBRFBP_rhs(T2, [Xf1(1); Xf1(2); Xf1(3); Xf1(4)], mu); 
    rhs_phi2=PBRFBP_rhs(T3, [Xf2(1); Xf2(2); Xf2(3); Xf2(4)], mu); 
    rhs_phi3=PBRFBP_rhs(T4, [Xf3(1); Xf3(2); Xf3(3); Xf3(4)], mu);

    Q11=-STMf1*rhs1+(2/3)*rhs_phi1;
    Q12=-(2/3)*STMf2*rhs2+(1/3)*rhs_phi2;
    Q13=-(1/3)*STMf3*rhs3;
    Q41=(1/3)*rhs_phi1;
    Q42=-(1/3)*STMf2*rhs2+(2/3)*rhs_phi2;
    Q43=-(2/3)*STMf3*rhs3+rhs_phi3;

    %Composition of GradCeq
    gradceq1=[STMf1, -eye(4), zeros(4), zeros(4), Q11, Q41];
    gradceq2=[zeros(4), STMf2, -eye(4), zeros(4), Q12, Q42];
    gradceq3=[zeros(4), zeros(4), STMf3, -eye(4), Q13, Q43];
    gradceq4=[R1, zeros(2,4), zeros(2,4), zeros(2,4), zeros(2)];
    gradceq5=[zeros(2,4), zeros(2,4), zeros(2,4), R4, zeros(2)];
    gradceq=[gradceq1; gradceq2; gradceq3; gradceq4; gradceq5];
    gradceq=gradceq';
end
end