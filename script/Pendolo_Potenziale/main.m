%%
clc
clear all
close all

%OPTIMIZATION MAIN

%Remark: in case of Keplerian potential if eps=0 we have the
%    case of system without repulsive potentials
%       in case of Mollifier potential it is enough to put
%    C=0 to avoid repulsion
%       in case of Vn it is enough to put n=0 to avoid
%    repulsion

global M m g l yp0 xp0 T xF thetaF x0 theta0 eps N c n k


%Costant parameters:
g=9.81;     %gravity accelleration, [m/s^2]
l=1;        %Length of the Pendulum [m]
M=1;        %Mass of the cart [kg]
m=0.01;     %Mass of the pendulum [kg]
T=1.5;        %Final time [s]
yp0=1;    %y of the obstacle
xp0=0;      %x of the obstacle
eps=0.1;      %Eps parameter for keplerian case, eps=0 in case of not potential
%eps=1;    %Eps parameter for Mollifier case
%c=1;      %Constant of mollifier, put C=0 in case of not potential
%k=0.1;    %Parameter of the Mollifier
%n=1000;      %Parameter in case of Vn, put n=0 in case of not repulsion
%k=0.1;    %Parameter in case of Vn
xF=2;       %Final configuration of x
thetaF=0;   %Final configuartion of theta
x0=-2;       %Intial configuration of x
theta0=0;   %Final configuration of theta

%Time discretization:
N=1000;
tspan=linspace(0,T,N);



%Optimization:
pOpt0=[5;0];    %intial guess

tic;
options = ...
    optimoptions('fsolve','Algorithm','trust-region-dogleg',...
        'Display','iter','MaxFunEvals',1000);
pOpt=fsolve(@funerr,pOpt0,options);




%Final computation of the dynamical of the system using pOpt:
K=[pOpt(1);pOpt(2)];
q0=[x0;0;0;0;theta0;0;K(1);K(2)];
[t,q]=ode45(@(t,q) System(t,q),tspan,q0);
toc



xOpt=q(:,1);        %optimal x
thetaOpt=q(:,5);    %optimal theta
y=l*cos(thetaOpt);  %y of pendulum

%Compute the control:
u=q(:,3).*(M+m*sin(q(:,5)).^2)-...
    m*sin(q(:,5)).*(l*q(:,6).^2-g*cos(q(:,5)));


%Plots:

figure(1)
plot(t,xOpt,'b',t,thetaOpt,'k',linspace(0,T,length(t)),...
     yp0*ones(length(t),1),'--k')
legend('x','theta','obstacle (y)')
title('x and theta solutions of equations')
xlabel('t')


figure(2)
plot(t,y,'b',linspace(0,T,length(t)),yp0*ones(length(t),1),'--k')
legend('y','obstacle (y)')
title('y of the Pendulum')
xlabel('t')


figure(3)
plot(t,u,'b',linspace(0,T,length(t)),0*ones(length(t),1),'--k')
axis([0,T,-20,100])
legend('control zoom')
title('Zoom of the Control')
xlabel('t')


figure(4)
plot(t,xOpt,'b',t,y,'k',linspace(0,T,length(t)),...
     yp0*ones(length(t),1),'--k',linspace(0,T,length(t)),...
     xp0*ones(length(t),1),'--r')
legend('x','y','obstacle (y)','obstacle (x)')
title('x of the Cart and y of the Pendulum')
xlabel('t')


figure(5)
plot(t,u,'k')
title('Control u')
xlabel('t')
ylabel('u')