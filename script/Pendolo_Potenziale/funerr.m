function F=funerr(K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function funerr to be minimize via optimization technique, it represents
%the difference between the computed final configurations and the precribed
%ones
%INPUT:
%K...parameter of the optimization
%OUTPUT:
%F...the function that we want to minimize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global M m g l yp0 xp0 T xF thetaF x0 theta0 eps N

tt=linspace(0,T,N);%time on which the system is integrated
q0=[x0;0;0;0;theta0;0;K(1);K(2)];%initial data

[t,qk]=ode45(@(t,q) System(t,q),tt,q0); %solve the system

%Compute the function F:
Ffx=xF-qk(end,1);
Fft=thetaF-qk(end,5);
F=[Ffx;Fft];
end