clear 
close all
clc

%% First Example

% Function g and h
g = @(x) eye(2) ; 
h =@(x ,xO  , pO) norm(x - xO,2) - pO ;
gradh =@(x ,xO) (x-xO)/norm(x-xO,2) ; 

% Initial and final points
p0 = [0,0]';
pg = [5,5]';

% Obstacle
pO = [3, 2.8]'; % center
rO = 0.7;  % radius
r1 = 0.7 ; 
theta_span = linspace(0,2*pi,100);
% Variable used to plot the circle
xO = pO(1) + rO * cos(theta_span);
yO = pO(2) + rO * sin(theta_span);

% Parameters
Ts = 0.01;
rate = rateControl(1/Ts);  % control rate
goal_tol = 0.05;
K = 5;  % control gain
alpha_values = [ 0.5 , 1 , 2 , 20 ,100 ] ; 

subplot(121)

counter = 1 ; 

%Creates a cela array that cointans in every the trajectory for a specific value of alpha
xout_all = cell(length(alpha_values) ,1 ) ; 

for i = 1:length(alpha_values)
    %Inizialising the variables for every alpha value 

    x= p0;
    goal_dist = norm(x-pg,2);
    alpha = alpha_values(i) ; 
    xout = [] ; 
    %Control Loop 
    while goal_dist >= goal_tol
      
        u_des = K*(pg-x)/norm(pg-x,2);  % desired input
        
        phi =  gradh(x,pO)' * u_des + alpha * h(x,pO,r1) ;
      
        if (phi < 0)
             u_safe = - (g(x) * gradh(x,pO))/(norm(g(x)*gradh(x,xO))) * phi  ;
             u = u_safe + u_des ; 
        else
            u = u_des ; 
        end
        
        %Trajectory 
        x = x + Ts * g(x) * [u(1), u(2)]';  % system (Explicit Euler) 
        xout(: , counter) = x ; 
        
        counter=counter +1 ; 
    
        goal_dist = norm(x-pg,2);  % update goal dist
        
        %{
        clf 
        %Scenarion
        plot(p0(1),p0(2),'bo','LineWidth',1.2), hold on,
        plot(pg(1),pg(2),'rx','LineWidth',1.2), hold on, 
        %plot(xO,yO,'k-','LineWidth',1.2 ), hold on,
        fill(xO,yO , 'black')
        scatter(x(1),x(2),'blue','filled','MarkerEdgeColor','k')
        axis equal
        
        %Wait for sample time
        waitfor(rate);
        %}

    end
   
    xout_all{i}=xout ; 
end
 hold on ; 
plot(p0(1),p0(2),'bx','LineWidth',1.2)
plot(pg(1),pg(2),'rx','LineWidth',1.2) 
fill(xO,yO, 'black' , 'Facealpha' , 0.5)
%scatter(x(1),x(2),'blue','filled','MarkerEdgeColor','k')
axis equal

x = linspace (p0(1) , pg(1) ) ; 
y = linspace (p0(2) , pg(2) ) ; 
plot (x,y  ,'LineStyle','--' , 'LineWidth', 1.1)

colors_1 = ['r', 'g', 'b' , 'c' , 'k']; % Colori per le traiettorie
for i = 1:length(xout_all)
    trajectory = xout_all{i};
    plot(trajectory(1, :), trajectory(2, :), 'Color', colors_1(i), 'LineWidth', 1.1);
end
for i=1:length(alpha_values)
legend_alpha{i} = sprintf ("alpha = %0.1f " , alpha_values(i)); 
end

legend ("Starting Point","Goal Point" , "Obstacle "  , "Desired trjectory", legend_alpha{:} ,  'Location','northwest')
axis equal ; 
xlabel("X[m]");
ylabel("Y[m") ; 




%% CBF-APF


% Initial and final points
p0 = [0,0]';
pg = [5,5]';

% Obstacle
pO = [3, 2.8]'; % center
rO = 0.7;  % radius
theta_span = linspace(0,2*pi,100);
xO = pO(1) + rO * cos(theta_span);
yO = pO(2) + rO * sin(theta_span);

% Parameters
Ts = 0.02;
%rate = rateControl(1/Ts);  % control rate
goal_tol = 0.05;
K = 5;  % control gain
alpha_1 =  3 ; 
Katt = 2 ; 
Krep = 2 ; 
rho_values = [0.7 , 1.2 ] ; 
sigma0 = 10^(-1) ;

 % Function definition 
 Uatt = @(x) 0.5*Katt*norm(x-pg,2)^2 ; 
 gradUatt = @(x) Katt * (x-pg) ; 
 
 rho = @(x,pO,rO) norm(x - pO,2) - rO ; 
 
 Urep = @(x,pO,rO) 0.5*Krep * (1/rho(x,pO,rO) - 1/rO )^2 ;
 gradUrep = @(x,pO,rO)  Krep*(1/rho(x,pO,rO) - 1/rO )*(-1/(rho(x,pO,rO)^2) * (x-pO/norm(x - pO,2))) ; 
 
 g= @(x) eye(2) ; 

subplot(122)

xout_all_APF = cell(length(rho_values) ,1 ) ; 

for i = 1:length(rho_values)
    
    rho0 =rho_values(i); 
    counter = 1 ; 
    x = p0;  %Agent state
    goal_dist = norm(x-pg,2);
    xout = [] ; 
    %Control Loop

    while goal_dist >= goal_tol
         u_des_APF = - gradUatt(x) ;   % desired input
       
        if (rho(x,pO,rO) < rho0)
            h_APF = 1/(1+Urep(x,pO,rO)) - sigma0 ; 
            gradh_APF = - gradUrep(x,pO,rO) / ((1+Urep(x,pO,rO))^2) ; 
        else
            h_APF = 1-sigma0 ; 
            gradh_APF = zeros(2,1) ; 
        end

        phi_APF =  gradh_APF' * u_des_APF + alpha_1*h_APF ; 
        
        if (phi_APF >= 0 )
            u=u_des_APF ; 
        else
            u_safe_APF= - (gradh_APF) / (norm(gradh_APF,2)^2) * phi_APF ; 
            u = u_safe_APF + u_des_APF ; 
        end    

        %Trajectory 
        x = x + Ts * g(x) * [u(1), u(2)]';  % system (Explicit Euler) 
        xout(: , counter) = x ; 
        
        counter=counter +1 ; 
    
        goal_dist = norm(x-pg,2);  % update goal dist
        
        %{
        clf 
        %Scenarion
        plot(p0(1),p0(2),'bo','LineWidth',1.2), hold on,
        plot(pg(1),pg(2),'rx','LineWidth',1.2), hold on, 
        %plot(xO,yO,'k-','LineWidth',1.2 ), hold on,
        fill(xO,yO , 'black')
        scatter(x(1),x(2),'blue','filled','MarkerEdgeColor','k')
        axis equal
        
        %Wait for sample time
        waitfor(rate);
        %}

    end
   xout_all_APF{i}=xout ; 
end
 hold on ; 
plot(p0(1),p0(2),'bx','LineWidth',1.2)
plot(pg(1),pg(2),'rx','LineWidth',1.2) 
fill(xO,yO, 'black' , 'Facealpha' , 0.5)
%scatter(x(1),x(2),'blue','filled','MarkerEdgeColor','k')
axis equal

x = linspace (p0(1) , pg(1) ) ; 
y = linspace (p0(2) , pg(2) ) ; 
plot (x,y  ,'LineStyle','--' , 'LineWidth', 1.1)

colors = ['r' , 'g' , 'b' , 'magenta'] ; 

for i = 1:length(rho_values)
    trajectory = xout_all_APF{i};
    plot(trajectory(1, :), trajectory(2, :),'Color',colors(i),'LineWidth', 1.1);
end

for i= 1: length(rho_values)
legend_rho{i} =  sprintf("rho = %0.1f" , rho_values(i)) ; 
end

legend ("Starting Point","Goal Point" , "Obstacle "  , "Desired trjectory", legend_rho{:},  'Location','northwest')
axis equal ;
xlabel("X[m]");
ylabel("Y[m") ; 
