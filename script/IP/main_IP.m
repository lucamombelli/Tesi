%% Inverse Pendulum on a Cart Limbo

clearvars;
close all;
clc;

%% Functions for simulation
function [sysfunStruct,state,sym_par,sys_par] = system_dynamics()

    % Variables and parameters
    syms x1 x2 x3 x4 real % (x,theta,dx,dtheta)
    syms u real % control input
    syms m l M g mu_d real % (ip mass, rod length, cart mass, gravity, friction coeff)
    syms t real % time
    
    % State vector & parameters
    sym_par = [M, m, l, g, mu_d];
    sys_par = [1, 0.1, 1, 9.81, 0.1];
    state = [x1, x2, x3, x4].';

    % Dynamics
    A_mat = [1, 0, 0, 0;
             0, 1, 0, 0;
             0, 0, (M+m), -m*l*cos(x2);
             0, 0, -m*l*cos(x2), m*l^2];
    
    b_vec = [x3;
             x4;
             -m*l*x4^2 * sin(x2)+u-mu_d*x3 ;
             m*g*l*sin(x2)];

    % Symbolic system 
    sys = A_mat \ b_vec;  % symbolic system
    f_sym = subs(sys, u, 0);  % symbolic drift
    g_sym = jacobian(sys, u);  % symbolic input map
    
    f_num_sym = subs(f_sym, sym_par, sys_par);
    g_num_sym = subs(g_sym, sym_par, sys_par);
    
    % Conversion to "function handle"
    f_fun = matlabFunction(f_num_sym, 'Vars', {state});  % drift
    g_fun = matlabFunction(g_num_sym, 'Vars', {state});  % input map
    
    % Output structure
    sysfunStruct.system = sys;
    sysfunStruct.symb_drift = f_sym;
    sysfunStruct.symb_input = g_sym;
    sysfunStruct.num_sym_drift = f_num_sym;
    sysfunStruct.num_sym_input = g_num_sym;
    sysfunStruct.fun_drift = f_fun;
    sysfunStruct.fun_input = g_fun;
end

function [A,B] = linearized_dynamics(f_sym,g_sym,state,sym_par,sys_par)
    
    % Jacobian w.r.t the state
    A = double(subs(jacobian(f_sym, state), [state; sym_par.'], [zeros(4,1); sys_par.']));

    % Jacpbian w.r.t the input
    B = double(subs(g_sym, [state; sym_par.'], [zeros(4,1); sys_par.']));
end

function [cbfObsStruct,cbfThetaStruct] = control_barrier_functions(state,obstacle,sym_par,sys_par,alpha_1,gamma_1,theta_crit,sysfunStruct)
    
    % Unpack the state (symbolic)
    x1 = state(1);
    x2 = state(2);
    x3 = state(3);
    x4 = state(4);

    % Unpack the parameters
    l = sym_par(3);
    f_num_sym = sysfunStruct.num_sym_drift;
    
    % Squared distance system - obstacle center
    r2 = x1^2 + l^2 + obstacle.height^2 - 2*x1*l*sin(x2) - 2*obstacle.height*l*cos(x2);
    r2_num = subs(r2,sym_par,sys_par);
    
    % FO-CBF for obstacle avoidance
    h0_sym = r2_num - obstacle.radius^2;  % symbolic constraint
    grad_h0_sym = gradient(h0_sym,state);  % symbolic gradient

    h1_sym = dot(grad_h0_sym,f_num_sym) + gamma_1*h0_sym ;  % symbolic CBF constraint
    grad_h1_sym = gradient(h1_sym,state);  % symbolic CBF gradient
    
    h1_fun = matlabFunction(h1_sym,'Vars',{state});  % CBF constraint (function handle)
    grad_h1_fun = matlabFunction(grad_h1_sym,'Vars',{state});  % CBF gradient (function handle)

    % FO-CBF for critic angle
    b0_sym = theta_crit^2 - x2^2;
    grad_b0_sym = gradient(b0_sym,state);

    Lf_b0_sym = grad_b0_sym.'*f_num_sym;
    b1_sym = Lf_b0_sym + alpha_1*b0_sym;
    grad_b1_sym = gradient(b1_sym,state);
    
    b1_fun = matlabFunction(b1_sym, 'Vars', {state});
    grad_b1_fun = matlabFunction(grad_b1_sym, 'Vars', {state});

    % Output structures
    cbfObsStruct.obs_cons = h1_fun;
    cbfObsStruct.obs_cons_gradient = grad_h1_fun;
    cbfThetaStruct.theta_cons = b1_fun;
    cbfThetaStruct.theta_cons_gradient = grad_b1_fun;
end

%% Obstacle 
% d = 1 works with gamma_1 = gamma_2 = 8-25 and epsilon = 0.1 - 0.49
% d = 0.9 works with gamma_1 = gamma_2 = 10-   and epsilon = 
obs.height = 1; % height
obs.radius = 10^(-1); % radius

%% System dynamics 
u_max = inf;  % lower bound
u_min = -inf;  % upper bound

[sysfunStruct,state,sym_par,sys_par] = system_dynamics();

sys = sysfunStruct.system;  % symbolic system
f_sym = sysfunStruct.symb_drift;  % symbolic drift
g_sym = sysfunStruct.symb_input;  % symbolic input 
f_num_sym = sysfunStruct.num_sym_drift;
g_num_sym = sysfunStruct.num_sym_input; 
f_fun = sysfunStruct.fun_drift;  % drift (function handle)
g_fun = sysfunStruct.fun_input;  % input (function handle)

%% LQR control 
% Q penalizes the error with respect to the target (default Q = diag([100, 500, 1, 1])) 
% R represents the weight matrix on the control action 
% - R large then and more expensive to use the control 
% - R small, control becomes more aggressive
[A,B] = linearized_dynamics(f_sym,g_sym,state,sym_par,sys_par);
Q = diag([10, 100, 1, 1]);
R = 1.0;
[K_lqr,P_lqr,~] = lqr(A,B,Q,R);

%% CBFs - obstacle and critic angle
theta_crit = pi/2;  % do be estimated
alpha_1 = 60;  % class-K function (slope of a straight line)
gamma_1 = 60;  % class-K function (slope of a straight line)

[cbfObsStruct,cbfThetaStruct] = control_barrier_functions(state,obs,sym_par,sys_par,alpha_1,gamma_1,theta_crit,sysfunStruct);

h1_fun = cbfObsStruct.obs_cons;  % constraint function (obstacle)
grad_h1_fun = cbfObsStruct.obs_cons_gradient;  % constraint gradient (obstacle)
b1_fun = cbfThetaStruct.theta_cons;  % constraint function (critic theta)
grad_b1_fun = cbfThetaStruct.theta_cons_gradient;  % constraint gradient (critic theta)

%% Simulation loop
x0 = -3;  % initial cart position
xg = 4;  % target cart position
X0 = [x0, 0, 0, 0]';  % initial state
Tf = 8;  % final time [s]

% Parameters of SO-CBFs
alpha_2 = 60; % class-K function (slope of a straight line)
gamma_2 = 60 ; % (lower bound = 8) class-K function (slope of a straight line)
[tout_rk,xout_rk,u_out_rk] = IP_simulation_rk(X0, f_fun, g_fun, b1_fun, grad_b1_fun, ...
                                  h1_fun, grad_h1_fun, gamma_2, alpha_2, u_min, u_max, K_lqr, xg, Tf);

%% Plots
plot_IP(tout_rk, xout_rk, u_out_rk, u_max, u_min, xg, theta_crit, sys_par, obs)
%plot_safe_set(sys_par, obs, xg, X0)
%phase_portrait(xout_rk)

%{
%% CBFs - obstacle and critic angle Euler 
theta_crit = pi/2;  % do be estimated
alpha_1 = 60;  % class-K function (slope of a straight line)
gamma_1 = 35;  % class-K function (slope of a straight line)

[cbfObsStruct,cbfThetaStruct] = control_barrier_functions(state,obs,sym_par,sys_par,alpha_1,gamma_1,theta_crit,sysfunStruct);

h1_fun = cbfObsStruct.obs_cons;  % constraint function (obstacle)
grad_h1_fun = cbfObsStruct.obs_cons_gradient;  % constraint gradient (obstacle)
b1_fun = cbfThetaStruct.theta_cons;  % constraint function (critic theta)
grad_b1_fun = cbfThetaStruct.theta_cons_gradient;  % constraint gradient (critic theta)

%% Simulation loop
x0 = -3;  % initial cart position
xg = 4;  % target cart position
X0 = [x0, 0, 0, 0]';  % initial state
Tf = 20;  % final time [s]

% Parameters of SO-CBFs
alpha_2 = 60; % class-K function (slope of a straight line)
gamma_2 = 45 ; % (lower bound = 8) class-K function (slope of a straight line)
[tout,xout,u_out] = IP_simulation(X0, f_fun, g_fun, b1_fun, grad_b1_fun, ...
                                  h1_fun, grad_h1_fun, gamma_2, alpha_2, u_min, u_max, K_lqr, xg, Tf);

%% Plots
plot_IP(tout, xout, u_out, u_max, u_min, xg, theta_crit, sys_par, obs)
%plot_safe_set(sys_par, obs, xg, X0)
%phase_portrait(xout_rk)
%}
%% Pendulum animation
sim = 0;
if(sim == 1)
    fig_anim = figure('Name','Pendulum Animation','NumberTitle','off', ...
                  'Position',[100,100,800,600]);

    i_end = size(xout_rk,1);
    %j_target = find(xout(:,1)>= xg ,1) ;
    for i_frame = 1:10:i_end
        % figure(fig_anim);
        set(0,'CurrentFigure',fig_anim);
        clf;
        cart_x = xout_rk(i_frame ,1 );  % Cart position
        pend_ang = xout_rk(i_frame ,2 );  % Pendulum angle
        IP_animation(cart_x, pend_ang , obs ,sys_par , xg , X0); 
        pause(0.01);
    end
end

plot_safe_set(sys_par, obs, xg, X0)



