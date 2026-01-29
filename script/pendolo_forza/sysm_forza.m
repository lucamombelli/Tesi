clearvars;
close all;
clc;
%% PARAMETRI E VARIABILI
syms x v_x
syms theta v_theta
syms u
syms m l M g
syms t

x_target_pos = 5;

% Parametri per l'ostacolo 
d = 1; % Altezza 
epsilon = 0.1 ; %Raggio dell'ostacolo

% Boundend Control
u_max = 300;
u_min = -300;

parametri_simbolici = [M, m, l, g];
valori_numerici     = [1, 0.1, 1, 9.81];
state = [x, theta, v_x, v_theta].';

%% DINAMICA
A_mat = [1, 0, 0, 0;
    0, 1, 0, 0;
    0, 0, (M+m),             -m*l*cos(theta);
    0, 0, -m*l*cos(theta),    m*l^2];

b_vec = [v_x;
    v_theta;
    -m*l*v_theta^2 * sin(theta)+u;
    m*g*l*sin(theta)];

sys = A_mat \ b_vec;
f_sym = subs(sys, u, 0);
g_sym = jacobian(sys, u);

f_num_sym = subs(f_sym, parametri_simbolici, valori_numerici);
g_num_sym = subs(g_sym, parametri_simbolici, valori_numerici);

f_fun = matlabFunction(f_num_sym, 'Vars', {state});
g_fun = matlabFunction(g_num_sym, 'Vars', {state});

sys_num = subs(sys, parametri_simbolici, valori_numerici);
dynamics_fun = matlabFunction(sys_num, 'Vars', {t, state, u});

%% LQR

A_lin = double(subs(jacobian(f_sym, state), [state; parametri_simbolici.'], [zeros(4,1); valori_numerici.']));
B_lin = double(subs(g_sym, [state; parametri_simbolici.'], [zeros(4,1); valori_numerici.']));

% PESI LQR
% Q penalizza l'errore rispetto al target
%Q = diag([100, 500, 1, 1]);

Q = zeros(4,4);
Q(1,1) = 50;
Q(2,2) = 100 ;
Q(3,3) = 1 ;
Q(4,4)  = 1 ;

R = 50; % R rappresenta la matrice di peso sull'azione del controllo 
% - R grande allora e piu costoso usare il controllo 
% - R piccolo , il controllo diventa piu aggressivo 

[K_lqr , P_lqr , ~ ] = lqr(A_lin, B_lin, Q, R);
%% 4. BARRIERE (CBF)
%barriera h per l'ostacolo
k1 =25;

r2 =x^2 +l^2 +d^2 - 2*x*l*sin(theta)-2*d*l*cos(theta);
r2_num = subs(r2 , parametri_simbolici , valori_numerici) ;
h0_sym = r2_num - epsilon^2 ;
grad_h0_sym = gradient(h0_sym, state);
h1_sym = dot(grad_h0_sym,f_num_sym)+k1*h0_sym ;
grad_h1_sym = gradient(h1_sym,state);

h1_fun = matlabFunction(h1_sym, 'Vars', {state});
grad_h1_fun = matlabFunction(grad_h1_sym, 'Vars', {state});


%Barriera B per l'angolo critico
a1 = 10; %a1 =80

theta_critico = atan(u_max / (9.81 * (valori_numerici(1) + valori_numerici(2))));
%b0_sym = theta_critico- sqrt(theta^2 + 1e-6);
b0_sym = theta_critico^2 - theta^2;
grad_b0_sym = gradient(b0_sym, state);
Lf_b0_sym = grad_b0_sym.' * f_num_sym;
b1_sym = Lf_b0_sym + a1 * b0_sym;
grad_b1_sym = gradient(b1_sym, state);

b1_fun = matlabFunction(b1_sym, 'Vars', {state});
grad_b1_fun = matlabFunction(grad_b1_sym, 'Vars', {state});

%% BARRIERA E ANGOLO

y0 = [ - 2; 0 ; 0 ; 0] ;
tstar = 15 ;

%Parametri CBF
a2 = 20; 
k2 = 50; % lower bound = 5 

[tout , yout , u_out ] = cbf_obs_angle(y0, f_fun, g_fun, b1_fun, grad_b1_fun, h1_fun , grad_h1_fun , k2,a2, u_min, u_max, K_lqr, x_target_pos,tstar);


plot_pendolo(tout , yout ,u_out ,u_max,u_min,x_target_pos ,theta_critico , valori_numerici , epsilon ,d )

%% ANIMAZIONE PENDOLO

fig_anim = figure('Name', 'Animazione Pendolo', 'NumberTitle', 'off', ...
                  'Position', [100, 100, 800, 600]);


i_end=size(yout,1);
%j_target = find(yout(:,1)>= x_target_pos ,1) ;
for i_frame = 1:10:i_end
    %figure(fig_anim);
    set(0, 'CurrentFigure', fig_anim);
    clf;
    cart_x   = yout(i_frame ,1 );  % Cart position
    pend_ang = yout(i_frame ,2 );  % Pendulum angle
    IP_Animation(cart_x, pend_ang , d , epsilon,valori_numerici , x_target_pos); 
    pause(0.001);
end

%% ONLY OBSTACLE
%{
k2_clf = 1 ; 
a2_clf = 12 ; 
tstar = 4 ; 
[t_obs , y_obs , u_out ] = cbf_obs(y0, f_fun, g_fun, b1_fun, grad_b1_fun, h1_fun , grad_h1_fun , k2_clf,a2_clf, u_min, u_max, K_lqr, x_target_pos,tstar,V_fun,grad_V_fun);


plot_pendolo(t_obs , y_obs ,u_out ,u_max,u_min,x_target_pos ,theta_critico , valori_numerici , epsilon ,d)


%% Controlo Lyapunov Function


target_state = [x_target_pos ; 0 ; 0 ; 0];

err = state - target_state;
V_sym = 0.5 * err.' * P_lqr * err;

% Calcola i gradienti simbolici
grad_V_sym = gradient(V_sym, state);
% Crea le funzioni Matlab
V_fun = matlabFunction(V_sym, 'Vars', {state});
grad_V_fun = matlabFunction(grad_V_sym, 'Vars', {state});

%% PROVA UTILIZZANDO SOLO LE CONTROL LYAPUNOV FUNCTION 

k2_clf = 1 ; 
a2_clf = 12 ; 
tstar = 4 ; 
[t_obs , y_obs , u_out] = cbf_obs(y0, f_fun, g_fun, b1_fun, grad_b1_fun, h1_fun , grad_h1_fun , k2_clf,a2_clf, u_min, u_max, K_lqr, x_target_pos,tstar,V_fun,grad_V_fun);


plot_pendolo(t_obs , y_obs ,u_out ,u_max,u_min,x_target_pos ,theta_critico , valori_numerici , epsilon ,d)

%}
%% SAFE SET 


figure;
% Definizione della funzione anonima h0(x, theta)
% Parametri: l = 5; d = 3; epsilon = 0.1;
% Formula: r^2 - epsilon^2
mia = @(x, theta) x.^2 + 4^2 + 3^2 - 2*4*x.*sin(theta) - 2*4*3*cos(theta) - 0.01;

% Plot della curva di livello 0 (il bordo dell'ostacolo nello spazio degli stati)
fimplicit(mia, [-5 5 -pi/2 pi/2], 'LineWidth', 2, 'Color', 'r');

grid on;
xlabel('Posizione Carrello x [m]');
ylabel('Angolo Pendolo \theta [rad]');
title('Bordo del Safe Set Geometrico (h_0 = 0)');
legend('Confine Ostacolo');

% Aggiungi una linea per vedere dove sei tu adesso (es. x=-2, theta=0)
%hold on;
%plot(-2, 0, 'bo', 'MarkerFaceColor', 'b');
%text(-2, 0.2, ' Start');
