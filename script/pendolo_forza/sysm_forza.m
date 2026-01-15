clearvars;
close all;
clc;
%% PARAMETRI E VARIABILI
syms x v_x
syms theta v_theta
syms u
syms m l M g
syms t

x_target_pos = 3;


d = 0.5; % Altezza dell'ostacolo
u_max = 100;
u_min = -100;


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

y_vec = [x; theta; v_x; v_theta];
f_fun = matlabFunction(f_num_sym, 'Vars', {y_vec});
g_fun = matlabFunction(g_num_sym, 'Vars', {y_vec});

sys_num = subs(sys, parametri_simbolici, valori_numerici);
dynamics_fun = matlabFunction(sys_num, 'Vars', {t, y_vec, u});

%% LQR

A_lin = double(subs(jacobian(f_sym, state), [state; parametri_simbolici.'], [zeros(4,1); valori_numerici.']));
B_lin = double(subs(g_sym, [state; parametri_simbolici.'], [zeros(4,1); valori_numerici.']));

% PESI LQR
% Q penalizza l'errore rispetto al target
%Q = diag([100, 500, 1, 1]);

Q = zeros(4,4);
Q(1,1) = 1000;
Q(2,2) = 10 ;
Q(3,3) =  10 ;
Q(4,4)  = 10 ;

R = 50;

K_lqr = lqr(A_lin, B_lin, Q, R);
%% 4. BARRIERE (CBF)
%barriera h per l'ostacolo
epsilon = 0.3 ;
k1 = 1;
k2 = 2;

r2 =x^2 +l^2 +d^2 - 2*x*l*sin(theta)-2*d*l*cos(theta);
r2_num = subs(r2 , parametri_simbolici , valori_numerici) ;
h0_sym = r2_num - epsilon ;
grad_h0_sym = gradient(h0_sym, state);
h1_sym = dot(grad_h0_sym,f_num_sym)+k1*h0_sym ;
grad_h1_sym = gradient(h1_sym,state);

h1_fun = matlabFunction(h1_sym, 'Vars', {y_vec});
grad_h1_fun = matlabFunction(grad_h1_sym, 'Vars', {y_vec});


%Barriera B per l'angolo critico
a1 = 10;
a2 = 10;

theta_critico = atan(u_max / (9.81 * (valori_numerici(1) + valori_numerici(2))));

b0_sym = theta_critico- sqrt(theta^2 + 1e-2);
b0_sym = theta_critico^2 - theta^2;
grad_b0_sym = gradient(b0_sym, state);
Lf_b0_sym = grad_b0_sym.' * f_num_sym;
b1_sym = Lf_b0_sym + a1 * b0_sym;
grad_b1_sym = gradient(b1_sym, state);

b1_fun = matlabFunction(b1_sym, 'Vars', {y_vec});
grad_b1_fun = matlabFunction(grad_b1_sym, 'Vars', {y_vec});


%% Controlo Lyapunov Function


target_state = [x_target_pos ; 0 ; 0 ; 0];

V = 0.5*norm(state - target_state,2)^2 ;
V_fun = matlabFunction(V,'Vars',{y_vec});
grad_V = state-target_state;
grad_V_fun = matlabFunction(grad_V , 'Vars' , {y_vec});


%% SOLO ANGOLO

%condizioni inziali
y0_angle = [ -2 ; -0.5 ; 0 ; 0] ;

[tout_angle , yout_angle , u_angle , u_nom ] = cbf_angle(y0_angle, f_fun, g_fun, b1_fun, grad_b1_fun,a2, u_min, u_max, K_lqr, x_target_pos,1);

% PLOT
figure();

subplot(4,1,1);
plot(tout_angle, yout_angle(:,1), 'LineWidth', 2);
hold on;
yline(x_target_pos, 'g--', 'Target');
xline(0,'w--','Obstacle')
grid on;
ylabel('Posizione x [m]'); title(['Posizione (Target: ' num2str(x_target_pos) 'm)']);

subplot(4,1,2);
plot(tout_angle, yout_angle(:,2), 'r', 'LineWidth', 2);
hold on;
yline(double(theta_critico), 'k--');
yline(-double(theta_critico), 'k--');
grid on; ylabel('Theta [rad]'); title('Angolo Pendolo');

subplot(4,1,3);
hold on,
plot(tout_angle, u_angle, 'b', 'LineWidth', 1.5);
plot(tout_angle ,u_nom)
yline(u_max, 'r--'); yline(u_min, 'r--');
grid on; ylabel('Controllo u [N]'); xlabel('Tempo [s]');
hold off ;

subplot(4,1,4)
y_angle = cos(yout_angle(:,2)) ;
plot(tout_angle,y_angle , 'LineWidth',1.4)

%% ANIMAZIONE PENDOLO ( SOLO ANGOLO ) 
figure; 

i_angle = size(yout_angle,1 );
for i_frame = 1:10:i_angle
    clf;
    cart_x   = yout_angle(i_frame ,1 );  % Cart position
    pend_ang = yout_angle(i_frame ,2 );  % Pendulum angle
    IP_Animation_angle(cart_x, pend_ang , valori_numerici); 
    pause(0.02);
end


%% BARRIERA E ANGOLO

y0 = [ - 2 ; 0 ; 0 ; 0] ;

[tout , yout , u_out , count] = cbf_obs_angle(y0, f_fun, g_fun, b1_fun, grad_b1_fun, h1_fun , grad_h1_fun , k2,a2, u_min, u_max, K_lqr, x_target_pos,30,V_fun,grad_V_fun);

if(count ~= 0)
    fprintf('Problem is infeasible  %i times \n ' , count),
end

plot_pendolo(tout , yout ,u_out ,u_max,u_min,x_target_pos ,theta_critico , valori_numerici , epsilon ,d )

%% Calcolo distanza dall'ostacolo
x_punta = yout(:,1) - valori_numerici(3)*sin(yout(:,2));
y_punta = valori_numerici(3)*cos(yout(:,2));

% Distanza dal centro dell'ostacolo (0, d)
dist_from_obstacle = sqrt((x_punta).^2 + (y_punta - d).^2);
clearance = dist_from_obstacle - epsilon;
% Trova violazioni
violations = find(clearance < 0);
if ~isempty(violations)
    fprintf('VIOLAZIONE OSTACOLO rilevata!\n');
    fprintf('   Prima violazione a t=%.3f s\n', tout(violations(1)));
    fprintf('   Penetrazione massima: %.4f m\n', min(clearance));
else
    fprintf('✓ Nessuna violazione ostacolo\n');
    fprintf('  Margine minimo: %.4f m\n', min(clearance));
end

%% ANIMAZIONE PENDOLO
figure
i_end = size(yout,1 );
j_target = find(yout(:,1)>= x_target_pos ,1) ;
for i_frame = 1:10:j_target+200
    clf;
    cart_x   = yout(i_frame ,1 );  % Cart position
    pend_ang = yout(i_frame ,2 );  % Pendulum angle
    IP_Animation(cart_x, pend_ang , d , epsilon,valori_numerici , x_target_pos); 
    pause(0.02);
end

%% PROVA STABILIZZAZIONE

%[t_stab , y_stab , u_out , count] = cbf_then_stabilize(y0, f_fun, g_fun, b1_fun, grad_b1_fun, h1_fun , grad_h1_fun , k2,a2, u_min, u_max, K_lqr, x_target_pos,4,V_fun,grad_V_fun);


%plot_pendolo(t_stab , y_stab ,u_out ,u_max,u_min,x_target_pos ,theta_critico , valori_numerici , epsilon ,d)
