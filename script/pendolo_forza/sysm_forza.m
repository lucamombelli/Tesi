clear all;
close all;
clc;

%% PARAMETRI E VARIABILI
syms x v_x
syms theta v_theta
syms u
syms m l M g
syms t

x_target_pos = 3;


d = 0.8; % Altezza dell'ostacolo
u_max = 400;
u_min = -400;


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
%C = ctrb(A_lin , B_lin) ;
%Q = diag([1, 100, 1, 1]);
%Q = C'*C;
Q = zeros(4,4);
Q(1,1) = 1000;
Q(2,2) = 10 ;
Q(3,3) =  10 ;
Q(4,4)  = 10 ;

R = 100;

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

[tout_angle , yout_angle , u_angle , u_nom ] = cbf_angle(y0_angle, f_fun, g_fun, b1_fun, grad_b1_fun,a2, u_min, u_max, K_lqr, x_target_pos,2);

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



%% BARRIERA E ANGOLO

y0 = [ - 1 ; 0.3 ; 0 ; 0] ;

[tout , yout , u_out , count] = cbf_obs_angle(y0, f_fun, g_fun, b1_fun, grad_b1_fun, h1_fun , grad_h1_fun , k2,a2, u_min, u_max, K_lqr, x_target_pos,4,V_fun,grad_V_fun);

if(count ~= 0)
    disp('Problem Infeasible at some instant'),
end

% Plot
j_obs = find(yout(:,1)>= 0 , 1) ; % inidice per cui trovo l'ostacolo
j_target = find(yout(:,1)>= x_target_pos ,1) ; % indice per cui ho superato x_target_pos

%Figure di tutta la simulazione
figure();
subplot(3,1,1);
plot(tout(), yout(:,1), 'LineWidth', 2);
hold on;
yline(x_target_pos, 'g--', 'Target');
grid on;
ylabel('x[m]'); title(['Posizione (Target: ' num2str(x_target_pos) 'm)']);
hold off ;

subplot(3,1,2);
plot(tout, yout(:,2), 'r', 'LineWidth', 2);
hold on;
yline(double(theta_critico), 'k--');
yline(-double(theta_critico), 'k--');
grid on; ylabel('Theta [rad]'); title('Angolo Pendolo');

subplot(3,1,3);
plot(tout, u_out, 'b', 'LineWidth', 1.5);
yline(u_max, 'r--'); yline(u_min, 'r--');
grid on; ylabel('Controllo u [N]'); xlabel('Tempo [s]');
title('Controllore Ottimo')
mesh = linspace(0,2*pi , 1000);
%{
subplot(4,1,4)
y = valori_numerici(3)*cos(yout(:,2)) ;
hold on ;
plot(tout,y , 'LineWidth',1)
scatter(tout(j_obs),y(j_obs),26)
plot(cos(mesh)*epsilon + tout(j_obs) , sin(mesh)*epsilon+d,'LineWidth',1,'Color','r')
axis equal ;
title('Posizione y')
hold off ;
%}

% ---------------------------------------------------------------------------------------------------
%Figure che si fa quando (e se) il pendolo passa per il target
if(j_target)
    figure();
    subplot(4,1,1);
    plot(tout(1:j_target), yout(1:j_target,1), 'LineWidth', 2);
    hold on;
    scatter(tout(j_obs) , yout(j_obs,1),26)
    yline(x_target_pos, 'g--', 'Target');
    grid on;
    ylabel('x[m]'); title(['Posizione (Target: ' num2str(x_target_pos) 'm)']);
    hold off ;

    subplot(4,1,2);
    plot(tout(1:j_target), yout(1:j_target,2), 'r', 'LineWidth', 2);
    hold on;
    yline(double(theta_critico), 'k--');
    yline(-double(theta_critico), 'k--');
    grid on; ylabel('Theta [rad]'); title('Angolo Pendolo');

    subplot(4,1,3);
    plot(tout(1:j_target), u_out(1:j_target), 'b', 'LineWidth', 1.5);
    yline(u_max, 'r--'); yline(u_min, 'r--');
    grid on; ylabel('Controllo u [N]'); xlabel('Tempo [s]');
    title('Controllore Ottimo')
    mesh = linspace(0,2*pi , 1000);

    subplot(4,1,4)
    %Sto plottando l'altezza del pendolo e anche l'ostacolo , da
    %controllare non so sicuro che sia giusto
    height = valori_numerici(3)*cos(yout(1:j_target,2)) ;
    hold on ;
    plot(tout(1:j_target),height , 'LineWidth',1)
    scatter(tout(j_obs), height(j_obs),26)
    plot(cos(mesh)*epsilon + tout(j_obs) , sin(mesh)*epsilon+d,'LineWidth',1,'Color','r')
    axis equal ;
    title('Posizione y')
    hold off ;
else
    disp('Il pendolo non ha raggiunto la target position entro il tempo della simulazione')
end

%------------------------------------------
% PLOT SPAZIALE
figure;

% Posizione punta pendolo nello spazio
x_punta = yout(:,1) - valori_numerici(3)*sin(yout(:,2));
y_punta = valori_numerici(3)*cos(yout(:,2));

% Plot traiettoria
plot(x_punta, y_punta, 'b-', 'LineWidth', 2);
hold on;

% Ostacolo in (0, d)
theta_circle = linspace(0, 2*pi, 100);
x_circle = epsilon * cos(theta_circle);
y_circle = epsilon * sin(theta_circle) + d;
plot(x_circle, y_circle, 'r-', 'LineWidth', 2);
fill(x_circle, y_circle, 'r', 'FaceAlpha', 0.3);
plot(0, d, 'rx', 'MarkerSize', 10, 'LineWidth', 3);

% Traiettoria del carrello (sulla rotaia a y=0)
plot(yout(:,1), zeros(size(yout(:,1))), 'k--', 'LineWidth', 1.5);

% Punti critici
if ~isempty(j_obs)
    scatter(x_punta(j_obs), y_punta(j_obs), 100, 'g', 'filled');
    text(x_punta(j_obs), y_punta(j_obs), ' ← x=0', 'FontSize', 12);
end

axis equal;
grid on;
xlabel('Posizione x [m]');
ylabel('Altezza y [m]');
title('Vista Spaziale: Traiettoria Punta Pendolo');
legend('Punta pendolo', 'Ostacolo', '', 'Centro ostacolo', 'Carrello', 'Passaggio x=0');
xlim([-2, x_target_pos]);
ylim([-0.5, valori_numerici(3)+0.3]);
hold off;


%% Calcolo distanza dall'ostacolo
x_punta = yout(:,1) - valori_numerici(3)*sin(yout(:,2));
y_punta = valori_numerici(3)*cos(yout(:,2));

% Distanza dal centro dell'ostacolo (0, d)
dist_from_obstacle = sqrt((x_punta).^2 + (y_punta - d).^2);
clearance = dist_from_obstacle - epsilon;
% Trova violazioni
violations = find(clearance < 0);
if ~isempty(violations)
    fprintf('⚠️ VIOLAZIONE OSTACOLO rilevata!\n');
    fprintf('   Prima violazione a t=%.3f s\n', tout(violations(1)));
    fprintf('   Penetrazione massima: %.4f m\n', min(clearance));
else
    fprintf('✓ Nessuna violazione ostacolo\n');
    fprintf('  Margine minimo: %.4f m\n', min(clearance));
end

%% PROVA STABILIZZAZIONE

[t_stab , y_stab , u_out , count] = cbf_then_stabilize(y0, f_fun, g_fun, b1_fun, grad_b1_fun, h1_fun , grad_h1_fun , k2,a2, u_min, u_max, K_lqr, x_target_pos,4,V_fun,grad_V_fun);
%Figure di tutta la simulazione
figure();
subplot(2,1,1);
plot(t_stab, y_stab(:,1), 'LineWidth', 2);
hold on;
yline(x_target_pos, 'g--', 'Target');
grid on;
ylabel('x[m]'); title(['Posizione (Target: ' num2str(x_target_pos) 'm)']);
hold off ;

subplot(2,1,2);
plot(t_stab, y_stab(:,2), 'r', 'LineWidth', 2);
hold on;
yline(double(theta_critico), 'k--');
yline(-double(theta_critico), 'k--');
grid on; ylabel('Theta [rad]'); title('Angolo Pendolo');


% PLOT SPAZIALE
figure;

% Posizione punta pendolo nello spazio
x_punta_stab = y_stab(:,1) - valori_numerici(3)*sin(y_stab(:,2));
y_punta_stab = valori_numerici(3)*cos(y_stab(:,2));

% Plot traiettoria
plot(x_punta_stab, y_punta_stab, 'b-', 'LineWidth', 2);
hold on;

% Ostacolo in (0, d)
theta_circle = linspace(0, 2*pi, 100);
x_circle = epsilon * cos(theta_circle);
y_circle = epsilon * sin(theta_circle) + d;
plot(x_circle, y_circle, 'r-', 'LineWidth', 2);
fill(x_circle, y_circle, 'r', 'FaceAlpha', 0.3);
plot(0, d, 'rx', 'MarkerSize', 10, 'LineWidth', 3);

% Traiettoria del carrello (sulla rotaia a y=0)
plot(y_stab(:,1), zeros(size(y_stab(:,1))), 'k--', 'LineWidth', 1.5);

% Punti critici
if ~isempty(j_obs)
    scatter(x_punta_stab(j_obs), y_punta_stab(j_obs), 100, 'g', 'filled');
    text(x_punta_stab(j_obs), y_punta_stab(j_obs), ' ← x=0', 'FontSize', 12);
end

axis equal;
grid on;
xlabel('Posizione x [m]');
ylabel('Altezza y [m]');
title('Vista Spaziale: Traiettoria Punta Pendolo');
legend('Punta pendolo', 'Ostacolo', '', 'Centro ostacolo', 'Carrello', 'Passaggio x=0');
xlim([-2, x_target_pos]);
ylim([-0.5, valori_numerici(3)+0.3]);
hold off;
%% Funzioni Ausiliare (Pulizia del codice )
function [tout , yout , u_applied , u_nom] = cbf_angle(y, f_fun, g_fun, b1_fun, grad_b1_fun,a2, u_min, u_max, K, x_tgt,tstar)

dt = 0.001;
tspan = 0:dt:tstar ;
yout = zeros(4,size(tspan,2)) ;
u_applied  = zeros(1,size(tspan,2))  ;
u_nom  = zeros(1,size(tspan,2))  ;
n = 1 ;
yout(:,n) = y ;
options = optimoptions('quadprog', 'Display', 'off');

for n=1:length(tspan)-1
    x_curr = yout(:,n);
    f_val = f_fun(yout(:,n));
    g_val = g_fun(yout(:,n));

    grad_b1 = grad_b1_fun(yout(:,n));
    val_b1  = b1_fun(yout(:,n));

    Lf_b1 = dot(grad_b1, f_val);
    Lg_b1 = dot(grad_b1, g_val);
    %val = yout(:,n) ;
    %u_nom(1,n) = (x_tgt-val(1)) ;
    %u_nom(1,n) = 1.2 ;

    y_des = [x_tgt;0;0;0];
    u_nom = - K * (x_curr - y_des) ;


    H = 1;
    f_qp = -u_nom;
    size(f_qp);
    A = -Lg_b1;
    b = Lf_b1 + a2 * val_b1;
    lb = u_min;
    ub = u_max;

    [z , ~ , flag] = quadprog(H, f_qp, A, b, [], [], lb, ub, [], options);

    if flag == -2
        error('Problem is infeasible')
    else
        u_applied(:,n) = z ;
    end

    yout(:,n+1) = yout(:,n)+ dt * (f_val + g_val*u_applied(:,n)) ;
end
tout = tspan ;
yout = yout.';
end


function [tout , yout , u_applied , counter_inf] = cbf_obs_angle(y, f_fun, g_fun, b1_fun, grad_b1_fun,h1_fun , grad_h1_fun , k2 ,a2, u_min, u_max, K, x_tgt,tstar,V_fun,grad_V_fun)

dt = 0.001;
tspan = 0:dt:tstar ;
yout = zeros(4,size(tspan,2)) ;
u_applied  = zeros(1,size(tspan,2))  ;
n = 1 ;
yout(:,n) = y ;
options = optimoptions('quadprog', 'Display', 'off' );
counter_inf = 0 ;
%Integrazione utilizzando Eulero-Esplicito
for n=1:length(tspan)-1
    x_curr = yout(:,n);

    f_val = f_fun(x_curr);
    g_val = g_fun(x_curr);

    grad_h1 = grad_h1_fun(x_curr);
    val_h1  = h1_fun(x_curr);
    Lf_h1   = dot(grad_h1, f_val);
    Lg_h1   = dot(grad_h1, g_val);

    grad_b1 = grad_b1_fun(x_curr);
    val_b1  = b1_fun(x_curr);
    Lf_b1   = dot(grad_b1, f_val);
    Lg_b1   = dot(grad_b1, g_val);

    grad_v = grad_V_fun(x_curr);
    v = V_fun(x_curr);
    Lf_v = dot(grad_v,f_val);
    Lg_v = dot(grad_v,g_val);
    la = 3;

    k2_adapt = k2;
    a2_adapt = a2;
    if val_h1 < 0.2
        k2_adapt = k2 * 0.6;
    end
    if val_b1 < 0.2
        a2_adapt = a2 * 0.6;
    end



    y_des = [x_tgt;0;0;0];
    u_nom = - K * (x_curr - y_des) ;



    H = diag([1,1e6 , 1e10 ]);
    f_qp = [-u_nom; 0; 0];
    A = [ -Lg_h1,  -1 , 0 ;
        -Lg_b1, 0 , -1  ;
        %Lg_v , 0 , 0 ,-1
        ];
    b = [ Lf_h1 + k2_adapt * val_h1;
        Lf_b1 + a2_adapt * val_b1 ];
    %-Lf_v - la*v];

    lb = [u_min; 0;0];
    ub = [u_max; inf ; inf;];


    [z, ~, flag] = quadprog(H, f_qp, A, b, [], [], lb, ub, [], options);

    if flag == -2
        warning("Problem Infeasible")
        counter_inf = counter_inf + 1;
        % Fallback: usa LQR saturato
        u_emergency = max(u_min, min(u_max, u_nom));
        u_applied(:,n) = u_emergency;
    else
        u_applied(:,n) = z(1);
        % Avvisa se le slack sono attive
        if z(2) > 1e-3 || z(3) > 1e-3
            %fprintf('Slack attive: [%.4f, %.4f] @ t=%.3f\n', z(2), z(3), tspan(n));
        end
    end
    yout(:,n+1) = yout(:,n)+ dt * (f_val + g_val*u_applied(:,n)) ;

end
tout = tspan ;
yout = yout.';
end



function [tout, yout, u_applied, counter_inf, phase_switch_idx] = cbf_then_stabilize(y, f_fun, g_fun, b1_fun, grad_b1_fun, h1_fun, grad_h1_fun, k2, a2, u_min, u_max, K, x_tgt, tstar, V_fun, grad_V_fun)

dt = 0.001;
max_iter = ceil(tstar / dt);
yout = zeros(4, max_iter);
u_applied = zeros(1, max_iter);
tout = zeros(1, max_iter);

yout(:,1) = y;
tout(1) = 0;
n = 1;

options = optimoptions('quadprog', 'Display', 'off');
counter_inf = 0;

% Soglia per cambio fase (SOLO POSIZIONE)
position_tolerance = 0.05;  % Quando |x - x_target| < 0.05, passa alla fase 2

phase = 1;  % 1 = navigazione, 2 = stabilizzazione
phase_switch_idx = [];

% LOOP PRINCIPALE
while n < max_iter && tout(n) < tstar
    x_curr = yout(:,n);
    
    % Controlla SOLO la posizione per il cambio fase
    position_error = abs(x_curr(1) - x_tgt);
    
    if phase == 1 && position_error < position_tolerance
        phase = 2;
        phase_switch_idx = n;
        fprintf('✓ Posizione target raggiunta a t=%.3f s\n', tout(n));
        fprintf('  Inizio stabilizzazione completa dello stato...\n');
    end
    
    f_val = f_fun(x_curr);
    g_val = g_fun(x_curr);
    
    if phase == 1
        % FASE 1: NAVIGAZIONE CON CBF
        % Obiettivo: raggiungere x_target evitando ostacoli
        
        grad_h1 = grad_h1_fun(x_curr);
        val_h1 = h1_fun(x_curr);
        Lf_h1 = dot(grad_h1, f_val);
        Lg_h1 = dot(grad_h1, g_val);
        
        grad_b1 = grad_b1_fun(x_curr);
        val_b1 = b1_fun(x_curr);
        Lf_b1 = dot(grad_b1, f_val);
        Lg_b1 = dot(grad_b1, g_val);
        
        % Tuning adattivo
        k2_adapt = k2;
        a2_adapt = a2;
        if val_h1 < 0.2
            k2_adapt = k2 * 0.6;
        end
        if val_b1 < 0.2
            a2_adapt = a2 * 0.6;
        end
        
        % Controllo nominale LQR verso il target
        y_des = [x_tgt; 0; 0; 0];
        u_nom = -K * (x_curr - y_des);
        
        % QP con vincoli CBF
        H = diag([1, 1e8, 1e10]);
        f_qp = [-u_nom; 0; 0];
        A = [-Lg_h1, -1,  0;
             -Lg_b1,  0, -1];
        b_qp = [Lf_h1 + k2_adapt * val_h1;
                Lf_b1 + a2_adapt * val_b1];
        lb = [u_min; 0; 0];
        ub = [u_max; inf; inf];
        
        [z, ~, flag] = quadprog(H, f_qp, A, b_qp, [], [], lb, ub, [], options);
        
        if flag == -2
            counter_inf = counter_inf + 1;
            u_emergency = max(u_min, min(u_max, u_nom));
            u_applied(n) = u_emergency;
            if counter_inf == 1
                warning('Problema infeasible @ t=%.3f, usando fallback LQR', tout(n));
            end
        else
            u_applied(n) = z(1);
            if z(2) > 1e-3 || z(3) > 1e-4
                % fprintf('Slack: [%.4f, %.4f] @ t=%.3f\n', z(2), z(3), tout(n));
            end
        end
        
    else
        % FASE 2: STABILIZZAZIONE (SOLO LQR)
        % Obiettivo: stabilizzare tutto lo stato [x, θ, v_x, v_θ] → [x_tgt, 0, 0, 0]
        
        y_des = [x_tgt; 0; 0; 0];
        error_state = x_curr - y_des;
        u_nom = -K * error_state;
        
        % Semplice saturazione (niente CBF)
        u_applied(n) = max(u_min, min(u_max, u_nom));
        
        % Criterio di stop: sistema completamente stabilizzato
        state_norm = norm(x_curr - x_tgt,inf);
        if state_norm < 1e-3
            fprintf('✓ Sistema completamente stabilizzato a t=%.3f s\n', tout(n));
            fprintf('  Stato finale: x=%.4f, θ=%.4f, v_x=%.4f, v_θ=%.4f\n', ...
                    x_curr(1), x_curr(2), x_curr(3), x_curr(4));
            n = n + 1;  % Includi l'ultimo punto
            break;
        end
    end
    
    % Integrazione Eulero esplicito
    yout(:,n+1) = yout(:,n) + dt * (f_val + g_val * u_nom);
    tout(n+1) = tout(n) + dt;
    n = n + 1;
end

% Taglia gli array alla lunghezza effettiva
tout = tout(1:n);
yout = yout(:, 1:n)';
u_applied = u_applied(1:n)';

%% STATISTICHE FINALI
fprintf('\n═══════════════════════════════════════\n');
fprintf('STATISTICHE SIMULAZIONE\n');
fprintf('═══════════════════════════════════════\n');
fprintf('Durata totale:     %.3f s\n', tout(end));
fprintf('Iterazioni totali: %d\n', n);

if ~isempty(phase_switch_idx)
    fprintf('\nFASE 1 (Navigazione CBF):\n');
    fprintf('  Durata: %.3f s\n', tout(phase_switch_idx));
    fprintf('  Problemi infeasible: %d\n', counter_inf);
    
    fprintf('\nFASE 2 (Stabilizzazione):\n');
    fprintf('  Durata: %.3f s\n', tout(end) - tout(phase_switch_idx));
    fprintf('  Stato finale:\n');
    fprintf('    x     = %.4f m (target: %.4f)\n', yout(end,1), x_tgt);
    fprintf('    θ     = %.4f rad\n', yout(end,2));
    fprintf('    v_x   = %.4f m/s\n', yout(end,3));
    fprintf('    v_θ   = %.4f rad/s\n', yout(end,4));
else
    fprintf('\n⚠ Target non raggiunto entro il tempo simulato\n');
    fprintf('  Distanza finale: %.4f m\n', abs(yout(end,1) - x_tgt));
end
fprintf('═══════════════════════════════════════\n\n');

end