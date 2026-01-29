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
        
        % Controllo nominale LQR verso il target
        y_des = [x_tgt; 0; 0; 0];
        u_nom = -K * (x_curr - y_des);

        % QP con vincoli CBF
        H = diag([1, 1e10, 1e10]);
        f_qp = [-u_nom; 0; 0];
        A = [-Lg_h1, 0,  0;
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

% STATISTICHE FINALI
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

end