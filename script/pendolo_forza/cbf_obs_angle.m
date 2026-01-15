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