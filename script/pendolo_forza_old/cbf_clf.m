function [tout , yout , u_applied ] = cbf_clf(y0, f_fun, g_fun, b1_fun, grad_b1_fun,h1_fun , grad_h1_fun , k2 ,a2, u_min, u_max, K, x_tgt,tstar,V_fun,grad_V_fun)

dt = 0.001;
tspan = 0:dt:tstar ;
yout = zeros(4,size(tspan,2)) ;
u_applied  = zeros(1,size(tspan,2))  ;
n = 1 ;
yout(:,n) = y0 ;
options = optimoptions('quadprog', 'Display', 'off' );
count_inf = 0 ;
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
    la = 0.3;

    alpha2 = @(h) k2 * h + 0.5 * h^3;

    %y_des = [x_tgt;0;0;0];
    %u_nom = - K * (x_curr - y_des) ;

    H = diag([eps, 0 , 0 , eps]);
    f_qp = [ 0  ; 0 ; 0 ; 0 ];
    A = [ -Lg_h1 , 0  , 0  ,0 ;
        %-Lg_b1, 0 , -1 ,0   ;
     Lg_v , 0 , 0 ,1 ];
    b = [Lf_h1 + 0.5*alpha2(val_h1);
        %Lf_b1 + a2 * val_b1;
    -Lf_v - la*v];

    lb = [u_min; 0;0;0 ];
    ub = [u_max; inf;inf; inf ];


    [z, ~, flag] = quadprog(H, f_qp, A, b, [], [], lb, ub, [], options);

    if flag == -2
        fprintf('Problem Infeasible at %f \n' , tspan(n));
        count_inf = count_inf + 1;
        % Fallback: usa LQR saturato

        K_brake = K;
        target_here = [x_curr(1); 0; 0; 0]; % Target = Dove sono ora
        u_brake = -K_brake * (x_curr - target_here);


        u_emergency = max(u_min, min(u_max, u_brake));
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

if(counter_inf ~= 0)
    fprintf('Problema Infeasible in %i istanti ' , counter_inf);
end


end