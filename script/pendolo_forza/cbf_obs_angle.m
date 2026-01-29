function [tout , yout , u_applied] = cbf_obs_angle(y, f_fun, g_fun, b1_fun, grad_b1_fun,h1_fun , grad_h1_fun , k2 ,a2, u_min, u_max, K, x_tgt,tstar)

dt = 0.001;
tspan = 0:dt:tstar ;
yout = zeros(4,size(tspan,2)) ;
u_applied  = zeros(1,size(tspan,2))  ;
n = 1 ;
yout(:,n) = y ;
options = optimoptions('quadprog', 'Display', 'off' );
counter_inf = 0 ;
counter_no_direction = 0;

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

    y_des = [x_tgt;0;0;0];
    u_nom = - K * (x_curr - y_des) ;

    H = diag([1e-2,1e12,1e14]);
    f_qp = [ -u_nom ; 0 ; 0 ];
    A = [ -Lg_h1 , -1  , 0   ;
        -Lg_b1, 0 , -1   ];
    b = [Lf_h1 + k2 * val_h1;
        Lf_b1 + a2 * val_b1];

    lb = [u_min; 0;0  ];
    ub = [u_max; inf;inf];

    [z, ~, flag] = quadprog(H, f_qp, A, b, [], [], lb, ub, [], options);

    if flag == -2
        fprintf('Problem Infeasible at %f \n' , tspan(n));
        counter_inf = counter_inf + 1;
        A 
        b 
    %{
        K_brake = K;
        target_here = [x_curr(1); 0; 0; 0]; % Target = Dove sono ora
        u_brake = -K_brake * (x_curr - target_here);


        u_emergency = max(u_min, min(u_max, u_brake));
        u_applied(:,n) = u_emergency;
    %}
    else
        u_applied(:,n) = z(1);
        % Avvisa se le slack sono attive
        if z(2) > 1e-3 || z(3) > 1e-3
            %fprintf('Slack attive: [%.4f, %.4f] @ t=%.3f\n', z(2), z(3), tspan(n));
        end
    end
    yout(:,n+1) = yout(:,n)+ dt * (f_val + g_val*u_applied(:,n)) ;
    
    if flag == -8 
        counter_no_direction = counter_no_direction + 1 ; 
        %fprintf("Unable to compute a step direction \n")
        fprintf("%f\n at time %i",z(1),tspan(n));
    end

end
if(counter_inf ~= 0)
    fprintf('Problem is infeasible  %i times \n ' , counter_inf),
else
     fprintf("Senza problemi!( No infeasible) \n");
end

if(counter_no_direction ~= 0)
    fprintf('Unable to compute a step direction  %i times \n ' , counter_no_direction),
else
    fprintf("Senza problemi! ( No unable to compute ) \n");
end



tout = tspan ;
yout = yout.';
end