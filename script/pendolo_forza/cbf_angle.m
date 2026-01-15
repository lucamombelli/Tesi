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
    u_nom = 1.2 ;

    y_des = [x_tgt;0;0;0];
    %u_nom = - K * (x_curr - y_des) ;


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
    %eulero Implicito
    yout(:,n+1) = yout(:,n)+ dt * (f_val + g_val*u_applied(:,n)) ;
end
tout = tspan ;
yout = yout.';
end
