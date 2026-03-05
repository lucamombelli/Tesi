function [tout,xout,u_applied] = IP_simulation(X, f_fun, g_fun, b1_fun, grad_b1_fun, ...
                                               h1_fun , grad_h1_fun , k2, a2, u_min, u_max, K, xg, Tf)
    
    % Parameters
    dt = 0.001;  % time-step
    tVec = 0:dt:Tf;  % time vector
    xout = zeros(4,numel(tVec));  % memory preallocation
    u_applied  = zeros(1,numel(tVec));  % memory preallocation
    xout(:,1) = X ;
    options = optimoptions('quadprog','Display','off');
    counter_inf = 0 ;
    counter_no_direction = 0;
    
    % Control loop
    for n = 1:length(tVec)-1
        XC = xout(:,n);  % current state
        
        % System components
        f_val = f_fun(XC);  % drift vector field
        g_val = g_fun(XC);  % input map
    
        % SO-CBFs (obstacle)
        val_h1 = h1_fun(XC);  % cbf
        grad_h1 = grad_h1_fun(XC);  % gradient
        Lf_h1 = dot(grad_h1,f_val);  % directional derivative along f(x)
        Lg_h1 = dot(grad_h1,g_val);  % directional derivative along g(x)
        
        % SO-CBFs (critical angle)
        val_b1 = b1_fun(XC);  % cbf
        grad_b1 = grad_b1_fun(XC);  % gradient
        Lf_b1 = dot(grad_b1,f_val); % directional derivative along f(x)
        Lg_b1 = dot(grad_b1,g_val); % directional derivative along g(x)
        
        % Target state and control
        XG = [xg, 0, 0, 0]';  
        u_lqr = -K*(XC-XG) ;
        
        % Constrained optimization with "quadprog"
        % weight if we don't turn a hard costrain is 1e12
        H = diag([1e-1, 1e12, 1e14]); % weight matrix
        f_qp = [-u_lqr, 0, 0]';  % linear term

        % Inequality constraint of the form Au <= b
        A = [-Lg_h1, -1, 0; ...
             -Lg_b1, 0, -1];  % left-hand-side
        b = [Lf_h1 + k2*val_h1, Lf_b1 + a2*val_b1]';  % rght-hand-side
        
       
        % Bounded control 
        lb = [u_min, 0, 0]';  % lower bound
        ub = [u_max, inf, inf]';  % upper bound
        
        % Quadprog
        [z,~,flag] = quadprog(H, f_qp, A, b, [], [], lb, ub, [], options);
        
        % Optimization flags
        if (flag == -2)
            error('Problem Infeasible at %f \n' , tVec(n));
        else
            u_applied(:,n) = z(1);
            % Avvisa se le slack sono attive
            if z(2) > 1e-3 || z(3) > 1e-3
                %fprintf('Slack attive: [%.4f, %.4f] @ t=%.3f\n', z(2), z(3), tspan(n));
            end
        end
        
        % Explicit Euler (metti RK se riesci) 
        xout(:,n+1) = xout(:,n) + dt*(f_val+g_val*u_applied(:,n)) ;
     
        if (flag == -8) 
            counter_no_direction = counter_no_direction + 1 ; 
        end
    
    end
    if (counter_inf ~= 0)
        fprintf('Problem is infeasible  %i times \n ' , counter_inf),
    else
         fprintf("No problem! (Feasible) \n");
    end
    
    if (counter_no_direction ~= 0)
        fprintf('Unable to compute a step direction  %i times \n ' , counter_no_direction),
    else
        fprintf("No problem! (Computed) \n");
    end
    tout = tVec ;
    xout = xout.';
end