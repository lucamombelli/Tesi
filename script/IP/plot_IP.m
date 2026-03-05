function plot_IP(tout, xout, u_out, u_max, u_min, xG, theta_crit, numerical_param, obs)

    j_obs = find(xout(:,1)>=0,1) ; % index where the obstacle can be found
    j_target = find(xout(:,1)>=xG,1) ; % index for which I passed the target
    
    Fig1=figure();
    subplot(3,1,1);
    plot(tout(),xout(:,1),'LineWidth',2);
    hold on;
    yline(xG,'g--','Target');
    xline(tout(j_obs),'r--','Obstacle Center')
    grid on;
    ylabel('$x$ [m]','Interpreter','latex','FontSize',12); 
    title(['Position (Target: ' num2str(xG) 'm)']);
    hold off ;
    
    subplot(3,1,2);
    plot(tout,xout(:,2),'r','LineWidth',2);
    hold on;
    yline(double(theta_crit),'k--');
    yline(-double(theta_crit),'k--');
    xline(tout(j_obs),'r--','Obstacle Center')
    grid on; 
    ylabel('$\theta$ [rad]','Interpreter','latex','FontSize',12); 
    title('IP angle');
    
    subplot(3,1,3);
    plot(tout,u_out,'b','LineWidth',1.5);
    %yline(u_max,'r--'); yline(u_min,'r--');
    xline(tout(j_obs),'r--','Obstacle Center')
    grid on; ylabel('u [N]'); xlabel('Time [s]');
    title('Optimal control')
    
    if(~j_target)
        figure();
        subplot(2,1,1)
        plot(tout,xout(:,3),'LineWidth',1.5) ;
        xline(tout(j_target),'r--','Obstacle Center')
        title('x velocity')
        grid on;
        subplot(2,1,2)
        plot(tout,xout(:,4),'LineWidth',1.5) ;
        xline(tout(j_target),'r--','Obstacle Center')
        title('theta velocity')
        grid on ;
    end

    % Plot
    Fig2 = figure;
    
    % Pendulum tip position
    x_tip = xout(:,1) - numerical_param(3)*sin(xout(:,2));
    y_tip = numerical_param(3)*cos(xout(:,2));
    maximum = max(x_tip);

    % Trajectory
    plot(x_tip,y_tip,'b-','LineWidth',2);
    hold on;
    
    % Obstacle
    theta_circle = linspace(0,2*pi,100);
    x_circle = obs.radius * cos(theta_circle);
    y_circle = obs.radius * sin(theta_circle) + obs.height;
    plot(x_circle,y_circle,'r-','LineWidth',2);
    fill(x_circle,y_circle,'r','FaceAlpha',0.7);
    
    % Trajectory
    plot(xout(:,1),zeros(size(xout(:,1))),'LineWidth',1.5,'Color','black');
    
    % Critical points
    xline(xG,'LineWidth',2,'LineStyle','--','Color','g','Label','Target')
    
    axis equal;
    grid on;
    xlabel('Position [m]','Interpreter','latex','FontSize',12);
    ylabel('Height [m]','Interpreter','latex','FontSize',12);
    title('Inverse Pendulum tip trajectory');
    legend('Pendulum tip', 'Obstacle','Location','southeast','Interpreter','latex','FontSize',12);
        
    %Controllo della distanza dall'ostacolo
    if maximum <= 10
        xlim([-4, maximum]);
    else
        xlim([-4, 10]);
    end
    ylim([-0.5, numerical_param(3)+0.5]);
    hold off;
    
    x_tip = xout(:,1) - numerical_param(3)*sin(xout(:,2));
    y_tip = numerical_param(3)*cos(xout(:,2));
    
    % Distanza dal centro dell'ostacolo (0, d)
    dist_from_obstacle = sqrt((x_tip).^2 + (y_tip - obs.height).^2);
    clearance = dist_from_obstacle - obs.radius;
    % Trova violazioni
    violations = find(clearance < 0);
    if ~isempty(violations)
        fprintf('CONSTRAINT VIOLATION DETECTED!\n');
        fprintf('   First violation at t=%.3f s\n', tout(violations(1)));
        fprintf('   Maximal penetration: %.4f m\n', min(clearance));
    else
        fprintf('✓ No obstacle violation\n');
        fprintf('  Minimum margin: %.4f m\n', min(clearance));
    end

%saveas(Fig1,'Fig1_eu_3','png');
%saveas(Fig2 , 'Fig2_eu_3','png');


end

