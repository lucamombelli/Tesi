function phase_portrait(xout)
    
    % State variables for phase portrait
    theta = xout(:,2);  % angle
    dtheta = xout(:,4);  % angular velocity
    n = numel(theta);  % number of samples
    
    % Phase portrait (target equilibrium)
    figure;
    xline(0,'--','Color','white','LineWidth',1.2)
    hold on
    yline(0,'--','Color','white','LineWidth',1.2)
    hold on
    plot(theta(round(n/4):end),dtheta(round(n/4):end),'b','LineWidth',1.5)
    title('Phase portrait at target','Interpreter','latex','FontSize',12)
    xlabel('$\theta$ (rad)','Interpreter','latex','FontSize',12);
    ylabel('$\dot{\theta}$ (rad/s)','Interpreter','latex','FontSize',12);
    xlim([-6e-3,6e-3])
    ylim([-0.020, 0.020])
    grid on
end