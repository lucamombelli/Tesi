function plot_safe_set(parameters, obs, xG, x0)
    
    % Parameters
    l = parameters(3);
    dy = obs.height;
    eps = obs.radius;
    dx = 0 ;
    
    figure;
    hold on;
    
    % h0(x1,x2) definition
    h0 = @(x1, x2) x1.^2 + l^2 + dy^2 + dx^2 -2*x1*dx + 2*dx*l*sin(x2) - 2*l*x1.*sin(x2) - 2*dy*l*cos(x2) - eps^2;
    fimplicit(h0,[x0(1)-1,xG+1,-pi/2,pi/2],'LineWidth',2,'Color','r');
    
    grid on;
    axis tight;
    xlabel('$x$ [m]','Interpreter','latex','FontSize',12);
    ylabel('$\theta$ [rad]','Interpreter','latex','FontSize',12);
    
    scatter(x0(1),x0(2),60,'blue','filled','DisplayName','Start');
    scatter(xG,0,60,'green','filled','DisplayName','Target');
  
    legend('Obstacle','Start','Target','Location','southeast','Interpreter','latex','FontSize',12);
    
    % Imposto la "leggenda" dell'asse y in radianti
    set(gca,'YTick',-pi/2:pi/4:pi/2); 
    set(gca,'YTickLabel',{'-\pi/2','-\pi/4','0','\pi/4','\pi/2'});
    set(gca,'TickLabelInterpreter','tex');
    
    hold off;
end