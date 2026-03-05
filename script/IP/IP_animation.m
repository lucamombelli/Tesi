function IP_animation(x1, x2, obs, parameters, xG , X0)
    
    % IP on cart parameters
    W  = 0; % width of cart
    H  = 0; % hight of Cart
    L  = parameters(3); % length of pendulum  
    wr = 0; % right wheel 
    
    % Position coordination
    y = H/2 + wr/2;
    w1x = x1 -0.9*W/2;
    w1y = 0;
    w2x = x1 + 0.9*W/2 - wr;
    w2y = 0;
    
    % Position of pendulum 
    px = x1 - L*sin(x2);
    py = y + L*cos(x2);
    
    % Position of the obstacle
    mesh = linspace(0,2*pi,100);
    x_obs = obs.radius*cos(mesh) ; 
    y_obs = obs.radius *sin(mesh) + obs.height ;
    
    % Plot
    base = plot([-15 15],[0 0],'k','LineWidth',2); % base line
    hold on;
    pendulum = plot([x1 px],[y py],'b','LineWidth',2); % Pendulum rod 
    p_cir = viscircles([px py],0.01,'Color',[1 0.1 0.1],'LineWidth',2.5); % Pendulum Trajectory(Circle)
    obsacle = fill(x_obs,y_obs,'r') ;
    target = xline(xG ,'Color','g','LineWidth',0.5,'LineStyle','--','Label','Target Positon');
    star = xline(X0(1) ,'Color','b','LineWidth',0.5,'LineStyle','--','Label','Start Positon');
    xlabel('X (m)');
    ylabel('Y (m)');
    title('Inverted pendulum')
    axis(gca,'equal');
    xlim([X0(1)-1 xG+1]);
    ylim([-0.5 1.5 ]);
    grid on;   
    hold off;
end