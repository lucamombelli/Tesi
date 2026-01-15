function animation =  IP_Animation(x,th,d ,epsilon, valori_numerici,tgt)

persistent j
if isempty(j)
    j = 0;
end

j = j+1;

W  = 0; % width of cart
H  = 0; % hight of Cart
L  = valori_numerici(3); % length of pendulum  
wr = 0; % right wheel 

% Position coordination
y = H/2+wr/2;
w1x = x -0.9*W/2;
w1y = 0;
w2x = x+0.9*W/2-wr;
w2y = 0;

% position of pendulum 
px = x - L*sin(th);
py = y + L*cos(th);

%Position of the obstacle
mesh = linspace(0,2*pi,100);
x_obs = epsilon*cos(mesh) ; 
y_obs = epsilon *sin(mesh) + d ;

base = plot([-15 15],[0 0],'k','LineWidth',2); % base line
hold on;
%cart = rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',0.1,'FaceColor',[1 0.1 0.1],'EdgeColor',[1 1 1]);
%left_wheel  = rectangle('Position',[w1x,w1y,wr,wr],'Curvature',1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
%right_wheel = rectangle('Position',[w2x,w2y,wr,wr],'Curvature',1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
    
pendulum = plot([x px],[y py],'b','LineWidth',2.5); % Pendulum rod 
p_cir = viscircles([px py],0.02,'Color',[1 0.1 0.1],'LineWidth',2.5); % Pendulum Trajectory(Circle)
%p_cir1 = viscircles([x y],0.02,'Color','w','LineWidth',0.2); % center of Cart
obstacle = fill(x_obs,y_obs,'r') ;
target = scatter(tgt , 0 , 50 , 'g','filled' );


xlabel('X (m)');
ylabel('Y (m)');
title('Inverted pendulum')
axis(gca,'equal');
xlim([-5 5]);
ylim([-1 1.5 ]);
grid on;   
   
   