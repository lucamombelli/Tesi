function plot_pendolo(tout , yout ,u_out ,u_max,u_min,x_target_pos ,theta_critico , valori_numerici ,obs)

j_obs = find(yout(:,1)>= 0 , 1) ; % inidice per cui trovo l'ostacolo
j_target = find(yout(:,1)>= x_target_pos ,1) ; % indice per cui ho superato x_target_pos

figure();
subplot(3,1,1);
plot(tout(), yout(:,1), 'LineWidth', 2);
hold on;
yline(x_target_pos, 'g--', 'Target');
xline(tout(j_obs) , 'r--' , 'Obstacle Center')
grid on;
ylabel('x[m]'); title(['Posizione (Target: ' num2str(x_target_pos) 'm)']);
hold off ;

subplot(3,1,2);
plot(tout, yout(:,2), 'r', 'LineWidth', 2);
hold on;
yline(double(theta_critico), 'k--');
yline(-double(theta_critico), 'k--');
xline(tout(j_obs) , 'r--' , 'Obstacle Center')
grid on; ylabel('Theta [rad]'); title('Angolo Pendolo');

subplot(3,1,3);
plot(tout, u_out, 'b', 'LineWidth', 1.5);
yline(u_max, 'r--'); yline(u_min, 'r--');
xline(tout(j_obs) , 'r--' , 'Obstacle Center')
grid on; ylabel('Controllo u [N]'); xlabel('Tempo [s]');
title('Controllore Ottimo')

if(~j_target)
    figure();
    subplot(2,1,1)
    plot(tout , yout(:,3) , 'LineWidth',1.5) ;
    xline(tout(j_target), 'r--' , 'Obstacle Center')
    title('Velocita x')
    grid on;
    subplot(2,1,2)
    plot(tout , yout(:,4) , 'LineWidth',1.5) ;
    xline(tout(j_target), 'r--' , 'Obstacle Center')
    title('Velocita theta')
    grid on ;
end
% PLOT SPAZIALE
figure;

% Posizione punta pendolo nello spazio
x_punta = yout(:,1) - valori_numerici(3)*sin(yout(:,2));
y_punta = valori_numerici(3)*cos(yout(:,2));
massimo = max(x_punta) ;
% Plot traiettoria
plot(x_punta, y_punta, 'b-', 'LineWidth', 2);
hold on;

% Ostacolo
theta_circle = linspace(0, 2*pi, 100);
x_circle = obs.epsilon * cos(theta_circle);
y_circle = obs.epsilon * sin(theta_circle) + obs.d;
plot(x_circle, y_circle, 'r-', 'LineWidth', 2);
fill(x_circle, y_circle, 'r', 'FaceAlpha', 0.7);

% Traiettoria del carrello
plot(yout(:,1), zeros(size(yout(:,1))), 'LineWidth', 1.5 , 'Color','black');

% Punti critici
xline(x_target_pos,'LineWidth', 2 ,'LineStyle','--' ,'Color','g','Label','Target')

axis equal;
grid on;
xlabel('Posizione x [m]');
ylabel('Altezza y [m]');
title('Vista Spaziale: Traiettoria Punta Pendolo');
legend('Punta pendolo', 'Ostacolo');
if massimo <= 10
    xlim([-2, massimo]);
else
    xlim([-2, 10]);
end
ylim([-0.5, valori_numerici(3)+0.5]);
hold off;

x_punta = yout(:,1) - valori_numerici(3)*sin(yout(:,2));
y_punta = valori_numerici(3)*cos(yout(:,2));

% Distanza dal centro dell'ostacolo (0, d)
dist_from_obstacle = sqrt((x_punta).^2 + (y_punta - obs.d).^2);
clearance = dist_from_obstacle - obs.epsilon;
% Trova violazioni
violations = find(clearance < 0);
if ~isempty(violations)
    fprintf('VIOLAZIONE OSTACOLO rilevata!\n');
    fprintf('   Prima violazione a t=%.3f s\n', tout(violations(1)));
    fprintf('   Penetrazione massima: %.4f m\n', min(clearance));
else
    fprintf('✓ Nessuna violazione ostacolo\n');
    fprintf('  Margine minimo: %.4f m\n', min(clearance));
end