function plot_safe_set(valori_numerici , obs , tgt , y0 )

l = valori_numerici(3);
d = obs.d;
eps = obs.epsilon;

figure;
hold on;

% Definiamo la funzione h0
mia = @(x, theta) x.^2 + l^2 + d^2 - 2*l*x.*sin(theta) - 2*d*l*cos(theta) - eps^2;

fimplicit(mia, [y0(1)-1, tgt+1, -pi/2, pi/2], 'LineWidth', 2, 'Color', 'r');

grid on;
axis tight;
xlabel('Posizione Carrello x [m]');
ylabel('Angolo Pendolo \theta [rad]');

scatter(y0(1), y0(2), 60, 'blue', 'filled', 'DisplayName', 'Start');
scatter(tgt, 0, 60, 'green', 'filled', 'DisplayName', 'Target');

legend('Ostacolo', 'Start', 'Target');

% Imposto la "leggenda" dell'asse y in radianti
set(gca, 'YTick', -pi/2:pi/4:pi/2); 
set(gca, 'YTickLabel', {'-\pi/2','-\pi/4'  ,'0','\pi/4', '\pi/2'});
set(gca, 'TickLabelInterpreter', 'tex');

hold off;
end