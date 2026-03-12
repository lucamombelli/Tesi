clear; clc;
% 1. Definisci le variabili simboliche (stato)
syms x theta vx vtheta real
state = [x;theta;vx;vtheta];
% Parametri fisici 
syms M m l g real
sym_par = [M, m, l, g ]';
sys_par = [1, 0.1, 1, 9.81 ]';

% Parametri fisici 
f = [ vx ; 
      vtheta;
      (-m*l*sin(theta)*vtheta^2+g*m*cos(theta)*sin(theta)) / (M+m*sin(theta)^2)
      (-m*l*cos(theta)*sin(theta)*vtheta^2 + g*(M+m)*sin(theta) ) / (l*(M+m*sin(theta)^2)) ; 
        ] ;

g = [ 0 ; 0 ; 1/(M+m*sin(theta)^2) ; cos(theta)/(l*(M+m*sin(theta)^2))] ; 

lie_brackets = @(f,g) simplify(jacobian(g,state)*f - jacobian(f,state)*g) ;

b = lie_brackets(f,g);
b1 = lie_brackets(f,b);
b2 = lie_brackets(f,b1);
b3 = lie_brackets(f,b2);
b4=lie_brackets(g,f);
b5=lie_brackets(b4,g);

D = [f,g,b,b1,b2 ,b3,b4,b5];
D_num = subs(D, [state ; sym_par] , [ zeros(4,1) ; sys_par]);
rank(D_num) 

%We want to obtain a system in the form dx/dt = Ax+Bu
A_sym = jacobian(f, state);
A = subs(A_sym, [state ; sym_par] , [ zeros(4,1) ; sys_par]);
B_sym = g ; 
B = subs(B_sym, [state ; sym_par] , [ zeros(4,1) ; sys_par]);
% Controllability iff rank(crt(A,B) ) == 4 
if (rank(ctrb(A,B)) == 4 ) 
    fprintf('The system is small time local controllable \n');
else
    fprintf('Porco Dio');
end