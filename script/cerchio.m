close all ;
clear all ; 

r = sqrt(5) ;
t = -r : 0.001:r ;
f1 = sqrt( 5 - t.^2 ) ;
f2 = - sqrt(5 - t.^2 ) ;

[ x , y ] = meshgrid( -3 : 3 ) ; 
t1 = -4 : 0.001:4


figure (1) 
hold on ; 
c_plot=plot (t, f1 ,'b' ,  t , f2 , 'b'   ) ; 
axis equal ; 
xlabel('x') 
ylabel('y') 
grid on ; 
q1=quiver ( x , y , x , y.^2 , 'r' )  ; 
q2=quiver( x , y  , 2.*x , 2.*y ) ;
q3 = plot( t1 , - (t1.^2).^(1/3)) 
legend ([c_plot(1) , q1 , q2 , q3 ] , "x^2+y^2=5 " ," (x,y^2) " , "(2x,2y)" , "y = -x^{2/3}"  )
print('../immagini/circonferenza.png ', '-dpng' , '-r300'  ); % Salva come PNG a 300 DPI
hold off ;

figure (2) 
hold on ; 
c_plot=plot (t, f1 ,'b' ,  t , f2 , 'b'   ) ; 
axis equal ; 
xlabel('x') 
ylabel('y') 
grid on ; 
q1=quiver ( x , y , -x , - y.^3, 'r' )  ; 
q2=quiver( x , y  , 2.*x , 2.*y ) ;
legend ([c_plot(1) , q1 , q2 , q3] , "x^2+y^2=5 " ," (-x,-y^3) " , "(2x,2y)"   )
print('../immagini/circonferenza_1.png ', '-dpng' , '-r300'  ); % Salva come PNG a 300 DPI
hold off ;



