using CairoMakie
using  LaTeXStrings
# 1. Definisci la griglia spaziale (i limiti del tuo grafico)
x_vals = -2.0:0.05:2.0
theta_vals = -(2*pi)/3:0.05: (2*pi)/3
v_theta_vals = - pi:0.05:pi 

R = 0.1
l = 1
M = 1 
m = 0.1
dy = 1 
dx = 0 
alpha1 = 40 ; 
theta_critico = pi/2 ; 

set_theme!(theme_latexfonts())

f(x, y) =  x^2 + l^2 + dx^2 +dy^2 -2( x*l*sin(y)+x*dx-dx*l*sin(y)+dy*l*cos(y)) - R^2 

Z = [f(x, y) for x in x_vals, y in theta_vals]

fig = Figure(size = (600, 400))
ax = Axis(fig[1, 1], 
    title = L"Regione Ammissibile S_1", 
    xlabel = L"x \text{ [m]}", 
    ylabel = L"\theta \text{ [rad]}",
    yticks = (
        [-pi, -pi/2, 0, pi/2, pi],                      # Posizioni numeriche
        [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]    # Etichette in LaTeX
    )
)


# f(x,theta) >= 0 
contourf!(ax, x_vals, theta_vals, Z, 
    levels = [0.0, maximum(Z) + 0.5], # Il range da colorare
   colormap = (:greens , 0.4 ) ,      # Scegli il colore della regione
)

# 6. f(x,θ) = 0
contour!(ax, x_vals, theta_vals, Z, 
    levels = [0.0], 
    color = :red, 
    linewidth = 2
)

fig 
save("/Users/luca/Tesi/immagini/obstacle_cbf.pdf", fig)
println("Grafico salvato con successo come obstacle_cbf.pdf!")
# ---------------------------------------------------------------------------------------------------
b0(x , y) = theta_critico^2 - y^2
B0 = [b0(y) for x in x_vals ,y in theta_vals] 

fig1= Figure(size = (600, 400))
ax = Axis(fig1[1, 1], 
    title = "Regione Ammissibile", 
    xlabel = L"x \text{ [m]}", 
    ylabel = L"\theta \text{ [rad/s]}",
        
     yticks = (
        [ -(2*pi)/3, -pi/2, 0, pi/2,  +(2*pi)/3],                      # Posizioni numeriche
        [L"-\frac{2 \pi}{3}", L"-\frac{\pi}{2}", L"0", L"\frac{\pi}{2}", L"\frac{2\pi}{3}"]    # Etichette in LaTeX
    )
)

contourf!(ax, x_vals,theta_vals ,B0, 
    levels = [0.0, maximum(B) + 0.5], 
    colormap = (:algae , 0.4 ) ,        
)

contour!(ax, x_vals,theta_vals ,B0, 
    levels = [0.0], 
    color = :red, 
    linewidth = 2
)

fig1 
save("/Users/luca/Tesi/immagini/angle_cbf.pdf", fig1)
println("Grafico salvato con successo come angle_velocities_pendolo.pdf!")

# ---------------------------------------------------------------------------------------------------
b(y , dy) = -2*dy*y + alpha1*(theta_critico^2 - y^2) 
B = [b(y, dy) for x in theta_vals, y in v_theta_vals] 

fig1= Figure(size = (600, 400))
ax = Axis(fig1[1, 1], 
    title = "Regione Ammissibile", 
    xlabel = L"\theta \text{ [rad]}", 
    ylabel = L"v_\theta \text{ [rad/s]}",
    
    xticks = (
        [-pi, -pi/2, 0, pi/2, pi],                      # Posizioni numeriche
        [L"-\pi", L"-\frac{\pi}{2}", L"0", L"\frac{\pi}{2}", L"\pi"]    # Etichette in LaTeX
    ),
    
     yticks = (
        [-pi, -pi/2, 0, pi/2, pi],                      # Posizioni numeriche
        [L"-\pi", L"-\frac{\pi}{2}", L"0", L"\frac{\pi}{2}", L"\pi"]    # Etichette in LaTeX
    )
)

contourf!(ax, theta_vals, v_theta_vals, B, 
    levels = [0.0, maximum(B) + 0.5], 
    colormap = (:algae , 0.4 ) ,        
)

contour!(ax, theta_vals, v_theta_vals, B, 
    levels = [0.0], 
    color = :red, 
    linewidth = 2
)

fig1 
save("/Users/luca/Tesi/immagini/angle_velocities_cbf.pdf", fig1)
println("Grafico salvato con successo come angle_pendolo.pdf!")

# ---------------------------------------------------------------------------------------------------