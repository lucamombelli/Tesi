using CairoMakie
using  LaTeXStrings

set_theme!(theme_latexfonts())

x_val = 0:0.01:5;
y_val = 0:0.01:5;

x_obs = 3 
y_obs = 2.8

r_obs = 0.5
h(x,y) = sqrt((x-x_obs)^2 + (y-y_obs)^2) - r_obs 

Z = [h(x, y) for x in x_val, y in y_val]

fig = Figure(size = (600, 400))
ax = Axis(fig[1, 1], 
    title = L"Regione\ Ammissibile \ S", 
    xlabel = L"x \text{ [m]}", 
    ylabel = L"y \text{ [m]}",
    aspect = DataAspect() 
)


# f(x,theta) >= 0 
contourf!(ax, x_val, y_val, Z, 
    levels = [0.0, maximum(Z) + 0.5], # Il range da colorare
   colormap = (:greens , 0.4 ) ,      # Scegli il colore della regione
)

# 6. f(x,θ) = 0
contour!(ax, x_val, y_val, Z, 
    levels = [0.0], 
    color = :red, 
    linewidth = 2
)

fig 
save("/Users/luca/Tesi/immagini/example_cbf.pdf", fig)
println("Grafico salvato con successo come example_cbf.pdf!")

