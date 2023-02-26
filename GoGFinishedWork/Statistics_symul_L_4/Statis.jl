using Plots
using Distributions
using Random
using StatsPlots

Numar_aruncari = 1000;

# Generam aruncarile cap -> 1, pajura -> 0
N_esantion = 0:Numar_aruncari;
p_true = 0.5;
Random.seed!(12);
data = rand(Bernoulli(p_true), last(N_esantion));

# Introducem bias - metoda simpla/barbara de a introduce un bias care mai apoi sa se corecteze (se poate dezvolta)
#data = fill(1, 50)
#data = vcat(data, rand(Bernoulli(p_true), 50))

prior = Beta(1,1);

anim = Animation()
for (i, N) in enumerate(N_esantion)

    cap = sum(data[1:i-1])
    pajura = N - cap

    updated = Beta(prior.α + cap, prior.β + pajura)

    plot(updated,
    framestyle =:box,
    yaxis=([], false),
    grid = false,
    size = (1280,900),  
    title = "Ajustarea Bayesiană după $N aruncări",
    xlabel = "Probabilitatea de a nimeri cap", 
    ylabel = "", 
    legend = nothing,
    xlim = (0,1),
    xticks = [0.0, 0.2, 0.35, 0.5, 0.65, 0.8, 1.0],
    fill = 0, 
    color =:green,
    α = 0.5, w = 5);
    vline!([p_true], w = 3);

    frame(anim)
end;

gif(anim, "CoinTossAveraging.gif", fps = 10)

# Pentru viitor: comparatie pentru acelasi set de date -> inferenta Bayesiana vs o inferenta descriptiva simpla cu distributie Gauss,
# (μ = sum(data) & σ cu formula consacrata la fiecare iteratie) pentru comparare -> caz cinstit vs caz cu bias (outlieri)