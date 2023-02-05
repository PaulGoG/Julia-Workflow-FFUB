using Plots
using LaTeXStrings
using LsqFit

# Cod de calcul pentru fitul parametrului Hubble

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru
const c = 2.998 * 1e8;

# Galaxii studiate: NGC 3627. NGC 3368, NGC 3147

function z(λ, λ₀, σ_λ)
    return [(λ - λ₀)/λ₀, sqrt(σ_λ^2 * (λ^2 + λ₀^2)/λ₀^4)]
end
function z_med(z, σ_z)
    Suma_Z = 0.0
    Suma_σ² = 0.0
    for i in eachindex(z)
        if z[i] != -1
            Suma_Z += z[i]
            Suma_σ² =+ σ_z[i]^2
        end
    end
    return [Suma_Z/length(z[z .!= -1]), sqrt(Suma_σ²)/length(z[z .!= -1])]
end
function distanta(d_et, s_et, s, σ_s_et, σ_s)
    return [d_et*s_et/s, sqrt((d_et*σ_s_et/s)^2 + (d_et*σ_s*s_et/s^2)^2)]
end
function fitare(v, d, σ_v, σ_d)
    model(t, p) = p[1].*t
    p0 = [1.0]
    wt = zeros(length(σ_v))
    for i in eachindex(σ_v)
        if σ_v[i] != 0 && σ_d[i] != 0
            wt[i] = 1/(σ_v[i]*σ_d[i])^2
        elseif σ_v[i] != 0
            wt[i] = 1/σ_v[i]^2
        elseif σ_d[i] != 0
            wt[i] = 1/σ_d[i]^2
        else
            wt[i] = 1.0
        end
    end
    fit = curve_fit(model, d, v, wt, p0)
    H_0 = last(fit.param)
    σ_H_0 = last(standard_errors(fit))
    return [H_0, σ_H_0]
end

# Ca-K, Ca-H, H-α, H-β, H-δ
λ₀ = [3933.7, 3968.5, 6562.8, 4861.3, 4101.7]; #Angstrom
σ_λ = 3.3/2;

d_1 = 10.73; #MPsec
σ_d_1 = 0;
s_1 = 5; #cm
σ_s_1 = 0.75;
λ_1 = [3941.67, 3980.0, 6578.34, 4873.33, 4113.3];
z_1 = [z(λ_1[i], λ₀[i], σ_λ)[1] for i in eachindex(λ₀)]; #adimensional
σ_z_1 = [z(λ_1[i], λ₀[i], σ_λ)[2] for i in eachindex(λ₀)];
z_mediu_1 = z_med(z_1, σ_z_1)[1];
σ_z_mediu_1 = z_med(z_1, σ_z_1)[2];

s_2 = 3.5; #cm
σ_s_2 = 0.8;
d_2 = distanta(d_1, s_1, s_2, σ_s_1, σ_s_2)[1]; #MPsec
σ_d_2 = distanta(d_1, s_1, s_2, σ_s_1, σ_s_2)[2];
λ_2 = [3944.995, 3980, 0.0, 0.0, 4116.33];
z_2 = [z(λ_2[i], λ₀[i], σ_λ)[1] for i in eachindex(λ₀)]; #adimensional
σ_z_2 = [z(λ_2[i], λ₀[i], σ_λ)[2] for i in eachindex(λ₀)];
z_mediu_2 = z_med(z_2, σ_z_2)[1];
σ_z_mediu_2 = z_med(z_2, σ_z_2)[2];

s_3 = 1.8; #cm
σ_s_3 = 0.5;
d_3 = distanta(d_1, s_1, s_3, σ_s_1, σ_s_3)[1]; #MPsec
σ_d_3 = distanta(d_1, s_1, s_3, σ_s_1, σ_s_3)[2];
λ_3 = [3970, 4006.66, 6620.0, 4916.66, 4136.66];
z_3 = [z(λ_3[i], λ₀[i], σ_λ)[1] for i in eachindex(λ₀)]; #adimensional
σ_z_3 = [z(λ_3[i], λ₀[i], σ_λ)[2] for i in eachindex(λ₀)];
z_mediu_3 = z_med(z_3, σ_z_3)[1];
σ_z_mediu_3 = z_med(z_3, σ_z_3)[2];

v = [0.0, z_mediu_1, z_mediu_2, z_mediu_3] .* c .* 1e-3; #km/s
d = [0.0, d_1, d_2, d_3]; 
σ_v = [0.0, σ_z_mediu_1, σ_z_mediu_2, σ_z_mediu_3] .* c .* 1e-3;
σ_d = [0.0, σ_d_1, σ_d_2, σ_d_3];

H_0 = fitare(v, d, σ_v, σ_d);

plt = scatter(
    d,
    v,
    xerr = σ_d, 
    yerr = σ_v,
    size = (1000, 950),
    dpi = 600,
    minorgrid = true,
    framestyle = :box,
    ylims = (0.0, last(v)+3*last(σ_v)),
    xlims = (0.0, last(d)+1.1*last(σ_d)),
    xlabel = "Distanța galaxiei față de observator [MPsc]",
    ylabel = "Viteza galaxiei față de observator [km/s]",
    title = latexstring("Determinarea parametrului Hubble \$\\mathrm{H}_0\$"),
    legend = false,
    color = :black
);
x = 0.0:last(d) + 3*last(σ_d);
plot!(
    plt,
    x,
    x .* H_0[1],
    ribbon = 3*H_0[2],
    color = :red
);
annotate!(
    plt,
    d[2]*1.05,
    v[2] - 5*σ_v[2],
    "NGC 3627"
);
annotate!(
    plt,
    d[3],
    v[3] - 3*σ_v[3],
    "NGC 3368"
);
annotate!(
    plt,
    d[4],
    v[4] - 3*σ_v[4],
    "NGC 3147"
);
annotate!(
    plt,
    d[2],
    H_0[1]*maximum(x)/2,
    latexstring("\$\\mathrm{H}_0\$  = $(round(H_0[1], digits = 3)) \$\\pm\$ $(round(H_0[2], digits = 5))  km/s*MPsc")
);
display(plt);
#savefig("Hubble.png");