using Plots; plotlyjs()

x=y = collect(-1000:1:1000.0)
Kek = tan(20*pi/180)
Q = zeros(length(x),length(y))
for i in 1:length(x)
    for j in 1:length(y)
        if y[j]/x[i] <= Kek && y[j]/x[i] > 0 && abs(y[j]) > 0
            Q[i,j] = atan(y[j]/x[i])
        end
    end
end
ok = Tuple.(findall(x->x!=0.0,Q))
a = getindex.(ok,1)
b = getindex.(ok,2)
c = zeros(length(a))
for i in 1:length(x)
    for j in 1:length(y)
        if findfirst(x->x==tuple(i,j),ok) !== nothing
            c[findfirst(x->x==tuple(i,j),ok)] = Q[i,j]
        end
    end
end
plot(x[a],y[b],c)
heatmap(x[a],y[b],c)
d = hcat(x[a],y[b])
for i in 1:1000
    q = d[rand(Tuple(collect(1:length(c)))),:]
    if abs(q[1]) <= 100
        println(q)
    end
end