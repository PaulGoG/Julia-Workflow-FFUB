# y'(x) = sqrt(abs(sin(y(x)))) & y(0) = 0

using Plots

function f(y)
    return sqrt(abs(sin(y)));
end

function RungeKutta4(y, dx)
    k1 = f(y);
    k2 = f(y + dx*k1/2);
    k3 = f(y + dx*k2/2);
    k4 = f(y + dx*k3);

    return y + (1/6)*dx*(k1+k2*2+k3*2+k4);
end

N = 1000;
a = 0;
b = 100;
dx = abs(b-a)/(N-1);
x = collect(a:dx:b);
y = zeros(N);

y[1] = 1;

for i in 1:(N-1)
    y[i+1] = RungeKutta4(y[i], dx)
end

plot(x,y)