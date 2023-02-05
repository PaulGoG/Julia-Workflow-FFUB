#=
y'(x) = z(x)
z'(x) = -sin(y(x))*y(x)
y(0) = 1 & z(0) = 0
=#

using Plots

function f(y, z)
    return [z, -sin(y)*y];
end

function RungeKutta4(y, z, dx)
    k1 = f(y, z);
    k2 = f(y + dx*k1[1]/2, z + dx*k1[2]/2);
    k3 = f(y + dx*k2[1]/2, z + dx*k2[2]/2);
    k4 = f(y + dx*k3[1], z + dx*k3[2]);

    new_y = y + (1/6)*dx*(k1[1] + k2[1]*2 + k3[1]*2 + k4[1]);
    new_z = z + (1/6)*dx*(k1[2] + k2[2]*2 + k3[2]*2 + k4[2]);

    return [new_y, new_z];
end

N = 1000;
a = 0;
b = 100;
dx = abs(b-a)/(N-1);
x = collect(a:dx:b);
y = zeros(N);
z = zeros(N);

y[1] = 1;
z[1] = 0;

for i in 1:(N-1)
    y[i+1] = RungeKutta4(y[i], z[i], dx)[1]
    z[i+1] = RungeKutta4(y[i], z[i], dx)[2]
end

plot(x, y)