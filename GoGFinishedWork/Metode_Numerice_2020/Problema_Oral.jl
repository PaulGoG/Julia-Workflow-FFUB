#Numerical method for finding root of function f(x)

using Plots

function f(x)
    x^2 - 4x + 4 - log(x)
end

@eval macro $(:do)(block, when::Symbol, condition) # Am creeat instructiunea do while deoarece Julia nu o are
           when ≠ :when && error("@do expected `when` got `$s`")
           quote
               let
                   $block
                   while $condition
                       $block
                   end
               end
           end |> esc
       end

a = 2
b = 4
Tol = 10^-6
n = 0
clearconsole()
@do begin
    global x = (a+b)/2
    global n+=1
    println(x)
    if f(a)*f(x) < 0
         global b = x
    else
         global a = x
    end
end when b - a > Tol
println("Solutia aproximativa dupa un numar de " ,n, " iteratii este " ,x,)

plot(f, 0, 5, framestyle =:box, title = "Graficul functiei f(x) = x^2 - 4x + 4 - ln(x)",
label = nothing)
xlabel!("x")
ylabel!("f(x)")
hline!([0], linestyle = :dot, color = :black, label = nothing)
scatter!([x],[0], label = nothing)
