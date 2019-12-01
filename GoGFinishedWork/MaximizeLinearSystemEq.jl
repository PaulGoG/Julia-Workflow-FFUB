using JuMP
using GLPK

model = Model(with_optimizer(GLPK.Optimizer))

@variable(model, 0 <= x1)
@variable(model, 0 <= x2)

@objective(model, Max, x1 + x2)

@constraint(model, con1, x1 + x2 <= 3)
@constraint(model, con2, -x1 + 3*x2 <= 1)
@constraint(model, con3, x2 <= 3)

optimize!(model)
termination_status(model)
primal_status(model)
dual_status(model)

objective_value(model)
value(x1)
value(x2)
