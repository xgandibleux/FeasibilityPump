# ==============================================================================
# Calcul de la relaxation linéaire du SPA

function computeLinearRelaxSPA( A::Array{Int,2}, c::Array{Int,1} )

  nbctr, nbvar = size(A) 
  model = Model(HiGHS.Optimizer)
  @variable(model, 0.0 <= x[1:nbvar] <= 1.0 )
  @objective(model, Min, sum((c[i])*x[i] for i in 1:nbvar))  
  @constraint(model, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
  set_silent(model)
  optimize!(model)
  return objective_value(model), value.(x)

end

# ==============================================================================
# calcul de la solution optimale en {0,1} pour permettre de se situer

function computeSPA( A::Array{Int,2}, c::Array{Int,1} )

  nbctr, nbvar = size(A) 
  model = Model(HiGHS.Optimizer)
  @variable(model, x[1:nbvar], Bin )
  @objective(model, Min, sum((c[i])*x[i] for i in 1:nbvar))  
  @constraint(model, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
  set_silent(model)
  optimize!(model)
  return objective_value(model), value.(x)

end
