# ==============================================================================
# Projete xTilde sur le polyedre X du SPA avec norme-L1
# version FP 2005

function Δ2SPA( A::Array{Int,2}, xTilde::Array{Int,1} )

    nbctr,nbvar = size(A)
    idxTilde0, idxTilde1 = split01(xTilde)

    proj = Model(HiGHS.Optimizer)
    @variable(proj, 0.0 <= x[1:length(xTilde)] <= 1.0 )
    @objective(proj, Min, sum(x[i] for i in idxTilde0) + sum((1-x[i]) for i in idxTilde1) )
    @constraint(proj, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
    set_silent(proj)
    optimize!(proj)
    return objective_value(proj), value.(x)
end


# ==============================================================================
# projecte la solution entiere correspondant au generateur k et test d'admissibilite
function projectingSolution!( vg::Vector{tGenerateur}, A::Array{Int,2}, c::Array{Int,1} )

    k=1

    # Projete la solution entiere sur le polytope X 
    fPrj, vg[k].sPrj.x = Δ2SPA(A,vg[k].sInt.x)

    # Nettoyage de vg[k].sPrj.x (reconditionne les 0 et 1 et arrondi les autres valeurs)
    nettoyageSolution!(vg[k].sPrj.x)

    # Teste si la projection est admissible
    if isInteger(vg[k].sPrj.x)

        # sauvegarde de la solution entiere admissible obtenue
        vg[k].sInt.x    = deepcopy(vg[k].sPrj.x)
        vg[k].sInt.y[1] = sum(c[i] for i in findall(v -> v > 0, vg[k].sPrj.x))
        #@show findall(v -> v > 0, vg[k].sPrj.x)
        #@show c[5], c[382], c[1498], sum([c[5], c[382], c[1498]])
        #@show vg[k].sInt.y[1]
        #brol
        vg[k].sFea = true
        @printf("  → projecting...  feasible    : \n")

    else

        vg[k].sFea = false
        @printf("  → projecting...  xxxxxxxx    : \n")

    end

end