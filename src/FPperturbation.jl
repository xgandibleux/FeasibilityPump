# ==============================================================================
# perturbe une solution binaire

function perturbSolution!( vg::Vector{tGenerateur}, c1::Array{Int,1} )

    k=1
    # nombre de variables maximum a considerer
    T = 20
    # nombre effectif de variables a flipper
    TT = rand(T/2:3*T/2)

    # liste des candidats (valeur, indice) et tri decroissant
    nbvar = length(vg[k].sInt.x)
    candidats=[( abs( vg[k].sPrj.x[i] - vg[k].sInt.x[i] ) , i ) for i=1:nbvar]
    sort!(candidats, rev=true, by = x -> x[1])
    #sort!(candidats,  by = x -> x[1])

    #@show vg[k].sPrj.x
    #@show vg[k].sInt.x
    #@show candidats
    #@show TT
    #@show nbvar

    i = 1
    while (i<= nbvar) && (i<=TT)
        j=candidats[i][2]
#        @show candidats[i][2]
        if vg[k].sInt.x[j] == 0
            vg[k].sInt.x[j] = 1
            vg[k].sInt.y[1] = vg[k].sInt.y[1] + c1[j]
        else
            vg[k].sInt.x[j] = 0
            vg[k].sInt.y[1] = vg[k].sInt.y[1] - c1[j]
        end
        i+=1
    end
    @printf("      perturbing...  : y= %8d \n", convert(Int64,round(vg[k].sInt.y[1], digits=3)) )
    return nothing
end
