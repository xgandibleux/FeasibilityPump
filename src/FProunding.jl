# ==============================================================================
# rend binaire une solution en arrondissant a 0 ou 1 ses variables fractionnaires 

function roundingSolution!(vg::Vector{tGenerateur}, c::Array{Int,1})

    k=1
    nbvar = length(vg[k].sInt.x)
    vg[k].sInt.y[1] = 0

    # fixe les variables 
    nbVarNonEntiere = 0
    for i in 1:nbvar
        if  isapprox(vg[k].sPrj.x[i] , 0.0, atol=1e-3)
            # variable "entiere" a 0
            vg[k].sInt.x[i] = 0
        elseif isapprox(vg[k].sPrj.x[i] , 1.0, atol=1e-3)
            # variable "entiere" a 1
            #@show vg[k].sPrj.x[i]
            vg[k].sInt.x[i] = 1
            vg[k].sInt.y[1] += c[i]
        else
            nbVarNonEntiere += 1
            if vg[k].sPrj.x[i] < 0.5
                #@show vg[k].sPrj.x[i]
                # variable fractionnaire => fixee a x[i]=0
                vg[k].sInt.x[i] = 0
            else
                # variable fractionnaire => fixee x[i]=1
                #@show vg[k].sPrj.x[i]
                vg[k].sInt.x[i] = 1
                vg[k].sInt.y[1] += c[i]
            end
        end
    end

    #@show findall(v -> v > 0, vg[k].sInt.x)
    @printf("  → rounding...    #varRounded : %4d    →    y= %5d  \n", nbVarNonEntiere, vg[k].sInt.y[1])

end