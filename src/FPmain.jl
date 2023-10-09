# ==============================================================================

println("""\nAlgorithm "Feasibility Pump" --------------------------------\n""")

const verbose = true

println("-) Activate the required packages\n")
using JuMP, HiGHS, Printf, Random 
verbose ? println("  Fait \n") : nothing


# ==============================================================================

include("FPdatastructures.jl") # types, datastructures and global variables specially defined for FP
include("FPparsers.jl")        # parsers of (bi-objective) instances 
include("FPjumpModels.jl")     # JuMP models for computing relaxed optima of the SPA
include("FProunding.jl")       # Startegies for rounding a LP-solution to a 01-solution
include("FPprojection.jl")     # JuMP models for computing the projection on the polytope of the SPA
include("FPperturbation.jl")   # routines dealing with the perturbation of a solution when a cycle is detected


# ==============================================================================
# Elabore 2 ensembles d'indices selon que xTilde[i] vaut 0 ou 1

function split01(xTilde::Array{Int,1})
   indices0 = (Int64)[]
   indices1 = (Int64)[]
   for i=1:length(xTilde)
       if xTilde[i] == 0
           push!(indices0,i)
       else
           push!(indices1,i)
       end
    end

   return indices0, indices1
end


# ==============================================================================
# test si une solution est admissible en verifiant si sa relaxation lineaire
# conduit a une solution entiere

function isInteger(x::Vector{Float64})
    admissible = true
    i=1
    while admissible && i<=length(x)
        if round(x[i], digits=3)!=0.0 && round(x[i], digits=3)!=1.0
            admissible = false
        end
        i+=1
    end
    return admissible
end


# ==============================================================================
# nettoyage des valeurs des variables d'une solution x relachee sur [0,1]

function nettoyageSolution!(x::Vector{Float64})
    nbvar=length(x)
    for i in 1:nbvar
        if     round(x[i], digits=3) == 0.0
                   x[i] = 0.0
        elseif round(x[i], digits=3) == 1.0
                   x[i] = 1.0
        else
                   x[i] = round(x[i], digits=3)
        end
    end
end


# ==============================================================================
# predicat : verifie si une solution entiere est realisable
function isFeasible(vg::Vector{tGenerateur}, k::Int64)
    #verbose && vg[k].sFea == true ? println("   feasible") : nothing
    return (vg[k].sFea == true)
end


# ==============================================================================
# predicat : verifie si le nombre d'essai maximum a ete tente
function isFinished(nIT::Int64, maxnIT::Int64)
#    verbose && nIT > maxnIT ? println("   maxnIT") : nothing
    return (nIT > maxnIT)
end


# ==============================================================================
# predicat : verifie si le budget de calcul maximum a ete consomme
function isTimeout(temps, maxTime)
#    verbose && time()- temps > maxTime ? println("   maxTime") : nothing
    return (time()- temps > maxTime)
end


# ==============================================================================
# point d'entree principal

function FP( fname::String, maxnIT::Int64, maxTime::Int64 )

    # =========================================================================

    # affiche la resolution a traiter -----------------------------------------
    println("\n0) instance = $fname | maxnIT = $maxnIT | maxTime = $maxTime\n") 

    # chargement de l'instance numerique --------------------------------------
    c1, c2, A = loadInstance2SPA(fname) # instance numerique de SPA (biobjectif->utiliser c1)
    nbctr,nbvar = size(A)
    nbobj = 1

    # calcul de la solution optimale en {0,1} pour permettre de se situer
    #fSPA, xSPA = computeSPA(A,c1)
    #@printf("   f_IP  = %8.2f   \n",round( fSPA, digits=3))
    #println(" |  ",  findall(v -> v > 0, value.(xSPA)))


    # =========================================================================

    # ligne 1: initialise le nombre d'iteration et calcule la valeur optimale relachee
    nIT = 0

    # declenche le compteur de temps ecoule
    temps = time()

    fRL, xfRL = computeLinearRelaxSPA(A, c1)
    nettoyageSolution!(xfRL)
    @printf("   f_LP  = %8.2f   \n",fRL)

    # ligne 2: test si la solution relachee est entiere -----------------------
    if isInteger(xfRL)
        println("\n  Integer! ")
        println("   nIT    =  ",nIT,"  |  Tps=", round(time()- temps, digits=4))
        #println("   x[i]   =  1 for i∈", findall(v -> v > 0, xfRL))
        println("   y=f(x) =  ",  fRL,"\n")
        return fRL, findall(v -> v > 0, xfRL)
    end
        
    # allocation de memoire pour la structure de donnees ----------------------
    vg = allocateDatastructure(1, nbvar, nbobj)

    # place la solution relachee dans la stru de donnees  ---------------------
    k = 1   
    ajouterX0!(vg, (tSolution{Float64})(xfRL,[fRL]))
    vg[k].sFea = false

    t2::Bool = false
    t3::Bool = false

    # initialise l'historique des solutions entieres non admissibles visitees
    H = []

    # ligne 3: rounding solution : met a jour sInt dans vg --------------------
    roundingSolution!(vg,c1)  
    
    # ajoute a l'historique la solution ---------------------------------------
    push!(H, findall(v -> v > 0, vg[k].sInt.x))

    # ligne 4: iteration de "pump" tant qu'une condition d'arret n'est pas rencontree
    while  !(t2=isFinished(nIT, maxnIT)) && !(t3=isTimeout(temps, maxTime))

        # ligne 5: passe a l'iteration suivante et calcule la projection
        nIT+=1

        # projecting solution : met a jour sPrj, sInt, sFea dans vg -------
        projectingSolution!(vg,A,c1)
        
        println("    t=",nIT,"  |  Tps=", round(time()- temps, digits=4))

        if isFeasible(vg,k)
            println("\n  Integer! ")
            #println("   x[i]   =  1 for i∈", findall(v -> v > 0, vg[k].sInt.x))
            println("   y=f(x) =  ",  vg[k].sInt.y[1])
            return vg[k].sInt.y[1], findall(v -> v > 0, vg[k].sInt.x)
        end

        # ligne 6 : rounding solution : met a jour sInt dans vg et test detection cycle sur solutions entieres 
        roundingSolution!(vg,c1)
        cycle = findall(v -> v > 0, vg[k].sInt.x) in H
        if (cycle == false)
            # rien a faire
        else
            println("      >>> CYCLE <<<")
            # perturbe la solution -------------------------------------
            perturbSolution!(vg,c1)
        end

        # archive dans l'historique la solution 
        push!(H, findall(v -> v > 0, vg[k].sInt.x))

    end

    println("\n  Failed ")
    if  t2
        println("    stop: maxnIT \n")
    elseif t3
        println("    stop: maxTime \n")
    end

end

# ==============================================================================

#@time FP("sppaa02.txt", 20, 20)
#@time FP("sppnw03.txt", 20, 20) #pb glpk
#@time FP("sppnw10.txt", 20, 20)
#@time FP("didactic3.txt", 5, 10)
#@time FP("didactic5.txt", 5, 10)
#@time FP("sppnw12.txt", 11, 20)

@time FP("didactic3.txt", 5, 10)

println("\n=======================================================")

for i in [1,2,3,5]
    fname = "sppaa0"*string(i)*".txt"
    @time FP(fname, 10, 20)
    println("\n=======================================================")
end


for i in 1:43
    if i<10
        fname = "sppnw0"*string(i)*".txt"
        @time FP(fname, 10, 20)
    else
        fname = "sppnw"*string(i)*".txt"
        @time FP(fname, 10, 20)
    end
    println("\n=======================================================")
end


#=
fname = "rail4284.txt"
@time FP(fname, 5, 10)
=#

nothing
