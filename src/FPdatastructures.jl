# types ------------------------------------------------------------------------

# type corresponding to a solution
mutable struct tSolution{T}
    x :: Vector{T}                # vector variables x (1..n)
    y :: Vector{T}                # vector outcomes  y (1..p)
end

# type corresponding to a point generator
mutable struct tGenerateur
    sRel :: tSolution{Float64}    # initial relaxed solution
    sInt :: tSolution{Int64}      # integer solution
    sPrj :: tSolution{Float64}    # projected solution
    sFea :: Bool                  # indicate if sInt is feasible or not
end


# type grouping the lists of points for dysplaying purposes
mutable struct tListDisplay
    xLf1  :: Vector{Float64};  yLf1  :: Vector{Float64} # liste des points (x,y) relaches
    xLf2  :: Vector{Float64};  yLf2  :: Vector{Float64} # liste des points (x,y) relaches
    xL    :: Vector{Float64};  yL    :: Vector{Float64} # liste des points (x,y) relaches
    XInt  :: Vector{Int64};    YInt  :: Vector{Int64}   # liste des points (x,y) entiers
    XProj :: Vector{Float64};  YProj :: Vector{Float64} # liste des points (x,y) projetes
    XFeas :: Vector{Int64};    YFeas :: Vector{Int64}   # liste des points (x,y) admissibles
    XPert :: Vector{Int64};    YPert :: Vector{Int64}   # liste des points (x,y) perturbes
end


# ==============================================================================
# Initialisation structure donnees contenant tous les generateurs

function allocateDatastructure(nbgen::Int64, nbvar::Int64, nbobj::Int64)

    verbose ? println("  → Allocation memoire pour ",nbgen," run(s) FP\n") : nothing

    vg = Vector{tGenerateur}(undef, nbgen)
    for k = 1:nbgen
        vg[k] = tGenerateur(  tSolution{Float64}(zeros(Float64,nbvar),zeros(Float64,nbobj)),
                              tSolution{Int64}(zeros(Int64,nbvar),zeros(Int64,nbobj)),
                              tSolution{Float64}(zeros(Float64,nbvar),zeros(Float64,nbobj)),
                              false
                            )
    end
    return vg
end


# ==============================================================================
# sauvegarde d'une solution relachee initiale

function ajouterX0!(vg::Vector{tGenerateur}, s::tSolution{Float64})
    k=1
    vg[k].sRel = deepcopy(s) # sauvegarde la solution relachee \bar{x}^{k}
    vg[k].sPrj = deepcopy(s) # sauvegarde la premiere projection \bar{x}^{k,0}
    return nothing
end


# ==============================================================================
# sauvegarde d'une solution entiere (arrondie ou perturbee) 

function ajouterXtilde!(vg::Vector{tGenerateur}, x::Vector{Int64}, y::Vector{Int64})
    k=1
    vg[k].sInt.x = deepcopy(x)
    vg[k].sInt.y = deepcopy(y)
    return nothing
end


# ==============================================================================
# sauvegarde d'une solution fractionnaire (projetee) 

function ajouterXbar!(vg::Vector{tGenerateur}, x::Vector{Float64}, y::Vector{Float64})
    k=1
    vg[k].sPrj.x = deepcopy(x)
    vg[k].sPrj.y = deepcopy(y)
    return nothing
end
