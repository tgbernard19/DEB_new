"""
Maturity formulations allocate a fraction of 
resources to maturity and reproduction.
"""
abstract type AbstractMaturity end

"""
    maturity!(o, u)
    maturity!(f, o, u)

Allocates reserve drain due to maturity maintenance.
Stores in M state variable if it exists.
"""
function maturity! end
maturity!(o, u) = maturity!(maturity_pars(o), o, u)
maturity!(f::Nothing, o, u) = nothing

"""
    Maturity(j_E_mat_mai, κmat, threshold)

A maturity model seperated to make maturity modeling optional.

$(FIELDDOCTABLE)
"""
@kwdef struct Maturity{TJ,TF,TM} <: AbstractMaturity
    j_E_mat_mai::TJ = 0.001
    κmat::TF = 0.05
    threshold::TM = 1.0
end

# n_N_M::MoMo      | 0.05            | mol*mol^-1      | [0.0, 1.0]   | _    | "N/C use for maturity"
# w_M::GMo         | 25.0            | g*mol^-1        | [0.0, 1.0]   | _    | "Mol-weight of shoot maturity reserve:"

maturity!(f::Maturity, o, u) = begin
    seedset = u[:M] > f.threshold ? u[:M] : zero(u[:M]) 
    mat = κmat(f) * E_ctb(o)
    flux(o)[:M,:gro] = mat - seedset * unit(mat)/unit(seedset)
    mat_mai = f.j_E_mat_mai * tempcorrection(o) * u[:V]
    drain = mat + mat_mai
    reserve_drain!(o, Val(:mat), drain)
end

"""
    κmat(maturity_pars::Maturity)

Kappa parameter for maturity.
"""
κmat(maturity_pars::Maturity) = maturity_pars.κmat
