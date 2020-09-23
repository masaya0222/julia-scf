module DIIS
using JuliaSCF.Mole
using LinearAlgebra
export Diis, insert, return_fock
mutable struct Diis
    ind::Int
    store_max::Int
    store_num::Int
    mo_num::Int
    occ_num::Int
    store_fock::Vector{Matrix{Float64}}
    store_error::Vector{Vector{Float64}}
    function Diis(mol::Molecule)
        new(1, 5, 0, mol.basis_num, mol.occ_num, Matrix{Float64}[], Vector{Float64}[])
    end
end

function return_error(D::Diis, fock_ao::Matrix{Float64}, coef::Matrix{Float64})
    fock_mo = coef'*fock_ao*coef
    error_vector = fock_mo[begin:D.occ_num,begin+D.occ_num:end][:]
    return error_vector
end

function insert(D::Diis, fock::Matrix{Float64}, coef::Matrix{Float64})
    if D.store_num < D.store_max
        push!(D.store_fock, copy(fock))
        push!(D.store_error, return_error(D, fock, coef))
        D.store_num += 1
    else
        D.store_fock[D.ind] = copy(fock)
        D.store_error[D.ind] = return_error(D, fock, coef)
    end
    D.ind = mod(D.ind+1, D.store_max)+1
end

function return_fock(D::Diis)
    B = ones(Float64, (D.store_num+1, D.store_num+1))*-1
    B[D.store_num+1,D.store_num+1] = 0.0
    for i = 1:D.store_num
        for j = 1:D.store_num
            B[i,j] = D.store_error[i]'*D.store_error[j]
        end
    end
    b = zeros(Float64, D.store_num+1)
    b[D.store_num+1] = -1.0
    w = b'*pinv(B)
    new_fock = sum([D.store_fock[i]*w[i] for i = 1:D.store_num])
    return new_fock
end
end