module UHF

using LinearAlgebra
using JuliaSCF.Mole
using JuliaSCF.Integral
using JuliaSCF.DIIS
using JuliaSCF.HF

export uhf, run_uhf!

mutable struct uhf
    mol::Molecule
    ovlp_ao::Matrix{Float64}
    kin_ao::Matrix{Float64}
    v1e_ao::Matrix{Float64}
    v2e_ao::Array{Float64,4}
    hcore::Matrix{Float64}
    mo_num::Int
    coeff_a::Matrix{Float64}
    coeff_b::Matrix{Float64}
    mo_ene_a::Vector{Float64}
    mo_ene_b::Vector{Float64}
    dense_ao_a::Matrix{Float64}
    dense_ao_b::Matrix{Float64}
    fock_ao_a::Matrix{Float64}
    fock_ao_b::Matrix{Float64}
    elec_ene::Float64
    nuc_ene::Float64
    total_ene::Float64
    diis_a::Diis
    diis_b::Diis

    function uhf(mol::Molecule, init=true::Bool)
        ovlp_ = Int1e_ovlp.get_ovlp(mol)
        kin_ = Int1e_kin.get_kin(mol)
        v1e_ = Int1e_nuc.get_v1e(mol)
        v2e_ = Int2e.get_v2e(mol)
        hcore_ = kin_+v1e_
        mo_num_ = mol.basis_num
        mo_ene_a_, coef_a_ = eigenh(hcore_,ovlp_)
        mo_ene_b_, coef_b_ = eigenh(hcore_,ovlp_)
 
        dm = init_guess(mol) 
        dense_ao_a_ = zero(dm) #get_dense_ao(mo_num_, coef_a_, mol.occ_num, mol.occ/2)
        dense_ao_a_[begin:div(mo_num_,2),begin:div(mo_num_,2)] = dm[begin:div(mo_num_,2),begin:div(mo_num_,2)]
        dense_ao_b_ = zero(dm) #get_dense_ao(mo_num_, coef_b_, mol.occ_num, mol.occ/2)
        dense_ao_b_[begin+div(mo_num_,2):end,begin+div(mo_num_,2):end] = dm[begin+div(mo_num_,2):end,begin+div(mo_num_,2):end] 

        fock_ao_a_, fock_ao_b_ = get_fock_ao(hcore_, mo_num_, dense_ao_a_, dense_ao_b_, v2e_)
        
        new(mol, ovlp_, kin_, v1e_, v2e_, hcore_, mo_num_, coef_a_, coef_b_,
            mo_ene_a_, mo_ene_b_, dense_ao_a_, dense_ao_b_, fock_ao_a_, fock_ao_b_, Inf, get_nuc_ene(mol), Inf, Diis(mol), Diis(mol))
        
    end

end

function get_dense_ao(len::Int, coef::Matrix{Float64}, occ_num::Int, occ::Vector{Float64})
    D = zeros(Float64, (len, len))
    for p = 1:len
        for q = 1:len
            for a = 1:occ_num
                D[p,q] += occ[a]*coef[p,a]*coef[q,a]
            end
            #D[p,q] = 2*dot(coef[p,begin:occ_num]',coef[q,begin:occ_num])
        end
    end
    return D
end

function get_dense_ao!(F::uhf)
    F.dense_ao_a = get_dense_ao(F.mo_num, F.coeff_a, F.mol.occ_num, F.mol.occ/2)
    F.dense_ao_b = get_dense_ao(F.mo_num, F.coeff_b, F.mol.occ_num, F.mol.occ/2)
    return F.dense_ao_a, F.dense_ao_b
end

function get_fock_ao(hcore::Matrix{Float64}, mo_num::Int, dense_ao_a::Matrix{Float64}, dense_ao_b::Matrix{Float64}, v2e_ao::Array{Float64, 4})
    f_a = copy(hcore)
    f_b = copy(hcore)
    for mu = 1:mo_num, nu = 1:mo_num, p = 1:mo_num, q = 1:mo_num
        f_a[mu,nu] += (dense_ao_a[p,q]+dense_ao_b[p,q])*v2e_ao[mu,nu,q,p]-dense_ao_a[p,q]*v2e_ao[mu,p,q,nu]
        f_b[mu,nu] += (dense_ao_a[p,q]+dense_ao_b[p,q])*v2e_ao[mu,nu,q,p]-dense_ao_b[p,q]*v2e_ao[mu,p,q,nu]
    end
    
    return f_a, f_b
end

function get_fock_ao!(F::uhf)
    F.fock_ao_a, F.fock_ao_b= get_fock_ao(F.hcore, F.mo_num, F.dense_ao_a, F.dense_ao_b, F.v2e_ao)
    return F.fock_ao_a, F.fock_ao_b
end

function get_nuc_ene(mol::Molecule)
    potential = 0.0
    for i = 1:length(mol.atoms)
        for j = i+1:length(mol.atoms)
            distance = norm(mol.atoms[i].cord-mol.atoms[j].cord)
            if distance <= 1e-10
                println("atom is arranged in very close")
            end
            potential += mol.atoms[i].atomic_num*mol.atoms[j].atomic_num/distance
        end
    end
    return potential
end

function converged(F::uhf, f_a::Matrix{Float64}, f_b::Matrix{Float64}, threshold::Float64)
    error_vector_a = f_a[begin:F.mol.occ_num,begin+F.mol.occ_num:end][:]
    error_a = norm(error_vector_a)

    error_vector_b = f_b[begin:F.mol.occ_num,begin+F.mol.occ_num:end][:]
    error_b = norm(error_vector_a)
    
    #@show error
    return error_b < threshold && error_b < threshold
end
function eigenh(f::Matrix{Float64}, s::Matrix{Float64})
    e,c = eigen(f,s)
    c = c[:,sortperm(e)]
    sort!(e)
    n = c'*s*c
    for i = 1:length(e)
        c[:,i] /= sqrt(n[i,i])
    end

    for i = 1:length(e)
        m = argmax(abs.(c[:,i]))
        if c[m,i] < 0
            c[:,i] *= -1
        end
    end
    
    return e,c
end

function run_uhf!(F::uhf)
    max_iteration = 100
    diis_on = true
    if F.mo_num -F.mol.occ_num <= 1
        diis_on = false
    end
    num = 0
    for i = 1:max_iteration
        #@show i
        get_fock_ao!(F)
        if diis_on
            insert(F.diis_a, F.fock_ao_a, F.coeff_a)
            insert(F.diis_b, F.fock_ao_b, F.coeff_b)
            F.fock_ao_a = return_fock(F.diis_a)
            F.fock_ao_b = return_fock(F.diis_b)
        end
        fock_mo_a = transpose(F.coeff_a)*F.fock_ao_a*F.coeff_a
        fock_mo_b = transpose(F.coeff_b)*F.fock_ao_b*F.coeff_b
        F.mo_ene_a, F.coeff_a = eigenh(F.fock_ao_a, F.ovlp_ao)
        F.mo_ene_b, F.coeff_b = eigenh(F.fock_ao_b, F.ovlp_ao)
        old_dense_a = copy(F.dense_ao_a)
        old_dense_b = copy(F.dense_ao_b)
        get_dense_ao!(F)
        if isapprox(old_dense_a, F.dense_ao_a) && isapprox(old_dense_b, F.dense_ao_b)#converged(F, fock_mo_a, fock_mo_b, 1e-7) && i > 1
            println("Converged")
            @show i
            break
        end
        num = i
    end
    if num == max_iteration
        println("Not Convergence")
    end
    #@show F.mo_ene
    F.elec_ene = 1/2*tr(F.dense_ao_a*(F.hcore+F.fock_ao_a)+F.dense_ao_b*(F.hcore+F.fock_ao_b))
    #println("HF elec energy")
    #println(F.elec_ene)
    F.total_ene = F.elec_ene + F.nuc_ene
    @show F.total_ene
    return F.total_ene
end

end
