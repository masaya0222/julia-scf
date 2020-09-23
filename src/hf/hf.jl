module HF

using LinearAlgebra
using JuliaSCF.Mole
using JuliaSCF.Integral
using JuliaSCF.DIIS

export hf, run!

mutable struct hf
    mol::Molecule
    ovlp_ao::Matrix{Float64}
    kin_ao::Matrix{Float64}
    v1e_ao::Matrix{Float64}
    v2e_ao::Array{Float64,4}
    hcore::Matrix{Float64}
    mo_num::Int
    coeff::Matrix{Float64}
    mo_ene::Vector{Float64}
    dense_ao::Matrix{Float64}
    fock_ao::Matrix{Float64}
    elec_ene::Float64
    nuc_ene::Float64
    total_ene::Float64
    diis::Diis


    function hf(mol::Molecule, init=true::Bool)
        ovlp_ = Int1e_ovlp.get_ovlp(mol)
        kin_ = Int1e_kin.get_kin(mol)
        v1e_ = Int1e_nuc.get_v1e(mol)
        v2e_ = Int2e.get_v2e(mol)
        hcore_ = kin_+v1e_
        mo_num_ = mol.basis_num
        mo_ene_, coef_ = eigenh(hcore_,ovlp_)
        if init
            dense_ao_ = init_guess(mol)
        else
            dense_ao_ = get_dense_ao(mo_num_, coef_, mol.occ_num, mol.occ)
        end
        fock_ao_ = get_fock_ao(kin_+v1e_,mo_num_, dense_ao_, v2e_)

        new(mol, ovlp_, kin_, v1e_, v2e_, kin_+v1e_, mo_num_, coef_, 
            mo_ene_, dense_ao_, fock_ao_, Inf, get_nuc_ene(mol), Inf, Diis(mol))
    end
end

function init_guess(mol::Molecule)
    dm = zeros(Float64, (mol.basis_num, mol.basis_num))
    ind = 1
    for i = 1:length(mol.atoms)
        a = mol.atoms[i]
        mol_atom = Molecule([atom(a.symbol, [0,0, 0,0,0.0])], mol.basis_name)
        mol_atom.occ, mol_atom.occ_num = make_atom_occ(mol_atom)
        kernel = hf(mol_atom, false)
        run_silence!(kernel)
        dm[ind:ind+mol_atom.basis_num-1,ind:ind+mol_atom.basis_num-1,] = kernel.dense_ao
        ind += mol_atom.basis_num
    end
    return dm
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

function get_dense_ao!(F::hf)
    F.dense_ao = get_dense_ao(F.mo_num, F.coeff, F.mol.occ_num, F.mol.occ)
    return F.dense_ao
end

function get_fock_ao(hcore::Matrix{Float64}, mo_num::Int, dense_ao::Matrix{Float64}, v2e_ao::Array{Float64, 4})
    f = copy(hcore)
    for mu = 1:mo_num, nu = 1:mo_num, p = 1:mo_num, q = 1:mo_num
        f[mu,nu] += dense_ao[p,q]*(v2e_ao[mu,nu,q,p]-0.5*v2e_ao[mu,p,q,nu])
    end
    
    return f
end

function get_fock_ao!(F::hf)
    F.fock_ao = get_fock_ao(F.hcore, F.mo_num, F.dense_ao, F.v2e_ao)
    return F.fock_ao
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

function converged(F::hf, f::Matrix{Float64}, threshold::Float64)
    error_vector = f[begin:F.mol.occ_num,begin+F.mol.occ_num:end][:]
    error = norm(error_vector)
    #@show error
    return error < threshold
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

function run_silence!(F::hf)
    max_iteration = 100
    for i = 1:max_iteration
        F.fock_ao = get_fock_ao!(F)
        insert(F.diis, F.fock_ao, F.coeff)
        F.fock_ao = return_fock(F.diis)
        fock_mo = transpose(F.coeff)*F.fock_ao*F.coeff
        F.mo_ene, F.coeff = eigenh(F.fock_ao, F.ovlp_ao)
        get_dense_ao!(F)
        if converged(F, fock_mo, 1e-7)
            break
        end
    end
    F.elec_ene = 1/2*tr(F.dense_ao*(F.hcore+F.fock_ao))
    F.total_ene = F.elec_ene + F.nuc_ene
    return F.total_ene
end

function run!(F::hf)
    max_iteration = 100
    for i = 1:max_iteration
        #@show i
        F.fock_ao = get_fock_ao!(F)
        insert(F.diis, F.fock_ao, F.coeff)
        F.fock_ao = return_fock(F.diis)
        fock_mo = transpose(F.coeff)*F.fock_ao*F.coeff
        F.mo_ene, F.coeff = eigenh(F.fock_ao, F.ovlp_ao)
        get_dense_ao!(F)
        if converged(F, fock_mo, 1e-7)
            break
        end
    end
    #@show F.mo_ene
    F.elec_ene = 1/2*tr(F.dense_ao*(F.hcore+F.fock_ao))
    #println("HF elec energy")
    #println(F.elec_ene)
    F.total_ene = F.elec_ene + F.nuc_ene
    @show F.total_ene
    return F.total_ene
end
#=
i = 4.0
m1 = Molecule([atom("O",[0.0,0.0,-i/2]),atom("O",[0.0,0.0,i/2])],"sto3g")
#kernel1 = hf(m1)
#run!(kernel1)

using PyCall
@pyimport pyscf
@pyimport numpy
X = 0.52918
mol_H2 = pyscf.gto.Mole()
mol_H2.build(atom="O 0 0 $(-i/2*X); O 0 0 $(i/2*X) ", basis="sto3g")

kernel1 = hf(m1)
kernel1.dense_ao = pyscf.scf.hf.get_init_guess(mol_H2)
kernel1.v2e_ao = mol_H2.intor("int2e")
kernel1.hcore = pyscf.scf.hf.get_hcore(mol_H2)
run!(kernel1)

m = pyscf.scf.RHF(mol_H2)
m.kernel()

@show pyscf.scf.hf.energy_tot(m,dm = kernel1.dense_ao, h1e = kernel1.hcore, vhf = pyscf.scf.hf.get_veff(mol_H2, kernel1.dense_ao))
@show m.mo_energy
=#
end