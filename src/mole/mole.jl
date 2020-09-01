module Mole
#include("../basis/tools.jl")
include("element_data.jl")
using JuliaSCF.Tools
using .Element: ELEMENTS_PROTON

export orb_detail, atom, Molecule

struct orb_detail 
    cord::Vector{Float64}
    orb_l::Int
    α_array::Vector{Float64}
    d_array::Vector{Vector{Float64}}
end

struct atom
    symbol::String
    atomic_num::Int
    cord::Vector{Float64}
    function atom(_symbol::String, _cord::Vector{Float64})
        new(_symbol, ELEMENTS_PROTON[_symbol], _cord)
    end
end

mutable struct Molecule
    atoms::Vector{atom}
    charge::Int
    basis_name::String
    basis::Vector{orb_detail}
    basis_num::Int
    elec_num::Int
    can_rhf::Bool
    occ::Vector{Int}
    occ_num::Int
    rotate_coef::Vector{Array{Float64, 4}}

    function Molecule(_atoms::Vector{atom}, _basis_name::String, _charge=0::Int)
        basis = Vector{orb_detail}()
        all_basis = Tools.get_basis(_basis_name)
        for _atom in _atoms
            append!(basis, format_basis(_atom, all_basis[_atom.symbol]))
        end
        _basis_num = count_basis(basis)
        _elec_num = sum([ELEMENTS_PROTON[_atom.symbol] for _atom = _atoms]) -_charge
        _occ, _occ_num = make_occ(_basis_num, _elec_num)
        _rotate_coef = make_rotate_coef(basis)
        new(_atoms, _charge, _basis_name, basis, _basis_num,_elec_num, true, _occ, _occ_num, _rotate_coef)
    end
end

function make_occ(basis_num::Int, elec_num::Int, )
    occ = zeros(Int, basis_num)
    occ_num = div(elec_num,2)
    for i = 1:occ_num
        occ[i] = 2
    end
    if elec_num % 2 == 1
        occ_num += 1
        occ[occ_num] = 1
    end
    return (occ, occ_num)
end

function format_basis(_atom::atom, _orb_info_vec::Vector{Tools.orb_info})
    convert_l = Dict([('S', 0), ('P', 1), ('D', 2)])
    res_s = Vector{orb_detail}()
    res_p = Vector{orb_detail}()
    res_d = Vector{orb_detail}()
    res_spd = [res_s, res_p, res_d]
    for b = _orb_info_vec
        if b.orb_l =="SP"
            for (i,l) = enumerate(b.orb_l)
                push!(res_spd[convert_l[l]+1], orb_detail(_atom.cord, convert_l[l], b.α_array, [b.d_array[i]]))
            end
        else
            push!(res_spd[convert_l[b.orb_l[1]]+1], orb_detail(_atom.cord, convert_l[b.orb_l[1]], b.α_array, b.d_array))
        end
    end
    res = Vector{orb_detail}()
    for r = res_spd
        if !isempty(r)
            append!(res, r)
        end
    end
    return res
end

function count_basis(_basis::Vector{orb_detail})
    num = 0

    for b = _basis
        num += (b.orb_l*2+1)*length(b.d_array)
    end
    return num
end
function make_rotate_coef(_basis::Vector{orb_detail})
    max_l = 0
    res = Vector{Array{Float64, 4}}()
    for b = _basis
        max_l = max(max_l, b.orb_l)
    end

    comb = zeros(Int64, (max_l+1, max_l+1))
    for i = 0:max_l, j = 0:max_l
        comb[i+1,j+1] = binomial(i,j)
    end

    for l = 0:max_l
        C = zeros(Float64, (l+1, div(l,2)+1, div(l,2)+1, l+1))
        for m_ = 0:l, t = 0:div(l,2), u = 0:div(l,2), v = 0:l
            if l >=t >=u && l-t>=m_+t && m_>=v
                C[m_+1,t+1,u+1,v+1] = (-2(t%2)+1.0)*(1/4^t)*comb[l+1,t+1]*comb[l-t+1,m_+t+1]*comb[t+1,u+1]*comb[m_+1,v+1]
            end
        end
        push!(res, C)
    end
    return res
end

end