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
    cord::Vector{Float64}
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

    function Molecule(_atoms::Vector{atom}, _basis_name::String, _charge=0::Int)
        basis = Vector{orb_detail}()
        all_basis = Tools.get_basis(_basis_name)
        for _atom in _atoms
            append!(basis, format_basis(_atom, all_basis[_atom.symbol]))
        end
        _basis_num = count_basis(basis)
        _elec_num = sum([ELEMENTS_PROTON[_atom.symbol] for _atom = _atoms]) -_charge
        _occ, _occ_num = make_occ(_basis_num, _elec_num)
        new(_atoms, _charge, _basis_name, basis, _basis_num,_elec_num, true, _occ, _occ_num)
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

end