module Tools
export orb_info, get_basis

mutable struct orb_info
    orb_l::String
    α_array::Array{Float64,1}
    d_array::Array{Array{Float64, 1}, 1}
    orb_info(orb_l::AbstractString, d_num::Int) = new(orb_l, Array{Float64, 1}[], Array{Array{Float64, 1}, 1}[[] for _= 1:d_num])
end

function get_basis(basisname::String)
    bs_dir = @__DIR__
    path = bs_dir *"/"*basisname
    l = []
    open("$(path)", "r") do fin
        l = [strip(s) for s = readlines(fin)]
    end
    basis = Dict{String, Array{orb_info, 1}}()
    flag = false
    i = 1
    while i +1 <= length(l)
        i += 1
        if l[i] == "BASIS \"ao basis\" PRINT"
            flag = true
            continue
        end
        if !flag || l[i] =="END"
            continue
        end
        if startswith(l[i], "#")
            i += 1
        end
        symbol = split(l[i])
        if !(symbol[1] in keys(basis))
            basis[symbol[1]] = Array{orb_info, 1}[]
        end
        alpha_d_num = length(split(l[i+1]))
        orb_instance = orb_info(symbol[2], alpha_d_num-1)
        while i + 1 <= length(l)
            i += 1
            l1 = split(l[i])
            if startswith(l1[1], "#") || l1[1] == "END" || isletter(l1[1][1])
                i -= 1
                break
            end
            for (j, s) = enumerate(split(l[i]))
                s = parse(Float64, s)
                if j ==1
                    push!(orb_instance.α_array, s)
                else
                    push!(orb_instance.d_array[j-1], s)
                end 
            end
        end
        push!(basis[symbol[1]], orb_instance)
    end
    return basis
end

#basis = get_basis("sto3g")
#@show basis["H"]
end

