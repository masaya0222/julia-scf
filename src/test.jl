module test
    using JuliaSCF.Mole
    using JuliaSCF.Integral
    using SpecialFunctions
    m = Molecule([atom("I",[0,0,-0.7]) ,atom("I",[0,0,+0.7]),atom("I",[0,0,-1.4]) ,atom("I",[0,0,+1.4])], "sto3g")
    @time S = Int1e_ovlp.get_ovlp(m)
    
    @show SpecialFunctions.gamma(2)
    v = Vector{Array{Int,2}}()
    a = [1 2;3 4]
    push!(v,a)
    @show v

end