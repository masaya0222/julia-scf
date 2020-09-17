module test
    using JuliaSCF.Mole
    using JuliaSCF.Integral
    using SpecialFunctions

    @show gamma

    a = [1,2,3,4,5]
    b  = view(a,2:length(a))
    @show b
    b[1] = 3
    @show b 
    @show a
    m = Molecule([atom("C",[0.0,0.0,-0.7]),atom("C",[0.0,0.0,0.7]),atom("C",[0.0,-0.7,0.0]),atom("C",[0.0,0.7,0.0]), atom("C",[-0.7,0.0,0.0]),atom("C",[0.7,0.0,0.0]), atom("C",[0.0,0.0,-1.4]),atom("C",[0.0,0.0,1.4])], "sto3g")
    @time v = Int2e.get_v2e(m)
    @show size(v)
end