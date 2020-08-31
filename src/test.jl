module test
    using JuliaSCF.Mole
    using JuliaSCF.Integral
    m = Molecule([atom("H",[0,0,-0.7]) ,atom("H",[0,0,+0.7])], "sto3g")
    S = Int1e_ovlp.get_ovlp(m)
    @show S
end