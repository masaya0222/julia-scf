module test
    using JuliaSCF.Mole
    using JuliaSCF.Integral
    using SpecialFunctions
    
    #m = Molecule([atom("C",[0.0,0.0,-0.7]), atom("C",[0.0,0.0,0.7]), atom("C",[0.0,0.0,1.4]), atom("C",[0.0,1.4,0.0]), atom("C",[1.4,0.0,0.0])], "sto3g")
    m = Molecule([atom("I",[0.0,0.0,-0.7]),atom("I",[0.0,0.0,0.7])],"sto3g")
    #@time v = Int2e.get_v2e(m)
    @time v1 = Int2e.get_v2e_single(m)
    using PyCall
    @pyimport pyscf
    X = 0.52918
    mol = pyscf.gto.Mole()
    mol.build(atom="I 0 0 $(-0.7*X); I 0 0 $(0.7*X)", basis="sto3g")
    #mol.build(atom="C 0 0 $(-0.7*X); C 0 0 $(0.4*X); C 0 $(0.4*X) 0; C 0 $(-0.4*X) 0; C $(-0.4*X) 0 0", basis="sto3g")
    @time v = mol.intor("int2e")
    
    
end