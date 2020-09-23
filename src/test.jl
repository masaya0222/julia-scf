module test
    using JuliaSCF.Mole
    using JuliaSCF.Integral
    using SpecialFunctions
    using LinearAlgebra
    using JuliaSCF.HF
    using Plots

    
    #m = Molecule([atom("C",[0.0,0.0,-0.7]), atom("C",[0.0,0.0,0.7]), atom("C",[0.0,0.0,1.4]), atom("C",[0.0,1.4,0.0]), atom("C",[1.4,0.0,0.0])], "sto3g")
    m = Molecule([atom("C",[0.0,0.0,-0.7]),atom("C",[0.0,0.0,0.7])],"ccpvdz")
    #@time v = Int2e.get_v2e_single(m)
    #@time v = Int2e.get_v2e(m)
    
    #=
    using PyCall
    @pyimport pyscf
    X = 0.52918
    mol = pyscf.gto.Mole()
    mol.build(atom="C 0 0 $(-0.8*X); C 0 0 $(0.9*X); C 0 $(-0.4*X) 0; C 0 $(0.7*X) 0", basis="ccpvdz")
    =#
    #mol.build(atom="C 0 0 $(-0.7*X); C 0 0 $(0.4*X); C 0 $(0.4*X) 0; C 0 $(-0.4*X) 0; C $(-0.4*X) 0 0", basis="sto3g")
    #@time v = mol.intor("int2e")
    gr()
    
    
    enes1 = []
    enes2 = []
    dis = []
    #=
    for i in 1.0:0.1:5
        @show i
        #=
        m1 = Molecule([atom("O",[0.0,0.0,-i/2]),atom("O",[0.0,0.0,i/2])],"sto3g")
        kernel1 = hf(m1)
        ene1 = run!(kernel1)
        =#
        mol = pyscf.gto.Mole()
        mol.build(atom="O 0 0 $(-i/2*X); O 0 0 $(i/2*X)", basis="sto3g")
        m = pyscf.scf.RHF(mol)
        ene1 = m.kernel()
        append!(enes1, ene1)
        
        
        #=
        m2 = Molecule([atom("O",[0.0,0.0,-i/2]),atom("O",[0.0,0.0,i/2])],"ccpvdz")
        kernel2 = hf(m2)
        ene2 = run!(kernel2)
        =#
        mol = pyscf.gto.Mole()
        mol.build(atom="O 0 0 $(-i/2*X); O 0 0 $(i/2*X)", basis="ccpvdz")
        m = pyscf.scf.RHF(mol)
        ene2 = m.kernel()
        
        append!(enes2, ene2)

        append!(dis, i)
    end
    
    labels = ["sto3g" "ccpVDZ"]

    display(plot(dis, [enes1 enes2], label = labels, dpi=500))
    =#
    #savefig("O2.png")
    #display(plot(dis, ex))
    
   
    
    
    #plot(enes)
    ##@show enes
    #m = Molecule([atom("H",[0.0,0.0,-0.7]),atom("H",[0.0,0.0,0.7])],"sto3g")
    #kernel = hf(m)
    #run!(kernel)
        
end