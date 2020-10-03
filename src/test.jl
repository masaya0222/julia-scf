module test
    using JuliaSCF.Mole
    using JuliaSCF.Integral
    using SpecialFunctions
    using LinearAlgebra
    using JuliaSCF.HF
    using JuliaSCF.UHF
    using Plots
    using LaTeXStrings
    using Lints
    X = 0.52918
    #=
    @lints begin
        mol = Lints.Molecule([1,1], [[0.0,0.0,-1.4*X],
                                        [0.0,0.0,1.4*X]])
        bas = Lints.BasisSet("sto-3g", mol)
        S = Lints.make_S(bas)
        @show size(S)
        @show ERI = Lints.make_ERI4(bas)
        @show size(ERI)
    end
    =#
    #m = Molecule([atom("C",[0.0,0.0,-0.7]), atom("C",[0.0,0.0,0.7]), atom("C",[0.0,0.0,1.4]), atom("C",[0.0,1.4,0.0]), atom("C",[1.4,0.0,0.0])], "sto3g")
    #m = Molecule([atom("C",[0.0,0.0,-0.7]),atom("C",[0.0,0.0,0.7])],"6-31g")
    #m = Molecule([atom("O",[0.0,0.0,0.0]),atom("H",[1.0,0.0,0.0]), atom("H",[0.0,1.0,0.0])],"ccpvdz")
    #@time kernel1 = HF.hf(m)
    #@time ene1 = run!(kernel1)
        
    #@time v = Int2e.get_v2e_single(m)
    #@time v = Int2e.get_v2e(m)
    #@show isapprox(v,ERI)
    
    using PyCall
    @pyimport pyscf
    X = 0.52918
    mol = pyscf.gto.Mole()
    
    mol.build(atom="H 0 0 $(-1.4*X); H 0 0 $(1.4*X)", basis="sto-3g")
    
    #mol.build(atom="C 0 0 $(-0.7*X); C 0 0 $(0.4*X); C 0 $(0.4*X) 0; C 0 $(-0.4*X) 0; C $(-0.4*X) 0 0", basis="sto3g")
    #@show  v = mol.intor("int2e")
    #size(v)
    
    gr()
    
    
    enes1 = []
    enes2 = []
    dis = []
    
    #H_ene = -0.466581849557275 #sto3g
    #H_ene = -0.49823291072907 #6-31g
    for i in 1.5:0.1:4
    #for i = 3.0
        @show i
        
        m1 = Molecule([atom("Li",[0.0,0.0,-i/2]),atom("Li",[0.0,0.0,i/2])],"6-31g")
        kernel1 = hf(m1)
        ene1 = run_rhf!(kernel1)
        #v = Int2e.get_v2e(m1)
        #@show v[1,1,1,:]
        append!(enes1, ene1)
        
        
        mol = pyscf.gto.Mole()
        mol.build(atom="Li 0 0 $(-i*X/2); Li 0 0 $(i*X/2)", basis="6-31g")
        m = pyscf.scf.RHF(mol)
        ene2 = m.kernel()
        append!(enes2, ene2)
        #v1 = mol.intor("int2e")
        #@show v1[1,1,1,:]
        append!(dis, i)
    end
    
    label = ["HF (6-31g)" "UHF (6-31g)"]
    xlabel=L"R (a.u.)"
    ylabel=L"E(H_2)-2E(H) (a.u.)"
    
    display(plot(dis, [enes1 enes2], label = label))#, xlabel=xlabel, ylabel=ylabel))#, dpi=500))
    
    #savefig("H2_6-31g.png")
    
    #display(plot(dis, ex))
    
   
    
    
    #plot(enes)
    ##@show enes
    #m = Molecule([atom("H",[0.0,0.0,-0.7]),atom("H",[0.0,0.0,0.7])],"sto3g")
    #kernel = hf(m)
    #run!(kernel)
     
end