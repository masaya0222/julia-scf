using Test
using PyCall
using JuliaSCF.Mole
using JuliaSCF.Integral

@pyimport pyscf

X = 0.52918
mol_H2 = pyscf.gto.Mole()
mol_H2.build(atom="H 0 0 $(-0.7*X); H 0 0 $(0.7*X) ", basis="sto3g")

mol_HLi = pyscf.gto.Mole()
mol_HLi.build(atom="H 0 0 $(-0.7*X); Li 0 0 $(0.7*X) ", basis="sto3g")

mol_KH = pyscf.gto.Mole()
mol_KH.build(atom="K 0 0 $(0*X); H 0 0 $(1*X) ", basis="sto3g")
(x1, y1, z1) = (0.4, -0.1, 1.0)
(x2, y2, z2) = (-0.3, 1.2, 1.1)
(x3, y3, z3) = (-0.5, 0.6, -0.1)
mol_ScH2 = pyscf.gto.Mole()
mol_ScH2.build(atom="Sc $(x1*X) $(y1*X) $(z1*X);H $(x2*X) $(y2*X) $(z2*X);H $(x3*X) $(y3*X) $(z3*X)", basis="sto3g", charge=+1)

mol_IH = pyscf.gto.Mole()
mol_IH.build(atom="I $(x1*X) $(y1*X) $(z1*X);H $(x2*X) $(y2*X) $(z2*X)", basis="sto3g")


m_H2 = Molecule([atom("H",[0,0,-0.7]), atom("H",[0,0,+0.7])], "sto3g")
m_HLi = Molecule([atom("H",[0,0,-0.7]), atom("Li",[0,0,+0.7])], "sto3g")
m_KH = Molecule([atom("K",[0,0,0.0]), atom("H",[0,0,1.0])], "sto3g")
m_ScH2 = Molecule([atom("Sc",[x1,y1,z1]), atom("H",[x2,y2,z2]), atom("H",[x3,y3,z3])], "sto3g")
m_IH = Molecule([atom("I",[x1,y1,z1]), atom("H",[x2,y2,z2])], "sto3g")

@testset "Overlap" begin
    @testset "Sab" begin
        function test_Sab()
            basis_a = orb_detail([0, 0,-0.7], 0, [3.425250914, 0.6239137298, 0.168855404], [[0.1543289673, 0.5353281423, 0.4446345422]])
            basis_b = orb_detail([0, 0, 0.7], 0, [3.425250914, 0.6239137298, 0.168855404], [[0.1543289673, 0.5353281423, 0.4446345422]])
            S_aa = Int1e_ovlp.cont_Sij(basis_a, basis_a)
            S_ab = Int1e_ovlp.cont_Sij(basis_a, basis_b)
            S_ba = Int1e_ovlp.cont_Sij(basis_b, basis_a)
            S_bb = Int1e_ovlp.cont_Sij(basis_b, basis_b)

            @test isapprox([S_aa[1], S_ab[1], S_ba[1], S_bb[1]],[1.0, 0.6593182, 0.6593182, 1.0, ]; rtol =1e-6)
        end
        test_Sab()
    end
    @testset "Slm1" begin
        function test_Slm1()
            S1 = mol_H2.intor("int1e_ovlp")
            basis_a = orb_detail([0, 0,-0.7], 0, [3.425250914, 0.6239137298, 0.168855404], [[0.1543289673, 0.5353281423, 0.4446345422]])
            basis_b = orb_detail([0, 0, 0.7], 0, [3.425250914, 0.6239137298, 0.168855404], [[0.1543289673, 0.5353281423, 0.4446345422]])
            S_aa = Int1e_ovlp.S_lm(basis_a, basis_a)
            S_ab = Int1e_ovlp.S_lm(basis_a, basis_b)
            S_ba = Int1e_ovlp.S_lm(basis_b, basis_a)
            S_bb = Int1e_ovlp.S_lm(basis_b, basis_b)
            
            @test isapprox([S_aa[1] S_ab[1]; S_ba[1] S_bb[1]],S1; rtol =1e-5)
        end
        test_Slm1()
    end
    @testset "get_ovlp1" begin
        function test_get_ovlp1()
            S1 = Int1e_ovlp.get_ovlp(m_HLi)
            S = mol_HLi.intor("int1e_ovlp")
            @test isapprox(S1, S;rtol = 1e-4, atol = 1e-15)
        end
        test_get_ovlp1()
    end
    @testset "get_ovlp2" begin
        function test_get_ovlp2()
            S1 = Int1e_ovlp.get_ovlp(m_KH)
            S = mol_KH.intor("int1e_ovlp")
            @test isapprox(S1, S;rtol = 1e-4, atol = 1e-15)
        end
        test_get_ovlp2()
    end
    @testset "get_ovlp3" begin
        function test_get_ovlp3()
            S1 = Int1e_ovlp.get_ovlp(m_ScH2)
            S = mol_ScH2.intor("int1e_ovlp")
            
            @test isapprox(S1, S;rtol = 1e-4, atol = 1e-15)
        end
        test_get_ovlp3()
    end
    @testset "get_ovlp4" begin
        function test_get_ovlp4()
            S1 = Int1e_ovlp.get_ovlp(m_IH)
            S = mol_IH.intor("int1e_ovlp")
            @test isapprox(S1, S;rtol = 1e-5, atol = 1e-15)
        end
        test_get_ovlp4()
    end

    @testset "Multiplication2" begin
        @test 1 * 1 == 1
        @test 2 * 3 == 6
    end
end

