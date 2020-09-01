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


m_H2 = Molecule([atom("H",[0,0,-0.7]), atom("H",[0,0,+0.7])], "sto3g")
m_HLi = Molecule([atom("H",[0,0,-0.7]), atom("Li",[0,0,+0.7])], "sto3g")
m_KH = Molecule([atom("K",[0,0,0.0]), atom("H",[0,0,1.0])], "sto3g")
m_ScH2 = Molecule([atom("Sc",[x1,y1,z1]), atom("H",[x2,y2,z2]), atom("H",[x3,y3,z3])], "sto3g")

@testset "Kin" begin
    @testset "get_kin1" begin
        function test_get_kin1()
            T = mol_HLi.intor("int1e_kin")
            T1 = Int1e_kin.get_kin(m_HLi)
            @test isapprox(T1, T;rtol = 1e-4, atol = 1e-15)
        end
        test_get_kin1()
    end
    @testset "get_kin2" begin
        function test_get_kin2()
            T = mol_KH.intor("int1e_kin")
            T1 = Int1e_kin.get_kin(m_KH)
            @test isapprox(T1, T;rtol = 1e-4, atol = 1e-15)
        end
        test_get_kin2()
    end
    @testset "get_kin3" begin
        function test_get_kin3()
            T = mol_ScH2.intor("int1e_kin")
            T1 = Int1e_kin.get_kin(m_ScH2)
            @test isapprox(T1, T;rtol = 1e-4, atol = 1e-15)
        end
        test_get_kin3()
    end
end