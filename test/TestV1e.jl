using Test
using PyCall
using JuliaSCF.Mole
using JuliaSCF.Integral

@pyimport pyscf

X = 0.52918
mol_H2 = pyscf.gto.Mole()
mol_H2.build(atom="H 0 0 $(-0.7*X); H 0 0 $(0.7*X) ", basis="sto3g")

mol_Li2 = pyscf.gto.Mole()
mol_Li2.build(atom="Li 0 0 $(-0.7*X); Li 0 0 $(0.7*X) ", basis="sto3g")

(x1, y1, z1) = (0.4, -0.1, 1.0)
(x2, y2, z2) = (-0.3, 1.2, 1.1)

mol_IH = pyscf.gto.Mole()
mol_IH.build(atom="I $(x1*X) $(y1*X) $(z1*X);H $(x2*X) $(y2*X) $(z2*X)", basis="sto3g")


m_H2 = Molecule([atom("H",[0,0,-0.7]), atom("H",[0,0,+0.7])], "sto3g")
m_Li2 = Molecule([atom("Li",[0,0,-0.7]), atom("Li",[0,0,+0.7])], "sto3g")
m_IH = Molecule([atom("I",[x1,y1,z1]), atom("H",[x2,y2,z2])], "sto3g")

@testset "Int1e_nuc" begin
    @testset "get_v1e1" begin
        function test_get_v1e1()
            V = mol_H2.intor("int1e_nuc")
            V1 = Int1e_nuc.get_v1e(m_H2)
            @test isapprox(V1, V;rtol = 1e-5, atol = 1e-15)
        end
        test_get_v1e1()
    end
    @testset "get_v1e2" begin
        function test_get_v1e2()
            V = mol_Li2.intor("int1e_nuc")
            V1 = Int1e_nuc.get_v1e(m_Li2)
            @test isapprox(V1, V;rtol = 1e-4, atol = 1e-15)
        end
        test_get_v1e2()
    end
    @testset "get_v1e3" begin
        function test_get_v1e3()
            V = mol_IH.intor("int1e_nuc")
            V1 = Int1e_nuc.get_v1e(m_IH)
            @test isapprox(V1, V;rtol = 1e-4, atol = 1e-15)
        end
        test_get_v1e3()
    end
end
