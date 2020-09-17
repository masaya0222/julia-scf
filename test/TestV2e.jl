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

mol_I2 = pyscf.gto.Mole()
mol_I2.build(atom="I $(x1*X) $(y1*X) $(z1*X);I $(x2*X) $(y2*X) $(z2*X)", basis="sto3g")


m_H2 = Molecule([atom("H",[0,0,-0.7]), atom("H",[0,0,+0.7])], "sto3g")
m_Li2 = Molecule([atom("Li",[0,0,-0.7]), atom("Li",[0,0,+0.7])], "sto3g")
m_I2 = Molecule([atom("I",[x1,y1,z1]), atom("I",[x2,y2,z2])], "sto3g")

@testset "Int2e" begin
    @testset "get_v2e1" begin
        function test_get_v2e1()
            V = mol_H2.intor("int2e")
            V1 = Int2e.get_v2e(m_H2)
            @test isapprox(V1, V; rtol=5e-6, atol=1e-15)
        end
        test_get_v2e1()
    end
    @testset "get_v2e2" begin
        function test_get_v2e2()
            V = mol_Li2.intor("int2e")
            V1 = Int2e.get_v2e(m_Li2)
            @test isapprox(V1, V; rtol=1e-5, atol=1e-15)
        end
        test_get_v2e2()
    end
    @testset "get_v2e3" begin
        function test_get_v2e3()
            V = mol_I2.intor("int2e")
            V1 = Int2e.get_v2e(m_I2)
            
            @test isapprox(V1, V; rtol=1e-5, atol=1e-15)
        end
        test_get_v2e3()
    end
end