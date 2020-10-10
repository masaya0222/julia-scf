using LinearAlgebra
function hello()
    println("hello from test1")
    println(@__FILE__)
end

mutable struct car
    name::String
    price::Int
    car() = new()
end
a = car()
a.name = "vox"

a = Vector{Vector{Float64}}()
b = [0.0,1.0,2.0]
push!(a, b)
@show a
b[1] = 3.0
@show a

