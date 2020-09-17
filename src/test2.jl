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
