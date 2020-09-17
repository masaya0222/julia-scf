mutable struct prim_data 
    F::Vector{Float64}
    U::Matrix{Float64}
    twozeta_a::Float64
    twozeta_b::Float64
    twozeta_c::Float64
    twozeta_d::Float64
    oo2z::Float64
    oo2n::Float64
    oo2zn::Float64
    poz::Float64
    pon::Float64
    oo2p::Float64
    ss_r12_ss::Float64
    prim_data() = new()
end

mutable struct  Libint_t
    int_stack::Vector{Float64}
    PrimQuartet::Vector{prim_data}
    AB::Vector{Float64}
    CD::Vector{Float64}
    vrr_classes::Matrix{Int}
    vrr_stack::Int
    Libint_t() = new()
end
