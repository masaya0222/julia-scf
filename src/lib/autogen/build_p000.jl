include("libints_header.jl")

# These machine-generated functions compute a quartet of (ps|0s) integrals 

function build_p000!(Data::prim_data, vp::Vector{Float64}, vp_ind::Int, I0::Vector{Float64}, I0_ind::Int, I1::Vector{Float64}, I1_ind::Int, I2::Vector{Float64}, I2_ind::Int, I3::Vector{Float64}, I3_ind::Int, I4::Vector{Float64}, I4_ind::Int)
	lpoz = Data.poz
	lpon = Data.pon
	U00 = Data.U[1,1]
	U01 = Data.U[1,2]
	U02 = Data.U[1,3]
	U40 = Data.U[5,1]
	U41 = Data.U[5,2]
	U42 = Data.U[5,3]


	vp[vp_ind+1] = U00*I0[I0_ind+1] + U40*I1[I1_ind+1]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+1] + U41*I1[I1_ind+1]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+1] + U42*I1[I1_ind+1]
	vp_ind += 1
end
# Total number of FLOPs = 9
