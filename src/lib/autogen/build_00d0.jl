include("libints_header.jl")

# These machine-generated functions compute a quartet of (0s|ds) integrals 

function build_00d0!(Data::prim_data, vp::Vector{Float64}, vp_ind::Int, I0::Vector{Float64}, I0_ind::Int, I1::Vector{Float64}, I1_ind::Int, I2::Vector{Float64}, I2_ind::Int, I3::Vector{Float64}, I3_ind::Int, I4::Vector{Float64}, I4_ind::Int)
	lpoz = Data.poz
	lpon = Data.pon
	oneo2n = 1.0*Data.oo2n
	U20 = Data.U[3,1]
	U21 = Data.U[3,2]
	U22 = Data.U[3,3]
	U50 = Data.U[6,1]
	U51 = Data.U[6,2]
	U52 = Data.U[6,3]


	vp[vp_ind+1] = U20*I0[I0_ind+1] + U50*I1[I1_ind+1] + (oneo2n)*(I2[I2_ind+1] - (lpon)*I3[I3_ind+1])
	vp_ind += 1
	vp[vp_ind+1] = U20*I0[I0_ind+2] + U50*I1[I1_ind+2]
	vp_ind += 1
	vp[vp_ind+1] = U20*I0[I0_ind+3] + U50*I1[I1_ind+3]
	vp_ind += 1
	vp[vp_ind+1] = U21*I0[I0_ind+2] + U51*I1[I1_ind+2] + (oneo2n)*(I2[I2_ind+1] - (lpon)*I3[I3_ind+1])
	vp_ind += 1
	vp[vp_ind+1] = U21*I0[I0_ind+3] + U51*I1[I1_ind+3]
	vp_ind += 1
	vp[vp_ind+1] = U22*I0[I0_ind+3] + U52*I1[I1_ind+3] + (oneo2n)*(I2[I2_ind+1] - (lpon)*I3[I3_ind+1])
	vp_ind += 1
end
# Total number of FLOPs = 30
