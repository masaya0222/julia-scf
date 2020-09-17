include("libints_header.jl")

# These machine-generated functions compute a quartet of (ps|gs) integrals 

function build_p0g0!(Data::prim_data, vp::Vector{Float64}, vp_ind::Int, I0::Vector{Float64}, I0_ind::Int, I1::Vector{Float64}, I1_ind::Int, I2::Vector{Float64}, I2_ind::Int, I3::Vector{Float64}, I3_ind::Int, I4::Vector{Float64}, I4_ind::Int)
	lpoz = Data.poz
	lpon = Data.pon
	oneo2zn = 1.0*Data.oo2zn
	twoo2zn = 2.0*Data.oo2zn
	threeo2zn = 3.0*Data.oo2zn
	fouro2zn = 4.0*Data.oo2zn
	U00 = Data.U[1,1]
	U01 = Data.U[1,2]
	U02 = Data.U[1,3]
	U40 = Data.U[5,1]
	U41 = Data.U[5,2]
	U42 = Data.U[5,3]


	vp[vp_ind+1] = U00*I0[I0_ind+1] + U40*I1[I1_ind+1] + (fouro2zn)*I4[I4_ind+1]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+2] + U40*I1[I1_ind+2] + (threeo2zn)*I4[I4_ind+2]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+3] + U40*I1[I1_ind+3] + (threeo2zn)*I4[I4_ind+3]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+4] + U40*I1[I1_ind+4] + (twoo2zn)*I4[I4_ind+4]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+5] + U40*I1[I1_ind+5] + (twoo2zn)*I4[I4_ind+5]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+6] + U40*I1[I1_ind+6] + (twoo2zn)*I4[I4_ind+6]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+7] + U40*I1[I1_ind+7] + (oneo2zn)*I4[I4_ind+7]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+8] + U40*I1[I1_ind+8] + (oneo2zn)*I4[I4_ind+8]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+9] + U40*I1[I1_ind+9] + (oneo2zn)*I4[I4_ind+9]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+10] + U40*I1[I1_ind+10] + (oneo2zn)*I4[I4_ind+10]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+11] + U40*I1[I1_ind+11]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+12] + U40*I1[I1_ind+12]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+13] + U40*I1[I1_ind+13]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+14] + U40*I1[I1_ind+14]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+15] + U40*I1[I1_ind+15]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+1] + U41*I1[I1_ind+1]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+2] + U41*I1[I1_ind+2] + (oneo2zn)*I4[I4_ind+1]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+3] + U41*I1[I1_ind+3]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+4] + U41*I1[I1_ind+4] + (twoo2zn)*I4[I4_ind+2]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+5] + U41*I1[I1_ind+5] + (oneo2zn)*I4[I4_ind+3]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+6] + U41*I1[I1_ind+6]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+7] + U41*I1[I1_ind+7] + (threeo2zn)*I4[I4_ind+4]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+8] + U41*I1[I1_ind+8] + (twoo2zn)*I4[I4_ind+5]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+9] + U41*I1[I1_ind+9] + (oneo2zn)*I4[I4_ind+6]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+10] + U41*I1[I1_ind+10]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+11] + U41*I1[I1_ind+11] + (fouro2zn)*I4[I4_ind+7]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+12] + U41*I1[I1_ind+12] + (threeo2zn)*I4[I4_ind+8]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+13] + U41*I1[I1_ind+13] + (twoo2zn)*I4[I4_ind+9]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+14] + U41*I1[I1_ind+14] + (oneo2zn)*I4[I4_ind+10]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+15] + U41*I1[I1_ind+15]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+1] + U42*I1[I1_ind+1]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+2] + U42*I1[I1_ind+2]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+3] + U42*I1[I1_ind+3] + (oneo2zn)*I4[I4_ind+1]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+4] + U42*I1[I1_ind+4]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+5] + U42*I1[I1_ind+5] + (oneo2zn)*I4[I4_ind+2]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+6] + U42*I1[I1_ind+6] + (twoo2zn)*I4[I4_ind+3]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+7] + U42*I1[I1_ind+7]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+8] + U42*I1[I1_ind+8] + (oneo2zn)*I4[I4_ind+4]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+9] + U42*I1[I1_ind+9] + (twoo2zn)*I4[I4_ind+5]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+10] + U42*I1[I1_ind+10] + (threeo2zn)*I4[I4_ind+6]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+11] + U42*I1[I1_ind+11]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+12] + U42*I1[I1_ind+12] + (oneo2zn)*I4[I4_ind+7]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+13] + U42*I1[I1_ind+13] + (twoo2zn)*I4[I4_ind+8]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+14] + U42*I1[I1_ind+14] + (threeo2zn)*I4[I4_ind+9]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+15] + U42*I1[I1_ind+15] + (fouro2zn)*I4[I4_ind+10]
	vp_ind += 1
end
# Total number of FLOPs = 195