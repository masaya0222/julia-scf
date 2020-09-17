include("libints_header.jl")

# These machine-generated functions compute a quartet of (fs|ps) integrals 

function build_f0p0!(Data::prim_data, vp::Vector{Float64}, vp_ind::Int, I0::Vector{Float64}, I0_ind::Int, I1::Vector{Float64}, I1_ind::Int, I2::Vector{Float64}, I2_ind::Int, I3::Vector{Float64}, I3_ind::Int, I4::Vector{Float64}, I4_ind::Int)
	lpoz = Data.poz
	lpon = Data.pon
	oneo2zn = 1.0*Data.oo2zn
	oneo2z = 1.0*Data.oo2z
	twoo2z = 2.0*Data.oo2z
	U00 = Data.U[1,1]
	U01 = Data.U[1,2]
	U02 = Data.U[1,3]
	U40 = Data.U[5,1]
	U41 = Data.U[5,2]
	U42 = Data.U[5,3]


	vp[vp_ind+1] = U00*I0[I0_ind+1] + U40*I1[I1_ind+1] + (twoo2z)*(I2[I2_ind+1] - (lpoz)*I3[I3_ind+1]) + (oneo2zn)*I4[I4_ind+1]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+2] + U40*I1[I1_ind+2] + (twoo2z)*(I2[I2_ind+2] - (lpoz)*I3[I3_ind+2])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+3] + U40*I1[I1_ind+3] + (twoo2z)*(I2[I2_ind+3] - (lpoz)*I3[I3_ind+3])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+4] + U40*I1[I1_ind+4] + (oneo2z)*(I2[I2_ind+4] - (lpoz)*I3[I3_ind+4]) + (oneo2zn)*I4[I4_ind+2]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+5] + U40*I1[I1_ind+5] + (oneo2z)*(I2[I2_ind+5] - (lpoz)*I3[I3_ind+5])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+6] + U40*I1[I1_ind+6] + (oneo2z)*(I2[I2_ind+6] - (lpoz)*I3[I3_ind+6])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+7] + U40*I1[I1_ind+7] + (oneo2z)*(I2[I2_ind+7] - (lpoz)*I3[I3_ind+7]) + (oneo2zn)*I4[I4_ind+3]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+8] + U40*I1[I1_ind+8] + (oneo2z)*(I2[I2_ind+8] - (lpoz)*I3[I3_ind+8])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+9] + U40*I1[I1_ind+9] + (oneo2z)*(I2[I2_ind+9] - (lpoz)*I3[I3_ind+9])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+10] + U40*I1[I1_ind+10] + (oneo2zn)*I4[I4_ind+4]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+11] + U40*I1[I1_ind+11]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+12] + U40*I1[I1_ind+12]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+13] + U40*I1[I1_ind+13] + (oneo2zn)*I4[I4_ind+5]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+14] + U40*I1[I1_ind+14]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+15] + U40*I1[I1_ind+15]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+16] + U40*I1[I1_ind+16] + (oneo2zn)*I4[I4_ind+6]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+17] + U40*I1[I1_ind+17]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+18] + U40*I1[I1_ind+18]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+10] + U41*I1[I1_ind+10] + (twoo2z)*(I2[I2_ind+4] - (lpoz)*I3[I3_ind+4])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+11] + U41*I1[I1_ind+11] + (twoo2z)*(I2[I2_ind+5] - (lpoz)*I3[I3_ind+5]) + (oneo2zn)*I4[I4_ind+4]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+12] + U41*I1[I1_ind+12] + (twoo2z)*(I2[I2_ind+6] - (lpoz)*I3[I3_ind+6])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+13] + U41*I1[I1_ind+13] + (oneo2z)*(I2[I2_ind+7] - (lpoz)*I3[I3_ind+7])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+14] + U41*I1[I1_ind+14] + (oneo2z)*(I2[I2_ind+8] - (lpoz)*I3[I3_ind+8]) + (oneo2zn)*I4[I4_ind+5]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+15] + U41*I1[I1_ind+15] + (oneo2z)*(I2[I2_ind+9] - (lpoz)*I3[I3_ind+9])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+16] + U41*I1[I1_ind+16]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+17] + U41*I1[I1_ind+17] + (oneo2zn)*I4[I4_ind+6]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+18] + U41*I1[I1_ind+18]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+16] + U42*I1[I1_ind+16] + (twoo2z)*(I2[I2_ind+7] - (lpoz)*I3[I3_ind+7])
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+17] + U42*I1[I1_ind+17] + (twoo2z)*(I2[I2_ind+8] - (lpoz)*I3[I3_ind+8])
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+18] + U42*I1[I1_ind+18] + (twoo2z)*(I2[I2_ind+9] - (lpoz)*I3[I3_ind+9]) + (oneo2zn)*I4[I4_ind+6]
	vp_ind += 1
end
# Total number of FLOPs = 182
