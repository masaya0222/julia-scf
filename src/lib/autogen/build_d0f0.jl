include("libints_header.jl")

# These machine-generated functions compute a quartet of (ds|fs) integrals 

function build_d0f0!(Data::prim_data, vp::Vector{Float64}, vp_ind::Int, I0::Vector{Float64}, I0_ind::Int, I1::Vector{Float64}, I1_ind::Int, I2::Vector{Float64}, I2_ind::Int, I3::Vector{Float64}, I3_ind::Int, I4::Vector{Float64}, I4_ind::Int)
	lpoz = Data.poz
	lpon = Data.pon
	oneo2zn = 1.0*Data.oo2zn
	twoo2zn = 2.0*Data.oo2zn
	threeo2zn = 3.0*Data.oo2zn
	oneo2z = 1.0*Data.oo2z
	U00 = Data.U[1,1]
	U01 = Data.U[1,2]
	U02 = Data.U[1,3]
	U40 = Data.U[5,1]
	U41 = Data.U[5,2]
	U42 = Data.U[5,3]


	vp[vp_ind+1] = U00*I0[I0_ind+1] + U40*I1[I1_ind+1] + (oneo2z)*(I2[I2_ind+1] - (lpoz)*I3[I3_ind+1]) + (threeo2zn)*I4[I4_ind+1]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+2] + U40*I1[I1_ind+2] + (oneo2z)*(I2[I2_ind+2] - (lpoz)*I3[I3_ind+2]) + (twoo2zn)*I4[I4_ind+2]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+3] + U40*I1[I1_ind+3] + (oneo2z)*(I2[I2_ind+3] - (lpoz)*I3[I3_ind+3]) + (twoo2zn)*I4[I4_ind+3]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+4] + U40*I1[I1_ind+4] + (oneo2z)*(I2[I2_ind+4] - (lpoz)*I3[I3_ind+4]) + (oneo2zn)*I4[I4_ind+4]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+5] + U40*I1[I1_ind+5] + (oneo2z)*(I2[I2_ind+5] - (lpoz)*I3[I3_ind+5]) + (oneo2zn)*I4[I4_ind+5]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+6] + U40*I1[I1_ind+6] + (oneo2z)*(I2[I2_ind+6] - (lpoz)*I3[I3_ind+6]) + (oneo2zn)*I4[I4_ind+6]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+7] + U40*I1[I1_ind+7] + (oneo2z)*(I2[I2_ind+7] - (lpoz)*I3[I3_ind+7])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+8] + U40*I1[I1_ind+8] + (oneo2z)*(I2[I2_ind+8] - (lpoz)*I3[I3_ind+8])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+9] + U40*I1[I1_ind+9] + (oneo2z)*(I2[I2_ind+9] - (lpoz)*I3[I3_ind+9])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+10] + U40*I1[I1_ind+10] + (oneo2z)*(I2[I2_ind+10] - (lpoz)*I3[I3_ind+10])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+11] + U40*I1[I1_ind+11] + (threeo2zn)*I4[I4_ind+7]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+12] + U40*I1[I1_ind+12] + (twoo2zn)*I4[I4_ind+8]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+13] + U40*I1[I1_ind+13] + (twoo2zn)*I4[I4_ind+9]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+14] + U40*I1[I1_ind+14] + (oneo2zn)*I4[I4_ind+10]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+15] + U40*I1[I1_ind+15] + (oneo2zn)*I4[I4_ind+11]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+16] + U40*I1[I1_ind+16] + (oneo2zn)*I4[I4_ind+12]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+17] + U40*I1[I1_ind+17]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+18] + U40*I1[I1_ind+18]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+19] + U40*I1[I1_ind+19]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+20] + U40*I1[I1_ind+20]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+21] + U40*I1[I1_ind+21] + (threeo2zn)*I4[I4_ind+13]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+22] + U40*I1[I1_ind+22] + (twoo2zn)*I4[I4_ind+14]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+23] + U40*I1[I1_ind+23] + (twoo2zn)*I4[I4_ind+15]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+24] + U40*I1[I1_ind+24] + (oneo2zn)*I4[I4_ind+16]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+25] + U40*I1[I1_ind+25] + (oneo2zn)*I4[I4_ind+17]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+26] + U40*I1[I1_ind+26] + (oneo2zn)*I4[I4_ind+18]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+27] + U40*I1[I1_ind+27]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+28] + U40*I1[I1_ind+28]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+29] + U40*I1[I1_ind+29]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+30] + U40*I1[I1_ind+30]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+11] + U41*I1[I1_ind+11] + (oneo2z)*(I2[I2_ind+1] - (lpoz)*I3[I3_ind+1])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+12] + U41*I1[I1_ind+12] + (oneo2z)*(I2[I2_ind+2] - (lpoz)*I3[I3_ind+2]) + (oneo2zn)*I4[I4_ind+7]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+13] + U41*I1[I1_ind+13] + (oneo2z)*(I2[I2_ind+3] - (lpoz)*I3[I3_ind+3])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+14] + U41*I1[I1_ind+14] + (oneo2z)*(I2[I2_ind+4] - (lpoz)*I3[I3_ind+4]) + (twoo2zn)*I4[I4_ind+8]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+15] + U41*I1[I1_ind+15] + (oneo2z)*(I2[I2_ind+5] - (lpoz)*I3[I3_ind+5]) + (oneo2zn)*I4[I4_ind+9]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+16] + U41*I1[I1_ind+16] + (oneo2z)*(I2[I2_ind+6] - (lpoz)*I3[I3_ind+6])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+17] + U41*I1[I1_ind+17] + (oneo2z)*(I2[I2_ind+7] - (lpoz)*I3[I3_ind+7]) + (threeo2zn)*I4[I4_ind+10]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+18] + U41*I1[I1_ind+18] + (oneo2z)*(I2[I2_ind+8] - (lpoz)*I3[I3_ind+8]) + (twoo2zn)*I4[I4_ind+11]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+19] + U41*I1[I1_ind+19] + (oneo2z)*(I2[I2_ind+9] - (lpoz)*I3[I3_ind+9]) + (oneo2zn)*I4[I4_ind+12]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+20] + U41*I1[I1_ind+20] + (oneo2z)*(I2[I2_ind+10] - (lpoz)*I3[I3_ind+10])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+21] + U41*I1[I1_ind+21]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+22] + U41*I1[I1_ind+22] + (oneo2zn)*I4[I4_ind+13]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+23] + U41*I1[I1_ind+23]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+24] + U41*I1[I1_ind+24] + (twoo2zn)*I4[I4_ind+14]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+25] + U41*I1[I1_ind+25] + (oneo2zn)*I4[I4_ind+15]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+26] + U41*I1[I1_ind+26]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+27] + U41*I1[I1_ind+27] + (threeo2zn)*I4[I4_ind+16]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+28] + U41*I1[I1_ind+28] + (twoo2zn)*I4[I4_ind+17]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+29] + U41*I1[I1_ind+29] + (oneo2zn)*I4[I4_ind+18]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+30] + U41*I1[I1_ind+30]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+21] + U42*I1[I1_ind+21] + (oneo2z)*(I2[I2_ind+1] - (lpoz)*I3[I3_ind+1])
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+22] + U42*I1[I1_ind+22] + (oneo2z)*(I2[I2_ind+2] - (lpoz)*I3[I3_ind+2])
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+23] + U42*I1[I1_ind+23] + (oneo2z)*(I2[I2_ind+3] - (lpoz)*I3[I3_ind+3]) + (oneo2zn)*I4[I4_ind+13]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+24] + U42*I1[I1_ind+24] + (oneo2z)*(I2[I2_ind+4] - (lpoz)*I3[I3_ind+4])
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+25] + U42*I1[I1_ind+25] + (oneo2z)*(I2[I2_ind+5] - (lpoz)*I3[I3_ind+5]) + (oneo2zn)*I4[I4_ind+14]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+26] + U42*I1[I1_ind+26] + (oneo2z)*(I2[I2_ind+6] - (lpoz)*I3[I3_ind+6]) + (twoo2zn)*I4[I4_ind+15]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+27] + U42*I1[I1_ind+27] + (oneo2z)*(I2[I2_ind+7] - (lpoz)*I3[I3_ind+7])
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+28] + U42*I1[I1_ind+28] + (oneo2z)*(I2[I2_ind+8] - (lpoz)*I3[I3_ind+8]) + (oneo2zn)*I4[I4_ind+16]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+29] + U42*I1[I1_ind+29] + (oneo2z)*(I2[I2_ind+9] - (lpoz)*I3[I3_ind+9]) + (twoo2zn)*I4[I4_ind+17]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+30] + U42*I1[I1_ind+30] + (oneo2z)*(I2[I2_ind+10] - (lpoz)*I3[I3_ind+10]) + (threeo2zn)*I4[I4_ind+18]
	vp_ind += 1
end
# Total number of FLOPs = 372
