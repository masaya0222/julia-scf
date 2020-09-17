include("libints_header.jl")

# These machine-generated functions compute a quartet of (gs|gs) integrals 

function build_g0g0!(Data::prim_data, vp::Vector{Float64}, vp_ind::Int, I0::Vector{Float64}, I0_ind::Int, I1::Vector{Float64}, I1_ind::Int, I2::Vector{Float64}, I2_ind::Int, I3::Vector{Float64}, I3_ind::Int, I4::Vector{Float64}, I4_ind::Int)
	lpoz = Data.poz
	lpon = Data.pon
	oneo2zn = 1.0*Data.oo2zn
	twoo2zn = 2.0*Data.oo2zn
	threeo2zn = 3.0*Data.oo2zn
	fouro2zn = 4.0*Data.oo2zn
	oneo2z = 1.0*Data.oo2z
	twoo2z = 2.0*Data.oo2z
	threeo2z = 3.0*Data.oo2z
	U00 = Data.U[1,1]
	U01 = Data.U[1,2]
	U02 = Data.U[1,3]
	U40 = Data.U[5,1]
	U41 = Data.U[5,2]
	U42 = Data.U[5,3]


	vp[vp_ind+1] = U00*I0[I0_ind+1] + U40*I1[I1_ind+1] + (threeo2z)*(I2[I2_ind+1] - (lpoz)*I3[I3_ind+1]) + (fouro2zn)*I4[I4_ind+1]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+2] + U40*I1[I1_ind+2] + (threeo2z)*(I2[I2_ind+2] - (lpoz)*I3[I3_ind+2]) + (threeo2zn)*I4[I4_ind+2]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+3] + U40*I1[I1_ind+3] + (threeo2z)*(I2[I2_ind+3] - (lpoz)*I3[I3_ind+3]) + (threeo2zn)*I4[I4_ind+3]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+4] + U40*I1[I1_ind+4] + (threeo2z)*(I2[I2_ind+4] - (lpoz)*I3[I3_ind+4]) + (twoo2zn)*I4[I4_ind+4]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+5] + U40*I1[I1_ind+5] + (threeo2z)*(I2[I2_ind+5] - (lpoz)*I3[I3_ind+5]) + (twoo2zn)*I4[I4_ind+5]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+6] + U40*I1[I1_ind+6] + (threeo2z)*(I2[I2_ind+6] - (lpoz)*I3[I3_ind+6]) + (twoo2zn)*I4[I4_ind+6]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+7] + U40*I1[I1_ind+7] + (threeo2z)*(I2[I2_ind+7] - (lpoz)*I3[I3_ind+7]) + (oneo2zn)*I4[I4_ind+7]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+8] + U40*I1[I1_ind+8] + (threeo2z)*(I2[I2_ind+8] - (lpoz)*I3[I3_ind+8]) + (oneo2zn)*I4[I4_ind+8]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+9] + U40*I1[I1_ind+9] + (threeo2z)*(I2[I2_ind+9] - (lpoz)*I3[I3_ind+9]) + (oneo2zn)*I4[I4_ind+9]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+10] + U40*I1[I1_ind+10] + (threeo2z)*(I2[I2_ind+10] - (lpoz)*I3[I3_ind+10]) + (oneo2zn)*I4[I4_ind+10]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+11] + U40*I1[I1_ind+11] + (threeo2z)*(I2[I2_ind+11] - (lpoz)*I3[I3_ind+11])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+12] + U40*I1[I1_ind+12] + (threeo2z)*(I2[I2_ind+12] - (lpoz)*I3[I3_ind+12])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+13] + U40*I1[I1_ind+13] + (threeo2z)*(I2[I2_ind+13] - (lpoz)*I3[I3_ind+13])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+14] + U40*I1[I1_ind+14] + (threeo2z)*(I2[I2_ind+14] - (lpoz)*I3[I3_ind+14])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+15] + U40*I1[I1_ind+15] + (threeo2z)*(I2[I2_ind+15] - (lpoz)*I3[I3_ind+15])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+16] + U40*I1[I1_ind+16] + (twoo2z)*(I2[I2_ind+16] - (lpoz)*I3[I3_ind+16]) + (fouro2zn)*I4[I4_ind+11]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+17] + U40*I1[I1_ind+17] + (twoo2z)*(I2[I2_ind+17] - (lpoz)*I3[I3_ind+17]) + (threeo2zn)*I4[I4_ind+12]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+18] + U40*I1[I1_ind+18] + (twoo2z)*(I2[I2_ind+18] - (lpoz)*I3[I3_ind+18]) + (threeo2zn)*I4[I4_ind+13]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+19] + U40*I1[I1_ind+19] + (twoo2z)*(I2[I2_ind+19] - (lpoz)*I3[I3_ind+19]) + (twoo2zn)*I4[I4_ind+14]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+20] + U40*I1[I1_ind+20] + (twoo2z)*(I2[I2_ind+20] - (lpoz)*I3[I3_ind+20]) + (twoo2zn)*I4[I4_ind+15]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+21] + U40*I1[I1_ind+21] + (twoo2z)*(I2[I2_ind+21] - (lpoz)*I3[I3_ind+21]) + (twoo2zn)*I4[I4_ind+16]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+22] + U40*I1[I1_ind+22] + (twoo2z)*(I2[I2_ind+22] - (lpoz)*I3[I3_ind+22]) + (oneo2zn)*I4[I4_ind+17]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+23] + U40*I1[I1_ind+23] + (twoo2z)*(I2[I2_ind+23] - (lpoz)*I3[I3_ind+23]) + (oneo2zn)*I4[I4_ind+18]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+24] + U40*I1[I1_ind+24] + (twoo2z)*(I2[I2_ind+24] - (lpoz)*I3[I3_ind+24]) + (oneo2zn)*I4[I4_ind+19]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+25] + U40*I1[I1_ind+25] + (twoo2z)*(I2[I2_ind+25] - (lpoz)*I3[I3_ind+25]) + (oneo2zn)*I4[I4_ind+20]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+26] + U40*I1[I1_ind+26] + (twoo2z)*(I2[I2_ind+26] - (lpoz)*I3[I3_ind+26])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+27] + U40*I1[I1_ind+27] + (twoo2z)*(I2[I2_ind+27] - (lpoz)*I3[I3_ind+27])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+28] + U40*I1[I1_ind+28] + (twoo2z)*(I2[I2_ind+28] - (lpoz)*I3[I3_ind+28])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+29] + U40*I1[I1_ind+29] + (twoo2z)*(I2[I2_ind+29] - (lpoz)*I3[I3_ind+29])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+30] + U40*I1[I1_ind+30] + (twoo2z)*(I2[I2_ind+30] - (lpoz)*I3[I3_ind+30])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+31] + U40*I1[I1_ind+31] + (twoo2z)*(I2[I2_ind+31] - (lpoz)*I3[I3_ind+31]) + (fouro2zn)*I4[I4_ind+21]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+32] + U40*I1[I1_ind+32] + (twoo2z)*(I2[I2_ind+32] - (lpoz)*I3[I3_ind+32]) + (threeo2zn)*I4[I4_ind+22]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+33] + U40*I1[I1_ind+33] + (twoo2z)*(I2[I2_ind+33] - (lpoz)*I3[I3_ind+33]) + (threeo2zn)*I4[I4_ind+23]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+34] + U40*I1[I1_ind+34] + (twoo2z)*(I2[I2_ind+34] - (lpoz)*I3[I3_ind+34]) + (twoo2zn)*I4[I4_ind+24]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+35] + U40*I1[I1_ind+35] + (twoo2z)*(I2[I2_ind+35] - (lpoz)*I3[I3_ind+35]) + (twoo2zn)*I4[I4_ind+25]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+36] + U40*I1[I1_ind+36] + (twoo2z)*(I2[I2_ind+36] - (lpoz)*I3[I3_ind+36]) + (twoo2zn)*I4[I4_ind+26]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+37] + U40*I1[I1_ind+37] + (twoo2z)*(I2[I2_ind+37] - (lpoz)*I3[I3_ind+37]) + (oneo2zn)*I4[I4_ind+27]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+38] + U40*I1[I1_ind+38] + (twoo2z)*(I2[I2_ind+38] - (lpoz)*I3[I3_ind+38]) + (oneo2zn)*I4[I4_ind+28]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+39] + U40*I1[I1_ind+39] + (twoo2z)*(I2[I2_ind+39] - (lpoz)*I3[I3_ind+39]) + (oneo2zn)*I4[I4_ind+29]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+40] + U40*I1[I1_ind+40] + (twoo2z)*(I2[I2_ind+40] - (lpoz)*I3[I3_ind+40]) + (oneo2zn)*I4[I4_ind+30]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+41] + U40*I1[I1_ind+41] + (twoo2z)*(I2[I2_ind+41] - (lpoz)*I3[I3_ind+41])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+42] + U40*I1[I1_ind+42] + (twoo2z)*(I2[I2_ind+42] - (lpoz)*I3[I3_ind+42])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+43] + U40*I1[I1_ind+43] + (twoo2z)*(I2[I2_ind+43] - (lpoz)*I3[I3_ind+43])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+44] + U40*I1[I1_ind+44] + (twoo2z)*(I2[I2_ind+44] - (lpoz)*I3[I3_ind+44])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+45] + U40*I1[I1_ind+45] + (twoo2z)*(I2[I2_ind+45] - (lpoz)*I3[I3_ind+45])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+46] + U40*I1[I1_ind+46] + (oneo2z)*(I2[I2_ind+46] - (lpoz)*I3[I3_ind+46]) + (fouro2zn)*I4[I4_ind+31]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+47] + U40*I1[I1_ind+47] + (oneo2z)*(I2[I2_ind+47] - (lpoz)*I3[I3_ind+47]) + (threeo2zn)*I4[I4_ind+32]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+48] + U40*I1[I1_ind+48] + (oneo2z)*(I2[I2_ind+48] - (lpoz)*I3[I3_ind+48]) + (threeo2zn)*I4[I4_ind+33]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+49] + U40*I1[I1_ind+49] + (oneo2z)*(I2[I2_ind+49] - (lpoz)*I3[I3_ind+49]) + (twoo2zn)*I4[I4_ind+34]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+50] + U40*I1[I1_ind+50] + (oneo2z)*(I2[I2_ind+50] - (lpoz)*I3[I3_ind+50]) + (twoo2zn)*I4[I4_ind+35]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+51] + U40*I1[I1_ind+51] + (oneo2z)*(I2[I2_ind+51] - (lpoz)*I3[I3_ind+51]) + (twoo2zn)*I4[I4_ind+36]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+52] + U40*I1[I1_ind+52] + (oneo2z)*(I2[I2_ind+52] - (lpoz)*I3[I3_ind+52]) + (oneo2zn)*I4[I4_ind+37]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+53] + U40*I1[I1_ind+53] + (oneo2z)*(I2[I2_ind+53] - (lpoz)*I3[I3_ind+53]) + (oneo2zn)*I4[I4_ind+38]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+54] + U40*I1[I1_ind+54] + (oneo2z)*(I2[I2_ind+54] - (lpoz)*I3[I3_ind+54]) + (oneo2zn)*I4[I4_ind+39]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+55] + U40*I1[I1_ind+55] + (oneo2z)*(I2[I2_ind+55] - (lpoz)*I3[I3_ind+55]) + (oneo2zn)*I4[I4_ind+40]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+56] + U40*I1[I1_ind+56] + (oneo2z)*(I2[I2_ind+56] - (lpoz)*I3[I3_ind+56])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+57] + U40*I1[I1_ind+57] + (oneo2z)*(I2[I2_ind+57] - (lpoz)*I3[I3_ind+57])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+58] + U40*I1[I1_ind+58] + (oneo2z)*(I2[I2_ind+58] - (lpoz)*I3[I3_ind+58])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+59] + U40*I1[I1_ind+59] + (oneo2z)*(I2[I2_ind+59] - (lpoz)*I3[I3_ind+59])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+60] + U40*I1[I1_ind+60] + (oneo2z)*(I2[I2_ind+60] - (lpoz)*I3[I3_ind+60])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+61] + U40*I1[I1_ind+61] + (oneo2z)*(I2[I2_ind+61] - (lpoz)*I3[I3_ind+61]) + (fouro2zn)*I4[I4_ind+41]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+62] + U40*I1[I1_ind+62] + (oneo2z)*(I2[I2_ind+62] - (lpoz)*I3[I3_ind+62]) + (threeo2zn)*I4[I4_ind+42]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+63] + U40*I1[I1_ind+63] + (oneo2z)*(I2[I2_ind+63] - (lpoz)*I3[I3_ind+63]) + (threeo2zn)*I4[I4_ind+43]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+64] + U40*I1[I1_ind+64] + (oneo2z)*(I2[I2_ind+64] - (lpoz)*I3[I3_ind+64]) + (twoo2zn)*I4[I4_ind+44]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+65] + U40*I1[I1_ind+65] + (oneo2z)*(I2[I2_ind+65] - (lpoz)*I3[I3_ind+65]) + (twoo2zn)*I4[I4_ind+45]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+66] + U40*I1[I1_ind+66] + (oneo2z)*(I2[I2_ind+66] - (lpoz)*I3[I3_ind+66]) + (twoo2zn)*I4[I4_ind+46]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+67] + U40*I1[I1_ind+67] + (oneo2z)*(I2[I2_ind+67] - (lpoz)*I3[I3_ind+67]) + (oneo2zn)*I4[I4_ind+47]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+68] + U40*I1[I1_ind+68] + (oneo2z)*(I2[I2_ind+68] - (lpoz)*I3[I3_ind+68]) + (oneo2zn)*I4[I4_ind+48]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+69] + U40*I1[I1_ind+69] + (oneo2z)*(I2[I2_ind+69] - (lpoz)*I3[I3_ind+69]) + (oneo2zn)*I4[I4_ind+49]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+70] + U40*I1[I1_ind+70] + (oneo2z)*(I2[I2_ind+70] - (lpoz)*I3[I3_ind+70]) + (oneo2zn)*I4[I4_ind+50]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+71] + U40*I1[I1_ind+71] + (oneo2z)*(I2[I2_ind+71] - (lpoz)*I3[I3_ind+71])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+72] + U40*I1[I1_ind+72] + (oneo2z)*(I2[I2_ind+72] - (lpoz)*I3[I3_ind+72])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+73] + U40*I1[I1_ind+73] + (oneo2z)*(I2[I2_ind+73] - (lpoz)*I3[I3_ind+73])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+74] + U40*I1[I1_ind+74] + (oneo2z)*(I2[I2_ind+74] - (lpoz)*I3[I3_ind+74])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+75] + U40*I1[I1_ind+75] + (oneo2z)*(I2[I2_ind+75] - (lpoz)*I3[I3_ind+75])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+76] + U40*I1[I1_ind+76] + (oneo2z)*(I2[I2_ind+76] - (lpoz)*I3[I3_ind+76]) + (fouro2zn)*I4[I4_ind+51]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+77] + U40*I1[I1_ind+77] + (oneo2z)*(I2[I2_ind+77] - (lpoz)*I3[I3_ind+77]) + (threeo2zn)*I4[I4_ind+52]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+78] + U40*I1[I1_ind+78] + (oneo2z)*(I2[I2_ind+78] - (lpoz)*I3[I3_ind+78]) + (threeo2zn)*I4[I4_ind+53]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+79] + U40*I1[I1_ind+79] + (oneo2z)*(I2[I2_ind+79] - (lpoz)*I3[I3_ind+79]) + (twoo2zn)*I4[I4_ind+54]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+80] + U40*I1[I1_ind+80] + (oneo2z)*(I2[I2_ind+80] - (lpoz)*I3[I3_ind+80]) + (twoo2zn)*I4[I4_ind+55]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+81] + U40*I1[I1_ind+81] + (oneo2z)*(I2[I2_ind+81] - (lpoz)*I3[I3_ind+81]) + (twoo2zn)*I4[I4_ind+56]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+82] + U40*I1[I1_ind+82] + (oneo2z)*(I2[I2_ind+82] - (lpoz)*I3[I3_ind+82]) + (oneo2zn)*I4[I4_ind+57]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+83] + U40*I1[I1_ind+83] + (oneo2z)*(I2[I2_ind+83] - (lpoz)*I3[I3_ind+83]) + (oneo2zn)*I4[I4_ind+58]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+84] + U40*I1[I1_ind+84] + (oneo2z)*(I2[I2_ind+84] - (lpoz)*I3[I3_ind+84]) + (oneo2zn)*I4[I4_ind+59]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+85] + U40*I1[I1_ind+85] + (oneo2z)*(I2[I2_ind+85] - (lpoz)*I3[I3_ind+85]) + (oneo2zn)*I4[I4_ind+60]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+86] + U40*I1[I1_ind+86] + (oneo2z)*(I2[I2_ind+86] - (lpoz)*I3[I3_ind+86])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+87] + U40*I1[I1_ind+87] + (oneo2z)*(I2[I2_ind+87] - (lpoz)*I3[I3_ind+87])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+88] + U40*I1[I1_ind+88] + (oneo2z)*(I2[I2_ind+88] - (lpoz)*I3[I3_ind+88])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+89] + U40*I1[I1_ind+89] + (oneo2z)*(I2[I2_ind+89] - (lpoz)*I3[I3_ind+89])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+90] + U40*I1[I1_ind+90] + (oneo2z)*(I2[I2_ind+90] - (lpoz)*I3[I3_ind+90])
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+91] + U40*I1[I1_ind+91] + (fouro2zn)*I4[I4_ind+61]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+92] + U40*I1[I1_ind+92] + (threeo2zn)*I4[I4_ind+62]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+93] + U40*I1[I1_ind+93] + (threeo2zn)*I4[I4_ind+63]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+94] + U40*I1[I1_ind+94] + (twoo2zn)*I4[I4_ind+64]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+95] + U40*I1[I1_ind+95] + (twoo2zn)*I4[I4_ind+65]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+96] + U40*I1[I1_ind+96] + (twoo2zn)*I4[I4_ind+66]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+97] + U40*I1[I1_ind+97] + (oneo2zn)*I4[I4_ind+67]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+98] + U40*I1[I1_ind+98] + (oneo2zn)*I4[I4_ind+68]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+99] + U40*I1[I1_ind+99] + (oneo2zn)*I4[I4_ind+69]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+100] + U40*I1[I1_ind+100] + (oneo2zn)*I4[I4_ind+70]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+101] + U40*I1[I1_ind+101]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+102] + U40*I1[I1_ind+102]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+103] + U40*I1[I1_ind+103]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+104] + U40*I1[I1_ind+104]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+105] + U40*I1[I1_ind+105]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+106] + U40*I1[I1_ind+106] + (fouro2zn)*I4[I4_ind+71]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+107] + U40*I1[I1_ind+107] + (threeo2zn)*I4[I4_ind+72]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+108] + U40*I1[I1_ind+108] + (threeo2zn)*I4[I4_ind+73]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+109] + U40*I1[I1_ind+109] + (twoo2zn)*I4[I4_ind+74]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+110] + U40*I1[I1_ind+110] + (twoo2zn)*I4[I4_ind+75]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+111] + U40*I1[I1_ind+111] + (twoo2zn)*I4[I4_ind+76]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+112] + U40*I1[I1_ind+112] + (oneo2zn)*I4[I4_ind+77]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+113] + U40*I1[I1_ind+113] + (oneo2zn)*I4[I4_ind+78]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+114] + U40*I1[I1_ind+114] + (oneo2zn)*I4[I4_ind+79]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+115] + U40*I1[I1_ind+115] + (oneo2zn)*I4[I4_ind+80]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+116] + U40*I1[I1_ind+116]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+117] + U40*I1[I1_ind+117]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+118] + U40*I1[I1_ind+118]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+119] + U40*I1[I1_ind+119]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+120] + U40*I1[I1_ind+120]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+121] + U40*I1[I1_ind+121] + (fouro2zn)*I4[I4_ind+81]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+122] + U40*I1[I1_ind+122] + (threeo2zn)*I4[I4_ind+82]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+123] + U40*I1[I1_ind+123] + (threeo2zn)*I4[I4_ind+83]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+124] + U40*I1[I1_ind+124] + (twoo2zn)*I4[I4_ind+84]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+125] + U40*I1[I1_ind+125] + (twoo2zn)*I4[I4_ind+85]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+126] + U40*I1[I1_ind+126] + (twoo2zn)*I4[I4_ind+86]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+127] + U40*I1[I1_ind+127] + (oneo2zn)*I4[I4_ind+87]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+128] + U40*I1[I1_ind+128] + (oneo2zn)*I4[I4_ind+88]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+129] + U40*I1[I1_ind+129] + (oneo2zn)*I4[I4_ind+89]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+130] + U40*I1[I1_ind+130] + (oneo2zn)*I4[I4_ind+90]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+131] + U40*I1[I1_ind+131]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+132] + U40*I1[I1_ind+132]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+133] + U40*I1[I1_ind+133]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+134] + U40*I1[I1_ind+134]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+135] + U40*I1[I1_ind+135]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+136] + U40*I1[I1_ind+136] + (fouro2zn)*I4[I4_ind+91]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+137] + U40*I1[I1_ind+137] + (threeo2zn)*I4[I4_ind+92]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+138] + U40*I1[I1_ind+138] + (threeo2zn)*I4[I4_ind+93]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+139] + U40*I1[I1_ind+139] + (twoo2zn)*I4[I4_ind+94]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+140] + U40*I1[I1_ind+140] + (twoo2zn)*I4[I4_ind+95]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+141] + U40*I1[I1_ind+141] + (twoo2zn)*I4[I4_ind+96]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+142] + U40*I1[I1_ind+142] + (oneo2zn)*I4[I4_ind+97]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+143] + U40*I1[I1_ind+143] + (oneo2zn)*I4[I4_ind+98]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+144] + U40*I1[I1_ind+144] + (oneo2zn)*I4[I4_ind+99]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+145] + U40*I1[I1_ind+145] + (oneo2zn)*I4[I4_ind+100]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+146] + U40*I1[I1_ind+146]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+147] + U40*I1[I1_ind+147]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+148] + U40*I1[I1_ind+148]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+149] + U40*I1[I1_ind+149]
	vp_ind += 1
	vp[vp_ind+1] = U00*I0[I0_ind+150] + U40*I1[I1_ind+150]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+91] + U41*I1[I1_ind+91] + (threeo2z)*(I2[I2_ind+46] - (lpoz)*I3[I3_ind+46])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+92] + U41*I1[I1_ind+92] + (threeo2z)*(I2[I2_ind+47] - (lpoz)*I3[I3_ind+47]) + (oneo2zn)*I4[I4_ind+61]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+93] + U41*I1[I1_ind+93] + (threeo2z)*(I2[I2_ind+48] - (lpoz)*I3[I3_ind+48])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+94] + U41*I1[I1_ind+94] + (threeo2z)*(I2[I2_ind+49] - (lpoz)*I3[I3_ind+49]) + (twoo2zn)*I4[I4_ind+62]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+95] + U41*I1[I1_ind+95] + (threeo2z)*(I2[I2_ind+50] - (lpoz)*I3[I3_ind+50]) + (oneo2zn)*I4[I4_ind+63]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+96] + U41*I1[I1_ind+96] + (threeo2z)*(I2[I2_ind+51] - (lpoz)*I3[I3_ind+51])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+97] + U41*I1[I1_ind+97] + (threeo2z)*(I2[I2_ind+52] - (lpoz)*I3[I3_ind+52]) + (threeo2zn)*I4[I4_ind+64]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+98] + U41*I1[I1_ind+98] + (threeo2z)*(I2[I2_ind+53] - (lpoz)*I3[I3_ind+53]) + (twoo2zn)*I4[I4_ind+65]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+99] + U41*I1[I1_ind+99] + (threeo2z)*(I2[I2_ind+54] - (lpoz)*I3[I3_ind+54]) + (oneo2zn)*I4[I4_ind+66]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+100] + U41*I1[I1_ind+100] + (threeo2z)*(I2[I2_ind+55] - (lpoz)*I3[I3_ind+55])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+101] + U41*I1[I1_ind+101] + (threeo2z)*(I2[I2_ind+56] - (lpoz)*I3[I3_ind+56]) + (fouro2zn)*I4[I4_ind+67]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+102] + U41*I1[I1_ind+102] + (threeo2z)*(I2[I2_ind+57] - (lpoz)*I3[I3_ind+57]) + (threeo2zn)*I4[I4_ind+68]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+103] + U41*I1[I1_ind+103] + (threeo2z)*(I2[I2_ind+58] - (lpoz)*I3[I3_ind+58]) + (twoo2zn)*I4[I4_ind+69]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+104] + U41*I1[I1_ind+104] + (threeo2z)*(I2[I2_ind+59] - (lpoz)*I3[I3_ind+59]) + (oneo2zn)*I4[I4_ind+70]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+105] + U41*I1[I1_ind+105] + (threeo2z)*(I2[I2_ind+60] - (lpoz)*I3[I3_ind+60])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+106] + U41*I1[I1_ind+106] + (twoo2z)*(I2[I2_ind+61] - (lpoz)*I3[I3_ind+61])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+107] + U41*I1[I1_ind+107] + (twoo2z)*(I2[I2_ind+62] - (lpoz)*I3[I3_ind+62]) + (oneo2zn)*I4[I4_ind+71]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+108] + U41*I1[I1_ind+108] + (twoo2z)*(I2[I2_ind+63] - (lpoz)*I3[I3_ind+63])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+109] + U41*I1[I1_ind+109] + (twoo2z)*(I2[I2_ind+64] - (lpoz)*I3[I3_ind+64]) + (twoo2zn)*I4[I4_ind+72]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+110] + U41*I1[I1_ind+110] + (twoo2z)*(I2[I2_ind+65] - (lpoz)*I3[I3_ind+65]) + (oneo2zn)*I4[I4_ind+73]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+111] + U41*I1[I1_ind+111] + (twoo2z)*(I2[I2_ind+66] - (lpoz)*I3[I3_ind+66])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+112] + U41*I1[I1_ind+112] + (twoo2z)*(I2[I2_ind+67] - (lpoz)*I3[I3_ind+67]) + (threeo2zn)*I4[I4_ind+74]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+113] + U41*I1[I1_ind+113] + (twoo2z)*(I2[I2_ind+68] - (lpoz)*I3[I3_ind+68]) + (twoo2zn)*I4[I4_ind+75]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+114] + U41*I1[I1_ind+114] + (twoo2z)*(I2[I2_ind+69] - (lpoz)*I3[I3_ind+69]) + (oneo2zn)*I4[I4_ind+76]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+115] + U41*I1[I1_ind+115] + (twoo2z)*(I2[I2_ind+70] - (lpoz)*I3[I3_ind+70])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+116] + U41*I1[I1_ind+116] + (twoo2z)*(I2[I2_ind+71] - (lpoz)*I3[I3_ind+71]) + (fouro2zn)*I4[I4_ind+77]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+117] + U41*I1[I1_ind+117] + (twoo2z)*(I2[I2_ind+72] - (lpoz)*I3[I3_ind+72]) + (threeo2zn)*I4[I4_ind+78]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+118] + U41*I1[I1_ind+118] + (twoo2z)*(I2[I2_ind+73] - (lpoz)*I3[I3_ind+73]) + (twoo2zn)*I4[I4_ind+79]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+119] + U41*I1[I1_ind+119] + (twoo2z)*(I2[I2_ind+74] - (lpoz)*I3[I3_ind+74]) + (oneo2zn)*I4[I4_ind+80]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+120] + U41*I1[I1_ind+120] + (twoo2z)*(I2[I2_ind+75] - (lpoz)*I3[I3_ind+75])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+121] + U41*I1[I1_ind+121] + (oneo2z)*(I2[I2_ind+76] - (lpoz)*I3[I3_ind+76])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+122] + U41*I1[I1_ind+122] + (oneo2z)*(I2[I2_ind+77] - (lpoz)*I3[I3_ind+77]) + (oneo2zn)*I4[I4_ind+81]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+123] + U41*I1[I1_ind+123] + (oneo2z)*(I2[I2_ind+78] - (lpoz)*I3[I3_ind+78])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+124] + U41*I1[I1_ind+124] + (oneo2z)*(I2[I2_ind+79] - (lpoz)*I3[I3_ind+79]) + (twoo2zn)*I4[I4_ind+82]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+125] + U41*I1[I1_ind+125] + (oneo2z)*(I2[I2_ind+80] - (lpoz)*I3[I3_ind+80]) + (oneo2zn)*I4[I4_ind+83]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+126] + U41*I1[I1_ind+126] + (oneo2z)*(I2[I2_ind+81] - (lpoz)*I3[I3_ind+81])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+127] + U41*I1[I1_ind+127] + (oneo2z)*(I2[I2_ind+82] - (lpoz)*I3[I3_ind+82]) + (threeo2zn)*I4[I4_ind+84]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+128] + U41*I1[I1_ind+128] + (oneo2z)*(I2[I2_ind+83] - (lpoz)*I3[I3_ind+83]) + (twoo2zn)*I4[I4_ind+85]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+129] + U41*I1[I1_ind+129] + (oneo2z)*(I2[I2_ind+84] - (lpoz)*I3[I3_ind+84]) + (oneo2zn)*I4[I4_ind+86]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+130] + U41*I1[I1_ind+130] + (oneo2z)*(I2[I2_ind+85] - (lpoz)*I3[I3_ind+85])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+131] + U41*I1[I1_ind+131] + (oneo2z)*(I2[I2_ind+86] - (lpoz)*I3[I3_ind+86]) + (fouro2zn)*I4[I4_ind+87]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+132] + U41*I1[I1_ind+132] + (oneo2z)*(I2[I2_ind+87] - (lpoz)*I3[I3_ind+87]) + (threeo2zn)*I4[I4_ind+88]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+133] + U41*I1[I1_ind+133] + (oneo2z)*(I2[I2_ind+88] - (lpoz)*I3[I3_ind+88]) + (twoo2zn)*I4[I4_ind+89]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+134] + U41*I1[I1_ind+134] + (oneo2z)*(I2[I2_ind+89] - (lpoz)*I3[I3_ind+89]) + (oneo2zn)*I4[I4_ind+90]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+135] + U41*I1[I1_ind+135] + (oneo2z)*(I2[I2_ind+90] - (lpoz)*I3[I3_ind+90])
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+136] + U41*I1[I1_ind+136]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+137] + U41*I1[I1_ind+137] + (oneo2zn)*I4[I4_ind+91]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+138] + U41*I1[I1_ind+138]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+139] + U41*I1[I1_ind+139] + (twoo2zn)*I4[I4_ind+92]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+140] + U41*I1[I1_ind+140] + (oneo2zn)*I4[I4_ind+93]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+141] + U41*I1[I1_ind+141]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+142] + U41*I1[I1_ind+142] + (threeo2zn)*I4[I4_ind+94]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+143] + U41*I1[I1_ind+143] + (twoo2zn)*I4[I4_ind+95]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+144] + U41*I1[I1_ind+144] + (oneo2zn)*I4[I4_ind+96]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+145] + U41*I1[I1_ind+145]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+146] + U41*I1[I1_ind+146] + (fouro2zn)*I4[I4_ind+97]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+147] + U41*I1[I1_ind+147] + (threeo2zn)*I4[I4_ind+98]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+148] + U41*I1[I1_ind+148] + (twoo2zn)*I4[I4_ind+99]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+149] + U41*I1[I1_ind+149] + (oneo2zn)*I4[I4_ind+100]
	vp_ind += 1
	vp[vp_ind+1] = U01*I0[I0_ind+150] + U41*I1[I1_ind+150]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+136] + U42*I1[I1_ind+136] + (threeo2z)*(I2[I2_ind+76] - (lpoz)*I3[I3_ind+76])
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+137] + U42*I1[I1_ind+137] + (threeo2z)*(I2[I2_ind+77] - (lpoz)*I3[I3_ind+77])
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+138] + U42*I1[I1_ind+138] + (threeo2z)*(I2[I2_ind+78] - (lpoz)*I3[I3_ind+78]) + (oneo2zn)*I4[I4_ind+91]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+139] + U42*I1[I1_ind+139] + (threeo2z)*(I2[I2_ind+79] - (lpoz)*I3[I3_ind+79])
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+140] + U42*I1[I1_ind+140] + (threeo2z)*(I2[I2_ind+80] - (lpoz)*I3[I3_ind+80]) + (oneo2zn)*I4[I4_ind+92]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+141] + U42*I1[I1_ind+141] + (threeo2z)*(I2[I2_ind+81] - (lpoz)*I3[I3_ind+81]) + (twoo2zn)*I4[I4_ind+93]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+142] + U42*I1[I1_ind+142] + (threeo2z)*(I2[I2_ind+82] - (lpoz)*I3[I3_ind+82])
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+143] + U42*I1[I1_ind+143] + (threeo2z)*(I2[I2_ind+83] - (lpoz)*I3[I3_ind+83]) + (oneo2zn)*I4[I4_ind+94]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+144] + U42*I1[I1_ind+144] + (threeo2z)*(I2[I2_ind+84] - (lpoz)*I3[I3_ind+84]) + (twoo2zn)*I4[I4_ind+95]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+145] + U42*I1[I1_ind+145] + (threeo2z)*(I2[I2_ind+85] - (lpoz)*I3[I3_ind+85]) + (threeo2zn)*I4[I4_ind+96]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+146] + U42*I1[I1_ind+146] + (threeo2z)*(I2[I2_ind+86] - (lpoz)*I3[I3_ind+86])
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+147] + U42*I1[I1_ind+147] + (threeo2z)*(I2[I2_ind+87] - (lpoz)*I3[I3_ind+87]) + (oneo2zn)*I4[I4_ind+97]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+148] + U42*I1[I1_ind+148] + (threeo2z)*(I2[I2_ind+88] - (lpoz)*I3[I3_ind+88]) + (twoo2zn)*I4[I4_ind+98]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+149] + U42*I1[I1_ind+149] + (threeo2z)*(I2[I2_ind+89] - (lpoz)*I3[I3_ind+89]) + (threeo2zn)*I4[I4_ind+99]
	vp_ind += 1
	vp[vp_ind+1] = U02*I0[I0_ind+150] + U42*I1[I1_ind+150] + (threeo2z)*(I2[I2_ind+90] - (lpoz)*I3[I3_ind+90]) + (fouro2zn)*I4[I4_ind+100]
	vp_ind += 1
end
# Total number of FLOPs = 1575
