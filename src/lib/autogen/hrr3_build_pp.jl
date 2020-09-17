# These machine-generated functions compute a quartet of |pp) integrals 

function hrr3_build_pp!(CD::Vector{Float64},vp::Vector{Float64}, vp_ind::Int, I0::Vector{Float64}, I0_ind::Int, I1::Vector{Float64}, I1_ind::Int, ab_num::Int)
	CD0 = CD[1]
	CD1 = CD[2]
	CD2 = CD[3]
	for ab = 1:ab_num
		vp[vp_ind+1] = I0[I0_ind+1] + CD0*I1[I1_ind+1]
		vp_ind += 1
		vp[vp_ind+1] = I0[I0_ind+2] + CD1*I1[I1_ind+1]
		vp_ind += 1
		vp[vp_ind+1] = I0[I0_ind+3] + CD2*I1[I1_ind+1]
		vp_ind += 1
		vp[vp_ind+1] = I0[I0_ind+2] + CD0*I1[I1_ind+2]
		vp_ind += 1
		vp[vp_ind+1] = I0[I0_ind+4] + CD1*I1[I1_ind+2]
		vp_ind += 1
		vp[vp_ind+1] = I0[I0_ind+5] + CD2*I1[I1_ind+2]
		vp_ind += 1
		vp[vp_ind+1] = I0[I0_ind+3] + CD0*I1[I1_ind+3]
		vp_ind += 1
		vp[vp_ind+1] = I0[I0_ind+5] + CD1*I1[I1_ind+3]
		vp_ind += 1
		vp[vp_ind+1] = I0[I0_ind+6] + CD2*I1[I1_ind+3]
		vp_ind += 1
		I0_ind += 6
		I1_ind += 3
	end
end
# Total number of FLOPs = 18 * ab_num 
