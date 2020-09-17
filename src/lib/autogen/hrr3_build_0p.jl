# These machine-generated functions compute a quartet of |0p) integrals 

function hrr3_build_0p!(CD::Vector{Float64},vp::Vector{Float64}, vp_ind::Int, I0::Vector{Float64}, I0_ind::Int, I1::Vector{Float64}, I1_ind::Int, ab_num::Int)
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
		I0_ind += 3
		I1_ind += 1
	end
end
# Total number of FLOPs = 6 * ab_num 
