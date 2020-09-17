# These machine-generated functions compute a quartet of |pp) integrals 

function hrr1_build_pp!(AB::Vector{Float64},vp::Vector{Float64}, vp_ind::Int, I0::Vector{Float64}, I0_ind::Int, I1::Vector{Float64}, I1_ind::Int, cd_num::Int)
	AB0 = AB[1]
	AB1 = AB[2]
	AB2 = AB[3]

	i0_ind = I0_ind
	i1_ind = I1_ind
	for cd = 1:cd_num
		vp[vp_ind+1] = I0[i0_ind+1] + AB0*I1[i1_ind+1]
		vp_ind += 1; i0_ind += 1; i1_ind += 1
	end
	i0_ind = I0_ind + 1*cd_num
	i1_ind = I1_ind
	for cd = 1:cd_num
		vp[vp_ind+1] = I0[i0_ind+1] + AB1*I1[i1_ind+1]
		vp_ind += 1; i0_ind += 1; i1_ind += 1
	end
	i0_ind = I0_ind + 2*cd_num
	i1_ind = I1_ind
	for cd = 1:cd_num
		vp[vp_ind+1] = I0[i0_ind+1] + AB2*I1[i1_ind+1]
		vp_ind += 1; i0_ind += 1; i1_ind += 1
	end
	i0_ind = I0_ind + 1*cd_num
	i1_ind = I1_ind + 1*cd_num
	for cd = 1:cd_num
		vp[vp_ind+1] = I0[i0_ind+1] + AB0*I1[i1_ind+1]
		vp_ind += 1; i0_ind += 1; i1_ind += 1
	end
	i0_ind = I0_ind + 3*cd_num
	i1_ind = I1_ind + 1*cd_num
	for cd = 1:cd_num
		vp[vp_ind+1] = I0[i0_ind+1] + AB1*I1[i1_ind+1]
		vp_ind += 1; i0_ind += 1; i1_ind += 1
	end
	i0_ind = I0_ind + 4*cd_num
	i1_ind = I1_ind + 1*cd_num
	for cd = 1:cd_num
		vp[vp_ind+1] = I0[i0_ind+1] + AB2*I1[i1_ind+1]
		vp_ind += 1; i0_ind += 1; i1_ind += 1
	end
	i0_ind = I0_ind + 2*cd_num
	i1_ind = I1_ind + 2*cd_num
	for cd = 1:cd_num
		vp[vp_ind+1] = I0[i0_ind+1] + AB0*I1[i1_ind+1]
		vp_ind += 1; i0_ind += 1; i1_ind += 1
	end
	i0_ind = I0_ind + 4*cd_num
	i1_ind = I1_ind + 2*cd_num
	for cd = 1:cd_num
		vp[vp_ind+1] = I0[i0_ind+1] + AB1*I1[i1_ind+1]
		vp_ind += 1; i0_ind += 1; i1_ind += 1
	end
	i0_ind = I0_ind + 5*cd_num
	i1_ind = I1_ind + 2*cd_num
	for cd = 1:cd_num
		vp[vp_ind+1] = I0[i0_ind+1] + AB2*I1[i1_ind+1]
		vp_ind += 1; i0_ind += 1; i1_ind += 1
	end
end
# Total number of FLOPs = 18 * cd_num 
