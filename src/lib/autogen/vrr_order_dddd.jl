include("libints_header.jl")
include("vrr_header.jl")

# Computes quartets necessary to compute (dd|dd) integrals

function vrr_order_dddd!(Libint::Libint_t, Data::prim_data)
	vrr_stack = Libint.vrr_stack
	build_p000!(Data, Libint.int_stack, vrr_stack+0, Data.F, 3, Data.F, 4, [0.0], -1, [0.0], -1, [0.0], -1)
	build_p000!(Data, Libint.int_stack, vrr_stack+3, Data.F, 2, Data.F, 3, [0.0], -1, [0.0], -1, [0.0], -1)
	build_d000!(Data, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+0, Data.F, 2, Data.F, 3, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+12, Data.F, 3, Data.F, 4, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+15, Data.F, 2, Data.F, 3, [0.0], -1, [0.0], -1, [0.0], -1)
	build_p0p0!(Data, Libint.int_stack, vrr_stack+18, Libint.int_stack, vrr_stack+15, Libint.int_stack, vrr_stack+12, [0.0], -1, [0.0], -1, Data.F, 3)
	build_00p0!(Data, Libint.int_stack, vrr_stack+27, Data.F, 1, Data.F, 2, [0.0], -1, [0.0], -1, [0.0], -1)
	build_p0p0!(Data, Libint.int_stack, vrr_stack+30, Libint.int_stack, vrr_stack+27, Libint.int_stack, vrr_stack+15, [0.0], -1, [0.0], -1, Data.F, 2)
	build_00p0!(Data, Libint.int_stack, vrr_stack+39, Data.F, 4, Data.F, 5, [0.0], -1, [0.0], -1, [0.0], -1)
	build_p0p0!(Data, Libint.int_stack, vrr_stack+42, Libint.int_stack, vrr_stack+12, Libint.int_stack, vrr_stack+39, [0.0], -1, [0.0], -1, Data.F, 4)
	build_d0p0!(Data, Libint.int_stack, vrr_stack+51, Libint.int_stack, vrr_stack+18, Libint.int_stack, vrr_stack+42, Libint.int_stack, vrr_stack+15, Libint.int_stack, vrr_stack+12, Libint.int_stack, vrr_stack+0)
	build_d0p0!(Data, Libint.int_stack, vrr_stack+69, Libint.int_stack, vrr_stack+30, Libint.int_stack, vrr_stack+18, Libint.int_stack, vrr_stack+27, Libint.int_stack, vrr_stack+15, Libint.int_stack, vrr_stack+3)
	build_f0p0!(Data, Libint.int_stack, vrr_stack+87, Libint.int_stack, vrr_stack+69, Libint.int_stack, vrr_stack+51, Libint.int_stack, vrr_stack+30, Libint.int_stack, vrr_stack+18, Libint.int_stack, vrr_stack+6)
	build_00d0!(Data, Libint.int_stack, vrr_stack+0, Libint.int_stack, vrr_stack+15, Libint.int_stack, vrr_stack+12, Data.F, 2, Data.F, 3, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+27, Libint.int_stack, vrr_stack+15, Data.F, 1, Data.F, 2, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+117, Libint.int_stack, vrr_stack+12, Libint.int_stack, vrr_stack+39, Data.F, 3, Data.F, 4, [0.0], -1)
	build_p0d0!(Data, Libint.int_stack, vrr_stack+123, Libint.int_stack, vrr_stack+0, Libint.int_stack, vrr_stack+117, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+12)
	build_p0d0!(Data, Libint.int_stack, vrr_stack+141, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+0, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+15)
	build_d0d0!(Data, Libint.int_stack, vrr_stack+159, Libint.int_stack, vrr_stack+141, Libint.int_stack, vrr_stack+123, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+0, Libint.int_stack, vrr_stack+18)
	build_00p0!(Data, Libint.int_stack, vrr_stack+18, Data.F, 0, Data.F, 1, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+21, Libint.int_stack, vrr_stack+18, Libint.int_stack, vrr_stack+27, Data.F, 0, Data.F, 1, [0.0], -1)
	build_p0d0!(Data, Libint.int_stack, vrr_stack+195, Libint.int_stack, vrr_stack+21, Libint.int_stack, vrr_stack+6, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+27)
	build_d0d0!(Data, Libint.int_stack, vrr_stack+213, Libint.int_stack, vrr_stack+195, Libint.int_stack, vrr_stack+141, Libint.int_stack, vrr_stack+21, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+30)
	tmp = vrr_stack + 213
	target_ptr = Libint.vrr_classes[3,3]
	for i = 1:36
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
	build_00p0!(Data, Libint.int_stack, vrr_stack+30, Data.F, 5, Data.F, 6, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+33, Libint.int_stack, vrr_stack+39, Libint.int_stack, vrr_stack+30, Data.F, 4, Data.F, 5, [0.0], -1)
	build_p0d0!(Data, Libint.int_stack, vrr_stack+249, Libint.int_stack, vrr_stack+117, Libint.int_stack, vrr_stack+33, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+39)
	build_d0d0!(Data, Libint.int_stack, vrr_stack+267, Libint.int_stack, vrr_stack+123, Libint.int_stack, vrr_stack+249, Libint.int_stack, vrr_stack+0, Libint.int_stack, vrr_stack+117, Libint.int_stack, vrr_stack+42)
	build_f0d0!(Data, Libint.int_stack, vrr_stack+303, Libint.int_stack, vrr_stack+159, Libint.int_stack, vrr_stack+267, Libint.int_stack, vrr_stack+141, Libint.int_stack, vrr_stack+123, Libint.int_stack, vrr_stack+51)
	build_f0d0!(Data, Libint.int_stack, vrr_stack+363, Libint.int_stack, vrr_stack+213, Libint.int_stack, vrr_stack+159, Libint.int_stack, vrr_stack+195, Libint.int_stack, vrr_stack+141, Libint.int_stack, vrr_stack+69)
	tmp = vrr_stack + 363
	target_ptr = Libint.vrr_classes[4,3]
	for i = 1:60
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
	build_00f0!(Data, Libint.int_stack, vrr_stack+195, Libint.int_stack, vrr_stack+0, Libint.int_stack, vrr_stack+117, Libint.int_stack, vrr_stack+15, Libint.int_stack, vrr_stack+12, [0.0], -1)
	build_00f0!(Data, Libint.int_stack, vrr_stack+42, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+0, Libint.int_stack, vrr_stack+27, Libint.int_stack, vrr_stack+15, [0.0], -1)
	build_00f0!(Data, Libint.int_stack, vrr_stack+52, Libint.int_stack, vrr_stack+117, Libint.int_stack, vrr_stack+33, Libint.int_stack, vrr_stack+12, Libint.int_stack, vrr_stack+39, [0.0], -1)
	build_p0f0!(Data, Libint.int_stack, vrr_stack+423, Libint.int_stack, vrr_stack+195, Libint.int_stack, vrr_stack+52, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+117)
	build_p0f0!(Data, Libint.int_stack, vrr_stack+453, Libint.int_stack, vrr_stack+42, Libint.int_stack, vrr_stack+195, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+0)
	build_d0f0!(Data, Libint.int_stack, vrr_stack+483, Libint.int_stack, vrr_stack+453, Libint.int_stack, vrr_stack+423, Libint.int_stack, vrr_stack+42, Libint.int_stack, vrr_stack+195, Libint.int_stack, vrr_stack+123)
	build_00f0!(Data, Libint.int_stack, vrr_stack+123, Libint.int_stack, vrr_stack+21, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+18, Libint.int_stack, vrr_stack+27, [0.0], -1)
	build_p0f0!(Data, Libint.int_stack, vrr_stack+543, Libint.int_stack, vrr_stack+123, Libint.int_stack, vrr_stack+42, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+6)
	build_d0f0!(Data, Libint.int_stack, vrr_stack+573, Libint.int_stack, vrr_stack+543, Libint.int_stack, vrr_stack+453, Libint.int_stack, vrr_stack+123, Libint.int_stack, vrr_stack+42, Libint.int_stack, vrr_stack+141)
	tmp = vrr_stack + 573
	target_ptr = Libint.vrr_classes[3,4]
	for i = 1:60
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
	build_00p0!(Data, Libint.int_stack, vrr_stack+27, Data.F, 6, Data.F, 7, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+133, Libint.int_stack, vrr_stack+30, Libint.int_stack, vrr_stack+27, Data.F, 5, Data.F, 6, [0.0], -1)
	build_00f0!(Data, Libint.int_stack, vrr_stack+139, Libint.int_stack, vrr_stack+33, Libint.int_stack, vrr_stack+133, Libint.int_stack, vrr_stack+39, Libint.int_stack, vrr_stack+30, [0.0], -1)
	build_p0f0!(Data, Libint.int_stack, vrr_stack+633, Libint.int_stack, vrr_stack+52, Libint.int_stack, vrr_stack+139, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+33)
	build_d0f0!(Data, Libint.int_stack, vrr_stack+663, Libint.int_stack, vrr_stack+423, Libint.int_stack, vrr_stack+633, Libint.int_stack, vrr_stack+195, Libint.int_stack, vrr_stack+52, Libint.int_stack, vrr_stack+249)
	build_f0f0!(Data, Libint.int_stack, vrr_stack+723, Libint.int_stack, vrr_stack+483, Libint.int_stack, vrr_stack+663, Libint.int_stack, vrr_stack+453, Libint.int_stack, vrr_stack+423, Libint.int_stack, vrr_stack+267)
	build_f0f0!(Data, Libint.int_stack, vrr_stack+823, Libint.int_stack, vrr_stack+573, Libint.int_stack, vrr_stack+483, Libint.int_stack, vrr_stack+543, Libint.int_stack, vrr_stack+453, Libint.int_stack, vrr_stack+159)
	tmp = vrr_stack + 823
	target_ptr = Libint.vrr_classes[4,4]
	for i = 1:100
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
	build_00g0!(Data, Libint.int_stack, vrr_stack+543, Libint.int_stack, vrr_stack+195, Libint.int_stack, vrr_stack+52, Libint.int_stack, vrr_stack+0, Libint.int_stack, vrr_stack+117, [0.0], -1)
	build_00g0!(Data, Libint.int_stack, vrr_stack+558, Libint.int_stack, vrr_stack+42, Libint.int_stack, vrr_stack+195, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+0, [0.0], -1)
	build_00g0!(Data, Libint.int_stack, vrr_stack+249, Libint.int_stack, vrr_stack+52, Libint.int_stack, vrr_stack+139, Libint.int_stack, vrr_stack+117, Libint.int_stack, vrr_stack+33, [0.0], -1)
	build_p0g0!(Data, Libint.int_stack, vrr_stack+923, Libint.int_stack, vrr_stack+543, Libint.int_stack, vrr_stack+249, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+52)
	build_p0g0!(Data, Libint.int_stack, vrr_stack+968, Libint.int_stack, vrr_stack+558, Libint.int_stack, vrr_stack+543, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+195)
	build_d0g0!(Data, Libint.int_stack, vrr_stack+1013, Libint.int_stack, vrr_stack+968, Libint.int_stack, vrr_stack+923, Libint.int_stack, vrr_stack+558, Libint.int_stack, vrr_stack+543, Libint.int_stack, vrr_stack+423)
	build_00g0!(Data, Libint.int_stack, vrr_stack+423, Libint.int_stack, vrr_stack+123, Libint.int_stack, vrr_stack+42, Libint.int_stack, vrr_stack+21, Libint.int_stack, vrr_stack+6, [0.0], -1)
	build_p0g0!(Data, Libint.int_stack, vrr_stack+1103, Libint.int_stack, vrr_stack+423, Libint.int_stack, vrr_stack+558, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+42)
	build_d0g0!(Data, Libint.int_stack, vrr_stack+1148, Libint.int_stack, vrr_stack+1103, Libint.int_stack, vrr_stack+968, Libint.int_stack, vrr_stack+423, Libint.int_stack, vrr_stack+558, Libint.int_stack, vrr_stack+453)
	tmp = vrr_stack + 1148
	target_ptr = Libint.vrr_classes[3,5]
	for i = 1:90
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
	build_00p0!(Data, Libint.int_stack, vrr_stack+558, Data.F, 7, Data.F, 8, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+561, Libint.int_stack, vrr_stack+27, Libint.int_stack, vrr_stack+558, Data.F, 6, Data.F, 7, [0.0], -1)
	build_00f0!(Data, Libint.int_stack, vrr_stack+149, Libint.int_stack, vrr_stack+133, Libint.int_stack, vrr_stack+561, Libint.int_stack, vrr_stack+30, Libint.int_stack, vrr_stack+27, [0.0], -1)
	build_00g0!(Data, Libint.int_stack, vrr_stack+558, Libint.int_stack, vrr_stack+139, Libint.int_stack, vrr_stack+149, Libint.int_stack, vrr_stack+33, Libint.int_stack, vrr_stack+133, [0.0], -1)
	build_p0g0!(Data, Libint.int_stack, vrr_stack+423, Libint.int_stack, vrr_stack+249, Libint.int_stack, vrr_stack+558, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+139)
	build_d0g0!(Data, Libint.int_stack, vrr_stack+1238, Libint.int_stack, vrr_stack+923, Libint.int_stack, vrr_stack+423, Libint.int_stack, vrr_stack+543, Libint.int_stack, vrr_stack+249, Libint.int_stack, vrr_stack+633)
	build_f0g0!(Data, Libint.int_stack, vrr_stack+1328, Libint.int_stack, vrr_stack+1013, Libint.int_stack, vrr_stack+1238, Libint.int_stack, vrr_stack+968, Libint.int_stack, vrr_stack+923, Libint.int_stack, vrr_stack+663)
	build_f0g0!(Data, Libint.int_stack, vrr_stack+1478, Libint.int_stack, vrr_stack+1148, Libint.int_stack, vrr_stack+1013, Libint.int_stack, vrr_stack+1103, Libint.int_stack, vrr_stack+968, Libint.int_stack, vrr_stack+483)
	tmp = vrr_stack + 1478
	target_ptr = Libint.vrr_classes[4,5]
	for i = 1:150
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
	build_g0d0!(Data, Libint.int_stack, vrr_stack+923, Libint.int_stack, vrr_stack+363, Libint.int_stack, vrr_stack+303, Libint.int_stack, vrr_stack+213, Libint.int_stack, vrr_stack+159, Libint.int_stack, vrr_stack+87)
	tmp = vrr_stack + 923
	target_ptr = Libint.vrr_classes[5,3]
	for i = 1:90
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
	build_g0f0!(Data, Libint.int_stack, vrr_stack+0, Libint.int_stack, vrr_stack+823, Libint.int_stack, vrr_stack+723, Libint.int_stack, vrr_stack+573, Libint.int_stack, vrr_stack+483, Libint.int_stack, vrr_stack+303)
	tmp = vrr_stack + 0
	target_ptr = Libint.vrr_classes[5,4]
	for i = 1:150
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
	build_g0g0!(Data, Libint.int_stack, vrr_stack+150, Libint.int_stack, vrr_stack+1478, Libint.int_stack, vrr_stack+1328, Libint.int_stack, vrr_stack+1148, Libint.int_stack, vrr_stack+1013, Libint.int_stack, vrr_stack+723)
	tmp = vrr_stack + 150
	target_ptr = Libint.vrr_classes[5,5]
	for i = 1:225
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
end

