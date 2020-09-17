include("libints_header.jl")
include("vrr_header.jl")

# Computes quartets necessary to compute (dp|dp) integrals

function vrr_order_dpdp!(Libint::Libint_t, Data::prim_data)
	vrr_stack = Libint.vrr_stack
	build_p000!(Data, Libint.int_stack, vrr_stack+0, Data.F, 2, Data.F, 3, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+3, Data.F, 2, Data.F, 3, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+6, Data.F, 1, Data.F, 2, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+9, Data.F, 3, Data.F, 4, [0.0], -1, [0.0], -1, [0.0], -1)
	build_p0p0!(Data, Libint.int_stack, vrr_stack+12, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+9, [0.0], -1, [0.0], -1, Data.F, 3)
	build_p0p0!(Data, Libint.int_stack, vrr_stack+21, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+3, [0.0], -1, [0.0], -1, Data.F, 2)
	build_d0p0!(Data, Libint.int_stack, vrr_stack+30, Libint.int_stack, vrr_stack+21, Libint.int_stack, vrr_stack+12, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+0)
	build_00d0!(Data, Libint.int_stack, vrr_stack+48, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+9, Data.F, 2, Data.F, 3, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+54, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+3, Data.F, 1, Data.F, 2, [0.0], -1)
	build_p0d0!(Data, Libint.int_stack, vrr_stack+60, Libint.int_stack, vrr_stack+54, Libint.int_stack, vrr_stack+48, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+3)
	build_00p0!(Data, Libint.int_stack, vrr_stack+0, Data.F, 0, Data.F, 1, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+78, Libint.int_stack, vrr_stack+0, Libint.int_stack, vrr_stack+6, Data.F, 0, Data.F, 1, [0.0], -1)
	build_p0d0!(Data, Libint.int_stack, vrr_stack+84, Libint.int_stack, vrr_stack+78, Libint.int_stack, vrr_stack+54, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+6)
	build_00p0!(Data, Libint.int_stack, vrr_stack+102, Data.F, 4, Data.F, 5, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+105, Libint.int_stack, vrr_stack+9, Libint.int_stack, vrr_stack+102, Data.F, 3, Data.F, 4, [0.0], -1)
	build_p0d0!(Data, Libint.int_stack, vrr_stack+111, Libint.int_stack, vrr_stack+48, Libint.int_stack, vrr_stack+105, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+9)
	build_d0d0!(Data, Libint.int_stack, vrr_stack+129, Libint.int_stack, vrr_stack+60, Libint.int_stack, vrr_stack+111, Libint.int_stack, vrr_stack+54, Libint.int_stack, vrr_stack+48, Libint.int_stack, vrr_stack+12)
	build_d0d0!(Data, Libint.int_stack, vrr_stack+165, Libint.int_stack, vrr_stack+84, Libint.int_stack, vrr_stack+60, Libint.int_stack, vrr_stack+78, Libint.int_stack, vrr_stack+54, Libint.int_stack, vrr_stack+21)
	tmp = vrr_stack + 165
	target_ptr = Libint.vrr_classes[3,3]
	for i = 1:36
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
	build_00f0!(Data, Libint.int_stack, vrr_stack+12, Libint.int_stack, vrr_stack+48, Libint.int_stack, vrr_stack+105, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+9, [0.0], -1)
	build_00f0!(Data, Libint.int_stack, vrr_stack+201, Libint.int_stack, vrr_stack+54, Libint.int_stack, vrr_stack+48, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+3, [0.0], -1)
	build_p0f0!(Data, Libint.int_stack, vrr_stack+211, Libint.int_stack, vrr_stack+201, Libint.int_stack, vrr_stack+12, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+48)
	build_00f0!(Data, Libint.int_stack, vrr_stack+241, Libint.int_stack, vrr_stack+78, Libint.int_stack, vrr_stack+54, Libint.int_stack, vrr_stack+0, Libint.int_stack, vrr_stack+6, [0.0], -1)
	build_p0f0!(Data, Libint.int_stack, vrr_stack+251, Libint.int_stack, vrr_stack+241, Libint.int_stack, vrr_stack+201, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+54)
	build_00p0!(Data, Libint.int_stack, vrr_stack+0, Data.F, 5, Data.F, 6, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+102, Libint.int_stack, vrr_stack+0, Data.F, 4, Data.F, 5, [0.0], -1)
	build_00f0!(Data, Libint.int_stack, vrr_stack+48, Libint.int_stack, vrr_stack+105, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+9, Libint.int_stack, vrr_stack+102, [0.0], -1)
	build_p0f0!(Data, Libint.int_stack, vrr_stack+281, Libint.int_stack, vrr_stack+12, Libint.int_stack, vrr_stack+48, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+105)
	build_d0f0!(Data, Libint.int_stack, vrr_stack+311, Libint.int_stack, vrr_stack+211, Libint.int_stack, vrr_stack+281, Libint.int_stack, vrr_stack+201, Libint.int_stack, vrr_stack+12, Libint.int_stack, vrr_stack+111)
	build_d0f0!(Data, Libint.int_stack, vrr_stack+371, Libint.int_stack, vrr_stack+251, Libint.int_stack, vrr_stack+211, Libint.int_stack, vrr_stack+241, Libint.int_stack, vrr_stack+201, Libint.int_stack, vrr_stack+60)
	tmp = vrr_stack + 371
	target_ptr = Libint.vrr_classes[3,4]
	for i = 1:60
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
	build_f0d0!(Data, Libint.int_stack, vrr_stack+431, Libint.int_stack, vrr_stack+165, Libint.int_stack, vrr_stack+129, Libint.int_stack, vrr_stack+84, Libint.int_stack, vrr_stack+60, Libint.int_stack, vrr_stack+30)
	tmp = vrr_stack + 431
	target_ptr = Libint.vrr_classes[4,3]
	for i = 1:60
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
	build_f0f0!(Data, Libint.int_stack, vrr_stack+0, Libint.int_stack, vrr_stack+371, Libint.int_stack, vrr_stack+311, Libint.int_stack, vrr_stack+251, Libint.int_stack, vrr_stack+211, Libint.int_stack, vrr_stack+129)
	tmp = vrr_stack + 0
	target_ptr = Libint.vrr_classes[4,4]
	for i = 1:100
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
end

