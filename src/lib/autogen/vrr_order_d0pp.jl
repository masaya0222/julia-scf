include("libints_header.jl")
include("vrr_header.jl")

# Computes quartets necessary to compute (d0|pp) integrals

function vrr_order_d0pp!(Libint::Libint_t, Data::prim_data)
	vrr_stack = Libint.vrr_stack
	build_p000!(Data, Libint.int_stack, vrr_stack+0, Data.F, 1, Data.F, 2, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+3, Data.F, 1, Data.F, 2, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+6, Data.F, 0, Data.F, 1, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+9, Data.F, 2, Data.F, 3, [0.0], -1, [0.0], -1, [0.0], -1)
	build_p0p0!(Data, Libint.int_stack, vrr_stack+12, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+9, [0.0], -1, [0.0], -1, Data.F, 2)
	build_p0p0!(Data, Libint.int_stack, vrr_stack+21, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+3, [0.0], -1, [0.0], -1, Data.F, 1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+30, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+9, Data.F, 1, Data.F, 2, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+36, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+3, Data.F, 0, Data.F, 1, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+42, Data.F, 3, Data.F, 4, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+45, Libint.int_stack, vrr_stack+9, Libint.int_stack, vrr_stack+42, Data.F, 2, Data.F, 3, [0.0], -1)
	build_p0d0!(Data, Libint.int_stack, vrr_stack+51, Libint.int_stack, vrr_stack+30, Libint.int_stack, vrr_stack+45, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+9)
	build_p0d0!(Data, Libint.int_stack, vrr_stack+69, Libint.int_stack, vrr_stack+36, Libint.int_stack, vrr_stack+30, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+3)
	build_d0p0!(Data, Libint.int_stack, vrr_stack+87, Libint.int_stack, vrr_stack+21, Libint.int_stack, vrr_stack+12, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+0)
	tmp = vrr_stack + 87
	target_ptr = Libint.vrr_classes[3,2]
	for i = 1:18
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
	build_d0d0!(Data, Libint.int_stack, vrr_stack+105, Libint.int_stack, vrr_stack+69, Libint.int_stack, vrr_stack+51, Libint.int_stack, vrr_stack+36, Libint.int_stack, vrr_stack+30, Libint.int_stack, vrr_stack+12)
	tmp = vrr_stack + 105
	target_ptr = Libint.vrr_classes[3,3]
	for i = 1:36
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
end

