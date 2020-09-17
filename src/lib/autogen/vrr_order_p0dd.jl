include("libints_header.jl")
include("vrr_header.jl")

# Computes quartets necessary to compute (p0|dd) integrals

function vrr_order_p0dd!(Libint::Libint_t, Data::prim_data)
	vrr_stack = Libint.vrr_stack
	build_00p0!(Data, Libint.int_stack, vrr_stack+0, Data.F, 1, Data.F, 2, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+3, Data.F, 2, Data.F, 3, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+0, Libint.int_stack, vrr_stack+3, Data.F, 1, Data.F, 2, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+12, Data.F, 0, Data.F, 1, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+15, Libint.int_stack, vrr_stack+12, Libint.int_stack, vrr_stack+0, Data.F, 0, Data.F, 1, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+21, Data.F, 3, Data.F, 4, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+24, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+21, Data.F, 2, Data.F, 3, [0.0], -1)
	build_00f0!(Data, Libint.int_stack, vrr_stack+30, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+24, Libint.int_stack, vrr_stack+0, Libint.int_stack, vrr_stack+3, [0.0], -1)
	build_00f0!(Data, Libint.int_stack, vrr_stack+40, Libint.int_stack, vrr_stack+15, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+12, Libint.int_stack, vrr_stack+0, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+12, Data.F, 4, Data.F, 5, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+50, Libint.int_stack, vrr_stack+21, Libint.int_stack, vrr_stack+12, Data.F, 3, Data.F, 4, [0.0], -1)
	build_00f0!(Data, Libint.int_stack, vrr_stack+56, Libint.int_stack, vrr_stack+24, Libint.int_stack, vrr_stack+50, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+21, [0.0], -1)
	build_00g0!(Data, Libint.int_stack, vrr_stack+66, Libint.int_stack, vrr_stack+30, Libint.int_stack, vrr_stack+56, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+24, [0.0], -1)
	build_00g0!(Data, Libint.int_stack, vrr_stack+50, Libint.int_stack, vrr_stack+40, Libint.int_stack, vrr_stack+30, Libint.int_stack, vrr_stack+15, Libint.int_stack, vrr_stack+6, [0.0], -1)
	build_p0d0!(Data, Libint.int_stack, vrr_stack+81, Libint.int_stack, vrr_stack+15, Libint.int_stack, vrr_stack+6, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+0)
	tmp = vrr_stack + 81
	target_ptr = Libint.vrr_classes[2,3]
	for i = 1:18
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
	build_p0f0!(Data, Libint.int_stack, vrr_stack+99, Libint.int_stack, vrr_stack+40, Libint.int_stack, vrr_stack+30, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+6)
	tmp = vrr_stack + 99
	target_ptr = Libint.vrr_classes[2,4]
	for i = 1:30
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
	build_p0g0!(Data, Libint.int_stack, vrr_stack+129, Libint.int_stack, vrr_stack+50, Libint.int_stack, vrr_stack+66, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+30)
	tmp = vrr_stack + 129
	target_ptr = Libint.vrr_classes[2,5]
	for i = 1:45
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
end

