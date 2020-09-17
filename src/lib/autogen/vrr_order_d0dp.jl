include("libints_header.jl")
include("vrr_header.jl")

# Computes quartets necessary to compute (d0|dp) integrals

function vrr_order_d0dp!(Libint::Libint_t, Data::prim_data)
	vrr_stack = Libint.vrr_stack
	build_00p0!(Data, Libint.int_stack, vrr_stack+0, Data.F, 2, Data.F, 3, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+3, Data.F, 1, Data.F, 2, [0.0], -1, [0.0], -1, [0.0], -1)
	build_p0p0!(Data, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+0, [0.0], -1, [0.0], -1, Data.F, 2)
	build_00d0!(Data, Libint.int_stack, vrr_stack+15, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+0, Data.F, 1, Data.F, 2, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+21, Data.F, 0, Data.F, 1, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+24, Libint.int_stack, vrr_stack+21, Libint.int_stack, vrr_stack+3, Data.F, 0, Data.F, 1, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+30, Data.F, 3, Data.F, 4, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+33, Libint.int_stack, vrr_stack+0, Libint.int_stack, vrr_stack+30, Data.F, 2, Data.F, 3, [0.0], -1)
	build_p0d0!(Data, Libint.int_stack, vrr_stack+39, Libint.int_stack, vrr_stack+15, Libint.int_stack, vrr_stack+33, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+0)
	build_p0d0!(Data, Libint.int_stack, vrr_stack+57, Libint.int_stack, vrr_stack+24, Libint.int_stack, vrr_stack+15, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+3)
	build_00f0!(Data, Libint.int_stack, vrr_stack+75, Libint.int_stack, vrr_stack+15, Libint.int_stack, vrr_stack+33, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+0, [0.0], -1)
	build_00f0!(Data, Libint.int_stack, vrr_stack+85, Libint.int_stack, vrr_stack+24, Libint.int_stack, vrr_stack+15, Libint.int_stack, vrr_stack+21, Libint.int_stack, vrr_stack+3, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+3, Data.F, 4, Data.F, 5, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+95, Libint.int_stack, vrr_stack+30, Libint.int_stack, vrr_stack+3, Data.F, 3, Data.F, 4, [0.0], -1)
	build_00f0!(Data, Libint.int_stack, vrr_stack+101, Libint.int_stack, vrr_stack+33, Libint.int_stack, vrr_stack+95, Libint.int_stack, vrr_stack+0, Libint.int_stack, vrr_stack+30, [0.0], -1)
	build_p0f0!(Data, Libint.int_stack, vrr_stack+111, Libint.int_stack, vrr_stack+75, Libint.int_stack, vrr_stack+101, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+33)
	build_p0f0!(Data, Libint.int_stack, vrr_stack+141, Libint.int_stack, vrr_stack+85, Libint.int_stack, vrr_stack+75, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+15)
	build_d0d0!(Data, Libint.int_stack, vrr_stack+171, Libint.int_stack, vrr_stack+57, Libint.int_stack, vrr_stack+39, Libint.int_stack, vrr_stack+24, Libint.int_stack, vrr_stack+15, Libint.int_stack, vrr_stack+6)
	tmp = vrr_stack + 171
	target_ptr = Libint.vrr_classes[3,3]
	for i = 1:36
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
	build_d0f0!(Data, Libint.int_stack, vrr_stack+207, Libint.int_stack, vrr_stack+141, Libint.int_stack, vrr_stack+111, Libint.int_stack, vrr_stack+85, Libint.int_stack, vrr_stack+75, Libint.int_stack, vrr_stack+39)
	tmp = vrr_stack + 207
	target_ptr = Libint.vrr_classes[3,4]
	for i = 1:60
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
end

