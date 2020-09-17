include("libints_header.jl")
include("vrr_header.jl")

# Computes quartets necessary to compute (p0|pp) integrals

function vrr_order_p0pp!(Libint::Libint_t, Data::prim_data)
	vrr_stack = Libint.vrr_stack
	build_00p0!(Data, Libint.int_stack, vrr_stack+0, Data.F, 1, Data.F, 2, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+3, Data.F, 0, Data.F, 1, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+6, Data.F, 2, Data.F, 3, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+9, Libint.int_stack, vrr_stack+0, Libint.int_stack, vrr_stack+6, Data.F, 1, Data.F, 2, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+15, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+0, Data.F, 0, Data.F, 1, [0.0], -1)
	build_p0p0!(Data, Libint.int_stack, vrr_stack+21, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+0, [0.0], -1, [0.0], -1, Data.F, 1)
	tmp = vrr_stack + 21
	target_ptr = Libint.vrr_classes[2,2]
	for i = 1:9
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
	build_p0d0!(Data, Libint.int_stack, vrr_stack+30, Libint.int_stack, vrr_stack+15, Libint.int_stack, vrr_stack+9, [0.0], -1, [0.0], -1, Libint.int_stack, vrr_stack+0)
	tmp = vrr_stack + 30
	target_ptr = Libint.vrr_classes[2,3]
	for i = 1:18
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
end

