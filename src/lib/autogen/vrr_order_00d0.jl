include("libints_header.jl")
include("vrr_header.jl")

# Computes quartets necessary to compute (00|d0) integrals

function vrr_order_00d0!(Libint::Libint_t, Data::prim_data)
	vrr_stack = Libint.vrr_stack
	build_00p0!(Data, Libint.int_stack, vrr_stack+0, Data.F, 1, Data.F, 2, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00p0!(Data, Libint.int_stack, vrr_stack+3, Data.F, 0, Data.F, 1, [0.0], -1, [0.0], -1, [0.0], -1)
	build_00d0!(Data, Libint.int_stack, vrr_stack+6, Libint.int_stack, vrr_stack+3, Libint.int_stack, vrr_stack+0, Data.F, 0, Data.F, 1, [0.0], -1)
	tmp = vrr_stack + 6
	target_ptr = Libint.vrr_classes[1,3]
	for i = 1:6
		Libint.int_stack[target_ptr+i] += Libint.int_stack[tmp+i]
	end
end

