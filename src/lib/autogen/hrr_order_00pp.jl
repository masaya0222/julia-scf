include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_00pp.jl")

# Computes quartets of (00|pp) integrals

function hrr_order_00pp!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[1,2] = 0
	Libint.vrr_classes[1,3] = 3
	for i = 1:9
		int_stack[i] = 0
	end

	Libint.vrr_stack = 9
	for i = 1:num_prim_comb
		vrr_order_00pp!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (00|pp) ---
	hrr3_build_pp!(Libint.CD, int_stack, 9, int_stack, 3, int_stack, 0, 1)
	return 9
end
