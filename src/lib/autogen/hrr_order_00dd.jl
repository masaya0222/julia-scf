include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_00dd.jl")

# Computes quartets of (00|dd) integrals

function hrr_order_00dd!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[1,3] = 0
	Libint.vrr_classes[1,4] = 6
	Libint.vrr_classes[1,5] = 16
	for i = 1:31
		int_stack[i] = 0
	end

	Libint.vrr_stack = 31
	for i = 1:num_prim_comb
		vrr_order_00dd!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (00|dp) ---
	hrr3_build_dp!(Libint.CD, int_stack, 31, int_stack, 6, int_stack, 0, 1)
# --- compute (00|fp) ---
	hrr3_build_fp!(Libint.CD, int_stack, 49, int_stack, 16, int_stack, 6, 1)
# --- compute (00|dd) ---
	hrr3_build_dd!(Libint.CD, int_stack, 79, int_stack, 49, int_stack, 31, 1)
	return 79
end
