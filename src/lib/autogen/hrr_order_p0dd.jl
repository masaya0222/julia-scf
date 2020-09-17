include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_p0dd.jl")

# Computes quartets of (p0|dd) integrals

function hrr_order_p0dd!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[2,3] = 0
	Libint.vrr_classes[2,4] = 18
	Libint.vrr_classes[2,5] = 48
	for i = 1:93
		int_stack[i] = 0
	end

	Libint.vrr_stack = 93
	for i = 1:num_prim_comb
		vrr_order_p0dd!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (p0|dp) ---
	hrr3_build_dp!(Libint.CD, int_stack, 93, int_stack, 18, int_stack, 0, 3)
# --- compute (p0|fp) ---
	hrr3_build_fp!(Libint.CD, int_stack, 147, int_stack, 48, int_stack, 18, 3)
# --- compute (p0|dd) ---
	hrr3_build_dd!(Libint.CD, int_stack, 237, int_stack, 147, int_stack, 93, 3)
	return 237
end
