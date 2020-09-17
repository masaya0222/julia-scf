include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_d0dd.jl")

# Computes quartets of (d0|dd) integrals

function hrr_order_d0dd!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[3,3] = 0
	Libint.vrr_classes[3,4] = 36
	Libint.vrr_classes[3,5] = 96
	for i = 1:186
		int_stack[i] = 0
	end

	Libint.vrr_stack = 186
	for i = 1:num_prim_comb
		vrr_order_d0dd!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (d0|dp) ---
	hrr3_build_dp!(Libint.CD, int_stack, 186, int_stack, 36, int_stack, 0, 6)
# --- compute (d0|fp) ---
	hrr3_build_fp!(Libint.CD, int_stack, 294, int_stack, 96, int_stack, 36, 6)
# --- compute (d0|dd) ---
	hrr3_build_dd!(Libint.CD, int_stack, 474, int_stack, 294, int_stack, 186, 6)
	return 474
end
