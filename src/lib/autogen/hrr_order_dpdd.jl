include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_dpdd.jl")

# Computes quartets of (dp|dd) integrals

function hrr_order_dpdd!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[3,3] = 0
	Libint.vrr_classes[3,4] = 36
	Libint.vrr_classes[3,5] = 96
	Libint.vrr_classes[4,3] = 186
	Libint.vrr_classes[4,4] = 246
	Libint.vrr_classes[4,5] = 346
	for i = 1:496
		int_stack[i] = 0
	end

	Libint.vrr_stack = 496
	for i = 1:num_prim_comb
		vrr_order_dpdd!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (d0|dp) ---
	hrr3_build_dp!(Libint.CD, int_stack, 496, int_stack, 36, int_stack, 0, 6)
# --- compute (d0|fp) ---
	hrr3_build_fp!(Libint.CD, int_stack, 604, int_stack, 96, int_stack, 36, 6)
# --- compute (d0|dd) ---
	hrr3_build_dd!(Libint.CD, int_stack, 784, int_stack, 604, int_stack, 496, 6)
# --- compute (f0|dp) ---
	hrr3_build_dp!(Libint.CD, int_stack, 496, int_stack, 246, int_stack, 186, 10)
# --- compute (f0|fp) ---
	hrr3_build_fp!(Libint.CD, int_stack, 1000, int_stack, 346, int_stack, 246, 10)
# --- compute (f0|dd) ---
	hrr3_build_dd!(Libint.CD, int_stack, 0, int_stack, 1000, int_stack, 496, 10)
# --- compute (dp|dd) ---
	hrr1_build_dp!(Libint.AB, int_stack, 1000, int_stack, 0, int_stack, 784, 36)
	return 1000
end
