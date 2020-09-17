include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_dddd.jl")

# Computes quartets of (dd|dd) integrals

function hrr_order_dddd!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[3,3] = 0
	Libint.vrr_classes[3,4] = 36
	Libint.vrr_classes[3,5] = 96
	Libint.vrr_classes[4,3] = 186
	Libint.vrr_classes[4,4] = 246
	Libint.vrr_classes[4,5] = 346
	Libint.vrr_classes[5,3] = 496
	Libint.vrr_classes[5,4] = 586
	Libint.vrr_classes[5,5] = 736
	for i = 1:961
		int_stack[i] = 0
	end

	Libint.vrr_stack = 961
	for i = 1:num_prim_comb
		vrr_order_dddd!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (d0|dp) ---
	hrr3_build_dp!(Libint.CD, int_stack, 961, int_stack, 36, int_stack, 0, 6)
# --- compute (d0|fp) ---
	hrr3_build_fp!(Libint.CD, int_stack, 1069, int_stack, 96, int_stack, 36, 6)
# --- compute (d0|dd) ---
	hrr3_build_dd!(Libint.CD, int_stack, 1249, int_stack, 1069, int_stack, 961, 6)
# --- compute (f0|dp) ---
	hrr3_build_dp!(Libint.CD, int_stack, 961, int_stack, 246, int_stack, 186, 10)
# --- compute (f0|fp) ---
	hrr3_build_fp!(Libint.CD, int_stack, 1465, int_stack, 346, int_stack, 246, 10)
# --- compute (f0|dd) ---
	hrr3_build_dd!(Libint.CD, int_stack, 0, int_stack, 1465, int_stack, 961, 10)
# --- compute (dp|dd) ---
	hrr1_build_dp!(Libint.AB, int_stack, 1465, int_stack, 0, int_stack, 1249, 36)
# --- compute (g0|dp) ---
	hrr3_build_dp!(Libint.CD, int_stack, 961, int_stack, 586, int_stack, 496, 15)
# --- compute (g0|fp) ---
	hrr3_build_fp!(Libint.CD, int_stack, 2113, int_stack, 736, int_stack, 586, 15)
# --- compute (g0|dd) ---
	hrr3_build_dd!(Libint.CD, int_stack, 360, int_stack, 2113, int_stack, 961, 15)
# --- compute (fp|dd) ---
	hrr1_build_fp!(Libint.AB, int_stack, 2113, int_stack, 360, int_stack, 0, 36)
# --- compute (dd|dd) ---
	hrr1_build_dd!(Libint.AB, int_stack, 0, int_stack, 2113, int_stack, 1465, 36)
	return 0
end
