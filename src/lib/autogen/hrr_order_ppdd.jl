include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_ppdd.jl")

# Computes quartets of (pp|dd) integrals

function hrr_order_ppdd!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[2,3] = 0
	Libint.vrr_classes[2,4] = 18
	Libint.vrr_classes[2,5] = 48
	Libint.vrr_classes[3,3] = 93
	Libint.vrr_classes[3,4] = 129
	Libint.vrr_classes[3,5] = 189
	for i = 1:279
		int_stack[i] = 0
	end

	Libint.vrr_stack = 279
	for i = 1:num_prim_comb
		vrr_order_ppdd!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (p0|dp) ---
	hrr3_build_dp!(Libint.CD, int_stack, 279, int_stack, 18, int_stack, 0, 3)
# --- compute (p0|fp) ---
	hrr3_build_fp!(Libint.CD, int_stack, 333, int_stack, 48, int_stack, 18, 3)
# --- compute (p0|dd) ---
	hrr3_build_dd!(Libint.CD, int_stack, 423, int_stack, 333, int_stack, 279, 3)
# --- compute (d0|dp) ---
	hrr3_build_dp!(Libint.CD, int_stack, 279, int_stack, 129, int_stack, 93, 6)
# --- compute (d0|fp) ---
	hrr3_build_fp!(Libint.CD, int_stack, 531, int_stack, 189, int_stack, 129, 6)
# --- compute (d0|dd) ---
	hrr3_build_dd!(Libint.CD, int_stack, 0, int_stack, 531, int_stack, 279, 6)
# --- compute (pp|dd) ---
	hrr1_build_pp!(Libint.AB, int_stack, 531, int_stack, 0, int_stack, 423, 36)
	return 531
end
