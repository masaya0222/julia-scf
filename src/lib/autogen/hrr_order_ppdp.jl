include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_ppdp.jl")

# Computes quartets of (pp|dp) integrals

function hrr_order_ppdp!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[2,3] = 0
	Libint.vrr_classes[2,4] = 18
	Libint.vrr_classes[3,3] = 48
	Libint.vrr_classes[3,4] = 84
	for i = 1:144
		int_stack[i] = 0
	end

	Libint.vrr_stack = 144
	for i = 1:num_prim_comb
		vrr_order_ppdp!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (p0|dp) ---
	hrr3_build_dp!(Libint.CD, int_stack, 144, int_stack, 18, int_stack, 0, 3)
# --- compute (d0|dp) ---
	hrr3_build_dp!(Libint.CD, int_stack, 198, int_stack, 84, int_stack, 48, 6)
# --- compute (pp|dp) ---
	hrr1_build_pp!(Libint.AB, int_stack, 306, int_stack, 198, int_stack, 144, 18)
	return 306
end
