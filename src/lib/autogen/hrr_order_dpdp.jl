include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_dpdp.jl")

# Computes quartets of (dp|dp) integrals

function hrr_order_dpdp!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[3,3] = 0
	Libint.vrr_classes[3,4] = 36
	Libint.vrr_classes[4,3] = 96
	Libint.vrr_classes[4,4] = 156
	for i = 1:256
		int_stack[i] = 0
	end

	Libint.vrr_stack = 256
	for i = 1:num_prim_comb
		vrr_order_dpdp!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (d0|dp) ---
	hrr3_build_dp!(Libint.CD, int_stack, 256, int_stack, 36, int_stack, 0, 6)
# --- compute (f0|dp) ---
	hrr3_build_dp!(Libint.CD, int_stack, 364, int_stack, 156, int_stack, 96, 10)
# --- compute (dp|dp) ---
	hrr1_build_dp!(Libint.AB, int_stack, 544, int_stack, 364, int_stack, 256, 18)
	return 544
end
