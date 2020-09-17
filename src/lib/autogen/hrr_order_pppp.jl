include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_pppp.jl")

# Computes quartets of (pp|pp) integrals

function hrr_order_pppp!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[2,2] = 0
	Libint.vrr_classes[2,3] = 9
	Libint.vrr_classes[3,2] = 27
	Libint.vrr_classes[3,3] = 45
	for i = 1:81
		int_stack[i] = 0
	end

	Libint.vrr_stack = 81
	for i = 1:num_prim_comb
		vrr_order_pppp!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (p0|pp) ---
	hrr3_build_pp!(Libint.CD, int_stack, 81, int_stack, 9, int_stack, 0, 3)
# --- compute (d0|pp) ---
	hrr3_build_pp!(Libint.CD, int_stack, 108, int_stack, 45, int_stack, 27, 6)
# --- compute (pp|pp) ---
	hrr1_build_pp!(Libint.AB, int_stack, 0, int_stack, 108, int_stack, 81, 9)
	return 0
end
