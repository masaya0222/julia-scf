include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_p0pp.jl")

# Computes quartets of (p0|pp) integrals

function hrr_order_p0pp!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[2,2] = 0
	Libint.vrr_classes[2,3] = 9
	for i = 1:27
		int_stack[i] = 0
	end

	Libint.vrr_stack = 27
	for i = 1:num_prim_comb
		vrr_order_p0pp!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (p0|pp) ---
	hrr3_build_pp!(Libint.CD, int_stack, 27, int_stack, 9, int_stack, 0, 3)
	return 27
end
