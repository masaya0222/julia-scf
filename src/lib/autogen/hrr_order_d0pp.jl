include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_d0pp.jl")

# Computes quartets of (d0|pp) integrals

function hrr_order_d0pp!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[3,2] = 0
	Libint.vrr_classes[3,3] = 18
	for i = 1:54
		int_stack[i] = 0
	end

	Libint.vrr_stack = 54
	for i = 1:num_prim_comb
		vrr_order_d0pp!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (d0|pp) ---
	hrr3_build_pp!(Libint.CD, int_stack, 54, int_stack, 18, int_stack, 0, 6)
	return 54
end
