include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_p0dp.jl")

# Computes quartets of (p0|dp) integrals

function hrr_order_p0dp!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[2,3] = 0
	Libint.vrr_classes[2,4] = 18
	for i = 1:48
		int_stack[i] = 0
	end

	Libint.vrr_stack = 48
	for i = 1:num_prim_comb
		vrr_order_p0dp!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (p0|dp) ---
	hrr3_build_dp!(Libint.CD, int_stack, 48, int_stack, 18, int_stack, 0, 3)
	return 48
end
