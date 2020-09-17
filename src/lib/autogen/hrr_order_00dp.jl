include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_00dp.jl")

# Computes quartets of (00|dp) integrals

function hrr_order_00dp!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[1,3] = 0
	Libint.vrr_classes[1,4] = 6
	for i = 1:16
		int_stack[i] = 0
	end

	Libint.vrr_stack = 16
	for i = 1:num_prim_comb
		vrr_order_00dp!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (00|dp) ---
	hrr3_build_dp!(Libint.CD, int_stack, 16, int_stack, 6, int_stack, 0, 1)
	return 16
end
