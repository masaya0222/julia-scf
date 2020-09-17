include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_ppd0.jl")

# Computes quartets of (pp|d0) integrals

function hrr_order_ppd0!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[2,3] = 0
	Libint.vrr_classes[3,3] = 18
	for i = 1:54
		int_stack[i] = 0
	end

	Libint.vrr_stack = 54
	for i = 1:num_prim_comb
		vrr_order_ppd0!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (pp|d0) ---
	hrr1_build_pp!(Libint.AB, int_stack, 54, int_stack, 18, int_stack, 0, 6)
	return 54
end
