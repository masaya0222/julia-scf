include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_d0d0.jl")

# Computes quartets of (d0|d0) integrals

function hrr_order_d0d0!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[3,3] = 0
	for i = 1:36
		int_stack[i] = 0
	end

	Libint.vrr_stack = 36
	for i = 1:num_prim_comb
		vrr_order_d0d0!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (d0|d0) ---
	return 0
end
