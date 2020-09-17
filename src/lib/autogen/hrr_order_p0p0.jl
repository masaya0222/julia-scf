include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_p0p0.jl")

# Computes quartets of (p0|p0) integrals

function hrr_order_p0p0!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[2,2] = 0
	for i = 1:9
		int_stack[i] = 0
	end

	Libint.vrr_stack = 9
	for i = 1:num_prim_comb
		vrr_order_p0p0!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (p0|p0) ---
	return 0
end
