include("libints_header.jl")
include("hrr_header.jl")
include("vrr_order_00p0.jl")

# Computes quartets of (00|p0) integrals

function hrr_order_00p0!(Libint::Libint_t, num_prim_comb::Int)
	Data = Libint.PrimQuartet
	Data_ind = 1
	int_stack = Libint.int_stack
	Libint.vrr_classes[1,2] = 0
	for i = 1:3
		int_stack[i] = 0
	end

	Libint.vrr_stack = 3
	for i = 1:num_prim_comb
		vrr_order_00p0!(Libint, Data[Data_ind])
		Data_ind += 1
	end
# --- compute (00|p0) ---
	return 0
end
