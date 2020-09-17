module Libint

include("libints_header.jl")
include("hrr_order.jl")

build_eri = Array{Function, 4}(undef, 3, 3, 3, 3)
libint_stack_size = zeros(Int, 3)
function init_libint_base()
	build_eri[1,1,2,1] = hrr_order_00p0!
	build_eri[1,1,3,1] = hrr_order_00d0!
	build_eri[1,1,2,2] = hrr_order_00pp!
	build_eri[1,1,3,2] = hrr_order_00dp!
	build_eri[1,1,3,3] = hrr_order_00dd!
	build_eri[2,1,2,1] = hrr_order_p0p0!
	build_eri[2,1,3,1] = hrr_order_p0d0!
	build_eri[2,1,2,2] = hrr_order_p0pp!
	build_eri[2,1,3,2] = hrr_order_p0dp!
	build_eri[2,1,3,3] = hrr_order_p0dd!
	build_eri[3,1,3,1] = hrr_order_d0d0!
	build_eri[3,1,2,2] = hrr_order_d0pp!
	build_eri[3,1,3,2] = hrr_order_d0dp!
	build_eri[3,1,3,3] = hrr_order_d0dd!
	build_eri[2,2,3,1] = hrr_order_ppd0!
	build_eri[2,2,2,2] = hrr_order_pppp!
	build_eri[2,2,3,2] = hrr_order_ppdp!
	build_eri[2,2,3,3] = hrr_order_ppdd!
	build_eri[3,2,3,2] = hrr_order_dpdp!
	build_eri[3,2,3,3] = hrr_order_dpdd!
	build_eri[3,3,3,3] = hrr_order_dddd!
	libint_stack_size[1] = 1
	libint_stack_size[2] = 222
	libint_stack_size[3] = 3193
end
init_libint_base()

function init_libint(libint::Libint_t, max_am::Int)
	libint.int_stack = zeros(Float64, libint_stack_size[max_am+1])
	libint.vrr_classes = zeros(Int, (7,7))
end
export build_eri, Libint_t, prim_data, init_libint
end
