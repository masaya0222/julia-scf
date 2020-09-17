import constant
from constant import hash, io
import hrr_build
import order
import vrr_build
import numpy as np

lib = constant.lib
new_am = constant.new_am
am = constant.am
libint = "libints.jl"

hrr_build.hrr_build()
f = open(lib+"/"+libint, 'w')
f.write('module Libint\n\
\n\
include("libints_header.jl")\n\
include("hrr_order.jl")\n\
\n\
build_eri = Array{Function, 4}(undef, %d, %d, %d, %d)\n\
libint_stack_size = zeros(Int, %d)\n\
function init_libint_base()\n'%(am+1,am+1,am+1,am+1,am+1))

libint_stack_size = np.ones(constant.MAX_AM//2+1)

order.emit_order(libint_stack_size)
vrr_build.vrr_build()
hrr_build.hrr_build()

new_am = constant.new_am
am_leter = constant.am_letter
for la in range(new_am+1):
    lb_max = la//2
    lb_min = 0 if la <= new_am//2 else la - new_am//2
    lc_min = 1 if la == 0 else la
    for lb in range(lb_min,lb_max+1):
        for lc in range(lc_min, new_am+1):
            ld_max = lc//2
            ld_min = 0 if lc <= new_am//2 else lc - new_am//2
            for ld in range(ld_min, ld_max+1):
                hrr_function_name = "hrr_order_"+am_leter[la-lb]+am_leter[lb]+am_leter[lc-ld]+am_leter[ld]+"!"
                f.write(f"\tbuild_eri[{la-lb+1},{lb+1},{lc-ld+1},{ld+1}] = {hrr_function_name}\n")
for l in range(am+1):
    f.write("\tlibint_stack_size[%d] = %d\n"%(l+1,libint_stack_size[l]))
f.write("end\n")
f.write("init_libint_base()\n\n")
f.write("function init_libint(libint::Libint_t, max_am::Int)\n")#, max_num_prim_quartets::Int)\n")
f.write("\tlibint.int_stack = zeros(Float64, libint_stack_size[max_am+1])\n")
#f.write("\tlibint.PrimQuartet = [prim_data() for i = 1:max_num_prim_quartets]\n")
f.write("\tlibint.vrr_classes = zeros(Int, (7,7))\n")
#f.write("\tfor i = 1:max_num_prim_quartets\n")
#f.write("\t\tlibint.PrimQuartet[i].U = zeros(Float64,(6,3))\n")
#f.write("\tend\n")
f.write("end\n")
f.write("export build_eri, Libint_t, prim_data, init_libint\n")
f.write("end\n")
f.close()
libints_header = "libints_header.jl"
f = open(lib+"/"+libints_header,'w')
f.write('mutable struct prim_data \n\
    F::Vector{Float64}\n\
    U::Matrix{Float64}\n\
    twozeta_a::Float64\n\
    twozeta_b::Float64\n\
    twozeta_c::Float64\n\
    twozeta_d::Float64\n\
    oo2z::Float64\n\
    oo2n::Float64\n\
    oo2zn::Float64\n\
    poz::Float64\n\
    pon::Float64\n\
    oo2p::Float64\n\
    ss_r12_ss::Float64\n\
    prim_data() = new()\n\
end\n\
\n\
mutable struct  Libint_t\n\
    int_stack::Vector{Float64}\n\
    PrimQuartet::Vector{prim_data}\n\
    AB::Vector{Float64}\n\
    CD::Vector{Float64}\n\
    vrr_classes::Matrix{Int}\n\
    vrr_stack::Int\n\
    Libint_t() = new()\n\
end\n')
f.close()



