import constant
from constant import hash, io
import numpy as np

k1 = np.zeros(constant.new_am)
k2 = np.zeros(constant.new_am)
k3 = np.zeros(constant.new_am)

def vrr_build():
    global k1, k2, k3
    new_am = constant.opt_am
    old_am = constant.old_am
    am_letter = constant.am_letter
    number = constant.number
    lib = constant.lib
    am_in = [0,0]
    am = np.zeros((2,3),dtype=int)
    k4 = ["lpoz", "lpon"]
    k1_suff = "o2z"
    k2_suff = "o2zn"
    k3_suff = "o2n"
    max1, max2 = 0, 0
    k1 = [number[i] + k1_suff for i in range(1, new_am+1)]
    k2 = [number[i] + k2_suff for i in range(1, new_am+1)]
    k3 = [number[i] + k3_suff for i in range(1, new_am+1)]
    k1max, k2max, k3max = (0,0,0)

    vrr_header = open(lib+"/"+"vrr_header.jl",'w')

    for la in range(new_am+1):
        lc_min = 0 if (la >= old_am + 1) else old_am + 1
        lc_max = new_am
        for lc in range(lc_min, lc_max+1):
            am_in[0] = la
            am_in[1] = lc
            if la == 0:
                a = 1
                k2max = la
                k3max = lc - 1
            else:
                a = 0
                k2max = lc
                k1max = la - 1
            foo = 5
            if (a == 0):
                foo = 4
            function_name = "build_%s0%s0"%(am_letter[la], am_letter[lc])
            code_name = function_name + ".jl"
            
            vrr_header.write('include("%s")\n'%code_name)
            
            code = open(lib+"/"+code_name, 'w')
            t1 = t2 = t3 = t4 = 0
            code.write('include("libints_header.jl")\n\n')
            code.write("# These machine-generated functions compute a quartet of (%cs|%cs) integrals \n\n"%(am_letter[la],am_letter[lc]))

            code.write("function %s!(Data::prim_data, vp::Vector{Float64}, vp_ind::Int, I0::Vector{Float64}, I0_ind::Int, I1::Vector{Float64}, I1_ind::Int, I2::Vector{Float64}, I2_ind::Int, I3::Vector{Float64}, I3_ind::Int, I4::Vector{Float64}, I4_ind::Int)\n"%function_name)
            code.write("\tlpoz = Data.poz\n")
            code.write("\tlpon = Data.pon\n")
            define_localv(a, foo, k1max, k2max, k3max, code)
            code.write("\n")

            FLOP_counter = 0

            for i in range(am_in[0]+1):
                am[0][0] = am_in[0] - i
                for j in range(i+1):
                    am[0][1] = i - j
                    am[0][2] = j
                    for k in range(am_in[1]+1):
                        am[1][0] = am_in[1] - k
                        for l in range(k+1):
                            am[1][1] = k - l
                            am[1][2] = l

                            if (am[a][2]): b = 2
                            if (am[a][1]): b = 1
                            if (am[a][0]): b = 0

                            am[a][b] = am[a][b] -1
                            am_in[a] = am_in[a] - 1
                            t2 = hash(am, am_in)
                            code.write("\tvp[vp_ind+1] = U%d%d*I0[I0_ind+%d] + U%d%d*I1[I1_ind+%d]"%(a*2, b, t2+1, foo, b, t2+1))
                            FLOP_counter += 3
                            if (am[a][b]):
                                am[a][b] = am[a][b] - 1
                                am_in[a] = am_in[a] - 1
                                t3 = hash(am, am_in)
                                code.write(" + (%s)*(I2[I2_ind+%d] - (%s)*I3[I3_ind+%d])"%((k1[am[a][b]] if a == 0 else k3[am[a][b]]), t3+1, k4[a], t3+1))
                                max1 = max(max1, am[a][b]+1)
                                am[a][b] = am[a][b] + 1
                                am_in[a] = am_in[a] + 1
                                FLOP_counter += 4
                            
                            
                            if (am[a^1][b]):
                               am[a^1][b] = am[a^1][b] - 1
                               am_in[a^1] = am_in[a^1] - 1
                               t4 = hash(am, am_in)
                               code.write(" + (%s)*I4[I4_ind+%d]"%(k2[am[a^1][b]], t4+1))
                               max2 = max(max2, am[a^1][b]+1)
                               am[a^1][b] = am[a^1][b] + 1
                               am_in[a^1] = am_in[a^1] + 1
                               FLOP_counter += 2
                            
                            code.write("\n")
                            code.write("\tvp_ind += 1\n")
                            am[a][b] = am[a][b] + 1
                            am_in[a] = am_in[a] + 1

                            t1 += 1

            code.write("end\n")
            code.write("# Total number of FLOPs = %d\n"%FLOP_counter)

            code.close()
            print("Done with %s"%code_name)
    vrr_header.close()

                               

def define_localv(a, foo, k1max, k2max, k3max, code):
    global k1, k2, k3
    for i in range(k2max):
        code.write("\t%s = %.1lf*Data.oo2zn\n"%(k2[i], float(i+1)))
    if (a == 0):
        for i in range(k1max):
            code.write("\t%s = %.1lf*Data.oo2z\n"%(k1[i], float(i+1)))
    else:
        for i in range(k3max):
            code.write("\t%s = %.1lf*Data.oo2n\n"%(k3[i], float(i+1)))
    
    code.write("\tU%d0 = Data.U[%d,1]\n"%(a*2, a*2+1))
    code.write("\tU%d1 = Data.U[%d,2]\n"%(a*2, a*2+1))
    code.write("\tU%d2 = Data.U[%d,3]\n"%(a*2, a*2+1))
    code.write("\tU%d0 = Data.U[%d,1]\n"%(foo, foo+1))
    code.write("\tU%d1 = Data.U[%d,2]\n"%(foo, foo+1))
    code.write("\tU%d2 = Data.U[%d,3]\n\n"%(foo, foo+1))


