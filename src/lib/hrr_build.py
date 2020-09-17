import constant
from constant import hash, io
import numpy as np
def hrr_build():
    new_am = constant.new_am
    am_letter = constant.am_letter
    lib = constant.lib
    am_in = [0,0]
    am = np.zeros((2,3),dtype=int)

    hrr_header = open(lib+"/"+"hrr_header.jl",'w')

    for lc in range(new_am+1):
        ld_max = lc if lc//2 + 1 > lc else lc//2 + 1
        for ld in range(1,ld_max+1):
            am_in[0] = lc-ld
            am_in[1] = ld

            split = 0

            function_name = "hrr3_build_"+am_letter[am_in[0]]+am_letter[am_in[1]]
            code_name = function_name+".jl"

            hrr_header.write('include("%s")\n'%code_name)

            code = open(lib+"/"+code_name, 'w')
            code.write("# These machine-generated functions compute a quartet of |%s%s) integrals \n\n" %(am_letter[am_in[0]],am_letter[am_in[1]]))
            code.write("function %s!(CD::Vector{Float64},vp::Vector{Float64}, vp_ind::Int, I0::Vector{Float64}, I0_ind::Int, I1::Vector{Float64}, I1_ind::Int, ab_num::Int)\n"%(function_name))

            #cuur_count = 0

            code.write("\tCD0 = CD[1]\n")
            code.write("\tCD1 = CD[2]\n")
            code.write("\tCD2 = CD[3]\n")
            nl = (am_in[1]*(am_in[1]+1))//2
            i0_step = (am_in[0]+2)*(am_in[0]+3)*nl//2
            i1_step = (am_in[0]+1)*(am_in[0]+2)*nl//2
            code.write("\tfor ab = 1:ab_num\n")

            FLOP_counter = 0

            for p in range(am_in[0]+1):
                am[0][0] = am_in[0]-p
                for q in range(p+1):
                    am[0][1] = p - q
                    am[0][2] = q
                    for r in range(am_in[1]+1):
                        am[1][0] = am_in[1] - r
                        for s in range(r+1):
                            am[1][1] = r-s
                            am[1][2] = s

                            if am[1][0]:
                                xyz = 0
                            elif am[1][1]:
                                xyz = 1
                            else:
                                xyz = 2
                            am[0][xyz] += 1
                            am_in[0] += 1
                            am[1][xyz] -= 1
                            am_in[1] -= 1
                            t0 = hash(am, am_in)
                            am[0][xyz] -= 1
                            am_in[0] -= 1
                            t1 = hash(am, am_in)
                            am[1][xyz] += 1
                            am_in[1] += 1

                            code.write("\t\tvp[vp_ind+1] = I0[I0_ind+%d] + CD%d*I1[I1_ind+%d]\n"
                            %(t0+1, xyz, t1+1))
                            code.write("\t\tvp_ind += 1\n")
                            FLOP_counter+=2
            if split == 0:
                code.write("\t\tI0_ind += %d\n"%(i0_step))
                code.write("\t\tI1_ind += %d\n"%(i1_step))
                code.write("\tend\nend\n")
            code.write("# Total number of FLOPs = %d * ab_num \n"%FLOP_counter)
            code.close()
            print("Done with %s"%code_name)

            la = lc - ld
            lb = ld
            am_in[0] = la
            am_in[1] = lb

            function_name = "hrr1_build_"+am_letter[am_in[0]]+am_letter[am_in[1]]
            code_name = function_name+".jl"

            hrr_header.write('include("%s")\n'%code_name)

            code = open(lib+"/"+code_name,'w')
            code.write("# These machine-generated functions compute a quartet of |%s%s) integrals \n\n" %(am_letter[am_in[0]],am_letter[am_in[1]]))
            code.write("function %s!(AB::Vector{Float64},vp::Vector{Float64}, vp_ind::Int, I0::Vector{Float64}, I0_ind::Int, I1::Vector{Float64}, I1_ind::Int, cd_num::Int)\n"%(function_name))

            code.write("\tAB0 = AB[1]\n")
            code.write("\tAB1 = AB[2]\n")
            code.write("\tAB2 = AB[3]\n")
            code.write("\n")

            #nj = (lb*(lb+1))//2
            FLOP_counter = 0
            for p in range(am_in[0]+1):
                am[0][0] = am_in[0] - p
                for q in range(p+1):
                    am[0][1] = p - q
                    am[0][2] = q

                    for r in range(am_in[1]+1):
                        am[1][0] = am_in[1] - r
                        for s in range(r+1):
                            am[1][1] = r - s
                            am[1][2] = s
                            if am[1][0]:
                                xyz = 0
                            elif am[1][1]:
                                xyz = 1
                            else:
                                xyz = 2
                            
                            am[0][xyz] += 1
                            am_in[0] += 1
                            am[1][xyz] -= 1
                            am_in[1] -= 1
                            t0 = hash(am, am_in)
                            am[0][xyz] -= 1
                            am_in[0] -= 1
                            t1 = hash(am, am_in)
                            am[1][xyz] += 1
                            am_in[1] += 1

                            if t0:
                                code.write("\ti0_ind = I0_ind + %d*cd_num\n"%t0)
                            else:
                                code.write("\ti0_ind = I0_ind\n")
                            if t1:
                                code.write("\ti1_ind = I1_ind + %d*cd_num\n"%t1)
                            else:
                                code.write("\ti1_ind = I1_ind\n")
                            
                            code.write("\tfor cd = 1:cd_num\n")
                            code.write("\t\tvp[vp_ind+1] = I0[i0_ind+1] + AB%d*I1[i1_ind+1]\n"%xyz)
                            code.write("\t\tvp_ind += 1; i0_ind += 1; i1_ind += 1\n")
                            code.write("\tend\n")
                            FLOP_counter += 2

            code.write("end\n")
            code.write("# Total number of FLOPs = %d * cd_num \n" %FLOP_counter)
            code.close
            print("Done with %s"%code_name)
    hrr_header.close()

        