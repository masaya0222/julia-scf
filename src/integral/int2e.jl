module Int2e
using JuliaSCF.Mole
using JuliaSCF.Integral.Int1e_nuc: Fn
using JuliaSCF.Lib.Libint

export get_v2e, get_v2e_single

function INT_NCART(am::Int)
    return div((am+1)*(am+2),2)
end

function cont_V2(basis_a::orb_detail, basis_b::orb_detail, basis_c::orb_detail, basis_d::orb_detail)
    w_fact = [1,1,3,15,105]
    am1 = basis_a.orb_l; am2 = basis_b.orb_l
    am3 = basis_c.orb_l; am4 = basis_d.orb_l
    am = am1 + am2 + am3 + am4
    a1s = basis_a.α_array; a2s = basis_b.α_array
    a3s = basis_c.α_array; a4s = basis_d.α_array
    da = basis_a.d_array; db = basis_b.d_array
    dc = basis_c.d_array; dd = basis_d.d_array
    A = basis_a.cord; B = basis_b.cord
    C = basis_c.cord; D = basis_d.cord
    AB2 = (A[1]-B[1])^2+(A[2]-B[2])^2+(A[3]-B[3])^2
    CD2 = (C[1]-D[1])^2+(C[2]-D[2])^2+(C[3]-D[3])^2

    nprim1 = length(a1s); nprim2 = length(a2s) 
    nprim3 = length(a3s); nprim4 = length(a4s)

    libint_ = Libint_t()
    init_libint(libint_, max(am1,am2,am3,am4))
    libint_.PrimQuartet = [prim_data() for i =1:(nprim1*nprim2*nprim3*nprim4)]
    for i = 1:(nprim1*nprim2*nprim3*nprim4)
        libint_.PrimQuartet[i].U = zeros(Float64, (6,3))    
    end
    fjt = zeros(Float64, am+1)
    libint_.AB = A - B
    libint_.CD = C - D
    size = INT_NCART(am1)*INT_NCART(am2)*INT_NCART(am3)*INT_NCART(am4)
    contV = Array{Vector{Float64},4}(undef, length(da),length(db),length(dc),length(dd))
    
    
    nprim = 0
    c1s = da[1];c2s = db[1]; c3s = dc[1]; c4s = dd[1]
    for p1 = 1:nprim1
        a1 = a1s[p1]
        c1 = c1s[p1]
        Na = (2*a1/pi)^(3/4)*sqrt(((4*a1)^am1)/w_fact[am1+1])
        c1 *= Na
        for p2 = 1:nprim2
            a2 = a2s[p2]
            c2 = c2s[p2]
            zeta = a1 + a2
            ooz = 1/zeta
            oo2z = 1/(2*zeta)
            P = (a1*A+a2*B)*ooz
            PA = P - A
            PB = P - B
            Nb = (2*a2/pi)^(3/4)*sqrt(((4*a2)^am2)/w_fact[am2+1])
            c2 *= Nb

            Sab = (pi*ooz)^(3/2)*exp(-a1*a2*ooz*AB2)*c1*c2
            for p3 = 1:nprim3
                a3 = a3s[p3]
                c3 = c3s[p3]
                Nc = (2*a3/pi)^(3/4)*sqrt(((4*a3)^am3)/w_fact[am3+1])
                c3 *= Nc
                for p4 = 1:nprim4
                    nprim += 1
                    a4 = a4s[p4]
                    c4 = c4s[p4]
                    Nd = (2*a4/pi)^(3/4)*sqrt(((4*a4)^am4)/w_fact[am4+1])
                    c4 *= Nd
                    nu = a3 + a4
                    oon = 1/nu
                    oo2n = 1/(2*nu)
                    oo2zn = 1/(2*(zeta+nu))
                    rho = (zeta*nu)/(zeta+nu)
                    oo2rho = 1/(2*rho)
                    Q = (a3*C+a4*D)*oon
                    QC = Q - C
                    QD = Q - D
                    PQ = P - Q
                    PQ2 = PQ[1]^2+PQ[2]^2+PQ[3]^2
                    W = (zeta*P+nu*Q)/(zeta+nu)
                    WP = W - P
                    WQ = W - Q
                    for i = 1:3
                        libint_.PrimQuartet[nprim].U[1,i] = PA[i]
                        libint_.PrimQuartet[nprim].U[2,i] = PB[i]
                        libint_.PrimQuartet[nprim].U[3,i] = QC[i]
                        libint_.PrimQuartet[nprim].U[4,i] = QD[i]
                        libint_.PrimQuartet[nprim].U[5,i] = WP[i]
                        libint_.PrimQuartet[nprim].U[6,i] = WQ[i]
                    end
                    libint_.PrimQuartet[nprim].oo2z = oo2z
                    libint_.PrimQuartet[nprim].oo2n = oo2n
                    libint_.PrimQuartet[nprim].oo2zn = oo2zn
                    libint_.PrimQuartet[nprim].poz = rho * ooz
                    libint_.PrimQuartet[nprim].pon = rho * oon
                    libint_.PrimQuartet[nprim].oo2p = oo2rho
                    T = rho*PQ2
                    
                    fjt = [Fn(n,T) for n = 0:am]
                    Scd = (pi*oon)^(3/2)*exp(-a3*a4*oon*CD2)*c3*c4
                    val = 2*sqrt(rho/pi)*Sab*Scd
                    libint_.PrimQuartet[nprim].F = fjt*val
                end
            end
        end
    end
    if am == 0
        temp = 0.0
        for i = 1:nprim
            temp += libint_.PrimQuartet[i].F[1]
        end
        return [temp]
    else
        target_ints = build_eri[am1+1,am2+1,am3+1,am4+1](libint_,nprim)
        return libint_.int_stack[begin+target_ints:target_ints+size]
    end
end

function V2e_lm(basis_a::orb_detail, basis_b::orb_detail, basis_c::orb_detail, basis_d::orb_detail, C_a::Array{Float64,4}, C_b::Array{Float64,4}, C_c::Array{Float64,4}, C_d::Array{Float64,4})
    V2e_abcd = cont_V2(basis_a, basis_b, basis_c, basis_d)
    
    I1 = basis_a.orb_l; J1 = basis_b.orb_l
    I2 = basis_c.orb_l; J2 = basis_d.orb_l
    a = basis_a.α_array; b = basis_b.α_array
    c = basis_c.α_array; d = basis_d.α_array
    da = basis_a.d_array; db = basis_b.d_array
    dc = basis_c.d_array; dd = basis_d.d_array
    
    max_l = max(I1,J1,I2,J2)
    fact = [factorial(i) for i = 0:2*max_l]
    index =[[[0]],  # for s
            [[2,1], [0]], # for p 0 = index[1][1][0], 1 = index[1][0][1], 2 = index[1][0][0]
            [[5,4,3],[2,1],[0]] # for d 0 = index[2][2][0], 1 = index[2][1][1], 2 = index[2][1][0], 3 = index[2][0][2], 4 = index[2][0][1], 5 = index[2][0][2]
            ]

    V2e_mamb = zeros(Float64, (2I1+1,2J1+1,2I2+1,2J2+1))
    
    for i=0:2I1, j=0:2J1, k=0:2I2, l=0:2J2
        ma = i-I1; mb = j-J1; mc = k-I2; md = l-J2
        ma_ = abs(ma); mb_ = abs(mb); mc_ = abs(mc); md_ = abs(md)
        Nma = (1/(2^ma_*fact[I1+1]))*sqrt(2*fact[I1+ma_+1]*fact[I1-ma_+1]/(2-(ma!=0)))
        Nmb = (1/(2^mb_*fact[J1+1]))*sqrt(2*fact[J1+mb_+1]*fact[J1-mb_+1]/(2-(mb!=0)))
        Nmc = (1/(2^mc_*fact[I2+1]))*sqrt(2*fact[I2+mc_+1]*fact[I2-mc_+1]/(2-(mc!=0)))
        Nmd = (1/(2^md_*fact[J2+1]))*sqrt(2*fact[J2+md_+1]*fact[J2-md_+1]/(2-(md!=0)))
        for ta=0:div(I1-ma_,2), tb=0:div(J1-mb_,2), tc=0:div(I2-mc_,2), td=0:div(J2-md_,2)
            for ua=0:ta, ub=0:tb, uc=0:tc, ud=0:td
                for va=0:div(ma_,2)
                    for vb=0:div(mb_,2)
                        for vc=0:div(mc_,2)
                            for vd=0:div(md_,2)
                                f = 1-2*((va+vb+vc+vd)%2)
                                va_ = va; vb_ = vb; vc_ = vc; vd_ = vd
                                if ma < 0
                                    va_ += 1/2
                                    if va_ > div(ma_-1,2)+1/2
                                        break
                                    end
                                end
                                if mb < 0
                                    vb_ += 1/2
                                    if vb_ > div(mb_-1,2)+1/2
                                        break
                                    end
                                end
                                if mc < 0
                                    vc_ += 1/2
                                    if vc_ > div(mc_-1,2)+1/2
                                        break
                                    end
                                end
                                if md < 0
                                    vd_ += 1/2
                                    if vd_ > div(md_-1,2)+1/2
                                        break
                                    end
                                end
                                pow_xa = Int(floor(2*ta+ma_-2*(ua+va_)))
                                pow_xb = Int(floor(2*tb+mb_-2*(ub+vb_)))
                                pow_xc = Int(floor(2*tc+mc_-2*(uc+vc_)))
                                pow_xd = Int(floor(2*td+md_-2*(ud+vd_)))
                                pow_ya = Int(floor(2*(ua+va_)))
                                pow_yb = Int(floor(2*(ub+vb_)))
                                pow_yc = Int(floor(2*(uc+vc_)))
                                pow_yd = Int(floor(2*(ud+vd_)))
                                tmp = 0
                                tmp += index[I1+1][pow_xa+1][pow_ya+1]*INT_NCART(J2)*INT_NCART(I2)*INT_NCART(J1)
                                tmp += index[J1+1][pow_xb+1][pow_yb+1]*INT_NCART(J2)*INT_NCART(I2)
                                tmp += index[I2+1][pow_xc+1][pow_yc+1]*INT_NCART(J2)
                                tmp += index[J2+1][pow_xd+1][pow_yd+1]

                                
                                V2e_mamb[i+1,j+1,k+1,l+1] += f*C_a[ma_+1,ta+1,ua+1,Int(floor(2*va_))+1] * 
                                                            C_b[mb_+1,tb+1,ub+1,Int(floor(2*vb_))+1] * 
                                                            C_c[mc_+1,tc+1,uc+1,Int(floor(2*vc_))+1] * 
                                                            C_d[md_+1,td+1,ud+1,Int(floor(2*vd_))+1] * 
                                                            V2e_abcd[tmp+1]
    
                            end
                        end
                    end
                end
            end
        end
        
        V2e_mamb[i+1,j+1,k+1,l+1] *= Nma*Nmb*Nmc*Nmd
        
    end
    return V2e_mamb
end

function get_v2e(mol::Molecule)
    basis = mol.basis_
    num = 0
    V2e = zeros(Float64, (mol.basis_num, mol.basis_num, mol.basis_num, mol.basis_num))
    basis_len = length(basis)
    check = zeros(Bool, (basis_len, basis_len, basis_len, basis_len))
    change = [[0],[1,2,0],[0,1,2,3,4]]
    Tasks = Vector{Vector{Int}}(undef, basis_len^4)
    ind_i = 1
    I1 = 0;J1 = 0;I2 = 0; J2 = 0
    for i = 1:basis_len
        ind_j = 1
        for j = 1:basis_len
            ind_k = 1
            for k = 1:basis_len
                ind_l = 1
                for l = 1:basis_len
                    I1 = basis[i].orb_l; J1 = basis[j].orb_l
                    I2 = basis[k].orb_l; J2 = basis[l].orb_l
                    if !check[i,j,k,l] && I1 >= J1 && I2 >= J2 && I2+J2 >= I1+J1
                        num += 1
                        check[i,j,k,l] = check[i,j,l,k] = check[j,i,k,l] = check[j,i,l,k] = true
                        check[k,l,i,j] = check[k,l,j,i] = check[l,k,i,j] = check[l,k,j,i] = true
                        Tasks[num] = [i,j,k,l,I1,J1,I2,J2]                   
                    end
                    ind_l += 2J2+1
                end
                ind_k += 2I2+1
            end
            ind_j += 2J1+1
        end
        ind_i += 2I1+1
    end
    results = Vector{Array{Float64, 4}}(undef, num)
    Threads.@threads for ind = 1:num
        results[ind] = V2e_lm(basis[Tasks[ind][1]],basis[Tasks[ind][2]],basis[Tasks[ind][3]],basis[Tasks[ind][4]],mol.rotate_coef[Tasks[ind][5]+1],mol.rotate_coef[Tasks[ind][6]+1],mol.rotate_coef[Tasks[ind][7]+1],mol.rotate_coef[Tasks[ind][8]+1])
    end
    
    check = zeros(Bool, (basis_len, basis_len, basis_len, basis_len))
    num = 0
    ind_i = 1
    I1 = 0;J1 = 0;I2 = 0; J2 = 0
    for i = 1:basis_len
        ind_j = 1
        for j = 1:basis_len
            ind_k = 1
            for k = 1:basis_len
                ind_l = 1
                for l = 1:basis_len
                    I1 = basis[i].orb_l; J1 = basis[j].orb_l
                    I2 = basis[k].orb_l; J2 = basis[l].orb_l
                    if !check[i,j,k,l] && I1 >= J1 && I2 >= J2 && I2+J2 >= I1+J1
                        num += 1
                        check[i,j,k,l] = check[i,j,l,k] = check[j,i,k,l] = check[j,i,l,k] = true
                        check[k,l,i,j] = check[k,l,j,i] = check[l,k,i,j] = check[l,k,j,i] = true
                        
                        V2elm = results[num]
                        
                        for ma=0:2I1, mb=0:2J1, mc=0:2I2, md=0:2J2
                            ans = V2elm[ma+1,mb+1,mc+1,md+1]
                            ind_1 = ind_i+change[I1+1][ma+1]
                            ind_2 = ind_j+change[J1+1][mb+1]
                            ind_3 = ind_k+change[I2+1][mc+1]
                            ind_4 = ind_l+change[J2+1][md+1]
                            V2e[ind_1,ind_2,ind_3,ind_4] = ans

                            V2e[ind_1,ind_2,ind_4,ind_3] = ans
                            V2e[ind_2,ind_1,ind_3,ind_4] = ans
                            V2e[ind_2,ind_1,ind_4,ind_3] = ans

                            V2e[ind_3,ind_4,ind_1,ind_2] = ans
                            V2e[ind_3,ind_4,ind_2,ind_1] = ans
                            V2e[ind_4,ind_3,ind_1,ind_2] = ans
                            V2e[ind_4,ind_3,ind_2,ind_1] = ans
                        end                
                    end
                    ind_l += 2J2+1
                end
                ind_k += 2I2+1
            end
            ind_j += 2J1+1
        end
        ind_i += 2I1+1
    end
    return V2e
end

function get_v2e_single(mol::Molecule)
    basis = mol.basis_
    V2e = zeros(Float64, (mol.basis_num, mol.basis_num, mol.basis_num, mol.basis_num))
    basis_len = length(basis)
    check = zeros(Bool, (basis_len, basis_len, basis_len, basis_len))
    change = [[0],[1,2,0],[0,1,2,3,4]]
    ind_i = 1
    I1 = 0;J1 = 0;I2 = 0; J2 = 0
    for i = 1:basis_len
        ind_j = 1
        for j = 1:basis_len
            ind_k = 1
            for k = 1:basis_len
                ind_l = 1
                for l = 1:basis_len
                    I1 = basis[i].orb_l; J1 = basis[j].orb_l
                    I2 = basis[k].orb_l; J2 = basis[l].orb_l
                    if !check[i,j,k,l] && I1 >= J1 && I2 >= J2 && I2+J2 >= I1+J1
                        check[i,j,k,l] = check[i,j,l,k] = check[j,i,k,l] = check[j,i,l,k] = true
                        check[k,l,i,j] = check[k,l,j,i] = check[l,k,i,j] = check[l,k,j,i] = true
                        V2elm = V2e_lm(basis[i],basis[j],basis[k],basis[l],mol.rotate_coef[I1+1],mol.rotate_coef[J1+1],mol.rotate_coef[I2+1],mol.rotate_coef[J2+1])
                        
                        for ma=0:2I1, mb=0:2J1, mc=0:2I2, md=0:2J2
                            ans = V2elm[ma+1,mb+1,mc+1,md+1]
                            ind_1 = ind_i+change[I1+1][ma+1]
                            ind_2 = ind_j+change[J1+1][mb+1]
                            ind_3 = ind_k+change[I2+1][mc+1]
                            ind_4 = ind_l+change[J2+1][md+1]
                            V2e[ind_1,ind_2,ind_3,ind_4] = ans

                            V2e[ind_1,ind_2,ind_4,ind_3] = ans
                            V2e[ind_2,ind_1,ind_3,ind_4] = ans
                            V2e[ind_2,ind_1,ind_4,ind_3] = ans

                            V2e[ind_3,ind_4,ind_1,ind_2] = ans
                            V2e[ind_3,ind_4,ind_2,ind_1] = ans
                            V2e[ind_4,ind_3,ind_1,ind_2] = ans
                            V2e[ind_4,ind_3,ind_2,ind_1] = ans
                        
                        end
                        
                    end
                    ind_l += 2J2+1
                end
                ind_k += 2I2+1
            end
            ind_j += 2J1+1
        end
        ind_i += 2I1+1
    end
    return V2e
end
    
end