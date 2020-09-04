module Int2e
using JuliaSCF.Mole
using JuliaSCF.Integral.Int1e_nuc: R_tuv, Fn, E_ij_t

function g_abcd(I1::Int, J1::Int, I2::Int, J2::Int, Ra::Vector{Float64}, Rb::Vector{Float64}, Rc::Vector{Float64}, Rd::Vector{Float64}, ai::Float64, bi::Float64, ci::Float64, di::Float64)
    p = ai + bi; q = ci + di
    α = p*q/(p+q)
    Rp = [(ai*Ra[i]+bi*Rb[i])/p for i = 1:3]
    Rq = [(ci*Rc[i]+di*Rd[i])/q for i = 1:3]

    Rtuv = R_tuv(I1+J1+I2+J2, Rp, Rq, α)
    gabcd = zeros(Float64, (I1+1,J1+1,I1+1,J1+1,I2+1,J2+1,I2+1,J2+1))

    Ejit1 = E_ij_t(I1, J1, Ra[1], Rb[1], ai, bi)
    Eklu1 = E_ij_t(I1, J1, Ra[2], Rb[2], ai, bi)
    Emnv1 = E_ij_t(I1, J1, Ra[3], Rb[3], ai, bi)

    Ejit2 = E_ij_t(I2, J2, Rc[1], Rd[1], ci, di)
    Eklu2 = E_ij_t(I2, J2, Rc[2], Rd[2], ci, di)
    Emnv2 = E_ij_t(I2, J2, Rc[3], Rd[3], ci, di)

    for i1=0:I1, j1=0:J1, k1=0:I1-i1, l1=0:J1-j1
        m1 = I1-i1-k1; n1 = J1-j1-l1
        for i2=0:I2, j2=0:J2, k2=0:I2-i2, l2=0:J2-j2
            m2 = I2-i2-k2; n2 = J2-j2-l2
            ans1 = 0.0
            for t1=0:i1+j1, u1=0:k1+l1, v1=0:m1+n1
                Etuv = Ejit1[i1+1,j1+1,t1+1]*Eklu1[k1+1,l1+1,u1+1]*Emnv1[m1+1,n1+1,v1+1]
                ans2 = 0.0
                for t2=0:i2+j2, u2=0:k2+l2, v2=0:m2+n2
                    f = 1.0-2*((t2+u2+v2)%2)
                    ans2+=f*Ejit2[i2+1,j2+1,t2+1]*Eklu2[k2+1,l2+1,u2+1]*Emnv2[m2+1,n2+1,v2+1]*Rtuv[t1+t2+1,u1+u2+1,v1+v2+1]
                end
                ans1 += ans2*Etuv
            end 
            gabcd[i1+1,j1+1,k1+1,l1+1,i2+1,j2+1,k2+1,l2+1] = ans1*(2*pi^(5/2))/(p*q*sqrt(p+q))
        end
     end
    return gabcd
end

function cont_V2(basis_a::orb_detail, basis_b::orb_detail, basis_c::orb_detail, basis_d::orb_detail)
    w_fact = [1,1,3,15,105]
    I1 = basis_a.orb_l; J1 = basis_b.orb_l
    I2 = basis_c.orb_l; J2 = basis_d.orb_l
    a = basis_a.α_array; b = basis_b.α_array
    c = basis_c.α_array; d = basis_d.α_array
    da = basis_a.d_array; db = basis_b.d_array
    dc = basis_c.d_array; dd = basis_d.d_array
    Ra = basis_a.cord; Rb = basis_b.cord
    Rc = basis_c.cord; Rd = basis_d.cord

    gabcd = [[[[g_abcd(I1,J1,I2,J2,Ra,Rb,Rc,Rd,ai,bi,ci,di) for di=d] for ci=c] for bi=b] for ai=a]
    contV = zeros(Float64, (length(da),length(db),length(dc),length(dd),I1+1,J1+1,I1+1,J1+1,I2+1,J2+1,I2+1,J2+1))
    for ind_a=1:length(da), ind_b=1:length(db), ind_c=1:length(dc), ind_d=1:length(dd) 
        for i1=0:I1,j1=0:J1,k1=0:I1-i1,l1=0:J1-j1
            m1 = I1-i1-k1; n1 = J1-j1-l1
            for i2=0:I2,j2=0:J2,k2=0:I2-i2,l2=0:J2-j2
                m2 = I2-i2-k2; n2 = J2-j2-l2
                ans = 0.0
                for p1=1:length(a), q1=1:length(b), p2=1:length(c), q2=1:length(d)
                    Na = (2*a[p1]/pi)^(3/4)*sqrt(((4*a[p1])^I1)/w_fact[I1+1])
                    Nb = (2*b[q1]/pi)^(3/4)*sqrt(((4*b[q1])^J1)/w_fact[J1+1])
                    Nc = (2*c[p2]/pi)^(3/4)*sqrt(((4*c[p2])^I2)/w_fact[I2+1])
                    Nd = (2*d[q2]/pi)^(3/4)*sqrt(((4*d[q2])^J2)/w_fact[J2+1])
                    ans += da[ind_a][p1]*db[ind_b][q1]*dc[ind_c][p2]*dd[ind_d][q2]*gabcd[p1][q1][p2][q2][i1+1,j1+1,k1+1,l1+1,i2+1,j2+1,k2+1,l2+1]*Na*Nb*Nc*Nd
                end
                contV[ind_a,ind_b,ind_c,ind_d,i1+1,j1+1,k1+1,l1+1,i2+1,j2+1,k2+1,l2+1] = ans
            end
        end 
    end
    return contV
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

    V2e_mamb = zeros(Float64, (length(da),length(db),length(dc),length(dd),2I1+1,2J1+1,2I2+1,2J2+1))

    for i=0:2I1, j=0:2J1, k=0:2I2, l=0:2J2
        ma = i-I1; mb = j-J1; mc = k-I2; md = l-J2
        ma_ = abs(ma); mb_ = abs(mb); mc_ = abs(mc); md_ = abs(md)
        Nma = (1/(2^ma_*fact[I1+1]))*sqrt(2*fact[I1+ma_+1]*fact[I1-ma_+1]/(2-(ma!=0)))
        Nmb = (1/(2^mb_*fact[J1+1]))*sqrt(2*fact[J1+mb_+1]*fact[J1-mb_+1]/(2-(mb!=0)))
        Nmc = (1/(2^mc_*fact[I2+1]))*sqrt(2*fact[I2+mc_+1]*fact[I2-mc_+1]/(2-(mc!=0)))
        Nmd = (1/(2^md_*fact[J2+1]))*sqrt(2*fact[J2+md_+1]*fact[J2-md_+1]/(2-(md!=0)))
        for ta=0:div(I1-ma_,2), tb=0:div(J1-mb_,2), tc=0:div(I2-mc_,2), td=0:div(J2-md_,2)
            for ua=0:ta, ub=0:tb, uc=0:tc, ud=0:td
                for va=0:div(ma_,2), vb=0:div(mb_,2), vc=0:div(mc_,2), vd=0:div(md_,2)
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

                    for ind_a=1:length(da), ind_b=1:length(db) ,ind_c=1:length(dc), ind_d=1:length(dd)
                        V2e_mamb[ind_a,ind_b,ind_c,ind_d,i+1,j+1,k+1,l+1] += f*C_a[ma_+1,ta+1,ua+1,Int(floor(2*va_))+1] * 
                                                                               C_b[mb_+1,tb+1,ub+1,Int(floor(2*vb_))+1] * 
                                                                               C_c[mc_+1,tc+1,uc+1,Int(floor(2*vc_))+1] * 
                                                                               C_d[md_+1,td+1,ud+1,Int(floor(2*vd_))+1] * 
                                                                               V2e_abcd[ind_a,ind_b,ind_c,ind_d,pow_xa+1,pow_xb+1,pow_ya+1,pow_yb+1,pow_xc+1,pow_xd+1,pow_yc+1,pow_yd+1]
                    end
                end
            end
        end
        for ind_a=1:length(da),ind_b=1:length(db),ind_c=1:length(dc),ind_d=1:length(dd)
            V2e_mamb[ind_a,ind_b,ind_c,ind_d,i+1,j+1,k+1,l+1] *= Nma*Nmb*Nmc*Nmd
        end
    end
    return V2e_mamb
end

function get_v2e(mol::Molecule)
    basis = mol.basis
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
                    if !check[i,j,k,l]
                        check[i,j,k,l] = check[i,j,l,k] = check[j,i,k,l] = check[j,i,l,k] = true
                        check[k,l,i,j] = check[k,l,j,i] = check[l,k,i,j] = check[l,k,j,i] = true
                        
                        V2elm = V2e_lm(basis[i],basis[j],basis[k],basis[l],mol.rotate_coef[I1+1],mol.rotate_coef[J1+1],mol.rotate_coef[I2+1],mol.rotate_coef[J2+1])

                        for ind_a=0:length(basis[i].d_array)-1, ind_b=0:length(basis[j].d_array)-1, ind_c=0:length(basis[k].d_array)-1, ind_d=0:length(basis[l].d_array)-1 
                            for ma=0:2I1, mb=0:2J1, mc=0:2I2, md=0:2J2
                                ans = V2elm[ind_a+1,ind_b+1,ind_c+1,ind_d+1,ma+1,mb+1,mc+1,md+1]
                                ind_1 = ind_i+ind_a*(2I1+1)+change[I1+1][ma+1]
                                ind_2 = ind_j+ind_b*(2J1+1)+change[J1+1][mb+1]
                                ind_3 = ind_k+ind_c*(2I2+1)+change[I2+1][mc+1]
                                ind_4 = ind_l+ind_d*(2J2+1)+change[J2+1][md+1]
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
                    end
                    ind_l += (2J2+1)*length(basis[l].d_array)
                end
                ind_k += (2I2+1)*length(basis[k].d_array)
            end
            ind_j += (2J1+1)*length(basis[j].d_array)
        end
        ind_i += (2I1+1)*length(basis[i].d_array)
    end
    return V2e
end
#m = Molecule([atom("Br",[0.0,0.0,-0.7]),atom("Br",[0.0,0.0,0.7])], "sto3g")
#m = Molecule([atom("C",[0.0,0.0,-0.7]),atom("C",[0.0,0.0,0.7]),atom("C",[0.0,0.7,0.0]),atom("C",[0.0,-0.7,0.0]),atom("C",[0.7,0.0,0.0])], "ccpvdz")
#m = Molecule([atom("C",[0.0,0.0,-0.7]),atom("C",[0.0,0.0,0.7]),atom("C",[0.0,-0.7,0.0]),atom("C",[0.0,0.7,0.0]),atom("C",[-0.7,0.0,0.0]),atom("C",[0.7,0.0,0.0])], "sto3g")

#@time v = get_v2e(m)
#@show size(v)
#=
using PyCall
@pyimport pyscf
X = 0.52918
mol = pyscf.gto.Mole()
mol.build(atom="C 0 0 $(-0.7*X); C 0 0 $(0.7*X);C 0 $(-0.7*X) 0; C 0 $(0.7*X) 0;C $(0.7*X) 0 0; C $(-0.7*X) 0 0; C $(-1.4*X) 0 0;C $(-0.7*X) 0 0", basis="ccpvdz")
@time v = mol.intor("int2e")
=#
#@show v

end