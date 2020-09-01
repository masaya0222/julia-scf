module Int1e_nuc
using JuliaSCF.Mole
using SpecialFunctions

export get_v1e

function V_ijklmn(I::Int, J::Int, Ra::Vector{Float64}, Rb::Vector{Float64}, Rc_list::Vector{Vector{Float64}}, Zc_list::Vector{Int}, ai::Float64, bi::Float64)
    p = ai + bi
    Rp = [(ai*Ra[i] + bi*Rb[i])/p for i = 1:3]
    Rtuv = zeros(Float64, (I+J+1,I+J+1, I+J+1))
    for (Rc, Zc) = zip(Rc_list, Zc_list)
        r = R_tuv(I+J, Rp, Rc, p)
        for t = 0:I+J, u = 0:I+J-t, v = 0:I+J-t-u
            Rtuv[t+1,u+1,v+1] += Zc * r[t+1,u+1,v+1]
        end
    end

    Eijt = E_ij_t(I, J, Ra[1], Rb[1], ai, bi)
    Eklu = E_ij_t(I, J, Ra[2], Rb[2], ai, bi)
    Emnv = E_ij_t(I, J, Ra[3], Rb[3], ai, bi)
    Vijklmn = zeros(Float64, (I+1, J+1, I+1, J+1))
    for i = 0:I, j = 0:J, k = 0:I-i, l = 0:J-j
        m = I-i-k
        n = J-j-l
        for t = 0:i+j, u = 0:k+l, v = 0:m+n
            Vijklmn[i+1,j+1,k+1,l+1] += Eijt[i+1,j+1,t+1]*Eklu[k+1,l+1,u+1]*Emnv[m+1,n+1,v+1]*Rtuv[t+1,u+1,v+1]
        end
        Vijklmn[i+1,j+1,k+1,l+1] *= -2*pi/p
    end 
    return Vijklmn
end

function R_tuv(IJ::Int64, Rp::Vector{Float64}, Rc::Vector{Float64}, p::Float64)
    Xpc = Rp[1] - Rc[1]
    Ypc = Rp[2] - Rc[2]
    Zpc = Rp[3] - Rc[3]
    Rpc_2 = Xpc^2 + Ypc^2 + Zpc^2
    R = [zeros(Float64, (IJ-N+1, IJ-N+1, IJ-N+1)) for N = 0:IJ]

    for n = 0:IJ
        R[n+1][1,1,1] += (-2*p)^n *Fn(n, p*Rpc_2)
    end
    for n = IJ-1:-1:0, t = 0:IJ-n, u = 0:IJ-n-t, v = 0:IJ-n-t-u
        if t >= 1
            if t >= 2
                R[n+1][t+1,u+1,v+1] += (t-1) * R[n+2][t-1,u+1,v+1]
            end
            R[n+1][t+1,u+1,v+1] += Xpc*R[n+2][t,u+1,v+1]
        elseif u >= 1
            if u >=2
                R[n+1][t+1,u+1,v+1] += (u-1)*R[n+2][t+1,u-1,v+1]
            end
            R[n+1][t+1,u+1,v+1] += Ypc*R[n+2][t+1,u,v+1]
        elseif v >= 1
            if v >= 2
                R[n+1][t+1,u+1,v+1] += (v-1)*R[n+2][t+1,u+1,v-1]
            end
            R[n+1][t+1,u+1,v+1] += Zpc*R[n+2][t+1,u+1,v]
        end
    end
    return R[1]
end

function Fn(n::Int, x::Float64)
    if x <= 0.1
        result = 1/(2*n+1)
        result_k = 1
        result_x = 1
        for k = 1:6
            result_x *= (-x)
            result_k *= k
            result += result_x/(result_k*(2*n+2*k+1))
        end
        return result
    end
    value_a = SpecialFunctions.gamma(n + 1/2)
    value_b = SpecialFunctions.gamma_inc(n + 1/2, x,0)[1]
    value_c = 2*x^(n + 1/2)
    return value_a*value_b/value_c
end

function E_ij_t(I::Int, J::Int, Ax::Float64, Bx::Float64, ai::Float64, bi::Float64)
    p = ai + bi
    mu = ai*bi/p
    Xab = Ax - Bx
    E = zeros(Float64, (I+1,J+1,I+J+1))
    E[1,1,1] = exp(-mu*Xab^2)
    for j = 0:J, i = 0:I
        if i == 0 && j ==0
            continue
        end
        if i == 0
            for t = 0:i+j
                if 0 <= t - 1
                    E[i+1,j+1,t+1] += (1/(2p))*E[i+1,j,t]
                end
                if t <= i + j -1
                    E[i+1,j+1,t+1] += (ai/p)*Xab*E[i+1,j,t+1]
                end
                if t + 1 <= i + j -1
                    E[i+1,j+1,t+1] += (t+1)*E[i+1,j,t+2]
                end
            end
        else
            for t = 0:i+j
                if 0 <= t - 1
                    E[i+1,j+1,t+1] += (1/(2p))*E[i,j+1,t]
                end
                if t <= i + j - 1
                    E[i+1,j+1,t+1] += (-bi/p)*Xab*E[i,j+1,t+1]
                end
                if t + 1 <= i + j - 1
                    E[i+1,j+1,t+1] += (t+1)*E[i,j+1,t+2]
                end
            end
        end
    end 
    return E
end

function cont_V1e(basis_a::orb_detail, basis_b::orb_detail, Rc_list::Vector{Vector{Float64}}, Zc_list::Vector{Int})
    w_fact = [1,1,3,15,105]
    I = basis_a.orb_l;J = basis_b.orb_l
    a = basis_a.α_array; b = basis_b.α_array
    da = basis_a.d_array; db = basis_b.d_array
    Vijklmn = [[V_ijklmn(I,J,basis_a.cord,basis_b.cord,Rc_list, Zc_list,ai,bi) for bi = b] for ai = a]

    contV = zeros(Float64,(length(da), length(db),I+1,J+1,I+1,J+1))
    for ind_a = 1:length(da), ind_b = 1:length(db)
        for i = 0:I, j = 0:J, k = 0:I-i, l = 0:J-j
            m = I-i-k
            n = J-j-l
            for p = 1:length(a), q = 1:length(b)
                ans = da[ind_a][p]*db[ind_b][q]*Vijklmn[p][q][i+1,j+1,k+1,l+1]
                Na = (2*a[p]/pi)^(3/4)*sqrt(((4*a[p])^I)/w_fact[I+1])
                Nb = (2*b[q]/pi)^(3/4)*sqrt(((4*b[q])^J)/w_fact[J+1])
                ans *= Na * Nb

                contV[ind_a,ind_b,i+1,j+1,k+1,l+1] += ans
            end
        end
    end
    return contV
end

function V1e_lm(basis_a::orb_detail, basis_b::orb_detail, Rc_list::Vector{Vector{Float64}}, Zc_list::Vector{Int}, C_a::Array{Float64,4}, C_b::Array{Float64,4})
    V1e_ab = cont_V1e(basis_a, basis_b, Rc_list, Zc_list)

    la = basis_a.orb_l; lb = basis_b.orb_l
    a = basis_a.α_array; b = basis_b.α_array
    da = basis_a.d_array; db = basis_b.d_array

    max_l = max(basis_a.orb_l, basis_b.orb_l)
    fact = [factorial(i) for i = 0:2*max_l]

    V1e_mamb = zeros(Float64, (length(da), length(db), 2*la+1, 2*lb+1))
    for i = 0:2la, j = 0:2lb
        ma = i - la; mb = j - lb
        ma_ = abs(ma); mb_ = abs(mb)
        Nma = (1/(2^ma_*fact[la+1]))*sqrt(2*fact[la+ma_+1]*fact[la-ma_+1]/(2-(ma!=0)))
        Nmb = (1/(2^mb_*fact[lb+1]))*sqrt(2*fact[lb+mb_+1]*fact[lb-mb_+1]/(2-(mb!=0)))
        for ta = 0:div(la-ma_,2), tb = 0:div(lb-mb_,2), ua = 0:ta, ub = 0:tb
            flag = false
            for va = 0:div(ma_,2), vb = 0:div(mb_,2)
                f = -2*((va+vb)%2)+1
                va_ = va
                vb_ = vb
                if ma < 0
                    va_ += 1/2
                    if va_ > div(ma_-1,2)+1/2
                        flag = true
                        break
                    end
                end
                if mb < 0
                    vb_ += 1/2
                    if vb_ > div(mb_-1,2)+1/2
                        flag = true
                        break
                    end
                end
                pow_xa = Int(floor(2*ta+ma_-2*(ua+va_)))
                pow_xb = Int(floor(2*tb+mb_-2*(ub+vb_)))
                pow_ya = Int(floor(2*(ua+va_)))
                pow_yb = Int(floor(2*(ub+vb_)))
                for ind_a = 1:length(da), ind_b = 1:length(db)
                    V1e_mamb[ind_a, ind_b,i+1, j+1] += f*C_a[ma_+1,ta+1,ua+1,Int(floor(2*va_))+1]*C_b[mb_+1,tb+1,ub+1,Int(floor(2*vb_))+1]*V1e_ab[ind_a, ind_b, pow_xa+1,pow_xb+1,pow_ya+1,pow_yb+1]
                end

                if flag
                    break
                end
            end
        end
        for ind_a = 1:length(da), ind_b = 1:length(db)
            V1e_mamb[ind_a, ind_b,i+1,j+1] *= Nma*Nmb
        end 
    end
    return V1e_mamb
end

function get_v1e(mol::Molecule)
    basis = mol.basis
    Rc_list = Vector{Vector{Float64}}()
    Zc_list = Vector{Int}()
    for atom = mol.atoms
        push!(Rc_list, atom.cord)
        push!(Zc_list, atom.atomic_num)
    end
    V1e = zeros(Float64, (mol.basis_num, mol.basis_num))
    basis_len = length(basis)
    check = zeros(Bool, (basis_len, basis_len))
    ind_i = 1
    change = [[0], [1, 2, 0], [0, 1, 2, 3, 4]]
    la = 0; lb = 0
    for i = 1:basis_len
        ind_j = 1
        for j = 1:basis_len
            la = basis[i].orb_l
            lb = basis[j].orb_l
            
            if !check[i,j]
                check[i,j] = check[j,i] = true
                V1elm = V1e_lm(basis[i], basis[j], Rc_list, Zc_list, mol.rotate_coef[la+1], mol.rotate_coef[lb+1])
                
                for ind_a = 0:length(basis[i].d_array)-1, ind_b = 0:length(basis[j].d_array)-1
                    for k = 0:2*la, l = 0:2*lb
                        if i == j
                            V1e[ind_i + ind_a*(2la+1) + change[la+1][k+1], ind_j + ind_b*(2lb+1) + change[lb+1][l+1]] = V1elm[ind_a+1,ind_b+1,k+1,l+1]
                        else
                            V1e[ind_i + ind_a*(2la+1) + change[la+1][k+1], ind_j + ind_b*(2lb+1) + change[lb+1][l+1]] = V1elm[ind_a+1,ind_b+1,k+1,l+1]
                            V1e[ind_j + ind_b*(2lb+1) + change[lb+1][l+1], ind_i + ind_a*(2la+1) + change[la+1][k+1]] = V1elm[ind_a+1,ind_b+1,k+1,l+1]
                        end
                    end
                end
            end
            ind_j += (2lb+1)*length(basis[j].d_array)
        end
        ind_i += (2la+1)*length(basis[i].d_array)
    end
    return V1e
end

end