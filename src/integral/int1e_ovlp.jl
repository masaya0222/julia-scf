module Int1e_ovlp
using JuliaSCF.Mole

export get_ovlp

function S_ij(I::Int, J::Int, Ax::Float64, Bx::Float64, ai::Float64, bi::Float64)
    p = ai + bi
    mu = ai*bi/p
    Xab = Ax - Bx
    S = zeros(Float64, (I+1, J+1))
    S[1,1] = sqrt(pi/p)*exp(-mu*Xab^2)
    for i = 1:I
        S[i+1,1] = (-bi/p)*Xab*S[i,1]
        if i !=1
            S[i+1,1] += (1/(2*p))*((i-1)*S[i-1][1])
        end
    end
    for j = 1:J ,i = 1:I+1
        S[i, j+1] = (ai/p)*Xab*S[i,j]
        if i != 1
            S[i, j+1] += (1/(2*p))*((i-1)*S[i-1,j])
        end
        if j != 1
            S[i, j+1] += (1/(2*p))*((j-1)*S[i, j-1])
        end 
    end
    return S
end

function cont_Sij(basis_a::orb_detail, basis_b::orb_detail)
    w_fact = [1,1,3,15,105]
    I = basis_a.orb_l; J = basis_b.orb_l
    a = basis_a.α_array; b = basis_b.α_array
    da = basis_a.d_array; db = basis_b.d_array
    Sij = [[S_ij(I, J, basis_a.cord[1],basis_b.cord[1], ai, bi) for bi = b] for ai = a]
    Skl = [[S_ij(I, J, basis_a.cord[2],basis_b.cord[2], ai, bi) for bi = b] for ai = a]
    Smn = [[S_ij(I, J, basis_a.cord[3],basis_b.cord[3], ai, bi) for bi = b] for ai = a]

    Sab = zeros(Float64,(length(da), length(db),I+1, J+1, I+1, J+1))
    for ind_a = 1:length(da), ind_b = 1:length(db)
        for i = 0:I, j = 0:J, k = 0:I-i,l = 0:J-j
            m = I-i-k
            n = J-j-l
            for p = 1:length(basis_a.α_array), q = 1:length(basis_b.α_array)
                ans = da[ind_a][p]*db[ind_b][q]*Sij[p][q][i+1,j+1]*Skl[p][q][k+1,l+1]*Smn[p][q][m+1,n+1]
                Na = (2*a[p]/pi)^(3/4)*sqrt(((4*a[p])^I)/w_fact[I+1])
                Nb = (2*b[q]/pi)^(3/4)*sqrt(((4*b[q])^J)/w_fact[J+1])
                ans *= Na *Nb
                Sab[ind_a, ind_b,i+1,j+1,k+1,l+1] += ans
            end
        end
    end
    return Sab
end

function S_lm(basis_a::orb_detail, basis_b::orb_detail, C_a::Array{Float64,4}, C_b::Array{Float64,4})
    Sab = cont_Sij(basis_a, basis_b)
    
    la = basis_a.orb_l; lb = basis_b.orb_l
    a = basis_a.α_array; b = basis_b.α_array
    da = basis_a.d_array; db = basis_b.d_array

    max_l = max(basis_a.orb_l, basis_b.orb_l)
    fact = [factorial(i) for i = 0:2*max_l]
    
    S_mamb = zeros(Float64, (length(da), length(db), 2*basis_a.orb_l+1, 2*basis_b.orb_l+1))
    
    for i = 0:2la, j = 0:2lb
        ma = i-la; mb = j-lb
        ma_ = abs(ma);mb_ = abs(mb)
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
                    S_mamb[ind_a, ind_b,i+1, j+1] += f*C_a[ma_+1,ta+1,ua+1,Int(floor(2*va_))+1]*C_b[mb_+1,tb+1,ub+1,Int(floor(2*vb_))+1]*Sab[ind_a, ind_b, pow_xa+1,pow_xb+1,pow_ya+1,pow_yb+1]
                end

                if flag
                    break
                end
            end
        end
        for ind_a = 1:length(da), ind_b = 1:length(db)
            S_mamb[ind_a, ind_b,i+1,j+1] *= Nma*Nmb
        end
    end
    return S_mamb    

end

function get_ovlp(mol::Molecule)
    basis = mol.basis
    S = zeros(Float64, (mol.basis_num, mol.basis_num))
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
                Slm = S_lm(basis[i], basis[j], mol.rotate_coef[la+1], mol.rotate_coef[lb+1])
                for ind_a = 0:length(basis[i].d_array)-1, ind_b = 0:length(basis[j].d_array)-1
                    for k = 0:2*la, l = 0:2*lb
                        S[ind_i + ind_a*(2la+1) + change[la+1][k+1], ind_j + ind_b*(2lb+1) + change[lb+1][l+1]] = Slm[ind_a+1,ind_b+1,k+1,l+1]
                        S[ind_j + ind_b*(2lb+1) + change[lb+1][l+1], ind_i + ind_a*(2la+1) + change[la+1][k+1]] = S[ind_i + ind_a*(2la+1) + change[la+1][k+1], ind_j + ind_b*(2lb+1) + change[lb+1][l+1]]
                    end
                end
            end
            ind_j += (2lb+1)*length(basis[j].d_array)
        end
        ind_i += (2la+1)*length(basis[i].d_array)
    end
    return S
end

end