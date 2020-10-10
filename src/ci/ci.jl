module CI
    using JuliaSCF.Mole
    using JuliaSCF.HF
    using OffsetArrays
    using LinearAlgebra

    export ci, run_ci!

    mutable struct ci
        mol::Molecule
        n::Int
        N::Int
        all::Int
        #C0::Matrix{Float64}
        C0::Vector{Float64}
        c_dav::Vector{Float64}
        #sigma0::Matrix{Float64}
        sigma0::Vector{Float64}
        H0::Vector{Float64}
        hcore::Matrix{Float64}
        v2e_mo::Array{Float64, 4}
        graph::AbstractArray{Int, 2}
        list::Vector{Vector{Int}}
        level_ind::Vector{Vector{Int}}
        E_0::Float64
        Ene::Float64
        nuc_ene::Float64

        function ci(kernel::hf)
            h = kernel.coeff'*kernel.hcore*kernel.coeff
            v2e_mo_ = zero(kernel.v2e_ao)
            coeff = kernel.coeff
            n = kernel.mol.basis_num
            N = div(kernel.mol.elec_num, 2)
            single = N*(n-N)
            double = div(N*(N-1),2)*div((n-N)*(n-N-1),2)
            all = 1+N*(n-N)+div(N*(N-1),2)*div((n-N)*(n-N-1),2)
            v2e_trans = reshape(kernel.v2e_ao, (n*n,n*n))
            coeff_kron = kron(coeff, coeff)
            v2e_trans = coeff_kron'*v2e_trans*coeff_kron
            v2e_mo_ = reshape(v2e_trans, (n,n,n,n))

            graph = make_graph(n,N)
            level_ind, list = make_list(n,N,graph)
            C0 = zeros(Float64, all+single*(1+single)+double)
            C0[1] = 1.0
            sigma = zeros(Float64, all+single*(1+single)+double)
            for k = 1:n, l = 1:n
                for j = 1:n
                    h[k,l] -= (1/2)*v2e_mo_[k,j,j,l]
                end
            end
            E0 = kernel.elec_ene
            nuc_ene = kernel.nuc_ene
            new(kernel.mol, n, N, all, C0, [0.0], sigma, [0.0], h, v2e_mo_, graph, list, level_ind, E0, E0, kernel.nuc_ene)
        end
    end

    mutable struct ci_box
        c::Vector{Vector{Float64}}
        sigma::Vector{Vector{Float64}}
        function ci_box()
            new(Vector{Vector{Float64}}(), Vector{Vector{Float64}}())
        end
    end

    function make_graph(n::Int, N::Int)
        level = 2 #excitation level. now CISD
        graph = OffsetArray(zeros(Int, (n+1,N+1)), 0:n, 0:N)
        graph[0,0] = 1
        for i = 1:n
            if i <= N
                ind = max(0,i-level)
            else
                ind = max(N-level,i-n+N,0)
            end
            for j = ind:min(i,N)
                if j == 0
                    graph[i,j] += graph[i-1,j]
                else
                    graph[i,j] += graph[i-1,j-1] + graph[i-1,j]
                end
            end
        end
        return graph
    end

    function use_graph(list::Vector{Int}, graph::AbstractArray{Int,2})
        result = 1
        elec_num = 0
        for (index, occ) in enumerate(list)
            if occ == 1
                elec_num += 1
                result += graph[index-1,elec_num]
            end
        end
        return result
    end

    function make_list(n::Int, N::Int, graph::AbstractArray{Int,2})
        level = 2 #excitation level. now CISD
        s_num = N*(n-N)
        d_num = div(N*(N-1),2)*div((n-N)*(n-N-1),2)
        list = Vector{Vector{Int}}(undef,1+s_num+d_num)
        # ground
        vec = zeros(Int, n)
        for ind in 1:N vec[ind] = 1 end
        ind = use_graph(vec, graph)
        list[ind] = vec
        ground = [ind]
        # s-level
        s_level = Vector{Int}()
        for i = 1:N, j = N+1:n
            vec = zeros(Int, n)
            for ind in 1:N vec[ind] = 1 end
            vec[i] = 0; vec[j] = 1
            ind = use_graph(vec, graph)
            list[ind] = vec
            push!(s_level, ind)
        end
        # d-level
        d_level = Vector{Int}()
        for i = 1:N, k = N+1:n
            for j = i+1:N, l = k+1:n
                vec = zeros(Int, n)
                for ind in 1:N vec[ind] = 1 end
                vec[i] = 0; vec[k] = 1
                vec[j] = 0; vec[l] = 1
                ind = use_graph(vec, graph)
                list[ind] = vec
                push!(d_level, ind)
            end
        end
        return [ground, s_level, d_level], list
    end
    
    function find_excitaion(n::Int, N::Int, vec::Vector{Int})
        occ = zeros(Int, N); occ_ind = 1
        unocc = zeros(Int, n-N); unocc_ind = 1
        for (ind, o) in enumerate(vec)
            if o == 1
                occ[occ_ind] = ind
                occ_ind+=1
            else
                unocc[unocc_ind] = ind
                unocc_ind+=1
            end
        end
        accum = accumulate(+, vec)
        
        ex = Vector{Tuple{Int,Int}}()
        sgn = Vector{Int}()
        sd = Vector{Bool}()
        for o in occ
            push!(ex, (o,o))
            push!(sgn, 1)
            if accum[N] >= N-2
                push!(sd, true)
            else
                push!(sd, false)
            end
            for u in unocc
                push!(ex, (u,o))
                parity1 = accum[o]-1
                parity2 = accum[u]
                if o < u
                    parity2 -= 1
                end
                push!(sgn,1-2*mod(parity1 - parity2,2))
                if (accum[N]-(o<=N ? 1 : 0)+(u<=N ? 1 : 0)) >= N-2
                    push!(sd, true)
                else
                    push!(sd, false)
                end
            end
        end
        return ex, sgn, sd
    end

    function find_excitaion2(n::Int, N::Int, vec::Vector{Int})
        occ = zeros(Int, N); occ_ind = 1
        unocc = zeros(Int, n-N); unocc_ind = 1
        for (ind, o) in enumerate(vec)
            if o == 1
                occ[occ_ind] = ind
                occ_ind+=1
            else
                unocc[unocc_ind] = ind
                unocc_ind+=1
            end
        end
        accum = accumulate(+, vec)
        
        ex = Vector{Tuple{Int,Int}}()
        sgn = Vector{Int}()
        sd = Vector{Bool}()
        for o in occ
            if accum[N] >= N-2
                push!(ex, (o,o))
                push!(sgn, 1)
                push!(sd, true)
            end
            for u in unocc
                if (accum[N]-(o<=N ? 1 : 0)+(u<=N ? 1 : 0)) >= N-2
                    push!(ex, (u,o))
                    parity1 = accum[o]-1
                    parity2 = accum[u]
                    if o < u
                        parity2 -= 1
                    end
                    push!(sgn,1-2*mod(parity1 - parity2,2))
                
                    push!(sd, true)
                end
            end
        end
        return ex, sgn, sd
    end

    function find_occ(n::Int, N::Int, vec::Vector{Int})
        occ = zeros(Int, N); occ_ind = 1
        unocc = zeros(Int, n-N); unocc_ind = 1
        for (ind, o) in enumerate(vec)
            if o == 1
                occ[occ_ind] = ind
                occ_ind+=1
            else
                unocc[unocc_ind] = ind
                unocc_ind+=1
            end
        end
        return occ, unocc
    end

    function Sigma(kernel::ci, C::Matrix{Float64})
        sigma = zeros(Float64, (kernel.all,kernel.all))
        # sigma1 
        for level = 0:2
            for I_beta in kernel.level_ind[level+1]
                F = zeros(Float64, kernel.all)
                I_beta_ket =  kernel.list[I_beta]
                exitations1, sgns1, checks1 = find_excitaion(kernel.n, kernel.N, I_beta_ket)
                for (ex_kl, sgn_kl, check1) in zip(exitations1, sgns1, checks1)
                    K_beta_ket = copy(I_beta_ket)
                    K_beta_ket[ex_kl[2]] = 0; K_beta_ket[ex_kl[1]] = 1
                    if check1
                        F[use_graph(K_beta_ket, kernel.graph)] += sgn_kl*kernel.hcore[ex_kl[1],ex_kl[2]]
                    end
                    exitations2, sgns2, checks2 = find_excitaion2(kernel.n, kernel.N, K_beta_ket)
                    for (ex_ij, sgn_ij, check2) in zip(exitations2, sgns2, checks2)
                        if !check2
                            continue
                        end
                        J_beta_ket = copy(K_beta_ket)
                        J_beta_ket[ex_ij[2]] = 0; J_beta_ket[ex_ij[1]] = 1
                        F[use_graph(J_beta_ket, kernel.graph)] += (1/2)*sgn_kl*sgn_ij*kernel.v2e_mo[ex_ij[1],ex_ij[2],ex_kl[1],ex_kl[2]]
                    end
                end
        
                for I_alpha in kernel.level_ind[1]
                    sigma[I_alpha, I_beta] = F'*C[I_alpha, :]
                end
                if level <= 1
                    for I_alpha in kernel.level_ind[2]
                        sigma[I_alpha, I_beta] = F'*C[I_alpha, :]
                    end
                end
                if level == 0
                    for I_alpha in kernel.level_ind[3]
                        sigma[I_alpha, I_beta] = F'*C[I_alpha, :]
                    end
                end
            end
        end
        #sigma2
        for level = 0:2
            for I_alpha in kernel.level_ind[level+1]
                F = zeros(Float64, kernel.all)
                I_alpha_ket =  kernel.list[I_alpha]
                exitations1, sgns1, checks1 = find_excitaion(kernel.n, kernel.N, I_alpha_ket)
                for (ex_kl, sgn_kl, check1) in zip(exitations1, sgns1, checks1)
                    K_alpha_ket = copy(I_alpha_ket)
                    K_alpha_ket[ex_kl[2]] = 0; K_alpha_ket[ex_kl[1]] = 1
                    if check1
                        F[use_graph(K_alpha_ket, kernel.graph)] += sgn_kl*kernel.hcore[ex_kl[1],ex_kl[2]]
                    end
                    exitations2, sgns2, checks2 = find_excitaion2(kernel.n, kernel.N, K_alpha_ket)
                    for (ex_ij, sgn_ij, check2) in zip(exitations2, sgns2, checks2)
                        if !check2
                            continue
                        end
                        J_alpha_ket = copy(K_alpha_ket)
                        J_alpha_ket[ex_ij[2]] = 0; J_alpha_ket[ex_ij[1]] = 1
                        F[use_graph(J_alpha_ket, kernel.graph)] += (1/2)*sgn_kl*sgn_ij*kernel.v2e_mo[ex_ij[1],ex_ij[2],ex_kl[1],ex_kl[2]]
                    end
                end
        
                for I_beta in kernel.level_ind[1]
                    sigma[I_alpha, I_beta] += F'*C[:, I_beta]
                end
                if level <= 1
                    for I_beta in kernel.level_ind[2]
                        sigma[I_alpha, I_beta] += F'*C[:, I_beta]
                    end
                end
                if level == 0
                    for I_beta in kernel.level_ind[3]
                        sigma[I_alpha, I_beta] += F'*C[:, I_beta]
                    end
                end
            end
        end
        
        #sigma 3
        for alpha_level = 0:2
            for I_alpha in kernel.level_ind[alpha_level+1]
                I_alpha_ket =  kernel.list[I_alpha]
                exitations1, sgns1, checks1 = find_excitaion2(kernel.n, kernel.N, I_alpha_ket)
                for (ex_kl, sgn_kl, check1) in zip(exitations1, sgns1, checks1)
                    if !check1
                        continue
                    end
                    J_alpha_ket = copy(I_alpha_ket)
                    J_alpha_ket[ex_kl[2]] = 0; J_alpha_ket[ex_kl[1]] = 1
                    J_alpha = use_graph(J_alpha_ket, kernel.graph)

                    I_beta_list = copy(kernel.level_ind[begin])
                    if alpha_level <= 1
                        union!(I_beta_list, kernel.level_ind[begin+1])
                    end
                    if alpha_level == 0
                        union!(I_beta_list, kernel.level_ind[begin+2])
                    end
                    for I_beta in I_beta_list
                        I_beta_ket =  kernel.list[I_beta]
                        exitations2, sgns2, checks2 = find_excitaion2(kernel.n, kernel.N, I_beta_ket)
                        for (ex_ij, sgn_ij, check2) in zip(exitations2, sgns2, checks2)
                            if !check2
                                continue
                            end
                            J_beta_ket = copy(I_beta_ket)
                            J_beta_ket[ex_ij[2]] = 0; J_beta_ket[ex_ij[1]] = 1
                            J_beta = use_graph(J_beta_ket, kernel.graph)
                            sigma[I_alpha, I_beta] += sgn_ij*sgn_kl*kernel.v2e_mo[ex_ij[1],ex_ij[2],ex_kl[1],ex_kl[2]]*C[J_alpha, J_beta]
                        end
                        
                    end
                end
            end
        end
        return sigma
    end

    function transform(kernel::ci, sigma::Matrix{Float64})
        trans_sigma = Vector{Float64}()
        for alpha_level= 1:3
            for I_alpha in kernel.level_ind[alpha_level]
                if alpha_level <= 3
                    for I_beta in kernel.level_ind[1]
                        push!(trans_sigma,sigma[I_alpha, I_beta])
                    end
                end
                if alpha_level <= 2
                    for I_beta in kernel.level_ind[2]
                        push!(trans_sigma,sigma[I_alpha, I_beta])
                    end
                end
                if alpha_level <= 1
                    for I_beta in kernel.level_ind[3]
                        push!(trans_sigma,sigma[I_alpha, I_beta])
                    end
                end
            end
        end
        return trans_sigma
    end

    function reverse_transform(kernel::ci, trans_sigma::Vector{Float64})
        sigma = zeros(Float64, (kernel.all, kernel.all))
        Ind = 1
        for alpha_level= 1:3
            for I_alpha in kernel.level_ind[alpha_level]
                if alpha_level <= 3
                    for I_beta in kernel.level_ind[1]
                        sigma[I_alpha, I_beta] = trans_sigma[Ind]
                        Ind += 1
                    end
                end
                if alpha_level <= 2
                    for I_beta in kernel.level_ind[2]
                        sigma[I_alpha, I_beta] = trans_sigma[Ind]
                        Ind += 1
                    end
                end
                if alpha_level <= 1
                    for I_beta in kernel.level_ind[3]
                        sigma[I_alpha, I_beta] = trans_sigma[Ind]
                        Ind += 1
                    end
                end
            end
        end
        return sigma
    end

    function make_H0(kernel::ci)
        H0 = Vector{Float64}()
        for alpha_level= 1:3
            for I_alpha in kernel.level_ind[alpha_level]
                I_alpha_ket = kernel.list[I_alpha]
                alpha_occ, alpha_unocc = find_occ(kernel.n, kernel.N, I_alpha_ket)
                h2 = sum([kernel.hcore[i,i] for i in alpha_occ])
                tmp = 0.0
                for i in alpha_occ
                    for j in alpha_unocc
                        tmp += kernel.v2e_mo[i,j,j,i]
                    end
                    for j in alpha_occ
                        tmp += kernel.v2e_mo[i,i,j,j]
                    end
                end
                h2 += tmp/2
                if alpha_level <= 3
                    for I_beta in kernel.level_ind[1]
                        I_beta_ket = kernel.list[I_beta]
                        beta_occ, beta_unocc = find_occ(kernel.n, kernel.N, I_beta_ket)

                        h1 = sum([kernel.hcore[i,i] for i in beta_occ])
                        tmp = 0.0
                        for i in beta_occ
                            for j in beta_unocc
                                tmp += kernel.v2e_mo[i,j,j,i]
                            end
                            for j in beta_occ
                                tmp += kernel.v2e_mo[i,i,j,j]
                            end
                        end
                        h1 += tmp/2
                        h3 = 0.0
                        for i in alpha_occ, j in beta_occ
                            h3 += kernel.v2e_mo[i,i,j,j]
                        end
                        push!(H0, h1+h2+h3)
                    end
                end
                if alpha_level <= 2
                    for I_beta in kernel.level_ind[2]
                        I_beta_ket = kernel.list[I_beta]
                        beta_occ, beta_unocc = find_occ(kernel.n, kernel.N, I_beta_ket)
        
                        h1 = sum([kernel.hcore[i,i] for i in beta_occ])
                        tmp = 0.0
                        for i in beta_occ
                            for j in beta_unocc
                                tmp += kernel.v2e_mo[i,j,j,i]
                            end
                            for j in beta_occ
                                tmp += kernel.v2e_mo[i,i,j,j]
                            end
                        end
                        h1 += tmp/2
                        h3 = 0.0
                        for i in alpha_occ, j in beta_occ
                            h3 += kernel.v2e_mo[i,i,j,j]
                        end
                        push!(H0, h1+h2+h3)
                    end
                end
                if alpha_level <= 1
                    for I_beta in kernel.level_ind[3]
                        I_beta_ket = kernel.list[I_beta]
                        beta_occ, beta_unocc = find_occ(kernel.n, kernel.N, I_beta_ket)
                
                        h1 = sum([kernel.hcore[i,i] for i in beta_occ])
                        tmp = 0.0
                        for i in beta_occ
                            for j in beta_unocc
                                tmp += kernel.v2e_mo[i,j,j,i]
                            end
                            for j in beta_occ
                                tmp += kernel.v2e_mo[i,i,j,j]
                            end
                        end
                        h1 += tmp/2
                        h3 = 0.0
                        for i in alpha_occ, j in beta_occ
                            h3 += kernel.v2e_mo[i,i,j,j]
                        end
                        push!(H0, h1+h2+h3)
                    end
                end
            end
        end
        return H0
    end

    function projection(kernel::ci, vec::Vector{Float64})
        return vec - dot(kernel.trans_C0, vec)*kernel.trans_C0
    end

    function make_c_new(kernel::ci)
        res = kernel.sigma0 - kernel.E_0*kernel.C0
        for (i,h) in enumerate(kernel.H0)
            if h-kernel.E_0 == 0.0
                res[i] /= 1e-17
                continue
            end
            res[i] /= (h-kernel.E_0)
        end
        return -res
    end
    #=
    function make_ene(kernel::ci)
        ene = kernel.E_0
        ene += 2*dot(kernel.c_dav, projection(kernel, kernel.trans_sigma0))
        project_c_dav = projection(kernel, kernel.c_dav)
        ene += dot(kernel.c_dav, projection(kernel, transform(kernel, Sigma(kernel, reverse_transform(kernel, project_c_dav)))))
        ene /= 1+dot(kernel.c_dav, project_c_dav)
        @show ene
        return ene
    end
    =#
    function make_ene(kernel::ci, box::ci_box)
        c = box.c
        sigma = box.sigma
        n = length(c)
        H = zeros(Float64, (n,n))
        S = zeros(Float64, (n,n))
        for i = 1:n, j = 1:n
            H[i,j] = dot(c[i],sigma[j])
            S[i,j] = dot(c[i], c[j])
        end
        vals, vecs = eigen(H,S)
        coeff = vecs[:,1]
        trans_C0 = zero(c[1])
        trans_sigma = zero(sigma[1])
        for i = 1:n
            trans_C0 += coeff[i]*c[i]
            trans_sigma += coeff[i]*sigma[i]
        end
        return vals[1], trans_C0, trans_sigma
    end

    function remake_C0(kernel::ci)
        projection_c_dav = projection(kernel, kernel.c_dav)
        new_C0 = kernel.trans_C0 + projection_c_dav
        new_C0 /= sqrt(1+dot(kernel.c_dav, projection_c_dav))
        return reverse_transform(kernel, new_C0)
    end

    function make_H(kernel::ci)
        n = kernel.n
        N = kernel.N
        s = (n-N)*N
        d = div((n-N)*(n-N-1),2)*div(N*(N-1),2)

        H = zeros(Float64, (1*(1+s+d)+s*(1+s)+d*1,1*(1+s+d)+s*(1+s)+d*1))
        
        ind = 1
        for alpha_level= 1:3
            for I_alpha in kernel.level_ind[alpha_level]
                I_alpha_ket = kernel.list[I_alpha]
                alpha_occ, alpha_unocc = find_occ(kernel.n, kernel.N, I_alpha_ket)
                if alpha_level <= 3
                    for I_beta in kernel.level_ind[1]
                        I_beta_ket = kernel.list[I_beta]
                        make_subH(kernel, H, ind, I_alpha_ket, I_beta_ket)
                        ind += 1
                    end
                end
                if alpha_level <= 2
                    for I_beta in kernel.level_ind[2]
                        I_beta_ket = kernel.list[I_beta]
                        make_subH(kernel, H, ind, I_alpha_ket, I_beta_ket)
                        ind += 1
                    end
                end
                if alpha_level <= 1
                    for I_beta in kernel.level_ind[3]
                        I_beta_ket = kernel.list[I_beta]
                        make_subH(kernel, H, ind, I_alpha_ket, I_beta_ket)
                        ind += 1
                    end
                end
            end
        end
        return H
    end

    function make_subH(kernel::ci, H::Matrix, ind::Int, I_alpha_ket::Vector{Int}, I_beta_ket::Vector{Int})
        ind1 = 1
        for alpha_level= 1:3
            for J_alpha in kernel.level_ind[alpha_level]
                J_alpha_ket = kernel.list[J_alpha]
                if alpha_level <= 3
                    for J_beta in kernel.level_ind[1]
                        J_beta_ket = kernel.list[J_beta]
                        #sig1
                        if I_alpha_ket == J_alpha_ket
                            exitations1, sgns1 = find_excitaion(kernel.n, kernel.N, I_beta_ket)
                            for (ex_kl, sgnkl) in zip(exitations1, sgns1)
                                ket_copy = copy(I_beta_ket)
                                ket_copy[ex_kl[2]] = 0; ket_copy[ex_kl[1]] = 1
                                if ket_copy == J_beta_ket 
                                    H[ind, ind1] += sgnkl*kernel.hcore[ex_kl[1],ex_kl[2]]
                                end
                                exitations2, sgns2 = find_excitaion(kernel.n, kernel.N, J_beta_ket)
                                for (ex_ji, sgnji) in zip(exitations2, sgns2)
                                    bra_copy = copy(J_beta_ket)
                                    bra_copy[ex_ji[2]] = 0; bra_copy[ex_ji[1]] = 1
                                    if bra_copy == ket_copy 
                                        H[ind, ind1] += (1/2)*sgnkl*sgnji*kernel.v2e_mo[ex_ji[2],ex_ji[1],ex_kl[1],ex_kl[2]]
                                    end
                                end
                            end
                        end
                        #sig2
                        if I_beta_ket == J_beta_ket
                            exitations1, sgns1 = find_excitaion(kernel.n, kernel.N, I_alpha_ket)
                            for (ex_kl, sgnkl) in zip(exitations1, sgns1)
                                ket_copy = copy(I_alpha_ket)
                                ket_copy[ex_kl[2]] = 0; ket_copy[ex_kl[1]] = 1
                                if ket_copy == J_alpha_ket
                                    H[ind, ind1] += sgnkl*kernel.hcore[ex_kl[1],ex_kl[2]]
                                end
                                exitations2, sgns2 = find_excitaion(kernel.n, kernel.N, J_alpha_ket)
                                for (ex_ji, sgnji) in zip(exitations2, sgns2)
                                    bra_copy = copy(J_alpha_ket)
                                    bra_copy[ex_ji[2]] = 0; bra_copy[ex_ji[1]] = 1
                                    if bra_copy == ket_copy
                                        H[ind, ind1] += (1/2)*sgnkl*sgnji*kernel.v2e_mo[ex_ji[2],ex_ji[1],ex_kl[1],ex_kl[2]]
                                    end
                                end
                            end
                        end
                        
                        #sig3
                        exitations1, sgns1 = find_excitaion(kernel.n, kernel.N, I_beta_ket)
                        for (ex_ij, sgnij) in zip(exitations1, sgns1)
                            ket_copy1 = copy(I_beta_ket)
                            ket_copy1[ex_ij[2]] = 0; ket_copy1[ex_ij[1]] = 1
                            if ket_copy1 == J_beta_ket
                                exitations2, sgns2 = find_excitaion(kernel.n, kernel.N, I_alpha_ket)
                                for (ex_kl, sgnkl) in zip(exitations2, sgns2)
                                    ket_copy2 = copy(I_alpha_ket)
                                    ket_copy2[ex_kl[2]] = 0; ket_copy2[ex_kl[1]] = 1
                                    if ket_copy2 == J_alpha_ket
                                        H[ind, ind1] += sgnij*sgnkl*kernel.v2e_mo[ex_ij[1],ex_ij[2],ex_kl[1], ex_kl[2]]
                                    end
                                end
                            end
                        end
                        ind1 += 1
                    end
                end
                if alpha_level <= 2
                    for J_beta in kernel.level_ind[2]
                        J_beta_ket = kernel.list[J_beta]
                        #sig1
                        if I_alpha_ket == J_alpha_ket
                            exitations1, sgns1 = find_excitaion(kernel.n, kernel.N, I_beta_ket)
                            for (ex_kl, sgnkl) in zip(exitations1, sgns1)
                                ket_copy = copy(I_beta_ket)
                                ket_copy[ex_kl[2]] = 0; ket_copy[ex_kl[1]] = 1
                                if ket_copy == J_beta_ket
                                    H[ind, ind1] += sgnkl*kernel.hcore[ex_kl[1],ex_kl[2]]
                                end
                                exitations2, sgns2 = find_excitaion(kernel.n, kernel.N, J_beta_ket)
                                for (ex_ji, sgnji) in zip(exitations2, sgns2)
                                    bra_copy = copy(J_beta_ket)
                                    bra_copy[ex_ji[2]] = 0; bra_copy[ex_ji[1]] = 1
                                    if bra_copy == ket_copy
                                        H[ind, ind1] += (1/2)*sgnkl*sgnji*kernel.v2e_mo[ex_ji[2],ex_ji[1],ex_kl[1],ex_kl[2]]
                                    end
                                end
                            end
                        end
                        
                        #sig2
                        if I_beta_ket == J_beta_ket
                            exitations1, sgns1 = find_excitaion(kernel.n, kernel.N, I_alpha_ket)
                            for (ex_kl, sgnkl) in zip(exitations1, sgns1)
                                ket_copy = copy(I_alpha_ket)
                                ket_copy[ex_kl[2]] = 0; ket_copy[ex_kl[1]] = 1
                                if ket_copy == J_alpha_ket
                                    H[ind, ind1] += sgnkl*kernel.hcore[ex_kl[1],ex_kl[2]]
                                end
                                exitations2, sgns2 = find_excitaion(kernel.n, kernel.N, J_alpha_ket)
                                for (ex_ji, sgnji) in zip(exitations2, sgns2)
                                    bra_copy = copy(J_alpha_ket)
                                    bra_copy[ex_ji[2]] = 0; bra_copy[ex_ji[1]] = 1
                                    if bra_copy == ket_copy
                                        H[ind, ind1] += (1/2)*sgnkl*sgnji*kernel.v2e_mo[ex_ji[2],ex_ji[1],ex_kl[1],ex_kl[2]]
                                    end
                                end
                            end
                        end
                        #sig3
                        exitations1, sgns1 = find_excitaion(kernel.n, kernel.N, I_beta_ket)
                        for (ex_ij, sgnij) in zip(exitations1, sgns1)
                            ket_copy1 = copy(I_beta_ket)
                            ket_copy1[ex_ij[2]] = 0; ket_copy1[ex_ij[1]] = 1
                            if ket_copy1 == J_beta_ket
                                exitations2, sgns2 = find_excitaion(kernel.n, kernel.N, I_alpha_ket)
                                for (ex_kl, sgnkl) in zip(exitations2, sgns2)
                                    ket_copy2 = copy(I_alpha_ket)
                                    ket_copy2[ex_kl[2]] = 0; ket_copy2[ex_kl[1]] = 1
                                    if ket_copy2 == J_alpha_ket
                                        H[ind, ind1] += sgnij*sgnkl*kernel.v2e_mo[ex_ij[1],ex_ij[2],ex_kl[1], ex_kl[2]]
                                    end
                                end
                            end
                        end
                        ind1 += 1
                    end
                end
                if alpha_level <= 1
                    for J_beta in kernel.level_ind[3]
                        J_beta_ket = kernel.list[J_beta]
                        #sig1
                        if I_alpha_ket == J_alpha_ket
                            exitations1, sgns1 = find_excitaion(kernel.n, kernel.N, I_beta_ket)
                            for (ex_kl, sgnkl) in zip(exitations1, sgns1)
                                ket_copy = copy(I_beta_ket)
                                ket_copy[ex_kl[2]] = 0; ket_copy[ex_kl[1]] = 1
                                if ket_copy == J_beta_ket
                                    H[ind, ind1] += sgnkl*kernel.hcore[ex_kl[1],ex_kl[2]]
                                end
                                exitations2, sgns2 = find_excitaion(kernel.n, kernel.N, J_beta_ket)
                                for (ex_ji, sgnji) in zip(exitations2, sgns2)
                                    bra_copy = copy(J_beta_ket)
                                    bra_copy[ex_ji[2]] = 0; bra_copy[ex_ji[1]] = 1
                                    if bra_copy == ket_copy
                                        H[ind, ind1] += (1/2)*sgnkl*sgnji*kernel.v2e_mo[ex_ji[2],ex_ji[1],ex_kl[1],ex_kl[2]]
                                    end
                                end
                            end
                        end
                        
                        #sig2
                        if I_beta_ket == J_beta_ket
                            exitations1, sgns1 = find_excitaion(kernel.n, kernel.N, I_alpha_ket)
                            for (ex_kl, sgnkl) in zip(exitations1, sgns1)
                                ket_copy = copy(I_alpha_ket)
                                ket_copy[ex_kl[2]] = 0; ket_copy[ex_kl[1]] = 1
                                if ket_copy == J_alpha_ket
                                    H[ind, ind1] += sgnkl*kernel.hcore[ex_kl[1],ex_kl[2]]
                                end
                                exitations2, sgns2 = find_excitaion(kernel.n, kernel.N, J_alpha_ket)
                                for (ex_ji, sgnji) in zip(exitations2, sgns2)
                                    bra_copy = copy(J_alpha_ket)
                                    bra_copy[ex_ji[2]] = 0; bra_copy[ex_ji[1]] = 1
                                    if bra_copy == ket_copy
                                        H[ind, ind1] += (1/2)*sgnkl*sgnji*kernel.v2e_mo[ex_ji[2],ex_ji[1],ex_kl[1],ex_kl[2]]
                                    end
                                end
                            end
                        end
                        #sig3
                        exitations1, sgns1 = find_excitaion(kernel.n, kernel.N, I_beta_ket)
                        for (ex_ij, sgnij) in zip(exitations1, sgns1)
                            ket_copy1 = copy(I_beta_ket)
                            ket_copy1[ex_ij[2]] = 0; ket_copy1[ex_ij[1]] = 1
                            if ket_copy1 == J_beta_ket
                                exitations2, sgns2 = find_excitaion(kernel.n, kernel.N, I_alpha_ket)
                                for (ex_kl, sgnkl) in zip(exitations2, sgns2)
                                    ket_copy2 = copy(I_alpha_ket)
                                    ket_copy2[ex_kl[2]] = 0; ket_copy2[ex_kl[1]] = 1
                                    if ket_copy2 == J_alpha_ket
                                        H[ind, ind1] += sgnij*sgnkl*kernel.v2e_mo[ex_ij[1],ex_ij[2],ex_kl[1], ex_kl[2]]
                                    end
                                end
                            end
                        end
                        ind1 += 1
                    end
                end
            end
        end
    end

    function run_ci!(kernel::ci)
        init_ene = kernel.E_0
        kernel.H0 = make_H0(kernel)
        box = ci_box()
        
        push!(box.c, copy(kernel.C0))
        kernel.sigma0 = transform(kernel, Sigma(kernel, reverse_transform(kernel, kernel.C0)))
        push!(box.sigma, copy(kernel.sigma0))
        converged = false
        max_iteration = 10
        iter = 0
        while !converged && iter < max_iteration 
            kernel.c_dav = make_c_new(kernel)
            if norm(kernel.c_dav) <= 1e-10
                converged = true
                break
            end
            push!(box.c, copy(kernel.c_dav))
            c_dav_sigma = transform(kernel, Sigma(kernel, reverse_transform(kernel, kernel.c_dav)))
            push!(box.sigma, copy(c_dav_sigma))
            E_old = kernel.E_0
            kernel.E_0, kernel.C0, kernel.sigma0 = make_ene(kernel, box)
            if abs(E_old-kernel.E_0) <= 1e-6
                converged = true
            end
            iter += 1
        end
        if converged
            println("converged!!")
            @show iter
        else
            println("Not converged...")
        end
        kernel.Ene = kernel.nuc_ene + kernel.E_0
        println("Correlation energy")
        corr = kernel.E_0-init_ene
        println(kernel.E_0-init_ene)
        @show kernel.Ene
        return kernel.Ene, corr
    end

    # this is function for debug
    function direct_run(kernel::ci, k::hf)
        H = make_H(kernel)
        vals, vecs = eigen(H)
        @show vals[1]+k.nuc_ene
        #@show vecs[:,1]
        @show vals[1]-k.elec_ene
    end
    #=
    m1 = Molecule([atom("H",[0.0,0.0,-0.7]),atom("Li",[0.0,0.0,0.7])],"ccpvdz")
    k = hf(m1)
    ene1 = run_rhf!(k)
    
    k3 = ci(k)
    run_ci!(k3)

    using PyCall
    @pyimport pyscf
    X = 0.52918
    mol = pyscf.gto.Mole()
    mol.build(atom="H 0 0 $(-0.7*X); Li 0 0 $(0.7*X)", basis="ccpvdz")
    mf = mol.HF().run()
    mycc = mf.CISD().run()
    =#
end