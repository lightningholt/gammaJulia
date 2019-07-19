
#This finds the linear Approximation of a 2D SYSTEM
function gammaLinApprox(N, rcpt_types, fs, vv_t, c, J0, i2e)
    Nthetas = round(Int32, N/2)
    W = J0

    #systematic paramters
    kk = 0.04
    nn = 2
    cons = length(c)

    fs = fs/1000
    # r_star = r_t[end,:,:]

    if length(size(vv_t)) > 2
        vv = vv_t[end, :,:]
    else
        vv = vv_t
    end

    r_star =  kk*max.([sum(vv[1:2:end,:], dims=1); sum(vv[2:2:end,:], dims=1)], zeros(N, cons)).^nn

    SpectE = zeros(length(fs), cons)
    #SpectI = zeros(length(fs), cons)

    for cc in 1:cons
        Phi = Diagonal(nn*kk.^(1/nn)*r_star[:,cc].^(1-1/nn))

        NoiseCov = ones(N)

        SpectorE = zeros(length(fs))
        #SpectorI = zeros(length(fs))

        eE = [1; 0]
        eI = [0; 1]

        #time constants
        tauE = 15
        tau_ratio = 1
        tauI = tauE/tau_ratio

        #multiple synaptic types
        t_scale = 1
        tauNMDA = 100*t_scale
        tauAMPA = 3*t_scale
        tauGABA = 5*t_scale
        nmdaRatio = 0.1

        NoiseNMDAratio = 0
        NoiseTau = 1*t_scale

        if rcpt_types > 1

            tauS = [tauAMPA, tauNMDA, tauGABA]
            tauSvec = kron(tauS, ones(N,1)); #vector time-scales

            Wrcpt = zeros(N, N, length(tauS));
            Wrcpt[:,:, 1] = (1-nmdaRatio)*[W[:,1] zeros(N, Nthetas)];
            Wrcpt[:,:,3] = [zeros(N,Nthetas)  W[:, 1+Nthetas]];
            Wrcpt[:,:, 2] = nmdaRatio*[W[:,1] zeros(N, Nthetas)];

            Wtot = zeros(N*rcpt_types, N*rcpt_types);
            Wtot[1:N, 1:N] = Wrcpt[:,:,1]
            Wtot[N+1:N+N, N+1:N+N] = Wrcpt[:,:,2]
            Wtot[2*N+1:3*N, 2*N+1:3*N] = Wrcpt[:,:,3]
        else
            Wrcpt = W;
            Wtot = W;
            tauSvec = tau;
        end


        if rcpt_types == 1
            for ff = 1:length(fs)
                Gf = conj((-1j*2*pi*fs[ff]*tau+1).^(-1)- W*Phi);
                vecE = Gf*eE;
                vecI = Gf*eI;

                SpectorE[ff,:] = vecE'*(NoiseCov.*conj(vecE));
                SpectorI[ff, :] = vecI'*(NoiseCov.*conj(vecI));

            end
        else
            J = [Wrcpt[:,:,1]*Phi; Wrcpt[:,:,2]*Phi; Wrcpt[:,:,3]*Phi]
            J = -I + repeat(J, 1, rcpt_types)
            Jacob = Diagonal(kron(inv.(tauS), ones(N)))*J
            JacobLambs, Jacobvec = eigen(Jacob) #eigenvalues and vectors of Jacob

            JacobLambsHz = 1000*JacobLambs/(2*pi)
            println(JacobLambsHz)

            eE = kron(ones(rcpt_types), eE)
            eI = kron(ones(rcpt_types), eI)
            ind = 0

            for ff in fs
                ind += 1

                Gf = -1im*2*pi*ff*Diagonal(kron(tauS, ones(N)))- J
                #Green's function

                vecE = Gf\eE
                vecE = (1-NoiseNMDAratio)*vecE[1:N]+ NoiseNMDAratio*vecE[N+1:2*N]
                SpectorE[ind] = (vecE'*(NoiseCov.*vecE)) * 2 * NoiseTau/abs(-1im *2 * pi * ff * NoiseTau  + 1)^2
            end
        end

        SpectE[:, cc] = SpectorE*2/1000


    end

    return SpectE
end


#This finds the linear Approximation of a 2D SYSTEM
function gammaTrackedLinApprox(N, rcpt_types, fs, vv, c, J0, i2e)
    #NOTE in this case vv is the final voltage from steady_Tracker_v
    Nthetas = round(Int32, N/2)
    W = J0


    #systematic paramters
    kk = 0.04
    nn = 2
    cons = length(c)

    fs = fs/1000


    r_star =  kk*max.([sum(vv[1:2:end,:], dims=1); sum(vv[2:2:end,:], dims=1)], zeros(N, cons)).^nn
    rs = nn * kk.^(1/nn) * r_star.^(1-1/nn)

    if cons == 1
        Phi(rr) = [rr[1] 0; 0 rr[2]]
    end

    # SpectE = zeros(length(fs), cons)
    #SpectI = zeros(length(fs), cons)

    eE = [1; 0]
    eI = [0; 1]

    #time constants
    tauE = 15
    tau_ratio = 1
    tauI = tauE/tau_ratio

    #multiple synaptic types
    t_scale = 1
    tauNMDA = 100*t_scale
    tauAMPA = 3*t_scale
    tauGABA = 5*t_scale
    nmdaRatio = 0.1

    NoiseNMDAratio = 0
    NoiseTau = 1*t_scale

    if rcpt_types > 1

        tauS = [tauAMPA, tauNMDA, tauGABA]
        tauSvec = kron(tauS, ones(1,N)); #vector time-scales

        #Tracker doesn't like this version
        # Wtot = [(1-nmdaRatio)*W[:,1] zeros(N, N*rcpt_types-1);
        #     zeros(N,N) nmdaRatio*W[:,1] zeros(N, N+1);
        #     zeros(N, N*rcpt_types-1) W[:,2]]

        Wtot = [(1-nmdaRatio)*W[1,1] 0 0 0 0 0;
           (1-nmdaRatio)*W[2,1] 0 0 0 0 0;
           0 0 nmdaRatio*W[1,1] 0 0 0;
           0 0 nmdaRatio*W[2,1] 0 0 0;
           0 0 0 0 0 W[1,2];
           0 0 0 0 0 W[2,2]]

        # J = - I + Wtot * kron(ones(rcpt_types, rcpt_type), nn*kk)
        # J = [-1+(1-nmdaRatio)*W[1,1]*rs[1,1] 0 0 0 0 0 -1+(1-nmdaRatio)*W[1,1]*rs[1,2] 0 0 0 0 0 -1+(1-nmdaRatio)*W[1,1]*rs[1,3] 0 0 0 0 0 -1+(1-nmdaRatio)*W[1,1]*rs[1,4] 0 0 0 0 0;
        #        (1-nmdaRatio)*W[2,1]*rs[1,1] -1 0 0 0 0 (1-nmdaRatio)*W[2,1]*rs[1,2] -1 0 0 0 0 (1-nmdaRatio)*W[2,1]*rs[1,3] -1 0 0 0 0 (1-nmdaRatio)*W[2,1]*rs[1,4] -1 0 0 0 0;
        #        0 0 nmdaRatio*W[1,1]*rs[1,1]-1 0 0 0 0 0 nmdaRatio*W[1,1]*rs[1,2]-1 0 0 0 0 0 nmdaRatio*W[1,1]*rs[1,3]-1 0 0 0 0 0 nmdaRatio*W[1,1]*rs[1,4]-1 0 0 0;
        #        0 0 nmdaRatio*W[2,1]*rs[1,1] -1 0 0 0 0 nmdaRatio*W[2,1]*rs[1,2] -1 0 0 0 0 nmdaRatio*W[2,1]*rs[1,3] -1 0 0 0 0 nmdaRatio*W[2,1]*rs[1,4] -1 0 0;
        #        0 0 0 0 -1 W[1,2]*rs[2,1] 0 0 0 0 -1 W[1,2]*rs[2,2] 0 0 0 0 -1 W[1,2]*rs[2,3] 0 0 0 0 -1 W[1,2]*rs[2,4];
        #        0 0 0 0 0 W[2,2]*rs[2,1]-1 0 0 0 0 0 W[2,2]*rs[2,2]-1 0 0 0 0 0 W[2,2]*rs[2,3]-1 0 0 0 0 0 W[2,2]*rs[2,4]-1]


    else
        Wrcpt = W;
        Wtot = W;
        tauSvec = tau;
    end

    #=
    # Phi = [r_star[1,1] 0 r_star[1, 2] 0 r_star[1,3] 0 r_star[1,4]; 0  r_star[2,1] 0  r_star[2,2] 0  r_star[2,3] 0  r_star[2,4]]
    # J = reshape(J, N*rcpt_types, N*rcpt_types, cons)
    J = -I + Wtot*kron(ones(rcpt_types, rcpt_types), nn*kk^(1/nn)*Phi(r_star).^(1-1/nn))

    Gf = -1im * 2 * pi * kron(fs, Diagonal(kron(tauS, ones(N)))) - kron(ones(length(fs)), J)

    vecE = Gf' \ kron(ones(rcpt_types), eE)

    kern = [1-NoiseNMDAratio; 0;  NoiseNMDAratio; 0; 0; 0]
    kern = repeat(kern, length(fs), N*rcpt_types);

    Gf = kern.*Gf

    # Gf = permutedims(reshape(kern, N*rcpt_types, N*rcpt_types, length(fs)), [3,1,2])
    Gf = permutedims(reshape(transpose(Gf .* kern), N*rcpt_types, N*rcpt_types, length(fs)), [2,1, 3])
    =#

    PS(NoiseNMDAratio, NoiseTau, fff, cc, r_star, Wtot,nn, kk, rcpt_types, tauS) = traceGreen(NoiseNMDAratio, NoiseTau, fff, makeGf(fff, cc, r_star, Wtot, nn, kk, rcpt_types, tauS))

    SpectE = zeros(length(fs), cons)

    SpectE = PS.(NoiseNMDAratio, NoiseTau, fs, 4, r_star, Wtot, nn, kk, rcpt_types, tauS)

    # for cc = 1:cons
    #     ind = 0
    #     for ff in fs
    #         ind += 1
    #         SpectE[ind, cc] = PS(NoiseNMDAratio, NoiseTau, ff, cc, r_star, Wtot, nn, kk, rcpt_types, tauS)
    #     end
    # end

    SpectE = 2*SpectE/1000

    # SpectE = PS.(NoiseNMDAratio, NoiseTau, fs, 1:cons, r_star, Wtot,nn, kk, rcpt_types, tauS)

    #=

    SpectE = sum(abs2, Gf, dims=[ 1, 2])

    # SpecE = greenInvPS(fs, tauS, Wtot, r_star)
    =#

    #=
    for cc in 1:cons

        # make a N*rcpt_types x N*rcpt_types with blocks of Phi everywhere
        # Phi = kron(ones(rcpt_types, rcpt_types), Diagonal(nn*kk.^(1/nn)*r_star[:,cc].^(1-1/nn)))


        NoiseCov = ones(N)

        SpectorE = zeros(length(fs))
        #SpectorI = zeros(length(fs))

        if rcpt_types == 1
            for ff = 1:length(fs)
                Gf = conj((-1j*2*pi*fs[ff]*tau+1).^(-1)- W*Phi);
                vecE = Gf*eE;
                vecI = Gf*eI;

                SpectorE[ff,:] = vecE'*(NoiseCov.*conj(vecE));
                SpectorI[ff, :] = vecI'*(NoiseCov.*conj(vecI));

            end
        else
            #J = [Wrcpt[:,:,1]*Phi; Wrcpt[:,:,2]*Phi; Wrcpt[:,:,3]*Phi]

            # J = -I + Wtot*Phi

            #J = -I + repeat(J, 1, rcpt_types)

            Jcon = J[:,:, cc]

            # Jacob = Diagonal(kron(inv.(tauS), ones(N)))*Jcon.data
            # JacobLambs, Jacobvec = eigen(Jacob) #eigenvalues and vectors of Jacob
            #
            # JacobLambsHz = 1000*JacobLambs/(2*pi)
            # println(JacobLambsHz)

            eE = kron(ones(rcpt_types), eE)
            eI = kron(ones(rcpt_types), eI)
            ind = 0

            Gf = -1im * 2 * pi * kron(fs, Diagonal(kron(tauS, ones(N)))) - kron(ones(length(fs)), Jcon)
            reGF = permutedims(reshape(transpose(Gf), N*rcpt_type, N*rcpt_type, length(fs)), [2,1,3])

            for ff in fs
                ind += 1

                # Gf = -1im*2*pi*ff*Diagonal(kron(tauS, ones(N)))- J
                Gf = -1im*2*pi*ff*Diagonal(kron(tauS, ones(N)))- Jcon
                #Green's function

                vecE = Gf\eE
                vecE = (1-NoiseNMDAratio)*vecE[1:N]+ NoiseNMDAratio*vecE[N+1:2*N]
                SpectorE[ind] = (vecE'*(NoiseCov.*vecE)) * 2 * NoiseTau/abs(-1im *2 * pi * ff * NoiseTau  + 1)^2
            end
        end

        SpectE[:, cc] = SpectorE*2/1000


    end
    =#

    return SpectE

end


function greenInvPS(fs, tauS, W, r_star)

    rcpt_types = length(tauS)
    nn = 2
    kk = 0.04
    NoiseNMDAratio = 0
    NoiseTau = 1* 100/maximum(tauS) # 100/max(tauS) should give me the time_scale

    N = 2

    cons = size(r_star, 2)

    # if rcpt_types >= 1
    #     N = round(Int64, size(W,1)/rcpt_types)
    # else
    #     N = 2
    #     rcpt_types = 1 #so I don't break things later
    # end

    eE = [1; 0]
    eE = kron(ones(rcpt_types), eE)

    NoiseCov = ones(N)

    println(size(W))
    println(rcpt_types)
    println(r_star)

    SpectorE = zeros(length(fs), cons)

    for cc = 1:cons

        J = -I + W * kron(ones(rcpt_types, rcpt_types), nn*kk^(1/nn)*Diagonal(r_star[:,cc]).^(1-1/nn))
        ind = 0

        for ff in fs

            ind += 1
            Gf = -1im * 2 * pi * ff * Diagonal(kron(tauS, ones(N))) -  J
            # Gf = -1im * 2 * pi * kron(fs, Diagonal(kron(tauS, ones(N)))) - kron(ones(length(fs)), J)


            vecE = Gf \ eE
            vecE = (1-NoiseNMDAratio)*vecE[1:N]+ NoiseNMDAratio*vecE[N+1:2*N]
            SpectorE[ind, cc] = (vecE'*(NoiseCov.*vecE)) * 2 * NoiseTau/abs(-1im *2 * pi * ff * NoiseTau  + 1)^2
        end
    end

    return SpectorE
end


function traceGreen(NoiseNMDAratio, NoiseTau, ff, Gf)
    N = 2
    rcpt_types = 3
    eE = kron(ones(rcpt_types), [1; 0])


    #=
    Gf = inv(Gf)

    eE = [1 - NoiseNMDAratio; 0; NoiseNMDAratio; 0; 0; 0]


    SpectE = eE' * (Gf' * Gf) * eE * 2 * NoiseTau/abs(-1im *2 * pi * fs * NoiseTau  + 1)^2
    =#
    NoiseCov = ones(N)
    vecE = Gf\ eE
    vecE = (1-NoiseNMDAratio)*vecE[1:N]+ NoiseNMDAratio*vecE[N+1:2*N]

    PS_fc = (vecE'*(NoiseCov.*vecE)) * 2 * NoiseTau/abs(-1im *2 * pi * ff * NoiseTau  + 1)^2

    return PS_fc
end


    # sGf, indGf = findmax(size(Gf))
    #
    # if sGf > N*rcpt_types
    #     if indGf == 1
    #         Gf = permutedims(reshape(transpose(Gf), N*rcpt_types, N*rcpt_types, sGf)
    #     elseif indGf == 2
    #         Gf = permutedims(reshape(Gf, N*rcpt_types, N*rcpt_types, sGf))
    #     end
    # end


function makeGf(ff, cc, r_star, Wtot, nn, kk, rcpt_types, tauS)
    makePhi(r) = [r[1] 0; 0 r[2]]
    phi = nn*kk^(1/nn)*makePhi(r_star[:,cc]).^(1-1/nn)

    J = -I + Wtot*kron(ones(rcpt_types, rcpt_types), phi)
    Gf = -1im * 2 * pi * kron(ff, Diagonal(kron(tauS, ones(N)))) - kron(ones(length(ff)),  J)
    # tauMat = [tauS[1] 0 0 0 0 0; 0 tauS[1] 0 0 0 0; 0 0 tauS[2] 0 0 0; 0 0 0 tauS[2] 0 0; 0 0 0 0 tauS[3] 0; 0 0 0 0 0 tauS[3]]
    # Gf = -1im * 2 * pi * ff * tauMat  - J
    return Gf
end


function sparseLinApprox(N, rcpt_types, fs, vv_t, c, J0, i2e)
    Nthetas = round(Int32, N/2)
    # J0 = param(J0)
    W = J0


    #systematic paramters
    kk = 0.04
    nn = 2
    cons = length(c)

    fs = fs/1000


    r_star =  kk*max.([sum(v1[1:2:end,:], dims=1); sum(v1[2:2:end,:], dims=1)], zeros(N, cons)).^nn
    rs = nn * kk.^(1/nn) * r_star.^(1-1/nn)
    
        #time constants
    tauE = 15
    tau_ratio = 1
    tauI = tauE/tau_ratio

    #multiple synaptic types
    t_scale = 1
    tauNMDA = 100*t_scale
    tauAMPA = 3*t_scale
    tauGABA = 5*t_scale
    nmdaRatio = 0.1

    NoiseNMDAratio = 0
    NoiseTau = 1*t_scale
    NoiseCov = [1; 1]

    if rcpt_types > 1

        tauS = [tauAMPA, tauNMDA, tauGABA]
        tauSvec = kron(tauS, ones(1,N)); #vector time-scales
    #      Wtot = [(1-nmdaRatio)*W[1,1] 0 0 0 0 0;
    #            (1-nmdaRatio)*W[2,1] 0 0 0 0 0;
    #            0 0 nmdaRatio*W[1,1] 0 0 0;
    #            0 0 nmdaRatio*W[2,1] 0 0 0;
    #            0 0 0 0 0 W[1,2];
    #            0 0 0 0 0 W[2,2]]
         Wtot = [(1-nmdaRatio)*Jee 0 0 0 0 0;
           (1-nmdaRatio)*Jie 0 0 0 0 0;
           0 0 nmdaRatio*Jee 0 0 0;
           0 0 nmdaRatio*Jie 0 0 0;
           0 0 0 0 0 -Jei;
           0 0 0 0 0 -Jii]

    else
        Wrcpt = W;
        Wtot = W;
        tauSvec = tau;
    end    
    
    eE = [1; 0]
    Phi(rr) = Diagonal(vec(rr))
    # J = -I + Wtot* kron(ones(rcpt_types, rcpt_types), Phi)
    J = [CuArray(-I + Wtot * kron(ones(rcpt_types, rcpt_types), Phi(rs[:,cc]))) for cc in 1:cons]
    eE = kron(ones(rcpt_types), eE)

    Gf = [CuArray(-1im * 2 * pi * ff * Diagonal(kron(tauS, ones(N))) - J[cc]) for cc in 1:cons for ff in fs];
    cuE = [CuArray(eE) for cc in 1:cons for ff in fs];

    fscons = repeat(fs, cons)
    
    SpectE = invToSpect(Gf, cuE, fscons, NoiseCov, NoiseNMDAratio, NoiseTau)
    
    return SpectE
end


function invToSpect(Gf, cuE, fscons, NoiseCov, NoiseNMDAratio, NoiseTau)
    info, invGf = CuArrays.CUBLAS.matinv_batched(Gf)
    
    kernel = let invGf = invGf
        function (i)
            x = invGf[i] * cuE[i]
            y = (1-NoiseNMDAratio)*x[1:N] + NoiseNMDAratio*x[N+1:2*N]
#             y = invGf[i]*(cuE[i].* CuArray([(1-NoiseNMDAratio); 0; NoiseNMDAratio; 0; 0; 0]))
            
            return y' * (CuArray(NoiseCov).*y) * 2 * NoiseTau/abs(-1im *2 * pi * fscons[i] * NoiseTau  + 1)^2
        end
    end
    
    return mapreduce(kernel, vcat, CuArray(axes(Gf, 1)); init=CuArray{Float32, 1}())
end