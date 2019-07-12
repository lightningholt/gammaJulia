
#This finds the linear Approximation of a 2D SYSTEM
function gammaLinApprox(N, rcpt_types, fs, r_t, c, J0, i2e)
    Nthetas = round(Int32, N/2)
    W = J0

    #systematic paramters
    kk = 0.04
    nn = 2
    cons = length(c)

    fs = fs/1000

    r_star = r_t[end,:,:]
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

    SpectE = zeros(length(fs), cons)
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
        tauSvec = kron(tauS, ones(N,1)); #vector time-scales

        #Tracker doesn't like this version
        # Wtot = [(1-nmdaRatio)*W[:,1] zeros(N, N*rcpt_types-1);
        #     zeros(N,N) nmdaRatio*W[:,1] zeros(N, N+1);
        #     zeros(N, N*rcpt_types-1) W[:,2]]
        J = [-1+(1-nmdaRatio)*W[1,1]*rs[1,1] 0 0 0 0 0 -1+(1-nmdaRatio)*W[1,1]*rs[1,2] 0 0 0 0 0 -1+(1-nmdaRatio)*W[1,1]*rs[1,3] 0 0 0 0 0 -1+(1-nmdaRatio)*W[1,1]*rs[1,4] 0 0 0 0 0;
               (1-nmdaRatio)*W[2,1]*rs[1,1] -1 0 0 0 0 (1-nmdaRatio)*W[2,1]*rs[1,2] -1 0 0 0 0 (1-nmdaRatio)*W[2,1]*rs[1,3] -1 0 0 0 0 (1-nmdaRatio)*W[2,1]*rs[1,4] -1 0 0 0 0;
               0 0 nmdaRatio*W[1,1]*rs[1,1]-1 0 0 0 0 0 nmdaRatio*W[1,1]*rs[1,2]-1 0 0 0 0 0 nmdaRatio*W[1,1]*rs[1,3]-1 0 0 0 0 0 nmdaRatio*W[1,1]*rs[1,4]-1 0 0 0;
               0 0 nmdaRatio*W[2,1]*rs[1,1] -1 0 0 0 0 nmdaRatio*W[2,1]*rs[1,2] -1 0 0 0 0 nmdaRatio*W[2,1]*rs[1,3] -1 0 0 0 0 nmdaRatio*W[2,1]*rs[1,4] -1 0 0;
               0 0 0 0 -1 W[1,2]*rs[2,1] 0 0 0 0 -1 W[1,2]*rs[2,2] 0 0 0 0 -1 W[1,2]*rs[2,3] 0 0 0 0 -1 W[1,2]*rs[2,4];
               0 0 0 0 0 W[2,2]*rs[2,1]-1 0 0 0 0 0 W[2,2]*rs[2,2]-1 0 0 0 0 0 W[2,2]*rs[2,3]-1 0 0 0 0 0 W[2,2]*rs[2,4]-1]


    else
        Wrcpt = W;
        Wtot = W;
        tauSvec = tau;
    end

    # Phi = [r_star[1,1] 0 r_star[1, 2] 0 r_star[1,3] 0 r_star[1,4]; 0  r_star[2,1] 0  r_star[2,2] 0  r_star[2,3] 0  r_star[2,4]]
    J = reshape(J, N*rcpt_types, N*rcpt_types, cons)

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

    return SpectE
end
