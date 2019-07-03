
#This finds the linear Approximation of a 2D SYSTEM
function gammaLinApprox(N, rcpt_types, fs, r_t, c, J0, i2e)
    Nthetas = round(Int32, N/2)
    W = J0

    #systematic paramters
    kk = 0.4
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

            eE = kron(ones(rcpt_types), eE)
            eI = kron(ones(rcpt_types), eI)
            ind = 0

            for ff in fs
                ind += 1
                println(ff)

                Gf = -1im*2*pi*ff*Diagonal(kron(tauS, ones(N)))- J
                #Green's function

                vecE = Gf\eE
                vecE = (1-NoiseNMDAratio)*vecE[1:N]+ NoiseNMDAratio*vecE[N+1:2*N]
                SpectorE[ind] = vecE'*(NoiseCov.*vecE) * 2 * NoiseTau/abs(-1im *2 * pi*ff * NoiseTau  + 1)^2
            end
        end

        SpectE[:, cc] = SpectorE


    end

    return SpectE
end
