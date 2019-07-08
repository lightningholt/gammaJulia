using LinearAlgebra
using Statistics

#this function uses an Euler Method to find the firing rates
function steady_rates(N, rcpt_types, t, c, J0, i2e)

    Nthetas = round(Int32, N/2)

    W = J0
    println(det(W))

    #systematic parameters
    kk = 0.04
    nn = 2
    cons = length(c)

    #define the inputs
    if rcpt_types > 1
        # g = [ones(N); zeros(N*2)] # Dimensions are [AMPA E I, NMDA E I, GABA E I]

        g = [1; i2e; zeros(N*2)] # Dimensions are [AMPA E I, NMDA E I, GABA E I]
        # g[2:2:end] = i2e .* g[2:2:end]
        # g[2] = i2e .* g[2]
    else
        g = ones(N, rcpt_types) # dimensions are E,I
        g[2,:] = i2e .* g[2,:]
    end

    I_spont = zeros(N*rcpt_types)
    I_mod = zeros(N*rcpt_types)

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

        #define receptor W

        ### Something weird is happening here. Gotta look into this more
        # Wrcpt = zeros(N, N, length(tauS));
        # Wrcpt[:,:, 1] = (1-nmdaRatio) .* [W[:,1] zeros(N)];
        # Wrcpt[:,:,3] = [zeros(N)  W[:, 1+Nthetas]];
        # Wrcpt[:,:, 2] = nmdaRatio*[W[:,1] zeros(N)];

        # Wtot = zeros(N*rcpt_types, N*rcpt_types);
        # Wtot[1:N, 1:N] = Wrcpt[:,:,1]
        # Wtot[N+1:N+N, N+1:N+N] = Wrcpt[:,:,2]
        # Wtot[2*N+1:3*N, 2*N+1:3*N] = Wrcpt[:,:,3]

        Wtot = zeros(N*rcpt_types, N*rcpt_types);
        Wtot[1:N, 1:N] = (1-nmdaRatio) .* [W[:,1] zeros(N)];
        Wtot[N+1:N+N, N+1:N+N] =  nmdaRatio*[W[:,1] zeros(N)];
        Wtot[2*N+1:3*N, 2*N+1:3*N] = [zeros(N)  W[:, 1+Nthetas]];

        #
        # for type = 0:rcpt_types-1
        #     Wtot[(1:N)+type*N, (1:N)+type*N] = Wrcpt[:,:, type+1];
        # end


    else
        tauSvec = tau;
        Wrcpt = W;
        Wtot = W;
    end

    #Time runs
    Nt = length(t)
    totalT = t[end]
    dt = mean(diff(t))
    dt2 = sqrt(dt)

    #Results variables
    v = zeros(N * rcpt_types, cons); #membrane potential
    r_starcons = zeros(N, cons);

    vv_t = zeros(Nt, N * rcpt_types, cons);
    vv1 = zeros(N*rcpt_types,1); #initial membrane
    r_t = zeros(Nt, N, cons); #firing rate
    tt_c = zeros(Nt, cons); #save out the time each simulation run

    rect_powerLaw(vv) = max([sum(vv[1:2:end]), sum(vv[2:2:end])], zeros(2)).^nn
    Conv = true #whether dvdt vanishes as t-> Nt

    #Euler Algorithm for dvdt
    for cc = 1:cons
        I_total = I_spont + I_mod + c[cc].*g

        dvdt(vv) = inv.(tauSvec)*dt .* (-vv + kk*Wtot*kron(ones(rcpt_types), rect_powerLaw(vv))+ I_total)
        vt = zeros(Nt, N*rcpt_types)
        v1 = zeros(N*rcpt_types)

        for tt = 1:Nt
            dv = dvdt(v1)
            v1 = dv + v1
            vv_t[tt, :, cc] = v1

            r_t[tt, :, cc] = kk*rect_powerLaw(v1)

            if tt >= Nt - 1000
                itr = maximum(abs.(dv))

                if itr > 0.01
                    Conv = false
                end
            end

        end

    end

    return  r_t, vv_t, Conv
end
