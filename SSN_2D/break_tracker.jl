using Flux, CuArrays
include("rates_2D.jl")


function invToSpect(Gf, cuE, fscons, NoiseCov, NoiseNMDAratio, NoiseTau)
    #info, invGf = CuArrays.CUBLAS.matinv_batched(Gf)

    invGf = inv(Gf)[2]

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


N = 2
Nthetas = round(Int32, N/2)
rcpt_types = 3
t = 0:0.1:5000
fs = 0:1:100

c = [0 25 50 100]

Jee = param(1.7)
Jei = param(1.525) #param(-0.5)
Jie = param(1.7) #param(2.7)
Jii = param(0.5) #param(-1.5)
i2e = param(0.6)

J0 = [Jee -Jei; Jie -Jii]

# things work though This
v1, conv = steady_Tracker_v(N, rcpt_types, t, c, J0, i2e)

#now showing track_LinApprox code
kk = 0.04
nn = 2
cons = length(c)

fs = fs/1000
r_star =  kk*max.([sum(v1[1:2:end,:], dims=1); sum(v1[2:2:end,:], dims=1)], zeros(N, cons)).^nn
#derivative of f(r_star), or equivalently f'(r_star) the input to Phi.
rs = nn * kk.^(1/nn) * r_star.^(1-1/nn)

t_scale = 1
tauNMDA = 100*t_scale
tauAMPA = 3*t_scale
tauGABA = 5*t_scale
nmdaRatio = 0.1

NoiseNMDAratio = 0
NoiseTau = 1*t_scale
NoiseCov = [1; 1]
tauS = [tauAMPA, tauNMDA, tauGABA]
tauSvec = kron(tauS, ones(1,N)); #vector time-scales
Wtot = [(1-nmdaRatio)*Jee 0 0 0 0 0;
(1-nmdaRatio)*Jie 0 0 0 0 0;
0 0 nmdaRatio*Jee 0 0 0;
0 0 nmdaRatio*Jie 0 0 0;
0 0 0 0 0 -Jei;
0 0 0 0 0 -Jii]

eE = [1; 0]
Phi(rr) = [rr[1] 0; 0 rr[2]]

J = [CuArray(-I + Wtot * kron(ones(rcpt_types, rcpt_types), Phi(rs[:,cc]))) for cc in 1:cons]
eE = kron(ones(rcpt_types), eE)

Gf = [CuArray(-1im * 2 * pi * ff * Diagonal(kron(tauS, ones(N))) - J[cc]) for cc in 1:cons for ff in fs];
cuE = [CuArray(eE) for cc in 1:cons for ff in fs];

fscons = repeat(fs, cons)
SpectE = invToSpect(Gf, cuE, fscons, NoiseCov, NoiseNMDAratio, NoiseTau)
