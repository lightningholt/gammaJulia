using Statistics

include("rates_2D.jl")
include("linearApprox.jl")

#function to find the loss, or mean square error between the expect/ideal Power Spectra
#and the actual linear approximated power spectra.
# function PS_Loss(Jee, Jei, Jie, Jii, i2e, idealSpect)
function PS_Loss(J0, i2e, idealSpect)
    N = 2 #number of neurons
    rcpt_types = 3
    t = 0:0.1:5000
    fs = 0:1:100
    c = [0 25 50 100]

    if J0[3] > 0
        error("Jei should be negative")
    end

    if J0[4] > 0
        error("Jii should be negative")
    end

    # J0 = [Jee Jei; Jie Jii]

    v_final, Conv = steady_Tracker_v(N, rcpt_types, t, c, J0, i2e)

    # println(v_final)

    SpectE = gammaTrackedLinApprox(N, rcpt_types, fs, v_final, c, J0, i2e)

    SpectE = SpectE./mean(SpectE)

    idealSpect = idealSpect./mean(idealSpect)

    LL = sum(abs.(SpectE - idealSpect))

    return LL, Conv, SpectE
end
