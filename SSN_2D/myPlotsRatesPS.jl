using Plots

#function to plot rates and Power spectra the same way and in the way that I like
#this function takes cc (contrasts) and r_t (the rates in time) as inputs
function myplot_rates(cc, r_t)
    gr()

    r_star = r_t[end,:,:]

    plot(cc', r_star', lc=[:red :blue], label=["Exc", "Inh"], lw=5, legend=:topleft, title="SS firing rates", xlabel="Contrasts", ylabel="Rates (Hz)")
end

#function to plot rates and Power spectra the same way and in the way that I like
#this function takes fs (frequencies) and spect (Power Spectra E) as inputs
function myplot_PS(fs, spect)
    gr()

    plot(fs, spect, lc=[:black :blue :green :red], label=["C = 0", "C = 25", "C = 50", "C = 100"], lw=2.5, legend=:top, xlabel="Frequencies (Hz)", ylabel="Power", title="Power Spectra")
end


function myplot_ps_rates(cc, r_t, fs, spect)
    gr()
    pR = plot(cc', r_star', lc=[:red :blue], label=["Exc", "Inh"], lw=5, legend=:topleft, title="SS firing rates", xlabel="Contrasts", ylabel="Rates (Hz)")
    pPS = plot(fs, spect, lc=[:black :blue :green :red], label=["C = 0", "C = 25", "C = 50", "C = 100"], lw=2.5, legend=:top, xlabel="Frequencies (Hz)", ylabel="Power", title="Power Spectra")
    plot(pR, pPS, layout=(1,2), widths=[0.3, 0.7])
end
