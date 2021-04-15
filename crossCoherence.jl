#==
crossCoherence.jl

calculates the average cross-spectral coherence at the
fundamental frequencies of the two time series. Outputs mean coherence
and the coherence at each fundamental frequency (in cohereStats), as
well as the power/frequency functions for each TS and the cross-spectral
coherence analysis.

Example Syntax:

samplerate = 100
windowSize = 256
windowOverlap = .5

crossCoherence(timeseries.p1_ts, timeseries.p2_ts, samplerate, windowSize, windowOverlap)

Tehran J. Davis, 2020

==#


using DSP, StatsBase, FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots


function crossCoherence(ts1, ts2, samplerate, windowSize, windowOverlap)

    tapering=rectangular

    ## Z-score normalization
    ts1 = StatsBase.zscore(ts1)
    ts2 = StatsBase.zscore(ts2)



    ## FTT parameters
    window = DSP.Windows.hanning(windowSize) #|> Int
    n_overlap = windowSize * windowOverlap |> Int


    ## Calculate the spectrum

    FP1 = DSP.welch_pgram(ts1, windowSize, Int(windowOverlap*windowSize))
    FP2 = DSP.welch_pgram(ts2, windowSize, Int(windowOverlap*windowSize))

    ## Calculate coherence measures
    𝘾 = coherence([ts1 ts2], samplerate, windowSize; tapering=tapering, smoothing=hannSmoother, tril=true)


    x_i = findall(FP1.power[2:end] .== maximum(FP1.power[2:end]))[1]
    y_i = findall(FP2.power[2:end] .== maximum(FP2.power[2:end]))[1]

    coherenceStats = [((𝘾.y[x_i] + 𝘾.y[y_i])/2)[2], 𝘾.y[x_i][2], 𝘾.y[y_i][2]]

    return coherenceStats
end
