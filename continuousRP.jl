 using DSP, StatsBase

function ContinuousRP(ts1, ts2)

    # center the data
    ts1 = ts1 .- mean(ts1)
    ts2 = ts2 .- mean(ts2)

    # hilbert transform
    h1 = DSP.hilbert(ts1) |> imag
    h2 = DSP.hilbert(ts2) |> imag
    num = h2.*ts1 - ts2.*h1
    denom = ts2.*ts1 + h2.*h1

    # time series of CRP in radians
    radians = atan.(num,denom)

    # convert to degrees
    pRP = rad2deg.(radians)

    # circular stats in radians
    wgh = ones(size(radians))
    wr = wgh'*exp.(1im*radians)
    meanRP = angle(wr) |> rad2deg
    rvRP = abs(wr)/sum(wgh)
    sdRP = sqrt(2*(1-rvRP))

    (meanRP = meanRP, rvRP = rvRP, sdRP = sdRP, pRP = pRP)

    # Dict("meanRP" => meanRP,
    #     "rvRP" => rvRP,
    #     "sdRP" => sdRP,
    #     "pRP" => pRP)


end
