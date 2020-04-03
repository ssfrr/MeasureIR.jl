"""
Estimate the period of `x`, from an initial guess `L`. The result will be the
peak of the continuous cross-correlation between `0.5L` and `1.5L`, i.e. the
period is not constrained to be an integer.

The algorithm is to first compute a discrete cross-correlation and oversample by
4x to account for inter-sample peaks. We use the peak of the oversampled
cross-correlation as a starting point for a numerical optimization of the
average-squared-distance function in the frequency domain.
"""
function getperiod(x::AbstractVector, L; max_ppm=500_000)
    maxD = L * max_ppm ÷ 1_000_000

    # create two signals that are offset by L (our expected period length), and
    # zero-pad by to do a linear (rather than circular) convolution. We set up
    # the padding to pre-shift so the cross-correlation starts start at their
    # left-most lag (`-maxD`), and the `maxD`th sample will actually be L lag)
    # we also add extra padding to account for time aliasing from the fractional
    # delay. Otherwise it can wrap around. 512 samples is enough for a sinc
    # function to decay more than 60dB

    # note this could be problematic if the signal is temporally compact, e.g.
    # an impulse train. Innacuracy in `L` could chop off an impulse from one or
    # the other of these.
    x1 = [zeros(maxD); x[L+1:end]; zeros(512)]
    x2 = [x[1:(end-L)]; zeros(maxD); zeros(512)]

    m = length(x1)
    X1 = rfft(x1)
    X2 = rfft(x2)
    M = length(X1)

    # do a discrete cross-correlation and upsample. The maximum here should be
    # a pretty good estimate of the true peak. Should be valid from 1:8L, which
    # represents lags from -L:L (upsampled by 4)
    # not sure why `xc` wasn't being type-inferrred`
    xc = resample(irfft(X1 .* conj.(X2), m)[1:2maxD+1], 4)::Vector{Float64}
    D_est1 = (findmax(xc)[2]-1)/4-maxD

    # note if `m` were guaranteed to be even the upper bound here would be π
    ω = range(0, (M-1)*2π/m, length=M)

    # error function for  frequency-domain delay estimation. This can be
    # interpreted as a frequency-domain ASDF (average-squared distance function)
    # as a function of the continuous delay
    err(D) = sum(abs2(X1[i] - X2[i] * exp(-im*ω[i]*(D+maxD)))
                 for i in 1:M)

    # opt = optimize(D->err(first(D)), [D_est1], NewtonTrustRegion(), autodiff=:forward)
    opt = optimize(err, D_est1-0.5, D_est1+0.5)
    if !converged(opt)
        throw(ErrorException("Delay estimation failed to converge"))
    end

    # return the delay and the normalized minimum error
    norm = 2sqrt(sum(abs2, X1)*sum(abs2, X2))
    L+minimizer(opt)[1], minimum(opt)/norm
end
