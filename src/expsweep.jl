struct ExpSweep{AT} <: IRMeasurement
    sig::AT
    w1::Float64
    w2::Float64
    nonlinear::Bool
end
ExpSweep(sig, w1, w2, nl) = ExpSweep(sig, Float64(w1), Float64(w2), nl)

"""
    expsweep(samples, minfreq=0.0025, maxfreq=π; options...)
    expsweep(seconds::Time, samplerate, minfreq=20, maxfreq=samplerate/2; options...)

Create an impulse response measurement using an exponential sinusoid sweep from
`minfreq` to `maxfreq`, over `samples` samples, following the procedure from
[1]. If no samplerate is given then frequencies are given as rad/sample. The
default minimum frequency of 0.0025 rad/sample corresponds to about 19Hz at
48kHz.

## Options

- `fadein`
- `fadeout`
- `fade`
- `optimize`
- `func`
- `nonlinear`

[1]: Farina, Angelo, Simultaneous measurement of impulse response and distortion
     with a swept-sine technique, Audio Engineering Society Convention 108 (2000).
"""
function expsweep(L, minfreq=0.0025, maxfreq=π;
        # fade times chosen heuristically to be pretty reasonable
        # we could definitely be smarter here, especially for small L
        fadein=round(Int, min(L/2, 200/minfreq)),
        fadeout=round(Int, min(L/4, 800/maxfreq)),
        fade=nothing,
        optimize=true, func=sin,
        nonlinear=false)
    if optimize
        minfreq = _optimizew1(minfreq, maxfreq, L)
    end
    if fade !== nothing
        fadein = fade
        fadeout = fade
    end
    sig = _expsweep(func, L, minfreq, maxfreq)
    winin = 0.5-0.5*cos.(linspace(0, π, fadein))
    winout = 0.5-0.5*cos.(linspace(π, 0, fadeout))
    sig[1:fadein] .*= winin
    sig[(end-fadeout+1):end] .*= winout
    ExpSweep(sig, minfreq, maxfreq, nonlinear)
end

function expsweep(t::Time, samplerate,
                  minfreq=20, maxfreq=striphz(samplerate)/2;
                  options...)
    sr = striphz(samplerate)
    w1 = striphz(minfreq) * 2π / sr
    w2 = striphz(maxfreq) * 2π / sr
    L = Int(stripsec(t) * sr)
    expsweep(L, w1, w2; options...)
end

function _expsweep(func, L, minfreq, maxfreq)
    # switch to notation from [1]
    w1 = minfreq
    w2 = maxfreq
    T = L
    t = 0:(T-1)
    @. func(w1*T /
           log(w2/w1) *
           (exp(t/T * log(w2/w1))-1))
end

"""
Compute a new starting freq w < w1 that gives us an integer number
of total cycles.
"""
function _optimizew1(w1, w2, L)
    # scale to [0-1]
    f1 = w1 * L
    f2 = w2 * L
    # compute the total (integer) number of cycles we want
    M = floor((f2 - f1) /
              (2π*log(f2/f1)))

    a(f) = (f2-f)/log(f2/f) - 2π*M
    fnew = find_zero(a, f1)

    fnew/L
end

stimulus(m::ExpSweep) = [m.sig; zeros(length(m.sig))]
function analyze(m::ExpSweep, response::AbstractArray)
    T = length(m.sig)
    w1 = m.w1
    w2 = m.w2
    t = 0:(T-1)
    # compute the normalizing factor from the derivative of the frequency -
    # the faster the frequency is changing the less time the signal spends
    # on that frequency
    amp = @. w1/T*log(w2/w1)*exp(t/T*log(w2/w1))
    invfilt = amp .* m.sig
    ir = mapslices(response, 1) do v
        xcorr(v, invfilt)
    end

    # keep the full IR (including non-causal parts representing nonlinearities)
    # if m.nonlinear is true
    m.nonlinear ? ir[T:3T, :] : ir[2T:3T, :]
end
