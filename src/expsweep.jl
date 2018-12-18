struct ExpSweep{AT, SR} <: IRMeasurement
    sig::AT
    w1::Float64
    w2::Float64
    gain::Float64
    samplerate::SR
end

# the default constructor doesn't auto-convert the non-parametric arguments
ExpSweep(sig::T, w1, w2, g, sr::SR) where {T, SR} = ExpSweep{T, SR}(sig, w1, w2, g, sr)

"""
    expsweep(length; options...)

Create an impulse response measurement using an exponential sinusoid sweep from
`minfreq` to `maxfreq`, following the procedure from [1]. Durations and
frequencies can be given with or without Unitful units. With the default
`samplerate` of `1` durations can be interpreted as samples and frequencies are
in cycles/sample. The default minimum frequency of 0.0004 cycles/sample
corresponds to about 19Hz at 48kHz.

## Options

- `samplerate=1`
- `minfreq=0.0004`
- `maxfreq=0.5`
$optiondoc_gain
- `fadein`
- `fadeout`
- `fade`
- `optimize`
- `func`

## Examples

```julia
using MeasureIR

# generate a 10000-sample sweep
expsweep(10000)

# generate a 10s sweep at 48kHz samplerate
using Unitful: s, kHz
expsweep(10s; samplerate=48kHz)
```

[1]: Farina, Angelo, Simultaneous measurement of impulse response and distortion
with a swept-sine technique, Audio Engineering Society Convention 108 (2000).
"""
# TODO: review Novak: https://ant-novak.com/posts/research/2015-10-30_JAES_Swept/
# it was recommended by Gordon
function expsweep(length;
        samplerate = 1,
        minfreq=0.0004*samplerate,
        maxfreq=0.5*samplerate,
        # default fades are computed below
        fadein=nothing,
        fadeout=nothing,
        fade=nothing,
        gain=expsweep_gain,
        # for long sweeps we don't need long silence
        prepad=0,
        optimize=true,
        func=sin)

    w1 = uconvert(NoUnits, minfreq * 2π / samplerate)
    w2 = uconvert(NoUnits, maxfreq * 2π / samplerate)
    L = round(Int, uconvert(NoUnits, length*samplerate))

    if fade !== nothing
        if fadein !== nothing || fadeout !== nothing
            throw(ArgumentError("`fade` argument cannot be used with `fadein` or `fadeout`"))
        end
        fadein = fade
        fadeout = fade
    end

    # default fade times chosen heuristically to be pretty reasonable
    # we could definitely be smarter here, especially for small L
    if fadein === nothing
        f1 = round(Int, min(L/2, 200/w1))
    else
        f1 = round(Int, uconvert(NoUnits, fadein*samplerate))
    end

    if fadeout === nothing
        f2 = round(Int, min(L/4, 800/w2))
    else
        f2 = round(Int, uconvert(NoUnits, fadeout*samplerate))
    end

    if optimize
        w1 = _optimizew1(w1, w2, L)
    end
    sig = _expsweep(func, L, w1, w2)
    winin = 0.5 .- 0.5*cos.(range(0, π; length=f1))
    winout = 0.5 .- 0.5*cos.(range(π, 0; length=f2))
    sig[1:f1] .*= winin
    sig[(end-f2+1):end] .*= winout
    ExpSweep(sig, w1, w2, gain, samplerate)
end

function _expsweep(func, L, minfreq, maxfreq)
    # switch to notation from [1]
    w1 = minfreq
    w2 = maxfreq
    T = L
    t = 0:(T-1)
    @. func(w1*T /
           log(w2/w1) * (exp(t/T * log(w2/w1))-1))
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

# SampledSignals doesn't accept unitful values for samplerates, so strip it
stimulus(m::ExpSweep{AT, <:Frequency}) where AT = SampleBuf(m.sig * m.gain,
                                       uconvert(NoUnits, m.samplerate*s))
stimulus(m::ExpSweep{AT, <:Real}) where AT = SampleBuf(m.sig * m.gain, m.samplerate)

"""
    function analyze(m::ExpSweep, response::AbstractArray; noncausal=false)

When the measurement is a sweep, you can pass the `noncausal=true` keyword
argument to include the noncausal part of the impulse response, which carries
energy from nonlinear distortion in the system. In this case the zero-lag part
of the impulse response starts at index `L÷2+1`, where `L` is the length of the
IR.
"""
function _analyze(m::ExpSweep, response::AbstractArray; noncausal=false)
    # we don't truncate the response here because typically sweep measurements
    # are taken without repeating, so there's no real errors we'd catch by
    # truncating.
    T = length(m.sig)
    w1 = m.w1
    w2 = m.w2
    t = 0:(T-1)
    # compute the normalizing factor from the derivative of the frequency -
    # the faster the frequency is changing the less time the signal spends
    # at that frequency so we need to boost the amplitude to compensate
    amp = @. w1/T*log(w2/w1)*exp(t/T*log(w2/w1))
    invfilt = amp .* m.sig
    # calculate an extra scaling factor to get the center of the passband to
    # unity. At some point it might make sense to look at a range instead of
    # just the center frequency, in case we catch it in the middle of some
    # ripple. Generally things should be pretty flat if there's a fade in/out
    # TODO: seems kinda silly to compute the whole FFT just to look at the
    # center frequency. We could just do the dot product of that freq
    roundtrip = xcorr(m.sig * m.gain, invfilt)
    centeridx = round(Int, (w1 + (w2-w1)/2)/π*(length(roundtrip)/2))
    scale = abs(rfft(roundtrip)[centeridx])

    ir = mapslices(response, dims=1) do v
        xcorr(v, invfilt) ./ scale
    end

    # keep the full IR (including non-causal parts representing nonlinearities)
    # if m.nonlinear is true
    zeroidx = size(ir,1)÷2+1
    if noncausal
        timeslice(ir, zeroidx-T:zeroidx+T)
    else
        # nonlinearity detection is not working right now, but it seems like a
        # useful feature to add in the future. One issue is that we need a way
        # to estimate the expected "noise floor" of the noncausal part of the
        # IR, which is _not_ the same as the noise floor of the response. It
        # definitely helps to deconvolve just the silence part and use that to
        # estimate the IR noise floor. The other issue is that the overall
        # energy in the noncausal part just doesn't change that much, because
        # all the extra energy is pretty concentrated and it gets averaged out.
        # probably to get this to work well we either need to explicitly look
        # for peaks, and possibly look at fixed delays from the main peak that
        # correspond to harmonics.
        # deconvsilence = mapslices(resp, 1) do v
        #     xcorr(v[1:prepadding(m)], invfilt)
        # end
        # nf = noisefloor(m, deconvsilence)

        # we don't check all the way up to the center because if the impulse is
        # bandlimited it'll have energy directly preceeding it. In reality
        # how far back we should go depends on the bandwidth of the response,
        # so this number may need to be tuned.
        # fortunately this is only to heuristically decide whether we print a
        # warning, so we're never going to get to 100% accuracy anyways...
        # preimpulse = 100
        # noncausalpower = sum(ir[zeroidx-T:zeroidx-preimpulse].^2) /
        #                  (2T+1-preimpulse)
        # if any(noncausalpower .> nf .* 1.1 .+ sqrt(eps()))
        #     @warn "Energy in noncausal IR is above noise floor. Check for nonlinearity"
        # end
        timeslice(ir, zeroidx:zeroidx+T)
    end
end
