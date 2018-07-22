struct GolaySequence{AT<:AbstractVector} <: IRMeasurement
    A::AT
    B::AT
    prepad::Int
    amp::Float64
    upsample::Int
end

prepadding(g::GolaySequence) = g.prepad

"""
    golay(samples; options...)
    golay(seconds::Time, samplerate; options...)
    golay(seconds::Time, samplerate::Frequency; options...)

Create an IR measurement using a complimentary Golay sequence, assuming that the
system being measured has a response less than length `samples`. The actual
length of the generated stimulus sequence might be greater than `samples`.

You can also specify the length with a duration in seconds along with a
sampling rate. The sampling rate can either be specified as a raw number or
a unitful frequency.

## Options

- `amp::AbstractFloat`: Amplitude of the binary signal. Defaults to 1/2.2, which
  provides enough headroom so that the analog stimulus has an maximum amplitude
  of about 1. This means that if you set the gain of your playback system such
  that a full-scale sinusoid does not distort, the default golay code should not
  distort either.
- `prepad::Int`: Defaults to `samples`. Puts a period of silence at the
  beginning of the measurement that can be used to estimate the noise floor.
- `upsample::Integer`: Defaults to 1, which creates a standard binary golay code
  with energy up to nyquist. Setting this to a higher number will create a
  bandlimited version, so `upsample=2` will only have energy up to 1/2 the
  nyquist frequency. Often the upper frequencies of the impulse response are
  time-variant, so they may not be accurately measurable anyways. Decreasing the
  bandwidth can also help avoid nonlinearities due to slew-rate limiting in the
  playback system.

## Example

```julia
using Plots: plot
using MeasureIR: golay, stimulus, analyze
using Unitful: s, Hz

meas = golay(2s, 48kHz)

# generate the test stimuli suitable for playback
stim = stimulus(meas)

# create a synthetic impulse response and simulate the system. This is where
# you'd normally play the stimuli through your system and record the response
irsim = 1./exp.(0:0.1:9.9) .* (rand(100) .- 0.5)
output = conv(stim, irsim)

# analyze to reconstruct the impulse response
ir = analyze(meas, output)

plot([irsim[1:100], ir[1:100]], labels=["Convolved IR", "Measured IR"])
```
"""
function golay end

function golay(L; amp=1/2.2, prepad=L, upsample=1)
    L = nextpow2(ceil(Int, L/upsample))
    A, B = _golay(L)
    GolaySequence(A, B, prepad, amp, upsample)
end

golay(t::Time, samplerate::Frequency; options...) = golay(Int(t*samplerate); options...)
golay(t::Time, samplerate; options...) = golay(Int(t/(1s) * samplerate); options...)

Base.:(==)(g1::GolaySequence, g2::GolaySequence) = g1.A == g2.A && g1.B == g2.B

function stimulus(sig::GolaySequence)
    stim = [sig.amp .* sig.A;
            zeros(length(sig.A));
            sig.amp .* sig.B;
            zeros(length(sig.B))]
    if sig.upsample != 1
        # force upsampling ratio to be a rational (DSP.jl issue #211)
        stim = filt(FIRFilter(sig.upsample//1), stim)
    end

    [zeros(sig.prepad); stim]
end

# naive implementation, allocates a lot, but runs reasonably fast
function _golay(L)
    if L == 1
        [1.0], [1.0]
    else
        a, b = _golay(L÷2)
        [a; b], [a; -b]
    end
end

function analyze(sig::GolaySequence, response::AbstractArray)
    # stimuli has 2 measurements, each with a length-L golay sequence
    # and length-L silence.
    L = Int(length(sig.A) * sig.upsample)
    np = noisefloor(sig, response)
    # we chop off the silence part - there's nothing causal for us there
    respA = @views response[sig.prepad+(1:2L), :]
    @views respB = truncresponse(response[sig.prepad+2L+1:end, :], 2L, np)
    if sig.upsample != 1
        # we can use the full-bandwidth zero-stuffed signal for the
        # cross-correlation, so that any weirdness caused by our upsampling
        # procedure only affects us once, and we get a nice bandlimited
        # impulse in the pure loopback case. There's an argument for using the
        # bandlimited version in that frequency content above the cutoff must
        # have been caused by nonlinearity, but for now we'll go for the clean
        # impulse. There's also some kind of strange issue where there ends up
        # being a bunch of pre-ringing garbage when using the bandlimited
        # version. Not sure what that's about. Also this way it doesn't cut off
        # the beginning of the bandlimited impulse when there's no delay.
        A = zerostuff(sig.A, Int(sig.upsample)) * sig.upsample
        B = zerostuff(sig.B, Int(sig.upsample)) * sig.upsample
        # A = filt(FIRFilter(sig.upsample), sig.A)
        # B = filt(FIRFilter(sig.upsample), sig.B)
    else
        A = sig.A
        B = sig.B
    end
    # run the cross-correlation, chopping off the non-causal part and the part
    # that would be corrupted if the IR is too long
    irA = mapslices(respA, 1) do v
        xcorr(v, A)[2L:(3L-1)] ./ (2L * sig.amp)
    end
    irB = mapslices(respB, 1) do v
        xcorr(v, B)[2L:(3L-1)] ./ (2L * sig.amp)
    end
    irA + irB
end

# this is no longer used because the brickwall cutoff causes unacceptable
# ringing in the time-domain signal, causing aliasing in the deconvolution
"""
upsamples/interpolates `A` treating it as periodic.
"""
function circ_upsample(A, ratio)
    L = length(A)
    L2 = L*ratio
    zpad = L2÷2 - L÷2
    spec = [rfft(A); zeros(Complex128, zpad)]
    if iseven(L)
        # in the original spectrum the nyquist bin was doubled up because of
        # periodicity. Now the energy is split between two bins.
        spec[L÷2+1] *= 0.5
    end

    irfft(spec, L2) * ratio
end

"""
upsample by zero-stuffing in between samples
"""
function zerostuff(A, ratio)
    A2 = zeros(length(A)*ratio)
    A2[1:ratio:end] .= A
    A2
end
