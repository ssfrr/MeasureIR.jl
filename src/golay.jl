# TODO: we could get away with removing samplerate from here and storing
# the A and B signals as SampleBufs
struct GolaySequence{AT<:AbstractVector, SR} <: IRMeasurement
    A::AT
    B::AT
    decay::Int
    gain::Float64
    bandlimit::Int
    samplerate::SR
end

GolaySequence(A::AT, B::AT, decay, gain, bandlimit, sr::SR) where {AT, SR} =
    GolaySequence{AT, SR}(A, B, decay, gain, bandlimit, sr)

"""
    golay(length; options...)

Create an IR measurement using a complimentary Golay sequence following the
method from [1], assuming that the system being measured has a response less
than length `samples`. The actual length of the generated stimulus sequence
might be greater than `samples`.

You can also specify the length with a duration in seconds along with a
sampling rate. The sampling rate can either be specified as a raw number or
a unitful frequency.

TODO: update length docs

## Options

$optiondoc_gain
- `bandlimit::Integer`: Defaults to 1, which creates a standard binary golay code
  with energy up to nyquist. Setting this to a higher number will create a
  bandlimited version, so `bandlimit=2` will only have energy up to 1/2 the
  nyquist frequency. Often the upper frequencies of the impulse response are
  time-variant, so they may not be accurately measurable anyways. Decreasing the
  bandwidth can also help avoid nonlinearities due to slew-rate limiting in the
  playback system.
- `decay`: Specify how long you expect the system decay length to be. The
  simulus will wait this long between the A and B parts. Note that this should
  include any delay introduced by the system or measurement process.

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


[1]: Berdahl, Edgar J. and Julius O. Smith III, "Transfer Function Measurement
Toolbox", REALSIMPLE Project. Released 2008-06-05 under the Creative Commons
License (Attribution 2.5) Center for Computer Research in Music and Acoustics
(CCRMA), Stanford University. Accessed 2018 from
https://ccrma.stanford.edu/realsimple/imp_meas/Golay_Code_Theory.html
"""
function golay end

function golay(length; decaylength=length, samplerate=1, gain=golay_gain, bandlimit=1)
    L = ceil(Int, uconvert(NoUnits, length*samplerate/bandlimit))
    dec = ceil(Int, uconvert(NoUnits, decaylength*samplerate))
    L = nextpow(2, L)
    A, B = _golay(L)
    GolaySequence(A, B, dec, gain, bandlimit, samplerate)
end

Base.:(==)(g1::GolaySequence, g2::GolaySequence) = g1.A == g2.A && g1.B == g2.B

function stimulus(sig::GolaySequence)
    stim = [sig.gain .* sig.A;
            zeros(sig.decay);
            sig.gain .* sig.B]
    if sig.bandlimit != 1
        # force upsampling ratio to be a rational (DSP.jl issue #211)
        # TODO: this is a bit broken because it doesn't handle the warm-up
        # in the convolution. but `resample` doesn't give you exactly the
        # right number of samples, so this needs to be looked into.
        # TODO: also I need to think about how we're handling the decay time
        # here, i.e. whether it shuold be before of after the resampling
        stim = filt(FIRFilter(sig.bandlimit//1), stim)
    end

    stim
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

function _analyze(sig::GolaySequence, response::AbstractArray)
    # stimuli has 2 measurements, each with a length-L golay sequence. They are
    # separated by `sig.decay` samples of silence.
    L = Int(length(sig.A) * sig.bandlimit)
    # np = noisefloor(sig, response)
    # we chop off the silence part - there's nothing causal for us there
    respA = timeslice(response, 1 : L+sig.decay)
    respB = timeslice(response, L+sig.decay+1 : size(response, 1))
    # @views respB = truncresponse(
    #         timeslice(response, L+sig.decay+1 : size(response, 1)), 2L, np)
    if sig.bandlimit != 1
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
        A = zerostuff(sig.A, Int(sig.bandlimit)) * sig.bandlimit
        B = zerostuff(sig.B, Int(sig.bandlimit)) * sig.bandlimit
        # A = filt(FIRFilter(sig.bandlimit), sig.A)
        # B = filt(FIRFilter(sig.bandlimit), sig.B)
    else
        A = sig.A
        B = sig.B
    end
    # run the cross-correlation, chopping off the non-causal part and the part
    # that would be corrupted if the IR is too long
    irA = mapslices(respA, dims=1) do v
        # warning - xcorr pads the inputs to equal length, so the total result
        # length should be 2*maxlen-1
        xc = xcorr(v, A) / (2L * sig.gain)
        time0 = (length(xc)+1)÷2-1
        xc[time0 .+ (1:sig.decay)]
    end
    irB = mapslices(respB, dims=1) do v
        xc = xcorr(v, B) / (2L * sig.gain)
        time0 = (length(xc)+1)÷2-1
        xc[time0 .+ (1:sig.decay)]
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
