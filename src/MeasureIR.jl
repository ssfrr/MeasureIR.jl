module MeasureIR

using Compat: @warn
# using Unitful: Time
import Unitful

export stimulus, analyze
export golay

abstract type IRMeasurement end

struct GolaySequence{AT<:AbstractVector} <: IRMeasurement
    A::AT
    B::AT
end

"""
    golay(L)

Create an IR measurement using a complimentary Golay sequence, assuming that the
system being measured has a response less than length `L`. The actual length of
the generated sequence might be greater than L.

## Example

```julia
using Plots: plot
using MeasureIR: golay, stimulus, analyze

meas = golay(4096)

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

function golay(L)
    ispow2(L) || L = nextpow2(L)
    GolaySequence(_golay(L)...)
end

golay(t::Time, samplerate::Frequency) = golay(Int(t*samplerate))
golay(t::Time, samplerate) = golay(Int(t/(1s) * samplerate))

Base.:(==)(g1::GolaySequence, g2::GolaySequence) = g1.A == g2.A && g1.B == g2.B

"""
    stimulus(sig)

Generate the stimulus signal for the given measurement process. This is the
signal that should be run through your test system to measure the impulse
response. After recording the response you can use the `analyze` function to
generate the impulse response.
"""
stimulus(sig::GolaySequence) = [sig.A;
                                zeros(length(sig.A));
                                sig.B;
                                zeros(length(sig.B))]

# naive implementation, allocates a lot, but runs reasonably fast
function _golay(L)
    if L == 1
        [1], [1]
    else
        a, b = _golay(L÷2)
        [a; b], [a; -b]
    end
end

"""
    analyze(sig::IRTestSignal, response::AbstractArray)

Convert a given stimulus response into an impulse response. Generally you will
first create the test signal with e.g. `golay`, then use `stimulus` to give you
the actual stimulus, then pass the original measurement object and the measures
response to `analyze`.

If `response` is a 2D array it is treated as an NxC multichannel response with
`N` frames and `C` channels.
"""
function analyze(sig::GolaySequence, response::AbstractArray)
    # stimuli has 2 measurements, each with a length-L golay sequence
    # and length-L silence.
    L = length(sig.A)
    if size(response, 1) > 4L
        trunc = false
        for ch in 1:size(response, 2)
            if !isapprox(response[(4L+1):end, ch],
                         zeros(size(response, 1)-4L),
                         atol=sqrt(eps()))
                trunc=true
            end
        end
        trunc && @warn "nonzero samples past end of analysis window. Check your test signal is long enough for the response"
    end
    # run the cross-correlation, chopping off the non-causal part and the part
    # that would be corrupted if the IR is too long
    irA = mapslices(response[1:2L, :], 1) do v
        xcorr(v, sig.A)[2L:(3L-1)] ./ 2L
    end
    irB = mapslices(response[(2L+1):4L, :], 1) do v
        xcorr(v, sig.B)[2L:(3L-1)] ./ 2L
    end
    irA + irB
end

end # module
