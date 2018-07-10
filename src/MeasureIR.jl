module MeasureIR

using Compat: @warn
export analyze, generate, golay

abstract type IRTestSignal end

struct GolaySequence{AT<:AbstractVector} <: IRTestSignal
    A::AT
    B::AT
end

"""
    golay(L)

Create a test signal using a complimentary Golay sequence, assuming that the
system being measured has a response less than length `L`. `L` should be a power
of two.

## Example

```julia
sig = golay(4096)

# generate the actual test signal suitable for playback
output = generate(sig)

# create a synthetic impulse response
ir = 1./exp.(0:100) .* rand(100)

analyze(sig, conv(output, ir)) ≈ ir # should be true
```
"""
function golay(L)
    if !ispow2(L)
        L2 = nextpow2(L)
        @warn "golay($L): should be power of two, upgrading to $L2"
        L = L2
    end
    GolaySequence(_golay(L)...)
end

"""
    generate(sig)

Generate the stimulus signal for the given measurement process. This is the
signal that should be run through your test system to measure the impulse
response. After recording the response you can use the `analyze` function to
generate the impulse response.
"""
generate(sig::GolaySequence) = [sig.A;
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
first create the test signal with e.g. `golay`, then use `generate` to give you
the actual stimulus, then pass the original stimulus and the measures response
to `analyze`.

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
