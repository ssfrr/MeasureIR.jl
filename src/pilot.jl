function pilot(length, freq; samplerate=1)
    L = round(Int, uconvert(NoUnits, length*samplerate))
    dθ = uconvert(NoUnits, 2π*freq/samplerate)

    cos.((0:L-1) * dθ)
end

function pilotfilt(sig, ω)
    # DSP.jl normalizes to half-cycle/sample, so that 1 is the nyquist,
    # but the rest of this code uses 1 cycle/sample (so 0.5 is nyquist)
    bpf = digitalfilter(Bandpass(2ω*0.97, 2ω/0.97),
                        FIRWindow(;transitionwidth=2ω*0.1))
    filt(bpf, sig)[length(bpf):end]
end

"""
Compute reasonable onset and offset thresholds to find a pulse of length L
within `sig`.
"""
function pulsethreshs(sig, L)
    if length(sig) < 4L
        @warn "`pulsethreshs` assumes `length(sig) > 4L` to estimate noise statistics"
    end
    lfrac = L/length(sig)
    qs = quantile(sig, [0.75, 1-lfrac/2])
    qs[1] .+ (0.55, 0.45) .* (qs[2]-qs[1])
end

"""
Find a pilot tone in the given signal. Note that the total signal length should
be at least 4 times longer than the expected length of the tone, so we can
estimate noise statistics.
"""
function findpilot(sig, freq, len; samplerate=1)
    ω = freqnorm(freq, samplerate)
    L = durnorm(len, samplerate)
    filtered = pilotfilt(sig, ω)
    env = abs2.(hilbert(filtered))
    onthresh, offthresh = pulsethreshs(env, L)
    # lpf = gaussian(round(Int, 0.25*L), 0.15)
    lpf = digitalfilter(Lowpass(1/4L),
                        FIRWindow(;transitionwidth=30/L))
    smoothed = filt(lpf, env)
    onset = findfirst((x->x>onthresh), smoothed)
    onset === nothing && return nothing, nothing
    offset = onset+findfirst(x->x<offthresh, smoothed[onset+1:end])
    filtdelay = length(lpf)÷2
    onset -= filtdelay
    offset === nothing && return onset, nothing
    offset -= filtdelay
    if !isapprox(offset-onset, L; rtol=0.5)
        @warn "Found pilot with length $(offset-onset) instead of $L"
    end

    # only pull the middle 80% of the tone to account for any other onset/offset
    # irregularities
    slack = 0.05(offset-onset)
    round.(Int, (onset+slack, offset-slack))
end

"""
Compute the average phase difference of `sig` (which should correspond to the
frequency). The frequency is given in cycles/sample
"""
function avgphasediff(sig)
    # we remove the ends to get rid of the ringing from the hilbert transform
    an = hilbert(sig)[2000:end-2000]
    diff = @view(an[2:end]) ./ @view an[1:end-1]
    mean(angle, diff)/2π
end

"""
Takes the given signal and pilot frequency and returns the resample factor
necessary to synchronize the signal.
"""
function pilotsync(sig, freq; samplerate=1)
    ω = freqnorm(freq, samplerate)
    # bpf = gaussian(N, 0.20) .* cos.((0:N-1) * ω)
    # filtered = filt(bpf, sig)[length(bpf):end]
    filtered = pilotfilt(sig, ω)
    ω̂₁ = avgphasediff(filtered)
    @info "frequency estimated from hilbert transform: $(ω̂₁)"

    # estimate amplitude, frequency, and phase
    model(n, p) = @. p[1] * cos(2π*p[2]*n + p[3])

    fit = curve_fit(model, (0:length(filtered)-1), filtered,
                    [1.0, ω̂₁, 0.0], autodiff=:forwarddiff)

    # @info "frequency estimated from lsq fit: $(fit.param[2]/2π*samplerate)"
    @info "fit amplitude: $(fit.param[1])"
    @info "fit frequency: $(fit.param[2] * samplerate)"
    @info "fit phase (rad): $(fit.param[3])"

    fit.param[2] / ω
end


"""
    phaselock(x; alpha=0.05)

Simple PLL (phase-locked loop). Assumes `x` to be a real or complex sinusoidal
signal. Returns a tuple of (freq, phase, err, x2), where each is a vector with
the same length as `x`. `x2` will be real or complex to match `x`.

`alpha` controls how agressive the PLL is in correcting its frequency and phase
estimates. Lower values are more stable but converge more slowly.

## Example

```
using Plots

in_p = 0:0.1:500
x = cos.(in_p) + randn(length(in_p))*0.3

f, p, err, x2 = phaselock(x)
```

Algorithm ported from C code at http://liquidsdr.org/blog/pll-simple-howto/
Released under the MIT License, Copyright (c) 2007 - 2016 Joseph Gaeddert
"""
function phaselock(x::Vector{<:Complex}; alpha=0.05)
    beta = 0.5*alpha^2
    freq = 0.0;
    phase = 0.0;
    freqs = zeros(length(x))
    phases = zeros(length(x))
    errs = zeros(length(x))
    out = zeros(Complex{Float64}, length(x))

    for i in eachindex(x)
        out[i] = exp(im*phase)
        err = angle(x[i]*conj(out[i]))

        freq += beta * err
        phase += freq + alpha * err

        freqs[i] = freq/2π
        phases[i] = phase/2π
        errs[i] = err
    end

    (freqs, phases, errs, out)
end

function phaselock(x::Vector{<:Real}; kwargs...)
    freqs, phases, errs, out = phaselock(hilbert(x); kwargs...)

    (freqs, phases, errs, real.(out))
end
