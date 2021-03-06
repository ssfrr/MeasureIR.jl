"""
    truncresponse(response, L, noisepower; warn=true)
Truncate the given array to the given length, printing a warning if the
truncated tail has more average power than the given noise floor.

`noisepower` can be a scalar or a row vector with one entry per channel

Warning can be disabled with the kwarg `warn=false`
"""
function truncresponse(response::AbstractArray, L, noisepower=eps(); warn=true)
    RL = size(response, 1)
    if RL > L
        truncenergy = sum(response[L+1:end, :].^2; dims=1)
        headroom = 1.05 # 5% over noise
        atol = sqrt(eps()) # absolute tolerance (for when noise floor = 0.0)
        if warn && any(truncenergy .> (atol + noisepower * headroom) * (RL-L))
            @warn string("nonsilent samples past end of analysis window. Check ",
                         "your test signal is long enough for the response")
            # Base.show_backtrace(STDERR, backtrace())
        end
        response = timeslice(response, 1:L)
    end

    response
end

# internal functions to do the actual truncation while keeping the
# dimensionality
timeslice(response::AbstractVector, I) = response[I]
timeslice(response::AbstractMatrix, I) = response[I, :]

# used to help parse arguments
striphz(f) = f
striphz(f::Frequency{T}) where T = T(f/Hz)
stripsec(t) = t
stripsec(t::Time{T}) where T = T(t/s)

energy(v) = sum(x->x^2, v)
power(v) = energy(v)/length(v)


"""
Normalize a frequency with a given samplerate to cycles/sample.
"""
freqnorm(f, sr) = uconvert(NoUnits, f/sr)

"""
Normalize a duration with a given samplerate to integer samples (rounded)
"""
durnorm(dur, sr) = round(Int, uconvert(NoUnits, dur*sr))


"""
    snr(mixture, signal)

Returns the signal-to-noise power ratio assuming that the mixture contains a
scaled version of the signal.
"""
function snr(mixture, signal)
    mixcoef = dot(mixture, signal) / dot(signal, signal)
    proj = mixcoef * signal

    energy(proj)/energy(mixture-proj)
end
