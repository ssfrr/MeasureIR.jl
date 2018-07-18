struct Impulse <: IRMeasurement
    len::Int
end

"""
    impulse(samples)
    impulse(seconds::Time, samplerate)
    impulse(seconds::Time, samplerate::Frequency)

Create an impulse response measurement using a literal digital impulse (a single
sample). This is typically not used in practice because the limited energy in
the stimulus signal gives poor signal-to-noise in the measured impulse response.
It is included in this package mostly as a reference.
"""
impulse(L) = Impulse(L)
impulse(t::Time, samplerate::Frequency) = impulse(Int(t*samplerate))
impulse(t::Time, samplerate) = impulse(Int(t/(1s) * samplerate))

stimulus(m::Impulse) = [1.0; zeros(m.len-1)]
analyze(m::Impulse, response::AbstractArray) = response
