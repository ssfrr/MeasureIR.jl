struct Impulse <: IRMeasurement
    len::Int
    prepad::Int
end

prepadding(i::Impulse) = i.prepad

"""
    impulse(samples)
    impulse(seconds::Time, samplerate)
    impulse(seconds::Time, samplerate::Frequency)

Create an impulse response measurement using a literal digital impulse (a single
sample). This is typically not used in practice because the limited energy in
the stimulus signal gives poor signal-to-noise in the measured impulse response.
It is included in this package mostly as a reference.
"""
impulse(L, prepad=L) = Impulse(L, prepad)
impulse(t::Time, samplerate::Frequency) = impulse(Int(t*samplerate))
impulse(t::Time, samplerate) = impulse(Int(t/(1s) * samplerate))

stimulus(m::Impulse) = [zeros(m.prepad); 1.0; zeros(m.len-1)]
analyze(m::Impulse, response::AbstractMatrix) = response[m.prepad+1:end, :]
analyze(m::Impulse, response::AbstractVector) = response[m.prepad+1:end]
