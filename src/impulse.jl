struct Impulse <: IRMeasurement
    len::Int
    prepad::Int
    gain::Float64
end

prepadding(i::Impulse) = i.prepad

"""
    impulse(samples; options...)
    impulse(seconds::Time, samplerate; options...)
    impulse(seconds::Time, samplerate::Frequency; options...)

Create an impulse response measurement using a literal digital impulse (a single
sample). This is typically not used in practice because the limited energy in
the stimulus signal gives poor signal-to-noise in the measured impulse response.
It is included in this package mostly as a reference.

## Options
- `gain` - Defaults to 0.89125 (-1dB). Set the gain applied to the stimulus.
- `prepad` - Defaults to `samples`. The number of samples of silence that
  preceeds the stimulus.
"""
impulse(L; prepad=L, gain=db2amp(-1)) = Impulse(L, prepad, gain)
impulse(t::Time, samplerate::Frequency; kwargs...) = impulse(Int(t*samplerate); kwargs...)
impulse(t::Time, samplerate; kwargs...) = impulse(Int(t/(1s) * samplerate); kwargs...)

stimulus(m::Impulse) = [zeros(m.prepad); m.gain; zeros(m.len-1)]
_analyze(m::Impulse, response::AbstractArray) =
        timeslice(response, m.prepad+1:size(response, 1)) ./ m.gain
