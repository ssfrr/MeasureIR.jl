struct RPMS{AT} <: IRMeasurement
    sig::AT
    prepad::Int
    gain::Float64
end

"""
    rpms(samples; options...)

Generate periodic noise using the random phase multisine medhod [1], by starting
with a flat spectrum and randomizing the phases.

## Options
- `prepad` - Defaults to `samples`, the number of zeros that preceed the stimulus.
- `gain` - Defaults to 0.231439, which scales the stimulus such that 99.99% of
  the analog signal should be below -1dBFS. See `scaling.jl` in the MeasureIR
  repository for details on how this was computed.
"""
function rpms(L; prepad=L, gain=0.231439)
    N = L÷2+1
    X = ones(Complex128, N)
    # set phases to be uniform random in [-π,π]
    @. X *= exp(im*π*(rand()*2-1))
    # DC component should be real
    X[1] = 1
    x = irfft(X, L) * sqrt(L)
    RPMS(x, prepad, gain)
end

# TODO: add trials support
stimulus(s::RPMS) = [zeros(s.prepad); s.sig * s.gain; zeros(length(s.sig))]
prepadding(s::RPMS) = s.prepad

function _analyze(s::RPMS, response::AbstractArray)
    mapslices(response, 1) do v
        xcorr(v[s.prepad+1:end], s.sig)[end÷2+1:end] / length(s.sig) / s.gain
    end
end
