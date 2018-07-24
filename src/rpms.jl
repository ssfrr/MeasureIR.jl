struct RPMS{AT} <: IRMeasurement
    sig::AT
    prepad::Int
    gain::Float64
end

"""
    rpms(samples; options...)

Generate periodic noise using the random phase multisine medhod [1], by starting
with a flat spectrum and randomizing the phases. Note that this is,
nondeterministic, so make sure you use the same measurement instance to generate
and analyze your stimulus.

## Options
- `prepad` - Defaults to `samples`, the number of zeros that preceed the stimulus.
- `gain` - Defaults to 0.231439, which scales the stimulus such that 99.99% of
  the analog signal should be below -1dBFS. See `scaling.jl` in the MeasureIR
  repository for details on how this was computed.

[1]: Mateljan, I., Signal selection for the room acoustics measurement, 1999,
IEEE Workshop on Applications of Signal Processing to Audio and Acoustics,
pp. 199–202.
"""
function rpms(L; prepad=L, gain=0.231439)
    N = L÷2+1
    X = ones(Complex128, N)
    # set phases to be uniform random in [0,2π]
    @. X *= exp(im*2π*rand())
    # DC component should be zero
    X[1] = 0.0
    if iseven(L)
        # Nyquist component should be real (it needs to equal its complex
        # conjugate)
        X[N] = 1.0
    end
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
