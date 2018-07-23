struct SchroederNoise{AT} <: IRMeasurement
    sig::AT
    gain::Float64
    prepad::Int
end

"""
Genereate schroeder noise by starting with a flat spectrum and
randomizing the phases.
"""
function schroeder(L; prepad=L)
    # we use our own RNG and seed every time so the result is deterministic.
    mt = MersenneTwister(0)
    N = L÷2+1
    X = ones(Complex128, N)
    # set phases to be uniform random in [-π,π]
    @. X *= exp(im*π*(rand(mt)*2-1))
    # DC component should be real
    X[1] = 1
    x = irfft(X, L)
    gain = 1/maximum(abs.(x))
    SchroederNoise(x, gain, prepad)
end

# TODO: add trials support
stimulus(s::SchroederNoise) = [zeros(s.prepad); s.sig * s.gain; zeros(length(s.sig))]
prepadding(s::SchroederNoise) = s.prepad

function _analyze(s::SchroederNoise, response::AbstractArray)
    mapslices(response, 1) do v
        xcorr(v[s.prepad+1:end], s.sig)[end÷2+1:end] / s.gain
    end
end
