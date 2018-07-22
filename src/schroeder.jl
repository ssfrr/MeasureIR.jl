struct SchroederNoise{AT} <: IRMeasurement
    sig::AT
    prepad::Int
end

"""
Genereate schroeder noise by starting with a flat spectrum and
randomizing the phases.
"""
function schroeder(L; prepad=L)
    # TODO: check amplitude scaling - seems very small
    N = L÷2+1
    X = ones(Complex128, N)
    # set phases to be uniform random in [-π,π]
    @. X *= exp(im*π*(rand()*2-1))
    # DC component should be real
    X[1] = 1
    x = irfft(X, L)
    SchroederNoise(x, prepad)
end

# TODO: add trials support
stimulus(s::SchroederNoise) = [zeros(s.prepad); s.sig]
prepadding(s::SchroederNoise) = s.prepad

function _analyze(s::SchroederNoise, response::AbstractArray)
    L = length(s.sig)
    mapslices(response, 1) do v
        xcorr(v[s.prepad+1:end], s.sig)[end÷2+1:end]
    end
end
