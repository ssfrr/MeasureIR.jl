struct SchroederNoise{AT} <: IRMeasurement
    sig::AT
end

"""
Genereate schroeder noise by starting with a flat spectrum and
randomizing the phases.
"""
function schroeder(L)
    N = L÷2+1
    X = ones(Complex128, N)
    # set phases to be uniform random in [-π,π]
    @. X *= exp(im*π*(rand()*2-1))
    # DC component should be real
    X[1] = 1
    x = irfft(X, L)
    SchroederNoise(x)
end

stimulus(s::SchroederNoise) = s.sig

function analyze(s::SchroederNoise, response::AbstractArray)
    L = length(s.sig)
    mapslices(response, 1) do v
        xcorr(v, s.sig)[L:end]
    end
end
