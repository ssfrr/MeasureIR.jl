using MeasureIR

m = golay()

"""
    chunkminmax(data, nchunks)

Break `data` into `nchunks` chunks, and return a 3-element `Tuple` of
where the 1st item is a `Range` with the indices of each chunk, and the
second and third are the lower and upper bounds of each chunk, respectively.

## Example

```julia
sig = randn(1_000_000) .* cos.(linspace(0,2π, 1_000_000))

i, l, u = chunkminmax(sig, 1000)
plot(i, l, fillrange=u)
```
"""
function chunkminmax(data, nchunks)
    N = size(data, 1)
    C = size(data, 2)
    nchunks < N || return 1:N, data, data
    lower = zeros(nchunks, C)
    upper = zeros(nchunks, C)

    chunksize = floor(Int, N/nchunks)
    @show nchunks
    for i in 1:nchunks
        offset = (i-1) * chunksize
        n = min(chunksize, N-offset)
        @show i
        @show offset
        lower[i, :] = minimum(view(data, (1:n) + offset, :), 1)
        upper[i, :] = maximum(view(data, (1:n) + offset, :), 1)
    end

    (0:(nchunks-1))*chunksize+1, lower, upper
end

function sigplot(sig, sr=1)
    i,l,u = chunkminmax(sig, 1000)
    plot(i/sr, l, fillrange=u, alpha=0.8)
end

chunkminmax(rand(20), 6)


plot(rand(10))


using Interpolations

A = zeros(20)
A = rand(20)
A[10] = 1
itp = interpolate(A, BSpline(Cubic(Line())), OnGrid())

x = linspace(1, 20, 500)
plot(x, itp[x], label="itp")
# plot!(x, x->sinc(x-10), label="sinc")
plot!(x, x->cubic(A, x, true), label="cubic")
scatter!(A, label="A")

# code adapted from http://paulbourke.net/miscellaneous/interpolation/
function cubic(A, i, catmullrom=true)
    y0 = i < 2 ? A[1] : A[floor(Int, i)-1]
    y1 = i < 1 ? A[1] : A[floor(Int, i)]
    y2 = i > length(A) ? A[end] : A[ceil(Int, i)]
    y3 = i > length(A)-1 ? A[end] : A[ceil(Int, i)+1]

    i -= floor(i)
    if catmullrom
        # variant where the derivatives are made to be equal for adjacent
        # regions
        a0 = -0.5*y0 + 1.5*y1 - 1.5*y2 + 0.5*y3
        a1 = y0 - 2.5*y1 + 2*y2 - 0.5*y3
        a2 = -0.5*y0 + 0.5*y2
        a3 = y1
    else
        a0 = y3 - y2 - y0 + y1
        a1 = y0 - y1 - a0
        a2 = y2 - y0
        a3 = y1
    end

    a0*i^3 + a1*i^2 + a2*i + a3
end

png("/home/sfr/Desktop/spline.png")

function measure(measurement)
    # wait so we don't get the sound of the keypress
    sleep(1)
    str = PortAudioStream(2, 2; synced=true)
    sig = stimulus(measurement)
    @async write(str, sig)
    resp = read(str, length(sig))
    close(str)
    resp
end

resp = run(swp)

ir = analyze(swp, resp)
@which analyze(swp, resp)
methods(analyze)

using Plots
plot(resp)
plot()

using Plots, PortAudio, MeasureIR
using Unitful: s, Hz, kHz

swp = expsweep(5*44100; nonlinear=true)
swp = golay(1*44100)
swp = mls(44100)
resp = convert(Array{Float64}, measure(swp))
ir = analyze(swp, resp)
sigplot(resp, 44100)
sigplot(ir, 44100)
size(ir)
plot(linspace(0,22050, size(ir,1)÷2+1), abs.(rfft(ir[:, 2])), legend=false)



testir = [zeros(100); (rand(100)*2-1) .* e.^(linspace(0,-10,100))]
plot(testir)

g = mls(44100)
stim = stimulus(g)
ir = analyze(g, conv(stim, testir))
plot(ir[1:500]*2.2)
plot!(testir)
