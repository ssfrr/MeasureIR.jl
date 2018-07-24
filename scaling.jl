# Compute the scaling factors for each measurement technique to ensure that
# analog signal is 99.9% below -1dBFS. Random gaussian noise is technically
# unbounded, so it's not useful to find a normalization level that will ensure
# we never clip, but we can at least find a reasonable level where any clipping
# should have a minimal impact.

"""
    percentilethresh(x, p)

Find the threshold where `p` percent of `x` will be below it.
"""
function percentilethresh(x, p)
    xs = sort(abs.(x))
    xs[ceil(Int, p/100*length(x))]
end

# expsweep and impulse are naturally limited, so no need to compute a scaling
# factor for them

using MeasureIR

rawsig(g::MeasureIR.GolaySequence) = [g.A; g.B]
rawsig(g::MeasureIR.IRMeasurement) = g.sig

data = Dict()
for f in [expsweep, golay, mls, schroeder]
    fname = split(string(f),".")[2]
    data[f] = Float64[]
    # we could just use a single long one, but it's also interesting to see
    # how they change at different signal lengths. (try plotting them)
    for N in 2.^(12:20)
        m = f(N)
        thresh = percentilethresh(resample(rawsig(m), 8//1), 99.99)
        println("$fname - $N: $thresh")
        push!(data[f], thresh)
    end
end

using DSP

ref = db2amp(-1)
println("Impulse Gain: ", ref)
println("ExpSweep Gain: ", ref / data[expsweep][end])
println("Golay Gain: ", ref / data[golay][end])
println("MLS Gain: ", ref / data[mls][end])
println("Schroeder Gain: ", ref / data[schroeder][end])
