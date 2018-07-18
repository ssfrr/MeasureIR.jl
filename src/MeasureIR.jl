module MeasureIR

using Compat: @warn
using Unitful: s, Hz, Time, Frequency
using SampledSignals: SampleBuf
using Roots: find_zero
using DSP: hanning, FIRFilter, filt

export stimulus, analyze
export expsweep, golay, impulse

"""
Abstract supertype for impulse response measurements. Subtypes should define a
constructor that allows the user to specify their system's expected decay time.

Subtypes should also define a `stimulus(sig::MyMeasurementType)` function that
returns a vector suitable for playing directly into the system, as well as a
`analyze(sig::MyMeasurementType, response)` that takes the system's response
and generates an impulse response.
"""
abstract type IRMeasurement end

"""
    stimulus(sig)

Generate the stimulus signal for the given measurement process. This is the
signal that should be run through your test system to measure the impulse
response. After recording the response you can use the `analyze` function to
generate the impulse response.
"""
function stimulus end

"""
    analyze(sig::IRMeasurement, response::AbstractArray)

Convert a given stimulus response into an impulse response. Generally you will
first create the test signal with e.g. `golay`, then use `stimulus` to give you
the actual stimulus, then pass the original measurement object and the measures
response to `analyze`.

If `response` is a 2D array it is treated as an NxC multichannel response with
`N` frames and `C` channels.
"""
function analyze end

# used to help parse arguments
striphz(f) = f
striphz(f::Frequency{T}) where T = T(f/Hz)
stripsec(t) = t
stripsec(t::Time{T}) where T = T(t/s)

include("golay.jl")
include("impulse.jl")
include("expsweep.jl")

end # module
