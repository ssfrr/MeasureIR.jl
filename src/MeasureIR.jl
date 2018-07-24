module MeasureIR

using Compat: @warn
using Unitful: s, Hz, Time, Frequency
using SampledSignals: SampleBuf, samplerate
using Roots: find_zero
using DSP: hanning, FIRFilter, filt, db2amp

export stimulus, analyze, noisefloor, prepadding, snr
export expsweep, golay, mls, rpms, impulse

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
    prepadding(meas::IRMeasurement)

Get the number of samples of prepadding silence for the given measurement.
"""
function prepadding end

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

"""
    noisefloor(sig::IRMeasurement, response::AbstractArray)

Estimate the noise floor in the given signal using the prepadded silence
specified by the measurement. Throws an error if the measurement has no
prepadding.

For a multichannel response returns a row vector with power per channel, or
returns a single number if the input is a Vector.
"""
function noisefloor(sig::IRMeasurement, response::AbstractArray)
    N = prepadding(sig)
    N > 0 || throw(ArgumentError("Can't estimate noise without prepadding"))
    _noisefloor(response, N)
end

_noisefloor(response::AbstractVector, N) = sum(response[1:N, :].^2, 1) / N
_noisefloor(response::AbstractMatrix, N) = sum(response[1:N].^2) / N

include("options.jl")
include("util.jl")
include("golay.jl")
include("impulse.jl")
include("expsweep.jl")
include("mls.jl")
include("rpms.jl")

# workaround needed because DSP.jl doesn't handle SampleBufs
# and Float32s well. We dispatch to an internal _analyze to avoid method
# ambiguity issues if each measurement type defined its own `analyze` method.
analyze(sig::IRMeasurement, response::AbstractArray; kwargs...) =
        _analyze(sig, response; kwargs...)
analyze(sig::IRMeasurement, response::SampleBuf; kwargs...) =
        SampleBuf(_analyze(sig, Float64.(response.data); kwargs...), samplerate(response))


end # module
