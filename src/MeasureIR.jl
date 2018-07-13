module MeasureIR

using Compat: @warn
import Unitful
using Unitful: s, Hz
using SampledSignals: SampleBuf

export stimulus, analyze
export golay

abstract type IRMeasurement end

include("golay.jl")

end # module
