const default_prepad=5*48000
# gain levels set from running scaling.jl
const impulse_gain = db2amp(-1)
const expsweep_gain = db2amp(-1)
const mls_gain = 0.38911592741900447
const golay_gain = 0.4143208008205667
const rpms_gain = 0.2316778907504581

const optiondoc_gain =
"""
- `gain::AbstractFloat`: Gain applied to the signal to reduce the chance of
  distortion. The default provides enough headroom so that the 99.99% of the
  analog stimulus is below -1dbFS. This means that if you set the gain of your
  playback system such that a full-scale sinusoid does not distort, this signal
  should not distort either.
"""

const optiondoc_prepad =
"""
- `prepad::Int`: Defaults to `$default_prepad`. Puts a period of silence at the
  beginning of the measurement that can be used to estimate the noise floor.
"""
