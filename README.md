<p align="center">
  <img width="450" src="https://ssfrr.github.io/MeasureIR.jl/readme_ir.png" alt="Simulated Impulse Response Plot">
</p>

# MeasureIR

[![Build Status](https://travis-ci.org/ssfrr/MeasureIR.jl.svg?branch=master)](https://travis-ci.org/ssfrr/MeasureIR.jl) [![Coverage Status](https://coveralls.io/repos/ssfrr/MeasureIR.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/ssfrr/MeasureIR.jl?branch=master) [![codecov.io](http://codecov.io/github/ssfrr/MeasureIR.jl/coverage.svg?branch=master)](http://codecov.io/github/ssfrr/MeasureIR.jl?branch=master)

MeasureIR is a Julia library for measuring and analyzing [impulse responses](https://en.wikipedia.org/wiki/Impulse_response) (IRs). Impulse responses are generally captured by playing some kind of test signal through the system under test, and then analyzing the response to extract the IR.

The simplest test signal is an actual impulse, like a gunshot, balloon pop, or hand clap. In this case there is no analysis necessary, because the recorded response is a direct impulse response. The downside to this simple approach is that the total energy in the impulse often can't be very large without overdriving the measurement equipment and causing nonlinearities.

For this reason it is more common to use other signals like pseudo-random noise or sine sweeps to characterize a system, as they can spread the energy over a longer time frame, which increases the signal-to-noise ratio of the measurement.

Measuring an impulse response takes place in the following steps:

1. Create a measurement (currently only `golay` is supported)
2. Generate a test signal for the measurement using `stimulus(m)`, where `m` is a measurement. This signal is a single-channel time-domain signal that could be played directly through a speaker, or saved to a file for later measurement.
3. Convolve the test signal with your system. This could be by playing the signal through a speaker into a room you're measuring, or using Julia's built-in `conv` function for testing. The result of this step should be a (possibly multichannel) response signal.
4. Analyze the system response to generate the impulse response. The form of this is `analyze(m, response)` where `m` is your measurement object and `response` is the measured output of your system in response to the stimuli.


## Example

```julia
using Plots: plot
using MeasureIR: golay, stimulus, analyze

meas = golay(4096)

# generate the test stimuli suitable for playback
stim = stimulus(meas)

# create a synthetic impulse response and simulate the system. This is where
# you'd normally play the stimuli through your system and record the response
irsim = 1./exp.(0:0.1:9.9) .* (rand(100) .- 0.5)
output = conv(stim, irsim)

# analyze to reconstruct the impulse response
ir = analyze(meas, output)

plot([irsim[1:100], ir[1:100]], labels=["Convolved IR", "Measured IR"])
```
