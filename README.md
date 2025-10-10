# gwsignal3gwsurr
gwsignal implementation to use waveform models from gwsurrogate

# installation
1. clone this repo 
```bash
git clone git@github.com:varenya27/gwsignal4gwsurr.git
```
2. install gwsignal4gwsurr 
```bash
python -m pip install .            # option 1
python -m pip install --editable . # option 2
```
# usage for (parallel) bilby pe
1. check that a class has been initialized for the necessary waveform model (eg. does `gwsurr.py` contain `NRHybSur3dq8_gwsurr` as a class?)
2. check that a wrapper has been implemented for the necessary waveform model (eg. does `NRHybSur3dq8_wrapper.py` exist?)
3. check that an appropriate `parameter_conversion` function has been implemented in `utils_bilby.py`
4. change the waveform arguments in the ini file:
```bash
################################################################################
## Waveform arguments
################################################################################

waveform-generator=gwsignal4gwsurr.utils_bilby.get_waveform_generator
waveform-approximant=NRSurr
frequency-domain-source-model=gwsignal4gwsurr.NRHybSur3dq8_wrapper.NRHybSur3dq8_wrapper
waveform-arguments-dict ={
  reference-frequency:20., 
  f-min:20.,
}
conversion-function=gwsignal4gwsurr.utils_bilby.parameter_conversion

```
5. run `parallel_bilby_generation` etc

# implemented models
1. NRHybSur3dq8
