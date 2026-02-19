# gwsignal4gwsurr
gwsignal implementation to do PE runs using waveform models from gwsurrogate

# Installation
1. Clone this repo 
```bash
git clone git@github.com:varenya27/gwsignal4gwsurr.git
```
2. Install gwsignal4gwsurr 
```bash
python -m pip install .            # option 1
python -m pip install --editable . # option 2
```

# Using this interface to do PE with `bilby`
Add to the .ini file:
```bash
################################################################################
## Waveform arguments
################################################################################

waveform-generator=gwsignal4gwsurr.waveform_generator.SurrogateWaveformGenerator
waveform-approximant=NRHybSur3dq8
waveform-arguments-dict ={
  reference-frequency:20.,
}
```

# File tree
```
.
├── examples
│   └── GW150914.ini          ==> example bilby .ini file to do PE on GW150914 with NRSur7dq4
├── gwsignal4gwsurr
│   ├── gwsurr_playground.py  ==> playground for new/experimental models
│   ├── gwsurr.py             ==> finalized models get their own class here
│   ├── __init__.py
│   └── utils_bilby.py        ==> bilby facing code goes here


```

# Implemented models
1. NRHybSur3dq8
2. NRSur7dq4
