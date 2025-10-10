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
4. change the waveform arguments in the ini file (eg. [see here](examples/GW150914.ini)):
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
5. generate submission and data files 
```bash
parallel_bilby_generation GW150914.ini
```
6. update the submission files in `outdir/submit/analysis_GW150914_0.sh`(for custom SLURM stuff); replace the stuff above the `mpirun` command; eg.
```bash
#!/bin/bash
#SBATCH --job-name=0_GW150914
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=24:00:00
#SBATCH --output=outdir/log_data_analysis/0_GW150914_%j.log
#SBATCH -p cpu,umd-cscdr-cpu,cpu-preempt

source ~/miniforge3/etc/profile.d/conda.sh
conda activate igwn-py310
```
# implemented models
1. NRHybSur3dq8
