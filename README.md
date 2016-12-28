# RECOAnalyzer
#

## SETUP (using CMSSW_8_1_0_pre15):

```
#### architechture set up
# For bash shell 
SCRAM_ARCH=slc6_amd64_gcc530
# For csh/tcsh shell:
setenv SCRAM_ARCH slc6_amd64_gcc530

#### CMSSW release set up
cmsrel CMSSW_8_1_0_pre15
cd CMSSW_8_1_0_pre15/src
cmsenv

## Tutorial software setup 
mkdir ECFAProj
cd ECFAProj
git clone https://github.com/drankincms/RECOAnalyzer.git
scram b -j8
cd RECOAnalyzer

```

## RUN:
```
cd python
cmsRun runEff.py
```
## PLOT:
```
root -l -b -q EffPlot.C
```
