# ProbabilisticInferenceForPromoterActivity

## General description
Transcriptional regulation can be tight and efficient. For example, with only ~10 repressors, many E. coli cells only contain few copies of proteins under the control of the lac operon. Fluorescent protein signals to monitor such a low level of promoter activity will be very close to that of the auto fluorescence. We developed a probabilistic algorithm to robustly infer the promoter activity in the low expression regime.

## Publication
This probabilistic inference algorithm was initially developed for 

Chang Chang, Mayra Garcia-Alcala, Leonor Saiz, Jose M.G. Vilar,  Philippe Cluzel. Robustness of DNA Looping Across Multiple Divisions in Individual Bacteria. bioRxiv, 2021, https://doi.org/10.1101/2021.12.28.474367

As suggested by the Referee, we decided to publish this algorithm and the codes as a separate paper.

## Sample datasets and demo
Sample datasets are available in the ExampleData folder. "MultipleLoops_data.tar" contains the time series for the Loops strain, and "NoLoop_data.tar" is for the No-loop strain.

As an example to run this algorithm, 
* Download all codes files and add them to the search path of MATLAB.
* Unzip one of the sample datasets. Take "MultipleLoops_data.tar" as an example. In the terminal of MATLAB, enter

```matlab
Inference('CL-54gamma', 'probabilistic');
```

One may change the setting for the algorithm in "Inference/Main/MethodParameter.m". For example, for "NoLoop_data.tar", we can change the setting to

```matlab
sigma.signal = 2.5;
sigma.derivative = [];
sigma.method = 'signal-frame-of-reference';
```

as a loose threshold.
