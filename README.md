# Sac
Sac is a state-of-the-art lossless audio compression model

Lossless audio compression is a complex problem, because PCM data is highly non-stationary and uses high sample resolution (typically >=16bit). That's why classic context modelling suffers from context dilution problems. Sac employs a simple OLS-NLMS predictor per frame including bias correction. Prediction residuals are encoded using a sophisticated bitplane coder including SSE and various forms of probability estimations. Meta-parameters of the predictor are optimized via binary search (or [DDS](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005WR004723)) on by-frame basis. This results in a highly asymmetric codec design. We throw a lot of muscles at the problem and archive only little gains - by practically predicting noise. 

I've tried adding neural networks, wavelet decomposition and predictor mixing in various ways, but so far with not much success. 
Future improvements have to mix several predictions on a bit-level by remapping residuals to probability distributions. 

This program wouldn't exist without the help from the following people (in no particular order):

Matt Mahoney, Dmitry Shkarin, Eugene D. Shelwien, Florin Ghido

## Benchmark
|Program|Parameters|
|:-|:-|
|SAC v0.6.3|--high --optimize high --sparse-pcm|
|OFR v5.100|--preset max|
|paq8px_v181fix1|-6|

Numbers are bits per sample (smaller is better)

| Name  | SAC | OFR | paq8px |
|:---|---:|---:|---:|
|ATrain|7.058|7.156|7.376|
|BeautySlept|7.539|7.790|7.846|
|chanchan|9.696|9.778|9.740|
|death2|5.091|5.465|5.224|
|experiencia|10.872|10.915|10.985|
|female_speech|4.385|4.498|4.691|
|FloorEssence|9.162|9.409|9.506|
|ItCouldBeSweet|8.208|8.310|8.362|
|Layla|9.558|9.571|9.742|
|LifeShatters|10.771|10.808|10.890|
|macabre|9.007|9.026|9.266|
|male_speech|4.292|4.256|4.532|
|SinceAlways|10.348|10.409|10.479|
|thear1|11.388|11.400|11.496|
|TomsDiner|7.001|7.108|7.087|
|velvet|9.833|9.990|10.059|
|*Mean*|**8.388**|8.493|8.580|
