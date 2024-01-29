# Sac
Sac is a state-of-the-art lossless audio compression model

Lossless audio compression is a complex problem, because PCM data is highly non-stationary and uses high sample resolution (typically >=16bit). That's why classic context modelling suffers from context dilution problems. Sac employs a simple OLS-NLMS predictor per frame including bias correction. Prediction residuals are encoded using a sophisticated bitplane coder including SSE and various forms of probability estimations. Meta-parameters of the predictor are optimized via binary search (or [DDS](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005WR004723)) on by-frame basis. This results in a highly asymmetric codec design. We throw a lot of muscles at the problem and archive only little gains - by practically predicting noise. 

I've tried adding neural networks, wavelet decomposition and predictor mixing in various ways, but so far with not much success. 
Future improvements have to mix several predictions on a bit-level by remapping residuals to probability distributions. 

This program wouldn't exist without the help from the following people (in no particular order):

Matt Mahoney, Dmitry Shkarin, Eugene D. Shelwien, Florin Ghido

## Benchmark
|Program|Parameters|Source|
|:-|:-|:-|
|SAC v0.6.7|--veryhigh|open|
|OFR v5.100|--preset max|closed|
|paq8px_v181fix1|-6|open|
|MAC v10.44|--insane|open|

Numbers are bits per sample (smaller is better)

| Name  | SAC | OFR | paq8px | MAC |
|:---|---:|---:|---:|---:|
|ATrain|7.044|7.156|7.376|7,269|
|BeautySlept|7.463|7.790|7.846|8,464|
|chanchan|9.698|9.778|9.740|9,951|
|death2|5.089|5.465|5.224|6,213|
|experiencia|10.872|10.915|10.985|11,005|
|female_speech|4.387|4.498|4.691|5,190|
|FloorEssence|9.117|9.409|9.506|9,537|
|ItCouldBeSweet|8.201|8.310|8.362|8,531|
|Layla|9.551|9.571|9.742|9,783|
|LifeShatters|10.771|10.808|10.890|10,838|
|macabre|9.001|9.026|9.266|9,172|
|male_speech|4.289|4.256|4.532|5,255|
|SinceAlways|10.345|10.409|10.479|10,522|
|thear1|11.386|11.400|11.496|11,451|
|TomsDiner|6.999|7.108|7.087|7,432|
|velvet|9.810|9.990|10.059|10,461|
|*Mean*|**8.377**|8.493|8.580|8,817|

