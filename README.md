# Sac
Sac is a state-of-the-art lossless audio compression model

Lossless audio compression is a complex problem, because PCM data is highly non-stationary and uses high sample resolution (typically >=16bit). That's why classic context modelling suffers from context dilution problems. Sac employs a simple OLS-NLMS predictor per frame including bias correction. Prediction residuals are encoded using a sophisticated bitplane coder including SSE and various forms of probability estimations. Meta-parameters of the predictor are optimized with [DDS](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005WR004723) on by-frame basis. This results in a highly asymmetric codec design. 

This program wouldn't exist without the help from the following people (in no particular order):

Matt Mahoney, Dmitry Shkarin, Eugene D. Shelwien, Florin Ghido, Grzegorz Ulacha

## Technical Features
* Input: wav file with 1-16 bit sample size, mono/stereo, PCM
* Output: sac file including all input metadata
* Decoded wav file is bit for bit identical to input wav file
* MD5 of raw PCM values
 
## Benchmark
|Program|Parameters|Source|
|:-|:-|:-|
|SAC v0.6.9|--veryhigh|open|
|OFR v5.100|--preset max|closed|
|paq8px_v208fix1|-6|open|
|MAC v10.44|--insane|open|

Numbers are bits per sample (smaller is better)

| Name  | SAC | OFR | paq8px | MAC |
|:---|---:|---:|---:|---:|
|ATrain|7.038|7.156|7.353|7,269|
|BeautySlept|7.450|7.790|7.826|8,464|
|chanchan|9.695|9.778|9.723|9,951|
|death2|5.085|5.465|5.215|6,213|
|experiencia|10.871|10.915|10.963|11,005|
|female_speech|4.383|4.498|4.708|5,190|
|FloorEssence|9.110|9.409|9.488|9,537|
|ItCouldBeSweet|8.196|8.310|8.330|8,531|
|Layla|9.542|9.571|9.725|9,783|
|LifeShatters|10.771|10.808|10.868|10,838|
|macabre|9.001|9.026|9.249|9,172|
|male_speech|4.292|4.256|4.509|5,255|
|SinceAlways|10.344|10.409|10.455|10,522|
|thear1|11.383|11.400|11.474|11,451|
|TomsDiner|6.997|7.108|7.057|7,432|
|velvet|9.768|9.990|10.030|10,461|
|*Mean*|**8.370**|8.493|8.561|8,817|

