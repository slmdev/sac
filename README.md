# Sac
Sac is a state-of-the-art lossless audio compression model

Lossless audio compression is a complex problem, because PCM data is highly non-stationary and uses high sample resolution (typically >=16bit). That's why classic context modelling suffers from context dilution problems. Sac employs a simple OLS-NLMS predictor per frame including bias correction. Prediction residuals are encoded using a sophisticated bitplane coder including SSE and various forms of probability estimations. Meta-parameters of the predictor are optimized with [DDS](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005WR004723) on by-frame basis. This results in a highly asymmetric codec design. 

This program wouldn't exist without the help from the following people (in no particular order):

Matt Mahoney, Dmitry Shkarin, Eugene D. Shelwien, Florin Ghido, Grzegorz Ulacha

## Technical features
* Input: wav file with 1-16 bit sample size, mono/stereo, pcm
* Output: sac file including all input metadata
* Decoded wav file is bit for bit identical to input wav file
* MD5 of raw pcm values

## Technical limitations
Sac uses fp64 for many internal calculations. The change of compiler options or (cpu-)platform might effect the output. Use at your own risk and for testing purposes only.
 
## Benchmarks
**Sac v0.7.7**

16 files (51.014.742 bytes) parallel on i7-13700H.

**Asymmetric encoding profiles** - bits per sample (bps) is mean bps over all files
|Profile|Size|Enc-time|Dec-time|bps|
|:-|:-|:-|:-|:-|
|--normal|26.812.983|00:00:44|00:00:33|8.441|
|--high|26.710.610|00:04:01|00:00:55|8.408|
|--veryhigh|26.643.674|00:33:37|00:00:56|8.386|
|--best|26.571.729|05:45:41|00:01:03|8.364|

&nbsp;

**Comparison with other lossless audio codecs**
|Program|Parameters|Source|
|:-|:-|:-|
|Sac v0.7.7|--best|open|
|OFR v5.100|--preset max|closed|
|paq8px_v208fix1|-6|open|
|MAC v10.44|--insane|open|

Numbers are bits per sample (bps)
| Name  | Sac | OFR | paq8px | MAC |
|:---|---:|---:|---:|---:|
|ATrain|7.026|7.156|7.353|7,269|
|BeautySlept|7.420|7.790|7.826|8,464|
|chanchan|9.691|9.778|9.723|9,951|
|death2|5.075|5.465|5.215|6,213|
|experiencia|10.866|10.915|10.963|11,005|
|female_speech|4.375|4.498|4.708|5,190|
|FloorEssence|9.086|9.409|9.488|9,537|
|ItCouldBeSweet|8.199|8.310|8.330|8,531|
|Layla|9.530|9.571|9.725|9,783|
|LifeShatters|10.766|10.808|10.868|10,838|
|macabre|8.987|9.026|9.249|9,172|
|male_speech|4.278|4.256|4.509|5,255|
|SinceAlways|10.343|10.409|10.455|10,522|
|thear1|11.382|11.400|11.474|11,451|
|TomsDiner|6.994|7.108|7.057|7,432|
|velvet|9.803|9.990|10.030|10,461|
|*Mean*|**8.364**|8.493|8.561|8,817|

