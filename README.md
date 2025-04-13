![Version](https://img.shields.io/github/release/slmdev/sac)
![Repo size](https://img.shields.io/github/repo-size/slmdev/sac)
![License](https://img.shields.io/github/license/slmdev/sac)

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
**Sac v0.7.18**

16 files (44.1Khz, stereo) - 51.014.742 bytes - parallel on i7-13700H

**Asymmetric encoding profiles** - bits per sample (bps) is mean bps over all files
|Profile||Size|Enc-time|Dec-time|bps|
|:-|-:|:-|:-|:-|:-|
|FLAC -8|100.0%|29.793.060|00:00:00|00:00:00|9.385|
|--normal|89.6%|26.691.217|00:00:22|00:00:20|8.401|
|--high|89.1%|26.546.317|00:03:42|00:00:50|8.352|
|--veryhigh|88.9%|26.480.793|00:17:38|00:00:53|8.330|
|--extrahigh|88.8%|26.453.344|00:48:46|00:00:52|8.322|
|--best|88.6%|26.397.759|04:37:00|00:01:00|8.303|

&nbsp;

**Comparison with other lossless audio codecs**
|Program|Parameters|Source|
|:-|:-|:-|
|Sac v0.7.18|--best|open|
|OFR v5.100|--preset max|closed|
|paq8px_v208fix1|-6|open|
|MP4ALS RM23|-b -p -z3|open|
|MAC v10.44|-c5000|open|
|FLAC v1.4.3|-8|open|

Numbers are bits per sample (bps)
| Name  | Sac | OFR | paq8px | MP4ALS | MAC | FLAC |
|:---|---:|---:|---:|---:|---:|---:|
|ATrain|6.991|7.156|7.353|7.232|7.269|7.962|
|BeautySlept|7.089|7.790|7.826|8.305|8.464|10.134|
|chanchan|9.658|9.778|9.723|9.886|9.951|10.206|
|death2|5.036|5.465|5.215|6.660|6.213|6.092|
|experiencia|10.825|10.915|10.963|10.992|11.005|11.428|
|female_speech|4.356|4.498|4.708|4.711|5.190|5.364|
|FloorEssence|9.034|9.409|9.488|9.509|9.537|10.201|
|ItCouldBeSweet|8.182|8.310|8.330|8.396|8.531|9.064|
|Layla|9.454|9.571|9.725|9.691|9.783|10.415|
|LifeShatters|10.748|10.808|10.868|10.836|10.838|11.194|
|macabre|8.968|9.026|9.249|9.076|9.172|10.043|
|male_speech|4.229|4.256|4.509|4.813|5.255|5.674|
|SinceAlways|10.329|10.409|10.455|10.473|10.522|11.254|
|thear1|11.362|11.400|11.474|11.425|11.451|11.783|
|TomsDiner|6.991|7.108|7.057|7.268|7.432|8.572|
|velvet|9.595|9.990|10.030|10.212|10.461|10.770|
|*Mean*|**8.303**|8.493|8.561|8.718|8.817|9.385|
|||||||
|Enc-time|04:37:00|00:00:09|00:05:37|00:00:15|00:00:01|00:00:00|
|Dec-time|00:01:00|00:00:02|00:05:40|00:00:13|00:00:01|00:00:00|

