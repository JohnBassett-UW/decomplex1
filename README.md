# decomplex1
A fork of deMULTIplex by Chris Mcginnis
Modified by John Bassett

# Description

Most cell and nuclei demultiplexing algorithms isolate signal by modeling the distribution of noise and removing it. The deMULTIplex algorithm provides the capability to optimize signal to noise by using thresholding in the signal isolation process. 

decomplex1 is a version of deMULTIplex which takes advantage of signal thresholding to demultiplex noisy data by altering the the thresholding parameters for the barcode signal distribution.

## Background

Single nucleus sequencing introduced the ability to sequence frozen tissue samples. In 2019, Gaublomme et Al. showed that nuclei hashing with barcoded antibodies works robustly for nuclei isolated from fresh frozen tissue. 

With the introduction of tumor slice culture models, primary patient tissue samples can be cultured for days befor being frozen. We have found that these frozen tumor slices can provide a significant challenge to nuclei hashing by dramatically increasing the background signal.

To improve the cost-efficiency of sequencing tumor slice culture experiments, we make modifications to deMULTIplex that improve sample demultiplexing of noisy data sets 

## Differences in signal distribution

Hashing whole cells and fresh frozen tissue, produces little background and signal is easily isolated. The figure below shows a sample distribution of hashing data from whole cells compared with data that we generated from tumor slice cultures. 

![Sample Distributions](/figures/SampleDistributions.pdf)

Whole cell hashing produces very little background and the vast majority of data points (blue line) contain zero reads from contaminating hashes.  

By contrast, in slice culture hashing experiments we see that the majority of data points fall in the background peak. This background also contains a signficantly greater number of reads compared with the signal.

## Changes to deMULTIplex

deMULTIplex generates thresholds by identifying the two highest modes in the dataset (blue and red lines in the image above). Noteably, nuceli from tumor slices have no second mode and therefore require an alternative method to thresholding. 
 
