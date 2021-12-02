# decomplex1
A fork of deMULTIplex by Chris Mcginnis  
Modified by John Bassett

# Description

Most cell and nuclei demultiplexing algorithms isolate signal by modeling the distribution of noise and removing it. The deMULTIplex algorithm provides the capability to optimize signal to noise by using thresholding in the signal isolation process. 

decomplex1 is a version of deMULTIplex which takes advantage of signal thresholding to demultiplex noisy data by altering the the thresholding parameters for the barcode signal distribution.

## Background

Single nucleus sequencing introduced the ability to sequence frozen tissue samples. In 2019, Gaublomme et Al. showed that nuclei hashing with barcoded antibodies works robustly for nuclei isolated from fresh frozen tissue. 

With the introduction of tumor slice culture models, primary patient tissue samples can be cultured for days before being frozen. We have found that these frozen tumor slices can provide a significant challenge to nuclei hashing by dramatically increasing the background signal.

To improve the cost-efficiency of sequencing tumor slice culture experiments, we make modifications to deMULTIplex that improve sample demultiplexing of noisy data sets 

## Differences in signal distribution

Hashing whole cells and fresh frozen tissue, produces little background and signal is easily isolated. The figure below shows a sample distribution of hashing data from whole cells compared with data that we generated from tumor slice cultures. 

![Sample Distributions](/Figures/SampleDistributions.jpg)



In our slice culture hashing experiments the majority of data points fall in the background peak. This background also contains a signficantly greater number of reads compared with the signal.

deMULTIplex generates thresholds by identifying the two largest modes in the dataset (blue and red lines in the image above). Noteably, nucelei from tumor slices have no second mode and therefore require an alternative method to thresholding.

## Alternate Thresholding

Instead of locating relative maxima, decomplex1 uses an estimate of the multiplet rate. This can be calculated based on the number of cells that have been loaded into the 10x controller. The actual quantity of doublets called will depend on the quantile which is optimized to call singlets. To adequatly perform a full range of quantile sweeps, it is beneficial to slightly underestimate the number of multiplets.


![Sample Thresholding](/Figures/EstMultiplets.jpg)

## Results

When alternative thresholding is used, the vast majority of barcoded nuclei can be correctly called. There remains a population of cells that are sufficiently obscured by the background signal that they cannot be called. Note that demultiplexing by barcoding calls doublets differently than by genotype. This is becuase barcoding allows for identification of homotypic doublets where as genotyping allows for identification of heterotypic doublets. It is always beneficial to use both methods, but when the samples are all the same genotype, such as in slice culture experiments, barcoding is the only feasibly method.

![Results](/Figures/Results.jpg)