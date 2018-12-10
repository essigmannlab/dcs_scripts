# Mutational Analysis

Draws mutational spectra for given .mutpos files (outputs of DCS pipeline); generates the following files:
* 2 spectra for each sample's total mutations: normalized ("-prop") and unnormalized ("-freq")
* 2 spectra for each sample's unique mutations: normalized ("-prop") and unnormalized ("-freq")
* 2 spectra for combined total mutations of all samples
* 2 spectra for combined unique mutations of all samples
* CSV files for each spectrum generated
* .mutpos files for each sample's unique mutations
* .mutpos files for each combination of mutations
* Intermediate files (`jf_counter.jf` and `ref_counts.txt`) for Jellyfish kmer counts
