# Selfing ML Sequences

This repository contains a Python script and an R script for estimating and visualizing selfing rates from multilocus sequence data.

## Installation

Clone the repository and make sure you have Python 3 installed, along with the necessary packages to run the script. You will also need R to visualize the results with the `plot_selfing_sequences.R` script.

## Usage

### Running the Python Script

The **`selfing_ML_sequences.py`** script takes a multifasta file and a species name as input and generates a tidy output file.

#### Command

```bash
python3 selfing_ML_sequences.py inputFile speciesName > outputFile
```

#### Example  
```bash
python3 selfing_ML_sequences.py Hibiscus_populations.fasta Hibiscus_laevis > selfing_Hibiscus_Hibiscus_laevis.txt
```
  
### Input File Format  
The input file is a multifasta file where each sequence ID follows this format:  
```shell
>{locus name}|{species name}|{individual name}|{Allele_1 or Allele_2}
ATGTCA....
```

The four fields in the sequence ID are separated by a pipe | and represent:
1. Locus name  
2. Species name  
3. Individual name  
4. Allele (Allele_1 or Allele_2)  


#### Example Data  
```objectivec
>1|Hibiscus_moscheutos|M10|Allele_1
CTTTCTATCACTGTTTCGTATGNGGGGCAGTAGAAGTATGCTGCCGAGTACCTTTCACTCTCTTGGTTGGCAATCACTCGGTGTGT
>1|Hibiscus_moscheutos|M10|Allele_2
CTTTCTATCACTGTTTCGTATGNGGGGCAGTAGAAGTATGCTGCCGAGTACCTTTCACTCTCTTGGTTGGCAATCACTCGGTGTGT
```

 
