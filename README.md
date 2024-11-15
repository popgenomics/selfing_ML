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
  
### Input File Format  
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

#### Output File Format  
The output file is in tidy format, making it easy to process with R's tidyverse package. It contains four columns:  
1. Species name
2. Individual name
3. Explored selfing rate values (default: a grid of 21 points between 0 and 1)
4. Log likelihood of the selfing rate

| Species          | Individual | Selfing Rate | Log-Likelihood      |
|------------------|------------|--------------|---------------------|
| Hibiscus_laevis  | L41        | 0.0          | -588.6670054354053  |
| Hibiscus_laevis  | L41        | 0.05         | -581.9224496416417  |
| Hibiscus_laevis  | L41        | 0.1          | -575.2568833051852  |
| Hibiscus_laevis  | L41        | 0.15         | -568.6800090923874  |
| Hibiscus_laevis  | L41        | 0.2          | -562.2066187379005  |
| Hibiscus_laevis  | L41        | 0.25         | -555.857834600882   |
| Hibiscus_laevis  | L41        | 0.3          | -549.6628670156902  |
| Hibiscus_laevis  | L41        | 0.35         | -543.661510900775   |
  
### Plotting the Results  
To graphically represent the results, you can use the plot_selfing_sequences.R script.  
  
#### Grant Execution Permissions  
```bash
chmod +x plot_selfing_sequences.R
```

#### Running the R Script
```bash
Rscript plot_selfing_sequences.R selfing_Hibiscus_Hibiscus_laevis.txt
```

This script will generate a PDF file with the plot:  

![alt text](https://github.com/popgenomics/selfing_ML/selfing_Hibiscus_Hibiscus_laevis.png)
 
