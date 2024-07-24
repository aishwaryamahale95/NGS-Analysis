# NGS-analysis
Python program to analyze a vcf file

## Overview
This script is written to analyze VCF files to perform the following tasks:

1.  Count the total number of variants per chromosome.
2.  Count the number of heterozygous and homozygous variants per chromosome.
3.  Identify the genotype with the largest indel.
4.  Outputs variants supported by less than 15% of the reads.

## Requirements
1.  Python 3.x
2.  'argparse' module (parse command line input)
3.  'os' module (output to the directory)
4.  'vcf' module (can be installed via pip: `pip install PyVCF' or 'conda install -c bioconda pyvcf')
## Usage
### Command Line Arguments
'vcf_file_path': Path to the input VCF file.

### Running a test file
python script.py test.vcf

### Running the script
To run the script, use the following command in your terminal:

python script.py <path_to_vcf_file>

## Output
The script will generate an output folder named 'vcf_analysis_results' which will include the output file 'vcf_analysis_output.txt'. This output file will contain the results for the vcf file analysis tasks.
