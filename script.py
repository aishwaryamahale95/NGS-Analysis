import argparse
import os
import vcf


# Read the vcf file
def vcf_file_read(lines):
    try:
        vcf_read = vcf.Reader(open(lines, 'r'))
        lines = [line for line in vcf_read]
        return lines
    except Exception as e:
        print(f"Error!! Please check the VCF file: {e}")
        return []

# 1.1 Count the total number of variants in the VCF file
def variants_count(lines, output_file):

    #total_variants = len(lines) #count the length of the variants
    #print("Total number of variants:", total_variants)

    # Count the total number of variants per chromosome
    variant_count = {}

    #Iterate through each line in vcf file and get the chromosome identifier
    for line in lines:  
        chrom = line.CHROM
        
        if chrom not in variant_count:  #check if the chromosome is not the dictionary variant_count
            variant_count[chrom] = 0
        variant_count[chrom] += 1
    output_file.write(f"Total number of variants per chromosome:\n")
    #print("Total number of variants per chromosome:")

    for chrom, count in variant_count.items():
        #print(f"Chromosome {chrom}: {count}")
        output_file.write(f"Chromosome {chrom}: {count}\n")


# 1.2 Count the total number of heterozygous and homozygous variants per chromosome
def variants_per_chromosome(lines, output_file):

    heterozygous_counts = {}
    homozygous_counts = {}

    for line in lines:
        chrom = line.CHROM
        genotype = line.samples[0].data.GT  #update the genotype for the line
        #print(genotype)
        

        if chrom not in heterozygous_counts:
            heterozygous_counts[chrom] = 0
        if chrom not in homozygous_counts:
            homozygous_counts[chrom] = 0

        if genotype == '0/1' or genotype == '1/0':      # Heterozygous variants (0/1 or 1/0)
            heterozygous_counts[chrom] += 1
        elif genotype == '1/1' or genotype == '0/0':    # Homozygous variants (1/1 or 0/0)
            homozygous_counts[chrom] += 1

    output_file.write(f"Heterozygous counts per chromosome:\n")
    #print("Heterozygous counts per chromosome:")
    
    for chrom, count in heterozygous_counts.items():    #Iterate through heterozygous_counts dictionary
        output_file.write(f"Chromosome {chrom}: {count}\n")
        #print(f"Chromosome {chrom}: {count}")

    output_file.write(f"Homozygous counts per chromosome:\n")
    #print("Homozygous counts per chromosome:")
    for chrom, count in homozygous_counts.items():      #Iterate through homozygous counts dictionary
        output_file.write(f"Chromosome {chrom}: {count}\n")
        #print(f"Chromosome {chrom}: {count}")

#2 Genotype with the largest indel in the VCF file
def largest_indel(lines, output_file):

    largest_indel_size = 0
    largest_indel = None

    #Iterate through each line and add the alleles from the REF and ALT columns respectively 
    for line in lines:
        ref = line.REF
        alt = line.ALT 
        if len(ref) != len(alt):  # This means that indels are present
            indel_size = abs(len(ref) - len(str(alt)))
            if indel_size > largest_indel_size:     #check if the indel size is larger than the largest indel size
                largest_indel_size = indel_size
                largest_indel = line.samples[0].data.GT     #update the genotype for the largest indel 
                
    if largest_indel:
        output_file.write(f"Genotype of the largest indel: {largest_indel}\n")
        #print(f"Genotype of the largest indel: {largest_indel}")
    else:
        output_file.write(f"No indels found in the VCF file.\n")
        #print("No indels found in the VCF file.")

#3 Output variants that are supported only by <15% of the reads
#Define a function to calculate percentage of reads supporting ALT allele

def supported_variants(lines, output_file):
    variant = []
    for line in lines:
        for sample in line.samples:
            if sample['DP'] > 0:  #check if there is coverage (read depth > 0)
                #print(sample['DP'])
                alt_alleles = sample['AD'][1] #since the AD[1] refers to the alternate alleles
                total_reads = sample['DP'] 
                if alt_alleles / total_reads < 0.15:
                    variant.append(line)

    #check for variants with less than 15% reads
    if variant:
        output_file.write(f"Variants supported only by less than 15% of the reads:\n")
        #print("Variants supported only by less than 15% of the reads:")
        for var in variant:
            output_file.write(f"{var}\n")
            #print(f"{var}\n")
    else:
        output_file.write(f"No variants supported by less than 15% of the reads found.\n")
        #print("No variants supported by less than 15% of the reads found.")


def main():
    parser = argparse.ArgumentParser(description='Process a VCF file.') #parsing the input vcf file
    parser.add_argument('vcf_file_path', type=str, help='Provide path to the VCF file')
    args = parser.parse_args()

    program_dir = os.path.dirname(os.path.abspath(__file__))     
    output_folder = os.path.join(program_dir, 'vcf_analysis_results')       #the output folder will be saved in the same folder as script.py file
    os.makedirs(output_folder, exist_ok = True)
    output_file = os.path.join(output_folder, 'vcf_analysis_output.txt')

    #Read VCF file and extract each line
    with open(output_file, 'w') as output_file:
        lines = vcf_file_read(args.vcf_file_path)

    #Calling the defined functions
        variants_count(lines, output_file)
        variants_per_chromosome(lines, output_file)
        largest_indel(lines, output_file)
        supported_variants(lines, output_file)
            
if __name__ == "__main__":
    main()
