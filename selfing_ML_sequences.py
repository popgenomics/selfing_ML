import argparse
import numpy as np
from scipy.special import logsumexp
from Bio.SeqIO import parse
from collections import defaultdict

# get the filename to read
parser = argparse.ArgumentParser(description="Read a genotype file and extract genotype list and p_list.")

# input file name
parser.add_argument("filename", type=str, help="Path to the input file")

# focal species name
parser.add_argument("species_name", type=str, help="Name of the species")

# parse the arguments
args = parser.parse_args()

# get arguments
filename = args.filename
species_name = args.species_name

# Define the range of selfing rates to explore
s_space = np.linspace(0, 1, 21)
s_space = [ round(i, 2) for i in s_space ]

# Reads the fastafile
def parse_fasta_2(filename, species_name, max_locus=50000):
    """
    Reads a FASTA file and extracts sequences corresponding to the specified species.

    Args:
        filename (str): The path to the FASTA file to read.
        species_name (str): The name of the species to filter.
        max_locus (int): The maximum number of locus to record.

    Returns:
        dict: A dictionary where the key is the locus and individual, and the value is a list containing two sequences (Allele_1 and Allele_2).
    """

    # Dictionnary to store the sequences
    fasta_dict = defaultdict(lambda: {"Allele_1": None, "Allele_2": None})

    # Read the FASTA
    for record in parse(filename, "fasta"):
        # Extract information from the header
        header = record.id.split("|")
        if len(header) != 4:
            print(f"En-tête mal formatée : {record.id}")
            continue

        locus, species, individual, allele = header

        # Filter the focal species
        if species != species_name:
            continue

        # Check the name of the allele and add the corresponding sequence
        if allele == "Allele_1":
            fasta_dict[(locus, individual)]["Allele_1"] = str(record.seq)
        elif allele == "Allele_2":
            fasta_dict[(locus, individual)]["Allele_2"] = str(record.seq)

    # Convert the dictionary to only keep individuals with two alleles
    result_dict = {}
    locus_count = 0

    for (locus, individual), seqs in fasta_dict.items():
        if seqs["Allele_1"] is not None and seqs["Allele_2"] is not None:
            result_dict[(locus, individual)] = [seqs["Allele_1"], seqs["Allele_2"]]
            locus_count += 1

        # Stop if we reach the maximum number of locus
        if locus_count >= max_locus:
#            print(f"Maximum de {max_locus} locus atteint. Arrêt de la lecture.")
            break

    return result_dict

def parse_fasta(filename, species_name):
	"""
	Reads a FASTA file and extracts sequences corresponding to the specified species.

	Args:
		filename (str): The path to the FASTA file to read.
		species_name (str): The name of the species to filter.

	Returns:
		dict: A dictionary where the key is the locus and individual, and the value is a list containing two sequences (Allele_1 and Allele_2).
	"""

	# Dictionnary to store the sequences
	fasta_dict = defaultdict(lambda: {"Allele_1": None, "Allele_2": None})

	# Read the FASTA
	for record in parse(filename, "fasta"):
		# Extract information from the header
		header = record.id.split("|")
		if len(header) != 4:
			print(f"En-tête mal formatée : {record.id}")
			continue

		locus, species, individual, allele = header

		# Filter the focal species 
		if species != species_name:
			continue

		# Check the name of the allele and add the corresponding sequence
		if allele == "Allele_1":
			fasta_dict[(locus, individual)]["Allele_1"] = str(record.seq)
		elif allele == "Allele_2":
			fasta_dict[(locus, individual)]["Allele_2"] = str(record.seq)

	# Convert the dictionnary to only keep individuals with two alleles
	result_dict = {
		(locus, individual): [seqs["Allele_1"], seqs["Allele_2"]]
		for (locus, individual), seqs in fasta_dict.items()
		if seqs["Allele_1"] is not None and seqs["Allele_2"] is not None
	}

	return result_dict


# Reads the data file and returns 2 lists: genotypes and allele frequencies
def read_data_file(filename):
	"""
	Reads the file and extracts the genotype list and p_list.

	Args:
		filename (str): Path to the input file.

	Returns:
		genotype_list (list): List of genotypes (strings).
		p_list (list): List of allele frequencies (floats).
	"""
	with open(filename, 'r') as file:
		lines = file.readlines()

	# Ignore the first line (comment)
	genotype_line = lines[1].strip()
	p_line = lines[2].strip()

	# Create lists from the lines
	genotype_list = genotype_line.split(',')
	p_list = [float(x) for x in p_line.split(',')]

	return genotype_list, p_list


def get_genotypes(result, locus):
	"""
	Analyzes a given locus to identify polymorphic positions and returns genotypes and allele frequencies.

	Args:
		result (dict): Dictionary containing sequences, as returned by the parse_fasta function.
		locus (str): Name of the locus to analyze.

	Returns:
		tuple: A list of alternative allele frequencies (p_list) and a dictionary of genotypes (genotype_dic).
	"""

	# Filter the dictionnary to only keep sequences of the specified locus
	sequences = [
		(individual, seqs[0], seqs[1])
		for (locus_name, individual), seqs in result.items()
		if locus_name == locus
	]

	if not sequences:
		print(f"Locus '{locus}' not find in the data.")
		return [], {}

	# Check the length of the sequences because they have to be aligned
	seq_length = len(sequences[0][1])
	for _, allele_1, allele_2 in sequences:
		if len(allele_1) != seq_length or len(allele_2) != seq_length:
			raise ValueError("Sequences of individuals are not aligned or have different lengths.")

	# Get the positions that are polymorphic
	polymorphic_positions = []
	ref_alt_bases = []

	for i in range(seq_length):
		# Get the bases of the two alleles for all individuals at position 'i'
		bases = set()
		for _, allele_1, allele_2 in sequences:
			if allele_1[i] in "ATGC":
				bases.add(allele_1[i])
			if allele_2[i] in "ATGC":
				bases.add(allele_2[i])

		# Check if the position is bi-allelic
		if len(bases) == 2:
			polymorphic_positions.append(i)
			ref_base, alt_base = list(bases)
			ref_alt_bases.append((ref_base, alt_base))

	# Produce p_list et genotype_dic
	p_list = [] # Contains the frequency of 'A' allele for different polymorphic positions.
	genotype_dic = defaultdict(list) # Contains the individual diploid genotypes for different polymorphic positions.

	for pos, (ref_base, alt_base) in zip(polymorphic_positions, ref_alt_bases):
		# Count the occurence of alleles to compute the allele frequency
		alt_count = 0
		total_count = 0

		for individual, allele_1, allele_2 in sequences:
			base1 = allele_1[pos]
			base2 = allele_2[pos]

			# Ignore the missing data ('N')
			if base1 not in "ATGC" or base2 not in "ATGC":
				continue

			if base1 == alt_base:
				alt_count += 1
			if base2 == alt_base:
				alt_count += 1

			total_count += 2

		# Compute the frequency of the alternative allele
		if total_count == 0:
			freq_alt = 0.0  # No available allele
		else:
			freq_alt = alt_count / total_count
		p_list.append(freq_alt)

		# Determine the genotype for each individual
		for individual, allele_1, allele_2 in sequences:
			base1 = allele_1[pos]
			base2 = allele_2[pos]

			if base1 not in "ATGC" or base2 not in "ATGC":
				genotype = "NA"  # Genotype is missing
			elif base1 == ref_base and base2 == ref_base:
				genotype = "00"
			elif base1 == alt_base and base2 == alt_base:
				genotype = "11"
			elif (base1 == ref_base and base2 == alt_base) or (base1 == alt_base and base2 == ref_base):
				genotype = "10"
			else:
				genotype = "NA"  # Unexpected case, assumed to be 'missing'

			genotype_dic[individual].append(genotype)

	# Check that all individuals share the same number of called genotypes
	num_polymorphic_sites = len(polymorphic_positions)
	for individual in genotype_dic:
		while len(genotype_dic[individual]) < num_polymorphic_sites:
			genotype_dic[individual].append("NA")

	return p_list, dict(genotype_dic)


def get_polymorphic_loci(result, max_loci=10000):
	"""
	Returns a list of loci containing at least one polymorphic site, with an optional maximum limit on the number of loci.

	Args:
		result (dict): Dictionary containing sequences, as returned by the parse_fasta function.
		max_loci (int): Maximum number of polymorphic loci to return (default is 5000).

	Returns:
		list: List of locus names that have at least one polymorphic site.
	"""
	polymorphic_loci = []

    # The list of uniq loci
	loci_list = list({locus for (locus, _) in result.keys()})

    # Loop over loci to check the presence of polymorphic sites
	for locus_i in loci_list:
		p_list, genotype_dic = get_genotypes(result, locus_i)

		if len(p_list) > 0:
			polymorphic_loci.append(locus_i)

        # Stop if we found a 'maximum' number of polymorphic loci
		if len(polymorphic_loci) >= max_loci:
#			print(f"Nombre maximum de loci polymorphiques atteint : {max_loci}")
			break

	return polymorphic_loci


# Compute the probability to observe a given genotype 11/10/00 as a function of the frequence fA
# Genotype AA is coded as 11.
# Genotype aa is coded as 00.
def P_geno(genotype, fA, s, epsilon=1e-10):
	"""
	Calculate the log-probability of a given genotype based on allele frequency (fA)
	and selfing rate (s) with numerical stability improvements.
	"""
	if genotype=='NA':
		return(0)

	if not (0 <= fA <= 1):
		raise ValueError("Allele frequency fA must be between 0 and 1.")
	if not (0 <= s <= 1):
		raise ValueError("Selfing rate s must be between 0 and 1.")

	# Compute the allele frequencies and the inbreeding coefficient F
	p = fA
	q = 1 - p
	F = s / (2 - s)
	
	# Calculate genotype probabilities with log1p for stability
	with np.errstate(divide='ignore'):
		P_11 = np.log(p**2 + p * q * F + epsilon) # p² + pqF
		P_10 = np.log(2 * p * q * (1 - F) + epsilon) # 2pq.(1-F)
		P_00 = np.log(q**2 + p * q * F + epsilon) # q² + pqF

	# Select the appropriate log-probability based on the genotype
	if genotype == '11':
		return P_11
	elif genotype == '00':
		return P_00
	elif genotype in ['10', '01']:  # Handle both heterozygote notations
		return P_10
	else:
		raise ValueError("Genotype must be '11', '10', or '00'.")


# Calculate the log-likelihood as the sum of P_geno over all polymorphic sites
# for a given individual and a given selfing rate s
def calculate_log_likelihood(result, s, loci_list):
	# get the list of individuals
	individuals_list = list({individual for (_, individual) in result.keys()})

	# store the LogLik for all individuals
	logLik_dict = defaultdict(float)

	# loop over all sites, of all loci, of all individuals
	for locus_i in loci_list:
		p_list, genotype_dic = get_genotypes(result, locus_i)

		# if some polymorphic sites were found in the locus_i
		nPolSite = len(p_list)
		if nPolSite>0:
			for site_i in range(nPolSite):
				for ind_i in genotype_dic.keys():
					logLik_dict[ind_i] += P_geno(genotype=genotype_dic[ind_i][site_i], fA=p_list[site_i], s=s, epsilon=1e-10)

	return(logLik_dict)

# Read the fasta dataset
result = parse_fasta_2(filename, species_name)
polymorphic_loci = get_polymorphic_loci(result)
#print('Number of polymorphic_loci: {nP}'.format(nP=len(polymorphic_loci)))

logLik_dict_final = defaultdict(list)

for s_i in s_space:
#	print('Compute the likelihood of the data for s={s_i}'.format(s_i=s_i))
	# get the logLik for all individuals, for a given value of 's'
	res_s_i = calculate_log_likelihood(result=result, s=s_i, loci_list=polymorphic_loci)
	for ind_i in res_s_i:
		logLik_dict_final[ind_i].append(res_s_i[ind_i])

print('species\tindividual\tselfing_rate\tlogLik')
for ind_i in res_s_i:
	for i in range(len(s_space)):
		print('{species}\t{ind}\t{s}\t{LL}'.format(species=species_name, ind=ind_i, s=s_space[i], LL=logLik_dict_final[ind_i][i]))

