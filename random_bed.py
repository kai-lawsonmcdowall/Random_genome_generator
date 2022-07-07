"""
#Essentially a pythonic rewrite of bedtools.random to generate a random bed file as a control. 

There are three main functions to create our random bed file: 

chrom_size_dictionary: 
----------------------
this takes in a chrom sizes file in the format of: 

chr1    0   121258967
chr2    0   232323643

The chr at the beginning is optional, however, it does require the 0 column for the moment. 
This generates a list of dictionaries based on the chrom sizes file, each dictionary contains the name
of the chromosome, it's length, and the probability of it being selected at random

* chrom_picker(): 
-----------------
this essentially creates a list of length n based on your dictionary, it generates this list of chromosomes,
importantly, it selects the chromosomes the chromosomes with appropriate probabilities (e.g. chr1 is the largest
chromosome, so we expect to see it more than chr16, for example)

* break_generator(): 
--------------------
the major function, this takes the dictionary outputted from chrom_size_dictionary() and chrom_picker() and uses them to generate
the bed file itself.
"""

#%%
def chrom_size_dictionary(filepath):
    from natsort import natsorted
    #additional functionality - add a conditional allowing start to not exist. 
    #shouldn't need to add a sort functionality

    #expect the file to be in the format of "chr1 0   size" and import to a pandas dataframe
    import pandas as pd
    chromosome_df = pd.read_csv(filepath, sep = '\t', header=None)
    chromosome_df.columns = ['name','start','length']

    # add chr to the first column (assuming it doesn't exist)
    if 'chr1' not in chromosome_df['name'].values :
        chromosome_df["name"] = 'chr' + chromosome_df["name"].astype(str)

    #sort your dataframe

    #total number of breaks which we will use to calculate probability
    total_bp = chromosome_df["length"].sum()

    # remove the "start" column which should only contain 0's anyway and insert the probability values
    nostart_chr_df = chromosome_df.drop("start",1)
    nostart_chr_df["prob"] = chromosome_df["length"]/total_bp

    #convert to dictionary
    chrom_length_and_prob = nostart_chr_df.to_dict('records')
    print(chrom_length_and_prob)

    chrom_size_dictionary.dict = chrom_length_and_prob

#%%
chrom_size_dictionary("test_chrom_sizes/hg19.chrom.sizes.bed")

#%%
#generates a list of chromosomes based on the probabilites stored in the dictionary. 
def chromosome_picker(n, chromosome_dictionary):
    
    import random
    import matplotlib.pyplot as plt
    import numpy as np
    
    chrs = []
    
    #create lists of probabilites and chromosomes from our list of dictionaries. 
    chromosomes =  [chrom['name'] for chrom in chromosome_dictionary]
    probabilites = [chrom['prob'] for chrom in chromosome_dictionary] 

    #selecting a chromsome at random based on our probabilites
    for i in range(n): 
        chrs.append(random.choices(chromosomes, weights=probabilites)[0])


    #correctly sorting the chromosomes    
    chr_indicies = [chromosomes.index(chr) + 1 for chr in chrs]
    #graphing to check that the chromosomes are being picked at random by graphing, the logic being the longer chromosomes should be picked more often.
    y = np.array(sorted(chr_indicies))
    plt.hist(y, bins = len(chromosomes)); #need to adjust the bins to show the chromosomes themselves.
    plt.show()

    chromosome_picker.output_chromosomes = chrs

# %%
chromosome_picker(91081245, chrom_size_dictionary.dict)

# %%
#requires the list of random chrs that we've generated and our original chromosome dictionary. 
def break_generator(chromosome_list, chromosome_dictionary):
    
    from natsort import natsorted
    import random
    import csv

    start_list = []
    end_list = [] 
    random_ID_list = []
    strand = []

    
    sorted_chromosome_list = natsorted(chromosome_list)

    #sorting the chromosome list

    for i in sorted_chromosome_list: 
        chrom_length = [chrom["length"] for chrom in chromosome_dictionary if chrom["name"] == i][0]

        start = random.randint(0, chrom_length)
        start_list.append(start)
        
        end = start + 1
        end_list.append(end)
        n=0
        
    for x in range(1, len(sorted_chromosome_list)+1):
        random_ID = "random_count_" + str(x)
        random_ID_list.append(random_ID)

        strand_binary = random.randint(0,1)
        if strand_binary == 0:
            strand.append("-")
        else:
            strand.append("+")

    print("writing out a file containing", len(sorted_chromosome_list), "breaks")

    
    data = zip(sorted_chromosome_list,start_list,end_list,random_ID_list,strand)
    with open('output.bed', 'w', newline='') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        for sorted_chromosome_list,start_list,end_list,random_ID_list,strand in data:
            tsv_output.writerow([sorted_chromosome_list,start_list,end_list,random_ID_list,strand])


# %%
break_generator(chromosome_picker.output_chromosomes,chrom_size_dictionary.dict, write_output="yes")

