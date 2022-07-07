# Random_genome_generator
a tool that generates a BED-format random genome file based on a genome size file you supply it
performs a similar function to bedtools random, except written entirely in python

requires a tab seperated input file in the format of

chr 0 end
e.g.
1	0	249250621


the output will look like

chr star  end ID  strand
1 342352  523623  random_sequence_n +

works with any genome so long as chr sizes are correct, allows you to vary the size of the random sequences, and contains a matplotlib function which shows you the number of random generated sequences for each chromosome, if these decrease gradually, with a spike at the x chromosome, then this indicates the probabilites are correct (i.e. it is generating these random frequencies across the genome in an expected manner).

Is a little slow and will likely optimize in the future (i.e. 100 million line file takes 15ish minutes to produce). 
