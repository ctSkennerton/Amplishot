# Amplishot: Amplicon-Shotgun
Currently microbial community profiling sudies rutienly use 454
pyrosequencing to generate 10,000 - 100,000 reads from particular 
variable regions of the 16S rRNA gene.  Unfortunately 454 pyrosequencing 
has a number of short falls such as homopolymer errors.  Furthermore
taxonomic resolution can be lost when using 454 pyrosequencing due
to the smaller fragment of the 16S rRNA gene that is being analyzed.
*Amplishot* combines amplification of the full 16S rRNA gene sequence
with *de novo* reconstruction of full-length 16S rRNA genes from specially
constructed "Amplishot" Illumina sequencing libraries or from
metagenomes. 
## Dependancies
* [Qiime](http://qiime.org) - tested only on version 1.6.0 & 1.8.0
*	[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) - tested with version 2.0.5
* [Pear](http://www.exelixis-lab.org/pear)
* [samtools](http://samtools.github.io/)
* [pyYAML](http://pyyaml.org/)
* numpy
* pysam
* [biom](http://biom-format.org/) - version 1 required (most likely version 2 will break things;
  only tested with version 1)


## OTU Clustering Dependancies
*	[cd-hit](http://cd-hit.org) - tested with 4.5.4
* uclust - now packaged with Qiime if using version 1.8.0 or higher.
  NOTE: uclust is not the same as usearch although google searches for
  the former will return the latter

### Assembly Dependancies
You must have one of the following
*	[phrap](http://www.phrap.org/) - tested with version 1.09518 (not
  currently recomented as it does not scale)
* [Ray](http://denovoassembler.sourceforge.net/) - tested with version 2.3.1
* [velvet](https://www.ebi.ac.uk/~zerbino/velvet/)

## Installation
You can either download the [latest source code](https://github.com/ctSkennerton/Amplishot)
 or a [particular version](https://github.com/ctSkennerton/Amplishot/tags) from github.
Once downloaded change into the Amplishot directory and run the following:
 
	 [sudo] python setup.py install
 
or if you do not have sudo on your computer use the `--prefix` option to change
the installation directory.

## Command-line interface  
Amplishot has a single executable called `amplishot`; you can see some basic help
by running the command `amplishot -h`.
The command-line options for Amplishot only offer a subset of the options that are 
available.  Most options are changed by using a configuration file.  The Amplishot 
configuration file is written in [YAML](http://yaml.org), which is a simple markup 
language; before you try to modify the configuration file it might be helpful to 
read up on the YAML [syntax](http://yaml.org/spec/1.1/).
Command-line options and a configuration file can be used in tandum.  Any options
specified on the command-line will **overwrite** the corresponding value in the
configuration file.  If changes have been made to a configuration by using 
command-line options, a new configuration file will be outputted to the global
output directory with a datetime signature so that no previous configuration details
are lost.  A new configuration file will not be outputted if there are no changes
to the current configuration set.

### Configuration options and their values
*	`threads`: Sets the number of threads/processes to use in the Amplishot pipeline. 
   	The value should be a single integer number (default: 1)
*	`log_level`: Changes the verbosity of logging messages.  The options from most 
   	verbose to least are: DEBUG, INFO, WARN, ERROR, CRITICAL (default: INFO)
*	`output_directory`: This is the directory where all results will be outputted.  
   	By default it is the current directory, symbolized by a '.' character
*	`input_raw_reads`: This must be a list of files raw Illumina sequencing read files
   	to input into Amplishot.  The format of the input reads must be a *list-of-lists*,
   	which can be added to the configuration in two ways:
        
		input_raw_reads:
	        - [/full/path/to/sample1.1.fq, /full/path/to/sample1.2.fq]
		    - [/full/path/to/sample2.1.fq, /full/path/to/sample2.2.fq]
		   
		input_raw_reads:
			- 
			    - /full/path/to/sample1.1.fq
				- /full/path/to/sample1.2.fq
			-
				- /full/path/to/sample2.1.fq
				- /full/path/to/sample2.2.fq

*	`aliases`: Use this option to set the sample names to be used in the output files.
	By default the filename is used without the file extension.  The form of the values
	must be a YAML list specified by either:
	
		aliases:
			- alias1
			- alias2
		
		aliases: [alias1, alias2]

* `skip_pairtigs`: specify true or false whether you would like to
  assemble paired reads first before mapping onto the reference 16S rDNA
  database. This option is highly recomended for samples that are from
  full metagenomes that will likely be mostly from non-rDNA source
*	`minimum_pairtig_length`: Specify the minimum length that pairtigs must be.  This 
	option has no effect if the `pairtig_read_files` option is set. (default: 350)
*	`pair_overlap_length`: The minimum number of nucleotides that two reads from a 
	pair must overlap by to generate a pairtig. (default: 30)
*	`mapper`: The name of the short read aligner used in Amplishot.  Currently only
	bowtie2 is implemented and therefore the only valid value for this option is
	bowtie
*	`mapper_database`: Give the **full** filepath to an index file generated by the
	short read mapper
*	`taxonomy_file`: the name of the file containing a mapping between the reference
	sequences and their taxon strings
*	`mapping_similarity_cutoffs`: A list of required similarity between a reference sequence
	and a pairtig. Reads will be segregated into a band of similarity
    and assembled separately in that band 
*	`taxon_coverage`: list of two integer numbers that determine whether there are 
	enough reads for assembly.  The first number must be the minimum coverage (vertical
	read depth) for a taxon; the second number is the number of bases that must 
	contain the minimum coverage. (default: [2, 1000])
*	`assembly_method`: *de novo* 16S reconstruction method.  The only valid
	values are `phrap`, `ray` and  `velvet`
*	`minimum_reconstruction_length`: minimum length of sequences that are defined 
	as 'full length' and used in taxonomic assignment.
*	`otu_clustering_method`: Currently the only valid value is `cdhit`
*	`otu_clustering_similarity`: the similarity used for clustering full-length 
	sequences from different samples into OTUs 
* `neighbours_file`: A file that calculate the phylogenetic distanse
  between two separate reference sequences

### Program related blocks
Some of the underlying programs used in Amplishot can be controlled precisely by
specifying a *block* in the configuration file containing options specific 
to that program.  Each of these blocks is specified with a key that is identical
to the program name; within each block are program specific key-value pairs.  
The program specific key-value pairs must be indented by 4 spaces ( **not tabs** ),
this indentation must be consistent throughout the entire configuration file.
Currently program related blocks are available for both the assembly and 
taxonomy assignment parts of Amplishot

#### Assembly 

##### Phrap
 Specify extra options
using the `phrap` key.  Any of the command-line options available in phrap 
(listed [here](http://www.phrap.org/phredphrap/phrap.html)) can be used as the keys
in the phrap block, however you must not add in the dash (-) prefix for the options.
For example to modify stringency of the assembly, you could change the scoring matrix:

	phrap:
		penalty: -9
		gap_ext: -11
		gap_init: -12
		minscore: 350

Just because you can do this does not mean that you should unless you know exactly
what you are doing or are experimenting when Amplishot is producing sub-standard
results.  The scoring matrix and other assembly parameters used in phrap have 
already been altered to generate accurate 16S assemblies, so the default
settings should work well.


#### Taxonomy Assignment
Taxonomy assignment is handled in Amplishot after the reconstruction of full-length
16S sequences has occurred.  There are a number of different methods for taxonomic
assignment that include some of those available in Qiime 1.6.0. The taxonomy 
assignment method is determined from the Amplishot configuration file with the
`assign_taxonomy_method` key.  By default the Bowtie2 taxonomy assigner is used.
The valid values for each classifier are shown below:
- `bowtie` for bowtie2 based assigner
- `blast` for Qiime blastall based assigner 
- `rdp` for RDP classifier
- `mothur` for Mothur classifier

#### Configuration File options

For all taxonomy assigners a special block can be given in the configuration file
for specific options.  The key to this block must be the same as the value of 
the `assign_taxonomy_method` key.  For example to use the blast taxon assigner the
following code could be added into the configuration file:

    assign_taxonomy_method: blast
    blast:
        evalue: 1e-50
        blast_db: /full/path/to/blast/database

#### Options specific to all taxon assigners
- `id_to_taxonomy_fp`: **Full** path to file containing a mapping between reference
sequences and their respective taxon strings.  By default all taxon assigners will
use the  value of the `taxonomy_file` key.  This option should only be used if 
different reference sequence set is being used for taxonomy assignment 
- `reference_sequences_fp`: **Full** path to file containing reference sequences

##### Bowtie
- `index`: **Full** file path to bowtie2 formatted index file for the reference sequences.  
By default the bow tie taxonomic assigner will use the value of the `mapper_database` 
key, however a different database file can be accessed  here for taxonomic assignment
- `threads`: specify the number of threads that bowtie can use 
- `percentId`: The minimum percent identity that a representative sequence must map with, 
any sequence below this threshold will not be given a taxonomy

##### Blast
- `blast_db`: **Full** file path to the file containing the blast database that must
be formatted using the `formatdb` utility for nucleotide sequences.  Do **not** add in
the file extensions usually associated with blast databases. e.g. `.nsq`, `.nin` etc.
- `evalue`: The maximum allowable e-value allowed for a given match.  If no match can
be found below this score, then a representative sequence will not be given a
taxonomic assignment.

##### Mothur
- `Confidence`: Minimum allowed confidence score for taxonomic assignment

##### RDP
- `Confidence`: Minimum allowed confidence score allowed for taxonomic assignment
- `max_memory`: Set the maximum memory allowed for the RDP java virtual machine
- `training_data_properties_fp`: **Full** path to a file containing pre-compiled 
training data.  
This option is overridden if both the `reference_sequences_fp` and 
`id_to_taxonomy_fp` keys are set.

### Example Configuration file

    ---
    threads: 5
    log_level: INFO
    minimum_pairtig_length: 350 # minimum length of the overlapped pairs
    pair_overlap_length: 30 # mimimum length of the overlap
    mapper: bowtie # program used for read mapping 
    mapping_similarity_cutoffs: [0.85, 0.90, 0.95, 0.98] # the sequence similarity required between the reference database and the reads
    taxon_coverage: [2, 1000] # list of two numbers. The first is the minimum coverage, the second is the number of bases that need to be covered
    assembly_method: ray # choose a genome assembler  
    minimum_reconstruction_length: 1000 # minimum length of sequences that we define as 'full length'
    otu_clustering_method: cdhit
    otu_clustering_similarity: 0.97 # the similarity used for clustering full-length sequences from different samples into OTUs
    read_mapping_percent: 0.90 # the percent identity that individual reads have to map with to be considered part of the reference
    assign_taxonomy_method: blast
    minimum_taxon_similarity: 0.90 # sequences that fall below this cutoff will be listed as no taxonomy
    blast_db: '/srv/whitlam/bio/db/gg/from_www.secongenome.com/2012_10/gg_12_10_otus/rep_set/99_otus.fasta'

### Tips for writing config files
Writing out the full file path names in the configuration file can be a
real pain.  However you can reduce the burden on yourself by taking
advantage of some of the advanced features in the `vim` text editor.
When in `INSERT` mode if you start typing a file path (like `~/`) and
then press CTRL-x CTRL-f, you'll get a popup menu of file paths!! You
can use this to quickly add in file names to your config file. 
