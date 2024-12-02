import click
import os
import os.path
import sys
import shutil
import random
import string
import itertools
import time
import pysam
import math
import glob
import logging

from pkg_resources import get_distribution
from subprocess import call, check_call
from .mgatkHelp import *
from ruamel.yaml import YAML
from itertools import repeat
from ruamel.yaml.scalarstring import SingleQuotedScalarString as sqs
from multiprocessing import Pool

# Setup logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@click.command()
@click.version_option()
@click.argument('mode', type=click.Choice(['bcall', 'call', 'tenx', 'check', 'support', 'remove-background']))
@click.option('--input', '-i', default=".", required=True, help='Input; either directory of singular .bam file; see documentation. REQUIRED.')
@click.option('--output', '-o', default="mgatk_out", help='Output directory for analysis required for `call` and `bcall`. Default = mgatk_out')
@click.option('--name', '-n', default="mgatk", help='Prefix for project name. Default = mgatk')
@click.option('--mito-genome', '-g', default="rCRS", required=True, help='Mitochondrial genome configuration. Choose hg19, hg38, mm10, (etc.) or a custom .fasta file; see documentation. Default = rCRS.')
@click.option('--ncores', '-c', default=1, help='Number of cores to run the main job in parallel.')
@click.option('--barcode-tag', '-bt', default="X", help='Read tag (generally two letters) to separate single cells; valid and required only in `bcall` mode.')
@click.option('--barcodes', '-b', default="", help='File path to barcodes that will be extracted; useful only in `bcall` mode. If none supplied, mgatk will learn abundant barcodes from the bam file (threshold defined by the -mb tag).')
@click.option('--min-barcode-reads', '-mb', default=1000, help='Minimum number of mitochondrial reads for a barcode to be genotyped; useful only in `bcall` mode; will not overwrite the `--barcodes` logic. Default = 1000.')
@click.option('--NHmax', default=1, help='Maximum number of read alignments allowed as governed by the NH flag. Default = 1.')
@click.option('--NMmax', default=4, help='Maximum number of paired mismatches allowed represented by the NM/nM tags. Default = 4.')
@click.option('--remove-duplicates', '-rd', is_flag=True, help='Remove duplicate (presumably PCR) reads')
@click.option('--umi-barcode', '-ub', default="", help='Read tag (generally two letters) to specify the UMI tag when removing duplicates for genotyping.')
@click.option('--handle-overlap', '-ho', is_flag=True, help='Only count each base in the overlap region between a pair of reads once')
@click.option('--low-coverage-threshold', '-lc', default=10, help='Variant count for each cell will be ignored below this when calculating VMR')
@click.option('--max-javamem', '-jm', default="8000m", help='Maximum memory for java for running duplicate removal per core. Default = 8000m.')
@click.option('--proper-pairs', '-pp', is_flag=True, help='Require reads to be properly paired.')
@click.option('--base-qual', '-q', default=0, help='Minimum base quality for inclusion in the genotype count. Default = 0.')
@click.option('--alignment-quality', '-aq', default=0, help='Minimum alignment quality to include read in genotype. Default = 0.')
@click.option('--emit-base-qualities', '-eb', is_flag=True, help='Output mean base quality per alt allele as part of the final output.')
@click.option('--nsamples', '-ns', default=7000, help='The number of samples / cells to be processed per iteration; Default = 7000. Supply 0 to try all.')
@click.option('--keep-samples', '-k', default="ALL", help='Comma separated list of sample names to keep; ALL (special string) by default. Sample refers to basename of .bam file')
@click.option('--ignore-samples', '-x', default="NONE", help='Comma separated list of sample names to ignore; NONE (special string) by default. Sample refers to basename of .bam file')
@click.option('--keep-temp-files', '-z', is_flag=True, help='Add this flag to keep all intermediate files.')
@click.option('--keep-qc-bams', '-qc', is_flag=True, help='Add this flag to keep the quality-controlled bams after processing.')
@click.option('--skip-R', '-sr', is_flag=True, help='Generate plain-text only output. Otherwise, this generates a .rds object that can be immediately read into R for downstream analysis.')
@click.option('--snake-stdout', '-so', is_flag=True, help='Write snakemake log to stdout rather than a file. May be necessary for certain HPC environments.')

def main(mode, input, output, name, mito_genome, ncores,
	barcode_tag, barcodes, min_barcode_reads,
	nhmax, nmmax, remove_duplicates, umi_barcode, handle_overlap, low_coverage_threshold,
	max_javamem, proper_pairs, base_qual, alignment_quality, emit_base_qualities,
	nsamples, keep_samples, ignore_samples,
	     keep_temp_files, keep_qc_bams, skip_r, snake_stdout):
	
	"""
	mgatk: a mitochondrial genome analysis toolkit. \n
	MODE = ['bcall', 'call', 'tenx'] \n
	See https://github.com/caleblareau/mgatk/wiki for more details.
	"""
	
	script_dir = os.path.dirname(os.path.realpath(__file__))
	cwd = os.getcwd()
	__version__ = get_distribution('mgatk').version
	logger.info("mgatk version %s" % __version__)
	
	# Verify dependencies	
	if remove_duplicates:
		check_software_exists("java")
	
	if (mode in ["call", "tenx", "bcall"]) and not skip_r:
		check_software_exists("R")
		check_R_packages(["data.table", "SummarizedExperiment", "GenomicRanges", "Matrix"])
		check_pip_packages(['matplotlib', 'pulp==2.7.0'])


	ncores = ncores if isinstance(ncores, int) else int(ncores)
	
	# Determine which genomes are available
	rawsg = glob.glob(script_dir + "/bin/anno/fasta/*.fasta")
	supported_genomes = [x.replace(script_dir + "/bin/anno/fasta/", "").replace(".fasta", "") for x in rawsg]  
		
	if mode in ["bcall", "tenx", "check"]:
		if barcode_tag == "X":
			logger.error('In `'+mode+'` mode, must specify a valid read tag ID (generally two letters).')
			sys.exit('ERROR: in `'+mode+'` mode, must specify a valid read tag ID (generally two letters).')
			
		# Input argument is assumed to be a .bam file
		filename, file_extension = os.path.splitext(input)
		if file_extension != ".bam":
			logger.error('In `'+mode+'` mode, the input should be an individual .bam file.')
			sys.exit('ERROR: in `'+mode+'` mode, the input should be an individual .bam file.')
		if not os.path.exists(input):
			logger.error('No file found called "' + input + '"; please specify a valid .bam file.')
			sys.exit('ERROR: No file found called "' + input + '"; please specify a valid .bam file.')
		if not os.path.exists(input + ".bai"):
			logger.info("Attempting to index: " + input)
			pysam.index(input)
			
			if not os.path.exists(input + ".bai"):
				logger.error('Index your input .bam file for `bcall` mode.')
				sys.exit('ERROR: index your input .bam file for `bcall` mode.')
		
		# Determine whether or not we have been supplied barcodes
		barcode_known = False
		if os.path.exists(barcodes) and barcodes != "":
			barcode_known = True
		else:
			if mode == "tenx":
				logger.error('Must specify a known barcode list with `tenx` mode')
				sys.exit(gettime() + 'Must specify a known barcode list with `tenx` mode')
			logger.info("Will determine barcodes with at least: " + str(min_barcode_reads) + " mitochondrial reads.")

	# Log the configuration details as a table
	config_details = [
		("Mode", mode),
		("Input", input),
		("Output", output),
		("Name", name),
		("Mitochondrial Genome", mito_genome),
		("Number of Cores", ncores),
		("Barcode Tag", barcode_tag),
		("Barcodes", barcodes),
		("Min Barcode Reads", min_barcode_reads),
		("NHmax", nhmax),
		("NMmax", nmmax),
		("Remove duplicates", remove_duplicates),
		("UMI Barcode", umi_barcode),
		("Handle Overlap", handle_overlap),
		("Low Coverage Threshold", low_coverage_threshold),
		("Max Java Memory", max_javamem),
		("Proper Pairs", proper_pairs),
		("Base Quality", base_qual),
		("Alignment Quality", alignment_quality),
		("Emit Base Qualities", emit_base_qualities),
		("Number of Samples", nsamples),
		("Keep Samples", keep_samples),
		("Ignore Samples", ignore_samples),
		("Keep Temp Files", keep_temp_files),
		("Keep QC Bams", keep_qc_bams),
		("Skip R", skip_r),
		("Snake Stdout", snake_stdout)
	]

	logger.info("Configuration for analysis:")
	for key, value in config_details:
		logger.info(f"{key:25}: {value}")
		
	# Remember that I started off as bcall as this will become overwritten
	wasbcall = False
	continue_check = True
	if mode in ["bcall", "tenx", "check"]:
					
		# Make temporary directory of inputs
		of = output; tf = of + "/temp"; bcbd = tf + "/barcoded_bams" # bcdb = barcoded bam directory
		folders = [of, tf, bcbd, of + "/final"]
		mkfolderout = [make_folder(x) for x in folders]
		# logger.info("Created necessary folders for processing: " + str(folders))
		
		# Handle fasta requirements
		fastaf, mito_chr, mito_length = handle_fasta_inference(mito_genome, supported_genomes, script_dir, mode, of)
		logger.info(f"Using fasta file: {fastaf}")
		idxs = pysam.idxstats(input).split("\n")
		
		# Handle common mtDNA reference genome errors
		bam_length = 0
		for i in idxs:
			if(i.split("\t")[0] == mito_chr):
				bam_length = int(i.split("\t")[1])
		
		if(mito_length == bam_length):
			pass
		#	logger.info("User specified mitochondrial genome matches .bam file")
		elif(bam_length == 16569):
		#	logger.info("User specified mitochondrial genome does NOT match .bam file; using rCRS instead (length == 16569)")
			fastaf, mito_chr, mito_length = handle_fasta_inference("rCRS", supported_genomes, script_dir, mode, of)
		elif(bam_length == 16571):
		#	logger.info("User specified mitochondrial genome does NOT match .bam file; using hg19 instead (length == 16571)")
			fastaf, mito_chr, mito_length = handle_fasta_inference("hg19", supported_genomes, script_dir, mode, of)
		else:
			sys.exit("User specified mitochondrial genome does NOT match .bam file; correctly specify reference genome or .fasta file")
		
		# Actually call the external script based on user input
		if(not barcode_known):
			barc_quant_file = of + "/final/barcodeQuants.tsv"
			passing_barcode_file = of + "/final/passingBarcodes.tsv"
			find_barcodes_py = script_dir + "/bin/python/find_barcodes.py"
			
			pycall = " ".join(['python', find_barcodes_py, input, bcbd, barcode_tag, str(min_barcode_reads), mito_chr, barc_quant_file, passing_barcode_file])
			logger.info(f"Running barcode finding script: {pycall}")
			os.system(pycall)
			barcodes = passing_barcode_file
			logger.info(f"Generated barcode files: {barc_quant_file}, {passing_barcode_file}")
		
		# logger.info("Finished processing barcodes for genotyping.")
		if(mode == "check" or mode == "tenx"):
			barcode_files = split_barcodes_file(barcodes, math.ceil(file_len(barcodes)/ncores), output)
		#	logger.info(f"Split barcodes into files for check/tenx mode: {barcode_files}")
			samples = [os.path.basename(os.path.splitext(sample)[0]) for sample in barcode_files] 
			samplebams = [of + "/temp/barcoded_bams/" + sample + ".bam" for sample in samples]
		#	logger.info(f"Sample BAM files: {samplebams}")
			
			if(umi_barcode == ""):
				umi_barcode = "XX"
			
			# Enact the split in a parallel manner if necessary
			if(mode == "tenx"):
				pool = Pool(processes=ncores)
				pmblah = pool.starmap(split_chunk_file, zip(barcode_files, repeat(script_dir), repeat(input), repeat(bcbd), repeat(barcode_tag), repeat(mito_chr), repeat(umi_barcode)))
				pool.close()
			#	logger.info("Completed parallel processing for tenx mode.")
				
			umi_barcode = "MU"
			continue_check = False
		
		# logger.info("Finished determining/splitting barcodes for genotyping.")

		# List the split BAM files and files created in the previous steps
		# logger.info("Listing split BAM files and files created in the previous steps:")
		# split_bam_files = glob.glob(bcbd + "/*.bam")
		# for bam_file in split_bam_files:
		#	logger.info(f"Split BAM file: {bam_file}")
		
		# created_files = glob.glob(of + "/final/*")
		# for created_file in created_files:
		#	logger.info(f"Created file: {created_file}")
			
	# -------------------------------
	# Determine samples for analysis
	# -------------------------------

	if(mode == "tenx"):
	
		# Make all of the output folders if necessary
		of = output; tf = of + "/temp"; qc = of + "/qc"; logs = of + "/logs"
		folders = [logs, of + "/logs/filterlogs", of + "/fasta", of + "/.internal",
			 of + "/.internal/parseltongue", of + "/.internal/samples", of + "/final", 
			 tf, tf + "/ready_bam", tf + "/temp_bam", tf + "/sparse_matrices", tf + "/quality",
			 qc, qc + "/quality", qc + "/depth"]

		mkfolderout = [make_folder(x) for x in folders]
		# logger.info("Created necessary folders for processing: " + str(folders))
		
		#-------------------
		# Handle .fasta file
		#-------------------			
		if(mode == "tenx"):
			# Logging		
			logf = open(output + "/logs" + "/base.mgatk.log", 'a')
		#	logger.info("Starting analysis with mgatk")
			
			if(mode == "call"):
				logger.info(nsamplesNote)

		if (remove_duplicates):
			make_folder(of + "/logs/rmdupslogs")
		#	logger.info("Created folder for remove duplicates logs")
	
		# Create internal README files 
		if not os.path.exists(of + "/.internal/README"):
			with open(of + "/.internal/README" , 'w') as outfile:
				outfile.write("This folder creates important (small) intermediate; don't modify it.\n\n")
		#	logger.info("Created .internal/README file")
		if not os.path.exists(of + "/.internal/parseltongue/README"):	
			with open(of + "/.internal/parseltongue/README" , 'w') as outfile:
				outfile.write("This folder creates intermediate output to be interpreted by Snakemake; don't modify it.\n\n")
		#	logger.info("Created .internal/parseltongue/README file")
		if not os.path.exists(of + "/.internal/samples/README"):
			with open(of + "/.internal" + "/samples" + "/README" , 'w') as outfile:
				outfile.write("This folder creates samples to be interpreted by Snakemake; don't modify it.\n\n")
		#	logger.info("Created .internal/samples/README file")
	
		# Set up sample bam plain text file
		for i in range(len(samples)):
			with open(of + "/.internal/samples/" + samples[i] + ".bam.txt" , 'w') as outfile:
				outfile.write(samplebams[i])
		#	logger.info(f"Created sample bam text file: {of + '/.internal/samples/' + samples[i] + '.bam.txt'}")
		
		# logger.info(f"Genotyping samples with {ncores} threads")
		
		# add sqs to get .yaml to play friendly https://stackoverflow.com/questions/39262556/preserve-quotes-and-also-add-data-with-quotes-in-ruamel
		dict1 = {'input_directory' : sqs(input), 'output_directory' : sqs(output), 'script_dir' : sqs(script_dir),
			'fasta_file' : sqs(fastaf), 'mito_chr' : sqs(mito_chr), 'mito_length' : sqs(mito_length), 'name' : sqs(name),
			'base_qual' : sqs(base_qual), 'remove_duplicates' : sqs(remove_duplicates), 'handle_overlap' : sqs(handle_overlap),
			'low_coverage_threshold' : sqs(low_coverage_threshold), 'barcode_tag' : sqs(barcode_tag), 'umi_barcode' : sqs(umi_barcode),
			'alignment_quality' : sqs(alignment_quality), 'emit_base_qualities' : sqs(emit_base_qualities),
			'proper_paired' : sqs(proper_pairs), 'NHmax' : sqs(nhmax), 'NMmax' : sqs(nmmax), 'max_javamem' : sqs(max_javamem)}
				
		y_s = of + "/.internal/parseltongue/snake.scatter.yaml"
		with open(y_s, 'w') as yaml_file:
			yaml=YAML()
			yaml.default_flow_style = False
			yaml.dump(dict1, yaml_file)
		# logger.info(f"Created YAML configuration file: {y_s}")
		
		cp_call = "cp " + y_s +  " " + logs + "/" + name + ".parameters.txt"
		os.system(cp_call)
		# logger.info(f"Copied YAML configuration to parameters file: {logs + '/' + name + '.parameters.txt'}")

		# Confirm with user to proceed
		# click.confirm(click.style('Proceed to genotyping?\n', fg='yellow'), abort=True)

		if(mode == "tenx"):
			import subprocess

			# Execute snakemake
			snake_log = logs + "/" + name + ".snakemake_tenx.log"

			snake_log_out = ""
			if not snake_stdout:
				snake_log_out = f' > {snake_log} 2>&1'

			snakecmd_tenx = f'snakemake --snakefile {script_dir}/bin/snake/Snakefile.tenx --cores {ncores} --config cfp="{y_s}" {snake_log_out}'

			# Log the Snakemake command with f-strings
			logger.info(f"Executing Snakemake command: {snakecmd_tenx}")

			# Execute the Snakemake command and capture output and errors
			try:
				result = subprocess.run(snakecmd_tenx, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			#	logger.info(f"Snakemake output: {result.stdout.decode('utf-8')}")
			#	logger.info(f"Snakemake errors: {result.stderr.decode('utf-8')}")
			except subprocess.CalledProcessError as e:
				logger.info(f"Snakemake failed with error: {str(e)}")
				logger.info(f"Snakemake output: {e.stdout.decode('utf-8')}")
				logger.info(f"{gettime()}Snakemake errors: {e.stderr.decode('utf-8')}")
				sys.exit(1)

	#-------
	# Gather
	#-------

	if(mode == "tenx"):

		# Make .rds file from the output
		Rcall = "Rscript " + script_dir + "/bin/R/toRDS.R " + output + "/final " + name
		os.system(Rcall)
		logger.info("Successfully created final output files")

	#--------
	# Cleanup
	#--------
	if(mode == "tenx"):
		if keep_qc_bams:
			logger.info("Final bams retained since --keep-qc-bams was specified.")
			dest = shutil.move(of + "/temp/ready_bam", of + "/qc_bam")    

		if keep_temp_files:
			logger.info("Temporary files not deleted since --keep-temp-files was specified.")
		else:
			shutil.rmtree(of+ "/fasta")
			shutil.rmtree(of + "/.internal")
			shutil.rmtree(of + "/temp")
		#	logger.info("Intermediate files successfully removed.")

		# Suspend logging
		logf.close()