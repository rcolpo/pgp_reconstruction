import argparse
import os
import os.path
import subprocess
import pickle
from datetime import datetime
from multiprocessing import freeze_support
import sys
import warnings

#enter to path two folders up. Important for the web application
try: current_path = os.path.dirname(os.path.realpath(__file__))
except: current_path = r'C:\phd\pgp\webApp\webApp\pgp_reconstruction\data\generated'
parent_path = os.path.abspath(os.path.join(current_path, '../..'))
sys.path.append(parent_path)


sys.path.append('/scratch/sbmlcomp/public2_html/webApp/webApp/python3.8-dependencies')
sys.path.append('/scratch/sbmlcomp/public2_html/webApp/webApp/python3.8-dependencies/lib/python3.8/site-packages/cplex-22.1.1.0-py3.8.egg')
sys.path.append('/scratch/sbmlcomp/public2_html/webApp/webApp/python3.8-dependencies/lib/python3.8/site-packages/cplex-22.1.1.0-py3.8.egg/cplex')
sys.path.append('/scratch/sbmlcomp/public2_html/webApp/webApp/python3.8-dependencies/lib/python3.8/site-packages')

from reframed import load_cbmodel
import cobra

import pgp_reconstruction
from pgp_reconstruction import config, project_dir
from pgp_reconstruction.reconstruction.findSoftConstraints import taxonomyBasedConstraints
from pgp_reconstruction.reconstruction.findOrfs import findOrfs
from pgp_reconstruction.reconstruction.prune_universal_model import prune_model
from pgp_reconstruction.reconstruction.scoring import reaction_scoring
from pgp_reconstruction.reconstruction.scoring import useReferenceModelData
from pgp_reconstruction.reconstruction.diamond import execute_diamond_blast, parse_diamond_output
from pgp_reconstruction.cli.download_missing_data import download_missing_files
from pgp_reconstruction.cli.util import saveProgressFile
from pgp_reconstruction.cli.util import first_run_check
from pgp_reconstruction.cli.util import loadConstraints


def maincall(inputFileName, outputfile=None, diamond_args=None, verbose=True, constraintsFilePath=None, reference=None, referenceScore=None, universalFile=None, gapfillMinimumMedia=None):


	if verbose: print('\nPreparing to reconstruct model. ' + str(datetime.now()) + '\n')

	if outputfile:
		model_id = os.path.splitext(os.path.basename(outputfile))[0]
	else:
		model_id = os.path.splitext(os.path.basename(inputFileName))[0]
		outputfile = os.path.splitext(inputFileName)[0] + '.xml'
		
	folder = os.path.split(inputFileName)[0]
	if folder: os.chdir(folder)
	

	outputfolder = os.path.abspath(os.path.dirname(outputfile))

	if not os.path.exists(outputfolder):
		try:
			os.makedirs(outputfolder)
		except:
			errorMessage = 'Unable to create output folder: ' + outputfolder
			raise ValueError(errorMessage)

	saveProgressFile(3, outputfolder)


	if universalFile:
		#open universal model providad by user
		reframedModel = load_cbmodel(universalFile, flavor='bigg')
		cobraModel = cobra.io.read_sbml_model(universalFile)
		
		#keep cobraModel annotation field always as a list
		for cobraObject in [cobraModel.reactions, cobraModel.metabolites]:
			for rxn in cobraObject:
				for db in rxn.annotation:
					if type(rxn.annotation[db]) == type(''):
						rxn.annotation[db] = [rxn.annotation[db]]
						
		#replace rxn.compartments by rxn.compartment.
		for rxn in cobraModel.reactions:
			compartmentSet = set()
			for met in rxn.metabolites:
				compartmentSet.add(met.compartment)
			rxn.compartment = list(compartmentSet)
	
	else:
		try:
			try:
				#try opening files saved as pickle. (much more efficient)
				pickle_file_path = os.path.join(project_dir, 'data/generated', 'cobraUniversalModel.pickle')
				with open(pickle_file_path, 'rb') as f:
					cobraModel = pickle.load(f)
					
				pickle_file_path = os.path.join(project_dir, 'data/generated', 'reframedUniversalModel.pickle')
				with open(pickle_file_path, 'rb') as f:
					reframedModel = pickle.load(f)
			
			except:
				#try opening universal model directely from its .XML file
				universe = os.path.join(project_dir, 'data/generated', 'universalRheaUnidirecional.xml')
				reframedModel = load_cbmodel(universe, flavor='bigg')
				cobraModel = cobra.io.read_sbml_model(universe)
				
				#keep cobraModel annotation field always as a list
				for cobraObject in [cobraModel.reactions, cobraModel.metabolites]:
					for rxn in cobraObject:
						for db in rxn.annotation:
							if type(rxn.annotation[db]) == type(''):
								rxn.annotation[db] = [rxn.annotation[db]]
								
				#replace rxn.compartments by rxn.compartment.
				for rxn in cobraModel.reactions:
					compartmentSet = set()
					for met in rxn.metabolites:
						compartmentSet.add(met.compartment)
					rxn.compartment = list(compartmentSet)
								
				with open(os.path.join(project_dir, 'data/generated', 'cobraUniversalModel.pickle'), 'wb') as handle:
					pickle.dump(cobraModel, handle, protocol=4)
				with open(os.path.join(project_dir, 'data/generated', 'reframedUniversalModel.pickle'), 'wb') as handle:
					pickle.dump(reframedModel, handle, protocol=4)
							
		except IOError:
			raise IOError(f'Failed to load universe model: {universe}\n')
		

	#load constraints
	if constraintsFilePath: constraintsFromFile = loadConstraints(constraintsFilePath, cobraModel, reframedModel)
	else: constraintsFromFile = {'soft':dict(),'hard':dict()}

	#create constraintsFromTaxonomy based on taxonomy
	constraintsFromTaxonomy, rxnsInTaxonomyConstraints, taxoOfTarget, rheaWithSameSyn, rxnsFromUniprot = taxonomyBasedConstraints(inputFileName, cobraModel)
	saveProgressFile(4, outputfolder)

	#se for DNA, roda prodigal para encontrar ORFs e traduzir as sequencias.
	inputfileNew, geneAndProteinNamePerSeqId = findOrfs(inputFileName)
	saveProgressFile(8, outputfolder)

	#parse_diamond_output
	filesList = os.listdir()
	diamond_db = project_dir + config.get('generated', 'diamond_db')
	
	if model_id + '-Diamond.tsv' in filesList: blast_output = model_id + '-Diamond.tsv'
	elif inputFileName + '-Diamond.tsv' in filesList: blast_output = inputFileName + '-Diamond.tsv'
	else: blast_output = model_id.split('-model')[0] + '-Diamond.tsv'
	
	if blast_output not in filesList or os.path.getsize(blast_output) == 0:
		
		if verbose:  print('Running diamond: ' + str(datetime.now()) + '\n')
		exit_code = execute_diamond_blast(inputfileNew, 'protein', blast_output, diamond_db, diamond_args, verbose)

		if exit_code is None:
			raise ValueError('Unable to run diamond (make sure diamond is available in your PATH).')

		if exit_code != 0:
			errorMessage = 'Failed to run diamond.'
			if diamond_args is not None: errorMessage += ' Incorrect diamond args? Please check documentation or use default args.'
			raise ValueError(errorMessage)
			
		if verbose:  print('Finished diamond: ' + str(datetime.now()) + '\n')
			
	saveProgressFile(30, outputfolder)

	diamondResult = parse_diamond_output(blast_output)
	
	rxnsScores, rheaIdToGene, bestMatchPerRead, singleRxnsInGenes = reaction_scoring(diamondResult, geneAndProteinNamePerSeqId, cobraModel, reframedModel, constraintsFromTaxonomy, constraintsFromFile, rxnsInTaxonomyConstraints, rxnsFromUniprot, outputfolder, verbose)
	
	saveProgressFile(40, outputfolder)
	
	if rxnsScores is None:
		raise ValueError('The input genome did not match sufficient genes/reactions in the database.')
	
	useReferenceModelData(reference, referenceScore, cobraModel, rxnsScores)
	
	if verbose: print('All in place! Starting to reconstruct model. ' + str(datetime.now()) + '\n')


	model = prune_model(reframedModel, cobraModel, rxnsScores, constraintsFromFile, rheaIdToGene, bestMatchPerRead, 
				taxoOfTarget, rheaWithSameSyn, singleRxnsInGenes, gapfillMinimumMedia, outputfolder, outputfile)

	if model is None:
		saveProgressFile("Failed to build model. ", outputfolder)
		raise ValueError("Failed to build model.")
		
	saveProgressFile(99, outputfolder)


def main():
	"""
	This is the main function that sets up the argument parser and initiates the metabolic model reconstruction process.
	"""
	
	parser = argparse.ArgumentParser(
		description="Reconstruct a metabolic model using 'Pathway-Guided Pruning Reconstruction'",
		formatter_class=argparse.RawTextHelpFormatter
	)

	parser.add_argument('input', metavar='INPUT', nargs='+', help="Input file. Acceptable formats: protein fasta, DNA fasta, GenBank, or Prokka annotation file.")

	parser.add_argument('--diamond-args', help="Additional arguments for running the diamond tool. (e.g., 'more-sensitive')")
	parser.add_argument('--constraints', help="Path to the constraints file (should be a .txt file).")
	parser.add_argument('--reference', help="Path to a manually curated model of a close reference species.")
	parser.add_argument('--referenceScore', default=0.1, help="Score given to reactions from reference model." )
	parser.add_argument('-o', '--output', dest='output', help="Path to the SBML output file.")
	parser.add_argument('-q', '--quiet', action='store_false', dest='verbose', help="Disable verbose mode.")
	parser.add_argument('--updateDB', action='store_true', help="Check for a more recent version of the databases used by the tool.")
	parser.add_argument('--universe-file', dest='universalFile', help="You can use your own universal model (SBML format)")
	parser.add_argument('--symbiosisDependent', action='store_false', help="Indicate if the organism can NOT grow while isolated.")

	args = parser.parse_args()

	# Check if diamond file exists in the data folder
	firstRun = first_run_check(args.updateDB)
	
	if firstRun:
		print('\n########\nThis was pgp_reconstruction first run. Files were included. Please start the application again for normal usage. If you keep seeing this message, manually download the missing files from:\n https://files.ufz.de/~umb-pgp_reconstruction-01/ \n########\n')
		return

	if len(args.input) > 1:
		parser.error('Can only accept one input per run. If your file name has spaces, try using using double quotes ( " ) instead of single quotes ( \' ), or replace the white space by underscore signs.')

	try:
		maincall(
			inputFileName=args.input[0],
			outputfile=args.output,
			diamond_args=args.diamond_args,
			verbose=args.verbose,
			constraintsFilePath=args.constraints,
			reference=args.reference,
			referenceScore=args.referenceScore,
			universalFile=args.universalFile,
			gapfillMinimumMedia=args.symbiosisDependent
		)
	except Exception as e:
		print(f"An error occurred: {e}", file=sys.stderr)
		sys.exit(1)


if __name__ == '__main__':
	freeze_support()
	main()

