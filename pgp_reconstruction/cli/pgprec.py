from pgp_reconstruction import config, project_dir
from pgp_reconstruction.reconstruction.findSoftConstraints import taxonomyBasedConstraints
from pgp_reconstruction.reconstruction.findOrfs import findOrfs
from pgp_reconstruction.reconstruction.prune_universal_model import prune_model
from pgp_reconstruction.reconstruction.scoring import reaction_scoring
from pgp_reconstruction.reconstruction.diamond import run_blast, load_diamond_results
from pgp_reconstruction.cli import download_missing_files
from reframed import load_cbmodel
from reframed.io.sbml import sanitize_id
import argparse
import cobra
import os
import os.path
from multiprocessing import Pool
from glob import glob
import subprocess
import pickle
import datetime
import sys
from multiprocessing import freeze_support

def first_run_check():

	download_missing_files()

	diamond_db = project_dir + config.get('generated', 'diamond_db')
	if not os.path.exists(diamond_db):
		print("Running diamond for the first time, please wait while we build the internal database...")
		fasta_file = project_dir + config.get('generated', 'fasta_file')
		cmd = ['diamond', 'makedb', '--in', fasta_file, '-d', diamond_db[:-5]]
		try:
			exit_code = subprocess.call(cmd)
		except OSError:
			raise ValueError('Unable to run diamond with the command "' + cmd + '"\nMake sure diamond is installed and available in your PATH.')
		else:
			if exit_code != 0:
				raise ValueError('Failed to run diamond (wrong arguments).')


def loadConstraints(soft, constraintsFromFile, cobraModel, reframedModel, sufix=''):
	#receive a file with ID of reactions and metabolites. IDs might be from ChEBI, KEGG, BIGG, SEED or Rhea. 
	#finds IDs in universal model, and asign a score
	
	#create copy of annotation with ID lower case 
	cobraModelRxns = dict()
	for rxn in cobraModel.reactions:
		cobraModelRxns[rxn.id.lower()] = dict()
		for db in rxn.annotation:
			cobraModelRxns[rxn.id.lower()][db] = set()
			for annotationId in rxn.annotation[db]:
				cobraModelRxns[rxn.id.lower()][db].add(annotationId.lower())
	
	#create copy of annotation with ID lower case 
	cobraModelMets = dict()
	for met in cobraModel.metabolites:
		cobraModelMets[met.id.lower()] = dict()
		for db in met.annotation:
			cobraModelMets[met.id.lower()][db] = set()
			for annotationId in met.annotation[db]:
				cobraModelMets[met.id.lower()][db].add(annotationId.lower())		
	
	#constraintsFromFile = {'soft':dict(),'hard':dict()}
	constraintsDict = {'metabolites':{'soft':dict(),'hard':dict()}, 'reactions':{'soft':dict(),'hard':dict()}}
	with open(soft) as file:
		for line in file:
			line = line.strip()
			if not line: continue
			lineSplit = line.lower().split('\t')
			#error handling
			if len(lineSplit) != 4 and len(lineSplit) != 3:
				sys.exit('ERROR: Constraint regarding metabolites should have 4 columns, separated by tabs. For example: "M_na    Soft    1    Media". Constraint regarding reactions should have 3 columns, separated by tabs. For example: "R_R09640    Soft    3".')
			if len(lineSplit[0]) <= 2:
				sys.exit('ERROR: First column should contain the IDs of reactions or metabolites. If reaction, it should start with "R_"; If metabolite, it should start with "M_". We accept IDs from BIGG, SEED and ChEBI. For example: "M_na".\nInput given: ' + str(line))
			if lineSplit[0][:2] == 'm_':
				if len(lineSplit) != 4:
					sys.exit('ERROR: Constraint regarding metabolites should have 4 columns, separated by tabs. For example: "M_na    Soft    1    Media".\nInput given: ' + str(line))
				if lineSplit[3] != 'media' and lineSplit[3] != 'product':
					sys.exit('ERROR: It should be informed if the constraint is regarding the consuption (using keyword "Media") or production (using keyword "Product"). For example: "M_na    Soft    1    Media".\nInput given: ' + str(line))
			if lineSplit[0][:2] == 'r_' and len(lineSplit) != 3:
				sys.exit('ERROR: Constraint regarding reactions should have 3 columns, separated by tabs. For example: "R_R09640    Soft    3".\nInput given: ' + str(line))
			if lineSplit[1] != 'soft' and lineSplit[1] != 'hard':
				sys.exit('ERROR: In the second columns, it should be informed is the constraint is "soft" (model will try to satisfy the condition) or "hard" (model will necessarely satisfy the condition). For example: "R_R09640    Soft    3".\nInput given: ' + str(line))
			#checking if score is number
			try: a = float(lineSplit[2])
			except: sys.exit('ERROR: In the third column it should be informed the reaction score in the optimization problem. If metabolite, the score will be atributed to the reactions producing the metabolite. If hard constraint, only the score sign will be taken into consideration. For example: "R_R09640    Soft    3".\nInput given: ' + str(line))
			
				
			#store data
			if lineSplit[0][:2] == 'm_':
				if lineSplit[1] == 'soft': constraintsDict['metabolites']['soft'][lineSplit[0][2:]] = {'score':float(lineSplit[2]), 'mediaOrProduct': lineSplit[3]}
				elif lineSplit[1] == 'hard': constraintsDict['metabolites']['hard'][lineSplit[0][2:]] = {'score':float(lineSplit[2]), 'mediaOrProduct': lineSplit[3]}
			elif lineSplit[0][:2] == 'r_':
				if lineSplit[1] == 'soft': constraintsDict['reactions']['soft'][lineSplit[0][2:]] = {'score':float(lineSplit[2])}
				elif lineSplit[1] == 'hard': constraintsDict['reactions']['hard'][lineSplit[0][2:]] = {'score':float(lineSplit[2])}
	
	
	
	##for mets
	#remove container from met ID and create a dict with the original ID. 
	constraintsMets = dict()
	for constraintType in constraintsDict['metabolites']:
		for metId in constraintsDict['metabolites'][constraintType]:
			if metId[-2:] == '_e' or metId[-2:] == '_c': metId2 = metId[2:-2]
			else: metId2 = metId
			constraintsMets[metId2] = metId
	
	sink = dict()
	exchange = dict()
	changeBound = set()
	#identify reactions producing (for media) and consuming (for product) the metabolite.
	for metId2 in constraintsMets:
		
		if constraintsMets[metId2] + '_e' in cobraModel.metabolites:
			met = cobraModel.metabolites.get_by_id(constraintsMets[metId2] + '_e')
			toIterate = [met]
		elif constraintsMets[metId2] + '_c' in cobraModel.metabolites:
			met = cobraModel.metabolites.get_by_id(constraintsMets[metId2] + '_c')
			toIterate = [met]
		else:
			toIterate = cobraModel.metabolites
	
		for met in toIterate:
			for db in ['bigg', 'chebi', 'kegg', 'seed']:
				#if db in met.annotation and metId in met.annotation[db]:
				if db in cobraModelMets[met.id.lower()] and metId2 in cobraModelMets[met.id.lower()][db]:
					
					for rxn in met.reactions:
					
						#if met.compartment != 'e' and not rxn.id.startswith('SK'): continue
					
						if len(rxn.metabolites) == 1:
						
							if met.compartment == 'c' and rxn.id.startswith('SK'):
								if met.id[:-2] not in sink: sink[met.id[:-2]] = set()
								sink[met.id[:-2]].add(rxn.id+sufix)
							else:
								if met.id[:-2] not in exchange: exchange[met.id[:-2]] = set()
								exchange[met.id[:-2]].add(rxn.id+sufix)
						
							if met in rxn.products:
								for constraintType in ['soft', 'hard']:
									if constraintsMets[metId2] in constraintsDict['metabolites'][constraintType]:
										if rxn.upper_bound > 0 and constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['mediaOrProduct'] == 'media': 
											constraintsFromFile['forward'][constraintType][rxn.id+ sufix] = constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['score']
										if rxn.lower_bound < 0 and constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['mediaOrProduct'] == 'product': 
											constraintsFromFile['reverse'][constraintType][rxn.id+ sufix] = constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['score']
										if met.compartment == 'c' and rxn.id.startswith('SK'):
											if constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['mediaOrProduct'] == 'media': 
												if rxn.upper_bound <= 0: changeBound.add(rxn.id)
												else: constraintsFromFile['forward'][constraintType][rxn.id+ sufix] = constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['score']
											if constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['mediaOrProduct'] == 'product':
												if rxn.lower_bound >= 0: changeBound.add(rxn.id)
												else: constraintsFromFile['reverse'][constraintType][rxn.id+sufix] = constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['score']
											
							elif met in rxn.reactants: 
								for constraintType in ['soft', 'hard']:
									if constraintsMets[metId2] in constraintsDict['metabolites'][constraintType]:
										if rxn.lower_bound < 0 and constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['mediaOrProduct'] == 'media': 
											constraintsFromFile['reverse'][constraintType][rxn.id+ sufix] = constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['score']
										if rxn.upper_bound > 0 and constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['mediaOrProduct'] == 'product':
											
											constraintsFromFile['forward'][constraintType][rxn.id+ sufix] = constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['score']
										if met.compartment == 'c' and rxn.id.startswith('SK'):
											if constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['mediaOrProduct'] == 'media': 
												if rxn.lower_bound >= 0: 
													if constraintType == 'hard' and constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['score'] < 0: continue # it is asking to do something that the model alreday can not do
													changeBound.add(rxn.id)
												else: constraintsFromFile['reverse'][constraintType][rxn.id+ sufix] = constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['score']
											if constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['mediaOrProduct'] == 'product':
												if rxn.upper_bound <= 0: 
													if constraintType == 'hard' and constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['score'] < 0: continue # it is asking to do something that the model alreday can not do
													changeBound.add(rxn.id)
												else: constraintsFromFile['forward'][constraintType][rxn.id+ sufix] = constraintsDict['metabolites'][constraintType][constraintsMets[metId2]]['score']

	for direction in ['forward', 'reverse']:
		#remove sink option if there is exchange
		for metId in sink:
			if metId in exchange:
				for rxnId in sink[metId]:
					if rxnId in constraintsFromFile[direction]['soft']: 
						del constraintsFromFile[direction]['soft'][rxnId]
						changeBound.discard(rxnId)
					if rxnId in constraintsFromFile[direction]['hard']: 
						del constraintsFromFile[direction]['hard'][rxnId]
						changeBound.discard(rxnId)
					

	#change sink rxns bounds if necessary, to make it with go to the direction asked in the constraint
	for cobraId in changeBound:
		reframedId = 'R_'+cobraId.replace('-','__45__').replace('.','__46__').replace('+','__43__')
		rxn = cobraModel.reactions.get_by_id(cobraId)
		if rxn.lower_bound == 0:
			rxn.lower_bound = -1000
			rxn.upper_bound = 0
			reframedModel.reactions[reframedId].lb = -1000
			reframedModel.reactions[reframedId].ub = 0
		else:
			rxn.lower_bound = 0
			rxn.upper_bound = 1000
			reframedModel.reactions[reframedId].lb = 0
			reframedModel.reactions[reframedId].ub = 1000
						
	##for reactions
	for direction in ['forward', 'reverse']:
		constraintsDictForAnnotation = {'soft':set(),'hard':set()}
		for constraintType in ['soft', 'hard']:
			for rxnId in constraintsDict['reactions'][constraintType]:
				sucess = 0
				if rxnId in cobraModel.reactions: 
					constraintsFromFile[direction][constraintType][rxnId[2:]+ sufix] = constraintsDict['reactions'][constraintType][rxnId]['score']
					sucess = 1
				if rxnId in cobraModel.reactions: 
					constraintsFromFile[direction][constraintType][rxnId+ sufix] = constraintsDict['reactions'][constraintType][rxnId]['score']
					sucess = 1
				if rxnId in cobraModel.reactions: 
					constraintsFromFile[direction][constraintType][rxnId+ sufix] = constraintsDict['reactions'][constraintType][rxnId]['score']
					sucess = 1
				if sucess == 0: constraintsDictForAnnotation[constraintType].add(rxnId)

		for rxn in cobraModel.reactions:
			for db in ['bigg', 'chebi', 'kegg', 'seed']:
				if db not in rxn.annotation: continue 
				for constraintType in ['soft', 'hard']:
					for rxnId in constraintsDictForAnnotation[constraintType]:
						if rxnId in cobraModelRxns[rxn.id.lower()][db]:
							constraintsFromFile[direction][constraintType][rxn.id+ sufix] = constraintsDict['reactions'][constraintType][rxnId]['score']


def maincall(inputFileName, outputfile=None, diamond_args=None,
		 verbose=True, constraints=None, reference=None):

	
	if verbose: print('\nPreparing to reconstruct model. ' + str(datetime.datetime.now()) + '\n')
		
	
	pickle_file_path = os.path.join(project_dir, 'data/generated', 'uniprotToRheaRxns.pickle')
	with open(pickle_file_path, 'rb') as f:
		uniprotToRheaRxns = pickle.load(f)

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
			print('Unable to create output folder:', outputfolder)
			return


	try:
		try:
			#try opening files saved as pickle. (more efficient)
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
							
			with open('c:/Users/colpoama/AppData/Local/Programs/Python/Python37/Lib/site-packages/pgp_reconstruction/data/generated/cobraUniversalModel.pickle', 'wb') as handle:
				pickle.dump(cobraModel, handle, protocol=4)
			with open('c:/Users/colpoama/AppData/Local/Programs/Python/Python37/Lib/site-packages/pgp_reconstruction/data/generated/reframedUniversalModel.pickle', 'wb') as handle:
				pickle.dump(reframedModel, handle, protocol=4)
						
	except IOError:
		raise IOError(f'Failed to load universe model: {universe}\n')
		
		

	#load constraints
	constraintsFromFile = {'forward':{'soft':dict(),'hard':dict()},'reverse':{'soft':dict(),'hard':dict()}}
	if constraints: loadConstraints(constraints, constraintsFromFile, cobraModel, reframedModel)

	#create soft_constraints based on taxonomy
	soft_constraints, rxnsInSoftSyn, taxoOfTarget, rheaWithSameSyn = taxonomyBasedConstraints(inputFileName, cobraModel)

	#se for DNA, roda prodigal para encontrar ORFs e traduzir as sequencias.
	inputfileNew, geneAndProteinNamePerSeqId = findOrfs(inputFileName)

	filesList = os.listdir()


	#load_diamond_results
	if verbose: print('Running diamond...')
	diamond_db = project_dir + config.get('generated', 'diamond_db')
	
	if model_id + '-Diamond.tsv' in filesList: blast_output = model_id + '-Diamond.tsv'
	else: blast_output = model_id_formated + '-Diamond.tsv'
	
	if blast_output not in filesList or os.path.getsize(blast_output) == 0:
		exit_code = run_blast(inputfileNew, 'protein', blast_output, diamond_db, diamond_args, verbose)

		if exit_code is None:
			print('Unable to run diamond (make sure diamond is available in your PATH).')
			return

		if exit_code != 0:
			print('Failed to run diamond.')
			if diamond_args is not None:
				print('Incorrect diamond args? Please check documentation or use default args.')
			return

	diamondResult = load_diamond_results(blast_output)
	
	rxnsScores, rheaIdToGene, rxnsCobraToreframedModel, soft_constraints_pathways, bestPerReadSimplified, softPositiveNewRxn, relevantRxnsPerGeneInModel = reaction_scoring(diamondResult, geneAndProteinNamePerSeqId, cobraModel, reframedModel, uniprotToRheaRxns, soft_constraints, rxnsInSoftSyn, verbose)
	
	if rxnsScores is None:
		print('The input genome did not match sufficient genes/reactions in the database.')
		return
	

	#include in rxnsScores reactions from reference model
	if reference:
		try:
			reframedModelReference = cobra.io.read_sbml_model(reference)
			
			for cobraObject in [reframedModelReference.reactions, reframedModelReference.metabolites]:
				for rxn in cobraObject:
					for db in rxn.annotation:
						if type(rxn.annotation[db]) == type(''):
							rxn.annotation[db] = [rxn.annotation[db]]
			
		except IOError:
			raise IOError(f'Failed to load reference model')
		
		rxnWithGenes = set()
		rxnBoundary = set()
		rxnWithoutGenes = set()
		
		for rxn in reframedModelReference.reactions:

			if rxn.genes:
				for db in rxn.annotation:
					for synId in rxn.annotation[db]: rxnWithGenes.add(synId)
				rxnWithGenes.add(rxn.id)
					
			elif rxn in reframedModelReference.boundary:
				for db in rxn.annotation:
					for synId in rxn.annotation[db]: rxnBoundary.add(synId)
				rxnBoundary.add(rxn.id)

			else:
				for db in rxn.annotation:
					for synId in rxn.annotation[db]: rxnWithoutGenes.add(synId)
				rxnWithoutGenes.add(rxn.id)
					
		rxnsFromReference = {'rxnWithGenes':rxnWithGenes, 'rxnBoundary':rxnBoundary, 'rxnWithoutGenes':rxnWithoutGenes}
		
		print('\n')
		print('len(rxnWithGenes) = ' + str(len(rxnWithGenes)))
		print('len(rxnBoundary) = ' + str(len(rxnBoundary)))
		print('len(rxnWithoutGenes) = ' + str(len(rxnWithoutGenes)))
		print('\n')
		
	else: rxnsFromReference = {'rxnWithGenes':set(), 'rxnBoundary':set(), 'rxnWithoutGenes':set()}

	
	if verbose: print('All in place! Starting to reconstruct model. ' + str(datetime.datetime.now()) + '\n')


	model = prune_model(reframedModel, cobraModel, rxnsScores, rheaIdToGene, uniprotToRheaRxns, inputFileName, rxnsCobraToreframedModel, soft_constraints_pathways, bestPerReadSimplified, 
				taxoOfTarget, rheaWithSameSyn, softPositiveNewRxn, constraintsFromFile, rxnsFromReference, relevantRxnsPerGeneInModel)

	if model is None:
		print("Failed to build model.")
		return


def main():

	parser = argparse.ArgumentParser(description="Reconstruct a metabolic model using 'Pathway-Guided Pruning Reconstruction'",
									 formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('input', metavar='INPUT', nargs='+',
						help="Input (protein fasta file, dna fasta file, GenBank file and Prokka annotation file can be used as input.\n"
						)

	parser.add_argument('--diamond-args', help="Additional arguments for running diamond")
	parser.add_argument('-o', '--output', dest='output', help="SBML output file (or output folder if -r is used)")
	parser.add_argument('-v', '--verbose', action='store_true', dest='verbose', help="Switch to verbose mode")
	parser.add_argument('--constraints', help="Constraints file")
	parser.add_argument('--reference', help="Manually curated model of a close reference species.")

	args = parser.parse_args()
		
	
	#check if diamond file exists on data folder
	first_run_check()

	if len(args.input) > 1:
		parser.error('Can only accept one input per run. If your file name has spaces, try using using double quotes ( " ) instead of single quotes ( \' ), or replace the white space by underscore signs.')

	maincall(
		inputFileName=args.input[0],
		outputfile=args.output,
		diamond_args=args.diamond_args,
		verbose=args.verbose,
		constraints=args.constraints,
		reference=args.reference,
	)


if __name__ == '__main__':
	freeze_support()
	main()
