import os
from pgp_reconstruction import config, project_dir
from pgp_reconstruction.cli.download_missing_data import download_missing_files

import os
import subprocess
import pickle
import sys

def saveProgressFile(progressValue, outputfolder):
	#Crete a file with the task's approximate progress. It is useful only for the web application. So the web application can use the same package version as in GitHub.
	
	progress_path = os.path.join(outputfolder, 'progress.txt')
	
	with open(progress_path, "w") as f:
		f.write(str(progressValue))
	
	return
	
	
def first_run_check(updateDB):

	if updateDB:
		#it can be called even if not first call
		firstRun = download_missing_files()
		if firstRun: return firstRun

	diamond_db = project_dir + config.get('generated', 'diamond_db')
	if not os.path.exists(diamond_db):
	
		#only tries to connect to internet to download files while diamond_db is not found
		firstRun = download_missing_files()
		if firstRun: return firstRun
	
		print("Running diamond for the first time, please wait while we build the internal database...")
		fasta_file = project_dir + config.get('generated', 'fasta_file')
		cmd = ['diamond', 'makedb', '--in', fasta_file, '-d', diamond_db[:-5]]
		cmdStr = ''
		for i in cmd: cmdStr += i + ' ' 
		
		try:
			exit_code = subprocess.call(cmd)
		except OSError:
			raise ValueError('Unable to run diamond with the command "' + cmdStr + '"\nMake sure diamond is installed and available in your PATH.')
		else:
			if exit_code != 0:
				raise ValueError('Failed to run diamond (wrong arguments).')
				
	return 0


def loadConstraints(constraintsFilePath, cobraModel, reframedModel):
	#receive a file with ID of reactions and metabolites. IDs might be from ChEBI, KEGG, BIGG, SEED or Rhea. 
	#finds IDs in universal model, and asign a score
	

	pickle_file_path = os.path.join(project_dir, 'data/generated', 'keggModules.pickle')
	with open(pickle_file_path, 'rb') as f:
		keggModules = pickle.load(f)
		
	pickle_file_path = os.path.join(project_dir, 'data/generated', 'biocycPathways.pickle')
	with open(pickle_file_path, 'rb') as f:
		biocycPathways = pickle.load(f)
		
	allPaths = set()
	for i in keggModules: allPaths.add(i.lower())
	for i in biocycPathways: allPaths.add(i.lower())
		
		
	constraintsDict = {'metabolites':{'soft':dict(),'hard':dict()}, 'reactions':{'soft':dict(),'hard':dict()}, 'pathways':{'soft':dict(),'hard':dict()}}
	with open(constraintsFilePath) as file:
		_ = next(file)  # Skip the first line
		for line in file:
			line = line.strip()
			if not line: continue
			lineSplit = line.split('\t')
			#error handling
			if len(lineSplit) != 4:
				sys.exit('ERROR: Constraint should have 4 columns, separated by tabs. For example: "M_na    Soft    1    Media".\nInput given: ' + str(line))
			if lineSplit[3] == 'taxonomy': continue
			
			itemId = lineSplit[0]
			constraintType = lineSplit[1].lower()
			itemScore = lineSplit[2]
			group = lineSplit[3].lower()
			
			
			if len(itemId) <= 2:
				sys.exit('ERROR: The first column should contain the IDs of reactions or metabolites. If reaction, it should start with "R_"; If metabolite, it should start with "M_". We accept IDs from BIGG, SEED, and ChEBI. For example: "M_na".\nInput given: ' + str(line))
			if itemId[:2].lower() == 'm_':
				if group != 'media' and group != 'product':
					sys.exit('ERROR: It should be informed whether the constraint is regarding consumption (using the keyword "Media") or production (using the keyword "Product"). For example: "M_na    Soft    1    Media".\nInput given: ' + str(line))
			if group == 'media' or group == 'product':
				if itemId[:2].lower() != 'm_':
					sys.exit('ERROR: the ID on the first column does not match the format of a metabolite ID, but in the last column you informed it is a consumed/produced metabolite. Metabolite IDs should start with "M_". Example of correct input: "M_61988_e    Soft    0.1    Media".\nInput given: ' + str(line))
			if itemId[:2].lower() == 'r_':
				if group != 'reaction':
					sys.exit('ERROR: the id on the first column seems to be of a reaction, but this is not the information given in the last column. Example of correct input: "R_R09640    Soft    3    Reaction".\nInput given: ' + str(line))
			if group[:2] == 'reaction':
				if itemId[:2].lower() != 'r_':
					sys.exit('ERROR: the ID on the first column does not match the format of a reaction ID, but in the last column you informed it is a reaction. Reactions IDs should start with "R_". Example of correct input: "R_R09640    Soft    3    Reaction".\nInput given: ' + str(line))
			if group == 'pathway':
				if itemId.lower() not in allPaths:
					sys.exit('ERROR: pathway ID was not recognized either as a MetaCyc pathway or as a KEGG module.\nInput given: ' + str(line))
			if constraintType != 'soft' and constraintType != 'hard':
				sys.exit('ERROR: In the second column, it should be informed if the constraint is "soft" (model will try to satisfy the condition) or "hard" (model will necessarely satisfy the condition). For example: "R_R09640    Soft    3    Reaction".\nInput given: ' + str(line))
			#checking if score is number
			try: _ = float(itemScore)
			except: sys.exit('ERROR: In the third column, it should be informed the reaction score. If metabolite, the score will be attributed to the reactions producing/consuming the metabolite. If hard constraint, only the score sign will be taken into consideration to decide if the reaction should be included or avoided. For example: "R_R09640    Soft    3    Reaction".\nInput given: ' + str(line))
			
			
			
			#store data
			if itemId[:2].lower() == 'm_':
				if constraintType == 'soft': 
					if itemId[2:] in constraintsDict['metabolites']['soft'] and group != constraintsDict['metabolites']['soft'][itemId[2:]][mediaOrProduct]:
						constraintsDict['metabolites']['soft'][itemId[2:]][mediaOrProduct] = 'media and product'
					else:
						constraintsDict['metabolites']['soft'][itemId[2:]] = {'score':float(itemScore), 'mediaOrProduct': group}
				elif constraintType == 'hard': 
					constraintsDict['metabolites']['hard'][itemId[2:]] = {'score':float(itemScore), 'mediaOrProduct': group}
			
			elif itemId[:2].lower() == 'r_':
				if constraintType == 'soft': constraintsDict['reactions']['soft'][itemId[2:]] = {'score':float(itemScore)}
				elif constraintType == 'hard': constraintsDict['reactions']['hard'][itemId[2:]] = {'score':float(itemScore)}
			
			elif group == 'pathway':
				if constraintType == 'soft': constraintsDict['pathways']['soft'][itemId.lower()] = {'score':float(itemScore)}
				elif constraintType == 'hard': constraintsDict['pathways']['hard'][itemId.lower()] = {'score':float(itemScore)}
	
	
	#find the exchange reactions producing or consuming the metabolite
	constraintsFromFile = {'soft':dict(),'hard':dict()}
	for constraintType in constraintsDict['metabolites']:
		for metId in constraintsDict['metabolites'][constraintType]:
			metInModel = None
			mediaOrProduct = constraintsDict['metabolites'][constraintType][metId]['mediaOrProduct']
			score = constraintsDict['metabolites'][constraintType][metId]['score']
			if metId in cobraModel.metabolites:
				metInModel = cobraModel.metabolites.get_by_id(metId)
			else:
				sucess = 0
				for met in cobraModel.metabolites:
					for db in ['bigg', 'chebi', 'kegg', 'seed']:
						if db not in met.annotation: continue
						if metId in met.annotation[db] or metId + '_e' in met.annotation[db]:
							sucess = 1
							break
					if sucess == 1: break
				if sucess == 1: metInModel = met
			
			if not metInModel: continue
			
			for rxn in metInModel.reactions:
				if len(rxn.metabolites) != 1: continue
				
				sucess = 0
				if mediaOrProduct == 'product' and ((rxn.lower_bound < 0 and rxn.metabolites[metInModel] > 0) or (rxn.upper_bound > 0 and rxn.metabolites[metInModel] < 0)): #producing
					constraintsFromFile[constraintType][rxn.id] = score
					sucess = 1
				
				if mediaOrProduct == 'media' and ((rxn.upper_bound > 0 and rxn.metabolites[metInModel] > 0) or (rxn.lower_bound < 0 and rxn.metabolites[metInModel] < 0)): #producing
					constraintsFromFile[constraintType][rxn.id] = score
					sucess = 1
			
				#if there are no constraints for both directions, makes the oposite direction imposible
				if sucess == 1 and constraintsDict['metabolites'][constraintType][metId]['mediaOrProduct'] != 'media and product':
				
					rxnReverse = None
					if 'forwardTemp' in rxn.id:
						rxnId = rxn.id.replace('forwardTemp', 'reverseTemp')
						if rxnId in cobraModel.reactions:
							rxnReverse = cobraModel.reactions.get_by_id(rxnId)
					if 'reverseTemp' in rxn.id:
						rxnId = rxn.id.replace('reverseTemp', 'forwardTemp')
						if rxnId in cobraModel.reactions:
							rxnReverse = cobraModel.reactions.get_by_id(rxnId)
					if rxnReverse:
						rxnReverse.lower_bound = 0
						rxnReverse.upper_bound = 0
						reframedId = 'R_'+rxnReverse.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
						reframedModel.reactions[reframedId].lb = 0
						reframedModel.reactions[reframedId].ub = 0
	
	
	#faz o mesmo para as reacoes
	for constraintType in constraintsDict['reactions']:
		for rxnId in constraintsDict['reactions'][constraintType]:
			
			score = constraintsDict['reactions'][constraintType][rxnId]['score']
			inModel = list()
			if rxnId in cobraModel.reactions:
				rxn = cobraModel.reactions.get_by_id(rxnId)
				inModel.append(rxn)
			else:
			
				for rxn in cobraModel.reactions:
					if rxn.id == rxnId + '_reverseTemp': inModel.append(rxn)
					elif rxn.id == rxnId + '_forwardTemp' : inModel.append(rxn)
					
				if not inModel:
			
					sucess = 0
					for db in ['rhea', 'bigg', 'kegg', 'seed', '']:
						for rxn in cobraModel.reactions:
							if db not in rxn.annotation: continue
							if rxnId in rxn.annotation[db] or rxnId.replace('_reverseTemp','').replace('_forwardTemp','') in rxn.annotation[db]:
								inModel.append(rxn)
						if inModel: break
			
			for rxn in inModel:
				constraintsFromFile[constraintType][rxn.id] = score
				
				
	#faz o mesmo para as pathways
	for constraintType in constraintsDict['pathways']:
		for pathwayIdLower in constraintsDict['pathways'][constraintType]:
			score = constraintsDict['pathways'][constraintType][pathwayIdLower]['score']
			#look for reactions translatable to kegg ids
			for pathwayId in keggModules:
				if pathwayId.lower() != pathwayIdLower: continue
				for eachProcess in keggModules[pathwayId]['RxnsInvolved']:
					for eachSet in eachProcess:
						for eachKeggId in eachSet:
							for rxn in cobraModel.reactions:
								if 'kegg' not in rxn.annotation: continue
								if eachKeggId in rxn.annotation['kegg']:
									if rxn.id in constraintsFromFile[constraintType]:
										#check if only entry has the same sign as the new entry
										if constraintsFromFile[constraintType][rxn.id] * score > 0:
											constraintsFromFile[constraintType][rxn.id] += 1
										else: pass
									else:
										constraintsFromFile[constraintType][rxn.id] = score
								
			#look for reactions translatable to metacyc ids
			for pathwayId in biocycPathways:
				if pathwayId.lower() != pathwayIdLower: continue
				for eachBiocycId in biocycPathways[pathwayId]['RxnsInvolved']:
					for rxn in cobraModel.reactions:
						if 'metacyc' not in rxn.annotation: continue
						if eachBiocycId in rxn.annotation['metacyc']:
							if rxn.id in constraintsFromFile[constraintType]:
								#check if only entry has the same sign as the new entry
								if constraintsFromFile[constraintType][rxn.id] * score > 0:
									constraintsFromFile[constraintType][rxn.id] += 1
								else: pass
							else:
								constraintsFromFile[constraintType][rxn.id] = score
								
	
	#manipulate hard constraints to remove same reaction with two directions from hard constraints
	hardInconsistent = dict()
	for cobraId in constraintsFromFile['hard']:
		cobraIdSimple = cobraId.replace('forwardTemp', '').replace('reverseTemp', '')
		if cobraIdSimple in hardInconsistent: hardInconsistent[cobraIdSimple] += 1
		else: hardInconsistent[cobraIdSimple] = 1
	for cobraIdSimple in hardInconsistent:
		if hardInconsistent[cobraIdSimple] > 1:
			del constraintsFromFile['hard'][cobraIdSimple+'forwardTemp']
			del constraintsFromFile['hard'][cobraIdSimple+'reverseTemp']
	
	
	return constraintsFromFile