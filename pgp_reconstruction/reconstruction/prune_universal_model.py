
import pickle
from pgp_reconstruction import project_dir
import cobra
from datetime import datetime
import os
import copy
import numpy as np
import re
from itertools import combinations


from pgp_reconstruction.reconstruction.makeEssentialRxns import makeEssentialGenesEssential
from pgp_reconstruction.cli.util import saveProgressFile
from pgp_reconstruction.reconstruction.scoring import findPathways

from reframed.solvers import solver_instance
from reframed.solvers.solver import VarType
from reframed.solvers.solution import Status

def includePathways(cobraModel, biocycPathways, rxnsPerModules):
		
	pickle_file_path = os.path.join(project_dir, 'data/generated', 'keggModules.pickle')
	with open(pickle_file_path, 'rb') as f:
		keggModules = pickle.load(f)
	
	dbToRxnsInModel = dict()
	for rxn in cobraModel.reactions:
		for db in ['metacyc', 'kegg']:
			if db not in rxn.annotation: continue
			for rxnId in rxn.annotation[db]:
				if rxnId not in dbToRxnsInModel:
					dbToRxnsInModel[rxnId] = set()
				dbToRxnsInModel[rxnId].add(rxn)
	
	minPathResult = findPathways(list(dbToRxnsInModel.keys()))
	
	#map biocyc/kegg rxns to rxn In model
	for pathway in minPathResult:
		if pathway not in minPathResult: continue
		
		for pathway in biocycPathways:
			for rxnId in biocycPathways[pathway]['RxnsInvolved']:
				if rxnId not in dbToRxnsInModel: continue
				for rxn in dbToRxnsInModel[rxnId]:
					if 'MetaCyc pathways' not in rxn.annotation: rxn.annotation['MetaCyc pathways'] = list()
					rxn.annotation['MetaCyc pathways'].append(biocycPathways[pathway]['name'])

		for pathway in rxnsPerModules:
			for rxnId in rxnsPerModules[pathway]:
				if rxnId not in dbToRxnsInModel: continue
				for rxn in dbToRxnsInModel[rxnId]:
					if 'KEGG pathways' not in rxn.annotation: rxn.annotation['KEGG pathways'] = list()
					rxn.annotation['KEGG pathways'].append(keggModules[pathway]['name'])	
	
	

def initSolver(reframedModel, positiveInObjective, constraintsFromFile=None, cobraModel=None, firstOptimization=0, notIncludeReframed=[], eps=1e-3, bigM=1000, min_growth=0.00001):

	#solver5 = initSolver(reframedModel, positiveInSolution|positiveDiffusion) # essentialRxnsInSpec = set()

	#multiple optimization problems are solved. To avoid repetition, a function is used to create the problem.


	solver = solver_instance(reframedModel)
	
	solver.neg_vars = []
	solver.pos_vars = []
	for rxnId in reframedModel.reactions: #reactions:
	
		if rxnId in notIncludeReframed: continue
	
		y_r, y_f = 'yr_' + rxnId, 'yf_' + rxnId
	
		if reframedModel.reactions[rxnId].lb < 0:
			solver.add_variable(y_r, 0, 1, vartype=VarType.BINARY, update=False)
			solver.neg_vars.append(y_r)
		if reframedModel.reactions[rxnId].ub > 0:
			solver.add_variable(y_f, 0, 1, vartype=VarType.BINARY, update=False)
			solver.pos_vars.append(y_f)
	solver.update()
	
	
	if firstOptimization:
		#force flux on rxns with genes, with positive score, and with flux when using FBA, to reduce solution space
	
		with cobraModel as model:
			
			boundaryIds = {rxnCobra.id for rxnCobra in model.boundary}
		
			for rxnCobra in model.reactions:
			
				if len(rxnCobra.compartments) == 2: continue
				if rxnCobra.id in boundaryIds: continue
				if rxnCobra.genes: continue
				if rxnCobra.lower_bound < 0: rxnCobra.lower_bound = -0.1
				if rxnCobra.upper_bound > 0: rxnCobra.upper_bound = 0.1
	
			fba = model.optimize()
		
		revelantRxnsWithFlux = list()
		for rxnCobra in cobraModel.reactions:
			if abs(fba.fluxes[rxnCobra.id]) > 0 and rxnCobra.genes: 
				if 'R_' + rxnCobra.id not in positiveInObjective: positiveInObjective.add('R_' + rxnCobra.id)
				revelantRxnsWithFlux.append(rxnCobra.id.replace('_forwardTemp','').replace('_reverseTemp',''))
		
		for rnxId in revelantRxnsWithFlux:
			constraintsEssentials = dict()
			if 'R_' + rnxId in reframedModel.reactions:
				y_r, y_f = 'yr_' + 'R_' + rnxId, 'yf_' + 'R_' + rnxId
				if y_r in solver.neg_vars: constraintsEssentials[y_r] = 1
				if y_f in solver.pos_vars: constraintsEssentials[y_f] = 1
			if 'R_' + rnxId + '_forwardTemp' in reframedModel.reactions:
				y_r, y_f = 'yr_' + 'R_' + rnxId + '_forwardTemp', 'yf_' + 'R_' + rnxId + '_forwardTemp'
				if y_r in solver.neg_vars: constraintsEssentials[y_r] = 1
				if y_f in solver.pos_vars: constraintsEssentials[y_f] = 1
			if 'R_' + rnxId + '_reverseTemp' in reframedModel.reactions:
				y_r, y_f = 'yr_' + 'R_' + rnxId + '_reverseTemp', 'yf_' + 'R_' + rnxId + '_reverseTemp'
				if y_r in solver.neg_vars: constraintsEssentials[y_r] = 1
				if y_f in solver.pos_vars: constraintsEssentials[y_f] = 1

			if constraintsEssentials:
				solver.add_constraint('forceRelevant_' + rnxId, constraintsEssentials, '=', 1, update=False)
		solver.update()
		
		
	if constraintsFromFile: hardConstraints = {'R_'+rxnId.replace('-','__45__').replace('.','__46__').replace('+','__43__') for rxnId in constraintsFromFile['hard']}
	else: hardConstraints = set()
	
	biomass = reframedModel.biomass_reaction
	#upper and lowe flux constraints:
	for rxnId in reframedModel.reactions:
	
		y_r, y_f = 'yr_' + rxnId, 'yf_' + rxnId
	
		if reframedModel.reactions[rxnId].lb < 0 or reframedModel.reactions[rxnId].ub > 0:

			if rxnId in positiveInObjective:
				#allow reverse and forward simultaneously if the score is negative
				
				if y_r in solver.neg_vars and '_reverseTemp' in y_r and y_f.replace('_reverseTemp','_forwardTemp') in solver.pos_vars: #when universal model is split
					solver.add_constraint('avoidTwoDirections_' + rxnId, {y_r: 1, 'yf_' + rxnId.replace('_reverseTemp','_forwardTemp'): 1}, '<', 1, update=False) 
				
				if y_r in solver.neg_vars and y_f in solver.pos_vars:
					solver.add_constraint('avoidTwoDirections_' + rxnId, {y_r: 1, y_f: 1}, '<', 1, update=False)
		
			#force flux on biomass
			if rxnId == reframedModel.biomass_reaction or rxnId in hardConstraints:
				if y_f in solver.pos_vars:
					solver.add_constraint('forceBiomass_' + rxnId, {y_f: 1}, '=', 1, update=False)
				elif y_r in solver.neg_vars:
					solver.add_constraint('forceBiomass_' + rxnId, {y_r: 1}, '=', 1, update=False)
				#solver.add_constraint('min_growth' + rxnId, {rxnId: 1}, '>', min_growth, update=False)
					
			#flux constraints
			if y_f in solver.pos_vars:
				solver.add_constraint('lb_pos_' + rxnId, {rxnId: 1, y_f: -eps}, '>', 0, update=False)
				solver.add_constraint('ub_pos_' + rxnId, {rxnId: 1, y_f: -bigM}, '<', 0, update=False)
			if y_r in solver.neg_vars:
				solver.add_constraint('lb_neg_' + rxnId, {rxnId: 1, y_r: bigM}, '>', 0, update=False)
				solver.add_constraint('ub_neg_' + rxnId, {rxnId: 1, y_r: eps}, '<', 0, update=False)		
	solver.update()
		
	return solver
	
def ensureFluxOnBiomass(activeRxnsreframedId3, possibleToRemove, reframedModel):

	positiveInObjective =  activeRxnsreframedId3
	
	solver5 = initSolver(reframedModel, positiveInObjective)

	objective = {}
	for rxnId in reframedModel.reactions:
	
		y_r, y_f = 'yr_' + rxnId, 'yf_' + rxnId
		
		if y_f in solver5.pos_vars:
			if rxnId in activeRxnsreframedId3: objective[y_f] = 10
			else: objective[y_f] = -1
			
		if y_r in solver5.neg_vars:
			if rxnId in activeRxnsreframedId3: objective[y_r] = 10
			else: objective[y_r] = -1
		
	
	solver5.update()
	solver5.set_objective(linear=objective, minimize=False)
	
	print('\nStart solving 5: ' + str(datetime.now()))
	solution5 = solver5.solve(emphasis=2, timelimit=600)
	
	if solution5.status != Status.OPTIMAL:
		print('\nStarting to solve 3.2: ' + str(datetime.now()))
		solution5 = solver5.solve(timelimit=3600)
	
	print('\nSolved 5: ' + str(datetime.now()))

	return solution5
	
	
def includeGenesRules(cobraModel, rheaIdToGene):

	##adiciona informacao de genes
	rxnsAndGenes = dict()
	for rxn in cobraModel.reactions:
		
		for eachDict in rheaIdToGene['bestMatchPerRead']:
			if 'R_' + rxn.id in eachDict['rxns']:
				if rxn.id not in rxnsAndGenes:
					rxnsAndGenes[rxn.id] = set()
				rxnsAndGenes[rxn.id].add(eachDict['gene'].replace(',','.').replace('(','').replace(')',''))
				
		if rxn.id not in rxnsAndGenes:
			bestMatch = dict()
			for eachDict in rheaIdToGene['notBestPerRead']:
				if 'R_' + rxn.id in eachDict['rxns']:
					if round(eachDict['score'],3) not in bestMatch:
						bestMatch[round(eachDict['score'],3)] = set()
					bestMatch[round(eachDict['score'],3)].add(eachDict['gene'])
			if bestMatch:
				rxnsAndGenes[rxn.id] = bestMatch[max(bestMatch)]
					

	for rxnId in rxnsAndGenes:
		ruleString = '('
		for uniprotId in rxnsAndGenes[rxnId]:
			if uniprotId.isdigit(): uniprotId = 'gene-' + uniprotId
			if ')' in uniprotId: uniprotId = uniprotId.replace(')', '-')
			if '(' in uniprotId: uniprotId = uniprotId.replace('(', '-')
			if ',' in uniprotId: uniprotId = uniprotId.replace(',', '-')
			ruleString += uniprotId + ' or '
		ruleString = ruleString[:-4] + ')'
		cobraModel.reactions.get_by_id(rxnId).gene_reaction_rule = ruleString
		
		
	genesToRxns = dict()
	for rxnId in rxnsAndGenes:
		for gene in rxnsAndGenes[rxnId]:
			if gene not in genesToRxns: genesToRxns[gene] = list()
			genesToRxns[gene].append(rxnId)
		
	essentialRxnsIds_Teoretical = set()
	for rxnId in rxnsAndGenes:
		essentialRxnsIds_Teoretical.add(rxnId)
	
	return genesToRxns, essentialRxnsIds_Teoretical
	
	

def saveReframedModel(reframedModel, solution3, passiveDiffusion, cobraModel, rheaIdToGene, rxnsScores, outputfile, spontaneousRheaRxns, biocycPathways, rxnsPerModules, outputfolder):
	#status = saveReframedModel(reframedModel, solution1, passiveDiffusion, cobraModel, rheaIdToGene, rxnsScores, outputfile, spontaneousRheaRxns)
	
	print('\nFinalizing model\n')


	#identify reactions with flux according to optimization
	activeRxns3 = set()
	activeRxnsreframedId3 = set()
	for rxnId in reframedModel.reactions:
		if ('yr_' + rxnId in solution3.values and solution3.values['yr_' + rxnId] > 1e-5) or ('yf_' + rxnId in solution3.values and solution3.values['yf_' + rxnId] > 1e-5) or (abs(solution3.values[rxnId]) > 1e-5):
			activeRxns3.add(rxnId[2:].replace('__45__','-').replace('__46__','.').replace('__43__','+').replace('_forwardTemp','').replace('_reverseTemp',''))
			activeRxnsreframedId3.add(rxnId)
				
	#check if solution allow for biomass production. Sometimes, because the solution is not optimal, there will be no biomass
	canGrow = 0
	with cobraModel as model:

		for rxn in model.reactions:
		
			rxnId = rxn.id.replace('_forwardTemp','').replace('_reverseTemp','')
			if rxnId not in activeRxns3:
				if rxn.upper_bound > 0: rxn.upper_bound = 0.01
				if rxn.lower_bound < 0: rxn.lower_bound = -0.01
			
		fba = model.optimize()
		
		if fba.status == 'optimal' and fba.objective_value != 0: canGrow = 1
		else: canGrow = 0
		
	if canGrow == 0:
	
		reframedModelTemp = copy.deepcopy(reframedModel)
		cobraModelTemp = copy.deepcopy(cobraModel)
		
		#reduce flux in rxns not from solution, to make reactions from solution more likely to be active.
		for rxn in cobraModelTemp.reactions:
		
			rxnId = rxn.id.replace('_forwardTemp','').replace('_reverseTemp','')
			if rxnId not in activeRxns3:
				if rxn.upper_bound > 0: rxn.upper_bound = 0.01
				if rxn.lower_bound < 0: rxn.lower_bound = -0.01
			
		fba = cobraModelTemp.optimize()
			
		fluxOnFba = set()
		for rxn in cobraModelTemp.reactions:
			rxnId = rxn.id.replace('_forwardTemp','').replace('_reverseTemp','')
			if fba.fluxes[rxn.id] != 0:
				fluxOnFba.add(rxnId)
				

		#remove as reacoes sem fluxo
		rxnsToRemove = set()
		possibleToRemove = set()
		for rxn in cobraModelTemp.reactions:
		
			rxnId = rxn.id.replace('_forwardTemp','').replace('_reverseTemp','')
			
			if rxnId in activeRxns3: pass
			elif rxnId in fluxOnFba: possibleToRemove.add('R_'+rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__'))
			else: rxnsToRemove.add(rxn)
		

		#deactivate reactions in rxnsToRemove
		for rxn in rxnsToRemove:
			
			rxn.lower_bound = 0
			rxn.upper_bound = 0
		
			idInreframed = 'R_'+rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
			reframedModelTemp.reactions[idInreframed].lb = 0
			reframedModelTemp.reactions[idInreframed].ub = 0
			
		#sometimes, there is no solution on solution5
		solution5 = ensureFluxOnBiomass(activeRxnsreframedId3, possibleToRemove, reframedModelTemp)
		
		activeRxns5 = set()
		if solution5.status == Status.OPTIMAL:		
			for rxnId in reframedModelTemp.reactions:
				if ('yr_' + rxnId in solution5.values and solution5.values['yr_' + rxnId] > 1e-5) or ('yf_' + rxnId in solution5.values and solution5.values['yf_' + rxnId] > 1e-5) or (abs(solution5.values[rxnId]) > 1e-5):
					if rxnId.startswith('R_'): rxnId = rxnId[2:]
					activeRxns5.add(rxnId.replace('__45__','-').replace('__46__','.').replace('__43__','+').replace('_forwardTemp','').replace('_reverseTemp',''))
					
		activeRxns3 = activeRxns3|activeRxns5
			
	saveProgressFile(96, outputfolder)
			
	##cria um novo modelo apenas com as reacoes com fluxo em solution5 e solution3. Eh mais rapido que deletar as reacoes que nao funcionam
	#find metabolites in reactions on solution5 and solution3
	metsToInclude = set()
	rxnsToInclude = set()
	for rxn in cobraModel.reactions:
		if rxn.id.replace('_forwardTemp','').replace('_reverseTemp','') in activeRxns3:
			rxnsToInclude.add(rxn.id)
			for met in rxn.metabolites:
				metsToInclude.add(met.id)
	
	
	#check		
	passiveDiffusionInC = {metId.replace('_e','_c') for metId in passiveDiffusion}
	passiveToInclude = passiveDiffusionInC.intersection(metsToInclude)
	passiveAlreadyIncluded = passiveDiffusionInC.intersection(metsToInclude)
	passiveToInclude = passiveToInclude-passiveAlreadyIncluded
	passiveDiffusionInE = {metId.replace('_c','_e') for metId in passiveToInclude}
	metsToInclude = metsToInclude|passiveDiffusionInE
	
	
	#find transport and exchange reactions
	if passiveDiffusionInE:
		for rxn in cobraModel.reactions:
			if rxn.id.replace('_forwardTemp','').replace('_reverseTemp','') in rxnsToInclude: continue
			for metId in passiveDiffusionInE:
				#looking for exchange
				if len(rxn.metabolites) == 1:
					sucesso = 0
					for met in rxn.metabolites: 
						if met.id == metId: sucesso = 1
					if sucesso == 1: rxnsToInclude.add(rxn.id)
				if len(rxn.compartment) == 2:
					sucesso = 0
					metsInRxn = {met.id for met in rxn.metabolites}
					if metId in metsInRxn:
						if metsInRxn.intersection(metsToInclude): sucesso = 1
					if sucesso == 1: rxnsToInclude.add(rxn.id)


	finalModel = cobra.Model('Model from ')
	
	#metabolitos
	for met in cobraModel.metabolites:
		if met.id not in metsToInclude: continue
		my_metabolite = cobra.Metabolite(met.id)
		my_metabolite.name = met.name
		my_metabolite.formula = met.formula
		my_metabolite.charge = met.charge
		my_metabolite.compartment = met.compartment
		my_metabolite.annotation = met.annotation
		_ = finalModel.add_metabolites([my_metabolite])
		
	#reacoes
	for rxn in cobraModel.reactions:
		if rxn.id not in rxnsToInclude: continue

		rxnIdToInclude = rxn.id.replace('_forwardTemp','').replace('_reverseTemp','')

		#if the reaction was already included, just change the bounds
		if rxnIdToInclude in finalModel.reactions:
			rxnInFinal = finalModel.reactions.get_by_id(rxnIdToInclude)
			if rxn.lower_bound < 0: rxnInFinal.lower_bound = -1000
			if rxn.upper_bound > 0: rxnInFinal.upper_bound = 1000
			continue
		
		my_reaction = cobra.Reaction(rxnIdToInclude)
		my_reaction.name = rxn.name
		my_reaction.annotation = rxn.annotation
		_ = finalModel.add_reactions([my_reaction])
		my_reaction.reaction = rxn.reaction
		my_reaction.objective_coefficient = rxn.objective_coefficient
		my_reaction.gene_reaction_rule = rxn.gene_reaction_rule
		my_reaction.lower_bound = rxn.lower_bound
		my_reaction.upper_bound = rxn.upper_bound

	#identify pathways and include information in each reaction
	includePathways(cobraModel, biocycPathways, rxnsPerModules)

	saveProgressFile(97, outputfolder)
	
	for sufix in [''] + list(range(1,100)):
		if os.path.isfile(os.path.splitext(os.path.basename(outputfile))[0] + str(sufix) + ".xml"): continue
		else: 
			cobra.io.write_sbml_model(finalModel, os.path.splitext(os.path.basename(outputfile))[0] + str(sufix) + ".xml")
			break
	return 1
		
		
def lowerScoreFromPromiscousEnzymes(cobraModel, biocycPathways, rxnsPerModules, rxnFromSolutionWithFlux, singleRxnsInGenes, rxnsInSol1, rheaIdToGene, rxnsScores):
	#some enzymes catalize multiple reactions. Maybe, just one of these reactions is important of the organism. This function make negative the score of reactions not present in almost complete pathways, or with not flux in pFBA. The score is only changed if such reaction is, in the ezyme-reaction association, always togther with a reaction associated as relevant.
	

	dbToRxnsInModel = dict()
	for rxn in cobraModel.reactions:
		if 'R_' + rxn.id not in rxnsInSol1: continue
		for db in ['metacyc', 'kegg']:
			if db not in rxn.annotation: continue
			for rxnId in rxn.annotation[db]:
				if rxnId not in dbToRxnsInModel:
					dbToRxnsInModel[rxnId] = set()
				dbToRxnsInModel[rxnId].add(rxn.id)
	
	metacycPlusKeggInModel = set(dbToRxnsInModel.keys())
	minPathResult = findPathways(metacycPlusKeggInModel)
	

	#check completness of all big biocyc pathway
	rxnsInPathways = set()
	for pathwayBiocyc in minPathResult:
		if pathwayBiocyc not in biocycPathways: continue
		synInModel = biocycPathways[pathwayBiocyc]['RxnsInvolved'].intersection(metacycPlusKeggInModel)
		if len(synInModel) / len(biocycPathways[pathwayBiocyc]['RxnsInvolved']) >= 0.5:
			for synId in synInModel: 
				for rxnId in dbToRxnsInModel[synId]:
					rxnsInPathways.add(rxnId)
				
	#check completness of all big kegg pathway
	for pathwayKegg in minPathResult:
		if pathwayKegg not in rxnsPerModules: continue
		synInModel = rxnsPerModules[pathwayKegg].intersection(metacycPlusKeggInModel)
		if len(synInModel) / len(rxnsPerModules[pathwayKegg]) >= 0.5:
			for synId in synInModel: 
				for rxnId in dbToRxnsInModel[synId]:
					rxnsInPathways.add(rxnId)

	
	toKeep = rxnsInPathways|rxnFromSolutionWithFlux|singleRxnsInGenes.intersection(rxnsInSol1)
	
	allRxnsFromGenes = set()
	for eachDict in rheaIdToGene['bestMatchPerRead']:
			for rxnId in eachDict['rxns']: allRxnsFromGenes.add(rxnId)
	

	for rxnId in allRxnsFromGenes - toKeep:
		sucess = 1
		for eachDict in rheaIdToGene['bestMatchPerRead']:
			if rxnId not in eachDict['rxns']: continue
			noPriority = eachDict['rxns'] - toKeep
			if not noPriority or len(noPriority) == eachDict['rxns']: 
				sucess = 0
				break
			
		#the score should only be changed if the reaction is together, in every gene, with reactions present in toKeep
		if sucess == 1 and rxnsScores[rxnId] > 0: 
			rxnsScores[rxnId] = -0.1
			
	return toKeep
	
def keepOneRheaSyn(rheaWithSameSyn, toKeep, rxnsScores, rheaIdToGene):
	#because the way the universal model was created, some rhea reactions can be synonims of other rhea rxns.Try to keep only one with a positive score, if possible.
	
	positiveScore = {rxnId for rxnId in rxnsScores if rxnsScores[rxnId] > 0}
	
	for rxnsSet in rheaWithSameSyn:
		toKeepInSet = rxnsSet.intersection(toKeep)
		if toKeepInSet and len(toKeepInSet) != len(rxnsSet):
			notKeepInSet = rxnsSet - toKeepInSet
			notKeepPositive = notKeepInSet.intersection(positiveScore)
			if not notKeepPositive: continue
		
			for rxnId in notKeepPositive:
				sucess = 1
				for eachDict in rheaIdToGene['bestMatchPerRead']:
					if rxnId not in eachDict['rxns']: continue
					if not toKeepInSet.intersection(eachDict['rxns']):
						sucess = 0
						break
						
				if sucess == 1 and rxnsScores[rxnId] > 0: 
					rxnsScores[rxnId] = -0.1
					
def findBiologComposition():
	
	import re

	# Initialize an empty dictionary to store CHEBI, KEGG, and MetaCyc IDs against Media IDs
	media_dict = {}

	# Open the file and read it line by line
	fileAndSource = {'Biolog_PM1.txt': 'carbon_1', 'Biolog_PM2.txt': 'carbon_2', 'Biolog_PM3.txt': 'nitrogen', 'Biolog_PM4.txt': 'phosphorus_and_sulfur'}

	for mediaFile in fileAndSource:
		source = fileAndSource[mediaFile]
		media_dict[source] = {}
		current_media_id = None  # Variable to hold the current media ID
		
		midiaPath = os.path.join(project_dir, 'data/generated', mediaFile)
		with open(midiaPath, 'r') as f:
			for line in f:
				# Strip leading and trailing whitespace
				line = line.strip()
				
				# Detect if line indicates a new media ID
				if line.startswith('-'):
					current_media_id = line.lstrip('-').strip()
					media_dict[source][current_media_id] = {'chebi': [], 'kegg': [], 'metacyc': []}
				
				# Search for lines containing "CHEBI:"
				if line.lower().startswith('chebi:'):
					ids = re.findall(r'\b\d+\b', line)
					media_dict[source][current_media_id]['chebi'].extend(ids)
				
				# Search for lines containing "KEGG Compound:"
				if line.lower().startswith('kegg compound:') or line.lower().startswith('kegg:'):
					ids = re.findall(r'[A-Za-z]\d+', line.split(':',1)[-1])
					media_dict[source][current_media_id]['kegg'].extend(ids)
				
				# Search for lines containing "BioCyc:"
				if line.lower().startswith('biocyc:'):
					for meta in line.strip().split(':',1)[-1].split(' '):
						if 'META:' in meta: media_dict[source][current_media_id]['metacyc'].append(meta.replace('META:',''))
	
			
	media_dict['phosphorus'] = dict()
	for letter in ['A','B','C','D','E']:
		for number in range(1,13):
			media_dict['phosphorus'][letter+str(number)] = dict(media_dict['phosphorus_and_sulfur'][letter+str(number)])

	media_dict['sulfur'] = dict()
	for letter in ['F','G','H']:
		for number in range(1,13):
			media_dict['sulfur'][letter+str(number)] = dict(media_dict['phosphorus_and_sulfur'][letter+str(number)])

	del media_dict['phosphorus_and_sulfur']				
						
	return media_dict
				

def mediaToExchangeRxns(cobraModel, media_dict):
	
	dbToRxns = dict()
	negativeMediaRxns = {'only Cl':set(),'only Na':set(),'only K':set(),'only Mg':set(),'only Fe':set()}
	otherMediaRxns = {'only O':set(), 'only Light':set(), 'nitrogen gas':set(), 'carbon dioxide':set()}
	for rxn in cobraModel.boundary:
		if rxn.compartments != {'e'}: continue
		if len(rxn.metabolites) != 1: continue
		for met in rxn.metabolites: break
			
		#look for reactions producing, instead of consuming
		if (rxn.metabolites[met] < 0 and rxn.lower_bound < 0) or (rxn.metabolites[met] > 0 and rxn.upper_bound > 0):
		
			elementsSum = 0
			if 'C' in met.elements: elementsSum += 1
			if 'N' in met.elements: elementsSum += 1
			if 'P' in met.elements: elementsSum += 1
			if 'S' in met.elements: elementsSum += 1
			if elementsSum > 1: continue
			
			if 'Cl' in met.elements and len(met.elements) == 1: negativeMediaRxns['only Cl'].add(rxn.id)
			if 'Mg' in met.elements and len(met.elements) == 1: negativeMediaRxns['only Mg'].add(rxn.id)
			if 'Fe' in met.elements and len(met.elements) == 1: negativeMediaRxns['only Fe'].add(rxn.id)
			if 'Na' in met.elements and len(met.elements) == 1: negativeMediaRxns['only Na'].add(rxn.id)
			if 'K' in met.elements and len(met.elements) == 1: negativeMediaRxns['only K'].add(rxn.id)
			if 'O' in met.elements and len(met.elements) == 1: otherMediaRxns['only O'].add(rxn.id)
			
			for dbSym in ['metacyc', 'chebi', 'kegg']:
				if dbSym in met.annotation:
					for synId in met.annotation[dbSym]:	
						if synId not in dbToRxns: dbToRxns[synId] = set()
						dbToRxns[synId].add(rxn.id)
						
				
	mediaToRxns = dict()
	for source in media_dict:
		mediaToRxns[source] = dict()
		for media_id in media_dict[source]:
			mediaToRxns[source][media_id] = set()	
			for db in media_dict[source][media_id]:
				for synID in media_dict[source][media_id][db]:
					if synID not in dbToRxns: continue
					for rxnId in dbToRxns[synID]: mediaToRxns[source][media_id].add(rxnId)
			
	try: otherMediaRxns['only Light'] = dbToRxns['Light']
	except: pass
	
	try:
		for rxnId in dbToRxns['NITROGEN-MOLECULE']:
			if '17997' in rxnId: otherMediaRxns['nitrogen gas'].append(rxnId)
	except: pass

	try:
		for rxnId in dbToRxns['NITROGEN-MOLECULE']:
			if '17997' in rxnId: otherMediaRxns['carbon dioxide'].append(rxnId)
	except: pass
					
	return mediaToRxns, negativeMediaRxns, otherMediaRxns
							
			
	
def findNegativeBiologMedia(cobraModel):
	
	media_dict = findBiologComposition()
	mediaToRxns, negativeMediaRxns, otherMediaRxns = mediaToExchangeRxns(cobraModel, media_dict)
	
	sourcesInBIolog = {'carbon':{'α- D-Glucose':'C9','D-Fructose':'C7','Sucrose':'D11','Lactose':'D9','D-Mannose':'A11','D-Mannitol':'B11','D-Sorbitol':'B2','Acetic acid':'C8','Pyruvic acid':'H8','Glycerol':'B3'},
	'nitrogen':{'Ammonia':'A2',	'Nitrate':'A4'},
	'phosphorus':{'Phosphate':'A2','Pyrophosphate':'A3'},
	'sulfur':{'Sulfate':'F2', 'Thiosulfate':'F3', 'Tetrathionate':'F4'}}
	
	sourcesToRxns = dict()
	for source in sourcesInBIolog:
		sourcesToRxns[source] = set()
		for moleculeName in sourcesInBIolog[source]:
			mediaId = sourcesInBIolog[source][moleculeName]
			if source == 'carbon': source1 = 'carbon_1'
			else: source1 = source
			for rxnId in mediaToRxns[source1][mediaId]:
				sourcesToRxns[source].add(rxnId)
	
	for rxn in otherMediaRxns['carbon dioxide']: sourcesToRxns['carbon'].add(rxn)
	del otherMediaRxns['carbon dioxide']

	for rxn in otherMediaRxns['nitrogen gas']: sourcesToRxns['nitrogen'].add(rxn)
	del otherMediaRxns['nitrogen gas']
				
	return sourcesToRxns, negativeMediaRxns, otherMediaRxns


def generate_combinations(toVerifyRxns, lenList):

	indexToRxn = dict()
	count = 0
	for rxn in toVerifyRxns:
		count += 1
		indexToRxn[count] = rxn

	numbers = indexToRxn.keys()
	combinationsRxns = []
	
	#for lenList in reversed(range(1, len(toVerifyRxns) + 1)):
	for combo in combinations(numbers, lenList):
		combinationsRxns.append(list(combo))
	
	return indexToRxn, combinationsRxns

def prune_model(reframedModel, cobraModel, rxnsScores, constraintsFromFile, rheaIdToGene, bestMatchPerRead,
				taxoOfTarget, rheaWithSameSyn, singleRxnsInGenes, gapfillMinimumMedia, outputfolder, outputfile,
				min_growth=0.1, eps=1e-3, bigM=1e3, default_score=-6.0, solver=None, minimumMedia=False):


	pickle_file_path = os.path.join(project_dir, 'data/generated', 'spontaneousRedundant.pickle')
	with open(pickle_file_path, 'rb') as f:
		spontaneousRheaRxns = pickle.load(f)
		
	pickle_file_path = os.path.join(project_dir, 'data/generated', 'rheaRxnsAndTax.pickle')
	with open(pickle_file_path, 'rb') as f:
		rheaRxnsAndTax = pickle.load(f)
		
	pickle_file_path = os.path.join(project_dir, 'data/generated', 'chebiMets.pickle')
	with open(pickle_file_path, 'rb') as f:
		chebiMets = pickle.load(f)
		
	pickle_file_path = os.path.join(project_dir, 'data/generated', 'biocycPathways.pickle')
	with open(pickle_file_path, 'rb') as f:
		biocycPathways = pickle.load(f)
		
	pickle_file_path = os.path.join(project_dir, 'data/generated', 'rxnsPerModules.pickle')
	with open(pickle_file_path, 'rb') as f:
		rxnsPerModules = pickle.load(f)
		
	saveProgressFile(41, outputfolder)

	base_score = min(rxnsScores.values())
	default_score = base_score -2

	#find metabolites with biological function:
	withWithBiologicalRole = set()
	fundamentalRelationship = set()
	metabolitesRelationship = set()
	for chebiId in chebiMets:
		if 'relationship' in chebiMets[chebiId] and chebiMets[chebiId]['relationship']:
			if 'fundamental metabolite' in chebiMets[chebiId]['relationship']: fundamentalRelationship.add(chebiId)
			if 'metabolite' in chebiMets[chebiId]['relationship']: metabolitesRelationship.add(chebiId)
			withWithBiologicalRole.add(chebiId)

	genesToRxns, essentialRxnsIds_Teoretical = includeGenesRules(cobraModel, rheaIdToGene)	
		
	rheaRxnsAndSyn = dict()
	biomass = reframedModel.biomass_reaction
	for rxn in cobraModel.reactions:
		idInreframed = 'R_' + rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
		rheaRxnsAndSyn[idInreframed] = {idInreframed[2:]}
		if 'rhea' in rxn.annotation:
			for rheaRxn in rxn.annotation['rhea']:
				rheaRxnsAndSyn[idInreframed].add(rheaRxn)

	rxnsScores[biomass] = 50
	
	saveProgressFile(47, outputfolder)

	#facilitates the inclusion of reactions of transport by passive diffusion 
	passiveDiffusion = set()
	for rxnCobra in cobraModel.reactions:
	
		idInreframed = 'R_' + rxnCobra.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
			
		if idInreframed in rxnsScores: continue
			
		if 'biomass metabolite demand' in rxnCobra.name:
			rxnsScores[idInreframed] = -0.1
			
		elif ' sink' in rxnCobra.name:
			rxnsScores[idInreframed] = -8
	
		elif rxnCobra.id.startswith('EX_'): 
			if '--' in rxnCobra.reaction: left_side, right_side = rxnCobra.reaction.split('--')
			if '=' in rxnCobra.reaction: left_side, right_side = rxnCobra.reaction.split('=')
			
			left_matches = [metId.replace('_e','') for metId in re.findall(r'\d+_e', left_side)]
			right_matches = [metId.replace('_e','') for metId in re.findall(r'\d+_e', right_side)]
			
			if (left_matches and not right_matches) and rxnCobra.lower_bound < 0:
				rxnsScores[idInreframed] = -0.1
			elif (right_matches and not left_matches) and rxnCobra.upper_bound > 0:
				rxnsScores[idInreframed] = -0.1
			else:
				rxnsScores[idInreframed] = -1
				
			#find metabolites that potentially can be passively transported 
			for met in rxnCobra.metabolites:
				if 'R' not in met.elements and met.formula_weight < 50 and met.charge == 0: passiveDiffusion.add(met.id)
		
		elif 'acidDissociation_' in idInreframed: rxnsScores[idInreframed] = -0.01


	min_growth=1
	minScore = min(rxnsScores.values())
	minGraterThanZero, minLowerThanZero = float(minScore), float(minScore)
	if minLowerThanZero > 0: minLowerThanZero = 0
	if minGraterThanZero < 0: minGraterThanZero = 0
	
	
	for idInreframed in reframedModel.reactions:
		
		rxnIdOriginal = idInreframed.replace('__45__','-').replace('__46__','.').replace('__43__','+')
		if rxnIdOriginal.startswith('R_'): rxnIdOriginal = rxnIdOriginal[2:]
		rxnIdOriginal = rxnIdOriginal.split('_')[0]
		
		if rxnIdOriginal in spontaneousRheaRxns and (idInreframed not in rxnsScores): 
			rxnsScores[idInreframed] = minScore
		
		if idInreframed in rxnsScores: continue
		
		rxnCobra = cobraModel.reactions.get_by_id(idInreframed[2:].replace('__45__','-').replace('__46__','.').replace('__43__','+'))
		
		#define score for transport reactions
		if len(rxnCobra.compartment) == 2:
			#makes easier for mets to leave the cell than to enter the cell
			if '--' in rxnCobra.reaction: left_side, right_side = rxnCobra.reaction.split('--')
			if '=' in rxnCobra.reaction: left_side, right_side = rxnCobra.reaction.split('=')
			
			left_matches = [metId.replace('_e','') for metId in re.findall(r'\d+_e', left_side)]
			right_matches = [metId.replace('_e','') for metId in re.findall(r'\d+_e', right_side)]
			
			if (left_matches and not right_matches) and rxnCobra.lower_bound < 0:
				rxnsScores[idInreframed] = -0.1
			elif (right_matches and not left_matches) and rxnCobra.upper_bound > 0:
				rxnsScores[idInreframed] = -0.1
			else:			
				if len(rxnCobra.metabolites) > 2: rxnsScores[idInreframed] = base_score -0.5 #give preference to more complex transport reactions
				else: rxnsScores[idInreframed] = base_score -1
			continue
	
		#gives default score to reactions
		rxnsScores[idInreframed] = default_score
		
	for idInreframed in reframedModel.reactions:
		#if it is not a Rhea reaction, but a biocyc translations, decrease score
		if not -0.01 <= rxnsScores[idInreframed] < 0: continue # to avoid flipping signs
		if 'metacyc' in rxnCobra.annotation: rxnsScores[idInreframed] += 0.01
		if 'kegg' in rxnCobra.annotation: rxnsScores[idInreframed] += 0.01
		

	if taxoOfTarget:
		#check taxonomical level of reactions
		taxoOfTarget.reverse()
		rxnsPerTax = dict()
		for tax in taxoOfTarget: rxnsPerTax[tax] = set()
		
		for idInreframed in reframedModel.reactions:
			for eachTxLevel in taxoOfTarget:
				sucess = 0
				for rheaSyn in rheaRxnsAndSyn[idInreframed]:
					if rheaSyn in rheaRxnsAndTax:
						if eachTxLevel.lower() in rheaRxnsAndTax[rheaSyn]: sucess = 1
				if sucess == 1:
					rxnsPerTax[eachTxLevel].add(idInreframed)
		
		
		#increse score of reactions identified in organisms sharing the same taxanomical lavel as the target organism.
		scoreToIncrese = list(np.arange(0.01, 0.5+0.01, 0.5/len(taxoOfTarget)))
		scoreToIncrese.reverse()
		idChanged = set()
		for eachTxLevel in taxoOfTarget:
			for rxnId in rxnsPerTax[eachTxLevel]:
				if rxnsScores[rxnId] > base_score: continue
				if rxnId in idChanged: continue
				elif rxnsScores[rxnId] + scoreToIncrese[taxoOfTarget.index(eachTxLevel)] > base_score -0.1: rxnsScores[rxnId] += base_score -0.1
				else: rxnsScores[rxnId] += scoreToIncrese[taxoOfTarget.index(eachTxLevel)]
	
	
	#reactions using mets with biological relevance, receive a boost in its score.
	for rxn in cobraModel.reactions:
		idInreframed = 'R_'+rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
		rxnWithRelevantMet = 0
		for met in rxn.metabolites:
			if 'chebi' not in met.annotation: continue
			for chebiId in met.annotation['chebi']:
				if chebiId in withWithBiologicalRole: rxnWithRelevantMet = 1
		if rxnWithRelevantMet == 1 and not -0.01 <= rxnsScores[idInreframed] < 0:
			rxnsScores[idInreframed] += 0.01
				
	
	positiveInObjective = set()
	for rxnId in reframedModel.reactions:
		if rxnId in rxnsScores and rxnsScores[rxnId] >= 0: positiveInObjective.add(rxnId)


	solver1 = initSolver(reframedModel, positiveInObjective, constraintsFromFile, cobraModel=cobraModel, firstOptimization=1)
	

	#creates objective function
	objective = {}
	for rxnId in reframedModel.reactions:
		
		y_r, y_f = 'yr_' + rxnId, 'yf_' + rxnId

		sucess = 0
		if y_f in solver1.pos_vars:
			if rxnsScores[rxnId] < 0: objective[y_f] = round(rxnsScores[rxnId],3)*3
			elif rxnsScores[rxnId] == 0: objective[y_f] = -0.0001
			else: objective[y_f] = round(rxnsScores[rxnId],3)

		if y_r in solver1.neg_vars:
			if rxnsScores[rxnId] < 0: objective[y_r] = round(rxnsScores[rxnId],3)*3
			elif rxnsScores[rxnId] == 0: objective[y_r] = -0.0001
			else: objective[y_r] = round(rxnsScores[rxnId],3)
	
	
	solver1.set_objective(linear=objective, minimize=False)
	print('\nStarting to solve 1: ' + str(datetime.now()))
	solution1 = solver1.solve(emphasis=1, timelimit=600)
	#solution1 = solver1.solve(emphasis=1)
	
	saveProgressFile(60, outputfolder)
	
	if solution1.status == Status.OPTIMAL: 
		print('Solved 1: ' + str(datetime.now()) + '\n')
	elif solution1.status == Status.UNKNOWN:
		errorMessage = 'It was not possible to solve the first optimization problem within the maximum time limit. The genome of the organism may be contaminated, incomplete or have the wrong taxonomical assigment. Alternatively, try including additional constraints (like providing information on the culture medium) or increasing the maximum processing time.'
		saveProgressFile("Failed to build model. " + errorMessage, outputfolder)
		raise Exception(errorMessage)
	else:
		errorMessage = 'It was not possible to find a solution to your problem. If you are using aditional constraints, try to relex them.'
		saveProgressFile("Failed to build model. " + errorMessage, outputfolder)
		raise Exception(errorMessage)


	#adjust scores based on first optimization	
	
	#makeEssentialGenesEssential specific for target organism
	deleted, rxnFromSolutionWithFlux = makeEssentialGenesEssential(solution1, bestMatchPerRead, rheaIdToGene, cobraModel, reframedModel, singleRxnsInGenes)
	

	cobraModel.remove_reactions(deleted)
	toDeletereframed = ['R_' + rxnId.replace('-','__45__').replace('.','__46__').replace('+','__43__') for rxnId in deleted]
	reframedModel.remove_reactions(toDeletereframed)
	
	saveProgressFile(74, outputfolder)


	#searching for the pathways present in solution1
	rxnsInSol1 = set()
	rxnsCobraInSol1 = set()
	sink = list()
	for rxnId in reframedModel.reactions:
		if ('yr_' + rxnId in solution1.values and solution1.values['yr_' + rxnId] > 1e-5) or ('yf_' + rxnId in solution1.values and solution1.values['yf_' + rxnId] > 1e-5) or (abs(solution1.values[rxnId]) > 1e-5):
			rxnsInSol1.add(rxnId)
			if rxnId.startswith('R_SK_'): sink.append(rxnId)
			rxnsCobraInSol1.add(rxnId[2:].replace('__45__','-').replace('__46__','.').replace('__43__','+'))
			if '_forwardTemp' in rxnId: rxnsInSol1.add(rxnId.replace('_forwardTemp','_reverseTemp'))
			if '_reverseTemp' in rxnId: rxnsInSol1.add(rxnId.replace('_reverseTemp','_forwardTemp'))
		

	toKeep = lowerScoreFromPromiscousEnzymes(cobraModel, biocycPathways, rxnsPerModules, rxnFromSolutionWithFlux, singleRxnsInGenes, rxnsInSol1, rheaIdToGene, rxnsScores)
	keepOneRheaSyn(rheaWithSameSyn, toKeep, rxnsScores, rheaIdToGene)
		
		
	#gapfill for minimum media
	if gapfillMinimumMedia:
	
		print('\nLine 916: ' + str(datetime.now()))
		
		sourcesToRxns, negativeMediaRxns, otherMediaRxns = findNegativeBiologMedia(cobraModel)
		
		simplifiedlMedia = set()
		for mediaType in [sourcesToRxns, negativeMediaRxns, otherMediaRxns]:
			for media in mediaType:
				for rxnId in mediaType[media]: simplifiedlMedia.add(rxnId)
		
		rxnsScoresMinimumMedia = dict(rxnsScores)
		
		allIntake = set()
		for rxn in cobraModel.boundary:
			for met in rxn.metabolites: break
				
			reframedRxnId = 'R_'+rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
				
			#look for reactions producing, instead of consuming
			if rxn in cobraModel.boundary:
				if (rxn.metabolites[met] < 0 and rxn.lower_bound < 0) or (rxn.metabolites[met] > 0 and rxn.upper_bound > 0):
					
					if rxn.id in simplifiedlMedia: 
						rxnsScoresMinimumMedia[reframedRxnId] = 0.1
					else: 
						rxnsScoresMinimumMedia[reframedRxnId] = -1
						allIntake.add(rxn)
					
		#save the boundaries
		allIntake_bounds = dict()
		for rxn in allIntake:
			allIntake_bounds[rxn] = {'lb':rxn.lower_bound ,'ub':rxn.upper_bound}
			
		#find essential exchangeRxns
		allIntake_toKeep = set()
		for rxn in allIntake:

			rxn.lower_bound, rxn.upper_bound = 0, 0
	
			fba_slim = cobraModel.slim_optimize()
			
			rxn.lower_bound, rxn.upper_bound = allIntake_bounds[rxn]['lb'], allIntake_bounds[rxn]['ub']
			if fba_slim == 0: allIntake_toKeep.add(rxn)
			
		
		for rxn in allIntake - allIntake_toKeep:
			for met in rxn.metabolites: break
			
			elementsSum = 0
			if 'C' in met.elements: elementsSum += 1
			if 'N' in met.elements: elementsSum += 1
			if 'P' in met.elements: elementsSum += 1
			if 'S' in met.elements: elementsSum += 1
			
			if elementsSum > 1: 
				rxn.lower_bound, rxn.upper_bound = 0, 0
			
		
		fba_slim = cobraModel.slim_optimize()
		if abs(fba_slim) > 1e-4:
		
			fba = cobraModel.optimize()
			#remove exchange rxns without flux on FBA
			allIntake_toRemove = set()
			for rxn in allIntake:
				if fba.fluxes[rxn.id] == 0: allIntake_toRemove.add(rxn)
			
			
			#verify essential exchange rxns
			toVerifyRxns = allIntake - allIntake_toRemove
			for rxn in toVerifyRxns:

				lb, ub = rxn.lower_bound, rxn.upper_bound
				rxn.lower_bound, rxn.upper_bound = 0, 0
		
				fba_slim = cobraModel.slim_optimize()
				
				rxn.lower_bound, rxn.upper_bound = lb, ub
				if fba_slim == 0: allIntake_toKeep.add(rxn)
				

			#simulate the removal of the rxns in toVerifyRxns. The more that can be removed, the better.
			toVerifyRxns = allIntake - (allIntake_toRemove|allIntake_toKeep)
			
			for rxn in allIntake_toRemove:
				rxn.lower_bound, rxn.upper_bound = 0, 0
	
			sucess = 0
			for lenList in reversed(range(1, len(toVerifyRxns) + 1)):
			
				indexToRxn, combinationsRxns = generate_combinations(toVerifyRxns, lenList)
				
				for combination in combinationsRxns:
			
					for rxnIndex in combination:
						indexToRxn[rxnIndex].lower_bound = 0
						indexToRxn[rxnIndex].upper_bound = 0
			
					fba_slim = cobraModel.slim_optimize()

					for rxnIndex in combination:
						indexToRxn[rxnIndex].lower_bound = allIntake_bounds[indexToRxn[rxnIndex]]['lb']
						indexToRxn[rxnIndex].upper_bound = allIntake_bounds[indexToRxn[rxnIndex]]['ub']
						
					if abs(fba_slim) > 1e-4: 
						sucess = 1
						break
				if sucess == 1: break
			

			for rxnIndex in combination:
				allIntake_toRemove.add(indexToRxn[rxnIndex])
				

			for rxn in allIntake_bounds:
				rxn.lower_bound = allIntake_bounds[rxn]['lb']
				rxn.upper_bound = allIntake_bounds[rxn]['ub']
			
			
			notIncludeCobra = {rxn.id for rxn in allIntake_toRemove}
			notIncludeReframed = {'R_'+ rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__') for rxn in allIntake_toRemove}
			
			with cobraModel as model:
			
				for rxnCobra in model.reactions:
					if rxnCobra.id in notIncludeCobra:
						if rxnCobra.lower_bound < 0: rxnCobra.lower_bound = 0
						if rxnCobra.upper_bound > 0: rxnCobra.upper_bound = 0
		
				fba = model.optimize()
				pfba_solution = cobra.flux_analysis.pfba(model)
			

			revelantRxnsWithFlux = set()
			allIntakeWithFlux = set()
			exchangeRxns = list()
			for rxn in cobraModel.reactions:
				if abs(pfba_solution.fluxes[rxn.id]) > 0: 
					if rxn.id in simplifiedlMedia:
						allIntakeWithFlux.add(rxn.id)
					else:
						if 'R_' + rxn.id not in positiveInObjective: positiveInObjective.add('R_' + rxn.id)
						revelantRxnsWithFlux.add('R_'+ rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__'))
						
			
			#reduce rxns in sourcesToRxns. there are multiple carbon sources. Only one should remain. Give priority to the ones from solution 1, and those with flux on FBA
			for mediaType in [sourcesToRxns, negativeMediaRxns, otherMediaRxns]:
				for media in mediaType:
					inter = mediaType[media].intersection(rxnsCobraInSol1)
					if not inter: inter = mediaType[media]
					inter2 = inter.intersection(allIntakeWithFlux)
					if not inter2: inter2 = inter
					if inter2 and inter2 != mediaType[media]:
						for rxnCobraId in mediaType[media] - inter2:
							reframedRxnId = 'R_'+rxnCobraId.replace('-','__45__').replace('.','__46__').replace('+','__43__')
							rxnsScoresMinimumMedia[reframedRxnId] = -0.5
						mediaType[media] = inter2
					
					
			print('\nLine 1020: ' + str(datetime.now()))
			
			#find all positive, to avoid double direction on these reactions
			positiveInObjective = set()
			for rxnId in reframedModel.reactions:
				if rxnId in rxnsScoresMinimumMedia and rxnsScoresMinimumMedia[rxnId] >= 0: positiveInObjective.add(rxnId)
			
			solver2 = initSolver(reframedModel, positiveInObjective, cobraModel=cobraModel, notIncludeReframed=notIncludeReframed)
			
			#creates objective function
			maximum = max(rxnsScoresMinimumMedia.values())
			objective = {}
			#maximum = max(rxnsScoresMinimumMedia.values())
			for rxnId in reframedModel.reactions:
				
				y_r, y_f = 'yr_' + rxnId, 'yf_' + rxnId

				if y_f in solver2.pos_vars:
					if rxnsScoresMinimumMedia[rxnId] < 0:
						#if rxnId in rxnsInSol1: objective[y_f] = round(rxnsScoresMinimumMedia[rxnId],3)/2
						#else: objective[y_f] = round(rxnsScoresMinimumMedia[rxnId],3)*2
						if rxnId in revelantRxnsWithFlux and rxnId in rxnsInSol1: objective[y_f] = -1
						elif rxnId in revelantRxnsWithFlux: objective[y_f] = -1.1
						elif rxnId in rxnsInSol1: objective[y_f] = -1.2
						else: objective[y_f] = -1.5
						
					elif rxnsScoresMinimumMedia[rxnId] == 0: objective[y_f] = -0.0001
					else:
						#objective[y_f] = rxnsScoresMinimumMedia[rxnId]/maximum
						if rxnId in revelantRxnsWithFlux and rxnId in rxnsInSol1: objective[y_f] = 0.2
						elif rxnId in revelantRxnsWithFlux: objective[y_f] = 0.18
						elif rxnId in rxnsInSol1: objective[y_f] = 0.15
						else: objective[y_f] = 0.05

				if y_r in solver2.neg_vars:
					if rxnsScoresMinimumMedia[rxnId] < 0:
						#if rxnId in rxnsInSol1: objective[y_r] = round(rxnsScoresMinimumMedia[rxnId],3)/2
						#else: objective[y_r] = round(rxnsScoresMinimumMedia[rxnId],3)*2
						if rxnId in revelantRxnsWithFlux and rxnId in rxnsInSol1: objective[y_r] = -1
						elif rxnId in revelantRxnsWithFlux: objective[y_r] = -1.1
						elif rxnId in rxnsInSol1: objective[y_r] = -1.2
						else: objective[y_r] = -1.5
					elif rxnsScoresMinimumMedia[rxnId] == 0: objective[y_r] = -0.0001
					else: 
						#objective[y_r] = rxnsScoresMinimumMedia[rxnId]/maximum
						if rxnId in revelantRxnsWithFlux and rxnId in rxnsInSol1: objective[y_r] = 0.2
						elif rxnId in revelantRxnsWithFlux: objective[y_r] = 0.18
						elif rxnId in rxnsInSol1: objective[y_r] = 0.15
						else: objective[y_r] = 0.05
			
						
			solver2.update()
			
			solver2.set_objective(linear=objective, minimize=False)
			
			print('\nStarting to solve 2.1: ' + str(datetime.now()))
			solution2 = solver2.solve(emphasis=1, timelimit=600)
			print('\nSolved 2.1: ' + str(datetime.now()))
			
			inSol2 = set()
			if solution2.status == Status.OPTIMAL:
				#make positive the score of reactions participating in this solution.
				for rxn in cobraModel.reactions:
					
					
					reframedRxnId = 'R_'+rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
				
					if ('yr_' + reframedRxnId in solution2.values and solution2.values['yr_' + reframedRxnId] > 1e-5) or ('yf_' + reframedRxnId in solution2.values and solution2.values['yf_' + reframedRxnId] > 1e-5) or (reframedRxnId in solution2.values and abs(solution2.values[reframedRxnId]) > 1e-5):
					
						inSol2.add(rxn.id)
						if '_forwardTemp' in rxn.id: inSol2.add(rxn.id.replace('_forwardTemp','_reverseTemp'))
						if '_reverseTemp' in rxn.id: inSol2.add(rxn.id.replace('_reverseTemp','_forwardTemp'))
					
						if rxnsScores[reframedRxnId] < 0: rxnsScores[reframedRxnId] = -0.01
					
			#find essential rxns in solution 2
			for rxn in cobraModel.reactions:
				if rxn.id not in inSol2: 
					rxn.lower_bound, rxn.upper_bound = 0, 0
					continue
			
			essentialInSolution2 = set()
			for rxn in cobraModel.reactions:
			
				if rxn.lower_bound == 0 and rxn.upper_bound == 0: continue

				lb, ub = rxn.lower_bound, rxn.upper_bound
				rxn.lower_bound, rxn.upper_bound = 0, 0
		
				fba_slim = cobraModel.slim_optimize()
				
				rxn.lower_bound, rxn.upper_bound = lb, ub
				if fba_slim == 0: essentialInSolution2.add(rxn)
				
			
			for rxn in essentialInSolution2:
				reframedRxnId = 'R_'+rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
				if rxnsScores[reframedRxnId] < 0: 
					rxnsScores[reframedRxnId] = 0.1
			
			
		
		for rxn in allIntake_bounds:
			rxn.lower_bound = allIntake_bounds[rxn]['lb']
			rxn.upper_bound = allIntake_bounds[rxn]['ub']
	

	#starting solution 3
	positiveInObjective = set()
	for rxnId in reframedModel.reactions:
		if rxnId in rxnsScores and rxnsScores[rxnId] >= 0: positiveInObjective.add(rxnId)
	

	solver3 = initSolver(reframedModel, positiveInObjective, constraintsFromFile)
	
	
	#creates objective function
	objective = {}
	for rxnId in reframedModel.reactions:
		
		y_r, y_f = 'yr_' + rxnId, 'yf_' + rxnId

		if y_f in solver3.pos_vars:
			if rxnsScores[rxnId] < 0:
				if rxnId in rxnsInSol1: objective[y_f] = round(rxnsScores[rxnId],3)/2
				else: objective[y_f] = round(rxnsScores[rxnId],3)*2
			elif rxnsScores[rxnId] == 0: objective[y_f] = -0.0001
			else:
				if rxnId in rxnsInSol1: objective[y_f] = round(rxnsScores[rxnId],3)*2
				else: objective[y_f] = round(rxnsScores[rxnId],3)

		if y_r in solver3.neg_vars:
			if rxnsScores[rxnId] < 0:
				if rxnId in rxnsInSol1: objective[y_r] = round(rxnsScores[rxnId],3)/2
				else: objective[y_r] = round(rxnsScores[rxnId],3)*2
			elif rxnsScores[rxnId] == 0: objective[y_r] = -0.0001
			else: 
				if rxnId in rxnsInSol1: objective[y_r] = round(rxnsScores[rxnId],3)*2
				else: objective[y_r] = round(rxnsScores[rxnId],3)
				
	solver3.update()
	
	solver3.set_objective(linear=objective, minimize=False)
	
	#solution3 = solver3.solve(timelimit=600)
	print('\nStarting to solve 3.1: ' + str(datetime.now()))
	solution3 = solver3.solve(timelimit=3600)
	if solution3.status != Status.OPTIMAL:
		print('\nStarting to solve 3.2: ' + str(datetime.now()))
		solution3 = solver3.solve(emphasis=1, timelimit=3600)
	print('\nSolved 3.1: ' + str(datetime.now()))
		
	saveProgressFile(94, outputfolder)
	
	if solution3.status == Status.OPTIMAL:
		status = saveReframedModel(reframedModel, solution3, passiveDiffusion, cobraModel, rheaIdToGene, rxnsScores, outputfile, spontaneousRheaRxns, biocycPathways, rxnsPerModules, outputfolder)
		return 1
	if solution3.status == Status.UNKNOWN:
		print('\nIt was not possible to solve the first optimization problem within the maximum time limit. The genome of the organism may be contaminated or incomplete. Alternatively, try including additional constraints or increasing the maximum processing time.')
	else:
		print('\nIt was not possible to find a solution to your problem. If you are using aditional constraints, try to relex them.')
		
	print('\nUsing solution of previous optimization')
		
	status = saveReframedModel(reframedModel, solution1, passiveDiffusion, cobraModel, rheaIdToGene, rxnsScores, outputfile, spontaneousRheaRxns, biocycPathways, rxnsPerModules, outputfolder)
	return 0
		
		

