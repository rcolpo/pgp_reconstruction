import pickle
from pgp_reconstruction import project_dir
import cobra
import datetime
import os
import copy
import numpy as np
import re

from pgp_reconstruction.reconstruction.makeEssentialRxns import makeEssentialGenesEssential

from reframed.solvers import solver_instance
from reframed.solvers.solver import VarType
from reframed.solvers.solution import Status



def initSolver(reframedModel, positiveInObjective, constraintsFromFile, cobraModel=None, firstOptimization=0, eps=1e-3, bigM=1000, min_growth=0.1):

	#solver5 = initSolver(reframedModel, positiveInSolution|positiveDiffusion, constraintsFromFile) # essentialRxnsInSpec = set()

	#multiple optimization problems are solved. To avoid repetition, a function is used to create the problem.

	solver = solver_instance(reframedModel)
	
	solver.neg_vars = []
	solver.pos_vars = []
	for rxnId in reframedModel.reactions: #reactions:
	
		y_r, y_f = 'yr_' + rxnId, 'yf_' + rxnId
	
		if reframedModel.reactions[rxnId].lb < 0:
			solver.add_variable(y_r, 0, 1, vartype=VarType.BINARY, update=False)
			solver.neg_vars.append(y_r)
		if reframedModel.reactions[rxnId].ub > 0:
			solver.add_variable(y_f, 0, 1, vartype=VarType.BINARY, update=False)
			solver.pos_vars.append(y_f)
	solver.update()
	
	
	if firstOptimization:
	
		#force flux on essentialRxns, to reduce solution space
	
		#try to make flux on essentialRxns more desirable
		with cobraModel as model:
			
			boundaryIds = {rxnCobra.id for rxnCobra in model.boundary}
		
			for rxnCobra in model.reactions:
			
				if len(rxnCobra.compartments) == 2: continue
				if rxnCobra.id in boundaryIds: continue
				if rxnCobra.genes: continue
				if rxnCobra.lower_bound < 0: rxnCobra.lower_bound = -0.1
				if rxnCobra.upper_bound > 0: rxnCobra.upper_bound = 0.1
	
			fba = model.optimize()
		
		essentialRxnsWithFlux = set()
		for rxnCobra in cobraModel.reactions:
			if abs(fba.fluxes[rxnCobra.id]) > 0 and rxnCobra.genes: 
				if 'R_' + rxnCobra.id not in positiveInObjective: positiveInObjective.append('R_' + rxnCobra.id)
				essentialRxnsWithFlux.add(rxnCobra.id.replace('_forwardTemp',''.replace('_reverseTemp','')))
	
		#for rxnId in reframedModel.reactions:
		
		for rnxId in essentialRxnsWithFlux:
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
				solver.add_constraint('forceEssential_' + rnxId, constraintsEssentials, '=', 1, update=False)
	solver.update()
	
	biomass = reframedModel.biomass_reaction
	#upper and lowe flux constraints:
	for rxnId in reframedModel.reactions:
	
		y_r, y_f = 'yr_' + rxnId, 'yf_' + rxnId
	
		if reframedModel.reactions[rxnId].lb < 0 or reframedModel.reactions[rxnId].ub > 0:

			if rxnId in positiveInObjective:# and (reframedModel.reactions[rxnId].lb < 0 and reframedModel.reactions[rxnId].ub > 0):
				#allow reverse and forward simultaneously if the score is negative
				
				if y_r in solver.neg_vars and '_reverseTemp' in y_r and y_f.replace('_reverseTemp','_forwardTemp') in solver.pos_vars: #when universal model is split
					solver.add_constraint('avoidTwoDirections_' + rxnId, {y_r: 1, 'yf_' + rxnId.replace('_reverseTemp','_forwardTemp'): 1}, '<', 1, update=False) 
				
				if y_r in solver.neg_vars and y_f in solver.pos_vars:
					solver.add_constraint('avoidTwoDirections_' + rxnId, {y_r: 1, y_f: 1}, '<', 1, update=False)
		
			#force flux on biomass
			if rxnId == reframedModel.biomass_reaction:
				if y_f in solver.pos_vars:
					solver.add_constraint('forceBiomass_' + rxnId, {y_f: 1}, '=', 1, update=False)
				if y_r in solver.neg_vars:
					solver.add_constraint('forceBiomass_' + rxnId, {y_r: 1}, '=', 1, update=False)
				solver.add_constraint('min_growth' + rxnId, {rxnId: 1}, '>', min_growth, update=False)
					
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

	constraintsFromFile = {'forward':{'soft':dict(),'hard':dict()},'reverse':{'soft':dict(),'hard':dict()}}
	positiveInObjective =  activeRxnsreframedId3
	
	solver5 = initSolver(reframedModel, positiveInObjective, constraintsFromFile)

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
	
	print('\nStart solving 5: ' + str(datetime.datetime.now()))
	solution5 = solver5.solve(emphasis=2, timelimit=600)
	
	if solution5.status != Status.OPTIMAL:
		print('\nStarting to solve 3.2: ' + str(datetime.datetime.now()))
		solution5 = solver5.solve(timelimit=3600)
	
	print('\nSolved 5: ' + str(datetime.datetime.now()))
	
	with open('solution5.pickle', 'wb') as handle:
			pickle.dump(solution5, handle, protocol=4)

	return solution5
	
	
def includeGenesRules(cobraModel, rheaIdToGene):

	##adiciona informacao de genes
	rxnsAndGenes = dict()
	for rxn in cobraModel.reactions:
		
		for eachDict in rheaIdToGene['bestPerReadSimplified']:
			if 'R_' + rxn.id in eachDict['rxns']:
				if rxn.id not in rxnsAndGenes:
					rxnsAndGenes[rxn.id] = set()
				rxnsAndGenes[rxn.id].add(eachDict['gene'])
				
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
			ruleString += uniprotId + ' or '
		ruleString = ruleString[:-4] + ')'
		cobraModel.reactions.get_by_id(rxnId).gene_reaction_rule = ruleString
		
		
	genesToRxns = dict()
	for rxnId in rxnsAndGenes:
		for gene in rxnsAndGenes[rxnId]:
			if gene not in genesToRxns: genesToRxns[gene] = list()
			genesToRxns[gene].append(rxnId)
	

	essentialRxns_Teoretical = set()
	for rxnId in rxnsAndGenes:
		essentialRxns_Teoretical.add(cobraModel.reactions.get_by_id(rxnId))
		
	essentialRxnsIds_Teoretical = set()
	for rxnId in rxnsAndGenes:
		essentialRxnsIds_Teoretical.add(rxnId)
	
	return genesToRxns, essentialRxns_Teoretical, essentialRxnsIds_Teoretical
	
	

def savereframeddModel(reframedModel, solution3, passiveDiffusion, cobraModel, rheaIdToGene, uniprotToRheaRxns, rxnsScores, inputFileName, spontaneousRheaRxns, constraintsFromFile):
	#status = savereframeddModel(reframedModel, solution1, passiveDiffusion, cobraModel, rheaIdToGene, uniprotToRheaRxns, rxnsScores, inputFileName, spontaneousRheaRxns, constraintsFromFile)
	
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
				if len(rxn.compartments) == 2:
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

	finalModel.slim_optimize()

	
	for sufix in [''] + list(range(1,100)):
		if os.path.isfile(os.path.splitext(os.path.basename(inputFileName))[0] + str(sufix) + ".xml"): continue
		else: 
			cobra.io.write_sbml_model(finalModel, os.path.splitext(os.path.basename(inputFileName))[0] + str(sufix) + ".xml")
			break
	return 1
		


def prune_model(reframedModel, cobraModel, rxnsScores, rheaIdToGene, uniprotToRheaRxns, inputFileName, rxnsCobraToreframedModel, soft_constraints_pathways, bestPerReadSimplified,
				taxoOfTarget, rheaWithSameSyn, softPositiveNewRxn, constraintsFromFile, rxnsFromReference, relevantRxnsPerGeneInModel,
				min_growth=0.1, min_atpm=0.1, eps=1e-3, bigM=1e3, default_score=-6.0,
				soft_score=10.0, ref_score=0.0, solver=None, debug_output=None, rsink = None, minimumMedia=False):
				
				


	""" Reconstruct a metabolic reframedModel using the reframedMe approach.

	Args:
		reframedModel (CBModel): universal reframedModel
		reaction_scores (pandas.DataFrame): reaction rxnsScores
		outputfile (str): write reframedModel to SBML file (optional)
		flavor (str): SBML flavor ('cobra' or 'fbc2', optional)
		default_score (float): penalty for non-annotated intracellular reactions (default: -1.0)
		soft_score (float): score for soft constraints (default: 1.0)
		init_env (Environment): initialize final reframedModel with given Environment (optional)

	Returns:
		CBModel: reconstructed reframedModel
	"""
		

	pickle_file_path = os.path.join(project_dir, 'data/generated', 'spontaneousRedundant.pickle')
	with open(pickle_file_path, 'rb') as f:
		spontaneousRheaRxns = pickle.load(f)
		
	pickle_file_path = os.path.join(project_dir, 'data/generated', 'superPathsDict.pickle')
	with open(pickle_file_path, 'rb') as f:
		superPathsDict = pickle.load(f)
		
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
		
	
		

	metsFromPseudomonas = {}

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

	genesToRxns, essentialRxns_Teoretical, essentialRxnsIds_Teoretical = includeGenesRules(cobraModel, rheaIdToGene)	
		
	rheaRxnsAndSyn = dict()
	biomass = reframedModel.biomass_reaction
	for rxn in cobraModel.reactions:
		idInreframed = 'R_' + rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
		rheaRxnsAndSyn[idInreframed] = {idInreframed[2:]}
		if 'rhea' in rxn.annotation:
			for rheaRxn in rxn.annotation['rhea']:
				rheaRxnsAndSyn[idInreframed].add(rheaRxn)

	rxnsScores[biomass] = 50

	#facilitates the inclusion of reactions of transport by passive diffusion 
	passiveDiffusion = set()
	for rxnCobra in cobraModel.reactions:
	
		idInreframed = 'R_' + rxnCobra.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
			
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
				if met.formula_weight < 50 and met.charge == 0: passiveDiffusion.add(met.id)
		
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
		if len(rxnCobra.compartments) == 2:
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
	
	
	##include soft constraints in objective
	#format: constraintsFromFile = {'forward':{'soft':dict(),'hard':dict()},'reverse':{'soft':dict(),'hard':dict()}}
	for cobraId in constraintsFromFile['forward']['soft']:
		rxnId = 'R_'+cobraId.replace('-','__45__').replace('.','__46__').replace('+','__43__')
		if reframedModel.reactions[rxnId].ub > 0:
			if rxnsScores[rxnId] > 0 and constraintsFromFile['forward']['soft'][cobraId] > 0 and rxnsScores[rxnId] > constraintsFromFile['forward']['soft'][cobraId]: continue
			if rxnsScores[rxnId] < 0 and constraintsFromFile['forward']['soft'][cobraId] < 0 and rxnsScores[rxnId] < constraintsFromFile['forward']['soft'][cobraId]: continue
			rxnsScores[rxnId] = constraintsFromFile['forward']['soft'][cobraId]
			sucesso = 1
	for cobraId in constraintsFromFile['forward']['hard']:
		rxnId = 'R_'+cobraId.replace('-','__45__').replace('.','__46__').replace('+','__43__')
		if reframedModel.reactions[rxnId].ub > 0:
			if constraintsFromFile['forward']['hard'][cobraId] < 0: rxnsScores[rxnId] = -100
			if constraintsFromFile['forward']['hard'][cobraId] > 0: rxnsScores[rxnId] = 100
			sucesso = 1	
	for cobraId in constraintsFromFile['reverse']['soft']:
		rxnId = 'R_'+cobraId.replace('-','__45__').replace('.','__46__').replace('+','__43__')
		if reframedModel.reactions[rxnId].lb < 0:
			if rxnsScores[rxnId] > 0 and constraintsFromFile['reverse']['soft'][cobraId] > 0 and rxnsScores[rxnId] > constraintsFromFile['reverse']['soft'][cobraId]: continue
			if rxnsScores[rxnId] < 0 and constraintsFromFile['reverse']['soft'][cobraId] < 0 and rxnsScores[rxnId] < constraintsFromFile['reverse']['soft'][cobraId]: continue
			rxnsScores[rxnId] = constraintsFromFile['reverse']['soft'][cobraId]
			sucesso = 1
	for cobraId in constraintsFromFile['reverse']['hard']:
		rxnId = 'R_'+cobraId.replace('-','__45__').replace('.','__46__').replace('+','__43__')
		if reframedModel.reactions[rxnId].lb < 0:
			if constraintsFromFile['reverse']['hard'][cobraId] < 0: rxnsScores[rxnId] = -100
			if constraintsFromFile['reverse']['hard'][cobraId] > 0: rxnsScores[rxnId] = 100
			sucesso = 1
	
	
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
	

	#use rxnsFromReference to increase scores
	for rxn in cobraModel.reactions:
		idInreframed = 'R_'+rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
		for db in rxn.annotation:
		
			if rxnsFromReference['rxnWithGenes'].intersection(rxn.annotation[db]):
				if idInreframed not in rxnsScores or rxnsScores[idInreframed] < 0.1: 
					print(rxn.id + ' ' + str(rxnsFromReference['rxnWithGenes'].intersection(rxn.annotation[db])))
					rxnsScores[idInreframed] = 0.1
				
			if rxnsFromReference['rxnWithoutGenes'].intersection(rxn.annotation[db]):
				if rxnsScores[idInreframed] < -0.1: 
					rxnsScores[idInreframed] = -0.1
				
			#if rxnsFromReference['rxnBoundary'].intersection(rxn.annotation[db]):
			#	if rxnsScores[idInreframed] < 0: 
			#		rxnsScores[idInreframed] = rxnsScores[idInreframed]/2
	

	print('\nLine 626: ' + str(datetime.datetime.now()))

	positiveInObjective = []
	solver1 = initSolver(reframedModel, positiveInObjective, constraintsFromFile, cobraModel=cobraModel, firstOptimization=1)
	
	#creates objective function
	objective = {}
	for rxnId in reframedModel.reactions:
		
		y_r, y_f = 'yr_' + rxnId, 'yf_' + rxnId

		if y_f in solver1.pos_vars:
			if rxnsScores[rxnId] < 0: objective[y_f] = round(rxnsScores[rxnId],3)*3
			elif rxnsScores[rxnId] == 0: objective[y_f] = -0.0001
			else: objective[y_f] = round(rxnsScores[rxnId],3)

		if y_r in solver1.neg_vars:
			if rxnsScores[rxnId] < 0: objective[y_r] = round(rxnsScores[rxnId],3)*3
			elif rxnsScores[rxnId] == 0: objective[y_r] = -0.0001
			else: objective[y_r] = round(rxnsScores[rxnId],3)
	
	
	solver1.set_objective(linear=objective, minimize=False)
	print('\nStarting to solve 1: ' + str(datetime.datetime.now()))
	solution1 = solver1.solve(emphasis=1, timelimit=300)#if does not find in 5 minutes, give up
	
	if solution1.status == Status.OPTIMAL: 
		print('Solved 1: ' + str(datetime.datetime.now()) + '\n')
	elif solution1.status == Status.UNKNOWN:
		print('\nIt was not possible to solve the first optimization problem within the maximum time limit. The genome of the organism may be contaminated, incomplete or have the wring taxonomical assigment. Alternatively, try including additional constraints (like providing information on the culture medium) or increasing the maximum processing time.')
		return
	else:
		print('\nIt was not possible to find a solution to your problem. If you are using aditional constraints, try to relex them.')
		return

	with open('solution1.pickle', 'wb') as handle:
		pickle.dump(solution1, handle, protocol=4)

	
	#start second optimization	
	
	#makeEssentialGenesEssential specific for target organism
	deleted = makeEssentialGenesEssential(solution1, bestPerReadSimplified, rheaIdToGene, uniprotToRheaRxns, cobraModel, reframedModel, relevantRxnsPerGeneInModel)
	cobraModel.remove_reactions(deleted)
	toDeletereframed = ['R_' + rxnId.replace('-','__45__').replace('.','__46__').replace('+','__43__') for rxnId in deleted]
	reframedModel.remove_reactions(toDeletereframed)

	#searching for the pathways present in solution1
	rxnsInSol1 = set()
	for rxnId in reframedModel.reactions:
		if ('yr_' + rxnId in solution1.values and solution1.values['yr_' + rxnId] > 1e-5) or ('yf_' + rxnId in solution1.values and solution1.values['yf_' + rxnId] > 1e-5) or (abs(solution1.values[rxnId]) > 1e-5):
			rxnsInSol1.add(rxnId.replace('_forwardTemp','').replace('_reverseTemp',''))
		
	#starting solution 3
	positiveInObjective = set()
	for rxnId in reframedModel.reactions:
		if rxnId in rxnsScores and rxnsScores[rxnId] >= 0: positiveInObjective.add(rxnId)
	

	#solver3 = initSolver(reframedModel, positiveInObjective, constraintsFromFile, essentialRxnsInSpec=essentialRxnsInSpec.intersection(rxnsInSol2OriginalIds))
	solver3 = initSolver(reframedModel, positiveInObjective, constraintsFromFile)
	
	
	#creates objective function
	objective = {}
	for rxnId in reframedModel.reactions:
		
		y_r, y_f = 'yr_' + rxnId, 'yf_' + rxnId

		if y_f in solver3.pos_vars:
			if rxnsScores[rxnId] < 0:
				if rxnId.replace('_forwardTemp','').replace('_reverseTemp','') in rxnsInSol1: objective[y_f] = round(rxnsScores[rxnId],3)/2
				else: objective[y_f] = round(rxnsScores[rxnId],3)
			elif rxnsScores[rxnId] == 0: objective[y_f] = -0.0001
			else:
				if rxnId.replace('_forwardTemp','').replace('_reverseTemp','') in rxnsInSol1: objective[y_f] = round(rxnsScores[rxnId],3)*2
				else: objective[y_f] = round(rxnsScores[rxnId],3)

		if y_r in solver3.neg_vars:
			if rxnsScores[rxnId] < 0:
				if rxnId.replace('_forwardTemp','').replace('_reverseTemp','') in rxnsInSol1: objective[y_r] = round(rxnsScores[rxnId],3)/2
				else: objective[y_r] = round(rxnsScores[rxnId],3)
			elif rxnsScores[rxnId] == 0: objective[y_r] = -0.0001
			else: 
				if rxnId.replace('_forwardTemp','').replace('_reverseTemp','') in rxnsInSol1: objective[y_r] = round(rxnsScores[rxnId],3)*2
				else: objective[y_r] = round(rxnsScores[rxnId],3)
				
	print('\nLine 651: ' + str(datetime.datetime.now()))
				
	solver3.update()
	
	solver3.set_objective(linear=objective, minimize=False)
	
	print('\nStarting to solve 3: ' + str(datetime.datetime.now()))
	
	solution3 = solver3.solve(timelimit=600)
	if solution3.status != Status.OPTIMAL:
		print('\nStarting to solve 3.2: ' + str(datetime.datetime.now()))
		solution3 = solver3.solve(emphasis=1, timelimit=3600)

	with open('solution3.pickle', 'wb') as handle:
		pickle.dump(solution3, handle, protocol=4)
		
	
	if solution3.status == Status.OPTIMAL:
		#send solution3 to create new model
		status = savereframeddModel(reframedModel, solution3, passiveDiffusion, cobraModel, rheaIdToGene, uniprotToRheaRxns, rxnsScores, inputFileName, spontaneousRheaRxns, constraintsFromFile)
		return 1
	if solution3.status == Status.UNKNOWN:
		print('\nIt was not possible to solve the first optimization problem within the manixum time limit. The genome of the organism may be contaminated or incomplete. Alternatively, try including additional constraints or increasing the maximum processing time.')
	else:
		print('\nIt was not possible to find a solution to your problem. If you are using aditional constraints, try to relex them.')
		
	print('\nUsing solution of previous optimization')
		
	status = savereframeddModel(reframedModel, solution1, passiveDiffusion, cobraModel, rheaIdToGene, uniprotToRheaRxns, rxnsScores, inputFileName, spontaneousRheaRxns, constraintsFromFile)
	return 0
		
		

