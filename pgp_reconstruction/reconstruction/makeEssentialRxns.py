import cobra
import numpy as np
import itertools
import os

from reframed.solvers.solution import Status
from pgp_reconstruction import project_dir
	
	
def includeGenesRules(cobraModel, rheaIdToGene, uniprotToRheaRxns):

	#para uma solucao da otimizacao
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
	
	



def findEssencialGenes(bestPerReadSimplified, reframedModel):
	#findEssencialGenes(bestPerReadSimplified, cobraModel)

	specGenAndFreq = dict()

	for file in ['deg_annotation_p.csv', 'deg_annotation_e.csv', 'deg_annotation_a.csv']:

		with open(os.path.join(project_dir, 'data/input', file), encoding="utf8") as f:
			fileLines = f.readlines()
		
		for line in fileLines:
			line = line.replace(';;',';"-";')
			genes = line.lower().split('";"')[2]
			species = line.split('";"')[7]
			if ' ' not in species: continue
			species = species.split(' ')[0]
			
			genesList = genes.replace('"','').split('/')
			for gene in genesList:
				if gene == '-': continue
				if species not in specGenAndFreq: specGenAndFreq[species] = dict()
				if gene not in specGenAndFreq[species]: specGenAndFreq[species][gene] = 0
				specGenAndFreq[species][gene] += 1
				
	essentialGenes = set()
	for species in specGenAndFreq:
		for gene in specGenAndFreq[species]:
			if specGenAndFreq[species][gene] <= 2: continue
			essentialGenes.add(gene)
			

	essentialRxnsInSpec = set()
	essentialGenesInSpec = set()
	for dictInList in bestPerReadSimplified:
		if dictInList['gene'].lower() in essentialGenes:
			essentialGenesInSpec.add(dictInList['gene'].lower())
			for rxnId in dictInList['rxns']:
				essentialRxnsInSpec.add(rxnId)
			
	return essentialGenesInSpec
		

def remove_high_outliers(data):
	q1 = np.percentile(data, 25)
	q3 = np.percentile(data, 75)
	iqr = q3 - q1
	upper_bound = q3 + iqr

	return [x for x in data if x <= upper_bound]




def makeEssentialGenesEssential(solution1, bestPerReadSimplified, rheaIdToGene, uniprotToRheaRxns, cobraModel, reframedModel, relevantRxnsPerGeneInModel):
	

	deleted = []
	if solution1.status == Status.OPTIMAL:


		genesToRxns, essentialRxns_Teoretical, essentialRxnsIds_Teoretical = includeGenesRules(cobraModel, rheaIdToGene, uniprotToRheaRxns)
		essentialGenesInSpec = findEssencialGenes(bestPerReadSimplified, cobraModel)
		
		rxnsInSol1 = set()
		for rxnId in reframedModel.reactions:
			if ('yr_' + rxnId in solution1.values and solution1.values['yr_' + rxnId] > 0) or ('yf_' + rxnId in solution1.values and solution1.values['yf_' + rxnId] > 0) or (abs(solution1.values[rxnId]) > 0):
				#cobraRxn = rxnId[2:].replace('__45__','-').replace('__46__','.').replace('__43__','+').replace('_forwardTemp','').replace('_reverseTemp','')
				cobraRxn = rxnId[2:].replace('__45__','-').replace('__46__','.').replace('__43__','+')
				rxnsInSol1.add(cobraRxn)
				rxnsInSol1.add(cobraRxn.replace('_forwardTemp','_reverseTemp'))
				rxnsInSol1.add(cobraRxn.replace('_reverseTemp','_forwardTemp'))
				
		rxnsWithGoodEvidence = set()		
		for i in relevantRxnsPerGeneInModel:
			rxnsWithGoodEvidence.add(rxnId[2:].replace('__45__','-').replace('__46__','.').replace('__43__','+'))
		

		##Create smaller model, to speedUp process
		mets = set()
		rxns = set()
		for rxn in cobraModel.reactions:
			if rxn.id in rxnsInSol1:
				rxns.add(rxn.id)
				for met in rxn.metabolites:
					mets.add(met.id)

		cobraModel2 = cobra.Model('Model from ')
		cobraModel2.tolerance = 1e-04
		
		#metabolitos
		for met in cobraModel.metabolites:
			if met.id not in mets: continue
			my_metabolite = cobra.Metabolite(met.id)
			my_metabolite.name = met.name
			my_metabolite.formula = met.formula
			my_metabolite.charge = met.charge
			my_metabolite.compartment = met.compartment
			my_metabolite.annotation = met.annotation
			_ = cobraModel2.add_metabolites([my_metabolite])
			
		#reacoes
		for rxn in cobraModel.reactions:
			if rxn.id not in rxns: continue
			
			my_reaction = cobra.Reaction(rxn.id)
			my_reaction.name = rxn.name
			my_reaction.annotation = rxn.annotation
			my_reaction.gene_reaction_rule = rxn.gene_reaction_rule
			_ = cobraModel2.add_reactions([my_reaction])
			my_reaction.reaction = rxn.reaction
			my_reaction.objective_coefficient = rxn.objective_coefficient
			my_reaction.lower_bound = rxn.lower_bound
			my_reaction.upper_bound = rxn.upper_bound

		cobraModel2 = cobra.flux_analysis.fastcc(cobraModel2)
		
		#save flux values. to use later
		initialFlux = dict()
		initialFluxIds = dict()
		for rxn in cobraModel2.reactions:
			initialFlux[rxn] = {'lb': rxn.lower_bound,'ub': rxn.upper_bound}
			initialFluxIds[rxn.id] = {'lb': rxn.lower_bound,'ub': rxn.upper_bound}
			

		#do not wich to minimize boundaries or transport
		boundaries = set(rxn.id for rxn in cobraModel.boundary)
		transport = set(rxn.id for rxn in cobraModel.reactions if len(rxn.compartments) == 2)
		acidDissociation = set(rxn.id for rxn in cobraModel.reactions if 'acidDissociation' in rxn.id)

		conter = 0
		while True:
			#find essential reactions:
			conter += 1
			print('In loop ' + str(conter) + '\n')
			activeRxns = set(rxn.id for rxn in cobraModel2.reactions)
			blockedInitial = cobra.flux_analysis.find_blocked_reactions(cobraModel2, list(set(cobraModel2.reactions) - set(cobraModel2.boundary)), processes=1, zero_cutoff=1e-4)
			singleDeletion = cobra.flux_analysis.single_reaction_deletion(cobraModel2, activeRxns)
			
			essentialRxns = set()
			for index, row in singleDeletion.iterrows():
				#to deal of different versions of cobraPy
				try:
					if abs(row['growth']) < 1e-5: essentialRxns.add(list(index)[0])
				except:
					if abs(row['growth']) < 1e-5: essentialRxns.add(list(row['ids'])[0])
					
			
			#exclude duplicates, to speed up the process of finding essential rxns
			rxnsList = list()
			geneList = list()
			for gene in genesToRxns:
				
				if gene not in essentialGenesInSpec: continue
				if essentialRxns.intersection(genesToRxns[gene]): continue
				if not activeRxns.intersection(genesToRxns[gene]): continue
				
				if activeRxns.intersection(genesToRxns[gene]) not in rxnsList: 
					rxnsList.append(activeRxns.intersection(genesToRxns[gene]))
					geneList.append(gene)
				
			genesToRxns2 = dict()
			for count, eachSet1 in enumerate(rxnsList):
				toTnclude = 1
				for eachSet2 in rxnsList: 
					if eachSet1 == eachSet2: continue
					if eachSet1.issubset(eachSet2): toTnclude = 0
				if toTnclude == 1: 
					genesToRxns2[geneList[count]] = eachSet1
			
			#execute simple deletion, to find potential secundary route that should not exist
			
			essentialPerGene = dict()
			for gene in genesToRxns2:
				
				for rxnId in genesToRxns2[gene]:
					rxn = cobraModel2.reactions.get_by_id(rxnId)
					rxn.lower_bound = 0
					rxn.upper_bound = 0
					
					
				fbanterior = cobraModel2.optimize()
				
				if fbanterior.status == 'optimal' and fbanterior.objective_value > 0.00001:
					
					potentialToDel = set()
					for rxn in cobraModel2.reactions:
						if rxn.id in fbanterior.fluxes and (fbanterior.fluxes[rxn.id] > 0): potentialToDel.add(rxn.id)
					
					toRemove = set()
					singleDeletion = cobra.flux_analysis.single_reaction_deletion(cobraModel2, potentialToDel - (essentialRxnsIds_Teoretical|rxnsWithGoodEvidence|essentialRxns|boundaries|transport|acidDissociation), processes=1)
					for index, row in singleDeletion.iterrows():
						try:
							if abs(row['growth']) < 1e-5: toRemove.add(cobraModel2.reactions.get_by_id(list(index)[0]))
						except:
							if abs(row['growth']) < 1e-5: toRemove.add(cobraModel2.reactions.get_by_id(list(row['ids'])[0]))

					if toRemove:
						print('len(toRemove) = ' + str(len(toRemove)))
						essentialPerGene[gene] = toRemove
					
				#restore fluxes
				for rxn in initialFlux:
					if rxn.lower_bound != initialFlux[rxn]['lb']: rxn.lower_bound = initialFlux[rxn]['lb']
					if rxn.upper_bound != initialFlux[rxn]['ub']: rxn.upper_bound = initialFlux[rxn]['ub']
					
			if not essentialPerGene: break
					
			essentialPerGeneLen = list()
			for gene in essentialPerGene:
				essentialPerGeneLen.append(len(essentialPerGene[gene]))
			goodLen = remove_high_outliers(essentialPerGeneLen)
			
			filtered_essentialPerGene = {k: v for k, v in essentialPerGene.items() if len(v) in goodLen}
					

			rxnsSolutionList = list()
			for rxnSet in filtered_essentialPerGene.values():
				for rxn in rxnSet: rxnsSolutionList.append(rxn)
				
			
			#check which of the possible reactions to delete result in less blocked reactions
			numberOfBlocked = dict()
			for rxn in set(rxnsSolutionList):
				with cobraModel2 as model:
					
					model.reactions.get_by_id(rxn.id).lower_bound=0
					model.reactions.get_by_id(rxn.id).upper_bound=0
					
					blocked = cobra.flux_analysis.find_blocked_reactions(model, list(set(model.reactions) - set(model.boundary)), processes=1, zero_cutoff=1e-4)
					numberOfBlocked[rxn.id] = blocked
					
			onePergene = list()
			for gene in filtered_essentialPerGene:
				tempSet = set()
				for rxn in filtered_essentialPerGene[gene]:
					if len(numberOfBlocked[rxn.id]) == 1:
						tempSet.add(rxn.id)
				if tempSet and tempSet not in onePergene: onePergene.append(tempSet)
			
			sucesso = 0
			combinations = list(itertools.product(*onePergene))
			for combination in combinations:
				if not combination: continue
				with cobraModel2 as model:
					for rxnId in combination:
						model.reactions.get_by_id(rxnId).lower_bound=0
						model.reactions.get_by_id(rxnId).upper_bound=0
				
					blocked = cobra.flux_analysis.find_blocked_reactions(model, list(set(model.reactions) - set(model.boundary)), processes=1, zero_cutoff=1e-4)
					
					if len(blocked) - len(blockedInitial) <= len(combination): 
						print('Deleting: ' + str(combination))
						sucesso = 1
						toDelete = combination
						break
			
			if sucesso == 1:
				print(str(toDelete) + ' was deleted to make genes, essential.' )
			
			if sucesso == 0:

				#delete reaction which removal blokes the least  number of reactions
				rxnInSolFreq = dict()
				for rxn in set(rxnsSolutionList):
					if rxnsSolutionList.count(rxn) not in rxnInSolFreq: rxnInSolFreq[rxnsSolutionList.count(rxn)] = list()
					rxnInSolFreq[rxnsSolutionList.count(rxn)].append(rxn)
					
				possibleToDelete =  rxnInSolFreq[max(rxnInSolFreq)]
				
				#check which of the possible reactions to delete result in less blocked reactions
				blockedPerSolution = list()
				for rxn in possibleToDelete:
					with cobraModel2 as model:
						
						model.reactions.get_by_id(rxn.id).lower_bound=0
						model.reactions.get_by_id(rxn.id).upper_bound=0
						
						blocked = cobra.flux_analysis.find_blocked_reactions(model, processes=1, zero_cutoff=1e-4)
						blockedPerSolution.append(blocked)
					
				minIndex = 0
				for count, eachSet in enumerate(blockedPerSolution):
					if count == minIndex: continue
					if len(eachSet) < len(blockedPerSolution[minIndex]): minIndex = count
				
				toDelete = [possibleToDelete[minIndex].id]
				
				print(str(toDelete) + ' was deleted to make ' + str(max(rxnInSolFreq)) + ' genes, essential.' )
				
			#remove reactions
			deleted += toDelete
			cobraModel2.remove_reactions(toDelete)
			
	print('deleted = ' + str(deleted))
	return deleted
