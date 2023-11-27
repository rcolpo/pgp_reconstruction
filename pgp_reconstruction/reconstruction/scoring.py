import pickle
from pgp_reconstruction import project_dir
from statistics import median
import sys
import os
from datetime import datetime
from pgp_reconstruction.cli.util import saveProgressFile

try:from pgp_reconstruction.dependencies.MinPath_master.MinPath import MinPathMain
except: pass


def useReferenceModelData(reference, referenceScore, cobraModel, rxnsScores):
	#include in rxnsScores reactions from reference model
	if reference:
		try:
			cobraModelReference = cobra.io.read_sbml_model(reference)
			
			for cobraObject in [cobraModelReference.reactions, cobraModelReference.metabolites]:
				for rxn in cobraObject:
					for db in rxn.annotation:
						if type(rxn.annotation[db]) == type(''):
							rxn.annotation[db] = [rxn.annotation[db]]
			
		except IOError:
			raise IOError(f'Failed to load reference model')
		
		rxnWithGenes = set()
		rxnBoundary = set()
		rxnWithoutGenes = set()
		
		for rxn in cobraModelReference.reactions:

			if rxn.genes:
				for db in rxn.annotation:
					for synId in rxn.annotation[db]: rxnWithGenes.add(synId)
				rxnWithGenes.add(rxn.id)
					
			elif rxn in cobraModelReference.boundary:
				for db in rxn.annotation:
					for synId in rxn.annotation[db]: rxnBoundary.add(synId)
				rxnBoundary.add(rxn.id)

			else:
				for db in rxn.annotation:
					for synId in rxn.annotation[db]: rxnWithoutGenes.add(synId)
				rxnWithoutGenes.add(rxn.id)
					
		rxnsFromReference = {'rxnWithGenes':rxnWithGenes, 'rxnBoundary':rxnBoundary, 'rxnWithoutGenes':rxnWithoutGenes}
		
		
		#use rxnsFromReference to increase scores
		for rxn in cobraModel.reactions:
			idInreframed = 'R_'+rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
			for db in rxn.annotation:
			
				if rxnsFromReference['rxnWithGenes'].intersection(rxn.annotation[db]):
					if idInreframed not in rxnsScores or rxnsScores[idInreframed] < 0.1: 
						print(rxn.id + ' ' + str(rxnsFromReference['rxnWithGenes'].intersection(rxn.annotation[db])))
						rxnsScores[idInreframed] = referenceScore
					
				if rxnsFromReference['rxnWithoutGenes'].intersection(rxn.annotation[db]):
					if referenceScore > 0:
						if rxnsScores[idInreframed] < -0.1: rxnsScores[idInreframed] = -0.1
						else: rxnsScores[idInreframed] = referenceScore

def findPrioritary(top50Simplified, categoriesCount, rxnsInTaxonomyConstraints, rxnsFromUniprot, rxnsInPathsFromSoft, rxnsInPathsNaive, subunits, geneAndProteinNamePerSeqId, swiss90tremble50SeqInfo, swissProtIds, firstloop = 1):
	#first check if there is any match from rxnsInTaxonomyConstraints
	filtered1 = list()
	
	#try to find the missing subunits
	if not filtered1 and firstloop == 1 and subunits:

		for eachDict in top50Simplified:
			if 'subunit' not in swiss90tremble50SeqInfo[eachDict['uniprotEntry']]['title'] and 'component' not in swiss90tremble50SeqInfo[eachDict['uniprotEntry']]['title']:
				sucesso = 0
				for eachDict in top50Simplified:
					subunitType = ''
					for subunitString in ['subunit', 'component']:
						if subunitString in swiss90tremble50SeqInfo[eachDict['uniprotEntry']]['title']:
							#to extract the subunit part
							if subunitString in swiss90tremble50SeqInfo[eachDict['uniprotEntry']]['title']:
							
								titleProvisional = swiss90tremble50SeqInfo[eachDict['uniprotEntry']]['title']
							
								if titleProvisional.count(subunitString) > 1: 
									titleProvisional = titleProvisional.replace(subunitString,'',1)
							
								if '-'+subunitString in titleProvisional.lower():
									for word in titleProvisional.lower().split(' '):
										if '-'+subunitString in word: break
									subunitType = word.replace('-'+subunitString,'')
								else:
								
									if '/'+subunitString in titleProvisional.lower():
										toReplace = titleProvisional.lower().split('/'+subunitString)[0]
										titleProvisional = titleProvisional.replace(toReplace.split(' ')[-1]+'/','')
									
									if subunitString+' ' in titleProvisional.lower():
										subunitType = titleProvisional.lower().split(subunitString+' ')[1]
										subunitType = subunitType.split(' ')[0]
									elif ' '+subunitString in titleProvisional.lower():
										subunitType = titleProvisional.lower().split(' '+subunitString)[0]
										subunitType = subunitType.split(' ')[-1]
					
					if subunitType:
						for subunit in subunits:
							if subunit['rxns'].intersection(eachDict['rxns']) and subunit['gene'] != eachDict['gene'] and subunit['subunit'] != subunitType:
								categoriesCount['subunits'] += 1
								filtered1.append(eachDict)
	

	if not filtered1 and firstloop == 1:
		for eachDict in top50Simplified:
			if eachDict['rxns'].intersection(rxnsFromUniprot):
				categoriesCount['rxnsFromUniprot'] += 1
				filtered1.append(eachDict)
	
	if not filtered1 and firstloop == 1:
		for eachDict in top50Simplified:
			if eachDict['rxns'].intersection(rxnsInTaxonomyConstraints):
				categoriesCount['fromMinPath'] += 1
				filtered1.append(eachDict)
				
	maxScore = 0
	for eachDict in top50Simplified:
		if eachDict['score'] >= maxScore: maxScore = eachDict['score']

	if not filtered1 and firstloop == 1:
		for eachDict in top50Simplified:
			if eachDict['score'] >=  maxScore*0.7 and eachDict['score']*0.7 > 100:
				if eachDict['rxns'].intersection(rxnsInPathsNaive):
					categoriesCount['fromNaivePaths'] += 1
					filtered1.append(eachDict)
				
	if not filtered1:
		for eachDict in top50Simplified:
			if eachDict['score'] >=  maxScore*0.8 and eachDict['score']*0.8 > 100:
				if eachDict['rxns'].intersection(rxnsInPathsFromSoft):
					categoriesCount['fromSoftConstraints'] += 1
					filtered1.append(eachDict)
	
	if not filtered1:
		#give priority to genes selected on annotatetion file
		for eachDict in top50Simplified:
			if eachDict['source_gene'] in geneAndProteinNamePerSeqId and geneAndProteinNamePerSeqId[eachDict['source_gene']]['gene'] and geneAndProteinNamePerSeqId[eachDict['source_gene']]['protein name 1'] and (geneAndProteinNamePerSeqId[eachDict['source_gene']]['gene'].lower() == swiss90tremble50SeqInfo[eachDict['uniprotEntry']]['gene'].lower() or swiss90tremble50SeqInfo[eachDict['uniprotEntry']]['title'].lower() in geneAndProteinNamePerSeqId[eachDict['source_gene']]['protein name 1'].lower()):
				categoriesCount['fromAnotation'] += 1
				filtered1.append(eachDict)
			elif eachDict['source_gene'] in geneAndProteinNamePerSeqId and geneAndProteinNamePerSeqId[eachDict['source_gene']]['gene'] and geneAndProteinNamePerSeqId[eachDict['source_gene']]['protein name 2'] and (swiss90tremble50SeqInfo[eachDict['uniprotEntry']]['title'].lower() in geneAndProteinNamePerSeqId[eachDict['source_gene']]['protein name 2'].lower()):
				categoriesCount['fromAnotation'] += 1
				filtered1.append(eachDict)
		if filtered1 == top50Simplified: filtered1 = list()
		
		
		#give priority to swissProtIds
		if not filtered1:
			for eachDict in top50Simplified:
				if eachDict['uniprotEntry'] in swissProtIds:
					filtered1.append(eachDict)
					categoriesCount['fromSwissProtIds'] += 1
				
	if not filtered1: filtered1 = top50Simplified
	
	return filtered1


def findBestPerRead(top60, swissProtIds, rxnsInTaxonomyConstraints, rxnsFromUniprot, rxnsInPathsFromSoft, rxnsInPathsNaive, subunits, swiss90tremble50SeqInfo, geneAndProteinNamePerSeqId, categoriesCount):
		
	#find me maximum score to cut out set of reactions < 0.5*maxScore
	maxScore = max([i['score'] for i in top60])
			
	top50Simplified = list()
	for eachDict in top60:
		if eachDict['score'] >= maxScore*0.5:
			top50Simplified.append(eachDict)
		
		
	#da prioridade ao match com maior intersecao com reacoes anotadas pelo biocyc/kegg, quando ha dois matchs, e os dois tem scores proximos
	filteredByPriority = findPrioritary(top50Simplified, categoriesCount, rxnsInTaxonomyConstraints, rxnsFromUniprot, rxnsInPathsFromSoft, rxnsInPathsNaive, subunits, geneAndProteinNamePerSeqId, swiss90tremble50SeqInfo, swissProtIds, firstloop = 1)
	while True:
		filteredByPriority2 = findPrioritary(filteredByPriority, categoriesCount, rxnsInTaxonomyConstraints, rxnsFromUniprot, rxnsInPathsFromSoft, rxnsInPathsNaive, subunits, geneAndProteinNamePerSeqId, swiss90tremble50SeqInfo, swissProtIds, firstloop = 0)
		if filteredByPriority2 == filteredByPriority: break
		filteredByPriority = filteredByPriority2
	
	maximumScore = 0
	maximumScoreIndex = 0
	for count, eachList in enumerate(filteredByPriority):
		if eachList['score'] > maximumScore:
			maximumScore = eachList['score']
			maximumScoreIndex = count

	#separate in two lists without intersection
	bestMatch = top50Simplified[maximumScoreIndex]
	del top50Simplified[maximumScoreIndex]
	
	notBestMatch = list()
	for eachDict in top50Simplified:
		if not eachDict['rxns'].issubset(bestMatch['rxns']):
			notBestMatch.append(eachDict)
			
	return bestMatch, notBestMatch
		
def changeScoreOfSoft(synId, rxnsScores, reframedModel, dbToCobraId, score):

	if synId not in dbToCobraId: return
	for rxnCobraId in dbToCobraId[synId]:
		rxnReframedId = 'R_'+rxnCobraId.replace('-','__45__').replace('.','__46__').replace('+','__43__')
		if rxnReframedId in reframedModel.reactions:	
			if rxnReframedId not in rxnsScores: 
				rxnsScores[rxnReframedId] = score
			else:
				if rxnsScores[rxnReframedId] < score: rxnsScores[rxnReframedId] = score
				elif rxnsScores[rxnReframedId] < 0 and rxnsScores[rxnReframedId] < -0.2: rxnsScores[rxnReframedId] += 0.1
	
	
def findPathways(metacycRxnsInModel, pathwaysToInclude=set()):

	pickle_file_path = os.path.join(project_dir, 'data/generated', 'biocycPathways.pickle')
	with open(pickle_file_path, 'rb') as f:
		biocycPathways = pickle.load(f)
		
	pickle_file_path = os.path.join(project_dir, 'data/generated', 'rxnsPerModules.pickle')
	with open(pickle_file_path, 'rb') as f:
		rxnsPerModules = pickle.load(f)

	paraEscrever = ""
	contador = 1
	for eachRxn in metacycRxnsInModel:
		paraEscrever += "read" + str(contador) + "	" + eachRxn + "\n"
		contador += 1
	paraEscrever = paraEscrever[:-1]
	humannRxns = open('biocycRxnsInModel.tsv', 'w')
	_ = humannRxns.write(paraEscrever)
	humannRxns.close()
	
	#creating file to be used as a mapping by the MinPath
	superPath = set()
	toWrite = "#Metacyc pathway and reactions mapping file\n#Pathway	ReactionID"
	for eachPath in biocycPathways:
		if pathwaysToInclude and eachPath not in pathwaysToInclude: continue
		if 'superpathway' in biocycPathways[eachPath]['name'].lower(): 
			superPath.add(eachPath)
			continue
		if len(biocycPathways[eachPath]['RxnsInvolved']) < 3: continue
		for eachRxn in biocycPathways[eachPath]['RxnsInvolved']:
			toWrite += "\n" + eachPath + '	' + eachRxn
	
	for eachPath in rxnsPerModules:
		if pathwaysToInclude and eachPath not in pathwaysToInclude: continue
		for eachRxn in rxnsPerModules[eachPath]:
			toWrite += "\n" + eachPath + '	' + eachRxn
	
	mapping = open('metcycRxnPathwayMap.tsv', 'w')
	_ = mapping.write(toWrite)
	mapping.close()
	
	#find the minimum set of pathways thta explains the reactions
	MinPathResult1 = MinPathMain(anyfile = 'biocycRxnsInModel.tsv', mapfile = 'metcycRxnPathwayMap.tsv')
	
	#find the minimum set of pathways thta explains the reactions, including superpathways
	toWrite = "#Metacyc pathway and reactions mapping file\n#Pathway	ReactionID"
	for eachPath in biocycPathways:
		if pathwaysToInclude and eachPath not in pathwaysToInclude: continue
		if eachPath in superPath or eachPath in MinPathResult1:
			for eachRxn in biocycPathways[eachPath]['RxnsInvolved']:
				toWrite += "\n" + eachPath + '	' + eachRxn
	
	mapping = open('metcycRxnPathwayMap.tsv', 'w')
	_ = mapping.write(toWrite)
	mapping.close()
	
	MinPathResult2 = MinPathMain(anyfile = 'biocycRxnsInModel.tsv', mapfile = 'metcycRxnPathwayMap.tsv')
	
	minPathResult = MinPathResult1|MinPathResult2
	
	return minPathResult

def rheaToIdInModel(cobraModel):

	dbToCobraId = dict()
	for cobraRxn in cobraModel.reactions: 
		if 'rhea' in cobraRxn.annotation:
			for rheaId in cobraRxn.annotation['rhea']:
				if rheaId not in dbToCobraId: dbToCobraId[int(rheaId)] = list()
				dbToCobraId[int(rheaId)].append(cobraRxn.id)
		else:		
			if 'kegg' in cobraRxn.annotation:
				for keggId in cobraRxn.annotation['kegg']:
					if keggId not in dbToCobraId: dbToCobraId[keggId] = list()
					dbToCobraId[keggId].append(cobraRxn.id)
			if 'metacyc' in cobraRxn.annotation:
				for metacycId in cobraRxn.annotation['metacyc']:
					if metacycId not in dbToCobraId: dbToCobraId[metacycId] = list()
					dbToCobraId[metacycId].append(cobraRxn.id)
			
	return dbToCobraId


	

def reaction_scoring(diamondResult, geneAndProteinNamePerSeqId, cobraModel, reframedModel, constraintsFromTaxonomy, constraintsFromFile, rxnsInTaxonomyConstraints, rxnsFromUniprot, outputfolder, verbose):
	
	

	""" Calculate reaction scores using new eggnog output.
	Args:
		diamondResult (pandas.DataFrame): gene diamondResult results
		gprs (pandas.DataFrame): BiGG GPR rules
		spontaneous_score (float): score to give to spontaneous reactions (default: 0.0)

	Returns:
		pandas.DataFrame: reaction scores
	"""

	pickle_file_path = os.path.join(project_dir, 'data/generated', 'biocycPathways.pickle')
	with open(pickle_file_path, 'rb') as f:
		biocycPathways = pickle.load(f)

	pickle_file_path = os.path.join(project_dir, 'data/generated', 'swissProtIds.pickle')
	print(pickle_file_path)
	with open(pickle_file_path, 'rb') as f:
		swissProtIds = pickle.load(f)
		
	pickle_file_path = os.path.join(project_dir, 'data/generated', 'swiss90tremble50SeqInfo.pickle')
	with open(pickle_file_path, 'rb') as f:
		swiss90tremble50SeqInfo = pickle.load(f)
		
	pickle_file_path = os.path.join(project_dir, 'data/generated', 'rxnsPerModules.pickle')
	with open(pickle_file_path, 'rb') as f:
		rxnsPerModules = pickle.load(f)
		
	
	saveProgressFile(34, outputfolder)
	
	
	if verbose: print('\nStarting reactions scoring ' + str(datetime.now()) + '\n')


	#cria lista de traducoes, do rhea para kegg e metacyc
	metacycToRheaDict = dict()
	metacycToModelIdDict = dict()
	modelIdToDbDict = dict()
	rheaToMetacycDict = dict()
	rheaIdToReframedId = dict()
	for rxn in cobraModel.reactions:
	
		idInReframed = 'R_' + rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
		if 'rhea' not in rxn.annotation: 
			
			for db in ['kegg', 'metacyc']:
				
				if db not in rxn.annotation: continue
			
				for rxnIdInDb in rxn.annotation[db]:
					if idInReframed not in rheaToMetacycDict: rheaToMetacycDict[idInReframed] = set()
					rheaToMetacycDict[idInReframed].add(rxnIdInDb)
					
					if rxnIdInDb not in metacycToModelIdDict: metacycToModelIdDict[rxnIdInDb] = set()
					metacycToModelIdDict[rxnIdInDb].add(idInReframed)

					if idInReframed not in modelIdToDbDict: modelIdToDbDict[idInReframed] = set()
					modelIdToDbDict[idInReframed].add(rxnIdInDb)
		
		else:
			for rheaId in rxn.annotation['rhea']:
			
				if rheaId not in rheaIdToReframedId: rheaIdToReframedId[int(rheaId)] = set()
				rheaIdToReframedId[int(rheaId)].add(idInReframed)
				
				for db in ['kegg', 'metacyc']:
					
					if db not in rxn.annotation: continue
				
					for rxnIdInDb in rxn.annotation[db]:
						if rxnIdInDb not in metacycToRheaDict: metacycToRheaDict[rxnIdInDb] = set()
						if rxnIdInDb not in metacycToModelIdDict: metacycToModelIdDict[rxnIdInDb] = set()
						if rheaId not in rheaToMetacycDict: rheaToMetacycDict[rheaId] = set()
						if idInReframed not in modelIdToDbDict: modelIdToDbDict[idInReframed] = set()
						metacycToRheaDict[rxnIdInDb].add(int(rheaId))
						metacycToModelIdDict[rxnIdInDb].add(idInReframed)
						rheaToMetacycDict[rheaId].add(rxnIdInDb)
						modelIdToDbDict[idInReframed].add(rxnIdInDb)
					
	soft_constraints_pathways = constraintsFromTaxonomy['metacyc']['pathways']|constraintsFromTaxonomy['kegg']['pathways']
	
	
	rxnsInPathsFromSoft = set()
	for eachPath in soft_constraints_pathways:
		if eachPath in rxnsPerModules:
			for eachRxn in rxnsPerModules[eachPath]:
				if eachRxn not in metacycToRheaDict: continue
				for rheaRxn in metacycToRheaDict[eachRxn]: rxnsInPathsFromSoft.add(rheaRxn)
		if eachPath in biocycPathways:
			for eachRxn in biocycPathways[eachPath]['RxnsInvolved']:
				if eachRxn not in metacycToRheaDict: continue
				for rheaRxn in metacycToRheaDict[eachRxn]: rxnsInPathsFromSoft.add(rheaRxn)
			
	#find reactions annotated by uniprot for the specie
	dbToCobraId = rheaToIdInModel(cobraModel)
			
			
	rxnsInTaxonomyConstraints = rxnsInTaxonomyConstraints
	rxnsInPathsNaive = set() # initializing variable.
			

	subunits = list()
	
	diamondResultDict = diamondResult.to_dict('records') # dict is faster than pandas in this application
	
	categoriesCount = {'subunits': 0, 'fromMinPath': 0, 'fromSoftConstraints': 0, 'fromNaivePaths':0, 'fromAnotation': 0, 'fromSwissProtIds':0, 'rxnsFromUniprot': 0}

	for repetition in [1,2]: #repeats the scoring two times. In the second time, favors reactions missing on almost complet pathways, and subunits of protein complexes. 


		#select one target_gene per "source_gene"
		currentRead = ''
		bestPerRead = list()
		notBestPerReadOdered = list()
		for eachMatch in diamondResultDict:
		
			if not swiss90tremble50SeqInfo[eachMatch['target_gene']]['gene']: continue
		
			if eachMatch['source_gene'] != currentRead: 
				#quando comeca o alinhamento de uma nova sequencia de proteina
				
				if currentRead and top60:
				
					notBestPerReadOdered.extend(top60)
					
					bestMatch, notBestMatch = findBestPerRead(top60, swissProtIds, rxnsInTaxonomyConstraints, rxnsFromUniprot, rxnsInPathsFromSoft, rxnsInPathsNaive, subunits, swiss90tremble50SeqInfo, geneAndProteinNamePerSeqId, categoriesCount) # encontra o melhor match entre os os resultados mais frequentes.
					if bestMatch: bestPerRead.append(bestMatch)
					for match in notBestMatch: notBestPerReadOdered.append(match)
					
				currentRead = eachMatch['source_gene']
				top60 = list()
			
			#penalize proteins with weak evidance of existence. level 1: protein evidence; level 2: transcript evidence; level 3: homolugos; level 4: predict; level 5: uncertain
			if swiss90tremble50SeqInfo[eachMatch['target_gene']]['evidence level'] == 1: newSocre = eachMatch['score']
			elif swiss90tremble50SeqInfo[eachMatch['target_gene']]['evidence level'] == 2: newSocre = eachMatch['score']/1.1
			elif swiss90tremble50SeqInfo[eachMatch['target_gene']]['evidence level'] == 3: newSocre = eachMatch['score']/1.2
			elif swiss90tremble50SeqInfo[eachMatch['target_gene']]['evidence level'] == 4: newSocre = eachMatch['score']/2
			elif swiss90tremble50SeqInfo[eachMatch['target_gene']]['evidence level'] == 5: newSocre = eachMatch['score']/4
			else: 
				print(swiss90tremble50SeqInfo[eachMatch['target_gene']]['evidence level'])
				newSocre = eachMatch['score']
					
			rxnsInGene = {i for i in swiss90tremble50SeqInfo[eachMatch['target_gene']]['rxns'] if i in dbToCobraId}
			if not rxnsInGene: continue
						
			if not top60:
				top60.append({'rxns':rxnsInGene,'count':1, 'uniprotEntry': eachMatch['target_gene'], 'source_gene': eachMatch['source_gene'], 'score':newSocre, 'gene': swiss90tremble50SeqInfo[eachMatch['target_gene']]['gene']})
			else:
				sucesso = 0
				for eachDict in top60:
					if eachDict['rxns'] == rxnsInGene:
						sucesso = 1
						eachDict['count'] += 1
						if eachDict['score'] < newSocre:
							eachDict['score'] = newSocre
							eachDict['gene'] = swiss90tremble50SeqInfo[eachMatch['target_gene']]['gene']
							eachDict['uniprotEntry'] = eachMatch['target_gene']
						break
				if sucesso == 0:
					top60.append({'rxns':rxnsInGene, 'count':1, 'uniprotEntry': eachMatch['target_gene'], 'source_gene': eachMatch['source_gene'], 'score':newSocre, 'gene': swiss90tremble50SeqInfo[eachMatch['target_gene']]['gene']})
			
		if top60:
			bestMatch, notBestMatch = findBestPerRead(top60, swissProtIds, rxnsInTaxonomyConstraints, rxnsFromUniprot, rxnsInPathsFromSoft, rxnsInPathsNaive, subunits, swiss90tremble50SeqInfo, geneAndProteinNamePerSeqId, categoriesCount)
			if bestMatch: bestPerRead.append(bestMatch)
			for match in notBestMatch: notBestPerReadOdered.append(match)


		#sort bestPerRead, from high score to low score
		scores = list() 
		for read in bestPerRead: scores.append(read['score'])
		scores.sort(reverse=True)
		
		bestPerReadSorted = list()
		for score in scores:
			for read in bestPerRead:
				if read['score'] == score:
					break
			bestPerReadSorted.append(read)
			bestPerRead.remove(read)
			
		#simplify bestPerReadSorted, by removing duplicates
		bestMatchPerRead = list()
		for read1 in bestPerReadSorted:
			sucesso = 0
			for read2 in bestMatchPerRead:
				if read2['rxns'] == read1['rxns']:
					sucesso = 1
					break
			if sucesso == 0:
				bestMatchPerRead.append(read1)
			
		#find pathways in organism. pathwaysPresent will be used outside this loop
		#encontra as reacoes 
		metacycRxnsInModel = set()
		for eachDict in bestMatchPerRead:
			for rheaId in eachDict['rxns']:
				rheaId = str(rheaId)
				if rheaId not in rheaToMetacycDict: 
					continue
				for rxnSyn in rheaToMetacycDict[rheaId]:
					metacycRxnsInModel.add(rxnSyn)
				
		minPathResult = findPathways(metacycRxnsInModel)
		
		#find the pathwas in constraintsFromTaxonomy
		if len(soft_constraints_pathways) > 100: pathwaysPresent = minPathResult.intersection(soft_constraints_pathways)
		else: pathwaysPresent = minPathResult
			
			
		if repetition == 1:
		
			#find proteins that are subunits of protein complexes
			subunits = list()
			uniprotSubunitEntries = set()
			for eachDict in bestMatchPerRead:
				if eachDict['uniprotEntry'] in uniprotSubunitEntries: continue
				subunitType = ''
				for subunitString in ['subunit', 'component']:
					if subunitString in swiss90tremble50SeqInfo[eachDict['uniprotEntry']]['title']:
						#to extract the subunit part
						if subunitString in swiss90tremble50SeqInfo[eachDict['uniprotEntry']]['title']:
						
							titleProvisional = swiss90tremble50SeqInfo[eachDict['uniprotEntry']]['title']
							if titleProvisional.count(subunitString) > 1: 
								titleProvisional = titleProvisional.replace(subunitString,'',1)
						
							if '-'+subunitString in titleProvisional.lower():
								for word in titleProvisional.lower().split(' '):
									if '-'+subunitString in word: break
								subunitType = word.replace('-'+subunitString,'')
							else:
							
								if '/'+subunitString in titleProvisional.lower():
									toReplace = titleProvisional.lower().split('/'+subunitString)[0]
									titleProvisional = titleProvisional.replace(toReplace.split(' ')[-1]+'/','')
								
								if subunitString+' ' in titleProvisional.lower():
									subunitType = titleProvisional.lower().split(subunitString+' ')[1]
									subunitType = subunitType.split(' ')[0]
								elif ' '+subunitString in titleProvisional.lower():
									subunitType = titleProvisional.lower().split(' '+subunitString)[0]
									subunitType = subunitType.split(' ')[-1]
				if subunitType:
					uniprotSubunitEntries.add(eachDict['uniprotEntry'])
					subunits.append({'rxns': eachDict['rxns'] , 'subunit': subunitType, 'gene':eachDict['gene'], 'uniprotEntry':eachDict['uniprotEntry']})
		
		
		#looking for rxns missing pathways present in the model
		rxnsPresentInPathways = set()
		for pathway in pathwaysPresent:
			if pathway in biocycPathways:
				if len(biocycPathways[pathway]['RxnsInvolved'].intersection(metacycRxnsInModel)) <= 3: continue
				rxnsPresentInPathways = rxnsPresentInPathways | biocycPathways[pathway]['RxnsInvolved']
			if pathway in rxnsPerModules:
				if len(rxnsPerModules[pathway].intersection(metacycRxnsInModel)) <= 3: continue
				rxnsPresentInPathways = rxnsPresentInPathways | rxnsPerModules[pathway]
		
		for metacycRxn in rxnsPresentInPathways:
			if metacycRxn not in metacycToRheaDict: rxnsInTaxonomyConstraints.add(metacycRxn)# continue
			else:
				for rheaRxn in metacycToRheaDict[metacycRxn]: rxnsInTaxonomyConstraints.add(rheaRxn)
				
				
		#find all pathways associated to all reactions matched
		rxnsInMappedPaths = set()
		for eachPath in biocycPathways:
			if eachPath in biocycPathways:
				if biocycPathways[eachPath]['RxnsInvolved'].intersection(metacycRxnsInModel):
					for rxnId in biocycPathways[eachPath]['RxnsInvolved']: rxnsInMappedPaths.add(rxnId)
			if eachPath in rxnsPerModules:
				if rxnsPerModules[eachPath].intersection(metacycRxnsInModel):
					for rxnId in rxnsPerModules[eachPath]: rxnsInMappedPaths.add(rxnId)
					rxnsPresentInPathways = rxnsPresentInPathways | rxnsPerModules[eachPath]
		
		for metacycRxn in rxnsInMappedPaths:
			if metacycRxn not in metacycToRheaDict: rxnsInPathsNaive.add(metacycRxn)# continue
			else:
				for rheaRxn in metacycToRheaDict[metacycRxn]: rxnsInPathsNaive.add(rheaRxn)
				
		saveProgressFile(36, outputfolder)
	saveProgressFile(38, outputfolder)



	#replace rhea IDs by IDs in reframed model
	rxnsInBest = set()
	for eachList in bestMatchPerRead:
		rxnsSet = set()
		for rxnIdSyn in eachList['rxns']:
			if rxnIdSyn not in rheaIdToReframedId: continue
			for rxnIdInModel in rheaIdToReframedId[rxnIdSyn]:
				rxnsSet.add(rxnIdInModel)
				rxnsInBest.add(rxnIdInModel)
		eachList['rxns'] = rxnsSet
	

	notBestPerRead = list()
	rxnsInNotBest = set()
	for eachList in notBestPerReadOdered:
		
		rxnsSet = set()
		for rxnIdSyn in eachList['rxns']:
			if rxnIdSyn in rxnsInBest: continue
			if rxnIdSyn not in rheaIdToReframedId: continue
			for rxnIdInModel in rheaIdToReframedId[rxnIdSyn]:
				rxnsSet.add(rxnIdInModel)
				rxnsInNotBest.add(rxnIdInModel)
				
		if rxnsSet: 
			eachList['rxns'] = rxnsSet
			notBestPerRead.append(eachList)
	

	#to retrive association between rxns and genes
	rheaIdToGene = {'bestMatchPerRead':bestMatchPerRead, 'notBestPerRead':notBestPerRead}


	#comeca a trabalhar apenas com reacoes presentes no modelo universal
	medianValue = median([eachDict['score'] for eachDict in bestMatchPerRead])


	make100Zero = 100/medianValue #100 might lead to wrong results. increased to 150. Maybe 200 to Trembl and 100 to SwissProt?
	make200Zero = 200/medianValue #100 might lead to wrong results. increased to 150. Maybe 200 to Trembl and 100 to SwissProt?

	#normalize
	for eachDict in bestMatchPerRead:
		if eachDict['uniprotEntry'] in swissProtIds: eachDict['scoreNormalized'] = round((eachDict['score']/medianValue) - make100Zero, 4)
		else: eachDict['scoreNormalized'] = round((eachDict['score']/medianValue) - make200Zero, 4)
		
		if eachDict['scoreNormalized'] > 0: eachDict['scoreNormalized'] = eachDict['scoreNormalized']*5
		elif eachDict['scoreNormalized'] < 0: eachDict['scoreNormalized'] = eachDict['scoreNormalized']/5
		elif eachDict['scoreNormalized'] == 0: eachDict['scoreNormalized'] = 0.01
	
	
	#finds subunits 
	moreThanTwoSubunits = set()
	for eachDictIndex1 in range(len(subunits)-1):
		for eachDictIndex2 in range(eachDictIndex1+1, len(subunits)):
			if not subunits[eachDictIndex1]['gene'] or not subunits[eachDictIndex2]['gene']: continue
			if subunits[eachDictIndex1]['rxns'].intersection(subunits[eachDictIndex2]['rxns']) and subunits[eachDictIndex1]['subunit'] != subunits[eachDictIndex2]['subunit'] and subunits[eachDictIndex1]['gene'] != subunits[eachDictIndex2]['gene']:
				moreThanTwoSubunits.add(subunits[eachDictIndex1]['uniprotEntry'])
				moreThanTwoSubunits.add(subunits[eachDictIndex2]['uniprotEntry'])
	

	#finding rhea without a complex
	rxnsScores = dict() #dict to store score of rhea reactions
	for eachDict in bestMatchPerRead:
		if eachDict['uniprotEntry'] in uniprotSubunitEntries:
			for rheaId in eachDict['rxns']:
				if eachDict['uniprotEntry'] in moreThanTwoSubunits:
					if rheaId not in rxnsScores: rxnsScores[rheaId] = eachDict['scoreNormalized']
					elif rxnsScores[rheaId] < eachDict['scoreNormalized']: rxnsScores[rheaId] = eachDict['scoreNormalized']
				else:
					if rheaId not in rxnsScores: 
						if eachDict['scoreNormalized'] < 0: rxnsScores[rheaId] = eachDict['scoreNormalized']
						else: rxnsScores[rheaId] = 0.01
		else:
			for rheaId in eachDict['rxns']:
				if rheaId not in rxnsScores: rxnsScores[rheaId] = eachDict['scoreNormalized']
				elif rxnsScores[rheaId] < eachDict['scoreNormalized']: rxnsScores[rheaId] = eachDict['scoreNormalized']
	
	
	#cut off very high values
	rxnBelonginToPaths = rxnsInPathsNaive|rxnsInTaxonomyConstraints
	for rheaId in rxnsScores:
		if rxnsScores[rheaId] > 0 and rheaId in rxnBelonginToPaths: rxnsScores[rheaId] = rxnsScores[rheaId]*1.5
		
		if rxnsScores[rheaId] > 10: rxnsScores[rheaId] = 10
		if rxnsScores[rheaId] > 0.1 and rheaId not in rxnBelonginToPaths: rxnsScores[rheaId] = rxnsScores[rheaId]/1.5
		

	minScoreBase = min(rxnsScores.values())
	
	if minScoreBase > 0: minScore = -0.1
	else: minScore = minScoreBase - 0.1
	
	#reactions not associated with the best match of each protein sequence
	for rheaId in rxnsInNotBest - rxnsInBest:
		rxnsScores[rheaId] = minScoreBase -0.5
	

	if constraintsFromTaxonomy:
		for rxnId in constraintsFromTaxonomy['kegg']['1 intersection']:
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, dbToCobraId, 1)
	
		for rxnId in constraintsFromTaxonomy['kegg']['1']:
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, dbToCobraId, -0.001)
		for rxnId in constraintsFromTaxonomy['kegg']['0.9 to 0.99']:
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, dbToCobraId, minScore - 0.1)
		for rxnId in constraintsFromTaxonomy['kegg']['0.8 to 0.89']:
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, dbToCobraId, minScore - 0.3)
		for rxnId in constraintsFromTaxonomy['kegg']['0.7 to 0.79']:
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, dbToCobraId, minScore - 0.5)
		for rxnId in constraintsFromTaxonomy['kegg']['0.6 to 0.69']:
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, dbToCobraId, minScore - 0.8)

		for rxnId in constraintsFromTaxonomy['metacyc']['1']:
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, dbToCobraId, -0.01)
		for rxnId in constraintsFromTaxonomy['metacyc']['0.9 to 0.99']:
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, dbToCobraId, minScore - 0.3)
		for rxnId in constraintsFromTaxonomy['metacyc']['0.8 to 0.89']:
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, dbToCobraId, minScore - 0.6)
		for rxnId in constraintsFromTaxonomy['metacyc']['0.7 to 0.79']:
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, dbToCobraId, minScore - 0.85)	
	
	
	##include soft constraints in objective
	for cobraId in constraintsFromFile['soft']:
		rxnId = 'R_'+cobraId.replace('-','__45__').replace('.','__46__').replace('+','__43__')
		
		if rxnId in rxnsScores:
			if constraintsFromFile['soft'][cobraId] > -1: #if greater than -1, means the user wants to allow the reaction to be present. Could be to force the reaction to be present (greater than 0), or to use the reaction as first option on gapfill (between -1 and 0)
				if constraintsFromFile['soft'][cobraId] > rxnsScores[rxnId]:
					rxnsScores[rxnId] = constraintsFromFile['soft'][cobraId]
			else: #if lower than -1, the user does not want the reaction, even if there is ma match in uniprot 
				rxnsScores[rxnId] = constraintsFromFile['soft'][cobraId]
		else:
			rxnsScores[rxnId] = constraintsFromFile['soft'][cobraId]
			
	for cobraId in constraintsFromFile['hard']:
		rxnId = 'R_'+cobraId.replace('-','__45__').replace('.','__46__').replace('+','__43__')
		rxnsScores[rxnId] = (constraintsFromFile['hard'][cobraId]/abs(constraintsFromFile['hard'][cobraId]))*100 # rxnsScores[rxnId] = 100, keeping the sign of the original constraint
		
	#include transport reactions related soft_constraints
	for constraintsType in ['soft', 'hard']:
		for cobraId in constraintsFromFile[constraintsType]:
		
			if not cobraId.startswith('EX_'): continue
		
			rxnOriginal = cobraModel.reactions.get_by_id(cobraId)
			for met in rxnOriginal.metabolites: break

			transport = set()
			exchange = set()
			other = set()
			
			for rxn in met.reactions:
				
				if rxnOriginal == rxn: continue
				
				if len(rxn.compartments) == 2: transport.add(rxn)
				elif (rxn.products and not rxn.reactants) or (not rxn.reactants and rxn.products): exchange.add(rxn)
				else: other.add(rxn)
				
			if not transport and other:
				for rxnOther in other:
					for met in rxnOther.metabolites:
						for rxn in met.reactions:
							if len(rxn.compartments) == 2: transport.add(rxn)
							
			for rxn in transport:
				rxnId = 'R_'+rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
				if rxnId in rxnsScores: continue
				if len(rxn.metabolites) > 2: rxnsScores[rxnId] = -0.01
				else: rxnsScores[rxnId] = -0.02
	
		
	#include reactions annotated by uniprot in the species	
	for rheaId in rxnsFromUniprot:
		if rheaId not in rheaIdToReframedId: continue
		for rxnId in rheaIdToReframedId[rheaId]:
			if rxnId in rxnsScores:
				if rxnsScores[rxnId] < minScore: rxnsScores[rxnId] = minScore
			else:
				rxnsScores[rxnId] = minScore

	#find the reactions missing in pathways that are almost complete. These reactions will receive the minimum score
	metacycPositiveScore = set()
	for modelId in rxnsInBest|rxnsInNotBest:
		if modelId not in modelIdToDbDict: continue
		for rxnSyn in modelIdToDbDict[modelId]: metacycPositiveScore.add(rxnSyn)
	
	
	#check completness of all big biocyc pathway
	rxnsMissingInPathways = set()
	for pathwayBiocyc in set(pathwaysPresent).intersection(biocycPathways):
		emModelo = biocycPathways[pathwayBiocyc]['RxnsInvolved'].intersection(metacycPositiveScore)
		if len(emModelo) / len(biocycPathways[pathwayBiocyc]['RxnsInvolved']) >= 0.5:
			rxnsMissingInPathways = rxnsMissingInPathways | (biocycPathways[pathwayBiocyc]['RxnsInvolved'] - metacycPositiveScore)
				
	#check completness of all big kegg pathway
	for pathwayKegg in set(pathwaysPresent).intersection(rxnsPerModules):
		emModelo = rxnsPerModules[pathwayKegg].intersection(metacycPositiveScore)
		if len(emModelo) / len(rxnsPerModules[pathwayKegg]) >= 0.5:
			rxnsMissingInPathways = rxnsMissingInPathways | (rxnsPerModules[pathwayKegg] - metacycPositiveScore)
	
	rxnsRheaFaltandoPathways = set()
	for rxnMetacycId in rxnsMissingInPathways:
		if rxnMetacycId in metacycToModelIdDict:
			for modelId in metacycToModelIdDict[rxnMetacycId]:
				rxnsRheaFaltandoPathways.add(modelId)
	
	minScore = min(rxnsScores.values())
	for rheaId in rxnsRheaFaltandoPathways:
		if rheaId in rxnsScores: 
			if rxnsScores[rheaId] < 0: rxnsScores[rheaId]=rxnsScores[rheaId]/2
		else: rxnsScores[rheaId] = minScore
		
	

	#find the reactions that are always alone in the gene-reactions association
	rxnsSimultaneous = dict()
	for eachDict in bestMatchPerRead:
		rxnsSet = set()
		for rxnId in eachDict['rxns']:
			rxnsSet.add(rxnId.replace('_forwardTemp','').replace('_reverseTemp',''))
		for rxnId in rxnsSet:
			if rxnId not in rxnsSimultaneous or len(rxnsSet) > rxnsSimultaneous[rxnId]: rxnsSimultaneous[rxnId] = len(rxnsSet)
			
	singleRxnsInGenes = set()
	for rxnId in rxnsSimultaneous:
		if rxnsSimultaneous[rxnId] != 1: continue
		if rxnId in reframedModel.reactions: singleRxnsInGenes.add(rxnId)
		if rxnId + '_forwardTemp' in reframedModel.reactions: singleRxnsInGenes.add(rxnId + '_forwardTemp')
		if rxnId + '_reverseTemp' in reframedModel.reactions: singleRxnsInGenes.add(rxnId + '_reverseTemp')
		

	return rxnsScores, rheaIdToGene, bestMatchPerRead, singleRxnsInGenes
	
