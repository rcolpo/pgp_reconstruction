import pickle
from pgp_reconstruction import project_dir
from statistics import median
import sys
import os
import datetime

try:from pgp_reconstruction.dependencies.MinPath_master.MinPath import MinPathMain
except: pass

def findPrioritary(top50Simplified, categoriesCount, rxnsInPathsFromMinPath, rxnsInPathsFromSoft, rxnsInPathsNaive, subunits, geneAndProteinNamePerSeqId, sprotSeqAndGenDict, swissProtIds, firstloop = 1):
	#first check if there is any match from rxnsInPathsFromMinPath
	filtered1 = list()
	
	#try to find the missing subunits
	if not filtered1 and firstloop == 1 and subunits:

		for eachDict in top50Simplified:
			if 'subunit' not in sprotSeqAndGenDict[eachDict['uniprotEntry']]['title'] and 'component' not in sprotSeqAndGenDict[eachDict['uniprotEntry']]['title']:
				sucesso = 0
				for eachDict in top50Simplified:
					subunitType = ''
					for subunitString in ['subunit', 'component']:
						if subunitString in sprotSeqAndGenDict[eachDict['uniprotEntry']]['title']:
							#to extract the subunit part
							if subunitString in sprotSeqAndGenDict[eachDict['uniprotEntry']]['title']:
							
								titleProvisional = sprotSeqAndGenDict[eachDict['uniprotEntry']]['title']
							
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
			if eachDict['rxns'].intersection(rxnsInPathsFromMinPath):
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
				
	if not filtered1 and firstloop == 1:
		for eachDict in top50Simplified:
			if eachDict['score'] >=  maxScore*0.8 and eachDict['score']*0.8 > 100:
				if eachDict['rxns'].intersection(rxnsInPathsFromSoft):
					categoriesCount['fromSoftConstraints'] += 1
					filtered1.append(eachDict)
	
	if not filtered1:
		#give priority to genes selected on annotatetion file
		for eachDict in top50Simplified:
			if eachDict['source_gene'] in geneAndProteinNamePerSeqId and geneAndProteinNamePerSeqId[eachDict['source_gene']]['gene'] and geneAndProteinNamePerSeqId[eachDict['source_gene']]['protein name 1'] and (geneAndProteinNamePerSeqId[eachDict['source_gene']]['gene'].lower() == sprotSeqAndGenDict[eachDict['uniprotEntry']]['gene'].lower() or sprotSeqAndGenDict[eachDict['uniprotEntry']]['title'].lower() in geneAndProteinNamePerSeqId[eachDict['source_gene']]['protein name 1'].lower()):
				categoriesCount['fromAnotation'] += 1
				filtered1.append(eachDict)
			elif eachDict['source_gene'] in geneAndProteinNamePerSeqId and geneAndProteinNamePerSeqId[eachDict['source_gene']]['gene'] and geneAndProteinNamePerSeqId[eachDict['source_gene']]['protein name 2'] and (sprotSeqAndGenDict[eachDict['uniprotEntry']]['title'].lower() in geneAndProteinNamePerSeqId[eachDict['source_gene']]['protein name 2'].lower()):
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


def findBestPerRead(top60, swissProtIds, rxnsInPathsFromMinPath, rxnsInPathsFromSoft, rxnsInPathsNaive, subunits, sprotSeqAndGenDict, geneAndProteinNamePerSeqId, categoriesCount):
	#find best match and alternative matches of each read

	freqRxnDict = dict()
	for eachDict in top60:
		for rxn in eachDict['rxns']:
			if rxn not in freqRxnDict: freqRxnDict[rxn] = 0
			freqRxnDict[rxn] += eachDict['count']

	top60v2 = list()
	for eachDict in top60:
		for rxn in freqRxnDict:
			if rxn in eachDict['rxns'] and freqRxnDict[rxn] < 10 and rxn not in rxnsInPathsFromMinPath and rxn not in rxnsInPathsFromSoft and rxn not in rxnsInPathsNaive:
				eachDict['rxns'].remove(rxn)
		if eachDict['rxns']:
			top60v2.append(eachDict)

	if top60v2:
		freqDict = dict()
		for eachDict in top60v2:
			if eachDict['count'] not in freqDict:
				freqDict[eachDict['count']] = list()
			freqDict[eachDict['count']].append(eachDict['score'])
		
		for eachFreq in freqDict:
			freqDict[eachFreq].sort(reverse=True)
			
		freqList = list(freqDict.keys())
		freqList.sort(reverse=True)
		
		top60Sorted = list()
		freqSum = sum(freqList)
		sumCount = 0
		#remove less frequent reactions sets to avoid wrong diamondResult on Trembl
		for eachFreq in freqList:
			if sumCount > freqSum*0.95: break
			for eachScore in freqDict[eachFreq]:
				if sumCount > freqSum*0.95: break
				for eachDict in top60v2:					
					if eachDict['count'] == eachFreq and eachDict['score'] == eachScore:
						top60Sorted.append(eachDict)
						sumCount += eachFreq
						break
			
		#find me maximum score to cut out set of reactions < 0.5*maxScore
		maxScore = max([i['score'] for i in top60Sorted])
				
		top50Simplified = list()
		for eachDict in top60Sorted:
			if eachDict['score'] >= maxScore*0.5:
				top50Simplified.append(eachDict)
	else: return 0, 0
		
		
	#da prioridade ao match com maior intersecao com reacoes anotadas pelo biocyc/kegg, quando ha dois matchs, e os dois tem scores proximos
	filteredByPriority = findPrioritary(top50Simplified, categoriesCount, rxnsInPathsFromMinPath, rxnsInPathsFromSoft, rxnsInPathsNaive, subunits, geneAndProteinNamePerSeqId, sprotSeqAndGenDict, swissProtIds, firstloop = 1)
	while True:
		filteredByPriority2 = findPrioritary(filteredByPriority, categoriesCount, rxnsInPathsFromMinPath, rxnsInPathsFromSoft, rxnsInPathsNaive, subunits, geneAndProteinNamePerSeqId, sprotSeqAndGenDict, swissProtIds)
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
		
def changeScoreOfSoft(rxnId, rxnsScores, reframedModel, score, softPositiveNewRxn):
	idInModel1 = 'R_'+rxnId.replace('-','__45__').replace('.','__46__').replace('+','__43__')
	for possibleId in [idInModel1]:
		if possibleId in reframedModel.reactions:	
			if possibleId not in rxnsScores: 
				rxnsScores[possibleId] = score
				if score >= 0: softPositiveNewRxn.add(rxnId)
			else:
				if rxnsScores[possibleId] < score: rxnsScores[possibleId] = score
				else: rxnsScores[possibleId] += 0.2

def reaction_scoring(diamondResult, geneAndProteinNamePerSeqId, cobraModel, reframedModel, uniprotToRheaRxns, soft_constraints, rxnsInPathsFromMinPath, verbose):
	
	""" Calculate reaction scores using new eggnog output.
	Args:
		diamondResult (pandas.DataFrame): gene diamondResult results
		gprs (pandas.DataFrame): BiGG GPR rules
		spontaneous_score (float): score to give to spontaneous reactions (default: 0.0)

	Returns:
		pandas.DataFrame: reaction scores
	"""

	pickle_file_path = os.path.join(project_dir, 'data/generated', 'keggModules.pickle')
	with open(pickle_file_path, 'rb') as f:
		keggModules = pickle.load(f)

	pickle_file_path = os.path.join(project_dir, 'data/generated', 'swissProtIds.pickle')
	with open(pickle_file_path, 'rb') as f:
		swissProtIds = pickle.load(f)
		
	pickle_file_path = os.path.join(project_dir, 'data/generated', 'biocycPathways.pickle')
	with open(pickle_file_path, 'rb') as f:
		biocycPathways = pickle.load(f)
		
	pickle_file_path = os.path.join(project_dir, 'data/generated', 'sprotSeqAndGenDict.pickle')
	with open(pickle_file_path, 'rb') as f:
		sprotSeqAndGenDict = pickle.load(f)
		
	pickle_file_path = os.path.join(project_dir, 'data/generated', 'superPathsDict.pickle')
	with open(pickle_file_path, 'rb') as f:
		superPathsDict = pickle.load(f)
		
	pickle_file_path = os.path.join(project_dir, 'data/generated', 'rxnsPerModules.pickle')
	with open(pickle_file_path, 'rb') as f:
		rxnsPerModules = pickle.load(f)

	
	if verbose: print('\nStarting reactions scoring ' + str(datetime.datetime.now()) + '\n')

	#cria lista de traducoes, do rhea para kegg e metacyc
	metacycToRheaDict = dict()
	rheaToMetacycDict = dict()
	rheaIdToreframedId = dict()
	for rxn in cobraModel.reactions:
	
		idInreframed = 'R_' + rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
		if 'rhea' not in rxn.annotation: 
			
			for db in ['kegg', 'metacyc']:
				
				if db not in rxn.annotation: continue
			
				for rxnKeggId in rxn.annotation[db]:
					if idInreframed not in rheaToMetacycDict: rheaToMetacycDict[idInreframed] = set()
					rheaToMetacycDict[idInreframed].add(rxnKeggId)
					
					if rxnKeggId not in metacycToRheaDict: metacycToRheaDict[rxnKeggId] = set()
					metacycToRheaDict[rxnKeggId].add(idInreframed)
					#if 'rhea' not in rxn.annotation: continue
					#for rheaSyn in rxn.annotation['rhea']: metacycToRheaDict[rxnKeggId].add(rheaSyn)
		
		else:
			for rheaId in rxn.annotation['rhea']:
			
				if rheaId not in rheaIdToreframedId: rheaIdToreframedId[rheaId] = set()
				rheaIdToreframedId[rheaId].add(idInreframed)
				
				for db in ['kegg', 'metacyc']:
					
					if db not in rxn.annotation: continue
				
					for rxnKeggId in rxn.annotation[db]:
						if rxnKeggId not in metacycToRheaDict: metacycToRheaDict[rxnKeggId] = set()
						if rheaId not in rheaToMetacycDict: rheaToMetacycDict[rheaId] = set()
						if idInreframed not in rheaToMetacycDict: rheaToMetacycDict[idInreframed] = set()
						metacycToRheaDict[rxnKeggId].add(idInreframed)
						metacycToRheaDict[rxnKeggId].add(rheaId)
						rheaToMetacycDict[rheaId].add(rxnKeggId)
						rheaToMetacycDict[idInreframed].add(rxnKeggId)


	#rewrite keggModules to have same formating as biocyc pathways
	keggModules2 = dict()
	for moduleId in keggModules:
		keggModules2[moduleId] = {'RxnsInvolved':set(),'MetsInvolved':dict()}
		for rxnsInvolvedListIndex in range(len(keggModules[moduleId]['RxnsInvolved'])):
			for rxnsInvolvedSet in keggModules[moduleId]['RxnsInvolved'][rxnsInvolvedListIndex]:
				for rxnKeggId in rxnsInvolvedSet:
					keggModules2[moduleId]['RxnsInvolved'].add(rxnKeggId)
					keggModules2[moduleId]['MetsInvolved'][rxnKeggId] = keggModules[moduleId]['MetsInvolved'][rxnsInvolvedListIndex]
					
	soft_constraints_pathways = soft_constraints['metacyc']['pathways']|soft_constraints['kegg']['pathways']
	

	rxnsInPathsFromSoft = set()
	for eachPath in soft_constraints_pathways:
		if eachPath not in biocycPathways: continue
		for eachRxn in biocycPathways[eachPath]['RxnsInvolved']:
			if eachRxn not in metacycToRheaDict: continue
			for rheaRxn in metacycToRheaDict[eachRxn]: rxnsInPathsFromSoft.add(rheaRxn)
	rxnsInPathsFromMinPath = set() # reinitializing variable. Not sure if it is a good idea to use it during first scoring run
	rxnsInPathsNaive = set() # reinitializing variable. Not sure if it is a good idea to use it during first scoring run
			

	subunits = list()
	
	diamondResultDict = diamondResult.to_dict('records') # dict is faster than pandas in this application
	
	categoriesCount = {'subunits': 0, 'fromMinPath': 0, 'fromSoftConstraints': 0, 'fromNaivePaths':0, 'fromAnotation': 0, 'fromSwissProtIds':0}

	for repetition in [1,2]: #repeats the scoring two times. In the second time, favors reactions missing on almost complet pathways, and subunits of protein complexes. 


		#select one target_gene per "source_gene"
		currentRead = ''
		bestPerRead = list()
		notBestPerReadOdered = list()
		#for entryNumber in range(len(diamondResult)):
		#for entryNumber in diamondResultDict['source_gene']:
		for eachMatch in diamondResultDict:
		
			if not sprotSeqAndGenDict[eachMatch['target_gene']]['gene']: continue
		
			if eachMatch['source_gene'] != currentRead:
				#quando comeca o alinhamento de uma nova sequencia de proteina
				
				if currentRead and top60:
				
					#if currentRead == 'MNAFOKBH_02063': break
					
					#quando tem dados da sequencia anterior
					#break
					notBestPerReadOdered.extend(top60)
					
					bestMatch, notBestMatch = findBestPerRead(top60, swissProtIds, rxnsInPathsFromMinPath, rxnsInPathsFromSoft, rxnsInPathsNaive, subunits, sprotSeqAndGenDict, geneAndProteinNamePerSeqId, categoriesCount) # encontra o melhor match entre os os resultados mais frequentes.
					if bestMatch:
						bestPerRead.append(bestMatch)
						#notBestPerReadOdered.extend(notBestMatch)
					
					
				currentRead = eachMatch['source_gene']
				top60 = list()
			
			rheaSet = set()
			if eachMatch['target_gene'] not in uniprotToRheaRxns: continue
			for rheaId in uniprotToRheaRxns[eachMatch['target_gene']]:
				rheaSet.add(rheaId)
			
			#penalize proteins with weak evidance of existence. level 1: protein evidence; level 2: transcript evidence; level 3: homolugos; level 4: predict; level 5: uncertain
			if sprotSeqAndGenDict[eachMatch['target_gene']]['Existence evidence level'] == 3:
				newSocre = eachMatch['score']/1.1
			elif sprotSeqAndGenDict[eachMatch['target_gene']]['Existence evidence level'] == 4:
				newSocre = eachMatch['score']/2
			elif sprotSeqAndGenDict[eachMatch['target_gene']]['Existence evidence level'] == 5:
				newSocre = eachMatch['score']/4
			else: newSocre = eachMatch['score']
					
						
			if not top60:
				top60.append({'rxns':rheaSet,'count':1, 'uniprotEntry': eachMatch['target_gene'], 'source_gene': eachMatch['source_gene'], 'score':newSocre, 'gene': sprotSeqAndGenDict[eachMatch['target_gene']]['gene']})
			else:
				sucesso = 0
				for eachDict in top60:
					if eachDict['rxns'] == rheaSet:
						eachDict['count'] += 1
						if eachDict['score'] < newSocre:
							eachDict['score'] = newSocre
							eachDict['gene'] = sprotSeqAndGenDict[eachMatch['target_gene']]['gene']
							eachDict['uniprotEntry'] = eachMatch['target_gene']
						sucesso = 1
						break
				if sucesso == 0:
					top60.append({'rxns':rheaSet,'count':1, 'uniprotEntry': eachMatch['target_gene'], 'source_gene': eachMatch['source_gene'], 'score':newSocre, 'gene': sprotSeqAndGenDict[eachMatch['target_gene']]['gene']})
			
		if eachMatch:
			bestMatch, notBestMatch = findBestPerRead(top60, swissProtIds, rxnsInPathsFromMinPath, rxnsInPathsFromSoft, rxnsInPathsNaive, subunits, sprotSeqAndGenDict, geneAndProteinNamePerSeqId, categoriesCount)
			if bestMatch:
				bestPerRead.append(bestMatch)
				notBestPerReadOdered = notBestPerReadOdered + notBestMatch



		#sort bestPerRead, from high score to low score
		scores = list()
		for read in bestPerRead:
			scores.append(read['score'])
		scores.sort(reverse=True)
		
		bestPerReadSorted = list()
		for score in scores:
			for read in bestPerRead:
				if read['score'] == score:
					break
			bestPerReadSorted.append(read)
			bestPerRead.remove(read)
			
		#simplify bestPerReadSorted, by removing duplicates
		bestPerReadOrdered = list()
		for read1 in bestPerReadSorted:
			sucesso = 0
			for read2 in bestPerReadOrdered:
				if read2['rxns'] == read1['rxns']:
					sucesso = 1
					break
			if sucesso == 0:
				bestPerReadOrdered.append(read1)
			
			
		if repetition == 1 or repetition == 2:
			#find proteins that are subunits of protein complexes
			subunits = list()
			uniprotSubunitEntries = set()
			for eachDict in bestPerReadOrdered:
				if eachDict['uniprotEntry'] in uniprotSubunitEntries: continue
				subunitType = ''
				for subunitString in ['subunit', 'component']:
					if subunitString in sprotSeqAndGenDict[eachDict['uniprotEntry']]['title']:
						#to extract the subunit part
						if subunitString in sprotSeqAndGenDict[eachDict['uniprotEntry']]['title']:
						
							titleProvisional = sprotSeqAndGenDict[eachDict['uniprotEntry']]['title']
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
		
		
		
		#find pathways in organism. pathwaysPresent will be used outside this loop
		#encontra as reacoes 
		metacycRxnsInModel = set()
		for eachDict in bestPerReadOrdered:
			for rheaId in eachDict['rxns']:
				if rheaId not in rheaToMetacycDict: 
					#rheaId
					continue
				for rxnSyn in rheaToMetacycDict[rheaId]:
					metacycRxnsInModel.add(rxnSyn)
				
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
		toWrite = "#Metacyc pathway and reactions mapping file\n#Pathway	ReactionID"
		for eachPath in biocycPathways:
			if eachPath in superPathsDict: continue
			if len(biocycPathways[eachPath]['RxnsInvolved']) < 3: continue
			for eachRxn in biocycPathways[eachPath]['RxnsInvolved']:
				toWrite += "\n" + eachPath + '	' + eachRxn
		
		for eachPath in rxnsPerModules:
			for eachRxn in rxnsPerModules[eachPath]:
				toWrite += "\n" + eachPath + '	' + eachRxn
		
		mapping = open('metcycRxnPathwayMap.tsv', 'w')
		_ = mapping.write(toWrite)
		mapping.close()
		
		#find the minimum set of pathways thta explains the reactions
		print('encontra o conjunto minimo de pathways')
		MinPathResult = MinPathMain(anyfile = 'biocycRxnsInModel.tsv', mapfile = 'metcycRxnPathwayMap.tsv')
		
		#find the minimum set of pathways thta explains the reactions, including superpathways
		
		#creating file to be used as a mapping by the MinPath
		toWrite = "#Metacyc pathway and reactions mapping file\n#Pathway	ReactionID"
		for eachPath in biocycPathways:
			if eachPath in superPathsDict or eachPath in MinPathResult:
				for eachRxn in biocycPathways[eachPath]['RxnsInvolved']:
					toWrite += "\n" + eachPath + '	' + eachRxn
		
		mapping = open('metcycRxnPathwayMap.tsv', 'w')
		_ = mapping.write(toWrite)
		mapping.close()
		
		
		print('encontra o conjunto minimo de pathways')
		MinPathResult2 = MinPathMain(anyfile = 'biocycRxnsInModel.tsv', mapfile = 'metcycRxnPathwayMap.tsv')
		
		
		#find the pathwas in soft_constraints
		if len(soft_constraints_pathways) > 100: pathwaysPresent = (MinPathResult|MinPathResult2).intersection(soft_constraints_pathways)
		else: pathwaysPresent = (MinPathResult|MinPathResult2)
			
		if repetition == 1:
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
				if metacycRxn not in metacycToRheaDict: rxnsInPathsFromMinPath.add(metacycRxn)# continue
				else:
					for rheaRxn in metacycToRheaDict[metacycRxn]: rxnsInPathsFromMinPath.add(rheaRxn)
					
					
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


	##identify reactions that are very likely present in the model: because the gene is only responsible for a single reaction; because the pathway has been identified; because there is flux on the reactions when doing FBA.
	#find genes with only one reaction, and list these reactions
	singleRxnsInGenes = set()
	moreThanOneRxnsInGenes = set()
	for eachList in bestPerReadOrdered:
		if len(eachList['rxns']) == 1: 
			for rxnIdSyn in eachList['rxns']: singleRxnsInGenes.add(rxnIdSyn)
		else:
			for rxnIdSyn in eachList['rxns']: moreThanOneRxnsInGenes.add(rxnIdSyn)
	allWithGenes = moreThanOneRxnsInGenes|singleRxnsInGenes
		
	#do FBA, to see which reactions with genes have flux
	rxnWithGenesAndFlux = set()
	with cobraModel as model:
		
		boundaryIds = {rxnCobra.id for rxnCobra in model.boundary}
	
		for rxnCobra in model.reactions:
			if len(rxnCobra.compartment) == 2: continue
			if rxnCobra.id in boundaryIds: continue
			if 'rhea' in rxnCobra.annotation:
				if allWithGenes.intersection(rxnCobra.annotation['rhea']): continue
			if 'kegg' in rxnCobra.annotation:
				if allWithGenes.intersection(rxnCobra.annotation['kegg']): continue
			if 'metacyc' in rxnCobra.annotation:
				if allWithGenes.intersection(rxnCobra.annotation['metacyc']): continue
			if rxnCobra.lower_bound < 0: rxnCobra.lower_bound = -0.1
			if rxnCobra.upper_bound > 0: rxnCobra.upper_bound = 0.1

		fba = model.optimize()
		
		for rxnCobra in cobraModel.reactions:
			if abs(fba.fluxes[rxnCobra.id]) > 0 and (rxnCobra.lower_bound < -0.1 or  rxnCobra.upper_bound > 0.1):
				if 'rhea' in rxnCobra.annotation:
					for i in rxnCobra.annotation['rhea']: rxnWithGenesAndFlux.add(i)
				if 'kegg' in rxnCobra.annotation:
					for i in rxnCobra.annotation['kegg']: rxnWithGenesAndFlux.add(i)
				if 'metacyc' in rxnCobra.annotation:
					for i in rxnCobra.annotation['metacyc']: rxnWithGenesAndFlux.add(i)
					
	relevantRxnsPerGene = singleRxnsInGenes|rxnsPresentInPathways|rxnWithGenesAndFlux
	
	#find rxns associated with genes, where one of the reaction has strong evidence. Consider the other reactions as probableArtefact. At least one of the reactions per gene has a positive score.
	artefacts = set()
	for eachList in bestPerReadOrdered:
		if len(eachList['rxns']) <= 1: continue
		probableArtefact = eachList['rxns'] - relevantRxnsPerGene
		if probableArtefact:
			if probableArtefact and eachList['rxns'] == probableArtefact: eachList['score'] = eachList['score']/2
			else: 
				for i in probableArtefact: artefacts.add(i)
				
	
	bestPerReadSimplified = list()
	rxnsInBest = set()
	for eachList in bestPerReadOrdered:
		
		rxnsSet = set()
		for rxnIdSyn in eachList['rxns']:
			if rxnIdSyn in rheaIdToreframedId:
				rxnsSet = rxnsSet|rheaIdToreframedId[rxnIdSyn]
				
		if rxnsSet: 
			eachList['rxns'] = rxnsSet
			bestPerReadSimplified.append(eachList)
			for i in rxnsSet: rxnsInBest.add(i)
			
	
	notBestPerRead = list()
	rxnsInNotBest = set()
	for eachList in notBestPerReadOdered:
		
		rxnsSet = set()
		for rxnId in eachList['rxns']:
			if rxnId in rheaIdToreframedId:
				rxnsSet = rxnsSet|rheaIdToreframedId[rxnId]
				
		if rxnsSet: 
			eachList['rxns'] = rxnsSet
			notBestPerRead.append(eachList)
			for i in rxnsSet: rxnsInNotBest.add(i)
			
	#to retrive association between rxns and genes
	rheaIdToGene = {'bestPerReadSimplified':bestPerReadSimplified, 'notBestPerRead':notBestPerRead}
	
			
	#comeca a trabalhar apenas com reacoes presentes no modelo universal
	medianValue = median([eachDict['score'] for eachDict in bestPerReadSimplified])

	make100Zero = 100/medianValue #100 might lead to wrong results. increased to 150. Maybe 200 to Trembl and 100 to SwissProt?
	make200Zero = 200/medianValue #100 might lead to wrong results. increased to 150. Maybe 200 to Trembl and 100 to SwissProt?

	#normalize
	for eachDict in bestPerReadSimplified:
		if eachDict['uniprotEntry'] in swissProtIds: eachDict['scoreNormalized'] = round((eachDict['score']/medianValue) - make100Zero, 4)
		else: eachDict['scoreNormalized'] = round((eachDict['score']/medianValue) - make200Zero, 4)
		
		if eachDict['scoreNormalized'] > 0: eachDict['scoreNormalized'] = eachDict['scoreNormalized']*1.5
		elif eachDict['scoreNormalized'] < 0: eachDict['scoreNormalized'] = eachDict['scoreNormalized']/10
		elif eachDict['scoreNormalized'] == 0: eachDict['scoreNormalized'] = 0.01
		
	uniprotScoreNormalized = dict()
	for eachDict in bestPerReadSimplified:
		uniprotScoreNormalized[eachDict['gene']] = eachDict['scoreNormalized']
	
	
	#finds subunits 
	maisDeDuasSubunidades = set()
	for eachDictIndex1 in range(len(subunits)-1):
		for eachDictIndex2 in range(eachDictIndex1+1, len(subunits)):
			if not subunits[eachDictIndex1]['gene'] or not subunits[eachDictIndex2]['gene']: continue
			if subunits[eachDictIndex1]['rxns'].intersection(subunits[eachDictIndex2]['rxns']) and subunits[eachDictIndex1]['subunit'] != subunits[eachDictIndex2]['subunit'] and subunits[eachDictIndex1]['gene'] != subunits[eachDictIndex2]['gene']:
				maisDeDuasSubunidades.add(subunits[eachDictIndex1]['uniprotEntry'])
				maisDeDuasSubunidades.add(subunits[eachDictIndex2]['uniprotEntry'])
	

	#finding rhea without a complex
	rxnsScores = dict() #dict to store score of rhea reactions
	for eachDict in bestPerReadSimplified:
		if eachDict['uniprotEntry'] in uniprotSubunitEntries:
			for rheaId in eachDict['rxns']:
				if eachDict['uniprotEntry'] in maisDeDuasSubunidades:
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
	for rheaId in rxnsScores:
		if rxnsScores[rheaId] > 10: rxnsScores[rheaId] = 10
		
	#reduce score of possible artefacts
	for rheaId in artefacts:
		if rheaId not in rheaIdToreframedId: continue
		for idInModel in rheaIdToreframedId[rheaId]:
			if rxnsScores[idInModel] > 0: rxnsScores[idInModel] = -0.01
				
			
	minScoreBase = min(rxnsScores.values())
	
	#reactions not associated with the best match of each protein sequence
	for rheaId in rxnsInNotBest - rxnsInBest:
		rxnsScores[rheaId] = minScoreBase -0.1
	
	if minScoreBase > 0: minScore = -0.1
	else: minScore = minScoreBase - 0.1
	

	rxnsCobraToreframedModel = dict()
	softPositiveNewRxn = set()
	
	
	if soft_constraints:
		for rxnId in soft_constraints['kegg']['1 intersection']:
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, 1, softPositiveNewRxn)
	
		for rxnId in soft_constraints['kegg']['1'] - soft_constraints['kegg']['1 intersection']:
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, -0.001, softPositiveNewRxn)
		for rxnId in soft_constraints['kegg']['0.9 to 0.99'] - (soft_constraints['kegg']['1']|soft_constraints['kegg']['1 intersection']):
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, minScore - 0.1, softPositiveNewRxn)
		for rxnId in soft_constraints['kegg']['0.8 to 0.89'] - (soft_constraints['kegg']['1']|soft_constraints['kegg']['1 intersection']|soft_constraints['kegg']['0.9 to 0.99']):
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, minScore - 0.3, softPositiveNewRxn)

		for rxnId in soft_constraints['metacyc']['1'] - soft_constraints['kegg']['1']:
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, -0.01, softPositiveNewRxn)
		for rxnId in soft_constraints['metacyc']['0.9 to 0.99'] - (soft_constraints['metacyc']['1']):
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, minScore - 0.3, softPositiveNewRxn)
		for rxnId in soft_constraints['metacyc']['0.8 to 0.89'] - (soft_constraints['metacyc']['1']|soft_constraints['metacyc']['0.9 to 0.99']):
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, minScore - 0.6, softPositiveNewRxn)
		

	#find the reactions missing in pathways that are almost complete. These reactions will receive the minimum score
	metacycPositiveScore = set()
	for rheaId in rxnsInBest|rxnsInNotBest:
		if rheaId not in rheaToMetacycDict: continue
		for rxnSyn in rheaToMetacycDict[rheaId]: metacycPositiveScore.add(rxnSyn)
	
	
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
		if rxnMetacycId in metacycToRheaDict:
			for rxnRheaId in metacycToRheaDict[rxnMetacycId]:
				rxnsRheaFaltandoPathways.add(rxnRheaId)
				
	rxnsMissingInPathways = rxnsRheaFaltandoPathways
	
	minScore = min(rxnsScores.values())
	for rheaId in rxnsMissingInPathways:
		if rheaId in rxnsScores: 
			if rxnsScores[rheaId] < 0: rxnsScores[rheaId]=rxnsScores[rheaId]/2
		else: rxnsScores[rheaId] = minScore
		
	
	relevantRxnsPerGeneInModel = set()
	relevantRxnsPerGene = singleRxnsInGenes|rxnsPresentInPathways
	for rxnCobra in cobraModel.reactions:
		if 'rhea' in rxnCobra.annotation:
			if relevantRxnsPerGene.intersection(rxnCobra.annotation['rhea']): relevantRxnsPerGeneInModel.add('R_'+rxnCobra.id.replace('-','__45__').replace('.','__46__').replace('+','__43__'))
		if 'kegg' in rxnCobra.annotation:
			if relevantRxnsPerGene.intersection(rxnCobra.annotation['kegg']): relevantRxnsPerGeneInModel.add('R_'+rxnCobra.id.replace('-','__45__').replace('.','__46__').replace('+','__43__'))
		if 'metacyc' in rxnCobra.annotation:
			if relevantRxnsPerGene.intersection(rxnCobra.annotation['metacyc']): relevantRxnsPerGeneInModel.add('R_'+rxnCobra.id.replace('-','__45__').replace('.','__46__').replace('+','__43__'))
	

	return rxnsScores, rheaIdToGene, rxnsCobraToreframedModel, soft_constraints_pathways, bestPerReadSimplified, softPositiveNewRxn, relevantRxnsPerGeneInModel
	
