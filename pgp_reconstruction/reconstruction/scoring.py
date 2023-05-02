import pickle
from pgp_reconstruction import project_dir
from pgp_reconstruction.dependencies.MinPath_master.MinPath import MinPathMain
from statistics import median
import sys
import os
import datetime


def findBestPerRead(top60, swissProtIds, rxnsInSoftSyn, subunits, sprotSeqAndGenDict, geneAndProteinNamePerSeqId):
	#find best match and alternative matches of each read

	freqRxnDict = dict()
	for eachDict in top60:
		for rxn in eachDict['rxns']:
			if rxn not in freqRxnDict: freqRxnDict[rxn] = 0
			freqRxnDict[rxn] += eachDict['count']

	top60v2 = list()
	for eachDict in top60:
		for rxn in freqRxnDict:
			if rxn in eachDict['rxns'] and freqRxnDict[rxn] < 10 and rxn not in rxnsInSoftSyn:
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
		
		
		#da prioridade ao match com maior intersecao com reacoes anotadas pelo biocyc/kegg, quando ha dois matchs, e os dois tem scores proximos
		swissProtId = list()
		prokkaMatching = list()
		intersectionCount = list()
		freqCount = list()
		for eachDict in top50Simplified:
			if eachDict['score'] >=  maxScore*0.5:
			
				#if input is prokka annotation, and gene matched is the one annotated by prokka. 
				if eachDict['query_gene'] in geneAndProteinNamePerSeqId and geneAndProteinNamePerSeqId[eachDict['query_gene']]['gene'] and geneAndProteinNamePerSeqId[eachDict['query_gene']]['protein name 1'] and (geneAndProteinNamePerSeqId[eachDict['query_gene']]['gene'].lower() == sprotSeqAndGenDict[eachDict['uniprotEntry']]['gene'].lower() or sprotSeqAndGenDict[eachDict['uniprotEntry']]['title'].lower() in geneAndProteinNamePerSeqId[eachDict['query_gene']]['protein name 1'].lower()):
					prokkaMatching.append(1)
				elif eachDict['query_gene'] in geneAndProteinNamePerSeqId and geneAndProteinNamePerSeqId[eachDict['query_gene']]['gene'] and geneAndProteinNamePerSeqId[eachDict['query_gene']]['protein name 2'] and (sprotSeqAndGenDict[eachDict['uniprotEntry']]['title'].lower() in geneAndProteinNamePerSeqId[eachDict['query_gene']]['protein name 2'].lower()):
					prokkaMatching.append(1)
				else: prokkaMatching.append(0)
			
				freqCount.append(len(eachDict['rxns'].intersection(rxnsInSoftSyn)))
				intersectionCount.append(len(eachDict['rxns'].intersection(rxnsInSoftSyn))/len(eachDict['rxns']))
				if eachDict['uniprotEntry'] in swissProtIds: swissProtId.append(1)
				else: swissProtId.append(0)
			else:
				freqCount.append(0)
				intersectionCount.append(0)
				swissProtId.append(0)


		fromProkka = 0
		fromSwiss = 0
		if sum(prokkaMatching) > 0:
			fromProkka = 1
			for valueIndex in range(len(prokkaMatching)):
				if not prokkaMatching[valueIndex]:
					swissProtId[valueIndex] = 0
					freqCount[valueIndex] = 0
					intersectionCount[valueIndex] = 0
		
		if sum(swissProtId) > 0:
			fromSwiss = 1
			for valueIndex in range(len(swissProtId)):
				if not swissProtId[valueIndex]:
					freqCount[valueIndex] = 0
					intersectionCount[valueIndex] = 0
						
		for valueIndex in range(len(freqCount)):
			if freqCount[valueIndex] != max(freqCount):
				freqCount[valueIndex] = 0
				intersectionCount[valueIndex] = 0		
				
		for valueIndex in range(len(intersectionCount)):
			if intersectionCount[valueIndex] != max(intersectionCount):
				freqCount[valueIndex] = 0
				intersectionCount[valueIndex] = 0
						
		#take the bestMatch the highest score between what is left
		maximumScore = 0
		maximumScoreIndex = 0
		if sum(freqCount):
			for valueIndex in range(len(freqCount)):
				if not freqCount[valueIndex]: continue
				if top50Simplified[valueIndex]['score'] > maximumScore:
					maximumScore = top50Simplified[valueIndex]['score']
					maximumScoreIndex = valueIndex
		else:
			prokkaInterSwissProt = list()
			for valueIndex in range(len(swissProtId)):
				if swissProtId[valueIndex] == 1 and prokkaMatching[valueIndex] == 1: prokkaInterSwissProt.append(valueIndex)
		
			if not prokkaInterSwissProt:
				for valueIndex in range(len(swissProtId)):
					if sum(swissProtId) and not swissProtId[valueIndex]: continue
					if top50Simplified[valueIndex]['score'] > maximumScore:
						maximumScore = top50Simplified[valueIndex]['score']
						maximumScoreIndex = valueIndex
			else:
				for valueIndex in prokkaInterSwissProt:
					if top50Simplified[valueIndex]['score'] > maximumScore:
						maximumScore = top50Simplified[valueIndex]['score']
						maximumScoreIndex = valueIndex
			
			
		#separate in two lists without intersection
		bestMatch = top50Simplified[maximumScoreIndex]
		del top50Simplified[maximumScoreIndex]
		
		notBestMatch = list()
		for eachDict in top50Simplified:
			if not eachDict['rxns'].issubset(bestMatch['rxns']):
				notBestMatch.append(eachDict)
				
				
		#verifica se algum dos que estao em notBestMatch eh uma subunit que faltou no processo anterior. Se for, troca bestMatch
		
		if subunits and ('subunit' not in sprotSeqAndGenDict[eachDict['uniprotEntry']]['title'] and 'component' not in sprotSeqAndGenDict[eachDict['uniprotEntry']]['title']):
		
			sucesso = 0
			for eachDict in notBestMatch:
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
							sucesso = 1
							break
				if sucesso == 1:
					break
			if sucesso == 1:
				notBestMatch.append(bestMatch)
				bestMatch = eachDict
				notBestMatch.remove(eachDict)
			
		return bestMatch, notBestMatch, fromProkka, fromSwiss
	else: return 0, 0, 0, 0
		
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

def reaction_scoring(diamondResult, geneAndProteinNamePerSeqId, cobraModel, reframedModel, uniprotToRheaRxns, soft_constraints, rxnsInSoftSyn, verbose):
	
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
	rxnsInSoftSyn = set() # reinitializing variable. Not sure if it is a good idea to use it during first scoring run
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
					

	subunits = list()
	soft_constraints_pathways = soft_constraints['metacyc']['pathways']|soft_constraints['kegg']['pathways']
	TotalProkka, TotalSwiss = 0,0
	
	diamondResultDict = diamondResult.to_dict('records') # dict is faster than pandas in this application
	
	for repetition in [1,2]: #repeats the scoring two times. In the second time, favors reactions missing on almost complet pathways, and subunits of protein complexes. 


		#select one BiGG_gene per "query_gene"
		currentRead = ''
		bestPerRead = list()
		notBestPerReadOdered = list()
		#for entryNumber in range(len(diamondResult)):
		#for entryNumber in diamondResultDict['query_gene']:
		for eachMatch in diamondResultDict:
		
			if not sprotSeqAndGenDict[eachMatch['BiGG_gene']]['gene']: continue
		
			if eachMatch['query_gene'] != currentRead:
				#quando comeca o alinhamento de uma nova sequencia de proteina
				
				if currentRead and top60:
				
					#if currentRead == 'MNAFOKBH_02063': break
					
					#quando tem dados da sequencia anterior
					#break
					notBestPerReadOdered.extend(top60)
					
					bestMatch, notBestMatch, fromProkka, fromSwiss = findBestPerRead(top60, swissProtIds, rxnsInSoftSyn, subunits, sprotSeqAndGenDict, geneAndProteinNamePerSeqId) # encontra o melhor match entre os os resultados mais frequentes.
					if bestMatch:
						bestPerRead.append(bestMatch)
						#notBestPerReadOdered.extend(notBestMatch)
						TotalProkka += fromProkka
						TotalSwiss += fromSwiss
					
					
				currentRead = eachMatch['query_gene']
				top60 = list()
			
			rheaSet = set()
			if eachMatch['BiGG_gene'] not in uniprotToRheaRxns: continue
			for rheaId in uniprotToRheaRxns[eachMatch['BiGG_gene']]:
				rheaSet.add(rheaId)
			
			#penalize proteins with weak evidance of existence. level 1: protein evidence; level 2: transcript evidence; level 3: homolugos; level 4: predict; level 5: uncertain
			if sprotSeqAndGenDict[eachMatch['BiGG_gene']]['Existence evidence level'] == 3:
				newSocre = eachMatch['score']/1.1
			elif sprotSeqAndGenDict[eachMatch['BiGG_gene']]['Existence evidence level'] == 4:
				newSocre = eachMatch['score']/2
			elif sprotSeqAndGenDict[eachMatch['BiGG_gene']]['Existence evidence level'] == 5:
				newSocre = eachMatch['score']/4
			else: newSocre = eachMatch['score']
					
						
			if not top60:
				top60.append({'rxns':rheaSet,'count':1, 'uniprotEntry': eachMatch['BiGG_gene'], 'query_gene': eachMatch['query_gene'], 'score':newSocre, 'gene': sprotSeqAndGenDict[eachMatch['BiGG_gene']]['gene']})
			else:
				sucesso = 0
				for eachDict in top60:
					if eachDict['rxns'] == rheaSet:
						eachDict['count'] += 1
						if eachDict['score'] < newSocre:
							eachDict['score'] = newSocre
							eachDict['gene'] = sprotSeqAndGenDict[eachMatch['BiGG_gene']]['gene']
							eachDict['uniprotEntry'] = eachMatch['BiGG_gene']
						sucesso = 1
						break
				if sucesso == 0:
					top60.append({'rxns':rheaSet,'count':1, 'uniprotEntry': eachMatch['BiGG_gene'], 'query_gene': eachMatch['query_gene'], 'score':newSocre, 'gene': sprotSeqAndGenDict[eachMatch['BiGG_gene']]['gene']})
			
		if eachMatch:
			bestMatch, notBestMatch, fromProkka, fromSwiss = findBestPerRead(top60, swissProtIds, rxnsInSoftSyn, subunits, sprotSeqAndGenDict, geneAndProteinNamePerSeqId)
			if bestMatch:
				bestPerRead.append(bestMatch)
				notBestPerReadOdered = notBestPerReadOdered + notBestMatch
				TotalProkka += fromProkka
				TotalSwiss += fromSwiss



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
		
		
		if repetition == 1:

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
				if metacycRxn not in metacycToRheaDict: rxnsInSoftSyn.add(metacycRxn)# continue
				else:
					for rheaRxn in metacycToRheaDict[metacycRxn]: rxnsInSoftSyn.add(rheaRxn)


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
			if len(rxnCobra.compartments) == 2: continue
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
				
			
	print('\n' + 'TotalProkka = ' + str(TotalProkka))
	print('TotalSwiss = ' + str(TotalSwiss)+ '\n')
	
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
	
	
	if 1==1:
		#if soft_constraints:
		'''
		#boundaries
		alterado = set()
		for rxn in cobraModel.boundary:
			if rxn.compartments != {'e'}: continue
			idInModel = 'R_'+rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
			for met in rxn.metabolites:
				
				if met.id[-2] in soft_constraints['User constraint']['To Include'] : rxnsScores[idInModel] = 10
				elif met.id in soft_constraints['User constraint']['To Include'] : rxnsScores[idInModel] = 10
				elif met.id[-2] in soft_constraints['User constraint']['To avoid'] : rxnsScores[idInModel] = -10
				elif met.id in soft_constraints['User constraint']['To avoid'] : rxnsScores[idInModel] = -10
				elif met.id[-2] in soft_constraints['metacyc met']['1'] and (idInModel not in rxnsScores or rxnsScores[idInModel] < -0.1): rxnsScores[idInModel] = 0
				elif met.id[-2] in soft_constraints['metacyc met']['0.9 to 0.99'] and (idInModel not in rxnsScores or rxnsScores[idInModel] < -0.2): rxnsScores[idInModel] = minScore - 0.2
				elif met.id[-2] in soft_constraints['metacyc met']['0.8 to 0.89'] and (idInModel not in rxnsScores or rxnsScores[idInModel] < -0.3): rxnsScores[idInModel] = minScore -0.3
				elif met.id[-2] in soft_constraints['metacyc met']['0.7 to 0.79'] and (idInModel not in rxnsScores or rxnsScores[idInModel] < -0.5): rxnsScores[idInModel] = minScore -0.5
				elif met.id[-2] in soft_constraints['metacyc met']['0.6 to 0.69'] and (idInModel not in rxnsScores or rxnsScores[idInModel] < -0.7): rxnsScores[idInModel] = minScore -0.7
				elif met.id[-2] in soft_constraints['metacyc met']['0.5 to 0.59'] and (idInModel not in rxnsScores or rxnsScores[idInModel] < -0.8): rxnsScores[idInModel] = minScore -0.8
				elif met.id[-2] in soft_constraints['metacyc met']['lower than 0.5'] and (idInModel not in rxnsScores or rxnsScores[idInModel] < -1): rxnsScores[idInModel] = minScore -1
		'''
		for rxnId in soft_constraints['kegg']['1 intersection']:
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, 1, softPositiveNewRxn)
	
		for rxnId in soft_constraints['kegg']['1'] - soft_constraints['kegg']['1 intersection']:
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, -0.001, softPositiveNewRxn)
		for rxnId in soft_constraints['kegg']['0.9 to 0.99'] - (soft_constraints['kegg']['1']|soft_constraints['kegg']['1 intersection']):
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, minScore - 0.1, softPositiveNewRxn)
		for rxnId in soft_constraints['kegg']['0.8 to 0.89'] - (soft_constraints['kegg']['1']|soft_constraints['kegg']['1 intersection']|soft_constraints['kegg']['0.9 to 0.99']):
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, minScore - 0.3, softPositiveNewRxn)
		#for rxnId in soft_constraints['kegg']['0.7 to 0.79'] - (soft_constraints['kegg']['1']|soft_constraints['kegg']['1 intersection']|soft_constraints['kegg']['0.9 to 0.99']|soft_constraints['kegg']['0.8 to 0.89']):
		#	changeScoreOfSoft(rxnId, rxnsScores, reframedModel, minScore - 0.5, softPositiveNewRxn)
		#for rxnId in soft_constraints['kegg']['0.6 to 0.69'] - (soft_constraints['kegg']['1']|soft_constraints['kegg']['1 intersection']|soft_constraints['kegg']['0.9 to 0.99']|soft_constraints['kegg']['0.8 to 0.89']|soft_constraints['kegg']['0.7 to 0.79']):
		#	changeScoreOfSoft(rxnId, rxnsScores, reframedModel, minScore - 0.7, softPositiveNewRxn)
		#for rxnId in soft_constraints['kegg']['0.5 to 0.59'] - (soft_constraints['kegg']['1']|soft_constraints['kegg']['1 intersection']|soft_constraints['kegg']['0.9 to 0.99']|soft_constraints['kegg']['0.8 to 0.89']|soft_constraints['kegg']['0.7 to 0.79']|soft_constraints['kegg']['0.6 to 0.69']):
		#	changeScoreOfSoft(rxnId, rxnsScores, reframedModel, minScore - 1, softPositiveNewRxn)

		for rxnId in soft_constraints['metacyc']['1'] - soft_constraints['kegg']['1']:
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, -0.01, softPositiveNewRxn)
		for rxnId in soft_constraints['metacyc']['0.9 to 0.99'] - (soft_constraints['metacyc']['1']):
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, minScore - 0.3, softPositiveNewRxn)
		for rxnId in soft_constraints['metacyc']['0.8 to 0.89'] - (soft_constraints['metacyc']['1']|soft_constraints['metacyc']['0.9 to 0.99']):
			changeScoreOfSoft(rxnId, rxnsScores, reframedModel, minScore - 0.6, softPositiveNewRxn)
		#for rxnId in soft_constraints['metacyc']['0.7 to 0.79'] - (soft_constraints['metacyc']['1']|soft_constraints['metacyc']['0.9 to 0.99']|soft_constraints['metacyc']['0.8 to 0.89']):
		#	changeScoreOfSoft(rxnId, rxnsScores, reframedModel, minScore - 1, softPositiveNewRxn)
		#for rxnId in soft_constraints['metacyc']['0.6 to 0.69'] - (soft_constraints['metacyc']['1']|soft_constraints['metacyc']['0.9 to 0.99']|soft_constraints['metacyc']['0.8 to 0.89']|soft_constraints['metacyc']['0.7 to 0.79']):
		#	changeScoreOfSoft(rxnId, rxnsScores, reframedModel, minScore - 1.5, softPositiveNewRxn)
		#for rxnId in soft_constraints['metacyc']['0.5 to 0.59'] - (soft_constraints['metacyc']['1']|soft_constraints['metacyc']['0.9 to 0.99']|soft_constraints['metacyc']['0.8 to 0.89']|soft_constraints['metacyc']['0.7 to 0.79']|soft_constraints['metacyc']['0.6 to 0.69']):
		#	changeScoreOfSoft(rxnId, rxnsScores, reframedModel, minScore - 2, softPositiveNewRxn)
		

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
	

	'''
	#faz com que metabolitos do controle negativo do biolog sejam mais facilmente incluidos no meio externo
	minimumMediaIds = ['ca2', 'cl', 'cobalt2', 'cu2', 'fe2', 'fe3', 'h', 'h2o', 'k', 'mg2', 'mobd', 'na1', 'nh4', 'ni2', 'pi', 'so4', 'zn2', 'o2']
	#mahirMedia = ['photon', '17632', '29101', '43474', '29103', '16189', '18420', '17996', '16947']
	mahirMedia = []
	minimumMediaFacultativeIds = ['glc__D', 'glyc', 'lac__D', 'lac__L']
	
	minimumMediaRxns = dict()
	for met in cobraModel.metabolites:
		if met.compartment == 'e':
			for metId in minimumMediaIds+minimumMediaFacultativeIds+mahirMedia:
				if metId not in minimumMediaRxns:
					minimumMediaRxns[metId] = {'exchange':{'products':set(), 'reactants':set()},'transport':{'products':set(), 'reactants':set()}}
				
				for db in ['bigg', 'chebi']:
					if db in met.annotation and metId in met.annotation[db]:
						if metId in minimumMediaIds+minimumMediaFacultativeIds+mahirMedia:
							for rxn in met.reactions:
								if len(rxn.metabolites) == 1:
									if met in rxn.products: minimumMediaRxns[metId]['exchange']['products'].add(rxn)
									elif met in rxn.reactants: minimumMediaRxns[metId]['exchange']['reactants'].add(rxn)
								if len(rxn.compartments) > 1:
									if met in rxn.products: minimumMediaRxns[metId]['transport']['products'].add(rxn)
									elif met in rxn.reactants: minimumMediaRxns[metId]['transport']['reactants'].add(rxn)
	
	for metId in minimumMediaRxns:
		for rxnType in ['exchange', 'transport']:
			if len(minimumMediaRxns[metId][rxnType]) == 1:
				for rxn in minimumMediaRxns[metId][rxnType]['reactants']|minimumMediaRxns[metId][rxnType]['products']:
					idInModel = 'R_' + rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
					if rxn in minimumMediaRxns[metId][rxnType]['products']: rxnsScores[idInModel] = -10
					else:
						if metId in minimumMediaIds+mahirMedia:
							if idInModel not in rxnsScores or rxnsScores[idInModel] < 1: rxnsScores[idInModel] = 5
						if metId in minimumMediaFacultativeIds:
							if idInModel not in rxnsScores or rxnsScores[idInModel] < 0.01: rxnsScores[idInModel] = 0.01
			if len(minimumMediaRxns[metId][rxnType]) > 1:
				for rxn in minimumMediaRxns[metId][rxnType]['reactants']|minimumMediaRxns[metId][rxnType]['products']:
					idInModel = 'R_' + rxn.id.replace('-','__45__').replace('.','__46__').replace('+','__43__')
					if rxn in minimumMediaRxns[metId][rxnType]['products']: rxnsScores[idInModel] = -10
					else:
						if metId in minimumMediaIds+mahirMedia:
							if idInModel not in rxnsScores or rxnsScores[idInModel] < 1: rxnsScores[idInModel] = 5
						if metId in minimumMediaFacultativeIds:
							if idInModel not in rxnsScores or rxnsScores[idInModel] < 0.01: rxnsScores[idInModel] = 0.01
	'''

	return rxnsScores, rheaIdToGene, rxnsCobraToreframedModel, soft_constraints_pathways, bestPerReadSimplified, softPositiveNewRxn, relevantRxnsPerGeneInModel
	
