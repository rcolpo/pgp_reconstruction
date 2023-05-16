#procura o core metabolismo de uma especie, e cria um arquivo com as reacoes do kegg e biocyc encontradas em organismos filogeneticamente proximos.
import os
import pickle
from pgp_reconstruction import project_dir
import re
import sys



def findAllMetsAndRxnsTranslations(cobraModel):
	#find RXNs translated to kegg and biocyc
	dbToRheaRxns = {'kegg':dict(), 'metacyc':dict()}
	for cobraRxn in cobraModel.reactions: 
		for eachDb in ['metacyc', 'kegg']:
			if eachDb in cobraRxn.annotation:
				for biocycRxnSyn in cobraRxn.annotation[eachDb]:
					idInDb = biocycRxnSyn.replace('META:','')
					if idInDb not in dbToRheaRxns[eachDb]:
						dbToRheaRxns[eachDb][idInDb] = set()
					dbToRheaRxns[eachDb][idInDb].add(cobraRxn.id.replace('_forwardTemp','').replace('_reverseTemp','').replace('_bidirectionalCopy',''))
	
	#find METs translated to kegg and biocyc
	dbToRheaMets = {'kegg':dict(), 'metacyc':dict()}
	for met in cobraModel.metabolites:
		for eachDb in ['kegg', 'metacyc']:
			if eachDb in met.annotation:
				for keggId in met.annotation[eachDb]:
					keggId = keggId.replace('META:','')
					if keggId not in dbToRheaMets:
						dbToRheaMets[eachDb][keggId] = list()
					dbToRheaMets[eachDb][keggId].append(met.id[:-2]) #tirar informacao de compartimento
					

	rheaWithSameSyn = list()
	for eachDb in dbToRheaRxns:
		for keggId in dbToRheaRxns[eachDb]:
			if len(dbToRheaRxns[eachDb][keggId]) > 1:
				rheaWithSameSyn.append(dbToRheaRxns[eachDb][keggId])
				
	return dbToRheaRxns, dbToRheaMets, rheaWithSameSyn
	
	

def findKeggRxnsInReferenceOrgs(taxonomicalLevelPerDb, dbToRheaRxns, allSpecKegg):
	#findKeggRxnsInReferenceOrgs(organismNameKeggList, dbToRheaRxns, allSpecKegg)
	
	keggAllSpec = dict()
	keggAllRxns = set()
	pathwaysKegg = set()
	for taxLavel in taxonomicalLevelPerDb['Kegg']:
		#encontra as reacoes anotadas pelo KEGG presentes no organismo
		
		keggAllSpec[taxLavel] = list()
		
		for organismNameKegg in taxonomicalLevelPerDb['Kegg'][taxLavel]:
		
			rheaInKeggSpec = set()
			pathwaysKeggSet = set()
			for pathway in allSpecKegg[organismNameKegg]:
				pathwaysKeggSet.add(pathway)
				for rxnKegg in allSpecKegg[organismNameKegg][pathway]:
					if rxnKegg not in dbToRheaRxns['kegg']:
						continue
					else:
						for rxnRheaId in dbToRheaRxns['kegg'][rxnKegg]:
							rheaInKeggSpec.add(rxnRheaId)
							keggAllRxns.add(rxnRheaId)
			keggAllSpec[taxLavel].append(rheaInKeggSpec)
			pathwaysKegg = pathwaysKegg|pathwaysKeggSet
		if not keggAllSpec[taxLavel]: del keggAllSpec[taxLavel]
		
	keggRxnsFreq = dict()
	for rxnId in keggAllRxns:
		freq = list()
		for taxLevel in keggAllSpec:
			count = 0
			for specList in keggAllSpec[taxLevel]:
				if rxnId in specList: count += 1
			if count == 0: freq.append(0)
			else: freq.append(count/len(keggAllSpec[taxLevel]))
		keggRxnsFreq[rxnId] = sum(freq)/len(freq)
			
		
	keggRxnsFreqSlots = {'lower than 0.5': set(), '0.5 to 0.59': set(), '0.6 to 0.69': set(), '0.7 to 0.79': set(), '0.8 to 0.89': set(), '0.9 to 0.99': set(), '1': set()}
	
	for eachRxn in keggRxnsFreq:
		if keggRxnsFreq[eachRxn] == 1: keggRxnsFreqSlots['1'].add(eachRxn)
		elif keggRxnsFreq[eachRxn] >= 0.9:  keggRxnsFreqSlots['0.9 to 0.99'].add(eachRxn)
		elif keggRxnsFreq[eachRxn] >= 0.8:  keggRxnsFreqSlots['0.8 to 0.89'].add(eachRxn)
		elif keggRxnsFreq[eachRxn] >= 0.7:  keggRxnsFreqSlots['0.7 to 0.79'].add(eachRxn)
		elif keggRxnsFreq[eachRxn] >= 0.6:  keggRxnsFreqSlots['0.6 to 0.69'].add(eachRxn)
		elif keggRxnsFreq[eachRxn] >= 0.5:  keggRxnsFreqSlots['0.5 to 0.59'].add(eachRxn)
		else: keggRxnsFreqSlots['lower than 0.5'].add(eachRxn)
			
	return keggRxnsFreqSlots, pathwaysKegg


def findBiocycRxnsInReferenceOrgs(taxonomicalLevelPerDb, dbToRheaRxns, dbToRheaMets, allSpecBioCyc, biocycPathways):
	
	biocycRxnsPerTax = dict()
	biocycMetsPerTax = dict()
	biocycAllRxns = set()
	biocycAllMets = set()
	pathwaysBiocyc = set()
	
	for taxLavel in taxonomicalLevelPerDb['Biocyc']:
		#encontra as reacoes anotadas pelo biocyc presentes no organismo
		
		biocycRxnsPerTax[taxLavel] = list()
		biocycMetsPerTax[taxLavel] = list()
		
		for organismNameBiocyc in taxonomicalLevelPerDb['Biocyc'][taxLavel]:
		
			pathwaysBiocycSet = set(allSpecBioCyc[organismNameBiocyc].keys()).intersection(biocycPathways)
			for path in pathwaysBiocycSet: pathwaysBiocyc.add(path)
			
			rheaRxnsInBiocycSpec = set()
			for rxnbiocyc in allSpecBioCyc[organismNameBiocyc]['rxnsInPathways']|allSpecBioCyc[organismNameBiocyc]['notInPathwaysNotTransportRxns']:
				if rxnbiocyc not in dbToRheaRxns['metacyc']: continue
				else:
					#rxnbiocyc
					for rxnRheaId in dbToRheaRxns['metacyc'][rxnbiocyc]:
						rheaRxnsInBiocycSpec.add(rxnRheaId)
						biocycAllRxns.add(rxnRheaId)
						
			rheaMetsInBiocycSpec = set()
			for met in allSpecBioCyc[organismNameBiocyc]['exchangedMets']|allSpecBioCyc[organismNameBiocyc]['demandedMets']:
				if met not in dbToRheaMets['metacyc']: continue
				else:
					for metRheaId in dbToRheaMets['metacyc'][met]:
						rheaMetsInBiocycSpec.add(metRheaId)
						biocycAllMets.add(metRheaId)
						
			biocycMetsPerTax[taxLavel].append(rheaMetsInBiocycSpec)
			biocycRxnsPerTax[taxLavel].append(rheaRxnsInBiocycSpec)
			
		if not biocycRxnsPerTax[taxLavel]: del biocycRxnsPerTax[taxLavel]
		
	biocycRxnsFreq = dict()
	for rxnId in biocycAllRxns:
		freq = list()
		for taxLevel in biocycRxnsPerTax:
			count = 0
			for specList in biocycRxnsPerTax[taxLevel]:
				if rxnId in specList: count += 1
			
			freq.append(count/len(biocycRxnsPerTax[taxLevel]))
		biocycRxnsFreq[rxnId] = sum(freq)/len(freq)
		
	biocycMetsFreq = dict()
	for metId in biocycAllMets:
		freq = list()
		for taxLevel in biocycMetsPerTax:
			count = 0
			for specList in biocycMetsPerTax[taxLevel]:
				if metId in specList: count += 1
			
			if count == 0: freq.append(0)
			else: freq.append(count/len(biocycMetsPerTax[taxLevel]))
		biocycMetsFreq[metId] = sum(freq)/len(freq)
		
		
	
	biocycRxnsFreqSlots = {'lower than 0.5': set(), '0.5 to 0.59': set(), '0.6 to 0.69': set(), '0.7 to 0.79': set(), '0.8 to 0.89': set(), '0.9 to 0.99': set(), '1': set()}
	for eachRxn in biocycRxnsFreq:
		if biocycRxnsFreq[eachRxn] == 1: biocycRxnsFreqSlots['1'].add(eachRxn)
		elif biocycRxnsFreq[eachRxn] >= 0.9:  biocycRxnsFreqSlots['0.9 to 0.99'].add(eachRxn)
		elif biocycRxnsFreq[eachRxn] >= 0.8:  biocycRxnsFreqSlots['0.8 to 0.89'].add(eachRxn)
		elif biocycRxnsFreq[eachRxn] >= 0.7:  biocycRxnsFreqSlots['0.7 to 0.79'].add(eachRxn)
		elif biocycRxnsFreq[eachRxn] >= 0.6:  biocycRxnsFreqSlots['0.6 to 0.69'].add(eachRxn)
		elif biocycRxnsFreq[eachRxn] >= 0.5:  biocycRxnsFreqSlots['0.5 to 0.59'].add(eachRxn)
		else: biocycRxnsFreqSlots['lower than 0.5'].add(eachRxn)
		
	biocycMetsFreqSlots = {'lower than 0.5': set(), '0.5 to 0.59': set(), '0.6 to 0.69': set(), '0.7 to 0.79': set(), '0.8 to 0.89': set(), '0.9 to 0.99': set(), '1': set()}
	for eachMet in biocycMetsFreq:
		if biocycMetsFreq[eachMet] == 1: biocycMetsFreqSlots['1'].add(eachMet)
		elif biocycMetsFreq[eachMet] >= 0.9:  biocycMetsFreqSlots['0.9 to 0.99'].add(eachMet)
		elif biocycMetsFreq[eachMet] >= 0.8:  biocycMetsFreqSlots['0.8 to 0.89'].add(eachMet)
		elif biocycMetsFreq[eachMet] >= 0.7:  biocycMetsFreqSlots['0.7 to 0.79'].add(eachMet)
		elif biocycMetsFreq[eachMet] >= 0.6:  biocycMetsFreqSlots['0.6 to 0.69'].add(eachMet)
		elif biocycMetsFreq[eachMet] >= 0.5:  biocycMetsFreqSlots['0.5 to 0.59'].add(eachMet)
		else: biocycMetsFreqSlots['lower than 0.5'].add(eachMet)
			
	return biocycRxnsFreqSlots, biocycMetsFreqSlots, pathwaysBiocyc
	
def findSimilarSpec(specName):
	#encontra especies no KEGG Organisms e Biocyc que sao taxonomicamente proximos ao organismo de interesse
	
	pickle_file_TaxoNcbi = os.path.join(project_dir, 'data/generated', 'spectAndTaxoNcbi.pickle')
	pickle_file_allSpecBioCyc = os.path.join(project_dir, 'data/generated', 'allSpecBioCyc.pickle')
	pickle_file_allSpecKegg = os.path.join(project_dir, 'data/generated', 'allSpecKegg.pickle')
	with open(pickle_file_TaxoNcbi, 'rb') as f:
		spectAndTaxoNcbi = pickle.load(f)
		
	specName = specName.replace('-', ' ').replace('_', ' ').replace('  ', ' ')
	if ' ' in specName: options = [specName.lower(), specName.split(' ')[0].lower(), specName.split(' ')[1].lower()]
	else: options = [specName.lower()]


	#find taxonomy classification of modelled organism
	#se encontrou uma especie e sua taxonomia - o nome do arquivo eh provavelmente de uma especia
	taxo = 0
	taxoOfTarget = 0
	for eachOption in options:
		for specKey in spectAndTaxoNcbi:
			if re.search(r'\b'+ eachOption +r'\b', specKey):
				taxo = spectAndTaxoNcbi[specKey]
				if eachOption not in taxo: taxo.append(eachOption)
				break
		if taxo: break

	taxonomicalLevelPerDb = {'Kegg': dict(), 'Biocyc': dict()}
	
	if not taxo:
		return taxonomicalLevelPerDb, list()
	
	
	if os.path.exists(pickle_file_allSpecBioCyc) and os.path.exists(pickle_file_allSpecKegg):
	
		with open(pickle_file_allSpecBioCyc, 'rb') as f:
			allSpecBioCyc = pickle.load(f)
			
		with open(pickle_file_allSpecKegg, 'rb') as f:
			allSpecKegg = pickle.load(f)
	
		#encontra especies baseado no primeiro nome de especias filogeneticamente proximas
		taxo.reverse()
		
		for tax in taxo:

			specSet = set()
			for specKey in spectAndTaxoNcbi:
				if tax in spectAndTaxoNcbi[specKey]:
					specSet.add(specKey.split(' ')[0])
					
			if specSet:
				keggCont, biocycCount = 0, 0
				for specName in specSet:
					taxonomicalLevelPerDb['Kegg'][specName] = list()
					taxonomicalLevelPerDb['Biocyc'][specName] = list()
					for specBiocyc in allSpecBioCyc:
						if specBiocyc.lower().startswith(specName):
							taxonomicalLevelPerDb['Biocyc'][specName].append(specBiocyc)
							biocycCount += 1
					for specBiocyc in allSpecKegg:
						if specBiocyc.lower().startswith(specName):
							taxonomicalLevelPerDb['Kegg'][specName].append(specBiocyc)
							keggCont += 1
							
			if keggCont > 15 and biocycCount > 15 and keggCont + biocycCount > 50: 
				taxo.reverse()
				taxoOfTarget = taxo[:taxo.index(tax)+1]
				break
		if 'cellular organisms' in taxoOfTarget: taxoOfTarget.remove('cellular organisms')
	else:
		if 'cellular organisms' in taxo: taxo.remove('cellular organisms')
		taxoOfTarget = taxo
				
	return taxonomicalLevelPerDb, taxoOfTarget
					

def taxonomyBasedConstraints(fileNameAndPath, cobraModel):

	#fileNameAndPath = 'C:/phd/Bangor/Aziz/Isolated/A25/Vagococcus acidifermentans.fa'
	print('fileNameAndPath = ' + str(fileNameAndPath))
	fileName = os.path.basename(fileNameAndPath)
	specName = os.path.splitext(fileName)[0].replace(' - ORFs','')
	if fileNameAndPath == fileName: filePath = os.getcwd()
	else: filePath = fileNameAndPath.replace(fileName,'')

	
	rxnsInSoftSyn = set()
	#organismNameKeggList, organismNameBiocycList, taxo = findSimilarSpec(specName) #{'Escherichia coli': {'Kegg': ['Escherichia coli ABU 83972 (2012)'], 'Biocyc': ['Escherichia coli K-12 substr.MG1655 (Tier 1) | 385 pathways; 1989 reactions']}}
	taxonomicalLevelPerDb, taxoOfTarget = findSimilarSpec(specName) #{'Kegg': {'tetragenococcus': ['Tetragenococcus koreensis (2018)', 'Tetragenococcus halophilus YJ1 (2019)'}, 'Biocyc': {'tetragenococcus': ['Tetragenococcus solitarius NBRC 100494 (Tier 3) | 163 pathways; 766 reactions']}}
	
	dbToRheaRxns, dbToRheaMets, rheaWithSameSyn = findAllMetsAndRxnsTranslations(cobraModel)
	
	#check if 'allSpecBioCyc' exists. 
	try:
		pickle_file_path = os.path.join(project_dir, 'data/generated', 'allSpecBioCyc.pickle')
		with open(pickle_file_path, 'rb') as f:
			allSpecBioCyc = pickle.load(f)
			
		pickle_file_path = os.path.join(project_dir, 'data/generated', 'allSpecKegg.pickle')
		with open(pickle_file_path, 'rb') as f:
			allSpecKegg = pickle.load(f)
			
		pickle_file_path = os.path.join(project_dir, 'data/generated', 'biocycPathways.pickle')
		with open(pickle_file_path, 'rb') as f:
			biocycPathways = pickle.load(f)
	except:
		soft_constraints = {'kegg': {'lower than 0.5': set(), '0.5 to 0.59': set(), '0.6 to 0.69': set(), '0.7 to 0.79': set(), '0.8 to 0.89': set(), '0.9 to 0.99': set(), '1': set(), '1 intersection':set(),'pathways':set()}, 'metacyc': {'lower than 0.5': set(), '0.5 to 0.59': set(), '0.6 to 0.69': set(), '0.7 to 0.79': set(), '0.8 to 0.89': set(), '0.9 to 0.99': set(), '1 intersection':set(),'1':set(),'pathways':set()}, 'metacyc met': {'lower than 0.5': set(), '0.5 to 0.59': set(), '0.6 to 0.69': set(), '0.7 to 0.79': set(), '0.8 to 0.89': set(), '0.9 to 0.99': set(), '1 intersection':set(),'1':set(),'pathways':set()}, 'User constraint':{'To Include': set(), 'To avoid': set()}}
		return soft_constraints, rxnsInSoftSyn, taxoOfTarget, rheaWithSameSyn
		

	
	keggRxnsFreqSlots, pathwaysKegg = findKeggRxnsInReferenceOrgs(taxonomicalLevelPerDb, dbToRheaRxns, allSpecKegg)
	
	biocycRxnsFreqSlots, biocycMetsFreqSlots, pathwaysBiocyc = findBiocycRxnsInReferenceOrgs(taxonomicalLevelPerDb, dbToRheaRxns, dbToRheaMets, allSpecBioCyc, biocycPathways)

	#organiza os resultados
	keggRxnsFreqSlots['1 intersection'] = keggRxnsFreqSlots['1'].intersection(biocycRxnsFreqSlots['1'])
	keggRxnsFreqSlots['1'] = keggRxnsFreqSlots['1'] - keggRxnsFreqSlots['1 intersection']
	keggRxnsFreqSlots['pathways'] = pathwaysKegg
	biocycRxnsFreqSlots['1 intersection'] = keggRxnsFreqSlots['1 intersection'] #a intersecao eh igual nos dois
	biocycRxnsFreqSlots['1'] = biocycRxnsFreqSlots['1'] - biocycRxnsFreqSlots['1 intersection']
	biocycRxnsFreqSlots['pathways'] = pathwaysBiocyc
	
	
	soft_constraints = {'kegg': keggRxnsFreqSlots, 'metacyc': biocycRxnsFreqSlots, 'metacyc met': biocycMetsFreqSlots, 'User constraint':{'To Include': set(), 'To avoid': set()}}
	#salvando arquivo ara usar como constraint
	softConstraintsTable = ''
	for eachKey in soft_constraints:
		for freqSlot in soft_constraints[eachKey]:
			for rheaId in soft_constraints[eachKey][freqSlot]:
				if rheaId.startswith('R_'): softConstraintsTable += rheaId + '\t' + eachKey + '\t' + freqSlot + '\n'
				else: softConstraintsTable += 'R_'+rheaId + '\t' + eachKey + '\t' + freqSlot + '\n'
				
	softConstraintsTable = softConstraintsTable[:-1]		

	fileS = open(filePath + '\\' + specName + ' - TaxoConstraints.tsv', 'w')
	fileS.write(softConstraintsTable)
	fileS.close()
	

	#cria set com todas as reacoes encontradas no modelo de acordo com o biocyc e kegg, assim como as suas traducoes
	
	rxnsInSoft = soft_constraints['kegg']['1 intersection']|soft_constraints['kegg']['1']|soft_constraints['kegg']['0.9 to 0.99']|soft_constraints['kegg']['0.8 to 0.89']|soft_constraints['kegg']['0.7 to 0.79']|soft_constraints['kegg']['0.6 to 0.69']|soft_constraints['kegg']['0.5 to 0.59']|soft_constraints['metacyc']['1 intersection']|soft_constraints['metacyc']['1']|soft_constraints['metacyc']['0.9 to 0.99']|soft_constraints['metacyc']['0.8 to 0.89']|soft_constraints['metacyc']['0.7 to 0.79']|soft_constraints['metacyc']['0.6 to 0.69']|soft_constraints['metacyc']['0.5 to 0.59']
	
	
	for rxn in cobraModel.reactions:
		if rxn.id.replace('_forwardTemp','').replace('_reverseTemp','').replace('_bidirectionalCopy','') in rxnsInSoft:
			rxnsInSoftSyn.add(rxn.id.replace('_forwardTemp','').replace('_reverseTemp','').replace('_bidirectionalCopy',''))
			if 'rhea' not in rxn.annotation: continue
			for rxnId in rxn.annotation['rhea']:
				rxnsInSoftSyn.add(rxnId)
	
	

	if softConstraintsTable:
		return soft_constraints, rxnsInSoftSyn, taxoOfTarget, rheaWithSameSyn
	else: 
		sys.exit("It was not possible to identify closely related organisms. The file name should be the organism name (e.g. Escherichia coli.fasta), or the taxonomical classification (e.g. 'Bacteroides.fasta' or 'Bacteroidia.fasta' or 'Cyanobacteria.fasta')")
		
