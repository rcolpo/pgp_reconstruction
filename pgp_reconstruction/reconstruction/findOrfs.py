#check if file is DNA or amino acid. If DNA, find ORFs with prodigal and translate sequences.


import os
import re
from pgp_reconstruction import project_dir
				

def writeFasta(inputFileName):
	#receive as input genbank file, extract proteins and write fasta
	
	f = open(inputFileName, 'r')
	assemblyFileLines = f.readlines()
	f.close()

	sequencesDict = dict()
	copyingTranslation = 0
	copyingProduct = 0
	translation = ''
	product = ''
	gene = ''
	ecNumber = ''
	for line in assemblyFileLines:
	
			if line.strip().startswith('/product="'):
				copyingProduct = 1
				lineT = line.strip().replace('/product="','')
				product += lineT.replace('"','')
				
				if lineT[-1] == '"':
					copyingProduct = 0
					
			elif copyingProduct == 1:
				lineT = line.strip()
				product += lineT.replace('"','')
				if lineT[-1] == '"':
					copyingProduct = 0
					
	
			if line.strip().startswith('/protein_id="'):
				protein_id = line.strip().replace('/protein_id="','')
				protein_id = protein_id.split(':')[-1].replace('"','')
				translation = ''
				
			if line.strip().startswith('/locus_tag="'):
				protein_id = line.strip().replace('/locus_tag="','')
				protein_id = protein_id.split(':')[-1].replace('"','')
				translation = ''
				
			if line.strip().startswith('/gene="'):
				gene = line.strip().replace('/gene="','').replace('"','')
				
			if line.strip().startswith('/EC_number="'):
				ecNumber = line.strip().replace('/EC_number="','').replace('"','')
				
			if line.strip().startswith('/translation="'):
				copyingTranslation = 1
				lineT = line.strip().replace('/translation="','')
				translation += lineT.replace('"','')
				
				if lineT[-1] == '"':
					copyingTranslation = 0
					sequencesDict[protein_id] = dict()
					sequencesDict[protein_id]['seq'] = translation
					sequencesDict[protein_id]['product'] = product
					sequencesDict[protein_id]['gene'] = gene
					sequencesDict[protein_id]['ecNumber'] = ecNumber
					translation = ''
					product = ''
					gene = ''
					ecNumber = ''

				
			elif copyingTranslation == 1:

				lineT = line.strip()
				translation += lineT.replace('"','')
			
				if lineT[-1] == '"':
					copyingTranslation = 0
					sequencesDict[protein_id] = dict()
					sequencesDict[protein_id]['seq'] = translation
					sequencesDict[protein_id]['product'] = product
					sequencesDict[protein_id]['gene'] = gene
					sequencesDict[protein_id]['ecNumber'] = ecNumber
					translation = ''
					product = ''
					gene = ''
					ecNumber = ''

	toWrite = ''
	for protein_id in sequencesDict:
	
		toWrite += '>'+protein_id
		if sequencesDict[protein_id]['product']:
			toWrite += ' ' + sequencesDict[protein_id]['product']
		if sequencesDict[protein_id]['gene']:
			toWrite += ' GN=' + sequencesDict[protein_id]['gene']
		if sequencesDict[protein_id]['ecNumber']:
			toWrite += ' EC=' + sequencesDict[protein_id]['ecNumber']		
		toWrite += '\n'
			
		chunks = [sequencesDict[protein_id]['seq'][i:i+60] for i in range(0, len(sequencesDict[protein_id]['seq']), 60)]
		for chunk in chunks:
			toWrite += chunk + '\n'
		
	toWrite = toWrite[:-1]
					
	fileName = os.path.splitext(os.path.basename(inputFileName))[0]
	
	fileLocation = inputFileName.split(fileName)[0]
	if fileLocation: fileLocation += '/'
	
	inputFileNameToWrite = fileLocation + fileName + '.fasta'
					
	f = open(inputFileNameToWrite, 'w')
	f.write(toWrite)
	f.close()
	
	return inputFileNameToWrite
				
def findOrfs(inputFileName):

	#identifica se eh DNA ou proteina. Se for DNA, identifica ORFS e traduz para proteina.
	
	#check if file is from genbank. if genbank, create fasta file for alingment
	if inputFileName.endswith('.gb') or inputFileName.endswith('.gbk') or inputFileName.endswith('.genbank') or inputFileName.endswith('.gbff'):
		inputFileName = writeFasta(inputFileName)
	
	dnaLetters = {'A', 'C', 'T', 'G', 'N'}
	aminoLetters = {'G', 'C', 'R', 'I', 'Y', 'D', 'Q', 'L', 'P', 'A', 'S', 'N', 'E', 'H', 'K', 'W', 'V', 'M', 'F', 'T'}
	dnaCount = 0
	aminoCount = 0
	
	#check if file is DNA or proteins
	f = open(inputFileName,'r')
	for line in f.readlines():
		if line and line[0] != '>': 
			for letter in line:
				if letter in dnaLetters:
					dnaCount += 1
				if letter in aminoLetters:
					aminoCount += 1
			if aminoCount != 0:
				break
	f.close()
	
	#if it is dna sequence, then uses prodigal to identify ORFs
	prodigal = 0
	if dnaCount == aminoCount:
		fileName = os.path.basename(inputFileName)
		filePath = inputFileName.replace(fileName,'')
		fileNameAndExtension = os.path.splitext(fileName)
	
		prodigalFile = ''
		for fileName in os.path.join(project_dir, 'dependencies'):
			if fileName in 'prodigal': prodigalFile = fileName
		if not prodigalFile: 
			sys.exit('Could not locate prodigal file to identify ORF in input genome.')
			
		prodigalPath = os.path.join(project_dir, 'dependencies', prodigalFile)
	
		os.system(prodigalPath + ' -i "' + inputFileName + '" -a "' + filePath + fileNameAndExtension[0] + ' - ORFs.faa"')
		prodigal = 1
	
		inputFileNameNew = filePath + fileNameAndExtension[0] + ' - ORFs.faa'
	else: inputFileNameNew = inputFileName
	

	#tenta idendificar gene and protein name
	geneAndProteinNamePerSeqId = dict()
	f = open(inputFileName,'r')
	for line in f.readlines():
		if line and line[0] == '>':
			line = line.strip()
			proteinId = line.split(' ')[0][1:]
			geneAndProteinNamePerSeqId[proteinId] = {'gene':'', 'protein name 1': '', 'protein name 2': ''}
			if prodigal == 0:
				if 'hypothetical protein' in line: continue
				
				if '[gene=' in line:
					#se for genBank annotation
					for eachSection in line.split(']'):
						if '[gene=' not in eachSection: continue
						a = re.search('(?<=\[gene\=).*', eachSection)
						geneAndProteinNamePerSeqId[proteinId]['gene'] = a.group().lower()
						
					for eachSection in line.split(']'):
						if '[protein=' not in eachSection: continue
						a = re.search('(?<=\[protein\=).*', eachSection)
						geneAndProteinNamePerSeqId[proteinId]['protein name 1'] = a.group().lower()
						
				elif ' GN=' in line:
					#se for genBank annotation
					geneAndProteinNamePerSeqId[proteinId]['gene'] = line.split(' GN=')[1].split(' ')[0].lower()
					geneAndProteinNamePerSeqId[proteinId]['protein name 1'] = line[1:].split(' GN=')[0].lower()[len(proteinId)+1:]
				else:
					#se for prokka annotation
					title = line.replace('>'+proteinId+' ', '')
					geneAndProteinNamePerSeqId[proteinId]['gene'] = title.split(' ')[-1].lower()
					geneAndProteinNamePerSeqId[proteinId]['protein name 1'] = title.replace(' '+title.split(' ')[-1], '').lower()
					geneAndProteinNamePerSeqId[proteinId]['protein name 2'] = title.lower()
					
	f.close()
	
	return inputFileNameNew, geneAndProteinNamePerSeqId