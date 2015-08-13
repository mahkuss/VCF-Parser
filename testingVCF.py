import numpy as np
import pandas as pd
import argparse
import sys
import os.path
import VCFparser
from collections import OrderedDict
from datetime import datetime

parser = argparse.ArgumentParser(description='Calculate pairwise concordance rates between variants in two VCF files.')
parser.add_argument('experimental_vcf', help='VCF file containing samples you want to check')
parser.add_argument('source_vcf', help='VCF file containing reference samples you want to compare to')
parser.add_argument('-o', '--output', help='Directory to write output to',
					default=os.getcwd())

args = parser.parse_args()
print('output directory: {}'.format(args.output))
try:
	log_fp = os.path.join( args.output, 'run.log' )
	log = open(log_fp, 'w')
	log.write( 'testingVCF.py Log File\n\n' )
	log.flush()
except Exception as e:
	sys.exit('Error: Could not write to specified output directory.')


VERBOSE = False
DEBUG = False
DUMP_CHROM = 'chr1'

if DEBUG:
	DEBUG_DUMP = '/home/rob/projects/c9orf72_ko_mouse/scratch/DEBUGchr1_testingVCF.out'
	debugDump = open(DEBUG_DUMP, 'w')


newVCF = VCFparser.VCF( args.experimental_vcf )
newVCF2 = VCFparser.VCF( args.source_vcf )


chromosomeOutput = []

wtf = [] # debugging list to collect oddities that pop up along the way

class ChromosomeMatchData:
	'''
	Class to hold and manipulate allele match data on a single chromosome for all samples in two separate VCF files.
	'''

	def __init__(self, chrom, expVCF, sourceVCF):
		'''
		Initializer function to create container object for all data related to VCF entries from the same chromosome
		'''

		#heterozygous_position_matches = { sample:[0] * (newVCF2.num_samples + 1) for sample in newVCF.samples }
		#homozygous_position_matches = { sample:[0] * (newVCF2.num_samples + 1) for sample in newVCF.samples }
		self.chromosome = chrom
		self.numRecordsInExpFile = 0 # total number of records in expFile corresponding to this chromosome
		self.numRecordsInSourceFile = 0 # total number of records in sourceFile corresponding to this chromosome
		self.numMatchingRecords = 0 # total number of records found in both input files

		self.heterozygousMatches = [0] * newVCF.num_samples # number of heterozygous variants in expFile 
															# for which there is an entry in sourceFile
		self.homozygousMatches = [0] * newVCF.num_samples 	# number of homozygous variants in expFile
															# for which there is an entry in sourceFile
		self.heterozygousSameAsRef = [0] * newVCF.num_samples  # number of self.heterozygousMatches for which
																# there isn't actually a variant in a sample
																# i.e., pulled along for the ride by a variant
																# called in another sample
		self.homozygousSameAsRef = [0] * newVCF.num_samples 	# number of self.homozygousMatches for which
																# there isn't actually a variant in a given sample
																# i.e., pulled along for the ride by a variant
																# called in another sample
		self.heterozygousNoCall = [0] * newVCF.num_samples  # number of self.heterozygousMatches for which
															# there is no base call information for a sample
		self.homozygousNoCall = [0] * newVCF.num_samples 	# number of self.homozygousMatches for which
															# there is no base call information for a sample

		self.compoundHeterozygous = [0] * newVCF.num_samples 	# count of compound heterozygous variants identified
																# in records with a matching record in both files

		self.sourceBaseCalled = { 
								  SNPtype: {sample:np.array( [0] * (newVCF2.num_samples + 1) ) for sample in newVCF.samples} \
								  for SNPtype in ['homozygous', 'heterozygous'] 
								 } # number of times a base was called for each sample in records that match across
								   # both input files; the + 1 is because the number of sites that don't match any
								   # sourceFile sample genome is stored in the -1 index of the homo or het _position_matches
								   # dictionary and therefore the -1 index of this array will be updated with the appropriate
								   # denominator for that calculation when the chromosome is summarized by call to self.sumarizeOutput()

		self.sourceBaseDiffFromRef = { 
								  		SNPtype: {sample:np.array( [0] * newVCF2.num_samples) for sample in newVCF.samples} \
								  		for SNPtype in ['homozygous', 'heterozygous'] 
								 	 } 	# number of times the base called in a sample for a record matching an expFile entry is
										# actually a variant relative to the reference

		# initialize OrderedDict objects to store the number of matching variants between any expFile sample and
		# any sourceFile sample; will use the last position of each list to record number of times there was no
		# sourceFile sample for which a matching base was called in an expFile sample 
		self.heterozygous_position_matches = OrderedDict()
		self.homozygous_position_matches = OrderedDict()
		for sample in expVCF.samples:
			self.heterozygous_position_matches[sample] = [0] * (sourceVCF.num_samples + 1)
			self.homozygous_position_matches[sample] = [0] * (sourceVCF.num_samples + 1)
			# num_sammples + 1 so that -1 index will hold number of variants that did not match any source genotype


	def _scanSourceFile(self, expBase, sample, sourceBaseCalls, ref, SNPtype, boolUpdateSourceRef):
		'''
		Private class method to scan through sourceFile base calls for a given expFile
		sample SNP and update the relevant position_matches dictionary that holds
		counts for that class of SNP (het vs homo, passed in through the argument
		SNPtype). This function is called exclusively by countMatches().
		'''

		updateSourceCalled = [0] * (len(sourceBaseCalls) + 1) 	# see comment on self.sourceBaseCalled in __init__ for
																# explanation of the +1 (lines 81 -88)
		updateSourceRefList = [0] * len(sourceBaseCalls)

		# define the index into chromData position matches (het vs homo) based on value of SNPtype
		if SNPtype == 'heterozygous':
			positionMatchDict = self.heterozygous_position_matches
		elif SNPtype == 'homozygous':
			positionMatchDict = self.homozygous_position_matches

		matchFound = False
		for ind_j in range( len( sourceBaseCalls ) ):
			for sourceBase in sourceBaseCalls[ ind_j ]: # check against each allele (homozygous matches will count as 2)

				# update counts for base being called at position and base being different than reference genome
				updateSourceCalled[ ind_j ] += int(sourceBase in ['A', 'C', 'T', 'G'])
				updateSourceRefList[ ind_j ] += int(sourceBase in ['A', 'C', 'T', 'G'] and sourceBase != ref)
				
				if sourceBase not in  ['A', 'C', 'T', 'G', 'ND']:
				 	wtf.append( sourceBase ) # if it's not a base then wtf is it?

				# check for a match at this allele
				if sourceBase == expBase:
					matchFound = True
					positionMatchDict[ sample ][ ind_j ] += 1

		if not matchFound: # if matchFound has not tripped to True then the variant in expFile is unique from all sourceFile genomes
			positionMatchDict[ sample ][-1] += 1

		self.sourceBaseCalled[ SNPtype ][ sample ] += np.array( updateSourceCalled )
		self.sourceBaseDiffFromRef[ SNPtype ][ sample ] += np.array( updateSourceRefList )


	def countMatches(self, expBaseCalls, sourceBaseCalls, ref):
		'''
		Class method to update position match tabulation between experimental file and 
		source file at a given position, segregated by homozygous and heterozygous SNPs. 
		Takes the base calls for a given position from each file as arguments. For each 
		heterozygous or homozygous base call in expBaseCalls it will add one to the corresponding 
		position_match dicitonary index for	any sourceBaseCalls samples that match. Therefore, 
		under current implementation, a homozygous match will add 1 and a heterozygous match
		will add 1 *to each source sample that matches EACH base*. Note that the code will check 
		that expBaseCall is actually different than the genome reference base at that position to 
		prevent counting of non-SNPs that hitch a ride with true SNPs in other samples in a multi-
		sample VCF file. If the expBaseCalls base is not found in any sample in sourceBaseCalls, 1 
		is added to the total in the -1 index of the relevant position_matches dictionary.
		'''


		## TO DO: 	
		##			-consider adding a check that a variant is truly a SNP in case polluted VCF fed into program

		for ind_i, expSample in enumerate( self.homozygous_position_matches.keys() ):
			if expBaseCalls[ ind_i ][0] == expBaseCalls[ ind_i ][1]: # test for homozygosity
				expBase = expBaseCalls[ ind_i ][0]
				if expBase == ref: # then it's not really a SNP -- it got pulled in by a variant in another sample
					self.homozygousSameAsRef[ ind_i ] += 1
					continue

				elif expBase == 'ND': # it's not called
				# elif expBaseCalls[ ind_i ][ 0 ] not in ['A', 'C', 'T', 'G']: # it's not called
					self.homozygousNoCall[ ind_i ] += 1 
				
				else:
					self.homozygousMatches[ ind_i ] += 1
					self._scanSourceFile( expBase, expSample, sourceBaseCalls, ref,
									'homozygous', ind_i ) # tabulates source file when last argument = 0
														# -3 shifts index to check when MouseWild is homo variant


			else: # now handle the heterozygous SNP case
				if expBaseCalls[ ind_i ][0] != ref: # test for compound heterozygote
					self.compoundHeterozygous[ ind_i ] += 1
					basesToCheck = [ base for base in expBaseCalls[ ind_i ] ]
				else:
					basesToCheck =[ expBaseCalls[ ind_i ][1] ]

				for expBase in basesToCheck:
					if expBase == ref: # then it's not really a SNP -- it got pulled in by a variant in another sample
						self.heterozygousSameAsRef[ ind_i ] += 1
						continue
					elif expBase == 'ND':
						self.heterozygousNoCall[ ind_i ] += 1 
						# wtf.append(expBase) # if it's not a base then wtf is it?
					else:
						self.heterozygousMatches[ ind_i ] += 1
						self._scanSourceFile( expBase, expSample, sourceBaseCalls, ref,
										'heterozygous', ind_i + 1 ) # ind + 1 effectively turns off checking for het variants



	def summarizeOutput(self): # consider changing this to NOT add the pandas dataframe as a class attribute but rather
								# simply generate the summary on the fly when called
		'''
		Function to summarize chromosome data container and pretty print output at end of run.
		'''

		logprint('\tGenerating summary for %s' %self.chromosome)
		
		expSamples = newVCF.samples
		sourceSamples = [ sample for sample in newVCF2.samples ]
		sourceSamples.append('No_Match')
		
		### THESE TABLES CURRENTLY CALCULATED BUT NOT USED ANYWHERE ####
		homoSourceCalled = pd.DataFrame(self.sourceBaseCalled[ 'homozygous' ], index=sourceSamples)
		hetSourceCalled = pd.DataFrame(self.sourceBaseCalled[ 'heterozygous' ], index=sourceSamples)

		homoSourceDiffFromRef = pd.DataFrame(self.sourceBaseDiffFromRef[ 'homozygous' ], index=sourceSamples[:-1])
		hetSourceDiffFromRef = pd.DataFrame(self.sourceBaseDiffFromRef[ 'heterozygous' ], index=sourceSamples[:-1])
		##### END OF TABLES CALCULATED BUT NOT YET USED #####

		homoSummary = pd.DataFrame(index=sourceSamples)
		hetSummary = pd.DataFrame(index=sourceSamples)

		homoNoMatch = []
		hetNoMatch = []
		homoNumChecked = []
		hetNumChecked = []

		# update the -1 index of self.sourceBaseCalled['heterozygous'] and ['homozygous'] to hold the number of
		# het and homozygous matches (i.e., self.homozygousMatches, self.heterozygousMatches); this is the appropriate
		# denominator for the calculation of percent of variants not matching any sourceFile sample genome
		for ind_i, expSample in enumerate(expSamples):
			self.sourceBaseCalled['homozygous'][ expSample ][ -1 ] = self.homozygousMatches[ ind_i ]
			self.sourceBaseCalled['heterozygous'][ expSample ][ -1 ] = self.heterozygousMatches[ ind_i ]

		for expSample in expSamples:
			matches = self.homozygous_position_matches[ expSample ]
			percentMatch = np.array( matches ) / self.sourceBaseCalled[ 'homozygous' ][ expSample ].astype(float) * 100
			homoSummary[ expSample ] = percentMatch
			homoNoMatch.append( float(self.homozygous_position_matches[ expSample ][ -1 ]) )
			homoNumChecked.append( np.array(self.homozygous_position_matches[ expSample ]).sum() )
			
			matches = self.heterozygous_position_matches[ expSample ]
			percentMatch = np.array( matches ) / self.sourceBaseCalled[ 'heterozygous' ][ expSample ].astype(float) * 100
			hetSummary[ expSample ] = percentMatch
			hetNoMatch.append( float(self.heterozygous_position_matches[ expSample ][ -1 ]) )
			hetNumChecked.append( np.array(self.heterozygous_position_matches[ expSample ]).sum() )

		self.homoSummary = homoSummary
		self.hetSummary = hetSummary
		self.homoNoMatch = np.array(homoNoMatch) / np.array(self.homozygousMatches)
		self.hetNoMatch = np.array(hetNoMatch) / np.array(self.heterozygousMatches)

		logprint('\tSummary for %s complete\n' %self.chromosome)


	def writeCSV(self):
		'''
		Class method to write results to a csv file.
		'''
		if type(self.homoSummary) == pd.core.frame.DataFrame:
			logprint('\tWriting homozygous percent matches for %s to csv' %self.chromosome)
			self.homoSummary.to_csv( os.path.join( args.output, self.chromosome + '_homozygous.csv' ) )
		if type(self.hetSummary) == pd.core.frame.DataFrame:
			logprint('\tWriting heterozygous percent matches for %s to csv' %self.chromosome)
			self.hetSummary.to_csv( os.path.join( args.output, self.chromosome + '_heterozygous.csv' ) )



## END CLASS DEFINITION ##


def skipHeader(vcf, vcfIterator):
	'''
	Simple helper function to advance past the header lines in the 2 VCF files.
	'''
	lineNum = 0
	while lineNum < vcf.header_line_count:
		junk = vcfIterator.readline()
		lineNum += 1


def homozygousMatchTable(chromosomeOutput):
	return( pd.DataFrame( { chrom.chromosome:pd.Series(chrom.homozygousMatches, index=newVCF.samples) for chrom in chromosomeOutput } ) )

def heterozygousMatchTable(chromosomeOutput):
	return( pd.DataFrame( { chrom.chromosome:pd.Series(chrom.heterozygousMatches, index=newVCF.samples) for chrom in chromosomeOutput } ) )

def homozygousSameAsRefTable(chromosomeOutput):
	return( pd.DataFrame( { chrom.chromosome:pd.Series(chrom.homozygousSameAsRef, index=newVCF.samples) for chrom in chromosomeOutput } ) )

def heterozygousSameAsRefTable(chromosomeOutput):
	return( pd.DataFrame( { chrom.chromosome:pd.Series(chrom.heterozygousSameAsRef, index=newVCF.samples) for chrom in chromosomeOutput } ) )

def homozygousNoCallTable(chromosomeOutput):
	return( pd.DataFrame( { chrom.chromosome:pd.Series(chrom.homozygousNoCall, index=newVCF.samples) for chrom in chromosomeOutput } ) )

def heterozygousNoCallTable(chromosomeOutput):
	return( pd.DataFrame( { chrom.chromosome:pd.Series(chrom.heterozygousNoCall, index=newVCF.samples) for chrom in chromosomeOutput } ) )

def getChromosomeMatchSummary(chromosomeOutput):
	dat = { 
			chrom.chromosome: pd.Series( 
											{
												'ExpFile': chrom.numRecordsInExpFile,
												'SourceFile': chrom.numRecordsInSourceFile,
												'Matching': chrom.numMatchingRecords
											}
										) 
			for chrom in chromosomeOutput
			}

	return( pd.DataFrame( dat ) )			


def writeChromosomeSummaryTable( func, fname):
	'''
	Helper function to write csv files of chromosome-wide summary objects created by their 
	respective summary functions (e.g., homozygousMatchTable() ).
	'''

	logprint('\t' + fname + '.csv')
	newTable = func(chromosomeOutput)
	newTable.to_csv( os.path.join( args.output, fname + '.csv' ) )


def wrapup(finalLinecounts = False):
	'''
	Function to process last analyzed chromosome and generate chromosome-wide summary tables.
	Writes files to directory specified in global variable OUTDIR. To be called when either 
	expFile or sourceFile has been exhausted.
	'''

	if not sourceFileExhausted: # if souceFileExhausted, this chromosome has already been summarized
		currentChrom.summarizeOutput()
		currentChrom.writeCSV()
		# append the last currentChrom object to chromosomeOutput list
		chromosomeOutput.append( currentChrom )

	if finalLinecounts:
		logprint( '\n' + addTimestamp( 'Analysis complete. Summarizing results across all analyzed contigs.') )
		logprint( '\n' + addTimestamp( 'Generating summary tables across chromosomes.') )
		
		writeChromosomeSummaryTable( homozygousMatchTable, 'homozygousMatchingRecords' )
		writeChromosomeSummaryTable( heterozygousMatchTable, 'heterozygousMatchingRecords' )
		writeChromosomeSummaryTable( homozygousSameAsRefTable, 'homozygousSameAsRef' )
		writeChromosomeSummaryTable( heterozygousSameAsRefTable, 'heterozygousSameAsRef' )
		writeChromosomeSummaryTable( homozygousNoCallTable, 'homozygousNoCall' )
		writeChromosomeSummaryTable( heterozygousNoCallTable, 'heterozygousNoCall' )
		writeChromosomeSummaryTable( getChromosomeMatchSummary, 'TotalRecordsByChromosome' )
		
		allCountedRecords = getChromosomeMatchSummary( chromosomeOutput ).sum(1)
		expFileCountedRecords = allCountedRecords['ExpFile']
		sourceFileCountedRecords = allCountedRecords['SourceFile']
		expTotal = newVCF.header_line_count + expFileCountedRecords + skippedRecords
		sourceTotal = newVCF2.header_line_count + sourceFileCountedRecords + sourceSkippedRecords

		logprint( '\n' + addTimestamp( 'Number of lines analyzed:\n') )
		logprint('\tInput VCF 1:')
		logprint('\t\tHeader lines: %s' %newVCF.header_line_count)
		logprint('\t\tRecords on analyzed contigs: %s' %expFileCountedRecords)
		logprint('\t\tSkipped records: %s' %skippedRecords)
		logprint('\t\tTotal: %s\n' %expTotal)
		logprint('\tInput VCF 2:')
		logprint('\t\tHeader lines: %s' %newVCF2.header_line_count)
		logprint('\t\tRecords on analyzed contigs: %s' %sourceFileCountedRecords)
		logprint('\t\tSkipped records: %s' %sourceSkippedRecords)
		logprint('\t\tTotal: %s\n' %sourceTotal)

	if DEBUG:
		debugDump.close()	

def addTimestamp(string):
	'''
	Function to reformat a string with a timestamp in front prior to printing.
	'''

	timestamp = datetime.now().strftime('%m-%d-%Y %H:%M:%S')
	return('[%s] %s' %(timestamp, string))

def logprint(string):
	'''
	Print string and write it to the log file.
	'''

	print( string )
	log.write( string + '\n' )
	log.flush()

#########

count = 0
lines_checked = 0
skippedRecords = 0
sourceSkippedRecords = 0
sourceFileExhausted = False # for possible future use if decide to continue counting remaining
							# expFile lines after exhausting sourceFile

currentChrom = False

# if VCF v4.1 or higher, determine intersection of contigs for the two VCF files and
# restrict analysis to just these contigs
logprint( addTimestamp( 'Finding contigs common to both input VCF files...' ) )
expContigs = newVCF.contigs
sourceContigs = newVCF2.contigs

if ( expContigs == set() ) or ( sourceContigs == set() ):
	logprint('Unable to determine contigs from VCF headers.')
	logprint('Ensure that input VCF files have contig identifiers in headers.')
	f.close()
	sys.exit('Unable to determine contigs from VCF headers.')
else:
	chromosomesToSearch = expContigs.intersection( sourceContigs )
	if chromosomesToSearch == set([]):
		logprint('No matching contigs found between input VCF files.')
		logprint('Check that input VCF files share common contigs.')
		f.close()
		sys.exit('No matching contigs found between input VCF files; check that input VCF files share common contigs.')	
	logprint( '\nFile 1 Contigs: {}'.format( ', '.join( expContigs ) ) )
	logprint( '\nFile 2 Contigs: {}'.format( ', '.join( sourceContigs ) ) )
	logprint( '\n' + addTimestamp( 'Common contigs to be searched:' ) )
	for chrom in sorted(chromosomesToSearch):
		logprint('\t%s' %chrom)

with open(newVCF.filename, 'r') as expFile, open(newVCF2.filename, 'r') as sourceFile:
	skipHeader(newVCF, expFile)
	skipHeader(newVCF2, sourceFile)

	currentChrom = None # initalize value to None so it can be populated by first valid expFile entry in for loop below
	workingOnChrom = None # variable to hold name of current chromosome being checked; initialize here because used in logical below

	sourceEntry = VCFparser.VCFentry( sourceFile.readline() )
	logprint( '\n' + addTimestamp('Starting analysis ...') )

	for expLine in expFile:
		match = True # boolean flag will be tripped to false if line read from sourceFile is past current expFile position
		lines_checked += 1
		
		if lines_checked % 1000 == 0:
			logprint( addTimestamp('Processed %s lines' %lines_checked) )
		expEntry = VCFparser.VCFentry( expLine )
		
		if VERBOSE:
			print(expEntry.CHROM, expEntry.POS)

		if expEntry.CHROM not in chromosomesToSearch: # make sure variant is from a chromosome we need to be checking
			skippedRecords += 1
			continue

		if sourceFileExhausted:
			if sourceEntry.CHROM == currentChrom.chromosome: 
				currentChrom.numRecordsInExpFile += 1
				continue

			elif expEntry.CHROM != workingOnChrom:
				skippedRecords += 1
				continue

		if currentChrom == None:
			workingOnChrom = expEntry.CHROM
			currentChrom = ChromosomeMatchData( expEntry.CHROM, newVCF, newVCF2 )	
			if sourceEntry.CHROM == currentChrom.chromosome: 
				currentChrom.numRecordsInSourceFile += 1

		elif expEntry.CHROM != workingOnChrom:
			logprint( addTimestamp('Completed analyzing %s\n' %workingOnChrom) )
			chromosomeOutput.append( currentChrom )
			currentChrom.summarizeOutput()
			currentChrom.writeCSV()
			currentChrom = ChromosomeMatchData( expEntry.CHROM, newVCF, newVCF2 )
			workingOnChrom = currentChrom.chromosome
			logprint( '\n' + addTimestamp('Starting analysis of %s' %workingOnChrom) )

		currentChrom.numRecordsInExpFile += 1

		# find entry in source that matches current exp line; if no such line exists, catch this
		# and respond by moving on to the next chromosome
		while ( sourceEntry.CHROM, sourceEntry.POS ) != ( expEntry.CHROM, expEntry.POS ):
			if sourceFileExhausted:
				match = False
				break

			if VERBOSE:
				print('Exp - %s:%s' %(expEntry.CHROM, expEntry.POS))
				print('Source - %s:%s' %(sourceEntry.CHROM, sourceEntry.POS))
			
			if (sourceEntry.POS > expEntry.POS) and (sourceEntry.CHROM == expEntry.CHROM):
				match = False
				if VERBOSE:
					print('WENT TOO FAR ON CURRENT CHROMOSOME!')
				break

			elif ( sourceEntry.CHROM != expEntry.CHROM ):
			# have seen all sourceFile records for this chromosome without a match so need to
			# advance expFile to next chromosome; trip match to False and break from while loop
			# as this statement will continue to evaluate True until we have read in first record
			# from next chromosome in expFile through the outer for loop iterator
				if ( len(chromosomeOutput) > 0 ) and ( sourceEntry.CHROM == chromosomeOutput[ -1 ].chromosome ):
					# source file is on previous chromosome and needs to catch up
					if VERBOSE:
						print('SOURCE FILE CHROM IS BEHIND AND NEEDS TO ADVANCE')
					while sourceEntry.CHROM != expEntry.CHROM:
						sourceLine = sourceFile.readline()
						if not sourceLine:
							logprint( '\n' + addTimestamp('Reached end of source VCF\n') )
							wrapup()
							sourceFileExhausted = True
							match = False
							logprint( '\n' + addTimestamp('Counting remaining experimental VCF records ...\n'))
							break

							if DEBUG:
								debugDump.close()

						sourceEntry = VCFparser.VCFentry( sourceLine )
						currentChrom.numRecordsInSourceFile += 1

						if VERBOSE:
							print('NEXT ENTRY IS FROM CHROM %s' %sourceEntry.CHROM)
						# if sourceEntry.CHROM == None:
						# 	sys.exit('Exiting: End of source file or misformated line in source file.')

					# now we are at first entry of next chromosome; check for match and exit while loop
					if VERBOSE:
						print('CHROMOSOMES BACK IN SYNCH')
					if ( sourceEntry.CHROM, sourceEntry.POS ) != ( expEntry.CHROM, expEntry.POS ):
						if VERBOSE:
							print('BUT NOT A MATCH HERE')
						match = False
						break
				
				else: # source file has advanced to next chromosome
					if VERBOSE:
						print('SOURCE FILE CHROMOSOME IS AHEAD OF EXP FILE!')
					match = False
					break


			# this is the fall-through case; if we reach here then we are on the correct chromosome, but at
			# an earlier position than the current expFile entry; try to read the next line and continue
			else:
				try:
					sourceLine = sourceFile.readline()
					if not sourceLine:
						logprint( '\n' + addTimestamp('Reached end of source VCF\n') )
						wrapup()
						sourceFileExhausted = True
						match = False
						logprint( '\n' + addTimestamp('Counting remaining experimental VCF records ...\n'))
						break

					sourceEntry = VCFparser.VCFentry( sourceLine )
					currentChrom.numRecordsInSourceFile += 1

				except Exception as e:
					logprint(sourceLine)
					logprint( str(e) )
					wrapup()
					f.close()
					sys.exit()

		if VERBOSE:
			print('checking for match')
		
		if match:
			count += 1
			currentChrom.numMatchingRecords += 1
			if VERBOSE:
				print('Match = %s' %match)
				print('Match # %s' %count)
				print( expEntry )
				print( sourceEntry )
				print( expEntry.BASE_CALLS )
				print( sourceEntry.BASE_CALLS )
			assert sourceEntry.REF == expEntry.REF
			currentChrom.countMatches( expEntry.BASE_CALLS, sourceEntry.BASE_CALLS, sourceEntry.REF )
			if VERBOSE:
				print( 'Heterozygous:' )
				print( currentChrom.heterozygous_position_matches )
				print( currentChrom.heterozygousMatches )
				print( 'no call: %s' %currentChrom.heterozygousNoCall )
				print( 'Homozygous:' )
				print( currentChrom.homozygous_position_matches )
				print( currentChrom.homozygousMatches )
				print( 'no call: %s' %currentChrom.homozygousNoCall )

				print( currentChrom.sourceBaseCalled )
				print( currentChrom.sourceBaseDiffFromRef )
				print( wtf )

			if DEBUG:
				if expEntry.CHROM == DUMP_CHROM:
					debugDump.write( str(expEntry) )
					debugDump.write( str(sourceEntry) )

	logprint( '\n' + addTimestamp('Reached end of expFile.') )
	if not sourceFileExhausted:
		try:
			sourceLine = sourceFile.readline()
		except:
			sourceLine = None # should catch the case where sourceFile had just read in its last line when expFile exhausted
		
		while sourceLine:
				sourceSkippedRecords += 1
				sourceLine = sourceFile.readline()

wrapup(finalLinecounts=True)
logprint( '\n' + addTimestamp('RUN COMPLETE') )