from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.stats as stats
import scipy.ndimage as ndimage

#filename = argv[1]
#OUTDIR = argv[2]

## 	FOR TESTING ONLY WILL CREATE FAKE VCF CLASS OBJECT TO WORK WITH ##
class TESTING:
	def __init__(self):
		self.num_samples = 4
		self.filename = '/home/rob/data/RNAseq/Naoki/wgs/VCFstats_filter_out.vcf'
		self.header_line_count = 94
self = TESTING()
###################################################################



# Initialize global variables
BIN_SIZE = 5e5
CHROM_DATA = '/home/rob/data/RNAseq/Naoki/neo_cassette/mm10chromInfo.txt'

# load in chromosome lengths to dictionary
ch_length = dict()

with open( CHROM_DATA, 'r' ) as f:
	for line in f:
		chrom = line.strip().split( '\t' )[:2]
		ch_length[ chrom[0] ] = int( chrom[1] )

# build bins and create dictionaries to hold number of homozygous and heterozygous snps in each bin
# data structure is a list of dicts, one for each sample in the vcf, with each dict holding two more
# dicts with keys "het" and "homo" for heterozygous and homozygous snps, respectively
# within the het and homo dicts will be keys corresponding to each chromosome which then hold a list
# of integers representing snp counts indexed by chromosomal coordinate start of the corresponding bin

snp = []
#totalGenes = dict()
#normalized = dict()

for sample in range( self.num_samples ):
	new_dict = {'heterozygous': dict(), 'homozygous_alt': dict() }
	for key, value in ch_length.items():
		new_dict[ 'heterozygous' ][ key ] = [0] * int( (value // BIN_SIZE) + 1 )
		new_dict[ 'homozygous_alt' ][ key ] = [0] * int( (value // BIN_SIZE) + 1 )
	#	totalGenes[key] = [0] * int( (value // BIN_SIZE) + 1 )
	#	normalized[key] = [0] * int( (value // BIN_SIZE) + 1 )

	snp.append( new_dict )

for key in snp[ 0 ][ 'heterozygous' ].keys():
	print( key + ': ' + str( len( snp[ 0 ][ 'heterozygous' ][ key ] ) ) )

print( len(snp) )

lines = 0

with open(self.filename, 'r') as f:
	for ind, line in enumerate(f):
		if ind < self.header_line_count:
			continue
		else:
			lines += 1
			data = VCFentry( line )
			if data.FILTER == 'PASS' or data.FILTER == '.':
				for sample in range(self.num_samples):
					try:
						# print(data.GENOTYPE_CLASS[ sample ])
						# print(data.CHROM)
						# print( int(data.POS // BIN_SIZE ))
						# print( snp[ sample ][ data.GENOTYPE_CLASS[ sample ] ][ data.CHROM ][ int( data.POS // BIN_SIZE ) ] )
						snp[ sample ][ data.GENOTYPE_CLASS[ sample ] ][ data.CHROM ][ int( data.POS // BIN_SIZE ) ] += 1 #sample will need to be argument specifying which column in multi-sample vcf
						# print( snp[ sample ][ data.GENOTYPE_CLASS[ sample ] ][ data.CHROM ][ int( data.POS // BIN_SIZE ) ] )
					except KeyError:
						print('genotype class not in index: %s', data.GENOTYPE_CLASS[ sample ])
						continue

print( 'lines: %s' %lines )
print( snp[ 0 ]['heterozygous']['chr5'])
print( snp[ 1 ]['heterozygous']['chr5'])
print( snp[ 2 ]['heterozygous']['chr5'])
print( snp[ 3 ]['heterozygous']['chr5'])




print
		try:
			hetSNP[chromosome][index] += int(q_val < SIG_CUTOFF)
		except IndexError:
			print(chromosome)
			print(index)
			print(locus)
			print(len(hetSNP[chromosome]))
			print(q_val)
			print(ensID)
		except:
			print("ERROR")


		totalGenes[chromosome][index] += 1
			#else:
			#	hetSNP[chromosome][index] = int(q_val < SIG_CUTOFF)
			#	totalGenes[chromosome][index] = 1

for key in totalGenes.keys():
	for bin in range(len(totalGenes[key])):
		if totalGenes[key][bin] != 0:
			normalized[key][bin] = hetSNP[key][bin] / totalGenes[key][bin]
		else:
			normalized[key][bin] = np.nan			

print(str(lines) + ' lines processed.')
print(str(status_ok + status_not_ok) + ' lines reported on desired conditions (WT v Homo)')
print("Status OK: " + str(status_ok))
print("Status NOT OK: " + str(status_not_ok))
print('\n')

clean = {}
target_bin = []
genes_per_window = []
for key in normalized.keys():
	clean[key] = []
	for i in range(len(normalized[key])):
		genes_per_window.append(totalGenes[key][i])

		if normalized[key][i] is not np.nan:
			print(str(i) + ' -- ' + 'bin: ' + str(int(i * BIN_SIZE)) + '\t' + 'percent DE: ' + str(round(normalized[key][i] * 100, 2)) + '%' +'\t' + str(totalGenes[key][i]) + ' genes.' )
			clean[key].append(normalized[key][i])

			if totalGenes[key][i] >= 8 and totalGenes[key][i] <= 12:
				target_bin.append(normalized[key][i])

		else:
			print(str(i) + ' -- ' + 'bin: ' + str(int(i * BIN_SIZE)) + '\t' + 'percent DE: NA')
	# try:

# for key in normalized.keys():
# 	try:
# 		fig = plt.figure()
# 		s = fig.add_subplot(111)
# 		# x,y = np.histogram(clean[key])

# 		# s.plot(x,y)
# 		s.hist(clean[key])
# 		fig.savefig(str(OUTDIR) + '/distribution_%s.png'%key)
# 		fig.close()
# 	except:
# 			print(clean[key])
# 			print('ERROR %s' % key)

fig = plt.figure()
s = fig.add_subplot(111)
weights = np.ones_like(target_bin)/len(target_bin)
# s.hist(target_bin, weights=weights, bins=np.arange(0,1,0.05))
# print(len(target_bin))
x,y = np.histogram(target_bin, weights=weights, bins=np.arange(0,1,0.05))
xs = ndimage.filters.gaussian_filter(x, 0.1, mode="constant")
xfit = stats.gaussian_kde(target_bin)
yf  = np.arange(0,1,0.01)
xf = np.array([xfit.evaluate(i) for i in yf])
xf = xf/xf.sum()
# print(len(xfit),len(y))
# s.plot(y[:-1],x,'r')
s.plot(yf,xf,'g')
s.plot([0.2,0.2],[0,0.03],'r')
s.plot([0.3,0.3],[0,0.03],'b')
plt.xlabel("Fraction Differentially Expressed Genes")
plt.ylabel("Frequency")
fig.savefig(str(OUTDIR) + '/target_bin_distribution.png')

fig = plt.figure()
s = fig.add_subplot(111)
weights = np.ones_like(genes_per_window) / len(genes_per_window)
s.hist(genes_per_window, weights=weights, bins=range(min(genes_per_window),max(genes_per_window)+5,5))
plt.xlabel('Number of Expressed Genes in %s bp Window'%BIN_SIZE)
plt.ylabel('Frequency')
fig.savefig(str(OUTDIR) + '/genes_per_window_distribution.png')

sort_order = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', 'X', 'Y', 'MT']

fig = plt.figure()
i = 0
for key in sorted(normalized, key=lambda x: sort_order.index(x), reverse=True):
	i += 1
	s = fig.add_subplot(len(normalized.keys()),1,i)
	s.imshow([normalized[key]]*4)
	s.set_ylabel(key)
	s.axes.get_xaxis().set_ticks([])
	s.axes.get_yaxis().set_ticks([])
fig.savefig(str(OUTDIR) + '/all_map.png')

	
	
# input()
# output_list = []
# sort_order = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', 'X', 'Y', 'MT']
# for item in counts.items():
# 	pct = round(item[1][0] / item[1][1] * 100, 2)
# 	output_list.append(str(item[0]) + ': ' + str(item[1][0]) + ' / ' + str(item[1][1]) + ' ... ' + str(pct) +'%')
# 	output_list.sort(key = lambda x: sort_order.index(x[:x.index(':')]))
# for line in output_list:
# 	print(line)

