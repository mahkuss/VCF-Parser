"""
Python tool for parsing VCF files.

Author: Rob Moccia
"""

import numpy as np
import re

class VCFentry:
	"""
	Class object to parse and represent a single entry in a VCF file.
	"""
	def __init__(self, entry):
		fields = entry.strip().split('\t')
		self.CHROM, self.POS, self.ID, self.REF, self.ALT, self.QUAL, self.FILTER, self.rawINFO =\
		fields[:8]

		self.INFO = { key: float(value) if ',' not in value else (float(x) for x in value.split(',')) for key, value in [ stat.split('=') for stat in self.rawINFO.split(';') ] } 
		self.POS = int(self.POS)
		self.QUAL = float(self.QUAL)

		def genotypeClass(genotype):
			"""
			Helper function in VCFentry __init__ method to determine if variant is het, homo ref, or homo alt.
			"""
			if '/' in genotype: #unphased
				allele1, allele2 = genotype.split('/')
			elif '|' in genotype: # phased
				allele1, allele2 = genotype.split('|')
			else:
				return(None) # input to function is not a valid genotype

			if allele1 != allele2:
				return('heterozygous')
			elif allele1 == '0':
				return ('homozygous_ref')
			elif allele1 == '.':
				return('no_data')
			else:
				return('homozygous_alt')

		if len(fields) > 8:
			self.FORMAT = fields[8]
			self.GENOTYPE = [ genotype for genotype in fields[9:] ]
			self.GENOTYPE_CLASS = [ genotypeClass( genotype.split(':')[0] ) for genotype in fields[9:] ]

		if ( len(self.REF) > 1 ) or ( len(self.ALT) > 1 ):
			self.TYPE = 'indel'
		else:
			self.TYPE = 'snp'

	def __str__(self):
		return '\t'.join([ self.CHROM, self.POS, self.ID, self.REF, self.ALT, self.QUAL,\
			self.FILTER, self.INFO, self.FORMAT, '\t'.join( [ genotype for genotype in self.GENOTYPE ] ) ] ) + '\n'


class VCF:
	"""
	Parse a VCF file.
	"""
	def __init__(self, filename):
		self.filename = filename
		self.header = []
		self.header_line_count = 0
		with open(filename, 'r') as vcf:
			for line in vcf:
				if line[0] == '#':
					self.header.append(line)
					self.header_line_count += 1
				else:
					self.num_samples = len(line.strip().split('\t')) - 9
					self.format = VCFentry(line).FORMAT
					self.info = [ key for key in VCFentry(line).INFO.keys() ]
					break

		try:
			self.samples = self.header[-1].strip().split('\t')[9:]

		except:
			print('Warning: Could not find sample names in header. Header may be missing from supplied VCF.')
			self.samples = ['NA'] * self.num_samples


	def split(self, writeSNP=True, writeINDEL=True):
		"""
		Class method to split a .vcf file into separate files containing SNPs and INDELs, respectively
		"""

		if writeSNP:
			snp_filename = self.filename.split('.vcf')[0] + 'SNP.vcf'
			snp_out = open(snp_filename, 'w')
			for line in self.header:
				snp_out.write(line)

		if writeINDEL:
			indel_filename = self.filename.split('.vcf')[0] + 'INDEL.vcf'
			indel_out = open(indel_filename, 'w')
			for line in self.header:
				indel_out.write(line)

		with open(self.filename, 'r') as vcf:
			for ind, line in enumerate(vcf):
				if ind < self.header_line_count:
					continue
				elif ( VCFentry(line).TYPE == 'snp' ) and writeSNP:
					snp_out.write(line)
				elif writeINDEL:
					indel_out.write(line)

		if writeSNP:
			snp_out.close()

		if writeINDEL:
			indel_out.close()


	def filter(self, output, quality=0, variant_type=['snp', 'indel'], genotype=False, chrom=False, columns=(), coords=(0,1e10) ):
		"""
		Class method to generate a new vcf file containing only calls passing a specified filter.
		If genotype as specified (choices: homozygous, heterozygous) you must also specify a range of columns
		(as a tuple) to check for multi-sample vcf files. By default, all columns will be checked. Current behavior
		for genotype=heterozyous will output variants if ANY sample is heterozygous at that locus. For genotype=
		homozygous, the line is written to output only if ALL samples are homozgous at that site.
		"""

		number_checked = 0
		number_passed = 0

		with open(output, 'w') as outfile:
			for line in self.header:
				outfile.write(line)

			with open(self.filename, 'r') as vcf:
				for ind, line in enumerate(vcf):
					if ind < self.header_line_count:
						#print('%s: header line' %ind)
						continue

					else:
						current_line = VCFentry(line)
						number_checked += 1
						# check if chromosome filter set and, if so, is current line from that chromosome?
						if chrom and current_line.CHROM not in chrom:
							#print('%s: failed chrom filter' %ind)
							continue

						if ( current_line.POS < coords[0] ) or ( current_line.POS > coords[1] ):
							#print('%s: failed coordinates filter' %ind)
							continue

						# check if genotype filter set and, if so, is current line of the correct genotype?
						if genotype:
							#print('%s: checking genotype' %ind)
							if columns == ():
								columns = (0, len(self.samples))
							print('columns: %s,%s' %(columns[0], columns[1]))
							current_line_genotype = [ current_line.GENOTYPE[ column ].split(':')[0] for column in range(columns[0], columns[1]) ]
							for entry in current_line_genotype:
								print(entry)
								if entry.split('/')[0] != entry.split('/')[1]:
									current_line_genotype = 'heterozygous'
									print('heterozygous')
									break
								else:
									print('homozygous')
									current_line_genotype = 'homozygous'
							#print('genotype = ' + current_line_genotype)
							if current_line_genotype != genotype:
								#print('%s: failed genotype test' %ind)
								continue

						# check that current line passes the quality and variant_type filters
						if ( current_line.QUAL >= quality ) and ( current_line.TYPE in variant_type ):
							# if you get this far, current line passes all filters and should be written to output file
							#print('%s: passed all filters, writing to outfile' %ind)
							outfile.write(line)
							number_passed += 1

		print('Variants checked: %s' %number_checked)
		print('Variants passed: %s' %number_passed)


	def qualityDistribution(self):
		"""
		Class method to plot distribution of quality scores in vcf file.
		"""

		import matplotlib.pyplot as plt

		quality_scores = []

		with open(self.filename, 'r') as vcf:
			for ind, line in enumerate(vcf):
				if ind < self.header_line_count:
					continue

				else:
					quality_scores.append( VCFentry(line).QUAL )

		plt.hist(filter(lambda x: x < 2000, quality_scores), bins=100)
		plt.title('Quality Score Distribution')
		plt.xlabel('Quality Score')
		plt.ylabel('Frequency')
		plt.show()

	def randomSample(self, num, output=False):
		"""
		Class method to generate a specified number of randomly samples calls from a VCF file.
		This can be useful for manually checking quality of post-filter calls in IGV.
		"""

		import random

		# first find number of lines in VCF file
		num_lines = sum( 1 for line in open(self.filename) if line.rstrip() )
		
		# randomly select specified number of lines out of file, not including header
		return_lines = random.sample( xrange( self.header_line_count + 1, num_lines + 1 ), num )

		# iterate through input vcf file and print the randomly selected lines (or write to file if output path provided)
		if output:
			outfile = open(output, 'w')
		with open(self.filename, 'r') as vcf:
			for ind, line in enumerate(vcf):
				if ind < self.header_line_count:
					if output:
						outfile.write(line)
					else:
						continue

				elif ind in return_lines:
					if output:
						outfile.write(line)
					else:
						print(line)

		if output:
			outfile.close()


	def VCFstats(self, quality=0, chrom='ALL', hfilter=False):
		"""
		VCF class method to report statistics on contents of VCF file. Reports separate counts of heterozygous and
		homozygous variants above a specified quality threshold broken down by sample in a multi-sample VCF. Also 
		report mean and sd of quality scores for each category.
		"""

		# initialize dictionaries to hold data
		dict_categories = ['heterozygous', 'homozygous_ref', 'homozygous_alt', 'no_data', 'analyzed']

		data = []

		failed_filter = 0

		for i in range( self.num_samples ):
			data.append( dict() )
			for cat in dict_categories:
				data[i][cat] = 0
				data[i]['qual'] = {'heterozygous': [], 'homozygous_alt': []}


		# read VCF line by line
		with open(self.filename, 'r') as vcf:
			for ind, line in enumerate(vcf):
				if ind < self.header_line_count:
					continue

				else:
					variants = VCFentry( line )

					if chrom != 'ALL' and variants.CHROM != chrom:
						failed_filter += 1
						continue

					if hfilter:
						rules = { rule[0].strip():( rule[1].strip(), float(rule[2].strip() ) ) 
									for rule in [ re.split('(==|>|<)', arg ) for arg in hfilter ] }
						pass_rule = True
						try:
							for key, value in rules.items():
								check = variants.INFO[ key ]
								if ( value[0] == '==' and check != value[1] ) or\
								( value[0] == '>' and check > value[1] ) or\
								(value[0] == '<' and check < value[1] ):
									pass_rule = False
									# print(variants.INFO)
									# print('failed: %s %s %s' %(key, value[0], value[1]))
									break
						except KeyError:
							print('WARNING: Invalid filter rule ignored on line %s: %s %s %s' %(ind, key, value[0], value[1]))
							continue
						if not pass_rule:
							failed_filter += 1
							continue

					if variants.QUAL >= quality:
					# for each line loop from 0 to number of samples
						for sample in range(self.num_samples):

						# at each step in loop: 
						#	-increment an accumulator for homozygous ref, homozygous alt, and heterozygous calls
						# 	 indexed by sample and chromosome; maintain separate count for 'no_data'
						#	-add current QUAL to either het or homo list
							data[ sample ][ variants.GENOTYPE_CLASS[ sample ] ] += 1
							if variants.GENOTYPE_CLASS[ sample ] not in ['no_data', 'homozygous_ref', None]:
								data[ sample ][ 'analyzed' ] += 1
								data[ sample ][ 'qual' ][ variants.GENOTYPE_CLASS[ sample ] ].append( variants.QUAL )

					else:
						failed_filter += 1

		# calculate mean and sd of qual score categories
		for idx in range( len(data) ):
			print(idx)
			print('Analyzed: %s' %data[idx]['analyzed'])
			print('Homozygous Reference: %s' %data[idx]['homozygous_ref'])
			print('Homozygous Alt: %s' %data[idx]['homozygous_alt'])
			print('Homozygous Quality: %s, %s' %( np.mean(data[idx]['qual']['homozygous_alt']), np.std(data[idx]['qual']['homozygous_alt']) ) ) 
			print('Heterozygous: %s' %data[idx]['analyzed'])
			print('Heterozygous Quality: %s, %s' %( np.mean(data[idx]['qual']['heterozygous']), np.std(data[idx]['qual']['heterozygous']) ) )
		print('Failed Filter: %s' %failed_filter)

		return data
		# format and report output
