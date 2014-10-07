"""
Python tool for parsing VCF files.
"""

class VCFentry:
	"""
	Class object to parse and represent a single entry in a VCF file.
	"""
	def __init__(self, entry):
		fields = entry.strip().split('\t')
		self.CHROM, self.POS, self.ID, self.REF, self.ALT, self.QUAL, self.FILTER, self.INFO =\
		fields[:8]

		self.POS = int(self.POS)
		self.QUAL = float(self.QUAL)

		if len(fields) > 8:
			self.FORMAT = fields[8]
			self.GENOTYPE = [ genotype for genotype in fields[9:] ]

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
		return_lines = random.sample( xrange( self.header_line_count + 1, num_lines ), num )
		print(return_lines)

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
