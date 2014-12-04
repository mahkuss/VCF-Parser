vcf = VCF('/home/rob/data/RNAseq/Naoki/wgs/c9orfKOMouse.vcf')
SNPstats = vcf.VCFstats(mode='snp')
SNPstats400 = vcf.VCFstats(quality=400, mode='snp')
#SNPstats400 = vcf.VCFstats(quality=400, chrom='chr4')

def filt(x):
	return x < 2000

fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, subplot_kw={'xlabel':'Quality', 'ylabel':'Count'})

axes[0,0].hist(filter(filt, SNPstats400[3]['qual']['heterozygous']), bins=100)
axes[0,1].hist(filter(filt, SNPstats400[0]['qual']['heterozygous']), bins=100)
axes[1,0].hist(filter(filt, SNPstats400[1]['qual']['heterozygous']), bins=100)
axes[1,1].hist(filter(filt, SNPstats400[2]['qual']['heterozygous']), bins=100)

axes[0,0].set_title('WT', fontsize=18)
axes[0,1].set_title('Het', fontsize=18)
axes[1,0].set_title('Homo', fontsize=18)
axes[1,1].set_title('NeoHomo', fontsize=18)

fig.suptitle('Heterozygous SNPs', fontsize=22)
fig.show()


fig2, axes2 = plt.subplots(2, 2, sharex=True, sharey=True, subplot_kw={'xlabel':'Quality', 'ylabel':'Count'})

axes2[0,0].hist(filter(filt, SNPstats400[3]['qual']['homozygous_alt']), bins=100)
axes2[0,1].hist(filter(filt, SNPstats400[0]['qual']['homozygous_alt']), bins=100)
axes2[1,0].hist(filter(filt, SNPstats400[1]['qual']['homozygous_alt']), bins=100)
axes2[1,1].hist(filter(filt, SNPstats400[2]['qual']['homozygous_alt']), bins=100)

axes2[0,0].set_title('WT', fontsize=18)
axes2[0,1].set_title('Het', fontsize=18)
axes2[1,0].set_title('Homo', fontsize=18)
axes2[1,1].set_title('NeoHomo', fontsize=18)

fig2.suptitle('Homozygous SNPs', fontsize=22)

fig2.show()


from pandas import Series, DataFrame
import pandas as pd

chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
				'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']

het_data = { 'WT': [], 'Het': [], 'Homo': [], 'NeoHomo': [] }
homo_data = { 'WT': [], 'Het': [], 'Homo': [], 'NeoHomo': [] }
for chromosome in chromosomes:
	data = vcf.VCFstats(hfilter=['QD < 2', 'FS > 60', 'MQ < 40', 'MQRankSum < -12.5', 'ReadPosRankSum < -8.0'], mode='snp', chrom=chromosome)
	het_data['WT'].append( data[3]['heterozygous'] )
	het_data['Het'].append( data[0]['heterozygous'] )
	het_data['Homo'].append( data[1]['heterozygous'] )
	het_data['NeoHomo'].append( data[2]['heterozygous'] )

	homo_data['WT'].append( data[3]['homozygous_alt'] )
	homo_data['Het'].append( data[0]['homozygous_alt'] )
	homo_data['Homo'].append( data[1]['homozygous_alt'] )
	homo_data['NeoHomo'].append( data[2]['homozygous_alt'] )

het_calls = DataFrame(het_data, index=chromosomes)
homo_calls = DataFrame(homo_data, index=chromosomes)

hetplot = het_calls.plot(kind='bar')
homoplot = homo_calls.plot(kind='bar')