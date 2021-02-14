"""
Multi-genome oligo design for MPRA
"""

#################
# CONFIG FILE:
#################

BED = config["target_bed"]
GENOME1 = config["source_genome"]
GTF = config["gtf"]

GENOMES = config["all_genomes"]
NAMES = config["names"]
FILTERS = config["filter_parameters"]

N_CONTROLS = config["n_controls"]

UPS5 = config["universal_priming_seq_fiveprime"]
UPS3 = config["universal_priming_seq_threeprime"]
RE = config["re_sites"]

#################
# WORKFLOW SETUP:
#################

TARGETS = "final_sequences.fa"

rule all:
	input:
		TARGETS

#################
# PIPELINE:
#################

# Tile 200mer windows at 2x coverage

rule create_tiles:
	input:
		bed=BED,
		genome=GENOME1,
		gtf=GTF
	output:
		merge="1_tiles/out_merged.bed",
		tiles="1_tiles/out_tiles.bed",
		strand="1_tiles/out_tiles.strand.bed",
		fa="1_tiles/out_tiles.fa"
	shell:
		"""
		# merge input bed within 400 bp (the maximum that would be joined by adjacent 200mers)
		bedtools merge -d 400 -i {input.bed} > {output.merge}

		# expand merged intervals such that they are multiples of 100...
		# int(($2+$3)/2) is the midpoint of the interval
		# int(($3-$2)/100+1)*100 is the interval size rounded up to the nearest 100
		# divide this quantity by two and add 100 to get the margin from midpoint
		# keep only 200mers, and name intervals with the coordinate

		awk -v OFS='\t' '{{ print $1, int(($2+$3)/2) - (int(($3-$2)/100+1)*100)/2 - 100, int(($2+$3)/2) + (int(($3-$2)/100+1)*100)/2 + 100 }}' {output.merge} |\
		awk -v OFS='\t' '{{if ($2<0) $2=0; print}}' |\
		bedtools makewindows -b stdin -w 200 -s 100 |\
		awk -v OFS='\t' '{{if ($3-$2==200) print}}' |\
		awk -v OFS='\t' '{{$4=$1":"$2"-"$3; print}}'> {output.tiles}

		# assign to strand of nearest transcribed feature (if tied, then choose at random)
		bedtools closest -a {output.tiles} -b {input.gtf} |\
		awk -v OFS='\t' '{{print $1,$2,$3,$4,0,$11}}' |\
		awk -v OFS='\t' -v seed=$RANDOM 'BEGIN {{srand(seed)}}; {{$(NF+1)=int(rand()*100000000)}}1' |\
		sort -k1,1 -k2,2n -k7,7n | sort -k1,1 -k2,2n -u | cut -f1-6 > {output.strand}
		
		# extract FASTA
		bedtools getfasta -s -name -fi {input.genome} -bed {output.strand} -fo {output.fa}
		"""

# BLAT tiles to all genomes of interest and filter in two rounds

NAMED_GENOMES = dict(zip(NAMES, GENOMES))
def get_genome(wildcards):
        return NAMED_GENOMES.get(wildcards.name)

NAMED_FILTERS = dict(zip(NAMES, FILTERS))
def get_filter(wildcards):
        return NAMED_FILTERS.get(wildcards.name)

rule blat:
	input:
		"1_tiles/out_tiles.fa",
	params:
		genome=get_genome,
		filter=get_filter
	output:
		align="2_blat/out_align.{name}.blast8",
		log="2_blat/out_filt_log.{name}.txt",
		filt="2_blat/out_filt.{name}.blast8"
	shell:
		"""
		module load blat
		blat -noHead -out=blast8 {params.genome} {input} {output.align}
		python scripts/doubleFiltBlat.py {output.align} {params.filter} {output.log} > {output.filt}
		"""

# Find tiles common to all genomes

rule common:
	input:
		expand("2_blat/out_filt.{name}.blast8", name=NAMES)
	output:
		tmp1="3_common/out_common.tmp1",
		tmp2="3_common/out_common.tmp2",
		comm="3_common/out_common.txt"
	shell:
		"""
		# IDs from first genome
		awk '{{print $1}}' {input[1]} | sort -u | grep -v quer > {output.comm}		

		# iterate over all genomes and find intersection
		for i in {input}; do
			awk '{{print $1}}' $i | sort -u > {output.tmp1}
			comm -12 {output.comm} {output.tmp1} > {output.tmp2}
			cp {output.tmp2} {output.comm}
		done
		"""

# Intersect filtered alignments with common tiles, fix coordinates, and extract FASTA

NAME2NAME = dict(zip(NAMES, NAMES))
def get_name(wildcards):
        return NAME2NAME.get(wildcards.name)

rule intersect:
	input:
		filt="2_blat/out_filt.{name}.blast8",
		ids="3_common/out_common.txt"
	params:
		name=get_name,
		genome=get_genome
	output:
		int="3_common/out_int.{name}.blast8",
		bed="3_common/out_int.{name}.bed",
		fa="3_common/out_int.{name}.fa"
	shell:
		"""
		grep -f {input.ids} {input.filt} > {output.int}
		python scripts/blast8_to_oligoBED.py {output.int} 190 210 200 {params.name} | sort -k1,1 -k2,2n > {output.bed}
		bedtools getfasta -s -name -fi {params.genome} -bed {output.bed} -fo stdout | sed 's/::.*$//g' > {output.fa}
		"""

# Combine FASTAs from all genomes and dedup; deduplicated entries are combined with a pipe

rule combine:
	input:
		expand("3_common/out_int.{name}.fa", name=NAMES)
	output:
		"3_common/out_all.fa"
	shell:
		"""
		cat {input} |\
		awk '{{ if ($0 !~ />/) {{print toupper($0)}} else {{print $0}} }}' |\
		awk '{{if(NR%2){{printf("%s\t",$1)}}else{{print}}}}' |\
		sort -k2,2 |\
		awk '{{if($2==SEQ){{gsub(">","",$1);ID=ID"|"$1}} else{{if(SEQ!=""){{printf("%s\\n%s\\n", ID,SEQ);}}SEQ=$2;ID=$1;}} }}END{{printf("%s\\n%s\\n", ID,SEQ)}}' > {output}
		"""


# Sample random source tiles and scramble to create negative controls

rule controls:
	input:
		ids="3_common/out_common.txt",
		genome=GENOME1
	params:
		N_CONTROLS
	output:
		bed="4_controls/out_random.bed",
		fa="4_controls/out_random.fa",
		ctrl="4_controls/out_controls.fa"
	shell:
		"""
		shuf -n {params} {input.ids} | sed 's/:/\t/g' | sed 's/-/\t/g' > {output.bed}
		bedtools getfasta -bed {output.bed} -fi {input.genome} -fo {output.fa}
		cat {output.fa} | python scripts/scrambleFasta.py |\
		awk '{{ if ($0 !~ />/) {{print toupper($0)}} else {{print $0}} }}' > {output.ctrl}
		"""

# Combine, add universal priming sequence, and remove sequences containing RE cutsites

rule final:
	input:
		blat="3_common/out_all.fa",
		ctrl="4_controls/out_controls.fa",
	params:
		ups5=UPS5,
                ups3=UPS3,
		re=RE
	output:
		ups="5_final/out_final.txt",
		tmp="5_final/out_final.tmp",
		final="final_sequences.fa"
	shell:
		"""
		# add UPS		
		cat {input.blat} {input.ctrl} | awk -v ups5={params.ups5} -v ups3={params.ups3} '{{ if ($0 !~ "^>") {{ $0=ups5$0ups3 }} }} {{print $0}}' |\
		awk '{{if(NR%2){{printf("%s\t",$1)}}else{{print}} }}' > {output.ups}
		
		# remove sequences containing RE sites
		for re in {params.re}; do
			grep -v $re {output.ups} > {output.tmp}
			cp {output.tmp} {output.ups}
		done

		sed 's/\t/\\n/g' {output.ups} > {output.final}
		"""
