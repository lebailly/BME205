PCEM := /projects/compbio/experiments/models.97
PCEM-SCRIPTS := ${PCEM}/scripts
PCB := /projects/compbio/bin
PCB-SCRIPTS := /projects/compbio/bin/scripts

default: Python2_f10.ids $(foreach x,Python2_f10_1 Python2_f10_2,${x}.ids ${x}.seqs)

MASTER ?= ${PCEM}/select-sets/dunbrack/2013/cullpdb_pc90_res1.6_R0.25_d131004_chains5130.ids 

ifdef REDEFINE

Python2_f10.ids: ${MASTER}
	${PCEM-SCRIPTS}/keep-ids ${MASTER} \
		< ${PCEM}/indexes/t06.ids \
		> $@

endif

%_1.ids %_2.ids: %.ids
	${PCEM-SCRIPTS}/partition-ids  <$^ \
		-rand_seed 3154765 \
		-num_groups 2 \
		-prefix $*_
	sort -o $*_1.ids $*_1.ids
	sort -o $*_2.ids $*_2.ids
	

%.x-seqs %.seqs : %.ids
	date
	-rm -f tmp-$*.ids
	-grep -v "_" $*.ids > tmp-$*.ids
	-chmod g+w tmp-$*.ids
	${PCEM-SCRIPTS}/get-seqs -ids tmp-$*.ids \
		-seq stride-mixed.seq tmp-$*.seqs \
		-seq stride-mixed.str2 tmp-$*.str2s \
		-seq stride-mixed.near-backbone-11 tmp-$*.near-backbone-11s
	${PCEM-SCRIPTS}/mask-lowercase < tmp-$*.seqs > tmp-$*.x-seqs
	mv -f tmp-$*.seqs $*.seqs
	mv -f tmp-$*.str2s $*.str2s
	mv -f tmp-$*.near-backbone-11s $*.near-backbone-11s
	mv -f tmp-$*.x-seqs $*.x-seqs
	rm -f tmp-$*.ids
	date

SCOP-IDS = /projects/compbio/data/scop/scop.ids

%-in-scop.ids: %.ids ${SCOP-IDS}
	${PCEM-SCRIPTS}/keep-ids $< < ${SCOP-IDS} \
	| uniq \
	> $@


# Set CODING-ARGS to --alphabet=ACDEFGHIKLMNPQRSTVWYZ to get coding costs
# with alphabet restricted to amino acids
CODING-ARGS ?=

COUNT-KMERS ?= ~karplus/bme205/f13/hw4-Markov/karplus/count-kmers
CODING-COST ?= ~karplus/bme205/f13/hw4-Markov/karplus/coding-cost ${CODING-ARGS}


define COUNT_RULE
%.count${order}:	%.seqs
	${COUNT-KMERS} -o${order} < $$< > $$@

endef

$(foreach order,0 1 2 3 4, $(eval $(call COUNT_RULE)))

count_tables:  Python2_f10_1.count0 Python2_f10_1.count1 \
	Python2_f10_2.count0 Python2_f10_2.count1 

CODING_FILES ?= $(foreach x,1 2,$(foreach y,1 2,$(foreach order,0 1 2 3 4,Python2_f10_$x.code_$y.${order})))

define CODING_RULE
%.code_${train}.${order}:  Python2_f10_${train}.count${order} %.seqs
	${CODING-COST} $$*.seqs < $$< > $$@

endef

$(foreach train,1 2,$(foreach order,0 1 2 3 4, $(eval $(call CODING_RULE))))

TINY_CODING_FILES ?= $(foreach order,0 1 2 3 4,tiny.code${order})

tiny.code%: tiny.count% tiny.seqs
	${CODING-COST} tiny.seqs < $< > $@

detected-alphabet-codings: count_tables
	mkdir -p detected-alphabet
	-rm ${CODING_FILES} ${TINY_CODING_FILES}
	${MAKE} CODING-ARGS= ${CODING_FILES} ${TINY_CODING_FILES}
	mv ${CODING_FILES}  ${TINY_CODING_FILES} detected-alphabet
AA-alphabet-codings: count_tables
	mkdir -p AA-alphabet
	-rm ${CODING_FILES} ${TINY_CODING_FILES}
	${MAKE} CODING-ARGS=--alphabet=ACDEFGHIKLMNPQRSTVWY ${CODING_FILES}  ${TINY_CODING_FILES}
	mv ${CODING_FILES}  ${TINY_CODING_FILES} AA-alphabet


