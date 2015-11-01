#!/bin/sh
DIR=/research/3rd_imp_data/users/bhoom/1kg/
for CHR in {1..22}; do
for i in {1..50}; do
    if [[ -e /research/3rd_imp_data/users/bhoom/1kg/imputed/aa/imputed_aa_chr${CHR}_pc${i}.gz ]]; then
	time sh prob2dose.sh /research/3rd_imp_data/users/bhoom/1kg/imputed/aa/imputed_aa_chr${CHR}_pc${i}.gz aa/aa${CHR}_pc${i} aa/aa.sample | tee aa/aa${CHR}_pc${i}.log
    else
	echo "../imputed/aa/imputed_aa_chr${CHR}_pc${i}.gz not found" >> aa/aa_files_not_found
    fi
done
done

for CHR in {1..22}; do
for i in {1..50}; do
    if [[ -e /research/3rd_imp_data/users/bhoom/1kg/imputed/ea/imputed_ea_chr${CHR}_pc${i}.gz ]]; then
	sh prob2dose.sh /research/3rd_imp_data/users/bhoom/1kg/imputed/ea/imputed_ea_chr${CHR}_pc${i}.gz ea/ea${CHR}_pc${i} ea/ea.sample  | tee ea/ea${CHR}_pc${i}.log
    else
	echo "../imputed/ea/imputed_ea_chr${CHR}_pc${i}.gz not found" >> ea/ea_files_not_found
    fi
done
done
