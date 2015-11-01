#!/bin/sh
DIR=/research/3rd_imp_data/users/bhoom/1kg/
for CHR in {15}; do
for i in {12..50}; do
    if [[ -e /research/3rd_imp_data/users/bhoom/1kg/imputed/ea/imputed_ea_chr${CHR}_pc${i}.gz ]]; then
	time sh prob2dose.sh /research/3rd_imp_data/users/bhoom/1kg/imputed/ea/imputed_ea_chr${CHR}_pc${i}.gz ea/ea${CHR}_pc${i} ea/ea.sample | tee ea/ea${CHR}_pc${i}.log
    else
	echo "../imputed/ea/imputed_ea_chr${CHR}_pc${i}.gz not found" >> ea/ea_files_not_found
    fi
done
done

for CHR in {16..22}; do
for i in {1..50}; do
    if [[ -e /research/3rd_imp_data/users/bhoom/1kg/imputed/ea/imputed_ea_chr${CHR}_pc${i}.gz ]]; then
	sh prob2dose.sh /research/3rd_imp_data/users/bhoom/1kg/imputed/ea/imputed_ea_chr${CHR}_pc${i}.gz ea/ea${CHR}_pc${i} ea/ea.sample  | tee ea/ea${CHR}_pc${i}.log
    else
	echo "../imputed/ea/imputed_ea_chr${CHR}_pc${i}.gz not found" >> ea/ea_files_not_found
    fi
done
done
