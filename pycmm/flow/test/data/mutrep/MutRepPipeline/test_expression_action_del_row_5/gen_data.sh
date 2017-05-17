#!/bin/bash

annotated_vcf="/proj/b2011117/private/databases/table_annovar/TYRCA/TYRCA_annotated.vcf.gz"
output_vcf="input.vcf"

tabix -h $annotated_vcf 1:1-1 > $output_vcf
tabix $annotated_vcf 1:1268000-1268000 >> $output_vcf
tabix $annotated_vcf 1:1334052-1334052 >> $output_vcf
tabix $annotated_vcf 1:1458900-1458900 >> $output_vcf
tabix $annotated_vcf 1:1857306-1857306 >> $output_vcf
tabix $annotated_vcf 1:3669356-3669356 >> $output_vcf
tabix $annotated_vcf 1:3683159-3683159 >> $output_vcf
tabix $annotated_vcf 1:47904668-47904668 >> $output_vcf
tabix $annotated_vcf 12:21176135-21176135 >> $output_vcf
tabix $annotated_vcf 17:16635980-16635980 >> $output_vcf
tabix $annotated_vcf 17:16675944-16675944 >> $output_vcf
tabix $annotated_vcf 17:17039561-17039561 >> $output_vcf

bgzip -f $output_vcf
tabix -p vcf $output_vcf.gz
