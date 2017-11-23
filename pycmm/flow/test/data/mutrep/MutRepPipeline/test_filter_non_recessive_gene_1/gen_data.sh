#!/bin/bash

annotated_vcf="/proj/b2011117/private/databases/table_annovar/WES294/WES294_annotated.vcf.gz"
output_vcf="input.vcf"

tabix -h $annotated_vcf 1:1-1 > $output_vcf
tabix $annotated_vcf 2:179408713-179411011 >> $output_vcf
tabix $annotated_vcf 6:31324205-31324601 >> $output_vcf
tabix $annotated_vcf 6:31740760-31740763 >> $output_vcf
tabix $annotated_vcf 7:100550486-100550486 >> $output_vcf
tabix $annotated_vcf 19:17392630-17392631 >> $output_vcf
tabix $annotated_vcf 19:23159486-23159486 >> $output_vcf
tabix $annotated_vcf 19:32083223-32083250 >> $output_vcf

bgzip -f $output_vcf
tabix -p vcf $output_vcf.gz
