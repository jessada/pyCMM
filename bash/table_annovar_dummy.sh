#!/bin/bash

module load annovar

cmd="table_annovar.pl $@"
eval "$cmd"
