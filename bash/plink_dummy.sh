#!/bin/bash

module load plink/1.07

cmd="plink $@"
eval "$cmd"
