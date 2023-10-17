#/bin/bash
##dic need -r, while file no -r 
for k in `ls`; do cp slurm.sh ${k}/; done
for k in `ls`; do cp hyper.py ${k}/; done