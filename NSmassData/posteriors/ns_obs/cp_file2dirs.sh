#/bin/bash
##dic need -r, while file no -r 
for k in `ls`; do cp -r mian_dic ${k}/; done