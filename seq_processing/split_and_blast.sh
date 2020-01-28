awk -v size=200 -v pre=splitted -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' $1

for f in `ls splitted*` ; do

  blastn -remote -db nt -max_target_seqs 40 -outfmt "6 std qlen qcovs sgi sseq ssciname staxid" -out $f.blastout -qcov_hsp_perc 90 -perc_identity 80 -query $f

sleep 2m

done

cat *blastout >> $1.blasthits
