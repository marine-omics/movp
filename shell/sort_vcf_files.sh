for f in *.vcf;do f1=${f%.vcf}; f2=${f1#freebayes.}; echo $f2 | awk -F "[:-]" '{print $1,$2,$3}';done | sort -k2,2 -k4,4n | awk '{printf("freebayes.%s:%s-%s.vcf\n",$1,$2,$3)}' 
