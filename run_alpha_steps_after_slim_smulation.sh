## declare the name of the file:
n=$1

##CUT TO MS FORMAT
sed '1,/Starting run/d'< ${n}.out > ${n}.out.ms

##RUN MSTATSPOP
./mstatspop -f ms -i ${n}.out.ms -o 0 -N 2 25 5 -G 1 -l 30000 -r 1 -n name_scaffold.txt -m mask_neutral.txt > ${n}.out.ms.neutral_statistics.txt
./mstatspop -f ms -i ${n}.out.ms -o 0 -N 2 25 5 -G 1 -l 30000 -r 1 -n name_scaffold.txt -m mask_functional.txt > ${n}.out.ms.functional_statistics.txt

##MAKE TABLE asymptotic-MK
grep 'Divergence\[0' ${n}.out.ms.functional_statistics.txt | tr '\t' ' ' | cut -d ' ' -f18 | sed '$d' > ${n}.out.ms.functional_divergence.txt #functional divergence
grep 'Divergence\[0' ${n}.out.ms.neutral_statistics.txt | tr '\t' ' ' | cut -d ' ' -f18 | sed '$d' > ${n}.out.ms.neutral_divergence.txt #neutral divergence
echo "freq\tsfs_funct\tsfs_neutral" > ${n}.SFS.table.txt #include a header.
grep 'fr\[0,' ${n}.out.ms.functional_statistics.txt | tr '\t' '\n' | cut -d ' ' -f2 > ${n}.out.ms.functional_statistics.txt.SFS.txt #create a column with functional SFS.
grep 'fr\[0,' ${n}.out.ms.neutral_statistics.txt | tr '\t' '\n' | cut -d ' ' -f2 > ${n}.out.ms.neutral_statistics.txt.SFS.txt #create a column with neutral SFS.
paste freq_col.txt ${n}.out.ms.functional_statistics.txt.SFS.txt ${n}.out.ms.neutral_statistics.txt.SFS.txt > ${n}.cols.txt #join the frequencies and the SFS for functional and neutral.
perl -ane 'print if $F[2]' ${n}.cols.txt > ${n}.cols_.txt #eliminate rows where syn is zero 
cat  ${n}.cols_.txt >> ${n}.SFS.table.txt #join header with SFS columns.

rm ${n}.out.ms.functional_statistics.txt.SFS.txt
rm ${n}.out.ms.neutral_statistics.txt.SFS.txt
rm ${n}.cols.txt
rm ${n}.cols_.txt

##RUN asymptotic-MK
div=$(cat ${n}.out.ms.functional_divergence.txt)
div0=$(cat ${n}.out.ms.neutral_divergence.txt) 
curl -F"d=${div}" -F"d0=${div0}" -F"xlow=0.05" -F"xhigh=0.95" -F"datafile=@${n}.SFS.table.txt" -o "${n}.asymptotic-MK_full.html" http://benhaller.com/cgi-bin/R/asymptoticMK_run.html

##MAKE MK from Theta estimates:
echo "ThetaWatt ThetaTaj ThetaFuLi ThetaFW" > ${n}.out.ms.functional_variability.txt #make a header.
echo "ThetaWatt ThetaTaj ThetaFuLi ThetaFW" > ${n}.out.ms.neutral_variability.txt #make a header.
grep 'Theta(' ${n}.out.ms.functional_statistics.txt | tr '\t' ' ' | cut -d ' ' -f4,6,8,10 >> ${n}.out.ms.functional_variability.txt #functional variability
grep 'Theta(' ${n}.out.ms.neutral_statistics.txt | tr '\t' ' ' | cut -d ' ' -f4,6,8,10 >> ${n}.out.ms.neutral_variability.txt #neutral variability
R --vanilla  --args $n < ./run_plots_Theta_alpha.R
