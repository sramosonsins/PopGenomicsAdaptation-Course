## declare an array variable
declare -a slim_filenames=("slim_template-SNMR0.slim" "slim_template-SNMRh.slim" "slim_template-SSR0.slim" "slim_template-SSRh.slim" "slim_template-PSL.slim" "slim_template-PSM.slim" "slim_template-PSH.slim" "slim_template-SNM-EXP.slim" "slim_template-SNM-RED.slim" "slim_template-DELF.slim" "slim_template-DELG.slim" "slim_template-DELF-EXP.slim" "slim_template-PSL-DELF.slim" "slim_template-PSM-DELF.slim" "slim_template-PSH-DELF.slim")

## RUN SLIM
for n in "${slim_filenames[@]}"
do
   slim -t -m $n > ${n}.out &
done

##RUN ALL THE STEPS
#for n in "${slim_filenames[@]}"
#do
#  sh sh run_alpha_steps_after_slim_smulation.sh $n
#done
