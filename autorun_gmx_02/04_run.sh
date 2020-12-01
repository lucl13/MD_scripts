gro=$1    
tpr=$2     
path=3.PETase-MHETase/1.no-substrate/


papp_cloud scp  dnlu@cg12:/home/dnlu/project/pengxue/"$path"/*"$gro"*   ./
wait
echo "\"step1 "$gro" gro file is download\""


nohup bash  03_auto_run-"$tpr".sh >03_auto_run-"$tpr".dat  2>&1 &
wait
echo "\"step2 "$tpr" tpr file generation is done\""

papp_cloud scp ./"*"$tpr".tpr"  dnlu@cg12:/home/dnlu/project/pengxue/"$path"
wait
echo "\"step3 "$tpr" tpr files uploaded to the server\""

exit
