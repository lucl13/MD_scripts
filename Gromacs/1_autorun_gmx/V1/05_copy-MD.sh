gro=$1    
path=3.PETase-MHETase/1.no-substrate/


papp_cloud scp  dnlu@cg12:/home/dnlu/project/pengxue/"$path"/*"$gro"*   ./
wait
echo "\"step1 "$gro" gro file is download\""


exit
