names=('2p9h_R-IPT')


for name in $names; do
	#statements
	for (( i = 7; i < 16; i++ )); do
		#statements
		 { echo -e "1 | 13 \n"; echo -e "name 18 protein_lig \n"; echo -e "! 18 \n"; echo -e "name 19 envir \n";  echo -e "q"; }  | gmx make_ndx -f ${i}_${name}_em.gro -o ${i}_${name}_index.ndx
	done
done

