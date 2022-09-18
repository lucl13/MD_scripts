for (( i = 1; i < 51; i++ )); do
	name="$i"_PETase-MHETase

#nvt:
	#gmx grompp -f nvt.mdp -c "$i"_em.gro -r "$i"_em.gro -p "$i"_topol.top -o "$name"-nvt.tpr

#npt:
	gmx grompp -f npt.mdp -c "$name"-nvt.gro -r "$name"-nvt.gro -t "$name"-nvt.cpt  -p "$i"_topol.top -o "$name"-npt.tpr

	
#md:
	#gmx grompp -f md.mdp -c "$name"-npt.gro -r "$name"-npt.gro -t "$name"-npt.cpt -p "$i"_topol.top  -o "$name"-md.tpr 
	

	wait
	
done  

