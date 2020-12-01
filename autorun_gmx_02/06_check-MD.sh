# !/bin/bash

for (( i = 1; i < 51; i++ )); do
	gmx check -f "$i"_PETase-MHETase-md.xtc
done

