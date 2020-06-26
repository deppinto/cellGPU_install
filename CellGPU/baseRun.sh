#$ -S /bin/sh
#$ -cwd
#$ -j y
#$ -l virtual_free=1G -l h_vmem=1G

cp $SGE_O_WORKDIR/* $TMPDIR
cd $TMPDIR
(time ./voronoi.out -a "{number}" -b "{tSteps}" -c "{GPU}" -d "{initSteps}" -e "{deltaT}" -f "{perimeter}" -g "{area}" -h "{velocity}" -i "{SubInt}" -j "{Tau}" -k "{IncVal}" -l "{dx}" -m "{sigma}" -n "{boxsize}" -o "{CellRadius}" -p "{CellClusterSize}") &> time.txt
cp * $SGE_O_WORKDIR/
rm *
