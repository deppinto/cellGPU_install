# !/bin/bash   

cd /home/diogo/CellGPU/CellGPU_Substrate2/
make
cd /home/diogo/CellGPU/

SCRIPT=5

mkdir -p scripts"$SCRIPT"

START=1
END=7
SEED_VAR=1

for ((ii=$START; ii<=$END; ii++))
do

JOBS=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $17}'`
#echo "$JOBS" 
for ((CURRJOB=$START; CURRJOB<=JOBS; CURRJOB++))
do

#echo "$SEED_VAR"
a=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $1}'`
b=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $2}'`
c=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $3}'`
d=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $4}'`
e=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $5}'`
f=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $6}'`
g=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $7}'`
h=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $8}'`
i=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $9}'`
j=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $10}'`
k=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $11}'`
l=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $12}'`
m=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $13}'`
n=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $14}'`
o=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $15}'`
p=`cat dados.txt | sed -n ${ii}','${ii}'p' | awk '{print $16}'`

cd scripts"$SCRIPT"/
mkdir -p Job_"$SEED_VAR"/
cd ..

cp -r /home/diogo/CellGPU/CellGPU_Substrate2/voronoi.out /home/diogo/CellGPU/scripts"$SCRIPT"/Job_"$SEED_VAR"

cp baseRun.sh ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{number}/$a/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{tSteps}/$b/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{GPU}/$c/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{initSteps}/$d/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{deltaT}/$e/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{perimeter}/$f/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{area}/$g/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{velocity}/$h/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{SubInt}/$i/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{Tau}/$j/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{IncVal}/$k/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{dx}/$l/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{sigma}/$m/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{boxsize}/$n/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{CellRadius}/$o/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{CellClusterSize}/$p/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh
sed -i "s/{Job}/$SEED_VAR/g" ./scripts"$SCRIPT"/Job_"$SEED_VAR"/run.sh

#sleep 1s

if [ $(find ./scripts"$SCRIPT"/Job_"$SEED_VAR"/ -name "*nc") ]; then 
  echo "File is found in Job_"$SEED_VAR""
else
  echo "File not found in Job_"$SEED_VAR""
  cd scripts"$SCRIPT"
  cd Job_"$SEED_VAR"
  qsub -V -l qs=true run.sh
  cd /home/diogo/CellGPU/
fi

#echo ./scripts/Job_"$SEED_VAR"/Submit__St${s}_T${z}_Iv${r}_P${p}_A${a}_V${v}.submit

((SEED_VAR++))

done
done
