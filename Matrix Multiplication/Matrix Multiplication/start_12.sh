touch results_n=12.txt
rm results_n=12.txt
echo "-----NP = 12------"
for ((n=0;n<5;n++))
do
  mpirun -hostfile hostfiles/6pcs -np 12 ./matrixmult < input.txt >> results.txt
done
echo "Finished"
