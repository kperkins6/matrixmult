touch results_n=8.txt
rm results_n=8.txt
echo "-----NP = 8------"
for ((n=0;n<5;n++))
do
  mpirun -hostfile hostfiles/8pcs -ppn 1 ./matrixmult < input.txt >> results.txt
done
echo "Finished"
