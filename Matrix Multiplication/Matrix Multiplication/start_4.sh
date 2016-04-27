touch results_n=4.txt
rm results_n=4.txt
echo "-----NP = 4------"
for ((n=0;n<5;n++))
do
  mpirun -hostfile hostfiles/4pcs -ppn 1 ./matrixmult < input.txt >> results.txt
done
echo "Finished"
