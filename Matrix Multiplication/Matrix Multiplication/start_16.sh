touch results_n=14.txt
rm results_n=14.txt
echo "-----NP = 16------"
for ((n=0;n<5;n++))
do
  mpirun -hostfile hostfiles/8pcs -np 16 ./matrixmult < input.txt >> results.txt
done
echo "Finished"
