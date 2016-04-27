touch results_n=20.txt
rm results_n=20.txt
echo "-----NP ==20------"
for ((n=0;n<5;n++))
do
  mpirun -hostfile hostfiles/10pcs -np 20 ./matrixmult < input.txt >> results.txt
done
echo "Finished"
