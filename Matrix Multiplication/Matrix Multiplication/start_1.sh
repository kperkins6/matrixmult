touch results_n=1.txt
rm results_n=1.txt
echo "-----NP = 1------"
for ((n=0;n<5;n++))
do
  mpirun -np 1 ./matrixmult < input.txt >> results_n=1.txt
done
echo "Finished"
