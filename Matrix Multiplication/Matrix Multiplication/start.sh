touch results.txt
rm results.txt
echo "-----NP = 1------"
echo "-----NP = 4------"
echo "-----NP = 8------"
echo "-----NP = 12------"
for ((n=0;n<5;n++))
do
  mpirun -np 12 ./matrixmult < input.txt >> results.txt
done
echo "-----NP = 14------"
for ((n=0;n<5;n++))
do
  mpirun -np 14 ./matrixmult < input.txt >> results.txt
done
echo "-----NP = 20------"
for ((n=0;n<5;n++))
do
  mpirun -np 20 ./matrixmult < input.txt >> results.txt
done
