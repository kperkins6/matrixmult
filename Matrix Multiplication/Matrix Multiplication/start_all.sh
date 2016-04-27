touch results.txt
rm results.txt
echo "-----NP = 20------"
for ((n=0;n<5;n++))
do
mpirun -hostfile hostfiles/10pcs -np 20 ./matrixmult < input.txt >> results.txt
done
echo "-----NP = 16------"
for ((n=0;n<5;n++))
do
mpirun -hostfile hostfiles/8pcs -np 16 ./matrixmult < input.txt >> results.txt
done
echo "-----NP = 12------"
for ((n=0;n<5;n++))
do
mpirun -hostfile hostfiles/6pcs -np 12 ./matrixmult < input.txt >> results.txt
done
echo "-----NP = 8------"
for ((n=0;n<5;n++))
do
mpirun -hostfile hostfiles/8pcs -ppn 1 ./matrixmult < input.txt >> results.txt
done
echo "-----NP = 4------"
for ((n=0;n<5;n++))
do
mpirun -hostfile hostfiles/4pcs -ppn 1 ./matrixmult < input.txt >> results.txt
done
echo "-----NP = 1------"
for ((n=0;n<5;n++))
do
  mpirun -np 1 ./matrixmult < input.txt >> results.txt
done
