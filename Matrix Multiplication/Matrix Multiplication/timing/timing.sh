#!/bin/bash


if [ $# == 0 ]; then
  echo "use flags -ijk -ikj or -kij to run tests on those forms"
  echo "or use the flag -all to run all forms (warning: will be long)"
else
  if [ $1 == "-ikj" ] || [ $1 == "-all" ]; then 
    rm results/ikj20.out
    for i in `seq 1 5`;
    do
      touch ./results/ikj20.out
      cat hostfiles/inputikj | mpirun -hostfile hostfiles/10pcs ../matrixmult >> ./results/ikj20.out
    done

    rm results/ikj16.out
    for i in `seq 1 5`;
    do
      touch ./results/ikj16.out
      cat hostfiles/inputikj | mpirun -hostfile hostfiles/8pcs ../matrixmult >> ./results/ikj16.out
    done

    rm results/ikj12.out
    for i in `seq 1 5`;
    do
      touch ./results/ikj12.out
      cat hostfiles/inputikj | mpirun -hostfile hostfiles/6pcs ../matrixmult >> ./results/ikj12.out
    done

    rm results/ikj8.out
    for i in `seq 1 5`;
    do
      touch ./results/ikj8.out
      cat hostfiles/inputikj | mpirun -hostfile hostfiles/8pcs -ppn 1 ../matrixmult >> ./results/ikj8.out
    done

    rm results/ikj4.out
    for i in `seq 1 5`;
    do
      touch ./results/ikj4.out
      cat hostfiles/inputikj | mpirun -hostfile hostfiles/4pcs -ppn 1 ../matrixmult >> ./results/ikj4.out
    done

    rm results/ikj1.out
    for i in `seq 1 5`;
    do
      touch ./results/ikj1.out
      cat hostfiles/inputikj | mpirun -n 1 ../matrixmult >> ./results/ikj1.out
    done
  fi


  ##################################################################################################
  ##
  ## 
  ##
  ##################################################################################################

  if [ $1 == "-kij" ] || [ $1 == "-all" ]; then 
    rm results/kij20.out
    for i in `seq 1 5`;
    do
      touch ./results/kij20.out
      cat hostfiles/inputkij | mpirun -hostfile hostfiles/10pcs ../matrixmult >> ./results/kij20.out
    done

    rm results/kij16.out
    for i in `seq 1 5`;
    do
      touch ./results/kij16.out
      cat hostfiles/inputkij | mpirun -hostfile hostfiles/8pcs ../matrixmult >> ./results/kij16.out
    done

    rm results/kij12.out
    for i in `seq 1 5`;
    do
      touch ./results/kij12.out
      cat hostfiles/inputkij | mpirun -hostfile hostfiles/6pcs ../matrixmult >> ./results/kij12.out
    done

    rm results/kij8.out
    for i in `seq 1 5`;
    do
      touch ./results/kij8.out
      cat hostfiles/inputkij | mpirun -hostfile hostfiles/8pcs -ppn 1 ../matrixmult >> ./results/kij8.out
    done

    rm results/kij4.out
    for i in `seq 1 5`;
    do
      touch ./results/kij4.out
      cat hostfiles/inputkij | mpirun -hostfile hostfiles/4pcs -ppn 1 ../matrixmult >> ./results/kij4.out
    done

    rm results/kij1.out
    for i in `seq 1 5`;
    do
      touch ./results/kij1.out
      cat hostfiles/inputkij | mpirun -n 1 ../matrixmult >> ./results/kij1.out
    done
  fi


  ##################################################################################################
  ##
  ## 
  ##
  ##################################################################################################

  if [ $1 == "-ijk" ] || [ $1 == "-all" ]; then 
    rm results/ijk20.out
    for i in `seq 1 5`;
    do
      touch ./results/ijk20.out
      cat hostfiles/inputijk | mpirun -hostfile hostfiles/10pcs ../matrixmult >> ./results/ijk20.out
    done
    
      rm results/ijk16.out
    for i in `seq 1 5`;
    do
      touch ./results/ijk16.out
      cat hostfiles/inputijk | mpirun -hostfile hostfiles/8pcs ../matrixmult >> ./results/ijk16.out
    done
    
      rm results/ijk12.out
    for i in `seq 1 5`;
    do
      touch ./results/ijk12.out
      cat hostfiles/inputijk | mpirun -hostfile hostfiles/6pcs ../matrixmult >> ./results/ijk12.out
    done
    
      rm results/ijk8.out
    for i in `seq 1 5`;
    do
      touch ./results/ijk8.out
      cat hostfiles/inputijk | mpirun -hostfile hostfiles/8pcs -ppn 1 ../matrixmult >> ./results/ijk8.out
    done
    
      rm results/ijk4.out
    for i in `seq 1 5`;
    do
      touch ./results/ijk4.out
      cat hostfiles/inputijk | mpirun -hostfile hostfiles/4pcs -ppn 1 ../matrixmult >> ./results/ijk4.out
    done
    
      rm results/ijk1.out
    for i in `seq 1 5`;
    do
      touch ./results/ijk1.out
      cat hostfiles/inputijk | mpirun -n 1 ../matrixmult >> ./results/ijk1.out
    done
  fi
fi

#mpirun -hostfile 2pcs -ppn 1 matrixmult >> script.out
#mpirun -hostfile 3pcs -ppn 1 matrixmult >> script.out
#mpirun -hostfile 5pcs -ppn 1 matrixmult >> script.out
#mpirun -hostfile 6pcs -ppn 1 matrixmult >> script.out
#mpirun -hostfile 7pcs -ppn 1 matrixmult >> script.out
#mpirun -hostfile 9pcs -ppn 1 matrixmult >> script.out
#mpirun -hostfile 10pcs -ppn 1 matrixmult >> script.out
#mpirun -hostfile 11pcs -ppn 1 matrixmult >> script.out
#mpirun -hostfile 7pcs matrixmult >> script.out
#mpirun -hostfile 9pcs matrixmult >> script.out
