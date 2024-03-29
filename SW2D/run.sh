#!/bin/bash
CPLUS_INCLUDE_PATH=MainCore
export CPLUS_INCLUDE_PATH
echo $CPLUS_INCLUDE_PATH
#export MAINCOREENV="MainCore"
#echo $MAINCOREENV
echo Compiling and Running Shallow Water 2D programm with $1.cpp as main program

if [ -n $2 ]
    then 
	export OMP_NUM_THREADS=$2
fi

g++ $1.cpp $CPLUS_INCLUDE_PATH/Raschet.cpp\
 $CPLUS_INCLUDE_PATH/Prepare_raschet.cpp\
  $CPLUS_INCLUDE_PATH/Problem_Definition.cpp\
   $CPLUS_INCLUDE_PATH/technical.cpp\
   $CPLUS_INCLUDE_PATH/Constants.cpp\
    $CPLUS_INCLUDE_PATH/Work_With_Techplot.cpp\
     $CPLUS_INCLUDE_PATH/BoundaryConditions.cpp\
      $CPLUS_INCLUDE_PATH/Transport_Problem.cpp\
       $CPLUS_INCLUDE_PATH/Numerical_scheme_Step_parallel.cpp\
        $CPLUS_INCLUDE_PATH/Time_control.cpp\
		 $CPLUS_INCLUDE_PATH/IO_system.cpp\
          $CPLUS_INCLUDE_PATH/Force.cpp\
    	  ephemeris/Ephemeris.cpp\
           ephemeris/Calendar.cpp\
       -o $1_compiled -std=c++11 -fopenmp -lm
./$1_compiled
