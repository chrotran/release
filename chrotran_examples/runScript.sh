#!/bin/bash

date
pwd
#

for i in "$@"
do
case $i in
    -t=*|--testdir=*)
    TESTDIR="${i#*=}"
    ;;
    -f=*|--filename=*)
    FILENAME="${i#*=}"
    ;;
    -np=*|--nproc=*)
    NPROC="${i#*=}"
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
done

# Absolute path for chrotran_example directory
MYDIR=$PWD

# Absolute path for chrotran executable
PFLOTRAN_EXE=$MYDIR/../src/pflotran/chrotran

INPUT_FILE=${FILENAME}.in
OUTPUT_FILE=${FILENAME}.log

cd ${MYDIR}/${TESTDIR}

if [ -f "${OUTPUT_FILE}" ]; then
  echo "log file exists!"
  echo "deleting that file"
  rm -rf ${OUTPUT_FILE}
fi

if [ "${NPROC}" -gt 1 ]
then
echo "mpirun -np ${NPROC} ${PFLOTRAN_EXE} -pflotranin ${INPUT_FILE} >& ${OUTPUT_FILE}"
mpirun -np ${NPROC} ${PFLOTRAN_EXE} -pflotranin ${INPUT_FILE}  >& ${OUTPUT_FILE}
else
echo "${PFLOTRAN_EXE} -pflotranin ${INPUT_FILE} >& ${OUTPUT_FILE}"
${PFLOTRAN_EXE} -pflotranin ${INPUT_FILE} >& ${OUTPUT_FILE}
fi
