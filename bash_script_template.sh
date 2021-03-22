#!/bin/bash



# get the number of input parameters
NPARAM=$#

# if user enters less than the required number of arguments, print the usage
if [ $NPARAM -lt X ]
then
	echo ""
	echo "USAGE :: "
	echo "./template.sh ARG1 ARG2 ... "
	echo "ARG1 is BLANK"
	echo "ARG2 is BLANK"
	echo ""
	exit
fi

# make sure you're not on the head node
HEADNODE="consign.hpc.local"

if [ "$HOSTNAME" == "$HEADNODE" ]
then
  echo "dont run scripts on the head node!"
  exit 1
fi


# assign inputs to variable names
# $1 is the first variable entered, $2 is the second, etc
VAR1=$1
VAR2=$2


# STOP if any of the input files are not found
if [ ! -e $VAR1 ]
then
	echo " $VAR1 :: file not found."
	exit 1
fi

# setup a command
CMD="exe -opt1 $VAR1
	 -opt2 $VAR2
	 -out  $OUT"

# run the command
eval $CMD > logfile.log

# if the command failed to run, report an error
if [ $? -ne 0 ]
then
	echo "$CMD failed. Aborting."
	exit 1
fi
