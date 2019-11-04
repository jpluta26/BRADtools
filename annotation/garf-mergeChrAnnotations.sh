#!/bin/bash
# pluta 10/14/19
# step 3 of garfield annotation process: join all the separate annotations into a single file
# this uses the output of garfieldCustomAnnotation.R


NPARAM=$#

if [ $NPARAM -lt 1 ]
then
        echo "USAGE :: "
        echo "mergeChrAnnotations.sh CHR"
        echo "CHR is 1..22 or X"
fi

HEADNODE="consign.hpc.local"

if [ "$HOSTNAME" == "$HEADNODE" ]
then
        echo "dont run scripts on the head node"
        exit 1
fi

CHR=$1

# multivariate equivalent of join
# taken from https://stackoverflow.com/questions/10726471/join-multiple-files
multijoin() {
        join_rec() {
            if [ $# -eq 1 ]; then
                join - "$1"
            else
                f=$1; shift
                join - "$f" | join_rec "$@"
            fi
        }

        if [ $# -le 2 ]; then
            join "$@"
        else
            f1=$1; f2=$2; shift 2
            join "$f1" "$f2" | join_rec "$@"
        fi
}



BEDSUFFIX="bed"$CHR
BEDS=(`ls *.$BEDSUFFIX`)
multijoin ${BEDS[@]} > chr${CHR}

# tmp1 holds the snp id, which is included in every file
# tmp2 holds the annotation information
more chr${CHR} | awk '{print $1}' > tmp1
more chr${CHR} | awk 'BEGIN{OFS="";} {$1=""; print $0}' > tmp2

paste -d " " tmp1 tmp2 > chr${CHR}

rm tmp1
rm tmp2
