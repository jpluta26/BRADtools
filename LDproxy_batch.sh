
#!/bin/bash

# get the number of input parameters
NPARAM=$#

# if user enters less than the required number of arguments, print the usage
if [ $NPARAM -lt 1 ]
then
	echo ""
	echo "USAGE :: "
	echo "./LDproxy_batch.sh SNP_LIST [TOKEN] [POP] [WINDOW] [R2_THRESHOLD] [GENOME_BUILD]"
	echo "    SNP_LIST      : A file containing a list of SNP IDs (one per line)."
	echo "    TOKEN         : Your LDlink token. Get from https://ldlink.nih.gov/?tab=apiaccess"
	echo "    POP           : The population code (default: ALL)."
	echo "    WINDOW        : The window size in base pairs (default: 50000 N.B. NIH default: 500000)."
	echo "    R2_THRESHOLD  : The minimum r2 threshold (default: 0.8)."
	echo "    GENOME_BUILD  : Genome build, one of: grch37, grch38, grch38_high_coverage (default: grch38)."
	echo ""
	exit 1
fi

# make sure you're not on the head node
HEADNODE="consign.hpc.local"

if [ "$HOSTNAME" == "$HEADNODE" ]
then
  echo "Don't run scripts on the head node"
  exit 1
fi

# assign inputs to variable names; use defaults if not provided
SNP_LIST=$1
TOKEN=$2
POP=${3:-"ALL"}
WINDOW=${4:-"50000"}
R2_THRESHOLD=${5:-"0.8"}
GENOME_BUILD=${6:-"grch38"}

# STOP if the SNP list file is not found
if [ ! -e "$SNP_LIST" ]
then
	echo "$SNP_LIST :: file not found."
	exit 1
fi


# STOP if the token is not found
if [ ! -e "$TOKEN" ]
then
	echo "Get your token from https://ldlink.nih.gov/?tab=apiaccess"
	exit 1
fi

# Validate the genome build option
if [[ "$GENOME_BUILD" != "grch37" && "$GENOME_BUILD" != "grch38" && "$GENOME_BUILD" != "grch38_high_coverage" ]]
then
    echo "Invalid genome_build specified. Options are: grch37, grch38, grch38_high_coverage."
    exit 1
fi

# Process each SNP in the input SNP list file
while read -r SNP; do
    # Skip empty lines
    if [ -z "$SNP" ]; then
        continue
    fi

    echo "Processing SNP: $SNP"

    # Construct the LDproxy API URL for the current SNP
    URL="https://ldlink.nih.gov/LDlinkRest/ldproxy?var=${SNP}&pop=${POP}&r2_d=r2&window=${WINDOW}&genome_build=${GENOME_BUILD}&token=${TOKEN}"

    # Call the API, filter the output directly, and save the filtered output.
    # The awk command assumes the first line is a header and that column 7 contains the LD (r2) value.
    curl -s -k -X GET "$URL" | awk -v r2="$R2_THRESHOLD" 'NR==1 || $7 > r2' > "${SNP}_LDgroup.txt"
    
    if [ $? -ne 0 ]; then
        echo "Error processing SNP $SNP. Skipping."
    else
        echo "Filtered results saved to ${SNP}_filtered.txt"
    fi
done < "$SNP_LIST"
