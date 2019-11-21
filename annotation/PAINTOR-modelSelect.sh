#!/bin/sh

# pluta 10/4/19

# run this on macbook

# get the baseline model baye's factor
CMD="~/PAINTOR_V3.0/PAINTOR -input input.files -in . -out . -Zhead z -LDname ld -Gname enrich.Base -Lname BF.Base -enumerate 2"
eval $CMD

if [ $? -ne 0 ]
then
        echo "$CMD"
        echo "command failed"
        echo "PAINTOR requires- non-singular LD matrix; spaces instead of tabs"
        exit 1
fi

if [ -e univariate.stats ]
then
        rm univariate.stats
fi

if [ ! -e annotation.names ]
then
  echo "annotation.names is missing - run addFunctionalAnnotation.R first"
  exit 1
fi

if [ ! -s annotation.names ]
then
        echo "annotation.names is empty, something went wrong in the previous step"
        exit 1
fi

# run univariate model for each feature
# LRT between baseline model and feature to determine significance
for i in `more annotation.names | tr "," "\n"`
do
        ~/PAINTOR_V3.0/PAINTOR -input input.files -in . -out . -Zhead z -LDname ld -Gname enrich.${i} -Lname BF.${i} -enumerate 2 -annotations $i
        CMD="Rscript ../LRT.BF.R $i >> univariate.stats"
        eval $CMD

        if [ $? -ne 0 ]
        then
                echo "$CMD"
                echo "command failed!"
                exit 1
        fi
done


 more univariate.stats | tr -d '"' | awk -F "," '$2 < 0.05 {print $1}' > sig.features

 NSIGFEATURE=`wc -l sig.features | awk '{print $1}'`

 if [ $NSIGFEATURE -gt 0 ]
 then
        echo "n.sig.features = $NSIGFEATURE"
 elif [ $NSIGFEATURE -eq 0 ]
then
        echo "no significant features founds"
        exit
 else
        echo "n.sig.features = $NSIGFEATURE, something went wrong"
 fi

# next, run multivariate model. dont necessarily need to restrict to statistically significant features; even a combination
# of non-significant terms can be enough information to choose a causal variant
# the -enumerate value in the multivariate model should match that in the baseline model, otherwise BF wont be in the
# same scale
