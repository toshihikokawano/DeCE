#!/bin/sh

source parameter.sh
outorg=output.org

for testcase in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
  ./decetest.sh $testcase
  diff ./$outbase/ENDFout$testcase.out ./$outorg/ENDFout$testcase.out
done
