#!/bin/sh
source parameter.sh

###############################
# chargediscrete
#    charged particle emissions to discrete levels given in MT >= 600

#chargediscrete=false
chargediscrete=true

###############################
# withconversion
#    discrete gamma transition table in MF12 includes gamma-ray probabilities

#withconversion=false
withconversion=true


if [ ! -f $toolbase/$decemf6 ]; then
    echo "compiling $decemf6"
    pushd $toolbase ; make $decemf6 ; popd
fi
if [ ! -f $toolbase/$decemf12 ]; then
    echo "compiling $decemf12"
    pushd $toolbase ; make $decemf12 ; popd
fi


echo "read data file and create a temporal file, $workfile"
if [ $chargediscrete = "true" ]; then
    echo "charged particle discrete transitions included"
    $dece $deceoption -o $workfile $template < $inbase/ENDFdataread_chargediscrete.dece
else
    echo "general case"
    $dece $deceoption -o $workfile $template < $inbase/ENDFdataread.dece
fi


echo "produce MF6, and store in $work6"
cat /dev/null > $workmf6
if [ $chargediscrete = "true" ]; then
    for mt in 16  22  28  91 102 649 699 749 799 849;do
        echo "make MF6 $mt from data/CoHSpectrum.dat"
        $toolbase/$decemf6 -t $mt -e data/CoHSpectrum_chargediscrete.dat -f $workfile >> $workmf6
    done
else
    for mt in 16  22  28  91 102 103 107;do
        echo "make MF6 $mt from data/CoHSpectrum.dat"
        $toolbase/$decemf6 -t $mt -e data/CoHSpectrum.dat -f $workfile >> $workmf6
    done
fi


echo "produce MF12 and MF14, and store in $workmf12"
if [ $withconversion = "true" ]; then
    $toolbase/$decemf12 data/gbranch_withconversion.dat $workfile > $workmf12
else
    $toolbase/$decemf12 data/gbranch.dat $workfile > $workmf12
fi


echo "assemble all results, $workfile, $workmf6, and $workmf12"
if [ $chargediscrete = "true" ]; then
    $dece -v -o $outfile $workfile <<EOF
libread 6  16 "$workmf6"
libread 6  22 "$workmf6"
libread 6  28 "$workmf6"
libread 6  91 "$workmf6"
libread 6 102 "$workmf6"
libread 6 649 "$workmf6"
libread 6 699 "$workmf6"
libread 6 749 "$workmf6"
libread 6 799 "$workmf6"
libread 6 849 "$workmf6"
multilibread 12  51  90 "$workmf12"
multilibread 12 600 640 "$workmf12"
multilibread 12 650 690 "$workmf12"
multilibread 12 700 740 "$workmf12"
multilibread 12 750 790 "$workmf12"
multilibread 12 800 840 "$workmf12"
multilibread 14  51  90 "$workmf12"
multilibread 14 600 640 "$workmf12"
multilibread 14 650 690 "$workmf12"
multilibread 14 700 740 "$workmf12"
multilibread 14 750 790 "$workmf12"
multilibread 14 800 840 "$workmf12"
set LineNumber
tpid "DeCE Example include MT600,800"
EOF
else
    $dece -v -o $outfile $workfile <<EOF
libread 6  16 "$workmf6"
libread 6  22 "$workmf6"
libread 6  28 "$workmf6"
libread 6  91 "$workmf6"
libread 6 102 "$workmf6"
libread 6 103 "$workmf6"
libread 6 107 "$workmf6"
multilibread 12 51 90 "$workmf12"
multilibread 14 51 90 "$workmf12"
set LineNumber
tpid "DeCE Example"
EOF
fi
