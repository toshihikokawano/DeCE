#!/bin/sh
source parameter.sh

chargediscrete=true

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
    echo "charged particle discrete transisions included"
    $dece $deceoption -o $workfile $template < $inbase/ENDFdataread_chargediscrete.dece
else
    echo "general case"
    $dece $deceoption -o $workfile $template < $inbase/ENDFdataread.dece
fi

exit

echo "produce MF6, and store in $work6"
cat /dev/null > $workmf6
if [ $chargediscrete = "true" ]; then
    for mt in 16  22  28  102 649 849;do
        echo "make MF6 $mt from data/CoHSpectrum.dat"
        $toolbase/$decemf6 -t $mt -e data/CoHSpectrum_chargediscrete.dat -f $workfile >> $workmf6
    done
else
    for mt in 16  22  28  102 103 107;do
        echo "make MF6 $mt from data/CoHSpectrum.dat"
        $toolbase/$decemf6 -t $mt -e data/CoHSpectrum.dat -f $workfile >> $workmf6
    done
fi


echo "produce MF12 and MF14, and store in $workmf12"
$toolbase/$decemf12 data/gbranch.dat $workfile > $workmf12


echo "assemble all results, $workfile, $workmf6, and $workmf12"
$dece -o $outfile $workfile <<EOF
libread 6 16 "$workmf6"
libread 6 22 "$workmf6"
libread 6 28 "$workmf6"
libread 6 102 "$workmf6"
multilibread 12 51 90 "$workmf12"
multilibread 14 51 90 "$workmf12"
EOF

if [ $chargediscrete = "true" ]; then
    $dece -o $outfile $workfile <<EOF
libread 6 649 "$workmf6"
libread 6 849 "$workmf6"
EOF
else
    $dece -o $outfile $workfile <<EOF
libread 6 103 "$workmf6"
libread 6 107 "$workmf6"
EOF
fi
