#!/bin/sh
source parameter.sh

#-------------------------------------------------------------------------------
#  Option setting
#

# chargediscrete
#    charged particle emissions to discrete levels given in MT >= 600
#    otherwise included in MF = 103 - 107

chargediscrete=true

# discretegammaspec
#    include discrete gamma-ray spectrum in MF6 rather than MF12

discretegammaspec=true

# withconversion
#    discrete gamma transition table constains gamma-ray probabilities
#    need to prepare different data file, see gbranch.dat and gbranch_withconversion.dat

withconversion=false

# resonanceangdist
#    generate Legendre coefficients from given resonance parameters
#    need MF2 and MF/MT=4/2 first, then data inside RRR will be replaced

resonanceangdist=false


#-------------------------------------------------------------------------------
#  Compile tools if not exist
#

if [ ! -f $toolbase/$decemf6 ]; then
    echo "compiling $decemf6"
    pushd $toolbase ; make $decemf6 ; popd
fi
if [ ! -f $toolbase/$decemf12 ]; then
    echo "compiling $decemf12"
    pushd $toolbase ; make $decemf12 ; popd
fi


#-------------------------------------------------------------------------------
#  Read CoH calculated data
#

echo "read data file and create a temporal file, $workfile"

if [ $chargediscrete = "true" ]; then
    if [ $discretegammaspec = "true" ]; then
        echo "charged particle discrete transitions included"
        echo "discrete gamma-ray in MF6"
        sed -e "s/DATADIR/${datadir}/g"   \
            -e "s/DISCRETEGAMMASPEC/\"${datadir}\/gbranch.dat\"/g" \
            ${inbase}/ENDFdataread_chargediscrete.dece > $deceinput
    else
        echo "charged particle discrete transitions included"
        echo "discrete gamma-ray in MF12"
        sed -e "s/DATADIR/${datadir}/g"   \
            -e "s/DISCRETEGAMMASPEC//g" \
            ${inbase}/ENDFdataread_chargediscrete.dece > $deceinput
    fi
else
    if [ $discretegammaspec = "true" ]; then
        echo "charged particle cross sections lumped"
        echo "discrete gamma-ray in MF6"
        sed -e "s/DATADIR/${datadir}/g"   \
            -e "s/DISCRETEGAMMASPEC/\"${datadir}\/gbranch.dat\"/g" \
            ${inbase}/ENDFdataread.dece > $deceinput
    else
        echo "charged particle cross sections lumped"
        echo "discrete gamma-ray in MF12"
        sed -e "s/DATADIR/${datadir}/g"   \
            -e "s/DISCRETEGAMMASPEC//g" \
            ${inbase}/ENDFdataread.dece > $deceinput
    fi
fi

$dece $deceoption -o $workfile $template < $deceinput


#-------------------------------------------------------------------------------
#  Make MF6 section
#

echo "produce MF6, and store in $work6"

cat /dev/null > $workmf6
if [ $chargediscrete = "true" ]; then
    for mt in 16  17  22  28  32  33  34  37  41  42  44  45  91 102 108 111 112 117 649 699 749 799 849;do
        echo "make MF6 $mt from ${datadir}/CoHSpectrum.dat"
        $toolbase/$decemf6 -t $mt -e ${datadir}/CoHSpectrum_chargediscrete.dat -f $workfile >> $workmf6
    done
else
    for mt in 16  17  22  28  32  33  34  37  41  42  44  45  91 102 103 104 105 106 107 108 111 112 117;do
        echo "make MF6 $mt from ${datadir}/CoHSpectrum.dat"
        $toolbase/$decemf6 -t $mt -e ${datadir}/CoHSpectrum.dat -f $workfile >> $workmf6
    done
fi


#-------------------------------------------------------------------------------
#  Make MF12 and MF14 sections if needed
#

if [ $discretegammaspec = "false" ]; then
    echo "produce MF12 and MF14, and store in $workmf12"
    if [ $withconversion = "true" ]; then
        $toolbase/$decemf12 ${datadir}/gbranch_withconversion.dat $workfile > $workmf12
    else
        $toolbase/$decemf12 ${datadir}/gbranch.dat $workfile > $workmf12
    fi
fi



#-------------------------------------------------------------------------------
#  Assemble all pieces
#

echo "assemble all results, $workfile, $workmf6, and $workmf12"

if [ $discretegammaspec = "true" ]; then
    $dece -v -o $outfile $workfile <<EOF
libread  6 0 "$workmf6"
EOF
else
    if [ $chargediscrete = "true" ]; then
        $dece -v -o $outfile $workfile <<EOF
libread  6 0 "$workmf6"
libread 12 0 "$workmf12"
EOF
    else
        $dece -v -o $outfile $workfile <<EOF
libread 6  0 "$workmf6"
multilibread 12 51 90 "$workmf12"
multilibread 14 51 90 "$workmf12"
EOF
    fi
fi



#-------------------------------------------------------------------------------
#  Replace scattering angular distributions by resonances
#

if [ $resonanceangdist = "true" ]; then
    cp $outfile $workfile
    $dece -v -o $outfile $workfile <<EOF
resonanceangdist 200
EOF
fi


#-------------------------------------------------------------------------------
#  Finally Tape ID and line numbers, and delete work files
#

cp $outfile $workfile
$dece -v -o $outfile $workfile <<EOF
#set LineNumber
tpid $tapeid
EOF


rm $workfile
rm $deceinput

