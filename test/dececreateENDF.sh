#!/bin/sh

source parameter.sh

if [ ! -f $toolbase/$decemf6 ]; then
    echo "compiling $decemf6"
    pushd $toolbase ; make $decemf6 ; popd
fi
if [ ! -f $toolbase/$decemf12 ]; then
    echo "compiling $decemf12"
    pushd $toolbase ; make $decemf12 ; popd
fi

echo "read data file and create a temporal file $workfile"

$dece $deceoption -o $workfile $template <<EOF

echo "read total cross section from external library"
libread 3 1 "ENDFMF3MT1.evl"

echo "remove from 20 to 150 MeV data"
delpoint 3 1 2.000001e+7 1.5e+8

echo "read elastic scattering from file"
read 3 2 "data/CoHCrossSection.dat"

echo "read (n,gamm) reaction"
read 3 102 "data/CoHCrossSection.dat"

echo "read (n,2n)"
read 3 16 "data/CoHParticleProduction.dat"

echo "read (n,n alpha)"
read 3 22 "data/CoHParticleProduction.dat"

echo "read (n,n p)"
read 3 28 "data/CoHParticleProduction.dat"

echo "read (n,p)"
read 3 103 "data/CoHParticleProduction.dat"

echo "read (n,alpha)"
read 3 107 "data/CoHParticleProduction.dat"

echo "read inelastic scattering cross sections"
multiread 3 51 91 "data/CoHLevelExcite1.dat"

echo "construct MF4"
make4

echo "read elastic scattering angular distribution"
angdist 4 2 "data/CoHLegendreCoefficient1.dat"

echo "read inelastic scattering angular distributions"
multiangdist 6 51 90 "data/CoHLegendreCoefficient1.dat"

echo "sum all the partial cross sections"
calc 3 = 4 +  16
calc 3 = 3 +  22
calc 3 = 3 +  28
calc 3 = 3 + 102
calc 3 = 3 + 103
calc 3 = 3 + 107

echo "re-calculate elastic by subtracting non-elastic from total"
calc 2 = 1 - 3

echo "re-calculate total again"
calc 1 = 2 + 3
EOF


echo "produce MF6, and store in $work6"
cat /dev/null > $workmf6
for mt in 16  22  28  102 103 107;do
    echo "make MF6 $mt from data/CoHSpectrum.dat"
    $toolbase/$decemf6 -t $mt -e data/CoHSpectrum.dat -f $workfile >> $workmf6
done


echo "produce MF12 and MF14, and store in $workmf12"
$toolbase/$decemf12 data/gbranch.dat $workfile > $workmf12


echo "assemble all results, $workfile, $workmf6, and $workmf12"
$dece -o $outfile $workfile <<EOF
echo "read MF6"
libread 6 16 "$workmf6"
libread 6 22 "$workmf6"
libread 6 28 "$workmf6"
libread 6 102 "$workmf6"
libread 6 103 "$workmf6"
libread 6 107 "$workmf6"

echo "read MF12"
multilibread 12 51 90 "$workmf12"

echo "read MF14"
multilibread 14 51 90 "$workmf12"

EOF
