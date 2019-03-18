#!/bin/sh

if [ $# -ne 1 ]; then
    echo 'decetest.sh N'
    echo '   where N = 1 ... 16'
    exit 1
fi


testcase=$1
source parameter.sh

echo '# TEST CASE: ' $testcase

###################################################
if [ $testcase = "1" ]; then
    echo '# Use a template file, which contains resonance parameters only'
    echo '# Read in external data from files and create a new ENDF file'
    echo '# Extract the created subsections'

###################################################
elif [ $testcase = "2" ]; then
    echo '# Read two levels of inelastic scattering cross sections'
    echo '# and create sum of them in MT3'

###################################################
elif [ $testcase = "3" ]; then
    echo '# Read all inelastic scattering cross sections'
    echo '# Create total inelastic scattering cross section (MT4)'

###################################################
elif [ $testcase = "4" ]; then
    echo '# Read total cross section, and copy to elastic'
    echo '# Divide elastic scattering cross section by 2'

###################################################
elif [ $testcase = "5" ]; then
    echo '# Delete several sections, or whole MF'

###################################################
elif [ $testcase = "6" ]; then
    echo '# Read angular distribution files and create MF4 and MF6'

###################################################
elif [ $testcase = "7" ]; then
    echo '# Import data from an ENDF-6 formatted file'
    echo '# Then remove data points above 20 MeV'

###################################################
elif [ $testcase = "8" ]; then
    echo '# Read (n,p) cross section from a file'
    echo '# Copy into temopral MT95 to MT99 sections'
    echo '# Distort the excitation function by Fermi-function and Gaussian'
    echo '# Multiply by a factor, and normalized to a data point at fixed energy'

###################################################
elif [ $testcase = "9" ]; then
    echo '# Read inelastic scattering cross sections from a file'
    echo '# Create total inelastic scattering, MT4'
    echo '# Apply the Fermi-function to MT4,'
    echo '# and re-distribute the change to all partial cross sections, MT51 - 91'

###################################################
elif [ $testcase = "10" ]; then
    echo '# Read total cross section from an evaluated file'
    echo '# Change the interpolation law '

###################################################
elif [ $testcase = "11" ]; then
    echo '# Read inelastic scattering cross sections from a file'
    echo '# Create total inelastic scattering, MT4'
    echo '# Check threshold energy of MT4'

###################################################
elif [ $testcase = "12" ]; then
    echo '# Read file that contains MF1 MT455 and 456 only'
    echo '# Create total nu-bar, MT452'
    template=ENDFMF1MT455.evl

###################################################
elif [ $testcase = "13" ]; then
    echo '# Read elastic scattering angular distribution'
    echo '# First, print Legendre coefficients,'
    echo '# then differential cross section at equi-angle points'

###################################################
elif [ $testcase = "14" ]; then
    echo '# Read two MF3 sections from file'
    echo '# then calculate the ratio and product of these sections'

###################################################
elif [ $testcase = "15" ]; then
    echo '# Reconstruct pointwise cross section'
    echo '# from resolved resonance parameters'
    echo '# This may take some time'

###################################################
elif [ $testcase = "16" ]; then
    echo '# Reconstruct scattering angular distributions'
    echo '# at every 1 keV interval from resolved resonance parameters'

###################################################
elif [ $testcase = "17" ]; then
    echo '# Reconstruct elastic scattering cross section'
    echo '# by summing all the partial cross sections, and'
    echo '# subtract it from the total cross section'

###################################################
else
    echo '# currently testcase shoud be less than 17'
    exit;
fi

$dece $deceoption -o $outbase/ENDFout$testcase.evl $template < $inbase/input$testcase.dece > $outbase/ENDFout$testcase.out 2> $outbase/ENDFout$testcase.log
cat $outbase/ENDFout$testcase.log
