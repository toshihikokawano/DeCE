#!/bin/sh

if [ $# -ne 1 ]; then
    echo 'decetest.sh N'
    echo '   where N = 1 ... 13'
    exit 1
fi


testcase=$1
source parameter.sh

###################################################
if [ $testcase = "1" ]; then
    echo '# DESCRIPTION:'
    echo '# Use a template file, which contains resonance parameters only'
    echo '# Read in external data from files and create a new ENDF file'
    echo '# Extract the created subsections'
    echo ''
    $dece $deceoption -o $outbase$testcase.evl $template < $inbase$testcase.dece
fi

###################################################
if [ $testcase = "2" ]; then
    echo '# DESCRIPTION:'
    echo '# Read two levels of inelastic scattering cross sections'
    echo '# and create sum of them in MT3'
    echo ''
    $dece $deceoption -o $outbase$testcase.evl $template < $inbase$testcase.dece
fi

###################################################
if [ $testcase = "3" ]; then
    echo '# DESCRIPTION:'
    echo '# Read all inelastic scattering cross sections'
    echo '# Create total inelastic scattering cross section (MT4)'
    echo ''
    $dece $deceoption -o $outbase$testcase.evl $template < $inbase$testcase.dece
fi

###################################################
if [ $testcase = "4" ]; then
    echo '# DESCRIPTION:'
    echo '# Read total cross section, and copy to elastic'
    echo '# Divide elastic scattering cross section by 2'
    echo ''
    $dece $deceoption -o $outbase$testcase.evl $template < $inbase$testcase.dece
fi

###################################################
if [ $testcase = "5" ]; then
    echo '# DESCRIPTION:'
    echo '# Delete several sections, or whole MF'
    echo ''
    $dece $deceoption -o $outbase$testcase.evl $template < $inbase$testcase.dece
fi

###################################################
if [ $testcase = "6" ]; then
    echo '# DESCRIPTION:'
    echo '# Read angular distribution files and create MF4 and MF6'
    echo ''
    $dece $deceoption -o $outbase$testcase.evl $template < $inbase$testcase.dece
fi

###################################################
if [ $testcase = "7" ]; then
    echo '# DESCRIPTION:'
    echo '# Import data from an ENDF-6 formatted file'
    echo '# Then remove data points above 20 MeV'
    echo ''
    $dece $deceoption -o $outbase$testcase.evl $template < $inbase$testcase.dece
fi

###################################################
if [ $testcase = "8" ]; then
    echo '# DESCRIPTION:'
    echo '# Read (n,p) cross section from a file'
    echo '# Copy into temopral MT95 to MT99 sections'
    echo '# Distort the excitation function by Fermi-function and Gaussian'
    echo '# Multiply by a factor, and normalized to a data point at fixed energy'
    echo ''
    $dece $deceoption -o $outbase$testcase.evl $template < $inbase$testcase.dece
fi

###################################################
if [ $testcase = "9" ]; then
    echo '# DESCRIPTION:'
    echo '# Read inelastic scattering cross sections from a file'
    echo '# Create total inelastic scattering, MT4'
    echo '# Apply the Fermi-function to MT4,'
    echo '# and re-distribute the change to all partial cross sections, MT51 - 91'
    echo ''
    $dece $deceoption -o $outbase$testcase.evl $template < $inbase$testcase.dece
fi

###################################################
if [ $testcase = "10" ]; then
    echo '# DESCRIPTION:'
    echo '# Read total cross section from an evaluated file'
    echo '# Change the interpolation law '
    echo ''
    $dece $deceoption -o $outbase$testcase.evl $template < $inbase$testcase.dece
fi

###################################################
if [ $testcase = "11" ]; then
    echo '# DESCRIPTION:'
    echo '# Read inelastic scattering cross sections from a file'
    echo '# Create total inelastic scattering, MT4'
    echo '# Check threshold energy of MT4'
    echo ''
    $dece $deceoption -o $outbase$testcase.evl $template < $inbase$testcase.dece
fi

###################################################
if [ $testcase = "12" ]; then
    echo '# DESCRIPTION:'
    echo '# Read file that contains MF1 MT455 and 456 only'
    echo '# Create total nu-bar, MT452'
    echo ''
    $dece $deceoption -o $outbase$testcase.evl ENDFMF1MT455.evl < $inbase$testcase.dece
fi

###################################################
if [ $testcase = "13" ]; then
    echo '# DESCRIPTION:'
    echo '# Read elastic scattering angular distribution'
    echo '# First, print Legendre coefficients,'
    echo '# then differential cross section at equi-angle points'
    $dece $deceoption -o $outbase$testcase.evl $template < $inbase$testcase.dece
fi

