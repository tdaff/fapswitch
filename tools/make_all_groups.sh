#!/bin/bash

##
# Demonstration script to generate the basic library of
# functional groups. Note that the actual library was
# hand generated and these will geometries will be different
##


PROG=./align_functional_group.py

#rm *.lib *.svg *.mol

set -o noglob

$PROG -t -n Hydrogen -s H -m -c Alkyl [H]
$PROG -t -n Fluoro -s F -m -c Halide F
$PROG -t -n Chloro -s Cl -m -c Halide Cl
$PROG -t -n Bromo -s Br -m -c Halide Br
$PROG -t -n Iodo -s I -m -c Halide I
$PROG -t -n Amine -s NH2 -m -c "N Containing" N
$PROG -t -n Methyl -s Me -m -c "Alkyl" C
$PROG -t -n Ethyl -s Et -m -c "Alkyl" CC
$PROG -t -n Propyl -s Pr -m -c "Alkyl" CCC
$PROG -t -n Nitro -s NO2 -m -c "N Containing" [N+]\(=O\)[O-]
$PROG -t -n Aldehyde -s HCO -m -c "O Containing" C=O
$PROG -t -n Carboxylate -s COOH -m -c "O Containing" C\(=O\)O
$PROG -t -n Alcohol -s OH -m -c "O Containing" O
$PROG -t -n Nitrile -s CN -m -c "N Containing" C#N
$PROG -t -n Phenyl -s Ph -m -c "Aromatic" c1ccccc1
$PROG -t -n Methoxy -s OMe -m -c "Ether" OC
$PROG -t -n Ethoxy -s OEt -m -c "Ether" OCC
$PROG -t -n PropylEther -s OPr -m -c "Ether" OCCC
$PROG -t -n MethylAmine -s NHMe -m -c "N Containing" NC
$PROG -t -n SulphonicAcid -s SO3H -m -c "S Containing" S\(=O\)\(=O\)O
$PROG -t -n Ethene -s CHCH2 -m -c "Hydrocarbyl" C=C
$PROG -t -n Ethyne -s CCH -m -c "Hydrocarbyl" C#C
$PROG -t -n PendantMethylAmine -s MeNH2 -m -c "N Containing" CN
$PROG -t -n Oxyethene -s OEte -m -c "O Containing" OC=C
$PROG -t -n Oxypropene -s OPre -m -c "O Containing" OCC=C
$PROG -t -n TriFluoroMethyl -s CF3 -m -c "Alkyl Halide" C\(F\)\(F\)F
$PROG -t -n Amide -s CONH2 -m -c "N Containing" C\(=O\)N
$PROG -t -n Ethanimine -s CHNH -m -c "N Containing" C=N

