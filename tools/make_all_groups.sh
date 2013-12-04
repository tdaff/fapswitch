#!/bin/bash

##
# Demonstration script to generate the basic library of
# functional groups. Note that the actual library was
# hand generated and these will geometries will be different
##


PROG=./align_functional_group.py

#rm *.lib *.svg *.mol

set -o noglob

#$PROG -t -n Hydrogen -s H -m -c Alkyl [H]
$PROG -t -n fluoro -s F -m -c Halide F
$PROG -t -n chloro -s Cl -m -c Halide Cl
$PROG -t -n bromo -s Br -m -c Halide Br
$PROG -t -n iodo -s I -m -c Halide I
$PROG -t -n amino -s NH2 -m -c "N Containing" N
$PROG -t -n methyl -s Me -m -c "Alkyl" C
$PROG -t -n ethyl -s Et -m -c "Alkyl" CC
$PROG -t -n propyl -s Pr -m -c "Alkyl" CCC
$PROG -t -n nitro -s NO2 -m -c "N Containing" [N+]\(=O\)[O-]
$PROG -t -n formyl -s HCO -m -c "O Containing" C=O
$PROG -t -n carboxy -s COOH -m -c "O Containing" C\(=O\)O
$PROG -t -n hydroxy -s OH -m -c "O Containing" O
$PROG -t -n cyano -s CN -m -c "N Containing" C#N
$PROG -t -n phenyl -s Ph -m -c "Aromatic" c1ccccc1
$PROG -t -n methoxy -s OMe -m -c "Ether" OC
$PROG -t -n ethoxy -s OEt -m -c "Ether" OCC
$PROG -t -n propoxy -s OPr -m -c "Ether" OCCC
$PROG -t -n methylamino -s NHMe -m -c "N Containing" NC
$PROG -t -n sulfo -s SO3H -m -c "S Containing" S\(=O\)\(=O\)O
$PROG -t -n vinyl -s CHCH2 -m -c "Hydrocarbyl" C=C
$PROG -t -n ethynyl -s CCH -m -c "Hydrocarbyl" C#C
$PROG -t -n aminomethyl -s MeNH2 -m -c "N Containing" CN
$PROG -t -n vinyloxy -s OEte -m -c "O Containing" OC=C
$PROG -t -n allyloxy -s OPre -m -c "O Containing" OCC=C
$PROG -t -n trifluoromethyl -s CF3 -m -c "Alkyl Halide" C\(F\)\(F\)F
$PROG -t -n carbamoyl -s CONH2 -m -c "N Containing" C\(=O\)N
$PROG -t -n methanimidoyl -s CHNH -m -c "N Containing" C=N

