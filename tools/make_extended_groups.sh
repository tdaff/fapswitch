#!/bin/bash

##
# Generate some common and more exotic functional groups
# from their smiles. These are included in the
# library 
##


PROG=./align_functional_group.py

rm *.flib *.svg *.mol *.html

set -o noglob

$PROG -t -n fluorocarbonyl -s COF -m C\(=O\)F
$PROG -t -n chlorocarbonyl -s COCl -m C\(=O\)Cl
$PROG -t -n bromocarbonyl -s COBr -m C\(=O\)Br
$PROG -t -n iodocarbonyl -s COI -m C\(=O\)I
$PROG -t -n fluoromethyl -s MeF -m CF
$PROG -t -n chloromethyl -s MeCl -m CCl
$PROG -t -n bromomethyl -s MeBr -m CBr
$PROG -t -n iodomethyl -s MeI -m CI
$PROG -t -n difluoromethyl -s MeF2 -m C\(F\)F
$PROG -t -n dichloromethyl -s MeCl2 -m C\(Cl\)Cl
$PROG -t -n dibromomethyl -s MeBr2 -m C\(Br\)Br
$PROG -t -n diiodomethyl -s MeI2 -m C\(I\)I
$PROG -t -n trichloromethyl -s MeCl3 -m C\(Cl\)\(Cl\)Cl
$PROG -t -n tribromomethyl -s MeBr3 -m C\(Br\)\(Br\)Br
$PROG -t -n triiodomethyl -s MeI3 -m C\(I\)\(I\)I

# Boron containing

$PROG -t -n 2H-azaborinin-1-yl -s AzBorn1 -m N1C=CC=CB1
$PROG -t -n 1,2-dihydroazaborinin-6-yl -s AzBorn6 C1=CC=CBN1
$PROG -t -n 1,2-dihydroazaborinin-5-yl -s AzBorn5 C1=CNBC=C1
$PROG -t -n 1,2-dihydroazaborinin-4-yl -s AzBorn4 C1=CBNC=C1
$PROG -t -n 1,2-dihydroazaborinin-3-yl -s AzBorn3 C1=CC=CNB1
$PROG -t -n 1H-azaborinin-2-yl -s AzBorn2 -m B1C=CC=CN1
$PROG -t -n boronooxy -s Boox -m OB\(O\)O
$PROG -t -n borol-1-yl -s Borl1 -m B1C=CC=C1
$PROG -t -n 1H-borol-2-yl -s Borl2 -m C1=CC=CB1
$PROG -t -n 1H-borol-3-yl -s Borl3 -m C1=CBC=C1
$PROG -t -n borono -s BOH2 -m B\(O\)O
$PROG -t -n dimethoxyboranyloxy -s DMeBanox -m OB\(OC\)\(OC\)
$PROG -t -n [hydroxy\(methoxy\)boranyl]oxy -s MeBanox -m OB\(O\)OC
$PROG -t -n hydroxy\(methoxy\)boranyl -s BOHOMe -m B\(OC\)O
$PROG -t -n methoxycarbonyl -s COOMe -m C\(=O\)OC
$PROG -t -n acetyl -s COMe -m C\(=O\)C
$PROG -t -n cyanato -s OCN OC#N
$PROG -t -n isocyanato -s NCO -m N=C=O
$PROG -t -n isothiocyanato -s NCS -m N=C=S
$PROG -t -n thiocyanato -s SCN -m SC#N
$PROG -t -n 1-methylenepropyl -s MeePr -m C\(=C\)CC
$PROG -t -n 1-methylprop-1-enyl -s MePre -m C\(C\)=CC
$PROG -t -n isopropenyl -s iPre -m C\(=C\)C
$PROG -t -n 1-methylallyl -s MeAll -m C\(C\)C=C
$PROG -t -n isopropyl -s iPr -m C\(C\)C

# Nitrogen aromatic

$PROG -t -n triazol-1-yl -s Trzl1 N1N=NC=C1
$PROG -t -n 1H-triazol-4-yl -s Trzl4 C1=CNN=N1
$PROG -t -n 1H-triazol-5-yl -s Trzl5 C1=CN=NN1
$PROG -t -n 1,2,4-oxadiazol-5-yl -s 124Oxdzl5 -m C1=NC=NO1
$PROG -t -n 1,2,4-oxadiazol-3-yl -s 124Oxdzl3 -m C1=NOC=N1
$PROG -t -n 1H-1,2,4-triazol-3-yl -s 124Trzl3 C1=NNC=N1
$PROG -t -n 1H-1,2,4-triazol-5-yl -s 124Trzl5 C1=NC=NN1
$PROG -t -n 1,2,4-triazol-1-yl -s 124Trzl1 N1N=CN=C1
$PROG -t -n 5-amino-1,3,4-oxadiazol-2-yl -s Am134Oxdzl C1=NN=C\(N\)O1
$PROG -t -n 1,3,4-oxadiazol-2-yl -s 134Oxdzl -m C1=NN=CO1
$PROG -t -n 4,6-diamino-1,3,5-triazin-2-yl -s DAmTzn c1nc\(N\)nc\(N\)n1
$PROG -t -n imidazol-1-yl -s Imidaz1 N1C=NC=C1
$PROG -t -n 1H-imidazol-2-yl -s Imidaz2 C1=NC=CN1
$PROG -t -n 1H-imidazol-4-yl -s Imidaz4 C1=CNC=N1
$PROG -t -n 1H-imidazol-5-yl -s Imidaz5 C1=CN=CN1
$PROG -t -n isoxazol-3-yl -s IsOxzl3 -m C1=NOC=C1
$PROG -t -n isoxazol-4-yl -s IsOxzl4 -m C1=CON=C1
$PROG -t -n isoxazol-5-yl -s IsOxzl5 -m C1=CC=NO1
$PROG -t -n isoxazolidin-2-yl -s IsOxzldn2 -m N1CCCO1
$PROG -t -n isoxazolidin-3-yl -s IsOxzldn3 -m C1NOCC1
$PROG -t -n isoxazolidin-4-yl -s IsOxzldn4 -m C1CONC1
$PROG -t -n isoxazolidin-5-yl -s IsOxzldn5 -m C1CCNO1
$PROG -t -n 5-nitroimidazol-1-yl -s NO2Imidaz1 N1C\([N+]\([O-]\)=O\)=CN=C1
$PROG -t -n 5-nitro-1H-imidazol-2-yl -s NO2Imidaz2 C1=NC=C\([N+]\([O-]\)=O\)N1
$PROG -t -n 5-nitro-1H-imidazol-4-yl -s NO2Imidaz4 C1=C\([N+]\([O-]\)=O\)NC=N1
$PROG -t -n oxazol-2-yl -s Oxzl2 -m C1=NC=CO1
$PROG -t -n oxazol-4-yl -s Oxzl4 -m C1=COC=N1
$PROG -t -n oxazol-5-yl -s Oxzl5 -m C1=CN=CO1
$PROG -t -n oxazolidin-2-yl -s Oxzldn2 -m C1OCCN1
$PROG -t -n oxazolidin-3-yl -s Oxzldn3 -m N1CCOC1
$PROG -t -n oxazolidin-4-yl -s Oxzldn4 -m C1NCOC1
$PROG -t -n oxazolidin-5-yl -s Oxzldn5 -m C1CNCO1
$PROG -t -n 4,5-dihydrooxazol-2-yl -s DHOxzl2 -m C1=NCCO1
$PROG -t -n 4,5-dihydrooxazol-4-yl -s DHOxzl4 -m C1N=COC1
$PROG -t -n 4,5-dihydrooxazol-5-yl -s DHOxzl5 -m C1OC=NC1
$PROG -t -n pyrrol-1-yl -s Pyrl1 -m N1C=CC=C1
$PROG -t -n 1H-pyrrol-2-yl -s Pyrl2 -m C1=CC=CN1
$PROG -t -n 1H-pyrrol-3-yl -s Pyrl3 -m C1=CNC=C1
$PROG -t -n 1,3,5-triazin-2-yl -s TAzn c1ncncn1

# Nitrogen containing

$PROG -t -n [4-hydroxy-4-oxo-but-2-enoyl]amino -s AmMale -m NC\(=O\)C=CC\(=O\)O
$PROG -t -n dicarboxyamino -s DCoxAm N\(C\(=O\)O\)C\(=O\)O
$PROG -t -n aziridin-1-yl -s AzIrdn -m N1CC1
$PROG -t -n diazenyl -s NNH -m N=N
$PROG -t -n dimethylamino -s NMe2 -m N\(C\)C
$PROG -t -n 1-hydroxyethylamino -s EtOHAm -m NC\(O\)C
$PROG -t -n methylazo -s NNMe -m N=NC
$PROG -t -n nitrooxy -s NO3 O[N+]\([O-]\)\(=O\)
$PROG -t -n nitroso -s NO N=O

# Oxygen containing

$PROG -t -n carboxyoxy -s OCOOH OC\(=O\)O
$PROG -t -n hydroperoxy -s OOH -m OO
$PROG -t -n methoxycarbonyloxy -s OCOOMe -m OC\(=O\)OC
$PROG -t -n methylperoxy -s OOMe -m OOC
$PROG -t -n nitroperoxy -s OONOO OO[N+]\(=O\)[O-]
$PROG -t -n nitrosooxy -s ONO ON=O
$PROG -t -n prop-1-enoxy -s OPr1e -m OC=CC
#$PROG -t -n propoxy -s OPr -m OCCC

# Phosphorous containing

$PROG -t -n dimethylphosphanyl -s PMe2 -m P\(C\)C
$PROG -t -n methylphosphanyl -s PHMe -m PC
$PROG -t -n [hydroxy\(methoxy\)phosphoryl]oxy -s MeHPO4 -m OP\(O\)\(=O\)OC
$PROG -t -n phosphonooxy -s PO4H2 OP\(O\)\(=O\)O
$PROG -t -n phosphono -s PO3H2 P\(O\)\(O\)=O

# Sulphur containing

$PROG -t -n disulfanyl -s SSH -m SS
$PROG -t -n methyldisulfanyl -s SSMe -m SSC
$PROG -t -n methylsulfanyl -s SMe -m SC
$PROG -t -n methylsulfinyl -s SOMe S\(=O\)C
$PROG -t -n ethanethioyl -s CSMe -m C\(=S\)C
$PROG -t -n aminosulfinyl -s SONH2 S\(=O\)N
$PROG -t -n hydroxysulfinyl -s SOOH S\(=O\)O
$PROG -t -n hydrosulfinyl -s SHO S\([H]\)\(=O\)
$PROG -t -n sulfamoyl -s SO2NH2 S\(=O\)\(=O\)N
$PROG -t -n methanethioyl -s CHS -m C=S
$PROG -t -n thiazol-1-yl -s Thzl1 -m S1C=NC=C1
$PROG -t -n thiazol-2-yl -s Thzl2 -m C1=NC=CS1
$PROG -t -n thiazol-4-yl -s Thzl4 -m C1=CSC=N1
$PROG -t -n thiazol-5-yl -s Thzl5 -m C1=CN=CS1
$PROG -t -n thiazol-3-ium-3-yl -s Thzlum3 -m [N+]1=CSC=C1
$PROG -t -n sulfanyl -s SH S\([H]\)

# Mixed groups

$PROG -t -n 2-bromo-3-furyl -s 2Br3Fur -m c1c\(Br\)occ1
$PROG -t -n 5-bromo-3-furyl -s 5Br3Fur -m c1cc\(Br\)oc1
$PROG -t -n 5-bromo-2-furyl -s 5Br2Fur -m c1ccc\(Br\)o1
$PROG -t -n 2-bromo-3-thienyl -s 2Br3Thi -m c1c\(Br\)scc1
$PROG -t -n 5-bromo-3-thienyl -s 5Br3Thi -m c1cc\(Br\)sc1
$PROG -t -n 5-bromo-2-thienyl -s 5Br2Thi -m c1ccc\(Br\)s1
$PROG -t -n 2-chloro-3-furyl -s 2Cl3Fur -m c1c\(Cl\)occ1
$PROG -t -n 5-chloro-3-furyl -s 5Cl3Fur -m c1cc\(Cl\)oc1
$PROG -t -n 5-chloro-2-furyl -s 5Cl2Fur -m c1ccc\(Cl\)o1
$PROG -t -n 2-chloro-3-thienyl -s 2Cl3Thi -m c1c\(Cl\)scc1
$PROG -t -n 5-chloro-3-thienyl -s 5Cl3Thi -m c1cc\(Cl\)sc1
$PROG -t -n 5-chloro-2-thienyl -s 5Cl2Thi -m c1ccc\(Cl\)s1
$PROG -t -n 2-ethyl-3-furyl -s 2Et3Fur -m c1c\(CC\)occ1
$PROG -t -n 5-ethyl-3-furyl -s 5Et3Fur -m c1cc\(CC\)oc1
$PROG -t -n 5-ethyl-2-furyl -s 5Et2Fur -m c1oc\(CC\)cc1
$PROG -t -n 2-ethyl-3-thienyl -s 2Et3Thi -m c1c\(CC\)scc1
$PROG -t -n 5-ethyl-3-thienyl -s 5Et3Thi -m c1cc\(CC\)sc1
$PROG -t -n 5-ethyl-2-thienyl -s 5Et2Thi -m c1sc\(CC\)cc1
$PROG -t -n 2-fluoro-3-furyl -s 2F3Fur -m c1c\(F\)occ1
$PROG -t -n 5-fluoro-3-furyl -s 5F3Fur -m c1cc\(F\)oc1
$PROG -t -n 5-fluoro-2-furyl -s 5F2Fur -m c1ccc\(F\)o1
$PROG -t -n 2-fluoro-3-thienyl -s 2F3Thi -m c1c\(F\)scc1
$PROG -t -n 5-fluoro-3-thienyl -s 5F3Thi -m c1cc\(F\)sc1
$PROG -t -n 5-fluoro-2-thienyl -s 5F2Thi -m c1ccc\(F\)s1
$PROG -t -n 2-methyl-3-furyl -s 2Me3Fur -m c1c\(C\)occ1
$PROG -t -n 5-methyl-3-furyl -s 5Me3Fur -m c1cc\(C\)oc1
$PROG -t -n 5-methyl-2-furyl -s 5Me2Fur -m c1oc\(C\)cc1
$PROG -t -n 2-methyl-3-thienyl -s 2Me3Thi -m c1c\(C\)scc1
$PROG -t -n 5-methyl-3-thienyl -s 5Me3Thi -m c1cc\(C\)sc1
$PROG -t -n 5-methyl-2-thienyl -s 5Me2Thi -m c1sc\(C\)cc1
$PROG -t -n 3-bromo-2-furyl -s 3Br2Fur -m c1occc1\(Br\)
$PROG -t -n 4-bromo-3-furyl -s 4Br3Fur -m c1c\(Br\)coc1
$PROG -t -n 4-bromo-2-furyl -s 4Br2Fur -m c1cc\(Br\)co1
$PROG -t -n 3-bromo-2-thienyl -s 3Br2Thi -m c1sccc1\(Br\)
$PROG -t -n 4-bromo-3-thienyl -s 4Br3Thi -m c1c\(Br\)csc1
$PROG -t -n 4-bromo-2-thienyl -s 4Br2Thi -m c1cc\(Br\)cs1
$PROG -t -n 3-chloro-2-furyl -s 3Cl2Fur -m c1occc1\(Cl\)
$PROG -t -n 4-chloro-3-furyl -s 4Cl3Fur -m c1c\(Cl\)coc1
$PROG -t -n 4-chloro-2-furyl -s 4Cl2Fur -m c1cc\(Cl\)co1
$PROG -t -n 3-chloro-2-thienyl -s 3Cl2Thi -m c1sccc1\(Cl\)
$PROG -t -n 4-chloro-3-thienyl -s 4Cl3Thi -m c1c\(Cl\)csc1
$PROG -t -n 4-chloro-2-thienyl -s 4Cl2Thi -m c1cc\(Cl\)cs1
$PROG -t -n 3-fluoro-2-furyl -s 3F2Fur -m c1occc1\(F\)
$PROG -t -n 4-fluoro-3-furyl -s 4F3Fur -m c1c\(F\)coc1
$PROG -t -n 4-fluoro-2-furyl -s 4F2Fur -m c1cc\(F\)co1
$PROG -t -n 3-fluoro-2-thienyl -s 3F2Thi -m c1sccc1\(F\)
$PROG -t -n 4-fluoro-3-thienyl -s 4F3Thi -m c1c\(F\)csc1
$PROG -t -n 4-fluoro-2-thienyl -s 4F2Thi -m c1cc\(F\)cs1
$PROG -t -n N-hydroxy-C-methyl-carbonimidoyl -s AcAdxm -m C\(C\)=NO
$PROG -t -n hydroxyiminomethyl -s Adxm -m C=NO
$PROG -t -n N-hydroxycarbamimidoyl -s AmAdxm -m C\(N\)=NO
$PROG -t -n C-fluoro-N-hydroxy-carbonimidoyl -s FAdxm -m C\(F\)=NO
$PROG -t -n C-chloro-N-hydroxy-carbonimidoyl -s ClAdxm -m C\(Cl\)=NO
$PROG -t -n C-bromo-N-hydroxy-carbonimidoyl -s BrAdxm -m C\(Br\)=NO
$PROG -t -n C-ethyl-N-hydroxy-carbonimidoyl -s PrAdxm -m C\(CC\)=NO
$PROG -t -n 2-\(2-furyl\)-3-furyl -s 22Fur3Fur -m c1c\(c2ccco2\)occ1
$PROG -t -n 5-\(2-furyl\)-3-furyl -s 52Fur3Fur -m c1coc\(c2ccco2\)c1
$PROG -t -n 5-\(2-furyl\)-2-furyl -s 52Fur2Fur -m c1ccc\(c2ccco2\)o1
$PROG -t -n 2-\(2-thienyl\)-3-thienyl -s 22Thi3Thi -m c1c\(c2cccs2\)scc1
$PROG -t -n 5-\(2-thienyl\)-3-thienyl -s 52Thi3Thi -m c1csc\(c2cccs2\)c1
$PROG -t -n 5-\(2-thienyl\)-2-thienyl -s 52Thi2Thi -m c1ccc\(c2cccs2\)s1
$PROG -t -n 2-\(5-bromo-2-furyl\)-3-furyl -s 25Br2Fur3Fur -m c1c\(c2ccc\(Br\)o2\)occ1
$PROG -t -n 5-\(5-bromo-2-furyl\)-3-furyl -s 55Br2Fur3Fur -m c1cc\(c2ccc\(Br\)o2\)oc1
$PROG -t -n 5-\(5-bromo-2-furyl\)-2-furyl -s 55Br2Fur2Fur -m c1ccc\(c2ccc\(Br\)o2\)o1
$PROG -t -n 2-\(5-bromo-2-thienyl\)-3-thienyl -s 25Br2Thi3Thi -m c1c\(c2ccc\(Br\)s2\)scc1
$PROG -t -n 5-\(5-bromo-2-thienyl\)-3-thienyl -s 55Br2Thi3Thi -m c1cc\(c2ccc\(Br\)s2\)sc1
$PROG -t -n 5-\(5-bromo-2-thienyl\)-2-thienyl -s 55Br2Thi2Thi -m c1ccc\(c2ccc\(Br\)s2\)s1
$PROG -t -n 2-furyl -s 2Fur -m C1=CC=CO1
$PROG -t -n 3-furyl -s 3Fur -m C1=COC=C1
$PROG -t -n 2-amino-3-furyl -s 2Am3Fur -m c1c\(N\)occ1
$PROG -t -n 5-amino-3-furyl -s 5Am3Fur -m c1coc\(N\)c1
$PROG -t -n 5-amino-2-furyl -s 5Am2Fur -m C1=CC=C\(N\)O1
$PROG -t -n pyrazin-2-yl -s Pyraz -m c1nccnc1
$PROG -t -n hydroxysulfanyl -s SOH -m SO
$PROG -t -n 2-thienyl -s 2Thi -m C1=CC=CS1
$PROG -t -n 3-thienyl -s 3Thi -m C1=CSC=C1
$PROG -t -n 2-amino-3-thienyl -s 2Am3Thi -m c1c\(N\)scc1
$PROG -t -n 5-amino-3-thienyl -s 5Am3Thi -m c1csc\(N\)c1
$PROG -t -n 5-amino-2-thienyl -s 5Am2Thi -m C1=CC=C\(N\)S1

# Most common from Ertl, not duplicating already included groups

$PROG -t -n 2,4-dichlorophenyl -s 24DClPh -m c1c\(Cl\)cc\(Cl\)cc1
$PROG -t -n butyl -s Bu -m CCCC
$PROG -t -n isobutyl -s iBu -m CC\(C\)C
$PROG -t -n allyl -s All -m CC=C
$PROG -t -n acetoxy -s OAc -m OC\(=O\)C
$PROG -t -n acetamido -s NCOMe NC\(=O\)C
$PROG -t -n anilino -s NPh -m Nc1ccccc1
$PROG -t -n benzyloxy -s OBnz -m OCc1ccccc1
$PROG -t -n benzyl -s Bnz -m Cc1ccccc1
$PROG -t -n 4-bromophenyl -s 4BrPh -m c1ccc\(Br\)cc1
$PROG -t -n 2-chlorophenyl -s 2ClPh -m c1c\(Cl\)cccc1
$PROG -t -n 3-chlorophenyl -s 3ClPh -m c1cc\(Cl\)ccc1
$PROG -t -n 4-chlorophenyl -s 4ClPh -m c1ccc\(Cl\)cc1
$PROG -t -n cyclohexyl -s CHex -m C1CCCCC1
$PROG -t -n ethoxycarbonyl -s COOEt -m C\(=O\)OCC
$PROG -t -n 4-fluorophenyl -s 4FPh -m c1ccc\(F\)cc1
$PROG -t -n hydroxymethyl -s COH -m CO
$PROG -t -n 4-methoxyphenyl -s 4OMePh -m c1ccc\(OC\)cc1
$PROG -t -n o-tolyl -s oTol -m c1c\(C\)cccc1
$PROG -t -n p-tolyl -s pTol -m c1ccc\(C\)cc1
$PROG -t -n morpholino -s Morph -m N1CCOCC1
$PROG -t -n 4-nitrophenyl -s 4NO2Ph -m c1ccc\([N+]\(=O\)[O-]\)cc1
$PROG -t -n benzoyl -s COPh -m C\(=O\)c1ccccc1
$PROG -t -n tert-butyl -s tBu -m C\(C\)\(C\)\(C\)
$PROG -t -n 3,4,5-trihydroxy-6-\(hydroxymethyl\)tetrahydropyran-2-yl -s AGlutol -m C1O[C@@H]\(CO\)[C@H]\(O\)[C@@H]\(O\)[C@@H]1\(O\)
$PROG -t -n 3-guanidinopropyl -s PrNCNN CCCNC\(=N\)N
$PROG -t -n \(4-hydroxyphenyl\)methyl -s CPh4OH -m Cc1ccc\(O\)cc1
$PROG -t -n carboxymethyl -s AcOH -m CC\(=O\)O
$PROG -t -n sec-butyl -s sBu -m C\(C\)CC
$PROG -t -n diethylamino -s NEt2 -m N\(CC\)CC
$PROG -t -n 1-hydroxyethyl -s 1OHEt -m C\(O\)C
$PROG -t -n carbamimidoyl -s Amidn C\(N\)=N
$PROG -t -n guanidino -s NCNN NC\(=N\)N
$PROG -t -n 1H-indol-3-yl -s 3Indl -m c1c\(C=CC=C2\)c2nc1
$PROG -t -n \(2-amino-2-oxo-ethyl\) -s MeCONH2 CC\(=O\)N
$PROG -t -n 4-aminobutyl -s BuNH2 -m CCCCN
$PROG -t -n 2-carboxyethyl -s PrOOH -m CCC\(=O\)O


