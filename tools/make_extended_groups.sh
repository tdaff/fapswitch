#!/bin/bash

##
# Generate some common and more exotic functional groups
# from their smiles. These are included in the
# library 
##


PROG=./align_functional_group.py

#rm *.flib *.svg *.mol

set -o noglob

$PROG -t -n AcylFluoride -s COF C\(=O\)F
$PROG -t -n AcylChloride -s COCl C\(=O\)Cl
$PROG -t -n AcylBromide -s COBr C\(=O\)Br
$PROG -t -n AcylIodide -s COI C\(=O\)I
$PROG -t -n MethylFluoride -s MeF CF
$PROG -t -n MethylChloride -s MeCl CCl
$PROG -t -n MethylBromide -s MeBr CBr
$PROG -t -n MethylIodide -s MeI CI
$PROG -t -n DiFluoroMethyl -s MeF2 C\(F\)F
$PROG -t -n DiChloroMethyl -s MeCl2 C\(Cl\)Cl
$PROG -t -n DiBromoMethyl -s MeBr2 C\(Br\)Br
$PROG -t -n DiIodoMethyl -s MeI2 C\(I\)I
$PROG -t -n TriChloroMethyl -s MeCl3 C\(Cl\)\(Cl\)Cl
$PROG -t -n TriBromoMethyl -s MeBr3 C\(Br\)\(Br\)Br
$PROG -t -n TriIodoMethyl -s MeI3 C\(I\)\(I\)I

# Boron containing

$PROG -t -n AzaBorine:1 -s AzBorn:1 N1C=CC=CB1
$PROG -t -n AzaBorine:2 -s AzBorn:2 C1=CC=CBN1
$PROG -t -n AzaBorine:3 -s AzBorn:3 C1=CNBC=C1
$PROG -t -n AzaBorine:4 -s AzBorn:4 C1=CBNC=C1
$PROG -t -n AzaBorine:5 -s AzBorn:5 C1=CC=CNB1
$PROG -t -n AzaBorine:6 -s AzBorn:6 B1C=CC=CN1
$PROG -t -n Borate -s Bate OB\(O\)O
$PROG -t -n Borole:1 -s Borl:1 B1C=CC=C1
$PROG -t -n Borole:2 -s Borl:2 C1=CC=CB1
$PROG -t -n Borole:3 -s Borl:3 C1=CBC=C1
$PROG -t -n Borono -s BOH2 B\(O\)O
$PROG -t -n DiMethylBorate -s DMeBate OB\(OC\)\(OC\)
$PROG -t -n MethylBorate -s MeBate OB\(O\)OC
$PROG -t -n MethylHydroxyBorino -s BOHOMe B\(OC\)O
$PROG -t -n MethylEster -s COOMe C\(=O\)OC
$PROG -t -n MethylKetone -s COMe C\(=O\)C
$PROG -t -n Cyanate -s OCN OC#N
$PROG -t -n IsoCyanate -s NCO N=C=O
$PROG -t -n IsoThioCyanate -s NCS N=C=S
$PROG -t -n ThioCyanate -s SCN SC#N
$PROG -t -n 2But1ene -s 2Bt1e C\(=C\)CC
$PROG -t -n 2But2ene -s 2Bt2e C\(C\)=CC
$PROG -t -n 2Pro1ene -s 2Pr1e C\(=C\)C
$PROG -t -n 3But1ene -s 3Bt1e C\(C\)C=C
$PROG -t -n IsoPropyl -s iPr C\(C\)C

# Nitrogen aromatic

$PROG -t -n 123Triazole:1 -s 123Trzl:1 N1N=NC=C1
$PROG -t -n 123Triazole:4 -s 123Trzl:4 C1=CNN=N1
$PROG -t -n 123Triazole:5 -s 123Trzl:5 C1=CN=NN1
$PROG -t -n 124Oxadiazole:3 -s 124Oxdzl:3 C1=NC=NO1
$PROG -t -n 124Oxadiazole:5 -s 124Oxdzl:5 C1=NOC=N1
$PROG -t -n 124Triazole:1 -s 124Trzl:1 C1=NNC=N1
$PROG -t -n 124Triazole:3 -s 124Trzl:3 C1=NC=NN1
$PROG -t -n 124Triazole:5 -s 124Trzl:5 N1N=CN=C1
$PROG -t -n 134AminoOxadiazole -s 134AmOxdzl C1=NN=C\(N\)O1
$PROG -t -n 134Oxadiazole -s 134Oxdzl C1=NN=CO1
$PROG -t -n DiAminoTriazine -s DAmTzn c1nc\(N\)nc\(N\)n1
$PROG -t -n Imidazole:1 -s Imdaz:1 N1C=NC=C1
$PROG -t -n Imidazole:2 -s Imdaz:2 C1=NC=CN1
$PROG -t -n Imidazole:4 -s Imdaz:4 C1=CNC=N1
$PROG -t -n Imidazole:5 -s Imdaz:5 C1=CN=CN1
$PROG -t -n IsOxazole:3 -s IsOxzl:3 C1=NOC=C1
$PROG -t -n IsOxazole:4 -s IsOxzl:4 C1=CON=C1
$PROG -t -n IsOxazole:5 -s IsOxzl:5 C1=CC=NO1
$PROG -t -n IsOxazolidine:2 -s IsOxzldn:2 N1CCCO1
$PROG -t -n IsOxazolidine:3 -s IsOxzldn:3 C1NOCC1
$PROG -t -n IsOxazolidine:4 -s IsOxzldn:4 C1CONC1
$PROG -t -n IsOxazolidine:5 -s IsOxzldn:5 C1CCNO1
$PROG -t -n NO2Imidazole:1 -s NO2Imdaz:1 N1C\([N+]\([O-]\)=O\)=CN=C1
$PROG -t -n NO2Imidazole:2 -s NO2Imdaz:2 C1=NC=C\([N+]\([O-]\)=O\)N1
$PROG -t -n NO2Imidazole:4 -s NO2Imdaz:4 C1=C\([N+]\([O-]\)=O\)NC=N1
$PROG -t -n Oxazole:2 -s Oxzl:2 C1=NC=CO1
$PROG -t -n Oxazole:4 -s Oxzl:4 C1=COC=N1
$PROG -t -n Oxazole:5 -s Oxzl:5 C1=CN=CO1
$PROG -t -n Oxazolidine:2 -s Oxzldn:2 C1OCCN1
$PROG -t -n Oxazolidine:3 -s Oxzldn:3 N1CCOC1
$PROG -t -n Oxazolidine:4 -s Oxzldn:4 C1NCOC1
$PROG -t -n Oxazolidine:5 -s Oxzldn:5 C1CNCO1
$PROG -t -n Oxazoline:2 -s Oxzln:2 C1=NCCO1
$PROG -t -n Oxazoline:4 -s Oxzln:4 C1N=COC1
$PROG -t -n Oxazoline:5 -s Oxzln:5 C1OC=NC1
$PROG -t -n Pyrrole:1 -s Pyrl:1 N1C=CC=C1
$PROG -t -n Pyrrole:2 -s Pyrl:2 C1=CC=CN1
$PROG -t -n Pyrrole:3 -s Pyrl:3 C1=CNC=C1
$PROG -t -n TriAzine -s TAzn c1ncncn1

# Nitrogen containing

$PROG -t -n AminoMaleic -s AmMale NC\(=O\)C=CC\(=O\)O
$PROG -t -n AminoMalonic -s AmMalon N\(C\(=O\)O\)C\(=O\)O
$PROG -t -n Aziridine -s AzIrdn N1CC1
$PROG -t -n Azo -s NNH N=N
$PROG -t -n DiMethylAmine -s NMe2 N\(C\)C
$PROG -t -n HemiAminal -s HmAmal NC\(O\)C
$PROG -t -n MethylAzo -s NNMe N=NC
$PROG -t -n Nitrate -s NO3 O[N+]\([O-]\)\(=O\)
$PROG -t -n Nitrosyl -s NO N=O

# Oxygen containing

$PROG -t -n Carbonate -s OCOOH OC\(=O\)O
$PROG -t -n Hydroperoxy -s OOH OO
$PROG -t -n MethylCarbonate -s OCOOMe OC\(=O\)OC
$PROG -t -n MethylPeroxy -s OOMe OOC
$PROG -t -n NitroPeroxy -s OONOOH OO[N+]\(=O\)[O-]
$PROG -t -n Nitrosooxy -s ONO ON=O
$PROG -t -n Pro1eneEther -s OPr1e OC=CC
$PROG -t -n PropylEther -s OPr OCCC

# Phosphorous containing

$PROG -t -n DiMethylPhosphanyl -s PMe2 P\(C\)C
$PROG -t -n MethylPhosphanyl -s PHMe PC
$PROG -t -n MethylPhosphate -s MeHPO4 OP\(O\)\(=O\)OC
$PROG -t -n Phosphate -s PO4H2 OP\(O\)\(=O\)O
$PROG -t -n Phosphono -s PO3H2 P\(O\)\(O\)=O

# Sulphur containing

$PROG -t -n DiSulfide -s SSH SS
$PROG -t -n MethylDisulfide -s SSMe SSC
$PROG -t -n MethylSulfide -s SMe SC
$PROG -t -n MethylSulfinyl -s SOMe S\(=O\)C
$PROG -t -n MethylThioyl -s CSMe C\(=S\)C
$PROG -t -n Sulfinamide -s SONH2 S\(=O\)N
$PROG -t -n Sulfino -s SOOH S\(=O\)O
$PROG -t -n Sulfinyl -s SOH S\([H]\)\(=O\)
$PROG -t -n Sulfonamide -s SO2NH2 S\(=O\)\(=O\)N
$PROG -t -n Thial -s CHS C=S
$PROG -t -n Thiazole:1 -s Thzl:1 S1C=NC=C1
$PROG -t -n Thiazole:2 -s Thzl:2 C1=NC=CS1
$PROG -t -n Thiazole:4 -s Thzl:4 C1=CSC=N1
$PROG -t -n Thiazole:5 -s Thzl:5 C1=CN=CS1
$PROG -t -n Thiazolium:1 -s Thzlum:1 S1C=[N+]C=C1
$PROG -t -n Thiazolium:2 -s Thzlum:2 C1=[N+]C=CS1
$PROG -t -n Thiazolium:3 -s Thzlum:3 [N+]1=CSC=C1
$PROG -t -n Thiazolium:4 -s Thzlum:4 C1=CSC=[N+]1
$PROG -t -n Thiazolium:5 -s Thzlum:5 C1=C[N+]=CS1
$PROG -t -n Thiol -s SH S\([H]\)

# Mixed groups

$PROG -t -n 2BromoFuran:3 -s 2BrFur:3 c1c\(Br\)occ1
$PROG -t -n 2BromoFuran:4 -s 2BrFur:4 c1cc\(Br\)oc1
$PROG -t -n 2BromoFuran:5 -s 2BrFur:5 c1ccc\(Br\)o1
$PROG -t -n 2BromoThiophene:3 -s 2BrThphn:3 c1c\(Br\)scc1
$PROG -t -n 2BromoThiophene:4 -s 2BrThphn:4 c1cc\(Br\)sc1
$PROG -t -n 2BromoThiophene:5 -s 2BrThphn:5 c1ccc\(Br\)s1
$PROG -t -n 2ChloroFuran:3 -s 2ClFur:3 c1c\(Cl\)occ1
$PROG -t -n 2ChloroFuran:4 -s 2ClFur:4 c1cc\(Cl\)oc1
$PROG -t -n 2ChloroFuran:5 -s 2ClFur:5 c1ccc\(Cl\)o1
$PROG -t -n 2ChloroThiophene:3 -s 2ClThphn:3 c1c\(Cl\)scc1
$PROG -t -n 2ChloroThiophene:4 -s 2ClThphn:4 c1cc\(Cl\)sc1
$PROG -t -n 2ChloroThiophene:5 -s 2ClThphn:5 c1ccc\(Cl\)s1
$PROG -t -n 2EthylFuran:3 -s 2EtFur:3 c1c\(CC\)occ1
$PROG -t -n 2EthylFuran:4 -s 2EtFur:4 c1cc\(CC\)oc1
$PROG -t -n 2EthylFuran:5 -s 2EtFur:5 c1oc\(CC\)cc1
$PROG -t -n 2EthylThiophene:3 -s 2EtThphn:3 c1c\(CC\)scc1
$PROG -t -n 2EthylThiophene:4 -s 2EtThphn:4 c1cc\(CC\)sc1
$PROG -t -n 2EthylThiophene:5 -s 2EtThphn:5 c1sc\(CC\)cc1
$PROG -t -n 2FluoroFuran:3 -s 2FFur:3 c1c\(F\)occ1
$PROG -t -n 2FluoroFuran:4 -s 2FFur:4 c1cc\(F\)oc1
$PROG -t -n 2FluoroFuran:5 -s 2FFur:5 c1ccc\(F\)o1
$PROG -t -n 2FluoroThiophene:3 -s 2FThphn:3 c1c\(F\)scc1
$PROG -t -n 2FluoroThiophene:4 -s 2FThphn:4 c1cc\(F\)sc1
$PROG -t -n 2FluoroThiophene:5 -s 2FThphn:5 c1ccc\(F\)s1
$PROG -t -n 2MethylFuran:3 -s 2MeFur:3 c1c\(C\)occ1
$PROG -t -n 2MethylFuran:4 -s 2MeFur:4 c1cc\(C\)oc1
$PROG -t -n 2MethylFuran:5 -s 2MeFur:5 c1oc\(C\)cc1
$PROG -t -n 2MethylThiophene:3 -s 2MeThphn:3 c1c\(C\)scc1
$PROG -t -n 2MethylThiophene:4 -s 2MeThphn:4 c1cc\(C\)sc1
$PROG -t -n 2MethylThiophene:5 -s 2MeThphn:5 c1sc\(C\)cc1
$PROG -t -n 3BromoFuran:2 -s 3BrFur:2 c1occc1\(Br\)
$PROG -t -n 3BromoFuran:4 -s 3BrFur:4 c1c\(Br\)coc1
$PROG -t -n 3BromoFuran:5 -s 3BrFur:5 c1cc\(Br\)co1
$PROG -t -n 3BromoThiophene:2 -s 3BrThphn:2 c1sccc1\(Br\)
$PROG -t -n 3BromoThiophene:4 -s 3BrThphn:4 c1c\(Br\)csc1
$PROG -t -n 3BromoThiophene:5 -s 3BrThphn:5 c1cc\(Br\)cs1
$PROG -t -n 3ChloroFuran:2 -s 3ClFur:2 c1occc1\(Cl\)
$PROG -t -n 3ChloroFuran:4 -s 3ClFur:4 c1c\(Cl\)coc1
$PROG -t -n 3ChloroFuran:5 -s 3ClFur:5 c1cc\(Cl\)co1
$PROG -t -n 3ChloroThiophene:2 -s 3ClThphn:2 c1sccc1\(Cl\)
$PROG -t -n 3ChloroThiophene:4 -s 3ClThphn:4 c1c\(Cl\)csc1
$PROG -t -n 3ChloroThiophene:5 -s 3ClThphn:5 c1cc\(Cl\)cs1
$PROG -t -n 3FluoroFuran:2 -s 3FFur:2 c1occc1\(F\)
$PROG -t -n 3FluoroFuran:4 -s 3FFur:4 c1c\(F\)coc1
$PROG -t -n 3FluoroFuran:5 -s 3FFur:5 c1cc\(F\)co1
$PROG -t -n 3FluoroThiophene:2 -s 3FThphn:2 c1sccc1\(F\)
$PROG -t -n 3FluoroThiophene:4 -s 3FThphn:4 c1c\(F\)csc1
$PROG -t -n 3FluoroThiophene:5 -s 3FThphn:5 c1cc\(F\)cs1
$PROG -t -n AcetAldoxime -s AcAdxm C\(C\)=NO
$PROG -t -n Aldoxime -s Adxm C=NO
$PROG -t -n AminoAldoxime -s AmAdxm C\(N\)=NO
$PROG -t -n FluoroAldoxime -s FAdxm C\(F\)=NO
$PROG -t -n ChloroAldoxime -s ClAdxm C\(Cl\)=NO
$PROG -t -n BromoAldoxime -s BrAdxm C\(Br\)=NO
$PROG -t -n PropionAldoxime -s PrAdxm C\(CC\)=NO
$PROG -t -n BiFuran:3 -s BiFur:3 c1c\(c2ccco2\)occ1
$PROG -t -n BiFuran:4 -s BiFur:4 c1coc\(c2ccco2\)c1
$PROG -t -n BiFuran:5 -s BiFur:5 c1ccc\(c2ccco2\)o1
$PROG -t -n BiThiophene:3 -s BiThphn:3 c1c\(c2cccs2\)scc1
$PROG -t -n BiThiophene:4 -s BiThphn:4 c1csc\(c2cccs2\)c1
$PROG -t -n BiThiophene:5 -s BiThphn:5 c1ccc\(c2cccs2\)s1
$PROG -t -n BromoBiFuran:3 -s BrBiFur:3 c1c\(c2ccc\(Br\)o2\)occ1
$PROG -t -n BromoBiFuran:4 -s BrBiFur:4 c1cc\(c2ccc\(Br\)o2\)oc1
$PROG -t -n BromoBiFuran:5 -s BrBiFur:5 c1ccc\(c2ccc\(Br\)o2\)o1
$PROG -t -n BromoBiThiophene:3 -s BrBiThphn:3 c1c\(c2ccc\(Br\)s2\)scc1
$PROG -t -n BromoBiThiophene:4 -s BrBiThphn:4 c1cc\(c2ccc\(Br\)s2\)sc1
$PROG -t -n BromoBiThiophene:5 -s BrBiThphn:5 c1ccc\(c2ccc\(Br\)s2\)s1
$PROG -t -n Furan:2 -s Fur:2 C1=CC=CO1
$PROG -t -n Furan:3 -s Fur:3 C1=COC=C1
$PROG -t -n FuranAmine:3 -s FurAm:3 c1c\(N\)occ1
$PROG -t -n FuranAmine:4 -s FurAm:4 c1coc\(N\)c1
$PROG -t -n FuranAmine:5 -s FurAm:5 C1=CC=C\(N\)O1
$PROG -t -n Pyrazine -s Pyraz c1nccnc1
$PROG -t -n SulfenicAcid -s Senic SO
$PROG -t -n Thiophene:2 -s Thphn:2 C1=CC=CS1
$PROG -t -n Thiophene:3 -s Thphn:2 C1=CSC=C1
$PROG -t -n Thiopheneamine:3 -s ThphnAm:3 c1c\(N\)scc1
$PROG -t -n Thiopheneamine:4 -s ThphnAm:4 c1csc\(N\)c1
$PROG -t -n Thiopheneamine:5 -s ThphnAm:5 C1=CC=C\(N\)S1

# Most common from Ertl no included before

$PROG -t -n 1_3_DiChloroBenzene:4 -s Ph13Cl2:4 c1c\(Cl\)cc\(Cl\)cc1
$PROG -t -n 1Butane -s 1Btn CCCC
$PROG -t -n 2MethylPropyl -s 2MePr CC\(C\)C
$PROG -t -n 3Prop1ene -s 3Pr1e CC=C
$PROG -t -n Acetyl -s OAc OC\(=O\)C
$PROG -t -n Amide:N -s CONH2:N NC\(=O\)C
$PROG -t -n AminoPhenyl:N -s AmPh:N Nc1ccccc1
$PROG -t -n BenzylEster -s BnzlEst OCc1ccccc1
$PROG -t -n Benzyl -s Bnzl Cc1ccccc1
$PROG -t -n BromoBenzene:4 -s PhBr:4 c1ccc\(Br\)cc1
$PROG -t -n ChloroBenzene:2 -s PhCl:2 c1c\(Cl\)cccc1
$PROG -t -n ChloroBenzene:3 -s PhCl:3 c1cc\(Cl\)ccc1
$PROG -t -n ChloroBenzene:4 -s PhCl:4 c1ccc\(Cl\)cc1
$PROG -t -n CycloHexane -s CycHex C1CCCCC1
$PROG -t -n EthylEster -s EtEst C\(=O\)OCC
$PROG -t -n FluoroBenzene:4 -s PhF:4 c1ccc\(F\)cc1
$PROG -t -n HydroxylMethyl -s COH CO
$PROG -t -n MethOxyBenzene:4 -s PhOMe:4 c1ccc\(OC\)cc1
$PROG -t -n MethylBenzene:2 -s PhMe:2 c1c\(C\)cccc1
$PROG -t -n MethylBenzene:4 -s PhMe:4 c1ccc\(C\)cc1
$PROG -t -n Morpholine -s Morph:N N1CCOCC1
$PROG -t -n NitroBenzene -s NO2Ph:4 c1ccc\([N+]\(=O\)[O-]\)cc1
$PROG -t -n PhenylKetone -s COPh C\(=O\)c1ccccc1
$PROG -t -n tert-Butyl -s tBut C\(C\)\(C\)\(C\)
$PROG -t -n 1,5-Anhydroglucitol:6 -s 15AG:6 C1O[C@@H]\(CO\)[C@H]\(O\)[C@@H]\(O\)[C@@H]1\(O\)
$PROG -t -n 1-PropylGuanidine -s PrNCNN CCCNC\(=N\)N
$PROG -t -n 4MethylPhenol -s CPh4OH Cc1ccc\(O\)cc1
$PROG -t -n AceticAcid -s AcOH CC\(=O\)O
$PROG -t -n Butane:3 -s Btn:3 C\(C\)CC
$PROG -t -n DiEthylAmine:N -s NEt2:N N\(CC\)CC
$PROG -t -n Ethanol:1 -s EtOH:1 C\(O\)C
$PROG -t -n FormAmidine -s Amidn C\(N\)=N
$PROG -t -n Guanidine -s NCNN NC\(=N\)N
$PROG -t -n Indole:3 -s Indl:3 c1c\(C=CC=C2\)c2nc1
$PROG -t -n MethylAmide:2 -s MeCONH2 CC\(=O\)N
$PROG -t -n n-ButylAmine:4 -s BuNH2:4 CCCCN
$PROG -t -n PropanoicAcid:3 -s PrOOH:3 CCC\(=O\)O


