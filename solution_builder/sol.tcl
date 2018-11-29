set sel [atomselect top all]
set selo [atomselect top {name O}]
set selh [atomselect top {name H}]
set selna [atomselect top {name Na}]
set selcl [atomselect top {name Cl}]
$selh set mass 1.0079401
$selh set charge .42379999
$selo set mass 15.999400
$selo set charge -0.84759998
$selna set mass 22.989999
$selna set charge 1
$selcl set charge -1
$selcl set mass 35.4500008
set minmax [measure minmax $sel -withradii]
set box [vecscale 1.1 [vecsub [lindex $minmax 1] [lindex $minmax 0]]]
pbc set $box
pbc box -center origin
pbc set {19.79388 19.79388 19.79388}
pbc box -center origin
topo guessbonds
# for {set i 0} {$i < 156} {incr i} {
# 	topo addbond [expr 3*$i] [expr 3*$i+1]
# 	topo addbond [expr 3*$i] [expr 3*$i+2]
# }

topo guessangles





set sel [atomselect top all]
pbc set {19.79388 19.79388 19.79388}
pbc warp
set gec [measure center $sel weight mass]
$sel moveby [vecscale -1.0 $gec]
pbc box -center origin
# topo writelammpsdata solsystemtcl.data full