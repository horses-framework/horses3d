#!/bin/bash
#
#  horses2plt regression test. Every case converts the 3-element Nodes_3 mesh
#  (mesh order 1, Gauss nodes) and compares the .tec file against reference/.
#
#     Nodes_3.hsol    : solution order 1 (same as the mesh)
#     Nodes_3_P2.hsol : solution order 2, obtained with horsesConverter
#                       (meshInterpolation of Nodes_3.hsol onto Nodes_3.hmesh)
#
cd "$(dirname "$0")"
TEST_DIR=$(pwd)
HORSES2PLT=$TEST_DIR/../../../bin/horses2plt
FILES=$TEST_DIR/../horsesConverter/Files
WORK=$TEST_DIR/RESULTS

rm -rf "$WORK"
mkdir -p "$WORK"
cp "$FILES/Nodes_3.hmesh" "$FILES/Nodes_3.hsol" "$TEST_DIR/Nodes_3_P2.hsol" "$WORK"

failed=0

#  compare <reference> <output>
compare() {
   awk '
      NR == FNR { ref[FNR] = $0; nref = FNR; next }
      {
         nout = FNR
         if ($1 ~ /^[A-Z]/) {
            if ($0 != ref[FNR]) { print "   line " FNR ": expected \"" ref[FNR] "\", got \"" $0 "\""; bad = 1 }
            next
         }
         nr = split(ref[FNR], r)
         if (nr != NF) { print "   line " FNR ": expected " nr " columns, got " NF; bad = 1; next }
         for (c = 1; c <= NF; c++) {
            d = $c - r[c]; if (d < 0) d = -d
            a = r[c];      if (a < 0) a = -a
            if (d > 1.0e-8 * a + 1.0e-12) { print "   line " FNR ", column " c ": expected " r[c] ", got " $c; bad = 1 }
         }
      }
      END {
         if (nout != nref) { print "   expected " nref " lines, got " nout; bad = 1 }
         exit bad
      }' "$1" "$2"
}

#  run_case <name> <solution file> [horses2plt flags...]
run_case() {
   local name=$1 sol=$2
   shift 2
   echo "--- $name"
   if ! (cd "$WORK" && "$HORSES2PLT" Nodes_3.hmesh "$sol" "$@" > "$name.log" 2>&1); then
      echo "   FAILED: horses2plt exited with an error (see RESULTS/$name.log)"
      failed=1
      return
   fi
   mv "$WORK/${sol%.hsol}.tec" "$WORK/$name.tec"
   if compare "$TEST_DIR/reference/$name.tec" "$WORK/$name.tec"; then
      echo "   PASSED"
   else
      echo "   FAILED: output differs from reference/$name.tec"
      failed=1
   fi
}

run_case SameOrder              Nodes_3.hsol
run_case MeshOrderDiffers       Nodes_3_P2.hsol
run_case LowerOrderHomogeneous  Nodes_3_P2.hsol --output-order=1
run_case LowerOrderGauss        Nodes_3_P2.hsol --output-order=1 --output-basis=Gauss

exit $failed
