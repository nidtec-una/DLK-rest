#!/bin/bash
set -e

GRID_M=${1:-50}
GRID_N=${2:-50}
OUT_DIR="results"

mkdir -p "$OUT_DIR"
rm -f "$OUT_DIR"/history_*.csv

echo "=========================================================================="
echo " EXPERIMENTOS PETSc: DLK-GMRES vs Nativos (Malla ${GRID_M}x${GRID_N})"
echo "=========================================================================="
printf "%-22s | %-12s | %-14s | %-10s\n" "Configuración Solver" "Iteraciones" "Residuo Final" "Tiempo (s)"
echo "--------------------------------------------------------------------------"

run_case() {
    NAME=$1
    shift
    ARGS="$@"
    HIST_FILE="$OUT_DIR/history_${NAME}.csv"
    OUTPUT=$(./test_dlk -m "$GRID_M" -n "$GRID_N" $ARGS -history_file "$HIST_FILE")
    
    ITS=$(echo "$OUTPUT" | grep "Iteraciones" | awk -F':' '{print $2}' | tr -d ' ')
    RES=$(echo "$OUTPUT" | grep "Residuo" | awk -F':' '{print $2}' | tr -d ' ')
    TIME=$(echo "$OUTPUT" | grep "Tiempo" | awk -F':' '{print $2}' | awk '{print $1}')
    
    printf "%-22s | %-12s | %-14s | %-10s\n" "$NAME" "$ITS" "$RES" "$TIME"
}

run_case "DLK_LQR_m20" -ksp_type dlk -ksp_dlk_m0 20 -ksp_dlk_controller lqr
run_case "DLK_PD_m20"  -ksp_type dlk -ksp_dlk_m0 20 -ksp_dlk_controller pd
run_case "GMRES_m20"   -ksp_type gmres -ksp_gmres_restart 20
run_case "GMRES_m30"   -ksp_type gmres -ksp_gmres_restart 30
run_case "BiCGStab"    -ksp_type bcgs

echo "--------------------------------------------------------------------------"
echo ""
echo "Generando gráfica comparativa..."
python3 plot_convergence.py "$OUT_DIR"/history_*.csv -o "$OUT_DIR/convergence_dlk_vs_gmres.png" -t "PETSc: DLK-GMRES Adaptativo vs Nativos (${GRID_M}x${GRID_N})"

echo ""
echo "=========================================================================="
echo " [✓] Experimento finalizado con éxito."
echo " [✓] Gráfica guardada en: $OUT_DIR/convergence_dlk_vs_gmres.png"
echo "=========================================================================="
