#!/bin/bash
set -e

MATRIX_BIN="data/bin/sherman5.bin"
OUT_DIR="results"

# Asegurar que el archivo binario de sherman5 exista
if [ ! -f "$MATRIX_BIN" ]; then
    echo "Convirtiendo coleccion de matrices a $MATRIX_BIN..."
    python3 data/convert_mat_to_bin.py
fi

mkdir -p "$OUT_DIR"
rm -f "$OUT_DIR"/history_sherman5_*.csv

echo "=========================================================================================="
echo " BENCHMARK SHERMAN5 (N = 3312, nnz = 20793): Solvers DLK vs Baselines Nativos"
echo "=========================================================================================="
printf "%-22s | %-12s | %-14s | %-10s\n" "Configuración Solver" "Iteraciones" "Residuo Final" "Tiempo (s)"
echo "------------------------------------------------------------------------------------------"

run_case() {
    NAME=$1
    shift
    ARGS="$@"
    HIST_FILE="$OUT_DIR/history_sherman5_${NAME}.csv"
    OUTPUT=$(./test_dlk -matrix "$MATRIX_BIN" $ARGS -pc_type ilu -tol 1e-9 -history_file "$HIST_FILE")
    
    ITS=$(echo "$OUTPUT" | grep "Iteraciones" | awk -F':' '{print $2}' | tr -d ' ')
    RES=$(echo "$OUTPUT" | grep "Residuo" | awk -F':' '{print $2}' | tr -d ' ')
    TIME=$(echo "$OUTPUT" | grep "Tiempo" | awk -F':' '{print $2}' | awk '{print $1}')
    
    printf "%-22s | %-12s | %-14s | %-10s\n" "$NAME" "$ITS" "$RES" "$TIME"
}

run_case "DLK_LQR_LGMRES"  -ksp_type dlk_lgmres -ksp_dlk_m0 30 -ksp_dlk_controller lqr -ksp_dlk_aug 2
run_case "DLK_LQR_GMRES"   -ksp_type dlk -ksp_dlk_m0 30 -ksp_dlk_controller lqr
run_case "DLK_PD_GMRES"    -ksp_type dlk -ksp_dlk_m0 30 -ksp_dlk_controller pd
run_case "LGMRES_m30"      -ksp_type lgmres -ksp_lgmres_restart 30 -ksp_lgmres_augment 2
run_case "GMRES_m30"       -ksp_type gmres -ksp_gmres_restart 30
run_case "BiCGStab"        -ksp_type bcgs

echo "------------------------------------------------------------------------------------------"
echo ""
echo "Generando gráfica de convergencia para sherman5..."
python3 plot_convergence.py "$OUT_DIR"/history_sherman5_*.csv -o "$OUT_DIR/sherman5_convergence.png" -t "PETSc: Convergencia en Matriz Sherman5 (N = 3312)"

echo ""
echo "=========================================================================================="
echo " [✓] Benchmark de sherman5 finalizado exitosamente."
echo " [✓] Gráfica guardada en: $OUT_DIR/sherman5_convergence.png"
echo "=========================================================================================="
