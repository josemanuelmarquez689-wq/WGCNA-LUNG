#!/usr/bin/env bash
# Corre todo el analisis en orden y junta las figuras finales en figures/.
# Uso: bash run_all.sh      (tiempo total ~1.5 h)
set -euo pipefail
cd "$(dirname "$0")"
R=/opt/R/current/bin/Rscript
mkdir -p logs figures
echo "[$(date '+%F %T')] 04 figuras de la red";              $R scripts/04_figuras_red.R               > logs/04_figuras_red.log 2>&1
echo "[$(date '+%F %T')] 03 conservacion del modulo";        NUCLEOS=20 $R scripts/03_conservacion_modulo_jaccard.R > logs/03_conservacion.log 2>&1
echo "[$(date '+%F %T')] 01 integridad global";              NUCLEOS=20 N_REP=30 $R scripts/01_robustez_global_integridad.R > logs/01_integridad.log 2>&1
echo "[$(date '+%F %T')] 02 control del mismo modulo";       NUCLEOS=24 N_REP=30 $R scripts/02_control_mismo_modulo.R > logs/02_control_modulo.log 2>&1
cp results/conservacion/Heatmap_Conservacion_Modulo.pdf            figures/Fig4_Conservacion_Modulo_Jaccard.pdf
cp results/integridad/Heatmap_Eficiencia_Global_Relativa.pdf       figures/Fig5_Eficiencia_Global_Relativa.pdf
cp results/integridad/Heatmap_Componente_Gigante_Relativo.pdf      figures/Fig6_Componente_Gigante_Relativo.pdf
cp results/integridad/Heatmap_Dano_vs_Azar_Eficiencia.pdf          figures/Fig7_Dano_vs_Azar_Red.pdf
cp results/integridad/Heatmap_Dano_vs_Azar_MismoModulo.pdf         figures/Fig8_Dano_vs_Azar_MismoModulo.pdf
cp results/integridad/Heatmap_Razon_Dano_Hubs_vs_MismoModulo.pdf   figures/Fig9_Razon_Dano_Hubs_vs_MismoModulo.pdf
echo "[$(date '+%F %T')] LISTO"
