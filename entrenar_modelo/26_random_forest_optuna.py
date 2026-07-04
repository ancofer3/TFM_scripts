import pandas as pd
import polars as pl
import optuna
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn.model_selection import cross_val_score
from pathlib import Path
import joblib

mt = pd.read_csv("megatabla_7.tsv", sep = "\t")
mt["RIN_TS"] = mt["RIN"].astype(str) + mt["TS"].astype(str)
rin_ts_combs = mt[["RIN_TS"]].drop_duplicates()
n_combs = mt.drop_duplicates(["associated_transcript", "RIN_TS"]).groupby("associated_transcript").size()
mt["n_combs"] = mt["associated_transcript"].map(n_combs)
n_combs_max = mt["n_combs"].max()
mt_filt = mt[mt["n_combs"] == n_combs_max]
total_counts = mt_filt.groupby(['associated_transcript', 'RIN_TS'])['counts'].sum().reset_index(name='total_transcript_counts')
ism_counts = mt_filt[mt_filt['structural_category'] == 'incomplete-splice_match'].groupby(['associated_transcript', 'RIN_TS'])['counts'].sum().reset_index(name='ism_counts')

ISMs = pd.merge(total_counts, ism_counts, on=['associated_transcript', 'RIN_TS'], how='left')
ISMs['ism_counts'] = ISMs['ism_counts'].fillna(0) # Si no hay ISMs, es 0
ISMs['porc_ISM'] = (ISMs['ism_counts'] / ISMs['total_transcript_counts']) * 100
cols_estaticas = ['associated_gene','associated_transcript', 'ref_length', 'ref_exons', 
                  'strand','counts_transcript','TS', 'RIN', 
                  'seqbatch', 'time','dip_test_TSS', 'dip_test_TSS_pval',
                  'exon_junction_density', 'porc_GC', 'CPM_transcript',
                  'length_CDS', 'porc_GC_CDS', 'length_3pUTR', 'porc_GC_3pUTR',
                  'length_5pUTR', 'porc_GC_5pUTR', 'RIN_TS', 'n_combs']
mt_estatico = mt_filt[cols_estaticas].drop_duplicates(subset=["associated_transcript","RIN_TS"])
ISMs = pd.merge(ISMs, mt_estatico, on=['associated_transcript','RIN_TS'], how='left')

cols_dinamicas = ['length','diff_to_gene_TSS','diff_to_gene_TTS','CV_diff_to_gene_TSS',
                  'CV_diff_to_gene_TTS','Cov_D1','Cov_D2', 'Cov_D3',
                  'Cov_D4', 'Cov_D5', 'Cov_D6', 'Cov_D7', 
                  'Cov_D8', 'Cov_D9','Cov_D10']
# Para cada columna de la que queremos sacar la media ponderada
for col in cols_dinamicas:
       mt_filt[f'{col}_ponderada'] = mt_filt[col]*mt_filt["counts"]
instrucciones = {f'{col}_ponderada' :'sum' for col in cols_dinamicas}
instrucciones["counts"] = 'sum'
dinamicas_grouped = mt_filt.groupby(["associated_transcript","RIN_TS"]).agg(instrucciones).reset_index()
for col in cols_dinamicas:
    dinamicas_grouped[f'avg_{col}'] = dinamicas_grouped[f'{col}_ponderada'] / dinamicas_grouped['counts']
cols_finales = ['associated_transcript', 'RIN_TS', 'counts'] + [f'avg_{col}' for col in cols_dinamicas]
df_dinamicas_limpio = dinamicas_grouped[cols_finales]
df_dinamicas_limpio.rename(columns={'counts': 'total_counts'}, inplace=True)
df_final = pd.merge(ISMs,df_dinamicas_limpio, on=["associated_transcript","RIN_TS"], how="left")
df_final['strand'] = df_final['strand'].map({'+': 1, '-': 0})
cols_drop = ["associated_transcript", "associated_gene", "RIN_TS","porc_ISM","ism_counts","seqbatch","time","n_combs","total_counts","TS"]

X_test = df_final[df_final["RIN"]==9.9].drop(columns=cols_drop)
X_train= df_final[df_final["RIN"]!=9.9].drop(columns=cols_drop)
y_test = df_final[df_final["RIN"]==9.9]["porc_ISM"]
y_train = df_final[df_final["RIN"]!=9.9]["porc_ISM"]
def objective(trial):
    # Definimos el espacio de búsqueda
    n_estimators = trial.suggest_int('n_estimators', 100, 1000)
    max_depth = trial.suggest_int('max_depth', 5, 50)
    min_samples_split = trial.suggest_int('min_samples_split', 2, 10)
    
    # Creamos el modelo con los parámetros sugeridos por Optuna
    model = RandomForestRegressor(
        n_estimators=n_estimators,
        max_depth=max_depth,
        min_samples_split=min_samples_split,
        random_state=42,
        n_jobs=1
    )
    
    # Usamos cross-validation para evaluar
    score = cross_val_score(model, X_train, y_train, cv=3, scoring='neg_mean_absolute_error').mean()
    return score

# Ejecutamos la optimización
study = optuna.create_study(direction='maximize') # maximize porque neg_MAE es negativo
study.optimize(objective, n_trials=50,n_jobs=32)
fig = optuna.visualization.plot_optimization_history(study)
fig.write_image("optimization_history.png")
fig2=optuna.visualization.plot_param_importances(study).show()
fig2.write_image("param_importances.png")
# Tras encontrar el mejor modelo:
mejor_modelo = RandomForestRegressor(**study.best_params, random_state=42)
pred=mejor_modelo.fit(X_train, y_train)
predicciones = mejor_modelo.predict(X_test)
print("MAE:", mean_absolute_error(y_test, pred))
print("RMSE:", mean_squared_error(y_test, pred, squared=False))
print("R²:", r2_score(y_test, pred))

df_resultados = df_final[df_final["RIN"]==9.9].copy()
df_resultados["prediccion"] = pred
# Guardar el modelo entrenado
joblib.dump(mejor_modelo, 'mejor_modelo_ism.pkl')
# Guardar los datos finales (X_test, y_test, predicciones) en un CSV o parquet
df_resultados.to_csv('resultados_modelo.csv', index=False)

import matplotlib.pyplot as plt

plt.scatter(y_test, predicciones, alpha=0.5)
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--')
plt.xlabel("Real %ISM (RIN 9.9)")
plt.ylabel("Predicción %ISM")
plt.title("Rendimiento del Modelo: Predicción vs Realidad")
plt.savefig("prediccion_vs_real.png")
