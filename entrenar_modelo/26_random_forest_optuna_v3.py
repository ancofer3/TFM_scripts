import pandas as pd
import polars as pl
import optuna
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn.model_selection import cross_val_score
from pathlib import Path
import joblib

mt = pd.read_csv("megatabla_8.tsv", sep = "\t")

mt["RIN_TS"] = mt["RIN"].astype(str) + mt["TS"].astype(str)
rin_ts_combs = mt[["RIN_TS"]].drop_duplicates()
n_combs = mt.drop_duplicates(["associated_transcript", "RIN_TS"]).groupby("associated_transcript").size()
mt["n_combs"] = mt["associated_transcript"].map(n_combs)
n_combs_max = mt["n_combs"].max()
mt_filt = mt.copy()
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
cols_codons=[]
for a in ["A","C","T","G"]:
    for b in ["A","C","T","G"]:
        for c in ["A","C","T","G"]:
            cols_codons.append(f"codon_{a}{b}{c}")
cols_estaticas.extend(cols_codons)
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


X_train = df_final[(df_final["RIN"] != 9.9) & (df_final["n_combs"] == n_combs_max)].drop(columns=cols_drop)
y_train = df_final[(df_final["RIN"] != 9.9) & (df_final["n_combs"] == n_combs_max)]["porc_ISM"]

# Test 1: Transcritos de RIN 9.9 que estaban en el set de entrenamiento
X_test_vistos = df_final[(df_final["RIN"] == 9.9) & (df_final["n_combs"] == n_combs_max)].drop(columns=cols_drop)
y_test_vistos = df_final[(df_final["RIN"] == 9.9) & (df_final["n_combs"] == n_combs_max)]["porc_ISM"]

# Test 2: Transcritos de RIN 9.9 NUEVOS (no estaban en todas las muestras)
X_test_nuevos = df_final[(df_final["RIN"] == 9.9) & (df_final["n_combs"] < n_combs_max)].drop(columns=cols_drop)
y_test_nuevos = df_final[(df_final["RIN"] == 9.9) & (df_final["n_combs"] < n_combs_max)]["porc_ISM"]
def objective(trial):
    # Definimos donde queremos que busque
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

study = optuna.create_study(direction='maximize') # maximizamos porque neg_MAE es negativo
study.optimize(objective, n_trials=50,n_jobs=32)

# Tras encontrar el mejor modelo:
mejor_modelo = RandomForestRegressor(**study.best_params,random_state=42, n_jobs=-1)
mejor_modelo.fit(X_train, y_train)

predicciones_vistos = mejor_modelo.predict(X_test_vistos)
print("\n--- TEST: TRANSCRITOS QUE YA ESTABAN EN EL TRAINING SET")
print("MAE:", mean_absolute_error(y_test_vistos, predicciones_vistos))
print("RMSE:", mean_squared_error(y_test_vistos, predicciones_vistos))
print("R²:", r2_score(y_test_vistos, predicciones_vistos))

predicciones_nuevos = mejor_modelo.predict(X_test_nuevos)
print("\n--- TEST: TRANSCRITOS QUE NO ESTABAN EN EL TRAINING SET")
print("MAE:", mean_absolute_error(y_test_nuevos, predicciones_nuevos))
print("RMSE:", mean_squared_error(y_test_nuevos, predicciones_nuevos))
print("R²:", r2_score(y_test_nuevos, predicciones_nuevos))

# Guardar el modelo entrenado
joblib.dump(mejor_modelo, 'mejor_modelo_ism_codons.pkl')
# Guardar los datos finales (X_test, y_test, predicciones) en un CSV o parquet
# Unir las predicciones a los dataframes originales para guardarlos
df_resultados_vistos = X_test_vistos.copy()
df_resultados_vistos['porc_ISM_real'] = y_test_vistos
df_resultados_vistos['prediccion'] = predicciones_vistos
df_resultados_vistos.to_csv('resultados_modelo_vistos.csv', index=False)

df_resultados_nuevos = X_test_nuevos.copy()
df_resultados_nuevos['porc_ISM_real'] = y_test_nuevos
df_resultados_nuevos['prediccion'] = predicciones_nuevos
df_resultados_nuevos.to_csv('resultados_modelo_nuevos.csv', index=False)

import matplotlib.pyplot as plt

n=0
for ax, i in zip(axes,[(y_test_nuevos,predicciones_nuevos),(y_test_vistos,predicciones_vistos)]):
    ax.scatter(i[0], i[1], alpha=0.5)
    ax.plot([i[0].min(), i[0].max()], [i[0].min(), i[0].max()], 'r--')
    ax.set_xlabel("Real %ISM (RIN 9.9)")
    ax.set_ylabel("Prediction %ISM")
    ax.set_title("New transcripts" if n==0 else "Known transcripts")
    n+=1
fig.suptitle("Model performance: Prediction vs Reality")
plt.savefig("scatter_performances_codons.png")

import seaborn as sns
importancias = mejor_modelo.feature_importances_
nombres_columnas = X_test_vistos.columns

# 2. Crear un DataFrame y ordenarlo de mayor a menor
df_importancias = pd.DataFrame({
    'Feature': nombres_columnas,
    'Importance': importancias
}).sort_values(by='Importance', ascending=False)

# Mostrar el top en texto por consola
print("--- FEATURE IMPORTANCES ---")
print(df_importancias)

# 3. Plotear un gráfico de barras horizontal (es el formato más legible)
plt.figure(figsize=(10, 6))
sns.set_theme(style="whitegrid")

# Usamos un paleta de colores degradada para que quede profesional
sns.barplot(
    x='Importance', 
    y='Feature', 
    data=df_importancias, 
    color="#fc8d59" # O 'viridis', 'crest', etc.
)
plt.yticks(fontsize=8)
plt.title('Feature Importances for ISM Percentage Prediction', fontsize=14, pad=15)
plt.xlabel('Relative Importance (Mean Decrease in Impurity)', fontsize=12)
plt.ylabel('Transcript / Sample Features', fontsize=12)
plt.tight_layout()
plt.savefig("Feature_importances_codons.png")
