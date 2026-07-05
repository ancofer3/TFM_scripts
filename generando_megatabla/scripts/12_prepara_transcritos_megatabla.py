'''
------------------------------------------------------------------------
Prepara transcritos

Fichero: 01_prepara_transcritos.py

Propósito del programa:
La idea final es hacer una megatabla con muchas columnas con info 
asociada a todos los transcritos de referencia que hay para cada 
muestra. Podemos empezar sacando una tabla en la que pongamos todos los 
transcritos de referencia de cada muestra (sin repetidos) con sus
counts, su gen asociado, diff to TSS (y su CV), diff to TTS (y su CV),
su long media de sus reads y añadiendo dos columnas con el TS y el RIN 
de la muestra a la que pertenecen. Luego, iremos añadiendo columnas a 
la tabla poco a poco. Por razones de control posterior, me he complicado
un poco la vida para agrupar las reads, en vez de por transcrito, por 
transcrito y numero de exones que cubre la read

Versión 3
Autores: Andrés Colomer Fernández
Fecha de última modificación: 20/06/2026

------------------------------------------------------------------------
'''
import pandas as pd
import glob as glob
import matplotlib.pyplot as plt
from diptest import diptest
import os

def cv(x):
    if x.mean() == 0 or x.std() == 0:
        return 0
    else:
        return x.std() / x.mean()
    
def covered_exons(row):
    diff = row["diff_to_gene_TSS"]
    if row["ref_length"] <= abs(diff):
        return 0
    if row["strand"] == "+":
        lengths = list(row["length_exones"])
    else:
        lengths = list(row["length_exones"])[::-1]
    n_exons = row["nExons"]
    acum = 0
    if diff >= 0:
        return n_exons
    else:
        for i, v in enumerate(lengths, start=1):
            acum += v
            # buffer de 50pb
            if (acum) + 50 >= abs(diff):
                return n_exons-i+1
        return 1
            
class Muestra_input:
    def __init__(self, path:str, cols:list, ref):
        # Abriremos el fichero 
        self.tabla = pd.read_csv(path, usecols = cols, sep = "\t")
        # Nos quitamos todos los novel
        self.tabla = self.tabla[self.tabla["associated_transcript"] != "novel"]
        self.tabla = self.tabla.merge(how="left", right=ref, left_on="associated_transcript", right_on="transcript")
        self.tabla["covered_exons"] = self.tabla.apply(covered_exons,axis=1)
        print(self.tabla.isnull().sum())
    def escribe_fich(self,path: str):
        # Guardamos la tabla en un fichero tsv
        self.tabla.to_csv(path, sep='\t', encoding='utf-8', index=False)
    def unicos(self, columna):
        # Realmente aqui lo que querriamos es sacar todos los FSMs 
        return self.tabla.drop_duplicates(subset=columna).reset_index(drop=True)
    def sacaMedia(self,group_by,col_objetivo):
        return self.tabla.groupby(by = group_by)[col_objetivo].mean()
    def sacaCV(self, col):
        cv_values = self.tabla.groupby("associated_transcript")[col].apply(cv)
        return cv_values.reset_index(name=f"CV_{col}")
    def sacaCounts(self,columna):
        counts = self.tabla[self.tabla["structural_category"] == "full-splice_match"].groupby(columna).size()
        return counts.reset_index(name="counts")
    def sacaDipTest(self):
        dip_values = self.tabla.groupby("associated_transcript")["diff_to_gene_TSS"].apply(lambda x: diptest(x)[0])
        return dip_values.reset_index(name="dip_test_TSS")
    def sacaDipTestpval(self):
        dip_pval_values = self.tabla.groupby("associated_transcript")["diff_to_gene_TSS"].apply(lambda x: diptest(x)[1])
        return dip_pval_values.reset_index(name="dip_test_TSS_pval")
    def splitea(self):
        FSMs = self.tabla[self.tabla["structural_category"] == "full-splice_match"]
        ISMs = self.tabla[self.tabla["structural_category"] == "incomplete-splice_match"]
        return FSMs.copy(), ISMs.copy()
    
def SacaRINyTS(path):
    fields = os.path.basename(path).split("_")
    TS = fields[0]
    RIN = fields[2]
    return RIN, TS


def main():
    muestra_tiempo={"TS12_RIN_9.9": 0,
                    "TS10_RIN_9.8": 0,
                    "TS11_RIN_9.7": 0,
                    "TS10_RIN_9.6": 0.5,
                    "TS11_RIN_9.6": 0.5,
                    "TS10_RIN_9.3": 1,
                    "TS11_RIN_9.3": 1,
                    "TS11_RIN_8.9": 3,
                    "TS11_RIN_8.8": 4,
                    "TS11_RIN_8.7": 6,
                    "TS10_RIN_8.4": 3,
                    "TS12_RIN_8.2": 6,
                    "TS10_RIN_7.7": 6,
                    "TS12_RIN_7.3": 8,
                    "TS12_RIN_7.2": 8,
                    }
    muestra_seqbatch={"TS12_RIN_9.9": 2,
                      "TS10_RIN_9.8": 1,
                      "TS11_RIN_9.7": 4,
                      "TS10_RIN_9.6": 1,
                      "TS11_RIN_9.6": 4,
                      "TS10_RIN_9.3": 1,
                      "TS11_RIN_9.3": 2,
                      "TS11_RIN_8.9": 4,
                      "TS11_RIN_8.8": 3,
                      "TS11_RIN_8.7": 3,
                      "TS10_RIN_8.4": 1,
                      "TS12_RIN_8.2": 3,
                      "TS10_RIN_7.7": 1,
                      "TS12_RIN_7.3": 2,
                      "TS12_RIN_7.2": 4,
                    }
    dfs_muestras = []
    # Incluimos información del genePred 
    ref = pd.read_csv("refAnnotation_.genePred", sep="\t", header = None)
    cols = ["transcript", "chrom", "strand", "txStart", "TxEnd", "cdsStart", "cdsEnd","nExons", "Exon_starts","Exon_ends","11","associated_gene","13", "14", "15"]
    ref.columns = cols
    longs_exones = []
    for starts_str, ends_str in zip(ref["Exon_starts"],ref["Exon_ends"]):
        starts = starts_str.split(",")
        ends = ends_str.split(",")
        longs = []
        for s,e in zip(starts,ends):
            if e != "" and s != "":
                longs.append(int(e)-int(s))
        longs_exones.append(longs)
    ref["length_exones"] = longs_exones
    ref = ref[["transcript","chrom","length_exones","strand","nExons"]]
    # Y ahora ya vamos classification por classification
    for i in glob.glob("inputs/*"):
        print(i)
        RIN, TS = SacaRINyTS(i)
        muestra = "_".join(os.path.basename(i).split("_")[0:3])
        seq_batch = muestra_seqbatch[muestra]
        time = muestra_tiempo[muestra]
        tabla = Muestra_input(path=i,
                              cols = ["isoform",
                                    "associated_transcript",
                                    "associated_gene",
                                    "structural_category",
                                    "diff_to_gene_TSS", 
                                    "diff_to_gene_TTS",
                                    "length",
                                    "ref_length",
                                    "ref_exons"],
                              ref=ref)
        # Por un lado nos quedamos con el numero de counts del transcrito:
        tabla.tabla["counts_transcript"] = tabla.tabla.groupby("associated_transcript")["isoform"].transform("count")
        # Y nos quitamos los que no tengan mínimo 5 counts
        tabla.tabla = tabla.tabla[tabla.tabla["counts_transcript"] >= 5]
        # El Hartigan's dip test lo sacamos a nivel de transcrito
        tabla.tabla["dip_test_TSS"] = tabla.tabla.groupby("associated_transcript")["diff_to_gene_TSS"].transform(lambda x: diptest(x)[0])
        tabla.tabla["dip_test_TSS_pval"] = tabla.tabla.groupby("associated_transcript")["diff_to_gene_TSS"].transform(lambda x: diptest(x)[1])

        # Dividimos FSMs y ISMs
        FSMs, ISMs = tabla.splitea()
        print("FSMs NaNs originales:", FSMs["covered_exons"].isnull().sum())
        print("ISMs NaNs originales:", ISMs["covered_exons"].isnull().sum())
        n=1
        dfs_limpios = []
        for df in [FSMs,ISMs]:
            grupos = ["associated_transcript"]
            if n > 1:
                grupos = ["associated_transcript","covered_exons"]
            df["counts"] = df.groupby(grupos)["isoform"].transform("count")
            df = df[df["counts"] >= 4].copy()
            df["CV_diff_to_gene_TSS"] = df.groupby(grupos)["diff_to_gene_TSS"].transform(cv)
            df["CV_diff_to_gene_TTS"] = df.groupby(grupos)["diff_to_gene_TTS"].transform(cv)
            df["diff_to_gene_TSS"] = df.groupby(grupos)["diff_to_gene_TSS"].transform("mean")
            df["diff_to_gene_TTS"] = df.groupby(grupos)["diff_to_gene_TTS"].transform("mean")
            df["length"]= df.groupby(grupos)["length"].transform("mean")
            
            df = df.drop_duplicates(grupos)
            df = df.drop(columns=["transcript","length_exones","nExons","isoform","chrom"])
            dfs_limpios.append(df)
            n += 1
        
        #counts = tabla.sacaCounts("associated_transcript")
        #filas_filt = counts[counts["counts"] >= 4]["associated_transcript"]
        #tabla.tabla = tabla.tabla[tabla.tabla["associated_transcript"].isin(filas_filt)]
        tabla_unicos = pd.concat(dfs_limpios)
        tabla_unicos["TS"] = TS
        tabla_unicos["RIN"] = RIN
        tabla_unicos["seqbatch"] = seq_batch
        tabla_unicos["time"] = time
        dfs_muestras.append(tabla_unicos)
        
    df_final = pd.concat(dfs_muestras)
    print("Valores nulos en df_final:")
    print(df_final.isnull().sum())
    print("Long df:", len(df_final))
    df_final = df_final.dropna()
    print("Long df sin NAs:", len(df_final))
    df_final.to_csv("megatablas/megatabla_1.tsv", sep="\t", encoding="utf-8", index=False)
    
if __name__ == "__main__":
    main()



