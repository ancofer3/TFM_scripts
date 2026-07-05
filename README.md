# TFM_scripts

Este repositorio contiene los scripts y cuadernos de código utilizados para el desarrollo de mi Trabajo de Fin de Máster (TFM). El proyecto se centra en el análisis de datos y el desarrollo de modelos de aprendizaje automático (Machine Learning) para [añadir breve descripción de la finalidad predictiva, por ejemplo: predecir la expresión génica, identificar variantes estructurales o clasificar muestras clínicas].

## 📂 Estructura del Repositorio

El flujo de trabajo del proyecto está organizado en las siguientes carpetas principales:

* **`generando_megatabla/scripts/`**: Scripts encargados del preprocesamiento y la integración de datos. Aquí se realiza la limpieza, filtrado y ensamblaje de los conjuntos de datos biológicos crudos para generar una "megatabla" consolidada que servirá como matriz de características para los modelos.
* **`entrenar_modelo/`**: Contiene el código fuente diseñado para la configuración, entrenamiento y validación de los modelos predictivos a partir de los datos previamente procesados.
* **`notebooks/`**: Cuadernos de Jupyter (Jupyter Notebooks) utilizados para el análisis exploratorio de datos (EDA), la evaluación de métricas de rendimiento y la visualización interactiva de los resultados.

## 🛠️ Requisitos Previos e Instalación

Para ejecutar los scripts y notebooks de este proyecto, es necesario tener un entorno con Python 3.x y las librerías de análisis de datos y bioinformática correspondientes. 

Se recomienda utilizar un entorno virtual (como `conda` o `venv`):

```bash
# Ejemplo para clonar el repositorio
git clone [https://github.com/ancofer3/TFM_scripts.git](https://github.com/ancofer3/TFM_scripts.git)
cd TFM_scripts

# Instalar dependencias (si decides añadir un requirements.txt)
pip install -r requirements.txt
