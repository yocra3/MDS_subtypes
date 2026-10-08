import pandas as pd

# Cargar archivo de Excel
cun_data = pd.read_excel('data/Mielodisplasias.xlsx', header = 0)

# Eliminar '_x000D_' del texto antes de dividir
texto_limpio = cun_data['DOCHEMATO'][0].replace('_x000D_', '')

# Dividir una cadena de texto que contiene \n en una lista y eliminar cadenas vacías
ej = [linea for linea in texto_limpio.split('\n') if linea != '']
