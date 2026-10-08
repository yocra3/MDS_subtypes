import pandas as pd
import re
# Cargar archivo de Excel
cun_data = pd.read_excel('data/datosCUN/Mielodisplasias.xlsx', header = 0)

# Eliminar '_x000D_' del texto antes de dividir
texto_limpio = cun_data['DOCHEMATO'][0].replace('_x000D_', '')

# Dividir una cadena de texto que contiene \n en una lista y eliminar cadenas vacías
ej = [linea for linea in texto_limpio.split('\n') if linea not in ['', '  ']]


for text in cun_data['DOCHEMATO']:
    print("DIAGNÓSTICO PRINCIPAL" in text)

texto_limpio2 = cun_data['DOCHEMATO'][3].replace('_x000D_', '')
print(texto_limpio2.replace('\n', '\r\n'))


## MEDULOGRAMA
def extraer_mielograma(texto: str) -> dict:
    if pd.isna(texto):
        return {}
    texto = texto.replace('_x000D_', '')
    lineas = [l.strip() for l in texto.split('\n') if l.strip()]
    en_seccion = False
    contenido = []
    for linea in lineas:
        if "Mielograma" in linea:
            en_seccion = True
            continue
        if "Morfología" in linea:
            en_seccion = False
        if en_seccion:
            cols = linea.split("\t \t")    
            contenido.append(";".join(cols))
    bloque = ";".join(contenido)
    # separa por comas o punto y coma para facilitar el parseo
    partes = re.split(r'[;,]', bloque)
    patron = re.compile(r'([^\t]+)\t\s*(-?\d+(?:[\.,]\d+)?)\s*%')
    res = {}
    for p in partes:
        m = patron.search(p)
        if m:
            nombre = m.group(1).strip().lower().replace('  ', ' ')
            valor = float(m.group(2).replace(',', '.'))
            res[nombre] = valor
    return res


mielograma_df = cun_data['DOCMEDULOGRAMA'].apply(extraer_mielograma).apply(pd.Series)

### Asumir no blastos es 0
mielograma_df['blastos'] = mielograma_df['blastos'].fillna(0)

## Cariotipo
car2 = cun_data['DOCHCARIOTIPO'][2].replace('_x000D_', '')
print(car2.replace('\n', '\r\n'))

for text in cun_data['DOCHCARIOTIPO']:
    print("Cariotipo" in text, "ISCN" in text)


def extraer_cariotipo(texto: str) -> dict:
    if pd.isna(texto):
        return {}
    texto = texto.replace('_x000D_', '')
    lineas = [l.strip() for l in texto.split('\n') if l.strip()]
    en_seccion = False
    contenido = []
    patron = re.compile(r'^\d')
    patron2 = re.compile(r'([^:]+):\s*(\d[\d\w\s\-\.,%\(\)]*)')
    for linea in lineas:
        if "Cariotipo" in linea or "ISCN" in linea:
            en_seccion = True
            m1 = patron2.search(linea)
            car = m1.group(2) if m1 else ""
            contenido.append(car)
            continue
        if "Morfología" in linea:
            en_seccion = False
        if en_seccion:
            m = patron.search(linea)
            if m:
                contenido.append(linea)
            else:
                en_seccion = False
                break
    return("/".join(contenido))

cariotipos = cun_data['DOCHCARIOTIPO'].apply(extraer_cariotipo)

## Genetica
gen = cun_data['DOCPANELNGS'][5].replace('_x000D_', '')
print(gen.replace('\n', '\r\n'))

for text in cun_data['DOCPANELNGS']:
    if isinstance(text, str):
        print("RESULTADOS CON RELEVANCIA CLINICA" in text)
    else:
        print("No NGS")


ngs_headers = ['Gen', 'CambioNucleotido', 'CambioAA', 'Posicion_hg19', 
               'TipoVariante', 'Categoria', 'VAF', 'Profundidad', 'BasesDatos']
def extraer_ngs(texto: str) -> dict:
    if pd.isna(texto):
        return pd.DataFrame(columns=ngs_headers)
    texto = texto.replace('_x000D_', '')
    lineas = [l.strip() for l in texto.split('\n') if l.strip()]
    en_seccion = False
    filas_tabla = []
    for linea in lineas:
        if linea == "RESULTADOS CON RELEVANCIA CLINICA" or linea == "RESULTADOS":
            en_seccion = True
            continue
        if "INTERPRETACIÓN CLÍNICA" in linea:
            en_seccion = False
            break
        if en_seccion:
            # Detectar si es una línea de tabla (contiene tabulaciones)
            if '\t' in linea:
                columnas = [col.strip() for col in linea.split('\t') if col.strip()]
                 # Si la longitud coincide con encabezados, es una fila de datos
                if len(columnas) == len(ngs_headers):
                    filas_tabla.append(columnas)
    # Crear DataFrame si hay datos
    if filas_tabla:
        df = pd.DataFrame(filas_tabla, columns=ngs_headers)
        return df
    else:
        return pd.DataFrame()

ngs_list = [extraer_ngs(text) for text in cun_data['DOCPANELNGS']]


ngs_vec = [";".join(df['Gen'].tolist()) if not df.empty else "" for df in ngs_list]

## Parse included genes
ngs_headers_method = ['GEN', 'CROMOSOMA', 'TRANSCRITO', 'EXONES']
ngs_headers_method2 = ['GEN', 'TRANSCRITO', 'EXONES']

def extraer_ngs_base(texto: str) -> dict:
    if pd.isna(texto):
        return pd.DataFrame(columns=ngs_headers_method)
    texto = texto.replace('_x000D_', '')
    lineas = [l.strip() for l in texto.split('\n') if l.strip()]
    en_seccion = False
    filas_tabla = []
    for linea in lineas:
        if linea == "METODOLOGÍA":
            en_seccion = True
            continue
        if linea == "LIMITACIONES DEL ESTUDIO" or  linea == "PARÁMETROS DE CALIDAD":
            en_seccion = False
            break
        if en_seccion:
            # Detectar si es una línea de tabla (contiene tabulaciones)
            if '\t' in linea:
                columnas = [col.strip() for col in linea.split('\t') if col.strip()]
                 # Si la longitud coincide con encabezados, es una fila de datos
                if len(columnas) == len(ngs_headers_method) or len(columnas) == len(ngs_headers_method2):
                    if columnas != ngs_headers_method and columnas != ngs_headers_method2:  # Evitar agregar la fila de encabezados
                        filas_tabla.append(columnas)
    # Crear DataFrame si hay datos
    if filas_tabla:
        if len(filas_tabla[0]) == len(ngs_headers_method):
            df = pd.DataFrame(filas_tabla, columns=ngs_headers_method)
        else:
            df = pd.DataFrame(filas_tabla, columns=ngs_headers_method2)
        return df
    else:
        return pd.DataFrame(columns=ngs_headers_method)
    
ngs_method_list = [extraer_ngs_base(text) for text in cun_data['DOCPANELNGS']]


## Preparar df final
cun_data['AGE'] = cun_data['FDIAG'].dt.year - cun_data['FNAC'].dt.year
cun_data['OS_STATUS'] = cun_data['FFALLE'].notna().astype(int)
cun_data['OS_YEARS'] = cun_data.apply(
    lambda row: (row['FFALLE'] - row['FDIAG']).days 
                if pd.notna(row['FFALLE']) 
                else (row['FULTIMOCONTACTO'] - row['FDIAG']).days,
    axis=1
)/365.25

cun_data['CARIOTIPO'] = cariotipos
cun_data['GENES_NGS'] = ngs_vec
cun_data['BM_BLASTS'] = mielograma_df['blastos']

cun_final = cun_data[['NH', 'FDIAG', 'SEXO', 'FNAC', 'AGE', 'CARIOTIPO', 'GENES_NGS',
                      'HB',  'LEUCOCITOS','NEUTOFILOS',  'PLAQUETAS',  'BM_BLASTS',  'MONOCITOS',
                       'OS_STATUS', 'OS_YEARS']]

cun_final.to_csv('data/datosCUN/cun_parseo_final.tsv', index=False, sep='\t')