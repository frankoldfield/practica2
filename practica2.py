import regex as re
import urllib.request
import sys
def quitacorte(diana):
    patron_corte = r'\^'
    er_corte = re.compile(patron_corte)
    corte_res = er_corte.search(diana)
    if corte_res:
        corte_pos = corte_res.start()
        return corte_pos
    else:
        return 0

def sust(m):#Esta función sirve para poder sustituir de golpe todas las coincidencias
    listabus = ['R', 'Y', 'M', 'K', 'S', 'W', 'B', 'D', 'H', 'V', 'N'] #Esta es la lista de caracteres a sustituir
    listasus = ['[AG]', '[CT]', '[AC]', '[GT]', '[CG]', '[AT]', '[CGT]', '[AGT]', '[ACT]', '[ACG]', '[ACGT]'] #Y esta es la lista con los caracteres que sustituyen a los otros
    listaret = '' #A esta variable vamos a ir adjuntando los caracteres que vamos a sustituir
    for i in m.group(0): #Recorremos todas las coincidencias
        listaret = listaret + listasus[listabus.index(i)]
    return listaret #Devolvemos la cadena completada



genes = open('ALL_C_genes_DNA.txt') #Comienza con una linea en blanco y entre cada bloque de texto hay dos líneas en blanco. En esta variable almacenamos el documento all_c_genes_DNA.txt
lineasgenes = genes.readlines() #Aqui guardamos las lineas del documento
for linea in lineasgenes: #Vamos leyendo todas las lineas del documento
    print(linea)
enzimas = open('link_bionet.txt')#Hay que eliminar las 10 primeras lineas
lineasenzimas = enzimas.readlines() #Lineas del documento link_bionet.txt
patron_enzima = r'[A-Z]([A-Za-z]|\d){1,}' #Expresión regular que describe las enzimas
er_enzima = re.compile(patron_enzima) #Compilamos esta expresión regular
patron_diana = r'([ATCGRYMKSWBDHVN^]{4,20})(?=$)' #Expresión regular que describe la diana
er_diana = re.compile(patron_diana) #La compilamos
Diccionario = {} #Iniciamos el diccionario donde vamos a almacenar las enzimas y sus dianas
for linea in lineasenzimas: #Recorremos todas las lineas del documento que contiene las enzimas
    enzima='' #Inicializamos la variable enzima
    diana = '' #Inicializamos la variable diana
    enzima_res = er_enzima.search(linea) #Este es un booleano que nos dice si hay una enzima en la linea que estamos recorriendo
    if enzima_res: #Si hay una enzima en la linea que estamos recorriendo
        enzima = linea[enzima_res.start():enzima_res.end()] #Entonces enzima es igual a lo que cumple el patrón de las enzimas en la línea que hemos recorrido
    diana_res = er_diana.search(linea)
    if diana_res:#Si hay una diana en la linea que estamos recorriendo
        diana = linea[diana_res.start():diana_res.end()] #La diana será el trozo de texto que encaja con la expresión regular de las dianas en esa línea
        patron_sustituye = r'[RYMKSWBDHVN]' #Este es el patrón de todos los caracteres que tenemos que sustituir
        er_sustituye = re.compile(patron_sustituye) #Lo compilamos
        diana = er_sustituye.sub(sust,diana) #Aqui invocamos la función sust, para sustituir todos los caracteres que tenemos que sustituir
        #print(diana)
        posicion = quitacorte(diana) #Aqui invocamos la función quitacorte que nos dice la posición en la que se halla el carácter ^
        patron_corte = r'\^' #El patrón con el carácter ^
        er_corte = re.compile(patron_corte) #Lo compilamos
        diana = er_corte.sub('',diana) #Ahora le quitamos el ^ a las dianas
        #print(diana)
        #print(posicion)
    Diccionario[enzima] = diana #Añadimos la entrada de esta línea al diccionario, la clave es la encima y el contenido es la diana
for k in Diccionario.keys(): #Esto sirve símplemente para recorrer el diccionario
    print('%s tiene valor %s' % (k,Diccionario[k]))
#Proyecto final asignatura  de Autómatas y Lenguajes Formales, bioinformática

#Primero importamos los archivos necesarios de las bases de datos
#A continuación procesamos la base de datos para generar un diccionario de enzimas cuyas entradas son sus
#dianas de reconocimiento
#Para almacenar las dianas de reconocimiento antes tenemos que tratarlas

#Después tendremos que hacer una función que permita seleccionar una serie de dianas de reconocimiento
#y poder aplicarlas a una secuencia de ADN dada, devolviendo el mapa de dianas

