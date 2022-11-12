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

enzimas = open('link_bionet.txt')#Hay que eliminar las 10 primeras lineas
lineasenzimas = enzimas.readlines() #Lineas del documento link_bionet.txt

patron_enzima = r'[A-Z]([A-Za-z]|\d){1,}' #Expresión regular que describe las enzimas
er_enzima = re.compile(patron_enzima) #Compilamos esta expresión regular
patron_diana = r'([ATCGRYMKSWBDHVN^]{4,20})(?=$)' #Expresión regular que describe la diana
er_diana = re.compile(patron_diana) #La compilamos
patron_sustituye = r'[RYMKSWBDHVN]' #Este es el patrón de todos los caracteres que tenemos que sustituir
er_sustituye = re.compile(patron_sustituye) #Lo compilamos
patron_corte = r'\^' #El patrón con el carácter ^
er_corte = re.compile(patron_corte) #Lo compilamos

patron_nombreCadenaADN = r'(?<=\>)(\.|\p{L}|\d)+(?= +)'
er_nombreCadenaADN = re.compile(patron_nombreCadenaADN)
nombreCadenaADN = ''
patron_cadenaADN = r'(?<=^)(([ATCG]+ )+)?[ATCG]+'
er_cadenaADN = re.compile(patron_cadenaADN)
cadenaADN = ''
patron_juntarCadena = r'( |\n|\r|\t)'
er_juntarCadena = re.compile(patron_juntarCadena)

Dic_enzimas = {} #Iniciamos el diccionario donde vamos a almacenar las enzimas y sus dianas
Dic_ADN = {} #Iniciamos el diccionario donde vamos a almacenar las cadenas de ADN

enzima = "" #Inicializamos la variable enzima
diana = "" #Inicializamos la variable diana
enzima_anterior = " " #Aqui iniciamos la variable enzima_anterior
posicion = ""

for linea in lineasenzimas: #Recorremos todas las lineas del documento que contiene las enzimas

    enzima_res = er_enzima.search(linea) #Este es un booleano que nos dice si hay una enzima en la linea que estamos recorriendo
    if enzima_res: #Si hay una enzima en la linea que estamos recorriendo
        enzima = linea[enzima_res.start():enzima_res.end()] #Entonces enzima es igual a lo que cumple el patrón de las enzimas en la línea que hemos recorrido

    diana_res = er_diana.search(linea)
    if diana_res:#Si hay una diana en la linea que estamos recorriendo
        diana = linea[diana_res.start():diana_res.end()] #La diana será el trozo de texto que encaja con la expresión regular de las dianas en esa línea
        posicion = quitacorte(diana) #Aqui invocamos la función quitacorte que nos dice la posición en la que se halla el carácter ^
        diana = er_sustituye.sub(sust,diana) #Aqui invocamos la función sust, para sustituir todos los caracteres que tenemos que sustituir
        diana = er_corte.sub('',diana) #Ahora le quitamos el ^ a las dianas

    if (enzima == enzima_anterior): #Comprobamos si la enzima es la misma que la anterior
        diana = "(" + str(Dic_enzimas.get(enzima_anterior)[0]) + "|" + diana + ")" #Modificamos la entrada de la diana concatenando lo necesario (hay que hacer casting de la entrada del diccionario)
    Dic_enzimas[enzima] = [diana,posicion] #Añadimos la entrada de esta línea al diccionario, la clave es la encima y el contenido es la diana
    enzima_anterior=enzima #Aqui guardamos la enzima para que se compruebe si en la siguiente línea vuelve a aparecer la misma
#Este for es el que crea el diccionario de las enzimas con sus dianas
#TENGO QUE AÑADIR EL RECORRIDO PARA CONVERTIR EN EXPRESIONES REGULARES COMPILADAS LAS CADENAS

#for k in Dic_enzimas.keys(): #Esto sirve símplemente para recorrer el diccionario de las enzimas
#    print('%s tiene valor %s' % (k, Dic_enzimas[k]))

for linea in lineasgenes:
    nombreCadenaADN_res = er_nombreCadenaADN.search(linea)
    cadenaADN_res = er_cadenaADN.search(linea)
    if nombreCadenaADN_res:#if la linea contiene un nombre de cadena
        cadenaADN = er_juntarCadena.sub('', cadenaADN)  # entonces guardamos la cadena de ADN, la formatamo
        # print(cadenaADN)
        Dic_ADN[nombreCadenaADN] = cadenaADN  # Y metemos en el diccionario la entrada
        cadenaADN = ' '  # reseteamos la variable de cadena de adn a una cadena vacía
        nombreCadenaADN = linea[nombreCadenaADN_res.start():nombreCadenaADN_res.end()]#entonces guardamos el nombre de la cadena
    elif cadenaADN_res:#elif si la linea contiene cadena de adn
        #print(linea[cadenaADN_res.start():cadenaADN_res.end()])
        cadenaADN = cadenaADN+linea[cadenaADN_res.start():cadenaADN_res.end()]#la concatenamos a la variable cadena de adn existente

for k in Dic_ADN.keys(): #Esto sirve símplemente para recorrer el diccionario de las enzimas
    print('%s tiene valor %s' % (k, Dic_ADN[k]))


#Proyecto final asignatura  de Autómatas y Lenguajes Formales, bioinformática

#Primero importamos los archivos necesarios de las bases de datos
#A continuación procesamos la base de datos para generar un diccionario de enzimas cuyas entradas son sus
#dianas de reconocimiento
#Para almacenar las dianas de reconocimiento antes tenemos que tratarlas

#Después tendremos que hacer una función que permita seleccionar una serie de dianas de reconocimiento
#y poder aplicarlas a una secuencia de ADN dada, devolviendo el mapa de dianas

