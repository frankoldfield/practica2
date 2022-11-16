import regex as re
import urllib.request
import ssl

Dic_ADN = {}  # Iniciamos el diccionario donde vamos a almacenar las cadenas de ADN
Dic_enzimas = {}  # Iniciamos el diccionario donde vamos a almacenar las enzimas y sus dianas
patron_corte = r'\^' #Este es el patrón que utilizaremos para buscar el carácter ^ para eliminarlo de las dianas
er_corte = re.compile(patron_corte) #Lo compilamos aqui fuera para que al iterarse varias veces la función quitacorte no tenga que compilarse un número de veces innecesario

def quitacorte(diana): #Esta función es la que utilizamos para quitar el carácter ^, y nos devuelve la posición en la que estaba
    corte_res = er_corte.search(diana)#Este booleano es cierto si se encuentra el patrón ^ en la cadena diana
    if corte_res: #En caso de que se encuentre entonces:
        corte_pos = corte_res.start() #La posición en la que estaba es el inicio de la cadena
        return corte_pos #La devolvemos como argumento de la función
    else:
        return 0 #Si no, devolvemos la posición 0, como indica el programa

def sust(m):#Esta función sirve para poder sustituir de golpe todas las coincidencias
    listabus = ['R', 'Y', 'M', 'K', 'S', 'W', 'B', 'D', 'H', 'V', 'N'] #Esta es la lista de caracteres a sustituir
    listasus = ['[AG]', '[CT]', '[AC]', '[GT]', '[CG]', '[AT]', '[CGT]', '[AGT]', '[ACT]', '[ACG]', '[ACGT]'] #Y esta es la lista con los caracteres que sustituyen a los otros
    listaret = '' #A esta variable vamos a ir adjuntando los caracteres que vamos a sustituir
    for i in m.group(0): #Recorremos todas las coincidencias
        listaret = listaret + listasus[listabus.index(i)]
    return listaret #Devolvemos la cadena completada

def leerEnzimas():
    patron_enzima = r'[A-Z]([A-Za-z]|\d){1,}'  # Expresión regular que describe las enzimas
    er_enzima = re.compile(patron_enzima)  # Compilamos esta expresión regular
    patron_diana = r'([ATCGRYMKSWBDHVN^]{4,20})(?=$)'  # Expresión regular que describe las dianas
    er_diana = re.compile(patron_diana)  # La compilamos
    patron_sustituye = r'[RYMKSWBDHVN]'  # Este es el patrón de todos los caracteres que tenemos que sustituir
    er_sustituye = re.compile(patron_sustituye)  # Lo compilamos
    enzima = ""  # Inicializamos la variable enzima
    diana = ""  # Inicializamos la variable diana
    enzima_anterior = " "  # Aqui inicializamos la variable enzima_anterior
    posicion = 0
    print("Cargando bionet...")
    try:
        url = 'https://aulavirtual.um.es/access/content/group/1896_G_2022_N_N/PRACTICAS/PRACTICA%202/link_bionet.txt'
        context = ssl.create_default_context()
        context.set_ciphers("DEFAULT")
        link = urllib.request.urlopen(url, context=context)
        contador = 0
        for linea in link:
            linea = linea.decode().strip()
            if contador > 9:
                enzima_res = er_enzima.search(linea)  # Este es un booleano que nos dice si hay una enzima en la linea que estamos recorriendo
                if enzima_res:  # Si hay una enzima en la linea que estamos recorriendo
                    enzima = linea[enzima_res.start():enzima_res.end()]  # Entonces enzima es igual a lo que cumple el patrón de las enzimas en la línea que hemos recorrido
                diana_res = er_diana.search(linea)
                if diana_res:  # Si hay una diana en la linea que estamos recorriendo
                    diana = linea[diana_res.start():diana_res.end()]  # La diana será el trozo de texto que encaja con la expresión regular de las dianas en esa línea
                    posicion = quitacorte(diana)  # Aqui invocamos la función quitacorte que nos dice la posición en la que se halla el carácter ^
                    diana = er_sustituye.sub(sust,diana)  # Aqui invocamos la función sust, para sustituir todos los caracteres que tenemos que sustituir
                    diana = er_corte.sub('', diana)  # Ahora le quitamos el ^ a las dianas
                if (enzima == enzima_anterior):  # Comprobamos si la enzima es la misma que la anterior
                    diana = str(Dic_enzimas.get(enzima_anterior)[0]) + "|" + diana  # Modificamos la entrada de la diana concatenando lo necesario (hay que hacer casting de la entrada del diccionario)
                Dic_enzimas[enzima] = [diana,posicion]  # Añadimos la entrada de esta línea al diccionario, la clave es la encima y el contenido es la diana
                enzima_anterior = enzima  # Aqui guardamos la enzima para que se compruebe si en la siguiente línea vuelve a aparecer la misma
            else:
                contador += 1
    except IOError as e:
        print('link_bionet no disponible:', e)

    #for linea in lineasenzimas:  # Recorremos todas las lineas del documento que contiene las enzimas


    for k in Dic_enzimas.keys():  # Este bucle sirve para recorrer el diccionario de las enzimas
        Dic_enzimas[k][0] = re.compile("(" + Dic_enzimas.get(k)[0] + ")") #Compilamos todas las er de las enzimas del diccionario de golpe
        print('%s tiene valor %s' % (k, Dic_enzimas[k]))

def leerGenes():
    print("Cargando All_C_genes_DNA.txt...")
    genes = open('ALL_C_genes_DNA.txt')  # Comienza con una linea en blanco y entre cada bloque de texto hay dos líneas en blanco. En esta variable almacenamos el documento all_c_genes_DNA.txt
    lineasgenes = genes.readlines()  # Aqui guardamos las lineas del documento
    patron_nombreCadenaADN = r'(?<=\>)(\.|\p{L}|\d)+(?= +)' #Este patrón lo utilizaremos para buscar los nombres de las cadenas de Adn
    er_nombreCadenaADN = re.compile(patron_nombreCadenaADN) #Lo compilamos
    nombreCadenaADN = '' #Inicializamos la variable que va a contener en cada iteración el nombre de la cadena de ADN
    patron_numNucleotidos = r'\d+(?= nt)' #Este es el patrón que utilizaremos para buscar el número de nucleótidos
    er_numNucleotidos = re.compile(patron_numNucleotidos) #La expresión regular ya compilada
    numNucleotidos = 0 #Inicializamos la variable del número de nucleótidos
    patron_cadenaADN = r'(?<=^)(([ATCG]+ )+)?[ATCG]+' #Expresión regular que va a buscar las cadenas de ADN
    er_cadenaADN = re.compile(patron_cadenaADN) #La compilamos
    cadenaADN = '' #Inicializamos la variable que va a contener la cadena de ADN en cada iteración
    patron_juntarCadena = r'( |\n|\r|\t)' #Esta expresión regular la vamos a utilizar para sustituir los caracteres no deseados de las cadenas de ADN
    er_juntarCadena = re.compile(patron_juntarCadena) #La compilamos
    try:
        url = 'https://aulavirtual.um.es/access/content/group/1896_G_2022_N_N/PRACTICAS/PRACTICA%202/All_C_genes_DNA.txt'
        context = ssl.create_default_context()
        context.set_ciphers("DEFAULT")
        link = urllib.request.urlopen(url, context=context)
        for linea in link:
            linea = linea.decode().strip()
            nombreCadenaADN_res = er_nombreCadenaADN.search(linea) #Este es el booleano que nos va a indicar en cada caso si se ha encontrado el nombre de un gen
            cadenaADN_res = er_cadenaADN.search(linea) #Este es el booleano que nos va a indicar en cada caso si se ha encontrado una cadena de ADN en la línea
            numNucleotidos_res = er_numNucleotidos.search(linea) #Este es el booleano que nos va a indicar en cada caso si se ha encontrado un número de nucleótidos en la línea
            if nombreCadenaADN_res:#Si la linea contiene un nombre de cadena
                cadenaADN = er_juntarCadena.sub('', cadenaADN)  # entonces guardamos la cadena de ADN, la formatamos
                Dic_ADN[nombreCadenaADN] = [cadenaADN,numNucleotidos]  # Y metemos en el diccionario la entrada
                cadenaADN = ' '  # reseteamos la variable de cadena de adn a una cadena vacía
                nombreCadenaADN = linea[nombreCadenaADN_res.start():nombreCadenaADN_res.end()]#entonces guardamos el nombre de la cadena
                numNucleotidos = linea[numNucleotidos_res.start():numNucleotidos_res.end()]#guardamos el número de nucleótidos de ese gen
            elif cadenaADN_res:#elif si la linea contiene cadena de adn
                cadenaADN = cadenaADN+linea[cadenaADN_res.start():cadenaADN_res.end()]#la concatenamos a la variable cadena de adn existente
        cadenaADN = er_juntarCadena.sub('', cadenaADN) #Tenemos que realizar esta última iteración fuera del bucle porque el último gen no se almacena, ya que hemos hecho que almacenarlo en el diccionario se realice \\
        Dic_ADN[nombreCadenaADN] = [cadenaADN, numNucleotidos] #Cuando encontramos el siguiente gen
    except IOError as e:
        print('All_C_genes_DNA no disponible:', e)
    #for k in Dic_ADN.keys():  # Esto sirve símplemente para recorrer el diccionario de los genes
        #print('%s tiene valor %s' % (k, Dic_ADN[k]))

def mapaDeDianas(gen):
    cadenaADN = str(Dic_ADN.get(gen)) #Necesitamos la cadena de ADN como string para imprimirla y para buscar sobre ella
    enzimaIntroducida = input("Introduzca una enzima de reconocimiento") #Esta es la enzima que pedimos al usuario
    er_enzimaIntroducida = re.compile(enzimaIntroducida)#Tenemos que generar una expresión regular con la entrada por si no fuese una enzima contenida en el diccionario
    listadianas = [] #Esta será la lista de dianas sobre las que tenemos que realizar el mapa de dianas del gen introducido
    mapaDianas = [] #Aqui iremos iterando el mapa de dianas de cada diana
    if enzimaIntroducida == '': #Si la enzima introducida es la cadena vacía
        return #Salimos de la función
    else:
        if enzimaIntroducida in Dic_enzimas.keys(): #Si encontramos que coindide la enzima introducida con una de las claves del diccionario
            listadianas.append( enzimaIntroducida)  # Entonces se añade a la lista de dianas (No hacemos una variable separada porque sería malgastar espacio, simplemente luego no le añadimos mas entradas)
        else: #En caso de que no se encuentre la enzima en el diccionario
            for k in Dic_enzimas: #Recorremos de nuevo el diccionario
                if er_enzimaIntroducida.fullmatch(k): #Si la clave leída tiene coincidencia con la expresión regular introducida
                    listadianas.append(k) #Se añade a la lista esta enzima (La de la clave del diccionario)
    if listadianas == []:
        print("ERROR: No se ha encontrado ninguna enzima que cumpla con la expresión regular introducida")
        return
    print('Enzima >> ' + str(enzimaIntroducida)) #Imprimimos la enzima que ha introducido el usuario
    for l in listadianas: #Recorremos la lista de dianas
        diana = Dic_enzimas[l][0] #Cogemos la expresión regular de la diana
        posicionCorte = Dic_enzimas[l][1] #Cogemos la posición de corte de la diana
        for r in diana.finditer(cadenaADN): #Ahora iteramos la cadena de ADN buscando coincidencias con la ER de la diana escogida
            mapaDianas.append(r.start()+posicionCorte-2) #Añadimos la posición de corte al mapa de dianas (Ponemos -2 porque tenemos que descontar las posiciones 0 de las variables utilizadas)
        if not mapaDianas == []:
            print(l + ' # ' + str(mapaDianas)) #Imprimimos el nombre de la enzima y su mapa de dianas HpyAII   C.Ykr231ORF3053P
        mapaDianas = [] #Reiniciamos la variable del mapa de dianas / EL PROBLEMA QUE TENGO ES QUE TOMA LO INTRODUCIDO TANTO COMO DIANA COMO ER
def programa():
    gen = input("Introduzca el nombre del gen deseado") #Le pedimos al usuario el nombre del gen que quiere consultar
    if gen == '': #En caso de que sea la cadena vacía
        return #Salimos del programa
    if not (gen in Dic_ADN): #En caso de que no esté en el diccionario
        print("El gen introducido no se encuentra en el diccionario de genes") #Se lo decimos al usuario
        return #Y salimos del programa
    else: #En otro caso
        print(Dic_ADN.get(gen)[0]) #Imprimimos la cadena de nucleótidos del gen
        print(Dic_ADN.get(gen)[1] + " nucleótidos") #Y el número de nucleótidos
        mapaDeDianas(gen) #Y llamamos a la función mapaDeDianas pasándole como argumento el nombre del gen


leerEnzimas()
leerGenes()
programa()

#TENGO QUE AÑADIR EL RECORRIDO PARA CONVERTIR EN EXPRESIONES REGULARES COMPILADAS LAS CADENAS
#TENGO QUE HACER QUE LAS PRIMERAS LÍNEAS NO VÁLIDAS NO SE LEAN

#Proyecto final asignatura  de Autómatas y Lenguajes Formales, bioinformática

#Primero importamos los archivos necesarios de las bases de datos
#A continuación procesamos la base de datos para generar un diccionario de enzimas cuyas entradas son sus
#dianas de reconocimiento
#Para almacenar las dianas de reconocimiento antes tenemos que tratarlas

#Después tendremos que hacer una función que permita seleccionar una serie de dianas de reconocimiento
#y poder aplicarlas a una secuencia de ADN dada, devolviendo el mapa de dianas
