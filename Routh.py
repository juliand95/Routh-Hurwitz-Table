#//====================================================================================================================
#// UNIVERSIDAD PEDAGOGICA Y TECNOLOGICA DE COLOMBIA
#// FACULTAD SECCIONAL SOGAMOSO
#// ESCUELA DE INGENIERIA ELECTRONICA
#// MICROPROCESADORES
#// JULIAN DARIO BRAVO TUAY
#// TABLA DE ROUTH-HURWITZ
#// GENERA LA TABLA DE ROUTH-HURWITZ PARA UN POLINOMIO DE ORDEN N INGRESADO POR TECLADO
#// FECHA: 27/09/2020
#//====================================================================================================================

import numpy as np
import sys
ep = sys.float_info.epsilon     #Se declara epsilon
print("Genera la tabla de Routh-Hurwitz para un polinomio de orden n \n")
print("Se debe ingresar los coeficientes del polinomio de orden mayor a menor \n\n")
Tam = int(input("Ingrese el orden del polinomio:"))
Tam_Final = Tam + 1;            #Se le suma uno para ingresar todos los coeficientes
Polinomio = []                  #Se declara el vector donde se guardaran los coeficientes del polinomio
for i in range(Tam_Final):
    Polinomio.append(float(input(f"Grado del polinomio {Tam-i} : "))) #Se ingresan los coeficientes del polinomio y se guardan en el vector        
Polinomio = Polinomio[::-1]     #Se inverte el orden del vector del polinomio 

#Se declaran dos vectores para guarda los coeficientes pares e impares
Par = []
Impar = []

print()
#Se genera un ciclo for para verificar que coeficiente es para e impar para guardar
#los coeficientes en su respectivo vector.
for pos,i in enumerate(Polinomio):
    if pos%2==0:            #verifica si es par cada posicion
        Par.append(i)       #Se guardan los coeficientes pares
    else:
        Impar.append(i)     #Se guardan los coeficientes impares

Par = Par[::-1]             #Se invierte el orden del vector par
Impar = Impar[::-1]         #Se invierte el orden del vector impar

for pos,i in enumerate(Polinomio):
    if pos%2==0:                    #verifica si es par cada posicion
        Matriz = [Par, Impar]       #Se concatena los vectores par e impar 
    else:
        Matriz = [Impar, Par]       #Se concatena los vectores impar y par 
print("La matriz a resolver es: \n")
##******************************************************************************
if len(Par)==len(Impar):                    #Se verifica si la longitud de los coeficientes pares e impares son iguales
    #print("OPCION 1")
    Par.append(0)                           #se agrega un cero al vector par
    Impar.append(0)                         #se agrega un cero al vector impar
    M0 = np.zeros((len(Par),), dtype=int)   #Se crea un vector de ceros de longitud igual al de par o impar 
        
    for pos,i in enumerate(Polinomio):
        if pos%2==0:                                                #verifica si es par cada posicion
            Matriz = [Par, Impar]                                   #Se concatena los vectores par e impar 
            Matriz_Ceros = np.zeros((Tam_Final-2,len(Par)))         #Se crea una matriz de ceros 
            Matriz = np.concatenate((Matriz,Matriz_Ceros),axis=0)   #Se concatena la matriz principal con la matriz de ceros
            Matriz = Matriz.astype(float)                           #Se hace que los valores de la matriz sean de tipo float
        else:
            Matriz = [Impar, Par]                                   #Se concatena los vectores impar y par
            Matriz_Ceros = np.zeros((Tam_Final-2,len(Impar)))       #Se crea una matriz de ceros 
            Matriz = np.concatenate((Matriz,Matriz_Ceros),axis=0)   #Se concatena la matriz principal con la matriz de ceros
            Matriz = Matriz.astype(float)                           #Se hace que los valores de la matriz sean de tipo float
    print(Matriz)
    print()
    print("La matriz  resuelta de Routh-Hurwitz es: \n")
    #Se inicializan variables
    f = 0       
    f1 = 1
    c1 = 1
    c3 = 0
    f3 = 2
    
    #Se genera el for que calcula los coeficientes de las siguientes filas para generar la tabla de routh
    for i in range(len(Matriz)-2):              #Recorre las filas de la matriz
        for j in range(len(Matriz[i])-1):       #Recorre las columnas de la matriz
            det = -((Matriz[f][0] * Matriz[f1][c1])-(Matriz[f1][0] * Matriz[f][c1]))/Matriz[f1][0]      #Calcula cada coeficiente
            c1 = c1 + 1                         #se incrementa las columnas 
            Matriz[f3][c3] = round(det, 2)      #Guarda el coeficiente en su posicion correspondiente dentro de la matriz
            c3 = c3 + 1                         #se incrementa las columnas 
            det = 0                             #se reinicia para volver a capturar el nuevo resultado del coeficiente
            #Verifica si genera alguna fila completa de ceros
        if (Matriz[f3] == M0).all():
            print("Se genero una fila de ceros: ",Matriz[f3])
            print()
            P = Matriz[f3-1]                #se extrae los coeficientes de la fila inmediatamente anterior para derivar
            print("polinomio a derivar:",P)
            print()
            fila = f3-1
            pot = Tam-fila                  #se calcula la potencia a la cual se deriva los coeficientes
            print("De orden:",pot)
            print()
            dev = []                        #Se crea vector para guardar los coeficientes derivados
            #Se crea un ciclo for para calcular la derivada de cada coeficiente y guardarlos en un vector
            for i in range(len(P)):
                deri = pot*P[i]             #calcula la derivada de cada coeficiente del vector (Se multiplica la potencia por el coeficiente)
                pot = pot - 2               #Se resta 2 a la potencia para poder calcular la derivada del siguiente coeficiente (es generico)
                dev.append(deri)            #Se guarda el resultado de la derivada en el vector 
            print("La derivada es ",dev)
            Matriz[f3] = dev                #Se reemplaza la fila que se genero de ceros por los nuevos valores que se obtuvo de la derivada
            #Verifica si se genera un cero en la primera columna
        if Matriz[f3][0] == 0:
            print("Aparece un cero \n")
            Matriz[f3][0] = ep  #Se reemplaza el cero por epsilon

        f = f + 1       #se incrementa las filas una vez termina el ciclo de las columnas
        f1 = f1 + 1     #se incrementa las filas una vez termina el ciclo de las columnas
        f3 = f3 +1      #se incrementa las filas una vez termina el ciclo de las columnas
        c1 = 1          #se reinicia a su posicion actual
        c3 = 0          #se reinicia a su posicion actual

    print(Matriz)
##******************************************************************************    
elif len(Par)>=len(Impar):
    #print("OPCION 2")
    Par.append(0)
    Impar.append(0)
    Impar.append(0)
    print(Par)
    print(Impar)
    
    M0 = np.zeros((len(Par),), dtype=int)
    
    for pos,i in enumerate(Polinomio):
        if pos%2==0:                                                #verifica si es par cada posicion
            Matriz = [Par, Impar]                                   #Se concatena los vectores par e impar 
            Matriz_Ceros = np.zeros((Tam_Final-2,len(Par)))         #Se crea una matriz de ceros 
            Matriz = np.concatenate((Matriz,Matriz_Ceros),axis=0)   #Se concatena la matriz principal con la matriz de ceros
            Matriz = Matriz.astype(float)                           #Se hace que los valores de la matriz sean de tipo float
        else:
            Matriz = [Impar, Par]                                   #Se concatena los vectores impar y par
            Matriz_Ceros = np.zeros((Tam_Final-2,len(Impar)))       #Se crea una matriz de ceros 
            Matriz = np.concatenate((Matriz,Matriz_Ceros),axis=0)   #Se concatena la matriz principal con la matriz de ceros
            Matriz = Matriz.astype(float)                           #Se hace que los valores de la matriz sean de tipo float
    print(Matriz)
    print()
    print("La matriz  resuelta de Routh-Hurwitz es: \n")
    #Se inicializan variables
    f = 0
    f1 = 1
    c1 = 1
    c3 = 0
    f3 = 2
    
    #Se genera el for que calcula los coeficientes de las siguientes filas para generar la tabla de routh
    for i in range(len(Matriz)-2):              #Recorre las filas de la matriz
        for j in range(len(Matriz[i])-1):       #Recorre las columnas de la matriz
            det = -((Matriz[f][0] * Matriz[f1][c1])-(Matriz[f1][0] * Matriz[f][c1]))/Matriz[f1][0]      #Calcula cada coeficiente
            c1 = c1 + 1                         #se incrementa las columnas 
            Matriz[f3][c3] = round(det, 2)      #Guarda el coeficiente en su posicion correspondiente dentro de la matriz
            c3 = c3 + 1                         #se incrementa las columnas 
            det = 0                             #se reinicia para volver a capturar el nuevo resultado del coeficiente
            #Verifica si genera alguna fila completa de ceros
        if (Matriz[f3] == M0).all():
            print("Se genero una fila de ceros: ",Matriz[f3])
            print()
            P = Matriz[f3-1]                #se extrae los coeficientes de la fila inmediatamente anterior para derivar
            print("polinomio a derivar:",P)
            print()
            fila = f3-1
            pot = Tam-fila                  #se calcula la potencia a la cual se deriva los coeficientes
            print("De orden:",pot)
            print()
            dev = []                        #Se crea vector para guardar los coeficientes derivados
            #Se crea un ciclo for para calcular la derivada de cada coeficiente y guardarlos en un vector
            for i in range(len(P)):
                deri = pot*P[i]         #calcula la derivada de cada coeficiente del vector (Se multiplica la potencia por el coeficiente)
                pot = pot - 2           #Se resta 2 a la potencia para poder calcular la derivada del siguiente coeficiente (es generico)
                dev.append(deri)        #Se guarda el resultado de la derivada en el vector 
            print("La derivada es ",dev)
            Matriz[f3] = dev            #Se reemplaza la fila que se genero de ceros por los nuevos valores que se obtuvo de la derivada
            #Verifica si se genera un cero en la primera columna
        if Matriz[f3][0] == 0:
            print("Aparece un cero \n")
            Matriz[f3][0] = ep          #Se reemplaza el cero por epsilon

        f = f + 1       #se incrementa las filas una vez termina el ciclo de las columnas
        f1 = f1 + 1     #se incrementa las filas una vez termina el ciclo de las columnas
        f3 = f3 +1      #se incrementa las filas una vez termina el ciclo de las columnas
        c1 = 1          #se reinicia a su posicion actual
        c3 = 0          #se reinicia a su posicion actual
    
    print(Matriz)

##    pass
##******************************************************************************
elif len(Impar)>=len(Par):
    pass
