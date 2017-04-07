# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 15:49:54 2017

@author: Amalik
"""

import wave
import binascii

NomFichier = "3Bonjours2.wav"
Monson = wave.open(NomFichier,'r')	# instanciation de l'objet Monson
Monson.setpos(44)
for i in range(100):
    a = Monson.readframes(1)
    print(a[0],a[1],a[2],a[3])
    b = a[0]*(2**8)+a[1]
    b = ((b/((2**16-1)))*2)-1
    print("b",b)
print("\nNombre de canaux :",Monson.getnchannels())
print("Taille d'un échantillon (en octets):",Monson.getsampwidth())
print("Fréquence d'échantillonnage (en Hz):",Monson.getframerate())
print("Nombre d'échantillons :",Monson.getnframes())
print("Type de compression :",Monson.getcompname())

echDebut = 0
echFin = 2

print("\nN° échantillon	Contenu")

Monson.setpos(echDebut)
plage = echFin - echDebut + 1
for i in range(0,plage):
    print(Monson.tell(),'\t\t',binascii.hexlify(Monson.readframes(4)))

#Monson.close()

a = int("33fd38fc",16)-((2**16-1)/2)
print(a)
