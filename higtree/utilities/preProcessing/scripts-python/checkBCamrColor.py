#!/usr/bin/env python
#coding: utf-8

import os
import os.path
import sys
import numpy as np
from numpy import *
import mmap
import glob

#Verificar instalação do termcolor sudo apt-get install python-termcolor
from termcolor import colored

def checkDir(filesDir, s2):
	print 'function-------'
	if os.path.exists(filesDir):
		for filename in glob.glob(os.path.join(filesDir, '*_bc')):
			#print filename
			cont = 0
			cont2 = 0
			for amrFile in glob.glob(os.path.join(filename, '*.amr')):
				#print amrFile
				path_amr = os.path.join(filename, amrFile)
				f = open(path_amr)
				#s = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
				#if s.find(str(stringValues)) != -1:
				#if s.find('   0.0000000000    0.0000000000    0.3000000000    0.4000000000    0.0000000000    0.1000000000 ') != -1:
				first_line = f.readline()
				sA = np.array(first_line.split())
				s3 = sA.astype(np.float)
				#s3[-1] = s3[-1].strip()

				#print s3
				#print s3[0]
				#print type(s3)

				t = np.array_equal(s2,s3)

				if (t):
					print 'True'
					print colored(amrFile, 'red', attrs=['reverse', 'blink'])
					cont = cont + 1
				cont2 = cont2+1
				f.close()
			if cont == 0:
				print 'Nenhum arquivo amr encontrado...\n'
				pass
			print 'Quantidade de amrs files verificados = ', cont2
	else:
		print '\nNão existe o diretório: ', filesDir


baseDir = os.getcwd()

print '\n ----- Verificação arquivos amrs para as condições de contorno\n'

s = raw_input('Digite o nome do diretório a ser verificado: ')
print 'Confirmando nome do diretório: ', s

print baseDir

if not os.path.exists(str(baseDir) +'/'+ str(s)):
	print '\nDiretório não existe: ', s
	print '\n --- Digite novamente'
	sys.exit('Error')

filesDir = str(baseDir) +'/'+ str(s) +'/amrs/'

print filesDir

dim = int(input ('Digite a dimensão do problema [3] or [2]: '))

if (dim == 2) or (dim == 3):
	print 'Valores dos pontos das condições de contorno à serem verificadas no domínio: xb xu yb yu zb zu'

	#xl = raw_input('Digite o valor para x inferior: ')

	while True:
		try:
			xb = float(input ('\tDigite valor para x inferior: '))
			xu = float(input ('\tDigite valor para x superior: '))
			yb = float(input ('\tDigite valor para y inferior: '))
			yu = float(input ('\tDigite valor para y superior: '))
			if(dim == 3):
				zb = float(input ('\tDigite valor para z inferior: '))
				zu = float(input ('\tDigite valor para z superior: '))
			else:
				pass
		except ValueError:
			print('When I ask for a number, give me a number. Come on!')
		else:
			if (dim == 2):
				stringValues = '   '+str(xb)+'    '+str(xu)+'    '+str(yb)+'    '+str(yu)+' '
			if (dim == 3):
				stringValues = '   '+str(xb)+'    '+str(xu)+'    '+str(yb)+'    '+str(yu)+'    '+str(zb)+'    '+str(zu)+' '
			sA1 = np.array(stringValues.split())
			s2 = sA1.astype(np.float)
			print s2
			#print s2[0]
			#print type(s2)
			print'---call function...\n'
			checkDir(filesDir, s2)
			break
else:
	print 'Dimensão diferente de 2 ou 3'










#f = open('/media/jonas/Dropbox/projICMC/paper/ns-3d_test02/amrs/Channel_3d_bc/mesh-channel3D-bc-1.amr')
#s = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
#if s.find('   0.0000000000    0.0000000000    0.3000000000    0.4000000000    0.0000000000    0.1000000000 ') != -1:
#	print 'True'
#else:
#	print 'False'



