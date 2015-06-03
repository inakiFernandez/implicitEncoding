import matplotlib.pyplot as plt
import random
import numpy as np
import copy
import sys
from PIL import Image

#image processing
#from scipy import ndimage
from scipy import misc

dictSet=dict()

patterns=list()
patterns.append([[True,False,False,False,False]])
patterns.append([[True,True,True,False,False],[True,True,False,False,True],[True,False,False,True,True],[True,False,True,True,False]])
patterns.append([[True,True,False,True,False],[True,False,True,False,True]])
patterns.append([[True,True,False,False,False],[True,False,True,False,False],[True,False,False,True,False],[True,False,False,False,True]])
patterns.append([[True,True,True,True,False],[True,True,True,False,True],[True,False,True,True,True],[True,True,False,True,True]])
patterns.append([[True,True,True,True,True]])
patterns.append([[False,False,False,False,False]])

patternsDict=dict()
cnt=0
for i in range(len(patterns)):
	for j in range(len(patterns[i])):
		patt=patterns[i][j]
		patternsDict[str(patt)]=[i+1,j+1,-cnt]
		cnt=cnt+1

#####################
########TILE#########
#####################
class Cell:
	def __init__(self,idI,i,j):
		#initialization
		#tiles neighbors N E S W 0 1 2 3
		self.neighbors=[None]*4
		# True site False non site
		self.state= False
		#distance
		self.d1 = np.inf
		#init identifier
		self.idI=idI
		self.idInit = set([idI])
		#compute identifier
		self.idC = set(self.idInit)
		#area identifier
		self.idArea = set(self.idInit)
		#ij
		self.i = i
		self.j = j
		self.bis = 0
		self.pattern = [0,0,0]
		#direction ids
		self.ids = [set(self.idInit)]*5
		self.conflict = None

	#####getters
	def nothing(self):
		a=None

	#integer row i on the grid
	def getI(self):
		return self.i
	
	#integer column j on the grid
	def getJ(self):
		return self.j

	#integer 2-tuple list : von neumann neighbors (i,j)
	def getVonNeumannN(self):
		return self.neighbors

	#integer set : set of cell id
	def getID(self,default=4):
		if default==4:
			return self.idC
		else:
			if self.ids[default] is not None:
				return self.ids[default]
			else:
				return self.idC
		#if default == 4:
		#	return self.idC
		#else:
		#	return self.ids[default]

	#integer intrinsic id
	def getIDI(self):
		return self.idI

	#integer set: set of intrinsic id
	def getIDInit(self):
		return self.idInit

	#integer set : set of area id (one number if min, multiple numbers if union)
	def getIDArea(self):
		return self.idArea

	#string computed id : set -> representation -> hash -> key dictionary -> index
	def getID_str(self):
		has = hash(repr(self.getID()))
		if not str(has) in dictSet:
			dictSet[str(has)]=len(dictSet)+1
		return dictSet[str(has)]

	#string computed id string : set -> list -> ordered -> string -> key dictionary -> index
	def getID_str2(self):
		l=list(self.getID())
		l.sort()
		strL = '#'.join(str(l))
		if not strL in dictSet:
			dictSet[strL]=len(dictSet)+1
		return dictSet[strL]

	#integer 1-norm to the closest site
	def getD1(self):
		return self.d1

	#boolean site cell or not
	def isSite(self):
		return self.state

	#interger bisector state 0 :not 1:thin 2:thick
	def getBisector(self):
		return self.bis

	####setters
	def setVonNeumannN(self,nTab):
		self.neighbors = nTab

	def setID(self, idI, default=5):
		#self.idC=set(idI)
		if default == 5:
		 	self.idC = set(idI)
		 	self.ids[default-1]=set(idI)
		 	self.ids[default-2]=set(idI)
		 	self.ids[default-3]=set(idI)
		 	self.ids[default-4]=set(idI)
		 	self.ids[default-5]=set(idI)
		else:
		 	self.ids[default]=set(idI)

	def setIDArea(self, idA):
		self.idArea = set(idA)

	def setD1(self, d1):
		self.d1 = d1

	def setSite(self,boolean):
		self.state = boolean
		#if self.state == True:
			#self.setD1(0)
			#cself.setID(self.getIDInit())

	def getConflict(self):
		return self.conflict

	def getPattern(self):
		return self.pattern

	def getPattern_str(self):
		return self.pattern[2]

	def setPattern(self,class_value):
		self.pattern = class_value

	def setBisector(self,value):
		self.bis = value

	def setConflict(self,value):
		self.conflict = value

	def drawSiteAndBisector(self):
		if self.isSite():
			return 2
 		return self.bis

	def printCell(self):
		return "state:%s d1:%s id:%s idI: %s i:%s j:%s pattern:%s idN:%s" % (self.state,self.d1,self.idC,self.idInit, self.i, self.j, self.pattern,self.ids)

	def duplicate(self):
		cop = Cell(self.getIDI(),self.getI(),self.getJ())
		cop.setSite(self.isSite())
		cop.setD1(self.getD1())
		#cop.setID(copy.deepcopy(self.getID()))
		cop.idC = copy.deepcopy(self.idC)
		cop.ids = copy.deepcopy(self.ids)
		cop.setConflict(self.getConflict())
		#cop.setID(copy.deepcopy(self.getID(0)),0)
		#cop.setID(copy.deepcopy(self.getID(1)),1)
		#cop.setID(copy.deepcopy(self.getID(2)),2)
		#cop.setID(copy.deepcopy(self.getID(3)),3)
		#cop.setID(copy.deepcopy(self.getID(4)),4)
		
		cop.setIDArea(copy.deepcopy(self.getIDArea()))
		cop.setVonNeumannN(copy.deepcopy(self.getVonNeumannN()))
		cop.setBisector(self.getBisector())
		cop.setPattern(copy.deepcopy(self.getPattern()))
		return cop

##FUNCTIONS

#generate grid of tiles
def generate_grid(cell_type,lig,col):
	cellsIJ = [[cell_type(i*col+j+1,i,j) for j in range(col)] for i in range(lig)]
	for i in range(lig):
		for j in range(col):
			#print(cellsIJ[i][j].printCell())
			d=cellsIJ[i][j]
			if i>0: #north neighbor
				d.neighbors[0] = (i-1,j)#((i-1)*col)+j
			if j>0: #west neighbor
				d.neighbors[3] = (i,j-1)#(i*col) + (j-1)
			if j < col-1: #east neighbor
				d.neighbors[1] = (i,j+1)#(i*col) + (j+1)
			if i < lig-1: # south neighbor
				d.neighbors[2] = (i+1,j)#((i+1)*col)+j
	return cellsIJ

def draw_grid(cells,fun):
	lig = len(cells)
	col = len(cells[0])
	fig, ax = plt.subplots()
	image = [[fun(cells[i][j]) for j in range(col)] for i in range(lig)]
	ax.imshow(image, cmap=plt.cm.gray_r, interpolation='nearest')
	ax.set_title('Cellular Automaton')
	plt.show()
	raw_input('\npress return to continue\n')
	plt.close()


def draw_neighbors_cell(sizeLig,sizeCol,lig,col):
	fig, ax = plt.subplots()
	image = [[0 for j in range(sizeCol)] for i in range(sizeLig)]
	d = cellsIJ[lig][col]
	print("\n" + "neighbor %s %s"%(lig,col))
	print("\t" + d.printCell())
	neighbors = d.getVonNeumannN()
	for neighbor in neighbors:
		if neighbor is not None:
			nI, nJ = neighbor[0], neighbor[1]
			neighbor = cellsIJ[nI][nJ]
			print("\t\t" +neighbor.printCell())
			image[nI][nJ] = 1
	ax.imshow(image, cmap=plt.cm.gray, interpolation='nearest')
	ax.set_title('Cellular Automaton')
	plt.show(block=False)
	raw_input('\npress return to continue\n')
	plt.close()

def draw_neighbors_cell2(sizeLig,sizeCol):
	fig, ax = plt.subplots()
	image = [[0 for j in range(sizeCol)] for i in range(sizeLig)]
	im = ax.imshow(image, cmap=plt.cm.gray, interpolation='nearest')
	plt.show(block=False)
	raw_input('\npress return to continue\n')
	for i in range(sizeLig):
		for j in range(sizeCol):
			ax.set_title('neighbors %s %s' % (i,j))
			image = None
			image = [[0 for k in range(sizeCol)] for l in range(sizeLig)]
			d = cellsIJ[i][j]
			neighbors = d.getVonNeumannN()
			for neighbor in neighbors:
				if neighbor is not None:
					nI, nJ = neighbor[0], neighbor[1]
					image[nI][nJ] = 1
					ax.clear()
					ax.matshow(image, cmap=plt.cm.gray, interpolation='nearest')

			plt.show(block=False)
			raw_input('\npress return to continue\n')		
			
	plt.close(fig)


def copyAutomaton(cells,cell_type):
	lig = len(cells)
	col = len(cells[0])
	automaton=None
	automaton=[[cell_type(i*col+j,i,j) for j in range(col)] for i in range(lig)]
	for i in range(lig):
		for j in range(col):
			automaton[i][j] = cells[i][j].duplicate()
	return automaton

#TODO
def simple_synchronism(cells,funC,funD,funP,cm=plt.cm.gray):
	lig = len(cells)
	col = len(cells[0])
	#init drawing
	fig, ax = plt.subplots()

	while True:
		#cellsIJ T copy
		snapshot = copyAutomaton(cells,Cell)
		#work on the copied automaton and update the next state automaton
		for i in range(lig):
			for j in range(col):
				funC(i,j,snapshot,cells)

		#printing
		imagePrint = [[ funP(cells[l][k]) for k in range(col) ] for l in range(lig) ]
		print(np.matrix(imagePrint))
		
		#drawing
		imageDraw = [[ funD(cells[l][k]) for k in range(col) ] for l in range(lig) ]
		ax.clear()
		#ax.matshow(imageDraw, cmap=plt.cm.gray, interpolation='nearest')
		ax.matshow(imageDraw, cmap=cm, interpolation='nearest')
		
		plt.show(block=False)
		raw_input('\npress return to continue\n')	
	
	plt.close(fig)

def simple_asynchronism(cells,funC,funD,funP,cm=plt.cm.gray):
	sizeLig = len(cells)
	sizeCol = len(cells[0])
	
	#init drawing
	fig, ax = plt.subplots()
	imageDraw = [[ funD(cells[l][k]) for k in range(sizeCol) ] for l in range(sizeLig) ]
	imagePrint = [[ funP(cells[l][k]) for k in range(sizeCol) ] for l in range(sizeLig) ]
	cnt = 0
	while True:
		#random execution
		randI = random.randint(0,sizeLig-1)
		randJ = random.randint(0,sizeCol-1)
		#print("cell %s %s"%(randI,randJ))
		funC(randI,randJ,cells,cells)
		cnt = cnt+1
		
		imagePrint[randI][randJ] = funP(cells[randI][randJ])
		imageDraw[randI][randJ] = funD(cells[randI][randJ])
		if cnt % 10000000 == 0:
			#printing
			print(np.matrix(imagePrint))
		
			#drawing
			ax.clear()
			ax.matshow(imageDraw, cmap=cm, interpolation='nearest')
			plt.show(block=False)
			raw_input('\npress return to continue\n')		

	plt.close(fig)

def computePatterns(cells):
	lig = len(cells)
	col = len(cells[0])
	for i in range(lig):
		for j in range(col):
			patternFound=False
			cell = cells[i][j]
			vN = cell.getVonNeumannN()
			localPat = [False]*5
			if cell.isSite():
				localPat[0] = cell.isSite()
				for k in range(len(vN)):
					if vN[k] is not None:
						vNK = vN[k]
						localPat[k+1] = cells[vNK[0]][vNK[1]].isSite()
				pat_Cl_Nb=patternsDict[str(localPat)]
				print(pat_Cl_Nb)
				cell.setPattern(pat_Cl_Nb)
	return cells

################
def use_case_test():
	cellsIJ = generate_grid(Cell,20,20)
	s1=cellsIJ[1][18]
	s2=cellsIJ[18][1]
	s1.setSite(True)
	s2.setSite(True)
	return cellsIJ

#72 36
def bis_case_1():
	cellsIJ = generate_grid(Cell,8,8)
	s1=cellsIJ[6][1]
	s2=cellsIJ[2][5]
	s1.setSite(True)
	s2.setSite(True)
	return cellsIJ

def bis_case_2():
	cellsIJ = generate_grid(Cell,8,8)
	s1=cellsIJ[6][1]
	s2=cellsIJ[1][6]
	s1.setSite(True)
	s2.setSite(True)
	return cellsIJ

def bis_case_3():
	cellsIJ = generate_grid(Cell,8,8)
	s1=cellsIJ[6][2]
	s2=cellsIJ[0][4]
	s1.setSite(True)
	s2.setSite(True)
	return cellsIJ

def bis_case_4():
	cellsIJ = generate_grid(Cell,8,8)
	s1=cellsIJ[6][2]
	s2=cellsIJ[0][5]
	s1.setSite(True)
	s2.setSite(True)
	return cellsIJ

def bis_case_5():
	cellsIJ = generate_grid(Cell,8,8)
	s1=cellsIJ[6][2]
	s2=cellsIJ[1][4]
	s1.setSite(True)
	s2.setSite(True)
	return cellsIJ

def bis_case_6():
	cellsIJ = generate_grid(Cell,8,8)
	s1=cellsIJ[6][3]
	s2=cellsIJ[1][4]
	s1.setSite(True)
	s2.setSite(True)
	return cellsIJ

#create first map
def map1():
	#tine num = (x+(y)*lig) y*lig + x
	lig = 32
	col = 32
	cellsIJ=generate_grid(Cell,lig,col)
	tab = []
	#boolean_mobile = False
	
	#map 32*32
	#obstacle 1
	
	#first rectangle
	for i in range(9):
			for j in range(2):
				#obstacle(4*lig+4+i+j*lig)
				tab.append(4*lig+4+i+j*lig)
	
	#second rectangle
	for i in range(2):
		for j in range(4):
			#obstacle(6*lig+4+i+j*lig)
			tab.append(6*lig+4+i+j*lig)
	
	#obstacle 2
	for i in range(9):
		for j in range(3):
			#obstacle(4*lig+18+i+j*lig)	
			tab.append(4*lig+18+i+j*lig)
	
	#obstacle 3
	for i in range(7):
		for j in range(2):
			#obstacle(9*lig+9+i+j*lig)
			tab.append(9*lig+9+i+j*lig)
	for i in range(5):
			#obstacle(11*lig+11+i)
			tab.append(11*lig+11+i)
	for i in range(4):
			#obstacle(12*lig+12+i)
			tab.append(12*lig+12+i)
	for i in range(3):
			#obstacle(13*lig+13+i)
			tab.append(13*lig+13+i)
	
	#obstacle 4
	for i in range(7):
		for j in range(3):
			#obstacle(10*lig+21+i+j*lig)		
			tab.append(10*lig+21+i+j*lig)
	for i in range(5):
			#obstacle(13*lig+21+i)
			tab.append(13*lig+21+i)
	for i in range(3):
			#obstacle(14*lig+22+i)
			tab.append(14*lig+22+i)
	#for i in range(1):
		#	obstacle(15*lig+23+i)

	#obstacle 5
	for i in range(4):
		for j in range(3):
			#obstacle(16*lig+18+i+j*lig)
			tab.append(16*lig+18+i+j*lig)

	#obstacle 6
	for i in range(8):
		for j in range(3):
			#obstacle(25*lig+4+i+j*lig)		
			tab.append(25*lig+4+i+j*lig)

	#obstacle 7
	for i in range(9):
		for j in range(6):
			#obstacle(22*lig+19+i+j*lig)
			tab.append(22*lig+19+i+j*lig)

	#obstacle 8
	for i in range(4):
		#obstacle(13*lig+5+i)
		tab.append(13*lig+5+i)

	for i in range(7):
		for j in range(7):
			#obstacle(14*lig+4+i+j*lig)
			tab.append(14*lig+4+i+j*lig)
	for i in range(4):
		#obstacle(21*lig+5+i)
		tab.append(21*lig+5+i)
	for i in range(2):
		for j in range(3):
			#obstacle(16*lig+11+i+j*lig)
			tab.append(16*lig+11+i+j*lig)

	#obstacle(15*lig+10+i)	
	#obstacle(19*lig+10+i)
	tab.append(15*lig+10+i)
	tab.append(19*lig+10+i)
	
	#obstacle border
	for i in range(32):
		#obstacle(i)
		#obstacle(31*lig+i)
		tab.append(i)
		tab.append(31*lig+i)		
		if i !=0 and i!=31:
			#obstacle(32*i)
			#obstacle(31+32*i)
			tab.append(32*i)
			tab.append(31+32*i)

	for t in tab:
			a=t/col
			b=t%col
			s=cellsIJ[a][b]
			s.setSite(True)
	
	return cellsIJ

# #create second map
def map2():
	lig = 32
	col = 32
	#boolean_mobile = False
	cellsIJ=generate_grid(Cell,lig,col)
	tab = []
	#map 27*27
	#obstacle border
	for i in range(32):
		#obstacle(i)
		#obstacle(lig+i)
		#obstacle(31*lig+i)
		#obstacle(30*lig+i)
		#obstacle(29*lig+i)
		#obstacle(28*lig+i)
		tmp = [(i),(lig+i),(31*lig+i),(30*lig+i),(29*lig+i),(28*lig+i)]
		tab.extend(tmp) 
		if i !=0 and i!=1 and i!=31 and i!=30 and i!=29 and i!=28:
			#obstacle(32*i)
			#obstacle(32*i+1)
			#obstacle(32*i+2)
			#obstacle(32*i+3)
			#obstacle(31+32*i)
			#obstacle(30+32*i)
			tmp = [(32*i),(32*i+1),(32*i+2),(32*i+3),(31+32*i),(30+32*i)]
			tab.extend(tmp) 
		
		for i in range(10):
			for j in range(2):
				#obstacle((26+j)*lig+4+i)
				tab.append((26+j)*lig+4+i) 
		
		for i in range(5):
			for j in range(4):
				#obstacle((22+j)*lig+9+i)
				tab.append((22+j)*lig+9+i)
		
		for i in range(3):
			for j in range(5):
				#obstacle((17+j)*lig+9+i)
				tab.append((17+j)*lig+9+i)

		for i in range(6):
			for j in range(18):
				#obstacle((7+j)*lig+17+i)
				tab.append((7+j)*lig+17+i)

		for i in range(2):
			for j in range(11):
				#obstacle((2+j)*lig+28+i)
				tab.append((2+j)*lig+28+i)

	for t in tab:
			a=t/col
			b=t%col
			s=cellsIJ[a][b]
			s.setSite(True)
	
	return cellsIJ

#create third map
def map3():
	lig = 30
	col = 30
	#boolean_mobile = False
	cellsIJ=generate_grid(Cell,lig,col)
	tab = []

	#map 27*27
	for i in range(30):
		#obstacle(i)
		#obstacle(lig+i)
		#obstacle(2*lig+i)
		#obstacle(29*lig+i)
		#obstacle(27*lig+i)
		#obstacle(28*lig+i)
		tmp = [i,lig+i,2*lig+i,29*lig+i,27*lig+i,28*lig+i]
		tab.extend(tmp)

		if i !=0 and i!=1 and i!=2 and i!=29 and i!=27 and i!=28:
			#obstacle(30*i)
			#obstacle(30*i+1)
			#obstacle(30*i+2)
			#obstacle(30*i+lig)
			#obstacle(30*i+1+lig)
			#obstacle(30*i+2+lig)
			tmp = [30*i,30*i+1,30*i+2,30*i+lig,30*i+1+lig,30*i+2+lig]
			tab.extend(tmp)
			
			if i not in [0,1,2,]:
				#obstacle(29+30*i)
				#obstacle(28+30*i)
				#obstacle(27+30*i)
				tmp = [29+30*i,28+30*i,27+30*i]
				tab.extend(tmp)
			
			if i not in [3,4,5]:				
				#obstacle(23+30*i)
				#obstacle(22+30*i)
				#obstacle(21+30*i)
				tmp = [23+30*i,22+30*i,21+30*i]
				tab.extend(tmp)
				
		for i in range(9):
			for j in range(3):
				#obstacle((6+j)*lig+12+i)
				tab.append((6+j)*lig+12+i)
		for i in range(3):
			for j in range(18):
				#obstacle((3+j)*lig+6+i)
				tab.append((3+j)*lig+6+i)
		for i in range(9):
			for j in range(3):
				#obstacle((18+j)*lig+9+i)
				tab.append((18+j)*lig+9+i)
		for i in range(3):
			for j in range(6):
				#obstacle((12+j)*lig+15+i)
				tab.append((12+j)*lig+15+i)
		for i in range(3):
			for j in range(3):
				#obstacle((12+j)*lig+12+i)
				tab.append((12+j)*lig+12+i)

	for t in tab:
			a=t/col
			b=t%col
			s=cellsIJ[a][b]
			s.setSite(True)
	
	return cellsIJ

#create third map
def map4():
	lig = 30
	col = 30
	cellsIJ=generate_grid(Cell,lig,col)
	tab = []
	#map 27*27
	for i in range(30):
		#obstacle(i)
		#obstacle(lig+i)
		#obstacle(2*lig+i)
		#obstacle(29*lig+i)
		#obstacle(28*lig+i)
		#obstacle(27*lig+i)	
		tmp = [i,lig+i,2*lig+i,29*lig+i,28*lig+i,27*lig+i]
		tab.extend(tmp)

		if i !=0 and i!=1 and i!=2 and i!=29 and i!=28 and i!=27:
		 	#obstacle(30*i)
		 	#obstacle(30*i+1)
		 	#obstacle(30*i+2)
		 	#obstacle(30*i+29)
		 	#obstacle(30*i+28)
		 	#obstacle(30*i+27)
		 	tmp = [30*i,30*i+1,30*i+2,30*i+29,30*i+28,30*i+27]
		 	tab.extend(tmp)

		 	if i not in [3,4,5]:	
				tmp = [23+30*i,22+30*i,21+30*i]
				tab.extend(tmp)

		for i in range(9):
		 	for j in range(3):
		 		#obstacle((6+j)*lig+12+i)
				tab.append(((6+j)*lig+12+i))

		for i in range(3): #PB
		 	for j in range(3):
		 		#obstacle((12+j)*lig+9+i)
				tab.append(((12+j)*lig+9+i))
		 		
		for i in range(9):
		 	for j in range(3):
		 		#obstacle((20+j)*lig+9+i)
				tab.append(((20+j)*lig+9+i))
		
		for i in range(3):
		 	for j in range(6):
		 		#obstacle((12+j)*lig+15+i)
				tab.append(((12+j)*lig+15+i))

		for i in range(3):
		 	for j in range(3):
		 		#obstacle((12+j)*lig+12+i)
				tab.append(((12+j)*lig+12+i))

		for j in range(2):
		 	for i in range(3):
		 		#obstacle((18+j)*lig+15+i)
				tab.append(((18+j)*lig+15+i))

	for t in tab:
		a=t/col
		b=t%col
		s=cellsIJ[a][b]
		s.setSite(True)
	
	return cellsIJ

#create third map
def map5():
	lig=79
	col=28
	cellsIJ=generate_grid(Cell,lig,col)
	tab = []

	#boolean_mobile = False
	#border of map
	
	#for i in range(col):
		#obstacle(i)
		#obstacle(i+(lig*(lig-1)))
	#	tab.append(i)
	#	tab.append(i+(lig*(lig-1)))

	#for i in range(lig):
		#obstacle(lig*i)
		#obstacle((lig*i)+col-1)
	#	tab.append(lig*i)
	#	tab.append((lig*i)+col-1)

	#map 80*80
	#obstacle(26+20*lig)
	#obstacle(25+78*lig)
	#obstacle(21+42*lig)
	#obstacle(21+74*lig)
	#obstacle(19+12*lig)
	#obstacle(19+21*lig)
	#obstacle(14+25*lig)
	#obstacle(13+56*lig)
	#obstacle(9+14*lig)
	#obstacle(7+39*lig)
	#obstacle(7+70*lig)
	#obstacle(5+10*lig)
	#obstacle(4+25*lig)
	#obstacle(3+69*lig)
	#obstacle(1+78*lig)

	tmp = [26+20*lig,25+78*lig,21+42*lig,21+74*lig,19+12*lig,19+21*lig,14+25*lig,13+56*lig,9+14*lig,7+39*lig,7+70*lig,5+10*lig,4+25*lig,3+69*lig,1+78*lig]
	tab.extend(tmp)
	for t in tab:
		a=t/lig
		b=t%lig
		s=cellsIJ[a][b]
		s.setSite(True)
	return cellsIJ
	

# #create sixth map
# def map6(lig,col):
# 	boolean_mobile = False
# 	#border of map
# 	for i in range(col):	
# 		obstacle(i)
# 	for i in range(lig):
# 		obstacle((col*(lig-1))+i)
# 	for i in range(col-2):	
# 		obstacle(lig+i*lig)
# 	for i in range(col-2):
# 		obstacle(2*lig-1+i*lig)

# 	for i in range(lig):
# 		obstacle(i+i*lig)
# 		if i != lig-1:
# 			obstacle(1+i+i*lig)

# #create seventh map
# def map7(lig, col):
# 	boolean_mobile = False
# 	#border of map
# 	for i in range(col):	
# 		obstacle(i)
# 	for i in range(col):	
# 		obstacle(i)
# 	for i in range(lig):
# 		obstacle((col*(lig-1))+i)
# 	for i in range(col-2):	
# 		obstacle(lig+i*lig)
# 	for i in range(col-2):
# 		obstacle(2*lig-1+i*lig)
# 	for i in range(lig):
# 		obstacle(i+i*lig)

def map8():
	lig = 20
	col = 20
	cellsIJ=generate_grid(Cell,lig,col)
	#border of map
	for i in range(col):	
		#obstacle(i)
		s=cellsIJ[0][i]
		s.setSite(True)
		s=cellsIJ[lig-1][i]
		s.setSite(True)
		s=cellsIJ[i][0]
		s.setSite(True)
		s=cellsIJ[i][col-1]
		s.setSite(True)

	#incomplete square
	#first point
	for i in range(col/4):
		#obstacle(col*(3*lig/8) + (3*col/8) +i*lig)
		a=(col*(3*lig/8) + (3*col/8) +i*lig)/col
		b=(col*(3*lig/8) + (3*col/8) +i*lig)%col
		s=cellsIJ[a][b]
		s.setSite(True)
	
	for i in range(col/4):
		#obstacle(col*(5*lig/8) + (3*col/8) +i)
		a=(col*(5*lig/8) + (3*col/8) +i)/col
		b=(col*(5*lig/8) + (3*col/8) +i)%col
		s=cellsIJ[a][b]
		s.setSite(True)
	
	for i in range(col/8 + 1):
		#obstacle(col*(3*lig/8) + (3*col/8) +i)
		a=(col*(3*lig/8) + (3*col/8) +i)/col
		b=(col*(3*lig/8) + (3*col/8) +i)%col
		s=cellsIJ[a][b]
		s.setSite(True)

	for i in range(col/8 + 1):
		#obstacle(col*(5*lig/8) + (5*col/8) -(i*lig))
		a=(col*(5*lig/8) + (5*col/8) -(i*lig))/col
		b=(col*(5*lig/8) + (5*col/8) -(i*lig))%col
		s=cellsIJ[a][b]
		s.setSite(True)
	return cellsIJ

def map0():
	lig = 32
	col = 32
	cellsIJ=generate_grid(Cell,lig,col)
	#borders
	for i in range(col):	
		s=cellsIJ[0][i]
		s.setSite(True)
		s=cellsIJ[lig-1][i]
		s.setSite(True)
		s=cellsIJ[i][0]
		s.setSite(True)
		s=cellsIJ[i][col-1]
		s.setSite(True)
	return cellsIJ

def mapTest():
	lig = 20
	col = 20
	cellsIJ=generate_grid(Cell,lig,col)
	#borders
	
	for i in range(col):	
		s=cellsIJ[0][i]
		s.setSite(True)
		s=cellsIJ[lig-1][i]
		s.setSite(True)
		s=cellsIJ[i][0]
		s.setSite(True)
		s=cellsIJ[i][col-1]
		s.setSite(True)
	
	#Case 4 5
	"""
	s=cellsIJ[lig/2][col/2]
	s.setSite(True)
	s=cellsIJ[lig/2+1][col/2]
	s.setSite(True)
	s=cellsIJ[lig/2+2][col/2]
	s.setSite(True)
	"""

	#Case2
	"""
	s=cellsIJ[lig/2][col/2]
	s.setSite(True)
	s=cellsIJ[lig/2][col/2+1]
	s.setSite(True)
	s=cellsIJ[lig/2][col/2+2]
	s.setSite(True)
	s=cellsIJ[lig/2][col/2+3]
	s.setSite(True)
	s=cellsIJ[lig/2][col/2+4]
	s.setSite(True)
	s=cellsIJ[lig/2][col/2+5]
	s.setSite(True)
	s=cellsIJ[lig/2-1][col/2]
	s.setSite(True)
	s=cellsIJ[lig/2-2][col/2]
	s.setSite(True)
	s=cellsIJ[lig/2-3][col/2]
	s.setSite(True)
	"""
	
	#Case3
	"""
	s=cellsIJ[lig/2][col/2]
	s.setSite(True)
	s=cellsIJ[lig/2][col/2+1]
	s.setSite(True)
	s=cellsIJ[lig/2][col/2+2]
	s.setSite(True)
	s=cellsIJ[lig/2+1][col/2+1]
	s.setSite(True)
	"""
	"""
	#Case1
	s=cellsIJ[lig/2][col/2-3]
	s.setSite(True)
	s=cellsIJ[lig/2][col/2-2]
	s.setSite(True)
	s=cellsIJ[lig/2][col/2-1]
	s.setSite(True)
	s=cellsIJ[lig/2][col/2]
	s.setSite(True)
	s=cellsIJ[lig/2][col/2+1]
	s.setSite(True)
	s=cellsIJ[lig/2][col/2+3]
	s.setSite(True)
	s=cellsIJ[lig/2][col/2+2]
	s.setSite(True)
	
	s=cellsIJ[lig/2+1][col/2-3]
	s.setSite(True)
	s=cellsIJ[lig/2+1][col/2-2]
	s.setSite(True)
	s=cellsIJ[lig/2+1][col/2-1]
	s.setSite(True)
	s=cellsIJ[lig/2+1][col/2]
	s.setSite(True)
	s=cellsIJ[lig/2+1][col/2+1]
	s.setSite(True)
	s=cellsIJ[lig/2+1][col/2+2]
	s.setSite(True)
	s=cellsIJ[lig/2+1][col/2+3]
	s.setSite(True)
	"""

	#lines
	max=20
	#for i in range(-max,max):
	#	s=cellsIJ[lig/2+i][col/2]
		#s.setSite(True)
	
	max=20
	#for i in range(-max,max):
	#	s=cellsIJ[lig/2][col/2+i]
		#s.setSite(True)

	#rectangle
	#max=lig/5
	#for i in range(max/2):
	#	for j in range(max):
	#		s=cellsIJ[lig/4+i][3*col/4+j]
			#s.setSite(True)

	#triangle
	max=4
	for i in range(max):
		for j in range(i):
			s=cellsIJ[lig/2+i][col/2+j]
			s.setSite(True)

	#triangle
	max=lig/5
	for i in range(max):
		for j in range(i):
			s=cellsIJ[7*lig/10+i][9*col/10-j]
			#s.setSite(True)

	#trapeze
	max=lig/5
	for i in range(max):
		for j in range(i+1):
			if i!=0 and i!=max-1:
				s=cellsIJ[lig/4+i][col/4+j]
				#s.setSite(True)
			if i!=0 and j!=0:
				s=cellsIJ[lig/4+i][col/4+j-1]
				#s.setSite(True)

	#square
	max=lig/10
	for i in range(max):
		for j in range(max):
			s=cellsIJ[lig/4+i][3*col/4+j]
			#s.setSite(True)

	for i in range(lig):# for every pixel:
		for j in range(col):
			if random.random() > 0.995:
				s=cellsIJ[i][j]
				#s.setSite(True)

	return cellsIJ

def cellsFromImg():
	#img = Image.new( 'RGB', (255,255), "black") # create a new black image
	#pixels = img.load() # create the pixel map
	#lig = img.size[0]
	#col = img.size[1]
	lig = 100
	col = 100
	cellsIJ=generate_grid(Cell,lig,col)
	for i in range(lig):# for every pixel:
		for j in range(col):
			if random.random() > 0.99:
				s=cellsIJ[i][j]
				s.setSite(True)
			    #pixels[i,j] = (i, j, 100) # set the colour accordingly
				#img.show()
	return cellsIJ

def read_imag(fileName):
	# file = open(fileName,"rb")
	# data = file.read()
	# file.close()
	# print(data)	
	
	#img = Image.new( 'RGB', (255,255), "black") # create a new black image
	#pixels = img.load() # create the pixel map
	#lig = img.size[0]
	#col = img.size[1]
	#img.show()
	#for i in range(lig):# for every pixel:
	#	for j in range(col):
			#pixels[i,j] = (i, j, 100) # set the colour accordingly
			
	#return cellsIJ
	#im = Image.open(fileName)
	#im.rotate(45).show()

	data = misc.imread(fileName, cmap=plt.cm.gray, interpolation='nearest')
	#l=misc.data()
	#plt.imshow(data, cmap=plt.cm.gray)


def read_img1(fileName):
	data = misc.imread(fileName)
	plt.imshow(data,cmap=plt.cm.gray)
	plt.show()

def read_img2(fileName):#PIL
	im = Image.open(fileName)
	(w,h) = im.size
	im = im.resize((w,h))#/10 IS GOOD
	(w,h) = im.size
	data = im.getdata()
	print(w,h)
	cellsIJ=generate_grid(Cell,h,w)
	for i in range(h):# for every pixel:
		for j in range(w):
			#print(i,j)
			#print(im.getpixel((i,j)))
			if im.getpixel((j,i)) == 0 or im.getpixel((j,i)) == (0,0,0) :
				s=cellsIJ[i][j]
				s.setSite(True)
	#im.show()
	return cellsIJ

