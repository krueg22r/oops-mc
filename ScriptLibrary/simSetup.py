#A program that creates the topology and coord data files for running a simulation
#At this point, it creates a box of given dimension filled with identical polymers w/config 
#specified in an xyz file. It does this by calling program packmol as an intermediate step


#!/usr/local/bin/python3.4
from subprocess import check_output

def main(): 
	polFile = input().split()[1] #name of xyz file with polymer structure
	inpFile = input().split()[1] #name of input file for my prog that this prog will write
	topFile = input().split()[1] #name of topology file for y prog that this prog will write
	nPol = int(input().split()[1]) #number of polymers we eventually want in box
	spacing = input().split()[1] #min particle spacing in file
	xLen = float(input().split()[1])
	yLen = float(input().split()[1])
	zLen = float(input().split()[1])
	#read number beads per polymer from single polymer file
	readPol = open(polFile,'r')
	nBead = int(readPol.readline())
	#print("nBead = "+str(nBead))
	#variables we'll want if this is a charge system
	nPlus = 0
	nMinus = 0
	chainCharge = 0
	chargeFrac = 1
	plusFile = 'sdfsd' 
	minusFile = 'sdlfj' 
	chargeFlag = int(input().split()[1]) #do we have charges in this run? if so, read in vals for all of the above stuff
	if(chargeFlag == 1): 
		nPlus = int(input().split()[1])
		nMinus = int(input().split()[1])
		plusFile = input().split()[1]
		minusFile = input().split()[1]
		chainCharge = int(input().split()[1])
		chargeFrac = float(input().split()[1])
	packInp = open('in.packmol','w') #open file for writing input to packmol
	#write packmol input file (in.packmol) 
	print("tolerance "+spacing,file=packInp)
	print("output packOut.xyz",file=packInp)
	print("filetype xyz",file=packInp)
	print("structure "+polFile,file=packInp)
	print("  number "+str(nPol),file=packInp)
	print("  inside box "+spacing+" "+spacing+" "+spacing+" "+str(xLen-2)+" "+str(yLen-2)+" "+str(zLen-2),file=packInp)
	print("end structure",file=packInp)
	#now, if we need ions also: 
	if(chargeFlag==1): 
		print("structure "+plusFile,file=packInp)
		print("  number "+str(nPlus),file=packInp)
		print("  inside box "+spacing+" "+spacing+" "+spacing+" "+str(xLen-2)+" "+str(yLen-2)+" "+str(zLen-2),file=packInp)
		print("end structure",file=packInp)
		print("structure "+minusFile,file=packInp)
		print("  number "+str(nMinus),file=packInp)
		print("  inside box "+spacing+" "+spacing+" "+spacing+" "+str(xLen-2)+" "+str(yLen-2)+" "+str(zLen-2),file=packInp)
		print("end structure",file=packInp)
	packInp.close()
	#execute packmol 
	with open("in.packmol","rb") as file: 
		packmolPath = "/Users/rachelkrueger/packmol/packmol"
		output = check_output(["packmolPath"], stdin=file)
		#change above to your path to packmol!
	#open packmol's output xyz file for reading
	xyz = open('packOut.xyz','r')	
	#open what will be my program's input file
	dataFile = open(inpFile, 'w')
	line = xyz.readline(); 
	line = line.strip('\n')
	print(line,file=dataFile)
	line = xyz.readline();
	line = line.strip('\n')
	print(line,file=dataFile)
	for i in range(0,nPol): 
		for j in range(0,nBead): 
			coords = xyz.readline()
			coords = coords.strip('\n')
			print(str(i) + "  " + coords,end='  ',file=dataFile)
			if(chargeFlag==1): 
				if(j % int(1/chargeFrac) == 0): 
					print(str(chainCharge),file=dataFile)
				else: 
					print("0",file=dataFile)
			else: 
				#print("executing the non-charge case!")
				print(" 0",file=dataFile)
	#finally, open what will be program's topology file
	top = open(topFile, 'w') 
	nTot = nPol * nBead + nPlus + nMinus
	print("nTot: "+str(nTot),file=top)
	nMol = nPol + nPlus + nMinus
	print("nMol: "+str(nMol),file=top)
	print("xLen: "+str(xLen),file=top)
	print("yLen: "+str(yLen),file=top)
	print("zLen: "+str(zLen),file=top)
	print("pBeadTranslate: 0.35",file=top)
	print("pComTranslate: 0.25",file=top)
	print("pPivot: 0.15",file=top)
	print("pCrankshaft: 0.25",file=top)
	print("bonds:",file=top)
	for i in range(0,nPol):
		for j in range(0,nBead-1): 
			print(str(i)+" "+str(j)+" "+str(j+1),file=top)
	print("-1",file=top)
	print("angles:",file=top)
	print("-1",file=top)
	print("dihedrals",file=top)
	print("-1",file=top)

main()
