def displayPotentialFromUg(folder=None, potSize=([256,256]), Ncells=([1,1]), showImag=False, doFlip=False, thresh = None, beta = None, InputUg=None, InputBeams = None, returnComplexUg=None, scaling=1,scaleImag=0):
    
	import os
	# import pyQEDA as pq #note this python source file/jupyter notebook should not be called 'pyQEDA.py' as this is reserved for the pyQEDA module as will cause a conflict if named so
	import numpy as np
	import math
	import matplotlib.pyplot as plt
	from copy import deepcopy
	currentPath = os.getcwd()
	

	def displayPotential_2D(UgList=None,InputBeams=None,potSize=([256,256]),Ncells=([1,1]),showImag=True,doFlip=False, thresh=thresh, beta = beta):
		# This function displays the real-space potential corrsponding to a set of
		# structure factors.  
		# Input: Ncells = number of unit cells
		# potSize: size of potential map in pixels
    
		NxMid = np.int(np.floor(potSize[1]/2))
		NyMid = np.int(np.floor(potSize[0]/2))
		potMap = np.zeros(potSize,dtype=complex)

		for j in range(len(UgList)):
			print(UgList[j])
			potMap[NyMid+InputBeams[j,0]*Ncells[1],NxMid+InputBeams[j,1]*Ncells[0]] = UgList[j]
		potMap = np.fft.ifft2(np.fft.ifftshift(potMap))

		if showImag:
			plt.imshow(potMap.real)
			plt.title(f"Potential Map Real - {folder}\n\n")
			plt.colorbar()
			plt.show()
			plt.imshow(potMap.imag)
			plt.title(f"Potential Map Imaginary - {folder}\n\n")
			plt.colorbar()
			plt.show()
			plt.imshow(np.abs(potMap))
			plt.title(f"Potential Map Abs - {folder}\n\n")
			plt.colorbar()
			plt.show()

		# Do charge flipping
        # Since charge/potential flipping only makes sense when we have a purely real potential, we will do two things:
        # a) drive the imaginary part of the potential towards zero by multiplying it with 0.5*beta
        # b) apply charge/potential flipping to the real part by multiplying those parts of it that are below 
        #    threshold by -beta. 
		if doFlip:      
			potMap = potMap.real+1j*scaleImag*beta*potMap.imag # this reduces the imaginary part
			tthresh = thresh*np.max(np.abs(potMap))
			ind = np.nonzero(potMap.real < tthresh)
            # the following expression flips the real part and keeps the imaginary part.
			potMap[ ind ] =  -beta*potMap[ind].real + 1j*potMap[ind].imag 

			#print(potMap[ind].shape)
			if showImag:
				plt.imshow(potMap.real)
				plt.title(f"Potential Map Charge Flipped - Real - {folder}\n\n")
				plt.colorbar()
				plt.show()
				plt.imshow(potMap.imag)
				plt.title(f"Potential Map Charge Flipped - Imaginary - {folder}\n\n")
				plt.colorbar()
				plt.show()
				plt.imshow(np.abs(potMap))
				plt.title(f"Potential Map Charge Flipped Abs - {folder}\n\n")
				plt.colorbar()
				plt.show()

			pmf = np.fft.fftshift(np.fft.fft2(potMap))	

			for j in range(len(UgList)):
				UgList[j] = pmf[NyMid+InputBeams[j,0]*Ncells[1],NxMid+InputBeams[j,1]*Ncells[0]]

			max_val = np.max(UgList)
			max_idx = np.argmax(UgList)
			print(f'max Ug {max_idx, max_val}')
		return UgList, potMap



	def displayPotentialFromUgFile(filepath=folder,potSize=potSize,Ncells=Ncells,showImag=showImag,doFlip=doFlip, thresh=thresh, beta=beta):
		beamfilename = '_BeamList.txt'
		beamfile = filepath + beamfilename
		print(beamfile)
		# UgReal = np.loadtxt(beamfile,skiprows=1, usecols= range(4, 5)  )
		# UgImag = np.loadtxt(beamfile,skiprows=1, usecols= range(5, 6)  )
		UgFileName = '_UgMasterList.txt'
		UgFile = filepath + UgFileName
		Ugimport = np.loadtxt(UgFile,skiprows=1,usecols= range(1, 3))
		Ug = np.array(range(int(Ugimport.size/2)), dtype=complex)
		Ug.real = Ugimport[:,0]
		Ug.imag = Ugimport[:,1]
		beams = np.asarray(np.loadtxt(beamfile,skiprows=1, usecols= range(1, 4)  ), dtype=int ) 
		Ugcopy = deepcopy(Ug)	
		UgList, potMap = displayPotential_2D(UgList=Ugcopy,InputBeams=beams,potSize=potSize,Ncells=Ncells,showImag=showImag,doFlip=doFlip, thresh=thresh, beta=beta)		
		if doFlip:
			print(f'\nDifference pre and post Ug charge flipping\n{Ugcopy[0]} -> {UgList[0]} ')

		
		return UgList, potMap, beams

	def displayPotentialFromUgVar(InputUg=InputUg, InputBeams=None, potSize=potSize,Ncells=Ncells,showImag=showImag,doFlip=doFlip, thresh=thresh, beta=beta):
		Ugcopy = deepcopy(InputUg)	
		UgList, potMap = displayPotential_2D(InputUg,InputBeams=InputBeams,potSize=potSize,Ncells=Ncells,showImag=showImag,doFlip=doFlip, thresh=thresh, beta=beta)		
		if doFlip:
			print(f'\nDifference pre and post Ug charge flipping\n{Ugcopy[0]} -> {UgList[0]} ')
			print(f'\nMax Difference =\n{np.max(UgList-Ugcopy)}')
		return UgList, potMap

	
	#MAIN #if a folder is selected then the Ug are read driectly from this folder. Otherwise, they are read from the input array.
	if folder!=None:
		filepath = "/Users/chenjie/pyQEDA/pyQEDA_simulations/" + folder + "/" + folder + '_experimental_data' +'/'
		[UgList, potMap, beams] = displayPotentialFromUgFile(filepath)
		if returnComplexUg:
			return UgList, potMap
		else:
			return [(np.ravel([UgList.real,UgList.imag],'F')), potMap, beams]
	else:
		print(f'\n\nInputUg into Potential = {InputUg}')
		UgComplex = np.zeros(np.int( (InputUg.size -1)/2),dtype=complex)
		UgComplex.real = InputUg[0:-1:2] # excludes last element which is thickness
		UgComplex.imag = InputUg[1:-1:2] 
		[UgList, potMap] = displayPotentialFromUgVar(UgComplex / scaling , InputBeams)
		
		
		if returnComplexUg:
			return UgList, potMap, InputBeams
		else:
			UgList = UgList * scaling
			UgList = (np.ravel([UgList.real,UgList.imag] ,'F'))
			UgList = np.append([UgList],[InputUg[-1]])
			return UgList, potMap, InputBeams

	
