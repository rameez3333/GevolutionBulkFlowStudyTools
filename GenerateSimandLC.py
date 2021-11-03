import os, sys
import numpy as np
from glob import glob


JNAME = 'GevoBF_'

JINDEX = 1


stlines = open('settingsSnaptemplate.ini').readlines()
lclines = open('settingsLCtemplate.ini').readlines()

#stlines = open('testtemplatesnap.ini').readlines()
#lclines = open('testtemplatelc.ini').readlines()

#gauge = 'Newtonian'
#gauge = 'Poisson'

def FindBiggestFlowinFile(fname, rad=250., angtol=30.):
    harr = np.genfromtxt(fname, skip_header=18).transpose()
    I, X, Y, Z, VX, VY, VZ, Mvir, M200b, M200c = harr[0], harr[8], harr[9],harr[10],harr[11],harr[12],harr[13], harr[2], harr[26], harr[27]
    VTOT = np.sqrt(VX*VX + VY*VY + VZ*VZ)
    XMAX, XMIN = np.max(X), np.min(X)
    YMAX, YMIN = np.max(Y), np.min(Y)
    ZMAX, ZMIN = np.max(Z), np.min(Z)
    print(len(X))
    print(XMAX, XMIN, YMAX, YMIN, ZMAX, ZMIN)
    print(np.max(VTOT), np.min(VTOT))
    centerX, centerY, centerZ = (XMAX + XMIN)/2. , (YMAX + YMIN)/2. , (ZMAX + ZMIN)/2.
    rads = np.sqrt((X - centerX)**2. + (Y - centerY)**2. + (Z - centerZ)**2. )
    indices = np.arange(len(rads))
    cindices = indices[rads<150]
    
    def FindBiggestCorrelatedFlow(indices, rad=100., angtol = 30.):
        BVels = []
        Angs = []
        for ind in indices:
            index = int(ind)
            Xnow, Ynow, Znow = X[index], Y[index], Z[index]
            VXnow, VYnow, VZnow = VX[index], VY[index], VZ[index]
        #print Xnow, Ynow, Znow, VXnow, VYnow, VZnow
            radcut = (np.sqrt(np.power(X-Xnow, 2) + np.power(Y-Ynow, 2) + np.power(Z-Znow, 2)) < rad)
            nbulkhaloes = len(X[radcut])
            VbulkX, VbulkY, VbulkZ = np.sum(VX[radcut])/float(nbulkhaloes), np.sum(VY[radcut])/float(nbulkhaloes), np.sum(VZ[radcut])/float(nbulkhaloes)
            Speedbulk = np.sqrt(VbulkX*VbulkX + VbulkY*VbulkY + VbulkZ*VbulkZ)
            cosangHalBF = (VXnow*VbulkX + VYnow*VbulkY + VZnow*VbulkZ)/(np.sqrt(VXnow*VXnow + VYnow*VYnow + VZnow*VZnow)*Speedbulk)
            print(ind, Speedbulk, np.rad2deg(np.arccos(cosangHalBF)))
            BVels.append(Speedbulk) 
            Angs.append(np.rad2deg(np.arccos(cosangHalBF)))
        BVels = np.asarray(BVels) 
        Angs = np.asarray(Angs)
        MBF = np.max(BVels[Angs<angtol])
        index = indices[BVels==MBF][0]
        Xnow, Ynow, Znow = X[index], Y[index], Z[index]
        VXnow, VYnow, VZnow = VX[index], VY[index], VZ[index]
        radcut = (np.sqrt(np.power(X-Xnow, 2) + np.power(Y-Ynow, 2) + np.power(Z-Znow, 2)) < rad)
        nbulkhaloes = len(X[radcut])
        VbulkX, VbulkY, VbulkZ = np.sum(VX[radcut])/float(nbulkhaloes), np.sum(VY[radcut])/float(nbulkhaloes), np.sum(VZ[radcut])/float(nbulkhaloes)
        Speedbulk = np.sqrt(VbulkX*VbulkX + VbulkY*VbulkY + VbulkZ*VbulkZ)
        cosangHalBF = (VXnow*VbulkX + VYnow*VbulkY + VZnow*VbulkZ)/(np.sqrt(VXnow*VXnow + VYnow*VYnow + VZnow*VZnow)*Speedbulk)
        print('Biggest Correlated Flow', index, Speedbulk, np.rad2deg(np.arccos(cosangHalBF)), np.asarray([VXnow, VYnow, VZnow]), np.asarray([VbulkX, VbulkY, VbulkZ]))
        return index, Speedbulk, np.rad2deg(np.arccos(cosangHalBF)), np.asarray([VXnow, VYnow, VZnow]), np.asarray([VbulkX, VbulkY, VbulkZ]), np.asarray([Xnow, Ynow, Znow])
    
    return FindBiggestCorrelatedFlow(cindices, rad=rad, angtol=angtol)







def DoASet(JINDEX, theory = 'GR', gauge='Poisson', SEED=str(int(np.random.choice(range(800))))):
    TJNAME = JNAME+theory+'_'+gauge+'_'+str(JINDEX) +'.ini'
    print 'Using seed', SEED, 'Jobname', JNAME+str(JINDEX)
    outdir = theory+'_'+gauge+'_'+str(JINDEX)
    fout = open(TJNAME, "w")
    for line in stlines:
        fout.write(line.replace('__SEED__', SEED).replace('__TH__', theory).replace('__ODIR__', outdir).replace('__GGE__', gauge))
    fout.close()
    os.system('mkdir moutputs/' + outdir)
    print 'Now running gevolution just for the snap'
    os.system('mpirun -np 16 ./gevolution -n 4 -m 4 -s ' + TJNAME)
    print 'Now running rockstar on the snap'
    os.system('./rockstar -c quickstart_snapshot.cfg moutputs/' + outdir+ '/lcdm_snap000_cdm ')
    os.system('mkdir moutputs/' + outdir+'/snap')
    os.system('mv halos_0.0.* moutputs/' + outdir+'/snap/')
    ind, SB, halbang, Velhal, Velbulk, BV = FindBiggestFlowinFile('moutputs/' + outdir+'/snap/halos_0.0.ascii')
    
    print 'Now summarizing bulk flow'
    
    
    TJNameOld = TJNAME
    
    fout=open('moutputs/' + outdir+'/snap/BulkFlowSummary.txt', 'w')
    fout.write(str(ind) + ' | ' +str(SB) + ' | ' + str(halbang) +'\n')
    fout.write(str(Velhal[0]) + ' | ' +str(Velhal[1]) + ' | ' + str(Velhal[2]) +'\n')
    fout.write(str(Velbulk[0]) + ' | ' +str(Velbulk[1]) + ' | ' + str(Velbulk[2]) +'\n')
    fout.close()
    
    VH = Velhal/np.sqrt(np.sum(Velhal*Velhal))*3
    
    LCX, LCY, LCZ = str(BV[0]), str(BV[1]), str(BV[2])
    VHX, VHY, VHZ = str(VH[0]), str(VH[1]), str(VH[2])
    
    TJNAME = JNAME+theory+'_'+ gauge +'_LC'+str(JINDEX) +'.ini'
    
    #TJNAME = JNAME+theory+'_LCtest'+str(JINDEX) +'.ini'
    
    print 'Generating lightcone at the vertex', LCX, LCY, LCZ, 'in the direction', VHX, VHY, VHZ
    
    fout = open(TJNAME, "w")
    for line in lclines:
        fout.write(line.replace('__SEED__', SEED).replace('__TH__', theory).replace('__ODIR__', outdir).replace('_X_', LCX).replace('_Y_', LCY).replace('_Z_', LCZ).replace('_DX_', VHX).replace('_DY_', VHY).replace('_DZ_', VHZ).replace('__GGE__', gauge) )
    fout.close()
    
    os.system('rm moutputs/' + outdir+'/*')
    
    print 'Now running gevolution again for the lightcone'
    os.system('mpirun -np 16 ./gevolution -n 4 -m 4 -s ' + TJNAME)
    
    print 'Now running lccat and lcmap'
    
    lcflist = glob('moutputs/' + outdir+'/lcdm_lightcone_*.map' )
    ind = [int(j.split('_')[-2]) for j in lcflist]
    
    
    print 'Now doing lcmap'
    
    print ' ./lcmap -s '+ TJNAME + ' -d 0.5 ' 
    
    
    #raw_input('Lets see whats going on ')
    
    os.system(' ./lcmap -s '+ TJNAME + ' -d 0.5 ' )
    
    
    
    cs = str(np.min(ind)) + '-' + str(np.max(ind))
    
    print 'Now doing lccat'
    
    print ' ./lccat -s '+ TJNAME + ' -c ' + cs
    
    os.system(' ./lccat -s '+ TJNAME + ' -c ' + cs)
    
    
    print 'Finished processing lccat'
    

    
    print 'Now running rockstar on the LC'
    os.system('./rockstar -c quickstart_lc.cfg moutputs/' + outdir+ '/lcdm_lightcone_cdm ')
    
    os.system('mkdir moutputs/' + outdir+'/lightcone')
    os.system('mv halos_0.0.* moutputs/' + outdir+'/lightcone/')
    
    print 'Now moving the Job files'
    os.system('mv ' + TJNameOld + ' moutputs/Jobs/')
    os.system('mv ' + TJNAME + ' moutputs/Jobs/')
    
    

#def DoASetBatch()


   
gauges = ['Poisson']#, 'Newtonian']
theories = ['GR', 'Newton']    
   
for i in range(42, 50):
    seed = str(int(np.random.choice(range(1300, 1800))))
    for t in theories:
        for g in gauges:
            DoASet(i, theory=t, gauge=g, SEED=seed)
