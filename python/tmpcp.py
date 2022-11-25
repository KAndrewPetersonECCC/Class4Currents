import os
import subprocess
import random
        
def tmpcpold(file, tmpdir='/space/hall1/sitestore/eccc/mrd/rpnenv/dpe000/tmpdir'):
     subprocess.call([ 'mkdir', '-p', tmpdir])
     subprocess.call([ 'scp', file, tmpdir+'/'])
     bile=os.path.basename(file)
     nile=tmpdir + '/' + bile
     return nile 

def tmprsync(file, tmpdir='/space/hall1/sitestore/eccc/mrd/rpnenv/dpe000/tmpdir'):
     subprocess.call([ 'mkdir', '-p', tmpdir])
     subprocess.call([ 'rsync', '-avz', file, tmpdir+'/'])
     bile=os.path.basename(file)
     nile=tmpdir + '/' + bile
     return nile 

def tmpcp(file, tmpdir='/space/hall1/sitestore/eccc/mrd/rpnenv/dpe000/tmpdir',scp='sscp'):
     subprocess.call([ 'mkdir', '-p', tmpdir])
     subprocess.call([ scp, file, tmpdir+'/'])
     bile=os.path.basename(file)
     nile=tmpdir + '/' + bile
     return nile 

def gunzip(file):    


    if ( file[-3:] == '.gz' ):
        nile = file[:-3]
        bile=os.path.basename(nile)
        ddir=os.path.dirname(nile)
        tile=ddir+'/tmp'+str(random.randint(0,100)).zfill(2)+'.'+bile+'.gz'
        if ( os.path.isfile(nile) ):
            print('file exists already')
            if ( os.path.isfile(file) ):
                subprocess.call([ 'mv', file, tile] )
                subprocess.call([ 'gunzip', tile])
                tile=tile[:-3]
                if ( os.path.isfile(tile) ):
                    rc = subprocess.call([ 'diff', tile, nile])
                    if ( rc == 0 ):
                        print('Removing duplicate', tile)
                        subprocess.call([ 'rm', tile])
                    elif ( os.path.isfile(tile) ):
                        print('Clobbering existing', nile)
                        subprocess.call([ 'mv', tile, nile])
        else:
            if ( os.path.isfile(file) ):
                subprocess.call([ 'gunzip', file])
    else:
        nile = file
        print('WARNING '+file+' does not end in .gz :: NO CHANGE')
         
    return nile

def tmpcpgz(file, tmpdir='/space/hall1/sitestore/eccc/mrd/rpnenv/dpe000/tmpdir'):
    nile = tmpcp(file, tmpdir=tmpdir)
    nnile = gunzip(nile)
    return nnile
       
def tmprm(file):
    subprocess.call([ 'rm', file ])
    return     

def cmcarc(file, arguements=["-v", "-t"], rm=False):
    command=['cmcarc']
    command.extend(arguements)
    command.extend(["-f", file])
    print(command)
    rc=subprocess.call(command)
    if ( rm ):
        subprocess.call(['rm', '-f', file])
    return rc

def cdcmcarc(file, dir, arguements=["-v", "-x"],rm=False):
    owd=os.getcwd()
    os.chdir(dir)
    cwd=os.getcwd()
    rc=cmcarc(file, arguements=arguements, rm=rm)
    os.chdir(owd)
    rwd=os.getcwd()
    #print(owd, cwd, rwd)
    #print(os.listdir(cwd))
    return rc
