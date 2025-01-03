import numpy as np

#NOTE:  I DO NOT KNOW THE PROPER WAY TO DEAL WITH BOUNDARIES YET.
#       WILL APPLY STENCIL ONLY OVER UNMASKED POINTS.
#       Normalization is by the sum of the unmasked point stencils.  

F=[]
F.append([2,1])                                        # p=0
F.append([10, 4, -1])                                # p=1
F.append([44, 15, -6, 1])                        #p=2
F.append([186, 56, -28, 8, -1])                        #p=3
F.append([772, 210, -120, 45, -10, 1])                #p=4
F.append([3172, 792, -495, 220, -66, 12, -1])        #p=5
F.append([12952, 3003, -2002, 1001, -364, 91, -14, 1])                                #p=6
F.append([52666, 11440, -8008, 4368, -1820, 560, -120, 16, -1])                        #p=7
F.append([213524, 43758, -31824, 18564, -8568, 3060, -816, 153, -18, 1])        #p=8
F.append([863820, 167960, -125970, 77520, -38760, 15504, -4845, 1140, -190, 20, -1])        #p=9

def shapiro(AX, p, npass=1, cyclic=True):    
    # AX is a masked? 1-D array
    # p is order of filter

    if ( not isinstance(p, int) ):
        raise Exception('p must be integer '+str(p) )
    if ( p > 9 ): 
        raise Exception('p must < 10, '+str(p) )
    f = F[p]
    
    if ( npass > 1):
        BX = AX.copy()
        for ipass in range(npass):
            BX=shapiro(BX,p)
    else:
        if ( AX.ndim == 2 ):
            ny, nx = AX.shape
            BX = AX.copy()
            if ( isinstance(cyclic, list) ):
                xcyclic=cyclic[0]
                ycyclic=cyclic[1]
            else:
                xcyclic=cyclic
                ycyclic=cyclic
            #print 'X direction'
            for iy in range(ny):
                BX[iy, :] = shapiro(AX[iy,:],p, cyclic=xcyclic)
            #print 'Y direction'
            for ix in range(nx):
                BX[:,ix] = shapiro(BX[:,ix], p, cyclic=ycyclic)
        elif ( AX.ndim > 2 ):
            BX = []
        else:        
            ni = len(AX)
            if ( cyclic ):    
                BD = np.zeros(ni).astype('float')
                if ( isinstance(AX, np.ma.core.MaskedArray) ):
                    AD = AX.data
                    AM = 1-AX.mask.astype(int)   # Returns masked (T) values as 0.   
                else:
                    AD = AX[:]
                    AM = np.ones(ni)

                for ii in range(ni):
                    Btmp = float(AD[ii]) * f[0] * AM[ii]
                    Ntmp = f[0] * AM[ii]
                    for ip in range(1,p+2):
                        # ii-ip
                        Btmp = Btmp + float(AD[(ii-ip)%ni])*f[ip]*AM[(ii-ip)%ni]
                        Ntmp = Ntmp+f[ip]*AM[(ii-ip)%ni]
                        # ii+ip
                        Btmp = Btmp + float(AD[(ii+ip)%ni])*f[ip]*AM[(ii+ip)%ni]
                        Ntmp = Ntmp + f[ip]*AM[(ii+ip)%ni]
                
                    #print 'ii', ii, Btmp, Ntmp
                    if ( Ntmp == 0 ):
                        BD[ii] = 0
                    else:
                        BD[ii] = AM[ii]*Btmp / Ntmp                
                if ( isinstance(AX, np.ma.core.MaskedArray) ):
                    BX = np.ma.masked_array(BD, mask=AX.mask)
                else:
                    BX = BD[:]
            else: # not cyclic
                # Add ghostzone and do again
                
                BX = AX.copy()
                if ( not isinstance(AX, np.ma.core.MaskedArray) ): 
                    BX = np.ma.masked_array(BX, mask=False)
                    BUF = np.ma.masked_array(np.zeros(p+1), mask=True)
                BX = np.ma.masked_array( np.append(BX.data, BUF.data), np.append(BX.mask, BUF.mask) )
                BX = shapiro(BX, p, cyclic=True)
                BX = BX[:ni]
                if ( not isinstance(AX, np.ma.core.MaskedArray) ): 
                    BX = BX.data                

    return BX

def shapiro2D(FLD,npass=1, wtc=None):
    #  !DESCRIPTION:  Smooth a 2D field
    #
    #  FLD is smoothed by niter applications of a 
    #  9-point averaging stencil:
    #
    #           1/16 alpha -- 2/16 alpha -- 1/16 alpha
    #           |                  |                 |
    #           2/16 alpha --     wtc    -- 2/16 alpha
    #           |                  |                 |
    #           1/16 alpha -- 2/16 alpha -- 1/16 alpha
    #
    #  Elements of the stencil are multiplied by the local cell thickness h_i
    #  and wtc = 1 - [ sum over eight surrounding cells (h_i * wt_i)]/h_c, where
    #  h_c is the thickness of the center cell.
    #  Land cells are given zero thickness so they are treated the same as other 
    #  points but given no weight in the averaging proceedure.
    #  Note that only the center weight is modified so that exchanges between
    #  cells remain balanced.  The normal practice of changing all weights
    #  when one or more cells are land points can result in different weights being 
    #  used when the cell is giving (i.e. on the boundaries of the stencil) from 
    #  when it is receiving (i.e., at the center of the cell) and this can result
    #  in spatial smoothing changing the net heat and salt sources and sinks. 
    #  The present approach avoids this problem.
    #
    #  For a function of the form cos(kx)*cos(ly), a single pass of this 
    #  spatial smoother in a region with uniform cell thicknesses will have 
    #  a gain of approximately 1 - 0.25*alpha*((k dx)**2 + (l dy)**2).
    #  For a feature with lambda_x/dx = lambda_y/dy = 2*n (k*dx=l*dy=pi/n) 
    #  the gain is 1 - 0.5*alpha*(pi/n)**2  --  for a single pass of the filter.
    #
    #  For the recursive filter applied to the mean (omega_c = 0) climatological
    #  conditions, the complex gain for variations at frequency omega is approximately
    #  1/[1 + 0.5*(alpha/Kappa)*(pi/n)**2 + i * omega*dt/(2*Kappa)].
    #  If alpha = M_alpha * Kappa, the complex gain is
    #  1/[1 + 0.5*M_alpha*(pi/n)**2 + (-1)^0.5 * omega*dt/(2*Kappa)].
    #  At very low frequencies, the gain reduces to 1/[1+0.5*M_alpha*(pi/n)**2].
    #  The reduction in this limit is due entirely to the spatial filter.
    #  Transient eddies with time scales short compared to 1/Kappa are damped 
    #  much less since the filtered nudge has small amplitude.
    #  The reduced gain in this case is due to the time filter.
    #
    #  For a feature with a diameter of 10*dx (n=10), and using M_alpha = 100, 
    #  the gain at  omega = 0 is approximately 1/6. With a diameter of 5*dx 
    #  the gain with M_alpha = 100 reduces to approximately 1/21.
    #  
    #  Fortran version: /fs/ssm/eccc/mrd/rpn/OCEAN/master_4.0/concepts-tools_5.0.2-intel-2016.1.156_ubuntu-14.04-amd64-64/src/src/filter/nemo_filter.F90
    #
    #-----------------------------------------------------------------------
    #

    if ( FLD.ndim < 2 ):
        raise Exception('Requires at least 2D field, nd = '+str(FLD.dim) )
    elif ( FLD.ndim > 2 ):
        # Recursively Apply Filter
        nz, nx, ny = FLD.shape
        NFLD=FLD
        for kk in range(nz):
            NFLD[kk, :, :], __ = shapiro2D(FLD[kk,:,:], npass=npass)
    else:
        if ( npass > 1 ):
            NFLD = FLD
            for ipass in range(npass):
                NFLD, wtc = shapiro2D(NFLD, wtc=wtc)
        else:  # Do the 1-pass
            # Depending on where field is coming from x/y dimensions may be flipped.
            # Choose conventions for standard file (but no special x or y treatment).    
            nx, ny = FLD.shape
        
            # Calulate Central Weight (16 - sum surrounding weights)
        
            if ( isinstance(FLD, np.ma.core.MaskedArray) ):
                FD = FLD.data
                FM = 1-FLD.mask.astype(int)   # Returns masked (T) values as 0.   
            else:
                FD = FLD[:,:]
                FM = np.ones( (nx, ny) )

            # Calulate Central Weight (16 - sum surrounding weights)
            # Coded to pass wtc to next pass.
            if ( isinstance(wtc, type(None)) ):
                wtc = np.zeros((nx,ny))
                for jj in range(1,ny-1):  
                    for ii in range(1,nx-1):
                        wtctemp = (  
                                     1*FM[ii-1,jj+1] + 2*FM[ii,jj+1] + 1*FM[ii+1,jj+1]  
                                   + 2*FM[ii-1,jj]   + 0             + 2*FM[ii+1,jj] 
                                   + 1*FM[ii-1,jj-1] + 2*FM[ii,jj-1] + 1*FM[ii+1,jj-1]
                                   )
                        wtc[ii,jj] = 16 - wtctemp
                            
            FFT = np.zeros((nx,ny))
            for jj in range(1,ny-1):  
                for ii in range(1,nx-1):
                    FIJ = float(  
                                 1*FD[ii-1,jj+1]*FM[ii-1,jj+1] +          2*FD[ii,jj+1]*FM[ii,jj+1] + 1*FD[ii+1,jj+1]*FM[ii+1,jj+1] 
                               + 2*FD[ii-1,jj]  *FM[ii-1,jj]   + wtc[ii,jj]*FD[ii,jj]  *FM[ii,jj]   + 2*FD[ii+1,jj]  *FM[ii+1,jj]
                               + 1*FD[ii-1,jj-1]*FM[ii-1,jj-1] +          2*FD[ii,jj-1]*FM[ii,jj-1] + 1*FD[ii+1,jj-1]*FM[ii+1,jj-1]
                               ) / 16.0 
                    FFT[ii,jj] = FIJ


            #print 'MAX/MIN', np.min(wtc[1:-1,1:-1]), np.max(wtc[1:-1,1:-1])
            #print 'max/min', np.min(FFT), np.max(FFT)
            #print 'max/min', np.min(FD), np.max(FD)
            #print 'max/min', np.min(FLD), np.max(FLD)                
            # APPLY THE NORTH FOLD CONDITION (SHOULD BE MODULE)
            FFT[1:,-1] = FFT[-1:0:-1,-3]
            # APPLY CYCLIC WITH BUFFER ZONES in X-direction (SHOULD BE MODULE)
            FFT[0, :] = FFT[nx-2,:]
            FFT[nx-1,:] = FFT[1,:]

            if ( isinstance(FLD, np.ma.core.MaskedArray) ):
                NFLD = np.ma.masked_array( FFT, mask=FLD.mask )
            else:
                NFLD = FFT
            #print 'max/min', np.min(NFLD), np.max(NFLD)                

    return NFLD, wtc

def SHAPIRO_RETURN(L, N, n=1, dx=1):
    SINE = np.sin( np.pi * dx / L )
    R = ( 1 - np.square(SINE)**(n) )**N
    return R

def LENGTH(max, iter=1):
    LENGTH = np.arange(2*iter, max, iter)
    return LENGTH

def SHAPIRO_CUTOFF(L, N, n=1, dx=1, cut=0.01):
    R =SHAPIRO_RETURN(L, N, n=n, dx=dx)
    PASS = np.where( ( R > cut ) & ( L > 2*dx ) )
    CUTT = np.where( R < cut )
    MIN_KEPTINN = np.min(L[PASS[0]])
    MAX_REMOVED = np.max(L[CUTT[-1]])
    return MAX_REMOVED, MIN_KEPTINN
