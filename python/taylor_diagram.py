import matplotlib.pyplot as plt
import math
import numpy as np
import scipy.stats

def calc_error_triangle(varobs, varfld, rcorr, sqrerr):

    if ( sqrerr ): calc_sqrerr=sqrerr
    if ( rcorr ): calc_rcorr=rcorr
    if ( not sqrerr ):
        calc_sqrerr = varobs + varfld - 2*math.sqrt(varobs * varfld)*rcorr
    if ( not rcorr ):
        calc_rcorr = ( varobs + varfld - sqrerr ) / 2.0 / math.sqrt(varobs*varfld)
    if ( sqrerr ): calc_sqrerr=sqrerr
    if ( rcorr ): calc_rcorr=rcorr
    # IS IT POSSIBLE TO CALC var_obs/var_fld if unknown?   
    residual =  varobs + varfld - 2*math.sqrt(varobs * varfld)*calc_rcorr - calc_sqrerr
    return calc_rcorr, calc_sqrerr, residual

def addpoint_taylor_diagram( varobs, varfld, rcorr, sqrerr, axe, color='k', label='', align=['left', 'bottom'], addpt=True):
    check_rcorr, check_sqrerr, residual = calc_error_triangle(varobs, varfld, rcorr, sqrerr)
    print('RESIDUAL', residual)
    if ( not sqrerr ): sqrerr = check_sqrerr
    if ( not rcorr ): rcorr = check_rcorr
    
    norm_varfld = varfld / varobs
    stdfld = math.sqrt(norm_varfld)
    
    cosphi, sinphi = cossin(rcorr)
    
    Xpt = stdfld * cosphi
    Ypt = stdfld * sinphi

    if ( addpt ):  
        axe.plot( Xpt, Ypt, marker='o', markerfacecolor=color, markersize=10, fillstyle='full')
        axe.text( Xpt, Ypt, label, horizontalalignment=align[0], verticalalignment=align[1])
    return Xpt, Ypt

def cossin(rcorr):
    cosphi = rcorr
    sinphi = math.sqrt(1.0-rcorr**2)
    return cosphi, sinphi    

def create_taylor(maxlim=2.0):
    xrange = [0, 2.0]
    yrange = [0, 2.0]
    if ( isinstance(maxlim, float) or isinstance(maxlim, int) ): 
       xrange[1] = maxlim
       yrange[1] = maxlim
    elif ( len(maxlim) == 2 ):
       xrange=maxlim
       yrange=maxlim
    elif( len(maxlim) == 4 ):
       xrange=maxlim[0:2]
       yrange=maxlim[2:4]
    
    fig, axe = plt.subplots()
    axe.set_xlim(xrange)
    axe.set_ylim(yrange)
    for radius in [0.25, 0.5, 1.0, 1.5]:
        circle0 = plt.Circle((0,0), radius, color='black', fill=False)
        circle1 = plt.Circle((1,0), radius, color='green', fill=False)
        axe.add_patch(circle0)
        axe.add_patch(circle1)
    return fig, axe
    
def add_spokes(axe):
    for rcorr in np.arange(0.0, 1.0, 0.1):
        radius=np.arange(0.0, 2.0, 0.1)
        cosphi, sinphi = cossin(rcorr)
        axe.plot( radius*cosphi, radius*sinphi, color='cyan' )
        
def make_taylor_figure( list_of_points, list_of_labels, list_of_colors=None, pltfile='taylor', maxlim=2.0, addvalues=False ):
    align_default = ['left', 'bottom']
    # EACH LIST ELEMENT ( varobs, varfld, rcorr, sqrerr)
    if ( maxlim == 'auto' ):
        fig, axe = create_taylor()
    else:
        fig, axe = create_taylor(maxlim=maxlim)
        align = align_default
    add_spokes(axe)

    if ( not list_of_colors ):
        list_of_colors = ['r', 'b', 'g', 'c', 'm' ]
    ncolors = len(list_of_colors) 
    minX = 100
    maxX = -100
    minY = 100
    maxY = -100  
    XYPTS=[]  
    for ipoint, points in enumerate(list_of_points):
        varobs, varfld, rcorr, sqrerr = points
        Xpt, Ypt = addpoint_taylor_diagram( varobs, varfld, rcorr, sqrerr, axe, addpt=False)
        XYPTS.append((Xpt, Ypt))
        if ( np.isfinite(Xpt) ):
          minX = np.min([minX, Xpt])
          maxX = np.max([maxX, Xpt])
        if ( np.isfinite(Ypt) ):
          minY = np.min([minY, Ypt])
          maxY = np.max([maxY, Ypt])
    midX = maxX - minX
    midY = maxY - minY

    for ipoint, points in enumerate(list_of_points):
        Xpt, Ypt = XYPTS[ipoint]
        clr = list_of_colors[ipoint%ncolors]
        label=list_of_labels[ipoint]
        varobs, varfld, rcorr, sqrerr = points
        align = align_default
        if ( maxlim == 'auto' ):
            if ( Xpt > midX ): align = ['right', align_default[1]]
        if ( addvalues ):
           rcor_str = "%.4f" % rcorr
           rmse_str = "%.4f" % math.sqrt(sqrerr)
           label=label+' '+'rcorr='+rcor_str+'/rmse='+rmse_str
        Xpt, Ypt = addpoint_taylor_diagram( varobs, varfld, rcorr, sqrerr, axe, color=clr, label=label, align=align)
    if ( maxlim == 'auto' ):
        axe.set_xlim([minX-(maxX-minX)/10.0, maxX+(maxX-minX)/10.0])
        axe.set_ylim([minY-(maxY-minY)/10.0, maxY+(maxY-minY)/10.0])
    fig.savefig(pltfile+'.pdf')
    fig.savefig(pltfile+'.png')
    plt.close(fig)
    return
    
def generate_data():
    N=100
    r = np.array([ [1, 0.9, 0.7, 0.5], [0.9, 1, 0.3, 0.1], [0.7, 0.3, 1, -0.1 ], [0.5, 0.1, -0.1, 1] ] )
    rng = np.random.default_rng()
    Y = rng.multivariate_normal((0,0,0,0), r, size=N)
    X0 = Y[:,0]
    X1 = Y[:,1]
    X2 = Y[:,2]
    X3 = Y[:,3]
    return X0, X1, X2, X3

def variance(X):
    XB = np.mean(X)
    VR = np.mean( np.square(X-XB) )
    return VR

def sqrerr(X, Y):
    sqrerr = np.mean( np.square(X-Y) )
    return sqrerr

def generate_points():
    X0, X1, X2, X3 = generate_data()
    V0 = variance(X0)
    V1 = variance(X1)
    V2 = variance(X2)
    V3 = variance(X3)
    
    E1 = sqrerr(X0, X1)
    E2 = sqrerr(X0, X2)
    E3 = sqrerr(X0, X3)
    
    r1, __, R = calc_error_triangle(V0, V1, None, E1); print('1', r1, scipy.stats.pearsonr(X0, X1)[0], R)
    r2, __, R = calc_error_triangle(V0, V2, None, E2); print('2', r2, scipy.stats.pearsonr(X0, X2)[0], R)
    r3, __, R = calc_error_triangle(V0, V3, None, E3); print('3', r3, scipy.stats.pearsonr(X0, X3)[0], R)

    return (V0, V1, r1, E1), (V0, V2, r2, E2), (V0, V3, r3, E3)
    
    
    
