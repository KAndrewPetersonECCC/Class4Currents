import numpy as np
import datetime
import os.path
import subprocess

def convert_datelist_strint(dates):
    if ( isinstance(dates, datetime.datetime) ):
       date_str, date_int = convert_date_strint(dates) 
       return date_str, date_int
    dates_str = []
    dates_int = []
    for date in dates:
        date_str, date_int = convert_date_strint(date)
        dates_str.append(date_str)
        dates_int.append(date_int)
    return dates_str, dates_int

def convert_date_strint(date):
    if ( isinstance(date, list) ):
       dates_str,dates_int = convert_datelist_strint(date)
       return dates_str, date_int
    if ( not isinstance(date, datetime.datetime) ):
       return None, None
    date_str=date.strftime("%Y%m%d%H")
    date_int=int(date_str)
    return date_str, date_int
    
def convert_strint_datelist(dates_str):
    if ( ( isinstance(dates_str, str) ) or ( isinstance(dates_str, int) ) or ( isinstance(dates_str, np.int64) ) ):
       date = convert_strint_date(dates_str) 
       return date
    dates = []
    for date_str in dates_str:
        date = convert_strint_date(date_str)
        dates.append(date)
    return dates

def convert_strint_date(date_str):
    if ( isinstance(date_str, int) ):
        date_str=str(date_str)
    if ( isinstance(date_str, np.int64) ):
        date_str=str(date_str)
    if ( isinstance(date_str, str) and ( len(date_str) == 8 ) ):
        date=datetime.datetime.strptime(date_str, '%Y%m%d')  
    elif ( isinstance(date_str, str) and ( len(date_str) == 10 ) ):
        date=datetime.datetime.strptime(date_str, '%Y%m%d%H')  
    elif ( isinstance(date_str, str) and ( len(date_str) == 12 ) ):
        date=datetime.datetime.strptime(date_str, '%Y%m%d%H%M')  
    elif ( isinstance(date_str, str) and ( len(date_str) == 14 ) ):
        date=datetime.datetime.strptime(date_str, '%Y%m%d%H%M%S')  
    elif ( isinstance(date_str, str) and ( len(date_str) == 20 ) ):
        date=datetime.datetime.strptime(date_str, '%Y%m%d%H%M%S%f')  
    else:
        date=None
    return date
    
def read_file(file='datafile.dat'):
    if ( not os.path.isfile(file) ):
        print('File does not exist, ', file)
        print('Create file')
        with open(file, mode='a'): pass
        return [], []
    try:
        alldata = np.loadtxt(file, unpack=True) 
        try: 
            nv, nt = alldata.shape
        except:
            nv = alldata.size
            nt = 1
            alldata = np.reshape(alldata, (nv,nt))
        dates = alldata[0,:]
        datas = (alldata[1:,:])
        dates = dates.astype(int)
    except:
        datas = np.array([])
        dates = np.array([])
    return dates, datas
    
def read_date_and_boolean(file='datafile.dat'):
    if ( not os.path.isfile(file) ):
        print('File does not exist, ', file)
        return [], []
    try:
        alldata = np.loadtxt(file, unpack=True) 
        try: 
            nv, nt = alldata.shape
        except:
            nv = alldata.size
            nt = 1
            alldata = np.reshape(alldata, (nv,nt))
        dates = alldata[0,:]
        boole = (alldata[1:,:])
        dates = dates.astype(int)
        boole = boole.astype(bool)
    except:
        boole = np.array([])
        dates = np.array([])
    return dates, boole
    
    
def write_file(dates, datas, file='datafile.dat'):
    with open(file, 'w') as f:
        nv, nt = datas.shape
        for idate in range(nt):
            if ( len(str(dates[idate])) == 8 ):
                f.write("%8d" % dates[idate])
            if ( len(str(dates[idate])) == 10 ):
                f.write("%10d" % dates[idate])
            if ( len(str(dates[idate])) == 12 ):
                f.write("%12d" % dates[idate])
            if ( len(str(dates[idate])) == 14 ):
                f.write("%14d" % dates[idate])
            for iv in range(nv):
                f.write(" %g" % datas[iv, idate])
            f.write("\r\n")
    return     
    
def write_date_and_boolean(dates, boole, file='datafile.dat'):
    with open(file, 'w') as f:
        nv, nt = boole.shape
        for idate in range(nt):
            if ( len(str(dates[idate])) == 8 ):
                f.write("%8d" % dates[idate])
            if ( len(str(dates[idate])) == 10 ):
                f.write("%10d" % dates[idate])
            if ( len(str(dates[idate])) == 12 ):
                f.write("%12d" % dates[idate])
            if ( len(str(dates[idate])) == 14 ):
                f.write("%14d" % dates[idate])
            for iv in range(nv):
                f.write(" %1d" % boole[iv, idate].astype(int))
            f.write("\r\n")
    return     
    
def write_2dates_file(dates1, dates2, datas, file='datafile.dat'):
    f = open(file, 'w')
    if ( datas.ndim == 1 ):
        nt=1
        nv=len(datas)
        datas=np.reshape(datas, (nv, nt))
    else:
        nv, nt = datas.shape
    for idate in range(nt):
        f.write("%10d  %10d" % ( dates1[idate], dates2[idate] ) )
        for iv in range(nv):
            f.write(" %g" % datas[iv, idate])
        f.write("\r\n")
    f.close()
    return     
    
def write_integer_to_file(n, file='datafile.dat'):
    with open(file, 'w') as f:
        f.write("%3d" % n)
        f.write("\r\n")
    return     

def read_integer_from_file(file='datafile.dat'):
    if ( not os.path.isfile(file) ):
        print('File does not exist, ', file)
        return []
    f = open(file, 'r')
    n=int(f.read(3))
    f.close()
    return n    

def write_intlist_to_file(int_list, file='datafile.dat'):
    f = open(file, 'w')
    for n in int_list:
        f.write("%6d" % n)
    f.write("\r\n")
    f.close()
    return     

def read_intlist_from_file(n_len, file='datafile.dat'):
    if ( not os.path.isfile(file) ):
        print('File does not exist, ', file)
        int_list=[]
        for i in range(n_len):
            int_list.append(0)
        return int_list
        
    f = open(file, 'r')
    int_list=[]
    for i in range(n_len):
        s = f.read(6)
        try:
            n = int(s)
        except:
            n = 0
        int_list.append(n)
    f.close()
    return int_list   

def add_to_file(date, data, file='datafile.dat', tmpfile='Null'):
    wrtfile=file
    if ( tmpfile != 'Null'): wrtfile=tmpfile
    if ( ( len(str(date)) < 8 ) ):
        raise Exception('date must be: CCYYMMDD')
    dates, datas = read_file(file=file)
    nv = len(data)
    if ( len(dates) != 0 ):
        ndv, ndt = datas.shape
        #print('DATA SHAPES', nv, ndv, ndt)
        if ( nv > ndv ):  # Add columns of zero to fill out data.
            ne = nv-ndv
            zeros = np.zeros((ne,ndt))
            datas=np.append(datas, zeros, axis=0)
        if ( nv < ndv ):
            print("Not enough data to fill", nv, ndv)
            return
    if ( len(dates) == 0 ):
        dates = np.atleast_1d(date)
        datas = np.reshape(np.array(data), (nv,1))
    else:
        iin = np.where(dates == date )
        if ( iin[0].size != 0 ):
            for iv in range(nv):
                datas[iv, iin] = data[iv]
        else:
            dates=np.append(dates, date)
            datas=np.append(datas, np.reshape(data, (nv,1)), axis=1)
    #print('SHAPE', datas.shape)
    if ( dates.size > 1 ):
        isort = np.argsort(dates)
        dates = dates[isort]
        for iv in range(nv):
            datas[iv, :] = datas[iv, isort]
        
    write_file(dates, datas, file=wrtfile)
    if ( wrtfile != file ):
        subprocess.call(['cp', wrtfile, file])
    return
    
def add_to_boolean_file(date, TorF, file='datafile.dat', tmpfile='Null'):
    if ( isinstance(TorF, list) ):
        for ie, element in enumerate(TorF):
            TorF[ie] = bool(element)
    elif ( isinstance(TorF, np.ndarray) ):
        TorF = TorF.astype(bool)
    else: 
        TorF = [bool(TorF)]
        
    wrtfile=file
    if ( tmpfile != 'Null'): wrtfile=tmpfile
    if ( ( len(str(date)) < 8 ) ):
        raise Exception('date must be: CCYYMMDD')
    dates, boole = read_date_and_boolean(file=file)
    nv = len(TorF)
    if ( len(dates) != 0 ):
        ndv, ndt = boole.shape
        if ( nv > ndv ):  # Add columns of zero to fill out data.
            ne = nv-ndv
            zeros = np.zeros((ne,ndt)).astype(bool)
            boole=np.append(datas, zeros, axis=0)
    if ( len(dates) == 0 ):
        dates = np.atleast_1d(date)
        boole = np.reshape(np.array(TorF), (nv,1))
    else:
        iin = np.where(dates == date )
        if ( iin[0].size != 0 ):
            for iv in range(nv):
                boole[iv, iin] = TorF[iv]
        else:
            dates=np.append(dates, date)
            boole=np.append(boole, np.reshape(TorF, (nv,1)), axis=1)
    #print('SHAPE', datas.shape)
    if ( dates.size > 1 ):
        isort = np.argsort(dates)
        dates = dates[isort]
        for iv in range(nv):
            boole[iv, :] = boole[iv, isort]
        
    write_date_and_boolean(dates, boole, file=wrtfile)
    if ( wrtfile != file ):
        subprocess.call(['cp', wrtfile, file])
    return

def get_from_boolean_file(date, file):
    dates, boole = read_date_and_boolean(file=file)
    boole = np.array(boole).astype(bool)
    iin = np.where(dates == date )
    if ( iin[0].size != 0 ):
        TorF = np.squeeze(boole[:,iin][:])
        if ( TorF.shape == () ): TorF=bool(TorF)
    else:
        TorF = None
    return TorF
    
def missing_from_file(file='datafile.dat', date_range=[datetime.datetime(2019,5,1), datetime.datetime(2019,6,30)]):
    dates, datas = read_file(file)
    missing = []
    if ( not isinstance(dates, np.ndarray ) ):
        print('File is empty or does not exist')
        print(os.path.isfile(file))
        return missing
     
    this_date = date_range[0]
    # Do not go past current date either.
    maxi_date = min([ date_range[1], datetime.datetime.now()])
    while ( this_date < maxi_date ):
        date_str = this_date.strftime("%Y%m%d%H")
        date_int = int(date_str)
        if (date_int not in dates):
            #print(date_int)
            missing.append(date_int)
        this_date = this_date + datetime.timedelta(days=1) 

    return missing 
     
     
