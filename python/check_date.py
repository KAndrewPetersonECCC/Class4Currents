import datetime

def check_date(date, outtype=str, dtlen=8):
    if ( (outtype==str) or (outtype==int) ):
      if ( isinstance(date, datetime.datetime) or isinstance(date, datetime.date) ):
        if ( dtlen == 8 ):  
          datestr=date.strftime("%Y%m%d")
        elif ( dtlen == 10 ): 
          datestr=date.strftime("%Y%m%d%H")
        elif ( dtlen == 12 ): 
          datestr=date.strftime("%Y%m%d%H%M")
        elif ( dtlen == 14 ): 
          datestr=date.strftime("%Y%m%d%H%M%S")
        else:
          datestr=date.strftime("%Y%m%d")
      if ( isinstance(date, int) ):  datestr=str(date)
      if ( isinstance(date, str) ): datestr=date
      if ( len(datestr) < dtlen ):
          datestr=datestr+str(0).zfill(dtlen-len(datestr))
      if ( len(datestr) > dtlen ):
          datestr=datestr[:dtlen]
    if ( outtype ==int ): datestr=int(datestr)
    if ( outtype== datetime.datetime ):
      if ( isinstance(date, int) ): date=str(date)
      if ( isinstance(date, str) and ( len(date) == 8 ) ):
        datestr=datetime.datetime.strptime(date, '%Y%m%d')  
      elif ( isinstance(date, str) and ( len(date) == 10 ) ):
        datestr=datetime.datetime.strptime(date, '%Y%m%d%H')  
      elif ( isinstance(date, str) and ( len(date) == 12 ) ):
        datestr=datetime.datetime.strptime(date, '%Y%m%d%H%M')  
      elif ( isinstance(date, str) and ( len(date) == 14 ) ):
        datestr=datetime.datetime.strptime(date, '%Y%m%d%H%M%S')  
      if ( isinstance(date, datetime.datetime) ): datestr=date
      if ( isinstance(date, datetime.date) ): datestr=datetime.datetime(*date.timetuple()[:4])
    return datestr
