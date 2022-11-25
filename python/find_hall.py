import subprocess

def get_host():
    host=subprocess.check_output('hostname')
    if ( not isinstance(host, str) ): host=host.decode('utf-8')
    return host   
def find_hall():
    host=get_host()
    #print('HOST = '+ host)
    hall='null'
    if ( host[0:3] == 'cs1' ): hall='hall1'
    if ( host[0:3] == 'cs2' ): hall='hall2'
    if ( host[0:3] == 'cs3' ): hall='hall3'
    if ( host[0:3] == 'cs4' ): hall='hall4'
    if ( host[6:10] == 'ppp1' ): hall='hall1'
    if ( host[6:10] == 'ppp2' ): hall='hall2'
    if ( host[5:9] == 'ppp3' ): hall='hall3'
    if ( host[5:9] == 'ppp4' ): hall='hall4'
    if ( host[6:10] == 'ppp3' ): hall='hall3'
    if ( host[6:10] == 'ppp4' ): hall='hall4'
    if ( host[0:5] == 'hpcr3' ): hall='hall3'
    if ( host[0:5] == 'hpcr4' ): hall='hall4'
    if ( host[0:4] == 'ppp5' ): hall='hall5'
    if ( host[0:4] == 'ppp6' ): hall='hall6'
    return hall

def get_main_host(hall=None):
    if ( hall == None ):
        hall = find_hall()
    host_num = hall[4]
    host='eccc-ppp'+host_num
    return host

def get_ppp(hall=None):
    ppp=get_main_host(hall=hall)
    return ppp

def get_site(hall=None):
    if ( hall == None ):
        hall = find_hall()
    host_num = hall[4]
    site='site'+host_num
    return site

