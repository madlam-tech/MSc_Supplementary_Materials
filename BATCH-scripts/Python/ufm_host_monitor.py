#!/usr/bin/python
"""
This is a sample of ufm api using ufmsdk

We also recommend to implement the parse_input and the Usage functions, though 
is not mandatory (they are marked with comment)

@copyright:
    Copyright (C) Mellanox Technologies Ltd. 2001-2011.  ALL RIGHTS RESERVED.
    This software product is a proprietary product of Mellanox Technologies Ltd.
    (the "Company") and all right, title, and interest in and to the software product,
    including all associated intellectual property rights, are and shall
    remain exclusively with the Company.

    This software product is governed by the End User License Agreement
    provided with the software product.

@author:
    Izzat Kukhon
@date:
    25 6 2012
@changed: (major changes only)

"""

__docformat__ = "javadoc"



#python imports
import sys
import string
import os
import getopt

#UFM imports
sys.path.append(os.path.dirname(sys.path[0]))

from ufmsdk import infratools
from ufmsdk import ufmapi

server = None
user = 'admin'
password = '123456'
version = 1.1
ufm_api = None
ret_code = 0
logger = infratools.initLogger(os.path.basename(sys.argv[0]), 1, True, os.path.basename(sys.argv[0]) + ".log", False)

ib_attr = [ 'Infiniband_Normalized_CBW' ]

def Usage ():
    """
    Prints Usage Guide of the script.

    """
    print
    print "Usage: %s -s UFM_REMOTE_SERVER [-u USER] [-p PASSWORD] [-v] [-h]" % os.path.basename(sys.argv[0])
    print
    print "Options:"
    print "     -s UFM REMOTE SERVER       - Connect to remote UFM server"
    print "     [-u USER]                  - User to connect to UFM server"
    print "     [-p PASSWORD]              - Password to connect to UFM server"
    print "     [-v]                       - Show version "
    print "     [-h]                       - Show this help "
    print

    sys.exit(1)

if __name__ == '__main__':
    try:
       opts, args = getopt.getopt(sys.argv[1:], "l:s:u:p:h:v", ["help"])
    except getopt.error :
       raise Usage()
    try:
        for opt, arg in opts :
            if opt == '-h':
                Usage()
            elif opt == '-u':
                user = arg
            elif opt == '-p':
                password = arg
            elif opt == '-v':
                print os.path.basename(sys.argv[0]), "Version" , version
                sys.exit(2)
            elif opt == '-s':
                server = arg
        if (server is None or server ==" "):
                Usage()
        
        severMap = dict ()
        severMap ["Info"] = "1"
        severMap ["Minor"] = "2"
        severMap ["Warning"] = "4"
        severMap ["Critical"] = "6"
        
        #connect to the UFM Server
        ufm_api = ufmapi.UFMAPI(user=user, password=password, server=server)
        
        site = ufm_api.getSite('default')
        if site.isIB:
            attr = ib_attr
        else:
            attr = eth_attr     
           
        devices = ufm_api.getDevices('default')  

        devicesDict = dict()
        for d in devices:
            devicesDict [str(d.name)] = d

        alldata = ufm_api.monitorDataOnce(high_cls='Site', low_cls='Device', objects=['Grid.default'], attributes=attr , functions=['RAW'])     
        for data in alldata:
            name = data.dname.split(" ")[0]
            severity = devicesDict [str(data.object_name.split(".")[2])].severity.__str__()
            sys.stdout.write("\nBEGIN_NODE " + name +"\n")
            for couple in data.couple:
                sys.stdout.write("\tIB_Congestion" + " " + str(couple.value) +"\n")
                sys.stdout.write("\tUFM_Severity" + " " + severMap[str(severity)].__str__() +"\n")
                sys.stdout.write("=========================================" +"\n")
            sys.stdout.write(" " +"\n")  
        
    except Exception, exc:
        logger.error("main_function throw unhandled exception: %s" % exc)
        ret_code = 1
            
    if ret_code == 0:
        logger.info("CMU-UFM Integration script:  " + os.path.basename(sys.argv[0]) + " completed successfully")
    elif ret_code == 1:
        logger.error("CMU-UFM Integration script: " + os.path.basename(sys.argv[0]) + " completed with errors")
    else:
        logger.error("CMU-UFM Integration script: " + os.path.basename(sys.argv[0]) + " completed with unknown return code")
    sys.exit(ret_code)
