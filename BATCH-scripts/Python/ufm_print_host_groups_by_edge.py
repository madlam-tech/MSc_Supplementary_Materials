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
severMap = None
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
        
        
        #connect to the UFM Server
        ufm_api = ufmapi.UFMAPI(user=user, password=password, server=server)        
           
        servers = ufm_api.getServers('default')
        serversDict = dict()
        for s in servers:
            serversDict [str(s.name)] = s
        
        switches = ufm_api.getSwitches('default')
        switchesDict = dict()
        for s in switches:
            switchesDict [str(s.name)] = s
        
        switchHosts = dict () 
               
        links = ufm_api.getLinks('default')  
        for l in links:
            host = None
            try:
                host = serversDict [l.srcguid]
            except:
                try:
                    host = serversDict [l.destguid]
                except:
                    pass
                
            if (host == None):
                # inter switch link
                continue
            
            switch = None
            try:
                switch = switchesDict [l.srcguid]
            except:
                try:
                    switch = switchesDict [l.destguid]
                except:
                    pass
        
            switchname = "IB_switch_"+switch.name
            
            if (switchHosts.get(switchname) == None):
                switchHosts [switchname] = []
            switchHosts [switchname].append (host.dname)    
            
        CMU_PATH="/opt/cmu"
        CMU_ADD_USER_GROUP="cmu_add_user_group"
        CMU_ADD_TO_USER_GROUP="cmu_add_to_user_group"
        CMU_DEL_USER_GROUP="cmu_del_user_group"
        for switchname,hostsList in switchHosts.iteritems():
            sys.stdout.write(CMU_PATH+"/bin/"+CMU_DEL_USER_GROUP+" "+switchname+"\n")
            sys.stdout.write(CMU_PATH+"/bin/"+CMU_ADD_USER_GROUP+" "+switchname+"\n")
            sys.stdout.write(CMU_PATH + "/bin/" + CMU_ADD_TO_USER_GROUP + " -t " + switchname+" ")
            for hostname in hostsList:
		host=hostname.split()[0]
                sys.stdout.write(host+" ")
            sys.stdout.write("\n")
    except Exception, exc:
        print exc.__dict__
        logger.error("main_function throw unhandled exception: %s" % exc)
        ret_code = 1
            
    if ret_code == 0:
        logger.info("CMU-UFM Integration script:  " + os.path.basename(sys.argv[0]) + " completed successfully")
    elif ret_code == 1:
        logger.error("CMU-UFM Integration script: " + os.path.basename(sys.argv[0]) + " completed with errors")
    else:
        logger.error("CMU-UFM Integration script: " + os.path.basename(sys.argv[0]) + " completed with unknown return code")
    sys.exit(ret_code)
    
