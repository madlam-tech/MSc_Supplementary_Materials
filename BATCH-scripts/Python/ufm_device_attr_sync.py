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
    20 6 2012
@changed: (major changes only)

"""

__docformat__ = "javadoc"


#python imports
import sys
import os
import getopt

#UFM imports
sys.path.append(os.path.dirname(sys.path[0]))
from ufmsdk import infratools
from ufmsdk import ufmapi

sys.path.append(os.path.dirname(sys.path[0]))
logger = infratools.initLogger(os.path.basename(sys.argv[0]), 1, True, os.path.basename(sys.argv[0]) + ".log", False)

server = None  # Null for localhost
user = 'admin'
password = '123456'
version = 1.1
ufm_api = None
attrDict = None
attr_file = None
ret_code = 0
def Usage ():
    """
    Prints Usage Guide of the script.

    """
    print
    print "Usage: %s  -s UFM REMOTE SERVER  [-u USER] [-p PASSWORD] [-v] [-h]" % os.path.basename(sys.argv[0])
    print
    print "Options:"
    print "     -s UFM REMOTE SERVER       - Connect to remote UFM server"
    print "     -i FILE                    - Device Properties File"
    print "     [-u USER]                  - User to connect to UFM server"
    print "     [-p PASSWORD]              - Password to connect to UFM server"
    print "     [-v]                       - Show version "
    print "     [-h]                       - Show this help "
    print
    sys.exit(1)

def getDeviceIdentificationProperties(site_name, device_name):
    """
    Returns the identification properties of the given device.
    @param site_name:
        site name.
    @type site_name: 
        C{str}
    @param device_name:
        device name.
    @type device_name: 
        C{str}
    @return: 
        True - operation success.
    @rtype:
        C{bool}
    """
    global ret_code
    try:
        identificationProperties = ufm_api.getDeviceIdentificationProperties(site_name, device_name)
        if(identificationProperties is not None):
            logger.info("\nSuccess. Identification Properties for device with name '" + device_name + "' returned successfully.")
        else:
            logger.error("\nFailed. Identification Properties for device with name '" + device_name + "' not returned.")
            ret_code = 1
    except Exception, exc:
        logger.error("main_function throw unhandled exception: %s" % exc)
        ret_code = 1
    return identificationProperties

class deviceAttrs:
    def __init__(self, name):
        self.name = str(name).strip()
        self.is_manual_name = False
        self.manual_name = None
        self.is_manual_ip = False
        self.manual_ip = None
        
    def update (self, IdentificationProp):    
        self.is_manual_ip = IdentificationProp.is_manual_ip
        self.is_manual_name = IdentificationProp.is_manual_name
        self.manual_ip = IdentificationProp.manual_ip
        self.manual_name = IdentificationProp.manual_name
        
def deviceAttrEq (da1, da2):
    return (da1.is_manual_name == da2.is_manual_name) \
        and (da1.manual_name == da2.manual_name) \
        and (da1.is_manual_ip == da2.is_manual_ip) \
        and (da1.manual_ip == da2.manual_ip)
        
def loadIdProperties (attr_file):
    try:
        f = open(attr_file, "r")
        lines = f.readlines()
    except:
        lines = None
    finally:
        f.close ()       

    if ((lines == None) or (len(lines) == 0)):
        return 1
    
    for line in lines:
        line = line.strip()
        cells = line.split (" ")
        hostname = cells[0]
        ip = cells[1]
        guid = cells[2].replace("0x","",1)
        da = deviceAttrs (str(guid).strip())
        da.is_manual_ip = True
        da.is_manual_name = True
        da.manual_ip = ip
        da.manual_name = hostname
        attrDict [str(guid).strip()] = da
      
    return 0
    
def queryIdProperties (name):
    try:
        da = attrDict [name]
    except:
        da = deviceAttrs (name)    
    return da

def updateIdProperties (dev, IdentificationProp):
    
    name = str(dev.name).strip()
    da1 = deviceAttrs (name)
    da1.update (IdentificationProp)  
    update = False

    da2 = queryIdProperties (name)
    
    if (deviceAttrEq (da1, da2) == True):
        return
    
    if (da2.is_manual_name == True and da2.manual_name != None):
        IdentificationProp.is_manual_name = da2.is_manual_name
        IdentificationProp.manual_name = da2.manual_name
        update = True
        
    if (da2.is_manual_ip == True and da2.manual_ip != None):
        IdentificationProp.is_manual_ip = da2.is_manual_ip
        IdentificationProp.manual_ip = da2.manual_ip
        update = True
                
    if (IdentificationProp.manual_ip == None):
        IdentificationProp.manual_ip = "0.0.0.0"
        update = True
        
    if (update == True):
        print "Updating attributes for node %s"%name             
        ufm_api.updateDeviceManualIdentificationProperties ('default', dev.name, IdentificationProp)



if __name__ == '__main__':
    try:
        
        opts, args = getopt.getopt(sys.argv[1:], "s:u:i:p:h:v", ["help"])
    
    except getopt.error :
        raise Usage()
    try:
        for opt, arg in opts :
            if opt in ('-h', '--help'):
                Usage()
            elif opt == '-u':
                user = arg
            elif opt == '-i':
                attr_file = arg
            elif opt == '-p':
                password = arg
            elif opt == '-v':
                print os.path.basename(sys.argv[0]), "Version" , version
                sys.exit(2)
            elif opt == '-s':
                server = arg
        if (server is None or server == " " or attr_file is None):
                Usage()
        
        attrDict = dict ()  
            
        if (loadIdProperties(attr_file) != 0):
            print "Cannot load the device attributes file %s"%attr_file
            sys.exit(1)
        
      # connect to UFM server
        ufm_api = ufmapi.UFMAPI(user=user, password=password, server=server)
       
        #prints info for switches in the site
        for dev in ufm_api.getServers('default'):
            if (dev.systype != 'VM' and dev.systype != 'VMM'):
                devIdProp = getDeviceIdentificationProperties ('default', dev.name)
                updateIdProperties (dev, devIdProp)
            
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
    
