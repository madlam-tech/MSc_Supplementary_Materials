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
    24 6 2012
@changed: (major changes only)

"""

__docformat__ = "javadoc"

#python imports
import sys
import string
import os
import getopt
import time
from xml.etree import ElementTree

#UFM imports
from ufmsdk import infratools
from ufmsdk import ufmapi

sys.path.append(os.path.dirname(sys.path[0]))


from suds.client import Client
from suds.transport.https import HttpAuthenticated

server = None
report_def = None
user = 'admin'
password = '123456'
version = 1.1
ret_code = 0
logger = infratools.initLogger(os.path.basename(sys.argv[0]), 1, True, os.path.basename(sys.argv[0]) + ".log", False)
ufm_api = None

def Usage ():
    print
    print "Usage: %s -s UFM REMOTE SERVER -i REPORT_CONF_FILE [-u USER] [-p PASSWORD] [-v] [-h]" % os.path.basename(sys.argv[0])
    print
    print "Options:"
    print "     -s UFM REMOTE SERVER       - Connect to remote UFM server"
    print "     -i REPORT_CONF_FILE        - Report Configuration XML"
    print "     [-u USER]                  - User to connect to UFM server"
    print "     [-p PASSWORD]              - Password to connect to UFM server"
    print "     [-v]                       - Show version "
    print "     [-h]                       - Show this help "
    print

    sys.exit(1)


def xmltodict(element):
#    if not isinstance(element, ElementTree.Element):
#        raise ValueError("must pass xml.etree.ElementTree.Element object")

    def xmltodict_handler(parent_element):
        result = dict()
        for element in parent_element:
            if len(element):
                obj = xmltodict_handler(element)
            else:
                obj = element.text

            if result.get(element.tag):
                if hasattr(result[element.tag], "append"):
                    result[element.tag].append(obj)
                else:
                    result[element.tag] = [result[element.tag], obj]
            else:
                result[element.tag] = obj
        return result

    return {element.tag: xmltodict_handler(element)}

def xmlfiletodict(filename):
    return xmltodict(ElementTree.parse(filename).getroot())

if __name__ == '__main__':
    try:
       opts, args = getopt.getopt(sys.argv[1:], "l:s:u:i:p:hv", ["help"])
    
    except getopt.error :
       raise Usage()
    try:
        for opt, arg in opts :
            if opt == '-h':
                Usage()
            elif opt == '-u':
                user = arg
            elif opt == '-i':
    	           report_def = arg
            elif opt == '-p':
                password = arg
            elif opt == '-v':
                print os.path.basename(sys.argv[0]), "Version" , version
                sys.exit(2)
            elif opt == '-s':
                server = arg
        if (server is None or server ==" " or report_def is None or report_def == " "):
            Usage()
            
        #connect to the UFM Server
        ufm_api = ufmapi.UFMAPI(user=user, password=password, server=server)
              
        logger.info("start report")
    
        report_def_dict = xmlfiletodict(report_def)
        report_id = ufm_api.startReport(report_def_dict['report_res'])
        logger.info("Report ID: " + str(report_id))
        logger.info("Waiting for report data")
        
        i = 0
        delay=10
        timeout=300.0/delay
        report_is_completed = False
        time.sleep(1)
        while ((i < timeout) and not(report_is_completed)):
            report_data = ufm_api.getReportData(report_id)
            if report_data[4] == "Processing":
                logger.info("Report is being processed. Waiting for %d seconds" % delay)
                time.sleep(delay)
            elif (len(report_data) > 5 and report_data[5] == "Completed"):
                logger.info("Report has been completed!")
                html = report_data[4]
                report_is_completed = True
                break
        if (not (report_is_completed)):
            logger.info("Report has not completed in %d seconds. Exit on timeout!" % (timeout * delay))        
                
        logger.info("Done")

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

