#!/usr/bin/python
"""
Created on Sep 18, 2017

@author: lioros
@copyright:
        Copyright (C) Mellanox Technologies Ltd. 2001-2017.  ALL RIGHTS RESERVED.

        This software product is a proprietary product of Mellanox Technologies Ltd.
        (the "Company") and all right, title, and interest in and to the software product,
        including all associated intellectual property rights, are and shall
        remain exclusively with the Company.

        This software product is governed by the End User License Agreement
        provided with the software product.
"""
import sys
import os
import time
import ufm_slurm_utils
import argparse
import subprocess
import logging

user = ufm_slurm_utils.read_conf_file("ufm_server_user")
password = ufm_slurm_utils.read_conf_file("ufm_server_pass")
login_details = {"-u":user, "-p":password}
sdk_ip = None
env_name = "slurm_env"
device_list = []
script_name = os.path.basename(sys.argv[0])
log_name = ufm_slurm_utils.read_conf_file("log_file_name")
parser = None
args = None

def parse_args():
    global args

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--job_id", action="store",
                                dest="job_id", default=None,
                                help="Job ID")
    args = parser.parse_args(sys.argv[1:])

if __name__ == '__main__':

     # connect to UFM server
    try:
        parse_args()
        format_str = "%(asctime)-15s  UFM-SLURM-Prolog JobID:" + args.job_id + "  %(levelname)-5s  : %(message)s"
        logging.basicConfig(format=format_str, level=logging.INFO, filename=log_name)
        logging.info("Start JobID: " + args.job_id)
        server, msg = ufm_slurm_utils.getUfmIP()
        if not server:
            if msg:
                logging.error("Failed to connect to UFM: " + str(msg))
            else:
                logging.error("Failed to connect to UFM.")
            sys.exit(1)
        else:
            ls_name = 'slurm_job_{0}'.format(args.job_id)
            logging.info("connecting to server ... %s", server)
            if ufm_slurm_utils.IsUfmRunning(server):
                logging.info("UFM: %s is running.." % server)
            else:
                logging.error("UFM: %s is not running..exiting" % server)
                sys.exit(1)
    except Exception as exc:
            logging.error("Can't connect %s : %s", (server, exc))
            sys.exit(1)

# # Get GUID of Nodes
    sdk_ip = None
    if ufm_slurm_utils.isRestSdkInstalled(sdk_ip):
        logging.info("UFM REST SDK on %s is installed.." % 'SLURM Controller machine.')
    else:
        logging.error("UFM REST SDK on %s is not installed..exiting" % 'SLURM Controller machine.')
        sys.exit(1)
    nodes_names = ufm_slurm_utils.getJobNodesName()
    nodes = nodes_names.splitlines()
    if not nodes:
        logging.error("Couldn't get nodes of the job.. existing..")
        sys.exit(1)
    hosts_guids = ufm_slurm_utils.getHostnameHostGuidDict(server)
    if not hosts_guids:
        logging.error("Couldn't get GUID of nodes. existing..")
        sys.exit(1)
    logging.info("The Job Nodes are:")
    for node in nodes:
        found = False
        if node in hosts_guids.keys():
            for guid in hosts_guids[node]:
                if ufm_slurm_utils.isSystemExist(server, guid, sdk_ip, login_details=login_details):
                    device_list.append(guid)
                    found = True      
                    break
            if found:
                logging.info(node)
            else:
                logging.error(node + " guid is not found in UFM fabric. It could not be added to the logical server")
        else:
            logging.warning(node + " is not in the same UFM fabric. It could not be added to the logical server.")
    if not device_list:
	logging.error("No GUIDS of nodes are found to add")
	sys.exit(1)
# # Create environment
    # If environment does not exist, create it
    if ufm_slurm_utils.getEnvironment(server, env_name, sdk_ip, login_details=login_details) is not None:
            logging.info('Environment %s already exists' % env_name)
    else:
        logging.info('Creating environment %s ...' % env_name)
        ret_code = ufm_slurm_utils.createEnvironment(server, env_name, sdk_ip, login_details=login_details)
    if ufm_slurm_utils.getEnvironment(server, env_name, sdk_ip, login_details=login_details) is None:
        logging.error("Failed to create enironment: %s" % env_name)
        sys.exit(1)

# # Create Logical Server

    try:
        ufm_slurm_utils.deleteLogicalServer(server, env_name, ls_name, sdk_ip, login_details=login_details)
        logging.info('Creating logical server: %s ...' % ls_name)
        if ufm_slurm_utils.createLogicalServer(server, env_name, ls_name, sdk_ip, login_details=login_details):
            logging.info('LS %s is created.' % ls_name)
        else:
            logging.info('Failed to create LS: %s.' % ls_name)
            sys.exit(1)
    except Exception, exc:
        logging.error("Error in creating LS %s: %s" % (ls_name, str(exc)))
        sys.exit(1)


# # Assign nodes to Logical Server
    try:
        logging.info("Allocate nodes to LS")
        if not device_list:
            logging.error("Error in allocate nodes to LS %s: No nodes related to UFM server are found." % ls_name)
            sys.exit(1)
        ufm_slurm_utils.allocateDevicesToLS(server, env_name, ls_name, device_list, sdk_ip, login_details=login_details)
        is_all_exist, lst = ufm_slurm_utils.isAllDevicesExistInLS(server, env_name, ls_name, device_list, sdk_ip, login_details=login_details)
        if is_all_exist:
            logging.info("Success. The allocation of the nodes to LS '" + ls_name + "' passed")
        else:
            logging.error("Not all the nodes are allocated to the logical server.")
            if lst:
                logging.error("The following nodes are not added:")
                for dev in lst:
                    logging.error(dev)
    except Exception, exc:
        logging.error(
        "Error in allocating nodes to LS. %s" % str(exc))
        ret_code = 1



