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
import argparse
import time
import ufm_slurm_utils
import logging
class UfmSlurmEpilog():

    def runSetup(self):
        self.log_name = ufm_slurm_utils.read_conf_file("log_file_name")
        self.script_name = os.path.basename(sys.argv[0])
	self.parse_args()
        format_str = "%(asctime)-15s  UFM-SLURM-Epilog JobID:" + self.args.job_id + "  %(levelname)-5s  : %(message)s"
        logging.basicConfig(format=format_str, level=logging.INFO, filename=self.log_name)
        self.user = ufm_slurm_utils.read_conf_file("ufm_server_user")
        self.password = ufm_slurm_utils.read_conf_file("ufm_server_pass")
        self.login_details = {"-u":self.user, "-p":self.password}
        self.ufm_server = None
        self.env_name = "slurm_env"
        self.sdk_ip = None


    def runTest(self):
        try:

            self.connect_to_ufm()
            if not ufm_slurm_utils.isRestSdkInstalled(self.sdk_ip):
                logging.error("UFM REST SDK on %s is not installed..exiting" % self.sdk_ip)
                sys.exit(1)
            ###### removing Logical server
            ls_name = 'slurm_job_{0}'.format(self.args.job_id)
            logging.info("Removing LS..")
            if ufm_slurm_utils.isLogicalServerExist(self.ufm_server, self.env_name, ls_name, self.sdk_ip, login_details=self.login_details):
                if ufm_slurm_utils.deleteLogicalServer(self.ufm_server, self.env_name, ls_name, self.sdk_ip, login_details=self.login_details):
                    logging.info("Logical Server: %s is removed successfully" % ls_name)
                else:
                    logging.error("LS: %s is not removed" % ls_name)
            else:
                logging.error("LS: %s is not exist to be removed" % ls_name)
        except Exception as exc:
            logging.error("Failed to remove LS:%s" % ls_name)
        finally:
            logging.info("## Done ##")

    def getJobInfo(self, job_id, attribute, attribute_level=1):
        command = "scontrol show job %s | grep ' %s' | awk '{print $%s}' | awk -F'=' '{print $2}'" % (job_id, attribute, attribute_level)
        _, res, _ = ufm_slurm_utils.runCommandStatus(command)
        return res

    def parse_args(self):
        self.parser = argparse.ArgumentParser()

        self.parser.add_argument("-i", "--job_id", action="store",
                                 dest="job_id", default=None,
                                 help="Job ID")

        self.args = self.parser.parse_args(sys.argv[1:])

    def connect_to_ufm(self):
        try:
            self.ufm_server, msg = ufm_slurm_utils.getUfmIP()
            print "server:", self.ufm_server, " msg:", msg
            if not self.ufm_server:
                if msg:
                    logging.error("Failed to connect to UFM: " + str(msg))
                else:
                    logging.error("Failed to connect to UFM.")
                    sys.exit(1)
            else:
                logging.info("connecting to server ... %s", self.ufm_server)
                if ufm_slurm_utils.IsUfmRunning(self.ufm_server):
                    logging.info("UFM: %s is running.." % self.ufm_server)
                else:
                    logging.error("UFM: %s is not running..exiting" % self.ufm_server)
                    sys.exit(1)
                self.sdk_ip = None
        except Exception as exc:
                logging.error("Can't connect %s : %s" % (self.ufm_server, exc))
                sys.exit(1)

if __name__ == '__main__':
    ret_code = 0
    try:
        test_name = os.path.basename(sys.argv[0])
        test = UfmSlurmEpilog()
        test.runSetup()
        test.runTest()
    except Exception as exc:
        print('Exception: ' + str(exc))
        print(test_name + " completed with errors")
        ret_code = 1
    sys.exit(ret_code)
