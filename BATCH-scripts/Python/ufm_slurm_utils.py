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
import os
import re
import socket
import subprocess
import sys
import platform
import warnings
from IPy import IP
import datetime
import json
import time

try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import paramiko
except ImportError:
    # Failed to import paramiko.
    pass
def read_conf_file(key):
    file = open('ufm_slurm.conf', 'r')
    confs = file.read()
    # match = re.search(r'%s.*=\s*(.*)$'%key,confs)
    match = re.search(r'%s.*=(.*)' % key, confs)
    if match:
      return match.groups()[0]
    else:
      return None


ufm_machine_user = read_conf_file("ufm_machine_user")
ufm_machine_pass = read_conf_file("ufm_machine_pass")
SDK_BASE_DIR = "/opt/ufmrestsdk"




def isLocalHost(device, verbose=False):
    if device is None:
        return True
    if ("127" in str(device)):
        return True
    try:
        localhost_ip = socket.gethostbyaddr(socket.gethostname())[2]
        device_ip = socket.gethostbyaddr(device)[2]
        if device_ip == localhost_ip:
            return True
        else:
            return False
    except Exception, exc:
        return False


def runCommand(command, device=None, verbose=True, timeout=None):
    """
    Run Shell command on local or remote device
    Return command exit code.
    """
    if not isLocalHost(device):
        returncode = None
        conn = None
        connectionError = False

        try:
            conn = connect_to_remote(device, ufm_machine_user, '')
        except Exception:
            connectionError = True

        if not conn:
            connectionError = True

        if connectionError:
            try:
                conn = connect_to_remote(device, ufm_machine_user,
                                         ufm_machine_pass)
                connectionError = False
            except Exception:
                pass

        if not conn:
            connectionError = True

        if not connectionError:
            command = command + " 2>&1"
            if timeout:
                cmd = "timeout {0} {1}\n".format(timeout, command)
                _, out, _ = conn.exec_command(cmd, timeout=timeout)
            else:
                _, out, _ = conn.exec_command(command)
            returncode = out.channel.recv_exit_status()
            conn.close()
        else:
            returncode, _, _ = runSshCommand(command, device, verbose)

        return returncode

    else:
        with open(os.devnull, 'w') as tempf:
            command = command + " 2>&1"

            if timeout:
                command = "timeout {0} {1}\n".format(timeout, command)
            print("runCommand : " + str(command))
            proc = subprocess.Popen(command, shell=True,
                                    stdout=tempf, stderr=tempf)
            proc.communicate()
            return proc.returncode
def _run_cmd_v2(command, verbose=True):
    """
    Run Shell command.
    Direct output to subprocess.PIPE.
    Return command exit code, stdout and stderr.
    """
    proc = subprocess.Popen(command, shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    return (proc.returncode, stdout.strip(), stderr.strip())

def runCommandStatus(command, device=None, verbose=True):
    """
    Run Shell command on local or remote device
    Return command rd, status, stderr.
    """
    if not isLocalHost(device):
        returncode, stdout, stderr, conn = (None,) * 4
        connectionError = False

        try:
            conn = connect_to_remote(device, ufm_machine_user, '')
        except Exception:
            connectionError = True

        if not conn:
            connectionError = True

        if connectionError:
            try:
                conn = connect_to_remote(device, ufm_machine_user, ufm_machine_pass)
                connectionError = False
            except Exception:
                pass

        if not conn:
            connectionError = True

        if not connectionError:
            _, out, err = conn.exec_command(command)
            stdout = out.read()
            stderr = err.read()
            returncode = out.channel.recv_exit_status()
            conn.close()

        else:
            returncode, stdout, stderr = runSshCommand(
                command, device, verbose)
        return (returncode, stdout, stderr)

    else:
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        return (proc.returncode, str(stdout).strip(), str(stderr).strip())

def runSshCommand(command, device=None, verbose=True, logg=None):
    """
    Run Shell command on local or remote device
    Return command rd, status, stderr.
    """

    if not isLocalHost(device):
        returncode, stdout, stderr, conn = (None,) * 4
        connectionError = False

        user_name = ufm_machine_user
        SSH_CMD = "ssh {user_name}@{device} ""{command}"""
        command = SSH_CMD.format(**locals())

    proc = subprocess.Popen(command, shell=True, close_fds=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    return (proc.returncode, str(stdout).strip(), str(stderr).strip())



def connect_to_remote(remote_host, user, passwd):

    """
    connect to remote host over ssh.
    """
    try:
        conn = paramiko.SSHClient()
        conn.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        conn.connect(remote_host, username=user, password=passwd)
        return conn
    except Exception, exc:
        raise exc

def getUfmHostName(guid, ufm_server="127.0.0.1"):
    dict = getHostnameHostGuidDict(ufm_server)
    host_name = None
    for device in dict:
        if dict[device]['guid'] in guid or (dict[device]['host_guid'] and dict[device]['host_guid'] in guid):
            host_name = device
        else:
            for port in dict[device]['ports']:
                if port['port_guid'] in guid:
                    host_name = device
                break
    if host_name:
        print "UFM Hostname : ", host_name
        return host_name
    else:
        print "No UFM Hostname found."
        return None

def isFileExist(file_name, device=None):
    if not device:
        if os.path.exists(file_name):
            return True
        else:
            return False
    else:
        rc = runCommand("ls " + file_name, device)
        if rc == 0:
            return True
        else:
            return False

def haStatusLocal(device=None):
    """
    Show UFM (HA) local.
    """
    if haInstalled(device=device):
        status = ufmStatus(device=device)
        return status.split("========================================")[1]


def _parseHaStatus(status):
    """
    Parse ha status to list of tuples.
    Hand this function the output of haStatusRemote/haStatusLocal.
    """
    ret = {}
    for line in status.splitlines()[2:]:
        try:
            tokens = line.split()
            if len(tokens) == 2:
                ret[tokens[0]] = tokens[1]
            elif len(tokens) == 3:
                ret[' '.join(tokens[:2])] = tokens[2]
            elif len(tokens) == 4:
                ret[' '.join(tokens[:3])] = tokens[2]
        except Exception:
            pass
    return ret

def _ufmhaRunning(device=None):
    """
    Return True if ufmha is running according to its status, False otherwise.
    Prereq : device have to be in ha mode. Use haInstalled to guarantee that.
    """
    status = haStatusLocal(device)
    dict = _parseHaStatus(status)
    if not dict.get('UFM Server'):
        return False
    if dict.get('UFM Server') != 'Running':
        return False
    return True

def _ufmdRunning(device=None):
    """
    Return True if ufmd is running according to its status, False otherwise.
    """
    status = ufmStatus(device=device)
    if re.search(r'Running', status, re.IGNORECASE):
        return True
    return False

def ufmStatus(device=None):
    """
    Reutrn UFM status (string).
    """
    if haInstalled(device=device):
        command = "/etc/init.d/ufmha status"
    else:
        command = "/etc/init.d/ufmd status"

    if device:
        returncode, stdout, stderr = runCommandStatus(
            command, device)
    else:
        returncode, stdout, stderr = runCommandStatus(command)
    return stdout

def haInstalled(device=None):
    """
    Return True if UFM is installed in high-availability mode, False otherwise.
    """
    indicator = '/opt/ufm/indicators/ufm_ha'
    if device:
        return isFileExist(indicator, device)
    else:
        return os.path.exists(indicator)

def IsUfmRunning(device=None):
    """
    Check if UFM is running on a given device.
    """
    if haInstalled(device=device):
        return _ufmhaRunning(device=device)
    else:
        return _ufmdRunning(device=device)

def getHostnameHostGuidDict(ufm_server):
    hosts_dict = {}

    cmd = "ibhosts"
    _, result, _ = runCommandStatus(cmd, ufm_server)

    ibhosts = result.splitlines()
    for ibhost in ibhosts:
        if ibhost.startswith('Ca'):
            m = re.search('(.*): ' + '(\S*)' + ' ports(.*) "' + '(\S*)' +
                          '(.*)', ibhost)
            host_guid = m.group(2).replace('0x', '')
            if "4036E" in ibhost:
                continue
            host_name = m.group(4).split(".")[0].replace('"', '')
            if 'bridge' in host_name:
                continue
            if host_name not in hosts_dict:
                hosts_dict[host_name]=[]
            hosts_dict[host_name].append(host_guid)
    return hosts_dict

def _run_cmd_v2_remote(command, remote_host, verbose=False):
    """
    Run Shell command on remote host.
    Return command exit code, stdout, stderr.
    """
    returncode, stdout, stderr, conn = (None,) * 4
    connectionError = False

    if verbose:
        print("_run_cmd_v2_remote (%s) : %s" % (remote_host, command))
    try:
        conn = connect_to_remote(remote_host, ufm_machine_user, '')
    except Exception, exc:
        connectionError = True

    if not conn:
        connectionError = True

    if connectionError:
        try:
            conn = connect_to_remote(remote_host, ufm_machine_user, ufm_machine_pass)
            connectionError = False
        except Exception:
            pass

    if not conn:
        connectionError = True

    if not connectionError:
        _, out, err = conn.exec_command(command)
        stdout = out.read()
        stderr = err.read()
        returncode = out.channel.recv_exit_status()
        conn.close()

    return (returncode, stdout, stderr)

def getInternalIP(interface, device):
    command = """ip addr show %s | grep 'inet\\b' | awk '{print $2}' | cut -d/ -f1 | awk 'FNR ==1 {print}'""" % interface
    _, init_add_result, _ = _run_cmd_v2_remote(command, device)
    init_add_result = init_add_result.replace('\n', '')
    print("in getInternalIP: init_add_result:" + str(init_add_result))
    return init_add_result

def getJobNodesName():
    command = "scontrol show hostname $SLURM_JOB_NODELIST"
    _, result, _ = runCommandStatus(command)
    return result

def getDevicesInLS(ufm_server, env_name, ls_name, ufm_sdk_server, login_details):
    if not ufm_sdk_server:
        ufm_sdk_server = ufm_server
    query = "get_all"
    sdk_params = {"environment": env_name, "logical_server":ls_name}
    _, out, _ = runSdk(ufm_server, "computes", query, sdk_server=ufm_sdk_server, sdk_params=sdk_params, login_details=login_details)
    return json.loads(out, strict=False)


def isAllDevicesExistInLS(ufm_server, env_name, ls_name, device_names_added, ufm_sdk_server, login_details):
    device_list = getDevicesInLS(ufm_server, env_name, ls_name, ufm_sdk_server, login_details=login_details)
    device_list = [dev["name"] for dev in device_list]
    not_included = []
    if device_list:
    	for dev in device_names_added:
            if dev not in device_list:
                not_included.append(dev)
        if not_included:
            return False, not_included
    else:
	    return False, device_names_added
    return True, None

def getLogicalServer(ufm_server, env_name, ls_name, ufm_sdk_server, login_details):
    if not ufm_sdk_server:
        ufm_sdk_server = ufm_server
    query = "get" if ls_name else "get_all"
    sdk_params = {"environment": env_name}
    if ls_name:
        sdk_params["name"] = ls_name
    _, out, _ = runSdk(ufm_server, "servers", query, sdk_server=ufm_sdk_server, sdk_params=sdk_params, login_details=login_details)
    if "not found" not in out.lower():
        return json.loads(out, strict=False)

def isLogicalServerExist(ufm_server, env_name, ls_name, ufm_sdk_server, login_details):
    if not ufm_sdk_server:
        ufm_sdk_server = ufm_server
    out = getLogicalServer(ufm_server, env_name, ls_name, ufm_sdk_server=ufm_sdk_server, login_details=login_details)
    if out:
        return True
    return False

def createLogicalServer(ufm_server, env_name, ls_name, ufm_sdk_server, login_details):
    if not ufm_sdk_server:
        ufm_sdk_server = ufm_server
    _, out, _ = runSdk(ufm_server, "servers", "create", ufm_sdk_server, sdk_params={"environment":env_name, "name":ls_name}, login_details=login_details)
    return out


def deleteLogicalServer(ufm_server, env_name, ls_name, ufm_sdk_server, login_details):
    if not ufm_sdk_server:
        ufm_sdk_server = ufm_server
    if isLogicalServerExist(ufm_server, env_name, ls_name, ufm_sdk_server, login_details=login_details):
        _, out, _ = runSdk(ufm_server, "servers", "delete", sdk_server=ufm_sdk_server, sdk_params={"environment":env_name, "name":ls_name}, login_details=login_details)
        if "Error" in out:
            return False
        else:
            return True
    else:
        return False


def allocateDevicesToLS(ufm_server, env_name, ls_name, device_list, ufm_sdk_server, login_details):
    if not ufm_sdk_server:
        ufm_sdk_server = ufm_server
    _, ret_code, _ = runSdk(ufm_server, "servers", "allocate-computes-manually", ufm_sdk_server, sdk_params={"environment":env_name, "name":ls_name, "computes":device_list}, login_details=login_details)
    if ret_code:
        return True
    return False

def getEnvironment(ufm_server, env_name, ufm_sdk_server, login_details):
    if not ufm_sdk_server:
        ufm_sdk_server = ufm_server
    _, out, _ = runSdk(ufm_server, "environments", "get", sdk_server=ufm_sdk_server, sdk_params={"name":env_name}, login_details=login_details)
    if not "not found" in out.lower():
        return json.loads(out, strict=False)

def createEnvironment(ufm_server, env_name, ufm_sdk_server, login_details):
    if not ufm_sdk_server:
        ufm_sdk_server = ufm_server
    _, out, _ = runSdk(ufm_server, "environments", "create", ufm_sdk_server, sdk_params={"name":env_name}, login_details=login_details)
    return out

def isSystemExist(ufm_server, system_guid_name, ufm_sdk_server, login_details):
    if not ufm_sdk_server:
        ufm_sdk_server = ufm_server
    out = getSystem(ufm_server, system_guid_name, ufm_sdk_server=ufm_sdk_server, login_details=login_details)
    if out:
        return True
    return False

def getSystem(ufm_server, system_guid_name, ufm_sdk_server, login_details):
    if not ufm_sdk_server:
        ufm_sdk_server = ufm_server
    _, out, _ = runSdk(ufm_server, "systems", "get", sdk_server=ufm_sdk_server, sdk_params={"name":system_guid_name}, login_details=login_details)
    if not "not found" in out.lower():
        return json.loads(out, strict=False)

def getUfmIP():
    try:
        ufm_manual_ip = read_conf_file("ufm_server")
        if ufm_manual_ip:
            try:
                IP(ufm_manual_ip)
                return ufm_manual_ip, None
            except Exception as ex:
                return None, 'Error in parsing manual UMF IP. ' + str(ex)
        else:
            sm_guid = smPortGuid()
            ufm_server_hostname = None
            if sm_guid:
                ufm_server_hostname = getUfmHostName(sm_guid, "127.0.0.1")
            if ufm_server_hostname:
                ufm_server = getInternalIP("eth0", ufm_server_hostname)
                return ufm_server, None
            else:
                return None, "There is no discoverd/running UFM on the fabric." + "sm_guid: " + str(sm_guid) + " ufm_hostname: " + str(ufm_server_hostname)
    except Exception as ex:
            return None, ex

def runSdk(ufm_server,
           sdk_name,
           sdk_option,
           sdk_server=None,
           login_details=None,
           protocol="http",
           sdk_params={},
           sdk_payload={},
           extra_options={}):

    """
    @param ufm_server: UFM server IP address

    @param sdk_name: Name of UFM SDK, please refer to
                     available SDK names in constants
                     file: sdk_constants.Sdk

    @param sdk_server: Server to run SDK script

    @param login_details: Username and Password info
                          as Dictionary:
                          {
                          "-u": "user",
                          "-p": "pass"
                          }

    @param protocol: can be either http or https, by
                     default it takes http.

    @param sdk_option: Each SDK has a set of available
                       options pre-defined in constants
                       file. Please refer to sdk_constants
                       to see available options for each
                       SDK.

    @param sdk_params: parameters required by certain
                       options such as "get" option you
                       need to specify group_id or maybe
                       event_id etc.

                       Examples:
                       1) One parameter with one value:
                          {"job_id": "7"}
                       2) One Parameter with multiple values:
                          {"element_id": ["5", "11", ...]}

    @param sdk_payload: payload required by certain options
                        such as "add" option.

                        Examples:
                        {"elementName":"group1",
                         "description":"testing groups"}
    """
    if sdk_name == "servers":
        sdk_path = os.path.join(SDK_BASE_DIR, "servers.pyc")
    elif sdk_name == "environments":
        sdk_path = os.path.join(SDK_BASE_DIR, "environments.pyc")
    elif sdk_name == "systems":
        sdk_path = os.path.join(SDK_BASE_DIR, "systems.pyc")
    elif sdk_name == "computes":
        sdk_path = os.path.join(SDK_BASE_DIR, "computes.pyc")
    user = read_conf_file("ufm_server_user")
    password = read_conf_file("ufm_server_pass")
    login_details = {"-u":user, "-p":password}
    login_info = ""
    for key, value in login_details.iteritems():
        login_info += " %s %s" % (key, value)

    sdk_command_line = "python {sdk_path} -s {ufm_server} {login_info}"\
        " -r {protocol} {sdk_option}"\
        .format(**locals())

    # Extract SDK Parameters
    if sdk_params:

        params_list = []
        all_params = ""
        for key, value in sdk_params.iteritems():
            if type(value) is list:
                params_list.append("--%s=%s" % (key, ",".join(value)))
            else:
                params_list.append("--%s=%s" % (key, value))
        if len(params_list) > 1:
            all_params = " ".join(params_list)
        else:
            all_params = params_list[0]

        print("SDK Parameters Parsed: %s" % all_params)
        sdk_command_line += " " + all_params

    # Extract SDK Payload
    if sdk_payload:
        payload = " --payload=\""
        payload_list = []
        all_payload = ""
        for key, value in sdk_payload.iteritems():
            if type(value) is list:
                payload_list.append("%s=%s" % (key, ",".join(value)))
            else:
                payload_list.append("%s=%s" % (key, value))
        if len(payload_list) > 1:
            all_payload = "&".join(payload_list)
        else:
            all_payload = payload_list[0]
        payload += all_payload
        payload += "\""
        sdk_command_line += payload

    # Extract Extra Options
    if extra_options:
        extra_options_info = ""
        for key, value in extra_options.iteritems():
            parsed_val = "\"{0}\"".format(value)
            extra_options_info += " %s %s" % (key, parsed_val)
        sdk_command_line += extra_options_info

    return runCommandStatus(sdk_command_line,
                                        None,
                                        verbose=True)

def parseSdkOutput(sdk_output):
    """
    @summary: Parses and extracts data from UFM SDK
              Output string (a.k.a execution log)

    @return: If sdk_output was parsed successfully
             returns a tuple of (status_code, result)
             If parsing fails returns None.
    """
    print("[*] Parsing SDK output ...")
    status_code = None
    sdk_result = None

    search_result = re.search("\[\*\][\s\S]*results:((?:[\S\s]*\n)+)", sdk_output)

    if search_result is not None:
        result = search_result.group(1)
        result = result.replace("=" * 70, "")

        # Getting Status Code from SDK Output
        search_for_status = re.search("Error\s+(\d\d\d)\s+(.*)",
                                      result)

        if search_for_status is not None:
            status_code = search_for_status.group(1)

        # Getting Result from SDK Output
        search_for_result = re.search(">>[\s\S]*HTTP\s*response\s*text:((?:[\S\s]*\n)+)",
                                      result)
        if search_for_result is not None:
            sdk_result = search_for_result.group(1)
    return status_code, sdk_result

def sminfo(device=None):
    """
    Run sminfo command.
    """
    command = "/sbin/sminfo"
    if device:
        return runCommandStatus(command, device)
    else:
        return runCommandStatus(command)


def parse_sminfo(sminfo):
    """
    Parse sminfo output.
    If sminfo had en error, return None.
    If parsing failed, return None.
    Otherwise, return the sm portguid.
    """
    returncode, stdout, stderr = sminfo
    if returncode != 0:
        return None
    if "mad_rpc: _do_madrpc failed" in stdout:
        return None
    if "mad_rpc: _do_madrpc failed" in stderr:
        return None
    if "sminfo: iberror: failed: query" in stdout:
        return None
    if "sminfo: iberror: failed: query" in stderr:
        return None
    match = re.search("sminfo: sm lid [\d]+ sm guid ([\w\d]+),", stdout)
    if match:
        return match.group(1)

def _run_method_with_timeout(method, args=(), repeat=True,
                             checkSuccess=None, checkFailure=None,
                             timeout=120, retry_after=10):
    """
    Run method with timeout.
    method - function to be running.
    args - the arguments for the method.
    checkSuccess - callback to determine success.
    checkFailure - callback to determine failure.
    retry_after - amount of time to wait before retrying.
    """
    first = True
    current_time = datetime.datetime.now()
    start_time = datetime.datetime.now()
    while current_time - start_time < datetime.timedelta(seconds=timeout):
        if first or repeat:
            retCode = method(*args)
        if first:
            first = False
        if checkSuccess and checkSuccess(retCode):
            return 0
        if checkFailure and checkFailure(retCode):
            return 1
        time.sleep(retry_after)
        current_time = datetime.datetime.now()

def smPortGuid(timeOut=5, device=None):
    """Return sm portGuid."""
    print ("in smport guid, dev=", device)
    retCode = _run_method_with_timeout(sminfo,
                                       args=(device,),
                                       repeat=True,
                                       checkSuccess=lambda x: parse_sminfo(
                                           x) != None,
                                       timeout=timeOut,
                                       retry_after=10)
    if retCode == 0:
        r = sminfo(device)
        return parse_sminfo(r)
    return

def isRestSdkInstalled(ufm_server):
    ufmrestsdk_dir = "/opt/ufmrestsdk"
    uninstall_script = "%s/uninstall.sh" % ufmrestsdk_dir
    return isFileExist(uninstall_script, ufm_server)

def getSlurmLogicalServers(ufm_server, env_name):
    all_ls = getLogicalServer(ufm_server, env_name)
    slurm_ls_names = [ls['name'] for ls in all_ls if ls['name'].startswith('slurm_job_')]
    return slurm_ls_names

def getRunningSLurmJobsID():
    command = """squeue| awk '{if ($5=="R") print $1}'"""
    _, result, _ = runCommandStatus(command)
    jobs = result.splitlines()
    return jobs
