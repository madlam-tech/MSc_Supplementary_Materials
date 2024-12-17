def host_allowlist(handler, host):
    return True
c.ServerProxy.host_allowlist = host_allowlist
