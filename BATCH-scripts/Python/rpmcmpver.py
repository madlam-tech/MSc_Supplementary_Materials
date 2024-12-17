"""
@copyright:
        Copyright (C) Mellanox Technologies Ltd. 2001-2012.  ALL RIGHTS RESERVED.

        This software product is a proprietary product of Mellanox Technologies Ltd.
        (the "Company") and all right, title, and interest in and to the software product,
        including all associated intellectual property rights, are and shall
        remain exclusively with the Company.

        This software product is governed by the End User License Agreement
        provided with the software product.
"""
__docformat__ = "javadoc"

import re  
import sys
import getopt

isalnum = re.compile('[^a-zA-Z0-9]')
usage_msg = "Usage: rpmcmpver ver1 ver2"   
      
def rpmvercmp(a, b):  
    # If they're the same, we're done  
    if a==b: return 0  
  
    def _gen_segments(val):  
        """ 
        Generator that splits a string into segments. 
        e.g., '2xFg33.+f.5' => ('2', 'xFg', '33', 'f', '5') 
        """  
        val = isalnum.split(val)  
        for dot in val:  
            res = ''  
            for s in dot:  
                if not res:  
                    res += s  
                elif (res.isdigit() and s.isdigit()) or \
                   (res.isalpha() and s.isalpha()):  
                    res += s  
                else:  
                    if res:  
                        yield res  
                    res = s  
            if res:  
                yield res  
      
    ver1, ver2 = a, b  
      
    # Get rid of the release number  
    ver1_rel, ver2_rel = None, None  
    if '-' in ver1: ver1, ver1_rel = ver1.rsplit('-')  
    if '-' in ver2: ver2, ver2_rel = ver2.rsplit('-')  
      
    l1, l2 = map(_gen_segments, (ver1, ver2))  
    while l1 and l2:  
        # Get the next segment; if none exists, done  
        try: s1 = l1.next()  
        except StopIteration: s1 = None  
        try: s2 = l2.next()  
        except StopIteration: s2 = None  
      
        if s1 is None and s2 is None: break  
        if (s1 and not s2): return 1  
        if (s2 and not s1): return -1  
  
        # Check for type mismatch  
        if s1.isdigit() and not s2.isdigit(): return 1  
        if s2.isdigit() and not s1.isdigit(): return -1  
  
        # Cast as ints if possible  
        if s1.isdigit(): s1 = int(s1)  
        if s2.isdigit(): s2 = int(s2)  
   
        rc = cmp(s1, s2)  
        if rc: return rc  
      
    # If we've gotten this far, check release numbers  
    if ver1_rel is not None and ver2_rel is not None:  
        return rpmvercmp(ver1_rel, ver2_rel)  
      
    return 0 


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
            if len(args) != 2:
                raise Usage(usage_msg) 
        except getopt.error, msg:
             raise Usage(usage_msg)
        rc = rpmvercmp(args[0], args[1])
        if rc == -1:
            rc = 2
        return rc
    except Usage, err:
        print err
        return 10

# exit 1 : ver1 > ver2
# exit 0 : ver1 == ver2
# exit 2 : ver1 < ver2 
if __name__ == "__main__":
    sys.exit(main())
    
    
    