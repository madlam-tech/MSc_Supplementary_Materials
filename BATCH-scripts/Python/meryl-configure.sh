#!/bin/sh


#  Path to Canu.

bin="/scale_wlg_persistent/filesets/opt_nesi/CS400_centos7_bdw/Canu/2.2-GCC-9.2.0/bin"

#  Report paths.

echo ""
echo "Found perl:"
echo "  " `which perl`
echo "  " `perl --version | grep version`
echo ""
echo "Found java:"
echo "  " `which /opt/nesi/CS400_centos7_bdw/Java/17/bin/java`
echo "  " `/opt/nesi/CS400_centos7_bdw/Java/17/bin/java -showversion 2>&1 | head -n 1`
echo ""
echo "Found canu:"
echo "  " $bin/canu
echo "  " `$bin/canu -version`
echo ""


#  Environment for any object storage.

export CANU_OBJECT_STORE_CLIENT=
export CANU_OBJECT_STORE_CLIENT_UA=
export CANU_OBJECT_STORE_CLIENT_DA=
export CANU_OBJECT_STORE_NAMESPACE=
export CANU_OBJECT_STORE_PROJECT=




/scale_wlg_persistent/filesets/opt_nesi/CS400_centos7_bdw/Canu/2.2-GCC-9.2.0/bin/meryl -C k=22 threads=4 memory=2 \
  count segment=1/01 ../../BC10.seqStore \
> BC10.ms22.config.01.out 2>&1
exit 0
