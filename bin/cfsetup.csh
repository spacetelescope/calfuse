# To execute, type "source cfsetup.csh"
#
set CF_DIR="${HOME}/calfuse"
set CF_VERSION="v3.2"
setenv PATH "${CF_DIR}/${CF_VERSION}/bin:${PATH}"
setenv CF_CALDIR "${CF_DIR}/${CF_VERSION}/calfiles"
setenv CF_PARMDIR "${CF_DIR}/${CF_VERSION}/parmfiles"
setenv CF_IDLDIR "${CF_DIR}/${CF_VERSION}/idl"
#
if ( `uname -s` !~ Darwin ) then
    if ( $?LD_LIBRARY_PATH ) then
	setenv LD_LIBRARY_PATH "${CF_DIR}/${CF_VERSION}/lib:${LD_LIBRARY_PATH}"
    else
	setenv LD_LIBRARY_PATH "${CF_DIR}/${CF_VERSION}/lib"
    endif
else
    if ( $?DYLD_LIBRARY_PATH ) then
	setenv DYLD_LIBRARY_PATH "${CF_DIR}/${CF_VERSION}/lib:${DYLD_LIBRARY_PATH}"
    else
	setenv DYLD_LIBRARY_PATH "${CF_DIR}/${CF_VERSION}/lib"
    endif
endif

set tcsh_path=`which tcsh`
if ( $tcsh_path != '/usr/local/bin/tcsh' ) then
    cp ${CF_DIR}/${CF_VERSION}/bin/calfuse.csh ${CF_DIR}/${CF_VERSION}/bin/calfuse_sv.csh
    sed -e "s#/usr/local/bin/tcsh#${tcsh_path}#g" ${CF_DIR}/${CF_VERSION}/bin/calfuse_sv.csh > ${CF_DIR}/${CF_VERSION}/bin/calfuse.csh
endif

rm -f ${CF_DIR}/${CF_VERSION}/bin/calfuse
ln -sf ${CF_DIR}/${CF_VERSION}/bin/calfuse.csh ${CF_DIR}/${CF_VERSION}/bin/calfuse
