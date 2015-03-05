# To execute, type "source cfsetup.sh"
#
CF_DIR="${HOME}/calfuse"
CF_VERSION="v3.2"
export PATH="${CF_DIR}/${CF_VERSION}/bin:${PATH}"
export LD_LIBRARY_PATH="${CF_DIR}/${CF_VERSION}/lib:${LD_LIBRARY_PATH}"
export DYLD_LIBRARY_PATH="${CF_DIR}/${CF_VERSION}/lib:${DYLD_LIBRARY_PATH}"
export CF_CALDIR="${CF_DIR}/${CF_VERSION}/calfiles"
export CF_PARMDIR="${CF_DIR}/${CF_VERSION}/parmfiles"
export CF_IDLDIR="${CF_DIR}/${CF_VERSION}/idl"
export CF_HISTDIR=`pwd`

sh_path=`which sh`
if [ $sh_path != '/bin/sh' ]; then
    cp ${CF_DIR}/${CF_VERSION}/bin/calfuse.sh ${CF_DIR}/${CF_VERSION}/bin/calfuse_sv.sh
    sed -e "s#/bin/sh#${sh_path}#g" ${CF_DIR}/${CF_VERSION}/bin/calfuse_sv.sh > ${CF_DIR}/${CF_VERSION}/bin/calfuse.sh
fi

rm -f ${CF_DIR}/${CF_VERSION}/bin/calfuse
ln -sf ${CF_DIR}/${CF_VERSION}/bin/calfuse.sh ${CF_DIR}/${CF_VERSION}/bin/calfuse
