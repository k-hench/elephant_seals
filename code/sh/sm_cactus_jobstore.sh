readonly CURRENT_DIR=$( pwd )
readonly CACTUS_IMAGE=$( echo "${snakemake_params}" | awk '{print $1}')
#readonly JOBSTORE_SZ=$( echo "${snakemake_params}" | awk '{print $2}')
readonly JOBSTORE_IMAGE=${snakemake_output}
readonly SEQFILE=$( echo "${snakemake_params}" | awk '{print $2}')
readonly SEQNAME=${SEQFILE##*/}
readonly RUN_ID=${SEQNAME%.txt}
readonly CACTUS_SCRATCH=../results/cactus/scratch/${RUN_ID}

echo "${CURRENT_DIR}" &> "${snakemake_log[0]}"
echo "==================" &>> "${snakemake_log[0]}"
echo "${SEQNAME}" &>> "${snakemake_log[0]}"
echo "==================" &>> "${snakemake_log[0]}"
echo "${CACTUS_IMAGE}" &>> "${snakemake_log[0]}"
echo "==================" &>> "${snakemake_log[0]}"
echo ${JOBSTORE_SZ} &>> "${snakemake_log[0]}"

restart=''
mkdir -p -m 777 ${CACTUS_SCRATCH}/upper ${CACTUS_SCRATCH}/work
touch ${JOBSTORE_IMAGE}
# truncate -s ${JOBSTORE_SZ} "${JOBSTORE_IMAGE}"
# apptainer exec --bind $(pwd) ${CACTUS_IMAGE} mkfs.ext3 -d ${CACTUS_SCRATCH} "${JOBSTORE_IMAGE}"

mkdir -m 700 -p ${CACTUS_SCRATCH}/tmp/steps-output
mkdir -p ../results/cactus/cactus_wd

cd ${CACTUS_SCRATCH}/tmp/steps-output
ln -s ${CURRENT_DIR}/../data/genomes/*.fa.gz ./