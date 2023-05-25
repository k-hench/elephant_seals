readonly CACTUS_IMAGE=$( echo "${snakemake_params}" | awk '{print $1}')
readonly SEQFILE=$( echo "${snakemake_params}" | awk '{print $2}')
readonly CACTUS_CORES=$( echo "${snakemake_params}" | awk '{print $3}')
readonly JOBSTORE_IMAGE="${snakemake_input}"
readonly SEQNAME=${SEQFILE##*/}
readonly RUN_ID=${SEQNAME%.txt}
readonly OUTPUTHAL=${RUN_ID}.hal
readonly CACTUS_SCRATCH=../results/cactus/scratch/${RUN_ID}

echo "store: ""${JOBSTORE_IMAGE}" &> "${snakemake_log[0]}"
echo "==================" &>> "${snakemake_log[0]}"
echo "file: " ${SEQFILE} &>> "${snakemake_log[0]}"
echo "==================" &>> "${snakemake_log[0]}"
echo "img: "${CACTUS_IMAGE} &>> "${snakemake_log[0]}"
echo "==================" &>> "${snakemake_log[0]}"
echo "cores: "${CACTUS_CORES} &>> "${snakemake_log[0]}"

CURRENT_DIR=$(pwd)
PARENT_DIR=${CURRENT_DIR%/*}

#   --overlay ${JOBSTORE_IMAGE} \
apptainer exec --cleanenv \
  --fakeroot --overlay ${CACTUS_SCRATCH} \
  --bind ${CACTUS_SCRATCH}/tmp:/tmp,${PARENT_DIR} \
  --env PYTHONNOUSERSITE=1 \
  ${CACTUS_IMAGE} \
  cactus-prepare \
  $(pwd)/${SEQFILE} \
  --cactusOptions '--maskMode none' \
  --defaultCores ${CACTUS_CORES} \
  --outDir /tmp/steps-output \
  --outSeqFile ${SEQFILE} \
  --outHal /tmp/steps-output/${OUTPUTHAL} \
  --jobStore /tmp/js > ${snakemake_output}