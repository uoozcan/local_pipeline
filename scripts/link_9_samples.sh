
#!/bin/bash



# Source directory

RAW_FASTQ_DIR="../raw_fastq"

INPUT_DIR="../pipeline_input"



# Function to create symlinks

create_link() {

    local SOURCE=$1

    local TARGET_NAME=$2

    

    if [ -f "${SOURCE}" ]; then

        ln -sf ${SOURCE} ${INPUT_DIR}/${TARGET_NAME}

        echo "✓ Linked: ${TARGET_NAME}"

    else

        echo "✗ WARNING: Source file not found: ${SOURCE}"

    fi

}



echo "Creating links for 9 samples..."



create_link "${RAW_FASTQ_DIR}/FH_2908_4_R1_GATCAG_L003_R1_001.trimmed.fastq.gz" "FH_2908_4_R1.fastq.gz"

create_link "${RAW_FASTQ_DIR}/FH_2908_4_R1_GATCAG_L003_R2_001.trimmed.fastq.gz" "FH_2908_4_R2.fastq.gz"



create_link "${RAW_FASTQ_DIR}/FH_2986_2_R1_ACTTGA_L003_R1_001.trimmed.fastq.gz" "FH_2986_2_R1.fastq.gz"

create_link "${RAW_FASTQ_DIR}/FH_2986_2_R1_ACTTGA_L003_R2_001.trimmed.fastq.gz" "FH_2986_2_R2.fastq.gz"



create_link "${RAW_FASTQ_DIR}/FH_3884_2_R1_S160_R1_001.trimmed.fastq.gz" "FH_3884_2_S160_R1.fastq.gz"

create_link "${RAW_FASTQ_DIR}/FH_3884_2_R1_S160_R2_001.trimmed.fastq.gz" "FH_3884_2_S160_R2.fastq.gz"



create_link "${RAW_FASTQ_DIR}/FH_3884_2_R1_S1_R1_001.trimmed.fastq.gz" "FH_3884_2_S1_R1.fastq.gz"

create_link "${RAW_FASTQ_DIR}/FH_3884_2_R1_S1_R2_001.trimmed.fastq.gz" "FH_3884_2_S1_R2.fastq.gz"



create_link "${RAW_FASTQ_DIR}/FH_4978_2_R1_CTTGTA_L004_R1_001.trimmed.fastq.gz" "FH_4978_2_R1.fastq.gz"

create_link "${RAW_FASTQ_DIR}/FH_4978_2_R1_CTTGTA_L004_R2_001.trimmed.fastq.gz" "FH_4978_2_R2.fastq.gz"



create_link "${RAW_FASTQ_DIR}/FH_5143_3_R1_GTGAAA_L008_R1_001.trimmed.fastq.gz" "FH_5143_3_R1.fastq.gz"

create_link "${RAW_FASTQ_DIR}/FH_5143_3_R1_GTGAAA_L008_R2_001.trimmed.fastq.gz" "FH_5143_3_R2.fastq.gz"



create_link "${RAW_FASTQ_DIR}/FH_5210_2_R1_TGACCA_L006_R1_001.trimmed.fastq.gz" "FH_5210_2_R1.fastq.gz"

create_link "${RAW_FASTQ_DIR}/FH_5210_2_R1_TGACCA_L006_R2_001.trimmed.fastq.gz" "FH_5210_2_R2.fastq.gz"



create_link "${RAW_FASTQ_DIR}/FH_5217_2_R1_ACAGTG_L006_R1_001.trimmed.fastq.gz" "FH_5217_2_R1.fastq.gz"

create_link "${RAW_FASTQ_DIR}/FH_5217_2_R1_ACAGTG_L006_R2_001.trimmed.fastq.gz" "FH_5217_2_R2.fastq.gz"



create_link "${RAW_FASTQ_DIR}/FH_5218_2_R1_GCCAAT_L006_R1_001.trimmed.fastq.gz" "FH_5218_2_R1.fastq.gz"

create_link "${RAW_FASTQ_DIR}/FH_5218_2_R1_GCCAAT_L006_R2_001.trimmed.fastq.gz" "FH_5218_2_R2.fastq.gz"



echo ""

echo "Done! Sample count:"

ls ${INPUT_DIR}/*_R1.fastq.gz | wc -l

echo "samples linked (should be 9)"



echo ""

echo "Samples:"

ls ${INPUT_DIR}/*_R1.fastq.gz | xargs -n1 basename | sed 's/_R1.fastq.gz//'

