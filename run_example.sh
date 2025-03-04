# nextflow run main.nf --input example/piv_samplesheet.csv --output example_output -profile docker "$1" 
nextflow run main.nf --input para_samplesheet.csv --run_kraken2 --kraken2_db ../adar_test_Data/krak --output para_out --skip_dedup -profile docker "$1"
