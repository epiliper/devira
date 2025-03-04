# nextflow run main.nf --input example/piv_samplesheet.csv --output example_output -profile docker "$1" 
nextflow run main.nf --input srr_samplesheet.csv --run_kraken2 --kraken2_db ../adar_test_Data/krak --output out --skip_dedup -profile docker "$1" -process.echo
