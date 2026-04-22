include { KRAKEN2 } from '../modules/kraken2'

// classify reads to a kraken2 database containing host reads. Return only the reads that haven't been classified.
workflow REMOVE_HOST {
        take:
        reads_meta      // tuple val(meta), path(reads)
        krak_db_host    // path to kraken2 database

        main:
        KRAKEN2(reads_meta, krak_db_host, false, false)

        emit:
        reads =  KRAKEN2.out.unclassified_reads_fastq
        report = KRAKEN2.out.report
    }

