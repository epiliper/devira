process MULTIQC {
    container 'quay.io/staphb/multiqc:1.22.3'
    tag "OUTDIR: ${output_dir}"
    label = 'process_single' 

    input:
    val run_name
    val output_dir 
    path(fastp_files)

    output:
    path "./${run_name}_multiqc.html", emit: report
    """
    echo "FastP files:"
    echo "${fastp_files}"
    
    multiqc -f \\
        -o . \\
        -n "${run_name}_multiqc.html" \\
        -m fastp \\
        ${fastp_files}
    """
}
