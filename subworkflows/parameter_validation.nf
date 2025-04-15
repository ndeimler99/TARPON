/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BARCODE_HAMMING_CHECK } from "../bin/process.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow validate_parameters {
    
    main:
        Pinguscript.ping_start(nextflow, workflow, params)
        parameters_passed = true
        print("Checking parameters")

        // Checking to see if report.html already exists in outdir to prevent accidental overwrite of data
        try {
            file("${params.outdir}/report.html", checkIfExists:true)
            if (!params.overwrite_outdir) {
                parameters_passed = false
                println("Out Directory Already Exists, Please Provide New Out Directory Name or Allow Overwriting of Pre-existing directory")
            }
        }
        catch (Exception e) {
            
        }

        // checking to see if input files/directory exist
        try {
            file(params.input, checkIfExists:true)
        }
        catch (Exception e) {
            parameters_passed = false
            println("Error - Input File or Directory Does not Exist")
        }

        if (params.c_strand_only && params.strand_comparison) {
            parameters_passed = false
            println("C Strand Only can not be specific with Strand Comparison")
        }

        // either adaptor sequence (simplex) or sample_file (multiplex) must be present otherwise telomere ends are not able to be accurately identified
        if (params.capture_probe_sequence == "" && params.sample_file == "") {
            parameters_passed = false
            println ("Adaptor Sequence and Sample File cannot both be empty")
        }

        // If a sample file is specified, it must be a valid file
        try {
            if (params.sample_file != ""){
                file(params.sample_file, checkIfExists:true)
            }
        }
        catch (Exception e) {
            parameters_passed = false
            println("Error - Sample File not Found")
        }

        //check to ensure barcodes hamming distance is greater than the number of allowable errors in the barcode
        if (params.sample_file != ''){
            try {
                barcode_check = BARCODE_HAMMING_CHECK(file(params.sample_file))
            }
            catch (Exception e){
                parameters_passed = false
                println "Supplied Barcode Sequences are too Similar for Demultiplexing with ${params.barcode_errors} Errors Allowed. Please reduce error amount."
            }
        }

        if (params.fast_basecalled && params.pod5_dir == ""){
            parameters_passed = false
            println "Pod5 Directory must be set if initial reads were basecalled using the fast dorado models"
        }

    emit:
        passed = parameters_passed
}
