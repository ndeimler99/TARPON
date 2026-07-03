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

        available = Runtime.runtime.availableProcessors()
        if (params.threads as int > available){
            println("Specified number of threads ${params.threads} is greater than the number of available Processors ${available}")
            println("Number of threads is reset to ${available}")
        }

        if (params.inheritance_mode){
            try {
                file(params.maternal_telobam, checkIfExists:true)
            }
            catch (Exception e) {
                parameters_passed = false
                println("Error - Maternal Telobam must be provided")
            }

            try {
                file(params.paternal_telobam, checkIfExists:true)
            }
            catch (Exception e) {
                parameters_passed = false
                println("Error - Paternal Telobam must be provided")
            }

            try {
                file(params.offspring_telobam, checkIfExists:true)
            }
            catch (Exception e) {
                parameters_passed = false
                println("Error - Offspring Telobam must be provided")
            }

            try {
                file(params.maternal_stats, checkIfExists:true)
            }
            catch (Exception e) {
                parameters_passed = false
                println("Error - Maternal Stats File must be provided")
            }

            try {
                file(params.paternal_stats, checkIfExists:true)
            }
            catch (Exception e) {
                parameters_passed = false
                println("Error - Paternal Stats File must be provided")
            }

            try {
                file(params.offspring_stats, checkIfExists:true)
            }
            catch (Exception e) {
                parameters_passed = false
                println("Error - Offspring Stats File must be provided")
            }
        }
        else if (params.cluster_concordance) {
            try {
                file(params.concordance_sample_file, checkIfExists: true)
            }
            catch (Exception e) {
                parameters_passed = false
                println("Error - Cluster Concordance File Must be Provided")
            }
        }
        else{
            if (params.recluster_only){
                try {
                    file(params.input, checkIfExists:true)
                }
                catch (Exception e) {
                    parameters_passed = false
                    println("Error - Telomeric Bam File Does not Exist")
                }

                try {
                    file(params.telomeric_stats, checkIfExists:true)
                }
                catch (Exception e) {
                    parameters_passed = false
                    println("Error - Telomeric Stats File Does not Exist")
                }
            }
            else if (params.alignment) {
                try {
                    file(params.reference, checkIfExists:true)
                }
                catch (Exception e) {
                    parameters_passed = false
                    println("Error - Reference File Does not Exist")
                }
            }
            else {
                
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

                if (params.repeat == params.mutant){
                    println("Mutant repeat can not be the same as wild type repeat")
                    parameters_passed = false
                }
                if (params.c_strand_only && params.strand_comparison) {
                    parameters_passed = false
                    println("C Strand Only can not be specific with Strand Comparison")
                }

                // either adaptor sequence (simplex) or sample_file (multiplex) must be present otherwise telomere ends are not able to be accurately identified
                if (params.capture_probe_sequence == "" && params.sample_file == "" && !params.no_capture_probe) {
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
            }
        }
    emit:
        passed = parameters_passed
}
