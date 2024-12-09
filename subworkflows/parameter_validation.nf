include { BARCODE_HAMMING_CHECK } from "../bin/process.nf"

workflow validate_parameters {
    
    main:
        parameters_passed = true
        print("Checking parameters")

        try {
            file("${params.outdir}/report.html", checkIfExists:true)
            if (!params.overwrite_outdir) {
                parameters_passed = false
                println("Out Directory Already Exists, Please Provide New Out Directory Name or Allow Overwriting of Pre-existing directory")
            }
        }
        catch (Exception e) {
            
        }

        try {
            file(params.input_file, checkIfExists:true)
        }
        catch (Exception e) {
            parameters_passed = false
            println("Error - Input File or Directory Does not Exist")
        }

        Pinguscript.ping_start(nextflow, workflow, params)

        if (params.adaptor_sequence == "" && params.sample_file == "") {
            parameters_passed = false
            println ("Adaptor Sequence and Sample File cannot both be empty")
        }

        if (!params.plot_telo_length && !params.plot_vrr_length) {
            parameters_passed = false
            println ("Either VRR or Telomere Length Must be Set")
        }
        try {
            if (params.sample_file != ""){
                file(params.sample_file, checkIfExists:true)
            }
        }
        catch (Exception e) {
            parameters_passed = false
            println("Error - Sample File not Found")
        }

        //check to make barcodes are more different than barcode errors
        // try catch is currently non-functional will simply cause pipeline error
        if (params.sample_file != ''){
            try {
                barcode_check = BARCODE_HAMMING_CHECK(file(params.sample_file))
            }
            catch (Exception e){
                parameters_passed = false
                println "Supplied Barcode Sequences are too Similar for Demultiplexing with ${params.barcode_errors} Errors Allowed. Please reduce error amount."
            }
        // }
        }

    emit:
        passed = parameters_passed
}