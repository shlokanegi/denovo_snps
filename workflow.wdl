version 1.0

workflow vcfeval_denovos_wf {
    meta {
	    author: "Shloka Negi"
        email: "shnegi@ucsc.edu"
        description: "Workflow built over vcfeval with custom-filtering to generate annotated denovo candidates in a trio"
    }

    parameter_meta {
        FAMILY: "Family name"
        CHILD_VCF: "Child's long-read harmonized VCF. Gzipped"
        FATHER_VCF: "Father's long-read harmonized VCF. Gzipped"
        MOTHER_VCF: "Mother's long-read harmonized VCF. Gzipped"
        FATHER_BAM: "Sorted and indexed father's BAM with the long reads. Optional. If present, denovos will be re-evaluated with read support. Helps filtering false-positives out"
        FATHER_BAM_INDEX: "Index for father's BAM"
        MOTHER_BAM: "Sorted and indexed mother's BAM with the long reads. Optional. If present, denovos will be re-evaluated with read support. Helps filtering false-positives out"
        MOTHER_BAM_INDEX: "Index for mother's BAM"
        REF: "Reference genome fasta"
        LOW_COMPLEXITY_BED: "GIAB low-complexity regions in GRCh38"
        SNPEFF_DB: "SNPeff annotation bundle (zip file)"
        SNPEFF_DB_NAME: "Name of SNPeff annotation bundle. E.g. GRCh38.105"
        CLINVAR_VCF: "VCF with all variants in ClinVar and their clinical significance (CLNSIG INFO field). Must be sorted, bgzipped, and indexed."
        CLINVAR_VCF_INDEX: "Index for CLINVAR_VCF (.tbi file)."
        DBNSFP_DB: "dbNSFP annotation bundle (block-gzip)"
        DBNSFP_DB_INDEX: "Index for DBNSFP_DB"
        GNOMAD_VCF: "VCF with all variants in gnomAD and their allele frequency (AF INFO field). Must be sorted, bgzipped, and indexed."
        GNOMAD_VCF_INDEX: "Index for GNOMAD_VCF (.tbi file)."
    }

    input {
        String FAMILY
        File CHILD_VCF
        File FATHER_VCF
        File MOTHER_VCF
        File? FATHER_BAM
        File? FATHER_BAM_INDEX
        File? MOTHER_BAM
        File? MOTHER_BAM_INDEX
        File REF
        File? LOW_COMPLEXITY_BED
        File? SNPEFF_DB
        String? SNPEFF_DB_NAME
        File? CLINVAR_VCF
        File? CLINVAR_VCF_INDEX
        File? DBNSFP_DB
        File? DBNSFP_DB_INDEX
        File? GNOMAD_VCF
        File? GNOMAD_VCF_INDEX
    }

    ## Runs vcfeval on child and parents' VCFs (trio) to generate total denovo variants per family. 
    ## Filters denovos by variant QUAL threshold.
    call run_vcfeval {
        input:
        child_vcf=CHILD_VCF,
        father_vcf=FATHER_VCF,
        mother_vcf=MOTHER_VCF,
        ref=REF,
        family=FAMILY
    }

    File denovos_vcf = run_vcfeval.vcf

    ## Restrict to denovos outside low-complexity regions.
    if(defined(LOW_COMPLEXITY_BED)) {
        call subset_denovos_by_region {
            input:
            input_vcf=denovos_vcf,
            lc_bed=select_first([LOW_COMPLEXITY_BED])
        }
    }

    File denovos_subset_vcf = select_first([subset_denovos_by_region.vcf, denovos_vcf])
    
    ## Effective denovo variant filtering using parent's VCFs to remove possible false-positives
    call run_filtering {
        input:
        input_vcf=denovos_subset_vcf,
        mother_vcf=MOTHER_VCF,
        father_vcf=FATHER_VCF
    }
    
    # File denovos_filtered_vcf = run_filtering.vcf
    File denovos_flagged_snps_vcf = run_filtering.flagged_snps

    ## Annotate and filter rare de novo snps
    if (defined(GNOMAD_VCF) && defined(GNOMAD_VCF_INDEX)) {
    call annotate_denovo_snps {
        input:
        input_vcf=denovos_flagged_snps_vcf,
        gnomad_vcf=GNOMAD_VCF,
        gnomad_vcf_index=GNOMAD_VCF_INDEX
        }
    }

    File denovos_rare_vcf = select_first([annotate_denovo_snps.vcf, denovos_flagged_snps_vcf])

    ## Annotate rare denovos with SnpEff, ClinVar and dbNSFP
    if(defined(SNPEFF_DB) && defined(SNPEFF_DB_NAME)){
        call annotate_with_snpeff {
            input:
            input_vcf=denovos_rare_vcf,
            snpeff_db=select_first([SNPEFF_DB]),
            db_name=select_first([SNPEFF_DB_NAME])
        }
    }

    File denovos_rare_snpeff_vcf = select_first([annotate_with_snpeff.vcf, denovos_rare_vcf])

    ## Annotate SNPs with presence in ClinVar
    # NOTE: first filter variants to keep those with high/moderate impact or with predicted loss of function (speeds up DB matching a lot)
    if (defined(CLINVAR_VCF) && defined(CLINVAR_VCF_INDEX) && defined(DBNSFP_DB) && defined(DBNSFP_DB_INDEX)){
        call subset_annotate_smallvars_with_db {
            input:
            input_vcf=denovos_rare_snpeff_vcf,
            clinvar_vcf=select_first([CLINVAR_VCF]),
            clinvar_vcf_index=select_first([CLINVAR_VCF_INDEX]),
            dbnsfp_db = select_first([DBNSFP_DB]),
            dbnsfp_db_index = select_first([DBNSFP_DB_INDEX])
        }
    }

    File denovos_rare_snpeff_clinvar_dbnsfp_vcf = select_first([subset_annotate_smallvars_with_db.vcf, denovos_rare_snpeff_vcf])

    ## Effective denovo variant filtering using parent's BAMs to flag false-positives
    if (defined(FATHER_BAM) && defined(FATHER_BAM_INDEX) && defined(MOTHER_BAM) && defined(MOTHER_BAM_INDEX)){
        call validate_denovos {
            input:
            input_vcf=denovos_rare_snpeff_clinvar_dbnsfp_vcf,
            mom_bam=select_first([MOTHER_BAM]),
            mom_bam_index=select_first([MOTHER_BAM_INDEX]),
            dad_bam=select_first([FATHER_BAM]),
            dad_bam_index=select_first([FATHER_BAM_INDEX]),
            reference_fasta=REF
        }
    }

    File final_vcf = select_first([validate_denovos.vcf, denovos_rare_snpeff_clinvar_dbnsfp_vcf])


    output {
        # File denovos_all = denovos_vcf                                      # all denovos genome-wide
        # File denovos_subset = denovos_subset_vcf                            # denovos in high-confidence regions
        # File denovos_filtered = denovos_filtered_vcf                        # denovos in high-confidence regions post-filtering
        # File denovos_flagged_snps = denovos_flagged_snps_vcf                # filtered denovo SNPs not called in both parents
        # File? denovos_snps_annotated = annotate_denovo_snps.annotated_vcf   # filtered denovo SNPs annotated with gnomad AF
        File denovos_snps_rare = denovos_rare_vcf                             # rare snp denovos
        File denovos_final_vcf = final_vcf
    }
}

task run_vcfeval {
    input {
        String family
        File child_vcf
        File father_vcf
        File mother_vcf
        File ref
        Int memSizeGB = 48
        Int threadCount = 32
        Int diskSizeGB = 5*round((size(child_vcf, "GB") + size(father_vcf, "GB")) + size(mother_vcf, "GB")) + 20
    }

    command <<<
        set -eux -o pipefail

        ## Preprocess VCFs
            # Normalize VCFs
            # Only keep PASS variants
            # Remove HOM REF variants and variants "missing" GT calls
            # Remove poor QUAL variants (<20)
            # Remove very long variants (>30 bps) to make sure we are not dealing with SVs
        zcat ~{child_vcf} | bcftools norm -m -any --threads ~{threadCount} | bcftools filter -e 'FILTER!="PASS"' | bcftools filter -e '(GT="0/0"| GT="0|0"| GT="./.")' | bcftools filter -i '(STRLEN(REF)<30 & STRLEN(ALT)<30)' -Oz -o child.filtered.vcf.gz
        zcat ~{father_vcf} | bcftools norm -m -any --threads ~{threadCount} | bcftools filter -e 'FILTER!="PASS"' | bcftools filter -e '(GT="0/0"| GT="0|0"| GT="./.")' | bcftools filter -i '(STRLEN(REF)<30 & STRLEN(ALT)<30)' -Oz -o father.filtered.vcf.gz
        zcat ~{mother_vcf} | bcftools norm -m -any --threads ~{threadCount} | bcftools filter -e 'FILTER!="PASS"' | bcftools filter -e '(GT="0/0"| GT="0|0"| GT="./.")' | bcftools filter -i '(STRLEN(REF)<30 & STRLEN(ALT)<30)' -Oz -o mother.filtered.vcf.gz

        ## Generate index
        tabix -p vcf child.filtered.vcf.gz
        tabix -p vcf father.filtered.vcf.gz
        tabix -p vcf mother.filtered.vcf.gz

        ## link the database VCF to make sure their indexes can be found
        ln -s child.filtered.vcf.gz child.vcf.gz
        ln -s father.filtered.vcf.gz father.vcf.gz
        ln -s mother.filtered.vcf.gz mother.vcf.gz
        ln -s child.filtered.vcf.gz.tbi child.vcf.gz.tbi
        ln -s father.filtered.vcf.gz.tbi father.vcf.gz.tbi
        ln -s mother.filtered.vcf.gz.tbi mother.vcf.gz.tbi

        ## Generate RTG Sequence data file for given reference genome
        rtg format -o GRCh38.sdf ~{ref}
        ## Run vcfeval
        rtg vcfeval --squash-ploidy -T ~{threadCount} -b child.vcf.gz -c father.vcf.gz -o child_not_father -t GRCh38.sdf        # false-negative variants -> present in child, but not in father
        rtg vcfeval --squash-ploidy -T ~{threadCount} -b child.vcf.gz -c mother.vcf.gz -o child_not_mother -t GRCh38.sdf        # false-negative variants -> present in child, but not in mother
        rtg vcfeval -T ~{threadCount} -b child_not_father/fn.vcf.gz -c child_not_mother/fn.vcf.gz -o denovos_out -t GRCh38.sdf  # true-positive variants -> present in child, but not in father and mother (denovos)

        ## Post-filtering of denovos by variant quality threshold
        zcat denovos_out/tp.vcf.gz | bcftools filter -e 'QUAL<20' -Oz -o ~{family}.denovos.vcf.gz

    >>>

	output {
		File vcf = "~{family}.denovos.vcf.gz"
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/shnegi/vcfeval@sha256:69e19ad005782b31767e75925e8f524fa0fe62ad1b65001f94a10bf3653e8fb0"
        preemptible: 1
    }
}

task subset_denovos_by_region {
    input {
        File input_vcf
        File lc_bed
        Int memSizeGB = 48
        Int threadCount = 16
        Int diskSizeGB = round(5*(size(input_vcf, "GB") + size(lc_bed, "GB"))) + 20
    }

    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")

    command <<<
        set -eux -o pipefail

        ## Filter denovos inside low-complexity regions in GRCh38
        bcftools view --threads ~{threadCount} -T ^~{lc_bed} -Oz -o ~{basen}.subset.vcf.gz ~{input_vcf}
        
    >>>

    output {
		File vcf = "~{basen}.subset.vcf.gz"
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/bcftools@sha256:f3a74a67de12dc22094e299fbb3bcd172eb81cc6d3e25f4b13762e8f9a9e80aa"
        preemptible: 1
    }
}

task run_filtering {
    input {
        File input_vcf
        File mother_vcf
        File father_vcf
        Int memSizeGB = 48
        Int threadCount = 32
        Int diskSizeGB = 5*round((size(input_vcf, "GB") + size(father_vcf, "GB")) + size(mother_vcf, "GB")) + 20
    }

    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")

    command <<<
        set -eux -o pipefail
        
        ## Apply variant filters to non-svim_asm variants, to remove possible false-positive denovos.
        zcat ~{input_vcf} | bcftools view --threads ~{threadCount} -i '((FORMAT/GQ>=20) & (FORMAT/DP>=20) & (FORMAT/GT!="1/1"))' | bcftools sort -Oz -o ~{basen}.filter.vcf.gz
        tabix -p vcf ~{basen}.filter.vcf.gz
        ln -s ~{basen}.filter.vcf.gz filtered.vcf.gz
        ln -s ~{basen}.filter.vcf.gz.tbi filtered.vcf.gz.tbi

        ## Index parent's VCFs
        tabix -p vcf ~{mother_vcf}; tabix -p vcf ~{father_vcf}
        ln -s ~{mother_vcf} mother.vcf.gz
        ln -s ~{mother_vcf}.tbi mother.vcf.gz.tbi
        ln -s ~{father_vcf} father.vcf.gz
        ln -s ~{father_vcf}.tbi father.vcf.gz.tbi

        ## Filtering denovos not called in both parents
        bcftools annotate -a mother.vcf.gz -m -MOM_MISSING filtered.vcf.gz -o itm.vcf.gz -Oz; tabix -p vcf itm.vcf.gz
        bcftools annotate -a father.vcf.gz -m -DAD_MISSING itm.vcf.gz -o ~{basen}.filter.flagged.vcf.gz -Oz; rm -rf itm.vcf.gz*
        
        ## Keep denovos (snps and indels) matching the filtering pattern
        bcftools view -v snps -i '(INFO/MOM_MISSING==1 & INFO/DAD_MISSING==1)' ~{basen}.filter.flagged.vcf.gz | bcftools sort -Oz -o ~{basen}.filter.flagged.snps.vcf.gz

    >>>

    output {
		File vcf = "~{basen}.filter.vcf.gz"
        File flagged_snps = "~{basen}.filter.flagged.snps.vcf.gz"
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/bcftools@sha256:f3a74a67de12dc22094e299fbb3bcd172eb81cc6d3e25f4b13762e8f9a9e80aa"
        preemptible: 1
    }
}

task annotate_denovo_snps {
    input {
        File input_vcf
        File? gnomad_vcf
        File? gnomad_vcf_index
        Int memSizeGB = 8
        Int threadCount = 32
        Int diskSizeGB = 5*round(size(input_vcf, "GB") + size(gnomad_vcf, "GB") + size(gnomad_vcf_index, "GB")) + 20
    }

    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")

    command <<<
    set -eux -o pipefail

    ## Index input VCF before transferring annotation
    tabix -p vcf ~{input_vcf}
    
    ## link the database VCF to make sure their indexes can be found
    ln -s ~{input_vcf} input.vcf.gz
    ln -s  ~{input_vcf}.tbi input.vcf.gz.tbi
    ln -s ~{gnomad_vcf} gnomad.vcf.bgz
    ln -s ~{gnomad_vcf_index} gnomad.vcf.bgz.tbi

    ## Annotate de novos with frequency in gnomAD
    bcftools annotate --threads ~{threadCount} -a gnomad.vcf.bgz -c AF input.vcf.gz -o ~{basen}.annotated.vcf.gz -Oz
    zcat ~{basen}.annotated.vcf.gz | bcftools view -i 'AF<0.001' | bcftools sort -o ~{basen}.rare.vcf.gz -Oz
    
    >>>

    output {
        File vcf = "~{basen}.rare.vcf.gz"
        File annotated_vcf = "~{basen}.annotated.vcf.gz"
    }
    
    runtime {
        memory: memSizeGB + " GB"
        cpu: 1
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/bcftools@sha256:f3a74a67de12dc22094e299fbb3bcd172eb81cc6d3e25f4b13762e8f9a9e80aa"
        preemptible: 1
    }
}

task annotate_with_snpeff {
    input {
        File input_vcf
        File snpeff_db
        String db_name
        Int memSizeGB = 16
        Int threadCount = 2
        Int diskSizeGB = 5*round(size(input_vcf, "GB") + size(snpeff_db, 'GB')) + 20
    }

    Int snpeffMem = if memSizeGB < 6 then 2 else memSizeGB - 4
    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")
    
    command <<<
    set -eux -o pipefail
    
    unzip ~{snpeff_db}
    
    zcat ~{input_vcf} | snpEff -Xmx~{snpeffMem}g -nodownload -no-intergenic \
                               -dataDir "${PWD}/data" ~{db_name} | gzip > ~{basen}.snpeff.vcf.gz
    >>>
    
    output {
        File vcf = "~{basen}.snpeff.vcf.gz"
    }
    
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/snpeff@sha256:7ac091da707f5d63f307eef4ee57c3f0e94eed49f86bbdace3d4be3a514ed410"
        preemptible: 1
    }
}

task subset_annotate_smallvars_with_db {
    input {
        File input_vcf
        File clinvar_vcf
        File clinvar_vcf_index
        File dbnsfp_db
        File dbnsfp_db_index
        Int memSizeGB = 16
        Int threadCount = 2
        Int diskSizeGB = 5*round(size(clinvar_vcf, 'GB') + size(dbnsfp_db, 'GB')) + 30
    }

    Int snpsiftMem = if memSizeGB < 6 then 2 else memSizeGB - 4
    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")
    
    command <<<
        set -eux -o pipefail

        ## link the database VCF to make sure their indexes can be found
        ln -s ~{clinvar_vcf} clinvar.vcf.bgz
        ln -s ~{clinvar_vcf_index} clinvar.vcf.bgz.tbi
        ln -s ~{dbnsfp_db} dbnsfp.txt.gz
        ln -s ~{dbnsfp_db_index} dbnsfp.txt.gz.tbi
        
        ## filter variants to keep those with high/moderate impact or with predicted loss of function
        # zcat ~{input_vcf} | SnpSift -Xmx1g filter "(ANN[*].IMPACT has 'HIGH') | (ANN[*].IMPACT has 'MODERATE') | ((exists LOF[*].PERC) & (LOF[*].PERC > 0.9))" | gzip  > smallvars.vcf.gz

        ## annotate IDs with clinvar IDs and add the CLNSIG INFO field
        zcat ~{input_vcf} | SnpSift -Xmx~{snpsiftMem}g annotate -info CLNSIG -v clinvar.vcf.bgz | gzip > smallvars.clinvar.vcf.gz
        
        ## annotate IDs with dbNSFP prediction scores and conservation scores
        zcat smallvars.clinvar.vcf.gz > smallvars.clinvar.vcf
        SnpSift -Xmx~{snpsiftMem}g dbnsfp -v -db dbnsfp.txt.gz -f GERP++_RS,CADD_raw,CADD_phred,MetaRNN_score,MetaRNN_pred,ALFA_Total_AF smallvars.clinvar.vcf | gzip > ~{basen}.clinvar.dbnsfp.vcf.gz
    >>>

    output {
        File vcf = "~{basen}.clinvar.dbnsfp.vcf.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/snpsift@sha256:049babfac841d15a92d8febfc10a25f5aa109c9fe6670af35ea79583a1c78402"
        preemptible: 1
    }
}

task validate_denovos {
    input {
        File input_vcf
        File mom_bam
        File mom_bam_index
        File dad_bam
        File dad_bam_index
        File reference_fasta
        Int memSizeGB = 30
        Int diskSizeGB = 5*round(size(input_vcf, "GB") + size(mom_bam, 'GB') + size(dad_bam, 'GB') + size(reference_fasta, 'GB')) + 30
    }

    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")

    command <<<
        set -eux -o pipefail

        ## link BAM files and indexes to make sure the indexes are found
        ln -s ~{mom_bam} mom.bam
        ln -s ~{mom_bam_index} mom.bam.bai
        ln -s ~{dad_bam} dad.bam
        ln -s ~{dad_bam_index} dad.bam.bai

        ## Validate and flag true positive de novos
        python3 /opt/scripts/dnvval.py -v ~{input_vcf} -mbam mom.bam -dbam dad.bam -r ~{reference_fasta} -o ~{basen}.dvval.vcf
        cat ~{basen}.dvval.vcf | bgzip > ~{basen}.dvval.vcf.gz

    >>>

    output {
        File vcf = "~{basen}.dnvval.vcf.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/shnegi/dvval@sha256:89233c11e7bb20f16294a7e7d56ce366741cbdbd0166bdaf833661bdb397e46d"
        preemptible: 1
    }


}