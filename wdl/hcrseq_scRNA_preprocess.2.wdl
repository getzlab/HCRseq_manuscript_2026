workflow hcrseq_scRNA_quantify{
  call preprocess
}

task preprocess{
  File h5_file
  File bam
  File bai
  File lesion_info
  File pathway_info
  String outstem
  Int? disk_GB = 200

  command {
    set -euo pipefail
	
    hcrseq scrna quantify --h5_file ${h5_file} \
			--bam ${bam} \
			--lesion_info ${lesion_info} \
			--pathway_info ${pathway_info} \
			--outstem ${outstem}
                           
  }


  output {
    File repair_annotated_h5ad="${outstem}.h5ad"

  }

  runtime {
    docker: "gcr.io/broad-getzlab-fmhcrsparc/hcrseq:v0.3"
    memory: "16GB"
    disks: "local-disk " + disk_GB + " HDD"
    cpu: "4"
    preemptible: "2"
  }


}