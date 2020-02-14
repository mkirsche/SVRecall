 javac -cp /work-zfs/mschatz1/pipelines/soft/git/Jasmine/src src/*.java
java -cp /work-zfs/mschatz1/pipelines/soft/git/Jasmine/src:src Rescreen merged_vcf=/scratch/groups/mschatz1/mkirsche/jasminepost/merged_clean_500_precise.vcf bam_file_list=/scratch/groups/mschatz1/mkirsche/jasminepost/bamfilelist.txt out_file=test.vcf sniffles_path=/work-zfs/mschatz1/pipelines/soft/sniffles samtools_path=/work-zfs/mschatz1/pipelines/soft/samtools | tee log.txt

