#!/bin/bash
set -e

ENV=envs/bleties_env

case "$1" in

  "get_software")
    mkdir -p opt
    cd opt
    # faFilter from KentUtils suite
    wget --no-clobber http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faFilter
    wget --no-clobber http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSomeRecords
    # ART read simulator for short reads
    wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
    tar -xzf artbinmountrainier2016.06.05linux64.tgz
    # PBSim2 read simulator for long reads
    git clone git@github.com:yukiteruono/pbsim2.git # commit e71f7892
    # Minimap2 mapper
    git clone https://github.com/lh3/minimap2
    # Seqtk
    git clone git@github.com:lh3/seqtk.git
    ;;

  "download_data")
    mkdir -p ref
    cd ref
    wget --no-clobber https://paramecium.i2bc.paris-saclay.fr/files/Paramecium/tetraurelia/51/sequences/ptetraurelia_mac_51.fa
    wget --no-clobber https://paramecium.i2bc.paris-saclay.fr/files/Paramecium/tetraurelia/51/sequences/ptetraurelia_mac_51_with_ies.fa
    ;;

  "filter_long_contigs")
    # Filter long contigs in MAC+IES assembly
    opt/faFilter -minSize=100000 -maxN=100 ref/ptetraurelia_mac_51_with_ies.fa ref/ptetraurelia_mac_51_with_ies.min100k_max100N.fa
    # Find corresponding contigs in the MAC-only assembly
    grep '>' ref/ptetraurelia_mac_51_with_ies.min100k_max100N.fa | sed 's/>//' | sed 's/_with_IES//' > ref/ctg_list
    opt/faSomeRecords ref/ptetraurelia_mac_51_with_ies.min100k_max100N.fa ref/ctg_list ref/ptetraurelia_mac_51.min100k_max100N.fa
    ;;

  "sim_pb_clr")
    # MAC reads
    opt/pbsim2/bin/pbsim \
      --prefix sim_pb_clr/ptet_mac_ies_sim_clr \
      --id-prefix ptet_mac_ies_sim_clr \
      --depth 50 \
      --length-min 100 \
      --length-max 1000000 \
      --seed 12345 \
      --hmm_model opt/pbsim2/data/P6C4.model \
      --length-mean 9000 \
      --length-sd 7000 \
      ref/ptetraurelia_mac_51_with_ies.min100k_max100N.fa &> sim_pb_clr/pbsim2_mac_ies.log
    # MAC+IES reads
    opt/pbsim2/bin/pbsim \
      --prefix sim_pb_clr/ptet_mac_sim_clr \
      --id-prefix ptet_mac_sim_clr \
      --depth 50 \
      --length-min 100 \
      --length-max 1000000 \
      --seed 12345 \
      --hmm_model opt/pbsim2/data/P6C4.model \
      --length-mean 9000 \
      --length-sd 7000 \
      ref/ptetraurelia_mac_51.min100k_max100N.fa &> sim_pb_clr/pbsim2_mac.log
    # Combine and gzip
    cat sim_pb_clr/ptet_mac_ies_sim_clr_*.fastq | gzip > sim_pb_clr/ptet_mac_ies_sim_clr.fq.gz
    cat sim_pb_clr/ptet_mac_sim_clr_*.fastq | gzip > sim_pb_clr/ptet_mac_sim_clr.fq.gz

  "mix_sim_pb_clr_reads")
    # downsample, rename, and combine reads
    # 5:45 mixture of MAC+IES:MAC
    opt/seqtk/seqtk sample -s 12345 sim_pb_clr/ptet_mac_ies_sim_clr.fq.gz 0.1 > sim_pb_clr/ptet_sim_clr.mix05.fq
    opt/seqtk/seqtk sample -s 12345 sim_pb_clr/ptet_mac_sim_clr.fq.gz 0.9 >> sim_pb_clr/ptet_sim_clr.mix05.fq
    gzip sim_pb_clr/ptet_sim_clr.mix05.fq
    # 10:40 mixture of MAC+IES:MAC
    opt/seqtk/seqtk sample -s 12345 sim_pb_clr/ptet_mac_ies_sim_clr.fq.gz 0.2 > sim_pb_clr/ptet_sim_clr.mix10.fq
    opt/seqtk/seqtk sample -s 12345 sim_pb_clr/ptet_mac_sim_clr.fq.gz 0.8 >> sim_pb_clr/ptet_sim_clr.mix10.fq
    gzip sim_pb_clr/ptet_sim_clr.mix10.fq
    # 20:30 mixture of MAC+IES:MAC
    opt/seqtk/seqtk sample -s 12345 sim_pb_clr/ptet_mac_ies_sim_clr.fq.gz 0.4 > sim_pb_clr/ptet_sim_clr.mix20.fq
    opt/seqtk/seqtk sample -s 12345 sim_pb_clr/ptet_mac_sim_clr.fq.gz 0.6 >> sim_pb_clr/ptet_sim_clr.mix20.fq
    gzip sim_pb_clr/ptet_sim_clr.mix20.fq
    # 40:10 mixture of MAC+IES:MAC
    opt/seqtk/seqtk sample -s 12345 sim_pb_clr/ptet_mac_ies_sim_clr.fq.gz 0.8 > sim_pb_clr/ptet_sim_clr.mix40.fq
    opt/seqtk/seqtk sample -s 12345 sim_pb_clr/ptet_mac_sim_clr.fq.gz 0.2 >> sim_pb_clr/ptet_sim_clr.mix40.fq
    gzip sim_pb_clr/ptet_sim_clr.mix40.fq
    # alias MAC-only library as mix00
    ln --relative -s sim_pb_clr/ptet_mac_sim_clr.fq.gz sim_pb_clr/ptet_sim_clr.mix00.fq.gz
    ;;

  "sim_ont"
    # MAC reads
    opt/pbsim2/bin/pbsim \
      --prefix sim_ont/ptet_mac_sim_ontR95 \
      --id-prefix ptet_mac_sim_ontR95 \
      --depth 50 \
      --length-min 100 \
      --length-max 1000000 \
      --seed 12345 \
      --difference-ratio 23:31:46 \
      --hmm_model opt/pbsim2/data/R95.model \
      --length-mean 9000 \
      --length-sd 7000 \
      ref/ptetraurelia_mac_51.min100k_max100N.fa &> sim_ont/pbsim2_mac_ontR95.log
    # MAC+IES reads
      opt/pbsim2/bin/pbsim \
        --prefix sim_ont/ptet_mac_ies_sim_ontR95 \
        --id-prefix ptet_mac_ies_sim_ontR95 \
        --depth 50 \
        --length-min 100 \
        --length-max 1000000 \
        --seed 12345 \
        --difference-ratio 23:31:46 \
        --hmm_model opt/pbsim2/data/R95.model \
        --length-mean 9000 \
        --length-sd 7000 \
        ref/ptetraurelia_mac_51_with_ies.min100k_max100N.fa &> sim_ont/pbsim2_mac_ies_ontR95.log
    # Combine and gzip
    cat sim_ont/ptet_mac_ies_sim_ontR95_*.fastq | gzip > sim_ont/ptet_mac_ies_sim_ontR95.fq
    cat sim_ont/ptet_mac_sim_ontR95_*.fastq | gzip > sim_ont/ptet_mac_sim_ontR95.fq
    ;;
  
  "mix_sim_ont_reads")
    # downsample, rename, and combine reads
    # 5:45 mixture of MAC+IES:MAC
    opt/seqtk/seqtk sample -s 12345 sim_ont/ptet_mac_ies_sim_ontR95.fq.gz 0.1 > sim_ont/ptet_sim_ontR95.mix05.fq
    opt/seqtk/seqtk sample -s 12345 sim_ont/ptet_mac_sim_ontR95.fq.gz 0.9 >> sim_ont/ptet_sim_ontR95.mix05.fq
    gzip sim_ont/ptet_sim_ontR95.mix05.fq
    # 10:40 mixture of MAC+IES:MAC
    opt/seqtk/seqtk sample -s 12345 sim_ont/ptet_mac_ies_sim_ontR95.fq.gz 0.2 > sim_ont/ptet_sim_ontR95.mix10.fq
    opt/seqtk/seqtk sample -s 12345 sim_ont/ptet_mac_sim_ontR95.fq.gz 0.8 >> sim_ont/ptet_sim_ontR95.mix10.fq
    gzip sim_ont/ptet_sim_ontR95.mix10.fq
    # 20:30 mixture of MAC+IES:MAC
    opt/seqtk/seqtk sample -s 12345 sim_ont/ptet_mac_ies_sim_ontR95.fq.gz 0.4 > sim_ont/ptet_sim_ontR95.mix20.fq
    opt/seqtk/seqtk sample -s 12345 sim_ont/ptet_mac_sim_ontR95.fq.gz 0.6 >> sim_ont/ptet_sim_ontR95.mix20.fq
    gzip sim_ont/ptet_sim_ontR95.mix20.fq
    # 40:10 mixture of MAC+IES:MAC
    opt/seqtk/seqtk sample -s 12345 sim_ont/ptet_mac_ies_sim_ontR95.fq.gz 0.8 > sim_ont/ptet_sim_ontR95.mix40.fq
    opt/seqtk/seqtk sample -s 12345 sim_ont/ptet_mac_sim_ontR95.fq.gz 0.2 >> sim_ont/ptet_sim_ontR95.mix40.fq
    gzip sim_ont/ptet_sim_ontR95.mix40.fq
    # alias MAC-only library as mix00
    ln --relative -s sim_ont/ptet_mac_sim_ontR95.fq.gz sim_ont/ptet_sim_ontR95.mix00.fq.gz
  ;;
  
  "map_pb")
    # PacBio CLR
    # Different coverage ratios of MAC+IES:MAC mixtures
    REF=ref/ptetraurelia_mac_51.min100k_max100N.fa
    for i in 00 05 10 20 40; do
        opt/minimap2/minimap2 -ax map-pb -t 16 \
          -o mapping_sim_pb/ptet_sim_clr.mix${i}.sam \
          $REF sim_pb_clr/ptet_sim_clr.mix${i}.fq.gz
        # convert SAM to BAM, sort and index
        samtools view -b -@16 --reference $REF \
          -o mapping_sim_pb/ptet_sim_clr.mix${i}.bam \
          mapping_sim_pb/ptet_sim_clr.mix${i}.sam
        samtools sort -@16 --reference $REF \
          -o mapping_sim_pb/ptet_sim_clr.mix${i}.sort.bam \
          mapping_sim_pb/ptet_sim_clr.mix${i}.bam
        samtools index -@4 mapping_sim_pb/ptet_sim_clr.mix${i}.sort.bam
    done
    ;;
  
  "map_ont")
    # ONT R9.5
    # Different coverage ratios of MAC+IES:MAC mixtures
    REF=ref/ptetraurelia_mac_51.min100k_max100N.fa
    for i in 00 05 10 20 40; do
        opt/minimap2/minimap2 -ax map-ont -t 16 \
          -o mapping_sim_ont/ptet_sim_ontR95.mix${i}.sam \
          $REF sim_ont/ptet_sim_ontR95.mix${i}.fq.gz
        # convert SAM to BAM, sort and index
        samtools view -b -@16 --reference $REF \
          -o mapping_sim_ont/ptet_sim_ontR95.mix${i}.bam \
          mapping_sim_ont/ptet_sim_ontR95.mix${i}.sam
        samtools sort -@16 --reference $REF \
          -o mapping_sim_ont/ptet_sim_ontR95.mix${i}.sort.bam \
          mapping_sim_ont/ptet_sim_ontR95.mix${i}.bam
        samtools index -@4 mapping_sim_ont/ptet_sim_ontR95.mix${i}.sort.bam
    done
    ;;

  "bleties_pb")
    source activate $ENV
    bleties --version
    mkdir -p bleties_pb
    for i in 00 05 10 20 40; do
      bleties --log bleties_pb/ptet_sim_clr.mix${i}.milraa_subreads.log \
        milraa \
        --bam mapping_sim_pb/ptet_sim_clr.mix${i}.sort.bam \
        --ref ref/ptetraurelia_mac_51.min100k_max100N.fa \
        --min_break_coverage 5 --min_del_coverage 5 \
        --type subreads \
        --dump --out bleties_pb/ptet_sim_clr.mix${i}.milraa_subreads
    done
    ;;

  "bleties_ont")
    source activate $ENV
    bleties --version
    mkdir -p bleties_ont
    for i in 00 05 10 20 40; do
      bleties --log bleties_ont/ptet_sim_ontR95.mix${i}.milraa_subreads.log \
        milraa \
        --bam mapping_sim_ont/ptet_sim_ontR95.mix${i}.sort.bam \
        --ref ref/ptetraurelia_mac_51.min100k_max100N.fa \
        --min_break_coverage 5 --min_del_coverage 5 \
        --type subreads \
        --dump --out bleties_ont/ptet_sim_ontR95.mix${i}.milraa_subreads
    done
    ;;

  *)
    echo "Subcommands"
    echo "  get_software"
    echo "  download_data"
    echo "  filter_long_contigs"
    echo "  sim_pb_clr"
    echo "  mix_sim_pb_clr_reads"
    echo "  sim_ont"
    echo "  mix_sim_ont_reads"
    echo "  map_pb"
    echo "  map_ont"
    echo "  bleties_pb"
    echo "  bleties_ont"
    exit 1

esac
