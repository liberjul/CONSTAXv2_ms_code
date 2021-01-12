#!/bin/bash -login

### Script to download databases, partition into exclusive datasets, and perform tests
### Assumes conda-installed CONSTAX
### Preferred memory allotment is 128 GB for SILVA, which likely requires a HPC cluster

### Environment variables
DIR=$(conda list | head -n 1 | rev | cut -d" " -f1 | rev | cut -d: -f1)
VERSION=2.0.3
BUILD=0
CONSTAXPATH=$DIR"/pkgs/constax-$VERSION-$BUILD/opt/constax-$VERSION"

### SILVA Bacteria
# Download data from SILVA and extract
SILVA_DATA=../../data/classification/silva
curl https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_tax_silva.fasta.gz > $SILVA_DATA/SILVA_138_SSURef_tax_silva.fasta.gz
gunzip $SILVA_DATA/SILVA_138_SSURef_tax_silva.fasta.gz

# Select records for Bacteria
python $CONSTAXPATH/fasta_select_by_keyword.py -i $SILVA_DATA/SILVA_138_SSURef_tax_silva.fasta -o $SILVA_DATA/SILVA_138_SSURef_tax_silva_bact.fasta -k "Bacteria;"
mkdir $SILVA_DATA/SILVA_Bacteria_tf
mv $SILVA_DATA/SILVA_138_SSURef_tax_silva_bact.fasta $SILVA_DATA/SILVA_Bacteria_tf/SILVA_138_SSURef_tax_silva_bact.fasta
# Create the formatted database files
python $CONSTAXPATH/FormatRefDB.py -d $SILVA_DATA/SILVA_Bacteria_tf/SILVA_138_SSURef_tax_silva_bact.fasta -t $SILVA_DATA/SILVA_Bacteria_tf -f SILVA -p $CONSTAXPATH

mkdir silva_parts
cd silva_parts

# Partition
sed 's/\t/;/' $SILVA_DATA/SILVA_Bacteria_tf/SILVA_138_SSURef_tax_silva_bact__RDP_taxonomy_headers.txt > $SILVA_DATA/SILVA_138_SSURef_tax_silva_bact__RDP_taxonomy_semis.txt
cut -d";" -f1-8 $SILVA_DATA/SILVA_138_SSURef_tax_silva_bact__RDP_taxonomy_semis.txt > $SILVA_DATA/SILVA_138_SSURef_tax_silva_bact__RDP_taxonomy_semis_cut.txt
for k in {0..4}
do # Partition according to Edgar 2016
  python partition.py -t $SILVA_DATA/SILVA_138_SSURef_tax_silva_bact__RDP_taxonomy_semis_cut.txt \
  -d $SILVA_DATA/SILVA_Bacteria_tf/SILVA_138_SSURef_tax_silva_bact.fasta \
  -f SILVA \
  -q $SILVA_DATA/silva_partition_"$k"_query \
  -r $SILVA_DATA/silva_partition_"$k"_ref

  # Subsample, because SILVA is too big to reasonably use full database
  python sub_sample_queries.py -d $SILVA_DATA/silva_partition_"$k"_query_fam.fasta -q $SILVA_DATA/silva_partition_"$k"_query_sub_fam.fasta -f SILVA -n 1000
  python sub_sample_queries.py -d $SILVA_DATA/silva_partition_"$k"_query_gen.fasta -q $SILVA_DATA/silva_partition_"$k"_query_sub_gen.fasta -f SILVA -n 1000
done
# Create simulated amplicons using in-silico PCR
for k in {0..4}
do
  for r in gen fam
  do
    python primer_matching.py -i $SILVA_DATA/silva_partition_"$k"_query_sub_1k_"$r".fasta \
      -o $SILVA_DATA/silva_partition_"$k"_query_"$r"_v4.fasta \
      -f GTGYCAGCMGCCGCGGTAA -r GGACTACNVGGGTWTCTAAT -m 3 # updated 515f 806R

    python primer_matching.py -i $SILVA_DATA/silva_partition_"$k"_query_sub_1k_"$r".fasta \
      -o $SILVA_DATA/silva_partition_"$k"_query_"$r"_v3-4.fasta \
      -f CCTACGGGNGGCWGCAG -r GACTACHVGGGTWTCTAAT -m 3 # updated 357wF and 785R from Herlemann et al. ISME 2011
  done
done
# Classify
for k in {0..4}
do
  for r in gen fam
  do
    # full length
    $CONSTAXPATH/constax.sh -c 0.8 -b -t -i $SILVA_DATA/silva_partition_"$k"_query_sub_1k_"$r".fasta \
      -n 16 -d $SILVA_DATA/silva_partition_"$k"_ref_"$r".fasta \
      -f $SILVA_DATA/tf_p"$k"_"$r"_1k -x $SILVA_DATA/tax_p"$k"_"$r"_1k -o $SILVA_DATA/out_p"$k"_"$r"_1k_cons --mem 128000 -m 20 \
      --conservative # conservative voting rule preferred for SILVA database
    # v4 region
    $CONSTAXPATH/constax.sh -c 0.8 -b -i $SILVA_DATA/silva_partition_"$k"_query_"$r"_v4.fasta \
      -n 16 -d $SILVA_DATA/silva_partition_"$k"_ref_"$r".fasta \
      -f $SILVA_DATA/tf_p"$k"_"$r"_1k -x $SILVA_DATA/tax_p"$k"_"$r"_v4 -o $SILVA_DATA/out_p"$k"_"$r"_v4_cons --mem 128000 -m 20 \
      --conservative
    # v3-4 regions
    $CONSTAXPATH/constax.sh -c 0.8 -b -i $SILVA_DATA/silva_partition_"$k"_query_"$r"_v3-4.fasta \
      -n 16 -d $SILVA_DATA/silva_partition_"$k"_ref_"$r".fasta \
      -f $SILVA_DATA/tf_p"$k"_"$r"_1k -x $SILVA_DATA/tax_p"$k"_"$r"_v3-4 -o $SILVA_DATA/out_p"$k"_"$r"_v3-4_cons --mem 128000 -m 20 \
      --conservative
    CONF=0.8; EVALUE=1; MAX_HITS=20; P_IDEN=0.8; FORMAT=SILVA
    for i in "1k" "v4" "v3-4"
    do
      OUTPUT=$SILVA_DATA/out_p"$k"_"$r"_"$i"
      mkdir $OUTPUT
      TAX=$SILVA_DATA/tax_p"$k"_"$r"_"$i"_cons
      TFILES=$SILVA_DATA/tf_p"$k"_"$r"_1k
      DB=$SILVA_DATA/silva_partition_"$k"_ref_"$r".fasta
      python $CONSTAXPATH/CombineTaxonomy.py -c $CONF -o "$OUTPUT/" -x "$TAX/" -b -e $EVALUE -m $MAX_HITS -p $P_IDEN -f $FORMAT -d $DB -t $TFILES -s False
    done
  done
done

# Parameter optimization for confidence threshold
for conf in 0.2 0.4 0.6 0.8 1.0
do
  for r in gen fam
  do
    for k in {0..4}
    do
      echo silva_partition_"$k"
      $CONSTAXPATH/constax.sh -c $conf -b -i $SILVA_DATA/silva_partition_"$k"_query_sub_1k_"$r".fasta \
        -n 16 --mem 10000 \
        -d $SILVA_DATA/silva_partition_"$k"_ref_"$r".fasta \
        -f $SILVA_DATA/tf_p"$k"_"$r"_1k -x $SILVA_DATA/tax_param -o $SILVA_DATA/out_param -m 20 --conservative
      mv $SILVA_DATA/out_param/combined_taxonomy.txt $SILVA_DATA/out_param/combined_taxonomy_"$r"_p"$k"_c"$conf"_m20_cons.txt
      mv $SILVA_DATA/out_param/consensus_taxonomy.txt $SILVA_DATA/out_param/consensus_taxonomy_"$r"_p"$k"_c"$conf"_m20_cons.txt
    done
  done
done
# Parameter optimization for max_hits parameter
for m in 1 3 5 10
do
  for r in gen fam
  do
    for k in {0..4}
    do
      echo silva_partition_"$k"
      $CONSTAXPATH/constax.sh -c 0.8 -b -i $SILVA_DATA/silva_partition_"$k"_query_sub_1k_"$r".fasta \
        -n 16 --mem 10000 \
        -d $SILVA_DATA/silva_partition_"$k"_ref_"$r".fasta \
        -f $SILVA_DATA/tf_p"$k"_"$r"_1k -x $SILVA_DATA/tax_param_m -o $SILVA_DATA/out_param_m --mem -m $m \
        --conservative
      mv $SILVA_DATA/out_param_m/combined_taxonomy.txt $SILVA_DATA/out_param/combined_taxonomy_"$r"_p"$k"_c0.8_m"$m"_cons.txt
      mv $SILVA_DATA/out_param_m/consensus_taxonomy.txt $SILVA_DATA/out_param/consensus_taxonomy_"$r"_p"$k"_c0.8_m"$m"_cons.txt
    done
  done
done

# Create true labels data for classification tests
mkdir query_dbs
for k in {0..4}
do
  for r in gen fam
  do
    python $CONSTAXPATH/FormatRefDB.py -d $SILVA_DATA/silva_partition_"$k"_query_sub_1k_"$r".fasta -t $SILVA_DATA/query_dbs -f SILVA -p $CONSTAXPATH
    for i in "v4" "v3-4"
    do
      python $CONSTAXPATH/FormatRefDB.py -d $SILVA_DATA/silva_partition_"$k"_query_"$r"_"$i".fasta -t $SILVA_DATA/query_dbs -f SILVA -p $CONSTAXPATH
    done
  done
done
cd ..

# Calculate classification performances
python get_metrics_silva_reg.py -d $SILVA_DATA
python get_metrics_silva_params.py -d $SILVA_DATA

# Run mothur for SILVA
for k in {0..4}
do
  for r in gen fam
  do
    grep "^>" $SILVA_DATA/silva_partition_"$k"_ref_"$r".fasta | sed 's/>//' | sed 's/ Bact/\tBact/' | sed 's/ /_/g' > \
      $SILVA_DATA/silva_partition_"$k"_ref_"$r".tax
  done
done
mothur mothur_class_seqs_silva_wang_knn_reg.txt
python tax_to_df_mothur_silva.py -d $SILVA_DATA
python tax_df_to_metrics_mothur_silva.py -d $SILVA_DATA

# Classification counts
$CONSTAXPATH/constax.sh -c 0.8 -b -t \
  -i  $SILVA_DATA/otus_16s.fasta \
  -n 16 \
  -d $SILVA_DATA/SILVA_Bacteria_tf/SILVA_138_SSURef_tax_silva_bact.fasta \
  -f $SILVA_DATA/SILVA_Bacteria_tf \
  -x $SILVA_DATA/tax_sil -o $SILVA_DATA/out_sil --mem 128000 -m 20 \
  --conservative
mv $SILVA_DATA/out_sil/combined_taxonomy.txt $SILVA_DATA/out_sil/combined_taxonomy_sil_blast.txt

### UNITE Fungi
UNITE_DATA=../../data/classification/unite
#Download and extract from UNITE
curl https://files.plutof.ut.ee/public/orig/E7/28/E728E2CAB797C90A01CD271118F574B8B7D0DAEAB7E81193EB89A2AC769A0896.gz > $UNITE_DATA/sh_general_release_04.02.2020.tar.gz
tar -xzvf $UNITE_DATA/sh_general_release_04.02.2020.tar.gz
mkdir $UNITE_DATA/UNITE_Fungi_tf
mv $UNITE_DATA/sh_general_release_04.02.2020/sh_general_release_dynamic_04.02.2020.fasta $UNITE_DATA/UNITE_Fungi_tf/sh_general_release_dynamic_04.02.2020.fasta
python $CONSTAXPATH/FormatRefDB.py -d $UNITE_DATA/UNITE_Fungi_tf/sh_general_release_dynamic_04.02.2020.fasta -t $UNITE_DATA/UNITE_Fungi_tf -f UNITE -p $CONSTAXPATH

sed 's/\t/;/' $UNITE_DATA/UNITE_Fungi_tf/sh_general_release_dynamic_04.02.2020__RDP_taxonomy_headers.txt > $UNITE_DATA/sh_general_release_dynamic_04.02.2020__RDP_taxonomy_semis.txt
cut -d';' -f1-8 $UNITE_DATA/sh_general_release_dynamic_04.02.2020__RDP_taxonomy_semis.txt > $UNITE_DATA/sh_general_release_dynamic_04.02.2020__RDP_taxonomy_semis_cut.txt
for k in {0..4}
do
  echo unite_partition_"$k"
  python partition.py -t $UNITE_DATA/sh_general_release_dynamic_04.02.2020__RDP_taxonomy_semis_cut.txt \
  -d $UNITE_DATA/UNITE_Fungi_tf/sh_general_release_dynamic_04.02.2020.fasta \
  -f UNITE \
  -q $UNITE_DATA/unite_partition_"$k"_query \
  -r $UNITE_DATA/unite_partition_"$k"_ref
done

# Use ITSx to extract its1 and 2 regions. Will not capture all, but works better than in-silico PCR
conda install -c bioconda itsx
for k in {0..4}
do
  for r in gen fam
  do
    ITSx -i $UNITE_DATA/unite_partition_"$k"_query_"$i".fasta \
      -o $UNITE_DATA/unite_partition_"$k"_query_"$i"_itsx \
      --cpu 16
  done
done

for k in {0..4}
do
  for r in gen fam
  do
    $CONSTAXPATH/constax.sh -c 0.8 -b -t -i $UNITE_DATA/unite_partition_"$k"_query_"$r".fasta \
      -n 16 -d $UNITE_DATA/unite_partition_"$k"_ref_"$r".fasta \
      -f $UNITE_DATA/tf_p"$k"_"$r" -x $UNITE_DATA/tax_p"$k"_"$r" -o $UNITE_DATA/out_p"$k"_"$r" --mem 32000 -m 5
    $CONSTAXPATH/constax.sh -c 0.8 -b -i $UNITE_DATA/unite_partition_"$k"_query_"$r"_itsx.ITS1.fasta \
      -n 16 -d $UNITE_DATA/unite_partition_"$k"_ref_"$r".fasta \
      -f $UNITE_DATA/tf_p"$k"_"$r" -x $UNITE_DATA/tax_p"$k"_"$r"_its1 -o $UNITE_DATA/out_p"$k"_"$r"_its1 --mem 32000 -m 5
    $CONSTAXPATH/constax.sh -c 0.8 -b -i $UNITE_DATA/unite_partition_"$k"_query_"$r".fasta \
      -n 16 -d $UNITE_DATA/unite_partition_"$k"_ref_"$r".fasta \
      -f $UNITE_DATA/tf_p"$k"_"$r" -x $UNITE_DATA/tax_p"$k"_"$r"_its2 -o $UNITE_DATA/out_p"$k"_"$r"_its2 --mem 32000 -m 5

    $CONSTAXPATH/constax.sh -c 0.8 -t -i $UNITE_DATA/unite_partition_"$k"_query_"$r".fasta \
      -n 16 -d $UNITE_DATA/unite_partition_"$k"_ref_"$r".fasta \
      -f $UNITE_DATA/tf_p"$k"_"$r"_u -x $UNITE_DATA/tax_p"$k"_"$r"_u -o $UNITE_DATA/out_p"$k"_"$r"_u --mem 32000 -m 5
    $CONSTAXPATH/constax.sh -c 0.8 -i $UNITE_DATA/unite_partition_"$k"_query_"$r"_itsx.ITS1.fasta \
      -n 16 -d $UNITE_DATA/unite_partition_"$k"_ref_"$r".fasta \
      -f $UNITE_DATA/tf_p"$k"_"$r"_u -x $UNITE_DATA/tax_p"$k"_"$r"_its1_u -o $UNITE_DATA/out_p"$k"_"$r"_its1_u --mem 32000 -m 5
    $CONSTAXPATH/constax.sh -c 0.8 -i $UNITE_DATA/unite_partition_"$k"_query_"$r".fasta \
      -n 16 -d $UNITE_DATA/unite_partition_"$k"_ref_"$r".fasta \
      -f $UNITE_DATA/tf_p"$k"_"$r"_u -x $UNITE_DATA/tax_p"$k"_"$r"_its2_u -o $UNITE_DATA/out_p"$k"_"$r"_its2_u --mem 32000 -m 5

    CONF=0.8; EVALUE=1; MAX_HITS=5; P_IDEN=0.8; FORMAT=UNITE
    for i in "" "_its1" "_its2"
    do
      OUTPUT=$UNITE_DATA/out_p"$k"_"$r"_u"$i"_cons
      mkdir $OUTPUT
      TAX=$UNITE_DATA/tax_p"$k"_"$r"_u"$i"
      TFILES=$UNITE_DATA/tf_p"$k"_"$r"_u
      DB=$UNITE_DATA/unite_partition_"$k"_ref_"$r".fasta
      python $CONSTAXPATH/CombineTaxonomy.py -c $CONF -o "$OUTPUT/" -x "$TAX/" -e $EVALUE -m $MAX_HITS -p $P_IDEN -f $FORMAT -d $DB -t $TFILES -s True
    done
  done
done

for conf in 0.2 0.4 0.6 0.8 1.0
do
  for r in gen fam
  do
    for k in {0..4}
    do
      for m in 5 20
      do
        echo unite_partition_"$k"
        /mnt/research/bonito_lab/CONSTAX_May2020/constax.sh -c $conf -b -i unite_partition_"$k"_query_sub_1k_"$r".fasta \
          -n $SLURM_CPUS_PER_TASK \
          -d unite_partition_"$k"_ref_"$r".fasta \
          -f tf_p"$k"_"$r" -x tax_p"$k"_"$r" -o out_param --mem $SLURM_MEM_PER_NODE -m $m
        mv out_param/combined_taxonomy.txt out_param/combined_taxonomy_"$r"_p"$k"_c"$conf"_m"$m".txt
        mv out_param/consensus_taxonomy.txt out_param/consensus_taxonomy_"$r"_p"$k"_c"$conf"_m"$m".txt
      done
    done
  done
done

for m in 1 3 10
do
  for r in gen fam
  do
    for k in {0..4}
    do
      echo unite_partition_"$k"
      $CONSTAXPATH/constax.sh -c 0.8 -b -i $UNITE_DATA/unite_partition_"$k"_query_sub_1k_"$r".fasta \
        -n $SLURM_CPUS_PER_TASK \
        -d $UNITE_DATA/unite_partition_"$k"_ref_"$r".fasta \
        -f $UNITE_DATA/tf_p"$k"_"$r" -x $UNITE_DATA/tax_p"$k"_"$r" -o $UNITE_DATA/out_param --mem $SLURM_MEM_PER_NODE -m $m
      mv $UNITE_DATA/out_param/combined_taxonomy.txt out_param/combined_taxonomy_"$r"_p"$k"_c0.8_m"$m".txt
      mv $UNITE_DATA/out_param/consensus_taxonomy.txt out_param/consensus_taxonomy_"$r"_p"$k"_c0.8_m"$m".txt
    done
  done
done

mkdir $UNITE_DATA/query_dbs
for k in {0..4}
do
  for r in gen fam
  do
    python $CONSTAXPATH/FormatRefDB.py -d $UNITE_DATA/unite_partition_"$k"_query_gen.fasta -t $UNITE_DATA/query_dbs -f UNITE -p $CONSTAXPATH
    for i in 1 2
    do
      python $CONSTAXPATH/FormatRefDB.py -d $UNITE_DATA/unite_partition_"$k"_query_"$r"_itsx.ITS"$i".fasta -t $UNITE_DATA/query_dbs -f UNITE -p $CONSTAXPATH
    done
  done
done

python get_metrics_unite_reg_blast_cons.py $UNITE_DATA
python get_metrics_unite_reg_utax_cons.py $UNITE_DATA
python get_metrics_unite_params.py $UNITE_DATA

# Run mothur for UNITE
# Download mothur release from UNITE
curl https://files.plutof.ut.ee/public/orig/B2/A2/B2A239FA96894583366BA7958877A59B62561C448F7D4E74F1BC50C0B5B8FD73.gz > $UNITE_DATA/sh_mothur_release_04.02.2020.tar.gz
tar -xzvf $UNITE_DATA/sh_mothur_release_04.02.2020.tar.gz
# Create partition files in mothur format
for i in {0..4}
do
  for r in gen fam
  do
    python to_mothur_format.py -i $UNITE_DATA/unite_partition_"$i"_ref_"$r".fasta \
      --itax $UNITE_DATA/sh_mothur_release_04.02.2020/UNITEv8_sh_99.tax --otax $UNITE_DATA/uni_parts_"$i"_"$r".tax --ialign $UNITE_DATA/sh_mothur_release_04.02.2020/UNITEv8_sh_99.fasta --oalign $UNITE_DATA/uni_parts_"$i"_"$r".align
  done
done
# Classify with mothur
mothur mothur_class_seqs_unite_wang_knn_its.txt
python tax_to_df_mothur_unite.py -d $UNITE_DATA
python tax_df_to_metrics_mothur_unite.py -d $UNITE_DATA

# Classification counts
$CONSTAXPATH/constax.sh -c 0.8 -b -t \
  -i  $UNITE_DATA/otus_ITS.fasta \
  -n 16 \
  -d $UNITE_DATA/UNITE_Fungi_tf/sh_general_release_dynamic_04.02.2020.fasta \
  -f $UNITE_DATA/UNITE_Fungi_tf \
  -x $UNITE_DATA/tax_uni -o $UNITE_DATA/out_uni --mem 32000 -m 5
mv $UNITE_DATA/out_uni/combined_taxonomy.txt $UNITE_DATA/out_uni/combined_taxonomy_uni_blast.txt
