"""
Porting the wrapper script into Snakemake for debugging

=head1 SYNOPSIS

  wrapper_phage_contigs_sorter_iPlant.pl --fasta sequences.fa

Required Arguments:

  -f|--fna       Fasta file of contigs

Options:

  -d|--dataset   Code dataset (DEFAULT "VIRSorter")
  --cp           Custom phage sequence
  --db           Either "1" (DEFAULT Refseqdb) or "2" (Viromedb)
  --wdir         Working directory (DEFAULT cwd)
  --ncpu         Number of CPUs (default: 4)
  --virome       Add this flag to enable virome decontamination mode, for datasets
                 mostly viral to force the use of generic metrics instead of
                 calculated from the whole dataset. (default: off)
  --data-dir     Path to "virsorter-data" directory (DEFAULT /data/virsorter-data)
  --diamond      Use diamond (in "--more-sensitive" mode) instead of blastp.
                 Diamond is much faster than blastp and may be useful for adding
		 many custom phages, or for processing extremely large Fasta files.
		 Unless you specify --diamond, VirSorter will use blastp.
  --keep-db      Specifying this flag keeps the new HMM and BLAST databases created
                 after adding custom phages. This is useful if you have custom phages
		 that you want to be included in several different analyses and want
		 to save the database and point VirSorter to it in subsequent runs.
                 By default, this is off, and you should only specify this flag if
		 you're SURE you need it.
  --help         Show help and exit

"""
import os
from datetime import datetime
from Bio import SeqIO

start_time = datetime.now()

snakefile_path = os.path.dirname(os.path.abspath(workflow.snakefile))
input_file = config['input_file']
code_dataset = config.get('dataset', 'Virsorter')
choice_database = int(config.get('db', 1))
tag_virome = config.get('virome', False)
data_dir = config['data_dir']
wdir = config.get('wdir', '.')
diamond = config.get('diamond', False)
custom_phage = config.get('cp', None)
keepdb = config.get('keep_db', False)

if not os.path.exists(input_file):
    raise Exception("Missing input file: " + input_file)

if choice_database not in [1,2,3]:
    raise Exception('choice_database must be 1, 2, or 3, not "{}"'.format(choice_database))

print('\n'.join(map(lambda t: "%-15s: %s" % t,
   [
    ('Bin',           snakefile_path),
    ('Dataset',       code_dataset),
    ('Input file',    input_file),
    ('Db',            choice_database),
    ('Working dir',   wdir),
    ('Custom phages', custom_phage),
    ('Data dir',      data_dir),
    ('blastp',        'diamond' if diamond else 'blastp'),
])))

if diamond:
    print("This VirSorter run uses DIAMOND (Buchfink et al., Nature Methods 2015) instead of blastp.")
if tag_virome:
    print( "WARNING: THIS WILL BE A VIROME DECONTAMINATION RUN")

# Need 2 databases
# PCs from Refseq (phages) or PCs from Refseq+Viromes
# PFAM (27.0)

script_dir         = os.path.join(snakefile_path, 'Scripts')

if tag_virome:
    readme_file    = os.path.join(data_dir, 'VirSorter_Readme_viromes.txt')
else:
    readme_file    = os.path.join(data_dir, 'VirSorter_Readme.txt')

generic_ref_file   = os.path.join(data_dir, 'Generic_ref_file.refs')

if choice_database == 2:
    dir_Phage_genes    = os.path.join(data_dir, 'Phage_gene_catalog_plus_viromes')
    ref_phage_clusters = os.path.join(data_dir,
        'Phage_gene_catalog_plus_viromes', 'Phage_Clusters_current.tab')
elif choice_database == 3:
    dir_Phage_genes    = os.path.join(data_dir, 'euk-virus')
    # ??? what goes here?  I don't have this file (sic)
    ref_phage_clusters = os.path.join(data_dir,
        'euk-virus', 'Phage_Clusters_current.tab')
else:
    dir_Phage_genes    = os.path.join(data_dir,'Phage_gene_catalog')
    ref_phage_clusters = os.path.join(data_dir,
                             'Phage_gene_catalog', 'Phage_Clusters_current.tab')

db_PFAM_a = os.path.join(data_dir, 'PFAM_27', 'Pfam-A.hmm')
db_PFAM_b = os.path.join(data_dir, 'PFAM_27', 'Pfam-B.hmm')

## SETTING UP THE WORKING DIRECTORY
log_dir = os.path.join(wdir, 'logs')

wildcard_constraints:
    r=r'\d+',
    table_type=r'(clusters|unclustered)'

rule all:
    input:
        os.path.join(wdir, 'Readme.txt'),
        os.path.join(wdir, 'Predicted_viral_sequences.fasta')

# cp fasta file in the wdir (and rename)
fastadir = os.path.join(wdir, 'fasta')
rule rename_fasta:
    input: input_file
    output: os.path.join(wdir, 'fasta', 'input_sequences.fna')
    run:
        with open(output[0], 'wt') as output_handle:
            for record in SeqIO.parse(input[0], 'fasta'):
                new_id = code_dataset + "_" + record.id
                output_handle.write(">{}\n{}\n".format(new_id,
                                                       str(record.seq)))


# detect circular, predict genes on contigs and extract proteins, as well
# as filtering on size (nb genes) and/or circular
rule step_0_5:
    input: rules.rename_fasta.output
    output:
        fasta_contigs_nett=os.path.join(fastadir,
                                        code_dataset + "_nett_filtered.fasta"),
        fasta_file_prots=os.path.join(fastadir, code_dataset + "_prots.fasta"),
        predict_file=os.path.join(fastadir, code_dataset + '_mga_final.predict')
        # ignored outputs:
        # VST_nett.fasta
        # VST_circu.list
        # VST_special_circus_mga.predict
        # VST_mga.predict
        # VST_contigs_circu_temp.fasta
    log: os.path.join(log_dir, "step_0_5.out")
    params:
        # at least 2 genes on the contig
        nb_gene_th=config.get('nb_gene_th', 2),
        path_script_step_1=os.path.join(script_dir,
                    "Step_1_contigs_cleaning_and_gene_prediction.pl")
    shell:
        "{params.path_script_step_1} {code_dataset} {fastadir} {input} \
            {params.nb_gene_th} > {log} 2>&1"

# Match against PFAM, once for all
# compare to PFAM a then b (hmmsearch)
rule step_0_7:
    input: rules.step_0_5.output.fasta_file_prots
    output:
        hmmsearch_pfam=os.path.join(wdir, 'Contigs_prots_vs_PFAM{ab}.tab'),
        hmmsearch_pfam_bis=os.path.join(wdir, 'Contigs_prots_vs_PFAM{ab}.out')
    log: os.path.join(log_dir, "step_0_7.out")
    threads: config.get('threads_hmmer', 2)
    params:
        pfam=lambda w: db_PFAM_a if w.ab=='a' else db_PFAM_b
    shell:
        "hmmsearch --tblout {output.hmmsearch_pfam} --cpu {threads} \
         -o {output.hmmsearch_pfam_bis} --noali {params.pfam} \
         {input} > {log} 2>&1 "


# Snakemake can't do a dynamic DAG, so we'll do all 9 revisions, but
# short circuit rules once we're done
max_revisions = config.get('max_revisions', 10)

# First revision, we just import the Refseq database
if custom_phage:
    rule add_custom_phage:
        input: custom_phage
        output: os.path.join(wdir, 'r_0', 'db', 'Pool_unclustered.faa')
        log: os.path.join(log_dir, "add_custom_phage.out")
        threads: config.get('threads_blast', 10)
        shell:
            """{script_dir}/Step_first_add_custom_phage_sequence.pl {input} \
                {dir_Phage_genes} r_0/db {threads} > {log} 2>&1
            """
else:
    rule copy_std_phage:
        output: os.path.join(wdir, 'r_0', 'db', 'Pool_unclustered.faa')
        shell:
            """cp {dir_Phage_genes}/* {wdir}/r_0/db/
            """

rule initial_blank_tab:
    output: os.path.join(wdir, 'r_0', 'Contigs_prots_vs_Phage_{table_type}.tab')
    shell: "touch {output}"

## Clustering of the new prots with the unclustered
script_new_cluster = os.path.join(script_dir, "Step_0_make_new_clusters.pl")

rule step_1_1:
    input:
        prots=rules.step_0_5.output.fasta_file_prots,
        prev_unclusterd=lambda w: os.path.join(wdir, \
                                               'r_' + str(int(w.r)-1), \
                                               'db', \
                                               'Pool_unclustered.faa'),
        prots_to_cluster=lambda w: os.path.join(wdir, \
                                           'r_' + str(int(w.r)-1), \
                                           code_dataset + '_new_prot_list.csv')
    output:
        new_db_profil=os.path.join(wdir, 'r_{r}', 'Pool_clusters.hmm'),
        pool_db=os.path.join(wdir, 'r_{r}', 'db', 'Pool_new_unclustered'),
        new_unclusterd=os.path.join(wdir, 'r_{r}', 'db', 'Pool_unclustered.faa'),
    log: os.path.join(log_dir, "step_1_1_r_{r}.log")
    params:
        rdir='r_{r}',
        diamond="diamond" if diamond else "",
    threads: config.get('threads_blast', 10)
    run:
        if os.path.getsize(input.prots_to_cluster) > 0:
            shell("""
                {script_new_cluster} {params.rdir} {input.prots} \
                 {input.prev_unclusterd} {input.prots_to_cluster} {threads} \
                {diamond} > {log} 2>&1
                  """)
        else:
            shell("touch {output}")

rule step_1_2:
    input:
        prots=rules.step_0_5.output.fasta_file_prots,
        profil=rules.step_1_1.output.new_db_profil,
        prots_to_cluster=lambda w: os.path.join(wdir, \
                                           'r_' + str(int(w.r)-1), \
                                           code_dataset + '_new_prot_list.csv'),
    output:
        hmmsearch_new=os.path.join(wdir, 'r_{r}', \
                                   'Contigs_prots_vs_New_clusters.tab'),
        hmmsearch_bis_new=os.path.join(wdir, 'r_{r}', \
                                       'Contigs_prots_vs_New_clusters.out'),
    log: os.path.join(log_dir, 'step_1_2_r_{r}.log')
    threads: config.get('threads_hmmer', 2)
    run:
        if os.path.getsize(input.prots_to_cluster) > 0:
            with open(input.profil[0]) as DB:
                for line in DB:
                    if line.startswith('NAME'):
                        new_clusters = True
                        break
                else:
                    new_clusters = False

            if new_clusters:
                shell("hmmsearh --tblout {output.hmmsearch_new} --cpu {threads} \
                   -o {output.hmmsearch_bis_new} --noali {input.profil} \
                   {input.prots} > {log} 2>&1")
            else:
                shell("touch {output}")
        else:
            shell("touch {output}")

rule pool_tab_outputs:
    input:
        prev_tab_pooled=lambda w: os.path.join(wdir, \
                                     'r_' + str(int(w.r)-1), \
                                     'Contigs_prots_vs_Phage_{}.tab'.format(w.table_type)),
        tab_new=os.path.join(wdir, 'r_{r}', \
                                   'Contigs_prots_vs_New_{table_type}.tab'),
    output:
        tab_pooled=os.path.join(wdir, 'r_{r}', \
                         'Contigs_prots_vs_Phage_{table_type}.tab')
    shell: """
        if [ -s {input.tab_new} ]; then
            cat {input} > {output}
        else
            ln -s $(realpath  --relative-to=$(dirname {output}) {input.prev_tab_pooled}) {output}
        fi
        """

rule step_1_3:
    input:
        prots=rules.step_0_5.output.fasta_file_prots,
        blastable_unclustered=os.path.join(wdir, 'r_{r}', 'db', \
                                           'Pool_new_unclustered'),
        prots_to_cluster=lambda w: os.path.join(wdir, \
                                           'r_' + str(int(w.r)-1), \
                                           code_dataset + '_new_prot_list.csv'),
    output:
        blast_new_unclustered=os.path.join(wdir, 'r_{r}', \
                             'Contigs_prots_vs_New_unclustered.tab')
    log: os.path.join(log_dir, 'step_1_3_r_{r}.log')
    threads: config.get('threads_blast', 10)
    run:
        if os.path.getsize(input.prots_to_cluster) > 0:
            if diamond:
                shell("diamond blastp --query {input.prots} \
                       --db {input.blastable_unclustered} \
                       --out {output.blast_new_unclustered} --threads {threads} \
                       --outfmt 6 -b 2 --more-sensitive -k 500 --evalue 0.001 \
                       > {log} 2>&1 ")
            else:
                shell("blastp -query {input.prots} -db {input.blastable_unclustered} \
                       -out {output.blast_new_unclustered} -num_threads {threads} \
                       -outfmt 6 -evalue 0.001 > {log} 2>&1 ")
        else:
            shell("touch {output}")


script_merge_annot = \
    os.path.join(script_dir, "Step_2_merge_contigs_annotation.pl")

## Complete the aff
rule step_2:
    input:
        predict_file=rules.step_0_5.output.predict_file,
        pooled_hmm=os.path.join(wdir, 'r_{r}', \
                         'Contigs_prots_vs_Phage_clusters.tab'),
        pooled_blast=os.path.join(wdir, 'r_{r}', \
                         'Contigs_prots_vs_Phage_unclustered.tab'),
        pfama=os.path.join(wdir, 'Contigs_prots_vs_PFAMa.tab'),
        pfamb=os.path.join(wdir, 'Contigs_prots_vs_PFAMb.tab'),
    output:
        affi=os.path.join(wdir, 'r_{r}', code_dataset + "_affi-contigs.csv")
    log: os.path.join(log_dir, 'step_2_r_{r}.log')
    params:
        prots_to_cluster=lambda w: os.path.join(wdir, \
                                           'r_' + str(int(w.r)-1), \
                                           code_dataset + '_new_prot_list.csv'),
        prev_affi=lambda w: os.path.join(wdir, 'r_' + str(int(w.r)-1), \
                                     code_dataset + "_affi-contigs.csv")
    shell:
        """
        if [ -s {params.prots_to_cluster} -o ! -e {params.prots_to_cluster} ]
        then
           {script_merge_annot} {input.predict_file} {input.pooled_hmm} \
            {input.pooled_blast} {input.pfama} {input.pfamb} \
            {ref_phage_clusters} {output.affi} > {log} 2>&1
        else
            ln -s $(realpath  --relative-to=$(dirname {output.affi}) {params.prev_affi}) {output.affi}
        fi
        """

## This generate a csv table including the map of each contig, with PFAM
#and Viral PCs annotations, as well as strand and length of genes


script_detect = os.path.join(script_dir, "Step_3_highlight_phage_signal.pl")
## Complete the summary
rule step_3:
    input:
        affi=rules.step_2.output.affi,
    output:
        phage_fragments=os.path.join(wdir, 'r_{r}', \
                                     code_dataset + '_phage-signal.csv')
        # other output: VST_affi-contigs.refs
    params:
        prots_to_cluster=lambda w: os.path.join(wdir, \
                                           'r_' + str(int(w.r)-1), \
                                           code_dataset + '_new_prot_list.csv'),
        prev_frags=lambda w: os.path.join(wdir, 'r_' + str(int(w.r)-1), \
                                     code_dataset + '_phage-signal.csv')
    log: os.path.join(log_dir, 'step_3_r_{r}.out')
    threads: config.get('threads_blast', 10)
    shell:
        """
        if [ -s {params.prots_to_cluster} -o ! -e {params.prots_to_cluster} ]
        then
            {script_detect} {input.affi} {output.phage_fragments} {threads} \
            > {log} 2>&1
        else
            ln -s $(realpath  --relative-to=$(dirname {output.phage_fragments}) {params.prev_frags}) {output.phage_fragments}
        fi
        """

# Decide which contigs are entirely viral and which are prophages, and
# which of both of these categories are phage enough to be added to the
# databases

script_summary = os.path.join(script_dir, "Step_4_summarize_phage_signal.pl")
rule step_4:
    input:
        affi=rules.step_2.output.affi,
        phage_fragments=rules.step_3.output.phage_fragments,
    output:
        rev_global=os.path.join(wdir, 'r_{r}', \
                                code_dataset + '_global-phage-signal.csv'),
        new_prots=os.path.join(wdir, 'r_{r}', \
                                code_dataset + "_new_prot_list.csv"),
    log: os.path.join(log_dir, 'step_4_r_{r}.out')
    params:
        prots_to_cluster=lambda w: os.path.join(wdir, \
                                           'r_' + str(int(w.r)-1), \
                                           code_dataset + '_new_prot_list.csv'),
        prev_global=lambda w: os.path.join(wdir, \
                                  'r_' + str(int(w.r)-1), \
                                  code_dataset + '_global-phage-signal.csv'),
    shell:
        """
        if [ -e {params.prev_global} ]; then
            cp {params.prev_global} {output.rev_global}
        fi
        if [ -s {params.prots_to_cluster} -o ! -e {params.prots_to_cluster} ]
        then
            {script_summary} {input.affi} {input.phage_fragments} \
               {output.rev_global} {output.new_prots} > {log} 2>&1
        fi
        touch {output.new_prots}
        """

# cop files from last iter into wdir
rule get_final_file:
    input: os.path.join(wdir, 'r_{}'.format(max_revisions), \
                        code_dataset + "_{suffix}.csv")
    output: os.path.join(wdir, code_dataset + "_{suffix}.csv")
    shell: "cp {input} {output}"

# Last step -> extract all sequences as fasta files and gb
script_generate_output = os.path.join(script_dir, 'Step_5_get_phage_fasta-gb.pl')
rule step_5:
    input:
        summary=os.path.join(wdir, code_dataset + "_global-phage-signal.csv"),
        last_affi=os.path.join(wdir, \
                               'r_{}'.format(max_revisions), \
                               code_dataset + "_phage-signal.csv"),
        affi_contigs=os.path.join(wdir, code_dataset + "_affi-contigs.csv"),
        fasta_contigs=os.path.join(fastadir, \
                                    code_dataset + "_nett_filtered.fasta"),
        fasta_prot=os.path.join(fastadir, code_dataset + "_prots.fasta"),
    output:
        phage=expand(os.path.join(wdir,
                                  'Predicted_viral_sequences',
                                  "{cd}_cat-{N}.fasta"),
                     N=[1,2,3], cd=code_dataset),
        prophage=expand(os.path.join(wdir,
                                     'Predicted_viral_sequences',
                                     "{cd}_prophage_cat-{N}.fasta"),
                        N=[4,5,6], cd=code_dataset)
    log: os.path.join(log_dir, 'step_5.out')
    shell: "{script_generate_output} {code_dataset} {wdir} > {log} 2>&1 "

rule collect_viral_fasta:
    input: rules.step_5.output.phage + rules.step_5.output.prophage
    output: os.path.join(wdir, 'Predicted_viral_sequences.fasta')
    shell: "cat {input} {output}"


#TODO: Is there anything below that needs to be done?
"""
#`mv $fastadir $wdir/Fasta_files`;

# We put all results from Hmmsearch and BLAST files in a separate directory
my $store_database_comparison = catdir($wdir, "tab_files");
mkpath($store_database_comparison) unless -d $store_database_comparison;
safe_mv($out_hmmsearch, $store_database_comparison);

# `mv $out_hmmsearch_bis $store_database_comparison/`;
safe_mv($out_blast_unclustered, $store_database_comparison);
safe_mv($out_hmmsearch_pfama, $store_database_comparison);
safe_mv($out_hmmsearch_pfama_bis, $store_database_comparison);
safe_mv($out_hmmsearch_pfamb, $store_database_comparison);
safe_mv($out_hmmsearch_pfamb_bis, $store_database_comparison);

#`mv error.log $log_dir`;
#`mv formatdb.log $log_dir`;
#
#my $final_error_log = catfile($log_dir, 'Virsorter_stderr_log');
#`mv log_err $final_error_log`;
#my $final_out_log = catfile($log_dir, 'Virsorter_stdout_log');
#`mv log_out $final_out_log`;

# We put all the files linked to the metric computation in a new directory
my $store_metric_files = catdir($wdir, 'Metric_files');

if (!-d $store_metric_files) {
    mkpath($store_metric_files);
}

safe_mv($out_file_affi, "$store_metric_files/VIRSorter_affi-contigs.tab");
my $out_file_affi_ref = catdir($wdir, $code_dataset . "_affi-contigs.refs");
safe_mv($out_file_affi_ref, $store_metric_files);
safe_mv($out_file_phage_fragments, "$store_metric_files/VIRSorter_phage_signal.tab");

safe_mv($new_prots_to_cluster, $store_metric_files);
"""

#TODO: Set up the readme file
# And we customize and add the readme file in the output directory
rule readme:
    input: rules.collect_viral_fasta.output
    output: os.path.join(wdir, 'Readme.txt')
    run:
        end_time = datetime.now()  #.isoformat(' ', 'seconds')
        with open(output[0], 'wt') as OUT:
            OUT.write( "VirSorter parameters used :\n")
            OUT.write( "--> Fasta file mined for viral sequences : $input_file\n")
            OUT.write( "--> Viral database used : ")

            if choice_database == 2:
                OUT.write(
                    "Viromes : all bacterial and archaeal virus genomes in Refseq,"
                    " as of January 2014, plus non-redundant predicted genes from viral"
                    " metagenomes (including seawater, freshwater, and human-related"
                    " samples)"
                )
            elif choice_database == 3:
                OUT.write( "Eukaryotic")
            else:
                OUT.write( "RefseqABVir (all bacterial and archaeal virus genomes"
                           " in Refseq, as of January 2014)")

            OUT.write("\n")

            if custom_phage is None:
                OUT.write( "--> No custom reference sequence was added to the database")
            else:
                OUT.write( "--> Custom reference sequences from fasta file $custom_phage"
                           " were added to the database")

            OUT.write("\n")

            if tag_virome:
                OUT.write(
                    "VirSorter was run with the in the 'Virome Decontamination' mode:"
                    " overall metrics for microbial sequences were not evaluated from the"
                    " complete dataset, but instead pre-computed values based on bacterial"
                    " and archaeal genomes from Refseq were used."
                )

            OUT.write("\n")
            OUT.write( "This VirSorter computation began on {}\n".format(start_time.isoformat(' ', 'seconds')))
            OUT.write( "This VirSorter computation finished on {}\n".format(end_time.isoformat(' ', 'seconds')))

            shell("cat {readme_file} >> {output}")
