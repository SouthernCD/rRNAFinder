from collections import OrderedDict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio import GFF
from toolbiox.lib.common.genome.seq_base import reverse_complement, read_fasta_by_faidx, get_seq_index_ignore_gap
from toolbiox.lib.common.os import rmdir, cmd_run
from toolbiox.lib.common.math.interval import section
from toolbiox.api.common.genome.blast import outfmt5_read
from hugep2g.src.genblasta import GenBlastAJob
from toolbiox.lib.common.fileIO import tsv_file_dict_parse
import os

script_dir_path = os.path.split(os.path.realpath(__file__))[0]

ath_query_rRNA_seq = script_dir_path+"/Ath_rRNA/rRNA.fa"
ath_query_rRNA_bed = script_dir_path+"/Ath_rRNA/rRNA.bed"

ath_query_rRNA_seq = os.path.abspath(ath_query_rRNA_seq)
ath_query_rRNA_bed = os.path.abspath(ath_query_rRNA_bed)


def rRNAFinder_main(args):
    args.fasta_file = os.path.abspath(args.fasta_file)

    if args.query_rRNA_seq:
        args.query_rRNA_seq = os.path.abspath(args.query_rRNA_seq)
    else:
        args.query_rRNA_seq = ath_query_rRNA_seq

    if args.query_rRNA_bed:
        args.query_rRNA_bed = os.path.abspath(args.query_rRNA_bed)
    else:
        args.query_rRNA_bed = ath_query_rRNA_bed

    args.output_dir = os.path.abspath(args.output_dir)
    if not args.log_file is None:
        args.log_file = os.path.abspath(args.log_file)

    # get rRNA unit

    # genblasta
    query_seq_dir = read_fasta_by_faidx(args.query_rRNA_seq)
    target_seq_dir = read_fasta_by_faidx(args.fasta_file)

    query_seq_list = [query_seq_dir[i] for i in query_seq_dir]
    target_seq_list = (target_seq_dir[i] for i in target_seq_dir)

    genblasta_job = GenBlastAJob(
        query_seq_list, target_seq_list, args.output_dir, blast_type='n2n')
    genblasta_job.run()

    genblasta_out_dir = {i[0]: i[1] for i in genblasta_job.parse_output(True)}

    best_genblast_hit = genblasta_out_dir[query_seq_list[0].seqname][0]

    # extract best subject hit
    best_hit_range = best_genblast_hit.rangeB

    s_file = genblasta_job.work_dir + "/best_gba_hit_flank.fa"

    s_flank_start = max(best_hit_range.start - 1000, 1)
    s_flank_end = min(best_hit_range.end + 1000,
                      target_seq_dir[best_hit_range.chr_id].seqs_length())
    s_seq = str(
        target_seq_dir[best_hit_range.chr_id].faidx[s_flank_start - 1:s_flank_end])
    s_base = s_flank_start
    if best_hit_range.strand == "-":
        s_seq = reverse_complement(s_seq)
        s_base = s_flank_end

    with open(s_file, 'w') as s:
        s.write(">%s\n%s" % ("best_gba_hit_flank", s_seq))

    # blastn
    blast_out_file = genblasta_job.work_dir + "/best_gba_hit_flank.bls"

    cmd_string = "makeblastdb -in %s -dbtype nucl" % s_file
    cmd_run(cmd_string, cwd=genblasta_job.work_dir, silence=True)
    cmd_string = "blastn -query %s -db %s -out %s -outfmt 5 -task blastn -evalue 1e-10 " % (
        args.query_rRNA_seq, s_file, blast_out_file)
    cmd_run(cmd_string, cwd=genblasta_job.work_dir, silence=True)

    blastn_output = outfmt5_read(blast_out_file)
    hsp_list = blastn_output[query_seq_list[0].seqname].hit[0].hsp

    query_str_dict = tsv_file_dict_parse(args.query_rRNA_bed, delimiter=",", key_col='name',
                                         fieldnames=['name', 'contig', 'start', 'end'])

    # get each rRNA gene
    rRNA_gene_file = args.output_dir + "/rRNA.gene.fa"

    with open(rRNA_gene_file, 'w') as f:
        type_seq_dict = {}
        for type_name in ['18S', '5.8S', '25S']:
            q_start = int(query_str_dict[type_name]['start'])
            q_end = int(query_str_dict[type_name]['end'])

            # get type overlap hsp
            over_hsp = [j for j in
                        [(section((q_start, q_end), (i.Hsp_query_from, i.Hsp_query_to)), i) for i in hsp_list] if
                        j[0][0]]

            if len(over_hsp) == 0:
                raise ValueError("can't find %s" % type_name)

            # get longest overlap hsp
            section_out, type_hsp = sorted(
                over_hsp, key=lambda x: x[0][1][1] - x[0][1][0], reverse=True)[0]

            deta = section_out[1]

            alig_start = get_seq_index_ignore_gap(
                deta[0], type_hsp.Hsp_qseq, type_hsp.Hsp_query_from, '-')
            alig_end = get_seq_index_ignore_gap(
                deta[1], type_hsp.Hsp_qseq, type_hsp.Hsp_query_from, '-')

            sub_type_seq = type_hsp.Hsp_hseq[alig_start -
                                             1:alig_end].replace('-', '')

            type_seq_dict[type_name] = sub_type_seq

            f.write(">%s\n%s\n" % (type_name, sub_type_seq))

    # blastn again
    gene_vs_flank_bls = genblasta_job.work_dir + "/gene_vs_flank.bls"

    cmd_string = "blastn -query %s -db %s -out %s -outfmt 5 -task blastn -evalue 1e-10 " % (
        rRNA_gene_file, s_file, gene_vs_flank_bls)
    cmd_run(cmd_string, cwd=genblasta_job.work_dir, silence=True)

    blastn_output = outfmt5_read(gene_vs_flank_bls)

    type_range = {}
    for type_name in ['18S', '5.8S', '25S']:
        top_hit = blastn_output[type_name].hit[0].hsp[0]
        type_range[type_name] = (top_hit.Hsp_hit_from, top_hit.Hsp_hit_to)

    type_range['ITS1'] = (type_range['18S'][1] + 1, type_range['5.8S'][0] - 1)
    type_range['ITS2'] = (type_range['5.8S'][1] + 1, type_range['25S'][0] - 1)

    type_range_dict = OrderedDict()
    for i in ['18S', 'ITS1', '5.8S', 'ITS2', '25S']:
        type_range_dict[i] = type_range[i]

    # output
    flank_seq = read_fasta_by_faidx(s_file)['best_gba_hit_flank']

    unit_fasta_file = args.output_dir + "/rRNA.unit.fa"
    gene_type_fasta_file = args.output_dir + "/rRNA.gene.fa"
    unit_gff_file = args.output_dir + "/rRNA.unit.gff"

    rna_unit = Seq(flank_seq.sub(
        type_range_dict['18S'][0], type_range_dict['25S'][1], "+", False))
    rec = SeqRecord(rna_unit, "rRNA_unit", description='')

    with open(unit_fasta_file, 'w') as f:
        f.write(rec.format('fasta'))

    type_noflank_range_dict = OrderedDict()
    flank_start = type_range_dict['18S'][0]
    with open(gene_type_fasta_file, 'w') as f:
        for i in ['18S', 'ITS1', '5.8S', 'ITS2', '25S']:
            type_noflank_range_dict[i] = (
                type_range_dict[i][0] - flank_start + 1, type_range_dict[i][1] - flank_start + 1)

            type_seq = Seq(flank_seq.sub(
                type_range_dict[i][0], type_range_dict[i][1], "+", False))
            type_rec = SeqRecord(type_seq, i, description='')
            f.write(type_rec.format('fasta'))

    rec.features = []
    qualifiers = {
        "source": "rRNAFinder",
        "ID": 'rRNA'
    }
    top_feature = SeqFeature(FeatureLocation(type_noflank_range_dict['18S'][0] - 1,
                                             type_noflank_range_dict['25S'][1], strand=1),
                             type="rRNA", qualifiers=qualifiers)
    top_feature.sub_features = []
    for i in ['18S', 'ITS1', '5.8S', 'ITS2', '25S']:

        sub_qualifiers = {
            "source": "rRNAFinder",
            "ID": i
        }
        if i == 'ITS1' or i == 'ITS2':
            type_tmp = 'intron'
        else:
            type_tmp = 'rRNA'

        sub_feature = SeqFeature(FeatureLocation(type_noflank_range_dict[i][0] - 1,
                                                 type_noflank_range_dict[i][1], strand=1),
                                 type=type_tmp, qualifiers=sub_qualifiers)

        top_feature.sub_features.append(sub_feature)

    rec.features.append(top_feature)

    with open(unit_gff_file, "w") as f:
        GFF.write([rec], f)

    rmdir(genblasta_job.work_dir)
