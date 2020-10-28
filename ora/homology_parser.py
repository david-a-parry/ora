from ora.alignments import get_align_pos, align_pos_to_amino_acid
from collections import namedtuple, defaultdict
from ora.uniprot_lookups import get_uniprot_features, feat_fields

ParalogLookup = namedtuple("ParalogLookup", "gene protein position")
uniprot_lookups = set()
header_fields = ['Query_Symbol', 'Query_Gene', 'Query_Protein', 'Query_Pos',
                 'Query_AA', 'Homolog_Symbol', 'Homolog_Gene',
                 'Homolog_Protein', 'Orthology_Type', 'Species', 'Percent_ID',
                 'Percent_Pos', 'Homolog_Pos', 'Homolog_AA']


def parse_homology_data(homs, positions, logger, protein_id=None,
                        skip_paralogs=False, output_all_orthologs=False):
    '''
        From a list of homologies (generated from Ensembl REST look-ups
        or ora.local_lookups) return results detailing homologies and any
        Uniprot features overlapping given positions.

        Args:
             homs: list of homology results as generated by
                   EnsemblRestQueries.get_homologies method.

             positions:
                   list of protein positions in query protein to check.

             logger:
                   logger object from caller.

             protein_id:
                   restrict lookups to source proteins matching this
                   protein ID

            skip_paralogs:
                   do not perform lookups for paralogs

            output_all_orthologs:
                   set 'should_output' to True for all ortholog results
                   even if no overlapping Uniprot features are identified.

    '''
    paralogs = defaultdict(list)
    results = []
    logger.info("Got {:,} homologs".format(len(homs)))
    for pos in positions:
        for i in range(len(homs)):
            s_id = homs[i]['source']['id']
            t_id = homs[i]['target']['id']
            s_protein = homs[i]['source']['protein_id']
            t_protein = homs[i]['target']['protein_id']
            s_seq = homs[i]['source']['align_seq']
            t_seq = homs[i]['target']['align_seq']
            p = get_align_pos(s_seq, pos)
            if (s_protein, pos) not in uniprot_lookups:
                ufeats = get_uniprot_features(s_protein, pos, pos)
                if ufeats:
                    for f in ufeats:
                        result = dict(
                            query_gene=s_id,
                            query_symbol=homs[i]['source'].get('symbol'),
                            query_protein=s_protein,
                            query_pos=pos,
                            query_aa=s_seq[p],
                            homolog_gene=s_id,
                            homolog_symbol=homs[i]['source'].get('symbol'),
                            homolog_protein=s_protein,
                            orthology_type="self",
                            species=homs[i]['source']['species'],
                            percent_id=100,
                            percent_pos=100,
                            homolog_pos=pos,
                            homolog_aa=s_seq[p],
                            query_species=homs[i]['source']['species'],
                            features=f,
                            should_output=True)
                        results.append(result)
                uniprot_lookups.add((s_protein, pos))
            if protein_id is not None and s_protein != protein_id:
                logger.info("Skipping paralog lookup for protein {}\n".format(
                    protein_id))
                return results, paralogs
            hom_type = homs[i]['type']
            if skip_paralogs and 'paralog' in hom_type:
                continue
            o, aa = align_pos_to_amino_acid(t_seq, p) if p > 0 else (-1, '-')
            if 'paralog' in hom_type and o > 0:
                paralogs[t_id].append(ParalogLookup(t_id, t_protein, o))
            if o > 0 and (t_protein, o) not in uniprot_lookups:
                result_template = dict(
                                    query_gene=s_id,
                                    query_protein=s_protein,
                                    query_symbol=homs[i]['source'].get(
                                        'symbol'),
                                    query_pos=pos,
                                    query_aa=s_seq[p],
                                    homolog_gene=t_id,
                                    homolog_protein=t_protein,
                                    homolog_symbol=homs[i]['target'].get(
                                        'symbol'),
                                    orthology_type=hom_type,
                                    species=homs[i]['target']['species'],
                                    percent_id=homs[i]['target']['perc_id'],
                                    percent_pos=homs[i]['target']['perc_pos'],
                                    homolog_pos=o,
                                    query_seq=s_seq,
                                    homolog_seq=t_seq,
                                    query_species=homs[i]['source']['species'],
                                    homolog_aa=aa,
                                    should_output=False)
                ufeats = get_uniprot_features(t_protein, o, o)
                if ufeats:
                    for f in ufeats:
                        result = result_template.copy()
                        result['features'] = f
                        result['should_output'] = True
                        results.append(result)
                else:
                    result = result_template.copy()
                    result['features'] = dict((x, '-') for x in feat_fields)
                    result['should_output'] = output_all_orthologs
                    results.append(result)
                uniprot_lookups.add((t_protein, o))
    if not skip_paralogs:
        logger.info("Got {} paralog sequences".format(len(paralogs)))
    return results, paralogs


def check_paralog_lookups(original_ensg, output_results, results):
    '''
        If we have results with a paralog of our original query gene as
        the query gene we should output an alignment of our original
        query gene and the paralog if we aren't already doing so. This
        returns results for our original query gene and these paralogs
        if not already in our output_results.
    '''
    # All pairwise comparisons in our output list
    pair_comps = set((r['query_gene'], r['homolog_gene']) for r in
                     output_results)
    # Pairwise comparisons where the query gene is paralogous to our query
    para_comps = set((r['query_gene'], r['homolog_gene']) for r in
                     output_results if r['query_gene'] != original_ensg)
    extra = list()
    for pair in para_comps:
        if (original_ensg, pair[0]) not in pair_comps:
            extra.extend(r for r in results if r['query_gene'] == original_ensg
                         and r['homolog_gene'] == pair[0])
    return extra
