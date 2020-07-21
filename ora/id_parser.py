import re

_ccds_re = re.compile(r'CCDS\d+(\.\d+)*')
_ensg_re = re.compile(r'ENS\w*G\d{11}(\.\d+)*')
_enst_re = re.compile(r'ENS\w*T\d{11}(\.\d+)*')
_ensp_re = re.compile(r'ENS\w*P\d{11}(\.\d+)*')
_entrez_re = re.compile(r'^\d+$')
_refrna_re = re.compile(r'[XN][MR]_\d+(\.\d+)*')
_refpro_re = re.compile(r'[XN]P_\d+(\.\d+)*')
_uniprot_re = re.compile(r'''^[OPQ][0-9][A-Z0-9]{3}[0-9]$|
                            ^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$''',
                         re.X)


def parse_id(id):
    '''
        For given identifier, returns a dict with the following keys:

            identifier_type: string representing best guess at identifier type

            is_transcript:   boolean indicating whether identifier type
                             represents transcript sequence

            is_protein:      boolean indicating whether identifier type
                             represents protein sequence

            is_ensembl_id:   boolean indicating whether identifier comes from
                             Ensembl

    '''
    result = dict(identifier_type=None,
                  is_transcript=False,
                  is_protein=False,
                  is_ensembl_id=False)
    if _ensg_re.match(id):
        result['is_ensembl_id'] = True
        result['identifier_type'] = "Ensembl Gene ID"
    elif _enst_re.match(id):
        result['is_transcript'] = True
        result['is_ensembl_id'] = True
        result['identifier_type'] = "Ensembl Transcript ID"
    elif _ensp_re.match(id):
        result['is_ensembl_id'] = True
        result['is_protein'] = True
        result['identifier_type'] = "Ensembl Protein ID"
    elif _refrna_re.match(id):
        result['is_transcript'] = True
        result['identifier_type'] = "RefSeq mRNA ID"
    elif _refpro_re.match(id):
        result['is_protein'] = True
        result['identifier_type'] = "RefSeq Protein ID"
    elif _uniprot_re.match(id):
        result['is_transcript'] = True
        result['is_protein'] = True
        result['identifier_type'] = "Uniprot ID"
    elif _ccds_re.match(id):
        result['is_transcript'] = True
        result['identifier_type'] = "CCDS ID"
    elif _entrez_re.match(id):
        result['identifier_type'] = "Entrez Gene ID"
    else:
        result['identifier_type'] = "Gene Symbol or other"
    return result
