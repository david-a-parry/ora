import csv
import gzip
import logging
import os
from collections import defaultdict

logger = logging.getLogger("Uniprot Lookups")
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
ch = logging.StreamHandler()
ch.setLevel(logger.level)
ch.setFormatter(formatter)
logger.addHandler(ch)

_sprot_prefix = os.path.join(os.path.dirname(__file__),
                             "data",
                             "sprot.")
_ens_id_tab = _sprot_prefix + 'ens_ids.tab.gz'
_features_tab = _sprot_prefix + 'features.tab.gz'
_isoform_seq_tab = _sprot_prefix + 'isoform_seqs.tab.gz'

ensp2uniprot = dict()
uniprot2feats = defaultdict(list)
feature_lookups = dict()


def get_uniprot_features(ensp, start, stop):
    '''
        Get Uniprot features from an Ensembl protein identifier and
        amino acid position.
    '''
    if not ensp2uniprot:
        _features_from_uniprot()
        _uniprot_from_ensp()
    ufeats = feature_lookups.get(ensp, feats_from_ensp(ensp))
    if ufeats is None:
        logger.debug("Could not find Uniprot features for {}".format(ensp))
        return None
    pos_feats = []
    for f in ufeats:
        if int(f['Start']) <= stop and int(f['Stop']) >= start:
            pos_feats.append(f)
    return pos_feats


def feats_from_ensp(ensp):
    u_info = ensp2uniprot.get(ensp)
    if u_info is None:
        logger.debug("Could not find Uniprot ID for {}".format(ensp))
        return None
    if u_info['displayed']:
        return uniprot2feats[u_info['id']]
    logger.warn("Found non-displayed Uniprot isoform {} for {} ".format(
        u_info['isoform'],
        ensp) + " - can not retrieve features for this isoform.")
    return None


def _features_from_uniprot():
    feats = list()
    found = False
    with gzip.open(_features_tab, 'rt') as tabfile:
        reader = csv.DictReader(tabfile, delimiter='\t')
        for f in ['UniprotID', 'Feature', 'Start', 'Stop', 'Description']:
            if f not in reader.fieldnames:
                raise ValueError("Invalid header for {}".format(_features_tab))
            for row in reader:
                uniprot2feats[row['UniprotID']].append(row)


def _uniprot_from_ensp():
    with gzip.open(_ens_id_tab, 'rt') as tabfile:
        reader = csv.DictReader(tabfile, delimiter='\t')
        for f in ['UniprotID', 'Protein', 'Isoform', 'Displayed']:
            if f not in reader.fieldnames:
                raise ValueError("Invalid header for {}".format(_ens_id_tab))
            for row in reader:
                ensp2uniprot[row['Protein']] = dict(
                    id=row['UniprotID'],
                    isoform=row['Isoform'],
                    displayed=bool(int(row['Displayed'])))
