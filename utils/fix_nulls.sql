UPDATE gene_member SET canonical_member_id = NULL WHERE canonical_member_id = '\N';
UPDATE gene_member SET display_label = NULL WHERE display_label = '\N';
UPDATE homology_member SET is_high_confidence = NULL WHERE is_high_confidence = '\N';
