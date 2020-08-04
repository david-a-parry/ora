CREATE TABLE `homology_member` (
    `homology_id` integer  NOT NULL
    ,  `gene_member_id` integer  NOT NULL
    ,  `seq_member_id` integer  DEFAULT NULL
    ,  `cigar_line` mediumtext DEFAULT NULL
    ,  `perc_cov` float  DEFAULT 0
    ,  `perc_id` float  DEFAULT 0
    ,  `perc_pos` float  DEFAULT 0
    ,  `description` text  NOT NULL
    ,  `is_high_confidence` integer DEFAULT NULL
    ,  PRIMARY KEY (`homology_id`,`gene_member_id`)
);
CREATE TABLE `sequence` (
      `sequence_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
    ,  `length` integer NOT NULL
    ,  `sequence` longtext NOT NULL
);
CREATE TABLE `gene_member` (
      `gene_member_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
    ,  `stable_id` varchar(128) NOT NULL
    ,  `version` integer DEFAULT 0
    ,  `taxon_id` integer  NOT NULL
    ,  `biotype_group` text  NOT NULL DEFAULT 'coding'
    ,  `canonical_member_id` integer  DEFAULT NULL
    ,  `description` text DEFAULT NULL
    ,  `display_label` varchar(128) DEFAULT NULL
    ,  UNIQUE (`stable_id`)
);

CREATE INDEX "idx_homology_member_gene_member_id" ON "homology_member" (`gene_member_id`);
CREATE INDEX "idx_homology_member_seq_member_id" ON "homology_member" (`seq_member_id`);
CREATE INDEX "idx_gene_member_canonical_member_id" ON "gene_member" (`canonical_member_id`);
