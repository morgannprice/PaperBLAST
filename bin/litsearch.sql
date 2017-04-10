
/* No PRIMARY KEY because may have duplicate entries */
CREATE TABLE GenePaper(
	geneId TEXT,
        queryTerm TEXT,
	pmcId TEXT, /* including the leading PMC */
	pmId TEXT,
	doi TEXT,
        title TEXT,
	authors TEXT,
	journal TEXT,
	year TEXT,
        isOpen INT /* missing for GeneRIF links */
);
CREATE INDEX 'GeneToPaper' ON GenePaper ('geneId' ASC);

CREATE TABLE Snippet(
	geneId TEXT,
        queryTerm TEXT,
        pmcId TEXT,
        pmId TEXT,
        snippet TEXT
);
CREATE INDEX 'GeneToSnippet' ON Snippet ('geneId' ASC, 'queryTerm' ASC, 'pmcId' ASC, 'pmId' ASC);

CREATE TABLE PaperAccess(
	pmcId TEXT,
        pmId TEXT,
        access TEXT,
        PRIMARY KEY (pmcId,pmId)
);

/* Note that some genes may be in Gene but have no links */
CREATE TABLE Gene(
	geneId TEXT NOT NULL,
	organism TEXT NOT NULL,
        protein_length INT NOT NULL,
        desc TEXT, /* may be missing */
        PRIMARY KEY (geneId)
);

CREATE TABLE UniProt(
	acc TEXT NOT NULL,
        desc TEXT NOT NULL,
        organism TEXT NOT NULL,
        comment TEXT NOT NULL,
        protein_length INT NOT NULL,
        PRIMARY KEY (acc)
);

CREATE TABLE GeneRIF(
	/* geneId is actually a fully specified protein id like YP_006960813.1 */
	geneId TEXT NOT NULL,
        pmcId TEXT NOT NULL,
        pmId TEXT NOT NULL,
        comment TEXT NOT NULL
);
CREATE INDEX 'GeneToRIF' ON GeneRIF ('geneId' ASC, 'pmcId' ASC, 'pmId' ASC);

CREATE TABLE EcoCyc(
       /* The protein_id will be something like 1-PFK-MONOMER
          In contrast, the protein_id in SeqToDuplicate or in the fasta file will be
          gnl|ECOLI|1-PFK-MONOMER */
	protein_id TEXT NOT NULL,
        protein_name TEXT NOT NULL,
        bnumber TEXT NOT NULL,
        desc TEXT NOT NULL,
        protein_length INT NOT NULL,
        PRIMARY KEY (protein_id)
);

CREATE TABLE EcoCycToPubMed(
	protein_id TEXT NOT NULL,
        pmId TEXT NOT NULL,
        PRIMARY KEY (protein_id, pmId)
);

/* From an (arbitrarily chosen) id with a unique sequence, to the identical sequence(s).
   Singleton sequences (with no duplicates) will not be listed. */
CREATE TABLE SeqToDuplicate(
	sequence_id TEXT NOT NULL,
        duplicate_id TEXT NOT NULL,
        PRIMARY KEY (sequence_id, duplicate_id)
);
CREATE INDEX 'DupToId' ON SeqToDuplicate ('duplicate_id' ASC, 'sequence_id' ASC);
