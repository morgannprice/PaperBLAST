
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
        isOpen INT
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

CREATE TABLE EcoCyc(
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

# From an (arbitrarily chosen) id with a unique sequence, to the identical sequence(s).
# Singleton sequences (with no duplicates) will not be listed.
CREATE TABLE SeqToDuplicate(
	sequence_id TEXT NOT NULL,
        duplicate_id TEXT NOT NULL,
        PRIMARY KEY (sequence_id, duplicate_id)
);
