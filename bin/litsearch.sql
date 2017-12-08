
/* This table includes text-mined links and GeneRIF links, but not curated links.
 * It has no PRIMARY KEY because it can have duplicate entries.
 */
CREATE TABLE GenePaper(
	geneId TEXT,
        queryTerm TEXT,
	pmcId TEXT, /* PubMedCentral identifier, including the leading PMC */
	pmId TEXT, /* PubMed identifier */
	doi TEXT, /* digital object identifier */
        title TEXT,
	authors TEXT,
	journal TEXT,
	year TEXT,
        isOpen INT /* missing for GeneRIF links */
);
CREATE INDEX 'GeneToPaper' ON GenePaper ('geneId' ASC);

/* Note that some genes may be in Gene but have no links */
CREATE TABLE Gene(
	geneId TEXT NOT NULL,
	organism TEXT NOT NULL,
        protein_length INT NOT NULL,
        desc TEXT, /* may be missing */
        PRIMARY KEY (geneId)
);

/* This table stores snippets extracted from the text of articles */
CREATE TABLE Snippet(
	geneId TEXT,
        queryTerm TEXT,
        pmcId TEXT,
        pmId TEXT,
        snippet TEXT
);
CREATE INDEX 'GeneToSnippet' ON Snippet ('geneId' ASC, 'queryTerm' ASC, 'pmcId' ASC, 'pmId' ASC);

/* This table stores comments from GeneRIF about the gene-paper links */
CREATE TABLE GeneRIF(
	/* geneId is actually a fully specified protein id like YP_006960813.1 */
	geneId TEXT NOT NULL,
        pmcId TEXT NOT NULL,
        pmId TEXT NOT NULL,
        comment TEXT NOT NULL
);
CREATE INDEX 'GeneToRIF' ON GeneRIF ('geneId' ASC, 'pmcId' ASC, 'pmId' ASC);

CREATE TABLE PaperAccess(
	pmcId TEXT,
        pmId TEXT,
        access TEXT,
        PRIMARY KEY (pmcId,pmId)
);

/* The curated proteins */
CREATE TABLE CuratedGene(
	/* The sequence id for these entries should be db::protId */
	db TEXT NOT NULL,
        protId TEXT NOT NULL,
        id2 TEXT NOT NULL, /* may be empty */
        name TEXT NOT NULL, /* may be empty */
        desc TEXT NOT NULL,
        organism TEXT NOT NULL, /* may be empty */
        protein_length INT NOT NULL,
        comment TEXT NOT NULL, /* may be empty */
        PRIMARY KEY (db,protId)
);

/* Link curated proteins to papers */
CREATE TABLE CuratedPaper(
	db TEXT NOT NULL,
        protId TEXT NOT NULL,
        pmId TEXT NOT NULL,
        PRIMARY KEY (db,protId,pmId)
);

/* This table supports a size reduction of the protein database
   by identifying identical sequences.

   It maps an (arbitrarily chosen) id with a unique sequence, which appears in
   the BLAST database, to the identical sequence(s), which do not.
   Singleton sequences (with no duplicates) are not listed.
*/
CREATE TABLE SeqToDuplicate(
	sequence_id TEXT NOT NULL,
        duplicate_id TEXT NOT NULL,
        PRIMARY KEY (sequence_id, duplicate_id)
);
CREATE INDEX 'DupToId' ON SeqToDuplicate ('duplicate_id' ASC, 'sequence_id' ASC);
