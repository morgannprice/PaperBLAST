
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
        snippet TEXT
);
CREATE INDEX 'GeneToSnippet' ON Snippet ('geneId' ASC, 'queryTerm' ASC, 'pmcId' ASC);

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
