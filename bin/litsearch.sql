
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

/* These tables describe characterized sequences with functional features
   such as active site residues, substrate binding residues, or
   other functional residues from mutagenesis experiments.    

   These sequences are stored in a separate BLAST database,
   with ids like db:protId or db:protId:chain
   (and without removal of any redundant sequences)
*/

CREATE TABLE Site(
  db TEXT, /* either SwissProt or PDB */
  id TEXT,
  chain TEXT,
  ligandId TEXT,
  ligandChain TEXT,
  type TEXT,
  posFrom INT,
  posTo INT,
  pdbFrom TEXT,
  pdbTo TEXT,
  comment TEXT,
  pmIds TEXT
);

CREATE INDEX 'IdToSite' ON Site (db,id,chain,ligandId);

CREATE TABLE HasSites(
  db TEXT,
  id TEXT,
  chain TEXT,
  id2 TEXT,
  length INT,
  desc TEXT,
  PRIMARY KEY(db,id,chain)
);

CREATE TABLE SeqHasSite (
  seqHash TEXT, /* the md5 hash of the sequence */
  seqLength INT,
  db TEXT,
  id TEXT,
  chain TEXT,
  PRIMARY KEY (seqHash, seqLength, db, id, chain)
);

CREATE TABLE PdbLigand(
  ligandId TEXT,
  ligandName TEXT,
  PRIMARY KEY (ligandId)
);

/* This stores the descriptions for the clustered PDB entires in pdbClust.faa, used by pdbBlast.cgi */
CREATE TABLE PdbClustInfo(
  id TEXT,
  chain TEXT,
  desc TEXT,
  PRIMARY KEY (id,chain)
);
