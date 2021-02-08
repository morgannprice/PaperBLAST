/* Describe all curated proteins and their sequences.
   CuratedInfo and CuratedSeq include only proteins with experimental evidence
   about their function.
   "Curated2" describes curated annotations with no experimental evidence.
*/
CREATE TABLE CuratedInfo(
  /* To handle curation of the identical sequence in different databases,
     curatedIds is a comma-delimited composite such as "SwissProt::P0C079,ecocyc::EG10836-MONOMER".
     Similarly, descs is delimited with ";;" (in the same order as the ids)
  */
  curatedIds TEXT PRIMARY KEY,
  seqLength INT NOT NULL,
  descs TEXT NOT NULL
);

CREATE TABLE CuratedSeq(
  curatedIds TEXT PRIMARY KEY,
  seq TEXT NOT NULL
);

CREATE TABLE CuratedPFam(
  curatedIds TEXT NOT NULL,
  hmmName TEXT NOT NULL,
  hmmAcc TEXT NOT NULL,
  evalue REAL NOT NULL,
  bits REAL NOT NULL,
  seqFrom INT NOT NULL,
  seqTo INT NOT NULL,
  seqLen INT NOT NULL,
  hmmFrom INT NOT NULL,
  hmmTo INT NOT NULL,
  hmmLen INT NOT NULL
);
CREATE INDEX 'PFamByIds' on CuratedPFam ('curatedIds');
CREATE INDEX 'CuratedByPFam' on CuratedPFam ('hmmAcc');

CREATE TABLE Curated2(
  protId TEXT PRIMARY KEY,
  desc TEXT NOT NULL,
  seq TEXT NOT NULL
);

/* Only heteromeric proteins are included */
CREATE TABLE Hetero(
  curatedId TEXT PRIMARY KEY, /* i.e., "BRENDA::P38540" -- just one id */
  comment TEXT
);



  
  
