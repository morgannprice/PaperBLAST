/* Describe all curated proteins and their sequences.
   CuratedInfo and CuratedSeq include only proteins with experimental evidence
   about their function.
   "Curated2" describes curated annotations with no experimental evidence.
*/
CREATE TABLE CuratedInfo(
  /* To handle curation of the identical sequence in different databases,
     curatedIds is a comma-delimited composite such as "SwissProt::P0C079,ecocyc::EG10836-MONOMER".
     Similarly, descs, id2s, and orgs are delimited with ";; " (in the same order as the ids)
     orgs and ids2 may be missing, in which case they are the empty string
  */
  curatedIds TEXT PRIMARY KEY,
  seqLength INT NOT NULL,
  descs TEXT NOT NULL,
  id2s TEXT NOT NULL, /* empty string if missing */
  orgs TEXT NOT NULL /* empty string if missing */
);

/* Make it possible to go from an individual curated id to the combination curatedIds.
   Trivial cases (with curatedId = curatedIds) are included
*/
CREATE TABLE CuratedIdToIds(
  curatedId TEXT PRIMARY KEY,
  curatedIds TEXT NOT NULL
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
  curatedIds TEXT PRIMARY KEY,
  comment TEXT
);

/* Links of reactions to compounds */
CREATE Table CompoundInReaction(
  rxnId TEXT NOT NULL, /* metacyc:metacycRxnId or RHEA:number */
  rxnLocation TEXT NTO NULL,
  cmpId TEXT NOT NULL, /* a metacyc compound id */
  cmpDesc TEXT NOT NULL,
  side INT NOT NULL, /* -1 for left side, +1 for right side */
  coefficient TEXT NOT NULL, /* usually a number, but sometimes an expression like 'n' or 'n+1' */
  compartment TEXT NOT NULL, /* the empty string, NIL, CYTOSOL, In, OUT, MEMBRANE, MIDDLE, PM-BAC-NEG */
  PRIMARY KEY (rxnId,cmpId,side)
);

CREATE TABLE EnzymeForReaction(
  curatedIds TEXT NOT NULL,
  rxnId TEXT NOT NULL,
  enzDesc TEXT NOT NULL,
  PRIMARY KEY (curatedIds,rxnId)
);
CREATE INDEX 'EnzymeByReaction' on EnzymeForReaction (rxnId,curatedIds);

/* Links of transporters to substrates (TCDB only) */
CREATE TABLE TransporterSubstrate(
  curatedIds TEXT NOT NULL,
  substrate TEXT NOT NULL
);
