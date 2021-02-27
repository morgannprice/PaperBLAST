/* A GapMind gaps database includes:
   All of the candidates for all of the steps
   The scoring of each step
   The scoring of each rule, including best paths
   Similarities to markers (for known gaps)
*/

CREATE TABLE Candidate(
  orgId TEXT NOT NULL,
  pathwayId TEXT NOT NULL,
  stepId TEXT NOT NULL,
  locusId TEXT NOT NULL,
  sysName TEXT NOT NULL, /* also known as locus tags */
  desc TEXT NOT NULL,
  /* locusId2 and sysName2 are empty except for split candidates */
  locusId2 TEXT NOT NULL,
  sysName2 TEXT NOT NULL,
  desc2 TEXT NOT NULL,
  score INT NOT NULL, /* 0, 1, or 2 */
  blastBits REAL NOT NULL,
  curatedIds TEXT NOT NULL, /* curatedIds, or curated2:protId, or uniprot:uniprotId */
  identity REAL NOT NULL,
  blastCoverage REAL NOT NULL,
  blastScore INT NOT NULL, /* 0, 1, or 2 */
  hmmId TEXT NOT NULL,
  hmmDesc TEXT NOT NULL,
  hmmBits REAL NOT NULL,
  hmmCoverage REAL NOT NULL,
  hmmScore INT NOT NULL, /* 0, 1, or 2 */
  otherIds TEXT NOT NULL, /* a curatedIds */
  otherBits REAL NOT NULL,
  otherIdentity REAL NOT NULL,
  otherCoverage REAL NOT NULL,
  PRIMARY KEY (orgId, pathwayId, stepId, locusId, locusId2)
);
CREATE INDEX 'CandidateByLocus' on Candidate (orgId, locusId);
CREATE INDEX 'CandidateByLocus2' on Candidate (orgId, locusId2);

CREATE TABLE StepScore(
  orgId TEXT NOT NULL,
  pathwayId TEXT NOT NULL,
  stepId TEXT NOT NULL,
  onBestPath INT NOT NULL, /* 0 or 1 */
  /* the best scoring candidate (if any) */
  score INT NOT NULL, /* 0, 1, or 2, or empty if there are no candidates at all */
  locusId TEXT NOT NULL, /* comma delimited if this is a split hit */
  sysName TEXT NOT NULL, /* ditto */
  /* the 2nd-best scoring candidate (if any) */
  score2 INT NOT NULL,
  locusId2 TEXT NOT NULL,
  sysName2 TEXT NOT NULL,
  PRIMARY KEY (orgId, pathwayId, stepId)
);

CREATE TABLE RuleScore(
  orgId TEXT NOT NULL,
  pathwayId TEXT NOT NULL,
  ruleId TEXT NOT NULL,
  /* number of high, medium, or low-confidence steps */
  nHi INT NOT NULL,
  nMed INT NOT NULL,
  nLo INT NOT NULL,
  score REAL NOT NULL,
  expandedPath TEXT NOT NULL, /* space-delimited list of steps used */
  path TEXT NOT NULL, /* space-delimited list of rules/steps used */
  path2 TEXT NOT NULL, /* 2nd-best path */
  PRIMARY KEY (orgId, pathwayId, ruleId)
);

CREATE TABLE RequirementNotMet(
  orgId TEXT NOT NULL,
  pathwayId TEXT NOT NULL,
  ruleId TEXT NOT NULL, /* or "all" */
  requiredPathwayId TEXT NOT NULL,
  requiredRuleId TEXT NOT NULL, /* or "all", or empty if requiredStepId is set */
  requiredStepId TEXT NOT NULL,
  isNot INT NOT NULL, /* 0 or 1 */
  comment TEXT NOT NULL,
  PRIMARY KEY (orgId, pathwayId, ruleId)
);

CREATE TABLE MarkerSimilarity(
  orgId TEXT NOT NULL,
  hitOrgId TEXT NOT NULL,
  identity REAL NOT NULL, /* up to 100% */
  nMarkers INT NOT NULL,
  PRIMARY KEY (orgId, hitOrgId)
);
