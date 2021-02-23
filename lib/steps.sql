/* A GapMind steps database includes:
   Descriptions of all rules and steps (parsed from Step files)
   Proteins and families that match those steps (as computed by gapquery.pl)
   Requirements (links between pathways)
   Known gaps
*/

/* Describe the pathways. Each pathway corresponds to a steps file,
   except the special pathwayId all, which describes the set of pathways.
*/
CREATE TABLE Pathway(
  pathwayId TEXT PRIMARY KEY,  /* i.e., arg for arginine synthesis */
  desc TEXT NOT NULL
);

CREATE TABLE Rule(
  pathwayId TEXT NOT NULL,
  /* Each pathway defines the rule all to achieve the entire pathway.
     Another common ruleId is "transport", for uptake of a catabolized compound */
  ruleId TEXT NOT NULL,
  PRIMARY KEY (pathwayId, ruleId)
);

/* A rule can be satisified by any of its instances.  The instances
   are in the order in the steps file; this also puts the rules in the
   correct order so that dependencies come first */
CREATE TABLE RuleInstance(
  pathwayId TEXT NOT NULL,
  ruleId TEXT NOT NULL,
  instanceId INT NOT NULL, /* increasing in order in the steps file */
  PRIMARY KEY (pathwayId, ruleId, instanceId)
);

/* Each instance has 1 or more components, either a subrule or a step */
CREATE TABLE InstanceComponent(
  pathwayId TEXT NOT NULL,
  ruleId TEXT NOT NULL,
  instanceId INT NOT NULL,
  componentId INT NOT NULL, /* increasing in order in the steps file */
  /* subRuleId or stepId will be empty */
  subRuleId TEXT NOT NULL,
  stepId TEXT NOT NULL,
  PRIMARY KEY (pathwayId, ruleId, instanceId, componentId)
);

CREATE TABLE Step(
  pathwayId TEXT NOT NULL,
  stepId TEXT NOT NULL,
  desc TEXT NOT NULL,
  PRIMARY KEY (pathwayId, stepId)
);

/* A step is defined by various terms or identifiers */
CREATE TABLE StepPart(
  pathwayId TEXT NOT NULL,
  stepId TEXT NOT NULL,
  partId INT NOT NULL, /* increasing in order in the steps file */
  /* partType must be EC, curated, uniprot, term, hmm, ignore, ignore_other, or ignore_hmm */
  partType TEXT NOT NULL,
  value TEXT NOT NULL,
  PRIMARY KEY (pathwayId, stepId, partId)
);  

/* A step is implemented as a list of queries. These don't
   correspond to StepPart in a 1:1 or 1:many way */
CREATE TABLE StepQuery(
  pathwayId TEXT NOT NULL,
  stepId TEXT NOT NULL,
  queryId INT NOT NULL, /* increasing in order */
  /* queryType must be curated, curated2, hmm, uniprot, or ignore.
     (curated2 refers to curated annotations of unchracterized proteins;
      ignore means that similarity to a characterized protein that is
      not part of the step does not lower the score) */
  queryType TEXT NOT NULL,
  /* Use empty instead of NULL for missing fields */
  curatedIds TEXT NOT NULL,
  uniprotId TEXT NOT NULL,
  protId TEXT NOT NULL, /* for curated2 */
  hmmId TEXT NOT NULL, /* i.e., PF02965 */
  hmmFileName TEXT NOT NULL, /* i.e., PF02965.17.hmm */
  /* The description is redundant with information in the curated database for
     queryType = curated, curated2, or ignore */
  desc TEXT NOT NULL,
  seq TEXT NOT NULL, /* not used for hmm */
  /* The HMM definitions are not included */
  PRIMARY KEY (pathwayId, stepId, queryId)
);

/* Store the actual models */
CREATE TABLE HMM(
  hmmId TEXT PRIMARY KEY,
  hmm BLOB NOT NULL /* HMMer format */
);

CREATE TABLE Requirement(
  pathwayId TEXT NOT NULL,
  ruleId TEXT NOT NULL, /* or "all" */
  requiredPathwayId TEXT NOT NULL,
  requiredRuleId TEXT NOT NULL, /* or "all", or empty if requiredStepId is set */
  requiredStepId TEXT NOT NULL,
  isNot INT NOT NULL, /* 0 or 1 */
  comment TEXT NOT NULL,
  PRIMARY KEY (pathwayId, ruleId, requiredPathwayId, requiredRuleId, requiredStepId)
);

/* Describe known gaps, that is, cases where a step is missing but the organism
   performs the pathway. These optionally have curated classifications and descriptions */
CREATE TABLE KnownGap(
  gdb TEXT NOT NULL, /* which genome database */
  gid TEXT NOT NULL, /* genome identifier */
  genomeName TEXT NOT NULL, /* for description */
  pathwayId TEXT NOT NULL,
  stepId TEXT NOT NULL,
  /* If not-empty, should be spurious (i.e., gene model is missing), diverged, or novel. */
  gapClass TEXT NOT NULL,
  comment TEXT NOT NULL,
  PRIMARY KEY (gdb, gid, pathwayId, stepId)
);

CREATE TABLE KnownGapMarker(
  gdb TEXT NOT NULL,
  gid TEXT NOT NULL,
  markerId TEXT NOT NULL,
  markerSeq TEXT NOT NULL, /* amino acid sequence of that marker gene */
  PRIMARY KEY (gdb, gid, markerId)
);
