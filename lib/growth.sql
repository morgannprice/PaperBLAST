/* GapMind can display growth information related to whether a pathway
   is present or not. 
*/

CREATE TABLE Growth(
  orgId TEXT NOT NULL,
  pathwayId TEXT NOT NULL,
  /* -1 for no growth, 0 for ambiguous, 1 for growth.
     Usually the row would be omitted if there is no data
  */
  status INT NOT NULL,
  comment TEXT NOT NULL, /* could be empty, or might give quantification */
  URL TEXT NOT NULL,     /* could be empty */
  PRIMARY KEY (orgId,pathwayId)
);
CREATE INDEX 'GrowthByPathawy' on Growth (pathwayId,orgId);
