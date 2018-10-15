/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asnvalid.h13";
static AsnValxNode avnx[11] = {
    {20,"none" ,0,0.0,&avnx[1] } ,
    {20,"info" ,1,0.0,&avnx[2] } ,
    {20,"warning" ,2,0.0,&avnx[3] } ,
    {20,"error" ,3,0.0,&avnx[4] } ,
    {20,"reject" ,4,0.0,&avnx[5] } ,
    {20,"fatal" ,5,0.0,NULL } ,
    {2,NULL,0,0.0,NULL } ,
    {3,NULL,2,0.0,NULL } ,
    {2,NULL,0,0.0,NULL } ,
    {2,NULL,1,0.0,NULL } ,
    {2,NULL,0,0.0,NULL } };

static AsnType atx[29] = {
  {401, "Severity-level" ,1,0,0,0,0,1,0,0,NULL,&atx[1],&avnx[0],0,&atx[2]} ,
  {310, "ENUMERATED" ,0,10,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {402, "Field-rule" ,1,0,0,0,0,1,0,0,NULL,&atx[9],&atx[3],0,&atx[10]} ,
  {0, "field-name" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[5]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "match-expression" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[6]} ,
  {0, "required" ,128,2,0,0,1,0,0,0,&avnx[6],&atx[7],NULL,0,&atx[8]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "severity" ,128,3,0,0,1,0,0,0,&avnx[7],&atx[0],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {403, "Field-set" ,1,0,0,0,0,1,0,0,NULL,&atx[12],&atx[11],0,&atx[13]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {404, "Dependent-field-rule" ,1,0,0,0,0,1,0,0,NULL,&atx[9],&atx[14],0,&atx[18]} ,
  {0, "match-name" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[15]} ,
  {0, "value-constraint" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[16]} ,
  {0, "other-fields" ,128,2,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[17]} ,
  {0, "disallowed-fields" ,128,3,0,1,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {405, "Dependent-field-set" ,1,0,0,0,0,1,0,0,NULL,&atx[12],&atx[19],0,&atx[20]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[13],NULL,0,NULL} ,
  {406, "Comment-rule" ,1,0,0,0,0,1,0,0,NULL,&atx[9],&atx[21],0,&atx[27]} ,
  {0, "prefix" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[22]} ,
  {0, "updated" ,128,1,0,0,1,0,0,0,&avnx[8],&atx[7],NULL,0,&atx[23]} ,
  {0, "fields" ,128,2,0,0,0,0,0,0,NULL,&atx[10],NULL,0,&atx[24]} ,
  {0, "require-order" ,128,3,0,0,1,0,0,0,&avnx[9],&atx[7],NULL,0,&atx[25]} ,
  {0, "allow-unlisted" ,128,4,0,0,1,0,0,0,&avnx[10],&atx[7],NULL,0,&atx[26]} ,
  {0, "dependent-rules" ,128,5,0,1,0,0,0,0,NULL,&atx[18],NULL,0,NULL} ,
  {407, "Comment-set" ,1,0,0,0,0,1,0,0,NULL,&atx[12],&atx[28],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[20],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-Structured-comment-validation" , "asnvalid.h13",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-Structured-comment-validation
*
**************************************************/

#define SEVERITY_LEVEL &at[0]

#define FIELD_RULE &at[2]
#define FIELD_RULE_field_name &at[3]
#define FIELD_RULE_match_expression &at[5]
#define FIELD_RULE_required &at[6]
#define FIELD_RULE_severity &at[8]

#define FIELD_SET &at[10]
#define FIELD_SET_E &at[11]

#define DEPENDENT_FIELD_RULE &at[13]
#define DEPENDENT_FIELD_RULE_match_name &at[14]
#define FIELD_RULE_value_constraint &at[15]
#define FIELD_RULE_other_fields &at[16]
#define FIELD_RULE_disallowed_fields &at[17]

#define DEPENDENT_FIELD_SET &at[18]
#define DEPENDENT_FIELD_SET_E &at[19]

#define COMMENT_RULE &at[20]
#define COMMENT_RULE_prefix &at[21]
#define COMMENT_RULE_updated &at[22]
#define COMMENT_RULE_fields &at[23]
#define COMMENT_RULE_require_order &at[24]
#define COMMENT_RULE_allow_unlisted &at[25]
#define COMMENT_RULE_dependent_rules &at[26]

#define COMMENT_SET &at[27]
#define COMMENT_SET_E &at[28]
