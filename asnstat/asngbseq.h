/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asngbseq.h10";
static AsnType atx[128] = {
  {401, "GBSet" ,1,0,0,0,0,0,0,0,NULL,&atx[22],&atx[1],0,&atx[2]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {402, "GBSeq" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[3],0,&atx[21]} ,
  {0, "locus" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[5]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "length" ,128,1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[7]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "strandedness" ,128,2,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[8]} ,
  {0, "moltype" ,128,3,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[9]} ,
  {0, "topology" ,128,4,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[10]} ,
  {0, "division" ,128,5,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[11]} ,
  {0, "update-date" ,128,6,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[12]} ,
  {0, "create-date" ,128,7,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[13]} ,
  {0, "update-release" ,128,8,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[14]} ,
  {0, "create-release" ,128,9,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[15]} ,
  {0, "definition" ,128,10,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[16]} ,
  {0, "primary-accession" ,128,11,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[17]} ,
  {0, "entry-version" ,128,12,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[18]} ,
  {0, "accession-version" ,128,13,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[19]} ,
  {0, "other-seqids" ,128,14,0,1,0,0,0,0,NULL,&atx[22],&atx[20],0,&atx[23]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[21],NULL,0,NULL} ,
  {403, "GBSeqid" ,1,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[25]} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "secondary-accessions" ,128,15,0,1,0,0,0,0,NULL,&atx[22],&atx[24],0,&atx[26]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[25],NULL,0,NULL} ,
  {404, "GBSecondary-accn" ,1,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[29]} ,
  {0, "project" ,128,16,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[27]} ,
  {0, "keywords" ,128,17,0,1,0,0,0,0,NULL,&atx[22],&atx[28],0,&atx[30]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[29],NULL,0,NULL} ,
  {405, "GBKeyword" ,1,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[36]} ,
  {0, "segment" ,128,18,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[31]} ,
  {0, "source" ,128,19,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[32]} ,
  {0, "organism" ,128,20,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[33]} ,
  {0, "taxonomy" ,128,21,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[34]} ,
  {0, "references" ,128,22,0,1,0,0,0,0,NULL,&atx[22],&atx[35],0,&atx[53]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[36],NULL,0,NULL} ,
  {406, "GBReference" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[37],0,&atx[56]} ,
  {0, "reference" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[38]} ,
  {0, "position" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[39]} ,
  {0, "authors" ,128,2,0,1,0,0,0,0,NULL,&atx[22],&atx[40],0,&atx[42]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[41],NULL,0,NULL} ,
  {412, "GBAuthor" ,1,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[47]} ,
  {0, "consortium" ,128,3,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[43]} ,
  {0, "title" ,128,4,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[44]} ,
  {0, "journal" ,128,5,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[45]} ,
  {0, "xref" ,128,6,0,1,0,0,0,0,NULL,&atx[22],&atx[46],0,&atx[51]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {413, "GBXref" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[48],0,&atx[60]} ,
  {0, "dbname" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[49]} ,
  {0, "id" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "pubmed" ,128,7,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[52]} ,
  {0, "remark" ,128,8,0,1,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "comment" ,128,23,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[54]} ,
  {0, "comment-set" ,128,24,0,1,0,0,0,0,NULL,&atx[22],&atx[55],0,&atx[66]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[56],NULL,0,NULL} ,
  {407, "GBComment" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[57],0,&atx[68]} ,
  {0, "type" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[58]} ,
  {0, "paragraphs" ,128,1,0,0,0,0,0,0,NULL,&atx[22],&atx[59],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[60],NULL,0,NULL} ,
  {414, "GBCommentParagraph" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[61],0,&atx[63]} ,
  {0, "items" ,128,0,0,0,0,0,0,0,NULL,&atx[22],&atx[62],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[63],NULL,0,NULL} ,
  {415, "GBCommentItem" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[64],0,&atx[72]} ,
  {0, "value" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[65]} ,
  {0, "url" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "struc-comments" ,128,25,0,1,0,0,0,0,NULL,&atx[22],&atx[67],0,&atx[76]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[68],NULL,0,NULL} ,
  {408, "GBStrucComment" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[69],0,&atx[81]} ,
  {0, "name" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[70]} ,
  {0, "items" ,128,1,0,0,0,0,0,0,NULL,&atx[22],&atx[71],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[72],NULL,0,NULL} ,
  {416, "GBStrucCommentItem" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[73],0,&atx[86]} ,
  {0, "tag" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[74]} ,
  {0, "value" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[75]} ,
  {0, "url" ,128,2,0,1,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "primary" ,128,26,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[77]} ,
  {0, "source-db" ,128,27,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[78]} ,
  {0, "database-reference" ,128,28,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[79]} ,
  {0, "feature-table" ,128,29,0,1,0,0,0,0,NULL,&atx[22],&atx[80],0,&atx[104]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[81],NULL,0,NULL} ,
  {409, "GBFeature" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[82],0,&atx[106]} ,
  {0, "key" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[83]} ,
  {0, "location" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[84]} ,
  {0, "intervals" ,128,2,0,1,0,0,0,0,NULL,&atx[22],&atx[85],0,&atx[94]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[86],NULL,0,NULL} ,
  {417, "GBInterval" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[87],0,&atx[99]} ,
  {0, "from" ,128,0,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[88]} ,
  {0, "to" ,128,1,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[89]} ,
  {0, "point" ,128,2,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[90]} ,
  {0, "iscomp" ,128,3,0,1,0,0,0,0,NULL,&atx[91],NULL,0,&atx[92]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "interbp" ,128,4,0,1,0,0,0,0,NULL,&atx[91],NULL,0,&atx[93]} ,
  {0, "accession" ,128,5,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "operator" ,128,3,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[95]} ,
  {0, "partial5" ,128,4,0,1,0,0,0,0,NULL,&atx[91],NULL,0,&atx[96]} ,
  {0, "partial3" ,128,5,0,1,0,0,0,0,NULL,&atx[91],NULL,0,&atx[97]} ,
  {0, "quals" ,128,6,0,1,0,0,0,0,NULL,&atx[22],&atx[98],0,&atx[102]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[99],NULL,0,NULL} ,
  {418, "GBQualifier" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[100],0,&atx[118]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[101]} ,
  {0, "value" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "xrefs" ,128,7,0,1,0,0,0,0,NULL,&atx[22],&atx[103],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {0, "feature-set" ,128,30,0,1,0,0,0,0,NULL,&atx[22],&atx[105],0,&atx[110]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[106],NULL,0,NULL} ,
  {410, "GBFeatureSet" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[107],0,&atx[114]} ,
  {0, "annot-source" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[108]} ,
  {0, "features" ,128,1,0,0,0,0,0,0,NULL,&atx[22],&atx[109],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[81],NULL,0,NULL} ,
  {0, "sequence" ,128,31,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[111]} ,
  {0, "contig" ,128,32,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[112]} ,
  {0, "alt-seq" ,128,33,0,1,0,0,0,0,NULL,&atx[22],&atx[113],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[114],NULL,0,NULL} ,
  {411, "GBAltSeqData" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[115],0,&atx[41]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[116]} ,
  {0, "items" ,128,1,0,1,0,0,0,0,NULL,&atx[22],&atx[117],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[118],NULL,0,NULL} ,
  {419, "GBAltSeqItem" ,1,0,0,0,0,0,0,0,NULL,&atx[50],&atx[119],0,NULL} ,
  {0, "interval" ,128,0,0,1,0,0,0,0,NULL,&atx[86],NULL,0,&atx[120]} ,
  {0, "isgap" ,128,1,0,1,0,0,0,0,NULL,&atx[91],NULL,0,&atx[121]} ,
  {0, "gap-length" ,128,2,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[122]} ,
  {0, "gap-type" ,128,3,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[123]} ,
  {0, "gap-linkage" ,128,4,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[124]} ,
  {0, "gap-comment" ,128,5,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[125]} ,
  {0, "first-accn" ,128,6,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[126]} ,
  {0, "last-accn" ,128,7,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[127]} ,
  {0, "value" ,128,8,0,1,0,0,0,0,NULL,&atx[4],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-GBSeq" , "asngbseq.h10",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = NULL;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-GBSeq
*
**************************************************/

#define GBSET &at[0]
#define GBSET_E &at[1]

#define GBSEQ &at[2]
#define GBSEQ_locus &at[3]
#define GBSEQ_length &at[5]
#define GBSEQ_strandedness &at[7]
#define GBSEQ_moltype &at[8]
#define GBSEQ_topology &at[9]
#define GBSEQ_division &at[10]
#define GBSEQ_update_date &at[11]
#define GBSEQ_create_date &at[12]
#define GBSEQ_update_release &at[13]
#define GBSEQ_create_release &at[14]
#define GBSEQ_definition &at[15]
#define GBSEQ_primary_accession &at[16]
#define GBSEQ_entry_version &at[17]
#define GBSEQ_accession_version &at[18]
#define GBSEQ_other_seqids &at[19]
#define GBSEQ_other_seqids_E &at[20]
#define GBSEQ_secondary_accessions &at[23]
#define GBSEQ_secondary_accessions_E &at[24]
#define GBSEQ_project &at[26]
#define GBSEQ_keywords &at[27]
#define GBSEQ_keywords_E &at[28]
#define GBSEQ_segment &at[30]
#define GBSEQ_source &at[31]
#define GBSEQ_organism &at[32]
#define GBSEQ_taxonomy &at[33]
#define GBSEQ_references &at[34]
#define GBSEQ_references_E &at[35]
#define GBSEQ_comment &at[53]
#define GBSEQ_comment_set &at[54]
#define GBSEQ_comment_set_E &at[55]
#define GBSEQ_struc_comments &at[66]
#define GBSEQ_struc_comments_E &at[67]
#define GBSEQ_primary &at[76]
#define GBSEQ_source_db &at[77]
#define GBSEQ_database_reference &at[78]
#define GBSEQ_feature_table &at[79]
#define GBSEQ_feature_table_E &at[80]
#define GBSEQ_feature_set &at[104]
#define GBSEQ_feature_set_E &at[105]
#define GBSEQ_sequence &at[110]
#define GBSEQ_contig &at[111]
#define GBSEQ_alt_seq &at[112]
#define GBSEQ_alt_seq_E &at[113]

#define GBSEQID &at[21]

#define GBSECONDARY_ACCN &at[25]

#define GBKEYWORD &at[29]

#define GBREFERENCE &at[36]
#define GBREFERENCE_reference &at[37]
#define GBREFERENCE_position &at[38]
#define GBREFERENCE_authors &at[39]
#define GBREFERENCE_authors_E &at[40]
#define GBREFERENCE_consortium &at[42]
#define GBREFERENCE_title &at[43]
#define GBREFERENCE_journal &at[44]
#define GBREFERENCE_xref &at[45]
#define GBREFERENCE_xref_E &at[46]
#define GBREFERENCE_pubmed &at[51]
#define GBREFERENCE_remark &at[52]

#define GBCOMMENT &at[56]
#define GBCOMMENT_type &at[57]
#define GBCOMMENT_paragraphs &at[58]
#define GBCOMMENT_paragraphs_E &at[59]

#define GBSTRUCCOMMENT &at[68]
#define GBSTRUCCOMMENT_name &at[69]
#define GBSTRUCCOMMENT_items &at[70]
#define GBSTRUCCOMMENT_items_E &at[71]

#define GBFEATURE &at[81]
#define GBFEATURE_key &at[82]
#define GBFEATURE_location &at[83]
#define GBFEATURE_intervals &at[84]
#define GBFEATURE_intervals_E &at[85]
#define GBFEATURE_operator &at[94]
#define GBFEATURE_partial5 &at[95]
#define GBFEATURE_partial3 &at[96]
#define GBFEATURE_quals &at[97]
#define GBFEATURE_quals_E &at[98]
#define GBFEATURE_xrefs &at[102]
#define GBFEATURE_xrefs_E &at[103]

#define GBFEATURESET &at[106]
#define GBFEATURESET_annot_source &at[107]
#define GBFEATURESET_features &at[108]
#define GBFEATURESET_features_E &at[109]

#define GBALTSEQDATA &at[114]
#define GBALTSEQDATA_name &at[115]
#define GBALTSEQDATA_items &at[116]
#define GBALTSEQDATA_items_E &at[117]

#define GBAUTHOR &at[41]

#define GBXREF &at[47]
#define GBXREF_dbname &at[48]
#define GBXREF_id &at[49]

#define GBCOMMENTPARAGRAPH &at[60]
#define GBCOMMENTPARAGRAPH_items &at[61]
#define GBCOMMENTPARAGRAPH_items_E &at[62]

#define GBCOMMENTITEM &at[63]
#define GBCOMMENTITEM_value &at[64]
#define GBCOMMENTITEM_url &at[65]

#define GBSTRUCCOMMENTITEM &at[72]
#define GBSTRUCCOMMENTITEM_tag &at[73]
#define GBSTRUCCOMMENTITEM_value &at[74]
#define GBSTRUCCOMMENTITEM_url &at[75]

#define GBINTERVAL &at[86]
#define GBINTERVAL_from &at[87]
#define GBINTERVAL_to &at[88]
#define GBINTERVAL_point &at[89]
#define GBINTERVAL_iscomp &at[90]
#define GBINTERVAL_interbp &at[92]
#define GBINTERVAL_accession &at[93]

#define GBQUALIFIER &at[99]
#define GBQUALIFIER_name &at[100]
#define GBQUALIFIER_value &at[101]

#define GBALTSEQITEM &at[118]
#define GBALTSEQITEM_interval &at[119]
#define GBALTSEQITEM_isgap &at[120]
#define GBALTSEQITEM_gap_length &at[121]
#define GBALTSEQITEM_gap_type &at[122]
#define GBALTSEQITEM_gap_linkage &at[123]
#define GBALTSEQITEM_gap_comment &at[124]
#define GBALTSEQITEM_first_accn &at[125]
#define GBALTSEQITEM_last_accn &at[126]
#define GBALTSEQITEM_value &at[127]
