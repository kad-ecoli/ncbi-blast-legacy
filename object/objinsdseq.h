#ifndef _objinsdseq_ 
#define _objinsdseq_ 

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif


#ifdef __cplusplus
extern "C" { /* } */
#endif


/**************************************************
*
*    Generated objects for Module INSD-INSDSeq
*    Generated using ASNCODE Revision: 6.17 at May 26, 2010 12:38 PM
*
**************************************************/

NLM_EXTERN Boolean LIBCALL
objinsdseqAsnLoad PROTO((void));


/**************************************************
*
*    INSDSet
*
**************************************************/
typedef struct struct_INSDSeq INSDSet;
typedef struct struct_INSDSeq PNTR INSDSetPtr;
#define INSDSetNew() INSDSeqNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN INSDSetPtr LIBCALL INSDSetFree PROTO ((INSDSetPtr ));
NLM_EXTERN INSDSetPtr LIBCALL INSDSetNew PROTO (( void ));
NLM_EXTERN INSDSetPtr LIBCALL INSDSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDSetAsnWrite PROTO (( INSDSetPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    INSDSeq
*
**************************************************/
typedef struct struct_INSDSeq {
   struct struct_INSDSeq PNTR next;
   Uint4 OBbits__;
   CharPtr   locus;
   Int4   length;
   CharPtr   strandedness;
   CharPtr   moltype;
   CharPtr   topology;
   CharPtr   division;
   CharPtr   update_date;
   CharPtr   create_date;
   CharPtr   update_release;
   CharPtr   create_release;
   CharPtr   definition;
   CharPtr   primary_accession;
   CharPtr   entry_version;
   CharPtr   accession_version;
   ValNodePtr   other_seqids;
   ValNodePtr   secondary_accessions;
   CharPtr   project;
   ValNodePtr   keywords;
   CharPtr   segment;
   CharPtr   source;
   CharPtr   organism;
   CharPtr   taxonomy;
   struct struct_INSDReference PNTR   references;
   CharPtr   comment;
   struct struct_INSDComment PNTR   comment_set;
   struct struct_INSDStrucComment PNTR   struc_comments;
   CharPtr   primary;
   CharPtr   source_db;
   CharPtr   database_reference;
   struct struct_INSDFeature PNTR   feature_table;
   struct struct_INSDFeatureSet PNTR   feature_set;
   CharPtr   sequence;
   CharPtr   contig;
   struct struct_INSDAltSeqData PNTR   alt_seq;
} INSDSeq, PNTR INSDSeqPtr;


NLM_EXTERN INSDSeqPtr LIBCALL INSDSeqFree PROTO ((INSDSeqPtr ));
NLM_EXTERN INSDSeqPtr LIBCALL INSDSeqNew PROTO (( void ));
NLM_EXTERN INSDSeqPtr LIBCALL INSDSeqAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDSeqAsnWrite PROTO (( INSDSeqPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDReference
*
**************************************************/
typedef struct struct_INSDReference {
   struct struct_INSDReference PNTR next;
   Uint4 OBbits__;
   CharPtr   reference;
   CharPtr   position;
   ValNodePtr   authors;
   CharPtr   consortium;
   CharPtr   title;
   CharPtr   journal;
   struct struct_INSDXref PNTR   xref;
#define OB__INSDReference_pubmed 0

   Int4   pubmed;
   CharPtr   remark;
} INSDReference, PNTR INSDReferencePtr;


NLM_EXTERN INSDReferencePtr LIBCALL INSDReferenceFree PROTO ((INSDReferencePtr ));
NLM_EXTERN INSDReferencePtr LIBCALL INSDReferenceNew PROTO (( void ));
NLM_EXTERN INSDReferencePtr LIBCALL INSDReferenceAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDReferenceAsnWrite PROTO (( INSDReferencePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDComment
*
**************************************************/
typedef struct struct_INSDComment {
   struct struct_INSDComment PNTR next;
   Uint4 OBbits__;
   CharPtr   type;
   struct struct_INSDCommentParagraph PNTR   paragraphs;
} INSDComment, PNTR INSDCommentPtr;


NLM_EXTERN INSDCommentPtr LIBCALL INSDCommentFree PROTO ((INSDCommentPtr ));
NLM_EXTERN INSDCommentPtr LIBCALL INSDCommentNew PROTO (( void ));
NLM_EXTERN INSDCommentPtr LIBCALL INSDCommentAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDCommentAsnWrite PROTO (( INSDCommentPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDStrucComment
*
**************************************************/
typedef struct struct_INSDStrucComment {
   struct struct_INSDStrucComment PNTR next;
   Uint4 OBbits__;
   CharPtr   name;
   struct struct_INSDStrucCommentItem PNTR   items;
} INSDStrucComment, PNTR INSDStrucCommentPtr;


NLM_EXTERN INSDStrucCommentPtr LIBCALL INSDStrucCommentFree PROTO ((INSDStrucCommentPtr ));
NLM_EXTERN INSDStrucCommentPtr LIBCALL INSDStrucCommentNew PROTO (( void ));
NLM_EXTERN INSDStrucCommentPtr LIBCALL INSDStrucCommentAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDStrucCommentAsnWrite PROTO (( INSDStrucCommentPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDFeature
*
**************************************************/
typedef struct struct_INSDFeature {
   struct struct_INSDFeature PNTR next;
   Uint4 OBbits__;
   CharPtr   key;
   CharPtr   location;
   struct struct_INSDInterval PNTR   intervals;
   CharPtr   operator__;
#define OB__INSDFeature_partial5 0

   Uint1   partial5;
#define OB__INSDFeature_partial3 1

   Uint1   partial3;
   struct struct_INSDQualifier PNTR   quals;
   struct struct_INSDXref PNTR   xrefs;
} INSDFeature, PNTR INSDFeaturePtr;


NLM_EXTERN INSDFeaturePtr LIBCALL INSDFeatureFree PROTO ((INSDFeaturePtr ));
NLM_EXTERN INSDFeaturePtr LIBCALL INSDFeatureNew PROTO (( void ));
NLM_EXTERN INSDFeaturePtr LIBCALL INSDFeatureAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDFeatureAsnWrite PROTO (( INSDFeaturePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDFeatureSet
*
**************************************************/
typedef struct struct_INSDFeatureSet {
   struct struct_INSDFeatureSet PNTR next;
   Uint4 OBbits__;
   CharPtr   annot_source;
   struct struct_INSDFeature PNTR   features;
} INSDFeatureSet, PNTR INSDFeatureSetPtr;


NLM_EXTERN INSDFeatureSetPtr LIBCALL INSDFeatureSetFree PROTO ((INSDFeatureSetPtr ));
NLM_EXTERN INSDFeatureSetPtr LIBCALL INSDFeatureSetNew PROTO (( void ));
NLM_EXTERN INSDFeatureSetPtr LIBCALL INSDFeatureSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDFeatureSetAsnWrite PROTO (( INSDFeatureSetPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDAltSeqData
*
**************************************************/
typedef struct struct_INSDAltSeqData {
   struct struct_INSDAltSeqData PNTR next;
   Uint4 OBbits__;
   CharPtr   name;
   struct struct_INSDAltSeqItem PNTR   items;
} INSDAltSeqData, PNTR INSDAltSeqDataPtr;


NLM_EXTERN INSDAltSeqDataPtr LIBCALL INSDAltSeqDataFree PROTO ((INSDAltSeqDataPtr ));
NLM_EXTERN INSDAltSeqDataPtr LIBCALL INSDAltSeqDataNew PROTO (( void ));
NLM_EXTERN INSDAltSeqDataPtr LIBCALL INSDAltSeqDataAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDAltSeqDataAsnWrite PROTO (( INSDAltSeqDataPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDXref
*
**************************************************/
typedef struct struct_INSDXref {
   struct struct_INSDXref PNTR next;
   Uint4 OBbits__;
   CharPtr   dbname;
   CharPtr   id;
} INSDXref, PNTR INSDXrefPtr;


NLM_EXTERN INSDXrefPtr LIBCALL INSDXrefFree PROTO ((INSDXrefPtr ));
NLM_EXTERN INSDXrefPtr LIBCALL INSDXrefNew PROTO (( void ));
NLM_EXTERN INSDXrefPtr LIBCALL INSDXrefAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDXrefAsnWrite PROTO (( INSDXrefPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDCommentParagraph
*
**************************************************/
typedef struct struct_INSDCommentParagraph {
   struct struct_INSDCommentParagraph PNTR next;
   Uint4 OBbits__;
   struct struct_INSDCommentItem PNTR   items;
} INSDCommentParagraph, PNTR INSDCommentParagraphPtr;


NLM_EXTERN INSDCommentParagraphPtr LIBCALL INSDCommentParagraphFree PROTO ((INSDCommentParagraphPtr ));
NLM_EXTERN INSDCommentParagraphPtr LIBCALL INSDCommentParagraphNew PROTO (( void ));
NLM_EXTERN INSDCommentParagraphPtr LIBCALL INSDCommentParagraphAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDCommentParagraphAsnWrite PROTO (( INSDCommentParagraphPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDCommentItem
*
**************************************************/
typedef struct struct_INSDCommentItem {
   struct struct_INSDCommentItem PNTR next;
   Uint4 OBbits__;
   CharPtr   value;
   CharPtr   url;
} INSDCommentItem, PNTR INSDCommentItemPtr;


NLM_EXTERN INSDCommentItemPtr LIBCALL INSDCommentItemFree PROTO ((INSDCommentItemPtr ));
NLM_EXTERN INSDCommentItemPtr LIBCALL INSDCommentItemNew PROTO (( void ));
NLM_EXTERN INSDCommentItemPtr LIBCALL INSDCommentItemAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDCommentItemAsnWrite PROTO (( INSDCommentItemPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDStrucCommentItem
*
**************************************************/
typedef struct struct_INSDStrucCommentItem {
   struct struct_INSDStrucCommentItem PNTR next;
   Uint4 OBbits__;
   CharPtr   tag;
   CharPtr   value;
   CharPtr   url;
} INSDStrucCommentItem, PNTR INSDStrucCommentItemPtr;


NLM_EXTERN INSDStrucCommentItemPtr LIBCALL INSDStrucCommentItemFree PROTO ((INSDStrucCommentItemPtr ));
NLM_EXTERN INSDStrucCommentItemPtr LIBCALL INSDStrucCommentItemNew PROTO (( void ));
NLM_EXTERN INSDStrucCommentItemPtr LIBCALL INSDStrucCommentItemAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDStrucCommentItemAsnWrite PROTO (( INSDStrucCommentItemPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDInterval
*
**************************************************/
typedef struct struct_INSDInterval {
   struct struct_INSDInterval PNTR next;
   Uint4 OBbits__;
#define OB__INSDInterval_from 0

   Int4   from;
#define OB__INSDInterval_to 1

   Int4   to;
#define OB__INSDInterval_point 2

   Int4   point;
#define OB__INSDInterval_iscomp 3

   Uint1   iscomp;
#define OB__INSDInterval_interbp 4

   Uint1   interbp;
   CharPtr   accession;
} INSDInterval, PNTR INSDIntervalPtr;


NLM_EXTERN INSDIntervalPtr LIBCALL INSDIntervalFree PROTO ((INSDIntervalPtr ));
NLM_EXTERN INSDIntervalPtr LIBCALL INSDIntervalNew PROTO (( void ));
NLM_EXTERN INSDIntervalPtr LIBCALL INSDIntervalAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDIntervalAsnWrite PROTO (( INSDIntervalPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDQualifier
*
**************************************************/
typedef struct struct_INSDQualifier {
   struct struct_INSDQualifier PNTR next;
   Uint4 OBbits__;
   CharPtr   name;
   CharPtr   value;
} INSDQualifier, PNTR INSDQualifierPtr;


NLM_EXTERN INSDQualifierPtr LIBCALL INSDQualifierFree PROTO ((INSDQualifierPtr ));
NLM_EXTERN INSDQualifierPtr LIBCALL INSDQualifierNew PROTO (( void ));
NLM_EXTERN INSDQualifierPtr LIBCALL INSDQualifierAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDQualifierAsnWrite PROTO (( INSDQualifierPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    INSDAltSeqItem
*
**************************************************/
typedef struct struct_INSDAltSeqItem {
   struct struct_INSDAltSeqItem PNTR next;
   Uint4 OBbits__;
   struct struct_INSDInterval PNTR   interval;
#define OB__INSDAltSeqItem_isgap 0

   Uint1   isgap;
#define OB__INSDAltSeqItem_gap_length 1

   Int4   gap_length;
   CharPtr   gap_type;
   CharPtr   gap_linkage;
   CharPtr   gap_comment;
   CharPtr   first_accn;
   CharPtr   last_accn;
   CharPtr   value;
} INSDAltSeqItem, PNTR INSDAltSeqItemPtr;


NLM_EXTERN INSDAltSeqItemPtr LIBCALL INSDAltSeqItemFree PROTO ((INSDAltSeqItemPtr ));
NLM_EXTERN INSDAltSeqItemPtr LIBCALL INSDAltSeqItemNew PROTO (( void ));
NLM_EXTERN INSDAltSeqItemPtr LIBCALL INSDAltSeqItemAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL INSDAltSeqItemAsnWrite PROTO (( INSDAltSeqItemPtr , AsnIoPtr, AsnTypePtr));

#ifdef __cplusplus
/* { */ }
#endif

#endif /* _objinsdseq_ */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

