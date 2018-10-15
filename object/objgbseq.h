#ifndef _objgbseq_ 
#define _objgbseq_ 

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
*    Generated objects for Module NCBI-GBSeq
*    Generated using ASNCODE Revision: 6.17 at May 26, 2010 12:37 PM
*
**************************************************/

NLM_EXTERN Boolean LIBCALL
objgbseqAsnLoad PROTO((void));


/**************************************************
*
*    GBSet
*
**************************************************/
typedef struct struct_GBSeq GBSet;
typedef struct struct_GBSeq PNTR GBSetPtr;
#define GBSetNew() GBSeqNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN GBSetPtr LIBCALL GBSetFree PROTO ((GBSetPtr ));
NLM_EXTERN GBSetPtr LIBCALL GBSetNew PROTO (( void ));
NLM_EXTERN GBSetPtr LIBCALL GBSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBSetAsnWrite PROTO (( GBSetPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    GBSeq
*
**************************************************/
typedef struct struct_GBSeq {
   struct struct_GBSeq PNTR next;
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
   struct struct_GBReference PNTR   references;
   CharPtr   comment;
   struct struct_GBComment PNTR   comment_set;
   struct struct_GBStrucComment PNTR   struc_comments;
   CharPtr   primary;
   CharPtr   source_db;
   CharPtr   database_reference;
   struct struct_GBFeature PNTR   feature_table;
   struct struct_GBFeatureSet PNTR   feature_set;
   CharPtr   sequence;
   CharPtr   contig;
   struct struct_GBAltSeqData PNTR   alt_seq;
} GBSeq, PNTR GBSeqPtr;


NLM_EXTERN GBSeqPtr LIBCALL GBSeqFree PROTO ((GBSeqPtr ));
NLM_EXTERN GBSeqPtr LIBCALL GBSeqNew PROTO (( void ));
NLM_EXTERN GBSeqPtr LIBCALL GBSeqAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBSeqAsnWrite PROTO (( GBSeqPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBReference
*
**************************************************/
typedef struct struct_GBReference {
   struct struct_GBReference PNTR next;
   Uint4 OBbits__;
   CharPtr   reference;
   CharPtr   position;
   ValNodePtr   authors;
   CharPtr   consortium;
   CharPtr   title;
   CharPtr   journal;
   struct struct_GBXref PNTR   xref;
#define OB__GBReference_pubmed 0

   Int4   pubmed;
   CharPtr   remark;
} GBReference, PNTR GBReferencePtr;


NLM_EXTERN GBReferencePtr LIBCALL GBReferenceFree PROTO ((GBReferencePtr ));
NLM_EXTERN GBReferencePtr LIBCALL GBReferenceNew PROTO (( void ));
NLM_EXTERN GBReferencePtr LIBCALL GBReferenceAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBReferenceAsnWrite PROTO (( GBReferencePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBComment
*
**************************************************/
typedef struct struct_GBComment {
   struct struct_GBComment PNTR next;
   Uint4 OBbits__;
   CharPtr   type;
   struct struct_GBCommentParagraph PNTR   paragraphs;
} GBComment, PNTR GBCommentPtr;


NLM_EXTERN GBCommentPtr LIBCALL GBCommentFree PROTO ((GBCommentPtr ));
NLM_EXTERN GBCommentPtr LIBCALL GBCommentNew PROTO (( void ));
NLM_EXTERN GBCommentPtr LIBCALL GBCommentAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBCommentAsnWrite PROTO (( GBCommentPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBStrucComment
*
**************************************************/
typedef struct struct_GBStrucComment {
   struct struct_GBStrucComment PNTR next;
   Uint4 OBbits__;
   CharPtr   name;
   struct struct_GBStrucCommentItem PNTR   items;
} GBStrucComment, PNTR GBStrucCommentPtr;


NLM_EXTERN GBStrucCommentPtr LIBCALL GBStrucCommentFree PROTO ((GBStrucCommentPtr ));
NLM_EXTERN GBStrucCommentPtr LIBCALL GBStrucCommentNew PROTO (( void ));
NLM_EXTERN GBStrucCommentPtr LIBCALL GBStrucCommentAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBStrucCommentAsnWrite PROTO (( GBStrucCommentPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBFeature
*
**************************************************/
typedef struct struct_GBFeature {
   struct struct_GBFeature PNTR next;
   Uint4 OBbits__;
   CharPtr   key;
   CharPtr   location;
   struct struct_GBInterval PNTR   intervals;
   CharPtr   operator__;
#define OB__GBFeature_partial5 0

   Uint1   partial5;
#define OB__GBFeature_partial3 1

   Uint1   partial3;
   struct struct_GBQualifier PNTR   quals;
   struct struct_GBXref PNTR   xrefs;
} GBFeature, PNTR GBFeaturePtr;


NLM_EXTERN GBFeaturePtr LIBCALL GBFeatureFree PROTO ((GBFeaturePtr ));
NLM_EXTERN GBFeaturePtr LIBCALL GBFeatureNew PROTO (( void ));
NLM_EXTERN GBFeaturePtr LIBCALL GBFeatureAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBFeatureAsnWrite PROTO (( GBFeaturePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBFeatureSet
*
**************************************************/
typedef struct struct_GBFeatureSet {
   struct struct_GBFeatureSet PNTR next;
   Uint4 OBbits__;
   CharPtr   annot_source;
   struct struct_GBFeature PNTR   features;
} GBFeatureSet, PNTR GBFeatureSetPtr;


NLM_EXTERN GBFeatureSetPtr LIBCALL GBFeatureSetFree PROTO ((GBFeatureSetPtr ));
NLM_EXTERN GBFeatureSetPtr LIBCALL GBFeatureSetNew PROTO (( void ));
NLM_EXTERN GBFeatureSetPtr LIBCALL GBFeatureSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBFeatureSetAsnWrite PROTO (( GBFeatureSetPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBAltSeqData
*
**************************************************/
typedef struct struct_GBAltSeqData {
   struct struct_GBAltSeqData PNTR next;
   Uint4 OBbits__;
   CharPtr   name;
   struct struct_GBAltSeqItem PNTR   items;
} GBAltSeqData, PNTR GBAltSeqDataPtr;


NLM_EXTERN GBAltSeqDataPtr LIBCALL GBAltSeqDataFree PROTO ((GBAltSeqDataPtr ));
NLM_EXTERN GBAltSeqDataPtr LIBCALL GBAltSeqDataNew PROTO (( void ));
NLM_EXTERN GBAltSeqDataPtr LIBCALL GBAltSeqDataAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBAltSeqDataAsnWrite PROTO (( GBAltSeqDataPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBXref
*
**************************************************/
typedef struct struct_GBXref {
   struct struct_GBXref PNTR next;
   Uint4 OBbits__;
   CharPtr   dbname;
   CharPtr   id;
} GBXref, PNTR GBXrefPtr;


NLM_EXTERN GBXrefPtr LIBCALL GBXrefFree PROTO ((GBXrefPtr ));
NLM_EXTERN GBXrefPtr LIBCALL GBXrefNew PROTO (( void ));
NLM_EXTERN GBXrefPtr LIBCALL GBXrefAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBXrefAsnWrite PROTO (( GBXrefPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBCommentParagraph
*
**************************************************/
typedef struct struct_GBCommentParagraph {
   struct struct_GBCommentParagraph PNTR next;
   Uint4 OBbits__;
   struct struct_GBCommentItem PNTR   items;
} GBCommentParagraph, PNTR GBCommentParagraphPtr;


NLM_EXTERN GBCommentParagraphPtr LIBCALL GBCommentParagraphFree PROTO ((GBCommentParagraphPtr ));
NLM_EXTERN GBCommentParagraphPtr LIBCALL GBCommentParagraphNew PROTO (( void ));
NLM_EXTERN GBCommentParagraphPtr LIBCALL GBCommentParagraphAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBCommentParagraphAsnWrite PROTO (( GBCommentParagraphPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBCommentItem
*
**************************************************/
typedef struct struct_GBCommentItem {
   struct struct_GBCommentItem PNTR next;
   Uint4 OBbits__;
   CharPtr   value;
   CharPtr   url;
} GBCommentItem, PNTR GBCommentItemPtr;


NLM_EXTERN GBCommentItemPtr LIBCALL GBCommentItemFree PROTO ((GBCommentItemPtr ));
NLM_EXTERN GBCommentItemPtr LIBCALL GBCommentItemNew PROTO (( void ));
NLM_EXTERN GBCommentItemPtr LIBCALL GBCommentItemAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBCommentItemAsnWrite PROTO (( GBCommentItemPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBStrucCommentItem
*
**************************************************/
typedef struct struct_GBStrucCommentItem {
   struct struct_GBStrucCommentItem PNTR next;
   Uint4 OBbits__;
   CharPtr   tag;
   CharPtr   value;
   CharPtr   url;
} GBStrucCommentItem, PNTR GBStrucCommentItemPtr;


NLM_EXTERN GBStrucCommentItemPtr LIBCALL GBStrucCommentItemFree PROTO ((GBStrucCommentItemPtr ));
NLM_EXTERN GBStrucCommentItemPtr LIBCALL GBStrucCommentItemNew PROTO (( void ));
NLM_EXTERN GBStrucCommentItemPtr LIBCALL GBStrucCommentItemAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBStrucCommentItemAsnWrite PROTO (( GBStrucCommentItemPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBInterval
*
**************************************************/
typedef struct struct_GBInterval {
   struct struct_GBInterval PNTR next;
   Uint4 OBbits__;
#define OB__GBInterval_from 0

   Int4   from;
#define OB__GBInterval_to 1

   Int4   to;
#define OB__GBInterval_point 2

   Int4   point;
#define OB__GBInterval_iscomp 3

   Uint1   iscomp;
#define OB__GBInterval_interbp 4

   Uint1   interbp;
   CharPtr   accession;
} GBInterval, PNTR GBIntervalPtr;


NLM_EXTERN GBIntervalPtr LIBCALL GBIntervalFree PROTO ((GBIntervalPtr ));
NLM_EXTERN GBIntervalPtr LIBCALL GBIntervalNew PROTO (( void ));
NLM_EXTERN GBIntervalPtr LIBCALL GBIntervalAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBIntervalAsnWrite PROTO (( GBIntervalPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBQualifier
*
**************************************************/
typedef struct struct_GBQualifier {
   struct struct_GBQualifier PNTR next;
   Uint4 OBbits__;
   CharPtr   name;
   CharPtr   value;
} GBQualifier, PNTR GBQualifierPtr;


NLM_EXTERN GBQualifierPtr LIBCALL GBQualifierFree PROTO ((GBQualifierPtr ));
NLM_EXTERN GBQualifierPtr LIBCALL GBQualifierNew PROTO (( void ));
NLM_EXTERN GBQualifierPtr LIBCALL GBQualifierAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBQualifierAsnWrite PROTO (( GBQualifierPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GBAltSeqItem
*
**************************************************/
typedef struct struct_GBAltSeqItem {
   struct struct_GBAltSeqItem PNTR next;
   Uint4 OBbits__;
   struct struct_GBInterval PNTR   interval;
#define OB__GBAltSeqItem_isgap 0

   Uint1   isgap;
#define OB__GBAltSeqItem_gap_length 1

   Int4   gap_length;
   CharPtr   gap_type;
   CharPtr   gap_linkage;
   CharPtr   gap_comment;
   CharPtr   first_accn;
   CharPtr   last_accn;
   CharPtr   value;
} GBAltSeqItem, PNTR GBAltSeqItemPtr;


NLM_EXTERN GBAltSeqItemPtr LIBCALL GBAltSeqItemFree PROTO ((GBAltSeqItemPtr ));
NLM_EXTERN GBAltSeqItemPtr LIBCALL GBAltSeqItemNew PROTO (( void ));
NLM_EXTERN GBAltSeqItemPtr LIBCALL GBAltSeqItemAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GBAltSeqItemAsnWrite PROTO (( GBAltSeqItemPtr , AsnIoPtr, AsnTypePtr));

#ifdef __cplusplus
/* { */ }
#endif

#endif /* _objgbseq_ */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

