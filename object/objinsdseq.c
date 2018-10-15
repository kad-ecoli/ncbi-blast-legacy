#include <asn.h>

#define NLM_GENERATED_CODE_PROTO

#include <objinsdseq.h>

static Boolean loaded = FALSE;

#include <asninsdseq.h>

#ifndef NLM_EXTERN_LOADS
#define NLM_EXTERN_LOADS {}
#endif

NLM_EXTERN Boolean LIBCALL
objinsdseqAsnLoad(void)
{

   if ( ! loaded) {
      NLM_EXTERN_LOADS

      if ( ! AsnLoad ())
      return FALSE;
      loaded = TRUE;
   }

   return TRUE;
}



/**************************************************
*    Generated object loaders for Module INSD-INSDSeq
*    Generated using ASNCODE Revision: 6.17 at May 26, 2010 12:38 PM
*
**************************************************/


/**************************************************
*
*    INSDSetFree()
*
**************************************************/
NLM_EXTERN 
INSDSetPtr LIBCALL
INSDSetFree(INSDSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) INSDSeqFree);
   return NULL;
}


/**************************************************
*
*    INSDSetAsnRead()
*
**************************************************/
NLM_EXTERN 
INSDSetPtr LIBCALL
INSDSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   INSDSetPtr ptr;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* INSDSet ::= (self contained) */
      atp = AsnReadId(aip, amp, INSDSET);
   } else {
      atp = AsnLinkType(orig, INSDSET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) INSDSeqAsnRead, (AsnOptFreeFunc) INSDSeqFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = INSDSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    INSDSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
INSDSetAsnWrite(INSDSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, INSDSET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) INSDSeqAsnWrite, aip, atp, INSDSET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    INSDSeqNew()
*
**************************************************/
NLM_EXTERN 
INSDSeqPtr LIBCALL
INSDSeqNew(void)
{
   INSDSeqPtr ptr = MemNew((size_t) sizeof(INSDSeq));

   return ptr;

}


/**************************************************
*
*    INSDSeqFree()
*
**************************************************/
NLM_EXTERN 
INSDSeqPtr LIBCALL
INSDSeqFree(INSDSeqPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> locus);
   MemFree(ptr -> strandedness);
   MemFree(ptr -> moltype);
   MemFree(ptr -> topology);
   MemFree(ptr -> division);
   MemFree(ptr -> update_date);
   MemFree(ptr -> create_date);
   MemFree(ptr -> update_release);
   MemFree(ptr -> create_release);
   MemFree(ptr -> definition);
   MemFree(ptr -> primary_accession);
   MemFree(ptr -> entry_version);
   MemFree(ptr -> accession_version);
   AsnGenericBaseSeqOfFree(ptr -> other_seqids ,ASNCODE_PTRVAL_SLOT);
   AsnGenericBaseSeqOfFree(ptr -> secondary_accessions ,ASNCODE_PTRVAL_SLOT);
   MemFree(ptr -> project);
   AsnGenericBaseSeqOfFree(ptr -> keywords ,ASNCODE_PTRVAL_SLOT);
   MemFree(ptr -> segment);
   MemFree(ptr -> source);
   MemFree(ptr -> organism);
   MemFree(ptr -> taxonomy);
   AsnGenericUserSeqOfFree(ptr -> references, (AsnOptFreeFunc) INSDReferenceFree);
   MemFree(ptr -> comment);
   AsnGenericUserSeqOfFree(ptr -> comment_set, (AsnOptFreeFunc) INSDCommentFree);
   AsnGenericUserSeqOfFree(ptr -> struc_comments, (AsnOptFreeFunc) INSDStrucCommentFree);
   MemFree(ptr -> primary);
   MemFree(ptr -> source_db);
   MemFree(ptr -> database_reference);
   AsnGenericUserSeqOfFree(ptr -> feature_table, (AsnOptFreeFunc) INSDFeatureFree);
   AsnGenericUserSeqOfFree(ptr -> feature_set, (AsnOptFreeFunc) INSDFeatureSetFree);
   MemFree(ptr -> sequence);
   MemFree(ptr -> contig);
   AsnGenericUserSeqOfFree(ptr -> alt_seq, (AsnOptFreeFunc) INSDAltSeqDataFree);
   return MemFree(ptr);
}


/**************************************************
*
*    INSDSeqAsnRead()
*
**************************************************/
NLM_EXTERN 
INSDSeqPtr LIBCALL
INSDSeqAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   INSDSeqPtr ptr;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* INSDSeq ::= (self contained) */
      atp = AsnReadId(aip, amp, INSDSEQ);
   } else {
      atp = AsnLinkType(orig, INSDSEQ);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = INSDSeqNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == INSDSEQ_locus) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> locus = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_length) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> length = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_strandedness) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> strandedness = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_moltype) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> moltype = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_topology) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> topology = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_division) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> division = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_update_date) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> update_date = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_create_date) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> create_date = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_update_release) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> update_release = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_create_release) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> create_release = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_definition) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> definition = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_primary_accession) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> primary_accession = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_entry_version) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> entry_version = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_accession_version) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> accession_version = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_other_seqids) {
      ptr -> other_seqids = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> other_seqids == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_secondary_accessions) {
      ptr -> secondary_accessions = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> secondary_accessions == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_project) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> project = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_keywords) {
      ptr -> keywords = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> keywords == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_segment) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> segment = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_source) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> source = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_organism) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> organism = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_taxonomy) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> taxonomy = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_references) {
      ptr -> references = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) INSDReferenceAsnRead, (AsnOptFreeFunc) INSDReferenceFree);
      if (isError && ptr -> references == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_comment) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> comment = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_comment_set) {
      ptr -> comment_set = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) INSDCommentAsnRead, (AsnOptFreeFunc) INSDCommentFree);
      if (isError && ptr -> comment_set == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_struc_comments) {
      ptr -> struc_comments = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) INSDStrucCommentAsnRead, (AsnOptFreeFunc) INSDStrucCommentFree);
      if (isError && ptr -> struc_comments == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_primary) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> primary = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_source_db) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> source_db = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_database_reference) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> database_reference = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_feature_table) {
      ptr -> feature_table = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) INSDFeatureAsnRead, (AsnOptFreeFunc) INSDFeatureFree);
      if (isError && ptr -> feature_table == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_feature_set) {
      ptr -> feature_set = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) INSDFeatureSetAsnRead, (AsnOptFreeFunc) INSDFeatureSetFree);
      if (isError && ptr -> feature_set == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_sequence) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> sequence = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_contig) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> contig = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSEQ_alt_seq) {
      ptr -> alt_seq = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) INSDAltSeqDataAsnRead, (AsnOptFreeFunc) INSDAltSeqDataFree);
      if (isError && ptr -> alt_seq == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = INSDSeqFree(ptr);
   goto ret;
}



/**************************************************
*
*    INSDSeqAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
INSDSeqAsnWrite(INSDSeqPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, INSDSEQ);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> locus != NULL) {
      av.ptrvalue = ptr -> locus;
      retval = AsnWrite(aip, INSDSEQ_locus,  &av);
   }
   av.intvalue = ptr -> length;
   retval = AsnWrite(aip, INSDSEQ_length,  &av);
   if (ptr -> strandedness != NULL) {
      av.ptrvalue = ptr -> strandedness;
      retval = AsnWrite(aip, INSDSEQ_strandedness,  &av);
   }
   if (ptr -> moltype != NULL) {
      av.ptrvalue = ptr -> moltype;
      retval = AsnWrite(aip, INSDSEQ_moltype,  &av);
   }
   if (ptr -> topology != NULL) {
      av.ptrvalue = ptr -> topology;
      retval = AsnWrite(aip, INSDSEQ_topology,  &av);
   }
   if (ptr -> division != NULL) {
      av.ptrvalue = ptr -> division;
      retval = AsnWrite(aip, INSDSEQ_division,  &av);
   }
   if (ptr -> update_date != NULL) {
      av.ptrvalue = ptr -> update_date;
      retval = AsnWrite(aip, INSDSEQ_update_date,  &av);
   }
   if (ptr -> create_date != NULL) {
      av.ptrvalue = ptr -> create_date;
      retval = AsnWrite(aip, INSDSEQ_create_date,  &av);
   }
   if (ptr -> update_release != NULL) {
      av.ptrvalue = ptr -> update_release;
      retval = AsnWrite(aip, INSDSEQ_update_release,  &av);
   }
   if (ptr -> create_release != NULL) {
      av.ptrvalue = ptr -> create_release;
      retval = AsnWrite(aip, INSDSEQ_create_release,  &av);
   }
   if (ptr -> definition != NULL) {
      av.ptrvalue = ptr -> definition;
      retval = AsnWrite(aip, INSDSEQ_definition,  &av);
   }
   if (ptr -> primary_accession != NULL) {
      av.ptrvalue = ptr -> primary_accession;
      retval = AsnWrite(aip, INSDSEQ_primary_accession,  &av);
   }
   if (ptr -> entry_version != NULL) {
      av.ptrvalue = ptr -> entry_version;
      retval = AsnWrite(aip, INSDSEQ_entry_version,  &av);
   }
   if (ptr -> accession_version != NULL) {
      av.ptrvalue = ptr -> accession_version;
      retval = AsnWrite(aip, INSDSEQ_accession_version,  &av);
   }
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> other_seqids ,ASNCODE_PTRVAL_SLOT, aip, INSDSEQ_other_seqids, INSDSEQ_other_seqids_E);
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> secondary_accessions ,ASNCODE_PTRVAL_SLOT, aip, INSDSEQ_secondary_accessions, INSDSEQ_secondary_accessions_E);
   if (ptr -> project != NULL) {
      av.ptrvalue = ptr -> project;
      retval = AsnWrite(aip, INSDSEQ_project,  &av);
   }
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> keywords ,ASNCODE_PTRVAL_SLOT, aip, INSDSEQ_keywords, INSDSEQ_keywords_E);
   if (ptr -> segment != NULL) {
      av.ptrvalue = ptr -> segment;
      retval = AsnWrite(aip, INSDSEQ_segment,  &av);
   }
   if (ptr -> source != NULL) {
      av.ptrvalue = ptr -> source;
      retval = AsnWrite(aip, INSDSEQ_source,  &av);
   }
   if (ptr -> organism != NULL) {
      av.ptrvalue = ptr -> organism;
      retval = AsnWrite(aip, INSDSEQ_organism,  &av);
   }
   if (ptr -> taxonomy != NULL) {
      av.ptrvalue = ptr -> taxonomy;
      retval = AsnWrite(aip, INSDSEQ_taxonomy,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> references, (AsnWriteFunc) INSDReferenceAsnWrite, aip, INSDSEQ_references, INSDSEQ_references_E);
   if (ptr -> comment != NULL) {
      av.ptrvalue = ptr -> comment;
      retval = AsnWrite(aip, INSDSEQ_comment,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> comment_set, (AsnWriteFunc) INSDCommentAsnWrite, aip, INSDSEQ_comment_set, INSDSEQ_comment_set_E);
   AsnGenericUserSeqOfAsnWrite(ptr -> struc_comments, (AsnWriteFunc) INSDStrucCommentAsnWrite, aip, INSDSEQ_struc_comments, INSDSEQ_struc_comments_E);
   if (ptr -> primary != NULL) {
      av.ptrvalue = ptr -> primary;
      retval = AsnWrite(aip, INSDSEQ_primary,  &av);
   }
   if (ptr -> source_db != NULL) {
      av.ptrvalue = ptr -> source_db;
      retval = AsnWrite(aip, INSDSEQ_source_db,  &av);
   }
   if (ptr -> database_reference != NULL) {
      av.ptrvalue = ptr -> database_reference;
      retval = AsnWrite(aip, INSDSEQ_database_reference,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> feature_table, (AsnWriteFunc) INSDFeatureAsnWrite, aip, INSDSEQ_feature_table, INSDSEQ_feature_table_E);
   AsnGenericUserSeqOfAsnWrite(ptr -> feature_set, (AsnWriteFunc) INSDFeatureSetAsnWrite, aip, INSDSEQ_feature_set, INSDSEQ_feature_set_E);
   if (ptr -> sequence != NULL) {
      av.ptrvalue = ptr -> sequence;
      retval = AsnWrite(aip, INSDSEQ_sequence,  &av);
   }
   if (ptr -> contig != NULL) {
      av.ptrvalue = ptr -> contig;
      retval = AsnWrite(aip, INSDSEQ_contig,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> alt_seq, (AsnWriteFunc) INSDAltSeqDataAsnWrite, aip, INSDSEQ_alt_seq, INSDSEQ_alt_seq_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    INSDReferenceNew()
*
**************************************************/
NLM_EXTERN 
INSDReferencePtr LIBCALL
INSDReferenceNew(void)
{
   INSDReferencePtr ptr = MemNew((size_t) sizeof(INSDReference));

   return ptr;

}


/**************************************************
*
*    INSDReferenceFree()
*
**************************************************/
NLM_EXTERN 
INSDReferencePtr LIBCALL
INSDReferenceFree(INSDReferencePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> reference);
   MemFree(ptr -> position);
   AsnGenericBaseSeqOfFree(ptr -> authors ,ASNCODE_PTRVAL_SLOT);
   MemFree(ptr -> consortium);
   MemFree(ptr -> title);
   MemFree(ptr -> journal);
   AsnGenericUserSeqOfFree(ptr -> xref, (AsnOptFreeFunc) INSDXrefFree);
   MemFree(ptr -> remark);
   return MemFree(ptr);
}


/**************************************************
*
*    INSDReferenceAsnRead()
*
**************************************************/
NLM_EXTERN 
INSDReferencePtr LIBCALL
INSDReferenceAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   INSDReferencePtr ptr;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* INSDReference ::= (self contained) */
      atp = AsnReadId(aip, amp, INSDREFERENCE);
   } else {
      atp = AsnLinkType(orig, INSDREFERENCE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = INSDReferenceNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == INSDREFERENCE_reference) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> reference = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDREFERENCE_position) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> position = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDREFERENCE_authors) {
      ptr -> authors = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> authors == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDREFERENCE_consortium) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> consortium = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDREFERENCE_title) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> title = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDREFERENCE_journal) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> journal = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDREFERENCE_xref) {
      ptr -> xref = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) INSDXrefAsnRead, (AsnOptFreeFunc) INSDXrefFree);
      if (isError && ptr -> xref == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDREFERENCE_pubmed) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> pubmed = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDREFERENCE_remark) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> remark = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = INSDReferenceFree(ptr);
   goto ret;
}



/**************************************************
*
*    INSDReferenceAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
INSDReferenceAsnWrite(INSDReferencePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, INSDREFERENCE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> reference != NULL) {
      av.ptrvalue = ptr -> reference;
      retval = AsnWrite(aip, INSDREFERENCE_reference,  &av);
   }
   if (ptr -> position != NULL) {
      av.ptrvalue = ptr -> position;
      retval = AsnWrite(aip, INSDREFERENCE_position,  &av);
   }
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> authors ,ASNCODE_PTRVAL_SLOT, aip, INSDREFERENCE_authors, INSDREFERENCE_authors_E);
   if (ptr -> consortium != NULL) {
      av.ptrvalue = ptr -> consortium;
      retval = AsnWrite(aip, INSDREFERENCE_consortium,  &av);
   }
   if (ptr -> title != NULL) {
      av.ptrvalue = ptr -> title;
      retval = AsnWrite(aip, INSDREFERENCE_title,  &av);
   }
   if (ptr -> journal != NULL) {
      av.ptrvalue = ptr -> journal;
      retval = AsnWrite(aip, INSDREFERENCE_journal,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> xref, (AsnWriteFunc) INSDXrefAsnWrite, aip, INSDREFERENCE_xref, INSDREFERENCE_xref_E);
   if (ptr -> pubmed || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> pubmed;
      retval = AsnWrite(aip, INSDREFERENCE_pubmed,  &av);
   }
   if (ptr -> remark != NULL) {
      av.ptrvalue = ptr -> remark;
      retval = AsnWrite(aip, INSDREFERENCE_remark,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    INSDCommentNew()
*
**************************************************/
NLM_EXTERN 
INSDCommentPtr LIBCALL
INSDCommentNew(void)
{
   INSDCommentPtr ptr = MemNew((size_t) sizeof(INSDComment));

   return ptr;

}


/**************************************************
*
*    INSDCommentFree()
*
**************************************************/
NLM_EXTERN 
INSDCommentPtr LIBCALL
INSDCommentFree(INSDCommentPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> type);
   AsnGenericUserSeqOfFree(ptr -> paragraphs, (AsnOptFreeFunc) INSDCommentParagraphFree);
   return MemFree(ptr);
}


/**************************************************
*
*    INSDCommentAsnRead()
*
**************************************************/
NLM_EXTERN 
INSDCommentPtr LIBCALL
INSDCommentAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   INSDCommentPtr ptr;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* INSDComment ::= (self contained) */
      atp = AsnReadId(aip, amp, INSDCOMMENT);
   } else {
      atp = AsnLinkType(orig, INSDCOMMENT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = INSDCommentNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == INSDCOMMENT_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDCOMMENT_paragraphs) {
      ptr -> paragraphs = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) INSDCommentParagraphAsnRead, (AsnOptFreeFunc) INSDCommentParagraphFree);
      if (isError && ptr -> paragraphs == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = INSDCommentFree(ptr);
   goto ret;
}



/**************************************************
*
*    INSDCommentAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
INSDCommentAsnWrite(INSDCommentPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, INSDCOMMENT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> type != NULL) {
      av.ptrvalue = ptr -> type;
      retval = AsnWrite(aip, INSDCOMMENT_type,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> paragraphs, (AsnWriteFunc) INSDCommentParagraphAsnWrite, aip, INSDCOMMENT_paragraphs, INSDCOMMENT_paragraphs_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    INSDStrucCommentNew()
*
**************************************************/
NLM_EXTERN 
INSDStrucCommentPtr LIBCALL
INSDStrucCommentNew(void)
{
   INSDStrucCommentPtr ptr = MemNew((size_t) sizeof(INSDStrucComment));

   return ptr;

}


/**************************************************
*
*    INSDStrucCommentFree()
*
**************************************************/
NLM_EXTERN 
INSDStrucCommentPtr LIBCALL
INSDStrucCommentFree(INSDStrucCommentPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> name);
   AsnGenericUserSeqOfFree(ptr -> items, (AsnOptFreeFunc) INSDStrucCommentItemFree);
   return MemFree(ptr);
}


/**************************************************
*
*    INSDStrucCommentAsnRead()
*
**************************************************/
NLM_EXTERN 
INSDStrucCommentPtr LIBCALL
INSDStrucCommentAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   INSDStrucCommentPtr ptr;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* INSDStrucComment ::= (self contained) */
      atp = AsnReadId(aip, amp, INSDSTRUCCOMMENT);
   } else {
      atp = AsnLinkType(orig, INSDSTRUCCOMMENT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = INSDStrucCommentNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == INSDSTRUCCOMMENT_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSTRUCCOMMENT_items) {
      ptr -> items = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) INSDStrucCommentItemAsnRead, (AsnOptFreeFunc) INSDStrucCommentItemFree);
      if (isError && ptr -> items == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = INSDStrucCommentFree(ptr);
   goto ret;
}



/**************************************************
*
*    INSDStrucCommentAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
INSDStrucCommentAsnWrite(INSDStrucCommentPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, INSDSTRUCCOMMENT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> name != NULL) {
      av.ptrvalue = ptr -> name;
      retval = AsnWrite(aip, INSDSTRUCCOMMENT_name,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> items, (AsnWriteFunc) INSDStrucCommentItemAsnWrite, aip, INSDSTRUCCOMMENT_items, INSDSTRUCCOMMENT_items_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    INSDFeatureNew()
*
**************************************************/
NLM_EXTERN 
INSDFeaturePtr LIBCALL
INSDFeatureNew(void)
{
   INSDFeaturePtr ptr = MemNew((size_t) sizeof(INSDFeature));

   return ptr;

}


/**************************************************
*
*    INSDFeatureFree()
*
**************************************************/
NLM_EXTERN 
INSDFeaturePtr LIBCALL
INSDFeatureFree(INSDFeaturePtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> key);
   MemFree(ptr -> location);
   AsnGenericUserSeqOfFree(ptr -> intervals, (AsnOptFreeFunc) INSDIntervalFree);
   MemFree(ptr -> operator__);
   AsnGenericUserSeqOfFree(ptr -> quals, (AsnOptFreeFunc) INSDQualifierFree);
   AsnGenericUserSeqOfFree(ptr -> xrefs, (AsnOptFreeFunc) INSDXrefFree);
   return MemFree(ptr);
}


/**************************************************
*
*    INSDFeatureAsnRead()
*
**************************************************/
NLM_EXTERN 
INSDFeaturePtr LIBCALL
INSDFeatureAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   INSDFeaturePtr ptr;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* INSDFeature ::= (self contained) */
      atp = AsnReadId(aip, amp, INSDFEATURE);
   } else {
      atp = AsnLinkType(orig, INSDFEATURE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = INSDFeatureNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == INSDFEATURE_key) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> key = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDFEATURE_location) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> location = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDFEATURE_intervals) {
      ptr -> intervals = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) INSDIntervalAsnRead, (AsnOptFreeFunc) INSDIntervalFree);
      if (isError && ptr -> intervals == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDFEATURE_operator) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> operator__ = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDFEATURE_partial5) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> partial5 = av.boolvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDFEATURE_partial3) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> partial3 = av.boolvalue;
      ptr -> OBbits__ |= 1<<1;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDFEATURE_quals) {
      ptr -> quals = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) INSDQualifierAsnRead, (AsnOptFreeFunc) INSDQualifierFree);
      if (isError && ptr -> quals == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDFEATURE_xrefs) {
      ptr -> xrefs = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) INSDXrefAsnRead, (AsnOptFreeFunc) INSDXrefFree);
      if (isError && ptr -> xrefs == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = INSDFeatureFree(ptr);
   goto ret;
}



/**************************************************
*
*    INSDFeatureAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
INSDFeatureAsnWrite(INSDFeaturePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, INSDFEATURE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> key != NULL) {
      av.ptrvalue = ptr -> key;
      retval = AsnWrite(aip, INSDFEATURE_key,  &av);
   }
   if (ptr -> location != NULL) {
      av.ptrvalue = ptr -> location;
      retval = AsnWrite(aip, INSDFEATURE_location,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> intervals, (AsnWriteFunc) INSDIntervalAsnWrite, aip, INSDFEATURE_intervals, INSDFEATURE_intervals_E);
   if (ptr -> operator__ != NULL) {
      av.ptrvalue = ptr -> operator__;
      retval = AsnWrite(aip, INSDFEATURE_operator,  &av);
   }
   if (ptr -> partial5 || (ptr -> OBbits__ & (1<<0) )){   av.boolvalue = ptr -> partial5;
      retval = AsnWrite(aip, INSDFEATURE_partial5,  &av);
   }
   if (ptr -> partial3 || (ptr -> OBbits__ & (1<<1) )){   av.boolvalue = ptr -> partial3;
      retval = AsnWrite(aip, INSDFEATURE_partial3,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> quals, (AsnWriteFunc) INSDQualifierAsnWrite, aip, INSDFEATURE_quals, INSDFEATURE_quals_E);
   AsnGenericUserSeqOfAsnWrite(ptr -> xrefs, (AsnWriteFunc) INSDXrefAsnWrite, aip, INSDFEATURE_xrefs, INSDFEATURE_xrefs_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    INSDFeatureSetNew()
*
**************************************************/
NLM_EXTERN 
INSDFeatureSetPtr LIBCALL
INSDFeatureSetNew(void)
{
   INSDFeatureSetPtr ptr = MemNew((size_t) sizeof(INSDFeatureSet));

   return ptr;

}


/**************************************************
*
*    INSDFeatureSetFree()
*
**************************************************/
NLM_EXTERN 
INSDFeatureSetPtr LIBCALL
INSDFeatureSetFree(INSDFeatureSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> annot_source);
   AsnGenericUserSeqOfFree(ptr -> features, (AsnOptFreeFunc) INSDFeatureFree);
   return MemFree(ptr);
}


/**************************************************
*
*    INSDFeatureSetAsnRead()
*
**************************************************/
NLM_EXTERN 
INSDFeatureSetPtr LIBCALL
INSDFeatureSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   INSDFeatureSetPtr ptr;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* INSDFeatureSet ::= (self contained) */
      atp = AsnReadId(aip, amp, INSDFEATURESET);
   } else {
      atp = AsnLinkType(orig, INSDFEATURESET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = INSDFeatureSetNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == INSDFEATURESET_annot_source) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> annot_source = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDFEATURESET_features) {
      ptr -> features = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) INSDFeatureAsnRead, (AsnOptFreeFunc) INSDFeatureFree);
      if (isError && ptr -> features == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = INSDFeatureSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    INSDFeatureSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
INSDFeatureSetAsnWrite(INSDFeatureSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, INSDFEATURESET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> annot_source != NULL) {
      av.ptrvalue = ptr -> annot_source;
      retval = AsnWrite(aip, INSDFEATURESET_annot_source,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> features, (AsnWriteFunc) INSDFeatureAsnWrite, aip, INSDFEATURESET_features, INSDFEATURESET_features_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    INSDAltSeqDataNew()
*
**************************************************/
NLM_EXTERN 
INSDAltSeqDataPtr LIBCALL
INSDAltSeqDataNew(void)
{
   INSDAltSeqDataPtr ptr = MemNew((size_t) sizeof(INSDAltSeqData));

   return ptr;

}


/**************************************************
*
*    INSDAltSeqDataFree()
*
**************************************************/
NLM_EXTERN 
INSDAltSeqDataPtr LIBCALL
INSDAltSeqDataFree(INSDAltSeqDataPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> name);
   AsnGenericUserSeqOfFree(ptr -> items, (AsnOptFreeFunc) INSDAltSeqItemFree);
   return MemFree(ptr);
}


/**************************************************
*
*    INSDAltSeqDataAsnRead()
*
**************************************************/
NLM_EXTERN 
INSDAltSeqDataPtr LIBCALL
INSDAltSeqDataAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   INSDAltSeqDataPtr ptr;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* INSDAltSeqData ::= (self contained) */
      atp = AsnReadId(aip, amp, INSDALTSEQDATA);
   } else {
      atp = AsnLinkType(orig, INSDALTSEQDATA);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = INSDAltSeqDataNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == INSDALTSEQDATA_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDALTSEQDATA_items) {
      ptr -> items = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) INSDAltSeqItemAsnRead, (AsnOptFreeFunc) INSDAltSeqItemFree);
      if (isError && ptr -> items == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = INSDAltSeqDataFree(ptr);
   goto ret;
}



/**************************************************
*
*    INSDAltSeqDataAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
INSDAltSeqDataAsnWrite(INSDAltSeqDataPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, INSDALTSEQDATA);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> name != NULL) {
      av.ptrvalue = ptr -> name;
      retval = AsnWrite(aip, INSDALTSEQDATA_name,  &av);
   }
   AsnGenericUserSeqOfAsnWrite(ptr -> items, (AsnWriteFunc) INSDAltSeqItemAsnWrite, aip, INSDALTSEQDATA_items, INSDALTSEQDATA_items_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    INSDXrefNew()
*
**************************************************/
NLM_EXTERN 
INSDXrefPtr LIBCALL
INSDXrefNew(void)
{
   INSDXrefPtr ptr = MemNew((size_t) sizeof(INSDXref));

   return ptr;

}


/**************************************************
*
*    INSDXrefFree()
*
**************************************************/
NLM_EXTERN 
INSDXrefPtr LIBCALL
INSDXrefFree(INSDXrefPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> dbname);
   MemFree(ptr -> id);
   return MemFree(ptr);
}


/**************************************************
*
*    INSDXrefAsnRead()
*
**************************************************/
NLM_EXTERN 
INSDXrefPtr LIBCALL
INSDXrefAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   INSDXrefPtr ptr;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* INSDXref ::= (self contained) */
      atp = AsnReadId(aip, amp, INSDXREF);
   } else {
      atp = AsnLinkType(orig, INSDXREF);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = INSDXrefNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == INSDXREF_dbname) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> dbname = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDXREF_id) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> id = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = INSDXrefFree(ptr);
   goto ret;
}



/**************************************************
*
*    INSDXrefAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
INSDXrefAsnWrite(INSDXrefPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, INSDXREF);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> dbname != NULL) {
      av.ptrvalue = ptr -> dbname;
      retval = AsnWrite(aip, INSDXREF_dbname,  &av);
   }
   if (ptr -> id != NULL) {
      av.ptrvalue = ptr -> id;
      retval = AsnWrite(aip, INSDXREF_id,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    INSDCommentParagraphNew()
*
**************************************************/
NLM_EXTERN 
INSDCommentParagraphPtr LIBCALL
INSDCommentParagraphNew(void)
{
   INSDCommentParagraphPtr ptr = MemNew((size_t) sizeof(INSDCommentParagraph));

   return ptr;

}


/**************************************************
*
*    INSDCommentParagraphFree()
*
**************************************************/
NLM_EXTERN 
INSDCommentParagraphPtr LIBCALL
INSDCommentParagraphFree(INSDCommentParagraphPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr -> items, (AsnOptFreeFunc) INSDCommentItemFree);
   return MemFree(ptr);
}


/**************************************************
*
*    INSDCommentParagraphAsnRead()
*
**************************************************/
NLM_EXTERN 
INSDCommentParagraphPtr LIBCALL
INSDCommentParagraphAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   INSDCommentParagraphPtr ptr;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* INSDCommentParagraph ::= (self contained) */
      atp = AsnReadId(aip, amp, INSDCOMMENTPARAGRAPH);
   } else {
      atp = AsnLinkType(orig, INSDCOMMENTPARAGRAPH);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = INSDCommentParagraphNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == INSDCOMMENTPARAGRAPH_items) {
      ptr -> items = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) INSDCommentItemAsnRead, (AsnOptFreeFunc) INSDCommentItemFree);
      if (isError && ptr -> items == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = INSDCommentParagraphFree(ptr);
   goto ret;
}



/**************************************************
*
*    INSDCommentParagraphAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
INSDCommentParagraphAsnWrite(INSDCommentParagraphPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, INSDCOMMENTPARAGRAPH);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   AsnGenericUserSeqOfAsnWrite(ptr -> items, (AsnWriteFunc) INSDCommentItemAsnWrite, aip, INSDCOMMENTPARAGRAPH_items, INSDCOMMENTPARAGRAPH_items_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    INSDCommentItemNew()
*
**************************************************/
NLM_EXTERN 
INSDCommentItemPtr LIBCALL
INSDCommentItemNew(void)
{
   INSDCommentItemPtr ptr = MemNew((size_t) sizeof(INSDCommentItem));

   return ptr;

}


/**************************************************
*
*    INSDCommentItemFree()
*
**************************************************/
NLM_EXTERN 
INSDCommentItemPtr LIBCALL
INSDCommentItemFree(INSDCommentItemPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> value);
   MemFree(ptr -> url);
   return MemFree(ptr);
}


/**************************************************
*
*    INSDCommentItemAsnRead()
*
**************************************************/
NLM_EXTERN 
INSDCommentItemPtr LIBCALL
INSDCommentItemAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   INSDCommentItemPtr ptr;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* INSDCommentItem ::= (self contained) */
      atp = AsnReadId(aip, amp, INSDCOMMENTITEM);
   } else {
      atp = AsnLinkType(orig, INSDCOMMENTITEM);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = INSDCommentItemNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == INSDCOMMENTITEM_value) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> value = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDCOMMENTITEM_url) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> url = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = INSDCommentItemFree(ptr);
   goto ret;
}



/**************************************************
*
*    INSDCommentItemAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
INSDCommentItemAsnWrite(INSDCommentItemPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, INSDCOMMENTITEM);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> value != NULL) {
      av.ptrvalue = ptr -> value;
      retval = AsnWrite(aip, INSDCOMMENTITEM_value,  &av);
   }
   if (ptr -> url != NULL) {
      av.ptrvalue = ptr -> url;
      retval = AsnWrite(aip, INSDCOMMENTITEM_url,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    INSDStrucCommentItemNew()
*
**************************************************/
NLM_EXTERN 
INSDStrucCommentItemPtr LIBCALL
INSDStrucCommentItemNew(void)
{
   INSDStrucCommentItemPtr ptr = MemNew((size_t) sizeof(INSDStrucCommentItem));

   return ptr;

}


/**************************************************
*
*    INSDStrucCommentItemFree()
*
**************************************************/
NLM_EXTERN 
INSDStrucCommentItemPtr LIBCALL
INSDStrucCommentItemFree(INSDStrucCommentItemPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> tag);
   MemFree(ptr -> value);
   MemFree(ptr -> url);
   return MemFree(ptr);
}


/**************************************************
*
*    INSDStrucCommentItemAsnRead()
*
**************************************************/
NLM_EXTERN 
INSDStrucCommentItemPtr LIBCALL
INSDStrucCommentItemAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   INSDStrucCommentItemPtr ptr;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* INSDStrucCommentItem ::= (self contained) */
      atp = AsnReadId(aip, amp, INSDSTRUCCOMMENTITEM);
   } else {
      atp = AsnLinkType(orig, INSDSTRUCCOMMENTITEM);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = INSDStrucCommentItemNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == INSDSTRUCCOMMENTITEM_tag) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> tag = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSTRUCCOMMENTITEM_value) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> value = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDSTRUCCOMMENTITEM_url) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> url = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = INSDStrucCommentItemFree(ptr);
   goto ret;
}



/**************************************************
*
*    INSDStrucCommentItemAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
INSDStrucCommentItemAsnWrite(INSDStrucCommentItemPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, INSDSTRUCCOMMENTITEM);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> tag != NULL) {
      av.ptrvalue = ptr -> tag;
      retval = AsnWrite(aip, INSDSTRUCCOMMENTITEM_tag,  &av);
   }
   if (ptr -> value != NULL) {
      av.ptrvalue = ptr -> value;
      retval = AsnWrite(aip, INSDSTRUCCOMMENTITEM_value,  &av);
   }
   if (ptr -> url != NULL) {
      av.ptrvalue = ptr -> url;
      retval = AsnWrite(aip, INSDSTRUCCOMMENTITEM_url,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    INSDIntervalNew()
*
**************************************************/
NLM_EXTERN 
INSDIntervalPtr LIBCALL
INSDIntervalNew(void)
{
   INSDIntervalPtr ptr = MemNew((size_t) sizeof(INSDInterval));

   return ptr;

}


/**************************************************
*
*    INSDIntervalFree()
*
**************************************************/
NLM_EXTERN 
INSDIntervalPtr LIBCALL
INSDIntervalFree(INSDIntervalPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> accession);
   return MemFree(ptr);
}


/**************************************************
*
*    INSDIntervalAsnRead()
*
**************************************************/
NLM_EXTERN 
INSDIntervalPtr LIBCALL
INSDIntervalAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   INSDIntervalPtr ptr;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* INSDInterval ::= (self contained) */
      atp = AsnReadId(aip, amp, INSDINTERVAL);
   } else {
      atp = AsnLinkType(orig, INSDINTERVAL);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = INSDIntervalNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == INSDINTERVAL_from) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> from = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDINTERVAL_to) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> to = av.intvalue;
      ptr -> OBbits__ |= 1<<1;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDINTERVAL_point) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> point = av.intvalue;
      ptr -> OBbits__ |= 1<<2;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDINTERVAL_iscomp) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> iscomp = av.boolvalue;
      ptr -> OBbits__ |= 1<<3;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDINTERVAL_interbp) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> interbp = av.boolvalue;
      ptr -> OBbits__ |= 1<<4;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDINTERVAL_accession) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> accession = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = INSDIntervalFree(ptr);
   goto ret;
}



/**************************************************
*
*    INSDIntervalAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
INSDIntervalAsnWrite(INSDIntervalPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, INSDINTERVAL);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> from || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> from;
      retval = AsnWrite(aip, INSDINTERVAL_from,  &av);
   }
   if (ptr -> to || (ptr -> OBbits__ & (1<<1) )){   av.intvalue = ptr -> to;
      retval = AsnWrite(aip, INSDINTERVAL_to,  &av);
   }
   if (ptr -> point || (ptr -> OBbits__ & (1<<2) )){   av.intvalue = ptr -> point;
      retval = AsnWrite(aip, INSDINTERVAL_point,  &av);
   }
   if (ptr -> iscomp || (ptr -> OBbits__ & (1<<3) )){   av.boolvalue = ptr -> iscomp;
      retval = AsnWrite(aip, INSDINTERVAL_iscomp,  &av);
   }
   if (ptr -> interbp || (ptr -> OBbits__ & (1<<4) )){   av.boolvalue = ptr -> interbp;
      retval = AsnWrite(aip, INSDINTERVAL_interbp,  &av);
   }
   if (ptr -> accession != NULL) {
      av.ptrvalue = ptr -> accession;
      retval = AsnWrite(aip, INSDINTERVAL_accession,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    INSDQualifierNew()
*
**************************************************/
NLM_EXTERN 
INSDQualifierPtr LIBCALL
INSDQualifierNew(void)
{
   INSDQualifierPtr ptr = MemNew((size_t) sizeof(INSDQualifier));

   return ptr;

}


/**************************************************
*
*    INSDQualifierFree()
*
**************************************************/
NLM_EXTERN 
INSDQualifierPtr LIBCALL
INSDQualifierFree(INSDQualifierPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> name);
   MemFree(ptr -> value);
   return MemFree(ptr);
}


/**************************************************
*
*    INSDQualifierAsnRead()
*
**************************************************/
NLM_EXTERN 
INSDQualifierPtr LIBCALL
INSDQualifierAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   INSDQualifierPtr ptr;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* INSDQualifier ::= (self contained) */
      atp = AsnReadId(aip, amp, INSDQUALIFIER);
   } else {
      atp = AsnLinkType(orig, INSDQUALIFIER);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = INSDQualifierNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == INSDQUALIFIER_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDQUALIFIER_value) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> value = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = INSDQualifierFree(ptr);
   goto ret;
}



/**************************************************
*
*    INSDQualifierAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
INSDQualifierAsnWrite(INSDQualifierPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, INSDQUALIFIER);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> name != NULL) {
      av.ptrvalue = ptr -> name;
      retval = AsnWrite(aip, INSDQUALIFIER_name,  &av);
   }
   if (ptr -> value != NULL) {
      av.ptrvalue = ptr -> value;
      retval = AsnWrite(aip, INSDQUALIFIER_value,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    INSDAltSeqItemNew()
*
**************************************************/
NLM_EXTERN 
INSDAltSeqItemPtr LIBCALL
INSDAltSeqItemNew(void)
{
   INSDAltSeqItemPtr ptr = MemNew((size_t) sizeof(INSDAltSeqItem));

   return ptr;

}


/**************************************************
*
*    INSDAltSeqItemFree()
*
**************************************************/
NLM_EXTERN 
INSDAltSeqItemPtr LIBCALL
INSDAltSeqItemFree(INSDAltSeqItemPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   INSDIntervalFree(ptr -> interval);
   MemFree(ptr -> gap_type);
   MemFree(ptr -> gap_linkage);
   MemFree(ptr -> gap_comment);
   MemFree(ptr -> first_accn);
   MemFree(ptr -> last_accn);
   MemFree(ptr -> value);
   return MemFree(ptr);
}


/**************************************************
*
*    INSDAltSeqItemAsnRead()
*
**************************************************/
NLM_EXTERN 
INSDAltSeqItemPtr LIBCALL
INSDAltSeqItemAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   INSDAltSeqItemPtr ptr;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* INSDAltSeqItem ::= (self contained) */
      atp = AsnReadId(aip, amp, INSDALTSEQITEM);
   } else {
      atp = AsnLinkType(orig, INSDALTSEQITEM);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = INSDAltSeqItemNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == INSDALTSEQITEM_interval) {
      ptr -> interval = INSDIntervalAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDALTSEQITEM_isgap) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> isgap = av.boolvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDALTSEQITEM_gap_length) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> gap_length = av.intvalue;
      ptr -> OBbits__ |= 1<<1;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDALTSEQITEM_gap_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> gap_type = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDALTSEQITEM_gap_linkage) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> gap_linkage = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDALTSEQITEM_gap_comment) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> gap_comment = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDALTSEQITEM_first_accn) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> first_accn = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDALTSEQITEM_last_accn) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> last_accn = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INSDALTSEQITEM_value) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> value = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = INSDAltSeqItemFree(ptr);
   goto ret;
}



/**************************************************
*
*    INSDAltSeqItemAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
INSDAltSeqItemAsnWrite(INSDAltSeqItemPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! objinsdseqAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, INSDALTSEQITEM);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> interval != NULL) {
      if ( ! INSDIntervalAsnWrite(ptr -> interval, aip, INSDALTSEQITEM_interval)) {
         goto erret;
      }
   }
   if (ptr -> isgap || (ptr -> OBbits__ & (1<<0) )){   av.boolvalue = ptr -> isgap;
      retval = AsnWrite(aip, INSDALTSEQITEM_isgap,  &av);
   }
   if (ptr -> gap_length || (ptr -> OBbits__ & (1<<1) )){   av.intvalue = ptr -> gap_length;
      retval = AsnWrite(aip, INSDALTSEQITEM_gap_length,  &av);
   }
   if (ptr -> gap_type != NULL) {
      av.ptrvalue = ptr -> gap_type;
      retval = AsnWrite(aip, INSDALTSEQITEM_gap_type,  &av);
   }
   if (ptr -> gap_linkage != NULL) {
      av.ptrvalue = ptr -> gap_linkage;
      retval = AsnWrite(aip, INSDALTSEQITEM_gap_linkage,  &av);
   }
   if (ptr -> gap_comment != NULL) {
      av.ptrvalue = ptr -> gap_comment;
      retval = AsnWrite(aip, INSDALTSEQITEM_gap_comment,  &av);
   }
   if (ptr -> first_accn != NULL) {
      av.ptrvalue = ptr -> first_accn;
      retval = AsnWrite(aip, INSDALTSEQITEM_first_accn,  &av);
   }
   if (ptr -> last_accn != NULL) {
      av.ptrvalue = ptr -> last_accn;
      retval = AsnWrite(aip, INSDALTSEQITEM_last_accn,  &av);
   }
   if (ptr -> value != NULL) {
      av.ptrvalue = ptr -> value;
      retval = AsnWrite(aip, INSDALTSEQITEM_value,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}

