/*  objseq.h
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* File Name:  objseq.h
*
* Author:  James Ostell
*   
* Version Creation Date: 4/1/91
*
* $Revision: 6.12 $
*
* File Description:  Object manager interface for module NCBI-Seq
*
* Modifications:  
* --------------------------------------------------------------------------
* Date	   Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*/

#ifndef _NCBI_Seq_
#define _NCBI_Seq_

#ifndef _ASNTOOL_
#include <asn.h>
#endif
#ifndef _NCBI_General_
#include <objgen.h>
#endif
#ifndef _NCBI_Seqloc_
#include <objloc.h>
#endif
#ifndef _NCBI_Pub_
#include <objpub.h>
#endif
#ifndef _NCBI_Seqalign_
#include <objalign.h>
#endif
#ifndef _NCBI_Pubdesc_
#include <objpubd.h>     /* separated out to avoid typedef order problems */
#endif
#ifndef _NCBI_Seqfeat_
#include <objfeat.h>       /* include organism for now */
#endif
#ifndef _NCBI_Seqres_
#include <objres.h>
#endif
#ifndef _NCBI_SeqBlock_
#include <objblock.h>
#endif
#ifndef _NCBI_SeqCode_
#include <objcode.h>
#endif

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef ValNode SeqEntry, FAR *SeqEntryPtr;

/*****************************************************************************
*
*   loader
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqAsnLoad PROTO((void));

/*****************************************************************************
*
*   internal structures for NCBI-Seq objects
*
*****************************************************************************/

/*****************************************************************************
*
*   SeqAnnot - Sequence annotations
*
*****************************************************************************/

/*****************************************************************************
*
*   TextAnnotId
*
*****************************************************************************/
typedef struct textannotid {
    CharPtr name,
        accession,
        release;
	Int2 version;             /* INT2_MIN (ncbilcl.h) = not set */
} TextAnnotId, PNTR TextAnnotIdPtr;

NLM_EXTERN TextAnnotIdPtr LIBCALL TextAnnotIdNew PROTO((void));
NLM_EXTERN Boolean      LIBCALL TextAnnotIdAsnWrite PROTO((TextAnnotIdPtr tsip, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN TextAnnotIdPtr LIBCALL TextAnnotIdAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN TextAnnotIdPtr LIBCALL TextAnnotIdFree PROTO((TextAnnotIdPtr tsip));

/*****************************************************************************
*
*   AnnotId is a choice using an ValNode, most types in data.ptrvalue
*      except integers, in data.intvalue
*   choice:
*   0 = not set
    1 = local Object-id ,      -- local use
	2 = ncbi INTEGER ,
	3 = general Dbtag ,        -- any database
	4 = other Textannot-id
*
*****************************************************************************/

typedef ValNode AnnotId, PNTR AnnotIdPtr;

#define ANNOTID_NOT_SET ( (Uint1)0)
#define ANNOTID_LOCAL ( (Uint1)1)
#define ANNOTID_NCBI ( (Uint1)2)
#define ANNOTID_GENERAL ((Uint1)3)
#define ANNOTID_OTHER ((Uint1)4)

NLM_EXTERN Boolean	 LIBCALL AnnotIdAsnWrite PROTO((AnnotIdPtr anp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN AnnotIdPtr LIBCALL AnnotIdAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN AnnotIdPtr LIBCALL AnnotIdFree PROTO((AnnotIdPtr anp));

/*****************************************************************************
*
*   These routines process sets or sequences of AnnotId's
*
*****************************************************************************/
NLM_EXTERN Boolean  LIBCALL AnnotIdSetAsnWrite PROTO((AnnotIdPtr anp, AsnIoPtr aip, AsnTypePtr settype, AsnTypePtr elementtype));
NLM_EXTERN AnnotIdPtr LIBCALL AnnotIdSetAsnRead PROTO((AsnIoPtr aip, AsnTypePtr settype, AsnTypePtr elementtype));
NLM_EXTERN AnnotIdPtr LIBCALL AnnotIdSetFree PROTO((AnnotIdPtr anp));

/*****************************************************************************
*
*   AlignDef
*     used to describe a set of alignments in a SeqAnnot
*     align_type
*        ALIGN_TYPE_REF   set of pairwise alignments to one "master" sequence
*        ALIGN_TYPE_ALT   set of alternate alignments of the same seqs
*        ALIGN_TYPE_BLOCKS  set of (unordered) aligned blocks
*    ids    SeqIds used in alignments
*        currently only used for ALIGN_TYPE_REF, where the one SeqId of
*           the "master" sequence is given (one SeqId)
*
*****************************************************************************/
typedef struct aligndef {
	Uint1 align_type;
	SeqIdPtr ids;
} AlignDef, PNTR AlignDefPtr;

#define ALIGN_TYPE_REF  1
#define ALIGN_TYPE_ALT  2
#define ALIGN_TYPE_BLOCKS 3
#define ALIGN_TYPE_OTHER 255

NLM_EXTERN AlignDefPtr LIBCALL AlignDefNew PROTO((void));
NLM_EXTERN Boolean     LIBCALL AlignDefAsnWrite PROTO((AlignDefPtr adp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN AlignDefPtr LIBCALL AlignDefAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN AlignDefPtr LIBCALL AlignDefFree PROTO((AlignDefPtr adp));

/*****************************************************************************
*
*   AnnotDescr uses an ValNode with choice = 
    1 = name VisibleString ,         -- a name for this annotation
    2 = title VisibleString ,        -- a title for this annotation
    3 = comment VisibleString ,      -- a more extensive comment
    4 = pub Pubdesc                 -- a reference to the publication
    5 = user UserObject             -- user structured data object
	6 = create-date Date            -- date entry created
	7 = update-date Date            -- date of last update
	8 = src Seq-id                  -- SeqId of entry this came from
	9 = align Align-def             -- Align-def pointer
	10= region Seq-loc              -- Seq-loc of region covered by this annot

    types with * use data.intvalue.  Other use data.ptrvalue
    ** uses a chain of ValNodes which use data.intvalue for enumerated type
*
*****************************************************************************/
#define Annot_descr_name 1
#define Annot_descr_title 2
#define Annot_descr_comment 3
#define Annot_descr_pub 4
#define Annot_descr_user 5
#define Annot_descr_create_date 6
#define Annot_descr_update_date 7
#define Annot_descr_src 8
#define Annot_descr_align 9
#define Annot_descr_region 10

typedef ValNode AnnotDescr, PNTR AnnotDescrPtr;

#define AnnotDescrNew SeqDescrNew
NLM_EXTERN Boolean    LIBCALL AnnotDescrAsnWrite PROTO((AnnotDescrPtr anp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN AnnotDescrPtr LIBCALL AnnotDescrAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN AnnotDescrPtr LIBCALL AnnotDescrFree PROTO((AnnotDescrPtr anp));

typedef ValNode AnnotDesc, PNTR AnnotDescPtr;

NLM_EXTERN Boolean    LIBCALL AnnotDescAsnWrite PROTO((AnnotDescPtr anp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN AnnotDescPtr LIBCALL AnnotDescAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN AnnotDescPtr LIBCALL AnnotDescFree PROTO((AnnotDescPtr anp));
NLM_EXTERN Int2 LIBCALL AnnotDescLabel PROTO((AnnotDescPtr anp, CharPtr buffer, Int2 buflen, Boolean content));


typedef struct seqannot {
    AnnotIdPtr id;
    Uint1 db;
    CharPtr name;
    AnnotDescrPtr desc;
    Uint1 type;             /* 1=ftable, 2=align, 3=graph 4=ids 5=locs */
    Pointer data;
    struct seqannot PNTR next;
	GatherIndex idx;        /* internal gather/objmgr tracking fields */
} SeqAnnot, PNTR SeqAnnotPtr;

NLM_EXTERN SeqAnnotPtr LIBCALL SeqAnnotNew PROTO((void));
NLM_EXTERN Boolean     LIBCALL SeqAnnotAsnWrite PROTO((SeqAnnotPtr sap, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN SeqAnnotPtr LIBCALL SeqAnnotAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN SeqAnnotPtr LIBCALL SeqAnnotFree PROTO((SeqAnnotPtr sap));
NLM_EXTERN Int2 LIBCALL SeqAnnotLabel PROTO((SeqAnnotPtr sap, CharPtr buffer, Int2 buflen, Uint1 content));

/*****************************************************************************
*
*   Sets of SeqAnnots
*
*****************************************************************************/
NLM_EXTERN Boolean     LIBCALL SeqAnnotSetAsnWrite PROTO((SeqAnnotPtr sap, AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element));
NLM_EXTERN SeqAnnotPtr LIBCALL SeqAnnotSetAsnRead PROTO((AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element));

/** SeqAnnotSetExtraCheck looks for extra Features from SeqEntryAsnOut
    and adds them to the stream when the annot pointer is NULL
*********/
NLM_EXTERN Boolean LIBCALL SeqAnnotSetExtraCheck (AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element);


/*****************************************************************************
*
*   SeqHist
*
*****************************************************************************/
typedef struct seqhist {
	SeqAlignPtr assembly;
	DatePtr replace_date;
	SeqIdPtr replace_ids;
	DatePtr replaced_by_date;
	SeqIdPtr replaced_by_ids;
	Boolean deleted;
	DatePtr deleted_date;
} SeqHist, PNTR SeqHistPtr;

NLM_EXTERN SeqHistPtr LIBCALL SeqHistNew PROTO((void));
NLM_EXTERN Boolean    LIBCALL SeqHistAsnWrite PROTO((SeqHistPtr shp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN SeqHistPtr LIBCALL SeqHistAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN SeqHistPtr LIBCALL SeqHistFree PROTO((SeqHistPtr shp));


/**************************************************
*
*    LinkageEvidence
*
**************************************************/
typedef struct struct_Linkage_evidence {
   struct struct_Linkage_evidence PNTR next;
   Int4   type;
} LinkageEvidence, PNTR LinkageEvidencePtr;


NLM_EXTERN LinkageEvidencePtr LIBCALL LinkageEvidenceFree PROTO ((LinkageEvidencePtr ));
NLM_EXTERN LinkageEvidencePtr LIBCALL LinkageEvidenceNew PROTO (( void ));
NLM_EXTERN LinkageEvidencePtr LIBCALL LinkageEvidenceAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL LinkageEvidenceAsnWrite PROTO (( LinkageEvidencePtr , AsnIoPtr, AsnTypePtr));

/**************************************************
*
*    SeqGap
*
**************************************************/
typedef struct struct_Seq_gap {
   Int4        type;
   Int4        linkage;
   ValNodePtr  linkage_evidence;
} SeqGap, PNTR SeqGapPtr;


NLM_EXTERN SeqGapPtr LIBCALL SeqGapFree PROTO ((SeqGapPtr ));
NLM_EXTERN SeqGapPtr LIBCALL SeqGapNew PROTO (( void ));
NLM_EXTERN SeqGapPtr LIBCALL SeqGapAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SeqGapAsnWrite PROTO (( SeqGapPtr , AsnIoPtr, AsnTypePtr));


/*****************************************************************************
*
*   Delta-seq is just a ValNode where
*     vnp->choice = 1, vnp->data.ptrvalue = SeqLocPtr
*     vnp->choice = 2, vnp->data.ptrvalue = SeqLitPtr
*
*   SeqLit is a Seq-literal
*
*****************************************************************************/
typedef ValNode DeltaSeq, PNTR DeltaSeqPtr;

NLM_EXTERN Boolean    LIBCALL DeltaSeqAsnWrite PROTO((DeltaSeqPtr dsp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN DeltaSeqPtr LIBCALL DeltaSeqAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN DeltaSeqPtr LIBCALL DeltaSeqFree PROTO((DeltaSeqPtr dsp));

/*** this read/write/free a linked list of DeltaSeq *******/

NLM_EXTERN Boolean    LIBCALL DeltaSeqSetAsnWrite PROTO((DeltaSeqPtr dsp, AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element));
NLM_EXTERN DeltaSeqPtr LIBCALL DeltaSeqSetAsnRead PROTO((AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element));
NLM_EXTERN DeltaSeqPtr LIBCALL DeltaSeqSetFree PROTO((DeltaSeqPtr dsp));

/*** SeqDataPtr is an opaque pointer that can hold a byte store or a SeqGap - need to cast pointer to actual type ***/

struct SeqData;
typedef struct SeqData PNTR SeqDataPtr;

/*** here is the SeqLit ***********************************
*	if seq_data is NULL and length != 0, represents a gap of known length.
*   if seq_data is NULL and length = 0, represents a gap of unknown length
*   if seq_data	is NULL and length != 0, and fuzz is type lim, class unk,
*      then represents a gap of unknown length, but which has been spread
*      evenly across a whole delta seq to match BioseqPtr->length. See
*      SpreadGapsInDeltaSeq() for support to do this. Useful when the length
*      of a partially sequenced clone is know, but not specific gap lengths
***********************************************************/

typedef struct seqlit {
	Int4 length;
	IntFuzzPtr fuzz;
	Uint1 seq_data_type;    /* see Bioseq below */
	SeqDataPtr seq_data;    /* was ByteStorePtr */ 
} SeqLit, PNTR SeqLitPtr;

NLM_EXTERN SeqLitPtr LIBCALL SeqLitNew PROTO((void));
NLM_EXTERN Boolean    LIBCALL SeqLitAsnWrite PROTO((SeqLitPtr slp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN SeqLitPtr LIBCALL SeqLitAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN SeqLitPtr LIBCALL SeqLitFree PROTO((SeqLitPtr slp));

/*****************************************************************************
*
*   Bioseq.
*     Inst is incorporated within Bioseq for efficiency
*     seq_data_type
*       0 = not set
*       1 = IUPACna
*       2 = IUPACaa
*       3 = NCBI2na
*       4 = NCBI4na
*       5 = NCBI8na
*       6 = NCBIpna
*       7 = NCBI8aa
*       8 = NCBIeaa
*       9 = NCBIpaa
*      11 = NCBIstdaa
*      12 = Seq-gap
*   seq_ext_type
*       0 = none
*       1 = seg-ext
*       2 = ref-ext
*       3 = map-ext
*       4 = delta-ext
*
*****************************************************************************/
#define Seq_code_iupacna 1
#define Seq_code_iupacaa 2
#define Seq_code_ncbi2na 3
#define Seq_code_ncbi4na 4
#define Seq_code_ncbi8na 5
#define Seq_code_ncbipna 6
#define Seq_code_ncbi8aa 7
#define Seq_code_ncbieaa 8
#define Seq_code_ncbipaa 9
#define Seq_code_iupacaa3 10
#define Seq_code_ncbistdaa 11
#define Seq_code_gap 12

#define Seq_repr_virtual 1
#define Seq_repr_raw    2
#define Seq_repr_seg    3
#define Seq_repr_const  4
#define Seq_repr_ref    5
#define Seq_repr_consen 6
#define Seq_repr_map    7
#define Seq_repr_delta  8
#define Seq_repr_other  255

#define Seq_mol_dna     1
#define Seq_mol_rna     2
#define Seq_mol_aa      3
#define Seq_mol_na      4
#define Seq_mol_other   255
#define ISA_na(x)  ((x==1)||(x==2)||(x==4))
#define ISA_aa(x)  (x == 3)

typedef struct bioseq {
    SeqIdPtr id;          /* Seq-ids */
    ValNodePtr descr;              /* Seq-descr */
    Uint1 repr,
        mol;
    Int4 length;            /* -1 if not set */
    IntFuzzPtr fuzz;
    Uint1 topology,
        strand,
        seq_data_type,      /* as in Seq_code_type above */
        seq_ext_type;
    SeqDataPtr seq_data;    /* was ByteStorePtr */ 
    Pointer seq_ext;
    SeqAnnotPtr annot;
	SeqHistPtr hist;
	Pointer omdp;           /* internal objmgrdataptr to speed up indexed functions */
	GatherIndex idx;        /* internal gather/objmgr tracking fields */
	SeqEntryPtr seqentry;   /* internal seqentry that points to this bioseq */
} Bioseq, PNTR BioseqPtr;

NLM_EXTERN BioseqPtr LIBCALL BioseqNew PROTO((void));
NLM_EXTERN Boolean   LIBCALL BioseqAsnWrite PROTO((BioseqPtr bsp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN BioseqPtr LIBCALL BioseqAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN BioseqPtr LIBCALL BioseqFree PROTO((BioseqPtr bsp));
NLM_EXTERN Int2 LIBCALL BioseqLabel PROTO((BioseqPtr bsp, CharPtr buffer, Int2 buflen, Uint1 content));

/*****************************************************************************
*
*   BioseqFreeComponents(bsp)
*       Frees data associated with the Bioseq.
*   	Does not free SeqIds
*   	Does not free Bioseq itself
*   	Used by SeqMgr for caching out
*   	Called by BioseqFree()
*
*****************************************************************************/
NLM_EXTERN BioseqPtr LIBCALL BioseqFreeComponents PROTO((BioseqPtr bsp));

NLM_EXTERN Boolean LIBCALL BioseqInstAsnWrite PROTO((BioseqPtr bsp, AsnIoPtr aip, AsnTypePtr orig));
NLM_EXTERN Boolean LIBCALL BioseqInstAsnRead PROTO((BioseqPtr bsp, AsnIoPtr aip, AsnTypePtr orig));

NLM_EXTERN Boolean LIBCALL SeqDataAsnWrite PROTO((SeqDataPtr sdp, Uint1 seqtype, AsnIoPtr aip, AsnTypePtr orig));
NLM_EXTERN Boolean LIBCALL SeqDataAsnWriteXML PROTO((SeqDataPtr sdp, Uint1 seqtype, AsnIoPtr aip, AsnTypePtr orig,Int4 seqlen));
NLM_EXTERN SeqDataPtr LIBCALL SeqDataAsnRead PROTO((AsnIoPtr aip, AsnTypePtr orig, Uint1Ptr typeptr, Int4 length));
NLM_EXTERN SeqDataPtr LIBCALL SeqDataFree PROTO((SeqDataPtr sd, Uint1 seqtype));

/*****************************************************************************
*
*   Initialize bioseq and seqcode tables and default numbering
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL BioseqLoad PROTO((void));

/*****************************************************************************
*
*   BioseqAsnRead Options
*
*****************************************************************************/

typedef struct op_objseq {
	SeqIdPtr sip;          /* seq id to find */
	Boolean found_it;      /* set to TRUE when BioseqAsnRead matches sip */
	Boolean load_by_id;    /* if TRUE, load only if sip matches */
} Op_objseq, PNTR Op_objseqPtr;

/* types for AsnIoOption OP_NCBIOBJSEQ */
#define BIOSEQ_CHECK_ID 1    /* match Op_objseq.sip */
#define CHECK_EXTRA_DESC 2   /* used by SeqEntryAsnOut */
#define CHECK_EXTRA_FEAT 3   /* used by SeqEntryAsnOut */

/*****************************************************************************
*
*   SeqDescr uses an ValNode with choice = 
    1 = * mol-type GIBB-mol ,          -- type of molecule
    2 = ** modif SET OF GIBB-mod ,      -- modifiers
    3 = * method GIBB-method ,         -- sequencing method
    4 = name VisibleString ,         -- a name for this sequence
    5 = title VisibleString ,        -- a title for this sequence
    6 = org Org-ref ,                -- if all from one organism
    7 = comment VisibleString ,      -- a more extensive comment
    8 = num Numbering ,              -- a numbering system
    9 = maploc Dbtag ,               -- map location of this sequence
    10 = pir PIR-block ,              -- PIR specific info
    11 = genbank GB-block ,           -- GenBank specific info
    12 = pub Pubdesc                 -- a reference to the publication
    13 = region VisibleString        -- name for this region of sequence
    14 = user UserObject             -- user structured data object
    15 = sp SP-block                 -- SWISSPROT specific info

	<in release 4, neighbors is replaced with dbxref -- OBSOLETE     >
    <16 = neighbors Entrez-link       -- links to sequence neighbors >

	16 = dbxref Dbtag                -- xref to other databases
	17 = embl EMBL-block             -- EMBL specific info
	18 = create-date Date            -- date entry created
	19 = update-date Date            -- date of last update
	20 = prf PRF-block 				 -- PRF specific information
	21 = pdb PDB-block              -- PDB specific information
	22 = het Heterogen              -- cofactor, etc associated but not bound
	23 = source BioSource
	24 = molinfo MolInfo
    types with * use data.intvalue.  Other use data.ptrvalue
    ** uses a chain of ValNodes which use data.intvalue for enumerated type
*
*****************************************************************************/
#define Seq_descr_mol_type 1
#define Seq_descr_modif 2
#define Seq_descr_method 3
#define Seq_descr_name 4
#define Seq_descr_title 5
#define Seq_descr_org 6
#define Seq_descr_comment 7
#define Seq_descr_num 8
#define Seq_descr_maploc 9
#define Seq_descr_pir 10
#define Seq_descr_genbank 11
#define Seq_descr_pub 12
#define Seq_descr_region 13
#define Seq_descr_user 14
#define Seq_descr_sp 15
#define Seq_descr_dbxref 16
#define Seq_descr_embl 17
#define Seq_descr_create_date 18
#define Seq_descr_update_date 19
#define Seq_descr_prf 20
#define Seq_descr_pdb 21
#define Seq_descr_het 22
#define Seq_descr_source 23
#define Seq_descr_molinfo 24
#define Seq_descr_modelev 25

#define SEQDESCR_MAX 26 /* size of array needed for seqdescr filter parameters */

/* SeqDescr/SeqDesc will start to be allocated as ObjValNodes with vnp->extended
   set to 1 and EXTRA_OBJMGR_FIELDS fields available with a cast pointer */

typedef ValNode SeqDescr, FAR *SeqDescrPtr;

NLM_EXTERN Boolean    LIBCALL SeqDescrAsnWrite PROTO((SeqDescrPtr anp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN SeqDescrPtr LIBCALL SeqDescrAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN SeqDescrPtr LIBCALL SeqDescrFree PROTO((SeqDescrPtr anp));

   /*** these replace the SeqDescrEAsnWrite() function and add others */

typedef ValNode SeqDesc, FAR *SeqDescPtr;

NLM_EXTERN Boolean    LIBCALL SeqDescAsnWrite PROTO((SeqDescPtr anp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN SeqDescPtr LIBCALL SeqDescAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN SeqDescPtr LIBCALL SeqDescFree PROTO((SeqDescPtr anp));
NLM_EXTERN Int2 LIBCALL SeqDescLabel PROTO((SeqDescPtr vnp, CharPtr buffer, Int2 buflen, Boolean content));

/** SeqDescrExtraCheck looks for extra Descriptors from SeqEntryAsnOut
    and adds them to the stream when the descr pointer is NULL
*********/
NLM_EXTERN Boolean LIBCALL SeqDescrExtraCheck (AsnIoPtr aip, AsnTypePtr atp);

/*****************************************************************************
*
*   Pubdesc, Numbering, MolInfo, BioSource types defined in objpubd.h
*
*****************************************************************************/

#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif


