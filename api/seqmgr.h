/*  seqmgr.h
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
* File Name:  seqmgr.h
*
* Author:  James Ostell
*   
* Version Creation Date: 9/94
*
* $Revision: 6.73 $
*
* File Description:  Manager for Bioseqs and BioseqSets
*
* Modifications:  
* --------------------------------------------------------------------------
* Date	   Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*/

#ifndef _NCBI_SeqMgr_
#define _NCBI_SeqMgr_

#ifndef _NCBI_ObjMgr_
#include <objmgr.h>		   /* the object manager interface */
#endif

#ifndef _NCBI_Seqset_
#include <objsset.h>		   /* the object loader interface */
#endif

#ifndef __NLM_THR__
#include <ncbithr.h>
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

/*****************************************************************************
*
*   Sequence Management Functions
*
*****************************************************************************/
	/** callbacks for data management **/

/*****************************************************************************
*
*   SeqMgr manipulates the "registry" of Bioseqs in memory
*       assigns "EntityID" to SeqEntrys loaded in memory..
*   		only top level SeqEntry gets an EntityID
*   		caching and locking also done on top level SeqEntry
*
*****************************************************************************/

#define BSF_TEMP 1      /* for BioseqFetch functions */

/* smp is really a SeqMgrPtr, but had to be Pointer to satisfy compiler */

typedef BioseqPtr (LIBCALLBACK * BSFetchTop)
		PROTO((SeqIdPtr sip, Uint1 ld_type));

typedef BioseqPtr (LIBCALLBACK * BSFetch) PROTO((SeqIdPtr sip, Pointer data));

typedef Int4 (LIBCALLBACK * SIDPreCacheFunc) (SeqEntryPtr sep, Boolean components, Boolean locations, Boolean products, Boolean alignments, Boolean history, Boolean inference, Boolean others);
typedef Int4 (LIBCALLBACK * SeqLenLookupFunc) (Int4 gi);
typedef CharPtr (LIBCALLBACK * AccnVerLookupFunc) (Int4 gi);
typedef SeqIdPtr (LIBCALLBACK * SeqIdSetLookupFunc) (Int4 gi);

typedef struct seqidindexelement {
	CharPtr str;               /* PRINTID_FASTA_SHORT string */
	ObjMgrDataPtr omdp;             /* the omdp containing the Bioseq */
} SeqIdIndexElement, PNTR SeqIdIndexElementPtr;

typedef struct seqidindexblock {
	SeqIdIndexElement sid[100];
	struct seqidindexblock PNTR next;
} SeqIdIndexBlock, PNTR SeqIdIndexBlockPtr;

typedef struct smscope {       /* for setting scope by thread */
	TNlmThread thr;            /* the thread the scope is valid for */
	SeqEntryPtr SEscope;       /* scope for that thread */
} SMScope, PNTR SMScopePtr;

typedef struct seqmng {        /* functions for sequence data management */
	SMScopePtr scope;
	Int2 total_scope,              /* sizeof scope */
		num_scope;                 /* current number */
	BSFetchTop bsfetch;            /* BioseqFetch into memory */
	Pointer TopData;               /* user data for BSFetchTop function */
	Boolean fetch_on_lock;         /* call fetch when locking? */
	Int4 NonIndexedBioseqCnt,      /* number of Bioseqs in NonIndexedBioseq */
		NonIndexedBioseqNum;       /* size of NonIndexedBioseq */
	BioseqPtr PNTR NonIndexedBioseq; /* Bioseqs waiting for SeqId indexing */
	Int4 BioseqIndexCnt,           /* number of elements used in BioseqIndex */
		BioseqIndexNum;            /* size of BioseqIndex */
	SeqIdIndexElementPtr PNTR BioseqIndex;  /* pointers to index elements */
	SeqIdIndexBlockPtr BioseqIndexData;    /* what BioseqIndex points to */
	Boolean is_write_locked;
	Int4 hold_indexing;      /* set by SeqMgrHoldIndexing */
	SIDPreCacheFunc seq_id_precache_func;
	SeqLenLookupFunc seq_len_lookup_func;
	AccnVerLookupFunc accn_ver_lookup_func;
	SeqIdSetLookupFunc seq_id_set_lookup_func;
} SeqMgr, PNTR SeqMgrPtr;

/**** All replaced in Object Manager ************/
/************************************************/
#define SM_BIOSEQ OBJ_BIOSEQ
#define SM_BIOSEQSET OBJ_BIOSEQSET

#define SeqMgrConnect(a,b,c,d) ObjMgrConnect(a,b,c,d)

/*****************************************************************************
*
*   SeqMgrAdd(type, data)
*   	adds a Bioseq or BioseqSet to the sequence manager
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrAdd PROTO((Uint2 type, Pointer data));
/*****************************************************************************
*
*   SeqMgrDelete(type, data)
*   	deletes a Bioseq or BioseqSet from the sequence manager
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrDelete PROTO((Uint2 type, Pointer data));

/*****************************************************************************
*
*   SeqMgrHoldIndexing(Boolean hold)
*       stops sequence indexing to allow bulk loading if hold = TRUE
*       starts it when hold = FALSE;
*       uses a counter so you must call it the same number of times
*        with TRUE as with FALSE
*       when the counter decrements to 0, it will index what it has.
*
*****************************************************************************/
NLM_EXTERN void LIBCALL SeqMgrHoldIndexing PROTO((Boolean hold));

/*****************************************************************************
*
*   SeqMgrAddToBioseqIndex(bsp)
*   	Indexes a BioseqPtr by SeqId(s)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrAddToBioseqIndex PROTO((BioseqPtr bsp));

/*****************************************************************************
*
*   SeqMgrDeleteDeleteFromBioseqIndex(bsp)
*   	Removes index on BioseqPtr SeqIds
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrDeleteFromBioseqIndex PROTO((BioseqPtr bsp));

/*****************************************************************************
*
*   SeqMgrReplaceInBioseqIndex(bsp)
*   	Replaces index on BioseqPtr SeqIds
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrReplaceInBioseqIndex PROTO((BioseqPtr bsp));

/*****************************************************************************
*
*   SeqMgrDeleteIndexesInRecord (sep)
*   	Bulk removal of SeqId index on entire entity prior to its deletion
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrDeleteIndexesInRecord (SeqEntryPtr sep);

/*****************************************************************************
*
*   SeqMgrClearBioseqIndex()
*   	Clears entire SeqId index for all entities
*
*****************************************************************************/
NLM_EXTERN void SeqMgrClearBioseqIndex (void);


NLM_EXTERN Boolean LIBCALL SeqMgrSeqEntry PROTO((Uint1 type, Pointer data, SeqEntryPtr sep));
NLM_EXTERN SeqEntryPtr LIBCALL SeqMgrGetSeqEntryForData PROTO((Pointer data));
NLM_EXTERN Int2 LIBCALL SeqMgrGetEntityIDForSeqEntry PROTO((SeqEntryPtr sep));
NLM_EXTERN SeqEntryPtr LIBCALL SeqMgrGetSeqEntryForEntityID PROTO((Int2 id));

/*****************************************************************************
*
*   SeqIdFetch functions
*     Convert between id types
*        Look first in memory
*        Then call registered OBJ_SEQID,OBJ_SEQID proceedures to satisfy
*        EntrezBioseqFetchEnable supports these for the Entrez Interface
*
*****************************************************************************/

/*****************************************************************************
*
*   GetSeqIdForGI(Int4)
*     returns the SeqId for a GI
*     returns NULL if can't find it
*     The returned SeqId is allocated. Caller must free it.
*
*****************************************************************************/
NLM_EXTERN SeqIdPtr LIBCALL GetSeqIdForGI PROTO((Int4 gi));

/*****************************************************************************
*
*   GetSeqIdSetForGI(Int4)
*     returns the chain of all SeqIds for a GI
*     returns NULL if can't find it
*     The returned SeqId chain is allocated. Caller must free it with SeqIdSetFree.
*
*****************************************************************************/
NLM_EXTERN SeqIdPtr LIBCALL GetSeqIdSetForGI PROTO((Int4 gi));


/*****************************************************************************
*
*   GetGIForSeqId(SeqIdPtr)
*     returns the GI for a SeqId
*     returns 0 if can't find it
*
*****************************************************************************/
NLM_EXTERN Int4 LIBCALL GetGIForSeqId PROTO((SeqIdPtr sid));

/*****************************************************************************
*
*   MakeReversedSeqIdString(sid, buf, len)
*     Prints FASTA_SHORT style in upper case reverse order for fast binary searches
*
*****************************************************************************/
NLM_EXTERN Boolean MakeReversedSeqIdString (SeqIdPtr sid, CharPtr buf, size_t len);

/*****************************************************************************
*
*   SeqMgrLinkSeqEntry(sep, parenttype, parentptr)
*      connects all component seq-entries by traversing the linked list
*        all calling SeqMgrConnect and SeqMgrSeqEntry appropriately
*        if parenttype != 0, then assumes seqentry is contained in parentptr
*           and should be connected to it
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrLinkSeqEntry PROTO((SeqEntryPtr sep, Uint2 parenttype, Pointer parentptr));

/*****************************************************************************
*
*   ClearBioseqFindCache()
*   	frees internal omdp and se caches which can thwart detection of colliding IDs
*
*****************************************************************************/
NLM_EXTERN void ClearBioseqFindCache (void);

/*****************************************************************************
*
*   SeqMgrFreeCache()
*   	frees all cached SeqEntrys
*   	returns FALSE if any errors occurred
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrFreeCache PROTO((void));

/********************************************************************************
*
*   BioseqReload (omdp, lockit)
*     reloads the cached SeqEntry at top of omdp
*     if (lockit) locks the record
*
*     returns NULL on failure
*     returns omdp of (possibly new) top level ObjMgrData containing the reloaded
*      data from the old omdp. Also returns NULL if omdp does not have a Bioseq
*      fetch function attached to it for reload.
*
*********************************************************************************/
NLM_EXTERN ObjMgrDataPtr LIBCALL BioseqReload PROTO((ObjMgrDataPtr omdp, Boolean lockit));

/*****************************************************************************
*
*   Selection Functions for data objects based on SeqLoc
*      See also general selection in objmgr.h
*
*****************************************************************************/

/*****************************************************************************
*
*   SeqMgrSelect(region)
*      region is a SeqLocPtr
*          It can only apply to one Bioseq
*          selected area will be extreme left and right ends
*          fuzziness is ignored
*      if something else selected, deselects it first, then selects requested
*        item
*      to select without deselecting something else, use SeqMgrAlsoSelect()
*      returns TRUE if item is now currently selected, FALSE if not
*      "region" is always copied. Caller is responsible for managment of
*         SeqLoc that is passed in.
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrSelect PROTO((SeqLocPtr region));
NLM_EXTERN Boolean LIBCALL SeqMgrAlsoSelect PROTO((SeqLocPtr region));

/*****************************************************************************
*
*   SeqMgrDeSelect(region)
*   	if this item was selected, then deselects and returns TRUE
*   	else returns FALSE
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrDeSelect PROTO((SeqLocPtr region));

/*****************************************************************************
*
*   SeqMgrSetColor(region, rgb)
*      region is a SeqLocPtr
*          It can only apply to one Bioseq
*          colored area will be extreme left and right ends
*          fuzziness is ignored
*      "region" is always copied. Caller is responsible for managment of
*         SeqLoc that is passed in.
*      rgb is a Uint1[3] array with RGB values
*         it is always copied so Caller is responsible for memory management
*         of passed in object
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrSetColor PROTO((SeqLocPtr region, Uint1Ptr rgb));

/************************************************/
/************************************************/


/*****************************************************************************
*
*   Return the current SeqMgr
*   	Initialize if not done already
*       This function will become obsolete
*
*****************************************************************************/
NLM_EXTERN SeqMgrPtr LIBCALL SeqMgrGet (void);

/*****************************************************************************
*
*   SeqMgrReadLock()
*   	Initialize if not done already
*       A thread can have only one read or write lock at a time
*       Many threads can have read locks
*       Only one thread can have a write lock
*       No other threads may have read locks if a write lock is granted
*       If another thread holds a write lock, this call blocks until write
*          is unlocked.
*
*****************************************************************************/
NLM_EXTERN SeqMgrPtr LIBCALL SeqMgrReadLock (void);

/*****************************************************************************
*
*   SeqMgrWriteLock
*   	Initialize if not done already
*       A thread can have only one read or write lock at a time
*       Many threads can have read locks
*       Only one thread can have a write lock
*       No other threads may have read locks if a write lock is granted
*       If another thread holds a read or write lock, this call blocks until write
*          is unlocked.
*
*****************************************************************************/
NLM_EXTERN SeqMgrPtr LIBCALL SeqMgrWriteLock (void);

/*****************************************************************************
*
*  SeqMgrUnlock()
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrUnlock (void);

/*****************************************************************************
*
*   SeqMgrSetBSFetchTop (fetch, data)
*   	sets the BSFetchTop routine to "fetch"
*       "data" pointer will be sent to function
*       returns previous value
*       set to NULL to turn off all fetching for that type
*
*       Current value can be called directly as BioseqFetch();
*   	Default is
*   		1) looks in memory
*   		2) looks locally if LocalBSFetch is set
*			3) looks remotely if RemoteBSFetch is set
*
*****************************************************************************/
NLM_EXTERN BSFetchTop LIBCALL SeqMgrSetBSFetchTop PROTO((BSFetchTop fetch, Pointer data));

/*****************************************************************************
*
*   SeqMgrSetFetchOnLock(value)
*   	if value = TRUE, manager will try to fetch the bioseq if not in
*          memory, when BioseqLock is called
*   	if FALSE, BioseqLock will only look in memory
*       returns previous value of fetch_on_lock
*       default is to fetch_on_lock
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqMgrSetFetchOnLock PROTO((Boolean value));



NLM_EXTERN BioseqPtr LIBCALL BioseqLock PROTO((BioseqPtr bsp));
NLM_EXTERN Boolean LIBCALL BioseqUnlock PROTO((BioseqPtr bsp));

/*****************************************************************************
*
*   BioseqFind(sid)
*   	Finds a Bioseq in memory based on SeqId
*   	Will also restore a Bioseq that has been cached out by SeqMgr
*
*****************************************************************************/
NLM_EXTERN BioseqPtr LIBCALL BioseqFind PROTO((SeqIdPtr sip));

/*****************************************************************************
*
*   BioseqFindCore(sid)
*   	Finds a Bioseq in memory based on SeqId when only "core" elements needed
*   	Will NOT restore a Bioseq that has been cached out by SeqMgr
*       This function is for use ONLY by functions that only need the parts
*         of the Bioseq left when cached out. This includes the SeqId chain,
*         and non-pointer components of the Bioseq.
*
*****************************************************************************/
NLM_EXTERN BioseqPtr LIBCALL BioseqFindCore PROTO((SeqIdPtr sip));

/*****************************************************************************
*
*   BioseqFindSpecial(sid)
*   	Finds a Bioseq in memory based on SeqId when only "core" elements needed
*   	Will NOT restore a Bioseq that has been cached out by SeqMgr
*       This function does not use the bioseq_cache mechanism, and is for
*         the validator to check for IdOnMultipleBioseqs.
*
*****************************************************************************/
NLM_EXTERN BioseqPtr LIBCALL BioseqFindSpecial (SeqIdPtr sip);

/*****************************************************************************
*
*   BioseqFindEntity(sid, itemIDptr)
*   	Finds a Bioseq in memory based on SeqId
*   	Will NOT restore a Bioseq that has been cached out by SeqMgr
*       returns EntityID if found, otherwise 0
*       itemIDptr is set to the value for itemID in ObjMgr functions
*       itemtype is OBJ_BIOSEQ of course
*
*****************************************************************************/
NLM_EXTERN Uint2 LIBCALL BioseqFindEntity PROTO((SeqIdPtr sip, Uint4Ptr itemIDptr));


/*****************************************************************************
*
*   BioseqLockById(sid)
*   	Like BioseqFind, except will also try to fetch the bioseq from
*         outside storage if not in memory already. Will cache out data
*   	  loaded this way if memory gets too full
*         Calls BioseqFetch to do the fetch
*
*****************************************************************************/
NLM_EXTERN BioseqPtr LIBCALL BioseqLockById PROTO((SeqIdPtr sid));
NLM_EXTERN Boolean LIBCALL BioseqUnlockById PROTO((SeqIdPtr sid));

NLM_EXTERN BioseqPtr LIBCALL BioseqFetch PROTO((SeqIdPtr sid, Uint1 ld_type));

#define BSFETCH_TEMP 1	   /* load called by software.. temporary use */
#define BSFETCH_STD  0	   /* load called by user, must be freed by user */

/*****************************************************************************
*
*   SeqEntry Management Functions
*
*****************************************************************************/
/*****************************************************************************
*
*   SeqEntrySetScope(sep)
*   	scopes global seqentry searches to sep
*       setting sep=NULL, opens scope to all seqentries in memory
*       returns the current scope
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr LIBCALL SeqEntrySetScope PROTO((SeqEntryPtr sep));

/*****************************************************************************
*
*   SeqEntryGetScope(sep)
*       returns the current scope or NULL if none set
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr LIBCALL SeqEntryGetScope PROTO((void));

NLM_EXTERN SeqEntryPtr LIBCALL SeqEntryFind PROTO((SeqIdPtr sip));

/*****************************************************************************
*
*   Context routines for Bioseqs in Seq-entrys
*      Context is the chain of Seqentries leading to the bioseq.
*      context[count-1] is SeqEntry for bsp itself
*      If Bioseq not in a Seqentry, count is 0 and bcp->se may be used
*        if a fake Seqentry is convenient.
*
*****************************************************************************/
#define BIOSEQCONTEXTMAX 20

typedef struct bioseqcontxt {
	BioseqPtr bsp;           /* the Bioseq in question */
	Int2 count;              /* number of elements in context */
	Boolean hit;             /* used by BioseqContextNew and ..GetSeqFeat */
	SeqEntryPtr context[BIOSEQCONTEXTMAX];  /* array of SeqEntryPtr (last is count -1) */
	ValNode se;             /* used for a tempory SeqEntryPtr when only a Bioseq */
	SeqFeatPtr sfp;          /* current sfp */
	SeqAnnotPtr sap;         /* current sap */
	Int2 sftype,             /* SeqFeat type to look for */
		in;					 /* 0=location, 1=product, 2=either */
} BioseqContext, PNTR BioseqContextPtr;

NLM_EXTERN BioseqContextPtr LIBCALL BioseqContextNew PROTO((BioseqPtr bsp));
NLM_EXTERN BioseqContextPtr LIBCALL BioseqContextFree PROTO((BioseqContextPtr bcp));
/*****************************************************************************
*
*   BioseqContextGetSeqDescr(bcp, type, curr, SeqEntryPtr PNTR sep)
*       returns pointer to the next SeqDescr of this type
*       type gives type of Seq-descr
*       if (type == 0)
*          get them all
*       curr is NULL or previous node of this type found
*       moves up from bsp
*		if (sep != NULL) sep set to SeqEntryPtr containing the SeqDescr.
*
*****************************************************************************/
NLM_EXTERN ValNodePtr LIBCALL BioseqContextGetSeqDescr PROTO((BioseqContextPtr bcp, Int2 type, ValNodePtr curr, SeqEntryPtr PNTR the_sep));
NLM_EXTERN CharPtr LIBCALL BioseqContextGetTitle PROTO((BioseqContextPtr bcp));
/*****************************************************************************
*
*   BioseqContextGetSeqFeat(bcp, type, curr, sapp, in)
*       returns pointer to the next Seq-feat of this type
*       type gives type of Seq-descr
*       if (type == 0)
*          get them all
*       curr is NULL or previous node of this type found
*       moves up from bsp
*   	if (sapp != NULL) is filled with SeqAnnotPtr containing the SeqFeat
*   	in:
*   		0 = sfp->location only
*   		1 = sfp->product only
*   		2 = either of above
*
*****************************************************************************/
NLM_EXTERN SeqFeatPtr LIBCALL BioseqContextGetSeqFeat PROTO((BioseqContextPtr bcp, Int2 type, SeqFeatPtr curr, SeqAnnotPtr PNTR sapp, Int2 in));

/*** works like BioseqContextGetSeqFeat, but a SeqEntryPtr and (optionally)
     a Bioseq will do ****************************************************/

/*****************************************************************************
*
*   SeqEntryGetSeqFeat(sep, type, curr, sapp)
*       returns pointer to the next Seq-feat of this type
*       type gives type of SeqFeat
*       if (type == 0)
*          get them all
*       curr is NULL or previous node of this type found
*       moves up from bsp
*   	if (sapp != NULL) is filled with SeqAnnotPtr containing the SeqFeat
*       if (bsp != NULL) then for that Bioseq match on location by "in"
*   	in:
*   		0 = sfp->location only
*   		1 = sfp->product only
*   		2 = either of above
*
*****************************************************************************/
NLM_EXTERN SeqFeatPtr LIBCALL SeqEntryGetSeqFeat PROTO((SeqEntryPtr sep, Int2 type, SeqFeatPtr curr, SeqAnnotPtr PNTR sapp, Int2 in, BioseqPtr bsp));

/*****************************************************************************
*
*   SpreadGapsInDeltaSeq(BioseqPtr bsp)
*      bsp must be a delta seq
*      function counts deltas with known lengths ( = known_len)
*               counts deltas which are gaps of unknown length ( = unk_count)
*                  these can delta of length 0, delta with fuzz = lim (unk),
*                    or SEQLOC_NULL
*               converts all unknown gaps to delta with fuzz = lim(unk)
*               sets length of all unknown gaps to
*                  (bsp->length - known_len)/unk_count
*                  any reminder spread over first few gaps
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SpreadGapsInDeltaSeq PROTO((BioseqPtr bsp));

/*****************************************************************************
*
*   CountGapsInDeltaSeq(BioseqPtr bsp, &num_segs, &num_gaps, &known_residues, &num_gaps_faked, CharPtr buf, Int2 buflen)
*      bsp must be a delta seq
*      function counts deltas and returns a profile
*          num_segs = total number of segments
*          num_gaps = total number of segments representing gaps
*          known_residues = number of real residues in the sequence (not gaps)
*          num_gaps_faked = number of gaps where real length is not known, but where
*                           a length was guessed by spreading the total gap length
*                           out over all gaps evenly.
*
*      NOTE: any of these pointers except bsp can be NULL
*
*      returns TRUE if values in argument were set.
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL CountGapsInDeltaSeq PROTO((BioseqPtr bsp, Int4Ptr num_segs, Int4Ptr num_gaps, Int4Ptr known_residues, Int4Ptr num_gaps_faked, CharPtr buf, Int4 buflen));

/*****************************************************************************
*
*   IsNonGappedLiteral(BioseqPtr bsp)
*      Returns TRUE if bsp is a delta seq is composed only of Seq-lits with
*      actual sequence data.  These are now made to allow optimal compression
*      of otherwise raw sequences with runs of ambiguous bases.
*
*****************************************************************************/
NLM_EXTERN Boolean IsNonGappedLiteral (BioseqPtr bsp);

/*****************************************************************************
*
*   GetUniGeneIDForSeqId(SeqIdPtr)
*     returns the UniGene ID for a SeqId
*     returns 0 if can't find it, or not a legal unigene id
*     This only applies to genomes division of entrez
*     These serve as temporary placeholders until NCBI establishes
*       a stable long-term ID system for these sequence clusters
*
*     The clusters begin with 1,000,000 and are grouped by organism
*
*****************************************************************************/

NLM_EXTERN Int4 LIBCALL GetUniGeneIDForSeqId PROTO((SeqIdPtr sip));

/*****************************************************************************
*
*   FetchFromSeqIdGiCache(gi, sipp)
*   FetchFromGiSeqIdCache(sip, gip)
*   RecordInSeqIdGiCache(gi, sip)
*   FreeSeqIdGiCache()
*     Internal functions to cache gi - SeqId associations
*
*****************************************************************************/
NLM_EXTERN Boolean FetchFromSeqIdGiCache (Int4 gi, SeqIdPtr PNTR sipp);
NLM_EXTERN Boolean FetchFromGiSeqIdCache (SeqIdPtr sip, Int4Ptr gip);

NLM_EXTERN void RecordInSeqIdGiCache (Int4 gi, SeqIdPtr sip);
NLM_EXTERN void FreeSeqIdGiCache (void);


/*****************************************************************************
*
*   BioseqExtra extensions to object manager to preindex features for rapid retrieval
*
*   Public functions moved to explore.h
*
*****************************************************************************/

/* the following structures are not frequently used directly by applications */

typedef struct smfeatitem {
  SeqFeatPtr   sfp;          /* freed when TL_CACHED, later will implement reassignment when reloaded */
  SeqAnnotPtr  sap;          /* SeqAnnot containing SeqFeat, same reap/reload criteria as above */
  BioseqPtr    bsp;          /* Bioseq on which this feature is indexed */
  CharPtr      label;        /* featdef content label */
  Int4         left;         /* extreme left on bioseq (first copy spanning origin is < 1) */
  Int4         right;        /* extreme right on bioseq (second copy spanning origin is > length) */
  Int4Ptr      ivals;        /* array of start/stop pairs */
  Int2         numivals;     /* number of start/stop pairs in ivals array */
  Int4         dnaStop;      /* last stop on protein mapped to DNA coordinate for flatfile */
  Boolean      partialL;     /* left end is partial */
  Boolean      partialR;     /* right end is partial */
  Boolean      farloc;       /* location has an accession not packaged in entity */
  Boolean      bad_order;    /* location is out of order - possibly trans-spliced */
  Boolean      mixed_strand; /* location has mixed strands - possibly trans-spliced */
  Boolean      ts_image;     /* trans-spliced image on another small chromosome, ignore packaging error */
  Uint1        strand;       /* strand (mapped to segmented bioseq if segmented) */
  Uint1        subtype;      /* featdef subtype */
  Uint4        itemID;       /* storing itemID so no need to gather again */
  Boolean      ignore;       /* ignore this second copy of a feature spanning the origin */
  Uint4        index;        /* position index needed for SeqMgrGetDesiredFeature */
  Int4         overlap;      /* for xxxByPos, index of leftmost candidate that overlaps this */
} SMFeatItem, PNTR SMFeatItemPtr;

typedef struct smfeatblock {
  struct smfeatblock PNTR  next;   /* pointer to next block of chunks */
  Uint2                    index;  /* latest offset within this block */
  SMFeatItemPtr            data;   /* allocated block for this chunk */
} SMFeatBlock, PNTR SMFeatBlockPtr;

typedef struct segpartsmap {
  struct segpartsmap PNTR  next;           /* pointer to next block of chunks */
  SeqLocPtr                slp;            /* allocated copy of seqLoc for part */
  CharPtr                  seqIdOfPart;    /* reverse upper case string of seqID of part */
  BioseqPtr                parentBioseq;   /* bioseq pointer for segmented parent */
  Int4                     cumOffset;      /* offset of part in segmented bioseq */
  Int4                     from;
  Int4                     to;
  Uint1                    strand;
  Uint4                    itemID;         /* OBJ_BIOSEQ_SEG itemID */
} SMSeqIdx, PNTR SMSeqIdxPtr;

typedef struct smdescitem {
  SeqDescrPtr  sdp;         /* freed when TL_CACHED, later will implement reassignment when reloaded */
  SeqEntryPtr  sep;         /* SeqEntry containing SeqDescr, same reap/reload criteria as above */
  Uint4        itemID;      /* storing itemID so no need to gather again */
  Uint4        index;       /* position index needed for SeqMgrGetDesiredDescriptor */
  Uint2        level;       /* packaging level - 0 is on Bioseq itself */
  Uint1        seqdesctype; /* seqdesc subtype */
} SMDescItem, PNTR SMDescItemPtr;

typedef struct smfiditem {
  SMFeatItemPtr  feat;
  CharPtr        fid;       /* string with numeric or alpha local feature ID */
} SMFidItem, PNTR SMFidItemPtr;

typedef struct bioseqextra {
  BioseqPtr           bsp;
  ObjMgrDataPtr       omdp;
  SeqFeatPtr          protFeat;        /* protein feature on whole protein bioseq gives name */
  SeqFeatPtr          cdsOrRnaFeat;    /* cds or rna whose product points to this bioseq */
  ValNodePtr          prodlisthead;    /* all features whose product points to this bioseq */

  SMFeatBlockPtr      featlisthead;    /* linked list of SMFeatItem chunks, arrays point to elements */
  SMFeatBlockPtr      featlisttail;    /* current block in linked list of SMFeatItem chunks */

  ValNodePtr          desclisthead;    /* linked list of ValNodes pointing to SMDescItem structures */

  SMDescItemPtr PNTR  descrsByID;      /* array of all descriptors on bioseq in original itemID order */
  SMDescItemPtr PNTR  descrsBySdp;     /* array of all features on bioseq sorted by SeqDescrPtr */
  SMDescItemPtr PNTR  descrsByIndex;   /* array of all features on bioseq sorted by order of presentation */

  AnnotDescPtr PNTR   annotDescByID;   /* array of all AnnotDesc (on entity) in original itemID order */

  SeqAlignPtr PNTR    alignsByID;      /* array of all alignments (on entity) in original itemID order */

  SMFeatItemPtr PNTR  featsByID;       /* array of all features on bioseq in original itemID order */
  SMFeatItemPtr PNTR  featsBySfp;      /* array of all features on bioseq sorted by SeqFeatPtr */
  SMFeatItemPtr PNTR  featsByPos;      /* array of all features on bioseq sorted by location */
  SMFeatItemPtr PNTR  featsByRev;      /* array of all features on bioseq sorted by reverse location */
  SMFeatItemPtr PNTR  featsByLabel;    /* array of all features on bioseq sorted by label */

  SMFeatItemPtr PNTR  genesByPos;      /* subset of featsByPos array containing only gene features */
  SMFeatItemPtr PNTR  mRNAsByPos;      /* subset of featsByPos array containing only mRNA features */
  SMFeatItemPtr PNTR  CDSsByPos;       /* subset of featsByPos array containing only CDS features */
  SMFeatItemPtr PNTR  pubsByPos;       /* subset of featsByPos array containing only publication features */
  SMFeatItemPtr PNTR  orgsByPos;       /* subset of featsByPos array containing only biosource features */
  SMFeatItemPtr PNTR  operonsByPos;    /* subset of featsByPos array containing only operon features */
  SMFeatItemPtr PNTR  genesByLocusTag; /* array of gene features sorted by locus_tag */

  SMFidItemPtr PNTR   featsByFeatID;   /* array of features sorted by feature ID string */

  BioseqPtr           parentBioseq;    /* segmented parent of this raw part all packaged together */
  SMSeqIdxPtr         segparthead;     /* linked list to speed mapping from parts to segmented bioseq */

  SMSeqIdxPtr PNTR    partsByLoc;      /* array of parts on segmented bioseq sorted by location */
  SMSeqIdxPtr PNTR    partsBySeqId;    /* array of parts on segmented bioseq sorted by reverse uppercase seqID */

  Int4                numdescs;        /* number of elements in descrsByID, descrsBySdp, and descrsByIndex arrays */
  Int4                numannotdesc;    /* number of elements in annotDescByID array */
  Int4                numaligns;       /* number of elements in alignsByID array */
  Int4                numfeats;        /* number of elements in featsByID, featsBySfp and featsByPos arrays */
  Int4                numgenes;        /* number of elements in genesByPos array */
  Int4                nummRNAs;        /* number of elements in mRNAsByPos array */
  Int4                numCDSs;         /* number of elements in CDSsByPos array */
  Int4                numpubs;         /* number of elements in pubsByPos array */
  Int4                numorgs;         /* number of elements in orgsByPos array */
  Int4                numoperons;      /* number of elements in operonsByPos array */
  Int4                numfids;         /* number of elements in featsByFeatID array */

  Int4                numsegs;         /* number of segments in partslist array */

  Int4                min;             /* used for finding best protein feature */
  Uint1               processed;       /* also used for finding best protein feature */
  Uint4               bspItemID;       /* for bioseq explore functions */
  Uint4               bspIndex;        /* for bioseq explore functions */
  Int2                blocksize;       /* size of SMFeatBlock.data array to avoid wasting space */
                                       /* additional fields to map between genome record and parts,
                                          genomic DNA and mRNA, and mRNA and protein */
} BioseqExtra, PNTR BioseqExtraPtr;

/* the following functions are not frequently called by applications */

/*****************************************************************************
*
*   Bioseq extra functions for reapextra, reloadextra, and freeextra take an
*     ObjMgrDataPtr as a parameter, and are only called by the object manager,
*     not the application program
*
*****************************************************************************/

NLM_EXTERN Pointer LIBCALLBACK SeqMgrReapBioseqExtraFunc PROTO((Pointer data));
NLM_EXTERN Pointer LIBCALLBACK SeqMgrReloadBioseqExtraFunc PROTO((Pointer data));
NLM_EXTERN Pointer LIBCALLBACK SeqMgrFreeBioseqExtraFunc PROTO((Pointer data));

/*****************************************************************************
*
*   SeqMgrFindSMFeatItemPtr and SeqMgrFindSMFeatItemByID return SMFeatItemPtr
*     to access internal fields, passing entityID and not bsp uses list attached
*     to top of entity containing index to all feature itemIDs regardless of
*     what bioseq they are indexed on
*   SeqMgrGetDesiredFeature in explore.h is the preferred public function
*   SeqMgrGetSfpProductList returns linked list of features whose sfp->product
*     points to the given bioseq
*
*****************************************************************************/

NLM_EXTERN SMFeatItemPtr LIBCALL SeqMgrFindSMFeatItemPtr PROTO((SeqFeatPtr sfp));
NLM_EXTERN SMFeatItemPtr LIBCALL SeqMgrFindSMFeatItemByID PROTO((Uint2 entityID, BioseqPtr bsp, Uint4 itemID));
NLM_EXTERN ValNodePtr LIBCALL SeqMgrGetSfpProductList (BioseqPtr bsp);

/*****************************************************************************
*
*   SeqMgrMapPartToSegmentedBioseq can speed up sequtil's CheckPointInBioseq
*     for indexed part bioseq to segmented bioseq mapping
*
*****************************************************************************/

NLM_EXTERN Int4 LIBCALL SeqMgrMapPartToSegmentedBioseq PROTO((BioseqPtr in, Int4 pos, BioseqPtr bsp, SeqIdPtr sip, BoolPtr flip_strand));

/*****************************************************************************
*
*   GenomePartToSegmentMap used for mapping of part positions not used on a given contig
*
*****************************************************************************/

NLM_EXTERN SMSeqIdxPtr GenomePartToSegmentMap (BioseqPtr in, BioseqPtr bsp, Int4 from, Int4 to);

/*****************************************************************************
*
*   TrimLocInSegment takes a location on an indexed far segmented part and trims
*     trims it to the region referred to by the parent segmented or delta bioseq.
*
*****************************************************************************/

NLM_EXTERN SeqLocPtr TrimLocInSegment (
  BioseqPtr master,
  SeqLocPtr location,
  BoolPtr p5ptr,
  BoolPtr p3ptr
);

/*****************************************************************************
*
*   SeqMgrIndexAlignments called by SeqMgrIndexFeatures, can be called separately
*
*****************************************************************************/

NLM_EXTERN void LIBCALL SeqMgrIndexAlignments (Uint2 entityID);

/*****************************************************************************
*
*   SeqMgrFindAnnotDescByID and SeqMgrFindSeqAlignByID uses new indexes to speed
*     lookup of AnnotDescPtr and SeqAlignPtr, respectively
*
*****************************************************************************/

NLM_EXTERN AnnotDescPtr LIBCALL SeqMgrFindAnnotDescByID (Uint2 entityID, Uint4 itemID);
NLM_EXTERN SeqAlignPtr LIBCALL SeqMgrFindSeqAlignByID PROTO((Uint2 entityID, Uint4 itemID));

/*****************************************************************************
*
*   LockFarComponents takes a top SeqEntryPtr and locks the far Bioseq components of
*     any segmented or delta sequences.  It turns a ValNode list of locked BioseqPtrs.
*   LockFarComponentsEx takes a top SeqEntryPtr and locks far Bioseqs that are either
*     components of any segmented or delta sequences, pointed to by feature locations,
*     or pointed to by feature products.  It turns a ValNode list of locked BioseqPtrs.
*   UnlockFarComponents takes the ValNode list of locked BioseqPtrs, unlocks each
*     Bioseq, frees the ValNode list, and returns NULL.
*   AdvcLockFarComponents is redesigned for better efficiency, can be multithreaded.
*
*****************************************************************************/

NLM_EXTERN ValNodePtr LockFarComponents (SeqEntryPtr sep);

NLM_EXTERN ValNodePtr LockFarComponentsEx (SeqEntryPtr sep, Boolean components, Boolean locations, Boolean products, SeqLocPtr loc);

NLM_EXTERN ValNodePtr UnlockFarComponents (ValNodePtr bsplist);

NLM_EXTERN ValNodePtr AdvcLockFarComponents (
  SeqEntryPtr sep,
  Boolean components,
  Boolean locations,
  Boolean products,
  SeqLocPtr loc,
  Boolean usethreads
);

/*****************************************************************************
*
*   LockFarAlignmentBioseqs finds Bioseqs in an alignment that could be
*   cached out and locks them.  It returns a ValNode list of the Bioseqs 
*   that can be unlocked with UnlockFarComponents.
*
*****************************************************************************/
NLM_EXTERN ValNodePtr LockFarAlignmentBioseqs (SeqAlignPtr salp);

/*****************************************************************************
*
*   SeqMgrSetPreCache registers the GiToSeqID precache function
*   LookupFarSeqIDs calls any registered function to preload the cache
*
*****************************************************************************/

NLM_EXTERN void LIBCALL SeqMgrSetPreCache (SIDPreCacheFunc func);

NLM_EXTERN Int4 LookupFarSeqIDs (
  SeqEntryPtr sep,
  Boolean components,
  Boolean locations,
  Boolean products,
  Boolean alignments,
  Boolean history,
  Boolean inference,
  Boolean others
);

/*****************************************************************************
*
*   SeqMgrSetLenFunc registers the GiToSeqLen lookup function
*   SeqMgrSetAccnVerFunc registers the GiToAccnVer lookup function
*   SeqMgrSetSeqIdSetFunc registers the GiToSeqIdSet lookup function
*
*****************************************************************************/

NLM_EXTERN void LIBCALL SeqMgrSetLenFunc (SeqLenLookupFunc func);

NLM_EXTERN void LIBCALL SeqMgrSetAccnVerFunc (AccnVerLookupFunc func);

NLM_EXTERN void LIBCALL SeqMgrSetSeqIdSetFunc (SeqIdSetLookupFunc func);

/*****************************************************************************
*
*   SeqEntryAsnOut (SeqEntryPtr sep, SeqIdPtr sip, Int2 retcode, AsnIoPtr aipout)
*
*      Takes top level SeqEntryPtr (sep) in memory, finds the Bioseq for sip,
*         then writes the relevant part of the SeqEntry into aipout depending
*         on retcode.. The SeqEntry in memory is unchanged.
*
*      retcode sets maximum complexity to return by values:
*         0 = return the whole blob
*         1 = return just the Bioseq and relevant descriptors and features
*         2 = return containing Seg-set if any
*         3 = return containing Nuc-prot set if any
*         4 = return containing Pub-set if any (this no longer used really)
*
******************************************************************************/
NLM_EXTERN Boolean SeqEntryAsnOut (SeqEntryPtr sep, SeqIdPtr sip,
                                    Int2 retcode, AsnIoPtr aipout);


NLM_EXTERN ObjMgrDataPtr SeqMgrGetOmdpForBioseq (BioseqPtr bsp);

NLM_EXTERN Pointer SeqMgrGetExtraDataForOmdp (ObjMgrDataPtr omdp);

NLM_EXTERN void SeqMgrRedoDescriptorIndexes (Uint2 entityID, Pointer ptr);
NLM_EXTERN void SeqMgrRedoFeatByLabelIndexes (Uint2 entityID, Pointer ptr);


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

