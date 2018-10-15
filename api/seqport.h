/*  seqport.h
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
* File Name:  seqport.h
*
* Author:  James Ostell
*   
* Version Creation Date: 7/13/91
*
* $Revision: 6.64 $
*
* File Description:  Ports onto Bioseqs
*
* Modifications:  
* --------------------------------------------------------------------------
* Date	   Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*/

#ifndef _NCBI_Seqport_
#define _NCBI_Seqport_

#include <sequtil.h>

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
*   SeqPort
*       will attach only to a Bioseq (SeqPortNew) or to a Seq-loc
*          (SeqPortNewByLoc) in any selected alphabet
*       You can then treat the sequence or location as a single contiguous
*          piece. You can Seek (SeqPortSeek) to any location. You can
*          SeqPortGetResidue, which get the residue at the current position
*          and seeks to the next residue. You can read a buffer of residues.
*   
*       Special characters are returned from SeqPortGetResidue
*   		SEQPORT_EOF   (end of sequence reached)
*           SEQPORT_VIRT  (hit a virtual sequence or gap)
*           INVALID_RESIDUE (residue not valid in original Bioseq)
*   		SEQPORT_EOS   (end of segment, not normally seen)
*
*       Some defined values are provided for the Int4 values passed as
*           start or stop
*       FIRST_RESIDUE   0	(first residue of sequence)
*       LAST_RESIDUE    -1  (last residue of sequence.. interpreted as
*                              length - 1)
*       APPEND_RESIDUE  -2  (interpreted as length.. off the end of the
*                              sequence. Not valid for SeqPort.. only
*                              used by editing functions )
*            
*          
*
*****************************************************************************/
#define SEQPORT_EOF 253       /* end of sequence data */
#define SEQPORT_EOS 252       /* end of segment */
#define SEQPORT_VIRT 251      /* skipping virtual sequence or gap */
#define IS_residue(x) (x <= 250)

#define FIRST_RESIDUE 0
#define LAST_RESIDUE -1
#define APPEND_RESIDUE -2

typedef struct spcache {
	Int2 ctr, total;
	Uint1 buf[100];
} SPCache, PNTR SPCachePtr;

typedef struct spcacheq {
	Int2 ctr, total;
	Char buf[400];
} SPCacheQ, PNTR SPCacheQPtr;

typedef struct seqport {
    BioseqPtr bsp;             /* 1 seqentry per port */
	Boolean locked;              /* TRUE if Lock function used */
    Int4 start, stop,            /* region of bsp covered */
        curpos,                  /* current position 0-(totlen-1) */
        totlen,                  /* total length of covered region */
		bytepos;                 /* current byte position in bsp->data */
    NumberingPtr currnum;        /* current numbering info */
    Uint1 strand,                /* as in seqloc */
        lastmsg;                 /* used by SeqPortRead() */
    Boolean is_circle ,          /* go around the end of a circle? */
            is_seg ,             /* return EOS at the end of segments? */
            do_virtual,          /* deliver 'N''X' over virtual seqs */
            gapIsZero,           /* deliver 0 for ncbi4na over virtual seqs - also needs do_virtual */
            eos,                 /* set when comp strand tries to back off */
			isa_null,            /* TRUE if seqport represents a NULL location */
			isa_virtual,         /* represents a virtual interval or Bioseq */
			backing;             /* signal to SeqPortSeek for backing up a layered SeqPort */
    SeqMapTablePtr smtp;         /* for mapping to requested alphabet */
    SeqCodeTablePtr sctp;        /* for getting symbols */
    Uint1 newcode,               /* requested output code */
          oldcode;               /* current input seq code (0 if not raw) */
    Uint1 byte,                    /* current byte in buf */
        bc,                      /* value to start bitctr */
        bitctr,                  /* current shift */
        lshift,                  /* amount to left shift on decompact */
        rshift,                  /* amount to right shift residue value */
        mask;                    /* mask for compact byte */
    struct seqport PNTR curr ,   /* current active seqport if seg or ref */
        PNTR segs,                    /* segments if seg or ref */
        PNTR next;                    /* if part of a segment chain */
	SPCachePtr cache;
	SPCacheQPtr cacheq;          /* used instead of cache for ncbi2na or ncbi4na to iupacna fasta lookup */
	ByteStorePtr bp;             /* used by both raw and delta seq pieces */
} SeqPort, PNTR SeqPortPtr;

/*****************************************************************************
*
*  Structure used in SeqPort DNA Compression
*
*****************************************************************************/
typedef struct SPCompress {
  Uint1Ptr buffer;  /* Buffer with 2na DNA sequence */
  Int4 type;        /* Type of stored sequence */
  Int4 residues;    /* Number of residues in buffer */
  Int4 used;        /* Number of bytes used in buffer */
  Int4 allocated;   /* Number of bytes allocated in buffer */ 
  Uint4Ptr lbytes;  /* Ambiguity information */
} SPCompress, PNTR SPCompressPtr;


NLM_EXTERN SeqPortPtr SeqPortNew PROTO((BioseqPtr bsp, Int4 start, Int4 stop, Uint1 strand, Uint1 code));
NLM_EXTERN SeqPortPtr SeqPortNewByLoc PROTO((SeqLocPtr seqloc, Uint1 code));
NLM_EXTERN SeqPortPtr SeqPortFree PROTO((SeqPortPtr spp));
NLM_EXTERN Int4 SeqPortTell PROTO((SeqPortPtr spp));
NLM_EXTERN Int2 SeqPortSeek PROTO((SeqPortPtr spp, Int4 offset, Int2 origin));
NLM_EXTERN Int4 SeqPortLen PROTO((SeqPortPtr spp));
NLM_EXTERN Uint1 LIBCALL SeqPortGetResidue PROTO((SeqPortPtr spp));
NLM_EXTERN Int2 LIBCALL SeqPortRead PROTO((SeqPortPtr spp, BytePtr buf, Int2 len));
NLM_EXTERN Uint1 GetGapCode PROTO((Uint1 seqcode));
NLM_EXTERN Boolean LIBCALL SeqPortSetUpFields PROTO((SeqPortPtr spp, Int4 start, Int4 stop, Uint1 strand, Uint1 newcode));
NLM_EXTERN Boolean LIBCALL SeqPortSetUpAlphabet PROTO((SeqPortPtr spp, Uint1 curr_code, Uint1 newcode));

/*******************************************************************************
*	
*   SeqPortStream (bsp, flags, userdata, proc)
*   SeqPortStreamInt (bsp, start, stop, strand, flags, userdata, proc)
*   SeqPortStreamLoc (slp, flags, userdata, proc)
*   SeqPortStreamLit (lit, flags, userdata, proc)
*       Efficient functions to stream through sequence
*
********************************************************************************/

typedef void (LIBCALLBACK *SeqPortStreamProc) (
  CharPtr sequence,
  Pointer userdata
);

typedef unsigned long StreamFlgType;

#define STREAM_EXPAND_GAPS      1
#define GAP_TO_SINGLE_DASH      2
#define EXPAND_GAPS_TO_DASHES   3

#define KNOWN_GAP_AS_PLUS       4

#define SUPPRESS_VIRT_SEQ       8
#define STREAM_VIRT_AS_PLUS    16

#define STREAM_CORRECT_INVAL   64

#define STREAM_ALLOW_NEG_GIS  128 /* for internal use only by NCBI ID group */

#define STREAM_HTML_SPANS     256 /* show span tags at begining of each line */

#define STREAM_ALL_FASTA_IDS  512 /* in FASTA streamer, show all Seq-ids */

NLM_EXTERN Int4 SeqPortStream (
  BioseqPtr bsp,
  StreamFlgType flags,
  Pointer userdata,
  SeqPortStreamProc proc
);

NLM_EXTERN Int4 SeqPortStreamInt (
  BioseqPtr bsp,
  Int4 start,
  Int4 stop,
  Uint1 strand,
  StreamFlgType flags,
  Pointer userdata,
  SeqPortStreamProc proc
);

NLM_EXTERN Int4 SeqPortStreamLoc (
  SeqLocPtr slp,
  StreamFlgType flags,
  Pointer userdata,
  SeqPortStreamProc proc
);

NLM_EXTERN Int4 SeqPortStreamLit (
  SeqLitPtr lit,
  StreamFlgType flags,
  Pointer userdata,
  SeqPortStreamProc proc
);

/*******************************************************************************
*	
*   StreamCacheSetup (bsp, slp, flags, scp)
*   StreamCacheGetResidue (scp)
*   StreamCacheSetPosition (scp, pos)
*       SeqPort functional replacement implemented on top of SeqPortStreams
*
********************************************************************************/

typedef struct streamcache {
  BioseqPtr     bsp;
  SeqLocPtr     slp;
  Char          buf [4004];
  Int2          ctr;
  Int2          total;
  Int4          offset;
  Int4          length;
  StreamFlgType flags;
  Boolean       failed;
} StreamCache, PNTR StreamCachePtr;

NLM_EXTERN Boolean StreamCacheSetup (
  BioseqPtr bsp,
  SeqLocPtr slp,
  StreamFlgType flags,
  StreamCache PNTR scp
);

NLM_EXTERN Uint1 StreamCacheGetResidue (
  StreamCache PNTR scp
);

NLM_EXTERN Boolean StreamCacheSetPosition (
  StreamCache PNTR scp,
  Int4 pos
);

/*
the following functions are for quick alphabet expansion, and require buffers
allocated with 4-byte or 2-byte alignment, because they cast 2 or 4 bytes into
Uint2 or Uint4 for fast integer copying.
*/

NLM_EXTERN Uint4Ptr LIBCALL MapNa2ByteToIUPACString PROTO((Uint1Ptr bytep, Uint4Ptr buf, Int4 total));
NLM_EXTERN Uint2Ptr LIBCALL MapNa4ByteToIUPACString PROTO((Uint1Ptr bytep, Uint2Ptr buf, Int4 total));
NLM_EXTERN Uint2Ptr LIBCALL MapNa4ByteToIUPACplusGapString PROTO((Uint1Ptr bytep, Uint2Ptr buf, Int4 total));
NLM_EXTERN Uint2Ptr LIBCALL MapNa2ByteToNa4String PROTO((Uint1Ptr bytep, Uint2Ptr buf, Int4 total));
NLM_EXTERN Uint4Ptr LIBCALL MapNa2ByteTo4BitString PROTO((Uint1Ptr bytep, Uint4Ptr buf, Int4 total));
NLM_EXTERN Uint2Ptr LIBCALL MapNa4ByteTo4BitString PROTO((Uint1Ptr bytep, Uint2Ptr buf, Int4 total));


/*****************************************************************************
*
*   SeqPortSetValues(spp)
*      Copies the values is_circle, is_seg, and do_virtual from spp to
*        any dependent SeqPortPtrs it contains. This is necessary for segmented
*        reference, or delta types of Bioseqs and on SeqPortNewByLoc()
*
*      SeqPortSet_... functions call this function
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqPortSetValues (SeqPortPtr spp);

/*****************************************************************************
*
*   SeqPortSet_is_circle(spp, value)
*      if (value) is TRUE, then SeqPort will go around the ends of a circular
*        molecule without stopping.
*      Default is FALSE
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqPortSet_is_circle (SeqPortPtr spp, Boolean value);

/*****************************************************************************
*
*   SeqPortSet_is_seg(spp, value)
*      if (value) is TRUE, then SeqPort will return SEQPORT_EOS whenever it
*         crosses a segment boundary in the SeqPort. When is_seg = TRUE,
*         SEQPORT_VIRT will NOT be returned on virtual or NULL segments. Instead
*         SEQPORT_EOS will be received only as the SeqPort passes over those
*         segments.
*      Default is FALSE, SeqPort will NEVER return SEQPORT_EOS
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqPortSet_is_seg (SeqPortPtr spp, Boolean value);

/*****************************************************************************
*
*   SeqPortSet_do_virtual(spp, value)
*      if (value) is TRUE, then SeqPort will return the appropriate ambiguity
*         character (e.g. "N" or "X") for the length of a virtual Bioseq or delta
*         gap segment. It will still return SEQPORT_VIRT for a "NULL" segment
*         (ie. gap of unknown length).
*         
*      Default is FALSE. In this case SeqPort will return a single SEQPORT_VIRT
*         when encountering a virtual Bioseq, just as for a "NULL" segment.
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqPortSet_do_virtual (SeqPortPtr spp, Boolean value);
NLM_EXTERN Boolean LIBCALL SeqPortSet_do_virtualEx (SeqPortPtr spp, Boolean value, Boolean gapIsZero);

/*****************************************************************************
*
*   BioseqHash(bsp)
*   	Computes a (almost) unique hash code for a bioseq
*
*****************************************************************************/
NLM_EXTERN Uint4 BioseqHash PROTO((BioseqPtr bsp));

/*****************************************************************************
*
*   ProteinFromCdRegion(sfp, include_stop)
*   	produces a ByteStorePtr containing the protein sequence in
*   ncbieaa code for the CdRegion sfp.  If include_stop, will translate
*   through stop codons.  If NOT include_stop, will stop at first stop
*   codon and return the protein sequence NOT including the terminating
*   stop.  Supports reading frame, alternate genetic codes, and code breaks
*   in the CdRegion. Removes trailing "X" on partial translation.
*
* if no explict partial at either end, but feature is
* annotated as partial, then guess should use internal
* amino acid code
*
*****************************************************************************/
NLM_EXTERN ByteStorePtr ProteinFromCdRegion PROTO(( SeqFeatPtr sfp, Boolean include_stop));
NLM_EXTERN ByteStorePtr ProteinFromCdRegionEx PROTO((SeqFeatPtr sfp, Boolean include_stop, Boolean remove_trailingX));
NLM_EXTERN ByteStorePtr ProteinFromCdRegionExEx PROTO((SeqFeatPtr sfp, Boolean include_stop, Boolean remove_trailingX, BoolPtr altStartP, Boolean farProdFetchOK));
NLM_EXTERN ByteStorePtr ProteinFromCdRegionExWithTrailingCodonHandling PROTO((SeqFeatPtr sfp, Boolean include_stop, Boolean remove_trailingX, Boolean no_stop_at_end_of_complete_cds));

/*****************************************************************************
*
*   Uint1 AAForCodon (Uint1Ptr codon, CharPtr codes)
*   	codon is 3 values in ncbi4na code
*       codes is the geneic code array to use
*          MUST have 'X' as unknown amino acid
*
*****************************************************************************/
NLM_EXTERN Uint1 AAForCodon PROTO((Uint1Ptr codon, CharPtr codes));

/*****************************************************************************
*
*   Uint1 IndexForCodon (codon, code)
*   	returns index into genetic codes codon array, give 3 bases of the
*       codon in any alphabet
*       returns INVALID_RESIDUE on failure
*   
*****************************************************************************/
NLM_EXTERN Uint1 IndexForCodon PROTO((Uint1Ptr codon, Uint1 code));

/*****************************************************************************
*
*   Boolean CodonForIndex (index, code, codon)
*   	Fills codon (3 Uint1 array) with codon corresponding to index,
*       in sequence alphabet code.
*       Index is the Genetic code index.
*       returns TRUE on success.
*
*****************************************************************************/
NLM_EXTERN Boolean CodonForIndex PROTO((Uint1 index, Uint1 code, Uint1Ptr codon));

/*****************************************************************************
*
*   Int2 GetFrameFromLoc (slp)
*   	returns 1,2,3 if can find the frame
*   	0 if not
*
*****************************************************************************/
NLM_EXTERN Uint1 GetFrameFromLoc PROTO((SeqLocPtr slp));

/******************************************************************
*
*	dnaLoc_to_aaLoc(sfp, dna_loc, merge, frame, allowTerminator)
*	map a SeqLoc on the DNA sequence
*       to a Seq-loc in the	protein sequence
*       through a CdRegion feature
*   if (merge) adjacent intervals on the amino acid sequence
*      are merged into one. This should be the usual case.
*   We try to report the frame if the caller provides a suitable pointer
*   If allowTerminator, can map the termination codon as a legal location
*
******************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL dnaLoc_to_aaLoc(SeqFeatPtr sfp, SeqLocPtr dna_loc, Boolean merge, Int4Ptr frame, Boolean allowTerminator);

/******************************************************************
*
*	productLoc_to_locationLoc(sfp, product_loc)
*	map a SeqLoc on the product sequence
*       to a Seq-loc in the location sequence
*       through a feature
*
*       this more general function is now called by
*          aaLoc_to_dnaLoc()
*
******************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL productLoc_to_locationLoc(SeqFeatPtr sfp, SeqLocPtr product_loc);

/******************************************************************
*
*	aaLoc_to_dnaLoc(sfp, aa_loc)
*	map a SeqLoc on the amino acid sequence
*       to a Seq-loc in the	DNA sequence
*       through a CdRegion feature
*
******************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL aaLoc_to_dnaLoc(SeqFeatPtr sfp, SeqLocPtr aa_loc);

/******************************************************************
*
*	aaFeatLoc_to_dnaFeatLoc(sfp, aa_loc)
*	map a SeqLoc on the amino acid sequence
*       to a Seq-loc in the	DNA sequence
*       through a CdRegion feature
*
*       uses aaLoc_to_dnaLoc() but does additional checks to
*       extend dnaLoc at either end to compensate for positions in
*       the dna which do not corresspond to the amino acid sequence
*       (partial codons which are not translated).
*
******************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL aaFeatLoc_to_dnaFeatLoc(SeqFeatPtr sfp, SeqLocPtr aa_loc);

/******************************************************************
*
*	productInterval_to_locationIntervals(sfp, aa_start, aa_stop)
*	map the amino acid sequence to a chain of Seq-locs in the 
*	DNA sequence through a CdRegion feature
*
******************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL productInterval_to_locationIntervals (SeqFeatPtr sfp, Int4 aa_start, Int4 aa_stop, Boolean aa_partialn);

/*-------------- BioseqRevComp () ---------------------------*/
/***********************************************************************
*   BioseqRevComp:   Takes the nucleic acid sequence from Bioseq
*	Entry and gives the reverse complement sequence in place
*       Does not change features.
************************************************************************/
NLM_EXTERN Boolean LIBCALL BioseqRevComp (BioseqPtr bsp);


/*-------------- BioseqComplement () ---------------------------*/
/***********************************************************************
*   BioseqComplement:   Takes the nucleic acid sequence from Bioseq
*	Entry and gives the complement sequence in place
*       Does not change features.
************************************************************************/
NLM_EXTERN Boolean LIBCALL BioseqComplement (BioseqPtr bsp);

           
/*-------------- BioseqReverse () ---------------------------*/
/***********************************************************************
*   BioseqReverse:   Takes nucleic acid sequence from Bioseq Entry and 
*	reverses the whole sequence in place
*       Does not change features.
************************************************************************/
NLM_EXTERN Boolean LIBCALL BioseqReverse (BioseqPtr bsp);


/*-------------- ContigRevComp () ---------------------------*/
/***********************************************************************
*   ContigRevComp:   Reverse complement segmented or delta bioseq
************************************************************************/
NLM_EXTERN Boolean LIBCALL ContigRevComp (BioseqPtr bsp);

/*****************************************************************************
*
*  SPCompressNew(void); - allocated memory for SPCompress structure
*
*****************************************************************************/
NLM_EXTERN SPCompressPtr SPCompressNew(void);

/*****************************************************************************
*
*  SPCompressFree(SPCompressPtr spc); -  free SPCompress structure
*
*****************************************************************************/
NLM_EXTERN void SPCompressFree(SPCompressPtr spc);

/*****************************************************************************
*
*   SPCompressDNA(SeqPortPtr spp);
*       converts a ncbi4na taken from spp into ncbi2na
*       buffer stored inside SPCompress structue together
*       with ambiguity information
*       returns pointer SPCompress structure or NULL if error
*
*       NOTE: In this function we do not know - what is length
*             of sequence to compress. Terminated flag for this 
*             function is SEQPORT_EOF returned from spp.
*
*****************************************************************************/
NLM_EXTERN SPCompressPtr SPCompressDNA(SeqPortPtr spp);

/*****************************************************************************
*
*   SPRebuildDNA(SPCompressPtr spc);
*       translates spc ncbi2na encoding buffer into
*       spc ncbi4na encoding buffer with rebuild ambiguities
*
*       spc - must be valid SPCompress structure returned
*       from SPCompressDNA() function in ncbi2na encoding
*
*****************************************************************************/
NLM_EXTERN Boolean SPRebuildDNA(SPCompressPtr spc);

/*****************************************************************************
*
*   ComposeCodonsRecognizedString (trna, buf, buflen);
*       Copies codon recognized string to buf, returns number of codons
*
*****************************************************************************/

NLM_EXTERN Int2 ComposeCodonsRecognizedString (tRNAPtr trna, CharPtr buf, size_t buflen);

/*****************************************************************************
*
*   TransTableNew (Int2 genCode);
*       Initializes TransTable finite state machine for 6-frame translation
*       and open reading frame search, allowing nucleotide ambiguity characters
*
*****************************************************************************/

typedef struct fsatranstable {
  Int2     genCode;
  Char     ncbieaa [65];
  Char     sncbieaa [65];
  Uint2    nextBase [3376];
  Char     aminoAcid [3376] [2];
  Char     orfStart [3376] [2];
  Uint1    basesToIdx [256];
} TransTable, PNTR TransTablePtr;

/* allocate 6-frame finite state translation table and initialize with indicated genetic code */
NLM_EXTERN TransTablePtr TransTableNew (Int2 genCode);
NLM_EXTERN TransTablePtr TransTableFree (TransTablePtr tbl);
NLM_EXTERN void TransTableFreeAll (void);

#define TTBL_TOP_STRAND  0
#define TTBL_BOT_STRAND  1

#define TTBL_ATG_STATE  48
#define TTBL_CAT_STATE 229

/* macros for using finite state machine for 6-frame translation */
#define NextCodonState(tbl,cur,ch) (tbl->nextBase [(int) (Uint2) cur] + tbl->basesToIdx [(int) (Uint1) ch])
#define GetCodonResidue(tbl,cur,stnd) (tbl->aminoAcid [(int) (Uint2) cur] [stnd])
#define GetStartResidue(tbl,cur,stnd) (tbl->orfStart [(int) (Uint2) cur] [stnd])
#define IsOrfStart(tbl,cur,stnd) ((Boolean) (GetStartResidue(tbl,cur,stnd) == 'M'))
#define IsAmbigStart(tbl,cur,stnd) ((Boolean) (GetStartResidue(tbl,cur,stnd) == 'X'))
#define IsAnyStart(tbl,cur,stnd) ((Boolean) (GetStartResidue(tbl,cur,stnd) != '-'))
#define IsOrfStop(tbl,cur,stnd) ((Boolean) (GetCodonResidue(tbl,cur,stnd) == '*'))
#define IsATGStart(tbl,cur,stnd) ((Boolean) (IsOrfStart(tbl,cur,stnd) && (stnd ? (cur == TTBL_CAT_STATE) : (cur == TTBL_ATG_STATE))))
#define IsAltStart(tbl,cur,stnd) ((Boolean) (IsOrfStart(tbl,cur,stnd) && (stnd ? (cur != TTBL_CAT_STATE) : (cur != TTBL_ATG_STATE))))

typedef void (LIBCALLBACK *TransTableMatchProc) (Int4 position, Char residue, Boolean atgStart, Boolean altStart, Boolean orfStop, Int2 frame, Uint1 strand, Pointer userdata);

/* convenience function calls user callback for each strand of entire bioseq */

NLM_EXTERN void TransTableProcessBioseq (
  TransTablePtr tbl,
  TransTableMatchProc matchProc,
  Pointer userdata,
  BioseqPtr bsp
);

/* trans table translation functions can be passed cds feature or individual parameters */

NLM_EXTERN ByteStorePtr TransTableTranslateCdRegion (
  TransTablePtr  PNTR tblptr,
  SeqFeatPtr cds,
  Boolean include_stop,
  Boolean remove_trailingX,
  Boolean no_stop_at_end_of_complete_cds
);

NLM_EXTERN ByteStorePtr TransTableTranslateCdRegionEx (
  TransTablePtr  PNTR tblptr,
  SeqFeatPtr cds,
  Boolean include_stop,
  Boolean remove_trailingX,
  Boolean no_stop_at_end_of_complete_cds,
  BoolPtr altStartP,
  Boolean farProdFetchOK
);

NLM_EXTERN ByteStorePtr TransTableTranslateSeqLoc (
  TransTablePtr  PNTR tblptr,
  SeqLocPtr location,
  Int2 genCode,
  Uint1 frame,
  Boolean include_stop,
  Boolean remove_trailingX
);

/* returns string of bases to translate */

NLM_EXTERN CharPtr ReadCodingRegionBases (
  SeqLocPtr location,
  Int4 len,
  Uint1 frame,
  Int4Ptr totalP
);

/* allow reuse of translation tables by saving as AppProperty, avoids unnecessary initializations */

NLM_EXTERN TransTablePtr PersistentTransTableByGenCode (
  Int2 genCode
);

NLM_EXTERN TransTablePtr PersistentTransTableByCdRegion (
  SeqFeatPtr cds
);

NLM_EXTERN ValNodePtr MakeCodeBreakList (
  SeqLocPtr cdslocation,
  Int4 len,
  CodeBreakPtr cbp,
  Uint1 frame
);

/*****************************************************************************
*
*   SeqSearch
*       Initializes SeqSearch finite state machine for sequence searching
*       Based on Practical Algorithms for Programmers by Binstock and Rex
*
*****************************************************************************/

struct SeqSearch;
typedef struct SeqSearch* SeqSearchPtr;

typedef void (LIBCALLBACK *SeqSearchMatchProc) (Int4 position, CharPtr name, CharPtr pattern, Int2 cutSite, Uint1 strand, Pointer userdata);

/* create empty nucleotide sequence search finite state machine */

NLM_EXTERN SeqSearchPtr SeqSearchNew (
  SeqSearchMatchProc matchproc,
  Pointer userdata
);

/*
   add nucleotide pattern or restriction site to sequence search finite state
   machine, expands using ambiguity codes R = A and G, H = A, C and T, etc.
*/

typedef unsigned long SearchFlgType;

#define SEQ_SEARCH_JUST_TOP_STRAND  1
#define SEQ_SEARCH_EXPAND_PATTERN   2
#define SEQ_SEARCH_ALLOW_MISMATCH   4

NLM_EXTERN void SeqSearchAddNucleotidePattern (
  SeqSearchPtr tbl,
  CharPtr name,
  CharPtr pattern,
  Int2 cutSite,
  SearchFlgType flags
);

/* program passes each character in turn to finite state machine */

NLM_EXTERN void SeqSearchProcessCharacter (
  SeqSearchPtr tbl,
  Char ch
);

/* convenience function calls SeqSearchProcessCharacter for entire bioseq */

NLM_EXTERN void SeqSearchProcessBioseq (
  SeqSearchPtr tbl,
  BioseqPtr bsp
);

/* reset state and position to allow another run with same search patterns */

NLM_EXTERN void SeqSearchReset (
  SeqSearchPtr tbl
);

/* clean up sequence search finite state machine allocated memory */

NLM_EXTERN SeqSearchPtr SeqSearchFree (
  SeqSearchPtr tbl
);


/*****************************************************************************
*
*   ProtSearch
*       Initializes ProtSearch finite state machine for sequence searching
*       Based on Practical Algorithms for Programmers by Binstock and Rex
*
*****************************************************************************/

struct ProtSearch;
typedef struct ProtSearch* ProtSearchPtr;

typedef void (LIBCALLBACK *ProtSearchMatchProc) (Int4 position, CharPtr name, CharPtr pattern, Pointer userdata);

/* create empty protein sequence search finite state machine */

NLM_EXTERN ProtSearchPtr ProtSearchNew (
  ProtSearchMatchProc matchproc,
  Pointer userdata
);

/*
   add protein pattern to protein sequence search finite state machine,
   expands using ambiguity codes B = D and N, Z = E and Q, etc.
*/

NLM_EXTERN void ProtSearchAddProteinPattern (
  ProtSearchPtr tbl,
  CharPtr name,
  CharPtr pattern,
  SearchFlgType flags
);

/* program passes each character in turn to finite state machine */

NLM_EXTERN void ProtSearchProcessCharacter (
  ProtSearchPtr tbl,
  Char ch
);

/* convenience function calls ProtSearchProcessCharacter for entire bioseq */

NLM_EXTERN void ProtSearchProcessBioseq (
  ProtSearchPtr tbl,
  BioseqPtr bsp
);


/* reset state and position to allow another run with same search patterns */

NLM_EXTERN void ProtSearchReset (
  ProtSearchPtr tbl
);

/* clean up sequence search finite state machine allocated memory */

NLM_EXTERN ProtSearchPtr ProtSearchFree (
  ProtSearchPtr tbl
);


/*****************************************************************************
*
*  Convenience functions for genome processing use BioseqLockById to get sequence
*  record (perhaps with phrap quality score graphs) so fetching from some network
*  or local server must be enabled, or sequences must already be in memory.
*
*****************************************************************************/

NLM_EXTERN CharPtr GetSequenceByBsp (
  BioseqPtr bsp
);

NLM_EXTERN CharPtr GetSequenceByIdOrAccnDotVer (
  SeqIdPtr sip,
  CharPtr accession,
  Boolean is_na
);

NLM_EXTERN CharPtr GetSequenceByFeature (
  SeqFeatPtr sfp
);

NLM_EXTERN CharPtr GetDNAbyAccessionDotVersion (
  CharPtr accession
);

NLM_EXTERN BytePtr GetScoresbyAccessionDotVersion (
  CharPtr accession,
  Int4Ptr bsplength
);

NLM_EXTERN BytePtr GetScoresbySeqId (
  SeqIdPtr sip,
  Int4Ptr bsplength
);

/*****************************************************************************
*
*   ConvertNsToGaps
*       Assumes string of Ns means a gap of known length
*
*****************************************************************************/

NLM_EXTERN void ConvertNsToGaps (
  BioseqPtr bsp,
  Pointer userdata
);

/**************************************************************
*
*  Returns a protein molecular weight for a SeqLoc
*    If it cannot calculate the value it returns -1.0
*    If sequence contains X, J, or O it fails
*
***************************************************************/
NLM_EXTERN FloatHi MolWtForLoc (SeqLocPtr slp);

NLM_EXTERN FloatHi MolWtForBsp (BioseqPtr bsp);

NLM_EXTERN FloatHi MolWtForStr (CharPtr str);



NLM_EXTERN Boolean LIBCALL ReverseSeqData (Uint1 seqtype, Int4 seqlen, SeqDataPtr sdp);
NLM_EXTERN Boolean ComplementSeqData (Uint1 seqtype, Int4 seqlen, SeqDataPtr sdp);


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
