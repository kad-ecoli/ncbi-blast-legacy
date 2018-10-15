/*  seqport.c
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
* File Name:  seqport.c
*
* Author:  James Ostell
*   
* Version Creation Date: 7/13/91
*
* $Revision: 6.186 $
*
* File Description:  Ports onto Bioseqs
*
* Modifications:  
* --------------------------------------------------------------------------
* Date       Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*/

/** for ErrPostEx() ****/

static char *this_module = "ncbiapi";
#define THIS_MODULE this_module
static char *this_file = __FILE__;
#define THIS_FILE this_file

/**********************/


#include <seqport.h>
#include <edutil.h>    /* for SeqLoc creation functions */
#include <gather.h>    /* for SeqLocOffset function */
#include <sqnutils.h>
#include <explore.h>   /* for BioseqFindFromSeqLoc function */
#include <subutil.h>
#include <tofasta.h>   /* for FastaSeqLineEx function */
#include <salutil.h> 
#include <alignmgr2.h> /* for correcting alignments when converting to delta */


NLM_EXTERN Boolean LIBCALL SeqPortAdjustLength (SeqPortPtr spp);

/*****************************************************************************
*
*   Fast mapping arrays
*
*****************************************************************************/

static Uint1Ptr Na2toIUPAC = NULL;
static Uint1Ptr Na4toIUPAC = NULL;
static Uint1Ptr Na4toIUPACplusGap = NULL;
static Uint1Ptr Na2toNa4 = NULL;
static Uint1Ptr Na2to4Bit = NULL;
static Uint1Ptr Na4to4Bit = NULL;
static TNlmMutex seqport_mutex = NULL;


/*****************************************************************************
*
*   MapNa2ByteToIUPACString and MapNa4ByteToIUPACString now copy directly to
*     the expanded character buffer for efficiency
*
*****************************************************************************/

static void InitNa2toIUPAC (void)

{
  Int2  base [4], index, j;
  Char  convert [4] = {'A', 'C', 'G', 'T'};
  Int4  ret;
  Uint1Ptr Na2toIUPAC_local = NULL;

  ret = NlmMutexLockEx (&seqport_mutex);  /* protect this section */
  if (ret) {
    ErrPostEx (SEV_FATAL, 0, 0, "MapNa2ByteToIUPACString mutex failed [%ld]", (long) ret);
    return;
  }

  if (Na2toIUPAC == NULL) {
    Na2toIUPAC_local = MemNew (sizeof (Uint1) * 1024);

    if (Na2toIUPAC_local != NULL) {
      for (base [0] = 0; base [0] < 4; (base [0])++) {
        for (base [1] = 0; base [1] < 4; (base [1])++) {
          for (base [2] = 0; base [2] < 4; (base [2])++) {
            for (base [3] = 0; base [3] < 4; (base [3])++) {
              index = 4 * (base [0] * 64 + base [1] * 16 + base [2] * 4 + base [3]);
              for (j = 0; j < 4; j++) {
                Na2toIUPAC_local [index + j] = convert [(base [j])];
              }
            }
          }
        }
      }
    }
    Na2toIUPAC = Na2toIUPAC_local;
  }

  NlmMutexUnlock (seqport_mutex);
}

NLM_EXTERN Uint4Ptr LIBCALL MapNa2ByteToIUPACString (Uint1Ptr bytep, Uint4Ptr buf, Int4 total)

{
  Uint4Ptr  bp;
  Uint1     byte;
  Int2      index;
  Int4      k;
  Uint4Ptr  ptr;

  if (bytep == NULL || buf == NULL) return buf;
  ptr = buf;

  /* initialize array if not yet set (first time function is called) */

  if (Na2toIUPAC == NULL) {
    InitNa2toIUPAC ();
  }

  if (Na2toIUPAC == NULL) return buf;

  /* now return 4 character string for each compressed byte */

  for (k = 0; k < total; k++) {
    byte = *bytep;
    bytep++;
    index = 4 * byte;
    bp = (Uint4Ptr) (Na2toIUPAC + index);
    /* copy 4 bytes at a time */
    /*
    for (j = 0; j < 4; j++) {
      *ptr = *bp;
      ptr++;
      bp++;
    }
    */
    *ptr = *bp;
    ptr++;
  }

  return ptr;
}

static void InitNa4toIUPAC (void)

{
  Int2  base [2], index, j;
  Char  convert [16] = {'N', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
                        'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
  Int4  ret;
  Uint1Ptr Na4toIUPAC_local = NULL;

  ret = NlmMutexLockEx (&seqport_mutex);  /* protect this section */
  if (ret) {
    ErrPostEx (SEV_FATAL, 0, 0, "MapNa4ByteToIUPACString mutex failed [%ld]", (long) ret);
    return;
  }

  if (Na4toIUPAC == NULL) {
    Na4toIUPAC_local = MemNew (sizeof (Uint1) * 512);

    if (Na4toIUPAC_local != NULL) {
      for (base [0] = 0; base [0] < 16; (base [0])++) {
        for (base [1] = 0; base [1] < 16; (base [1])++) {
          index = 2 * (base [0] * 16 + base [1]);
          for (j = 0; j < 2; j++) {
            Na4toIUPAC_local [index + j] = convert [(base [j])];
          }
        }
      }
    }
    Na4toIUPAC = Na4toIUPAC_local;
  }

  NlmMutexUnlock (seqport_mutex);
}

NLM_EXTERN Uint2Ptr LIBCALL MapNa4ByteToIUPACString (Uint1Ptr bytep, Uint2Ptr buf, Int4 total)

{
  Uint2Ptr  bp;
  Uint1     byte;
  Int2      index;
  Int4      k;
  Uint2Ptr  ptr;

  if (bytep == NULL || buf == NULL) return buf;
  ptr = buf;

  /* initialize array if not yet set (first time function is called) */

  if (Na4toIUPAC == NULL) {
    InitNa4toIUPAC ();
  }

  if (Na4toIUPAC == NULL) return buf;

  /* now return 2 character string for each compressed byte */

  for (k = 0; k < total; k++) {
    byte = *bytep;
    bytep++;
    index = 2 * byte;
    bp = (Uint2Ptr)  (Na4toIUPAC + index);
    /* copy 2 bytes at a time */
    /*
    for (j = 0; j < 2; j++) {
      *ptr = *bp;
      ptr++;
      bp++;
    }
    */
    *ptr = *bp;
    ptr++;
  }

  return ptr;
}

static void InitNa4toIUPACplusGap (void)

{
  Int2  base [2], index, j;
  Char  convert [16] = {'-', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
                        'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
  Int4  ret;
  Uint1Ptr Na4toIUPACplusGap_local = NULL;

  ret = NlmMutexLockEx (&seqport_mutex);  /* protect this section */
  if (ret) {
    ErrPostEx (SEV_FATAL, 0, 0, "MapNa4ByteToIUPACplusGapString mutex failed [%ld]", (long) ret);
    return;
  }

  if (Na4toIUPACplusGap == NULL) {
    Na4toIUPACplusGap_local = MemNew (sizeof (Uint1) * 512);

    if (Na4toIUPACplusGap_local != NULL) {
      for (base [0] = 0; base [0] < 16; (base [0])++) {
        for (base [1] = 0; base [1] < 16; (base [1])++) {
          index = 2 * (base [0] * 16 + base [1]);
          for (j = 0; j < 2; j++) {
            Na4toIUPACplusGap_local [index + j] = convert [(base [j])];
          }
        }
      }
    }
    Na4toIUPACplusGap = Na4toIUPACplusGap_local;
  }

  NlmMutexUnlock (seqport_mutex);
}

NLM_EXTERN Uint2Ptr LIBCALL MapNa4ByteToIUPACplusGapString (Uint1Ptr bytep, Uint2Ptr buf, Int4 total)

{
  Uint2Ptr  bp;
  Uint1     byte;
  Int2      index;
  Int4      k;
  Uint2Ptr  ptr;

  if (bytep == NULL || buf == NULL) return buf;
  ptr = buf;

  /* initialize array if not yet set (first time function is called) */

  if (Na4toIUPACplusGap == NULL) {
    InitNa4toIUPACplusGap ();
  }

  if (Na4toIUPACplusGap == NULL) return buf;

  /* now return 2 character string for each compressed byte */

  for (k = 0; k < total; k++) {
    byte = *bytep;
    bytep++;
    index = 2 * byte;
    bp = (Uint2Ptr)  (Na4toIUPACplusGap + index);
    /* copy 2 bytes at a time */
    /*
    for (j = 0; j < 2; j++) {
      *ptr = *bp;
      ptr++;
      bp++;
    }
    */
    *ptr = *bp;
    ptr++;
  }

  return ptr;
}

static void InitNa2toNa4 (void)

{
  Int2   pair [2], index, j;
  Uint1  convert [16] = {17,  18,  20,  24,  33,  34,  36,  40,
                         65,  66,  68,  72, 129, 130, 132, 136};
  Int4   ret;
  Uint1Ptr Na2toNa4_local = NULL;

  ret = NlmMutexLockEx (&seqport_mutex);  /* protect this section */
  if (ret) {
    ErrPostEx (SEV_FATAL, 0, 0, "MapNa2ByteToIUPACString mutex failed [%ld]", (long) ret);
    return;
  }

  if (Na2toNa4 == NULL) {
    Na2toNa4_local = MemNew (sizeof (Uint1) * 512);

    if (Na2toNa4_local != NULL) {
      for (pair [0] = 0; pair [0] < 16; (pair [0])++) {
        for (pair [1] = 0; pair [1] < 16; (pair [1])++) {
          index = 2 * (pair [0] * 16 + pair [1]);
          for (j = 0; j < 2; j++) {
            Na2toNa4_local [index + j] = convert [(pair [j])];
          }
        }
      }
    }
    Na2toNa4 = Na2toNa4_local;
  }

  NlmMutexUnlock (seqport_mutex);
}

NLM_EXTERN Uint2Ptr LIBCALL MapNa2ByteToNa4String (Uint1Ptr bytep, Uint2Ptr buf, Int4 total)

{
  Uint2Ptr  bp;
  Uint1     byte;
  Int2      index;
  Int4      k;
  Uint2Ptr  ptr;

  if (bytep == NULL || buf == NULL) return buf;
  ptr = buf;

  /* initialize array if not yet set (first time function is called) */

  if (Na2toNa4 == NULL) {
    InitNa2toNa4 ();
  }

  if (Na2toNa4 == NULL) return buf;

  /* now return 2 character byte for each compressed byte */

  for (k = 0; k < total; k++) {
    byte = *bytep;
    bytep++;
    index = 2 * byte;
    bp = (Uint2Ptr)  (Na2toNa4 + index);
    /* copy 2 bytes at a time */
    /*
    for (j = 0; j < 2; j++) {
      *ptr = *bp;
      ptr++;
      bp++;
    }
    */
    *ptr = *bp;
    ptr++;
  }

  return ptr;
}

static void InitNa2to4Bit (void)

{
  Int2  base [4], index, j;
  Uint1  convert [4] = {1, 2, 4, 8};
  Int4  ret;
  Uint1Ptr Na2to4Bit_local = NULL;

  ret = NlmMutexLockEx (&seqport_mutex);  /* protect this section */
  if (ret) {
    ErrPostEx (SEV_FATAL, 0, 0, "MapNa2ByteTo4BitString mutex failed [%ld]", (long) ret);
    return;
  }

  if (Na2to4Bit == NULL) {
    Na2to4Bit_local = MemNew (sizeof (Uint1) * 1024);

    if (Na2to4Bit_local != NULL) {
      for (base [0] = 0; base [0] < 4; (base [0])++) {
        for (base [1] = 0; base [1] < 4; (base [1])++) {
          for (base [2] = 0; base [2] < 4; (base [2])++) {
            for (base [3] = 0; base [3] < 4; (base [3])++) {
              index = 4 * (base [0] * 64 + base [1] * 16 + base [2] * 4 + base [3]);
              for (j = 0; j < 4; j++) {
                Na2to4Bit_local [index + j] = convert [(base [j])];
              }
            }
          }
        }
      }
    }
    Na2to4Bit = Na2to4Bit_local;
  }

  NlmMutexUnlock (seqport_mutex);
}

NLM_EXTERN Uint4Ptr LIBCALL MapNa2ByteTo4BitString (Uint1Ptr bytep, Uint4Ptr buf, Int4 total)

{
  Uint4Ptr  bp;
  Uint1     byte;
  Int2      index;
  Int4      k;
  Uint4Ptr  ptr;

  if (bytep == NULL || buf == NULL) return buf;
  ptr = buf;

  /* initialize array if not yet set (first time function is called) */

  if (Na2to4Bit == NULL) {
    InitNa2to4Bit ();
  }

  if (Na2to4Bit == NULL) return buf;

  /* now return 4 byte string for each compressed byte */

  for (k = 0; k < total; k++) {
    byte = *bytep;
    bytep++;
    index = 4 * byte;
    bp = (Uint4Ptr) (Na2to4Bit + index);
    /* copy 4 bytes at a time */
    /*
    for (j = 0; j < 4; j++) {
      *ptr = *bp;
      ptr++;
      bp++;
    }
    */
    *ptr = *bp;
    ptr++;
  }

  return ptr;
}

static void InitNa4to4Bit (void)

{
  Int2  base [2], index, j;
  Char  convert [16] = {15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  Int4  ret;
  Uint1Ptr Na4to4Bit_local = NULL;

  ret = NlmMutexLockEx (&seqport_mutex);  /* protect this section */
  if (ret) {
    ErrPostEx (SEV_FATAL, 0, 0, "MapNa4ByteToIUPACString mutex failed [%ld]", (long) ret);
    return;
  }

  if (Na4to4Bit == NULL) {
    Na4to4Bit_local = MemNew (sizeof (Uint1) * 512);

    if (Na4to4Bit_local != NULL) {
      for (base [0] = 0; base [0] < 16; (base [0])++) {
        for (base [1] = 0; base [1] < 16; (base [1])++) {
          index = 2 * (base [0] * 16 + base [1]);
          for (j = 0; j < 2; j++) {
            Na4to4Bit_local [index + j] = convert [(base [j])];
          }
        }
      }
    }
    Na4to4Bit = Na4to4Bit_local;
  }

  NlmMutexUnlock (seqport_mutex);
}

NLM_EXTERN Uint2Ptr LIBCALL MapNa4ByteTo4BitString (Uint1Ptr bytep, Uint2Ptr buf, Int4 total)

{
  Uint2Ptr  bp;
  Uint1     byte;
  Int2      index;
  Int4      k;
  Uint2Ptr  ptr;

  if (bytep == NULL || buf == NULL) return buf;
  ptr = buf;

  /* initialize array if not yet set (first time function is called) */

  if (Na4to4Bit == NULL) {
    InitNa4to4Bit ();
  }

  if (Na4to4Bit == NULL) return buf;

  /* now return 2 character string for each compressed byte */

  for (k = 0; k < total; k++) {
    byte = *bytep;
    bytep++;
    index = 2 * byte;
    bp = (Uint2Ptr)  (Na4to4Bit + index);
    /* copy 2 bytes at a time */
    /*
    for (j = 0; j < 2; j++) {
      *ptr = *bp;
      ptr++;
      bp++;
    }
    */
    *ptr = *bp;
    ptr++;
  }

  return ptr;
}

/*****************************************************************************
*
*   SeqPort Routines
*
*****************************************************************************/

/*****************************************************************************
*
*   SeqPortFree(spp)
*
*****************************************************************************/
NLM_EXTERN SeqPortPtr SeqPortFree (SeqPortPtr spp)

{
    SeqPortPtr tspp, nextspp;

    if (spp == NULL)
        return NULL;

    if (spp->locked)              /* locked during access */
        BioseqUnlock(spp->bsp);   /* make available for freeing */

    tspp = spp->segs;
    while (tspp != NULL)
    {
        nextspp = tspp->next;
        SeqPortFree(tspp);
        tspp = nextspp;
    }

    MemFree(spp->cache);
    MemFree (spp->cacheq);

    MemFree(spp);

    return NULL;
}

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
NLM_EXTERN Boolean LIBCALL SeqPortSetValues (SeqPortPtr spp)
{
    SeqPortPtr tmp;

    if (spp == NULL)
        return FALSE;

    for (tmp = spp->segs; tmp != NULL; tmp = tmp->next)
    {
        tmp->is_circle = spp->is_circle;
        tmp->is_seg = spp->is_seg;
        tmp->do_virtual = spp->do_virtual;
        tmp->gapIsZero = spp->gapIsZero;

        if (tmp->segs != NULL)
            SeqPortSetValues(tmp);
    }

    return TRUE;
}


NLM_EXTERN Boolean LIBCALL SeqPortSet_is_circle (SeqPortPtr spp, Boolean value)
{
    if (spp == NULL)
        return FALSE;
    spp->is_circle = value;
    return SeqPortSetValues(spp);
}

NLM_EXTERN Boolean LIBCALL SeqPortSet_is_seg (SeqPortPtr spp, Boolean value)
{
    if (spp == NULL)
        return FALSE;
    spp->is_seg = value;
    return SeqPortSetValues(spp);
}

/**************************************************************
*
*  This function adjusts the length of seqport to remove virtual
*   segments or add them back as needed
*
**************************************************************/
NLM_EXTERN Boolean LIBCALL SeqPortAdjustLength (SeqPortPtr spp)
{
    SeqPortPtr tmp;
    Int4 len = 0;

    if (spp == NULL)
        return FALSE;


    if (spp->isa_virtual)
    {
        if (spp->do_virtual)
            spp->totlen = spp->stop - spp->start + 1;
        else
            spp->totlen = 0;
        if (spp->totlen == 0)
            spp->isa_null = TRUE;
        else
            spp->isa_null = FALSE;
    }
    else if (spp->segs != NULL)
    {
        for (tmp = spp->segs; tmp != NULL; tmp = tmp->next)
        {
            SeqPortAdjustLength (tmp);
            len += tmp->totlen;
        }
        spp->totlen = len;
    }
    else if (! spp->isa_null)
        spp->totlen = spp->stop - spp->start + 1;
    spp->curpos = -1;  /* reset to unused */

    return TRUE;

}

NLM_EXTERN Boolean LIBCALL SeqPortSet_do_virtualEx (SeqPortPtr spp, Boolean value, Boolean gapIsZero)
{
    Boolean do_it = FALSE, has_virtual=FALSE;
    SeqPortPtr tmp;

    if (spp == NULL)
        return FALSE;
    
    if (spp->isa_virtual == TRUE)
        has_virtual = TRUE;
    if (spp->do_virtual != value)
        do_it = TRUE;
    if (spp->gapIsZero != gapIsZero)
        do_it = TRUE;
    for (tmp = spp->segs; tmp != NULL; tmp = tmp->next)
    {
        if (tmp->isa_virtual == TRUE)
            has_virtual = TRUE;
        if (tmp->do_virtual != value)
            do_it = TRUE;
        if (tmp->gapIsZero != gapIsZero)
            do_it = TRUE;
    }

    if (! do_it)   /* no change needed */
        return TRUE;

    
    spp->do_virtual = value;
    spp->gapIsZero = gapIsZero;
    SeqPortSetValues(spp);
    if (has_virtual)   /* have to check the SeqPort */
    {
        SeqPortAdjustLength(spp);
        SeqPortSeek(spp, 0, SEEK_SET);
    }

    return TRUE;
}


NLM_EXTERN Boolean LIBCALL SeqPortSet_do_virtual (SeqPortPtr spp, Boolean value)
{
    return SeqPortSet_do_virtualEx (spp, value, FALSE);
}

NLM_EXTERN Boolean LIBCALL SeqPortSetUpFields (SeqPortPtr spp, Int4 start, Int4 stop, Uint1 
strand, Uint1 newcode)
{
    if (spp == NULL) return FALSE;
    spp->start = start;
    spp->stop = stop;
    spp->strand = strand;
    spp->curpos = -1;    /* not set */
    spp->totlen = stop - start + 1;
    spp->newcode = newcode;
    spp->sctp = SeqCodeTableFind(newcode);

    return TRUE;
}
NLM_EXTERN Boolean LIBCALL SeqPortSetUpAlphabet(SeqPortPtr spp, Uint1 curr_code, Uint1 
newcode)
{
    if (spp == NULL) return FALSE;

        spp->oldcode = curr_code;
        spp->sctp = SeqCodeTableFind(curr_code);

        switch (curr_code)
        {
            case Seq_code_ncbi2na:
                spp->bc = 4;            /* bit shifts needed */
                spp->rshift = 6;
                spp->lshift = 2;
                spp->mask = 192;
                break;
            case Seq_code_ncbi4na:
                spp->bc = 2;
                spp->rshift = 4;
                spp->lshift = 4;
                spp->mask = 240;
                break;
            default:
                spp->bc = 1;
                spp->mask = 255;
                break;
        }

        if ((newcode) && (newcode != curr_code))    /* conversion alphabet */
        {
            if ((spp->smtp = SeqMapTableFind(newcode, curr_code)) != NULL)
                spp->sctp = SeqCodeTableFind(newcode);
        }

        return TRUE;
}

/*****************************************************************************
*
*   SeqPortNew(bsp, start, stop, strand, newcode)
*       if bsp == NULL, creates an empty port
*       see objloc.h for strand defines
*
*****************************************************************************/
NLM_EXTERN SeqPortPtr SeqPortNew (BioseqPtr bsp, Int4 start, Int4 stop, Uint1 strand, Uint1 
newcode)

{
    SeqPortPtr spp, spps, sppcurr = NULL, spprev, prev, curr = NULL;
    Uint1 curr_code, repr, tstrand = 0;
    SeqLocPtr the_segs = NULL, currseg;
    Int4 len, ctr, tlen = 0, tfrom = 0, tto = 0, xfrom, xto, tstart, tstop;
    Char errbuf[41], idbuf[41];
    ValNode fake;
    Boolean done, started;
    BioseqPtr tbsp;
    ValNodePtr currchunk;  /* can be a SeqLoc or an element of a Delta Seq 
*/
    Boolean do_multi_loc, cycle2;
    SeqLitPtr slitp = NULL;
    SeqIdPtr tsip;

    spp = (SeqPortPtr) MemNew(sizeof(SeqPort));
    errbuf[0] = '\0';

    if (bsp == NULL)     /* a NULL section */
        return spp;

    spp->bsp = bsp;                    /* get ready for error 
msgs */
    SeqIdWrite(SeqIdFindBest(bsp->id, 0), errbuf, PRINTID_FASTA_SHORT, 40);
    len = BioseqGetLen(bsp);
    if (start < 0)
    {
        ErrPostEx(SEV_ERROR, 0,0  ,
                 "SeqPortNew: %s start (%ld)< 0", errbuf, 
(long)start);
        MemFree(spp);
        return NULL;
    }
    if (start >= len)
    {
        ErrPostEx(SEV_ERROR,0,0,
                 "SeqPortNew: %s start(%ld) >= len(%ld)",
                    errbuf, (long)start, (long)len);
        MemFree(spp);
        return NULL;
    }
    if (stop == LAST_RESIDUE)
        stop = len - 1;
    else if (stop < start)
    {
        ErrPostEx(SEV_ERROR,0,0,
                 "SeqPortNew: %s stop(%ld) < start(%ld)",
                    errbuf, (long)stop, (long)start);
        MemFree(spp);
        return NULL;
    }
    else if (stop >= len)
    {
        ErrPostEx(SEV_ERROR,0,0,
                 "SeqPortNew: %s stop(%ld) >= len(%ld)",
                    errbuf, (long)stop, (long)len);
        MemFree(spp);
        return NULL;
    }

    SeqPortSetUpFields (spp, start,stop, strand, newcode);

    spp->currnum = BioseqGetSeqDescr(bsp, Seq_descr_num, NULL);
    if (spp->currnum == NULL)   /* no numbering set */
        spp->currnum = NumberingDefaultGet();   /* use default */

    repr = Bioseq_repr(bsp);
    if ((repr == Seq_repr_virtual) ||    /* virtual sequence */
        (repr == Seq_repr_map ))         /* map sequence */
    {
        spp->isa_virtual = TRUE;
        spp->curpos = 0;
    }
    else if ((repr == Seq_repr_seg) ||   /* segmented */
        (repr == Seq_repr_ref) ||     /* reference */
        (repr == Seq_repr_delta))     /* delta */
    {
        spp->oldcode = 0;        /* no code, not raw */

        if (repr == Seq_repr_seg)  /* segmented */
        {
            fake.choice = SEQLOC_MIX;   /* make SEQUENCE OF Seq-loc, 
into one */
            fake.data.ptrvalue = bsp->seq_ext;
            fake.next = NULL;
            the_segs = (SeqLocPtr)&fake;
        }
        else if (repr == Seq_repr_ref)        /* reference: is a Seq-loc 
*/
            the_segs = (SeqLocPtr)bsp->seq_ext;

        if (repr == Seq_repr_delta)   /* chain of deltas to follow */
            currchunk = (ValNodePtr)(bsp->seq_ext);
        else                          /* seqlocs */
            currchunk = (ValNodePtr)SeqLocFindNext(the_segs, NULL);

        currseg = NULL;
        ctr = 0;
        done = FALSE;
        started = FALSE;
        while ((! done) && (currchunk != NULL))
        {
            do_multi_loc = FALSE;
            cycle2 = TRUE;     /* only really needed for complicated 
delta seq locs */
            currseg = NULL;
            if (repr == Seq_repr_delta)
            {
                if (currchunk->choice == 1)  /* it's a SeqLocPtr 
*/
                {
                    currseg = 
(SeqLocPtr)(currchunk->data.ptrvalue);
                    if (! IS_one_loc(currseg, FALSE)) /* 
don't do complicated cases here */
                    {
                        do_multi_loc = TRUE;
                        currseg = 
SeqLocFindNext((SeqLocPtr)(currchunk->data.ptrvalue), NULL);
                    }
                }
                else                         /* it's a SeqLitPtr 
*/
                {
                    currseg = NULL;
                    slitp = 
(SeqLitPtr)(currchunk->data.ptrvalue);
                    tlen = slitp->length;
                    tstrand = Seq_strand_plus;
                    tfrom = 0;
                    tto = tlen - 1;
                }
            }
            else
                currseg = (SeqLocPtr)currchunk;

            while (cycle2)   /* normally once, except for 
complicated delta locs */
            {
                if (currseg != NULL)   /* for segs and deltas of 
type loc */
                {
                    tlen = SeqLocLen(currseg);
                    tstrand = SeqLocStrand(currseg);
                    tfrom = SeqLocStart(currseg);
                    tto = SeqLocStop(currseg);
                }
        
                if (! started)
                {
                    if ((ctr + tlen - 1) >= start)
                    {
                        tstart = start - ctr;
                        started = TRUE;
                    }
                    else
                        tstart = -1;
                }
                else
                    tstart = 0;

                if (tstart >= 0)   /* have a start */
                {
                    if ((ctr + tlen - 1) >= stop)
                    {
                        done = TRUE;   /* hit the end */
                        tstop = ((ctr + tlen - 1) - 
stop);
                    }
                    else
                        tstop = 0;

                    if (tstrand == Seq_strand_minus)
                    {
                        xfrom = tfrom + tstop;
                        xto = tto - tstart;
                    }
                    else
                    {
                        xfrom = tfrom + tstart;
                        xto = tto - tstop;
                    }

                    if (currseg != NULL)    /* working off locs */
                    {
                        if (currseg->choice == SEQLOC_NULL)
                        {
                            tbsp = NULL;
                            spps = SeqPortNew(tbsp, xfrom, xto, tstrand, newcode);
                            spps->isa_null = TRUE;
                        }
                        else
                        {
                            tsip = SeqLocId(currseg);
                            tbsp = BioseqLockById(tsip);
                            if (tbsp != NULL)
                                spps = SeqPortNew(tbsp, xfrom, xto, tstrand, newcode);
                            else
                            {
                                spps = NULL;
                                if (tsip != NULL)
                                    SeqIdWrite(tsip, idbuf, PRINTID_FASTA_SHORT, 40);
                                else
                                    StringMove(idbuf,"seqid=NULL");
                                ErrPostEx(SEV_ERROR,0,0,
                                      "SeqPortNew: %s could not find component %s", 
                                      errbuf, idbuf);
                                return SeqPortFree(spp);
                            }
                        }

                    }
                    else
                    {
                        spps = (SeqPortPtr) MemNew(sizeof(SeqPort));
                        SeqPortSetUpFields (spps, xfrom, 
xto, tstrand, newcode);
                        SeqPortSetUpAlphabet(spps, 
slitp->seq_data_type, newcode);
                        if (slitp->seq_data != NULL)
                            spps->bp = (ByteStorePtr)
slitp->seq_data;
                        else
                        {
                            spps->isa_virtual = TRUE;
                                if (slitp->length == 0)
                                spps->isa_null = TRUE;
                            else
                            {   /* default for delta gaps */
                                spps->do_virtual = TRUE;
                            }

                        }
                    }
        
                    if (spps == NULL)
                    {
                        ErrPostEx(SEV_ERROR,0,0,
                         "SeqPortNew: %s unexpected null during recursion", 
                                 errbuf);
                        return SeqPortFree(spp);
                    }

                    if (currseg != NULL)
                        spps->locked = TRUE;

                    if (sppcurr == NULL)
                        spp->segs = spps;
                    else
                        sppcurr->next = spps;
                    sppcurr = spps;
                }

                ctr += tlen;

                if (! do_multi_loc)
                    cycle2 = FALSE;
                else
                {
                    currseg = 
SeqLocFindNext((SeqLocPtr)(currchunk->data.ptrvalue), currseg);
                    if (currseg == NULL)
                        cycle2 = FALSE;
                }
            }

            if (repr == Seq_repr_delta)
                currchunk = currchunk->next;
            else
                currchunk = SeqLocFindNext(the_segs, currchunk);
        }
        if (strand == Seq_strand_minus)  /* reverse seqport order */
        {
            prev = spp->segs;
            spprev = spp->segs;
            spp->segs = NULL;
            sppcurr = NULL;
            while (prev != NULL)
            {
                curr = spprev;
                prev = NULL;
                while (curr->next != NULL)  /* end of chain */
                {
                    prev = curr;
                    curr = curr->next;
                }
                if (prev != NULL)
                    prev->next = NULL;
                if (sppcurr == NULL)
                    spp->segs = curr;
                else
                    sppcurr->next = curr;
                sppcurr = curr;
            }
            curr->next = NULL;   /* last one in chain */
        }
        spp->curr = spp->segs;

          if (! started)   /* nothing found */
          {
               ErrPostEx(SEV_ERROR,0,0,"SeqPortNew: no data found for %s", 
                           errbuf);
             return SeqPortFree(spp);
          }
    }
    else if ((repr == Seq_repr_raw) ||   /* sequence not by reference */
        (repr == Seq_repr_const))
    {
        curr_code = BioseqGetCode(bsp);

        SeqPortSetUpAlphabet(spp, curr_code, newcode);
        spp->bp = (ByteStorePtr) bsp->seq_data;

     /* allocate fast lookup caches for 2na or 4na to iupacna or 4na conversion */

        if ((newcode == Seq_code_iupacna || newcode == Seq_code_ncbi4na) &&
            (curr_code == Seq_code_ncbi2na || curr_code == Seq_code_ncbi4na)) {
            spp->cacheq = (SPCacheQPtr) MemNew (sizeof (SPCacheQ));
        }

    }

    SeqPortAdjustLength (spp);
    SeqPortSeek(spp, 0, SEEK_SET);
    return spp;
}

/*****************************************************************************
*
*   SeqPortNewByLoc(loc, code)
*       builds a new seqport based on a SeqLoc
*
*****************************************************************************/
NLM_EXTERN SeqPortPtr SeqPortNewByLoc (SeqLocPtr loc, Uint1 code)

{
    BioseqPtr bsp = NULL;
    SeqPortPtr spp = NULL, sppcurr, spps;
    Int4 start = 0, stop = 0;
    Uint1 strand = Seq_strand_unknown;
    SeqLocPtr currloc = NULL;
    CharPtr locptr, currlocptr;

    if (loc == NULL)
        return spp;

                   /* get the needed components */

    switch (loc->choice)
    {
        case SEQLOC_INT:      /* int */
        case SEQLOC_PNT:      /* pnt */
        case SEQLOC_PACKED_PNT:      /* packed-pnt   */
            start = SeqLocStart(loc);
            stop = SeqLocStop(loc);
            strand = SeqLocStrand(loc);
        case SEQLOC_WHOLE:      /* whole */
            bsp = BioseqLockById(SeqLocId(loc));  /* need the bioseq 
now */
            if (bsp == NULL)
                return NULL;    /* can't do it */
    }



    switch (loc->choice)
    {
        case SEQLOC_EMPTY:      /* empty */
        case SEQLOC_EQUIV:     /* equiv */
        case SEQLOC_BOND:      /* bond */
            break;

        case SEQLOC_NULL:      /* null */
            spp = SeqPortNew(NULL, FIRST_RESIDUE, LAST_RESIDUE, 0, 
code);
            spp->isa_null = TRUE;
            break;

        case SEQLOC_WHOLE:      /* whole */
            spp = SeqPortNew(bsp, FIRST_RESIDUE, LAST_RESIDUE, 0, code);
            if (spp != NULL)
                spp->locked = TRUE;
            else
                BioseqUnlock(bsp);
            break;

        case SEQLOC_INT:      /* int */
        case SEQLOC_PNT:      /* pnt */
        case SEQLOC_PACKED_PNT:      /* packed-pnt   */
            spp = SeqPortNew(bsp, start, stop, strand, code);
            if (spp != NULL)
                spp->locked = TRUE;
            else
                BioseqUnlock(bsp);
            break;

        case SEQLOC_PACKED_INT:      /* packed seqint */
        case SEQLOC_MIX:      /* mix */
            spp = (SeqPortPtr) MemNew(sizeof(SeqPort));
            spp->totlen = SeqLocLen(loc);
            spp->start = 0;
            spp->stop = spp->totlen - 1;
            spp->curpos = -1;    /* not set */
            spp->currnum = NULL;   /* use numbering from parts */
            currloc = NULL;
            sppcurr = NULL;
            while ((currloc = SeqLocFindNext(loc, currloc)) != NULL)
            {
                spps = SeqPortNewByLoc(currloc, code);
                if (spps == NULL)
                {
                    locptr = SeqLocPrint(loc);
                    currlocptr = SeqLocPrint(currloc);
                    ErrPostEx(SEV_ERROR, 0,0  ,
        "SeqPortNewByLoc unexpected null during recursion [loc=%s][curr=%s]",
                    locptr, currlocptr);
                    MemFree(locptr);
                    MemFree(currlocptr);
                    SeqPortFree(spp);
                    return NULL;
                }
                if (sppcurr == NULL)
                    spp->segs = spps;
                else
                    sppcurr->next = spps;
                sppcurr = spps;
            }
            spp->curr = spp->segs;
            break;
        case SEQLOC_FEAT:
            ErrPostEx(SEV_ERROR, 0,0  ,
                 "SeqLocNewByLoc: Seq-loc.feat not supported");
            break;
    }

    SeqPortAdjustLength (spp);
    SeqPortSeek(spp, 0, SEEK_SET);

    return spp;
}

/*****************************************************************************
*
*   SeqPortSeek(spp, offset, origin)
*       works like fseek()
*           returns 0 on success   (weird but true)
*           non-zero on fail
*       uses coordinates 0-(len - 1)  no matter what region seqport covers
*       
*
*****************************************************************************/
static void ClearQCache (SeqPortPtr spp, Int4 sp)

{
    SPCacheQPtr spcpq;

    spcpq = spp->cacheq;
    if (spcpq != NULL) {
        spcpq->ctr = 0;
        spcpq->total = 0; /* clear out cache parameters to force new read */
    }
    spp->curpos = sp;
    spp->byte = SEQPORT_EOF;
}

NLM_EXTERN Int2 SeqPortSeek (SeqPortPtr spp, Int4 offset, Int2 origin)

{
    Int4 sp, curpos, left, pos, lim, diff;
    Boolean plus_strand;
    Uint1 the_byte, the_residue;
    Int2 bitctr;
    SeqPortPtr curspp;
    Uint1Ptr buf;
    SPCachePtr spcp;

    if (spp == NULL)
        return 1;

    spp->eos = FALSE;   /* unset flag set when moving off segment */

                                /* get position as positive offset from 0 */
    if (spp->strand == Seq_strand_minus)
        plus_strand = FALSE;
    else
        plus_strand = TRUE;

    sp = spp->curpos;    /* current offset, 0 - (totlen - 1)  */
    switch (origin)
    {
        case SEEK_SET:
            spp->backing = FALSE;  /* reset.. not backing */
            if ((offset > spp->totlen) || (offset < 0)) {
                ClearQCache (spp, sp);
                return 1;
            }
            sp = offset;
            break;
        case SEEK_CUR:
            if (((sp + offset) > spp->totlen) ||
                ((sp + offset) < 0 ))
            {
        /** check for reverse complement backing **/
        if ((sp + offset < 0) && (offset == -2)) 
        {
            if (spp->backing == 1)
            {
                ClearQCache(spp, -1);
                spp->eos = TRUE; /* note backing off segment */
                return 0;
            }
            if (spp->curpos == -1)  /* not set */
                return 0;
        }

                if (! spp->is_circle) {
                    ClearQCache (spp, sp);
                    return 1;
                }
            }
            else
                sp += offset;
            if (spp->is_circle)
            {
                while (sp >= spp->totlen)   /* circle adjustments */
                    sp -= spp->totlen;
                while (sp < 0)
                    sp += spp->totlen;
            }
            break;
        case SEEK_END:
            if ((ABS(offset) > spp->totlen) || (offset > 0)) {
                ClearQCache (spp, sp);
                return 1;
            }
            sp = spp->totlen + offset;
            break;
        default:
            ClearQCache (spp, sp);
            return 1;
    }

    if (sp == spp->curpos)     /* already in right position */
        return 0;

    if (sp == spp->totlen)    /* seek to EOF */
    {
        spp->curpos = sp;
        spp->byte = SEQPORT_EOF;    /* set to nothing */
        ClearQCache (spp, sp);
        return 0;
    }

    if (spp->oldcode)       /* has data, is raw or const type */
    {

        /* if 2na or 4na to iupacna, now only need fast lookup caches */
        if (spp->cacheq != NULL) {
            ClearQCache (spp, sp);
            return 0; /* bypass remaining code */
        }

    /* original code using cache direct from byte store */

        if (spp->cache == NULL)     /* allocate a cache */
            spp->cache = (SPCachePtr)MemNew(sizeof(SPCache));
        spcp = spp->cache;
        buf = spcp->buf;

        if (plus_strand)
        {
            curpos = sp + spp->start;
            pos = curpos / (Int4) (spp->bc);
            lim = spp->stop / (Int4) (spp->bc);
            diff = lim - pos + 1;
            if (diff > 100)
            {
                diff = 100;
                lim = pos + diff - 1;
            }
            BSSeek(spp->bp, pos, SEEK_SET);
            spcp->total = (Int2) BSRead(spp->bp, (VoidPtr)buf, 
diff);
            spcp->ctr = 0;
            spp->bytepos = lim;
        }
        else
        {
            curpos = spp->stop - sp;
            pos = curpos / (Int4) (spp->bc);
            lim = spp->start / (Int4) (spp->bc);
            diff = pos - lim + 1;
            if (diff > 100)
            {
                diff = 100;
                lim = pos - diff + 1;
            }
            BSSeek(spp->bp, lim, SEEK_SET);
            spcp->total = (Int2) BSRead(spp->bp, (VoidPtr)buf, 
diff);
            spcp->ctr = (Int2)(diff - 1);
            spp->bytepos = lim;
        }
        left = curpos % (Int4) (spp->bc);
        the_byte = spcp->buf[spcp->ctr];
        if ((plus_strand) || (spp->bc == 1))
            the_residue = the_byte;
        else        /* reverse compressed bit orders */
        {
            left = spp->bc - 1 - left;
            the_residue = 0;
            bitctr = spp->bc;
            while (bitctr)
            {
                the_residue |= the_byte & spp->mask;
                bitctr--;
                if (bitctr)
                {
                    the_residue >>= spp->lshift;
                    the_byte <<= spp->lshift;
                }
            }
        }
        bitctr = spp->bc;
        while (left)
        {
            the_residue <<= spp->lshift;
            left--; bitctr--;
        }
        spp->byte = the_residue;
        spp->bitctr = (Uint1) bitctr;
        spp->curpos = sp;
        return 0;
    }
    else if ((spp->isa_virtual) || (spp->isa_null))   /* virtual or NULL */
    {
        spp->curpos = sp;
        return 0;
    }
    else                    /* segmented, reference sequences */
    {
 
        if (spp->backing == 1)  /* check for backing off segment */
        {
            if ((spp->curr->curpos == 1) &&
                (! spp->curr->backing))  /* yup */
            {
                spp->curr->curpos = -1;  /* just set the flag */
                spp->curpos -= 2;
                return 0;                /* no eos needed, -1 
will do */
            }
        }

        curpos = 0;
        curspp = spp->segs;
        if (curspp == NULL) return 1;
        while ((curpos + curspp->totlen) <= sp)
        {
            curpos += curspp->totlen;
            curspp = curspp->next;
            if (curspp == NULL)
                return 1;
        }
        if (plus_strand)
            curpos = sp - curpos;
        else
            curpos = (curspp->totlen - 1) - (sp - curpos);
        curspp->backing = spp->backing;
        if (! SeqPortSeek(curspp, curpos, SEEK_SET))
        {
            curspp->backing = FALSE;
            spp->curr = curspp;
            spp->curpos = sp;
            return 0;
        }
        else
        {
            curspp->backing = FALSE;
            return 1;
        }
    }
}

/*****************************************************************************
*
*   Int4 SeqPortTell(spp)
*
*****************************************************************************/
NLM_EXTERN Int4 SeqPortTell (SeqPortPtr spp)

{
    if (spp == NULL)
        return -1L;

    return spp->curpos;
}

/*****************************************************************************
*
*   SeqPortGetResidue(spp)
*       returns residue at current location in requested codeing
*       SEQPORT_EOF = end of file
*
*****************************************************************************/
static Uint1 LIBCALL SeqPortQuickGetResidue (SeqPortPtr spp, SPCacheQPtr spcpq, Boolean plus_strand)

{
  Uint1    bytes [100];
  Int4     curpos, pos, lim, diff;
  CharPtr  ptr;
  Uint1    residue = INVALID_RESIDUE;
  Int2     total, i, j;

  if (spp == NULL || spcpq == NULL) return INVALID_RESIDUE;

  if (spp->curpos == spp->totlen) return SEQPORT_EOF;

  if (spp->curpos < spp->totlen) {

    if (spcpq->ctr >= spcpq->total) {

      /* read next buffer of bytes */

      if (plus_strand) {

        curpos = spp->curpos + spp->start;
        pos = curpos / (Int4) (spp->bc);
        lim = spp->stop / (Int4) (spp->bc);
        diff = lim - pos + 1;
        if (diff > 100) {
          diff = 100;
          lim = pos + diff - 1;
        }
        BSSeek (spp->bp, pos, SEEK_SET);
        total = (Int2) BSRead (spp->bp, (VoidPtr) bytes, diff);
        spp->bytepos = lim;

      } else {

        curpos = spp->stop - spp->curpos;
        pos = curpos / (Int4) (spp->bc);
        lim = spp->start / (Int4) (spp->bc);
        diff = pos - lim + 1;
        if (diff > 100) {
          diff = 100;
          lim = pos - diff + 1;
        }
        BSSeek (spp->bp, lim, SEEK_SET);
        total = (Int2) BSRead (spp->bp, (VoidPtr) bytes, diff);
        spp->bytepos = lim;

      }

      /* buffer is not null terminated, so uses special copy function */

      ptr = spcpq->buf;

      if (spp->newcode == Seq_code_iupacna) {
        if (spp->oldcode == Seq_code_ncbi2na) {
          ptr = (CharPtr) MapNa2ByteToIUPACString (bytes, (Uint4Ptr) ptr, total);
        } else if (spp->oldcode == Seq_code_ncbi4na) {
          ptr = (CharPtr) MapNa4ByteToIUPACString (bytes, (Uint2Ptr) ptr, total);
        }
      } else if (spp->newcode == Seq_code_ncbi4na) {
        if (spp->oldcode == Seq_code_ncbi2na) {
          ptr = (CharPtr) MapNa2ByteTo4BitString (bytes, (Uint4Ptr) ptr, total);
        } else if (spp->oldcode == Seq_code_ncbi4na) {
          ptr = (CharPtr) MapNa4ByteTo4BitString (bytes, (Uint2Ptr) ptr, total);
        }
      }

      spcpq->total = ptr - spcpq->buf;
      spcpq->ctr = 0;

      /* deal with end conditions */

      if (plus_strand) {
        spcpq->ctr += (curpos % (Int4) (spp->bc));
        if (lim == (spp->stop / (Int4) (spp->bc))) {
          diff = (spp->stop + 1) % (Int4) (spp->bc);
          if (diff > 0) {
            spcpq->total -= (Int4) (spp->bc) - diff;
          }
        }
      } else {
        if (pos == (curpos / (Int4) (spp->bc))) {
          diff = (curpos + 1) % (Int4) (spp->bc);
          if (diff > 0) {
            spcpq->total -= (Int4) (spp->bc) - diff;
          }
        }
        if (lim == (spp->start / (Int4) (spp->bc))) {
          spcpq->ctr += (spp->start) % (Int4) (spp->bc);
        }

        /* reverse complement */

        for (i = spcpq->ctr, j = spcpq->total - 1; i < j; i++, j--) {
          residue = spcpq->buf [i];
          spcpq->buf [i] = spcpq->buf [j];
          spcpq->buf [j] = residue;
        }
        for (i = spcpq->ctr; i < spcpq->total; i++) {
          residue = spcpq->buf [i];
          spcpq->buf [i] = SeqCodeTableComp (spp->sctp, residue);
        }

      }

    }

    /* now get residue directly from uncompressed buffer */

    residue = spcpq->buf [spcpq->ctr];
    spcpq->ctr++;
  }

  spp->curpos++;

  return residue;
}

NLM_EXTERN Uint1 LIBCALL SeqPortGetResidue (SeqPortPtr spp)

{
    Uint1 residue = INVALID_RESIDUE, the_byte, the_residue, the_code;
    Boolean plus_strand = TRUE, moveup;
    Int2 bitctr, index;
    Int4 pos, lim, diff;
    SPCachePtr spcp;
    SeqPortPtr tmp, prev;
    SPCacheQPtr spcpq;

    if (spp != NULL)
        spp->backing = FALSE;  /* clear it on read */

    if (spp != NULL && spp->cacheq != NULL && spp->curpos < spp->totlen) {
        spcpq = spp->cacheq;
        if (spcpq->ctr < spcpq->total) {
            residue = spcpq->buf [spcpq->ctr];
            spcpq->ctr++;
            spp->curpos++;
            return residue;
        }
    }

    if ((spp == NULL) || ((spp->bp == NULL) && (spp->oldcode)))
        return SEQPORT_EOF;

    if (spp->isa_null) { /* NULL interval */
        spp->eos = TRUE; /* moving off the segment */
        return SEQPORT_VIRT;
    }

    if (spp->eos)       /* end of reverse complement spp */
        return SEQPORT_EOF;

    if (spp->curpos == spp->totlen)
    {
        if (spp->is_circle)
        {
            SeqPortSeek(spp, 0, SEEK_SET);  /* go to start */
            if (spp->is_seg)   /* give EOS? */
                return SEQPORT_EOS;
        }
        else
            return SEQPORT_EOF;         /* EOF really */
    }
 
    if (spp->curpos == -1)        /* backed off end */
    {
        if (spp->is_circle)
        {
            SeqPortSeek(spp, -1, SEEK_END);  /* go to end */
            if (spp->is_seg)   /* give EOS? */
                return SEQPORT_EOS;
        }
        else
            return SEQPORT_EOF;         /* EOF really */
    }

    if (spp->strand == Seq_strand_minus)
        plus_strand = FALSE;

    if (spp->oldcode)    /* its a raw or const sequence */
    {

    /* separate function for quick lookup to avoid cluttering old code */
        if (spp->cacheq != NULL) {
            return SeqPortQuickGetResidue (spp, spp->cacheq, plus_strand);
        }

        residue = spp->byte & spp->mask;
        residue >>= spp->rshift;
        spp->byte <<= spp->lshift;
        spp->bitctr--;
        if (spp->curpos < (spp->totlen - 1))  /* curpos not incremented yet */
        {
            if (spp->bitctr == 0)
            {
                spcp = spp->cache;
                if (! plus_strand) /* need previous byte */
                {
                    spcp->ctr--;
                    if (spcp->ctr < 0)
                    {
                        pos = spp->bytepos - 1;
                        lim = spp->start / 
(Int4)(spp->bc);
                        diff = pos - lim + 1;
                        if (diff > 100)
                        {
                            diff = 100;
                            lim = pos - 100 + 1;
                        }
                        BSSeek(spp->bp, lim, SEEK_SET);
                        spcp->total = 
(Int2)BSRead(spp->bp, (VoidPtr)(spcp->buf), diff);
                        spcp->ctr = (Int2)(diff - 1);
                        spp->bytepos = lim;
                    }
                }
                else                /* need next 
byte */
                {
                    spcp->ctr++;
                    if (spcp->ctr >= spcp->total)
                    {
                        pos = spp->bytepos + 1;
                        lim = spp->stop / 
(Int4)(spp->bc);
                        diff = lim - pos + 1;
                        if (diff > 100)
                        {
                            diff = 100;
                            lim = pos + diff - 1;
                        }
                        BSSeek(spp->bp, pos, SEEK_SET);
                        spcp->total = 
(Int2)BSRead(spp->bp, (VoidPtr)(spcp->buf), diff);
                        spcp->ctr = 0;
                        spp->bytepos = lim;
                    }
                }
                the_byte = spcp->buf[spcp->ctr];

                if ((plus_strand) || (spp->bc == 1))
                    the_residue = the_byte;
                else        /* reverse compressed bit orders */
                {
                    the_residue = 0;
                    bitctr = spp->bc;
                    while (bitctr)
                    {
                        the_residue |= the_byte & spp->mask;
                        bitctr--;
                        if (bitctr)
                        {
                            the_residue >>= spp->lshift;
                            the_byte <<= spp->lshift;
                        }
                    }       
                }
                spp->byte = the_residue;
                spp->bitctr = spp->bc;
            }
        }

        if (spp->smtp == NULL)   /* no conversion, check now */
        {
            index = (Int2)residue - (Int2)(spp->sctp->start_at);
            if ((index < 0) || (index >= (Int2)(spp->sctp->num)))
                residue = INVALID_RESIDUE;
            else if (*(spp->sctp->names[index]) == '\0')
                residue = INVALID_RESIDUE;
        }
    }
    else if (spp->isa_virtual)  /* virtual */
    {
        if (spp->do_virtual)
        {
            if (spp->newcode)
                the_code = spp->newcode;
            else
                the_code = spp->oldcode;
            if (spp->gapIsZero && the_code == Seq_code_ncbi4na) {
                residue = 0;
            } else {
                residue = GetGapCode (the_code);
            }
            spp->curpos++;
            return residue;
        }
        else
        {
            spp->curpos++;
            return SEQPORT_VIRT;
        }
    }
    else              /* segmented or reference sequence */
    {
        residue = SeqPortGetResidue(spp->curr);
        while (! IS_residue(residue))
        {
            /* spp->curr->eos = FALSE;  just in case was set */
            moveup = FALSE;

            switch (residue)
            {
                case SEQPORT_VIRT:
                case SEQPORT_EOS:
                    if (spp->curr->segs == NULL)  /* this 
did not come up a layer */
                        moveup = TRUE;
                    break;
                case SEQPORT_EOF:
                    moveup = TRUE;
                    break;
                default:
                    break;
            }

            if (moveup)
            {
                if ((spp->curr->curpos == -1) && (! 
spp->curr->eos))   /* moving backwards, many layers deep */
                {
                    prev = NULL;
                    for (tmp = spp->segs; tmp != spp->curr; 
tmp = tmp->next)
                        prev = tmp;
                    if (prev != NULL)
                        spp->curr = prev;
                    else if (spp->is_circle)  /* go to end 
*/
                    {
                        for (tmp = spp->segs; tmp->next 
!= NULL; tmp = tmp->next)
                            continue;
                        spp->curr = tmp;
                    }
                    else
                        return SEQPORT_EOF;

                    if (! plus_strand)
                        SeqPortSeek(spp->curr, 0, 
SEEK_SET);
                    else if (! (spp->curr->isa_null))
                        SeqPortSeek(spp->curr, -1, 
SEEK_END);
                    else
                        spp->curr->curpos = -1;   /* 
flag the null for next time around */
                }
                else                           /* moving 
forwards */
                {
                    if (spp->curr->next != NULL)
                        spp->curr = spp->curr->next;
                    else if (spp->is_circle)
                        spp->curr = spp->segs;
                    else
                        return SEQPORT_EOF;

                    if (plus_strand)
                        SeqPortSeek(spp->curr, 0, 
SEEK_SET);
                    else
                        SeqPortSeek(spp->curr, -1, 
SEEK_END);
                }

                if (spp->is_seg)
                    return SEQPORT_EOS;
            }

            if ((residue == SEQPORT_VIRT) || (residue == 
INVALID_RESIDUE))
                return residue;
            residue = SeqPortGetResidue(spp->curr);
        }

        if (! plus_strand)
        {
            spp->curr->backing++;     /* signal we are backing 
up */
            if (SeqPortSeek(spp->curr, -2, SEEK_CUR))  /* back up to "next" */
                spp->curr->eos = TRUE;

        }
    }
    
    if (spp->smtp != NULL)
        residue = SeqMapTableConvert(spp->smtp, residue);

    if (! plus_strand)
        residue = SeqCodeTableComp(spp->sctp, residue);

    spp->curpos++;
    return residue;
}

/*****************************************************************************
*
*   GetGapCode(seqcode)
*       returns code to use for virtual sequence residues for sequence
*         code seqcode
*       returns INVALID_RESIDUE if seqcode invalid
*
*****************************************************************************/
NLM_EXTERN Uint1 GetGapCode (Uint1 seqcode)
{
    Uint1 residue = INVALID_RESIDUE;
    
    switch (seqcode)
    {
        case Seq_code_iupacna:
            residue = 'N';
            break;
        case Seq_code_iupacaa:
        case Seq_code_ncbieaa:
            residue = 'X';
            break;
        case Seq_code_ncbi2na:    /* there isn't ambiguity */
            break;
        case Seq_code_ncbi8na:
        case Seq_code_ncbi4na:
            residue = 15;
            break;
        case Seq_code_iupacaa3:  /* no 1 letter character */
        case Seq_code_ncbipna:
        case Seq_code_ncbipaa:
            break;
        case Seq_code_ncbistdaa:
            residue = 21;
            break;

    }

    return residue;
}


/*****************************************************************************
*
*   SeqPortRead(spp, buf, len)
*       returns bytes read
*       if returns a negative number, then ABS(return value) gives the
*         same codes as SeqPortGetResidue for EOF or EOS
*
*****************************************************************************/
NLM_EXTERN Int2 LIBCALL SeqPortRead (SeqPortPtr spp, Uint1Ptr buf, Int2 len)

{
    Int2 ctr = 0;
    Int4 loopmax;
    Uint1 retval;
    SPCacheQPtr spcpq;

    if ((spp == NULL) || (buf == NULL) || (len <= 0))
        return 0;

    if (spp->lastmsg)    /* previous EOF or EOS saved */
    {
        ctr = spp->lastmsg;
        spp->lastmsg = 0;
        ctr *= -1;
        return ctr;
    }

    spcpq = spp->cacheq;
    while (ctr < len) {
        loopmax = 0;
        if (spcpq != NULL && spp->curpos < spp->totlen && spcpq->ctr < spcpq->total) {
            loopmax = MIN ((spp->totlen - spp->curpos), (spcpq->total - spcpq->ctr));
            loopmax = MIN (loopmax, (Int4) (len - ctr));
        }
        /* loopmax saves multiple comparisons, speeds up significantly */
        if (loopmax > 0) {
            while (loopmax > 0) {
                retval = spcpq->buf [spcpq->ctr];
                spcpq->ctr++;
                spp->curpos++;
                loopmax--;
                if (IS_residue (retval)) {
                    *buf = retval;
                    buf++;
                    ctr++;
                } else {
                    if (! ctr)   /* first one */
                    {
                        ctr = retval;   /* send return as negative number */
                        ctr *= -1;
                        return ctr;
                    } else {
                        spp->lastmsg = retval;
                        return ctr;
                    }
                }
            }
        } else {
            retval = SeqPortGetResidue(spp);
            if (IS_residue(retval))
            {
                *buf = retval;
                buf++;
                ctr++;
            }
            else
            {
                if (! ctr)   /* first one */
                {
                    ctr = retval;   /* send return as negative number */
                    ctr *= -1;
                    return ctr;
                }
                else
                {
                    spp->lastmsg = retval;
                    return ctr;
                }
            }
        }
    }
    return ctr;
}

/*******************************************************************************
*    
*   SeqPortStream (bsp, flags, userdata, proc)
*   SeqPortStreamInt (bsp, start, stop, strand, flags, userdata, proc)
*   SeqPortStreamLoc (slp, flags, userdata, proc)
*       Efficient functions to stream through sequence
*
********************************************************************************/

/* structure for passing common arguments internal functions */

typedef struct streamdata {
  StreamFlgType      flags;
  Pointer            userdata;
  SeqPortStreamProc  proc;
  Uint1              letterToComp [256];
  CharPtr            tmp;
  Boolean            failed;
  Int2               depth;
  SeqEntryPtr        scope;
} StreamData, PNTR StreamDataPtr;

/* prototype for main internal recursive processing function */

static Int4 SeqPortStreamWork (
  BioseqPtr bsp,
  Int4 start,
  Int4 stop,
  Uint1 strand,
  StreamDataPtr sdp
);

#define STREAM_GAP_MASK (STREAM_EXPAND_GAPS | GAP_TO_SINGLE_DASH | EXPAND_GAPS_TO_DASHES)

static Int4 SeqPortStreamGap (
  Int4 length,
  Boolean is_na,
  Boolean is_virt,
  Boolean is_known,
  StreamDataPtr sdp
)

{
  Char     buf [4004];
  Char     ch, gapchar = '-';
  Boolean  expand_gaps, many_dashes, many_pluses, single_dash;
  Int4     len;

  if (sdp == NULL) return 0;

  many_pluses = FALSE;
  if (is_virt) {
    if ((sdp->flags & SUPPRESS_VIRT_SEQ) != 0) return 0;
    if ((sdp->flags & STREAM_VIRT_AS_PLUS) != 0) {
      many_pluses = TRUE;
      gapchar = '+';
    }
  } else if (is_known) {
    if ((sdp->flags & KNOWN_GAP_AS_PLUS) != 0) {
      many_pluses = TRUE;
      gapchar = '+';
    }
  }

  expand_gaps = (Boolean) ((sdp->flags & STREAM_GAP_MASK) == STREAM_EXPAND_GAPS);
  single_dash = (Boolean) ((sdp->flags & STREAM_GAP_MASK) == GAP_TO_SINGLE_DASH);
  many_dashes = (Boolean) ((sdp->flags & STREAM_GAP_MASK) == EXPAND_GAPS_TO_DASHES);

  /* if both gap flags are false, ignore gap */

  if ((! expand_gaps) && (! single_dash) && (! many_dashes)) return 0;

  if (single_dash) {

    /* if only indicating gap presence, send one gap character, return 0 count */

    buf [0] = gapchar;
    buf [1] = '\0';

    sdp->proc (buf, sdp->userdata);

    return 0;
  }

  /* if not single dash to mark any gap, need at least one base or residue */

  if (length < 1) return 0;

  if (many_dashes || many_pluses) {
    ch = gapchar;
  } else if (is_na) {
    ch = 'N';
  } else {
    ch = 'X';
  }

  len = MIN (length, 4000L);
  MemSet ((Pointer) buf, ch, len);
  buf [(int) (Int2) len] = '\0';

  for (len = length; len > 0; len -= 4000L) {

    /* on last loop, send only partial buffer */

    if (len < 4000L) {
      buf [(int) (Int2) len] = '\0';
    }

    sdp->proc (buf, sdp->userdata);
  }

  /* return number of N or X or gap characters sent */

  return length;
}

static Uint1Ptr LIBCALL MapNa8ByteToIUPACString (
  Uint1Ptr bytep,
  Uint1Ptr buf,
  Int4 total,
  Uint1 badchar,
  SeqMapTablePtr smtp,
  StreamDataPtr sdp
)

{
  Uint1     ch;
  Int4      k;
  Uint1Ptr  ptr;
  Uint1     residue;

  if (bytep == NULL || buf == NULL || sdp == NULL) return buf;
  ptr = buf;

  for (k = 0; k < total; k++) {
    residue = *bytep;
    if (smtp != NULL) {
      ch = SeqMapTableConvert (smtp, residue);
      if (ch == INVALID_RESIDUE && (Boolean) ((sdp->flags & STREAM_CORRECT_INVAL) != 0)) {
        *ptr = badchar;
      } else {
        *ptr = ch;
      }
    } else {
      *ptr = residue;
    }
    bytep++;
    ptr++;
  }

  return ptr;
}

static Int4 SeqPortStreamBlock (
  ByteStorePtr bs,
  Int4 blk,
  Int4 compress,
  Uint1 alphabet,
  Uint1 badchar,
  SeqMapTablePtr smtp,
  Int4 start,
  Int4 stop,
  Boolean revcomp,
  StreamDataPtr sdp
)

{
  Uint4     uncomp [1001]; /* 4000 characters + extra for end-of-string null byte */
  Uint4     compr [251];   /* 1000 bytes + extra for safety */
  CharPtr   buf;
  Uint1Ptr  bytes;
  /*
  Char      buf [4004];
  Uint1     bytes [1004];
  */
  Char      ch;
  Int4      count = 0, cumulative, total;
  Int2      from, to;
  Boolean   many_dashes, single_dash;
  CharPtr   nd, ptr, str, tmp;

  if (bs == NULL || sdp == NULL) return 0;

  /* Uint4 arrays ensure 4-byte address alignment by the compiler, no need for & since array is pointer */

  buf = (CharPtr) uncomp;
  bytes = (Uint1Ptr) compr;

  BSSeek (bs, blk, SEEK_SET);

  total = BSRead (bs, (VoidPtr) bytes, 1000L);
  if (total < 1) return 0;

  ptr = buf;
  switch (alphabet) {
    case Seq_code_ncbi2na :
      ptr = (CharPtr) MapNa2ByteToIUPACString (bytes, (Uint4Ptr) ptr, total);
      break;
    case Seq_code_ncbi4na :
      single_dash = (Boolean) ((sdp->flags & STREAM_GAP_MASK) == GAP_TO_SINGLE_DASH);
      many_dashes = (Boolean) ((sdp->flags & STREAM_GAP_MASK) == EXPAND_GAPS_TO_DASHES);
      if (single_dash || many_dashes) {
        ptr = (CharPtr) MapNa4ByteToIUPACplusGapString (bytes, (Uint2Ptr) ptr, total);
      } else {
        ptr = (CharPtr) MapNa4ByteToIUPACString (bytes, (Uint2Ptr) ptr, total);
      }
      break;
    default :
      ptr = (CharPtr) MapNa8ByteToIUPACString (bytes, (Uint1Ptr) ptr, total, badchar, smtp, sdp);
      break;
  }
  *ptr = '\0';

  cumulative = blk * compress;

  /* deal with end conditions */

  total = ptr - buf;

  /* check for bsp->length > actual raw data */

  if (start > cumulative + total) {
    sdp->failed = TRUE;
    return 0;
  }

  from = 0;
  if (start > cumulative && start < cumulative + total) {
    from += start - cumulative;
  }

  if (stop < cumulative + total) {
    to = (Int2) (stop - cumulative + 1);
    buf [to] = '\0';
  }

  str = buf + from;

  if (revcomp) {

    /* reverse string first - middle base not touched, so cannot also complement here */

    nd = str;
    while (*nd != '\0') {
      nd++;
    }
    nd--;

    tmp = str;
    while (nd > tmp) {
      ch = *nd;
      *nd = *tmp;
      *tmp = ch;
      nd--;
      tmp++;
    }

    /* now complement every base in string */

    nd = str;
    ch = *nd;
    while (ch != '\0') {
      *nd = sdp->letterToComp [(int) (Uint1) ch];
      nd++;
      ch = *nd;
    }

  }

  /* send characters to stream callback */

  sdp->proc (str, sdp->userdata);

  /* return number of characters sent */

  tmp = str;
  while (*tmp != '\0') {
    count++;
    tmp++;
  }

  return count;
}

static Int4 SeqPortStreamRaw (
  BioseqPtr bsp,
  Int4 start,
  Int4 stop,
  Uint1 strand,
  StreamDataPtr sdp
)

{
  Uint1           alphabet, code;
  Char            badchar;
  ByteStorePtr    bs;
  Int4            blk, compress, count = 0, from, to;
  Boolean         is_na, revcomp = FALSE;
  SeqMapTablePtr  smtp = NULL;

  if (bsp == NULL || sdp == NULL) return 0;
  if (bsp->repr != Seq_repr_raw && bsp->repr != Seq_repr_const) return 0;

  is_na = (Boolean) ISA_na (bsp->mol);

  if (bsp->seq_data_type == Seq_code_gap) {

    /* support for new Seq-data.gap */

    count += SeqPortStreamGap (stop - start + 1, is_na, FALSE, FALSE, sdp);

    return count;
  }

  /* otherwise Seq-data is a byte store */

  bs = (ByteStorePtr) bsp->seq_data;
  if (bs == NULL) return 0;

  alphabet = bsp->seq_data_type;

  if (strand == Seq_strand_minus && is_na) {
    revcomp = TRUE;
  }

  /* setup code conversion and decompression parameters */

  if (is_na) {
    code = Seq_code_iupacna;
    badchar = 'N';
  } else {
    code = Seq_code_ncbieaa;
    badchar = 'X';
  }

  switch (alphabet) {
    case Seq_code_ncbi2na :
      compress = 4;
      break;
    case Seq_code_ncbi4na :
      compress = 2;
      break;
    default :
      compress = 1;
      break;
  }

  if (code != alphabet) {
    smtp = SeqMapTableFind (code, alphabet);
    if (smtp == NULL) return 0;
  }

  /* calculate bytestore block addresses in chunks of 1000 */

  from = ((start / compress) / 1000L) * 1000L;
  to = ((stop / compress) / 1000L) * 1000L;

  /* process sequential blocks of sequence */

  if (revcomp) {

    for (blk = to; blk >= from; blk -= 1000) {
      count += SeqPortStreamBlock (bs, blk, compress, alphabet, badchar, smtp, start, stop, TRUE, sdp);
    }

  } else {

    for (blk = from; blk <= to; blk += 1000) {
      count += SeqPortStreamBlock (bs, blk, compress, alphabet, badchar, smtp, start, stop, FALSE, sdp);
    }
  }

  return count;
}

static Int4 SeqPortStreamSeqLit (
  SeqLitPtr slitp,
  Boolean is_na,
  Int4 start,
  Int4 stop,
  Uint1 strand,
  StreamDataPtr sdp
)

{
  Bioseq   bsq;
  Int4     count = 0;
  Boolean  is_known = TRUE;

  if (slitp == NULL || sdp == NULL) return 0;

  /* ignore gaps of unknown length */

  if (slitp->length < 1) return 0;

  if (slitp->seq_data == NULL || slitp->seq_data_type == Seq_code_gap) {

    /* literal without sequence data is a virtual gap, also handle new gap type */

    if (slitp->fuzz != NULL) {
      is_known = FALSE;
    }

    count += SeqPortStreamGap (stop - start + 1, is_na, FALSE, is_known, sdp);

    return count;
  }

  /* otherwise fake a Bioseq with the literal as its data */

  MemSet ((Pointer) &bsq, 0, sizeof (Bioseq));

  bsq.repr = Seq_repr_raw;
  if (is_na) {
    bsq.mol = Seq_mol_dna;
  } else {
    bsq.mol = Seq_mol_aa;
  }
  bsq.seq_data_type = slitp->seq_data_type;
  bsq.seq_data = slitp->seq_data;
  bsq.length = slitp->length;

  /* call SeqPortStreamRaw to handle sequence data in the byte store */

  count += SeqPortStreamRaw (&bsq, start, stop, strand, sdp);

  return count;
}

static Int2     stream_retry_attempts = 0;
static Boolean  stream_retry_count_set = FALSE;

static Int2     stream_retry_sleep = 0;
static Boolean  stream_retryp_sleep_set = FALSE;

static Int4 SeqPortStreamSeqLoc (
  SeqLocPtr slp,
  Int4 start,
  Int4 stop,
  Uint1 strand,
  StreamDataPtr sdp,
  SeqIdPtr parentID
)

{
  BioseqPtr    bsp;
  Char         buf [64];
  Int4         count = 0;
  SeqEntryPtr  oldscope = NULL;
  Char         pid [64];
  SeqIdPtr     sip;
#ifdef OS_UNIX
  Int2         attempts;
  CharPtr      str;
  int          val = 0;
#endif

  if (slp == NULL || sdp == NULL) return 0;

  if (start < 0 || stop < 0) return 0;

  sip = SeqLocId (slp);
  if (sip == NULL) return 0;

  if (sip->choice == SEQID_GI && sip->data.intvalue <= 0 &&
      (Boolean) ((sdp->flags & STREAM_ALLOW_NEG_GIS) == 0)) {
    
    /* gi 0 is always a data error, just report and bail */
    /* negative gi sometimes used in-house, allow if flag set */

    SeqIdWrite (sip, buf, PRINTID_FASTA_SHORT, sizeof (buf) - 1);
    if (parentID != NULL) {
      SeqIdWrite (parentID, pid, PRINTID_FASTA_LONG, sizeof (pid) - 1);
      ErrPostEx (SEV_ERROR, 0, 0, "SeqPortStream ignoring Bioseq %s component of %s", buf, pid);
    } else {
      ErrPostEx (SEV_ERROR, 0, 0, "SeqPortStream ignoring Bioseq %s", buf);
    }
    sdp->failed = TRUE;
    return 0;
  }

  oldscope = SeqEntrySetScope (sdp->scope);
  bsp = BioseqLockById (sip);
  SeqEntrySetScope (oldscope);

#ifdef OS_UNIX
  if (bsp == NULL) {

    /* number of retries and sleep between retries now configured by environment variable */

    if (! stream_retry_count_set) {
      str = (CharPtr) getenv ("SEQPORT_STREAM_FETCH_ATTEMPTS");
      if (StringDoesHaveText (str)) {
        if (sscanf (str, "%d", &val) == 1) {
          stream_retry_attempts = (Uint2) val;
        }
      }
      stream_retry_count_set = TRUE;
    }

    if (! stream_retryp_sleep_set) {
      str = (CharPtr) getenv ("SEQPORT_STREAM_RETRY_SLEEP");
      if (StringDoesHaveText (str)) {
        if (sscanf (str, "%d", &val) == 1) {
          stream_retry_sleep = (Uint2) val;
        }
      }
      stream_retryp_sleep_set = TRUE;
    }

    /* retry failed fetch attempt up to specified limit */

    if (stream_retry_attempts > 1) {
      attempts = 1;
      while (bsp == NULL && attempts < stream_retry_attempts) {
        if (stream_retry_sleep > 0) {
          sleep (stream_retry_sleep);
        }

        oldscope = SeqEntrySetScope (sdp->scope);
        bsp = BioseqLockById (sip);
        SeqEntrySetScope (oldscope);
        attempts++;
      }
      if (bsp != NULL) {
        SeqIdWrite (sip, buf, PRINTID_FASTA_SHORT, sizeof (buf) - 1);
        if (parentID != NULL) {
          SeqIdWrite (parentID, pid, PRINTID_FASTA_LONG, sizeof (pid) - 1);
          ErrPostEx (SEV_WARNING, 0, 0,
                     "SeqPortStream loaded Bioseq %s component of %s after %d attempts",
                     buf, pid, (int) attempts);
        } else {
          ErrPostEx (SEV_WARNING, 0, 0,
                     "SeqPortStream loaded Bioseq %s after %d attempts",
                     buf, (int) attempts);
        }
      }
    }
  }
#endif

  if (bsp == NULL) {
    SeqIdWrite (sip, buf, PRINTID_FASTA_SHORT, sizeof (buf) - 1);
    if (parentID != NULL) {
      SeqIdWrite (parentID, pid, PRINTID_FASTA_LONG, sizeof (pid) - 1);
      ErrPostEx (SEV_ERROR, 0, 0, "SeqPortStream failed to load Bioseq %s component of %s", buf, pid);
    } else {
      ErrPostEx (SEV_ERROR, 0, 0, "SeqPortStream failed to load Bioseq %s", buf);
    }
    sdp->failed = TRUE;
    return 0;
  }

  count = SeqPortStreamWork (bsp, start, stop, strand, sdp);

  BioseqUnlock (bsp);

  return count;
}

/* structure for processing components in forward or reverse direction */

typedef struct streamobj {
  SeqLocPtr  slp;
  SeqLitPtr  slitp;
  Int4       from;
  Int4       to;
  Uint1      strand;
} StreamObj, PNTR StreamObjPtr;

static StreamObjPtr StreamObjNew (
  SeqLocPtr slp,
  SeqLitPtr slitp,
  Int4 from,
  Int4 to,
  Uint1 strand
)

{
  StreamObjPtr  sop;

  sop = (StreamObjPtr) MemNew (sizeof (StreamObj));
  if (sop == NULL) return NULL;

  sop->slp = slp;
  sop->slitp = slitp;
  sop->from = from;
  sop->to = to;
  sop->strand = strand;

  return sop;
}

static Int4 SeqPortStreamDelta (
  BioseqPtr bsp,
  Int4 start,
  Int4 stop,
  Uint1 strand,
  StreamDataPtr sdp
)

{
  Int4          count = 0, cumulative, from, to, len;
  DeltaSeqPtr   dsp;
  ValNodePtr    head = NULL, last = NULL, vnp;
  Boolean       is_na;
  Boolean       revcomp = FALSE;
  SeqLitPtr     slitp;
  SeqLocPtr     slp;
  StreamObjPtr  sop;

  if (bsp == NULL || sdp == NULL) return 0;

  is_na = (Boolean) ISA_na (bsp->mol);

  if (strand == Seq_strand_minus && is_na) {
    revcomp = TRUE;
  }

  /* build linked list in forward or reverse order, depending upon input strand */

  for (dsp = (DeltaSeqPtr) bsp->seq_ext, cumulative = 0;
       dsp != NULL && cumulative <= stop;
       dsp = dsp->next, cumulative += len) {

    len = 0;

    switch (dsp->choice) {

      case 1 :
        slp = (SeqLocPtr) dsp->data.ptrvalue;
        if (slp == NULL) continue;

        if (slp->choice == SEQLOC_NULL) continue;

        from = SeqLocStart (slp);
        to = SeqLocStop (slp);
        strand = SeqLocStrand (slp);

        if (from < 0 || to < 0) continue;

        len = to - from + 1;

        if (cumulative + len <= start) continue;

        /* adjust from and to if not using entire interval */

        if (strand == Seq_strand_minus) {

          if (start > cumulative) {
            to -= start - cumulative;
          }

          if (stop < cumulative + len) {
            from += cumulative + len - stop - 1;
          }

        } else {

          if (start > cumulative) {
            from += start - cumulative;
          }

          if (stop < cumulative + len) {
            to -= cumulative + len - stop - 1;
          }
        }

        if (revcomp) {
          if (strand == Seq_strand_minus) {
            strand = Seq_strand_plus;
          } else {
            strand = Seq_strand_minus;
          }
        }

        sop = StreamObjNew (slp, NULL, from, to, strand);
        if (sop == NULL) continue;

        if (revcomp) {

          vnp = ValNodeAddPointer (NULL, 0, (Pointer) sop);
          vnp->next = head;
          head = vnp;

        } else {

          vnp = ValNodeAddPointer (&last, 0, (Pointer) sop);
          if (head == NULL) {
            head = vnp;
          }
          last = vnp;
        }
        break;

      case 2 :
        slitp = (SeqLitPtr) dsp->data.ptrvalue;
        if (slitp == NULL) continue;

        from = 0;
        to = slitp->length - 1;
        strand = Seq_strand_plus;

        if (from < 0 || to < 0) continue;

        len = to - from + 1;

        if (cumulative + len <= start) continue;

        /* adjust from and to if not using entire interval */

        if (start > cumulative) {
          from += start - cumulative;
        }

        if (stop < cumulative + len) {
          to -= cumulative + len - stop - 1;
        }

        if (revcomp) {
          if (strand == Seq_strand_minus) {
            strand = Seq_strand_plus;
          } else {
            strand = Seq_strand_minus;
          }
        }

        sop = StreamObjNew (NULL, slitp, from, to, strand);
        if (sop == NULL) continue;

        if (revcomp) {

          vnp = ValNodeAddPointer (NULL, 0, (Pointer) sop);
          vnp->next = head;
          head = vnp;

        } else {

          vnp = ValNodeAddPointer (&last, 0, (Pointer) sop);
          if (head == NULL) {
            head = vnp;
          }
          last = vnp;
        }
        break;

      default :
        break;
    }
  }

  /* process components in correct order */

  for (vnp = head; vnp != NULL; vnp = vnp->next) {

    sop = (StreamObjPtr) vnp->data.ptrvalue;
    if (sop == NULL) continue;

    if (sop->slp != NULL) {

      count += SeqPortStreamSeqLoc (sop->slp, sop->from, sop->to, sop->strand, sdp, bsp->id);

    } else if (sop->slitp != NULL) {

      count += SeqPortStreamSeqLit (sop->slitp, is_na, sop->from, sop->to, sop->strand, sdp);
    }
  }

  /* free control list */

  ValNodeFreeData (head);

  return count;
}

static Int4 SeqPortStreamSeg (
  BioseqPtr bsp,
  Int4 start,
  Int4 stop,
  Uint1 strand,
  StreamDataPtr sdp
)

{
  Int4          count = 0, cumulative, from, to, len;
  ValNodePtr    head = NULL, last = NULL, vnp;
  Boolean       is_na;
  Boolean       revcomp = FALSE;
  SeqLocPtr     slp;
  StreamObjPtr  sop;

  if (bsp == NULL || sdp == NULL) return 0;

  is_na = (Boolean) ISA_na (bsp->mol);

  if (strand == Seq_strand_minus && is_na) {
    revcomp = TRUE;
  }

  /* build linked list in forward or reverse order, depending upon input strand */

  for (slp = (SeqLocPtr) bsp->seq_ext, cumulative = 0;
       slp != NULL && cumulative <= stop;
       slp = slp->next, cumulative += len) {

    len = 0;

    if (slp->choice == SEQLOC_NULL) continue;

    from = SeqLocStart (slp);
    to = SeqLocStop (slp);
    strand = SeqLocStrand (slp);

    if (from < 0 || to < 0) continue;

    len = to - from + 1;

    if (cumulative + len <= start) continue;

    /* adjust from and to if not using entire interval */

    if (strand == Seq_strand_minus) {

      if (start > cumulative) {
        to -= start - cumulative;
      }

      if (stop < cumulative + len) {
        from += cumulative + len - stop - 1;
      }

    } else {

      if (start > cumulative) {
        from += start - cumulative;
      }

      if (stop < cumulative + len) {
        to -= cumulative + len - stop - 1;
      }
    }

    if (revcomp) {
      if (strand == Seq_strand_minus) {
        strand = Seq_strand_plus;
      } else {
        strand = Seq_strand_minus;
      }
    }

    sop = StreamObjNew (slp, NULL, from, to, strand);
    if (sop == NULL) continue;

    if (revcomp) {

      vnp = ValNodeAddPointer (NULL, 0, (Pointer) sop);
      vnp->next = head;
      head = vnp;

    } else {

      vnp = ValNodeAddPointer (&last, 0, (Pointer) sop);
      if (head == NULL) {
        head = vnp;
      }
      last = vnp;
    }
  }

  /* process components in correct order */

  for (vnp = head; vnp != NULL; vnp = vnp->next) {

    sop = (StreamObjPtr) vnp->data.ptrvalue;
    if (sop == NULL) continue;

    if (sop->slp != NULL) {

      count += SeqPortStreamSeqLoc (sop->slp, sop->from, sop->to, sop->strand, sdp, bsp->id);
    }
  }

  /* free control list */

  ValNodeFreeData (head);

  return count;
}

/* SeqPortStreamWork calls appropriate representation-specific function */

static Int4 SeqPortStreamWork (
  BioseqPtr bsp,
  Int4 start,
  Int4 stop,
  Uint1 strand,
  StreamDataPtr sdp
)

{
  Int4  count = 0;

  if (bsp == NULL || sdp == NULL) return 0;

  /* start and stop position reality checks */

  if (start < 0) {
    start = 0;
  }
  if (stop < 0) {
    stop = bsp->length - 1;
  }

  /* if start or stop are beyond sequence length, set failed flag */

  if (start >= bsp->length || stop >= bsp->length) {
    sdp->failed = TRUE;
    return 0;
  }

  if (start > stop) return 0;

  /* stack depth overflow check for recursively-defined sequence instances */

  (sdp->depth)++;

  if (sdp->depth > 20) {
    ErrPostEx (SEV_ERROR, 0, 0, "SeqPortStreamWork stack depth overflow");
    sdp->failed = TRUE;
    return 0;
  }

  /* call appropriate stream function */

  switch (bsp->repr) {

    case Seq_repr_virtual :
      count += SeqPortStreamGap (stop - start + 1, ISA_na (bsp->mol), TRUE, FALSE, sdp);
      break;

    case Seq_repr_raw :
    case Seq_repr_const :
      count += SeqPortStreamRaw (bsp, start, stop, strand, sdp);
      break;

    case Seq_repr_seg :
      if (bsp->seq_ext_type == 1) {
        count += SeqPortStreamSeg (bsp, start, stop, strand, sdp);
      }
      break;

    case Seq_repr_delta :
      if (bsp->seq_ext_type == 4) {
        count += SeqPortStreamDelta (bsp, start, stop, strand, sdp);
      }
      break;

    default :
      break;
  }

  /* restore stack depth value */

  (sdp->depth)--;

  return count;
}

/* default callback for copying to allocated buffer */

static void LIBCALLBACK SaveStreamSequence (
  CharPtr sequence,
  Pointer userdata
)

{
  CharPtr       tmp;
  CharPtr PNTR  tmpp;

  tmpp = (CharPtr PNTR) userdata;
  tmp = *tmpp;

  tmp = StringMove (tmp, sequence);

  *tmpp = tmp;
}

/* SeqPortStreamSetup creates revcomp table, calls SeqPortStreamWork  */

static Int4 SeqPortStreamSetup (
  BioseqPtr bsp,
  Int4 start,
  Int4 stop,
  Uint1 strand,
  SeqLocPtr loc,
  SeqLitPtr lit,
  StreamFlgType flags,
  Pointer userdata,
  SeqPortStreamProc proc
)

{
  Char        ch, lttr;
  CharPtr     complementBase = " TVGH  CD  M KN   YSAABW R ";
  Int4        count = 0, from, to;
  Uint2       entityID;
  Int2        i;
  Boolean     is_na;
  StreamData  sd;
  SeqLocPtr   slp;

  if (bsp == NULL && loc == NULL && lit == NULL) return 0;
  if (proc == NULL && userdata == NULL) return 0;

  MemSet ((Pointer) &sd, 0, sizeof (StreamData));

  sd.flags = flags;
  sd.userdata = userdata;
  sd.proc = proc;
  sd.tmp = NULL;
  sd.failed = FALSE;
  sd.depth = 0;

  /* if NULL callback, copy into allocated userdata string */

  if (proc == NULL) {
    sd.proc = SaveStreamSequence;
    sd.tmp = userdata;
    sd.userdata = &(sd.tmp);
  }

  /* set up nucleotide complementation lookup table */

  for (i = 0; i < 256; i++) {
    sd.letterToComp [i] = '\0';
  }
  for (ch = 'A', i = 1; ch <= 'Z'; ch++, i++) {
    lttr = complementBase [i];
    if (lttr != ' ') {
      sd.letterToComp [(int) (Uint1) ch] = lttr;
    }
  }
  for (ch = 'a', i = 1; ch <= 'z'; ch++, i++) {
    lttr = complementBase [i];
    if (lttr != ' ') {
      sd.letterToComp [(int) (Uint1) ch] = lttr;
    }
  }

  /* commence streaming */

  if (bsp != NULL) {

    entityID = ObjMgrGetEntityIDForPointer (bsp);
    sd.scope = GetTopSeqEntryForEntityID (entityID);

    count += SeqPortStreamWork (bsp, start, stop, strand, &sd);

  } else if (loc != NULL) {

    sd.scope = SeqEntryGetScope ();

    slp = SeqLocFindNext (loc, NULL);
    while (slp != NULL) {

      from = SeqLocStart (slp);
      to = SeqLocStop (slp);
      strand = SeqLocStrand (slp);

      if (from < 0 || to < 0) {
        sd.failed = TRUE;
        return -1;
      }

      count += SeqPortStreamSeqLoc (slp, from, to, strand, &sd, NULL);

      slp = SeqLocFindNext (loc, slp);
    }

  } else if (lit != NULL) {

    is_na = TRUE;
    switch (lit->seq_data_type) {
      case Seq_code_iupacaa :
      case Seq_code_ncbi8aa :
      case Seq_code_ncbieaa :
      case Seq_code_ncbipaa :
      case Seq_code_iupacaa3 :
      case Seq_code_ncbistdaa :
        is_na = FALSE;
        break;
      default :
        break;
    }

    count += SeqPortStreamSeqLit (lit, is_na, 0, lit->length - 1, Seq_strand_plus, &sd);
  }

  /* return number of bases or residues streamed to callback */

  if (sd.failed) {
    if (count  < 1) return -1;
    return -count;
  }

  return count;
}

/* public functions all call SeqPortStreamSetup */

NLM_EXTERN Int4 SeqPortStream (
  BioseqPtr bsp,
  StreamFlgType flags,
  Pointer userdata,
  SeqPortStreamProc proc
)

{
  return SeqPortStreamSetup (bsp, 0, -1, Seq_strand_unknown, NULL, NULL, flags, userdata, proc);
}

NLM_EXTERN Int4 SeqPortStreamInt (
  BioseqPtr bsp,
  Int4 start,
  Int4 stop,
  Uint1 strand,
  StreamFlgType flags,
  Pointer userdata,
  SeqPortStreamProc proc
)

{
  return SeqPortStreamSetup (bsp, start, stop, strand, NULL, NULL, flags, userdata, proc);
}

NLM_EXTERN Int4 SeqPortStreamLoc (
  SeqLocPtr slp,
  StreamFlgType flags,
  Pointer userdata,
  SeqPortStreamProc proc
)

{
  return SeqPortStreamSetup (NULL, 0, 0, 0, slp, NULL, flags, userdata, proc);
}

NLM_EXTERN Int4 SeqPortStreamLit (
  SeqLitPtr lit,
  StreamFlgType flags,
  Pointer userdata,
  SeqPortStreamProc proc
)

{
  return SeqPortStreamSetup (NULL, 0, 0, 0, NULL, lit, flags, userdata, proc);
}

/*******************************************************************************
*    
*   StreamCacheSetup (bsp, slp, flags, scp)
*   StreamCacheGetResidue (scp)
*   StreamCacheSetPosition (scp, pos)
*       SeqPort functional replacement implemented on top of SeqPortStreams
*
********************************************************************************/

NLM_EXTERN Boolean StreamCacheSetup (
  BioseqPtr bsp,
  SeqLocPtr slp,
  StreamFlgType flags,
  StreamCache PNTR scp
)

{
  if (bsp == NULL && slp == NULL) return FALSE;
  if (scp == NULL) return FALSE;

  MemSet ((Pointer) scp, 0, sizeof (StreamCache));

  if (bsp != NULL) {
    scp->bsp = bsp;
    scp->length = bsp->length;
  } else {
    scp->slp = slp;
    scp->length = SeqLocLen (slp);
  }
  scp->flags = flags;

  return TRUE;
}

static Boolean StreamCacheRefreshBuffer (
  StreamCache PNTR scp
)

{
  Bioseq         bsq;
  Int4           count;
  StreamFlgType  flags;
  SeqLocPtr      loc;
  SeqLoc         sl;
  SeqLocPtr      slp;
  Int4           stop;

  if (scp == NULL) return FALSE;

  if (scp->ctr >= scp->total) {
    scp->offset += (Int4) scp->total;
    scp->ctr = 0;
    scp->total = 0;

    MemSet ((Pointer) &(scp->buf), 0, sizeof (scp->buf));

    if (scp->offset < 0 || scp->offset >= scp->length) return FALSE;

    stop = MIN (scp->offset + 4000L, scp->length);

    flags = scp->flags;
    if ((flags & STREAM_GAP_MASK) == GAP_TO_SINGLE_DASH || (flags & STREAM_GAP_MASK) == 0) {
      /* if expand_gaps_to_dashes not equal to gaps_to_single_dash + stream_gap_mask, need to clear other bits first */
      flags |= EXPAND_GAPS_TO_DASHES;
    }
    if ((flags & SUPPRESS_VIRT_SEQ) != 0) {
      flags ^= SUPPRESS_VIRT_SEQ;
      flags |= STREAM_VIRT_AS_PLUS;
    }

    if (scp->bsp != NULL) {

      count = SeqPortStreamInt (scp->bsp, scp->offset, stop - 1, Seq_strand_plus,
                                flags, (Pointer) &(scp->buf), NULL);
      if (count < 0) {
        scp->failed = TRUE;
      }

    } else if (scp->slp != NULL) {

      slp = scp->slp;
      MemSet ((Pointer) &bsq, 0, sizeof (Bioseq));
      MemSet ((Pointer) &sl, 0, sizeof (SeqLoc));
      bsq.repr = Seq_repr_seg;
      bsq.mol = Seq_mol_na;
      bsq.seq_ext_type = 1;
      bsq.length = SeqLocLen (slp);
      bsq.seq_ext = &sl;
      if (slp->choice == SEQLOC_MIX || slp->choice == SEQLOC_PACKED_INT) {
        loc = (SeqLocPtr) slp->data.ptrvalue;
        if (loc != NULL) {
          sl.choice = loc->choice;
          sl.data.ptrvalue = (Pointer) loc->data.ptrvalue;
          sl.next = loc->next;
        }
      } else {
        sl.choice = slp->choice;
        sl.data.ptrvalue = (Pointer) slp->data.ptrvalue;
        sl.next = NULL;
      }

      SeqPortStreamInt (&bsq, scp->offset, stop - 1, Seq_strand_plus,
                        flags, (Pointer) &(scp->buf), NULL);
    }

    scp->total = StringLen (scp->buf);
  }

  return TRUE;
}

NLM_EXTERN Uint1 StreamCacheGetResidue (
  StreamCache PNTR scp
)

{
  Uint1  residue = '\0';

  if (scp == NULL) return '\0';

  if (scp->ctr >= scp->total) {
    if (! StreamCacheRefreshBuffer (scp)) return '\0';
  }

  if (scp->ctr < scp->total) {
    residue = scp->buf [(int) scp->ctr];
    (scp->ctr)++;

    if (residue == '-') {
    
      if ((scp->flags & STREAM_GAP_MASK) == 0) {
        while (residue == '-') {
          if (scp->ctr >= scp->total) {
            if (! StreamCacheRefreshBuffer (scp)) return '\0';
          }

          while (scp->ctr < scp->total && residue == '-') {
            residue = scp->buf [(int) scp->ctr];
            (scp->ctr)++;
          }
        }
        if (residue == '-') return '\0';

      } else if ((scp->flags & STREAM_GAP_MASK) == GAP_TO_SINGLE_DASH) {

        while (residue == '-') {
          if (scp->ctr >= scp->total) {
            if (! StreamCacheRefreshBuffer (scp)) return '-';
          }

          while (scp->ctr < scp->total && residue == '-') {
            residue = scp->buf [(int) scp->ctr];
            if (residue != '-') return '-';
            (scp->ctr)++;
          }
        }
      }

    } else if (residue == '+') {

      if ((scp->flags & SUPPRESS_VIRT_SEQ) != 0) {
        while (residue == '+') {
          if (scp->ctr >= scp->total) {
            if (! StreamCacheRefreshBuffer (scp)) return '\0';
          }

          while (scp->ctr < scp->total && residue == '+') {
            residue = scp->buf [(int) scp->ctr];
            (scp->ctr)++;
          }
        }
        if (residue == '+') return '\0';
      }
    }
  }

  return residue;
}

NLM_EXTERN Boolean StreamCacheSetPosition (
  StreamCache PNTR scp,
  Int4 pos
)

{
  if (scp == NULL) return FALSE;

  if (scp->offset <= pos && scp->offset + (Int4) scp->total >= pos) {
    scp->ctr = (Int2) (pos - scp->offset);
    return TRUE;
  }

  scp->ctr = 0;
  scp->total = 0;
  scp->offset = pos;

  if (scp->offset < 0 || scp->offset >= scp->length) {
    scp->offset = 0;
    return FALSE;
  }

  return TRUE;
}

/*******************************************************************************
*    
*    ProteinFromCdRegionEx ( SeqFeatPtr sfp, Boolean include_stop, Boolean remove_trailingX)
*        replacement for old ProteinFromCdRegionEx, but using TransTableTranslateCdRegion. 
*
********************************************************************************/

NLM_EXTERN ByteStorePtr ProteinFromCdRegionExEx (SeqFeatPtr sfp, Boolean include_stop, Boolean remove_trailingX, BoolPtr altStartP, Boolean farProdFetchOK)

{
  ByteStorePtr   bs;
  CdRegionPtr    crp;
  Int2           genCode = 0;
  Char           str [32];
  Boolean        tableExists = FALSE;
  TransTablePtr  tbl = NULL;
  ValNodePtr     vnp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return NULL;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL) return NULL;

  /* find genetic code */

  if (crp->genetic_code != NULL) {
    vnp = (ValNodePtr) crp->genetic_code->data.ptrvalue;
    while (vnp != NULL) {
      if (vnp->choice == 2) {
        genCode = (Int2) vnp->data.intvalue;
      }
      vnp = vnp->next;
    }
  }

  if (genCode == 7) {
    genCode = 4;
  } else if (genCode == 8) {
    genCode = 1;
  } else if (genCode == 0) {
    genCode = 1;
  }

  /* set app property name for storing desired FSA */

  sprintf (str, "TransTableFSAforGenCode%d", (int) genCode);

  /* get FSA for desired genetic code if it already exists */

  tbl = (TransTablePtr) GetAppProperty (str);
  tableExists = (Boolean) (tbl != NULL);

  bs = TransTableTranslateCdRegionEx (&tbl, sfp, include_stop, remove_trailingX,
                                      FALSE, altStartP, farProdFetchOK);

  /* save FSA in genetic code-specific app property name */

  if (! tableExists) {
    SetAppProperty (str, (Pointer) tbl);
  }

  return bs;
}

NLM_EXTERN ByteStorePtr ProteinFromCdRegionEx (SeqFeatPtr sfp, Boolean include_stop, Boolean remove_trailingX)

{
  return ProteinFromCdRegionExEx (sfp, include_stop, remove_trailingX, NULL, TRUE);
}

NLM_EXTERN ByteStorePtr ProteinFromCdRegionExWithTrailingCodonHandling
(
  SeqFeatPtr sfp,
  Boolean include_stop,
  Boolean remove_trailingX,
  Boolean no_stop_at_end_of_complete_cds
)

{
  ByteStorePtr   bs;
  CdRegionPtr    crp;
  Int2           genCode = 0;
  Char           str [32];
  Boolean        tableExists = FALSE;
  TransTablePtr  tbl = NULL;
  ValNodePtr     vnp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return NULL;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL) return NULL;

  /* find genetic code */

  if (crp->genetic_code != NULL) {
    vnp = (ValNodePtr) crp->genetic_code->data.ptrvalue;
    while (vnp != NULL) {
      if (vnp->choice == 2) {
        genCode = (Int2) vnp->data.intvalue;
      }
      vnp = vnp->next;
    }
  }

  if (genCode == 7) {
    genCode = 4;
  } else if (genCode == 8) {
    genCode = 1;
  } else if (genCode == 0) {
    genCode = 1;
  }

  /* set app property name for storing desired FSA */

  sprintf (str, "TransTableFSAforGenCode%d", (int) genCode);

  /* get FSA for desired genetic code if it already exists */

  tbl = (TransTablePtr) GetAppProperty (str);
  tableExists = (Boolean) (tbl != NULL);

  bs = TransTableTranslateCdRegion (&tbl, sfp, include_stop, remove_trailingX,
                                    no_stop_at_end_of_complete_cds);

  /* save FSA in genetic code-specific app property name */

  if (! tableExists) {
    SetAppProperty (str, (Pointer) tbl);
  }

  return bs;
}

NLM_EXTERN Uint1 AAForCodon (Uint1Ptr codon, CharPtr codes);

/*****************************************************************************
*
*   ProteinFromCdRegion(sfp, include_stop)
*       produces a ByteStorePtr containing the protein sequence in
*   ncbieaa code for the CdRegion sfp.  If include_stop, will translate
*   through stop codons.  If NOT include_stop, will stop at first stop
*   codon and return the protein sequence NOT including the terminating
*   stop.  Supports reading frame, alternate genetic codes, and code breaks
*   in the CdRegion. Removes trailing "X" from partial translation.
*
*****************************************************************************/
NLM_EXTERN ByteStorePtr ProteinFromCdRegion(SeqFeatPtr sfp, Boolean include_stop)
{
    return ProteinFromCdRegionEx(sfp, include_stop, TRUE);
}


/* old version of ProteinFromCdRegionEx no longer compiled (below) */

#if 0
/*******************************************************************************
*    
*    ProteinFromCdRegionEx( SeqFeatPtr sfp, Boolean include_stop, Boolean remove_trailingX)
*        same behavior as ProteinFromCdRegion, but another Boolean remove_trailingX
*    specifies whether trailing X's should be removed. 
*
********************************************************************************/

NLM_EXTERN ByteStorePtr Old_ProteinFromCdRegionEx (SeqFeatPtr sfp, Boolean include_stop, Boolean remove_trailingX)
{
    SeqPortPtr spp = NULL;
    ByteStorePtr bs = NULL;
    Uint1 residue = 0;
    Int4 pos1, pos2, pos, len;
    Int4Ptr the_breaks = NULL;
    Uint1Ptr the_residues = NULL;
    Int2 num_code_break = 0, use_break;
    SeqLocPtr tmp;
    Int2 i;
    Uint1 codon[3], aa;
    CdRegionPtr crp;
    ValNodePtr vnp;
    GeneticCodePtr gcp;
    CharPtr vals, codes;
    CodeBreakPtr cbp;
    Boolean bad_base, no_start, check_start, got_stop;
    Uint2 part_prod = 0, part_loc = 0;
    Boolean incompleteLastCodon;

    if ((sfp == NULL) || (sfp->data.choice != 3))
        return NULL;

    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    len = SeqLocLen(sfp->location);

    num_code_break = 0;
    if (crp->code_break != NULL)
    {
        cbp = crp->code_break;
        while (cbp != NULL)
        {
            num_code_break++;
            cbp = cbp->next;
        }
        the_breaks = (Int4Ptr) MemNew((size_t)(num_code_break * sizeof(Int4)));
        the_residues = (Uint1Ptr) MemNew((size_t)(num_code_break * sizeof(Uint1)));

        num_code_break = 0;
        cbp = crp->code_break;
        while (cbp != NULL)
        {
            pos1 = INT4_MAX;
            pos2 = -10;
            tmp = NULL;
            while ((tmp = SeqLocFindNext(cbp->loc, tmp)) != NULL)
            {
                pos = GetOffsetInLoc(tmp, sfp->location, 
SEQLOC_START);
                if (pos < pos1)
                    pos1 = pos;
                pos = GetOffsetInLoc(tmp, sfp->location, 
SEQLOC_STOP);
                if (pos > pos2)
                    pos2 = pos;
            }
            pos = pos2 - pos1; /* codon length */
            if (pos == 2 || (pos >= 0 && pos <= 1 && pos2 == len - 1))   /*  a codon */
            /* allowing a partial codon at the end */
            {
                the_breaks[num_code_break] = pos1;
                the_residues[num_code_break] = (Uint1) 
cbp->aa.value.intvalue;
                num_code_break++;
            }
            else
            {
                ErrPost(CTX_NCBIOBJ, 1, "Invalid Code-break.loc");
            }

            cbp = cbp->next;
        }
    }

    gcp = NULL;
    if (crp->genetic_code != NULL)
    {
        vnp = (ValNodePtr)(crp->genetic_code->data.ptrvalue);
        while ((vnp != NULL) && (gcp == NULL))
        {
            switch (vnp->choice)
            {
            case 1:   /* name */
                gcp = GeneticCodeFind(0, 
(CharPtr)vnp->data.ptrvalue);
                break;
            case 2:   /* id */
                gcp = GeneticCodeFind(vnp->data.intvalue, NULL);
                break;
            case 3:   /* ncbieaa */
            case 6:   /* sncbieaa */
            case 4:   /* ncbi8aa */
            case 5:      /* ncbistdaa */
            case 7:   /* sncbi8aa */
            case 8:   /* sncbistdaa */
            default:
                break;
            }
            vnp = vnp->next;
        }
    }
    if (gcp == NULL)
        gcp = GeneticCodeFind(1, NULL);   /* use universal */
    if (gcp == NULL)
        goto erret;

    vals = NULL;
    codes = NULL;
    for (vnp = (ValNodePtr)gcp->data.ptrvalue; vnp != NULL; vnp = vnp->next)
    {
        if (vnp->choice == 6)   /* sncbieaa */
            vals = (CharPtr)vnp->data.ptrvalue;
        else if (vnp->choice == 3)  /* ncbieaa */
            codes = (CharPtr)vnp->data.ptrvalue;
    }
    if (codes == NULL)
        goto erret;

    no_start = FALSE;
    part_loc = SeqLocPartialCheck(sfp->location);
    part_prod = SeqLocPartialCheck(sfp->product);
    if ((part_loc & SLP_START) || (part_prod & SLP_START))
        no_start = TRUE;

    if ((vals == NULL) || (no_start) || (crp->frame > 1))  /* no special 
starts */
    {
        vals = codes;
        check_start = FALSE;
    }
    else
        check_start = TRUE;

    spp = SeqPortNewByLoc(sfp->location, Seq_code_ncbi4na);
    if (spp == NULL)
        goto erret;

    /* len = SeqLocLen(sfp->location); - saved above */    /* size of coding region */
    len /= 3;                           /* size of 
protein */
    len += 1;                           /* allow 
partial codon at end */
    bs = BSNew(len);
    if (bs == NULL)
        goto erret;

    if (crp->frame == 2)     /* skip partial first codon */
        pos = 1;
    else if (crp->frame == 3)
        pos = 2;
    else
        pos = 0;
    SeqPortSeek(spp, pos, SEEK_SET);
    got_stop = FALSE;

    incompleteLastCodon = FALSE;

    do
    {
        use_break = -1;
        for (i = 0; i < num_code_break; i++)
        {
            if (pos == the_breaks[i])
            {
                use_break = i;
                i = num_code_break;
            }
        }

        bad_base = FALSE;
        for (i = 0; i < 3; i++)
        {
            residue = SeqPortGetResidue(spp);
            if (residue == SEQPORT_VIRT || residue == SEQPORT_EOS) {
                /* skip past null NULL in seqport, get next - JK */
                residue = SeqPortGetResidue(spp);
            }
            if (residue == SEQPORT_EOF)
                break;
            if (residue == INVALID_RESIDUE)
                bad_base = TRUE;
            codon[i] = residue;
        }
        if (! i)   /* no bases */
            break;
        while (i < 3)      /* incomplete last codon */
        {
            codon[i] = 15;   /* N */
            i++;
            incompleteLastCodon = TRUE;
        }

        pos += 3;
        if (use_break >= 0)
            aa = the_residues[use_break];
        else if (bad_base)
            aa = 'X';
        else
        {
            aa = AAForCodon(codon, vals);
            if (check_start)   /* first codon on possibly complete 
CDS */
            {
                if (aa == '-')   /* invalid start */
                {
                    /* if no explict partial at either end, but 
feature is */
                    /* annotated as partial, then guess should 
use internal */
                    /* amino acid code */

                    if ((! ((part_loc & SLP_STOP) || 
(part_prod & SLP_STOP))) &&
                        (sfp->partial))
                        aa = AAForCodon(codon, codes);  
/* get internal aa */
                }
                check_start = FALSE;
            }
        }

        if ((! include_stop) && (aa == '*'))
        {
            got_stop = TRUE;
            break;
        }

        BSPutByte(bs, (Int2)aa);

        vals = codes;     /* not a start codon anymore */

    } while (residue != SEQPORT_EOF);

    if ((! got_stop) && incompleteLastCodon) {
        BSSeek(bs, -1, SEEK_END);  /* remove last X if incomplete last codon */
        aa = (Uint1)BSGetByte(bs);
        if ((aa == 'X') && (BSLen(bs)))
        {
            BSSeek(bs, -1, SEEK_END);
            BSDelete(bs, 1);
            BSSeek(bs, -1, SEEK_END);
        }
    }
    if ((! got_stop) && remove_trailingX)   /* only remove trailing X on partial CDS */
    {
        BSSeek(bs, -1, SEEK_END);  /* back up to last residue */
        aa = (Uint1)BSGetByte(bs);
        while ((aa == 'X') && (BSLen(bs)))
        {
            BSSeek(bs, -1, SEEK_END);
            BSDelete(bs, 1);
            BSSeek(bs, -1, SEEK_END);
            aa = (Uint1)BSGetByte(bs);
        }
    }

    if (! BSLen(bs)) goto erret;

ret:
    SeqPortFree(spp);
    MemFree(the_breaks);
    MemFree(the_residues);
    return bs;
erret:
    bs = BSFree(bs);
    goto ret;
}
#endif

/* old version of ProteinFromCdRegionEx no longer compiled (above) */


/*****************************************************************************
*
*   Uint1 AAForCodon (Uint1Ptr codon, CharPtr codes)
*       codon is 3 values in ncbi4na code
*       codes is the geneic code array to use
*          MUST have 'X' as unknown amino acid
*
*****************************************************************************/
NLM_EXTERN Uint1 AAForCodon (Uint1Ptr codon, CharPtr codes)
{
    register Uint1 aa = 0, taa;
    register int i, j, k, index0, index1, index2;
    static Uint1 mapping[4] = { 8,     /* T in ncbi4na */
                                2,     /* C */
                                1,     /* A */
                                4 };   /* G */


    for (i = 0; i < 4; i++)
    {
        if (codon[0] & mapping[i])
        {
            index0 = i * 16;
            for (j = 0; j < 4; j++)
            {
                if (codon[1] & mapping[j])
                {
                    index1 = index0 + (j * 4);
                    for (k = 0; k < 4; k++)
                    {
                        if (codon[2] & mapping[k])
                        {
                            index2 = index1 + k;
                            taa = codes[index2];
                            if (! aa)
                                aa = taa;
                            else
                            {
                                if (taa != aa)
                                {
                                    aa = 
'X';
                                    break;
                                }
                            }
                        }
                        if (aa == 'X')
                            break;
                    }
                }
                if (aa == 'X')
                    break;
            }
        }
        if (aa == 'X')
            break;
    }
    return aa;
}

static    Uint1 codon_xref [4] = {   /* mapping from NCBI2na to codon codes */
        2,  /* A */
        1,  /* C */
        3,  /* G */
        0 }; /* T */

/*****************************************************************************
*
*   Uint1 IndexForCodon (codon, code)
*       returns index into genetic codes codon array, give 3 bases of the
*       codon in any alphabet
*       returns INVALID_RESIDUE on failure
*   
*****************************************************************************/
NLM_EXTERN Uint1 IndexForCodon (Uint1Ptr codon, Uint1 code)
{
    Int2 i, j;
    SeqMapTablePtr smtp;
    Uint1 residue, index = 0;

    smtp = SeqMapTableFind(Seq_code_ncbi2na, code);
    if (smtp == NULL) return INVALID_RESIDUE;

    for (i=0, j=16; i < 3; i++, j /= 4)
    {
        residue = SeqMapTableConvert(smtp, codon[i]);
        if (residue > 3) return INVALID_RESIDUE;
        residue = codon_xref[residue];
        index += (Uint1)(residue * j);
    }

    return index;
}

/*****************************************************************************
*
*   Boolean CodonForIndex (index, code, codon)
*       Fills codon (3 Uint1 array) with codon corresponding to index,
*       in sequence alphabet code.
*       Index is the Genetic code index.
*       returns TRUE on success.
*
*****************************************************************************/
NLM_EXTERN Boolean CodonForIndex (Uint1 index, Uint1 code, Uint1Ptr codon)
{
    Int2 i, j, k;
    SeqMapTablePtr smtp;
    Uint1 residue;
    
    if (codon == NULL) return FALSE;
    if (index > 63) return FALSE;

    smtp = SeqMapTableFind(code, Seq_code_ncbi2na);
    if (smtp == NULL) return FALSE;

    for (i = 0, j = 16; i < 3; i++, j /= 4)
    {
        residue = (Uint1)((Int2)index / j);
        index -= (Uint1)(residue * j);
        for (k = 0; k < 4; k++)
        {
            if (codon_xref[k] == residue)
            {
                residue = (Uint1)k;
                break;
            }
        }
        residue = SeqMapTableConvert(smtp, residue);
        codon[i] = residue;
    }

    return TRUE;
}

/*----------- GetFrameFromLoc()-----------------*/

/*****************************************************************************
*
*   Int2 GetFrameFromLoc (slp)
*       returns 1,2,3 if can find the frame
*       0 if not
*
*****************************************************************************/
NLM_EXTERN Uint1 GetFrameFromLoc (SeqLocPtr slp)
{
    Uint1 frame = 0;
    SeqLocPtr curr, last;
    Boolean is_partial;
    SeqIntPtr sip;
    SeqPntPtr spp;

    if (slp == NULL)
        return frame;

    curr = SeqLocFindNext(slp, NULL);

    is_partial = FALSE;
    switch (curr->choice)
    {
        case SEQLOC_INT:
            sip = (SeqIntPtr)curr->data.ptrvalue;
            if (sip->strand == Seq_strand_minus)
            {
                if (sip->if_to != NULL)
                    is_partial = TRUE;
            }
            else if (sip->if_from != NULL)
                is_partial = TRUE;
            break;
        case SEQLOC_PNT:
            spp = (SeqPntPtr)curr->data.ptrvalue;
            if (spp->fuzz != NULL)
                is_partial = TRUE;
            break;
        default:
            return frame;
    }
        

    if (! is_partial)
        return (Int2) 1;    /* complete 5' end, it's frame 1 */

    is_partial = FALSE;
    last = curr;
    while ((curr = SeqLocFindNext(slp, last)) != NULL)
        last = curr;

    switch (last->choice)
    {
        case SEQLOC_INT:
            sip = (SeqIntPtr) last->data.ptrvalue;
            if (sip->strand == Seq_strand_minus)
            {
                if (sip->if_from != NULL)
                    return frame;
            }
            else if (sip->if_to != NULL)
                return frame;
            break;
        case SEQLOC_PNT:
            spp = (SeqPntPtr) last->data.ptrvalue;
            if (spp->fuzz != NULL)
                return frame;
            break;
        default:
            return frame;
    }

                      /* have complete last codon, get frame 
from length */
    frame = (Uint1)(SeqLocLen(slp) % 3);
    if (frame == 0)
        frame = 1;
    else if (frame == 1)
        frame = 2;
    else
        frame = 3;

    return frame;
}

static Boolean add_fuzziness_to_loc (SeqLocPtr slp, Boolean less)
{
    IntFuzzPtr ifp;
    SeqIntPtr sint;
    SeqPntPtr spnt;    

    sint = NULL;
    spnt = NULL;

    if(slp->choice == SEQLOC_INT)
        sint = (SeqIntPtr) slp->data.ptrvalue;
    else
    {
        if(slp->choice == SEQLOC_PNT)
            spnt = (SeqPntPtr) slp->data.ptrvalue;
        else
            return FALSE;
    }
    ifp = IntFuzzNew();
    ifp->choice = 4;
    ifp->a = less ? 2 : 1;

    if(spnt != NULL)
        spnt->fuzz = ifp;
    else if (sint != NULL)
    {
        if(less)
            sint->if_from = ifp;
        else
            sint->if_to = ifp;
    }

    return TRUE;
}


static Boolean load_fuzz_to_DNA(SeqLocPtr dnaLoc, SeqLocPtr aaLoc, Boolean 
first)
{
    Uint1 strand;
    SeqPntPtr spnt;
    SeqIntPtr sint;
    IntFuzzPtr ifp;
    Boolean load, less;

    load = FALSE;
    strand = SeqLocStrand(aaLoc);
    if(aaLoc->choice == SEQLOC_INT)
    {
        sint = (SeqIntPtr) aaLoc->data.ptrvalue;
        if((first && strand != Seq_strand_minus ) || 
            (!first && strand == Seq_strand_minus))    /*the first 
Seq-loc*/
        {
            ifp = sint->if_from;
            if(ifp && ifp->choice == 4 )
                load = (ifp->a == 2);
        }
        else
        {
            ifp = sint->if_to;
            if(ifp && ifp->choice == 4)
                load = (ifp->a == 1);
        }
    }
    else if(aaLoc->choice == SEQLOC_PNT)
    {
        spnt = (SeqPntPtr) aaLoc->data.ptrvalue;
        ifp = spnt->fuzz;
        if(ifp && ifp->choice == 4)
        {
            if(first)
                load = (ifp->a == 2);
            else
                load = (ifp->a == 1);
        }
    }

    if(load)
    {
        if(SeqLocStrand(dnaLoc) == Seq_strand_minus)
            less = (first == FALSE);
        else
            less = first;
        add_fuzziness_to_loc (dnaLoc, less);
        return TRUE;
    }
    else
        return FALSE;
}    

/******************************************************************
*
*    aaLoc_to_dnaLoc(sfp, aa_loc)
*    map a SeqLoc on the amino acid sequence
*       to a Seq-loc in the    DNA sequence
*       through a CdRegion feature
*       
*       This now calls the more general productLoc_to_locationLoc(sfp, productLoc)
*
******************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL aaLoc_to_dnaLoc(SeqFeatPtr sfp, SeqLocPtr aa_loc)
{
    return productLoc_to_locationLoc(sfp, aa_loc);
}

/******************************************************************
*
*    aaLoc_to_dnaLoc(sfp, productLoc)
*    map a SeqLoc on the product sequence
*       to a Seq-loc in the    location sequence
*       through a feature.
*
*       if the feature is a CdRegion, converts by modulo 3
*       to support aaLoc_to_dnaLoc() function
*
******************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL productLoc_to_locationLoc(SeqFeatPtr sfp, SeqLocPtr productLoc)
{
    SeqLocPtr head = NULL, slp, tmp, next;
    Int4 product_start, product_stop;
    SeqBondPtr sbp;
    ValNode vn;
    Boolean is_cdregion = FALSE;
  Boolean partial5, partial3;

    if ((sfp == NULL) || (productLoc == NULL)) return head;
    if (sfp->data.choice == 3) is_cdregion = TRUE;
    if (sfp->product == NULL) return head;
    if (! (SeqIdForSameBioseq(SeqLocId(productLoc), SeqLocId(sfp->product))))
        return head;

    if (productLoc->choice == SEQLOC_BOND)   /* fake this one in */
    {
        sbp = (SeqBondPtr)(productLoc->data.ptrvalue);
        tmp = productInterval_to_locationIntervals(sfp, sbp->a->point, sbp->a->point, FALSE);
        if (sbp->b == NULL)  /* one point in bond */
            return tmp;

        SeqLocAdd(&head, tmp, TRUE, FALSE);
        tmp = productInterval_to_locationIntervals(sfp, sbp->b->point, sbp->b->point, FALSE);
        if (tmp == NULL)
            return head;

        vn.choice = SEQLOC_NULL;  /* make a mix with an internal NULL */
        vn.next = NULL;
        vn.data.ptrvalue = NULL;

        SeqLocAdd(&head, &vn, TRUE, TRUE);  /* copy it in */
        SeqLocAdd(&head, tmp, TRUE, FALSE); /* put real 3 base int in */

        goto ret;
    }

  CheckSeqLocForPartial (productLoc, &partial5, &partial3);
    slp = NULL;
    while ((slp = SeqLocFindNext(productLoc, slp)) != NULL)
    {
        product_start = SeqLocStart(slp);
        product_stop = SeqLocStop(slp);
        if ((product_start >= 0) && (product_stop >= 0))
        {
           tmp = productInterval_to_locationIntervals(sfp, product_start, product_stop, partial5);
           if(tmp != NULL)
            load_fuzz_to_DNA(tmp, slp, TRUE);
           while (tmp != NULL)
           {
               next = tmp->next;
               tmp->next = NULL;
               if(next == NULL)
                load_fuzz_to_DNA(tmp, slp, FALSE);
               SeqLocAdd(&head, tmp, TRUE, FALSE);
               tmp = next;
           }
        } else if (slp->choice == SEQLOC_NULL) {
            vn.choice = SEQLOC_NULL;  /* make a mix with an internal NULL */
            vn.next = NULL;
            vn.data.ptrvalue = NULL;
            SeqLocAdd(&head, &vn, TRUE, TRUE); 
        }
    }
ret:               
    return SeqLocPackage(head);
}

/******************************************************************
*
*       aaFeatLoc_to_dnaFeatLoc(sfp, aa_loc)
*       map a SeqLoc on the amino acid sequence
*       to a Seq-loc in the     DNA sequence
*       through a CdRegion feature
*
*       uses aaLoc_to_dnaLoc() but does additional checks to
*       extend dnaLoc at either end to compensate for positions in
*       the dna which do not corresspond to the amino acid sequence
*       (partial codons which are not translated).
*
******************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL aaFeatLoc_to_dnaFeatLoc(SeqFeatPtr sfp,
                                                     SeqLocPtr aa_loc)
{
    SeqLocPtr dnaLoc = NULL;
    Uint2 dnaPartial;
    Int4 aaPos;
    SeqLocPtr tmp1 = NULL, tmp2 = NULL, tmp;
    SeqIdPtr sip;
    CdRegionPtr crp;
    SeqIntPtr sp1, sp2;
    BioseqPtr bsp;
  Boolean   aa_partialn, aa_partialc;

    dnaLoc = aaLoc_to_dnaLoc(sfp, aa_loc);
    if (dnaLoc == NULL) return dnaLoc;

    if (! sfp->partial)  /* no partial checks needed */
        return dnaLoc;


  CheckSeqLocForPartial (aa_loc, &aa_partialn, &aa_partialc);
    crp = (CdRegionPtr)(sfp->data.value.ptrvalue);

    aaPos = SeqLocStart(aa_loc);
    if ((! aaPos) && (crp->frame > 1) && aa_partialn)   /* using first amino acid */
    {
        tmp1 = SeqLocFindNext(sfp->location, NULL);
        tmp2 = SeqLocFindNext(dnaLoc, NULL);

        if ((tmp1->choice == SEQLOC_INT) &&
                         (tmp2->choice == SEQLOC_INT))
        {
            sp1 = (SeqIntPtr)(tmp1->data.ptrvalue);
            sp2 = (SeqIntPtr)(tmp2->data.ptrvalue);
            if (sp1->strand ==  Seq_strand_minus)
            {
                sp2->to = sp1->to;  /* add partial codon */
            }
            else
            {
                sp2->from = sp1->from;
            }
        }
    }

    dnaPartial = SeqLocPartialCheck(sfp->location);
    if ((dnaPartial & SLP_STOP) && aa_partialc)   /* missing 3' end of cdregion */
    {
        sip = SeqLocId(aa_loc);
        bsp = BioseqFindCore(sip);
        if (bsp != NULL)
        {
            aaPos = SeqLocStop(aa_loc);
            if (aaPos == (bsp->length - 1)) /* last amino acid */
            {
                tmp = NULL;
                while ((tmp = SeqLocFindNext(sfp->location,tmp)) != NULL)
                {
                    tmp1 = tmp;
                }
                tmp = NULL;
                while ((tmp = SeqLocFindNext(dnaLoc,tmp)) != NULL)
                {
                    tmp2 = tmp;
                }
            
                if (tmp1 != NULL && tmp2 != NULL && (tmp1->choice == SEQLOC_INT) &&
                    (tmp2->choice == SEQLOC_INT))
                {
                    sp1 = (SeqIntPtr)(tmp1->data.ptrvalue);
                    sp2 = (SeqIntPtr)(tmp2->data.ptrvalue);
                    if (sp1->strand ==  Seq_strand_minus)
                    {
                        sp2->from = sp1->from;  /* add partial codon */
                    }
                    else
                    {
                        sp2->to = sp1->to;
                    }
                }
            }
    
        }
    }
    return dnaLoc;
}


static SeqLocPtr 
NucLocFromProtInterval 
(SeqFeatPtr cds, 
 Int4 prot_start, 
 Int4 prot_stop, 
 Boolean n_partial)
{
  CdRegionPtr crp;
  Int4        aa_before = 0, nt_this, prev_nt = 0, part_codon;
  SeqLocPtr   result = NULL;
  SeqLocPtr   slp = NULL; /* used for iterating through locations in the coding region */
  SeqLocPtr   loc; /* used for creating interval on NT sequence */
  Boolean     first_loc = TRUE;
  Int4        cds_int_start, cds_int_stop, cds_int_len;
  Int4        frame_start = 0;
  Int4        aa_int_start = 0, aa_int_stop = 0, aa_len, this_aa, aa_needed, aa_unneeded, aa_accumulated = 0;
  Int4        aa_from_this_interval;
  Uint1       strand;

  if (cds == NULL || cds->data.choice != SEQFEAT_CDREGION || prot_start < 0 || prot_stop < prot_start) {
    return NULL;
  }

  crp = (CdRegionPtr) cds->data.value.ptrvalue;
  if (crp == NULL) {
    return NULL;
  }
  if (crp->frame > 1) {
    frame_start = crp->frame - 1;
  }

  aa_len = prot_stop - prot_start + 1;

  while((slp = SeqLocFindNext(cds->location, slp)) != NULL) {
    cds_int_len = SeqLocLen (slp);
    cds_int_start = SeqLocStart (slp);
    cds_int_stop = SeqLocStop (slp);
    strand = SeqLocStrand (slp);

    if (first_loc) {
      if (strand == Seq_strand_minus) {
        cds_int_stop -= frame_start;
      } else {
        cds_int_start += frame_start;
      }
      cds_int_len -= frame_start;
    }

    /* calculate the number of NT that "count" for this interval -
      * don't include the NT in a partial codon at the beginning of
      * of the feature, but do include NT from a partial codon at
      * the end of the previous interval.
      */
    nt_this = cds_int_len + prev_nt;
    part_codon = nt_this % 3;
    nt_this -= part_codon;
    
    /* calculate how many AA are covered by this interval */
    this_aa = nt_this / 3;

    if (aa_before + this_aa >= prot_start) {

      /* figure out whether to take all of this interval, or just part of it */
      aa_from_this_interval = this_aa;

      /* 5' end (left for plus strand, right for minus) */
      if (aa_before < prot_start) {
        /* skip some at the beginning */
        aa_unneeded = prot_start - aa_before;
        aa_from_this_interval -= aa_unneeded;

        if (strand == Seq_strand_minus) {
          aa_int_stop = cds_int_stop + prev_nt - (3 * aa_unneeded);
        } else {
          aa_int_start = cds_int_start - prev_nt + (3 * aa_unneeded);
        }          
      } else {
        /* start at the beginning */
        if (strand == Seq_strand_minus) {
          aa_int_stop = cds_int_stop;
          if (first_loc) {
            if (n_partial) {
              /* put frame shift back in, if first loc and n-partial */
              aa_int_stop += frame_start;
            } else if (aa_before == prot_start) {
              /* starts in this interval, but after "remainder" of previous codon */
              aa_int_stop -= prev_nt;
            }
          }
        } else {
          aa_int_start = cds_int_start;
          if (first_loc) {
            if (n_partial) {
              /* put frame shift back in, if first loc and n-partial */
              aa_int_start -= frame_start;
            } else if (aa_before == prot_start) {
              /* starts in this interval, but after "remainder" of previous codon */
              aa_int_start += prev_nt;
            }
          }
        }
      }
      
      /* 3' end (right for plus strand, left for minus) */
      if (aa_accumulated + aa_from_this_interval < aa_len) {
        if (strand == Seq_strand_minus) {
          aa_int_start = cds_int_start;
        } else {
          aa_int_stop = cds_int_stop;
        }
      } else {
        /* just take the part that we need */
        aa_needed = aa_len - aa_accumulated;
        aa_unneeded = aa_from_this_interval - aa_needed;

        if (strand == Seq_strand_minus) {
          aa_int_start = cds_int_start + part_codon + (3 * aa_unneeded);
        } else {
          aa_int_stop = cds_int_stop - part_codon - (3 * aa_unneeded);
        }
        aa_from_this_interval -= aa_unneeded;
      }

      /* note - if aa_int_start > aa_int_stop, that means we eliminated
       * both ends of the interval.
       */
      if (aa_int_start <= aa_int_stop) {
        /* aa_accumulated now includes the number of complete codons that have
        * been accounted for (not counting a partial codon at the end of this
        * interval, if any
        */
        aa_accumulated += aa_from_this_interval;

        /* add interval to result */
            loc = SeqLocIntNew(aa_int_start, aa_int_stop, strand, SeqLocId(slp));
            SeqLocAdd(&result, loc, TRUE, FALSE);
      }
    }

    first_loc = FALSE;
    aa_before += this_aa;
    prev_nt = part_codon;

    if (aa_before > prot_stop) {
      break;
    }
  }

  return result;  
}


static SeqLocPtr NaLocFromNaInterval (SeqFeatPtr sfp, Int4 product_start, Int4 product_stop)
{
  SeqLocPtr slp = NULL;
  SeqLocPtr location_loc, loc;            /*for the sfp.location location*/

  Boolean is_end;            /**is the end for process reached?**/
  Int4 p_start=0, p_stop=0;        /**product sequence start & stop in defined
                    corresponding sfp.product **/
  Int4 cur_pos;            /**current sfp.product sequence position in process**/
  Int4 product_len;        /**length of the sfp.product **/

  Int4 d_start, d_stop;        /*the start and the stop of the sfp.location sequence*/
  Int4 offset;            /*offset from the start of the current exon*/
  Int4 aa_len;
  Uint1 strand;
  Int4 p_end_pos;    /*the end of the product sequence in the current loc*/

  cur_pos= product_start;
  product_len = 0;
  is_end = FALSE;
  p_start = 0;
  slp = NULL;
  location_loc= NULL;
  while(!is_end && ((slp = SeqLocFindNext(sfp->location, slp))!=NULL))
  {
      product_len += SeqLocLen(slp);
    p_stop = product_len - 1;

      p_end_pos = p_stop;

      if(p_stop >= product_stop)
      {
        p_stop = product_stop;        /**check if the end is reached**/
        is_end = TRUE;
      }

      if(p_stop >= cur_pos)    /*get the exon*/
      {
        offset = cur_pos - p_start;

          strand = SeqLocStrand(slp);
          if(strand == Seq_strand_minus)
            d_start = SeqLocStop(slp) - offset;
          else
            d_start = SeqLocStart(slp) + offset;

          d_stop = d_start;

      aa_len = MIN(p_stop, product_stop) - cur_pos +1;

      if(strand == Seq_strand_minus)
          {
              if(aa_len >= 0)
              {
                  d_stop -= (aa_len - 1);
              }
              else
        {
                  ++d_stop;
        }
                                
              d_stop = MAX(d_stop, SeqLocStart(slp));
              loc = SeqLocIntNew(d_stop, d_start, strand, SeqLocId(slp));
          }
          else
          {
              if(aa_len >= 0)
              {
                  d_stop += (aa_len - 1);
              }
              else
                  --d_stop;
                  
              d_stop = MIN(d_stop, SeqLocStop(slp));
              loc = SeqLocIntNew(d_start, d_stop, strand, SeqLocId(slp));
          }
          SeqLocAdd(&location_loc, loc, TRUE, FALSE);

          cur_pos = p_stop+1;
    }

    p_start = p_stop +1;

  }/**end of while(slp && !is_end) **/

  return location_loc;
}

/******************************************************************
*
*    productInterval_to_locationIntervals(sfp, product_start, product_stop)
*    map the amino acid sequence to a chain of Seq-locs in the 
*    DNA sequence through a CdRegion feature
*
******************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL 
productInterval_to_locationIntervals
(SeqFeatPtr sfp,
 Int4 product_start,
 Int4 product_stop, 
 Boolean aa_partialn)
{

  if (sfp->data.choice == SEQFEAT_CDREGION) {
    return NucLocFromProtInterval (sfp, product_start, product_stop, aa_partialn);
  } else {
    return NaLocFromNaInterval (sfp, product_start, product_stop);
  }

}

static Boolean load_fuzz_to_DNA PROTO((SeqLocPtr dnaLoc, SeqLocPtr aaLoc, 
Boolean first));
/******************************************************************
*
*    dnaLoc_to_aaLoc(sfp, location_loc, merge)
*    map a SeqLoc on the DNA sequence
*       to a Seq-loc in the    protein sequence
*       through a CdRegion feature
*   if (merge) adjacent intervals on the amino acid sequence
*      are merged into one. This should be the usual case.
*
******************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL dnaLoc_to_aaLoc(SeqFeatPtr sfp, SeqLocPtr location_loc, Boolean 
merge, Int4Ptr frame, Boolean allowTerminator)
{
    SeqLocPtr aa_loc = NULL, loc;
    CdRegionPtr crp;
    Int4 product_len, end_pos, frame_offset;
    GatherRange gr;
    Int4 a_left = 0, a_right, last_aa = -20, aa_from, aa_to;
  Int4 cds_left, cds_right;
    SeqLocPtr slp;
    Int2 cmpval;
    SeqIdPtr aa_sip;
    BioseqPtr bsp;
  Boolean partial5, partial3;
  Uint1 strand;

    if ((sfp == NULL) || (location_loc == NULL)) return aa_loc;
    if (sfp->data.choice != 3) return aa_loc;
    if (sfp->product == NULL) return aa_loc;

    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if(crp == NULL) return aa_loc;

    /* location_loc must be equal or contained in feature */
    cmpval = SeqLocCompare(location_loc, sfp->location);
    if (! ((cmpval == SLC_A_IN_B) || (cmpval == SLC_A_EQ_B)))
        return aa_loc;

    aa_sip = SeqLocId(sfp->product);
    if (aa_sip == NULL) return aa_loc;
    bsp = BioseqLockById(aa_sip);
    if (bsp == NULL) return aa_loc;
    end_pos = bsp->length - 1;
    BioseqUnlock(bsp);

    if(crp->frame == 0)
        frame_offset = 0;
    else
        frame_offset = (Int4)crp->frame-1;

  cds_left = SeqLocStart (sfp->location);
  cds_right = SeqLocStop (sfp->location);


    slp = NULL;
    product_len = 0;
    loc = NULL;
    while ((slp = SeqLocFindNext(sfp->location, slp))!=NULL)
    {
    if (SeqLocOffset(location_loc, slp, &gr, 0))
      {
            SeqLocOffset(slp, location_loc, &gr, 0);

            a_left = gr.left + product_len;
            a_right = gr.right + product_len;
      if (frame_offset > 0) {
        a_left -= frame_offset;
        a_right -= frame_offset;
      }

            if (a_left < 0)
      {
        CheckSeqLocForPartial (slp, &partial5, &partial3);
        strand = SeqLocStrand (slp);
        if ((partial5 && strand != Seq_strand_minus) || (partial3 && strand == Seq_strand_minus)) {
                  a_left = gr.left;
        } else {
          a_left += 3;
        }
      }
      if (a_right > (bsp->length) * 3 - 1 && !allowTerminator) {
        CheckSeqLocForPartial (slp, &partial5, &partial3);
        strand = SeqLocStrand (slp);
        if (partial3 && a_right == bsp->length * 3) {
          /* it's ok, leave it alone */
        } else if ((partial5 && strand != Seq_strand_minus) || (partial3 && strand == Seq_strand_minus)) {
                  a_right = (bsp->length * 3) - 1;
        } else {
          a_right -= 3;
        }
      }

            aa_from = a_left / 3;
            aa_to = a_right / 3;

            if (aa_to > end_pos && !allowTerminator)
                aa_to = end_pos;

            if (merge)
            {
                if (aa_from <= last_aa)  /* overlap due to codons */
                    aa_from = last_aa+1;  /* set up to merge */
            }

      /* NOTE - if a_left is not <= a_right, then a correction for frame may have
       * caused the location to not actually be mappable to the protein sequence.
       */
            if ((aa_from <= aa_to || (allowTerminator && aa_from == aa_to + 1)) && a_left <= a_right)
            {
                if(loc != NULL)
                {
                    if(aa_loc == NULL)
                        load_fuzz_to_DNA(loc, location_loc, TRUE);
                    SeqLocAdd(&aa_loc, loc, merge, FALSE);
                }
                loc = SeqLocIntNew(aa_from, aa_to, 0, aa_sip);
                last_aa = aa_to;
            }
    }

    product_len += SeqLocLen(slp);        
    }

    if(loc != NULL)
    {
        if(aa_loc == NULL)
            load_fuzz_to_DNA(loc, location_loc, TRUE);
        load_fuzz_to_DNA(loc, location_loc, FALSE);
        SeqLocAdd(&aa_loc, loc, merge, FALSE);
    }
    if (frame != NULL)
        *frame = a_left % 3;

    return SeqLocPackage(aa_loc);
}

/*****************************************************************************
*
*   BioseqHash(bsp)
*       Computes a (almost) unique hash code for a bioseq
*
*****************************************************************************/
NLM_EXTERN Uint4 BioseqHash (BioseqPtr bsp)
{
    Uint4 hashval = 0;
    SeqPortPtr spp;
    Uint1 code;
    Int2 residue;

    if (bsp == NULL) return hashval;

    if (ISA_na(bsp->mol))
        code = Seq_code_iupacna;
    else
        code = Seq_code_ncbieaa;

    spp = SeqPortNew(bsp, 0, -1, 0, code);
    if (spp == NULL) return hashval;

    while ((residue = SeqPortGetResidue(spp)) != SEQPORT_EOF)
    {
        hashval *= 1103515245;
        hashval += (Uint4)residue + 12345;
    }

    SeqPortFree(spp);

    return hashval;
}


/*-------------- BioseqRevComp () ---------------------------*/
/***********************************************************************
*   BioseqRevComp:   Takes the nucleic acid sequence from Bioseq
*    Entry and gives the reverse complement sequence in place
*       Does not change features.
************************************************************************/
NLM_EXTERN Boolean LIBCALL BioseqRevComp (BioseqPtr bsp)
{
    Boolean retval;

    retval = BioseqReverse (bsp);
    if (retval)
        retval = BioseqComplement(bsp);
    return retval;
}

NLM_EXTERN Boolean ComplementSeqData (Uint1 seqtype, Int4 seqlen, SeqDataPtr sdp)
{
    SeqCodeTablePtr sctp;
    ByteStorePtr    bysp;
    long            readbyte, bslen;
    Uint1           byte = 0, byte_to, newbyte = 0, residue;
    Uint1           comp, bitctr, mask, lshift, rshift, bc;
    
    if (seqtype == Seq_code_gap) return FALSE;

    bysp = (ByteStorePtr) sdp;
    if (bysp == NULL)
    {
        ErrPostEx(SEV_ERROR,0,0, "Error:  no sequence data\n");
        return FALSE;      
    }

    if ((sctp = SeqCodeTableFind (seqtype)) == NULL)
    {
        ErrPostEx(SEV_ERROR,0,0, "Can't open table\n");
        return FALSE;
    }
    switch (seqtype)        /*determine type of base encoding*/
    {
        case Seq_code_ncbi2na:
            bc = 4;
            rshift = 6;
            lshift = 2;
            mask = 192;
            break;

        case Seq_code_ncbi4na:
            bc = 2;
            rshift = 4;
            lshift = 4;
            mask = 240;
            break;

                case Seq_code_iupacna:
                case Seq_code_ncbi8na:
            bc = 1;
            rshift = 0;
            lshift = 0;
            mask = 255;
            break;
    case Seq_code_iupacaa:
    case Seq_code_ncbi8aa:
    case Seq_code_ncbieaa:
    case Seq_code_ncbipaa:
    case Seq_code_iupacaa3:
    case Seq_code_ncbistdaa:             /* ignore amino acid */
      ErrPostEx(SEV_ERROR,0,0, "Error:  cannot complement aa ; No ->mol flag on Bioseq\n");
      return FALSE;
    case Seq_code_ncbipna:
      ErrPostEx(SEV_WARNING,0,0, "Error: Don't yet know how to complement profile\n");
            return FALSE;
        default:
            return FALSE;
    }

    bslen = BSLen(bysp);
    bitctr = 0;
    readbyte = 0;

    while (readbyte < bslen)
    {
        if (!bitctr)
        {                /*get new byte*/
            BSSeek (bysp, readbyte, SEEK_SET);
            newbyte = byte_to = byte = residue = 0;
            byte = (Uint1)BSGetByte (bysp);
            bitctr = bc;
            readbyte++;
        }

        for (; bitctr; bitctr--)
        {
            residue = byte & mask;    /*mask out all but one base*/
            residue >>= rshift;
            byte <<= lshift;

            comp = SeqCodeTableComp (sctp, residue); /*get 
complement*/

            newbyte <<= lshift;
            byte_to = newbyte;
            newbyte = (comp | byte_to);    /*put complements 
together*/

        }

        if (readbyte)            /*put back byte with comps*/
        {
            BSSeek (bysp, readbyte-1, SEEK_SET);
            BSPutByte (bysp, newbyte);
        }
    }
    return TRUE;
  
}


static Boolean DeltaBioseqComplement (BioseqPtr bsp)
{
  DeltaSeqPtr dsp;
  SeqLitPtr   slip;
  Boolean     rval = FALSE;
  
  if (bsp == NULL || bsp->repr != Seq_repr_delta)
  {
    return rval;
  }
  
  dsp = (DeltaSeqPtr) bsp->seq_ext;
  while (dsp != NULL)
  {
    if (dsp->choice != 2) 
    {
      ErrPostEx(SEV_ERROR,0,0, "Error: Can't complement delta sequences with far locs\n");
      return FALSE;  
    }
    dsp = dsp->next;
  }
  rval = TRUE;
  dsp = (DeltaSeqPtr) bsp->seq_ext;
  while (dsp != NULL)
  {
    slip = (SeqLitPtr) dsp->data.ptrvalue;
    /* complement data */
    if (slip->seq_data != NULL)
    {
      rval &= ComplementSeqData (slip->seq_data_type, slip->length, slip->seq_data);
    }
    dsp = dsp->next;
  }
  return rval;
}


/*-------------- BioseqComplement () ---------------------------*/
/***********************************************************************
*   BioseqComplement:   Takes the nucleic acid sequence from Bioseq
*    Entry and gives the complement sequence in place
*       Does not change features.
************************************************************************/
NLM_EXTERN Boolean LIBCALL BioseqComplement (BioseqPtr bsp)
{
    Boolean         rval = FALSE;
    
  if (bsp == NULL)
  {
    ErrPostEx(SEV_ERROR,0,0, "Error: not a BioseqPtr\n");
    rval = FALSE;  
  }
  else if (ISA_aa(bsp->mol)) 
  {
    ErrPostEx(SEV_ERROR,0,0, "Error:  cannot complement aa\n");
        rval = FALSE;  
  }
  else if (bsp->repr == Seq_repr_delta)
  {
    rval = DeltaBioseqComplement (bsp);
  }
  else if (bsp->repr == Seq_repr_raw)
  {
    rval = ComplementSeqData (bsp->seq_data_type, bsp->length, bsp->seq_data);
  }
  else
  {
    ErrPostEx(SEV_ERROR,0,0, "Error: not a raw or delta sequence\n");
        rval = FALSE;  
  }
  return rval;                        

} /* BioseqComplement */


NLM_EXTERN Boolean LIBCALL ReverseSeqData (Uint1 seqtype, Int4 seqlen, SeqDataPtr sdp)
{
    ByteStorePtr     bysp1, bysp2 = '\0';
    long         readbyte, bslen = 0;
    Int4         count = 0;
    Uint1     byte = 0, byte2, byte_to = 0, byte_to2, newbyte = 0;
    Uint1        newbyte2, finalbyte, residue, residue2, bitctr, bc2 = 0;
    Uint1         bitctr2, mask, mask2, lshift, rshift, bc = 0, jagged;
    
    if (seqtype == Seq_code_gap) return FALSE;

    bysp1 = (ByteStorePtr) sdp;

  if (bysp1 == NULL)
  {
    ErrPostEx(SEV_ERROR,0,0, "Error:  No sequence data\n");
    return FALSE;
  }

    switch (seqtype){
        case Seq_code_ncbi2na:        /*bitshifts needed*/
            mask = 192;
            mask2 = 3;
            lshift = 2;
            rshift = 6;
                        jagged = seqlen%4;
            switch (jagged)    /*change if jagged last byte*/
            {
                case 1:
                    bc = 1;
                    bc2 = 3;
                    break;
                case 2:
                    bc = 2;
                    bc2 = 2;
                    break;
                case 3:
                    bc = 3;
                    bc2 = 1;
                    break;
                default:
                    bc = 4;
                    bc2 = 0;
                    break;
            }
            break;
        case Seq_code_ncbi4na:
            mask = 240;
            mask2 = 15;
            lshift = 4;
            rshift = 4;
                        jagged = seqlen%2;
            switch (jagged)
            {
                case 1:
                    bc = 1;
                    bc2 = 1;
                    break;
                default:
                    bc = 2;
                    bc2 = 0;
                    break;
            }
            break;
                case Seq_code_iupacna:
                case Seq_code_ncbi8na:

                case Seq_code_iupacaa:
                case Seq_code_ncbi8aa:
                case Seq_code_ncbieaa:
                case Seq_code_ncbistdaa:
            bc = 1;
                        bc2 = 0;
            rshift = 0;
            lshift = 0;
                        jagged = 0;
            mask = 255;
                        mask2 = 0;
            break;
                case Seq_code_ncbipaa:
                case Seq_code_iupacaa3:
                    ErrPostEx(SEV_ERROR,0,0, "Error:  cannot  reverse %s protein alphabet",(int)seqtype);
                    return FALSE;
                case Seq_code_ncbipna:
                    ErrPostEx(SEV_WARNING,0,0, "Error: Don't yet know how to reverse profile\n");
        default:        /*ignores amino acid sequence*/
            return FALSE;
    }
    bysp2 = BSDup(bysp1);
    bslen = BSLen (bysp1);
    bitctr = bitctr2 = 0;
    readbyte = 0;
    count = 0;

    if (!jagged)            /*no jagged last byte*/
    {
        while ((readbyte != BSLen(bysp1)))
        {
            count = rshift;
            if (!bitctr)        /*get new byte*/
            {
                newbyte = byte_to = byte = residue = 0;
                BSSeek (bysp2, --bslen, SEEK_SET);
                byte = (Uint1)BSGetByte (bysp2);
                bitctr = bc;
                readbyte++;
            }

            for (;bitctr; bitctr--)
            {
                residue = byte & mask;
                residue >>= count;
                byte <<= lshift;
                count = count - lshift;
    
                newbyte = (residue | byte_to);
                byte_to = newbyte;
            }

            BSSeek (bysp1, readbyte-1, SEEK_SET);
            BSPutByte (bysp1, newbyte);

        }
    }
    else                /*jagged last byte*/
    {
        /*Gets two bytes prior to loop*/
        newbyte = newbyte2 = byte_to = byte_to2 = 0;
        byte2 = residue = residue2 = 0;
        BSSeek (bysp2, bslen-2, SEEK_SET);
        byte2 = (Uint1) BSGetByte (bysp2);    /*byte closer to beginning*/
        byte = (Uint1) BSGetByte (bysp2);
        bitctr = bc;
        bitctr2 = bc2;
        bslen = bslen - 2;
        readbyte = 1;

        while (readbyte != BSLen(bysp1))
        {
            count = rshift;
            if (!bitctr)        /*when needed gets another 
byte*/
            {
                newbyte = newbyte2 = byte_to = byte_to2 = 0;
                byte2 = finalbyte = residue = residue2 = 0;
                BSSeek (bysp2, --bslen, SEEK_SET);
                byte2 = (Uint1) BSGetByte (bysp2);
                bitctr = bc;
                bitctr2 = bc2;
                ++readbyte;
            }
            for (; bitctr; bitctr--)
            {
                residue = byte & mask;        /*reverses 1st 
byte*/
                residue >>= count;
                byte <<= lshift;
                byte_to = newbyte;
                newbyte = (residue | byte_to);
                count = count - lshift;
            }
            for (; bitctr2; bitctr2--)
            {
                residue2 = byte2 & mask2;   /*reverses 2nd */
                byte2 >>= lshift;        /*partially to 
join*/
                newbyte2 <<= lshift;        /*with the 1st*/
                byte_to2 = newbyte2;
                newbyte2 = (residue2 | byte_to2);
            }
            newbyte <<= (8 - (bc*lshift));    /*joins 1st & 2nd 
bytes*/
            finalbyte = (newbyte | newbyte2);
            byte2 <<= (bc2 * lshift);
            byte = byte2;

            BSSeek (bysp1, readbyte-1, SEEK_SET);
            BSPutByte (bysp1, finalbyte);
        }
    }
    BSFree(bysp2);
    return TRUE;
} /* ReverseSeqData */


static Boolean DeltaBioseqReverse (BioseqPtr bsp)
{
  DeltaSeqPtr dsp, next_dsp, newchain = NULL;
  SeqLitPtr   slip;
  Boolean     rval = FALSE;
  
  if (bsp == NULL || bsp->repr != Seq_repr_delta)
  {
    return rval;
  }
  
  dsp = (DeltaSeqPtr) bsp->seq_ext;
  while (dsp != NULL)
  {
    if (dsp->choice != 2) 
    {
      ErrPostEx(SEV_ERROR,0,0, "Error: Can't reverse delta sequences with far locs\n");
      return FALSE;  
    }
    dsp = dsp->next;
  }
  
  dsp = (DeltaSeqPtr) bsp->seq_ext;  
  rval = TRUE;
  while (dsp != NULL)
  {
    slip = (SeqLitPtr) dsp->data.ptrvalue;
    /* reverse data */
    if (slip->seq_data != NULL)
    {
      rval &= ReverseSeqData (slip->seq_data_type, slip->length, slip->seq_data);
    }
    
    /* reverse the chain */
    next_dsp = dsp->next;
    dsp->next = newchain;
    newchain = dsp;
    
    dsp = next_dsp;
  }
  bsp->seq_ext = newchain;
  return rval;
}
           
/*-------------- BioseqReverse () ---------------------------*/
/***********************************************************************
*   BioseqReverse:   Takes nucleic acid sequence from Bioseq Entry and 
*    reverses the whole sequence in place
*       Does not change features.
************************************************************************/
NLM_EXTERN Boolean LIBCALL BioseqReverse (BioseqPtr bsp)
{
    Boolean       rval;
    
  if (bsp == NULL)
  {
    ErrPostEx(SEV_ERROR,0,0, "Error: not a BioseqPtr\n");
    rval = FALSE;  
  }
  else if (bsp->repr == Seq_repr_delta)
  {
    rval = DeltaBioseqReverse (bsp);
  }                       
  else if (bsp->repr == Seq_repr_raw)
  {
    rval = ReverseSeqData (bsp->seq_data_type, bsp->length, bsp->seq_data);
  }
  else
  {
    ErrPostEx(SEV_ERROR,0,0, "Error: not a raw or delta sequence\n");
    rval = FALSE;  
  }
                        
    return rval;
} /* BioseqReverse */

#define SPC_BUFF_CHUNK 1024

/*****************************************************************************
*
*  ContigRevComp
*
*****************************************************************************/

static Boolean SegRevComp (BioseqPtr bsp)

{
  ValNodePtr  head = NULL;
  Int4        from, to, tmp;
  Boolean     partial5, partial3;
  SeqIntPtr   sintp;
  SeqLocPtr   slp;
  ValNode     vn;
  ValNodePtr  vnp;

  MemSet ((Pointer) &vn, 0, sizeof (ValNode));
  vn.choice = SEQLOC_MIX;
  vn.data.ptrvalue = bsp->seq_ext;

  /* get each location component */

  slp = SeqLocFindNext (&vn, NULL);
  while (slp != NULL) {

    /* copy component, reversing strand */

    vnp = NULL;
    if (slp->choice == SEQLOC_NULL) {

      vnp = ValNodeAddPointer (NULL, SEQLOC_NULL, NULL);

    } else if (slp->choice == SEQLOC_INT) {

      sintp = (SeqIntPtr) slp->data.ptrvalue;
      if (sintp != NULL) {
        CheckSeqLocForPartial (slp, &partial5, &partial3);
        from = sintp->from;
        to = sintp->to;
        if (sintp->strand != Seq_strand_minus) {
          tmp = from;
          from = to;
          to = tmp;
        }
        vnp = AddIntervalToLocation (NULL, sintp->id, from, to, partial3, partial5);
      }

    }

    /* save in new list in reverse order */

    if (vnp != NULL) {
      vnp->next = head;
      head = vnp;
    }

    slp = SeqLocFindNext (&vn, slp);
  }

  if (head == NULL) return FALSE;

  bsp->seq_ext = SeqLocSetFree ((ValNodePtr) bsp->seq_ext);
  bsp->seq_ext = head;

  bsp->hist = SeqHistFree (bsp->hist);

  return TRUE;
}

static Boolean DeltaRevComp (BioseqPtr bsp)

{
  DeltaSeqPtr  dsp, dspnext;
  ValNodePtr   head = NULL;
  Int4         from, to, tmp;
  Boolean      partial5, partial3;
  SeqIntPtr    sintp;
  SeqLocPtr    slp;
  SeqLitPtr    slitp, slip;
  ValNodePtr   vnp;

  for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp != NULL; dsp = dsp->next) {
    vnp = NULL;

    if (dsp->choice == 1) {

      slp = (SeqLocPtr) dsp->data.ptrvalue;
      if (slp != NULL) {

        if (slp->choice == SEQLOC_NULL) {

          vnp = ValNodeAddPointer (NULL, SEQLOC_NULL, NULL);

        } else if (slp->choice == SEQLOC_INT) {

          sintp = (SeqIntPtr) slp->data.ptrvalue;
          if (sintp != NULL) {
            CheckSeqLocForPartial (slp, &partial5, &partial3);
            from = sintp->from;
            to = sintp->to;
            if (sintp->strand != Seq_strand_minus) {
              tmp = from;
              from = to;
              to = tmp;
            }
            vnp = AddIntervalToLocation (NULL, sintp->id, from, to, partial3, partial5);
          }
        }
      }

    } else if (dsp->choice == 2) {

      slitp = (SeqLitPtr) dsp->data.ptrvalue;
      if (slitp != NULL && slitp->seq_data == NULL) {
        slip = SeqLitNew ();
        if (slip != NULL) {
          slip->length = slitp->length;
          /* not copying fuzz */
          slip->seq_data_type = slitp->seq_data_type;
          vnp = ValNodeAddPointer (NULL, 2, (Pointer) slip);
        }
      } else {
        ValNodeFree (head);
        return FALSE;
      }
    }

    /* save in new list in reverse order */

    if (vnp != NULL) {
      vnp->next = head;
      head = vnp;
    }
  }

  if (head == NULL) return FALSE;

  dsp = (DeltaSeqPtr) bsp->seq_ext;
  while (dsp != NULL) {
    dspnext = dsp->next;
    dsp->next = NULL;
    DeltaSeqFree (dsp);
    dsp = dsp->next;
  }
  bsp->seq_ext = head;

  bsp->hist = SeqHistFree (bsp->hist);

  return TRUE;
}

NLM_EXTERN Boolean LIBCALL ContigRevComp (BioseqPtr bsp)

{
  if (bsp == NULL) {
    ErrPostEx (SEV_ERROR, 0, 0, "ContigRevComp: empty BioseqPtr");
    return FALSE;
  }

  if (bsp->repr == Seq_repr_seg && bsp->seq_ext_type == 1 && bsp->seq_ext != NULL) {
    return SegRevComp (bsp);
  }
  if (bsp->repr == Seq_repr_delta && bsp->seq_ext_type == 4 && bsp->seq_ext != NULL) {
    return DeltaRevComp (bsp);
  }

  ErrPostEx (SEV_ERROR, 0, 0, "ContigRevComp: not a segmented or delta BioseqPtr");
  return FALSE;
}

/*****************************************************************************
*
*  SPCompressNew(void); - allocated memory for SPCompress structure
*
*****************************************************************************/
NLM_EXTERN SPCompressPtr SPCompressNew(void)
{
    SPCompressPtr spc;
    
    spc = (SPCompressPtr) MemNew(sizeof(SPCompress));
    spc->buffer = (Uint1Ptr) MemNew(SPC_BUFF_CHUNK);
    spc->allocated = SPC_BUFF_CHUNK;
    spc->residues = 0;
    spc->lbytes = NULL;
    
    return spc;
}
/*****************************************************************************
*
*  SPCompressFree(SPCompressPtr spc); -  free SPCompress structure
*
*****************************************************************************/
NLM_EXTERN void SPCompressFree(SPCompressPtr spc)
{

  MemFree(spc->buffer);
  MemFree(spc->lbytes);
  MemFree(spc);

}
/*****************************************************************************
*
*  Int4 SPCompressRead (Pointer data, Uint1Ptr buf, Int4 length);
*        Hook read-function for SPCompressDNA()
*
*****************************************************************************/
static Int4 SPCompressRead (Pointer data, Uint1Ptr buf, Int4 length);
static Int4 SPCompressRead (Pointer data, Uint1Ptr buf, Int4 length)
{
  SeqPortPtr spp;
  Uint1 residue = 0;
  Int4 total_read=0, index=0;

  Boolean second = FALSE;

  spp = (SeqPortPtr) data;
  MemSet(buf, 0, length);  /* Clear buffer first */

  while (index < length && (residue=SeqPortGetResidue(spp)) != SEQPORT_EOF) {
    if (IS_residue(residue)) {
      if(second) {
        buf[index] += residue;
        index++;
        second = FALSE;
      } else {
        residue <<= 4;
        buf[index] += residue;
        second = TRUE;
      }
      total_read++;
    } else if (residue == SEQPORT_VIRT) { /* No sequence, return NULL. */
      continue;
    } else {
      ErrPostEx(SEV_WARNING, 0, 0,"[Bad residue]\n");
      return -1;
    }
  }
  return total_read;
}

/*****************************************************************************
*
*  Int4 SPCompressWrite (Pointer data, Uint1Ptr buf, Int4 length);
*        Hook write-function for SPCompressDNA()
*
*****************************************************************************/
static Int4 SPCompressWrite (Pointer data, Uint1Ptr buf, Int4 length);
static Int4 SPCompressWrite (Pointer data, Uint1Ptr buf, Int4 length)
{
  SPCompressPtr spc;
  spc = (SPCompressPtr) data;
  
  if((spc->used + length) >= spc->allocated) {
    spc->allocated += SPC_BUFF_CHUNK;
    spc->buffer = (Uint1Ptr)Realloc(spc->buffer, 
                                    spc->allocated); 
  }
  
  if((MemCpy(spc->buffer + spc->used, buf, length)) == NULL)
    return -1;
  
  spc->used += length;
  
  return length;
}

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
NLM_EXTERN Boolean SPRebuildDNA(SPCompressPtr spc)
{
    ByteStorePtr bsp, bsp_plain;
    Int4 residues;

    if(spc == NULL || spc->type != Seq_code_ncbi2na)
        return FALSE;
    
    residues = (spc->used-1)*4 + (spc->buffer[spc->used-1] & 0x3);
    bsp = BSNew(spc->used);
    BSWrite(bsp, spc->buffer, spc->used);
    
    if((bsp_plain = BSConvertSeq(bsp, Seq_code_ncbi4na, 
                                 Seq_code_ncbi2na, residues)) == NULL) {
        return FALSE;
    }
    
    BSRebuildDNA_4na(bsp_plain, spc->lbytes);
    
    spc->buffer = (Uint1Ptr) Realloc(spc->buffer, residues/2+1);
    BSRead(bsp_plain, spc->buffer, residues/2+1);
    spc->type = Seq_code_ncbi4na;
    spc->residues = residues;
    BSFree(bsp_plain);
    
    return TRUE;
}

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
NLM_EXTERN SPCompressPtr SPCompressDNA(SeqPortPtr spp)
{
  SPCompressPtr spc;
  
  if (spp == NULL || spp->newcode != Seq_code_ncbi4na)
    return NULL;
    
  spc = SPCompressNew();
  if(!GenericCompressDNA((VoidPtr) spp, (VoidPtr) spc, 
                         (Uint4) -1, /* Length of sequence unknown */
                         SPCompressRead, 
                         SPCompressWrite, 
                         &spc->lbytes
                         )) {
    return NULL;
  }
  spc->type = Seq_code_ncbi2na;
  return spc;
}

/*****************************************************************************
*
*   ComposeCodonsRecognizedString (trna, buf, buflen);
*       Copies codon recognized string to buf, returns number of codons
*
*****************************************************************************/

static int LIBCALLBACK SortCodonByName (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1;
  CharPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      str1 = (CharPtr) vnp1->data.ptrvalue;
      str2 = (CharPtr) vnp2->data.ptrvalue;
      if (str1 != NULL && str2 != NULL) {
        return StringICmp (str1, str2);
      } else {
        return 0;
      }
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

static Uint1 MakeDegenerateBase (Uint1 ch1, Uint1 ch2, Uint1Ptr chrToInt, CharPtr intToChr)

{
  Uint1  idx;

  idx = chrToInt [(int) ch1] | chrToInt [(int) ch2];
  return intToChr [(int) idx];
}

NLM_EXTERN Int2 ComposeCodonsRecognizedString (tRNAPtr trna, CharPtr buf, size_t buflen)

{
  Char          ch;
  Uint1         chrToInt [256];
  Uint1         codon [4];
  Int2          count = 0;
  ValNodePtr    head, next, vnp;
  Int2          k;
  Uint1         i, j;
  CharPtr       intToChr = "?ACMGRSVUWYHKDBN";
  CharPtr       prefix, ptr, str1, str2;
  Pointer PNTR  prev;

  if (trna == NULL || buf == NULL || buflen < 25) return 0;

  *buf = '\0';
  codon [3] = '\0';
  head = NULL;

  for (j = 0; j < 6; j++) {
    if (trna->codon [j] < 64) {
      if (CodonForIndex (trna->codon [j], Seq_code_iupacna, codon)) {
        for (k = 0; k < 3; k++) {
          if (codon [k] == 'T') {
            codon [k] = 'U';
          }
        }
        ValNodeCopyStr (&head, 0, (CharPtr) codon);
      }
    }
  }

  head = ValNodeSort (head, SortCodonByName);

  if (head == NULL) return 0;

  for (k = 0; k < 256; k++) {
    chrToInt [k] = 0;
  }
  for (i = 1; i < 16; i++) {
    ch = intToChr [i];
    chrToInt [(int) ch] = i;
  }

  count = ValNodeLen (head);
  str1 = (CharPtr) head->data.ptrvalue;
  vnp = head->next;
  prev = (Pointer PNTR) &(head->next);
  while (vnp != NULL) {
    next = vnp->next;
    str2 = (CharPtr) vnp->data.ptrvalue;
    if (str1 != NULL && str2 != NULL &&
        str1 [0] == str2 [0] && str1 [1] == str2 [1]) {
      str1 [2] = MakeDegenerateBase (str1 [2], str2 [2], chrToInt, intToChr);
      *prev = next;
      vnp->next = NULL;
      ValNodeFreeData (vnp);
    } else {
      str1 = str2;
      prev = (Pointer PNTR) &(vnp->next);
    } 
    vnp = next;
  }

  for (vnp = head, ptr = buf, i = 0, prefix = NULL; vnp != NULL;
       vnp = vnp->next, prefix = ", ", i++) {
    ptr = StringMove (ptr, prefix);
    ptr = StringMove (ptr, (CharPtr) vnp->data.ptrvalue);
  }

  ValNodeFreeData (head);
  return count;
}

/*****************************************************************************
*
*   TransTableNew (Int2 genCode);
*       Initializes TransTable finite state machine for 6-frame translation
*       and open reading frame search, allowing nucleotide ambiguity characters
*
*****************************************************************************/

static Boolean SetGenCode (Int2 genCode, CharPtr PNTR ncbieaa, CharPtr PNTR sncbieaa)

{
  GeneticCodePtr  codes;
  GeneticCodePtr  gcp;
  Int4            id;
  ValNodePtr      vnp;

  if (ncbieaa == NULL || sncbieaa == NULL) return FALSE;

  codes = GeneticCodeTableLoad ();
  if (codes == NULL) return FALSE;
  for (gcp = codes; gcp != NULL; gcp = gcp->next) {
    id = 0;
    *ncbieaa = NULL;
    *sncbieaa = NULL;
    for (vnp = (ValNodePtr) gcp->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
      switch (vnp->choice) {
        case 2 :
          id = vnp->data.intvalue;
          break;
        case 3 :
          *ncbieaa = (CharPtr) vnp->data.ptrvalue;
          break;
        case 6 :
          *sncbieaa = (CharPtr) vnp->data.ptrvalue;
          break;
        default :
          break;
      }
    }
    if (genCode == id) return TRUE;
  }

  return FALSE;
}

typedef enum {
  BASE_A = 0,  /* A    */
  BASE_C,      /* C    */
  BASE_G,      /* G    */
  BASE_T,      /* T    */
  BASE_M,      /* AC   */
  BASE_R,      /* AG   */
  BASE_W,      /* AT   */
  BASE_S,      /* CG   */
  BASE_Y,      /* CT   */
  BASE_K,      /* GT   */
  BASE_V,      /* ACG  */
  BASE_H,      /* ACT  */
  BASE_D,      /* AGT  */
  BASE_B,      /* CGT  */
  BASE_N       /* ACGT */
} BaseCode;

NLM_EXTERN TransTablePtr TransTableNew (Int2 genCode)

{
  Char     ch, tpaa, btaa, tporf, btorf;
  Char     charToBase [16] = "ACGTMRWSYKVHDBN";
  Int2     fournaToBase [16] = {
             BASE_N, BASE_A, BASE_C, BASE_M, BASE_G, BASE_R, BASE_S, BASE_V,
             BASE_T, BASE_W, BASE_Y, BASE_H, BASE_K, BASE_D, BASE_B, BASE_N};
  Int2     expansions [75] = {
             BASE_A, -1,     -1,     -1,     -1,
             BASE_C, -1,     -1,     -1,     -1,
             BASE_G, -1,     -1,     -1,     -1,
             BASE_T, -1,     -1,     -1,     -1,
             BASE_A, BASE_C, -1,     -1,     -1,
             BASE_A, BASE_G, -1,     -1,     -1,
             BASE_A, BASE_T, -1,     -1,     -1,
             BASE_C, BASE_G, -1,     -1,     -1,
             BASE_C, BASE_T, -1,     -1,     -1,
             BASE_G, BASE_T, -1,     -1,     -1,
             BASE_A, BASE_C, BASE_G, -1,     -1,
             BASE_A, BASE_C, BASE_T, -1,     -1,
             BASE_A, BASE_G, BASE_T, -1,     -1,
             BASE_C, BASE_G, BASE_T, -1,     -1,
             BASE_A, BASE_C, BASE_G, BASE_T, -1};
  Boolean  goOn;
  Int2     i, j, k, st, nx, cd;
  Int2     p, q, r, x, y, z;
  Uint1    ui;
  Int2     codonidx [4] = {2, 1, 3, 0};  /* in genetic code table, T = 0, C = 1, A = 2, G = 3, */
  Int2     complidx [4] = {0, 3, 1, 2};  /* and index = (base1 * 16) + (base2 * 4) + base3 */
  CharPtr  ncbieaa = NULL, sncbieaa = NULL;
  TransTablePtr  tbl;

  tbl = (TransTablePtr) MemNew (sizeof (TransTable));
  if (tbl == NULL) return NULL;
  MemSet ((Pointer) tbl, 0, sizeof (TransTable));

  if (genCode == 7) {
    genCode = 4;
  } else if (genCode == 8) {
    genCode = 1;
  } else if (genCode == 0) {
    genCode = 1;
  }

  if ((! SetGenCode (genCode, &ncbieaa, &sncbieaa)) || ncbieaa == NULL || sncbieaa == NULL) {
    ncbieaa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    sncbieaa = "---M---------------M---------------M----------------------------";
  }

  tbl->genCode = genCode;
  StringNCpy_0 (tbl->ncbieaa, ncbieaa, sizeof (tbl->ncbieaa));
  StringNCpy_0 (tbl->sncbieaa, sncbieaa, sizeof (tbl->sncbieaa));

  /* table to convert any ASCII character to BASE_x integer from 0 through 14 */
  for (i = 0; i < 256; i++) {
    tbl->basesToIdx [i] = BASE_N;
  }

  /* map iupacna alphabet to BaseCode */
  for (ui = BASE_A; ui <= BASE_N; ui++) {
    ch = charToBase [ui];
    tbl->basesToIdx [(int) ch] = ui;
    ch = TO_LOWER (ch);
    tbl->basesToIdx [(int) ch] = ui;
  }
  tbl->basesToIdx [(int) 'U'] = BASE_T;
  tbl->basesToIdx [(int) 'u'] = BASE_T;
  tbl->basesToIdx [(int) 'X'] = BASE_N;
  tbl->basesToIdx [(int) 'x'] = BASE_N;

  /* also map ncbi4na alphabet to BaseCode */
  for (i = 0; i < 16; i++) {
    tbl->basesToIdx [(int) i] = fournaToBase [i];
  }

  /* add tbl->basesToIdx [(int) ch] to tbl->nextBase [state] to get next state */

  /* treat state 0 as already having seen NN, avoiding single and double letter states */
  tbl->nextBase [0] = 3361;

  /* states 1 through 3375 are triple letter states (AAA, AAC, ..., NNT, NNN) */
  for (i = BASE_A, st = 1; i <= BASE_N; i++) {
    for (j = BASE_A, nx = 1; j <= BASE_N; j++) {
      for (k = BASE_A; k <= BASE_N; k++, st++, nx += 15) {
        tbl->nextBase [st] = nx;
      }
    }
  }

  /* tbl->aminoAcid [state] [strand] contains amino acid encoded by state */

  /* initialize all states to return unknown amino acid X */
  for (st = 0; st < 3376; st++) {
    tbl->aminoAcid [st] [TTBL_TOP_STRAND] = 'X';
    tbl->aminoAcid [st] [TTBL_BOT_STRAND] = 'X';
    tbl->orfStart [st] [TTBL_TOP_STRAND] = '-';
    tbl->orfStart [st] [TTBL_BOT_STRAND] = '-';
  }

  /* lookup amino acid for each codon in genetic code table */
  for (i = BASE_A, st = 1; i <= BASE_N; i++) {
    for (j = BASE_A; j <= BASE_N; j++) {
      for (k = BASE_A; k <= BASE_N; k++, st++) {
        /* st = 225 * i + 15 * j + k + 1; */

        tpaa = '\0';
        btaa = '\0';
        tporf = '\0';
        btorf = '\0';
        goOn = TRUE;

        /* expand ambiguous IJK nucleotide symbols into component bases XYZ */
        for (p = i * 5, x = expansions [p]; x != -1 && goOn; p++, x = expansions [p]) {
          for (q = j * 5, y = expansions [q]; y != -1 && goOn; q++, y = expansions [q]) {
            for (r = k * 5, z = expansions [r]; z != -1 && goOn; r++, z = expansions [r]) {

              /* lookup amino acid for codon XYZ */
              cd = 16 * codonidx [x] + 4 * codonidx [y] + codonidx [z];
              ch = ncbieaa [cd];
              if (tpaa == '\0') {
                tpaa = ch;
              } else if (tpaa != ch) {
                /* allow Asx (Asp or Asn) and Glx (Glu or Gln) and Xle (Leu or Ile) */
                if ((tpaa == 'B' || tpaa == 'D' || tpaa == 'N') && (ch == 'D' || ch == 'N')) {
                  tpaa = 'B';
                } else if ((tpaa == 'Z' || tpaa == 'E' || tpaa == 'Q') && (ch == 'E' || ch == 'Q')) {
                  tpaa = 'Z';
                } else if ((tpaa == 'J' || tpaa == 'I' || tpaa == 'L') && (ch == 'I' || ch == 'L')) {
                  tpaa = 'J';
                } else {
                  tpaa = 'X';
                }
              }
              /* and translation start flag on top strand */
              ch = sncbieaa [cd];
              if (tporf == '\0') {
                tporf = ch;
              } else if (tporf != ch) {
                tporf = 'X'; /* was '-' */
              }

              /* lookup amino acid for complement of reversed ZYX */
              cd = 16 * complidx [z] + 4 * complidx [y] + complidx [x];
              ch = ncbieaa [cd];
              if (btaa == '\0') {
                btaa = ch;
              } else if (btaa != ch) {
                /* allow Asx (Asp or Asn) and Glx (Glu or Gln) and Xle (Leu or Ile) */
                if ((btaa == 'B' || btaa == 'D' || btaa == 'N') && (ch == 'D' || ch == 'N')) {
                  btaa = 'B';
                } else if ((btaa == 'Z' || btaa == 'E' || btaa == 'Q') && (ch == 'E' || ch == 'Q')) {
                  btaa = 'Z';
                } else if ((btaa == 'J' || btaa == 'I' || btaa == 'L') && (ch == 'I' || ch == 'L')) {
                  btaa = 'J';
                } else {
                  btaa = 'X';
                }
              }
              /* and translation start flag on bottom strand */
              ch = sncbieaa [cd];
              if (btorf == '\0') {
                btorf = ch;
              } else if (btorf != ch) {
                btorf = 'X'; /* was '-' */
              }

              /* drop out of loop as soon as answer is known */
              if (tpaa == 'X' && btaa == 'X' && tporf == 'X' && btorf == 'X') { /* was '-' for orfs */
                goOn = FALSE;
              }
            }
          }
        }

        /* assign amino acid */
        tbl->aminoAcid [st] [TTBL_TOP_STRAND] = tpaa;
        tbl->aminoAcid [st] [TTBL_BOT_STRAND] = btaa;

        /* assign orf start */
        tbl->orfStart [st] [TTBL_TOP_STRAND] = tporf;
        tbl->orfStart [st] [TTBL_BOT_STRAND] = btorf;
      }
    }
  }

  /* finite state machine for 6-frame translation and ORF search is now initialized */
  return tbl;
}

NLM_EXTERN TransTablePtr TransTableFree (TransTablePtr tbl)

{
  return MemFree (tbl);
}

NLM_EXTERN void TransTableFreeAll (void)

{
    Int2           genCode;
    Char           str [32];
    TransTablePtr  tbl;
    
    for (genCode = 1; genCode < 40; genCode++) {
        sprintf (str, "TransTableFSAforGenCode%d", (int) genCode);
        tbl = (TransTablePtr) GetAppProperty (str);
        if (tbl != NULL) {
            SetAppProperty (str, NULL);
            TransTableFree (tbl);
        }
    }
    return;
}

/* convenience function does translation for entire bioseq */

NLM_EXTERN void TransTableProcessBioseq (
  TransTablePtr tbl,
  TransTableMatchProc matchProc,
  Pointer userdata,
  BioseqPtr bsp
)

{
  Boolean     altStart, atgStart, orfStop;
  Byte        bases [400];
  Char        ch;
  Int2        ctr, frame, i, j, state;
  Int4        position;
  Uint1       residue;
  SeqPortPtr  spp;

  if (tbl == NULL || matchProc == NULL || bsp == NULL) return;

  if (! ISA_na (bsp->mol)) return;

  spp = SeqPortNew (bsp, 0, -1, 0, Seq_code_iupacna);
  if (spp == NULL) return;

  if (bsp->repr == Seq_repr_delta) {
    SeqPortSet_do_virtual (spp, TRUE);
  }

  /* read first block of bases, reality check on length */

  ctr = SeqPortRead (spp, bases, sizeof (bases));
  if (ctr < 6) {
    SeqPortFree (spp);
    return;
  }

  state = 0;
  position = 0;
  frame = 0;

  i = 0;
  residue = (Uint1) bases [i];

  /* prime finite state machine with first two bases */

  for (j = 0; j < 2 && residue != SEQPORT_EOF; j++) {
    if (IS_residue (residue)) {
      state = NextCodonState (tbl, state, residue);
    }
    i++;
    residue = (Uint1) bases [i];
  }

  /* loop on all remaining bases */

  while (residue != SEQPORT_EOF) {
    if (IS_residue (residue)) {
      state = NextCodonState (tbl, state, residue);

      /* get amino acid for codon on top strand */

      ch = GetCodonResidue (tbl, state, TTBL_TOP_STRAND);
      atgStart = IsATGStart (tbl, state, TTBL_TOP_STRAND);
      altStart = IsAltStart (tbl, state, TTBL_TOP_STRAND);
      orfStop = IsOrfStop (tbl, state, TTBL_TOP_STRAND);
      matchProc (position, ch, atgStart, altStart, orfStop, frame, Seq_strand_plus, userdata);

      /* get amino acid for codon on top strand */

      ch = GetCodonResidue (tbl, state, TTBL_BOT_STRAND);
      atgStart = IsATGStart (tbl, state, TTBL_BOT_STRAND);
      altStart = IsAltStart (tbl, state, TTBL_BOT_STRAND);
      orfStop = IsOrfStop (tbl, state, TTBL_BOT_STRAND);
      matchProc (position, ch, atgStart, altStart, orfStop, frame, Seq_strand_minus, userdata);

      /* advance base position, also keep track of frame */

      position++;
      frame++;
      if (frame > 2) {
        frame = 0;
      }
    }

    /* increment base counter */

    i++;
    if (i >= ctr) {
      i = 0;

      /* read next block of bases */

      ctr = SeqPortRead (spp, bases, sizeof (bases));
      if (ctr < 0) {
        bases [0] = -ctr;
      } else if (ctr < 1) {
        bases [0] = SEQPORT_EOF;
      }
    }
    residue = (Uint1) bases [i];
  }

  SeqPortFree (spp);
}

/* trans table translation function can be passed cds feature or individual parameters */

typedef struct readcdsdata {
  CharPtr  tmp;
  size_t   frame;
  Int4     max;
  Boolean  overflow;
} ReadCdsData, PNTR ReadCdsPtr;

/* callback allows skipping one or two bases at beginning */

static void LIBCALLBACK SaveCdsBases (
  CharPtr sequence,
  Pointer userdata
)

{
  Char        ch;
  CharPtr     from, to;
  unsigned int len;
  Int4        max;
  ReadCdsPtr  rcp;

  rcp = (ReadCdsPtr) userdata;

  if (rcp->frame > 0) {
    len = 0;
    ch = sequence [len];
    while (ch != '\0' && len <= rcp->frame) {
      len++;
      ch = sequence [len];
    }
    /* len = StringLen (sequence); */
    if (rcp->frame >= len) {

      /* unusual locations can have fewer bases in the first segments than the frame, so just decrement */

      rcp->frame -= len;
      return;
    }
  }

  /* rcp->tmp = StringMove (rcp->tmp, sequence + rcp->frame); */

  from = sequence + rcp->frame;
  to = rcp->tmp;
  max = rcp->max;

  ch = *from;
  while (ch != '\0' && max > 0) {
    *to = ch;
    to++;
    from++;
    ch = *from;
    max--;
  }
  *to = '\0';
  if (ch != '\0') {
    rcp->overflow = TRUE;
  }

  rcp->tmp = to;
  rcp->max = max;

  rcp->frame = 0;
}

NLM_EXTERN CharPtr ReadCodingRegionBases (SeqLocPtr location, Int4 len, Uint1 frame, Int4Ptr totalP)

{
  CharPtr      bases, txt;
  Int4         mod;
  ReadCdsData  rcd;
  /*
  Int2         actual, cnt;
  BioseqPtr    bsp;
  Int4         mod, position;
  SeqIdPtr     sip;
  SeqLocPtr    slp;
  SeqPortPtr   spp;
  */

  bases = MemNew ((size_t) (len + 6));
  if (bases == NULL) return NULL;

  rcd.tmp = bases;
  rcd.max = len;
  rcd.overflow = FALSE;

  /* adjust start position */

  if (frame == 2) {
    rcd.frame = 1;
  } else if (frame == 3) {
    rcd.frame = 2;
  } else {
    rcd.frame = 0;
  }

  SeqPortStreamLoc (location, STREAM_EXPAND_GAPS, (Pointer) &rcd, SaveCdsBases);

  txt = rcd.tmp;

  if (rcd.overflow) {
    ErrPostEx (SEV_ERROR, 0, 0, "ReadCodingRegionBases overflow caught");
  }

#if 0
  spp = SeqPortNewByLoc (location, Seq_code_iupacna);
  if (spp == NULL) {
    MemFree (bases);
    return NULL;
  }

  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL) {
    sip = SeqLocId (slp);
    if (sip != NULL) {
      bsp = BioseqFind (sip);
      if (bsp != NULL) {
        if (bsp->repr == Seq_repr_delta || bsp->repr == Seq_repr_virtual) {
          SeqPortSet_do_virtual (spp, TRUE);
        }
      }
    }
    slp = SeqLocFindNext (location, slp);
  }

  /* adjust start position */

  if (frame == 2) {
    position = 1;
  } else if (frame == 3) {
    position = 2;
  } else {
    position = 0;
  }
  SeqPortSeek (spp, position, SEEK_SET);
  len -= position;

  /* read nucleotides into temporary buffer */

  cnt = (Int2) MIN (len, 32000L);
  txt = bases;
  actual = 1;
  while (cnt > 0 && len > 0 && actual > 0) {
    actual = SeqPortRead (spp, (BytePtr) txt, cnt);
    if (actual < 0) {
      actual = -actual;
      if (actual == SEQPORT_VIRT || actual == SEQPORT_EOS) {
        actual = 1; /* ignore, keep going */
      } else if (actual == SEQPORT_EOF) {
        actual = 0; /* stop */
      }
    } else if (actual > 0) {
      len -= actual;
      txt += actual;
      cnt = (Int2) MIN (len, 32000L);
    }
  }

  SeqPortFree (spp);
#endif

  /* pad incomplete last codon with Ns */

  len = StringLen (bases);
  if (len > 0) {
    mod = len % 3;
    if (mod == 1) {
      txt = StringMove (txt, "NN");
    } else if (mod == 2) {
      txt = StringMove (txt, "N");
    }
  }
  if (totalP != NULL) {
    *totalP = len;
  }

  return bases;
}

NLM_EXTERN ValNodePtr MakeCodeBreakList (SeqLocPtr cdslocation, Int4 len, CodeBreakPtr cbp, Uint1 frame)

{
  Int4        adjust = 0, pos, pos1, pos2;
  SeqLocPtr   tmp;
  ValNodePtr  vnphead = NULL;

  if (cdslocation == NULL || cbp == NULL) return NULL;

  if (frame == 2) {
    adjust = 1;
  } else if (frame == 3) {
    adjust = 2;
  } else {
    adjust = 0;
  }

  while (cbp != NULL) {
    pos1 = INT4_MAX;
    pos2 = -10;
    tmp = NULL;

    while ((tmp = SeqLocFindNext (cbp->loc, tmp)) != NULL) {
      pos = GetOffsetInLoc (tmp, cdslocation, SEQLOC_START);
      if (pos < pos1) {
        pos1 = pos;
      }
      pos = GetOffsetInLoc (tmp, cdslocation, SEQLOC_STOP);
      if (pos > pos2) {
        pos2 = pos;
      }
    }

    pos = pos2 - pos1; /* codon length */
    /* allow partial codon at the end */
    if (pos == 2 || (pos >= 0 && pos <= 1 && pos2 == len - 1)) {
      pos1 -= adjust;
      ValNodeAddInt (&vnphead, (Int2) cbp->aa.value.intvalue, (Int4) (pos1 / 3));
    }

    cbp = cbp->next;
  }

  vnphead = ValNodeSort (vnphead, SortByIntvalue);

  return vnphead;
}

static ByteStorePtr TransTableTranslateCommon (
  TransTablePtr  PNTR tblptr,
  SeqLocPtr location,
  SeqLocPtr product,
  Boolean partial,
  Int2 genCode,
  Uint1 frame,
  CodeBreakPtr code_break,
  Boolean include_stop,
  Boolean remove_trailingX,
  Boolean no_stop_at_end_of_complete_cds,
  BoolPtr altStartP,
  Boolean farProdFetchOK
)

{
  Char           aa;
  Int2           j, state = 0;
  Boolean        bad_base, no_start, check_start, got_stop,
                 incompleteLastCodon, use_break, is_first;
  CharPtr        bases, txt, protseq;
  ByteStorePtr   bs;
  ValNodePtr     codebreakhead = NULL, vnp;
  TransTablePtr  localtbl = NULL, tbl;
  Uint2          part_prod = 0, part_loc = 0;
  Int4           dnalen, protlen, total, k, p, q;
  Uint1          residue = 0;

  /* if table pointer not passed in from calling stack, use local table */

  if (tblptr == NULL) {
    tblptr = &localtbl;
  }

  if (location == NULL) return NULL;
  dnalen = SeqLocLen (location);
  if (dnalen < 1) return NULL;

  /* adjust for obsolete genetic code numbers */

  if (genCode == 7) {
    genCode = 4;
  } else if (genCode == 8) {
    genCode = 1;
  } else if (genCode == 0) {
    genCode = 1;
  }

  /* can store table for reuse on calling function's stack, replace if code is changing */

  tbl = *tblptr;
  if (tbl != NULL && genCode != tbl->genCode) {
    tbl = TransTableFree (tbl);
    *tblptr = tbl;
  }
  if (tbl == NULL) {
    tbl = TransTableNew (genCode);
    *tblptr = tbl;
  }
  if (tbl == NULL) return NULL;

  /* read bases, pad last codon with Ns, get total base count without padding */

  bases = ReadCodingRegionBases (location, dnalen, frame, &total);
  if (bases == NULL) {
    TransTableFree (localtbl);
    return NULL;
  }

  /* reality check on length */

  if (StringLen (bases) < 3) {
    MemFree (bases);
    TransTableFree (localtbl);
    return NULL;
  }

  /* process code breaks into list of aa (choice) and protein offset (data.intvalue) */

  codebreakhead = MakeCodeBreakList (location, dnalen, code_break, frame);

  no_start = FALSE;
  part_loc = SeqLocPartialCheck (location);
  part_prod = SeqLocPartialCheckEx (product, farProdFetchOK);
  if ((part_loc & SLP_START) /* || (part_prod & SLP_START) */) {
    no_start = TRUE;
  }
  if (StringHasNoText (tbl->sncbieaa) || no_start || frame > 1) {
    check_start = FALSE;
  } else {
    check_start = TRUE;
  }

  /* size of protein, allow partial codon at end */

  protlen = dnalen;
  protlen /= 3;
  protlen += 1;

  protseq = (CharPtr) MemNew ((size_t) protlen + 2);
  if (protseq == NULL) {
    MemFree (bases);
    ValNodeFree (codebreakhead);
    TransTableFree (localtbl);
    return NULL;
  }

  bs = BSNew (protlen);
  if (bs == NULL) {
    MemFree (bases);
    MemFree (protseq);
    ValNodeFree (codebreakhead);
    TransTableFree (localtbl);
    return NULL;
  }

  got_stop = FALSE;
  incompleteLastCodon = FALSE;
  is_first = TRUE;
  use_break = FALSE;
  state = 0;

  k = 0;
  p = 0;
  q = 0;
  txt = bases;
  residue = (Uint1) *txt;

  if (altStartP != NULL) {
    *altStartP = FALSE;
  }

  /* loop through all codons */

  while (residue != '\0') {
    for (j = 0, bad_base = FALSE; j < 3; j++, k++, txt++, residue = (Uint1) *txt) {
      if (IS_residue (residue)) {
        state = NextCodonState (tbl, state, residue);
      } else {
        state = NextCodonState (tbl, state, 'N');
        bad_base = TRUE;
      }
    }

    for (vnp = codebreakhead; vnp != NULL && vnp->data.intvalue != p; vnp = vnp->next) continue;
    use_break = (Boolean) (vnp != NULL);

    if (use_break) {
      aa = (Char) vnp->choice;
    } else if (bad_base) {
      aa = 'X';
    } else if (is_first && check_start) {

      /* ambiguous start codon that MAY be an initiator now translated to ambiguous X amino acid */
      aa = GetStartResidue (tbl, state, TTBL_TOP_STRAND);
      if (aa == '-') {
        if ((! ((part_loc & SLP_STOP) || (part_prod & SLP_STOP))) && (partial)) {
          aa = GetCodonResidue (tbl, state, TTBL_TOP_STRAND);
        }
      } else {
        if (altStartP != NULL) {
          if (IsAltStart (tbl, state, TTBL_TOP_STRAND)) {
            *altStartP = TRUE;
          }
        }
      }
    } else {

      aa = GetCodonResidue (tbl, state, TTBL_TOP_STRAND);
    }
    is_first = FALSE;

    if (aa == '*'
      && (! include_stop
        || (no_stop_at_end_of_complete_cds && ! partial && *(txt + 1) == 0)))
    {
      got_stop = TRUE;
      residue = '\0'; /* signal end of loop */

    } else {

      if (q < protlen) { /* protect against accidental buffer overflow */
        protseq [q] = aa;
      }
      q++;
      /*
      BSPutByte (bs, (Int2) aa);
      */
    }

    /* advance protein position for code break test */

    p++;
  }

  if (q > protlen) {
    ErrPostEx (SEV_ERROR, 0, 0, "TransTableTranslate - %ld characters written, %ld characters expected", (long) q, (long) protlen);
  }

  if (k > total) {
    incompleteLastCodon = TRUE;
  }

  if ((! got_stop) && incompleteLastCodon && q > 0) {
    aa = protseq [q - 1];
    if ((aa == 'X' /* || aa == 'B' || aa == 'Z' */) && q > 0) {
      q--;
    }
#if 0
    BSSeek (bs, -1, SEEK_END);  /* remove last X if incomplete last codon */
    aa = (Char) BSGetByte (bs);
    if ((aa == 'X' /* || aa == 'B' || aa == 'Z' */) && BSLen (bs) > 0) {
      BSSeek (bs, -1, SEEK_END);
      BSDelete (bs, 1);
      BSSeek (bs, -1, SEEK_END);
    }
#endif
  }

  if ((! got_stop) && remove_trailingX && q > 0) { /* only remove trailing X on partial CDS */
    aa = protseq [q - 1];
    while ((aa == 'X' /* || aa == 'B' || aa == 'Z' */) && q > 0) {
      q--;
      aa = protseq [q - 1];
    }
#if 0
    BSSeek (bs, -1, SEEK_END);  /* back up to last residue */
    aa = (Char) BSGetByte (bs);
    while ((aa == 'X' /* || aa == 'B' || aa == 'Z' */) && BSLen (bs) > 0) {
      BSSeek (bs, -1, SEEK_END);
      BSDelete (bs, 1);
      BSSeek (bs, -1, SEEK_END);
      aa = (Char) BSGetByte (bs);
    }
#endif
  }

  BSWrite (bs, (Pointer) protseq, q);

  if (BSLen (bs) < 1) {
    bs = BSFree (bs);
  }

  /* clean up temporarily allocated memory */

  MemFree (bases);
  MemFree (protseq);
  ValNodeFree (codebreakhead);

  /* free local table, if allocated */

  TransTableFree (localtbl);

  return bs;
}

/* public functions for trans table translation */

NLM_EXTERN ByteStorePtr TransTableTranslateSeqLoc (
  TransTablePtr  PNTR tblptr,
  SeqLocPtr location,
  Int2 genCode,
  Uint1 frame,
  Boolean include_stop,
  Boolean remove_trailingX
)

{
  return TransTableTranslateCommon (tblptr, location, NULL, FALSE, genCode,
                                    frame, NULL, include_stop,
                                    remove_trailingX, FALSE, NULL, TRUE);
}

NLM_EXTERN ByteStorePtr TransTableTranslateCdRegionEx (
  TransTablePtr  PNTR tblptr,
  SeqFeatPtr cds,
  Boolean include_stop,
  Boolean remove_trailingX,
  Boolean no_stop_at_end_of_complete_cds,
  BoolPtr altStartP,
  Boolean farProdFetchOK
)

{
  CdRegionPtr  crp;
  Int2         genCode = 0;
  ValNodePtr   vnp;
  Boolean      partial5, partial3;

  if (cds == NULL || cds->data.choice != SEQFEAT_CDREGION) return NULL;
  crp = (CdRegionPtr) cds->data.value.ptrvalue;
  if (crp == NULL) return NULL;

  /* set genCode variable from genetic_code parameter, if id choice is used */

  if (crp->genetic_code != NULL) {
    vnp = (ValNodePtr) crp->genetic_code->data.ptrvalue;
    while (vnp != NULL) {
      if (vnp->choice == 2) {
        genCode = (Int2) vnp->data.intvalue;
      }
      vnp = vnp->next;
    }
  }
  CheckSeqLocForPartial (cds->location, &partial5, &partial3);

  return TransTableTranslateCommon (tblptr, cds->location, cds->product, partial3,
                                    genCode, crp->frame, crp->code_break,
                                    include_stop, remove_trailingX,
                                    no_stop_at_end_of_complete_cds, altStartP, farProdFetchOK);
}

NLM_EXTERN ByteStorePtr TransTableTranslateCdRegion (
  TransTablePtr  PNTR tblptr,
  SeqFeatPtr cds,
  Boolean include_stop,
  Boolean remove_trailingX,
  Boolean no_stop_at_end_of_complete_cds
)

{
  return TransTableTranslateCdRegionEx (tblptr, cds, include_stop, remove_trailingX,
                                        no_stop_at_end_of_complete_cds, NULL, TRUE);
}

/* allow reuse of translation tables by saving as AppProperty */

static TransTablePtr  PersistentTransTableCommon (
  SeqFeatPtr cds,
  Int2 genCode
)

{
  CdRegionPtr    crp;
  Char           str [32];
  TransTablePtr  tbl = NULL;
  ValNodePtr     vnp;

  if (cds != NULL && cds->data.choice == SEQFEAT_CDREGION) {
    crp = (CdRegionPtr) cds->data.value.ptrvalue;
    if (crp != NULL && crp->genetic_code != NULL) {
      vnp = (ValNodePtr) crp->genetic_code->data.ptrvalue;
      while (vnp != NULL) {
        if (vnp->choice == 2) {
          genCode = (Int2) vnp->data.intvalue;
        }
        vnp = vnp->next;
      }
    }
  }

  if (genCode == 7) {
    genCode = 4;
  } else if (genCode == 8) {
    genCode = 1;
  } else if (genCode == 0) {
    genCode = 1;
  }

  /* set app property name for storing desired FSA */

  sprintf (str, "TransTableFSAforGenCode%d", (int) genCode);

  /* get FSA for desired genetic code if it already exists */

  tbl = (TransTablePtr) GetAppProperty (str);

  /* if not already exists, save FSA in genetic code-specific app property name */

  if (tbl == NULL) {
    tbl = TransTableNew (genCode);
    SetAppProperty (str, (Pointer) tbl);
  }

  return tbl;
}

NLM_EXTERN TransTablePtr PersistentTransTableByGenCode (
  Int2 genCode
)

{
  return PersistentTransTableCommon (NULL, genCode);
}

NLM_EXTERN TransTablePtr PersistentTransTableByCdRegion (
  SeqFeatPtr cds
)

{
  return PersistentTransTableCommon (cds, 0);
}

/*****************************************************************************
*
*   SeqSearch
*       Initializes SeqSearch finite state machine for sequence searching
*       Based on Practical Algorithms for Programmers by Binstock and Rex
*
*****************************************************************************/

/* general purpose DNA sequence search finite state machine */

typedef struct seqpattern {
  CharPtr           name;
  CharPtr           pattern;
  Int2              cutSite;
  Uint1             strand;
  struct seqpattern * next;
} SeqPatternItem, PNTR SeqPatternPtr;

typedef struct seqmatch {
  CharPtr         name;
  CharPtr         pattern;
  Int2            cutSite;
  Uint1           strand;
  struct seqmatch * next;
} SeqMatchItem, PNTR SeqMatchPtr;

typedef struct seqstate {
  Int2         onfailure;
  Int2         transitions [15]; /* order is ACGTMRWSYKVHDBN */
  SeqMatchPtr  matches;
} SeqStateItem, PNTR SeqStatePtr;

typedef struct SeqSearch {
  SeqStatePtr         stateArray;
  SeqPatternPtr       patternList;
  Int4                maxPatLen;
  Int2                maxState;
  Int2                highState;
  Int2                currentState;
  Int4                currentPos;
  Boolean             primed;
  SeqSearchMatchProc  matchproc;
  Pointer             userdata;
  Uint1               letterToIdx [256];
  Uint1               letterToComp [256];
} SeqSearchData;

#define FAIL_STATE -1

/* returns next state given current state and next character */

static Int2 SeqSearchGotoState (
  SeqSearchPtr tbl,
  Int2 state,
  Char ch,
  Boolean zeroFailureReturnsZero
)

{
  int          index;
  Int2         newstate;
  SeqStatePtr  sp;

  sp = &(tbl->stateArray [(int) state]);
  index = tbl->letterToIdx [(int) (Uint1) ch];
  newstate = sp->transitions [index];

  if (newstate != 0) return newstate;

  if (state == 0 && zeroFailureReturnsZero) return 0;

  return FAIL_STATE;
}

/* returns state to check next if current pattern broken */

static Int2 SeqSearchFailState (
  SeqSearchPtr tbl,
  Int2 state
)

{
  SeqStatePtr  sp;

  sp = &(tbl->stateArray [(int) state]);
  return sp->onfailure;
}

/* add a single character transition from one state to another */

static void SeqSearchAddTransition (
  SeqSearchPtr tbl,
  Int2 oldState,
  Char ch,
  Int2 newState
)

{
  int          index;
  SeqStatePtr  sp;

  sp = &(tbl->stateArray [(int) oldState]);
  index = tbl->letterToIdx [(int) (Uint1) ch];
  sp->transitions [index] = newState;
}

/* given state should report a successful match */

static void SeqSearchAddOutput (
  SeqSearchPtr tbl,
  Int2 state,
  CharPtr name,
  CharPtr pattern,
  Int2 cutSite,
  Uint1 strand
)

{
  SeqMatchPtr  mp;
  SeqStatePtr  sp;

  sp = &(tbl->stateArray [(int) state]);
  for (mp = sp->matches; mp != NULL; mp = mp->next) {
    if (StringCmp (name, mp->name) == 0) return;
  }

  mp = (SeqMatchPtr) MemNew (sizeof (SeqMatchItem));
  if (mp == NULL) return;

  mp->name = StringSave (name);
  mp->pattern = StringSave (pattern);
  mp->cutSite = cutSite;
  mp->strand = strand;

  mp->next = sp->matches;
  sp->matches = mp;
}

/* add one nucleotide sequence pattern to the finite state machine */

static Int2 SeqSearchEnterNucWord (
  SeqSearchPtr tbl,
  Int2 highState,
  Int2 maxState,
  CharPtr name,
  CharPtr pattern,
  Int2 cutSite,
  Uint1 strand
)

{
  Char     ch;
  Int2     next, patLen, state;
  CharPtr  ptr;

  state = 0;
  next = 0;

  patLen = StringLen (pattern);

  /* try to overlay beginning of pattern onto existing table */

  for (ptr = pattern, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    next = SeqSearchGotoState (tbl, state, ch, FALSE);
    if (next == FAIL_STATE) break;
    state = next;
  }

  /* now create new states for remaining characters in pattern */

  for ( ; ch != '\0'; ptr++, ch = *ptr) {
    highState++;
    SeqSearchAddTransition (tbl, state, ch, highState);
    state = highState;
  }

  /* at end of pattern record match information */

  SeqSearchAddOutput (tbl, state, name, pattern, cutSite, strand);

  return highState;
}

/* FIFO queue and other functions for building failure states */

static void SeqSearchQueueAdd (
  Int2Ptr queue,
  Int2 qbeg,
  Int2 val
)

{
  Int2  q;

  q = queue [qbeg];
  if (q == 0) {
    queue [qbeg] = val;
  } else {
    for ( ; queue [q] != 0; q = queue [q]) continue;
    queue [q] = val;
  }
  queue [val] = 0;
}

static void SeqSearchFindFail (
  SeqSearchPtr tbl,
  Int2 state,
  Int2 newState,
  Char ch
)

{
  SeqMatchPtr  mp;
  Int2         next;
  SeqStatePtr  sp;

  /* traverse existing failure path */

  while ((next = SeqSearchGotoState (tbl, state, ch, TRUE)) == FAIL_STATE) {
    state = SeqSearchFailState (tbl, state);
  }

  /* add new failure state */

  sp = &(tbl->stateArray [(int) newState]);
  sp->onfailure = next;

  /* add matches of substring at new state */

  sp = &(tbl->stateArray [(int) next]);
  for (mp = sp->matches; mp != NULL; mp = mp->next) {
    SeqSearchAddOutput (tbl, newState, mp->name, mp->pattern,
                        mp->cutSite, mp->strand);
  }
}

static void SeqSearchComputeFail (
  SeqSearchPtr tbl,
  Int2Ptr queue
)

{
  CharPtr      charToNuc = "ACGTMRWSYKVHDBN";
  Char         ch;
  Int2         qbeg, r, s, state;
  int          index;
  SeqStatePtr  sp;

  qbeg = 0;
  queue [0] = 0;

  /* queue up states reached directly from state 0 (depth 1) */

  sp = &(tbl->stateArray [0]);
  for (index = 0; index < 15; index++) {
    s = sp->transitions [index];
    if (s == 0) continue;
    sp->onfailure = 0;
    SeqSearchQueueAdd (queue, qbeg, s);
  }

  while (queue [qbeg] != 0) {
    r = queue [qbeg];
    qbeg = r;

    /* depth 1 states beget depth 2 states, etc. */

    sp = &(tbl->stateArray [r]);
    for (index = 0; index < 15; index++) {
      ch = charToNuc [index];
      s = sp->transitions [index];
      if (s == 0) continue;
      SeqSearchQueueAdd (queue, qbeg, s);

      /*
         Search for nucleotide sequences GTCGAC and TCATGA

         State   Substring   Transitions   Failure
           2       GT          C ->   3       7
           3       GTC         G ->   4       ?
           ...
           7       T           C ->   8       0
           8       TC          A ->   9

         For example, r = 2 (GT), if 'C' would go to s = 3 (GTC).
         From previous computation, 2 (GT) fails to 7 (T).  So we
         are not in a pattern starting with GT, but we may be in
         a pattern starting with the next character after G, or T.
         Thus, check state 7 (T) for any transitions using 'C'.
         Since 7 (T) 'C' -> 8 (TC), therefore set fail [3] -> 8.
      */

      state = SeqSearchFailState (tbl, r);
      SeqSearchFindFail (tbl, state, s, ch);
    }
  }
}

/* on first character, populate state transition table */

static void SeqSearchPrimeStateArray (
  SeqSearchPtr tbl
)

{
  Int2           highState, maxState;
  SeqPatternPtr  pp;
  Int2Ptr        queue;
  SeqStatePtr    stateArray;

  if (tbl == NULL || tbl->primed || tbl->patternList == NULL) return;

  for (maxState = 1, pp = tbl->patternList; pp != NULL; pp = pp->next) {
    maxState += StringLen (pp->pattern);
  }

  if (maxState > 4000) {
    Message (MSG_POST, "FiniteStateSearch cannot handle %d states", (int) maxState);
    return;
  }

  stateArray = (SeqStatePtr) MemNew (sizeof (SeqStateItem) * (size_t) maxState);
  queue = (Int2Ptr) MemNew (sizeof (Int2) * maxState);

  if (stateArray == NULL || queue == NULL) {
    MemFree (stateArray);
    MemFree (queue);
    Message (MSG_POST, "SequenceSearch unable to allocate buffers");
    return;
  }

  tbl->stateArray = stateArray;
  tbl->maxState = maxState;

  for (highState = 0, pp = tbl->patternList; pp != NULL; pp = pp->next) {
    highState = SeqSearchEnterNucWord (tbl, highState, maxState, pp->name,
                                       pp->pattern, pp->cutSite, pp->strand);
  }

  SeqSearchComputeFail (tbl, queue);

  MemFree (queue);

  tbl->highState = highState;
  tbl->currentState = 0;
  tbl->currentPos = 0;
  tbl->primed = TRUE;
}

/* for testing, print summary of transition table */

/*
static void PrintSeqSearchTable (
  SeqSearchPtr tbl,
  FILE *fp
)

{
  Int2         i;
  SeqMatchPtr  mp;
  SeqStatePtr  sp;
  Int2         state;

  if (tbl == NULL || fp == NULL) return;
  if (! tbl->primed) {
    SeqSearchPrimeStateArray (tbl);
  }
  if (tbl->stateArray == NULL) return;
  if (tbl->highState > 99) return;

  fprintf (fp, "State Fail A  C  G  T  M  R  W  S  Y  K  V  H  D  B  N\n");

  for (state = 0; state <= tbl->highState; state++) {
    sp = &(tbl->stateArray [(int) state]);
    fprintf (fp, " %3d  %3d", (int) state, (int) sp->onfailure);

    for (i = 0; i < 15; i++) {
      if (sp->transitions [i] != 0) {
        fprintf (fp, "%3d", (int) sp->transitions [i]);
      } else {
        fprintf (fp, "   ");
      }
    }

    for (mp = sp->matches; mp != NULL; mp = mp->next) {
      fprintf (fp, " %s", mp->name);
    }

    fprintf (fp, "\n");
  }
}
*/

/* create empty nucleotide sequence search finite state machine */

NLM_EXTERN SeqSearchPtr SeqSearchNew (
  SeqSearchMatchProc matchproc,
  Pointer userdata
)

{
  CharPtr       charToNuc = "ACGTMRWSYKVHDBN";
  Char          ch, lttr;
  CharPtr       complementBase = " TVGH  CD  M KN   YSAABW R ";
  Int2          i;
  Uint1         k;
  SeqSearchPtr  tbl;

  if (matchproc == NULL) return NULL;
  tbl = (SeqSearchPtr) MemNew (sizeof (SeqSearchData));
  if (tbl == NULL) return NULL;

  tbl->stateArray = NULL;
  tbl->patternList = NULL;
  tbl->maxPatLen = 0;
  tbl->maxState = 0;
  tbl->highState = 0;
  tbl->currentState = 0;
  tbl->currentPos = 0;
  tbl->matchproc = matchproc;
  tbl->userdata = userdata;
  tbl->primed = FALSE;

  /* initialize table to convert character to transition index from 0 (A) to 14 (N) */

  for (i = 0; i < 256; i++) {
    tbl->letterToIdx [i] = 14;
  }
  for (k = 0; k < 15; k++) {
    ch = charToNuc [k];
    tbl->letterToIdx [(int) ch] = k;
    ch = TO_LOWER (ch);
    tbl->letterToIdx [(int) ch] = k;
  }
  tbl->letterToIdx [(int) 'U'] = tbl->letterToIdx [(int) 'T'];
  tbl->letterToIdx [(int) 'u'] = tbl->letterToIdx [(int) 'T'];
  tbl->letterToIdx [(int) 'X'] = tbl->letterToIdx [(int) 'N'];
  tbl->letterToIdx [(int) 'x'] = tbl->letterToIdx [(int) 'N'];

  /* initialize table to convert character to complement character */

  for (i = 0; i < 256; i++) {
    tbl->letterToComp [i] = '\0';
  }
  for (ch = 'A', i = 1; ch <= 'Z'; ch++, i++) {
    lttr = complementBase [i];
    if (lttr != ' ') {
      tbl->letterToComp [(int) (Uint1) ch] = lttr;
    }
  }
  for (ch = 'a', i = 1; ch <= 'z'; ch++, i++) {
    lttr = complementBase [i];
    if (lttr != ' ') {
      tbl->letterToComp [(int) (Uint1) ch] = lttr;
    }
  }

  return tbl;
}

/* table to expand ambiguity letter to all matching nucleotide letters */

static CharPtr  nucExpandList [26] = {
  "A",
  "CGT",
  "C",
  "AGT",
  "",
  "",
  "G",
  "ACT",
  "",
  "",
  "GT",
  "",
  "AC",
  "ACGT",
  "",
  "",
  "",
  "AG",
  "CG",
  "T",
  "T",
  "ACG",
  "AT",
  "",
  "CT",
  ""
};

/* recursive function to expand and store appropriate individual patterns */

static void StoreSeqPattern (
  SeqSearchPtr tbl,
  CharPtr name,
  CharPtr str,
  Int2 cutSite,
  Uint1 strand
)

{
  Int4           patLen;
  SeqPatternPtr  pp;

  pp = (SeqPatternPtr) MemNew (sizeof (SeqPatternItem));
  if (pp == NULL) return;

  pp->name = StringSave (name);
  pp->pattern = StringSave (str);
  pp->cutSite = cutSite;
  pp->strand = strand;

  pp->next = tbl->patternList;
  tbl->patternList = pp;
  patLen = StringLen (str);
  if (patLen > tbl->maxPatLen) {
    tbl->maxPatLen = patLen;
  }
}

static void ExpandSeqPattern (
  SeqSearchPtr tbl,
  CharPtr name,
  CharPtr pattern,
  Int2 cutSite,
  Uint1 strand,
  size_t patLen,
  CharPtr str,
  Uint2 position,
  SearchFlgType flags
)

{
  Char     ch, lttr;
  Uint2     idx;
  CharPtr  ptr;

  if (position < patLen) {

    if ((Boolean) ((flags & SEQ_SEARCH_EXPAND_PATTERN) != 0)) {

      /* given ambiguity letter, get index into nucExpandList */

      ch = pattern [position];
      idx = ch - 'A';
      ptr = nucExpandList [idx];

      /* put every ACGT letter at current position, recurse for next position */

      for (lttr = *ptr; lttr != '\0'; ptr++, lttr = *ptr) {
        str [position] = lttr;
        ExpandSeqPattern (tbl, name, pattern, cutSite, strand,
                          patLen, str, position + 1, flags);
      }

    } else {

      /* if matching ambiguity characters in sequence, do not expand each base */

      str [position] = pattern [position];
      ExpandSeqPattern (tbl, name, pattern, cutSite, strand,
                        patLen, str, position + 1, flags);
    }

    /* do not run into pattern storage section of code located below */

    return;
  }

  /* when position reaches pattern length, store one fully expanded string */

  StoreSeqPattern (tbl, name, str, cutSite, strand);

  if ((Boolean) ((flags & SEQ_SEARCH_ALLOW_MISMATCH) == 0)) return;

  for (idx = 0; idx < patLen; idx++) {
    ch = str [idx];

    /* put N at every position if a single mismatch is allowed */

    str [idx] = 'N';

    StoreSeqPattern (tbl, name, str, cutSite, strand);

    /* now restore proper character, go on to put N in next position */

    str [idx] = ch;
  }
}

/* add restriction site to sequence search finite state machine */

NLM_EXTERN void SeqSearchAddNucleotidePattern (
  SeqSearchPtr tbl,
  CharPtr name,
  CharPtr pattern,
  Int2 cutSite,
  SearchFlgType flags
)

{
  Char     ch, comp [128], pat [128], str [128];
  Int2     i, j;
  size_t   len;
  Uint1    strand;
  Boolean  symmetric = TRUE;

  if (tbl == NULL || StringHasNoText (name) || StringHasNoText (pattern)) return;

  StringNCpy_0 (pat, pattern, sizeof (pat));
  TrimSpacesAroundString (pat);

  len = StringLen (pat);

  /* upper case working copy of pattern string */

  for (i = 0; i < len; i++) {
    ch = pat [i];
    pat [i] = TO_UPPER (ch);
  }

  /* reverse complement pattern to see if it is symetrical */

  for (i = 0, j = len - 1; i < len; i++, j--) {
    ch = pat [i];
    comp [j] = tbl->letterToComp [(int) (Uint1) ch];
  }
  comp [len] = '\0';
  symmetric = (Boolean) (StringICmp (pat, comp) == 0);

  if (symmetric) {
    strand = Seq_strand_both;
  } else {
    strand = Seq_strand_plus;
  }

  /* record expansion of entered pattern */

  MemSet ((Pointer) str, 0, sizeof (str));
  ExpandSeqPattern (tbl, name, pat, cutSite, strand,
                    len, str, 0, flags);

  if (symmetric) return;
  if ((Boolean) ((flags & SEQ_SEARCH_JUST_TOP_STRAND) != 0)) return;

  /* record expansion of reverse complement of asymmetric pattern */

  MemSet ((Pointer) str, 0, sizeof (str));
  ExpandSeqPattern (tbl, name, comp, len - cutSite, Seq_strand_minus,
                    len, str, 0, flags);
}

/* program passes each character in turn to finite state machine */

static void SeqSearchProcessCharacterEx (
  SeqSearchPtr tbl,
  Char ch,
  Int4 length
)

{
  Int2         curr, next;
  SeqMatchPtr  mp;
  Int4         patLen;
  SeqStatePtr  sp;

  if (tbl == NULL) return;
  if (! tbl->primed) {
    SeqSearchPrimeStateArray (tbl);
  }
  if (tbl->stateArray == NULL) return;

  curr = tbl->currentState;

  /* loop through failure states until match or back to state 0 */

  while ((next = SeqSearchGotoState (tbl, curr, ch, TRUE)) == FAIL_STATE) {
    curr = SeqSearchFailState (tbl, curr);
  }

  tbl->currentState = next;
  (tbl->currentPos)++;

  /*
     States while traversing search sequence containing EcoRI site (GAATTC)
                                                        ------
     AAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTCGAGCTCGGTACCCGGGGATCCTC
     00100010001000100110012000001211200000111000012100012345612100011000001111200000
                                                        *
  */

  /* report any matches at current state to callback function */

  sp = &(tbl->stateArray [(int) next]);
  for (mp = sp->matches; mp != NULL; mp = mp->next) {

    /* for circular sequences, prevent multiple reports of patterns */

    patLen = StringLen (mp->pattern);
    if (tbl->currentPos - patLen < length) {
      tbl->matchproc (tbl->currentPos - patLen,
                      mp->name, mp->pattern, mp->cutSite,
                      mp->strand, tbl->userdata);
    }
  }
}

NLM_EXTERN void SeqSearchProcessCharacter (
  SeqSearchPtr tbl,
  Char ch
)

{
  SeqSearchProcessCharacterEx (tbl, ch, INT4_MAX);
}

/* convenience function calls SeqSearchProcessCharacter for entire nucleotide bioseq */

typedef struct seqsrchdata {
  SeqSearchPtr  tbl;
  Int4          length;
} SeqSrchData, PNTR SeqSrchPtr;

static void LIBCALLBACK SearchSeqProc (
  CharPtr sequence,
  Pointer userdata
)

{
  Char        ch;
  CharPtr     ptr;
  SeqSrchPtr  ssp;

  ssp = (SeqSrchPtr) userdata;

  ptr = sequence;
  ch = *ptr;
  while (ch != '\0') {
    ch = TO_UPPER (ch);
    SeqSearchProcessCharacterEx (ssp->tbl, ch, ssp->length);
    ptr++;
    ch = *ptr;
  }
}

NLM_EXTERN void SeqSearchProcessBioseq (
  SeqSearchPtr tbl,
  BioseqPtr bsp
)

{
  SeqSrchData  ssd;

  SeqSearchReset (tbl);

  if (tbl == NULL || bsp == NULL) return;

  if (! ISA_na (bsp->mol)) return;

  ssd.tbl = tbl;
  ssd.length = bsp->length;

  SeqPortStream (bsp, STREAM_EXPAND_GAPS, (Pointer) &ssd, SearchSeqProc);

  /* for circular molecules, check for patterns spanning origin */

  if (bsp->topology == TOPOLOGY_CIRCULAR && bsp->length > tbl->maxPatLen) {
    SeqPortStreamInt (bsp, 0, tbl->maxPatLen, Seq_strand_plus, STREAM_EXPAND_GAPS, (Pointer) &ssd, SearchSeqProc);
  }

  SeqSearchReset (tbl);
}

/* reset state and position to allow another run with same search patterns */

NLM_EXTERN void SeqSearchReset (
  SeqSearchPtr tbl
)

{
  if (tbl == NULL) return;

  tbl->currentState = 0;
  tbl->currentPos = 0;
}

/* clean up sequence search finite state machine allocated memory */

static SeqPatternPtr FreePatternList (
  SeqPatternPtr pp
)

{
  SeqPatternPtr  next;

  while (pp != NULL) {
    next = pp->next;
    pp->next = NULL;
    MemFree (pp->name);
    MemFree (pp->pattern);
    MemFree (pp);
    pp = next;
  }

  return NULL;
}

static SeqMatchPtr FreeMatchList (
  SeqMatchPtr mp
)

{
  SeqMatchPtr  next;

  while (mp != NULL) {
    next = mp->next;
    mp->next = NULL;
    MemFree (mp->name);
    MemFree (mp->pattern);
    MemFree (mp);
    mp = next;
  }

  return NULL;
}

NLM_EXTERN SeqSearchPtr SeqSearchFree (
  SeqSearchPtr tbl
)

{
  Int2  maxState, state;

  if (tbl == NULL) return NULL;

  maxState = tbl->maxState;

  for (state = 0; state < maxState; state++) {
    FreeMatchList (tbl->stateArray [state].matches);
  }

  FreePatternList (tbl->patternList);

  MemFree (tbl->stateArray);
  return MemFree (tbl);
}

/*

static CharPtr testseq =
 "AAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCGRATYCCCGGGTACCGAGCTATYCCGAATTCGAGCTCGGTACCCGGGGATCCTCGANTTCATTCGAPTTCCAGTC";

static void MatchProc (Int4 position, CharPtr name, CharPtr pattern,
                       Int2 cutSite, Uint1 strand, Pointer userdata)

{
  Message (MSG_POST, "Name '%s', Pattern '%s', Position %ld",
           name, pattern, (long) position);
}


extern void TestSeqSearch (void);
extern void TestSeqSearch (void)

{
  Char          ch;
  CharPtr       ptr;
  SeqSearchPtr  tbl;

  tbl = SeqSearchNew (MatchProc, NULL);
  if (tbl == NULL) return;

  SeqSearchAddNucleotidePattern (tbl, "AmbiG", "GRATYC", 1, SEQ_SEARCH_EXPAND_PATTERN);
  SeqSearchAddNucleotidePattern (tbl, "ExacT", "GRAT", 1, SEQ_SEARCH_JUST_TOP_STRAND);

  for (ptr = testseq, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    SeqSearchProcessCharacter (tbl, ch);
  }

  SeqSearchFree (tbl);
}

*/

/*****************************************************************************
*
*   ProtSearch
*       Initializes ProtSearch finite state machine for sequence searching
*       Based on Practical Algorithms for Programmers by Binstock and Rex
*
*****************************************************************************/

/* general purpose protein sequence search finite state machine */

typedef struct protpattern {
  CharPtr           name;
  CharPtr           pattern;
  struct protpattern * next;
} ProtPatternItem, PNTR ProtPatternPtr;

typedef struct protmatch {
  CharPtr          name;
  CharPtr          pattern;
  struct protmatch * next;
} ProtMatchItem, PNTR ProtMatchPtr;

typedef struct protstate {
  Int2          onfailure;
  Int2          transitions [27]; /* order is ABCDEFGHIJKLMNOPQRSTUVWXYZ */
  ProtMatchPtr  matches;
} ProtStateItem, PNTR ProtStatePtr;

typedef struct ProtSearch {
  ProtStatePtr         stateArray;
  ProtPatternPtr       patternList;
  Int4                 maxPatLen;
  Int2                 maxState;
  Int2                 highState;
  Int2                 currentState;
  Int4                 currentPos;
  Boolean              primed;
  ProtSearchMatchProc  matchproc;
  Pointer              userdata;
  Uint1                letterToIdx [256];
} ProtSearchData;

#define FAIL_STATE -1

/* returns next state given current state and next character */

static Int2 ProtSearchGotoState (
  ProtSearchPtr tbl,
  Int2 state,
  Char ch,
  Boolean zeroFailureReturnsZero
)

{
  int           index;
  Int2          newstate;
  ProtStatePtr  sp;

  sp = &(tbl->stateArray [(int) state]);
  index = tbl->letterToIdx [(int) (Uint1) ch];
  newstate = sp->transitions [index];

  if (newstate != 0) return newstate;

  if (state == 0 && zeroFailureReturnsZero) return 0;

  return FAIL_STATE;
}

/* returns state to check next if current pattern broken */

static Int2 ProtSearchFailState (
  ProtSearchPtr tbl,
  Int2 state
)

{
  ProtStatePtr  sp;

  sp = &(tbl->stateArray [(int) state]);
  return sp->onfailure;
}

/* add a single character transition from one state to another */

static void ProtSearchAddTransition (
  ProtSearchPtr tbl,
  Int2 oldState,
  Char ch,
  Int2 newState
)

{
  int           index;
  ProtStatePtr  sp;

  sp = &(tbl->stateArray [(int) oldState]);
  index = tbl->letterToIdx [(int) (Uint1) ch];
  sp->transitions [index] = newState;
}

/* given state should report a successful match */

static void ProtSearchAddOutput (
  ProtSearchPtr tbl,
  Int2 state,
  CharPtr name,
  CharPtr pattern
)

{
  ProtMatchPtr  mp;
  ProtStatePtr  sp;

  sp = &(tbl->stateArray [(int) state]);
  for (mp = sp->matches; mp != NULL; mp = mp->next) {
    if (StringCmp (name, mp->name) == 0) return;
  }

  mp = (ProtMatchPtr) MemNew (sizeof (ProtMatchItem));
  if (mp == NULL) return;

  mp->name = StringSave (name);
  mp->pattern = StringSave (pattern);

  mp->next = sp->matches;
  sp->matches = mp;
}

/* add one protein sequence pattern to the finite state machine */

static Int2 ProtSearchEnterProtWord (
  ProtSearchPtr tbl,
  Int2 highState,
  Int2 maxState,
  CharPtr name,
  CharPtr pattern
)

{
  Char     ch;
  Int2     next, patLen, state;
  CharPtr  ptr;

  state = 0;
  next = 0;

  patLen = StringLen (pattern);

  /* try to overlay beginning of pattern onto existing table */

  for (ptr = pattern, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    next = ProtSearchGotoState (tbl, state, ch, FALSE);
    if (next == FAIL_STATE) break;
    state = next;
  }

  /* now create new states for remaining characters in pattern */

  for ( ; ch != '\0'; ptr++, ch = *ptr) {
    highState++;
    ProtSearchAddTransition (tbl, state, ch, highState);
    state = highState;
  }

  /* at end of pattern record match information */

  ProtSearchAddOutput (tbl, state, name, pattern);

  return highState;
}

/* FIFO queue and other functions for building failure states */

static void ProtSearchQueueAdd (
  Int2Ptr queue,
  Int2 qbeg,
  Int2 val
)

{
  Int2  q;

  q = queue [qbeg];
  if (q == 0) {
    queue [qbeg] = val;
  } else {
    for ( ; queue [q] != 0; q = queue [q]) continue;
    queue [q] = val;
  }
  queue [val] = 0;
}

static void ProtSearchFindFail (
  ProtSearchPtr tbl,
  Int2 state,
  Int2 newState,
  Char ch
)

{
  ProtMatchPtr  mp;
  Int2          next;
  ProtStatePtr  sp;

  /* traverse existing failure path */

  while ((next = ProtSearchGotoState (tbl, state, ch, TRUE)) == FAIL_STATE) {
    state = ProtSearchFailState (tbl, state);
  }

  /* add new failure state */

  sp = &(tbl->stateArray [(int) newState]);
  sp->onfailure = next;

  /* add matches of substring at new state */

  sp = &(tbl->stateArray [(int) next]);
  for (mp = sp->matches; mp != NULL; mp = mp->next) {
    ProtSearchAddOutput (tbl, newState, mp->name, mp->pattern);
  }
}

static void ProtSearchComputeFail (
  ProtSearchPtr tbl,
  Int2Ptr queue
)

{
  CharPtr       charToProt = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  Char          ch;
  Int2          qbeg, r, s, state;
  int           index;
  ProtStatePtr  sp;

  qbeg = 0;
  queue [0] = 0;

  /* queue up states reached directly from state 0 (depth 1) */

  sp = &(tbl->stateArray [0]);
  for (index = 0; index < 26; index++) {
    s = sp->transitions [index];
    if (s == 0) continue;
    sp->onfailure = 0;
    ProtSearchQueueAdd (queue, qbeg, s);
  }

  while (queue [qbeg] != 0) {
    r = queue [qbeg];
    qbeg = r;

    /* depth 1 states beget depth 2 states, etc. */

    sp = &(tbl->stateArray [r]);
    for (index = 0; index < 26; index++) {
      ch = charToProt [index];
      s = sp->transitions [index];
      if (s == 0) continue;
      ProtSearchQueueAdd (queue, qbeg, s);

      state = ProtSearchFailState (tbl, r);
      ProtSearchFindFail (tbl, state, s, ch);
    }
  }
}

/* on first character, populate state transition table */

static void ProtSearchPrimeStateArray (
  ProtSearchPtr tbl
)

{
  Int2            highState, maxState;
  ProtPatternPtr  pp;
  Int2Ptr         queue;
  ProtStatePtr    stateArray;

  if (tbl == NULL || tbl->primed || tbl->patternList == NULL) return;

  for (maxState = 1, pp = tbl->patternList; pp != NULL; pp = pp->next) {
    maxState += StringLen (pp->pattern);
  }

  if (maxState > 4000) {
    Message (MSG_POST, "FiniteStateSearch cannot handle %d states", (int) maxState);
    return;
  }

  stateArray = (ProtStatePtr) MemNew (sizeof (ProtStateItem) * (size_t) maxState);
  queue = (Int2Ptr) MemNew (sizeof (Int2) * maxState);

  if (stateArray == NULL || queue == NULL) {
    MemFree (stateArray);
    MemFree (queue);
    Message (MSG_POST, "SequenceSearch unable to allocate buffers");
    return;
  }

  tbl->stateArray = stateArray;
  tbl->maxState = maxState;

  for (highState = 0, pp = tbl->patternList; pp != NULL; pp = pp->next) {
    highState = ProtSearchEnterProtWord (tbl, highState, maxState, pp->name,
                                         pp->pattern);
  }

  ProtSearchComputeFail (tbl, queue);

  MemFree (queue);

  tbl->highState = highState;
  tbl->currentState = 0;
  tbl->currentPos = 0;
  tbl->primed = TRUE;
}

/* for testing, print summary of transition table */

/*
static void PrintProtSearchTable (
  ProtSearchPtr tbl,
  FILE *fp
)

{
  Int2         i;
  ProtMatchPtr  mp;
  ProtStatePtr  sp;
  Int2         state;

  if (tbl == NULL || fp == NULL) return;
  if (! tbl->primed) {
    ProtSearchPrimeStateArray (tbl);
  }
  if (tbl->stateArray == NULL) return;
  if (tbl->highState > 99) return;

  fprintf (fp, "State Fail A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z\n");

  for (state = 0; state <= tbl->highState; state++) {
    sp = &(tbl->stateArray [(int) state]);
    fprintf (fp, " %3d  %3d", (int) state, (int) sp->onfailure);

    for (i = 0; i < 26; i++) {
      if (sp->transitions [i] != 0) {
        fprintf (fp, "%3d", (int) sp->transitions [i]);
      } else {
        fprintf (fp, "   ");
      }
    }

    for (mp = sp->matches; mp != NULL; mp = mp->next) {
      fprintf (fp, " %s", mp->name);
    }

    fprintf (fp, "\n");
  }
}
*/

/* create empty protein sequence search finite state machine */

NLM_EXTERN ProtSearchPtr ProtSearchNew (
  ProtSearchMatchProc matchproc,
  Pointer userdata
)

{
  CharPtr        charToProt = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  Char           ch;
  Int2           i;
  ProtSearchPtr  tbl;

  if (matchproc == NULL) return NULL;
  tbl = (ProtSearchPtr) MemNew (sizeof (ProtSearchData));
  if (tbl == NULL) return NULL;

  tbl->stateArray = NULL;
  tbl->patternList = NULL;
  tbl->maxPatLen = 0;
  tbl->maxState = 0;
  tbl->highState = 0;
  tbl->currentState = 0;
  tbl->currentPos = 0;
  tbl->matchproc = matchproc;
  tbl->userdata = userdata;
  tbl->primed = FALSE;

  /* initialize table to convert character to transition index from 0 (A) to 25 (Z) */

  for (i = 0; i < 256; i++) {
    tbl->letterToIdx [i] = 23;
  }
  for (i = 0; i < 26; i++) {
    ch = charToProt [i];
    tbl->letterToIdx [(int) ch] = i;
    ch = TO_LOWER (ch);
    tbl->letterToIdx [(int) ch] = i;
  }

  return tbl;
}

/* table to expand ambiguity letter to all matching nucleotide letters */

static CharPtr  protExpandList [26] = {
  "A",
  "DN",
  "C",
  "D",
  "E",
  "F",
  "G",
  "H",
  "I",
  "IL",
  "K",
  "L",
  "M",
  "N",
  "O",
  "P",
  "Q",
  "R",
  "S",
  "T",
  "U",
  "V",
  "W",
  "ACDEFGHIKLMNOPQRSTUVWY",
  "Y",
  "EQ"
};

/* recursive function to expand and store appropriate individual patterns */

static void StoreProtPattern (
  ProtSearchPtr tbl,
  CharPtr name,
  CharPtr str
)

{
  Int4            patLen;
  ProtPatternPtr  pp;

  pp = (ProtPatternPtr) MemNew (sizeof (ProtPatternItem));
  if (pp == NULL) return;

  pp->name = StringSave (name);
  pp->pattern = StringSave (str);

  pp->next = tbl->patternList;
  tbl->patternList = pp;
  patLen = StringLen (str);
  if (patLen > tbl->maxPatLen) {
    tbl->maxPatLen = patLen;
  }
}

static void ExpandProtPattern (
  ProtSearchPtr tbl,
  CharPtr name,
  CharPtr pattern,
  size_t patLen,
  CharPtr str,
  Int2 position,
  SearchFlgType flags
)

{
  Char     ch, lttr;
  Int2     idx;
  CharPtr  ptr;

  if (position < patLen) {

    if ((Boolean) ((flags & SEQ_SEARCH_EXPAND_PATTERN) != 0)) {

      /* given ambiguity letter, get index into protExpandList */

      ch = pattern [position];
      idx = ch - 'A';
      ptr = protExpandList [idx];

      /* put every unambiguous amino acid letter at current
         position, recurse for next position */

      for (lttr = *ptr; lttr != '\0'; ptr++, lttr = *ptr) {
        str [position] = lttr;
        ExpandProtPattern (tbl, name, pattern, patLen, str, position + 1, flags);
      }

    } else {

      /* if matching ambiguity characters in sequence, do not expand each base */

      str [position] = pattern [position];
      ExpandProtPattern (tbl, name, pattern, patLen, str, position + 1, flags);
    }

    /* do not run into pattern storage section of code located below */

    return;
  }

  /* when position reaches pattern length, store one fully expanded string */

  StoreProtPattern (tbl, name, str);

  if ((Boolean) ((flags & SEQ_SEARCH_ALLOW_MISMATCH) == 0)) return;

  for (idx = 0; idx < patLen; idx++) {
    ch = str [idx];

    /* put X at every position if a single mismatch is allowed */

    str [idx] = 'X';

    StoreProtPattern (tbl, name, str);

    /* now restore proper character, go on to put X in next position */

    str [idx] = ch;
  }
}

/* add protein to sequence search finite state machine */

NLM_EXTERN void ProtSearchAddProteinPattern (
  ProtSearchPtr tbl,
  CharPtr name,
  CharPtr pattern,
  SearchFlgType flags
)

{
  Char     ch, pat [128], str [128];
  Int2     i;
  size_t   len;

  if (tbl == NULL || StringHasNoText (name) || StringHasNoText (pattern)) return;

  StringNCpy_0 (pat, pattern, sizeof (pat));
  TrimSpacesAroundString (pat);

  len = StringLen (pat);

  /* upper case working copy of pattern string */

  for (i = 0; i < len; i++) {
    ch = pat [i];
    pat [i] = TO_UPPER (ch);
  }

  /* record expansion of entered pattern */

  MemSet ((Pointer) str, 0, sizeof (str));
  ExpandProtPattern (tbl, name, pat, len, str, 0, flags);
}

/* program passes each character in turn to finite state machine */

static void ProtSearchProcessCharacterEx (
  ProtSearchPtr tbl,
  Char ch,
  Int4 length
)

{
  Int2          curr, next;
  ProtMatchPtr  mp;
  Int4          patLen;
  ProtStatePtr  sp;

  if (tbl == NULL) return;
  if (! tbl->primed) {
    ProtSearchPrimeStateArray (tbl);
  }
  if (tbl->stateArray == NULL) return;

  curr = tbl->currentState;

  /* loop through failure states until match or back to state 0 */

  while ((next = ProtSearchGotoState (tbl, curr, ch, TRUE)) == FAIL_STATE) {
    curr = ProtSearchFailState (tbl, curr);
  }

  tbl->currentState = next;
  (tbl->currentPos)++;

  /* report any matches at current state to callback function */

  sp = &(tbl->stateArray [(int) next]);
  for (mp = sp->matches; mp != NULL; mp = mp->next) {

    /* for circular sequences, prevent multiple reports of patterns */

    patLen = StringLen (mp->pattern);
    if (tbl->currentPos - patLen < length) {
      tbl->matchproc (tbl->currentPos - patLen,
                      mp->name, mp->pattern, tbl->userdata);
    }
  }
}

NLM_EXTERN void ProtSearchProcessCharacter (
  ProtSearchPtr tbl,
  Char ch
)

{
  ProtSearchProcessCharacterEx (tbl, ch, INT4_MAX);
}

/* convenience function calls ProtSearchProcessCharacter for entire protein bioseq */

typedef struct protsrchdata {
  ProtSearchPtr  tbl;
  Int4          length;
} ProtSrchData, PNTR ProtSrchPtr;

static void LIBCALLBACK SearchProtProc (
  CharPtr sequence,
  Pointer userdata
)

{
  Char        ch;
  CharPtr     ptr;
  ProtSrchPtr  ssp;

  ssp = (ProtSrchPtr) userdata;

  ptr = sequence;
  ch = *ptr;
  while (ch != '\0') {
    ch = TO_UPPER (ch);
    ProtSearchProcessCharacterEx (ssp->tbl, ch, ssp->length);
    ptr++;
    ch = *ptr;
  }
}

NLM_EXTERN void ProtSearchProcessBioseq (
  ProtSearchPtr tbl,
  BioseqPtr bsp
)

{
  ProtSrchData  ssd;

  ProtSearchReset (tbl);

  if (tbl == NULL || bsp == NULL) return;

  if (! ISA_aa (bsp->mol)) return;

  ssd.tbl = tbl;
  ssd.length = bsp->length;

  SeqPortStream (bsp, STREAM_EXPAND_GAPS, (Pointer) &ssd, SearchProtProc);

  ProtSearchReset (tbl);
}

/* reset state and position to allow another run with same search patterns */

NLM_EXTERN void ProtSearchReset (
  ProtSearchPtr tbl
)

{
  if (tbl == NULL) return;

  tbl->currentState = 0;
  tbl->currentPos = 0;
}

/* clean up sequence search finite state machine allocated memory */

static ProtPatternPtr FreeProtPatternList (
  ProtPatternPtr pp
)

{
  ProtPatternPtr  next;

  while (pp != NULL) {
    next = pp->next;
    pp->next = NULL;
    MemFree (pp->name);
    MemFree (pp->pattern);
    MemFree (pp);
    pp = next;
  }

  return NULL;
}

static ProtMatchPtr FreeProtMatchList (
  ProtMatchPtr mp
)

{
  ProtMatchPtr  next;

  while (mp != NULL) {
    next = mp->next;
    mp->next = NULL;
    MemFree (mp->name);
    MemFree (mp->pattern);
    MemFree (mp);
    mp = next;
  }

  return NULL;
}

NLM_EXTERN ProtSearchPtr ProtSearchFree (
  ProtSearchPtr tbl
)

{
  Int2  maxState, state;

  if (tbl == NULL) return NULL;

  maxState = tbl->maxState;

  for (state = 0; state < maxState; state++) {
    FreeProtMatchList (tbl->stateArray [state].matches);
  }

  FreeProtPatternList (tbl->patternList);

  MemFree (tbl->stateArray);
  return MemFree (tbl);
}

/*

static CharPtr testseq =
 "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN";

static void MatchProc (Int4 position, CharPtr name, CharPtr pattern, Pointer userdata)

{
  Message (MSG_POST, "Name '%s', Pattern '%s', Position %ld",
           name, pattern, (long) position);
}


extern void TestProtSearch (void);
extern void TestProtSearch (void)

{
  Char           ch;
  CharPtr        ptr;
  ProtSearchPtr  tbl;

  tbl = ProtSearchNew (MatchProc, NULL);
  if (tbl == NULL) return;

  ProtSearchAddProteinPattern (tbl, "AmbiG", "GRATYC", 1, SEQ_SEARCH_EXPAND_PATTERN);

  for (ptr = testseq, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    ProtSearchProcessCharacter (tbl, ch);
  }

  ProtSearchFree (tbl);
}

*/

/*****************************************************************************
*
*  Convenience functions for genome processing use BioseqLockById to get sequence
*  record (perhaps with phrap quality score graphs) so fetching from some network
*  or local server must be enabled, or sequences must already be in memory.
*
*****************************************************************************/

NLM_EXTERN CharPtr GetSequenceByFeature (SeqFeatPtr sfp)

{
  Int4        len;
  CharPtr     str = NULL;
  /*
  Int2        actual, cnt;
  BioseqPtr   bsp;
  SeqPortPtr  spp;
  CharPtr     str = NULL, txt;
  */

  if (sfp == NULL) return NULL;
  len = SeqLocLen (sfp->location);
  if (len > 0 && len < MAXALLOC) {
    str = MemNew (sizeof (Char) * (len + 2));
    if (str != NULL) {
      SeqPortStreamLoc (sfp->location, STREAM_EXPAND_GAPS, (Pointer) str, NULL);
    
#if 0
      spp = SeqPortNewByLoc (sfp->location, Seq_code_iupacna);
      if (spp != NULL) {

        bsp = BioseqFindFromSeqLoc (sfp->location);
        if (bsp != NULL) {
          if (bsp->repr == Seq_repr_delta || bsp->repr == Seq_repr_virtual) {
            SeqPortSet_do_virtual (spp, TRUE);
          }
        }

        cnt = (Int2) MIN (len, 32000L);
        txt = str;
        actual = 1;

        while (cnt > 0 && len > 0 && actual > 0) {
          actual = SeqPortRead (spp, (BytePtr) txt, cnt);
          if (actual < 0) {
            actual = -actual;
            if (actual == SEQPORT_VIRT || actual == SEQPORT_EOS) {
              actual = 1; /* ignore, keep going */
            } else if (actual == SEQPORT_EOF) {
              actual = 0; /* stop */
            }
          } else if (actual > 0) {
            len -= actual;
            txt += actual;
            cnt = (Int2) MIN (len, 32000L);
          }
        }

        SeqPortFree (spp);
      }
#endif
    }
  }
    
  return str;
}

NLM_EXTERN CharPtr GetSequenceByBsp (BioseqPtr bsp)

{
  CharPtr     str = NULL;
  /*
  Int2        actual, cnt;
  Uint1       code = Seq_code_iupacna;
  Int4        len;
  SeqPortPtr  spp;
  CharPtr     str = NULL, txt;
  */

  if (bsp == NULL || bsp->length >= MAXALLOC) return NULL;

  str = MemNew (sizeof (Char) * (bsp->length + 2));
  if (str == NULL) return NULL;

  SeqPortStream (bsp, STREAM_EXPAND_GAPS, (Pointer) str, NULL);

#if 0
  if (ISA_aa (bsp->mol)) {
      code = Seq_code_iupacaa;
  }
  spp = SeqPortNew (bsp, 0, -1, 0, code);
  if (spp == NULL) {
    MemFree (str);
    return NULL;
  }

  if (bsp->repr == Seq_repr_delta || bsp->repr == Seq_repr_virtual) {
    SeqPortSet_do_virtual (spp, TRUE);
  }

  len = bsp->length;
  cnt = (Int2) MIN (len, 32000L);
  txt = str;
  actual = 1;

  while (cnt > 0 && len > 0 && actual > 0) {
    actual = SeqPortRead (spp, (BytePtr) txt, cnt);
    if (actual < 0) {
      actual = -actual;
      if (actual == SEQPORT_VIRT || actual == SEQPORT_EOS) {
        actual = 1; /* ignore, keep going */
      }
    } else if (actual > 0) {
      len -= actual;
      txt += actual;
      cnt = (Int2) MIN (len, 32000L);
    }
  }

  SeqPortFree (spp);
#endif

  return str;
}

NLM_EXTERN CharPtr GetSequenceByIdOrAccnDotVer (SeqIdPtr sip, CharPtr accession, Boolean is_na)

{
  BioseqPtr  bsp;
  SeqIdPtr   deleteme = NULL;
  CharPtr    str = NULL;

  if (sip == NULL) {
    if (StringHasNoText (accession)) return NULL;
    sip = SeqIdFromAccessionDotVersion (accession);
    deleteme = sip; /* allocated seqid, so must later delete it */
  }
  if (sip == NULL) return NULL;

  bsp = BioseqLockById (sip);
  SeqIdFree (deleteme);
  if (bsp == NULL) return NULL;

  if ((ISA_na (bsp->mol) && is_na) || (ISA_aa (bsp->mol) && (! is_na))) {
    if (bsp->length < MAXALLOC) {
      str = GetSequenceByBsp (bsp);
    }
  }

  BioseqUnlock (bsp);
  return str;
}

/* original convenience function now calls more advanced version that can get proteins */

NLM_EXTERN CharPtr GetDNAbyAccessionDotVersion (CharPtr accession)

{
  return GetSequenceByIdOrAccnDotVer (NULL, accession, TRUE);
}


static void FixGapLength (BioseqPtr bsp, Int4 offset, Int4 diff)
{
  CharPtr     extra_ns;
  SeqLocPtr   slp;
  ValNodePtr  align_annot_list, vnp;
  SeqAnnotPtr sanp;
  
  if (bsp == NULL || bsp->id == NULL || diff == 0) return;
  
  align_annot_list = FindAlignSeqAnnotsForBioseq (bsp);

  if (diff > 0)
  {
    extra_ns = (CharPtr)MemNew ((diff + 1) * sizeof (Char));
    if (extra_ns != NULL)
    {
      MemSet (extra_ns, 'N', diff);
      extra_ns [diff] = 0;
      insertchar (extra_ns, offset, bsp->id, bsp->mol, FALSE);
    }
      slp = SeqLocIntNew (offset, offset + diff - 1, Seq_strand_plus, bsp->id);
    for (vnp = align_annot_list; vnp != NULL; vnp = vnp->next)
    {
      sanp = vnp->data.ptrvalue;
      if (sanp != NULL && sanp->type == 2)
      {
        sanp->data = SeqAlignInsertByLoc (slp, sanp->data);
      }
    }
    SeqLocFree (slp);
  }
  else
  {
      slp = SeqLocIntNew (offset, offset - diff - 1, Seq_strand_plus, bsp->id);
    SeqDeleteByLoc (slp, TRUE, FALSE);
 
    for (vnp = align_annot_list; vnp != NULL; vnp = vnp->next)
    {
      sanp = vnp->data.ptrvalue;
      if (sanp != NULL && sanp->type == 2)
      {
        sanp->data = SeqAlignDeleteByLoc (slp, sanp->data);
      }
    }

    SeqLocFree (slp);
  }
}

static Int4 AddSeqLitData (CharPtr str, ValNodePtr PNTR seq_ext)
{
  Int4      len;
  SeqLitPtr slp;
   
  if (StringHasNoText (str)) {
    return 0;
  }
  len = StringLen (str);  
  slp = (SeqLitPtr) MemNew (sizeof (SeqLit));
  if (slp != NULL) {
    slp->length = len;
    ValNodeAddPointer (seq_ext, (Int2) 2, (Pointer) slp);
    slp->seq_data = (SeqDataPtr) BSNew (slp->length);
    slp->seq_data_type = Seq_code_iupacna;
    AddBasesToByteStore ((ByteStorePtr) slp->seq_data, str);
  }
  return len;
}

static Boolean IsGapUnknown (Int4 gap_len,
                             Int4 unknown_gap_size, 
                             Int4 known_gap_size, 
                             Boolean unknown_greater_than_or_equal, 
                             Boolean known_greater_than_or_equal)
{
  Boolean    make_unknown_size = FALSE;
  
  if (gap_len == 0)
  {
    make_unknown_size = FALSE;
  }
  else if (gap_len == unknown_gap_size)
  {
    make_unknown_size = TRUE;
  }
  else if (gap_len == known_gap_size)
  {
    make_unknown_size = FALSE;
  }
  else if (gap_len > unknown_gap_size && unknown_greater_than_or_equal)
  {
    if (!known_greater_than_or_equal)
    {
      make_unknown_size = TRUE;
    }
    else if (unknown_gap_size > known_gap_size)
    {
      make_unknown_size = TRUE;
    }
    else if (gap_len < known_gap_size)
    {
      make_unknown_size = TRUE;
    }
  }
  return make_unknown_size;
}


static Int4 AddGap(Int4 gap_len, 
                  Boolean    make_unknown_size,
                   BioseqPtr bsp,
                   Int4      len,
                   ValNodePtr PNTR seq_ext)
{
  Int4       added_len = 0;
  SeqLitPtr  slp;
  IntFuzzPtr ifp;
    
  if (gap_len > 0) {
    slp = (SeqLitPtr) MemNew (sizeof (SeqLit));
    if (slp != NULL) {
      slp->length = gap_len;
      ValNodeAddPointer (seq_ext, (Int2) 2, (Pointer) slp);
      if (make_unknown_size) {
        ifp = IntFuzzNew ();
        ifp->choice = 4;
        slp->fuzz = ifp;
        if (slp->length != 100) {
          FixGapLength (bsp, len, 100 - slp->length);
             slp->length = 100;
        }
      }
      added_len += slp->length;
    }
  }
  return added_len;
}


static Boolean NeedToConvert 
(CharPtr  bases, 
 Boolean  unknown_greater_than_or_equal,
 Boolean     known_greater_than_or_equal,
 Int4        unknown_gap_size,
 Int4        known_gap_size)
{
  Int4        gap_len;
  CharPtr     cp;

  if (StringHasNoText(bases)) { 
      return FALSE;
  }
  cp = bases;
  while (*cp != '\0') {

    gap_len = StringSpn (cp, "N");
    if (gap_len > 0 ) {
      if ((gap_len == unknown_gap_size
            || (gap_len > unknown_gap_size && unknown_greater_than_or_equal)
            || gap_len == known_gap_size
            || (gap_len > known_gap_size && known_greater_than_or_equal))) {
        return TRUE;
      } else {
        cp += gap_len;
      }
    } else {
      gap_len = StringCSpn (cp, "N");
      cp += gap_len;
    }      
  }
  return FALSE;
}

/*****************************************************************************
*
*   ConvertNsToGaps
*       Assumes string of Ns means a gap of known length
*
*****************************************************************************/

NLM_EXTERN void ConvertNsToGaps (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CharPtr     bases, str, txt;
  Char        ch;
  Int4        len;
  ValNodePtr  seq_ext;
  Boolean     unknown_greater_than_or_equal = FALSE;
  Boolean     known_greater_than_or_equal = FALSE;
  Int4Ptr     gap_sizes;
  Int4        unknown_gap_size = 0;
  Int4        known_gap_size = 0;
  Int4        gap_len;
  Boolean     make_unknown_size;

  if (bsp == NULL || bsp->repr != Seq_repr_raw || ISA_aa (bsp->mol)) return;
  if (userdata == NULL)
  {
    known_greater_than_or_equal = TRUE;
  }
  else
  {
    gap_sizes = (Int4Ptr) userdata;
    unknown_gap_size = gap_sizes[0];
    known_gap_size = gap_sizes[1];
    if (unknown_gap_size < 0)
    {
      unknown_greater_than_or_equal = TRUE;
      unknown_gap_size = 0 - unknown_gap_size;
    }
    if (known_gap_size < 0)
    {
      known_greater_than_or_equal = TRUE;
      known_gap_size = 0 - known_gap_size;
    }
  }

  bases = GetSequenceByBsp (bsp);
  if (bases == NULL) return;

  if (!NeedToConvert(bases, unknown_greater_than_or_equal, known_greater_than_or_equal, unknown_gap_size, known_gap_size)) {
    MemFree (bases);
    return;
  }

  seq_ext = NULL;
  len = 0;

  txt = bases;
  str = txt;
  ch = *txt;
  
  while (*str != '\0') {

    gap_len = StringSpn (str, "N");
    if (gap_len > 0 ) {
      if ((gap_len == unknown_gap_size
            || (gap_len > unknown_gap_size && unknown_greater_than_or_equal)
            || gap_len == known_gap_size
            || (gap_len > known_gap_size && known_greater_than_or_equal))) {
        /* add any prior sequence data as literal */
        ch = *str;
        *str = '\0';
        len += AddSeqLitData (txt, &(seq_ext));
        *str = ch;        
        /* add a gap */
        make_unknown_size = IsGapUnknown (gap_len,
                                          unknown_gap_size, 
                                          known_gap_size, 
                                          unknown_greater_than_or_equal, 
                                          known_greater_than_or_equal);

        len += AddGap(gap_len, 
                      make_unknown_size,
                      bsp,
                      len,
                      &(seq_ext));
        txt = str + gap_len;
      }
      str += gap_len;
    } else {
      gap_len = StringCSpn (str, "N");
      str += gap_len;
    }      
  }
  /* at end, add last sequence data literal */
  len += AddSeqLitData (txt, &(seq_ext));
    
  MemFree (bases);

  bsp->seq_data = SeqDataFree (bsp->seq_data, bsp->seq_data_type);
  bsp->seq_data_type = 0;
  bsp->repr = Seq_repr_delta;
  bsp->seq_ext_type = 4;
  bsp->seq_ext = seq_ext;
  bsp->length = len;

  BioseqPack (bsp);
}


/* Protein Molecular Weight Section */

/* Values are A through Z order:
   B is really D or N, but they are close so is treated as D
   Z is really E or Q, but they are close so is treated as E
   X is hard to guess, so the calculation fails on X
   - and * are skipped
   water molecule is removed for in-peptide atom counts

  A  B  C  D  E  F  G  H   I   J   K   L  M  N   O  P  Q   R  S  T  U  V   W  X  Y  Z
*/
Uint1 C_atoms[26] =
{ 3, 4, 3, 4, 5, 9, 2, 6,  6,  6,  6,  6, 5, 4, 12, 5, 5,  6, 3, 4, 3, 5, 11, 0, 9, 5};
Uint1 H_atoms[26] =
{ 5, 5, 5, 5, 7, 9, 3, 7, 11, 11, 12, 11, 9, 6, 19, 7, 8, 12, 5, 7, 5, 9, 10, 0, 9, 7};
Uint1 N_atoms[26] =
{ 1, 1, 1, 1, 1, 1, 1, 3,  1,  1,  2,  1, 1, 2,  3, 1, 2,  4, 1, 1, 1, 1,  2, 0, 1, 1};
Uint1 O_atoms[26] =
{ 1, 3, 1, 3, 3, 1, 1, 1,  1,  1,  1,  1, 1, 2,  2, 1, 2,  1, 2, 2, 1, 1,  1, 0, 2, 3};
Uint1 S_atoms[26] =
{ 0, 0, 1, 0, 0, 0, 0, 0,  0,  0,  0,  0, 1, 0,  0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0};
Uint1 Se_atoms[26] =
{ 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0,  0, 0, 0,  0, 0, 0,  0, 0, 0, 1, 0,  0, 0, 0, 0};

/**************************************************************
*
*  Returns a protein molecular weight for a SeqLoc
*    If it cannot calculate the value it returns -1.0
*    If sequence contains X,B,U,*,orZ it fails
*
***************************************************************/
NLM_EXTERN FloatHi MolWtForLoc (SeqLocPtr slp)
{
    StreamCache sc;
    Int2 res;
    int residue;
    Int4    Ccnt,
        Hcnt,
        Ncnt,
        Ocnt,
        Scnt,
        Secnt;
    FloatHi retval = -1.0;

    if (slp == NULL) return retval;
    StreamCacheSetup (NULL, slp, 0, &sc);

    Ccnt = 0;  /* initialize counters */
    Hcnt = 2;  /* always start with water */
    Ocnt = 1;  /* H20 */
    Ncnt = 0;
    Scnt = 0;
    Secnt = 0;

    while ((res = StreamCacheGetResidue (&sc)) != '\0')
    {
        if (IS_LOWER (res)) {
            res = TO_UPPER (res);
        }
        if (IS_UPPER (res)) {
            residue = res - 'A';
            if (H_atoms[residue] == 0) { /* unsupported AA */
                return retval;    /* bail out */
            }
            Ccnt += C_atoms[residue];
            Hcnt += H_atoms[residue];
            Ncnt += N_atoms[residue];
            Ocnt += O_atoms[residue];
            Scnt += S_atoms[residue];
            Secnt += Se_atoms[residue];
        } else if (res != '-' && res != '*') {
            return retval;    /* bail out */
        }
    }

    retval = (12.01115 * Ccnt) + (1.0079 * Hcnt) +
         (14.0067 * Ncnt) + (15.9994 * Ocnt) +
         (32.064 * Scnt) + (78.96 * Secnt);

    return retval;
}

NLM_EXTERN FloatHi MolWtForBsp (BioseqPtr bsp)
{
    StreamCache sc;
    Int2 res;
    int residue;
    Int4    Ccnt,
        Hcnt,
        Ncnt,
        Ocnt,
        Scnt,
        Secnt;
    FloatHi retval = -1.0;

    if (bsp == NULL) return retval;
    if (! ISA_aa (bsp->mol)) return retval;
    StreamCacheSetup (bsp, NULL, 0, &sc);

    Ccnt = 0;  /* initialize counters */
    Hcnt = 2;  /* always start with water */
    Ocnt = 1;  /* H20 */
    Ncnt = 0;
    Scnt = 0;
    Secnt = 0;

    while ((res = StreamCacheGetResidue (&sc)) != '\0')
    {
        if (IS_LOWER (res)) {
            res = TO_UPPER (res);
        }
        if (IS_UPPER (res)) {
            residue = res - 'A';
            if (H_atoms[residue] == 0) { /* unsupported AA */
                return retval;    /* bail out */
            }
            Ccnt += C_atoms[residue];
            Hcnt += H_atoms[residue];
            Ncnt += N_atoms[residue];
            Ocnt += O_atoms[residue];
            Scnt += S_atoms[residue];
            Secnt += Se_atoms[residue];
        } else if (res != '-' && res != '*') {
            return retval;    /* bail out */
        }
    }

    retval = (12.01115 * Ccnt) + (1.0079 * Hcnt) +
         (14.0067 * Ncnt) + (15.9994 * Ocnt) +
         (32.064 * Scnt) + (78.96 * Secnt);

    return retval;
}

NLM_EXTERN FloatHi MolWtForStr (CharPtr str)
{
    Char res;
    int residue;
    Int4    Ccnt,
        Hcnt,
        Ncnt,
        Ocnt,
        Scnt,
        Secnt;
    FloatHi retval = -1.0;

    if (str == NULL) return retval;

    Ccnt = 0;  /* initialize counters */
    Hcnt = 2;  /* always start with water */
    Ocnt = 1;  /* H20 */
    Ncnt = 0;
    Scnt = 0;
    Secnt = 0;

    res = *str;
    while (res != '\0')
    {
        if (IS_LOWER (res)) {
            res = TO_UPPER (res);
        }
        if (IS_UPPER (res)) {
            residue = res - 'A';
            if (H_atoms[residue] == 0) { /* unsupported AA */
                return retval;    /* bail out */
            }
            Ccnt += C_atoms[residue];
            Hcnt += H_atoms[residue];
            Ncnt += N_atoms[residue];
            Ocnt += O_atoms[residue];
            Scnt += S_atoms[residue];
            Secnt += Se_atoms[residue];
        } else if (res != '-' && res != '*') {
            return retval;    /* bail out */
        }
        str++;
        res = *str;
    }

    retval = (12.01115 * Ccnt) + (1.0079 * Hcnt) +
         (14.0067 * Ncnt) + (15.9994 * Ocnt) +
         (32.064 * Scnt) + (78.96 * Secnt);

    return retval;
}

