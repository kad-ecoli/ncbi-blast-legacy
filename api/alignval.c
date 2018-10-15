/*  alignval.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  alignval.c
*
* Author:  Jian Ye, Colombe Chappey
*
* Version Creation Date:   6/3/99
*
* $Revision: 6.77 $
*
* File Description:  To validate sequence alignment.
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

 
#include <ncbi.h>
#include <seqmgr.h>
#include <objmgr.h>
#include <sequtil.h> 
#include <sqnutils.h>
#include <satutil.h>
#include <salsap.h>
#include <txalign.h>
#include <salpacc.h>
#include <alignval.h>
#include <valid.h>
#include <alignmgr2.h>


Uint1  jybitnum[8]={0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};

typedef struct saval {
  Boolean     message;
  Boolean     msg_success;
  Boolean     find_remote_bsp;
  Boolean     find_acc_bsp;
  Boolean     delete_salp;
  Boolean     delete_bsp;
  Boolean     retdel;
  Boolean     do_hist_assembly;
  ValNodePtr  ids;
  Uint2       entityID;
  Boolean     dirty;
} SaVal, PNTR SaValPtr;

typedef struct JY_error_msg {
        Uint1 level;/* corresponds to levels of ErrPostEx [none(0), info(1), war
n(2), error(3) and fatal(4)] */
        CharPtr msg;
} JYErrorMsg, *JYErrorMsgPtr;

/******************************************************************
***
*** Error Messaging
***    copies of the BLASt functions in blastpri.h
***    JYConstructErrorMessage = BlastConstructErrorMessage
***    JYErrorChainDestroy = BlastErrorChainDestroy
***
******************************************************************/ 

static ValNodePtr errorp = NULL;
#define BUFFER_LENGTH 512

static Uint2 AlignmentPercentIdentityEx (SeqAlignPtr salp, Boolean internal_gaps, Boolean internal_validation);

static ValNodePtr JYConstructErrorMessage (CharPtr function, CharPtr message, Uint1 level, ValNodePtr PNTR vnpp)
{
        Char buffer[BUFFER_LENGTH];
        CharPtr ptr;
        JYErrorMsgPtr error_msg;

        if (vnpp == NULL)
                return NULL;

        buffer[0] = NULLB;
        ptr = buffer;
        if (function != NULL)
        {
                sprintf(buffer, "%s: ", function);
                ptr = buffer + StringLen(buffer);
        }

        if (message != NULL)
        {
                sprintf(ptr, "%s", message);
        }

        error_msg = (JYErrorMsgPtr) MemNew(sizeof(JYErrorMsg));
        error_msg->msg = StringSave(buffer);
        error_msg->level = level;

        ValNodeAddPointer(vnpp, 0, error_msg);

        return *vnpp;
}

static ValNodePtr JYErrorChainDestroy (ValNodePtr vnp)

{
        ValNodePtr start = vnp;
        JYErrorMsgPtr error_msg;

        while (vnp)
        {
           error_msg = (JYErrorMsgPtr) vnp->data.ptrvalue;
           if (error_msg != NULL) {
              MemFree(error_msg->msg);
           }
           vnp->data.ptrvalue = MemFree(vnp->data.ptrvalue);
           vnp = vnp->next;
        }

        ValNodeFree(start);

        return NULL;
}
/******************************************************************
Output error message according to code defined in alignval.h.  
id refers to seqid of the sequence that causes the error 
and idcontext refers to other sequences in the same segment.  
Intvalue is used to indicate 1) the segment where the sequence 
with error is, or 2) the segtype in case of segtype error.  
Please note that not all errors report all three 
parameters(id, idcontext, Intvalue)
******************************************************************/ 

static Boolean useValErr = FALSE;
static Boolean useLockByID = FALSE;
static ValidStructPtr useVsp = NULL;

static BioseqPtr AlignValBioseqLockById (SeqIdPtr sid)

{
  Int4 old_sev;
  BioseqPtr bsp = NULL;

  if (useLockByID) {
    old_sev = ErrSetMessageLevel (SEV_WARNING);
    bsp = BioseqLockById (sid);
    ErrSetMessageLevel ((ErrSev) old_sev);
  } else {
    bsp = BioseqFindCore (sid);
  }
  return bsp;
}

static Boolean AlignValBioseqUnlock (BioseqPtr bsp)

{
  if (useLockByID) {
    return BioseqUnlock (bsp);
  } else {
    return TRUE;
  }
}

NLM_EXTERN void CDECL  ValidErr VPROTO((ValidStructPtr vsp, int severity, int code1, int code2, const char *fmt, ...));

/*****************************************************************
*  get the approximate sequence coordinate for an alignment segment
*  sip == NULL -> get alignment coordinate
*****************************************************************/
static Int4 valmsggetseqpos(SeqAlignPtr sap, Int4 segment, SeqIdPtr sip)
{
   Int4          c;
   DenseDiagPtr  ddp;
   DenseSegPtr   dsp;
   Boolean       found;
   Int4          i;
   Int4          j;
   Int4          pos;
   PackSegPtr    psp;
   Uint1Ptr      seqpresence;
   SeqIdPtr      sip_tmp;
   SeqLocPtr     slp;
   StdSegPtr     ssp;

   if (sap == NULL || sap->segs == NULL) {
      return -1;
   } else if (segment == 0) {
      return 0;
   }
   if (sap->segtype == SAS_DENSEG)
   {
      dsp = (DenseSegPtr)sap->segs;
      if (sip == NULL)
      {
         pos = 0;
         for (c=0; c<segment; c++)
         {
            pos += dsp->lens[c];
         }
         return pos;
      }
      sip_tmp = dsp->ids;
      i = 0;
      found = FALSE;
      while (!found && sip_tmp != NULL)
      {
         if (SeqIdComp(sip, sip_tmp) == SIC_YES)
            found = TRUE;
         else
         {
            sip_tmp = sip_tmp->next;
            i++;
         }
      }
      if (!found || i>dsp->dim || segment > dsp->numseg)
         return -1;
      pos = 0;
      for (c=0; c<segment; c++)
      {
         if ((j = dsp->starts[(dsp->dim*c)+i])>0)
            pos=j;
      }
      return pos;
   } else if (sap->segtype == SAS_DENDIAG)
   {
      ddp = (DenseDiagPtr)sap->segs;
      pos = 0;
      for (c=0; c<segment; c++)
      {
         pos += ddp->len;
         ddp = ddp->next;
         if (ddp == NULL)
            return -1;
      }
      if (sip == NULL)
         return pos;
      sip_tmp = ddp->id;
      i = 0;
      found = FALSE;
      while (!found && sip_tmp != NULL)
      {
         if (SeqIdComp(sip, sip_tmp) == SIC_YES)
            found = TRUE;
         else
         {
            sip_tmp = sip_tmp->next;
            i++;
         }
      }
      if (!found || i>ddp->dim)
         return -1;
      return (ddp->starts[i]);
   } else if (sap->segtype == SAS_STD)
   {
      ssp = (StdSegPtr)(sap->segs);
      pos = 0;
      for (c=0; c<segment-1; c++)
      {
         pos += SeqLocLen(ssp->loc);
         ssp = ssp->next;
         if (ssp == NULL)
            return -1;
      }
      if (sip == NULL)
         return pos;
      slp = ssp->loc;
      found = FALSE;
      while (!found && slp!=NULL)
      {
         sip_tmp = SeqLocId(slp);
         if (SeqIdComp(sip, sip_tmp) == SIC_YES)
            found = TRUE;
         else
            slp = slp->next;
      }
      if (!found)
         return -1;
      return (SeqLocStart(slp));
   } else if (sap->segtype == SAS_PACKED)
   {
      psp = (PackSegPtr)(sap->segs);
      if (segment > psp->numseg)
         return -1;
      if (sip == NULL)
      {
         pos = 0;
         for (c=0; c<segment; c++)
         {
            pos += psp->lens[c];
         }
         return pos;
      }
      sip_tmp = psp->ids;
      i = 0;
      found = FALSE;
      while (!found && sip_tmp != NULL)
      {
         if (SeqIdComp(sip, sip_tmp) == SIC_YES)
            found = TRUE;
         else
         {
            sip_tmp = sip_tmp->next;
            i++;
         }
      }
      if (!found || i>psp->dim)
         return -1;
      pos = 0;
      seqpresence = NULL;
      BSSeek(psp->present, 0, SEEK_SET);
      seqpresence=MemNew(BSLen(psp->present));
      if(!seqpresence)
         return -1;
      BSRead(psp->present, seqpresence, BSLen(psp->present));
      for (c=0; c<segment; c++)
      {
         if (seqpresence[(c*psp->numseg+i)/8]&jybitnum[(c*psp->numseg+i)%8])
            pos+=psp->lens[c];
      }
      return pos;
   } else
      return -1;
}


static BioseqPtr BioseqForAlignment (SeqAlignPtr salp)
{
  Int4 row, num_rows;
  BioseqPtr bsp = NULL;
  SeqIdPtr  sip;
  SeqEntryPtr oldscope;
  DenseDiagPtr ddp;
  
  oldscope = SeqEntrySetScope (NULL);
  /* NOTE - can't index DenseDiag chain during validation because we're examining the individual DenseDiags,
   * and indexing converts it to DenseSegs.
   */
  if (salp->segtype == SAS_DENDIAG && salp->segs != NULL) {
    ddp = (DenseDiagPtr) salp->segs;
    while (bsp == NULL && ddp != NULL) {
      for (sip = ddp->id; bsp == NULL && sip != NULL; sip = sip->next) {
        bsp = BioseqFind (sip);
        sip = sip->next;
      }
      ddp = ddp->next;
    }
  } else {
    AlnMgr2IndexSingleChildSeqAlign(salp);
    num_rows = AlnMgr2GetNumRows(salp);
    for (row = 1; row <= num_rows && bsp == NULL; row++) {
      sip = AlnMgr2GetNthSeqIdPtr(salp, row);
      bsp = BioseqFind(sip);
    }
  }
  SeqEntrySetScope (oldscope);  
  return bsp;
}


static void ValMessage (SeqAlignPtr salp, Int1 MessageCode, ErrSev errlevel, SeqIdPtr id, SeqIdPtr idcontext , Int4 Intvalue) 
{
  
  Char     buf[256], 
           buf3[64],
           string1[64],
           string2[552];
  GatherContextPtr gcp;
  Int4     pos;

  string1[0] = '\0';
  string2[0] = '\0';
  SeqIdWrite(id, buf, PRINTID_FASTA_LONG, sizeof(buf)-1);
  switch(MessageCode)
  {
    case Err_SeqId:
      sprintf(string1, "SeqId");
      sprintf(string2, "The sequence corresponding to SeqId %s could not be found", buf);
      break;

    case Err_Strand_Rev:
      pos = valmsggetseqpos(salp, Intvalue, id);
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Strand");
      sprintf(string2, "The strand labels for SeqId %s are inconsistent across the alignment; the first inconsistent region is the %ld(th) region, near sequence position %ld, context %s", buf, (long) Intvalue, (long) pos, buf3);
      break;

    case Err_Denseg_Len_Start:
      pos = valmsggetseqpos(salp, Intvalue, id);
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Start/Length");
      sprintf(string2, "There is a problem with sequence %s, in segment %ld (near sequence position %ld), context %s: the segment is too long or short or the next segment has an incorrect start position", buf, (long) Intvalue, (long) pos, buf3);
      break;

    case  Err_Start_Less_Than_Zero:
      pos = valmsggetseqpos(salp, Intvalue, id);
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Start");
      sprintf(string2, "Start point is less than zero in segment %ld (near sequence position %ld) for sequence ID: %s in the context of %s", (long) Intvalue, (long) pos, buf, buf3);
      break;

    case Err_Start_More_Than_Biolen:
      pos = valmsggetseqpos(salp, Intvalue, id);
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Start");
      sprintf(string2, "In sequence %s, segment %ld (near sequence position %ld) context %s, the alignment claims to contain residue coordinates that are past the end of the sequence.  Either the sequence is too short, or there are extra characters or formatting errors in the alignment", buf, (long) Intvalue, (long) pos, buf3);
      break;

    case Err_End_Less_Than_Zero:
      pos = valmsggetseqpos(salp, Intvalue, id);
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Length");
      sprintf(string2, "End point is less than zero in segment %ld (near position %d) for sequence ID: %s in the context of %s.  This could be a formatting error", (long) Intvalue, (int) pos,buf, buf3);
      break;

    case Err_End_More_Than_Biolen:
      pos = valmsggetseqpos(salp, Intvalue, id);
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Length");
      sprintf(string2, "In sequence %s, segment %ld (near sequence position %ld) context %s, the alignment claims to contain residue coordinates that are past the end of the sequence.  Either the sequence is too short, or there are extra characters or formatting errors in the alignment", buf, (long) Intvalue, (long) pos, buf3);
      break;

    case Err_Len_Less_Than_Zero:
      pos = valmsggetseqpos(salp, Intvalue, id);
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Length");
      sprintf(string2, "Segment length is less than zero in segment %ld (near sequence position %ld) for sequence ID: %s in the context of %s.  Look for extra characters in this segment or flanking segments", (long) Intvalue, (long) pos, buf, buf3); 
      break;

    case Err_Len_More_Than_Biolen:
      pos = valmsggetseqpos(salp, Intvalue, id);
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Length");
      sprintf(string2, "In sequence %s, segment %ld (near sequence position %ld) context %s, the alignment claims to contain residue coordinates that are past the end of the sequence.  Either the sequence is too short, or there are extra characters or formatting errors in the alignment", buf, (long) Intvalue, (long) pos, buf3);
      break; 
 
    case Err_Sum_Len_Start:
      pos = valmsggetseqpos(salp, Intvalue, id);
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Start");
      sprintf(string2, "In sequence %s, segment %ld (near sequence position %ld) context %s, the alignment claims to contain residue coordinates that are past the end of the sequence.  Either the sequence is too short, or there are extra characters or formatting errors in the alignment", buf, (long) Intvalue, (long) pos, buf3);
      break;

    case Err_SeqAlign_DimSeqId_Not_Match:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "SeqId");
      sprintf(string2, "The Seqalign has more or fewer ids than the number of rows in the alignment (context %s).  Look for possible formatting errors in the ids.", buf3);
      break;

    case Err_Segs_DimSeqId_Not_Match:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "SeqId");
      sprintf(string2, "In segment %ld, there are more or fewer rows than there are seqids (context %s).  Look for possible formatting errors in the ids.", (long) Intvalue, buf3);
      break;

    case Err_Fastalike:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Fasta");
      sprintf(string2, "This may be a fasta-like alignment for SeqId: %s in the context of %s", buf, buf3); 
      break;

    case Err_Null_Segs:
      sprintf(string1, "Segs");
      sprintf(string2, "This alignment is missing all segments.  This is a non-correctable error -- look for serious formatting problems.");
      break;

    case Err_Segment_Gap:
      pos = valmsggetseqpos(salp, Intvalue, id);
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Segs");
      sprintf(string2, "Segment %ld (near alignment position %ld) in the context of %s contains only gaps.  Each segment must contain at least one actual sequence -- look for columns with all gaps and delete them.", (long) Intvalue + 1, (long) pos, buf3);
      break;

    case Err_Segs_Dim_One:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Segs");
      sprintf(string2, "Segment %ld apparently has only one sequence.  Each portion of the alignment must have at least two sequences.  context %s", (long) Intvalue, buf3);
      break;

    case Err_SeqAlign_Dim_One:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Dim");
      sprintf(string2, "This seqalign apparently has only one sequence.  Each alignment must have at least two sequences.  context %s", buf3);
      break;

    case Err_Segtype :
      /* ignore new segtype warnings in genomic gpipe sequence */
      if (useValErr && useVsp != NULL && useVsp->is_gpipe_in_sep && useVsp->bsp_genomic_in_sep) return;
      sprintf(string1, "Segs");
      sprintf(string2, "This alignment has an undefined or unsupported Seqalign segtype %ld", (long) Intvalue);
      break;
      
    case Err_Pcnt_ID :
      sprintf(string1, "PercentIdentity");
      sprintf(string2, "This alignment has a percent identity of %d%%", Intvalue);
      break;

    case Err_Short_Aln:
      sprintf(string1, "ShortAln");
      sprintf(string2, "This alignment is shorter than at least one non-farpointer sequence.");
      break;

    case Err_Unexpected_Alignment_Type:
      sprintf(string1, "UnexpectedAlignmentType");
      sprintf (string2, "This is not a DenseSeg alignment.");
      break;

    default:
      break;
  }
  if (useValErr) {
    if (salp != NULL && useVsp != NULL) {
      gcp = useVsp->gcp;
      if (gcp != NULL) {
          gcp->entityID = salp->idx.entityID;
          gcp->itemID = salp->idx.itemID;
          gcp->thistype = salp->idx.itemtype;

        useVsp->bsp = BioseqForAlignment(salp);
        ValidErr (useVsp, errlevel, 6, MessageCode, "%s: %s", string1, string2);
      }
    }
    return;
  }
  if (StringLen(string1) > 0)
     errorp = JYConstructErrorMessage (string1, string2, errlevel, &errorp);
}

 
/******************************************************************
return the number of seqid
******************************************************************/ 
static Int2 CountSeqIdInSip (SeqIdPtr sip)
{
    Int2 numids=0;

     while(sip) 
       { 
     numids++;
     sip=sip->next;
       }
     return numids;
}

/*********************************************************/
static void delete_bioseqs (ValNodePtr ids, Uint2 entityID)
{
  SeqEntryPtr  sep_top;
  SeqEntryPtr  sep_del;
  ValNodePtr   vnp;
  SeqIdPtr     sip;
  SeqLocPtr    slp;
  BioseqPtr    bsp;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;

  if (ids == NULL)
     return;
  sep_top = GetTopSeqEntryForEntityID (entityID);
  SaveSeqEntryObjMgrData (sep_top, &omdptop, &omdata);
  GetSeqEntryParent (sep_top, &parentptr, &parenttype);

  vnp=ids;
  while (vnp!=NULL)
  {
     sip = (SeqIdPtr) vnp->data.ptrvalue;
     if (sip!=NULL) {
        slp = (SeqLocPtr)ValNodeNew (NULL);
        slp->choice = SEQLOC_WHOLE;
        slp->data.ptrvalue = sip;
        bsp = GetBioseqGivenSeqLoc (slp, entityID);
        if (bsp!=NULL) {
           sep_del=GetBestTopParentForData (entityID, bsp);
           RemoveSeqEntryFromSeqEntry (sep_top, sep_del, FALSE);
        }
        slp->data.ptrvalue = NULL;
        SeqLocFree (slp);
     }
     vnp=vnp->next;
  }
  SeqMgrLinkSeqEntry (sep_top, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep_top, omdptop, &omdata);
  RenormalizeNucProtSets (sep_top, TRUE);

  for (vnp=ids; vnp!=NULL; vnp=vnp->next) {
     SeqIdFree ((SeqIdPtr) vnp->data.ptrvalue);
     vnp->data.ptrvalue = NULL;
  }
  ValNodeFree (vnp);
  return;
}


/******************************************************************
validate a SeqId
******************************************************************/ 
static void ValidateSeqId (SeqIdPtr sip, SeqAlignPtr salp)
{
  SeqIdPtr  siptemp=NULL, sipnext;
  BioseqPtr bsp=NULL;
  
  for(siptemp=sip; siptemp!=NULL; siptemp=siptemp->next)
  {
    /*
    bsp = AlignValBioseqLockById(siptemp);
    if(!bsp)
        ValMessage (salp, Err_SeqId, SEV_ERROR, siptemp, NULL, 0);
    else
        AlignValBioseqUnlockById(siptemp);
    */
    sipnext = siptemp->next;
    siptemp->next = NULL;
    bsp = BioseqFindCore (siptemp);
    if (bsp == NULL && siptemp->choice == SEQID_LOCAL) {
        ValMessage (salp, Err_SeqId, SEV_ERROR, siptemp, NULL, 0);
    }
    siptemp->next = sipnext;
  }
  return;
}

/******************************************************************
return seqid for each seg.  
Note that a newly created seqid chain is returned for stdseg 
and you need to free the memory after you use it in this case
******************************************************************/ 
static SeqIdPtr SeqIdInAlignSegs(Pointer segs, Uint1 segtype, SeqAlignPtr salp)
{

  SeqIdPtr sip=NULL;
  StdSegPtr ssp;
  DenseDiagPtr ddp;
  DenseSegPtr dsp;
  PackSegPtr psp;
  SeqLocPtr slp=NULL, slptemp;

  if(!segs)
  {
      ValMessage (salp, Err_Null_Segs, SEV_ERROR, NULL, NULL, 0);
      return NULL;
  }
  if(segtype==1) 
  { /* DenseDiag */
      
      ddp=(DenseDiagPtr)segs;    
      sip=ddp->id;
  }
  else if (segtype==2)
  { /* DenseSeg */
      
      dsp = (DenseSegPtr) segs;
      sip=dsp->ids;
  }
  else if (segtype==3)
  { /* StdSeg */
      
      ssp = (StdSegPtr)segs;
      slp = ssp->loc;
      /*make a new linked list of SeqId*/
      for(slptemp=slp; slptemp!=NULL; slptemp=slptemp->next)
        AddSeqId(&sip, SeqLocId(slptemp));
      
  }
  else if(segtype==4)
  { /* Packed Seg. Optimal for editing alignments */
      
      psp = (PackSegPtr)segs;
      if (psp!=NULL)
        sip = psp->ids;
  }      
  return sip;
}

 
/******************************************************************
validate SeqId in sequence alignment
******************************************************************/ 
static void  ValidateSeqIdInSeqAlign (SeqAlignPtr salp)
{
  SeqIdPtr sip=NULL;
  Pointer segptr=NULL;
  DenseDiagPtr ddp=NULL, ddptemp;
  StdSegPtr    ssp=NULL, ssptemp;
 

  if(salp)
    {     
      segptr=salp->segs;
      if(!segptr)
    ValMessage (salp, Err_Null_Segs, SEV_ERROR, NULL, NULL, 0);
      else
    {

      /*densediag */
      if(salp->segtype==1)
        {
          /*cast to appropriate pointer*/
          ddp=(DenseDiagPtr)segptr;
          for(ddptemp=ddp; ddptemp!=NULL; ddptemp=ddptemp->next)
        {
          
          sip=SeqIdInAlignSegs((Pointer)ddptemp, salp->segtype, salp);    
          ValidateSeqId(sip, salp);
        }
        }
      
      /*Stdseg*/
      else if(salp->segtype==3)
        {
          /*cast to appropriate pointer*/
          ssp=(StdSegPtr)segptr;
          for(ssptemp=ssp; ssptemp!=NULL; ssptemp=ssptemp->next)
        {
          
          sip=SeqIdInAlignSegs((Pointer)ssptemp, salp->segtype, salp);    
          ValidateSeqId(sip, salp);
          /*free Seqid if sip is a new chain created by SeqIdinAlignSegs*/
          SeqIdSetFree(sip);
        }
        }
      
      /*Denseseg, Packseg*/
      else if(salp->segtype==2||salp->segtype==4)
        {
          
          sip=SeqIdInAlignSegs(segptr, salp->segtype, salp);    
          ValidateSeqId(sip, salp);
        } 
    }
    }
}

/******************************************************************
return true if  two sip are the same, false otherwise.  
Also return false if there is error in sip
******************************************************************/ 
static Boolean SeqIdCmp (SeqIdPtr sip1, SeqIdPtr sip2)
{
  Char buf1[256], buf2[256];
 
  if(!sip1||!sip2)
    return FALSE;

  SeqIdWrite(sip1, buf1, PRINTID_FASTA_LONG, 255);
  SeqIdWrite(sip2, buf2, PRINTID_FASTA_LONG, 255);
  return(!StringCmp(buf1, buf2));
 
}
 

/******************************************************************
return the strand for a seqloc with seqid=sip in a stdseg.  
Note, it returns 255 if null sip or ssp
******************************************************************/ 
static Uint1 SeqLocStrandForSipInStdSeg (SeqIdPtr sip, StdSegPtr ssp, SeqAlignPtr salp)
{
  SeqLocPtr slp, slptemp;
  Uint1     strand=0;
    
  if(!sip||!ssp)
    return (255);

  slp=ssp->loc;
  for(slptemp=slp; slptemp!=NULL; slptemp=slptemp->next)
  {
      if(SeqIdCmp(sip, SeqLocId(slptemp)))
    {
      strand=SeqLocStrand(slptemp);
      break;
    }
  }
  return strand;
}


/******************************************************************
check if the  strand is consistent in Stdseg
******************************************************************/ 
static void ValidateStrandInStdSeg(StdSegPtr ssp, SeqAlignPtr salp)
{
  SeqIdPtr     sip=NULL,  sip_inseg=NULL;
  Uint1           strand1=0, strand2=0;
  StdSegPtr    ssptemp, ssptemp2, ssptemp3;
  SeqLocPtr    slp, slptemp;
  ValNodePtr   FinishedSip=NULL, temp;
  Boolean      CheckedStatus;
  Int4         start_numseg=0, end_numseg=0;
  
  if(!ssp)
    ValMessage (salp, Err_Null_Segs, SEV_ERROR, NULL, NULL, 0);
  else
    for(ssptemp=ssp; ssptemp!=NULL; ssptemp=ssptemp->next)
      {
    sip_inseg=SeqIdInAlignSegs((Pointer)ssptemp, 3, salp);
    start_numseg++;
    slp=ssptemp->loc;
    for(slptemp=slp; slptemp!=NULL; slptemp=slptemp->next)
      {
        
        CheckedStatus=FALSE;
        sip=SeqLocId(slptemp);
        if(sip)
          {
        /*if a seqloc represented by a sip has been checked, set the checkedstatus flag to true so it will not be checked again*/
        for(temp=FinishedSip; temp!=NULL; temp=temp->next)
          {
            if(SeqIdCmp(sip, temp->data.ptrvalue))
              {
            CheckedStatus=TRUE;
            break;
              }
          }
        /*seqloc not checked yet*/
        if(!CheckedStatus)
          {
            
            /*keep a record of  checked sip*/
            ValNodeAddPointer(&FinishedSip, 0, sip);
            end_numseg=start_numseg;
            /*go through all segs to get at least two strand, if any, for this seqloc*/
            for(ssptemp2=ssptemp; ssptemp2!=NULL; ssptemp2=ssptemp2->next, end_numseg++)
              {
            /*get the first defined strand */
            strand1=SeqLocStrandForSipInStdSeg(sip, ssptemp2, salp);
            
            if(strand1!=0&&strand1!=255)
              {
                ssptemp2=ssptemp2->next;
                break;
              }
            
              }
            
            if(strand1!=0&&strand1!=255)
              /*continue to get next strand */
              for(ssptemp3=ssptemp2; ssptemp3!=NULL; ssptemp3=ssptemp3->next, end_numseg++)
            {
              strand2=SeqLocStrandForSipInStdSeg(sip, ssptemp3, salp);
              if(strand2==0||strand2==255)
                continue;
              
              if(strand2!=0&&strand2!=255)
                /*strand should be same for a given seq*/ 
                if(strand1!=strand2)
                  
                  ValMessage (salp, Err_Strand_Rev, SEV_ERROR, sip, sip_inseg, end_numseg+1);
            }
            }
          }
      }
    SeqIdSetFree(sip_inseg);
    
      }
  
  ValNodeFree(FinishedSip);
}
 
 
/******************************************************************
check if the  strand is consistent in Denseseg
******************************************************************/ 
static void ValidateStrandInPack_DenseSeg(Pointer segs, Uint1 segtype, SeqAlignPtr salp)
{ 
  DenseSegPtr dsp=NULL;
  PackSegPtr psp=NULL;
  Int4         numseg, aligndim, dimnumseg, i, j, m;
  SeqIdPtr     sip=NULL, siptemp;
  Uint1           strand1=0, strand2=0;
  Uint1Ptr strandptr=NULL;
        
  if(!segs)
  {
    ValMessage (salp, Err_Null_Segs, SEV_ERROR, NULL, NULL, 0);
  } 
  else if(segtype==2||segtype==4)
  {
    if(segtype==2)
    {
      dsp=(DenseSegPtr)segs;
      strandptr=dsp->strands;
      sip=dsp->ids;
      numseg=dsp->numseg;
      aligndim=dsp->dim;
    }     
    else if(segtype==4)
    {
      psp=(PackSegPtr)segs;
      strandptr=psp->strands;
      sip=psp->ids;
      numseg=psp->numseg;
      aligndim=psp->dim;
    }

    dimnumseg=numseg*aligndim;
    if(strandptr)
    {     
      /*go through id for each alignment sequence*/
      for(j=0; j<aligndim; j++)
      {
        /* first  strand value for each sequence*/ 
        strand1=strandptr[j];
        /* go through all strand values for each sequence*/  
        for(i=j+aligndim; i<dimnumseg; i=i+aligndim)
        {          
          strand2=strandptr[i];
          
          if(strand1==0||strand1==255)
          {
            strand1=strand2;
            continue;
          }
          
          /*skip undefined strand*/
          if(strand2!=0&&strand2!=255) 
          {
            /*strand should be same for a given seq*/ 
            if(strand1!=strand2)
            {
              /*find current seqid*/
            
              siptemp=sip;
              for(m=0; m<j&&siptemp!=NULL; m++)
              {
                siptemp=siptemp->next;
              }
              ValMessage (salp, Err_Strand_Rev, SEV_ERROR, siptemp, sip, i/aligndim+1);
            }
          }
        }
      }
    }
  }
}




/******************************************************************
check if the  strand is consistent in SeqAlignment of global 
or partial type
******************************************************************/ 
static void ValidateStrandinSeqAlign(SeqAlignPtr salp)
{
  StdSegPtr ssp=NULL ;
  
  if(salp)
    {
   
      /*Strands needs to be validated  in case of global or partial alignment*/ 
     
      /*denseseg or packseg*/
      if(salp->segtype==2||salp->segtype==4)
 
    ValidateStrandInPack_DenseSeg(salp->segs, salp->segtype, salp);

      /*stdseg*/
      else if(salp->segtype==3)
    {
      ssp=(StdSegPtr)salp->segs;
      ValidateStrandInStdSeg(ssp, salp);
    }
   } 
}



/******************************************************************
Make sure that, in Densediag alignment, segment length and 
start point is not less than zero, and  segment length is not greater 
than Bioseq length
******************************************************************/ 
static void ValidateSeqlengthInDenseDiag (DenseDiagPtr ddp, SeqAlignPtr salp)
{
  Int4Ptr      stptr=NULL; 
  DenseDiagPtr ddptemp;
  Int2         numseg, i;
  SeqIdPtr     sip=NULL, siptemp;
  Int4         bslen;
  BioseqPtr    bsp=NULL;
  

  if(!ddp)
    ValMessage (salp, Err_Null_Segs, SEV_ERROR, NULL, NULL, 0);
  else
    {
      for(ddptemp=ddp, numseg=0; ddptemp!=NULL; ddptemp=ddptemp->next, numseg++)
    {
      sip=ddp->id;
      stptr=ddptemp->starts;
      
      if(stptr)
        {
          for(i=0, siptemp=sip; i<ddptemp->dim; i++, siptemp=siptemp->next)
        {
          bsp=AlignValBioseqLockById(siptemp);
          if(bsp)
            {
              bslen=bsp->length; 
              AlignValBioseqUnlock (bsp);
              /*verify start*/
              if(stptr[i]<0)
            ValMessage (salp, Err_Start_Less_Than_Zero, SEV_ERROR, siptemp, sip , numseg);     
              if(stptr[i]>=bslen)
            ValMessage (salp, Err_Start_More_Than_Biolen, SEV_ERROR, siptemp, sip , numseg); 
              
              /*verify length*/
              
              if(ddptemp->len<0)
            ValMessage (salp, Err_Len_Less_Than_Zero, SEV_ERROR, siptemp, sip , numseg); 
              
              if(ddptemp->len+stptr[i]>bslen)
            ValMessage (salp, Err_Sum_Len_Start, SEV_ERROR, siptemp, sip , numseg);  
            }
        }
        }
    }
    }
}


/******************************************************************
return a new copy of len array in reversed order 
******************************************************************/ 
static Int4Ptr GetReverseLength (Int2 numseg, Int4Ptr lenptr)
{
  Int4Ptr lenptrtemp=NULL;
  Int2 p;
  
  if(!lenptr)
    return NULL;

  lenptrtemp=(Int4Ptr)MemNew(numseg*sizeof(Int4Ptr));
  if(!lenptrtemp)
  {
      ErrPostEx (SEV_ERROR, 0,0,  "Warning:insufficient memory");
      return NULL;
  }
  for(p=0; p<numseg; p++)    
    lenptrtemp[p]=lenptr[numseg-1-p];
  return lenptrtemp;

}

/******************************************************************
return a new copy of start array in reversed "numseg" order .  
Note that the relative position of starts in each numseg has not changed.  
Example:  original length={0, 0, 10, -1, 30, 10}, numseg=3, 
lens={10, 20, 40}, the reversed length={30, 10, 10, -1, 0, 0}
******************************************************************/ 
static Int4Ptr GetReverseStart(Int2 numseg, Int2 dim, Int4Ptr stptr)
{
  Int4Ptr stptrtemp=NULL;
  Int2 p, q;

  if(!stptr)
    return NULL;

  stptrtemp=(Int4Ptr)MemNew(numseg*dim*sizeof(Int4Ptr));
  if(!stptrtemp)
  {
      ErrPostEx (SEV_ERROR, 0,0,  "Warning:insufficient memory"); 
      return NULL; 
  }
  for(p=0; p<numseg; p++)
    for(q=0; q<dim; q++)
      stptrtemp[q+p*dim]=stptr[q+(numseg-1-p)*dim];

  return stptrtemp;
}

 

/******************************************************************
Make sure that, in Denseseg alignment, segment length and 
start point agrees each other and the sum of segment length 
is not greater than Bioseq length
******************************************************************/ 
static void ValidateSeqlengthInDenseSeg (DenseSegPtr dsp, SeqAlignPtr salp)
{

  Int4Ptr      lenptr=NULL, stptr=NULL, lenptrtemp=NULL, stptrtemp=NULL, lenptrtemp2=NULL, stptrtemp2=NULL;
  
  Int2         numseg, aligndim, i, j;
  SeqIdPtr     sip=NULL, siptemp;
  Int4         bslen = 0;
  BioseqPtr    bsp=NULL;

 if(!dsp)
   ValMessage (salp, Err_Null_Segs, SEV_ERROR, NULL, NULL, 0);
 else
    {
      numseg=dsp->numseg;
      aligndim=dsp->dim;     
      
      stptr=dsp->starts;
      lenptr=dsp->lens;
      sip=dsp->ids;
     
      if(stptr==NULL||lenptr==NULL)
    return;
      
  
      /*go through each sequence*/
      for(j=0, siptemp=sip; j<aligndim&&siptemp; j++, siptemp=siptemp->next)
    {
       
      lenptrtemp=lenptr;
      stptrtemp=stptr;
      /*if on minus strand, use reversed length and start array*/
      if(dsp->strands)
        {
          if(dsp->strands[j]==Seq_strand_minus)
        {
          if(!lenptrtemp2&&!stptrtemp2)
            {
              lenptrtemp2= GetReverseLength (numseg, lenptr);
              if (lenptrtemp2==NULL)
                 return;
              stptrtemp2= GetReverseStart (numseg, aligndim, stptr);
              if (stptrtemp2==NULL)
                 return;
            }
          lenptrtemp=lenptrtemp2;
          stptrtemp=stptrtemp2;
        }
        }

      bsp=AlignValBioseqLockById(siptemp);
      if(bsp!=NULL)
        {
          bslen=bsp->length;  
          AlignValBioseqUnlock (bsp);
        }

      /*go through each segment for a given sequence*/
      for(i=0; i<numseg; i++)
        {
       
          /*no need to verify if segment is not present*/
          if(stptrtemp[j+i*aligndim]!=-1)
        {
 
          /*length plus start should be equal to next start*/
          /*check a start if it's not the last one and the next start is not -1*/
          if(i!=numseg-1&&stptrtemp[j+(i+1)*aligndim]!=-1)
            {      
              
              if(stptrtemp[j+i*aligndim]+lenptrtemp[i]!=stptrtemp[j+(i+1)*aligndim]) 
            {
              if (dsp->strands)
                {
                  if(dsp->strands[j]==2)
                ValMessage (salp, Err_Denseg_Len_Start, SEV_ERROR, siptemp, sip , numseg-i); 
                  else
                ValMessage (salp, Err_Denseg_Len_Start, SEV_ERROR, siptemp, sip , i+1);    
                }
                          else
                ValMessage (salp, Err_Denseg_Len_Start, SEV_ERROR, siptemp, sip , i+1);
            }
            }
          /*check a start if it's not the last one and the next start is -1*/
          else if (i!=numseg-1&&stptrtemp[j+(i+1)*aligndim]==-1)
            {
              Int4 k=i+1;
              /*find the next start that is not last and not -1*/
              while(k<numseg&&stptrtemp[j+k*aligndim]==-1)
            k++;

              /*length plus start should be equal to the closest next start that is not -1*/           
     
              if(k<numseg&&stptrtemp[j+i*aligndim]+lenptrtemp[i]!=stptrtemp[j+k*aligndim])
            {
              if (dsp->strands)
                {
                  if(dsp->strands[j]==2)
                ValMessage (salp, Err_Denseg_Len_Start, SEV_ERROR, siptemp, sip , numseg-i); 
                  else
                 ValMessage (salp, Err_Denseg_Len_Start, SEV_ERROR, siptemp, sip , i+1); 
                }
                          else
                ValMessage (salp, Err_Denseg_Len_Start, SEV_ERROR, siptemp, sip , i+1);    
            }
            }
          
          
         /*make sure the start plus segment does not exceed total bioseq length*/ 
          if(bsp!=NULL)
            {
              
              if(stptrtemp[j+i*aligndim]+lenptrtemp[i]>bslen)
            if (dsp->strands)
              {
                if(dsp->strands[j]==2)
                  ValMessage (salp, Err_Sum_Len_Start, SEV_ERROR, siptemp, sip , numseg-1); 
                else
                  ValMessage (salp, Err_Sum_Len_Start, SEV_ERROR, siptemp, sip , i+1); 
              }
            else
              ValMessage (salp, Err_Sum_Len_Start, SEV_ERROR, siptemp, sip , i+1); 
            }
          
        }
                    
        }        
    }
    }        


 MemFree(lenptrtemp2);
 MemFree(stptrtemp2);
                  

}

/******************************************************************
Make sure that, in Seqloc of a Stdseg alignment, 
end point, start point and length are not less than zero, 
and are not greater than Bioseq length
******************************************************************/ 
static void ValidateSeqlengthInStdSeg (StdSegPtr ssp, SeqAlignPtr salp)
{ 
  StdSegPtr    ssptemp;
  Int2         numseg;
  SeqIdPtr     sip=NULL, siptemp;
  Int4         start, end, length, bslen;
  BioseqPtr    bsp=NULL;
  SeqLocPtr    slp=NULL, slptemp;

  if(!ssp) {
    ValMessage (salp, Err_Null_Segs, SEV_ERROR, NULL, NULL, 0);
  } else {
    for(ssptemp=ssp, numseg=0; ssptemp!=NULL; ssptemp=ssptemp->next, numseg++) { 
      /*get all seqid in current segment*/
      sip=SeqIdInAlignSegs((Pointer)ssptemp, 3, salp);         
      slp=ssptemp->loc;
      if(slp==NULL)
        return;
      for(slptemp=slp; slptemp!=NULL; slptemp=slptemp->next) { 
        siptemp=SeqLocId(slptemp);
        start=SeqLocStart(slptemp);
        end=SeqLocStop(slptemp);
        length=SeqLocLen(slptemp);
        
        bsp=AlignValBioseqLockById(siptemp);
        if(bsp) {
          bslen=bsp->length;
          AlignValBioseqUnlock (bsp);
     
          /*verify start*/
          if(start<0) {
            ValMessage (salp, Err_Start_Less_Than_Zero, SEV_ERROR, siptemp, sip , numseg+1);      
          }
            
          if(start>bslen-1) {
            ValMessage (salp, Err_Start_More_Than_Biolen, SEV_ERROR, siptemp, sip , numseg+1); 
          }
          
            /*verify end*/
          if(end<0) {
            ValMessage (salp, Err_End_Less_Than_Zero, SEV_ERROR, siptemp, sip , numseg+1); 
          }
          if(end>bslen-1) {
            ValMessage (salp, Err_End_More_Than_Biolen, SEV_ERROR, siptemp, sip , numseg+1); 
          }
                                  
          /*verify length*/
          if(length<0) {
            ValMessage (salp, Err_Len_Less_Than_Zero, SEV_ERROR, siptemp, sip , numseg+1); 
          }
            
          if(length>bslen) {
            ValMessage (salp, Err_Len_More_Than_Biolen, SEV_ERROR, siptemp, sip , numseg+1);  
          }
        
        }
      }
      /*free Seqid if sip is a new chain created by SeqIdinAlignSegs*/      
      SeqIdSetFree(sip);
    }
  }
}

/******************************************************************
validate the start and segment length in packseg
******************************************************************/ 
static void ValidateSeqlengthInPackSeg (PackSegPtr psp, SeqAlignPtr salp)
{
  Uint1Ptr     seqpresence=NULL;
  Int2         numseg, aligndim, i, j; 
  SeqIdPtr     sip=NULL, siptemp;
  Int4Ptr      stptr=NULL, lenptr=NULL; 
  BioseqPtr    bsp=NULL;
  Int4         bslen, seg_start;

  if(!psp)
    ValMessage (salp, Err_Null_Segs, SEV_ERROR, NULL, NULL, 0);
  else
    {
      numseg=psp->numseg;
      aligndim=psp->dim;           
      sip=psp->ids;
      stptr=psp->starts;
      lenptr=psp->lens;
    
      if(stptr&&lenptr)
    {
      if(psp->present)
        {
          BSSeek(psp->present, 0, SEEK_SET);
          seqpresence=MemNew(BSLen(psp->present));
          if(!seqpresence)
        {
          
          ErrPostEx (SEV_ERROR, 0,0,  "Warning:insufficient memory");
          return;
          
        }
          BSRead(psp->present, seqpresence, BSLen(psp->present));
          /*go through each sequence*/
          for(j=0, siptemp=sip; j<aligndim && siptemp != NULL; siptemp=siptemp->next, j++)
        {  
          bsp=AlignValBioseqLockById(siptemp);
          if(bsp)
            {
              bslen=bsp->length; 
              AlignValBioseqUnlock (bsp);
              seg_start=stptr[j];
              /*check start*/
              if(seg_start<0)
            ValMessage (salp, Err_Start_Less_Than_Zero, SEV_ERROR, siptemp, sip , 0);     
              if(seg_start>=bslen)
            ValMessage (salp, Err_Start_More_Than_Biolen, SEV_ERROR, siptemp, sip , 0);
              
              /*go through each segment*/
              for(i=0; i<numseg; i++)
            {
              /*if this segment is present*/
              if(seqpresence[(i*aligndim+j)/8]&jybitnum[(i*aligndim+j)%8])      
                {
                  /*check start plus seg length*/
                  seg_start=seg_start+lenptr[i];
                  if(seg_start>bslen)
                 ValMessage (salp, Err_Sum_Len_Start, SEV_ERROR, siptemp, sip, numseg);
                }
            }
            }
        }
        }
    }
    }
  MemFree(seqpresence);         
}

/******************************************************************
check segment length, start and end point in Denseseg, Densediag and Stdseg
******************************************************************/ 
static void  ValidateSeqlengthinSeqAlign (SeqAlignPtr salp)
{
   
  if (salp)
  { 
      if(salp->segtype==1)
    ValidateSeqlengthInDenseDiag ((DenseDiagPtr)salp->segs, salp);
      else if(salp->segtype==2)
    ValidateSeqlengthInDenseSeg ((DenseSegPtr)salp->segs, salp);
      else if(salp->segtype==3)
    ValidateSeqlengthInStdSeg ((StdSegPtr)salp->segs, salp);
      else if(salp->segtype==4)
    ValidateSeqlengthInPackSeg ((PackSegPtr)salp->segs, salp);
  }
}

/******************************************************************
check if # of seqid matches the dimensions, and 
if there is only one seqeuence in seqalign
******************************************************************/ 
static void ValidateDimSeqIds (SeqAlignPtr salp)
{
  SeqIdPtr sip=NULL;
  DenseDiagPtr ddp=NULL, ddptemp;
  StdSegPtr ssp=NULL, ssptemp;
  DenseSegPtr dsp=NULL;
  Int4 numseg=0;
  
 if(salp)
   {
     /*densediag */
     if(salp->segtype==1)
       {
     
     ddp=(DenseDiagPtr)salp->segs;
     if(!ddp)
       ValMessage (salp, Err_Null_Segs, SEV_ERROR, NULL, NULL, 0);
     else
       for(ddptemp=ddp, numseg=0; ddptemp!=NULL; ddptemp=ddptemp->next, numseg++)
         {
           sip=ddptemp->id;
           if(ddptemp->dim==1)
         ValMessage (salp, Err_Segs_Dim_One, SEV_ERROR, NULL, sip , numseg+1);
           if(ddptemp->dim!=CountSeqIdInSip(sip))          
         ValMessage (salp, Err_Segs_DimSeqId_Not_Match, SEV_ERROR, NULL, sip , numseg+1);
          
         }
       }
     
     /*denseseg, packseg */
     else if(salp->segtype==2||salp->segtype==4)
       {
     dsp=(DenseSegPtr) (salp->segs);
     if(!dsp)
       ValMessage (salp, Err_Null_Segs, SEV_ERROR, NULL, NULL, 0);
     else
       {
         sip=dsp->ids;
         if(dsp->dim==1)
           ValMessage (salp, Err_SeqAlign_Dim_One, SEV_ERROR, NULL, sip , 0); 
         if(dsp->dim!=CountSeqIdInSip(sip)) 
            ValMessage (salp, Err_SeqAlign_DimSeqId_Not_Match, SEV_ERROR, NULL, sip , 0); 
           
       }
       }
     
     /*stdseg */
     else if(salp->segtype==3)
       {
     
     ssp=(StdSegPtr)salp->segs;
     if(!ssp)
       ValMessage (salp, Err_Null_Segs, SEV_ERROR, NULL, NULL, 0);
     else
       for(ssptemp=ssp, numseg=0; ssptemp!=NULL; ssptemp=ssptemp->next, numseg++)
         {
           
           sip=SeqIdInAlignSegs((Pointer)ssptemp, 3, salp);
           if(ssptemp->dim==1)
         ValMessage (salp, Err_Segs_Dim_One, SEV_ERROR, NULL, sip , numseg+1);
           if(ssptemp->dim!=CountSeqIdInSip( sip)) 
         ValMessage (salp, Err_Segs_DimSeqId_Not_Match, SEV_ERROR, NULL, sip , numseg+1);
           /*free Seqid if sip is a new chain created by SeqIdinAlignSegs*/
           
           SeqIdSetFree(sip);
         }
       }
   }
}

/******************************************************************
return true if a sip is contained in a seg, or false if otherwise 
Note it returns FASLE for an empty seqloc.  
It also returns false if error in sip or ssp
******************************************************************/ 
static Boolean IsSipContainedInStdseg(SeqIdPtr sip, StdSegPtr ssp)
{
  SeqLocPtr slp, slptemp;
  
  if(!sip||!ssp)
    return FALSE;

  slp=ssp->loc;
  for(slptemp=slp; slptemp!=NULL; slptemp=slptemp->next)
    {
      if(slptemp->choice!=SEQLOC_EMPTY&&SeqIdCmp(sip, SeqLocId(slptemp)))
    return TRUE;
    }
  
  return FALSE;
}

static Int4 PercentStringMatch (CharPtr string1, CharPtr string2)
{
  Int4 len1, len2, min_len, k, max_len;
  Int4 num_match = 0;
  
  if (StringHasNoText (string1) || StringHasNoText (string2))
  {
      return 0;
  }
  len1 = StringLen (string1);
  len2 = StringLen (string2);
  
  if (len1 > len2)
  {
      min_len = len2;
      max_len = len1;
  }
  else
  {
      min_len = len1;
      max_len = len2;
  }
  
  for (k = 0; k < min_len; k++)
  {
      if (string1[k] == string2[k] || string1[k] == 'N' || string2[k] == 'N')
      {
        num_match++;
      }
  }
  return (100 * num_match) / min_len;
}

static Boolean CheckForPercentMatch (SeqIdPtr sip_list)
{
  SeqIdPtr  sip_temp, sip_next;
  BioseqPtr bsp;
  CharPtr   master_seq = NULL, this_seq = NULL;  
  
  if (sip_list == NULL) return FALSE;
  sip_next = sip_list->next;
  sip_list->next = NULL;
  bsp = BioseqFind (sip_list);
  if (bsp != NULL)
  {
    master_seq = GetSequenceByBsp (bsp);      
  }
  sip_list->next = sip_next;
  sip_temp = sip_next;
  if (bsp == NULL || master_seq == NULL) 
  {
      return FALSE;
  }
  
  for (sip_temp = sip_next; sip_temp != NULL; sip_temp = sip_next)
  {
      sip_next = sip_temp->next;
      sip_temp->next = NULL;
      
      bsp = BioseqFind (sip_temp);
      if (bsp != NULL)
      {
        this_seq = GetSequenceByBsp (bsp);
      } else {
        this_seq = NULL;
      }
      
      sip_temp->next = sip_next;
      if (bsp == NULL || StringHasNoText (this_seq) || PercentStringMatch (master_seq, this_seq) < 50)
      {
        MemFree (this_seq);
        return FALSE;
      }
      MemFree (this_seq);
  }
  return TRUE;
}


/******************************************************************
check if an alignment is FASTA-like.  
If all gaps are at the 3' ends with dimensions>2, it's FASTA-like
******************************************************************/ 
static Boolean Is_Fasta_Seqalign (SeqAlignPtr salp)
{

  SeqIdPtr    siptemp=NULL;
  DenseSegPtr dsp;
  Int4Ptr     startp;
  Boolean     gap;
  Int4        k;
  Int2        j;
  SeqIdPtr    bad_sip = NULL;
  
  /*check only global or partial type*/
  if(salp->type!=1&&salp->type!=3)
    return FALSE;

  if (salp->segtype != SAS_DENSEG) {
    ValMessage (salp, Err_Unexpected_Alignment_Type, SEV_ERROR, NULL, NULL, 0);
  } else {
    dsp = (DenseSegPtr) salp->segs;
    if(!dsp)
    {
      ValMessage (salp, Err_Null_Segs, SEV_ERROR, NULL, NULL, 0);
    }
    else
    {
      if(dsp->dim<=2)
      {
        return FALSE;
      }
      /* if any sequence has gaps at the 5' end or internal gaps, the entire
       * alignment is declared to be valid.
       * if the sequence contains no gaps at all or only 3' end gaps, check
       * sequences for matches against the first sequence - if more than half
       * of the nucleotides are matches, then call this not FASTA-like.
       */ 
      for (j=0, siptemp=dsp->ids; j<dsp->dim&&siptemp; j++, siptemp=siptemp->next)
      {
        gap=FALSE;
          
        for (k=0; k<dsp->numseg; k++)
        {
          startp=dsp->starts;
          
          /*if start value is -1, set gap flag to true*/
          if (startp[dsp->dim*k + j] < 0)
          {
            gap = TRUE;              
          }
          /*if a positive start value is found after the initial -1 start value, then it's not  fasta like, no need to check this sequence further */
          else if(gap)
          {
            if (bad_sip != NULL)
            {
              SeqIdFree (bad_sip);
            }
            return FALSE;              
          }
          /* if no positive start value is found after the initial -1 start value
           * (indicating that gaps exist only at the 5' end) or if no gaps
           * were found at all, flag this sequence as bad if it is the first found.
           */
          if(k==dsp->numseg-1)
          {
            if (bad_sip == NULL)
            {
              bad_sip = SeqIdDup (siptemp);
            }
          }
        }
      }
      if (bad_sip != NULL)
      {
        if (! CheckForPercentMatch (dsp->ids))
        {
          ValMessage (salp, Err_Fastalike, SEV_WARNING, bad_sip, dsp->ids, 0);
          SeqIdFree (bad_sip);
          return TRUE;        
        }
        SeqIdFree (bad_sip);
        return FALSE;
      }
    }
  }
  /*no fasta like sequence is found*/
  return FALSE;
  
}  
  
 

/******************************************************************
check if there is a gap for all sequence in a segment
******************************************************************/ 
static void Segment_Gap_In_SeqAlign(SeqAlignPtr salp)
{
  Int4Ptr      stptr=NULL;
  DenseSegPtr  dsp=NULL;
  DenseDiagPtr ddp=NULL, ddptemp;
  StdSegPtr    ssp=NULL, ssptemp;
  PackSegPtr   psp=NULL;
  Uint1Ptr     seqpresence=NULL;
  Int2         numseg, aligndim, i, j; 
  SeqIdPtr     sip=NULL;
  SeqLocPtr    slp=NULL, slptemp;
  

  if(salp)
    {
      /*densediag*/
      if(salp->segtype==1)
    {
      ddp=(DenseDiagPtr)salp->segs;
      if(!ddp)
        ValMessage (salp, Err_Null_Segs, SEV_ERROR, NULL, NULL, 0);
      else
        {
          for(ddptemp=ddp, numseg=0; ddptemp!=NULL; ddptemp=ddptemp->next, numseg++)
        {
          sip=ddptemp->id;
          /*empty segment*/
          if(ddptemp->dim==0)   
            ValMessage (salp, Err_Segment_Gap, SEV_ERROR, NULL, sip, numseg);
        }
        }
    }
 
   
      /*denseseg*/
     else if(salp->segtype==2)
    {
      dsp=(DenseSegPtr)salp->segs;
      if(!dsp)
        ValMessage (salp, Err_Null_Segs, SEV_ERROR, NULL, NULL, 0);
      else
        {
          numseg=dsp->numseg;
          aligndim=dsp->dim;           
          stptr=dsp->starts;
          sip=dsp->ids;
          
          if(stptr==NULL)
        return;
          
          /*go through each segment*/
          for(j=0; j<numseg; j++)
        {    
          /*go through each sequence */
          for(i=0; i<aligndim; i++)
            {
              
              if(stptr[j*aligndim+i]==-1)
            {  
              /*all starts are -1 in this segment*/
              if(i==aligndim-1)
                ValMessage (salp, Err_Segment_Gap, SEV_ERROR, NULL, sip, j);
            }
              /*at least one start that is not -1*/
              else
            break;
              
            }
        }
        }
    }

        /*stdseg*/
     else if(salp->segtype==3)
    {
      ssp=(StdSegPtr)salp->segs;
      if(!ssp)
        ValMessage (salp, Err_Null_Segs, SEV_ERROR, NULL, NULL, 0);
      else
        {
          /*go through each segment*/
          for(ssptemp=ssp, numseg=0; ssptemp!=NULL; ssptemp=ssptemp->next, numseg++)
        {
          sip=SeqIdInAlignSegs((Pointer)ssptemp, 3, salp);
          slp=ssptemp->loc;
          /*go through each sequence*/
          for(slptemp=slp; slptemp!=NULL; slptemp=slptemp->next)
            { 
              if(slptemp->choice==SEQLOC_EMPTY||slptemp->choice==SEQLOC_NULL) 
            {
              if(slptemp->next)   
                continue;
              /*all seqloc are empty*/ 
              else
                ValMessage (salp, Err_Segment_Gap, SEV_ERROR, NULL, sip, numseg);
            }
              /*at least one non-empty seqloc*/
              else
            break;
            }
          /*free Seqid if sip is a new chain created by SeqIdinAlignSegs*/
          SeqIdSetFree(sip);
 
        }
        }
    }
      /*packseg*/
      else if(salp->segtype==4)
    {
      psp=(PackSegPtr)salp->segs;
      if(!psp)
        ValMessage (salp, Err_Null_Segs, SEV_ERROR, NULL, NULL, 0);
      else
        {
          numseg=psp->numseg;
          aligndim=psp->dim;           
          sip=psp->ids;
          if(psp->present)
        {
          BSSeek(psp->present, 0, SEEK_SET);
          seqpresence=MemNew(BSLen(psp->present));
          if(!seqpresence)
            {
              ErrPostEx (SEV_ERROR, 0,0,  "Warning:insufficient memory");
              return;
              
            }
          BSRead(psp->present, seqpresence, BSLen(psp->present));
          
          /*go through each segment*/
          for(j=0; j<numseg; j++)
            {    
              /*go through each sequence */
              for(i=0; i<aligndim; i++)
            {
              /*check the presence of each sequence by determining the bit value in a byte (0, not present; otherwise present)*/
              if(!(seqpresence[(j*aligndim+i)/8]&jybitnum[(j*aligndim+i)%8]))      
                {
                  /*more sequence to go*/
                  if(i<aligndim-1)   
                continue;
                  /*no sequence is present in this segment*/
                  else if(i==aligndim-1)
                ValMessage (salp, Err_Segment_Gap, SEV_ERROR, NULL, sip, j);
                }
              /*at least one sequence is present*/
              else
                break;
            }
            }
        MemFree(seqpresence);
        }
        }
    
    }


    }
}


static Boolean IsAlignmentTPA (SeqAlignPtr salp)
{
  Boolean isTPA = FALSE;
  BioseqPtr bsp;
  SeqIdPtr  sip = NULL, tmp_sip;
  SeqEntryPtr oldscope;
  DenseDiagPtr ddp;
  StdSegPtr    ssp;

  if (salp == NULL) {
    return FALSE;
  }

  oldscope = SeqEntrySetScope (NULL);

  switch (salp->segtype) {
    case SAS_DENDIAG:
      /*densediag */
      for (ddp = (DenseDiagPtr) salp->segs; ddp != NULL && !isTPA; ddp = ddp->next) {
        for (sip = SeqIdInAlignSegs((Pointer)ddp, salp->segtype, salp);
             sip != NULL && !isTPA;
             sip = sip->next) {
          bsp = BioseqLockById(sip);
          isTPA = HasTpaUserObject(bsp);
          BioseqUnlock(bsp);
        }
      }
      break;
    case SAS_STD: 
      /*Stdseg*/
      for (ssp = (StdSegPtr) salp->segs; ssp != NULL && !isTPA; ssp = ssp->next) {
        sip = SeqIdInAlignSegs((Pointer)ssp, salp->segtype, salp);
        for (tmp_sip = sip;
             tmp_sip != NULL && !isTPA;
             tmp_sip = tmp_sip->next) {
          bsp = BioseqLockById(tmp_sip);
          isTPA = HasTpaUserObject(bsp);
          BioseqUnlock(bsp);
        }
      }
      /*free Seqid if sip is a new chain created by SeqIdinAlignSegs*/
      SeqIdSetFree(sip);
      break;
    case SAS_DENSEG:
    case SAS_PACKED:
      /*Denseseg, Packseg*/
      for (sip=SeqIdInAlignSegs(salp->segs, salp->segtype, salp);
           sip != NULL && !isTPA;
           sip = sip->next) {
        bsp = BioseqLockById(sip);
        isTPA = HasTpaUserObject(bsp);
        BioseqUnlock(bsp);
      }
      break;
  }

  SeqEntrySetScope (oldscope);  
  return isTPA;
}


static void CheckAlnSeqLens (SeqAlignPtr salp)
{
  Int4     aln_len, start, stop;
  Int4     num_rows, row;
  SeqIdPtr sip;
  BioseqPtr bsp;
  Boolean   is_shorter = FALSE;

  if (salp == NULL) return;

  aln_len =  AlnMgr2GetAlnLength(salp, FALSE);
  num_rows = AlnMgr2GetNumRows(salp);
  if (num_rows < 0) {
    return;
  }

  for (row = 1; row <= num_rows && !is_shorter; row++) {
    sip = AlnMgr2GetNthSeqIdPtr(salp, row);
    bsp = BioseqFind (sip);
    if (bsp != NULL && bsp->idx.entityID == salp->idx.entityID) {
      AlnMgr2GetNthSeqRangeInSA(salp, row, &start, &stop);
      if ((stop > start && stop < bsp->length - 1) || (start > stop && start > bsp->length - 1)) {
        is_shorter = TRUE;
      }
    }
    sip = SeqIdFree (sip);
  }
  if (is_shorter) {
    ValMessage (salp, Err_Short_Aln, SEV_INFO, NULL, NULL, 0);
  }
}

 
/******************************************************************
validate seqid, segment length, strand in Seqalignment for Denseseg, 
Densediag and Stdseg.  Also check if it's FASTA-like
******************************************************************/ 
static Boolean ValidateSeqAlignFunc (SeqAlignPtr salp, Boolean find_remote_bsp)
{
  Boolean   error=FALSE;
  Uint2     pcnt_identity;
  SeqAlignPtr salp_test;
  
  if(salp==NULL)
    return FALSE;

  /*validate if dimesion equals number of seqid*/     
  ValidateDimSeqIds (salp);
        
  if (find_remote_bsp) {
    ValidateSeqIdInSeqAlign (salp);
    ValidateSeqlengthinSeqAlign (salp);
  }
  /*validate strand*/
  ValidateStrandinSeqAlign (salp);
       
  /*validate Fasta like*/
  if (Is_Fasta_Seqalign (salp))
  {
      error = TRUE;
  }
      
  /*validate segment gap*/
  Segment_Gap_In_SeqAlign (salp);
  
  if (!IsAlignmentTPA(salp)) {
    if (salp->segtype == SAS_DENDIAG) {
      /* duplicate alignment, to prevent indexing from changing the original type */
      salp_test = SeqAlignDup (salp);
      pcnt_identity = AlignmentPercentIdentityEx (salp_test, FALSE, TRUE);
      salp_test = SeqAlignFree (salp_test);
    } else {
      pcnt_identity = AlignmentPercentIdentityEx (salp, FALSE, TRUE);
    }

    if (pcnt_identity < 50) {
      ValMessage (salp, Err_Pcnt_ID, SEV_WARNING, NULL, NULL, pcnt_identity);
    }

/*    CheckAlnSeqLens (salp); */
  }
  
  return error;
}


/******************************************************************
validate each alignment sequentially.  
This function will subject the seqalign to all validation functions
******************************************************************/ 
NLM_EXTERN Boolean ValidateSeqAlign (SeqAlignPtr salp, Uint2 entityID, Boolean message,
                         Boolean msg_success, Boolean find_remote_bsp,
                         Boolean delete_bsp, Boolean delete_salp, BoolPtr dirty)
{  
  SeqAlignPtr  pre,
               salptmp;
  SaVal        sv;
  SaValPtr     svp;
  ValNodePtr   vnp;
  JYErrorMsgPtr bemp;
  MsgAnswer    ans;
  Int2         err_count=0,
               salp_count=0;
  Boolean      retdel = FALSE; 

  if(salp!=NULL)
  {
     sv.message = message;
     sv.msg_success = msg_success;
     sv.find_remote_bsp = find_remote_bsp;
     sv.delete_salp = delete_salp;
     sv.delete_bsp = delete_bsp;
     sv.retdel = TRUE;
     sv.do_hist_assembly = FALSE;
     sv.ids = NULL;
     sv.entityID = entityID; 
     sv.dirty = FALSE;   
     svp = &sv;   
     pre=NULL;
     salptmp=salp; 
     while (salptmp)
     {
        salp_count++;
        if (salptmp->segtype == SAS_SPARSE) {
           ValMessage (salp, Err_Segtype, SEV_WARNING, NULL, NULL, salptmp->segtype);
        } else if (salptmp->segtype == SAS_SPLICED) {
           ValMessage (salp, Err_Segtype, SEV_WARNING, NULL, NULL, salptmp->segtype);
        }
        else if (salptmp->segtype==5)
        {
           ValidateSeqAlign ((SeqAlignPtr) (salptmp->segs), entityID, message, msg_success, find_remote_bsp, delete_bsp, delete_salp, &svp->dirty);
        } 
        else if (salptmp->segtype<1 || salptmp->segtype>4)
        {
           ValMessage (salp, Err_Segtype, SEV_ERROR, NULL, NULL, salptmp->segtype);
        }
        else {
              ValidateSeqAlignFunc (salptmp, svp->find_remote_bsp);
        }         
           if (errorp)
           {
              if(svp->message)
              {
                 for (vnp=errorp; vnp!=NULL; vnp=vnp->next)
                 {
                    bemp=(JYErrorMsgPtr)vnp->data.ptrvalue;
                    ErrPostEx ((ErrSev) bemp->level, 0, 0, bemp->msg);
                 }
              }
              errorp = JYErrorChainDestroy (errorp);
              if (svp->delete_salp)
              {
            if (pre==NULL) {
              salp=salptmp->next;
              salptmp->next = NULL;
              SeqAlignFree (salptmp);
              salptmp = salp;
            }
            else {
              pre->next = salptmp->next;
              salptmp->next = NULL;
              SeqAlignFree (salptmp);
              salptmp = pre->next;
            }
           }
              else {
                 salptmp = salptmp->next;
              }
              err_count++;
           svp->retdel=FALSE;
        }
           else {
              salptmp = salptmp->next;
           }
     }
     if (err_count==0 && svp->msg_success) {
        if (salp_count>1)
           ans = Message (MSG_OK, "Validation test of %d alignments succeeded", salp_count);
        else
           ans = Message (MSG_OK, "Validation test of the alignment succeeded");
     }
     if (dirty)
        *dirty = svp->dirty;
     retdel = svp->retdel;
  }   
  return retdel;
} 


/******************************************************************
call back function for REGISTER_ALIGNVALIDATION defined in sequin4.c.  
Starting point for seqalign validation if user clicked on 
SeqalignValidation under menu Filer/Alignment.  
Either individual alignment or alignment block 
should be highlighted for this validation to work
******************************************************************/ 

NLM_EXTERN Int2 LIBCALLBACK ValidateSeqAlignFromData (Pointer data)
{ 
 
  OMProcControlPtr  ompcp;
  SeqAlignPtr       salp=NULL;
  SeqAnnotPtr       sap=NULL;
  SeqEntryPtr       sep=NULL;
  
  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  
  switch(ompcp->input_itemtype)
    {
    case OBJ_BIOSEQ :
      sep = SeqMgrGetSeqEntryForData (ompcp->input_data);
      break;
    case OBJ_BIOSEQSET :
      sep = SeqMgrGetSeqEntryForData (ompcp->input_data);
      break;
      /*if clicked on alignment block*/
    case OBJ_SEQANNOT:
      sap=(SeqAnnotPtr) (ompcp->input_data);
      break;
      /*if clicked on individual alignment*/
    case OBJ_SEQALIGN:
      salp=(SeqAlignPtr) (ompcp->input_data);
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  
  ErrSetMessageLevel(SEV_ERROR);
  if(sap!=NULL)
  {
     salp=is_salp_in_sap(sap, 2);
     ValidateSeqAlign (salp, 0, TRUE, TRUE, TRUE, FALSE, FALSE, NULL);
  }
  if (salp!=NULL) {
     ValidateSeqAlign (salp, 0, TRUE, TRUE, TRUE, FALSE, FALSE, NULL);
  }
  if (sep!=NULL) {
     ValidateSeqAlignInSeqEntry (sep, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE);
  }
  return OM_MSG_RET_DONE;
}

static void ValidateSeqAlignInAnnot (SeqAnnotPtr sap, SaValPtr svp)

{
  SeqAlignPtr  salp;

  while (sap != NULL) {
    if (sap->type == 2) {
      salp = (SeqAlignPtr) sap->data;
      if (salp != NULL) {
        ValidateSeqAlign (salp, svp->entityID, svp->message, svp->msg_success, svp->find_remote_bsp, svp->delete_bsp, svp->delete_salp, &svp->dirty);
      }
    }
    sap = sap->next;
  }
}

static void ValidateSeqAlignInHist (SeqHistPtr hist, SaValPtr svp)

{
  SeqAlignPtr  salp;

  if (hist == NULL) return;
  salp = hist->assembly;
  /* ValidateSeqAlign will validate the entire chain */
  ValidateSeqAlign (salp, svp->entityID, svp->message, svp->msg_success, svp->find_remote_bsp, svp->delete_bsp, svp->delete_salp, &svp->dirty);
}

static void ValidateSeqAlignCallback (SeqEntryPtr sep, Pointer mydata,
                                          Int4 index, Int2 indent)
{
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  SaValPtr           svp;

  if (sep != NULL && sep->data.ptrvalue && mydata != NULL) {
     svp = (SaValPtr)mydata;
     if (IS_Bioseq(sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp!=NULL) {
           ValidateSeqAlignInAnnot (bsp->annot, svp);
           if (svp != NULL && svp->do_hist_assembly) {
              ValidateSeqAlignInHist (bsp->hist, svp);
           }
        }
     }   
     else if(IS_Bioseq_set(sep)) {
        bssp = (BioseqSetPtr)sep->data.ptrvalue;
        if (bssp!=NULL) {
           ValidateSeqAlignInAnnot (bssp->annot, svp);
        }
     }
  }
}



NLM_EXTERN Boolean ValidateSeqAlignInSeqEntry (SeqEntryPtr sep, Boolean message, 
                                 Boolean msg_success, Boolean find_remote_bsp, 
                                 Boolean delete_bsp, Boolean delete_salp,
                                 Boolean do_hist_assembly)
{
  SeqEntryPtr      sep_head;
  Uint2            entityID;
  SaVal            sv;
  Boolean          success=TRUE;

  entityID = ObjMgrGetEntityIDForChoice (sep);
  if (entityID > 0) {
     sep_head = GetTopSeqEntryForEntityID (entityID);
     if (sep_head != NULL) {
        sv.message = message;
        sv.msg_success = msg_success;
        sv.find_remote_bsp = find_remote_bsp;
        sv.find_acc_bsp = FALSE;
        sv.delete_salp = delete_salp;
        sv.delete_bsp = delete_bsp;
        sv.retdel = TRUE;
        sv.do_hist_assembly = do_hist_assembly;
        sv.ids = NULL;
        sv.entityID = entityID; 
        sv.dirty = FALSE;
        SeqEntryExplore (sep_head, (Pointer)&sv, ValidateSeqAlignCallback);
        if (sv.dirty) {
           ObjMgrSetDirtyFlag (entityID, TRUE);
           ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
        }
        success = sv.retdel;
     }
  }
  return success;
}


/* alignment validator private for regular validator */

NLM_EXTERN Boolean ValidateSeqAlignWithinValidator (ValidStructPtr vsp, SeqEntryPtr sep, Boolean find_remote_bsp, Boolean do_hist_assembly);

NLM_EXTERN Boolean ValidateSeqAlignWithinValidator (ValidStructPtr vsp, SeqEntryPtr sep, Boolean find_remote_bsp, Boolean do_hist_assembly)

{
  GatherContext  gc;
  Boolean        rsult;

  if (vsp == NULL || sep == NULL) return FALSE;
  useLockByID = vsp->farIDsInAlignments;
  useValErr = TRUE;
  useVsp = vsp;
  vsp->gcp = &gc;
  vsp->bssp = NULL;
  vsp->bsp = NULL;
  vsp->sfp = NULL;
  vsp->descr = NULL;
  MemSet ((Pointer) &gc, 0, sizeof (GatherContext));
  rsult = ValidateSeqAlignInSeqEntry (sep, FALSE, FALSE, find_remote_bsp, FALSE, FALSE, do_hist_assembly);
  useLockByID = TRUE;
  useValErr = FALSE;
  useVsp = NULL;
  return rsult;
}


/* PopulateSample and ReadFromAlignmentSample are utility functions for AlignmentPercentIdentity */
static void PopulateSample (Uint1Ptr seqbuf_list, Int4Ptr start_list, 
                            Int4 sample_len, BioseqPtr PNTR bsp_list,
                            Int4 row)
{
  Char ch;
  
  if (seqbuf_list == NULL || start_list == NULL || sample_len < 1 || row < 0 || bsp_list == NULL
      || bsp_list[row] == NULL || start_list[row] < 0 || start_list[row] >= bsp_list[row]->length) {
      return;
  }

  ch = *(seqbuf_list + (row + 1) * sample_len);

  SeqPortStreamInt (bsp_list[row],
                    start_list[row], 
                    MIN (start_list[row] + sample_len - 1, bsp_list[row]->length - 1), 
                    Seq_strand_plus,
                    STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL,
                    seqbuf_list + row * sample_len,
                    NULL);

  /* put back char overwritten by SeqPortStreamInt */
  *(seqbuf_list + (row + 1) * sample_len) = ch;

}


static Uint1 ComplementChar (Uint1 ch)
{
  if (ch == 'A') {
    return 'T';
  } else if (ch == 'T') {
    return 'A';
  } else if (ch == 'G') {
    return 'C';
  } else if (ch == 'C') {
    return 'G';
  } else {
    return ch;
  }
}

static Uint1 ReadFromAlignmentSample(Uint1Ptr seqbuf_list, Int4Ptr start_list, 
                                     Int4 sample_len, BioseqPtr PNTR bsp_list,
                                     Uint1Ptr strand_list,
                                     Int4 row, Int4 seq_pos)
{
  Uint1 ch = 0;
  
  if (seqbuf_list == NULL || start_list == NULL || sample_len < 1 || row < 0 || bsp_list == NULL
      || bsp_list[row] == NULL || seq_pos < 0 || seq_pos >= bsp_list[row]->length) {
      return 0;
  }
  
  if (seq_pos < start_list[row] || seq_pos >= start_list[row] + sample_len) {
    start_list[row] = (seq_pos / sample_len) * sample_len;
    PopulateSample (seqbuf_list, start_list, 
                    sample_len, bsp_list,
                    row);
  }
  ch = seqbuf_list[(row * sample_len) + seq_pos - start_list[row]];
  if (strand_list[row] == Seq_strand_minus) {
    ch = ComplementChar(ch);
  }
  return ch;
}

typedef struct ambchar {
  Char ambig_char;
  CharPtr match_list;
} AmbCharData, PNTR AmbCharPtr;

static const AmbCharData ambiguity_list[] = {
 { 'R', "AG" },
 { 'Y', "CT" },
 { 'M', "AC" },
 { 'K', "GT" },
 { 'S', "CG" },
 { 'W', "AT" },
 { 'H', "ACT" },
 { 'B', "CGT" },
 { 'V', "ACG" },
 { 'D', "AGT" }};

static const Int4 num_ambiguities = sizeof (ambiguity_list) / sizeof (AmbCharData);

static Char AmbiguousMatch (Char ch1, Char ch2)
{
  Int4 i;
  for (i = 0; i < num_ambiguities; i++) {
    if (ch1 == ambiguity_list[i].ambig_char
        && StringChr (ambiguity_list[i].match_list, ch2)) {
      return ch2;
    } else if (ch2 == ambiguity_list[i].ambig_char
        && StringChr (ambiguity_list[i].match_list, ch1)) {
      return ch1;
    }
  }
  return 0;
}


extern double *
GetAlignmentColumnPercentIdentities 
(SeqAlignPtr salp,
 Int4    start,
 Int4    stop,
 Boolean internal_gaps,
 Boolean internal_validation)
{
  Int4       aln_len, num_rows, row, col_count = 0;
  Int4       num_match;
  Int4       aln_pos, seq_pos, k;
  Uint1          row_ch;
  SeqEntryPtr    oldscope;
  SeqIdPtr PNTR  sip_list;
  BioseqPtr PNTR bsp_list;
  Uint1Ptr       strand_list;
  BoolPtr        start_gap, end_gap;
  Int4Ptr        start_list;
  Uint1Ptr       seqbuf_list;
  Int4           sample_len = 50;
  Int4           chars_appearing[5]; /* 0 is A, 1 is T, 2 is G, 3 is C, 4 is internal gap */
  Int4           max_app, total_app, i;
  double *       pct_ids;
  
  if (salp == NULL || start < 0 || stop < start) return NULL;
 
  AlnMgr2IndexSingleChildSeqAlign(salp);
  aln_len = AlnMgr2GetAlnLength(salp, FALSE);
  num_rows = AlnMgr2GetNumRows(salp);
  if (num_rows < 0) {
    Message (MSG_POSTERR, "AlnMgr2GetNumRows failed");
    return NULL;
  }

  pct_ids = (double *) MemNew (sizeof (double) * (stop - start + 1));
  MemSet (pct_ids, 0, sizeof (double) * (stop - start + 1));

  bsp_list = (BioseqPtr PNTR) MemNew (num_rows * sizeof (BioseqPtr));
  sip_list = (SeqIdPtr PNTR) MemNew (num_rows * sizeof(SeqIdPtr));
  strand_list = (Uint1Ptr) MemNew (num_rows * sizeof(Uint1));
  start_gap = (BoolPtr) MemNew (num_rows * sizeof(Boolean));
  end_gap = (BoolPtr) MemNew (num_rows * sizeof(Boolean));
  for (row = 1; row <= num_rows; row++) {
    sip_list[row - 1] = AlnMgr2GetNthSeqIdPtr(salp, row);
    strand_list[row - 1] = AlnMgr2GetNthStrand(salp, row);
    bsp_list[row - 1] = BioseqLockById(sip_list[row - 1]);
    if (bsp_list[row - 1] == NULL) {
      oldscope = SeqEntrySetScope (NULL);
      bsp_list[row - 1] = BioseqLockById(sip_list[row - 1]);
      SeqEntrySetScope(oldscope);
      if (bsp_list[row - 1] == NULL) {
        break;
      }
    }
    start_gap[row - 1] = TRUE;
    end_gap[row - 1] = FALSE;
  }
  
  if (row <= num_rows) {
    Message (MSG_POSTERR, "Unable to locate Bioseq in alignment");
    while (row >= 0) {
      sip_list[row] = SeqIdFree(sip_list[row]);
      BioseqUnlock(bsp_list[row]);
      row--;
    }
    sip_list = MemFree (sip_list);
    bsp_list = MemFree (bsp_list);
    start_gap = MemFree (start_gap);
    end_gap = MemFree (end_gap);
    return 0;
  }
  
  start_list = (Int4Ptr) MemNew (num_rows * sizeof(Int4));
  seqbuf_list = (Uint1Ptr) MemNew (num_rows * sample_len * sizeof(Uint1));
  for (row = 0; row < num_rows; row++) {
    start_list[row] = 0;
    PopulateSample (seqbuf_list, start_list, 
                    sample_len, bsp_list,
                    row);
  }
  
  num_match = 0;
  for (aln_pos = start; aln_pos < aln_len && aln_pos <= stop; aln_pos++) {
    /* init lists */
    MemSet (chars_appearing, 0, sizeof (chars_appearing));
    for (row = 1; row <= num_rows; row++) {
      if (end_gap[row - 1]) {
        continue;
      }
      seq_pos = AlnMgr2MapSeqAlignToBioseq(salp, aln_pos, row);
      if (seq_pos < 0) {
        if (start_gap[row - 1] || end_gap[row - 1]) {
          /* beginning/end gap - never counts against percent identity */
        } else {
          k = aln_pos + 1;
          while (k < aln_len && seq_pos < 0) {
            seq_pos = AlnMgr2MapSeqAlignToBioseq(salp, k, row);
            k++;
          }
          if (seq_pos < 0) {
            /* now in end_gap for this sequence */
            end_gap[row - 1] = TRUE;
          } else {
            /* internal gaps count against percent identity when specified */
            if (internal_gaps) {
              chars_appearing[4] ++;
            }
          }
        }
      } else {
        start_gap[row - 1] = FALSE;
        
        row_ch = ReadFromAlignmentSample(seqbuf_list, start_list, 
                                         sample_len, bsp_list, strand_list,
                                         row - 1, seq_pos);
        switch (row_ch) {
          case 'A':
            chars_appearing[0]++;
            break;
          case 'T':
            chars_appearing[1]++;
            break;
          case 'G':
            chars_appearing[2]++;
            break;
          case 'C':
            chars_appearing[3]++;
            break;
          default:
            /* we don't count ambiguity characters */
            break;
       }
      }
    }
    max_app = 0;
    total_app = 0;
    for (i = 0; i < 4; i++) {
      if (chars_appearing[i] > max_app) {
        max_app = chars_appearing[i];
      }
      total_app += chars_appearing[i];
    }
    /* add in internal gaps */
    total_app += chars_appearing[4];
    if (total_app > 0) {
      pct_ids[aln_pos - start] = (double) max_app / (double) total_app;
    }
    col_count++;
  }
  
  for (row = 0; row < num_rows; row++) {
    sip_list[row] = SeqIdFree(sip_list[row]);
    BioseqUnlock(bsp_list[row]);
  }
  sip_list = MemFree (sip_list);
  bsp_list = MemFree (bsp_list);
  start_gap = MemFree (start_gap);
  end_gap = MemFree (end_gap);
  start_list = MemFree (start_list);
  seqbuf_list = MemFree (seqbuf_list);
      
  return pct_ids;  
}


static Uint2 AlignmentPercentIdentityEx (SeqAlignPtr salp, Boolean internal_gaps, Boolean internal_validation)
{
  Int4       aln_len, num_rows, row, col_count = 0;
  Int4       num_match;
  Uint2      pcnt;
  Boolean    row_match;
  Int4       aln_pos, seq_pos, tmp;
  Uint1          seq_ch, row_ch, amb_match;
  SeqEntryPtr    oldscope;
  SeqIdPtr PNTR  sip_list;
  BioseqPtr PNTR bsp_list;
  Uint1Ptr       strand_list;
  Int4Ptr        start_list;
  Uint1Ptr       seqbuf_list;
  Int4           sample_len = 50;
  Int4Ptr        starts, stops;
  Int4           match_25 = 0;
  
  if (salp == NULL) return 0;
 
  AlnMgr2IndexSingleChildSeqAlign(salp);
  aln_len = AlnMgr2GetAlnLength(salp, FALSE);
  num_rows = AlnMgr2GetNumRows(salp);
  if (num_rows < 0) {
    if (! internal_validation) {
      Message (MSG_POSTERR, "AlnMgr2GetNumRows failed");
    }
    return 0;
  }
  bsp_list = (BioseqPtr PNTR) MemNew (num_rows * sizeof (BioseqPtr));
  sip_list = (SeqIdPtr PNTR) MemNew (num_rows * sizeof(SeqIdPtr));
  strand_list = (Uint1Ptr) MemNew (num_rows * sizeof(Uint1));
  starts = (Int4Ptr) MemNew (num_rows * sizeof (Int4));
  stops = (Int4Ptr) MemNew (num_rows * sizeof (Int4));
  for (row = 1; row <= num_rows; row++) {
    sip_list[row - 1] = AlnMgr2GetNthSeqIdPtr(salp, row);
    strand_list[row - 1] = AlnMgr2GetNthStrand(salp, row);
    bsp_list[row - 1] = BioseqLockById(sip_list[row - 1]);
    if (bsp_list[row - 1] == NULL) {
      oldscope = SeqEntrySetScope (NULL);
      bsp_list[row - 1] = BioseqLockById(sip_list[row - 1]);
      SeqEntrySetScope(oldscope);
      if (bsp_list[row - 1] == NULL) {
        break;
      }
    }
    /* get endpoints for each row */
    AlnMgr2GetNthSeqRangeInSA(salp, row, starts + row - 1, stops + row - 1);
    starts[row - 1] = AlnMgr2MapBioseqToSeqAlign (salp, starts[row - 1], row);
    stops[row - 1] = AlnMgr2MapBioseqToSeqAlign (salp, stops[row - 1], row);
    if (starts[row - 1] > stops[row - 1]) {
      tmp = starts[row - 1];
      starts[row - 1] = stops[row - 1];
      stops[row - 1] = tmp;
    }

  }
  
  if (row <= num_rows) {
    if (! internal_validation) {
      Message (MSG_POSTERR, "Unable to locate Bioseq in alignment");
    }
    while (row > 0) {
      sip_list[row - 1] = SeqIdFree(sip_list[row - 1]);
      BioseqUnlock(bsp_list[row - 1]);
      row--;
    }
    sip_list = MemFree (sip_list);
    bsp_list = MemFree (bsp_list);
    starts = MemFree (starts);
    stops = MemFree (stops);
    return 0;
  }
  
  start_list = (Int4Ptr) MemNew (num_rows * sizeof(Int4));
  seqbuf_list = (Uint1Ptr) MemNew ((num_rows * sample_len + 1) * sizeof(Uint1));
  for (row = 0; row < num_rows; row++) {
    start_list[row] = 0;
    PopulateSample (seqbuf_list, start_list, 
                    sample_len, bsp_list,
                    row);
  }
  
  num_match = 0;
  for (aln_pos = 0; aln_pos < aln_len; aln_pos++) {
    row_match = TRUE;
    seq_ch = 0;
    for (row = 1; row <= num_rows; row++) {
      if (aln_pos < starts[row - 1] || aln_pos > stops[row - 1]) {
        continue;
      }
      seq_pos = AlnMgr2MapSeqAlignToBioseq(salp, aln_pos, row);
      if (seq_pos < 0) {
        if (internal_gaps) {
          row_match = FALSE;
        }
      } else {        
        row_ch = ReadFromAlignmentSample(seqbuf_list, start_list, 
                                         sample_len, bsp_list, strand_list,
                                         row - 1, seq_pos);
        if (row_ch == 'N') {
          /* do nothing - Ns do not count against percent identity */
        } else if (seq_ch == 0) {
          seq_ch = row_ch;
        } else if (seq_ch != row_ch) {
          amb_match = AmbiguousMatch (seq_ch, row_ch);
          if (amb_match == 0) {
            row_match = FALSE;
          } else {
            seq_ch = amb_match;
          }
        }
      }
      if (!row_match) {
        break;
      }
    }
    if (row_match) {
      num_match++;
      match_25++;
    }
    col_count++;
    if (col_count % 25 == 0) {
      match_25 = 0;
    }
  }
  
  for (row = 0; row < num_rows; row++) {
    sip_list[row] = SeqIdFree(sip_list[row]);
    BioseqUnlock(bsp_list[row]);
  }
  sip_list = MemFree (sip_list);
  bsp_list = MemFree (bsp_list);
  starts = MemFree (starts);
  stops = MemFree (stops);
  start_list = MemFree (start_list);
  seqbuf_list = MemFree (seqbuf_list);
      
  if (col_count == 0) {
      pcnt = 0;
  } else {
      pcnt = (100 * num_match) / col_count;
  }
  return pcnt;
}

extern Uint2 AlignmentPercentIdentity (SeqAlignPtr salp, Boolean internal_gaps)
{
  return AlignmentPercentIdentityEx (salp, internal_gaps, FALSE);
}

extern Uint2 WeightedAlignmentPercentIdentity (SeqAlignPtr salp, Boolean internal_gaps)
{
  Int4       aln_len, num_rows, row, col_count = 0;
  Int4       num_match;
  Uint2      pcnt;
  Int4       aln_pos, seq_pos, k;
  Uint1          row_ch;
  SeqEntryPtr    oldscope;
  SeqIdPtr PNTR  sip_list;
  BioseqPtr PNTR bsp_list;
  Uint1Ptr       strand_list;
  BoolPtr        start_gap, end_gap;
  Int4Ptr        start_list;
  Uint1Ptr       seqbuf_list;
  Int4           sample_len = 50;
  Int4           chars_appearing[5]; /* 0 is A, 1 is T, 2 is G, 3 is C, 4 is internal gap */
  double         col_pct, col_pct_total = 0;
  Int4           max_app, total_app, i;
  
  if (salp == NULL) return 0;
 
  AlnMgr2IndexSingleChildSeqAlign(salp);
  aln_len = AlnMgr2GetAlnLength(salp, FALSE);
  num_rows = AlnMgr2GetNumRows(salp);
  if (num_rows < 0) {
    Message (MSG_POSTERR, "AlnMgr2GetNumRows failed");
    return 0;
  }
  bsp_list = (BioseqPtr PNTR) MemNew (num_rows * sizeof (BioseqPtr));
  sip_list = (SeqIdPtr PNTR) MemNew (num_rows * sizeof(SeqIdPtr));
  strand_list = (Uint1Ptr) MemNew (num_rows * sizeof(Uint1));
  start_gap = (BoolPtr) MemNew (num_rows * sizeof(Boolean));
  end_gap = (BoolPtr) MemNew (num_rows * sizeof(Boolean));
  for (row = 1; row <= num_rows; row++) {
    sip_list[row - 1] = AlnMgr2GetNthSeqIdPtr(salp, row);
    strand_list[row - 1] = AlnMgr2GetNthStrand(salp, row);
    bsp_list[row - 1] = BioseqLockById(sip_list[row - 1]);
    if (bsp_list[row - 1] == NULL) {
      oldscope = SeqEntrySetScope (NULL);
      bsp_list[row - 1] = BioseqLockById(sip_list[row - 1]);
      SeqEntrySetScope(oldscope);
      if (bsp_list[row - 1] == NULL) {
        break;
      }
    }
    start_gap[row - 1] = TRUE;
    end_gap[row - 1] = FALSE;
  }
  
  if (row <= num_rows) {
    Message (MSG_POSTERR, "Unable to locate Bioseq in alignment");
    while (row >= 0) {
      sip_list[row] = SeqIdFree(sip_list[row]);
      BioseqUnlock(bsp_list[row]);
      row--;
    }
    sip_list = MemFree (sip_list);
    bsp_list = MemFree (bsp_list);
    start_gap = MemFree (start_gap);
    end_gap = MemFree (end_gap);
    return 0;
  }
  
  start_list = (Int4Ptr) MemNew (num_rows * sizeof(Int4));
  seqbuf_list = (Uint1Ptr) MemNew (num_rows * sample_len * sizeof(Uint1));
  for (row = 0; row < num_rows; row++) {
    start_list[row] = 0;
    PopulateSample (seqbuf_list, start_list, 
                    sample_len, bsp_list,
                    row);
  }
  
  num_match = 0;
  for (aln_pos = 0; aln_pos < aln_len; aln_pos++) {
    /* init lists */
    MemSet (chars_appearing, 0, sizeof (chars_appearing));
    for (row = 1; row <= num_rows; row++) {
      if (end_gap[row - 1]) {
        continue;
      }
      seq_pos = AlnMgr2MapSeqAlignToBioseq(salp, aln_pos, row);
      if (seq_pos < 0) {
        if (start_gap[row - 1] || end_gap[row - 1]) {
          /* beginning/end gap - never counts against percent identity */
        } else {
          k = aln_pos + 1;
          while (k < aln_len && seq_pos < 0) {
            seq_pos = AlnMgr2MapSeqAlignToBioseq(salp, k, row);
            k++;
          }
          if (seq_pos < 0) {
            /* now in end_gap for this sequence */
            end_gap[row - 1] = TRUE;
          } else {
            /* internal gaps count against percent identity when specified */
            if (internal_gaps) {
              chars_appearing[4] ++;
            }
          }
        }
      } else {
        start_gap[row - 1] = FALSE;
        
        row_ch = ReadFromAlignmentSample(seqbuf_list, start_list, 
                                         sample_len, bsp_list, strand_list,
                                         row - 1, seq_pos);
        switch (row_ch) {
          case 'A':
            chars_appearing[0]++;
            break;
          case 'T':
            chars_appearing[1]++;
            break;
          case 'G':
            chars_appearing[2]++;
            break;
          case 'C':
            chars_appearing[3]++;
            break;
          default:
            /* we don't count ambiguity characters */
            break;
       }
      }
    }
    max_app = 0;
    total_app = 0;
    for (i = 0; i < 4; i++) {
      if (chars_appearing[i] > max_app) {
        max_app = chars_appearing[i];
      }
      total_app += chars_appearing[i];
    }
    if (total_app > 0) {
      col_pct = (double) max_app / (double) total_app;
      col_pct_total += col_pct;
    }
    col_count++;
  }
  
  for (row = 0; row < num_rows; row++) {
    sip_list[row] = SeqIdFree(sip_list[row]);
    BioseqUnlock(bsp_list[row]);
  }
  sip_list = MemFree (sip_list);
  bsp_list = MemFree (bsp_list);
  start_gap = MemFree (start_gap);
  end_gap = MemFree (end_gap);
  start_list = MemFree (start_list);
  seqbuf_list = MemFree (seqbuf_list);
      
  if (col_count == 0) {
      pcnt = 0;
  } else {
      pcnt = (100 * col_pct_total) / col_count;
  }
  return pcnt;
}

