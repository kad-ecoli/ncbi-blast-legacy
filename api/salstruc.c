/* ===========================================================================
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
* File Name:  salstruc.c
*
* Author:  Colombe Chappey
*
* Version Creation Date:   1/27/96
*
* ==========================================================================
*/
#include <salstruc.h>
#include <salsa.h>
#include <salutil.h>
#include <salsap.h>
#include <subutil.h>
#include <gather.h>
#include <satutil.h>
#include <sqnutils.h>
#include <alignmgr2.h>


NLM_EXTERN SelStructPtr BufferFree (SelStructPtr ssp)
{
  SelStructPtr   ssptmp, ssptmp2, next;
  SelEdStructPtr sep;
  ValNodePtr     vnp;

  for (ssptmp = ssp; ssptmp != NULL; ssptmp = ssptmp->next)
            ssptmp->region = NULL;
  ssptmp = ssp; 
  while (ssptmp != NULL)
  {
         next = ssptmp->next;
       ssptmp2 = ssptmp;
         if (ssptmp2->region != NULL)
         {
            sep = (SelEdStructPtr) ssptmp2->region;
            if (sep->region != NULL) 
               sep->region = SeqLocFree ((SeqLocPtr) sep->region);
            if (sep->data != NULL) {
               vnp = (ValNodePtr) sep->data;
               vnp->data.ptrvalue = NULL;
               sep->data = ValNodeFree (vnp);
            }
            MemFree (sep);
            ssptmp2->region = NULL;
         }
         MemFree (ssptmp2);
         ssptmp2 = NULL;
         ssptmp = next;
  }
  return NULL;
}

/*******************************************************************
***  
***    SetupDataBuffer
***    minbufferlength = (Int4) ( WINPERBUF * MAXCharLine * j );
*******************************************************************/
NLM_EXTERN EditAlignDataPtr SetupDataBuffer (EditAlignDataPtr adp)
{
  if (adp==NULL) 
     return NULL;
  if ( adp->seqnumber == 0 ) return NULL;
  adp->minbufferlength = TMP_BUFFERLENGTH;
  adp->bufferlength = (Int4) MIN (adp->length, adp->minbufferlength);
  return adp;
}

/*******************************************************************
***  
***    SetupDataPanel
***
*******************************************************************/
NLM_EXTERN EditAlignDataPtr SetupDataPanel (EditAlignDataPtr adp)
{
  Int4  j;
  Int4  lg;

  if (adp==NULL) 
     return NULL;
  if ( adp->seqnumber == 0 ) return NULL;
  lg = MAXLineWindow + 4;

  if ( adp->item_id != NULL ) MemFree (adp->item_id);
  adp->item_id = NULL;
  adp->item_id = (Uint4Ptr)MemNew ((size_t) (lg * sizeof(Uint4)));
  for (j=0; j<MAXLineWindow; j++) adp->item_id[j] = 0;

  if ( adp->seqEntity_id != NULL ) MemFree (adp->seqEntity_id);
  adp->seqEntity_id = NULL;
  adp->seqEntity_id =(Uint2Ptr)MemNew((size_t) (lg * sizeof(Uint2)) );
  for (j=0; j<MAXLineWindow; j++) adp->seqEntity_id[j] = 0;

  if ( adp->itemtype != NULL ) MemFree (adp->itemtype);
  adp->itemtype = NULL;
  adp->itemtype =(Uint2Ptr)MemNew ((size_t) (lg * sizeof(Uint2)));
  for (j=0; j<MAXLineWindow; j++) adp->itemtype[j] = 0;

  if ( adp->itemsubtype != NULL ) MemFree (adp->itemsubtype);
  adp->itemsubtype = NULL;
  adp->itemsubtype =(Uint2Ptr)MemNew((size_t) (lg * sizeof(Uint2)) );
  for (j=0; j<MAXLineWindow; j++) adp->itemsubtype[j] = 0;

  if ( adp->alignline != NULL ) MemFree (adp->alignline);
  adp->alignline = NULL;
  adp->alignline= (Uint4Ptr)MemNew((size_t) (lg * sizeof(Uint4)));
  for (j=0; j<MAXLineWindow; j++) adp->alignline[j] = 0;

  if ( adp->colonne != NULL ) MemFree (adp->colonne);
  adp->colonne = NULL;
  lg = adp->minbufferlength + adp->editbuffer + 4;
  adp->colonne = (Int4Ptr) MemNew ((size_t) (lg * sizeof(Int4)));
  for (j=0; j<adp->minbufferlength +adp->editbuffer; j++) adp->colonne[j] = -1;

  if (adp->item_id==NULL || adp->seqEntity_id==NULL || adp->itemsubtype==NULL 
  || adp->alignline==NULL || adp->colonne==NULL)  {
         adp->seqnumber = 0;
         return NULL;
  }
  return adp;
}

static Uint2 OBJ_ (Uint2 feattype)
{
  if ( feattype == FEATDEF_BAD ) return OBJ_BIOSEQ;
  return OBJ_SEQFEAT;
}

/***********************************************************
***
************************************************************/
typedef struct ccid {
  SeqIdPtr    sip;
  SeqEntryPtr sep;
  BioseqPtr   bsp;
} CcId, PNTR CcIdPtr;
 
typedef struct orgscan {
  ObjMgrPtr  omp;
  Int2       nuclCode;
  Int2       mitoCode;
  Boolean    mito;
  Char       taxname [64];
} OrgScan, PNTR OrgScanPtr;

static Boolean CC_OrgScanGatherFunc (GatherContextPtr gcp)

{
  BioSourcePtr   biop;
  ObjMgrTypePtr  omtp;
  OrgNamePtr     onp;
  OrgRefPtr      orp;
  OrgScanPtr     osp;
  ValNodePtr     sdp;
  SeqFeatPtr     sfp;
  Uint2          subtype;
  Int2           val;
  ValNodePtr     vnp;

  if (gcp == NULL || gcp->thisitem == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT  && gcp->thistype != OBJ_SEQDESC) return TRUE;

  osp = (OrgScanPtr) gcp->userdata;
  if (osp == NULL) return TRUE;

  subtype = 0;   
  if (gcp->thistype == OBJ_SEQFEAT  || gcp->thistype == OBJ_SEQDESC) {
    omtp = ObjMgrTypeFind (osp->omp, gcp->thistype, NULL, NULL);
    if (omtp == NULL) {
      return TRUE;
    }
    if (omtp->subtypefunc != NULL) {
      subtype = (*(omtp->subtypefunc)) (gcp->thisitem);
    }
  }  

  orp = NULL;
  biop = NULL;
  switch (gcp->thistype) {
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) gcp->thisitem;
      switch (subtype) {
        case FEATDEF_ORG :
          orp = (OrgRefPtr) sfp->data.value.ptrvalue;
          break;
        case FEATDEF_BIOSRC :
          biop = (BioSourcePtr) sfp->data.value.ptrvalue;
          break;
        default :
          break;
      }
      break;
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) gcp->thisitem;
      switch (subtype) {
        case Seq_descr_modif :
          vnp = (ValNodePtr) sdp->data.ptrvalue;
          while (vnp != NULL) {
            val = (Int2) vnp->data.intvalue;
            if (val == MODIF_mitochondrial || val == MODIF_kinetoplast) {
              osp->mito = TRUE;
            }  
            vnp = vnp->next;
          }
          break;
        case Seq_descr_org :
          orp = (OrgRefPtr) sdp->data.ptrvalue;
          break;
        case Seq_descr_source :
          biop = (BioSourcePtr) sdp->data.ptrvalue;
          break;
        default :
          break;
      }
      break;
    default :
      break;
  }
 
  if (orp == NULL && biop != NULL) {
    orp = biop->org;
    osp->mito = (Boolean) (biop->genome == 4 || biop->genome == 5 || biop->genome == 20);
  }
  if (orp != NULL) {
    StringNCpy_0 (osp->taxname, orp->taxname, sizeof (osp->taxname));
    onp = orp->orgname;
    if (onp != NULL) {
      osp->nuclCode = onp->gcode;
      osp->mitoCode = onp->mgcode;
    }
  }
    
  return TRUE;
}

static Int2 CC_SeqEntryOrEntityIDToGeneticCode (SeqEntryPtr sep, Uint2 entityID, BoolPtr mito, CharPtr taxname, size_t maxsize)
{
  GatherScope  gs;
  OrgScan      osp;

  if (mito != NULL) {
    *mito = FALSE;
  }
  if (taxname != NULL && maxsize) {
    *taxname = '\0';
  }
  osp.mito = FALSE;
  osp.nuclCode = 0;
  osp.mitoCode = 0;
  osp.omp = ObjMgrGet ();
  osp.taxname [0] = '\0';
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQDESC] = FALSE;
  if (sep != NULL) {
    gs.scope = sep;
    GatherSeqEntry (sep, (Pointer) &osp, CC_OrgScanGatherFunc, &gs);
  } else if (entityID) {
    GatherEntity (entityID, (Pointer) &osp, CC_OrgScanGatherFunc, &gs);
  }
  if (mito != NULL) {
    *mito = osp.mito;
  }
  if (taxname != NULL && maxsize) {
    StringNCpy_0 (taxname, osp.taxname, maxsize);
  }
  if (osp.mito) {
    return osp.mitoCode;
  } else {
    return osp.nuclCode;
  }
}


static void FindSeqEntryForSeqIdCallback (SeqEntryPtr sep, Pointer mydata,
                                          Int4 index, Int2 indent)
{
  BioseqPtr          bsp;
  SeqIdPtr           sip;
  CcIdPtr            cip;
 
  if (sep != NULL && sep->data.ptrvalue && mydata != NULL) {
     cip = (CcIdPtr)mydata;
     if (cip->sep==NULL && IS_Bioseq(sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp!=NULL) {
           sip = SeqIdFindBest(bsp->id, 0);
           if (SeqIdForSameBioseq(cip->sip, sip)) {
              cip->sep = sep;
              cip->bsp = bsp;
           }
        }
     }   
  }
}
 
 
static Int2 CC_SeqEntryToGeneticCode (Uint2 entityID, SeqIdPtr sip)
{
  SeqEntryPtr sep_head,
              sep;
  CcId        ci;
  Int2        genCode = 0;

  sep_head  = GetTopSeqEntryForEntityID (entityID);
  ci.sip = SeqIdDup (sip);
  ci.sep = NULL;
  SeqEntryExplore(sep_head,(Pointer)&ci, FindSeqEntryForSeqIdCallback);
  sep = ci.sep;
  SeqIdFree (ci.sip);
  if (sep!=NULL) {
/*
     genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);*/
     genCode = CC_SeqEntryOrEntityIDToGeneticCode (sep, 0, NULL, NULL, 0);
  }
  return genCode;
}

/************************************************** 
***  next_notemptyline: search for not empty line (not_empty_string) 
***      if not exists, go to next alignment line  
***  get_substring
***
***************************************************/
static CharPtr get_substring (CharPtr str, Int4 drw_start, Int4 drw_width)
{
  Int4            width;
  Int4            stringlens;
  CharPtr         strp;
  if (str == NULL ) return NULL; 
  stringlens = StringLen (str);
  if ( drw_start >= stringlens ) { return NULL; } 
  strp = str + drw_start;
  stringlens = StringLen (strp);
  if (stringlens == 0) return NULL; 
  width = MIN ((Int2) drw_width, (Int2) stringlens);
  if ( !not_empty_string (strp, width) ) return NULL;
  return strp;
}

NLM_EXTERN CharPtr next_notemptyline (ValNodePtr anp_list, ValNodePtr linebuff, Int2 numberalignline, Int2 *index, Int4 start, Int4 *drw_width, TextAlignBufPtr *tdp, AlignNodePtr *anp)
{
  TextAlignBufPtr tdptmp;
  ValNodePtr      vnp;
  ValNodePtr      vnpanp;
  CharPtr         str;
  Int4            width = *drw_width;
  Int2            j = 1;

  if (*index > numberalignline) return NULL;
  vnp = linebuff;
  tdptmp = (TextAlignBufPtr) vnp->data.ptrvalue;
  vnpanp = anp_list;
  if (j < *index ) {
         vnp = vnp->next;
         while (j < *index) {
           tdptmp = (TextAlignBufPtr) vnp->data.ptrvalue;
           if (OBJ_(tdptmp->feattype) == OBJ_BIOSEQ) vnpanp = vnpanp->next;
           j++;
           if (j == *index ) break;
           vnp = vnp->next;
         }
  }
  while ( vnp != NULL && j <= numberalignline ) 
  {
         str = get_substring (tdptmp->buf, start, width);
         if ( str != NULL ) break;
         vnp = vnp->next;
         if (vnp == NULL) break;
         tdptmp = (TextAlignBufPtr) vnp->data.ptrvalue;
         if (OBJ_(tdptmp->feattype) == OBJ_BIOSEQ) vnpanp = vnpanp->next;
         j++;
  } 
  if ( j > numberalignline && vnp == NULL ) return NULL;
  *index = j;
  *drw_width = width;
  *tdp = tdptmp;
  if (OBJ_(tdptmp->feattype) == OBJ_BIOSEQ) 
     *anp= (AlignNodePtr) vnpanp->data.ptrvalue;
  else 
     *anp = NULL;
  return str;
}

/******************************************************
***  GetPerCol
***    .seqline    : index of the sequences  (0, 1..)
***    .alignline  : index of alignment line (0, 1..)  
***    
*******************************************************/
static void GetPerCol (EditAlignDataPtr adp, Int4 hoffset)
{
  SeqAlignPtr  salp;
  CompSegPtr   dsp;
  Int4Ptr      lenp;
  AlignNodePtr anp;
  AlignSegPtr  segs;
  Int4         j = 0; 
  Int4         m;
  Int4         k, l;

  if (adp->minbufferlength == 0) return;
  if (adp->sap_align != NULL && (salp = (SeqAlignPtr)adp->sap_align->data) != NULL) {
         j = 0;
         while ( j < adp->bufferlength && salp != NULL )
         {
                dsp = (CompSegPtr) salp->segs;
                lenp = dsp->lens;
                k = m = 0; 
                while (k < dsp->numseg && m < hoffset) 
                {
                    for (l = 0; l < *lenp && m < hoffset; l++, m++) continue;
                    if (m == hoffset) break;
                    lenp++;
                    k++;
                }
                for ( ; k <= dsp->numseg && j < adp->bufferlength; k++, lenp++)
                {
                   for (l=0; l < *lenp && j < adp->bufferlength; l++, m++, j++)
                   {
                      adp->colonne[j] = m;
                   }
                }
                if ( j < adp->bufferlength && salp->next != NULL ) 
                {
                   for (l=0; l < adp->intersalpwidth && j < adp->bufferlength; l++, j++)
                      adp->colonne[j] = -1;
                   salp = salp->next;
                }
                else break;
         }
  }
  else if ( adp->anp_list != NULL ) 
  {
         anp = (AlignNodePtr) adp->anp_list->data.ptrvalue;
         segs = anp->segs;
         j=0;
         for (k=0, m=0; segs != NULL && m < hoffset; k++, segs=segs->next)
            for (l=0; l<(segs->gr.right-segs->gr.left+1) && m<hoffset; l++, m++)
                continue;
         for (; segs != NULL && j < adp->bufferlength; k++, segs=segs->next)
            for(l=0; l<(segs->gr.right-segs->gr.left+1) && j< adp->bufferlength; l++, m++, j
++)
                adp->colonne[j] = m;
  }
  else {
         return;
  }
  if ( j < adp->minbufferlength + adp->editbuffer && adp->colonne[j-1] > -1) {
         adp->colonne[j] = adp->colonne[j-1]+1;
         j++;
  }
  if ( j < adp->minbufferlength + adp->editbuffer )
         for (; j < adp->minbufferlength + adp->editbuffer; j++)
                adp->colonne[j] = -1;
  return ;
}

/*********************************************************************
***  is_feature_to_buffer
*********************************************************************/
NLM_EXTERN SelEdStructPtr is_feature_to_buffer (ValNodePtr vnphead, Uint4 bspitemID, Uint2 entityID, Int4 from, Int4 drw_width, SeqAlignPtr salp, Uint2 seqedit, ValNodePtr sqloc_list)
{
  SelEdStructPtr  cds;
  SeqIdPtr        sip;
  SeqLocPtr       slp;
  Int4            start, stop;
  Int4            start2, stop2;
  Int2            chklocp;
  Boolean         draw = FALSE;

  if (vnphead == NULL) return NULL;
  cds = (SelEdStructPtr) vnphead->data.ptrvalue;
  if ( cds == NULL ) 
     return NULL;
  sip = SeqLocId ((SeqLocPtr) cds->region);
  for (; cds != NULL && !draw; cds = cds->next) 
  {
     if (entityID == cds->entityID && bspitemID == cds->bsp_itemID)
     {
        if ( cds->regiontype == OM_REGION_SEQLOC && cds->region != NULL ) 
        {
           slp = (SeqLocPtr) cds->region;
           start2 = SeqLocStart(slp);
           chklocp =chkloc(sip, start2, sqloc_list, &start2);
           start= SeqCoordToAlignCoord(start2, sip, salp, 0, chklocp);
           stop2 = SeqLocStop(slp);
           chklocp =chkloc(sip, stop2, sqloc_list, &stop2);
           stop = SeqCoordToAlignCoord(stop2, sip, salp, 0, chklocp);
           if (start <= stop && start < from + drw_width && stop >= from)  {
                 draw = TRUE;
                 break;
           }
        }
     }
  }
  if (draw) 
     return cds;
  return NULL;
}

/******************************************************************/
NLM_EXTERN ByteStorePtr cds_to_pept (SeqLocPtr slp, Uint1 frame, Int2 gencode, Boolean include_stop)
{
  ByteStorePtr  bs;
  ValNodePtr    code;
  CdRegionPtr   crp;
  SeqFeatPtr    sfp;
  ValNodePtr    vnp;

  if (slp == NULL) return NULL;
  sfp = SeqFeatNew ();
  if (sfp == NULL) return NULL;
  sfp->data.choice = SEQFEAT_CDREGION;
  crp = CdRegionNew ();
  sfp->data.value.ptrvalue = (Pointer) crp;
  if (crp == NULL) {
         SeqFeatFree (sfp);
         return NULL;
  }
  crp->orf = FALSE;
  crp->conflict = FALSE;
  crp->frame = frame;
  crp->gaps = 0;
  crp->mismatch = 0;
  crp->stops = 0;
  code = ValNodeNew (NULL);
  if (code != NULL) {
         code->choice = 254;
         vnp = ValNodeNew (NULL);
         code->data.ptrvalue = vnp;
         if (vnp != NULL) {
            vnp->choice = 2;
            vnp->data.intvalue = (Int4) gencode;
         }
  }
  crp->genetic_code = code;
  crp->code_break = NULL;
  sfp->location = slp;
  bs = ProteinFromCdRegion (sfp, include_stop);
  if (bs == NULL) {
     sfp->location = NULL;
     SeqFeatFree (sfp);
     return NULL;
  }
  return bs;
}

/******************************************************************/
static CharPtr SeqAlignTranslate (SeqAlignPtr salp, Uint2 entityID, Int4 from, Int4 to, Uint1 codonbase, SeqIdPtr sip, Uint1 strand, ValNodePtr sqlocs)
{
  SeqLocPtr    slp;
  SeqIntPtr    sit;
  ByteStorePtr bsp;
  CharPtr      pep = NULL, 
               pepPtr = NULL;
  CharPtr      str = NULL, 
               strPtr = NULL;
  CharPtr      buffer = NULL, 
               bufferPtr = NULL;
  Int4         k, 
               slplens,
               strlens;
  Int4         start, stop,
               len;
  Int2         cb;  
  Int2         genCode;
  BioseqPtr    bsq;

  if ( salp == NULL ) return NULL; 
  if ( salp->segtype != COMPSEG ) {
     return NULL; 
  }
  start= (Int4) AlignCoordToSeqCoord (from, sip, salp, sqlocs, 0);
  stop = (Int4) AlignCoordToSeqCoord (to, sip, salp, sqlocs, 0);

  bsq = BioseqLockById (sip);
  if (bsq != NULL) {
     if (start + stop >= bsq->length) 
        stop = bsq->length - 1;
     BioseqUnlock (bsq);
  }
  else return NULL;

  slp = fuzz_loc (start, stop, strand, sip, TRUE, TRUE);
  if (slp == NULL) {
     return NULL; 
  }
  slplens = SeqLocLen (slp);
  if ( slplens < (Int4) (3 + codonbase) ) {
     return NULL; 
  }
  genCode = CC_SeqEntryToGeneticCode (entityID, sip);
  if (genCode == 0)
     genCode = Seq_code_ncbieaa; 

  sit = (SeqIntPtr) slp->data.ptrvalue;
  if (strand == Seq_strand_plus) 
  {
     sit->from = sit->from + codonbase;
     slplens = SeqLocLen (slp);
     cb = (Int2)(slplens % (Int4) 3); 
     if (cb == 1) {
          sit->to --;
     } 
  }
  else if (strand == Seq_strand_minus) {
     sit->to = sit->to - codonbase;
     slplens = SeqLocLen (slp);
     cb = (Int2)(slplens % (Int4) 3); 
     if (cb == 1 && sit->from >0) {
          sit->from --;
     } else if (cb == 2) {
          sit->from ++;
     } 
     if (cb == 0) codonbase = 0;
     else if (cb == 1) codonbase = 1;
     else if (cb == 2) codonbase = 2;
  }
  slplens = SeqLocLen (slp);
  if ( slplens >= 3 ) {
     bsp = cds_to_pept (slp, 1, genCode, TRUE);
     str = (CharPtr) BSMerge (bsp, NULL);
     BSFree (bsp);
     pep = (CharPtr)MemNew ((size_t) ((slplens +5) *sizeof(Char)));
     pep = emptystring (pep, (Int4)(slplens + 5));
     pep [slplens +3] = '\0';
     pepPtr = pep;
     *pepPtr = ' ';
     pepPtr += codonbase +1;
     strlens = 3*StringLen(str);
     if (slplens < strlens) {
        strlens=(Int4)(slplens/(Int4)3);
        str [strlens] ='\0';
     }
     if (strand == Seq_strand_minus)
          reverse_string (str);
     strlens = StringLen(str);
     strPtr = str;
     for (k = 0; k < strlens; k++, pepPtr += 3, strPtr++) {
          *pepPtr = *strPtr; 
     }
     MemFree (str);
     buffer = (CharPtr)MemNew ((size_t) ((to -from +5) *sizeof(Char)));
     buffer = emptystring (buffer, (Int4)(to -from +5));
     buffer [to -from +3] = '\0';
     bufferPtr = buffer;
     *bufferPtr = ' ';
     buffer = ReadBufferFromSap (pep, buffer, salp, sip, from, to, &len);
     MemFree (pep);
  }
  SeqLocFree (slp);
  return buffer; 
}

static CharPtr prot_to_putprot (CharPtr trans)
{
  CharPtr strp;
  Boolean in = TRUE;

  for (strp = trans; *strp != '\0'; strp++)
  {
         if (*strp == 'M') in = TRUE;
         else if (*strp == '*') {
            if (!in) *strp = ' ';
            else in = FALSE;
         }
         else if (in) *strp = '~';
         else *strp = ' ';
  }
  return trans;
}

static CharPtr prot_to_rputprot (CharPtr trans)
{
  CharPtr strp, strptmp;
  Boolean in = TRUE;
  Int4    j, k;

  strptmp = strp = trans; 
  while (*strptmp != '\0') {
     for (j=0; *strptmp != '\0'; strptmp++, j++) {
         if (*strptmp == '*') {
            for (k=0; k<=j && *strp != '\0'; strp++, k++) {
               if (*strp != '*' && *strp != 'M' ) *strp = ' ';
            }
            in = TRUE; 
            break;
         }
         if ( *strptmp == 'M' ) {
            if (in) {
               for (k=0; k<=j && *strp != '\0'; strp++, k++) {
                  if (*strp != 'M' && *strp != '*') *strp = '~';
               }
               in = FALSE; 
            }
            else {
               for (k=0; k<j && *strp != '\0'; strp++, k++) {
                  if (*strp != 'M' && *strp != '*') *strp = ' ';
               }
            }
            break;
         }
     }
     if (strptmp == '\0') break;
     strptmp++;
     j++;
  }
  if (in) {
     for (k=0; k<j && *strp != '\0'; strp++, k++) 
         if (*strp != 'M' && *strp != '*') *strp = '~';
  }
  else {
     for (k=0; k<j && *strp != '\0'; strp++, k++) 
         if (*strp != 'M' && *strp != '*') *strp = ' ';
  }
  strptmp = strp = trans; 
/*
  strp++;
  for(; *strp != '\0'; strp++, strptmp++)
     if (*strptmp == '*' && *strp != '~' ) *strptmp = ' ';
*/
  return trans;
}


static SelStructPtr addssp (SelStructPtr *ssphead, Uint2 ei, Uint2 choice, Pointer pt, Uint2 iID)
{
  SelStructPtr ssp, hssp;

  ssp = (SelStructPtr) MemNew (sizeof (SelStruct));
  if (ssp == NULL) 
     return *ssphead;
  ssp->entityID = ei;
  ssp->itemID = iID;
  ssp->itemtype = choice;
  ssp->regiontype = OM_REGION_SEQLOC;  
  ssp->region = (Pointer) pt;
  ssp->next = NULL;
  ssp->prev = NULL;
  hssp = *ssphead; 
  if (hssp == NULL) {
     *ssphead = ssp;
  }
  else {
     for (; hssp->next != NULL; hssp = hssp->next) continue;
     hssp->next = ssp;
     ssp->prev = hssp;
  }
  return ssp;
}

static SelStructPtr MakeRf (SelStructPtr buffer, Uint2 entityID, Uint2 itemID, Uint2 feattype, SeqIdPtr bspsip, Uint1 strand, Uint1 j, SeqAlignPtr salp, Int4 fromseq, Int4 toseq, EditAlignDataPtr adp, Int4 from_bufferstart, Int4 to_bufferstart, Uint2 typerf)
{
  SelEdStructPtr  rf;
  CharPtr         trans;

  trans = (CharPtr) SeqAlignTranslate (salp, entityID, fromseq, toseq, (Uint1)(j), bspsip, strand, adp->sqloc_list);
  if (trans != NULL) {
     if ( adp->prot_mode == PUTPROT ) {
        if (strand == Seq_strand_minus) trans = prot_to_rputprot (trans);
        else trans = prot_to_putprot (trans); 
     }
     rf=new_seledstruct(entityID, itemID, feattype, 0,itemID, from_bufferstart, to_bufferstart,  bspsip, strand, TRUE, NULL, (Pointer)trans, 0, 1);
     addssp (&(buffer), entityID, typerf, (Pointer) rf, itemID);
  }
  return buffer;
}


static SelStructPtr get_firstline (SelEdStructPtr sesp1, SelStructPtr buffer)
{
  SelEdStructPtr sesp;
  SelStructPtr   buf;

  if (buffer == NULL) {
         return NULL;
  }
  if (sesp1 == NULL) {
         return buffer;
  }
  for (buf=buffer; buf!=NULL; buf=buf->next) 
  {
         sesp = (SelEdStructPtr) buf->region;
         if (is_sameId (sesp1->entityID, sesp1->itemID, sesp1->itemtype, 255, sesp->entityID, sesp->itemID, sesp->itemtype, 255) ) 
            break;
  }
  if (buf == NULL) {
         return buffer;
  }
  return buf;
}

static Boolean has_complement (ValNodePtr params, Uint2 entityID, Uint4 itemID)
{
  ValNodePtr  vnp;
  SeqParamPtr prm;

  for (vnp = params; vnp != NULL; vnp = vnp->next)
  {
     prm = (SeqParamPtr) vnp->data.ptrvalue;
     if (prm->entityID == entityID && prm->itemID == itemID) {
        if ( prm->complement ) return TRUE;
     }
  }
  return FALSE;
}

static Boolean rf_on (ValNodePtr params, Uint2 entityID, Uint4 itemID, Uint2 rf)
{
  ValNodePtr  vnp;
  SeqParamPtr prm;

  for (vnp = params; vnp != NULL; vnp = vnp->next)
  {
     prm = (SeqParamPtr) vnp->data.ptrvalue;
     if (prm->entityID == entityID && prm->itemID == itemID) {
        return prm->rf[rf];
     }
  }
  return FALSE;
}

/*********************************************************************
***  arrange_buffer
*********************************************************************/
static void arrange_buffer (EditAlignDataPtr adp)
{
  TextAlignBufPtr tdp;
  AlignNodePtr    anp;
  SelEdStructPtr  ssp = NULL, 
                  bspssp = NULL;
  SelEdStructPtr  firstsesp = NULL;
  SeqIdPtr        bspsip;
  CharPtr         bufstr = NULL;
  Int4            fromseq = adp->bufferstart;
  Int4            toseq = adp->bufferstart + adp->bufferlength -1;
  Int4            from = 0;
  Int4            to = adp->bufferlength-1;
  Int4            length = adp->bufferlength;
  Int2            index = 1;
  Uint4           itemID;
  Uint2           entityID;
  Int2            j, k;

  SeqAlignPtr     salp = (SeqAlignPtr) adp->sap_align->data;
  SelEdStructPtr  drawfeat = NULL;
  ValNodePtr      vnpfeat;

  if (adp->firstssp != NULL) {
     firstsesp = SelEdStructDup ((SelEdStructPtr) adp->firstssp->region);
  }
  BufferFree (adp->buffer);
  adp->buffer = NULL;
  if (adp->draw_scale) {
         ssp =new_seledstruct(LINE0,EDITDEF_SCA,EDITDEF_SCA,0,LINE0,from, to, NULL, 0, FALSE, NULL, NULL, 0, 1);
         addssp (&(adp->buffer), 0, EDITDEF_SCA, (Pointer) ssp, EDITDEF_SCA);
  }
  if (adp->draw_bars) {
         ssp =new_seledstruct(LINE0,EDITDEF_SCB,EDITDEF_SCB,0,LINE0,from, to, NULL, 0, FALSE, NULL, NULL, 0, 1);
         addssp (&(adp->buffer), 0, EDITDEF_SCB, (Pointer) ssp, EDITDEF_SCB);
  }
  bufstr =next_notemptyline (adp->anp_list, adp->linebuff, (Int2)adp->numberalignline, &index, from, &length, &tdp, &anp);
  if ( index > adp->numberalignline ) {
         return;
  }
  while ( index <= adp->numberalignline && bufstr != NULL) 
  {
         if ( OBJ_(tdp->feattype) == OBJ_BIOSEQ ) {
            itemID = anp->bsp_itemID;
            if ( adp->input_format == OBJ_BIOSEQ ) {
               entityID = tdp->seqEntityID;
            }
            else if ( adp->input_format == OBJ_SEQALIGN ) {
               entityID = anp->seq_entityID;
            }
            bspsip = anp->sip;
/**
WriteLog ("ARRANGE %d %d   %d %d %d %d  %d %d %d %d  %d  \n", entityID, itemID, anp->entityID, anp->itemID, anp->seq_entityID, anp->bsp_itemID, tdp->entityID, tdp->itemID, tdp->seqEntityID, tdp->bsp_itemID, OBJ_(tdp->feattype));
**/
         } 
         else {
            itemID = tdp->itemID;
            entityID = tdp->entityID;
/**
WriteLog ("ARRANGFEATE %d %d   %d %d %d %d  %d   %d  \n", entityID, itemID, tdp->entityID, tdp->itemID, tdp->seqEntityID, tdp->bsp_itemID, tdp->feattype, OBJ_(tdp->feattype));
**/
         }
         if (OBJ_(tdp->feattype) == OBJ_BIOSEQ && (is_seqvisible(entityID, itemID, adp->seq_info)) )
         {
            bspssp = new_seledstruct (entityID, itemID, OBJ_(tdp->feattype), 0,itemID, 
                             from + adp->bufferstart, to + adp->bufferstart,  
                             bspsip, Seq_strand_plus, TRUE, NULL, (Pointer) tdp, 0, 1);
            addssp (&(adp->buffer), entityID, FEATDEF_BAD, (Pointer) bspssp, itemID);

            if (has_complement (adp->params, bspssp->entityID, itemID)) {
                addssp (&(adp->buffer), entityID, EDITDEF_CPL, (Pointer) bspssp, itemID);
            }
            for (j=0; j<3;j++) {
               if (rf_on (adp->params, bspssp->entityID, itemID, j))  
                  MakeRf (adp->buffer, entityID, itemID, OBJ_(tdp->feattype), bspsip, Seq_strand_plus, (Uint1)j, salp, fromseq, toseq, adp, from + adp->bufferstart, to + adp->bufferstart, (Uint2)(EDITDEF_RF1 + j));
            }
            for (j=3, k=2; j<6;j++, k--) {
               if (rf_on (adp->params, bspssp->entityID, itemID, j))  
                  MakeRf (adp->buffer, entityID, itemID, OBJ_(tdp->feattype), bspsip, Seq_strand_minus, (Uint1)k, salp, fromseq, toseq, adp, from + adp->bufferstart, to + adp->bufferstart, (Uint2)(EDITDEF_RF4 + k));
            }
            if ( adp->seqfeat != NULL ) 
            {
                vnpfeat = adp->seqfeat;
                for (; vnpfeat != NULL; vnpfeat = vnpfeat->next) 
                {
                   drawfeat = is_feature_to_buffer (vnpfeat, itemID, entityID, fromseq, length, salp, adp->input_format, adp->sqloc_list);
                   if (drawfeat != NULL) 
                   {
                      ssp = drawfeat;
                      addssp (&(adp->buffer), entityID, (Uint1)vnpfeat->choice, (Pointer) ssp, ssp->entityID);
                      if ( vnpfeat->choice == FEATDEF_CDS ) 
                      {
                         if ( drawfeat->data != NULL ) {
                            ssp = drawfeat;
                            addssp (&(adp->buffer), entityID, FEATDEF_TRSL,(Pointer) ssp, ssp->entityID);
                         }
                         tdp = TextAlignBufFind (adp->linebuff, drawfeat->entityID, drawfeat->itemID, drawfeat->itemtype);
                         if ( tdp != NULL )
                         {
                            ssp = new_seledstruct(tdp->entityID, tdp->itemID,  OBJ_SEQFEAT, 0,itemID, 
                                                  from + adp->bufferstart, to + adp->bufferstart, bspsip, Seq_strand_plus, TRUE, NULL, (Pointer) tdp, 0, 1);
                            addssp(&(adp->buffer), entityID, FEATDEF_PROT, (Pointer) ssp, tdp->itemID);
                         }
                      }
                   }
                }
            }
            if ( adp->feat != NULL ) {
                vnpfeat = adp->feat;
                for (; vnpfeat != NULL; vnpfeat = vnpfeat->next)
                {
                   drawfeat = is_feature_to_buffer(vnpfeat, itemID, entityID, fromseq, length, salp, adp->input_format, adp->sqloc_list);
                   if (drawfeat != NULL) 
                   {
                      ssp = drawfeat;
                      addssp (&(adp->buffer), entityID, (Uint1)vnpfeat->choice, (Pointer) ssp, ssp->entityID);
                      if (vnpfeat->choice == SEQFEAT_CDREGION 
                      && drawfeat->data != NULL) {
                         ssp = drawfeat;
                         addssp (&(adp->buffer), entityID, FEATDEF_TRSL, (Pointer) ssp, ssp->entityID);
                      }
                   }
                }
            }
         }
         index++;
         length = adp->bufferlength;
         bufstr = next_notemptyline (adp->anp_list, adp->linebuff, 
                     (Int2)adp->numberalignline, &index, from, &length, &tdp, &anp);
  }
  adp->buffertail = adp->buffer;
  if (adp->buffertail !=NULL)
     for (; adp->buffertail->next !=NULL; adp->buffertail =adp->buffertail->next)
        continue;
  adp->firstssp = get_firstline (firstsesp, adp->buffer);
  SelEdStructDel (firstsesp);
  return;
}


static Boolean update_fromalignnode (EditAlignDataPtr adp)
{
  AlignNodePtr  anp;
  ValNodePtr    curr;      /*for the list of AlignNodePtr*/
  ValNodePtr    list;      /*list of DrawText*/
  ValNodePtr    tdp_list;  /*list of DrawText*/
  TextAlignBufPtr tdp;
  Int4          start;
  Int4          p_stop=0;
  Uint2         entityID;
  SeqIdPtr      sip;
  ValNodePtr    vnp;
  SeqParamPtr   prm;
  SelEdStructPtr seq_info;
  Boolean       first_draw;
  Int2 j;

  if (adp->minbufferlength == 0) return FALSE;
  if (adp->linebuff != NULL) 
         adp->linebuff = (ValNodePtr) FreeTextAlignList(adp->linebuff);
  adp->linebuff = NULL;
  adp->numberalignline = 0;
  if ( adp->anp_list == NULL ) {
     ErrPostEx (SEV_ERROR, 0, 0, "fail in update_fromalignnode [1]");
  }
{
/***!!!!!!!!!!!!!!**/
SelEdStructPtr tmp;
  for (tmp=adp->seq_info; tmp!=NULL; tmp=tmp->next)
  {
     if (tmp->entityID==0)
        break;
  }
  first_draw=(Boolean)(tmp!=NULL);
}
  vnp = adp->params;
  for (curr = adp->anp_list, j=0; curr != NULL; curr = curr->next, j++)
  {
         anp = (AlignNodePtr) curr->data.ptrvalue;
         if ( anp == NULL ) {
                break;
         }
         start = (Int4) (adp->gr.left +adp->bufferstart);
         list = (ValNodePtr) ProcessTextAlignNode (anp, start, 
                start + adp->bufferlength + adp->editbuffer, &p_stop, NULL, 
                adp->visibleWidth, (Int1)0, (Uint4)0, NULL);  
         if ( list == NULL ) {
                break;
         }
         tdp_list = list;

         while (tdp_list != NULL)
	 {
                tdp = (TextAlignBufPtr) tdp_list->data.ptrvalue;
                if (adp->input_format == OBJ_BIOSEQ) {
                   entityID = tdp->seqEntityID;
                } else if ( adp->input_format == OBJ_SEQALIGN ) {
                   entityID = anp->seq_entityID;
                }
                if (adp->master.entityID == 0)
                {
                   if (adp->master.region != NULL) {
                      sip = SeqLocId((SeqLocPtr)adp->master.region);
                      if (SeqIdForSameBioseq(sip, anp->sip)) 
                      {
                         adp->master.entityID = entityID;
                         adp->master.itemID = anp->bsp_itemID;
                         adp->master.itemtype = OBJ_BIOSEQ;
                         adp->caret.entityID = entityID;
                         adp->caret.itemID = anp->bsp_itemID;
                         adp->caret.itemtype = OBJ_BIOSEQ;
                      }
                   }
                }
                if (vnp != NULL) {
                   prm = (SeqParamPtr) vnp->data.ptrvalue;
                   if (prm->entityID == 0) {
                      prm->entityID = entityID;
                      prm->itemID = anp->bsp_itemID;
                   }
                }
                if(tdp->buf != NULL)
                {
                       if (tdp->label == NULL) {
                          tdp->label=(CharPtr)MemNew((size_t)(64*sizeof(Char)));
                          tdp->label = StringSave ("unknown");
                       }
                       adp->numberalignline++;
                       ValNodeAddPointer (&adp->linebuff, 0, (Pointer) tdp);
                }
                if (first_draw)
                {
                   for (seq_info=adp->seq_info;seq_info!=NULL;seq_info=seq_info->next)
                   {
                      sip=SeqLocId((SeqLocPtr)seq_info->region);
                      if (SeqIdForSameBioseq(sip, anp->sip))
                      {
                         seq_info->entityID=entityID;
                         seq_info->itemID=anp->bsp_itemID;
                         seq_info->itemtype=OBJ_BIOSEQ;
                         break;
                      }
                   }
                }
                tdp_list = tdp_list->next;
         }
         if (vnp != NULL) vnp = vnp->next;
  }
  if (adp->numberalignline == 0)
  { 
         return FALSE;
  } 
  GetPerCol (adp, adp->bufferstart);
  return TRUE;
}

static Int4 get_bufferstart (EditAlignDataPtr adp)
{
  Int4 start, modstart = 0;

  if ( adp->minbufferlength == 0 ) return 0;
  if ( adp->hoffset < adp->minbufferlength/3 ) { 
         start = 0;
  }
  else if ( adp->length < adp->minbufferlength ) {
         start = 0;
  }
  else if ( adp->hoffset > adp->length- adp->minbufferlength) {
         start = adp->length - adp->minbufferlength;
  }
  else {
         start = adp->hoffset-(adp->minbufferlength/3);
  }
  if (start > 0) {
         modstart = start % adp->visibleWidth; 
  }
  if (modstart > 0) {
         start -= modstart;
  }
  return start;
}

static void count_feature_buf_line (ValNodePtr fnp_list, Int4 g_left, Int4 g_right, ValNodePtr PNTR feature_line)
{
        FeatNodePtr fnp;
        Int4 c_left, c_right;
        ValNodePtr vnp;
        Boolean found;

        if(fnp_list == NULL)
                return;
        
        while(fnp_list)
        {
           fnp = (FeatNodePtr)fnp_list->data.ptrvalue;
           c_left = fnp->extremes.left;
           c_right = fnp->extremes.right;
           if(!(c_left > g_right || c_right < g_left))
           {
                found = FALSE;
                for(vnp = *feature_line; vnp != NULL && !found; vnp = vnp->next)
                {
                        if(vnp->data.intvalue == (Int4)(fnp->itemID))
                                found = TRUE;
                }
                if(!found)
                        ValNodeAddInt(feature_line, 0, (Int4)(fnp->itemID));
           }
           fnp_list = fnp_list->next;
        }
}

static Int4 CountTextAlignNodeNum(AlignNodePtr anp, Int4 m_left, Int4 m_right)
{
        Int4 num_line = 0;
        Int4 g_left, g_right;

        AlignSegPtr asp;
        ValNodePtr feature_line, curr;  /*the number of lines for a feature*/

        g_left = anp->extremes.left;
        g_right = anp->extremes.right;
        if(m_left > g_right || m_right < g_left)
                return 0;

        num_line = 1;
        feature_line = NULL;

        /*process  the GAPs and the DIAGs segs*/
        for(asp = anp->segs; asp !=NULL; asp = asp->next)
        {
           g_left = asp->gr.left;
           g_right = asp->gr.right;
           if(!(g_left > m_right || g_right < m_left))
           {
              switch(asp->type)
              {  
                case GAP_SEG:
                   break;

                case REG_SEG:
                case DIAG_SEG:
                   g_left = MAX(m_left, g_left);
                   g_right = MIN(m_right, g_right);
                   count_feature_buf_line (asp->cnp, g_left, g_right, &feature_line);
                   break;
                default:
                   break;
              }
           }
           if(g_left > m_right)
                break;
        }
        if(feature_line != NULL)
        {
           for(curr = feature_line; curr != NULL; curr = curr->next)
                ++num_line;
           ValNodeFree(feature_line);
        }
 
        return num_line;
}

static Int4 count_total_line(ValNodePtr anp_list, Int4 line_len, Int4 left, Int4 right, Int4Ptr h_blocks)
{
        AlignNodePtr anp;
        Int4 c_start, c_stop;
        Int4 line_num = 0;
        ValNodePtr curr;
         
        if(anp_list == NULL)
                return line_num;
        if(h_blocks != NULL)
                *h_blocks = 0;
        anp = (AlignNodePtr)anp_list->data.ptrvalue;
        if(left == -1)
                left = anp->extremes.left;
        if(right == -1)
                right = anp->extremes.right;
        if(left > anp->extremes.right || right < anp->extremes.left)
                return line_num;
        left = MAX(left, anp->extremes.left);
        right = MIN(right, anp->extremes.right);
        if (left >= right)
                return line_num;
        c_start = left;
        while(c_start <=right)
        {
           c_stop = MIN(right, (c_start+line_len-1));
           for(curr = anp_list; curr != NULL; curr = curr->next)
           {  
              anp = (AlignNodePtr)curr->data.ptrvalue;
              line_num += CountTextAlignNodeNum(anp, c_start, c_stop);
           }
           if(h_blocks != NULL)
              ++(*h_blocks);
           c_start = c_stop+1;
        }
  return line_num;
}
 
static Int4 addline_perblock (EditAlignDataPtr adp, Int4 diffs)
{
  Int4        line = 0;
  ValNodePtr  vnp;
  SeqParamPtr prm;
  Int1        j;

  if (adp->draw_scale) line += (Int4) diffs;
  if (adp->draw_bars)  line += (Int4) diffs;
  for (vnp = adp->params; vnp != NULL; vnp = vnp->next)
  {
     prm = (SeqParamPtr) vnp->data.ptrvalue;
     if ( prm->complement ) line += (Int4) diffs;
     for (j=0; j<=6; j++) 
        if (prm->rf[j]) line += (Int4) diffs;
  }
  return line;
}  
 
static Int2 feat_linenum (Int4 slp_start, Int4 slp_stop, Int4 line_len, Int4 left,
Int4 right)
{
  Int4         modstart;
  Int4         modstop;
 
  slp_start = MAX (slp_start, left);
  modstart = slp_start % line_len;
  if ( modstart > 0) slp_start -= modstart;
 
  slp_stop = MIN (slp_stop, right);
  modstop = slp_stop % line_len;
  if ( modstop > 0) slp_stop += line_len;
 
  return (Int2)((slp_stop - slp_start) / line_len);
}
 
static Int4 CountFeatNum (ValNodePtr adpfeat, Int4 line_len, Int4 left, Int4 right){
  ValNodePtr   vnp;
  SelEdStructPtr sesp;
  SeqLocPtr    slp;
  Int4         line = 0;
 
  for (vnp = adpfeat; vnp != NULL; vnp = vnp->next)
  {
     sesp = (SelEdStructPtr) vnp->data.ptrvalue;
     for (; sesp != NULL; sesp = sesp->next) 
     {
        if (vnp->choice ==FEATDEF_CDS && sesp->regiontype ==OM_REGION_SEQLOC 
        && sesp->region !=NULL)
        {
           slp = (SeqLocPtr) sesp->region;;
           if (SeqLocStart(slp) > right || SeqLocStop(slp) < left) 
              line+=0; 
           else {
              line+= feat_linenum (SeqLocStart(slp), SeqLocStop(slp), line_len, left, right);
           }
        }
     }
  }
  return line;
}
 
static Int4 get_tot_line (EditAlignDataPtr adp, Int4 line_len, Int4 left, Int4 right)
{
  Int4         diffs, line = 0;

  if (left >= right )           
      return line;
  line=count_total_line (adp->anp_list, line_len, left, right, NULL);
  diffs = right - left;
  if ( (diffs % line_len) > 0 ) diffs += line_len;
  diffs = (Int4) (diffs / line_len);
  line += (Int4) addline_perblock (adp, diffs);
  line += (Int4) CountFeatNum (adp->feat, line_len, left, right);
  line += (Int4) CountFeatNum (adp->seqfeat, line_len, left, right);
  return line;
}

NLM_EXTERN void data_collect_arrange (EditAlignDataPtr adp, Boolean recollect)
{
  ValNodePtr         vnp = NULL;
  Int4               maxscroll;
  Int2               x;
  Boolean            goOn = FALSE;

  x = adp->seqnumber;
  if (adp->draw_scale) x++;
  if (adp->draw_bars) x++;
  maxscroll = adp->length * x / adp->visibleWidth;
  
  if (adp->int4value2 != -1) /* Feature line count is disabled when scrolling */
  {
  	adp->nlines = get_tot_line (adp, adp->visibleWidth, 0, adp->length-1);
  	adp->nlines = adp->nlines - adp->pnlLine +3;
  	adp->nlines = MAX ((Int4)0, adp->nlines);
  }

  if (recollect) 
  {
     adp->bufferstart = get_bufferstart (adp);
     adp->bufferlength= MIN ((Int4) (adp->length - adp->bufferstart), 
                             (Int4) adp->minbufferlength);
     if (abs((adp->bufferstart +adp->bufferlength ) -adp->length) < 100) 
     {
        adp->bufferlength = adp->length - adp->bufferstart;
     }
     goOn = update_fromalignnode (adp);

     if (adp->Cn3D_On)
{{
  TextAlignBufPtr tdp, 
                  tdp1;
  CharPtr         tmp;
  DenseDiagPtr    ddp;
  Int4            start, len, j;
  Int4            start1, stop1;
  SeqIdPtr        sip,
                  seqsip;
  SeqAlignPtr     salp = (SeqAlignPtr)(adp->sap_align->data);
  Int2        chklocp;

  if(adp->seqnumber > 1) {
     for (vnp=adp->linebuff; vnp!=NULL; vnp=vnp->next) {
        tdp = (TextAlignBufPtr) vnp->data.ptrvalue;
        for (tmp=tdp->buf; *tmp!='\0'; ++tmp) {
           *tmp = TO_LOWER (*tmp);
        }
     }
  }
/* yanli end, 7/19/1999 */

  if (adp->blocks!=NULL) {
     tdp1 = (TextAlignBufPtr) adp->linebuff->data.ptrvalue;
     ddp =(DenseDiagPtr)adp->blocks->segs;
     if (ddp != NULL) 
     {
        for (; ddp!=NULL; ddp=ddp->next) 
        {
           sip = ddp->id;
           start = *(ddp->starts); 
           len = ddp->len;
           if (start+len > adp->gr.left +adp->bufferstart) 
           {
              start = MAX(start, adp->gr.left +adp->bufferstart);
              start1=SeqCoordToAlignCoord (start, sip, salp, 0, 0);
              chklocp =chkloc (sip, start+len-1, adp->sqloc_list, &stop1);
              stop1 =SeqCoordToAlignCoord (stop1, sip, salp, 0, chklocp);
              for (vnp=adp->linebuff; vnp!=NULL; vnp=vnp->next) 
              {
                 sip = ddp->id;
                 tdp = (TextAlignBufPtr) vnp->data.ptrvalue;
                 seqsip=(SeqIdPtr)SeqIdFromAlignNode(adp->anp_list, tdp->seqEntityID, tdp->bsp_itemID, OBJ_BIOSEQ);
                 while (sip!=NULL)
                 {
                    if (SeqIdForSameBioseq(seqsip, sip) )
                    {
                       for (j=start1; j<=stop1; j++) 
                          if (tdp1->buf[j] != '-') 
                             tdp->buf[j] = TO_UPPER(tdp->buf[j]);
                       break;
                    }
                    sip=sip->next;
                 }
              }
           }
        }      
     }
  }
}}

     if ( !goOn ) {
        adp->firstssp = NULL;
        return;
     }
  }
  if (goOn || adp->firstssp == NULL) {
     arrange_buffer (adp);
  }
}

/***************************************************************
***
****************************************************************/
static void print_scale (FILE *fout, EditAlignDataPtr adp, Int4 hoffset, Int2 leftmargin, Int2 scalelength)
{
  Int4   scal;
  Char   str[128];
  Int4   j;

  for (j = 0; j < leftmargin; j++)
         fprintf (fout, " ");
  for ( j = hoffset; j < hoffset + scalelength; j++) {
         scal = (Int4)(adp->gr.left + 1 + j);
         if (scal % 10 == 0) 
         {
            sprintf (str, "%d", (int)scal);
            fprintf (fout, "     ");
            fprintf (fout, "%5s", str);
         }
         if (adp->columnpcell > 0 && j > hoffset)
            if ((Int2) j % (Int2) adp->columnpcell == 0) 
               fprintf (fout, " ");
  }
  fprintf (fout, "\n");
}

static void print_bars (FILE *fout, EditAlignDataPtr adp, Int4 hoffset, Int2 leftmargin, Int2 scalelength)
{
  Int4   scal;
  Int4   j;

  for (j = 0; j < leftmargin; j++)
         fprintf (fout, " ");
  for ( j = hoffset; j < hoffset + scalelength; j++)  {
         scal = (Int4)(adp->gr.left + 1 + j);
         if (scal % 10 == 0)  {
            fprintf (fout, "|");
         }
         else if (scal % 5 == 0)  {
            fprintf (fout, ".");
         }
         else { 
            fprintf (fout, " ");
         }
         if (adp->columnpcell > 0 && j > hoffset)
            if ((Int2) j % (Int2) adp->columnpcell == 0) fprintf (fout, " ");
  }
  fprintf (fout, "\n");
}

static void print_line (FILE *fout, CharPtr label, CharPtr txt, Int2 leftmargin, Uint1 columnpcell, SeqIdPtr gi_sip, Boolean  html_option, Boolean is_emptyline)
{
  Int4     strlens;
  Int4     j;
  CharPtr  txtp;
  Boolean  goOn=FALSE;

  txtp = txt;
  if (is_emptyline) {
     for (j=0; j< (Int4)StringLen(txt); j++, txtp++) 
        if (*txtp != '-') { 
           goOn=TRUE; break; }
     if (!goOn) return;
  }
  strlens = 0;
  if (label!=NULL)
     strlens = (Int2)StringLen(label);
  if (strlens > 0) {
     if (strlens > leftmargin)  
        label[leftmargin] = '\0';
     fprintf (fout, "%s", label);       
     if (html_option && gi_sip != NULL) 
        fprintf (fout, "</A>");
     for (j = strlens; j < leftmargin; j++)
        fprintf (fout, " ");
  }
  else {
     for (j = 0; j < leftmargin; j++)
        fprintf (fout, " ");
  }
  txtp = txt;
  strlens = StringLen(txt);
  for (j=0; j< strlens; j++, txtp++) {
         fprintf (fout, "%c", (char)(*txtp));
         if (columnpcell && j)
            if ((Int4) (j+1) % (Int4) columnpcell == 0) {
               fprintf (fout, " ");
            }
  }
  fprintf (fout, "\n");
}

static CharPtr restrict_todiff (CharPtr str1, CharPtr str0)
{
  size_t   j;

  if (str1 != NULL && str0 != NULL) {
     for (j=0; j<MIN(StringLen(str1), StringLen(str0)); j++) {
        if (str1[j] == str0[j]) str1[j] = '.';
     }
  }
  return str1;
}



/************************************************
****  get_master sequence  
************************************************/
NLM_EXTERN CharPtr get_master (ValNodePtr linebuff, Uint2 entityID, Uint4 itemID, Uint2 itemtype)
{
  ValNodePtr      vnp;
  TextAlignBufPtr tap;
  Uint2           tentityID;

  if ( linebuff == NULL ) return NULL;
  for (vnp = linebuff; vnp != NULL; vnp = vnp->next)
  {
         tap = (TextAlignBufPtr) vnp->data.ptrvalue;
         if ( tap != NULL)
         {
            if (OBJ_(tap->feattype) == OBJ_BIOSEQ) 
                 tentityID = tap->seqEntityID;
            else tentityID = tap->entityID;
            if (tentityID == entityID && tap->bsp_itemID == itemID && OBJ_(tap->feattype) == itemtype)
               break;
         }
  }
  if (vnp==NULL || tap == NULL) return NULL;
  if (tap->buf == NULL) return NULL; 
  return (tap->buf);
}


/************************************************************************
***  read_buffer_fromalignnode
*************************************************************************/
NLM_EXTERN Boolean read_buffer_fromalignnode (EditAlignDataPtr adp, ValNodePtr *linebuff, Int4 bufferstart, Int4 minbufferlength, Int2 *numberalignline)
{
  AlignNodePtr  anp;
  ValNodePtr    curr;      /*for the list of AlignNodePtr*/
  ValNodePtr    list;      /*list of DrawText*/
  ValNodePtr    tdp_list;  /*list of DrawText*/
  ValNodePtr    linebufftmp = NULL; 
  TextAlignBufPtr tdp;
  Int4          p_stop=0;

  *linebuff = NULL;
  *numberalignline = 0;
  if ( adp->anp_list == NULL ) {
         ErrPostEx (SEV_ERROR, 0, 0, "fail in read_buffer_fromalignnode [1]");
  }
  for(curr = adp->anp_list; curr !=NULL; curr = curr->next)
  {
         anp = (AlignNodePtr) curr->data.ptrvalue;
         if ( anp == NULL ) {
                ErrPostEx (SEV_ERROR, 0, 0, "fail in read_buffer_fromalignnode [2]");
                break;
         }
                /*generate the DrawText buffer*/
         list = (ValNodePtr) ProcessTextAlignNode (anp, 
                (Int4) (adp->gr.left +bufferstart), 
                (Int4) (adp->gr.left + bufferstart + minbufferlength), 
                &p_stop, NULL, adp->visibleWidth, (Int1)0, (Uint4)0, NULL);
         if ( list == NULL ) {
                ErrPostEx (SEV_ERROR, 0, 0, "fail in read_buffer_fromalignnode [3]");
                break;
         }
	 tdp_list = list;
         while ( tdp_list != NULL )
	 {
                tdp = (TextAlignBufPtr) tdp_list->data.ptrvalue;
                if(tdp->buf != NULL)
                {
                       if (tdp->label == NULL) {
                          tdp->label=(CharPtr)MemNew((size_t)(64*sizeof(Char)));
                          emptystring (tdp->label, 64);
                          StringNCpy_0 (tdp->label, "unknown", 7);
                       }
                       (*numberalignline)++;
                       ValNodeAddPointer (&linebufftmp, 0, (Pointer) tdp);
                }
                tdp_list->data.ptrvalue = NULL;  /*<<<<<<<<<NEW */
                tdp_list = tdp_list->next;
         }
         ValNodeFree (list);                     /*<<<<<<<<<NEW */
  }
  if (linebufftmp == NULL)
     return FALSE;
  *linebuff = linebufftmp;
  return TRUE;
}


static void AdjustAlignmentsInAnnot (SeqAnnotPtr PNTR p_sap, SeqLocPtr slp)
{
  SeqAnnotPtr sap, sap_prev = NULL, sap_next;
  SeqAlignPtr salp, salp_prev = NULL, salp_next;

  if (p_sap == NULL || *p_sap == NULL || slp == NULL) {
    return;
  }

  sap = *p_sap;
  while (sap != NULL) {
    sap_next = sap->next;
    if (sap->type == 2) {
      salp = (SeqAlignPtr) sap->data;
      salp_prev = NULL;
      while (salp != NULL) {
        salp_next = salp->next;
        if (SeqAlignDeleteByLoc  (slp, salp) == NULL) {
          if (salp_prev == NULL) {
            sap->data = salp_next;
          } else {
            salp_prev->next = salp_next;
          }
        } else {
          salp->saip = SeqAlignIndexFree (salp->saip);
          AlnMgr2IndexSeqAlign (salp);

          salp_prev = salp;
        }
        salp = salp_next;
      }
      if (sap->data == NULL) {
        if (sap_prev == NULL) {
          *p_sap = sap_next;
        } else {
          sap_prev->next = NULL;
        }
        sap->next = NULL;
        sap = SeqAnnotFree (sap);
      } else {
        sap_prev = sap;
      }
    } else {
      sap_prev = sap;
    }
    sap = sap_next;
  }
}


/**********************************
***
*** BioseqTrimN
***   truncates the nnn's at the extremities of a bioseq bsp
***   providing TopSeqEntry sep allows to modify the SeqAlign if any
***
***********************************/
NLM_EXTERN void SeqAlignDeleteByLocCallback (SeqEntryPtr sep, Pointer mydata,
                                          Int4 index, Int2 indent)
{
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  SeqLocPtr          slp;

  slp = (SeqLocPtr)(mydata);
  if (slp!=NULL && sep != NULL && sep->data.ptrvalue) {
     if (IS_Bioseq(sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp!=NULL) {
           AdjustAlignmentsInAnnot (&(bsp->annot), slp);
        }
     }
     else if(IS_Bioseq_set(sep)) {
        bssp = (BioseqSetPtr)sep->data.ptrvalue;
        if (bssp!=NULL) {
           AdjustAlignmentsInAnnot (&(bssp->annot), slp);
        }
     }
  }
}

NLM_EXTERN Boolean BioseqTrimN (BioseqPtr bsp, SeqEntryPtr sep)
{
  SeqIdPtr      sip;
  SeqLocPtr     slp1 = NULL,
                slp2 = NULL;
  CharPtr       str;
  Int4          j,
                lens;
  Boolean       truncate = FALSE;
  
  if (bsp==NULL)
     return FALSE;
  sip = bsp->id;
  str = load_seq_data (sip, -1, -1, FALSE, &lens);
  if (str != NULL) 
  {
     j = lens-1;
     while (j>0) {
        if (str[j] != 'n' && str[j] != 'N') 
           break;
        j--;
     }
     if (j<lens-1) 
     {
        slp1 = SeqLocIntNew (j+1, lens-1, Seq_strand_plus, sip);
        SeqDeleteByLoc (slp1, TRUE, FALSE);
        truncate = TRUE;
     }
     j=0;
     while (j<lens) {
        if (str[j] != 'n' && str[j] != 'N') 
           break;
        j++;
     }
     if (j>0) {
        slp2 = SeqLocIntNew (0, j-1, Seq_strand_plus, sip);
        SeqDeleteByLoc (slp2, TRUE, FALSE);
        truncate = TRUE;
     }
     if (slp1!=NULL) {
        if (sep!=NULL)
           SeqEntryExplore (sep, (Pointer)slp1, SeqAlignDeleteByLocCallback);
        ValNodeFree (slp1);
     }
     if (slp2!=NULL) {
        if (sep!=NULL)
           SeqEntryExplore (sep, (Pointer)slp2, SeqAlignDeleteByLocCallback);
        ValNodeFree (slp2);
     }
  }
  return truncate;
}


/**********************************************************
***  GetFeatureForEditor
***
***********************************************************/
static Boolean is_newfeat_static (ValNodePtr feathead, Uint2 eID, Uint2 subtype, SeqLocPtr slp)
{
  SelEdStructPtr   psp;
  ValNodePtr       vnp = NULL;

  if (feathead == NULL)  { 
         return TRUE;
  }
  for (vnp = feathead; vnp != NULL; vnp = vnp->next)
  {
        psp = (SelEdStructPtr) vnp->data.ptrvalue;
        if (psp->entityID == eID && psp->itemtype == subtype)
        {
           if (SeqLocCompare(slp, (SeqLocPtr) psp->region) == SLC_A_EQ_B) {
              return FALSE;
           }
        }
  }
  return TRUE;
}

NLM_EXTERN ValNodePtr AddFeatFunc (SelEdStructPtr feat, ValNodePtr *featlist, Uint2 itemsubtype)
{
  SelEdStructPtr   psp, 
                   prepsp, tmp;
  ValNodePtr       feathead;
  ValNodePtr       vnp = NULL;
  Int4             featstart, pspstart, pspnextstart;
  Int1             insert;
  Uint4            itemID;

  if (feat == NULL) 
         return *featlist;
  itemID = feat->itemID;

  feathead = *featlist;
  if (feathead == NULL)
  {
         feathead = ValNodeNew (NULL);
         if (feathead == NULL) {
                return feathead;
         }
         feathead->choice = (Uint1)itemsubtype;
         feathead->data.ptrvalue = (Pointer) feat;
         feat->prev = NULL;
         return feathead;
  }
  if (feathead != NULL) 
  {
         for (vnp = feathead; vnp != NULL; vnp = vnp->next)
         {
            if (vnp->choice == itemsubtype) 
            {
               psp = (SelEdStructPtr) vnp->data.ptrvalue;
               if (is_sameses (psp, feat))
               {
                 insert = 0;
                 featstart = SeqLocStart((SeqLocPtr)feat->region);
                 prepsp = NULL;
                 while (psp!= NULL) 
                 {
                   if (psp->next == NULL) {
                      if (itemID == psp->itemID) {
                         if (overlapp_ssp ((SeqLocPtr)feat->region, (SeqLocPtr) psp->region)
) { 
/*{
CharPtr str1, str2;
str1 = SeqLocPrint((SeqLocPtr) feat->region);
str2 = SeqLocPrint((SeqLocPtr) psp->region);
WriteLog("%d %d   %d %d>>> %s  >> %s\n", feat->entityID, feat->itemID, psp->entityID, psp->itemID, str1, str2);
}*/
                            insert = 99; /*break; */
                         }
                         pspstart = SeqLocStart((SeqLocPtr)psp->region);
                         if ( featstart < pspstart) { 
                            insert=-1; break; }
                         else { 
                            insert = +1; break; }
                      }
                      if (itemID < psp->itemID) { insert=-1; break; }
                      if (itemID > psp->itemID) { insert=+1; break; }
                   }
                   else {
                      pspstart = SeqLocStart((SeqLocPtr)psp->region);
                      pspnextstart = SeqLocStart((SeqLocPtr)psp->next->region);
                      if (itemID == psp->itemID ) { 
                         if (overlapp_ssp ((SeqLocPtr)feat->region, (SeqLocPtr) psp->region)) { 
/*{
CharPtr str1, str2;
str1 = SeqLocPrint((SeqLocPtr) feat->region);
str2 = SeqLocPrint((SeqLocPtr) psp->region);
WriteLog("%d %d   %d %d>>> %s  >> %s\n", feat->entityID, feat->itemID, psp->entityID, psp->itemID, str1, str2);
}*/
                            insert = 99; break; }
                         if (itemID == psp->itemID && featstart < pspstart) { 
                            insert=-1; break; 
                         }
                         if (itemID == psp->next->itemID) { 
                            if (overlapp_ssp ((SeqLocPtr)feat->region, (SeqLocPtr) psp->next->region)) { 
/*{
CharPtr str1, str2;
str1 = SeqLocPrint((SeqLocPtr) feat->region);
str2 = SeqLocPrint((SeqLocPtr) psp->region);
WriteLog("%d %d   %d %d>>> %s  >> %s\n", feat->entityID, feat->itemID, psp->entityID, psp->itemID, str1, str2);
}*/
                               insert = 99; break; }
                            if(featstart >pspstart && featstart <pspnextstart) {
                               insert=+1; break; 
                            }
                         }
                         if (itemID < psp->next->itemID) { 
                            if (featstart > pspstart ) {
                               insert=+1; break; 
                            }
                         }
                      }
                      if (itemID >psp->itemID && itemID ==psp->next->itemID) { 
                         if (overlapp_ssp((SeqLocPtr)feat->region, (SeqLocPtr) psp->next->region)) {
/*{
CharPtr str1, str2;
str1 = SeqLocPrint((SeqLocPtr) feat->region);
str2 = SeqLocPrint((SeqLocPtr) psp->region);
WriteLog("%d %d   %d %d>>> %s  >> %s\n", feat->entityID, feat->itemID, psp->entityID, psp->itemID, str1, str2);
}*/
                            insert=99; break; }
                         if (featstart < pspnextstart ) {
                            insert=+1; break; 
                         }
                      }
                      if (itemID <psp->itemID) { 
                         insert=-1; break; 
                      }
                   }
                   prepsp = psp;
                   psp = psp->next;
                 }
                 if (insert == 99) {
/* exons overlapp */
                    MemFree (feat);
                    feat=NULL;
                    return feathead;
                 }
                 if (insert < 0 && prepsp == NULL) {
                    feat->next = psp;
                    vnp->data.ptrvalue = (Pointer) feat;
                    psp->prev = feat;
                    feat->prev = NULL;
                 }
                 else if (insert < 0 && prepsp != NULL) {
                    tmp = prepsp->next;
                    prepsp->next = feat;
                    feat->next = tmp;
                    psp->prev = feat;
                    feat->prev = prepsp;
                 }
                 else  if (insert > 0 && psp->next == NULL) {
                    psp->next = feat;
                    feat->prev = psp;
                 }
                 else  if (insert > 0 && psp->next != NULL) {
                    tmp = psp->next;
                    psp->next = feat;
                    feat->next = tmp;
                    tmp->prev = feat;
                    feat->prev = prepsp;
                 }
                 else 
                    ErrPostEx (SEV_ERROR, 0, 0, "fail in MadeFeatProc [99]"); 
                 break;
               }
            }
         }
  }
  if (vnp == NULL) {
         vnp = ValNodeAddPointer (&feathead, 0, (Pointer) feat);
         vnp->choice = (Uint1)itemsubtype;
         return feathead;
  }
  return feathead;
}

/****************************************************************
*
*   satcollfunc()
*   callback function for collecting features on Sequence 
*   alignment. It recalculates the feature intervals based on 
*   the intervals in the aligned segments
*
****************************************************************/
typedef struct alignfeat
{
   ObjMgrPtr   omp;
   ValNodePtr  data;
   Int2        filter_level;
   Uint2       entityID;
   CollectSeqOptionPtr csop;
   BioseqPtr   bsp;

}AlignFeat, PNTR AlignFeatPtr;

static Boolean slpfeatcollfunc(GatherContextPtr gcp)
{
  CollectSeqOptionPtr csop;
  ObjMgrTypePtr  omtp;
  SeqFeatPtr     sfp;
  AlignFeatPtr   afp;
  SeqLocPtr      slp = NULL;
  SeqLocPtr      curr;
  SelEdStructPtr feat;
  SeqIdPtr       sip;
  BioseqPtr      bsp;
  CdRegionPtr    crp;
  Char           label[101];
  Int4           start, stop;
  Uint2          eID, bspID;
  Uint4          iID;
  Uint2          feat_subtype;   /*types defined by objfdef.h*/
  Int2           label_size;
  Uint1          strand;
  Uint1          codonstart;

  if(gcp->thistype != OBJ_SEQFEAT)
      return TRUE;
  afp = (AlignFeatPtr)(gcp->userdata);
  if(afp == NULL || afp->csop == NULL)
      return FALSE;
  if(afp->filter_level == gcp->seglevel+1)
      return TRUE;
  csop = afp->csop;
  if(csop->features == NULL)
      return FALSE;

  omtp = ObjMgrTypeFind (afp->omp, OBJ_SEQFEAT, NULL, NULL);
  if(omtp == NULL)
      return TRUE;

  feat_subtype = 0;
  if(omtp->subtypefunc !=NULL)
      feat_subtype = (*(omtp->subtypefunc)) (gcp->thisitem); 
   /*do not collect the current feature*/
  if(csop->features[feat_subtype] == FALSE)
      return TRUE;
  sfp = (SeqFeatPtr)gcp->thisitem;
/*
  label_size = MIN(100, csop->label_size);
*/
  label_size = 50;
  label [0] = '\0';
  if (omtp->labelfunc != NULL)
  {
        (*(omtp->labelfunc))(sfp, label, label_size, csop->flabel_format [feat_subtype] );
  }
  if(gcp->product)
     slp = sfp->product;
  else {
     slp = sfp->location;
     eID = gcp->entityID;
     iID = gcp->itemID;
     bspID = afp->entityID;
     if (sfp->data.choice == SEQFEAT_CDREGION) {
        crp = (CdRegionPtr) sfp->data.value.ptrvalue;
        codonstart = crp->frame;
     }
     else codonstart = 0;
     if (slp->choice == SEQLOC_PACKED_INT || slp->choice == SEQLOC_MIX)
     {
        strand = SeqLocStrand (slp);
        curr = NULL;
        while ((curr = SeqLocFindNext(slp, curr)) != NULL)
        {
           bsp = afp->bsp;
           if (bsp != NULL) {
                 start = GetOffsetInBioseq (curr, bsp, SEQLOC_LEFT_END);
                 stop = GetOffsetInBioseq (curr, bsp, SEQLOC_RIGHT_END);
                 if (start != -1 || stop != -1) 
                 {
                    start = MAX (SeqLocStart (curr), 0);
                    stop = MIN ((bsp->length - 1), SeqLocStop (curr));
                    sip = SeqIdFindBest (bsp->id, 0);
                    feat = new_seledstruct (eID, iID, OBJ_SEQFEAT, 0,bspID,  
                                  start, stop,  sip, strand, FALSE, label, NULL, 0 , codonstart);
                    afp->data=AddFeatFunc(feat, &(afp->data), feat_subtype);
                 }
           }
           codonstart = 1;
        }
     }
     else  {
        bsp = afp->bsp;
        if (bsp != NULL) {
           start = GetOffsetInBioseq (slp, bsp, SEQLOC_LEFT_END);
           stop = GetOffsetInBioseq (slp, bsp, SEQLOC_RIGHT_END);
           if (start != -1 || stop != -1) 
           {
                 start = MAX (SeqLocStart (slp), 0);
                 stop = MIN ((bsp->length - 1), SeqLocStop (slp));
                 sip = SeqIdFindBest (bsp->id, 0);
                 feat = new_seledstruct (eID, iID, OBJ_SEQFEAT, 0,bspID, 
                               start, stop,  sip, SeqLocStrand(slp), FALSE, label, NULL, 0, codonstart);
                 afp->data = AddFeatFunc (feat, &(afp->data), feat_subtype);
           }
        }
     }
  }
  return TRUE;
}

/******************************************************************
*
*    CollectFeatureForEditor (slp, anp, csop)
*   collect feature for the alignment
*   slp: the target Seq-loc
*   anp: the AlignNode belong to the target Seq-loc
*   csop: the option for gathering the features
*   
******************************************************************/
NLM_EXTERN ValNodePtr CollectFeatureForEditor (SeqLocPtr slp, ValNodePtr seqfeat, Uint2 seq_entityID, Uint4 bsp_itemID, Uint1 *featOrder, Boolean all_feat)
{
  CollectSeqOption cs_option;
  GatherScope      gs;
  AlignFeat        af;
  BioseqPtr        bsp;
  ValNodePtr       vnp;
  ValNodePtr       vnp2, next;
  SelEdStructPtr   psp;
  SeqEntryPtr      sep,
                   old;
  Int2             j;
  Boolean          show_feature, collect = FALSE;
   
  if(slp == NULL)
      return FALSE;
  if(seq_entityID == 0)
      return FALSE;
  cs_option.nointerval = FALSE;
  cs_option.slabel_format = PRINTID_FASTA_LONG;   /*PRINTID_TEXTID_ACCESSION;*/
  cs_option.seglevels = 0;
  for( j = 0; j < FEATDEF_ANY; ++j)   
  {
     if (all_feat)
        show_feature = TRUE;
     else
        show_feature = (Boolean)(featOrder[j] != 0);
     cs_option.features[j] = show_feature;
     if(show_feature) 
        collect = TRUE;
  }
  if(collect)
  {
     MemSet ((Pointer) (cs_option.flabel_format), OM_LABEL_CONTENT, (size_t) FEATDEF_ANY*sizeof(Uint1));
     bsp = BioseqLockById(SeqLocId(slp));

/*****/
sep =  GetTopSeqEntryForEntityID (seq_entityID);
old = SeqEntrySetScope (sep);
/****/
        
     MemSet((Pointer)&gs, 0, sizeof (GatherScope));
     gs.get_feats_location = TRUE;
     gs.get_feats_product =( bsp->mol == Seq_mol_aa);
     MemSet ((Pointer)(gs.ignore), (int)TRUE, (size_t) (OBJ_MAX) * sizeof (Boolean));
     gs.ignore[OBJ_SEQANNOT] = FALSE;
     gs.ignore[OBJ_SEQFEAT] = FALSE;
     gs.nointervals = FALSE;   
     gs.ignore_top = FALSE;
     gs.currlevel = 0;
     gs.offset = 0;               /*anp->extremes.left;*/
     gs.target = slp;
     gs.scope = sep;

     af.data = NULL;
     af.csop = &cs_option;
     af.omp = ObjMgrGet();
     af.filter_level = 0;
     af.entityID = bsp_itemID;
     af.bsp = bsp;
     GatherEntity (seq_entityID, (Pointer)(&af), slpfeatcollfunc, &gs);
     BioseqUnlock (bsp);

SeqEntrySetScope (old);

     if (seqfeat != NULL) {
        for (vnp=seqfeat; vnp->next!=NULL; vnp=vnp->next)
           continue;
        vnp2 = af.data;
        while (vnp2!=NULL) {
           next = vnp2->next;
           vnp2->next = NULL;
           psp = (SelEdStructPtr) vnp2->data.ptrvalue;
           if (is_newfeat_static (seqfeat, psp->entityID, 0, (SeqLocPtr)psp->region))
           {
              vnp->next = vnp2;
              vnp = vnp2;
           }
           else {
              SelEdStructDel (psp); 
              vnp2->data.ptrvalue = NULL;
              ValNodeFree (vnp2);
           }
           vnp2 = next;
        }
     }
     else seqfeat = af.data;
     return seqfeat;
  }
  return NULL;
}


/**********************************************************************
*** ShowAlignmentText
***
**********************************************************************/

static void print_firstlinePHYLIP (FILE *fout, Int2 seqnumber, Int4 length)
{
  fprintf (fout, "%7d%5d\n", seqnumber, length);
}

NLM_EXTERN void ShowAlignmentText (FILE *fout, EditAlignDataPtr adp, SelStructPtr ssp, Int2 leftmargin, Int4 printfrom, Int4 printto, Boolean html_option)
{
  ValNodePtr      linebuff=NULL;
  TextAlignBufPtr tdp;
  AlignNodePtr    anp;
  SeqIdPtr        bspsip;
  SeqAlignPtr     salp = (SeqAlignPtr) adp->sap_align->data;
  CharPtr         bufstr,
                  trans;
  CharPtr         masterbuf =NULL;
  Int4            width;
  Int4            widthtmp;
  Int4            from; 
  Uint4           itemID; 
  Uint2           entityID, itemtype;
  Int2            numberalignline = 0;
  Int2            index = 0;
  Int2            j, k;
  Boolean         firstline = TRUE;

  if (printfrom > adp->length || printto > adp->length) {
     Message(MSG_OK, "fail in ShowAlignment: %ld > %ld or %ld > %ld", (long)printfrom, (long)adp->length, (long)printto, (long)adp->length);
     return;
  }  
  if (adp->align_format == SALSA_PHYLIP)
     print_firstlinePHYLIP (fout, adp->seqnumber, adp->length);
  from = printfrom;
  width = (Int4)MIN ((Int4) adp->visibleWidth, (Int4)(printto - from +1));

  while ( read_buffer_fromalignnode (adp, &linebuff, from, width-1, &numberalignline) )
  {
     if (adp->charmode)
        masterbuf = get_master (linebuff, adp->master.entityID, adp->master.itemID, adp->master.itemtype);
     else masterbuf = NULL;

     if (adp->draw_scale && (adp->align_format != SALSA_PHYLIP)) {
         print_scale (fout, adp, from, leftmargin, (Int2)width);
     }
     if (adp->draw_bars && (adp->align_format != SALSA_PHYLIP)) {
         print_bars (fout, adp, from, leftmargin, (Int2)width);
     }
     index = 1;
     widthtmp = (Int4)width;
     bufstr = next_notemptyline (adp->anp_list, linebuff, numberalignline, 
                                &index, (Int4)0, &widthtmp, &tdp, &anp);
     while ( index <= numberalignline && bufstr != NULL) 
     {
         if ( OBJ_(tdp->feattype) == OBJ_BIOSEQ ) {
              itemID = anp->bsp_itemID;
              if (adp->input_format == OBJ_BIOSEQ) {
                 entityID = tdp->seqEntityID;
              }
              else {
                 entityID = anp->seq_entityID;
              }
              bspsip = anp->sip;
         } 
         else {
              itemID = tdp->itemID;
              entityID = tdp->entityID;
         }
         itemtype = OBJ_(tdp->feattype);
         if (not_empty_string (bufstr, width)) 
         {
            if (masterbuf != NULL && (adp->master.entityID != entityID || adp->master.itemID != itemID))
               bufstr = restrict_todiff (bufstr, masterbuf);
            if ((adp->align_format == SALSA_PHYLIP) && !firstline)
               print_line (fout, NULL, bufstr, leftmargin, adp->columnpcell, bspsip, html_option, FALSE);
            else
               print_line (fout, tdp->label, bufstr, leftmargin, adp->columnpcell, bspsip, html_option, FALSE);
         }
         if (OBJ_(tdp->feattype) == OBJ_BIOSEQ && not_empty_string (bufstr, width)) 
         {
            if (has_complement (adp->params, entityID, itemID)) {
               print_line (fout, "complement", complement_string(bufstr), leftmargin, adp->columnpcell, NULL, html_option, FALSE);
            }
            for (j=0; j<3;j++) {
               if (rf_on (adp->params, entityID, itemID, j))  {
                   trans = (CharPtr) SeqAlignTranslate (salp, entityID, from, from + widthtmp - 1,  (Uint1)(j), bspsip, Seq_strand_plus, adp->sqloc_list);
                   if ( adp->prot_mode == PUTPROT )
                      trans = prot_to_putprot (trans);
                   print_line (fout, "RF >", trans, leftmargin, adp->columnpcell, NULL, html_option, FALSE);
               }
            }
            for (j=3, k=2; j<6;j++, k--) {
               if (rf_on (adp->params, entityID, itemID, j))  {
                   trans = (CharPtr) SeqAlignTranslate (salp, entityID, from, from + widthtmp - 1,  (Uint1)k, bspsip, Seq_strand_minus, adp->sqloc_list);
                   if ( adp->prot_mode == PUTPROT )
                      trans = prot_to_rputprot (trans);
                   print_line (fout, "RF <", trans, leftmargin, adp->columnpcell, NULL, html_option, FALSE);
               }
            }
         }
         index++;
         widthtmp = width;
         bufstr = next_notemptyline (adp->anp_list, linebuff, numberalignline, 
                                     &index, (Int4)0, &widthtmp, &tdp, &anp);
     }
     firstline = FALSE;
     from += width;
     if (from >= printto) break;
     width = (Int4)MIN ((Int4) adp->visibleWidth, (Int4)(printto - from +1));
     fprintf (fout, "\n");
  }
  return;
}

NLM_EXTERN void showfastagap_fromalign (SeqAlignPtr salp, Int4 line, FILE *f)
{
  BioseqPtr bsp;
  CharPtr   str,
            strp;
  Char      buffer [255];
  CharPtr   bufp;
  SeqIdPtr  sip;
  Int4      len,
            al_len;
  Int4      from, to, offset,
            ncar = 0;
  Int2      j;
  Boolean   goOn = TRUE,
            is_prot;

  Char      strLog[120];

  if (salp == NULL)
     return;
  for (j=0; j<salp->dim; j++) {
     sip = SeqAlignId (salp, j);
     if (sip!=NULL) {
      bsp = BioseqLockById (sip);
      if (bsp!=NULL) {
        is_prot = (Boolean) (ISA_aa (bsp->mol));
        SeqIdWrite (bsp->id, strLog, PRINTID_FASTA_LONG, 120);
        BioseqUnlock (bsp);
        fprintf(f, ">%s\n", strLog);
        str = load_seq_data (sip, -1, -1, is_prot, &len);
        strp = str;
        al_len = SeqAlignLength (salp);
        from = 0;
        offset = (Int4) MIN (line, (al_len-from));
        to = offset -1;
        goOn = (Boolean) (to > from);
        while (goOn) {
/*           fprintf(f, "%5d  ", from); no more line numbers! SW 11/01 */
           bufp = ReadBufferFromSap (strp, buffer, salp, sip, from, to, &ncar);
           fprintf(f, "%s\n", bufp);
           strp+= ncar;
           from = to+1;
           offset = (Int4) MIN (line, (al_len-from));
           to = from +offset -1;
           goOn = (Boolean) (to > from);
        }
      }
     }
  }
}


