/*   sqnutil1.c
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
* File Name:  sqnutil1.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   9/2/97
*
* $Revision: 6.676 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#include <sqnutils.h>
#include <gather.h>
#include <subutil.h>
#include <objfdef.h>
#include <seqport.h>
#include <objproj.h>
/* #include <objmmdb1.h> */
#include <gbfeat.h>
#include <gbftdef.h>
#include <edutil.h>
#include <tofasta.h>
#include <parsegb.h>
#include <utilpars.h>
#include <validatr.h>
#include <explore.h>
#include <subutil.h>
#include <asn2gnbi.h>
#include <salpacc.h>
#include <alignmgr2.h>
#include <valid.h>

#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>

static int descr_insert_order [] = {
  Seq_descr_title,
  Seq_descr_source,
  Seq_descr_molinfo,
  Seq_descr_het,
  Seq_descr_pub,
  Seq_descr_comment,
  Seq_descr_name,
  Seq_descr_user,
  Seq_descr_maploc,
  Seq_descr_region,
  Seq_descr_num,
  Seq_descr_dbxref,
  Seq_descr_mol_type,
  Seq_descr_modif,
  Seq_descr_method,
  Seq_descr_org,
  Seq_descr_sp,
  Seq_descr_pir,
  Seq_descr_prf,
  Seq_descr_pdb,
  Seq_descr_embl,
  Seq_descr_genbank,
  Seq_descr_modelev,
  Seq_descr_create_date,
  Seq_descr_update_date,
  0
};

static void NormalizeDescriptorProc (
  SeqEntryPtr sep,
  Pointer data,
  Int4 index,
  Int2 indent
)

{
  BioseqPtr         bsp;
  BioseqSetPtr      bssp;
  /* arrays are SEQDESCR_MAX + 1, last slot stores unexpected descriptor numbers */
  SeqDescrPtr       first [SEQDESCR_MAX + 1];
  SeqDescrPtr       last [SEQDESCR_MAX + 1];
  int               i;
  int               idx;
  SeqDescrPtr PNTR  head = NULL;
  SeqDescrPtr PNTR  prev = NULL;
  SeqDescrPtr       next;
  SeqDescrPtr       sdp;

  if (sep == NULL) return;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    head = &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    head = &(bssp->descr);
  }
  if (head == NULL) return;

  MemSet ((Pointer) &first, 0, sizeof (first));
  MemSet ((Pointer) &last, 0, sizeof (last));

  prev = head;
  sdp = *prev;
  while (sdp != NULL) {
    next = sdp->next;

    *prev = sdp->next;
    sdp->next = NULL;

    idx = (int) sdp->choice;
    /* unexpected descriptor numbers go into last slot */
    if (idx <= 0 || idx >= SEQDESCR_MAX) {
      idx = SEQDESCR_MAX;
    }
    if (idx > 0 && idx <= SEQDESCR_MAX) {
      if (first [idx] == NULL) {
        first [idx] = sdp;
      }
      if (last [idx] != NULL) {
        (last [idx])->next = sdp;
      }
      last [idx] = sdp;
    }

    sdp = next;
  }

  for (i = 0; descr_insert_order [i] != 0; i++) {
    idx = descr_insert_order [i];
    sdp = first [idx];
    if (sdp == NULL) continue;
    ValNodeLink (head, sdp);
  }
}

NLM_EXTERN void NormalizeDescriptorOrder (
  SeqEntryPtr sep
)

{
  SeqEntryExplore (sep, NULL, NormalizeDescriptorProc);
}

NLM_EXTERN DatePtr DateAdvance (DatePtr dp, Uint1 monthsToAdd)

{
  if (dp == NULL) {
    dp = DateCurr ();
  }
  if (dp != NULL && dp->data [0] == 1 && dp->data [1] > 0) {
    while (monthsToAdd > 12) {
      monthsToAdd--;
      (dp->data [1])++;
    }
    if (dp->data [2] < 13 - monthsToAdd) {
      (dp->data [2]) += monthsToAdd;
    } else {
      (dp->data [1])++;
      (dp->data [2]) -= (12 - monthsToAdd);
    }
    if (dp->data [2] == 0) {
      dp->data [2] = 1;
    }
    if (dp->data [3] == 0) {
      switch (dp->data [2]) {
        case 4 :
        case 6 :
        case 9 :
        case 11 :
          dp->data [3] = 30;
          break;
        case 2 :
          dp->data [3] = 28;
          break;
        default :
          dp->data [3] = 31;
          break;
      }
    }
  }
  if (dp != NULL) {
    switch (dp->data [2]) {
      case 4 :
      case 6 :
      case 9 :
      case 11 :
        if (dp->data [3] > 30) {
          dp->data [3] = 30;
        }
        break;
      case 2 :
        if (dp->data [3] > 28) {
          dp->data [3] = 28;
        }
        break;
      default :
        if (dp->data [3] > 31) {
          dp->data [3] = 31;
        }
        break;
    }
  }
  return dp;
}

typedef struct orgscan {
  ObjMgrPtr     omp;
  Int2          nuclCode;
  Int2          mitoCode;
  Int2          pstdCode;
  Boolean       mito;
  Boolean       plastid;
  Char          taxname [196];
  BioSourcePtr  biop;
} OrgScan, PNTR OrgScanPtr;

static Boolean OrgScanGatherFunc (GatherContextPtr gcp)

{
  BioSourcePtr   biop;
  Boolean        doCodes = FALSE;
  Boolean        doMito = FALSE;
  Boolean        doTaxname = FALSE;
  Boolean        mito = FALSE;
  Int2           mitoCode = 0;
  Int2           nuclCode = 0;
  Int2           pstdCode = 0;
  ObjMgrTypePtr  omtp;
  OrgNamePtr     onp;
  OrgRefPtr      orp;
  OrgScanPtr     osp;
  ValNodePtr     sdp;
  SeqFeatPtr     sfp;
  Uint2          subtype;
  CharPtr        taxname = NULL;
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
              mito = TRUE;
              doMito = TRUE;
              /* osp->mito = TRUE; */
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
    mito = (Boolean) (biop->genome == GENOME_kinetoplast ||
                      biop->genome == GENOME_mitochondrion ||
                      biop->genome == GENOME_hydrogenosome);
    doMito = TRUE;
    /* osp->mito = (Boolean) (biop->genome == 4 || biop->genome == 5); */
  }
  if (orp != NULL) {
    taxname = orp->taxname;
    doTaxname = TRUE;
    /* StringNCpy_0 (osp->taxname, orp->taxname, sizeof (osp->taxname)); */
    onp = orp->orgname;
    if (onp != NULL) {
      nuclCode = onp->gcode;
      mitoCode = onp->mgcode;
      pstdCode = onp->pgcode;
      doCodes = TRUE;
      /* osp->nuclCode = onp->gcode;
      osp->mitoCode = onp->mgcode; */
    }
  }
  if (biop != NULL) {
    if (osp->biop == NULL || biop->is_focus) {
      osp->biop = biop;
      if (doMito) {
        osp->mito = mito;
      }
      osp->plastid = (Boolean) (biop->genome == GENOME_chloroplast ||
                                biop->genome == GENOME_chromoplast ||
                                biop->genome == GENOME_plastid ||
                                biop->genome == GENOME_cyanelle ||
                                biop->genome == GENOME_apicoplast ||
                                biop->genome == GENOME_leucoplast ||
                                biop->genome == GENOME_proplastid ||
                                biop->genome == GENOME_chromatophore);
      if (doCodes) {
        osp->nuclCode = nuclCode;
        osp->mitoCode = mitoCode;
        osp->pstdCode = pstdCode;
      }
      if (doTaxname) {
        StringNCpy_0 (osp->taxname, taxname, sizeof (osp->taxname));
      }
    }
  }

  return TRUE;
}

static Int2 SeqEntryOrEntityIDToGeneticCode (SeqEntryPtr sep, Uint2 entityID, BoolPtr mito,
                                             CharPtr taxname, size_t maxsize,
                                             BioSourcePtr PNTR biopp)

{
  GatherScope  gs;
  OrgScan      osp;

  if (mito != NULL) {
    *mito = FALSE;
  }
  if (taxname != NULL && maxsize > 0) {
    *taxname = '\0';
  }
  osp.mito = FALSE;
  osp.plastid = FALSE;
  osp.nuclCode = 0;
  osp.mitoCode = 0;
  osp.pstdCode = 0;
  osp.omp = ObjMgrGet ();
  osp.taxname [0] = '\0';
  osp.biop = NULL;
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
    GatherSeqEntry (sep, (Pointer) &osp, OrgScanGatherFunc, &gs);
  } else if (entityID > 0) {
    GatherEntity (entityID, (Pointer) &osp, OrgScanGatherFunc, &gs);
  }
  if (mito != NULL) {
    *mito = osp.mito;
  }
  if (taxname != NULL && maxsize > 0) {
    StringNCpy_0 (taxname, osp.taxname, maxsize);
  }
  if (biopp != NULL) {
    *biopp = osp.biop;
  }
  if (osp.plastid) {
    if (osp.pstdCode > 0) {
      return osp.pstdCode;
    } else {
      return 11;
    }
  } else if (osp.mito) {
    return osp.mitoCode;
  } else {
    return osp.nuclCode;
  }
}

NLM_EXTERN Int2 EntityIDToGeneticCode (Uint2 entityID, BoolPtr mito, CharPtr taxname, size_t maxsize)

{
  return SeqEntryOrEntityIDToGeneticCode (NULL, entityID, mito, taxname, maxsize, NULL);
}

NLM_EXTERN Int2 SeqEntryToGeneticCode (SeqEntryPtr sep, BoolPtr mito, CharPtr taxname, size_t maxsize)

{
  return SeqEntryOrEntityIDToGeneticCode (sep, 0, mito, taxname, maxsize, NULL);
}

NLM_EXTERN Int2 SeqEntryToBioSource (SeqEntryPtr sep, BoolPtr mito, CharPtr taxname, size_t maxsize, BioSourcePtr PNTR biopp)

{
  return SeqEntryOrEntityIDToGeneticCode (sep, 0, mito, taxname, maxsize, biopp);
}

NLM_EXTERN Boolean BioseqToGeneticCode (
  BioseqPtr bsp,
  Int2Ptr gencodep,
  BoolPtr mitop,
  BoolPtr plastidp,
  CharPtr taxnamep,
  size_t maxsize,
  BioSourcePtr PNTR biopp
)

{
  BioSourcePtr       biop = NULL;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  Int2               gencode = 0;
  Boolean            mito = FALSE;
  Int2               mitoCode = 0;
  Int2               nuclCode = 0;
  Int2               pstdCode = 0;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  Boolean            plastid = FALSE;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  CharPtr            taxname = NULL;

  if (bsp == NULL) return FALSE;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
  }

  if (biop == NULL) {
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
    if (sfp != NULL) {
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    }
  }

  if (biop == NULL) return FALSE;
  orp = biop->org;
  if (orp == NULL) return FALSE;

  taxname = orp->taxname;
  if (StringHasNoText (taxname)) return FALSE;

  onp = orp->orgname;
  if (onp != NULL) {
    nuclCode = onp->gcode;
    mitoCode = onp->mgcode;
    pstdCode = onp->pgcode;
  }

  mito = (Boolean) (biop->genome == GENOME_kinetoplast ||
                    biop->genome == GENOME_mitochondrion ||
                    biop->genome == GENOME_hydrogenosome);

  plastid = (Boolean) (biop->genome == GENOME_chloroplast ||
                       biop->genome == GENOME_chromoplast ||
                       biop->genome == GENOME_plastid ||
                       biop->genome == GENOME_cyanelle ||
                       biop->genome == GENOME_apicoplast ||
                       biop->genome == GENOME_leucoplast ||
                       biop->genome == GENOME_proplastid ||
                       biop->genome == GENOME_chromatophore);

  if (plastid) {
    if (pstdCode > 0) {
      gencode = pstdCode;
    } else {
      gencode = 11;
    }
  } else if (mito) {
    gencode = mitoCode;
  } else {
    gencode = nuclCode;
  }

  if (gencodep != NULL) {
    *gencodep = gencode;
  }
  if (mitop != NULL) {
    *mitop = mito;
  }
  if (plastidp != NULL) {
    *plastidp = plastid;
  }
  if (taxnamep != NULL && maxsize > 0) {
    StringNCpy_0 (taxnamep, taxname, maxsize);
  }
  if (biopp != NULL) {
    *biopp = biop;
  }

  return TRUE;
}


static Boolean FindBspItem (GatherContextPtr gcp)

{
  BioseqPtr  PNTR bspp;

  bspp = (BioseqPtr PNTR) gcp->userdata;
  if (bspp != NULL && gcp->thistype == OBJ_BIOSEQ) {
    *bspp = (BioseqPtr) gcp->thisitem;
  }
  return TRUE;
}

NLM_EXTERN BioseqPtr GetBioseqGivenIDs (Uint2 entityID, Uint4 itemID, Uint2 itemtype)

{
  BioseqPtr  bsp;

  bsp = NULL;
  if (entityID > 0 && itemID > 0 && itemtype == OBJ_BIOSEQ) {
    GatherItem (entityID, itemID, itemtype, (Pointer) (&bsp), FindBspItem);
  }
  return bsp;
}

NLM_EXTERN BioseqPtr GetBioseqGivenSeqLoc (SeqLocPtr slp, Uint2 entityID)

{
  BioseqPtr    bsp;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;

  if (slp == NULL) return NULL;
  bsp = NULL;
  sip = SeqLocId (slp);
  if (sip != NULL) {
    bsp = BioseqFind (sip);
  } else if (entityID > 0) {
    slp = SeqLocFindNext (slp, NULL);
    if (slp != NULL) {
      sip = SeqLocId (slp);
      if (sip != NULL) {
        bsp = BioseqFind (sip);
        if (bsp != NULL) {
          sep = GetBestTopParentForData (entityID, bsp);
          if (sep != NULL) {
            sep = FindNucSeqEntry (sep);
            if (sep != NULL && sep->choice == 1) {
              bsp = (BioseqPtr) sep->data.ptrvalue;
            }
          }
        }
      }
    }
  }
  return bsp;
}

typedef struct tripletdata {
    Uint2      entityID;
    Uint4      itemID;
    Uint2      itemtype;
    Pointer    lookfor;
} TripletData, PNTR TripletDataPtr;

static Boolean FindIDsFromPointer (GatherContextPtr gcp)

{
  TripletDataPtr  tdp;

  tdp = (TripletDataPtr) gcp->userdata;
  if (tdp != NULL && gcp->thisitem == tdp->lookfor) {
    tdp->entityID = gcp->entityID;
    tdp->itemID = gcp->itemID;
    tdp->itemtype = gcp->thistype;
  }
  return TRUE;
}

NLM_EXTERN Uint4 GetItemIDGivenPointer (Uint2 entityID, Uint2 itemtype, Pointer lookfor)

{
  GatherScope  gs;
  TripletData  td;

  if (entityID > 0 && itemtype > 0 && itemtype < OBJ_MAX && lookfor != NULL) {
    td.entityID = 0;
    td.itemID = 0;
    td.itemtype = 0;
    td.lookfor = lookfor;
    MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
    gs.seglevels = 1;
    gs.get_feats_location = FALSE;
    MemSet ((Pointer)(gs.ignore), (int)(FALSE), (size_t)(OBJ_MAX * sizeof(Boolean)));
    /* gs.ignore[itemtype] = FALSE; */
    GatherEntity (entityID, (Pointer) (&td), FindIDsFromPointer, &gs);
    if (td.entityID == entityID && td.itemID > 0 && td.itemtype == itemtype) {
      return td.itemID;
    }
  }
  return 0;
}

static void AddNucPart (BioseqPtr segseq, BioseqSetPtr parts, SeqEntryPtr addme)

{
  BioseqPtr    bsp;
  SeqLocPtr    slp;
  SeqEntryPtr  tmp;

  if (segseq == NULL || addme == NULL) return;
  if (addme->choice != 1 || addme->data.ptrvalue == NULL) return;
  bsp = (BioseqPtr) addme->data.ptrvalue;

  slp = ValNodeNew ((ValNodePtr) segseq->seq_ext);
  if (slp == NULL) return;
  if (segseq->seq_ext == NULL) {
    segseq->seq_ext = (Pointer) slp;
  }
  if (bsp->length >= 0) {
    segseq->length += bsp->length;
    slp->choice = SEQLOC_WHOLE;
    slp->data.ptrvalue = (Pointer) SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
  } else {
    slp->choice = SEQLOC_NULL;
    addme = SeqEntryFree (addme);
    return;
  }

  if (parts == NULL) {
    addme = SeqEntryFree (addme);
    return;
  }
  if (parts->seq_set != NULL) {
    tmp = parts->seq_set;
    while (tmp->next != NULL) {
      tmp = tmp->next;
    }
    tmp->next = addme;
  } else {
    parts->seq_set = addme;
  }
}

NLM_EXTERN void GetSeqEntryParent (SeqEntryPtr target, Pointer PNTR parentptr, Uint2Ptr parenttype)

{
  ObjMgrPtr      omp;
  ObjMgrDataPtr  omdp;

  if (parentptr == NULL || parenttype == NULL) return;
  *parenttype = 0;
  *parentptr = NULL;
  if (target == NULL || target->data.ptrvalue == NULL) return;
  omp = ObjMgrGet ();
  if (omp == NULL) return;
  omdp = ObjMgrFindByData (omp, target->data.ptrvalue);
  if (omdp == NULL) return;
  *parenttype = omdp->parenttype;
  *parentptr = omdp->parentptr;
}

NLM_EXTERN void SaveSeqEntryObjMgrData (SeqEntryPtr target, ObjMgrDataPtr PNTR omdptopptr, ObjMgrData PNTR omdataptr)

{
  ObjMgrPtr         omp;
  ObjMgrDataPtr  omdp, omdptop = NULL;

  if (target == NULL || omdptopptr == NULL || omdataptr == NULL) return;
  *omdptopptr = NULL;
  MemSet ((Pointer) omdataptr, 0, sizeof (ObjMgrData));
  omp = ObjMgrGet ();
  if (omp == NULL) return;
  omdp = ObjMgrFindByData (omp, target->data.ptrvalue);
  if (omdp == NULL) return;
  omdptop = ObjMgrFindTop (omp, omdp);
  if (omdptop == NULL) return;
  if (omdptop->EntityID == 0) return;
  *omdptopptr = omdptop;
  MemCopy ((Pointer) omdataptr, omdptop, sizeof (ObjMgrData));
  omdptop->userdata = NULL;
}

extern void ObjMgrRemoveEntityIDFromRecycle (Uint2 entityID, ObjMgrPtr omp);
extern void ObjMgrRecordOmdpByEntityID (Uint2 entityID, ObjMgrDataPtr omdp);
NLM_EXTERN void RestoreSeqEntryObjMgrData (SeqEntryPtr target, ObjMgrDataPtr omdptop, ObjMgrData PNTR omdataptr)

{
  ObjMgrPtr    omp;
  ObjMgrDataPtr omdp, omdpnew = NULL;

  if (target == NULL || omdptop == NULL || omdataptr == NULL) return;
  if (omdataptr->EntityID == 0) return;
  omp = ObjMgrGet ();
  if (omp == NULL) return;
  omdp = ObjMgrFindByData (omp, target->data.ptrvalue);
  if (omdp == NULL) return;
  omdpnew = ObjMgrFindTop (omp, omdp);
  if (omdpnew == NULL) return;
  if (omdpnew != omdptop) {
    omdpnew->EntityID = omdataptr->EntityID;
    omdptop->EntityID = 0;
    omdpnew->lockcnt = omdataptr->lockcnt;
    omdpnew->tempload = omdataptr->tempload;
    omdpnew->clipboard = omdataptr->clipboard;
    omdpnew->dirty = omdataptr->dirty;
    omdpnew->being_freed = omdataptr->being_freed;
    omdpnew->free = omdataptr->free;
    omdpnew->options = omdataptr->options;
    ObjMgrRemoveEntityIDFromRecycle (omdpnew->EntityID, omp);
    ObjMgrRecordOmdpByEntityID (omdpnew->EntityID, omdpnew);
  }
  omdpnew->userdata = omdataptr->userdata;
}

NLM_EXTERN void AddSeqEntryToSeqEntry (SeqEntryPtr target, SeqEntryPtr insert, Boolean relink)

{
  SeqEntryPtr    first;
  BioseqPtr      insertbsp;
  BioseqSetPtr   nuc_prot;
  Uint2          parenttype;
  Pointer        parentptr;
  BioseqSetPtr   parts;
  BioseqPtr      seg;
  BioseqSetPtr   segs;
  BioseqPtr      targetbsp;
  BioseqSetPtr   targetbssp;
  SeqEntryPtr    the_nuc;
  SeqEntryPtr    the_prt;
  SeqEntryPtr    tmp;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;

  if (target == NULL || insert == NULL) return;
  if (target->data.ptrvalue == NULL || insert->data.ptrvalue == NULL) return;

  if (relink) {
    SaveSeqEntryObjMgrData (target, &omdptop, &omdata);
    GetSeqEntryParent (target, &parentptr, &parenttype);
  }

  if (IS_Bioseq (target) && IS_Bioseq (insert)) {
    targetbsp = (BioseqPtr) target->data.ptrvalue;
    insertbsp = (BioseqPtr) insert->data.ptrvalue;
    if (ISA_na (targetbsp->mol)) {
      if (ISA_na (insertbsp->mol)) {

        seg = BioseqNew ();
        if (seg == NULL) return;
        seg->mol = targetbsp->mol;
        seg->repr = Seq_repr_seg;
        seg->seq_ext_type = 1;
        seg->length = 0;
        /* seg->id = MakeSeqID ("SEG_dna"); */
        /* seg->id = MakeNewProteinSeqId (NULL, NULL); */
        seg->id = MakeUniqueSeqID ("segseq_");
        SeqMgrAddToBioseqIndex (seg);

        the_nuc = SeqEntryNew ();
        if (the_nuc == NULL) return;
        the_nuc->choice = 1;
        the_nuc->data.ptrvalue = (Pointer) seg;

        segs = BioseqSetNew ();
        if (segs == NULL) return;
        segs->_class = 2;
        segs->seq_set = the_nuc;

        parts = BioseqSetNew ();
        if (parts == NULL) return;
        parts->_class = 4;

        tmp = SeqEntryNew ();
        if (tmp == NULL) return;
        tmp->choice = 2;
        tmp->data.ptrvalue = (Pointer) parts;
        the_nuc->next = tmp;

        first = SeqEntryNew ();
        if (first == NULL) return;
        first->choice = 1;
        first->data.ptrvalue = (Pointer) targetbsp;
        target->choice = 2;
        target->data.ptrvalue = (Pointer) segs;

        AddNucPart (seg, parts, first);
        AddNucPart (seg, parts, insert);

      } else if (ISA_aa (insertbsp->mol)) {

        nuc_prot = BioseqSetNew ();
        if (nuc_prot == NULL) return;
        nuc_prot->_class = 1;

        the_nuc = SeqEntryNew ();
        if (the_nuc == NULL) return;
        the_nuc->choice = 1;
        the_nuc->data.ptrvalue = (Pointer) targetbsp;
        target->choice = 2;
        target->data.ptrvalue = (Pointer) nuc_prot;
        nuc_prot->seq_set = the_nuc;

        the_nuc->next = insert;

      }
    } else if (ISA_aa (targetbsp->mol)) {
      if (ISA_na (insertbsp->mol)) {

        nuc_prot = BioseqSetNew ();
        if (nuc_prot == NULL) return;
        nuc_prot->_class = 1;

        the_prt = SeqEntryNew ();
        if (the_prt == NULL) return;
        the_prt->choice = 1;
        the_prt->data.ptrvalue = (Pointer) targetbsp;
        target->choice = 2;
        target->data.ptrvalue = (Pointer) nuc_prot;
        nuc_prot->seq_set = insert;

        the_prt->next = insert->next;
        insert->next = the_prt;

      }
    }
  } else if (IS_Bioseq_set (target)) {
    targetbssp = (BioseqSetPtr) target->data.ptrvalue;
    if (targetbssp->_class == 1 && IS_Bioseq (insert)) {
     insertbsp = (BioseqPtr) insert->data.ptrvalue;
     if (ISA_aa (insertbsp->mol)) {

        nuc_prot = targetbssp;
        if (nuc_prot->seq_set != NULL) {
          tmp = nuc_prot->seq_set;
          while (tmp->next != NULL) {
            tmp = tmp->next;
          }
          tmp->next = insert;
        } else {
          nuc_prot->seq_set = insert;
        }

      }
    } else if (targetbssp->_class == 2 && IS_Bioseq (insert)) {
      insertbsp = (BioseqPtr) insert->data.ptrvalue;
      if (ISA_na (insertbsp->mol)) {

        the_nuc = FindNucSeqEntry (target);
        if (the_nuc != NULL && the_nuc->next != NULL) {
          tmp = the_nuc->next;
          if (tmp->choice == 2 && tmp->data.ptrvalue != NULL) {
            parts = (BioseqSetPtr) tmp->data.ptrvalue;
            if (parts->_class == 4 && the_nuc->choice == 1) {
              seg = (BioseqPtr) the_nuc->data.ptrvalue;
              AddNucPart (seg, parts, insert);
            }
          }
        }

      } else if (ISA_aa (insertbsp->mol)) {

        nuc_prot = BioseqSetNew ();
        if (nuc_prot == NULL) return;
        nuc_prot->_class = 1;

        first = SeqEntryNew ();
        if (first == NULL) return;
        first->choice = 2;
        first->data.ptrvalue = (Pointer) targetbssp;
        target->choice = 2;
        target->data.ptrvalue = (Pointer) nuc_prot;
        nuc_prot->seq_set = first;

        first->next = insert;

      }
    } else if (targetbssp->_class == 7) {

      if (targetbssp->seq_set != NULL) {
        tmp = targetbssp->seq_set;
        while (tmp->next != NULL) {
          tmp = tmp->next;
        }
        tmp->next = insert;
      } else {
        targetbssp->seq_set = insert;
      }
    } else if ((targetbssp->_class >= BioseqseqSet_class_mut_set &&
                targetbssp->_class <= BioseqseqSet_class_eco_set) ||
               targetbssp->_class == BioseqseqSet_class_wgs_set ||
               targetbssp->_class == BioseqseqSet_class_small_genome_set) {

      if (targetbssp->seq_set != NULL) {
        tmp = targetbssp->seq_set;
        while (tmp->next != NULL) {
          tmp = tmp->next;
        }
        tmp->next = insert;
      } else {
        targetbssp->seq_set = insert;
      }

    } else if (targetbssp->_class == BioseqseqSet_class_gen_prod_set) {

      if (targetbssp->seq_set != NULL) {
        tmp = targetbssp->seq_set;
        while (tmp->next != NULL) {
          tmp = tmp->next;
        }
        tmp->next = insert;
      } else {
        targetbssp->seq_set = insert;
      }

    }
  }

  if (relink) {
    SeqMgrLinkSeqEntry (target, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (target, omdptop, &omdata);
  }
}

NLM_EXTERN void ReplaceSeqEntryWithSeqEntry (SeqEntryPtr target, SeqEntryPtr replaceWith, Boolean relink)

{
  Uint2          parenttype;
  Pointer        parentptr;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;

  if (target == NULL || replaceWith == NULL) return;

  if (relink) {
    SaveSeqEntryObjMgrData (target, &omdptop, &omdata);
    GetSeqEntryParent (target, &parentptr, &parenttype);
  }

  if (target->choice == 1) {
    BioseqFree ((BioseqPtr) target->data.ptrvalue);
  } else if (target->choice == 2) {
    BioseqSetFree ((BioseqSetPtr) target->data.ptrvalue);
  }
  target->choice = replaceWith->choice;
  target->data.ptrvalue = replaceWith->data.ptrvalue;
  MemFree (replaceWith);

  if (relink) {
    SeqMgrLinkSeqEntry (target, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (target, omdptop, &omdata);
  }
}

static void SeqEntryRemoveLoop (SeqEntryPtr sep, SeqEntryPtr del, SeqEntryPtr PNTR prev)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   next;

  while (sep != NULL) {
    next = sep->next;
    if (sep == del) {
      *prev = sep->next;
      sep->next = NULL;
      SeqEntryFree (sep);
    } else {
      prev = (SeqEntryPtr PNTR) &(sep->next);
      if (IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (bssp != NULL) {
          SeqEntryRemoveLoop (bssp->seq_set, del, &(bssp->seq_set));
        }
      }
    }
    sep = next;
  }
}

NLM_EXTERN void RemoveSeqEntryFromSeqEntry (SeqEntryPtr top, SeqEntryPtr del, Boolean relink)

{
  SeqEntryPtr    dummy;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;

  if (top == NULL || del == NULL) return;
  if (top->data.ptrvalue == NULL || del->data.ptrvalue == NULL) return;

  if (relink) {
    SaveSeqEntryObjMgrData (top, &omdptop, &omdata);
    GetSeqEntryParent (top, &parentptr, &parenttype);
  }

  dummy = NULL;
  SeqEntryRemoveLoop (top, del, &dummy);

  if (relink) {
    SeqMgrLinkSeqEntry (top, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (top, omdptop, &omdata);
  }
}


typedef struct commontitle {
  BioseqPtr bsp;
  SeqDescPtr sdp;
} CommonTitleData, PNTR CommonTitlePtr;


static CommonTitlePtr CommonTitleNew (BioseqPtr bsp, SeqDescPtr sdp)
{
  CommonTitlePtr c = (CommonTitlePtr) MemNew (sizeof (CommonTitleData));
  c->bsp = bsp;
  c->sdp = sdp;
  return c;
}


static CommonTitlePtr CommonTitleFree (CommonTitlePtr c)
{
  if (c != NULL) {
    c = MemFree (c);
  }
  return c;
}


static ValNodePtr CommonTitleListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = CommonTitleFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static void RemoveCommonTitles (ValNodePtr vnp, CharPtr common_title)
{
  CommonTitlePtr c;
  ObjValNodePtr  ovp;

  while (vnp != NULL) {
    c = vnp->data.ptrvalue;
    if (StringCmp (c->sdp->data.ptrvalue, common_title) == 0 && c->sdp->extended > 0) {
      ovp = (ObjValNodePtr) c->sdp;
      ovp->idx.deleteme = TRUE;
    }
    vnp = vnp->next;
  }
}


static int LIBCALLBACK SortCommonTitle (VoidPtr ptr1, VoidPtr ptr2)

{
  CommonTitlePtr c1;
  CommonTitlePtr c2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      c1 = (CommonTitlePtr) vnp1->data.ptrvalue;
      c2 = (CommonTitlePtr) vnp2->data.ptrvalue;
      if (c1 != NULL && c2 != NULL && c1->sdp != NULL && c2->sdp != NULL
          && c1->sdp->data.ptrvalue != NULL && c2->sdp->data.ptrvalue != NULL) {
        return StringCmp (c1->sdp->data.ptrvalue, c2->sdp->data.ptrvalue);
      }
    }
  }
  return 0;
}


static void CollectCommonTitle (BioseqPtr bsp, Pointer data)
{
  SeqDescPtr sdp;

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL) {
    return;
  }

  sdp = bsp->descr;
  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_title) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, 0, CommonTitleNew (bsp, sdp));
    }
    sdp = sdp->next;
  }
}


static CharPtr FindCommonTitleFromList (ValNodePtr list)
{
  ValNodePtr vnp;
  CommonTitlePtr c;
  Int4 num_common = 0, num_total, num_expected;
  CharPtr common_title;

  if (list == NULL) {
    return NULL;
  }
  num_total = ValNodeLen (list);
  if (num_total % 2 != 0 || num_total < 4) {
    return NULL;
  }
  num_expected = num_total / 2;

  c = list->data.ptrvalue;
  common_title = c->sdp->data.ptrvalue;
  num_common = 1;

  for (vnp = list->next; vnp != NULL; vnp = vnp->next) {
    c = (CommonTitlePtr) vnp->data.ptrvalue;
    if (StringCmp (common_title, c->sdp->data.ptrvalue) == 0) {
      num_common++;
    } else {
      num_common = 1;
      common_title = c->sdp->data.ptrvalue;
    }
  }
  if (num_common == num_expected) {
    return StringSave (common_title);
  }

  return NULL;
}


static void PromoteCommonTitlesSetCallback (BioseqSetPtr bssp, Pointer data)
{
  ValNodePtr list = NULL;
  CharPtr common_title = NULL;
  SeqDescrPtr sdp;
  Int4        num_member = 0;
  SeqEntryPtr s;
  CharPtr     set_title = NULL;

  if (bssp == NULL || !GetsDocsumTitle (bssp->_class)) {
    return;
  }

  VisitBioseqsInSet (bssp, &list, CollectCommonTitle);
  list = ValNodeSort (list, SortCommonTitle);

  common_title = FindCommonTitleFromList (list);
  if (common_title != NULL) {
    s = bssp->seq_set;
    while (s != NULL) {
      num_member++;
      s = s->next;
    }
    if (ValNodeLen (list) == num_member) {
      for (sdp = bssp->descr; sdp != NULL && set_title == NULL; sdp = sdp->next) {
        if (sdp->choice == Seq_descr_title) {
          set_title = sdp->data.ptrvalue;
        }
      }
      if (set_title != NULL
          && StringCmp (set_title, common_title) != 0) {
        /* don't remove, the seq titles just happen to be identical */
        common_title = MemFree (common_title);
      }
    }
  }
  if (common_title != NULL) {
    sdp = SeqDescrNew (NULL);
    sdp->choice = Seq_descr_title;
    sdp->data.ptrvalue = common_title;
    sdp->next = bssp->descr;
    bssp->descr = sdp;
    RemoveCommonTitles (list, common_title);
  }
  list = CommonTitleListFree(list);
}


NLM_EXTERN void PromoteCommonTitlesToSet (SeqEntryPtr sep)
{
  VisitSetsInSep (sep, NULL, PromoteCommonTitlesSetCallback);
}


NLM_EXTERN void DeleteMultipleTitles (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Boolean       hastitle;
  ValNodePtr    nextsdp;
  Pointer PNTR  prevsdp;
  ValNodePtr    sdp;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
    prevsdp = (Pointer PNTR) &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
  } else return;
  hastitle = FALSE;
  while (sdp != NULL) {
    nextsdp = sdp->next;
    if (sdp->choice == Seq_descr_title) {
      if (hastitle) {
        *(prevsdp) = sdp->next;
        sdp->next = NULL;
        SeqDescFree (sdp);
      } else {
        hastitle = TRUE;
        prevsdp = (Pointer PNTR) &(sdp->next);
      }
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
}

NLM_EXTERN void RenormalizeNucProtSets (SeqEntryPtr sep, Boolean relink)

{
  SeqAnnotPtr    annot;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ValNodePtr     descr;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;
  SeqAnnotPtr    sap, tmp_sap;
  SeqEntryPtr    seqentry;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (bssp->_class >= 13 && bssp->_class <= 16) ||
                         bssp->_class == BioseqseqSet_class_wgs_set ||
                         bssp->_class == BioseqseqSet_class_gen_prod_set ||
                         bssp->_class == BioseqseqSet_class_small_genome_set)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        RenormalizeNucProtSets (sep, relink);
      }
      return;
    }
    if (bssp != NULL && bssp->_class == 1) {
      seqentry = bssp->seq_set;
      if (seqentry != NULL && seqentry->next == NULL) {

        if (relink) {
          SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
          GetSeqEntryParent (sep, &parentptr, &parenttype);
        }

        descr = bssp->descr;
        bssp->descr = NULL;
        annot = bssp->annot;
        bssp->annot = NULL;

        sep->choice = seqentry->choice;
        sep->data.ptrvalue = seqentry->data.ptrvalue;
        seqentry->data.ptrvalue = NULL;
        bssp->seq_set = NULL;
        bssp->seqentry = NULL;
        MemFree (seqentry);
        BioseqSetFree (bssp);

        sap = NULL;
        if (IS_Bioseq (sep)) {
          bsp = (BioseqPtr) sep->data.ptrvalue;
          ValNodeLink (&(bsp->descr), descr);
          if (bsp->annot == NULL) {
            bsp->annot = annot;
            annot = NULL;
          } else {
            sap = bsp->annot;
          }
        } else if (IS_Bioseq_set (sep)) {
          bssp = (BioseqSetPtr) sep->data.ptrvalue;
          ValNodeLink (&(bssp->descr), descr);
          if (bssp->annot == NULL) {
            bssp->annot = annot;
            annot = NULL;
          } else {
            sap = bssp->annot;
          }
        }
        if (sap != NULL) {
          tmp_sap = sap;
          while (tmp_sap->next != NULL) {
            tmp_sap = tmp_sap->next;
          }
          tmp_sap->next = annot;
          MergeAdjacentAnnotsInList (sap);
        }

        DeleteMultipleTitles (sep, NULL, 0, 0);

        if (relink) {
          SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
          RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
        }
      }
    }
  }
}

NLM_EXTERN ValNodePtr ExtractBioSourceAndPubs (SeqEntryPtr sep)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ValNodePtr    descr;
  ValNodePtr    last;
  ValNodePtr    nextsdp;
  Pointer PNTR  prevsdp;
  ValNodePtr    sdp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return NULL;
  descr = NULL;
  last = NULL;
  sdp = NULL;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
    prevsdp = (Pointer PNTR) &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
  } else return NULL;
  while (sdp != NULL) {
    nextsdp = sdp->next;
    if (sdp->choice == Seq_descr_pub || sdp->choice == Seq_descr_source) {
      *(prevsdp) = sdp->next;
      sdp->next = NULL;
      if (descr == NULL) {
        descr = sdp;
        last = descr;
      } else if (last != NULL) {
        last->next = sdp;
        last = last->next;
      }
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
  return descr;
}

NLM_EXTERN void ReplaceBioSourceAndPubs (SeqEntryPtr sep, ValNodePtr descr)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ValNodePtr    last;
  Pointer PNTR  prevsdp;
  ValNodePtr    sdp;

  if (sep == NULL || descr == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
    prevsdp = (Pointer PNTR) &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
  } else return;
  last = descr;
  while (last->next != NULL) {
    last = last->next;
  }
  last->next = sdp;
  *(prevsdp) = descr;
}

typedef struct targetdata {
  BioseqPtr    bsp;
  SeqEntryPtr  nps;
  Boolean      skipGenProdSet;
} TargetData, PNTR TargetDataPtr;

static Boolean ReturnStackToItem (GatherContextPtr gcp)

{
  BioseqSetPtr   bssp;
  Int2           i;
  Uint2          itemtype;
  TargetDataPtr  tdp;

  if (gcp == NULL) return TRUE;
  tdp = (TargetDataPtr) gcp->userdata;
  if (tdp == NULL) return TRUE;
  if (gcp->gatherstack != NULL && gcp->numstack > 0) {
    for (i = 0; i < gcp->numstack; i++) {
      itemtype = gcp->gatherstack [i].itemtype;
      if (itemtype == OBJ_BIOSEQ || itemtype == OBJ_BIOSEQSET) {
        tdp->nps = SeqMgrGetSeqEntryForData (gcp->gatherstack [i].thisitem);
        if (gcp->gatherstack [i].itemtype == OBJ_BIOSEQSET) {
          bssp = (BioseqSetPtr) gcp->gatherstack [i].thisitem;
          if (bssp->_class != BioseqseqSet_class_genbank &&
              bssp->_class != BioseqseqSet_class_mut_set &&
              bssp->_class != BioseqseqSet_class_pop_set &&
              bssp->_class != BioseqseqSet_class_phy_set &&
              bssp->_class != BioseqseqSet_class_eco_set &&
              bssp->_class != BioseqseqSet_class_wgs_set &&
              bssp->_class != BioseqseqSet_class_small_genome_set &&
              (bssp->_class != BioseqseqSet_class_gen_prod_set ||
           (! tdp->skipGenProdSet))) {
            return FALSE;
          }
        } else if (gcp->gatherstack [i].itemtype == OBJ_BIOSEQ) {
          return FALSE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean GetStackToTarget (GatherContextPtr gcp)

{
  TargetDataPtr  tdp;

  if (gcp == NULL) return TRUE;
  tdp = (TargetDataPtr) gcp->userdata;
  if (tdp == NULL) return TRUE;
  if (gcp->thistype == OBJ_BIOSEQ) {
    if (tdp->bsp == (BioseqPtr) gcp->thisitem) {
      return ReturnStackToItem (gcp);
    }
  }
  return TRUE;
}

NLM_EXTERN SeqEntryPtr LIBCALL GetBestTopParentForDataEx (Uint2 entityID, BioseqPtr bsp, Boolean skipGenProdSet)

{
  BioseqSetPtr  bssp;
  BioseqSetPtr  parent;
  GatherScope   gs;
  TargetData    td;

  td.bsp = bsp;
  td.nps = NULL;
  td.skipGenProdSet = skipGenProdSet;
  if (entityID > 0 && bsp != NULL) {
    if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
      bssp = (BioseqSetPtr) bsp->idx.parentptr;
      if (bssp != NULL && bssp->_class == BioseqseqSet_class_parts && bssp->idx.parenttype == OBJ_BIOSEQSET) {
        parent = (BioseqSetPtr) bssp->idx.parentptr;
        if (parent != NULL && parent->_class == BioseqseqSet_class_segset) {
          bssp = parent;
        }
      }
      if (bssp != NULL && bssp->_class == BioseqseqSet_class_segset && bssp->idx.parenttype == OBJ_BIOSEQSET) {
        parent = (BioseqSetPtr) bssp->idx.parentptr;
        if (parent != NULL && parent->_class == BioseqseqSet_class_nuc_prot) {
          bssp = parent;
        }
      }
      if (bssp != NULL && bssp->seqentry != NULL) {
        if (bssp->_class == BioseqseqSet_class_nuc_prot ||
            bssp->_class == BioseqseqSet_class_segset ||
            bssp->_class == BioseqseqSet_class_parts) {
          return bssp->seqentry;
        }
      }
      if (bsp->seqentry != NULL) {
        return bsp->seqentry;
      }
    }
    MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
    gs.seglevels = 1;
    MemSet ((Pointer) (gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
    gs.ignore[OBJ_BIOSEQ] = FALSE;
    gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
    GatherEntity (entityID, (Pointer) &td, GetStackToTarget, &gs);
  }
  return td.nps;
}

NLM_EXTERN SeqEntryPtr LIBCALL GetBestTopParentForData (Uint2 entityID, BioseqPtr bsp)

{
  return GetBestTopParentForDataEx (entityID, bsp, FALSE);
}

NLM_EXTERN SeqEntryPtr LIBCALL GetBestTopParentForItemIDEx (Uint2 entityID, Uint4 itemID, Uint2 itemtype, Boolean skipGenProdSet)

{
  TargetData  td;

  td.bsp = NULL;
  td.nps = NULL;
  td.skipGenProdSet = skipGenProdSet;
  if (entityID > 0 && itemID > 0 && itemtype > 0) {
    GatherItem (entityID, itemID, itemtype, (Pointer) &td, ReturnStackToItem);
  }
  return td.nps;
}

NLM_EXTERN SeqEntryPtr LIBCALL GetBestTopParentForItemID (Uint2 entityID, Uint4 itemID, Uint2 itemtype)

{
  return GetBestTopParentForItemIDEx (entityID, itemID, itemtype, FALSE);
}

NLM_EXTERN SeqEntryPtr LIBCALL GetTopSeqEntryForEntityID (Uint2 entityID)

{
  ObjMgrDataPtr  omdp;
  SeqSubmitPtr   ssp;

  omdp = ObjMgrGetData (entityID);
  if (omdp != NULL) {
    switch (omdp->datatype) {
      case OBJ_SEQSUB :
        ssp = (SeqSubmitPtr) omdp->dataptr;
        if (ssp != NULL && ssp->datatype == 1) {
          return (SeqEntryPtr) ssp->data;
        }
        break;
      case OBJ_BIOSEQ :
        return (SeqEntryPtr) omdp->choice;
      case OBJ_BIOSEQSET :
        return (SeqEntryPtr) omdp->choice;
      default :
        break;
    }
  }
  return NULL;
}

NLM_EXTERN Boolean CheckSeqLocForPartialEx (SeqLocPtr location, BoolPtr p5ptr, BoolPtr p3ptr, Int4Ptr limptr)

{
  SeqLocPtr   firstSlp;
  IntFuzzPtr  ifp;
  SeqLocPtr   lastSlp;
  Int4        lim;
  Boolean     partial5;
  Boolean     partial3;
  SeqIntPtr   sip;
  SeqLocPtr   slp;
  SeqPntPtr   spp;

  partial5 = FALSE;
  partial3 = FALSE;
  lim = -1;
  if (location != NULL) {
    firstSlp = NULL;
    lastSlp = NULL;
    slp = SeqLocFindNext (location, NULL);
    while (slp != NULL) {
      if (firstSlp == NULL) {
        firstSlp = slp;
      }
      lastSlp = slp;
      slp = SeqLocFindNext (location, slp);
    }
    if (firstSlp != NULL) {
      if (firstSlp->choice == SEQLOC_INT && firstSlp->data.ptrvalue != NULL) {
        sip = (SeqIntPtr) firstSlp->data.ptrvalue;
        if (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev) {
          ifp = sip->if_to;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 1) {
            partial5 = TRUE;
          }
        } else {
          ifp = sip->if_from;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 2) {
            partial5 = TRUE;
          }
        }
      } else if (firstSlp->choice == SEQLOC_PNT && firstSlp->data.ptrvalue != NULL) {
        spp = (SeqPntPtr) firstSlp->data.ptrvalue;
        if (spp->strand == Seq_strand_minus || spp->strand == Seq_strand_both_rev) {
          ifp = spp->fuzz;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 1) {
            partial5 = TRUE;
          }
        } else {
          ifp = spp->fuzz;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 2) {
            partial5 = TRUE;
          }
        }
        ifp = spp->fuzz;
        if (ifp != NULL && ifp->choice == 4) {
          lim = ifp->a;
        }
      }
    }
    if (lastSlp != NULL) {
      if (lastSlp->choice == SEQLOC_INT && lastSlp->data.ptrvalue != NULL) {
        sip = (SeqIntPtr) lastSlp->data.ptrvalue;
        if (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev) {
          ifp = sip->if_from;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 2) {
            partial3 = TRUE;
          }
        } else {
          ifp = sip->if_to;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 1) {
            partial3 = TRUE;
          }
        }
      } else if (lastSlp->choice == SEQLOC_PNT && lastSlp->data.ptrvalue != NULL) {
        spp = (SeqPntPtr) lastSlp->data.ptrvalue;
        if (spp->strand == Seq_strand_minus || spp->strand == Seq_strand_both_rev) {
          ifp = spp->fuzz;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 2) {
            partial3 = TRUE;
          }
        } else {
          ifp = spp->fuzz;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 1) {
            partial3 = TRUE;
          }
        }
        ifp = spp->fuzz;
        if (ifp != NULL && ifp->choice == 4) {
          lim = ifp->a;
        }
      }
    }
  }
  if (p5ptr != NULL) {
    *p5ptr = partial5;
  }
  if (p3ptr != NULL) {
    *p3ptr = partial3;
  }
  if (limptr != NULL) {
    *limptr = lim;
  }
  return (Boolean) (partial5 || partial3 || lim == 3 || lim == 4);
}

NLM_EXTERN Boolean CheckSeqLocForPartial (SeqLocPtr location, BoolPtr p5ptr, BoolPtr p3ptr)

{
  return CheckSeqLocForPartialEx (location, p5ptr, p3ptr, NULL);
}

static void ConvertWholeToIntLoc (SeqLocPtr slp)
{
  BioseqPtr bsp;
  SeqIntPtr sip;
  
  if (slp == NULL || slp->choice != SEQLOC_WHOLE || slp->data.ptrvalue == NULL)
  {
    return;
  }
  bsp = BioseqFind (slp->data.ptrvalue);
  if (bsp == NULL)
  {
    return;
  }
  
  sip = SeqIntNew ();
  if (sip != NULL)
  {
    sip->from = 0;
    sip->to = bsp->length - 1;
    sip->id = SeqIdDup (SeqIdFindBest (bsp->id, 0));
    sip->strand = bsp->strand;
    slp->data.ptrvalue = SeqIdFree (slp->data.ptrvalue);
    slp->data.ptrvalue = sip;
    slp->choice = SEQLOC_INT;
  } 
}

NLM_EXTERN void SetSeqLocPartialEx (SeqLocPtr location, Boolean partial5, Boolean partial3, Int4 lim)

{
  SeqLocPtr   firstSlp;
  IntFuzzPtr  ifp;
  SeqLocPtr   lastSlp;
  SeqIntPtr   sip;
  SeqLocPtr   slp;
  SeqPntPtr   spp;

  if (location != NULL) {
    /* if whole, need to convert to int */
    if (partial5 || partial3)
    {
      ConvertWholeToIntLoc (location);
    }
  
    firstSlp = NULL;
    lastSlp = NULL;
    slp = SeqLocFindNext (location, NULL);
    while (slp != NULL) {
      if (firstSlp == NULL) {
        firstSlp = slp;
      }
      lastSlp = slp;
      slp = SeqLocFindNext (location, slp);
    }
    if (firstSlp != NULL) {
      if (firstSlp->choice == SEQLOC_INT && firstSlp->data.ptrvalue != NULL) {
        sip = (SeqIntPtr) firstSlp->data.ptrvalue;
        if (partial5) {
          ifp = IntFuzzNew ();
          if (ifp != NULL) {
            ifp->choice = 4;
            if (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev) {
              sip->if_to = IntFuzzFree (sip->if_to);
              sip->if_to = ifp;
              ifp->a = 1;
            } else {
              sip->if_from = IntFuzzFree (sip->if_from);
              sip->if_from = ifp;
              ifp->a = 2;
            }
          }
        } else {
          if (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev) {
            sip->if_to = IntFuzzFree (sip->if_to);
          } else {
            sip->if_from = IntFuzzFree (sip->if_from);
          }
        }
      } else if (firstSlp->choice == SEQLOC_PNT && firstSlp->data.ptrvalue != NULL) {
        spp = (SeqPntPtr) firstSlp->data.ptrvalue;
        if (partial5) {
          ifp = IntFuzzNew ();
          if (ifp != NULL) {
            ifp->choice = 4;
            if (spp->strand == Seq_strand_minus || spp->strand == Seq_strand_both_rev) {
              spp->fuzz = IntFuzzFree (spp->fuzz);
              spp->fuzz = ifp;
              ifp->a = 1;
            } else {
              spp->fuzz = IntFuzzFree (spp->fuzz);
              spp->fuzz = ifp;
              ifp->a = 2;
            }
          }
        } else if (lim == 3 || lim == 4) {
          ifp = IntFuzzNew ();
          if (ifp != NULL) {
            ifp->choice = 4;
            spp->fuzz = IntFuzzFree (spp->fuzz);
            spp->fuzz = ifp;
            ifp->a = lim;
          }
        } else {
          if (spp->strand == Seq_strand_minus || spp->strand == Seq_strand_both_rev) {
            spp->fuzz = IntFuzzFree (spp->fuzz);
          } else {
            spp->fuzz = IntFuzzFree (spp->fuzz);
          }
        }
      }
    }
    if (lastSlp != NULL) {
      if (lastSlp->choice == SEQLOC_INT && lastSlp->data.ptrvalue != NULL) {
        sip = (SeqIntPtr) lastSlp->data.ptrvalue;
        if (partial3) {
          ifp = IntFuzzNew ();
          if (ifp != NULL) {
            ifp->choice = 4;
            if (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev) {
              sip->if_from = IntFuzzFree (sip->if_from);
              sip->if_from = ifp;
              ifp->a = 2;
            } else {
              sip->if_to = IntFuzzFree (sip->if_to);
              sip->if_to = ifp;
              ifp->a = 1;
            }
          }
        } else {
          if (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev) {
            sip->if_from = IntFuzzFree (sip->if_from);
          } else {
            sip->if_to = IntFuzzFree (sip->if_to);
          }
        }
      } else if (lastSlp->choice == SEQLOC_PNT && lastSlp->data.ptrvalue != NULL) {
        spp = (SeqPntPtr) lastSlp->data.ptrvalue;
        if (partial3) {
          ifp = IntFuzzNew ();
          if (ifp != NULL) {
            ifp->choice = 4;
            if (spp->strand == Seq_strand_minus || spp->strand == Seq_strand_both_rev) {
              spp->fuzz = IntFuzzFree (spp->fuzz);
              spp->fuzz = ifp;
              ifp->a = 2;
            } else {
              spp->fuzz = IntFuzzFree (spp->fuzz);
              spp->fuzz = ifp;
              ifp->a = 1;
            }
          }
        } else if (lim == 3 || lim == 4) {
          ifp = IntFuzzNew ();
          if (ifp != NULL) {
            ifp->choice = 4;
            spp->fuzz = IntFuzzFree (spp->fuzz);
            spp->fuzz = ifp;
            ifp->a = lim;
          }
        } else {
          if (spp->strand == Seq_strand_minus || spp->strand == Seq_strand_both_rev) {
            spp->fuzz = IntFuzzFree (spp->fuzz);
          } else {
            spp->fuzz = IntFuzzFree (spp->fuzz);
          }
        }
      }
    }
  }
}

NLM_EXTERN void SetSeqLocPartial (SeqLocPtr location, Boolean partial5, Boolean partial3)

{
  SetSeqLocPartialEx (location, partial5, partial3, -1);
}

NLM_EXTERN ValNodePtr GetSeqLocPartialSet (SeqLocPtr location)

{
  ValNodePtr  head = NULL, last = NULL, vnp;
  Int4        lim;
  Boolean     noLeft;
  Boolean     noRight;
  SeqLocPtr   slp;
  Int4        val;

  if (location == NULL) return NULL;

  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL) {
    CheckSeqLocForPartialEx (slp, &noLeft, &noRight, &lim);
    val = 0;
    if (noLeft) {
      val |= 2;
    }
    if (noRight) {
      val |= 1;
    }
    if (lim == 3) {
      val |= 4;
    } else if (lim == 4) {
      val |= 8;
    }
    vnp = ValNodeAddInt (&last, 0, val);
    if (head == NULL) {
      head = vnp;
    }
    last = vnp;
    slp = SeqLocFindNext (location, slp);
  }

  return head;
}

NLM_EXTERN void SetSeqLocPartialSet (SeqLocPtr location, ValNodePtr vnp)

{
  Int4        lim;
  Boolean     noLeft;
  Boolean     noRight;
  SeqLocPtr   slp;
  Int4        val;

  if (location == NULL || vnp == NULL) return;

  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL && vnp != NULL) {
    val = (Int4) vnp->data.intvalue;
    noLeft = (Boolean) ((val & 2) != 0);
    noRight = (Boolean) ((val & 1) != 0);
    lim = -1;
    if ((val & 4) != 0) {
      lim = 3;
    } else if ((val & 8) != 0) {
      lim = 4;
    }
    SetSeqLocPartialEx (slp, noLeft, noRight, lim);
    slp = SeqLocFindNext (location, slp);
    vnp = vnp->next;
  }
}

/* KeyTag section */

NLM_EXTERN int LIBCALLBACK SortVnpByString (VoidPtr ptr1, VoidPtr ptr2)

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
      }
    }
  }
  return 0;
}

NLM_EXTERN int LIBCALLBACK SortVnpByStringCS (VoidPtr ptr1, VoidPtr ptr2)

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
        return StringCmp (str1, str2);
      }
    }
  }
  return 0;
}

NLM_EXTERN int LIBCALLBACK SortVnpByStringCI (VoidPtr ptr1, VoidPtr ptr2)

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
        return StringCmp (str1, str2);
      }
    }
  }
  return 0;
}

NLM_EXTERN int LIBCALLBACK SortVnpByStringCIUCFirst (VoidPtr ptr1, VoidPtr ptr2)

{
  int         comp;
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
        comp = StringICmp (str1, str2);
        if (comp != 0) return comp;
        return StringCmp (str1, str2);
      }
    }
  }
  return 0;
}

NLM_EXTERN int LIBCALLBACK SortVnpByStringCILCFirst (VoidPtr ptr1, VoidPtr ptr2)

{
  int         comp;
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
        comp = StringICmp (str1, str2);
        if (comp != 0) return comp;
        return StringCmp (str2, str1);
      }
    }
  }
  return 0;
}

NLM_EXTERN int LIBCALLBACK SortVnpByNaturalCS (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1, str2;
  ValNodePtr  vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;

  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;

  str1 = (CharPtr) vnp1->data.ptrvalue;
  str2 = (CharPtr) vnp2->data.ptrvalue;
  if (str1 == NULL || str2 == NULL) return 0;

  return NaturalStringCmp (str1, str2);
}

NLM_EXTERN int LIBCALLBACK SortVnpByNaturalCI (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1, str2;
  ValNodePtr  vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;

  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;

  str1 = (CharPtr) vnp1->data.ptrvalue;
  str2 = (CharPtr) vnp2->data.ptrvalue;
  if (str1 == NULL || str2 == NULL) return 0;

  return NaturalStringICmp (str1, str2);
}

NLM_EXTERN ValNodePtr UniqueValNode (ValNodePtr list)

{
  CharPtr       last;
  ValNodePtr    next;
  Pointer PNTR  prev;
  CharPtr       str;
  ValNodePtr    vnp;

  if (list == NULL) return NULL;
  last = (CharPtr) list->data.ptrvalue;
  vnp = list->next;
  prev = (Pointer PNTR) &(list->next);
  while (vnp != NULL) {
    next = vnp->next;
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringICmp (last, str) == 0) {
      vnp->next = NULL;
      *prev = next;
      ValNodeFreeData (vnp);
    } else {
      last = (CharPtr) vnp->data.ptrvalue;
      prev = (Pointer PNTR) &(vnp->next);
    }
    vnp = next;
  }

  return list;
}

NLM_EXTERN ValNodePtr UniqueStringValNodeCS (ValNodePtr list)

{
  CharPtr       last;
  ValNodePtr    next;
  Pointer PNTR  prev;
  CharPtr       str;
  ValNodePtr    vnp;

  if (list == NULL) return NULL;
  last = (CharPtr) list->data.ptrvalue;
  vnp = list->next;
  prev = (Pointer PNTR) &(list->next);
  while (vnp != NULL) {
    next = vnp->next;
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringCmp (last, str) == 0) {
      vnp->next = NULL;
      *prev = next;
      ValNodeFreeData (vnp);
    } else {
      last = (CharPtr) vnp->data.ptrvalue;
      prev = (Pointer PNTR) &(vnp->next);
    }
    vnp = next;
  }

  return list;
}

NLM_EXTERN ValNodePtr UniqueStringValNodeCI (ValNodePtr list)

{
  CharPtr       last;
  ValNodePtr    next;
  Pointer PNTR  prev;
  CharPtr       str;
  ValNodePtr    vnp;

  if (list == NULL) return NULL;
  last = (CharPtr) list->data.ptrvalue;
  vnp = list->next;
  prev = (Pointer PNTR) &(list->next);
  while (vnp != NULL) {
    next = vnp->next;
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringICmp (last, str) == 0) {
      vnp->next = NULL;
      *prev = next;
      ValNodeFreeData (vnp);
    } else {
      last = (CharPtr) vnp->data.ptrvalue;
      prev = (Pointer PNTR) &(vnp->next);
    }
    vnp = next;
  }

  return list;
}

NLM_EXTERN int LIBCALLBACK SortByIntvalue (VoidPtr ptr1, VoidPtr ptr2)

{
  Int4        val1;
  Int4        val2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  val1 = (Int4) vnp1->data.intvalue;
  val2 = (Int4) vnp2->data.intvalue;
  if (val1 > val2) {
    return 1;
  } else if (val1 < val2) {
    return -1;
  }
  return 0;
}

NLM_EXTERN ValNodePtr UniqueIntValNode (ValNodePtr list)

{
  Int4          curr, last;
  ValNodePtr    next;
  Pointer PNTR  prev;
  ValNodePtr    vnp;

  if (list == NULL) return NULL;
  last = (Int4) list->data.intvalue;
  vnp = list->next;
  prev = (Pointer PNTR) &(list->next);
  while (vnp != NULL) {
    next = vnp->next;
    curr = (Int4) vnp->data.intvalue;
    if (last == curr) {
      vnp->next = NULL;
      *prev = next;
      ValNodeFree (vnp);
    } else {
      last = (Int4) vnp->data.intvalue;
      prev = (Pointer PNTR) &(vnp->next);
    }
    vnp = next;
  }

  return list;
}

NLM_EXTERN int LIBCALLBACK SortByPtrvalue (VoidPtr ptr1, VoidPtr ptr2)

{
  VoidPtr     val1;
  VoidPtr     val2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  val1 = (VoidPtr) vnp1->data.ptrvalue;
  val2 = (VoidPtr) vnp2->data.ptrvalue;
  if (val1 > val2) {
    return 1;
  } else if (val1 < val2) {
    return -1;
  }
  return 0;
}

NLM_EXTERN ValNodePtr UniquePtrValNode (ValNodePtr list)

{
  VoidPtr       curr, last;
  ValNodePtr    next;
  Pointer PNTR  prev;
  ValNodePtr    vnp;

  if (list == NULL) return NULL;
  last = (VoidPtr) list->data.ptrvalue;
  vnp = list->next;
  prev = (Pointer PNTR) &(list->next);
  while (vnp != NULL) {
    next = vnp->next;
    curr = (VoidPtr) vnp->data.ptrvalue;
    if (last == curr) {
      vnp->next = NULL;
      *prev = next;
      ValNodeFree (vnp);
    } else {
      last = (VoidPtr) vnp->data.ptrvalue;
      prev = (Pointer PNTR) &(vnp->next);
    }
    vnp = next;
  }

  return list;
}

NLM_EXTERN void KeyTagInit (KeyTag PNTR ktp, ValNodePtr list)

{
  Int2          i;
  CharPtr PNTR  index;
  Int2          num;
  ValNodePtr    vnp;

  if (ktp == NULL || list == NULL) return;
  list = ValNodeSort (list, SortVnpByString);
  list = UniqueValNode (list);
  num = ValNodeLen (list);
  index = MemNew (sizeof (CharPtr) * (num + 1));

  for (vnp = list, i = 0; vnp != NULL && i < num; vnp = vnp->next, i++) {
    index [i] = (CharPtr) vnp->data.ptrvalue;
  }

  ktp->num = num;
  ktp->list = list;
  ktp->index = index;
}

NLM_EXTERN void KeyTagClear (KeyTag PNTR ktp)

{
  if (ktp == NULL) return;
  ktp->num = 0;
  ktp->list = ValNodeFreeData (ktp->list);
  ktp->index = MemFree (ktp->index);
}

NLM_EXTERN Int2 KeyFromTag (KeyTag PNTR ktp, CharPtr tag)

{
  Int2  L, R, mid, compare;

  if (ktp == NULL || ktp->list == NULL || ktp->index == NULL) return 0;
  if (tag == NULL) return 0;

  L = 0;
  R = ktp->num - 1;
  while (L < R) {
    mid = (L + R) / 2;
    compare = StringICmp (ktp->index [mid], tag);
    if (compare < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }
  if (StringICmp (ktp->index [R], tag) == 0) {
    return (R + 1);
  }

  return 0;
}

NLM_EXTERN CharPtr TagFromKey (KeyTag PNTR ktp, Int2 key)

{
  if (ktp == NULL || ktp->list == NULL || ktp->index == NULL) return 0;
  if (key < 1 || key > ktp->num) return 0;
  key--;
  return ktp->index [key];
}

/* begin PromoteXrefs section */

typedef struct geneextendlist {
  GeneRefPtr  grp;
  SeqLocPtr   slp;
  ObjMgrPtr   omp;
  Boolean     rsult;
  Char        label [41];
} GeneExtendList, PNTR GeneExtendPtr;

static Boolean GeneExtendFunc (GatherContextPtr gcp)

{
  BioseqPtr      bsp;
  GeneExtendPtr  gep;
  GeneRefPtr     grp;
  Boolean        hasNulls;
  ObjMgrTypePtr  omtp;
  SeqFeatPtr     sfp;
  SeqLocPtr      slp;
  Char           thislabel [41];

  if (gcp == NULL) return TRUE;

  gep = (GeneExtendPtr) gcp->userdata;
  if (gep == NULL ) return TRUE;

  thislabel [0] = '\0';

  if (gcp->thistype == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE && sfp->data.value.ptrvalue != NULL) {
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      omtp = ObjMgrTypeFind (gep->omp, gcp->thistype, NULL, NULL);
      if (omtp == NULL) {
        return TRUE;
      }
      if (omtp->labelfunc != NULL) {
        (*(omtp->labelfunc)) (gcp->thisitem, thislabel, 40, OM_LABEL_CONTENT);
      }
      if (thislabel [0] != '\0') {
        if (StringICmp (thislabel, gep->label) == 0) {
          if (SeqLocCompare (gep->slp, sfp->location) != SLC_NO_MATCH) {
            bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
            if (bsp != NULL) {
              slp = SeqLocMerge (bsp, sfp->location, gep->slp, TRUE, FALSE, FALSE);
              if (slp != NULL) {
                sfp->location = SeqLocFree (sfp->location);
                sfp->location = slp;
                if (bsp->repr == Seq_repr_seg) {
                  slp = SegLocToPartsEx (bsp, sfp->location, TRUE);
                  sfp->location = SeqLocFree (sfp->location);
                  sfp->location = slp;
                  hasNulls = LocationHasNullsBetween (sfp->location);
                  sfp->partial = (sfp->partial || hasNulls);
                }
                FreeAllFuzz (slp);
                gep->rsult = TRUE;
              }
            }
          }
          return FALSE;
        }
      }
    }
  }
  return TRUE;
}

/*
static Boolean ExtendGene (GeneRefPtr grp, SeqEntryPtr nsep, SeqLocPtr slp)

{
  GeneExtendList  gel;
  GatherScope     gs;
  ObjMgrTypePtr   omtp;
  SeqFeatPtr      sfp;

  if (grp == NULL || nsep == NULL || slp == NULL) return FALSE;
  gel.grp = grp;
  gel.slp = slp;
  gel.omp = ObjMgrGet ();
  gel.label [0] = '\0';
  gel.rsult = FALSE;
  omtp = ObjMgrTypeFind (gel.omp, OBJ_SEQFEAT, NULL, NULL);
  if (omtp != NULL && omtp->labelfunc != NULL) {
    sfp = SeqFeatNew ();
    if (sfp != NULL) {
      sfp->data.choice = SEQFEAT_GENE;
      sfp->data.value.ptrvalue = (Pointer) grp;
      (*(omtp->labelfunc)) ((Pointer) sfp, gel.label, 40, OM_LABEL_CONTENT);
      sfp->data.value.ptrvalue = NULL;
      SeqFeatFree (sfp);
    }
  }
  MemSet ((Pointer)(&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherSeqEntry (nsep, (Pointer) &gel, GeneExtendFunc, &gs);
  return gel.rsult;
}
*/

NLM_EXTERN void SetEmptyGeneticCodes (SeqAnnotPtr sap, Int2 genCode)

{
  CdRegionPtr     crp;
  GeneticCodePtr  gc;
  SeqFeatPtr      sfp;
  ValNodePtr      vnp;

  if (sap == NULL || sap->type != 1) return;
  for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
    if (sfp->data.choice == SEQFEAT_CDREGION) {
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL) {
        gc = crp->genetic_code;
        if (gc != NULL) {
          vnp = gc->data.ptrvalue;
          if (vnp != NULL && vnp->choice == 2) {
            vnp->data.intvalue = (Int4) genCode;
            /*
            if (vnp->data.intvalue == 0) {
              vnp->data.intvalue = (Int4) genCode;
            }
            */
          }
        }
      }
    }
  }
}

NLM_EXTERN void PromoteXrefsExEx (
  SeqFeatPtr sfp,
  BioseqPtr bsp,
  Uint2 entityID,
  Boolean include_stop,
  Boolean remove_trailingX,
  Boolean gen_prod_set,
  Boolean force_local_id,
  BoolPtr seq_fetch_failP
)

{
  Int2                 adv;
  ByteStorePtr         bs;
  BioseqSetPtr         bssp;
  Char                 ch;
  CharPtr              comment;
  CdRegionPtr          crp;
  Int2                 ctr = 1;
  ValNodePtr           descr;
  SeqFeatPtr           first;
  GBQualPtr            gbq;
  Int4                 i;
  Char                 id [128];
  SeqEntryPtr          last;
  Char                 lcl [128];
  BioseqPtr            mbsp;
  MolInfoPtr           mip;
  SeqEntryPtr          msep;
  SeqFeatXrefPtr       next;
  GBQualPtr            nextqual;
  SeqEntryPtr          old;
  ObjMgrDataPtr        omdptop;
  ObjMgrData           omdata;
  Uint2                parenttype;
  Pointer              parentptr;
  Boolean              partial5;
  Boolean              partial3;
  BioseqPtr            pbsp;
  SeqFeatXrefPtr PNTR  prev;
  GBQualPtr PNTR       prevqual;
  SeqFeatPtr           prot;
  CharPtr              protseq;
  ProtRefPtr           prp, prp2;
  SeqEntryPtr          psep;
  CharPtr              ptr;
  CharPtr              rnaseq;
  SeqEntryPtr          sep;
  SeqHistPtr           shp;
  SeqIdPtr             sip;
  SeqEntryPtr          target = NULL;
  Uint4                version = 0;
  long int             val;
  ValNodePtr           vnp;
  SeqFeatXrefPtr       xref;
  Boolean              ok_to_remove;
  /*
  DbtagPtr             dbt;
  SeqFeatPtr           gene;
  GeneRefPtr           grp;
  */

  if (seq_fetch_failP != NULL) {
    *seq_fetch_failP = FALSE;
  }

  if (sfp == NULL || bsp == NULL) return;

  /* set subtypes, used to find mRNA features for genomic product sets */

  first = sfp;
  while (sfp != NULL) {
    if (sfp->idx.subtype == 0) {
      sfp->idx.subtype = FindFeatDefType (sfp);
    }
    sfp = sfp->next;
  }

  /* no longer expand genes specified by qualifiers on other features (except repeat_region) */

  /*
  sfp = first;
  while (sfp != NULL) {
    prev = &(sfp->xref);
    xref = sfp->xref;
    while (xref != NULL) {
      next = xref->next;
      if (xref->data.choice == SEQFEAT_GENE &&
          sfp->data.choice != SEQFEAT_GENE &&
          sfp->idx.subtype != FEATDEF_repeat_region) {
        grp = (GeneRefPtr) xref->data.value.ptrvalue;
        if (grp != NULL && SeqMgrGeneIsSuppressed (grp)) {
        } else {
          xref->data.value.ptrvalue = NULL;
          if (grp != NULL) {
            sep = SeqMgrGetSeqEntryForData (bsp);
            if (ExtendGene (grp, sep, sfp->location)) {
              GeneRefFree (grp);
            } else {
              gene = CreateNewFeature (sep, NULL, SEQFEAT_GENE, NULL);
              if (gene != NULL) {
                gene->data.value.ptrvalue = (Pointer) grp;
                gene->location = SeqLocFree (gene->location);
                gene->location = AsnIoMemCopy (sfp->location,
                                               (AsnReadFunc) SeqLocAsnRead,
                                               (AsnWriteFunc) SeqLocAsnWrite);
                for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
                  dbt = (DbtagPtr) vnp->data.ptrvalue;
                  if (dbt == NULL) continue;
                  ValNodeAddPointer (&(gene->dbxref), 0, (Pointer) DbtagDup (dbt));
                }
              }
            }
          }
          *(prev) = next;
          xref->next = NULL;
          xref->data.choice = 0;
          SeqFeatXrefFree (xref);
        }
      } else {
        prev = &(xref->next);
      }
      xref = next;
    }
    sfp = sfp->next;
  }
  */

  /* expand mRNA features into cDNA product sequences */

  bssp = NULL;
  sep = NULL;
  last = NULL;
  if (gen_prod_set) {
    sep = GetTopSeqEntryForEntityID (entityID);
    if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp != NULL && bssp->seq_set != NULL) {
        last = bssp->seq_set;
        while (last->next != NULL) {
          last = last->next;
        }
      }
    }
  }

  if (gen_prod_set && sep != NULL && bssp != NULL && last != NULL) {
    target = sep;
    SaveSeqEntryObjMgrData (target, &omdptop, &omdata);
    GetSeqEntryParent (target, &parentptr, &parenttype);
    sfp = first;
    while (sfp != NULL) {
      if (sfp->data.choice == SEQFEAT_RNA &&
          /* sfp->idx.subtype != FEATDEF_tRNA && */
          sfp->product == NULL && (! sfp->pseudo)) {
        gbq = sfp->qual;
        prevqual = (GBQualPtr PNTR) &(sfp->qual);
        id [0] = '\0';
        sip = NULL;
        comment = NULL;
        while (gbq != NULL) {
          nextqual = gbq->next;
          if (StringICmp (gbq->qual, "transcript_id") == 0) {
            if (StringDoesHaveText (id) && StringDoesHaveText (gbq->val)) {
              ErrPostEx (SEV_WARNING, ERR_FEATURE_QualWrongThisFeat,
                         "RNA transcript_id %s replacing %s", gbq->val, id);
            }
            *(prevqual) = gbq->next;
            gbq->next = NULL;
            StringNCpy_0 (id, gbq->val, sizeof (id));
            GBQualFree (gbq);
          } else if (StringICmp (gbq->qual, "comment") == 0 &&
                     StringDoesHaveText (gbq->val)) {
            *(prevqual) = gbq->next;
            gbq->next = NULL;
            comment = StringSave (gbq->val);
            GBQualFree (gbq);
          } else {
            prevqual = (GBQualPtr PNTR) &(gbq->next);
          }
          gbq = nextqual;
        }
        if (! StringHasNoText (id)) {
          if (StringChr (id, '|') != NULL) {
            sip = SeqIdParse (id);
          } else if (force_local_id) {
            sprintf (lcl, "lcl|%s", id);
            sip = SeqIdParse (lcl);
          } else {
            adv = ValidateAccnDotVer (id);
            if (adv == 0 || adv == -5) {
              ptr = StringChr (id, '.');
              if (ptr != NULL) {
                *ptr = '\0';
                ptr++;
                if (sscanf (ptr, "%ld", &val) == 1) {
                  version = (Uint4) val;
                }
              }
              sip = SeqIdFromAccession (id, version, NULL);
            } else {
              sprintf (lcl, "lcl|%s", id);
              sip = SeqIdParse (lcl);
            }
          }
        }
        if (sip != NULL || sfp->idx.subtype == FEATDEF_mRNA) {
          rnaseq = GetSequenceByFeature (sfp);
          if (rnaseq == NULL && seq_fetch_failP != NULL) {
            *seq_fetch_failP = TRUE;
          }
          if (rnaseq != NULL) {
            i = (Int4) StringLen (rnaseq);
            bs = BSNew (i + 2);
            if (bs != NULL) {
              BSWrite (bs, (VoidPtr) rnaseq, (Int4) StringLen (rnaseq));
              mbsp = BioseqNew ();
              if (mbsp != NULL) {
                mbsp->repr = Seq_repr_raw;
                mbsp->mol = Seq_mol_rna;
                mbsp->seq_data_type = Seq_code_iupacna;
                mbsp->seq_data = (SeqDataPtr) bs;
                mbsp->length = BSLen (bs);
                BioseqPack (mbsp);
                bs = NULL;
                /*
                sep = GetTopSeqEntryForEntityID (entityID);
                */
                old = SeqEntrySetScope (sep);
                if (sip != NULL) {
                  mbsp->id = sip;
                } else if (sfp->idx.subtype == FEATDEF_mRNA) {
                  /* actually just making rapid unique ID for mRNA */
                  mbsp->id = MakeNewProteinSeqIdEx (sfp->location, NULL, NULL, &ctr);
                }
                CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
                SeqMgrAddToBioseqIndex (mbsp);
                SeqEntrySetScope (old);
                msep = SeqEntryNew ();
                if (msep != NULL) {
                  msep->choice = 1;
                  msep->data.ptrvalue = (Pointer) mbsp;
                  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) mbsp, msep);
                  mip = MolInfoNew ();
                  if (mip != NULL) {
                    switch (sfp->idx.subtype) {
                      case FEATDEF_preRNA :
                        mip->biomol = MOLECULE_TYPE_PRE_MRNA;
                        break;
                      case FEATDEF_mRNA :
                        mip->biomol = MOLECULE_TYPE_MRNA;
                        break;
                      case FEATDEF_tRNA :
                        mip->biomol = MOLECULE_TYPE_TRNA;
                        break;
                      case FEATDEF_rRNA :
                        mip->biomol = MOLECULE_TYPE_RRNA;
                        break;
                      case FEATDEF_snRNA :
                        mip->biomol = MOLECULE_TYPE_SNRNA;
                        break;
                      case FEATDEF_scRNA :
                        mip->biomol = MOLECULE_TYPE_SCRNA;
                        break;
                      case FEATDEF_otherRNA :
                        mip->biomol = MOLECULE_TYPE_TRANSCRIBED_RNA;
                        break;
                      case FEATDEF_snoRNA :
                        mip->biomol = MOLECULE_TYPE_SNORNA;
                        break;
                      case FEATDEF_ncRNA :
                        mip->biomol = MOLECULE_TYPE_NCRNA;
                        break;
                      case FEATDEF_tmRNA :
                        mip->biomol = MOLECULE_TYPE_TMRNA;
                        break;
                      default :
                        mip->biomol = 0;
                        break;
                    }
                    if (partial5 && partial3) {
                      mip->completeness = 5;
                    } else if (partial5) {
                      mip->completeness = 3;
                    } else if (partial3) {
                      mip->completeness = 4;
                    }
                    vnp = CreateNewDescriptor (msep, Seq_descr_molinfo);
                    if (vnp != NULL) {
                      vnp->data.ptrvalue = (Pointer) mip;
                    }
                  }
                  if (comment != NULL) {
                    vnp = CreateNewDescriptor (msep, Seq_descr_comment);
                    if (vnp != NULL) {
                      vnp->data.ptrvalue = (Pointer) comment;
                    }
                  }
                  /* add mRNA sequence to genomic product set */
                  last->next = msep;
                  last = msep;
                  SetSeqFeatProduct (sfp, mbsp);
                }
              }
            }
            rnaseq = MemFree (rnaseq);
          }
        }
      }
      sfp = sfp->next;
    }
    SeqMgrLinkSeqEntry (target, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (target, omdptop, &omdata);
  }

  /* expand coding region features into protein product sequences */

  last = NULL;
  sfp = first;
  while (sfp != NULL) {
    prev = &(sfp->xref);
    xref = sfp->xref;
    while (xref != NULL) {
      next = xref->next;
      if (xref->data.choice == SEQFEAT_PROT &&
          sfp->data.choice == SEQFEAT_CDREGION &&
          sfp->product == NULL && (! sfp->pseudo)) {
        prp = (ProtRefPtr) xref->data.value.ptrvalue;
        ok_to_remove = TRUE;
        if (prp != NULL) {
          crp = (CdRegionPtr) sfp->data.value.ptrvalue;
          if (crp != NULL) {
/**
            crp->frame = 0;
**/
            bs = ProteinFromCdRegionEx (sfp, include_stop, remove_trailingX);
            if (bs == NULL && seq_fetch_failP != NULL) {
              *seq_fetch_failP = TRUE;
            }
            if (bs != NULL) {
              protseq = BSMerge (bs, NULL);
              bs = BSFree (bs);
              if (protseq != NULL) {
                ptr = protseq;
                ch = *ptr;
                while (ch != '\0') {
                  *ptr = TO_UPPER (ch);
                  ptr++;
                  ch = *ptr;
                }
                i = (Int4) StringLen (protseq);
                if (i > 0 && protseq [i - 1] == '*') {
                  protseq [i - 1] = '\0';
                }
                bs = BSNew (i + 2);
                if (bs != NULL) {
                  ptr = protseq;
                  /*
                  if (protseq [0] == '-') {
                    ptr++;
                  }
                  */
                  BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
                }
                protseq = MemFree (protseq);
              }
              pbsp = BioseqNew ();
              if (pbsp != NULL) {
                pbsp->repr = Seq_repr_raw;
                pbsp->mol = Seq_mol_aa;
                pbsp->seq_data_type = Seq_code_ncbieaa;
                pbsp->seq_data = (SeqDataPtr) bs;
                pbsp->length = BSLen (bs);
                bs = NULL;
                sep = NULL;
                mbsp = NULL;
                if (gen_prod_set) {
                  gbq = sfp->qual;
                  prevqual = (GBQualPtr PNTR) &(sfp->qual);
                  id [0] = '\0';
                  sip = NULL;
                  while (gbq != NULL) {
                    nextqual = gbq->next;
                    if (StringICmp (gbq->qual, "transcript_id") == 0) {
                      if (StringDoesHaveText (id) && StringDoesHaveText (gbq->val)) {
                        ErrPostEx (SEV_WARNING, ERR_FEATURE_QualWrongThisFeat,
                                   "CDS transcript_id %s replacing %s", gbq->val, id);
                      }
                      *(prevqual) = gbq->next;
                      gbq->next = NULL;
                      StringNCpy_0 (id, gbq->val, sizeof (id));
                      GBQualFree (gbq);
                    } else if (StringICmp (gbq->qual, "secondary_accession") == 0) {
                      *(prevqual) = gbq->next;
                      gbq->next = NULL;
                      shp = ParseStringIntoSeqHist (NULL, gbq->val);
                      if (shp != NULL) {
                        pbsp->hist = shp;
                      }
                      GBQualFree (gbq);
                    } else {
                      prevqual = (GBQualPtr PNTR) &(gbq->next);
                    }
                    gbq = nextqual;
                  }
                  if (StringHasNoText (id)) {
                    Message (MSG_POSTERR, "No transcript_id on CDS - unable to create nuc-prot set");
                  } else {
                    if (StringChr (id, '|') != NULL) {
                      sip = SeqIdParse (id);
                    } else if (force_local_id) {
                      sprintf (lcl, "lcl|%s", id);
                      sip = SeqIdParse (lcl);
                    } else {
                      adv = ValidateAccnDotVer (id);
                      if (adv == 0 || adv == -5) {
                        ptr = StringChr (id, '.');
                        if (ptr != NULL) {
                          *ptr = '\0';
                          ptr++;
                          if (sscanf (ptr, "%ld", &val) == 1) {
                            version = (Uint4) val;
                          }
                        }
                        sip = SeqIdFromAccession (id, version, NULL);
                      } else {
                        sprintf (lcl, "lcl|%s", id);
                        sip = SeqIdParse (lcl);
                      }
                    }
                  }
                  mbsp = BioseqFind (sip);
                  SeqIdFree (sip);
                  if (mbsp != NULL) {
                    sep = SeqMgrGetSeqEntryForData (mbsp);
                  /*
                  } else {
                    sep = GetBestTopParentForDataEx (entityID, bsp, TRUE);
                  */
                  }
                } else {
                  sep = GetBestTopParentForData (entityID, bsp);
                }
                if (sep == NULL) {
                  Message (MSG_POSTERR, "No location for nuc-prot set for CDS - unable to create nuc-prot set");
                  pbsp = BioseqFree (pbsp);
                  ok_to_remove = FALSE;
                } else {
                  old = SeqEntrySetScope (sep);
                  gbq = sfp->qual;
                  prevqual = (GBQualPtr PNTR) &(sfp->qual);
                  id [0] = '\0';
                  sip = NULL;
                  while (gbq != NULL) {
                    nextqual = gbq->next;
                    if (StringICmp (gbq->qual, "protein_id") == 0) {
                      if (StringDoesHaveText (id) && StringDoesHaveText (gbq->val)) {
                                ErrPostEx (SEV_WARNING, ERR_FEATURE_QualWrongThisFeat,
                          "CDS protein_id %s replacing %s", gbq->val, id);
                      }
                      *(prevqual) = gbq->next;
                      gbq->next = NULL;
                      StringNCpy_0 (id, gbq->val, sizeof (id));
                      GBQualFree (gbq);
                    } else {
                      prevqual = (GBQualPtr PNTR) &(gbq->next);
                    }
                    gbq = nextqual;
                  }
                  if (! StringHasNoText (id)) {
                    if (StringChr (id, '|') != NULL) {
                      sip = SeqIdParse (id);
                    } else if (force_local_id) {
                      sprintf (lcl, "lcl|%s", id);
                      sip = SeqIdParse (lcl);
                    } else {
                      adv = ValidateAccnDotVer (id);
                      if (adv == 0 || adv == -5) {
                        ptr = StringChr (id, '.');
                        if (ptr != NULL) {
                          *ptr = '\0';
                          ptr++;
                          if (sscanf (ptr, "%ld", &val) == 1) {
                            version = (Uint4) val;
                          }
                        }
                        sip = SeqIdFromAccession (id, version, NULL);
                      } else {
                        sprintf (lcl, "lcl|%s", id);
                        sip = SeqIdParse (lcl);
                      }
                    }
                  }
                  if (sip != NULL) {
                    pbsp->id = sip;
                  } else {
                    pbsp->id = MakeNewProteinSeqIdEx (sfp->location, NULL, NULL, &ctr);
                  }
                  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
                  SeqMgrAddToBioseqIndex (pbsp);
                  SeqEntrySetScope (old);
                  psep = SeqEntryNew ();
                  if (psep != NULL) {
                    psep->choice = 1;
                    psep->data.ptrvalue = (Pointer) pbsp;
                    SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) pbsp, psep);
                    mip = MolInfoNew ();
                    if (mip != NULL) {
                      mip->biomol = 8;
                      mip->tech = 8;
                      if (partial5 && partial3) {
                        mip->completeness = 5;
                      } else if (partial5) {
                        mip->completeness = 3;
                      } else if (partial3) {
                        mip->completeness = 4;
                      }
                      vnp = CreateNewDescriptor (psep, Seq_descr_molinfo);
                      if (vnp != NULL) {
                        vnp->data.ptrvalue = (Pointer) mip;
                      }
                    }
                    /* the first protein may change the set/seq structure,
                    so goes through AddSeqEntryToSeqEntry */

                    if (gen_prod_set || last == NULL) {
                      descr = ExtractBioSourceAndPubs (sep);
                      AddSeqEntryToSeqEntry (sep, psep, TRUE);
                      ReplaceBioSourceAndPubs (sep, descr);
                      last = psep;
                    } else {
                      last->next = psep;
                      last = psep;
                    }
                    if (target == NULL) {
                      target = sep;
                      SaveSeqEntryObjMgrData (target, &omdptop, &omdata);
                      GetSeqEntryParent (target, &parentptr, &parenttype);
                    }
                    SetSeqFeatProduct (sfp, pbsp);
                    psep = SeqMgrGetSeqEntryForData (pbsp);
                    if (psep != NULL) {
                      last = psep;
                      prot = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
                      if (prot != NULL) {
                        prot->data.value.ptrvalue = (Pointer) prp;
                        SetSeqLocPartial (prot->location, partial5, partial3);
                        prot->partial = (Boolean) (partial5 || partial3);
                      }
                    }
                  }
                }
              }
            }
          }
        }
        if (ok_to_remove) {
          xref->data.value.ptrvalue = NULL;
          *(prev) = next;
          xref->next = NULL;
          xref->data.choice = 0;
          SeqFeatXrefFree (xref);
        } else {
          prev = &(xref->next);
        }
      } else {
        prev = &(xref->next);
      }
      xref = next;
    }
    sfp = sfp->next;
  }

  /* expand mat_peptide features with protein_id qualifiers into protein product sequences */

  last = NULL;
  sfp = first;
  while (sfp != NULL) {
    if (sfp->data.choice == SEQFEAT_PROT && sfp->product == NULL) {
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      gbq = sfp->qual;
      prevqual = (GBQualPtr PNTR) &(sfp->qual);
      id [0] = '\0';
      sip = NULL;
      while (gbq != NULL) {
        nextqual = gbq->next;
        if (StringICmp (gbq->qual, "protein_id") == 0) {
          if (StringDoesHaveText (id) && StringDoesHaveText (gbq->val)) {
            ErrPostEx (SEV_WARNING, ERR_FEATURE_QualWrongThisFeat,
                       "Protein protein_id %s replacing %s",
                       gbq->val, id);
          }
          *(prevqual) = gbq->next;
          gbq->next = NULL;
          StringNCpy_0 (id, gbq->val, sizeof (id));
          GBQualFree (gbq);
        } else {
          prevqual = (GBQualPtr PNTR) &(gbq->next);
        }
        gbq = nextqual;
      }
      if (! StringHasNoText (id)) {
        if (StringChr (id, '|') != NULL) {
          sip = SeqIdParse (id);
        } else if (force_local_id) {
          sprintf (lcl, "lcl|%s", id);
          sip = SeqIdParse (lcl);
        } else {
          adv = ValidateAccnDotVer (id);
          if (adv == 0 || adv == -5) {
            ptr = StringChr (id, '.');
            if (ptr != NULL) {
              *ptr = '\0';
              ptr++;
              if (sscanf (ptr, "%ld", &val) == 1) {
                version = (Uint4) val;
              }
            }
            sip = SeqIdFromAccession (id, version, NULL);
          } else {
            sprintf (lcl, "lcl|%s", id);
            sip = SeqIdParse (lcl);
          }
        }
      }
      if (sip != NULL) {
        protseq = GetSequenceByFeature (sfp);
        if (protseq == NULL && seq_fetch_failP != NULL) {
          *seq_fetch_failP = TRUE;
        }
        if (protseq != NULL) {
          i = (Int4) StringLen (protseq);
          bs = BSNew (i + 2);
          if (bs != NULL) {
            BSWrite (bs, (VoidPtr) protseq, (Int4) StringLen (protseq));
            pbsp = BioseqNew ();
            if (pbsp != NULL) {
              pbsp->repr = Seq_repr_raw;
              pbsp->mol = Seq_mol_aa;
              pbsp->seq_data_type = Seq_code_ncbieaa;
              pbsp->seq_data = (SeqDataPtr) bs;
              pbsp->length = BSLen (bs);
              bs = NULL;
              /*
              sep = GetTopSeqEntryForEntityID (entityID);
              */
              sep = GetBestTopParentForData (entityID, bsp);
              old = SeqEntrySetScope (sep);
              if (sip != NULL) {
                pbsp->id = sip;
              } else {
                pbsp->id = MakeNewProteinSeqIdEx (sfp->location, NULL, NULL, &ctr);
              }
              CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
              SeqMgrAddToBioseqIndex (pbsp);
              SeqEntrySetScope (old);
              psep = SeqEntryNew ();
              if (psep != NULL) {
                psep->choice = 1;
                psep->data.ptrvalue = (Pointer) pbsp;
                SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) pbsp, psep);
                mip = MolInfoNew ();
                if (mip != NULL) {
                  mip->biomol = MOLECULE_TYPE_PEPTIDE;
                  if (partial5 && partial3) {
                    mip->completeness = 5;
                  } else if (partial5) {
                    mip->completeness = 3;
                  } else if (partial3) {
                    mip->completeness = 4;
                  }
                  vnp = CreateNewDescriptor (psep, Seq_descr_molinfo);
                  if (vnp != NULL) {
                    vnp->data.ptrvalue = (Pointer) mip;
                  }
                }
                if (last == NULL) {
                  AddSeqEntryToSeqEntry (sep, psep, TRUE);
                  last = psep;
                } else {
                  last->next = psep;
                  last = psep;
                }
                SetSeqFeatProduct (sfp, pbsp);
                if (prp != NULL) {
                  prp2 = AsnIoMemCopy ((Pointer) prp,
                                       (AsnReadFunc) ProtRefAsnRead,
                                       (AsnWriteFunc) ProtRefAsnWrite);
                  if (prp2 != NULL) {
                    psep = SeqMgrGetSeqEntryForData (pbsp);
                    if (psep != NULL) {
                      prot = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
                      if (prot != NULL) {
                        prot->data.value.ptrvalue = prp2;
                        SetSeqLocPartial (prot->location, partial5, partial3);
                        prot->partial = (Boolean) (partial5 || partial3);
                      }
                    }
                  }
                }
              }
            }
          }
          protseq = MemFree (protseq);
        }
      }
    }
    sfp = sfp->next;
  }

  if (target != NULL) {
    SeqMgrLinkSeqEntry (target, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (target, omdptop, &omdata);
  }
}

NLM_EXTERN void PromoteXrefsEx (SeqFeatPtr sfp, BioseqPtr bsp, Uint2 entityID, Boolean include_stop,
                                Boolean remove_trailingX, Boolean gen_prod_set)

{
  PromoteXrefsExEx (sfp, bsp, entityID, include_stop, remove_trailingX, gen_prod_set, FALSE, NULL);
}

NLM_EXTERN void PromoteXrefs (SeqFeatPtr sfp, BioseqPtr bsp, Uint2 entityID)

{
  PromoteXrefsExEx (sfp, bsp, entityID, TRUE, FALSE, FALSE, FALSE, NULL);
}

/* begin BasicSeqEntryCleanup section */

static Boolean HasNoText (CharPtr str)

{
  Uchar  ch;    /* to use 8bit characters in multibyte languages */

  if (str != NULL) {
    ch = *str;
    while (ch != '\0') {
      if (ch > ' ') {
        return FALSE;
      }
      str++;
      ch = *str;
    }
  }
  return TRUE;
}

static Boolean AlreadyInVnpList (ValNodePtr head, ValNodePtr curr)

{
  if (head == NULL || curr == NULL) return FALSE;
  /* since we cannot sort these lists, must check against all previous entries */
  while (head != curr && head != NULL) {
    if (StringICmp (head->data.ptrvalue, curr->data.ptrvalue) == 0) return TRUE;
    head = head->next;
  }
  return FALSE;
}

NLM_EXTERN CharPtr TrimSpacesAndSemicolons (CharPtr str)

{
  CharPtr  amp;
  Uchar    ch;    /* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    if (ch != '\0' && (ch <= ' ' || ch == ';')) {
      while (ch != '\0' && (ch <= ' ' || ch == ';')) {
        ptr++;
        ch = *ptr;
      }
      while (ch != '\0') {
        *dst = ch;
        dst++;
        ptr++;
        ch = *ptr;
      }
      *dst = '\0';
    }
    amp = NULL;
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == '&') {
        amp = ptr;
        dst = NULL;
      } else if (ch <= ' ') {
        if (dst == NULL) {
          dst = ptr;
        }
        amp = NULL;
      } else if (ch == ';') {
        if (dst == NULL && amp == NULL) {
          dst = ptr;
        }
      } else {
        dst = NULL;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

NLM_EXTERN CharPtr TrimSpacesAndJunkFromEnds (
  CharPtr str,
  Boolean allowEllipsis
)

{
  Uchar    ch;    /* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  Boolean  isPeriod;
  Boolean  isTilde;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    if (ch != '\0' && ch <= ' ') {
      while (ch != '\0' && ch <= ' ') {
        ptr++;
        ch = *ptr;
      }
      while (ch != '\0') {
        *dst = ch;
        dst++;
        ptr++;
        ch = *ptr;
      }
      *dst = '\0';
    }
    dst = NULL;
    ptr = str;
    ch = *ptr;
    isPeriod = FALSE;
    isTilde = FALSE;
    while (ch != '\0') {
      if (ch <= ' ' || ch == '.' || ch == ',' || ch == '~' || ch == ';') {
        if (dst == NULL) {
          dst = ptr;
        }
        isPeriod = (Boolean) (isPeriod || ch == '.');
        isTilde = (Boolean) (isTilde || ch == '~');
      } else {
        dst = NULL;
        isPeriod = FALSE;
        isTilde = FALSE;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      /* allow one period at end */
      if (isPeriod) {
        *dst = '.';
        dst++;
        /* ellipsis are now okay */
        if (allowEllipsis && *dst == '.' && dst [1] == '.') {
          dst += 2;
        }
      } else if (isTilde) {
        /* allow double tilde at end */
        if (*dst == '~' && dst [1] == '~') {
          dst += 2;
        }
      }
      *dst = '\0';
    }
  }
  return str;
}

static CharPtr TrimSpacesSemicolonsAndCommas (CharPtr str)

{
  CharPtr  amp;
  Uchar    ch;    /* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    if (ch != '\0' && (ch <= ' ' || ch == ';' || ch == ',')) {
      while (ch != '\0' && (ch <= ' ' || ch == ';' || ch == ',')) {
        ptr++;
        ch = *ptr;
      }
      while (ch != '\0') {
        *dst = ch;
        dst++;
        ptr++;
        ch = *ptr;
      }
      *dst = '\0';
    }
    amp = NULL;
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == '&') {
        amp = ptr;
        dst = NULL;
      } else if (ch <= ' ') {
        if (dst == NULL) {
          dst = ptr;
        }
        amp = NULL;
      } else if (ch == ';') {
        if (dst == NULL && amp == NULL) {
          dst = ptr;
        }
      } else if (ch == ',') {
        if (dst == NULL) {
          dst = ptr;
        }
        amp = NULL;
      } else {
        dst = NULL;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

static CharPtr TrimSpacesSemicolonsAndCommasAndParens (CharPtr str)

{
  CharPtr  amp;
  Uchar    ch;    /* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    if (ch != '\0' && (ch <= ' ' || ch == ';' || ch == ',' || ch == '(')) {
      while (ch != '\0' && (ch <= ' ' || ch == ';' || ch == ',' || ch == '(')) {
        ptr++;
        ch = *ptr;
      }
      while (ch != '\0') {
        *dst = ch;
        dst++;
        ptr++;
        ch = *ptr;
      }
      *dst = '\0';
    }
    amp = NULL;
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == '&') {
        amp = ptr;
        dst = NULL;
      } else if (ch <= ' ') {
        if (dst == NULL) {
          dst = ptr;
        }
        amp = NULL;
      } else if (ch == ';') {
        if (dst == NULL && amp == NULL) {
          dst = ptr;
        }
      } else if (ch == ',') {
        if (dst == NULL) {
          dst = ptr;
        }
        amp = NULL;
      } else if (ch == ')') {
        if (dst == NULL) {
          dst = ptr;
        }
        amp = NULL;
      } else {
        dst = NULL;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

static CharPtr TrimFlankingQuotes (CharPtr str)

{
  size_t  len;

  if (str != NULL && str [0] != '\0') {
    len = StringLen (str);
    while (len > 0) {
      if (str [0] == '"' && str [len - 1] == '"') {
        str [0] = ' ';
        str [len - 1] = ' ';
      } else if (str [0] == '\'' && str [len - 1] == '\'') {
        str [0] = ' ';
        str [len - 1] = ' ';
      } else {
        return str;
      }
      TrimSpacesAroundString (str);
      len = StringLen (str);
    }
  }
  return str;
}

static void RemoveFlankingQuotes (CharPtr PNTR strp)

{
  if (strp == NULL) return;
  if (*strp == NULL) return;
  TrimFlankingQuotes (*strp);
  if (HasNoText (*strp)) {
    *strp = MemFree (*strp);
  }
}

static void RemoveFlankingQuotesList (ValNodePtr PNTR vnpp)

{
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       vnp;

  if (vnpp == NULL) return;
  prev = vnpp;
  vnp = *vnpp;
  while (vnp != NULL) {
    next = vnp->next;
    TrimFlankingQuotes (vnp->data.ptrvalue);
    if (HasNoText (vnp->data.ptrvalue) || AlreadyInVnpList (*vnpp, vnp)) {
      *prev = vnp->next;
      vnp->next = NULL;
      ValNodeFreeData (vnp);
    } else {
      prev = &(vnp->next);
    }
    vnp = next;
  }
}

static void CleanVisString (CharPtr PNTR strp)

{
  if (strp == NULL) return;
  if (*strp == NULL) return;
  TrimSpacesSemicolonsAndCommas (*strp);
  if (HasNoText (*strp)) {
    *strp = MemFree (*strp);
  }
}

static void CleanVisStringAndCompress (CharPtr PNTR strp)

{
  if (strp == NULL) return;
  if (*strp == NULL) return;
  TrimSpacesSemicolonsAndCommas (*strp);
  Asn2gnbkCompressSpaces (*strp);
  if (HasNoText (*strp)) {
    *strp = MemFree (*strp);
  }
}

static void CleanVisStringJunk (CharPtr PNTR strp)

{
  if (strp == NULL) return;
  if (*strp == NULL) return;
  TrimSpacesAndJunkFromEnds (*strp, TRUE);
  if (HasNoText (*strp)) {
    *strp = MemFree (*strp);
  }
}

static void CleanVisStringJunkAndCompress (CharPtr PNTR strp)

{
  if (strp == NULL) return;
  if (*strp == NULL) return;
  TrimSpacesAndJunkFromEnds (*strp, TRUE);
  Asn2gnbkCompressSpaces (*strp);
  if (HasNoText (*strp)) {
    *strp = MemFree (*strp);
  }
}

static void CleanDoubleQuote (CharPtr str)

{
  Char  ch;

  if (str == NULL) return;
  ch = *str;
  while (ch != '\0') {
    if (ch == '"') {
      *str = '\'';
    }
    str++;
    ch = *str;
  }
}

static CharPtr RemoveSpacesBetweenTildes (CharPtr str)

{
  Char     ch;
  CharPtr  dst;
 CharPtr  ptr;
  CharPtr  tmp;

  if (str == NULL || str [0] == '\0') return str;

  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    *dst = ch;
    dst++;
    ptr++;
    if (ch == '~') {
      tmp = ptr;
      ch = *tmp;
      while (ch != 0 && ch <= ' ') {
        tmp++;
        ch = *tmp;
      }
      if (ch == '~') {
        ptr = tmp;
      }
    }
    ch = *ptr;
  }
  *dst = '\0';

  return str;
}

static void CleanVisStringList (ValNodePtr PNTR vnpp)

{
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       vnp;

  if (vnpp == NULL) return;
  prev = vnpp;
  vnp = *vnpp;
  while (vnp != NULL) {
    next = vnp->next;
    TrimSpacesSemicolonsAndCommas (vnp->data.ptrvalue);
    if (HasNoText (vnp->data.ptrvalue) || AlreadyInVnpList (*vnpp, vnp)) {
      *prev = vnp->next;
      vnp->next = NULL;
      ValNodeFreeData (vnp);
    } else {
      prev = &(vnp->next);
    }
    vnp = next;
  }
}

static void CleanVisStringListAndCompress (ValNodePtr PNTR vnpp)

{
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       vnp;

  if (vnpp == NULL) return;
  prev = vnpp;
  vnp = *vnpp;
  while (vnp != NULL) {
    next = vnp->next;
    TrimSpacesSemicolonsAndCommas (vnp->data.ptrvalue);
    Asn2gnbkCompressSpaces (vnp->data.ptrvalue);
    if (HasNoText (vnp->data.ptrvalue) || AlreadyInVnpList (*vnpp, vnp)) {
      *prev = vnp->next;
      vnp->next = NULL;
      ValNodeFreeData (vnp);
    } else {
      prev = &(vnp->next);
    }
    vnp = next;
  }
}

static Boolean AlreadyInVnpListCaseSensitive (ValNodePtr head, ValNodePtr curr)

{
  if (head == NULL || curr == NULL) return FALSE;
  /* since we cannot sort these lists, must check against all previous entries */
  while (head != curr && head != NULL) {
    if (StringCmp (head->data.ptrvalue, curr->data.ptrvalue) == 0) return TRUE;
    head = head->next;
  }
  return FALSE;
}

static void CleanVisStringListCaseSensitive (ValNodePtr PNTR vnpp)

{
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       vnp;

  if (vnpp == NULL) return;
  prev = vnpp;
  vnp = *vnpp;
  while (vnp != NULL) {
    next = vnp->next;
    TrimSpacesSemicolonsAndCommas (vnp->data.ptrvalue);
    if (HasNoText (vnp->data.ptrvalue) || AlreadyInVnpListCaseSensitive (*vnpp, vnp)) {
      *prev = vnp->next;
      vnp->next = NULL;
      ValNodeFreeData (vnp);
    } else {
      prev = &(vnp->next);
    }
    vnp = next;
  }
}

static void CleanDoubleQuoteList (ValNodePtr vnp)

{
  while (vnp != NULL) {
    CleanDoubleQuote ((CharPtr) vnp->data.ptrvalue);
    vnp = vnp->next;
  }
}

static Boolean HandledGBQualOnGene (SeqFeatPtr sfp, GBQualPtr gbq)

{
  Int2        choice = 0;
  GeneRefPtr  grp;

  if (StringICmp (gbq->qual, "map") == 0) {
    choice = 2;
  } else if (StringICmp (gbq->qual, "allele") == 0) {
    choice = 3;
  } else if (StringICmp (gbq->qual, "locus_tag") == 0) {
    choice = 4;
  } else if (StringICmp (gbq->qual, "old_locus_tag") == 0) {
    choice = 5;
  }
  if (choice > 0) {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (grp == NULL) return FALSE;
    switch (choice) {
      case 2 :
        if (grp->maploc != NULL) return FALSE;
        if (StringHasNoText (gbq->val)) return FALSE;
        grp->maploc = StringSave (gbq->val);
        break;
      case 3 :
        if (StringHasNoText (gbq->val)) return FALSE;
        if (grp->allele != NULL) {
          if (StringICmp (gbq->val, grp->allele) == 0) return TRUE;
          return FALSE;
        }
        grp->allele = StringSave (gbq->val);
        break;
      case 4 :
        if (grp->locus_tag != NULL) return FALSE;
        if (StringHasNoText (gbq->val)) return FALSE;
        grp->locus_tag = StringSave (gbq->val);
        break;
      case 5 :
/* removed by indexer request */
/*        if (StringHasNoText (gbq->val)) return FALSE;
 *       if (grp->locus_tag != NULL) {
 *         if (StringICmp (gbq->val, grp->locus_tag) == 0) return TRUE;
 *         return FALSE;
 *       }
 */
        return FALSE;
        break;
      default :
        break;
    }
    return TRUE;
  }
  return FALSE;
}

/* code break parser functions from the flatfile parser */

static Uint1 GetQualValueAa (CharPtr qval)

{
   CharPtr  str, eptr, ptr;
   Uint1    aa;

    str = StringStr(qval, "aa:");
    if (str != NULL) {
        str += 3;
    } else {
        ErrPostEx (SEV_WARNING, ERR_QUALIFIER_InvalidDataFormat,
                   "bad transl_except %s", qval);
        str = StringStr(qval, ",");
        if (str != NULL) {
            str = StringStr(str, ":");
            if (str != NULL) {
              str++;
            }
        }
    }

    if (str == NULL) return (Uint1) 'X';

       while (*str == ' ')
           ++str;
       for (eptr = str; *eptr != ')' && *eptr != ' ' && *eptr != '\0';  eptr++) continue;

    ptr = TextSave(str, eptr-str);
    aa = ValidAminoAcid(ptr);
    MemFree(ptr);  

    return (aa);
}

static CharPtr SimpleValuePos (CharPtr qval)

{
   CharPtr bptr, eptr;

   if ((bptr = StringStr(qval, "(pos:")) == NULL) {
           return NULL;
   }
    
   bptr += 5;
   while (*bptr == ' ')
       ++bptr;
   for (eptr = bptr; *eptr != ',' && *eptr != '\0'; eptr++) continue;

   return (TextSave(bptr, eptr-bptr));
}

extern Boolean ParseAnticodon (SeqFeatPtr sfp, CharPtr val, Int4 offset);
extern Boolean ParseAnticodon (SeqFeatPtr sfp, CharPtr val, Int4 offset)

{
  Int4       diff;
  Int2       j;
  Boolean    locmap;
  int        num_errs;
  CharPtr    pos;
  Boolean    pos_range = FALSE;
  RnaRefPtr  rrp;
  SeqIntPtr  sintp;
  SeqIdPtr   sip;
  Boolean    sitesmap;
  SeqLocPtr  slp;
  SeqPntPtr  spp;
  Uint1      strand;
  Int4       temp;
  tRNAPtr    trp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return FALSE;
  if (StringHasNoText (val)) return FALSE;

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL) return FALSE;

  if (rrp->ext.choice == 0 && rrp->ext.value.ptrvalue == NULL) {
    rrp->ext.choice = 2;
    trp = (tRNAPtr) MemNew (sizeof (tRNA));
    rrp->ext.value.ptrvalue = (Pointer) trp;
    if (trp != NULL) {
      trp->aatype = 2;
      for (j = 0; j < 6; j++) {
        trp->codon [j] = 255;
      }
    }
  }
  if (rrp->ext.choice != 2) return FALSE;

  trp = (tRNAPtr) rrp->ext.value.ptrvalue;
  if (trp == NULL) return FALSE;
      
  /* find SeqId to use */
  sip = SeqLocId (sfp->location);
  if (sip == NULL) {
    slp = SeqLocFindNext (sfp->location, NULL);
    if (slp != NULL) {
      sip = SeqLocId (slp);
    }
  }
  if (sip == NULL) return FALSE;

  /* parse location */
  pos = SimpleValuePos (val);
  if (pos == NULL) {
    ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
               "anticodon parsing failed, %s, drop the anticodon", val);
    return FALSE;
  }

  trp->anticodon = Nlm_gbparseint (pos, &locmap, &sitesmap, &num_errs, sip);
  if (trp->anticodon == NULL) {
    ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
               "anticodon parsing failed, %s, drop the anticodon", pos);
    MemFree (pos);
    return FALSE;
  }

  if (trp->anticodon->choice == SEQLOC_PNT) {
    /* allow a single point */
    spp = trp->anticodon->data.ptrvalue;
    if (spp != NULL) {
      spp->point += offset;
    }
  }
  if (trp->anticodon->choice == SEQLOC_INT) {
    sintp = trp->anticodon->data.ptrvalue;
    if (sintp == NULL) {
      MemFree (pos);
      return FALSE;
    }
    sintp->from += offset;
    sintp->to += offset;
    if (sintp->from > sintp->to) {
      temp = sintp->from;
      sintp->from = sintp->to;
      sintp->to = temp;
    }
    sintp->strand = SeqLocStrand (sfp->location);
    strand = sintp->strand;
    diff = SeqLocStop(trp->anticodon) - SeqLocStart(trp->anticodon); /* SeqLocStop/Start does not do what you think */
    /*
    if ((diff != 2 && (strand != Seq_strand_minus)) ||
        (diff != -2 && (strand == Seq_strand_minus))) {
      pos_range = TRUE;
    }
    */
    if (diff != 2) {
      pos_range = TRUE;
    }
    if (num_errs > 0 || pos_range) {
      ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
                 "anticodon range is wrong, %s, drop the anticodon", pos);
      MemFree (pos);
      return FALSE;
    }
    if (SeqLocCompare (sfp->location, trp->anticodon) != SLC_B_IN_A) {
      ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
                 "/anticodon not in tRNA: %s", val);
      MemFree (pos);
      return FALSE;
    }
  }

  MemFree (pos);

  return TRUE;
}

extern Boolean ParseCodeBreak (SeqFeatPtr sfp, CharPtr val, Int4 offset);
extern Boolean ParseCodeBreak (SeqFeatPtr sfp, CharPtr val, Int4 offset)

{
  CodeBreakPtr  cbp;
  CdRegionPtr   crp;
  Int4          diff;
  CodeBreakPtr  lastcbp;
  Boolean       locmap;
  int           num_errs;
  CharPtr       pos;
  Boolean       pos_range = FALSE;
  SeqIntPtr     sintp;
  SeqIdPtr      sip;
  Boolean       sitesmap;
  SeqLocPtr     slp;
  SeqPntPtr     spp;
  Uint1         strand;
  Int4          temp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return FALSE;
  if (StringHasNoText (val)) return FALSE;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL) return FALSE;

  /* find SeqId to use */
  sip = SeqLocId (sfp->location);
  if (sip == NULL) {
    slp = SeqLocFindNext (sfp->location, NULL);
    if (slp != NULL) {
      sip = SeqLocId (slp);
    }
  }
  if (sip == NULL) return FALSE;

  cbp = CodeBreakNew ();
  if (cbp == NULL) return FALSE;
  cbp->aa.choice = 1; /* ncbieaa */
  cbp->aa.value.intvalue = (Int4) GetQualValueAa (val);

  /* parse location */
  pos = SimpleValuePos (val);
  if (pos == NULL) {
    ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
               "transl_except parsing failed, %s, drop the transl_except", val);
    return FALSE;
  }
  cbp->loc = Nlm_gbparseint (pos, &locmap, &sitesmap, &num_errs, sip);
  if (cbp->loc == NULL) {
    CodeBreakFree (cbp);
    ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
               "transl_except parsing failed, %s, drop the transl_except", pos);
    MemFree (pos);
    return FALSE;
  }
  if (cbp->loc->choice == SEQLOC_PNT) {
    /* allow a single point */
    spp = cbp->loc->data.ptrvalue;
    if (spp != NULL) {
      spp->point += offset;
    }
  }
  if (cbp->loc->choice == SEQLOC_INT) {
    sintp = cbp->loc->data.ptrvalue;
    if (sintp == NULL) {
      MemFree (pos);
      return FALSE;
    }
    sintp->from += offset;
    sintp->to += offset;
    if (sintp->from > sintp->to) {
      temp = sintp->from;
      sintp->from = sintp->to;
      sintp->to = temp;
    }
    sintp->strand = SeqLocStrand (sfp->location);
    strand = sintp->strand;
    diff = SeqLocStop(cbp->loc) - SeqLocStart(cbp->loc); /* SeqLocStop/Start does not do what you think */
    /*
    if ((diff != 2 && (strand != Seq_strand_minus)) ||
        (diff != -2 && (strand == Seq_strand_minus))) {
      pos_range = TRUE;
    }
    */
    if (diff != 2) {
      pos_range = TRUE;
    }
    if (num_errs > 0 || pos_range) {
      CodeBreakFree (cbp);
      ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
                 "transl_except range is wrong, %s, drop the transl_except", pos);
      MemFree (pos);
      return FALSE;
    }
    if (SeqLocCompare (sfp->location, cbp->loc) != SLC_B_IN_A) {
      CodeBreakFree (cbp);
      ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
                 "/transl_except not in CDS: %s", val);
      MemFree (pos);
      return FALSE;
    }
  }

  /* add to code break list */
  lastcbp = crp->code_break;
  if (lastcbp == NULL) {
    crp->code_break = cbp;
  } else {
     while (lastcbp->next != NULL) {
      lastcbp = lastcbp->next;
    }
    lastcbp->next = cbp;
  }
  MemFree (pos);
  return TRUE;
}

static Boolean CodonsAlreadyInOrder (tRNAPtr trp)

{
  Int2  i, j;

  if (trp == NULL) return TRUE;
  for (i = 0, j = 1; i < 5; i++, j++) {
    if (trp->codon [i] > trp->codon [j]) return FALSE;
  }
  return TRUE;
}

static int LIBCALLBACK SortCodons (VoidPtr ptr1, VoidPtr ptr2)

{
  Uint1  codon1, codon2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  codon1 = *((Uint1Ptr) ptr1);
  codon2 = *((Uint1Ptr) ptr2);
  if (codon1 > codon2) {
    return 1;
  } else if (codon1 < codon2) {
    return -1;
  }
  return 0;
}

static void UniqueCodons (tRNAPtr trp)

{
  Int2   i, j;
  Uint1  last = 255, next;

  if (trp == NULL) return;

  for (i = 0, j = 0; i < 6; i++) {
    next = trp->codon [i];
    if (next != last) {
      trp->codon [j] = next;
      last = next;
      j++;
    }
  }
  while (j < 6) {
    trp->codon [j] = 255;
    j++;
  }
}

static CharPtr  codonLetterExpand [] =
{
  "?", "A", "C", "AC",
  "G", "AG", "CG", "ACG",
  "T", "AT", "CT", "ACT",
  "GT", "AGT", "CGT", "ACGT",
  NULL
};

NLM_EXTERN Boolean ParseDegenerateCodon (tRNAPtr trp, Uint1Ptr codon)

{
  Uint1    ch;
  Uint1    chrToInt [256];
  Int2     k;
  Uint1    i, j;
  Uint1    idx;
  CharPtr  intToChr = "?ACMGRSVTWYHKDBN";
  CharPtr  ptr, str;

  if (trp == NULL || codon == NULL) return FALSE;

  for (i = 0; i < 2; i++) {
    ch = codon [i];
    if (ch != 'A' && ch != 'C' && ch != 'G' && ch != 'T') return FALSE;
  }

  for (k = 0; k < 256; k++) {
    chrToInt [k] = 0;
  }
  for (i = 1; i < 16; i++) {
    ch = intToChr [i];
    chrToInt [(int) ch] = i;
  }

  idx = chrToInt [(int) codon [2]];
  if (idx > 15) return FALSE;

  str = codonLetterExpand [idx];
  ptr = str;
  ch = *ptr;
  j = 0;
  codon [3] = '\0';
  while (ch != '\0' && j < 6) {
    codon [2] = ch;
    trp->codon [j] = IndexForCodon (codon, Seq_code_iupacna);
    ptr++;
    ch = *ptr;
    j++;
  }

  return TRUE;
}

static void CleanupTrna (SeqFeatPtr sfp, tRNAPtr trp)

{
  Uint1           aa = 0;
  Char            codon [16];
  Uint1           curraa;
  Uint1           from = 0;
  Int2            i;
  Int2            j;
  Boolean         justTrnaText;
  Boolean         okayToFree = TRUE;
  SeqMapTablePtr  smtp;
  CharPtr         str;
  Uint1           trpcodon [6];

  /* look for tRNA-OTHER with actual amino acid in comment */

  if (trp == NULL) return;

  if (sfp != NULL && sfp->comment != NULL && trp->codon [0] == 255) {
    codon [0] = '\0';
    if (StringNICmp (sfp->comment, "codon recognized: ", 18) == 0) {
      StringNCpy_0 (codon, sfp->comment + 18, sizeof (codon));
    } else if (StringNICmp (sfp->comment, "codons recognized: ", 19) == 0) {
      StringNCpy_0 (codon, sfp->comment + 19, sizeof (codon));
    }
    if (StringDoesHaveText (codon)) {
      if (StringLen (codon) > 3 && codon [3] == ';') {
        codon [3] = '\0';
        okayToFree = FALSE;
      }
      if (StringLen (codon) == 3) {
        for (i = 0; i < 3; i++) {
          if (codon [i] == 'U') {
            codon [i] = 'T';
          }
        }
        if (ParseDegenerateCodon (trp, (Uint1Ptr) codon)) {
          if (okayToFree) {
            sfp->comment = MemFree (sfp->comment);
          } else {
            str = StringSave (sfp->comment + 22);
            TrimSpacesAroundString (str);
            sfp->comment = MemFree (sfp->comment);
            if (StringHasNoText (str)) {
              str = MemFree (str);
            }
            sfp->comment = str;
          }
        }
      }
    }
  }

  if (! CodonsAlreadyInOrder (trp)) {
    StableMergeSort ((VoidPtr) &(trp->codon), 6, sizeof (Uint1), SortCodons);
  }
  UniqueCodons (trp);

  /* now always switch iupacaa to ncbieaa (was just for selenocysteine) */

  if (trp->aatype == 1 /* && trp->aa == 'U' */) {
    trp->aatype = 2;
  }

  if (sfp == NULL || sfp->comment == NULL) return;

  if (trp->aatype == 2) {
    aa = trp->aa;
  } else {
    switch (trp->aatype) {
      case 0 :
        from = 0;
        break;
      case 1 :
        from = Seq_code_iupacaa;
        break;
      case 2 :
        from = Seq_code_ncbieaa;
        break;
      case 3 :
        from = Seq_code_ncbi8aa;
        break;
      case 4 :
        from = Seq_code_ncbistdaa;
        break;
      default:
        break;
    }
    smtp = SeqMapTableFind (Seq_code_ncbieaa, from);
    if (smtp != NULL) {
      aa = SeqMapTableConvert (smtp, trp->aa);
    }
  }
  if (aa != 'X') {
    curraa = ParseTRnaString (sfp->comment, &justTrnaText, trpcodon, TRUE);
    if (aa == 0 && curraa != 0) {
      aa = curraa;
      trp->aa = curraa;
      trp->aatype = 2;
    }
    if (aa != 0 && aa == curraa) {
      if (justTrnaText) {
        for (j = 0; j < 6; j++) {
          if (trp->codon [j] == 255) {
            trp->codon [j] = trpcodon [j];
          }
        }
        if (StringCmp (sfp->comment, "fMet") != 0) {
          sfp->comment = MemFree (sfp->comment);
        }
      }
    }
    return;
  }
  aa = ParseTRnaString (sfp->comment, &justTrnaText, trpcodon, TRUE);
  if (aa == 0) return;
  trp->aa = aa;
  trp->aatype = 2;
  if (justTrnaText) {
    for (j = 0; j < 6; j++) {
      if (trp->codon [j] == 255) {
        trp->codon [j] = trpcodon [j];
      }
    }
    if (StringCmp (sfp->comment, "fMet") != 0) {
      sfp->comment = MemFree (sfp->comment);
    }
  }
}

NLM_EXTERN SeqFeatPtr LIBCALL GetBestProteinFeatureUnindexed (SeqLocPtr product)

{
  BioseqPtr    bsp;
  SeqFeatPtr   prot = NULL;
  SeqAnnotPtr  sap;
  SeqFeatPtr   tmp;
  ValNode      vn;

  if (product == NULL) return NULL;
  bsp = BioseqFindFromSeqLoc (product);
  if (bsp == NULL || bsp->repr != Seq_repr_raw) return NULL;
  vn.choice = SEQLOC_WHOLE;
  vn.data.ptrvalue = (Pointer) SeqIdFindBest (bsp->id, 0);
  vn.next = NULL;
  for (sap = bsp->annot; sap != NULL && prot == NULL; sap = sap->next) {
    if (sap->type == 1) {
      for (tmp = (SeqFeatPtr) sap->data; tmp != NULL && prot == NULL; tmp = tmp->next) {
        if (tmp->data.choice == SEQFEAT_PROT) {
          if (SeqLocCompare (tmp->location, &vn)) {
            /* find first protein feature packaged on and located on bioseq */
            prot = tmp;
          }
        }
      }
    }
  }
  return prot;
}

static void CleanupECNumber (CharPtr str)

{
  size_t len;

  len = StringLen (str);
  if (len < 1) return;
  if (str [len - 1] == '.') {
    str [len - 1] = ' ';
  }
  if (StringNICmp (str, "EC ", 3) == 0) {
    str [0] = ' ';
    str [1] = ' ';
  } else if (StringNICmp (str, "EC:", 3) == 0) {
    str [0] = ' ';
    str [1] = ' ';
    str [2] = ' ';
  }
  TrimSpacesAroundString (str);
}

static Boolean HandledGBQualOnCDS (SeqFeatPtr sfp, GBQualPtr gbq, ValNodePtr PNTR afterMe)

{
  Int2            choice = 0;
  CdRegionPtr     crp;
  Uint1           frame;
  ValNodePtr      gcp;
  ValNodePtr      prev;
  SeqFeatPtr      prot;
  ProtRefPtr      prp = NULL;
  Char            str [16];
  Int4            transl_table;
  int             val;
  ValNodePtr      vnp;
  SeqFeatXrefPtr  xref;

  if (StringICmp (gbq->qual, "product") == 0) {
    choice = 1;
  } else if (StringICmp (gbq->qual, "function") == 0) {
    choice = 2;
  } else if (StringICmp (gbq->qual, "EC_number") == 0) {
    choice = 3;
  } else if (StringICmp (gbq->qual, "prot_note") == 0) {
    choice = 4;
  }
  if (choice > 0) {
    prot = GetBestProteinFeatureUnindexed (sfp->product);
    if (prot != NULL) {
      prp = (ProtRefPtr) prot->data.value.ptrvalue;
    }
    if (prp == NULL) {
      /* otherwise make cross reference */
      xref = sfp->xref;
      while (xref != NULL && xref->data.choice != SEQFEAT_PROT) {
        xref = xref->next;
      }
      if (xref == NULL) {
        prp = ProtRefNew ();
        if (prp == NULL) return FALSE;
        xref = SeqFeatXrefNew ();
        if (xref == NULL) return FALSE;
        xref->data.choice = SEQFEAT_PROT;
        xref->data.value.ptrvalue = (Pointer) prp;
        xref->next = sfp->xref;
        sfp->xref = xref;
      }
      if (xref != NULL) {
        prp = (ProtRefPtr) xref->data.value.ptrvalue;
      }
    }
    if (prp == NULL) return FALSE;
    switch (choice) {
      case 1 :
        if (prot != NULL && prot->data.value.ptrvalue != NULL) {
          if (*afterMe == NULL) {
            /* if protein product exists, product gbqual becomes first name */
            vnp = ValNodeCopyStr (NULL, 0, gbq->val);
            if (vnp != NULL) {
              vnp->next = prp->name;
              prp->name = vnp;
            }
            *afterMe = vnp;
          } else {
            vnp = ValNodeCopyStr (NULL, 0, gbq->val);
            prev = *afterMe;
            if (vnp != NULL) {
              vnp->next = prev->next;
              prev->next = vnp;
            }
            *afterMe = vnp;
          }
        } else {
          /* if local xref, append to name */
          ValNodeCopyStr (&(prp->name), 0, gbq->val);
        }
        break;
      case 2 :
        ValNodeCopyStr (&(prp->activity), 0, gbq->val);
        break;
      case 3 :
        ValNodeCopyStr (&(prp->ec), 0, gbq->val);
        break;
      case 4 :
        if (prot == NULL) {
          return FALSE;
        } else {
          prot->comment = StringSave (gbq->val);
        }
        break;
      default :
        break;
    }
    return TRUE;
  }

  if (StringICmp (gbq->qual, "transl_except") == 0) {
    return ParseCodeBreak (sfp, gbq->val, 0);
  }

  if (StringICmp (gbq->qual, "codon_start") == 0) {
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if (crp != NULL) {
      frame = crp->frame;
      if (frame == 0) {
        StringNCpy_0 (str, gbq->val, sizeof (str));
        if (sscanf (str, "%d", &val) == 1) {
          if (val > 0 && val < 4) {
            crp->frame = (Uint1) val;
            return TRUE;
          }
        }
        frame = 1;
      }
      sprintf (str, "%d", (int) frame);
      if (StringICmp (str, gbq->val) == 0) {
        return TRUE;
      } else if (sfp->pseudo && sfp->product == NULL) {
        StringNCpy_0 (str, gbq->val, sizeof (str));
        if (sscanf (str, "%d", &val) == 1) {
          if (val > 0 && val < 4) {
            crp->frame = (Uint1) val;
            return TRUE;
          }
        }
      }
    }
  }

  if (StringICmp (gbq->qual, "transl_table") == 0) {
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if (crp != NULL) {
      transl_table = 0;
      gcp = crp->genetic_code;
      if (gcp != NULL) {
        for (vnp = gcp->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
          if (vnp->choice == 2 && vnp->data.intvalue != 0) {
            transl_table = vnp->data.intvalue;
          }
        }
        if (transl_table == 0) {
          transl_table = 1;
        }
        sprintf (str, "%ld", (long) transl_table);
        if (StringICmp (str, gbq->val) == 0) {
          return TRUE;
        }
      } else {
        StringNCpy_0 (str, gbq->val, sizeof (str));
        if (sscanf (str, "%d", &val) == 1) {
          vnp = ValNodeNew (NULL);
          if (vnp != NULL) {
            vnp->choice = 2;
            vnp->data.intvalue = (Int4) val;
            gcp = GeneticCodeNew ();
            if (gcp != NULL) {
              gcp->data.ptrvalue = vnp;
              crp->genetic_code = gcp;
              return TRUE;
            }
          }
        }
      }
    }
  }

  if (StringICmp (gbq->qual, "translation") == 0) {
    return TRUE;
  }

  return FALSE;
}


static Boolean HandledGBQualOnRNA (SeqFeatPtr sfp, GBQualPtr gbq, Boolean isEmblOrDdbj)

{
  Uint1      aa;
  BioseqPtr  bsp;
  Uint1      codon [6];
  Boolean    emptyRNA;
  Int4       from;
  Boolean    is_fMet = FALSE;
  Boolean    is_std_name = FALSE;
  Int2       j;
  Boolean    justTrnaText;
  size_t     len;
  CharPtr    name;
  CharPtr    ptr;
  RnaRefPtr  rrp;
  SeqIntPtr  sintp;
  SeqIdPtr   sip;
  CharPtr    str;
  Char       tmp [64];
  Int4       to;
  tRNAPtr    trp;
  long int   val;

  is_std_name = (Boolean) (StringICmp (gbq->qual, "standard_name") == 0);
  if (StringICmp (gbq->qual, "product") == 0 ||
      (is_std_name && (! isEmblOrDdbj) )) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rrp == NULL) return FALSE;
    if (rrp->type == 0) {
      rrp->type = 255;
    }
    if (rrp->type == 255 && is_std_name) return FALSE;
    if (rrp->ext.choice == 1) {
      name = (CharPtr) rrp->ext.value.ptrvalue;
      if (StringHasNoText (name)) {
        rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
        rrp->ext.choice = 0;
      }
    }
    if (rrp->ext.choice == 2) {
      trp = (tRNAPtr) rrp->ext.value.ptrvalue;
      if (trp != NULL) {
        if (trp->aatype == 0 && trp->aa == 0 && trp->anticodon == NULL) {
          emptyRNA = TRUE;
          for (j = 0; j < 6; j++) {
            if (trp->codon [j] != 255) {
              emptyRNA = FALSE;
            }
          }
          if (emptyRNA) {
            rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
            rrp->ext.choice = 0;
          }
        }
      }
    }
    if (rrp->type == 3 && rrp->ext.choice == 1) {
      name = (CharPtr) rrp->ext.value.ptrvalue;
      aa = ParseTRnaString (name, &justTrnaText, codon, FALSE);
      if (aa != 0) {
        is_fMet = (Boolean) (StringStr (name, "fMet") != NULL);
        rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
        trp = (tRNAPtr) MemNew (sizeof (tRNA));
        if (trp != NULL) {
          trp->aatype = 2;
          for (j = 0; j < 6; j++) {
            trp->codon [j] = 255;
          }
          if (justTrnaText) {
            for (j = 0; j < 6; j++) {
              trp->codon [j] = codon [j];
            }
          }
          trp->aa = aa;
          rrp->ext.choice = 2;
          rrp->ext.value.ptrvalue = (Pointer) trp;
          if (aa == 'M') {
            if (is_fMet) {
              if (sfp->comment == NULL) {
                sfp->comment = StringSave ("fMet");
              } else {
                len = StringLen (sfp->comment) + StringLen ("fMet") + 5;
                str = MemNew (sizeof (Char) * len);
                StringCpy (str, sfp->comment);
                StringCat (str, "; ");
                StringCat (str, "fMet");
                sfp->comment = MemFree (sfp->comment);
                sfp->comment = str;
              }
            }
          }
          CleanupTrna (sfp, trp);
        }
      }
    }
    if (rrp->type == 3 && rrp->ext.choice == 0) {
      AddQualifierToFeature (sfp, "product", gbq->val);
      return TRUE;
    }
    if (rrp->type == 3 && rrp->ext.choice == 2) {
      trp = (tRNAPtr) rrp->ext.value.ptrvalue;
      if (trp != NULL && trp->aatype == 2) {
        if (trp->aa == ParseTRnaString (gbq->val, NULL, NULL, FALSE)) {
          return TRUE;
        }
      }
    }
    if (rrp->ext.choice != 0 && rrp->ext.choice != 1) return FALSE;
    name = (CharPtr) rrp->ext.value.ptrvalue;
    if (! HasNoText (name)) {
      if (StringICmp (name, gbq->val) == 0) {
        return TRUE;
      }
      str = StringStr (gbq->val, "rDNA");
      if (str != NULL) {
        str [1] = 'R';
        if (StringICmp (name, gbq->val) == 0) {
          return TRUE;
        }
      }
      if (rrp->type == 255 || rrp->type == 8 || rrp->type == 9 || rrp->type == 10) {
        /* new convention follows ASN.1 spec comments, allows new RNA types */
        return FALSE;
      }
      /* subsequent /product now added to comment */
      if (sfp->comment == NULL) {
        sfp->comment = gbq->val;
        gbq->val = NULL;
      } else if (StringStr (gbq->val, sfp->comment) == NULL) {
        len = StringLen (sfp->comment) + StringLen (gbq->val) + 5;
        str = MemNew (sizeof (Char) * len);
        StringCpy (str, sfp->comment);
        StringCat (str, "; ");
        StringCat (str, gbq->val);
        sfp->comment = MemFree (sfp->comment);
        sfp->comment = str;
      }
      /* return FALSE; */
      return TRUE;
    }
    if (rrp->type == 8 || rrp->type == 9 || rrp->type == 10) {
      /* new convention follows ASN.1 spec comments, allows new RNA types */
      return FALSE;
    }
    if (rrp->ext.choice == 1 && rrp->ext.value.ptrvalue != NULL) {
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
    }
    if (rrp->ext.choice == 0 || rrp->ext.choice == 1) {
      rrp->ext.choice = 1;
      rrp->ext.value.ptrvalue = StringSave (gbq->val);
      return TRUE;
    }
  } else if (StringICmp (gbq->qual, "anticodon") == 0) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rrp == NULL) return FALSE;
    if (rrp->type == 0) {
      rrp->type = 255;
    }
    if (rrp->type == 3 && rrp->ext.choice == 0) {
      trp = (tRNAPtr) MemNew (sizeof (tRNA));
      if (trp != NULL) {
        rrp->ext.choice = 2;
        rrp->ext.value.ptrvalue = trp;
        for (j = 0; j < 6; j++) {
          trp->codon [j] = 255;
        }
      }
    }
    if (rrp->type == 3 && rrp->ext.choice == 2) {
      trp = (tRNAPtr) rrp->ext.value.ptrvalue;
      if (trp != NULL) {
        StringNCpy_0 (tmp, gbq->val, sizeof (tmp));
        ptr = StringStr (tmp, "(");
        if (ptr != NULL) {
          ptr = StringStr (ptr + 1, "pos");
          if (ptr != NULL) {
            ptr = StringStr (ptr + 3, ":");
          }
        }
        if (ptr != NULL) {
          str = ptr + 1;
          ptr = StringStr (str, "..");
          if (ptr != NULL) {
            *ptr = '\0';
            if (sscanf (str, "%ld", &val) == 1) {
              from = val - 1;
              str = ptr + 2;
              ptr = StringStr (str, ",");
              if (ptr != NULL) {
                *ptr = '\0';
                if (sscanf (str, "%ld", &val) == 1) {
                  to = val - 1;
                  sip = SeqLocId (sfp->location);
                  if (sip != NULL) {
                    bsp = BioseqFind (sip);
                    if (bsp != NULL) {
                      if (from >= 0 && from < bsp->length - 1) {
                        if (to >= 0 && to < bsp->length - 1) {
                          sintp = SeqIntNew ();
                          if (sintp != NULL) {
                            if (from > to) {
                              sintp->from = to;
                              sintp->to = from;
                              sintp->strand = Seq_strand_minus;
                            } else {
                              sintp->from = from;
                              sintp->to = to;
                              sintp->strand = Seq_strand_plus;
                            }
                            sintp->id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
                            trp->anticodon = ValNodeAddPointer (NULL, SEQLOC_INT, (Pointer) sintp);
                            if (trp->aatype == 0 && trp->aa == 0) {
                              ptr = StringStr (ptr + 1, "aa:");
                              if (ptr != NULL) {
                                str = ptr + 3;
                                ptr = StringStr (str, ")");
                                if (ptr != NULL) {
                                  *ptr = '\0';
                                  trp->aa = ParseTRnaString (str, NULL, NULL, FALSE);
                                  if (trp->aa != 0) {
                                    trp->aatype = 2;
                                  }
                                }
                              }
                            }
                            return TRUE;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return FALSE;
}

static Boolean HandledGBQualOnProt (SeqFeatPtr sfp, GBQualPtr gbq)

{
  Int2        choice = 0;
  ProtRefPtr  prp;
  ValNodePtr  vnp;

  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL) return FALSE;
  if (StringICmp (gbq->qual, "product") == 0) {
    choice = 1;
  } else if (StringICmp (gbq->qual, "function") == 0) {
    choice = 2;
  } else if (StringICmp (gbq->qual, "EC_number") == 0) {
    choice = 3;
  } else if (StringICmp (gbq->qual, "standard_name") == 0) {
    choice = 4;
  } else if (StringICmp (gbq->qual, "label") == 0) {
    choice = 5;
  } else if (StringICmp (gbq->qual, "allele") == 0) {
      choice = 6;
  }
  if (choice == 1 || choice == 4) {
    vnp = prp->name;
    if (vnp != NULL && (! HasNoText (vnp->data.ptrvalue))) return FALSE;
    ValNodeCopyStr (&(prp->name), 0, gbq->val);
    vnp = prp->name;
    if (vnp != NULL && prp->desc != NULL) {
      if (StringICmp (vnp->data.ptrvalue, prp->desc) == 0) {
        prp->desc = MemFree (prp->desc);
      }
    }
    return TRUE;
  } else if (choice == 2) {
    ValNodeCopyStr (&(prp->activity), 0, gbq->val);
    return TRUE;
  } else if (choice == 3) {
    ValNodeCopyStr (&(prp->ec), 0, gbq->val);
    return TRUE;
  } else if (choice == 5) {
    return FALSE; /* keep label gbqual only */
  } else if (choice == 6) {
      return FALSE;
  }

  if (StringICmp (gbq->qual, "experiment") == 0 ||
      StringICmp (gbq->qual, "inference") == 0) {
    return FALSE;
  }

  if (StringICmp (gbq->qual, "UniProtKB_evidence") == 0) {
    return FALSE;
  }

  return TRUE; /* all other gbquals not appropriate on protein features */
}

static Boolean HandledGBQualOnImp (SeqFeatPtr sfp, GBQualPtr gbq)

{
  Char        ch;
  ImpFeatPtr  ifp;
  Int4        len;
  CharPtr     ptr;

  if (StringICmp (gbq->qual, "rpt_unit") == 0) {
    if (HasNoText (gbq->val)) return FALSE;
    ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
    if (ifp == NULL) return FALSE;
    if (StringICmp (ifp->key, "repeat_region") != 0) return FALSE;
    len = SeqLocLen (sfp->location);
    if (len != (Int4) StringLen (gbq->val)) return FALSE;
    ptr = gbq->val;
    ch = *ptr;
    while (ch != '\0') {
      if (StringChr ("ACGTNacgtn", ch) == NULL) return FALSE;
      ptr++;
      ch = *ptr;
    }
    /* return TRUE; */
  }
  return FALSE;
}

static void CleanupRptUnit (GBQualPtr gbq)

{
  Char     ch;
  size_t   len;
  CharPtr  ptr;
  CharPtr  str;
  CharPtr  tmp;

  if (gbq == NULL) return;
  if (StringHasNoText (gbq->val)) return;
  len = StringLen (gbq->val) * 2 + 1;
  str = MemNew (sizeof (Char) * len);
  if (str == NULL) return;
  ptr = str;
  tmp = gbq->val;
  ch = *tmp;
  while (ch != '\0') {
    while (ch == '(' || ch == ')' || ch == ',') {
      *ptr = ch;
      ptr++;
      tmp++;
      ch = *tmp;
    }
    while (IS_WHITESP (ch)) {
      tmp++;
      ch = *tmp;
    }
    while (IS_DIGIT (ch)) {
      *ptr = ch;
      ptr++;
      tmp++;
      ch = *tmp;
    }
    if (ch == '.' || ch == '-') {
      while (ch == '.' || ch == '-') {
        tmp++;
        ch = *tmp;
      }
      *ptr = '.';
      ptr++;
      *ptr = '.';
      ptr++;
    }
    while (IS_WHITESP (ch)) {
      tmp++;
      ch = *tmp;
    }
    while (IS_DIGIT (ch)) {
      *ptr = ch;
      ptr++;
      tmp++;
      ch = *tmp;
    }
    while (IS_WHITESP (ch)) {
      tmp++;
      ch = *tmp;
    }
    if (ch == '\0' || ch == '(' || ch == ')' || ch == ',' || ch == '.' || IS_WHITESP (ch) || IS_DIGIT (ch)) {
    } else {
      MemFree (str);
      /* lower case the contents */
      ptr = gbq->val;
      ch = *ptr;
      while (ch != '\0') {
        if (IS_UPPER (ch)) {
          *ptr = TO_LOWER (ch);
        }
        ptr++;
        ch = *ptr;
      }
      return;
    }
  }
  *ptr = '\0';
  gbq->val = MemFree (gbq->val);
  gbq->val = str;
  /* and lower case the contents */
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_UPPER (ch)) {
      *ptr = TO_LOWER (ch);
    }
    ptr++;
    ch = *ptr;
  }
}

static void CleanupRptUnitSeq (GBQualPtr gbq)

{
  Char     ch;
  CharPtr  ptr;

  if (gbq == NULL) return;
  if (StringHasNoText (gbq->val)) return;
  /* lower case, and convert U to T */
  ptr = gbq->val;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_UPPER (ch)) {
      ch = TO_LOWER (ch);
      if (ch == 'u') {
        ch = 't';
      }
      *ptr = ch;
    }
    ptr++;
    ch = *ptr;
  }
}

static void CleanupRptUnitRange (GBQualPtr gbq)

{
  Char     ch;
  Int2     dashes = 0;
  Int2     dots = 0;
  size_t   len;
  CharPtr  ptr;
  CharPtr  str;
  CharPtr  tmp;

  if (gbq == NULL) return;
  if (StringHasNoText (gbq->val)) return;
  ptr = gbq->val;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == '-') {
      dashes++;
    } else if (ch == '.') {
      dots++;
    } else if (IS_DIGIT (ch)) {
      /* okay */
    } else return;
    ptr++;
    ch = *ptr;
  }

  if (dashes > 0 && dots == 0) {
    len = StringLen (gbq->val + dashes);
    str = (CharPtr) MemNew (sizeof (Char) * (len + 5));
    tmp = str;
    ptr = gbq->val;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == '-') {
        *tmp = '.';
        tmp++;
        *tmp = '.';
        tmp++;
      } else {
        *tmp = ch;
        tmp++;
      }
      ptr++;
      ch = *ptr;
    }
    gbq->val = MemFree (gbq->val);
    gbq->val = str;
  }
}

static void CleanupReplace (GBQualPtr gbq)

{
  Char     ch;
  CharPtr  ptr;

  if (gbq == NULL) return;
  if (StringHasNoText (gbq->val)) return;
  ptr = gbq->val;
  ch = *ptr;
  while (ch != '\0') {
    if (StringChr ("ACGTUacgtu", ch) == NULL) return;
    ptr++;
    ch = *ptr;
  }
  /* lower case, and convert U to T */
  ptr = gbq->val;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_UPPER (ch)) {
      ch = TO_LOWER (ch);
      if (ch == 'u') {
        ch = 't';
      }
      *ptr = ch;
    }
    ptr++;
    ch = *ptr;
  }
}

static CharPtr evCategoryPfx [] = {
  "",
  "COORDINATES: ",
  "DESCRIPTION: ",
  "EXISTENCE: ",
  NULL
};

static void CleanupInference (GBQualPtr gbq)

{
  Char     ch;
  CharPtr  colon;
  CharPtr  dst;
  Int2     j;
  size_t   len;
  CharPtr  ptr;
  CharPtr  skip;
  CharPtr  space;
  CharPtr  str;

  if (gbq == NULL) return;
  if (StringHasNoText (gbq->val)) return;

  str = gbq->val;
  space = NULL;
  colon = NULL;

  skip = NULL;
  for (j = 0; evCategoryPfx [j] != NULL; j++) {
    len = StringLen (evCategoryPfx [j]);
    if (StringNICmp (str, evCategoryPfx [j], len) != 0) continue;
    skip = str + len;
  }
  if (skip != NULL) {
    str = skip;
  }

  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    *dst = ch;
    if (ch == ' ') {
      if (space == NULL) {
        space = dst;
      }
    } else if (ch == ':') {
      if (space != NULL) {
        dst = space;
        *dst = ch;
      }
      space = NULL;
      colon = dst;
    } else {
      if (space != NULL && colon != NULL) {
        colon++;
        dst = colon;
        *dst = ch;
      }
      space = NULL;
      colon = NULL;
    }
    dst++;
    ptr++;
    ch = *ptr;
  }
  *dst = '\0';
}

static CharPtr evCategoryNoSpace [] = {
  "",
  "COORDINATES:",
  "DESCRIPTION:",
  "EXISTENCE:",
  NULL
};

static void RepairInference (GBQualPtr gbq)

{
  Int2     j;
  size_t   len;
  CharPtr  ptr;
  CharPtr  skip;
  CharPtr  str;

  if (gbq == NULL) return;
  if (StringHasNoText (gbq->val)) return;

  str = gbq->val;
  for (j = 0; evCategoryNoSpace [j] != NULL; j++) {
    len = StringLen (evCategoryNoSpace [j]);
    if (StringNICmp (str, evCategoryNoSpace [j], len) != 0) continue;
    if (StringNICmp (str, evCategoryPfx [j], len + 1) == 0) continue;
    /* need to repair */
    skip = str + len;
    ptr = MemNew (StringLen (skip) + 20);
    if (ptr == NULL) return;
    StringCpy (ptr, evCategoryPfx [j]);
    StringCat (ptr, skip);
    gbq->val = MemFree (gbq->val);
    gbq->val = ptr;
    return;
  }
}

static void CleanupConsSplice (GBQualPtr gbq)

{
  size_t   len;
  CharPtr  ptr;
  CharPtr  str;

  if (StringNICmp (gbq->val, "(5'site:", 8) != 0) return;
  ptr = StringStr (gbq->val, ",3'site:");
  if (ptr == NULL) return;
  len = StringLen (gbq->val) + 5;
  str = (CharPtr) MemNew (len);
  if (str == NULL) return;
  *ptr = '\0';
  ptr++;
  StringCpy (str, gbq->val);
  StringCat (str, ", ");
  StringCat (str, ptr);
  gbq->val = MemFree (gbq->val);
  gbq->val = str;
}

static Boolean ExpandParenGroup (GBQualPtr headgbq)

{
  Char       ch;
  GBQualPtr  lastgbq;
  size_t     len;
  Int2       nesting;
  GBQualPtr  newgbq;
  GBQualPtr  nextqual;
  CharPtr    ptr;
  CharPtr    str;
  CharPtr    tmp;

  nextqual = headgbq->next;
  lastgbq = headgbq;
  ptr = headgbq->val;
  tmp = StringSave (ptr + 1);
  len = StringLen (tmp);
  if (len > 0 && tmp [len - 1] == ')') {
    tmp [len - 1] = '\0';
  }
  str = tmp;
  nesting = 0;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == '(') {
      nesting++;
    } else if (ch == ')') {
      nesting--;
      if (nesting < 0) {
        MemFree (tmp);
        return FALSE;
      }
    } else if (ch == ',') {
      if (nesting < 0) {
        MemFree (tmp);
        return FALSE;
      }
    }
    ptr++;
    ch = *ptr;
  }
  while (! StringHasNoText (str)) {
    ptr = StringChr (str, ',');
    if (ptr == NULL) {
      ptr = StringRChr (str, ')');
    }
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }
    TrimSpacesAroundString (str);
    newgbq = GBQualNew ();
    if (newgbq != NULL) {
      newgbq->qual = StringSave (headgbq->qual);
      newgbq->val = StringSave (str);
      newgbq->next = nextqual;
      lastgbq->next = newgbq;
      lastgbq = newgbq;
    }
    str = ptr;
  }
  MemFree (tmp);
  return TRUE;
}

static Boolean IsBaseRange (CharPtr str)

{
  CharPtr   ptr;
  Char      tmp [32];
  long int  val;

  if (StringLen (str) > 25) return FALSE;
  StringNCpy_0 (tmp, str, sizeof (tmp));
  ptr = StringStr (tmp, "..");
  if (ptr == NULL) return FALSE;
  *ptr = '\0';
  if (StringHasNoText (tmp)) return FALSE;
  if (sscanf (tmp, "%ld", &val) != 1 || val < 1) return FALSE;
  ptr += 2;
  if (StringHasNoText (ptr)) return FALSE;
  if (sscanf (ptr, "%ld", &val) != 1 || val < 1) return FALSE;
  return TRUE;
}

static void ModernizeFeatureGBQuals (SeqFeatPtr sfp)

{
  GBQualPtr       gbq;
  size_t          len;
  GBQualPtr       nextqual;
  GBQualPtr PNTR  prevqual;
  CharPtr         str;
  Boolean         unlink;

  if (sfp == NULL) return;
  gbq = sfp->qual;
  prevqual = (GBQualPtr PNTR) &(sfp->qual);
  while (gbq != NULL) {
    CleanVisString (&(gbq->qual));
    CleanVisString (&(gbq->val));
    if (gbq->qual == NULL) {
      gbq->qual = StringSave ("");
    }
    if (StringIsJustQuotes (gbq->val)) {
      gbq->val = MemFree (gbq->val);
    }
    if (gbq->val == NULL) {
      gbq->val = StringSave ("");
    }
    nextqual = gbq->next;
    unlink = TRUE;
    if (StringICmp (gbq->qual, "rpt_unit_seq") == 0) {
      str = gbq->val;
      len = StringLen (str);
      if (len > 1 && *str == '{' && str [len - 1] == '}') {
        *str = '(';
        str [len - 1] = ')';
      }
      if (len > 1 && *str == '(' && str [len - 1] == ')' /* && StringChr (str + 1, '(') == NULL */) {
        if (ExpandParenGroup (gbq)) {
          nextqual = gbq->next;
          /* individual parsed out (xxx,xxx) qualifiers will be processed next, now get rid of original */
          unlink = TRUE;
        } else {
          unlink = FALSE;
        }
      } else {
        unlink = FALSE;
      }
    } else if (StringICmp (gbq->qual, "rpt_type") == 0 ||
        StringICmp (gbq->qual, "rpt_unit") == 0 ||
        StringICmp (gbq->qual, "rpt_unit_range") == 0 ||
        StringICmp (gbq->qual, "rpt_unit_seq") == 0 ||
        StringICmp (gbq->qual, "replace") == 0 ||
        StringICmp (gbq->qual, "compare") == 0 ||
        StringICmp (gbq->qual, "old_locus_tag") == 0 ||
        StringICmp (gbq->qual, "usedin") == 0) {
      str = gbq->val;
      len = StringLen (str);
      if (len > 1 && *str == '{' && str [len - 1] == '}') {
        *str = '(';
        str [len - 1] = ')';
      }
      if (len > 1 && *str == '(' && str [len - 1] == ')' && StringChr (str + 1, '(') == NULL) {
        if (ExpandParenGroup (gbq)) {
          nextqual = gbq->next;
          /* individual parsed out (xxx,xxx) qualifiers will be processed next, now get rid of original */
          unlink = TRUE;
        } else {
          unlink = FALSE;
        }
      } else {
        unlink = FALSE;
      }
    } else {
      unlink = FALSE;
    }
    if (unlink) {
      *(prevqual) = gbq->next;
      gbq->next = NULL;
      gbq->qual = MemFree (gbq->qual);
      gbq->val = MemFree (gbq->val);
      GBQualFree (gbq);
    } else {
      prevqual = (GBQualPtr PNTR) &(gbq->next);
    }
    gbq = nextqual;
  }
}


static void MendSatelliteQualifier (CharPtr PNTR satellite)
{
  Int4 microsatellite_len = StringLen ("microsatellite");
  Int4 minisatellite_len = StringLen ("minisatellite");
  Int4 satellite_len = StringLen ("satellite");
  Int4 type_len = 0;
  CharPtr new_qual, colon, src, dst;

  if (satellite == NULL || StringHasNoText (*satellite)) {
    return;
  }

  if (StringNCmp (*satellite, "microsatellite", microsatellite_len) == 0) {
    type_len = microsatellite_len;
  } else if (StringNCmp (*satellite, "minisatellite", minisatellite_len) == 0) {
    type_len = minisatellite_len;
  } else if (StringNCmp (*satellite, "satellite", satellite_len) == 0) {
    type_len = satellite_len;
  }
  
  if (type_len == 0) {
    new_qual = (CharPtr) MemNew (sizeof (Char) * (StringLen (*satellite) + satellite_len + 3));
    sprintf (new_qual, "satellite:%s", *satellite);
    *satellite = MemFree (*satellite);
    *satellite = new_qual;
  } else if (*(*satellite + type_len) == ' ') {
    *(*satellite + type_len) = ':';
  }

  /* remove spaces after colon */
  colon = StringChr (*satellite, ':');
  if (colon != NULL) {
    src = colon + 1;
    dst = colon + 1;
    while (*src == ' ') {
      src++;
    }
    while (*src != 0) {
      *dst = *src;
      dst++;
      src++;
    }
    *dst = 0;
  }
}


static void CleanupFeatureGBQuals (SeqFeatPtr sfp, Boolean isEmblOrDdbj)

{
  ValNodePtr      afterMe = NULL;
  Boolean         all_digits;
  Char            ch;
  DbtagPtr        db;
  GBQualPtr       gbq;
  GeneRefPtr      grp;
  ImpFeatPtr      ifp;
  size_t          len;
  GBQualPtr       nextqual;
  ObjectIdPtr     oip;
  GBQualPtr PNTR  prevqual;
  CharPtr         ptr;
  GBQualPtr       rpt_unit_range = NULL;
  GBQualPtr       rpt_unit_seq = NULL;
  CharPtr         str;
  CharPtr         tag;
  Boolean         unlink;
  ValNodePtr      vnp;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL) return;
  gbq = sfp->qual;
  prevqual = (GBQualPtr PNTR) &(sfp->qual);
  while (gbq != NULL) {
    CleanVisString (&(gbq->qual));
    CleanVisString (&(gbq->val));
    if (gbq->qual == NULL) {
      gbq->qual = StringSave ("");
    }
    if (StringIsJustQuotes (gbq->val)) {
      gbq->val = MemFree (gbq->val);
    }
    if (gbq->val == NULL) {
      gbq->val = StringSave ("");
    }
    if (StringICmp (gbq->qual, "replace") == 0) {
      if (sfp->data.choice == SEQFEAT_IMP) {
        ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
        if (ifp != NULL) {
          if (StringICmp (ifp->key, "variation") == 0 && gbq->val != NULL) {
            ptr = gbq->val;
            ch = *ptr;
            while (ch != '\0') {
              *ptr = TO_LOWER (ch);
              ptr++;
              ch = *ptr;
            }
          }
        }
      }
    }
    nextqual = gbq->next;
    unlink = TRUE;
    if (StringICmp (gbq->qual, "partial") == 0) {
      sfp->partial = TRUE;
    } else if (StringICmp (gbq->qual, "evidence") == 0) {
      /*
      if (StringICmp (gbq->val, "experimental") == 0) {
        if (sfp->exp_ev != 2) {
          sfp->exp_ev = 1;
        }
      } else if (StringICmp (gbq->val, "not_experimental") == 0) {
        sfp->exp_ev = 2;
      }
      */
    } else if (StringICmp (gbq->qual, "exception") == 0) {
      sfp->excpt = TRUE;
      if (! HasNoText (gbq->val)) {
        if (StringICmp (gbq->val, "TRUE") != 0) {
          if (sfp->except_text == NULL) {
            sfp->except_text = StringSaveNoNull (gbq->val);
          }
        }
      }
    } else if (StringICmp (gbq->qual, "note") == 0 ||
               StringICmp (gbq->qual, "notes") == 0 ||
               StringICmp (gbq->qual, "comment") == 0) {
      if (sfp->comment == NULL) {
        sfp->comment = gbq->val;
        gbq->val = NULL;
      } else {
        len = StringLen (sfp->comment) + StringLen (gbq->val) + 5;
        str = MemNew (sizeof (Char) * len);
        StringCpy (str, sfp->comment);
        StringCat (str, "; ");
        StringCat (str, gbq->val);
        sfp->comment = MemFree (sfp->comment);
        sfp->comment = str;
      }
    } else if (StringICmp (gbq->qual, "label") == 0) {
      if (StringICmp (gbq->val, FindKeyFromFeatDefType (sfp->idx.subtype, FALSE)) == 0) {
        /* skip label that is simply the feature key */
      } else if (sfp->comment == NULL || StringISearch (sfp->comment, gbq->qual) == NULL) {
        /* if label is not already in comment, append */
        len = StringLen (sfp->comment) + StringLen (gbq->val) + StringLen ("label: ") + 5;
        str = MemNew (sizeof (Char) * len);
        if (sfp->comment == NULL) {
          StringCpy (str, "label: ");
          StringCat (str, gbq->val);
          sfp->comment = str;
        } else {
          StringCpy (str, sfp->comment);
          StringCat (str, "; ");
          StringCat (str, "label: ");
          StringCat (str, gbq->val);
          sfp->comment = MemFree (sfp->comment);
          sfp->comment = str;
        }
      }
    } else if (StringICmp (gbq->qual, "db_xref") == 0) {
      tag = gbq->val;
      ptr = StringChr (tag, ':');
      if (ptr != NULL) {
        vnp = ValNodeNew (NULL);
        db = DbtagNew ();
        vnp->data.ptrvalue = db;
        *ptr = '\0';
        ptr++;
        db->db = StringSave (tag);
        oip = ObjectIdNew ();
        oip->str = StringSave (ptr);
        db->tag = oip;
        vnp->next = sfp->dbxref;
        sfp->dbxref = vnp;
      } else {
        /*
        db->db = StringSave ("?");
        oip = ObjectIdNew ();
        oip->str = StringSave (tag);
        db->tag = oip;
        vnp->next = sfp->dbxref;
        sfp->dbxref = vnp;
        */
        unlink = FALSE;
      }
    } else if (StringICmp (gbq->qual, "gdb_xref") == 0) {
      vnp = ValNodeNew (NULL);
      db = DbtagNew ();
      vnp->data.ptrvalue = db;
      db->db = StringSave ("GDB");
      oip = ObjectIdNew ();
      oip->str = StringSave (gbq->val);
      db->tag = oip;
      vnp->next = sfp->dbxref;
      sfp->dbxref = vnp;
    } else if (StringICmp (gbq->qual, "cons_splice") == 0) {
      /*
      CleanupConsSplice (gbq);
      unlink = FALSE;
      */
    } else if (StringICmp (gbq->qual, "replace") == 0) {
      CleanupReplace (gbq);
      unlink = FALSE;
    } else if (StringICmp (gbq->qual, "rpt_unit_seq") == 0) {
      CleanupRptUnitSeq (gbq);
      unlink = FALSE;
    } else if (StringICmp (gbq->qual, "rpt_unit_range") == 0) {
      CleanupRptUnitRange (gbq);
      unlink = FALSE;
    } else if (sfp->data.choice == SEQFEAT_GENE && HandledGBQualOnGene (sfp, gbq)) {
    } else if (sfp->data.choice == SEQFEAT_CDREGION && HandledGBQualOnCDS (sfp, gbq, &afterMe)) {
    } else if (sfp->data.choice == SEQFEAT_RNA && HandledGBQualOnRNA (sfp, gbq, isEmblOrDdbj)) {
    } else if (sfp->data.choice == SEQFEAT_PROT && HandledGBQualOnProt (sfp, gbq)) {
    } else if (sfp->data.choice == SEQFEAT_IMP && HandledGBQualOnImp (sfp, gbq)) {
    } else if (StringICmp (gbq->qual, "rpt_unit") == 0) {
      if (IsBaseRange (gbq->val)) {
        gbq->qual = MemFree (gbq->qual);
        gbq->qual = StringSave ("rpt_unit_range");
        unlink = FALSE;
      } else {
        gbq->qual = MemFree (gbq->qual);
        gbq->qual = StringSave ("rpt_unit_seq");
        unlink = FALSE;
      }
    } else if (StringICmp (gbq->qual, "EC_number") == 0) {
      CleanupECNumber (gbq->val);
      unlink = FALSE;
    } else if (StringICmp (gbq->qual, "pseudo") == 0) {
      sfp->pseudo = TRUE;
    } else if (StringICmp (gbq->qual, "ribosomal_slippage") == 0 ||
               StringICmp (gbq->qual, "ribosomal-slippage") == 0 ||
               StringICmp (gbq->qual, "ribosomal slippage") == 0) {
      sfp->excpt = TRUE;
      if (HasNoText (gbq->val)) {
        if (sfp->except_text == NULL) {
          sfp->except_text = StringSaveNoNull ("ribosomal slippage");
        }
      }
    } else if (StringICmp (gbq->qual, "trans_splicing") == 0 ||
               StringICmp (gbq->qual, "trans-splicing") == 0 ||
               StringICmp (gbq->qual, "trans splicing") == 0) {
      sfp->excpt = TRUE;
      if (HasNoText (gbq->val)) {
        if (sfp->except_text == NULL) {
          sfp->except_text = StringSaveNoNull ("trans-splicing");
        }
      }
    } else if (StringICmp (gbq->qual, "artificial_location") == 0 ||
               StringICmp (gbq->qual, "artificial-location") == 0 ||
               StringICmp (gbq->qual, "artificial location") == 0) {
      sfp->excpt = TRUE;
      if (HasNoText (gbq->val)) {
        if (sfp->except_text == NULL) {
          sfp->except_text = StringSaveNoNull ("artificial location");
        }
      }
    } else if (StringICmp (gbq->qual, "gene") == 0 && (! StringHasNoText (gbq->val))) {
      grp = GeneRefNew ();
      grp->locus = StringSave (gbq->val);
      xref = SeqFeatXrefNew ();
      xref->data.choice = SEQFEAT_GENE;
      xref->data.value.ptrvalue = (Pointer) grp;
      xref->specialCleanupFlag = TRUE; /* flag to test for overlapping gene later */
      xref->next = sfp->xref;
      sfp->xref = xref;
    } else if (sfp->data.choice != SEQFEAT_CDREGION && StringICmp (gbq->qual, "codon_start") == 0) {
      /* not legal on anything but CDS, so remove it */
    } else if (StringICmp (gbq->qual, "experiment") == 0 &&
               StringICmp (gbq->val, "experimental evidence, no additional details recorded") == 0) {
      /* remove default experiment string if instantiated */
    } else if (StringICmp (gbq->qual, "inference") == 0) {
      if (StringICmp (gbq->val, "non-experimental evidence, no additional details recorded") == 0) {
        /* remove default inference string if instantiated */
      } else {
        CleanupInference (gbq);
        RepairInference (gbq);
        unlink = FALSE;
      }
    } else if (StringICmp (gbq->qual, "transposon") == 0) {
      if (StringICmp (gbq->val, "class I integron") == 0 ||
          StringICmp (gbq->val, "class II integron") == 0 ||
          StringICmp (gbq->val, "class III integron") == 0 ||
          StringICmp (gbq->val, "class 1 integron") == 0 ||
          StringICmp (gbq->val, "class 2 integron") == 0 ||
          StringICmp (gbq->val, "class 3 integron") == 0) {
        len = StringLen ("integron") + StringLen (gbq->val) + 5;
        str = MemNew (sizeof (Char) * len);
        StringCpy (str, "integron");
        StringCat (str, ":");
        ptr = StringStr (gbq->val, " integron");
        if (ptr != NULL) {
          *ptr = '\0';
        }
        StringCat (str, gbq->val);
        gbq->val = MemFree (gbq->val);
        gbq->val = str;
        gbq->qual = MemFree (gbq->qual);
        gbq->qual = StringSave ("mobile_element");
        unlink = FALSE;
      } else {
        len = StringLen ("transposon") + StringLen (gbq->val) + 5;
        str = MemNew (sizeof (Char) * len);
        StringCpy (str, "transposon");
        StringCat (str, ":");
        StringCat (str, gbq->val);
        gbq->val = MemFree (gbq->val);
        gbq->val = str;
        gbq->qual = MemFree (gbq->qual);
        gbq->qual = StringSave ("mobile_element");
        unlink = FALSE;
      }
    } else if (StringICmp (gbq->qual, "insertion_seq") == 0) {
      len = StringLen ("insertion sequence") + StringLen (gbq->val) + 5;
      str = MemNew (sizeof (Char) * len);
      StringCpy (str, "insertion sequence");
      StringCat (str, ":");
      StringCat (str, gbq->val);
      gbq->val = MemFree (gbq->val);
      gbq->val = str;
      gbq->qual = MemFree (gbq->qual);
      gbq->qual = StringSave ("mobile_element");
      unlink = FALSE;
    } else if (StringCmp (gbq->qual, "satellite") == 0) {
      MendSatelliteQualifier(&(gbq->val));
      unlink = FALSE;
    } else {
      unlink = FALSE;
    }
    if (StringICmp (gbq->qual, "mobile_element") == 0) {
      if (StringStr (gbq->val, " :") == 0 || StringStr (gbq->val, ": ") == 0) {
        len = StringLen (gbq->val) + 5;
        ptr = StringChr (gbq->val, ':');
        if (ptr != NULL) {
          *ptr = '\0';
          ptr++;
          TrimSpacesAroundString (gbq->val);
          TrimSpacesAroundString (ptr);
          str = MemNew (sizeof (Char) * len);
          StringCpy (str, gbq->val);
          StringCat (str, ":");
          StringCat (str, ptr);
          gbq->val = MemFree (gbq->val);
          gbq->val = str;
        }
      }
    }
    
    if (StringICmp (gbq->qual, "mobile_element") == 0) {
      if (sfp->data.choice == SEQFEAT_IMP) {
        ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
        if (ifp != NULL) {
          if (StringICmp (ifp->key, "repeat_region") == 0 && gbq->val != NULL) {
            gbq->qual = MemFree (gbq->qual);
            gbq->qual = StringSave ("mobile_element_type");
            ifp->key = MemFree (ifp->key);
            ifp->key = StringSave ("mobile_element");
            sfp->idx.subtype = FEATDEF_mobile_element;
          }
        }
      }
    }
    if (StringICmp (gbq->qual, "estimated_length") == 0) {
      all_digits = TRUE;
      ptr = gbq->val;
      ch = *ptr;
      while (ch != '\0') {
        if (! IS_DIGIT (ch)) {
          all_digits = FALSE;
        }
        ptr++;
        ch = *ptr;
      }
      if (! all_digits) {
        if (StringICmp (gbq->val, "unknown") != 0) {
          MemFree (gbq->val);
          gbq->val = StringSave ("unknown");
        }
      }
    }

    if (sfp->data.choice == SEQFEAT_IMP) {
      ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
      if (ifp != NULL) {
        if (StringICmp (ifp->key, "conflict") == 0 ) {
          ifp->key = MemFree (ifp->key);
          ifp->key = StringSave ("misc_difference");
          sfp->idx.subtype = FEATDEF_misc_difference;
          len = StringLen (sfp->comment) + StringLen ("conflict") + 5;
          str = MemNew (sizeof (Char) * len);
          if (sfp->comment == NULL) {
            StringCpy (str, "conflict");
            sfp->comment = str;
          } else {
            StringCpy (str, "conflict; ");
            StringCat (str, sfp->comment);
            sfp->comment = MemFree (sfp->comment);
            sfp->comment = str;
          }
        }
      }
    }

    if (rpt_unit_seq != NULL) {
      CleanupRptUnit (rpt_unit_seq);
    }
    if (rpt_unit_range != NULL) {
      CleanupRptUnit (rpt_unit_range);
    }
    if (unlink) {
      *(prevqual) = gbq->next;
      gbq->next = NULL;
      gbq->qual = MemFree (gbq->qual);
      gbq->val = MemFree (gbq->val);
      GBQualFree (gbq);
    } else {
      prevqual = (GBQualPtr PNTR) &(gbq->next);
    }
    gbq = nextqual;
  }
}

static int LIBCALLBACK SortByGBQualKeyAndVal (VoidPtr ptr1, VoidPtr ptr2)

{
  int        compare;
  GBQualPtr  gbq1;
  GBQualPtr  gbq2;
  CharPtr    str1;
  CharPtr    str2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  gbq1 = *((GBQualPtr PNTR) ptr1);
  gbq2 = *((GBQualPtr PNTR) ptr2);
  if (gbq1 == NULL || gbq2 == NULL) return 0;
  str1 = (CharPtr) gbq1->qual;
  str2 = (CharPtr) gbq2->qual;
  if (str1 == NULL || str2 == NULL) return 0;
  compare = StringICmp (str1, str2);
  if (compare != 0) return compare;
  str1 = (CharPtr) gbq1->val;
  str2 = (CharPtr) gbq2->val;
  if (str1 == NULL || str2 == NULL) return 0;
  compare = StringICmp (str1, str2);
  return compare;
}

static Boolean GBQualsAlreadyInOrder (GBQualPtr list)

{
  int        compare;
  GBQualPtr  curr;
  GBQualPtr  next;

  if (list == NULL || list->next == NULL) return TRUE;
  curr = list;
  next = curr->next;
  while (next != NULL) {
    compare = StringICmp (curr->qual, next->qual);
    if (compare > 0) return FALSE;
    if (compare == 0) {
      compare = StringICmp (curr->val, next->val);
      if (compare > 0) return FALSE;
    }
    curr = next;
    next = curr->next;
  }
  return TRUE;
}

static GBQualPtr SortFeatureGBQuals (GBQualPtr list)

{
  size_t     count, i;
  GBQualPtr  gbq, PNTR head;

  if (list == NULL) return NULL;
  if (GBQualsAlreadyInOrder (list)) return list;

  for (gbq = list, count = 0; gbq != NULL; gbq = gbq->next, count++) continue;
  head = MemNew (sizeof (GBQualPtr) * (count + 1));

  for (gbq = list, i = 0; gbq != NULL && i < count; i++) {
    head [i] = gbq;
    gbq = gbq->next;
  }

  StableMergeSort (head, count, sizeof (GBQualPtr), SortByGBQualKeyAndVal);

  for (i = 0; i < count; i++) {
    gbq = head [i];
    gbq->next = head [i + 1];
  }

  list = head [0];
  MemFree (head);

  return list;
}

static void CleanupDuplicateGBQuals (GBQualPtr PNTR prevgbq)

{
  GBQualPtr  gbq;
  GBQualPtr  last = NULL;
  GBQualPtr  next;
  Boolean    unlink;

  if (prevgbq == NULL) return;
  gbq = *prevgbq;
  while (gbq != NULL) {
    next = gbq->next;
    unlink = FALSE;
    if (last != NULL) {
      if (StringICmp (last->qual, gbq->qual) == 0 &&
          StringICmp (last->val, gbq->val) == 0) {
        unlink = TRUE;
      }
    } else {
      last = gbq;
    }
    if (unlink) {
      *prevgbq = gbq->next;
      gbq->next = NULL;
      GBQualFree (gbq);
    } else {
      last = gbq;
      prevgbq = (GBQualPtr PNTR) &(gbq->next);
    }
    gbq = next;
  }
}

/* this identifies gbquals that should have been placed into special fields */

#define NUM_ILLEGAL_QUALS 14

/* StringICmp use of TO_UPPER means translation should go before transl_XXX */

static CharPtr illegalGbqualList [NUM_ILLEGAL_QUALS] = {
  "anticodon",
  "citation",
  "codon_start",
  "db_xref",
  "evidence",
  "exception",
  "gene",
  "note",
  "protein_id",
  "pseudo",
  "transcript_id",
  "transl_except",
  "transl_table",
  "translation",
};

static Int2 QualifierIsIllegal (CharPtr qualname)

{
  Int2  L, R, mid;

  if (qualname == NULL || *qualname == '\0') return FALSE;

  L = 0;
  R = NUM_ILLEGAL_QUALS - 1;

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (illegalGbqualList [mid], qualname) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (illegalGbqualList [R], qualname) == 0) {
    return TRUE;
  }

  return FALSE;
}

static void GbqualLink (GBQualPtr PNTR head, GBQualPtr qual)

{
  GBQualPtr  gbq;

  if (head == NULL || qual == NULL) return;
  gbq = *head;
  if (gbq != NULL) {
    while (gbq->next != NULL) {
      gbq = gbq->next;
    }
    gbq->next = qual;
  } else {
    *head = qual;
  }
}

static GBQualPtr SortIllegalGBQuals (GBQualPtr list)

{
  GBQualPtr  gbq, next, legal = NULL, illegal = NULL;

  gbq = list;
  while (gbq != NULL) {
    next = gbq->next;
    gbq->next = NULL;
    if (QualifierIsIllegal (gbq->qual)) {
      GbqualLink (&illegal, gbq);
    } else {
      GbqualLink (&legal, gbq);
    }
    gbq = next;
  }
  GbqualLink (&legal, illegal);
  return legal;
}

static Boolean IsSubString (CharPtr str1, CharPtr str2)

{
  Char    ch;
  size_t  len1, len2;

  len1 = StringLen (str1);
  len2 = StringLen (str2);
  if (len1 >= len2) return FALSE;
  if (StringNICmp (str1, str2, len1) != 0) return FALSE;
  ch = str2 [len1];
  if (IS_ALPHANUM (ch)) return FALSE;
  return TRUE;
}

static int LIBCALLBACK SortByOrgModSubtype (VoidPtr ptr1, VoidPtr ptr2)

{
  int        compare;
  OrgModPtr  omp1;
  OrgModPtr  omp2;
  CharPtr    str1;
  CharPtr    str2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  omp1 = *((OrgModPtr PNTR) ptr1);
  omp2 = *((OrgModPtr PNTR) ptr2);
  if (omp1 == NULL || omp2 == NULL) return 0;
  if (omp1->subtype > omp2->subtype) {
    return 1;
  } else if (omp1->subtype < omp2->subtype) {
    return -1;
  }
  str1 = (CharPtr) omp1->subname;
  str2 = (CharPtr) omp2->subname;
  if (str1 == NULL || str2 == NULL) return 0;
  compare = StringICmp (str1, str2);
  return compare;
}

static Boolean OrgModsAlreadyInOrder (OrgModPtr list)

{
  int        compare;
  OrgModPtr  curr;
  OrgModPtr  next;
  CharPtr    str1;
  CharPtr    str2;

  if (list == NULL || list->next == NULL) return TRUE;
  curr = list;
  next = curr->next;
  while (next != NULL) {
    if (curr->subtype > next->subtype) return FALSE;
    str1 = (CharPtr) curr->subname;
    str2 = (CharPtr) next->subname;
    compare = StringICmp (str1, str2);
    if (compare > 0) return FALSE;
    curr = next;
    next = curr->next;
  }
  return TRUE;
}

static OrgModPtr SortOrgModList (OrgModPtr list)

{
  size_t     count, i;
  OrgModPtr  omp, PNTR head;

  if (list == NULL) return NULL;
  if (OrgModsAlreadyInOrder (list)) return list;

  for (omp = list, count = 0; omp != NULL; omp = omp->next, count++) continue;
  head = MemNew (sizeof (OrgModPtr) * (count + 1));

  for (omp = list, i = 0; omp != NULL && i < count; i++) {
    head [i] = omp;
    omp = omp->next;
  }

  StableMergeSort (head, count, sizeof (OrgModPtr), SortByOrgModSubtype);

  for (i = 0; i < count; i++) {
    omp = head [i];
    omp->next = head [i + 1];
  }

  list = head [0];
  MemFree (head);

  return list;
}


static void RemoveSpaceBeforeAndAfterColon (CharPtr str)
{
  CharPtr pColon, cp, src, dst;

  if (StringHasNoText (str)) {
    return;
  }

  pColon = StringChr (str, ':');
  while (pColon != NULL) {
    cp = pColon - 1;
    while (cp > str && isspace (*cp)) {
      cp--;
    }
    if (cp < str || !isspace (*cp)) {
      cp++;
    }
    *cp = ':';
    dst = cp + 1;
    cp = pColon + 1;
    while (isspace (*cp)) {
      cp++;
    }
    src = cp;
    pColon = dst - 1;
    if (src != dst) {
      while (*src != 0) {
        *dst = *src;
        dst++; src++;
      }
      *dst = 0;
    }
    pColon = StringChr (pColon + 1, ':');
  }
}


static void CleanOrgModListEx (OrgModPtr PNTR ompp, CharPtr orpcommon)

{
  Char            ch;
  OrgModPtr       last = NULL;
  OrgModPtr       next;
  OrgModPtr       omp;
  OrgModPtr       omp_anamorph, omp_gb_anamorph, omp_other;
  OrgModPtr PNTR  prev;
  CharPtr         ptr;
  Boolean         redund;
  CharPtr         str;
  CharPtr         tmp;
  Boolean         unlink;

  if (ompp == NULL) return;
  prev = ompp;
  omp = *ompp;
  while (omp != NULL) {
    next = omp->next;
    unlink= FALSE;
    CleanVisStringAndCompress (&(omp->subname));
    TrimSpacesAndJunkFromEnds (omp->subname, FALSE);
    RemoveFlankingQuotes (&(omp->subname));
    CleanVisStringAndCompress (&(omp->attrib));
    if (omp->subtype == ORGMOD_common && StringICmp (omp->subname, orpcommon) == 0) {
      unlink = TRUE;
    } else if (last != NULL) {
      if (HasNoText (omp->subname)) {
        unlink = TRUE;
      } else if (last->subtype == omp->subtype &&
                 StringICmp (last->subname, omp->subname) == 0 ||
                 (last->subtype == omp->subtype &&
                 last->subtype == ORGMOD_other &&
                  StringStr (last->subname, omp->subname) != NULL)) {
        unlink = TRUE;
      } else if (last->subtype == omp->subtype &&
                 last->subtype == ORGMOD_other &&
                 IsSubString (last->subname, omp->subname)) {
        last->subname = MemFree (last->subname);
        last->subname = omp->subname;
        omp->subname = NULL;
        unlink = TRUE;
      }
    } else if (HasNoText (omp->subname) ||
               StringCmp (omp->subname, ")") == 0 ||
               StringCmp (omp->subname, "(") == 0) {
      unlink = TRUE;
    } else {
      last = omp;
    }
    if (unlink) {
      *prev = omp->next;
      omp->next = NULL;
      OrgModFree (omp);
    } else {
      last = omp;
      prev = &(omp->next);
    }
    omp = next;
  }

  for (omp = *ompp; omp != NULL; omp = omp->next) {
    if (omp->subtype != ORGMOD_specimen_voucher &&
        omp->subtype != ORGMOD_culture_collection &&
        omp->subtype != ORGMOD_bio_material) continue;
    if (StringHasNoText (omp->subname)) continue;
    RemoveSpaceBeforeAndAfterColon (omp->subname);
    ptr = StringStr (omp->subname, "::");
    if (ptr == NULL) continue;
    ptr++;
    tmp = ptr;
    tmp++;
    ch = *tmp;
    while (ch != '\0') {
      *ptr = ch;
      ptr++;
      tmp++;
      ch = *tmp;
    }
    *ptr = '\0';
  }

  omp_anamorph = NULL;
  omp_gb_anamorph = NULL;
  omp_other = NULL;
  redund = FALSE;

  for (omp = *ompp; omp != NULL; omp = omp->next) {
    if (omp->subtype == ORGMOD_anamorph) {
      omp_anamorph = omp;
    } else if (omp->subtype == ORGMOD_gb_anamorph) {
      omp_gb_anamorph = omp;
    } else if (omp->subtype == ORGMOD_other) {
      omp_other = omp;
    }
  }
  if (omp_other != NULL && StringNICmp (omp_other->subname, "anamorph:", 9) == 0) {
    ptr = omp_other->subname + 9;
    ch = *ptr;
    while (ch == ' ') {
      ptr++;
      ch = *ptr;
    }
    if (omp_anamorph != NULL) {
      str = omp_anamorph->subname;
      if (StringCmp (ptr, str) == 0) {
        redund = TRUE;
      }
    } else if (omp_gb_anamorph != NULL) {
      str = omp_gb_anamorph->subname;
      if (StringCmp (ptr, str) == 0) {
        redund = TRUE;
      }
    }
  }
  if (redund) {
    prev = ompp;
    omp = *ompp;
    while (omp != NULL) {
      next = omp->next;
      unlink= FALSE;
      if (omp == omp_other) {
        unlink= TRUE;
      }
      if (unlink) {
        *prev = omp->next;
        omp->next = NULL;
        OrgModFree (omp);
      } else {
        prev = &(omp->next);
      }
      omp = next;
    }
  }
}

NLM_EXTERN void CleanOrgModList (OrgModPtr PNTR ompp)

{
  CleanOrgModListEx (ompp, NULL);
}

static Boolean IsNoNameSubSource (SubSourcePtr ssp)

{
  if (ssp == NULL) return FALSE;

  return (Boolean) (ssp->subtype == SUBSRC_germline ||
                    ssp->subtype == SUBSRC_rearranged ||
                    ssp->subtype == SUBSRC_transgenic ||
                    ssp->subtype == SUBSRC_environmental_sample ||
                    ssp->subtype == SUBSRC_metagenomic);
}

static int LIBCALLBACK SortBySubSourceSubtype (VoidPtr ptr1, VoidPtr ptr2)

{
  int           compare;
  SubSourcePtr  ssp1;
  SubSourcePtr  ssp2;
  CharPtr       str1;
  CharPtr       str2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  ssp1 = *((SubSourcePtr PNTR) ptr1);
  ssp2 = *((SubSourcePtr PNTR) ptr2);
  if (ssp1 == NULL || ssp2 == NULL) return 0;
  if (ssp1->subtype > ssp2->subtype) {
    return 1;
  } else if (ssp1->subtype < ssp2->subtype) {
    return -1;
  }
  if (IsNoNameSubSource (ssp1)) return 0;
  str1 = (CharPtr) ssp1->name;
  str2 = (CharPtr) ssp2->name;
  if (str1 == NULL || str2 == NULL) return 0;
  compare = StringICmp (str1, str2);
  return compare;
}

static Boolean SubSourceAlreadyInOrder (SubSourcePtr list)

{
  int           compare;
  SubSourcePtr  curr;
  SubSourcePtr  next;
  CharPtr       str1;
  CharPtr       str2;

  if (list == NULL || list->next == NULL) return TRUE;
  curr = list;
  next = curr->next;
  while (next != NULL) {
    if (curr->subtype > next->subtype) return FALSE;
    if (curr->subtype == next->subtype) {
      if (! IsNoNameSubSource (curr)) {
        str1 = (CharPtr) curr->name;
        str2 = (CharPtr) next->name;
        compare = StringICmp (str1, str2);
        if (compare > 0) return FALSE;
      }
    }
    curr = next;
    next = curr->next;
  }
  return TRUE;
}

static SubSourcePtr SortSubSourceList (SubSourcePtr list)

{
  size_t        count, i;
  SubSourcePtr  ssp, PNTR head;

  if (list == NULL) return NULL;
  if (SubSourceAlreadyInOrder (list)) return list;

  for (ssp = list, count = 0; ssp != NULL; ssp = ssp->next, count++) continue;
  head = MemNew (sizeof (SubSourcePtr) * (count + 1));

  for (ssp = list, i = 0; ssp != NULL && i < count; i++) {
    head [i] = ssp;
    ssp = ssp->next;
  }

  StableMergeSort (head, count, sizeof (SubSourcePtr), SortBySubSourceSubtype);

  for (i = 0; i < count; i++) {
    ssp = head [i];
    ssp->next = head [i + 1];
  }

  list = head [0];
  MemFree (head);

  return list;
}

static CharPtr TrimParenthesesAndCommasAroundString (CharPtr str)

{
  Uchar    ch;    /* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && (ch < ' ' || ch == '(' || ch == ',')) {
      ptr++;
      ch = *ptr;
    }
    while (ch != '\0') {
      *dst = ch;
      dst++;
      ptr++;
      ch = *ptr;
    }
    *dst = '\0';
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch != ')' && ch != ',') {
        dst = NULL;
      } else if (dst == NULL) {
        dst = ptr;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

static CharPtr CombineSplitQual (CharPtr origval, CharPtr newval)

{
  size_t   len;
  CharPtr  str = NULL;

  if (StringStr (origval, newval) != NULL) return origval;
  len = StringLen (origval) + StringLen (newval) + 5;
  str = MemNew (sizeof (Char) * len);
  if (str == NULL) return origval;
  TrimParenthesesAndCommasAroundString (origval);
  TrimParenthesesAndCommasAroundString (newval);
  StringCpy (str, "(");
  StringCat (str, origval);
  StringCat (str, ",");
  StringCat (str, newval);
  StringCat (str, ")");
  /* free original string, knowing return value will replace it */
  MemFree (origval);
  return str;
}

static Uint1 LocationForPlastidText (CharPtr plastid_name)
{
  if (StringICmp (plastid_name, "chloroplast") == 0) {
    return GENOME_chloroplast;
  } else if (StringICmp (plastid_name, "chromoplast") == 0) {
    return GENOME_chromoplast;
  } else if (StringICmp (plastid_name, "kinetoplast") == 0) {
    return GENOME_kinetoplast;
  } else if (StringICmp (plastid_name, "plastid") == 0) {
    return GENOME_plastid;
  } else if (StringICmp (plastid_name, "apicoplast") == 0) {
    return GENOME_apicoplast;
  } else if (StringICmp (plastid_name, "leucoplast") == 0) {
    return GENOME_leucoplast;
  } else if (StringICmp (plastid_name, "proplastid") == 0) {
    return GENOME_proplastid;
  } else if (StringICmp (plastid_name, "chromatophore") == 0) {
    return GENOME_chromatophore;
  } else {
    return 0;
  }
}

NLM_EXTERN void StringToLower (CharPtr str)

{
  Char  ch;

  if (str == NULL) return;
  ch = *str;
  while (ch != '\0') {
    *str = TO_LOWER (ch);
    str++;
    ch = *str;
  }
}


static void CleanPCRPrimerSeq (CharPtr seq)
{
  CharPtr ptr, src, dst, tmp;
  Char    ch;
  Boolean in_brackets = FALSE;
  Int4    i;

  if (StringHasNoText (seq)) {
    return;
  }

  /* upper case sequence */
  ptr = seq;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_UPPER (ch)) {
      *ptr = TO_LOWER (ch);
    }
    ptr++;
    ch = *ptr;
  }
  /* remove any spaces in sequence outisde of <modified base> */
  src = seq;
  dst = seq;
  ch = *src;
  while (ch != '\0') {
    if (ch == '<') {
      in_brackets = TRUE;
      *dst = ch;
      dst++;
    } else if (ch == '>') {
      in_brackets = FALSE;
      *dst = ch;
      dst++;
    } else if (ch != ' ') {
      *dst = ch;
      dst++;
    } else if (in_brackets) {
      *dst = ch;
      dst++;
    }
    src++;
    ch = *src;
  }
  *dst = '\0';
  /* upper case modified base <OTHER> */
  ptr = seq;
  tmp = StringStr (ptr, "<other>");
  while (tmp != NULL) {
    ptr = tmp + 7;
    for (i = 1; i < 6; i++) {
      ch = tmp [i];
      tmp [i] = TO_UPPER (ch);
    }
    tmp = StringStr (ptr, "<other>");
  }
}


static void CleanupPCRPrimers (PCRPrimerPtr PNTR pppp)

{
  PCRPrimerPtr       next;
  PCRPrimerPtr PNTR  prev;
  PCRPrimerPtr       ppp;

  if (pppp == NULL) return;

  prev = pppp;
  ppp = *pppp;
  while (ppp != NULL) {
    next = ppp->next;

    CleanVisString (&(ppp->seq));
    CleanPCRPrimerSeq (ppp->seq);
    CleanVisString (&(ppp->name));

    if (ppp->seq == NULL && ppp->name == NULL) {
      *prev = next;
      ppp->next = NULL;
      PCRPrimerFree (ppp);
    } else {
      StringToLower (ppp->seq);
      prev = &(ppp->next);
    }

    ppp = next;
  }

  /* fix artifact caused by fwd/rev-primer-seq starting with colon, separating name and seq */

  ppp = *pppp;
  if (ppp == NULL) return;
  next = ppp->next;
  if (next == NULL) return;
  if (next->next != NULL) return;

  if (ppp->name != NULL && ppp->seq == NULL && next->name == NULL && next->seq != NULL) {
    ppp->seq = next->seq;
    next->seq = NULL;
    ppp->next = NULL;
    PCRPrimerFree (next);
  } else if (ppp->seq != NULL && ppp->name == NULL && next->seq == NULL && next->name != NULL) {
    ppp->name = next->name;
    next->name = NULL;
    ppp->next = NULL;
    PCRPrimerFree (next);
  }
}

static void CleanupPCRReactionSet (PCRReactionSetPtr PNTR prpp)

{
  PCRReactionSetPtr       next;
  PCRReactionSetPtr PNTR  prev;
  PCRReactionSetPtr       prp;

  if (prpp == NULL) return;

  prev = prpp;
  prp = *prpp;
  while (prp != NULL) {
    next = prp->next;

    CleanupPCRPrimers (&(prp->forward));
    CleanupPCRPrimers (&(prp->reverse));

    if (prp->forward == NULL && prp->reverse == NULL) {
      *prev = next;
      prp->next = NULL;
      PCRReactionFree (prp);
    } else {
      prev = &(prp->next);
    }

    prp = next;
  }
}

extern void CleanSubSourceList (SubSourcePtr PNTR sspp, Uint1 location)

{
  Char               ch;
  CharPtr            dst;
  Int2               i;
  Boolean            in_brackets = FALSE;
  SubSourcePtr       last = NULL;
  SubSourcePtr       next;
  SubSourcePtr PNTR  prev;
  CharPtr            ptr;
  CharPtr            src;
  SubSourcePtr       ssp;
  CharPtr            str;
  CharPtr            tmp;
  Boolean            unlink;
  /*
  SubSourcePtr       fwd_seq = NULL, rev_seq = NULL, fwd_name = NULL, rev_name = NULL;
  size_t             len;
  */

  if (sspp == NULL) return;
  prev = sspp;
  ssp = *sspp;
  while (ssp != NULL) {
    next = ssp->next;
    unlink= FALSE;
    if (! IsNoNameSubSource (ssp)) {
      CleanVisStringAndCompress (&(ssp->name));
      TrimSpacesAndJunkFromEnds (ssp->name, FALSE);
      RemoveFlankingQuotes (&(ssp->name));
    } else /* if (StringICmp (ssp->name, "TRUE") == 0) */ {
      ssp->name = MemFree (ssp->name);
      ssp->name = StringSave ("");
    }
    if (ssp->subtype == SUBSRC_country) {
      CleanVisStringJunk (&(ssp->name));
      if (StringICmp (ssp->name, "United States") == 0 ||
          StringICmp (ssp->name, "United States of America") == 0 ||
          StringICmp (ssp->name, "U.S.A.") == 0) {
        ssp->name = MemFree (ssp->name);
        ssp->name = StringSave ("USA");
      }
      if (StringNICmp (ssp->name, "United States:", 14) == 0) {
        str = ssp->name;
        str [0] = ' ';
        str [1] = ' ';
        str [2] = ' ';
        str [3] = ' ';
        str [4] = ' ';
        str [5] = ' ';
        str [6] = ' ';
        str [7] = ' ';
        str [8] = ' ';
        str [9] = ' ';
        str [10] = 'U';
        str [11] = 'S';
        str [12] = 'A';
        TrimSpacesAroundString (ssp->name);
      }
    } else if (ssp->subtype == SUBSRC_clone) {
      CleanVisStringJunk (&(ssp->name));
    }
    if (ssp->subtype == SUBSRC_fwd_primer_seq ||
        ssp->subtype == SUBSRC_rev_primer_seq) {
      if (ssp->name != NULL) {
        /* upper case sequence */
        ptr = ssp->name;
        ch = *ptr;
        while (ch != '\0') {
          if (IS_UPPER (ch)) {
            *ptr = TO_LOWER (ch);
          }
          ptr++;
          ch = *ptr;
        }
        /* remove any spaces in sequence outisde of <modified base> */
        src = ssp->name;
        dst = ssp->name;
        ch = *src;
        while (ch != '\0') {
          if (ch == '<') {
            in_brackets = TRUE;
            *dst = ch;
            dst++;
          } else if (ch == '>') {
            in_brackets = FALSE;
            *dst = ch;
            dst++;
          } else if (ch != ' ') {
            *dst = ch;
            dst++;
          } else if (in_brackets) {
            *dst = ch;
            dst++;
          }
          src++;
          ch = *src;
        }
        *dst = '\0';
        /* upper case modified base <OTHER> */
        ptr = ssp->name;
        tmp = StringStr (ptr, "<other>");
        while (tmp != NULL) {
          ptr = tmp + 7;
          for (i = 1; i < 6; i++) {
            ch = tmp [i];
            tmp [i] = TO_UPPER (ch);
          }
          tmp = StringStr (ptr, "<other>");
        }
      }
    }
    /*
    if (ssp->subtype == SUBSRC_fwd_primer_seq) {
      if (fwd_seq == NULL) {
        fwd_seq = ssp;
      } else {
        fwd_seq->name = CombineSplitQual (fwd_seq->name, ssp->name);
        unlink = TRUE;
      }
    }
    if (ssp->subtype == SUBSRC_rev_primer_seq) {
      if (rev_seq == NULL) {
        rev_seq = ssp;
      } else {
        rev_seq->name = CombineSplitQual (rev_seq->name, ssp->name);
        unlink = TRUE;
      }
    }
    if (ssp->subtype == SUBSRC_fwd_primer_name) {
      if (fwd_name == NULL) {
        fwd_name = ssp;
      } else {
        fwd_name->name = CombineSplitQual (fwd_name->name, ssp->name);
        unlink = TRUE;
      }
    }
    if (ssp->subtype == SUBSRC_rev_primer_name) {
      if (rev_name == NULL) {
        rev_name = ssp;
      } else {
        rev_name->name = CombineSplitQual (rev_name->name, ssp->name);
        unlink = TRUE;
      }
    }
    */
    CleanVisString (&(ssp->attrib));
    if (last != NULL) {
      if (HasNoText (ssp->name) && (! IsNoNameSubSource (ssp))) {
        unlink = TRUE;
      } else if (last->subtype == ssp->subtype &&
                 (IsNoNameSubSource (ssp) ||
                  StringICmp (last->name, ssp->name) == 0 ||
                  (last->subtype == SUBSRC_other &&
                   StringStr (last->name, ssp->name) != NULL))) {
        unlink = TRUE;
      } else if (last->subtype == ssp->subtype &&
                 last->subtype == SUBSRC_other &&
                 IsSubString (last->name, ssp->name)) {
        last->name = MemFree (last->name);
        last->name = ssp->name;
        ssp->name = NULL;
        unlink = TRUE;
      } else if (ssp->subtype == SUBSRC_plastid_name &&
                 location != 0
                 && location == LocationForPlastidText (ssp->name)) {
        unlink = TRUE;
      }
    } else if (HasNoText (ssp->name) && (! IsNoNameSubSource (ssp))) {
      unlink = TRUE;
    } else if (ssp->subtype == SUBSRC_plastid_name &&
               location != 0
               && location == LocationForPlastidText (ssp->name)) {
      unlink = TRUE;
    } else {
      last = ssp;
    }
    if (unlink) {
      *prev = ssp->next;
      ssp->next = NULL;
      SubSourceFree (ssp);
    } else {
      last = ssp;
      prev = &(ssp->next);
    }
    ssp = next;
  }
  /*
  if (fwd_seq != NULL) {
    if (StringChr (fwd_seq->name, ',') != NULL) {
      ptr = fwd_seq->name;
      len = StringLen (ptr);
      if (ptr [0] != '(' || ptr [len - 1] != ')') {
        TrimParenthesesAndCommasAroundString (fwd_seq->name);
        str = MemNew (sizeof (Char) * (len + 4));
        if (str != NULL) {
          StringCpy (str, "(");
          StringCat (str, fwd_seq->name);
          StringCat (str, ")");
          fwd_seq->name = MemFree (fwd_seq->name);
          fwd_seq->name = str;
        }
      }
    }
  }
  if (rev_seq != NULL) {
    if (StringChr (rev_seq->name, ',') != NULL) {
      ptr = rev_seq->name;
      len = StringLen (ptr);
      if (ptr [0] != '(' || ptr [len - 1] != ')') {
        TrimParenthesesAndCommasAroundString (rev_seq->name);
        str = MemNew (sizeof (Char) * (len + 4));
        if (str != NULL) {
          StringCpy (str, "(");
          StringCat (str, rev_seq->name);
          StringCat (str, ")");
          rev_seq->name = MemFree (rev_seq->name);
          rev_seq->name = str;
        }
      }
    }
  }
  if (fwd_name != NULL) {
    if (StringChr (fwd_name->name, ',') != NULL) {
      ptr = fwd_name->name;
      len = StringLen (ptr);
      if (ptr [0] != '(' || ptr [len - 1] != ')') {
        TrimParenthesesAndCommasAroundString (fwd_name->name);
        str = MemNew (sizeof (Char) * (len + 4));
        if (str != NULL) {
          StringCpy (str, "(");
          StringCat (str, fwd_name->name);
          StringCat (str, ")");
          fwd_name->name = MemFree (fwd_name->name);
          fwd_name->name = str;
        }
      }
    }
  }
  if (rev_name != NULL) {
    if (StringChr (rev_name->name, ',') != NULL) {
      ptr = rev_name->name;
      len = StringLen (ptr);
      if (ptr [0] != '(' || ptr [len - 1] != ')') {
        TrimParenthesesAndCommasAroundString (rev_name->name);
        str = MemNew (sizeof (Char) * (len + 4));
        if (str != NULL) {
          StringCpy (str, "(");
          StringCat (str, rev_name->name);
          StringCat (str, ")");
          rev_name->name = MemFree (rev_name->name);
          rev_name->name = str;
        }
      }
    }
  }
  */
}

extern void CleanSubSourcePrimers (SubSourcePtr PNTR sspp)

{
  SubSourcePtr       fwd_seq = NULL, rev_seq = NULL, fwd_name = NULL, rev_name = NULL;
  size_t             len;
  SubSourcePtr       next;
  SubSourcePtr PNTR  prev;
  CharPtr            ptr;
  SubSourcePtr       ssp;
  CharPtr            str;
  Boolean            unlink;

  if (sspp == NULL) return;
  prev = sspp;
  ssp = *sspp;
  while (ssp != NULL) {
    next = ssp->next;
    unlink= FALSE;
    if (ssp->subtype == SUBSRC_fwd_primer_seq) {
      if (fwd_seq == NULL) {
        fwd_seq = ssp;
      } else {
        fwd_seq->name = CombineSplitQual (fwd_seq->name, ssp->name);
        unlink = TRUE;
      }
    }
    if (ssp->subtype == SUBSRC_rev_primer_seq) {
      if (rev_seq == NULL) {
        rev_seq = ssp;
      } else {
        rev_seq->name = CombineSplitQual (rev_seq->name, ssp->name);
        unlink = TRUE;
      }
    }
    if (ssp->subtype == SUBSRC_fwd_primer_name) {
      if (fwd_name == NULL) {
        fwd_name = ssp;
      } else {
        fwd_name->name = CombineSplitQual (fwd_name->name, ssp->name);
        unlink = TRUE;
      }
    }
    if (ssp->subtype == SUBSRC_rev_primer_name) {
      if (rev_name == NULL) {
        rev_name = ssp;
      } else {
        rev_name->name = CombineSplitQual (rev_name->name, ssp->name);
        unlink = TRUE;
      }
    }
    if (unlink) {
      *prev = ssp->next;
      ssp->next = NULL;
      SubSourceFree (ssp);
    } else {
      prev = &(ssp->next);
    }
    ssp = next;
  }
  if (fwd_seq != NULL) {
    if (StringChr (fwd_seq->name, ',') != NULL) {
      ptr = fwd_seq->name;
      len = StringLen (ptr);
      if (ptr [0] != '(' || ptr [len - 1] != ')') {
        TrimParenthesesAndCommasAroundString (fwd_seq->name);
        str = MemNew (sizeof (Char) * (len + 4));
        if (str != NULL) {
          StringCpy (str, "(");
          StringCat (str, fwd_seq->name);
          StringCat (str, ")");
          fwd_seq->name = MemFree (fwd_seq->name);
          fwd_seq->name = str;
        }
      }
    }
  }
  if (rev_seq != NULL) {
    if (StringChr (rev_seq->name, ',') != NULL) {
      ptr = rev_seq->name;
      len = StringLen (ptr);
      if (ptr [0] != '(' || ptr [len - 1] != ')') {
        TrimParenthesesAndCommasAroundString (rev_seq->name);
        str = MemNew (sizeof (Char) * (len + 4));
        if (str != NULL) {
          StringCpy (str, "(");
          StringCat (str, rev_seq->name);
          StringCat (str, ")");
          rev_seq->name = MemFree (rev_seq->name);
          rev_seq->name = str;
        }
      }
    }
  }
  if (fwd_name != NULL) {
    if (StringChr (fwd_name->name, ',') != NULL) {
      ptr = fwd_name->name;
      len = StringLen (ptr);
      if (ptr [0] != '(' || ptr [len - 1] != ')') {
        TrimParenthesesAndCommasAroundString (fwd_name->name);
        str = MemNew (sizeof (Char) * (len + 4));
        if (str != NULL) {
          StringCpy (str, "(");
          StringCat (str, fwd_name->name);
          StringCat (str, ")");
          fwd_name->name = MemFree (fwd_name->name);
          fwd_name->name = str;
        }
      }
    }
  }
  if (rev_name != NULL) {
    if (StringChr (rev_name->name, ',') != NULL) {
      ptr = rev_name->name;
      len = StringLen (ptr);
      if (ptr [0] != '(' || ptr [len - 1] != ')') {
        TrimParenthesesAndCommasAroundString (rev_name->name);
        str = MemNew (sizeof (Char) * (len + 4));
        if (str != NULL) {
          StringCpy (str, "(");
          StringCat (str, rev_name->name);
          StringCat (str, ")");
          rev_name->name = MemFree (rev_name->name);
          rev_name->name = str;
        }
      }
    }
  }
}

/* if string starts with given prefix, return pointer to remaining text */

static CharPtr StringHasPrefix (CharPtr str, CharPtr pref, Boolean novalneeded, Boolean skippref)

{
  Char     ch;
  size_t   len;
  Char     tmp [64];
  CharPtr  val;

  if (StringHasNoText (str) || StringHasNoText (pref)) return NULL;
  len = StringLen (pref);
  StringNCpy_0 (tmp, pref, sizeof (tmp));
  if (StringNICmp (str, tmp, len) != 0) {
    /* try after replacing dash with underscore */
    val = tmp;
    ch = *val;
    while (ch != '\0') {
      if (ch == '-') {
        *val = '_';
      }
      val++;
      ch = *val;
    }
    if (StringNICmp (str, tmp, len) != 0) return NULL;
  }
  if (skippref) {
    val = str + len;
  } else {
    val = str;
  }
  if (StringHasNoText (val)) {
    if (novalneeded) return " ";
    return NULL;
  }
  ch = *(str + len);
  if (ch != '=' && ch != ' ' && ch != ':' && ch != '\0') return NULL;
  ch = *val;
  while (ch == '=' || ch == ' ' || ch == ':') {
    val++;
    ch = *val;
  }
  if (StringHasNoText (val)) return NULL;
  return val;
}


Nlm_QualNameAssoc current_orgmod_subtype_alist[] = {
  {" ",                   0},
  {"Acronym",            ORGMOD_acronym},
  {"Anamorph",           ORGMOD_anamorph},
  {"Authority",          ORGMOD_authority},
  {"Bio-material",       ORGMOD_bio_material},
  {"Biotype",            ORGMOD_biotype},
  {"Biovar",             ORGMOD_biovar},
  {"Breed",              ORGMOD_breed},
  {"Chemovar",           ORGMOD_chemovar},
  {"Common",             ORGMOD_common},
  {"Cultivar",           ORGMOD_cultivar},
  {"Culture-collection", ORGMOD_culture_collection},
  {"Ecotype",            ORGMOD_ecotype},
  {"Forma",              ORGMOD_forma},
  {"Forma-specialis",    ORGMOD_forma_specialis},
  {"Group",              ORGMOD_group},
  {"Host",               ORGMOD_nat_host},
  {"Isolate",            ORGMOD_isolate},
  {"Metagenome-source",  ORGMOD_metagenome_source},
  {"Pathovar",           ORGMOD_pathovar},
  {"Serogroup",          ORGMOD_serogroup},
  {"Serotype",           ORGMOD_serotype},
  {"Serovar",            ORGMOD_serovar},
  {"Specimen-voucher",   ORGMOD_specimen_voucher},
  {"Strain",             ORGMOD_strain},
  {"Subgroup",           ORGMOD_subgroup},
  {"Sub-species",        ORGMOD_sub_species},
  {"Substrain",          ORGMOD_substrain},
  {"Subtype",            ORGMOD_subtype},
  {"Synonym",            ORGMOD_synonym},
  {"Teleomorph",         ORGMOD_teleomorph},
  {"Type",               ORGMOD_type},
  {"Variety",            ORGMOD_variety},
  { NULL, 0 } };

Nlm_QualNameAssoc discouraged_orgmod_subtype_alist[] = {
  {"Old Lineage",      ORGMOD_old_lineage},
  {"Old Name",         ORGMOD_old_name},
  { NULL, 0 } };

Nlm_QualNameAssoc discontinued_orgmod_subtype_alist[] = {
  {"Dosage",           ORGMOD_dosage},
  { NULL, 0 } };


Nlm_NameNameAssoc orgmod_aliases[] = {
  {"Sub-species",   "subspecies", ORGMOD_sub_species},
  {"Host", "nat-host",   ORGMOD_nat_host},
  {"Host", "specific-host",   ORGMOD_nat_host},
  {"Substrain",   "Sub_strain", ORGMOD_substrain},
  { NULL, NULL, 0 } };

extern CharPtr GetOrgModQualName (Uint1 subtype)
{
  Int4 i;
  
  if (subtype == ORGMOD_other) {
    return "Note";
  }
  for (i = 0; current_orgmod_subtype_alist[i].name != NULL; i++) {
    if (current_orgmod_subtype_alist[i].value == subtype) {
      return current_orgmod_subtype_alist[i].name;
    }
  }
  for (i = 0; discouraged_orgmod_subtype_alist[i].name != NULL; i++) {
    if (discouraged_orgmod_subtype_alist[i].value == subtype) {
      return discouraged_orgmod_subtype_alist[i].name;
    }
  }

  for (i = 0; discontinued_orgmod_subtype_alist[i].name != NULL; i++) {
    if (discontinued_orgmod_subtype_alist[i].value == subtype) {
      return discontinued_orgmod_subtype_alist[i].name;
    }
  }

  return NULL;
}


extern void BioSourceHasOldOrgModQualifiers (BioSourcePtr biop, BoolPtr has_discouraged, BoolPtr has_discontinued)
{
  OrgModPtr mod;
  Boolean   discouraged = FALSE, discontinued = FALSE;
  Int4      i;

  if (biop != NULL && biop->org != NULL && biop->org->orgname != NULL) {
    mod = biop->org->orgname->mod;
    while (mod != NULL && (!discouraged || !discontinued)) {
      for (i = 0; discouraged_orgmod_subtype_alist[i].name != NULL && !discouraged; i++) {
        if (mod->subtype == discouraged_orgmod_subtype_alist[i].value) {
          discouraged = TRUE;
        }
      }
      for (i = 0; discontinued_orgmod_subtype_alist[i].name != NULL && !discontinued; i++) {
        if (mod->subtype == discontinued_orgmod_subtype_alist[i].value) {
          discontinued = TRUE;
        }
      }
      mod = mod->next;
    }
  }

  if (has_discouraged != NULL) {
    *has_discouraged = discouraged;
  }
  if (has_discontinued != NULL) {
    *has_discontinued = discontinued;
  }
}


static void StringHasOrgModPrefix (CharPtr str, CharPtr PNTR pval, Uint1Ptr p_subtypeval, Boolean skippref)
{
  Int2          i;
  CharPtr       val = NULL;
  Uint1         subtype_val = 0;
  
  for (i = 0; current_orgmod_subtype_alist[i].name != NULL && subtype_val == 0; i++) {
    if (current_orgmod_subtype_alist[i].value == ORGMOD_nat_host) continue;
    val = StringHasPrefix (str, current_orgmod_subtype_alist [i].name, FALSE, skippref);
    if (val != NULL) {
      subtype_val = current_orgmod_subtype_alist[i].value;
    }
  }  
  if (subtype_val == 0) {
    for (i = 0; orgmod_aliases[i].name != NULL && subtype_val == 0; i++) {
      if (orgmod_aliases[i].value == ORGMOD_nat_host) continue;
      val = StringHasPrefix (str, orgmod_aliases [i].alias, FALSE, skippref);
      if (val != NULL) {
        subtype_val = orgmod_aliases[i].value;
      }
    }
  }
  if (pval != NULL) {
    *pval = val;
  }
  if (p_subtypeval != NULL) {
    *p_subtypeval = subtype_val;
  }
}

static void OrpModToOrgMod (ValNodePtr PNTR vnpp, OrgModPtr PNTR ompp)

{
  Char        ch;
  ValNodePtr  next;
  Int2        numcommas;
  Int2        numspaces;
  OrgModPtr   omp;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     val;
  ValNodePtr  vnp;
  Uint1       subtype;

  if (vnpp == NULL || ompp == NULL) return;
  vnp = *vnpp;
  while (vnp != NULL) {
    next = vnp->next;
    str = (CharPtr) vnp->data.ptrvalue;
    val = NULL;
    subtype = 0;
    StringHasOrgModPrefix (str, &val, &subtype, TRUE);
    if (val != NULL) {
      numspaces = 0;
      numcommas = 0;
      ptr = str;
      ch = *ptr;
      while (ch != '\0') {
        if (ch == ' ') {
          numspaces++;
        } else if (ch == ',') {
          numcommas++;
        }
        ptr++;
        ch = *ptr;
      }
      if (numspaces > 4 || numcommas > 0) {
        val = NULL;
      }
    }
    if (val != NULL) {
      omp = OrgModNew ();
      if (omp != NULL) {
        omp->subtype = (Uint1) subtype;
        omp->subname = StringSave (val);
        omp->next = *ompp;
        *ompp = omp;
      }
      *vnpp = vnp->next;
      vnp->next = NULL;
      ValNodeFreeData (vnp);
    } else {
      vnpp = &(vnp->next);
    }
    vnp = next;
  }
}

Nlm_QualNameAssoc current_subsource_subtype_alist[] = {
  {" ",                      0},
  {"Cell-line",             SUBSRC_cell_line},
  {"Cell-type",             SUBSRC_cell_type},
  {"Chromosome",            SUBSRC_chromosome},
  {"Clone",                 SUBSRC_clone},
  {"Clone-lib",             SUBSRC_clone_lib},
  {"Collected-by",          SUBSRC_collected_by},
  {"Collection-date",       SUBSRC_collection_date},
  {"Country",               SUBSRC_country},
  {"Dev-stage",             SUBSRC_dev_stage},
  {"Endogenous-virus-name", SUBSRC_endogenous_virus_name},
  {"Environmental-sample",  SUBSRC_environmental_sample},
  {"Frequency",             SUBSRC_frequency},
  {"Genotype",              SUBSRC_genotype},
  {"Germline",              SUBSRC_germline},
  {"Haplogroup",            SUBSRC_haplogroup},
  {"Haplotype",             SUBSRC_haplotype},
  {"Identified-by",         SUBSRC_identified_by},
  {"Isolation-source",      SUBSRC_isolation_source},
  {"Lab-host",              SUBSRC_lab_host},
  {"Lat-Lon",               SUBSRC_lat_lon},
  {"Linkage-group",         SUBSRC_linkage_group},
  {"Map",                   SUBSRC_map},
  {"Mating-type",           SUBSRC_mating_type},
  {"Metagenomic",           SUBSRC_metagenomic},
  {"Plasmid-name",          SUBSRC_plasmid_name},
  {"Pop-variant",           SUBSRC_pop_variant},
  {"Rearranged",            SUBSRC_rearranged},
  {"Segment",               SUBSRC_segment},
  {"Sex",                   SUBSRC_sex},
  {"Subclone",              SUBSRC_subclone},
  {"Tissue-lib",            SUBSRC_tissue_lib},
  {"Tissue-type",           SUBSRC_tissue_type},
  {"Transgenic",            SUBSRC_transgenic},
  { NULL, 0 } };

Nlm_QualNameAssoc discouraged_subsource_subtype_alist[] = {
  {"Plastid-name",          SUBSRC_plastid_name},
  { NULL, 0 } };

Nlm_QualNameAssoc discontinued_subsource_subtype_alist[] = {
  {"Ins-seq-name",          SUBSRC_insertion_seq_name},
  {"Transposon-name",       SUBSRC_transposon_name},
  {"Fwd-PCR-primer-name",   SUBSRC_fwd_primer_name},
  {"Fwd-PCR-primer-seq",    SUBSRC_fwd_primer_seq},
  {"Rev-PCR-primer-name",   SUBSRC_rev_primer_name},
  {"Rev-PCR-primer-seq",    SUBSRC_rev_primer_seq},
  { NULL, 0 } };

Nlm_NameNameAssoc subsource_aliases[] = {
  {"Fwd-PCR-primer-name", "fwd-primer-name",    SUBSRC_fwd_primer_name},
  {"Fwd-PCR-primer-seq",  "fwd-primer-seq",     SUBSRC_fwd_primer_seq},
  {"Rev-PCR-primer-name", "rev-primer-name",    SUBSRC_rev_primer_name},
  {"Rev-PCR-primer-seq",  "rev-primer-seq",     SUBSRC_rev_primer_seq},
  {"Subclone",            "sub-clone",          SUBSRC_subclone},
  {"Lat-Lon",             "Lat-long",           SUBSRC_lat_lon},
  {"Lat-Lon",             "Latitude-Longitude", SUBSRC_lat_lon },
  { NULL, NULL, 0 } };

extern CharPtr GetSubsourceQualName (Uint1 subtype)
{
  Int4 i;
  
  if (subtype == SUBSRC_other) {
    return "Note";
  }
  for (i = 0; current_subsource_subtype_alist[i].name != NULL; i++) {
    if (current_subsource_subtype_alist[i].value == subtype) {
      return current_subsource_subtype_alist[i].name;
    }
  }

  for (i = 0; discouraged_subsource_subtype_alist[i].name != NULL; i++) {
    if (discouraged_subsource_subtype_alist[i].value == subtype) {
      return discouraged_subsource_subtype_alist[i].name;
    }
  }

  for (i = 0; discontinued_subsource_subtype_alist[i].name != NULL; i++) {
    if (discontinued_subsource_subtype_alist[i].value == subtype) {
      return discontinued_subsource_subtype_alist[i].name;
    }
  }

  return NULL;
}


extern void BioSourceHasOldSubSourceQualifiers (BioSourcePtr biop, BoolPtr has_discouraged, BoolPtr has_discontinued)
{
  SubSourcePtr ssp;
  Boolean   discouraged = FALSE, discontinued = FALSE;
  Int4      i;

  if (biop != NULL) {
    ssp = biop->subtype;
    while (ssp != NULL && (!discouraged || !discontinued)) {
      for (i = 0; discouraged_subsource_subtype_alist[i].name != NULL && !discouraged; i++) {
        if (ssp->subtype == discouraged_subsource_subtype_alist[i].value) {
          discouraged = TRUE;
        }
      }
      for (i = 0; discontinued_subsource_subtype_alist[i].name != NULL && !discontinued; i++) {
        if (ssp->subtype == discontinued_subsource_subtype_alist[i].value) {
          discontinued = TRUE;
        }
      }
      ssp = ssp->next;
    }
  }

  if (has_discouraged != NULL) {
    *has_discouraged = discouraged;
  }
  if (has_discontinued != NULL) {
    *has_discontinued = discontinued;
  }
}


static void StringHasSubSourcePrefix (CharPtr str, CharPtr PNTR pval, Uint1Ptr p_subtypeval, Boolean skippref)
{
  Int2          i;
  CharPtr       val = NULL;
  Uint1         subtype_val = 0;
  
  for (i = 0; current_subsource_subtype_alist[i].name != NULL && subtype_val == 0; i++) {
    val = StringHasPrefix (str, current_subsource_subtype_alist [i].name,
                           (Boolean) (current_subsource_subtype_alist[i].value == SUBSRC_germline ||
                                      current_subsource_subtype_alist[i].value == SUBSRC_rearranged ||
                                      current_subsource_subtype_alist[i].value == SUBSRC_transgenic ||
                                      current_subsource_subtype_alist[i].value == SUBSRC_environmental_sample ||
                                      current_subsource_subtype_alist[i].value == SUBSRC_metagenomic),
                           skippref);
    if (val != NULL) {
      subtype_val = current_subsource_subtype_alist[i].value;
    }
  }  
  if (subtype_val == 0) {
    for (i = 0; subsource_aliases[i].name != NULL && subtype_val == 0; i++) {
      val = StringHasPrefix (str, subsource_aliases [i].alias,
                             (Boolean) (subsource_aliases[i].value == SUBSRC_germline ||
                                        subsource_aliases[i].value == SUBSRC_rearranged ||
                                        subsource_aliases[i].value == SUBSRC_transgenic ||
                                        subsource_aliases[i].value == SUBSRC_environmental_sample ||
                                        subsource_aliases[i].value == SUBSRC_metagenomic),
                             skippref);
      if (val != NULL) {
        subtype_val = subsource_aliases[i].value;
      }
    }
  }
  if (pval != NULL) {
    *pval = val;
  }
  if (p_subtypeval != NULL) {
    *p_subtypeval = subtype_val;
  }
}

static void OrpModToSubSource (ValNodePtr PNTR vnpp, SubSourcePtr PNTR sspp)

{
  Char          ch;
  ValNodePtr    next;
  Int2          numcommas;
  Int2          numspaces;
  CharPtr       ptr;
  SubSourcePtr  ssp;
  CharPtr       str;
  CharPtr       val;
  ValNodePtr    vnp;
  Uint1         subtype_val = 0;

  if (vnpp == NULL || sspp == NULL) return;
  vnp = *vnpp;
  while (vnp != NULL) {
    next = vnp->next;
    str = (CharPtr) vnp->data.ptrvalue;
    val = NULL;
    subtype_val = 0;
    StringHasSubSourcePrefix (str, &val, &subtype_val, TRUE);

    if (val != NULL) {
      numspaces = 0;
      numcommas = 0;
      ptr = str;
      ch = *ptr;
      while (ch != '\0') {
        if (ch == ' ') {
          numspaces++;
        } else if (ch == ',') {
          numcommas++;
        }
        ptr++;
        ch = *ptr;
      }
      if (numspaces > 4 || numcommas > 0) {
        val = NULL;
      }
    }
    if (val != NULL) {
      ssp = SubSourceNew ();
      if (ssp != NULL) {
        ssp->subtype = subtype_val;
        ssp->name = StringSave (val);
        ssp->next = *sspp;
        *sspp = ssp;
      }
      *vnpp = vnp->next;
      vnp->next = NULL;
      ValNodeFreeData (vnp);
    } else {
      vnpp = &(vnp->next);
    }
    vnp = next;
  }
}

static void GbqualToOrpMod (GBQualPtr PNTR prevgbq, ValNodePtr PNTR vnpp)

{
  GBQualPtr  gbq;
  size_t     len;
  GBQualPtr  next;
  CharPtr    str;
  Boolean    unlink;
  CharPtr    val;
  Uint1      subtype_val;

  if (prevgbq == NULL) return;
  gbq = *prevgbq;
  while (gbq != NULL) {
    next = gbq->next;
    unlink = FALSE;
    str = gbq->qual;
    if (str != NULL) {
      val = NULL;
      subtype_val = 0;
      StringHasOrgModPrefix (str, &val, &subtype_val, FALSE);
      if (val == NULL) {
        subtype_val = 0;
        StringHasSubSourcePrefix (str, &val, &subtype_val, FALSE);

      }
      if (val != NULL) {
        len = StringLen (gbq->val);
        str = MemNew (sizeof (Char) * (len + 64));
        if (str != NULL) {
          StringCpy (str, val);
          StringCat (str, "=");
          StringCat (str, gbq->val);
          ValNodeAddStr (vnpp, 0, str);
          unlink = TRUE; 
        }
      }
    }
    if (unlink) {
      *prevgbq = gbq->next;
      gbq->next = NULL;
      GBQualFree (gbq);
    } else {
      prevgbq = (GBQualPtr PNTR) &(gbq->next);
    }
    gbq = next;
  }
}

#define IS_WHITESP(c) (((c) == ' ') || ((c) == '\n') || ((c) == '\r') || ((c) == '\t'))

static Boolean IsStringSingleToken (CharPtr str)

{
  Char  ch;

  if (StringHasNoText (str)) return FALSE;

  ch = *str;
  while (ch != '\0') {
    if (IS_WHITESP (ch)) return FALSE;
    str++;
    ch = *str;
  }

  return TRUE;
}

static CharPtr FindAnOrgMod (OrgNamePtr onp, Uint1 subtype)

{
  OrgModPtr  omp;

  if (onp == NULL || subtype == 0) return NULL;

  for (omp = onp->mod; omp != NULL; omp = omp->next) {
    if (omp->subtype != subtype) continue;
    if (StringHasNoText (omp->subname)) continue;
    return omp->subname;
  }

  return NULL;
}

static CharPtr FindASubSource (BioSourcePtr biop, Uint1 subtype)

{
  SubSourcePtr  ssp;

  if (biop == NULL || subtype == 0) return NULL;

  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype != subtype) continue;
    if (StringHasNoText (ssp->name)) continue;
    return ssp->name;
  }

  return NULL;
}

static CharPtr FindNextSingleTilde (CharPtr str)

{
  Char  ch;

  if (StringHasNoText (str)) return NULL;

  ch = *str;
  while (ch != '\0') {
    if (ch == ' ') {
      if (str [1] == '~') {
        str++;
        ch = *str;
        while (ch == '~') {
          str++;
          ch = *str;
        }
      } else {
        str++;
        ch = *str;
      }
    } else if (ch == '~') {
      if (str [1] != '~') return str;
      str++;
      ch = *str;
      while (ch == '~') {
        str++;
        ch = *str;
      }
    } else {
      str++;
      ch = *str;
    }
  }

  return NULL;
}

static ValNodePtr SplitAtSingleTilde (CharPtr strs)

{
  ValNodePtr  head = NULL;
  CharPtr     ptr, str, tmp;

  if (StringHasNoText (strs)) return NULL;

  tmp = StringSave (strs);
  str = tmp;

  while (StringDoesHaveText (str)) {
    ptr = FindNextSingleTilde (str);
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }
    TrimSpacesAroundString (str);
    ValNodeCopyStr (&head, 0, str);
    str = ptr;
  }

  MemFree (tmp);
  return head;
}

static CharPtr MergeTildeStrings (ValNodePtr head)

{
  size_t      len = 0;
  CharPtr     prefix = "", ptr, str;
  ValNodePtr  vnp;

  if (head == NULL) return NULL;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    len += StringLen (str) + 1;
  }
  if (len < 1) return NULL;

  ptr = MemNew (sizeof (Char) * (len + 2));
  if (ptr == NULL) return NULL;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    StringCat (ptr, prefix);
    StringCat (ptr, str);
    prefix = "~";
  }

  return ptr;
}

static void CleanupOrgModOther (BioSourcePtr biop, OrgNamePtr onp)

{
  ValNodePtr      head, vnp;
  OrgModPtr       next;
  OrgModPtr       omp;
  OrgModPtr PNTR  prev;
  CharPtr         str;
  Uint1           subtype_val;
  CharPtr         tmp;
  Boolean         unlink;
  CharPtr         val;

  if (biop == NULL || onp == NULL) return;

  prev = &(onp->mod);
  omp = onp->mod;
  while (omp != NULL) {
    next = omp->next;
    unlink= FALSE;
    if (omp->subtype == ORGMOD_other) {
      str = omp->subname;
      head = SplitAtSingleTilde (str);
      for (vnp = head; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        val = NULL;
        subtype_val = 0;
        StringHasOrgModPrefix (str, &val, &subtype_val, TRUE);
        if (val != NULL) {
          tmp = FindAnOrgMod (onp, subtype_val);
          if (tmp != NULL && StringICmp (tmp, val) == 0) {
            vnp->data.ptrvalue = NULL;
          }
        } else {
          subtype_val = 0;
          StringHasSubSourcePrefix (str, &val, &subtype_val, TRUE);
          if (val != NULL) {
            tmp = FindASubSource (biop, subtype_val);
            if (tmp != NULL && StringICmp (tmp, val) == 0) {
              vnp->data.ptrvalue = NULL;
            }
          }
        }
      }
      str = MergeTildeStrings (head);
      ValNodeFreeData (head);
      omp->subname = MemFree (omp->subname);
      omp->subname = str;
      if (StringHasNoText (str)) {
        unlink = TRUE;
      }
    }
    if (unlink) {
      *prev = omp->next;
      omp->next = NULL;
      OrgModFree (omp);
    } else {
      prev = &(omp->next);
    }
    omp = next;
  }
}

static void CleanupSubSourceOther (BioSourcePtr biop, OrgNamePtr onp)

{
  ValNodePtr         head, vnp;
  SubSourcePtr       next;
  SubSourcePtr PNTR  prev;
  SubSourcePtr       ssp;
  CharPtr            str;
  Uint1              subtype_val;
  CharPtr            tmp;
  Boolean            unlink;
  CharPtr            val;

  if (biop == NULL /* || onp == NULL */ ) return;

  prev = &(biop->subtype);
  ssp = biop->subtype;
  while (ssp != NULL) {
    next = ssp->next;
    unlink = FALSE;
    if (ssp->subtype == SUBSRC_other) {
      str = ssp->name;
      head = SplitAtSingleTilde (str);
      for (vnp = head; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        val = NULL;
        subtype_val = 0;
        StringHasOrgModPrefix (str, &val, &subtype_val, TRUE);
        if (val != NULL) {
          tmp = FindAnOrgMod (onp, subtype_val);
          if (tmp != NULL && StringICmp (tmp, val) == 0) {
            vnp->data.ptrvalue = NULL;
          }
        } else {
          subtype_val = 0;
          StringHasSubSourcePrefix (str, &val, &subtype_val, TRUE);
          if (val != NULL) {
            tmp = FindASubSource (biop, subtype_val);
            if (tmp != NULL && StringICmp (tmp, val) == 0) {
              vnp->data.ptrvalue = NULL;
            }
          }
        }
      }
      str = MergeTildeStrings (head);
      ValNodeFreeData (head);
      ssp->name = MemFree (ssp->name);
      ssp->name = str;
      if (StringHasNoText (str)) {
        unlink = TRUE;
      }
    }
    if (unlink) {
      *prev = ssp->next;
      ssp->next = NULL;
      SubSourceFree (ssp);
    } else {
      prev = &(ssp->next);
    }
    ssp = next;
  }
}

static int LIBCALLBACK SortDbxref (VoidPtr ptr1, VoidPtr ptr2)

{
  int          compare;
  DbtagPtr     dbt1;
  DbtagPtr     dbt2;
  ObjectIdPtr  oip1;
  ObjectIdPtr  oip2;
  CharPtr      str1;
  CharPtr      str2;
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  dbt1 = (DbtagPtr) vnp1->data.ptrvalue;
  dbt2 = (DbtagPtr) vnp2->data.ptrvalue;
  if (dbt1 == NULL || dbt2 == NULL) return 0;
  str1 = (CharPtr) dbt1->db;
  str2 = (CharPtr) dbt2->db;
  if (str1 == NULL || str2 == NULL) return 0;
  compare = StringICmp (str1, str2);
  if (compare != 0) return compare;
  oip1 = dbt1->tag;
  oip2 = dbt2->tag;
  if (oip1 == NULL || oip2 == NULL) return 0;
  str1 = oip1->str;
  str2 = oip2->str;
  if (str1 != NULL && str2 != NULL) {
    return StringICmp (str1, str2);
  } else if (str1 == NULL && str2 == NULL) {
    if (oip1->id > oip2->id) {
      return 1;
    } else if (oip1->id < oip2->id) {
      return -1;
    }
  } else if (str1 != NULL) {
    return 1;
  } else if (str2 != NULL) {
    return -1;
  }
  return 0;
}

static void FixNumericDbxrefs (ValNodePtr vnp)

{
  Char         ch;
  DbtagPtr     dbt;
  Boolean      isNum;
  Boolean      leadingzero;
  Boolean      notallzero;
  ObjectIdPtr  oip;
  CharPtr      ptr;
  long         val;

  while (vnp != NULL) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      oip = dbt->tag;
      if (oip != NULL) {
        ptr = oip->str;
        if (ptr != NULL) {
          leadingzero = FALSE;
          notallzero = FALSE;
          isNum = TRUE;
          ch = *ptr;
          if (ch == '0') {
            leadingzero = TRUE;
          }
          while (ch != '\0') {
            if ((! IS_DIGIT (ch)) && (! IS_WHITESP (ch))) {
              isNum = FALSE;
            } else if ('1'<= ch && ch <='9') {
              notallzero = TRUE;
            }
            ptr++;
            ch = *ptr;
          }
          if (isNum) {
            if (leadingzero && notallzero) {
              /* suppress conversion */
            } else if (sscanf (oip->str, "%ld", &val) == 1) {
              oip->id = (Int4) val;
              oip->str = MemFree (oip->str);
            }
          }
        }
      }
    }
    vnp = vnp->next;
  }
}

static void FixOldDbxrefs (ValNodePtr vnp)

{
  Boolean      all_digits;
  Char         ch;
  DbtagPtr     dbt;
  ObjectIdPtr  oip;
  CharPtr      ptr;
  CharPtr      str;
  CharPtr      tmp;
  ValNodePtr   vp2;

  while (vnp != NULL) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {

      TrimSpacesAroundString (dbt->db);
      oip = dbt->tag;
      if (oip != NULL && oip->str != NULL) {
        /*
        TrimSpacesAroundString (oip->str);
        */
        TrimSpacesSemicolonsAndCommasAndParens (oip->str);
      }

      if (StringICmp (dbt->db, "SWISS-PROT") == 0 &&
          StringCmp (dbt->db, "Swiss-Prot") != 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("Swiss-Prot");
      } else if (StringICmp (dbt->db, "SPTREMBL") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("TrEMBL");
      } else if (StringICmp (dbt->db, "SUBTILIS") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("SubtiList");
      } else if (StringICmp (dbt->db, "MGD") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("MGI");
      } else if (StringCmp (dbt->db, "cdd") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("CDD");
      } else if (StringCmp (dbt->db, "FlyBase") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("FLYBASE");
      } else if (StringCmp (dbt->db, "GENEDB") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("GeneDB");
      } else if (StringCmp (dbt->db, "GreengenesID") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("Greengenes");
      } else if (StringCmp (dbt->db, "HMPID") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("HMP");
      }
      if (StringICmp (dbt->db, "HPRD") == 0) {
        oip = dbt->tag;
        if (oip != NULL && StringDoesHaveText (oip->str)) {
          str = oip->str;
          if (StringNICmp (str, "HPRD_", 5) == 0) {
            str [0] = ' ';
            str [1] = ' ';
            str [2] = ' ';
            str [3] = ' ';
            str [4] = ' ';
            TrimSpacesAroundString (str);
          }
        }
      } else if (StringICmp (dbt->db, "MGI") == 0) {
        oip = dbt->tag;
        if (oip != NULL && StringDoesHaveText (oip->str)) {
          str = oip->str;
          if (StringNICmp (str, "MGI:", 4) == 0 || StringNICmp (str, "MGD:", 4) == 0) {
            str [0] = ' ';
            str [1] = ' ';
            str [2] = ' ';
            str [3] = ' ';
            TrimSpacesAroundString (str);
          } else if (StringNICmp (str, "J:", 2) == 0) {
            ptr = str + 2;
            ch = *ptr;
            all_digits = TRUE;
            while (ch != '\0') {
              if (! IS_DIGIT (ch)) {
                all_digits = FALSE;
              }
              ptr++;
              ch = *ptr;
            }
            if (all_digits) {
              oip->str = MemFree (oip->str);
              oip->str = StringSave ("");
            }
          }
        }
      }
      if (StringICmp (dbt->db, "Swiss-Prot") == 0 ||
          StringICmp (dbt->db, "SWISSPROT") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("UniProt/Swiss-Prot");
      } else if (StringICmp (dbt->db, "TrEMBL") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("UniProt/TrEMBL");
      } else if (StringICmp (dbt->db, "LocusID") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("GeneID");
      } else if (StringICmp (dbt->db, "MaizeDB") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("MaizeGDB");
      }
      if (StringICmp (dbt->db, "UniProt/Swiss-Prot") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("UniProtKB/Swiss-Prot");
      } else if (StringICmp (dbt->db, "UniProt/TrEMBL") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("UniProtKB/TrEMBL");
      } else if (StringICmp (dbt->db, "Genew") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("HGNC");
      } else if (StringICmp (dbt->db, "IFO") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("NBRC");
      } else if (StringICmp (dbt->db, "BHB") == 0 ||
          StringICmp (dbt->db, "BioHealthBase") == 0) {
        dbt->db = MemFree (dbt->db);
        dbt->db = StringSave ("IRD");
      }

      /* expand db_xrefs with colons inside tags */

      oip = dbt->tag;
      if (oip != NULL && oip->str != NULL) {
        ptr = StringChr (oip->str, ':');
        if (ptr != NULL) {
          if (StringHasNoText (ptr + 1)) {
            *ptr = '\0';
          } else {
            tmp = dbt->db;
            dbt = DbtagNew ();
            if (dbt != NULL) {
              oip = ObjectIdNew ();
              if (oip != NULL) {
                vp2 = ValNodeNew (NULL);
                if (vp2 != NULL) {
                  *ptr = '\0';
                  ptr++;
                  TrimSpacesAroundString (ptr);
                  dbt->db = StringSave (tmp);
                  oip->str = StringSave (ptr);
                  dbt->tag = oip;
                  vp2->data.ptrvalue = (Pointer) dbt;
                  vp2->next = vnp->next;
                  vnp->next = vp2;
                }
              }
            }
          }
        }
      }
    }

    vnp = vnp->next;
  }
}

static void CleanupDuplicateDbxrefs (ValNodePtr PNTR prevvnp)

{
  DbtagPtr     dbt;
  DbtagPtr     last = NULL;
  ValNodePtr   nextvnp;
  ObjectIdPtr  oip1;
  ObjectIdPtr  oip2;
  CharPtr      str1;
  CharPtr      str2;
  Boolean      unlink;
  ValNodePtr   vnp;

  if (prevvnp == NULL) return;
  vnp = *prevvnp;
  while (vnp != NULL) {
    nextvnp = vnp->next;
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      unlink = FALSE;
      if (last != NULL) {
        str1 = (CharPtr) dbt->db;
        str2 = (CharPtr) last->db;
        if (str1 != NULL && str2 != NULL && StringICmp (str1, str2) == 0) {
          oip1 = dbt->tag;
          oip2 = last->tag;
          if (oip1 != NULL && oip2 != NULL) {
            str1 = oip1->str;
            str2 = oip2->str;
            if (str1 != NULL && str2 != NULL) {
              if (StringICmp (str1, str2) == 0) {
                unlink = TRUE;
              }
            } else if (str1 == NULL && str2 == NULL) {
              if (oip1->id == oip2->id) {
                unlink = TRUE;
              }
            }
          }
        }
      } else {
        last = dbt;
      }
      if (unlink) {
        *prevvnp = vnp->next;
        vnp->next = NULL;
        DbtagFree (dbt);
        ValNodeFree (vnp);
      } else {
        last = dbt;
        prevvnp = (ValNodePtr PNTR) &(vnp->next);
      }
    }
    vnp = nextvnp;
  }
}

static void CleanupObsoleteDbxrefs (ValNodePtr PNTR prevvnp)

{
  DbtagPtr     dbt;
  ValNodePtr   nextvnp;
  ObjectIdPtr  oip;
  CharPtr      str;
  Boolean      unlink;
  ValNodePtr   vnp;

  if (prevvnp == NULL) return;
  vnp = *prevvnp;
  while (vnp != NULL) {
    nextvnp = vnp->next;
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      unlink = FALSE;
      str = (CharPtr) dbt->db;
      if (StringHasNoText (str) ||
          StringICmp (str, "PID") == 0 ||
          StringICmp (str, "PIDg") == 0 ||
          /*
          StringICmp (str, "PIDe") == 0 ||
          StringICmp (str, "PIDd") == 0 ||
          */
          /*
          StringICmp (str, "GI") == 0 ||
          */
          StringICmp (str, "NID") == 0) {
        unlink = TRUE;
      }
      oip = dbt->tag;
      if (oip == NULL) {
        unlink = TRUE;
      } else if (oip->str != NULL) {
        if (StringHasNoText (oip->str)) {
          unlink = TRUE;
        }
      } else if (oip->id == 0) {
        unlink = TRUE;
      }
      if (unlink) {
        *prevvnp = vnp->next;
        vnp->next = NULL;
        DbtagFree (dbt);
        ValNodeFree (vnp);
      } else {
        prevvnp = (ValNodePtr PNTR) &(vnp->next);
      }
    }
    vnp = nextvnp;
  }
}

static int LIBCALLBACK SortCits (VoidPtr ptr1, VoidPtr ptr2)

{
  int         compare;
  Char        label1 [128], label2 [128];
  ValNodePtr  ppr1, ppr2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  ppr1 = *((ValNodePtr PNTR) ptr1);
  ppr2 = *((ValNodePtr PNTR) ptr2);
  if (ppr1 == NULL || ppr2 == NULL) return 0;
  PubLabel (ppr1, label1, 127, OM_LABEL_CONTENT);
  PubLabel (ppr2, label2, 127, OM_LABEL_CONTENT);
  compare = StringICmp (label1, label2);
  return compare;
}

static Boolean CitGenTitlesMatch (ValNodePtr pub1, ValNodePtr pub2)

{
  CitGenPtr  cgp1, cgp2;

  if (pub1->choice == PUB_Gen) {
    cgp1 = (CitGenPtr) pub1->data.ptrvalue;
    if (cgp1->serial_number != -1 && pub1->next != NULL) {
      pub1 = pub1->next;
    }
  }
  if (pub2->choice == PUB_Gen) {
    cgp2 = (CitGenPtr) pub2->data.ptrvalue;
    if (cgp2->serial_number != -1 && pub2->next != NULL) {
      pub2 = pub2->next;
    }
  }

  if (pub1->choice != PUB_Gen || pub2->choice != PUB_Gen) return TRUE;
  cgp1 = (CitGenPtr) pub1->data.ptrvalue;
  cgp2 = (CitGenPtr) pub2->data.ptrvalue;
  if (cgp1->title == NULL || cgp2->title == NULL) return TRUE;
  if (StringCmp (cgp1->title, cgp2->title) != 0) return FALSE;
  return TRUE;
}

static void CleanupDuplicateCits (ValNodePtr PNTR prevvnp)

{
  Char        label1 [128], label2 [128];
  ValNodePtr  last = NULL;
  ValNodePtr  nextvnp;
  Boolean     unlink;
  ValNodePtr  vnp;

  if (prevvnp == NULL) return;
  vnp = *prevvnp;
  while (vnp != NULL) {
    nextvnp = vnp->next;
    unlink = FALSE;
    if (last != NULL) {
      PubLabelUnique (last, label1, 127, OM_LABEL_CONTENT, TRUE);
      PubLabelUnique (vnp, label2, 127, OM_LABEL_CONTENT, TRUE);
      if (StringCmp (label1, label2) == 0 && CitGenTitlesMatch (last, vnp)) {
        unlink = TRUE;
      }
    } else {
      last = vnp;
    }
    if (unlink) {
      *prevvnp = vnp->next;
      vnp->next = NULL;
      PubFree (vnp);
    } else {
      last = vnp;
      prevvnp = (ValNodePtr PNTR) &(vnp->next);
    }
    vnp = nextvnp;
  }
}

/* name processing code from Sequin editors */

NLM_EXTERN void FirstNameToInitials (CharPtr first, CharPtr inits, size_t maxsize)

{
  Char  ch;
  Uint2  i;

  if (inits != NULL && maxsize > 0) {
    inits [0] = '\0';
    if (first != NULL) {
      i = 0;
      ch = *first;
      while (ch != '\0' && i < maxsize) {
        while (ch != '\0' && (ch <= ' ' || ch == '-')) {
          first++;
          ch = *first;
        }
        if (IS_ALPHA (ch)) {
          inits [i] = ch;
          i++;
          first++;
          ch = *first;
        }
        while (ch != '\0' && ch > ' ' && ch != '-') {
          first++;
          ch = *first;
        }
        if (ch == '-') {
          inits [i] = ch;
          i++;
          first++;
          ch = *first;
        }
      }
      inits [i] = '\0';
    }
  }
}

static void StripPeriods (CharPtr str)

{
  Char     ch;
  CharPtr  dst;

  if (str != NULL) {
    dst = str;
    ch = *str;
    while (ch != '\0') {
      if (ch != '.') {
        *dst = ch;
        dst++;
      }
      str++;
      ch = *str;
    }
    *dst = '\0';
  }
}

static void TrimLeadingSpaces (CharPtr str)

{
  Char     ch;
  CharPtr  dst;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ch = *str;
    while (ch != '\0' && ch <= ' ') {
      str++;
      ch = *str;
    }
    while (ch != '\0') {
      *dst = ch;
      dst++;
      str++;
      ch = *str;
    }
    *dst = '\0';
  }
}

static void ExtractSuffixFromInitials (NameStdPtr nsp)

{
  Char     ch;
  Boolean  has_period = FALSE;
  size_t   len;
  CharPtr  str;

  str = nsp->names [4];
  ch = *str;
  while (ch != '\0') {
    if (ch == '.') {
      has_period = TRUE;
    }
    str++;
    ch = *str;
  }
  if (! has_period) return;
  str = nsp->names [4];
  len = StringLen (str);
  if (len >= 4 && StringCmp (str +  len - 3, "III") == 0) {
    str [len - 3] = '\0';
    nsp->names [5] = StringSave ("III");
  } else if (len >= 5 && StringCmp (str +  len - 4, "III.") == 0) {
    str [len - 4] = '\0';
    nsp->names [5] = StringSave ("III");
  } else if (len >= 3 && StringCmp (str +  len - 2, "Jr") == 0) {
    str [len - 2] = '\0';
    nsp->names [5] = StringSave ("Jr");
  } else if (len >= 4 && StringCmp (str +  len - 3, "2nd") == 0) {
    str [len - 3] = '\0';
    nsp->names [5] = StringSave ("II");
  } else if (len >= 3 && StringCmp (str +  len - 2, "IV") == 0) {
    str [len - 2] = '\0';
    nsp->names [5] = StringSave ("IV");
  } else if (len >= 4 && StringCmp (str +  len - 3, "IV.") == 0) {
    str [len - 3] = '\0';
    nsp->names [5] = StringSave ("IV");
  }
}

static CharPtr NameStdPtrToTabbedString (NameStdPtr nsp, Boolean fixInitials)

{
  Char   first [256];
  Char   frstinits [64];
  Char   initials [64];
  Int2   j;
  Char   last [256];
  Char   middle [128];
  Char   str [512];
  Char   suffix [64];

  if (nsp == NULL) return NULL;
  if (nsp->names [5] == NULL && nsp->names [4] != NULL) {
    ExtractSuffixFromInitials (nsp);
  }
  str [0] = '\0';
  StringNCpy_0 (first, nsp->names [1], sizeof (first));
  TrimSpacesAroundString (first);
  StringNCpy_0 (initials, nsp->names [4], sizeof (initials));
  StripPeriods (initials);
  TrimLeadingSpaces (initials);
  StringNCpy_0 (last, nsp->names [0], sizeof (last));
  TrimLeadingSpaces (last);
  StringNCpy_0 (middle, nsp->names [2], sizeof (middle));
  TrimLeadingSpaces (middle);
  if (StringCmp (initials, "al") == 0 &&
      StringCmp (last, "et") == 0 &&
      first [0] == '\0') {
    initials [0] = '\0';
    StringCpy (last, "et al.");
  }
  /*
  if (first [0] == '\0') {
    StringNCpy_0 (first, initials, sizeof (first));
    if (IS_ALPHA (first [0])) {
      if (first [1] == '-') {
        first [3] = '\0';
      } else {
        first [1] = '\0';
      }
    } else {
      first [0] = '\0';
    }
  }
  */
  frstinits [0] = '\0';
  FirstNameToInitials (first, frstinits, sizeof (frstinits) - 1);
  StripPeriods (first);
  TrimLeadingSpaces (first);
  if (first [0] != '\0') {
    StringCat (str, first);
  } else {
    /*
    StringCat (str, " ");
    */
  }
  StringCat (str, "\t");
  if (fixInitials) {
    j = 0;
    while (initials [j] != '\0' && TO_UPPER (initials [j]) == TO_UPPER (frstinits [j])) {
      j++;
    }
    if (initials [j] != '\0') {
      StringCat (str, initials + j);
    } else {
      /*
      StringCat (str, " ");
      */
    }
  } else if (initials [0] != '\0') {
    StringCat (str, initials);
  } else if (frstinits [0] != '\0') {
    StringCat (str, frstinits);
  }
  StringCat (str, "\t");
  StringCat (str, last);
  StringNCpy_0 (suffix, nsp->names [5], sizeof (suffix));
  StringCat (str, "\t");
  StripPeriods (suffix);
  TrimLeadingSpaces (suffix);
  if (suffix [0] != '\0') {
    StringCat (str, suffix);
  } else {
    /*
    StringCat (str, " ");
    */
  }
  StringCat (str, "\t");
  StringCat (str, middle);
  StringCat (str, "\n");
  return StringSave (str);
}

static CharPtr XtractTagListColumn (CharPtr source, Int2 col)

{
  Char     ch;
  size_t   count;
  CharPtr  ptr;
  CharPtr  str;

  if (source == NULL || source [0] == '\0' || col < 0) return NULL;

  ptr = source;
  ch = *ptr;
  while (col > 0 && ch != '\n' && ch != '\0') {
    while (ch != '\t' && ch != '\n' && ch != '\0') {
      ptr++;
      ch = *ptr;
    }
    if (ch == '\t') {
      ptr++;
      ch = *ptr;
    }
    col--;
  }

  count = 0;
  ch = ptr [count];
  while (ch != '\t' && ch != '\n' && ch != '\0') {
    count++;
    ch = ptr [count];
  }
  str = (CharPtr) MemNew(count + 1);
  if (str != NULL) {
    MemCpy (str, ptr, count);
  }
  return str;
}

static NameStdPtr TabbedStringToNameStdPtr (CharPtr txt, Boolean fixInitials)

{
  Char        ch;
  CharPtr     first;
  Char        initials [64];
  Int2        j;
  Int2        k;
  Char        last;
  Int2        len;
  NameStdPtr  nsp;
  Char        periods [128];
  CharPtr     str;
  Char        str1 [64];
  Char        suffix [80];

  if (txt == NULL) return NULL;
  nsp = NameStdNew ();
  if (nsp == NULL) return NULL;
  nsp->names [0] = XtractTagListColumn (txt, 2);
  TrimLeadingSpaces (nsp->names [0]);
  first = XtractTagListColumn (txt, 0);
  StripPeriods (first);
  nsp->names [1] = StringSave (first);
  TrimLeadingSpaces (nsp->names [1]);
  str1 [0] = '\0';
  if (fixInitials) {
    FirstNameToInitials (first, str1, sizeof (str1) - 1);
  }
  str = XtractTagListColumn (txt, 1);
  StringNCat (str1, str, sizeof (str1) - 1);
  MemFree (str);
  j = 0;
  k = 0;
  ch = str1 [j];
  while (ch != '\0') {
    if (ch != ' ') {
      initials [k] = ch;
      k++;
    }
    j++;
    ch = str1 [j];
  }
  initials [k] = '\0';
  periods [0] = '\0';
          j = 0;
          ch = initials [j];
          while (ch != '\0') {
            if (ch == ',') {
              initials [j] = '.';
            }
            j++;
            ch = initials [j];
          }
          str = StringStr (initials, ".ST.");
          if (str != NULL) {
            *(str + 2) = 't';
          }
  j = 0;
  k = 0;
  ch = initials [j];
  while (ch != '\0') {
    if (ch == '-') {
      periods [k] = ch;
      k++;
      j++;
      ch = initials [j];
    } else if (ch == '.') {
      j++;
      ch = initials [j];
            } else if (ch == ' ') {
              j++;
              ch = initials [j];
    } else {
      periods [k] = ch;
              last = ch;
      k++;
      j++;
      ch = initials [j];
              if (ch == '\0') {
                if (! (IS_LOWER (last))) {
                  periods [k] = '.';
                  k++;
                }
              /* } else if (ch == '.' && initials [j + 1] == '\0') { */
              } else if (! (IS_LOWER (ch))) {
                periods [k] = '.';
                k++;
              }
    }
  }
  if (k > 0 && periods [k - 1] != '.') {
    periods [k] = '.';
    k++;
  }
  periods [k] = '\0';
  nsp->names [4] = StringSave (periods);
  TrimLeadingSpaces (nsp->names [4]);
  str = XtractTagListColumn (txt, 3);
  StringNCpy_0 (str1, str, sizeof (str1));
  MemFree (str);
  j = 0;
  k = 0;
  ch = str1 [j];
  while (ch != '\0') {
    if (ch != ' ') {
      suffix [k] = ch;
      k++;
    }
    j++;
    ch = str1 [j];
  }
  suffix [k] = '\0';
  if (suffix [0] != '\0') {
    len = StringLen (suffix);
    if (len > 0 && suffix [len - 1] == '.') {
      suffix [len - 1] = '\0';
    }
    if (StringICmp (suffix, "1d") == 0) {
      StringCpy (suffix, "I");
    } else if (StringICmp (suffix, "1st") == 0) {
      StringCpy (suffix, "I");
    } else if (StringICmp (suffix, "2d") == 0) {
      StringCpy (suffix, "II");
    } else if (StringICmp (suffix, "2nd") == 0) {
      StringCpy (suffix, "II");
    } else if (StringICmp (suffix, "3d") == 0) {
      StringCpy (suffix, "III");
    } else if (StringICmp (suffix, "3rd") == 0) {
      StringCpy (suffix, "III");
    } else if (StringICmp (suffix, "4th") == 0) {
      StringCpy (suffix, "IV");
    } else if (StringICmp (suffix, "5th") == 0) {
      StringCpy (suffix, "V");
    } else if (StringICmp (suffix, "6th") == 0) {
      StringCpy (suffix, "VI");
    } else if (StringICmp (suffix, "Sr") == 0) {
      StringCpy (suffix, "Sr.");
    } else if (StringICmp (suffix, "Jr") == 0) {
      StringCpy (suffix, "Jr.");
    }
    /*
    len = StringLen (suffix);
    if (len > 0 && suffix [len - 1] != '.') {
      StringCat (suffix, ".");
    }
    */
    nsp->names [5] = StringSave (suffix);
    TrimLeadingSpaces (nsp->names [5]);
  }
  if (StringCmp (nsp->names [0], "et al") == 0) {
    nsp->names [0] = MemFree (nsp->names [0]);
    nsp->names [0] = StringSave ("et al.");
  }
  nsp->names [2] = XtractTagListColumn (txt, 4);
  TrimLeadingSpaces (nsp->names [2]);
  if (StringHasNoText (nsp->names [2])) {
    nsp->names [2] = MemFree (nsp->names [2]);
  }
  MemFree (first);
  return nsp;
}

static AffilPtr CleanAffil (AffilPtr afp)

{
  if (afp == NULL) return NULL;
  CleanVisStringJunkAndCompress (&(afp->affil));
  CleanVisStringJunkAndCompress (&(afp->div));
  CleanVisStringJunkAndCompress (&(afp->city));
  CleanVisStringJunkAndCompress (&(afp->sub));
  CleanVisStringJunkAndCompress (&(afp->country));
  CleanVisStringJunkAndCompress (&(afp->street));
  CleanVisStringJunkAndCompress (&(afp->email));
  CleanVisStringJunkAndCompress (&(afp->fax));
  CleanVisStringJunkAndCompress (&(afp->phone));
  CleanVisStringJunkAndCompress (&(afp->postal_code));
  if (afp->choice == 2) {
    if (StringCmp (afp->country, "U.S.A.") == 0) {
      afp->country = MemFree (afp->country);
      afp->country = StringSave ("USA");
    }
  }
  if (afp->affil == NULL &&
      afp->div == NULL &&
      afp->city == NULL &&
      afp->sub == NULL &&
      afp->country == NULL &&
      afp->street == NULL &&
      afp->email == NULL &&
      afp->fax == NULL &&
      afp->phone == NULL &&
      afp->postal_code == NULL) {
    afp = MemFree (afp);
  }
  return afp;
}

static void NormalizeAuthors (AuthListPtr alp, Boolean fixInitials)

{
  AuthorPtr        ap;
  CharPtr          initials;
  size_t           len;
  ValNodePtr       names;
  ValNodePtr       next;
  NameStdPtr       nsp;
  PersonIdPtr      pid;
  ValNodePtr PNTR  prev;
  CharPtr          str;
  Boolean          upcaseinits;
  ValNodePtr       vnp;
  Boolean          zap;

  if (alp == NULL) return;
  alp->affil = CleanAffil (alp->affil);

  if (alp->choice == 2 || alp->choice == 3) {
    for (vnp = alp->names; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      TrimSpacesAroundString (str);
      Asn2gnbkCompressSpaces (str);
    }
  }
  if (alp->choice != 1) return;

  prev = &(alp->names);
  names = alp->names;
  while (names != NULL) {
    next = names->next;
    zap = FALSE;
    ap = names->data.ptrvalue;
    if (ap != NULL) {
      pid = ap->name;
      if (pid == NULL) {
        /* continue */
      } else if (pid->choice == 2) {
        nsp = pid->data;
        if (nsp != NULL /* && nsp->names [4] != NULL */) {
          upcaseinits = FALSE;
          initials = nsp->names [4];
          if (StringLen (initials) > 0) {
            if (IS_UPPER (initials [0])) {
              upcaseinits = TRUE;
            }
          }
          str = NameStdPtrToTabbedString (nsp, fixInitials);
          pid->data = NameStdFree (nsp);
          nsp = TabbedStringToNameStdPtr (str, fixInitials);
          if (upcaseinits) {
            initials = nsp->names [4];
            if (StringLen (initials) > 0) {
              if (IS_LOWER (initials [0])) {
                initials [0] = TO_UPPER (initials [0]);
              }
            }
          }
          pid->data = nsp;
          MemFree (str);
          CleanVisString (&(nsp->names [0]));
          CleanVisString (&(nsp->names [1]));
          CleanVisString (&(nsp->names [2]));
          CleanVisString (&(nsp->names [3]));
          CleanVisString (&(nsp->names [4]));
          CleanVisString (&(nsp->names [5]));
          CleanVisString (&(nsp->names [6]));
          if (StringCmp (nsp->names [0], "et") == 0 &&
              (StringCmp (nsp->names [4], "al") == 0 ||
               StringCmp (nsp->names [4], "al.") == 0 ||
               StringCmp (nsp->names [4], "Al.") == 0) &&
              (StringHasNoText (nsp->names [1]) ||
               StringCmp (nsp->names [1], "a") == 0)) {
            nsp->names [4] = MemFree (nsp->names [4]);
            nsp->names [1] = MemFree (nsp->names [1]);
            nsp->names [0] = MemFree (nsp->names [0]);
            nsp->names [0] = StringSave ("et al.");
          }
          str = nsp->names [0];
          len = StringLen (str);
          if (len > 4 && StringHasNoText (nsp->names [5])) {
            if (StringCmp (str + len - 4, " Jr.") == 0 ||
                StringCmp (str + len - 4, " Sr.") == 0) {
              nsp->names [5] = StringSave (str + len - 3);
              str [len - 4] = '\0';
              TrimSpacesAroundString (str);
            }
          }
          str = nsp->names [4];
          len = StringLen (str);
          if (len > 4 && StringHasNoText (nsp->names [5])) {
            if (StringCmp (str + len - 4, ".Jr.") == 0 ||
                StringCmp (str + len - 4, ".Sr.") == 0) {
              nsp->names [5] = StringSave (str + len - 3);
              str [len - 3] = '\0';
              TrimSpacesAroundString (str);
            }
          }
          if (StringHasNoText (nsp->names [0]) &&
              StringHasNoText (nsp->names [1]) &&
              StringHasNoText (nsp->names [2]) &&
              StringHasNoText (nsp->names [3]) &&
              StringHasNoText (nsp->names [4]) &&
              StringHasNoText (nsp->names [5]) &&
              StringHasNoText (nsp->names [6])) {
            zap = TRUE;
          }
          /* last name is required, so zap if not present */
          if (StringHasNoText (nsp->names [0])) {
            zap = TRUE;
          }
        }
      } else if (pid->choice == 3 || pid->choice == 4 || pid->choice == 5) {
        TrimSpacesAroundString ((CharPtr) pid->data);
        if (StringHasNoText ((CharPtr) pid->data)) {
          zap = TRUE;
        }
      }
    }
    if (zap) {
      /* remove empty authors */
      *prev = names->next;
      names->next = NULL;
      AuthorFree (ap);
      ValNodeFree (names);
    } else {
      prev = &(names->next);
    }
    names = next;
  }
  /* if no remaining authors, put in default author for legal ASN.1 */
  if (alp->names == NULL) {
    names = ValNodeNew (NULL);
    if (names != NULL) {
      /*
      ap = AuthorNew ();
      if (ap != NULL) {
        pid = PersonIdNew ();
        if (pid != NULL) {
          pid->choice = 4;
          pid->data = (Pointer) StringSave ("?");
          ap->name = pid;
          names->choice = 1;
          names->data.ptrvalue = ap;
          alp->names = names;
        }
      }
      */
      names->choice = 3;
      names->data.ptrvalue = StringSave ("?");
      alp->names = names;
      alp->choice = 3;
    }
  }
}

static void StrStripSpaces (
  CharPtr str
)

{
  CharPtr  new_str;

  if (str == NULL) return;

  new_str = str;
  while (*str != '\0') {
    *new_str++ = *str;
    if (*str == ' ' || *str == '\t' || *str == '(') {
      for (str++; *str == ' ' || *str == '\t'; str++) continue;
      if (*str == ')' || *str == ',') {
        new_str--;
      }
    } else {
      str++;
    }
  }
  *new_str = '\0';
}

/* from utilpub.c */
static Boolean empty_citgen(CitGenPtr  cit)
{
    if (cit == NULL)
        return TRUE;
    if (cit->cit)
        return FALSE;
    if (cit->authors)
        return FALSE;
    if (cit->muid > 0)
        return FALSE;
    if (cit->journal)
        return FALSE;
    if (cit->volume)
        return FALSE;
    if (cit->issue)
        return FALSE;
    if (cit->pages)
        return FALSE;
    if (cit->date)
        return FALSE;
    if (cit->serial_number > 0)
        return FALSE;
    if (cit->title)
        return FALSE;
    if (cit->pmid > 0)
        return FALSE;
    return TRUE;
}

static void NormalizePubAuthors (ValNodePtr vnp, Boolean stripSerial, Boolean fixInitials)

{
  CitArtPtr    cap;
  CitBookPtr   cbp;
  CitGenPtr    cgp;
  CitPatPtr    cpp;
  CitSubPtr    csp;

  if (vnp == NULL) return;
  if (vnp->choice == PUB_PMid || vnp->choice == PUB_Muid) return;
  if (vnp->data.ptrvalue == NULL) return;
  switch (vnp->choice) {
    case PUB_Gen :
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      NormalizeAuthors (cgp->authors, fixInitials);
      break;
    case PUB_Sub :
      csp = (CitSubPtr) vnp->data.ptrvalue;
      NormalizeAuthors (csp->authors, fixInitials);
      break;
    case PUB_Article :
      cap = (CitArtPtr) vnp->data.ptrvalue;
      NormalizeAuthors (cap->authors, fixInitials);
      if (cap->from == 2 || cap->from == 3) {
        cbp = (CitBookPtr) cap->fromptr;
        if (cbp != NULL) {
          NormalizeAuthors (cbp->authors, fixInitials);
        }
      }
      break;
    case PUB_Book :
      cbp = (CitBookPtr) vnp->data.ptrvalue;
      NormalizeAuthors (cbp->authors, fixInitials);
      break;
    case PUB_Man :
      cbp = (CitBookPtr) vnp->data.ptrvalue;
      if (cbp->othertype == 2 && cbp->let_type == 3) {
        NormalizeAuthors (cbp->authors, fixInitials);
      }
      break;
    case PUB_Patent :
      cpp = (CitPatPtr) vnp->data.ptrvalue;
      NormalizeAuthors (cpp->authors, fixInitials);
      NormalizeAuthors (cpp->applicants, fixInitials);
      NormalizeAuthors (cpp->assignees, fixInitials);
      break;
    default :
      break;
  }
}

static void NormalizeAPub (ValNodePtr vnp, Boolean stripSerial, Boolean fixInitials)

{
  AffilPtr     affil;
  AuthListPtr  alp;
  CitArtPtr    cap;
  CitBookPtr   cbp;
  CitGenPtr    cgp;
  CitJourPtr   cjp;
  CitPatPtr    cpp;
  CitSubPtr    csp;
  ImprintPtr   imp;
  CharPtr      str;
  CharPtr      tmp;
  ValNodePtr   ttl;

  if (vnp == NULL) return;
  if (vnp->choice == PUB_PMid || vnp->choice == PUB_Muid) return;
  if (vnp->data.ptrvalue == NULL) return;
  imp = NULL;
  switch (vnp->choice) {
    case PUB_Gen :
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      if (stripSerial) {
        cgp->serial_number = -1; /* but does not remove if empty */
      }
      if (StringNICmp (cgp->cit, "unpublished", 11) == 0) {
        cgp->cit [0] = 'U';
        /* cgp->date = DateFree (cgp->date); */ /* remove date if unpublished */
        if (cgp->journal == NULL) {
          cgp->volume = MemFree (cgp->volume);
          cgp->issue = MemFree (cgp->issue);
          cgp->pages = MemFree (cgp->pages);
        }
      }
      TrimSpacesAroundString (cgp->cit);
      if (StringDoesHaveText (cgp->title)) {
        StrStripSpaces (cgp->title);
      }
      break;
    case PUB_Sub :
      csp = (CitSubPtr) vnp->data.ptrvalue;
      alp = csp->authors;
      imp = csp->imp;
      if (alp != NULL && alp->affil == NULL && imp != NULL && imp->pub != NULL) {
        alp->affil = imp->pub;
        imp->pub = NULL;
      }
      if (csp->date == NULL && imp != NULL && imp->date != NULL) {
        csp->date = imp->date;
        imp->date = NULL;
      }
      if (imp != NULL && imp->date == NULL) {
        csp->imp = ImprintFree (csp->imp);
      }
      if (alp != NULL && alp->affil != NULL) {
        affil = alp->affil;
        if (affil->choice == 1) {
          str = affil->affil;
          if (StringNICmp (str, "to the ", 7) == 0) {
            if (StringNICmp (str + 24, " databases", 10) == 0) {
              str += 34;
              if (*str == '.') {
                str++;
              }
              tmp = StringSaveNoNull (TrimSpacesAroundString (str));
              affil->affil = MemFree (affil->affil);
              affil->affil = tmp;
            }
          }
        }
        alp->affil = CleanAffil (alp->affil);
      }
      imp = csp->imp;
      break;
    case PUB_Article :
      cap = (CitArtPtr) vnp->data.ptrvalue;
      if (cap != NULL) {
        if (cap->from == 1) {
          cjp = (CitJourPtr) cap->fromptr;
          if (cjp != NULL) {
            imp = cjp->imp;
          }
        } else if (cap->from == 2 || cap->from == 3) {
          cbp = (CitBookPtr) cap->fromptr;
          if (cbp != NULL) {
            imp = cbp->imp;
          }
        }
        for (ttl = cap->title; ttl != NULL; ttl = ttl->next) {
          if (ttl->choice == Cit_title_name) {
            str = (CharPtr) ttl->data.ptrvalue;
            if (StringHasNoText (str)) continue;
            StrStripSpaces (str);
          }
        }
      }
      break;
    case PUB_Patent :
      cpp = (CitPatPtr) vnp->data.ptrvalue;
      if (cpp != NULL) {
        if (StringCmp (cpp->country, "USA") == 0) {
          cpp->country = MemFree (cpp->country);
          cpp->country = StringSave ("US");
        }
      }
      break;
    default :
      break;
  }
  if (imp != NULL) {
    CleanVisStringAndCompress (&(imp->volume));
    CleanVisStringAndCompress (&(imp->issue));
    CleanVisStringAndCompress (&(imp->pages));
    CleanVisStringAndCompress (&(imp->section));
    CleanVisStringAndCompress (&(imp->part_sup));
    CleanVisStringAndCompress (&(imp->language));
    CleanVisStringAndCompress (&(imp->part_supi));
  }
}

NLM_EXTERN void CleanUpPubdescAuthors (PubdescPtr pdp)

{
  Char             buf1 [121];
  Boolean          fixInitials = TRUE;
  Boolean          hasArt = FALSE;
  Boolean          hasUid = FALSE;
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       vnp;

  if (pdp == NULL) return;
  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Muid || vnp->choice == PUB_PMid) {
      if (vnp->data.intvalue > 0) {
        hasUid = TRUE;
      }
    } else if (vnp->choice == PUB_Article) {
      hasArt = TRUE;
    }
  }
  if (hasArt && hasUid) {
    fixInitials = FALSE;
  }
  prev = &(pdp->pub);
  vnp = pdp->pub;
  while (vnp != NULL) {
    next = vnp->next;
    PubLabelUnique (vnp, buf1, sizeof (buf1) - 1, OM_LABEL_CONTENT, TRUE);
    NormalizePubAuthors (vnp, TRUE, fixInitials);
    vnp = next;
  }
}

static void NormalizePubdesc (PubdescPtr pdp, Boolean stripSerial, Boolean doAuthors, ValNodePtr PNTR publist)

{
  ArticleIdPtr     aip;
  Int4             artpmid = 0;
  Char             buf1 [121];
  Char             buf2 [121];
  CitArtPtr        cap = NULL;
  CitGenPtr        cgp;
  CitJourPtr       cjp;
  Boolean          fixInitials = TRUE;
  Boolean          hasArt = FALSE;
  Boolean          hasUid = FALSE;
  ImprintPtr       imp;
  ValNodePtr       next;
  Int4             pmid = 0;
  ValNodePtr PNTR  prev;
  ValNodePtr       vnp;

  if (pdp == NULL) return;
  CleanVisString (&(pdp->comment));
  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Muid || vnp->choice == PUB_PMid) {
      if (vnp->data.intvalue > 0) {
        hasUid = TRUE;
      }
    } else if (vnp->choice == PUB_Article) {
      hasArt = TRUE;
    }
  }
  if (hasArt && hasUid) {
    fixInitials = FALSE;
  }
  prev = &(pdp->pub);
  vnp = pdp->pub;
  if (vnp != NULL && vnp->next == NULL && vnp->choice == PUB_Gen) {
    cgp = (CitGenPtr) vnp->data.ptrvalue;
    buf1 [0] = '\0';
    PubLabelUnique (vnp, buf1, sizeof (buf1) - 1, OM_LABEL_CONTENT, TRUE);
    if (doAuthors) {
      NormalizeAuthors (cgp->authors, fixInitials);
    }
    if (stripSerial) {
      cgp->serial_number = -1;
    }
    if (StringNICmp (cgp->cit, "unpublished", 11) == 0) {
      cgp->cit [0] = 'U';
      /* cgp->date = DateFree (cgp->date); */ /* remove date if unpublished */
      if (cgp->journal == NULL) {
        cgp->volume = MemFree (cgp->volume);
        cgp->issue = MemFree (cgp->issue);
        cgp->pages = MemFree (cgp->pages);
      }
    }
    TrimSpacesAroundString (cgp->cit);
    if (StringDoesHaveText (cgp->title)) {
      StrStripSpaces (cgp->title);
    }
    buf2 [0] = '\0';
    PubLabelUnique (vnp, buf2, sizeof (buf2) - 1, OM_LABEL_CONTENT, TRUE);
    if (StringCmp (buf1, buf2) != 0) {
      ValNodeCopyStr (publist, 1, buf1);
      ValNodeCopyStr (publist, 2, buf2);
    }
    return; /* but does not remove if empty and only element of Pub */
  }
  while (vnp != NULL) {
    next = vnp->next;
    buf1 [0] = '\0';
    PubLabelUnique (vnp, buf1, sizeof (buf1) - 1, OM_LABEL_CONTENT, TRUE);
    if (doAuthors) {
      NormalizePubAuthors (vnp, stripSerial, fixInitials);
    }
    NormalizeAPub (vnp, stripSerial, fixInitials);
    if (vnp->choice == PUB_Article) {
      cap = (CitArtPtr) vnp->data.ptrvalue;
      if (cap != NULL && cap->from == 1) {
        cjp = (CitJourPtr) cap->fromptr;
        if (cjp != NULL) {
          imp = cjp->imp;
          if (imp != NULL) {
            if (imp->pubstatus == PUBSTATUS_aheadofprint && imp->prepub != 2) {
              if (StringHasNoText (imp->volume) || StringHasNoText (imp->pages)) {
                imp->prepub = 2;
              }
            }
            if (imp->pubstatus == PUBSTATUS_aheadofprint && imp->prepub == 2) {
              if (StringDoesHaveText (imp->volume) && StringDoesHaveText (imp->pages)) {
                imp->prepub = 0;
              }
            }
            if (imp->pubstatus == PUBSTATUS_epublish && imp->prepub == 2) {
              imp->prepub = 0;
            }
          }
        }
        for (aip = cap->ids; aip != NULL; aip = aip->next) {
          if (aip->choice == ARTICLEID_PUBMED) {
            artpmid = aip->data.intvalue;
          }
        }
      }
    } else if (vnp->choice == PUB_PMid) {
      pmid = vnp->data.intvalue;
    }
    if (vnp->choice == PUB_Gen && empty_citgen ((CitGenPtr) vnp->data.ptrvalue)) {
      *prev = vnp->next;
      vnp->next = NULL;
      PubFree (vnp);
    } else {
      prev = &(vnp->next);
      buf2 [0] = '\0';
      PubLabelUnique (vnp, buf2, sizeof (buf2) - 1, OM_LABEL_CONTENT, TRUE);
      if (StringCmp (buf1, buf2) != 0) {
        ValNodeCopyStr (publist, 1, buf1);
        ValNodeCopyStr (publist, 2, buf2);
      }
    }
    vnp = next;
  }
  if (pmid == 0 && artpmid > 0) {
    ValNodeAddInt (&(pdp->pub), PUB_PMid, artpmid);
  } else if (pmid > 0 && artpmid == 0 && cap != NULL) {
    ValNodeAddInt (&(cap->ids), ARTICLEID_PUBMED, pmid);
  }
}

NLM_EXTERN void CleanUpPubdescBody (PubdescPtr pdp, Boolean stripSerial)

{
  if (pdp == NULL) return;
  NormalizePubdesc (pdp, stripSerial, FALSE, NULL);
}

static Boolean KeywordAlreadyInList (ValNodePtr head, CharPtr kwd)

{
  ValNodePtr  vnp;

  if (head == NULL || kwd == NULL) return FALSE;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    if (StringICmp ((CharPtr) vnp->data.ptrvalue, kwd) == 0) return TRUE;
  }

  return FALSE;
}

static Boolean CopyGeneXrefToGeneFeat (GeneRefPtr grp, GeneRefPtr grx)

{
  if (grp == NULL || grx == NULL) return FALSE;
  if (grx->db != NULL) {
    ValNodeLink (&(grp->db), grx->db);
    grx->db = NULL;
  }
  if (grx->locus == NULL && grx->allele == NULL &&
      grx->desc == NULL && grx->maploc == NULL &&
      grx->locus_tag == NULL && grx->db == NULL &&
      grx->syn == NULL) return TRUE;
  return FALSE;
}

static void HandleXrefOnGene (SeqFeatPtr sfp)

{
  GeneRefPtr           grp;
  GeneRefPtr           grx;
  SeqFeatXrefPtr       next;
  SeqFeatXrefPtr PNTR  prev;
  SeqFeatXrefPtr       xref;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE) return;
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return;
   prev = &(sfp->xref);
  xref = sfp->xref;
  while (xref != NULL) {
    next = xref->next;
    if (xref->data.choice == SEQFEAT_GENE) {
      grx = (GeneRefPtr) xref->data.value.ptrvalue;
      if (CopyGeneXrefToGeneFeat (grp, grx)) {
        *(prev) = next;
        xref->next = NULL;
        SeqFeatXrefFree (xref);
      } else {
        prev = &(xref->next);
      }
    } else {
      prev = &(xref->next);
    }
    xref = next;
  }
}

static void CopyProtXrefToProtFeat (ProtRefPtr prp, ProtRefPtr prx)

{
  ValNodePtr       curr;
  size_t           len;
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  CharPtr          str;

  if (prp == NULL || prx == NULL) return;

  if (prx->db != NULL) {
    ValNodeLink (&(prp->db), prx->db);
    prx->db = NULL;
  }

  prev = &(prx->name);
  curr = prx->name;
  while (curr != NULL) {
    next = curr->next;
    str = (CharPtr) curr->data.ptrvalue;
    if (! KeywordAlreadyInList (prp->name, str)) {
      ValNodeCopyStr (&(prp->name), 0, str);
      *(prev) = next;
      curr->next = NULL;
      curr->data.ptrvalue = NULL;
      ValNodeFree (curr);
    } else {
      prev = &(curr->next);
    }
    curr = next;
  }

  if (prp->desc == NULL) {
    prp->desc = prx->desc;
    prx->desc = NULL;
  } else if (prx->desc != NULL) {
    if (StringCmp (prx->desc, prp->desc) != 0) {
      len = StringLen (prp->desc) + StringLen (prx->desc) + 6;
      str = MemNew (len);
      if (str != NULL) {
        StringCpy (str, prp->desc);
        StringCat (str, "; ");
        StringCat (str, prx->desc);
        prp->desc = MemFree (prp->desc);
        prp->desc = str;
      }
    }
  }

  prev = &(prx->ec);
  curr = prx->ec;
  while (curr != NULL) {
    next = curr->next;
    str = (CharPtr) curr->data.ptrvalue;
    if (! KeywordAlreadyInList (prp->ec, str)) {
      ValNodeCopyStr (&(prp->ec), 0, str);
      *(prev) = next;
      curr->next = NULL;
      curr->data.ptrvalue = NULL;
      ValNodeFree (curr);
    } else {
      prev = &(curr->next);
    }
    curr = next;
  }

  prev = &(prx->activity);
  curr = prx->activity;
  while (curr != NULL) {
    next = curr->next;
    str = (CharPtr) curr->data.ptrvalue;
    if (! KeywordAlreadyInList (prp->activity, str)) {
      ValNodeCopyStr (&(prp->activity), 0, str);
      curr->data.ptrvalue = NULL;
    }
    *(prev) = next;
    curr->next = NULL;
    curr->data.ptrvalue = NULL;
    ValNodeFree (curr);
    curr = next;
  }
}

static Boolean InGpsGenomic (SeqFeatPtr sfp)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;

  if (sfp == NULL) return FALSE;
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return FALSE;
  if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) bsp->idx.parentptr;
    while (bssp != NULL) {
      if (bssp->_class == BioseqseqSet_class_nuc_prot) return FALSE;
      if (bssp->_class == BioseqseqSet_class_gen_prod_set) return TRUE;
      if (bssp->idx.parenttype != OBJ_BIOSEQSET) return FALSE;
      bssp = (BioseqSetPtr) bssp->idx.parentptr;
    }
  }
  return FALSE;
}

static void HandleXrefOnCDS (SeqFeatPtr sfp)

{
  SeqFeatXrefPtr       next;
  SeqFeatXrefPtr PNTR  prev;
  SeqFeatPtr           prot;
  ProtRefPtr           prp;
  ProtRefPtr           prx;
  SeqFeatXrefPtr       xref;

  if (sfp != NULL && sfp->product != NULL) {
    if (InGpsGenomic (sfp)) return;
    prot = GetBestProteinFeatureUnindexed (sfp->product);
    if (prot != NULL) {
      prp = (ProtRefPtr) prot->data.value.ptrvalue;
      if (prp != NULL) {
        prev = &(sfp->xref);
        xref = sfp->xref;
        while (xref != NULL) {
          next = xref->next;
          if (xref->data.choice == SEQFEAT_PROT) {
            prx = (ProtRefPtr) xref->data.value.ptrvalue;
            CopyProtXrefToProtFeat (prp, prx);
            *(prev) = next;
            xref->next = NULL;
            SeqFeatXrefFree (xref);
          } else {
            prev = &(xref->next);
          }
          xref = next;
        }
      }
    }
  }
}

static void CleanUserStrings (
  UserFieldPtr ufp,
  Pointer userdata
)

{
  ObjectIdPtr  oip;

  oip = ufp->label;
  if (oip != NULL && oip->str != NULL) {
    if (! StringHasNoText (oip->str)) {
      CleanVisString (&(oip->str));
    }
  }
  if (ufp->choice == 1) {
    if (! StringHasNoText ((CharPtr) ufp->data.ptrvalue)) {
      CleanVisString ((CharPtr PNTR) &(ufp->data.ptrvalue));
    }
  }
}

static void CleanUserFields (
  UserFieldPtr ufp,
  Pointer userdata
)

{
  ObjectIdPtr  oip;

  oip = ufp->label;
  if (oip != NULL && oip->str != NULL) {
    if (! StringHasNoText (oip->str)) {
      CleanVisString (&(oip->str));
    }
  }
  VisitUserFieldsInUfp (ufp, userdata, CleanUserStrings);
}


NLM_EXTERN void CleanStructuredComment (
  UserObjectPtr uop
)

{
  Boolean      genome_assembly_data = FALSE;
  UserFieldPtr ufp;
  CharPtr      str, core, new_str;

  if (uop == NULL || uop->type == NULL 
      || StringCmp (uop->type->str, "StructuredComment") != 0) {
    return;
  }

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    if (ufp->label != NULL 
        && ufp->choice == 1 
        && (str = (CharPtr) ufp->data.ptrvalue) != NULL) {
      if (StringCmp (ufp->label->str, "StructuredCommentPrefix") == 0) {
        core = StructuredCommentDbnameFromString(str);
        new_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (core) + 15));
        sprintf (new_str, "##%s-START##", core);
        str = MemFree (str);
        ufp->data.ptrvalue = new_str;
        if (StringCmp (core, "Genome-Assembly-Data") == 0) {
          genome_assembly_data = TRUE;
        }
        core = MemFree (core);
      } else if (StringCmp (ufp->label->str, "StructuredCommentSuffix") == 0) {
        core = StructuredCommentDbnameFromString(str);
        new_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (core) + 15));
        sprintf (new_str, "##%s-END##", core);
        str = MemFree (str);
        ufp->data.ptrvalue = new_str;
        if (StringCmp (core, "Genome-Assembly-Data") == 0) {
          genome_assembly_data = TRUE;
        }
        core = MemFree (core);
      }
    }
  }

  if (genome_assembly_data) {
    for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
      if (ufp->label != NULL 
          && ufp->choice == 1 
          && (str = (CharPtr) ufp->data.ptrvalue) != NULL) {
        if (StringCmp (ufp->label->str, "Finishing Goal") == 0 ||
            StringCmp (ufp->label->str, "Current Finishing Status") == 0) {
          if (StringCmp (str, "High Quality Draft") == 0) {
            ufp->data.ptrvalue = StringSave ("High-Quality Draft");
            str = MemFree (str);
          } else if (StringCmp (str, "Improved High Quality Draft") == 0) {
            ufp->data.ptrvalue = StringSave ("Improved High-Quality Draft");
            str = MemFree (str);
          } else if (StringCmp (str, "Annotation Directed") == 0) {
            ufp->data.ptrvalue = StringSave ("Annotation-Directed Improvement");
            str = MemFree (str);
          } else if (StringCmp (str, "Non-contiguous Finished") == 0) {
            ufp->data.ptrvalue = StringSave ("Noncontiguous Finished");
            str = MemFree (str);
          }
        }
      }
    }
  }
}


static void CleanUserObject (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr  oip;

  oip = uop->type;
  if (oip != NULL && oip->str != NULL) {
    if (! StringHasNoText (oip->str)) {
      CleanVisString (&(oip->str));
    }
  }
  VisitUserFieldsInUop (uop, userdata, CleanUserFields);
  CleanStructuredComment (uop);
}

static CharPtr bsecSiteList [] = {
  "", "active", "binding", "cleavage", "inhibit", "modifi",
  "glycosylation", "myristoylation", "mutagenized", "metal-binding",
  "phosphorylation", "acetylation", "amidation", "methylation",
  "hydroxylation", "sulfatation", "oxidative-deamination",
  "pyrrolidone-carboxylic-acid", "gamma-carboxyglutamic-acid",
  "blocked", "lipid-binding", "np-binding", "DNA-binding",
  "signal-peptide", "transit-peptide", "transmembrane-region",
  "nitrosylation", NULL
};

static CharPtr uninfStrings [] = {
  "signal",
  "transit",
  "peptide",
  "signal peptide",
  "signal-peptide",
  "signal_peptide",
  "transit peptide",
  "transit-peptide",
  "transit_peptide",
  "unnamed",
  "unknown",
  "putative",
  NULL
};

static Boolean InformativeString (CharPtr str)

{
  Int2  i;

  if (StringHasNoText (str)) return FALSE;

  for (i = 0; uninfStrings [i] != NULL; i++) {
    if (StringICmp (str, uninfStrings [i]) == 0) return FALSE;
  }

  return TRUE;
}

static void CleanUpExceptText (SeqFeatPtr sfp)

{
  ValNodePtr  head, vnp;
  size_t      len;
  CharPtr     prefix, ptr, str, tmp;

  if (sfp == NULL || sfp->except_text == NULL) return;
  if (StringStr (sfp->except_text, "ribosome slippage") == NULL &&
      StringStr (sfp->except_text, "trans splicing") == NULL &&
      StringStr (sfp->except_text, "alternate processing") == NULL &&
      StringStr (sfp->except_text, "non-consensus splice site") == NULL &&
      StringStr (sfp->except_text, "adjusted for low quality genome") == NULL) return;

  head = NULL;
  str = sfp->except_text;
  tmp = str;
  while (! StringHasNoText (tmp)) {
    ptr = StringChr (tmp, ',');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }
    TrimSpacesAroundString (tmp);
    ValNodeCopyStr (&head, 0, tmp);
    tmp = ptr;
  }
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    tmp = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (tmp)) continue;
    if (StringCmp (tmp, "ribosome slippage") == 0) {
      vnp->data.ptrvalue = MemFree (tmp);
      vnp->data.ptrvalue = StringSave ("ribosomal slippage");
    } else if (StringCmp (tmp, "trans splicing") == 0) {
      vnp->data.ptrvalue = MemFree (tmp);
      vnp->data.ptrvalue = StringSave ("trans-splicing");
    } else if (StringCmp (tmp, "alternate processing") == 0) {
      vnp->data.ptrvalue = MemFree (tmp);
      vnp->data.ptrvalue = StringSave ("alternative processing");
    } else if (StringCmp (tmp, "non-consensus splice site") == 0) {
      vnp->data.ptrvalue = MemFree (tmp);
      vnp->data.ptrvalue = StringSave ("nonconsensus splice site");
    } else if (StringCmp (tmp, "adjusted for low quality genome") == 0) {
      vnp->data.ptrvalue = MemFree (tmp);
      vnp->data.ptrvalue = StringSave ("adjusted for low-quality genome");
    }
  }

  len = 0;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    tmp = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (tmp)) continue;
    len += StringLen (tmp) + 2;
  }

  str = (CharPtr) MemNew (len + 2);
  if (str == NULL) return;

  prefix = "";
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    tmp = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (tmp)) continue;
    StringCat (str, prefix);
    StringCat (str, tmp);
    prefix = ", ";
  }

  sfp->except_text = MemFree (sfp->except_text);
  sfp->except_text = str;

  ValNodeFreeData (head);
}

static Boolean ExpandGeneSynCom (ValNodePtr headsyn)

{
  ValNodePtr  lastsyn;
  ValNodePtr  newsyn;
  ValNodePtr  nextsyn;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     tmp;

  str = (CharPtr) headsyn->data.ptrvalue;
  if (StringHasNoText (str)) return TRUE;
  if (StringChr (str, ',') == NULL) return FALSE;

  nextsyn = headsyn->next;
  lastsyn = headsyn;
  tmp = StringSave ((CharPtr) headsyn->data.ptrvalue);
  str = tmp;

  while (! StringHasNoText (str)) {
    ptr = StringChr (str, ',');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }
    TrimSpacesAroundString (str);
    newsyn = ValNodeNew (NULL);
    if (newsyn != NULL) {
      newsyn->data.ptrvalue = StringSave (str);
      newsyn->next = nextsyn;
      lastsyn->next = newsyn;
      lastsyn = newsyn;
    }
    str = ptr;
  }

  MemFree (tmp);
  return TRUE;
}

static Boolean ExpandGeneSynSem (ValNodePtr headsyn)

{
  ValNodePtr  lastsyn;
  ValNodePtr  newsyn;
  ValNodePtr  nextsyn;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     tmp;

  str = (CharPtr) headsyn->data.ptrvalue;
  if (StringHasNoText (str)) return TRUE;
  if (StringStr (str, "; ") == NULL) return FALSE;

  nextsyn = headsyn->next;
  lastsyn = headsyn;
  tmp = StringSave ((CharPtr) headsyn->data.ptrvalue);
  str = tmp;

  while (! StringHasNoText (str)) {
    ptr = StringStr (str, "; ");
    if (ptr != NULL) {
      ptr++;
      *ptr = '\0';
      ptr++;
    }
    TrimSpacesAroundString (str);
    newsyn = ValNodeNew (NULL);
    if (newsyn != NULL) {
      newsyn->data.ptrvalue = StringSave (str);
      newsyn->next = nextsyn;
      lastsyn->next = newsyn;
      lastsyn = newsyn;
    }
    str = ptr;
  }

  MemFree (tmp);
  return TRUE;
}

static void ExpandGeneSynList (GeneRefPtr grp)

{
  ValNodePtr       currsyn;
  ValNodePtr       nextsyn;
  ValNodePtr PNTR  prevsyn;

  if (grp == NULL || grp->syn == NULL) return;

  currsyn = grp->syn;
  prevsyn = &(grp->syn);
  while (currsyn != NULL) {
    if (ExpandGeneSynCom (currsyn)) {
      nextsyn = currsyn->next;
      *(prevsyn) = currsyn->next;
      currsyn->next = NULL;
      ValNodeFreeData (currsyn);
    } else {
      nextsyn = currsyn->next;
      prevsyn = (ValNodePtr PNTR) &(currsyn->next);
    }
    currsyn = nextsyn;
  }

  currsyn = grp->syn;
  prevsyn = &(grp->syn);
  while (currsyn != NULL) {
    if (ExpandGeneSynSem (currsyn)) {
      nextsyn = currsyn->next;
      *(prevsyn) = currsyn->next;
      currsyn->next = NULL;
      ValNodeFreeData (currsyn);
    } else {
      nextsyn = currsyn->next;
      prevsyn = (ValNodePtr PNTR) &(currsyn->next);
    }
    currsyn = nextsyn;
  }
}

typedef struct gosstruc {
  CharPtr       term;
  Char          goid [32];
  CharPtr       evidence;
  Int4          pmid;
  CharPtr       goref;
  UserFieldPtr  ufp;
} GosStruc, PNTR GosStrucPtr;

static int LIBCALLBACK SortVnpByGssp (VoidPtr ptr1, VoidPtr ptr2)

{
  int           compare;
  GosStrucPtr   gsp1, gsp2;
  ValNodePtr    vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  gsp1 = (GosStrucPtr) vnp1->data.ptrvalue;
  gsp2 = (GosStrucPtr) vnp2->data.ptrvalue;
  if (gsp1 == NULL || gsp2 == NULL) return 0;

  compare = StringICmp (gsp1->goid, gsp2->goid);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  compare = StringICmp (gsp1->term, gsp2->term);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  compare = StringICmp (gsp1->evidence, gsp2->evidence);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  if (gsp1->pmid == 0) return 1;
  if (gsp2->pmid == 0) return -1;
  if (gsp1->pmid > gsp2->pmid) {
    return 1;
  } else if (gsp1->pmid < gsp2->pmid) {
    return -1;
  }

  return 0;
}

static CharPtr bsecGoQualType [] = {
  "", "Process", "Component", "Function", NULL
};

static CharPtr bsecGoFieldType [] = {
  "", "text string", "go id", "pubmed id", "go ref", "evidence", NULL
};

static UserFieldPtr SortGoTerms (
  UserFieldPtr entryhead
)

{
  UserFieldPtr  entry, topufp, ufp, lastufp;
  CharPtr       evidence, goid, goref, textstr;
  Char          gid [32];
  GosStrucPtr   gsp, lastgsp;
  ValNodePtr    head = NULL, vnp;
  Int2          j;
  ObjectIdPtr   oip;
  Int4          pmid;

  if (entryhead == NULL) return entryhead;

  for (entry = entryhead; entry != NULL; entry = entry->next) {
    if (entry == NULL || entry->choice != 11) break;
    topufp = (UserFieldPtr)  entry->data.ptrvalue;
    if (topufp == NULL) continue;

    textstr = NULL;
    evidence = NULL;
    goid = NULL;
    goref = NULL;
    pmid = 0;
    for (ufp = topufp; ufp != NULL; ufp = ufp->next) {
      oip = ufp->label;
      if (oip == NULL) continue;
      for (j = 0; bsecGoFieldType [j] != NULL; j++) {
        if (StringICmp (oip->str, bsecGoFieldType [j]) == 0) break;
      }
      if (bsecGoFieldType [j] == NULL) continue;
      switch (j) {
        case 1 :
          if (ufp->choice == 1) {
            textstr = (CharPtr) ufp->data.ptrvalue;
          }
          break;
        case 2 :
          if (ufp->choice == 1) {
            goid = (CharPtr) ufp->data.ptrvalue;
          } else if (ufp->choice == 2) {
            sprintf (gid, "%ld", (long) (Int4) ufp->data.intvalue);
            goid = (CharPtr) gid;
          }
          break;
        case 3 :
          if (ufp->choice == 2) {
            pmid = (Int4) ufp->data.intvalue;
          }
          break;
        case 4 :
          if (ufp->choice == 1) {
            goref = (CharPtr) ufp->data.ptrvalue;
          }
          break;
        case 5 :
          if (ufp->choice == 1) {
            evidence = (CharPtr) ufp->data.ptrvalue;
          }
          break;
        default :
          break;
      }
    }

    if (StringDoesHaveText (textstr)) {
      gsp = (GosStrucPtr) MemNew (sizeof (GosStruc));
      if (gsp != NULL) {
        gsp->term = textstr;
        StringNCpy_0 (gsp->goid, goid, sizeof (gsp->goid));
        gsp->evidence = evidence;
        gsp->pmid = pmid;
        gsp->goref = goref;
        gsp->ufp = entry;
        ValNodeAddPointer (&head, 0, (Pointer) gsp);
      }
    }
  }

  if (head == NULL) return entryhead;
  head = ValNodeSort (head, SortVnpByGssp);

  entryhead = NULL;
  lastgsp = NULL;
  lastufp = NULL;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    gsp = (GosStrucPtr) vnp->data.ptrvalue;
    if (gsp == NULL || gsp->ufp == NULL) continue;
    if (lastgsp != NULL &&
        (StringICmp (gsp->term, lastgsp->term) == 0 || StringICmp (gsp->goid, lastgsp->goid) == 0) &&
         (gsp->pmid == lastgsp->pmid &&
          StringICmp (gsp->goref, lastgsp->goref) == 0 &&
          StringICmp (gsp->evidence, lastgsp->evidence) == 0)) {
      gsp->ufp->next = NULL;
      UserFieldFree (gsp->ufp);
    } else {
      if (lastufp != NULL) {
        lastufp->next = gsp->ufp;
      } else {
        entryhead = gsp->ufp;
      }
      lastufp = gsp->ufp;
      lastufp->next = NULL;
    }
    lastgsp = gsp;
  }

  ValNodeFreeData (head);

  return entryhead;
}

static void SortGoTermsUfp (
  UserFieldPtr ufp,
  Pointer userdata
)

{
  UserFieldPtr  entry;
  Int2          i;
  ObjectIdPtr   oip;
 
  if (ufp == NULL || ufp->choice != 11) return;
  oip = ufp->label;
  if (oip == NULL) return;
  for (i = 0; bsecGoQualType [i] != NULL; i++) {
    if (StringICmp (oip->str, bsecGoQualType [i]) == 0) break;
  }
  if (bsecGoQualType [i] == NULL) return;

  entry = ufp->data.ptrvalue;
  if (entry == NULL || entry->choice != 11) return;

  ufp->data.ptrvalue = SortGoTerms (entry);
}

static void SortGoTermsSfp (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr  oip;

  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "GeneOntology") == 0) {
    VisitUserFieldsInUop (uop, userdata, SortGoTermsUfp);
  }
}

static void CleanupGoTerms (
  UserFieldPtr entryhead
)

{
  UserFieldPtr  entry, topufp, ufp;
  CharPtr       goid, goref, str;
  Int2          j;
  ObjectIdPtr   oip;

  if (entryhead == NULL) return;

  for (entry = entryhead; entry != NULL; entry = entry->next) {
    if (entry == NULL || entry->choice != 11) break;
    topufp = (UserFieldPtr)  entry->data.ptrvalue;
    if (topufp == NULL) continue;

    goid = NULL;
    goref = NULL;
    for (ufp = topufp; ufp != NULL; ufp = ufp->next) {
      oip = ufp->label;
      if (oip == NULL) continue;
      for (j = 0; bsecGoFieldType [j] != NULL; j++) {
        if (StringICmp (oip->str, bsecGoFieldType [j]) == 0) break;
      }
      if (bsecGoFieldType [j] == NULL) continue;
      switch (j) {
        case 2 :
          if (ufp->choice == 1) {
            goid = (CharPtr) ufp->data.ptrvalue;
            if (goid != NULL && *goid != '\0') {
              if (StringNICmp (goid, "GO:", 3) == 0) {
                str = StringSave (goid + 3);
                ufp->data.ptrvalue = (Pointer) str;
                MemFree (goid);
              }
            }
          }
          break;
        case 4 :
          if (ufp->choice == 1) {
            goref = (CharPtr) ufp->data.ptrvalue;
            if (goref != NULL && *goref != '\0') {
              if (StringNICmp (goref, "GO_REF:", 7) == 0) {
                str = StringSave (goref + 7);
                ufp->data.ptrvalue = (Pointer) str;
                MemFree (goref);
              }
            }
          }
          break;
        default :
          break;
      }
    }
  }
}

static void CleanupGoTermsUfp (
  UserFieldPtr ufp,
  Pointer userdata
)

{
  UserFieldPtr  entry;
  Int2          i;
  ObjectIdPtr   oip;
 
  if (ufp == NULL || ufp->choice != 11) return;
  oip = ufp->label;
  if (oip == NULL) return;
  for (i = 0; bsecGoQualType [i] != NULL; i++) {
    if (StringICmp (oip->str, bsecGoQualType [i]) == 0) break;
  }
  if (bsecGoQualType [i] == NULL) return;

  entry = ufp->data.ptrvalue;
  if (entry == NULL || entry->choice != 11) return;

  CleanupGoTerms (entry);
}

static void CleanupGoTermsSfp (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr  oip;

  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "GeneOntology") == 0) {
    VisitUserFieldsInUop (uop, userdata, CleanupGoTermsUfp);
  }
}

static CharPtr CleanUpSgml (
  CharPtr str
)

{
  Int2     ascii_len;
  Char     buf [256];
  CharPtr  ptr;

  if (StringHasNoText (str)) return NULL;
  if (StringChr (str, '&') == NULL) return NULL;

  ascii_len = Sgml2AsciiLen (str);
  if (ascii_len + 2 >= sizeof (buf)) return NULL;

  buf [0] = '\0';
  Sgml2Ascii (str, buf, ascii_len + 1);
  if (StringHasNoText (buf)) return NULL;
  if (StringCmp (str, buf) == 0) return NULL;

  ptr = StringChr (buf, '<');
  if (ptr != NULL) {
    *ptr = ' ';
  }
  ptr = StringChr (buf, '>');
  if (ptr != NULL) {
    *ptr = ' ';
  }
  TrimSpacesAroundString (buf);
  Asn2gnbkCompressSpaces (buf);

  return StringSave (buf);
}

/* special exception for genome pipeline rRNA names */

static Boolean NotExceptedRibosomalName (
  CharPtr name
)

{
  Char     ch;
  CharPtr  str;

  str = StringStr (name, " ribosomal");
  if (str == NULL) return FALSE;

  str += 10;
  ch = *str;
  while (ch != '\0') {
    if (ch == ' ' || IS_DIGIT (ch)) {
      /* okay */
    } else {
      return TRUE;
    }
    str++;
    ch = *str;
  }

  return FALSE;
}

NLM_EXTERN void CleanupSubSourceOrgModOtherFeat (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioSourcePtr  biop;
  OrgNamePtr    onp = NULL;
  OrgRefPtr     orp;

  if (sfp == NULL) return;
  if (sfp->data.choice != SEQFEAT_BIOSRC) return;
  biop = (BioSourcePtr) sfp->data.value.ptrvalue;
  if (biop == NULL) return;
  orp = biop->org;
  if (orp != NULL) {
    onp = orp->orgname;
    if (orp != NULL) {
      CleanupOrgModOther (biop, onp);
    }
  }
  CleanupSubSourceOther (biop, onp);
}

NLM_EXTERN void CleanupSubSourceOrgModOtherDesc (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  BioSourcePtr  biop;
  OrgNamePtr    onp = NULL;
  OrgRefPtr     orp;

  if (sdp == NULL) return;
  if (sdp->choice != Seq_descr_source) return;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;
  orp = biop->org;
  if (orp != NULL) {
    onp = orp->orgname;
    if (orp != NULL) {
      CleanupOrgModOther (biop, onp);
    }
  }
  CleanupSubSourceOther (biop, onp);
}

typedef struct xmltable {
  CharPtr  code;
  size_t   len;
  CharPtr  letter;
} XmlTable, PNTR XmlTablePtr;

static XmlTable xmlunicodes [] = {
  { "&amp",     4, "&"},
  { "&apos",    5, "\'"},
  { "&gt",      3, ">"},
  { "&lt",      3, "<"},
  { "&quot",    5, "\""},
  { "&#13&#10", 8, ""},
  { "&#916",    5, "Delta"},
  { "&#945",    5, "alpha"},
  { "&#946",    5, "beta"},
  { "&#947",    5, "gamma"},
  { "&#952",    5, "theta"},
  { "&#955",    5, "lambda"},
  { "&#956",    5, "mu"},
  { "&#957",    5, "nu"},
  { "&#8201",   6, " "},
  { "&#8206",   6, ""},
  { "&#8242",   6, "'"},
  { "&#8594",   6, "->"},
  { "&#8722",   6, "-"},
  { "&#8710",   6, "delta"},
  { "&#64257",  7, "fi"},
  { "&#64258",  7, "fl"},
  { "&#65292",  7, ","},
  { NULL,       0, ""}
};

static CharPtr BSECDecodeXml (
  CharPtr str
)

{
  Char         ch, nxt;
  CharPtr      dst, ptr, src;
  Int2         i;
  size_t       len;
  XmlTablePtr  xtp;

  if (StringHasNoText (str)) return str;

  src = str;
  dst = str;
  ch = *src;
  while (ch != '\0') {
    if (ch == '&') {
      xtp = NULL;
      len = 1;
      for (i = 0; xmlunicodes [i].code != NULL; i++) {
        if (StringNICmp (src, xmlunicodes [i].code, xmlunicodes [i].len) == 0) {
          nxt = *(src +xmlunicodes [i].len);
          if (nxt == ';') {
            xtp = &(xmlunicodes [i]);
            len = xtp->len + 1;
            break;
          } else if (nxt == ' ' || nxt == '\0') {
            xtp = &(xmlunicodes [i]);
            len = xtp->len;
            break;
          }
        }
      }
      if (xtp != NULL) {
        if (StringLen (xtp->letter) > 0) {
          ptr = xtp->letter;
          ch = *ptr;
          while (ch != '\0') {
            *dst = ch;
            dst++;
            ptr++;
            ch = *ptr;
          }
        }
        src += len;
      } else {
        *dst = ch;
        dst++;
        src++;
      }
    } else {
      *dst = ch;
      dst++;
      src++;
    }
    ch = *src;
  }
  *dst = '\0';

  return str;
}

static void CleanupFeatureStrings (
  SeqFeatPtr sfp,
  Boolean isJscan,
  Boolean stripSerial,
  Boolean modernizeFeats,
  ValNodePtr PNTR publist
)

{
  Uint1         aa;
  BioSourcePtr  biop;
  Char          ch;
  Uint1         codon [6];
  GeneRefPtr    grp;
  ImpFeatPtr    ifp;
  Boolean       is_fMet = FALSE;
  Int2          j;
  Boolean       justTrnaText;
  size_t        len;
  CharPtr       name;
  OrgNamePtr    onp = NULL;
  OrgRefPtr     orp;
  PubdescPtr    pdp;
  ProtRefPtr    prp;
  CharPtr       ptr;
  RnaRefPtr     rrp;
  SubSourcePtr  ssp;
  CharPtr       str;
  CharPtr       suff;
  CharPtr       temp;
  Char          tmp [64];
  Boolean       trimming_junk;
  tRNAPtr       trp;
  CharPtr       val;
  ValNodePtr    vnp;
  RNAGenPtr     rgp;
  RNAQualPtr    rqp;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL) return;
  BSECDecodeXml (sfp->comment);
  CleanVisString (&(sfp->comment));
  len = StringLen (sfp->comment);
  if (len > 4) {
    if (StringCmp (sfp->comment + len - 3, ",..") == 0 ||
        StringCmp (sfp->comment + len - 3, ".,.") == 0 ||
        StringCmp (sfp->comment + len - 3, "..,") == 0 ||
        StringCmp (sfp->comment + len - 3, ",.,") == 0) {
      sfp->comment [len - 3] = '.';
      sfp->comment [len - 2] = '.';
      sfp->comment [len - 1] = '.';
    }
  }
  BSECDecodeXml (sfp->title);
  CleanVisString (&(sfp->title));
  CleanVisString (&(sfp->except_text));
  if (StringDoesHaveText (sfp->except_text)) {
    CleanUpExceptText (sfp);
  }
  CleanDoubleQuote (sfp->comment);
  if (StringCmp (sfp->comment, ".") == 0) {
    sfp->comment = MemFree (sfp->comment);
  }
  /*
  if (sfp->ext != NULL) {
    VisitUserObjectsInUop (sfp->ext, NULL, SortGoTermsSfp);
  }
  */
  if (sfp->ext != NULL) {
    VisitUserObjectsInUop (sfp->ext, NULL, CleanupGoTermsSfp);
  }

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->data.choice != SEQFEAT_PROT) continue;
    prp = (ProtRefPtr) xref->data.value.ptrvalue;
    if (prp == NULL) continue;
    RemoveFlankingQuotes (&(prp->desc));
    RemoveFlankingQuotesList (&(prp->name));
  }

  switch (sfp->data.choice) {
    case SEQFEAT_BOND :
    case SEQFEAT_PSEC_STR :
    case SEQFEAT_COMMENT:
      return;
    case SEQFEAT_SITE :
      for (j = 0; bsecSiteList [j] != NULL; j++) {
        StringNCpy_0 (tmp, bsecSiteList [j], sizeof (tmp));
        len = StringLen (tmp);
        if (StringNICmp (sfp->comment, tmp, len) == 0) {
          if (sfp->data.value.intvalue == 0 || sfp->data.value.intvalue == 255) {
            sfp->data.value.intvalue = j;
            if (StringHasNoText (sfp->comment + len) || StringICmp (sfp->comment + len, " site") == 0) {
              sfp->comment = MemFree (sfp->comment);
            }
          }
        } else {
          val = tmp;
          ch = *val;
          while (ch != '\0') {
            if (ch == '-') {
              *val = ' ';
            }
            val++;
            ch = *val;
          }
          if (StringNICmp (sfp->comment, tmp, len) == 0) {
            if (sfp->data.value.intvalue == 0 || sfp->data.value.intvalue == 255) {
              sfp->data.value.intvalue = j;
              if (StringHasNoText (sfp->comment + len) || StringICmp (sfp->comment + len, " site") == 0) {
                sfp->comment = MemFree (sfp->comment);
              }
            }
          }
        }
      }
      break;
    default :
      break;
  }
  if (sfp->data.value.ptrvalue == NULL) return;

  biop = NULL;
  orp = NULL;
  switch (sfp->data.choice) {
    case SEQFEAT_ORG :
      orp = (OrgRefPtr) sfp->data.value.ptrvalue;
      break;
    case SEQFEAT_BIOSRC :
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
      }
    default :
      break;
  }
  if (orp != NULL && sfp->qual != NULL) {
    GbqualToOrpMod (&(sfp->qual), &(orp->mod));
  }

  biop = NULL;
  orp = NULL;
  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (sfp->xref != NULL) {
        HandleXrefOnGene (sfp);
      }
      BSECDecodeXml (grp->locus);
      CleanVisString (&(grp->locus));
      /*
      if (isJscan && StringDoesHaveText (grp->locus)) {
        ptr = CleanUpSgml (grp->locus);
        if (ptr != NULL) {
          grp->locus = MemFree (grp->locus);
          grp->locus = StringSave (ptr);
        }
      }
      */
      CleanVisString (&(grp->allele));
      CleanVisString (&(grp->desc));
      CleanVisString (&(grp->maploc));
      CleanVisString (&(grp->locus_tag));
      ExpandGeneSynList (grp);
      /*
      if (isJscan && grp->syn != NULL) {
        for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          ptr = CleanUpSgml (str);
          if (ptr != NULL) {
            vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
            vnp->data.ptrvalue = StringSave (ptr);
          }
        }
      }
      */
      for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        BSECDecodeXml (str);
      }
      CleanVisStringListCaseSensitive (&(grp->syn));
      grp->syn = ValNodeSort (grp->syn, SortVnpByStringCS);
      grp->syn = UniqueStringValNodeCS (grp->syn);
      grp->syn = ValNodeSort (grp->syn, SortVnpByStringCILCFirst);
      CleanDoubleQuote (grp->locus);
      CleanDoubleQuote (grp->allele);
      CleanDoubleQuote (grp->desc);
      /*
      if (isJscan && StringDoesHaveText (grp->desc)) {
        ptr = CleanUpSgml (grp->desc);
        if (ptr != NULL) {
          grp->desc = MemFree (grp->desc);
          grp->desc = StringSave (ptr);
        }
      }
      */
      CleanDoubleQuote (grp->maploc);
      CleanDoubleQuote (grp->locus_tag);
      CleanDoubleQuoteList (grp->syn);
      FixOldDbxrefs (grp->db);
      FixNumericDbxrefs (grp->db);
      grp->db = ValNodeSort (grp->db, SortDbxref);
      CleanupDuplicateDbxrefs (&(grp->db));
      CleanupObsoleteDbxrefs (&(grp->db));
      /* now move grp->dbxref to sfp->dbxref */
      vnp = grp->db;
      grp->db = NULL;
      ValNodeLink ((&sfp->dbxref), vnp);
      if (grp->locus != NULL && grp->syn != NULL) {
        vnp = grp->syn;
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringCmp (grp->locus, str) == 0) {
          grp->syn = vnp->next;
          vnp->next = NULL;
          ValNodeFreeData (vnp);
        }
      }
      /*
      if (grp->locus != NULL && sfp->comment != NULL) {
        if (StringCmp (grp->locus, sfp->comment) == 0) {
          sfp->comment = MemFree (sfp->comment);
        }
      }
      */
      break;
    case SEQFEAT_ORG :
      orp = (OrgRefPtr) sfp->data.value.ptrvalue;
      break;
    case SEQFEAT_CDREGION :
      if (sfp->xref != NULL && sfp->product != NULL) {
        HandleXrefOnCDS (sfp);
      }
      break;
    case SEQFEAT_PROT :
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        CleanupECNumber (str);
      }
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        BSECDecodeXml (str);
      }
      BSECDecodeXml (prp->desc);
      CleanVisStringAndCompress (&(prp->desc));
      CleanVisStringListAndCompress (&(prp->name));
      CleanVisStringList (&(prp->ec));
      CleanVisStringList (&(prp->activity));
      CleanDoubleQuote (prp->desc);
      CleanDoubleQuoteList (prp->name);
      CleanDoubleQuoteList (prp->ec);
      CleanDoubleQuoteList (prp->activity);
      RemoveFlankingQuotes (&(prp->desc));
      RemoveFlankingQuotesList (&(prp->name));
      FixOldDbxrefs (prp->db);
      FixNumericDbxrefs (prp->db);
      prp->db = ValNodeSort (prp->db, SortDbxref);
      CleanupDuplicateDbxrefs (&(prp->db));
      CleanupObsoleteDbxrefs (&(prp->db));
      /* now move prp->dbxref to sfp->dbxref */
      vnp = prp->db;
      prp->db = NULL;
      ValNodeLink ((&sfp->dbxref), vnp);
      if (prp->processed != 3 && prp->processed != 4 &&
          prp->name == NULL && sfp->comment != NULL) {
        if (StringICmp (sfp->comment, "putative") != 0) {
          ValNodeAddStr (&(prp->name), 0, sfp->comment);
          sfp->comment = NULL;
        }
      }
      if (prp->processed == 3 || prp->processed == 4) {
        if (prp->name != NULL) {
          str = (CharPtr) prp->name->data.ptrvalue;
          if ((StringStr (str, "putative") != NULL ||
               StringStr (str, "put. ") != NULL) &&
              sfp->comment == NULL) {
            sfp->comment = StringSave ("putative");
          }
          if (! InformativeString (str)) {
            prp->name = ValNodeFreeData (prp->name);
          }
        }
      }
      if ((prp->processed == 1 || prp->processed == 2) && prp->name == NULL) {
        ValNodeCopyStr (&(prp->name), 0, "unnamed");
      }
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringICmp (str, "RbcL") == 0 || StringICmp (str, "rubisco large subunit") == 0) {
          vnp->data.ptrvalue = StringSave ("ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit");
          MemFree (str);
        } else if (StringICmp (str, "RbcS") == 0 || StringICmp (str, "rubisco small subunit") == 0) {
          vnp->data.ptrvalue = StringSave ("ribulose-1,5-bisphosphate carboxylase/oxygenase small subunit");
          MemFree (str);
        }
      }
      if (StringDoesHaveText (prp->desc)) {
        for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          if (StringCmp (prp->desc, str) == 0) {
            prp->desc = MemFree (prp->desc);
          }
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp->ext.choice == 1) {
        BSECDecodeXml ((CharPtr) rrp->ext.value.ptrvalue);
        CleanVisStringAndCompress ((CharPtr PNTR) &(rrp->ext.value.ptrvalue));
        CleanDoubleQuote ((CharPtr) rrp->ext.value.ptrvalue);
        RemoveFlankingQuotes ((CharPtr PNTR) &(rrp->ext.value.ptrvalue));
        if (rrp->ext.value.ptrvalue == NULL) {
          rrp->ext.choice = 0;
        } else if (rrp->type == 4) {
          name = (CharPtr) rrp->ext.value.ptrvalue;
          len = StringLen (name);
          if (len > 5) {
            if (len > 14 && StringNICmp (name + len - 14, " ribosomal rRNA", 14) == 0) {
            } else if (StringNICmp (name + len - 5, " rRNA", 5) == 0) {
              str = MemNew (len + 10);
              if (str != NULL) {
                StringNCpy (str, name, len - 5);
                StringCat (str, " ribosomal RNA");
                rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
                rrp->ext.value.ptrvalue = (Pointer) str;
              }
            } else if (StringNICmp (name + len - 5, "_rRNA", 5) == 0) {
              str = MemNew (len + 10);
              if (str != NULL) {
                StringNCpy (str, name, len - 5);
                StringCat (str, " ribosomal RNA");
                rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
                rrp->ext.value.ptrvalue = (Pointer) str;
              }
            }
          }
        } else if (rrp->type == 3) {
          name = (CharPtr) rrp->ext.value.ptrvalue;
          aa = ParseTRnaString (name, &justTrnaText, codon, FALSE);
          if (aa != 0) {
            is_fMet = (Boolean) (StringStr (name, "fMet") != NULL);
            rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
            trp = (tRNAPtr) MemNew (sizeof (tRNA));
            if (trp != NULL) {
              trp->aatype = 2;
              for (j = 0; j < 6; j++) {
                trp->codon [j] = 255;
              }
              if (justTrnaText) {
                for (j = 0; j < 6; j++) {
                  trp->codon [j] = codon [j];
                }
              }
              trp->aa = aa;
              rrp->ext.choice = 2;
              rrp->ext.value.ptrvalue = (Pointer) trp;
              CleanupTrna (sfp, trp);
            }
            if (is_fMet) {
              if (sfp->comment == NULL) {
                sfp->comment = StringSave ("fMet");
              } else {
                len = StringLen (sfp->comment) + StringLen ("fMet") + 5;
                str = MemNew (sizeof (Char) * len);
                StringCpy (str, sfp->comment);
                StringCat (str, "; ");
                StringCat (str, "fMet");
                sfp->comment = MemFree (sfp->comment);
                sfp->comment = str;
              }
            }
          }
        }
      } else if (rrp->ext.choice == 2) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        CleanupTrna (sfp, trp);
      } else if (rrp->type == 3 && (! StringHasNoText (sfp->comment))) {
        aa = ParseTRnaString (sfp->comment, &justTrnaText, codon, TRUE);
        if (aa != 0) {
          trp = (tRNAPtr) MemNew (sizeof (tRNA));
          if (trp != NULL) {
            trp->aatype = 2;
            for (j = 0; j < 6; j++) {
              trp->codon [j] = 255;
            }
            if (justTrnaText) {
              for (j = 0; j < 6; j++) {
                trp->codon [j] = codon [j];
              }
            }
            trp->aa = aa;
            rrp->ext.choice = 2;
            rrp->ext.value.ptrvalue = (Pointer) trp;
            if (justTrnaText) {
              if (StringCmp (sfp->comment, "fMet") != 0 &&
                  StringCmp (sfp->comment, "fMet tRNA") != 0 &&
                  StringCmp (sfp->comment, "fMet-tRNA") != 0) {
                sfp->comment = MemFree (sfp->comment);
              } else {
                sfp->comment = MemFree (sfp->comment);
                sfp->comment = StringSave ("fMet");
              }
            }
          }
        }
      }
      if (rrp->ext.choice == 3) {
        rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
        if (rgp != NULL) {
          CleanVisStringAndCompress (&(rgp->product));
          CleanDoubleQuote (rgp->product);
          RemoveFlankingQuotes (&(rgp->product));
          CleanVisStringAndCompress (&(rgp->_class));
          CleanDoubleQuote (rgp->_class);
          for (rqp = rgp->quals; rqp != NULL; rqp = rqp->next) {
            CleanVisStringAndCompress (&(rqp->qual));
            CleanDoubleQuote (rqp->qual);
            CleanVisStringAndCompress (&(rqp->val));
            CleanDoubleQuote (rqp->val);
          }
        }
      }
      if (rrp->ext.choice == 0 && sfp->comment != NULL && rrp->type == 4) {
        len = StringLen (sfp->comment);
        if (len > 15 && len < 20) {
          if (StringNICmp (sfp->comment + len - 15, "S ribosomal RNA", 15) == 0) {
            rrp->ext.choice = 1;
            rrp->ext.value.ptrvalue = sfp->comment;
            sfp->comment = NULL;
          }
        } else if (len > 6 && len < 20) {
          if (StringNICmp (sfp->comment + len - 6, "S rRNA", 6) == 0) {
            rrp->ext.choice = 1;
            rrp->ext.value.ptrvalue = sfp->comment;
            sfp->comment = NULL;
          }
        }
      }
/*
 * This section has been commented out based on a request by DeAnne Cravaritis.
 * If left in, this causes unexpected results when RNA comments are copied to 
 * the product name or vice versa.      
      if (rrp->ext.choice == 1 && rrp->ext.value.ptrvalue != NULL) {
        if (StringICmp ((CharPtr) rrp->ext.value.ptrvalue, sfp->comment) == 0) {
          sfp->comment = MemFree (sfp->comment);
        }
      }
*/      
      if (rrp->type == 4 && rrp->ext.choice == 1 ) {
        name = (CharPtr) rrp->ext.value.ptrvalue;
        len = StringLen (name);
        if (len > 5 && NotExceptedRibosomalName (name)) {
          suff = NULL;
          str = StringStr (name, " ribosomal");
          if (str != NULL) {
            suff = str + 10;
            ch = *suff;
            if (ch != '\0' && ch != ' ') {
              suff = NULL;
              str = NULL;
            }
          }
          if (str == NULL) {
            str = StringStr (name, " rRNA");
            if (str != NULL) {
              suff = str + 5;
              ch = *suff;
              if (ch != '\0' && ch != ' ') {
                suff = NULL;
                str = NULL;
              }
            }
          }
          if (suff != NULL && StringNICmp (suff, " RNA", 4) == 0) {
            suff += 4;
          }
          if (suff != NULL && StringNICmp (suff, " DNA", 4) == 0) {
            suff += 4;
          }
          if (suff != NULL && StringNICmp (suff, " ribosomal", 10) == 0) {
            suff += 10;
          }
          TrimSpacesAroundString (suff);
          if (str != NULL) {
            *str = '\0';
            len = StringLen (name);
            if (StringHasNoText (suff)) {
              suff = NULL;
            }
            if (suff != NULL) {
              len += StringLen (suff) + 2;
            }
            str = MemNew (len + 15);
            if (str != NULL) {
              StringCpy (str, name);
              StringCat (str, " ribosomal RNA");
              if (suff != NULL) {
                ch = *suff;
                if (ch != ',' && ch != ';') {
                  StringCat (str, " ");
                }
                StringCat (str, suff);
              }
              rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
              rrp->ext.value.ptrvalue = (Pointer) str;
            }
          }
        }
        name = (CharPtr) rrp->ext.value.ptrvalue;
        len = StringLen (name);
        if (len > 5) {
          ch = *name;
          while (ch != '\0' && (ch == '.' || (IS_DIGIT (ch)))) {
            name++;
            ch = *name;
          }
          /*
          if (ch == 's' && StringCmp (name, "s ribosomal RNA") == 0) {
            *name = 'S';
          }
          */
          if (ch == 's' && name [1] == ' ') {
            *name = 'S';
          }
        }
        StrStripSpaces ((CharPtr) rrp->ext.value.ptrvalue);
        name = (CharPtr) rrp->ext.value.ptrvalue;
        len = StringLen (name);
        if (len > 17) {
          if (StringNICmp (name + len - 17, "ribosomal RNA RNA", 17) == 0) {
            *(name + len - 4) = '\0';
          }
        }
        trimming_junk = TRUE;
        while (trimming_junk) {
          StrStripSpaces ((CharPtr) rrp->ext.value.ptrvalue);
          name = (CharPtr) rrp->ext.value.ptrvalue;
          ptr = StringStr (name, "ribosomal ribosomal");
          if (ptr != NULL) {
            suff = ptr + 19;
            *(ptr + 10) = '\0';
            temp = MemNew (StringLen (name) + StringLen (suff) + 2);
            TrimSpacesAroundString (suff);
            StringCpy (temp, name);
            if (suff [0] != ' ' && suff [0] != '\0') {
              StringCat (temp, " ");
            }
            StringCat (temp, suff);
            rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
            rrp->ext.value.ptrvalue = (Pointer) temp;
          } else {
            ptr = StringStr (name, "RNA RNA");
            if (ptr != NULL) {
              suff = ptr + 7;
              *(ptr + 4) = '\0';
              temp = MemNew (StringLen (name) + StringLen (suff) + 2);
              TrimSpacesAroundString (suff);
              StringCpy (temp, name);
              if (suff [0] != ' ' && suff [0] != '\0') {
                StringCat (temp, " ");
              }
              StringCat (temp, suff);
              rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
              rrp->ext.value.ptrvalue = (Pointer) temp;
            } else {
              ptr = StringStr (name, "ribosomal RNA ribosomal");
              if (ptr != NULL) {
                suff = ptr + 23;
                *(ptr + 14) = '\0';
                temp = MemNew (StringLen (name) + StringLen (suff) + 2);
                TrimSpacesAroundString (suff);
                StringCpy (temp, name);
                if (suff [0] != ' ' && suff [0] != '\0') {
                  StringCat (temp, " ");
                }
                StringCat (temp, suff);
                rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
                rrp->ext.value.ptrvalue = (Pointer) temp;
              } else {
                ptr = StringStr (name, "ribosomal rRNA");
                if (ptr != NULL) {
                  suff = ptr + 14;
                  *(ptr + 10) = '\0';
                  temp = MemNew (StringLen (name) + StringLen (" RNA") + StringLen (suff) + 2);
                  TrimSpacesAroundString (suff);
                  StringCpy (temp, name);
                  StringCat (temp, " RNA");
                  if (suff [0] != ' ' && suff [0] != '\0') {
                    StringCat (temp, " ");
                  }
                  StringCat (temp, suff);
                  rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
                  rrp->ext.value.ptrvalue = (Pointer) temp;
                } else {
                  ptr = StringStr (name, "RNA rRNA");
                  if (ptr != NULL) {
                    suff = ptr + 8;
                    *(ptr + 3) = '\0';
                    temp = MemNew (StringLen (name) + StringLen (suff) + 2);
                    TrimSpacesAroundString (suff);
                    StringCpy (temp, name);
                    if (suff [0] != ' ' && suff [0] != '\0') {
                      StringCat (temp, " ");
                    }
                    StringCat (temp, suff);
                    rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
                    rrp->ext.value.ptrvalue = (Pointer) temp;
                  } else {
                    trimming_junk = FALSE;
                  }
                }
              }
            }
          }
        }
        TrimSpacesAroundString ((CharPtr) rrp->ext.value.ptrvalue);
      }
      /*
      if (rrp->type == 2 && rrp->ext.choice == 0 && sfp->comment != NULL) {
        rrp->ext.choice = 1;
        rrp->ext.value.ptrvalue = sfp->comment;
        sfp->comment = NULL;
      }
      */
      if (rrp->type == 2 && rrp->ext.choice == 0 && sfp->comment != NULL) {
        len = StringLen (sfp->comment);
        if (len > 5) {
          if (StringNICmp (sfp->comment + len - 4, " RNA", 4) == 0 ||
              StringNICmp (sfp->comment + len - 5, " mRNA", 5) == 0) {
            rrp->ext.choice = 1;
            rrp->ext.value.ptrvalue = sfp->comment;
            sfp->comment = NULL;
          }
        }
      }
      if (rrp->type == 255 || rrp->type == 10) {
        name = GetRNARefProductString (rrp, NULL);
        if (StringICmp (name, "its1") == 0 || StringICmp (name, "its 1") == 0) {
          SetRNARefProductString (rrp, NULL, "internal transcribed spacer 1", ExistingTextOption_replace_old);
        } else if (StringICmp (name, "its2") == 0 || StringICmp (name, "its 2") == 0) {
          SetRNARefProductString (rrp, NULL, "internal transcribed spacer 2", ExistingTextOption_replace_old);
        } else if (StringICmp (name, "its3") == 0 || StringICmp (name, "its 3") == 0) {
          SetRNARefProductString (rrp, NULL, "internal transcribed spacer 3", ExistingTextOption_replace_old);
        }
        name = MemFree (name);
      }
      if ((rrp->type == 255 || rrp->type == 10) && rrp->ext.choice == 0 && sfp->comment != NULL) {
        if (StringICmp (sfp->comment, "internal transcribed spacer 1") == 0 ||
            StringICmp (sfp->comment, "internal transcribed spacer 2") == 0 ||
            StringICmp (sfp->comment, "internal transcribed spacer 3") == 0) {
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = sfp->comment;
          sfp->comment = NULL;
        }
      }
      break;
    case SEQFEAT_PUB :
      pdp = (PubdescPtr) sfp->data.value.ptrvalue;
      CleanDoubleQuote (pdp->comment);
      NormalizePubdesc (pdp, stripSerial, TRUE, publist);
      break;
    case SEQFEAT_SEQ :
      break;
    case SEQFEAT_IMP :
      ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
      CleanVisString (&(ifp->key));
      CleanVisString (&(ifp->loc));
      CleanVisString (&(ifp->descr));
      break;
    case SEQFEAT_REGION :
      CleanVisStringAndCompress ((CharPtr PNTR) &(sfp->data.value.ptrvalue));
      CleanDoubleQuote ((CharPtr) sfp->data.value.ptrvalue);
      if (sfp->data.value.ptrvalue == NULL) {
        sfp->data.choice = SEQFEAT_COMMENT;
      }
      break;
    case SEQFEAT_COMMENT :
      break;
    case SEQFEAT_BOND :
      break;
    case SEQFEAT_SITE :
      break;
    case SEQFEAT_RSITE :
      break;
    case SEQFEAT_USER :
      VisitAllUserObjectsInUop ((UserObjectPtr) sfp->data.value.ptrvalue, NULL, CleanUserObject);
      break;
    case SEQFEAT_TXINIT :
      break;
    case SEQFEAT_NUM :
      break;
    case SEQFEAT_PSEC_STR :
      break;
    case SEQFEAT_NON_STD_RESIDUE :
      break;
    case SEQFEAT_HET :
      break;
    case SEQFEAT_BIOSRC :
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      if (biop != NULL) {
        if (biop->genome == GENOME_virion) {
          biop->genome = GENOME_unknown;
        }
        orp = biop->org;
        if (orp != NULL) {
          CleanVisStringListAndCompress (&(orp->mod));
          OrpModToSubSource (&(orp->mod), &(biop->subtype));
          onp = orp->orgname;
          if (onp != NULL) {
            CleanupOrgModOther (biop, onp);
          }
        }
        biop->subtype = SortSubSourceList (biop->subtype);
        CleanSubSourceList (&(biop->subtype), biop->genome);
        CleanupSubSourceOther (biop, onp);
        biop->subtype = SortSubSourceList (biop->subtype);
        if (modernizeFeats) {
          ModernizePCRPrimers (biop);
        }
        CleanupPCRReactionSet (&(biop->pcr_primers));
        if (biop->genome == GENOME_unknown || biop->genome == GENOME_genomic) {
          for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
            if (ssp->subtype == SUBSRC_plasmid_name) {
              biop->genome = GENOME_plasmid;
            }
          }
        }
      }
      break;
    default :
      break;
  }
  if (orp != NULL) {
    CleanVisStringAndCompress (&(orp->taxname));
    CleanVisStringAndCompress (&(orp->common));
    CleanVisStringList (&(orp->mod));
    CleanVisStringList (&(orp->syn));
    FixOldDbxrefs (orp->db);
    FixNumericDbxrefs (orp->db);
    orp->db = ValNodeSort (orp->db, SortDbxref);
    CleanupDuplicateDbxrefs (&(orp->db));
    CleanupObsoleteDbxrefs (&(orp->db));
    onp = orp->orgname;
    while (onp != NULL) {
      CleanVisString (&(onp->attrib));
      CleanVisString (&(onp->lineage));
      CleanVisString (&(onp->div));
      OrpModToOrgMod (&(orp->mod), &(onp->mod));
      onp->mod = SortOrgModList (onp->mod);
      CleanOrgModListEx (&(onp->mod), orp->common);
      onp->mod = SortOrgModList (onp->mod);
      onp = onp->next;
    }
  }
}

static void CleanupDescriptorStrings (
  ValNodePtr sdp,
  Boolean stripSerial,
  Boolean modernizeFeats,
  ValNodePtr PNTR publist,
  Boolean isEmblOrDdbj
)

{
  BioSourcePtr  biop;
  EMBLBlockPtr  ebp;
  GBBlockPtr    gbp;
  OrgNamePtr    onp = NULL;
  OrgRefPtr     orp;
  PubdescPtr    pdp;
  SubSourcePtr  ssp;

  if (sdp == NULL) return;
  switch (sdp->choice) {
    case Seq_descr_mol_type :
    case Seq_descr_method :
      return;
    default :
      break;
  }
  if (sdp->data.ptrvalue == NULL) return;

  biop = NULL;
  orp = NULL;
  switch (sdp->choice) {
    case Seq_descr_mol_type :
      break;
    case Seq_descr_modif :
      break;
    case Seq_descr_method :
      break;
    case Seq_descr_name :
      CleanVisString ((CharPtr PNTR) &sdp->data.ptrvalue);
      if (sdp->data.ptrvalue == NULL) {
        sdp->data.ptrvalue = StringSave ("");
      }
      break;
    case Seq_descr_title :
      BSECDecodeXml ((CharPtr) sdp->data.ptrvalue);
      CleanVisStringAndCompress ((CharPtr PNTR) &sdp->data.ptrvalue);
      if (sdp->data.ptrvalue == NULL) {
        sdp->data.ptrvalue = StringSave ("");
      }
      break;
    case Seq_descr_org :
      orp = (OrgRefPtr) sdp->data.ptrvalue;
      break;
    case Seq_descr_comment :
      BSECDecodeXml ((CharPtr) sdp->data.ptrvalue);
      CleanVisStringJunk ((CharPtr PNTR) &sdp->data.ptrvalue);
      RemoveSpacesBetweenTildes ((CharPtr) sdp->data.ptrvalue);
      if (sdp->data.ptrvalue == NULL) {
        sdp->data.ptrvalue = StringSave ("");
      }
      break;
    case Seq_descr_num :
      break;
    case Seq_descr_maploc :
      break;
    case Seq_descr_pir :
      break;
    case Seq_descr_genbank :
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      CleanVisStringList (&(gbp->extra_accessions));
      gbp->extra_accessions = ValNodeSort (gbp->extra_accessions, SortVnpByString);
      gbp->extra_accessions = UniqueValNode (gbp->extra_accessions);
      if (isEmblOrDdbj) {
        CleanVisStringListCaseSensitive (&(gbp->keywords));
      } else {
        CleanVisStringList (&(gbp->keywords));
      }
      CleanVisStringJunk (&(gbp->source));
      if (StringCmp (gbp->source, ".") == 0) {
        gbp->source = MemFree (gbp->source);
      }
      CleanVisStringJunk (&(gbp->origin));
      if (StringCmp (gbp->origin, ".") == 0) {
        gbp->origin = MemFree (gbp->origin);
      }
      CleanVisString (&(gbp->date));
      CleanVisString (&(gbp->div));
      CleanVisString (&(gbp->taxonomy));
      break;
    case Seq_descr_pub :
      pdp = (PubdescPtr) sdp->data.ptrvalue;
      CleanDoubleQuote (pdp->comment);
      NormalizePubdesc (pdp, stripSerial, TRUE, publist);
      break;
    case Seq_descr_region :
      CleanVisString ((CharPtr PNTR) &sdp->data.ptrvalue);
      if (sdp->data.ptrvalue == NULL) {
        sdp->data.ptrvalue = StringSave ("");
      }
      break;
    case Seq_descr_user :
      VisitAllUserObjectsInUop ((UserObjectPtr) sdp->data.ptrvalue, NULL, CleanUserObject);
      break;
    case Seq_descr_sp :
      break;
    case Seq_descr_dbxref :
      break;
    case Seq_descr_embl :
      ebp = (EMBLBlockPtr) sdp->data.ptrvalue;
      CleanVisStringList (&(ebp->extra_acc));
      ebp->extra_acc = ValNodeSort (ebp->extra_acc, SortVnpByString);
      CleanVisStringListCaseSensitive (&(ebp->keywords));
      break;
    case Seq_descr_create_date :
      break;
    case Seq_descr_update_date :
      break;
    case Seq_descr_prf :
      break;
    case Seq_descr_pdb :
      break;
    case Seq_descr_het :
      break;
    case Seq_descr_source :
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop != NULL) {
        if (biop->genome == GENOME_virion) {
          biop->genome = GENOME_unknown;
        }
        orp = biop->org;
        if (orp != NULL) {
          CleanVisStringList (&(orp->mod));
          OrpModToSubSource (&(orp->mod), &(biop->subtype));
          onp = orp->orgname;
          if (onp != NULL) {
            CleanupOrgModOther (biop, onp);
          }
        }
        biop->subtype = SortSubSourceList (biop->subtype);
        CleanSubSourceList (&(biop->subtype), biop->genome);
        CleanupSubSourceOther (biop, onp);
        biop->subtype = SortSubSourceList (biop->subtype);
        if (modernizeFeats) {
          ModernizePCRPrimers (biop);
        }
        CleanupPCRReactionSet (&(biop->pcr_primers));
        if (biop->genome == GENOME_unknown || biop->genome == GENOME_genomic) {
          for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
            if (ssp->subtype == SUBSRC_plasmid_name) {
              biop->genome = GENOME_plasmid;
            }
          }
        }
      }
      break;
    case Seq_descr_molinfo :
      break;
    default :
      break;
  }
  if (orp != NULL) {
    CleanVisStringAndCompress (&(orp->taxname));
    CleanVisStringAndCompress (&(orp->common));
    CleanVisStringList (&(orp->mod));
    CleanVisStringList (&(orp->syn));
    FixOldDbxrefs (orp->db);
    FixNumericDbxrefs (orp->db);
    orp->db = ValNodeSort (orp->db, SortDbxref);
    CleanupDuplicateDbxrefs (&(orp->db));
    CleanupObsoleteDbxrefs (&(orp->db));
    onp = orp->orgname;
    while (onp != NULL) {
      CleanVisString (&(onp->attrib));
      CleanVisString (&(onp->lineage));
      CleanVisString (&(onp->div));
      OrpModToOrgMod (&(orp->mod), &(onp->mod));
      onp->mod = SortOrgModList (onp->mod);
      CleanOrgModListEx (&(onp->mod), orp->common);
      onp->mod = SortOrgModList (onp->mod);
      onp = onp->next;
    }
  }
}

static Int2 CheckForQual (GBQualPtr gbqual, CharPtr string_q, CharPtr string_v)

{
  GBQualPtr curq;

  for (curq = gbqual; curq; curq = curq->next) {
    if (StringCmp (string_q, curq->qual) == 0) {
      if (curq->val == NULL) {
        curq->val = StringSave (string_v);
        return 1;
      } 
      if (StringCmp (string_v, curq->val) == 0) return 1;
    }
  }
  return 0;
}
static GBQualPtr AddGBQual (GBQualPtr gbqual, CharPtr qual, CharPtr val)

{
  GBQualPtr curq;

  if (StringCmp (qual, "translation") == 0) {
    if (val == NULL)  return gbqual;
    if (*val == '\0') return gbqual;
  }
  if (gbqual) {
    if (CheckForQual (gbqual, qual, val) == 1) return gbqual;
    for (curq = gbqual; curq->next != NULL; curq = curq->next) continue;
    curq->next = GBQualNew ();
    curq = curq->next;
    if (val)
      curq->val = StringSave (val);
    curq->qual = StringSave (qual);
  } else {
    gbqual = GBQualNew ();
    gbqual->next = NULL;
    if (val)
      gbqual->val = StringSave (val);
    gbqual->qual = StringSave (qual);
  }
  return gbqual;
}

static void AddReplaceQual (SeqFeatPtr sfp, CharPtr p)

{
  CharPtr s, val;

  val = StringChr (p, '\"');
  if (val == NULL) return;
  val++;
  s = p + StringLen (p) - 1;
  if (*s != ')') return;
  for (s--; s > val && *s != '\"'; s--) continue;
  if (*s != '\"') return;
  *s = '\0';
  sfp->qual = (GBQualPtr) AddGBQual (sfp->qual, "replace", val);
  *s = '\"';
}

NLM_EXTERN Boolean SerialNumberInString (CharPtr str)

{
  Char     ch;
  Boolean  hasdigits;
  CharPtr  ptr;
  Boolean  suspicious = FALSE;

  if (str == NULL || StringHasNoText (str)) return FALSE;
  ptr = StringChr (str, '[');
  while ((! suspicious) && ptr != NULL) {
    hasdigits = FALSE;
    ptr++;
    ch = *ptr;
    while (IS_DIGIT (ch)) {
      hasdigits = TRUE;
      ptr++;
      ch = *ptr;
    }
    if (ch == ']' && hasdigits) {
      suspicious = TRUE;
    }
    if (! suspicious) {
      ptr = StringChr (ptr, '[');
    }
  }
  return suspicious;
}

/* now only strips serials for local, general, refseq, and 2+6 genbank ids */
static void CheckForSwissProtID (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  SeqIdPtr      sip;
  BoolPtr       stripSerial;
  TextSeqIdPtr  tsip;

  if (sep == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    stripSerial = (BoolPtr) mydata;
    if (stripSerial == NULL) return;
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      switch (sip->choice) {
        case SEQID_GIBBSQ :
        case SEQID_GIBBMT :
          *stripSerial = FALSE;
          break;
        case SEQID_EMBL :
        case SEQID_PIR :
        case SEQID_SWISSPROT :
        case SEQID_PATENT :
        case SEQID_DDBJ :
        case SEQID_PRF :
        case SEQID_PDB :
        case SEQID_TPE:
        case SEQID_TPD:
        case SEQID_GPIPE:
          *stripSerial = FALSE;
          break;
        case SEQID_GENBANK :
        case SEQID_TPG:
          tsip = (TextSeqIdPtr) sip->data.ptrvalue;
          if (tsip != NULL) {
            if (StringLen (tsip->accession) == 6) {
              *stripSerial = FALSE;
            }
          }
          break;
        case SEQID_NOT_SET :
        case SEQID_LOCAL :
        case SEQID_OTHER :
        case SEQID_GENERAL :
          break;
        default :
          break;
      }
    }
  }
}

static void CheckForEmblDdbjID (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr  bsp;
  BoolPtr    isEmblOrDdbj;
  SeqIdPtr   sip;

  if (sep == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    isEmblOrDdbj = (BoolPtr) mydata;
    if (isEmblOrDdbj == NULL) return;
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      switch (sip->choice) {
        case SEQID_EMBL :
        case SEQID_DDBJ :
          *isEmblOrDdbj = TRUE;
          break;
        default :
          break;
      }
    }
  }
}

static void CheckForJournalScanID (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr  bsp;
  BoolPtr    isJScan;
  SeqIdPtr   sip;

  if (sep == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    isJScan = (BoolPtr) mydata;
    if (isJScan == NULL) return;
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      switch (sip->choice) {
        case SEQID_GIBBSQ :
        case SEQID_GIBBMT :
        case SEQID_GIIM :
          *isJScan = TRUE;
          break;
        default :
          break;
      }
    }
  }
}

NLM_EXTERN void CleanUpSeqLoc (SeqLocPtr slp)

{
  BioseqPtr  bsp;
  SeqLocPtr  curr;
  SeqLocPtr  head;
  SeqLocPtr  last;
  SeqLocPtr  loc;
  SeqLocPtr  next;
  SeqIdPtr   sip;
  SeqIntPtr  sintp;
  SeqPntPtr  spp;
  Int4       swp;
  SeqLocPtr  tail;

  if (slp == NULL) return;

  if (slp->choice == SEQLOC_WHOLE) {
    sip = (SeqIdPtr) slp->data.ptrvalue;
    if (sip != NULL) {
      bsp = BioseqFind (sip);
      if (bsp != NULL) {
        sintp = SeqIntNew ();
        if (sintp != NULL) {
          sintp->from = 0;
          sintp->to = bsp->length - 1;
          sintp->id = sip; /* reuse existing slp->data.ptrvalue, no need to free */
          slp->choice = SEQLOC_INT;
          slp->data.ptrvalue = (Pointer) sintp;
        }
      }
    }
  }

  /* from < to for all intervals */
  loc = SeqLocFindNext (slp, NULL);
  while (loc != NULL) {
    if (loc->choice == SEQLOC_INT) {
      sintp = (SeqIntPtr) loc->data.ptrvalue;
      if (sintp != NULL) {
        if (sintp->from > sintp->to) {
          swp = sintp->from;
          sintp->from = sintp->to;
          sintp->to = swp;
        }
        if (sintp->strand == Seq_strand_both) {
          sintp->strand = Seq_strand_plus;
        } else if (sintp->strand == Seq_strand_both_rev) {
          sintp->strand = Seq_strand_minus;
        }
      }
    } else if (loc->choice == SEQLOC_PNT) {
      spp = (SeqPntPtr) loc->data.ptrvalue;
      if (spp != NULL) {
        if (spp->strand == Seq_strand_both) {
          spp->strand = Seq_strand_plus;
        } else if (spp->strand == Seq_strand_both_rev) {
          spp->strand = Seq_strand_minus;
        }
      }
    }
    loc = SeqLocFindNext (slp, loc);
  }

  if (slp->choice == SEQLOC_PACKED_INT) {
    loc = (SeqLocPtr) slp->data.ptrvalue;
    if (loc == NULL || loc->next != NULL) return;
    /* here seqloc_packed_int points to a single location element, so no need for seqloc_packed_int parent */
    slp->choice = loc->choice;
    slp->data.ptrvalue = (Pointer) loc->data.ptrvalue;
    MemFree (loc);
    return;
  }

  if (slp->choice != SEQLOC_MIX) return;
  loc = (SeqLocPtr) slp->data.ptrvalue;
  if (loc == NULL) return;

  if (loc->next != NULL) {
    /* check for null NULL at beginning */
    if (loc->choice == SEQLOC_NULL) {
      slp->data.ptrvalue = (Pointer) loc->next;
      loc->next = NULL;
      ValNodeFree (loc);
    }
    /* check for null NULL at end */
    loc = (SeqLocPtr) slp->data.ptrvalue;
    last = NULL;
    while (loc->next != NULL) {
      last = loc;
      loc = loc->next;
    }
    if (loc->choice == SEQLOC_NULL && last != NULL) {
      last->next = NULL;
      ValNodeFree (loc);
    }
  }

  loc = (SeqLocPtr) slp->data.ptrvalue;
  if (loc == NULL) return;

  if (loc->next == NULL) {
    /* here seqloc_mix points to a single location element, so no need for seqloc_mix parent */
    slp->choice = loc->choice;
    slp->data.ptrvalue = (Pointer) loc->data.ptrvalue;
    MemFree (loc);
    return;
  }

  /* check for nested seqloc_mix, remove nesting */
  curr = loc;
  last = NULL;
  while (curr != NULL) {
    next = curr->next;
    if (curr->choice == SEQLOC_MIX) {
      head = (SeqLocPtr) curr->data.ptrvalue;
      if (head != NULL) {
        tail = head;
        while (tail->next != NULL) {
          tail = tail->next;
        }
        if (last != NULL) {
          last->next = head;
        }
        tail->next = curr->next;
        curr->next = NULL;
        curr = MemFree (curr);
      }
    } else {
      last = curr;
    }
    curr = next;
  }
}

typedef struct cbloc {
  CodeBreakPtr  cbp;
  Int4          pos;
} CbLoc, PNTR CbLocPtr;

static int LIBCALLBACK SortByCodeBreakLoc (VoidPtr ptr1, VoidPtr ptr2)

{
  CbLocPtr  clp1;
  CbLocPtr  clp2;

  clp1 = (CbLocPtr) ptr1;
  clp2 = (CbLocPtr) ptr2;
  if (clp1 == NULL || clp2 == NULL) return 0;
  if (clp1->pos < clp2->pos) {
    return -1;
  } else if (clp1->pos > clp2->pos) {
    return 1;
  }
  return 0;
}

static CodeBreakPtr SortCodeBreaks (SeqFeatPtr sfp, CodeBreakPtr list)

{
  BioseqPtr     bsp;
  CodeBreakPtr  cbp;
  CbLocPtr      head;
  size_t        count, i;
  Boolean       out_of_order = FALSE;
  Int4          pos;
  SeqLocPtr     slp;

  if (sfp == NULL || list == NULL) return list;
  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) return list;

  for (cbp = list, count = 0; cbp != NULL; cbp = cbp->next, count++) continue;
  if (count < 2) return list;

  head = (CbLocPtr) MemNew (sizeof (CbLoc) * (count + 1));
  if (head == NULL) return list;

  for (cbp = list, i = 0; cbp != NULL && i < count; i++) {
    head [i].cbp = cbp;
    slp = dnaLoc_to_aaLoc (sfp, cbp->loc, TRUE, NULL, TRUE);
    head [i].pos = GetOffsetInBioseq (slp, bsp, SEQLOC_START) + 1;
    SeqLocFree (slp);
    cbp = cbp->next;
  }

  pos = head [0].pos;
  for (i = 1; i < count; i++) {
    if (head [i].pos < pos) {
      out_of_order = TRUE;
    }
    pos = head [i].pos;
  }

  if (out_of_order) {
    StableMergeSort (head, count, sizeof (CbLoc), SortByCodeBreakLoc);

    for (i = 0; i < count; i++) {
      cbp = head [i].cbp;
      cbp->next = head [i + 1].cbp;
    }

    list = head [0].cbp;
  }

  MemFree (head);

  return list;
}

static void CleanupDuplicatedCodeBreaks (CodeBreakPtr PNTR prevcbp)

{
  CodeBreakPtr  cbp;
  CodeBreakPtr  last = NULL;
  CodeBreakPtr  next;
  Boolean       unlink;

  if (prevcbp == NULL) return;
  cbp = *prevcbp;
  while (cbp != NULL) {
    next = cbp->next;
    unlink = FALSE;
    if (last != NULL) {
      if (SeqLocCompare (cbp->loc, last->loc) == SLC_A_EQ_B &&
          cbp->aa.choice == last->aa.choice &&
          cbp->aa.value.intvalue == last->aa.value.intvalue) {
        unlink = TRUE;
      }
    } else {
      last = cbp;
    }
    if (unlink) {
      *prevcbp = cbp->next;
      cbp->next = NULL;
      CodeBreakFree (cbp);
    } else {
      last = cbp;
      prevcbp = (CodeBreakPtr PNTR) &(cbp->next);
    }
    cbp = next;
  }
}


CharPtr ncrnaClassList[] = {
"antisense_RNA",
"autocatalytically_spliced_intron",
"hammerhead_ribozyme",
"ribozyme",
"RNase_P_RNA",
"RNase_MRP_RNA",
"telomerase_RNA",
"guide_RNA",
"rasiRNA",
"scRNA",
"siRNA",
"miRNA",
"piRNA",
"snoRNA",
"snRNA",
"SRP_RNA",
"vault_RNA",
"Y_RNA",
"other",
NULL};

Int4 NcrnaOTHER = sizeof (ncrnaClassList) / sizeof (CharPtr) - 1;


extern Boolean IsStringInNcRNAClassList (CharPtr str)
{
  CharPtr PNTR p;

  if (StringHasNoText (str)) return FALSE;
  for (p = ncrnaClassList; *p != NULL; p++)
  {
    if (StringICmp (str, *p) == 0)
    {
      return TRUE;
    }
  }
  return FALSE;
}


static void AddNonCopiedQual (SeqFeatPtr sfp, CharPtr qual, CharPtr class_val)
{
  GBQualPtr gbq;

  if (sfp == NULL || StringHasNoText (qual) || StringHasNoText (class_val)) 
  {
    return;
  }
  gbq = sfp->qual;
  while (gbq != NULL 
          && (StringCmp (gbq->qual, qual) != 0
              || StringCmp (gbq->val, class_val) != 0))
  {
    gbq = gbq->next;
  }
  if (gbq == NULL)
  {
    gbq = GBQualNew ();
    gbq->qual = StringSave (qual);
    gbq->val = StringSave (class_val);
    gbq->next = sfp->qual;
    sfp->qual = gbq;
  }

}


static CharPtr GetMiRNAProduct (CharPtr str)
{
  Int4    len;
  CharPtr product = NULL;

  if (StringHasNoText (str)) return NULL;
  if (StringNCmp (str, "miRNA ", 6) == 0)
  {
    product = StringSave (str + 6);
  }
  else if (StringNCmp (str, "microRNA ", 9) == 0)
  {
    product = StringSave (str + 9);
  }
  else
  {
    len = StringLen (str);
    if (len > 6 && StringCmp (str + len - 6, " miRNA") == 0
        && (len < 15 || StringCmp (str + len - 15, "precursor miRNA") != 0))
    {
      product = (CharPtr) MemNew (sizeof (Char) * (len - 5));
      StringNCpy (product, str, len - 6);
      product[len - 6] = 0;
    }
    else if (len > 9 && StringCmp (str + len - 9, " microRNA") == 0
             && (len < 18 || StringCmp (str + len - 18, "precursor microRNA") != 0))
    {
      product = (CharPtr) MemNew (sizeof (Char) * (len - 8));
      StringNCpy (product, str, len - 9);
      product[len - 9] = 0;
    }
  }
  return product;
}


static Boolean ConvertToNcRNA (SeqFeatPtr sfp)
{
  GBQualPtr gbq;
  RnaRefPtr rrp;
  Boolean was_converted = FALSE;
  CharPtr miRNAproduct = NULL;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL)
  {
    return FALSE;
  }
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp->type == 5 || rrp->type == 6 || rrp->type == 7)
  {
    if (rrp->type == 5)
    {
      AddNonCopiedQual (sfp, "ncRNA_class", "snRNA");
    }
    else if (rrp->type == 6)
    {
      AddNonCopiedQual (sfp, "ncRNA_class", "scRNA");
    }
    else if (rrp->type == 7)
    {
      AddNonCopiedQual (sfp, "ncRNA_class", "snoRNA");
    }
    if (rrp->ext.choice == 1)
    {
      AddNonCopiedQual (sfp, "product", rrp->ext.value.ptrvalue);
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
    }
    rrp->ext.choice = 1;
    rrp->ext.value.ptrvalue = StringSave ("ncRNA");
    rrp->type = 255;
    was_converted = TRUE;
  }
  else if (rrp->type == 255 && rrp->ext.choice == 1)
  {
    if (IsStringInNcRNAClassList (rrp->ext.value.ptrvalue)) 
    {
      AddNonCopiedQual (sfp, "ncRNA_class", rrp->ext.value.ptrvalue);
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
      rrp->ext.value.ptrvalue = StringSave ("ncRNA");
      was_converted = TRUE;
    }
    else if ((miRNAproduct = GetMiRNAProduct (rrp->ext.value.ptrvalue)) != NULL)
    {
      AddNonCopiedQual (sfp, "ncRNA_class", "miRNA");
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
      rrp->ext.value.ptrvalue = StringSave ("ncRNA");
      AddNonCopiedQual (sfp, "product", miRNAproduct);
      miRNAproduct = MemFree (miRNAproduct);
      was_converted = TRUE;
    }
    else if (StringCmp (rrp->ext.value.ptrvalue, "ncRNA") != 0
             && StringCmp (rrp->ext.value.ptrvalue, "tmRNA") != 0
             && StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") != 0)
    {
      AddNonCopiedQual (sfp, "product", rrp->ext.value.ptrvalue);
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
      rrp->ext.value.ptrvalue = StringSave ("misc_RNA");
      was_converted = TRUE;
    }
  }
  if (rrp->type == 255 && rrp->ext.choice == 0) {
    rrp->ext.choice = 1;
    rrp->ext.value.ptrvalue = StringSave ("misc_RNA");
  }
  if (rrp->type == 255 && rrp->ext.choice == 1 &&
      StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") == 0) {
    for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
      if (StringCmp (gbq->qual, "ncRNA_class") == 0) {
        rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
        rrp->ext.value.ptrvalue = StringSave ("ncRNA");
        was_converted = TRUE;
      } else if (StringCmp (gbq->qual, "tag_peptide") == 0) {
        rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
        rrp->ext.value.ptrvalue = StringSave ("tmRNA");
        was_converted = TRUE;
      }
    }
  }
  return was_converted;
}

static void ModernizeFeatureStrings (SeqFeatPtr sfp, Boolean isEmblOrDdbj)

{
  CharPtr      desc;
  GBQualPtr    gbq;
  CharPtr      name;
  ProtRefPtr   prp;
  RnaRefPtr    rrp;
  CharPtr      str;
  ValNodePtr   vnp;

  if (sfp == NULL) return;

  /* skip feature types that do not use data.value.ptrvalue */
  switch (sfp->data.choice) {
    case SEQFEAT_COMMENT:
    case SEQFEAT_BOND:
    case SEQFEAT_SITE:
    case SEQFEAT_PSEC_STR:
      return;
    default:
      break;
  }

  if (sfp->data.value.ptrvalue == NULL) return;

  switch (sfp->data.choice) {
    case SEQFEAT_PROT:
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      desc = prp->desc;
      if (! isEmblOrDdbj) {
        CleanVisStringList (&(prp->name));
        break;
      }
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (StringICmp (str, "RbcL") == 0 || StringICmp (str, "rubisco large subunit") == 0) {
          vnp->data.ptrvalue = StringSave ("ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit");
          str = MemFree (str);
          if (StringICmp (desc, "RbcL") == 0 || StringICmp (desc, "rubisco large subunit") == 0) {
            prp->desc = MemFree (prp->desc);
          }
        } else if (StringICmp (str, "RbcS") == 0 || StringICmp (str, "rubisco small subunit") == 0) {
          vnp->data.ptrvalue = StringSave ("ribulose-1,5-bisphosphate carboxylase/oxygenase small subunit");
          str = MemFree (str);
          if (StringICmp (desc, "RbcS") == 0 || StringICmp (desc, "rubisco small subunit") == 0) {
            prp->desc = MemFree (prp->desc);
          }
        } else if (StringCmp (desc, str) == 0) {
          prp->desc = MemFree (prp->desc);
        }
        if (StringStr (str, "ribulose") != NULL &&
            StringStr (str, "bisphosphate") != NULL &&
            StringStr (str, "methyltransferase") == NULL &&
            StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit") != 0 &&
            StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase small subunit") != 0) {
          if (StringICmp (str, "ribulose 1,5-bisphosphate carboxylase/oxygenase large subunit") == 0 ||
              StringICmp (str, "ribulose 1,5-bisphosphate carboxylase large subunit") == 0 ||
              StringICmp (str, "ribulose bisphosphate carboxylase large subunit") == 0 ||
              StringICmp (str, "ribulose-bisphosphate carboxylase large subunit") == 0 ||
              StringICmp (str, "ribulose-1,5-bisphosphate carboxylase large subunit") == 0 ||
              StringICmp (str, "ribulose-1,5-bisphosphate carboxylase, large subunit") == 0 ||
              StringICmp (str, "large subunit of ribulose-1,5-bisphosphate carboxylase/oxygenase") == 0 ||
              StringICmp (str, "ribulose-1,5-bisphosphate carboxylase oxygenase large subunit") == 0 ||
              StringICmp (str, "ribulose bisphosphate carboxylase large chain") == 0 ||
              StringICmp (str, "ribulose 1,5-bisphosphate carboxylase-oxygenase large subunit") == 0 ||
              StringICmp (str, "ribulose bisphosphate carboxylase oxygenase large subunit") == 0 ||
              StringICmp (str, "ribulose 1,5 bisphosphate carboxylase large subunit") == 0 ||
              StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase, large subunit") == 0 ||
              StringICmp (str, "large subunit of ribulose-1,5-bisphosphate carboxylase/oxgenase") == 0 ||
              StringICmp (str, "ribulose bisphosphate carboxylase/oxygenase large subunit") == 0 ||
              StringICmp (str, "ribulose-1,5-bisphosphate carboxylase oxygenase, large subunit") == 0 ||
              StringICmp (str, "ribulose 5-bisphosphate carboxylase, large subunit") == 0 ||
              StringICmp (str, "ribulosebisphosphate carboxylase large subunit") == 0 ||
              StringICmp (str, "ribulose bisphosphate large subunit") == 0 ||
              StringICmp (str, "ribulose 1,5 bisphosphate carboxylase/oxygenase large subunit") == 0 ||
              StringICmp (str, "ribulose 1,5-bisphosphate carboxylase/oxygenase large chain") == 0 ||
              StringICmp (str, "large subunit ribulose-1,5-bisphosphate carboxylase/oxygenase") == 0 ||
              StringICmp (str, "ribulose-bisphosphate carboxylase, large subunit") == 0 ||
              StringICmp (str, "ribulose-1, 5-bisphosphate carboxylase/oxygenase large-subunit") == 0) {
            vnp->data.ptrvalue = StringSave ("ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit");
            str = MemFree (str);
          }
        }
      }
      CleanVisStringList (&(prp->name));
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp->type == 255 && rrp->ext.choice == 1) {
        name = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringCmp (name, "misc_RNA") == 0) {
          for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
            if (StringCmp (gbq->qual, "product") != 0) continue;
            name = gbq->val;
            if (StringHasNoText (name)) continue;
            if (StringICmp (name, "its1") == 0 || StringICmp (name, "its 1") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 1");
            } else if (StringICmp (name, "its2") == 0 || StringICmp (name, "its 2") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 2");
            } else if (StringICmp (name, "its3") == 0 || StringICmp (name, "its 3") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 3");
            } else if (StringICmp (name, "Ribosomal DNA internal transcribed spacer 1") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 1");
            } else if (StringICmp (name, "Ribosomal DNA internal transcribed spacer 2") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 2");
            } else if (StringICmp (name, "Ribosomal DNA internal transcribed spacer 3") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 3");
            } else if (StringICmp (name, "internal transcribed spacer 1 (ITS1)") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 1");
            } else if (StringICmp (name, "internal transcribed spacer 2 (ITS2)") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 2");
            } else if (StringICmp (name, "internal transcribed spacer 3 (ITS3)") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 3");
            }
          }
        }
      }
      break;
    default:
      break;
  }
}

static Boolean IsFeatureCommentRedundant (SeqFeatPtr sfp)

{
  Uint1            aa;
  Choice           cbaa;
  CodeBreakPtr     cbp;
  CharPtr          comment;
  CdRegionPtr      crp;
  SeqFeatPtr       feat;
  Uint1            from;
  GBQualPtr        gbq;
  GeneRefPtr       grp;
  CharPtr          name;
  BioseqPtr        prod;
  ProtRefPtr       prp;
  Uint1            residue;
  RNAGenPtr        rgp;
  RNAQualPtr       rqp;
  RnaRefPtr        rrp;
  SeqAnnotPtr      sap;
  SeqCodeTablePtr  sctp;
  Uint1            seqcode;
  SeqIdPtr         sip;
  SeqMapTablePtr   smtp;
  CharPtr          str;
  tRNAPtr          trp;
  ValNodePtr       vnp;

  if (sfp == NULL) return FALSE;
  comment = sfp->comment;
  if (StringHasNoText (comment)) return FALSE;

  if (sfp->excpt && StringDoesHaveText (sfp->except_text)) {
    if (StringCmp (comment, sfp->except_text) == 0) return TRUE;
  }

  /* skip feature types that do not use data.value.ptrvalue */
  switch (sfp->data.choice) {
    case SEQFEAT_COMMENT:
    case SEQFEAT_BOND:
    case SEQFEAT_SITE:
    case SEQFEAT_PSEC_STR:
      return FALSE;
    default:
      break;
  }

  if (sfp->data.value.ptrvalue == NULL) return FALSE;

  switch (sfp->data.choice) {
    case SEQFEAT_GENE:
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (StringCmp (comment, grp->locus) == 0) return TRUE;
      if (StringCmp (comment, grp->desc) == 0) return TRUE;
      if (StringCmp (comment, grp->locus_tag) == 0) return TRUE;
      for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (StringCmp (comment, str) == 0) return TRUE;
      }
      break;
    case SEQFEAT_CDREGION:
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
        seqcode = 0;
        sctp = NULL;
        cbaa = cbp->aa;
        switch (cbaa.choice) {
          case 1 :
            seqcode = Seq_code_ncbieaa;
            break;
          case 2 :
            seqcode = Seq_code_ncbi8aa;
            break;
          case 3 :
            seqcode = Seq_code_ncbistdaa;
            break;
          default :
            break;
        }
        if (seqcode != 0) {
          sctp = SeqCodeTableFind (seqcode);
          if (sctp != NULL) {
            residue = cbaa.value.intvalue;
            if (residue != 42) {
              if (seqcode != Seq_code_ncbieaa) {
                smtp = SeqMapTableFind (seqcode, Seq_code_ncbieaa);
                residue = SeqMapTableConvert (smtp, residue);
              }
              if (residue == 'U') {
                if (StringCmp (comment, "selenocysteine") == 0) return TRUE;
              } else if (residue == 'O') {
                if (StringCmp (comment, "pyrrolysine") == 0) return TRUE;
              }
            }
          }
        }
      }
      if (sfp->product != NULL) {
        sip = SeqLocId (sfp->product);
        if (sip != NULL) {
          prod = BioseqFind (sip);
          if (prod != NULL) {
            for (sap = prod->annot; sap != NULL; sap = sap->next) {
              if (sap->type != 1) continue;
              for (feat = (SeqFeatPtr) sap->data; feat != NULL; feat = feat->next) {
                if (feat->data.choice != SEQFEAT_PROT) continue;
                prp = (ProtRefPtr) feat->data.value.ptrvalue;
                if (prp == NULL) continue;
                for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
                  str = (CharPtr) vnp->data.ptrvalue;
                  if (StringHasNoText (str)) continue;
                  if (StringCmp (comment, str) == 0) return TRUE;
                }
              }
            }
          }
        }
      }
      break;
    case SEQFEAT_PROT:
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (StringCmp (comment, str) == 0) return TRUE;
      }
      if (StringDoesHaveText (prp->desc)) {
        if (StringCmp (comment, prp->desc) == 0) return TRUE;
      }
      for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (StringCmp (comment, str) == 0) return TRUE;
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp->type == 255 && rrp->ext.choice == 1) {
        name = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringCmp (name, "misc_RNA") == 0) {
          for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
            if (StringCmp (gbq->qual, "product") != 0) continue;
            name = gbq->val;
            if (StringHasNoText (name)) continue;
            /*
            if (StringICmp (name, "internal transcribed spacer 1") == 0) {
              if (StringICmp (comment, "its1") == 0 || StringICmp (comment, "its 1") == 0) return TRUE;
            } else if (StringICmp (name, "internal transcribed spacer 2") == 0) {
              if (StringICmp (comment, "its2") == 0 || StringICmp (comment, "its 2") == 0) return TRUE;
            } else if (StringICmp (name, "internal transcribed spacer 3") == 0) {
              if (StringICmp (comment, "its3") == 0 || StringICmp (comment, "its 3") == 0) return TRUE;
            }
            */
          }
        }
      } else if (rrp->type == 3 && rrp->ext.choice == 2) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        if (trp != NULL) {
          aa = 0;
          if (trp->aatype == 2) {
            aa = trp->aa;
          } else {
            from = 0;
            switch (trp->aatype) {
              case 0 :
                from = 0;
                break;
              case 1 :
                from = Seq_code_iupacaa;
                break;
              case 2 :
                from = Seq_code_ncbieaa;
                break;
              case 3 :
                from = Seq_code_ncbi8aa;
                break;
              case 4 :
                from = Seq_code_ncbistdaa;
                break;
              default:
                break;
            }
            seqcode = Seq_code_ncbieaa;
            smtp = SeqMapTableFind (seqcode, from);
            if (smtp != NULL) {
              aa = SeqMapTableConvert (smtp, trp->aa);
              if (aa == 255 && from == Seq_code_iupacaa) {
                if (trp->aa == 'U') {
                  aa = 'U';
                } else if (trp->aa == 'O') {
                  aa = 'O';
                }
              }
            }
          }
          if (aa > 0 && aa != 255) {
            if (StringNCmp (comment, "aa: ", 4) == 0) {
              comment += 4;
            }
            residue = FindTrnaAA3 (comment);
            if (residue == aa) {
              if (aa == 'M' && StringICmp ("fMet", comment) == 0) return FALSE;
              return TRUE;
            }
            residue = FindTrnaAA (comment);
            if (residue == aa) return TRUE;
          }
        }
      } else if (rrp->ext.choice == 3) {
        rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
        if (rgp != NULL) {
          if (StringCmp (comment, rgp->product) == 0) return TRUE;
          if (StringCmp (comment, rgp->_class) == 0) return TRUE;
          for (rqp = rgp->quals; rqp != NULL; rqp = rqp->next) {
            if (StringCmp (comment, rqp->val) == 0) return TRUE;
          }
        }
      }
      break;
    default:
      break;
  }

  return FALSE;
}


static CharPtr ExtractSatelliteFromComment (CharPtr comment)
{
  CharPtr satellite_type = NULL, satellite_start = NULL;
  CharPtr satellite_qual = NULL;
  Int4    satellite_len, i;

  if (StringHasNoText (comment)) {
    return NULL;
  }

  if (StringNCmp (comment, "microsatellite", 14) == 0) { 
    satellite_type = "microsatellite";
    satellite_start = comment;
  } else if (StringNCmp (comment, "minisatellite", 13) == 0) {
    satellite_type = "minisatellite";
    satellite_start = comment;
  } else if (StringNCmp (comment, "satellite", 9) == 0) {
    satellite_type = "satellite";
    satellite_start = comment;
  }

  if (satellite_start == NULL) {
    return NULL;
  }

  satellite_len = StringLen (satellite_type);
  if (comment[satellite_len] == '\0') {
    satellite_qual = StringSave (satellite_type);
    *comment = 0;
  } else if (comment[satellite_len] == ';') {
    satellite_qual = StringSave (satellite_type);
    for (i = 0; i <= satellite_len; i++) {
      comment [i] = ' ';
    }
    TrimSpacesAroundString (comment);
  }
  if (comment != NULL && comment [0] == '~' && comment [1] != '~') {
    comment [0] = ' ';
    TrimSpacesAroundString (comment);
  }

  return satellite_qual;
}

static void DoModernizeRNAFields (SeqFeatPtr sfp)

{
  RNAQualSetPtr       nextrqp;
  RNAQualSetPtr PNTR  prevrqp;
  RNAGenPtr           rgp;
  RNAQualSetPtr       rqp;
  RnaRefPtr           rrp;
  CharPtr             str;
  Boolean             unlink;
  Int2                i;
  size_t              len;
  CharPtr             ncclass;
  CharPtr             product;
  CharPtr             tmp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return;

  ModernizeRNAFields (sfp);
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL) return;

  if (rrp->ext.choice == 1 && rrp->type == 10) {
    str = rrp->ext.value.ptrvalue;
    if (StringHasNoText (str)) return;

    rgp = (RNAGenPtr) MemNew (sizeof (RNAGen));
    if (rgp == NULL) return;
    rrp->ext.choice = 3;
    rrp->ext.value.ptrvalue = (Pointer) rgp;
    rgp->product = str;
  }

  if (rrp->ext.choice != 3) return;

  rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
  if (rgp == NULL) return;

  rqp = rgp->quals;
  prevrqp = (RNAQualSetPtr PNTR) &(rgp->quals);
  while (rqp != NULL) {
    nextrqp = rqp->next;
    unlink = FALSE;
    if (StringHasNoText (rqp->qual) || StringHasNoText (rqp->val)) {
      unlink = TRUE;
    }
    if (unlink) {
      *(prevrqp) = rqp->next;
      rqp->next = NULL;
      RNAQualFree (rqp);
    } else {
      prevrqp = (RNAQualSetPtr PNTR) &(rqp->next);
    }
    rqp = nextrqp;
  }

  if (rrp->type == 10 && StringDoesHaveText (rgp->product) && rgp->_class == NULL) {
    ncclass = rgp->product;
    for (i = 0; ncrnaClassList [i] != NULL; i++) {
      str = ncrnaClassList [i];
      if (StringHasNoText (str)) continue;
      len = StringLen (str);
      if (len < 1) continue;
      if (StringNICmp (ncclass, str, len) != 0) continue;
      if (ncclass [len] != ' ') continue;
      tmp = ncclass + len + 1;
      if (StringHasNoText (tmp)) continue;
      ncclass [len] = '\0';
      rgp->_class = StringSave (ncclass);
      product = StringSave (tmp);
      rgp->product = MemFree (rgp->product);
      rgp->product = product;
      TrimSpacesAroundString (rgp->_class);
      TrimSpacesAroundString (rgp->product);
      rrp->type = 8;
      sfp->idx.subtype = FEATDEF_ncRNA;
    }
  }

  if (rgp->quals != NULL) return;
  if (StringDoesHaveText (rgp->_class) || StringDoesHaveText (rgp->product)) return;

  rrp->ext.value.ptrvalue = NULL;
  rrp->ext.choice = 0;
  RNAGenFree (rgp);
}


static void FixncRNAClass (SeqFeatPtr sfp)
{
  RnaRefPtr rrp;
  RNAGenPtr rgp;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_ncRNA
      || (rrp = (RnaRefPtr) sfp->data.value.ptrvalue) == NULL
      || rrp->ext.choice != 3
      || (rgp = (RNAGenPtr) rrp->ext.value.ptrvalue) == NULL)
  {
    return;
  }

  if (StringICmp (rgp->_class, "antisense") == 0) {
    rgp->_class = MemFree (rgp->_class);
    rgp->_class = StringSave ("antisense_RNA");
  }
}


NLM_EXTERN void CleanUpSeqFeat (
  SeqFeatPtr sfp,
  Boolean isEmblOrDdbj,
  Boolean isJscan,
  Boolean stripSerial,
  Boolean modernizeFeats,
  ValNodePtr PNTR publist
)

{
  BioseqPtr     bsp;
  CodeBreakPtr  cbp;
  CdRegionPtr   crp;
  GBQualPtr     gbq;
  Boolean       emptyRNA;
  GeneRefPtr    grp;
  Boolean       hasNulls;
  SeqIdPtr      id;
  ImpFeatPtr    ifp;
  Int2          j;
  CharPtr       name;
  Boolean       partial5;
  Boolean       partial3;
  Uint1         processed;
  ProtRefPtr    prp;
  ValNodePtr    psp;
  RNAGenPtr     rgp;
  RNAQualPtr    rqp;
  RnaRefPtr     rrp;
  Uint1         rrptype;
  CharPtr       satellite_type;
  SeqIntPtr     sintp;
  SeqPntPtr     pntp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  CharPtr       str;
  Uint1         strand;
  tRNAPtr       trp;
  SeqFeatXrefPtr  xref, next, PNTR prevlink;

  if (sfp == NULL) return;
  crp = NULL;
  if (sfp->data.choice == SEQFEAT_IMP) {
    ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
    if (ifp != NULL) {
      if (ifp->loc != NULL) {
        str = StringStr (ifp->loc, "replace");
        if (str != NULL) {
          AddReplaceQual (sfp, str);
          ifp->loc = MemFree (ifp->loc);
        }
      }
      if (StringCmp (ifp->key, "CDS") == 0) {
        if (! isEmblOrDdbj) {
          sfp->data.value.ptrvalue = ImpFeatFree (ifp);
          sfp->data.choice = SEQFEAT_CDREGION;
          crp = CdRegionNew ();
          sfp->data.value.ptrvalue = crp;
          sfp->idx.subtype = FEATDEF_CDS;
        }
      } else if (StringCmp (ifp->key, "allele") == 0 ||
                 StringCmp (ifp->key, "mutation") == 0) {
        ifp->key = MemFree (ifp->key);
        ifp->key = StringSave ("variation");
        sfp->idx.subtype = FEATDEF_variation;
      } else if (StringCmp (ifp->key, "Import") == 0 ||
                 StringCmp (ifp->key, "virion") == 0) {
        ifp->key = MemFree (ifp->key);
        ifp->key = StringSave ("misc_feature");
        sfp->idx.subtype = FEATDEF_misc_feature;
      } else if (StringCmp (ifp->key, "repeat_unit") == 0 ) {
        ifp->key = MemFree (ifp->key);
        ifp->key = StringSave ("repeat_region");
        sfp->idx.subtype = FEATDEF_repeat_region;
      } else if (StringCmp (ifp->key, "misc_bind") == 0) {
        ifp->key = MemFree (ifp->key);
        ifp->key = StringSave ("misc_binding");
        sfp->idx.subtype = FEATDEF_misc_binding;
      } else if (StringCmp (ifp->key, "satellite") == 0 && (! isEmblOrDdbj)) {
        ifp->key = MemFree (ifp->key);
        ifp->key = StringSave ("repeat_region");
        sfp->idx.subtype = FEATDEF_repeat_region;
        gbq = GBQualNew ();
        if (gbq != NULL) {
          gbq->qual = StringSave ("satellite");
          gbq->val = ExtractSatelliteFromComment (sfp->comment);
          if (gbq->val == NULL) {
            gbq->val = StringSave ("satellite");
          }
          gbq->next = sfp->qual;
          sfp->qual = gbq;
        }
      } else if (StringHasNoText (ifp->loc)) {
        rrptype = 0;
        if (StringCmp (ifp->key, "precursor_RNA") == 0) {
          rrptype = 1;
        } else if (StringCmp (ifp->key, "mRNA") == 0) {
          rrptype = 2;
        } else if (StringCmp (ifp->key, "tRNA") == 0) {
          rrptype = 3;
        } else if (StringCmp (ifp->key, "rRNA") == 0) {
          rrptype = 4;
        } else if (StringCmp (ifp->key, "snRNA") == 0) {
          rrptype = 5;
        } else if (StringCmp (ifp->key, "scRNA") == 0) {
          rrptype = 6;
        } else if (StringCmp (ifp->key, "snoRNA") == 0) {
          rrptype = 7;
        } else if (StringCmp (ifp->key, "misc_RNA") == 0) {
          rrptype = 255;
        }
        if (rrptype != 0) {
          sfp->data.value.ptrvalue = ImpFeatFree (ifp);
          sfp->data.choice = SEQFEAT_RNA;
          rrp = RnaRefNew ();
          sfp->data.value.ptrvalue = rrp;
          rrp->type = rrptype;
          sfp->idx.subtype = FindFeatDefType (sfp);
        } else {
          processed = 0;
          if (StringCmp (ifp->key, "proprotein") == 0 || StringCmp (ifp->key, "preprotein") == 0) {
            processed = 1;
          } else if (StringCmp (ifp->key, "mat_peptide") == 0) {
            processed = 2;
          } else if (StringCmp (ifp->key, "sig_peptide") == 0) {
            processed = 3;
          } else if (StringCmp (ifp->key, "transit_peptide") == 0) {
            processed = 4;
          }
          if (processed != 0 || StringCmp (ifp->key, "Protein") == 0) {
            bsp = BioseqFind (SeqLocId (sfp->location));
            if (bsp != NULL && ISA_aa (bsp->mol)) {
              sfp->data.value.ptrvalue = ImpFeatFree (ifp);
              sfp->data.choice = SEQFEAT_PROT;
              prp = ProtRefNew ();
              sfp->data.value.ptrvalue = prp;
              prp->processed = processed;
              sfp->idx.subtype = FindFeatDefType (sfp);
            }
          }
        }
      }
      if (sfp->data.choice == SEQFEAT_IMP && StringCmp (ifp->key, "repeat_region") == 0 && (! isEmblOrDdbj)) {
        satellite_type = ExtractSatelliteFromComment (sfp->comment);
        if (satellite_type != NULL) {
          gbq = GBQualNew ();
          if (gbq != NULL) {
            gbq->qual = StringSave ("satellite");
            gbq->val = satellite_type;
            gbq->next = sfp->qual;
            sfp->qual = gbq;
          }
        }
      }
    }
  }
  if (crp != NULL && crp->frame == 0 && (! sfp->pseudo)) {
    crp->frame = GetFrameFromLoc (sfp->location);
  }
  if (sfp->data.choice == SEQFEAT_RNA) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rrp != NULL) {
      if (rrp->ext.choice == 1) {
        name = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringHasNoText (name)) {
          rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
          rrp->ext.choice = 0;
        }
      } else if (rrp->ext.choice == 2) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        if (trp != NULL) {
          if (trp->aatype == 0 && trp->aa == 0 && trp->anticodon == NULL) {
            emptyRNA = TRUE;
            for (j = 0; j < 6; j++) {
              if (trp->codon [j] != 255) {
                emptyRNA = FALSE;
              }
            }
            if (emptyRNA) {
              rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
              rrp->ext.choice = 0;
            }
          }
        }
      } else if (rrp->ext.choice == 3) {
        rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
        if (rgp != NULL) {
          if (StringHasNoText (rgp->_class) && StringHasNoText (rgp->product)) {
            emptyRNA = TRUE;
            for (rqp = rgp->quals; rqp != NULL; rqp = rqp->next) {
              if (StringDoesHaveText (rqp->qual) && StringDoesHaveText (rqp->val)) {
                emptyRNA = FALSE;
              }
            } 
            if (emptyRNA) {
              rrp->ext.value.ptrvalue = RNAGenFree (rrp->ext.value.ptrvalue);
              rrp->ext.choice = 0;
            }
          }
        }
      }
    }
  }
  ModernizeFeatureGBQuals (sfp);
  sfp->qual = SortFeatureGBQuals (sfp->qual);
  CleanupDuplicateGBQuals (&(sfp->qual));
  CleanupFeatureGBQuals (sfp, isEmblOrDdbj);
  sfp->qual = SortIllegalGBQuals (sfp->qual);
  CleanupFeatureStrings (sfp, isJscan, stripSerial, modernizeFeats, publist);
  FixOldDbxrefs (sfp->dbxref);
  FixNumericDbxrefs (sfp->dbxref);
  sfp->dbxref = ValNodeSort (sfp->dbxref, SortDbxref);
  CleanupDuplicateDbxrefs (&(sfp->dbxref));
  CleanupObsoleteDbxrefs (&(sfp->dbxref));
  psp = sfp->cit;
  if (psp != NULL && psp->data.ptrvalue) {
    psp->data.ptrvalue = ValNodeSort ((ValNodePtr) psp->data.ptrvalue, SortCits);
    CleanupDuplicateCits ((ValNodePtr PNTR) &(psp->data.ptrvalue));
  }
  CleanUpSeqLoc (sfp->location);
  strand = SeqLocStrand (sfp->location);
  id = SeqLocId (sfp->location);
  if (sfp->data.choice == SEQFEAT_GENE) {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (grp != NULL) {
      if (grp->pseudo) {
        sfp->pseudo = TRUE;
        grp->pseudo = FALSE;
      }
    }
  } else if (sfp->data.choice == SEQFEAT_CDREGION) {
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if (crp != NULL) {
      crp->code_break = SortCodeBreaks (sfp, crp->code_break);
      CleanupDuplicatedCodeBreaks (&(crp->code_break));
      for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
        CleanUpSeqLoc (cbp->loc);
        if (strand == Seq_strand_minus && id != NULL) {
          slp = cbp->loc;
          if (slp != NULL && slp->choice == SEQLOC_INT) {
            sip = SeqLocId (slp);
            if (sip != NULL && SeqIdComp (id, sip) == SIC_YES) {
              sintp = (SeqIntPtr) slp->data.ptrvalue;
              if (sintp != NULL) {
                sintp->strand = Seq_strand_minus;
              }
            }
          }
        }
      }
    }
  } else if (sfp->data.choice == SEQFEAT_RNA) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rrp != NULL) {
      if (rrp->pseudo) {
        sfp->pseudo = TRUE;
        rrp->pseudo = FALSE;
      }
    }
    if (rrp != NULL && rrp->ext.choice == 2) {
      trp = (tRNAPtr) rrp->ext.value.ptrvalue;
      if (trp != NULL && trp->anticodon != NULL) {
        CleanUpSeqLoc (trp->anticodon);
        if (strand == Seq_strand_minus && id != NULL) {
          slp = trp->anticodon;
          if (slp != NULL && slp->choice == SEQLOC_INT) {
            sip = SeqLocId (slp);
            if (sip != NULL && SeqIdComp (id, sip) == SIC_YES) {
              sintp = (SeqIntPtr) slp->data.ptrvalue;
              if (sintp != NULL) {
                sintp->strand = Seq_strand_minus;
              }
            }
          }
        }
      }
    }
    if (ConvertToNcRNA (sfp)) {
      sfp->idx.subtype = FindFeatDefType (sfp);
    }
    if (sfp->idx.subtype == FEATDEF_ncRNA) {
      FixncRNAClass (sfp);
    }       
  } else if (sfp->data.choice == SEQFEAT_REGION ||
             sfp->data.choice == SEQFEAT_SITE ||
             sfp->data.choice == SEQFEAT_BOND ||
             sfp->data.choice == SEQFEAT_PROT) {
    bsp = BioseqFind (SeqLocId (sfp->location));
    if (bsp != NULL && ISA_aa (bsp->mol)) {
      slp = SeqLocFindNext (sfp->location, NULL);
      while (slp != NULL) {
        if (slp->choice == SEQLOC_INT) {
          sintp = (SeqIntPtr) slp->data.ptrvalue;
          if (sintp != NULL) {
            if (sintp->strand != Seq_strand_unknown) {
              sintp->strand = Seq_strand_unknown;
            }
          }
        } else if (slp->choice == SEQLOC_PNT) {
          pntp = (SeqPntPtr) slp->data.ptrvalue;
          if (pntp->strand != Seq_strand_unknown) {
            pntp->strand = Seq_strand_unknown;
          }
        }
        slp = SeqLocFindNext (sfp->location, slp);
      }
    }
  }

  ModernizeFeatureStrings (sfp, isEmblOrDdbj);

  if (sfp->data.choice == SEQFEAT_GENE) {
    if (modernizeFeats) {
      ModernizeGeneFields (sfp);
    }
  }

  if (sfp->data.choice == SEQFEAT_RNA) {
    if (modernizeFeats) {
      DoModernizeRNAFields (sfp);
    }
  }

  if (IsFeatureCommentRedundant (sfp)) {
    sfp->comment = MemFree (sfp->comment);
  }

  /* sort and unique gbquals again after recent processing */
  sfp->qual = SortFeatureGBQuals (sfp->qual);
  CleanupDuplicateGBQuals (&(sfp->qual));
  sfp->qual = SortIllegalGBQuals (sfp->qual);

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  hasNulls = LocationHasNullsBetween (sfp->location);
  sfp->partial = (sfp->partial || partial5 || partial3 || (hasNulls && ! isEmblOrDdbj));

  prevlink = (SeqFeatXrefPtr PNTR) &(sfp->xref);
  xref = sfp->xref;
  while (xref != NULL) {
    next = xref->next;

    if (xref->id.choice == 0 && xref->data.choice == 0) {
      *prevlink = xref->next;
      xref->next = NULL;
      MemFree (xref);
    } else {
      prevlink = (SeqFeatXrefPtr PNTR) &(xref->next);
    }

    xref = next;
  }
}


static void CleanUpSeqGraph (SeqGraphPtr sgp)

{
  if (sgp == NULL) return;
  if (sgp->loc != NULL) {
    CleanUpSeqLoc (sgp->loc);
  }
}

static void RemoveZeroLengthSeqLits (BioseqPtr bsp)
{
  DeltaSeqPtr dsp, prev = NULL, dsp_next;
  SeqLitPtr slip;

  if (bsp == NULL || bsp->repr != Seq_repr_delta) {
    return;
  }

  for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp != NULL; dsp = dsp_next) {
    dsp_next = dsp->next;
    if (dsp->choice == 2 && (slip = (SeqLitPtr) (dsp->data.ptrvalue)) != NULL 
        && slip->length == 0 && slip->seq_data_type == 1
        && slip->seq_data != NULL) {
      if (prev == NULL) {
        bsp->seq_ext = dsp->next;
      } else {
        prev->next = dsp->next;
      }
      dsp->next = NULL;
      dsp = DeltaSeqFree (dsp);
    } else {
      prev = dsp;
    }
  }
}


static Boolean CleanUpSeqIdText (SeqIdPtr sip)
{
  ObjectIdPtr  oip;
  Boolean      rval = FALSE;

  if (sip == NULL) return FALSE;
  if (sip->choice == SEQID_LOCAL) {
    oip = (ObjectIdPtr) sip->data.ptrvalue;
    if (oip != NULL) {
      if (StringDoesHaveText (oip->str)) {
        if (isspace (oip->str[0]) || isspace (oip->str[StringLen (oip->str) - 1])) {
          TrimSpacesAroundString (oip->str);
          rval = TRUE;
        }
      }
    }
  }
  return rval;
}


static void CleanUpSeqId (
  SeqIdPtr sip,
  Pointer userdata
)

{
  CleanUpSeqIdText (sip);
}

static void CleanSeqIdInBioseq (BioseqPtr bsp, Pointer userdata)

{
  SeqIdPtr sip;
  Boolean  need_reindex = FALSE;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (CleanUpSeqIdText (sip)) {
      need_reindex = TRUE;
    }
  }
  if (need_reindex) {
    SeqMgrReplaceInBioseqIndex (bsp);
  }
}

static void CleanSeqIdInSeqFeat (SeqFeatPtr sfp, Pointer userdata)

{
  VisitSeqIdsInSeqFeat (sfp, NULL, CleanUpSeqId);
}

static void CleanSeqIdInSeqAlign (SeqAlignPtr sap, Pointer userdata)

{
  VisitSeqIdsInSeqAlign (sap, NULL, CleanUpSeqId);
}

static void CleanSeqIdInSeqGraph (SeqGraphPtr sgp, Pointer userdata)

{
  VisitSeqIdsInSeqGraph (sgp, NULL, CleanUpSeqId);
}

static void CleanSeqIdInSeqAnnot (SeqAnnotPtr annot, Pointer userdata)

{
  VisitSeqIdsInSeqAnnot (annot, NULL, CleanUpSeqId);
}

typedef struct npcounts {
  Int4     nucs;
  Int4     prots;
  Boolean  make_genbank;
} NPCounts, PNTR NPCountsPtr;

static void CountNucsAndProts (BioseqPtr bsp, Pointer userdata)

{
  NPCountsPtr  ncp;

  if (bsp == NULL) return;
  ncp = (NPCountsPtr) userdata;
  if (ncp == NULL) return;

  if (ISA_na (bsp->mol)) {
    (ncp->nucs)++;
  } else if (ISA_aa (bsp->mol)) {
    (ncp->prots)++;
  }
}

static void CheckInnerSets (BioseqSetPtr bssp, Pointer userdata)

{
  NPCountsPtr  ncp;

  if (bssp == NULL) return;
  ncp = (NPCountsPtr) userdata;
  if (ncp == NULL) return;

  if (bssp->_class == BioseqseqSet_class_segset || bssp->_class == BioseqseqSet_class_parts) return;
  ncp->make_genbank = TRUE;
}

static void FixBadSetClass (BioseqSetPtr bssp, Pointer userdata)

{
  NPCounts  nc;

  if (bssp == NULL) return;
  if (bssp->_class != BioseqseqSet_class_not_set && bssp->_class != BioseqseqSet_class_other) return;

  MemSet ((Pointer) &nc, 0, sizeof (NPCounts));
  VisitSequencesInSet (bssp, (Pointer) &nc, VISIT_MAINS, CountNucsAndProts);
  VisitSetsInSet (bssp, (Pointer) &nc, CheckInnerSets);
  if (nc.nucs == 1 && nc.prots > 0 && (! nc.make_genbank)) {
    bssp->_class = BioseqseqSet_class_nuc_prot;
  } else {
    bssp->_class = BioseqseqSet_class_genbank;
  }
}

static void RemoveDuplicateSeqIds (BioseqPtr bsp)

{
  SeqIdPtr sip, sip_cmp, sip_prev, sip_next;

  if (bsp == NULL) {
    return;
  }

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    sip_prev = sip;
    for (sip_cmp = sip->next; sip_cmp != NULL; sip_cmp = sip_next) {
      sip_next = sip_cmp->next;
      if (SeqIdComp (sip, sip_cmp) == SIC_YES) {
        sip_prev->next = sip_cmp->next;
        sip_cmp->next = NULL;
        sip_cmp = SeqIdFree (sip_cmp);
      } else {
        sip_prev = sip_cmp;
      }
    }
  }
}


static void BasicSeqEntryCleanupInternal (
  SeqEntryPtr sep,
  ValNodePtr PNTR publist,
  Boolean isEmblOrDdbj,
  Boolean isJscan,
  Boolean stripSerial
)

{
  BioSourcePtr  biop;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqDescrPtr   desc;
  Char          div [10];
  GBBlockPtr    gbp;
  MolInfoPtr    mip;
  OrgNamePtr    onp;
  OrgRefPtr     orp;
  SeqAnnotPtr   sap = NULL;
  ValNodePtr    sdp = NULL;
  SeqFeatPtr    sfp;
  SeqGraphPtr   sgp;
  SeqEntryPtr   tmp;

  if (sep == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    /* remove duplicate SeqIds on the same Bioseq */
    RemoveDuplicateSeqIds (bsp);

    /* repair damaged delta sequences */
    RemoveZeroLengthSeqLits (bsp);

    sap = bsp->annot;
    sdp = bsp->descr;
    desc = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
    if (desc != NULL && desc->choice == Seq_descr_molinfo) {
      mip = (MolInfoPtr) desc->data.ptrvalue;
      if (mip != NULL) {
        /* repair if bsp.mol is not-set */
        if (bsp->mol == 0) {
          switch (mip->biomol) {
            case MOLECULE_TYPE_GENOMIC :
              bsp->mol = Seq_mol_na;
              break;
            case MOLECULE_TYPE_PRE_MRNA :
            case MOLECULE_TYPE_MRNA :
            case MOLECULE_TYPE_RRNA :
            case MOLECULE_TYPE_TRNA :
            case MOLECULE_TYPE_SNRNA :
            case MOLECULE_TYPE_SCRNA :
            case MOLECULE_TYPE_CRNA :
            case MOLECULE_TYPE_SNORNA :
            case MOLECULE_TYPE_TRANSCRIBED_RNA :
            case MOLECULE_TYPE_NCRNA :
            case MOLECULE_TYPE_TMRNA :
              bsp->mol = Seq_mol_rna;
              break;
            case MOLECULE_TYPE_PEPTIDE :
              bsp->mol = Seq_mol_aa;
              break;
            case MOLECULE_TYPE_OTHER_GENETIC_MATERIAL :
              bsp->mol = Seq_mol_other;
              break;
            case MOLECULE_TYPE_GENOMIC_MRNA_MIX :
              bsp->mol = Seq_mol_na;
              break;
            default :
              break;
          }
        } else if (bsp->mol != Seq_mol_rna 
                   && (mip->biomol == MOLECULE_TYPE_CRNA || mip->biomol == MOLECULE_TYPE_MRNA)) {
          bsp->mol = Seq_mol_rna;
        }
      }
    }
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
      BasicSeqEntryCleanupInternal (tmp, publist, isEmblOrDdbj, isJscan, stripSerial);
    }
    sap = bssp->annot;
    sdp = bssp->descr;
  } else return;
  biop = NULL;
  orp = NULL;
  gbp = NULL;
  div [0] = '\0';
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        CleanUpSeqFeat (sfp, isEmblOrDdbj, isJscan, stripSerial, TRUE, publist);
        sfp = sfp->next;
      }
    } else if (sap->type == 3) {
      sgp = (SeqGraphPtr) sap->data;
      while (sgp != NULL) {
        CleanUpSeqGraph (sgp);
        sgp = sgp->next;
      }
    }
    sap = sap->next;
  }
  while (sdp != NULL) {
    switch (sdp->choice) {
      case Seq_descr_org :
        orp = (OrgRefPtr) sdp->data.ptrvalue;
        break;
      case Seq_descr_genbank :
        gbp = (GBBlockPtr) sdp->data.ptrvalue;
        break;
      case Seq_descr_source :
        biop = (BioSourcePtr) sdp->data.ptrvalue;
        if (biop != NULL) {
          orp = biop->org;
        }
        break;
      default :
        break;
    }
    CleanupDescriptorStrings (sdp, stripSerial, TRUE, publist, isEmblOrDdbj);
    sdp = sdp->next;
  }

  /* copy genbank block division into biosource, if necessary */

  if (orp != NULL && gbp != NULL) {
    StringNCpy_0 (div, gbp->div, sizeof (div));
    if (StringHasNoText (div)) return;
    onp = orp->orgname;
    while (onp != NULL) {
      if (StringHasNoText (onp->div)) {
        onp->div = MemFree (onp->div);
        onp->div = StringSaveNoNull (div);
      }
      onp = onp->next;
    }
  }
}

static void ReplaceCitOnFeat (CitGenPtr cgp, ValNodePtr publist)

{
  ValNodePtr  nxt;
  ValNodePtr  vnp;

  for (vnp = publist; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice != 1) continue;
    if (StringCmp (cgp->cit, (CharPtr) vnp->data.ptrvalue) == 0) {
      nxt = vnp->next;
      if (nxt != NULL && nxt->choice == 2) {
        cgp->cit = MemFree (cgp->cit);
        cgp->cit = StringSaveNoNull ((CharPtr) nxt->data.ptrvalue);
        if (cgp->cit != NULL) {
          if (StringNICmp (cgp->cit, "unpublished", 11) == 0) {
            cgp->cit [0] = 'U';
          }
        }
      }
      return;
    }
  }
}

static void ChangeCitsOnFeats (SeqFeatPtr sfp, Pointer userdata)

{
  CitGenPtr   cgp;
  ValNodePtr  ppr;
  ValNodePtr  psp;
  ValNodePtr  vnp;

  psp = sfp->cit;
  if (psp != NULL && psp->data.ptrvalue) {
    for (ppr = (ValNodePtr) psp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
      vnp = NULL;
      if (ppr->choice == PUB_Gen) {
        vnp = ppr;
      } else if (ppr->choice == PUB_Equiv) {
        for (vnp = (ValNodePtr) ppr->data.ptrvalue;
             vnp != NULL && vnp->choice != PUB_Gen;
             vnp = vnp->next) continue;
      }
      if (vnp != NULL && vnp->choice == PUB_Gen) {
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL && (! StringHasNoText (cgp->cit))) {
          ReplaceCitOnFeat (cgp, (ValNodePtr) userdata);
        }
      }
    }
  }
}

static Int4 GetPmidForMuid (ValNodePtr pairlist, Int4 muid)

{
  ValNodePtr  vnp;

  vnp = pairlist;
  while (vnp != NULL) {
    if (muid == vnp->data.intvalue) {
      vnp = vnp->next;
      if (vnp == NULL) return 0;
      return vnp->data.intvalue;
    } else {
      vnp = vnp->next;
      if (vnp == NULL) return 0;
      vnp = vnp->next;
    }
  }

  return 0;
}

static void ChangeFeatCitsToPmid (SeqFeatPtr sfp, Pointer userdata)

{
  Int4        muid = 0;
  Int4        pmid = 0;
  ValNodePtr  ppr;
  ValNodePtr  psp;
  ValNodePtr  vnp;

  psp = sfp->cit;
  if (psp != NULL && psp->data.ptrvalue) {
    for (ppr = (ValNodePtr) psp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
      vnp = NULL;
      if (ppr->choice == PUB_Muid) {
        vnp = ppr;
      } else if (ppr->choice == PUB_Equiv) {
        for (vnp = (ValNodePtr) ppr->data.ptrvalue;
             vnp != NULL && vnp->choice != PUB_Muid;
             vnp = vnp->next) continue;
      }
      if (vnp != NULL && vnp->choice == PUB_Muid) {
        muid = vnp->data.intvalue;
        if (muid != 0) {
          pmid = GetPmidForMuid ((ValNodePtr) userdata, muid);
          if (pmid != 0) {
            vnp->choice = PUB_PMid;
            vnp->data.intvalue = pmid;
          }
        }
      }
    }
  }
}

static void GetMuidPmidPairs (PubdescPtr pdp, Pointer userdata)

{
  Int4        muid = 0;
  Int4        pmid = 0;
  ValNodePtr  vnp;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Muid :
        muid = vnp->data.intvalue;
        break;
      case PUB_PMid :
        pmid = vnp->data.intvalue;
        break;
      default :
        break;
    }
  }
  if (muid == 0 || pmid == 0) return;
  ValNodeAddInt ((ValNodePtr PNTR) userdata, 0, muid);
  ValNodeAddInt ((ValNodePtr PNTR) userdata, 0, pmid);
}

static void FlattenPubSet (ValNodePtr PNTR prev)

{
  ValNodePtr  next;
  ValNodePtr  ppr;
  ValNodePtr  vnp;

  if (prev == NULL || *prev == NULL) return;
  ppr = *prev;
  while (ppr != NULL) {
    next = ppr->next;

    if (ppr->choice == PUB_Equiv) {
      vnp = (ValNodePtr) ppr->data.ptrvalue;
      if (vnp != NULL && vnp->next == NULL) {
        ppr->choice = vnp->choice;
        switch (vnp->choice) {
          case PUB_Muid :
          case PUB_PMid :
            ppr->data.intvalue = vnp->data.intvalue;
            break;
          default :
            ppr->data.ptrvalue = vnp->data.ptrvalue;
            break;
        }
        ValNodeFree (vnp);
      }
    }

    ppr = next;
  }
} 

static void FlattenDupInPubSet (ValNodePtr PNTR prev)

{
  ValNodePtr  next;
  ValNodePtr  nxt;
  ValNodePtr  ppr;
  ValNodePtr  vnp;

  if (prev == NULL || *prev == NULL) return;
  ppr = *prev;
  while (ppr != NULL) {
    next = ppr->next;

    if (ppr->choice == PUB_Equiv) {
      vnp = (ValNodePtr) ppr->data.ptrvalue;
      if (vnp != NULL) {
        nxt = vnp->next;
        if (nxt != NULL && nxt->next == NULL && vnp->choice == nxt->choice) {
          switch (vnp->choice) {
            case PUB_Muid :
            case PUB_PMid :
              if (vnp->data.intvalue == nxt->data.intvalue) {
                vnp->next = ValNodeFree (nxt);
              }
              break;
            default :
              break;
          }
        }
      }
    }

    ppr = next;
  }
} 

static void FlattenPubdesc (PubdescPtr pdp, Pointer userdata)

{
  FlattenPubSet (&(pdp->pub));
}

static void FlattenSfpCit (SeqFeatPtr sfp, Pointer userdata)

{
  ValNodePtr  psp;

  psp = sfp->cit;
  if (psp == NULL) return;
  FlattenDupInPubSet ((ValNodePtr PNTR) &(psp->data.ptrvalue));
  FlattenPubSet ((ValNodePtr PNTR) &(psp->data.ptrvalue));
}

typedef struct fastnode {
  ValNodePtr  head;
  ValNodePtr  tail;
} FastNode, PNTR FastNodePtr;

static void GetCitGenLabels (PubdescPtr pdp, Pointer userdata)

{
  Char             buf [121];
  CitGenPtr        cgp;
  FastNodePtr      labellist;
  ValNodePtr       tmp;
  ValNodePtr       vnp;
 
  if (pdp == NULL) return;
  labellist = (FastNodePtr) userdata;
  if (labellist == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice != PUB_Gen) continue;
    cgp = (CitGenPtr) vnp->data.ptrvalue;
    if (cgp == NULL) continue;
    if (cgp->cit == NULL && cgp->journal == NULL &&
        cgp->date == NULL && cgp->serial_number) continue;
    PubLabelUnique (vnp, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE);
    tmp = ValNodeCopyStr (&(labellist->tail), 0, buf);
    if (labellist->head == NULL) {
      labellist->head = tmp;
    }
    labellist->tail = tmp;
  }
}

static void ReplaceShortCitGenOnFeat (CitGenPtr cgp, ValNodePtr labellist)

{
  Char        buf [128];
  Char        ch;
  size_t      len1;
  size_t      len2;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     tmp;
  ValNodePtr  vnp;

  for (vnp = labellist; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    len1 = StringLen (cgp->cit);
    if (len1 < 2 || len1 > 120) continue;
    StringCpy (buf, cgp->cit);
    ptr = StringStr (buf, "Unpublished");
    if (ptr != NULL) {
      ptr += 11;
      *ptr = '\0';
      tmp = StringStr (cgp->cit, "Unpublished");
      if (tmp != NULL) {
        tmp += 11;
        ch = *tmp;
        while (ch == ' ') {
          tmp++;
          ch = *tmp;
        }
        StringCat (buf, tmp);
      }
    }
    len1 = StringLen (buf);
    if (buf [len1 - 1] != '>') continue;
    len1--;
    len2 = StringLen (str);
    if (len1 >= len2) continue;
    if (StringNCmp (str, buf, len1) == 0) {
      cgp->cit = MemFree (cgp->cit);
      cgp->cit = StringSaveNoNull (str);
      if (cgp->cit != NULL) {
        if (StringNICmp (cgp->cit, "unpublished", 11) == 0) {
          cgp->cit [0] = 'U';
        }
      }
      return;
    }
  }
}

static void UpdateShortFeatCits (SeqFeatPtr sfp, Pointer userdata)

{
  CitGenPtr   cgp;
  ValNodePtr  ppr;
  ValNodePtr  psp;
  ValNodePtr  vnp;

  psp = sfp->cit;
  if (psp != NULL && psp->data.ptrvalue) {
    for (ppr = (ValNodePtr) psp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
      vnp = NULL;
      if (ppr->choice == PUB_Gen) {
        vnp = ppr;
      } else if (ppr->choice == PUB_Equiv) {
        for (vnp = (ValNodePtr) ppr->data.ptrvalue;
             vnp != NULL && vnp->choice != PUB_Gen;
             vnp = vnp->next) continue;
      }
      if (vnp != NULL && vnp->choice == PUB_Gen) {
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL && (! StringHasNoText (cgp->cit))) {
          ReplaceShortCitGenOnFeat (cgp, (ValNodePtr) userdata);
        }
      }
    }
  }
}

NLM_EXTERN void BasicSeqAnnotCleanup (SeqAnnotPtr sap)

{
  SeqFeatPtr   sfp;
  SeqGraphPtr  sgp;

  if (sap == NULL) return;

  VisitSeqIdsInSeqAnnot (sap, NULL, CleanUpSeqId);

  if (sap->type == 1) {
    sfp = (SeqFeatPtr) sap->data;
    while (sfp != NULL) {
      CleanUpSeqFeat (sfp, FALSE, FALSE, TRUE, TRUE, NULL);
      sfp = sfp->next;
    }
  } else if (sap->type == 3) {
    sgp = (SeqGraphPtr) sap->data;
    while (sgp != NULL) {
      CleanUpSeqGraph (sgp);
      sgp = sgp->next;
    }
  }
}

/*
static CharPtr proteinOrganellePrefix [] = {
  NULL,
  NULL,
  "chloroplast",
  "chromoplast",
  "kinetoplast",
  "mitochondrion",
  "plastid",
  "macronuclear",
  "extrachromosomal",
  "plasmid",
  NULL,
  NULL,
  "cyanelle",
  "proviral",
  "virus",
  "nucleomorph",
  "apicoplast",
  "leucoplast",
  "protoplast",
  "endogenous virus",
  "hydrogenosome",
  "chromosome",
  "chromatophore"
};
*/

static CharPtr proteinOrganellePrefix [] = {
  NULL,
  NULL,
  "chloroplast",
  NULL,
  NULL,
  "mitochondrion",
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
};

static CharPtr TitleEndsInOrganism (
  CharPtr title,
  CharPtr organism,
  CharPtr organelle,
  CharPtr PNTR onlp
)

{
  int      genome;
  size_t   len1, len2, len3;
  CharPtr  onl, ptr, tmp;

  if (onlp != NULL) {
    *onlp = NULL;
  }
  if (StringHasNoText (title) || StringHasNoText (organism)) return NULL;
  len1 = StringLen (title);
  len2 = StringLen (organism);
  if (len2 + 4 >= len1) return NULL;

  tmp = title + len1 - len2 - 3;
  if (tmp [0] != ' ' || tmp [1] != '[' || tmp [len2 + 2] != ']') return NULL;
  if (StringNICmp (tmp + 2, organism, len2) != 0) return NULL;

  if (onlp != NULL) {
    len3 = len1 - len2 - 3;
    for (genome = GENOME_chloroplast; genome <= GENOME_chromatophore; genome++) {
      ptr = proteinOrganellePrefix [genome];
      if (ptr == NULL) continue;
      len2 = StringLen (ptr);
      if (len2 + 4 >= len3) continue;
      onl = title + len3 - len2 - 3;
      if (onl [0] != ' ' || onl [1] != '(' || onl [len2 + 2] != ')') continue;
      if (StringNICmp (onl + 2, ptr, len2) != 0) continue;
      *onlp = onl;
      break;
    }
  }

  return tmp;
}

static void AddPartialToProteinTitle (
  BioseqPtr bsp,
  Pointer userdata
)

{
  BioSourcePtr  biop;
  int           genome = 0;
  size_t        len;
  MolInfoPtr    mip;
  CharPtr       oldname = NULL;
  OrgModPtr     omp;
  OrgNamePtr    onp;
  CharPtr       organelle = NULL;
  OrgRefPtr     orp;
  Boolean       partial = FALSE;
  CharPtr       penult = NULL;
  SeqDescrPtr   sdp;
  SeqIdPtr      sip;
  CharPtr       str;
  CharPtr       suffix = NULL;
  CharPtr       taxname = NULL;
  CharPtr       title;
  CharPtr       tmp;
  SeqDescrPtr   ttl = NULL;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_SWISSPROT) return;
  }

  ttl = BioseqGetSeqDescr (bsp, Seq_descr_title, NULL);
  if (ttl == NULL || ttl->choice != Seq_descr_title) return;
  str = (CharPtr) ttl->data.ptrvalue;
  if (StringHasNoText (str)) return;

  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
  if (sdp != NULL && sdp->choice == Seq_descr_molinfo) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL && mip->completeness > 1 && mip->completeness < 6) {
      partial = TRUE;
    }
  }

  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_source, NULL);
  if (sdp != NULL && sdp->choice == Seq_descr_source) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL) {
      genome = biop->genome;
      if (genome >= GENOME_chloroplast && genome <= GENOME_chromatophore) {
        organelle = proteinOrganellePrefix [genome];
      }
      orp = biop->org;
      if (orp != NULL) {
        taxname = orp->taxname;
        if (StringNICmp (organelle, taxname, StringLen (organelle)) == 0) {
          organelle = NULL;
        }
        onp = orp->orgname;
        if (onp != NULL) {
          for (omp = onp->mod; omp != NULL; omp = omp->next) {
            if (omp->subtype == ORGMOD_old_name) {
              oldname = omp->subname;
            }
          }
        }
      }
    }
  }

  /* search for partial, must be just before parenthesized organelle or bracketed organism */
  tmp = StringSearch (str, ", partial [");
  if (tmp == NULL) {
    tmp = StringSearch (str, ", partial (");
  }

  /* find oldname or taxname in brackets at end of protein title */
  if (oldname != NULL && taxname != NULL) {
    suffix = TitleEndsInOrganism (str, oldname, organelle, &penult);
  }
  if (suffix == NULL && taxname != NULL) {
    suffix = TitleEndsInOrganism (str, taxname, organelle, &penult);
    if (suffix != NULL) {
      if (organelle == NULL && penult != NULL) {
      } else if (organelle != NULL && penult == NULL) {
      } else if (StringCmp (organelle, penult) != 0) {
      } else {
        /* bail if no need to change partial text (organelle) [organism name] */
        if (partial) {
          if (tmp != NULL) return;
        } else {
          if (tmp == NULL) return;
        }
      }
    }
  }

  /* do not change unless [genus species] was at the end */
  if (suffix == NULL) return;

  /* truncate bracketed info from end of title, will replace with current taxname */
  *suffix = '\0';
  suffix = taxname;

  /* truncate parenthesized info from just before bracketed taxname, will replace with current organelle */
  if (penult != NULL) {
    *penult = '\0';
  }

  /* if ", partial [/(" was indeed just before the [genus species] or (organelle), it will now be ", partial" */
  if (! partial && tmp != NULL && StringCmp (tmp, ", partial") == 0) {
    *tmp = '\0';
  }
  TrimSpacesAroundString (str);

  len = StringLen (str) + StringLen (organelle) + StringLen (suffix) + 20;
  title = MemNew (sizeof (Char) * len);
  if (title == NULL) return;

  StringCpy (title, str);
  if (partial && tmp == NULL) {
    StringCat (title, ", partial");
  }
  if (organelle != NULL) {
    StringCat (title, " (");
    StringCat (title, organelle);
    StringCat (title, ")");
  }
  if (suffix != NULL) {
    StringCat (title, " [");
    StringCat (title, suffix);
    StringCat (title, "]");
  }
  MemFree (str);
  ttl->data.ptrvalue = title;
}

NLM_EXTERN void CleanUpProteinTitles (SeqEntryPtr sep)

{
  if (sep == NULL) return;
  VisitBioseqsInSep (sep, NULL, AddPartialToProteinTitle);
}

NLM_EXTERN void BasicSeqEntryCleanup (SeqEntryPtr sep)

{
  AuthorPtr       ap;
  ContactInfoPtr  cip;
  CitSubPtr       csp;
  Uint2           entityID;
  Boolean         isEmblOrDdbj = FALSE;
  Boolean         isJscan = FALSE;
  FastNode        labelnode;
  ValNodePtr      pairlist = NULL;
  ValNodePtr      publist = NULL;
  SeqEntryPtr     oldscope;
  ObjMgrDataPtr   omdp;
  SubmitBlockPtr  sbp;
  SeqSubmitPtr    ssp;
  Boolean         stripSerial = TRUE;

  if (sep == NULL) return;

  /* InGpsGenomic needs idx fields assigned */

  entityID = SeqMgrGetEntityIDForSeqEntry (sep);
  AssignIDsInEntityEx (entityID, 0, NULL, NULL);

  /* HandleXrefOnCDS call to GetBestProteinFeatureUnindexed now scoped within record */

  oldscope = SeqEntrySetScope (sep);

  /* clean up spaces in local IDs */

  VisitBioseqsInSep (sep, NULL, CleanSeqIdInBioseq);
  VisitFeaturesInSep (sep, NULL, CleanSeqIdInSeqFeat);
  VisitAlignmentsInSep (sep, NULL, CleanSeqIdInSeqAlign);
  VisitGraphsInSep (sep, NULL, CleanSeqIdInSeqGraph);
  VisitAnnotsInSep (sep, NULL, CleanSeqIdInSeqAnnot);

  /* Fix Bioseq-sets with class 0 */

  VisitSetsInSep (sep, NULL, FixBadSetClass);

  /* removed unnecessarily nested Pub-equivs */

  VisitPubdescsInSep (sep, NULL, FlattenPubdesc);
  VisitFeaturesInSep (sep, NULL, FlattenSfpCit);

  SeqEntryExplore (sep, (Pointer) &stripSerial, CheckForSwissProtID);
  SeqEntryExplore (sep, (Pointer) &isEmblOrDdbj, CheckForEmblDdbjID);
  SeqEntryExplore (sep, (Pointer) &isJscan, CheckForJournalScanID);
#ifdef SUPPRESS_STRIP_SERIAL_DIFFERENCES
  stripSerial = FALSE;
#endif

  BasicSeqEntryCleanupInternal (sep, &publist, isEmblOrDdbj, isJscan, stripSerial);
  if (publist != NULL) {
    VisitFeaturesInSep (sep, (Pointer) publist, ChangeCitsOnFeats);
  }
  ValNodeFreeData (publist);

  /* now get muid/pmid pairs, update sfp->cits to pmids */

  VisitPubdescsInSep (sep, (Pointer) &pairlist, GetMuidPmidPairs);
  if (pairlist != NULL) {
    VisitFeaturesInSep (sep, (Pointer) pairlist, ChangeFeatCitsToPmid);
  }
  ValNodeFree (pairlist);

  labelnode.head = NULL;
  labelnode.tail = NULL;
  VisitPubdescsInSep (sep, (Pointer) &labelnode, GetCitGenLabels);
  if (labelnode.head != NULL) {
    VisitFeaturesInSep (sep, (Pointer) labelnode.head, UpdateShortFeatCits);
  }
  ValNodeFreeData (labelnode.head);

  SeqEntrySetScope (oldscope);

  /* also normalize authors on submit block citation */

  entityID = SeqMgrGetEntityIDForSeqEntry (sep);
  omdp = ObjMgrGetData (entityID);
  if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
    ssp = (SeqSubmitPtr) omdp->dataptr;
    if (ssp != NULL && ssp->datatype == 1) {
      sbp = ssp->sub;
      if (sbp != NULL) {
        csp = sbp->cit;
        if (csp != NULL) {
          NormalizeAuthors (csp->authors, TRUE);
        }
        cip = sbp->contact;
        if (cip != NULL) {
          ap = cip->contact;
          if (ap != NULL) {
            ap->affil = CleanAffil (ap->affil);
          }
        }
      }
    }
  }

  /*
  dynamically add missing partial to already instantiated protein
  titles, in between main title and bracketed organism name
  */

  VisitBioseqsInSep (sep, NULL, AddPartialToProteinTitle);
}

typedef struct bsecsmfedata {
  Int4  max;
  Int4  num_at_max;
} BsecSmfeData, PNTR BsecSmfePtr;

static Boolean LIBCALLBACK BsecSMFEProc (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)


{
  BsecSmfePtr  bsp;
  Int4         len;

  if (sfp == NULL || context == NULL) return TRUE;
  bsp = context->userdata;
  if (bsp == NULL) return TRUE;

  len = SeqLocLen (sfp->location);
  if (len < bsp->max) {
    bsp->max = len;
    bsp->num_at_max = 1;
  } else if (len == bsp->max) {
    (bsp->num_at_max)++;
  }

  return TRUE;
}

NLM_EXTERN void RemoveUnnecessaryGeneXrefs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BsecSmfeData         bsd;
  SeqFeatPtr           cds;
  Int2                 count;
  SeqFeatXrefPtr       curr, next;
  SeqMgrFeatContext    fcontext;
  GeneRefPtr           grp, grpx;
  SeqFeatXrefPtr PNTR  last;
  BioseqPtr            prd;
  SeqFeatPtr           sfpx;
  CharPtr              syn1, syn2;

  if (sfp == NULL || sfp->data.choice == SEQFEAT_GENE) return;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL || SeqMgrGeneIsSuppressed (grp)) return;

  grpx = NULL;
  sfpx = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
  if (sfpx != NULL) {
    if (sfpx->data.choice != SEQFEAT_GENE) return;
    grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
  } else {
    prd = BioseqFindFromSeqLoc (sfp->location);
    if (prd != NULL && ISA_aa (prd->mol)) {
      cds = SeqMgrGetCDSgivenProduct (prd, NULL);
      if (cds != NULL) {
        grpx = SeqMgrGetGeneXref (cds);
        if (grpx == NULL) {
          sfpx = SeqMgrGetOverlappingGene (cds->location, &fcontext);
          if (sfpx != NULL && sfpx->data.choice == SEQFEAT_GENE) {
            grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
          }
        }
      }
    }
  }
  if (grpx == NULL || SeqMgrGeneIsSuppressed (grp)) return;

  if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (grpx->locus_tag)) {
    if (StringICmp (grp->locus_tag, grpx->locus_tag) != 0) return;
  } else if (StringDoesHaveText (grp->locus) && StringDoesHaveText (grpx->locus)) {
    if (StringICmp (grp->locus, grpx->locus) != 0) return;
  } else if (grp->syn != NULL && grpx->syn != NULL) {
    syn1 = (CharPtr) grp->syn->data.ptrvalue;
    syn2 = (CharPtr) grpx->syn->data.ptrvalue;
    if (StringDoesHaveText (syn1) && StringDoesHaveText (syn2)) {
      if (StringICmp (syn1, syn2) != 0) return;
    }
  }

  MemSet ((Pointer) &bsd, 0, sizeof (BsecSmfeData));
  bsd.max = INT4_MAX;
  bsd.num_at_max = 0;
  count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_GENE,
                                           NULL, 0, LOCATION_SUBSET,
                                           (Pointer) &bsd, BsecSMFEProc);

  if (bsd.num_at_max < 2) {
    last = (SeqFeatXrefPtr PNTR) &(sfp->xref);
    curr = sfp->xref;
    while (curr != NULL) {
      next = curr->next;
      if (curr->data.choice == SEQFEAT_GENE) {
        *last = next;
        curr->next = NULL;
        SeqFeatXrefFree (curr);
      } else {
        last = &(curr->next);
      }
      curr = next;
    }
  }
}

static void SortSeqFeatFields (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CdRegionPtr  crp;
  ValNodePtr   psp;

  if (sfp == NULL) return;

  sfp->qual = SortFeatureGBQuals (sfp->qual);

  sfp->qual = SortIllegalGBQuals (sfp->qual);

  sfp->dbxref = ValNodeSort (sfp->dbxref, SortDbxref);

  psp = sfp->cit;
  if (psp != NULL && psp->data.ptrvalue) {
    psp->data.ptrvalue = ValNodeSort ((ValNodePtr) psp->data.ptrvalue, SortCits);
  }

  if (sfp->data.choice == SEQFEAT_CDREGION) {
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if (crp != NULL) {
      crp->code_break = SortCodeBreaks (sfp, crp->code_break);
    }
  }
}

static void SortBioSourceFields (
  BioSourcePtr biop,
  Pointer userdata
)

{
  OrgNamePtr  onp;
  OrgRefPtr   orp;

  if (biop == NULL) return;

  orp = biop->org;
  if (orp != NULL) {
    orp->db = ValNodeSort (orp->db, SortDbxref);

    for (onp = orp->orgname; onp != NULL; onp = onp->next) {
      onp->mod = SortOrgModList (onp->mod);
    }
  }

  biop->subtype = SortSubSourceList (biop->subtype);
}

NLM_EXTERN void SortSeqEntryQualifiers (
  SeqEntryPtr sep
)

{
  if (sep == NULL) return;

  VisitFeaturesInSep (sep, NULL, SortSeqFeatFields);
  VisitBioSourcesInSep (sep, NULL, SortBioSourceFields);
}

/* end BasicSeqEntryCleanup section */

NLM_EXTERN void ResynchCDSPartials (SeqFeatPtr sfp, Pointer userdata)

{
  SeqFeatPtr   bestprot;
  BioseqPtr    bsp;
  MolInfoPtr   mip;
  Boolean      partial5;
  Boolean      partial3;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;
  SeqLocPtr    slp;
  ValNodePtr   vnp;
  /* variables for logging */
  LogInfoPtr    lip;
  CharPtr orig_loc = NULL, new_loc;
  Char    id_buf[100];

  if (sfp->data.choice != SEQFEAT_CDREGION) return;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  sfp->partial = (Boolean) (sfp->partial || partial5 || partial3);
  /*
  slp = SeqLocFindNext (sfp->location, NULL);
  if (slp == NULL) return;
  */
  sip = SeqLocId (sfp->product);
  if (sip == NULL) return;
  bsp = BioseqFind (sip);
  if (bsp == NULL || !ISA_aa (bsp->mol) || bsp->repr != Seq_repr_raw) return;

  bestprot = SeqMgrGetBestProteinFeature (bsp, NULL);
  if (bestprot == NULL) {
    bestprot = GetBestProteinFeatureUnindexed (sfp->product);
  }
  if (bestprot == NULL) return;

  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) return;

  if ((lip = (LogInfoPtr)userdata) != NULL) {
    orig_loc = SeqLocPrintUseBestID (bestprot->location);
  }
  slp = NULL;
  sip = SeqLocId (bestprot->location);
  if (sip != NULL) {
    slp = WholeIntervalFromSeqId (sip);
  }
  if (slp == NULL) {
    slp = CreateWholeInterval (sep);
  }
  if (slp != NULL) {
    bestprot->location = SeqLocFree (bestprot->location);
    bestprot->location = slp;
  }
  SetSeqLocPartial (bestprot->location, partial5, partial3);
  bestprot->partial = sfp->partial;
  if (lip != NULL) {
    new_loc = SeqLocPrintUseBestID (bestprot->location);
    if (StringCmp (orig_loc, new_loc) != 0) {
      lip->data_in_log = TRUE;
      if (lip->fp != NULL) {
        fprintf (lip->fp, "Synchronized coding region partials for protein feature location at %s\n", orig_loc, new_loc);
      }
    }
    new_loc = MemFree (new_loc);
  }
  orig_loc = MemFree (orig_loc);
  vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
  id_buf[0] = 0;
  if (vnp == NULL) {
    vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
    if (vnp != NULL) {
      mip = MolInfoNew ();
      vnp->data.ptrvalue = (Pointer) mip;
      if (mip != NULL) {
        mip->biomol = 8; /* peptide */
        mip->tech = 13; /* concept-trans-author */
        if (lip != NULL) {
          if (lip->fp != NULL) {
            SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_FASTA_SHORT, sizeof (id_buf) - 1);
            fprintf (lip->fp, "Added MolInfo descriptor for %s\n", id_buf);
          }
          lip->data_in_log = TRUE;
        }
      }
    }
  }
  if (vnp != NULL) {
    mip = (MolInfoPtr) vnp->data.ptrvalue;
    if (mip != NULL) {
      if (partial5 && partial3) {
        if (mip->completeness != 5) {
          mip->completeness = 5;
          if (lip != NULL) {
            if (lip->fp != NULL) {
              if (id_buf[0] == 0) {
                SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_FASTA_SHORT, sizeof (id_buf) - 1);
              }
              fprintf (lip->fp, "Adjusted completeness for MolInfo descriptor on %s\n", id_buf);
              lip->data_in_log = TRUE;
            }
          }
        }
      } else if (partial5) {
        if (mip->completeness != 3) {
          mip->completeness = 3;
          if (lip != NULL) {
            if (lip->fp != NULL) {
              if (id_buf[0] == 0) {
                SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_FASTA_SHORT, sizeof (id_buf) - 1);
              }
              fprintf (lip->fp, "Adjusted completeness for MolInfo descriptor on %s\n", id_buf);
            }
            lip->data_in_log = TRUE;
          }
        }
      } else if (partial3) {
        if (mip->completeness != 4) {
          mip->completeness = 4;
          if (lip != NULL) {
            if (lip->fp != NULL) {
              if (id_buf[0] == 0) {
                SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_FASTA_SHORT, sizeof (id_buf) - 1);
              }
              fprintf (lip->fp, "Adjusted completeness for MolInfo descriptor on %s\n", id_buf);
            }
            lip->data_in_log = TRUE;
          }
        }
      } else if (sfp->partial) {
        if (mip->completeness != 2) {
          mip->completeness = 2;
          if (lip != NULL) {
            if (lip->fp != NULL) {
              if (id_buf[0] == 0) {
                SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_FASTA_SHORT, sizeof (id_buf) - 1);
              }
              fprintf (lip->fp, "Adjusted completeness for MolInfo descriptor on %s\n", id_buf);
            }
            lip->data_in_log = TRUE;
          }
        }
      } else {
        if (mip->completeness != 0) {
          mip->completeness = 0;
          if (lip != NULL) {
            if (lip->fp != NULL) {
              if (id_buf[0] == 0) {
                SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_FASTA_SHORT, sizeof (id_buf) - 1);
              }
              fprintf (lip->fp, "Adjusted completeness for MolInfo descriptor on %s\n", id_buf);
            }
            lip->data_in_log = TRUE;
          }
        }
      }
    }
  }
}


NLM_EXTERN Boolean ResynchCodingRegionPartialsEx (SeqEntryPtr sep, FILE *log_fp)

{
  LogInfoData lid;
  MemSet (&lid, 0, sizeof (LogInfoData));
  lid.fp = log_fp;
  VisitFeaturesInSep (sep, &lid, ResynchCDSPartials);
  return lid.data_in_log;
}

NLM_EXTERN void ResynchCodingRegionPartials (SeqEntryPtr sep)

{
  ResynchCodingRegionPartialsEx (sep, NULL);
}


NLM_EXTERN void ResynchMRNAPartials (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr    bsp;
  MolInfoPtr   mip;
  Boolean      partial5;
  Boolean      partial3;
  RnaRefPtr    rrp;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;
  SeqLocPtr    slp;
  ValNodePtr   vnp;

  if (sfp->data.choice != SEQFEAT_RNA) return;
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL || rrp->type != 2) return;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  sfp->partial = (Boolean) (sfp->partial || partial5 || partial3);
  slp = SeqLocFindNext (sfp->location, NULL);
  if (slp == NULL) return;
  sip = SeqLocId (sfp->product);
  if (sip == NULL) return;
  bsp = BioseqFind (sip);
  if (bsp != NULL && ISA_na (bsp->mol) && bsp->repr == Seq_repr_raw) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    if (sep == NULL) return;
    vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
    if (vnp == NULL) {
      vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
      if (vnp != NULL) {
        mip = MolInfoNew ();
        vnp->data.ptrvalue = (Pointer) mip;
        if (mip != NULL) {
          mip->biomol = 3; /* mRNA */
          mip->tech = 1; /* standard */
        }
      }
    }
    if (vnp != NULL) {
      mip = (MolInfoPtr) vnp->data.ptrvalue;
      if (mip != NULL) {
        if (partial5 && partial3) {
          mip->completeness = 5;
        } else if (partial5) {
          mip->completeness = 3;
        } else if (partial3) {
          mip->completeness = 4;
        } else if (sfp->partial) {
          mip->completeness = 2;
        } else {
          mip->completeness = 0;
        }
      }
    }
  }
}

NLM_EXTERN void ResynchMessengerRNAPartials (SeqEntryPtr sep)

{
  VisitFeaturesInSep (sep, NULL, ResynchMRNAPartials);
}

NLM_EXTERN void ResynchPeptidePartials (SeqFeatPtr sfp, Pointer userdata)

{
  SeqFeatPtr   bestprot;
  BioseqPtr    bsp;
  MolInfoPtr   mip;
  Boolean      partial5;
  Boolean      partial3;
  ProtRefPtr   prp;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;
  SeqLocPtr    slp;
  ValNodePtr   vnp;

  if (sfp->data.choice != SEQFEAT_PROT) return;
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL) return;
  if (prp->processed < 1 || prp->processed > 4) return;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  sfp->partial = (Boolean) (sfp->partial || partial5 || partial3);
  /*
  slp = SeqLocFindNext (sfp->location, NULL);
  if (slp == NULL) return;
  */
  sip = SeqLocId (sfp->product);
  if (sip == NULL) return;
  bsp = BioseqFind (sip);
  if (bsp != NULL && ISA_aa (bsp->mol) && bsp->repr == Seq_repr_raw) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    if (sep == NULL) return;
    bestprot = SeqMgrGetBestProteinFeature (bsp, NULL);
    if (bestprot == NULL) {
      bestprot = GetBestProteinFeatureUnindexed (sfp->product);
    }
    if (bestprot != NULL) {
      slp = NULL;
      sip = SeqLocId (bestprot->location);
      if (sip != NULL) {
        slp = WholeIntervalFromSeqId (sip);
      }
      if (slp == NULL) {
        slp = CreateWholeInterval (sep);
      }
      if (slp != NULL) {
        bestprot->location = SeqLocFree (bestprot->location);
        bestprot->location = slp;
      }
      SetSeqLocPartial (bestprot->location, partial5, partial3);
      bestprot->partial = sfp->partial;
    }
    vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
    if (vnp == NULL) {
      vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
      if (vnp != NULL) {
        mip = MolInfoNew ();
        vnp->data.ptrvalue = (Pointer) mip;
        if (mip != NULL) {
          mip->biomol = 8;
          mip->tech = 13;
        }
      }
    }
    if (vnp != NULL) {
      mip = (MolInfoPtr) vnp->data.ptrvalue;
      if (mip != NULL) {
        if (partial5 && partial3) {
          mip->completeness = 5;
        } else if (partial5) {
          mip->completeness = 3;
        } else if (partial3) {
          mip->completeness = 4;
        } else if (sfp->partial) {
          mip->completeness = 2;
        } else {
          mip->completeness = 0;
        }
      }
    }
  }
}

NLM_EXTERN void ResynchProteinPartials (SeqEntryPtr sep)

{
  VisitFeaturesInSep (sep, NULL, ResynchPeptidePartials);
}

/* SeqIdStripLocus removes the SeqId.name field if accession is set */

NLM_EXTERN SeqIdPtr SeqIdStripLocus (SeqIdPtr sip)

{
  TextSeqIdPtr  tip;

  if (sip != NULL) {
    switch (sip->choice) {
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
      case SEQID_OTHER :
      case SEQID_TPG:
      case SEQID_TPE:
      case SEQID_TPD:
      case SEQID_GPIPE:
        tip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tip != NULL) {
          if (! HasNoText (tip->accession)) {
            tip->name = MemFree (tip->name);
          }
        }
        break;
      default :
        break;
    }
  }
  return sip;
}

NLM_EXTERN SeqLocPtr StripLocusFromSeqLoc (SeqLocPtr location)

{
  SeqLocPtr      loc;
  SeqLocPtr      next;
  PackSeqPntPtr  psp;
  SeqBondPtr     sbp;
  SeqIntPtr      sinp;
  SeqIdPtr       sip;
  SeqLocPtr      slp;
  SeqPntPtr      spp;

  if (location == NULL) return NULL;
  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL) {
    next = SeqLocFindNext (location, slp);
    switch (slp->choice) {
      case SEQLOC_EMPTY :
      case SEQLOC_WHOLE :
        sip = (SeqIdPtr) slp->data.ptrvalue;
        SeqIdStripLocus (sip);
        break;
      case SEQLOC_INT :
        sinp = (SeqIntPtr) slp->data.ptrvalue;
        if (sinp != NULL) {
          SeqIdStripLocus (sinp->id);
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        loc = (SeqLocPtr) slp->data.ptrvalue;
        while (loc != NULL) {
          sip = SeqLocId (loc);
          SeqIdStripLocus (sip);
          loc = loc->next;
        }
        break;
      case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          spp = sbp->a;
          if (spp != NULL) {
            SeqIdStripLocus (spp->id);
          }
          spp = sbp->b;
          if (spp != NULL) {
            SeqIdStripLocus (spp->id);
          }
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL) {
          SeqIdStripLocus (spp->id);
        }
        break;
      case SEQLOC_PACKED_PNT :
        psp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (psp != NULL) {
          SeqIdStripLocus (psp->id);
        }
        break;
      default :
        break;
    }
    slp = next;
  }
  return location;
}

static void GetRidOfLocusCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  sap = NULL;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1 && sap->data != NULL) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        StripLocusFromSeqLoc (sfp->location);
        StripLocusFromSeqLoc (sfp->product);
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
}

NLM_EXTERN void GetRidOfLocusInSeqIds (Uint2 entityID, SeqEntryPtr sep)

{
  if (entityID < 1 && sep == NULL) return;
  if (entityID > 0 && sep == NULL) {
    sep = GetTopSeqEntryForEntityID (entityID);
  }
  if (sep == NULL) return;
  SeqEntryExplore (sep, NULL, GetRidOfLocusCallback);
}

/* Mac can now use static parse tables by using
   Make Strings Read-Only and Store Static Data in TOC
#ifdef OS_MAC
#define ASNLOAD_NEEDED 1
#endif
*/
#if defined(OS_DOS) || defined(WIN16)
#define ASNLOAD_NEEDED 1
#endif

static Boolean FileExists (CharPtr dirname, CharPtr subname, CharPtr filename)

{
  Char  path [PATH_MAX];

  StringNCpy_0 (path, dirname, sizeof (path));
  FileBuildPath (path, subname, NULL);
  FileBuildPath (path, NULL, filename);
  return (Boolean) (FileLength (path) > 0);
}

/*
static Boolean CheckAsnloadPath (CharPtr dirname, CharPtr subdir)

{
#ifdef ASNLOAD_NEEDED
  Char  fname [16];
  int   i;

  for (i = 60; i <= 99; ++i) {
    sprintf (fname, "asnmedli.l%02d", (int) i);
    if (FileExists (dirname, subdir, fname)) {
      return TRUE;
    }
  }
  return FALSE;
#else
  return TRUE;
#endif
}
*/

static Boolean CheckDataPath (CharPtr dirname, CharPtr subdir)

{
  if (FileExists (dirname, subdir, "seqcode.val")) return TRUE;
  return (Boolean) (FileExists (dirname, subdir, "objprt.prt"));
}

static Boolean CheckErrMsgPath (CharPtr dirname, CharPtr subdir)

{
  return (Boolean) (FileExists (dirname, subdir, "valid.msg"));
}

static void SetTransientPath (CharPtr dirname, CharPtr subname, CharPtr file,
                              CharPtr section, CharPtr type)

{
  Char  path [PATH_MAX];

  StringNCpy_0 (path, dirname, sizeof (path));
  FileBuildPath (path, subname, NULL);
  TransientSetAppParam (file, section, type, path);
}

NLM_EXTERN Boolean UseLocalAsnloadDataAndErrMsg (void)

{
  Boolean  dataFound;
  Char     path [PATH_MAX];
  Char     appPath[PATH_MAX];
  CharPtr  ptr;

  ProgramPath (appPath, sizeof (appPath));
  StrCpy(path, appPath);
  /* data a sibling of our application? */
  ptr = StringRChr (path, DIRDELIMCHR);
  if (ptr != NULL) {
    ptr++;
    *ptr = '\0';
  }
  dataFound = CheckDataPath (path, "data");
  if (! (dataFound)) {
  /* data an uncle of our application? */
    if (ptr != NULL) {
      ptr--;
      *ptr = '\0';
      ptr = StringRChr (path, DIRDELIMCHR);
      if (ptr != NULL) {
        ptr++;
        *ptr = '\0';
      }
      dataFound = CheckDataPath (path, "data");
    }
  }
#ifdef OS_UNIX_DARWIN
  if (! (dataFound) && IsApplicationPackage (appPath)) {
      /* is data inside our application within Contents/Resources? */
      StrCpy (path, appPath);
      FileBuildPath (path, "Contents", NULL);
      FileBuildPath (path, "Resources", NULL);
      dataFound = CheckDataPath (path, "data");
      if (! dataFound) {
        StrCpy (path, appPath);
        ptr = StringStr (path, "/ncbi/build/");
        if (ptr != NULL) {
          /* see if running under older Xcode 3 build environment */
          ptr [5] = '\0';
          dataFound = CheckDataPath (path, "data");
        }
      }
      if (! dataFound) {
        StrCpy (path, appPath);
        ptr = StringStr (path, "/ncbi/make/");
        if (ptr != NULL) {
          /* see if running under newer Xcode 3 build environment */
          ptr [5] = '\0';
          dataFound = CheckDataPath (path, "data");
        }
      }
      if (! dataFound) {
          StrCpy (path, appPath);
          ptr = StringStr (path, "/Library/Developer/");
          if (ptr != NULL) {
              /* see if running under Xcode 4 build environment */
              ptr [19] = '\0';
              dataFound = CheckDataPath (path, "data");
          }
      }
  }
#endif
  if (dataFound) {
    SetTransientPath (path, "asnload", "NCBI", "NCBI", "ASNLOAD");
    SetTransientPath (path, "data", "NCBI", "NCBI", "DATA");
    if (CheckErrMsgPath (path, "errmsg")) {
      SetTransientPath (path, "errmsg", "NCBI", "ErrorProcessing", "MsgPath");
      TransientSetAppParam ("NCBI", "ErrorProcessing", "EO_BEEP", "No");
    }
    return TRUE;
  }
  return FALSE;
}

typedef struct miscdata {
  SeqEntryPtr  sep;
  Int2         count;
  Int2         desired;
  Uint1        _class;
} MiscData, PNTR MiscDataPtr;

static void FindNthSeqEntryCallback (SeqEntryPtr sep, Pointer mydata,
                                     Int4 index, Int2 indent)

{
  MiscDataPtr  mdp;

  if (sep != NULL && mydata != NULL) {
    mdp = (MiscDataPtr) mydata;
    (mdp->count)++;
    if (mdp->count == mdp->desired) {
      mdp->sep = sep;
    }
  }
}

NLM_EXTERN SeqEntryPtr LIBCALL FindNthSeqEntry (SeqEntryPtr sep, Int2 seq)

{
  MiscData  md;

  md.sep = NULL;
  md.count = 0;
  md.desired = seq;
  if (sep != NULL) {
    SeqEntryExplore (sep, (Pointer) (&md), FindNthSeqEntryCallback);
  }
  return md.sep;
}

NLM_EXTERN SeqEntryPtr LIBCALL FindNthBioseq (SeqEntryPtr sep, Int2 seq)

{
  MiscData  md;

  md.sep = NULL;
  md.count = 0;
  md.desired = seq;
  if (sep != NULL) {
    BioseqExplore (sep, (Pointer) (&md), FindNthSeqEntryCallback);
  }
  return md.sep;
}

NLM_EXTERN SeqEntryPtr LIBCALL FindNthSequinEntry (SeqEntryPtr sep, Int2 seq)

{
  MiscData  md;

  md.sep = NULL;
  md.count = 0;
  md.desired = seq;
  if (sep != NULL) {
    SequinEntryExplore (sep, (Pointer) (&md), FindNthSeqEntryCallback);
  }
  return md.sep;
}

static void FindNucSeqEntryCallback (SeqEntryPtr sep, Pointer mydata,
                                     Int4 index, Int2 indent)

{
  BioseqPtr    bsp;
  MiscDataPtr  mdp;

  if (sep != NULL && sep->choice == 1 && mydata != NULL) {
    mdp = (MiscDataPtr) mydata;
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp != NULL && ISA_na (bsp->mol)) {
      if (mdp->sep == NULL) {
        mdp->sep = sep;
      }
    }
  }
}

NLM_EXTERN SeqEntryPtr LIBCALL FindNucSeqEntry (SeqEntryPtr sep)

{
  MiscData  md;

  md.sep = NULL;
  md.count = 0;
  md.desired = 0;
  if (sep != NULL) {
    BioseqExplore (sep, (Pointer) (&md), FindNucSeqEntryCallback);
  }
  return md.sep;
}

NLM_EXTERN BioseqPtr LIBCALL FindNucBioseq (SeqEntryPtr sep)

{
  BioseqPtr    nbsp;
  SeqEntryPtr  nsep;

  nsep = FindNucSeqEntry (sep);
  if (nsep == NULL) return NULL;
  if (! IS_Bioseq (nsep)) return NULL;
  nbsp = (BioseqPtr) nsep->data.ptrvalue;
  return nbsp;
}

static void FindBioseqSetByClassCallback (SeqEntryPtr sep, Pointer mydata,
                                          Int4 index, Int2 indent)

{
  BioseqSetPtr  bssp;
  MiscDataPtr   mdp;

  if (sep != NULL && sep->choice == 2 && mydata != NULL) {
    mdp = (MiscDataPtr) mydata;
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && bssp->_class == mdp->_class) {
      if (mdp->sep == NULL) {
        mdp->sep = sep;
      }
    }
  }
}

NLM_EXTERN SeqEntryPtr LIBCALL FindBioseqSetByClass (SeqEntryPtr sep, Uint1 _class)

{
  MiscData  md;

  md.sep = NULL;
  md.count = 0;
  md.desired = 0;
  md._class = _class;
  if (sep != NULL) {
    SeqEntryExplore (sep, (Pointer) (&md), FindBioseqSetByClassCallback);
  }
  return md.sep;
}

typedef struct kinddata {
  Boolean  hasNuc;
  Boolean  hasProt;
} KindData, PNTR KindPtr;

static void HasNucOrProtCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr  bsp;
  KindPtr    kptr;

  if (sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL && mydata != NULL) {
    kptr = (KindPtr) mydata;
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (ISA_na (bsp->mol)) {
      kptr->hasNuc = TRUE;
    } else if (ISA_aa (bsp->mol)) {
      kptr->hasProt = TRUE;
    }
  }
}

NLM_EXTERN Boolean LIBCALL SeqEntryHasNucs (SeqEntryPtr sep)

{
  KindData  kd;

  kd.hasNuc = FALSE;
  kd.hasProt = FALSE;
  if (sep != NULL) {
    BioseqExplore (sep, (Pointer) (&kd), HasNucOrProtCallback);
  }
  return kd.hasNuc;
}

NLM_EXTERN Boolean LIBCALL SeqEntryHasProts (SeqEntryPtr sep)

{
  KindData  kd;

  kd.hasNuc = FALSE;
  kd.hasProt = FALSE;
  if (sep != NULL) {
    BioseqExplore (sep, (Pointer) (&kd), HasNucOrProtCallback);
  }
  return kd.hasProt;
}

static Boolean CheckForAlignments (GatherContextPtr gcp)

{
  BoolPtr  boolptr;

  if (gcp == NULL) return TRUE;

  boolptr = (BoolPtr) gcp->userdata;
  if (boolptr == NULL ) return TRUE;

  switch (gcp->thistype) {
    case OBJ_SEQALIGN :
    case OBJ_SEQHIST_ALIGN :
      *boolptr = TRUE;
      return TRUE;
    default :
      break;
  }
  return TRUE;
}

NLM_EXTERN Boolean LIBCALL SeqEntryHasAligns (Uint2 entityID, SeqEntryPtr sep)

{
  GatherScope  gs;
  Boolean      rsult;

  rsult = FALSE;
  if (entityID == 0 || sep == NULL) return FALSE;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  MemSet((Pointer) (gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQALIGN] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQHIST] = FALSE;
  gs.ignore[OBJ_SEQHIST_ALIGN] = FALSE;
  gs.scope = sep;
  GatherEntity (entityID, (Pointer) (&rsult), CheckForAlignments, &gs);
  return rsult;
}

static void FindPowerBLASTAsnCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  AnnotDescrPtr  desc;
  ObjectIdPtr    oip;
  SeqAnnotPtr    sap;
  BoolPtr        rsult;

  if (sep == NULL || sep->data.ptrvalue == NULL || mydata == NULL) return;
  rsult = (BoolPtr) mydata;
  sap = (IS_Bioseq (sep)) ?
         ((BioseqPtr) sep->data.ptrvalue)->annot :
         ((BioseqSetPtr) sep->data.ptrvalue)->annot;
  while (sap != NULL) {
    if (sap->type == 2) {
      desc = NULL;
      while ((desc = ValNodeFindNext (sap->desc, desc, Annot_descr_user)) != NULL) {
        if (desc->data.ptrvalue != NULL) {
          oip = ((UserObjectPtr) desc->data.ptrvalue)->type;
          if (oip != NULL && StringCmp (oip->str, "Hist Seqalign") == 0) {
            *rsult = TRUE;
          }
        }
      }
    }
    sap = sap->next;
  }
}

NLM_EXTERN Boolean LIBCALL PowerBLASTASN1Detected (SeqEntryPtr sep)

{
  Boolean  rsult;

  rsult = FALSE;
  SeqEntryExplore (sep, (Pointer) &rsult, FindPowerBLASTAsnCallback);
  return rsult;
}

NLM_EXTERN SeqLocPtr CreateWholeInterval (SeqEntryPtr sep)

{
  BioseqPtr  bsp;
  SeqIntPtr  sip;
  SeqLocPtr  slp;

  slp = NULL;
  if (sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    slp = ValNodeNew (NULL);
    if (slp != NULL) {
      sip = SeqIntNew ();
      if (sip != NULL) {
        slp->choice = SEQLOC_INT;
        slp->data.ptrvalue = (Pointer) sip;
        sip->from = 0;
        sip->to = bsp->length - 1;
        if (ISA_na (bsp->mol)) {
          sip->strand = Seq_strand_plus;
        }
        sip->id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
      }
    }
  }
  return slp;
}

NLM_EXTERN SeqLocPtr WholeIntervalFromSeqId (SeqIdPtr sip)

{
  BioseqPtr  bsp;
  SeqIntPtr  sintp;
  SeqLocPtr  slp;

  if (sip == NULL) return NULL;
  bsp = BioseqFindCore (sip);
  if (bsp == NULL) return NULL;
  slp = ValNodeNew (NULL);
  if (slp == NULL) return NULL;
  sintp = SeqIntNew ();
  if (sintp == NULL) return NULL;
  slp->choice = SEQLOC_INT;
  slp->data.ptrvalue = (Pointer) sintp;
  sintp->from = 0;
  sintp->to = bsp->length - 1;
  if (ISA_na (bsp->mol)) {
    sintp->strand = Seq_strand_plus;
  }
  sintp->id = SeqIdStripLocus (SeqIdDup (sip));
  return slp;
}

NLM_EXTERN void FreeAllFuzz (SeqLocPtr location)

{
  SeqIntPtr  sip;
  SeqLocPtr  slp;

  if (location == NULL) return;
  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL) {
    if (slp->choice == SEQLOC_INT) {
      sip = (SeqIntPtr) slp->data.ptrvalue;
      if (sip != NULL) {
        sip->if_to = IntFuzzFree (sip->if_to);
        sip->if_from = IntFuzzFree (sip->if_from);
      }
    }
    slp = SeqLocFindNext (location, slp);
  }
}

NLM_EXTERN Boolean LocationHasNullsBetween (SeqLocPtr location)

{
  SeqLocPtr  slp;

  if (location == NULL) return FALSE;
  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL) {
    if (slp->choice == SEQLOC_NULL) return TRUE;
    slp = SeqLocFindNext (location, slp);
  }
  return FALSE;
}

NLM_EXTERN void NormalizeNullsBetween (SeqLocPtr location)

{
  SeqLocPtr  next, tmp, vnp;

  if (location == NULL) return;
  if (! LocationHasNullsBetween (location)) return;

  if (location->choice != SEQLOC_MIX) return;
  vnp = (ValNodePtr) location->data.ptrvalue;
  if (vnp == NULL) return;

  while (vnp != NULL && vnp->next != NULL) {
    next = vnp->next;
    if (vnp->choice != SEQLOC_NULL && next->choice != SEQLOC_NULL) {
      tmp = ValNodeNew (NULL);
      if (tmp != NULL) {
        tmp->choice = SEQLOC_NULL;
        tmp->next = vnp->next;
        vnp->next = tmp;
      }
    }
    vnp = next;
  }
}

NLM_EXTERN Uint1 FindFeatFromFeatDefType (Uint2 subtype)

{
  switch (subtype) {
    case FEATDEF_GENE :
      return SEQFEAT_GENE;
    case FEATDEF_ORG :
      return SEQFEAT_ORG;
    case FEATDEF_CDS :
      return SEQFEAT_CDREGION;
    case FEATDEF_PROT :
      return SEQFEAT_PROT;
    case FEATDEF_PUB :
      return SEQFEAT_PUB;
    case FEATDEF_SEQ :
      return SEQFEAT_SEQ;
    case FEATDEF_REGION :
      return SEQFEAT_REGION;
    case FEATDEF_COMMENT :
      return SEQFEAT_COMMENT;
    case FEATDEF_BOND :
      return SEQFEAT_BOND;
    case FEATDEF_SITE :
      return SEQFEAT_SITE;
    case FEATDEF_RSITE :
      return SEQFEAT_RSITE;
    case FEATDEF_USER :
      return SEQFEAT_USER;
    case FEATDEF_TXINIT :
      return SEQFEAT_TXINIT;
    case FEATDEF_NUM :
      return SEQFEAT_NUM;
    case FEATDEF_PSEC_STR :
      return SEQFEAT_PSEC_STR;
    case FEATDEF_NON_STD_RESIDUE :
      return SEQFEAT_NON_STD_RESIDUE;
    case FEATDEF_HET :
      return SEQFEAT_HET;
    case FEATDEF_BIOSRC :
      return SEQFEAT_BIOSRC;
    default :
      if (subtype >= FEATDEF_preRNA && subtype <= FEATDEF_otherRNA) {
        return SEQFEAT_RNA;
      }
      if (subtype == FEATDEF_snoRNA) {
        return SEQFEAT_RNA;
      }
      if (subtype >= FEATDEF_ncRNA && subtype <= FEATDEF_tmRNA) {
        return SEQFEAT_RNA;
      }
      if (subtype >= FEATDEF_preprotein && subtype <= FEATDEF_transit_peptide_aa) {
        return SEQFEAT_PROT;
      }
      if (subtype >= FEATDEF_IMP && subtype <= FEATDEF_site_ref) {
        return SEQFEAT_IMP;
      }
      if (subtype >= FEATDEF_gap && subtype <= FEATDEF_oriT) {
        return SEQFEAT_IMP;
      }
      if (subtype == FEATDEF_mobile_element) {
        return SEQFEAT_IMP;
      }
  }
  return 0;
}

NLM_EXTERN SeqIdPtr MakeSeqID (CharPtr str)

{
  CharPtr   buf;
  Int4      len;
  SeqIdPtr  sip;

  sip = NULL;
  if (str != NULL && *str != '\0') {
    if (StringChr (str, '|') != NULL) {
      sip = SeqIdParse (str);
    } else {
      len = StringLen (str) + 5;
      buf = (CharPtr) MemNew (sizeof (Char) * len);
      sprintf (buf, "lcl|%s", str);
      sip = SeqIdParse (buf);
      buf = MemFree (buf);
    }
  }
  return sip;
}

NLM_EXTERN SeqIdPtr MakeUniqueSeqID (CharPtr prefix)

{
    Char buf[40];
    CharPtr tmp;
    Int2 ctr;
    ValNodePtr newid;
    ObjectIdPtr oid;
    ValNode vn;
    TextSeqId tsi;
    ValNodePtr altid;
    size_t len;

    altid = &vn;
    vn.choice = SEQID_GENBANK;
    vn.next = NULL;
    vn.data.ptrvalue = &tsi;
    tsi.name = NULL;
    tsi.accession = NULL;
    tsi.release = NULL;
    tsi.version = INT2_MIN;

    len = StringLen (prefix);
    if (len > 0 && len < 32) {
        tmp = StringMove(buf, prefix);
    } else {
        tmp = StringMove(buf, "tmpseq_");
    }

    newid = ValNodeNew(NULL);
    oid = ObjectIdNew();
    oid->str = buf;   /* allocate this later */
    newid->choice = SEQID_LOCAL;
    newid->data.ptrvalue = oid;

    tsi.name = buf;   /* check for alternative form */

    for (ctr = 1; ctr < 32000; ctr++)
    {
        sprintf(tmp, "%d", (int)ctr);
        if ((BioseqFindCore(newid) == NULL) && (BioseqFindCore(altid) == NULL))
        {
            oid->str = StringSave(buf);
            return newid;
        }
    }

    return NULL;
}

NLM_EXTERN SeqIdPtr SeqIdFindWorst (SeqIdPtr sip)

{
  Uint1  order [NUM_SEQID];

  SeqIdBestRank (order, NUM_SEQID);
  order [SEQID_LOCAL] = 10;
  order [SEQID_GENBANK] = 5;
  order [SEQID_EMBL] = 5;
  order [SEQID_PIR] = 5;
  order [SEQID_SWISSPROT] = 5;
  order [SEQID_DDBJ] = 5;
  order [SEQID_PRF] = 5;
  order [SEQID_PDB] = 5;
  order [SEQID_TPG] = 5;
  order [SEQID_TPE] = 5;
  order [SEQID_TPD] = 5;
  order [SEQID_GPIPE] = 9;
  order [SEQID_NAMED_ANNOT_TRACK] = 9;
  order [SEQID_PATENT] = 10;
  order [SEQID_OTHER] = 8;
  order [SEQID_GENERAL] = 15;
  order [SEQID_GIBBSQ] = 15;
  order [SEQID_GIBBMT] = 15;
  order [SEQID_GIIM] = 20;
  order [SEQID_GI] = 20;
  return SeqIdSelect (sip, order, NUM_SEQID);
}

NLM_EXTERN SeqFeatPtr CreateNewFeature (SeqEntryPtr sep, SeqEntryPtr placeHere,
                             Uint1 choice, SeqFeatPtr useThis)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqFeatPtr    prev;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;

  if (sep == NULL || sep->choice != 1) return NULL;
  sfp = NULL;
  bsp = NULL;
  bssp = NULL;
  if (placeHere == NULL) {
    placeHere = sep;
  }
  if (placeHere != NULL && placeHere->data.ptrvalue != NULL) {
    if (placeHere->choice == 1) {
      bsp = (BioseqPtr) placeHere->data.ptrvalue;
      sap = bsp->annot;
      while (sap != NULL && (sap->name != NULL || sap->desc != NULL || sap->type != 1)) {
        sap = sap->next;
      }
      if (sap == NULL) {
        sap = SeqAnnotNew ();
        if (sap != NULL) {
          sap->type = 1;
          sap->next = bsp->annot;
          bsp->annot = sap;
        }
        sap = bsp->annot;
      }
    } else if (placeHere->choice == 2) {
      bssp = (BioseqSetPtr) placeHere->data.ptrvalue;
      sap = bssp->annot;
      while (sap != NULL && (sap->name != NULL || sap->desc != NULL || sap->type != 1)) {
        sap = sap->next;
      }
      if (sap == NULL) {
        sap = SeqAnnotNew ();
        if (sap != NULL) {
          sap->type = 1;
          sap->next = bssp->annot;
          bssp->annot = sap;
        }
        sap = bssp->annot;
      }
    } else {
      return NULL;
    }
    if (sap != NULL) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (useThis != NULL) {
        sfp = useThis;
      } else {
        sfp = SeqFeatNew ();
      }
      if (sap->data != NULL) {
        prev = sap->data;
        while (prev->next != NULL) {
          prev = prev->next;
        }
        prev->next = sfp;
      } else {
        sap->data = (Pointer) sfp;
      }
      if (sfp != NULL) {
        sfp->data.choice = choice;
        if (useThis == NULL) {
          sfp->location = CreateWholeInterval (sep);
        }
      }
    }
  }
  return sfp;
}

NLM_EXTERN SeqFeatPtr CreateNewFeatureOnBioseq (BioseqPtr bsp, Uint1 choice, SeqLocPtr slp)

{
  SeqEntryPtr  sep;
  SeqFeatPtr   sfp;

  if (bsp == NULL) return NULL;
  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) return NULL;
  sfp = CreateNewFeature (sep, NULL, choice, NULL);
  if (sfp == NULL) return NULL;
  if (slp != NULL) {
    sfp->location = SeqLocFree (sfp->location);
    sfp->location = AsnIoMemCopy (slp, (AsnReadFunc) SeqLocAsnRead,
                                  (AsnWriteFunc) SeqLocAsnWrite);
  }
  return sfp;
}

NLM_EXTERN ValNodePtr CreateNewDescriptor (SeqEntryPtr sep, Uint1 choice)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Uint1         _class;
  ValNodePtr    descr;
  SeqEntryPtr   seqentry;
  ValNodePtr    vnp;

  vnp = NULL;
  if (sep != NULL) {
    descr = NULL;
    vnp = NULL;
    bsp = NULL;
    bssp = NULL;
    seqentry = sep;
    while (seqentry != NULL) {
      if (seqentry->choice == 1) {
        bsp = (BioseqPtr) seqentry->data.ptrvalue;
        if (bsp != NULL) {
          descr = bsp->descr;
          vnp = SeqDescrNew (descr);
          if (descr == NULL) {
            bsp->descr = vnp;
          }
        }
        seqentry = NULL;
      } else if (seqentry->choice == 2) {
        bssp = (BioseqSetPtr) seqentry->data.ptrvalue;
        if (bssp != NULL) {
          _class = bssp->_class;
          if (_class == 7) {
            descr = bssp->descr;
            vnp = SeqDescrNew (descr);
            if (descr == NULL) {
              bssp->descr = vnp;
            }
            seqentry = NULL;
          } else if ((_class >= 5 && _class <= 8) || _class == 11 /* || _class == BioseqseqSet_class_gen_prod_set */) {
            seqentry = bssp->seq_set;
          } else {
            descr = bssp->descr;
            vnp = SeqDescrNew (descr);
            if (descr == NULL) {
              bssp->descr = vnp;
            }
            seqentry = NULL;
          }
        } else {
          seqentry = NULL;
        }
      } else {
        seqentry = NULL;
      }
    }
    if (vnp != NULL) {
      vnp->choice = choice;
    }
  }
  return vnp;
}

NLM_EXTERN ValNodePtr CreateNewDescriptorOnBioseq (BioseqPtr bsp, Uint1 choice)

{
  SeqEntryPtr  sep;

  if (bsp == NULL) return NULL;
  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) return NULL;
  return CreateNewDescriptor (sep, choice);
}

/* common functions to scan binary ASN.1 file of entire release as Bioseq-set */

static Int4 VisitSeqIdList (SeqIdPtr sip, Pointer userdata, VisitSeqIdFunc callback)

{
  Int4  index = 0;

  while (sip != NULL) {
    if (callback != NULL) {
      callback (sip, userdata);
    }
    index++;
    sip = sip->next;
  }
  return index;
}

NLM_EXTERN Int4 VisitSeqIdsInSeqLoc (SeqLocPtr slp, Pointer userdata, VisitSeqIdFunc callback)

{
  Int4           index = 0;
  SeqLocPtr      loc;
  PackSeqPntPtr  psp;
  SeqBondPtr     sbp;
  SeqIntPtr      sinp;
  SeqIdPtr       sip;
  SeqPntPtr      spp;

  if (slp == NULL) return index;

  while (slp != NULL) {
    switch (slp->choice) {
      case SEQLOC_NULL :
        break;
      case SEQLOC_EMPTY :
      case SEQLOC_WHOLE :
        sip = (SeqIdPtr) slp->data.ptrvalue;
        index += VisitSeqIdList (sip, userdata, callback);
        break;
      case SEQLOC_INT :
        sinp = (SeqIntPtr) slp->data.ptrvalue;
        if (sinp != NULL) {
          sip = sinp->id;
          index += VisitSeqIdList (sip, userdata, callback);
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL) {
          sip = spp->id;
          index += VisitSeqIdList (sip, userdata, callback);
        }
        break;
      case SEQLOC_PACKED_PNT :
        psp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (psp != NULL) {
          sip = psp->id;
          index += VisitSeqIdList (sip, userdata, callback);
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        loc = (SeqLocPtr) slp->data.ptrvalue;
        while (loc != NULL) {
          index += VisitSeqIdsInSeqLoc (loc, userdata, callback);
          loc = loc->next;
        }
        break;
      case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          spp = (SeqPntPtr) sbp->a;
          if (spp != NULL) {
            sip = spp->id;
            index += VisitSeqIdList (sip, userdata, callback);
          }
          spp = (SeqPntPtr) sbp->b;
          if (spp != NULL) {
            sip = spp->id;
            index += VisitSeqIdList (sip, userdata, callback);
          }
        }
        break;
      case SEQLOC_FEAT :
        break;
      default :
        break;
    }
    slp = slp->next;
  }

  return index;
}

NLM_EXTERN Int4 VisitSeqIdsInBioseq (BioseqPtr bsp, Pointer userdata, VisitSeqIdFunc callback)

{
  Int4  index = 0;

  if (bsp == NULL) return index;

  if (bsp->id != NULL) {
    index += VisitSeqIdList (bsp->id, userdata, callback);
  }

  return index;
}

NLM_EXTERN Int4 VisitSeqIdsInSeqFeat (SeqFeatPtr sfp, Pointer userdata, VisitSeqIdFunc callback)

{
  CodeBreakPtr  cbp;
  CdRegionPtr   crp;
  Int4          index = 0;
  RnaRefPtr     rrp;
  tRNAPtr       trp;

  if (sfp == NULL) return index;

  index += VisitSeqIdsInSeqLoc (sfp->location, userdata, callback);
  index += VisitSeqIdsInSeqLoc (sfp->product, userdata, callback);

  switch (sfp->data.choice) {
    case SEQFEAT_CDREGION :
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL) {
        for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
          index += VisitSeqIdsInSeqLoc (cbp->loc, userdata, callback);
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->ext.choice == 2) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        if (trp != NULL && trp->anticodon != NULL) {
          index += VisitSeqIdsInSeqLoc (trp->anticodon, userdata, callback);
        }
      }
      break;
    default :
      break;
  }

  return index;
}

NLM_EXTERN Int4 VisitSeqIdsInSeqAlign (SeqAlignPtr sap, Pointer userdata, VisitSeqIdFunc callback)

{
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  Int4          index = 0;
  SeqIdPtr      sip;
  SeqLocPtr     slp = NULL;
  StdSegPtr     ssp;

  if (sap == NULL) return index;

  if (sap->bounds != NULL) {
    sip = SeqLocId (sap->bounds);
    index += VisitSeqIdList (sip, userdata, callback);
  }

  if (sap->segs == NULL) return index;

  switch (sap->segtype) {
    case SAS_DENDIAG :
      ddp = (DenseDiagPtr) sap->segs;
      if (ddp != NULL) {
        for (sip = ddp->id; sip != NULL; sip = sip->next) {
          index += VisitSeqIdList (sip, userdata, callback);
        }
      }
      break;
    case SAS_DENSEG :
      dsp = (DenseSegPtr) sap->segs;
      if (dsp != NULL) {
        for (sip = dsp->ids; sip != NULL; sip = sip->next) {
          index += VisitSeqIdList (sip, userdata, callback);
        }
      }
      break;
    case SAS_STD :
      ssp = (StdSegPtr) sap->segs;
      for (slp = ssp->loc; slp != NULL; slp = slp->next) {
        sip = SeqLocId (slp);
        index += VisitSeqIdList (sip, userdata, callback);
      }
      break;
    case SAS_DISC :
      /* recursive */
      for (sap = (SeqAlignPtr) sap->segs; sap != NULL; sap = sap->next) {
        index += VisitSeqIdsInSeqAlign (sap, userdata, callback);
      }
      break;
    default :
      break;
  }

  return index;
}

NLM_EXTERN Int4 VisitSeqIdsInSeqGraph (SeqGraphPtr sgp, Pointer userdata, VisitSeqIdFunc callback)

{
  Int4      index = 0;
  SeqIdPtr  sip;

  if (sgp == NULL) return index;

  if (sgp->loc != NULL) {
    sip = SeqLocId (sgp->loc);
    index += VisitSeqIdList (sip, userdata, callback);
  }

  return index;
}

NLM_EXTERN Int4 VisitSeqIdsInSeqAnnot (SeqAnnotPtr annot, Pointer userdata, VisitSeqIdFunc callback)

{
  Int4         index = 0;
  SeqAlignPtr  sap;
  SeqFeatPtr   sfp;
  SeqGraphPtr  sgp;

  if (annot == NULL || annot->data == NULL) return index;

  switch (annot->type) {

    case 1 :
      for (sfp = (SeqFeatPtr) annot->data; sfp != NULL; sfp = sfp->next) {
        index += VisitSeqIdsInSeqFeat (sfp, userdata, callback);
      }
      break;

    case 2 :
      for (sap = (SeqAlignPtr) annot->data; sap != NULL; sap = sap->next) {
        index += VisitSeqIdsInSeqAlign (sap, userdata, callback);
      }
      break;

    case 3 :
      for (sgp = (SeqGraphPtr) annot->data; sgp != NULL; sgp = sgp->next) {
        index += VisitSeqIdsInSeqGraph (sgp, userdata, callback);
      }
      break;

    default :
      break;
  }

  return index;
}

NLM_EXTERN Int4 VisitUserFieldsInUfp (UserFieldPtr ufp, Pointer userdata, VisitUserFieldsFunc callback)

{
  UserFieldPtr  curr;
  Int4          index = 0;
  Boolean       nested = FALSE;

  if (ufp == NULL) return index;
  if (ufp->choice == 11) {
    for (curr = (UserFieldPtr) ufp->data.ptrvalue; curr != NULL; curr = curr->next) {
      index += VisitUserFieldsInUfp (curr, userdata,callback);
      nested = TRUE;
    }
  }
  if (! nested) {
    if (callback != NULL) {
      callback (ufp, userdata);
    }
    index++;
  }
  return index;
}

NLM_EXTERN Int4 VisitUserFieldsInUop (UserObjectPtr uop, Pointer userdata, VisitUserFieldsFunc callback)

{
  Int4          index = 0;
  UserFieldPtr  ufp;

  if (uop == NULL) return index;
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    if (callback != NULL) {
      callback (ufp, userdata);
    }
    index++;
  }
  return index;
}

/* Visits only unnested nodes */
NLM_EXTERN Int4 VisitUserObjectsInUop (UserObjectPtr uop, Pointer userdata, VisitUserObjectFunc callback)

{
  Int4           index = 0;
  Boolean        nested = FALSE;
  UserObjectPtr  obj;
  UserFieldPtr   ufp;

  if (uop == NULL) return index;
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    if (ufp->choice == 6) {
      obj = (UserObjectPtr) ufp->data.ptrvalue;
      index += VisitUserObjectsInUop (obj, userdata, callback);
      nested = TRUE;
    } else if (ufp->choice == 12) {
      for (obj = (UserObjectPtr) ufp->data.ptrvalue; obj != NULL; obj = obj->next) {
        index += VisitUserObjectsInUop (obj, userdata, callback);
      }
      nested = TRUE;
    }
  }
  if (! nested) {
    if (callback != NULL) {
      callback (uop, userdata);
    }
    index++;
  }
  return index;
}

NLM_EXTERN Int4 VisitAllUserObjectsInUop (UserObjectPtr uop, Pointer userdata, VisitUserObjectFunc callback)

{
  Int4           index = 0;
  UserObjectPtr  obj;
  UserFieldPtr   ufp;

  if (uop == NULL) return index;
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    if (ufp->choice == 6) {
      obj = (UserObjectPtr) ufp->data.ptrvalue;
      index += VisitAllUserObjectsInUop (obj, userdata, callback);
    } else if (ufp->choice == 12) {
      for (obj = (UserObjectPtr) ufp->data.ptrvalue; obj != NULL; obj = obj->next) {
        index += VisitAllUserObjectsInUop (obj, userdata, callback);
      }
    }
  }
  if (callback != NULL) {
    callback (uop, userdata);
  }
  index++;
  return index;
}

typedef struct uopdata {
  UserObjectPtr  rsult;
  CharPtr        tag;
} UopData, PNTR UopDataPtr;

static void FindUopProc (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr  oip;
  UopDataPtr   udp;

  if (uop == NULL || userdata == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  udp = (UopDataPtr) userdata;
  if (StringICmp (oip->str, udp->tag) != 0) return;
  udp->rsult = uop;
}

NLM_EXTERN UserObjectPtr FindUopByTag (UserObjectPtr top, CharPtr tag)

{
  UopData  ud;

  if (top == NULL || StringHasNoText (tag)) return NULL;
  ud.rsult = NULL;
  ud.tag = tag;
  VisitUserObjectsInUop (top, (Pointer) &ud, FindUopProc);
  return ud.rsult;
}

NLM_EXTERN UserObjectPtr CombineUserObjects (UserObjectPtr origuop, UserObjectPtr newuop)

{
  UserFieldPtr   prev = NULL;
  ObjectIdPtr    oip;
  UserFieldPtr   ufp;
  UserObjectPtr  uop;

  if (newuop == NULL) return origuop;
  if (origuop == NULL) return newuop;

  /* adding to an object that already chaperones at least two user objects */

  oip = origuop->type;
  if (oip != NULL && StringICmp (oip->str, "CombinedFeatureUserObjects") == 0) {

    for (ufp = origuop->data; ufp != NULL; ufp = ufp->next) {
      prev = ufp;
    }

    ufp = UserFieldNew ();
    oip = ObjectIdNew ();
    oip->id = 0;
    ufp->label = oip;
    ufp->choice = 6; /* user object */
    ufp->data.ptrvalue = (Pointer) newuop;

    /* link new set at end of list */

    if (prev != NULL) {
      prev->next = ufp;
    } else {
      origuop->data = ufp;
    }
    return origuop;
  }

  /* creating a new chaperone, link in first two user objects */

  uop = UserObjectNew ();
  oip = ObjectIdNew ();
  oip->str = StringSave ("CombinedFeatureUserObjects");
  uop->type = oip;

  ufp = UserFieldNew ();
  oip = ObjectIdNew ();
  oip->id = 0;
  ufp->label = oip;
  ufp->choice = 6; /* user object */
  ufp->data.ptrvalue = (Pointer) origuop;
  uop->data = ufp;
  prev = ufp;

  ufp = UserFieldNew ();
  oip = ObjectIdNew ();
  oip->id = 0;
  ufp->label = oip;
  ufp->choice = 6; /* user object */
  ufp->data.ptrvalue = (Pointer) newuop;
  prev->next = ufp;

  return uop;
}


static Int4 VisitDescriptorsProc (SeqDescrPtr descr, Pointer userdata, VisitDescriptorsFunc callback)

{
  Int4         index = 0;
  SeqDescrPtr  sdp;

  for (sdp = descr; sdp != NULL; sdp = sdp->next) {
    if (callback != NULL) {
      callback (sdp, userdata);
    }
    index++;
  }
  return index;
}

NLM_EXTERN Int4 VisitDescriptorsOnBsp (BioseqPtr bsp, Pointer userdata, VisitDescriptorsFunc callback)

{
  Int4  index = 0;

  if (bsp == NULL) return index;
  index += VisitDescriptorsProc (bsp->descr, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitDescriptorsOnSet (BioseqSetPtr bssp, Pointer userdata, VisitDescriptorsFunc callback)

{
  Int4  index = 0;

  if (bssp == NULL) return index;
  index += VisitDescriptorsProc (bssp->descr, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitDescriptorsInSet (BioseqSetPtr bssp, Pointer userdata, VisitDescriptorsFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  index += VisitDescriptorsProc (bssp->descr, userdata, callback);
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitDescriptorsInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitDescriptorsOnSep (SeqEntryPtr sep, Pointer userdata, VisitDescriptorsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitDescriptorsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitDescriptorsOnSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitDescriptorsInSep (SeqEntryPtr sep, Pointer userdata, VisitDescriptorsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitDescriptorsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitDescriptorsInSet (bssp, userdata, callback);
  }
  return index;
}


static Int4 VisitFeaturesProc (SeqAnnotPtr annot, Pointer userdata, VisitFeaturesFunc callback)

{
  Int4         index = 0;
  SeqAnnotPtr  sap;
  SeqFeatPtr   sfp;

  for (sap = annot; sap != NULL; sap = sap->next) {
    if (sap->type != 1) continue;
    for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
      if (callback != NULL) {
        callback (sfp, userdata);
      }
      index++;
    }
  }
  return index;
}

NLM_EXTERN Int4 VisitFeaturesOnSap (SeqAnnotPtr sap, Pointer userdata, VisitFeaturesFunc callback)

{
  Int4        index = 0;
  SeqFeatPtr  sfp;

  if (sap == NULL) return index;
  if (sap->type != 1) return index;
  for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
    if (callback != NULL) {
      callback (sfp, userdata);
    }
    index++;
  }
  return index;
}

NLM_EXTERN Int4 VisitFeaturesOnBsp (BioseqPtr bsp, Pointer userdata, VisitFeaturesFunc callback)

{
  Int4  index = 0;

  if (bsp == NULL) return index;
  index += VisitFeaturesProc (bsp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitFeaturesOnSet (BioseqSetPtr bssp, Pointer userdata, VisitFeaturesFunc callback)

{
  Int4  index = 0;

  if (bssp == NULL) return index;
  index += VisitFeaturesProc (bssp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitFeaturesInSet (BioseqSetPtr bssp, Pointer userdata, VisitFeaturesFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  index += VisitFeaturesProc (bssp->annot, userdata, callback);
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitFeaturesInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitFeaturesOnSep (SeqEntryPtr sep, Pointer userdata, VisitFeaturesFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitFeaturesOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitFeaturesOnSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitFeaturesInSep (SeqEntryPtr sep, Pointer userdata, VisitFeaturesFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitFeaturesOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitFeaturesInSet (bssp, userdata, callback);
  }
  return index;
}


static Int4 VisitAlignmentsOnDisc (Pointer segs, Pointer userdata, VisitAlignmentsFunc callback)

{
  Int4         index = 0;
  SeqAlignPtr  salp;

  for (salp = (SeqAlignPtr) segs; salp != NULL; salp = salp->next) {
    if (callback != NULL) {
      callback (salp, userdata);
    }
    index++;
    if (salp->segtype == SAS_DISC) {
      index += VisitAlignmentsOnDisc (salp->segs, userdata, callback);
    }
  }
  return index;
}

static Int4 VisitAlignmentsProc (SeqAnnotPtr annot, Pointer userdata, VisitAlignmentsFunc callback)

{
  Int4         index = 0;
  SeqAlignPtr  salp;
  SeqAnnotPtr  sap;

  for (sap = annot; sap != NULL; sap = sap->next) {
    if (sap->type != 2) continue;
    for (salp = (SeqAlignPtr) sap->data; salp != NULL; salp = salp->next) {
      if (callback != NULL) {
        callback (salp, userdata);
      }
      index++;
      if (salp->segtype == SAS_DISC) {
        index += VisitAlignmentsOnDisc (salp->segs, userdata, callback);
      }
    }
  }
  return index;
}

NLM_EXTERN Int4 VisitAlignmentsOnSap (SeqAnnotPtr sap, Pointer userdata, VisitAlignmentsFunc callback)

{
  Int4         index = 0;
  SeqAlignPtr  salp;

  if (sap == NULL) return index;
  if (sap->type != 2) return index;
  for (salp = (SeqAlignPtr) sap->data; salp != NULL; salp = salp->next) {
    if (callback != NULL) {
      callback (salp, userdata);
    }
    index++;
    if (salp->segtype == SAS_DISC) {
      index += VisitAlignmentsOnDisc (salp->segs, userdata, callback);
    }
  }
  return index;
}

NLM_EXTERN Int4 VisitAlignmentsOnBsp (BioseqPtr bsp, Pointer userdata, VisitAlignmentsFunc callback)

{
  Int4  index = 0;

  if (bsp == NULL) return index;
  index += VisitAlignmentsProc (bsp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitAlignmentsOnSet (BioseqSetPtr bssp, Pointer userdata, VisitAlignmentsFunc callback)

{
  Int4  index = 0;

  if (bssp == NULL) return index;
  index += VisitAlignmentsProc (bssp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitAlignmentsInSet (BioseqSetPtr bssp, Pointer userdata, VisitAlignmentsFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  index += VisitAlignmentsProc (bssp->annot, userdata, callback);
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitAlignmentsInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitAlignmentsOnSep (SeqEntryPtr sep, Pointer userdata, VisitAlignmentsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitAlignmentsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitAlignmentsOnSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitAlignmentsInSep (SeqEntryPtr sep, Pointer userdata, VisitAlignmentsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitAlignmentsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitAlignmentsInSet (bssp, userdata, callback);
  }
  return index;
}


static Int4 VisitGraphsProc (SeqAnnotPtr annot, Pointer userdata, VisitGraphsFunc callback)

{
  Int4         index = 0;
  SeqAnnotPtr  sap;
  SeqGraphPtr  sgp;

  for (sap = annot; sap != NULL; sap = sap->next) {
    if (sap->type != 3) continue;
    for (sgp = (SeqGraphPtr) sap->data; sgp != NULL; sgp = sgp->next) {
      if (callback != NULL) {
        callback (sgp, userdata);
      }
      index++;
    }
  }
  return index;
}

NLM_EXTERN Int4 VisitGraphsOnSap (SeqAnnotPtr sap, Pointer userdata, VisitGraphsFunc callback)

{
  Int4         index = 0;
  SeqGraphPtr  sgp;

  if (sap == NULL) return index;
  if (sap->type != 3) return index;
  for (sgp = (SeqGraphPtr) sap->data; sgp != NULL; sgp = sgp->next) {
    if (callback != NULL) {
      callback (sgp, userdata);
    }
    index++;
  }
  return index;
}

NLM_EXTERN Int4 VisitGraphsOnBsp (BioseqPtr bsp, Pointer userdata, VisitGraphsFunc callback)

{
  Int4  index = 0;

  if (bsp == NULL) return index;
  index += VisitGraphsProc (bsp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitGraphsOnSet (BioseqSetPtr bssp, Pointer userdata, VisitGraphsFunc callback)

{
  Int4  index = 0;

  if (bssp == NULL) return index;
  index += VisitGraphsProc (bssp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitGraphsInSet (BioseqSetPtr bssp, Pointer userdata, VisitGraphsFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  index += VisitGraphsProc (bssp->annot, userdata, callback);
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitGraphsInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitGraphsOnSep (SeqEntryPtr sep, Pointer userdata, VisitGraphsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitGraphsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitGraphsOnSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitGraphsInSep (SeqEntryPtr sep, Pointer userdata, VisitGraphsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitGraphsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitGraphsInSet (bssp, userdata, callback);
  }
  return index;
}


static Int4 VisitAnnotsProc (SeqAnnotPtr annot, Pointer userdata, VisitAnnotsFunc callback)

{
  Int4         index = 0;
  SeqAnnotPtr  sap;

  for (sap = annot; sap != NULL; sap = sap->next) {
    if (callback != NULL) {
      callback (sap, userdata);
    }
    index++;
  }
  return index;
}

NLM_EXTERN Int4 VisitAnnotsOnBsp (BioseqPtr bsp, Pointer userdata, VisitAnnotsFunc callback)

{
  Int4  index = 0;

  if (bsp == NULL) return index;
  index += VisitAnnotsProc (bsp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitAnnotsOnSet (BioseqSetPtr bssp, Pointer userdata, VisitAnnotsFunc callback)

{
  Int4  index = 0;

  if (bssp == NULL) return index;
  index += VisitAnnotsProc (bssp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitAnnotsInSet (BioseqSetPtr bssp, Pointer userdata, VisitAnnotsFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  index += VisitAnnotsProc (bssp->annot, userdata, callback);
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitAnnotsInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitAnnotsOnSep (SeqEntryPtr sep, Pointer userdata, VisitAnnotsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitAnnotsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitAnnotsOnSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitAnnotsInSep (SeqEntryPtr sep, Pointer userdata, VisitAnnotsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitAnnotsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitAnnotsInSet (bssp, userdata, callback);
  }
  return index;
}


static Int4 VisitAuthorsProc (AuthListPtr alp, Pointer userdata, VisitAuthorFunc callback)

{
  AuthorPtr    ap;
  Int4         index = 0;
  ValNodePtr   names;
  NameStdPtr   nsp;
  PersonIdPtr  pid;

  if (alp == NULL || alp->choice != 1) return index;

  for (names = alp->names; names != NULL; names = names->next) {
    ap = names->data.ptrvalue;
    if (ap == NULL) continue;
    pid = ap->name;
    if (pid == NULL || pid->choice != 2) continue;
    nsp = pid->data;
    if (nsp == NULL) continue;
    if (callback != NULL) {
      callback (nsp, userdata);
    }
    index++;
  }

  return index;
}

NLM_EXTERN Int4 VisitAuthorsInPub (PubdescPtr pdp, Pointer userdata, VisitAuthorFunc callback)

{
  CitArtPtr   cap;
  CitBookPtr  cbp;
  CitGenPtr   cgp;
  CitPatPtr   cpp;
  CitSubPtr   csp;
  Int4        index = 0;
  ValNodePtr  vnp;

  if (pdp == NULL) return index;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_PMid || vnp->choice == PUB_Muid) continue;
    if (vnp->data.ptrvalue == NULL) continue;
    switch (vnp->choice) {
      case PUB_Gen :
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        index += VisitAuthorsProc (cgp->authors, userdata, callback);
        break;
      case PUB_Sub :
        csp = (CitSubPtr) vnp->data.ptrvalue;
        index += VisitAuthorsProc (csp->authors, userdata, callback);
        break;
      case PUB_Article :
        cap = (CitArtPtr) vnp->data.ptrvalue;
        index += VisitAuthorsProc (cap->authors, userdata, callback);
        if (cap->from == 2 || cap->from == 3) {
          cbp = (CitBookPtr) cap->fromptr;
          if (cbp != NULL) {
            index += VisitAuthorsProc (cbp->authors, userdata, callback);
          }
        }
        break;
      case PUB_Book :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        index += VisitAuthorsProc (cbp->authors, userdata, callback);
        break;
      case PUB_Man :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp->othertype == 2 && cbp->let_type == 3) {
          index += VisitAuthorsProc (cbp->authors, userdata, callback);
        }
        break;
      case PUB_Patent :
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        index += VisitAuthorsProc (cpp->authors, userdata, callback);
        index += VisitAuthorsProc (cpp->applicants, userdata, callback);
        index += VisitAuthorsProc (cpp->assignees, userdata, callback);
        break;
      default :
        break;
    }
  }

  return index;
}


static Int4 VisitPubdescsProc (SeqDescrPtr descr, SeqAnnotPtr annot, Pointer userdata, VisitPubdescsFunc callback)

{
  Int4         index = 0;
  PubdescPtr   pdp;
  SeqAnnotPtr  sap;
  SeqDescrPtr  sdp;
  SeqFeatPtr   sfp;

  for (sdp = descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_pub) {
      pdp = (PubdescPtr) sdp->data.ptrvalue;
      if (pdp != NULL) {
        if (callback != NULL) {
          callback (pdp, userdata);
        }
        index++;
      }
    }
  }
  for (sap = annot; sap != NULL; sap = sap->next) {
    if (sap->type != 1) continue;
    for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
      if (sfp->data.choice == SEQFEAT_PUB) {
        pdp = (PubdescPtr) sfp->data.value.ptrvalue;
        if (pdp != NULL) {
          if (callback != NULL) {
            callback (pdp, userdata);
          }
          index++;
        }
      }
    }
  }
  return index;
}

NLM_EXTERN Int4 VisitPubdescsOnBsp (BioseqPtr bsp, Pointer userdata, VisitPubdescsFunc callback)

{
  Int4  index = 0;

  if (bsp == NULL) return index;
  index += VisitPubdescsProc (bsp->descr, bsp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitPubdescsOnSet (BioseqSetPtr bssp, Pointer userdata, VisitPubdescsFunc callback)

{
  Int4  index = 0;

  if (bssp == NULL) return index;
  index += VisitPubdescsProc (bssp->descr, bssp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitPubdescsInSet (BioseqSetPtr bssp, Pointer userdata, VisitPubdescsFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  index += VisitPubdescsProc (bssp->descr, bssp->annot, userdata, callback);
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitPubdescsInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitPubdescsOnSep (SeqEntryPtr sep, Pointer userdata, VisitPubdescsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitPubdescsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitPubdescsOnSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitPubdescsInSep (SeqEntryPtr sep, Pointer userdata, VisitPubdescsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitPubdescsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitPubdescsInSet (bssp, userdata, callback);
  }
  return index;
}


static Int4 VisitBioSourcesProc (SeqDescrPtr descr, SeqAnnotPtr annot, Pointer userdata, VisitBioSourcesFunc callback)

{
  BioSourcePtr  biop;
  Int4          index = 0;
  SeqAnnotPtr   sap;
  SeqDescrPtr   sdp;
  SeqFeatPtr    sfp;

  for (sdp = descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_source) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop != NULL) {
        if (callback != NULL) {
          callback (biop, userdata);
        }
        index++;
      }
    }
  }
  for (sap = annot; sap != NULL; sap = sap->next) {
    if (sap->type != 1) continue;
    for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
      if (sfp->data.choice == SEQFEAT_BIOSRC) {
        biop = (BioSourcePtr) sfp->data.value.ptrvalue;
        if (biop != NULL) {
          if (callback != NULL) {
            callback (biop, userdata);
          }
          index++;
        }
      }
    }
  }
  return index;
}

NLM_EXTERN Int4 VisitBioSourcesOnBsp (BioseqPtr bsp, Pointer userdata, VisitBioSourcesFunc callback)

{
  Int4  index = 0;

  if (bsp == NULL) return index;
  index += VisitBioSourcesProc (bsp->descr, bsp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitBioSourcesOnSet (BioseqSetPtr bssp, Pointer userdata, VisitBioSourcesFunc callback)

{
  Int4  index = 0;

  if (bssp == NULL) return index;
  index += VisitBioSourcesProc (bssp->descr, bssp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitBioSourcesInSet (BioseqSetPtr bssp, Pointer userdata, VisitBioSourcesFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  index += VisitBioSourcesProc (bssp->descr, bssp->annot, userdata, callback);
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitBioSourcesInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitBioSourcesOnSep (SeqEntryPtr sep, Pointer userdata, VisitBioSourcesFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitBioSourcesOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitBioSourcesOnSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitBioSourcesInSep (SeqEntryPtr sep, Pointer userdata, VisitBioSourcesFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitBioSourcesOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitBioSourcesInSet (bssp, userdata, callback);
  }
  return index;
}


NLM_EXTERN Int4 VisitBioseqsInSet (BioseqSetPtr bssp, Pointer userdata, VisitBioseqsFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitBioseqsInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitBioseqsInSep (SeqEntryPtr sep, Pointer userdata, VisitBioseqsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (callback != NULL) {
      callback (bsp, userdata);
    }
    index++;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitBioseqsInSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitSequencesInSet (BioseqSetPtr bssp, Pointer userdata, Int2 filter, VisitSequencesFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  if (bssp->_class == BioseqseqSet_class_parts) {
    if (filter != VISIT_PARTS) return index;
    filter = VISIT_MAINS;
  }
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitSequencesInSep (tmp, userdata, filter, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitSequencesInSep (SeqEntryPtr sep, Pointer userdata, Int2 filter, VisitSequencesFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (filter == VISIT_MAINS ||
        (filter == VISIT_NUCS && ISA_na (bsp->mol)) ||
        (filter == VISIT_PROTS && ISA_aa (bsp->mol))) {
      if (callback != NULL) {
        callback (bsp, userdata);
      }
      index++;
    }
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitSequencesInSet (bssp, userdata, filter, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitSetsInSet (BioseqSetPtr bssp, Pointer userdata, VisitSetsFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  if (callback != NULL) {
    callback (bssp, userdata);
  }
  index++;
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitSetsInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitSetsInSep (SeqEntryPtr sep, Pointer userdata, VisitSetsFunc callback)

{
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitSetsInSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitElementsInSep (SeqEntryPtr sep, Pointer userdata, VisitElementsFunc callback)

{
  BioseqSetPtr  bssp;
  Int4          index = 0;
  SeqEntryPtr   tmp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return index;
    if (bssp->_class == 7 ||
        (bssp->_class >= 13 && bssp->_class <= 16) ||
        bssp->_class == BioseqseqSet_class_wgs_set ||
        bssp->_class == BioseqseqSet_class_gen_prod_set ||
        bssp->_class == BioseqseqSet_class_small_genome_set) {
      for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
        index += VisitElementsInSep (tmp, userdata, callback);
      }
      return index;
    }
  }
  if (callback != NULL) {
    callback (sep, userdata);
  }
  index++;
  return index;
}

NLM_EXTERN Boolean IsPopPhyEtcSet (Uint1 _class)

{
  if (_class == BioseqseqSet_class_mut_set ||
      _class == BioseqseqSet_class_pop_set ||
      _class == BioseqseqSet_class_phy_set ||
      _class == BioseqseqSet_class_eco_set ||
      _class == BioseqseqSet_class_wgs_set ||
      _class == BioseqseqSet_class_small_genome_set) return TRUE;
  return FALSE;
}


static Int4 ScanBioseqSetReleaseInt (
  CharPtr inputFile,
  Boolean binary,
  Boolean compressed,
  Pointer userdata,
  ScanBioseqSetFunc callback,
  Boolean freesep,
  TNlmMutexPtr mutex
)

{
  AsnIoPtr      aip;
  AsnModulePtr  amp;
  AsnTypePtr    atp, atp_bss, atp_se;
  FILE          *fp;
  Int4          index = 0;
  SeqEntryPtr   sep;
#ifdef OS_UNIX
  Char          cmmd [256];
  CharPtr       gzcatprog;
  int           ret;
  Boolean       usedPopen = FALSE;
#endif
  if (StringHasNoText (inputFile) || callback == NULL) return index; 

#ifndef OS_UNIX
  if (compressed) {
    Message (MSG_ERROR, "Can only decompress on-the-fly on UNIX machines");
    return index;
  }
#endif

  amp = AsnAllModPtr ();
  if (amp == NULL) {
    Message (MSG_ERROR, "Unable to load AsnAllModPtr");
    return index;
  }

  atp_bss = AsnFind ("Bioseq-set");
  if (atp_bss == NULL) {
    Message (MSG_ERROR, "Unable to find ASN.1 type Bioseq-set");
    return index;
  }

  atp_se = AsnFind ("Bioseq-set.seq-set.E");
  if (atp_se == NULL) {
    Message (MSG_ERROR, "Unable to find ASN.1 type Bioseq-set.seq-set.E");
    return index;
  }

#ifdef OS_UNIX
  if (compressed) {
    gzcatprog = getenv ("NCBI_UNCOMPRESS_BINARY");
    if (gzcatprog != NULL) {
      sprintf (cmmd, "%s %s", gzcatprog, inputFile);
    } else {
      ret = system ("gzcat -h >/dev/null 2>&1");
      if (ret == 0) {
        sprintf (cmmd, "gzcat %s", inputFile);
      } else if (ret == -1) {
        Message (MSG_FATAL, "Unable to fork or exec gzcat in ScanBioseqSetRelease");
        return index;
      } else {
        ret = system ("zcat -h >/dev/null 2>&1");
        if (ret == 0) {
          sprintf (cmmd, "zcat %s", inputFile);
        } else if (ret == -1) {
          Message (MSG_FATAL, "Unable to fork or exec zcat in ScanBioseqSetRelease");
          return index;
        } else {
          Message (MSG_FATAL, "Unable to find zcat or gzcat in ScanBioseqSetRelease - please edit your PATH environment variable");
          return index;
        }
      }
    }
    fp = popen (cmmd, /* binary? "rb" : */ "r");
    usedPopen = TRUE;
  } else {
    fp = FileOpen (inputFile, binary? "rb" : "r");
  }
#else
  fp = FileOpen (inputFile, binary? "rb" : "r");
#endif
  if (fp == NULL) {
    Message (MSG_POSTERR, "FileOpen failed for input file '%s'", inputFile);
    return index;
  }

  aip = AsnIoNew (binary? ASNIO_BIN_IN : ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
  if (aip == NULL) {
    Message (MSG_ERROR, "AsnIoNew failed for input file '%s'", inputFile);
    return index;
  }

  atp = atp_bss;

  while ((atp = AsnReadId (aip, amp, atp)) != NULL) {
    if (atp == atp_se) {
      if (mutex != NULL) {
        NlmMutexLockEx (mutex);
      }
      SeqMgrHoldIndexing (TRUE);
      sep = SeqEntryAsnRead (aip, atp);
      SeqMgrHoldIndexing (FALSE);
      if (mutex != NULL) {
        NlmMutexUnlock (*mutex);
      }
      callback (sep, userdata);
      if (freesep) {
        SeqEntryFree (sep);
      }
      index++;
    } else {
      AsnReadVal (aip, atp, NULL);
    }
  }

  AsnIoFree (aip, FALSE);

#ifdef OS_UNIX
  if (usedPopen) {
    pclose (fp);
  } else {
    FileClose (fp);
  }
#else
  FileClose (fp);
#endif
  return index;
}

NLM_EXTERN Int4 ScanBioseqSetRelease (
  CharPtr inputFile,
  Boolean binary,
  Boolean compressed,
  Pointer userdata,
  ScanBioseqSetFunc callback
)

{
  return ScanBioseqSetReleaseInt (inputFile, binary, compressed, userdata, callback, TRUE, NULL);
}

static TNlmMutex  scan_bioseq_set_release_mutex = NULL;

NLM_EXTERN Int4 ScanBioseqSetReleaseMT (
  CharPtr inputFile,
  Boolean binary,
  Boolean compressed,
  Pointer userdata,
  ScanBioseqSetFunc callback
)

{
  return ScanBioseqSetReleaseInt (inputFile, binary, compressed, userdata, callback, FALSE, &scan_bioseq_set_release_mutex);
}

NLM_EXTERN SeqEntryPtr LIBCALL FreeScanSeqEntryMT (
  SeqEntryPtr sep
)

{
  if (sep == NULL) return NULL;

  NlmMutexLockEx (&scan_bioseq_set_release_mutex);

  SeqMgrHoldIndexing (TRUE);
  SeqEntryFree (sep);
  SeqMgrHoldIndexing (FALSE);

  NlmMutexUnlock (scan_bioseq_set_release_mutex);

  return NULL;
}

NLM_EXTERN Int4 ScanEntrezgeneSetRelease (
  CharPtr inputFile,
  Boolean binary,
  Boolean compressed,
  Pointer userdata,
  ScanEntrezgeneSetFunc callback
)

{
  AsnIoPtr       aip;
  AsnModulePtr   amp;
  AsnTypePtr     atp, atp_egs, atp_egse;
  EntrezgenePtr  egp;
  FILE           *fp;
  Int4           index = 0;
#ifdef OS_UNIX
  Char           cmmd [256];
  CharPtr        gzcatprog;
  int            ret;
  Boolean        usedPopen = FALSE;
#endif
  if (StringHasNoText (inputFile) || callback == NULL) return index; 

#ifndef OS_UNIX
  if (compressed) {
    Message (MSG_ERROR, "Can only decompress on-the-fly on UNIX machines");
    return index;
  }
#endif

  amp = AsnAllModPtr ();
  if (amp == NULL) {
    Message (MSG_ERROR, "Unable to load AsnAllModPtr");
    return index;
  }

  atp_egs = AsnFind ("Entrezgene-Set");
  if (atp_egs == NULL) {
    Message (MSG_ERROR, "Unable to find ASN.1 type Entrezgene-Set");
    return index;
  }

  atp_egse = AsnFind ("Entrezgene-Set.E");
  if (atp_egse == NULL) {
    Message (MSG_ERROR, "Unable to find ASN.1 type Entrezgene-Set.E");
    return index;
  }

#ifdef OS_UNIX
  if (compressed) {
    gzcatprog = getenv ("NCBI_UNCOMPRESS_BINARY");
    if (gzcatprog != NULL) {
      sprintf (cmmd, "%s %s", gzcatprog, inputFile);
    } else {
      ret = system ("gzcat -h >/dev/null 2>&1");
      if (ret == 0) {
        sprintf (cmmd, "gzcat %s", inputFile);
      } else if (ret == -1) {
        Message (MSG_FATAL, "Unable to fork or exec gzcat in ScanEntrezgeneSetRelease");
        return index;
      } else {
        ret = system ("zcat -h >/dev/null 2>&1");
        if (ret == 0) {
          sprintf (cmmd, "zcat %s", inputFile);
        } else if (ret == -1) {
          Message (MSG_FATAL, "Unable to fork or exec zcat in ScanEntrezgeneSetRelease");
          return index;
        } else {
          Message (MSG_FATAL, "Unable to find zcat or gzcat in ScanEntrezgeneSetRelease - please edit your PATH environment variable");
          return index;
        }
      }
    }
    fp = popen (cmmd, /* binary? "rb" : */ "r");
    usedPopen = TRUE;
  } else {
    fp = FileOpen (inputFile, binary? "rb" : "r");
  }
#else
  fp = FileOpen (inputFile, binary? "rb" : "r");
#endif
  if (fp == NULL) {
    Message (MSG_POSTERR, "FileOpen failed for input file '%s'", inputFile);
    return index;
  }

  aip = AsnIoNew (binary? ASNIO_BIN_IN : ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
  if (aip == NULL) {
    Message (MSG_ERROR, "AsnIoNew failed for input file '%s'", inputFile);
    return index;
  }

  atp = atp_egs;

  while ((atp = AsnReadId (aip, amp, atp)) != NULL) {
    if (atp == atp_egse) {
      egp = EntrezgeneAsnRead (aip, atp);
      callback (egp, userdata);
      EntrezgeneFree (egp);
      index++;
    } else {
      AsnReadVal (aip, atp, NULL);
    }
  }

  AsnIoFree (aip, FALSE);

#ifdef OS_UNIX
  if (usedPopen) {
    pclose (fp);
  } else {
    FileClose (fp);
  }
#else
  FileClose (fp);
#endif
  return index;
}



NLM_EXTERN int LIBCALL ObjectIdCompare (ObjectIdPtr a, ObjectIdPtr b)
{
  int rval = 0;
  Char buf[30];

    if (a == b) {
        rval = 0;
  } else if (a == NULL) {
    rval = -1;
  } else if (b == NULL) {
    rval = 1;
  } else if (a->str == NULL && b->str == NULL) {
    if (a->id < b->id) {
      rval = -1;
    } else if (a->id > b->id) {
      rval = 1;
    }
  } else if (a->str == NULL) {
    sprintf (buf, "%d", a->id);
    rval = StringCmp (buf, b->str);
  } else if (b->str == NULL) {
    sprintf (buf, "%d", b->id);
    rval = StringCmp (a->str, buf);
  } else {
    rval = StringCmp (a->str, b->str);
  }
  return rval; 
}


/*****************************************************************************
*
*   DbtagMatch(a, b)
*
*****************************************************************************/
NLM_EXTERN int LIBCALL DbtagCompare (DbtagPtr a, DbtagPtr b)
{
  int rval = 0;

    if (a == b) {
        rval = 0;
  } else if (a == NULL) {
    rval = -1;
  } else if (b == NULL) {
    rval = 1;
  } else if ((rval = StringICmp (a->db, b->db)) == 0) {
    rval = ObjectIdCompare (a->tag, b->tag);
  }
  return rval;
}


static int LIBCALLBACK SortVnpByDbtag (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      return DbtagCompare (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    }
  }
  return 0;
}

NLM_EXTERN int LIBCALL OrgModSetCompare (OrgModPtr mod1, OrgModPtr mod2)
{
  int rval = 0;

  while (mod1 != NULL && mod2 != NULL && rval == 0)
  {
    if (mod1->subtype < mod2->subtype)
    {
      rval = -1;
    }
    else if (mod1->subtype > mod2->subtype)
    {
      rval = 1;
    }
    else if ((rval = StringCmp (mod1->subname, mod2->subname)) == 0
           && (rval = StringCmp (mod1->attrib, mod2->attrib)) == 0)
    {
      mod1 = mod1->next;
      mod2 = mod2->next;
    }
  }

  if (rval == 0)
  {
    if (mod1 == NULL && mod2 == NULL)
    {
      rval = 0;
    }
    else if (mod1 == NULL)
    {
      rval = -1;
    }
    else if (mod2 == NULL)
    {
      rval = 1;
    }
  }
  return rval;
}


NLM_EXTERN int LIBCALL OrgNameCompare (OrgNamePtr onp1, OrgNamePtr onp2)
{
  int rval = 0;

  while (onp1 != NULL && onp2 != NULL && rval == 0)
  {
    if ((rval = OrgModSetCompare(onp1->mod, onp2->mod)) != 0
        || (rval = StringCmp (onp1->lineage, onp2->lineage)) != 0
        || (rval = StringCmp (onp1->div, onp2->div)) != 0
        || (rval = StringCmp (onp1->attrib, onp2->attrib)) != 0)
    {
      /* no further processing */
    }
    else if (onp1->choice < onp2->choice)
    {
      rval = -1;
    }
    else if (onp1->choice > onp2->choice)
    {
      rval = 1;
    }
    else if (onp1->gcode < onp2->gcode)
    {
      rval = -1;
    } 
    else if (onp1->gcode > onp2->gcode)
    {
      rval = 1;
    }
    else if (onp1->mgcode < onp2->mgcode)
    {
      rval = -1;
    } 
    else if (onp1->mgcode > onp2->mgcode)
    {
      rval = 1;
    }
    else if (onp1->pgcode < onp2->pgcode)
    {
      rval = -1;
    } 
    else if (onp1->pgcode > onp2->pgcode)
    {
      rval = 1;
    }
    onp1 = onp1->next;
    onp2 = onp2->next;
  }
  if (rval == 0) 
  {
    if (onp1 == NULL && onp2 == NULL) 
    {
      rval = 0;
    }
    else if (onp1 == NULL)
    {
      rval = -1;
    } 
    else if (onp2 == NULL) 
    {
      rval = 1;
    }
  }
  return rval;
}


/*****************************************************************************
*
*   OrgRefCompare (orp1, orp2)
*
*****************************************************************************/
NLM_EXTERN int LIBCALL OrgRefCompare (OrgRefPtr orp1, OrgRefPtr orp2)
{
  int rval = 0;
  if (orp1 == NULL && orp2 == NULL)
  {
    return 0;
  }
  else if (orp1 == NULL) 
  {
    return -1;
  }
  else if (orp2 == NULL)
  {
    return 1;
  }
  else if ((rval = StringCmp (orp1->taxname, orp2->taxname)) != 0) 
  {
    return rval;
  }
  else if ((rval = StringCmp (orp1->common, orp2->common)) != 0)
  {
    return rval;
  }
  else if ((rval = ValNodeCompare (orp1->syn, orp2->syn, SortVnpByString)) != 0) 
  {
    return rval;
  }
  else if ((rval = ValNodeCompare (orp1->db, orp2->db, SortVnpByDbtag)) != 0)
  {
    return rval;
  }
  else
  {
    rval = OrgNameCompare (orp1->orgname, orp2->orgname);
  }
  return rval;
}


static Boolean DoStringsMatch (CharPtr str1, CharPtr str2, Boolean case_sensitive)
{
  Boolean rval = FALSE;

  if (case_sensitive) {
    if (StringCmp (str1, str2) == 0) {
      rval = TRUE;
    }
  } else if (StringICmp (str1, str2) == 0) {
    rval = TRUE;
  }
  return rval;
}


static Boolean DoGBQualListsMatch (GBQualPtr gbq1, GBQualPtr gbq2, Boolean case_sensitive)
{
  Boolean rval = TRUE;

  while (rval && gbq1 != NULL && gbq2 != NULL) {
    if (!DoStringsMatch (gbq1->qual, gbq2->qual, case_sensitive)) {
      rval = FALSE;
    } else if (!DoStringsMatch (gbq1->val, gbq2->val, case_sensitive)) {
      rval = FALSE;
    } else {
      gbq1 = gbq1->next;
      gbq2 = gbq2->next;
    }
  }
  if (gbq1 != NULL || gbq2 != NULL) {
    rval = FALSE;
  }
  return rval;
}


static Boolean CheckBioseqForPartial (BioseqPtr bsp, BoolPtr partial5, BoolPtr partial3)
{
  SeqMgrDescContext context;
  SeqDescrPtr       sdp;
  MolInfoPtr        mip;
  Boolean           rval = FALSE;

  if (bsp == NULL) {
    return FALSE;
  }
  if (partial5 != NULL) {
    *partial5 = FALSE;
  }
  if (partial3 != NULL) {
    *partial3 = FALSE;
  }
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
  if (sdp != NULL && (mip = (MolInfoPtr) sdp->data.ptrvalue) != NULL) {
    /* partial 5 */
    if (mip->completeness == 3 || mip->completeness == 5) {
      if (partial5 != NULL) {
        *partial5 = TRUE;
      }
      rval = TRUE;
    }
    /* partial 3 */
    if (mip->completeness == 4 || mip->completeness == 5) {
      if (partial3 != NULL) {
        *partial3 = TRUE;
      }
      rval = TRUE;
    }
    if (mip->completeness == 2) {
      rval = TRUE;
    }
  }
  return rval;
}


static Boolean ProductsMatch (SeqLocPtr slp1, SeqLocPtr slp2, Boolean case_sensitive, Boolean ignore_partial)
{
  BioseqPtr bsp1, bsp2;
  Int2      ctr, pos1, pos2;
  Char      buf1[51];
  Char      buf2[51];
  Int4      len = 50;
  SeqFeatPtr sfp1, sfp2;
  SeqMgrFeatContext fcontext1, fcontext2;
  Boolean           partial5_1, partial5_2, partial3_1, partial3_2;

  if (slp1 == NULL && slp2 == NULL) {
    return TRUE;
  } else if (slp1 == NULL || slp2 == NULL) {
    return FALSE;
  } else if (SeqLocCompare (slp1, slp2) == SLC_A_EQ_B) {
    return TRUE;
  } else {
    bsp1 = BioseqFindFromSeqLoc (slp1);
    bsp2 = BioseqFindFromSeqLoc (slp2);
    if (bsp1 == NULL || bsp2 == NULL) {
      /* can't compare, assume they don't match */
      return FALSE;
    } else if (bsp1->length != bsp2->length) {
      return FALSE;
    } else {
      CheckBioseqForPartial (bsp1, &partial5_1, &partial3_1);
      CheckBioseqForPartial (bsp2, &partial5_2, &partial3_2);
      if (!ignore_partial 
          && ((partial5_1 && !partial5_2) || (!partial5_1 && partial5_2) 
          || (partial3_1 && !partial3_2) || (!partial3_1 && partial3_2))) {
        return FALSE;
      }
      /* check that translation sequences match */
      pos1 = 0;
      pos2 = 0;
      if (ignore_partial) {
        if (partial5_1 || partial5_2) {
          pos1++;
          pos2++;
        }
      }
      while (pos1 < bsp1->length && pos2 < bsp2->length) {
        ctr = SeqPortStreamInt (bsp1, pos1, MIN(pos1 + len - 1, bsp1->length - 1), Seq_strand_plus,
                            STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL,
                            (Pointer) buf1, NULL);
        ctr = SeqPortStreamInt (bsp2, pos2, MIN(pos2 + len - 1, bsp2->length - 1), Seq_strand_plus,
                            STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL,
                            (Pointer) buf2, NULL);
        if (StringNCmp (buf1, buf2, ctr) != 0) {
          return FALSE;
        }
        pos1 += len;
        pos2 += len;
      }

      /* now check that protein features match */
      sfp1 = SeqMgrGetNextFeature (bsp1, NULL, 0, 0, &fcontext1);
      sfp2 = SeqMgrGetNextFeature (bsp2, NULL, 0, 0, &fcontext2);
      while (sfp1 != NULL && sfp2 != NULL) {
        if (!DoFeaturesMatch (sfp1, sfp2, TRUE, case_sensitive, ignore_partial)) {
          return FALSE;
        }        
        sfp1 = SeqMgrGetNextFeature (bsp1, sfp1, SEQFEAT_PROT, 0, &fcontext1);
        sfp2 = SeqMgrGetNextFeature (bsp2, sfp2, SEQFEAT_PROT, 0, &fcontext2);
      }
      if (sfp1 != NULL || sfp2 != NULL) {
        return FALSE;
      } else {
        return TRUE;
      }
    }
  }
}


static Boolean DoLocationPartialsMatch (SeqLocPtr slp1, SeqLocPtr slp2)
{
  Boolean partial5_1, partial3_1, partial1;
  Boolean partial5_2, partial3_2, partial2;

  partial1 = CheckSeqLocForPartial (slp1, &partial5_1, &partial3_1);
  partial2 = CheckSeqLocForPartial (slp2, &partial5_2, &partial3_2);
  if ((partial1 && !partial2) || (!partial1 && partial2)) {
    return FALSE;
  }
  if ((partial5_1 && !partial5_2) || (!partial5_1 && partial5_2)) {
    return FALSE;
  }
  if ((partial3_1 && !partial3_2) || (!partial3_1 && partial3_2)) {
    return FALSE;
  }
  return TRUE;
}


static Boolean DoLocationsMatch (SeqLocPtr slp1, SeqLocPtr slp2, Boolean allow_different_sequences, Boolean ignore_partial)
{
  SeqLocPtr slp_tmp1, slp_tmp2;

  if (slp1 == NULL && slp2 == NULL) {
    return TRUE;
  } else if (slp1 == NULL || slp2 == NULL) {
    return FALSE;
  }

  if (!ignore_partial && !DoLocationPartialsMatch (slp1, slp2)) {
    return FALSE;
  }
  if (allow_different_sequences) {
    for (slp_tmp1 = SeqLocFindNext (slp1, NULL), slp_tmp2 = SeqLocFindNext (slp2, NULL);
         slp_tmp1 != NULL && slp_tmp2 != NULL;
         slp_tmp1 = SeqLocFindNext (slp1, slp_tmp1), slp_tmp2 = SeqLocFindNext (slp2, slp_tmp2)) {
      if (SeqLocStart (slp_tmp1) != SeqLocStart (slp_tmp2)
          || SeqLocStop (slp_tmp1) != SeqLocStop (slp_tmp2)
          || (!ignore_partial && !DoLocationPartialsMatch (slp_tmp1, slp_tmp2))) {
        return FALSE;
      }
    }
  } else if (SeqLocCompare (slp1, slp2) != SLC_A_EQ_B) {
    return FALSE;
  }
  return TRUE;
}


static Boolean DoCdRegionsMatch (CdRegionPtr crp1, CdRegionPtr crp2)
{
  if (crp1 == NULL && crp2 == NULL) {
    return TRUE;
  } else if (crp1 == NULL || crp2 == NULL) {
    return FALSE;
  } else if ((crp1->orf && !crp2->orf) || (!crp1->orf && crp2->orf)){
    return FALSE;
  } else if ((crp1->conflict && !crp2->conflict) || (!crp1->conflict && crp2->conflict)){
    return FALSE;
  } else if (crp1->gaps != crp2->gaps) {
    return FALSE;
  } else if (crp1->mismatch != crp2->mismatch) {
    return FALSE;
  } else if (crp1->stops != crp2->stops) {
    return FALSE;
  } else if ((crp1->genetic_code == NULL && crp2->genetic_code != NULL)
             || (crp1->genetic_code != NULL && crp2->genetic_code == NULL)
             || (crp1->genetic_code != NULL && crp2->genetic_code != NULL 
                 && !AsnIoMemComp (crp1->genetic_code, crp2->genetic_code, (AsnWriteFunc) GeneticCodeAsnWrite))) {
    return FALSE;
  } else if ((crp1->code_break == NULL && crp2->code_break != NULL)
             || (crp1->code_break != NULL && crp2->code_break == NULL)
             || (crp1->code_break != NULL && crp2->code_break != NULL 
                 && !AsnIoMemComp (crp1->code_break, crp2->code_break, (AsnWriteFunc) CodeBreakAsnWrite))) {
    return FALSE;
  } else if (crp1->frame != crp2->frame) {
    if ((crp1->frame == 0 || crp1->frame == 1) && (crp2->frame == 0 || crp2->frame == 1)) {
      /* both effectively frame 1, ignore this difference */
    } else {
      return FALSE;
    }
  }
  return TRUE;
}


static Boolean DoesSeqFeatDataMatch (ChoicePtr d1, ChoicePtr d2)
{
  if (d1 == NULL && d2 == NULL) {
    return TRUE;
  } else if (d1 == NULL || d2 == NULL) {
    return FALSE;
  } else if (d1->choice != d2->choice) {
    return FALSE;
  } else if (d1->choice == SEQFEAT_CDREGION) {
    return DoCdRegionsMatch(d1->value.ptrvalue, d2->value.ptrvalue);
  } else {
    return AsnIoMemComp(d1, d2, (AsnWriteFunc) SeqFeatDataAsnWrite);
  }
}


NLM_EXTERN Boolean DoFeaturesMatch (SeqFeatPtr sfp1, SeqFeatPtr sfp2, Boolean allow_different_sequences, Boolean case_sensitive, Boolean ignore_partial)
{
  if (sfp1 == NULL && sfp2 == NULL) {
    return TRUE;
  } else if (sfp1 == NULL || sfp2 == NULL) {
    return FALSE;
  } if (sfp1->data.choice != sfp2->data.choice) {
    return FALSE;
  } else if (sfp1->idx.subtype != sfp2->idx.subtype) {
    return FALSE;
  } else if (!ignore_partial && ((sfp1->partial && !sfp2->partial) || (!sfp1->partial && sfp2->partial))) {
    return FALSE;
  } else if ((sfp1->pseudo && !sfp2->pseudo) || (!sfp1->pseudo && sfp2->pseudo)) {
    return FALSE;
  } else if ((sfp1->excpt && !sfp2->excpt) || (!sfp1->excpt && sfp2->excpt)) {
    return FALSE;
  } else if (!DoLocationsMatch (sfp1->location, sfp2->location, allow_different_sequences, ignore_partial)) {
    return FALSE;
  } else if (!DoStringsMatch (sfp1->comment, sfp2->comment, case_sensitive)) {
    return FALSE;
  } else if (!DoStringsMatch (sfp1->title, sfp2->title, case_sensitive)) {
    return FALSE;
  } else if (sfp1->ext != NULL || sfp2->ext != NULL) {
    return FALSE;
  } else if (sfp1->exts != NULL || sfp2->exts != NULL) {
    return FALSE;
  } else if (!DoStringsMatch (sfp1->except_text, sfp2->except_text, case_sensitive)) {
    return FALSE;
  } else if (sfp1->exp_ev != sfp2->exp_ev) {
    return FALSE;
  } else if (!DoGBQualListsMatch (sfp1->qual, sfp2->qual, case_sensitive)) {
    return FALSE;
  } else if ((sfp1->cit != NULL || sfp2->cit != NULL) && PubMatch (sfp1->cit, sfp2->cit) != 0) {
    return FALSE;
  } else if (!DbxrefsMatch (sfp1->dbxref, sfp2->dbxref, case_sensitive)) {
    return FALSE;
  } else if (!DoesSeqFeatDataMatch(&(sfp1->data), &(sfp2->data))) {
    return FALSE;
  } else if (!XrefsMatch (sfp1->xref, sfp2->xref)) {
    return FALSE;
  } else if (!ProductsMatch (sfp1->product, sfp2->product, case_sensitive, ignore_partial)) {
    return FALSE;
  } else {
    return TRUE;
  }
}


NLM_EXTERN void CleanupStringsForOneDescriptor (SeqDescPtr sdp, SeqEntryPtr sep)
{
  Boolean stripSerial = FALSE;
  Boolean isEmblOrDdbj = FALSE;

  if (sdp == NULL) {
    return;
  }
  SeqEntryExplore (sep, (Pointer) &stripSerial, CheckForSwissProtID);
  SeqEntryExplore (sep, (Pointer) &isEmblOrDdbj, CheckForEmblDdbjID);

  if (sdp->choice == Seq_descr_pub) {
    FlattenPubdesc (sdp->data.ptrvalue, NULL);
  }

  CleanupDescriptorStrings (sdp, stripSerial, TRUE, NULL, isEmblOrDdbj);
}



NLM_EXTERN void CleanupOneSeqFeat (SeqFeatPtr sfp)
{
  Boolean           isEmblOrDdbj = FALSE;
  Boolean           isJscan = FALSE;
  Boolean           stripSerial = TRUE;
  ValNodePtr        publist = NULL;
  SeqEntryPtr       sep;

  if (sfp->idx.entityID == 0) {
    return;
  }
  sep = GetTopSeqEntryForEntityID (sfp->idx.entityID);

  SeqEntryExplore (sep, (Pointer) &stripSerial, CheckForSwissProtID);
  SeqEntryExplore (sep, (Pointer) &isEmblOrDdbj, CheckForEmblDdbjID);
  SeqEntryExplore (sep, (Pointer) &isJscan, CheckForJournalScanID);
  FlattenSfpCit (sfp, NULL);
  CleanUpSeqFeat (sfp, isEmblOrDdbj, isJscan, stripSerial, TRUE, &publist);

  if (publist != NULL) {
   ChangeCitsOnFeats (sfp, publist);
  }
  ValNodeFreeData (publist);
}

/* special cases for chloroplast genetic code until implemented in taxonomy database */

typedef struct pgorg {
  CharPtr  organism;
  Uint1    pgcode;
} PgOrg;

static PgOrg pgOrgList [] = {
  { "Chromera velia", 4 } ,
  { NULL, 0 }
};

typedef struct pglin {
  CharPtr  lineage;
  Uint1    pgcode;
} PgLin;

static PgLin pgLinList [] = {
  { "Eukaryota; Alveolata; Apicomplexa; Coccidia; ", 4 } ,
  { NULL, 0 }
};

NLM_EXTERN Uint1 GetSpecialPlastidGenCode (
  CharPtr taxname,
  CharPtr lineage
)

{
  Int2    i;
  size_t  max;
  Uint1   pgcode = 0;

  if (StringDoesHaveText (taxname)) {
    for (i = 0; pgOrgList [i].organism != NULL; i++) {
      if (StringICmp (taxname, pgOrgList [i].organism) != 0) continue;
      pgcode = pgOrgList [i].pgcode;
    }
  }

  if (StringDoesHaveText (lineage)) {
    for (i = 0; pgLinList [i].lineage != NULL; i++) {
      max = StringLen (pgLinList [i].lineage);
      if (StringNICmp (lineage, pgLinList [i].lineage, max) != 0) continue;
      pgcode = pgLinList [i].pgcode;
    }
  }

  if (pgcode == 11) {
    pgcode = 0;
  }

  return pgcode;
}


static void TrimStopsFromCompleteCodingRegionsCallback (SeqFeatPtr sfp, Pointer data)
{
  Boolean p5, p3;
  BioseqPtr protbsp;
  CharPtr   prot_str;
  Int4      len;
  /* variables for shortening protein features */
  SeqFeatPtr        prot_sfp;
  SeqMgrFeatContext fcontext;
  SeqIntPtr         sintp;
  /* variables for logging */
  LogInfoPtr lip;
  Char      id_buf[100];
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || sfp->product == NULL) {
    return;
  }

  CheckSeqLocForPartial (sfp->location, &p5, &p3);
  if (p3) {
    return;
  }

  protbsp = BioseqFindFromSeqLoc (sfp->product);
  if (protbsp == NULL) {
    return;
  }

  prot_str = GetSequenceByBsp (protbsp);
  if (prot_str == NULL || (len = StringLen (prot_str)) == 0
      || prot_str[len - 1] != '*') {
    prot_str = MemFree (prot_str);
    return;
  }

  BSSeek ((ByteStorePtr) protbsp->seq_data, -1, SEEK_END);
  BSDelete ((ByteStorePtr) protbsp->seq_data, 1);
  protbsp->length -= 1;
  prot_str = MemFree (prot_str);

  for (prot_sfp = SeqMgrGetNextFeature (protbsp, NULL, 0, 0, &fcontext);
       prot_sfp != NULL;
       prot_sfp = SeqMgrGetNextFeature (protbsp, prot_sfp, 0, 0, &fcontext)) {
    if (prot_sfp->location != NULL
        && prot_sfp->location->choice == SEQLOC_INT 
        && (sintp = (SeqIntPtr)prot_sfp->location->data.ptrvalue) != NULL) {
      if (sintp->to > protbsp->length - 1) {
        sintp->to = protbsp->length - 1;
      }
    }
  }

  lip = (LogInfoPtr) data;
  if (lip != NULL) {
    if (lip->fp != NULL) {
      SeqIdWrite (SeqIdFindBest (protbsp->id, SEQID_GENBANK), id_buf, PRINTID_FASTA_SHORT, sizeof (id_buf) - 1);
      fprintf (lip->fp, "Trimmed trailing * from %s\n", id_buf);
    }
    lip->data_in_log = TRUE;
  }
}


NLM_EXTERN Boolean TrimStopsFromCompleteCodingRegions (SeqEntryPtr sep, FILE *log_fp)
{
  LogInfoData lid;
  MemSet (&lid, 0, sizeof (LogInfoData));
  lid.fp = log_fp;
  VisitFeaturesInSep (sep, &lid, TrimStopsFromCompleteCodingRegionsCallback);
  return lid.data_in_log;
}


NLM_EXTERN void 
FixCapitalizationInTitle 
(CharPtr PNTR pTitle,
 Boolean      first_is_upper,
 ValNodePtr   org_names)
{
  if (pTitle == NULL) return;
  ResetCapitalization (first_is_upper, *pTitle);
  FixAbbreviationsInElement (pTitle);
  FixOrgNamesInString (*pTitle, org_names);
}


typedef struct structuredcommentconversion {
  Int4 num_converted;
  Int4 num_unable_to_convert;
} StructuredCommentConversionData, PNTR StructuredCommentConversionPtr;

static void CommentWithSpacesToStructuredCommentCallback (SeqDescPtr sdp, Pointer userdata)
{
  UserObjectPtr uop;
  CharPtr       str, start, stop;
  Int4          len;
  UserFieldPtr  ufp = NULL, prev_ufp = NULL;
  StructuredCommentConversionPtr sd;

  if (sdp == NULL || sdp->choice != Seq_descr_comment || StringHasNoText (sdp->data.ptrvalue)) {
    return;
  }

  uop = UserObjectNew ();
  uop->type = ObjectIdNew ();
  uop->type->str = StringSave ("StructuredComment");

  start = sdp->data.ptrvalue;
  while (*start != 0) {
    stop = start + StringCSpn (start, " ~");
    while (*stop != 0 && *stop != '~' && !isspace (*(stop + 1)) && *(stop + 1) != 0) {
      stop = stop + 1 + StringCSpn (stop + 1, " ~");
    }
    len = 1 + stop - start;
    str = (CharPtr) MemNew (sizeof (Char) * len);
    StringNCpy (str, start, len - 1);
    str[len - 1] = 0;
    if (ufp == NULL) {
      /* add new field */
      ufp = UserFieldNew ();
      if (prev_ufp == NULL) {
        uop->data = ufp;
      } else {
        prev_ufp->next = ufp;
      }
      ufp->label = ObjectIdNew ();
      ufp->label->str = str;
    } else {
      /* add value to last field */
      ufp->choice = 1;
      ufp->data.ptrvalue = str;
      prev_ufp = ufp;
      ufp = NULL;
    }
    if (*stop == 0) {
      start = stop;
    } else {
      start = stop + 1 + StringSpn (stop + 1, " ");
    }
  }

  if (prev_ufp == NULL) {
    uop = UserObjectFree (uop);
    return;
  }
  sd = (StructuredCommentConversionPtr) userdata;
  if (ufp == NULL) {
    sdp->data.ptrvalue = MemFree (sdp->data.ptrvalue);
    sdp->data.ptrvalue = uop;
    sdp->choice = Seq_descr_user;
    if (sd != NULL) {
      sd->num_converted++;
    }
  } else {
    uop = UserObjectFree (uop);
    if (sd != NULL) {
      sd->num_unable_to_convert++;
    }
  }
}


NLM_EXTERN Int4 ConvertCommentsWithSpacesToStructuredCommentsForSeqEntry (SeqEntryPtr sep)
{
  StructuredCommentConversionData sd;
  
  MemSet (&sd, 0, sizeof (StructuredCommentConversionData));
  VisitDescriptorsInSep (sep, &sd, CommentWithSpacesToStructuredCommentCallback);

  return sd.num_unable_to_convert;
}


NLM_EXTERN void RemoveFeatureLink (SeqFeatPtr sfp1, SeqFeatPtr sfp2)
{
  SeqFeatXrefPtr  xref, next, PNTR prevlink;
  ObjectIdPtr     oip;
  SeqFeatPtr      link_sfp;
  Char            buf [32];
  CharPtr         str = NULL;

  if (sfp1 == NULL) return;

  prevlink = (SeqFeatXrefPtr PNTR) &(sfp1->xref);
  xref = sfp1->xref;
  while (xref != NULL) {
    next = xref->next;
    link_sfp = NULL;

    if (xref->id.choice == 3) {
      oip = (ObjectIdPtr) xref->id.value.ptrvalue;
      if (oip != NULL) {
        if (StringDoesHaveText (oip->str)) {
          str = oip->str;
        } else {
          sprintf (buf, "%ld", (long) oip->id);
          str = buf;
        }
        link_sfp = SeqMgrGetFeatureByFeatID (sfp1->idx.entityID, NULL, str, NULL, NULL);
      }
    }
    if (link_sfp == sfp2) {
      *prevlink = xref->next;
      xref->next = NULL;
      MemFree (xref);
    } else {
      prevlink = (SeqFeatXrefPtr PNTR) &(xref->next);
    }

    xref = next;
  }
}


NLM_EXTERN void LinkTwoFeatures (SeqFeatPtr dst, SeqFeatPtr sfp)

{
  ChoicePtr       cp;
  ObjectIdPtr     oip;
  SeqFeatXrefPtr  xref, prev_xref, next_xref;
  SeqFeatPtr      old_match;

  if (dst == NULL || sfp == NULL) return;

  cp = &(dst->id);
  if (cp == NULL) return;
  if (cp->choice == 3) {
    /* don't create a duplicate xref, remove links to other features */
    xref = sfp->xref;
    prev_xref = NULL;
    while (xref != NULL) {
      next_xref = xref->next;
      if (xref->id.choice == 3 && xref->id.value.ptrvalue != NULL) {
        if (ObjectIdMatch (cp->value.ptrvalue, xref->id.value.ptrvalue)) {
          /* already have this xref */
          return;
        } else {
          old_match = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, NULL, xref, NULL);
          RemoveFeatureLink (sfp, old_match);
          RemoveFeatureLink (old_match, sfp);
        }
      } else {
        prev_xref = xref;
      }
      xref = next_xref;
    }

    oip = (ObjectIdPtr) cp->value.ptrvalue;
    if (oip != NULL) {
      oip = AsnIoMemCopy (oip, (AsnReadFunc) ObjectIdAsnRead,
                          (AsnWriteFunc) ObjectIdAsnWrite);
      if (oip != NULL) {
        xref = SeqFeatXrefNew ();
        if (xref != NULL) {
          xref->id.choice = 3;
          xref->id.value.ptrvalue = (Pointer) oip;
          xref->next = sfp->xref;
          sfp->xref = xref;
        }
      }
    }
  }
}


static void MakeFeatureXrefsFromProteinIdQualsCallback (SeqFeatPtr sfp, Pointer data)
{
  GBQualPtr gbq;
  SeqIdPtr sip;
  BioseqPtr pbsp;
  SeqFeatPtr cds;
  CharPtr    product;
  ProtRefPtr prp;
  SeqEntryPtr sep;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_mRNA || (sep = (SeqEntryPtr) data) == NULL) {
    return;
  }

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "protein_id") == 0 || StringICmp (gbq->qual, "orig_protein_id") == 0) {
      sip = CreateSeqIdFromText (gbq->val, sep);
      pbsp = BioseqFind (sip);
      cds = SeqMgrGetCDSgivenProduct (pbsp, NULL);
      if (cds != NULL) {
        LinkTwoFeatures (cds, sfp);
        LinkTwoFeatures (sfp, cds);
        product = GetRNAProductString(sfp, NULL);
        if (StringHasNoText (product)) {
          prp = GetProtRefForFeature (cds);
          if (prp != NULL && prp->name != NULL && !StringHasNoText (prp->name->data.ptrvalue)) {
            SetRNAProductString (sfp, NULL, prp->name->data.ptrvalue, ExistingTextOption_replace_old);
          }
        }
        product = MemFree (product);
      }
    }
  }
}


NLM_EXTERN void MakeFeatureXrefsFromProteinIdQuals (SeqEntryPtr sep)
{
  /* assign feature IDs, so that we can create xrefs that use them */
  AssignFeatureIDs (sep);

  VisitFeaturesInSep (sep, (Pointer) sep, MakeFeatureXrefsFromProteinIdQualsCallback);
}


static void MakeFeatureXrefsFromTranscriptIdQualsCallback (SeqFeatPtr sfp, Pointer data)
{
  GBQualPtr gbq;
  SeqIdPtr sip;
  BioseqPtr pbsp;
  SeqFeatPtr cds;
  CharPtr    product;
  ProtRefPtr prp;
  SeqEntryPtr sep;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_mRNA || (sep = (SeqEntryPtr) data) == NULL) {
    return;
  }

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "transcript_id") == 0 || StringICmp (gbq->qual, "orig_transcript_id") == 0) {
      sip = CreateSeqIdFromText (gbq->val, sep);
      pbsp = BioseqFind (sip);
      cds = SeqMgrGetCDSgivenProduct (pbsp, NULL);
      if (cds != NULL) {
        LinkTwoFeatures (cds, sfp);
        LinkTwoFeatures (sfp, cds);
        product = GetRNAProductString(sfp, NULL);
        if (StringHasNoText (product)) {
          prp = GetProtRefForFeature (cds);
          if (prp != NULL && prp->name != NULL && !StringHasNoText (prp->name->data.ptrvalue)) {
            SetRNAProductString (sfp, NULL, prp->name->data.ptrvalue, ExistingTextOption_replace_old);
          }
        }
        product = MemFree (product);
      }
    }
  }
}


NLM_EXTERN void MakeFeatureXrefsFromTranscriptIdQuals (SeqEntryPtr sep)
{
  /* assign feature IDs, so that we can create xrefs that use them */
  AssignFeatureIDs (sep);

  VisitFeaturesInSep (sep, (Pointer) sep, MakeFeatureXrefsFromTranscriptIdQualsCallback);
}


static void FinishHalfXrefsCallback (SeqFeatPtr sfp, Pointer data)
{
  SeqFeatPtr other;
  SeqFeatXrefPtr xref, xref_other;
  Boolean has_other_xref;

  if (sfp == NULL) {
    return;
  }

  xref = sfp->xref;
  while (xref != NULL) {
    if (xref->id.choice == 3) {
      other = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, NULL, xref, NULL);
      if (other != NULL) {
        xref_other = other->xref;
        has_other_xref = FALSE;
        while (xref_other != NULL && !has_other_xref) {
          if (xref_other->id.choice == 3) {
            has_other_xref = TRUE;
          }
          xref_other = xref_other->next;
        }
        if (!has_other_xref) {
          LinkTwoFeatures (sfp, other);
        }
      }
    }
    xref = xref->next;
  }
}


NLM_EXTERN void FinishHalfXrefs (SeqEntryPtr sep)
{
  VisitFeaturesInSep (sep, (Pointer) sep, FinishHalfXrefsCallback);
}


NLM_EXTERN Uint1 GetAaFromtRNA (tRNAPtr trp)
{
  Uint1 aa;
  Uint1 from;
  SeqMapTablePtr  smtp;

  if (trp == NULL) {
    return 0;
  }

  aa = 0;
  if (trp->aatype == 2) {
    aa = trp->aa;
  } else {
    from = 0;
    switch (trp->aatype) {
    case 0:
      from = 0;
      break;
    case 1:
      from = Seq_code_iupacaa;
      break;
    case 2:
      from = Seq_code_ncbieaa;
      break;
    case 3:
      from = Seq_code_ncbi8aa;
      break;
    case 4:
      from = Seq_code_ncbistdaa;
      break;
    default:
      break;
    }
    smtp = SeqMapTableFind (Seq_code_ncbieaa, from);
    if (smtp != NULL) {
      aa = SeqMapTableConvert (smtp, trp->aa);
    }
  }
  return aa;
}


NLM_EXTERN CharPtr GetCodesFortRNA (SeqFeatPtr sfp, Int2 *pCode)
{
  BioseqPtr       bsp;
  Int2            code = 0;
  GeneticCodePtr  gncp;
  ValNodePtr      vnp;
  CharPtr         codes = NULL;

  if (sfp == NULL) {
    return NULL;
  }

  /* find genetic code table */

  bsp = GetBioseqGivenSeqLoc (sfp->location, sfp->idx.entityID);
  BioseqToGeneticCode (bsp, &code, NULL, NULL, NULL, 0, NULL);

  gncp = GeneticCodeFind (code, NULL);
  if (gncp == NULL) {
    gncp = GeneticCodeFind (1, NULL);
    code = 1;
  }
  if (gncp != NULL) {
    for (vnp = (ValNodePtr) gncp->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice != 3) continue;
      codes = (CharPtr) vnp->data.ptrvalue;
      break;
    }
  }
  if (pCode != NULL) {
    *pCode = code;
  }
  return codes;
}


static Boolean DoesCodonMatchAminoAcid (Uint1 aa, Uint1 index, CharPtr codes)
{
  Uint1           taa;
  Boolean         rval = FALSE;
    
  if (aa == 0 || aa == 255 || codes == NULL) 
  {
    return TRUE;
  }
  taa = codes [index];

  if (taa == aa) 
  {
    rval = TRUE;
  }
  /* selenocysteine normally uses TGA (14), so ignore without requiring exception in record */
  else if (aa == 'U' && taa == '*' && index == 14) 
  {
    rval = TRUE;
  }
  /* pyrrolysine normally uses TAG (11) in archaebacteria, ignore without requiring exception */
  else if (aa == 'O' && taa == '*' && index == 11) {
    rval = TRUE;
  }
  /* TAA (10) is not yet known to be used for an exceptional amino acid, but the night is young */

  return rval;
}


static Boolean IsATGC (Char ch)
{
  if (ch == 'A' || ch == 'T' || ch == 'G' || ch == 'C') {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Char s_comp (Char ch)
{
  if (ch == 'A') {
    return 'T';
  } else if (ch == 'G') {
    return 'C';
  } else if (ch == 'C') {
    return 'G';
  } else if (ch == 'T') {
    return 'A';
  } else {
    return 'N';
  }
}


static CharPtr GetFlipCodonLoggingInfo (SeqFeatPtr sfp)
{
  SeqFeatPtr gene = NULL;
  GeneRefPtr grp = NULL;
  ValNode   vn;
  CharPtr txt = NULL;

  GetGeneInfoForFeature (sfp, &grp, &gene);
  if (grp != NULL && !StringHasNoText (grp->locus_tag)) {
    txt = StringSave (grp->locus_tag);
  } else {
    MemSet (&vn, 0, sizeof (ValNode));
    vn.choice = OBJ_SEQFEAT;
    vn.data.ptrvalue = sfp;
    txt = GetDiscrepancyItemText (&vn);
  }
  return txt;
}


static Int4 CountCodonsRecognized (tRNAPtr trp)
{
  Int4 num = 0, i;

  if (trp == NULL) {
    return 0;
  }
  for (i = 0; i < 6; i++) {
    if (trp->codon [i] < 64) {
      num++;
    }
  }
  return num;
}


static Int4 CountMatchingCodons (tRNAPtr trp, Uint1 aa, CharPtr codes)
{
  Int4 num = 0, i;

  if (trp == NULL) {
    return 0;
  }
  for (i = 0; i < 6; i++) {
    if (trp->codon [i] < 64) {
      if (DoesCodonMatchAminoAcid (aa, trp->codon[i], codes)) {
        num++;
      }
    }
  }

  return num;
}


static Int4 CountFlippableCodons (tRNAPtr trp, Uint1 aa, CharPtr codes, Int2 code)
{
  Int4 num = 0, i;
  Int2      index;
  Uint1     codon [4];
  Uint1     rcodon [4];

  if (trp == NULL) {
    return 0;
  }
	/* Note - it is important to set the fourth character in the codon array to NULL
		* because CodonForIndex only fills in the three characters of actual codon,
		* so if you StringCpy the codon array and the NULL character is not found after
		* the three codon characters, you will write in memory you did not intend to.
		*/
	codon [3] = 0;
  rcodon [3] = 0;
  for (i = 0; i < 6; i++) 
  {
    if (trp->codon [i] < 64 
        && !DoesCodonMatchAminoAcid (aa, trp->codon[i], codes) 
        && CodonForIndex (trp->codon [i], Seq_code_iupacna, codon)
        && IsATGC(codon[0])
        && IsATGC(codon[1])
        && IsATGC(codon[2])) 
    {
      rcodon[0] = s_comp(codon[2]);
      rcodon[1] = s_comp(codon[1]);
      rcodon[2] = s_comp(codon[0]);
      index = IndexForCodon (rcodon, code);
      if (index < 64 && DoesCodonMatchAminoAcid(aa, index, codes))
      {
        num++;
      }
    }
  }

  return num;
}


static Int4 FlipFlippableCodons (tRNAPtr trp, Uint1 aa, CharPtr codes, Int2 code)
{
  Int4 num = 0, i;
  Int2      index;
  Uint1     codon [4];
  Uint1     rcodon [4];

  if (trp == NULL) {
    return 0;
  }
	/* Note - it is important to set the fourth character in the codon array to NULL
		* because CodonForIndex only fills in the three characters of actual codon,
		* so if you StringCpy the codon array and the NULL character is not found after
		* the three codon characters, you will write in memory you did not intend to.
		*/
	codon [3] = 0;
  rcodon [3] = 0;
  for (i = 0; i < 6; i++) 
  {
    if (trp->codon [i] < 64 
        && !DoesCodonMatchAminoAcid (aa, trp->codon[i], codes) 
        && CodonForIndex (trp->codon [i], Seq_code_iupacna, codon)
        && IsATGC(codon[0])
        && IsATGC(codon[1])
        && IsATGC(codon[2])) 
    {
      rcodon[0] = s_comp(codon[2]);
      rcodon[1] = s_comp(codon[1]);
      rcodon[2] = s_comp(codon[0]);
      index = IndexForCodon (rcodon, code);
      if (index < 64 && DoesCodonMatchAminoAcid(aa, index, codes))
      {
        trp->codon[i] = index;
        num++;
      }
    }
  }

  return num;
}


static Boolean IgnoretRNACodonRecognized (SeqFeatPtr sfp)
{
  if (sfp == NULL 
      || StringISearch (sfp->except_text, "RNA editing") != NULL
      || StringISearch (sfp->except_text, "modified codon recognition") != NULL)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static void FlipCodonRecognizedCallback (SeqFeatPtr sfp, Pointer data)
{
  RnaRefPtr rrp;
  tRNAPtr   trp;
  Uint1     aa;  
  CharPtr   txt;
  LogInfoPtr lip;
  Int2       code = 0;
  CharPtr    codes = NULL;
  Int4       num_codons, num_match, num_flippable;

  if (IgnoretRNACodonRecognized(sfp)
      || sfp->idx.subtype != FEATDEF_tRNA
      || (rrp = (RnaRefPtr) sfp->data.value.ptrvalue) == NULL
      || rrp->ext.choice != 2
      || (trp = (tRNAPtr)(rrp->ext.value.ptrvalue)) == NULL)
  {
    return;
  }

  num_codons = CountCodonsRecognized (trp);
  if (num_codons == 0) {
    return;
  }

  lip = (LogInfoPtr) data;

  aa = GetAaFromtRNA (trp);

  /* find genetic code table */
  codes = GetCodesFortRNA (sfp, &code);

  if (codes == NULL) return;

  num_match = CountMatchingCodons (trp, aa, codes);
  if (num_codons == num_match) {
    return;
  } else if (num_codons > 1) {
    if (lip != NULL) 
    {
      if (lip->fp != NULL) 
      {
        /* text for log */
        txt = GetFlipCodonLoggingInfo (sfp);
        fprintf (lip->fp, "Unable to flip bad codon_recognized for %s\n", txt);
        txt = MemFree (txt);
      }
      lip->data_in_log = TRUE;
    }
  } else {
    num_flippable = CountFlippableCodons(trp, aa, codes, code);
    if (num_flippable == num_codons) {
      FlipFlippableCodons (trp, aa, codes, code);
    } else {
      if (lip != NULL) 
      {
        if (lip->fp != NULL) 
        {
          /* text for log */
          txt = GetFlipCodonLoggingInfo (sfp);
          fprintf (lip->fp, "Unable to flip bad codon_recognized for %s\n", txt);
          txt = MemFree (txt);
        }
        lip->data_in_log = TRUE;
      }
    }
  }
}


NLM_EXTERN void FlipCodonRecognizedInSeqEntry (SeqEntryPtr sep, LogInfoPtr lip)
{
  VisitFeaturesInSep (sep, lip, FlipCodonRecognizedCallback);
}


static void RemoveBadCodonRecognizedCallback (SeqFeatPtr sfp, Pointer data)
{
  RnaRefPtr rrp;
  tRNAPtr   trp;
  Int2      j, k;
  Uint1     aa;  
  Uint1     codon [4];
  Uint1     rcodon [4];
  CharPtr   txt;
  LogInfoPtr lip;
  Int2       code = 0;
  CharPtr    codes = NULL;
  Int4       num_codons, num_match;

  if (IgnoretRNACodonRecognized(sfp)
      || sfp->idx.subtype != FEATDEF_tRNA
      || (rrp = (RnaRefPtr) sfp->data.value.ptrvalue) == NULL
      || rrp->ext.choice != 2
      || (trp = (tRNAPtr)(rrp->ext.value.ptrvalue)) == NULL)
  {
    return;
  }

  num_codons = CountCodonsRecognized (trp);
  if (num_codons == 0) {
    return;
  }

  lip = (LogInfoPtr) data;

  aa = GetAaFromtRNA (trp);

  /* find genetic code table */
  codes = GetCodesFortRNA (sfp, &code);

  if (codes == NULL) return;

  num_match = CountMatchingCodons (trp, aa, codes);
  if (num_match == num_codons) {
    return;
  }

	/* Note - it is important to set the fourth character in the codon array to NULL
		* because CodonForIndex only fills in the three characters of actual codon,
		* so if you StringCpy the codon array and the NULL character is not found after
		* the three codon characters, you will write in memory you did not intend to.
		*/
	codon [3] = 0;
  rcodon [3] = 0;

  for (j = 0; j < 6; j++) 
  {
    if (trp->codon [j] < 64) 
    {
      if (DoesCodonMatchAminoAcid (aa, trp->codon[j], codes)) 
      {
        /* already ok - skip it */
      } 
      else if (CodonForIndex (trp->codon [j], Seq_code_iupacna, codon)
          && IsATGC(codon[0])
          && IsATGC(codon[1])
          && IsATGC(codon[2])) 
      {
        for (k = j + 1; k < 6; k++) 
        {
          trp->codon[k - 1] = trp->codon[k];
        }
        trp->codon[5] = 255;
        if (lip != NULL) 
        {
          if (lip->fp != NULL) 
          {
            /* text for log */
            txt = GetFlipCodonLoggingInfo (sfp);
            fprintf (lip->fp, "Removed codon_recognized '%s' for %s\n", codon, txt);
            txt = MemFree (txt);
          }
          lip->data_in_log = TRUE;
        }
        /* push index down, so we don't skip over a codon */
        j--;
      }
    }
  }
}


NLM_EXTERN void RemoveBadCodonRecognizedInSeqEntry (SeqEntryPtr sep, LogInfoPtr lip)
{
  VisitFeaturesInSep (sep, lip, RemoveBadCodonRecognizedCallback);
}


NLM_EXTERN void ReverseBioseqInAlignment (SeqAlignPtr salp, Pointer userdata)
{
  BioseqPtr bsp;
  SeqIdPtr  sip;
  Boolean   found = FALSE;
  Int4      order;
  
  if (salp == NULL || userdata == NULL) return;
  
  bsp = (BioseqPtr) userdata;
  
  for (sip = bsp->id; sip != NULL && ! found; sip = sip->next) 
  {
    order = SeqIdOrderInBioseqIdList(sip, SeqIdPtrFromSeqAlign (salp));
    if (order > 0) {
      AlnMgr2IndexSeqAlignEx(salp, FALSE);
      ReverseAlignmentStrand (salp, order);
      SeqAlignIndexFree(salp->saip);
      salp->saip = NULL;
      found = TRUE;
    }
  }
}


/* need to reverse the order of the segments and flip the strands */
NLM_EXTERN void FlipAlignment (SeqAlignPtr salp)
{
  DenseSegPtr dsp;
  Int4        row, seg, swap_start, swap_len, opp_seg;
  Score    swap_score;
  Uint1       swap_strand;
  
  if (salp == NULL || salp->segtype != SAS_DENSEG || salp->segs == NULL)
  {
    return;
  }
  
  dsp = (DenseSegPtr) salp->segs;
  if (dsp->strands == NULL) {
    dsp->strands = (Uint1Ptr) MemNew (dsp->numseg * dsp->dim * sizeof (Uint1));
    MemSet (dsp->strands, Seq_strand_plus, dsp->numseg * dsp->dim * sizeof (Uint1));
  }

  for (seg = 0; seg < dsp->numseg / 2; seg++) {
    /* swap segments to reverse order */
    opp_seg = dsp->numseg - 1 - seg;
    /* swap lens */
    swap_len = dsp->lens[seg];
    dsp->lens[seg] = dsp->lens[opp_seg];
    dsp->lens[opp_seg] = swap_len;
    /* swap scores */
    if (dsp->scores != NULL) {
      swap_score = dsp->scores[seg];
      dsp->scores[seg] = dsp->scores[opp_seg];
      dsp->scores[opp_seg] = swap_score;
    }
    for (row = 0; row < dsp->dim; row++) {
      /* swap strands */
      swap_strand = dsp->strands[dsp->dim * seg + row];
      dsp->strands[dsp->dim * seg + row] = dsp->strands[dsp->dim * opp_seg + row];
      dsp->strands[dsp->dim * opp_seg + row] = swap_strand;
      
      /* swap starts */
      swap_start = dsp->starts[dsp->dim * seg + row];
      dsp->starts[dsp->dim * seg + row] = dsp->starts[dsp->dim * opp_seg + row];
      dsp->starts[dsp->dim * opp_seg + row] = swap_start;      
    }
  }
  
  /* reverse segments */
  for (seg = 0; seg < dsp->numseg; seg++) {
    for (row = 0; row < dsp->dim; row++) {
      if (dsp->strands[dsp->dim * seg + row] == Seq_strand_minus) {
        dsp->strands[dsp->dim * seg + row] = Seq_strand_plus;
      } else {
        dsp->strands[dsp->dim * seg + row] = Seq_strand_minus;
      }
    }
  }
  SAIndex2Free2(salp->saip);
  salp->saip = NULL;
}


NLM_EXTERN void FlipEntireAlignmentIfAllSequencesFlipped (SeqAnnotPtr sap, Pointer userdata)
{
  SeqAlignPtr salp;
  ValNodePtr  vnp;
  BioseqPtr   bsp;
  SeqIdPtr    sip;
  Boolean     found;
  Int4 row, num_rows;

  if (sap == NULL || sap->type != 2 || userdata == NULL) return;
  salp = (SeqAlignPtr) sap->data;
  if (salp == NULL || salp->idx.deleteme) return;
  
  
  AlnMgr2IndexSingleChildSeqAlign(salp);
  num_rows = AlnMgr2GetNumRows(salp);
  for (row = 1; row <= num_rows; row++) {
    sip = AlnMgr2GetNthSeqIdPtr(salp, row);
    found = FALSE;
    vnp = (ValNodePtr)userdata;
    while (vnp != NULL && !found) {
      bsp = (BioseqPtr) vnp->data.ptrvalue;
      if (SeqIdOrderInBioseqIdList (sip, bsp->id) > 0) {
        found = TRUE;
      }
      vnp = vnp->next;
    }
    if (!found) return;
  }
  
  FlipAlignment(salp);      
}


NLM_EXTERN ValNodePtr ListSequencesWithAlignments (ValNodePtr bsp_list)
{
  BioseqPtr     bsp;
  ValNodePtr    vnp, aln_bsp = NULL;

  for (vnp = bsp_list; vnp != NULL; vnp = vnp->next) {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (bsp != NULL && IsBioseqInAnyAlignment (bsp, bsp->idx.entityID)) {
      ValNodeAddPointer (&aln_bsp, 0, bsp);
    }
  }
  return aln_bsp;
}


NLM_EXTERN void RevCompBioseqList (ValNodePtr bsp_list, 
                                   Uint2 entityID, 
                                   BioseqFunc func,
                                   Boolean revCompFeats,
                                   Boolean check_for_aln)
{
  SeqEntryPtr sep;
  BioseqPtr   bsp;
  ValNodePtr  vnp;

  sep = GetTopSeqEntryForEntityID (entityID);
  
  for (vnp = bsp_list; vnp != NULL; vnp = vnp->next) {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (func != NULL) {
      func (bsp);
      if (check_for_aln) {
        VisitAlignmentsInSep (sep, (Pointer) bsp, ReverseBioseqInAlignment);
      }
    }
    if (revCompFeats) {
      if (bsp->repr == Seq_repr_raw || bsp->repr == Seq_repr_const) {
        
        if (sep != NULL) {
          SeqEntryExplore (sep, (Pointer) bsp, RevCompFeats);
        }
      }
    }
  }
}


typedef struct bioseqinalignmentdata {
	Boolean   found;
	BioseqPtr lookingfor;
} BioseqInAlignmentData, PNTR BioseqInAlignmentPtr;

static Boolean IsBioseqInThisAlignment (SeqAlignPtr salp, BioseqPtr bsp)
{
  SeqIdPtr sip;
  Boolean found = FALSE;

  for (sip = bsp->id; sip != NULL && ! found; sip = sip->next) 
  {
    found = SeqAlignFindSeqId (salp, sip);
  }
  return found;
}

static void FindAlignmentCallback (SeqAnnotPtr sap, Pointer userdata)
{
  BioseqInAlignmentPtr biap;
  SeqAlignPtr          salp;

  if (sap == NULL || sap->type != 2 || userdata == NULL) 
  {
    return;
  }
  biap = (BioseqInAlignmentPtr) userdata;
  if (biap->found) return;
  salp = (SeqAlignPtr) sap->data;
  if (salp == NULL) return;
  biap->found = IsBioseqInThisAlignment (salp, biap->lookingfor);

}

NLM_EXTERN Boolean IsBioseqInAnyAlignment (BioseqPtr bsp, Uint2 input_entityID)
{
  SeqEntryPtr           topsep;
  BioseqInAlignmentData biad;

  topsep = GetTopSeqEntryForEntityID (input_entityID);
  biad.found = FALSE;
  biad.lookingfor = bsp;

  VisitAnnotsInSep (topsep, &biad, FindAlignmentCallback);
  return biad.found;
}


static void RemoveAlignmentsWithSequenceCallback (SeqAnnotPtr sap, Pointer userdata)
{
  SeqAlignPtr salp;
  SeqIdPtr    sip;

  if (sap == NULL || sap->type != 2 || userdata == NULL) return;
  salp = (SeqAlignPtr) sap->data;
  if (salp == NULL || salp->idx.deleteme) return;
  sip = (SeqIdPtr) userdata;
  while (sip != NULL && !sap->idx.deleteme) {
    if (FindSeqIdinSeqAlign (salp, sip)) {
  	  sap->idx.deleteme = TRUE;
  	}
  	sip = sip->next;
  }
}

NLM_EXTERN void RemoveAlignmentsWithSequence (BioseqPtr bsp, Uint2 input_entityID)
{
  SeqEntryPtr           topsep;

  if (bsp == NULL) return;
  topsep = GetTopSeqEntryForEntityID (input_entityID);

  VisitAnnotsInSep (topsep, bsp->id, RemoveAlignmentsWithSequenceCallback);
}


/* assumes locations on same Bioseq */
static Boolean OutOfOrder (SeqLocPtr slp_prev, SeqLocPtr slp_next)
{
  Uint1 strand_p, strand_n;
  Boolean rval = FALSE;
  Int4 start_p, start_n, stop_p, stop_n;

  if (slp_prev == NULL || slp_next == NULL) 
  {
    return FALSE;
  }

  strand_p = SeqLocStrand (slp_prev);
  strand_n = SeqLocStrand (slp_next);
  if (strand_p == Seq_strand_minus) 
  {
    if (strand_n != Seq_strand_minus) 
    {
      /* mixed strand, not necessarily out of order */
      rval = FALSE;
    } else {
      start_p = SeqLocStart (slp_prev);
      stop_p = SeqLocStop (slp_prev);
      start_n = SeqLocStart (slp_next);
      stop_n = SeqLocStop (slp_next);
      if (start_p < start_n || stop_p < stop_n) 
      {
        rval = TRUE;
      }
    }
  } else {
    if (strand_n == Seq_strand_minus)
    {
      /* mixed strand, not necessarily out of order */
      rval = FALSE;
    } else {
      start_p = SeqLocStart (slp_prev);
      stop_p = SeqLocStop (slp_prev);
      start_n = SeqLocStart (slp_next);
      stop_n = SeqLocStop (slp_next);
      if (start_p > start_n || stop_p > stop_n) 
      {
        rval = TRUE;
      }
    }
  }
  return rval;
}


/* assumes locations on same Bioseq and in order on same strand*/
static Boolean TooFarApartForTransSplicing (SeqLocPtr slp_prev, SeqLocPtr slp_next)
{
  Boolean rval = FALSE;
  Int4 start_n, start_p, stop_n, stop_p;

  if (slp_prev == NULL || slp_next == NULL) 
  {
    return FALSE;
  }

  if (SeqLocStrand (slp_prev) == Seq_strand_minus) 
  {
    start_p = SeqLocStart (slp_prev);
    stop_n = SeqLocStop (slp_next);
    if (start_p - stop_n > 10000) 
    {
      rval = TRUE;
    }
  } else {
    stop_p = SeqLocStop (slp_prev);
    start_n = SeqLocStart (slp_next);
    if (start_n - stop_p > 10000) 
    {
      rval = TRUE;
    }
  }
  return rval;
}


NLM_EXTERN SeqLocPtr MakeGeneLocForFeatureLoc (SeqLocPtr floc, Uint2 entityID, Boolean trans_spliced)
{
  /* in the age of small-set genomes, we're going to pretend that segmented sets do not exist.
   * A gene location for a feature location that includes multiple bioseqs should include
   * one interval per bioseq that covers all locations of the feature that occur on that bioseq.
   */

  SeqLocPtr slp_new = NULL, slp_tmp, slp_last = NULL, add_slp;
  SeqLocPtr PNTR pAddSlp = NULL;
  BioseqPtr bsp, last_bsp = NULL;
  Boolean partial5 = FALSE, partial3 = FALSE;
  Uint2   strand, last_strand = Seq_strand_plus;

  pAddSlp = &slp_new;
  for (slp_tmp = SeqLocFindNext (floc, NULL);
       slp_tmp != NULL;
       slp_tmp = SeqLocFindNext (floc, slp_tmp))
  {
    bsp = GetBioseqGivenSeqLoc (slp_tmp, entityID);
    strand = SeqLocStrand (slp_tmp);
    if (bsp != last_bsp || strand != last_strand 
        || (trans_spliced && OutOfOrder (slp_last, slp_tmp))
        || (trans_spliced && TooFarApartForTransSplicing(slp_last, slp_tmp))) {
      add_slp = SeqLocMerge (bsp, slp_tmp, NULL, TRUE, FALSE, FALSE);
      if (slp_last == NULL) {
        slp_new = add_slp;
      } else {
        slp_last->next = add_slp;
        pAddSlp = &(slp_last->next);
      }
      slp_last = add_slp;
      last_bsp = bsp;
      last_strand = strand;
    } else {
      add_slp = SeqLocMerge (bsp, *pAddSlp, slp_tmp, TRUE, FALSE, FALSE);
      *pAddSlp = SeqLocFree (*pAddSlp);
      *pAddSlp = add_slp;
      slp_last = add_slp;
    }
  }
  if (slp_new != NULL && slp_new->next != NULL) {
    slp_tmp = ValNodeNew (NULL);
    slp_tmp->choice = SEQLOC_MIX;
    slp_tmp->data.ptrvalue = slp_new;
    slp_new = slp_tmp;
  }
  if (slp_new != NULL) {
    CheckSeqLocForPartial (floc, &partial5, &partial3);
    SetSeqLocPartial (slp_new, partial5, partial3);
  }

  return slp_new;
}


/* code for resolving conflicting IDs */
typedef struct {
  CharPtr  oldStr;
  SeqIdPtr newSip;
} ReplaceIDStruct, PNTR ReplaceIDStructPtr;


/********************************************************************
*
* SeqLocReplaceLocalID
*   replaces the Seq-Id in a Seq-Loc (slp) with a new Seq-Id (new_sip)
*   only if the Seq-Id is a local one.
*
**********************************************************************/

static SeqLocPtr SeqLocReplaceLocalID (SeqLocPtr slp,
				       SeqIdPtr  new_sip)
{
  SeqLocPtr        curr;
  PackSeqPntPtr    pspp;
  SeqIntPtr        target_sit;
  SeqPntPtr        spp;
  SeqIdPtr         currId;

  switch (slp->choice) {
     case SEQLOC_PACKED_INT :
     case SEQLOC_MIX :
     case SEQLOC_EQUIV :
        curr = NULL;
        while ((curr = SeqLocFindNext (slp, curr)) != NULL) {
           curr = SeqLocReplaceLocalID (curr, new_sip);
        }
        break;
     case SEQLOC_PACKED_PNT :
        pspp = (PackSeqPntPtr) slp->data.ptrvalue;
        if ((pspp != NULL) && (pspp->id->choice == SEQID_LOCAL)) {
          SeqIdFree (pspp->id);
          pspp->id = SeqIdDup (new_sip);
        }
        break;
     case SEQLOC_EMPTY :
     case SEQLOC_WHOLE :
        currId = (SeqIdPtr) slp->data.ptrvalue;
	if (currId->choice == SEQID_LOCAL)
	  {
	    SeqIdFree (currId);
	    slp->data.ptrvalue = (Pointer) SeqIdDup (new_sip);
	  }
        break;
     case SEQLOC_INT :
        target_sit = (SeqIntPtr) slp->data.ptrvalue;
	if (target_sit->id->choice == SEQID_LOCAL)
	  {
	    SeqIdFree (target_sit->id);
	    target_sit->id = SeqIdDup (new_sip);
	  }
        break;
     case SEQLOC_PNT :
        spp = (SeqPntPtr)slp->data.ptrvalue;
	if (spp->id->choice == SEQID_LOCAL)
	  {
	    SeqIdFree(spp->id);
	    spp->id = SeqIdDup(new_sip);
	  }
        break;
     default :
        break;
  }
  return slp;
}

static void ReplaceIdForFeature (SeqFeatPtr sfp, SeqIdPtr sip)
{
  CdRegionPtr  crp;
  CodeBreakPtr cbp;
  RnaRefPtr    rrp;
  tRNAPtr      trp;

  if (sfp == NULL || sip == NULL) {
    return;
  }
  /* replace local ID in location */
  if (sfp->location != NULL) {
    SeqLocReplaceLocalID (sfp->location, sip);
  }

  /* also replace local ID in code breaks */
  if (sfp->data.choice == SEQFEAT_CDREGION
      && (crp = (CdRegionPtr)sfp->data.value.ptrvalue) != NULL
      && crp->code_break != NULL) {
    for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
      SeqLocReplaceLocalID (cbp->loc, sip);
    }
  }

  /* also replace local ID in anticodons */
  if (sfp->data.choice == SEQFEAT_RNA
      && (rrp = (RnaRefPtr) sfp->data.value.ptrvalue) != NULL
      && rrp->type == 3 && rrp->ext.choice == 2
      && (trp = (tRNAPtr) rrp->ext.value.ptrvalue) != NULL
      && trp->anticodon != NULL) {
    SeqLocReplaceLocalID (trp->anticodon, sip);
  }
}


static void ReplaceLocalIdOnLoc_callback (SeqFeatPtr sfp, Pointer userdata)
{
  SeqIdPtr     sip;

  if (sfp == NULL) {
    return;
  }

  sip = (SeqIdPtr) userdata;
  ReplaceIdForFeature (sfp, sip);
}


static void CheckFeatForNuclID_callback (SeqFeatPtr sfp, Pointer userdata)
{
  SeqIdPtr            featSip = NULL;
  ReplaceIDStructPtr  idsPtr;
  ObjectIdPtr         oip;
  Char                tmpIdStr [128];

  if (NULL == sfp)
    return;

  /* Get the old Seq Id and the new */
  /* one that it was changed to.    */
    
  idsPtr = (ReplaceIDStructPtr) userdata;
  if ((NULL == idsPtr)         ||
      (NULL == idsPtr->oldStr) ||
      (NULL == idsPtr->newSip))
    return;

  /* Get the location Seq ID for this CDS feature */
    
  featSip = SeqLocId (sfp->location);
  if (featSip == NULL) return;
  oip     = (ObjectIdPtr) featSip->data.ptrvalue;
    
  /* If the location Seq ID matches the old Seq Id */
  /* then change the location to point to the new. */
    
  if (NULL == oip->str) {
    sprintf (tmpIdStr, "%d", oip->id);
    if (StringCmp (tmpIdStr, idsPtr->oldStr) == 0) {
      ReplaceIdForFeature (sfp, idsPtr->newSip);
    }
  } else if (StringCmp (oip->str, idsPtr->oldStr) == 0){
    ReplaceIdForFeature (sfp, idsPtr->newSip);
  }
}


static void CheckFeatForProductID_callback (SeqFeatPtr sfp, Pointer userdata)
{
  SeqIdPtr            featSip = NULL;
  ReplaceIDStructPtr  idsPtr;
  ObjectIdPtr         oip;
  Char                tmpIdStr [128];

  if (NULL == sfp)
    return;

  if ((sfp->data.choice == SEQFEAT_CDREGION) &&
      (sfp->product != NULL)) {

    /* Get the old Seq Id and the new */
    /* one that it was changed to.    */
    
    idsPtr = (ReplaceIDStructPtr) userdata;
    if ((NULL == idsPtr)         ||
	(NULL == idsPtr->oldStr) ||
	(NULL == idsPtr->newSip))
      return;

    /* Get the product Seq ID for this CDS feature */
    
    featSip = SeqLocId (sfp->product);
    oip     = (ObjectIdPtr) featSip->data.ptrvalue;
    
    /* If the product Seq ID matches the old Seq Id */
    /* then change the product to point to the new. */
    
    if (NULL == oip->str) {
      sprintf (tmpIdStr, "%d", oip->id);
      if (StringCmp (tmpIdStr, idsPtr->oldStr) == 0)
	SeqLocReplaceLocalID (sfp->product, idsPtr->newSip);
    }
    if (StringCmp (oip->str, idsPtr->oldStr) == 0)
      SeqLocReplaceLocalID (sfp->product, idsPtr->newSip);
    
  }
}


static void ReplaceLocalID (BioseqPtr bsp,
			    SeqIdPtr sip,
			    CharPtr key,
			    Int2 count)

{
  ObjectIdPtr      oip;
  Char             str [64];
  Char             tmp [70];
  BioseqSetPtr     bssp = NULL;
  ReplaceIDStruct  ids;
  BioseqPtr        siblingBsp;
  SeqEntryPtr      sep;
  Int2             parentType;

  if (bsp == NULL || sip == NULL || StringHasNoText (key)) return;
  oip = (ObjectIdPtr) sip->data.ptrvalue;
  if (oip == NULL) return;

  /* Create the new ID string */

  StringNCpy_0 (str, key, sizeof (str));
  sprintf (tmp, "%s__%d", str, (int) count);

  /* Save the original SeqId for later passing */
  /* to CheckSetForNuclID_callback () and      */
  /* CheckSetForProductId_callback ().         */

  if (NULL != oip->str)
    ids.oldStr = StringSave (oip->str);
  else {
    ids.oldStr = (CharPtr) MemNew (32);
    sprintf (ids.oldStr, "%d", oip->id);
  }
    

  /* Update the Seq ID with the new string */

  oip->str = StringSave (tmp);
  ids.newSip = sip;
  SeqMgrReplaceInBioseqIndex (bsp);

  /* Replace the local ID on all the features of the bioseq */

  VisitFeaturesOnBsp (bsp, (Pointer) sip, ReplaceLocalIdOnLoc_callback);

  /* Check the parent (and grandparent, etc.) BioseqSet */
  /* for features that use the changed ID.              */

  parentType = bsp->idx.parenttype;
  if (parentType == OBJ_BIOSEQSET) 
    bssp = (BioseqSetPtr) bsp->idx.parentptr;

  while (parentType == OBJ_BIOSEQSET) {

    if ((bssp != NULL) && (bssp->_class == 1)) {
      
      /* Check features that are attached to */
      /* the parent set itself.              */
      
      if (ISA_na(bsp->mol))
	VisitFeaturesOnSet (bssp, (Pointer) &ids,
			    CheckFeatForNuclID_callback);
      else if (ISA_aa(bsp->mol))
	VisitFeaturesOnSet (bssp, (Pointer) &ids,
			    CheckFeatForProductID_callback);
      
      /* Check features that are attached to */
      /* other Bioseqs in the set.           */
      
      sep = bssp->seqentry;
      while (NULL != sep) {
	if (sep->choice == 1) { /* bioseq */
	  siblingBsp = (BioseqPtr) sep->data.ptrvalue;
	  if (ISA_na(bsp->mol))
	    VisitFeaturesOnBsp (siblingBsp, (Pointer) sip,
				CheckFeatForNuclID_callback);
	  else if (ISA_aa(bsp->mol))
	    VisitFeaturesOnBsp (siblingBsp, (Pointer) sip,
				CheckFeatForProductID_callback);
	}
	sep = sep->next;
      }
      
      sep = bssp->seq_set;
      while (NULL != sep) {
	if (sep->choice == 1) { /* bioseq */
	  siblingBsp = (BioseqPtr) sep->data.ptrvalue;
	  if (ISA_na(bsp->mol))
	    VisitFeaturesOnBsp (siblingBsp, (Pointer) sip,
				CheckFeatForNuclID_callback);
	  else if (ISA_aa(bsp->mol))
	    VisitFeaturesOnBsp (siblingBsp, (Pointer) sip,
				CheckFeatForProductID_callback);
	}
	sep = sep->next;
      }
    }
    parentType = bssp->idx.parenttype;
    bssp = (BioseqSetPtr) bssp->idx.parentptr;
  }

  /* Clean up before exiting */

  MemFree (ids.oldStr);

}


static void BuildLclTree (LclIdListPtr PNTR head, BioseqPtr bsp, CharPtr x, SeqIdPtr sip)

{
  Int2          comp;
  LclIdListPtr  idlist;

  if (*head != NULL) {
    idlist = *head;
    comp = StringICmp (idlist->key, x);
    if (comp < 0) {
      BuildLclTree (&(idlist->right), bsp, x, sip);
    } else if (comp > 0) {
      BuildLclTree (&(idlist->left), bsp, x, sip);
    } else {
      if (idlist->firstbsp != NULL && idlist->firstsip != NULL) {
        ReplaceLocalID (idlist->firstbsp, idlist->firstsip, x, 1);
        idlist->count = 2;
        idlist->firstbsp = NULL;
        idlist->firstsip = NULL;
      }
      ReplaceLocalID (bsp, sip, x, idlist->count);
      (idlist->count)++;
    }
  } else {
    idlist = MemNew (sizeof (LclIdList));
    if (idlist != NULL) {
      *head = idlist;
      idlist->firstbsp = bsp;
      idlist->firstsip = sip;
      idlist->count = 1;
      idlist->key = StringSave (x);
      idlist->left = NULL;
      idlist->right = NULL;
    }
  }
}

NLM_EXTERN void FreeLclTree (LclIdListPtr PNTR head)

{
  LclIdListPtr  idlist;

  if (head != NULL && *head != NULL) {
    idlist = *head;
    FreeLclTree (&(idlist->left));
    FreeLclTree (&(idlist->right));
    MemFree (idlist->key);
    MemFree (idlist);
  }
}


NLM_EXTERN void ResolveExistingIDsCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr          bsp;
  LclIdListPtr PNTR  head;
  SeqIdPtr           sip;
  Char               str [64];

  head = (LclIdListPtr PNTR) mydata;
  if (sep == NULL || head == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp != NULL) {
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_LOCAL) {
          SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
          BuildLclTree (head, bsp, str, sip);
        }
      }
    }
  }
}


static Boolean DoesIdListHaveLocal (SeqIdPtr sip)
{
  while (sip != NULL) {
    if (sip->choice == SEQID_LOCAL) {
      return TRUE;
    }
    sip = sip->next;
  }
  return FALSE;
}


static Boolean DoesSeqLocListHaveLocalId (SeqLocPtr slp)
{
  SeqLocPtr      loc;
  PackSeqPntPtr  psp;
  SeqBondPtr     sbp;
  SeqIntPtr      sinp;
  SeqIdPtr       sip;
  SeqPntPtr      spp;
  Boolean        has_local = FALSE;

  while (slp != NULL) {
    switch (slp->choice) {
      case SEQLOC_NULL :
        break;
      case SEQLOC_EMPTY :
      case SEQLOC_WHOLE :
        sip = (SeqIdPtr) slp->data.ptrvalue;
        has_local = DoesIdListHaveLocal (sip);
        break;
      case SEQLOC_INT :
        sinp = (SeqIntPtr) slp->data.ptrvalue;
        if (sinp != NULL) {
          sip = sinp->id;
          has_local = DoesIdListHaveLocal (sip);
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL) {
          sip = spp->id;
          has_local = DoesIdListHaveLocal (sip);
        }
        break;
      case SEQLOC_PACKED_PNT :
        psp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (psp != NULL) {
          sip = psp->id;
          has_local = DoesIdListHaveLocal (sip);
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        loc = (SeqLocPtr) slp->data.ptrvalue;
        while (loc != NULL && !has_local) {
          has_local = DoesSeqLocListHaveLocalId(loc);
          loc = loc->next;
        }
        break;
      case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          spp = (SeqPntPtr) sbp->a;
          if (spp != NULL) {
            sip = spp->id;
            has_local = DoesIdListHaveLocal (sip);
          }
          spp = (SeqPntPtr) sbp->b;
          if (spp != NULL) {
            sip = spp->id;
            has_local = DoesIdListHaveLocal (sip);
          }
        }
        break;
      case SEQLOC_FEAT :
        break;
      default :
        break;
    }
    slp = slp->next;
  }
  return FALSE;
}


static void SeqEntryHasAlignmentsWithLocalIDsCallback (SeqAnnotPtr sap, Pointer userdata)
{
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  PackSegPtr    psp;
  SeqAlignPtr   salp;
  StdSegPtr     ssp;
  Boolean       has_local = FALSE;
  BoolPtr     bp;

  if (sap == NULL || sap->type != 2 || userdata == NULL) return;
  salp = (SeqAlignPtr) sap->data;
  if (salp != NULL) 
  {
    switch (salp->segtype) {
      case SAS_DENDIAG :
        for (ddp = salp->segs; ddp != NULL && !has_local; ddp = ddp->next) {
          has_local = DoesIdListHaveLocal (ddp->id);
        }
        break;
      case SAS_DENSEG :
        dsp = salp->segs;
        if (dsp != NULL) {
          has_local = DoesIdListHaveLocal (dsp->ids);
        }
        break;
      case SAS_STD :
        for (ssp = salp->segs; ssp != NULL && !has_local; ssp = ssp->next) {
          has_local = DoesIdListHaveLocal (ssp->ids);
          if (!has_local) {
            has_local = DoesSeqLocListHaveLocalId (ssp->loc);
          }
        }
        break;
      case SAS_PACKED :
        psp = (PackSegPtr) salp->segs;
        if (psp != NULL) {
          has_local = DoesIdListHaveLocal (psp->ids);
        }
        break;
      default :
        break;
    }
  }

  bp = (BoolPtr) userdata;
  *bp |= has_local;
}


NLM_EXTERN Boolean HasAlignmentsWithLocalIDs (SeqEntryPtr sep)
{
  Boolean has_alignments = FALSE;

  VisitAnnotsInSep (sep, (Pointer) &has_alignments, SeqEntryHasAlignmentsWithLocalIDsCallback);

  return has_alignments;
}

/*
  EC number replacement - copied from Sequin, with protection
  multiple reads if no replacement file available
*/

typedef struct ecrepdata {
  CharPtr  before;
  CharPtr  after;
} EcRepData, PNTR EcRepPtr;

static ValNodePtr     ec_rep_list = NULL;
static EcRepPtr PNTR  ec_rep_data = NULL;
static Int4           ec_rep_len = 0;
static Boolean        ec_rep_read = FALSE;

static int LIBCALLBACK SortVnpByEcBefore (VoidPtr ptr1, VoidPtr ptr2)

{
  EcRepPtr    erp1, erp2;
  CharPtr     str1, str2;
  ValNodePtr  vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  erp1 = (EcRepPtr) vnp1->data.ptrvalue;
  erp2 = (EcRepPtr) vnp2->data.ptrvalue;
  if (erp1 == NULL || erp2 == NULL) return 0;
  str1 = erp1->before;
  str2 = erp2->before;
  if (str1 == NULL || str2 == NULL) return 0;
  return StringCmp (str1, str2);
}

static void SetupECReplacementTable (CharPtr file)

{
  EcRepPtr    erp;
  FileCache   fc;
  FILE        *fp = NULL;
  Int4        i;
  ValNodePtr  last = NULL;
  Char        line [512];
  Char        path [PATH_MAX];
  CharPtr     ptr;
  ErrSev      sev;
  CharPtr     str;
  ValNodePtr  vnp;

  if (ec_rep_data != NULL) return;
  if (ec_rep_read) return;

  if (FindPath ("ncbi", "ncbi", "data", path, sizeof (path))) {
    FileBuildPath (path, NULL, file);
    sev = ErrSetMessageLevel (SEV_ERROR);
    fp = FileOpen (path, "r");
    ErrSetMessageLevel (sev);
    if (fp != NULL) {
      FileCacheSetup (&fc, fp);
  
      str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      while (str != NULL) {
        if (StringDoesHaveText (str)) {
          ptr = StringChr (str, '\t');
          if (ptr != NULL) {
            *ptr = '\0';
            ptr++;
            /* only replace if a single destination number, not a split from one to many */
            if (StringChr (ptr, '\t') == NULL) {
              erp = (EcRepPtr) MemNew (sizeof (EcRepData));
              if (erp != NULL) {
                erp->before = StringSave (str);
                erp->after = StringSave (ptr);
                vnp = ValNodeAddPointer (&last, 0, (Pointer) erp);
                if (ec_rep_list == NULL) {
                  ec_rep_list = vnp;
                }
                last = vnp;
              }
            }
          }
        }
        str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
      }

      FileClose (fp);
      ec_rep_len = ValNodeLen (ec_rep_list);
      if (ec_rep_len > 0) {
        ec_rep_list = ValNodeSort (ec_rep_list, SortVnpByEcBefore);
        ec_rep_data = (EcRepPtr PNTR) MemNew (sizeof (EcRepPtr) * (ec_rep_len + 1));
        if (ec_rep_data != NULL) {
          for (vnp = ec_rep_list, i = 0; vnp != NULL; vnp = vnp->next, i++) {
            erp = (EcRepPtr) vnp->data.ptrvalue;
            ec_rep_data [i] = erp;
          }
        }
      }
    }
  }

  ec_rep_read = TRUE;
}

static CharPtr GetECReplacement (CharPtr str)

{
  EcRepPtr  erp;
  Int4      L, R, mid;
 
  if (StringHasNoText (str)) return NULL;

  L = 0;
  R = ec_rep_len - 1;
  while (L < R) {
    mid = (L + R) / 2;
    erp = ec_rep_data [(int) mid];
    if (erp != NULL && StringCmp (erp->before, str) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }
  erp = ec_rep_data [(int) R];
  if (erp != NULL && StringCmp (erp->before, str) == 0 && StringChr (erp->after, '\t') == NULL) return erp->after;

  return NULL;
}

static void UpdateProtEC (SeqFeatPtr sfp, Pointer userdata)

{
  Int4Ptr     countP;
  Int2        inf_loop_check;
  CharPtr     nxt;
  ProtRefPtr  prp;
  CharPtr     rep;
  CharPtr     str;
  ValNodePtr  vnp;
 
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT) return;
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL || prp->ec == NULL) return;
  countP = (Int4Ptr) userdata;

  for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    rep = GetECReplacement (str);
    if (rep == NULL) continue;
    nxt = rep;
    inf_loop_check = 0;
    while (nxt != NULL && inf_loop_check < 10) {
      rep = nxt;
      inf_loop_check++;
      nxt = GetECReplacement (rep);
    }
    vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
    vnp->data.ptrvalue = StringSave (rep);
    if (countP != NULL) {
      (*countP)++;
    }
  }
}

NLM_EXTERN Int4 UpdateReplacedECNumbers (SeqEntryPtr sep)

{
  Int4  count = 0;

  if (sep == NULL) return 0;

  SetupECReplacementTable ("ecnum_replaced.txt");
  if (ec_rep_data != NULL && ec_rep_len > 0) {
    VisitFeaturesInSep (sep, (Pointer) &count, UpdateProtEC);
  }

  return count;
}

static void DeleteBadProtEC (SeqFeatPtr sfp, Pointer userdata)

{
  Int4Ptr     countP;
  ProtRefPtr  prp;
  CharPtr     str;
  ValNodePtr  vnp;
 
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT) return;
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL || prp->ec == NULL) return;
  countP = (Int4Ptr) userdata;

  for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    if (ValidateECnumber (str) && (! ECnumberNotInList (str))) continue;
    vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
    if (countP != NULL) {
      (*countP)++;
    }
  }
}

NLM_EXTERN Int4 DeleteBadECNumbers (SeqEntryPtr sep)

{
  Int4  count = 0;

  if (sep == NULL) return 0;

  VisitFeaturesInSep (sep, (Pointer) &count, DeleteBadProtEC);

  return count;
}




