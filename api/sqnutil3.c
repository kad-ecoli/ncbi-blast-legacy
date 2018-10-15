/*   sqnutil3.c
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
* File Name:  sqnutil3.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   2/7/00
*
* $Revision: 6.805 $
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
#include <gbfeat.h>
#include <gbftdef.h>
#include <edutil.h>
#include <tofasta.h>
#include <parsegb.h>
#include <utilpars.h>
#include <validatr.h>
#include <explore.h>
#include <salsap.h>
#include <salutil.h>
#include <salpedit.h>
#include <alignmgr2.h>
#include <actutils.h>
#include <utilpub.h>
/* included for discrepancy report */
#include <asn2gnbk.h>
#include <asn2gnbp.h>
#include <valid.h>
#include <findrepl.h>

#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>

/* functions for associating CDS and parent mRNA using featureIDs */

NLM_EXTERN void ClearFeatIDs (
  SeqFeatPtr sfp
)

{
  if (sfp == NULL) return;
  SeqFeatIdFree (&sfp->id);
  sfp->id.choice = 0;
}

NLM_EXTERN void ClearFeatIDXrefs (
  SeqFeatPtr sfp
)

{
  SeqFeatXrefPtr  xref, next, PNTR prevlink;

  if (sfp == NULL) return;

  prevlink = (SeqFeatXrefPtr PNTR) &(sfp->xref);
  xref = sfp->xref;
  while (xref != NULL) {
    next = xref->next;

    if (xref->id.choice != 0) {
      SeqFeatIdFree (&xref->id);
      xref->id.choice = 0;
    }
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

static void SfpClearFeatIDs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  if (sfp == NULL) return;
  ClearFeatIDs (sfp);
  ClearFeatIDXrefs (sfp);
}

NLM_EXTERN void ClearFeatureIDs (
  SeqEntryPtr sep
)

{
  VisitFeaturesInSep (sep, NULL, SfpClearFeatIDs);
}

typedef struct idpair {
  Int4  before;
  Int4  after;
} IdPairData, PNTR IdPairPtr;

typedef struct fiddata {
  Int4       highestID;
  Int4       highestRef;
  Int4       offset;
  Int4       count;
  IdPairPtr  pairs;
} FidData, PNTR FidDataPtr;

static void FindHighestFeatID (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  FidDataPtr      fip;
  ObjectIdPtr     oip;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL) return;
  fip = (FidDataPtr) userdata;
  if (fip == NULL) return;

  if (sfp->id.choice == 3) {
    oip = (ObjectIdPtr) sfp->id.value.ptrvalue;
    if (oip != NULL) {
      if (oip->str == NULL) {
        if (oip->id >= fip->highestID) {
          fip->highestID = oip->id;
        }
      }
    }
  }

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice != 3) continue;
    oip = (ObjectIdPtr) xref->id.value.ptrvalue;
    if (oip != NULL) {
      if (oip->str == NULL) {
        if (oip->id >= fip->highestRef) {
          fip->highestRef = oip->id;
        }
      }
    }
  }
}

NLM_EXTERN Int4 FindHighestFeatureID (
  SeqEntryPtr sep
)

{
  FidData  fd;

  MemSet ((Pointer) &fd, 0, sizeof (FidData));
  fd.highestID = 0;
  fd.highestRef = 0;
  VisitFeaturesInSep (sep, (Pointer) &fd, FindHighestFeatID);
  return fd.highestID;
}

static void SfpAssignFeatIDs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  FidDataPtr   fip;
  ObjectIdPtr  oip;

  if (sfp == NULL) return;
  fip = (FidDataPtr) userdata;
  if (fip == NULL) return;

  if (sfp->id.choice == 3) return;
  oip = ObjectIdNew ();
  if (oip == NULL) return;

  (fip->highestID)++;
  oip->id = fip->highestID;

  sfp->id.value.ptrvalue = (Pointer) oip;
  sfp->id.choice = 3;
}

NLM_EXTERN void AssignFeatureIDs (
  SeqEntryPtr sep
)

{
  FidData  fd;

  MemSet ((Pointer) &fd, 0, sizeof (FidData));
  fd.highestID = 0;
  fd.highestRef = 0;
  VisitFeaturesInSep (sep, (Pointer) &fd, FindHighestFeatID);
  VisitFeaturesInSep (sep, (Pointer) &fd, SfpAssignFeatIDs);
}

static void SfpOffsetFeatIDs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  FidDataPtr   fip;
  ObjectIdPtr  oip;

  if (sfp == NULL) return;
  fip = (FidDataPtr) userdata;
  if (fip == NULL) return;

  if (sfp->id.choice == 3) {
    oip = (ObjectIdPtr) sfp->id.value.ptrvalue;
    if (oip != NULL) {
      if (oip->str == NULL) {
        oip->id += fip->offset;
      }
    }
  }
}

NLM_EXTERN void OffsetFeatureIDs (
  SeqEntryPtr sep,
  Int4 offset
)

{
  FidData  fd;

  MemSet ((Pointer) &fd, 0, sizeof (FidData));
  fd.offset = offset;
  VisitFeaturesInSep (sep, (Pointer) &fd, SfpOffsetFeatIDs);
}

static void SfpOffsetFeatIDXrefs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  FidDataPtr      fip;
  ObjectIdPtr     oip;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL) return;
  fip = (FidDataPtr) userdata;
  if (fip == NULL) return;

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice != 3) continue;
    oip = (ObjectIdPtr) xref->id.value.ptrvalue;
    if (oip != NULL) {
      if (oip->str == NULL) {
        oip->id += fip->offset;
      }
    }
  }
}

NLM_EXTERN void OffsetFeatureIDXrefs (
  SeqEntryPtr sep,
  Int4 offset
)

{
  FidData  fd;

  MemSet ((Pointer) &fd, 0, sizeof (FidData));
  fd.offset = offset;
  VisitFeaturesInSep (sep, (Pointer) &fd, SfpOffsetFeatIDXrefs);
}

static void SfpMakePairList (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  FidDataPtr   fip;
  Int4         idx;
  IdPairPtr    ipp;
  ObjectIdPtr  oip;

  if (sfp == NULL) return;
  fip = (FidDataPtr) userdata;
  if (fip == NULL) return;
  if (fip->pairs == NULL) return;

  if (sfp->id.choice != 3) return;
  oip = (ObjectIdPtr) sfp->id.value.ptrvalue;
  if (oip == NULL) return;

  idx = fip->highestID;
  ipp = &(fip->pairs [idx]);

  (fip->highestID)++;
  ipp->before = oip->id;
  ipp->after = fip->highestID;
}

static int LIBCALLBACK SortPairList (VoidPtr ptr1, VoidPtr ptr2)

{
  IdPairPtr  ipp1 = (IdPairPtr) ptr1;
  IdPairPtr  ipp2 = (IdPairPtr) ptr2;

  if (ipp1 == NULL || ipp2 == NULL) return 0;
  if (ipp1->before > ipp2->before) return 1;
  if (ipp1->before < ipp2->before) return -1;
  return 0;
}

static Int4 LookupNewFeatID (
  FidDataPtr fip,
  Int4 before
)

{
  IdPairPtr  ipp;
  Int4       L;
  Int4       mid;
  Int4       R;

  if (fip == NULL || fip->pairs == NULL || fip->count < 1) return 0;

  L = 0;
  R = fip->count - 1;
  while (L < R) {
    mid = (L + R) / 2;
    ipp = &(fip->pairs [mid]);
    if (ipp->before < before) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (R < fip->count) {
    ipp = &(fip->pairs [R]);
    if (ipp->before == before) return ipp->after;
  }

  return 0;
}

static void SfpReassignPairList (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  FidDataPtr      fip;
  ObjectIdPtr     oip;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL) return;
  fip = (FidDataPtr) userdata;
  if (fip == NULL) return;
  if (fip->pairs == NULL) return;

  if (sfp->id.choice == 3) {
    oip = (ObjectIdPtr) sfp->id.value.ptrvalue;
    if (oip != NULL) {
      if (oip->str == NULL) {
        oip->id = LookupNewFeatID (fip, oip->id);
      }
    }
  }

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice != 3) continue;
    oip = (ObjectIdPtr) xref->id.value.ptrvalue;
    if (oip != NULL) {
      if (oip->str == NULL) {
        oip->id = LookupNewFeatID (fip, oip->id);
      }
    }
  }
}

NLM_EXTERN void ReassignFeatureIDs (
  SeqEntryPtr sep
)

{
  Int4     count;
  FidData  fd;

  count = VisitFeaturesInSep (sep, NULL, NULL);
  if (count < 1) return;

  MemSet ((Pointer) &fd, 0, sizeof (FidData));
  fd.highestID = 0;
  fd.highestRef = 0;
  fd.count = count;
  fd.pairs = (IdPairPtr) MemNew (sizeof (IdPairData) * (count + 1));
  if (fd.pairs == NULL) return;

  VisitFeaturesInSep (sep, (Pointer) &fd, SfpMakePairList);

  StableMergeSort (fd.pairs, (size_t) count, sizeof (IdPairData), SortPairList);

  VisitFeaturesInSep (sep, (Pointer) &fd, SfpReassignPairList);

  MemFree (fd.pairs);
}

typedef struct vcmdata {
  Boolean     accounted_for;
  SeqFeatPtr  cds;
  SeqFeatPtr  mrna;
  SeqFeatPtr  partner;
} VcmData, PNTR VcmDataPtr;

typedef struct loopdata {
  Int2        count;
  SeqFeatPtr  cds;
  SeqFeatPtr  mrna;
} LoopData, PNTR LoopDataPtr;

static Boolean LIBCALLBACK GetSingleMrnaProc (
  SeqFeatPtr mrna,
  SeqMgrFeatContextPtr context
)

{
  LoopDataPtr  ldp;
  VcmDataPtr   vdp;

  ldp = (LoopDataPtr) context->userdata;

  vdp = (VcmDataPtr) mrna->idx.scratch;
  if (vdp != NULL && vdp->accounted_for) return TRUE;

  (ldp->count)++;
  ldp->mrna = mrna;

  return TRUE;
}

static void BspLinkCDSmRNAbyOverlap (
  BioseqPtr bsp,
  Pointer userdata
)

{
  Int2               count;
  SeqMgrFeatContext  fcontext;
  Boolean            goOn;
  Int4               id;
  LoopData           ld;
  ObjectIdPtr        oip;
  SeqFeatPtr         partner, sfp;
  VcmDataPtr         vdp;
  SeqFeatXrefPtr     xref;

  if (bsp == NULL || ISA_aa (bsp->mol)) return;

  /* add scratch structure to CDS and mRNA features */

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
  while (sfp != NULL) {
    sfp->idx.scratch = (Pointer) MemNew (sizeof (VcmData));
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_mRNA, &fcontext);
  while (sfp != NULL) {
    sfp->idx.scratch = (Pointer) MemNew (sizeof (VcmData));
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_mRNA, &fcontext);
  }

  /* loop through CDS features, finding single unused mRNA partner */

  goOn = TRUE;
  while (goOn) {
    goOn = FALSE;
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
    while (sfp != NULL) {
      vdp = (VcmDataPtr) sfp->idx.scratch;
      if (vdp != NULL && (! vdp->accounted_for)) {
        ld.count = 0;
        ld.cds = sfp;
        ld.mrna = NULL;
        if (sfp->excpt &&
            (StringISearch (sfp->except_text, "ribosomal slippage") != NULL ||
             StringISearch (sfp->except_text, "trans-splicing") != NULL)) {
          count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_mRNA, NULL, 0,
                                                   LOCATION_SUBSET, (Pointer) &ld,
                                                   GetSingleMrnaProc);
        } else {
          count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_mRNA, NULL, 0,
                                                   CHECK_INTERVALS, (Pointer) &ld,
                                                   GetSingleMrnaProc);
        }
        if (ld.count == 1 && ld.mrna != NULL) {
          vdp->accounted_for = TRUE;
          vdp->cds = ld.cds;
          vdp->mrna = ld.mrna;
          vdp->partner = ld.mrna;
          vdp = (VcmDataPtr) ld.mrna->idx.scratch;
          if (vdp != NULL) {
            vdp->accounted_for = TRUE;
            vdp->cds = ld.cds;
            vdp->mrna = ld.mrna;
            vdp->partner = ld.cds;
            goOn = TRUE;
          }
        }
      }
      sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
    }
  }

  /* assign xrefs between CDS and mRNA features */

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    vdp = (VcmDataPtr) sfp->idx.scratch;
    if (vdp != NULL && vdp->accounted_for) {
      partner = vdp->partner;
      if (partner != NULL && partner->id.choice == 3) {
        oip = (ObjectIdPtr) partner->id.value.ptrvalue;
        if (oip != NULL && oip->str == NULL) {
          id = oip->id;
          if (id > 0) {
            for (xref = sfp->xref; xref != NULL && xref->id.choice != 3; xref = xref->next) continue;
            if (xref != NULL) {
              oip = (ObjectIdPtr) xref->id.value.ptrvalue;
              if (oip != NULL) {
                if (oip->str != NULL) {
                  oip->str = MemFree (oip->str);
                }
                oip->id = id;
              }
            } else {
              xref = SeqFeatXrefNew ();
              if (xref != NULL) {
                oip = ObjectIdNew ();
                if (oip != NULL) {
                  oip->id = id;
                  xref->id.choice = 3;
                  xref->id.value.ptrvalue = (Pointer) oip;
                  xref->next = sfp->xref;
                  sfp->xref = xref;
                }
              }
            }
          }
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  /* free scratch structure in CDS and mRNA features */

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    if (sfp->idx.scratch != NULL) {
      sfp->idx.scratch = MemFree (sfp->idx.scratch);
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }
}

NLM_EXTERN void LinkCDSmRNAbyOverlap (
  SeqEntryPtr sep
)

{
  AssignFeatureIDs (sep);
  VisitBioseqsInSep (sep, NULL, BspLinkCDSmRNAbyOverlap);
}

static void BspLinkCDSmRNAbyLabel (
  BioseqPtr bsp,
  Pointer userdata
)

{
  SeqFeatPtr         cds, mrna;
  SeqMgrFeatContext  ccontext;
  SeqMgrFeatContext  mcontext;
  Int4               id;
  ObjectIdPtr        oip;
  SeqFeatXrefPtr     xref;

  if (bsp == NULL || ISA_aa (bsp->mol)) return;

  /* loop through CDS features, finding mRNA partner by label */

  cds = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &ccontext);
  while (cds != NULL) {
    if (StringDoesHaveText (ccontext.label)) {
      mrna = SeqMgrGetFeatureByLabel (bsp, ccontext.label, 0, FEATDEF_mRNA, &mcontext);
      if (mrna != NULL && StringCmp (ccontext.label, mcontext.label) == 0) {
        if (cds->id.choice == 3 && mrna->id.choice == 3) {

          /* assign xrefs between CDS and mRNA features */

          oip = (ObjectIdPtr) mrna->id.value.ptrvalue;
          if (oip != NULL && oip->str == NULL) {
            id = oip->id;
            if (id > 0) {
              for (xref = cds->xref; xref != NULL && xref->id.choice != 3; xref = xref->next) continue;
              if (xref != NULL) {
                oip = (ObjectIdPtr) xref->id.value.ptrvalue;
                if (oip != NULL) {
                  if (oip->str != NULL) {
                    oip->str = MemFree (oip->str);
                  }
                  oip->id = id;
                }
              } else {
                xref = SeqFeatXrefNew ();
                if (xref != NULL) {
                  oip = ObjectIdNew ();
                  if (oip != NULL) {
                    oip->id = id;
                    xref->id.choice = 3;
                    xref->id.value.ptrvalue = (Pointer) oip;
                    xref->next = cds->xref;
                    cds->xref = xref;
                  }
                }
              }
            }
          }

          oip = (ObjectIdPtr) cds->id.value.ptrvalue;
          if (oip != NULL && oip->str == NULL) {
            id = oip->id;
            if (id > 0) {
              for (xref = mrna->xref; xref != NULL && xref->id.choice != 3; xref = xref->next) continue;
              if (xref != NULL) {
                oip = (ObjectIdPtr) xref->id.value.ptrvalue;
                if (oip != NULL) {
                  if (oip->str != NULL) {
                    oip->str = MemFree (oip->str);
                  }
                  oip->id = id;
                }
              } else {
                xref = SeqFeatXrefNew ();
                if (xref != NULL) {
                  oip = ObjectIdNew ();
                  if (oip != NULL) {
                    oip->id = id;
                    xref->id.choice = 3;
                    xref->id.value.ptrvalue = (Pointer) oip;
                    xref->next = mrna->xref;
                    mrna->xref = xref;
                  }
                }
              }
            }
          }
        }
      }
    }
    cds = SeqMgrGetNextFeature (bsp, cds, SEQFEAT_CDREGION, 0, &ccontext);
  }
}

NLM_EXTERN void LinkCDSmRNAbyLabel (
  SeqEntryPtr sep
)

{
  AssignFeatureIDs (sep);
  VisitBioseqsInSep (sep, NULL, BspLinkCDSmRNAbyLabel);
}


static void MakeOneLink (
  SeqFeatPtr f1, 
  SeqFeatPtr f2
)

{
  ObjectIdPtr        oip;
  SeqFeatXrefPtr     xref;
  Int4               id;

  if (f1 == NULL || f2 == NULL || f1->id.choice != 3 || f2->id.choice != 3) {
    return;
  }

  oip = (ObjectIdPtr) f1->id.value.ptrvalue;
  if (oip != NULL && oip->str == NULL) {
    id = oip->id;
    if (id > 0) {
      for (xref = f2->xref; xref != NULL && xref->id.choice != 3; xref = xref->next) continue;
      if (xref != NULL) {
        oip = (ObjectIdPtr) xref->id.value.ptrvalue;
        if (oip != NULL) {
          if (oip->str != NULL) {
            oip->str = MemFree (oip->str);
          }
          oip->id = id;
        }
      } else {
        xref = SeqFeatXrefNew ();
        if (xref != NULL) {
          oip = ObjectIdNew ();
          if (oip != NULL) {
            oip->id = id;
            xref->id.choice = 3;
            xref->id.value.ptrvalue = (Pointer) oip;
            xref->next = f2->xref;
            f2->xref = xref;
          }
        }
      }
    }
  }
}


static void CreateReciprocalLink (
  SeqFeatPtr f1, 
  SeqFeatPtr f2
)

{
  if (f1 == NULL || f2 == NULL || f1->id.choice != 3 || f2->id.choice != 3) {
    return;
  }

  MakeOneLink (f1, f2);
  MakeOneLink (f2, f1);
}


static void LinkCDSmRNAbyLabelAndLocationCallback (
  BioseqPtr bsp, 
  Pointer userdata
)

{
  SMFeatItemPtr PNTR  array;
  BioseqExtraPtr      bspextra;
  Uint2               entityID;
  SMFeatItemPtr       feat;
  Int4                i, j, best_index, best_diff, diff;
  Int4                num;
  ObjMgrDataPtr       omdp;

  if (bsp == NULL) return;

  omdp = SeqMgrGetOmdpForBioseq (bsp);
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return;

  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return;
  array = bspextra->featsByLabel;
  num = bspextra->numfeats;
  if (array == NULL || num < 1) return;

  entityID = bsp->idx.entityID;
  if (entityID < 1) {
    entityID = ObjMgrGetEntityIDForPointer (omdp->dataptr);
  }

  /* labels are all grouped together - for each cds/mRNA in group of identical labels,
   * find match with best location.
   */
  for (i = 0; i < num - 1; i++) {
    feat = array [i];
    if (feat->sfp == NULL) {
      continue;
    } else if (feat->sfp->xref != NULL) {
      /* already assigned feat xref */
      continue;
    } else if (feat->sfp->idx.subtype != FEATDEF_CDS && feat->sfp->idx.subtype != FEATDEF_mRNA) {
      /* not interested in these feature types */
    } else {
      best_index = -1;
      for (j = i + 1; j < num && StringCmp (feat->label, array[j]->label) == 0; j++) {
        if (array[j]->sfp == NULL) {
          /* bad */
        } else if (array[j]->sfp->xref != NULL) {
          /* already assigned feat xref */
        } else if (feat->sfp->idx.subtype == FEATDEF_CDS) {
          if (array[j]->sfp->idx.subtype != FEATDEF_mRNA) {
            /* wrong feature type */
          } else if ((diff = SeqLocAinB (feat->sfp->location, array[j]->sfp->location)) < 0) {
            /* locations don't match */
          } else {
            if (best_index == -1) {
              /* don't have a best yet */
              best_index = j;
              best_diff = diff;
            } else if (diff < best_diff) {
              best_index = j;
              best_diff = diff;
            }
          }
        } else if (feat->sfp->idx.subtype == FEATDEF_mRNA) {
          if (array[j]->sfp->idx.subtype != FEATDEF_CDS) {
            /* wrong feature type */
          } else if ((diff = SeqLocAinB (array[j]->sfp->location, feat->sfp->location)) < 0) {
            /* locations don't match */
          } else {
            if (best_index == -1) {
              /* don't have a best yet */
              best_index = j;
              best_diff = diff;
            } else if (diff < best_diff) {
              best_index = j;
              best_diff = diff;
            }
          }
        }
      }
      if (best_index > -1) {
        CreateReciprocalLink (feat->sfp, array[best_index]->sfp);
      }
    }
  }
}


NLM_EXTERN void LinkCDSmRNAbyLabelAndLocation (
  SeqEntryPtr sep
)

{
  AssignFeatureIDs (sep);
  VisitBioseqsInSep (sep, NULL, LinkCDSmRNAbyLabelAndLocationCallback);
}


typedef struct ovpdata {
  SeqFeatPtr  sfp;
  Char        revstr [42];
} OvpData, PNTR OvpDataPtr;

static int LIBCALLBACK SortOvpByString (VoidPtr ptr1, VoidPtr ptr2)

{
  OvpDataPtr  odp1;
  OvpDataPtr  odp2;
  CharPtr     str1;
  CharPtr     str2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  odp1 = *((OvpDataPtr PNTR) ptr1);
  odp2 = *((OvpDataPtr PNTR) ptr2);
  if (odp1 == NULL || odp2 == NULL) return 0;
  str1 = odp1->revstr;
  str2 = odp2->revstr;
  if (str1 == NULL || str2 == NULL) return 0;
  return StringICmp (str1, str2);
}

static void FindProtBsp (BioseqPtr bsp, Pointer userdata)

{
  BioseqPtr PNTR  protP;

  if (bsp == NULL || ! (ISA_aa (bsp->mol))) return;
  protP = (BioseqPtr PNTR) userdata;
  *protP = bsp;
}

static void BspLinkCDSmRNAbyProduct (
  BioseqPtr bsp,
  Pointer userdata
)

{
  BioseqSetPtr       bssp;
  Char               buf [42];
  BioseqPtr          cdna, prot;
  SeqFeatPtr         cds, mrna, sfp;
  OvpDataPtr         PNTR cdsarray = NULL, PNTR mrnaarray = NULL;
  ValNodePtr         cdshead = NULL, mrnahead = NULL, vnp;
  int                compare;
  Uint2              entityID;
  SeqMgrFeatContext  fcontext;
  Int2               i, numcds, nummrna, L, R, mid;
  Int4               id;
  OvpDataPtr         odp;
  ObjectIdPtr        oip;
  SeqEntryPtr        sep;
  SeqIdPtr           sip;
  SeqFeatXrefPtr     xref;

  if (bsp == NULL || ISA_aa (bsp->mol)) return;

  numcds = 0;
  nummrna = 0;

  /* count CDS and mRNA features, make revstr from product SeqId */

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    switch (sfp->idx.subtype) {
      case FEATDEF_CDS :
        if (sfp->product != NULL) {
          numcds++;
          sip = SeqLocId (sfp->product);
          if (sip == NULL) break;
          MakeReversedSeqIdString (sip, buf, sizeof (buf) - 1);
          if (StringHasNoText (buf)) break;
          odp = (OvpDataPtr) MemNew (sizeof (OvpData));
          if (odp == NULL) break;
          odp->sfp = sfp;
          StringCpy (odp->revstr, buf);
          vnp = ValNodeAddPointer (NULL, 0, (Pointer) odp);
          if (vnp == NULL) break;
          vnp->next = cdshead;
          cdshead = vnp;
        }
        break;
      case FEATDEF_mRNA :
        if (sfp->product != NULL) {
          nummrna++;
          sip = SeqLocId (sfp->product);
          if (sip == NULL) break;
          MakeReversedSeqIdString (sip, buf, sizeof (buf) - 1);
          if (StringHasNoText (buf)) break;
          odp = (OvpDataPtr) MemNew (sizeof (OvpData));
          if (odp == NULL) break;
          odp->sfp = sfp;
          StringCpy (odp->revstr, buf);
          vnp = ValNodeAddPointer (NULL, 0, (Pointer) odp);
          if (vnp == NULL) break;
          vnp->next = mrnahead;
          mrnahead = vnp;
        }
        break;
      default :
        break;
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  if (numcds > 0 && nummrna > 0) {
    cdsarray = (OvpDataPtr PNTR) MemNew (sizeof (OvpDataPtr) * (numcds + 1));
    mrnaarray = (OvpDataPtr PNTR) MemNew (sizeof (OvpDataPtr) * (nummrna + 1));

    /* populate and sort arrays to search for feature by product SeqId */

    if (cdsarray != NULL && mrnaarray != NULL) {
      for (vnp = cdshead, i = 0; vnp != NULL; vnp = vnp->next, i++) {
        cdsarray [i] = (OvpDataPtr) vnp->data.ptrvalue;
      }
      for (vnp = mrnahead, i = 0; vnp != NULL; vnp = vnp->next, i++) {
        mrnaarray [i] = (OvpDataPtr) vnp->data.ptrvalue;
      }

      StableMergeSort (cdsarray, (size_t) numcds, sizeof (OvpDataPtr), SortOvpByString);
      StableMergeSort (mrnaarray, (size_t) nummrna, sizeof (OvpDataPtr), SortOvpByString);

      for (i = 0; i < nummrna; i++) {
        odp = (OvpDataPtr) mrnaarray [i];
        if (odp == NULL) continue;
        mrna = odp->sfp;
        if (mrna == NULL || mrna->product == NULL) continue;
        sip = SeqLocId (mrna->product);
        if (sip == NULL) continue;

        cdna = BioseqLockById (sip);
        if (cdna == NULL) continue;
        entityID = ObjMgrGetEntityIDForPointer (cdna);
        if (entityID < 1) continue;
        if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
          sep = GetTopSeqEntryForEntityID (entityID);
          if (sep == NULL) continue;
          AssignIDsInEntity (entityID, 0, NULL);
        }
        if (cdna->idx.parenttype == OBJ_BIOSEQSET) {
          bssp = (BioseqSetPtr) cdna->idx.parentptr;
          if (bssp == NULL) continue;
          if (bssp->_class == BioseqseqSet_class_nuc_prot) {
            prot = NULL;
            if (VisitBioseqsInSet (bssp, (Pointer) &prot, FindProtBsp) == 2) {
              for (sip = prot->id; sip != NULL; sip = sip->next) {
                MakeReversedSeqIdString (sip, buf, sizeof (buf) - 1);
    
                /* binary search */
    
                L = 0;
                R = numcds - 1;
                while (L < R) {
                  mid = (L + R) / 2;
                  odp = cdsarray [mid];
                  compare = StringCmp (odp->revstr, buf);
                  if (compare < 0) {
                    L = mid + 1;
                  } else {
                    R = mid;
                  }
                }
                odp = cdsarray [R];
                if (odp != NULL && StringCmp (odp->revstr, buf) == 0) {
                  cds = odp->sfp;
                  if (cds == NULL) continue;
    
                  /* make reciprocal feature ID xrefs */
    
                  if (cds->id.choice == 3) {
                    oip = (ObjectIdPtr) cds->id.value.ptrvalue;
                    if (oip != NULL && oip->str == NULL) {
                      id = oip->id;
                      if (id > 0) {
                        for (xref = mrna->xref; xref != NULL && xref->id.choice != 3; xref = xref->next) continue;
                        if (xref != NULL) {
                          oip = (ObjectIdPtr) xref->id.value.ptrvalue;
                          if (oip != NULL) {
                            if (oip->str != NULL) {
                              oip->str = MemFree (oip->str);
                            }
                            oip->id = id;
                          }
                        } else {
                          xref = SeqFeatXrefNew ();
                          if (xref != NULL) {
                            oip = ObjectIdNew ();
                            if (oip != NULL) {
                              oip->id = id;
                              xref->id.choice = 3;
                              xref->id.value.ptrvalue = (Pointer) oip;
                              xref->next = mrna->xref;
                              mrna->xref = xref;
                            }
                          }
                        }
                      }
                    }
                  }
    
                  if (mrna->id.choice == 3) {
                    oip = (ObjectIdPtr) mrna->id.value.ptrvalue;
                    if (oip != NULL && oip->str == NULL) {
                      id = oip->id;
                      if (id > 0) {
                        for (xref = cds->xref; xref != NULL && xref->id.choice != 3; xref = xref->next) continue;
                        if (xref != NULL) {
                          oip = (ObjectIdPtr) xref->id.value.ptrvalue;
                          if (oip != NULL) {
                            if (oip->str != NULL) {
                              oip->str = MemFree (oip->str);
                            }
                            oip->id = id;
                          }
                        } else {
                          xref = SeqFeatXrefNew ();
                          if (xref != NULL) {
                            oip = ObjectIdNew ();
                            if (oip != NULL) {
                              oip->id = id;
                              xref->id.choice = 3;
                              xref->id.value.ptrvalue = (Pointer) oip;
                              xref->next = cds->xref;
                              cds->xref = xref;
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
        BioseqUnlock (cdna);
      }
    }

    /* clean up */

    MemFree (cdsarray);
    MemFree (mrnaarray);
  }

  /* more cleanup */

  ValNodeFreeData (cdshead);
  ValNodeFreeData (mrnahead);
}

NLM_EXTERN void LinkCDSmRNAbyProduct (
  SeqEntryPtr sep
)

{
  AssignFeatureIDs (sep);
  VisitBioseqsInSep (sep, NULL, BspLinkCDSmRNAbyProduct);
}

NLM_EXTERN void StripFeatIDXrefAsnFilter (
  AsnIoPtr aip,
  AsnIoPtr aop
)

{
  AsnModulePtr    amp;
  AsnTypePtr      atp, atp_se, atp_sfx, atp_sfxe;
  DataVal         dv;
  Boolean         inxrefs;
  SeqFeatXrefPtr  xref;

  if (aip == NULL || aop == NULL) return;

  amp = AsnAllModPtr ();
  if (amp == NULL) return;
  atp_se = AsnFind ("Seq-entry");
  atp_sfx = AsnFind ("Seq-feat.xref");
  atp_sfxe = AsnFind ("Seq-feat.xref.E");
  if (atp_se == NULL || atp_sfx == NULL || atp_sfxe == NULL) return;

  inxrefs = FALSE;
  atp = atp_se;

  while ((atp = AsnReadId (aip, amp, atp)) != NULL) {
    if (atp == atp_sfxe) {
      xref = SeqFeatXrefAsnRead (aip, atp);
      if (xref->data.choice != 0) {
        if (! inxrefs) {
          inxrefs = TRUE;
          AsnOpenStruct (aop, atp_sfx, (Pointer) NULL);
        }
        SeqFeatXrefAsnWrite (xref, aop, atp);
      }
      SeqFeatXrefFree (xref);
    } else if (atp == atp_sfx) {
      AsnReadVal (aip, atp, &dv);
      /* only send struct as open and close item */
      AsnKillValue (atp, &dv);
    } else {
      if (inxrefs) {
        AsnCloseStruct (aop, atp_sfx, (Pointer) NULL);
        inxrefs = FALSE;
      }
      AsnReadVal (aip, atp, &dv);
      AsnWrite (aop, atp, &dv);
      AsnKillValue (atp, &dv);
    }
  }
}

NLM_EXTERN void StripSeqDataGapAsnFilter (
  AsnIoPtr aip,
  AsnIoPtr aop
)

{
  AsnModulePtr  amp;
  AsnTypePtr    atp, atp_se, atp_dsl;
  DataVal       dv;
  SeqLitPtr     slp;

  if (aip == NULL || aop == NULL) return;

  amp = AsnAllModPtr ();
  if (amp == NULL) return;
  atp_se = AsnFind ("Seq-entry");
  atp_dsl = AsnFind ("Delta-seq.literal");
  if (atp_se == NULL || atp_dsl == NULL) return;

  atp = atp_se;

  while ((atp = AsnReadId (aip, amp, atp)) != NULL) {
    if (atp == atp_dsl) {
      slp = SeqLitAsnRead (aip, atp);
      if (slp != NULL && slp->seq_data != NULL && slp->seq_data_type == Seq_code_gap) {
        slp->seq_data = SeqDataFree (slp->seq_data, slp->seq_data_type);
      }
      SeqLitAsnWrite (slp, aop, atp);
      SeqLitFree (slp);
    } else {
      AsnReadVal (aip, atp, &dv);
      AsnWrite (aop, atp, &dv);
      AsnKillValue (atp, &dv);
    }
  }
}

NLM_EXTERN void StripNewFeatMolInfoFieldsAsnFilter (
  AsnIoPtr aip,
  AsnIoPtr aop
)

{
  AsnModulePtr   amp;
  AsnTypePtr     atp, atp_se, atp_sf, atp_mi;
  DataVal        dv;
  MolInfoPtr     mip;
  SeqFeatPtr     sfp;
  ValNodePtr     ids;
  UserObjectPtr  exts;
  CharPtr        gbmoltype;

  if (aip == NULL || aop == NULL) return;

  amp = AsnAllModPtr ();
  if (amp == NULL) return;
  atp_se = AsnFind ("Seq-entry");
  atp_sf = AsnFind ("Seq-feat");
  atp_mi = AsnFind ("Mol-info");
  if (atp_se == NULL || atp_sf == NULL || atp_mi == NULL) return;

  atp = atp_se;

  while ((atp = AsnReadId (aip, amp, atp)) != NULL) {
    if (atp == atp_sf) {
      sfp = SeqFeatAsnRead (aip, atp);
      ids = sfp->ids;
      exts = sfp->exts;
      sfp->ids = NULL;
      sfp->exts = NULL;
      SeqFeatAsnWrite (sfp, aop, atp);
      sfp->ids = ids;
      sfp->exts = exts;
      SeqFeatFree (sfp);
    } else if (atp == atp_mi) {
      mip = MolInfoAsnRead (aip, atp);
      gbmoltype = mip->gbmoltype;
      mip->gbmoltype = NULL;
      MolInfoAsnWrite (mip, aop, atp);
      mip->gbmoltype = gbmoltype;
      MolInfoFree (mip);
    } else {
      AsnReadVal (aip, atp, &dv);
      AsnWrite (aop, atp, &dv);
      AsnKillValue (atp, &dv);
    }
  }
}

NLM_EXTERN void StripPCRPrimerAsnFilter (
  AsnIoPtr aip,
  AsnIoPtr aop
)

{
  AsnModulePtr       amp;
  AsnTypePtr         atp, atp_se, atp_bs;
  BioSourcePtr       biop;
  DataVal            dv;
  PCRReactionSetPtr  pcr_primers;

  if (aip == NULL || aop == NULL) return;

  amp = AsnAllModPtr ();
  if (amp == NULL) return;
  atp_se = AsnFind ("Seq-entry");
  atp_bs = AsnFind ("BioSource");
  if (atp_se == NULL || atp_bs == NULL) return;

  atp = atp_se;

  while ((atp = AsnReadId (aip, amp, atp)) != NULL) {
    if (atp == atp_bs) {
      biop = BioSourceAsnRead (aip, atp);
      pcr_primers = biop->pcr_primers;
      biop->pcr_primers = NULL;
      BioSourceAsnWrite (biop, aop, atp);
      biop->pcr_primers = pcr_primers;
      BioSourceFree (biop);
    } else {
      AsnReadVal (aip, atp, &dv);
      AsnWrite (aop, atp, &dv);
      AsnKillValue (atp, &dv);
    }
  }
}

NLM_EXTERN void StripOrgNamePgcodeAsnFilter (
  AsnIoPtr aip,
  AsnIoPtr aop
)

{
  AsnModulePtr  amp;
  AsnTypePtr    atp, atp_se, atp_on;
  DataVal       dv;
  OrgNamePtr    onp;

  if (aip == NULL || aop == NULL) return;

  amp = AsnAllModPtr ();
  if (amp == NULL) return;
  atp_se = AsnFind ("Seq-entry");
  atp_on = AsnFind ("OrgName");
  if (atp_se == NULL || atp_on == NULL) return;

  atp = atp_se;

  while ((atp = AsnReadId (aip, amp, atp)) != NULL) {
    if (atp == atp_on) {
      onp = OrgNameAsnRead (aip, atp);
      onp->pgcode = 0;
      OrgNameAsnWrite (onp, aop, atp);
      OrgNameFree (onp);
    } else {
      AsnReadVal (aip, atp, &dv);
      AsnWrite (aop, atp, &dv);
      AsnKillValue (atp, &dv);
    }
  }
}

static void AddGBQualToFeature (
  SeqFeatPtr sfp,
  CharPtr qual,
  CharPtr val
)

{
  GBQualPtr  gbq;

  if (sfp == NULL || StringHasNoText (qual) || StringHasNoText (val)) return;

  gbq = GBQualNew ();
  if (gbq == NULL) return;

  gbq->qual = StringSave (qual);
  gbq->val = StringSave (val);

  gbq->next = sfp->qual;
  sfp->qual = gbq;
}

NLM_EXTERN void StripGeneRnaPcrAsnFilter (
  AsnIoPtr aip,
  AsnIoPtr aop
)

{
  AsnModulePtr         amp;
  AsnTypePtr           atp, atp_se, atp_sf, atp_bs;
  BioSourcePtr         biop;
  DataVal              dv;
  GeneNomenclaturePtr  formal_name = NULL;
  GeneRefPtr           grp = NULL;
  RNAGenPtr            rgp = NULL;
  RNAQualPtr           rqp;
  RnaRefPtr            rrp = NULL;
  SeqFeatPtr           sfp;
  PCRReactionSetPtr    pcr_primers = NULL;

  if (aip == NULL || aop == NULL) return;

  amp = AsnAllModPtr ();
  if (amp == NULL) return;
  atp_se = AsnFind ("Seq-entry");
  atp_sf = AsnFind ("Seq-annot.data.ftable.E");
  atp_bs = AsnFind ("Seqdesc.source");
  if (atp_se == NULL || atp_sf == NULL || atp_bs == NULL) return;

  atp = atp_se;

  while ((atp = AsnReadId (aip, amp, atp)) != NULL) {
    if (atp == atp_bs) {
      biop = BioSourceAsnRead (aip, atp);
      pcr_primers = biop->pcr_primers;
      biop->pcr_primers = NULL;
      BioSourceAsnWrite (biop, aop, atp);
      if (pcr_primers != NULL) {
        pcr_primers = PCRReactionSetFree (pcr_primers);
      }
      BioSourceFree (biop);
    } else if (atp == atp_sf) {
      sfp = SeqFeatAsnRead (aip, atp);
      if (sfp->data.choice == SEQFEAT_GENE) {
        grp = (GeneRefPtr) sfp->data.value.ptrvalue;
        if (grp != NULL) {
          formal_name = grp->formal_name;
          grp->formal_name = NULL;
        }
      } else if (sfp->data.choice == SEQFEAT_RNA) {
        rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
        if (rrp != NULL) {
          if (rrp->type == 8 || rrp->type == 9 || rrp->type == 10) {
            if (rrp->ext.choice == 3) {
              rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
              rrp->ext.value.ptrvalue = NULL;
              rrp->ext.choice = 0;
            }
            if (rrp->ext.choice == 0) {
              rrp->ext.choice = 1;
              switch (rrp->type) {
                case 8 :
                  rrp->ext.value.ptrvalue = StringSave ("ncRNA");
                  break;
                case 9 :
                  rrp->ext.value.ptrvalue = StringSave ("tmRNA");
                  break;
                case 10 :
                  rrp->ext.value.ptrvalue = StringSave ("misc_RNA");
                  break;
                default :
                  break;
              }
            }
            if (rgp != NULL) {
              if (StringDoesHaveText (rgp->_class)) {
                AddGBQualToFeature (sfp, "ncRNA_class", rgp->_class);
              }
              if (StringDoesHaveText (rgp->product)) {
                AddGBQualToFeature (sfp, "product", rgp->product);
              }
              for (rqp = rgp->quals; rqp != NULL; rqp = rqp->next) {
                AddGBQualToFeature (sfp, rqp->qual, rqp->val);
              }
            }
          }
        }
      }
      SeqFeatAsnWrite (sfp, aop, atp);
      if (formal_name != NULL) {
        formal_name = GeneNomenclatureFree (formal_name);
      }
      if (rgp != NULL) {
        rgp = RNAGenFree (rgp);
      }
      SeqFeatFree (sfp);
    } else {
      AsnReadVal (aip, atp, &dv);
      AsnWrite (aop, atp, &dv);
      AsnKillValue (atp, &dv);
    }
  }
}

NLM_EXTERN void StripSeqFeatSupportAsnFilter (
  AsnIoPtr aip,
  AsnIoPtr aop
)

{
  AsnModulePtr       amp;
  AsnTypePtr         atp, atp_se, atp_sf;
  DataVal            dv;
  SeqFeatPtr         sfp;
  SeqFeatSupportPtr  support;

  if (aip == NULL || aop == NULL) return;

  amp = AsnAllModPtr ();
  if (amp == NULL) return;
  atp_se = AsnFind ("Seq-entry");
  atp_sf = AsnFind ("Seq-annot.data.ftable.E");
  if (atp_se == NULL || atp_sf == NULL) return;

  atp = atp_se;

  while ((atp = AsnReadId (aip, amp, atp)) != NULL) {
    if (atp == atp_sf) {
      sfp = SeqFeatAsnRead (aip, atp);
      support = sfp->support;
      sfp->support = NULL;
      SeqFeatAsnWrite (sfp, aop, atp);
      sfp->support = support;
      SeqFeatFree (sfp);
    } else {
      AsnReadVal (aip, atp, &dv);
      AsnWrite (aop, atp, &dv);
      AsnKillValue (atp, &dv);
    }
  }
}

/* CautiousSeqEntryCleanup section */

static Boolean EmptyOrNullString (CharPtr str)

{
  Char  ch;

  if (str == NULL) return TRUE;
  ch = *str;
  while (ch != '\0') {
    if (ch > ' ' && ch <= '~') return FALSE;
    str++;
    ch = *str;
  }
  return TRUE;
}

/* RemoveMultipleTitles currently removes FIRST title in chain */

static void RemoveMultipleTitles (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  SeqDescrPtr    descr = NULL;
  SeqDescrPtr    lasttitle = NULL;
  ObjValNodePtr  ovp;
  SeqDescrPtr    sdp;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    descr = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    descr = bssp->descr;
  } else return;
  for (sdp = descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_title) continue;
    if (lasttitle != NULL) {
      if (lasttitle->extended != 0) {
        ovp = (ObjValNodePtr) lasttitle;
        ovp->idx.deleteme = TRUE;
      }
      lasttitle = sdp;
    } else {
      lasttitle = sdp;
    }
  }
}

static void MakeBioSourceCopy (SeqEntryPtr sep, Pointer userdata)

{
  BioSourcePtr  biop;
  OrgRefPtr     master;
  OrgRefPtr     orp;
  SeqDescrPtr   sdp;

  master = (OrgRefPtr) userdata;
  sdp = SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL);
  if (sdp != NULL) return;
  biop = BioSourceNew ();
  if (biop == NULL) return;
  orp = OrgRefNew ();
  if (orp == NULL) return;
  biop->org = orp;
  orp->taxname = StringSave (master->taxname);
  orp->common = StringSave (master->common);
  sdp = CreateNewDescriptor (sep, Seq_descr_source);
  if (sdp == NULL) return;
  sdp->data.ptrvalue = (Pointer) biop;
}

static void ReplicatePopPhyMutSetBioSource (SeqEntryPtr sep)

{
  BioSourcePtr   biop;
  BioseqSetPtr   bssp;
  OrgRefPtr      orp;
  ObjValNodePtr  ovp;
  SeqDescrPtr    sdp;

  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;
  if (bssp->_class == 7 ||
      (bssp->_class >= 13 && bssp->_class <= 16)) {
    sdp = SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL);
    if (sdp == NULL) return;
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop == NULL) return;
    orp = biop->org;
    if (orp == NULL) return;
    VisitElementsInSep (sep, (Pointer) orp, MakeBioSourceCopy);
    if (sdp->extended != 0) {
      ovp = (ObjValNodePtr) sdp;
      ovp->idx.deleteme = TRUE;
    }
  }
}

static SeqFeatPtr BestCDS (SeqLocPtr loc, ValNodePtr cdslist)

{
  SeqFeatPtr  best_cds = NULL;
  SeqFeatPtr  cds;
  Int4        diff;
  Int4        min = INT4_MAX;
  ValNodePtr  vnp;

  if (loc == NULL || cdslist == NULL) return NULL;
  for (vnp = cdslist; vnp != NULL; vnp = vnp->next) {
    cds = (SeqFeatPtr) vnp->data.ptrvalue;
    diff = SeqLocAinB (loc, cds->location);
    if (diff >= 0) {
      if (diff < min) {
        min = diff;
        best_cds = cds;
      }
    }
  }
  return best_cds;
}

#define num_bond 5
static CharPtr feat_bond [num_bond] = {
  NULL,
  "disulfide bond",
  "thiolester bond",
  "xlink bond",
  "thioether bond"
};

#define num_site 27
static CharPtr feat_site [num_site] = {
  NULL, 
  "active", 
  "binding",
  "cleavage",
  "inhibit",
  "modified",
  "glycosylation",
  "myristoylation",
  "mutagenized",
  "metal-binding",
  "phosphorylation",
  "acetylation",
  "amidation",
  "methylation",
  "hydroxylation",
  "sulfatation",
  "oxidative-deamination",
  "pyrrolidone-carboxylic-acid",
  "gamma-carboxyglutamic-acid",
  "blocked",
  "lipid-binding",
  "np-binding",
  "dna-binding",
  "signal-peptide",
  "transit-peptide",
  "transmembrane-region",
  "nitrosylation"
};

static Int2 FindStr (CharPtr PNTR array, Int2 array_num, CharPtr str)

{
  Int2 i;

  for (i = 0; i < array_num; i++) {
    if (array [i] == NULL) continue;
    if (StringNCmp (str, array [i], StringLen (array [i])) == 0) return i;
  }
  return -1;
}

static SeqLocPtr fake_bond_loc (SeqLocPtr slp)

{
  SeqLocPtr loc, l, lnext, ldata;

  if (slp == NULL) return NULL;
  loc = MemNew (sizeof (SeqLoc));
  MemCopy (loc, slp, sizeof (SeqLoc));
  ldata = (SeqLocPtr) loc->data.ptrvalue;
  if (slp->choice != SEQLOC_MIX) return loc;
  for (l = ldata; l != NULL; l = lnext) {
    lnext = l->next;
    if (l->choice == SEQLOC_NULL) {
      ldata = remove_node (ldata, l);
    }
  }
  return loc;
}

static void ConvertImpFeatToProt (SeqFeatPtr feat, Pointer userdata)

{
  SeqFeatPtr  best_cds = NULL;
  Int2        bond = 0;
  BioseqPtr   bsp;
  ValNodePtr  cdslist;
  Uint1       choice = 0;
  Int4        frame;
  ImpFeatPtr  ifp;
  SeqLocPtr   loc;
  Uint1       processed = 0;
  ProtRefPtr  prp;
  SeqFeatPtr  sfp;
  SeqIdPtr    sip;
  Int2        site = 0;
  SeqLocPtr   slp;
  Uint1       subtype = 0;

  if (feat == NULL || feat->data.choice != SEQFEAT_IMP) return;
  ifp = (ImpFeatPtr) feat->data.value.ptrvalue;
  if (ifp == NULL) return;
  cdslist = (ValNodePtr) userdata;
  if (StringCmp (ifp->key, "mat_peptide") == 0) {
    processed = 2;
    choice = SEQFEAT_PROT;
    subtype = FEATDEF_mat_peptide_aa;
  } else if (StringCmp (ifp->key, "sig_peptide") == 0) {
    processed = 3;
    choice = SEQFEAT_PROT;
    subtype = FEATDEF_sig_peptide_aa;
  } else if (StringCmp (ifp->key, "transit_peptide") == 0) {
    processed = 4;
    choice = SEQFEAT_PROT;
    subtype = FEATDEF_transit_peptide_aa;
  } else if (StringCmp (ifp->key, "misc_feature") == 0 && feat->comment != NULL) {
    site = FindStr (feat_site, num_site, feat->comment);
    if (site != -1) {
      choice = SEQFEAT_SITE;
      subtype = FEATDEF_SITE;
    } else {
      bond = FindStr (feat_bond, num_bond, feat->comment);
      if (bond != -1) {
        choice = SEQFEAT_BOND;
        subtype = FEATDEF_BOND;
      }
    }
  }
  if (choice == 0) return;

  if (processed != 0 || site != 0) {
    best_cds = BestCDS (feat->location, cdslist);
  } else if (bond != 0) {
    loc = fake_bond_loc (feat->location);
    best_cds = BestCDS (loc, cdslist);
    SeqLocFree (loc);
  }
  if (best_cds == NULL) return;
  slp = dnaLoc_to_aaLoc (best_cds, feat->location, TRUE, &frame, FALSE);
  if (slp == NULL) return;
  sip = SeqLocId (best_cds->product);
  if (sip == NULL) return;
  bsp = BioseqLockById (sip);
  if (bsp == NULL) return;
  sfp = CreateNewFeatureOnBioseq (bsp, choice, slp);
  BioseqUnlock (bsp);
  if (sfp == NULL) return;

  sfp->partial = feat->partial;
  sfp->excpt = feat->excpt;
  sfp->exp_ev = feat->exp_ev;
  sfp->pseudo = feat->pseudo;

  sfp->comment = feat->comment;
  feat->comment = NULL;
  sfp->qual = feat->qual;
  feat->qual = NULL;
  sfp->title = feat->title;
  feat->title = NULL;
  sfp->ext = feat->ext;
  feat->ext = NULL;
  sfp->cit = feat->cit;
  feat->cit = NULL;

  sfp->xref = feat->xref;
  feat->xref = NULL;
  sfp->dbxref = feat->dbxref;
  feat->dbxref = NULL;
  sfp->except_text = feat->except_text;
  feat->except_text = NULL;

  if (choice == SEQFEAT_PROT) {
    prp = ProtRefNew ();
    sfp->data.value.ptrvalue = (Pointer) prp;
    if (prp != NULL) {
      prp->processed = processed;
    }
    switch (processed) {
    }
  } else if (choice == SEQFEAT_SITE) {
    sfp->data.value.intvalue = site;
  } else if (choice == SEQFEAT_BOND) {
    sfp->data.value.intvalue = bond;
  }
  sfp->idx.subtype = subtype;

  feat->idx.deleteme = TRUE;
}

static void GetListOfCDSs (SeqFeatPtr sfp, Pointer userdata)

{
  ValNodePtr PNTR  head;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;
  head = (ValNodePtr PNTR) userdata;
  ValNodeAddPointer (head, 0, sfp->data.value.ptrvalue);
}

static void ChangeImpFeatToProt (SeqEntryPtr sep)

{
  ValNodePtr  cdslist = NULL;

  VisitFeaturesInSep (sep, (Pointer) &cdslist, GetListOfCDSs);
  VisitFeaturesInSep (sep, (Pointer) cdslist, ConvertImpFeatToProt);
  ValNodeFree (cdslist);
}


NLM_EXTERN void MergeAdjacentAnnotsInList (SeqAnnotPtr sap)
{
  SeqAnnotPtr nextsap;
  SeqFeatPtr  sfp;

  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 1 && nextsap != NULL && nextsap->type == 1) {
      if (sap->id == NULL && nextsap->id == NULL &&
          sap->name == NULL && nextsap->name == NULL &&
          sap->db == 0 && nextsap->db == 0 &&
          sap->desc == NULL && nextsap->desc == NULL &&
          sap->data != NULL && nextsap->data != NULL) {
        sfp = (SeqFeatPtr) sap->data;
        while (sfp->next != NULL) {
          sfp = sfp->next;
        }
        sfp->next = (SeqFeatPtr) nextsap->data;
        nextsap->data = NULL;
        sap->next = nextsap->next;
        SeqAnnotFree (nextsap);
        nextsap = sap->next;
      }
    }
    sap = nextsap;
  }
}


static void MergeAdjacentAnnotsCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqAnnotPtr   sap;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  MergeAdjacentAnnotsInList (sap);
}

NLM_EXTERN Boolean PubIsEffectivelyEmpty (PubdescPtr pdp)

{
  ValNodePtr  vnp;

  if (pdp == NULL) return FALSE;
  vnp = pdp->pub;
  if (vnp != NULL && vnp->next == NULL && vnp->choice == PUB_Gen) {
    if (empty_citgen ((CitGenPtr) vnp->data.ptrvalue)) {
      return TRUE;
    }
  }
  return FALSE;
}

static void MarkEmptyDescsForCleanup (SeqDescrPtr sdp, Pointer userdata)

{
  GBBlockPtr     gbp;
  ObjValNodePtr  ovp;
  PubdescPtr     pdp;
  CharPtr        str;

  if (sdp == NULL || sdp->extended == 0) return;
  ovp = (ObjValNodePtr) sdp;
  if (sdp->choice == Seq_descr_title) {
    str = (CharPtr) sdp->data.ptrvalue;
    if (StringHasNoText (str)) {
      ovp->idx.deleteme = TRUE;
    }
  } else if (sdp->choice == Seq_descr_pub) {
    pdp = (PubdescPtr) sdp->data.ptrvalue;
    if (pdp == NULL) return;
    if (PubIsEffectivelyEmpty (pdp)) {
      ovp->idx.deleteme = TRUE;
    }
  } else if (sdp->choice == Seq_descr_genbank) {
    gbp = (GBBlockPtr) sdp->data.ptrvalue;
    if (gbp == NULL) return;
    /* gbp->source = MemFree (gbp->source); */
    /* gbp->origin = MemFree (gbp->origin); */
    gbp->taxonomy = MemFree (gbp->taxonomy);
    if (gbp->extra_accessions == NULL && gbp->source == NULL &&
        gbp->keywords == NULL && gbp->origin == NULL &&
        gbp->date == NULL && gbp->entry_date == NULL &&
        gbp->div == NULL && gbp->taxonomy == NULL) {
      ovp->idx.deleteme = TRUE;
    }
  }
}

static void MarkEmptyFeatsForCleanup (SeqFeatPtr sfp, Pointer userdata)

{
  GeneRefPtr  grp;
  PubdescPtr  pdp;
  ProtRefPtr  prp;
  ValNodePtr  vnp;

  if (sfp == NULL) return;
  if (sfp->data.choice == SEQFEAT_GENE && sfp->data.value.ptrvalue != NULL) {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (EmptyOrNullString (grp->locus)) {
      grp->locus = MemFree (grp->locus);
    }
    if (EmptyOrNullString (grp->allele)) {
      grp->allele = MemFree (grp->allele);
    }
    if (EmptyOrNullString (grp->desc)) {
      grp->desc = MemFree (grp->desc);
    }
    if (EmptyOrNullString (grp->maploc)) {
      grp->maploc = MemFree (grp->maploc);
    }
    if (EmptyOrNullString (grp->locus_tag)) {
      grp->locus_tag = MemFree (grp->locus_tag);
    }
    if (EmptyOrNullString (grp->locus) &&
        EmptyOrNullString (grp->allele) &&
        EmptyOrNullString (grp->desc) &&
        EmptyOrNullString (grp->maploc) &&
        EmptyOrNullString (grp->locus_tag) &&
        grp->db == NULL && grp->syn == NULL) {
      sfp->idx.deleteme = TRUE;
    }
  } else if (sfp->data.choice == SEQFEAT_PROT && sfp->data.value.ptrvalue != NULL) {
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    if (prp->processed != 3 && prp->processed != 4) {
      vnp = prp->name;
      if ((vnp == NULL || EmptyOrNullString ((CharPtr) vnp->data.ptrvalue)) &&
          EmptyOrNullString (prp->desc) &&
          prp->ec == NULL && prp->activity == NULL && prp->db == NULL) {
        sfp->idx.deleteme = TRUE;
      }
    }
  } else if (sfp->data.choice == SEQFEAT_PUB && sfp->data.value.ptrvalue != NULL) {
    pdp = (PubdescPtr) sfp->data.value.ptrvalue;
    if (PubIsEffectivelyEmpty (pdp)) {
      sfp->idx.deleteme = TRUE;
    }
  } else if (sfp->data.choice == SEQFEAT_COMMENT && EmptyOrNullString (sfp->comment)) {
    sfp->idx.deleteme = TRUE;
  }
}

static void ConvertPubFeatDescProc (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr      bsp;
  size_t         len;
  ObjValNodePtr  ovp;
  PubdescPtr     pdp;
  SeqDescPtr     sdp;
  SeqEntryPtr    sep;
  SeqIdPtr       sip;
  CharPtr        str;
  ValNode        vn;

  /* look for publication features */
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PUB) return;
  /* get bioseq by feature location */
  sip = SeqLocId (sfp->location);
  bsp = BioseqFind (sip);
  if (bsp == NULL) return;
  sip = SeqIdFindBest(bsp->id, 0);
  if (sip == NULL) return;
  vn.choice = SEQLOC_WHOLE;
  vn.extended = 0;
  vn.data.ptrvalue = (Pointer) sip;
  vn.next = NULL;
  /* is feature full length? */
  if (SeqLocCompare (sfp->location, &vn) != SLC_A_EQ_B) return;
  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) return;
  sdp = CreateNewDescriptor (sep, Seq_descr_pub);
  if (sdp == NULL) return;
  /* move publication from feature to descriptor */
  sdp->data.ptrvalue = sfp->data.value.ptrvalue;
  if (sdp->extended != 0) {
    ovp = (ObjValNodePtr) sdp;
    ovp->idx.subtype = Seq_descr_pub;
  }
  sfp->data.value.ptrvalue = NULL;
  /* flag old feature for removal */
  sfp->idx.deleteme = TRUE;
  /* move comment to remark */
  if (sfp->comment == NULL) return;
  pdp = (PubdescPtr) sdp->data.ptrvalue;
  if (pdp == NULL) return;
  if (pdp->comment == NULL) {
    pdp->comment = sfp->comment;
    sfp->comment = NULL;
  } else {
    len = StringLen (pdp->comment) + StringLen (sfp->comment) + 5;
    str = MemNew (sizeof (Char) * len);
    StringCpy (str, pdp->comment);
    StringCat (str, "; ");
    StringCat (str, sfp->comment);
    pdp->comment = MemFree (pdp->comment);
    pdp->comment = str;
  }
}

extern void ConvertSourceFeatDescProc (SeqFeatPtr sfp, Pointer userdata)

{
  BioSourcePtr   biop;
  BioseqPtr      bsp;
  SubSourcePtr   lastssp;
  ObjValNodePtr  ovp;
  SeqDescPtr     sdp;
  SeqEntryPtr    sep;
  SeqIdPtr       sip;
  SubSourcePtr   ssp;
  ValNode        vn;
  ValNodePtr     last_dbxref;

  /* look for biosource features */
  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC) return;
  /* get bioseq by feature location */
  sip = SeqLocId (sfp->location);
  bsp = BioseqFind (sip);
  if (bsp == NULL) return;
  sip = SeqIdFindBest(bsp->id, 0);
  if (sip == NULL) return;
  vn.choice = SEQLOC_WHOLE;
  vn.extended = 0;
  vn.data.ptrvalue = (Pointer) sip;
  vn.next = NULL;
  /* is feature full length? */
  if (SeqLocCompare (sfp->location, &vn) != SLC_A_EQ_B) return;
  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) return;
  sdp = CreateNewDescriptor (sep, Seq_descr_source);
  if (sdp == NULL) return;
  /* move biosource from feature to descriptor */
  sdp->data.ptrvalue = sfp->data.value.ptrvalue;
  if (sdp->extended != 0) {
    ovp = (ObjValNodePtr) sdp;
    ovp->idx.subtype = Seq_descr_source;
  }
  sfp->data.value.ptrvalue = NULL;
  /* flag old feature for removal */
  sfp->idx.deleteme = TRUE;
  /* move comment to subsource note */
  if (sfp->comment == NULL) return;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;
  ssp = SubSourceNew ();
  if (ssp == NULL) return;
  ssp->subtype = SUBSRC_other;
  ssp->name = sfp->comment;
  sfp->comment = NULL;
  /* link in at end, since BasicSeqEntry will have sorted this list */
  if (biop->subtype == NULL) {
    biop->subtype = ssp;
  } else {
    lastssp = biop->subtype;
    while (lastssp->next != NULL) {
      lastssp = lastssp->next;
    }
    lastssp->next = ssp;
  }

  /* move dbxrefs on feature to source */
  if (sfp->dbxref != NULL) {
    if (biop->org == NULL) {
      biop->org = OrgRefNew();
    }
    last_dbxref = biop->org->db;
    while (last_dbxref != NULL && last_dbxref->next != NULL) {
      last_dbxref = last_dbxref->next;
    }
    if (last_dbxref == NULL) {    
      biop->org->db = sfp->dbxref;
    } else {
      last_dbxref->next = sfp->dbxref;
    }
    sfp->dbxref = NULL;
  }
}

static void PromoteOrgRefDescToBioSource (SeqDescrPtr sdp, Pointer userdata)

{
  BioSourcePtr   biop;
  OrgRefPtr      orp;
  ObjValNodePtr  ovp;

  if (sdp->choice != Seq_descr_org) return;
  orp = (OrgRefPtr) sdp->data.ptrvalue;
  if (orp == NULL) return;
  biop = BioSourceNew ();
  if (biop == NULL) return;
  biop->org = orp;
  sdp->choice = Seq_descr_source;
  sdp->data.ptrvalue = (Pointer) biop;
  if (sdp->extended != 0) {
    ovp = (ObjValNodePtr) sdp;
    ovp->idx.subtype = Seq_descr_source;
  }
}

static void PromoteOrgRefFeatToBioSource (SeqFeatPtr sfp, Pointer userdata)

{
  BioSourcePtr  biop;
  OrgRefPtr     orp;

  if (sfp->data.choice != SEQFEAT_ORG) return;
  orp = (OrgRefPtr) sfp->data.value.ptrvalue;
  if (orp == NULL) return;
  biop = BioSourceNew ();
  if (biop == NULL) return;
  biop->org = orp;
  sfp->data.choice = SEQFEAT_BIOSRC;
  sfp->data.value.ptrvalue = (Pointer) biop;
  sfp->idx.subtype = FEATDEF_BIOSRC;
}

static void DeleteBadMarkedGeneXrefs (SeqFeatPtr sfp, Pointer userdata)

{
  SeqFeatXrefPtr       nextxref;
  SeqFeatXrefPtr PNTR  prevxref;
  Boolean              unlink;
  SeqFeatXrefPtr       xref;

  if (sfp == NULL) return;
  xref = sfp->xref;
  prevxref = (SeqFeatXrefPtr PNTR) &(sfp->xref);
  while (xref != NULL) {
    nextxref = xref->next;
    unlink = FALSE;
    if (xref->specialCleanupFlag && xref->data.choice == SEQFEAT_GENE) {
      if (SeqMgrGetOverlappingGene (sfp->location, NULL) != NULL) {
        unlink = TRUE;
      }
    }
    xref->specialCleanupFlag = FALSE;
    if (unlink) {
      *(prevxref) = xref->next;
      xref->next = NULL;
      SeqFeatXrefFree (xref);
    } else {
      prevxref = (SeqFeatXrefPtr PNTR) &(xref->next);
    }
    xref = nextxref;
  }
}

static void LookForMarkedGeneXrefs (SeqFeatPtr sfp, Pointer userdata)

{
  BoolPtr         hasMarkedGenes;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL || sfp->xref == NULL) return;
  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->specialCleanupFlag) {
      hasMarkedGenes = (BoolPtr) userdata;
      *hasMarkedGenes = TRUE;
      return;
    }
  }
}

NLM_EXTERN void CautiousSeqEntryCleanup (SeqEntryPtr sep, SeqEntryFunc taxfun, SeqEntryFunc taxmerge)

{
  /*
  Boolean      correct = FALSE;
  */
  Uint2        entityID;
  Boolean      hasMarkedGenes;
  ErrSev       lsev;
  ErrSev       msev;
  SeqEntryPtr  oldscope;
  /*
  Boolean      strip = TRUE;
  */
  Boolean      taxserver;

  if (sep == NULL) return;
  msev = ErrSetMessageLevel (SEV_MAX);
  lsev = ErrSetLogLevel (SEV_MAX);
  entityID = SeqMgrGetEntityIDForSeqEntry (sep);

  BasicSeqEntryCleanup (sep);

  VisitFeaturesInSep (sep, NULL, PromoteOrgRefFeatToBioSource);
  VisitDescriptorsInSep (sep, NULL, PromoteOrgRefDescToBioSource);

  oldscope = SeqEntrySetScope (sep);
  VisitFeaturesInSep (sep, NULL, ConvertSourceFeatDescProc);
  VisitFeaturesInSep (sep, NULL, ConvertPubFeatDescProc);
  SeqEntrySetScope (oldscope);

  VisitFeaturesInSep (sep, NULL, MarkEmptyFeatsForCleanup);
  VisitDescriptorsInSep (sep, NULL, MarkEmptyDescsForCleanup);
  DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);

  SeqEntryExplore (sep, NULL, MergeAdjacentAnnotsCallback);

  ChangeImpFeatToProt (sep);
  DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);

  VisitBioseqsInSep (sep, NULL, ExtendSingleGeneOnMRNA);

  ReplicatePopPhyMutSetBioSource (sep);
  SeqEntryExplore (sep, NULL, RemoveMultipleTitles);

  /* LoopSeqEntryToAsn3 section here */
  taxserver = (Boolean) (taxfun != NULL || taxmerge != NULL);

  /*
  if (correct) {
    SeqEntryExplore(sep, (Pointer)(&porg), CorrectSourceFeat);
  }
  */







  /* a few more things to do here */

  hasMarkedGenes = FALSE;
  VisitFeaturesInSep (sep, (Pointer) &hasMarkedGenes, LookForMarkedGeneXrefs);
  if (hasMarkedGenes) {
    SeqMgrIndexFeatures (entityID, NULL);
    VisitFeaturesInSep (sep, NULL, DeleteBadMarkedGeneXrefs);
    SeqMgrClearFeatureIndexes (entityID, NULL);
  }

  BasicSeqEntryCleanup (sep);

  AssignIDsInEntity (entityID, 0, NULL);

  ErrSetMessageLevel (msev);
  ErrSetLogLevel (lsev);
}

/*
static Int4 LoopSeqEntryToAsn3 (SeqEntryPtr sep, Boolean strip, Boolean correct, SeqEntryFunc taxfun, SeqEntryFunc taxmerge)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   oldscope;
  Int4          rsult;
  Boolean       taxserver;

  rsult = 0;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        rsult += LoopSeqEntryToAsn3 (sep, strip, correct, taxfun, taxmerge);
      }
      return rsult;
    }
  }
  oldscope = SeqEntrySetScope (sep);
  taxserver = (Boolean) (taxfun != NULL || taxmerge != NULL);
  rsult = SeqEntryToAsn3Ex (sep, strip, correct, taxserver, taxfun, taxmerge);
  SeqEntrySetScope (oldscope);
  return rsult;
}
  LoopSeqEntryToAsn3 (sep, TRUE, FALSE, taxfun, taxmerge);

*/

typedef struct featdefNameStruct {
  Uint1    type;
  CharPtr  name;
} FeatdefNameData, PNTR FeatdefNamePtr;

static FeatdefNameData featdefWithName [] = {
  { FEATDEF_10_signal ,          "-10_signal"         },
  { FEATDEF_35_signal ,          "-35_signal"         },
  { FEATDEF_3clip ,              "3'clip"             },
  { FEATDEF_3UTR ,               "3'UTR"              },
  { FEATDEF_5clip ,              "5'clip"             },
  { FEATDEF_5UTR ,               "5'UTR"              },
  { FEATDEF_attenuator ,         "attenuator"         },
  { FEATDEF_BOND ,               "Bond"               },
  { FEATDEF_CAAT_signal ,        "CAAT_signal"        },
  { FEATDEF_CDS ,                "CDS"                },
  { FEATDEF_CLONEREF ,           "CloneRef"           },
  { FEATDEF_PUB ,                "Cit"                },
  { FEATDEF_COMMENT ,            "Comment"            },
  { FEATDEF_conflict ,           "conflict"           },
  { FEATDEF_C_region ,           "C_region"           },
  { FEATDEF_D_loop ,             "D-loop"             },
  { FEATDEF_D_segment ,          "D_segment"          },
  { FEATDEF_enhancer ,           "enhancer"           },
  { FEATDEF_exon ,               "exon"               },
  { FEATDEF_gap ,                "gap"                },
  { FEATDEF_GC_signal ,          "GC_signal"          },
  { FEATDEF_GENE ,               "Gene"               },
  { FEATDEF_HET ,                "Het"                },
  { FEATDEF_iDNA ,               "iDNA"               },
  { FEATDEF_IMP ,                "Import"             },
  { FEATDEF_Imp_CDS ,            "Imp_CDS"            },
  { FEATDEF_intron ,             "intron"             },
  { FEATDEF_J_segment ,          "J_segment"          },
  { FEATDEF_LTR ,                "LTR"                },
  { FEATDEF_mat_peptide_aa ,     "mat_peptide"        },
  { FEATDEF_mat_peptide ,        "mat_peptide_nt"     },
  { FEATDEF_misc_binding ,       "misc_binding"       },
  { FEATDEF_misc_difference ,    "misc_difference"    },
  { FEATDEF_misc_feature ,       "misc_feature"       },
  { FEATDEF_misc_recomb ,        "misc_recomb"        },
  { FEATDEF_otherRNA ,           "misc_RNA"           },
  { FEATDEF_misc_signal ,        "misc_signal"        },
  { FEATDEF_misc_structure ,     "misc_structure"     },
  { FEATDEF_mobile_element ,     "mobile_element"     },
  { FEATDEF_modified_base ,      "modified_base"      },
  { FEATDEF_mRNA ,               "mRNA"               },
  { FEATDEF_NON_STD_RESIDUE ,    "NonStdRes"          },
  { FEATDEF_NUM ,                "Num"                },
  { FEATDEF_N_region ,           "N_region"           },
  { FEATDEF_ncRNA ,              "ncRNA"              },
  { FEATDEF_old_sequence ,       "old_sequence"       },
  { FEATDEF_operon ,             "operon"             },
  { FEATDEF_oriT ,               "oriT"               },
  { FEATDEF_polyA_signal ,       "polyA_signal"       },
  { FEATDEF_polyA_site ,         "polyA_site"         },
  { FEATDEF_preRNA ,             "precursor_RNA"      },
  { FEATDEF_preprotein ,         "preprotein"         },
  { FEATDEF_primer_bind ,        "primer_bind"        },
  { FEATDEF_prim_transcript ,    "prim_transcript"    },
  { FEATDEF_promoter ,           "promoter"           },
  { FEATDEF_PROT ,               "Protein"            },
  { FEATDEF_protein_bind ,       "protein_bind"       },
  { FEATDEF_RBS ,                "RBS"                },
  { FEATDEF_REGION ,             "Region"             },
  { FEATDEF_repeat_region ,      "repeat_region"      },
  { FEATDEF_repeat_unit ,        "repeat_unit"        },
  { FEATDEF_rep_origin ,         "rep_origin"         },
  { FEATDEF_rRNA ,               "rRNA"               },
  { FEATDEF_RSITE ,              "Rsite"              },
  { FEATDEF_satellite ,          "satellite"          },
  { FEATDEF_scRNA ,              "scRNA"              },
  { FEATDEF_PSEC_STR ,           "SecStr"             },
  { FEATDEF_sig_peptide_aa ,     "sig_peptide"        },
  { FEATDEF_sig_peptide ,        "sig_peptide_nt"     },
  { FEATDEF_SITE ,               "Site"               },
  { FEATDEF_site_ref ,           "Site-ref"           },
  { FEATDEF_snoRNA ,             "snoRNA"             },
  { FEATDEF_snRNA ,              "snRNA"              },
  { FEATDEF_source ,             "source"             },
  { FEATDEF_BIOSRC ,             "Src"                },
  { FEATDEF_stem_loop ,          "stem_loop"          },
  { FEATDEF_STS ,                "STS"                },
  { FEATDEF_S_region ,           "S_region"           },
  { FEATDEF_TATA_signal ,        "TATA_signal"        },
  { FEATDEF_terminator ,         "terminator"         },
  { FEATDEF_tmRNA ,              "tmRNA"              },
  { FEATDEF_transit_peptide_aa , "transit_peptide"    },
  { FEATDEF_transit_peptide ,    "transit_peptide_nt" },
  { FEATDEF_tRNA ,               "tRNA"               },
  { FEATDEF_TXINIT ,             "TxInit"             },
  { FEATDEF_unsure ,             "unsure"             },
  { FEATDEF_USER ,               "User"               },
  { FEATDEF_variation ,          "variation"          },
  { FEATDEF_VARIATIONREF ,       "VariationRef"       },
  { FEATDEF_virion ,             "virion"             },
  { FEATDEF_V_region ,           "V_region"           },
  { FEATDEF_V_segment ,          "V_segment"          },
  { FEATDEF_SEQ ,                "Xref"               }
};

NLM_EXTERN Uint1 FindFeatDefTypeFromKey (CharPtr key)

{
  Int2  L, R, mid;

  if (key == NULL || *key == '\0') return FEATDEF_BAD;

  L = 0;
  R = (sizeof (featdefWithName) / sizeof (FeatdefNameData)) - 1;

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (featdefWithName [mid].name, key) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (featdefWithName [R].name, key) == 0) {
    return featdefWithName [R].type;
  }

  return FEATDEF_BAD;
}

static CharPtr featurekeys [] = {
  "???" ,
  "Gene" ,
  "Org" ,
  "CDS" ,
  "Protein" ,
  "precursor_RNA" ,
  "mRNA" ,
  "tRNA" ,
  "rRNA" ,
  "snRNA" ,
  "scRNA" ,
  "misc_RNA" ,
  "Cit" ,
  "Xref" ,
  "Import" ,
  "allele" ,
  "attenuator" ,
  "C_region" ,
  "CAAT_signal" ,
  "CDS" ,
  "conflict" ,
  "D-loop" ,
  "D_segment" ,
  "enhancer" ,
  "exon" ,
  "GC_signal" ,
  "iDNA" ,
  "intron" ,
  "J_segment" ,
  "LTR" ,
  "mat_peptide" ,
  "misc_binding" ,
  "misc_difference" ,
  "misc_feature" ,
  "misc_recomb" ,
  "misc_RNA" ,
  "misc_signal" ,
  "misc_structure" ,
  "modified_base" ,
  "mutation" ,
  "N_region" ,
  "old_sequence" ,
  "polyA_signal" ,
  "polyA_site" ,
  "precursor_RNA" ,
  "prim_transcript" ,
  "primer_bind" ,
  "promoter" ,
  "protein_bind" ,
  "RBS" ,
  "repeat_region" ,
  "repeat_unit" ,
  "rep_origin" ,
  "S_region" ,
  "satellite" ,
  "sig_peptide" ,
  "source" ,
  "stem_loop" ,
  "STS" ,
  "TATA_signal" ,
  "terminator" ,
  "transit_peptide" ,
  "unsure" ,
  "V_region" ,
  "V_segment" ,
  "variation" ,
  "virion" ,
  "3'clip" ,
  "3'UTR" ,
  "5'clip" ,
  "5'UTR" ,
  "-10_signal" ,
  "-35_signal" ,
  "Site-ref" ,
  "Region" ,
  "Comment" ,
  "Bond" ,
  "Site" ,
  "Rsite" ,
  "User" ,
  "TxInit" ,
  "Num" ,
  "SecStr" ,
  "NonStdRes" ,
  "Het" ,
  "Src" ,
  "proprotein" ,
  "mat_peptide" ,
  "sig_peptide" ,
  "transit_peptide",
  "snoRNA",
  "gap",
  "operon",
  "oriT",
  "ncRNA",
  "tmRNA",
  "CloneRef",
  "VariationRef",
  "mobile_element"
};

NLM_EXTERN CharPtr FindKeyFromFeatDefType (Uint1 type, Boolean forGBFF)

{
  CharPtr  key;

  if (type < FEATDEF_GENE || type >= FEATDEF_MAX) {
    type = FEATDEF_BAD;
  }
  key = featurekeys [type];

  if (forGBFF) {
    if (type == FEATDEF_GENE) {
      key = "gene";
    } else if (type == FEATDEF_REGION ||
               type == FEATDEF_COMMENT ||
               type == FEATDEF_BOND ||
               type == FEATDEF_SITE) {
      key = "misc_feature";
    } else if (type == FEATDEF_VARIATIONREF) {
      key = "variation";
    }
  }

  return key;
}

/* tRNA codon index to codon string lookup table functions */

typedef struct gcCodonStruct {
  Uint1    index;
  CharPtr  codon;
} GcCodonData, PNTR GcCodonPtr;

static CharPtr    gcCodonStrings = NULL;
static GcCodonPtr codonGcIndex = NULL;

/* mapping from NCBI2na to codon codes */

static Uint1 codon_xref [4] = {
  2,  /* A */
  1,  /* C */
  3,  /* G */
  0   /* T */
};

static int LIBCALLBACK SortCodonByString (
  VoidPtr vp1,
  VoidPtr vp2
)

{
  int         compare;
  GcCodonPtr  gcp1 = vp1;
  GcCodonPtr  gcp2 = vp2;

  if (gcp1 == NULL || gcp2 == NULL) return 0;

  compare = StringICmp (gcp1->codon, gcp2->codon);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  return 0;
}

static void InitGcCodons (void)

{
  Uint1           codon [4], index;
  GcCodonPtr      codonGcIdx;
  CharPtr         gcCodonStr;
  Int2            i, j, k;
  int             idx, offset;
  CharPtr         ptr;
  Uint1           residue;
  SeqMapTablePtr  smtp;

  if (codonGcIndex != NULL && gcCodonStrings != NULL) return;

  gcCodonStr = (CharPtr) MemNew (sizeof (Char) * 256);
  if (gcCodonStr == NULL) return;
  codonGcIdx = (GcCodonPtr) MemNew (sizeof (GcCodonData) * 64);
  if (codonGcIdx == NULL) return;

  smtp = SeqMapTableFind (Seq_code_iupacna, Seq_code_ncbi2na);
  if (smtp == NULL) return;

  for (idx = 0; idx < 64; idx++) {
    index = (Uint1) idx;

    for (i = 0, j = 16; i < 3; i++, j /= 4) {
      residue = (Uint1) ((Int2) index / j);
      index -= (Uint1) (residue * j);
      for (k = 0; k < 4; k++) {
        if (codon_xref [k] == residue) {
          residue = (Uint1) k;
          break;
        }
      }
      residue = SeqMapTableConvert (smtp, residue);
      codon [i] = residue;
    }
    codon [3] = 0;

    offset = 4 * idx;
    ptr = gcCodonStr + offset;
    StringCpy (ptr, (CharPtr) codon);

    codonGcIdx [idx].index = (Uint1) idx;
    codonGcIdx [idx].codon = ptr;
  }

  StableMergeSort (codonGcIdx, (size_t) 64, sizeof (GcCodonData), SortCodonByString);

  gcCodonStrings = gcCodonStr;
  codonGcIndex = codonGcIdx;
}

NLM_EXTERN Uint1 CodonToGcIndex (CharPtr codon)

{
  Char  ch;
  Int2  i, L, R, mid;
  Char  tmp [4];

  if (codonGcIndex == NULL) {
    InitGcCodons ();
  }
  if (codonGcIndex == NULL) return 255;
  if (StringLen (codon) != 3) return 255;
  StringNCpy_0 (tmp, codon, sizeof (tmp));

  for (i = 0; i < 3; i++) {
    ch = tmp [i];
    ch = TO_UPPER (ch);
    if (ch == 'U') {
       ch = 'T';
    }
    tmp [i] = ch;
  }

  L = 0;
  R = 63;

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (codonGcIndex [mid].codon, tmp) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (codonGcIndex [R].codon, tmp) == 0) {
    return codonGcIndex [R].index;
  }

  return 255;
}

NLM_EXTERN CharPtr GcIndextoCodon (Uint1 index)

{
  int      offset;
  CharPtr  ptr;

  if (gcCodonStrings == NULL) {
    InitGcCodons ();
  }
  if (gcCodonStrings == NULL) return NULL;
  if (index > 63) return NULL;

  offset = 4 * index;
  ptr = gcCodonStrings + offset;

  return ptr;
}

static FloatHi GetCddBitScore (SeqFeatPtr sfp)

{
  ObjectIdPtr    oip;
  UserFieldPtr   ufp;
  UserObjectPtr  uop;

  if (sfp == NULL) return 0.0;
  uop = sfp->ext;
  if (uop == NULL) return 0.0;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "cddScoreData") != 0) return 0.0;
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (oip != NULL && StringICmp (oip->str, "bit_score") == 0) {
      if (ufp->choice == 3) {
        return ufp->data.realvalue;
      }
    }
  }
  return 0.0;
}

static Boolean FeatIsCDD (
  SeqFeatPtr sfp,
  FloatHi PNTR scoreP
)

{
  DbtagPtr    dbt;
  ValNodePtr  vnp;

  if (scoreP != NULL) {
    *scoreP = 0.0;
  }
  for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      if (StringCmp (dbt->db, "CDD") == 0 || StringCmp (dbt->db, "cdd") == 0) {
        if (scoreP != NULL) {
          *scoreP = GetCddBitScore (sfp);
        }
        return TRUE;
      }
    }
  }

  return FALSE;
}
static void BestCDDperBioseq (BioseqPtr bsp, Pointer userdata)

{
  SeqFeatPtr         best;
  SeqMgrFeatContext  context;
  FloatHi            currscore;
  Int4               right;
  SeqFeatPtr         sfp;
  FloatHi            topscore;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
  while (sfp != NULL) {
    if (context.featdeftype == FEATDEF_REGION && FeatIsCDD (sfp, &currscore)) {
      best = sfp;
      right = context.right;
      topscore = currscore;
      sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
      while (sfp != NULL && context.featdeftype == FEATDEF_REGION &&
             FeatIsCDD (sfp, &currscore) && context.left < right) {
        right = MAX (context.right, right);
        if (currscore <= topscore) {
          sfp->idx.deleteme = TRUE;
        } else {
          best->idx.deleteme = TRUE;
          best = sfp;
          topscore = currscore;
        }
        sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
      }
    } else {
      sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
    }
  }
}

NLM_EXTERN void LeaveBestCDD (SeqEntryPtr sep)

{
  Uint2  entityID;

  if (sep == NULL) return;
  entityID = ObjMgrGetEntityIDForChoice (sep);
  if (entityID < 1) return;

  if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
    SeqMgrIndexFeatures (entityID, NULL);
  }

  VisitBioseqsInSep (sep, NULL, BestCDDperBioseq);
  DeleteMarkedObjects (entityID, 0, NULL);

  SeqMgrClearFeatureIndexes (entityID, NULL);
}

static CharPtr CompressNonBases (CharPtr str)

{
  Char     ch;
  CharPtr  dst;
  CharPtr  ptr;

  if (str == NULL || str [0] == '\0') return NULL;

  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_ALPHA (ch)) {
      *dst = ch;
      dst++;
    }
    ptr++;
    ch = *ptr;
  }
  *dst = '\0';

  return str;
}

static void LIBCALLBACK SPStreamToRaw (
  CharPtr sequence,
  Pointer userdata
)

{
  ByteStorePtr  bs;
  Char          ch;
  size_t        len;
  CharPtr       tmp;

  bs = (ByteStorePtr) userdata;
  tmp = sequence;
  ch = *tmp;
  while (ch != '\0') {
    if (ch == '\n' || ch == '\r' || ch == '\t') {
      *tmp = ' ';
    } else {
      *tmp = TO_UPPER (ch);
    }
    tmp++;
    ch = *tmp;
  }
  TrimSpacesAroundString (sequence);
  CompressNonBases (sequence);

  len = StringLen (sequence);
  if (len < 1) return;
  BSWrite (bs, sequence, len * sizeof (Char));
}

NLM_EXTERN void SegOrDeltaBioseqToRaw (BioseqPtr bsp)

{
  ByteStorePtr  bs;

  if (bsp == NULL || (bsp->repr != Seq_repr_seg && bsp->repr != Seq_repr_delta)) return;
  if (! ISA_na (bsp->mol)) return;
  bs = BSNew (bsp->length);
  if (bs == NULL) return;

  SeqPortStream (bsp, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) bs, SPStreamToRaw);

  if (bsp->repr == Seq_repr_seg && bsp->seq_ext_type == 1) {
    bsp->seq_ext = SeqLocSetFree ((ValNodePtr) bsp->seq_ext);
    bsp->seq_ext_type = 0;
  } else if (bsp->repr == Seq_repr_delta && bsp->seq_ext_type == 4) {
    bsp->seq_ext = NULL; /* for now just NULL out */
    bsp->seq_ext_type = 0;
  }
  bsp->seq_data = SeqDataFree (bsp->seq_data, bsp->seq_data_type);
  bsp->seq_data = (SeqDataPtr) bs;
  bsp->length = BSLen (bs);
  bsp->repr = Seq_repr_raw;
  bsp->seq_data_type = Seq_code_iupacna;
}

typedef struct segtodelta
{
  ValNodePtr seq_ext;
  Int4       len;
  SeqIdPtr   master_sip;
  BioseqPtr  master_bsp;  
  Int4       num_segs_converted;
} SegToDeltaData, PNTR SegToDeltaPtr;


static ValNodePtr CombineDescriptorLists (ValNodePtr target, ValNodePtr insert)
{
  ValNodePtr combined_list = NULL;
  ValNodePtr vnp, vnp_next;
  ValNodePtr title_descr = NULL, prev_descr = NULL;
  CharPtr    combined_title;
  Int4       combined_title_len;
  
  if (target == NULL)
  {
    combined_list = insert;
  }
  else if (insert == NULL)
  {
    combined_list = target;
  }
  else
  {
    combined_list = target;
      for (vnp = target; vnp->next != NULL; vnp = vnp->next)
      {
        if (vnp->choice == Seq_descr_title)
        {
          title_descr = vnp;
        }
      }
      prev_descr = vnp;
      if (title_descr == NULL)
      {
        prev_descr->next = insert;
      }
      else
      {
        for (vnp = insert; vnp != NULL; vnp = vnp_next)
        {
          vnp_next = vnp->next;
          vnp->next = NULL;
          if (vnp->choice == Seq_descr_title)
          {
            /* combine with previous title */
            combined_title_len = StringLen (title_descr->data.ptrvalue)
                                + StringLen (vnp->data.ptrvalue)
                                + 3;
            combined_title = (CharPtr) MemNew (sizeof (Char) * combined_title_len);
            if (combined_title != NULL)
            {
              StringCpy (combined_title, title_descr->data.ptrvalue);
              StringCat (combined_title, "; ");
              StringCat (combined_title, vnp->data.ptrvalue);
              title_descr->data.ptrvalue = MemFree (title_descr->data.ptrvalue);
              title_descr->data.ptrvalue = combined_title;
            }
            ValNodeFreeData (vnp);
          }
          else
          {
            /* add to master list */
            prev_descr->next = vnp;
            prev_descr = vnp;
          }
        } 
      }
  }
  return combined_list;
}

static void MoveSegmentLocToMaster (SeqLocPtr slp, SegToDeltaPtr sdp)
{
  SeqIntPtr     sintp;
  SeqLocPtr     slp2;
  SeqPntPtr     spp;
  PackSeqPntPtr pspp;
  Int4          i;
  
  if (slp == NULL || sdp == NULL) return;
  
  switch (slp->choice)
  {
    case SEQLOC_WHOLE:
    case SEQLOC_EMPTY:
      slp->data.ptrvalue = SeqIdFree (slp->data.ptrvalue);
      slp->data.ptrvalue = SeqIdDup (sdp->master_sip);
      break;
    case SEQLOC_INT:
      sintp = (SeqIntPtr) slp->data.ptrvalue;
      if (sintp != NULL)
      {
        sintp->id = SeqIdFree (sintp->id);
        sintp->id = SeqIdDup (sdp->master_sip);
        sintp->from += sdp->len;
        sintp->to += sdp->len;
        /* strand stays the same */
      }
      break;
    case SEQLOC_PACKED_INT:
    case SEQLOC_MIX:    
    case SEQLOC_EQUIV:
            slp2 = (SeqLocPtr)slp->data.ptrvalue;
            while (slp2 != NULL)
            {
                MoveSegmentLocToMaster (slp2, sdp);
                slp2 = slp2->next;
            }
      break;
    case SEQLOC_PNT:
      spp = (SeqPntPtr) slp->data.ptrvalue;
      if (spp != NULL)
      {
        spp->id = SeqIdFree (spp->id);
        spp->id = SeqIdDup (sdp->master_sip);
        spp->point += sdp->len;
      }
      break;
    case SEQLOC_PACKED_PNT:
            pspp = (PackSeqPntPtr)slp->data.ptrvalue;
            while (pspp != NULL)
      {
        for (i = 0; i < pspp->used; i++)
        {
          pspp->pnts[i] += sdp->len;
        }
        pspp->id = SeqIdFree (pspp->id);
        pspp->id = SeqIdDup (sdp->master_sip);
        pspp = pspp->next;
      }
      break;
  }
}

static void MoveSegmentFeaturesToMaster (SeqFeatPtr sfp, Pointer userdata)

{
  SegToDeltaPtr   segdeltptr;

  if (sfp == NULL || userdata == NULL) return;
  
  segdeltptr = (SegToDeltaPtr) userdata;

  MoveSegmentLocToMaster (sfp->location, segdeltptr);
}

#if 0
static void AdjustAlignmentOffsetsForDeltaConversion (SeqAlignPtr salp, Int4Ptr offsets, BoolPtr is_gap, Int4 num_sets)
{
  DenseSegPtr dsp;
  Int4        aln_seg_num, j, index;
  
  if (salp == NULL || offsets == NULL) return;

  /* adjust alignment starts to match delta sequence coordinates */
  if (salp->segtype == 2)
  {
    dsp = (DenseSegPtr) (salp->segs);
    aln_seg_num = 0;
    for (j = 0; j < num_sets; j++)
    {
      if (!is_gap [j])
      {
        for (index = 0; index < dsp->numseg; index++)
        {
          if (dsp->starts [dsp->dim * index + aln_seg_num] != -1)
          {
            dsp->starts [dsp->dim * index + aln_seg_num] += offsets [j];  
          }
        }
        aln_seg_num++;
      }
    }
  }      
}
#endif

static SeqAnnotPtr CombineAnnots (SeqAnnotPtr target, SeqAnnotPtr insert, Int4 offset)
{
  SeqAnnotPtr combined_list = NULL;
  SeqAnnotPtr feature_sap = NULL;
  SeqAnnotPtr prev_sap = NULL;
  SeqAnnotPtr sap, next_sap;
  SeqFeatPtr  last_feat, first_feat;
  
  if (target == NULL)
  {
    combined_list = insert;
  }
  else if (insert == NULL)
  {
    combined_list = target;
  }
  else
  {
    combined_list = target;
    for (sap = target; sap != NULL; sap = sap->next)
    {
      if (sap->type == 1 && sap->name == NULL && sap->desc == NULL)
      {
        feature_sap = sap;
      }
      prev_sap = sap;
    }
    for (sap = insert; sap != NULL; sap = next_sap)
    {
      next_sap = sap->next;
      sap->next = NULL;
      if (sap->type == 1 && sap->name == NULL && sap->desc == NULL && feature_sap != NULL)
      {
        first_feat = (SeqFeatPtr) sap->data;
        if (first_feat != NULL)
        {
          for (last_feat = (SeqFeatPtr) feature_sap->data;
               last_feat != NULL && last_feat->next != NULL;
               last_feat = last_feat->next)
          {  
          }
          if (last_feat == NULL)
          {
            feature_sap->data = first_feat;    
          }
          else
          {
            last_feat->next = first_feat;
          }
        }
        sap->data = NULL;
        SeqAnnotFree (sap);
      }
      else
      {
        prev_sap->next = sap;
        prev_sap = sap;
      }
    }
  }
  return combined_list;
}

static Int4 AddGapSeqLit (ValNodePtr PNTR seq_ext)
{
  SeqLitPtr       slip;
  IntFuzzPtr      ifp;
  CharPtr         gap_chars = "NNNNNNNNNN"
                              "NNNNNNNNNN"
                              "NNNNNNNNNN"
                              "NNNNNNNNNN"
                              "NNNNNNNNNN"
                              "NNNNNNNNNN"
                              "NNNNNNNNNN"
                              "NNNNNNNNNN"
                              "NNNNNNNNNN"
                              "NNNNNNNNNN";
                              
  if (seq_ext == NULL) return 0;
                                
  slip = (SeqLitPtr) MemNew (sizeof (SeqLit));
  if (slip != NULL) {
    slip->length = 100;
    ValNodeAddPointer (seq_ext, (Int2) 2, (Pointer) slip);
    ifp = IntFuzzNew ();
    ifp->choice = 4;
      
    slip->fuzz = ifp;
    slip->seq_data = (SeqDataPtr) BSNew (slip->length);
    slip->seq_data_type = Seq_code_iupacna;
    AddBasesToByteStore ((ByteStorePtr) slip->seq_data, gap_chars);
    return 100;
  }
  return 0;
}

static Boolean LIBCALLBACK 
AddSegmentToDeltaSeq 
(SeqLocPtr slp,
 SeqMgrSegmentContextPtr context)

{
  SegToDeltaPtr   segdeltptr;
  SeqIdPtr        sip;
  BioseqPtr       bsp;
  CharPtr         bases;
  SeqLitPtr       slip;

  SeqLocPtr         loc;

  if (slp == NULL || context == NULL) return FALSE;
  segdeltptr = (SegToDeltaPtr) context->userdata;
  if (segdeltptr == NULL) return FALSE;

  sip = SeqLocId (slp);
  
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) 
  {
    return TRUE;
  }

  bsp = BioseqFind (sip);

  if (bsp == NULL)
  {
    return TRUE;
  }
  
  bases = GetSequenceByBsp (bsp);
  if (bases == NULL) 
  {
    bsp->idx.deleteme = TRUE;
    return TRUE;    
  }
  
  if (segdeltptr->seq_ext != NULL)
  {
    /* insert gap of unknown length between the previous segment
     * and this one.
     */
    segdeltptr->len += AddGapSeqLit (&(segdeltptr->seq_ext));
  }

  /* move descriptors to master_bsp */
  segdeltptr->master_bsp->descr = CombineDescriptorLists (segdeltptr->master_bsp->descr, bsp->descr);
  bsp->descr = NULL;
  
  /* move features to master_bsp */
  VisitFeaturesOnBsp (bsp, segdeltptr, MoveSegmentFeaturesToMaster);
  segdeltptr->master_bsp->annot = CombineAnnots (segdeltptr->master_bsp->annot, bsp->annot, segdeltptr->len);
  bsp->annot = NULL;
  
  slip = (SeqLitPtr) MemNew (sizeof (SeqLit));
  if (slip != NULL) 
  {
    slip->length = StringLen (bases);
    ValNodeAddPointer (&(segdeltptr->seq_ext), (Int2) 2, (Pointer) slip);
    slip->seq_data = (SeqDataPtr) BSNew (slip->length);
    slip->seq_data_type = Seq_code_iupacna;
    AddBasesToByteStore ((ByteStorePtr) slip->seq_data, bases);
    segdeltptr->len += slip->length;
  }

  segdeltptr->num_segs_converted ++;
  return TRUE;
}

static BioseqPtr GetDeltaSeqFromMasterSeg (BioseqPtr bsp)
{
  BioseqPtr      new_bsp;
  SegToDeltaData sdd;
  BioseqSetPtr   segset;
  
  if (bsp == NULL || bsp->repr != Seq_repr_seg 
      || bsp->seq_ext == NULL || bsp->seq_ext_type != 1) 
  {
    return NULL;
  }
  
  if (! ISA_na (bsp->mol)) return NULL;

  /* use SeqMgrExploreSegments to build a list of SeqLitPtr */
  sdd.seq_ext = NULL;
  sdd.len = 0;
  sdd.master_bsp = bsp;
  sdd.master_sip = bsp->id;
  sdd.num_segs_converted = 0;
  
  /* move descriptors and features from segset to master seg */
  if (bsp->idx.parenttype == OBJ_BIOSEQSET)
  {
    segset = (BioseqSetPtr) bsp->idx.parentptr;
    if (segset != NULL)
    {
      bsp->descr = CombineDescriptorLists (bsp->descr, segset->descr);
      segset->descr = NULL;
    }
  }  

  SeqMgrExploreSegments (bsp, (Pointer) &sdd, AddSegmentToDeltaSeq);
  
  new_bsp = BioseqNew ();
  new_bsp->descr = bsp->descr;
  bsp->descr = NULL;
  new_bsp->annot = bsp->annot;
  bsp->annot = NULL;
  new_bsp->seq_data = NULL;
  new_bsp->seq_data_type = 0;
  new_bsp->repr = Seq_repr_delta;
  new_bsp->seq_ext_type = 4;
  new_bsp->seq_ext = sdd.seq_ext;
  new_bsp->length = sdd.len;
  new_bsp->id = SeqIdDup (bsp->id); 
/*  new_bsp->id = MakeUniqueSeqID ("delta_"); */
  new_bsp->mol = bsp->mol;

  BioseqPack (new_bsp);  
  return new_bsp;
}

NLM_EXTERN void ConvertSegSetsToDeltaSequences (SeqEntryPtr sep)
{
  BioseqSetPtr  bssp;
  SeqEntryPtr   sub_sep, prev_sep, next_sep;
  ObjMgrDataPtr omdptop;
  ObjMgrData    omdata;
  Uint2         parenttype;
  Pointer       parentptr;
  SeqEntryPtr   new_sep;
  BioseqPtr     bsp, new_bsp = NULL;
  BioseqSetPtr  parent_set;
  
  if (sep == NULL || !IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp->_class == 2)
  {
    SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
    GetSeqEntryParent (sep, &parentptr, &parenttype);
  
    parent_set = (BioseqSetPtr)(bssp->idx.parentptr);
    prev_sep = NULL;
    for (sub_sep = bssp->seq_set; sub_sep != NULL && !IS_Bioseq (sub_sep); sub_sep = sub_sep->next)
    {
      prev_sep = sub_sep;
    }
    if (sub_sep != NULL)
    {
      bsp = sub_sep->data.ptrvalue;
      new_bsp = GetDeltaSeqFromMasterSeg (sub_sep->data.ptrvalue);
      new_sep = SeqEntryNew();
      new_sep->choice = 1;
      new_sep->data.ptrvalue = new_bsp;
            
      /* add new seq entry to parent set */
      AddSeqEntryToSeqEntry (parent_set->seqentry, new_sep, TRUE);

      /* remove segset */      
      bssp->idx.deleteme = TRUE;
    }
    SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
    DeleteMarkedObjects (0, OBJ_BIOSEQSET, parent_set);
    SeqMgrReplaceInBioseqIndex (new_bsp); 
  }
  else
  {
    for (sub_sep = bssp->seq_set; sub_sep != NULL; sub_sep = next_sep)
    {
      next_sep = sub_sep->next;
      ConvertSegSetsToDeltaSequences (sub_sep);
    }
  }
}

static PubMedFetchFunc pmf_pubfetch = NULL;

NLM_EXTERN void LIBCALL PubMedSetFetchFunc (PubMedFetchFunc func)

{
  pmf_pubfetch = func;
}

NLM_EXTERN PubmedEntryPtr LIBCALL GetPubMedForUid (Int4 uid)

{
  PubMedFetchFunc  func;

  if (uid < 1) return NULL;
  func = pmf_pubfetch;
  if (func == NULL) return NULL;
  return func (uid);
}

static Boolean IsTerminator (int c)
{
  if (c == '\n' || c == '\r') {
    return TRUE;
  } else {
    return FALSE;
  }
}

typedef struct bufferedread {
  CharPtr data;
  Int4    len;
  Int4    offset;
} BufferedReadData, PNTR BufferedReadPtr;

static BufferedReadPtr BufferedReadFree (BufferedReadPtr brp)
{
  if (brp == NULL) return NULL;
  if (brp->data != NULL) {
    MemFree (brp->data);
    brp->data = NULL;
  }
  brp->offset = 0;
  brp->len = 0;
  return NULL;
}

extern void FreeBufferedReadList (ValNodePtr vnp)
{
  if (vnp == NULL) return;
  FreeBufferedReadList (vnp->next);
  vnp->next = NULL;
  vnp->data.ptrvalue = BufferedReadFree ( (BufferedReadPtr)vnp->data.ptrvalue); 
  ValNodeFree (vnp);
}

/* three possible return codes:
 * 0 = no terminators seen at all
 * 1 = have terminator plus one character
 * 2 = last is terminator - need more characters
 */
static Int4 HasTerminator (ValNodePtr list, Int4 PNTR len)
{
  CharPtr      cp;
  ValNodePtr   vnp;
  BufferedReadPtr brp;

  if (len == NULL) return 0;
  *len = 0;
  if (list == NULL) return 0;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    if (vnp->data.ptrvalue == NULL) continue;
    brp = (BufferedReadPtr) vnp->data.ptrvalue;
    if (brp->data == NULL) continue;
    for (cp = brp->data + brp->offset; *cp != 0; cp++) {
      if (IsTerminator (*cp)) {
        if (* (cp + 1) != 0 || vnp->next != NULL) {
          return 1;
        } else {
          return 2;
        }
      } else { 
        (*len) ++;
      }
    }
  }
  return 0;
}

static CharPtr GetLineFromBuffer (ValNodePtr PNTR current_data, Int4 len)
{
  ValNodePtr      vnp, next_vnp;
  BufferedReadPtr brp;
  CharPtr         cp;
  CharPtr         new_line;
  Int4            ctr;
  Char            this_terminator;
  CharPtr         next_char;

  if (current_data == NULL || *current_data == NULL) return NULL;

  new_line = MemNew (len + 1);
  if (new_line == NULL) return NULL;

  ctr = 0;
  vnp = *current_data;
  while (vnp != NULL && ctr < len) {
    if ((brp = (BufferedReadPtr)vnp->data.ptrvalue) == NULL || brp->data == NULL) {
      next_vnp = vnp->next;
      vnp->next = NULL;
      vnp->data.ptrvalue = BufferedReadFree (brp);
      ValNodeFree (vnp);
      vnp = next_vnp;
    } else {
      if (ctr + brp->len <= len) {
        MemCpy (new_line + ctr, brp->data + brp->offset, brp->len);
        ctr += brp->len;
        next_vnp = vnp->next;
        vnp->next = NULL;
        vnp->data.ptrvalue = BufferedReadFree (brp);
        ValNodeFree (vnp);
        vnp = next_vnp;
      } else {
        MemCpy (new_line + ctr, brp->data + brp->offset, len - ctr);
        brp->offset += len - ctr;
        brp->len -= (len - ctr);
        ctr = len;
      }
    }
  }
  if (vnp != NULL) {
    brp = (BufferedReadPtr)vnp->data.ptrvalue;
    if (brp->len >= 0) {
      cp = brp->data + brp->offset;
      this_terminator = *cp;
      /* handle condition when last character in data is terminator */
      if (* (cp + 1) == 0) {
        next_vnp = vnp->next;
        vnp->next = NULL;
        vnp->data.ptrvalue = BufferedReadFree (brp);
        ValNodeFree (vnp);
        vnp = next_vnp;
        while (vnp != NULL && (brp = (BufferedReadPtr)vnp->data.ptrvalue) == NULL) {
          next_vnp = vnp->next;
          vnp->next = NULL;
          vnp->data.ptrvalue = BufferedReadFree (brp);
          ValNodeFree (vnp);
          vnp = next_vnp;
        }
        if (vnp == NULL) {
          *current_data = NULL;
          new_line [len] = 0;
          return new_line;
        } else {
          next_char = brp->data + brp->offset;
          if (IsTerminator (*next_char) && *next_char != this_terminator) {
            brp->offset ++;
            brp->len --;
            if (brp->len == 0) {
              next_vnp = vnp->next;
              vnp->next = NULL;
              vnp->data.ptrvalue = BufferedReadFree (brp);
              ValNodeFree (vnp);
              vnp = next_vnp;
            }
          }
        }
      } else {
        next_char = cp + 1;
        if (IsTerminator (*next_char) && *next_char != this_terminator) {
          brp->offset += 2;
          brp->len -= 2;
        } else {
          brp->offset ++;
          brp->len --;
        }
      }
      if (brp->len <= 0) {
        next_vnp = vnp->next;
        vnp->next = NULL;
        vnp->data.ptrvalue = BufferedReadFree (brp);
        ValNodeFree (vnp);
        vnp = next_vnp;
      }
    }
  }
  *current_data = vnp;
  new_line [len] = 0;
  return new_line;
}

#define READ_BUFFER_SIZE 5000

static ValNodePtr AddToBuffer (ValNodePtr current_data, FILE *fp)
{
  ValNodePtr vnp;
  BufferedReadPtr brp;

  vnp = ValNodeNew (current_data);
  if (vnp == NULL) return NULL;
 
  brp = (BufferedReadPtr) MemNew (sizeof (BufferedReadData));
  if (brp == NULL) return NULL;
  brp->data = MemNew (READ_BUFFER_SIZE);
  if (brp->data == NULL) return NULL;
  brp->offset = 0;
 
  brp->len = fread (brp->data, 1, READ_BUFFER_SIZE - 1, fp);
  *(char *)(brp->data + brp->len) = 0; 

  vnp->data.ptrvalue = brp;
  return vnp;
}

extern CharPtr MyFGetLine (FILE *fp, ValNodePtr PNTR current_data)
{
  Int4       terminator_status;
  Int4       data_len;
  ValNodePtr last_vnp;

  terminator_status = HasTerminator (*current_data, &data_len);
  while (!feof (fp) && terminator_status == 0) {
    last_vnp = AddToBuffer (*current_data, fp);
    if (*current_data == NULL) {
      *current_data = last_vnp;
    }
    terminator_status = HasTerminator (*current_data, &data_len);
  }

  if (!feof (fp) && terminator_status == 2) {
    AddToBuffer (*current_data, fp);
  }
  return GetLineFromBuffer (current_data, data_len);
} 

/* PCR_primer manipulation functions */

static ValNodePtr ParsePCRComponent (
  CharPtr strs
)

{
  ValNodePtr  head = NULL;
  size_t      len;
  CharPtr     ptr, str, tmp;

  if (StringHasNoText (strs)) return NULL;

  tmp = StringSave (strs);
  if (tmp == NULL) return NULL;

  str = tmp;
  len = StringLen (str);
  if (len > 1 && *str == '(' && str [len - 1] == ')' && StringChr (str + 1, '(') == NULL) {
    str [len - 1] = '\0';
    str++;
  }

  while (StringDoesHaveText (str)) {
    ptr = StringChr (str, ',');
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

NLM_EXTERN ValNodePtr ParsePCRStrings (
  CharPtr fwd_primer_seq,
  CharPtr rev_primer_seq,
  CharPtr fwd_primer_name,
  CharPtr rev_primer_name
)

{
  ValNodePtr  curr_fwd_name;
  ValNodePtr  curr_fwd_seq;
  ValNodePtr  curr_rev_name;
  ValNodePtr  curr_rev_seq;
  CharPtr     fwd_name;
  CharPtr     fwd_seq;
  CharPtr     rev_name;
  CharPtr     rev_seq;
  ValNodePtr  fwd_name_list = NULL;
  ValNodePtr  fwd_seq_list = NULL;
  ValNodePtr  rev_name_list = NULL;
  ValNodePtr  rev_seq_list = NULL;
  ValNodePtr  head = NULL;
  Boolean     okay;
  Int2        orig_order = 0;
  PcrSetPtr   psp;

  fwd_seq_list = ParsePCRComponent (fwd_primer_seq);
  rev_seq_list = ParsePCRComponent (rev_primer_seq);
  fwd_name_list = ParsePCRComponent (fwd_primer_name);
  rev_name_list = ParsePCRComponent (rev_primer_name);
    
  curr_fwd_seq = fwd_seq_list;
  curr_rev_seq = rev_seq_list;
  curr_fwd_name = fwd_name_list;
  curr_rev_name = rev_name_list;

  while (curr_fwd_seq != NULL || curr_rev_seq != NULL || curr_fwd_name != NULL || curr_rev_name != NULL) {
    fwd_seq = NULL;
    rev_seq = NULL;
    fwd_name = NULL;
    rev_name = NULL;
    okay = FALSE;

    if (curr_fwd_seq != NULL) {
      fwd_seq = (CharPtr) curr_fwd_seq->data.ptrvalue;
      curr_fwd_seq = curr_fwd_seq->next;
      okay = TRUE;
    }

    if (curr_rev_seq != NULL) {
      rev_seq = (CharPtr) curr_rev_seq->data.ptrvalue;
      curr_rev_seq = curr_rev_seq->next;
      okay = TRUE;
    }

    if (curr_fwd_name != NULL) {
      fwd_name = (CharPtr) curr_fwd_name->data.ptrvalue;
      curr_fwd_name = curr_fwd_name->next;
      okay = TRUE;
    }

    if (curr_rev_name != NULL) {
      rev_name = (CharPtr) curr_rev_name->data.ptrvalue;
      curr_rev_name = curr_rev_name->next;
      okay = TRUE;
    }

    if (okay) {
      psp = (PcrSetPtr) MemNew (sizeof (PcrSet));
      if (psp != NULL) {
        psp->fwd_seq = StringSaveNoNull (fwd_seq);
        psp->rev_seq = StringSaveNoNull (rev_seq);
        psp->fwd_name = StringSaveNoNull (fwd_name);
        psp->rev_name = StringSaveNoNull (rev_name);
        orig_order++;
        psp->orig_order = orig_order;
        ValNodeAddPointer (&head, 0, (Pointer) psp);
      }
    }
  }

  ValNodeFreeData (fwd_seq_list);
  ValNodeFreeData (rev_seq_list);
  ValNodeFreeData (fwd_name_list);
  ValNodeFreeData (rev_name_list);

  return head;
}

NLM_EXTERN ValNodePtr ParsePCRSet (
  BioSourcePtr biop
)

{
  CharPtr       fwd_primer_seq = NULL;
  CharPtr       rev_primer_seq = NULL;
  CharPtr       fwd_primer_name = NULL;
  CharPtr       rev_primer_name = NULL;
  SubSourcePtr  ssp;

  if (biop == NULL) return NULL;

  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_fwd_primer_seq) {
      fwd_primer_seq = ssp->name;
    } else if (ssp->subtype == SUBSRC_rev_primer_seq) {
      rev_primer_seq = ssp->name;
    } else if (ssp->subtype == SUBSRC_fwd_primer_name) {
      fwd_primer_name = ssp->name;
    } else if (ssp->subtype == SUBSRC_rev_primer_name) {
      rev_primer_name = ssp->name;
    }
  }

  return ParsePCRStrings (fwd_primer_seq, rev_primer_seq, fwd_primer_name, rev_primer_name);
}

NLM_EXTERN int LIBCALLBACK SortVnpByPCRSetSeq (VoidPtr ptr1, VoidPtr ptr2)

{
  int         compare;
  PcrSetPtr   psp1, psp2;
  ValNodePtr  vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  psp1 = (PcrSetPtr) vnp1->data.ptrvalue;
  psp2 = (PcrSetPtr) vnp2->data.ptrvalue;
  if (psp1 == NULL || psp2 == NULL) return 0;

  compare = StringICmp (psp1->fwd_seq, psp2->fwd_seq);
  if (compare != 0) return compare;

  compare = StringICmp (psp1->rev_seq, psp2->rev_seq);
  if (compare != 0) return compare;

  compare = StringICmp (psp1->fwd_name, psp2->fwd_name);
  if (compare != 0) return compare;

  compare = StringICmp (psp1->rev_name, psp2->rev_name);
  if (compare != 0) return compare;

  if (psp1->orig_order > psp2->orig_order) {
    return 1;
  } else if (psp1->orig_order < psp2->orig_order) {
    return -1;
  }

  return 0;
}

NLM_EXTERN ValNodePtr UniqueVnpByPCRSetSeq (ValNodePtr pset)

{
  PcrSetPtr     last;
  ValNodePtr    next;
  Pointer PNTR  prev;
  PcrSetPtr     psp;
  ValNodePtr    vnp;

  if (pset == NULL) return NULL;
  last = (PcrSetPtr) pset->data.ptrvalue;
  vnp = pset->next;
  prev = (Pointer PNTR) &(pset->next);
  while (vnp != NULL) {
    next = vnp->next;
    psp = (PcrSetPtr) vnp->data.ptrvalue;
    if (last != NULL && psp != NULL &&
        StringICmp (last->fwd_seq, psp->fwd_seq) == 0 &&
        StringICmp (last->rev_seq, psp->rev_seq) == 0 &&
        StringICmp (last->fwd_name, psp->fwd_name) == 0 &&
        StringICmp (last->rev_name, psp->rev_name) == 0) {
      vnp->next = NULL;
      *prev = next;
      MemFree (psp->fwd_seq);
      MemFree (psp->rev_seq);
      MemFree (psp->fwd_name);
      MemFree (psp->rev_name);
      ValNodeFreeData (vnp);
    } else {
      last = (PcrSetPtr) vnp->data.ptrvalue;
      prev = (Pointer PNTR) &(vnp->next);
    }
    vnp = next;
  }

  return pset;
}

NLM_EXTERN int LIBCALLBACK SortVnpByPCRSetOrder (VoidPtr ptr1, VoidPtr ptr2)

{
  PcrSetPtr   psp1, psp2;
  ValNodePtr  vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  psp1 = (PcrSetPtr) vnp1->data.ptrvalue;
  psp2 = (PcrSetPtr) vnp2->data.ptrvalue;
  if (psp1 == NULL || psp2 == NULL) return 0;

  if (psp1->orig_order > psp2->orig_order) {
    return 1;
  } else if (psp1->orig_order < psp2->orig_order) {
    return -1;
  }

  return 0;
}

static CharPtr CombinePCRItems (
  ValNodePtr list
)

{
  Int4        count;
  size_t      len;
  CharPtr     ptr;
  CharPtr     str;
  ValNodePtr  vnp;

  if (list == NULL) return NULL;
  count = ValNodeLen (list);
  if (count == 1) {
    ptr = (CharPtr) list->data.ptrvalue;
    return StringSaveNoNull (ptr);
  }

  len = 0;
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    ptr = (CharPtr) vnp->data.ptrvalue;
    if (ptr == NULL) continue;
    len += StringLen (ptr) + 1;
  }
  str = (CharPtr) MemNew (sizeof (Char) * (len + 4));
  if (str == NULL) return NULL;
  StringCpy (str, "(");

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    ptr = (CharPtr) vnp->data.ptrvalue;
    if (ptr == NULL) continue;
    StringCat (str, ptr);
    if (vnp->next != NULL) {
      StringCat (str, ",");
    }
  }

  StringCat (str, ")");
  return str;
}

NLM_EXTERN SubSourcePtr WritePCRSet (
  ValNodePtr pset
)

{
  ValNodePtr    fwd_name_list = NULL;
  ValNodePtr    fwd_seq_list = NULL;
  ValNodePtr    rev_name_list = NULL;
  ValNodePtr    rev_seq_list = NULL;
  SubSourcePtr  head = NULL;
  SubSourcePtr  last = NULL;
  PcrSetPtr     psp;
  SubSourcePtr  ssp;
  CharPtr       str;
  ValNodePtr    vnp;

  if (pset == NULL) return NULL;

  for (vnp = pset; vnp != NULL; vnp = vnp->next) {
    psp = (PcrSetPtr) vnp->data.ptrvalue;
    if (psp == NULL) continue;
    if (StringDoesHaveText (psp->fwd_seq)) {
      ValNodeCopyStr (&fwd_seq_list, 0, psp->fwd_seq);
    }
    if (StringDoesHaveText (psp->rev_seq)) {
      ValNodeCopyStr (&rev_seq_list, 0, psp->rev_seq);
    }
    if (StringDoesHaveText (psp->fwd_name)) {
      ValNodeCopyStr (&fwd_name_list, 0, psp->fwd_name);
    }
    if (StringDoesHaveText (psp->rev_name)) {
      ValNodeCopyStr (&rev_name_list, 0, psp->rev_name);
    }
  }

  str = CombinePCRItems (fwd_seq_list);
  if (str != NULL) {
    ssp = SubSourceNew ();
    ssp->subtype = SUBSRC_fwd_primer_seq;
    ssp->name = str;
    if (head == NULL) {
      head = ssp;
    }
    if (last != NULL) {
      last->next = ssp;
    }
    last = ssp;
  }

  str = CombinePCRItems (rev_seq_list);
  if (str != NULL) {
    ssp = SubSourceNew ();
    ssp->subtype = SUBSRC_rev_primer_seq;
    ssp->name = str;
    if (head == NULL) {
      head = ssp;
    }
    if (last != NULL) {
      last->next = ssp;
    }
    last = ssp;
  }

  str = CombinePCRItems (fwd_name_list);
  if (str != NULL) {
    ssp = SubSourceNew ();
    ssp->subtype = SUBSRC_fwd_primer_name;
    ssp->name = str;
    if (head == NULL) {
      head = ssp;
    }
    if (last != NULL) {
      last->next = ssp;
    }
    last = ssp;
  }

  str = CombinePCRItems (rev_name_list);
  if (str != NULL) {
    ssp = SubSourceNew ();
    ssp->subtype = SUBSRC_rev_primer_name;
    ssp->name = str;
    if (head == NULL) {
      head = ssp;
    }
    if (last != NULL) {
      last->next = ssp;
    }
    last = ssp;
  }

  return head;
}

NLM_EXTERN ValNodePtr FreePCRSet (
  ValNodePtr pset
)

{
  PcrSetPtr   psp;
  ValNodePtr  vnp;

  if (pset == NULL) return NULL;

  for (vnp = pset; vnp != NULL; vnp = vnp->next) {
    psp = (PcrSetPtr) vnp->data.ptrvalue;
    if (psp == NULL) continue;
    MemFree (psp->fwd_seq);
    MemFree (psp->rev_seq);
    MemFree (psp->fwd_name);
    MemFree (psp->rev_name);
  }

  return ValNodeFreeData (pset);
}

static ValNodePtr ParsePCRColonString (
  CharPtr strs
)

{
  ValNodePtr  head = NULL;
  size_t      len;
  CharPtr     ptr, str, tmp;

  if (StringHasNoText (strs)) return NULL;

  tmp = StringSave (strs);
  str = tmp;
  len = StringLen (str);
  if (len > 1 && StringChr (str, ':') != NULL) {
    while (StringDoesHaveText (str)) {
      ptr = StringChr (str, ':');
      if (ptr != NULL) {
        *ptr = '\0';
        ptr++;
      }
      TrimSpacesAroundString (str);
      ValNodeCopyStr (&head, 0, str);
      str = ptr;
    }
  } else {
    ValNodeCopyStr (&head, 0, str);
  }

  MemFree (tmp);
  return head;
}

static CharPtr FusePrimerNames (
  CharPtr first,
  CharPtr second
)

{
  size_t   len;
  CharPtr  str;

  if (first == NULL) return second;
  if (second == NULL) return first;

  len = StringLen (first) + StringLen (second) + 5;
  str = MemNew (len);
  if (str == NULL) return NULL;

  StringCpy (str, first);
  StringCat (str, ":");
  StringCat (str, second);

  return str;
}

static PCRPrimerPtr ModernizePCRPrimerHalf (
  CharPtr seq,
  CharPtr name
)

{
  CharPtr       curr_name = NULL, curr_seq = NULL, fused_name;
  PCRPrimerPtr  curr_primer = NULL, last_primer = NULL, primer_set = NULL;
  ValNodePtr    name_list, seq_list, name_vnp, seq_vnp;

  seq_list = ParsePCRColonString (seq);
  name_list = ParsePCRColonString (name);

  seq_vnp = seq_list;
  name_vnp = name_list;

  while (seq_vnp != NULL /* || name_vnp != NULL */) {
    if (seq_vnp != NULL) {
      curr_seq = (CharPtr) seq_vnp->data.ptrvalue;
      seq_vnp = seq_vnp->next;
    }
    if (name_vnp != NULL) {
      curr_name = (CharPtr) name_vnp->data.ptrvalue;
      name_vnp = name_vnp->next;
    } else {
      curr_name = NULL;
    }

    curr_primer = (PCRPrimerPtr) MemNew (sizeof (PCRPrimer));
    if (curr_primer != NULL) {
      curr_primer->seq = StringSaveNoNull (curr_seq);
      curr_primer->name = StringSaveNoNull (curr_name);

      if (primer_set == NULL) {
        primer_set = curr_primer;
      }
      if (last_primer != NULL) {
        last_primer->next = curr_primer;
      }
      last_primer = curr_primer;
    }
  }

  while (name_vnp != NULL && last_primer != NULL) {
    curr_name = (CharPtr) name_vnp->data.ptrvalue;
    fused_name = FusePrimerNames (last_primer->name, curr_name);
    MemFree (last_primer->name);
    last_primer->name = StringSaveNoNull (fused_name);
    name_vnp = name_vnp->next;
  }

  while (name_vnp != NULL && last_primer == NULL) {
    curr_name = (CharPtr) name_vnp->data.ptrvalue;
    curr_primer = (PCRPrimerPtr) MemNew (sizeof (PCRPrimer));
    if (curr_primer != NULL) {
      curr_primer->name = StringSaveNoNull (curr_name);

      if (primer_set == NULL) {
        primer_set = curr_primer;
      }
      if (last_primer != NULL) {
        last_primer->next = curr_primer;
      }
      last_primer = curr_primer;
    }
    name_vnp = name_vnp->next;
  }

  ValNodeFreeData (seq_list);
  ValNodeFreeData (name_list);

  return primer_set;
}

NLM_EXTERN void ModernizePCRPrimers (
  BioSourcePtr biop
)

{
  PCRReactionSetPtr  curr_reaction, last_reaction = NULL, reaction_set = NULL;
  PCRPrimerPtr       forward, reverse;
  PcrSetPtr          psp;
  ValNodePtr         pset, vnp;
  SubSourcePtr       nextssp;
  SubSourcePtr PNTR  prevssp;
  SubSourcePtr       ssp;
  Boolean            unlink;

  if (biop == NULL) return;
  /* if (biop->pcr_primers != NULL) return; */

  pset = ParsePCRSet (biop);
  if (pset == NULL) return;

  for (vnp = pset; vnp != NULL; vnp = vnp->next) {
    psp = (PcrSetPtr) vnp->data.ptrvalue;
    if (psp == NULL) continue;

    forward = ModernizePCRPrimerHalf (psp->fwd_seq, psp->fwd_name);
    reverse = ModernizePCRPrimerHalf (psp->rev_seq, psp->rev_name);

    if (forward != NULL || reverse != NULL) {

      curr_reaction = (PCRReactionSetPtr) MemNew (sizeof (PCRReactionSet));
      if (curr_reaction != NULL) {
        curr_reaction->forward = forward;
        curr_reaction->reverse = reverse;

        if (reaction_set == NULL) {
          reaction_set = curr_reaction;
        }
        if (last_reaction != NULL) {
          last_reaction->next = curr_reaction;
        }
        last_reaction = curr_reaction;
      }
    }
  }

  FreePCRSet (pset);

  if (reaction_set != NULL) {
    if (last_reaction != NULL) {
      /* merge with existing structured pcr_primers */
      last_reaction->next = biop->pcr_primers;
    }
    biop->pcr_primers = reaction_set;

    ssp = biop->subtype;
    prevssp = (SubSourcePtr PNTR) &(biop->subtype);
    while (ssp != NULL) {
      nextssp = ssp->next;
      unlink= FALSE;

      if (ssp->subtype == SUBSRC_fwd_primer_seq ||
          ssp->subtype == SUBSRC_rev_primer_seq ||
          ssp->subtype == SUBSRC_fwd_primer_name ||
          ssp->subtype == SUBSRC_rev_primer_name) {
        unlink = TRUE;
      }

      if (unlink) {
        *prevssp = ssp->next;
        ssp->next = NULL;
        SubSourceFree (ssp);
      } else {
        prevssp = (SubSourcePtr PNTR) &(ssp->next);
      }
      ssp = nextssp;
    }
  }
}


NLM_EXTERN void ModernizeRNAFields (
  SeqFeatPtr sfp
)

{
  GBQualPtr       gbq;
  GBQualPtr       nextqual;
  GBQualPtr PNTR  prevqual;
  RNAGenPtr       rgp;
  RNAQualPtr      rqp;
  RNAQualPtr      last_rqp = NULL;
  RnaRefPtr       rrp;
  CharPtr         str;
  Boolean         unlink;

  if (sfp == NULL) return;
  if (sfp->data.choice != SEQFEAT_RNA) return;

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL) return;

  if (rrp->type != 255 || rrp->ext.choice != 1) return;
  str = rrp->ext.value.ptrvalue;
  if (StringHasNoText (str)) return;

  if (StringCmp (str, "ncRNA") == 0) {
    rrp->type = 8;
  } else if (StringCmp (str, "tmRNA") == 0) {
    rrp->type = 9;
  } else if (StringCmp (str, "misc_RNA") == 0) {
    rrp->type = 10;
  } else return;
  rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
  rrp->ext.choice = 0;

  rgp = (RNAGenPtr) MemNew (sizeof (RNAGen));
  if (rgp == NULL) return;
  rrp->ext.choice = 3;
  rrp->ext.value.ptrvalue = (Pointer) rgp;

  gbq = sfp->qual;
  prevqual = (GBQualPtr PNTR) &(sfp->qual);
  while (gbq != NULL) {
    nextqual = gbq->next;
    unlink = FALSE;
    if (StringCmp (gbq->qual, "ncRNA_class") == 0) {
      rgp->_class = StringSaveNoNull (gbq->val);
      unlink = TRUE;
    } else if (StringCmp (gbq->qual, "product") == 0) {
      rgp->product = StringSaveNoNull (gbq->val);
      unlink = TRUE;
    } else if (StringCmp (gbq->qual, "tag_peptide") == 0) {
      rqp = (RNAQualPtr) MemNew (sizeof (RNAQual));
      if (rqp != NULL) {
        rqp->qual = StringSave (gbq->qual);
        rqp->val = StringSave (gbq->val);
        if (rgp->quals == NULL) {
          rgp->quals = rqp;
        }
        if (last_rqp != NULL) {
          last_rqp->next = rqp;
        }
        last_rqp = rqp;
        unlink = TRUE;
      }
    }
    if (unlink) {
      *(prevqual) = gbq->next;
      gbq->next = NULL;
      GBQualFree (gbq);
    } else {
      prevqual = (GBQualPtr PNTR) &(gbq->next);
    }
    gbq = nextqual;
  }
}

static DbtagPtr DbtagParse (
  CharPtr str
)

{
  Boolean      all_digits = TRUE;
  Char         ch;
  DbtagPtr     dbt;
  long         num;
  Int2         num_digits = 0;
  ObjectIdPtr  oip;
  CharPtr      ptr;
  CharPtr      tmp;

  if (StringHasNoText (str)) return NULL;
  ptr = StringChr (str, ':');
  if (ptr == NULL) return NULL;

  dbt = DbtagNew ();
  oip = ObjectIdNew ();
  if (dbt == NULL || oip == NULL) return NULL;

  if (ptr != NULL) {
    *ptr = '\0';
    ptr++;
  }

  dbt->db = StringSave (str);
  dbt->tag = oip;

  tmp = ptr;
  ch = *tmp;
  while (ch != '\0') {
    if (IS_DIGIT (ch)) {
      num_digits++;
    } else {
      all_digits = FALSE;
    }
    tmp++;
    ch = *tmp;
  }

  if (all_digits) {
    if (num_digits < 10 || (num_digits == 10 && StringCmp (ptr, "2147483647") <= 0)) {
      sscanf (ptr, "%ld", &num);
      oip->id = (Int4) num;
      return dbt;
    }
  }

  oip->str = StringSave (ptr);

  return dbt;
}

static void GetNomenclatureUOP (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr         oip;
  UserObjectPtr PNTR  uopp;

  if (uop == NULL || userdata == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "OfficialNomenclature") != 0) return;
  uopp = (UserObjectPtr PNTR) userdata;
  *uopp = uop;
}

NLM_EXTERN void ModernizeGeneFields (
  SeqFeatPtr sfp
)

{
  GeneNomenclaturePtr  gnp;
  GeneRefPtr           grp;
  ObjectIdPtr          oip;
  CharPtr              str;
  CharPtr              symbol = NULL, name = NULL, source = NULL;
  Uint2                status = 0;
  UserFieldPtr         ufp;
  UserObjectPtr        uop = NULL;
  UserObjectPtr        curr, next;
  UserObjectPtr PNTR   prev;

  if (sfp == NULL) return;
  if (sfp->data.choice != SEQFEAT_GENE) return;

  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return;

  if (grp->formal_name != NULL) return;

  if (sfp->ext == NULL) return;
  VisitUserObjectsInUop (sfp->ext, (Pointer) &uop, GetNomenclatureUOP);
  if (uop == NULL) return;

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (oip == NULL || oip->str == NULL) continue;
    if (StringICmp (oip->str, "Symbol") == 0) {
      if (ufp->choice == 1) {
        str = (CharPtr) ufp->data.ptrvalue;
        if (str != NULL) {
          symbol = str;
        }
      }
    } else if (StringICmp (oip->str, "Name") == 0) {
      if (ufp->choice == 1) {
        str = (CharPtr) ufp->data.ptrvalue;
        if (str != NULL) {
          name = str;
        }
      }
    } else if (StringICmp (oip->str, "DataSource") == 0) {
      if (ufp->choice == 1) {
        str = (CharPtr) ufp->data.ptrvalue;
        if (str != NULL) {
          source = str;
        }
      }
    } else if (StringICmp (oip->str, "Status") == 0) {
      if (ufp->choice == 1) {
        str = (CharPtr) ufp->data.ptrvalue;
        if (str != NULL) {
          if (StringICmp (str, "Official") == 0) {
            status = 1;
          } else if (StringICmp (str, "Interim") == 0) {
            status = 2;
          }
        }
      }
    }
  }
  if (symbol == NULL && name == NULL && source == NULL && status == 0) return;

  gnp = GeneNomenclatureNew ();
  if (gnp == NULL) return;

  gnp->status = status;
  gnp->symbol = StringSaveNoNull (symbol);
  gnp->name = StringSaveNoNull (name);
  gnp->source = DbtagParse (source);

  grp->formal_name = gnp;

  prev = (UserObjectPtr PNTR) &(sfp->ext);
  curr = sfp->ext;
  while (curr != NULL) {
    next = curr->next;
    if (uop == curr) {
      *(prev) = curr->next;
      curr->next = NULL;
      UserObjectFree (curr);
    } else {
      prev = (UserObjectPtr PNTR) &(curr->next);
    }
    curr = next;
  }
}



static void AddDefLinesToAlignmentSequences 
(TAlignmentFilePtr afp,
 SeqEntryPtr sep_head)
{
  BioseqSetPtr bssp;
  SeqEntryPtr  sep;
  Int4         index;
  ValNodePtr   sdp;
  CharPtr      new_title;
  Uint4         new_title_len;
  Int4         curr_seg;
  Int4         num_sets = 1;
  Boolean      one_defline_per_sequence = TRUE;
  Boolean      all_extra_empty;

 
  if (afp == NULL || sep_head == NULL || ! IS_Bioseq_set (sep_head))
  {
    return;
  }
  if ((afp->num_deflines == 0 || afp->deflines == NULL)
    && (afp->num_organisms == 0 || afp->organisms == NULL))
  {
    return;
  }
  bssp = sep_head->data.ptrvalue;
  
  /* find out if all of our deflines are real */
  if (afp->num_segments > 1 && afp->num_deflines == afp->num_sequences)
  {
    one_defline_per_sequence = FALSE;
    num_sets = afp->num_sequences / afp->num_segments;
    all_extra_empty = TRUE;
    for (curr_seg = num_sets; curr_seg < afp->num_deflines && all_extra_empty; curr_seg ++)
    {
      if (afp->deflines [curr_seg] != NULL)
      {
        all_extra_empty = FALSE;
      }
    }
      if (all_extra_empty)
      {
          one_defline_per_sequence = TRUE;
      }
  }
  
  for (sep = bssp->seq_set, index = 0;
       sep != NULL && (index < afp->num_deflines || index < afp->num_organisms);
       sep = sep->next, index++)
  {
    new_title_len = 0;
    /* get lengths for organisms for this sequence */
    
    if (afp->num_segments > 1 && afp->num_organisms == afp->num_sequences)
    {
      /* have one organism per segment, in which case use only the first one */
      curr_seg = index * afp->num_segments;
    }
    else
    { /* otherwise one organism per sequence */
      curr_seg = index;
    }
    if (curr_seg < afp->num_organisms)
    {
      new_title_len += StringLen (afp->organisms [curr_seg]) + 1;
    }

    /* get lengths for deflines for this sequence */
    if (! one_defline_per_sequence)
    { /* have one defline per segment, in which use only the first one */
      curr_seg = index * afp->num_segments;
    }
    else
    { /* otherwise one defline per sequence */
      curr_seg = index;    
    }
    if (curr_seg < afp->num_deflines) 
    {
      new_title_len += StringLen (afp->deflines [curr_seg]) + 1;
    }
      
    if (new_title_len > 0) {
      new_title = (CharPtr) MemNew (new_title_len);
      if (new_title == NULL) return;
      new_title [0] = 0;
      
      /* list organisms at beginning of new defline */
      if (afp->num_segments > 1 && afp->num_organisms == afp->num_sequences)
      { /* have one organism per segment, in which case use only first one */
        curr_seg = index * afp->num_segments;
      }
      else
      { /* otherwise one organism per sequence */
          curr_seg = index;
      }

      if (curr_seg < afp->num_organisms) {
        StringCat (new_title, afp->organisms [curr_seg]);
        if (new_title_len > StringLen (new_title) + 1)
        {
          StringCat (new_title, " ");
        }
      }
      
      if (!one_defline_per_sequence)
      { /* have one defline per segment, in which case all go to same sequence */
        curr_seg = index * afp->num_segments;
      }
      else 
      {
          curr_seg = index;
      }
      if (curr_seg < afp->num_deflines) 
      {
        StringCat (new_title, afp->deflines [curr_seg]);
      }

      sdp = CreateNewDescriptor (sep, Seq_descr_title);
      if (sdp != NULL) {
        sdp->data.ptrvalue = new_title;
      } else {
        MemFree (new_title);
      }
    }
  }
}

#if 0
static SeqEntryPtr 
MakeDeltaSetFromAlignment 
(SeqEntryPtr sep_list,
 TAlignmentFilePtr afp,
 Uint1 moltype,
 Int4  gap_length
 )
{
  BioseqPtr    bsp, deltabsp;
  SeqEntryPtr  this_list, last_sep, next_list, sep, nextsep;
  SeqEntryPtr  topsep, last_delta_sep;
  SeqIdPtr     sip;
  Int4         curr_seg;
  CharPtr      seqbuf;
  ValNodePtr   vnp;
  SeqLitPtr    slp;
  IntFuzzPtr   ifp;
  SeqEntryPtr  delta_list = NULL;
  
  delta_list = NULL;
  last_delta_sep = NULL;
  this_list = sep_list;
  while (this_list != NULL)
  {
    last_sep = this_list;
    curr_seg = 0;
    while (last_sep != NULL && curr_seg < afp->num_segments - 1)
    {
        last_sep = last_sep->next;
        curr_seg++;
    }
    if (last_sep == NULL) return NULL;
    next_list = last_sep->next;
    last_sep->next = NULL;
  
    bsp = (BioseqPtr)this_list->data.ptrvalue;
    if (bsp == NULL) return NULL;

    sip = SeqIdDup (bsp->id);
    vnp = ValNodeExtract (&(bsp->descr), Seq_descr_title);

    deltabsp = BioseqNew ();
    if (deltabsp == NULL) return NULL;
    deltabsp->repr = Seq_repr_delta;
    deltabsp->seq_ext_type = 4;
    deltabsp->mol = moltype;
    deltabsp->length = 0;

    topsep = SeqEntryNew ();
    if (topsep == NULL) return NULL;
    topsep->choice = 1;
    topsep->data.ptrvalue = (Pointer) deltabsp;

    for (sep = this_list; sep != NULL; sep = nextsep) {
      nextsep = sep->next;
      sep->next = NULL;

      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (bsp == NULL) continue;

      if (bsp->repr == Seq_repr_raw) {
        BioseqRawConvert (bsp, Seq_code_iupacna);
        seqbuf = BSMerge ((ByteStorePtr) bsp->seq_data, NULL);
        slp = (SeqLitPtr) MemNew (sizeof (SeqLit));
        if (slp == NULL) continue;

        slp->length = bsp->length;
        ValNodeAddPointer ((ValNodePtr PNTR) &(deltabsp->seq_ext), (Int2) 2, (Pointer) slp);
        slp->seq_data = BSNew (slp->length);
        slp->seq_data_type = Seq_code_iupacna;
        AddBasesToByteStore (slp->seq_data, seqbuf);
        MemFree(seqbuf);

        deltabsp->length += slp->length;

      } else if (bsp->repr == Seq_repr_virtual) {
        slp = (SeqLitPtr) MemNew (sizeof (SeqLit));
        if (slp == NULL) continue;
        slp->length = bsp->length;
        if (slp == NULL) continue;

        slp->length = bsp->length;
        ValNodeAddPointer ((ValNodePtr PNTR) &(deltabsp->seq_ext), (Int2) 2, (Pointer) slp);
        if (slp->length < 1) {
          slp->length = 0;
          ifp = IntFuzzNew ();
          ifp->choice = 4;
          slp->fuzz = ifp;
        }

        deltabsp->length += slp->length;
      }
      SeqEntryFree (sep);
      
      if (nextsep != NULL)
      {
        /* add gap */
        slp = (SeqLitPtr) MemNew (sizeof (SeqLit));
        if (slp == NULL) continue;
        slp->length = gap_length;
        ValNodeAddPointer ((ValNodePtr PNTR) &(deltabsp->seq_ext), (Int2) 2, (Pointer) slp);
        deltabsp->length += slp->length;        
      }
    }

    ValNodeLink (&(deltabsp->descr), vnp);
    deltabsp->id = sip;
    
    if (last_delta_sep == NULL)
    {
        delta_list = topsep;
    }
    else 
    {
        last_delta_sep->next = topsep;    
    }
    last_delta_sep = topsep;
    
    this_list = next_list;    
  }
  return delta_list;
}
#endif

static void RenameSegSet (SeqEntryPtr sep)
{
  BioseqSetPtr bssp, seg_bssp;
  SeqEntryPtr  seg_sep;
  BioseqPtr    main_bsp = NULL;
  BioseqPtr    seg_bsp = NULL;
  Char         new_id_str [255];
  
  if (sep == NULL || !IS_Bioseq_set (sep) || (bssp = sep->data.ptrvalue) == NULL
      || bssp->_class != BioseqseqSet_class_segset)
  {
      return;
  }
  
  sep = bssp->seq_set;
  while (sep != NULL && (seg_bsp == NULL || main_bsp == NULL))
  {
      if (IS_Bioseq (sep))
      {
        main_bsp = (BioseqPtr) sep->data.ptrvalue;
      }
      else if (IS_Bioseq_set (sep))
      {
        seg_bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (seg_bssp != NULL && seg_bssp->_class == BioseqseqSet_class_parts)
        {
            seg_sep = seg_bssp->seq_set;
            while (seg_sep != NULL && seg_bsp == NULL)
            {
              if (IS_Bioseq (seg_sep))
              {
                seg_bsp = seg_sep->data.ptrvalue;
            }
            seg_sep = seg_sep->next;
          }
        }
      }
      sep = sep->next;
  }
  if (main_bsp == NULL || seg_bsp == NULL)
  {
      return;
  }
  SeqIdWrite (seg_bsp->id, new_id_str, PRINTID_FASTA_SHORT, sizeof (new_id_str) - 7);
  StringCat (new_id_str, "_master");
  SeqIdFree (main_bsp->id);
  main_bsp->id = MakeSeqID (new_id_str);
}

static SeqEntryPtr 
MakeSegmentedSetFromAlignment 
(SeqEntryPtr       sep_list,
 TAlignmentFilePtr afp,
 Uint1             moltype,
 Int4Ptr           segs_per_set)
{
  SeqEntryPtr  this_list, last_sep, next_list, nextsep, last_segset;
  Int4         curr_seg;
  Int4         set_index = 0;
  
  this_list = sep_list;
  sep_list = NULL;
  last_segset = NULL;
  while (this_list != NULL)
  {
    last_sep = this_list;
    curr_seg = 0;
    while (last_sep != NULL && curr_seg < segs_per_set [set_index] - 1)
    {
      if (!IS_Bioseq (last_sep)) return NULL;
      last_sep = last_sep->next;
      curr_seg++;
    }
    if (last_sep == NULL) return NULL;
    next_list = last_sep->next;
    last_sep->next = NULL;
    
    last_sep = this_list->next;
    this_list->next = NULL;
    while (last_sep != NULL)
    {
      nextsep = last_sep->next;
      last_sep->next = NULL;
      AddSeqEntryToSeqEntry (this_list, last_sep, FALSE);
      last_sep = nextsep;
    }

    /* fix IDs for seg sets */    
    RenameSegSet (this_list);
    
    if (sep_list == NULL) 
    {
      sep_list = this_list;
    }
    else
    {
      last_segset->next = this_list;
    }
    last_segset = this_list;
    
    this_list = next_list;
    set_index++;
  }
  return sep_list;     
}


extern CharPtr AlignmentStringToSequenceString (CharPtr aln_str, Uint1 moltype)
{
  CharPtr cp_aln, cp_seq;
  Char    ch;
  CharPtr seq_str;
  
  if (aln_str == NULL) return NULL;
  seq_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (aln_str) + 1));
  if (seq_str == NULL) return NULL;
  cp_seq = seq_str;
  for (cp_aln = aln_str; *cp_aln != 0; cp_aln++)
  {
    ch = *cp_aln;
    ch = TO_UPPER (ch); 
    if ( ISA_na (moltype) ) 
    {
      if (ch == 'U') ch = 'T';
      if (ch == 'X') ch = 'N';
      if ( StringChr ("EFIJLOPQXZ-.*", ch) == NULL )  
      { 
        *cp_seq = ch;
        cp_seq++;
      }
    }
    else 
    {
      if ( StringChr("JO-.", ch) == NULL ) 
      {
        *cp_seq = ch;
        cp_seq++;
      } 
    }
  }
  *cp_seq = 0;
  return seq_str;
}

NLM_EXTERN SeqEntryPtr SequenceStringToSeqEntry (CharPtr str, SeqIdPtr sip, Uint1 mol_type)
{
  SeqEntryPtr  sep;
  BioseqPtr    bsp;
  ByteStorePtr bs;

  if (str == NULL || sip == NULL) return NULL;
  sep = SeqEntryNew ();
  if (sep == NULL) return NULL;
  bsp = BioseqNew ();
  if (bsp == NULL) 
  { 
    ValNodeFree (sep); 
    return NULL; 
  }
  sep->choice = 1;
  sep->data.ptrvalue = (Pointer) bsp;
  bsp->id = SeqIdDup (sip);
  bsp->id->next = NULL;
  /* Note - use SeqMgrReplaceInBioseqIndex instead of SeqMgrAddToBioseqIndex because 
   * BioseqNew called SeqMgrAddToBioseqIndex, without the IDs, and won't add it again.
   */
  SeqMgrReplaceInBioseqIndex (bsp);
  bsp->repr = Seq_repr_raw;
  if ( ISA_na (mol_type) ) 
  {
    bsp->mol = Seq_mol_na;
    bsp->seq_data_type = Seq_code_iupacna;
  } 
  else
  {
    bsp->mol = Seq_mol_aa;
    bsp->seq_data_type = Seq_code_ncbieaa;
  }
  bsp->length = StringLen (str);
  if ( bsp->length == 0 ) 
  {
    BioseqFree (bsp);
    ValNodeFree (sep);
    return NULL;
  }
  bs = BSNew (bsp->length);
  bsp->seq_data = (SeqDataPtr) bs;
  BSWrite (bs, str, bsp->length);

  return sep;
}

#if 0
static SeqEntryPtr MakeDeltaSeqsFromAlignmentSequences (TAlignmentFilePtr afp, Uint1 moltype, CharPtr PNTR seq_str)
{
  Int4            num_sets, next_start, k, index;
  SeqIdPtr        sip;
  SeqLitPtr       slip;
  SeqEntryPtr     sep_list = NULL, sep, sep_last = NULL;
  BioseqPtr       new_bsp;
  ValNodePtr      seq_ext = NULL;
  
  if (afp == NULL || seq_str == NULL) return NULL;
  
  num_sets = afp->num_sequences / afp->num_segments;
  for (k = 0; k < num_sets; k++)
  {
    sep = SeqEntryNew ();
    if (sep == NULL) return NULL;
    new_bsp = BioseqNew ();
    if (new_bsp == NULL) return NULL;
    sip = MakeSeqID (afp->ids [k * afp->num_segments]);
    new_bsp->id = sip;
    sep->choice = 1;
    sep->data.ptrvalue = new_bsp;
    SeqMgrAddToBioseqIndex (new_bsp);
    
    if (sep_last == NULL)
    {
      sep_list = sep;
    }
    else
    {
      sep_last->next = sep;
    }
    sep_last = sep;
    
    new_bsp->seq_data = NULL;
    new_bsp->seq_data_type = 0;
    new_bsp->repr = Seq_repr_delta;
    new_bsp->seq_ext_type = 4;
    new_bsp->mol = moltype;
    new_bsp->seq_ext = NULL;
    new_bsp->length = 0;
    next_start = (k + 1) * afp->num_segments;
    seq_ext = NULL;
    for (index = k * afp->num_segments; index < next_start; index++)
    {
      if (seq_ext != NULL)
      {
        /* insert gap of unknown length between the previous segment
         * and this one.
         */
        new_bsp->length += AddGapSeqLit (&seq_ext);
      }

      if (StringHasNoText (seq_str [index]))
      {
        /* add gap to represent missing sequence */
        new_bsp->length += AddGapSeqLit (&seq_ext);        
      }
      else
      {
        slip = (SeqLitPtr) MemNew (sizeof (SeqLit));
        if (slip != NULL) 
        {
          slip->length = StringLen (seq_str [index]);
          ValNodeAddPointer (&seq_ext, (Int2) 2, (Pointer) slip);
          slip->seq_data = BSNew (slip->length);
          slip->seq_data_type = Seq_code_iupacna;
          AddBasesToByteStore (slip->seq_data, seq_str [index]);
          new_bsp->length += slip->length;
        } 
      }
    }
    new_bsp->seq_ext = seq_ext;
    BioseqPack (new_bsp);        
  }
    
  return sep_list;
}
#endif

static SeqIdPtr GetFarPointerID (CharPtr id_str)
{
  CharPtr  tmp_id_str;
  CharPtr  cp_start, cp_end;
  Int4     len;
  SeqIdPtr sip;
  
  if (id_str == NULL)
  {
    return NULL;
  }

  cp_start = StringChr (id_str, '|');
  if (cp_start == NULL)
  {
    cp_start = id_str;
    len = StringLen (id_str);
  }
  else
  {
    cp_start++;
    cp_end = StringChr (cp_start, '|');
    if (cp_end == NULL)
    {
      len = StringLen (cp_start);
    }
    else
    {
      len = cp_end - cp_start;
    }
  }
  if (len == 0)
  {
    return NULL;
  }
  tmp_id_str = (CharPtr) MemNew ((len + 4) * sizeof (Char));
  if (tmp_id_str == NULL)
  {
    return NULL;
  }
  StringCpy (tmp_id_str, "acc");
  StringNCat (tmp_id_str, cp_start, len);
  tmp_id_str [len + 3] = 0;
  sip = MakeSeqID (tmp_id_str);
  MemFree (tmp_id_str);
  return sip;
}

static void ReplacePipesWithUnderscores (CharPtr seqid_str)
{
  CharPtr cp;
  
  if (seqid_str == NULL)
  {
    return;
  }
  
  cp = seqid_str;
  while (*cp != 0)
  {
    if (*cp == '|')
    {
      *cp = '_';
    }
    cp++;
  }
}

extern SeqEntryPtr MakeSequinDataFromAlignmentEx (TAlignmentFilePtr afp, Uint1 moltype, Boolean check_ids) 
{
  SeqIdPtr    PNTR sip_list;
  SeqIdPtr    PNTR sip_prev;
  SeqAnnotPtr sap = NULL;
  SeqAlignPtr salp_list, salp_last;
  ValNodePtr  PNTR seqvnp;
  SeqEntryPtr sep_list;
  SeqEntryPtr sep, sep_prev;
  SeqIdPtr    sip;
  ValNodePtr  vnp;
  Int4        index, curr_seg, num_sets;
  BioseqPtr   bsp;
  CharPtr     tmp_id_str;
  MsgAnswer   ans;
  Int4Ptr      segs_per_set = NULL;
  Int4Ptr      segs_per_aln = NULL;
  Boolean      found_empty_seg = FALSE;
  CharPtr      seq_data = NULL;

  if (afp == NULL) return NULL;
  
  if (afp->num_sequences == 0) return NULL;
  if (afp->num_segments < 1) return NULL;
  
  sip_list = (SeqIdPtr PNTR) MemNew (afp->num_segments * sizeof (SeqIdPtr));  
  sip_prev = (SeqIdPtr PNTR) MemNew (afp->num_segments * sizeof (SeqIdPtr));  
  seqvnp = (ValNodePtr PNTR) MemNew (afp->num_segments * sizeof (ValNodePtr));  
  segs_per_set = (Int4Ptr) MemNew (sizeof (Int4Ptr) * afp->num_sequences);
  segs_per_aln = (Int4Ptr) MemNew (sizeof (Int4Ptr) * afp->num_segments);
  if (sip_list == NULL || sip_prev == NULL || seqvnp == NULL
      || segs_per_set == NULL || segs_per_aln == NULL)
  {
    MemFree (sip_list);
    MemFree (sip_prev);
      MemFree (seqvnp);
      MemFree (segs_per_set);
      MemFree (segs_per_aln);
      return NULL;
  }
 
  for (curr_seg = 0; curr_seg < afp->num_segments; curr_seg ++)
  {
    sip_list [curr_seg] = NULL;
    sip_prev [curr_seg] = NULL;
      seqvnp [curr_seg] = NULL;
      segs_per_aln [curr_seg] = 0;
  }

  sep_list = NULL;
  sep_prev = NULL;
  curr_seg = 0;

  for (index = 0; index < afp->num_sequences; index++) {
    seq_data = AlignmentStringToSequenceString (afp->sequences [index], moltype);
    if (StringHasNoText (seq_data))
    {
      found_empty_seg = TRUE;
    }
    else
    {
      sip = MakeSeqID (afp->ids [index]);
      if (sip == NULL && StringChr (afp->ids [index], '|') != NULL)
      {
        ReplacePipesWithUnderscores (afp->ids [index]);
        sip = MakeSeqID (afp->ids [index]);
      }
      if (sip != NULL)
      {
        sip->next = SeqIdFree (sip->next);
      }
      if (check_ids && StringNCmp (afp->ids[index], "acc", 3) != 0)
      {
        bsp = BioseqFind (sip);
        if (bsp == NULL)
        {
          sip = SeqIdFree (sip);
          tmp_id_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (afp->ids [index]) + 4));
          sprintf (tmp_id_str, "gb|%s", afp->ids [index]);
          sip = MakeSeqID (tmp_id_str);
          MemFree (tmp_id_str);
          bsp = BioseqFind (sip);
        }
        if (bsp == NULL)
        {
          ans = Message (MSG_YN, "Can't find sequence %s in set - is this a far pointer?", afp->ids[index]);
          if (ans == ANS_YES)
          {
            sip = SeqIdFree (sip);
            sip = GetFarPointerID (afp->ids [index]);
          }
          else
          {
            sip = SeqIdFree (sip);
            sip = MakeSeqID (afp->ids [index]);
          }
          if (sip != NULL)
          {
            sip->next = SeqIdFree (sip->next);
          }
        }
      }

      sep = SequenceStringToSeqEntry (seq_data, sip, moltype);
      if (sep != NULL) {
        if (sep_list == NULL) {
          sep_list = sep;
        } else {
          sep_prev->next = sep;
        }
        sep_prev = sep;
        vnp = ValNodeNew (seqvnp[curr_seg]);
        if (seqvnp[curr_seg] == NULL) seqvnp[curr_seg] = vnp;
        vnp->data.ptrvalue = afp->sequences [index];
      
        /* only add SeqID to list if adding segment */
        if (sip_prev[curr_seg] == NULL) {
          sip_list[curr_seg] = sip;
        } else {
          sip_prev[curr_seg]->next = sip;
        }
        sip_prev[curr_seg] = sip;
        
        /* add to totals for this set and for this alignment */
        segs_per_set [index / afp->num_segments] ++;
        segs_per_aln [index % afp->num_segments] ++;
      }
    }
    seq_data = MemFree (seq_data);
    curr_seg ++;
    if (curr_seg >= afp->num_segments) 
    {
      curr_seg = 0;
    }
  }

  if (found_empty_seg)
  {
    Boolean   indexerVersion;
    MsgAnswer ans = ANS_YES;
    
    if (afp->num_segments > 1)
    {
      indexerVersion = (Boolean) (GetAppProperty ("InternalNcbiSequin") != NULL);
      if (indexerVersion)
      {
        ans = Message (MSG_YN, "This alignment of segmented sets contains a segment that is all gaps - do you wish to continue?");
      }
    }
    else
    {
      Message (MSG_ERROR, "This alignment contains a sequence that is all gaps.");
      ans = ANS_NO;
    }
    if (ans == ANS_NO)
    {
      for (curr_seg = 0; curr_seg < afp->num_segments; curr_seg ++)
      {
        ValNodeFree (seqvnp [curr_seg]);
      }
      MemFree (seqvnp);
      MemFree (sip_list);
      MemFree (sip_prev);
      MemFree (segs_per_set);
      MemFree (segs_per_aln);
      sep_list = SeqEntryFree (sep_list);
      return NULL;
    }
  }

  
  if (afp->num_segments == 1) 
  {
    sap = LocalAlignToSeqAnnotDimn (seqvnp[0], sip_list[0], NULL, afp->num_sequences,
                                    0, NULL, FALSE);
    sep_list = make_seqentry_for_seqentry (sep_list);
    SeqAlignAddInSeqEntry (sep_list, sap);      
  } 
  else 
  {
    sep_list = MakeSegmentedSetFromAlignment (sep_list, afp, moltype, segs_per_set);
    sep_list = make_seqentry_for_seqentry (sep_list);
    num_sets = afp->num_sequences / afp->num_segments;
    salp_list = NULL;
    salp_last = NULL;

    for (curr_seg = 0; curr_seg < afp->num_segments; curr_seg++)
    {      
      sap = LocalAlignToSeqAnnotDimn (seqvnp[curr_seg], sip_list[curr_seg], NULL, segs_per_aln [curr_seg],
                                    0, NULL, FALSE);
      if (sap != NULL)
      {
        SeqAlignAddInSeqEntry (sep_list, sap);
      }
    }
  }

  for (curr_seg = 0; curr_seg < afp->num_segments; curr_seg ++)
  {
    ValNodeFree (seqvnp [curr_seg]);
  }
  MemFree (seqvnp);
  MemFree (sip_list);
  MemFree (sip_prev);
  MemFree (segs_per_set);
  MemFree (segs_per_aln);

  AddDefLinesToAlignmentSequences (afp, sep_list);

  return sep_list;
}

extern SeqEntryPtr MakeSequinDataFromAlignment (TAlignmentFilePtr afp, Uint1 moltype) 
{
  return MakeSequinDataFromAlignmentEx (afp, moltype, FALSE);
}

/* Create sequences and alignment annotation */

/**********************************************************/
extern SeqEntryPtr make_seqentry_for_seqentry (SeqEntryPtr sep)
{
  SeqEntryPtr  sep1 = NULL,
               tmp;
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
 
  if (sep == NULL) return NULL;

  if (! IS_Bioseq (sep) && ! IS_Bioseq_set (sep)) {
    return sep;
  } else if (sep->next == NULL) {
    return sep;
  } else if ((bssp = BioseqSetNew ()) == NULL) {
    return sep;
  } else {
    bssp->_class = 14;
    bssp->seq_set = sep;
    sep1 = SeqEntryNew ();
    sep1->choice = 2;
    sep1->data.ptrvalue = bssp;
    SeqMgrLinkSeqEntry (sep1, 0, NULL);
          
    for (tmp = bssp->seq_set; tmp!=NULL; tmp=tmp->next) {
      if (IS_Bioseq(tmp)) {
        bsp = (BioseqPtr) tmp->data.ptrvalue;
        ObjMgrConnect (OBJ_BIOSEQ, (Pointer) bsp, OBJ_BIOSEQSET, (Pointer) bssp);
      }
    }
  }
  return sep1;
}

/* These two functions are used for removing mRNAs that overlap pseudo misc_feats
 * and marking genes that overlap pseudo misc_feats as pseudo.
 */
static void PseudoMiscFeatProcessingCallback (SeqFeatPtr sfp, Pointer userdata)
{
  SeqFeatPtr        gene, mRNA;
  SeqMgrFeatContext gcontext, mcontext;
  
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_misc_feature) return;
  /* we only want to process misc_feats if the pseudo flag is set or the 
   * comment contains the word "pseudogene".
   */
#if 0
  if (!sfp->pseudo && StringISearch (sfp->comment, "pseudogene") == NULL) return;
#endif

  gene = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
  mRNA = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_mRNA, NULL, 0, NULL,
                                      RANGE_MATCH, &mcontext);
  if (gene != NULL)
  {
      gene->pseudo = TRUE;
  }
  if (mRNA != NULL && mRNA->product == NULL) /* only delete mRNAs without products */
  {
      mRNA->idx.deleteme = TRUE;
  }  
}

extern void ProcessPseudoMiscFeatsForEntityID (Uint2 entityID)
{
  SeqEntryPtr sep;
  
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return;
  
  VisitFeaturesInSep (sep, (Pointer) NULL, PseudoMiscFeatProcessingCallback);
  DeleteMarkedObjects (entityID, 0, NULL);      
}

/* These three functions are used for converting pseudo CDSs to misc_features. */
NLM_EXTERN Boolean ConvertOnePseudoCDSToMiscFeatEx (SeqFeatPtr sfp, Boolean remove_product)
{
  BioseqPtr  bsp;
  SeqFeatPtr new_sfp;
  ImpFeatPtr ifp;
  
  if (sfp == NULL || (sfp->data.choice != SEQFEAT_CDREGION) || (! sfp->pseudo)) return FALSE;
  
  bsp = BioseqFindFromSeqLoc (sfp->location);  
  if (bsp == NULL) return FALSE;
  ifp = ImpFeatNew ();
  if (ifp == NULL) return FALSE;
  new_sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_IMP, sfp->location);
  if (new_sfp == NULL) 
  {
      ImpFeatFree (ifp);
      return FALSE;
  }
  new_sfp->data.value.ptrvalue = (Pointer) ifp;
  ifp->key = StringSave ("misc_feature");
  new_sfp->comment = sfp->comment;
  sfp->comment = NULL;
  new_sfp->qual = sfp->qual;
  sfp->qual = NULL;
  
  if (remove_product && sfp->product != NULL)
  {
      bsp = BioseqFindFromSeqLoc (sfp->product);
      sfp->product = SeqLocFree (sfp->product);
      bsp->idx.deleteme = TRUE;
  }
  sfp->idx.deleteme = TRUE;
  return TRUE;
}


extern Boolean ConvertOnePseudoCDSToMiscFeat (SeqFeatPtr sfp)
{
  return ConvertOnePseudoCDSToMiscFeatEx (sfp, TRUE);
}

static void ConvertPseudoCDSToMiscFeatCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ConvertOnePseudoCDSToMiscFeat (sfp);
}

extern void ConvertPseudoCDSToMiscFeatsForEntityID (Uint2 entityID)
{
  SeqEntryPtr sep;
  
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return;
  
  VisitFeaturesInSep (sep, (Pointer) NULL, ConvertPseudoCDSToMiscFeatCallback);
  DeleteMarkedObjects (entityID, 0, NULL);      
}

typedef struct alignmentforbsp
{
  BioseqPtr   bsp;
  SeqAlignPtr salp_list;
  SeqAlignPtr salp_last;
  ValNodePtr  seq_annot_list;
} AlignmentForBspData, PNTR AlignmentForBspPtr;

static void FindAlignmentsForBioseqCallback (SeqAnnotPtr sap, Pointer userdata)
{
  AlignmentForBspPtr   afbp;
  SeqAlignPtr          salp;
  SeqIdPtr             sip;
  Boolean              found = FALSE;

  if (sap == NULL || sap->type != 2 || userdata == NULL) 
  {
    return;
  }
  afbp = (AlignmentForBspPtr) userdata;
  if (afbp->bsp == NULL)
  {
    return;
  }
  salp = (SeqAlignPtr) sap->data;
  if (salp == NULL) return;
  for (sip = afbp->bsp->id; sip != NULL && !found; sip = sip->next)
  {
    if (SeqAlignFindSeqId (salp, sip))
    {
      salp = AlnMgr2DupAlnAndIndexes(salp);
      AlnMgr2IndexSeqAlign(salp);
      if (afbp->salp_last == NULL)
      {
        afbp->salp_list = salp; 
      }
      else
      {
        afbp->salp_last->next = salp;
      }
      afbp->salp_last = salp;
      found = TRUE;
    }
  }
}

extern SeqAlignPtr FindAlignmentsForBioseq (BioseqPtr bsp)
{
  SeqEntryPtr         topsep;
  AlignmentForBspData afbd;
  SeqLocPtr           slp;
  SeqIdPtr            sip;
  
  if (bsp == NULL) return NULL;
  topsep = GetTopSeqEntryForEntityID (bsp->idx.entityID);
  afbd.salp_list = NULL;
  afbd.salp_last = NULL;
  if (bsp->repr == Seq_repr_seg)
  {
    for (slp = bsp->seq_ext; slp != NULL; slp = slp->next)
    {
      sip = SeqLocId (slp);
      afbd.bsp = BioseqFind (sip);
      VisitAnnotsInSep (topsep, &afbd, FindAlignmentsForBioseqCallback);
    }
  }
  else
  {
    afbd.bsp = bsp;
    VisitAnnotsInSep (topsep, &afbd, FindAlignmentsForBioseqCallback);
  }
  
  return afbd.salp_list;
}

static void FindAlignSeqAnnotsForBioseqCallback (SeqAnnotPtr sap, Pointer userdata)
{
  AlignmentForBspPtr   afbp;
  SeqAlignPtr          salp;
  SeqIdPtr             sip;
  Boolean              found = FALSE;

  if (sap == NULL || sap->type != 2 || userdata == NULL) 
  {
    return;
  }
  afbp = (AlignmentForBspPtr) userdata;
  if (afbp->bsp == NULL)
  {
    return;
  }
  salp = (SeqAlignPtr) sap->data;
  if (salp == NULL) return;
  for (sip = afbp->bsp->id; sip != NULL && !found; sip = sip->next)
  {
    if (SeqAlignFindSeqId (salp, sip))
    {
      ValNodeAddPointer (&(afbp->seq_annot_list), 0, sap);
      found = TRUE;
    }
  }
}

extern ValNodePtr FindAlignSeqAnnotsForBioseq (BioseqPtr bsp)
{
  SeqEntryPtr         topsep;
  AlignmentForBspData afbd;
  SeqLocPtr           slp;
  SeqIdPtr            sip;
  
  if (bsp == NULL) return NULL;
  topsep = GetTopSeqEntryForEntityID (bsp->idx.entityID);
  afbd.salp_list = NULL;
  afbd.salp_last = NULL;
  afbd.seq_annot_list = NULL;
  if (bsp->repr == Seq_repr_seg)
  {
    for (slp = bsp->seq_ext; slp != NULL; slp = slp->next)
    {
      sip = SeqLocId (slp);
      afbd.bsp = BioseqFind (sip);
      VisitAnnotsInSep (topsep, &afbd, FindAlignSeqAnnotsForBioseqCallback);
    }
  }
  else
  {
    afbd.bsp = bsp;
    VisitAnnotsInSep (topsep, &afbd, FindAlignSeqAnnotsForBioseqCallback);
  }
  
  return afbd.seq_annot_list;
}

NLM_EXTERN void ChangeSeqIdToWorstID (SeqIdPtr sip)
{
  BioseqPtr       bsp;
  SeqIdPtr        id;
  Pointer         pnt;

  if (sip == NULL)
    return;
  bsp = BioseqFindCore (sip);
  if (bsp == NULL)
    return;
  id = SeqIdDup (SeqIdFindWorst (bsp->id));
  if (id == NULL)
    return;
  /* now remove SeqId contents to reuse SeqId valnode */
  pnt = sip->data.ptrvalue;
  switch (sip->choice) {
  case SEQID_LOCAL:            /* local */
    ObjectIdFree ((ObjectIdPtr) pnt);
    break;
  case SEQID_GIBBSQ:           /* gibbseq */
  case SEQID_GIBBMT:           /* gibbmt */
    break;
  case SEQID_GIIM:             /* giimid */
    GiimFree ((GiimPtr) pnt);
    break;
  case SEQID_GENBANK:          /* genbank */
  case SEQID_EMBL:             /* embl */
  case SEQID_PIR:              /* pir   */
  case SEQID_SWISSPROT:        /* swissprot */
  case SEQID_OTHER:            /* other */
  case SEQID_DDBJ:
  case SEQID_PRF:
  case SEQID_TPG:
  case SEQID_TPE:
  case SEQID_TPD:
  case SEQID_GPIPE:
    TextSeqIdFree ((TextSeqIdPtr) pnt);
    break;
  case SEQID_PATENT:           /* patent seq id */
    PatentSeqIdFree ((PatentSeqIdPtr) pnt);
    break;
  case SEQID_GENERAL:          /* general */
    DbtagFree ((DbtagPtr) pnt);
    break;
  case SEQID_GI:               /* gi */
    break;
  case SEQID_PDB:
    PDBSeqIdFree ((PDBSeqIdPtr) pnt);
    break;
  }
  sip->choice = id->choice;
  sip->data.ptrvalue = id->data.ptrvalue;
  SeqIdStripLocus (sip);
}

NLM_EXTERN void ChangeSeqLocToWorstID (SeqLocPtr slp)
{
  SeqLocPtr       loc;
  PackSeqPntPtr   psp;
  SeqBondPtr      sbp;
  SeqIntPtr       sinp;
  SeqIdPtr        sip;
  SeqPntPtr       spp;

  while (slp != NULL) {
    switch (slp->choice) {
    case SEQLOC_NULL:
      break;
    case SEQLOC_EMPTY:
    case SEQLOC_WHOLE:
      sip = (SeqIdPtr) slp->data.ptrvalue;
      ChangeSeqIdToWorstID (sip);
      break;
    case SEQLOC_INT:
      sinp = (SeqIntPtr) slp->data.ptrvalue;
      if (sinp != NULL) {
        sip = sinp->id;
        ChangeSeqIdToWorstID (sip);
      }
      break;
    case SEQLOC_PNT:
      spp = (SeqPntPtr) slp->data.ptrvalue;
      if (spp != NULL) {
        sip = spp->id;
        ChangeSeqIdToWorstID (sip);
      }
      break;
    case SEQLOC_PACKED_PNT:
      psp = (PackSeqPntPtr) slp->data.ptrvalue;
      if (psp != NULL) {
        sip = psp->id;
        ChangeSeqIdToWorstID (sip);
      }
      break;
    case SEQLOC_PACKED_INT:
    case SEQLOC_MIX:
    case SEQLOC_EQUIV:
      loc = (SeqLocPtr) slp->data.ptrvalue;
      while (loc != NULL) {
        ChangeSeqLocToWorstID (loc);
        loc = loc->next;
      }
      break;
    case SEQLOC_BOND:
      sbp = (SeqBondPtr) slp->data.ptrvalue;
      if (sbp != NULL) {
        spp = (SeqPntPtr) sbp->a;
        if (spp != NULL) {
          sip = spp->id;
          ChangeSeqIdToWorstID (sip);
        }
        spp = (SeqPntPtr) sbp->b;
        if (spp != NULL) {
          sip = spp->id;
          ChangeSeqIdToWorstID (sip);
        }
      }
      break;
    case SEQLOC_FEAT:
      break;
    default:
      break;
    }
    slp = slp->next;
  }
}

/* This function will remove DenDiag and pairwise alignments if they contain
 * the sequence identified by sip, otherwise it will remove the sequence from
 * the alignment.
 */
static SeqAlignPtr RemoveOneSequenceFromAlignment (SeqIdPtr sip, SeqAlignPtr salphead)
{
  Uint4       seqid_order;
  SeqIdPtr    tmpsip;
  SeqAlignPtr salp, salp_next, prev_salp, remove_salp, last_remove;
  
  if (!FindSeqIdinSeqAlign (salphead, sip)) return salphead;
  
  salp = salphead;
  prev_salp = NULL;
  remove_salp = NULL;
  last_remove = NULL;
  while (salp != NULL)
  {
    salp_next = salp->next;
    tmpsip = SeqIdPtrFromSeqAlign (salp);
    seqid_order = SeqIdOrderInBioseqIdList(sip, tmpsip);
    if (seqid_order == 0)
    {
      /* do nothing for this subalignment */
      prev_salp = salp;
    }
    else if (salp->dim == 2 || salphead->segtype ==1)
    {
      /* This is for a pairwise alignment or a DENDIAG alignment */
      if (prev_salp == NULL)
      {
          salphead = salp->next;
      }
      else
      {
          prev_salp->next = salp->next;
      }
      /* save the alignments that we want to free in a list and get rid of them
       * at the end - freeing them beforehand causes problems with listing the
       * IDs in the alignment.
       */
      salp->next = NULL;
      if (remove_salp == NULL)
      {
          remove_salp = salp;
      }
      else
      {
          last_remove->next = salp;
      }
      last_remove = salp;
    }
    else 
    {
      SeqAlignBioseqDeleteById (salphead, sip);  
      prev_salp = salp;
    }
    salp = salp_next;
  }
  /* Now we can free the alignment */
  SeqAlignFree (remove_salp);
  return salphead;
}

static void RemoveSequenceFromAlignmentsCallback (SeqAnnotPtr sap, Pointer userdata)
{
  SeqAlignPtr salp;
  SeqIdPtr    sip;

  if (sap == NULL || sap->type != 2 || userdata == NULL) return;
  salp = (SeqAlignPtr) sap->data;
  if (salp == NULL) return;
  sip = (SeqIdPtr) userdata;
  sap->data = RemoveOneSequenceFromAlignment (sip, salp);
  /* if we've deleted all of the alignments, get rid of the annotation as well */
  if (sap->data == NULL)
  {
      sap->idx.deleteme = TRUE;
  }
}

typedef struct checkforremovesequencefromalignments
{
  Boolean  found_problem;
  SeqIdPtr sip;
} CheckForRemoveSequenceFromAlignmentsData, PNTR CheckForRemoveSequenceFromAlignmentsPtr;

/* This is the callback function for looking for pairwise alignments.
 * If we delete the first sequence in a pairwise alignment, we end up deleting
 * the entire alignment because that sequence is paired with every other sequence.
 */
static void CheckForRemoveSequenceFromAlignmentsProblemsCallback (SeqAnnotPtr sap, Pointer userdata)
{
  CheckForRemoveSequenceFromAlignmentsPtr p;
  SeqAlignPtr salphead, salp;
  Uint4       seqid_order;
  SeqIdPtr    tmpsip;
  
  if (sap == NULL || sap->type != 2
      || (p = (CheckForRemoveSequenceFromAlignmentsPtr)userdata) == NULL
      || p->found_problem)
  {
      return;
  }
  salphead = (SeqAlignPtr) sap->data;
  if (salphead == NULL) return;
  
  if (!FindSeqIdinSeqAlign (salphead, p->sip))
  {
      return;
  }
  for (salp = salphead; salp != NULL; salp = salp->next)
  {
    tmpsip = SeqIdPtrFromSeqAlign (salp);
    seqid_order = SeqIdOrderInBioseqIdList(p->sip, tmpsip);
    if (seqid_order == 0)
    {
      continue;
    }
    else if (seqid_order == 1 && salp->dim == 2)
    {
      p->found_problem = TRUE;      
    }
  }
}

extern Boolean IsSequenceFirstInPairwise (SeqEntryPtr sep, SeqIdPtr sip)
{
  CheckForRemoveSequenceFromAlignmentsData data;
  
  if (sep == NULL || sip == NULL)
  {
    return FALSE;
  }
  
    data.sip = sip;
    data.found_problem = FALSE;
  
  VisitAnnotsInSep (sep, (Pointer) &data, CheckForRemoveSequenceFromAlignmentsProblemsCallback);
  return data.found_problem;
}

extern Boolean RemoveSequenceFromAlignments (SeqEntryPtr sep, SeqIdPtr sip)
{
  if (sep == NULL || sip == NULL)
  {
    return FALSE;
  }
  if (IsSequenceFirstInPairwise (sep, sip))
  {
    return FALSE;
  }
  VisitAnnotsInSep (sep, (Pointer) sip, RemoveSequenceFromAlignmentsCallback);
  return TRUE;
}

static CharPtr evCategoryPrefix [] = {
  "",
  "COORDINATES: ",
  "DESCRIPTION: ",
  "EXISTENCE: ",
  NULL
};

static CharPtr inferencePrefix [] = {
  "",
  "similar to sequence",
  "similar to AA sequence",
  "similar to DNA sequence",
  "similar to RNA sequence",
  "similar to RNA sequence, mRNA",
  "similar to RNA sequence, EST",
  "similar to RNA sequence, other RNA",
  "profile",
  "nucleotide motif",
  "protein motif",
  "ab initio prediction",
  "alignment",
  NULL
};

static Int2 ValidateInferenceAccession (CharPtr str, Char chr, Boolean fetchAccn, Boolean has_fetch_function)

{
  Int2      accnv, rsult;
  Boolean   is_insd = FALSE, is_refseq = FALSE;
  ErrSev    sev;
  SeqIdPtr  sip;
  CharPtr   tmp;

  if (StringHasNoText (str)) return EMPTY_INFERENCE_STRING;

  rsult = VALID_INFERENCE;

  tmp = StringChr (str, chr);
  if (tmp != NULL) {
    *tmp = '\0';
    tmp++;
    TrimSpacesAroundString (str);
    TrimSpacesAroundString (tmp);
    if (StringDoesHaveText (tmp)) {
      if (StringICmp (str, "INSD") == 0) {
        is_insd = TRUE;
      }
      if (StringICmp (str, "RefSeq") == 0) {
        is_refseq = TRUE;
      }
      if (is_insd || is_refseq) {
        if (StringLen (tmp) > 3) {
          if (tmp [2] == '_') {
            if (is_insd) {
              rsult = BAD_ACCESSION_TYPE;
            }
          } else {
            if (is_refseq) {
              rsult = BAD_ACCESSION_TYPE;
            }
          }
        }
        accnv = ValidateAccnDotVer (tmp);
        if (accnv == -5 || accnv == -6) {
          rsult = BAD_INFERENCE_ACC_VERSION;
        } else if (accnv != 0) {
          rsult = BAD_INFERENCE_ACCESSION;
        } else if (fetchAccn) {
          sip = SeqIdFromAccessionDotVersion (tmp);
          sev = ErrGetMessageLevel ();
          ErrSetMessageLevel (SEV_ERROR);
          if (has_fetch_function && GetGIForSeqId (sip) == 0) {
            rsult = ACC_VERSION_NOT_PUBLIC;
          }
          ErrSetMessageLevel (sev);
          SeqIdFree (sip);
        }
      }
    }
    if (StringChr (tmp, ' ') != NULL) rsult = SPACES_IN_INFERENCE;
  } else {
    rsult = SINGLE_INFERENCE_FIELD;
  }

  return rsult;
}

static Char NextColonOrVerticalBar (CharPtr ptr)

{
  Char  ch = '\0';

  if (ptr == NULL) return ch;

  ch = *ptr;
  while (ch != '\0') {
    if (ch == ':' || ch == '|') return ch;
    ptr++;
    ch = *ptr;
  }

  return ch;
}

NLM_EXTERN Int2 ValidateInferenceQualifier (CharPtr val, Boolean fetchAccn)

{
  Int2           best, j, rsult, tmprsult;
  Char           ch;
  Boolean        has_fetch_function, same_species;
  size_t         len;
  CharPtr        nxt, ptr, rest, skip, str;
  ObjMgrProcPtr  ompp = NULL;

  if (StringHasNoText (val)) return EMPTY_INFERENCE_STRING;

  skip = NULL;
  for (j = 0; evCategoryPrefix [j] != NULL; j++) {
    len = StringLen (evCategoryPrefix [j]);
    if (StringNICmp (val, evCategoryPrefix [j], len) != 0) continue;
    skip = val + len;
  }
  if (skip != NULL) {
    val = skip;
  }

  for (j = 0; inferencePrefix [j] != NULL; j++) {
    len = StringLen (inferencePrefix [j]);
    if (StringNICmp (val, inferencePrefix [j], len) != 0) continue;
    rest = val + len;
    best = j;
  }

  if (best < 0 || inferencePrefix [best] == NULL) return BAD_INFERENCE_PREFIX;

  if (rest == NULL) return BAD_INFERENCE_BODY;

  same_species = FALSE;
  ch = *rest;
  while (IS_WHITESP (ch)) {
    rest++;
    ch = *rest;
  }
  if (StringNICmp (rest, "(same species)", 14) == 0) {
    same_species = TRUE;
    rest += 14;
  }
  ch = *rest;
  while (IS_WHITESP (ch) || ch == ':') {
    rest++;
    ch = *rest;
  }

  if (StringHasNoText (rest)) return BAD_INFERENCE_BODY;

  rsult = VALID_INFERENCE;
  if (same_species && best > 7) {
    rsult = SAME_SPECIES_MISUSED;
  }

  has_fetch_function = FALSE;
  while ((ompp = ObjMgrProcFindNext(NULL, OMPROC_FETCH, OBJ_SEQID, OBJ_SEQID, ompp)) != NULL) {
    if ((ompp->subinputtype == 0) && (ompp->suboutputtype == SEQID_GI)) {
      has_fetch_function = TRUE;
    }
  }

  str = StringSave (rest);

  if (best >= 1 && best <= 7) {
    ptr = str;
    while (ptr != NULL) {
      nxt = StringChr (ptr, ',');
      if (nxt != NULL) {
        *nxt = '\0';
        nxt++;
      }
      tmprsult = ValidateInferenceAccession (ptr, ':', fetchAccn, has_fetch_function);
      if (tmprsult != VALID_INFERENCE) {
        rsult = tmprsult;
      }
      ptr = nxt;
    }
  } else if (best == 12) {
    tmprsult = VALID_INFERENCE;
    ptr = StringRChr (str, ':');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }
    while (ptr != NULL) {
      nxt = StringChr (ptr, ',');
      if (nxt != NULL) {
        *nxt = '\0';
        nxt++;
      }
      ch = NextColonOrVerticalBar (ptr);
      tmprsult = ValidateInferenceAccession (ptr, ch, fetchAccn, has_fetch_function);
      if (tmprsult != VALID_INFERENCE) {
        rsult = tmprsult;
      }
      ptr = nxt;
    }
  }

  MemFree (str);

  return rsult;
}

extern void MergeFeatureIntervalsToParts (SeqFeatPtr sfp, Boolean ordered)
{
  BioseqPtr     bsp;
  CodeBreakPtr  cbp;
  CdRegionPtr   crp;
  Boolean       noLeft;
  Boolean       noRight;
  RnaRefPtr     rrp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  tRNAPtr       trna;
  
  if (sfp == NULL || sfp->location == NULL) return;
  sip = SeqLocId (sfp->location);
  if (sip == NULL) return;
  bsp = BioseqFind (sip);
  if (bsp == NULL) return;
  if (ISA_aa (bsp->mol)) return;
  CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
  slp = SegLocToPartsEx (bsp, sfp->location, ordered);
  if (slp == NULL) return;
  sfp->location = SeqLocFree (sfp->location);
  sfp->location = slp;
  FreeAllFuzz (sfp->location);
  SetSeqLocPartial (sfp->location, noLeft, noRight);
  sfp->partial = (sfp->partial || noLeft || noRight);
  switch (sfp->data.choice) {
    case SEQFEAT_CDREGION :
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL && crp->code_break != NULL) {
        for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
          slp = SegLocToPartsEx (bsp, cbp->loc, ordered);
          if (slp != NULL) {
            cbp->loc = SeqLocFree (cbp->loc);
            cbp->loc = slp;
            FreeAllFuzz (cbp->loc);
          }
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->type == 3 && rrp->ext.choice == 2) {
        trna = rrp->ext.value.ptrvalue;
        if (trna != NULL && trna->anticodon != NULL) {
          slp = SegLocToPartsEx (bsp, trna->anticodon, ordered);
          if (slp != NULL) {
            trna->anticodon = SeqLocFree (trna->anticodon);
            trna->anticodon = slp;
            FreeAllFuzz (trna->anticodon);
          }
        }
      }
      break;
    default :
      break;
  }  
}

extern void ExtendSingleGeneOnMRNA (BioseqPtr bsp, Pointer userdata)

{
  MolInfoPtr        mip;
  SeqDescrPtr       sdp;
  Boolean           is_mrna = FALSE, is_master_seq = FALSE, has_nulls = FALSE;
  SeqFeatPtr        gene = NULL;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext context;
  Int4              num_cds = 0;
  Int4              num_mrna = 0;
  SeqIdPtr          sip;
  SeqLocPtr         slp;
  Boolean           partial5, partial3;
  BioSourcePtr      biop;
  OrgRefPtr         orp;
  BioseqSetPtr      bssp;

  if (bsp == NULL || bsp->length == 0
      || !ISA_na (bsp->mol)) {
    return;
  }

  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
  if (sdp != NULL) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL && mip->biomol == MOLECULE_TYPE_MRNA) {
      is_mrna = TRUE;
    }
  }
  if (!is_mrna) {
    return;
  }

  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_source, NULL);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL) {
      if (biop->origin == ORG_ARTIFICIAL) {
        orp = biop->org;
        if (orp != NULL) {
          if (StringICmp (orp->taxname, "synthetic construct") == 0) return;
        }
      }
    }
  }

  if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) bsp->idx.parentptr;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_segset) {
      is_master_seq = TRUE;
    }
  }
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context)) {
    if (sfp->data.choice == SEQFEAT_GENE) {
      /* skip this sequence if it has more than one gene */
      if (gene == NULL) {
        gene = sfp;
      } else {
        return; 
      }
    } else if (sfp->data.choice == SEQFEAT_CDREGION) {
      num_cds++;
      /* skip this sequence if it has more than one coding region */
      if (num_cds > 1 && !is_master_seq) {
        return;
      }
    } else if (sfp->idx.subtype == FEATDEF_mRNA) {
      num_mrna++;
      /* skip this sequence if it has more than one mRNA */
      if (num_mrna > 1) return;
    }
  }

  if (gene != NULL && gene->location != NULL) {
    slp = gene->location;
    if (slp->choice != SEQLOC_INT) {
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        /* skip this sequence if it is multi-interval and EMBL or DDBJ */
        if (sip->choice == SEQID_EMBL || sip->choice == SEQID_DDBJ) return;
      }
    }
  }

  if (gene != NULL && BioseqFindFromSeqLoc (gene->location) == bsp) {
    CheckSeqLocForPartial (gene->location, &partial5, &partial3);
    has_nulls = LocationHasNullsBetween (gene->location);
    /* gene should cover entire length of sequence */
    slp = SeqLocIntNew (0, bsp->length - 1, SeqLocStrand (gene->location), SeqIdFindBest (bsp->id, 0));
    SetSeqLocPartial (slp, partial5, partial3);
    gene->location = SeqLocFree (gene->location);
    gene->location = slp;
    if (is_master_seq) {
      MergeFeatureIntervalsToParts (gene, has_nulls);
    }
  }
}

/* Functions for the Discrepancy Report */


static Boolean IsProdBiomol (Uint1 biomol)
{
  if (biomol == MOLECULE_TYPE_MRNA
      || biomol == MOLECULE_TYPE_NCRNA 
      || biomol == MOLECULE_TYPE_RRNA
      || biomol == MOLECULE_TYPE_PRE_MRNA
      || biomol == MOLECULE_TYPE_TRNA)  {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean IsmRNASequenceInGenProdSet (BioseqPtr bsp)
{
  SeqMgrDescContext dcontext;
  BioseqSetPtr      bssp;
  SeqDescrPtr sdp;
  MolInfoPtr  mip;
  Boolean rval = FALSE;

  if (bsp != NULL && bsp->mol == Seq_mol_rna && bsp->idx.parentptr != NULL && bsp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) bsp->idx.parentptr;
    if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
      sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
      if (sdp != NULL && sdp->data.ptrvalue != NULL && sdp->choice == Seq_descr_molinfo) {
        mip = (MolInfoPtr) sdp->data.ptrvalue;
        rval = IsProdBiomol (mip->biomol);
      }
    } else if (bssp->_class == BioseqseqSet_class_nuc_prot && bssp->idx.parentptr != NULL && bssp->idx.parenttype == OBJ_BIOSEQSET) {
      bssp = (BioseqSetPtr) bssp->idx.parentptr;
      if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
        sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
        if (sdp != NULL && sdp->data.ptrvalue != NULL && sdp->choice == Seq_descr_molinfo) {
          mip = (MolInfoPtr) sdp->data.ptrvalue;
          rval = IsProdBiomol (mip->biomol);
        }
      }
    }
  }
  return rval;
}


typedef struct skipmrnafeaturesingenprodset {
  Pointer userdata;
  VisitFeaturesFunc callback;
} SkipmRNAFeaturesInGenProdSetData, PNTR SkipmRNAFeaturesInGenProdSetPtr;

static void VisitGenProdSetFeaturesCallback (SeqFeatPtr sfp, Pointer userdata)
{
  SkipmRNAFeaturesInGenProdSetPtr p;

  if (sfp == NULL || userdata == NULL) {
    return;
  }

  p = (SkipmRNAFeaturesInGenProdSetPtr) userdata;
  if (p->callback == NULL) {
    return;
  }

  if (!IsmRNASequenceInGenProdSet(BioseqFindFromSeqLoc (sfp->location))) {
    (p->callback) (sfp, p->userdata);
  }
}


extern void VisitGenProdSetFeatures (SeqEntryPtr sep, Pointer userdata, VisitFeaturesFunc callback)
{
  SkipmRNAFeaturesInGenProdSetData d;

  d.callback = callback;
  d.userdata = userdata;
  VisitFeaturesInSep (sep, &d, VisitGenProdSetFeaturesCallback);
}


extern ClickableItemPtr 
NewClickableItem 
(Uint4           clickable_item_type,
 CharPtr         description_fmt,
 ValNodePtr      item_list)
{
  ClickableItemPtr dip;
  
  dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  if (dip != NULL)
  {
    dip->clickable_item_type = clickable_item_type;
    dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (description_fmt) + 15));
    sprintf (dip->description, description_fmt, ValNodeLen (item_list));
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    dip->item_list = item_list;
    dip->subcategories = NULL;
    dip->expanded = FALSE;
    dip->level = 0;
  }
  return dip;  
}


extern ClickableItemPtr 
NewClickableItemNoList 
(Uint4           clickable_item_type,
 CharPtr         description)
{
  ClickableItemPtr dip;
  
  dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  if (dip != NULL)
  {
    dip->clickable_item_type = clickable_item_type;
    dip->description = StringSave (description);
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    dip->item_list = NULL;
    dip->subcategories = NULL;
    dip->expanded = FALSE;
    dip->level = 0;
  }
  return dip;  
}


extern ValNodePtr ClickableItemObjectListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;
  ObjValNodePtr ovn;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next= NULL;
    if (vnp->extended != 0) {
      ovn = (ObjValNodePtr) vnp;
      ovn->idx.scratch = FieldTypeListFree (ovn->idx.scratch);
    }
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


extern ValNodePtr ClickableItemObjectListCopy (ValNodePtr orig)
{
  ValNodePtr cpy = NULL, vnp_last = NULL, vnp;
  ObjValNodePtr ovn, ovn_cpy;

  while (orig != NULL) {
    if (orig->extended != 0) {
      ovn = (ObjValNodePtr) orig;
      ovn_cpy = (ObjValNodePtr) SeqDescrNew (NULL);
      ovn_cpy->idx.scratch = FieldTypeListCopy (ovn->idx.scratch);
      ovn_cpy->vn.choice = ovn->vn.choice;
      ovn_cpy->vn.data.ptrvalue = ovn->vn.data.ptrvalue;
      if (vnp_last == NULL) {
        cpy = (ValNodePtr) ovn_cpy;
      } else {
        vnp_last->next = (ValNodePtr) ovn_cpy;
      }
      vnp_last = (ValNodePtr) ovn_cpy;
    } else {
      vnp = ValNodeNew (NULL);
      vnp->choice = orig->choice;
      vnp->data.ptrvalue = orig->data.ptrvalue;
      if (vnp_last == NULL) {
        cpy = vnp;
      } else {
        vnp_last->next = vnp;
      }
      vnp_last = vnp;
    }
    orig = orig->next;
  }
  return cpy;
}

static ValNodePtr MakeObjectListWithFields (ValNodePtr item_list, ValNodePtr field_list)
{
  ValNodePtr vnp, vnp_last = NULL, extended_item_list = NULL;
  ObjValNodePtr ovn;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    ovn = (ObjValNodePtr) SeqDescrNew (NULL);
    ovn->vn.choice = vnp->choice;
    ovn->vn.data.ptrvalue = vnp->data.ptrvalue;
    ovn->idx.scratch = FieldTypeListCopy (field_list);
    if (vnp_last == NULL) {
      extended_item_list = (ValNodePtr)ovn;
    } else {
      vnp_last->next = (ValNodePtr)ovn;
    }
    vnp_last = (ValNodePtr)ovn;
  }
  return extended_item_list;
}


extern ClickableItemPtr ClickableItemFree (ClickableItemPtr cip)
{
  if (cip != NULL)
  {
    cip->description = MemFree (cip->description);
    if (cip->datafree_func != NULL)
    {
      (cip->datafree_func) (cip->callback_data);
    }
    cip->item_list = ClickableItemObjectListFree (cip->item_list);
  
    cip->subcategories = FreeClickableList (cip->subcategories);
    cip = MemFree (cip);
  }
  return cip;
}


extern ValNodePtr FreeClickableList (ValNodePtr list)
{
  ValNodePtr       list_next;
  
  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = ClickableItemFree (list->data.ptrvalue);
    list = ValNodeFree (list);
    list = list_next;
  }
  return list;
}


extern Boolean AnyDiscrepanciesChosen (ValNodePtr cip_list)
{
  ClickableItemPtr cip;
  Boolean          any_chosen = FALSE;

  while (cip_list != NULL && !any_chosen) {
    cip = (ClickableItemPtr) cip_list->data.ptrvalue;
    if (cip != NULL
        && (cip->chosen 
            || (cip->expanded && AnyDiscrepanciesChosen (cip->subcategories)))) {
      any_chosen = TRUE;
    }
    cip_list = cip_list->next;
  }
  return any_chosen;
}


NLM_EXTERN void ChooseAllDiscrepancies (ValNodePtr cip_list)
{
  ClickableItemPtr cip;

  while (cip_list != NULL) {
    cip = cip_list->data.ptrvalue;
    if (cip != NULL) {
      cip->chosen = TRUE;
    }
    cip_list = cip_list->next;
  }
}


NLM_EXTERN int LIBCALLBACK SortVnpByClickableItemChosen (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  ClickableItemPtr cip1, cip2;
  Int4        rval = 0;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  if (vnp1->data.ptrvalue == NULL || vnp2->data.ptrvalue == NULL) return 0;
  
  cip1 = vnp1->data.ptrvalue;
  cip2 = vnp2->data.ptrvalue;

  if (cip1->chosen && cip2->chosen) {
    rval = 0;
  } else if (!cip1->chosen && !cip2->chosen) {
    rval = 0;
  } else if (cip1->chosen) {
    rval = -1;
  } else if (cip2->chosen) {
    rval = 1;
  }
  return rval;
}


static CharPtr GetClickableItemDescription (ValNodePtr vnp)
{
  ClickableItemPtr cip;

  if (vnp == NULL || vnp->data.ptrvalue == NULL) return NULL;
  cip = (ClickableItemPtr) vnp->data.ptrvalue;
  return cip->description;
}


extern void ExpandClickableItemList (ValNodePtr vnp)
{
  ClickableItemPtr cip;
  
  while (vnp != NULL) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    cip->expanded = TRUE;
    ExpandClickableItemList (cip->subcategories);
    vnp = vnp->next;
  }
}


extern void ContractClickableItemList (ValNodePtr vnp)
{
  ClickableItemPtr cip;
  
  while (vnp != NULL) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    cip->expanded = FALSE;
    ContractClickableItemList (cip->subcategories);
    vnp = vnp->next;
  }
}


extern int LIBCALLBACK SortVnpByClickableItemDescription (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1;
  CharPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      str1 = GetClickableItemDescription (vnp1);
      str2 = GetClickableItemDescription (vnp2);
      if (str1 != NULL && str2 != NULL) {
        return StringICmp (str1, str2);
      }
    }
  }
  return 0;
}


/* utility functions for the discrepancy report tests */
static void ValNodeLinkCopy (ValNodePtr PNTR list1, ValNodePtr list2)
{
  if (list1 == NULL) return;
  while (list2 != NULL)
  {
    ValNodeAddPointer (list1, list2->choice, list2->data.ptrvalue);
    list2 = list2->next;
  }
}


static Boolean ValNodeStringListMatch (ValNodePtr vnp1, ValNodePtr vnp2)
{
  if (vnp1 == NULL && vnp2 == NULL)
  {
    return TRUE;
  }
  else if (vnp1 == NULL || vnp2 == NULL)
  {
    return FALSE;
  }
  else if (StringCmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue) != 0)
  {
    return FALSE;
  }
  else
  {
    return ValNodeStringListMatch (vnp1->next, vnp2->next);
  }
}


static ValNodePtr ItemListFromSubcategories (ValNodePtr subcategories)
{
  ValNodePtr       vnp;
  ClickableItemPtr cip;
  ValNodePtr       item_list = NULL;

  for (vnp = subcategories; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL) {
      ValNodeLink (&item_list, ClickableItemObjectListCopy(cip->item_list));
    }
  }
  return item_list;
}


static void ValNodeExtractMatch (ValNodePtr PNTR list, ValNodePtr match)
{
  ValNodePtr vnp_prev = NULL, vnp_next, vnp;

  if (list == NULL) return;

  for (vnp = *list; vnp != NULL; vnp = vnp_next) {
    vnp_next = vnp->next;
    if (vnp->choice == match->choice && vnp->data.ptrvalue == match->data.ptrvalue) {
      if (vnp_prev == NULL) {
        *list = vnp_next;
      } else {
        vnp_prev->next = vnp_next;
      }
      vnp->next = NULL;
      vnp = ValNodeFree (vnp);
    } else {
      vnp_prev = vnp;
    }
  }
}


NLM_EXTERN void RemoveDuplicateItems (ValNodePtr PNTR item_list)
{
  ValNodePtr vnp;

  if (item_list == NULL) {
    return;
  }
  for (vnp = *item_list; vnp != NULL && vnp->next != NULL; vnp = vnp->next) {
    ValNodeExtractMatch (&(vnp->next), vnp);
  }
}


extern GlobalDiscrepancyPtr 
GlobalDiscrepancyNew (CharPtr str, Uint1 data_choice, Pointer data)
{
  GlobalDiscrepancyPtr g;

  g = (GlobalDiscrepancyPtr) MemNew (sizeof (GlobalDiscrepancyData));
  g->str = StringSave (str);
  g->data_choice = data_choice;
  if (g->data_choice == 0) {
    g->data = StringSave (data);
  } else {
    g->data = data;
  }
  return g;
}


extern GlobalDiscrepancyPtr GlobalDiscrepancyFree (GlobalDiscrepancyPtr g)
{
  if (g != NULL) {
    g->str = MemFree (g->str);
    if (g->data_choice == 0) {
      g->data = MemFree (g->data);
    }
    g = MemFree (g);
  }
  return g;
}


extern ValNodePtr FreeGlobalDiscrepancyList (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = GlobalDiscrepancyFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


extern void ConvertGlobalDiscrepancyToText (GlobalDiscrepancyPtr g, Boolean use_feature_fmt, CharPtr filename)
{
  ValNode vn;
  ValNodePtr list_copy;

  if (g == NULL || g->data_choice == 0) return;
  
  MemSet (&vn, 0, sizeof (ValNode));
  vn.choice = g->data_choice;
  vn.data.ptrvalue = g->data;
  vn.next = NULL;

  g->data_choice = 0;
  if (use_feature_fmt) {  
    list_copy = ReplaceDiscrepancyItemWithFeatureTableStrings (&vn);
    g->data = list_copy->data.ptrvalue;
    list_copy = ValNodeFree (list_copy);
  } else {

    g->data = GetDiscrepancyItemTextEx (&vn, filename);
  }
}


extern void ConvertGlobalDiscrepancyListToText (ValNodePtr vnp, Boolean use_feature_fmt, CharPtr filename)
{
  while (vnp != NULL) {
    ConvertGlobalDiscrepancyToText (vnp->data.ptrvalue, use_feature_fmt, filename);
    vnp = vnp->next;
  }
}


extern ValNodePtr GetGlobalDiscrepancyItem (GlobalDiscrepancyPtr g)
{
  ValNodePtr rval = NULL;
  if (g != NULL) {
    rval = ValNodeNew (NULL);
    rval->choice = g->data_choice;
    if (rval->choice == 0) {
      rval->data.ptrvalue = StringSave (g->data);
    } else {
      rval->data.ptrvalue = g->data;
    }
  }
  return rval;
}


extern CharPtr GetGlobalDiscrepancyStr (GlobalDiscrepancyPtr g)
{
  CharPtr rval = NULL;
  if (g != NULL) {
    rval = g->str;
  }
  return rval;
}


NLM_EXTERN int LIBCALLBACK SortVnpByGlobalDiscrepancyString (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  GlobalDiscrepancyPtr g1, g2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      g1 = (GlobalDiscrepancyPtr) vnp1->data.ptrvalue;
      g2 = (GlobalDiscrepancyPtr) vnp2->data.ptrvalue;
      if (g1 != NULL && g2 != NULL && g1->str != NULL && g2->str != NULL) {
        return StringICmp (g1->str, g2->str);
      }
    }
  }
  return 0;
}


static Int4 CountDupGlobalDiscrepancy (ValNodePtr vnp)
{
  GlobalDiscrepancyPtr g1, g2;
  Int4                 num_dup = 1;

  if (vnp == NULL 
      || (g1 = (GlobalDiscrepancyPtr) vnp->data.ptrvalue) == NULL
      || StringHasNoText (g1->str)) {
    return 0;
  } else if (vnp->next == NULL) {
    return 1;
  }
  vnp = vnp->next;
  while (vnp != NULL
         && (g2 = (GlobalDiscrepancyPtr) vnp->data.ptrvalue) != NULL 
         && StringICmp (g1->str, g2->str) == 0) {
    num_dup++;
    vnp = vnp->next;
  }
  return num_dup;
}


extern ClickableItemPtr
ReportNonUniqueGlobalDiscrepancy 
(ValNodePtr vnp, 
 CharPtr    label_fmt,
 CharPtr    ind_cat_fmt,
 Uint4      clickable_item_type,
 Boolean    keep_top_category)

{
  Int4          num_dup, total_dup = 0, i;
  ValNodePtr       item_list;
  ClickableItemPtr cip = NULL;
  ValNodePtr       subcategories = NULL;
  CharPtr          str;

  while (vnp != NULL) {
    num_dup = CountDupGlobalDiscrepancy (vnp);
    if (num_dup > 1) {
      total_dup += num_dup;
      str = GetGlobalDiscrepancyStr (vnp->data.ptrvalue);
      if (str == NULL) str = "";
      item_list = NULL;
      i = num_dup;
      while (i > 0) {
        ValNodeLink (&item_list, GetGlobalDiscrepancyItem (vnp->data.ptrvalue));
        vnp = vnp->next;
        i--;
      }
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (ind_cat_fmt) + StringLen (str) + 15));
      sprintf (cip->description, ind_cat_fmt, num_dup, str);
      cip->clickable_item_type = clickable_item_type;
      cip->item_list = item_list;
      ValNodeAddPointer (&subcategories, 0, cip);
    } else {
      vnp = vnp->next;
    }
  }
  if (subcategories != NULL) {
    if (subcategories->next == NULL && !keep_top_category) {
      cip = subcategories->data.ptrvalue;
      subcategories = ValNodeFree (subcategories);
    } else {
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (label_fmt) + 15));
      sprintf (cip->description, label_fmt, total_dup);
      cip->clickable_item_type = clickable_item_type;
      cip->subcategories = subcategories;
    }
  }
  return cip;
}


static Boolean IsLocusTagFormatBad (CharPtr locus_tag)
{
  CharPtr cp;
  Boolean after_underscore = FALSE;
  
  if (StringHasNoText (locus_tag))
  {
    return FALSE;
  }
  
  cp = locus_tag;
  if (!isalpha (*cp))
  {
    return TRUE;
  }
  cp++;
  while (*cp != 0)
  {
    if (*cp == '_')
    {
      if (after_underscore)
      {
        return TRUE;
      }
      else
      {
        after_underscore = TRUE;
        if (*(cp + 1) == 0)
        {
          return TRUE;
        }
      }
    }
    else if (!isalpha (*cp) && !isdigit (*cp))
    {
      return TRUE;
    }
    cp++;
  }
  if (after_underscore)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


extern ClickableItemPtr ReportBadLocusTagFormat (ValNodePtr list)
{
  ValNodePtr vnp, item_list = NULL;
  ClickableItemPtr cip = NULL;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    if (IsLocusTagFormatBad (GetGlobalDiscrepancyStr (vnp->data.ptrvalue))) {
      ValNodeLink (&item_list, GetGlobalDiscrepancyItem (vnp->data.ptrvalue));
    }
  }
  if (item_list != NULL) {
    cip = NewClickableItem (DISC_GENE_LOCUS_TAG_BAD_FORMAT, "%d locus tags are incorrectly formatted.", item_list);
  }
  return cip;
}


static CharPtr GetGlobalDiscrepancyPrefix (GlobalDiscrepancyPtr g)
{
  CharPtr cp, prefix = NULL;
  Int4    len;

  if (g == NULL) return NULL;
  cp = StringChr (g->str, '_');
  if (cp != NULL) {
    len = cp - g->str;
    prefix = MemNew (sizeof (Char) * (len + 1));
    StringNCpy (prefix, g->str, len);
    prefix[len] = 0;
  }
  return prefix;
}


static Int4 CountDupGlobalDiscrepancyPrefix (ValNodePtr vnp)
{
  GlobalDiscrepancyPtr g1, g2;
  CharPtr              cp;
  Int4                 len;
  Int4                 num_dup = 1;

  if (vnp == NULL 
      || (g1 = (GlobalDiscrepancyPtr) vnp->data.ptrvalue) == NULL
      || StringHasNoText (g1->str)
      || (cp = StringChr (g1->str, '_')) == NULL) {
    return 0;
  } else if (vnp->next == NULL) {
    return 1;
  }
  len = cp - g1->str + 1;
  vnp = vnp->next;
  while (vnp != NULL
         && (g2 = (GlobalDiscrepancyPtr) vnp->data.ptrvalue) != NULL 
         && StringNCmp (g1->str, g2->str, len) == 0) {
    num_dup++;
    vnp = vnp->next;
  }
  return num_dup;
}


extern ValNodePtr ReportInconsistentGlobalDiscrepancyPrefixes
(ValNodePtr vnp, 
 CharPtr    label_fmt,
 Uint4      clickable_item_type)

{
  Int4          num_dup;
  CharPtr       prefix;
  ValNodePtr    disc_list = NULL;
  ClickableItemPtr cip;

  if (vnp == NULL) return NULL;

  num_dup = CountDupGlobalDiscrepancyPrefix (vnp);
  if (num_dup < ValNodeLen (vnp)) {
    while (vnp != NULL) {
      prefix = GetGlobalDiscrepancyPrefix (vnp->data.ptrvalue);
      num_dup = CountDupGlobalDiscrepancyPrefix (vnp);
      if (num_dup < 1) {
        vnp = vnp->next;
      } else if (prefix != NULL) {
        cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        cip->clickable_item_type = clickable_item_type;
        cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (label_fmt) + StringLen (prefix) + 15));
        sprintf (cip->description, label_fmt, num_dup, prefix);
        /* skip duplicates without printing */
        while (num_dup > 0) {
          ValNodeLink (&cip->item_list, GetGlobalDiscrepancyItem (vnp->data.ptrvalue));
          vnp = vnp->next;
          num_dup--;
        }
        prefix = MemFree (prefix);
        ValNodeAddPointer (&disc_list, 0, cip);
      } else {
        /* skip items without prefix */
        while (num_dup > 0) {
          vnp = vnp->next;
          num_dup--;
        }
      }
    }
  } 
  return disc_list;
}


extern ValNodePtr ReportInconsistentGlobalDiscrepancyStrings
(ValNodePtr vnp, 
 CharPtr    label_fmt,
 Uint4      clickable_item_type)

{
  Int4          num_dup;
  CharPtr       prefix;
  ValNodePtr    disc_list = NULL;
  ClickableItemPtr cip;

  if (vnp == NULL) return NULL;

  num_dup = CountDupGlobalDiscrepancy (vnp);
  if (num_dup < ValNodeLen (vnp)) {
    while (vnp != NULL) {
      prefix = GetGlobalDiscrepancyStr (vnp->data.ptrvalue);
      num_dup = CountDupGlobalDiscrepancy (vnp);
      if (num_dup < 1) {
        vnp = vnp->next;
      } else if (prefix != NULL) {
        cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        cip->clickable_item_type = clickable_item_type;
        cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (label_fmt) + StringLen (prefix) + 15));
        sprintf (cip->description, label_fmt, num_dup, prefix);
        /* skip duplicates without printing */
        while (num_dup > 0) {
          ValNodeLink (&cip->item_list, GetGlobalDiscrepancyItem (vnp->data.ptrvalue));
          vnp = vnp->next;
          num_dup--;
        }
        ValNodeAddPointer (&disc_list, 0, cip);
      } else {
        /* skip items without prefix */
        while (num_dup > 0) {
          vnp = vnp->next;
          num_dup--;
        }
      }
    }
  } 
  return disc_list;
}


extern ClickableItemPtr ReportMissingFields (ValNodePtr list, CharPtr label_fmt, Uint4 clickable_item_type)
{
  ClickableItemPtr cip;

  if (list == NULL) return NULL;

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (label_fmt) + 15));
  sprintf (cip->description, label_fmt, ValNodeLen (list));
  cip->clickable_item_type = clickable_item_type;
  while (list != NULL) {
    ValNodeLink (&(cip->item_list), GetGlobalDiscrepancyItem (list->data.ptrvalue));
    list = list->next;
  }
  return cip;  
}


extern Boolean GeneRefMatch (GeneRefPtr grp1, GeneRefPtr grp2)
{
  if (grp1 == NULL && grp2 == NULL)
  {
    return TRUE;
  }
  else if (grp1 == NULL || grp2 == NULL)
  {
    return FALSE;
  }
  else if (StringCmp (grp1->locus, grp2->locus) != 0
           || StringCmp (grp1->allele, grp2->allele) != 0
           || StringCmp (grp1->desc, grp2->desc) != 0
           || StringCmp (grp1->maploc, grp2->maploc) != 0
           || StringCmp (grp1->locus_tag, grp2->locus_tag) != 0
           || (grp1->pseudo && !grp2->pseudo)
           || (!grp1->pseudo && grp2->pseudo)
           || !ValNodeStringListMatch (grp1->db, grp2->db)
           || !ValNodeStringListMatch (grp1->syn, grp2->syn))
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

extern Boolean DbxrefsMatch (ValNodePtr vnp1, ValNodePtr vnp2, Boolean case_sensitive)
{
  Boolean rval = TRUE;

  while (rval && vnp1 != NULL && vnp2 != NULL) {
    if (DbtagMatchEx (vnp1->data.ptrvalue, vnp2->data.ptrvalue, case_sensitive)) {
      vnp1 = vnp1->next;
      vnp2 = vnp2->next;
    } else {
      rval = FALSE;
    }
  }
  if (vnp1 != NULL || vnp2 != NULL) {
    rval = FALSE;
  }
  return rval;
}


extern Boolean XrefsMatch (SeqFeatXrefPtr x1, SeqFeatXrefPtr x2)
{
  Boolean rval = TRUE;

  while (rval && x1 != NULL && x2 != NULL) {
    rval = AsnIoMemComp (x1, x2, (AsnWriteFunc) SeqFeatXrefAsnWrite);
    x1 = x1->next;
    x2 = x2->next;
  }
  if (x1 != NULL || x2 != NULL) {
    rval = FALSE;
  }
  return rval;
}


extern Boolean ProtRefMatch (ProtRefPtr prp1, ProtRefPtr prp2)
{
  if (prp1 == NULL && prp2 == NULL) {
    return TRUE;
  } else if (prp1 == NULL || prp2 == NULL) {
    return FALSE;
  } else if (!ValNodeStringListMatch (prp1->name, prp2->name)
             || StringCmp (prp1->desc, prp2->desc) != 0
             || !ValNodeStringListMatch (prp1->ec, prp2->ec)
             || !ValNodeStringListMatch (prp1->activity, prp2->activity)
             || !DbxrefsMatch (prp1->db, prp2->db, TRUE)
             || prp1->processed != prp2->processed) {
    return FALSE;
  } else {
    return TRUE;
  }
}


/* declarations for discrepancy tests */
extern void AddMissingAndSuperfluousGeneDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void AddDiscrepanciesForNonGeneLocusTags (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindMissingProteinIDs (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindCDSmRNAGeneLocationDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindCDSGeneProductConflicts (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindDuplicateGeneLocus (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void AddECNumberNoteDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindPseudoDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void AddJoinedFeatureDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void AddOverlappingGeneDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void AddContainedCodingRegionDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void AddRNACDSOverlapDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindShortContigs (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindNonmatchingContigSources (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindSuspectProductNames (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindSuspectPhrases (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindInconsistentSourceAndDefline (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindParticalCDSsInCompleteSequences (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindUnknownProteinsWithECNumbers (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindShortSequences (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void tRNACountFeaturesAndFindDups (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void tRNAFindBadLength (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void rRNACountFeaturesAndFindDups (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindRNAsWithoutProducts (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindTranslExceptNotes (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindCDSOverlappingtRNAs (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern int LIBCALLBACK SortVnpByClickableItemDescription (VoidPtr ptr1, VoidPtr ptr2);
extern void CountProteins (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindFeaturesOverlappingSrcFeatures (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
static void PercentNDiscrepanciesForSeqEntry (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
static void FindAdjacentPseudoGenes (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);

/* J. Chen */
static void ShowTranslExcept(ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
static void ShowCDsHavingGene(ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
static void TestDeflineExistence(ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list); 
static void TestMrnaOverlappingPseudo(ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list); 
static void RmvMrnaOverlappingPseudoGene(ValNodePtr item_list, Pointer data, LogInfoPtr lip);
/* J. Chen */


static void RmvMrnaOverlappingPseudoGene(ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr vnp, entityIDList = NULL;
  SeqFeatPtr sfp;
  SeqMgrFeatContext context;
  CharPtr           txt;

/*
  if (Message (MSG_OKC, "Are you sure you want to remove mRNA overlapping a pseudogene?") == ANS_CANCEL) {
    return;
  }
*/

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT) {
      sfp = (BioseqPtr) vnp->data.ptrvalue;
      sfp->idx.deleteme = TRUE;
      ValNodeAddInt (&entityIDList, 0, sfp->idx.entityID);
    }
  }

  for (vnp = entityIDList; vnp != NULL; vnp = vnp->next) {
    DeleteMarkedObjects (vnp->data.intvalue, 0, NULL);
    ObjMgrSetDirtyFlag (vnp->data.intvalue, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, vnp->data.intvalue, 0, 0);
  }

  entityIDList = ValNodeFree (entityIDList);

}  /* RmvMrnaOverlappingPseudoGene */




static void GetMrnaOverlappingPseudoGene(SeqFeatPtr sfp, Pointer userdata)
{
  GeneRefPtr      grp;
  SeqFeatPtr      gene_sfp = NULL;
  RnaRefPtr	  rna_rp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || userdata == NULL)
  {
    return;
  }
  rna_rp = (RnaRefPtr)sfp->data.value.ptrvalue;
  if (rna_rp->type != 2) return;  /* not a mRNA */
  
  gene_sfp = GetGeneForFeature(sfp);
  if (gene_sfp == NULL)
  {
    return;
  }

  if (gene_sfp->pseudo)
  {
    ValNodeAddPointer (userdata, OBJ_SEQFEAT, sfp);
  }
}  /* GetMrnaOverlappingPseudoGene */



static void TestMrnaOverlappingPseudoGene(ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr         pseudo_features = NULL, vnp;
  ClickableItemPtr   dip;
  CharPtr            bad_fmt = "%d Pseudogenes have overlapping mRNAs.";

  if (discrepancy_list == NULL) return;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &pseudo_features, GetMrnaOverlappingPseudoGene);
  }

  if (pseudo_features != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = TEST_MRNA_OVERLAPPING_PSEUDO_GENE;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (pseudo_features));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = pseudo_features;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
} /* TestMrnaOverlappingPseudoGene */



static void HasDefline(BioseqPtr bsp, Pointer userdata)
{
  SeqDescrPtr sdp;
  SeqMgrDescContext context;

  if (bsp == NULL || ISA_aa (bsp->mol) || userdata == NULL) {
    return;
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &context);
  if (sdp != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_BIOSEQ, bsp);
  }

}  /* HasDefline */




static void FindOneDefline (SeqEntryPtr sep, Pointer item_list)
{
   BioseqSetPtr bssp;
   SeqEntryPtr  tmp;

   if (IS_Bioseq(sep)) VisitBioseqsInSep (sep, item_list, HasDefline);
   else if (IS_Bioseq_set(sep)) {
            bssp = (BioseqSetPtr) sep->data.ptrvalue;
            for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
                  FindOneDefline(tmp, item_list);
                  if (item_list != NULL) break;
            }
   }
} /* FindOneDefline */





static void TestDeflineExistence(ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr  vnp, item_list = NULL;
  SeqEntryPtr sep;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = (SeqEntryPtr) vnp->data.ptrvalue;
    FindOneDefline(sep, &item_list);
    if (item_list != NULL) {
        ValNodeAddPointer (discrepancy_list, 0, 
	     NewClickableItem (TEST_DEFLINE_PRESENT, "Bioseq has definition line", item_list));
        break;
    }
  }
}  /* TestDeflineExistence */




/* J. Chen */
static void FindCDsHavingGeneName(SeqFeatPtr sfp, Pointer userdata)
{
  SeqFeatPtr          gene, protein_feat;
  BioseqPtr           protein_seq;
  SeqMgrFeatContext   fcontext;
  ValNodePtr          prot_nm;
  ProtRefPtr          prp;
  GeneRefPtr          gene_p;

   if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;

  gene = GetGeneForFeature (sfp);  /* one gene */
  if (gene == NULL) { /* no gene means no gene name */
    return;
  }
  gene_p = (GeneRefPtr)gene->data.value.ptrvalue;
  if (gene_p == NULL) return;
  if (gene_p->locus == NULL) return; /* no gene name */

  protein_seq = BioseqFindFromSeqLoc (sfp->product);
  protein_feat = SeqMgrGetNextFeature (protein_seq, 0, SEQFEAT_PROT, FEATDEF_PROT, &fcontext);
  if (protein_feat == NULL
                || (prp = (ProtRefPtr) protein_feat->data.value.ptrvalue) == NULL) return;

  for (prot_nm = prp->name;  prot_nm !=  NULL; prot_nm = prot_nm->next) {
    if (strstr(prot_nm->data.ptrvalue, "hypothetical protein") != NULL) {
        ValNodeAddPointer (userdata, OBJ_SEQFEAT, sfp);
    }
  }
}    /* FindCDsHavingGeneName */

   
  

  
/* Display hypothetic protein having a gene name: J. Chen */
static void ShowCDsHavingGene(ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{ 
   ValNodePtr       vnp, item_list;
   SeqEntryPtr      sep;
   ClickableItemPtr cip;
   CharPtr          show_CDs = "%d hypothetical coding regions have a gene name";

   for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
     sep = vnp->data.ptrvalue;
     item_list = NULL;
     VisitFeaturesInSep (sep, &item_list, FindCDsHavingGeneName);
     if (item_list != NULL) {
       cip = NewClickableItem (SHOW_HYPOTHETICAL_CDS_HAVING_GENE_NAME, show_CDs, item_list);
       ValNodeAddPointer (discrepancy_list, 0, cip);
     }
   }
}   /* ShowCDsHavingGene() */


/* Find code breaks in a coding region: J. Chen */
static Boolean CodingRegionHasCodeBreak(SeqFeatPtr sfp)
{
  CdRegionPtr  crp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION
      || (crp = (CdRegionPtr)sfp->data.value.ptrvalue) == NULL
      || crp->code_break == NULL)
  {
      return FALSE;
  }
  else return TRUE;
}

/* Find all coding regions that have a translation exception: J. Chen */
static void CodingRegionsHaveTranslExcept(SeqFeatPtr sfp, Pointer userdata)
{
  if (sfp != NULL && userdata != NULL && sfp->data.choice == SEQFEAT_CDREGION) {

    if (CodingRegionHasCodeBreak(sfp)) {
        ValNodeAddPointer (userdata, OBJ_SEQFEAT, sfp);
    } 
  }

} /* CodingRegionsHaveTranslExcept() */



/* Show the translation exceptions: J. Chen */
static void ShowTranslExcept(ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{ 
  ValNodePtr       vnp, item_list;
  SeqEntryPtr      sep;
  ClickableItemPtr cip;
  CharPtr          show_transl_except = "%d coding regions have a translation exception";

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    item_list = NULL;
    VisitFeaturesInSep (sep, &item_list, CodingRegionsHaveTranslExcept);
    if (item_list != NULL) {
      cip = NewClickableItem (SHOW_TRANSL_EXCEPT, show_transl_except, item_list);
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
  }
}  /* ShowTranslExceptInCDs */


/* functions for the missing and superfluous gene tests */
static Boolean GeneRefMatchForSuperfluousCheck (GeneRefPtr grp1, GeneRefPtr grp2)
{
  if (grp1 == NULL && grp2 == NULL)
  {
    return TRUE;
  }
  else if (grp1 == NULL || grp2 == NULL)
  {
    return FALSE;
  }
  else if (StringCmp (grp1->locus, grp2->locus) != 0
           || StringCmp (grp1->locus_tag, grp2->locus_tag) != 0
           || (grp1->pseudo && !grp2->pseudo)
           || (!grp1->pseudo && grp2->pseudo))
  {
    return FALSE;
  }
  else if (StringHasNoText (grp1->allele) 
           && !StringHasNoText (grp2->allele) 
           && StringCmp (grp1->allele, grp2->allele) != 0) 
  {
    return FALSE;
  }
  else if (StringHasNoText (grp1->desc) 
           && !StringHasNoText (grp2->desc) 
           && StringCmp (grp1->desc, grp2->desc) != 0) 
  {
    return FALSE;
  }
  else if (StringHasNoText (grp1->maploc) 
           && !StringHasNoText (grp2->maploc) 
           && StringCmp (grp1->maploc, grp2->maploc) != 0) 
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


static void ExtractGeneFromListByGeneRef (ValNodePtr PNTR list, GeneRefPtr grp)
{
  ValNodePtr prev = NULL, this_vnp, next_vnp;
  SeqFeatPtr gene_feat;
  
  if (list == NULL || grp == NULL)
  {
    return;
  }
  
  this_vnp = *list;
  while (this_vnp != NULL)
  {
    next_vnp = this_vnp->next;
    gene_feat = (SeqFeatPtr) this_vnp->data.ptrvalue;
    if (gene_feat != NULL && GeneRefMatchForSuperfluousCheck (gene_feat->data.value.ptrvalue, grp))
    {
      if (prev == NULL)
      {
        *list = next_vnp;
      }
      else
      {
        prev->next = next_vnp;
      }
      this_vnp->next = NULL;
      ValNodeFree (this_vnp);
    }
    else
    {
      prev = this_vnp;
    }
    this_vnp = next_vnp;    
  }
}


static void ExtractGeneFromListByGene (ValNodePtr PNTR list, SeqFeatPtr gene)
{
  ValNodePtr prev = NULL, this_vnp, next_vnp;
  
  if (list == NULL || gene == NULL)
  {
    return;
  }
  
  this_vnp = *list;
  while (this_vnp != NULL)
  {
    next_vnp = this_vnp->next;
    if (this_vnp->data.ptrvalue == gene)
    {
      if (prev == NULL)
      {
        *list = next_vnp;
      }
      else
      {
        prev->next = next_vnp;
      }
      this_vnp->next = NULL;
      ValNodeFree (this_vnp);
    }
    else
    {
      prev = this_vnp;
    }
    this_vnp = next_vnp;    
  }
}


static void 
CheckGenesForFeatureType 
(ValNodePtr PNTR features_without_genes,
 ValNodePtr PNTR superfluous_genes,
 BioseqPtr  bsp,
 Uint1      feature_type,
 Uint1      feature_subtype,
 Boolean    makes_gene_not_superfluous)
{
  SeqFeatPtr         sfp, gene_sfp;
  GeneRefPtr         grp;
  SeqMgrFeatContext  context;
  
  if (features_without_genes == NULL
      || superfluous_genes == NULL
      || bsp == NULL)
  {
    return;
  }
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, feature_type, feature_subtype, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, feature_type, feature_subtype, &context))
  {
    if (sfp->data.choice == SEQFEAT_GENE) {
      continue;
    }
    /* check for gene xref */
    grp = SeqMgrGetGeneXref (sfp);
    if (grp != NULL)
    {
      if (SeqMgrGeneIsSuppressed (grp))
      {
        ValNodeAddPointer (features_without_genes, OBJ_SEQFEAT, sfp);
      }
      else
      {
        ExtractGeneFromListByGeneRef (superfluous_genes, grp);
      }
    }
    else
    {
      gene_sfp = SeqMgrGetOverlappingGene (sfp->location, NULL);
      if (gene_sfp == NULL)
      {
        ValNodeAddPointer (features_without_genes, OBJ_SEQFEAT, sfp);
      }
      else if (makes_gene_not_superfluous)
      {
        ExtractGeneFromListByGene (superfluous_genes, gene_sfp);
      }
    }  
  }  
}

typedef struct misssupergenes
{
  ValNodePtr missing_list;
  ValNodePtr super_list;
  Boolean    any_genes;
} MissSuperGenesData, PNTR MissSuperGenesPtr;


static void FindMissingGenes (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr         sfp;
  SeqMgrFeatContext  context;
  ValNodePtr         features_without_genes = NULL;
  ValNodePtr         superfluous_genes = NULL;
  MissSuperGenesPtr  msgp;
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  msgp = (MissSuperGenesPtr) userdata;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, FEATDEF_GENE, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_GENE, FEATDEF_GENE, &context))
  {
    ValNodeAddPointer (&superfluous_genes, OBJ_SEQFEAT, sfp);
  }
  
  CheckGenesForFeatureType (&features_without_genes, &superfluous_genes, bsp, 
                            SEQFEAT_CDREGION, 0, TRUE);
  CheckGenesForFeatureType (&features_without_genes, &superfluous_genes, bsp, 
                            SEQFEAT_RNA, 0, TRUE);
  CheckGenesForFeatureType (&features_without_genes, &superfluous_genes, bsp, 
                            SEQFEAT_IMP, FEATDEF_RBS, FALSE);
  CheckGenesForFeatureType (&features_without_genes, &superfluous_genes, bsp, 
                            SEQFEAT_IMP, FEATDEF_exon, FALSE);
  CheckGenesForFeatureType (&features_without_genes, &superfluous_genes, bsp, 
                            SEQFEAT_IMP, FEATDEF_intron, FALSE);

  ValNodeLink (&(msgp->missing_list), features_without_genes);
  if (IsmRNASequenceInGenProdSet (bsp)) {
    superfluous_genes = ValNodeFree (superfluous_genes);
  } else {
    ValNodeLink (&(msgp->super_list), superfluous_genes);
  }
}


static void 
GetPseudoAndNonPseudoGeneList 
(ValNodePtr      super_list,
 ValNodePtr PNTR pseudo_list, 
 ValNodePtr PNTR non_pseudo_list)
{
  ValNodePtr vnp;
  SeqFeatPtr gene;
  GeneRefPtr grp;
  
  if (pseudo_list == NULL || non_pseudo_list == NULL)
  {
    return;
  }
  *pseudo_list = NULL;
  *non_pseudo_list = NULL;
  
  for (vnp = super_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == OBJ_SEQFEAT)
    {
      gene = (SeqFeatPtr) vnp->data.ptrvalue;
      if (gene != NULL && gene->data.choice == SEQFEAT_GENE)
      {
        grp = (GeneRefPtr) gene->data.value.ptrvalue;
        if (gene->pseudo || (grp != NULL && grp->pseudo))
        {
          ValNodeAddPointer (pseudo_list, OBJ_SEQFEAT, gene);
        }
        else
        {
          ValNodeAddPointer (non_pseudo_list, OBJ_SEQFEAT, gene);
        }
      }
    }
  }
}


static void 
GetFrameshiftAndNonFrameshiftGeneList 
(ValNodePtr      super_list,
 ValNodePtr PNTR frameshift_list, 
 ValNodePtr PNTR non_frameshift_list)
{
  ValNodePtr vnp;
  SeqFeatPtr gene;
  
  if (frameshift_list == NULL || non_frameshift_list == NULL)
  {
    return;
  }
  *frameshift_list = NULL;
  *non_frameshift_list = NULL;
  
  for (vnp = super_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == OBJ_SEQFEAT)
    {
      gene = (SeqFeatPtr) vnp->data.ptrvalue;
      if (gene != NULL 
          && (StringISearch (gene->comment, "frameshift") != NULL
              || StringISearch (gene->comment, "frame shift") != NULL)) 
      {
        ValNodeAddPointer (frameshift_list, OBJ_SEQFEAT, gene);
      }
      else
      {
        ValNodeAddPointer (non_frameshift_list, OBJ_SEQFEAT, gene);
      }
    }
  }
}


extern void AddMissingAndSuperfluousGeneDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip, pseudo_dip, non_pseudo_dip;
  CharPtr            missing_genes_fmt = "%d features have no genes.";
  CharPtr            extra_genes_fmt = "%d gene features are not associated with a CDS or RNA feature.";
  CharPtr            pseudo_extra_genes_fmt = "%d pseudo gene features are not associated with a CDS or RNA feature.";
  CharPtr            non_pseudo_frameshift_extra_genes_fmt = "%d non-pseudo gene features are not associated with a CDS or RNA feature and have frameshift in the comment.";
  CharPtr            non_pseudo_non_frameshift_extra_genes_fmt = "%d non-pseudo gene features are not associated with a CDS or RNA feature and do not have frameshift in the comment.";
  MissSuperGenesData msgd;
  ValNodePtr         non_pseudo_list = NULL, pseudo_list = NULL, vnp;
  ValNodePtr         non_frameshift_list = NULL, frameshift_list = NULL;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  msgd.missing_list = NULL;
  msgd.super_list = NULL;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &msgd, FindMissingGenes);
  }
  
  if (msgd.missing_list != NULL)
  {
    dip = NewClickableItem (DISC_GENE_MISSING, missing_genes_fmt, msgd.missing_list);
    if (dip != NULL)
    {
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
  
  if (msgd.super_list != NULL)
  {
    GetPseudoAndNonPseudoGeneList (msgd.super_list, &pseudo_list, &non_pseudo_list);
    GetFrameshiftAndNonFrameshiftGeneList (non_pseudo_list, &frameshift_list, &non_frameshift_list);
    non_pseudo_list = ValNodeFree (non_pseudo_list);
    dip = NewClickableItem (DISC_SUPERFLUOUS_GENE, extra_genes_fmt, msgd.super_list);
    if (dip != NULL)
    {
      ValNodeAddPointer (discrepancy_list, 0, dip);     
 
      if (frameshift_list != NULL) 
      {
        non_pseudo_dip = NewClickableItem (DISC_SUPERFLUOUS_GENE, non_pseudo_frameshift_extra_genes_fmt, frameshift_list);
        non_pseudo_dip->level = 1;
        ValNodeAddPointer (&(dip->subcategories), 0, non_pseudo_dip);
      }
      if (non_frameshift_list != NULL) 
      {
        non_pseudo_dip = NewClickableItem (DISC_SUPERFLUOUS_GENE, non_pseudo_non_frameshift_extra_genes_fmt, non_frameshift_list);
        non_pseudo_dip->level = 1;
        ValNodeAddPointer (&(dip->subcategories), 0, non_pseudo_dip);
      }
      if (pseudo_list != NULL) 
      {
        pseudo_dip = NewClickableItem (DISC_SUPERFLUOUS_GENE, pseudo_extra_genes_fmt, pseudo_list);
        pseudo_dip->level = 1;
        ValNodeAddPointer (&(dip->subcategories), 0, pseudo_dip);
      }
    }
  }  
}


static Boolean CommentHasPhrase (CharPtr comment, CharPtr phrase)
{
  CharPtr cp;
  Int4    len;

  if (StringHasNoText (comment) || StringHasNoText (phrase)) {
    return FALSE;
  }
  len = StringLen (phrase);
  cp = comment;
  while (cp != NULL) {
    if (StringNICmp (comment, phrase, len) == 0 && (*(cp + len) == ';' || *(cp + len) == 0)) {
      return TRUE;
    } else {
      cp = StringChr (cp, ';');
      if (cp != NULL) {
        cp++;
        cp += StringSpn (cp, " ");
      }
    }
  }
  return FALSE;
}


static Boolean IsOkSuperfluousGene (SeqFeatPtr sfp)
{
  GeneRefPtr grp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE) {
    return FALSE;
  } else if (CommentHasPhrase (sfp->comment, "coding region not determined")) {
    return TRUE;
  } else if ((grp = sfp->data.value.ptrvalue) == NULL) {
    return FALSE;
  } else if (CommentHasPhrase (grp->desc, "coding region not determined")) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void FindOncallerMissingGenes (BioseqPtr bsp, Pointer data)
{
  MissSuperGenesPtr msgp;
  ValNodePtr        features_without_genes = NULL, superfluous_genes = NULL;
  ValNodePtr        other = NULL, prev = NULL, vnp, vnp_next;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext context;

  if (bsp == NULL || ISA_aa (bsp->mol) || (msgp = (MissSuperGenesPtr) data) == NULL) {
    return;
  }

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, FEATDEF_GENE, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_GENE, FEATDEF_GENE, &context))
  {
    msgp->any_genes = TRUE;
    if (sfp->pseudo) {
      continue;
    }

    ValNodeAddPointer (&superfluous_genes, OBJ_SEQFEAT, sfp);
  }
  
  /* look for features without genes that we care about */
  CheckGenesForFeatureType (&features_without_genes, &superfluous_genes, bsp, 
                            SEQFEAT_CDREGION, 0, TRUE);
  CheckGenesForFeatureType (&features_without_genes, &superfluous_genes, bsp, 
                            SEQFEAT_RNA, FEATDEF_mRNA, TRUE);
  CheckGenesForFeatureType (&features_without_genes, &superfluous_genes, bsp, 
                            SEQFEAT_RNA, FEATDEF_tRNA, TRUE);

  /* all other feature types make genes not superfluous */
  CheckGenesForFeatureType (&other, &superfluous_genes, bsp, 
                            0, 0, TRUE);

  other = ValNodeFree (other);

  ValNodeLink (&(msgp->missing_list), features_without_genes);
  if (IsmRNASequenceInGenProdSet (bsp)) {
    superfluous_genes = ValNodeFree (superfluous_genes);
  } else {
    /* remove genes with explanatory comments/descriptions */
    for (vnp = superfluous_genes; vnp != NULL; vnp = vnp_next) {
      vnp_next = vnp->next;
      if (IsOkSuperfluousGene((SeqFeatPtr)vnp->data.ptrvalue)) {
        if (prev == NULL) {
          superfluous_genes = vnp->next;
        } else {
          prev->next = vnp->next;
        }
        vnp->next = NULL;
        vnp = ValNodeFree (vnp);
      } else {
        prev = vnp;
      }
    }

    ValNodeLink (&(msgp->super_list), superfluous_genes);
  }

}


static void OnCallerMissingAndSuperfluousGenes (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr   dip;
  CharPtr            missing_genes_fmt = "%d features have no genes.";
  CharPtr            extra_genes_fmt = "%d gene features are not associated with any feature and are not pseudo.";
  MissSuperGenesData msgd;
  ValNodePtr         vnp;
  SeqEntryPtr        oldscope;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  msgd.missing_list = NULL;
  msgd.super_list = NULL;
  msgd.any_genes = FALSE;
  
  oldscope = SeqEntryGetScope ();
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    SeqEntrySetScope (vnp->data.ptrvalue);
    VisitBioseqsInSep (vnp->data.ptrvalue, &msgd, FindOncallerMissingGenes);
  }
  SeqEntrySetScope (oldscope);
  
  if (msgd.any_genes) {
    if (msgd.missing_list != NULL)
    {
      dip = NewClickableItem (ONCALLER_GENE_MISSING, missing_genes_fmt, msgd.missing_list);
      if (dip != NULL)
      {
        ValNodeAddPointer (discrepancy_list, 0, dip);
      }
    }
    
    if (msgd.super_list != NULL)
    {
      dip = NewClickableItem (ONCALLER_SUPERFLUOUS_GENE, extra_genes_fmt, msgd.super_list);
      if (dip != NULL)
      {
        ValNodeAddPointer (discrepancy_list, 0, dip);      
      }
    }
  } else {
    msgd.missing_list = ValNodeFree (msgd.missing_list);
    msgd.super_list = ValNodeFree (msgd.super_list);
  }
}


/* test for missing or inconsistent protein IDs */
typedef struct prefixcheck 
{
  CharPtr prefix;
  ValNodePtr feature_list;
} PrefixCheckData, PNTR PrefixCheckPtr;


static ValNodePtr FreePrefixCheckList (ValNodePtr prefix_list)
{
  PrefixCheckPtr pcp;
  
  if (prefix_list == NULL)
  {
    return NULL;
  }
  
  prefix_list->next = FreePrefixCheckList (prefix_list->next);
  
  pcp = (PrefixCheckPtr) prefix_list->data.ptrvalue;
  if (pcp != NULL)
  {
    pcp->prefix = MemFree (pcp->prefix);
    pcp->feature_list = ValNodeFree (pcp->feature_list);
    pcp = MemFree (pcp);
  }
  prefix_list = ValNodeFree (prefix_list);
  return NULL;
}


static ClickableItemPtr InconsistentPrefix (PrefixCheckPtr pcp, CharPtr bad_fmt, DiscrepancyType disc_type)
{
  ClickableItemPtr dip = NULL;

  if (pcp == NULL || StringHasNoText (pcp->prefix) || pcp->feature_list == NULL)
  {
    return NULL;
  }
  dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  if (dip != NULL)
  {
    dip->clickable_item_type = disc_type;
    dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + StringLen (pcp->prefix)+ 15));
    sprintf (dip->description, bad_fmt, ValNodeLen (pcp->feature_list), pcp->prefix);
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    dip->item_list = pcp->feature_list;
    pcp->feature_list = NULL;
  }      
  return dip;
}


CharPtr discReportInconsistentLocusTagPrefixFmt = "%d features have locus tag prefix %s.";
CharPtr discReportInconsistentProteinIDPrefixFmt = "%d sequences have protein ID prefix %s.";
CharPtr discReportBadProteinIdFmt = "%d proteins have invalid IDs.";

static ClickableItemPtr InconsistentLocusTagPrefix (PrefixCheckPtr pcp)
{
  return InconsistentPrefix (pcp, discReportInconsistentLocusTagPrefixFmt, DISC_GENE_LOCUS_TAG_INCONSISTENT_PREFIX);
}


extern void FindProteinIDCallback (BioseqPtr bsp, Pointer userdata)
{
  ProtIdListsPtr pip;
  SeqIdPtr       sip;
  DbtagPtr       dbt = NULL;

  if (bsp == NULL || ! ISA_aa (bsp->mol) || userdata == NULL)
  {
    return;
  }

  pip = (ProtIdListsPtr) userdata;

  for (sip = bsp->id; sip != NULL && dbt == NULL; sip = sip->next)
  {
    if (sip->choice == SEQID_GENERAL)
    {
      dbt = (DbtagPtr) sip->data.ptrvalue;
      if (IsSkippableDbtag(dbt))
      {
        dbt = NULL;
      }
    }
  }
  if (dbt == NULL)
  {
    ValNodeAddPointer (&(pip->missing_gnl_list), 0, GlobalDiscrepancyNew (NULL, OBJ_BIOSEQ, bsp));
  }
  else
  {  
    ValNodeAddPointer (&(pip->gnl_list), 0, GlobalDiscrepancyNew (dbt->db, OBJ_BIOSEQ, bsp));
  }
}



extern void FindMissingProteinIDs (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip;
  ProtIdListsData  pid;
  ValNodePtr       vnp;
  
  if (discrepancy_list == NULL) return;
  
  MemSet (&pid, 0, sizeof (ProtIdListsData));
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &pid, FindProteinIDCallback);
  }
  
  if (pid.missing_gnl_list != NULL)
  {
    dip = ReportMissingFields (pid.missing_gnl_list, discReportBadProteinIdFmt, DISC_MISSING_PROTEIN_ID);
    if (dip != NULL) {
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
    pid.missing_gnl_list = FreeGlobalDiscrepancyList (pid.missing_gnl_list);
  }
  if (pid.gnl_list != NULL)
  {
    pid.gnl_list = ValNodeSort (pid.gnl_list, SortVnpByGlobalDiscrepancyString);
    ValNodeLink (discrepancy_list, 
                 ReportInconsistentGlobalDiscrepancyStrings (pid.gnl_list,
                                                             discReportInconsistentProteinIDPrefixFmt,
                                                             DISC_INCONSISTENT_PROTEIN_ID_PREFIX));
    pid.gnl_list = FreeGlobalDiscrepancyList (pid.gnl_list);
  } 
}


typedef struct locustagcheck
{
  ValNodePtr locus_tags_list;
  ValNodePtr missing_list;
  Boolean    exclude_dirsub;
} LocusTagCheckData, PNTR LocusTagCheckPtr;

static void GeneLocusTagDiscrepancyCallback (ValNodePtr item_list, Pointer userdata)
{
  Message (MSG_OK, "I could launch the editor for the individual gene...");
}

static Boolean IsBacterialBioSource (BioSourcePtr biop);

/* Not WGS, genome, or RefSeq */
static Boolean IsLocationDirSub (SeqLocPtr slp)
{
  SeqIdPtr sip;
  Boolean  rval = TRUE, is_complete = FALSE;
  BioseqPtr bsp;
  SeqDescrPtr sdp;
  SeqMgrDescContext context;
  GBBlockPtr gbp;
  MolInfoPtr mip;
  BioSourcePtr biop;
  ValNodePtr vnp;

  if (slp == NULL) {
    rval = FALSE;
  } else {
    sip = SeqLocId (slp);
    if (sip == NULL) {
      rval = FALSE;
    } else if (sip->choice == SEQID_OTHER) {
      rval = FALSE;
    } else {
      bsp = BioseqLockById (sip);
      if (bsp == NULL) {
        rval = TRUE;
      } else {
        rval = TRUE;
        for (sip = bsp->id; sip != NULL && !rval; sip = sip->next) {
          if (sip->choice == SEQID_OTHER) {
            rval = FALSE;
          }
        }
        /* look for WGS keyword */
        for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &context);
              sdp != NULL && rval;
              sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_genbank, &context)) {
          gbp = (GBBlockPtr) sdp->data.ptrvalue;
          for (vnp = gbp->keywords; vnp != NULL && rval; vnp = vnp->next) {
            if (StringICmp ((CharPtr)vnp->data.ptrvalue, "WGS") == 0) {
              rval = FALSE;
            }
          }       
        }
        for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
              sdp != NULL && rval;
              sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_molinfo, &context)) {
          mip = (MolInfoPtr) sdp->data.ptrvalue;
          if (mip->tech == MI_TECH_wgs) {
            rval = FALSE;
          }
          if (mip->completeness == 1) {
            is_complete = TRUE;
          }
        }
        /* is genome? (complete and bacterial)? */
        if (is_complete) {
          for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
              sdp != NULL && rval;
              sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &context)) {
            biop = (BioSourcePtr) sdp->data.ptrvalue;
            if (IsBacterialBioSource(biop)) {
              rval = FALSE;
            }
          }
        }
      }
      BioseqUnlock (bsp);
    }
  }
  return rval;
}


static void CheckGeneLocusTag (SeqFeatPtr sfp, Pointer userdata)
{
  GeneRefPtr         grp;
  LocusTagCheckPtr   ltcp;
  
  if (sfp == NULL || userdata == NULL || sfp->data.choice != SEQFEAT_GENE || sfp->data.value.ptrvalue == NULL)
  {
    return;
  }
  
  ltcp = (LocusTagCheckPtr) userdata;

  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp != NULL) {
    if (grp->pseudo) return;
    if (StringDoesHaveText (grp->locus_tag)) {
      ValNodeAddPointer (&(ltcp->locus_tags_list), 0, 
                          GlobalDiscrepancyNew (grp->locus_tag, OBJ_SEQFEAT, sfp));
    } else {
      if (!ltcp->exclude_dirsub || !IsLocationDirSub (sfp->location)) {
        ValNodeAddPointer (&(ltcp->missing_list), 0,
                            GlobalDiscrepancyNew (NULL, OBJ_SEQFEAT, sfp));
      }
    }
  }
}

static Boolean AlreadyInList (ValNodePtr vnp, SeqFeatPtr sfp)
{
  while (vnp != NULL && vnp->data.ptrvalue != sfp)
  {
    vnp = vnp->next;
  }
  if (vnp == NULL)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


static SeqFeatPtr GetNextGene (SeqFeatPtr sfp)
{
  BioseqPtr bsp;
  SeqFeatPtr sfp_next;
  SeqMgrFeatContext fcontext;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE) return NULL;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return NULL;
  /* initialize fcontext for search */
  sfp_next = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, FEATDEF_GENE, &fcontext);
  while (sfp_next != sfp && sfp_next != NULL)
  {
    sfp_next = SeqMgrGetNextFeature (bsp, sfp_next, SEQFEAT_GENE, FEATDEF_GENE, &fcontext);
  }
  if (sfp_next != sfp) return NULL;

  sfp_next = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_GENE, FEATDEF_GENE, &fcontext);
  return sfp_next;
}


static ValNodePtr FindValNodeForGlobalDiscrepancyFeature (ValNodePtr start_search, Int4 len_search, SeqFeatPtr sfp)
{
  GlobalDiscrepancyPtr g;

  while (start_search != NULL && len_search > 0) {
    g = (GlobalDiscrepancyPtr) start_search->data.ptrvalue;
    if (g != NULL && g->data_choice == OBJ_SEQFEAT && g->data == sfp) {
      return start_search;
    } else {
      start_search = start_search->next;
      len_search--;
    }
  }
  return NULL;
}


static ValNodePtr FindAdjacentGenesInSubList (ValNodePtr sub_list, Int4 list_len)
{
  GlobalDiscrepancyPtr g;
  SeqFeatPtr           sfp, sfp_next;
  ValNodePtr           vnp, found_match, adj_list = NULL;
  Int4                 len;
  
  vnp = sub_list;
  len = list_len;
  while (vnp != NULL && len > 0) {
    g = (GlobalDiscrepancyPtr) vnp->data.ptrvalue;
    if (g->data_choice == OBJ_SEQFEAT && g->data != NULL) {
      sfp = g->data;
      sfp_next = GetNextGene (sfp);
      if (sfp_next != NULL) {
        found_match = FindValNodeForGlobalDiscrepancyFeature (sub_list, list_len, sfp_next);
        if (found_match != NULL) {
          if (vnp->choice == 0) {
            ValNodeAddPointer (&adj_list, OBJ_SEQFEAT, sfp);
            vnp->choice = 1;
          }
          if (found_match->choice == 0) {
            ValNodeAddPointer (&adj_list, OBJ_SEQFEAT, sfp_next);
            found_match->choice = 1;
          }
        }
      }
    }
    vnp = vnp->next;
    len --;
  }
  /* set choices back to zero */
  for (vnp = sub_list, len = 0; vnp != NULL && len < list_len; vnp = vnp->next, len++) {
    vnp->choice = 0;
  }
  return adj_list;
}
 

extern ClickableItemPtr FindAdjacentDuplicateLocusTagGenes (ValNodePtr locus_tag_list)
{
  ValNodePtr       vnp, adjacent_list = NULL;
  ClickableItemPtr cip = NULL;
  Int4             num_dup;
  CharPtr          duplicate_adjacent_fmt = "%d genes are adjacent to another gene with the same locus tag.";

  vnp = locus_tag_list;
  while (vnp != NULL) {
    num_dup = CountDupGlobalDiscrepancy (vnp);
    if (num_dup > 1) {
      ValNodeLink (&adjacent_list, FindAdjacentGenesInSubList (vnp, num_dup));      
      while (num_dup > 0) {
        vnp = vnp->next;
        num_dup--;
      }
    } else {
      vnp = vnp->next;
    }
  }

  if (adjacent_list != NULL) {
    cip = NewClickableItem (DISC_GENE_DUPLICATE_LOCUS_TAG, duplicate_adjacent_fmt, adjacent_list);
  }
  return cip;
}


NLM_EXTERN ValNodePtr ValNodeDupStringList (ValNodePtr vnp)
{
  ValNodePtr cpy = NULL, last = NULL, tmp;

  while (vnp != NULL)
  {
    tmp = ValNodeNew (NULL);
    tmp->choice = vnp->choice;
    tmp->data.ptrvalue = StringSave (vnp->data.ptrvalue);
    if (last == NULL)
    {
      cpy = tmp;
    }
    else
    {
      last->next = tmp;
    }
    last = tmp;
    vnp = vnp->next;
  }
  return cpy;
}
      

NLM_EXTERN ValNodePtr ValNodeDupIntList (ValNodePtr vnp)
{
  ValNodePtr cpy = NULL, last = NULL, tmp;

  while (vnp != NULL)
  {
    tmp = ValNodeNew (NULL);
    tmp->choice = vnp->choice;
    tmp->data.intvalue = vnp->data.intvalue;
    if (last == NULL)
    {
      cpy = tmp;
    }
    else
    {
      last->next = tmp;
    }
    last = tmp;
    vnp = vnp->next;
  }
  return cpy;
}
      

NLM_EXTERN ValNodePtr FindBadLocusTagsInList (ValNodePtr list)
{
  ValNodePtr       bad_list = NULL, list_copy;
  ValNodePtr       vnp;
  
  list_copy = ValNodeDupStringList (list);
  list_copy = ValNodeSort (list_copy, SortVnpByString);

  for (vnp = list_copy; vnp != NULL; vnp = vnp->next) {
    /* look for badly formatted locus tags */
    if (IsLocusTagFormatBad (vnp->data.ptrvalue)) {
      ValNodeAddPointer(&bad_list, eLocusTagErrorBadFormat, StringSave (vnp->data.ptrvalue));
    }
  }
  list_copy = ValNodeFreeData (list_copy);
  return bad_list;  
}


CharPtr discReportDuplicateLocusTagFmt = "%d genes have duplicate locus tags.";
CharPtr discReportOneDuplicateLocusTagFmt = "%d genes have locus tag %s.";
CharPtr discReportMissingLocusTags = "%d genes have no locus tags.";

extern void AddDiscrepanciesForMissingOrNonUniqueGeneLocusTagsEx (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list, Boolean exclude_dirsub)
{
  LocusTagCheckData  ltcd;
  ClickableItemPtr dip = NULL, dip_sub;
  ValNodePtr         vnp;
  
  if (discrepancy_list == NULL) return;
  ltcd.locus_tags_list = NULL;
  ltcd.missing_list = NULL;
  ltcd.exclude_dirsub = exclude_dirsub;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitGenProdSetFeatures (vnp->data.ptrvalue, &ltcd, CheckGeneLocusTag);
  }
  
  if (ltcd.locus_tags_list != NULL) {
    ltcd.locus_tags_list = ValNodeSort (ltcd.locus_tags_list, SortVnpByGlobalDiscrepancyString);
    ltcd.missing_list = ValNodeSort (ltcd.missing_list, SortVnpByGlobalDiscrepancyString);
    
    if (ltcd.missing_list != NULL) {
      dip = ReportMissingFields (ltcd.missing_list, discReportMissingLocusTags, DISC_GENE_MISSING_LOCUS_TAG);
      if (dip != NULL) {
        ValNodeAddPointer (discrepancy_list, 0, dip);
      }
    }
    dip = ReportNonUniqueGlobalDiscrepancy (ltcd.locus_tags_list, 
                                            discReportDuplicateLocusTagFmt, 
                                            discReportOneDuplicateLocusTagFmt,
                                            DISC_GENE_DUPLICATE_LOCUS_TAG,
                                            FALSE);
    if (dip != NULL) {
      dip_sub = FindAdjacentDuplicateLocusTagGenes (ltcd.locus_tags_list);
      if (dip_sub != NULL) {
        ValNodeAddPointer (&(dip->subcategories), 0, dip_sub);
      }
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }

    /* inconsistent locus tags */
    ValNodeLink (discrepancy_list,
                 ReportInconsistentGlobalDiscrepancyPrefixes (ltcd.locus_tags_list,
                                                              discReportInconsistentLocusTagPrefixFmt,
                                                              DISC_GENE_LOCUS_TAG_INCONSISTENT_PREFIX));
    /* bad formats */
    dip = ReportBadLocusTagFormat (ltcd.locus_tags_list);
    if (dip != NULL) {
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }

  ltcd.locus_tags_list = FreeGlobalDiscrepancyList (ltcd.locus_tags_list);
  ltcd.missing_list = FreeGlobalDiscrepancyList (ltcd.missing_list);
}

extern void AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  AddDiscrepanciesForMissingOrNonUniqueGeneLocusTagsEx (discrepancy_list, sep_list, FALSE);
}


static void AddDiscrepancyForNonGeneLocusTag (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR    locus_tag_list;
  GBQualPtr          qual;
  
  if (sfp == NULL || userdata == NULL || sfp->data.choice == SEQFEAT_GENE)
  {
    return;
  }
  
  locus_tag_list = (ValNodePtr PNTR) userdata;
  
  for (qual = sfp->qual; qual != NULL; qual = qual->next)
  {
    if (StringICmp(qual->qual, "locus_tag") == 0) 
    {
      ValNodeAddPointer (locus_tag_list, OBJ_SEQFEAT, sfp);
      return;
    }
  }
}

extern void AddDiscrepanciesForNonGeneLocusTags (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr locus_tag_list = NULL, vnp;
  CharPtr    bad_fmt = "%d non-gene features have locus tags.";
  ClickableItemPtr dip;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &locus_tag_list, AddDiscrepancyForNonGeneLocusTag);
  }
  
  if (locus_tag_list != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_NON_GENE_LOCUS_TAG;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (locus_tag_list));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = locus_tag_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }  
}


static Boolean ShouldCollectPseudoGeneText (CharPtr str)
{
  if (StringHasNoText (str)) {
    return FALSE;
  } else if (StringSearch (str, "hypothetical") != NULL) {
    return FALSE;
  } else if (StringSearch (str, "transposase") != NULL) {
    return FALSE;
  } else {
    return TRUE;
  }
}

static ValNodePtr MakeTokens (CharPtr str, Char sep)
{
  CharPtr cp_end, cp, txt;
  Int4       len;
  ValNodePtr list = NULL;

  if (StringHasNoText (str)) {
    return NULL;
  }
  cp = str;
  cp_end = StringChr (cp, sep);
  while (cp_end != NULL) {
    len = cp_end - cp + 1;
    txt = (CharPtr) MemNew (sizeof (Char) * len);
    StringNCpy (txt, cp, len - 1);
    txt[len - 1] = 0;
    TrimSpacesAroundString (txt);
    if (StringHasNoText (txt)) {
      txt = MemFree (txt);
    } else {
      ValNodeAddPointer (&list, 0, txt);
    }
    cp = cp_end + 1;
    cp_end = StringChr (cp, sep);
  }
  txt = StringSave (cp);
  TrimSpacesAroundString (txt);
  if (StringHasNoText (txt)) {
    txt = MemFree (txt);
  } else {
    ValNodeAddPointer (&list, 0, txt);
  }
  return list;
}


/* Note - this function assumes that the lists have been sorted */
static CharPtr FindVnpStringMatches (ValNodePtr list1, ValNodePtr list2, Boolean case_sensitive)
{
  ValNodePtr vnp1, vnp2;
  CharPtr    rval = NULL, tmp;
  Int4       cmp;

  vnp1 = list1;
  vnp2 = list2;
  while (vnp1 != NULL && vnp2 != NULL) {
    if (case_sensitive) {
      cmp = StringCmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    } else {
      cmp = StringICmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    }
    if (cmp == 0) {
      if (rval == NULL) {
        rval = StringSave (vnp1->data.ptrvalue);
      } else {
        tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (rval) + StringLen (vnp1->data.ptrvalue) + 2));
        sprintf (tmp, "%s;%s", rval, (CharPtr) vnp1->data.ptrvalue);
        rval = MemFree (rval);
        rval = tmp;
      }
      vnp1 = vnp1->next;
      vnp2 = vnp2->next;
    } else if (cmp < 1) {
      vnp1 = vnp1->next;
    } else {
      vnp2 = vnp2->next;
    }
  }
  return rval;
}


static void ExtractVnpByStringSearch (ValNodePtr PNTR list, CharPtr search)
{
  ValNodePtr vnp, prev = NULL, vnp_next;
  if (list == NULL) return;

  for (vnp = *list; vnp != NULL; vnp = vnp_next) {
    vnp_next = vnp->next;
    if (StringSearch (vnp->data.ptrvalue, search) != NULL) {
      if (prev == NULL) {
        *list = vnp_next;
      } else {
        prev->next = vnp_next;
      }
      vnp->next = NULL;
      vnp = ValNodeFreeData (vnp);
    } else {
      prev = vnp;
    }
  }
}


static CharPtr GetGeneStringMatch (CharPtr str1, CharPtr str2)
{
  ValNodePtr list1, list2;
  CharPtr    rval;

  if (!ShouldCollectPseudoGeneText(str1) || ! ShouldCollectPseudoGeneText (str2)) {
    return NULL;
  }
  list1 = MakeTokens (str1, ';');
  list2 = MakeTokens (str2, ';');
  list1 = ValNodeSort (list1, SortVnpByString);
  list2 = ValNodeSort (list2, SortVnpByString);
  
  rval = FindVnpStringMatches (list1, list2, FALSE);
  list1 = ValNodeFreeData (list1);
  list2 = ValNodeFreeData (list1);
  return rval;
} 


static void FindAdjacentPseudoGenesOnBioseq (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr sfp, sfp_next;
  SeqMgrFeatContext fcontext;
  GeneRefPtr        grp, grp_next;
  CharPtr           fmt = "Adjacent pseudogenes have the same text: %s";
  CharPtr           match_txt;
  ClickableItemPtr  cip;

  if (bsp == NULL || ISA_aa (bsp->mol)) return;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
       sfp != NULL;
       sfp = sfp_next) 
  {
    sfp_next = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
    if (sfp_next != NULL && sfp->data.choice == SEQFEAT_GENE && sfp->pseudo
        && sfp_next->data.choice == SEQFEAT_GENE && sfp_next->pseudo
        && SeqLocStrand (sfp->location) == SeqLocStrand (sfp_next->location))
    {
      match_txt = GetGeneStringMatch (sfp->comment, sfp_next->comment);
      if (match_txt == NULL) 
      {
        grp = (GeneRefPtr) sfp->data.value.ptrvalue;
        grp_next = (GeneRefPtr) sfp_next->data.value.ptrvalue;
        if (grp != NULL && grp_next != NULL) {
          match_txt = GetGeneStringMatch (grp->locus, grp_next->locus);
          if (match_txt == NULL) {
            match_txt = GetGeneStringMatch (grp->desc, grp_next->desc);
          }
        }
      }
      if (match_txt != NULL)
      {
        cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (match_txt)));
        sprintf (cip->description, fmt, match_txt);
        cip->clickable_item_type = DISC_ADJACENT_PSEUDOGENE;
        ValNodeAddPointer (&(cip->item_list), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(cip->item_list), OBJ_SEQFEAT, sfp_next);
        ValNodeAddPointer ((ValNodePtr PNTR) userdata, 0, cip);
      }
    }
  }
}

static ClickableItemPtr DiscrepancyForPairs (Uint4 item_type, CharPtr bad_fmt, ValNodePtr item_list);


static ValNodePtr SubcategoriesForIdenticalClickableItemDescriptions (ValNodePtr discrepancy_list)
{
  ClickableItemPtr cip, cip_new;
  ValNodePtr       subcategories = NULL;
  ValNodePtr       vnp_start = NULL, vnp_prev, vnp;
  CharPtr          last_str = NULL;
  CharPtr          fmt = "%d genes: %s";

  if (discrepancy_list == NULL || discrepancy_list->next == NULL) return NULL;
  discrepancy_list = ValNodeSort (discrepancy_list, SortVnpByClickableItemDescription);

  vnp_start = discrepancy_list;
  vnp_prev = vnp_start;
  cip = (ClickableItemPtr) vnp_start->data.ptrvalue;
  last_str = cip->description;

  for (vnp = discrepancy_list->next; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (StringCmp (last_str, cip->description) != 0) {
      vnp_prev->next = NULL;     
      cip_new = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      cip_new->clickable_item_type = cip->clickable_item_type;
      cip_new->item_list = ItemListFromSubcategories (vnp_start);
      cip_new->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (last_str) + 15));
      sprintf (cip_new->description, fmt, ValNodeLen (cip_new->item_list), last_str);
      cip_new->subcategories = vnp_start;
      ValNodeAddPointer (&subcategories, 0, cip_new);
      vnp_start = vnp;
      last_str = cip->description;
    }
    vnp_prev = vnp;
  }
  cip_new = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  cip_new->clickable_item_type = cip->clickable_item_type;
  cip_new->item_list = ItemListFromSubcategories (vnp_start);
  cip_new->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (last_str) + 15));
  sprintf (cip_new->description, fmt, ValNodeLen (cip_new->item_list), last_str);
  cip_new->subcategories = vnp_start;
  ValNodeAddPointer (&subcategories, 0, cip_new);
  
  return subcategories;  
}


static void FindAdjacentPseudoGenes (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr pair_list = NULL, vnp, subcategories, item_list;
  SeqEntryPtr sep;
  ClickableItemPtr  cip;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitBioseqsInSep (sep, &pair_list, FindAdjacentPseudoGenesOnBioseq);
  }
  if (pair_list != NULL) {
    subcategories = SubcategoriesForIdenticalClickableItemDescriptions (pair_list);

    item_list = ItemListFromSubcategories (subcategories);
    cip = DiscrepancyForPairs (DISC_ADJACENT_PSEUDOGENE, "%d pseudogenes match an adjacent pseudogene's text", item_list);
    cip->subcategories = subcategories;
    ValNodeAddPointer (discrepancy_list, 0, cip); 
  }
}


static void FindBioseqsWithoutAnnotationCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;

  if (bsp == NULL || userdata == NULL) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  if (sfp == NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_BIOSEQ, bsp);
  }
}

  
static void FindBioseqsWithoutAnnotation (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL, vnp;
  SeqEntryPtr sep;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitBioseqsInSep (sep, &item_list, FindBioseqsWithoutAnnotationCallback);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_NO_ANNOTATION, "%d bioseqs have no features", item_list));
  }
}


/* For Influenza A viruses, the year of the collection date should be
 * the last number before the second set of parentheses in the tax name.
 * For Influenza B viruses, the year should be the last token in the
 * tax name (tokens are separated by / characters).
 */
static Boolean DoInfluenzaStrainAndCollectionDateMisMatch (BioSourcePtr biop)
{
  CharPtr cp;
  Int4    year = 0, coll_year;
  SubSourcePtr ssp;

  if (biop == NULL || biop->org == NULL) {
    return FALSE;
  }
  if (StringNCmp (biop->org->taxname, "Influenza A virus ", 18) == 0) {
    cp = StringChr (biop->org->taxname, '(');
    if (cp != NULL) {
      cp = StringChr (cp + 1, '(');
      if (cp != NULL) {
        cp--;
        while (isspace (*cp) && cp > biop->org->taxname) {
          cp--;
        }
        if (isdigit (*cp)) {
          while (cp > biop->org->taxname + 1 && isdigit (*(cp - 1))) {
            cp--;
          }
          if (isdigit (*cp)) {
            year = atoi (cp);
          }
        }
      }
    }
  } else if (StringNCmp (biop->org->taxname, "Influenza B virus ", 18) == 0) {
    cp = StringRChr (biop->org->taxname, '/');
    if (cp != NULL) {
      cp++;
      while (isspace (*cp)) {
        cp++;
      }
      if (isdigit (*cp)) {
        year = atoi (cp);
      }
    }
  } else {
    return FALSE;
  }

  if (year > 0) {
    for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
      if (ssp->subtype == SUBSRC_collection_date) {
        cp = StringRChr (ssp->name, '-');
        if (cp == NULL) {
          coll_year = atoi (ssp->name);
        } else if (!isdigit (*(cp + 1))) {
          return TRUE;
        } else {
          coll_year = atoi (cp + 1);
        }
        if (coll_year == year) {
          return FALSE;
        } else {
          return TRUE;
        }
      }
    }
  }
  return TRUE;
} 


static void FindBioSourceDescWithInfluenzaStrainCollectionDateMismatch (SeqDescrPtr sdp, Pointer data)
{
  if (data != NULL && sdp != NULL && sdp->choice == Seq_descr_source && DoInfluenzaStrainAndCollectionDateMisMatch (sdp->data.ptrvalue)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}


static void FindBioSourceFeatWithInfluenzaStrainCollectionDateMismatch (SeqFeatPtr sfp, Pointer data)
{
  if (data != NULL && sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC && DoInfluenzaStrainAndCollectionDateMisMatch (sfp->data.value.ptrvalue)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
}


static void FindInfluenzaStrainCollectionDateMismatches (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL, vnp;
  SeqEntryPtr sep;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitDescriptorsInSep (sep, &item_list, FindBioSourceDescWithInfluenzaStrainCollectionDateMismatch);
    VisitFeaturesInSep (sep, &item_list, FindBioSourceFeatWithInfluenzaStrainCollectionDateMismatch);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_INFLUENZA_DATE_MISMATCH, "%d influenza strains conflict with collection date", item_list));
  }
}


static Boolean PosIsAt3End (Int4 pos, SeqLocPtr slp)
{
  BioseqPtr bsp;
  SeqLocPtr tmp;
  Int4      seq_end = 0;

  if ((bsp = BioseqFindFromSeqLoc (slp)) == NULL) {
    return FALSE;
  } else if (pos == bsp->length - 1) {
    return TRUE;
  } else if (bsp->repr != Seq_repr_seg || bsp->seq_ext_type != 1) {
    return FALSE;
  } else {
    for (tmp = (SeqLocPtr)bsp->seq_ext; tmp != NULL; tmp = tmp->next) {
      seq_end += SeqLocLen (tmp);
      if (pos == seq_end - 1) {
        return TRUE;
      }
    }
    bsp = BioseqFind(SeqLocId (slp));
    if (pos == bsp->length -1) {
      return TRUE;
    }
  }
  return FALSE;
}
    

static Boolean PosIsAt5End (Int4 pos, SeqLocPtr slp)
{
  BioseqPtr bsp;
  Int4      seq_end = 0;

  if (slp == NULL) {
    return FALSE;
  } else if (pos == 0) {
    return TRUE;
  } else if ((bsp = BioseqFindFromSeqLoc (slp)) == NULL
             || bsp->repr != Seq_repr_seg || bsp->seq_ext_type != 1) {
    return FALSE;
  } else {
    for (slp = (SeqLocPtr)bsp->seq_ext; slp != NULL; slp = slp->next) {
      seq_end += SeqLocLen (slp);
      if (pos == seq_end) {
        return TRUE;
      }
    }
  }
  return FALSE;
}


static void FindShortIntronsCallback (SeqFeatPtr sfp, Pointer data)
{
  SeqLocPtr slp;
  Int4      last_start, last_stop, start, stop;
  Boolean   found_short = FALSE, partial5, partial3;
  Uint1     strand;

  if (sfp == NULL || data == NULL || IsPseudo (sfp)) {
    return;
  }
  if (sfp->idx.subtype == FEATDEF_intron) {
    if (SeqLocLen (sfp->location) < 11) {
      strand = SeqLocStrand (sfp->location);
      CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
      if (partial5 && strand != Seq_strand_minus 
          && PosIsAt5End (SeqLocStart (sfp->location), sfp->location)) {
        /* partial at end of sequence, ok */
      } else if (partial3 && strand == Seq_strand_minus 
          && PosIsAt5End (SeqLocStop (sfp->location), sfp->location)) {
        /* partial at end of sequence, ok */
      } else if (partial5 && strand == Seq_strand_minus 
          && PosIsAt3End (SeqLocStart (sfp->location), sfp->location)) {
        /* partial at end of sequence, ok */
      } else if (partial3 && strand != Seq_strand_minus
          && PosIsAt3End (SeqLocStop (sfp->location), sfp->location)) {
        /* partial at end of sequence, ok */
      } else {
        ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
      }
    }
  } else if (sfp->idx.subtype == FEATDEF_CDS && !sfp->excpt) {
    slp = SeqLocFindNext (sfp->location, NULL);
    last_start = SeqLocStart (slp);
    last_stop = SeqLocStop (slp);
    slp = SeqLocFindNext (sfp->location, slp);
    while (slp != NULL && !found_short) {
      start = SeqLocStart (slp);
      stop = SeqLocStop (slp);
      if (ABS (start - last_stop) < 11 || ABS (stop - last_start) < 11) {
        found_short = TRUE;
      }
      last_start = start;
      last_stop = stop;
      slp = SeqLocFindNext (sfp->location, slp);
    }
    if (found_short) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
    }
  }
}

    
extern void FindShortIntrons (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL, vnp, with_exception = NULL;
  SeqFeatPtr sfp;
  SeqEntryPtr sep;
  Boolean     any_no_exception = FALSE;
  ClickableItemPtr cip;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitFeaturesInSep (sep, &item_list, FindShortIntronsCallback);
  }
  if (item_list != NULL) {
    for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
      sfp = vnp->data.ptrvalue;
      if (sfp != NULL && sfp->excpt) {
        ValNodeAddPointer (&with_exception, OBJ_SEQFEAT, sfp);
      } else {
        any_no_exception = TRUE;
      }
    }
    if (any_no_exception) {
      cip = NewClickableItem (DISC_SHORT_INTRON, "%d introns are shorter than 10 nt", item_list);
      if (with_exception != NULL) {
        ValNodeAddPointer (&cip->subcategories, 0, NewClickableItem (DISC_SHORT_INTRON, "%d introns are shorter than 11 nt and have an exception", with_exception));
      }
      ValNodeAddPointer (discrepancy_list, 0, cip);
    } else {
      ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_SHORT_INTRON, "%d introns are shorter than 11 nt and have an exception", with_exception));
    }
  }
}


NLM_EXTERN void AddExceptionsToShortIntrons (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  CharPtr    txt;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT && (sfp = (SeqFeatPtr) vnp->data.ptrvalue) != NULL
        && StringStr (sfp->except_text, "low-quality sequence region") == NULL) {
      SetStringValue (&(sfp->except_text), "low-quality sequence region", ExistingTextOption_append_semi);
      sfp->excpt = TRUE;
      if (lip != NULL && lip->fp != NULL) {
        txt = GetDiscrepancyItemText (vnp);
        fprintf (lip->fp, "Added low-quality sequence region exception to %s\n", txt);
        txt = MemFree (txt);
      }
    }
  }
}



static Boolean StrandOk (Uint1 strand1, Uint1 strand2)
{
  if (strand1 == Seq_strand_minus && strand2 != Seq_strand_minus) {
    return FALSE;
  } else if (strand1 != Seq_strand_minus && strand2 == Seq_strand_minus) {
    return FALSE;
  } else {
    return TRUE;
  }
}


static Boolean IsMixedStrandGeneLocationOk (SeqLocPtr feat_loc, SeqLocPtr gene_loc)
{
  SeqLocPtr         gene_subloc, feat_subloc;
  Uint1             gene_strand, feat_strand;
  Int4              gene_start, gene_stop;
  Int4              feat_start, feat_stop;

  gene_subloc = SeqLocFindNext (gene_loc, NULL);
  feat_subloc = SeqLocFindNext (feat_loc, NULL);
  while (gene_subloc != NULL && feat_subloc != NULL) {

    gene_strand = SeqLocStrand (gene_subloc);
    feat_strand = SeqLocStrand (feat_subloc);
    if (!StrandOk (gene_strand, feat_strand)) {
      return FALSE;
    }

    gene_start = SeqLocStart (gene_subloc);
    gene_stop = SeqLocStop (gene_subloc);
    feat_start = SeqLocStart (feat_subloc);
    feat_stop = SeqLocStop (feat_subloc);

    if (gene_strand == Seq_strand_minus) {
      if (gene_stop != feat_stop) {
        return FALSE;
      }
      while (gene_start != feat_start && feat_subloc != NULL) {
        while ((feat_subloc = SeqLocFindNext (feat_loc, feat_subloc)) != NULL) {
          feat_strand = SeqLocStrand (feat_subloc);
          if (!StrandOk (gene_strand, feat_strand)) {
            return FALSE;
          }
          feat_start = SeqLocStart (feat_subloc);
          if (feat_start < gene_start) {
            return FALSE;
          } else if (feat_start == gene_start) {
            break;
          }
        }
      }

    } else {
      if (gene_start != feat_start) {
        return FALSE;
      }
      while (gene_stop != feat_stop && feat_subloc != NULL) {
        while ((feat_subloc = SeqLocFindNext (feat_loc, feat_subloc)) != NULL) {
          feat_strand = SeqLocStrand (feat_subloc);
          if (!StrandOk (gene_strand, feat_strand)) {
            return FALSE;
          }
          feat_stop = SeqLocStop (feat_subloc);
          if (feat_stop > gene_stop) {
            return FALSE;
          } else if (feat_stop == gene_stop) {
            break;
          }
        }
      }
    }
    if (feat_subloc == NULL) {
      return FALSE;
    }
    gene_subloc = SeqLocFindNext (gene_loc, gene_subloc);
    feat_subloc = SeqLocFindNext (feat_loc, feat_subloc);
  }
  

  if (gene_subloc != NULL || feat_subloc != NULL) 
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


static Boolean 
IsGeneLocationOk 
(SeqFeatPtr           feat,
 SeqMgrFeatContextPtr feat_context, 
 SeqFeatPtr           gene,
 SeqMgrFeatContextPtr gene_context,
 BioseqPtr            bsp)
{
  SeqFeatPtr        rbs_sfp;
  SeqMgrFeatContext rbs_context;
  
  if (feat_context == NULL || gene_context == NULL)
  {
    return FALSE;
  } 
  else if (feat_context->mixed_strand || gene_context->mixed_strand) 
  {
    /* special handling for trans-spliced */
    return IsMixedStrandGeneLocationOk (feat->location, gene->location);
  }
  else if ((feat_context->strand == Seq_strand_minus && gene_context->strand != Seq_strand_minus)
           || (feat_context->strand != Seq_strand_minus && gene_context->strand == Seq_strand_minus))
  {
    return FALSE;
  }
  else if (gene_context->left == feat_context->left && gene_context->right == feat_context->right)
  {
    return TRUE;
  }
  else if ((gene_context->strand == Seq_strand_minus && gene_context->left == feat_context->left)
           || (gene_context->strand != Seq_strand_minus && gene_context->right == feat_context->right))
  {
    /* find RBS to extend gene on 5' end */
    for (rbs_sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_RBS, &rbs_context);
         rbs_sfp != NULL;
         rbs_sfp = SeqMgrGetNextFeature (bsp, rbs_sfp, 0, FEATDEF_RBS, &rbs_context))
    {
      if (rbs_context.strand != gene_context->strand)
      {
        continue;
      }
      if (rbs_context.strand == Seq_strand_minus)
      {
        if (rbs_context.right == gene_context->right 
            && rbs_context.left >= feat_context->right)
        {
          return TRUE;
        }
      }
      else
      {
        if (rbs_context.left == gene_context->left
            && rbs_context.right <= feat_context->left)
        {
          return  TRUE;
        }
      }
    }
  }
  return FALSE;
}

static ClickableItemPtr GeneLocationDiscrepancy (Uint1 feature_type, SeqFeatPtr gene, SeqFeatPtr sfp)
{
  ClickableItemPtr cip;

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  MemSet (cip, 0, sizeof (ClickableItemData));
  cip->clickable_item_type = DISC_GENE_CDS_mRNA_LOCATION_CONFLICT;
  if (feature_type == SEQFEAT_CDREGION) {
    cip->description = StringSave ("Coding region location does not match gene location");
  } else if (feature_type == SEQFEAT_RNA) {
    cip->description = StringSave ("RNA feature location does not match gene location");
  } else {
    cip->description = StringSave ("Feature location does not match gene location");
  }
  ValNodeAddPointer (&cip->item_list, OBJ_SEQFEAT, sfp);
  ValNodeAddPointer (&cip->item_list, OBJ_SEQFEAT, gene);
  return cip;
}

static ClickableItemPtr MissingGeneXrefDiscrepancy (Uint1 feature_type, SeqFeatPtr sfp)
{
  ClickableItemPtr cip;

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  MemSet (cip, 0, sizeof (ClickableItemData));
  cip->clickable_item_type = DISC_GENE_CDS_mRNA_LOCATION_CONFLICT;
  if (feature_type == SEQFEAT_CDREGION) {
    cip->description = StringSave ("Coding region xref gene does not exist");
  } else if (feature_type == SEQFEAT_RNA) {
    cip->description = StringSave ("RNA feature xref gene does not exist");
  } else {
    cip->description = StringSave ("Feature xref gene does not exist");
  }
  ValNodeAddPointer (&cip->item_list, OBJ_SEQFEAT, sfp);
  return cip;
}


static void
CheckFeatureTypeForLocationDiscrepancies 
(BioseqPtr       bsp, 
 Uint1           feature_type,
 Uint1           exclude_featdef,
 ValNodePtr PNTR discrepancy_list)
{
  SeqMgrFeatContext context, gene_context;
  GeneRefPtr        grp;
  SeqFeatPtr        sfp, gene_sfp;
  Boolean           found_match;
  
  if (bsp == NULL || ISA_aa (bsp->mol) || discrepancy_list == NULL || IsmRNASequenceInGenProdSet(bsp))
  {
    return;
  }

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, feature_type, 0, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, feature_type, 0, &context))
  {
    grp = SeqMgrGetGeneXref (sfp);
    if (grp == NULL)
    {
      gene_sfp = SeqMgrGetOverlappingGene (sfp->location, &gene_context);
      if (gene_sfp != NULL && !IsGeneLocationOk (sfp, &context, gene_sfp, &gene_context, bsp) && sfp->idx.subtype != exclude_featdef)
      {
        ValNodeAddPointer (discrepancy_list, 0, GeneLocationDiscrepancy(feature_type, gene_sfp, sfp));
      }
    }
    else if (!SeqMgrGeneIsSuppressed (grp))
    {
      found_match = FALSE;
      for (gene_sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, FEATDEF_GENE, &gene_context);
           gene_sfp != NULL && ! found_match;
           gene_sfp = SeqMgrGetNextFeature (bsp, gene_sfp, SEQFEAT_GENE, FEATDEF_GENE, &gene_context))
      {
        if (GeneRefMatch (gene_sfp->data.value.ptrvalue, grp) && gene_context.strand == context.strand)
        {
          if (IsGeneLocationOk (sfp, &context, gene_sfp, &gene_context, bsp))
          {
            found_match = TRUE;
          }
          else if (sfp->idx.subtype != exclude_featdef)
          {
            ValNodeAddPointer (discrepancy_list, 0, GeneLocationDiscrepancy(feature_type, gene_sfp, sfp));
          }
        }
      }
      if (!found_match) {
        ValNodeAddPointer (discrepancy_list, 0, MissingGeneXrefDiscrepancy(feature_type, sfp));
      }
    }
  }
  
}


static Boolean HasLineage (BioSourcePtr biop, CharPtr lineage)
{
  CharPtr forced_lineage;

  forced_lineage = GetAppProperty ("ReportLineage");
  if (StringISearch (forced_lineage, lineage) != NULL) 
  {
    return TRUE;
  } 
  else if (StringHasNoText (forced_lineage) 
           && biop != NULL && biop->org != NULL && biop->org->orgname != NULL
           && StringISearch (biop->org->orgname->lineage, lineage) != NULL) 
  {
    return TRUE;
  } 
  else 
  {
    return FALSE;
  }
}


static Boolean BioseqHasLineage (BioseqPtr bsp, CharPtr lineage)
{
  SeqMgrDescContext context;
  SeqDescrPtr       sdp;
  BioSourcePtr      biop;
  CharPtr           forced_lineage;

  forced_lineage = GetAppProperty ("ReportLineage");
  if (!StringHasNoText (forced_lineage)) {
    if (StringISearch (forced_lineage, lineage) != NULL) 
    {
      return TRUE;
    } else {
      return FALSE;
    } 
  } else if (bsp == NULL) {
    return FALSE;
  }
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp == NULL || (biop = sdp->data.ptrvalue) == NULL
      || biop->org == NULL
      || biop->org->orgname == NULL
      || StringISearch (biop->org->orgname->lineage, lineage) == NULL) {
    return FALSE;
  } else {
    return TRUE;
  }

}


static Boolean IsEukaryoticBioSource (BioSourcePtr biop)
{
  return HasLineage(biop, "Eukaryota");
}


static Boolean IsViralBioSource (BioSourcePtr biop)
{
  return HasLineage(biop, "Viruses");
}


static Boolean IsBacterialBioSource (BioSourcePtr biop)
{
  return HasLineage(biop, "Bacteria");
}


static Boolean IsEukaryotic (BioseqPtr bsp)
{
  SeqMgrDescContext context;
  SeqDescrPtr       sdp;
  BioSourcePtr      biop;

  if (bsp == NULL) {
    return FALSE;
  }
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp == NULL || (biop = sdp->data.ptrvalue) == NULL
      || biop->genome == GENOME_mitochondrion
      || biop->genome == GENOME_chloroplast
      || biop->genome == GENOME_plastid
      || biop->genome == GENOME_apicoplast
      || !IsEukaryoticBioSource(biop)) {
    return FALSE;
  } else {
    return TRUE;
  }
}


static void CDSmRNAGeneLocationDiscrepanciesCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR   discrepancy_list;

  if (bsp == NULL || ! ISA_na (bsp->mol) || userdata == NULL)
  {
    return;
  }
  
  discrepancy_list = (ValNodePtr PNTR) userdata;
  
  if (IsEukaryotic (bsp)) {
    CheckFeatureTypeForLocationDiscrepancies (bsp, SEQFEAT_RNA, 0, discrepancy_list);
  } else {
    CheckFeatureTypeForLocationDiscrepancies (bsp, SEQFEAT_CDREGION, 0, discrepancy_list);
    CheckFeatureTypeForLocationDiscrepancies (bsp, SEQFEAT_RNA, 0, discrepancy_list);
  }
}

static ValNodePtr ValNodePointerDup (ValNodePtr vnp);

extern void FindCDSmRNAGeneLocationDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr         subcategories = NULL, feature_list = NULL, vnp;
  CharPtr            bad_fmt = "%d features have inconsistent gene locations.";
  ClickableItemPtr dip;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &subcategories, CDSmRNAGeneLocationDiscrepanciesCallback);
  }
  
  if (subcategories != NULL)
  {
    for (vnp = subcategories; vnp != NULL; vnp = vnp->next) {
      dip = vnp->data.ptrvalue;
      if (dip != NULL && dip->item_list != NULL) {
        ValNodeLink (&feature_list, ValNodePointerDup (dip->item_list));
      }
    }
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_GENE_CDS_mRNA_LOCATION_CONFLICT;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (subcategories));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->subcategories = subcategories;

      dip->item_list = feature_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }  
}


typedef struct cdsgeneproduct 
{
  ValNodePtr cds_list;
  CharPtr    gene_locus;
  CharPtr    product_name;
} CDSGeneProductData, PNTR CDSGeneProductPtr;


static ValNodePtr CDSGeneProductListFree (ValNodePtr cds_list)
{
  CDSGeneProductPtr cgpp;
  
  if (cds_list == NULL) {
    return cds_list;
  }
  
  cds_list->next = CDSGeneProductListFree(cds_list->next);
  
  cgpp = (CDSGeneProductPtr) cds_list->data.ptrvalue;
  if (cgpp != NULL) {
    cgpp->cds_list = ValNodeFree (cgpp->cds_list);
  }
  ValNodeFreeData (cds_list);
  return NULL;
}


static CharPtr GetGeneLabel (SeqFeatPtr sfp)
{
  GeneRefPtr grp;
  SeqFeatPtr gene_sfp;
  
  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL)
  {
    gene_sfp = SeqMgrGetOverlappingGene (sfp->location, NULL);
    if (gene_sfp != NULL)
    {
      grp = gene_sfp->data.value.ptrvalue;
    }
  }
  if (grp != NULL)
  {
    if (!StringHasNoText (grp->locus))
    {
      return grp->locus;
    }
  }
  return NULL;
}

static void FindCDSGeneProductConflictsCallback (SeqFeatPtr sfp, Pointer userdata)
{
  CDSGeneProductPtr cgpp, cgpp_compare;
  SeqMgrFeatContext context;
  ValNodePtr PNTR   cds_list;
  ValNodePtr        vnp;
  Boolean           found_match = FALSE;
  CharPtr           gene_label;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || userdata == NULL)
  {
    return;
  }
  
  sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, 0, 0, sfp, &context);
  if (sfp == NULL)
  {
    return;
  }
  
  gene_label = GetGeneLabel (sfp);
  if (StringHasNoText (gene_label)) return;

  cds_list = (ValNodePtr PNTR) userdata;

  if (*cds_list == NULL) {
    cgpp = (CDSGeneProductPtr) MemNew (sizeof (CDSGeneProductData));
    if (cgpp != NULL)
    {
      ValNodeAddPointer (&(cgpp->cds_list), OBJ_SEQFEAT, sfp);
      cgpp->gene_locus = gene_label;
      cgpp->product_name = StringSave (context.label);
      ValNodeAddPointer (cds_list, 0, cgpp);
    }
  } else {
    vnp = *cds_list;
    while (vnp != NULL && !found_match)
    {
      cgpp_compare = (CDSGeneProductPtr) vnp->data.ptrvalue;
      if (cgpp_compare != NULL 
          && StringCmp (cgpp_compare->gene_locus, gene_label) == 0
          && StringCmp (cgpp_compare->product_name, context.label) != 0)
      {
        found_match = TRUE;
        vnp->choice = 1;
        ValNodeAddPointer (&(cgpp_compare->cds_list), OBJ_SEQFEAT, sfp);
      }
      vnp = vnp->next;
    }
    if (!found_match) {
      cgpp = (CDSGeneProductPtr) MemNew (sizeof (CDSGeneProductData));
      if (cgpp != NULL)
      {
        ValNodeAddPointer (&(cgpp->cds_list), OBJ_SEQFEAT, sfp);
        cgpp->gene_locus = gene_label;
        cgpp->product_name = StringSave (context.label);
        ValNodeAddPointer (cds_list, 0, cgpp);
      }
    }
  }  
}

extern void FindCDSGeneProductConflicts (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr         cds_list = NULL, non_conflict = NULL, vnp;
  CDSGeneProductPtr  cgpp;
  CharPtr            bad_fmt = "%d coding regions have the same gene name as another coding region but a different product.";
  CharPtr            bad_cat_fmt = "%d coding regions have the same gene name(%s) as another coding region but a different product.";
  ClickableItemPtr dip;
  ValNodePtr         item_list = NULL, cds_vnp;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitGenProdSetFeatures (vnp->data.ptrvalue, &cds_list, FindCDSGeneProductConflictsCallback);
  }

  /* remove CDSs without conflicts */
  non_conflict = ValNodeExtractList (&cds_list, 0);
  non_conflict = CDSGeneProductListFree (non_conflict);
  
  /* for each item, replace structure used for search with just the feature */
  for (vnp = cds_list; vnp != NULL; vnp = vnp->next)
  {
    cgpp = (CDSGeneProductPtr) vnp->data.ptrvalue;
    if (cgpp != NULL)
    {
      for (cds_vnp = cgpp->cds_list; cds_vnp != NULL; cds_vnp = cds_vnp->next) {
          ValNodeAddPointer (&item_list, OBJ_SEQFEAT, cds_vnp->data.ptrvalue);
      }
      
      cgpp->product_name = MemFree (cgpp->product_name);
      
      dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      if (dip != NULL)
      {
        dip->clickable_item_type = DISC_GENE_PRODUCT_CONFLICT;
        dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + StringLen (cgpp->gene_locus) + 15));
        sprintf (dip->description, bad_cat_fmt, ValNodeLen (cgpp->cds_list), cgpp->gene_locus == NULL ? "" : cgpp->gene_locus);
        dip->item_list = cgpp->cds_list;
        cgpp->cds_list = NULL;

        vnp->choice = 0;
        vnp->data.ptrvalue = dip;      
      } else {
        cgpp->cds_list = ValNodeFree (cgpp->cds_list);
        vnp->choice = 0;
        vnp->data.ptrvalue = NULL;
      }
      /* note - we are not freeing gene_locus because we didn't make a copy */
      cgpp = MemFree (cgpp);
    }
  }
    
  if (cds_list != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_GENE_PRODUCT_CONFLICT;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (item_list));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = item_list;
      dip->subcategories = cds_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }  
}


static void FindDuplicateGeneLocusBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;
  GeneRefPtr        grp;
  ValNodePtr        locus_list = NULL;
  ClickableItemPtr  cip;
  Char              buf[255];
  CharPtr           fixed_fmt = "%%d genes have the same locus as another gene on %s", tmp_fmt;
  
  if (bsp == NULL || IsmRNASequenceInGenProdSet (bsp) || userdata == NULL) {
    return;
  }

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_GENE, 0, &fcontext)) {
    grp = sfp->data.value.ptrvalue;
    if (grp != NULL && !StringHasNoText (grp->locus)) {
      ValNodeAddPointer (&locus_list, 0, GlobalDiscrepancyNew (grp->locus, OBJ_SEQFEAT, sfp));
    }
  }
  locus_list = ValNodeSort (locus_list, SortVnpByGlobalDiscrepancyString);
  if (locus_list != NULL) {
    SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), buf, PRINTID_FASTA_SHORT, sizeof (buf) - 1);
    tmp_fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (fixed_fmt) + StringLen (buf)));
    sprintf (tmp_fmt, fixed_fmt, buf);
    cip = ReportNonUniqueGlobalDiscrepancy (locus_list,
                                            tmp_fmt,
                                            "%d genes have locus %s",
                                            DISC_GENE_DUPLICATE_LOCUS,
                                            TRUE);
    tmp_fmt = MemFree (tmp_fmt);
    if (cip != NULL) {
      if (cip->item_list == NULL) {
        cip->item_list = ItemListFromSubcategories (cip->subcategories);
      }
      ValNodeAddPointer ((ValNodePtr PNTR) userdata, 0, cip);
    }
    locus_list = FreeGlobalDiscrepancyList (locus_list);
  }
}


extern void FindDuplicateGeneLocus (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr       disc_list = NULL, item_list, vnp;
  CharPtr          bad_fmt = "%d genes have the same locus as another gene on the same Bioseq.";
  ClickableItemPtr cip;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &disc_list, FindDuplicateGeneLocusBioseqCallback);
  }
  
  if (disc_list != NULL)
  {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (cip != NULL)
    {
      item_list = ItemListFromSubcategories (disc_list);
      cip->clickable_item_type = DISC_GENE_DUPLICATE_LOCUS;
      cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (cip->description, bad_fmt, ValNodeLen (item_list));
      cip->callback_func = NULL;
      cip->datafree_func = NULL;
      cip->callback_data = NULL;
      cip->item_list = item_list;
      cip->subcategories = disc_list;
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
  }  
  
}


static void FindECNumberNotes (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR    ec_number_features;
  BioseqPtr          prot_bsp;
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         prot_sfp;
  ProtRefPtr         prp;
  ValNodePtr         vnp;
  
  if (sfp == NULL || userdata == NULL || StringHasNoText (sfp->comment))
  {
    return;
  }
  
  ec_number_features = (ValNodePtr PNTR) userdata;
  
  if (LookForECnumberPattern (sfp->comment))
  {
    ValNodeAddPointer (ec_number_features, OBJ_SEQFEAT, sfp);
  }
  else if (sfp->data.choice == SEQFEAT_CDREGION && sfp->product != NULL) 
  {
    prot_bsp = BioseqFindFromSeqLoc(sfp->product);
    prot_sfp = SeqMgrGetNextFeature(prot_bsp, NULL, SEQFEAT_PROT, FEATDEF_PROT, &fcontext);
    if (prot_sfp != NULL && prot_sfp->data.value.ptrvalue != NULL) {
      prp = (ProtRefPtr) prot_sfp->data.value.ptrvalue;
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        if (LookForECnumberPattern (vnp->data.ptrvalue)) {
          ValNodeAddPointer (ec_number_features, OBJ_SEQFEAT, sfp);
          return;
        }
      }
      if (LookForECnumberPattern (prp->desc)) {
        ValNodeAddPointer (ec_number_features, OBJ_SEQFEAT, sfp);
        return;
      }
    }
  }  
}

extern void AddECNumberNoteDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr ec_number_features = NULL, vnp;
  ClickableItemPtr dip;
  CharPtr            bad_fmt = "%d features have EC numbers in notes or products.";
  
  if (discrepancy_list == NULL) return;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &ec_number_features, FindECNumberNotes);
  }
  
  if (ec_number_features != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_EC_NUMBER_NOTE;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (ec_number_features));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = ec_number_features;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
  
}


static void FindPseudoDiscrepanciesCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR pseudo_features;
  GeneRefPtr      grp;
  SeqFeatPtr      gene_sfp = NULL;
  
  if (sfp == NULL || (sfp->data.choice != SEQFEAT_CDREGION && sfp->data.choice != SEQFEAT_RNA)
      || userdata == NULL)
  {
    return;
  }
  
  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL)
  {
    return;
  }
  
  gene_sfp = SeqMgrGetOverlappingGene (sfp->location, NULL);
  if (gene_sfp == NULL)
  {
    return;
  }
  
  if (sfp->pseudo && ! gene_sfp->pseudo)
  {
    pseudo_features = (ValNodePtr PNTR) userdata;
    ValNodeAddPointer (pseudo_features, OBJ_SEQFEAT, sfp);
    if (gene_sfp != NULL)
    {
      ValNodeAddPointer (pseudo_features, OBJ_SEQFEAT, gene_sfp);
    }
  }
}


extern void FindPseudoDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr pseudo_features = NULL, vnp;
  ClickableItemPtr dip;
  CharPtr            bad_fmt = "%d CDSs, RNAs, and genes have mismatching pseudos.";
  
  if (discrepancy_list == NULL) return;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &pseudo_features, FindPseudoDiscrepanciesCallback);
  }
  
  if (pseudo_features != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_PSEUDO_MISMATCH;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (pseudo_features));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = pseudo_features;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
  
}



static Boolean IsProtRefEmpty (ProtRefPtr prp)
{
  if (prp == NULL) {
    return TRUE;
  } else if (prp->name != NULL || prp->desc != NULL || prp->ec != NULL
    || prp->activity != NULL || prp->db != NULL || prp->processed != 0) {
    return FALSE;
  } else {
    return TRUE;
  }
}


NLM_EXTERN void OncallerToolPseudoDiscrepanciesFix (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr vnp, entityIDList = NULL;
  ValNode    vn;
  SeqFeatPtr sfp, mrna;
  CharPtr    feat_txt;
  SeqMgrFeatContext fcontext;
  SeqFeatXrefPtr xref, prev_xref, next_xref;
  ValNodePtr  next_name;
  ProtRefPtr  prp;
  RnaRefPtr   rrp;

  MemSet (&vn, 0, sizeof (ValNode));
  vn.choice = OBJ_SEQFEAT;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT) {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      ValNodeAddInt (&entityIDList, 0, sfp->idx.entityID);
      if (!sfp->pseudo) {
        if (lip != NULL && lip->fp != NULL) {
          feat_txt = GetDiscrepancyItemText (vnp);
          fprintf (lip->fp, "Added pseudo to %s\n", feat_txt);
          feat_txt = MemFree (feat_txt);
          lip->data_in_log = TRUE;
        }
        sfp->pseudo = TRUE;
      }
      if (sfp->data.choice == SEQFEAT_CDREGION) {
        mrna = SeqMgrGetOverlappingmRNA (sfp->location, &fcontext);
        if (mrna != NULL && !mrna->pseudo) {
          if (lip != NULL && lip->fp != NULL) {
            vn.data.ptrvalue = mrna;
            feat_txt = GetDiscrepancyItemText (&vn);
            fprintf (lip->fp, "Added pseudo to %s\n", feat_txt);
            feat_txt = MemFree (feat_txt);
            lip->data_in_log = TRUE;
          }
          mrna->pseudo = TRUE;
          /* move mRNA product to comment */
          if ((rrp = (RnaRefPtr) mrna->data.value.ptrvalue) != NULL
            && rrp->ext.choice == 1) {
            SetStringValue (&(mrna->comment), rrp->ext.value.ptrvalue, ExistingTextOption_append_semi);
            rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
            rrp->ext.choice = 0;
          }
        }
        
        /* move CDS protein name to comment */
        prev_xref = NULL;
        for (xref = sfp->xref; xref != NULL; xref = next_xref) {
          next_xref = xref->next;
          if (xref->data.choice == SEQFEAT_PROT
              && (prp = (ProtRefPtr) xref->data.value.ptrvalue) != NULL
              && prp->name != NULL
              && !StringHasNoText (prp->name->data.ptrvalue)) {
            SetStringValue (&(sfp->comment), prp->name->data.ptrvalue, ExistingTextOption_append_semi);
            prp->name->data.ptrvalue = MemFree (prp->name->data.ptrvalue);
            next_name = prp->name->next;
            prp->name->next = NULL;
            prp->name = ValNodeFreeData (prp->name);
            prp->name = next_name;
            if (IsProtRefEmpty(prp)) {
              if (prev_xref == NULL) {
                sfp->xref = next_xref;
              } else {
                prev_xref->next = next_xref;
              }
              xref->next = NULL;
              xref = SeqFeatXrefFree (xref);
            } else {
              prev_xref = xref;
            }
          } else {
            prev_xref = xref;
          }
        }            
      }
    }
  }

  entityIDList = ValNodeSort (entityIDList, SortByIntvalue);
  ValNodeUnique (&entityIDList, SortByIntvalue, ValNodeFree);
  
  for (vnp = entityIDList; vnp != NULL; vnp = vnp->next) {
    ObjMgrSetDirtyFlag (vnp->data.intvalue, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, vnp->data.intvalue, 0, 0);    
  }

}


static void FindJoinedLocations (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR joined_features;
  SeqMgrFeatContext context;
  SeqFeatPtr sfp;
  
  if (bsp == NULL || userdata == NULL || IsEukaryotic (bsp)) 
  {
    return;
  }

  joined_features = (ValNodePtr PNTR) userdata;

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context)) 
  {
    if (sfp->location != NULL 
        && (sfp->location->choice == SEQLOC_MIX || sfp->location->choice == SEQLOC_PACKED_INT)) 
    {
      ValNodeAddPointer (joined_features, OBJ_SEQFEAT, sfp);
    }
  }
}


static int CompareFeaturesByException (SeqFeatPtr sfp1, SeqFeatPtr sfp2)
{
  int         rval = 0;

  if (sfp1 == NULL || sfp2 == NULL) {
    return 0;
  }
  if (sfp1->excpt && !sfp2->excpt) {
    rval = -1;
  } else if (!sfp1->excpt && sfp2->excpt) {
    rval = 1;
  } else if (!sfp1->excpt && !sfp2->excpt) {
    rval = 0;
  } else {
    rval = StringICmp (sfp1->except_text, sfp2->except_text);
  }
  return rval;
}


static int LIBCALLBACK SortVnpByFeatureExceptionAndLocation (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  SeqFeatPtr  sfp1, sfp2;
  SeqMgrFeatContext fcontext1, fcontext2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL 
        && (sfp1 = vnp1->data.ptrvalue) != NULL
        && (sfp2 = vnp2->data.ptrvalue) != NULL) {
      rval = CompareFeaturesByException (sfp1, sfp2);
      if (rval == 0) {
        sfp1 = SeqMgrGetDesiredFeature (sfp1->idx.entityID, NULL, sfp1->idx.itemID, 0, sfp1, &fcontext1);
        sfp2 = SeqMgrGetDesiredFeature (sfp2->idx.entityID, NULL, sfp2->idx.itemID, 0, sfp2, &fcontext2);
        if (fcontext1.left < fcontext2.left) {
          rval = -1;
        } else if (fcontext1.left > fcontext2.left) {
          rval = 1;
        }
      }
    }
  }
  return rval;
}


extern void AddJoinedFeatureDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr joined_features = NULL, vnp, dup_list, subcat_item_list, vnp_next;
  ValNodePtr subcategories = NULL;
  SeqFeatPtr sfp;
  
  ClickableItemPtr dip;
  CharPtr            bad_fmt = "%d features have joined locations.";
  CharPtr            exception_fmt = "%d features have joined location but exception '%s'";
  
  if (discrepancy_list == NULL) return;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &joined_features, FindJoinedLocations);
  }
  
  if (joined_features != NULL)
  {    
    /* get list of joined features, separated into exception categories */
    dup_list = ValNodeCopyPtr (joined_features);
    dup_list = ValNodeSort (dup_list, SortVnpByFeatureExceptionAndLocation);
    subcat_item_list = dup_list;
    for (vnp = dup_list; vnp != NULL; vnp = vnp_next) {
      vnp_next = vnp->next;
      if (vnp_next == NULL || CompareFeaturesByException(vnp->data.ptrvalue, vnp_next->data.ptrvalue) != 0) {
        vnp->next = NULL;
        sfp = vnp->data.ptrvalue;
        if (!sfp->excpt) {
          ValNodeAddPointer (&subcategories, 0, 
                             NewClickableItem (DISC_JOINED_FEATURES, 
                                               "%d features have joined location but no exception",
                                               subcat_item_list));
        } else if (StringHasNoText (sfp->except_text)) {
          ValNodeAddPointer (&subcategories, 0, 
                             NewClickableItem (DISC_JOINED_FEATURES, 
                                               "%d features have joined location but a blank exception",
                                               subcat_item_list));
        } else {
          dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
          dip->clickable_item_type = DISC_JOINED_FEATURES;
          dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (exception_fmt) + StringLen (sfp->except_text) + 15));
          sprintf (dip->description, exception_fmt, ValNodeLen (subcat_item_list), sfp->except_text);
          dip->callback_func = NULL;
          dip->datafree_func = NULL;
          dip->callback_data = NULL;
          dip->item_list = subcat_item_list;
          ValNodeAddPointer (&subcategories, 0, dip);
        }
        subcat_item_list = vnp_next;
      }
    }

    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_JOINED_FEATURES;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (joined_features));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = joined_features;
      dip->subcategories = subcategories;

      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}


static void FindOverlappingGenes (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr         sfp, sfp_compare;
  SeqMgrFeatContext  context;
  ValNodePtr PNTR    overlapping_genes = NULL, non_overlap;
  ValNodePtr         gene_list = NULL, vnp, vnp_next;
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  overlapping_genes = (ValNodePtr PNTR) userdata;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, FEATDEF_GENE, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_GENE, FEATDEF_GENE, &context))
  {
    ValNodeAddPointer (&gene_list, 0, sfp);
  }
  
  for (vnp = gene_list; vnp != NULL && vnp->next != NULL; vnp = vnp->next)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    for (vnp_next = vnp->next; vnp_next != NULL; vnp_next = vnp_next->next)
    {
      sfp_compare = (SeqFeatPtr) vnp_next->data.ptrvalue;
      
      if (SeqLocStrand (sfp->location) != SeqLocStrand (sfp_compare->location))
      {
        continue;
      }
      
      if (SeqLocCompare (sfp->location, sfp_compare->location) != SLC_NO_MATCH)
      {
        vnp->choice = OBJ_SEQFEAT;
        vnp_next->choice = OBJ_SEQFEAT;
      }
    }
  }
  
  non_overlap = ValNodeExtractList (&gene_list, 0);
  non_overlap = ValNodeFree (non_overlap);
  ValNodeLink (overlapping_genes, gene_list);
  
}


extern void AddOverlappingGeneDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip;
  CharPtr            bad_fmt = "%d genes overlap another gene on the same strand.";
  ValNodePtr         overlapping_genes = NULL, vnp;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &overlapping_genes, FindOverlappingGenes);
  }
  
  if (overlapping_genes != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_OVERLAPPING_GENES;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (overlapping_genes));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = overlapping_genes;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}


typedef struct cdsoverlap 
{
  CharPtr    product_name;
  SeqFeatPtr sfp;  
  Int4       left;
  Int4       right;
} CDSOverlapData, PNTR CDSOverlapPtr;


static CDSOverlapPtr CDSOverlapNew (SeqFeatPtr sfp, CharPtr product_name, Int4 left, Int4 right)
{
  CDSOverlapPtr cop;
  
  cop = (CDSOverlapPtr) MemNew (sizeof (CDSOverlapData));
  if (cop != NULL)
  {
    cop->product_name = StringSave (product_name);
    cop->sfp = sfp;
    cop->left = left;
    cop->right = right;
  }
  return cop;
}


static ValNodePtr FreeCDSOverlapList (ValNodePtr vnp)
{
  CDSOverlapPtr cop;
  
  if (vnp != NULL)  
  {
    vnp->next = FreeCDSOverlapList (vnp->next);
    cop = (CDSOverlapPtr) vnp->data.ptrvalue;
    if (cop != NULL)
    {
      cop->product_name = MemFree (cop->product_name);
      cop = MemFree (cop);
      vnp->data.ptrvalue = NULL;
    }
    vnp = ValNodeFree (vnp);
  }
  return vnp;
}


static ValNodePtr FeatureListFromOverlapList (ValNodePtr vnp)
{
  ValNodePtr     feat_list = NULL;
  CDSOverlapPtr cop;
  
  while (vnp != NULL)
  {
    if (vnp->choice != 0 && vnp->data.ptrvalue != NULL)
    {
      cop = (CDSOverlapPtr) vnp->data.ptrvalue;
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, cop->sfp);
    }
    vnp = vnp->next;
  }
  return feat_list;
}


static CharPtr similar_product_words[] = 
{ "transposase",
  "integrase"
};

const int num_similar_product_words = sizeof (similar_product_words) / sizeof (CharPtr);

static CharPtr ignore_similar_product_words[] = 
{ "hypothetical protein",
  "phage",
  "predicted protein"
};

const int num_ignore_similar_product_words = sizeof (ignore_similar_product_words) / sizeof (CharPtr);


static Boolean OverlappingProductNameSimilar (CharPtr str1, CharPtr str2)
{
  Int4 i;
  Boolean str1_has_similarity_word = FALSE, str2_has_similarity_word = FALSE;
  
  if (StringHasNoText (str1) && StringHasNoText (str2))
  {
    return TRUE;
  }
  else if (StringHasNoText (str1) || StringHasNoText (str2))
  {
    return FALSE;
  }
  
  /* if both product names contain one of the special case similarity words,
   * the product names are similar. */
  for (i = 0; i < num_similar_product_words; i++)
  {
    if (StringISearch (str1, similar_product_words [i]) != NULL)
    {
      str1_has_similarity_word = TRUE;
    }
    if (StringISearch (str2, similar_product_words [i]) != NULL)
    {
      str2_has_similarity_word = TRUE;
    }
  }
  if (str1_has_similarity_word && str2_has_similarity_word)
  {
    return TRUE;
  }
  
  /* otherwise, if one of the product names contains one of special ignore similarity
   * words, the product names are not similar.
   */
  for (i = 0; i < num_ignore_similar_product_words; i++)
  {
    if (StringISearch (str1, ignore_similar_product_words[i]) != NULL
        || StringISearch (str2, ignore_similar_product_words[i]) != NULL)
    {
      return FALSE;
    }
  }
  
  if (StringICmp (str1, str2) == 0)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static void RemoveCodingRegionsWithSuppressionWords (ValNodePtr PNTR cds_list)
{
  FeatureFieldPtr    field;
  CharPtr            product;
  ValNodePtr         prev = NULL, vnp, vnp_next;
  SeqFeatPtr         cds;

  if (cds_list == NULL || *cds_list == NULL) {
    return;
  }

  field = FeatureFieldNew ();
  field->type = Macro_feature_type_cds;
  field->field = ValNodeNew (NULL);
  field->field->choice = FeatQualChoice_legal_qual;
  field->field->data.intvalue = Feat_qual_legal_product;

  for (vnp = *cds_list; vnp != NULL; vnp = vnp_next) {
    vnp_next = vnp->next;
    cds = (SeqFeatPtr) vnp->data.ptrvalue;
    product = GetQualFromFeature (cds, field, NULL);
    if (DoesStringContainPhrase (product, "ABC", TRUE, TRUE)
        || DoesStringContainPhrase (product, "transposon", FALSE, FALSE) 
        || DoesStringContainPhrase (product, "transposase", FALSE, FALSE)) {
      if (prev == NULL) {
        *cds_list = vnp_next;
      } else {
        prev->next = vnp_next;
      }
      vnp->next = NULL;
      vnp = ValNodeFree (vnp);
    } else {
      prev = vnp;
    }
    product = MemFree (product);

  }
  field = FeatureFieldFree (field);
}


static void FindOverlappingCDSs (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr         sfp;
  SeqMgrFeatContext  context;
  ValNodePtr PNTR    overlapping_cds = NULL, cds_list;
  ValNodePtr         overlap_list = NULL, vnp, vnp_next;
  CDSOverlapPtr      cop, cop_compare;
  Uint1              strand1, strand2;
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  overlapping_cds = (ValNodePtr PNTR) userdata;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, FEATDEF_CDS, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, FEATDEF_CDS, &context))
  {
    ValNodeAddPointer (&overlap_list, 0, CDSOverlapNew (sfp, context.label, context.left, context.right));
  }
  
  for (vnp = overlap_list; vnp != NULL && vnp->next != NULL; vnp = vnp->next)
  {
    cop = (CDSOverlapPtr) vnp->data.ptrvalue;
    if (cop == NULL)
    {
      continue;
    }
    for (vnp_next = vnp->next; vnp_next != NULL; vnp_next = vnp_next->next)
    {
      cop_compare = (CDSOverlapPtr) vnp_next->data.ptrvalue;
      if (cop_compare == NULL)
      {
        continue;
      }
      else if (cop_compare->left > cop->right)
      {
        break;
      }
      if (!OverlappingProductNameSimilar (cop->product_name, cop_compare->product_name))
      {
        continue;
      }
      strand1 = SeqLocStrand (cop->sfp->location);
      strand2 = SeqLocStrand (cop_compare->sfp->location);
      if ((strand1 == Seq_strand_minus && strand2 != Seq_strand_minus)
          || (strand1 != Seq_strand_minus && strand2 == Seq_strand_minus))
      {
        continue;
      }
      
      if (SeqLocCompare (cop->sfp->location, cop_compare->sfp->location) != SLC_NO_MATCH)
      {
        vnp->choice = OBJ_SEQFEAT;
        vnp_next->choice = OBJ_SEQFEAT;
      }
    }
  }
  
  cds_list = FeatureListFromOverlapList(overlap_list);
  /* remove features with suppression words */
  RemoveCodingRegionsWithSuppressionWords (&cds_list);

  if (cds_list != NULL)
  {
    ValNodeLink (overlapping_cds, cds_list);
  }
  overlap_list = FreeCDSOverlapList (overlap_list);
}

const CharPtr kOverlappingCDSNoteText = "overlaps another CDS with the same product name";
const CharPtr kOverlappingCDSNeedsNoteFmt = "%d coding regions overlap another coding region with a similar or identical name that do not have the appropriate note text";

static Boolean HasOverlapComment (SeqFeatPtr sfp)
{
  if (sfp == NULL || StringHasNoText (sfp->comment)) {
    return FALSE;
  }
  if (StringISearch (sfp->comment, "overlap") != NULL
      || StringISearch (sfp->comment, "frameshift") != NULL
      || StringISearch (sfp->comment, "frame shift") != NULL
      || StringISearch (sfp->comment, "extend") != NULL) {
    return TRUE;
  } else {
    return FALSE;
  }
}


extern void AddOverlappingCodingRegionDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip, sub;
  CharPtr            bad_fmt = "%d coding regions overlap another coding region with a similar or identical name.";
  ValNodePtr         overlapping_cds = NULL, vnp;
  ValNodePtr         comment_cds = NULL, no_comment_cds = NULL;
  SeqFeatPtr         cds;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &overlapping_cds, FindOverlappingCDSs);
  }
  
  if (overlapping_cds != NULL)
  {
    dip = NewClickableItem (DISC_OVERLAPPING_CDS, bad_fmt, overlapping_cds);

    if (dip != NULL)
    {
      ValNodeAddPointer (discrepancy_list, 0, dip);
      /* suppress coding regions that have ABC, transposon, or transposase in the name */
      /* create subcategories:
       * has overlap comment
       * does not have overlap comment
       */

      for (vnp = overlapping_cds; vnp != NULL; vnp = vnp->next)
      {
        if (vnp->choice == OBJ_SEQFEAT)
        {
          cds = vnp->data.ptrvalue;
          if (HasOverlapComment (cds)) {
            ValNodeAddPointer (&comment_cds, OBJ_SEQFEAT, cds);
          } else {
            ValNodeAddPointer (&no_comment_cds, OBJ_SEQFEAT, cds);
          }
        }
      }
      if (no_comment_cds != NULL)
      {
        sub = NewClickableItem (DISC_OVERLAPPING_CDS, kOverlappingCDSNeedsNoteFmt, no_comment_cds);
        if (sub != NULL) {
          ValNodeAddPointer (&dip->subcategories, 0, sub);
        }
      }
      if (comment_cds != NULL)
      {
        sub = NewClickableItem (DISC_OVERLAPPING_CDS, "%d coding regions overlap another coding region with a similar or identical name but have the appropriate note text", comment_cds);
        if (sub != NULL) {
          ValNodeAddPointer (&dip->subcategories, 0, sub);
        }
      }
    }
  }
}


NLM_EXTERN void MarkOverlappingCDSs (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  CharPtr    feat_txt;
  Boolean    has_title = FALSE;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT 
        && (sfp = vnp->data.ptrvalue) != NULL 
        && sfp->data.choice == SEQFEAT_CDREGION
        && StringISearch (sfp->comment, kOverlappingCDSNoteText) == NULL) {
      SetStringValue (&(sfp->comment), kOverlappingCDSNoteText, ExistingTextOption_append_semi);
      if (lip != NULL && lip->fp != NULL) {
        if (!has_title) {
          fprintf (lip->fp, "Added \"overlaps another CDS with the same product name\" to CDS note for overlapping CDSs with similar product names\n");
          has_title = TRUE;
        }

        feat_txt = GetDiscrepancyItemText (vnp);
        fprintf (lip->fp, "Added overlapping CDS note to %s", feat_txt);
        feat_txt = MemFree (feat_txt);
        lip->data_in_log = TRUE;
      }
    }
  }
  if (has_title) {
    fprintf (lip->fp, "\n");
  }
}


typedef struct twolists {
  ValNodePtr first_list;
  ValNodePtr second_list;
} TwoListsData, PNTR TwoListsPtr;

static void FindContainedCDSs (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr         sfp, sfp_compare;
  SeqMgrFeatContext  context;
  TwoListsPtr        two_lists;
  ValNodePtr         cds_list = NULL;
  ValNodePtr         contained_list_this_strand = NULL, contained_list_other_strand = NULL, vnp, vnp_next, last;
  ValNodePtr         last_this_strand = NULL, last_other_strand = NULL;
  Int2               loc_compare;
  Uint1              strand, strand_compare;
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  two_lists = (TwoListsPtr) userdata;
  last = NULL;
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, FEATDEF_CDS, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, FEATDEF_CDS, &context))
  {
    ValNodeAddPointer (&last, OBJ_SEQFEAT, sfp);
    if (cds_list == NULL) 
    {
      cds_list = last;
    }
  }
  
  last = NULL;
  for (vnp = cds_list; vnp != NULL && vnp->next != NULL; vnp = vnp->next)
  {    
    sfp = vnp->data.ptrvalue;
    strand = SeqLocStrand (sfp->location);
    for (vnp_next = vnp->next; vnp_next != NULL; vnp_next = vnp_next->next)
    {
      sfp_compare = vnp_next->data.ptrvalue;
      strand_compare = SeqLocStrand (sfp_compare->location);
      loc_compare = SeqLocCompare (sfp->location, sfp_compare->location);
      if (loc_compare == SLC_A_IN_B || loc_compare == SLC_B_IN_A || loc_compare == SLC_A_EQ_B)
      {
        if (StrandOk (strand, strand_compare)) {
          if (!AlreadyInList (contained_list_this_strand, sfp)) 
          {
            ValNodeAddPointer (&last_this_strand, OBJ_SEQFEAT, sfp);
            if (contained_list_this_strand == NULL) 
            {
              contained_list_this_strand = last_this_strand;
            }
          }
          if (!AlreadyInList (contained_list_this_strand, sfp_compare)) 
          {
            ValNodeAddPointer (&last_this_strand, OBJ_SEQFEAT, sfp_compare);
            if (contained_list_this_strand == NULL) 
            {
              contained_list_this_strand = last_this_strand;
            }
          }
        } else {
          if (!AlreadyInList (contained_list_other_strand, sfp)) 
          {
            ValNodeAddPointer (&last_other_strand, OBJ_SEQFEAT, sfp);
            if (contained_list_other_strand == NULL) 
            {
              contained_list_other_strand = last_other_strand;
            }
          }
          if (!AlreadyInList (contained_list_other_strand, sfp_compare)) 
          {
            ValNodeAddPointer (&last_other_strand, OBJ_SEQFEAT, sfp_compare);
            if (contained_list_other_strand == NULL) 
            {
              contained_list_other_strand = last_other_strand;
            }
          }
        }
      }
    }
  }
  cds_list = ValNodeFree (cds_list);
  
  ValNodeLink (&(two_lists->first_list), contained_list_this_strand);
  ValNodeLink (&(two_lists->second_list), contained_list_other_strand);
}


extern void AddContainedCodingRegionDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip;
  CharPtr          bad_fmt = "%d coding regions are completely contained in another coding region.";
  CharPtr          same_strand_fmt = "%d coding regions are completely contained in another coding region on the same strand.";
  CharPtr          other_strand_fmt = "%d coding regions are completely contained in another coding region, but on the opposite strand.";
  TwoListsData     two_lists;
  ValNodePtr       vnp, subcategories = NULL, item_list;
  

  if (discrepancy_list == NULL)
  {
    return;
  }

  MemSet (&two_lists, 0, sizeof (TwoListsData));
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &two_lists, FindContainedCDSs);
  }
  
  if (two_lists.first_list != NULL && two_lists.second_list != NULL) {
    ValNodeAddPointer (&subcategories, 0, NewClickableItem (DISC_CONTAINED_CDS, same_strand_fmt, two_lists.first_list));
    ValNodeAddPointer (&subcategories, 0, NewClickableItem (DISC_CONTAINED_CDS, other_strand_fmt, two_lists.second_list));
    item_list = ItemListFromSubcategories (subcategories);
    dip = NewClickableItem (DISC_CONTAINED_CDS, bad_fmt, item_list);
    dip->subcategories = subcategories;
    ValNodeAddPointer (discrepancy_list, 0, dip);
  } else if (two_lists.first_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_CONTAINED_CDS, same_strand_fmt, two_lists.first_list));
  } else if (two_lists.second_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_CONTAINED_CDS, other_strand_fmt, two_lists.second_list));
  }
}


typedef struct cdsrnaoverlap {
  ValNodePtr cds_in_rna;
  ValNodePtr rna_in_cds;
  ValNodePtr exact_match;
  ValNodePtr overlap_same_strand;
  ValNodePtr overlap_opp_strand;
  ValNodePtr overlap;
  ValNodePtr all;
} CDSRNAOverlapData, PNTR CDSRNAOverlapPtr;

static void FindCDSRNAOverlaps (BioseqPtr bsp, Pointer data)
{
  CDSRNAOverlapPtr  p;
  SeqFeatPtr        sfp, rna;
  ValNodePtr        rna_list = NULL, vnp;
  SeqMgrFeatContext fcontext;
  Int2              cmp;
  Uint1             strand1, strand2;

  if (bsp == NULL || data == NULL) return;

  p = (CDSRNAOverlapPtr) data;

  for (rna = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_RNA, 0, &fcontext);
       rna != NULL;
       rna = SeqMgrGetNextFeature (bsp, rna, SEQFEAT_RNA, 0, &fcontext))
  {
    if (rna->idx.subtype == FEATDEF_mRNA) continue;
    ValNodeAddPointer (&rna_list, OBJ_SEQFEAT, rna);
  }

  if (rna_list == NULL) return;

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext))
  {
    for (vnp = rna_list; vnp != NULL; vnp = vnp->next)
    {
      rna = (SeqFeatPtr) vnp->data.ptrvalue;
      cmp = SeqLocCompare (sfp->location, rna->location);
      if (cmp == SLC_A_EQ_B)
      {
        ValNodeAddPointer (&(p->exact_match), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(p->exact_match), OBJ_SEQFEAT, rna);
        ValNodeAddPointer (&(p->all), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(p->all), OBJ_SEQFEAT, rna);
      }
      else if (cmp == SLC_A_IN_B)
      {
        ValNodeAddPointer (&(p->cds_in_rna), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(p->cds_in_rna), OBJ_SEQFEAT, rna);
        ValNodeAddPointer (&(p->all), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(p->all), OBJ_SEQFEAT, rna);
      }
      else if (cmp == SLC_B_IN_A)
      {
        ValNodeAddPointer (&(p->rna_in_cds), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(p->rna_in_cds), OBJ_SEQFEAT, rna);
        ValNodeAddPointer (&(p->all), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(p->all), OBJ_SEQFEAT, rna);
      }
      else if (cmp != SLC_NO_MATCH)
      {
        strand1 = SeqLocStrand (sfp->location);
        strand2 = SeqLocStrand (rna->location);
        if ((strand1 == Seq_strand_minus && strand2 != Seq_strand_minus)
          || (strand2 == Seq_strand_minus && strand1 != Seq_strand_minus))
        {
          ValNodeAddPointer (&(p->overlap_opp_strand), OBJ_SEQFEAT, sfp);
          ValNodeAddPointer (&(p->overlap_opp_strand), OBJ_SEQFEAT, rna);
        }
        else
        {
          ValNodeAddPointer (&(p->overlap_same_strand), OBJ_SEQFEAT, sfp);
          ValNodeAddPointer (&(p->overlap_same_strand), OBJ_SEQFEAT, rna);
        }
        ValNodeAddPointer (&(p->overlap), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(p->overlap), OBJ_SEQFEAT, rna);
        ValNodeAddPointer (&(p->all), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(p->all), OBJ_SEQFEAT, rna);
      }
    }
  }
  rna_list = ValNodeFree (rna_list);
}


static ClickableItemPtr DiscrepancyForPairs (Uint4 item_type, CharPtr bad_fmt, ValNodePtr item_list)
{
  ClickableItemPtr dip;
  Int4             num_feat;

  if (StringHasNoText (bad_fmt)) return NULL;
  dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  dip->clickable_item_type = item_type;
  dip->item_list = item_list;
  num_feat = ValNodeLen (dip->item_list) / 2;
  dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
  sprintf (dip->description, bad_fmt, num_feat);
  return dip;
}


extern void AddRNACDSOverlapDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip, overlap_dip;
  ValNodePtr       vnp;
  CDSRNAOverlapData d;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  MemSet (&d, 0, sizeof (CDSRNAOverlapData));
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &d, FindCDSRNAOverlaps);
  }
  
  if (d.all != NULL)
  {
    dip = DiscrepancyForPairs (DISC_RNA_CDS_OVERLAP, "%d coding regions overlap RNA features", d.all);
    if (d.exact_match != NULL)
    {
      ValNodeAddPointer (&(dip->subcategories), 0, 
                         DiscrepancyForPairs (DISC_RNA_CDS_OVERLAP, 
                                              "%d coding region locations exactly match an RNA location",
                                              d.exact_match));
    }
    if (d.cds_in_rna != NULL)
    {
      ValNodeAddPointer (&(dip->subcategories), 0, 
                         DiscrepancyForPairs (DISC_RNA_CDS_OVERLAP,
                                              "%d coding regions are completely contained in RNAs",
                                              d.cds_in_rna));
    }
    if (d.rna_in_cds != NULL)
    {
      ValNodeAddPointer (&(dip->subcategories), 0, 
                         DiscrepancyForPairs (DISC_RNA_CDS_OVERLAP, 
                                              "%d coding regions completely contain RNAs",
                                              d.rna_in_cds));
    }
    if (d.overlap != NULL)
    {
      overlap_dip = DiscrepancyForPairs (DISC_RNA_CDS_OVERLAP, 
                                          "%d coding regions overlap RNAs (no containment)",
                                          d.overlap);
      if (d.overlap_same_strand != NULL) {
        ValNodeAddPointer (&(overlap_dip->subcategories), 0, 
                            DiscrepancyForPairs (DISC_RNA_CDS_OVERLAP, 
                                              "%d coding regions overlap RNAs on the same strand (no containment)",
                                              d.overlap_same_strand));
      }
      if (d.overlap_opp_strand != NULL) {
        ValNodeAddPointer (&(overlap_dip->subcategories), 0, 
                            DiscrepancyForPairs (DISC_RNA_CDS_OVERLAP, 
                                              "%d coding regions overlap RNAs on the opposite strand (no containment)",
                                              d.overlap_opp_strand));                                                   
      }
      ValNodeAddPointer (&(dip->subcategories), 0, overlap_dip);
    }
    
    ValNodeAddPointer (discrepancy_list, 0, dip);
  }
}


static void FindShortContigsCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR bioseq_list;
  
  if (bsp == NULL || !ISA_na (bsp->mol) || userdata == NULL || bsp->length >= 200) {
    return;
  }

  if (IsmRNASequenceInGenProdSet (bsp)) {
    return;
  }
  
  bioseq_list = (ValNodePtr PNTR) userdata;
  
  ValNodeAddPointer (bioseq_list, OBJ_BIOSEQ, bsp);
}

extern void FindShortContigs (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip;
  CharPtr            bad_fmt = "%d contigs are shorter than 200 nt.";
  ValNodePtr         bioseq_list = NULL, vnp;
  
  if (discrepancy_list == NULL) return;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &bioseq_list, FindShortContigsCallback);
  }
  
  if (bioseq_list != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_SHORT_CONTIG;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (bioseq_list));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = bioseq_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}


static void RemoveShortContigsWithoutAnnotation (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr vnp, entityIDList = NULL;
  BioseqPtr  bsp;
  SeqFeatPtr sfp;
  SeqMgrFeatContext context;
  CharPtr           txt;

  if (Message (MSG_OKC, "Are you sure you want to remove short contigs without annotation?") == ANS_CANCEL) {
    return;
  }

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_BIOSEQ) {
      bsp = (BioseqPtr) vnp->data.ptrvalue;
      if (bsp->annot == NULL) {
        sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
        if (sfp == NULL) {
          if (lip != NULL) {
            lip->data_in_log = TRUE;
            if (lip->fp != NULL) {
              txt = GetDiscrepancyItemText (vnp);
              fprintf (lip->fp, "Removed short contig without annotation: %s\n", txt);
              txt = MemFree (txt);
            }
          }
          bsp->idx.deleteme = TRUE;
          ValNodeAddInt (&entityIDList, 0, bsp->idx.entityID);
        }
      }
    }
  }

  entityIDList = ValNodeSort (entityIDList, SortByIntvalue);
  ValNodeUnique (&entityIDList, SortByIntvalue, ValNodeFree);
  
  for (vnp = entityIDList; vnp != NULL; vnp = vnp->next) {
    DeleteMarkedObjects (vnp->data.intvalue, 0, NULL);
    ObjMgrSetDirtyFlag (vnp->data.intvalue, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, vnp->data.intvalue, 0, 0);    
  }
  entityIDList = ValNodeFree (entityIDList);
}


static void FindShortSequencesCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR bioseq_list;
  BioseqSetPtr    bssp;
  
  if (bsp == NULL || !ISA_na (bsp->mol) || userdata == NULL || bsp->length >= 50
      || IsmRNASequenceInGenProdSet (bsp))
  {
    return;
  }
  
  if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) bsp->idx.parentptr;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_parts) {
      return;
    }
  }
  
  bioseq_list = (ValNodePtr PNTR) userdata;
  
  ValNodeAddPointer (bioseq_list, OBJ_BIOSEQ, bsp);
}

extern void FindShortSequences (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip;
  CharPtr            bad_fmt = "%d sequences are shorter than 50 nt.";
  ValNodePtr         bioseq_list = NULL, vnp;
  
  if (discrepancy_list == NULL) return;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &bioseq_list, FindShortSequencesCallback);
  }
  
  if (bioseq_list != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_SHORT_CONTIG;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (bioseq_list));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = bioseq_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}


static void FindShortProtSequencesCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR bioseq_list;
  BioseqSetPtr    bssp;
  SeqDescrPtr     sdp;
  SeqMgrDescContext context;
  MolInfoPtr        mip;
  
  if (bsp == NULL || !ISA_aa (bsp->mol) || userdata == NULL || bsp->length >= 50)
  {
    return;
  }
  
  if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) bsp->idx.parentptr;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_parts) {
      return;
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
  if (sdp != NULL && (mip = (MolInfoPtr) sdp->data.ptrvalue) != NULL
      && mip->completeness != 1) {
    return;
  }
  
  bioseq_list = (ValNodePtr PNTR) userdata;
  
  ValNodeAddPointer (bioseq_list, OBJ_BIOSEQ, bsp);
}


static void FindShortProtSequences (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip;
  CharPtr            bad_fmt = "%d protein sequences are shorter than 50 aa.";
  ValNodePtr         bioseq_list = NULL, vnp;
  
  if (discrepancy_list == NULL) return;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &bioseq_list, FindShortProtSequencesCallback);
  }
  
  if (bioseq_list != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = SHORT_PROT_SEQUENCES;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (bioseq_list));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = bioseq_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}


typedef struct sdpandbsp {
  SeqDescrPtr sdp;
  BioseqPtr bsp;
} SdpAndBspData, PNTR SdpAndBspPtr;

typedef struct biosrccheck 
{
  BioSourcePtr biop;
  ValNodePtr   sdp_list;
} BioSrcCheckData, PNTR BioSrcCheckPtr;

static ValNodePtr FreeBioSrcCheckList (ValNodePtr biosrc_list)
{
  BioSrcCheckPtr  bscp;
  
  if (biosrc_list == NULL)
  {
    return NULL;
  }
  
  biosrc_list->next = FreeBioSrcCheckList (biosrc_list->next);
  
  bscp = (BioSrcCheckPtr) biosrc_list->data.ptrvalue;
  if (bscp != NULL)
  {
    bscp->sdp_list = ValNodeFreeData (bscp->sdp_list);
    bscp = MemFree (bscp);
  }
  biosrc_list = ValNodeFree (biosrc_list);
  return NULL;
}


static void FindInconsistentSourcesCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR biosrc_list, vnp;
  SeqDescrPtr     sdp;
  BioSrcCheckPtr  bscp;
  Boolean         found = FALSE;
  SeqMgrDescContext context;
  SdpAndBspPtr      sabp;
  
  if (bsp == NULL || !ISA_na (bsp->mol) || userdata == NULL)
  {
    return;
  }
  
  biosrc_list = (ValNodePtr PNTR) userdata;
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp != NULL)
  {
    sabp = (SdpAndBspPtr) MemNew (sizeof (SdpAndBspData));
    sabp->sdp = sdp;
    sabp->bsp = bsp;

    for (vnp = *biosrc_list; vnp != NULL && !found; vnp = vnp->next)
    {
      bscp = (BioSrcCheckPtr) vnp->data.ptrvalue;
      if (bscp != NULL && BioSourceMatch (sdp->data.ptrvalue, bscp->biop))
      {
        ValNodeAddPointer (&(bscp->sdp_list), 0, sabp);
        found = TRUE;
      }
    }
    if (!found)
    {
      bscp = (BioSrcCheckPtr) MemNew (sizeof (BioSrcCheckData));
      if (bscp != NULL)
      {
        bscp->biop = sdp->data.ptrvalue;
        ValNodeAddPointer (&(bscp->sdp_list), 0, sabp);
        ValNodeAddPointer (biosrc_list, 0, bscp);
      }
    }
  }
}


static ClickableItemPtr InconsistentBiosrc (BioSrcCheckPtr bscp)
{
  ClickableItemPtr dip = NULL;
  CharPtr          bad_fmt = "%d contigs have identical sources that do not match another contig source.";
  ValNodePtr       vnp;
  SdpAndBspPtr     sabp;

  if (bscp == NULL || bscp->sdp_list == NULL)
  {
    return NULL;
  }
  
  dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  if (dip != NULL)
  {
    dip->clickable_item_type = DISC_INCONSISTENT_BIOSRC;
    dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
    sprintf (dip->description, bad_fmt, ValNodeLen (bscp->sdp_list));
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    for (vnp = bscp->sdp_list; vnp != NULL; vnp = vnp->next) {
      sabp = (SdpAndBspPtr) vnp->data.ptrvalue;
      ValNodeAddPointer (&(dip->item_list), OBJ_SEQDESC, sabp->sdp);
      ValNodeAddPointer (&(dip->item_list), OBJ_BIOSEQ, sabp->bsp);
    }
  }      
  return dip;
}


extern void FindNonmatchingContigSources (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip;
  ValNodePtr       biosrc_list = NULL, vnp, vnp_s, sub_list = NULL, item_list = NULL;
  BioSrcCheckPtr   bscp;
  SdpAndBspPtr     sabp;
  CharPtr          disc_fmt = "%d inconsistent contig sources";

  if (discrepancy_list == NULL) return;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &biosrc_list, FindInconsistentSourcesCallback);
  }
  
  if (biosrc_list != NULL && biosrc_list->next != NULL)
  {
    for (vnp = biosrc_list; vnp != NULL; vnp = vnp->next)
    {
      bscp = (BioSrcCheckPtr) vnp->data.ptrvalue;
      dip = InconsistentBiosrc (bscp);
      ValNodeAddPointer (&sub_list, 0, dip);

      /* add sdp and bsp to item list */
      for (vnp_s = bscp->sdp_list; vnp_s != NULL; vnp_s = vnp_s->next) {
        sabp = (SdpAndBspPtr) vnp_s->data.ptrvalue;
        ValNodeAddPointer (&item_list, OBJ_SEQDESC, sabp->sdp);
        ValNodeAddPointer (&item_list, OBJ_BIOSEQ, sabp->bsp);
      }
    }
  }
  biosrc_list = FreeBioSrcCheckList (biosrc_list);
  if (item_list != NULL) {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (dip, 0, sizeof (ClickableItemData));
    dip->clickable_item_type = DISC_INCONSISTENT_BIOSRC;
    dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (disc_fmt) + 15));
    sprintf (dip->description, disc_fmt, ValNodeLen (item_list) / 2);
    dip->item_list = item_list;
    dip->subcategories = sub_list;
    ValNodeAddPointer (discrepancy_list, 0, dip);
  }
}


static Boolean IsWordChar (Char ch)
{
  if (isalpha (ch) || isdigit (ch)) {
    return TRUE;
  } else {
    return FALSE;
  }
}


NLM_EXTERN Boolean DoesStringContainPhrase (CharPtr str, CharPtr phrase, Boolean case_sensitive, Boolean whole_word)
{
  CharPtr cp;
  Boolean rval = FALSE;
  Int4    len;

  if (StringHasNoText (str) || StringHasNoText (phrase)) {
    return FALSE;
  }

  if (case_sensitive) {
    cp = StringSearch (str, phrase);
  } else {
    cp = StringISearch (str, phrase);
  }

  if (cp != NULL) {
    if (whole_word) {
      while (cp != NULL && !rval) {
        len = StringLen (phrase);
        if ((cp == str || !IsWordChar (*(cp - 1)))
            && (cp [len] == 0 || !IsWordChar (cp [len]))) {
          rval = TRUE;
        } else {
          if (case_sensitive) {
            cp = StringSearch (cp + 1, phrase);
          } else {
            cp = StringISearch (cp + 1, phrase);
          }
        }
      }
    } else {
      rval = TRUE;
    }
  }
  return rval;
}

typedef Boolean (*SuspectProductNameSearchFunc) PROTO ((CharPtr, CharPtr));
typedef void (*SuspectProductNameReplaceFunc) PROTO ((CharPtr PNTR, CharPtr, CharPtr, SeqFeatPtr));

typedef enum {
  eSuspectNameType_None = 0,
  eSuspectNameType_Typo = 1,
  eSuspectNameType_QuickFix,
  eSuspectNameType_NoOrganelleForProkaryote,
  eSuspectNameType_MightBeNonfunctional,
  eSuspectNameType_Database,
  eSuspectNameType_RemoveOrganismName,
  eSuspectNameType_InappropriateSymbol,
  eSuspectNameType_EvolutionaryRelationship,
  eSuspectNameType_UseProtein,
  eSuspectNameType_Max
} ESuspectNameType;

static CharPtr suspect_name_category_names[] = {
  "Unknown category",
  "Typo",
  "Quick fix",
  "Organelles not appropriate in prokaryote",
  "Suspicious phrase; should this be nonfunctional?",
  "May contain database identifier more appropriate in note; remove from product name",
  "Remove organism from product name",
  "Possible parsing error or incorrect formatting; remove inappropriate symbols",
  "Implies evolutionary relationship; change to -like protein",
  "Add protein to the end of product name",
  "Unknown category"
};


static Boolean CategoryOkForBioSource (BioSourcePtr biop, ESuspectNameType name_type)
{
  if (name_type != eSuspectNameType_NoOrganelleForProkaryote) {
    return TRUE;
  } else if (!HasTaxonomyID (biop)) {
    return TRUE;
  } else if (IsEukaryoticBioSource(biop)) {
    return FALSE;
  } else {
    return TRUE;
  }
}


typedef struct suspectproductname {
  CharPtr pattern;
  SuspectProductNameSearchFunc search_func;
  ESuspectNameType fix_type;
  CharPtr replace_phrase;
  SuspectProductNameReplaceFunc replace_func;
} SuspectProductNameData, PNTR SuspectProductNamePtr;


static Boolean EndsWithPattern (CharPtr pattern, CharPtr search)
{
  Int4 phrase_len, len;
  
  phrase_len = StringLen (pattern);
  len = StringLen (search);

  if (len >= phrase_len && StringICmp (search + len - phrase_len, pattern) == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean StartsWithPattern (CharPtr pattern, CharPtr search)
{
  Int4 phrase_len, len;
  
  phrase_len = StringLen (pattern);
  len = StringLen (search);

  if (len >= phrase_len && StringNICmp (search, pattern, phrase_len) == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static CharPtr s_putative_replacements[] = {
  "possible", 
  "potential",
  "predicted",
  "probable",
  NULL
};


static Boolean StartsWithPutativeReplacement (CharPtr pattern, CharPtr search)
{
  Int4 i;

  for (i = 0; s_putative_replacements[i] != NULL; i++) {
    if (StartsWithPattern(s_putative_replacements[i], search)) {
      return TRUE;
    }
  }
  return FALSE;
}


static Boolean MayContainPlural (CharPtr pattern, CharPtr search)
{
  return StringMayContainPlural (search);
}



static Boolean ContainsBracketsOrParentheses (CharPtr pattern, CharPtr search)
{
  return ContainsNorMoreSetsOfBracketsOrParentheses (search, 1);
}


static Boolean ContainsTwoSetsOfBracketsOrParentheses (CharPtr pattern, CharPtr search)
{
  return ContainsNorMoreSetsOfBracketsOrParentheses (search, 2);
}


static Boolean EndsWithPunct (CharPtr pattern, CharPtr search)
{
  Int4 len;
  Char last_ch;

  len = StringLen (search);
  last_ch = search[len - 1];
  if (last_ch == '.' || last_ch == ',' || last_ch == '-'
    || last_ch == '_' || last_ch == ':' || last_ch == '/')
  {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean BeginsWithPunct (CharPtr pattern, CharPtr search)
{
  if (search == NULL) return FALSE;
  if (search[0] == '.' || search[0] == ',' || search[0] == '-'
      || search[0] == '_' || search[0] == ':' || search[0] == '/') {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean BeginsOrEndsWithQuotes (CharPtr pattern, CharPtr search)
{
  Int4 len;

  if (search == NULL) return FALSE;
  if (search[0] == '\'' || search[0] == '"') {
    return TRUE;
  } else {
    len = StringLen (search);
    if (search[len - 1] == '\'' || search[len - 1] == '"') {
      return TRUE;
    } else {
      return FALSE;
    }
  }
}


static Boolean ContainsUnknownName (CharPtr pattern, CharPtr search)
{
  if (StringISearch(search, pattern) != NULL
      && StringISearch (search, "protein of unknown function") == NULL
      && StringISearch (search, "domain of unknown function") == NULL) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean ContainsWholeWordCaseSensitive (CharPtr pattern, CharPtr search)
{
  if (DoesStringContainPhrase (search, pattern, TRUE, TRUE)) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean ContainsWholeWord (CharPtr pattern, CharPtr search)
{
  if (DoesStringContainPhrase(search, pattern, FALSE, TRUE)) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean IsSingleWord (CharPtr pattern, CharPtr search)
{
  if (StringICmp (search, pattern) == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static CharPtr s_weasels[] = {
  "hypothetical", 
  "probable", 
  "putative",
  NULL
};


static Boolean IsSingleWordOrWeaselPlusSingleWord (CharPtr pattern, CharPtr search)
{
  Int4 i, len;

  if (StringICmp (search, pattern) == 0) {
    return TRUE;
  } else {
    for (i = 0; s_weasels[i] != NULL; i++) {
      len = StringLen (s_weasels[i]);
      if (StringNICmp (search, s_weasels[i], len) == 0
          && StringCmp (search + len + StringSpn (search + len, " "), pattern) == 0) {
        return TRUE;
      }
    }
    return FALSE;
  }
}


static Boolean NormalSearch (CharPtr pattern, CharPtr search)
{
  if (StringISearch(search, pattern) != NULL) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean ThreeOrMoreNumbersTogether (CharPtr pattern, CharPtr search)
{
  return ContainsThreeOrMoreNumbersTogether (search);
}


static Boolean ContainsUnderscore (CharPtr pattern, CharPtr search)
{
  return StringContainsUnderscore (search);
}


static Boolean PrefixPlusNumbersOnly (CharPtr pattern, CharPtr search)
{
  return IsPrefixPlusNumbers (pattern, search);
}


static Boolean EndsWithFold (CharPtr pattern, CharPtr search)
{
  Int4 len;

  if (search == NULL) {
    return FALSE;
  }
  len = StringLen (search);
  if (len < 4) {
    return FALSE;
  }
  if (StringICmp (search + len - 4, "fold") == 0) {
    if (StringCmp (search + len - 4, "folD") == 0
        || StringCmp (search + len - 4, "FolD") == 0) {
      return FALSE;
    } else {
      return TRUE;
    }
  } else {
    return FALSE;
  }
}

static Boolean AllCapitalLetters (CharPtr pattern, CharPtr search)
{
  CharPtr cp;
  Boolean any_alpha = FALSE;

  if (search == NULL) {
    return FALSE;
  }
  cp = search;
  while (*cp != 0) {
    if (isalpha (*cp)) {
      any_alpha = TRUE;
      if (islower (*cp)) {
        return FALSE;
      }
    }
    ++cp;
  }
  return any_alpha;
}


static Boolean ContainsUnbalancedParentheses (CharPtr pattern, CharPtr search)
{
  return StringContainsUnbalancedParentheses (search);
}


static Boolean IsTooLong (CharPtr pattern, CharPtr search)
{
  if (StringISearch (search, "bifunctional") != NULL 
    || StringISearch (search, "multifunctional") != NULL) {
    return FALSE;
  } else if (StringLen (search) > 100) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static CharPtr s_pseudoweasels[] = {
  "pseudouridine",
  "pseudoazurin",
  "pseudouridylate",
  NULL};

static Boolean ContainsPseudo (CharPtr pattern, CharPtr search)
{
  CharPtr cp;
  Int4    i, len;
  Boolean bad_pseudo;

  if (search == NULL) {
    return FALSE;
  }
  cp = StringISearch (search, "pseudo");
  while (cp != NULL) {
    bad_pseudo = FALSE;
    for (i = 0; s_pseudoweasels[i] != NULL && !bad_pseudo; i++) {
      len = StringLen (s_pseudoweasels[i]);
      if (StringNCmp (cp, s_pseudoweasels[i], len) == 0) {
        bad_pseudo = TRUE;
        cp = StringISearch (cp + len, "pseudo");
      }
    }
    if (!bad_pseudo) {
      return TRUE;
    }
  }
  return FALSE;
}


static Boolean ContainsDoubleSpace (CharPtr pattern, CharPtr search)
{
  if (search == NULL) {
    return FALSE;
  }
  if (StringSearch (search, "  ") != NULL) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void SimpleReplaceFunc (CharPtr PNTR orig, CharPtr find, CharPtr replace, SeqFeatPtr sfp)
{
  FindReplaceString (orig, find, replace, FALSE, TRUE);
}


static void SimpleReplaceAnywhereFunc (CharPtr PNTR orig, CharPtr find, CharPtr replace, SeqFeatPtr sfp)
{
  FindReplaceString (orig, find, replace, FALSE, FALSE);
}


static void ReplaceWholeNameFunc (CharPtr PNTR orig, CharPtr find, CharPtr replace, SeqFeatPtr sfp)
{
  if (orig == NULL) {
    return;
  }
  if (IsSingleWordOrWeaselPlusSingleWord(find, *orig)) {
    *orig = MemFree (*orig);
    *orig = StringSave (replace);
  }
}


static void ReplaceWholeNameAddNoteFunc (CharPtr PNTR orig, CharPtr find, CharPtr replace, SeqFeatPtr sfp)
{
  if (orig == NULL) {
    return;
  }
  if (IsSingleWordOrWeaselPlusSingleWord(find, *orig)) {
    SetStringValue (&(sfp->comment), *orig, ExistingTextOption_append_semi);
    *orig = MemFree (*orig);
    *orig = StringSave (replace);
  }
}


static void ReplaceAtFront (CharPtr PNTR orig, CharPtr find, CharPtr replace, SeqFeatPtr sfp)
{
  Int4 orig_len, find_len, replace_len, new_len;
  CharPtr new_str;

  if (orig == NULL || find == NULL) {
    return;
  }

  orig_len = StringLen (*orig);
  find_len = StringLen (find);
  if (find_len > orig_len || StringNICmp (*orig, find, find_len) != 0) {
    return;
  }
  replace_len = StringLen (replace);

  new_len = orig_len + replace_len - find_len;
  new_str = (CharPtr) MemNew (sizeof (Char) * (new_len + 1));
  if (replace_len > 0) {
    StringCpy (new_str, replace);
  }
  StringCat (new_str, (*orig) + find_len);
  *orig = MemFree (*orig);
  *orig = new_str;
}


static void ReplaceAtEnd (CharPtr PNTR orig, CharPtr find, CharPtr replace, SeqFeatPtr sfp)
{
  Int4 orig_len, find_len, replace_len, new_len;
  CharPtr new_str;

  if (orig == NULL || find == NULL) {
    return;
  }

  orig_len = StringLen (*orig);
  find_len = StringLen (find);
  if (find_len > orig_len || StringICmp ((*orig) + orig_len - find_len, find) != 0) {
    return;
  }
  replace_len = StringLen (replace);

  new_len = orig_len + replace_len - find_len;
  new_str = (CharPtr) MemNew (sizeof (Char) * (new_len + 1));
  StringNCpy (new_str, *orig, orig_len - find_len);
  if (replace_len > 0) {
    StringCat (new_str, replace);
  }
  *(new_str + new_len) = 0;
  *orig = MemFree (*orig);
  *orig = new_str;
}


static void UsePutative (CharPtr PNTR orig, CharPtr find, CharPtr replace, SeqFeatPtr sfp)
{
  Int4 i;
  for (i = 0; s_putative_replacements[i] != NULL; i++) {
    ReplaceAtFront (orig, s_putative_replacements[i], "putative", sfp);
  }
}


static void RemoveBeginningAndEndingQuotes (CharPtr PNTR orig, CharPtr find, CharPtr replace, SeqFeatPtr sfp)
{
  CharPtr src, dst;
  Int4 len;

  if (orig == NULL || *orig == NULL || !BeginsOrEndsWithQuotes (NULL, *orig)) {
    return;
  }
  src = *orig;
  dst = *orig;
  if (*src == '\'' || *src == '"') {
    src++;
    while (*src != 0) {
      *dst = *src;
      dst++;
      src++;
    }
    *dst = 0;
  }
  len = StringLen (*orig);
  if ((*orig)[len - 1] == '\'' || (*orig)[len - 1] == '"') {
    (*orig)[len - 1] = 0;
  }
}


static void FixLongProduct (CharPtr PNTR orig, CharPtr find, CharPtr replace, SeqFeatPtr sfp)
{
  Int4 len, keep_len;
  if (orig == NULL || *orig == NULL || sfp == NULL || *orig == sfp->comment) {
    return;
  }
  len = StringLen (*orig);
  keep_len = StringCSpn (*orig, ",;(");
  if (keep_len < len) {
    SetStringValue (&(sfp->comment), *orig, ExistingTextOption_append_semi);
    *((*orig) + keep_len) = 0;
  }
}


static void HaemReplaceFunc (CharPtr PNTR orig, CharPtr find, CharPtr replace, SeqFeatPtr sfp)
{
  if (orig == NULL || *orig == NULL) {
    return;
  }

  FindReplaceString (orig, find, "heme", FALSE, TRUE);
  FindReplaceString (orig, find, "hem", FALSE, FALSE);
}


static CharPtr SummarizeSuspectPhraseFunc (SuspectProductNameSearchFunc s)
{
  if (s == NULL) {
    return "NULL function";
  } else  if (s == EndsWithPattern) {
    return "occurs at end of text";
  } else if (s == ContainsWholeWord) {
    return "contains phrase as whole word";
  } else if (s == StartsWithPattern) {
    return "occurs at beginning of text";
  } else if (s == ContainsWholeWordCaseSensitive) {
    return "contains phrase as whole word, case sensitive";
  } else if (s == IsSingleWord) {
    return "entire text matches (not case sensitive)";
  } else if (s == IsSingleWordOrWeaselPlusSingleWord) {
    return "entire text matches (not case sensitive) or text matches after weasel word";
  } else if (s == NormalSearch) {
    return "contains phrase anywhere, not case sensitive";
  } else if (s == ContainsDoubleSpace) {
    return "contains double space";
  } else if (s == PrefixPlusNumbersOnly) {
    return "entire product is prefix followed by numbers";
  } else if (s == IsTooLong) {
    return "longer than 50 characters";
  } else {
    return "special rules";
  }
}


static CharPtr SummarizeSuspectReplacementPhrase (SuspectProductNameReplaceFunc s, CharPtr replace_phrase)
{
  CharPtr phrase = NULL;
  CharPtr simple_fmt = "Replace with '%s' (whole word)";
  CharPtr simple_anywhere_fmt = "Replace with '%s'";
  CharPtr whole_fmt = "Replace entire product name with '%s'";
  CharPtr whole_note_fmt = "Move product name to note, use '%s' for product name";

  
  if (s == NULL) {
    return StringSave ("No replacement");
  } else if (s == SimpleReplaceFunc) {
    phrase = (CharPtr) MemNew (sizeof (Char) * (StringLen (simple_fmt) + StringLen (replace_phrase)));
    sprintf (phrase, simple_fmt, replace_phrase);
  } else if (s == SimpleReplaceAnywhereFunc) {
    phrase = (CharPtr) MemNew (sizeof (Char) * (StringLen (simple_anywhere_fmt) + StringLen (replace_phrase)));
    sprintf (phrase, simple_anywhere_fmt, replace_phrase);
  } else if (s == FixLongProduct) {
    phrase = StringSave ("Truncate at first comma or semicolon");
  } else if (s == UsePutative) {
    phrase = StringSave ("Replace with 'putative'");
  } else if (s == ReplaceWholeNameFunc) {
    phrase = (CharPtr) MemNew (sizeof (Char) * (StringLen (whole_fmt) + StringLen (replace_phrase)));
    sprintf (phrase, whole_fmt, replace_phrase);
  } else if (s == ReplaceWholeNameAddNoteFunc) {
    phrase = (CharPtr) MemNew (sizeof (Char) * (StringLen (whole_note_fmt) + StringLen (replace_phrase)));
    sprintf (phrase, whole_note_fmt, replace_phrase);
  } else if (s == ReplaceAtEnd || s == ReplaceAtFront) {
    phrase = (CharPtr) MemNew (sizeof (Char) * (StringLen (simple_anywhere_fmt) + StringLen (replace_phrase)));
    sprintf (phrase, simple_anywhere_fmt, replace_phrase);
  } else {
    phrase = StringSave ("Unknown replacement action");
  }
  return phrase;
}


static SuspectProductNameData suspect_product_terms[] = {
  { "beginning with period, comma, or hyphen" , BeginsWithPunct, eSuspectNameType_InappropriateSymbol, NULL, NULL } ,
  { "begins or ends with quotes", BeginsOrEndsWithQuotes, eSuspectNameType_QuickFix, NULL, RemoveBeginningAndEndingQuotes } ,
  { "binding" , EndsWithPattern, eSuspectNameType_UseProtein, NULL, NULL } ,
  { "domain", EndsWithPattern, eSuspectNameType_UseProtein, NULL, NULL } ,
  { "like" , EndsWithPattern, eSuspectNameType_UseProtein, NULL, NULL } ,
  { "motif" , EndsWithPattern, eSuspectNameType_UseProtein, NULL, NULL } ,
  { "related" , EndsWithPattern, eSuspectNameType_UseProtein, NULL, NULL } ,
  { "repeat", EndsWithPattern, eSuspectNameType_UseProtein, NULL, NULL } ,
  { "fold" , EndsWithFold, eSuspectNameType_UseProtein, NULL, NULL } ,
  { "Arabidopsis" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Aspergillus" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "B.subtilis" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Bacillus" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Bacteroides" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Campylobacter" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Chlamydial" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Chlamydomonas" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Drosophila" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "E.coli" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Escherichia" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Helicobacter" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Includes:" , ContainsWholeWord, eSuspectNameType_InappropriateSymbol, NULL, NULL } ,
  { "Jejuni" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Leishmania" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Marinococcus" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Mus musculus" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Mycobacterium" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Pestis" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Rhodobacter" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Salmonella" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Staphlococcal" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Staphlococcus" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Staphylococcus" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Streptococcus" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Subtilis" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Tuberculosis" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Typhimurium" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Yersinia" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "aminotransferasee" , ContainsWholeWord, eSuspectNameType_Typo , "aminotransferase", SimpleReplaceFunc } ,
  { "arginin " , ContainsWholeWord, eSuspectNameType_Typo , "arginine ", SimpleReplaceFunc } ,
  { "argininte" , ContainsWholeWord, eSuspectNameType_Typo , "arginine", SimpleReplaceFunc } ,
  { "aureus" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "bioin" , ContainsWholeWord, eSuspectNameType_Typo , "biotin", SimpleReplaceFunc } ,
  { "biosythesis" , ContainsWholeWord, eSuspectNameType_Typo , "biosynthesis", SimpleReplaceFunc } ,
  { "cerevisiae" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "chelatin" , ContainsWholeWord, eSuspectNameType_Typo , "chelating", SimpleReplaceFunc } ,
  { "coli" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "contain" , ContainsWholeWord, eSuspectNameType_None, NULL, NULL } ,
  { "deydrogenase" , ContainsWholeWord, eSuspectNameType_Typo, "dehydrogenase", SimpleReplaceFunc } ,
  { "diacyglycerol" , ContainsWholeWord, eSuspectNameType_Typo, "diacylglycerol", SimpleReplaceFunc } ,
  { "domainl", ContainsWholeWord, eSuspectNameType_Typo, "domain", SimpleReplaceFunc } ,
  { "enterica" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "exporte" , ContainsWholeWord, eSuspectNameType_Typo, "exported", SimpleReplaceFunc } ,
  { "familie" , ContainsWholeWord, eSuspectNameType_Typo, "family", SimpleReplaceFunc } ,
  { "gene" , ContainsWholeWord, eSuspectNameType_None, NULL, NULL } ,
  { "genes" , ContainsWholeWord, eSuspectNameType_None, NULL, NULL } ,
  { "glycin" , ContainsWholeWord, eSuspectNameType_Typo, "glycine", SimpleReplaceFunc } ,
  { "glycosy" , ContainsWholeWord, eSuspectNameType_Typo, "glucosyl", SimpleReplaceFunc } ,
  { "halophilus" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "hemaggltinin" , ContainsWholeWord, eSuspectNameType_Typo, "hemagglutinin", SimpleReplaceFunc } ,
  { "hexpeptide" , ContainsWholeWord, eSuspectNameType_Typo, "hexapeptide", SimpleReplaceFunc } ,
  { "histide" , ContainsWholeWord, eSuspectNameType_Typo, "histidine", SimpleReplaceFunc } ,
  { "homo" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "homocystein" , ContainsWholeWord, eSuspectNameType_Typo, "homocysteine", SimpleReplaceFunc } ,
  { "hyp domain protein" , IsSingleWord, eSuspectNameType_Typo, "hypothetical protein", SimpleReplaceFunc },
  { "hypot" , ContainsWholeWord, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "hypothe" , ContainsWholeWord, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "hypothet" , ContainsWholeWord, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "hypothetic" , ContainsWholeWord, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "hypothetica" , ContainsWholeWord, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "hypothetical domain protein" , IsSingleWord, eSuspectNameType_Typo, "hypothetical protein", SimpleReplaceFunc },
  { "inactivated derivative" , ContainsWholeWord, eSuspectNameType_None, NULL, NULL } ,
  { "initation" , ContainsWholeWord, eSuspectNameType_Typo, "initiation", SimpleReplaceFunc } ,
  { "invertion" , ContainsWholeWord, eSuspectNameType_Typo, "inversion", SimpleReplaceFunc } ,
  { "isomaerase" , ContainsWholeWord, eSuspectNameType_Typo, "isomerase", SimpleReplaceFunc } ,
  { "mobilisation" , ContainsWholeWord, eSuspectNameType_Typo, "mobilization", SimpleReplaceFunc } ,
  { "mouse" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "mutatrotase" , ContainsWholeWord, eSuspectNameType_Typo, "mutarotase", SimpleReplaceFunc } ,
  { "ncharacterized" , ContainsWholeWord, eSuspectNameType_Typo, "uncharacterized", SimpleReplaceFunc } ,
  { "ndoribonuclease" , ContainsWholeWord, eSuspectNameType_Typo, "endoribonuclease", SimpleReplaceFunc } ,
  { "niger" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "ntegral " , ContainsWholeWord, eSuspectNameType_Typo, "integral ", SimpleReplaceFunc } ,
  { "obalt" , ContainsWholeWord, eSuspectNameType_Typo, "cobalt", SimpleReplaceFunc } ,
  { "odule" , ContainsWholeWord, eSuspectNameType_None, NULL, NULL } ,
  { "orf, hyp" , IsSingleWord, eSuspectNameType_Typo, "hypothetical protein", SimpleReplaceFunc },
  { "orf, hypothetical" , IsSingleWord, eSuspectNameType_Typo, "hypothetical protein", SimpleReplaceFunc },
  { "oxidoreductasee" , ContainsWholeWord, eSuspectNameType_Typo, "oxidoreductase", SimpleReplaceFunc } ,
  { "oxidoredutase" , ContainsWholeWord, eSuspectNameType_Typo, "oxidoreductase", SimpleReplaceFunc } ,
  { "periplamic" , ContainsWholeWord, eSuspectNameType_Typo, "periplasmic", SimpleReplaceFunc } ,
  { "periplasmc" , ContainsWholeWord, eSuspectNameType_Typo, "periplasmic", SimpleReplaceFunc } ,
  { "phosphatidyltransferse" , ContainsWholeWord, eSuspectNameType_Typo, "phosphatidyltransferase", SimpleReplaceFunc } ,
  { "phosphopantethiene" , ContainsWholeWord, eSuspectNameType_Typo, "phosphopantetheine", SimpleReplaceFunc } ,
  { "pombe" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "portein" , ContainsWholeWord, eSuspectNameType_Typo, "protein", SimpleReplaceFunc } ,
  { "protei" , ContainsWholeWord, eSuspectNameType_Typo, "protein", SimpleReplaceFunc } ,
  { "protwin" , ContainsWholeWord, eSuspectNameType_Typo, "protein", SimpleReplaceFunc } ,
  { "pseudo" , ContainsWholeWord, eSuspectNameType_MightBeNonfunctional, NULL, NULL } ,
  { "pseudomonas" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "puter" , ContainsWholeWord, eSuspectNameType_Typo, "outer", SimpleReplaceFunc } ,
  { "pylori" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "rat" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "reductasee" , ContainsWholeWord, eSuspectNameType_Typo, "reductase", SimpleReplaceFunc } ,
  { "rsponse" , ContainsWholeWord, eSuspectNameType_Typo, "response", SimpleReplaceFunc } ,
  { "serovar" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "sigm" , ContainsWholeWord, eSuspectNameType_Typo, "sigma", NULL } ,
  { "sreptomyces" , ContainsWholeWord, eSuspectNameType_None, NULL, NULL } ,
  { "staphylococcal" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "start codon" , ContainsWholeWord, eSuspectNameType_MightBeNonfunctional, NULL, NULL } ,
  { "streptococcal" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "streptomyces" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "subsp" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "tetracenpmycin" , ContainsWholeWord, eSuspectNameType_Typo, "tetracenomycin", SimpleReplaceFunc } ,
  { "thaliana" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "thiamin/thiamin" , ContainsWholeWord, eSuspectNameType_Typo, "thiamin/thiamine", SimpleReplaceFunc } ,
  { "thioderoxin" , ContainsWholeWord, eSuspectNameType_Typo, "thioredoxin", SimpleReplaceFunc } ,
  { "threonin" , ContainsWholeWord, eSuspectNameType_Typo, "threonine", SimpleReplaceFunc } ,
  { "transcrIptional" , ContainsWholeWordCaseSensitive, eSuspectNameType_Typo, "transcriptional", SimpleReplaceFunc } ,
  { "transemembrane" , ContainsWholeWord, eSuspectNameType_Typo, "transmembrane", SimpleReplaceFunc } ,
  { "transferasee" , ContainsWholeWord, eSuspectNameType_Typo, "transferase", SimpleReplaceFunc } ,
  { "transmebrane" , ContainsWholeWord, eSuspectNameType_Typo, "transmembrane", SimpleReplaceFunc } ,
  { "unkn", IsSingleWord, eSuspectNameType_None, "hypothetical protein", SimpleReplaceFunc },
  { "unnamed" , ContainsWholeWord, eSuspectNameType_None, NULL, NULL } ,
  { "utilisation" , ContainsWholeWord, eSuspectNameType_Typo, "utilization", SimpleReplaceFunc } ,
  { "xenopus" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "yeast" , ContainsWholeWord, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "ypothetical" , ContainsWholeWord, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "ytochrome" , ContainsWholeWord, eSuspectNameType_Typo, "cytochrome", SimpleReplaceFunc } ,
  { "containing" , StartsWithPattern, eSuspectNameType_None, NULL, NULL } ,
  { "from" , StartsWithPattern, eSuspectNameType_None, NULL, NULL } ,
  { "CHC2 zinc finger" , IsSingleWord, eSuspectNameType_UseProtein, NULL, NULL } ,
  { "SWIM zinc finger" , IsSingleWord, eSuspectNameType_UseProtein, NULL, NULL } ,
  { "probable protein" , IsSingleWord, eSuspectNameType_None, NULL, NULL } ,
  { "protein" , IsSingleWord, eSuspectNameType_None, NULL, NULL } ,
  { "sodium" , IsSingleWord, eSuspectNameType_None, NULL, NULL } ,
  { "IS" , PrefixPlusNumbersOnly, eSuspectNameType_None, NULL, NULL } ,
  { "three or more numbers together, not after 'UPF' or 'DUF' or 'IS' and not followed by the word 'family' and not preceded by either 'cytochrome' or 'coenzyme'" , ThreeOrMoreNumbersTogether,
 eSuspectNameType_Database, NULL, NULL } ,
  { "all capital letters" , AllCapitalLetters, eSuspectNameType_InappropriateSymbol, NULL, NULL } ,
  { "#" , NormalSearch, eSuspectNameType_InappropriateSymbol, NULL, NULL } ,
  { ". " , NormalSearch, eSuspectNameType_InappropriateSymbol, NULL, NULL } ,
  { "=" , NormalSearch, eSuspectNameType_InappropriateSymbol, NULL, NULL } ,
  { "?" , NormalSearch, eSuspectNameType_InappropriateSymbol, NULL, NULL } ,
  { "%" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "Chloroplast" , NormalSearch, eSuspectNameType_NoOrganelleForProkaryote, NULL, NULL } ,
  { "ECOLI" , NormalSearch, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Fragment" , NormalSearch, eSuspectNameType_MightBeNonfunctional, NULL, NULL } ,
  { "Frameshift" , NormalSearch, eSuspectNameType_MightBeNonfunctional, NULL, NULL } ,
  { "Homolog" , NormalSearch, eSuspectNameType_EvolutionaryRelationship, NULL, NULL } ,
  { "Intein" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "Intiation" , NormalSearch, eSuspectNameType_Typo, "initiation", SimpleReplaceFunc } ,
  { "K potassium" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "K+ potassium" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "Mitochondrial" , NormalSearch, eSuspectNameType_NoOrganelleForProkaryote, NULL, NULL } ,
  { "No definition line found" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "Plasmodium" , NormalSearch, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "Portein" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "Related to" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "Similar to" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "Transemembrane" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "Transmebrane" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "\\-PA" , NormalSearch, eSuspectNameType_InappropriateSymbol, NULL, NULL } ,
  { "accessroy" , NormalSearch, eSuspectNameType_Typo, "accessory", SimpleReplaceFunc } ,
  { "aceytltranferase" , NormalSearch, eSuspectNameType_Typo, "acetyltransferase", SimpleReplaceFunc } ,
  { "active site" , NormalSearch, eSuspectNameType_UseProtein, NULL, NULL } ,
  { "adenylattransferase" , NormalSearch, eSuspectNameType_Typo, "adenylate transferase", SimpleReplaceFunc } ,
  { "adenylytransferase" , NormalSearch, eSuspectNameType_Typo, "adenylyltransferase", SimpleReplaceFunc } ,
  { "alternate protein name" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "aluminium" , NormalSearch, eSuspectNameType_Typo, "aluminum", SimpleReplaceFunc } ,
  { "aminopetidase" , NormalSearch, eSuspectNameType_Typo, "aminopeptidase", SimpleReplaceFunc } ,
  { "aparaginase" , NormalSearch, eSuspectNameType_Typo, "asparaginase", SimpleReplaceFunc } ,
  { "asparate" , NormalSearch, eSuspectNameType_Typo, "aspartate", SimpleReplaceFunc } ,
  { "authentic point mutation" , NormalSearch, eSuspectNameType_MightBeNonfunctional, NULL, NULL } ,
  { "bifunctional" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "bifunctionnal" , NormalSearch, eSuspectNameType_Typo, "bifunctional", SimpleReplaceFunc } ,
  { "bigenesis" , NormalSearch, eSuspectNameType_Typo, "biogenesis", SimpleReplaceFunc } ,
  { "biosyntesis" , NormalSearch, eSuspectNameType_Typo, "biosynthesis", SimpleReplaceFunc } ,
  { "bnding" , NormalSearch, eSuspectNameType_Typo, "binding", SimpleReplaceFunc } ,
  { "bos taurus" , NormalSearch, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "carboxilic" , NormalSearch, eSuspectNameType_Typo, "carboxylic", SimpleReplaceFunc } ,
  { "cell divisionFtsK/SpoIIIE" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "characteris" , NormalSearch, eSuspectNameType_Typo, "characteriz", SimpleReplaceAnywhereFunc } ,
  { "coantaining" , NormalSearch, eSuspectNameType_Typo, "containing", SimpleReplaceFunc } ,
  { "coenzye" , NormalSearch, eSuspectNameType_Typo, "coenzyme", SimpleReplaceFunc } ,
  { "componenet" , NormalSearch, eSuspectNameType_Typo, "component", SimpleReplaceFunc } ,
  { "componnent" , NormalSearch, eSuspectNameType_Typo, "component", SimpleReplaceFunc } ,
  { "consevered" , NormalSearch, eSuspectNameType_Typo, "conserved", SimpleReplaceFunc } ,
  { "containg" , NormalSearch, eSuspectNameType_Typo, "containing", SimpleReplaceFunc } ,
  { "cotaining" , NormalSearch, eSuspectNameType_Typo, "containing", SimpleReplaceFunc } ,
  { "degration" , NormalSearch, eSuspectNameType_Typo, "degradation", SimpleReplaceFunc } ,
  { "deletion" , NormalSearch, eSuspectNameType_MightBeNonfunctional, NULL, NULL } ,
  { "dependant" , NormalSearch, eSuspectNameType_Typo, "dependent", SimpleReplaceFunc } ,
  { "dimerisation" , NormalSearch, eSuspectNameType_Typo, "dimerization", SimpleReplaceFunc } ,
  { "dimerising" , NormalSearch, eSuspectNameType_Typo, "dimerizing", SimpleReplaceFunc } ,
  { "dioxyenase" , NormalSearch, eSuspectNameType_Typo, "dioxygenase", SimpleReplaceFunc } ,
  { "disulphide" , NormalSearch, eSuspectNameType_Typo, "disulfide", SimpleReplaceFunc } ,
  { "divison" , NormalSearch, eSuspectNameType_Typo, "division", SimpleReplaceFunc } ,
  { "domain domain" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "domain protein domain protein" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "domian" , NormalSearch, eSuspectNameType_Typo, "domain", SimpleReplaceFunc } ,
  { "dyhydrogenase" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "dyhydrogenase" , NormalSearch, eSuspectNameType_Typo, "dehydrogenase", SimpleReplaceFunc } ,
  { "enentioselective" , NormalSearch, eSuspectNameType_Typo, "enantioselective", SimpleReplaceFunc } ,
  { "facotr" , NormalSearch, eSuspectNameType_Typo, "factor", SimpleReplaceFunc } ,
  { "fagella", NormalSearch, eSuspectNameType_Typo, "flagella", SimpleReplaceFunc } ,
  { "family family" , NormalSearch, eSuspectNameType_Typo, "family", SimpleReplaceFunc } ,
  { "flageller" , NormalSearch, eSuspectNameType_Typo, "flagellar", SimpleReplaceFunc } ,
  { "frame shift" , NormalSearch, eSuspectNameType_MightBeNonfunctional, NULL, NULL } ,
  { "gIycerol" , NormalSearch, eSuspectNameType_Typo, "glycerol", SimpleReplaceFunc } ,
  { "glcosyl" , NormalSearch, eSuspectNameType_Typo, "glycosyl", SimpleReplaceFunc } ,
  { "glucosainyl" , NormalSearch, eSuspectNameType_Typo, "glucosaminyl", SimpleReplaceFunc } ,
  { "glutaminne" , NormalSearch, eSuspectNameType_Typo, "glutamine", SimpleReplaceFunc } ,
  { "golgi" , NormalSearch, eSuspectNameType_NoOrganelleForProkaryote, NULL, NULL } ,
  { "haem" , NormalSearch, eSuspectNameType_Typo, "heme", HaemReplaceFunc } ,
  { "haemagglutination" , NormalSearch, eSuspectNameType_Typo, "hemagglutination", SimpleReplaceFunc } ,
  { "heam" , NormalSearch, eSuspectNameType_Typo, "heme", HaemReplaceFunc } ,
  { "hemelysin" , NormalSearch, eSuspectNameType_Typo, "hemolysin", SimpleReplaceFunc } ,
  { "hemoglobine" , NormalSearch, eSuspectNameType_Typo, "hemoglobin", SimpleReplaceFunc } ,
  { "hexapaptide" , NormalSearch, eSuspectNameType_Typo, "hexapeptide", SimpleReplaceFunc } ,
  { "highly conserved" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "histadine" , NormalSearch, eSuspectNameType_Typo, "histidine", SimpleReplaceFunc } ,
  { "homeserine" , NormalSearch, eSuspectNameType_Typo, "homoserine", SimpleReplaceFunc } ,
  { "homo sapiens" , NormalSearch, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "hpothetical" , NormalSearch, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "human" , NormalSearch, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "hyphotetical" , NormalSearch, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "hyphotheical" , NormalSearch, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "hypotehtical" , NormalSearch, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "hypotethical" , NormalSearch, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "hypotetical" , NormalSearch, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "hypotheical" , NormalSearch, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "hypotheitcal" , NormalSearch, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "hypothetcial" , NormalSearch, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "hypothteical" , NormalSearch, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "hypothtical" , NormalSearch, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "hypthetical" , NormalSearch, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "hyptothetical" , NormalSearch, eSuspectNameType_Typo, "hypothetical", SimpleReplaceFunc } ,
  { "inductible" , NormalSearch, eSuspectNameType_Typo, "inducible", SimpleReplaceFunc } ,
  { "interrupt" , NormalSearch, eSuspectNameType_MightBeNonfunctional, NULL, NULL } ,
  { "isomerse" , NormalSearch, eSuspectNameType_Typo, "isomerase", SimpleReplaceFunc } ,
  { "majour" , NormalSearch, eSuspectNameType_Typo, "major", SimpleReplaceFunc } ,
  { "mambrane" , NormalSearch, eSuspectNameType_Typo, "membrane", SimpleReplaceFunc } ,
  { "meausure" , NormalSearch, eSuspectNameType_Typo, "measure", SimpleReplaceFunc } ,
  { "membranne" , NormalSearch, eSuspectNameType_Typo, "membrane", SimpleReplaceFunc } ,
  { "methlytransferase" , NormalSearch, eSuspectNameType_Typo, "methyltransferase", SimpleReplaceFunc } ,
  { "metylase" , NormalSearch, eSuspectNameType_Typo, "methylase", SimpleReplaceFunc } ,
  { "molibdenum" , NormalSearch, eSuspectNameType_Typo, "molybdenum", SimpleReplaceFunc } ,
  { "molybopterin" , NormalSearch, eSuspectNameType_Typo, "molybdopterin", SimpleReplaceFunc } ,
  { "molydopterin" , NormalSearch, eSuspectNameType_Typo, "molybdopterin", SimpleReplaceFunc } ,
  { "monooxigenase" , NormalSearch, eSuspectNameType_Typo, "monooxygenase", SimpleReplaceFunc } ,
  { "monoxyde" , NormalSearch, eSuspectNameType_Typo, "monoxide", SimpleReplaceFunc } ,
  { "monoxygenase" , NormalSearch, eSuspectNameType_Typo, "monooxygenase", SimpleReplaceFunc } ,
  { "mulitdrug" , NormalSearch, eSuspectNameType_Typo, "multidrug", SimpleReplaceFunc } ,
  { "multifunctional", NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "narrowly conserved" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "nickle" , NormalSearch, eSuspectNameType_Typo, "nickel", SimpleReplaceFunc } ,
  { "novel protein" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "nucelar" , NormalSearch, eSuspectNameType_Typo, "nuclear", SimpleReplaceFunc } ,
  { "nucleotydyl" , NormalSearch, eSuspectNameType_Typo, "nucleotidyl", SimpleReplaceFunc } ,
  { "nulcear" , NormalSearch, eSuspectNameType_Typo, "nuclear", SimpleReplaceFunc } ,
  { "open reading frame" , NormalSearch, eSuspectNameType_MightBeNonfunctional, NULL, NULL } ,
  { "or related" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "orphan protein" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "ortholog" , NormalSearch, eSuspectNameType_EvolutionaryRelationship, NULL, NULL } ,
  { "outers" , NormalSearch, eSuspectNameType_Typo, "outer", SimpleReplaceFunc } ,
  { "oxidoreducatse" , NormalSearch, eSuspectNameType_Typo, "oxidoreductase", SimpleReplaceFunc } ,
  { "oxidoreductasse" , NormalSearch, eSuspectNameType_Typo, "oxidoreductase", SimpleReplaceFunc } ,
  { "oxidoreduxtase" , NormalSearch, eSuspectNameType_Typo, "oxidoreductase", SimpleReplaceFunc } ,
  { "oxydase" , NormalSearch, eSuspectNameType_Typo, "oxidase", SimpleReplaceFunc } ,
  { "paralog" , NormalSearch, eSuspectNameType_EvolutionaryRelationship, NULL, NULL } ,
  { "partial" , NormalSearch, eSuspectNameType_MightBeNonfunctional, NULL, NULL } ,
  { "peptidodoglycan" , NormalSearch, eSuspectNameType_Typo, "peptidoglycan", SimpleReplaceFunc } ,
  { "periplsmic" , NormalSearch, eSuspectNameType_Typo, "periplasmic", SimpleReplaceFunc } ,
  { "phophate" , NormalSearch, eSuspectNameType_Typo, "phosphate", SimpleReplaceFunc } ,
  { "phopho" , NormalSearch, eSuspectNameType_Typo, "phospho", SimpleReplaceFunc } ,
  { "phophoserine" , NormalSearch, eSuspectNameType_Typo, "phosphoserine", SimpleReplaceFunc } ,
  { "phoshate" , NormalSearch, eSuspectNameType_Typo, "phosphate", SimpleReplaceFunc } ,
  { "phosphatransferase" , NormalSearch, eSuspectNameType_Typo, "phosphotransferase", SimpleReplaceFunc } ,
  { "phosphotase" , NormalSearch, eSuspectNameType_Typo, "phosphatase", SimpleReplaceFunc } ,
  { "posible" , NormalSearch, eSuspectNameType_Typo, "possible", SimpleReplaceFunc } ,
  { "presursor" , NormalSearch, eSuspectNameType_Typo, "precursor", SimpleReplaceFunc } ,
  { "probable putative" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "proein" , NormalSearch, eSuspectNameType_Typo, "protein", SimpleReplaceFunc } ,
  { "prortein" , NormalSearch, eSuspectNameType_Typo, "protein", SimpleReplaceFunc } ,
  { "proteine" , NormalSearch, eSuspectNameType_Typo, "protein", SimpleReplaceFunc } ,
  { "proteinn" , NormalSearch, eSuspectNameType_Typo, "protein", SimpleReplaceFunc } ,
  { "protien" , NormalSearch, eSuspectNameType_Typo, "protein", SimpleReplaceFunc } ,
  { "protrein" , NormalSearch, eSuspectNameType_Typo, "protein", SimpleReplaceFunc } ,
  { "prptein" , NormalSearch, eSuspectNameType_Typo, "protein", SimpleReplaceFunc } ,
  { "pseudogene" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "puatative" , NormalSearch, eSuspectNameType_Typo, "putative", SimpleReplaceFunc } ,
  { "puative" , NormalSearch, eSuspectNameType_Typo, "putative", SimpleReplaceFunc } ,
  { "putaitive" , NormalSearch, eSuspectNameType_Typo, "putative", SimpleReplaceFunc } ,
  { "putaitve" , NormalSearch, eSuspectNameType_Typo, "putative", SimpleReplaceFunc } ,
  { "putaive" , NormalSearch, eSuspectNameType_Typo, "putative", SimpleReplaceFunc } ,
  { "putataive" , NormalSearch, eSuspectNameType_Typo, "putative", SimpleReplaceFunc } ,
  { "putatitve" , NormalSearch, eSuspectNameType_Typo, "putative", SimpleReplaceFunc } ,
  { "putative orphan protein" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "putative probable" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "putative putative" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "putative, putative" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "putatuve" , NormalSearch, eSuspectNameType_Typo, "putative", SimpleReplaceFunc } ,
  { "putatve" , NormalSearch, eSuspectNameType_Typo, "putative", SimpleReplaceFunc } ,
  { "putatvie" , NormalSearch, eSuspectNameType_Typo, "putative", SimpleReplaceFunc } ,
  { "putayive" , NormalSearch, eSuspectNameType_Typo, "putative", SimpleReplaceFunc } ,
  { "putitive" , NormalSearch, eSuspectNameType_Typo, "putative", SimpleReplaceFunc } ,
  { "qlcohol" , NormalSearch, eSuspectNameType_Typo, "alcohol", SimpleReplaceFunc } ,
  { "recognised" , NormalSearch, eSuspectNameType_Typo, "recognized", SimpleReplaceFunc } ,
  { "regulatot" , NormalSearch, eSuspectNameType_Typo, "regulator", SimpleReplaceFunc } ,
  { "reponse" , NormalSearch, eSuspectNameType_Typo, "response", SimpleReplaceFunc } ,
  { "resistence" , NormalSearch, eSuspectNameType_Typo, "resistance", SimpleReplaceFunc } ,
  { "ribosimal" , NormalSearch, eSuspectNameType_Typo, "ribosomal", SimpleReplaceFunc } ,
  { "ribosoml" , NormalSearch, eSuspectNameType_Typo, "ribosomal", SimpleReplaceFunc } ,
  { "sapiens" , NormalSearch, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "serinr" , NormalSearch, eSuspectNameType_Typo, "serine", SimpleReplaceFunc } ,
  { "signalling" , NormalSearch, eSuspectNameType_Typo, "signaling", SimpleReplaceFunc } ,
  { "similar" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "simmilar" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "specfic" , NormalSearch, eSuspectNameType_Typo, "specific", SimpleReplaceFunc } ,
  { "sphaeroides" , NormalSearch, eSuspectNameType_RemoveOrganismName, NULL, NULL } ,
  { "spscific" , NormalSearch, eSuspectNameType_Typo, "specific", SimpleReplaceFunc } ,
  { "stabilisation" , NormalSearch, eSuspectNameType_Typo, "stabilization", SimpleReplaceFunc } ,
  { "subnit" , NormalSearch, eSuspectNameType_Typo, "subunit", SimpleReplaceFunc } ,
  { "suger" , NormalSearch, eSuspectNameType_Typo, "sugar", SimpleReplaceFunc } ,
  { "sulpho" , NormalSearch, eSuspectNameType_None, "sulfo", SimpleReplaceFunc } ,
  { "sulphur" , NormalSearch, eSuspectNameType_Typo, "sulfur", SimpleReplaceFunc } ,
  { "systhesis" , NormalSearch, eSuspectNameType_Typo, "synthesis", SimpleReplaceFunc } ,
  { "sythase" , NormalSearch, eSuspectNameType_Typo, "synthase", SimpleReplaceFunc } ,
  { "thiredoxin" , NormalSearch, eSuspectNameType_Typo, "thioredoxin", SimpleReplaceFunc } ,
  { "trancsriptional" , NormalSearch, eSuspectNameType_Typo, "transcription", SimpleReplaceFunc } ,
  { "tranferase" , NormalSearch, eSuspectNameType_Typo, "transferase", SimpleReplaceFunc } ,
  { "tranporter" , NormalSearch, eSuspectNameType_Typo, "transporter", SimpleReplaceFunc } ,
  { "transcirbed" , NormalSearch, eSuspectNameType_Typo, "transcribed", SimpleReplaceFunc } ,
  { "transcriptonal" , NormalSearch, eSuspectNameType_Typo, "transcriptional", SimpleReplaceFunc } ,
  { "transcritional" , NormalSearch, eSuspectNameType_Typo, "transcriptional", SimpleReplaceFunc } ,
  { "transebrane" , NormalSearch, eSuspectNameType_Typo, "transmembrane", SimpleReplaceFunc } ,
  { "transglycolase" , NormalSearch, eSuspectNameType_Typo, "transglycosylase", SimpleReplaceFunc } ,
  { "transorter" , NormalSearch, eSuspectNameType_Typo, "transporter", SimpleReplaceFunc } ,
  { "transpoase" , NormalSearch, eSuspectNameType_Typo, "transposase", SimpleReplaceFunc } ,
  { "transportor" , NormalSearch, eSuspectNameType_Typo, "transporter", SimpleReplaceFunc } ,
  { "transproter" , NormalSearch, eSuspectNameType_Typo, "transporter", SimpleReplaceFunc } ,
  { "transulfuration" , NormalSearch, eSuspectNameType_Typo, "transsulfuration", SimpleReplaceFunc } ,
  { "trnasporter" , NormalSearch, eSuspectNameType_Typo, "transporter", SimpleReplaceFunc } ,
  { "truncat" , NormalSearch, eSuspectNameType_MightBeNonfunctional, NULL, NULL } ,
  { "ttg start" , NormalSearch, eSuspectNameType_MightBeNonfunctional, NULL, NULL } ,
  { "tumour" , NormalSearch, eSuspectNameType_Typo, "tumor", SimpleReplaceFunc } ,
  { "typr" , NormalSearch, eSuspectNameType_Typo, "type", SimpleReplaceFunc } ,
  { "uncharacterized protein", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } ,
  { "uncharaterized" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "undecapaprenyl" , NormalSearch, eSuspectNameType_Typo, "undecaprenyl", SimpleReplaceFunc } ,
  { "unkown" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "utilising" , NormalSearch, eSuspectNameType_Typo, "utilizing", SimpleReplaceFunc } ,
  { "weakly conserved" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "widely conserved" , NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "|" , NormalSearch, eSuspectNameType_InappropriateSymbol, NULL, NULL } ,
  { "C term" , ProductContainsTerm, eSuspectNameType_MightBeNonfunctional, NULL, NULL } ,
  { "C-term" , ProductContainsTerm, eSuspectNameType_MightBeNonfunctional, NULL, NULL } ,
  { "N term" , ProductContainsTerm, eSuspectNameType_MightBeNonfunctional, NULL, NULL } ,
  { "N-term" , ProductContainsTerm, eSuspectNameType_MightBeNonfunctional, NULL, NULL } ,
  { "Two or more sets of brackets or parentheseis" , ContainsTwoSetsOfBracketsOrParentheses, eSuspectNameType_None, NULL, NULL } ,
  { "unknown" , ContainsUnknownName, eSuspectNameType_None, NULL, NULL } ,
  { "double space" , ContainsDoubleSpace, eSuspectNameType_None, NULL, NULL } ,
  { "COG" , ContainsWholeWordCaseSensitive, eSuspectNameType_Database, NULL, NULL } ,
  { "DUF" , ContainsWholeWordCaseSensitive, eSuspectNameType_Database, NULL, NULL } ,
  { "EST" , ContainsWholeWordCaseSensitive, eSuspectNameType_Database, NULL, NULL } ,
  { "FOG" , ContainsWholeWordCaseSensitive, eSuspectNameType_Database, NULL, NULL } ,
  { "UPF" , ContainsWholeWordCaseSensitive, eSuspectNameType_Database, NULL, NULL } ,
  { "_" , ContainsUnderscore, eSuspectNameType_Database, NULL, NULL } ,
  { "ending with period, comma, hyphen, underscore, colon, or forward slash" , EndsWithPunct, eSuspectNameType_InappropriateSymbol, NULL, NULL } ,
  { "PTS system" , IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_UseProtein, NULL, NULL } ,
  { "helix-turn-helix" , IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_UseProtein, NULL, NULL } ,
  { "transposase of" , IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_None, NULL, NULL } ,
  { "zinc finger" , IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_UseProtein, NULL, NULL } ,
  { "may contain a plural" , MayContainPlural, eSuspectNameType_None, NULL, NULL } ,
  { "unbalanced brackets or parentheses" , ContainsUnbalancedParentheses, eSuspectNameType_InappropriateSymbol, NULL, NULL } ,
  { "long product name that may contain descriptive information more appropriate in a note", IsTooLong, eSuspectNameType_QuickFix, NULL, NULL } ,
  { "Product name begins with possible, potential, predicted or probable.  Please use putative.", StartsWithPutativeReplacement, eSuspectNameType_QuickFix, "putative", UsePutative } ,

  { "CDS", NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "doubtful", NormalSearch, eSuspectNameType_None, NULL, NULL } ,
  { "alternate protein name", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "conser", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "conserve", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "conserved hypothetical", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "conserved hypothetical protein", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "conserved", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "domain family", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "domain of unknown function", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "domain protein", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "domain", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "doubtful CDS found within S. typhi pathogenicity island", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "factor", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "family protein", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "hypo", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "hypothetical ORF", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "hypothetical", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "hypothetical domain protein", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "No definition line found", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "orphan protein", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "ORF", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "orf, hyp", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "orf, hypothetical", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "peptide", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "precursor", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "probable", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "predicted", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "predicted protein", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "probable protein", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "protein containing", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "protein of unknown function", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "protein-containing", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "pseudo", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "putative conserved hypothetical", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "putative hypothetical", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "putative protein", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "putative", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "uncharacterized conserved protein", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "unnamed", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameFunc } , 
  { "o252", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "o252 protein", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Alanine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Arginine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Asparagine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Aspartic acid", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Cysteine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "DNA", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Glutamic acid", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Glutamine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Glycine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Histidine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Isoleucine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Leucine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Lysine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Methionine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "NAD", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "PASTA", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_UseProtein, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Phenylalanine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Proline", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "RNA", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Serine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Threonine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Tryptophan", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Tyrosine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "Valine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "adenine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "amino acid", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "barrel", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "carbon", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "citrate", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "cytosine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "finger", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "ggdef", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "guanine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "helium", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "helix", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "hydrogen", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "insertion sequence", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "iron", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "mRNA", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "membrane", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "ncRNA", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "nitrogen", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "oxygen", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "p-loop", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_UseProtein, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "peptide", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "phage", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "plasmid", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "purine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "rRNA", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "repeat", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "secreted", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "signal peptide", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_UseProtein, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "signal", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "subunit", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "tRNA", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "thymine", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "transport-associated", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "transposon", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "uracil", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } , 
  { "zinc", IsSingleWordOrWeaselPlusSingleWord, eSuspectNameType_QuickFix, "hypothetical protein", ReplaceWholeNameAddNoteFunc } 
};


const int num_suspect_product_terms = sizeof (suspect_product_terms) / sizeof (SuspectProductNameData);


static void FixSuspectProductNameTyposInOneFeature (SeqFeatPtr cds, LogInfoPtr lip, ESuspectNameType fix_type)
{
  Int4            k;
  ProtRefPtr      prp;
  ValNodePtr      vnp;
  CharPtr         tmp, desc;
  ValNode         vn;
  SeqFeatPtr      mrna;
  SeqMgrFeatContext context;
  RnaRefPtr         rrp;
  CharPtr           extra;
  CharPtr           and_associated_mrna = " and associated mRNA";

  if (cds == NULL || cds->data.choice != SEQFEAT_CDREGION || cds->data.value.ptrvalue == NULL
      || cds->product == NULL || (prp = GetProtRefForFeature(cds)) == NULL)
  {
    return;
  }


  
  for (k = 0; k < num_suspect_product_terms; k++)
  {    
    for (vnp = prp->name; vnp != NULL; vnp = vnp->next) 
    {
      if (suspect_product_terms[k].fix_type == fix_type
        && suspect_product_terms[k].replace_func != NULL
        && suspect_product_terms[k].search_func != NULL
        && (suspect_product_terms[k].search_func) (suspect_product_terms[k].pattern, vnp->data.ptrvalue)) 
      {
        if (lip != NULL && lip->fp != NULL) {
          tmp = StringSave ((CharPtr) vnp->data.ptrvalue);
          (suspect_product_terms[k].replace_func)(&tmp, 
			                                      suspect_product_terms[k].pattern, 
												  suspect_product_terms[k].replace_phrase,
												  cds);
          if (StringCmp (tmp, vnp->data.ptrvalue) != 0) {
            extra = "";
            mrna = SeqMgrGetOverlappingmRNA (cds->location, &context);
            if (mrna != NULL && mrna->data.choice == SEQFEAT_RNA
                && (rrp = mrna->data.value.ptrvalue) != NULL
                && rrp->ext.choice == 1
                && StringCmp (rrp->ext.value.ptrvalue, vnp->data.ptrvalue) == 0) {
                rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
                rrp->ext.value.ptrvalue = StringSave (tmp);
                extra = and_associated_mrna;
            }
            MemSet (&vn, 0, sizeof (ValNode));
            vn.choice = OBJ_SEQFEAT;
            vn.data.ptrvalue = cds;
            desc = GetDiscrepancyItemText (&vn);
            fprintf (lip->fp, "Changed '%s' to '%s' for %s%s\n", (CharPtr) vnp->data.ptrvalue, tmp, desc, extra);
            vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
            vnp->data.ptrvalue = tmp;
            tmp = NULL;
            desc = MemFree (desc);
            lip->data_in_log = TRUE;
          }
          tmp = MemFree (tmp);
        } else {
          tmp = (CharPtr) vnp->data.ptrvalue;
          (suspect_product_terms[k].replace_func)(&tmp, suspect_product_terms[k].pattern, suspect_product_terms[k].replace_phrase, cds);
          vnp->data.ptrvalue = tmp;
        }
        break;
      }
      /* only check the first name */
      if (!StringHasNoText (vnp->data.ptrvalue)) {
        break;
      }
    }
  }
}


static void FixSuspectProductNameTypos (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr vnp;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT) {
      FixSuspectProductNameTyposInOneFeature ((SeqFeatPtr) vnp->data.ptrvalue, lip, eSuspectNameType_Typo);
    }
  }
}


static void FixSuspectProductNameQuickFixes (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr vnp;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT) {
      FixSuspectProductNameTyposInOneFeature ((SeqFeatPtr) vnp->data.ptrvalue, lip, eSuspectNameType_QuickFix);
    }
  }
}

static void FindSuspectProductNamesCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR feature_list;
  Int4            k;
  ProtRefPtr      prp;
  ValNodePtr      vnp;
  BioseqPtr       bsp;
  SeqFeatPtr      cds;
  BioSourcePtr    biop = NULL;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT || sfp->data.value.ptrvalue == NULL
      || userdata == NULL)
  {
    return;
  }
  
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  feature_list = (ValNodePtr PNTR) userdata;

  /* add coding region rather than protein */
  if (sfp->idx.subtype == FEATDEF_PROT) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
      if (cds != NULL) {
        sfp = cds;
      }
      /* find BioSource, to check whether we want to run all categories */
      biop = GetBiopForBsp (bsp);
    }
  }
  
  for (k = 0; k < num_suspect_product_terms; k++)
  {
    if (!CategoryOkForBioSource(biop, suspect_product_terms[k].fix_type)) {
      continue;
    }
    for (vnp = prp->name; vnp != NULL; vnp = vnp->next) 
    {
      if (suspect_product_terms[k].search_func != NULL
        && (suspect_product_terms[k].search_func) (suspect_product_terms[k].pattern, vnp->data.ptrvalue)) 
      {
        ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
        break;
      }
      /* only check the first name */
      if (!StringHasNoText (vnp->data.ptrvalue)) {
        break;
      }
    }
  }
        
}


static ClickableItemPtr SuspectPhraseEx (Uint4 clickable_item_type, CharPtr phrase, Boolean quote_phrase, CharPtr feat_type, ValNodePtr feature_list)
{
  ClickableItemPtr dip = NULL;
  CharPtr          bad_fmt_quote = "%d %ss contain '%s'";
  CharPtr          bad_fmt_noquote = "%d %ss contain %s";
  CharPtr          bad_fmt;

  if (feature_list == NULL || phrase == NULL || StringHasNoText (feat_type))
  {
    return NULL;
  }
  if (quote_phrase) {
    bad_fmt = bad_fmt_quote;
  } else {
    bad_fmt = bad_fmt_noquote;
  }
  
  dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  if (dip != NULL)
  {
    dip->clickable_item_type = clickable_item_type;
    dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + StringLen (phrase) + StringLen (feat_type) + 15));
    sprintf (dip->description, bad_fmt, ValNodeLen (feature_list), feat_type, phrase);
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    dip->item_list = feature_list;
  }      
  return dip;
}


static ClickableItemPtr SuspectPhrase (Uint4 clickable_item_type, CharPtr phrase, CharPtr feat_type, ValNodePtr feature_list)
{
  return SuspectPhraseEx (clickable_item_type, phrase, TRUE, feat_type, feature_list);
}


static ClickableItemPtr SuspectPhraseEnd (Uint4 clickable_item_type, CharPtr phrase, CharPtr feat_type, ValNodePtr feature_list)
{
  ClickableItemPtr dip = NULL;
  CharPtr            bad_fmt = "%d %ss end with %s";

  if (feature_list == NULL || StringHasNoText (phrase) || StringHasNoText (feat_type))
  {
    return NULL;
  }
  
  dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  if (dip != NULL)
  {
    dip->clickable_item_type = clickable_item_type;
    dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + StringLen (phrase) + StringLen (feat_type) + 15));
    sprintf (dip->description, bad_fmt, ValNodeLen (feature_list), feat_type, phrase);
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    dip->item_list = feature_list;
  }      
  return dip;
}


static ClickableItemPtr SuspectPhraseStart (Uint4 clickable_item_type, CharPtr phrase, CharPtr feat_type, ValNodePtr feature_list)
{
  ClickableItemPtr dip = NULL;
  CharPtr            bad_fmt = "%d %ss start with %s";

  if (feature_list == NULL || StringHasNoText (phrase) || StringHasNoText (feat_type))
  {
    return NULL;
  }
  
  dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  if (dip != NULL)
  {
    dip->clickable_item_type = clickable_item_type;
    dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + StringLen (phrase) + StringLen (feat_type) + 15));
    sprintf (dip->description, bad_fmt, ValNodeLen (feature_list), feat_type, phrase);
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    dip->item_list = feature_list;
  }      
  return dip;
}


static Uint4 ClickableItemTypeForNameCat (Int4 k)
{
  if (k == eSuspectNameType_Typo) {
    return DISC_PRODUCT_NAME_TYPO;
  } else if (k == eSuspectNameType_QuickFix) {
    return DISC_PRODUCT_NAME_QUICKFIX;
  } else {
    return DISC_SUSPECT_PRODUCT_NAME;
  }
}

typedef struct suspectrulefeats {
  SuspectRuleSetPtr rule_list;
  ValNodePtr PNTR   feature_list;
  Int4 num_rules;
} SuspectRuleFeatsData, PNTR SuspectRuleFeatsPtr;


static void FindSuspectProductNamesWithRulesCallback (SeqFeatPtr sfp, Pointer userdata)
{
  SuspectRuleFeatsPtr srlist;
  SuspectRulePtr  rule;
  Int4            k;
  ProtRefPtr      prp;
  BioseqPtr       bsp;
  SeqFeatPtr      cds;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT || sfp->data.value.ptrvalue == NULL
      || (srlist = (SuspectRuleFeatsPtr)userdata) == NULL)
  {
    return;
  }
  
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;

  if (prp == NULL || prp->name == NULL) {
    return;
  }

  /* add coding region rather than protein */
  if (sfp->idx.subtype == FEATDEF_PROT) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
      if (cds != NULL) {
        sfp = cds;
      }
    }
  }
  
  for (k = 0, rule = srlist->rule_list; k < srlist->num_rules && rule != NULL; k++, rule = rule->next)
  {
    if (DoesStringMatchSuspectRule (prp->name->data.ptrvalue, sfp, rule)) 
    {
      ValNodeAddPointer (&(srlist->feature_list[k]), OBJ_SEQFEAT, sfp);
    }
  }
        
}


static void AutoFixSuspectProductRules (ValNodePtr item_list, Pointer userdata, LogInfoPtr lip)
{
  SuspectRulePtr rule;
  ValNodePtr     vnp;

  if ((rule = (SuspectRulePtr) userdata) == NULL || item_list == NULL) {
    return;
  }

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT) {
      if (ApplySuspectProductNameFixToFeature (rule, (SeqFeatPtr) vnp->data.ptrvalue, lip == NULL ? NULL : lip->fp)) {
        if (lip != NULL) {
          lip->data_in_log = TRUE;
        }
      }
    }
  }
}


static void 
FindSuspectProductNamesWithRules 
(ValNodePtr PNTR discrepancy_list,
 ValNodePtr sep_list,
 SuspectRuleSetPtr rule_list)
{
  SuspectRuleFeatsData srdata;
  SuspectRulePtr       rule;
  CharPtr              summ;
  CharPtr              fmt = "%d features %s";
  ValNodePtr PNTR      name_cat;
  ValNodePtr         master_list = NULL, vnp;
  Int4               k;
  ClickableItemPtr   dip, tdip = NULL;
  ValNodePtr         subcategories = NULL;
  Int4               num_cat = Fix_type_gene + 1;

  if (discrepancy_list == NULL) return;

  srdata.num_rules = CountSuspectRuleSet (rule_list);
  if (srdata.num_rules == 0) {
    return;
  }

  srdata.rule_list = rule_list;
  srdata.feature_list = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * srdata.num_rules);
  if (srdata.feature_list == NULL) return;

  name_cat = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * num_cat);

  /* initialize array for suspicious product names */
  for (k = 0; k < srdata.num_rules; k++)
  {
    srdata.feature_list[k] = NULL;
  }

  /* initialize named categories */
  for (k = 0; k < num_cat; k++) {
    name_cat[k] = NULL;
  }
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) 
  {
    VisitGenProdSetFeatures (vnp->data.ptrvalue, &srdata, FindSuspectProductNamesWithRulesCallback);
  }

  for (k = 0, rule = srdata.rule_list; k < srdata.num_rules && rule != NULL; k++, rule = rule->next)
  {
    if (srdata.feature_list[k] != NULL)
    {
      summ = SummarizeSuspectRuleEx(rule, TRUE);
      dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      dip->clickable_item_type = DISC_SUSPECT_PRODUCT_NAME;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (summ) + 15));
      sprintf (dip->description, fmt, ValNodeLen (srdata.feature_list[k]), summ);
      summ = MemFree (summ);
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = srdata.feature_list[k];
      if (rule->replace != NULL) {
        dip->autofix_func = AutoFixSuspectProductRules;
        dip->autofix_data = rule;
      }
      ValNodeAddPointer (&name_cat[rule->rule_type], 0, dip);
      ValNodeLinkCopy (&master_list, srdata.feature_list[k]);
    }
  }
  if (master_list != NULL)
  {
    for (k = 0; k < num_cat; k++) {
      if (name_cat[k] != NULL) {
        tdip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        MemSet (tdip, 0, sizeof (ClickableItemData));
        tdip->description = StringSave (SummarizeFixType(k));
        tdip->item_list = ItemListFromSubcategories (name_cat[k]);
        tdip->clickable_item_type = DISC_SUSPECT_PRODUCT_NAME;
        tdip->subcategories = name_cat[k];
        tdip->expanded = TRUE;
        ValNodeAddPointer (&subcategories, 0, tdip);
      }
    }
    dip = SuspectPhraseEx (DISC_SUSPECT_PRODUCT_NAME, "suspect phrase or characters", FALSE, "product_name", master_list);
    if (dip != NULL)
    {
      dip->subcategories = subcategories;      
      dip->expanded = TRUE;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }

  MemFree (srdata.feature_list);
  MemFree (name_cat);
}



static void FindSuspectProductNamesWithStaticList (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr PNTR   feature_list = NULL;
  ValNodePtr         master_list = NULL, vnp;
  Int4               k;
  ClickableItemPtr   dip, tdip = NULL;
  ValNodePtr         name_cat[eSuspectNameType_Max];
  ValNodePtr         subcategories = NULL;
  
  if (discrepancy_list == NULL) return;

  feature_list = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * num_suspect_product_terms);
  if (feature_list == NULL) return;

  MemSet (&name_cat, 0, sizeof (name_cat));
  
  /* initialize array for suspicious product names */
  for (k = 0; k < num_suspect_product_terms; k++)
  {
    feature_list[k] = NULL;
  }
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) 
  {
    VisitGenProdSetFeatures (vnp->data.ptrvalue, feature_list, FindSuspectProductNamesCallback);
  }
  
  for (k = 0; k < num_suspect_product_terms; k++)
  {
    if (feature_list[k] != NULL)
    {
      if (suspect_product_terms[k].search_func == EndsWithPattern) 
      {
        dip = SuspectPhraseEnd (ClickableItemTypeForNameCat(suspect_product_terms[k].fix_type), suspect_product_terms[k].pattern, "product name", feature_list[k]);
      }
      else if (suspect_product_terms[k].search_func == StartsWithPattern) 
      {
        dip = SuspectPhraseStart (ClickableItemTypeForNameCat(suspect_product_terms[k].fix_type), suspect_product_terms[k].pattern, "product name", feature_list[k]);
      }
      else 
      {
        dip = SuspectPhrase (ClickableItemTypeForNameCat(suspect_product_terms[k].fix_type), suspect_product_terms[k].pattern, "product name", feature_list[k]);
      }
      if (dip != NULL)
      {
        ValNodeAddPointer (&name_cat[suspect_product_terms[k].fix_type], 0, dip);
      }
      ValNodeLinkCopy (&master_list, feature_list[k]);
    }
  }
  if (master_list != NULL)
  {
    for (k = 0; k < eSuspectNameType_Max; k++) {
      if (name_cat[k] != NULL) {
        tdip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        MemSet (tdip, 0, sizeof (ClickableItemData));
        tdip->description = StringSave (suspect_name_category_names[k]);
        tdip->item_list = ItemListFromSubcategories (name_cat[k]);
        tdip->clickable_item_type = ClickableItemTypeForNameCat(suspect_product_terms[k].fix_type);
        tdip->subcategories = name_cat[k];
        tdip->expanded = TRUE;
        ValNodeAddPointer (&subcategories, 0, tdip);
      }
    }
    dip = SuspectPhraseEx (DISC_SUSPECT_PRODUCT_NAME, "suspect phrase or characters", FALSE, "product_name", master_list);
    if (dip != NULL)
    {
      dip->subcategories = subcategories;      
      dip->expanded = TRUE;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }

  MemFree (feature_list);
}


static SuspectRuleSetPtr s_SuspectProductRuleList = NULL;
static Boolean s_TriedToReadRules = FALSE;

extern void FindSuspectProductNames (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  Char rule_file[PATH_MAX];
  AsnIoPtr aip;
  

  if (s_SuspectProductRuleList == NULL && !s_TriedToReadRules) 
  {
    if (GetAppParam ("SEQUINCUSTOM", "SETTINGS", "PRODUCT_RULES_LIST", NULL, rule_file, sizeof (rule_file) - 1)
         || GetAppParam ("SEQUIN", "SETTINGS", "PRODUCT_RULES_LIST", NULL, rule_file, sizeof (rule_file) - 1))
    {    
      if ((aip = AsnIoOpen (rule_file, "r")) == NULL) {
        Message (MSG_ERROR, "Unable to read %s", rule_file);
      } else {
        if ((s_SuspectProductRuleList = SuspectRuleSetAsnRead (aip, NULL)) == NULL) {
          Message (MSG_ERROR, "Unable to read suspect product rules from %s", rule_file);
        }
        AsnIoClose (aip);
      }
    }
    if (s_SuspectProductRuleList == NULL) 
    {
      if (FindPath ("ncbi", "ncbi", "data", rule_file, sizeof (rule_file))) 
      {
        FileBuildPath (rule_file, NULL, "product_rules.prt");
        if ((aip = AsnIoOpen (rule_file, "r")) == NULL) {
          Message (MSG_ERROR, "Unable to read %s", rule_file);
        } else {
          if ((s_SuspectProductRuleList = SuspectRuleSetAsnRead (aip, NULL)) == NULL) {
            Message (MSG_ERROR, "Unable to read suspect product rules from %s", rule_file);
          }
          AsnIoClose (aip);
        }
      }
    }
    s_TriedToReadRules = TRUE;
  }
  if (s_SuspectProductRuleList == NULL) 
  {
    FindSuspectProductNamesWithStaticList(discrepancy_list, sep_list);
  } 
  else 
  {
    FindSuspectProductNamesWithRules(discrepancy_list, sep_list, s_SuspectProductRuleList);
  }
}


NLM_EXTERN Boolean IsProductNameOk (CharPtr product_name)
{
  Int4     k;
  Boolean  rval = TRUE;

  for (k = 0; k < num_suspect_product_terms && rval; k++)
  {    
    if (suspect_product_terms[k].search_func != NULL
      && (suspect_product_terms[k].search_func) (suspect_product_terms[k].pattern, product_name)) 
    {
      rval = FALSE;
    }
  }
  return rval;
}


NLM_EXTERN Boolean ReportProductNameProblems (CharPtr product_name, FILE *output_file, CharPtr prefix)
{
  Int4 k;
  Boolean any_problems = FALSE;
  CharPtr func_name;

  for (k = 0; k < num_suspect_product_terms; k++)
  {    
    if (suspect_product_terms[k].search_func != NULL
      && (suspect_product_terms[k].search_func) (suspect_product_terms[k].pattern, product_name)) 
    {
      if (suspect_product_terms[k].search_func == EndsWithPattern) {
        func_name = "Ends with";
      } else if (suspect_product_terms[k].search_func == StartsWithPattern) {
        func_name = "Starts with";
      } else {
        func_name = "Contains";
      }
      if (output_file) {
        if (prefix == NULL) {
          fprintf (output_file, "%s\t%s '%s'\n", product_name, func_name, suspect_product_terms[k].pattern);
        } else {
          fprintf (output_file, "%s\t%s\t%s '%s'\n", prefix, product_name, func_name, suspect_product_terms[k].pattern);
        }
      } else {
        if (prefix == NULL) {
          printf ("%s\t%s '%s'\n", product_name, func_name, suspect_product_terms[k].pattern);
        } else {
          printf ("%s\t%s\t%s '%s'\n", prefix, product_name, func_name, suspect_product_terms[k].pattern);
        }
      }
      any_problems = TRUE;
    }
  }
  return any_problems;
}


static CharPtr suspect_phrases[] = 
{
"fragment",
"frameshift",
"%",
"E-value",
"E value",
"Evalue"
};

const int num_suspect_phrases = sizeof (suspect_phrases) / sizeof (CharPtr);


static void FindSuspectPhrasesCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR feature_list;
  Int4            k;
  ProtRefPtr      prp;
  CharPtr         check_str = NULL;
  
  if (sfp == NULL || (sfp->data.choice != SEQFEAT_PROT && sfp->data.choice != SEQFEAT_CDREGION) || sfp->data.value.ptrvalue == NULL
      || userdata == NULL)
  {
    return;
  }
  
  if (sfp->data.choice == SEQFEAT_PROT) {  
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    check_str = prp->desc;
  } else if (sfp->data.choice == SEQFEAT_CDREGION) {
    check_str = sfp->comment;
  }
  if (StringHasNoText (check_str)) return;
  
  feature_list = (ValNodePtr PNTR) userdata;
  
  for (k = 0; k < num_suspect_phrases; k++)
  {
    if (StringISearch(check_str, suspect_phrases[k]) != NULL)
    {
      ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
      break;
    }
  }  
}

extern void FindSuspectPhrases (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr PNTR  feature_list = NULL;
  ValNodePtr       vnp, subcat = NULL;
  ClickableItemPtr dip;
  Int4             k;
  
  if (discrepancy_list == NULL) return;

  feature_list = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * num_suspect_phrases);
  for (k = 0; k < num_suspect_phrases; k++)
  {
    feature_list[k] = NULL;
  }

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, feature_list, FindSuspectPhrasesCallback);
  }

  for (k = 0; k < num_suspect_phrases; k++)
  {
    if (feature_list[k] != NULL) {
      dip = SuspectPhrase (DISC_SUSPECT_PHRASES, suspect_phrases[k], "cds comments or protein description", feature_list[k]);
      if (dip != NULL)
      {
        ValNodeAddPointer (&subcat, 0, dip);
      }
    }
  }

  if (subcat != NULL) 
  {
    dip = SuspectPhraseEx (DISC_SUSPECT_PRODUCT_NAME, "suspect phrases", FALSE, "cds comments or protein description", ItemListFromSubcategories (subcat));
    if (dip != NULL)
    {
      dip->subcategories = subcat;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }

  MemFree (feature_list);  
}


static void CheckForStr 
(CharPtr    check_str, 
 SeqFeatPtr sfp,
 ValNodePtr PNTR feature_list,
 CharPtr PNTR find_list,
 Int4         num_find)
{
  Int4 k;

  if (StringHasNoText (check_str) || sfp == NULL || feature_list == NULL || find_list == NULL) {
    return;
  }

  for (k = 0; k < num_find; k++)
  {
    if (StringISearch(check_str, find_list[k]) != NULL)
    {
      ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
      break;
    }
  }  

}


static CharPtr suspicious_note_phrases[] = 
{
 "characterised",
 "recognised",
 "characterisation",
 "localisation",
 "tumour",
 "uncharacterised",
 "oxydase",
 "colour",
 "localise",
 "faecal",
 "orthologue",
 "paralogue",
 "homolog",
 "homologue",
 "intronless gene"
};

const int num_suspicious_note_phrases = sizeof (suspicious_note_phrases) / sizeof (CharPtr);


static void FindSuspiciousNoteTextCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR feature_list;
  ProtRefPtr      prp;
  GeneRefPtr      grp;
  
  if (sfp == NULL || (feature_list = (ValNodePtr PNTR)userdata) == NULL) {
    return;
  }

  if (sfp->data.choice == SEQFEAT_GENE) {
    /* look in gene comment and gene description */
    CheckForStr (sfp->comment, sfp, feature_list, suspicious_note_phrases, num_suspicious_note_phrases);
    if ((grp = sfp->data.value.ptrvalue) != NULL) {
      CheckForStr (grp->desc, sfp, feature_list, suspicious_note_phrases, num_suspicious_note_phrases);
    }
  } else if (sfp->data.choice == SEQFEAT_CDREGION) {
    /* look in CDS comment */
    CheckForStr (sfp->comment, sfp, feature_list, suspicious_note_phrases, num_suspicious_note_phrases);
  } else if (sfp->idx.subtype == FEATDEF_PROT) {
    /* look in protein description */
    if ((prp = sfp->data.value.ptrvalue) != NULL) {
      CheckForStr (prp->desc, sfp, feature_list, suspicious_note_phrases, num_suspicious_note_phrases);
    }
  } else if (sfp->idx.subtype == FEATDEF_misc_feature) {
    /* look in misc_feature comment */
    CheckForStr (sfp->comment, sfp, feature_list, suspicious_note_phrases, num_suspicious_note_phrases);
  }
}


static void FindSuspiciousPhraseInNoteText (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr PNTR  feature_list = NULL;
  ValNodePtr       vnp, subcat = NULL;
  ClickableItemPtr dip;
  Int4             k;
  
  if (discrepancy_list == NULL) return;

  feature_list = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * num_suspicious_note_phrases);
  for (k = 0; k < num_suspicious_note_phrases; k++)
  {
    feature_list[k] = NULL;
  }

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, feature_list, FindSuspiciousNoteTextCallback);
  }

  for (k = 0; k < num_suspicious_note_phrases; k++)
  {
    if (feature_list[k] != NULL) {
      dip = SuspectPhrase (DISC_SUSPICIOUS_NOTE_TEXT, suspicious_note_phrases[k], "note text", feature_list[k]);
      if (dip != NULL)
      {
        ValNodeAddPointer (&subcat, 0, dip);
      }
    }
  }

  if (subcat != NULL) 
  {
    dip = SuspectPhraseEx (DISC_SUSPICIOUS_NOTE_TEXT, "suspicious phrases", FALSE, "note text", ItemListFromSubcategories (subcat));
    if (dip != NULL)
    {
      dip->subcategories = subcat;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }

  MemFree (feature_list);  
}


static void FindUnknownProteinsWithECNumbersCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ProtRefPtr prp;
  ValNodePtr PNTR feature_list;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT || sfp->data.value.ptrvalue == NULL || userdata == NULL) 
  {
    return;
  }
  
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp->name == NULL || prp->ec == NULL) return;

  if (StringISearch (prp->name->data.ptrvalue, "hypothetical protein") != NULL
      || StringISearch (prp->name->data.ptrvalue, "unknown protein") != NULL) 
  {
    feature_list = (ValNodePtr PNTR) userdata;
    ValNodeAddPointer (feature_list, OBJ_SEQFEAT, sfp);
  }
}


extern void FindUnknownProteinsWithECNumbers (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr         feature_list = NULL, vnp;
  ClickableItemPtr dip;
  
  if (discrepancy_list == NULL) return;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &feature_list, FindUnknownProteinsWithECNumbersCallback);
  }
  
  if (feature_list != NULL) {
    dip = NewClickableItem (DISC_EC_NUMBER_ON_HYPOTHETICAL_PROTEIN, "%d protein features have an EC number and a protein name of 'unknown protein' or 'hypothetical protein'", feature_list);
    ValNodeAddPointer (discrepancy_list, 0, dip);
  }  
}

static ClickableItemPtr InconsistentSourceDefline (SeqDescrPtr biop_sdp, SeqDescrPtr title_sdp)
{
  ClickableItemPtr dip = NULL;
  CharPtr            bad_fmt = "Organism description not found in definition line: %s.";
  BioSourcePtr       biop;
  CharPtr            desc = NULL;

  if (biop_sdp == NULL || title_sdp == NULL)
  {
    return NULL;
  }
  
  biop = (BioSourcePtr) biop_sdp->data.ptrvalue;
  if (biop != NULL && biop->org != NULL && !StringHasNoText (biop->org->taxname))
  {
    desc = biop->org->taxname;
  }
  else
  {
    desc = title_sdp->data.ptrvalue;
  }
  if (StringHasNoText (desc)) {
    return NULL;
  }
  
  dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  if (dip != NULL)
  {
    dip->clickable_item_type = DISC_INCONSISTENT_BIOSRC_DEFLINE;
    dip->description = (CharPtr)MemNew (StringLen (bad_fmt) + StringLen (desc));
    sprintf (dip->description, bad_fmt, desc);
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    dip->item_list = NULL;
    ValNodeAddPointer (&(dip->item_list), OBJ_SEQDESC, biop_sdp);
    ValNodeAddPointer (&(dip->item_list), OBJ_SEQDESC, title_sdp);
  }      
  return dip;
}


static void FindInconsistentSourceAndDeflineCallback (BioseqPtr bsp, Pointer userdata)
{
  ClickableItemPtr dip;
  ValNodePtr PNTR discrepancy_list;
  SeqDescrPtr        biop_sdp, title_sdp;
  SeqMgrDescContext  context;
  BioSourcePtr       biop;
  
  discrepancy_list = (ValNodePtr PNTR) userdata;
  if (bsp == NULL || discrepancy_list == NULL) return;
  
  biop_sdp = SeqMgrGetNextDescriptor(bsp, NULL, Seq_descr_source, &context);
  if (biop_sdp == NULL || biop_sdp->data.ptrvalue == NULL)
  {
    return;
  }
  biop = (BioSourcePtr) biop_sdp->data.ptrvalue;
  if (biop->org == NULL)
  {
    return;
  }
  if (StringHasNoText (biop->org->taxname)) 
  {
    return;
  }
  
  title_sdp = SeqMgrGetNextDescriptor(bsp, NULL, Seq_descr_title, &context);
  if (title_sdp == NULL) return;
  
  if (StringStr (title_sdp->data.ptrvalue, biop->org->taxname) == NULL)
  {
    dip = InconsistentSourceDefline (biop_sdp, title_sdp);
    if (dip != NULL)
    {
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}


extern void FindInconsistentSourceAndDefline (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{  
  ValNodePtr disc_pairs = NULL, vnp;
  CharPtr    bad_fmt = "%d sources do not match definition lines.";
  ClickableItemPtr dip;
  
  if (discrepancy_list == NULL) return;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &disc_pairs, FindInconsistentSourceAndDeflineCallback);
  }

  if (disc_pairs == NULL) 
  {
    return;
  }
  else if (disc_pairs->next == NULL)
  {
    ValNodeLink (discrepancy_list, disc_pairs);
  }
  else
  {
    dip = NewClickableItem (DISC_INCONSISTENT_BIOSRC_DEFLINE, bad_fmt, disc_pairs);
    dip->item_list = NULL;
    dip->subcategories = disc_pairs;
    
    ValNodeAddPointer (discrepancy_list, 0, dip);
  }
}


static void FindParticalCDSsInCompleteSequencesCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR    cds_list;
  SeqDescrPtr        molinfo_sdp;
  SeqMgrDescContext  context;
  SeqFeatPtr         cds;
  SeqMgrFeatContext  fcontext;
  MolInfoPtr         mip;
  Boolean            partial5, partial3;
  
  cds_list = (ValNodePtr PNTR) userdata;
  if (bsp == NULL || cds_list == NULL) return;
  
  molinfo_sdp = SeqMgrGetNextDescriptor(bsp, NULL, Seq_descr_molinfo, &context);
  if (molinfo_sdp == NULL || molinfo_sdp->data.ptrvalue == NULL)
  {
    return;
  }
  mip = (MolInfoPtr) molinfo_sdp->data.ptrvalue;
  if (mip->completeness != 1)
  {
    return;
  }
  
  cds = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, FEATDEF_CDS, &fcontext);
  while (cds != NULL) {
      CheckSeqLocForPartial (cds->location, &partial5, &partial3);
      if (cds->partial || partial5 || partial3) {
          ValNodeAddPointer (cds_list, OBJ_SEQFEAT, cds);
      }
      cds = SeqMgrGetNextFeature (bsp, cds, SEQFEAT_CDREGION, FEATDEF_CDS, &fcontext);
  }
}


extern void FindParticalCDSsInCompleteSequences (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{  
  ValNodePtr cds_list = NULL, vnp;
  CharPtr    bad_fmt = "%d partial CDSs in complete sequences.";
  ClickableItemPtr dip;
  
  if (discrepancy_list == NULL) return;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &cds_list, FindParticalCDSsInCompleteSequencesCallback);
  }

  if (cds_list == NULL) 
  {
    return;
  }
  else
  {
    dip = NewClickableItem (DISC_PARTIAL_CDS_IN_COMPLETE_SEQUENCE, bad_fmt, cds_list);
    dip->subcategories = NULL;
    
    ValNodeAddPointer (discrepancy_list, 0, dip);
  }
}

static Boolean RnaRefMatch (RnaRefPtr rrp1, RnaRefPtr rrp2)
{
  tRNAPtr tp1, tp2;
  Boolean rval = FALSE;

  if (rrp1 == NULL && rrp2 == NULL) {
    rval = TRUE;
  } else if (rrp1 == NULL || rrp2 == NULL) {
    rval = FALSE;
  } else if (rrp1->type != rrp2->type) {
    rval = FALSE;
  } else if (rrp1->ext.choice != rrp2->ext.choice) {
    return FALSE;
  } else {
    switch (rrp1->ext.choice) {
      case 0:
        rval = TRUE;
        break;
      case 1:
        if (StringCmp (rrp1->ext.value.ptrvalue, rrp2->ext.value.ptrvalue) == 0) {
          rval = TRUE;
        } else {
          rval = FALSE;
        }
        break;
      case 2:
        tp1 = rrp1->ext.value.ptrvalue;
        tp2 = rrp2->ext.value.ptrvalue;
        if (tp1 == NULL && tp2 == NULL) {
          rval = TRUE;
        } else if (tp1 == NULL || tp2 == NULL) {
          rval = FALSE;
        } else if (tp1->aa == tp2->aa) {
          rval = TRUE;
        } else {
          rval = FALSE;
        }
        break;
      default:
        rval = FALSE;
        break;
    }
  }
  return rval;
}

static void AddRNAMatch (SeqFeatPtr sfp, ValNodePtr PNTR puniq_list)
{
  ValNodePtr vnp, uniq_vnp;
  SeqFeatPtr sfp_match;
  RnaRefPtr rrp_match, rrp_find;
  Boolean   found_match = FALSE;

  if (sfp == NULL || sfp->data.value.ptrvalue == NULL || puniq_list == NULL) return;
  rrp_find = (RnaRefPtr) sfp->data.value.ptrvalue;

  uniq_vnp = *puniq_list;

  if (uniq_vnp == NULL) {
    vnp = ValNodeNew(NULL);
    vnp->choice = OBJ_SEQFEAT;
    vnp->data.ptrvalue = sfp;
    vnp->next = NULL;
    ValNodeAddPointer (puniq_list, 0, vnp);
    found_match = TRUE;
  }
  while (uniq_vnp != NULL && !found_match) {
    vnp = uniq_vnp->data.ptrvalue;
    if (vnp == NULL) {
      /* fill in empty list */
      ValNodeAddPointer (&vnp, OBJ_SEQFEAT, sfp);
      uniq_vnp->data.ptrvalue = vnp;
      found_match = TRUE;
    } else {
      sfp_match = vnp->data.ptrvalue;
      if (sfp_match != NULL && sfp_match->data.choice == SEQFEAT_RNA && sfp_match->data.value.ptrvalue != NULL) {
        rrp_match = sfp_match->data.value.ptrvalue;
        if (RnaRefMatch(rrp_match, rrp_find)) {
          ValNodeAddPointer (&vnp, OBJ_SEQFEAT, sfp);
          found_match = TRUE;
          /* set flag so we know this list has duplicates */
          uniq_vnp->choice = 1;
        } 
      }
      if (!found_match) {
        if (uniq_vnp->next == NULL) {
          /* add to end of list */
          uniq_vnp->next = ValNodeNew(NULL);
          uniq_vnp->next->next = NULL;
          uniq_vnp->next->choice = 0;
          vnp = ValNodeNew(NULL);
          vnp->choice = OBJ_SEQFEAT;
          vnp->data.ptrvalue = sfp;
          vnp->next = NULL;
          uniq_vnp->next->data.ptrvalue = vnp;
          found_match = TRUE;
        } else {
          uniq_vnp = uniq_vnp->next;
        }
      }
    }
  }
}
  

static void FindDupRNAsInList (ValNodePtr rna_list, ValNodePtr PNTR discrepancy_list, CharPtr label, CharPtr id_str)
{
  ValNodePtr vnp, uniq_list = NULL;
  ValNodePtr dup_list = NULL;
  CharPtr          dup_fmt = "%d %s features on %s have the same name (%s)";
  ClickableItemPtr cip;
  SeqFeatPtr       sfp;
  SeqMgrFeatContext fcontext;

  for (vnp = rna_list; vnp != NULL; vnp = vnp->next) {
    AddRNAMatch (vnp->data.ptrvalue, &uniq_list);
  }

  dup_list = ValNodeExtractList (&uniq_list, 1);

  for (vnp = uniq_list; vnp != NULL; vnp = vnp->next) {
    uniq_list->data.ptrvalue = ValNodeFree (uniq_list->data.ptrvalue);
  }
  uniq_list = ValNodeFree (uniq_list);

  for (vnp = dup_list; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->item_list = vnp->data.ptrvalue;
    sfp = cip->item_list->data.ptrvalue;
    sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &fcontext);
    cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (dup_fmt) + StringLen (label) + StringLen (id_str) + StringLen (fcontext.label) + 15));
    sprintf (cip->description, dup_fmt, ValNodeLen (cip->item_list), label, id_str, fcontext.label);
    if (sfp->idx.subtype == FEATDEF_tRNA) {
      cip->clickable_item_type = DISC_DUP_TRNA;
    } else {
      cip->clickable_item_type = DISC_DUP_RRNA;
    }
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
  dup_list = ValNodeFree (dup_list);

}

typedef struct desiredaa {
  Char    short_symbol;
  CharPtr long_symbol;
  Int4    num_expected;
} DesiredAAData, PNTR DesiredAAPtr;

static DesiredAAData desired_aaList [] = {
{'A', "Ala", 1 },
{'B', "Asx", 0 },
{'C', "Cys", 1 },
{'D', "Asp", 1 },
{'E', "Glu", 1 },
{'F', "Phe", 1 },
{'G', "Gly", 1 },
{'H', "His", 1 },
{'I', "Ile", 1 },
{'J', "Xle", 0 },
{'K', "Lys", 1 },
{'L', "Leu", 2 },
{'M', "Met", 1 },
{'N', "Asn", 1 },
{'P', "Pro", 1 },
{'Q', "Gln", 1 },
{'R', "Arg", 1 },
{'S', "Ser", 2 },
{'T', "Thr", 1 },
{'V', "Val", 1 },
{'W', "Trp", 1 },
{'X', "Xxx", 0 },
{'Y', "Tyr", 1 },
{'Z', "Glx", 0 },
{'U', "Sec", 0 },
{'O', "Pyl", 0 },
{'*', "Ter", 0 }
};

static void AddMissingtRNADiscrepancy (CharPtr str, ValNodePtr PNTR discrepancy_list, CharPtr id_str, BioseqPtr bsp)
{
  ClickableItemPtr cip;
  CharPtr          desc_fmt = "Sequence %s is missing trna-%s";

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  MemSet (cip, 0, sizeof (ClickableItemData));

  cip->clickable_item_type = DISC_COUNT_TRNA;
  ValNodeAddPointer (&(cip->item_list), OBJ_BIOSEQ, bsp);
  cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (desc_fmt) + StringLen (id_str) + StringLen (str)));
  sprintf (cip->description, desc_fmt, id_str, str);
  ValNodeAddPointer (discrepancy_list, 0, cip);
}

static void 
AddExtratRNADiscrepancy 
(CharPtr         str, 
 Int4            num, 
 ValNodePtr PNTR discrepancy_list, 
 CharPtr         id_str, 
 BioseqPtr       bsp,
 ValNodePtr      rna_list)
{
  ClickableItemPtr cip;
  SeqMgrFeatContext fcontext;
  CharPtr          desc_fmt = "Sequence %s has %d trna-%s features";
  ValNodePtr       vnp;
  SeqFeatPtr       sfp;

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  MemSet (cip, 0, sizeof (ClickableItemData));

  cip->clickable_item_type = DISC_COUNT_TRNA;
  ValNodeAddPointer (&(cip->item_list), OBJ_BIOSEQ, bsp);
  for (vnp = rna_list; vnp != NULL; vnp = vnp->next) {
    sfp = vnp->data.ptrvalue;
    sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &fcontext);
    if (StringSearch (fcontext.label, str) != NULL) {
      ValNodeAddPointer (&(cip->item_list), OBJ_SEQFEAT, sfp);
    }
  }
  cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (desc_fmt) + StringLen (id_str) + StringLen (str) + 15));
  sprintf (cip->description, desc_fmt, id_str, num, str);
  ValNodeAddPointer (discrepancy_list, 0, cip);
}

static void FindMissingRNAsInList (ValNodePtr rna_list, ValNodePtr PNTR discrepancy_list, CharPtr id_str, BioseqPtr bsp)
{
  ValNodePtr       vnp;
  SeqFeatPtr       sfp;
  SeqMgrFeatContext fcontext;
  Uint1            num;
  Int4Ptr          num_present;
  Uint1            i;

  num = sizeof (desired_aaList) / sizeof (DesiredAAData);

  num_present = (Int4Ptr) MemNew (sizeof (Int4) * num);
  MemSet (num_present, 0, sizeof (Int4) * num);

  for (vnp = rna_list; vnp != NULL; vnp = vnp->next) {
    sfp = vnp->data.ptrvalue;
    sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &fcontext);
    for (i = 0; i < num; i++) {
      if (StringSearch (fcontext.label, desired_aaList[i].long_symbol) != NULL) {
        num_present[i] ++;
        break;
      }
    }
  }
  for (i = 0; i < num; i++) {
    if (num_present[i] < desired_aaList[i].num_expected) {
      AddMissingtRNADiscrepancy (desired_aaList[i].long_symbol, discrepancy_list, id_str, bsp);
    } else if (num_present[i] > desired_aaList[i].num_expected) {
      AddExtratRNADiscrepancy (desired_aaList[i].long_symbol, num_present[i], discrepancy_list, id_str, bsp, rna_list);
    }
  }
}

typedef struct featcount {
  Uint1      featdeftype;
  ValNodePtr discrepancy_list;
} FeatCountData, PNTR FeatCountPtr;

static void RNACountFeaturesBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqMgrFeatContext fcontext;
  SeqMgrDescContext dcontext;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  ValNodePtr         feat_list = NULL;
  BioSourcePtr       biop;
  Boolean            run_test = FALSE;
  FeatCountPtr       fcp;
  CharPtr            count_fmt = "%d %s features found on %s";
  CharPtr            label;
  ClickableItemPtr   cip;
  Char        id_str[45];

  if (bsp == NULL || ISA_aa (bsp->mol) || userdata == NULL) {
    return;
  }

  fcp = (FeatCountPtr) userdata;

  /* look for Bioseq with organelle */
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  while (sdp != NULL && !run_test) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL 
        && (biop->genome == GENOME_plastid
            || biop->genome == GENOME_mitochondrion
            || biop->genome == GENOME_chloroplast)) {
      run_test = TRUE;
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext);
  }
  if (!run_test) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, fcp->featdeftype, &fcontext);
  while (sfp != NULL) {
    ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_RNA, fcp->featdeftype, &fcontext);
  }

  if (feat_list != NULL) {
    label = (CharPtr) FeatDefTypeLabel(feat_list->data.ptrvalue);
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_str, PRINTID_REPORT, sizeof (id_str) - 1);
    cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (count_fmt) + StringLen (label) + StringLen (id_str) + 15));
    sprintf (cip->description, count_fmt, ValNodeLen (feat_list), label, id_str);
    cip->item_list = feat_list;
    if (fcp->featdeftype == FEATDEF_tRNA) {
      cip->clickable_item_type = DISC_COUNT_TRNA;
    } else {
      cip->clickable_item_type = DISC_COUNT_RRNA;
    }
    ValNodeAddPointer (&(fcp->discrepancy_list), 0, cip);
    if (fcp->featdeftype == FEATDEF_tRNA) {
      FindMissingRNAsInList (feat_list, &(fcp->discrepancy_list), id_str, bsp);
    } else {
      FindDupRNAsInList (feat_list, &(fcp->discrepancy_list), label, id_str);
    }
  }
}


static BioseqPtr GetRNATestBioseq (ValNodePtr vp)
{
  ClickableItemPtr cip;
  ValNodePtr       vnp;
  BioseqPtr        bsp = NULL;
  SeqFeatPtr       sfp;

  if (vp == NULL || vp->data.ptrvalue == NULL) return NULL;
  cip = (ClickableItemPtr) vp->data.ptrvalue;
  for (vnp = cip->item_list; vnp != NULL && bsp == NULL; vnp = vnp->next) {
    if (vnp->data.ptrvalue == NULL) continue;
    if (vnp->choice == OBJ_SEQFEAT) {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      bsp = BioseqFindFromSeqLoc (sfp->location);
    }
  }
  return bsp;
}

static void AddRNANumList (ValNodePtr PNTR discrepancy_list, ValNodePtr list_start)
{
  ClickableItemPtr cip;
  CharPtr          cp;
  CharPtr          desc_fmt = "%d sequences have ";
  CharPtr          desc_str;
  Int4             copy_len, orig_len;
  ValNodePtr       vnp;
  BioseqPtr        bsp;

  if (discrepancy_list == NULL || list_start == NULL) return;
  desc_str = GetClickableItemDescription (list_start);
  cp = StringSearch (desc_str, " found on");
  if (cp != NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->subcategories = list_start;
    copy_len = cp - desc_str;
    cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (desc_fmt) + 15 + copy_len));
    sprintf (cip->description, desc_fmt, ValNodeLen (list_start));
    orig_len = StringLen (cip->description);
    StringNCat (cip->description, desc_str, copy_len);
    cip->description [orig_len + copy_len] = 0;

    for (vnp = list_start; vnp != NULL; vnp = vnp->next) {
      bsp = GetRNATestBioseq (vnp);
      if (bsp != NULL) {
        ValNodeAddPointer (&(cip->item_list), OBJ_BIOSEQ, bsp);
      }
    }
    ValNodeAddPointer (discrepancy_list, 0, cip);
  } else {
    ValNodeLink (discrepancy_list, list_start);
  }
}

static void RNACountFeaturesAndFindDups (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list, Uint1 featdeftype)
{
  ValNodePtr       vnp, list_start, list_last;
  FeatCountData    fcd;
  SeqEntryPtr      sep;
  CharPtr          cp, compare1;
  Int4             compare_len;

  fcd.featdeftype = featdeftype;
  fcd.discrepancy_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitBioseqsInSep (sep, &fcd, RNACountFeaturesBioseqCallback);
  }

  /* count how many Bioseqs have different numbers of features */
  fcd.discrepancy_list = ValNodeSort (fcd.discrepancy_list, SortVnpByClickableItemDescription);

  while (fcd.discrepancy_list != NULL) {
    list_start = fcd.discrepancy_list;
    compare1 = GetClickableItemDescription (list_start);
    cp = StringSearch (compare1, "found on");
    if (cp == NULL) {
      fcd.discrepancy_list = fcd.discrepancy_list->next;
      list_start->next = NULL;
      AddRNANumList (discrepancy_list, list_start);
    } else {
      compare_len = cp - compare1;
      list_last = list_start;
      vnp = list_start->next;
      while (vnp != NULL && StringNCmp (compare1, GetClickableItemDescription (vnp), compare_len) == 0) {
        list_last = vnp;
        vnp = vnp->next;
      }
      
      list_last->next = NULL;
      fcd.discrepancy_list = vnp;
      AddRNANumList (discrepancy_list, list_start);
    }
  }
}

extern void tRNACountFeaturesAndFindDups (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  RNACountFeaturesAndFindDups (discrepancy_list, sep_list, FEATDEF_tRNA);
}

extern void rRNACountFeaturesAndFindDups (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  RNACountFeaturesAndFindDups (discrepancy_list, sep_list, FEATDEF_rRNA);
}

/* do not count short tRNAs if they are partial */
static void CountShorttRNA (SeqFeatPtr sfp, Pointer data)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_tRNA || data == NULL || sfp->partial) return;

  if (SeqLocLen (sfp->location) < 50) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
}

static void CountLongtRNA (SeqFeatPtr sfp, Pointer data)
{
  SeqMgrFeatContext fcontext;
  Int4 len;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_tRNA || data == NULL) return;

  if ((len = SeqLocLen (sfp->location)) > 90) {
    if (len > 100
        || (sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, 0, 0, sfp, &fcontext)) == NULL
        || (StringCmp (fcontext.label, "Ser") != 0
            && StringCmp (fcontext.label, "Leu") != 0
            && StringCmp (fcontext.label, "Sec") != 0)) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
    }
  }
}

extern void tRNAFindBadLength (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr  vnp;
  SeqEntryPtr sep;
  ValNodePtr  too_short = NULL, too_long = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitFeaturesInSep (sep, &too_short, CountShorttRNA);
  }
  if (too_short != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_BADLEN_TRNA, "%d tRNAs are too short", too_short));
  }

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitFeaturesInSep (sep, &too_long, CountLongtRNA);
  }
  if (too_long != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_BADLEN_TRNA, "%d tRNAs are too long", too_long));
  }

}


static void FindRNAsWithoutProductsCallback (SeqFeatPtr sfp, Pointer data)
{
  ValNode field;
  FeatureFieldPtr ff;
  CharPtr str;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) {
    return;
  }

  if (sfp->idx.subtype == FEATDEF_otherRNA) {
    if (StringNICmp (sfp->comment, "contains ", 9) == 0
        || StringNICmp (sfp->comment, "may contain", 11) == 0) {
      return;
    }
  } else if (sfp->idx.subtype == FEATDEF_tmRNA) {
    /* don't require products for tmRNA */
    return;
  }

  ff = FeatureFieldNew ();
  ff->type = Macro_feature_type_any;
  ValNodeAddInt (&ff->field, FeatQualChoice_legal_qual, Feat_qual_legal_product);
  field.choice = FieldType_feature_field;
  field.data.ptrvalue = ff;
  field.next = NULL;

  str = GetFieldValueForObject (OBJ_SEQFEAT, sfp, &field, NULL);
  if (StringHasNoText (str)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
  str = MemFree (str);
  ff = FeatureFieldFree (ff);
}


static ClickableItemPtr PseudoAndNonPseudoClickableItem (Uint4 clickable_item_type, CharPtr format, ValNodePtr item_list)
{
  ValNodePtr pseudo_list = NULL, non_pseudo_list = NULL, vnp;
  CharPtr pseudo_fmt = " and are pseudo", non_pseudo_fmt = " and are not pseudo";
  ClickableItemPtr cip, pseudo_cip = NULL, non_pseudo_cip = NULL;

  if (item_list == NULL) {
    return NULL;
  }

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT) {
      if (IsPseudo(vnp->data.ptrvalue)) {
        ValNodeAddPointer (&pseudo_list, OBJ_SEQFEAT, vnp->data.ptrvalue);
      } else {
        ValNodeAddPointer (&non_pseudo_list, OBJ_SEQFEAT, vnp->data.ptrvalue);
      }
    }
  }
  if (pseudo_list != NULL) {
    pseudo_cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (pseudo_cip, 0, sizeof (ClickableItemData));
    pseudo_cip->clickable_item_type = clickable_item_type;
    pseudo_cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (format) + StringLen (pseudo_fmt) + 15));
    sprintf (pseudo_cip->description, format, ValNodeLen (pseudo_list));
    StringCat (pseudo_cip->description, pseudo_fmt);
    pseudo_cip->item_list = pseudo_list;
  }

  if (non_pseudo_list != NULL) {
    non_pseudo_cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (non_pseudo_cip, 0, sizeof (ClickableItemData));
    non_pseudo_cip->clickable_item_type = clickable_item_type;
    non_pseudo_cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (format) + StringLen (non_pseudo_fmt) + 15));
    sprintf (non_pseudo_cip->description, format, ValNodeLen (non_pseudo_list));
    StringCat (non_pseudo_cip->description, non_pseudo_fmt);
    non_pseudo_cip->item_list = non_pseudo_list;
  }

  if (pseudo_cip == NULL) {
    cip = non_pseudo_cip;
  } else if (non_pseudo_cip == NULL) {
    cip = pseudo_cip;
  } else {
    cip = NewClickableItem (clickable_item_type, format, item_list);
    ValNodeAddPointer (&(cip->subcategories), 0, non_pseudo_cip);
    ValNodeAddPointer (&(cip->subcategories), 0, pseudo_cip);
  }
  return cip;
}


extern void FindRNAsWithoutProducts (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr       vnp;
  SeqEntryPtr      sep;
  ValNodePtr       rna_list = NULL;
  ClickableItemPtr cip;
  SeqEntryPtr      oldscope;

  oldscope = SeqEntrySetScope (NULL);
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    SeqEntrySetScope (sep);
    VisitFeaturesInSep (sep, &rna_list, FindRNAsWithoutProductsCallback);
  }
  SeqEntrySetScope (oldscope);
  if (rna_list != NULL) {
   cip = PseudoAndNonPseudoClickableItem (DISC_RNA_NO_PRODUCT, "%d RNA features have no product", rna_list);
   ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


static void tRNASameStrandBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqMgrFeatContext fcontext;
  SeqMgrDescContext dcontext;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  ValNodePtr         feat_list = NULL;
  BioSourcePtr       biop;
  Boolean            run_test = FALSE;
  ValNodePtr PNTR    discrepancy_list;
  ClickableItemPtr   cip;

  Uint1              strand, this_strand;
  Boolean            mixed_strand = FALSE;

  if (bsp == NULL || ISA_aa (bsp->mol) || userdata == NULL) {
    return;
  }

  discrepancy_list = (ValNodePtr PNTR) userdata;

  /* look for Bioseq with organelle */
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  while (sdp != NULL && !run_test) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL 
        && (biop->genome == GENOME_plastid
            || biop->genome == GENOME_mitochondrion
            || biop->genome == GENOME_chloroplast)) {
      run_test = TRUE;
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext);
  }
  if (!run_test) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_RNA, FEATDEF_tRNA, &fcontext);
  while (sfp != NULL && !mixed_strand) {
    if (feat_list == NULL) {
      strand = SeqLocStrand (sfp->location);
    } else {
      this_strand = SeqLocStrand (sfp->location);
      if ((strand == Seq_strand_minus && this_strand != Seq_strand_minus)
          || (strand != Seq_strand_minus && this_strand == Seq_strand_minus)) {
        mixed_strand = TRUE;
      }
    }
    ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_RNA, FEATDEF_tRNA, &fcontext);
  }

  if (mixed_strand) {
    feat_list = ValNodeFree (feat_list);
  } else if (feat_list != NULL) {
    if (strand == Seq_strand_minus) {
      cip = NewClickableItem (DISC_STRAND_TRNA, "%d tRNAs on minus strand", feat_list);
    } else {
      cip = NewClickableItem (DISC_STRAND_TRNA, "%d tRNAs on plus strand", feat_list);
    }
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


extern void FindtRNAsOnSameStrand (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr       vnp;
  SeqEntryPtr      sep;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitBioseqsInSep (sep, discrepancy_list, tRNASameStrandBioseqCallback);
  }
}


typedef struct translnote {
  ValNodePtr transl_no_note;
  ValNodePtr note_no_transl;
  ValNodePtr transl_too_long;
} TranslNoteData, PNTR TranslNotePtr;

NLM_EXTERN Boolean CodingRegionHasTranslExcept (SeqFeatPtr sfp)
{
  CodeBreakPtr cbp;
  Int4         len, tmp_len;
  CdRegionPtr  crp;
  SeqLocPtr    slp;
  Int4         codon_start, codon_stop, pos, codon_length;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION
      || (crp = (CdRegionPtr)sfp->data.value.ptrvalue) == NULL
      || crp->code_break == NULL)
  {
      return FALSE;
  }

  len = SeqLocLen (sfp->location);
  tmp_len = len;

  if (crp->frame == 2) {
    tmp_len -= 1;
  } else if (crp->frame == 3) {
    tmp_len -= 2;
  } 
  if (tmp_len % 3 == 0) 
  {
      return FALSE;
  }
  for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next)
  {
    if (cbp->aa.choice != 1 || cbp->aa.value.intvalue != 42) {
      continue;
    }
    codon_start = INT4_MAX;
    codon_stop = -10;
    slp = NULL;
    while ((slp = SeqLocFindNext (cbp->loc, slp)) != NULL) {
      pos = GetOffsetInLoc (slp, sfp->location, SEQLOC_START);
      if (pos <= codon_start)
      {
        codon_start = pos;
        pos = GetOffsetInLoc (slp, sfp->location, SEQLOC_STOP);
        if (pos > codon_stop)
        {
          codon_stop = pos;
        }
        codon_length = codon_stop - codon_start;      /* codon length */
        if (codon_length >= 0 && codon_length <= 1 && codon_stop == len - 1)
        {                       /*  a codon */
          /* allowing a partial codon at the end */
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean TranslTooLong (SeqFeatPtr sfp)
{
  CodeBreakPtr cbp;
  CdRegionPtr  crp;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION
      || (crp = (CdRegionPtr)sfp->data.value.ptrvalue) == NULL
      || crp->code_break == NULL)
  {
      return FALSE;
  }

  for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next)
  {
    if (cbp->aa.choice == 1 
        && cbp->aa.value.intvalue == 42
        && SeqLocLen (cbp->loc) > 3) {
      return TRUE;
    }
  }
  return FALSE;
}

static void FindTranslNoNote (SeqFeatPtr sfp, Pointer userdata)
{
  TranslNotePtr tnp;
  CharPtr       note_txt = "TAA stop codon is completed by the addition of 3' A residues to the mRNA";

  if (sfp != NULL && userdata != NULL && sfp->data.choice == SEQFEAT_CDREGION) {
    tnp = (TranslNotePtr) userdata;
    if (CodingRegionHasTranslExcept (sfp)) {
      if (StringStr (sfp->comment, note_txt) == NULL) {
        ValNodeAddPointer (&(tnp->transl_no_note), OBJ_SEQFEAT, sfp);
      }
    } else if (StringStr (sfp->comment, note_txt) != NULL) {
      ValNodeAddPointer (&(tnp->note_no_transl), OBJ_SEQFEAT, sfp);
    }
    if (TranslTooLong(sfp)) {
      ValNodeAddPointer (&(tnp->transl_too_long), OBJ_SEQFEAT, sfp);
    }
  }
}

extern void FindTranslExceptNotes (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr       vnp;
  TranslNoteData   tnd;
  SeqEntryPtr      sep;
  ClickableItemPtr cip;
  CharPtr          transl_no_note_fmt = "%d features have a translation exception but no note";
  CharPtr          note_no_transl_fmt = "%d features have a note but not translation exception";
  CharPtr          transl_too_long_fmt = "%d features have translation exceptions longer than 3 bp";

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    tnd.transl_no_note = NULL;
    tnd.note_no_transl = NULL;
    tnd.transl_too_long = NULL;
    VisitFeaturesInSep (sep, &tnd, FindTranslNoNote);
    if (tnd.transl_no_note != NULL) {
      cip = NewClickableItem (DISC_TRANSL_NO_NOTE, transl_no_note_fmt, tnd.transl_no_note);
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
    if (tnd.note_no_transl != NULL) {
      cip = NewClickableItem (DISC_NOTE_NO_TRANSL, note_no_transl_fmt, tnd.note_no_transl);
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
    if (tnd.transl_too_long != NULL) {
      cip = NewClickableItem (DISC_TRANSL_TOO_LONG, transl_too_long_fmt, tnd.note_no_transl);
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
  }
}

static Boolean GetOverlappingTRNAs (BioseqPtr bsp, SeqLocPtr slp, Int4 loc_right, ValNodePtr PNTR list)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext  context;
  Boolean            found_any = FALSE;
  Uint1              slp_strand, rna_strand;

  if (bsp == NULL || slp == NULL || list == NULL) return FALSE;
  slp_strand = SeqLocStrand (slp);

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_RNA, FEATDEF_tRNA, &context);
       sfp != NULL && context.left <= loc_right;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_RNA, FEATDEF_tRNA, &context))
  {
    rna_strand = SeqLocStrand (sfp->location);
    if (((slp_strand == Seq_strand_minus && rna_strand == Seq_strand_minus)
         || (slp_strand != Seq_strand_minus && rna_strand != Seq_strand_minus))
        && SeqLocCompare (sfp->location, slp) != SLC_NO_MATCH) {
      ValNodeAddPointer (list, OBJ_SEQFEAT, sfp);
      found_any = TRUE;
    }
  }
  return found_any;
}

static void FindCDSOverlappingtRNAsBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr         sfp;
  SeqMgrFeatContext  context;
  ValNodePtr         subcategories = NULL;
  ValNodePtr PNTR    discrepancy_list;
  ValNodePtr         item_list, all_item_list = NULL;
  ValNodePtr         trna_list = NULL;
  ClickableItemPtr   cip;
  CharPtr            list_fmt = "%d coding regions have overlapping tRNAs";
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  discrepancy_list = (ValNodePtr PNTR) userdata;
  
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, FEATDEF_CDS, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, FEATDEF_CDS, &context))
  {
    
    item_list = NULL;
    trna_list = NULL;
    if (GetOverlappingTRNAs (bsp, sfp->location, context.right, &trna_list)) {
      ValNodeAddPointer (&item_list, OBJ_SEQFEAT, sfp);
      ValNodeLink (&item_list, trna_list);
      ValNodeLink (&all_item_list, ValNodePointerDup(item_list));
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      MemSet (cip, 0, sizeof (ClickableItemData));
      cip->item_list = item_list;
      cip->description = StringSave ("Coding region overlaps tRNAs");
      cip->clickable_item_type = DISC_CDS_OVERLAP_TRNA;
      ValNodeAddPointer (&subcategories, 0, cip);
    }
  }  
  if (subcategories != NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof(ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->clickable_item_type = DISC_CDS_OVERLAP_TRNA;
    cip->item_list = all_item_list;
    cip->subcategories = subcategories;
    cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (list_fmt) + 15));
    sprintf (cip->description, list_fmt, ValNodeLen (subcategories));
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


extern void FindCDSOverlappingtRNAs (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr       vnp, this_list = NULL;
  SeqEntryPtr      sep;
  ClickableItemPtr cip;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitBioseqsInSep (sep, &this_list, FindCDSOverlappingtRNAsBioseqCallback);
  }
  if (this_list != NULL) {
    cip = NewClickableItem (DISC_CDS_OVERLAP_TRNA, "%d Bioseqs have coding regions that overlap tRNAs", this_list);
    cip->subcategories = this_list;
    cip->item_list = ItemListFromSubcategories (this_list);
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


static void RemoveFeaturesInsideLocation (ValNodePtr PNTR item_list, SeqLocPtr slp)
{
  ValNodePtr vnp, vnp_prev = NULL, vnp_next;
  SeqFeatPtr sfp;
  Boolean    do_remove;
  Int2       cmp;

  if (item_list == NULL || slp == NULL) {
    return;
  }

  for (vnp = *item_list; vnp != NULL; vnp = vnp_next) {
    vnp_next = vnp->next;
    do_remove = FALSE;
    if (vnp->choice == OBJ_SEQFEAT) {
      sfp = vnp->data.ptrvalue;
      if (sfp == NULL) {
        do_remove = TRUE;
      } else {
        cmp = SeqLocCompare (sfp->location, slp);
        if (cmp == SLC_A_IN_B || cmp == SLC_A_EQ_B) {
          do_remove = TRUE;
        }
      }
    }
    if (do_remove) {
      if (vnp_prev == NULL) {
        *item_list = vnp->next;
      } else {
        vnp_prev->next = vnp->next;
      }
      vnp->next = NULL;
      vnp = ValNodeFree (vnp);
    } else {
      vnp_prev = vnp;
    }
  }
}


static void FindFeaturesOverlappingSrcFeaturesBioseqCallback (BioseqPtr bsp, Pointer data)
{
  ClickableItemPtr  cip;
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;
  ValNodePtr        this_list, src_vnp;

  if (bsp == NULL || data == NULL) return;


  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_BIOSRC, &fcontext);
  while (sfp != NULL)
  {
    this_list = ListFeaturesOverlappingLocation (bsp, sfp->location, 0, 0);
    RemoveFeaturesInsideLocation (&this_list, sfp->location);

    if (this_list != NULL)
    {
      cip = NewClickableItem (DISC_FEAT_OVERLAP_SRCFEAT, "%d features overlap a source feature", this_list);
      /* insert source feature at beginning of item list */
      src_vnp = ValNodeNew (NULL);     
      src_vnp->choice = OBJ_SEQFEAT;
      src_vnp->data.ptrvalue = sfp;
      src_vnp->next = cip->item_list;
      cip->item_list = src_vnp;
      ValNodeAddPointer ((ValNodePtr PNTR) data, 0, cip);
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_BIOSRC, &fcontext);
  }
}


extern void FindFeaturesOverlappingSrcFeatures (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr       vnp;
  SeqEntryPtr      sep;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitBioseqsInSep (sep, discrepancy_list, FindFeaturesOverlappingSrcFeaturesBioseqCallback);
  }
}


static void LIBCALLBACK CountNProc (CharPtr sequence, Pointer userdata)
{
  Int4Ptr p_i;
  CharPtr cp;

  if (sequence == NULL || userdata == NULL) return;
  p_i = (Int4Ptr) userdata;

  for (cp = sequence; *cp != 0; cp++)
  {
    if (*cp == 'N')
    {
      (*p_i) ++;
    }
  }
}


NLM_EXTERN FloatLo PercentNInBioseq (BioseqPtr bsp, Boolean include_gaps)
{
  Int4 num_n = 0;
  Int4 flags = 0;
  
  if (bsp->length == 0) return 0;

  if (include_gaps) {
    flags |= STREAM_EXPAND_GAPS;
  }

  SeqPortStream (bsp, flags, (Pointer) &num_n, CountNProc);

  return ((FloatLo)num_n * 100) / (FloatLo) bsp->length;
}


static Boolean IsDeltaSeqWithFarpointers (BioseqPtr bsp)
{
  DeltaSeqPtr dsp;
  Boolean rval = FALSE;

  if (bsp == NULL || bsp->repr != Seq_repr_delta) {
    return FALSE;
  }
  for (dsp = (DeltaSeqPtr) (bsp->seq_ext); dsp != NULL && !rval; dsp = dsp->next) {
    if (dsp->choice == 1) {
      rval = TRUE;
    }
  }
  return rval;
}


static void PercentNDiscrepancy (BioseqPtr bsp, Pointer userdata)
{
  FloatLo pct;

  if (bsp == NULL || ISA_aa (bsp->mol) || userdata == NULL || IsDeltaSeqWithFarpointers (bsp))
  {
    return;
  }

  pct = PercentNInBioseq (bsp, FALSE);
  if (pct > 5.0) 
  {
    ValNodeAddPointer ((ValNodePtr PNTR)userdata, OBJ_BIOSEQ, bsp);
  }
}


static void PercentNDiscrepanciesForSeqEntry (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  SeqEntryPtr      sep;
  ValNodePtr       vnp, list = NULL;
  CharPtr top_fmt = "%d sequences have > 5%% Ns";

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = (SeqEntryPtr) vnp->data.ptrvalue;
    VisitBioseqsInSep (sep, &list, PercentNDiscrepancy);
  }

  if (list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_PERCENTN, top_fmt, list));
  }
}


typedef struct basecountandnrun {
  ValNodePtr no_a;
  ValNodePtr no_t;
  ValNodePtr no_c;
  ValNodePtr no_g;
  ValNodePtr n_run;
} BaseCountAndNRunData, PNTR BaseCountAndNRunPtr;


typedef struct basecounts {
  Int4 num_a;
  Int4 num_t;
  Int4 num_g;
  Int4 num_c;
  Int4 n_run;
  Boolean has_n_run;
  Int4 n_run_start;
  Int4 pos;
  ValNodePtr run_locations;
} BaseCountsData, PNTR BaseCountsPtr;

typedef struct intervalpair {
  Int4 start;
  Int4 stop;
} IntervalPairData, PNTR IntervalPairPtr;

static IntervalPairPtr IntervalPairNew (Int4 start, Int4 stop)
{
  IntervalPairPtr i;

  i = (IntervalPairPtr) MemNew (sizeof (IntervalPairData));
  i->start = start;
  i->stop = stop;
  return i;
}


static IntervalPairPtr IntervalPairFree (IntervalPairPtr i)
{
  i = MemFree (i);
  return i;
}


static void LIBCALLBACK CountBaseProc (CharPtr sequence, Pointer userdata)
{
  BaseCountsPtr counts;
  CharPtr cp;

  if (sequence == NULL || userdata == NULL) return;
  counts = (BaseCountsPtr) userdata;

  for (cp = sequence; *cp != 0; cp++, counts->pos ++)
  {
    if (*cp == 'N')
    {
      if (counts->n_run == 0) {
        counts->n_run_start = counts->pos;
      }
      counts->n_run ++;
    }
    else
    {
      if (counts->n_run >= 10)    /* 20->10, per Larissa's request, by J. Chen  */
      {
        counts->has_n_run = TRUE;
        ValNodeAddPointer (&(counts->run_locations), 0, IntervalPairNew (counts->n_run_start, counts->pos - 1));
      }
      counts->n_run = 0;
      switch (*cp)
      {
        case 'A':
          counts->num_a++;
          break;
        case 'T':
          counts->num_t++;
          break;
        case 'G':
          counts->num_g++;
          break;
        case 'C':
          counts->num_c++;
          break;
      }
    }
  }
}


static CharPtr FormatIntervalListString (ValNodePtr interval_list)
{
  CharPtr interval = NULL, interval_fmt = "%d-%d", cp;
  IntervalPairPtr i;
  ValNodePtr vnp;
  Int4       num;

  num = ValNodeLen (interval_list);
  if (num > 0) {
    interval = (CharPtr) MemNew (sizeof (Char) * ((StringLen (interval_fmt) + 30) * num));
    cp = interval;
    for (vnp = interval_list; vnp != NULL; vnp = vnp->next) {
      i = (IntervalPairPtr) vnp->data.ptrvalue;
      if (i != NULL) {
        sprintf (cp, interval_fmt, i->start + 1, i->stop + 1);
        StringCat (cp, ", ");
        cp += StringLen (cp);
      }
    }
    /* remove terminal comma and space */
    interval[StringLen(interval) - 2] = 0;
  }

  return interval;
}


static void BaseCountAndNRunDiscrepancyForBioseq (BioseqPtr bsp, Pointer userdata)
{
  BaseCountsData base_counts;
  BaseCountAndNRunPtr errs;
  ClickableItemPtr    cip;
  CharPtr             fmt = "%s has runs of Ns at the following locations: %s";
  CharPtr             interval;
  Char                id_buf[255];

  if (bsp == NULL || ISA_aa (bsp->mol) || userdata == NULL || IsDeltaSeqWithFarpointers (bsp)) return;
  errs = (BaseCountAndNRunPtr) userdata;
  MemSet (&base_counts, 0, sizeof (BaseCountsData));
  SeqPortStream (bsp, 0, (Pointer) &base_counts, CountBaseProc);
  if (base_counts.num_a == 0) {
    ValNodeAddPointer (&(errs->no_a), OBJ_BIOSEQ, bsp);
  }
  if (base_counts.num_c == 0) {
    ValNodeAddPointer (&(errs->no_c), OBJ_BIOSEQ, bsp);
  }
  if (base_counts.num_t == 0) {
    ValNodeAddPointer (&(errs->no_t), OBJ_BIOSEQ, bsp);
  }
  if (base_counts.num_g == 0) {
    ValNodeAddPointer (&(errs->no_g), OBJ_BIOSEQ, bsp);
  }
  if (base_counts.n_run >= 10) {   /* 20->10: per Larissa's request, by J. Chen */
    ValNodeAddPointer (&(base_counts.run_locations), 0, IntervalPairNew (base_counts.n_run_start, base_counts.pos - 1));
    base_counts.has_n_run = TRUE;
 }
 if (base_counts.has_n_run) { 
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    cip->clickable_item_type = DISC_N_RUNS;
    ValNodeAddPointer (&(cip->item_list), OBJ_BIOSEQ, bsp);
    SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_REPORT, sizeof (id_buf) - 1);
    interval = FormatIntervalListString (base_counts.run_locations);
    cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (id_buf) + StringLen (interval) + 1));
    sprintf (cip->description, fmt, id_buf, interval);
    interval = MemFree (interval);
    base_counts.run_locations = ValNodeFreeData (base_counts.run_locations);
    ValNodeAddPointer (&(errs->n_run), 0, cip);
  }
}


static void BaseCountAndNRunDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  SeqEntryPtr sep;
  ValNodePtr  vnp, zero_base = NULL, zero_base_tot_list = NULL, item_list;
  BaseCountAndNRunData lists;
  ClickableItemPtr cip;

  MemSet (&lists, 0, sizeof (BaseCountAndNRunData));
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitBioseqsInSep (sep, &lists, BaseCountAndNRunDiscrepancyForBioseq);
  }

  if (lists.n_run != NULL) {
    item_list = ItemListFromSubcategories (lists.n_run);

    /* 20->10: per Larissa's request, by J. Chen */
    cip = NewClickableItem (DISC_N_RUNS, "%d sequences have runs of 10 or more Ns", item_list);

    cip->subcategories = lists.n_run;
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }

  if (lists.no_a != NULL) {
    ValNodeAddPointer (&zero_base, 0, NewClickableItem (DISC_ZERO_BASECOUNT, "%d sequences have no As", lists.no_a));
    ValNodeLinkCopy (&zero_base_tot_list, lists.no_a);
  }
  if (lists.no_t != NULL) {
    ValNodeAddPointer (&zero_base, 0, NewClickableItem (DISC_ZERO_BASECOUNT, "%d sequences have no Ts", lists.no_t));
    ValNodeLinkCopy (&zero_base_tot_list, lists.no_t);
  }
  if (lists.no_g != NULL) {
    ValNodeAddPointer (&zero_base, 0, NewClickableItem (DISC_ZERO_BASECOUNT, "%d sequences have no Gs", lists.no_g));
    ValNodeLinkCopy (&zero_base_tot_list, lists.no_g);
  }
  if (lists.no_c != NULL) {
    ValNodeAddPointer (&zero_base, 0, NewClickableItem (DISC_ZERO_BASECOUNT, "%d sequences have no Cs", lists.no_c));
    ValNodeLinkCopy (&zero_base_tot_list, lists.no_c);
  }
    
  if (zero_base_tot_list != NULL) {
    cip = NewClickableItem (DISC_ZERO_BASECOUNT, "%d sequences have a zero basecount for a nucleotide", zero_base_tot_list);
    cip->subcategories = zero_base;
    ValNodeAddPointer (discrepancy_list, 0, cip);
  } 
}


CharPtr discReportDuplicateProteinIDFmt = "%d coding regions have non-unique protein IDs";
CharPtr discReportOneDuplicateProteinIDFmt = "%d coding regions have protein ID %s";
CharPtr discReportMissingProteinIDFmt = "%d coding regions have missing protein IDs";
CharPtr discReportDuplicateTranscriptIdFmt = "%d mRNAs have non-unique transcript IDs";
CharPtr discReportOneDuplicateTranscriptIdFmt = "%d mRNAs have non-unique transcript ID %s";
CharPtr discReportMissingTranscriptIDFmt = "%d mRNAs have missing transcript IDs";


/* look for duplicate protein IDs and duplicate transcript IDs */
/* every coding region should have a protein ID and a transcript ID */
/* RNA should have a transcript ID to match. */
static void CheckGenProdSetBioseq (BioseqPtr bsp, GenProdSetDiscrepancyListsPtr lists)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;
  SeqIdPtr          sip;
  Char              buf [96];

  if (bsp == NULL || !ISA_na (bsp->mol) || lists == NULL) {
    return;
  }

  /* look for missing protein IDs and duplicate protein IDs on coding regions */
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext)) {
    if (sfp->product == NULL) {
      if (!sfp->pseudo) {
        ValNodeAddPointer (&(lists->missing_protein_id), 0,
                           GlobalDiscrepancyNew (NULL, OBJ_SEQFEAT, sfp));
      }
    } else {
      sip = SeqLocId (sfp->product);
      SeqIdWrite (sip, buf, PRINTID_REPORT, sizeof (buf) - 1);
      ValNodeAddPointer (&(lists->cds_product_list), 0,
                         GlobalDiscrepancyNew (buf, OBJ_SEQFEAT, sfp));
    }
  }

  /* look for missing transcript IDs and duplicate transcript IDs on mRNAs */
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_mRNA, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_mRNA, &fcontext)) {
    if (sfp->product == NULL) {
      ValNodeAddPointer (&(lists->missing_mrna_product), 0, 
                         GlobalDiscrepancyNew (NULL, OBJ_SEQFEAT, sfp));
    } else {
      sip = SeqLocId (sfp->product);
      SeqIdWrite (sip, buf, PRINTID_REPORT, sizeof (buf) - 1);
      ValNodeAddPointer (&(lists->mrna_product_list), 0,
                         GlobalDiscrepancyNew (buf, OBJ_SEQFEAT, sfp));
    }
  }
}


extern void CheckGenProdSetsInSeqEntry (SeqEntryPtr sep, GenProdSetDiscrepancyListsPtr lists)
{
  BioseqSetPtr bssp;

  if (sep == NULL || !IS_Bioseq_set (sep) || sep->data.ptrvalue == NULL || lists == NULL) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
    if (IS_Bioseq (bssp->seq_set)) {
      CheckGenProdSetBioseq(bssp->seq_set->data.ptrvalue, lists);
    }
  } else {
    sep = bssp->seq_set;
    while (sep != NULL) {
      CheckGenProdSetsInSeqEntry (sep, lists);
      sep = sep->next;
    }
  }
}


static void CheckListForGenProdSets (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr       vnp, disc_list = NULL;
  ClickableItemPtr cip;
  GenProdSetDiscrepancyListsData lists;

  MemSet (&lists, 0, sizeof (GenProdSetDiscrepancyListsData));
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    CheckGenProdSetsInSeqEntry (vnp->data.ptrvalue, &lists);
  }

  if (lists.missing_protein_id != NULL) {
    cip = ReportMissingFields (lists.missing_protein_id, discReportMissingProteinIDFmt, DISC_MISSING_GENPRODSET_PROTEIN);
    if (cip != NULL) {
      ValNodeAddPointer (&disc_list, 0, cip);
    }
    lists.missing_protein_id = FreeGlobalDiscrepancyList (lists.missing_protein_id);
  }

  if (lists.cds_product_list != NULL) {
    lists.cds_product_list = ValNodeSort (lists.cds_product_list, SortVnpByGlobalDiscrepancyString);
    cip = ReportNonUniqueGlobalDiscrepancy (lists.cds_product_list, 
                                            discReportDuplicateProteinIDFmt,
                                            discReportOneDuplicateProteinIDFmt,
                                            DISC_DUP_GENPRODSET_PROTEIN,
                                            FALSE);
    if (cip != NULL) {
      ValNodeAddPointer (&disc_list, 0, cip);
    }
    lists.cds_product_list = FreeGlobalDiscrepancyList (lists.cds_product_list);
  }
  
  
  if (lists.missing_mrna_product != NULL) {
    cip = ReportMissingFields (lists.missing_mrna_product, discReportMissingTranscriptIDFmt, DISC_MISSING_GENPRODSET_TRANSCRIPT_ID);
    if (cip != NULL) {
      ValNodeAddPointer (&disc_list, 0, cip);
    }
    lists.missing_mrna_product = FreeGlobalDiscrepancyList (lists.missing_mrna_product);
  }


  if (lists.mrna_product_list != NULL) {
    lists.mrna_product_list = ValNodeSort (lists.mrna_product_list, SortVnpByGlobalDiscrepancyString);
    cip = ReportNonUniqueGlobalDiscrepancy (lists.mrna_product_list, 
                                            discReportDuplicateTranscriptIdFmt,
                                            discReportOneDuplicateTranscriptIdFmt,
                                            DISC_DUP_GENPRODSET_TRANSCRIPT_ID,
                                            FALSE);
    if (cip != NULL) {
      ValNodeAddPointer (&disc_list, 0, cip);
    }
    lists.mrna_product_list = FreeGlobalDiscrepancyList (lists.mrna_product_list);
  }
  
    


  if (disc_list != NULL) {
    if (disc_list->next == NULL) {
      ValNodeLink (discrepancy_list, disc_list);
    } else {
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      cip->description = StringSave ("GenProdSet Errors");
      cip->subcategories = disc_list;
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
  }
}


static void CountProteinsBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  if (bsp != NULL && ISA_aa (bsp->mol) && userdata != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_BIOSEQ, bsp);
  }
}

extern void CountProteins (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr       vnp;
  SeqEntryPtr      sep;
  ValNodePtr       proteins;
  ClickableItemPtr cip;
  CharPtr          prot_count_fmt = "%d protein sequences in record";

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    proteins = NULL;
    VisitBioseqsInSep (sep, &proteins, CountProteinsBioseqCallback);
    if (proteins != NULL) {
      cip = NewClickableItem (DISC_COUNT_PROTEINS, prot_count_fmt, proteins);
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
  }
}


typedef struct missingviralqualsdata {
  ValNodePtr missing_collection_date;
  ValNodePtr missing_country;
  ValNodePtr missing_specific_host;
} MissingViralQualsData, PNTR MissingViralQualsPtr;


static void AddMissingViralQualsDiscrepancies (BioSourcePtr biop, Uint1 choice, Pointer data, MissingViralQualsPtr q)
{
  SubSourcePtr ssp;
  OrgModPtr    mod;
  Boolean has_collection_date = FALSE;
  Boolean has_country = FALSE;
  Boolean has_specific_host = FALSE;

  if (!IsViralBioSource(biop) || q == NULL) {
    return;
  }

  for (ssp = biop->subtype; ssp != NULL && (!has_collection_date || !has_country); ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_collection_date) {
      has_collection_date = TRUE;
    } else if (ssp->subtype == SUBSRC_country) {
      has_country = TRUE;
    }
  }

  for (mod = biop->org->orgname->mod; mod != NULL && !has_specific_host; mod = mod->next) {
    if (mod->subtype == ORGMOD_nat_host) {
      has_specific_host = TRUE;
    }
  }

  if (!has_collection_date) {
    ValNodeAddPointer (&(q->missing_collection_date), choice, data);
  }
  if (!has_country) {
    ValNodeAddPointer (&(q->missing_country), choice, data);
  }
  if (!has_specific_host) {
    ValNodeAddPointer (&(q->missing_specific_host), choice, data);
  }
}


static void FindMissingViralQualsFeatCallback (SeqFeatPtr sfp, Pointer data)
{
  if (sfp != NULL 
      && sfp->data.choice == SEQFEAT_BIOSRC
      && data != NULL) {
    AddMissingViralQualsDiscrepancies (sfp->data.value.ptrvalue, OBJ_SEQFEAT, sfp, (MissingViralQualsPtr) data);
  }
}


static void FindMissingViralQualsDescCallback (SeqDescrPtr sdp, Pointer data)
{
  if (sdp != NULL 
      && sdp->choice == Seq_descr_source
      && data != NULL) {
    AddMissingViralQualsDiscrepancies (sdp->data.ptrvalue, OBJ_SEQDESC, sdp, (MissingViralQualsPtr) data);
  }
}


static void FindMissingViralQuals (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr       vnp;
  SeqEntryPtr      sep;
  MissingViralQualsData q;
  ClickableItemPtr cip;
  ValNodePtr       subcategories = NULL, item_list;

  q.missing_collection_date = NULL;
  q.missing_country = NULL;
  q.missing_specific_host = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitFeaturesInSep (sep, &q, FindMissingViralQualsFeatCallback);
    VisitDescriptorsInSep (sep, &q, FindMissingViralQualsDescCallback);
  }
  if (q.missing_collection_date != NULL) {
    ValNodeAddPointer (&subcategories, 0, NewClickableItem (DISC_MISSING_VIRAL_QUALS, "%d virus organisms are missing suggested qualifier collection date", q.missing_collection_date));
  }
  if (q.missing_country != NULL) {
    ValNodeAddPointer (&subcategories, 0, NewClickableItem (DISC_MISSING_VIRAL_QUALS, "%d virus organisms are missing suggested qualifier country", q.missing_country));
  }
  if (q.missing_specific_host != NULL) {
    ValNodeAddPointer (&subcategories, 0, NewClickableItem (DISC_MISSING_VIRAL_QUALS, "%d virus organisms are missing suggested qualifier specific-host", q.missing_specific_host));
  }
  if (subcategories != NULL) {
    item_list = ItemListFromSubcategories (subcategories);
    RemoveDuplicateItems (&item_list);
    cip = NewClickableItem (DISC_MISSING_VIRAL_QUALS, "%d virus organisms are missing required qualifiers", item_list);
    cip->subcategories = subcategories;
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


typedef struct duplicatequal {
  Uint1 choice;
  Pointer data;
  CharPtr val;
  ValNodePtr qual;
} DuplicateQualData, PNTR DuplicateQualPtr;

static DuplicateQualPtr DuplicateQualNew (Uint1 choice, Pointer data, ValNodePtr qual)
{
  DuplicateQualPtr dq;
  SourceQualChoicePtr s;

  dq = (DuplicateQualPtr) MemNew (sizeof (DuplicateQualData));
  dq->choice = choice;
  dq->data = data;
  dq->qual = AsnIoMemCopy (qual, (AsnReadFunc) FieldTypeAsnRead, (AsnWriteFunc) FieldTypeAsnWrite);
  dq->val = GetFieldValueForObject (choice, data, dq->qual, NULL);
  if (StringHasNoText (dq->val) && dq->qual != NULL && dq->qual->choice == FieldType_source_qual
      && (s = (SourceQualChoicePtr) dq->qual->data.ptrvalue) != NULL
      && s->choice == SourceQualChoice_location) {
    dq->val = MemFree (dq->val);
    dq->val = StringSave ("genomic");
  }
  return dq;
}


static void AddFieldValueToDuplicateQual (DuplicateQualPtr dq, ValNodePtr qual)
{
  CharPtr new_val, tmp;
  if (dq == NULL || qual == NULL) {
    return;
  }

  tmp = GetFieldValueForObject (dq->choice, dq->data, qual, NULL);
  if (!StringHasNoText (tmp)) {
    new_val = (CharPtr) MemNew (sizeof (Char) * (StringLen (dq->val) + StringLen (tmp) + 2));
    StringCpy (new_val, dq->val);
    StringCat (new_val, " ");
    StringCat (new_val, tmp);
    dq->val = MemFree (dq->val);
    dq->val = new_val;
  }
  tmp = MemFree (tmp);
}


static int CompareDuplicateQual (DuplicateQualPtr dq1, DuplicateQualPtr dq2)
{
  int         rval = 0;

  if (dq1 != NULL && dq2 != NULL) {
    rval = CompareFieldTypes (dq1->qual, dq2->qual);
    if (rval == 0) {
      rval = StringCmp (dq1->val, dq2->val);
    }
  }
  return rval;
}


static DuplicateQualPtr DuplicateQualFree (DuplicateQualPtr dq)
{
  if (dq != NULL) {
    dq->qual = FieldTypeFree (dq->qual);
    dq->val = MemFree (dq->val);
    dq = MemFree (dq);
  }
  return dq;
}


static ValNodePtr DuplicateQualListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = DuplicateQualFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static int LIBCALLBACK SortVnpByDuplicateQualFieldTypeThenValue (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);

    if (vnp1->data.ptrvalue != NULL && vnp2->data.ptrvalue != NULL) {
      rval = CompareDuplicateQual (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    }
  }

  return rval;
}


static int LIBCALLBACK SortVnpByDuplicateQualObjectThenValue (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  DuplicateQualPtr dq1, dq2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);

    if (vnp1->data.ptrvalue != NULL && vnp2->data.ptrvalue != NULL) {
      dq1 = vnp1->data.ptrvalue;
      dq2 = vnp2->data.ptrvalue;
      if (dq1->choice < dq2->choice) {
        rval = -1;
      } else if (dq1->choice > dq2->choice) {
        rval = 1;
      } else if (dq1->data < dq2->data) {
        rval = -1;
      } else if (dq1->data > dq2->data) {
        rval = 1;
      } else {
        rval = StringCmp (dq1->val, dq2->val);
      }
    }
  }

  return rval;
}


static Boolean IsCollectedByQual (FieldTypePtr qual)
{
  SourceQualChoicePtr sq;

  if (qual == NULL || qual->choice != FieldType_source_qual || qual->data.ptrvalue == NULL) {
    return FALSE;
  }
  sq = (SourceQualChoicePtr) qual->data.ptrvalue;
  if (sq->choice != SourceQualChoice_textqual) {
    return FALSE;
  } else if (sq->data.intvalue == Source_qual_collected_by) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean IsIdentifiedByQual (FieldTypePtr qual)
{
  SourceQualChoicePtr sq;

  if (qual == NULL || qual->choice != FieldType_source_qual || qual->data.ptrvalue == NULL) {
    return FALSE;
  }
  sq = (SourceQualChoicePtr) qual->data.ptrvalue;
  if (sq->choice != SourceQualChoice_textqual) {
    return FALSE;
  } else if (sq->data.intvalue == Source_qual_identified_by) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void 
ReportSameValueMultipleQuals 
(ValNodePtr PNTR discrepancy_list,
 Uint1           choice,
 Pointer         data, 
 CharPtr         val, 
 ValNodePtr      qual_list,
 Uint4           item_type)
{
  ValNodePtr vnp, name_list = NULL;
  Boolean    do_report = FALSE;
  CharPtr    qual_name, fmt = "BioSource has value '%s' for these qualifiers: ";
  ClickableItemPtr cip;
  Int4       names_len = 0;

  if (discrepancy_list == NULL || StringHasNoText (val) || qual_list == NULL || qual_list->next == NULL) {
    return;
  }

  /* make sure we have quals that are not collected by and identified by */
  for (vnp = qual_list; vnp != NULL; vnp = vnp->next) {
    if (!IsCollectedByQual (vnp->data.ptrvalue) && !IsIdentifiedByQual (vnp->data.ptrvalue)) {
      do_report = TRUE;
    }
    qual_name = SummarizeFieldType (vnp->data.ptrvalue);
    names_len += StringLen (qual_name) + 2;
    ValNodeAddPointer (&name_list, 0, qual_name);
  }

  if (do_report) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    cip->clickable_item_type = item_type;
    ValNodeAddPointer (&(cip->item_list), choice, data);
    cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (val) + names_len + 1));
    sprintf (cip->description, fmt, val);
    for (vnp = name_list; vnp != NULL; vnp = vnp->next) {
      StringCat (cip->description, vnp->data.ptrvalue);
      if (vnp->next != NULL) {
        StringCat (cip->description, ", ");
      }
    }
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }

  name_list = ValNodeFreeData (name_list);
}


static Boolean LIBCALLBACK IsUnwantedSourceQualifier (ValNodePtr vnp)
{
  if (vnp == NULL) {
    return TRUE;
  } else if (vnp->choice != FieldType_source_qual) {
    return FALSE;
  }
  vnp = vnp->data.ptrvalue;
  if (vnp == NULL) {
    return TRUE;
  } else if (vnp->choice != SourceQualChoice_textqual) {
    return FALSE;
  } else if (vnp->data.intvalue == Source_qual_common_name
             || vnp->data.intvalue == Source_qual_lineage
             || vnp->data.intvalue == Source_qual_division
             || vnp->data.intvalue == Source_qual_old_name
             || vnp->data.intvalue == Source_qual_old_lineage
             || vnp->data.intvalue == Source_qual_gb_acronym
             || vnp->data.intvalue == Source_qual_gb_anamorph
             || vnp->data.intvalue == Source_qual_gb_synonym) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void AdjustSourceQualSampleFieldListForOnCallerTest (ValNodePtr PNTR qual_list, ValNodePtr object_list)
{
  ValNodePtr vnp, field;
  AECRSamplePtr sample;
  if (qual_list == NULL || *qual_list == NULL) {
    return;
  }

  ValNodePurge (qual_list, IsUnwantedSourceQualifier, FieldTypeListFree);

  vnp = ValNodeNew (NULL);
  vnp->choice = SourceQualChoice_location;
  vnp->data.intvalue = 0;
  field = ValNodeNew (NULL);
  field->choice = FieldType_source_qual;
  field->data.ptrvalue = vnp;
  sample = GetAECRSampleFromObjectList (object_list, field);
  if (sample != NULL && sample->num_found > 0) {
    ValNodeLink (qual_list, field);
  } else {
    field = FieldTypeFree (vnp);
  }
  sample = AECRSampleFree (sample);
}


static ValNodePtr SourceQualListForOnCallerTest (SeqEntryPtr sep, ValNodePtr object_list)
{
  ValNodePtr qual_list;

  qual_list = GetSourceQualSampleFieldList (sep);
  AdjustSourceQualSampleFieldListForOnCallerTest (&qual_list, object_list);
  return qual_list;
}


static ClickableItemPtr FindMultipleSourceQuals (ValNodePtr qual, ValNodePtr item_list)
{
  ClickableItemPtr cip = NULL;
  ValNodePtr vnp;
  StringConstraintPtr scp;
  CharPtr             str1, str2, qualname, fmt;
  CharPtr             has_multi_fmt = "%%d sources have multiple %s qualifiers";
  ValNodePtr          has_multi = NULL;
  ValNodePtr          src_choice;

  if (qual == NULL || item_list == NULL) {
    return NULL;
  }
  if (qual->choice == FieldType_source_qual 
      && (src_choice = qual->data.ptrvalue) != NULL 
      && src_choice->choice != SourceQualChoice_textqual) {
    return NULL;
  }

  scp = StringConstraintNew ();
  scp->not_present = TRUE;
  scp->match_location = String_location_equals;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    str1 = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, qual, NULL);
    if (str1 != NULL) {
      scp->match_text = str1;
      str2 = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, qual, scp);
      if (str2 != NULL) {
        ValNodeAddPointer (&has_multi, vnp->choice, vnp->data.ptrvalue);
        str2 = MemFree (str2);
      }
      str1 = MemFree (str1);
    }
  }

  if (has_multi != NULL) {
    qualname = SummarizeFieldType (qual);
    fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (has_multi_fmt) + StringLen (qualname)));
    sprintf (fmt, has_multi_fmt, qualname);
    cip = NewClickableItem (DISC_DUP_SRC_QUAL, fmt, has_multi);
    fmt = MemFree (fmt);
    qualname = MemFree (qualname);
  }
  return cip;
}


static ClickableItemPtr 
SourceQualProblemItem 
(ValNodePtr qual,
 ValNodePtr dup_list,
 ValNodePtr missing_list,
 ValNodePtr src_list,
 ValNodePtr unique_list,
 Uint4      item_type)
{
  ClickableItemPtr cip = NULL, cip_dup, cip_multi;
  CharPtr          some_missing_some_dup = "%s (some missing, some duplicate%s)";
  CharPtr          some_missing = "%s (some missing, all unique%s)";
  CharPtr          some_dup = "%s (all present, some duplicate%s)";
  CharPtr          good = "%s (all present, all unique%s)";
  CharPtr          some_missing_all_same = "%s (some missing, all same%s)";
  CharPtr          all_present_all_same = "%s (all present, all same%s)";
  CharPtr          some_multi = ", some multi";
  CharPtr          unique_fmt = "%%d sources have unique values for %s", unique_desc;
  CharPtr          fmt = NULL;
  CharPtr          qual_name;

  qual_name = SummarizeFieldType (qual);

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  cip->clickable_item_type = item_type;
  cip->item_list = NULL;
  cip->callback_func = NULL;
  cip->datafree_func = NULL;
  cip->callback_data = NULL;
  cip->chosen = 0;
  cip->expanded = FALSE;
  cip->level = 0;
  cip->subcategories = NULL;

  cip_multi = FindMultipleSourceQuals (qual, src_list);

  if (dup_list == NULL && missing_list == NULL) {
    fmt = good;
    cip->item_list = ValNodeCopyPtr (src_list);
  } else if (dup_list != NULL && missing_list != NULL) {
    if (dup_list->next == NULL
        && (cip_dup = dup_list->data.ptrvalue) != NULL
        && ValNodeLen (cip_dup->item_list) == ValNodeLen (src_list) - ValNodeLen (missing_list)) {
      fmt = some_missing_all_same;
    } else {
      fmt = some_missing_some_dup;
    }
    ValNodeLink (&(cip->subcategories), missing_list);
    ValNodeLink (&(cip->subcategories), dup_list);
  } else if (dup_list != NULL) {
    if (dup_list->next == NULL
        && (cip_dup = dup_list->data.ptrvalue) != NULL
        && ValNodeLen (cip_dup->item_list) == ValNodeLen (src_list)) {
      fmt = all_present_all_same;
    } else {      
      fmt = some_dup;
    }
    ValNodeLink (&(cip->subcategories), dup_list);
  } else if (missing_list != NULL) {
    fmt = some_missing;
    cip->subcategories = missing_list;
  }

  if (fmt != NULL) {
    cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (qual_name) + (cip_multi == NULL ? 0 : StringLen (some_multi))));
    sprintf (cip->description, fmt, qual_name, cip_multi == NULL ? "" : some_multi);
  }

  /* note - if we don't use unique_list, we need to free it */
  if (unique_list != NULL && (dup_list != NULL || missing_list != NULL)) {
    unique_list = ValNodeSort (unique_list, SortVnpByDiscrepancyItemText);
    unique_desc = (CharPtr) MemNew (sizeof (Char) * (StringLen (unique_fmt) + StringLen (qual_name)));
    sprintf (unique_desc, unique_fmt, qual_name);
    ValNodeAddPointer (&(cip->subcategories), 0, NewClickableItem (item_type, unique_desc, unique_list));
    unique_desc = MemFree (unique_desc);
  } else {
    unique_list = ValNodeFree (unique_list);
  }


  if (cip_multi != NULL) {
    ValNodeAddPointer (&(cip->subcategories), 0, cip_multi);
  }

  return cip;
}


static void FindRepeatedFieldValues (ValNodePtr PNTR discrepancy_list, ValNodePtr PNTR combo_list, Uint4 item_type)
{
  ValNodePtr repeated = NULL;
  DuplicateQualPtr dq1, dq2;
  ValNodePtr val_dup_list = NULL, item_list, vnp_c;
  ClickableItemPtr cip;

  if (discrepancy_list == NULL || combo_list == NULL || *combo_list == NULL) {
    return;
  }
  /* now look for repeated field values in individual organisms */
  *combo_list = ValNodeSort (*combo_list, SortVnpByDuplicateQualObjectThenValue);

  dq1 = (*combo_list)->data.ptrvalue;
  for (vnp_c = (*combo_list)->next; vnp_c != NULL; vnp_c = vnp_c->next) {
    dq2 = vnp_c->data.ptrvalue;
    if (dq1->choice == dq2->choice && dq1->data == dq2->data && StringCmp (dq1->val, dq2->val) == 0) {
      if (repeated == NULL) {
        ValNodeAddPointer (&repeated, 0, dq1->qual);
      }
      ValNodeAddPointer (&repeated, 0, dq2->qual);
    } else {
      if (repeated != NULL) {
        ReportSameValueMultipleQuals (&val_dup_list, dq1->choice, dq1->data, dq1->val, repeated, item_type);
        repeated = ValNodeFree (repeated);
      }
    }
    dq1 = dq2;
  }
  if (repeated != NULL) {
    ReportSameValueMultipleQuals (&val_dup_list, dq1->choice, dq1->data, dq1->val, repeated, item_type);
    repeated = ValNodeFree (repeated);
  }

  if (val_dup_list != NULL) {
    item_list = ItemListFromSubcategories (val_dup_list);
    RemoveDuplicateItems (&item_list);
    cip = NewClickableItem (item_type, "%d sources have two qualifiers with the same value", item_list);
    cip->subcategories = val_dup_list;
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


static void AddDiscrepanciesForSourceQualComboList (ValNodePtr PNTR discrepancy_list, ValNodePtr PNTR combo_list, ValNodePtr src_list, Uint4 item_type)
{
  ValNodePtr missing_for_qual = NULL, repeated = NULL, unique_list = NULL, dup_qual_list = NULL;
  ValNodePtr subcat = NULL;
  ValNodePtr vnp_c;
  DuplicateQualPtr dq1, dq2;
  ClickableItemPtr cip;
  CharPtr          missing_fmt = "%%d sources are missing %s";
  CharPtr          dup_fmt = "%%d sources have '%s' for %s";
  CharPtr          fmt, qual_name;
  Char             tmp[30];
  ErrSev           msev, lsev;

  if (combo_list == NULL || *combo_list == NULL || src_list == NULL) {
    return;
  }

  /* look for uniqueness across organisms */
  *combo_list = ValNodeSort (*combo_list, SortVnpByDuplicateQualFieldTypeThenValue);
  dq1 = (*combo_list)->data.ptrvalue;
  ValNodeAddPointer (&repeated, dq1->choice, dq1->data);
  for (vnp_c = (*combo_list)->next; vnp_c != NULL; vnp_c = vnp_c->next) {
    dq2 = vnp_c->data.ptrvalue;
    if (CompareDuplicateQual (dq1, dq2) != 0) {
      if (dq1->val == NULL || (StringHasNoText (dq1->val) && !IsNonTextFieldType (dq1->qual))) {
        repeated = ValNodeSort (repeated, SortVnpByDiscrepancyItemText);
        qual_name = SummarizeFieldType (dq1->qual);
        fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (missing_fmt) + StringLen (qual_name)));
        sprintf (fmt, missing_fmt, qual_name);
        ValNodeAddPointer (&missing_for_qual, 0, NewClickableItem (item_type, fmt, repeated));
        qual_name = MemFree (qual_name);
        fmt = MemFree (fmt);
        repeated = NULL;
      } else if (repeated != NULL && repeated->next != NULL) {
        repeated = ValNodeSort (repeated, SortVnpByDiscrepancyItemText);
        qual_name = SummarizeFieldType (dq1->qual);
        fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (dup_fmt) + StringLen (qual_name) + StringLen (dq1->val)));
        sprintf (fmt, dup_fmt, dq1->val, qual_name);
        ValNodeAddPointer (&dup_qual_list, 0, NewClickableItem (item_type, fmt, repeated));
        qual_name = MemFree (qual_name);
        fmt = MemFree (fmt);
        repeated = NULL;
      } else {
        ValNodeLink (&unique_list, repeated);
        repeated = NULL;
      }
      if (CompareFieldTypes (dq1->qual, dq2->qual) != 0) {
        dup_qual_list = ValNodeSort (dup_qual_list, SortVnpByDiscrepancyDescription);
        ValNodeReverse (&dup_qual_list);
        ValNodeAddPointer (&subcat, 0, SourceQualProblemItem (dq1->qual, dup_qual_list, missing_for_qual, src_list, unique_list, item_type));
        dup_qual_list = NULL;
        missing_for_qual = NULL;
        unique_list = NULL;
      }
    }
    ValNodeAddPointer (&repeated, dq2->choice, dq2->data);
    dq1 = dq2;
  }

  if (dq1->val == NULL || (StringHasNoText (dq1->val) && !IsNonTextFieldType (dq1->qual))) {
    repeated = ValNodeSort (repeated, SortVnpByDiscrepancyItemText);
    qual_name = SummarizeFieldType (dq1->qual);
    fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (missing_fmt) + StringLen (qual_name)));
    sprintf (fmt, missing_fmt, qual_name);
    ValNodeAddPointer (&missing_for_qual, 0, NewClickableItem (DISC_MISSING_SRC_QUAL, fmt, repeated));
    qual_name = MemFree (qual_name);
    fmt = MemFree (fmt);
    repeated = NULL;
  } else if (repeated != NULL && repeated->next != NULL) {
    repeated = ValNodeSort (repeated, SortVnpByDiscrepancyItemText);
    qual_name = SummarizeFieldType (dq1->qual);
    fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (dup_fmt) + StringLen (qual_name) + StringLen (dq1->val)));
    sprintf (fmt, dup_fmt, dq1->val, qual_name);
    ValNodeAddPointer (&dup_qual_list, 0, NewClickableItem (DISC_DUP_SRC_QUAL, fmt, repeated));
    qual_name = MemFree (qual_name);
    fmt = MemFree (fmt);
    repeated = NULL;
  } else {
    ValNodeLink (&unique_list, repeated);
    repeated = NULL;
  }
  dup_qual_list = ValNodeSort (dup_qual_list, SortVnpByDiscrepancyDescription);
  ValNodeReverse (&dup_qual_list);
  ValNodeAddPointer (&subcat, 0, SourceQualProblemItem (dq1->qual, dup_qual_list, missing_for_qual, src_list, unique_list, item_type));
  dup_qual_list = NULL;
  missing_for_qual = NULL;
  unique_list = NULL;

  if (subcat != NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->subcategories = subcat;
    cip->clickable_item_type = item_type;
    cip->description = StringSave ("Source Qualifier Report");

    msev = ErrSetMessageLevel (SEV_MAX);
    lsev = ErrSetLogLevel (SEV_MAX);
    if (GetAppParam ("SEQUINCUSTOM", "ONCALLERTOOL", "EXPAND_SRCQUAL_REPORT", NULL, tmp, sizeof (tmp) - 1)
        && StringICmp (tmp, "TRUE") == 0) {
      cip->expanded = TRUE; /* initially source qualifier report should be open */
    }
    ErrSetMessageLevel (msev);
    ErrSetLogLevel (lsev);
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }

}


static void CheckBioSourceQualsEx (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list, Boolean combine_seqentry_reports, Uint4 item_type)
{
  ValNodePtr src_list = NULL, combo_list = NULL, qual_list = NULL;
  ValNodePtr vnp, vnp_q, vnp_s;
  DuplicateQualPtr dq1;
  SeqEntryPtr      sep;

  if (combine_seqentry_reports) {
    src_list = NULL;
    for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
      ValNodeLink (&src_list, GetObjectListForFieldType (FieldType_source_qual, vnp->data.ptrvalue));
    }
    qual_list = GetSourceQualSampleFieldListForSeqEntryList (sep_list);
    AdjustSourceQualSampleFieldListForOnCallerTest (&qual_list, src_list);

    combo_list = NULL;
    /* get all values for all organisms */
    for (vnp_q = qual_list; vnp_q != NULL; vnp_q = vnp_q->next) {
      for (vnp_s = src_list; vnp_s != NULL; vnp_s = vnp_s->next) {
        dq1 = DuplicateQualNew (vnp_s->choice, vnp_s->data.ptrvalue, vnp_q);
        ValNodeAddPointer (&combo_list, 0, dq1);
      }
    }  
    AddDiscrepanciesForSourceQualComboList (discrepancy_list, &combo_list, src_list, item_type);
    /* now look for repeated field values in individual organisms */
    FindRepeatedFieldValues (discrepancy_list, &combo_list, item_type);

    combo_list = DuplicateQualListFree (combo_list);
    src_list = ValNodeFree (src_list);
    qual_list = FieldTypeListFree (qual_list);
  } else {
    for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
      sep = (SeqEntryPtr) vnp->data.ptrvalue;
      src_list = GetObjectListForFieldType (FieldType_source_qual, sep);
      qual_list = SourceQualListForOnCallerTest (sep, src_list); 
      combo_list = NULL;

      /* get all values for all organisms */
      for (vnp_q = qual_list; vnp_q != NULL; vnp_q = vnp_q->next) {
        for (vnp_s = src_list; vnp_s != NULL; vnp_s = vnp_s->next) {
          dq1 = DuplicateQualNew (vnp_s->choice, vnp_s->data.ptrvalue, vnp_q);
          ValNodeAddPointer (&combo_list, 0, dq1);
        }
      }
      AddDiscrepanciesForSourceQualComboList (discrepancy_list, &combo_list, src_list, item_type);
      /* now look for repeated field values in individual organisms */
      FindRepeatedFieldValues (discrepancy_list, &combo_list, item_type);
      combo_list = DuplicateQualListFree (combo_list);
      src_list = ValNodeFree (src_list);
      qual_list = FieldTypeListFree (qual_list);
    }
  }
}


extern void CheckBioSourceQuals (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  CheckBioSourceQualsEx (discrepancy_list, sep_list, FALSE, DISC_SRC_QUAL_PROBLEM);
}


extern void CheckBioSourceQualsAsnDisc (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  CheckBioSourceQualsEx (discrepancy_list, sep_list, TRUE, DISC_SOURCE_QUALS_ASNDISC);
}


typedef Boolean (*BioSourceTestFunc) PROTO ((BioSourcePtr));

typedef struct biosourcetest {
  BioSourceTestFunc func;
  ValNodePtr list;
} BioSourceTestData, PNTR BioSourceTestPtr;


static void BioSourceTestFeatCallback (SeqFeatPtr sfp, Pointer data)
{
  BioSourceTestPtr testdata;
  if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC
      && (testdata = (BioSourceTestPtr) data) != NULL
      && testdata->func != NULL
      && testdata->func (sfp->data.value.ptrvalue)) {
    ValNodeAddPointer (&(testdata->list), OBJ_SEQFEAT, sfp);
  }
}


static void BioSourceTestDescCallback (SeqDescrPtr sdp, Pointer data)
{
  BioSourceTestPtr testdata;
  if (sdp != NULL && sdp->choice == Seq_descr_source
      && (testdata = (BioSourceTestPtr) data) != NULL
      && testdata->func != NULL
      && testdata->func (sdp->data.ptrvalue)) {
    ValNodeAddPointer (&(testdata->list), OBJ_SEQDESC, sdp);
  }
}


static ValNodePtr RunBioSourceTest (SeqEntryPtr sep, BioSourceTestFunc func)
{
  BioSourceTestData data;

  data.func = func;
  data.list = NULL;
  VisitDescriptorsInSep (sep, &data, BioSourceTestDescCallback);
  VisitFeaturesInSep (sep, &data, BioSourceTestFeatCallback);
  return data.list;
}

typedef Boolean (*BioseqTestFunc) PROTO ((BioseqPtr));

typedef struct bioseqtest {
  BioseqTestFunc func;
  ValNodePtr list;
} BioseqTestData, PNTR BioseqTestPtr;


static void BioseqTestBioseqCallback (BioseqPtr bsp, Pointer data)
{
  BioseqTestPtr testdata;
  if (bsp != NULL
      && (testdata = (BioseqTestPtr) data) != NULL
      && testdata->func != NULL
      && testdata->func (bsp)) {
    ValNodeAddPointer (&(testdata->list), OBJ_BIOSEQ, bsp);
  }
}


static ValNodePtr RunBioseqTest (SeqEntryPtr sep, BioseqTestFunc func)
{
  BioseqTestData data;

  data.func = func;
  data.list = NULL;
  VisitBioseqsInSep (sep, &data, BioseqTestBioseqCallback);
  return data.list;
}




static Boolean HasAmplifiedWithSpeciesSpecificPrimerNote (BioSourcePtr biop)
{
  SubSourcePtr ssp;
  OrgModPtr    mod;
  Boolean      rval = FALSE;
  
  if (biop == NULL) {
    return FALSE;
  }
  for (ssp = biop->subtype; ssp != NULL && !rval; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_other 
        && StringCmp (ssp->name, "amplified with species-specific primers") == 0) {
      rval = TRUE;
    }
  }
  if (!rval && biop->org != NULL && biop->org->orgname != NULL) {
    for (mod = biop->org->orgname->mod; mod != NULL && !rval; mod = mod->next) {
      if (mod->subtype == ORGMOD_other
          && StringCmp (mod->subname, "amplified with species-specific primers") == 0) {
        rval = TRUE;
      }
    }
  }
  return rval;
}


static Boolean IsMissingRequiredClone (BioSourcePtr biop)
{
  Boolean needs_clone = FALSE;
  Boolean has_clone = FALSE;
  Boolean has_gel_band_isolate = FALSE;
  SubSourcePtr ssp;
  OrgModPtr    mod;

  if (biop == NULL || HasAmplifiedWithSpeciesSpecificPrimerNote(biop)) {
    return FALSE;
  }
  
  if (biop->org != NULL && StringISearch (biop->org->taxname, "uncultured") != NULL) {
    needs_clone = TRUE;
  }
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_environmental_sample) {
      needs_clone = TRUE;
    } else if (ssp->subtype == SUBSRC_clone) {
      has_clone = TRUE;
    }
  }

  if (needs_clone && !has_clone) {
    /* look for gel band isolate */
    if (biop->org != NULL && biop->org->orgname != NULL) {
      for (mod = biop->org->orgname->mod; mod != NULL && !has_gel_band_isolate; mod = mod->next) {
        if (mod->subtype == ORGMOD_isolate && StringISearch (mod->subname, "gel band") != NULL) {
          has_gel_band_isolate = TRUE;
        }
      }
    }
    if (has_gel_band_isolate) {
      needs_clone = FALSE;
    }
  }

  if (needs_clone && !has_clone) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void FindRequiredClones (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    ValNodeLink (&item_list, RunBioSourceTest (vnp->data.ptrvalue, IsMissingRequiredClone));
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_REQUIRED_CLONE, "%d biosources are missing required clone value", item_list));
  }
}


static Boolean IsMissingRequiredStrain (BioSourcePtr biop)
{
  OrgModPtr mod;

  if (biop == NULL || !IsBacterialBioSource(biop)
    || biop->org == NULL || biop->org->orgname == NULL) {
    return FALSE;
  }
  for (mod = biop->org->orgname->mod; mod != NULL; mod = mod->next) {
    if (mod->subtype == ORGMOD_strain) {
      return FALSE;
    }
  }
  return TRUE;
}


static void FindRequiredStrains (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    ValNodeLink (&item_list, RunBioSourceTest (vnp->data.ptrvalue, IsMissingRequiredStrain));
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_REQUIRED_STRAIN, "%d biosources are missing required strain value", item_list));
  }
}


static Boolean BacterialTaxShouldEndWithStrain (BioSourcePtr biop)
{
  OrgModPtr mod;
  Int4      tax_len, len;

  if (biop == NULL || !IsBacterialBioSource(biop)
      || biop->org == NULL || biop->org->orgname == NULL) {
    return FALSE;
  }
  tax_len = StringLen (biop->org->taxname);
  for (mod = biop->org->orgname->mod; mod != NULL; mod = mod->next) {
    if (mod->subtype == ORGMOD_strain) {
      len = StringLen (mod->subname);
      if (len > tax_len || StringCmp (biop->org->taxname + tax_len - len, mod->subname) != 0) {
        return TRUE;
      }
    }
  }
  return FALSE;
}


static void FindBacterialTaxStrainMismatch (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    ValNodeLink (&item_list, RunBioSourceTest (vnp->data.ptrvalue, BacterialTaxShouldEndWithStrain));
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_BACTERIAL_TAX_STRAIN_MISMATCH, "%d biosources have tax name/strain mismatch", item_list));
  }
}


static Boolean SpNotUncultured (BioSourcePtr biop)
{
  Int4 len;

  if (biop == NULL || biop->org == NULL || (len = StringLen(biop->org->taxname)) < 4
    || StringCmp (biop->org->taxname + len - 4, " sp.") != 0
    || StringNICmp (biop->org->taxname, "uncultured ", 11) == 0) {
    return FALSE;
  } else {
    return TRUE;
  }
}


static void FindSpNotUncultured (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    ValNodeLink (&item_list, RunBioSourceTest (vnp->data.ptrvalue, SpNotUncultured));
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_SP_NOT_UNCULTURED, "%d biosources have taxnames that end with ' sp.' but do not start with 'uncultured'", item_list));
  }
}


static void RetroviridaeDNACallback (BioseqPtr bsp, Pointer data)
{
  SeqMgrDescContext context;
  SeqDescrPtr       sdp;
  BioSourcePtr      biop;

  if (bsp == NULL || bsp->mol != Seq_mol_dna || data == NULL) {
    return;
  }
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp == NULL || (biop = sdp->data.ptrvalue) == NULL
      || biop->genome == GENOME_proviral
      || !HasLineage(biop, "Retroviridae")) {
    return;
  } else {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}


static void CheckRetroviridaeDNA (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, RetroviridaeDNACallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_RETROVIRIDAE_DNA, "%d Retroviridae biosources on DNA sequences are not proviral", item_list));
    item_list = NULL;
  }

}


static void MakeLocationProviral (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr   vnp;
  BioSourcePtr biop;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    biop = GetBioSourceFromObject (vnp->choice, vnp->data.ptrvalue);
    if (biop != NULL) {
      biop->genome = GENOME_proviral;
    }
  }
}


static void CheckForMapChromosomeConflictsCallback (BioseqPtr bsp, Pointer data)
{
  BioSourcePtr biop;
  SeqDescrPtr  sdp;
  SeqMgrDescContext context;
  SubSourcePtr ssp;
  Boolean has_map = FALSE, has_chromosome = FALSE;

  if (!IsEukaryotic (bsp) || ISA_aa (bsp->mol) || data == NULL) {
    return;
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp != NULL && (biop = sdp->data.ptrvalue) != NULL ) {
    for (ssp = biop->subtype; ssp != NULL && (!has_map || !has_chromosome); ssp = ssp->next) {
      if (ssp->subtype == SUBSRC_map) {
        has_map = TRUE;
      } else if (ssp->subtype == SUBSRC_chromosome) {
        has_chromosome = TRUE;
      }
    }
    if (has_map && !has_chromosome) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
    }
  }
}


static void CheckForMapChromosomeConflicts (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL, vnp;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, CheckForMapChromosomeConflictsCallback);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_MAP_CHROMOSOME_CONFLICT, "%d sources on eukaryotic sequences have map but not chromosome", item_list));
  }
}


static void CheckMoltypes (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, object_list, vnp_o, dq_list = NULL;
  DuplicateQualPtr dq1, dq2;
  ValNodePtr field, field2;
  BioseqPtr  bsp;
  ValNodePtr repeated, moltype_list;
  ClickableItemPtr cip;
  CharPtr          moltype_fmt = "%%d sequences have moltype %s";
  CharPtr          fmt;
  Boolean          any_errors = FALSE;

  vnp = ValNodeNew (NULL);
  vnp->choice = MolinfoField_molecule;
  vnp->data.intvalue = 0;

  field = ValNodeNew (NULL);
  field->choice = FieldType_molinfo_field;
  field->data.ptrvalue = vnp;

  vnp = ValNodeNew (NULL);
  vnp->choice = MolinfoField_mol_class;
  vnp->data.intvalue = 0;
  field2 = ValNodeNew (NULL);
  field2->choice = FieldType_molinfo_field;
  field2->data.ptrvalue = vnp;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    object_list = GetObjectListForFieldType (FieldType_molinfo_field, vnp->data.ptrvalue);
    for (vnp_o = object_list; vnp_o != NULL; vnp_o = vnp_o->next) {
      if (vnp_o->choice == OBJ_BIOSEQ && (bsp = vnp_o->data.ptrvalue) != NULL && !ISA_aa (bsp->mol)) {
        dq1 = DuplicateQualNew (vnp_o->choice, vnp_o->data.ptrvalue, field);
        if (StringHasNoText (dq1->val)) {
          dq1->val = MemFree (dq1->val);
          dq1->val = StringSave ("genomic");
        }
        AddFieldValueToDuplicateQual (dq1, field2);
        ValNodeAddPointer (&dq_list, 0, dq1);
      }
    }
    object_list = FreeObjectList (object_list);
    if (dq_list != NULL && dq_list->next != NULL) {
      dq_list = ValNodeSort (dq_list, SortVnpByDuplicateQualFieldTypeThenValue);
      dq1 = dq_list->data.ptrvalue;
      repeated = NULL;
      ValNodeAddPointer (&repeated, dq1->choice, dq1->data);
      moltype_list = NULL;
      for (vnp_o = dq_list->next; vnp_o != NULL; vnp_o = vnp_o->next) {
        dq2 = vnp_o->data.ptrvalue;
        if (StringCmp (dq1->val, dq2->val) != 0) {
          fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (moltype_fmt) + StringLen (dq1->val)));
          sprintf (fmt, moltype_fmt, dq1->val);
          ValNodeAddPointer (&moltype_list, 0, NewClickableItem (DISC_INCONSISTENT_MOLTYPES, fmt, repeated));
          fmt = MemFree (fmt);
          repeated = NULL;
        }
        ValNodeAddPointer (&repeated, dq2->choice, dq2->data);
        dq1 = dq2;
      }
      if (moltype_list == NULL) {
        repeated = ValNodeFree (repeated);
      } else {
        fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (moltype_fmt) + StringLen (dq1->val)));
        sprintf (fmt, moltype_fmt, dq1->val);
        ValNodeAddPointer (&moltype_list, 0, NewClickableItem (DISC_INCONSISTENT_MOLTYPES, fmt, repeated));
        fmt = MemFree (fmt);
        cip = NewClickableItem (DISC_INCONSISTENT_MOLTYPES, "%d sequences have inconsistent moltypes", ItemListFromSubcategories (moltype_list));
        cip->subcategories = moltype_list;
        ValNodeAddPointer (discrepancy_list, 0, cip);
        any_errors = TRUE;
      }
    }
    dq_list = DuplicateQualListFree (dq_list);
  }
  field = FieldTypeFree (field);
  field2 = FieldTypeFree (field2);
  if (!any_errors) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->clickable_item_type = DISC_INCONSISTENT_MOLTYPES;
    cip->description = StringSave ("Moltypes are consistent");
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


static Boolean CitSubMatchExceptDate (CitSubPtr csp1, CitSubPtr csp2)
{
  if (csp1 == NULL && csp2 == NULL) {
    return TRUE;
  } else if (csp1 == NULL || csp2 == NULL) {
    return FALSE;
  } else if (StringCmp (csp1->descr, csp2->descr) != 0
    || csp1->medium != csp2->medium) {
    return FALSE;
  } else if ((csp1->authors == NULL && csp2->authors != NULL)  
             || (csp1->authors != NULL && csp2->authors == NULL)
             || (csp1->authors != NULL && csp2->authors != NULL 
             && !AsnIoMemComp (csp1->authors, csp2->authors, (AsnWriteFunc) AuthListAsnWrite))) {
    return FALSE;
  } else if ((csp1->imp == NULL && csp2->imp != NULL)
             || (csp1->imp != NULL && csp2->imp == NULL)
             || (csp1->imp != NULL && csp2->imp != NULL
                 && !AsnIoMemComp (csp1->imp, csp2->imp, (AsnWriteFunc) ImprintAsnWrite))) {
    return FALSE;
  } else {
    return TRUE;
  }
}


static Boolean SubmitBlockMatchExceptDate (SubmitBlockPtr sb1, SubmitBlockPtr sb2)
{
  if (sb1 == NULL && sb2 == NULL) {
    return TRUE;
  } else if (sb1 == NULL || sb2 == NULL) {
    return FALSE;
  } else if (!AsnIoMemComp (sb1->contact, sb2->contact, (AsnWriteFunc) ContactInfoAsnWrite)) {
    return FALSE;
  } else if (!CitSubMatchExceptDate(sb1->cit, sb2->cit)) {
    return FALSE;
  } else if ((!sb1->hup && sb2->hup) || (sb1->hup && !sb2->hup)) {
    return FALSE;
  } else if (sb1->hup && !DateMatch (sb1->reldate, sb2->reldate, TRUE)) {
    return FALSE;
  } else if (sb1->subtype != sb2->subtype) {
    return FALSE;
  } else if (StringCmp (sb1->tool, sb2->tool) != 0) {
    return FALSE;
  } else if (StringCmp (sb1->user_tag, sb2->user_tag) != 0) {
    return FALSE;
  } else if (StringCmp (sb1->comment, sb2->comment) != 0) {
    return FALSE;
  } else {
    return TRUE;
  }
}


NLM_EXTERN SeqSubmitPtr FindSeqSubmitForSeqEntry (SeqEntryPtr sep)
{
  BioseqPtr bsp;
  BioseqSetPtr bssp;
  SeqSubmitPtr ssp = NULL;

  if (sep == NULL) {
    return NULL;
  }
  if (IS_Bioseq (sep)) {
    bsp = sep->data.ptrvalue;
    if (bsp != NULL && bsp->idx.parentptr != NULL && bsp->idx.parenttype == OBJ_SEQSUB) {
      ssp = bsp->idx.parentptr;
    }
  } else if (IS_Bioseq_set (sep)) {
    bssp = sep->data.ptrvalue;
    if (bssp != NULL && bssp->idx.parentptr != NULL && bssp->idx.parenttype == OBJ_SEQSUB) {
      ssp = bssp->idx.parentptr;
    }
  }
  return ssp;
}


static SubmitBlockPtr FindSubmitBlockForSeqEntry (SeqEntryPtr sep)
{
  SubmitBlockPtr sbp = NULL;
  SeqSubmitPtr ssp = NULL;

  ssp = FindSeqSubmitForSeqEntry (sep);
  if (ssp != NULL) {
    sbp = ssp->sub;
  }
  return sbp;
}


static void CheckSubmitBlockConflicts (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, vnp_m, vnp_s;
  ValNodePtr missing_list = NULL, match_lists = NULL, subcat = NULL, item_list;
  ClickableItemPtr cip;
  SeqEntryPtr sep;
  Boolean     has_any = FALSE, found_match;
  SubmitBlockPtr sbp;

  if (discrepancy_list == NULL || sep_list == NULL || sep_list->next == NULL) {
    return;
  }

  for (vnp = sep_list; vnp != NULL && !has_any; vnp = vnp->next) {
    if (FindSubmitBlockForSeqEntry (vnp->data.ptrvalue) != NULL) {
      has_any = TRUE;
    }
  }

  if (has_any) {
    for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
      sep = vnp->data.ptrvalue;
      if (sep != NULL) {
        sbp = FindSubmitBlockForSeqEntry (sep);
        if (sbp == NULL) {
          if (IS_Bioseq (sep)) {
            ValNodeAddPointer (&missing_list, OBJ_BIOSEQ, sep->data.ptrvalue);
          } else if (IS_Bioseq_set (sep)) {
            ValNodeAddPointer (&missing_list, OBJ_BIOSEQSET, sep->data.ptrvalue);
          }
        } else {
          found_match = FALSE;
          for (vnp_m = match_lists; vnp_m != NULL && !found_match; vnp_m = vnp_m->next) {
            vnp_s = vnp_m->data.ptrvalue;
            if (SubmitBlockMatchExceptDate(sbp, FindSubmitBlockForSeqEntry (vnp_s->data.ptrvalue))) {
              found_match = TRUE;
              ValNodeAddPointer (&vnp_s, 0, sep);
            }
          }
          if (!found_match) {
            vnp_s = ValNodeNew (NULL);
            vnp_s->choice = 0;
            vnp_s->data.ptrvalue = sep;
            ValNodeAddPointer (&match_lists, 0, vnp_s);
          }
        }
      }
    }
    if (missing_list != NULL || (match_lists != NULL && match_lists->next != NULL)) {
      if (missing_list != NULL) {
        ValNodeAddPointer (&subcat, 0, NewClickableItem (DISC_SUBMITBLOCK_CONFLICT, "%d records have no submit-block", missing_list));
      }
      if (match_lists != NULL) {
        for (vnp_m = match_lists; vnp_m != NULL; vnp_m = vnp_m->next) {
          item_list = NULL;
          for (vnp_s = vnp_m->data.ptrvalue; vnp_s != NULL; vnp_s = vnp_s->next) {
            sep = vnp_s->data.ptrvalue;
            if (IS_Bioseq (sep)) {
              ValNodeAddPointer (&item_list, OBJ_BIOSEQ, sep->data.ptrvalue);
            } else if (IS_Bioseq_set (sep)) {
              ValNodeAddPointer (&item_list, OBJ_BIOSEQSET, sep->data.ptrvalue);
            }
          }
          ValNodeAddPointer (&subcat, 0, NewClickableItem (DISC_SUBMITBLOCK_CONFLICT, "%d records have identical submit-blocks", item_list));
        }
      }
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      MemSet (cip, 0, sizeof (ClickableItemData));
      cip->clickable_item_type = DISC_SUBMITBLOCK_CONFLICT;
      cip->description = StringSave ("SubmitBlock Conflicts");
      cip->subcategories = subcat;
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
    for (vnp_m = match_lists; vnp_m != NULL; vnp_m = vnp_m->next) {
      vnp_m->data.ptrvalue = ValNodeFree (vnp_m->data.ptrvalue);
    }
    match_lists = ValNodeFree (match_lists);
  }
}


static PubdescPtr PubdescFromItem (ValNodePtr vnp)
{
  PubdescPtr pdp = NULL;
  SeqDescrPtr sdp;
  SeqFeatPtr sfp;

  if (vnp == NULL) {
    return NULL;
  }
  if (vnp->choice == OBJ_SEQDESC) {
    sdp = (SeqDescrPtr) vnp->data.ptrvalue;
    if (sdp != NULL && sdp->choice == Seq_descr_pub) {
      pdp = sdp->data.ptrvalue;
    }
  } else if (vnp->choice == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_PUB) {
      pdp = sfp->data.value.ptrvalue;
    }
  }
  return pdp;
}


static CitSubPtr CitSubFromPubEquiv (ValNodePtr pub)
{
  CitSubPtr csp = NULL;

  while (pub != NULL && csp == NULL) {
    if (pub->choice == PUB_Sub) {
      csp = pub->data.ptrvalue;
    } else if (pub->choice == PUB_Equiv) {
      csp = CitSubFromPubEquiv (pub->data.ptrvalue);
    }
    pub = pub->next;
  }
  return csp;
}


static CitSubPtr CitSubFromPubdesc (PubdescPtr pdp)
{
  CitSubPtr csp = NULL;

  if (pdp == NULL) {
    return NULL;
  } else {
    csp = CitSubFromPubEquiv (pdp->pub);
  }
  return csp;
}


static AffilPtr AffilFromCitSub (CitSubPtr csp)
{
  AffilPtr affil = NULL;
  if (csp != NULL && csp->authors != NULL ) {
    affil = csp->authors->affil;
  }
  return affil;
}

static int ComparePubAffilForItem (ValNodePtr vnp1, ValNodePtr vnp2)
{
  AffilPtr afp1, afp2;
  CharPtr  str1, str2;
  int rval = 0;

  if (vnp1 == NULL && vnp2 == NULL) {
    rval = 0;
  } else if (vnp1 == NULL) {
    rval = -1;
  } else if (vnp2 == NULL) {
    rval = 1;
  } else {
    afp1 = AffilFromCitSub (CitSubFromPubdesc (PubdescFromItem (vnp1)));
    afp2 = AffilFromCitSub (CitSubFromPubdesc (PubdescFromItem (vnp2)));
    str1 = GetFlatFileAffilString (afp1);
    str2 = GetFlatFileAffilString (afp2);
    rval = StringCmp (str1, str2);
    str1 = MemFree (str1);
    str2 = MemFree (str2);
  }
  return rval;
}

static int LIBCALLBACK SortVnpByPubAffil (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    rval = ComparePubAffilForItem (vnp1, vnp2);
  }

  return rval;
}


static void CollectCitSubPubsFeatCallback (SeqFeatPtr sfp, Pointer data)
{
  ValNode vn;

  if (data != NULL) {
    MemSet (&vn, 0, sizeof (ValNode));
    vn.choice = OBJ_SEQFEAT;
    vn.data.ptrvalue = sfp;
    vn.next = NULL;
    if (CitSubFromPubdesc (PubdescFromItem(&vn)) != NULL) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
    }
  }
}


static void CollectCitSubPubsDescCallback (SeqDescrPtr sdp, Pointer data)
{
  ValNode vn;

  if (data != NULL) {
    MemSet (&vn, 0, sizeof (ValNode));
    vn.choice = OBJ_SEQDESC;
    vn.data.ptrvalue = sdp;
    vn.next = NULL;
    if (CitSubFromPubdesc (PubdescFromItem(&vn)) != NULL) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
    }
  }
}


static void FindMismatchedCitSubAffiliations (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, cit_sub_list = NULL, repeated = NULL, subcat = NULL;
  CharPtr    summ1 = NULL, summ2, fmt, affil_fmt = "%%d CitSubs have affiliation %s";
  ClickableItemPtr cip;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitDescriptorsInSep (vnp->data.ptrvalue, &cit_sub_list, CollectCitSubPubsDescCallback);
    VisitFeaturesInSep (vnp->data.ptrvalue, &cit_sub_list, CollectCitSubPubsFeatCallback);
  }

  cit_sub_list = ValNodeSort (cit_sub_list, SortVnpByPubAffil);
  if (cit_sub_list != NULL && cit_sub_list->next != NULL) {
    summ1 = GetFlatFileAffilString (AffilFromCitSub (CitSubFromPubdesc (PubdescFromItem (cit_sub_list))));
    ValNodeAddPointer (&repeated, cit_sub_list->choice, cit_sub_list->data.ptrvalue);
    for (vnp = cit_sub_list->next; vnp != NULL; vnp = vnp->next) {
      summ2 = GetFlatFileAffilString (AffilFromCitSub (CitSubFromPubdesc (PubdescFromItem (vnp))));
      if (StringCmp (summ1, summ2) != 0) {
        repeated = ValNodeSort (repeated, SortVnpByDiscrepancyItemText);
        if (StringHasNoText (summ1)) {
          ValNodeAddPointer (&subcat, 0, NewClickableItem (DISC_CITSUBAFFIL_CONFLICT, "%d Cit-subs have no affiliation", repeated));
        } else {
          fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (affil_fmt) + StringLen (summ1)));
          sprintf (fmt, affil_fmt, summ1);
          ValNodeAddPointer (&subcat, 0, NewClickableItem (DISC_CITSUBAFFIL_CONFLICT, fmt, repeated));
          fmt = MemFree (fmt);
        }
        repeated = NULL;
      }
      ValNodeAddPointer (&repeated, vnp->choice, vnp->data.ptrvalue);
      summ1 = MemFree (summ1);
      summ1 = summ2;
    }
    repeated = ValNodeSort (repeated, SortVnpByDiscrepancyItemText);
    if (StringHasNoText (summ1)) {
      ValNodeAddPointer (&subcat, 0, NewClickableItem (DISC_CITSUBAFFIL_CONFLICT, "%d Cit-subs have no affiliation", repeated));
    } else {
      fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (affil_fmt) + StringLen (summ1)));
      sprintf (fmt, affil_fmt, summ1);
      ValNodeAddPointer (&subcat, 0, NewClickableItem (DISC_CITSUBAFFIL_CONFLICT, fmt, repeated));
      fmt = MemFree (fmt);
    }
    repeated = NULL;
  }

  if (subcat == NULL) {
    if (cit_sub_list == NULL) {
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      MemSet (cip, 0, sizeof (ClickableItemData));
      cip->clickable_item_type = DISC_CITSUBAFFIL_CONFLICT;
      cip->description = StringSave ("No citsubs were found!");
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
  } else if (subcat->next == NULL && !StringHasNoText (summ1)) {
    /* Make no report if all values match 
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->clickable_item_type = DISC_CITSUBAFFIL_CONFLICT;
    cip->description = StringSave ("All citsub affiliations match");
    ValNodeAddPointer (discrepancy_list, 0, cip); */
    subcat = FreeClickableList (subcat);
  } else {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->clickable_item_type = DISC_CITSUBAFFIL_CONFLICT;
    cip->description = StringSave ("Citsub affiliation conflicts found");
    cip->subcategories = subcat;
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
  summ1 = MemFree (summ1);
  cit_sub_list = ValNodeFree (cit_sub_list);
}


typedef struct haplotypesequence {
  CharPtr haplotype;
  CharPtr taxname;
  BioseqPtr bsp;
} HaplotypeSequenceData, PNTR HaplotypeSequencePtr;


static HaplotypeSequencePtr HaplotypeSequenceNew (CharPtr haplotype, CharPtr taxname, BioseqPtr bsp)
{
  HaplotypeSequencePtr h;

  h = (HaplotypeSequencePtr) MemNew (sizeof (HaplotypeSequenceData));
  h->haplotype = haplotype;
  h->taxname = taxname;
  h->bsp = bsp;
  return h;
}



static int CompareSubSequences (BioseqPtr bsp1, Int4 pos1, BioseqPtr bsp2, Int4 pos2, Int4 cmp_len, Boolean allow_Ndiff)
{
  int  rval = 0;
  Int4 buf_len = 49;
  Char buf1[50];
  Char buf2[50];
  CharPtr cp1, cp2;
  Int2 ctr;

  if (bsp1 == NULL && bsp2 == NULL) {
    return 0;
  } else if (bsp1 == NULL) {
    return -1;
  } else if (bsp2 == NULL) {
    return 1;
  }

  while (pos1 < bsp1->length && pos2 < bsp2->length && rval == 0 && cmp_len > 0) {
    ctr = SeqPortStreamInt (bsp1, pos1, MIN (MIN(pos1 + buf_len - 1, bsp1->length - 1), pos1 + cmp_len - 1), Seq_strand_plus,
                        STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL,
                        (Pointer) buf1, NULL);
    buf1[ctr] = 0;
    ctr = SeqPortStreamInt (bsp2, pos2, MIN (MIN(pos2 + buf_len - 1, bsp2->length - 1), pos2 + cmp_len - 1), Seq_strand_plus,
                        STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL,
                        (Pointer) buf2, NULL);
    buf2[ctr] = 0;

    cp1 = buf1;
    cp2 = buf2;
    while (*cp1 != 0 && *cp2 != 0 && rval == 0) {
      if (allow_Ndiff && (*cp1 == 'N' || *cp2 == 'N')) {
        /* ok - can continue */
      } else if (*cp1 == *cp2) {
        /* identical, can continue */
      } else if (*cp1 < *cp2) {
        rval = -1;
      } else {
        rval = 1;
      }
      ++cp1;
      ++cp2;
    }
    if (*cp1 == 0 && *cp2 != 0) {
      rval = -1;
    } else if (*cp1 != 0 && *cp2 == 0) {
      rval = 1;
    }
    pos1 += buf_len;
    pos2 += buf_len;
    cmp_len -= buf_len;
  }
  return rval;
}


static int CompareSequences (BioseqPtr bsp1, BioseqPtr bsp2, Boolean allow_Ndiff)
{
  int       rval = 0;
  
  if (bsp1 != NULL && bsp2 != NULL) {
    if (bsp1->length < bsp2->length) {
      rval = -1;
    } else if (bsp1->length > bsp2->length) {
      rval = 1;
    } else {
      rval = CompareSubSequences (bsp1, 0, bsp2, 0, bsp1->length, allow_Ndiff);
    }
  }
  return rval;
}


static Boolean SequencesHaveOverlap (BioseqPtr bsp1, BioseqPtr bsp2, Boolean allow_Ndiff)
{
  Int4 pct_overlap_required = 50;
  Int4 overlap_len, min_overlap_len;
  int  rval = -1;
  Int4 offset = 0;

  if (bsp1->length > bsp2->length) {
    min_overlap_len = (pct_overlap_required * bsp2->length) / 100;
  } else {
    min_overlap_len = (pct_overlap_required * bsp1->length) / 100;
  }
  while (rval != 0 && offset < bsp1->length - min_overlap_len) {
    overlap_len = MIN (bsp2->length - offset, bsp1->length - offset);
    rval = CompareSubSequences(bsp1, offset, bsp2, 0, overlap_len, allow_Ndiff);
    offset++;
  }
  offset = 0;
  while (rval != 0 && offset < bsp2->length - min_overlap_len) {
    overlap_len = MIN (bsp2->length - offset, bsp1->length - offset);
    rval = CompareSubSequences(bsp2, offset, bsp1, 0, overlap_len, allow_Ndiff);
    offset++;
  }
  if (rval == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean DoSequencesMatchForHaplotype (BioseqPtr bsp1, BioseqPtr bsp2, Boolean allow_Ndiff)
{
  Int4 diff;
  int  rval = -1;


  if (bsp1 == NULL && bsp2 == NULL) {
    return TRUE;
  } else if (bsp1 == NULL || bsp2 == NULL) {
    return FALSE;
  }

  if (bsp1->length == bsp2->length) {
    rval = CompareSubSequences (bsp1, 0, bsp2, 0, bsp1->length, allow_Ndiff);
  } else if (bsp1->length > bsp2->length) {
    diff = bsp1->length - bsp2->length;
    while (rval != 0 && diff >= 0) {
      rval = CompareSubSequences (bsp1, diff, bsp2, 0, bsp2->length, allow_Ndiff);
      diff--;
    }
  } else {
    diff = bsp2->length - bsp1->length;
    while (rval != 0 && diff >= 0) {
      rval = CompareSubSequences (bsp1, 0, bsp2, diff, bsp1->length, allow_Ndiff);
      diff--;
    }
  }
  if (rval != 0 && SequencesHaveOverlap(bsp1, bsp2, allow_Ndiff)) {
    rval = 0;
  }

  if (rval == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static int CompareHaplotypeThenSequence (HaplotypeSequencePtr a, HaplotypeSequencePtr b, Boolean allowNDiff)
{
  int   rval = 0;
  
  if (a != NULL && b != NULL) {
    rval = StringCmp (a->taxname, b->taxname);
    if (rval == 0) {
      rval = StringCmp (a->haplotype, b->haplotype);
      if (rval == 0) {
        rval = CompareSequences (a->bsp, b->bsp, allowNDiff);
      }
    }
  }
  return rval;
}


static int LIBCALLBACK SortVnpByHaplotypeThenSequence (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    rval = CompareHaplotypeThenSequence (vnp1->data.ptrvalue, vnp2->data.ptrvalue, FALSE);
  }

  return rval;
}


static int LIBCALLBACK SortVnpByHaplotypeThenSequenceAllowNDiff (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    rval = CompareHaplotypeThenSequence (vnp1->data.ptrvalue, vnp2->data.ptrvalue, TRUE);
  }

  return rval;
}


static int CompareSequenceThenHaplotype (HaplotypeSequencePtr a, HaplotypeSequencePtr b, Boolean allowNDiff)
{
  int   rval = 0;
  
  if (a != NULL && b != NULL) {
    rval = CompareSequences (a->bsp, b->bsp, allowNDiff);
    if (rval == 0) {
      rval = StringCmp (a->taxname, b->taxname);
      if (rval == 0) {
        rval = StringCmp (a->haplotype, b->haplotype);
      }
    }
  }
  return rval;
}


static int LIBCALLBACK SortVnpBySequenceThenHaplotype (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    rval = CompareSequenceThenHaplotype (vnp1->data.ptrvalue, vnp2->data.ptrvalue, FALSE);
  }

  return rval;
}


static int LIBCALLBACK SortVnpBySequenceThenHaplotypeAllowNDiff (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    rval = CompareSequenceThenHaplotype (vnp1->data.ptrvalue, vnp2->data.ptrvalue, TRUE);
  }

  return rval;
}



static void HaplotypeCollectionCallback (BioseqPtr bsp, Pointer data)
{
  SeqDescrPtr sdp;
  SeqMgrDescContext context;
  BioSourcePtr      biop;
  SubSourcePtr      ssp;
  CharPtr           taxname = NULL;

  if (bsp != NULL && data != NULL && !ISA_aa (bsp->mol)) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
    if (sdp != NULL) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop != NULL) {
        ssp = biop->subtype;
        while (ssp != NULL && ssp->subtype != SUBSRC_haplotype) {
          ssp = ssp->next;
        }
        if (ssp != NULL) {
          if (biop->org != NULL) {
            taxname = biop->org->taxname;
          }
          ValNodeAddPointer ((ValNodePtr PNTR) data, 0, HaplotypeSequenceNew (ssp->name, taxname, bsp));
        }
      }
    }
  }
}


static void 
ReportOneHaplotypeSequenceMismatch 
(ValNodePtr PNTR discrepancy_list, 
 CharPtr taxname, 
 CharPtr haplotype, 
 ValNodePtr mismatch_list,
 Boolean    allowNDiff)
{
  CharPtr  seq_mismatch_fmt = "%%d sequences have organism %s haplotype %s but the sequences do not match%s";
  CharPtr  allow_N_fmt = " (allowing N to match any)";
  CharPtr  strict_N_fmt = " (strict match)";
  Int4     fmt_len;
  CharPtr  fmt;

  fmt_len = StringLen (seq_mismatch_fmt) + StringLen (taxname) + StringLen (haplotype);
  if (allowNDiff) {
    fmt_len += StringLen (allow_N_fmt);
  } else {
    fmt_len += StringLen (strict_N_fmt);
  }

  fmt = (CharPtr) MemNew (sizeof (Char) * fmt_len);
  sprintf (fmt, seq_mismatch_fmt, taxname, haplotype, allowNDiff ? allow_N_fmt : strict_N_fmt);

  ValNodeAddPointer (discrepancy_list, 0, 
                     NewClickableItem (DISC_HAPLOTYPE_MISMATCH, fmt, mismatch_list));
  fmt = MemFree (fmt);
}


static void
ReportOneSequenceMatchHaplotypeMismatch
(ValNodePtr PNTR discrepancy_list, 
 ValNodePtr mismatch_list,
 Boolean    allowNDiff)
{
  CharPtr  hap_mismatch_strict_fmt = "%d sequences are identical (strict match) but have different haplotypes";
  CharPtr  hap_mismatch_allowN_fmt = "%d sequences are identical (allowing N to match any) but have different haplotypes";
  ValNodePtr src_qual, field_list, extended_item_list;

  src_qual = ValNodeNew (NULL);
  src_qual->choice = SourceQualChoice_textqual;
  src_qual->data.intvalue = Source_qual_haplotype;
  field_list = ValNodeNew (NULL);
  field_list->choice = FieldType_source_qual;
  field_list->data.ptrvalue = src_qual;
    
  extended_item_list = MakeObjectListWithFields (mismatch_list, field_list);
  mismatch_list = ValNodeFree (mismatch_list);
  field_list = FieldTypeListFree (field_list);

  ValNodeAddPointer (discrepancy_list, 0,
                     NewClickableItem (DISC_HAPLOTYPE_MISMATCH, 
                                       allowNDiff ? hap_mismatch_allowN_fmt : hap_mismatch_strict_fmt, extended_item_list));
}


static ValNodePtr ReportHaplotypeSequenceMismatchForList (ValNodePtr PNTR haplotype_sequence_list, Boolean allow_NDiff)
{
  ValNodePtr vnp_h;
  ValNodePtr same_list, subcat = NULL;
  HaplotypeSequencePtr h1, h2;
  Boolean  have_mismatch;

  if (haplotype_sequence_list == NULL || *haplotype_sequence_list == NULL) {
    return subcat;
  }

  /* first, look for same taxname, same haplotype, different sequence */
  *haplotype_sequence_list = ValNodeSort (*haplotype_sequence_list, allow_NDiff ? SortVnpByHaplotypeThenSequence : SortVnpByHaplotypeThenSequence);
  have_mismatch = FALSE;
  same_list = NULL;
  h1 = (*haplotype_sequence_list)->data.ptrvalue;
  for (vnp_h = (*haplotype_sequence_list)->next; vnp_h != NULL; vnp_h = vnp_h->next) {
    h2 = vnp_h->data.ptrvalue;
    if (StringCmp (h1->taxname, h2->taxname) == 0 && StringCmp (h1->haplotype, h2->haplotype) == 0) {
      if (same_list == NULL) {
        have_mismatch = FALSE;
        ValNodeAddPointer (&same_list, OBJ_BIOSEQ, h1->bsp);
      }
      ValNodeAddPointer (&same_list, OBJ_BIOSEQ, h2->bsp);
      if (!DoSequencesMatchForHaplotype (h1->bsp, h2->bsp, allow_NDiff)) {
        have_mismatch = TRUE;
      }
    } else {
      if (same_list != NULL) {
        /* add discrepancy report */
        if (have_mismatch) {
          ReportOneHaplotypeSequenceMismatch (&subcat, h1->taxname, h1->haplotype, same_list, allow_NDiff);
        } else {
          same_list = ValNodeFree (same_list);
        }
      }
      same_list = NULL;
    }
    h1 = h2;
  }
  if (same_list != NULL) {
    if (have_mismatch) {
      /* add discrepancy report */
      ReportOneHaplotypeSequenceMismatch (&subcat, h1->taxname, h1->haplotype, same_list, allow_NDiff);
    } else {
      same_list = ValNodeFree (same_list);
    }
  }

  /* now look for sequence that match but have different haplotypes */
  *haplotype_sequence_list = ValNodeSort (*haplotype_sequence_list, allow_NDiff ? SortVnpBySequenceThenHaplotypeAllowNDiff : SortVnpBySequenceThenHaplotype);
  same_list = NULL;
  have_mismatch = FALSE;
  h1 = (*haplotype_sequence_list)->data.ptrvalue;
  for (vnp_h = (*haplotype_sequence_list)->next; vnp_h != NULL; vnp_h = vnp_h->next) {
    h2 = vnp_h->data.ptrvalue;
    if (CompareSequences (h1->bsp, h2->bsp, allow_NDiff) == 0) {
      if (same_list == NULL) {
        ValNodeAddPointer (&same_list, OBJ_BIOSEQ, h1->bsp);
        have_mismatch = FALSE;
      }
      ValNodeAddPointer (&same_list, OBJ_BIOSEQ, h2->bsp);
      if (StringCmp (h1->haplotype, h2->haplotype) != 0) {
        have_mismatch = TRUE;
      }
    } else {
      if (same_list != NULL) {
        if (have_mismatch) {
          ReportOneSequenceMatchHaplotypeMismatch (&subcat, same_list, allow_NDiff);
        } else {
          same_list = ValNodeFree (same_list);
        }
        same_list = NULL;
      }
    }
    h1 = h2;
  }
  if (same_list != NULL) {
    if (have_mismatch) {
      ReportOneSequenceMatchHaplotypeMismatch (&subcat, same_list, allow_NDiff);
    } else {
      same_list = ValNodeFree (same_list);
    }
    same_list = NULL;
  }

  subcat = ValNodeSort (subcat, SortVnpByDiscrepancyDescription);
  ValNodeReverse (&subcat);

  return subcat;
}


static void ReportHaplotypeSequenceMismatch (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp;
  ValNodePtr haplotype_sequence_list, strict_match = NULL, nonstrict_match = NULL, subcat = NULL;
  CharPtr  mismatch_loose_fmt = "There are %d haplotype problems (loose match, allowing Ns to differ)";
  CharPtr  mismatch_strict_fmt = "There are %d haplotype problems (strict match)";
  ClickableItemPtr cip_main, cip_loose = NULL, cip_strict = NULL;

  if (discrepancy_list == NULL || sep_list == NULL) {
    return;
  }

  /* Note - analysis should be performed separately for each SeqEntry, rather than for the list as a whole */
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    haplotype_sequence_list = NULL;
    VisitBioseqsInSep (vnp->data.ptrvalue, &haplotype_sequence_list, HaplotypeCollectionCallback);

    nonstrict_match = ReportHaplotypeSequenceMismatchForList (&haplotype_sequence_list, TRUE);
    if (nonstrict_match != NULL) {
      cip_loose = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      cip_loose->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (mismatch_loose_fmt) + 15));
      sprintf (cip_loose->description, mismatch_loose_fmt, ValNodeLen (nonstrict_match));
      cip_loose->subcategories = nonstrict_match;
      ValNodeAddPointer (&subcat, 0, cip_loose);
    }

    strict_match = ReportHaplotypeSequenceMismatchForList (&haplotype_sequence_list, FALSE);
    if (strict_match != NULL) {
      cip_strict = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      cip_strict->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (mismatch_strict_fmt) + 15));
      sprintf (cip_strict->description, mismatch_strict_fmt, ValNodeLen (strict_match));
      cip_strict->subcategories = strict_match;
      ValNodeAddPointer (&subcat, 0, cip_strict);
    }
  }

  if (subcat != NULL) {
    cip_main = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip_main, 0, sizeof (ClickableItemData));
    cip_main->clickable_item_type = DISC_HAPLOTYPE_MISMATCH;
    cip_main->description = StringSave ("Haplotype Problem Report");
    cip_main->subcategories = subcat;
    ValNodeAddPointer (discrepancy_list, 0, cip_main);
  }
}


static Boolean IsGenomicDNASequence (BioseqPtr bsp)
{
  SeqMgrDescContext dcontext;
  SeqDescrPtr       sdp;
  Boolean           rval = FALSE;
  MolInfoPtr        mip;

  if (bsp == NULL) {
    rval = FALSE;
  } else if (bsp->mol != Seq_mol_dna) {
    rval = FALSE;
  } else {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
    if (sdp != NULL) {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL && mip->biomol == MOLECULE_TYPE_GENOMIC) {
        rval = TRUE;
      }
    }
  }
  return rval;
}


static void ReportFeatureMoltypeMismatchCallback (BioseqPtr bsp, Pointer data)
{
  SeqFeatPtr  sfp;
  SeqMgrFeatContext fcontext;

  if (bsp == NULL || data == NULL) {
    return;
  }

  if (!IsGenomicDNASequence(bsp)) {
    sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_rRNA, &fcontext);
    if (sfp != NULL) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQ, bsp);
    } else if ((sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_otherRNA, &fcontext)) != NULL) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQ, bsp);
    }
  }
}


static void ReportFeatureMoltypeMismatch (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp;
  ValNodePtr item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, ReportFeatureMoltypeMismatchCallback);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_FEATURE_MOLTYPE_MISMATCH, "%d sequences have rRNA or misc_RNA features but are not genomic DNA", item_list));
  }
}


/* change the sequences on which rRNA features are located to genomic DNA */
static void ChangeMoltypeToGenomicDNA (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr vnp;
  BioseqPtr  bsp;
  SeqDescrPtr sdp;
  SeqMgrDescContext dcontext;
  MolInfoPtr        mip;
  Char              id_txt[255];

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_BIOSEQ) {
      bsp = vnp->data.ptrvalue;
      if (bsp != NULL) {
        bsp->mol = Seq_mol_dna;
        sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
        if (sdp == NULL) {
          sdp = CreateNewDescriptorOnBioseq (bsp, Seq_descr_molinfo);
        }
        mip = (MolInfoPtr) sdp->data.ptrvalue;
        if (mip == NULL) {
          mip = MolInfoNew ();
          sdp->data.ptrvalue = mip;
        }
        mip->biomol = MOLECULE_TYPE_GENOMIC;
        bsp->strand = 0;
        if (lip != NULL && lip->fp != NULL) {
          SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
          fprintf (lip->fp, "Changed biomol for %s\n", id_txt);
          lip->data_in_log = TRUE;
        }
      }
    }
  }
}


const CharPtr kmRNAVariant = ", transcript variant ";
const CharPtr kCDSVariant = ", isoform ";

NLM_EXTERN Boolean ProductsMatchForRefSeq (CharPtr cds_str, CharPtr mrna_str)
{
  CharPtr join_mrna, join_cds;
  Int4    len;

  if (StringHasNoText (cds_str) || StringHasNoText (mrna_str)) {
    return FALSE;
  }

  join_mrna = StringStr (mrna_str, kmRNAVariant);
  if (join_mrna == NULL) {
    return FALSE;
  }
  join_cds = StringStr (cds_str, kCDSVariant);
  if (join_cds == NULL) {
    return FALSE;
  }
  len = join_mrna - mrna_str;
  if (len != join_cds - cds_str) {
    return FALSE;
  } else if (StringNCmp (cds_str, mrna_str, len) != 0) {
    return FALSE;
  }
  cds_str = join_cds + StringLen (kCDSVariant);
  mrna_str = join_mrna + StringLen (kmRNAVariant);
  if (StringCmp (cds_str, mrna_str) != 0) {
    return FALSE;
  } else {
    return TRUE;
  }
}


NLM_EXTERN SeqFeatPtr GetmRNAforCDS (SeqFeatPtr cds)
{
  SeqFeatPtr      mrna = NULL;
  SeqFeatXrefPtr  xref;
  SeqMgrFeatContext mcontext;

  /* first, check for mRNA identified by feature xref */
  for (xref = cds->xref; xref != NULL && mrna == NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      mrna = SeqMgrGetFeatureByFeatID (cds->idx.entityID, NULL, NULL, xref, NULL);
      if (mrna != NULL && mrna->idx.subtype != FEATDEF_mRNA) {
        mrna = NULL;
      }
    }
  }
  
  /* try by location if not by xref */
  if (mrna == NULL) {
    mrna = SeqMgrGetLocationSupersetmRNA (cds->location, &mcontext);
    if (mrna == NULL) {
      mrna = SeqMgrGetOverlappingmRNA (cds->location, &mcontext);
    }
  }
  return mrna;
}

typedef struct underlyingfeat {
  SeqFeatPtr orig_feat;
  ValNodePtr matching_features;
} UnderlyingFeatData, PNTR UnderlyingFeatPtr;

static Boolean LIBCALLBACK FindUnderlyingCDS (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)

{
  UnderlyingFeatPtr  uf;

  if (sfp == NULL || context == NULL) return TRUE;
  uf = context->userdata;
  if (uf == NULL) return TRUE;

  if (TestFeatOverlap(uf->orig_feat, sfp, CHECK_INTERVALS) >= 0) {
    ValNodeAddPointer (&(uf->matching_features), OBJ_SEQFEAT, sfp);
  }

  return TRUE;
}


NLM_EXTERN SeqFeatPtr GetCDSformRNA (SeqFeatPtr mrna)
{
  SeqFeatPtr      cds = NULL;
  SeqFeatXrefPtr  xref;
  Int2 count;
  UnderlyingFeatData uf;

  /* first, check for cds identified by feature xref */
  for (xref = mrna->xref; xref != NULL && cds == NULL; xref = xref->next) {
    if (xref->id.choice != 0) {
      cds = SeqMgrGetFeatureByFeatID (mrna->idx.entityID, NULL, NULL, xref, NULL);
      if (cds != NULL && cds->idx.subtype != FEATDEF_CDS) {
        cds = NULL;
      }
    }
  }
  
  /* try by location if not by xref */
  if (cds == NULL) {
    MemSet (&uf, 0, sizeof (UnderlyingFeatData));
    uf.orig_feat = mrna;
    count = SeqMgrGetAllOverlappingFeatures (mrna->location, FEATDEF_CDS, NULL, 0, 
                                             SIMPLE_OVERLAP, &uf, FindUnderlyingCDS);
    if (uf.matching_features != NULL) {
      cds = uf.matching_features->data.ptrvalue;
      uf.matching_features = ValNodeFree (uf.matching_features);
    }
  }
  return cds;
}


static void ReportCDSWithoutmRNACallback (BioseqPtr bsp, Pointer data)
{
  SeqMgrFeatContext fcontext;
  SeqMgrDescContext dcontext;
  SeqFeatPtr        sfp, mRNA;
  SeqDescrPtr       sdp;
  MolInfoPtr        mip;
  CharPtr           feat_product, mrna_product;
  ValNode           field;
  FeatureFieldPtr   ff;
  BioSourcePtr      biop;

  if (bsp == NULL || bsp->mol != Seq_mol_dna || data == NULL) {
    return;
  }
  
  if (!IsEukaryotic (bsp)) {
    return;
  }
  biop = GetBiopForBsp(bsp);
  if (biop != NULL && IsLocationOrganelle(biop->genome)) {
    return;
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp == NULL || sdp->data.ptrvalue == NULL) {
    return;
  }
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip->biomol != MOLECULE_TYPE_GENOMIC) {
    return;
  }

  ff = FeatureFieldNew ();
  ff->type = Macro_feature_type_any;
  ValNodeAddInt (&(ff->field), FeatQualChoice_legal_qual, Feat_qual_legal_product);
  field.choice = FieldType_feature_field;
  field.data.ptrvalue = ff;
  field.next = NULL;

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_CDS, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_CDS, &fcontext)) {
    if (IsPseudo (sfp)) {
      continue;
    }

    mRNA = GetmRNAforCDS(sfp);

    if (mRNA == NULL) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
    } else {
      feat_product = GetFieldValueForObject (OBJ_SEQFEAT, sfp, &field, NULL);
      mrna_product = GetFieldValueForObject (OBJ_SEQFEAT, mRNA, &field, NULL);
      if (StringCmp (feat_product, mrna_product) != 0 && !ProductsMatchForRefSeq(feat_product, mrna_product)) {
        ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
      }
      feat_product = MemFree (feat_product);
      mrna_product = MemFree (mrna_product);
    }
  }

  ff = FeatureFieldFree (ff);
}


static void ReportCDSWithoutmRNA (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp;
  ValNodePtr item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, ReportCDSWithoutmRNACallback);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_CDS_WITHOUT_MRNA, "%d coding regions do not have an mRNA", item_list));
  }
}


static ProtRefPtr FindBestProtRef (Uint2 entityID, SeqFeatPtr cds)

{
  SeqFeatPtr bestprot;
  
  if (cds == NULL) return NULL;
  bestprot = FindBestProtein (entityID, cds->product);
  if (bestprot != NULL) {
    return bestprot->data.value.ptrvalue;
  } else {
    return NULL;
  }
}


static SeqIntPtr CombineSeqInt (SeqIntPtr sint1, SeqIntPtr sint2)
{
  SeqIntPtr sint_combined = NULL;

  if (sint1 == NULL || sint2 == NULL) {
    return NULL;
  }
  sint_combined = SeqIntNew ();
  sint_combined->id = SeqIdDup (sint1->id);
  sint_combined->strand = sint1->strand;
  sint_combined->from = sint1->from;
  sint_combined->if_from = AsnIoMemCopy (sint1->if_from, (AsnReadFunc)IntFuzzAsnRead, (AsnWriteFunc) IntFuzzAsnWrite);
  sint_combined->to = sint2->to;
  sint_combined->if_to = AsnIoMemCopy (sint2->if_to, (AsnReadFunc)IntFuzzAsnRead, (AsnWriteFunc) IntFuzzAsnWrite);

  return sint_combined;
}


static SeqLocPtr CombineLocations (SeqLocPtr slp1, SeqLocPtr utr, BioseqPtr bsp)
{
  SeqLocPtr slp_combined = NULL, slp, tmp, slp_prev;
  Uint1     strand1, strand2;
  Int4         start1, start2, stop1, stop2;
  SeqIntPtr    sint_combined;

  if (slp1 == NULL || utr == NULL) {
    return NULL;
  }

  strand1 = SeqLocStrand (slp1);
  strand2 = SeqLocStrand (utr);
  if (strand1 == Seq_strand_minus && strand2 != Seq_strand_minus) {
    return NULL;
  } else if (strand1 != Seq_strand_minus && strand2 == Seq_strand_minus) {
    return NULL;
  }

  start1 = SeqLocStart (slp1);
  stop1 = SeqLocStop (slp1);
  start2 = SeqLocStart (utr);
  stop2 = SeqLocStop (utr);
  if (strand1 == Seq_strand_minus) {
    /* allow overlap for 3' UTR */
    if (stop2 >= start1 - 1 && stop2 < start1 + 3) {
      slp_combined = SeqLocMergeEx (bsp, slp1, utr, FALSE, FALSE, FALSE, FALSE);
      if (slp_combined != NULL && slp_combined->choice == SEQLOC_MIX) {
        slp = slp_combined->data.ptrvalue;
        while (slp != NULL && SeqLocStart (slp) != start1) {
          slp = slp->next;
        }
        /* if we have adjacent intervals at the point where the main loc ends, combine them */
        if (slp != NULL && slp->next != NULL && slp->choice == SEQLOC_INT && slp->next->choice == SEQLOC_INT) {
          sint_combined = CombineSeqInt (slp->next->data.ptrvalue, slp->data.ptrvalue);
          if (sint_combined != NULL) {
            tmp = slp->next;
            slp->next = slp->next->next;
            tmp->next = NULL;
            tmp = SeqLocFree (tmp);
            slp->data.ptrvalue = SeqIntFree (slp->data.ptrvalue);
            slp->data.ptrvalue = sint_combined;
          }
        }
      }
      /* no overlap for 5' UTR */
    } else if (start2 == stop1 + 1) {
      slp_combined = SeqLocMergeEx (bsp, utr, slp1, FALSE, FALSE, FALSE, FALSE);
      if (slp_combined != NULL && slp_combined->choice == SEQLOC_MIX) {
        slp = slp_combined->data.ptrvalue;
        slp_prev = NULL;
        while (slp != NULL && SeqLocStop (slp) != stop1) {
          slp_prev = slp;
          slp = slp->next;
        }
        if (slp != NULL && slp_prev != NULL && slp_prev->choice == SEQLOC_INT && slp->choice == SEQLOC_INT) {
          sint_combined = CombineSeqInt (slp->data.ptrvalue, slp_prev->data.ptrvalue);
          if (sint_combined != NULL) {
            slp_prev->next = slp->next;
            slp->next = NULL;
            slp = SeqLocFree (slp);
            slp_prev->data.ptrvalue = SeqIntFree (slp_prev->data.ptrvalue);
            slp_prev->data.ptrvalue = sint_combined;
          }
        }
      }
    }
  } else {
    /* allow overlap for 3' UTR */
    if (start2 > stop1 - 3 && start2 <= stop1 + 1) {
      slp_combined = SeqLocMergeEx (bsp, slp1, utr, FALSE, FALSE, FALSE, FALSE);
      if (slp_combined != NULL && slp_combined->choice == SEQLOC_MIX) {
        slp = slp_combined->data.ptrvalue;
        while (slp != NULL && SeqLocStop (slp) != stop1) {
          slp = slp->next;
        }
        /* if we have adjacent intervals at the point where the main loc ends, combine them */
        if (slp != NULL && slp->next != NULL && slp->choice == SEQLOC_INT && slp->next->choice == SEQLOC_INT) {
          sint_combined = CombineSeqInt (slp->data.ptrvalue, slp->next->data.ptrvalue);
          if (sint_combined != NULL) {
            tmp = slp->next;
            slp->next = slp->next->next;
            tmp->next = NULL;
            tmp = SeqLocFree (tmp);
            slp->data.ptrvalue = SeqIntFree (slp->data.ptrvalue);
            slp->data.ptrvalue = sint_combined;
          }
        }
      }
      /* no overlap for 5' UTR */
    } else if (stop2 == start1 - 1) {
      slp_combined = SeqLocMergeEx (bsp, utr, slp1, FALSE, FALSE, FALSE, FALSE);
      if (slp_combined != NULL && slp_combined->choice == SEQLOC_MIX) {
        slp = slp_combined->data.ptrvalue;
        slp_prev = NULL;
        while (slp != NULL && SeqLocStart (slp) != start1) {
          slp_prev = slp;
          slp = slp->next;
        }
        if (slp != NULL && slp_prev != NULL && slp_prev->choice == SEQLOC_INT && slp->choice == SEQLOC_INT) {
          sint_combined = CombineSeqInt (slp_prev->data.ptrvalue, slp->data.ptrvalue);
          if (sint_combined != NULL) {
            slp_prev->next = slp->next;
            slp->next = NULL;
            slp = SeqLocFree (slp);
            slp_prev->data.ptrvalue = SeqIntFree (slp_prev->data.ptrvalue);
            slp_prev->data.ptrvalue = sint_combined;
          }
        }
      }
    }
  }

  return slp_combined;
}


NLM_EXTERN SeqLocPtr GetmRNALocationFromCDSLocation (SeqLocPtr slp, Uint2 entityID)
{
  BioseqPtr bsp;
  SeqLocPtr slp_mrna = NULL, tmp;
  Uint1     strand;
  SeqFeatPtr utr5, utr3;
  SeqMgrFeatContext context;
  Int4              pos5, pos3;
  Boolean           found;
  Boolean           partial5 = TRUE, partial3 = TRUE;

  bsp = GetBioseqGivenSeqLoc (slp, entityID);
  strand = SeqLocStrand (slp);
  if (strand == Seq_strand_minus) {
    pos5 = SeqLocStop (slp);
    pos3 = SeqLocStart (slp);
  } else {
    pos5 = SeqLocStart (slp);
    pos3 = SeqLocStop (slp);
  }

  slp_mrna = AsnIoMemCopy ((Pointer) slp,
                                      (AsnReadFunc) SeqLocAsnRead,
                                      (AsnWriteFunc) SeqLocAsnWrite);  

  utr5 = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_5UTR, &context);
  if (strand == Seq_strand_minus) {
    while (utr5 != NULL && context.left < pos5 + 1) {
      utr5 = SeqMgrGetNextFeature (bsp, utr5, 0, FEATDEF_5UTR, &context);
    }
    if (context.left == pos5 + 1 && utr5 != NULL) {
      tmp = CombineLocations (slp_mrna, utr5->location, bsp);
      if (tmp != NULL) {
        slp_mrna = SeqLocFree (slp_mrna);
        slp_mrna = tmp;
        CheckSeqLocForPartial (utr5->location, &partial5, NULL);
      }
    } 
  } else {
    found = FALSE;
    while (utr5 != NULL && !found && context.left < pos5) {
      if (context.right == pos5 - 1) {
        tmp = CombineLocations (slp_mrna, utr5->location, bsp);
        if (tmp != NULL) {
          slp_mrna = SeqLocFree (slp_mrna);
          slp_mrna = tmp;
          CheckSeqLocForPartial (utr5->location, &partial5, NULL);
        }
        found = TRUE;
      } else {
        utr5 = SeqMgrGetNextFeature (bsp, utr5, 0, FEATDEF_5UTR, &context);
      }
    }
  }

  utr3 = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_3UTR, &context);
  if (strand == Seq_strand_minus) {
    found = FALSE;
    while (utr3 != NULL && !found && context.left < pos3 + 2) {
      if (context.right >= pos3 - 1 && context.right < pos3 + 2) {
        tmp = CombineLocations (slp_mrna, utr3->location, bsp);
        if (tmp != NULL) {
          slp_mrna = SeqLocFree (slp_mrna);
          slp_mrna = tmp;
          CheckSeqLocForPartial (utr3->location, NULL, &partial3);
        }
        found = TRUE;
      } else {
        utr3 = SeqMgrGetNextFeature (bsp, utr3, 0, FEATDEF_3UTR, &context);
      }
    }
  } else {
    found = FALSE;
    while (utr3 != NULL && !found && context.left < pos3 + 2) {
      if (context.left >= pos3 - 2 && context.left < pos3 + 2) {
        found = TRUE;
        tmp = CombineLocations (slp_mrna, utr3->location, bsp);
        if (tmp != NULL) {
          slp_mrna = SeqLocFree (slp_mrna);
          slp_mrna = tmp;
          CheckSeqLocForPartial (utr3->location, NULL, &partial3);
        }
      } else {
        utr3 = SeqMgrGetNextFeature (bsp, utr3, 0, FEATDEF_3UTR, &context);
      }
    }
  }

  SetSeqLocPartial (slp_mrna, partial5, partial3);
  return slp_mrna;
}


NLM_EXTERN void AddmRNAForCDS (SeqFeatPtr sfp)
{
  RnaRefPtr rrp;
  ProtRefPtr prp;
  ValNodePtr name;
  CharPtr    mRNAname = NULL;
  SeqFeatPtr rna, gene;
  SeqEntryPtr sep;
  Boolean     partial5, partial3;
  BioseqPtr   bsp;
  SeqMgrFeatContext fcontext;

  rrp = RnaRefNew ();
  if (rrp != NULL) {
    rrp->type = 2;
    prp = FindBestProtRef (sfp->idx.entityID, sfp);
    if (prp != NULL) {
      name = prp->name;
      if (name != NULL && !StringHasNoText (name->data.ptrvalue)) {
        mRNAname = StringSave (name->data.ptrvalue);
      } else if (!StringHasNoText (prp->desc)) {
        mRNAname = StringSave (prp->desc);
      }
    }
    if (mRNAname!= NULL) {
      rrp->ext.choice = 1;
      rrp->ext.value.ptrvalue = mRNAname;
    }
    rna = SeqFeatNew ();
    if (rna != NULL) {
      rna->data.choice = SEQFEAT_RNA;
      rna->data.value.ptrvalue = (Pointer) rrp;
      rna->location = GetmRNALocationFromCDSLocation (sfp->location, sfp->idx.entityID);
      CheckSeqLocForPartial (rna->location, &partial5, &partial3);
      rna->partial = (rna->partial || partial5 || partial3);
      bsp = GetBioseqGivenSeqLoc (rna->location, sfp->idx.entityID);
      if (bsp != NULL) {
        sep = SeqMgrGetSeqEntryForData (bsp);
        if (sep != NULL) {
          CreateNewFeature (sep, NULL, SEQFEAT_RNA, rna);
        } else {
          rna->next = sfp->next;
          sfp->next = rna;
        }
      } else {
        rna->next = sfp->next;
        sfp->next = rna;
      }
      /* if gene location matches mRNA exactly, make it partial on both ends */
      gene = SeqMgrGetOverlappingGene (rna->location, &fcontext);
      if (gene != NULL && SeqLocAinB (rna->location, gene->location) == 0) {
        SetSeqLocPartial (gene->location, TRUE, TRUE);
        gene->partial = TRUE;
      }
    }
  }
}


static void AddMissingmRNA (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr vnp;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT) {
      AddmRNAForCDS (vnp->data.ptrvalue);
    }
  }
}


static void ReportmRNAOnNonGenomicEukaryoticSequencesCallback (BioseqPtr bsp, Pointer data)
{
  SeqMgrFeatContext fcontext;
  SeqMgrDescContext dcontext;
  SeqFeatPtr        sfp;
  SeqDescrPtr       sdp;
  MolInfoPtr        mip;
  BioSourcePtr      biop;

  if (bsp == NULL || bsp->mol != Seq_mol_dna || data == NULL) {
    return;
  }
  
  if (!IsEukaryotic (bsp)) {
    return;
  }
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp == NULL || sdp->data.ptrvalue == NULL) {
    return;
  }
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip->biomol != MOLECULE_TYPE_GENOMIC) {
    return;
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL || sdp->data.ptrvalue == NULL) {
    return;
  }
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop->genome == GENOME_macronuclear || biop->genome == GENOME_unknown || biop->genome == GENOME_genomic) {
    return;
  }

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_RNA, FEATDEF_mRNA, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_RNA, FEATDEF_mRNA, &fcontext)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
}


static void ReportmRNAOnNonGenomicEukaryoticSequences (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp;
  ValNodePtr item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, ReportmRNAOnNonGenomicEukaryoticSequencesCallback);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_mRNA_ON_WRONG_SEQUENCE_TYPE, "%d mRNAs are located on eukaryotic sequences that do not have genomic or plasmid sources", item_list));
  }
}




/* exon and intron locations */
/* if an intron starts or stops between the start of one exon and the end of the next exon, it should
 * abut both exons.
 */

static ValNodePtr CompareIntronExonList (ValNodePtr exon_list, ValNodePtr intron_list)
{
  SeqFeatPtr        exon, next_exon, intron;
  ValNodePtr        vnp_e, vnp_i;
  Int4              exon_start, exon_stop, intron_start, intron_stop, next_exon_start, next_exon_stop;
  ValNodePtr        problem_list = NULL;

  if (exon_list != NULL && intron_list != NULL) {
    exon = exon_list->data.ptrvalue;
    exon_start = SeqLocStart (exon->location);
    exon_stop = SeqLocStop (exon->location);
    intron = intron_list->data.ptrvalue;
    intron_start = SeqLocStart (intron->location);
    intron_stop = SeqLocStop (intron->location);

    if (intron_start < exon_start) {
      if (intron_stop != exon_start - 1) {
        ValNodeAddPointer (&problem_list, OBJ_SEQFEAT, intron_list->data.ptrvalue);
        ValNodeAddPointer (&problem_list, OBJ_SEQFEAT, exon_list->data.ptrvalue);
      }
      vnp_i = intron_list->next;
      if (vnp_i != NULL) {
        intron = vnp_i->data.ptrvalue;
        intron_start = SeqLocStart (intron->location);
        intron_stop = SeqLocStop (intron->location);
      }
    } else {
      vnp_i = intron_list;
    }

    for (vnp_e = exon_list->next; vnp_e != NULL && vnp_i != NULL; vnp_e = vnp_e->next) {
      next_exon = vnp_e->data.ptrvalue;
      next_exon_start = SeqLocStart (next_exon->location);
      next_exon_stop = SeqLocStop (next_exon->location);
      while (vnp_i != NULL && intron_start < next_exon_start) {
        if (intron_start != exon_stop + 1 || intron_stop != next_exon_start - 1) {
          if (intron_start != exon_stop + 1) {
            ValNodeAddPointer (&problem_list, OBJ_SEQFEAT, exon);
          }
          ValNodeAddPointer (&problem_list, OBJ_SEQFEAT, intron);
          if (intron_stop != next_exon_start - 1) {
            ValNodeAddPointer (&problem_list, OBJ_SEQFEAT, next_exon);
          }
        }
        vnp_i = vnp_i->next;
        if (vnp_i != NULL) {
          intron = vnp_i->data.ptrvalue;
          intron_start = SeqLocStart (intron->location);
          intron_stop = SeqLocStop (intron->location);
        }
      }
      exon = next_exon;
      exon_start = next_exon_start;
      exon_stop = next_exon_stop;
    }
    if (vnp_i != NULL) {
      if (intron_start != exon_stop + 1) {
        ValNodeAddPointer (&problem_list, OBJ_SEQFEAT, exon);
        ValNodeAddPointer (&problem_list, OBJ_SEQFEAT, intron);
      }
    }

    RemoveDuplicateItems (&problem_list);
  }
  return problem_list;
}


static ValNodePtr GetFeatureListForGene (BioseqPtr bsp, SeqFeatPtr gene, Uint1 featdef)
{
  SeqFeatPtr        feat, feat_gene;
  SeqMgrFeatContext fcontext, gcontext;
  ValNodePtr        feat_list = NULL;
  GeneRefPtr        grp;

  for (feat = SeqMgrGetNextFeature (bsp, NULL, 0, featdef, &fcontext);
       feat != NULL;
       feat = SeqMgrGetNextFeature (bsp, feat, 0, featdef, &fcontext)) {
    if (gene == NULL) {
      /* collect all */
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, feat);
    } else if ((grp = SeqMgrGetGeneXref (feat)) == NULL) {
      /* find by overlap */
      feat_gene = SeqMgrGetOverlappingGene(feat->location, &gcontext);
      if (feat_gene == gene) {
        ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, feat);
      }
    } else if (!SeqMgrGeneIsSuppressed(grp) && GeneRefMatch(grp, gene->data.value.ptrvalue)) {
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, feat);
    }
  }
  return feat_list;
}


static void CheckIntronAndExonLocationsOnBioseq (BioseqPtr bsp, Pointer data)
{
  SeqFeatPtr        gene;
  SeqMgrFeatContext gcontext;
  ValNodePtr        exon_list = NULL, intron_list = NULL;
  ValNodePtr        problem_list = NULL;
  Char              id[255];
  CharPtr           fmt, problems_fmt = "%%d introns and exons have location conflicts on %s";

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL) {
    return;
  }

  gene = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &gcontext);
  if (gene == NULL) {
    /* no genes - just do all exons and introns present */
    exon_list = GetFeatureListForGene (bsp, NULL, FEATDEF_exon);
    intron_list = GetFeatureListForGene (bsp, NULL, FEATDEF_intron);
    problem_list = CompareIntronExonList (exon_list, intron_list);
    exon_list = ValNodeFree (exon_list);
    intron_list = ValNodeFree (intron_list);
  } else {
    while (gene != NULL) {
      if (StringICmp (gene->except_text, "trans-splicing") != 0) {
        exon_list = GetFeatureListForGene (bsp, gene, FEATDEF_exon);
        intron_list = GetFeatureListForGene (bsp, gene, FEATDEF_intron);
        ValNodeLink (&problem_list, CompareIntronExonList (exon_list, intron_list));
        exon_list = ValNodeFree (exon_list);
        intron_list = ValNodeFree (intron_list);
      }
      gene = SeqMgrGetNextFeature (bsp, gene, SEQFEAT_GENE, 0, &gcontext);
    }
  }

  if (problem_list != NULL) {
    SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id, PRINTID_REPORT, sizeof (id) - 1);
    fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (problems_fmt) + StringLen (id)));
    sprintf (fmt, problems_fmt, id);
    ValNodeAddPointer ((ValNodePtr PNTR) data, 0, NewClickableItem (DISC_EXON_INTRON_CONFLICT, fmt, problem_list));
    fmt = MemFree (fmt);
  }

}


static void CheckIntronAndExonLocations (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp;
  ValNodePtr disc_list = NULL;
  ClickableItemPtr cip;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &disc_list, CheckIntronAndExonLocationsOnBioseq);
  }
  if (disc_list != NULL) {
    cip = NewClickableItem (DISC_EXON_INTRON_CONFLICT, "%d introns and exons are incorrectly positioned", ItemListFromSubcategories (disc_list));
    cip->subcategories = disc_list;
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


typedef Boolean (*ValNodeExtractTestFunc) PROTO ((ValNodePtr, Pointer));

static ValNodePtr ValNodeExtractListByFunction (ValNodePtr PNTR list, ValNodeExtractTestFunc func, Pointer data)
{
  ValNodePtr vnp, vnp_prev = NULL, vnp_next, new_list = NULL;

  if (list == NULL || *list == NULL || func == NULL) {
    return NULL;
  }

  for (vnp = *list; vnp != NULL; vnp = vnp_next) {
    vnp_next = vnp->next;
    if (func (vnp, data)) {
      if (vnp_prev == NULL) {
        *list = vnp->next;
      } else {
        vnp_prev->next = vnp->next;
      }
      vnp->next = NULL;
      ValNodeLink (&new_list, vnp);
    } else {
      vnp_prev = vnp;
    }
  }
  return new_list;
}
  

typedef struct featurecountdata {
  BioseqPtr bsp;
  CharPtr   seq_id_txt;
  Int4  featdef;
  Int4  num_feats;
} FeatureCountData, PNTR FeatureCountPtr;


static FeatureCountPtr FeatureCountNew (BioseqPtr bsp, Int4 featdef)
{
  FeatureCountPtr f;

  f = (FeatureCountPtr) MemNew (sizeof (FeatureCountData));
  f->bsp = bsp;
  f->seq_id_txt = NULL;
  f->featdef = featdef;
  f->num_feats = 0;
  return f;
}


static FeatureCountPtr FeatureCountFree (FeatureCountPtr f)
{
  if (f != NULL) {
    f->seq_id_txt = MemFree (f->seq_id_txt);
    f = MemFree (f);
  }
  return f;
}


static ValNodePtr FeatureCountListFree (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = FeatureCountFree (list->data.ptrvalue);
    list = ValNodeFree (list);
    list = list_next;
  }
  return list;
}


static void SaveFeatureCountSequenceIds (ValNodePtr list, CharPtr filename)
{
  FeatureCountPtr f;
  ValNode vn;

  MemSet (&vn, 0, sizeof (ValNode));
  vn.choice = OBJ_BIOSEQ;
  vn.next = NULL;
  while (list != NULL) {
    f = (FeatureCountPtr) list->data.ptrvalue;
    if (f != NULL && f->bsp != NULL) {
      vn.data.ptrvalue = f->bsp;
      f->seq_id_txt = GetDiscrepancyItemTextEx (&vn, filename);
      f->bsp = NULL;
    }
    list = list->next;
  }
}


static ValNodePtr GetSequenceIdListFromFeatureCountList (ValNodePtr feat_count_list)
{
  ValNodePtr seq_list = NULL, vnp;
  FeatureCountPtr f;

  for (vnp = feat_count_list; vnp != NULL; vnp = vnp->next) {
    f = (FeatureCountPtr) vnp->data.ptrvalue;
    if (f != NULL && !StringHasNoText (f->seq_id_txt)) {
      ValNodeAddPointer (&seq_list, 0, StringSave (f->seq_id_txt));
    }
  }
  return seq_list;
}


static ValNodePtr GetFeatureTypesFromFeatureCounts (ValNodePtr feat_count_list)
{
  ValNodePtr vnp;
  ValNodePtr feat_type_list = NULL;
  FeatureCountPtr f;
  Int4 sort_countdown = 100;

  for (vnp = feat_count_list; vnp != NULL; vnp = vnp->next) {
    f = (FeatureCountPtr) vnp->data.ptrvalue;
    if (f != NULL) {
      ValNodeAddInt (&feat_type_list, 0, f->featdef);
    }
    sort_countdown--;
    if (sort_countdown == 0) {
      feat_type_list = ValNodeSort (feat_type_list, SortByIntvalue);
      ValNodeUnique (&feat_type_list, SortByIntvalue, ValNodeFree);
      sort_countdown = 100;
    }
  }

  feat_type_list = ValNodeSort (feat_type_list, SortByIntvalue);
  ValNodeUnique (&feat_type_list, SortByIntvalue, ValNodeFree);
  return feat_type_list;
}


static Int4 GetNumFeaturesInList (ValNodePtr list)
{
  FeatureCountPtr f;
  ValNodePtr vnp;
  Int4 num = 0;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    f = (FeatureCountPtr) vnp->data.ptrvalue;
    if (f != NULL) {
      num += f->num_feats;
    }
  }
  return num;
}


static Boolean FeatureCountHasFeatdef (ValNodePtr vnp, Pointer data)
{
  Int4 featdef;
  FeatureCountPtr f;

  if (vnp == NULL || data == NULL) {
    return FALSE;
  }

  featdef = *((Int4Ptr)data);
  f = (FeatureCountPtr) vnp->data.ptrvalue;
  if (f != NULL && f->featdef == featdef) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean FeatureCountHasNumFeats (ValNodePtr vnp, Pointer data)
{
  Int4 num_feats;
  FeatureCountPtr f;

  if (vnp == NULL || data == NULL) {
    return FALSE;
  }

  num_feats = *((Int4Ptr)data);
  f = (FeatureCountPtr) vnp->data.ptrvalue;
  if (f != NULL && f->num_feats == num_feats) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void InsertMissingFeatureCountsWithSeqIdTxt (ValNodePtr PNTR feat_count_list)
{
  ValNodePtr seq_list, feat_list, feat_seq_list, tmp_list, vnp, new_list = NULL;
  ValNodePtr v1, v2;
  Int4       featdef;
  FeatureCountPtr f;

  if (feat_count_list == NULL || *feat_count_list == NULL) {
    return;
  }

  seq_list = GetSequenceIdListFromFeatureCountList (*feat_count_list);
  seq_list = ValNodeSort (seq_list, SortVnpByString);
  feat_list = GetFeatureTypesFromFeatureCounts (*feat_count_list);

  for (vnp = feat_list; vnp != NULL; vnp = vnp->next) {
    featdef = vnp->data.intvalue;
    tmp_list = ValNodeExtractListByFunction (feat_count_list, FeatureCountHasFeatdef, &featdef);
    feat_seq_list = GetSequenceIdListFromFeatureCountList (tmp_list);
    feat_seq_list = ValNodeSort (feat_seq_list, SortVnpByString);
    v1 = seq_list;
    v2 = feat_seq_list;
    while (v1 != NULL) {
      if (v2 != NULL && StringCmp (v2->data.ptrvalue, v1->data.ptrvalue) == 0) {
        v2 = v2->next;
      } else {
        f = FeatureCountNew (NULL, featdef);
        f->seq_id_txt = StringSave (v1->data.ptrvalue);
        ValNodeAddPointer (&tmp_list, 0, f);
      }
      v1 = v1->next;
    }
    ValNodeLink (&new_list, tmp_list);
    feat_seq_list = ValNodeFreeData (feat_seq_list);
    tmp_list = NULL;
  }
  seq_list = ValNodeFreeData (seq_list);
  feat_list = ValNodeFree (feat_list);

  *feat_count_list = new_list;
}


static int CompareFeatureCounts (FeatureCountPtr f1, FeatureCountPtr f2)
{
  int rval = 0;

  if (f1 != NULL && f2 != NULL) {
    if (f1->featdef < f2->featdef) {
      rval = -1;
    } else if (f1->featdef > f2->featdef) {
      rval = 1;
    } else if (f1->num_feats < f2->num_feats) {
      rval = -1;
    } else if (f1->num_feats > f2->num_feats) {
      rval = 1;
    }
  }
  return rval;
}


static int LIBCALLBACK SortVnpFeatureCount (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);

    if (vnp1->data.ptrvalue != NULL && vnp2->data.ptrvalue != NULL) {
      rval = CompareFeatureCounts (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    }
  }

  return rval;
}


static void CountFeaturesOnSequenceCallback (BioseqPtr bsp, Pointer data)
{
  ValNodePtr featdef_list = NULL, vnp;
  SeqFeatPtr sfp;
  SeqMgrFeatContext context;
  FeatureCountPtr f;

  if (bsp == NULL || data == NULL) {
    return;
  }

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context)) {
    if (sfp->idx.subtype != FEATDEF_PROT) {
      ValNodeAddInt (&featdef_list, 0, sfp->idx.subtype);
    }
  }
  featdef_list = ValNodeSort (featdef_list, SortByIntvalue);

  if (featdef_list != NULL) {
    f = FeatureCountNew (bsp, featdef_list->data.intvalue);
    f->num_feats = 1;
    ValNodeAddPointer ((ValNodePtr PNTR) data, 0, f);
    vnp = featdef_list->next;
    while (vnp != NULL) {
      if (vnp->data.intvalue == f->featdef) {
        f->num_feats++;
      } else {
        f = FeatureCountNew (bsp, vnp->data.intvalue);
        f->num_feats = 1;
        ValNodeAddPointer ((ValNodePtr PNTR) data, 0, f);
      }
      vnp = vnp->next;
    }
    featdef_list = ValNodeFree (featdef_list);
  }
}


typedef struct missingcountsdata {
  ValNodePtr feat_type_list;
  ValNodePtr feat_count_list;
} MissingCountsData, PNTR MissingCountsPtr;

static void AddMissingFeatureCountsCallback (BioseqPtr bsp, Pointer data)
{
  MissingCountsPtr  m;
  ValNodePtr        vnp;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext context;
  Uint1             seqfeattype;

  if (bsp == NULL || data == NULL) {
    return;
  }

  m = (MissingCountsPtr) data;

  for (vnp = m->feat_type_list; vnp != NULL; vnp = vnp->next) {
    seqfeattype = FindFeatFromFeatDefType (vnp->data.intvalue);
    if ((seqfeattype == SEQFEAT_PROT && ISA_aa (bsp->mol))
        || (seqfeattype != SEQFEAT_PROT && !ISA_aa (bsp->mol))) {
      sfp = SeqMgrGetNextFeature (bsp, NULL, 0, vnp->data.intvalue, &context);
      if (sfp == NULL) {
        ValNodeAddPointer (&(m->feat_count_list), 0, FeatureCountNew (bsp, vnp->data.intvalue));
      }
    }
  }
}


static ClickableItemPtr AddFeatureCountReport (Int4 featdef, Int4 num, ValNodePtr bsp_list)
{
  ClickableItemPtr cip;
  CharPtr          fmt = "%d bioseqs have %d %s features";
  Int4             feature_type;
  CharPtr          feature_name;

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  cip->clickable_item_type = DISC_FEATURE_COUNT;
  cip->item_list = bsp_list;

  feature_type = GetFeatureTypeFromFeatdef (featdef);
  feature_name = GetFeatureNameFromFeatureType (feature_type);

  cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (feature_name) + 30));
  sprintf (cip->description, fmt, ValNodeLen (bsp_list), num, feature_name);

  cip->callback_func = NULL;
  cip->datafree_func = NULL;
  cip->callback_data = NULL;
  cip->subcategories = NULL;
  cip->expanded = FALSE;
  cip->level = 0;

  return cip;
}


static ClickableItemPtr AddFeatureTypeSummary (Int4 featdef, ValNodePtr disc_list, Int4 total)
{
  ClickableItemPtr cip;
  Int4             feature_type, len;
  CharPtr          feature_name;
  CharPtr          fmt = "%s: %d present%s";
  CharPtr          inconsistent = " (inconsistent)";

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  cip->clickable_item_type = DISC_FEATURE_COUNT;
  cip->item_list = NULL;

  feature_type = GetFeatureTypeFromFeatdef (featdef);
  feature_name = GetFeatureNameFromFeatureType (feature_type);

  len = StringLen (fmt) + StringLen (feature_name) + 15;
  if (disc_list != NULL) {
    len += StringLen (inconsistent);
    disc_list = ValNodeSort (disc_list, SortVnpByDiscrepancyDescription);
    ValNodeReverse (&disc_list);
  }
  cip->description = (CharPtr) MemNew (sizeof (Char) * len);
  sprintf (cip->description, fmt, feature_name, total, disc_list == NULL ? "" : inconsistent);

  cip->callback_func = NULL;
  cip->datafree_func = NULL;
  cip->callback_data = NULL;
  cip->subcategories = disc_list;
  cip->expanded = FALSE;
  cip->level = 0;

  return cip;
}


static void CountFeaturesOnSequences (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp;
  ValNodePtr feat_count_list = NULL;
  ValNodePtr bsp_list = NULL, num_list = NULL, type_list = NULL, ok_list = NULL;
  ClickableItemPtr cip;
  Int4             current_featdef, current_num, feat_total;
  MissingCountsData m;
  FeatureCountPtr f;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &feat_count_list, CountFeaturesOnSequenceCallback);
  }

  m.feat_type_list = GetFeatureTypesFromFeatureCounts (feat_count_list);
  m.feat_count_list = NULL;
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &m, AddMissingFeatureCountsCallback);
  }
  ValNodeLink (&feat_count_list, m.feat_count_list);
  m.feat_type_list = ValNodeFree (m.feat_type_list);

  feat_count_list = ValNodeSort (feat_count_list, SortVnpFeatureCount);
  if (feat_count_list != NULL) {
    f = feat_count_list->data.ptrvalue;
    current_featdef = f->featdef;
    current_num = f->num_feats;
    feat_total = current_num;
    ValNodeAddPointer (&bsp_list, OBJ_BIOSEQ, f->bsp);
    vnp = feat_count_list->next;
    while (vnp != NULL) {
      f = vnp->data.ptrvalue;
      if (f->featdef == current_featdef) {
        if (f->num_feats != current_num) {
          bsp_list = ValNodeSort (bsp_list, SortVnpByDiscrepancyItemText);
          ValNodeAddPointer (&num_list, 0, AddFeatureCountReport (current_featdef, current_num, bsp_list));
          bsp_list = NULL;
          current_featdef = f->featdef;
          current_num = f->num_feats;
        }
        feat_total += current_num;
        ValNodeAddPointer (&bsp_list, OBJ_BIOSEQ, f->bsp);
      } else {
        bsp_list = ValNodeSort (bsp_list, SortVnpByDiscrepancyItemText);
        ValNodeAddPointer (&num_list, 0, AddFeatureCountReport (current_featdef, current_num, bsp_list));
        bsp_list = NULL;
        if (num_list->next == NULL) {
          cip = AddFeatureTypeSummary (current_featdef, NULL, feat_total);
          cip->subcategories = num_list;
          ValNodeAddPointer (&ok_list, 0, cip);
        } else {
          ValNodeAddPointer (&type_list, 0, AddFeatureTypeSummary (current_featdef, num_list, feat_total));
        }
        num_list = NULL;
        current_featdef = f->featdef;
        current_num = f->num_feats;
        feat_total = current_num;
        ValNodeAddPointer (&bsp_list, OBJ_BIOSEQ, f->bsp);
      }
      vnp = vnp->next;
    }

    bsp_list = ValNodeSort (bsp_list, SortVnpByDiscrepancyItemText);
    ValNodeAddPointer (&num_list, 0, AddFeatureCountReport (current_featdef, current_num, bsp_list));
    bsp_list = NULL;
    if (num_list->next == NULL) {
      cip = AddFeatureTypeSummary (current_featdef, NULL, feat_total);
      cip->subcategories = num_list;
      ValNodeAddPointer (&ok_list, 0, cip);
    } else {
      ValNodeAddPointer (&type_list, 0, AddFeatureTypeSummary (current_featdef, num_list, feat_total));
    }

    ValNodeLink (&type_list, ok_list);

    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    cip->clickable_item_type = DISC_FEATURE_COUNT;
    cip->subcategories = type_list;
    cip->item_list = NULL;
    cip->description = StringSave ("Feature Counts");

    cip->callback_func = NULL;
    cip->datafree_func = NULL;
    cip->callback_data = NULL;
    cip->expanded = TRUE;
    cip->level = 0;

    ValNodeAddPointer (discrepancy_list, 0, cip);

    /* free list now that we are done with it */
    feat_count_list = FeatureCountListFree (feat_count_list);
  }
}


typedef struct taxnameconflict {
  CharPtr qual;
  CharPtr taxname;
  Uint1   obj_type;
  Pointer obj_data;
} TaxNameConflictData, PNTR TaxNameConflictPtr;


static TaxNameConflictPtr TaxNameConflictNew (CharPtr qual, CharPtr taxname, Uint1 obj_type, Pointer obj_data)
{
  TaxNameConflictPtr h;

  h = (TaxNameConflictPtr) MemNew (sizeof (TaxNameConflictData));
  h->qual = qual;
  h->taxname = taxname;
  h->obj_type = obj_type;
  h->obj_data = obj_data;
  return h;
}


static int LIBCALLBACK SortTaxNameConflict (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  TaxNameConflictPtr s1, s2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    s1 = vnp1->data.ptrvalue;
    s2 = vnp2->data.ptrvalue;
    if (s1 != NULL && s2 != NULL) {
      rval = StringICmp (s1->qual, s2->qual);
      if (rval == 0) {
        rval = StringICmp (s1->taxname, s2->taxname);
      }
    }
  }

  return rval;
}


static void CollectTaxnameConflictDiscrepancies 
(ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list, 
 VisitFeaturesFunc feat_callback, VisitDescriptorsFunc desc_callback,
 CharPtr qual_name, Uint4 item_type)
{
  ValNodePtr vnp;
  ValNodePtr st_list = NULL, spec_list = NULL, disc_list = NULL;
  TaxNameConflictPtr st1, st2 = NULL;
  Boolean               have_mismatch;
  CharPtr               spec_fmt = "%%d biosources have %s %s but do not have the same taxnames";
  CharPtr               top_fmt = "%%d BioSources have %s/taxname conflicts";
  CharPtr               fmt;
  ClickableItemPtr      cip;

  if (sep_list == NULL || discrepancy_list == NULL) {
    return;
  }

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitDescriptorsInSep (vnp->data.ptrvalue, &st_list, desc_callback);
    VisitFeaturesInSep (vnp->data.ptrvalue, &st_list, feat_callback);
  }

  if (st_list != NULL) {
    st_list = ValNodeSort (st_list, SortTaxNameConflict);
    st1 = st_list->data.ptrvalue;
    ValNodeAddPointer (&spec_list, st1->obj_type, st1->obj_data);
    have_mismatch = FALSE;
    for (vnp = st_list->next; vnp != NULL; vnp = vnp->next) {
      st2 = vnp->data.ptrvalue;
      if (StringICmp (st1->qual, st2->qual) == 0) {
        ValNodeAddPointer (&spec_list, st2->obj_type, st2->obj_data);
        if (StringICmp (st1->taxname, st2->taxname) != 0) {
          have_mismatch = TRUE;
        }
      } else {
        if (have_mismatch) {
          fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (spec_fmt) + StringLen (qual_name) + StringLen (st1->qual)));
          sprintf (fmt, spec_fmt, qual_name, st1->qual);
          ValNodeAddPointer (&disc_list, 0, NewClickableItem (item_type, fmt, spec_list));
          fmt = MemFree (fmt);
          spec_list = NULL;
        } 
        spec_list = ValNodeFree (spec_list);
        have_mismatch = FALSE;
        ValNodeAddPointer (&spec_list, st2->obj_type, st2->obj_data);
      }
      st1 = st2;
    }
    if (have_mismatch) {
      fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (spec_fmt) + StringLen (qual_name) + StringLen (st1->qual)));
      sprintf (fmt, spec_fmt, qual_name, st1->qual);
      ValNodeAddPointer (&disc_list, 0, NewClickableItem (item_type, fmt, spec_list));
      fmt = MemFree (fmt);
      spec_list = NULL;
    } 
    spec_list = ValNodeFree (spec_list);
    have_mismatch = FALSE;

    if (disc_list != NULL) {
      if (disc_list->next == NULL) {
        ValNodeLink (discrepancy_list, disc_list);
      } else {
        fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (top_fmt) + StringLen (qual_name)));
        sprintf (fmt, top_fmt, qual_name);
        cip = NewClickableItem (item_type, fmt, ItemListFromSubcategories (disc_list));
        fmt = MemFree (fmt);
        cip->subcategories = disc_list;
        ValNodeAddPointer (discrepancy_list, 0, cip);
      }
    }
    st_list = ValNodeFreeData (st_list);
  }

}


static Boolean s_StringHasVoucherSN (CharPtr str)
{
  if (DoesStringContainPhrase (str, "s.n.", TRUE, TRUE)) {
    return TRUE;
  } else if (DoesStringContainPhrase (str, "sn", TRUE, TRUE)) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void CollectSpecVoucherTaxnameCallback (Uint1 obj_type, Pointer obj_data, BioSourcePtr biop, ValNodePtr PNTR list)
{
  OrgModPtr mod;

  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL || list == NULL) {
    return;
  }

  mod = biop->org->orgname->mod;
  while (mod != NULL && (mod->subtype != ORGMOD_specimen_voucher || s_StringHasVoucherSN(mod->subname))) {
    mod = mod->next;
  }
  if (mod != NULL) {
    ValNodeAddPointer (list, 0, TaxNameConflictNew (mod->subname, biop->org->taxname, obj_type, obj_data));
  }
}


static void CollectSpecVoucherTaxnameFeat (SeqFeatPtr sfp, Pointer data)
{
  if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC) {
    CollectSpecVoucherTaxnameCallback (OBJ_SEQFEAT, sfp, sfp->data.value.ptrvalue, data);
  }
}


static void CollectSpecVoucherTaxnameDesc (SeqDescrPtr sdp, Pointer data)
{
  if (sdp != NULL && sdp->choice == Seq_descr_source) {
    CollectSpecVoucherTaxnameCallback (OBJ_SEQDESC, sdp, sdp->data.ptrvalue, data);
  }
}


static void CollectSpecVoucherTaxnameDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  CollectTaxnameConflictDiscrepancies (discrepancy_list, sep_list, 
                                       CollectSpecVoucherTaxnameFeat,
                                       CollectSpecVoucherTaxnameDesc,
                                       "specimen voucher",
                                       DISC_SPECVOUCHER_TAXNAME_MISMATCH);
}


static void CollectStrainTaxnameCallback (Uint1 obj_type, Pointer obj_data, BioSourcePtr biop, ValNodePtr PNTR list)
{
  OrgModPtr mod;

  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL || list == NULL) {
    return;
  }

  mod = biop->org->orgname->mod;
  while (mod != NULL && mod->subtype != ORGMOD_strain) {
    mod = mod->next;
  }
  if (mod != NULL) {
    ValNodeAddPointer (list, 0, TaxNameConflictNew (mod->subname, biop->org->taxname, obj_type, obj_data));
  }
}


static void CollectStrainTaxnameFeat (SeqFeatPtr sfp, Pointer data)
{
  if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC) {
    CollectStrainTaxnameCallback (OBJ_SEQFEAT, sfp, sfp->data.value.ptrvalue, data);
  }
}


static void CollectStrainTaxnameDesc (SeqDescrPtr sdp, Pointer data)
{
  if (sdp != NULL && sdp->choice == Seq_descr_source) {
    CollectStrainTaxnameCallback (OBJ_SEQDESC, sdp, sdp->data.ptrvalue, data);
  }
}


static void CollectStrainTaxnameDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  CollectTaxnameConflictDiscrepancies (discrepancy_list, sep_list, 
                                       CollectStrainTaxnameFeat,
                                       CollectStrainTaxnameDesc,
                                       "strain",
                                       DISC_STRAIN_TAXNAME_MISMATCH);
}


typedef struct partialconflictdata {
  ValNodePtr intron_list;
  ValNodePtr exon_list;
  ValNodePtr promoter_list;
  ValNodePtr RNA_list;
  ValNodePtr utr3_list;
  ValNodePtr utr5_list;
  ValNodePtr cds_list;
  ValNodePtr misc_feature_list; 
} PartialConflictData, PNTR PartialConflictPtr;


static ValNodePtr ReportPartialConflictsForFeatureType (BioseqPtr bsp, Int4 seqfeat, Int4 featdef, CharPtr label)
{
  SeqFeatPtr sfp, gene;
  SeqMgrFeatContext context;
  Boolean partial5, partial3, gene_partial5, gene_partial3;
  SeqLocPtr feat_loc, gene_loc;
  Int4 feat_start, feat_stop, gene_start, gene_stop; 
  Uint1 feat_strand, gene_strand;
  Boolean conflict5, conflict3;
  CharPtr conflict_both_fmt = "%s feature partialness conflicts with gene on both ends";
  CharPtr conflict_5_fmt = "%s feature partialness conflicts with gene on 5' end";
  CharPtr conflict_3_fmt = "%s feature partialness conflicts with gene on 3' end";
  CharPtr fmt;
  ClickableItemPtr cip;
  ValNodePtr disc_list = NULL;

  if (bsp == NULL || ISA_aa (bsp->mol) || label == NULL) {
    return NULL;
  }

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, seqfeat, featdef, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, seqfeat, featdef, &context)) {
    if (sfp->data.choice == SEQFEAT_GENE) {
      continue;
    }
    gene = GetGeneForFeature (sfp);
    if (gene != NULL) {
      feat_loc = SeqLocMerge (bsp, sfp->location, NULL, FALSE, FALSE, FALSE);
      gene_loc = SeqLocMerge (bsp, gene->location, NULL, FALSE, FALSE, FALSE);
      feat_strand = SeqLocStrand (feat_loc);
      if (feat_strand == Seq_strand_minus) {
        feat_start = SeqLocStop (feat_loc);
        feat_stop = SeqLocStart (feat_loc);
      } else {
        feat_start = SeqLocStart (feat_loc);
        feat_stop = SeqLocStop (feat_loc);
      }
      gene_strand = SeqLocStrand (gene_loc);
      if (gene_strand == Seq_strand_minus) {
        gene_start = SeqLocStop (gene_loc);
        gene_stop = SeqLocStart (gene_loc);
      } else {
        gene_start = SeqLocStart (gene_loc);
        gene_stop = SeqLocStop (gene_loc);
      }
      CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
      CheckSeqLocForPartial (gene->location, &gene_partial5, &gene_partial3);
      
      if (((partial5 && !gene_partial5) || (!partial5 && gene_partial5)) && feat_start == gene_start) {
        conflict5 = TRUE;
      } else {
        conflict5 = FALSE;
      }

      if (((partial3 && !gene_partial3) || (!partial3 && gene_partial3)) && feat_stop == gene_stop) {
        conflict3 = TRUE;
      } else {
        conflict3 = FALSE;
      }

      if (conflict5 || conflict3) {
        cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        cip->clickable_item_type = DISC_GENE_PARTIAL_CONFLICT;
        cip->subcategories = NULL;
        cip->item_list = NULL;
        ValNodeAddPointer (&(cip->item_list), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(cip->item_list), OBJ_SEQFEAT, gene);
        if (conflict5 && conflict3) {
          fmt = conflict_both_fmt;
        } else if (conflict5) {
          fmt = conflict_5_fmt;
        } else {
          fmt = conflict_3_fmt;
        }
        cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (label) + StringLen (fmt)));
        sprintf (cip->description, fmt, label);

        cip->callback_func = NULL;
        cip->datafree_func = NULL;
        cip->callback_data = NULL;
        cip->expanded = FALSE;
        cip->level = 0;
        ValNodeAddPointer (&disc_list, 0, cip);
      }
      feat_loc = SeqLocFree (feat_loc);
      gene_loc = SeqLocFree (gene_loc);
    }
  }
  return disc_list;
}
     

static ClickableItemPtr ClickableItemForNumberOfCategories (Uint4 clickable_item_type, CharPtr fmt, ValNodePtr subcategories)
{
  ClickableItemPtr cip;

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  cip->clickable_item_type = clickable_item_type;
  cip->subcategories = subcategories;
  cip->item_list = ItemListFromSubcategories (subcategories);
  cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + 15));
  sprintf (cip->description, fmt, ValNodeLen (subcategories));

  cip->callback_func = NULL;
  cip->datafree_func = NULL;
  cip->callback_data = NULL;
  cip->expanded = FALSE;
  cip->level = 0;
  return cip;
}


static void ReportPartialConflictsBioseqCallback (BioseqPtr bsp, Pointer data)
{
  PartialConflictPtr p;

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL) {
    return;
  }

  p = (PartialConflictPtr) data;

  ValNodeLink (&(p->intron_list), ReportPartialConflictsForFeatureType (bsp, 0, FEATDEF_intron, "intron"));
  ValNodeLink (&(p->exon_list), ReportPartialConflictsForFeatureType (bsp, 0, FEATDEF_exon, "exon"));
  ValNodeLink (&(p->promoter_list), ReportPartialConflictsForFeatureType (bsp, 0, FEATDEF_promoter, "promoter"));
  ValNodeLink (&(p->RNA_list), ReportPartialConflictsForFeatureType (bsp, SEQFEAT_RNA, 0, "RNA"));
  ValNodeLink (&(p->utr3_list), ReportPartialConflictsForFeatureType (bsp, 0, FEATDEF_3UTR, "3' UTR"));
  ValNodeLink (&(p->utr5_list), ReportPartialConflictsForFeatureType (bsp, 0, FEATDEF_5UTR, "5' UTR"));
  if (!IsEukaryotic(bsp)) {
    ValNodeLink (&(p->cds_list), ReportPartialConflictsForFeatureType (bsp, 0, FEATDEF_CDS, "coding region"));
  }
  ValNodeLink (&(p->misc_feature_list), ReportPartialConflictsForFeatureType (bsp, 0, FEATDEF_misc_feature, "misc_feature"));
}


static void ReportPartialConflicts (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, other_list = NULL, disc_list = NULL;
  PartialConflictData p;
  ClickableItemPtr cip;
  ValNodePtr item_list = NULL;

  if (discrepancy_list == NULL || sep_list == NULL) {
    return;
  }

  p.intron_list = NULL;
  p.cds_list = NULL;
  p.exon_list = NULL;
  p.misc_feature_list = NULL;
  p.promoter_list = NULL;
  p.RNA_list = NULL;
  p.utr3_list = NULL;
  p.utr5_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &p, ReportPartialConflictsBioseqCallback);
  }
  
  ValNodeLink (&other_list, p.exon_list);
  ValNodeLink (&other_list, p.intron_list);
  ValNodeLink (&other_list, p.promoter_list);
  ValNodeLink (&other_list, p.RNA_list);
  ValNodeLink (&other_list, p.utr3_list);
  ValNodeLink (&other_list, p.utr5_list);

  if (other_list != NULL) {
    cip = ClickableItemForNumberOfCategories (DISC_GENE_PARTIAL_CONFLICT, "%d features that are not coding regions or misc_features conflict with partialness of overlapping gene", other_list);
    ValNodeLink (&item_list, ItemListFromSubcategories (cip->subcategories));
    ValNodeAddPointer (&disc_list, 0, cip);
  }

  if (p.cds_list != NULL) {
    cip = ClickableItemForNumberOfCategories (DISC_GENE_PARTIAL_CONFLICT, "%d coding region locations conflict with partialness of overlapping gene", p.cds_list);
    ValNodeLink (&item_list, ItemListFromSubcategories (cip->subcategories));
    ValNodeAddPointer (&disc_list, 0, cip);
  }
  if (p.misc_feature_list != NULL) {
    cip = ClickableItemForNumberOfCategories (DISC_GENE_PARTIAL_CONFLICT, "%d misc_feature locations conflict with partialness of overlapping gene", p.misc_feature_list);
    ValNodeLink (&item_list, ItemListFromSubcategories (cip->subcategories));
    ValNodeAddPointer (&disc_list, 0, cip);
  }
  
  if (disc_list != NULL) {
    cip = DiscrepancyForPairs (DISC_GENE_PARTIAL_CONFLICT, "%d feature locations conflict with partialness of overlapping gene", item_list);
    cip->subcategories = disc_list;
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


typedef struct objfindbytext {
  CharPtr search_text;
  ValNodePtr item_list;
} ObjFindByTextData, PNTR ObjFindByTextPtr;


typedef struct objfindlistoftext {
  CharPtr PNTR search_items;
  ValNodePtr PNTR item_lists;
  Boolean whole_word;
} ObjFindListOfTextData, PNTR ObjFindListOfTextPtr;

static void RemoveTranslation (CharPtr str)
{
  CharPtr cp, cp_end;

  cp = StringSearch (str, "/translation=\"");
  if (cp != NULL) {
    cp_end = StringChr (cp + 14, '"');
    if (cp_end != NULL) {
      cp_end++;
      while (*cp_end != 0) {
        *cp = *cp_end;
        cp++;
        cp_end++;
      }
      *cp = 0;
    }
  }
}


static CharPtr GetTaxnameForObject (Uint2 entityID, Uint2 itemtype, Uint4 itemID)
{
  BioseqPtr        bsp = NULL;
  SeqFeatPtr       sfp;
  SeqDescrPtr      sdp;
  ObjValNodePtr    ovn;
  SeqMgrFeatContext fcontext;
  SeqMgrDescContext dcontext;
  BioSourcePtr      biop;
  CharPtr           taxname = NULL;

  switch (itemtype) {
    case OBJ_BIOSEQ:
      bsp =  GetBioseqGivenIDs (entityID, itemID, itemtype);
      break;
    case OBJ_SEQFEAT:
      sfp = SeqMgrGetDesiredFeature (entityID, NULL, itemID, 0, NULL, &fcontext);          
      if (sfp != NULL) {
        bsp = BioseqFindFromSeqLoc (sfp->location);
      }
      break;
    case OBJ_SEQDESC:
      sdp = SeqMgrGetDesiredDescriptor (entityID, NULL, itemID, 0, NULL, &dcontext);
      if (sdp != NULL && sdp->extended != 0) {
        ovn = (ObjValNodePtr) sdp;
        if (ovn->idx.parenttype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) ovn->idx.parentptr;
        } else if (ovn->idx.parenttype == OBJ_BIOSEQSET && ovn->idx.parentptr != NULL) {
          bsp = GetRepresentativeBioseqFromBioseqSet (ovn->idx.parentptr);
        }
      }
      break;
  }
  if (bsp != NULL) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
    if (sdp != NULL) {
      biop = sdp->data.ptrvalue;
      if (biop != NULL && biop->org != NULL) {
        taxname = biop->org->taxname;
      }
    }
  }
  return taxname;
}


static void FlatfileTextFind (
  CharPtr str,
  Pointer userdata,
  BlockType blocktype,
  Uint2 entityID,
  Uint2 itemtype,
  Uint4 itemID,
  Int4 left,
  Int4 right
)

{
  BioseqPtr        bsp;
  SeqFeatPtr       sfp;
  SeqDescrPtr      sdp;
  SeqMgrFeatContext fcontext;
  SeqMgrDescContext dcontext;
  ObjFindListOfTextPtr  obj;
  CharPtr               cpy;
  Int4                  i;
  CharPtr               taxname = NULL;
  Boolean               do_add;

  if (blocktype == SEQUENCE_BLOCK) return;
  /* don't spellcheck organism name or lineage */
  if (blocktype == ORGANISM_BLOCK) return;
  if (userdata == NULL) return;

  obj = (ObjFindListOfTextPtr) userdata;

  cpy = StringSave (str);

  if (blocktype == FEATURE_BLOCK) {
    RemoveTranslation (cpy);
  }

  for (i = 0; obj->search_items[i] != NULL; i++) {
    do_add = FALSE;
    if (DoesStringContainPhrase (cpy, obj->search_items[i], FALSE, obj->whole_word)) {
      if (taxname == NULL) {
        /* remove taxname */
        taxname = GetTaxnameForObject (entityID, itemtype, itemID);
        FindReplaceString (&cpy, taxname, "", FALSE, TRUE);
        if (DoesStringContainPhrase (cpy, obj->search_items[i], FALSE, obj->whole_word)) {
          do_add = TRUE;
        }
      } else {
        do_add = TRUE;
      }
    }
    if (do_add) {
      switch (itemtype) {
        case OBJ_BIOSEQ:
          bsp =  GetBioseqGivenIDs (entityID, itemID, itemtype);
          if (bsp != NULL) {
            ValNodeAddPointer (&(obj->item_lists[i]), OBJ_BIOSEQ, bsp);
          }
          break;
        case OBJ_SEQFEAT:
          sfp = SeqMgrGetDesiredFeature (entityID, NULL, itemID, 0, NULL, &fcontext);
          if (sfp != NULL) {
            ValNodeAddPointer (&(obj->item_lists[i]), OBJ_SEQFEAT, sfp);
          }
          break;
        case OBJ_SEQDESC:
          sdp = SeqMgrGetDesiredDescriptor (entityID, NULL, itemID, 0, NULL, &dcontext);
          if (sdp != NULL) {
            ValNodeAddPointer (&(obj->item_lists[i]), OBJ_SEQDESC, sdp);
          }
          break;
      }
    }
  }
  cpy = MemFree (cpy);
}


static void FindTextInFlatfileEx (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list, Uint4 item_type, CharPtr PNTR find_list, Boolean whole_word)
{
  XtraBlock       xtra;
  ObjFindListOfTextData od;
  ErrSev          level;
  Boolean         okay;
  SeqEntryPtr     oldscope;
  SeqEntryPtr     sep;
  CharPtr         find_fmt = "%%d objects contain %s", fmt;
  Int4            i, num = 0;
  ValNodePtr      vnp;

  if (discrepancy_list == NULL || sep_list == NULL) return;

  MemSet ((Pointer) &xtra, 0, sizeof (XtraBlock));
  xtra.ffwrite = FlatfileTextFind;
  xtra.userdata = (Pointer) &od;
  xtra.reindex = TRUE;
  level = ErrSetMessageLevel (SEV_MAX);

  od.whole_word = whole_word;
  od.search_items = find_list;
  for (i = 0; find_list[i] != NULL; i++) {
    num++;
  }

  od.item_lists = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * num);
  for (i = 0; find_list[i] != NULL; i++) {
    od.item_lists[i] = NULL;
  }

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    oldscope = SeqEntrySetScope (sep);
    okay = SeqEntryToGnbk (sep, NULL, GENBANK_FMT, SEQUIN_MODE, NORMAL_STYLE,
                        SHOW_CONTIG_FEATURES, 0, 0, &xtra, NULL);
    SeqEntrySetScope (oldscope);
  }
  for (i = 0; find_list[i] != NULL; i++) {
    if (od.item_lists[i] != NULL) {
      fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (find_fmt) + StringLen (find_list[i])));
      sprintf (fmt, find_fmt, find_list[i]);
      ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (item_type, fmt, od.item_lists[i]));
      od.item_lists[i] = NULL;
      fmt = MemFree (fmt);
    }
  }

  od.item_lists = MemFree (od.item_lists);

  ErrSetMessageLevel (level);
}


static CharPtr flatfile_find_list_oncaller[] = {
 "univeristy",
 "univerisity",
 "univercity",
 "uiniversity",
 "uinversity",
 "univesity",
 "uviversity",
 "putatvie",
 "putaitve",
 "protien",
 "simmilar",
 "Insitiute",
 "Instutite",
 "instute",
 "institue",
 "insitute",
 "insititute",
 "ribosoml",
 "transcirbed",
 "Agricultutral",
 "agriculturral",
 "resaerch",
 "charaterization",
 "clonging",
 "anaemia",
 "heam",
 "haem",
 "technlogy",
 "technolgy",
 "biotechnlogy",
 "biotechnolgy",
 "biotechology",
 "enviroment",
 "hypotetical",
 "puatative",
 "putaive",
 "putatitve",
 "putataive",
 "putatuve",
 "cotaining",
 "hypothteical",
 "hypotethical",
 "hypothetcial",
 "consevered",
 "haemagglutination",
 "indepedent",
 "reserch",
 "agricultral",
 "Bacilllus",
 "catalize",
 "subitilus",
 "P.R.Chian",
 "PRChian",
 "phylogentic",
 "pylogeny",
 "reseach",
 "ribossomal",
 "mithocon",
 "scencies",
 "scinece",
 "enivronment",
 "structual",
 "sulfer",
 NULL
};


static CharPtr flatfile_find_list_oncaller_wholeword[] = {
  "caputre",
  "casette",
  "chian",
  "cytochome",
  "diveristy",
  "genone",
  "muesum",
  "musuem",
  "nuclear shutting",
  "reserach",
  "transcirption",
  "unversity",
  "varent",
  NULL
};


static void FindTextInFlatfileOncaller (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  FindTextInFlatfileEx (discrepancy_list, sep_list, DISC_FLATFILE_FIND_ONCALLER, flatfile_find_list_oncaller, FALSE);
  FindTextInFlatfileEx (discrepancy_list, sep_list, DISC_FLATFILE_FIND_ONCALLER, flatfile_find_list_oncaller_wholeword, TRUE);
}


typedef struct replacepair {
  CharPtr find;
  CharPtr replace;
} ReplacePairData, PNTR ReplacePairPtr;


static ReplacePairData oncaller_tool_spell_fixes[] = {
  {"homologue", "homolog" },
  {"charaterization", "characterization"},
  {NULL, NULL}};

static void OncallerToolSpellFix (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr entityID_list = NULL, vnp;
  Uint2 entityID;
  Int4 i;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    entityID = GetEntityIdFromObject (vnp->choice, vnp->data.ptrvalue);
    if (entityID != 0) {
      ValNodeAddInt (&entityID_list, 0, entityID);
    }
  }
  entityID_list = ValNodeSort (entityID_list, SortByIntvalue);
  ValNodeUnique (&entityID_list, SortByIntvalue, ValNodeFree);
  for (vnp = entityID_list; vnp != NULL; vnp = vnp->next) {
    entityID = vnp->data.intvalue;
    for (i = 0; oncaller_tool_spell_fixes[i].find != NULL; i++) {
      FindReplaceInEntity (entityID, oncaller_tool_spell_fixes[i].find, oncaller_tool_spell_fixes[i].replace, 
                           FALSE, FALSE, TRUE,
                              FALSE, UPDATE_NEVER, NULL, NULL, NULL, FALSE, NULL, NULL);
    }
  }
}


static SuspectProductNameData cds_product_find[] = {
  { "-like", EndsWithPattern } ,
  { "pseudo", ContainsPseudo } ,
  { "fragment", ContainsWholeWord } ,
  { "similar", ContainsWholeWord } ,
  { "frameshift", ContainsWholeWord } ,
  { "partial", ContainsWholeWord } ,
  { "homolog", ContainsWholeWord } ,
  { "homologue", ContainsWholeWord } ,
  { "paralog", ContainsWholeWord } ,
  { "paralogue", ContainsWholeWord } ,
  { "ortholog", ContainsWholeWord } ,
  { "orthologue", ContainsWholeWord } ,
  { "gene", ContainsWholeWord } ,
  { "genes", ContainsWholeWord } ,
  { "related", ContainsWholeWord } ,
  { "terminus", ContainsWholeWord } ,
  { "N-terminus", ContainsWholeWord } ,
  { "C-terminus", ContainsWholeWord } ,
  { "characterised", ContainsWholeWord } ,
  { "recognised", ContainsWholeWord } ,
  { "characterisation", ContainsWholeWord } ,
  { "localisation", ContainsWholeWord } ,
  { "tumour", ContainsWholeWord } ,
  { "uncharacterised", ContainsWholeWord } ,
  { "oxydase", ContainsWholeWord } ,
  { "colour", ContainsWholeWord } ,
  { "localise", ContainsWholeWord } ,
  { "faecal", ContainsWholeWord } ,
  { "frame"}
};

const int num_cds_product_find = sizeof (cds_product_find) / sizeof (SuspectProductNameData);

static void FindCodingRegions (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR feature_list;
  Int4            k;
  ProtRefPtr      prp;
  ValNodePtr      vnp;
  BioseqPtr       bsp;
  SeqFeatPtr      cds;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT || sfp->data.value.ptrvalue == NULL
      || userdata == NULL)
  {
    return;
  }
  
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  feature_list = (ValNodePtr PNTR) userdata;

  /* add coding region rather than protein */
  if (sfp->idx.subtype == FEATDEF_PROT) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
      if (cds != NULL) {
        sfp = cds;
      }
    }
  }
  
  for (k = 0; k < num_cds_product_find; k++)
  {
    for (vnp = prp->name; vnp != NULL; vnp = vnp->next) 
    {
      if (cds_product_find[k].search_func != NULL
        && (cds_product_find[k].search_func) (cds_product_find[k].pattern, vnp->data.ptrvalue)) 
      {
        ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
        break;
      }
    }
  }
}


static void FindTextInCDSProduct (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr PNTR feature_list;
  ValNodePtr         master_list = NULL, vnp;
  Int4               k;
  ClickableItemPtr dip;
  ValNodePtr         subcategories = NULL;

  if (discrepancy_list == NULL || sep_list == NULL) {
    return;
  }

  feature_list = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * num_cds_product_find);
  if (feature_list == NULL) return;
  
  /* initialize array for suspicious product names */
  for (k = 0; k < num_cds_product_find; k++)
  {
    feature_list[k] = NULL;
  }

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, feature_list, FindCodingRegions);
  }
  
  for (k = 0; k < num_cds_product_find; k++)
  {
    if (feature_list[k] != NULL)
    {
      if (cds_product_find[k].search_func == EndsWithPattern) 
      {
        dip = SuspectPhraseEnd (DISC_CDS_PRODUCT_FIND, cds_product_find[k].pattern, "coding region product", feature_list[k]);
      }
      else if (cds_product_find[k].search_func == StartsWithPattern) 
      {
        dip = SuspectPhraseStart (DISC_CDS_PRODUCT_FIND, cds_product_find[k].pattern, "coding region product", feature_list[k]);
      }
      else 
      {
        dip = SuspectPhrase (DISC_CDS_PRODUCT_FIND, cds_product_find[k].pattern, "coding region product", feature_list[k]);
      }
      if (dip != NULL)
      {
        ValNodeAddPointer (&subcategories, 0, dip);
      }
      ValNodeLinkCopy (&master_list, feature_list[k]);
    }
  }
  
  if (master_list != NULL)
  {
    dip = SuspectPhraseEx (DISC_CDS_PRODUCT_FIND, "suspect phrase or characters", FALSE, "coding region product", master_list);
    if (dip != NULL)
    {
      dip->subcategories = subcategories;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }

  MemFree (feature_list);
}


static void FindDupDeflineCallback (BioseqPtr bsp, Pointer data)
{
  SeqDescrPtr sdp;

  if (bsp != NULL && data != NULL && !ISA_aa (bsp->mol)) {
    for (sdp = bsp->descr; sdp != NULL; sdp = sdp->next) {
      if (sdp->choice == Seq_descr_title) {
        ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
      }
    }
  }
}


static int LIBCALLBACK SortVnpByDescriptorText (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  SeqDescrPtr sdp1, sdp2;
  int rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      sdp1 = vnp1->data.ptrvalue;
      sdp2 = vnp2->data.ptrvalue;
      if (sdp1 != NULL && sdp2 != NULL) {
        rval = StringICmp (sdp1->data.ptrvalue, sdp2->data.ptrvalue);
      }
    }
  }
  return rval;
}


static void FindDupDeflines (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr title_list = NULL, vnp, repeated = NULL, subcat = NULL, unique = NULL;
  SeqDescrPtr sdp1, sdp2;
  ClickableItemPtr cip;
  CharPtr dup_fmt = "%d definition lines are identical";
  CharPtr unique_fmt = "%d definition lines are unique";
  Boolean any_errors = FALSE;

  if (discrepancy_list == NULL || sep_list == NULL) {
    return;
  }

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &title_list, FindDupDeflineCallback);
  }

  title_list = ValNodeSort (title_list, SortVnpByDescriptorText);
  if (title_list != NULL && title_list->next != NULL) {
    sdp1 = title_list->data.ptrvalue;
    ValNodeAddPointer (&repeated, OBJ_SEQDESC, sdp1);
    for (vnp = title_list->next; vnp != NULL; vnp = vnp->next) {
      sdp2 = vnp->data.ptrvalue;
      if (StringICmp (sdp1->data.ptrvalue, sdp2->data.ptrvalue) != 0) {
        if (repeated->next == NULL) {
          ValNodeLink (&unique, repeated);
        } else {
          ValNodeAddPointer (&subcat, 0, NewClickableItem (DISC_DUP_DEFLINE, dup_fmt, repeated));
        }
        repeated = NULL;
      }
      ValNodeAddPointer (&repeated, OBJ_SEQDESC, sdp2);
      sdp1 = sdp2;      
    }
    if (repeated->next == NULL) {
      ValNodeLink (&unique, repeated);
    } else {
      ValNodeAddPointer (&subcat, 0, NewClickableItem (DISC_DUP_DEFLINE, dup_fmt, repeated));
    }
    if (subcat != NULL) {
      if (unique != NULL) {
        ValNodeAddPointer (&subcat, 0, NewClickableItem (DISC_DUP_DEFLINE, unique_fmt, unique));
        unique = NULL;
      }
      if (subcat->next == NULL) {
        ValNodeLink (discrepancy_list, subcat);
      } else {
        cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        MemSet (cip, 0, sizeof (ClickableItemData));
        cip->clickable_item_type = DISC_DUP_DEFLINE;
        cip->description = StringSave ("Defline Problem Report");
        cip->subcategories = subcat;
        ValNodeAddPointer (discrepancy_list, 0, cip);
      }
      any_errors = TRUE;
    }
  }
  title_list = ValNodeFree (title_list);
  unique = ValNodeFree (unique);
  if (!any_errors) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->clickable_item_type = DISC_DUP_DEFLINE;
    cip->description = StringSave ("All deflines are unique");
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


static void CountNucBioseqCallback (BioseqPtr bsp, Pointer data)
{
  if (bsp != NULL && ISA_na (bsp->mol) && data != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQ, bsp);
  }
}


static void CountNucSeqs (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, bsp_list = NULL;

  if (discrepancy_list == NULL || sep_list == NULL) {
    return;
  }

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &bsp_list, CountNucBioseqCallback);
  }

  if (bsp_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_COUNT_NUCLEOTIDES, "%d nucleotide Bioseqs are present", bsp_list));
  }
}


static Boolean HasCultureCollectionForATCCStrain (OrgModPtr mods, CharPtr str)
{
  OrgModPtr mod;
  CharPtr   cp;
  Boolean   rval = FALSE;

  if (StringHasNoText (str)) {
    return TRUE;
  } else if (mods == NULL) {
    return FALSE;
  }

  for (mod = mods; mod != NULL && !rval; mod = mod->next) {
    if (mod->subtype == ORGMOD_culture_collection
      && StringNCmp (mod->subname, "ATCC:", 5) == 0) {
      cp = StringChr (str, ';');
      if (cp == NULL) {
        if (StringCmp (mod->subname + 5, str) == 0) {
          rval = TRUE;
        }
      } else if (StringNCmp (mod->subname + 5, str, cp - str) == 0) {
        rval = TRUE;
      }
    }
  }
  return rval;
}


static Boolean HasATCCStrainForCultureCollection (OrgModPtr mods, CharPtr str)
{
  OrgModPtr mod;
  CharPtr   cp;
  Boolean   rval = FALSE;

  if (StringHasNoText (str)) {
    return TRUE;
  } else if (mods == NULL) {
    return FALSE;
  }

  for (mod = mods; mod != NULL && !rval; mod = mod->next) {
    if (mod->subtype == ORGMOD_strain
      && StringNCmp (mod->subname, "ATCC ", 5) == 0) {
      cp = StringChr (mod->subname, ';');
      if (cp == NULL) {
        if (StringCmp (mod->subname + 5, str) == 0) {
          rval = TRUE;
        }
      } else if (StringNCmp (mod->subname + 5, str, cp - mod->subname - 5) == 0) {
        rval = TRUE;
      }
    }
  }
  return rval;
}


typedef Boolean (*CollectBioSourceTest) PROTO ((BioSourcePtr));

typedef struct collectbiosource {
  CollectBioSourceTest test_func;
  ValNodePtr pass_list;
  ValNodePtr fail_list;
} CollectBioSourceData, PNTR CollectBioSourcePtr;


static void CollectBioSourceDescCallback (SeqDescrPtr sdp, Pointer data)
{
  CollectBioSourcePtr cb;

  if (sdp != NULL && sdp->choice == Seq_descr_source 
      && (cb = (CollectBioSourcePtr)data) != NULL
      && cb->test_func != NULL) {
    if ((cb->test_func) (sdp->data.ptrvalue)) {
      ValNodeAddPointer (&(cb->pass_list), OBJ_SEQDESC, sdp);
    } else {
      ValNodeAddPointer (&(cb->fail_list), OBJ_SEQDESC, sdp);
    }
  }
}


static void CollectBioSourceFeatCallback (SeqFeatPtr sfp, Pointer data)
{
  CollectBioSourcePtr cb;

  if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC 
      && (cb = (CollectBioSourcePtr)data) != NULL
      && cb->test_func != NULL) {
    if ((cb->test_func)(sfp->data.value.ptrvalue)) {
      ValNodeAddPointer (&(cb->pass_list), OBJ_SEQFEAT, sfp);
    } else {
      ValNodeAddPointer (&(cb->fail_list), OBJ_SEQFEAT, sfp);
    }
  }
}


static ValNodePtr CollectBioSources (ValNodePtr sep_list, CollectBioSourceTest test_func, Boolean want_pass)
{
  CollectBioSourceData cb;
  ValNodePtr vnp;

  cb.test_func = test_func;
  cb.pass_list = NULL;
  cb.fail_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitDescriptorsInSep (vnp->data.ptrvalue, &cb, CollectBioSourceDescCallback);
    VisitFeaturesInSep (vnp->data.ptrvalue, &cb, CollectBioSourceFeatCallback);
  }

  if (want_pass) {
    cb.fail_list = ValNodeFree (cb.fail_list);
    return cb.pass_list;
  } else {
    cb.pass_list = ValNodeFree (cb.pass_list);
    return cb.fail_list;
  }
}


static Boolean IsATCCStrainInCultureCollectionForBioSource (BioSourcePtr biop)
{
  OrgModPtr mod;
  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) {
    return TRUE;
  }
  for (mod = biop->org->orgname->mod; mod != NULL; mod = mod->next) {
    if (mod->subtype == ORGMOD_strain && StringNCmp (mod->subname, "ATCC ", 5) == 0) {
      if (!HasCultureCollectionForATCCStrain(biop->org->orgname->mod, mod->subname + 5)) {
        return FALSE;
      }
    } else if (mod->subtype == ORGMOD_culture_collection && StringNCmp (mod->subname, "ATCC:", 5) == 0) {
      if (!HasATCCStrainForCultureCollection (biop->org->orgname->mod, mod->subname + 5)) {
        return FALSE;
      }
    }
  }
  
  return TRUE;
}


static void CheckATCCStrainCultureCollConflict (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr src_list = NULL;

  src_list = CollectBioSources (sep_list, IsATCCStrainInCultureCollectionForBioSource, FALSE);

  if (src_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DUP_DISC_ATCC_CULTURE_CONFLICT, "%d biosources have conflicting ATCC strain and culture collection values", src_list));
  }
}


static void AddATCCStrainToCultureColl (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  AECRParseActionPtr  parse;
  SourceQualPairPtr   pair;
  ValNodePtr          field_from, field_to, vnp;
  CharPtr             str1, str2, cp, new_str;

  parse = AECRParseActionNew ();

  parse->fields = ValNodeNew (NULL);
  parse->fields->choice = FieldPairType_source_qual;
  pair = SourceQualPairNew ();
  pair->field_from = Source_qual_strain;
  pair->field_to = Source_qual_culture_collection;
  parse->fields->data.ptrvalue = pair;

  parse->portion = TextPortionNew ();
  parse->portion->left_marker = ValNodeNew (NULL);
  parse->portion->left_marker = MakeTextTextMarker ("ATCC ");
  parse->portion->include_left = FALSE;
  parse->portion->right_marker = NULL;
  parse->portion->include_right = FALSE;
  parse->portion->inside = TRUE;
  parse->portion->case_sensitive = FALSE;
  parse->portion->whole_word = FALSE;

  parse->remove_from_parsed = FALSE;
  parse->remove_left = FALSE;
  parse->remove_right = FALSE;
  parse->existing_text = ExistingTextOption_add_qual;

  field_from = GetFromFieldFromFieldPair (parse->fields);
  field_to = GetToFieldFromFieldPair (parse->fields);

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    str1 = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_from, NULL);
    str2 = GetTextPortionFromString (str1, parse->portion);
    if (str2 != NULL) {
      cp = StringChr (str2, ';');
      if (cp != NULL) {
        *cp = 0;
      }
      new_str = (CharPtr) MemNew (sizeof (Char) * (5 + StringLen (str2) + 1));
      sprintf (new_str, "ATCC:%s", str2);
      SetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_to, NULL, new_str, parse->existing_text);
      new_str = MemFree (new_str);
    }
    str1 = MemFree (str1);
    str2 = MemFree (str2);
  }
  field_from = FieldTypeFree (field_from);
  field_to = FieldTypeFree (field_to);
  parse = AECRParseActionFree (parse);
}


static ReplacePairData us_state_abbrev_fixes[] = {
 {"AL", "Alabama"},
 {"AK", "Alaska"},
 {"AZ", "Arizona"},
 {"AR", "Arkansas"},
 {"CA", "California"},
 {"CO", "Colorado"},
 {"CT", "Connecticut"},
 {"DE", "Delaware"},
 {"FL", "Florida"},
 {"GA", "Georgia"},
 {"HI", "Hawaii"},
 {"ID", "Idaho"},
 {"IL", "Illinois"},
 {"IN", "Indiana"},
 {"IA", "Iowa"},
 {"KS", "Kansas"},
 {"KY", "Kentucky"},
 {"LA", "Louisiana"},
 {"ME", "Maine"},
 {"MD", "Maryland"},
 {"MA", "Massachusetts"},
 {"MI", "Michigan"},
 {"MN", "Minnesota"},
 {"MS", "Mississippi"},
 {"MO", "Missouri"},
 {"MT", "Montana"},
 {"NE", "Nebraska"},
 {"NV", "Nevada"},
 {"NH", "New Hampshire"},
 {"NJ", "New Jersey"},
 {"NM", "New Mexico"},
 {"NY", "New York"},
 {"NC", "North Carolina"},
 {"ND", "North Dakota"},
 {"OH", "Ohio"},
 {"OK", "Oklahoma"},
 {"OR", "Oregon"},
 {"PA", "Pennsylvania"},
 {"PR", "Puerto Rico"},
 {"RI", "Rhode Island"},
 {"SC", "South Carolina"},
 {"SD", "South Dakota"},
 {"TN", "Tennessee"},
 {"TX", "Texas"},
 {"UT", "Utah"},
 {"VT", "Vermont"},
 {"VA", "Virginia"},
 {"WA", "Washington"}, 
 {"WV", "West Virginia"},
 {"WI", "Wisconsin"},
 {"WY", "Wyoming"},
 {NULL, NULL}
};


static Boolean IsPubdescSubmit (PubdescPtr pdp)
{
  PubPtr pub;

  if (pdp == NULL) {
    return FALSE;
  }
  for (pub = pdp->pub; pub != NULL; pub = pub->next) {
    if (pub->choice == PUB_Sub) {
      return TRUE;
    }
  }
  return FALSE;
}


static void CollectPubsForUSAStateFeatCallback (SeqFeatPtr sfp, Pointer data)
{
  if (sfp != NULL 
      && sfp->data.choice == SEQFEAT_PUB 
      && data != NULL 
      && IsPubdescSubmit (sfp->data.value.ptrvalue)) {

    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
}


static void CollectPubsForUSAStateDescCallback (SeqDescrPtr sdp, Pointer data)
{
  if (sdp != NULL 
      && sdp->choice == Seq_descr_pub
      && data != NULL 
      && IsPubdescSubmit (sdp->data.ptrvalue)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}


static void CheckUSAStates (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, pub_list = NULL, item_list = NULL;
  ValNode field_c, field_s;
  CharPtr country, state;
  Boolean found;
  Int4    i;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitDescriptorsInSep (vnp->data.ptrvalue, &pub_list, CollectPubsForUSAStateDescCallback);
    VisitFeaturesInSep (vnp->data.ptrvalue, &pub_list, CollectPubsForUSAStateFeatCallback);
  }

  field_c.choice = FieldType_pub;
  field_c.data.intvalue = Publication_field_affil_country;
  field_c.next = NULL;

  field_s.choice = FieldType_pub;
  field_s.data.intvalue = Publication_field_affil_sub;
  field_s.next = NULL;

  for (vnp = pub_list; vnp != NULL; vnp = vnp->next) {
    country = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, &field_c, NULL);
    if (StringCmp (country, "USA") == 0) {
      state = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, &field_s, NULL);
      if (StringCmp (state, "Washington DC") == 0) {
        found = TRUE;
      } else if (state == NULL || !isupper (state[0]) || !isupper (state[1]) || state[2] != 0) {
        ValNodeAddPointer (&item_list, vnp->choice, vnp->data.ptrvalue);
      } else {
        found = FALSE;
        for (i = 0; us_state_abbrev_fixes[i].find != NULL && !found; i++) {
          if (StringCmp (us_state_abbrev_fixes[i].find, state) == 0) {
            found = TRUE;
          }
        }
        if (!found && StringCmp ("DC", state) == 0) {
          found = TRUE;
        }
        if (!found) {
          ValNodeAddPointer (&item_list, vnp->choice, vnp->data.ptrvalue);
        }
      }
      state = MemFree (state);
    }
    country = MemFree (country);
  }

  if (item_list != NULL) {
    item_list = ValNodeSort (item_list, SortVnpByDiscrepancyItemText);
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_USA_STATE, "%d cit-subs are missing state abbreviations", item_list));
  }
}


static void FixUSAStates (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNode    field_c, field_s;
  ValNodePtr vnp;
  CharPtr    state, country;
  Int4       i;

  field_c.choice = FieldType_pub;
  field_c.data.intvalue = Publication_field_affil_country;
  field_c.next = NULL;

  field_s.choice = FieldType_pub;
  field_s.data.intvalue = Publication_field_affil_sub;
  field_s.next = NULL;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    country = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, &field_c, NULL);
    if (StringCmp (country, "USA") == 0) {
      state = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, &field_s, NULL);
      for (i = 0; us_state_abbrev_fixes[i].find != NULL; i++) {
        if (StringICmp (us_state_abbrev_fixes[i].replace, state) == 0) {
          SetFieldValueForObject (vnp->choice, 
                                  vnp->data.ptrvalue, 
                                  &field_s, NULL, 
                                  us_state_abbrev_fixes[i].find, 
                                  ExistingTextOption_replace_old);
          break;
        }
      }
      state = MemFree (state);
    }
    country = MemFree (country);
  }
}


static void CheckForLinkerSequenceCallback (BioseqPtr bsp, Pointer data)
{
  SeqDescrPtr sdp;
  MolInfoPtr  mip;
  SeqMgrDescContext context;
  Int2              ctr;
  Char              buf1[50];
  CharPtr           cp;
  Int4              tail_len = 0;
  Boolean           found_linker = FALSE;

  if (bsp == NULL || bsp->mol != Seq_mol_rna || data == NULL) {
    return;
  }

  /* only inspect mRNA sequences */
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
  if (sdp == NULL || sdp->data.ptrvalue == NULL) {
    return;
  }
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip->biomol != MOLECULE_TYPE_MRNA) {
    return;
  }

  if (bsp->length < 30) {
    /* not long enough to have poly-a tail */
    return;
  }

  ctr = SeqPortStreamInt (bsp, bsp->length - 30, bsp->length - 1, Seq_strand_plus,
                        STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL,
                        (Pointer) buf1, NULL);
  buf1[ctr] = 0;
  cp = buf1;
  while (*cp != 0 && !found_linker) {
    if (*cp == 'A') {
      tail_len++;
    } else if (tail_len > 20) {
      found_linker = TRUE;
    } else {
      tail_len = 0;
    }
    cp++;
  }

  if (found_linker) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQ, bsp);
  }
}


static void CheckForLinkerSequence (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, CheckForLinkerSequenceCallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_POSSIBLE_LINKER, "%d bioseqs may have linker sequence after the poly-A tail", item_list));
  }
}


static Boolean IsMrnaSequence (BioseqPtr bsp)
{
  SeqDescrPtr sdp;
  MolInfoPtr  mip;
  SeqMgrDescContext dcontext;
  
  if (bsp == NULL || bsp->mol != Seq_mol_rna) {
    return FALSE;
  }

  /* only inspect mRNA sequences */
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp == NULL || sdp->data.ptrvalue == NULL) {
    return FALSE;
  }
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip->biomol != MOLECULE_TYPE_MRNA) {
    return FALSE;
  }
  return TRUE;
}


static void FindExonsOnMrnaCallback (BioseqPtr bsp, Pointer data)
{
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;

  /* only inspect mRNA sequences */
  if (!IsMrnaSequence(bsp) || data == NULL) {
    return;
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_exon, &fcontext);
  if (sfp != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQ, bsp);
  }
}


static void FindExonsOnMrna (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, FindExonsOnMrnaCallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_EXON_ON_MRNA, "%d mRNA bioseqs have exon features", item_list));
  }
}


NLM_EXTERN void RemoveExonsOnMrna (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr entityIDList = NULL, vnp;
  BioseqPtr  bsp;
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;

  if (item_list == NULL) {
    return;
  }

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_exon, &fcontext);
         sfp != NULL;
         sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_exon, &fcontext)) {
      sfp->idx.deleteme = TRUE;
    }
    ValNodeAddInt (&entityIDList, 0, bsp->idx.entityID);
  }

  entityIDList = ValNodeSort (entityIDList, SortByIntvalue);
  ValNodeUnique (&entityIDList, SortByIntvalue, ValNodeFree);
  
  for (vnp = entityIDList; vnp != NULL; vnp = vnp->next) {
    DeleteMarkedObjects (vnp->data.intvalue, 0, NULL);
    ObjMgrSetDirtyFlag (vnp->data.intvalue, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, vnp->data.intvalue, 0, 0);    
  }
  ValNodeFree (entityIDList);
}


static int LIBCALLBACK SortObjectListByPubTitleAndAuthors (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  int         rval = 0;
  ValNode     title_field, auth_field;
  CharPtr     title1, title2, auth1, auth2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      title_field.choice = FieldType_pub;
      title_field.data.intvalue = Publication_field_title;
      title_field.next = NULL;
      title1 = GetFieldValueForObject (vnp1->choice, vnp1->data.ptrvalue, &title_field, NULL);
      title2 = GetFieldValueForObject (vnp2->choice, vnp2->data.ptrvalue, &title_field, NULL);
      rval = StringCmp (title1, title2);
      title1 = MemFree (title1);
      title2 = MemFree (title2);
      if (rval == 0) {
        auth_field.choice = FieldType_pub;
        auth_field.data.intvalue = Publication_field_authors_initials;
        auth_field.next = NULL;
        auth1 = GetFieldValueForObject (vnp1->choice, vnp1->data.ptrvalue, &auth_field, NULL);
        auth2 = GetFieldValueForObject (vnp2->choice, vnp2->data.ptrvalue, &auth_field, NULL);
        rval = StringCmp (auth1, auth2);
        auth1 = MemFree (auth1);
        auth2 = MemFree (auth2);
      }
    }
  }
  return rval;
}


static void CollectPubsForTitleAuthorConflictsFeatCallback (SeqFeatPtr sfp, Pointer data)
{
  if (sfp != NULL 
      && sfp->data.choice == SEQFEAT_PUB 
      && data != NULL 
      && !IsPubdescSubmit (sfp->data.value.ptrvalue)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
}


static void CollectPubsForTitleAuthorConflictsDescCallback (SeqDescrPtr sdp, Pointer data)
{
  if (sdp != NULL 
      && sdp->choice == Seq_descr_pub
      && data != NULL 
      && !IsPubdescSubmit (sdp->data.ptrvalue)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}


static ClickableItemPtr MakeAuthorTitleItem (CharPtr title, CharPtr authors, ValNodePtr list)
{
  CharPtr    title_author_fmt = "%%d articles have title '%s' and author list '%s'";
  CharPtr    fmt;
  ClickableItemPtr cip;

  if (title == NULL) {
    title = "NULL";
  }

  fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (title_author_fmt) + StringLen (title) + StringLen (authors)));;
  sprintf (fmt, title_author_fmt, title, authors);
  cip = NewClickableItem (DISC_TITLE_AUTHOR_CONFLICT, fmt, list);
  fmt = MemFree (fmt);
  return cip;
}


static void CheckForTitleAuthorConflicts (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp;
  ValNodePtr pub_list = NULL, repeated = NULL, author_cluster, author_cluster_list = NULL;
  CharPtr    last_title = NULL, this_title, last_authors, this_authors;
  ValNode    title_field, auth_field;
  Boolean    author_conflict = FALSE;
  ValNodePtr disc_list = NULL;
  CharPtr    fmt, title_fmt = "%%d articles have title '%s' but do not have the same author list";
  ClickableItemPtr cip;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitDescriptorsInSep (vnp->data.ptrvalue, &pub_list, CollectPubsForTitleAuthorConflictsDescCallback);
    VisitFeaturesInSep (vnp->data.ptrvalue, &pub_list, CollectPubsForTitleAuthorConflictsFeatCallback);
  }

  pub_list = ValNodeSort (pub_list, SortObjectListByPubTitleAndAuthors);

  if (pub_list != NULL && pub_list->next != NULL) {
    title_field.choice = FieldType_pub;
    title_field.data.intvalue = Publication_field_title;
    title_field.next = NULL;
    auth_field.choice = FieldType_pub;
    auth_field.data.intvalue = Publication_field_authors_initials;
    auth_field.next = NULL;
    last_title = GetFieldValueForObject (pub_list->choice, pub_list->data.ptrvalue, &title_field, NULL);
    last_authors = GetFieldValueForObject (pub_list->choice, pub_list->data.ptrvalue, &auth_field, NULL);
    author_cluster = NULL;
    author_cluster_list = NULL;
    ValNodeAddPointer (&author_cluster, pub_list->choice, pub_list->data.ptrvalue);
    ValNodeAddPointer (&repeated, pub_list->choice, pub_list->data.ptrvalue);
    for (vnp = pub_list->next; vnp != NULL; vnp = vnp->next) {
      this_title = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, &title_field, NULL);
      this_authors = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, &auth_field, NULL);
      if (StringCmp (last_title, this_title) == 0) {
        ValNodeAddPointer (&repeated, vnp->choice, vnp->data.ptrvalue);
        if (StringCmp (last_authors, this_authors) != 0) {
          ValNodeAddPointer (&author_cluster_list, 0, MakeAuthorTitleItem (last_title, last_authors, author_cluster));
          author_cluster = NULL;
        }
        ValNodeAddPointer (&author_cluster, vnp->choice, vnp->data.ptrvalue);
      } else {
        ValNodeAddPointer (&author_cluster_list, 0, MakeAuthorTitleItem (last_title, last_authors, author_cluster));
        author_cluster = NULL;
        if (author_cluster_list != NULL && author_cluster_list->next != NULL) {
          fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (title_fmt) + StringLen (last_title)));
          sprintf (fmt, title_fmt, last_title);
          cip = NewClickableItem (DISC_TITLE_AUTHOR_CONFLICT, fmt, repeated);
          cip->item_list = ValNodeFree (cip->item_list);
          author_cluster_list = ValNodeSort (author_cluster_list, SortVnpByDiscrepancyDescription);
          ValNodeReverse (&author_cluster_list);
          cip->subcategories = author_cluster_list;
          author_cluster_list = NULL;
          ValNodeAddPointer (&disc_list, 0, cip);
          fmt = MemFree (fmt);
          repeated = NULL;
        } else {
          author_cluster_list = FreeClickableList (author_cluster_list);
          repeated = ValNodeFree (repeated);
        }
        author_conflict = FALSE;
        ValNodeAddPointer (&repeated, vnp->choice, vnp->data.ptrvalue);
        ValNodeAddPointer (&author_cluster, vnp->choice, vnp->data.ptrvalue);
      }
      last_title = MemFree (last_title);
      last_authors = MemFree (last_authors);
      last_title = this_title;
      last_authors = this_authors;
    }

    ValNodeAddPointer (&author_cluster_list, 0, MakeAuthorTitleItem (last_title, last_authors, author_cluster));
    author_cluster = NULL;
    if (author_cluster_list != NULL && author_cluster_list->next != NULL) {
      fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (title_fmt) + StringLen (last_title)));
      sprintf (fmt, title_fmt, last_title);
      cip = NewClickableItem (DISC_TITLE_AUTHOR_CONFLICT, fmt, repeated);
      cip->item_list = ValNodeFree (cip->item_list);
      author_cluster_list = ValNodeSort (author_cluster_list, SortVnpByDiscrepancyDescription);
      ValNodeReverse (&author_cluster_list);
      cip->subcategories = author_cluster_list;
      author_cluster_list = NULL;
      ValNodeAddPointer (&disc_list, 0, cip);
      fmt = MemFree (fmt);
      repeated = NULL;
    } else {
      repeated = ValNodeFree (repeated);
      author_cluster_list = FreeClickableList (author_cluster_list);
    }
    last_title = MemFree (last_title);
    last_authors = MemFree (last_authors);

    if (disc_list != NULL) {
      if (disc_list->next == NULL) {
        ValNodeLink (discrepancy_list, disc_list);
      } else{
        cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        MemSet (cip, 0, sizeof (ClickableItemData));
        cip->clickable_item_type = DISC_TITLE_AUTHOR_CONFLICT;
        cip->description = StringSave ("Publication Title/Author Inconsistencies");
        cip->subcategories = disc_list;
        ValNodeAddPointer (discrepancy_list, 0, cip);
      }
    }
  }
          
  pub_list = ValNodeFree (pub_list);
}

/* note, "de la" needs to be before "de" so that it will be skipped in its entirety */
static CharPtr ShortAuthorNames[] = {
"de la",
"del",
"de",
"da",
"du",
"dos",
"la",
"le",
"van",
"von",
"der",
"den",
"di",
NULL};


static Boolean IsNameCapitalizationOk (CharPtr str)
{
  CharPtr cp;
  Int4    i, len;
  Boolean need_cap = TRUE, rval = TRUE, found;
  Boolean needed_lower = FALSE, found_lower = FALSE;

  if (StringHasNoText(str)) {
    return TRUE;
  }

  cp = str;
  while (*cp != 0 && rval) {
    if (isalpha (*cp)) {
      if (!need_cap) {
        needed_lower = TRUE;
        if (!isupper(*cp)) {
          found_lower = TRUE;
        }
      }
      if (need_cap && !isupper (*cp)) {
        if (cp == str || *(cp - 1) == ' ') {
          /* check to see if this is a short name */
          found = FALSE;
          for (i = 0; ShortAuthorNames[i] != NULL && !found; i++) {
            len = StringLen (ShortAuthorNames[i]);
            if (StringNCmp (cp, ShortAuthorNames[i], len) == 0
                && *(cp + len) == ' ') {
              found = TRUE;
              cp += len - 1;
            }
          }
          if (!found) {
            rval = FALSE;
          }
        } else {
          rval = FALSE;
        }
      }
      need_cap = FALSE;
    } else {
      need_cap = TRUE;
    }
    cp++;
  }
  if (needed_lower && !found_lower) {
    rval = FALSE;
  }
  return rval;
}

static Boolean IsAuthorInitialsCapitalizationOk (CharPtr init)
{
  CharPtr cp;

  if (StringHasNoText (init)) {
    return TRUE;
  }

  cp = init;
  while (*cp != 0) {
    if (isalpha (*cp) && !isupper(*cp)) {
      return FALSE;
    }
    cp++;
  }
  return TRUE;
}


static void CheckAuthCapsAuthCallback (NameStdPtr nsp, Pointer userdata)
{
  BoolPtr pIsBad;

  if (nsp == NULL || (pIsBad = (BoolPtr)userdata) == NULL || *pIsBad) {
    return;
  }

  if (!IsNameCapitalizationOk (nsp->names[0])) {
    /* last name bad */
    *pIsBad = TRUE;
  } else if(!IsNameCapitalizationOk (nsp->names[1])) {
    /* first name bad */
    *pIsBad = TRUE;
  } else if(!IsAuthorInitialsCapitalizationOk (nsp->names[4])) {
    /* initials bad */
    *pIsBad = TRUE;
  }
}


static Boolean AreBadAuthCapsInPubdesc (PubdescPtr pubdesc)
{
  Boolean is_bad = FALSE;

  if (pubdesc == NULL) {
    return FALSE;
  }
  VisitAuthorsInPub (pubdesc, &is_bad, CheckAuthCapsAuthCallback);
  return is_bad;
}


static void CheckAuthCapsFeatCallback (SeqFeatPtr sfp, Pointer data)
{
  if (data != NULL && sfp != NULL && sfp->data.choice == SEQFEAT_PUB && AreBadAuthCapsInPubdesc (sfp->data.value.ptrvalue)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
}

static void CheckAuthCapsDescrCallback (SeqDescrPtr sdp, Pointer data)
{
  if (data != NULL && sdp != NULL && sdp->choice == Seq_descr_pub && AreBadAuthCapsInPubdesc (sdp->data.ptrvalue)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}


static Boolean AreAuthCapsOkInSubmitBlock (SubmitBlockPtr sbp)
{
  Boolean is_bad = FALSE;
  AuthorPtr    ap;
  PersonIdPtr  pid;
  ValNodePtr vnp;

  if (sbp == NULL || sbp->cit == NULL || sbp->cit->authors == NULL || sbp->cit->authors->choice != 1) {
    return TRUE;
  }
  for (vnp = sbp->cit->authors->names; vnp != NULL && !is_bad; vnp = vnp->next) {
    if ((ap = (AuthorPtr) vnp->data.ptrvalue) != NULL
        && (pid = (PersonIdPtr) ap->name) != NULL
        && pid->choice == 2) {
      CheckAuthCapsAuthCallback (pid->data, &is_bad);
    }
  }
  return !is_bad;
}


static void CheckAuthCaps (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp;
  ValNodePtr pub_list = NULL;
  SeqEntryPtr sep;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitDescriptorsInSep (sep, &pub_list, CheckAuthCapsDescrCallback);
    VisitFeaturesInSep (sep, &pub_list, CheckAuthCapsFeatCallback);
    if (!AreAuthCapsOkInSubmitBlock(FindSubmitBlockForSeqEntry (sep))) {
      if (IS_Bioseq (sep)) {
        ValNodeAddPointer (&pub_list, OBJ_BIOSEQ, sep->data.ptrvalue);
      } else if (IS_Bioseq_set (sep)) {
        ValNodeAddPointer (&pub_list, OBJ_BIOSEQSET, sep->data.ptrvalue);
      }
    }
  }

  if (pub_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_CHECK_AUTH_CAPS, "%d pubs have incorrect author capitalization", pub_list));
  }
}


static void FixAuthCapsAuthCallback (NameStdPtr nsp, Pointer userdata)
{
  CharPtr cp;

  if (nsp == NULL) {
    return;
  }

  FixCapitalizationInElement (&(nsp->names[0]), FALSE, FALSE, TRUE);
  FixCapitalizationInElement (&(nsp->names[1]), FALSE, FALSE, FALSE);
  /* Set initials to all caps */
  for (cp = nsp->names[4]; cp != NULL && *cp != 0; cp++)
  {
    *cp = toupper (*cp);
  }
}


static void FixAuthCapsInSubmitBlock (SubmitBlockPtr sbp)
{
  ValNodePtr vnp;
  AuthorPtr  ap;
  PersonIdPtr pid;

  if (sbp == NULL || sbp->cit == NULL || sbp->cit->authors == NULL || sbp->cit->authors->choice != 1) {
    return;
  }
  for (vnp = sbp->cit->authors->names; vnp != NULL; vnp = vnp->next) {
    if ((ap = (AuthorPtr) vnp->data.ptrvalue) != NULL
        && (pid = (PersonIdPtr) ap->name) != NULL
        && pid->choice == 2) {
      FixAuthCapsAuthCallback (pid->data, NULL);
    }
  }
}


static void FixAuthCaps (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;
  BioseqPtr bsp;
  BioseqSetPtr bssp;
  SeqSubmitPtr ssp;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case OBJ_SEQFEAT:
        sfp = vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PUB) {
          VisitAuthorsInPub (sfp->data.value.ptrvalue, NULL, FixAuthCapsAuthCallback);
        }
        break;
      case OBJ_SEQDESC:
        sdp = vnp->data.ptrvalue;
        if (sdp != NULL && sdp->choice == Seq_descr_pub) {
          VisitAuthorsInPub (sdp->data.ptrvalue, NULL, FixAuthCapsAuthCallback);
        }
        break;
      case OBJ_BIOSEQ:
        bsp = (BioseqPtr) vnp->data.ptrvalue;
        if (bsp != NULL && bsp->idx.parentptr != NULL && bsp->idx.parenttype == OBJ_SEQSUB) {
          ssp = bsp->idx.parentptr;
          if (ssp != NULL) {
            FixAuthCapsInSubmitBlock (ssp->sub);
          }
        }
        break;
      case OBJ_BIOSEQSET:
        bssp = (BioseqSetPtr) vnp->data.ptrvalue;
        if (bssp != NULL && bssp->idx.parentptr != NULL && bssp->idx.parenttype == OBJ_SEQSUB) {
          ssp = bssp->idx.parentptr;
          if (ssp != NULL) {
            FixAuthCapsInSubmitBlock (ssp->sub);
          }
        }
        break;
    }
  }
}


static CharPtr suspect_rna_product_names[] = 
{
  "gene",
  "genes"
};

const int num_suspect_rna_product_names = sizeof (suspect_rna_product_names) / sizeof (CharPtr);

static void CheckRNAProductsAndCommentsCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR feature_list;
  Int4            k;
  CharPtr         str;
  BoolPtr         phrase_in_product;
  
  if (sfp == NULL || (feature_list = (ValNodePtr PNTR) userdata) == NULL 
      || (sfp->idx.subtype != FEATDEF_rRNA && sfp->idx.subtype != FEATDEF_tRNA)) {
    return;
  }

  phrase_in_product = (BoolPtr) MemNew (sizeof (Boolean) * num_suspect_rna_product_names);
  for (k = 0; k < num_suspect_rna_product_names; k++) {
    phrase_in_product[k] = FALSE;
  }

  /* check product */
  str = GetRNAProductString (sfp, NULL);
  for (k = 0; k < num_suspect_rna_product_names; k++) {
    if (DoesStringContainPhrase (str, suspect_rna_product_names[k], FALSE, TRUE)) {
      ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
      phrase_in_product[k] = TRUE;
    }
  }
  str = MemFree (str);

  /* check comment */

  for (k = 0; k < num_suspect_rna_product_names; k++) {
    if (!phrase_in_product[k] && DoesStringContainPhrase (sfp->comment, suspect_rna_product_names[k], FALSE, TRUE)) {
      ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
    }
  }
  phrase_in_product = MemFree (phrase_in_product);
}


static void 
CheckForSuspectPhraseByList 
(ValNodePtr PNTR discrepancy_list,
 ValNodePtr sep_list,
 CharPtr PNTR phrase_list,
 Int4         num_phrases,
 VisitFeaturesFunc callback,
 Uint4             item_type,
 CharPtr           phrase_loc)
{
  ValNodePtr vnp, master_list = NULL, subcategories = NULL;
  ValNodePtr PNTR feature_list;
  ClickableItemPtr dip;
  SeqEntryPtr sep;
  Int4        k;

  feature_list = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * num_phrases);

  for (k = 0; k < num_phrases; k++) {
    feature_list[k] = NULL;
  }
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitFeaturesInSep (sep, feature_list, callback);
  }

  for (k = 0; k < num_phrases; k++)
  {
    if (feature_list[k] != NULL)
    {
      dip = SuspectPhrase (item_type, phrase_list[k], phrase_loc, feature_list[k]);
      if (dip != NULL)
      {
        ValNodeAddPointer (&subcategories, 0, dip);
      }
      ValNodeLinkCopy (&master_list, feature_list[k]);
    }
  }
  
  if (master_list != NULL)
  {
    dip = SuspectPhrase (item_type, "suspect phrase", phrase_loc, master_list);
    if (dip != NULL)
    {
      dip->subcategories = subcategories;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }

  MemFree (feature_list);
}


static void CheckRNAProductsAndComments (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  CheckForSuspectPhraseByList (discrepancy_list, sep_list, 
                               suspect_rna_product_names, num_suspect_rna_product_names,
                               CheckRNAProductsAndCommentsCallback,
                               DISC_CHECK_RNA_PRODUCTS_AND_COMMENTS,
                               "RNA product_name or comment");
}


static void CheckMicrosatelliteRepeatTypeCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR item_list;
  Boolean is_microsatellite = FALSE;
  Boolean is_tandem = FALSE;
  GBQualPtr qual;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_repeat_region
    || (item_list = (ValNodePtr PNTR) userdata) == NULL) {
    return;
  }

  for (qual = sfp->qual; qual != NULL && (!is_microsatellite || !is_tandem); qual = qual->next) {
    if (StringCmp (qual->qual, "satellite") == 0) {
      if (StringICmp (qual->val, "microsatellite") == 0
          || StringNICmp (qual->val, "microsatellite:", 15) == 0) {
        is_microsatellite = TRUE;
      }
    } else if (StringCmp (qual->qual, "rpt_type") == 0) {
      if (StringCmp (qual->val, "tandem") == 0) {
        is_tandem = TRUE;
      }
    }
  }

  if (is_microsatellite && !is_tandem) {
    ValNodeAddPointer (item_list, OBJ_SEQFEAT, sfp);
  }
}


static void CheckMicrosatelliteRepeatType (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;
  SeqEntryPtr sep;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitFeaturesInSep (sep, &item_list, CheckMicrosatelliteRepeatTypeCallback);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_MICROSATELLITE_REPEAT_TYPE, "%d microsatellites do not have a repeat type of tandem", item_list));
  }
}


static void AddRepeatTypeTandem (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  GBQualPtr qual;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT) {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      qual = GBQualNew ();
      qual->qual = StringSave ("rpt_type");
      qual->val = StringSave ("tandem");
      qual->next = sfp->qual;
      sfp->qual = qual;
    }
  }
}


static void CheckMitochondrionRequiredCallback (BioseqPtr bsp, Pointer data)
{
  SeqMgrFeatContext fcontext;
  SeqMgrDescContext dcontext;
  SeqFeatPtr   sfp;
  SeqDescrPtr  sdp;
  BioSourcePtr biop;
  Boolean      needs_mitochondrial = FALSE;

  if (bsp == NULL || data == NULL) {
    return;
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_D_loop, &fcontext);
  if (sfp == NULL) {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_misc_feature, &fcontext);
         sfp != NULL && !needs_mitochondrial;
         sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_misc_feature, &fcontext)) {
      if (StringISearch (sfp->comment, "control region") != NULL) {
        needs_mitochondrial = TRUE;
      }
    }
  } else {
    needs_mitochondrial = TRUE;
  }

  if (needs_mitochondrial) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
    if (sdp == NULL || (biop = sdp->data.ptrvalue) == NULL || biop->genome != GENOME_mitochondrion) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQ, bsp);
    }
  }
}


static void CheckMitochondrionRequired (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;
  SeqEntryPtr sep;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitBioseqsInSep (sep, &item_list, CheckMitochondrionRequiredCallback);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_MITOCHONDRION_REQUIRED, "%d bioseqs have D-loop or control region misc_feature, but are do not have mitochondrial source", item_list));
  }
}


static void MakeLocationMitochondrial (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr   vnp;
  BioSourcePtr biop;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    biop = GetBioSourceFromObject (vnp->choice, vnp->data.ptrvalue);
    if (biop != NULL) {
      biop->genome = GENOME_mitochondrion;
    }
  }
}


static Boolean DoesPubdescContainUnpubPubWithoutTitle (PubdescPtr pdp)
{
  ValNodePtr pub;
  Boolean    rval = FALSE;
  Int4       status;
  CharPtr    title;

  if (pdp == NULL) {
    return FALSE;
  }

  for (pub = pdp->pub; pub != NULL && !rval; pub = pub->next) {
    status = GetPubMLStatus (pub);
    if (status == Pub_type_unpublished) {
      title = GetPubFieldFromPub(pub, Publication_field_title, NULL);
      if (StringHasNoText (title) || StringICmp (title, "Direct Submission") == 0) {
        rval = TRUE;
      }
      title = MemFree (title);
    }
  }
  return rval;
}


static void FindUnpubPubsWithoutTitlesFeatCallback (SeqFeatPtr sfp, Pointer data)
{
  if (data != NULL && sfp != NULL && sfp->data.choice == SEQFEAT_PUB && DoesPubdescContainUnpubPubWithoutTitle (sfp->data.value.ptrvalue)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
}


static void FindUnpubPubsWithoutTitlesDescCallback (SeqDescrPtr sdp, Pointer data)
{
  if (data != NULL && sdp != NULL && sdp->choice == Seq_descr_pub && DoesPubdescContainUnpubPubWithoutTitle (sdp->data.ptrvalue)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}

static void FindUnpubPubsWithoutTitles (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;
  SeqEntryPtr sep;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitFeaturesInSep (sep, &item_list, FindUnpubPubsWithoutTitlesFeatCallback);
    VisitDescriptorsInSep (sep, &item_list, FindUnpubPubsWithoutTitlesDescCallback);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_UNPUB_PUB_WITHOUT_TITLE, "%d unpublished pubs have no title", item_list));
  }
}

static Boolean AreIntervalStrandsOk (SeqLocPtr gene, SeqLocPtr feat)
{
  SeqLocPtr sub_gene, sub_feat;
  Boolean   found_match;
  Boolean   found_bad = FALSE;
  Int2      cmp;
  Uint1     feat_strand, gene_strand;

  sub_feat = SeqLocFindNext (feat, NULL);
  while (sub_feat != NULL && !found_bad) {
    found_match = FALSE;
    sub_gene = SeqLocFindNext (gene, NULL);
    while (sub_gene != NULL && !found_match) {
      cmp = SeqLocCompare (sub_feat, sub_gene);
      if (cmp == SLC_A_IN_B || cmp == SLC_A_EQ_B) {
        found_match = TRUE;
        feat_strand = SeqLocStrand (sub_feat);
        gene_strand = SeqLocStrand (sub_gene);
        if (!StrandOk(feat_strand, gene_strand)) {
          found_bad = TRUE;
        }
      }
      sub_gene = SeqLocFindNext (gene, sub_gene);
    }
    sub_feat = SeqLocFindNext (feat, sub_feat);
  }
  return !found_bad;
}


static void CheckGeneFeatureStrandConflictsCallback (BioseqPtr bsp, Pointer data)
{
  SeqFeatPtr gene, sfp;
  SeqMgrFeatContext gene_context, fcontext;
  ClickableItemPtr cip;
  Boolean          is_error;

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL) {
    return;
  }

  for (gene = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &gene_context);
       gene != NULL;
       gene = SeqMgrGetNextFeature (bsp, gene, SEQFEAT_GENE, 0, &gene_context)) {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
         sfp != NULL && fcontext.left <= gene_context.right;
         sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext)) {
      if (sfp->data.choice == SEQFEAT_GENE) {
        continue;
      }
      if (sfp->idx.subtype == FEATDEF_primer_bind) {
        continue;
      }
      if (fcontext.left == gene_context.left || fcontext.right == gene_context.right) {
        is_error = FALSE;
        if (gene_context.mixed_strand) {
          /* trans-splicing - compare each interval */
          is_error = !AreIntervalStrandsOk(gene->location, sfp->location);
        } else if (!StrandOk (fcontext.strand, gene_context.strand)) {
          is_error = TRUE;
        }
        if (is_error) {
          cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
          MemSet (cip, 0, sizeof (ClickableItemData));
          cip->clickable_item_type = DISC_BAD_GENE_STRAND;
          cip->description = StringSave ("Gene and feature strands conflict");
          ValNodeAddPointer (&(cip->item_list), OBJ_SEQFEAT, gene);
          ValNodeAddPointer (&(cip->item_list), OBJ_SEQFEAT, sfp);
          ValNodeAddPointer ((ValNodePtr PNTR) data, 0, cip);
        }
      }
    }
  }    
}


static void CheckGeneFeatureStrandConflicts (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, disc_list = NULL;
  CharPtr    fmt = "%d feature locations conflict with gene location strands";
  ClickableItemPtr cip;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &disc_list, CheckGeneFeatureStrandConflictsCallback);
  }
  if (disc_list != NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->clickable_item_type = DISC_BAD_GENE_STRAND;
    cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + 15));
    sprintf (cip->description, fmt, ValNodeLen (disc_list));
    cip->subcategories = disc_list;
    cip->item_list = ItemListFromSubcategories (disc_list);
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


static void FindRBSWithoutGene (BioseqPtr bsp, Pointer data)
{
  SeqMgrFeatContext context, gcontext;
  SeqFeatPtr        sfp;

  if (bsp == NULL || data == NULL || ISA_aa (bsp->mol)) {
    return;
  }

  /* only report RBSs without genes if there are any genes on the sequence */
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &context);
  if (sfp == NULL) {
    return;
  }

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_RBS, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_RBS, &context)) {
    if (SeqMgrGetOverlappingGene (sfp->location, &gcontext) == NULL) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
    }
  }
}
    

static void CheckForRBSWithoutGene (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;
  CharPtr    fmt = "%d RBS features do not have overlapping genes";

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, FindRBSWithoutGene);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_RBS_WITHOUT_GENE, fmt, item_list));
  }
}


/* use integer with flags
 * 1 : found some with graphs
 * 2 : found some without graphs
 */
static void CheckForQualityScoresCallback (BioseqPtr bsp, Pointer data)
{
  SeqAnnotPtr sap;
  Int4Ptr     p_i;
  Int4        i;

  if (bsp == NULL || ISA_aa (bsp->mol) || (p_i = (Int4Ptr)data) == NULL) {
    return;
  }

  for (sap = bsp->annot; sap != NULL && sap->type != 3; sap = sap->next) {
  }
  i = *p_i;
  if (sap == NULL) {
    i |= 2;
  } else {
    i |= 1;
  }
  *p_i = i;
}


static void CheckForQualityScores (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  Int4 i = 0;
  ValNodePtr vnp;
  ClickableItemPtr cip;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &i, CheckForQualityScoresCallback);
  }

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  MemSet (cip, 0, sizeof (ClickableItemData));
  cip->clickable_item_type = DISC_QUALITY_SCORES;
  if (i == 1) {
    cip->description = StringSave ("Quality scores are present on all sequences.");
  } else if (i == 2) {
    cip->description = StringSave ("Quality scores are missing on all sequences.");
  } else if (i == 3) {
    cip->description = StringSave ("Quality scores are missing on some sequences.");
  }
  ValNodeAddPointer (discrepancy_list, 0, cip);
}


static void InternalTranscribedrRNACallback (SeqFeatPtr sfp, Pointer data)
{
  CharPtr product;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_rRNA || data == NULL) {
    return;
  }

  product = GetRNAProductString (sfp, NULL);

  if (StringISearch (product, "internal") != NULL
      || StringISearch (product, "transcribed") != NULL
      || StringISearch (product, "spacer") != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
}


static void InternalTranscribedSpacerrRNA (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL, vnp;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &item_list, InternalTranscribedrRNACallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_INTERNAL_TRANSCRIBED_SPACER_RRNA, "%d rRNA feature products contain 'internal', 'transcribed', or 'spacer'", item_list));
  }
}


static Int4 DistanceToUpstreamGap (Int4 pos, BioseqPtr bsp)
{
  Int4 last_gap = -1, offset = 0;
  DeltaSeqPtr dsp;
  SeqLitPtr slp;

  if (pos < 0 || bsp == NULL || bsp->repr != Seq_repr_delta) {
    return -1;
  }

  for (dsp = bsp->seq_ext; dsp != NULL; dsp = dsp->next) {
    if (dsp->choice == 1) {
      offset += SeqLocLen (dsp->data.ptrvalue);
    } else if (dsp->choice == 2) {
      slp = (SeqLitPtr) dsp->data.ptrvalue;
      offset += slp->length;
      if (IsDeltaSeqGap (dsp)) {
        last_gap = offset;
      }
    }
    if (offset > pos) {
      if (last_gap > -1) {
        return pos - last_gap;
      } else {
        return -1;
      }
    }
  }

  return -1;
}


static Int4 DistanceToDownstreamGap (Int4 pos, BioseqPtr bsp)
{
  Int4 offset = 0;
  DeltaSeqPtr dsp;
  SeqLitPtr slp;

  if (pos < 0 || bsp == NULL || bsp->repr != Seq_repr_delta) {
    return -1;
  }

  for (dsp = bsp->seq_ext; dsp != NULL; dsp = dsp->next) {
    if (dsp->choice == 1) {
      offset += SeqLocLen (dsp->data.ptrvalue);
    } else if (dsp->choice == 2) {
      slp = (SeqLitPtr) dsp->data.ptrvalue;
      if (IsDeltaSeqGap (dsp) && offset > pos) {
        return offset - pos - 1;
      } else {
        offset += slp->length;
      }
    }
  }

  return -1;
}


static Boolean CouldExtendLeft (BioseqPtr bsp, Int4 pos)
{
  Boolean   rval = FALSE;
  Int4      distance;

  if (pos == 0) {
    rval = FALSE;
  } else if (pos < 3) {
    rval = TRUE;
  } else if (bsp->repr == Seq_repr_delta) {
    /* wasn't close to the sequence end, but perhaps it is close to a gap */
    distance = DistanceToUpstreamGap (pos, bsp);
    if (distance > 0 && distance < 3) {
      rval = TRUE;
    }
  }
  return rval;
}


static Boolean CouldExtendRight (BioseqPtr bsp, Int4 pos)
{
  Boolean   rval = FALSE;
  Int4      distance;

  if (pos == bsp->length - 1) {
    rval = FALSE;
  } else if (pos > bsp->length - 4) {
    rval = TRUE;
  } else if (bsp->repr == Seq_repr_delta) {
    /* wasn't close to the sequence end, but perhaps it is close to a gap */
    distance = DistanceToDownstreamGap (pos, bsp);
    if (distance > 0 && distance < 3) {
      rval = TRUE;
    }
  }

  return rval;
}


NLM_EXTERN Int4 
Extend5PartialSeqIntToEndOrGap 
(SeqIntPtr sint,
 BioseqPtr bsp,
 Boolean   short_only)
{
  Int4      distance = 0;

  if (sint == NULL || bsp == NULL) {
    return FALSE;
  }

  if (sint->strand == Seq_strand_minus) {
    if (sint->if_to != NULL && sint->to != bsp->length - 1) {
      distance = DistanceToDownstreamGap (sint->to, bsp);
      if (distance == 1 || distance == 2 || (distance > -1 && !short_only)) {
        sint->to += distance;
      } else if (!short_only || sint->to > bsp->length - 4) {
        distance = bsp->length - 1 - sint->to;
        sint->to = bsp->length - 1;     
      } else {
        distance = 0;
      }
    }
  } else {
    if (sint->if_from != NULL && sint->from != 0) {
      distance = DistanceToUpstreamGap (sint->from, bsp);
      if (distance == 1 || distance == 2 || (distance > -1 && !short_only)) {
        sint->from -= distance;
      } else if (!short_only || sint->from < 3) {
        distance = sint->from;
        sint->from = 0;
      } else {
        distance = 0;
      }
    }
  }

  return distance;
}


NLM_EXTERN Int4 
Extend3PartialSeqIntToEndOrGap 
(SeqIntPtr sint,
 BioseqPtr bsp,
 Boolean   short_only)
{
  Int4      distance = 0;

  if (sint == NULL || bsp == NULL) {
    return FALSE;
  }

  if (sint->strand == Seq_strand_minus) {
    if (sint->if_from != NULL && sint->from != 0) {
      distance = DistanceToUpstreamGap (sint->from, bsp);
      if (distance == 1 || distance == 2 || (distance > -1 && !short_only)) {
        sint->from -= distance;
      } else if (!short_only || sint->from < 3) {
        distance = sint->from;
        sint->from = 0;
      } else {
        distance = 0;
      }
    }
  } else {
    if (sint->if_to != NULL && sint->to != bsp->length - 1) {
      distance = DistanceToDownstreamGap (sint->to, bsp);
      if (distance == 1 || distance == 2 || (distance > -1 && !short_only)) {
        sint->to += distance;
      } else if (!short_only || sint->to > bsp->length - 4) {
        distance = bsp->length - 1 - sint->to;
        sint->to = bsp->length - 1;     
      } else {
        distance = 0;
      }
    }
  }
  return distance;
}



static Boolean ExtendPartialSeqIntToEndOrGap (SeqIntPtr sint, BioseqPtr bsp)
{
  Boolean rval = FALSE;
  if (Extend5PartialSeqIntToEndOrGap (sint, bsp, TRUE) > 0) {
    rval = TRUE;
  }
  
  if (Extend3PartialSeqIntToEndOrGap (sint, bsp, TRUE) > 0) {
    rval = TRUE;
  }

  return rval;
}


NLM_EXTERN Int4 ExtendSeqLocToEndOrGap (SeqLocPtr slp, BioseqPtr bsp, Boolean end5)
{
  Int4 diff = 0;
  SeqLocPtr slp_index;

  if (slp == NULL || bsp == NULL) return 0;

  switch (slp->choice)
  {
    case SEQLOC_INT:
      if (end5) {
        diff = Extend5PartialSeqIntToEndOrGap (slp->data.ptrvalue, bsp, FALSE);
      } else {
        diff = Extend3PartialSeqIntToEndOrGap (slp->data.ptrvalue, bsp, FALSE);
      }
      break;
    case SEQLOC_MIX:
      case SEQLOC_PACKED_INT:
      if (end5) {
        /* take the first one */
        diff = ExtendSeqLocToEndOrGap (slp->data.ptrvalue, bsp, end5);
      } else {
        /* take the last one */
        for (slp_index = slp->data.ptrvalue; slp_index != NULL && slp_index->next != NULL; slp_index = slp_index->next) {
        }
        if (slp_index != NULL) {
          diff = ExtendSeqLocToEndOrGap (slp_index, bsp, end5);
        }
      }
      break;
  }

  return diff;
}


NLM_EXTERN SeqFeatPtr FindBestProtein (Uint2 entityID, SeqLocPtr product)

{
  SeqFeatPtr        sfp, bestprot = NULL;
  SeqMgrFeatContext context;
  BioseqPtr         bsp;
  SeqLocPtr         slp = NULL;
  
  if (product == NULL) return NULL;
  
  bsp = BioseqFindFromSeqLoc (product);
  if (bsp == NULL) return NULL;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_PROT, 0, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_PROT, 0, &context))
  {
    if (slp == NULL) 
    {
      bestprot = sfp;
      slp = sfp->location;
    } else if (SeqLocCompare (slp, sfp->location) == SLC_A_IN_B) {
      bestprot = sfp;
      slp = sfp->location;
    }
  }
  return bestprot;
}


NLM_EXTERN Boolean RetranslateOneCDS 
( SeqFeatPtr sfp,
  Uint2 entityID,
  Boolean include_stop,
  Boolean no_stop_at_end_of_complete_cds)

{
  SeqFeatPtr    bestprot;
  ByteStorePtr  bs;
  BioseqPtr     bsp;
  Char          ch;
  SeqFeatPtr    gene;
  GeneRefPtr    grp;
  SeqEntryPtr   master;
  MolInfoPtr    mip;
  SeqEntryPtr   old;
  Boolean       partial5;
  Boolean       partial3;
  CharPtr       prot;
  CharPtr       ptr;
  SeqEntryPtr   sep;
  SeqIdPtr      sip;
  ValNodePtr    vnp;
  ProtRefPtr    prp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return TRUE;

  /* bail on pseudo CDS */

  if (sfp->pseudo) return TRUE;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL) {
    if (grp->pseudo) return TRUE;
  } else {
    gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
    if (gene != NULL) {
      if (gene->pseudo) return TRUE;
      grp = (GeneRefPtr) gene->data.value.ptrvalue;
      if (grp != NULL && grp->pseudo) return TRUE;
    }
  }

  if (sfp->location == NULL) return TRUE;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);

  if (sfp->product == NULL) {
    master = NULL;
    old = NULL;
    bsp = GetBioseqGivenSeqLoc (sfp->location, entityID);
    if (bsp != NULL) {
      master = GetBestTopParentForData (entityID, bsp);
    }
    bsp = BioseqNew ();
    if (bsp != NULL) {
      bsp->mol = Seq_mol_aa;
      bsp->repr = Seq_repr_raw;
      bsp->seq_data_type = Seq_code_ncbieaa;
      bsp->length = 0;
      bsp->seq_data = (SeqDataPtr) BSNew (0);
      if (master != NULL) {
        old = SeqEntrySetScope (master);
      }
      bsp->id = MakeNewProteinSeqId (sfp->location, NULL);
      SeqMgrAddToBioseqIndex (bsp);
      if (master != NULL) {
        SeqEntrySetScope (old);
      }
      sep = SeqEntryNew ();
      if (sep != NULL) {
        sep->choice = 1;
        sep->data.ptrvalue = (Pointer) bsp;
        SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
      }
      SetSeqFeatProduct (sfp, bsp);
      if (master != NULL && sep != NULL) {
        AddSeqEntryToSeqEntry (master, sep, TRUE);
      }
    }
  }

  sip = SeqLocId (sfp->product);
  if (sip != NULL) {
    bsp = BioseqFind (sip);
    if (bsp != NULL && ISA_aa (bsp->mol) && bsp->repr == Seq_repr_raw) {
      bestprot = FindBestProtein (entityID, sfp->product);
      bs = ProteinFromCdRegionExWithTrailingCodonHandling (sfp,
                                              include_stop,
                                              FALSE,
                                              no_stop_at_end_of_complete_cds );
      if (bs == NULL) return TRUE;
      prot = BSMerge (bs, NULL);
      bs = BSFree (bs);
      if (prot == NULL) return TRUE;
      ptr = prot;
      ch = *ptr;
      while (ch != '\0') {
        *ptr = TO_UPPER (ch);
        ptr++;
        ch = *ptr;
      }
      bs = BSNew (1000);
      if (bs != NULL) {
        ptr = prot;
        /*
        if (prot [0] == '-') {
          ptr++;
        }
        */
        BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
      }
      bsp->repr = Seq_repr_raw;
      bsp->mol = Seq_mol_aa;
      bsp->seq_data = SeqDataFree (bsp->seq_data, bsp->seq_data_type);
      bsp->seq_data = (SeqDataPtr) bs;
      bsp->seq_data_type = Seq_code_ncbieaa;
      bsp->length = BSLen (bs);
      sep = SeqMgrGetSeqEntryForData (bsp);
      if (sep == NULL) return TRUE;
      if (bestprot == NULL)
      {
        bestprot = CreateNewFeature (sep, NULL, SEQFEAT_PROT, NULL);
        prp = ProtRefNew ();
        bestprot->data.value.ptrvalue = prp;
      }
      if (bestprot != NULL) {
        bestprot->location = SeqLocFree (bestprot->location);
        bestprot->location = CreateWholeInterval (sep);
        SetSeqLocPartial (bestprot->location, partial5, partial3);
        bestprot->partial = (partial5 || partial3);
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
          /*
          } else if (partial) {
            mip->completeness = 2;
          */
          } else {
            mip->completeness = 0;
          }
        }
      }
    }
  }
  return TRUE;
}


NLM_EXTERN Boolean ExtendPartialsToEndOrGap (SeqFeatPtr sfp)
{
  Boolean rval = FALSE;
  SeqIdPtr sip;
  SeqLocPtr slp, slp_start = NULL, slp_stop = NULL;
  BioseqPtr bsp;
  Uint1 strand = Seq_strand_unknown;
  Int4  end5 = 0, diff = 0;
  CdRegionPtr crp;
  SeqFeatPtr  gene = NULL, mrna = NULL;
  SeqMgrFeatContext mrna_context;

  if (sfp == NULL || (slp = sfp->location) == NULL) {
    return FALSE;
  }

  if (sfp->data.choice != SEQFEAT_GENE) {
    gene = GetGeneForFeature (sfp);
  }

  if (sfp->idx.subtype != FEATDEF_mRNA) {
    mrna = SeqMgrGetOverlappingmRNA (sfp->location, &mrna_context);
  }

  if (slp->choice == SEQLOC_INT) {
    slp_start = slp;
    slp_stop = slp;
    sip = SeqLocId (slp);
  } else if (slp->choice == SEQLOC_MIX) {
    if ((sip = SeqLocId (slp)) != NULL) /* can only process if all on one bioseq */ {
      slp_start = slp->data.ptrvalue;
      slp_stop = slp_start;
      while (slp_stop->next != NULL) {
        slp_stop = slp_stop->next;
      }
    }
  }
  if (slp_start == NULL || slp_stop == NULL || slp_start->choice != SEQLOC_INT || slp_stop->choice != SEQLOC_INT) {
    return FALSE;
  }

  bsp = BioseqFind (sip);
  if (bsp == NULL) {
    return FALSE;
  }

  if (sfp->data.choice == SEQFEAT_CDREGION) {
    strand = SeqLocStrand (slp);
    if (strand == Seq_strand_minus) {
      end5 = SeqLocStop (slp);
    } else {
      end5 = SeqLocStart (slp);
    }
  }

  rval = ExtendPartialSeqIntToEndOrGap (slp_start->data.ptrvalue, bsp);
  if (slp_stop != slp_start) {
    rval |= ExtendPartialSeqIntToEndOrGap (slp_stop->data.ptrvalue, bsp);
  }

  if (rval) {
    if (sfp->data.choice == SEQFEAT_CDREGION) {
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp == NULL) {
        crp = CdRegionNew ();
        sfp->data.value.ptrvalue = crp;
      }
      if (strand == Seq_strand_minus) {
        diff = SeqLocStop (sfp->location) - end5;
      } else {
        diff = end5 - SeqLocStart (sfp->location);
      }
      if (diff > 0) {
        switch (crp->frame) {
          case 0:
          case 1:
            crp->frame = diff + 1;
            break;
          case 2:
            if (diff == 1) {
              crp->frame = 1;
            } else if (diff == 2) {
              crp->frame = 3;
            }
            break;
          case 3:
            if (diff == 1) {
              crp->frame = 2;
            } else if (diff == 2) {
              crp->frame = 1;
            }
            break;
        }
      }
    }
    /* retranslate coding region */
    RetranslateOneCDS (sfp, sfp->idx.entityID, TRUE, TRUE);

    /* also extend overlapping gene */
    ExtendPartialsToEndOrGap (gene);

    /* and overlapping mRNA */
    ExtendPartialsToEndOrGap (mrna);
  }

  return rval;
}




static void FindExtendablePartialsCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr  sfp;
  SeqMgrFeatContext fcontext;
  Boolean partialL, partialR, partial5, partial3;

  if (bsp == NULL || ISA_aa (bsp->mol) || userdata == NULL) {
    return;
  }

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
        sfp != NULL;
        sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext)) {
    CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
    if (fcontext.strand == Seq_strand_minus) {
      partialL = partial3;
      partialR = partial5;
    } else {
      partialL = partial5;
      partialR = partial3;
    }
    if ((partialL && CouldExtendLeft (bsp, fcontext.left))
      || (partialR && CouldExtendRight (bsp, fcontext.right))) {
      ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_SEQFEAT, sfp);
    }
  }
}


extern void FindExtendablePartials (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL, vnp;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, FindExtendablePartialsCallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_PARTIAL_PROBLEMS, "%d features have partial ends that do not abut the end of the sequence or a gap, but could be extended by 3 or fewer nucleotides to do so", item_list));
  }
}


static Boolean IsNonExtendableLeft (BioseqPtr bsp, Int4 pos)
{
  Boolean   rval = TRUE;
  Int4      distance;

  if (pos < 3) {
    /* is either at the end or is within extending distance */
    return FALSE;
  } else if (bsp->repr == Seq_repr_delta) {
    /* wasn't close to the sequence end, but perhaps it is close to a gap */
    distance = DistanceToUpstreamGap (pos, bsp);
    if (distance > -1 && distance < 3) {
      rval = FALSE;
    }
  }
  return rval;
}


static Boolean IsNonExtendableRight (BioseqPtr bsp, Int4 pos)
{
  Boolean   rval = TRUE;
  Int4      distance;

  if (pos > bsp->length - 4) {
    /* is either at the end or is within extending distance */
    rval = FALSE;
  } else if (bsp->repr == Seq_repr_delta) {
    /* wasn't close to the sequence end, but perhaps it is close to a gap */
    distance = DistanceToDownstreamGap (pos, bsp);
    if (distance > -1 && distance < 3) {
      rval = FALSE;
    }
  }

  return rval;
}


static const CharPtr kNonExtendableException = "unextendable partial coding region";

static void FindBacterialNonExtendablePartialsCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqDescrPtr sdp;
  SeqFeatPtr  sfp;
  SeqMgrDescContext dcontext;
  SeqMgrFeatContext fcontext;
  BioSourcePtr biop;
  Boolean partialL, partialR, partial5, partial3;

  if (bsp == NULL || ISA_aa (bsp->mol) || userdata == NULL) {
    return;
  }

  /* only perform test if associated organism cannot be identified as eukaryote */
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL || (biop = sdp->data.ptrvalue) == NULL || !IsEukaryoticBioSource(biop)) {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
         sfp != NULL;
         sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext)) {
      CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
      /* skip feature if it already has the exception */
      if (StringISearch (sfp->except_text, kNonExtendableException) != NULL) {
        continue;
      }
      if (fcontext.strand == Seq_strand_minus) {
        partialL = partial3;
        partialR = partial5;
      } else {
        partialL = partial5;
        partialR = partial3;
      }
      if ((partialL && IsNonExtendableLeft (bsp, fcontext.left))
        || (partialR && IsNonExtendableRight (bsp, fcontext.right))) {
        ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_SEQFEAT, sfp);
      }
    }
  }
}


extern void FindBacterialNonExtendablePartials (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL, vnp;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, FindBacterialNonExtendablePartialsCallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_BACTERIAL_PARTIAL_NONEXTENDABLE_PROBLEMS, "%d features have partial ends that do not abut the end of the sequence or a gap, and cannot be extended by 3 or fewer nucleotides to do so", item_list));
  }
}


NLM_EXTERN void AddNonExtendableException (SeqFeatPtr sfp)
{
  CharPtr new_except;
  BioseqPtr bsp;

  if (sfp != NULL && StringISearch (sfp->except_text, kNonExtendableException) == NULL) {
    if (sfp->except_text == NULL) {
      sfp->except_text = StringSave (kNonExtendableException);
    } else {
      new_except = (CharPtr) MemNew (sizeof (Char) * (StringLen (sfp->except_text) + StringLen (kNonExtendableException) + 3));
      sprintf (new_except, "%s; %s", sfp->except_text, kNonExtendableException);
      sfp->except_text = MemFree (sfp->except_text);
      sfp->except_text = new_except;
    }
    sfp->excpt = TRUE;
    if (sfp->data.choice == SEQFEAT_CDREGION && sfp->product != NULL) {
      bsp = BioseqFindFromSeqLoc (sfp->product);
      if (bsp != NULL) {
        UpdateProteinTitle (bsp);
      }
    }
  }
}


static void FindBacterialNonExtendablePartialsWithExceptionsCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqDescrPtr sdp;
  SeqFeatPtr  sfp;
  SeqMgrDescContext dcontext;
  SeqMgrFeatContext fcontext;
  BioSourcePtr biop;
  Boolean partialL, partialR, partial5, partial3;

  if (bsp == NULL || ISA_aa (bsp->mol) || userdata == NULL) {
    return;
  }

  /* only perform test if associated organism cannot be identified as eukaryote */
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL || (biop = sdp->data.ptrvalue) == NULL || !IsEukaryoticBioSource(biop)) {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
         sfp != NULL;
         sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext)) {
      CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
      /* skip feature if it does not have the exception */
      if (StringISearch (sfp->except_text, kNonExtendableException) == NULL) {
        continue;
      }
      if (fcontext.strand == Seq_strand_minus) {
        partialL = partial3;
        partialR = partial5;
      } else {
        partialL = partial5;
        partialR = partial3;
      }
      if ((partialL && IsNonExtendableLeft (bsp, fcontext.left))
        || (partialR && IsNonExtendableRight (bsp, fcontext.right))) {
        ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_SEQFEAT, sfp);
      }
    }
  }
}


static void FindBacterialNonExtendablePartialsWithExceptions (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL, vnp;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, FindBacterialNonExtendablePartialsWithExceptionsCallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_BACTERIAL_PARTIAL_NONEXTENDABLE_EXCEPTION, "%d features have partial ends that do not abut the end of the sequence or a gap, and cannot be extended by 3 or fewer nucleotides to do so, but have the correct exception", item_list));
  }
}


NLM_EXTERN void FixBacterialNonExtendablePartials (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp, gene;
  Boolean    has_title = FALSE;
  CharPtr    orig_location, key;
  GeneRefPtr grp;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT && (sfp = vnp->data.ptrvalue) != NULL) {
      orig_location = SeqLocPrintUseBestID (sfp->location);
      if (sfp->data.choice == SEQFEAT_GENE) {
        gene = sfp;
      } else {
        gene = GetGeneForFeature (sfp);
      }
      AddNonExtendableException (sfp);
      if (gene == NULL) {
        grp = SeqMgrGetGeneXref (sfp);
      } else {
        grp = gene->data.value.ptrvalue;
      }

      key = StringSaveNoNull (FeatDefTypeLabel (sfp));

      if (lip != NULL && lip->fp != NULL) {
        if (!has_title) {
          fprintf (lip->fp, "Exceptions for extendable partials:\n");
          has_title = TRUE;
        }
        if (grp != NULL && !StringHasNoText (grp->locus_tag )) {
          fprintf (lip->fp, "Added exception to %s (%s) at %s\n", key == NULL ? "Unknown feature type" : key,
                                                                grp->locus_tag,
                                                                orig_location);
        } else {
          fprintf (lip->fp, "Added exception to %s at %s \n", key == NULL ? "Unknown feature type" : key,
                                                              orig_location);
        }
      }
      key = MemFree (key);
      if (lip != NULL) {
        lip->data_in_log = TRUE;
      }
      orig_location = MemFree (orig_location);
    }
  }
  if (has_title) {
    fprintf (lip->fp, "\n");
  }
}


NLM_EXTERN void FixExtendablePartials (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp, gene;
  CharPtr    orig_location, new_location, key;
  GeneRefPtr grp;
  Boolean    has_title = FALSE;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT && (sfp = vnp->data.ptrvalue) != NULL) {
      orig_location = SeqLocPrintUseBestID (sfp->location);
      if (sfp->data.choice == SEQFEAT_GENE) {
        gene = sfp;
      } else {
        gene = GetGeneForFeature (sfp);
      }
      if (ExtendPartialsToEndOrGap (sfp) && lip != NULL && lip->fp != NULL) {
        if (!has_title) {
          fprintf (lip->fp, "Extended Partial Features\n");
          has_title = TRUE;
        }
        new_location = SeqLocPrintUseBestID (sfp->location);
        if (gene == NULL) {
          grp = SeqMgrGetGeneXref (sfp);
        } else {
          grp = gene->data.value.ptrvalue;
        }

        key = StringSaveNoNull (FeatDefTypeLabel (sfp));

        if (grp != NULL && !StringHasNoText (grp->locus_tag )) {
          fprintf (lip->fp, "Extended %s (%s) from %s to %s\n", key == NULL ? "Unknown feature type" : key,
                                                                grp->locus_tag,
                                                                orig_location, new_location);
        } else {
          fprintf (lip->fp, "Extended %s %s to %s\n", key == NULL ? "Unknown feature type" : key,
                                                      orig_location, new_location);
        }
        key = MemFree (key);
        new_location = MemFree (new_location);
        lip->data_in_log = TRUE;
      }
      orig_location = MemFree (orig_location);
    }
  }
  if (has_title) {
    fprintf (lip->fp, "\n");
  }
}


static CharPtr suspect_rrna_product_names[] = 
{
"domain",
"partial",
"5s_rRNA",
"16s_rRNA",
"23s_rRNA"
};

const int num_suspect_rrna_product_names = sizeof (suspect_rrna_product_names) / sizeof (CharPtr);

static void FindSuspectrRNAProductsCallback (SeqFeatPtr sfp, Pointer data)
{
  Int4 k;
  CharPtr product;
  ValNodePtr PNTR feature_list;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_rRNA || (feature_list = (ValNodePtr PNTR)data) == NULL) {
    return;
  }

  product = GetRNAProductString (sfp, NULL);
  for (k = 0; k < num_suspect_rrna_product_names; k++) {
    if (DoesStringContainPhrase (product, suspect_rrna_product_names[k], FALSE, FALSE)) {
      ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
    }
  }
  product = MemFree (product);
}


static StringConstraintPtr MakeSimpleSearchConstraint (CharPtr search, Boolean whole_word)
{
  StringConstraintPtr scp;
  scp = StringConstraintNew();
  scp->match_text = StringSave (search);
  scp->whole_word = whole_word;
  return scp;
}


static SuspectRulePtr MakeSimpleSearchRule (CharPtr search, Boolean whole_word)
{
  SuspectRulePtr rule;

  rule = SuspectRuleNew();
  rule->find = ValNodeNew (NULL);
  rule->find->choice = SearchFunc_string_constraint;
  rule->find->data.ptrvalue = MakeSimpleSearchConstraint (search, whole_word);
  return rule;
}


static SuspectRuleSetPtr MakeSuspectrRNARules (void)
{
  SuspectRuleSetPtr rna_rules = NULL, last_rule = NULL, tmp;
  Int4 i;

  for (i = 0; i < num_suspect_rrna_product_names; i++) {
    tmp = MakeSimpleSearchRule (suspect_rrna_product_names[i], FALSE);
    if (last_rule == NULL) {
      rna_rules = tmp;
    } else {
      last_rule->next = tmp;
    }
    last_rule = tmp;
  }

  tmp = MakeSimpleSearchRule("8S", TRUE);
  tmp->except = ValNodeNew (NULL);
  tmp->except->choice = SearchFunc_string_constraint;
  tmp->except->data.ptrvalue = MakeSimpleSearchConstraint("5.8S", TRUE);
  if (last_rule == NULL) {
    rna_rules = tmp;
  } else {
    last_rule->next = tmp;
  }
  last_rule = tmp;

  return rna_rules;
}


static void FindSuspectrRNAProducts (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  SuspectRuleSetPtr rna_rules;
  ValNodePtr        subcat;
  ClickableItemPtr  cip;

  rna_rules = MakeSuspectrRNARules();

  while (sep_list != NULL) {
    subcat = GetSuspectRuleDiscrepancies (sep_list->data.ptrvalue, rna_rules, FEATDEF_rRNA, DISC_SUSPECT_RRNA_PRODUCTS);
    if (subcat != NULL) {
      cip = SuspectPhraseEx (DISC_SUSPECT_RRNA_PRODUCTS, "suspect phrase", FALSE, "rRNA product name", ItemListFromSubcategories (subcat));
      cip->subcategories = subcat;
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
    sep_list = sep_list->next;
  }
  rna_rules = SuspectRuleSetFree (rna_rules);
}


static CharPtr bad_misc_comment_phrases[] = {
  "catalytic intron"
};

const Int4 num_bad_misc_comment_phrases = sizeof (bad_misc_comment_phrases) / sizeof (CharPtr);

static void FindBadMiscFeaturesCallback (SeqFeatPtr sfp, Pointer data)
{
  ValNodePtr PNTR feature_list;
  Int4 k;

  if (sfp != NULL && sfp->idx.subtype == FEATDEF_misc_feature && (feature_list = (ValNodePtr PNTR)data) != NULL) {
    for (k = 0; k < num_bad_misc_comment_phrases; k++) {
      if (DoesStringContainPhrase (sfp->comment, bad_misc_comment_phrases[k], FALSE, FALSE)) {
        ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
      }
    }
  }
}


static void FindBadMiscFeatures (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  CheckForSuspectPhraseByList (discrepancy_list, sep_list, 
                               bad_misc_comment_phrases, num_bad_misc_comment_phrases,
                               FindBadMiscFeaturesCallback,
                               DISC_SUSPECT_MISC_FEATURES,
                               "misc_feature comment");
}


static Boolean HasParentheses (CharPtr cp)
{
  CharPtr p1;
  Int4    len;

  p1 = StringChr (cp, '(');
  if (p1 == NULL) {
    return FALSE;
  }
  len = StringLen (p1);
  if (*(p1 + len - 1) == ')') {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean HasMissingBacteriaStrain (BioSourcePtr biop)
{
  CharPtr cp;
  OrgModPtr mod;
  Boolean   found = FALSE;

  if (biop == NULL || biop->org == NULL) {
    return FALSE;
  }
  
  cp = StringSearch (biop->org->taxname, " sp. ");
  if (cp == NULL) {
    return FALSE;
  }
  cp += 5;
  if (StringHasNoText (cp) || HasParentheses (cp)) {
    return FALSE;
  }

  if (StringISearch (biop->org->taxname, "enrichment culture clone") != NULL) {
    /* ignore enrichment culture clones */
    return FALSE;
  }

  if (biop->org->orgname == NULL) {
    return FALSE;
  }

  if (!IsBacterialBioSource(biop)) {
    return FALSE;
  }

  for (mod = biop->org->orgname->mod; mod != NULL && !found; mod = mod->next) {
    if (mod->subtype == ORGMOD_strain && StringCmp (mod->subname, cp) == 0) {
      found = TRUE;
    }
  }
  return !found;
}


static void FindMissingBacteriaStrain (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr  item_list = NULL;

  item_list = CollectBioSources (sep_list, HasMissingBacteriaStrain, TRUE);

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_BACTERIA_MISSING_STRAIN, "%d bacterial biosources have taxname 'Genus sp. strain' but no strain", item_list));
  }
}


static Boolean IsBacterialIsolate (BioSourcePtr biop)
{
  OrgModPtr mod;
  Boolean has_bad_isolate = FALSE;

  if (biop == NULL 
      || !IsBacterialBioSource(biop)
      || biop->org == NULL 
      || biop->org->orgname == NULL 
      || biop->org->orgname->mod == NULL
      || HasAmplifiedWithSpeciesSpecificPrimerNote(biop)) {
    return FALSE;
  }
  
  for (mod = biop->org->orgname->mod; mod != NULL && !has_bad_isolate; mod = mod->next) {
    if (mod->subtype == ORGMOD_isolate
        && StringNICmp (mod->subname, "DGGE gel band", 13) != 0
        && StringNICmp (mod->subname, "TGGE gel band", 13) != 0
        && StringNICmp (mod->subname, "SSCP gel band", 13) != 0) {
      has_bad_isolate = TRUE;
    }
  }
  return has_bad_isolate;
}


static void FindBacteriaIsolate (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr  item_list = NULL;

  item_list = CollectBioSources (sep_list, IsBacterialIsolate, TRUE);

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_BACTERIA_SHOULD_NOT_HAVE_ISOLATE, "%d bacterial biosources have isolate", item_list));
  }
}


static void FindMetagenomic (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr  vnp, item_list = NULL, constraint = NULL, src_list, vnp_s;
  SeqEntryPtr sep;
  SourceConstraintPtr src;

  src = SourceConstraintNew ();
  src->field1 = ValNodeNew (NULL);
  src->field1->choice = SourceQualChoice_textqual;
  src->field1->data.intvalue = Source_qual_metagenomic;
  ValNodeAddPointer (&constraint, ConstraintChoice_source, src);

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = (SeqEntryPtr) vnp->data.ptrvalue;
    src_list = GetObjectListForFieldType (FieldType_source_qual, sep);
    for (vnp_s = src_list; vnp_s != NULL; vnp_s = vnp_s->next) {
      if (DoesObjectMatchConstraintChoiceSet (vnp_s->choice, vnp_s->data.ptrvalue, constraint)) {
        ValNodeAddPointer (&item_list, vnp_s->choice, vnp_s->data.ptrvalue);
      }
    }
    src_list = FreeObjectList (src_list);
  }
  constraint = ConstraintChoiceSetFree (constraint);

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_METAGENOMIC, "%d biosources have metagenomic qualifier", item_list));
  }
}


static void FindMetagenomeSource (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr  vnp, item_list = NULL, constraint = NULL, src_list, vnp_s;
  SeqEntryPtr sep;
  SourceConstraintPtr src;

  src = SourceConstraintNew ();
  src->field1 = ValNodeNew (NULL);
  src->field1->choice = SourceQualChoice_textqual;
  src->field1->data.intvalue = Source_qual_metagenome_source;
  ValNodeAddPointer (&constraint, ConstraintChoice_source, src);

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = (SeqEntryPtr) vnp->data.ptrvalue;
    src_list = GetObjectListForFieldType (FieldType_source_qual, sep);
    for (vnp_s = src_list; vnp_s != NULL; vnp_s = vnp_s->next) {
      if (DoesObjectMatchConstraintChoiceSet (vnp_s->choice, vnp_s->data.ptrvalue, constraint)) {
        ValNodeAddPointer (&item_list, vnp_s->choice, vnp_s->data.ptrvalue);
      }
    }
    src_list = FreeObjectList (src_list);
  }
  constraint = ConstraintChoiceSetFree (constraint);

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_METAGENOME_SOURCE, "%d biosources have metagenome_source qualifier", item_list));
  }
}


static void FindBacteriamRNACallback (BioseqPtr bsp, Pointer data)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext context;

  if (bsp == NULL || !BioseqHasLineage(bsp, "Bacteria") || data == NULL) {
    return;
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_mRNA, &context);
  if (sfp != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQ, bsp);
  }
}


static void FindBacteriamRNA (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr  vnp, item_list = NULL;
  SeqEntryPtr sep;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = (SeqEntryPtr) vnp->data.ptrvalue;
    VisitBioseqsInSep (sep, &item_list, FindBacteriamRNACallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_BACTERIA_SHOULD_NOT_HAVE_MRNA, "%d bacterial sequences have mRNA features", item_list));
  }
}

static void FindMissingDefinitionLinesCallback (BioseqPtr bsp, Pointer data)
{
  SeqDescrPtr sdp;
  SeqMgrDescContext context;

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL) {
    return;
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &context);
  if (sdp == NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQ, bsp);
  }
}


static void FindMissingDefinitionLines (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr  vnp, item_list = NULL;
  SeqEntryPtr sep;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = (SeqEntryPtr) vnp->data.ptrvalue;
    VisitBioseqsInSep (sep, &item_list, FindMissingDefinitionLinesCallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_MISSING_DEFLINES, "%d bioseqs have no definition line", item_list));
  }
}


static Boolean IsAffilMissingFromPubdesc (PubdescPtr pdp)
{
  CharPtr str;
  ValNodePtr pub;
  Boolean    rval = FALSE;

  if (pdp == NULL) {
    return FALSE;
  }

  for (pub = pdp->pub; pub != NULL && !rval; pub = pub->next) {
    if (pub->choice == PUB_Sub) {
      str = GetPubFieldFromPub (pub, Publication_field_affiliation, NULL);
      if (StringHasNoText (str)) {
        rval = TRUE;
      }
      str = MemFree (str);
    }
  }
  return rval;
}


static void FindMissingAffiliationsFeatCallback (SeqFeatPtr sfp, Pointer data)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PUB || data == NULL) {
    return;
  }
  if (IsAffilMissingFromPubdesc (sfp->data.value.ptrvalue)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_MAX + OBJ_SEQFEAT, sfp);
  }
}


static void FindMissingAffiliationsDescCallback (SeqDescPtr sdp, Pointer data)
{
  if (sdp == NULL || sdp->choice != Seq_descr_pub || data == NULL) {
    return;
  }
  if (IsAffilMissingFromPubdesc (sdp->data.ptrvalue)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_MAX + OBJ_SEQDESC, sdp);
  }
}


static void FindMissingAffiliations (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr  vnp, item_list = NULL;
  SeqEntryPtr sep;
  SeqSubmitPtr ssp;
  SubmitBlockPtr sbp;
  Boolean        add_object;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = (SeqEntryPtr) vnp->data.ptrvalue;
    add_object = FALSE;
    ssp = FindSeqSubmitForSeqEntry (sep);
    if (ssp != NULL) {
      sbp = ssp->sub;
      if (sbp == NULL) {
        add_object = TRUE;
      } else if (sbp->contact == NULL || sbp->contact->contact == NULL || sbp->contact->contact->affil == NULL
        || StringHasNoText (sbp->contact->contact->affil->affil)) {
        add_object = TRUE;
      } else if (sbp->cit == NULL || sbp->cit->authors == NULL || sbp->cit->authors->affil == NULL
        || StringHasNoText (sbp->cit->authors->affil->affil)) {
        add_object = TRUE;
      }
      if (add_object) {
        ValNodeAddPointer (&item_list, OBJ_SEQSUB, ssp);
      }
    }

    VisitFeaturesInSep (sep, &item_list, FindMissingAffiliationsFeatCallback);
    VisitDescriptorsInSep (sep, &item_list, FindMissingAffiliationsDescCallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_MISSING_AFFIL, "%d citsubs are missing affiliation", item_list));
  }
}


static CharPtr new_exceptions[] = 
{ 
  "annotated by transcript or proteomic data",
  "heterogeneous population sequenced",
  "low-quality sequence region",
  "unextendable partial coding region",
  NULL
};


static void FindCDSNewExceptionCallback (SeqFeatPtr sfp, Pointer data)
{
  Int4 i;

  if (sfp != NULL && data != NULL && sfp->data.choice == SEQFEAT_CDREGION) {
    for (i = 0; new_exceptions[i] != NULL; i++) {
      if (StringISearch (sfp->except_text, new_exceptions[i]) != NULL) {
        ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
        break;
      }
    }
  }
}


static void FindCDSNewException (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr  vnp, item_list = NULL;
  SeqEntryPtr sep;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = (SeqEntryPtr) vnp->data.ptrvalue;
    VisitFeaturesInSep (sep, &item_list, FindCDSNewExceptionCallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_CDS_HAS_NEW_EXCEPTION, "%d coding regions have new exceptions", item_list));
  }
}


typedef struct srcqualkeyword {
  Int4 srcqual;
  CharPtr keyword;
} SrcQualKeywordData, PNTR SrcQualKeywordPtr;

static SrcQualKeywordData srcqual_keywords[] = {
  { Source_qual_forma_specialis, " f. sp." } ,
  { Source_qual_forma, " f." } ,
  { Source_qual_sub_species, " subsp." } ,
  { Source_qual_variety, " var." } ,
  { Source_qual_pathovar, " pv." }
};

#define NUM_srcqual_keywords sizeof (srcqual_keywords) / sizeof (SrcQualKeywordData)


static Boolean IsTrinomialWithoutQualifier (BioSourcePtr biop)
{
  Int4 i, len;
  CharPtr cp, val;
  ValNode vn;
  Boolean rval = FALSE;

  if (biop == NULL || biop->org == NULL || StringHasNoText (biop->org->taxname)) {
    return FALSE;
  }

  /* ignore viruses */
  if (IsViralBioSource(biop)) {
    return FALSE;
  }

  for (i = 0; i < NUM_srcqual_keywords; i++) {
    if ((cp = StringISearch (biop->org->taxname, srcqual_keywords[i].keyword)) != NULL) {
      cp += StringLen (srcqual_keywords[i].keyword);
      while (isspace (*cp)) {
        cp++;
      }
      if (!StringHasNoText (cp)) {
        vn.next = NULL;
        vn.choice = SourceQualChoice_textqual;
        vn.data.intvalue = srcqual_keywords[i].srcqual;
        val = GetSourceQualFromBioSource (biop, &vn, NULL);
        len = StringLen (val);
        if (StringNCmp (cp, val, len) != 0) {
          rval = TRUE;
        }
      }
      break;
    }
  }
  return rval;
}


static void FindTrinomialWithoutQualifier (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr  item_list = NULL;

  item_list = CollectBioSources (sep_list, IsTrinomialWithoutQualifier, TRUE);

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_TRINOMIAL_SHOULD_HAVE_QUALIFIER, "%d trinomial sources lack corresponding qualifier", item_list));
  }
}


static CharPtr rRNATerms[] = {
"16S",
"18S",
"23S",
"26S",
"28S",
"small",
"large",
NULL };

static void FindShortrRNAsCallback (SeqFeatPtr sfp, Pointer data)
{
  ValNodePtr PNTR item_list;
  Int4 i;
  CharPtr    rrna_name;
  Boolean    is_bad = FALSE;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_rRNA || sfp->partial || (item_list = (ValNodePtr PNTR) data) == NULL) {
    return;
  }

  if (SeqLocLen (sfp->location) > 1000) {
    return;
  }

  rrna_name = GetRNAProductString(sfp, NULL);

  for (i = 0; rRNATerms[i] != NULL && !is_bad; i++) {
    if (StringISearch (rrna_name, rRNATerms[i]) != NULL) {
      is_bad = TRUE;
    }
  }

  rrna_name = MemFree (rrna_name);
  if (is_bad) {
    ValNodeAddPointer (item_list, OBJ_SEQFEAT, sfp);
  }
}


static void FindShortrRNAs (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL, vnp;


  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &item_list, FindShortrRNAsCallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_SHORT_RRNA, "%d rRNA features are too short", item_list));
  }
}


static void FindStandardNameCallback (SeqFeatPtr sfp, Pointer data)
{
  GBQualPtr q;

  if (sfp == NULL || data == NULL) {
    return;
  }

  for (q = sfp->qual; q != NULL; q = q->next) {
    if (StringCmp (q->qual, "standard_name") == 0) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
      return;
    }
  }
}


static void FindStandardName (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL, vnp;


  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &item_list, FindStandardNameCallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (ONCALLER_HAS_STANDARD_NAME, "%d features have standard_name qualifier", item_list));
  }
}


static Boolean DoAuthorityAndTaxnameConflict (BioSourcePtr biop)
{
  OrgModPtr mod;
  CharPtr   end;
  size_t    len;

  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL || StringHasNoText (biop->org->taxname)) {
    return FALSE;
  }

  for (mod = biop->org->orgname->mod; mod != NULL && mod->subtype != ORGMOD_authority; mod = mod->next) 
  {}

  if (mod == NULL) {
    return FALSE;
  }

  end = StringChr (biop->org->taxname, ' ');
  if (end != NULL) {
    end = StringChr (end + 1, ' ');
  }

  if (end == NULL) {
    len = StringLen (biop->org->taxname);
  } else {
    len = end - biop->org->taxname;
  }

  if (StringLen (mod->subname) < len || StringNCmp (mod->subname, biop->org->taxname, len) != 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void CheckAuthorityTaxnameConflict (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL;

  item_list = CollectBioSources (sep_list, DoAuthorityAndTaxnameConflict, TRUE);

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (ONCALLER_CHECK_AUTHORITY, "%d biosources have taxname/authority conflict", item_list));
  }
}


NLM_EXTERN AuthListPtr PNTR GetAuthListForPub (PubPtr the_pub)
{
  CitGenPtr  cgp;
  CitSubPtr  csp;
  CitArtPtr  cap;
  CitBookPtr cbp;
  CitPatPtr  cpp;

  if (the_pub == NULL)
  {
    return NULL;
  }
  switch (the_pub->choice) {
    case PUB_Gen :
      cgp = (CitGenPtr) the_pub->data.ptrvalue;
      return &(cgp->authors);
      break;
    case PUB_Sub :
      csp = (CitSubPtr) the_pub->data.ptrvalue;
      return &(csp->authors);
      break;
    case PUB_Article :
      cap = (CitArtPtr) the_pub->data.ptrvalue;
      return &(cap->authors);
      break;
    case PUB_Book :
    case PUB_Man :
      cbp = (CitBookPtr) the_pub->data.ptrvalue;
      return &(cbp->authors);
      break;
    case PUB_Patent :
      cpp = (CitPatPtr) the_pub->data.ptrvalue;
      return &(cpp->authors);
      break;
    default :
      break;
  }
  return NULL;
}


static Boolean AuthorIsConsortium (AuthorPtr ap)
{
  if (ap != NULL && ap->name != NULL && ap->name->choice == 5) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean AuthListHasConsortium (AuthListPtr auth_list)
{
  ValNodePtr  names;
  AuthorPtr   ap;

  if (auth_list == NULL || auth_list->choice != 1) {
    return FALSE;
  }
  for (names = auth_list->names; names != NULL; names = names->next) {
    ap = names->data.ptrvalue;
    if (AuthorIsConsortium(ap)) {
      return TRUE;
    }
  }
  return FALSE;
}


static Boolean PubEquivHasConsortium (ValNodePtr pub)
{
  ValNodePtr vnp;
  AuthListPtr PNTR p_auth;

  for (vnp = pub; vnp != NULL; vnp = vnp->next) {
    p_auth = GetAuthListForPub (vnp);
    if (p_auth != NULL && *p_auth != NULL && AuthListHasConsortium(*p_auth)) {
      return TRUE;
    }
  }
  return FALSE;
}


static Boolean PubHasConsortium (PubdescPtr pdp)
{
  if (pdp == NULL) {
    return FALSE;
  } else {
    return PubEquivHasConsortium (pdp->pub);
  }
}


static Boolean ContactInfoHasConsortium (ContactInfoPtr contact_info)
{
  if (contact_info != NULL && AuthorIsConsortium(contact_info->contact)) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void FindConsortiumsDescCallback (SeqDescrPtr sdp, Pointer data)
{
  if (sdp == NULL || sdp->choice != Seq_descr_pub || data == NULL) {
    return;
  }
  if (PubHasConsortium (sdp->data.ptrvalue)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}


static void FindConsortiumsFeatCallback (SeqFeatPtr sfp, Pointer data)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PUB || data == NULL) {
    return;
  }
  if (PubHasConsortium (sfp->data.value.ptrvalue)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
}


static void FindConsortiums (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL, vnp;
  SeqSubmitPtr ssp;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitDescriptorsInSep (vnp->data.ptrvalue, &item_list, FindConsortiumsDescCallback);
  }
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &item_list, FindConsortiumsFeatCallback);
  }
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    ssp = FindSeqSubmitForSeqEntry (vnp->data.ptrvalue);
    if (ssp != NULL && ssp->sub != NULL 
        && (ContactInfoHasConsortium (ssp->sub->contact) 
            || (ssp->sub->cit != NULL && AuthListHasConsortium(ssp->sub->cit->authors)))) {
      ValNodeAddPointer (&item_list, OBJ_SEQSUB, ssp);
    }
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (ONCALLER_CONSORTIUM, "%d publications/submitter blocks have consortium", item_list));
  }
}


static void RemoveConsortiumFromAuthList (AuthListPtr alp)
{
  ValNodePtr  names, prev = NULL, names_next;
  AuthorPtr   ap;

  if (alp == NULL) {
    return;
  }

  for (names = alp->names; names != NULL; names = names_next) 
  {
    names_next = names->next;
    ap = names->data.ptrvalue;
    if (ap->name->choice == 5)
    {
      ap->name->data = MemFree (ap->name->data);
      AuthorFree (ap);
      if (prev == NULL)
      {
        alp->names = names->next;
      }
      else
      {
        prev->next = names->next;
      }
      names->next = NULL;
      names = ValNodeFree (names);
    }
    else
    {
      prev = names;
    }
  }  
}


NLM_EXTERN void RemoveConsortiumFromPub (PubPtr pub)
{
  AuthListPtr PNTR p_auth_list;
  AuthListPtr alp;
  
  p_auth_list = GetAuthListForPub (pub);
  if (p_auth_list == NULL || (alp = *p_auth_list) == NULL) {
    return;
  }
  RemoveConsortiumFromAuthList (alp);  
}


static void RemoveConsortiumFromPubdesc (PubdescPtr pdp)
{
  ValNodePtr pub;

  if (pdp != NULL) {
    for (pub = pdp->pub; pub != NULL; pub = pub->next) {
      RemoveConsortiumFromPub (pub);
    }
  }
}


NLM_EXTERN void RemoveConsortiums (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  SeqDescPtr sdp;
  PubdescPtr pdp;
  SeqSubmitPtr ssp;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case OBJ_SEQFEAT:
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PUB) {
          pdp = sfp->data.value.ptrvalue;
          RemoveConsortiumFromPubdesc (pdp);
        }
        break;
      case OBJ_SEQDESC:
        sdp = (SeqDescPtr) vnp->data.ptrvalue;
        if (sdp != NULL && sdp->choice == Seq_descr_pub) {
          pdp = sdp->data.ptrvalue;
          RemoveConsortiumFromPubdesc (pdp);
        }
        break;
      case OBJ_SEQSUB:
        ssp = (SeqSubmitPtr) vnp->data.ptrvalue;
        if (ssp != NULL && ssp->sub != NULL) {
          if (ssp->sub->contact != NULL && AuthorIsConsortium (ssp->sub->contact->contact)) {
            ssp->sub->contact->contact = AuthorFree (ssp->sub->contact->contact);
          }
          if (ssp->sub->cit != NULL) {
            RemoveConsortiumFromAuthList(ssp->sub->cit->authors);
          }
        }
        break;
    }
  }
}


static Boolean MatchExceptSpaceColon (CharPtr str1, CharPtr str2)
{
  if (str1 == NULL && str2 == NULL) {
    return TRUE;
  }
  while ((str1 == NULL || *str1 != 0) && (str2 == NULL || *str2 != 0)) {
    if (str1 != NULL && (*str1 == ':' || isspace (*str1))) {
      str1++;
    } else if (str2 != NULL && (*str2 == ':' || isspace (*str2))) {
      str2++;
    } else if (str1 != NULL && str2 != NULL && *str1 != *str2) {
      return FALSE;
    } else if (str1 == NULL && *str2 != 0) {
      return FALSE;
    } else if (str2 == NULL && *str1 != 0) {
      return FALSE;
    } else {
      if (str1 != NULL && *str1 != 0) {
        str1++;
      }
      if (str2 != NULL && *str2 != 0) {
        str2++;
      }
    } 
  }
  if ((str1 != NULL && *str1 != 0) || (str2 != NULL && *str2 != 0)) {
    return FALSE;
  } else {
    return TRUE;
  }
}


static Boolean BioSourceHasConflictingStrainAndCultureCollectionValues (BioSourcePtr biop)
{
  OrgModPtr strain, culture;
  Boolean   has_conflict = FALSE, has_match = FALSE;
  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL || biop->org->orgname->mod == NULL) {
    return FALSE;
  }

  for (strain = biop->org->orgname->mod; strain != NULL && !has_match; strain = strain->next) {
    if (strain->subtype == ORGMOD_strain) {
      for (culture = biop->org->orgname->mod; culture != NULL && !has_match; culture = culture->next) {
        if (culture->subtype == ORGMOD_culture_collection) {
          if (MatchExceptSpaceColon(strain->subname, culture->subname)) {
            has_match = TRUE;
          } else {
            has_conflict = TRUE;
          }
        }
      }
    }
  }
  if (has_conflict && !has_match) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void FindStrainCultureCollectionMismatch (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL, extended_item_list = NULL;
  ValNodePtr field_list = NULL, field, src_qual;

  item_list = CollectBioSources (sep_list, BioSourceHasConflictingStrainAndCultureCollectionValues, TRUE);

  if (item_list != NULL) {
    src_qual = ValNodeNew (NULL);
    src_qual->choice = SourceQualChoice_textqual;
    src_qual->data.intvalue = Source_qual_strain;
    field = ValNodeNew (NULL);
    field->choice = FieldType_source_qual;
    field->data.ptrvalue = src_qual;
    field_list = field;
    src_qual = ValNodeNew (NULL);
    src_qual->choice = SourceQualChoice_textqual;
    src_qual->data.intvalue = Source_qual_culture_collection;
    field = ValNodeNew (field_list);
    field->choice = FieldType_source_qual;
    field->data.ptrvalue = src_qual;
    
    extended_item_list = MakeObjectListWithFields (item_list, field_list);
    item_list = ValNodeFree (item_list);
    field_list = FieldTypeListFree (field_list);
    ValNodeAddPointer (discrepancy_list, 0, 
                        NewClickableItem (ONCALLER_STRAIN_CULTURE_COLLECTION_MISMATCH, 
                        "%d organisms have conflicting strain and culture-collection values", 
                        extended_item_list));
  }
}


static Boolean HasMultiSrc (BioSourcePtr biop)
{
  OrgModPtr mod;

  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) {
    return FALSE;
  }

  for (mod = biop->org->orgname->mod; mod != NULL; mod = mod->next) {
    if ((mod->subtype == ORGMOD_strain || mod->subtype == ORGMOD_isolate)
      && (StringChr (mod->subname, ',') != NULL || StringChr (mod->subname, ';') != NULL)) {
      return TRUE;
    }
  }
  return FALSE;
}


static void FindMultiSrc (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL;

  item_list = CollectBioSources (sep_list, HasMultiSrc, TRUE);

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (ONCALLER_MULTISRC, "%d organisms have comma or semicolon in strain or isolate", item_list));
  }
}

static Boolean HasMultipleCultureCollection (BioSourcePtr biop)
{
  OrgModPtr mod;
  Boolean   has_one = FALSE;

  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) {
    return FALSE;
  }

  for (mod = biop->org->orgname->mod; mod != NULL; mod = mod->next) {
    if (mod->subtype == ORGMOD_culture_collection) {
      if (has_one) {
        return TRUE;
      } else {
        has_one = TRUE;
      }
    }
  }
  return FALSE;
}


static void FindMultipleCultureCollection (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL, extended_item_list = NULL, src_qual, field_list;

  item_list = CollectBioSources (sep_list, HasMultipleCultureCollection, TRUE);

  if (item_list != NULL) {
    src_qual = ValNodeNew (NULL);
    src_qual->choice = SourceQualChoice_textqual;
    src_qual->data.intvalue = Source_qual_culture_collection;
    field_list = ValNodeNew (NULL);
    field_list->choice = FieldType_source_qual;
    field_list->data.ptrvalue = src_qual;
    
    extended_item_list = MakeObjectListWithFields (item_list, field_list);
    field_list = FieldTypeListFree (field_list);
    item_list = ValNodeFree (item_list);

    ValNodeAddPointer (discrepancy_list, 0, 
                       NewClickableItem (ONCALLER_MULTIPLE_CULTURE_COLLECTION, 
                                         "%d organisms have multiple culture-collection qualifiers", 
                                         extended_item_list));
  }
}


static Boolean HasHumanHost(BioSourcePtr biop)
{
  OrgModPtr mod;

  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) {
    return FALSE;
  }

  for (mod = biop->org->orgname->mod; mod != NULL; mod = mod->next) {
    if (mod->subtype == ORGMOD_nat_host && DoesStringContainPhrase (mod->subname, "human", FALSE, TRUE)) {
      return TRUE;
    }
  }
  return FALSE;
}


static void FindHumanHosts (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL;

  item_list = CollectBioSources (sep_list, HasHumanHost, TRUE);

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_HUMAN_HOST, "%d organisms have 'human' host qualifiers", item_list));
  }
}


NLM_EXTERN void FixHumanHosts (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr   vnp, entityIDList = NULL;
  Uint2        entityID;
  BioSourcePtr biop;
  OrgModPtr    mod;
  Int4         num_changed = 0;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    biop = GetBioSourceFromObject (vnp->choice, vnp->data.ptrvalue);
    if (biop != NULL && biop->org != NULL && biop->org->orgname != NULL) {      
      for (mod = biop->org->orgname->mod; mod != NULL; mod = mod->next) {
        if (mod->subtype == ORGMOD_nat_host && DoesStringContainPhrase (mod->subname, "human", FALSE, TRUE)) {
          FindReplaceString (&(mod->subname), "human", "Homo sapiens", FALSE, TRUE);
          num_changed++;
          entityID = GetEntityIdFromObject (vnp->choice, vnp->data.ptrvalue);
          ValNodeAddInt (&entityIDList, 0, entityID);
        }
      }
    }
  }

  entityIDList = ValNodeSort (entityIDList, SortByIntvalue);
  ValNodeUnique (&entityIDList, SortByIntvalue, ValNodeFree);
  
  for (vnp = entityIDList; vnp != NULL; vnp = vnp->next) {
    ObjMgrSetDirtyFlag (vnp->data.intvalue, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, vnp->data.intvalue, 0, 0);    
  }

  if (num_changed > 0 && lip != NULL && lip->fp != NULL) {
    fprintf (lip->fp, "Changed %d host qualifiers from 'human' to 'Homo sapiens'\n", num_changed);
    lip->data_in_log = TRUE;
  }
}


static void FindSegSetsCallback (BioseqSetPtr bssp, Pointer data)
{
  if (bssp != NULL && data != NULL && bssp->_class == BioseqseqSet_class_segset) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQSET, bssp);
  }
}


static void FindSegSets (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL, vnp;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitSetsInSep (vnp->data.ptrvalue, &item_list, FindSegSetsCallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_SEGSETS_PRESENT, "%d segsets are present", item_list));
  }
}


static void FindNonWGSSetsCallback (BioseqSetPtr bssp, Pointer data)
{
  if (bssp != NULL && data != NULL
      && (bssp->_class == BioseqseqSet_class_eco_set
          || bssp->_class == BioseqseqSet_class_mut_set
          || bssp->_class == BioseqseqSet_class_phy_set
          || bssp->_class == BioseqseqSet_class_pop_set)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQSET, bssp);
  }
}


static void FindNonWGSSets (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr item_list = NULL, vnp;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitSetsInSep (vnp->data.ptrvalue, &item_list, FindNonWGSSetsCallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_NONWGS_SETS_PRESENT, "%d sets are of type eco, mut, phy or pop", item_list));
  }
}


static void ListAllFeatures (SeqFeatPtr sfp, Pointer data)
{
  ValNodePtr PNTR feature_list;
  ValNodePtr vnp, vnp_last = NULL, item_list;

  if (sfp != NULL 
      && sfp->idx.subtype != FEATDEF_gap
      && sfp->idx.subtype != FEATDEF_PROT
      && (feature_list = (ValNodePtr PNTR)data) != NULL) {
    for (vnp = *feature_list; vnp != NULL && vnp->choice != sfp->idx.subtype; vnp = vnp->next) {
      vnp_last = vnp;
    }
    if (vnp == NULL) {
      vnp = ValNodeNew (vnp_last);
      vnp->choice = sfp->idx.subtype;
      if (vnp_last == NULL) {
        *feature_list = vnp;
      }
    }
    item_list = vnp->data.ptrvalue;
    ValNodeAddPointer (&item_list, OBJ_SEQFEAT, sfp);
    vnp->data.ptrvalue = item_list;
  }
}


static void GetFeatureList  (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr feature_list = NULL, vnp, item_list;
  ClickableItemPtr cip;
  CharPtr fmt_fmt = "%%d %s features";
  CharPtr fmt;
  FeatDefPtr  curr;
  Uint1       key;
  CharPtr     label = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &feature_list, ListAllFeatures);
  }

  for (vnp = feature_list; vnp != NULL; vnp = vnp->next) {
    item_list = vnp->data.ptrvalue;
    label = NULL;
    curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
    while (curr != NULL && curr->featdef_key != vnp->choice) {
      curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
    }
    if (curr == NULL) {
      label = "unknown";
    } else {
      label = curr->typelabel;
    }
    fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt_fmt) + StringLen (label)));
    sprintf (fmt, fmt_fmt, label);
    cip = NewClickableItem (DISC_FEATURE_LIST, fmt, item_list);
    fmt = MemFree (fmt);
    vnp->choice = 0;
    vnp->data.ptrvalue = cip;
  }

  if (feature_list != NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->clickable_item_type = DISC_FEATURE_LIST;
    cip->description = StringSave ("Feature List");
    cip->item_list = NULL;
    cip->subcategories = feature_list;

    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}

static void FindBadBacterialGeneNamesCallback (BioseqPtr bsp, Pointer data)
{
  SeqDescrPtr sdp;
  SeqMgrDescContext dcontext;
  BioSourcePtr biop;
  SeqFeatPtr gene;
  SeqMgrFeatContext fcontext;
  GeneRefPtr grp;

  if (bsp == NULL || data == NULL) {
    return;
  }
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL || (biop = sdp->data.ptrvalue) == NULL || biop->org == NULL
      || biop->org->orgname == NULL
      || !IsBacterialBioSource (biop)) {
    return;
  }

  for (gene = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &fcontext);
       gene != NULL;
       gene = SeqMgrGetNextFeature (bsp, gene, SEQFEAT_GENE, 0, &fcontext)) {
    grp = (GeneRefPtr) gene->data.value.ptrvalue;
    if (grp != NULL && grp->locus != NULL && (!isalpha (*(grp->locus)) || !islower (*(grp->locus)))) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, gene);
    }
  }
}


typedef Boolean (*BadGeneNameTestFunc) PROTO ((CharPtr, CharPtr, SeqFeatPtr));

typedef struct badgenename {
  CharPtr pattern;
  BadGeneNameTestFunc func;
} BadGeneNameData, PNTR BadGeneNamePtr;

static Boolean GeneNameLongerThanTenChars (CharPtr pattern, CharPtr search, SeqFeatPtr sfp)
{
  if (StringLen (search) > 10) {
    return TRUE;
  } else {
    return FALSE;
  }
}

static Boolean GeneNameContainsPhrase (CharPtr pattern, CharPtr search, SeqFeatPtr sfp)
{
  if (StringISearch (search, pattern) != NULL) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean GeneNameHas4Numbers (CharPtr pattern, CharPtr search, SeqFeatPtr sfp)
{
  CharPtr cp;
  Int4    num_digits = 0;

  if (search == NULL) {
    return FALSE;
  }

  for (cp = search; *cp != 0 && num_digits < 4; cp++) {
    if (isdigit (*cp)) {
      ++num_digits;
    } else {
      num_digits = 0;
    }
  }
  if (num_digits >= 4) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static BadGeneNameData bad_gene_rules[] = {
  { "more than 10 characters", GeneNameLongerThanTenChars },
  { "putative", GeneNameContainsPhrase },
  { "fragment", GeneNameContainsPhrase },
  { "gene", GeneNameContainsPhrase },
  { "orf", GeneNameContainsPhrase },
  { "like", GeneNameContainsPhrase },
  { "4 or more consecutive numbers", GeneNameHas4Numbers }
};


static const Int4 kNumBadGeneRules = sizeof (bad_gene_rules) / sizeof (BadGeneNameData);

static void FindBadGeneNameCallback (SeqFeatPtr sfp, Pointer data)
{
  ValNodePtr PNTR feature_lists;
  GeneRefPtr grp;
  Int4 k;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE 
      || (grp = (GeneRefPtr) sfp->data.value.ptrvalue) == NULL
      || StringHasNoText (grp->locus)
      || (feature_lists = (ValNodePtr PNTR) data) == NULL) {
    return;
  }

  for (k = 0; k < kNumBadGeneRules; k++) {
    if (bad_gene_rules[k].func(bad_gene_rules[k].pattern, grp->locus, sfp)) {
      ValNodeAddPointer (feature_lists + k, OBJ_SEQFEAT, sfp);
    }
  }
}

  
static void FindBadGeneNames (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr PNTR feature_lists, vnp;
  ValNodePtr bad_bacterial_genes = NULL;
  ValNodePtr subcat = NULL;
  CharPtr fmt = "%d bacterial genes do not start with lowercase letters";
  Int4 k;
  ClickableItemPtr dip;

  feature_lists = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * kNumBadGeneRules);
  MemSet (feature_lists, 0, sizeof (ValNodePtr) * kNumBadGeneRules);
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, feature_lists, FindBadGeneNameCallback);
    VisitBioseqsInSep (vnp->data.ptrvalue, &bad_bacterial_genes, FindBadBacterialGeneNamesCallback);
  }

  if (bad_bacterial_genes != NULL) {
    ValNodeAddPointer (&subcat, 0, NewClickableItem (DISC_BAD_BACTERIAL_GENE_NAME, fmt, bad_bacterial_genes));
  }

  for (k = 0; k < kNumBadGeneRules; k++) {
    if (feature_lists[k] != NULL) {
      ValNodeAddPointer (&subcat, 0, SuspectPhraseEx(TEST_BAD_GENE_NAME, bad_gene_rules[k].pattern, FALSE, "gene", feature_lists[k]));
    }
  }
  feature_lists = MemFree (feature_lists);

  if (subcat == NULL) {
    /* do nothing */
  } else if (subcat->next == NULL) {
    ValNodeLink (discrepancy_list, subcat);
  } else {
    dip = SuspectPhraseEx (TEST_BAD_GENE_NAME, "suspect phrase or characters", FALSE, "gene", ItemListFromSubcategories (subcat));
    if (dip != NULL)
    {
      dip->subcategories = subcat;      
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}


static void MoveBadGeneNames (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  SeqFeatPtr sfp;
  GeneRefPtr grp;
  ValNodePtr vnp;
  Int4 num = 0;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT && (sfp = (SeqFeatPtr) vnp->data.ptrvalue) != NULL
        && sfp->data.choice == SEQFEAT_GENE
        && (grp = (GeneRefPtr) sfp->data.value.ptrvalue) != NULL
        && !StringHasNoText (grp->locus)) {
      SetStringValue (&(sfp->comment), grp->locus, ExistingTextOption_append_semi);
      grp->locus = MemFree (grp->locus);
      num++;
    }
  }
  if (num > 0 && lip != NULL) {
    lip->data_in_log = TRUE;
    if (lip->fp != NULL) {
      fprintf (lip->fp, "Moved %d bad gene names to gene comment.\n", num);
    }
  }
}


typedef struct bioseqsetclassnameclassval {
  CharPtr class_name;
  Uint1   class_val;
} BioseqSetClassNameClassValData, PNTR BioseqSetClassNameClassValPtr;

static BioseqSetClassNameClassValData bioseqsetclassname_classval[] = {
  {"  ",                     BioseqseqSet_class_not_set},
  {"Nuc-prot",               BioseqseqSet_class_nuc_prot},
  {"Segset",                 BioseqseqSet_class_segset},
  {"Conset",                 BioseqseqSet_class_conset},
  {"Parts",                  BioseqseqSet_class_parts},
  {"Gibb",                   BioseqseqSet_class_gibb},
  {"GI",                     BioseqseqSet_class_gi},
  {"Genbank",                BioseqseqSet_class_genbank},
  {"PIR",                    BioseqseqSet_class_pir},
  {"Pubset",                 BioseqseqSet_class_pub_set},
  {"Equiv",                  BioseqseqSet_class_equiv},
  {"Swissprot",              BioseqseqSet_class_swissprot},
  {"PDB-entry",              BioseqseqSet_class_pdb_entry},
  {"Mut-set",                BioseqseqSet_class_mut_set},
  {"Pop-set",                BioseqseqSet_class_pop_set},
  {"Phy-set",                BioseqseqSet_class_phy_set},
  {"Eco-set",                BioseqseqSet_class_eco_set},
  {"Gen-prod-set",           BioseqseqSet_class_gen_prod_set},
  {"WGS-set",                BioseqseqSet_class_wgs_set},
  {"Small-genome-set",       BioseqseqSet_class_small_genome_set},
  {"Other",                  BioseqseqSet_class_other}};

#define NUM_bioseqsetclassname_classval sizeof (bioseqsetclassname_classval) / sizeof (BioseqSetClassNameClassValData)


NLM_EXTERN CharPtr GetSetClassName (Uint1 class_val)
{
  Int4 i;

  for (i = 0; i < NUM_bioseqsetclassname_classval; i++) {
    if (bioseqsetclassname_classval[i].class_val == class_val) {
      return bioseqsetclassname_classval[i].class_name;
    }
  }
  return NULL;
}


static void ReorganizeDoubleGenBankSets (BioseqSetPtr parent_bssp) 
{
  SeqEntryPtr this_sep;
  BioseqSetPtr target_bssp;
  SeqAnnotPtr  sap_last;

  if (parent_bssp == NULL) {
    return;
  }

  for (this_sep = parent_bssp->seq_set; this_sep != NULL; this_sep = this_sep->next) {
    if (IS_Bioseq_set (this_sep)) {
      target_bssp = this_sep->data.ptrvalue;
      if (target_bssp == NULL) {
        continue;
      }
      if (parent_bssp->_class == BioseqseqSet_class_genbank 
          && target_bssp->_class == BioseqseqSet_class_genbank) {
        ValNodeLink (&(parent_bssp->seq_set), target_bssp->seq_set);
        target_bssp->seq_set = NULL;
        ValNodeLink (&(parent_bssp->descr), target_bssp->descr);          
        target_bssp->descr = NULL;
        sap_last = parent_bssp->annot;
        while (sap_last != NULL && sap_last->next != NULL) {
          sap_last = sap_last->next;
        }
        if (sap_last == NULL) {
          parent_bssp->annot = target_bssp->annot;
        } else {
          sap_last->next = target_bssp->annot;
        }
        target_bssp->annot = NULL;
        target_bssp->idx.deleteme = TRUE;
      } else {
        ReorganizeDoubleGenBankSets (target_bssp);
      }
    }
  }
}


static void RemoveDoubleGenBankSets (BioseqSetPtr parent_bssp, Uint2 entityID)
{
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  Uint2             top_parenttype;
  Pointer           top_parentptr;
  SeqEntryPtr       top_sep;

  if (parent_bssp == NULL) {
    return;
  }

  top_sep = GetTopSeqEntryForEntityID (entityID);
  if (top_sep == NULL) return;
  SaveSeqEntryObjMgrData (top_sep, &omdptop, &omdata);
  GetSeqEntryParent (top_sep, &top_parentptr, &top_parenttype);

  ReorganizeDoubleGenBankSets (parent_bssp);

  SeqMgrLinkSeqEntry (top_sep, top_parenttype, top_parentptr);
    
  SeqMgrClearFeatureIndexes (entityID, NULL);
  SeqMgrIndexFeatures (entityID, NULL);

  RestoreSeqEntryObjMgrData (top_sep, omdptop, &omdata);
  
  SeqMgrClearFeatureIndexes (entityID, NULL);
  SeqMgrIndexFeatures (entityID, NULL);

  DeleteMarkedObjects (entityID, 0, NULL);
  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);    

}


NLM_EXTERN void FixNonWGSSets (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr vnp;
  BioseqSetPtr  bssp;
  CharPtr       class_name;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_BIOSEQSET) {
      bssp = vnp->data.ptrvalue;
      if (bssp != NULL) {
        class_name = GetSetClassName (bssp->_class);

        bssp->_class = BioseqseqSet_class_genbank;
        if (bssp->idx.parenttype == OBJ_BIOSEQSET) {
          RemoveDoubleGenBankSets (bssp->idx.parentptr, bssp->idx.entityID);
        }
        if (lip != NULL && lip->fp != NULL) {
          fprintf (lip->fp, "Bioseq-set class changed from %s to genbank\n", class_name == NULL ? "unknown" : class_name);
          lip->data_in_log = TRUE;
        }
      }
    }
  }
}


static void FindMismatchedCommentsCallback (SeqDescrPtr sdp, Pointer data)
{
  if (sdp != NULL && sdp->choice == Seq_descr_comment && data != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}


NLM_EXTERN void FindMismatchedComments (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr comment_list = NULL, cat_list = NULL, vnp_prev = NULL, vnp;
  CharPtr    curr_val = NULL;
  SeqDescrPtr sdp;
  ClickableItemPtr cip = NULL;
  CharPtr sub_fmt = "%%d comments contain %s";
  CharPtr fmt;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitDescriptorsInSep (vnp->data.ptrvalue, &comment_list, FindMismatchedCommentsCallback);
  }

  if (comment_list == NULL) {
    return;
  }

  comment_list = ValNodeSort (comment_list, SortVnpByObject);
  for (vnp = comment_list; vnp != NULL; vnp = vnp->next) {
    sdp = (SeqDescrPtr) vnp->data.ptrvalue;
    if (curr_val != NULL && StringCmp (curr_val, sdp->data.ptrvalue) != 0) {
      if (vnp_prev != NULL) {
        vnp_prev->next = NULL;
      }
      fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (sub_fmt) + StringLen ((CharPtr) sdp->data.ptrvalue)));
      sprintf (fmt, sub_fmt, (CharPtr) sdp->data.ptrvalue);
      cip = NewClickableItem (DISC_MISMATCHED_COMMENTS, fmt, comment_list);
      ValNodeAddPointer (&cat_list, 0, cip);
      comment_list = vnp;
    }
    curr_val = (CharPtr) sdp->data.ptrvalue;
    vnp_prev = vnp;
  }
  if (cat_list != NULL) {
    fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (sub_fmt) + StringLen ((CharPtr) sdp->data.ptrvalue)));
    sprintf (fmt, sub_fmt, (CharPtr) sdp->data.ptrvalue);
    cip = NewClickableItem (DISC_MISMATCHED_COMMENTS, fmt, comment_list);
    ValNodeAddPointer (&cat_list, 0, cip);
    comment_list = NULL;
  }

  if (cat_list != NULL && cat_list->next != NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->clickable_item_type = DISC_MISMATCHED_COMMENTS;
    cip->description = StringSave ("Mismatched comments were found");
    cip->subcategories = cat_list;
    cat_list = NULL;
    cip->item_list = ItemListFromSubcategories (cip->subcategories);
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }

  comment_list = ValNodeFree (comment_list);
  cat_list = FreeClickableList (cat_list);
}


static void FixMismatchedCommentsCallback (SeqDescrPtr sdp, Pointer data)
{
  if (sdp != NULL && sdp->choice == Seq_descr_comment && data != NULL) {
    sdp->data.ptrvalue = MemFree (sdp->data.ptrvalue);
    sdp->data.ptrvalue = StringSave ((CharPtr) data);
  }
}


NLM_EXTERN void FixMismatchedComments (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr entityIDList = NULL, vnp;
  SeqDescrPtr sdp;
  CharPtr     new_val;
  ObjValNodePtr ovp;
  SeqEntryPtr   sep;

  if (item_list == NULL) {
    return;
  }

  sdp = (SeqDescrPtr) item_list->data.ptrvalue;
  if (sdp == NULL || sdp->data.ptrvalue == NULL) {
    return;
  }

  new_val = StringSave (sdp->data.ptrvalue);

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    sdp = (SeqDescrPtr) vnp->data.ptrvalue;
    if (sdp->extended) {
      ovp = (ObjValNodePtr) sdp;
      ValNodeAddInt (&entityIDList, 0, ovp->idx.entityID);
    }
  }

  entityIDList = ValNodeSort (entityIDList, SortByIntvalue);
  ValNodeUnique (&entityIDList, SortByIntvalue, ValNodeFree);
  
  for (vnp = entityIDList; vnp != NULL; vnp = vnp->next) {
    sep = GetTopSeqEntryForEntityID (vnp->data.intvalue);
    VisitDescriptorsInSep (sep, new_val, FixMismatchedCommentsCallback);
    ObjMgrSetDirtyFlag (vnp->data.intvalue, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, vnp->data.intvalue, 0, 0);    
  }
  if (lip != NULL && lip->fp != NULL) {
    fprintf (lip->fp, "Replaced all coments with '%s'\n", new_val);
    lip->data_in_log = TRUE;
  }
  new_val = MemFree (new_val);
}


static void FindOrderedLocationsCallback (SeqFeatPtr sfp, Pointer data)
{
  if (sfp != NULL && data != NULL && LocationHasNullsBetween(sfp->location)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
}


NLM_EXTERN void FindOrderedLocations (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp;
  ValNodePtr feat_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &feat_list, FindOrderedLocationsCallback);
  }

  if (feat_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, 
                       NewClickableItem (ONCALLER_ORDERED_LOCATION, "%d features have ordered locations", feat_list));
  }
}


static void FixOrderedLocationsCallback (SeqFeatPtr sfp, Pointer data)
{
  SeqLocPtr slp_prev = NULL, slp, slp_next;
  Boolean changed_loc = FALSE;
  CharPtr orig, repl;
  LogInfoPtr lip;

  if (sfp != NULL && sfp->location != NULL && sfp->location->choice == SEQLOC_MIX) {
    orig = SeqLocPrint (sfp->location);
    for (slp = sfp->location->data.ptrvalue; slp != NULL; slp = slp_next) {
      slp_next = slp->next;
      if (slp->choice == SEQLOC_NULL) {
        if (slp_prev == NULL) {
          sfp->location->data.ptrvalue = slp_next;
        } else {
          slp_prev->next = slp_next;
        }
        slp->next = NULL;
        slp = SeqLocFree (slp);
        changed_loc = TRUE;
      } else {
        slp_prev = slp;
      }
    }
    if (changed_loc && (lip = (LogInfoPtr) data) != NULL && lip->fp != NULL) {
      repl = SeqLocPrint (sfp->location);
      fprintf (lip->fp, "Changed location from %s to %s", orig, repl);
      repl = MemFree (repl);
      lip->data_in_log = TRUE;
    }
    orig = MemFree (orig);
  }
}


NLM_EXTERN void FixOrderedLocations (ValNodePtr item_list, Pointer data, LogInfoPtr lip)
{
  ValNodePtr entityIDList = NULL, vnp;
  SeqFeatPtr sfp;

  if (item_list == NULL) {
    return;
  }

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    ValNodeAddInt (&entityIDList, 0, sfp->idx.entityID);
    FixOrderedLocationsCallback (vnp->data.ptrvalue, lip);
  }

  entityIDList = ValNodeSort (entityIDList, SortByIntvalue);
  ValNodeUnique (&entityIDList, SortByIntvalue, ValNodeFree);
  
  for (vnp = entityIDList; vnp != NULL; vnp = vnp->next) {
    ObjMgrSetDirtyFlag (vnp->data.intvalue, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, vnp->data.intvalue, 0, 0);    
  }
}


static void FindCommentDescriptorsCallback(SeqDescrPtr sdp, Pointer data)
{
  if (sdp != NULL && sdp->choice == Seq_descr_comment && data != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}


NLM_EXTERN void FindCommentDescriptors (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp;
  ValNodePtr desc_list = NULL;
  Boolean    all_same = TRUE;
  SeqDescPtr sdp1, sdp2;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitDescriptorsInSep (vnp->data.ptrvalue, &desc_list, FindCommentDescriptorsCallback);
  }

  if (desc_list != NULL) {
    sdp1 = desc_list->data.ptrvalue;
    vnp = desc_list->next;
    while (vnp != NULL && all_same) {
      sdp2 = vnp->data.ptrvalue;
      if (StringCmp (sdp1->data.ptrvalue, sdp2->data.ptrvalue) != 0) {
        all_same = FALSE;
      }
      vnp = vnp->next;
    }
    ValNodeAddPointer (discrepancy_list, 0, 
                       NewClickableItem (ONCALLER_COMMENT_PRESENT, 
                       all_same ? "%d comment descriptors were found (all same)" : "%d comment descriptors were found (some different)",
                       desc_list));
  }
}


static void FindTitlesOnSetsCallback (BioseqSetPtr bssp, Pointer data)
{
  SeqDescPtr sdp;
  ClickableItemPtr cip;

  if (bssp == NULL || data == NULL) {
    return;
  }
  for (sdp = bssp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_title) {
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      cip->description = StringSave ((CharPtr) sdp->data.ptrvalue);
      cip->clickable_item_type = ONCALLER_DEFLINE_ON_SET;
      ValNodeAddPointer (&(cip->item_list), OBJ_SEQDESC, sdp);
      ValNodeAddPointer ((ValNodePtr PNTR) data, 0, cip);
    }
  }
}


NLM_EXTERN void FindTitlesOnSets (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp;
  ValNodePtr title_list = NULL, item_list;
  ClickableItemPtr cip;
  Char             tmp[30];

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitSetsInSep (vnp->data.ptrvalue, &title_list, FindTitlesOnSetsCallback);
  }

  if (title_list != NULL) {
    item_list = ItemListFromSubcategories(title_list);
    cip = NewClickableItem (ONCALLER_DEFLINE_ON_SET, "%d titles on sets were found", item_list);
    cip->subcategories = title_list;
    if (GetAppParam ("SEQUINCUSTOM", "ONCALLERTOOL", "EXPAND_DEFLINE_ON_SET", NULL, tmp, sizeof (tmp) - 1)
        && StringICmp (tmp, "TRUE") == 0) {
      cip->expanded = TRUE;
    }
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


static void FindInconsistentHIVRNACallback (BioseqPtr bsp, Pointer data)
{
  SeqDescrPtr bsdp, msdp;
  SeqMgrDescContext context;
  BioSourcePtr biop;
  MolInfoPtr   mip;

  if (bsp == NULL || data == NULL || bsp->mol != Seq_mol_rna) {
    return;
  }

  bsdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (bsdp == NULL 
      || (biop = (BioSourcePtr) bsdp->data.ptrvalue) == NULL 
      || biop->genome == GENOME_unknown
      || biop->org == NULL) {
    return;
  }
  if (StringICmp (biop->org->taxname, "Human immunodeficiency virus") != 0 
      && StringICmp (biop->org->taxname, "Human immunodeficiency virus 1") != 0
      && StringICmp (biop->org->taxname, "Human immunodeficiency virus 2") != 0) {
    return;
  }

  if (biop->genome == GENOME_genomic) {
    msdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
    if (msdp != NULL && (mip = (MolInfoPtr) msdp->data.ptrvalue) != NULL && mip->biomol == MOLECULE_TYPE_GENOMIC) {
      return;
    }
  }
  ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQ, bsp);
}


static void FindInconsistentHIVRNA (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp;
  ValNodePtr item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, FindInconsistentHIVRNACallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, 
                       NewClickableItem (ONCALLER_HIV_RNA_INCONSISTENT, "%d HIV RNA bioseqs have inconsistent location/moltype", item_list));
  }
}


typedef struct bspprojectid {
  Int4 projectID;
  BioseqPtr bsp;
} BspProjectIdData, PNTR BspProjectIdPtr;


static BspProjectIdPtr BspProjectIdNew (BioseqPtr bsp, Int4 projectID)
{
  BspProjectIdPtr b;

  b = (BspProjectIdPtr) MemNew (sizeof (BspProjectIdData));
  b->projectID = projectID;
  b->bsp = bsp;
  return b;
}


static BspProjectIdPtr BspProjectIdFree (BspProjectIdPtr b)
{
  if (b != NULL) {
    b = MemFree (b);
  }
  return b;
}


static ValNodePtr BspProjectIdListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = BspProjectIdFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static int LIBCALLBACK SortVnpByBspProjectId (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  BspProjectIdPtr b1, b2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);

    if (vnp1 != NULL &&& vnp2 != NULL) {
      b1 = (BspProjectIdPtr) vnp1->data.ptrvalue;
      b2 = (BspProjectIdPtr) vnp2->data.ptrvalue;
      if (b1 != NULL && b2 != NULL) {
        if (b1->projectID < b2->projectID) {
          rval = -1;
        } else if (b1->projectID > b2->projectID) {
          rval = 1;
        }
      }
    }
  }

  return rval;
}


static void FindProjectIdSequenceCallback (BioseqPtr bsp, Pointer data)
{
  Int4 projectID;

  if (bsp == NULL || data == NULL) {
    return;
  }

  projectID = GetGenomeProjectID (bsp);
  if (projectID > 0) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, bsp->mol, BspProjectIdNew (bsp, projectID));
  }
}


static Boolean AllProjectIdsInListSame (ValNodePtr list)
{
  BspProjectIdPtr bid;
  Boolean         rval = TRUE;
  Int4            first_id;

  if (list == NULL || list->next == NULL) {
    return TRUE;
  }

  bid = (BspProjectIdPtr) list->data.ptrvalue;
  first_id = bid->projectID;
  list = list->next;
  while (list != NULL && rval) {
    bid = (BspProjectIdPtr) list->data.ptrvalue;
    if (first_id != bid->projectID) {
      rval = FALSE;
    }
    list = list->next;
  }
  return rval;
}


static void AddProjectIdSequencesFromList (ValNodePtr PNTR discrepancy_list, ValNodePtr list)
{
  ValNodePtr       vnp, subcat_items = NULL, item_list = NULL, subcat = NULL;
  BspProjectIdPtr  b;
  CharPtr          fmt = "%%d %s sequences have project ID %d";
  CharPtr          all_fmt = "%%d %s sequences have project IDs (%s)";
  Char             format[150];
  Int4             last_project_id = 0;
  ClickableItemPtr cip;
  Boolean          all_same;

  if (list == NULL) {
    return;
  }
  all_same = AllProjectIdsInListSame (list);
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    b = (BspProjectIdPtr) vnp->data.ptrvalue;
    if (b->projectID != last_project_id && last_project_id > 0) {
        sprintf (format, fmt, ISA_aa (b->bsp->mol) ? "protein" : "nucleotide", last_project_id); 
      ValNodeAddPointer (&subcat, 0, 
                         NewClickableItem (TEST_HAS_PROJECT_ID, format, subcat_items));
      subcat_items = NULL;
    }
    ValNodeAddPointer (&subcat_items, OBJ_BIOSEQ, b->bsp);
    ValNodeAddPointer (&item_list, OBJ_BIOSEQ, b->bsp);
    last_project_id = b->projectID;
  }
  if (last_project_id > 0) {
    sprintf (format, fmt, ISA_aa (b->bsp->mol) ? "protein" : "nucleotide", last_project_id); 
    ValNodeAddPointer (&subcat, 0, 
                       NewClickableItem (TEST_HAS_PROJECT_ID, format, subcat_items));
    subcat_items = NULL;
  }

  if (item_list != NULL) {
      sprintf (format, all_fmt, ISA_aa (b->bsp->mol) ? "protein" : "nucleotide", all_same ? "all same" : "some different");
    cip = NewClickableItem (TEST_HAS_PROJECT_ID, format, item_list);
    cip->subcategories = subcat;
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


static void FindProjectIdSequences (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp;
  ValNodePtr id_list = NULL, subcat = NULL;
  ValNodePtr prot_list = NULL;
  Int4       num_seq = 0;
  ClickableItemPtr cip;
  CharPtr          all_fmt_same = "%d sequences have project IDs (all same)";
  CharPtr          all_fmt_diff = "%d sequences have project IDs (some different)";
  CharPtr          all_fmt;
  Boolean          all_same = TRUE;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &id_list, FindProjectIdSequenceCallback);
  }

  num_seq = ValNodeLen (id_list);
  all_same = AllProjectIdsInListSame (id_list);
  id_list = ValNodeSort (id_list, SortVnpByBspProjectId);
  prot_list = ValNodeExtractList (&id_list, 3);

  AddProjectIdSequencesFromList (&subcat, id_list);
  AddProjectIdSequencesFromList (&subcat, prot_list);

  if (subcat != NULL) {
    if (id_list == NULL) {
      ValNodeLink (discrepancy_list, subcat);
    } else if (prot_list == NULL) {
      ValNodeLink (discrepancy_list, subcat);
    } else {
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      cip->clickable_item_type = TEST_HAS_PROJECT_ID;
      cip->subcategories = subcat;
      cip->expanded = 1;
      if (all_same) {
        all_fmt = all_fmt_same;
      } else {
        all_fmt = all_fmt_diff;
      }
      cip->description = (CharPtr) MemNew (sizeof (CharPtr) * (StringLen (all_fmt) + 15));
      sprintf (cip->description, all_fmt, num_seq);
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
  }
  id_list = BspProjectIdListFree(id_list);
  prot_list = BspProjectIdListFree(id_list);
}


static void FindSeqWithStructuredComments (BioseqPtr bsp, Pointer data)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  Uint1             num_present = 0;
  UserObjectPtr     uop;

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL) {
    return;
  }

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &context)) {
    if ((uop = (UserObjectPtr) sdp->data.ptrvalue) != NULL
        && uop->type != NULL
        && StringICmp (uop->type->str, "StructuredComment") == 0) {
      num_present++;
    }
  }
  ValNodeAddPointer ((ValNodePtr PNTR) data, num_present, bsp);
}


static void FindMissingStructuredComments (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr count_list = NULL;
  ValNodePtr tmp_list = NULL;
  ValNodePtr vnp;
  CharPtr    fmt;
  CharPtr    num_fmt = "%%d sequences have %d structured comments";
  ClickableItemPtr cip;
  ValNodePtr subcat = NULL;
  Uint1      orig_choice;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &count_list, FindSeqWithStructuredComments);
  }

  if (count_list == NULL) {
    return;
  }

  tmp_list = ValNodeExtractList (&count_list, 0);
  if (tmp_list == NULL) {
    /* no sequences have 0 */
    tmp_list = ValNodeExtractList (&count_list, count_list->choice);
  }
  if (count_list == NULL) {
    /* all sequences have same number of structured comments, no report */
    tmp_list = ValNodeFree (tmp_list);
  } else {
    while (tmp_list != NULL) {
      orig_choice = tmp_list->choice;
      for (vnp = tmp_list; vnp != NULL; vnp = vnp->next) {
        vnp->choice = OBJ_BIOSEQ;
      }
      fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (num_fmt) + 15));
      sprintf (fmt, num_fmt, orig_choice);
      cip = NewClickableItem (ONCALLER_MISSING_STRUCTURED_COMMENTS, fmt, tmp_list);
      fmt = MemFree (fmt);
      ValNodeAddPointer (&subcat, 0, cip);
      if (count_list == NULL) {
        tmp_list = NULL;
      } else {
        tmp_list = ValNodeExtractList (&count_list, count_list->choice);
      }
    }
    if (subcat->next == NULL) {
      subcat = FreeClickableList (subcat);
    } else {
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      MemSet (cip, 0, sizeof (ClickableItemData));
      cip->clickable_item_type = ONCALLER_MISSING_STRUCTURED_COMMENTS;
      cip->subcategories = subcat;
      cip->description = StringSave ("Sequences have different numbers of structured comments");
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
  }
}


static void MissingGenomeAssemblyStructuredCommentCallback (BioseqPtr bsp, Pointer data)
{
  SeqDescrPtr sdp;
  SeqMgrDescContext dcontext;
  Boolean found = FALSE;
  UserObjectPtr uop;
  UserFieldPtr ufp;

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL) {
    return;
  }
  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
       sdp != NULL && !found;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext)) {
    if ((uop = (UserObjectPtr) sdp->data.ptrvalue) != NULL
        && uop->type != NULL
        && StringICmp (uop->type->str, "StructuredComment") == 0) {
      for (ufp = uop->data; ufp != NULL && !found; ufp = ufp->next) {
        if (StringICmp (ufp->label->str, "StructuredCommentPrefix") == 0) {
          if (ufp->choice == 1 && StringICmp (ufp->data.ptrvalue, "##Genome-Assembly-Data-START##") == 0) {
            found = TRUE;
          }
          break;
        }
      }
    }
  }
  if (!found) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQ, bsp);
  }
}


static void FindMissingGenomeAssemblyStructuredComments (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, MissingGenomeAssemblyStructuredCommentCallback);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (MISSING_GENOMEASSEMBLY_COMMENTS, "%d bioseqs are missing GenomeAssembly structured comments", item_list));
  }
}


static void FindCDSWithCDDXrefCallback (SeqFeatPtr sfp, Pointer data)
{
  ValNodePtr vnp;
  DbtagPtr dbtag;
  Boolean  has_cdd_xref = FALSE;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || data == NULL) {
    return;
  }

  for (vnp = sfp->dbxref; vnp != NULL && !has_cdd_xref; vnp = vnp->next) {
    if ((dbtag = (DbtagPtr) vnp->data.ptrvalue) != NULL && StringICmp (dbtag->db, "CDD") == 0) {
      has_cdd_xref = TRUE;
    }
  }

  if (has_cdd_xref) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
}


static void FindCDSWithCDDXref (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &item_list, FindCDSWithCDDXrefCallback);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_CDS_HAS_CDD_XREF, "%d features have CDD Xrefs", item_list));
  }
}


static void LIBCALLBACK CountUnusualNTProc (CharPtr sequence, Pointer userdata)
{
  Int4Ptr p_i;
  CharPtr cp;

  if (sequence == NULL || userdata == NULL) return;
  p_i = (Int4Ptr) userdata;

  for (cp = sequence; *cp != 0; cp++)
  {
    if (*cp != 'N' && *cp != 'A' && *cp != 'T' && *cp != 'G' && *cp != 'C')
    {
      (*p_i) ++;
    }
  }
}


static void FindUnusualNTCallback (BioseqPtr bsp, Pointer data)
{
  Int4 num_bad = 0;
  Int4 flags = 0;

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL) {
    return;
  }

  SeqPortStream (bsp, flags, (Pointer) &num_bad, CountUnusualNTProc);
  if (num_bad > 0) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQ, bsp);
  }

}


static void FindUnusualNT (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, FindUnusualNTCallback);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_UNUSUAL_NT, "%d sequences contain nucleotides that are not ATCG or N", item_list));
  }
}


typedef struct qualityinterval {
  Int4 start;
  Int4 pos;
  Int4 num_ns;
  FloatLo min_pct;
  Int4 min_length;
  Boolean found_interval;
} QualityIntervalData, PNTR QualityIntervalPtr;


static void LIBCALLBACK FindLowQualityIntervalProc (CharPtr sequence, Pointer userdata)
{
  QualityIntervalPtr p_i;
  CharPtr cp;
  Int4    len;

  if (sequence == NULL || userdata == NULL) return;
  p_i = (QualityIntervalPtr) userdata;

  for (cp = sequence; *cp != 0; cp++)
  {
    if (*cp != 'A' && *cp != 'T' && *cp != 'G' && *cp != 'C') {
      if (p_i->start == -1) {
        /* start new interval if we aren't already in one */
        p_i->start = p_i->pos;
        p_i->num_ns = 1;
      } else {
        /* add to number of ns in this interval */
        p_i->num_ns++;
      }
    } else {
      if (p_i->start > -1) {
        /* if we are already in an interval, see if we should continue to be */
        len = p_i->pos - p_i->start;
        if ((FloatLo) p_i->num_ns / (FloatLo) len >= p_i->min_pct) {
          /* yes */
        } else {
          /* no */
          /* is the interval long enough to qualify? */
          if (len >= p_i->min_length) {
            p_i->found_interval = TRUE;
          }
          /* reset for next interval */
          p_i->start = -1;
          p_i->num_ns = 0;
        }
      }
    }
    p_i->pos ++;
  }
}


static void FindLowQualityRegionsCallback (BioseqPtr bsp, Pointer data)
{
  QualityIntervalData q;

  Int4 flags = 0;

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL) {
    return;
  }
  MemSet (&q, 0, sizeof (QualityIntervalData));
  q.start = -1;
  q.min_pct = 0.25;
  q.min_length = 30;

  SeqPortStream (bsp, flags, (Pointer) &q, FindLowQualityIntervalProc);
  /* check final interval, in case the end of the sequence is low quality */
  if (q.start > -1 && q.pos - q.start >= q.min_length) {
    q.found_interval = TRUE;
  }

  if (q.found_interval) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQ, bsp);
  }

}


static void FindLowQualityRegions (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, FindLowQualityRegionsCallback);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_LOW_QUALITY_REGION, "%d sequences contains low quality region", item_list));
  }
}


NLM_EXTERN Boolean IsLocationOrganelle (Uint1 genome)
{
  if (genome == GENOME_chloroplast
      || genome == GENOME_chromoplast
      || genome == GENOME_kinetoplast
      || genome == GENOME_mitochondrion
      || genome == GENOME_cyanelle
      || genome == GENOME_nucleomorph
      || genome == GENOME_apicoplast
      || genome == GENOME_leucoplast
      || genome == GENOME_proplastid
      || genome == GENOME_hydrogenosome
      || genome == GENOME_plastid
      || genome == GENOME_chromatophore) {
    return TRUE;
  } else {
    return FALSE;
  }
}

static void FindOrganelleNotGenomicCallback(BioseqPtr bsp, Pointer data)
{
  SeqDescPtr sdp;
  SeqMgrDescContext context;
  MolInfoPtr    mip;
  BioSourcePtr  biop;

  if (bsp == NULL || ISA_aa(bsp->mol) || data == NULL) {
    return;
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
  if (sdp == NULL || (mip = (MolInfoPtr) sdp->data.ptrvalue) == NULL) {
    return;
  } else if ((mip->biomol == MOLECULE_TYPE_GENOMIC || mip->biomol == 0) && bsp->mol == Seq_mol_dna) {
    return;
  }
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp != NULL && (biop = (BioSourcePtr) sdp->data.ptrvalue) != NULL
      && IsLocationOrganelle(biop->genome)) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}


static void FindOrganelleNotGenomic (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, FindOrganelleNotGenomicCallback);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_ORGANELLE_NOT_GENOMIC, "%d non-genomic sequences are organelles", item_list));
  }
}


static Boolean HasUnculturedNonOrganelleName (CharPtr taxname)
{
  if (StringCmp (taxname, "uncultured organism") == 0
      || StringCmp (taxname, "uncultured microorganism") == 0
      || StringCmp (taxname, "uncultured bacterium") == 0
      || StringCmp (taxname, "uncultured archaeon") == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static CharPtr kIntergenicSpacerNames[] = {
  "trnL-trnF intergenic spacer",
  "trnH-psbA intergenic spacer",
  "trnS-trnG intergenic spacer",
  "trnF-trnL intergenic spacer",
  "psbA-trnH intergenic spacer",
  "trnG-trnS intergenic spacer",
  NULL};

static Boolean HasIntergenicSpacerName(CharPtr str)
{
  Int4 i;
  Boolean rval = FALSE;

  for (i = 0; kIntergenicSpacerNames[i] != NULL && !rval; i++) {
    if (StringISearch (str, kIntergenicSpacerNames[i]) != NULL) {
      rval = TRUE;
    }
  }
  return rval;
}


static void FindUnwantedSpacersCallback(BioseqPtr bsp, Pointer data)
{
  SeqDescPtr sdp;
  SeqMgrDescContext context;
  BioSourcePtr  biop;
  SeqMgrFeatContext fcontext;
  SeqFeatPtr sfp;

  if (bsp == NULL || data == NULL) {
    return;
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp == NULL || (biop = (BioSourcePtr) sdp->data.ptrvalue) == NULL
      || biop->genome == GENOME_chloroplast || biop->genome == GENOME_plastid) {
    return;
  }
  /* shouldn't be uncultured non-organelle */
  if (biop != NULL && biop->org != NULL && HasUnculturedNonOrganelleName(biop->org->taxname)) {
    return;
  }

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_misc_feature, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_misc_feature, &fcontext)) {
    if (HasIntergenicSpacerName(sfp->comment)) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
    }
  }
}


static void FindUnwantedSpacers (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, FindUnwantedSpacersCallback);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_UNWANTED_SPACER, "%d suspect intergenic spacer notes not organelle", item_list));
  }
}


static SuspectRuleSetPtr OrganelleRules = NULL;
static Boolean OrganelleRuleReadAttempted = FALSE;

static SuspectRuleSetPtr ReadOrganelleRules(void)
{
  AsnIoPtr     aip;
  Char         buf [PATH_MAX];
  SuspectRuleSetPtr   rule_list;

  if (! FindPath("ncbi", "ncbi", "data", buf, sizeof (buf)))
  {
    Message (MSG_POSTERR, "Failed to find organelle product rules");
    return NULL;
  }

  StringCat(buf, "organelle_products.prt");

  aip = AsnIoOpen (buf, "r");
  if (aip == NULL) {
    Message (MSG_POSTERR, "Unable to open %s", buf);
    return NULL;
  }

  rule_list = SuspectRuleSetAsnRead (aip, NULL);
  if (rule_list == NULL) {
    Message (MSG_POSTERR, "Unable to read organelle product rule list from %s.", buf);
  }

  AsnIoClose (aip);
  return rule_list;
}


typedef struct findorganelleproducts {
  SuspectRuleSetPtr rule_list;
  ValNodePtr item_list;
} FindOrganelleProductsData, PNTR FindOrganelleProductsPtr;

static void FindOrganelleProductsCallback(BioseqPtr bsp, Pointer data)
{
  SeqDescPtr sdp;
  SeqMgrDescContext context;
  BioSourcePtr  biop;
  SeqMgrFeatContext fcontext, pcontext;
  SeqFeatPtr sfp, protsfp;
  ProtRefPtr prp;
  SuspectRulePtr rule;
  FindOrganelleProductsPtr fop;
  Boolean match;
  BioseqPtr protbsp;

  if (bsp == NULL || (fop = (FindOrganelleProductsPtr)data) == NULL) {
    return;
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp == NULL || (biop = (BioSourcePtr) sdp->data.ptrvalue) == NULL
      || biop->genome == GENOME_mitochondrion
      || biop->genome == GENOME_chloroplast 
      || biop->genome == GENOME_plastid) {
    return;
  }

  /* source should not be bacterial or viral */
  if (biop != NULL && biop->org != NULL && biop->org->orgname != NULL) {
    if (IsBacterialBioSource (biop) || IsViralBioSource(biop)) {
      return;
    }
  }

  /* shouldn't be uncultured non-organelle */
  if (biop != NULL && biop->org != NULL && HasUnculturedNonOrganelleName(biop->org->taxname)) {
    return;
  }

  /* look for misc_features */
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_misc_feature, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_misc_feature, &fcontext)) {
    if (StringNICmp (sfp->comment, "contains ", 9) == 0) {
      match = FALSE;
      for (rule = fop->rule_list; rule != NULL && !match; rule = rule->next) {
        match = DoesStringMatchSuspectRule (sfp->comment, sfp, rule);
      }
      if (match) {
        ValNodeAddPointer (&(fop->item_list), OBJ_SEQFEAT, sfp);
      }
    }
  }

  /* also look for coding regions */
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_CDS, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_CDS, &fcontext)) {
    protbsp = BioseqFindFromSeqLoc (sfp->product);
    protsfp = SeqMgrGetNextFeature (protbsp, NULL, 0, FEATDEF_PROT, &pcontext);
    if (protsfp != NULL && (prp = (ProtRefPtr) protsfp->data.value.ptrvalue) != NULL
      && prp->name != NULL) {
      match = FALSE;
      for (rule = fop->rule_list; rule != NULL && !match; rule = rule->next) {
        match = DoesStringMatchSuspectRule (prp->name->data.ptrvalue, sfp, rule);
      }
      if (match) {
        ValNodeAddPointer (&(fop->item_list), OBJ_SEQFEAT, sfp);
      }
    }
  }
}


static void FindOrganelleProducts(ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp;
  FindOrganelleProductsData fd;

  if (!OrganelleRuleReadAttempted) {
    OrganelleRules = ReadOrganelleRules();
    OrganelleRuleReadAttempted = TRUE;
  }
  if (OrganelleRules == NULL) {
    return;
  }

  MemSet (&fd, 0, sizeof (FindOrganelleProductsData));
  fd.rule_list = OrganelleRules;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &fd, FindOrganelleProductsCallback);
  }
  if (fd.item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_ORGANELLE_PRODUCTS, "%d suspect products not organelle", fd.item_list));
  }
}


static void FindBadMrnaQualCallback (BioseqPtr bsp, Pointer data)
{
  SeqDescPtr sdp;
  SeqMgrDescContext context;
  BioSourcePtr biop;
  SubSourcePtr ssp;
  Boolean found = FALSE;

  if (!IsMrnaSequence(bsp) || data == NULL) {
    return;
  }
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp == NULL || (biop = (BioSourcePtr) sdp->data.ptrvalue) == NULL) {
    return;
  }

  for (ssp = biop->subtype; ssp != NULL && !found; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_germline || ssp->subtype == SUBSRC_rearranged) {
      found = TRUE;
    }
  }
  if (found) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}


static void FindBadMrnaQual (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, FindBadMrnaQualCallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_BAD_MRNA_QUAL, "%d mRNA sequences have germline or rearranged qualifier", item_list));
  }
}


/* A warning when environmental sample qualifier is present and the organism name 
 * does not contain 'uncultured' or 'enrichment culture' or 'metagenome' or 'unidentified'
 * and the source does not have note (orgmod or subsrc) 
 * 'amplified with species-specific primers' 
 *  and the /metagenomic-source qualifier is not used
 */
static Boolean HasUnnecessaryEnvironmental(BioSourcePtr biop)
{
  SubSourcePtr ssp;
  OrgModPtr    mod;
  Boolean found = FALSE;
  Boolean has_note = FALSE;
  Boolean has_metagenomic = FALSE;

  if (biop == NULL) {
    return FALSE;
  }

  for (ssp = biop->subtype; ssp != NULL && !has_note && !has_metagenomic; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_environmental_sample) {
      found = TRUE;
    } else if (ssp->subtype == SUBSRC_other && StringISearch (ssp->name, "amplified with species-specific primers") != NULL) {
      has_note = TRUE;
    } else if (ssp->subtype == SUBSRC_metagenomic) {
      has_metagenomic = TRUE;
    }
  }

  if (!found || has_note || has_metagenomic) {
    return FALSE;
  }
  if (biop->org != NULL) {
    if (StringISearch (biop->org->taxname, "uncultured") != NULL
        || StringISearch (biop->org->taxname, "enrichment culture") != NULL
        || StringISearch (biop->org->taxname, "metagenome") != NULL
        || StringISearch (biop->org->taxname, "environmental sample") != NULL
        || StringISearch (biop->org->taxname, "unidentified") != NULL) {
      return FALSE;
    }
    if (biop->org->orgname != NULL) {
      for (mod = biop->org->orgname->mod; mod != NULL && !has_note; mod = mod->next) {
        if (mod->subtype == ORGMOD_other && StringISearch (mod->subname, "amplified with species-specific primers") != NULL) {
          has_note = TRUE;
        }
      }
      if (has_note) {
        return FALSE;
      }
    }
  }
  return TRUE;
}


static void FindUnnecessaryEnvironmental (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    ValNodeLink (&item_list, RunBioSourceTest (vnp->data.ptrvalue, HasUnnecessaryEnvironmental));
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_UNNECESSARY_ENVIRONMENTAL, "%d biosources have unnecessary environmental qualifier", item_list));
  }
}


static void FindUnnecessaryVirusGeneCallback(BioseqPtr bsp, Pointer data)
{
  BioSourcePtr biop;
  SeqMgrFeatContext context;
  SeqFeatPtr sfp;

  if (bsp == NULL || data == NULL || ISA_aa(bsp->mol)) {
    return;
  }

  biop = GetBiopForBsp(bsp);
  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) {
    return;
  }
  if (HasLineage (biop, "Picornaviridae")
      || HasLineage (biop, "Potyviridae")
      || HasLineage (biop, "Flaviviridae")
      || HasLineage (biop, "Togaviridae")) {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &context);
         sfp != NULL;
         sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_GENE, 0, &context)) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
    }
  }
}


static void FindUnnecessaryVirusGene (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, FindUnnecessaryVirusGeneCallback);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_UNNECESSARY_VIRUS_GENE, "%d virus genes need to be removed", item_list));
  }
}


typedef struct isunwanted {
  Boolean has_sat_feat;
  Boolean has_non_sat_feat;
  Boolean has_rearranged;
} IsUnwantedData, PNTR IsUnwantedPtr;


static Boolean IsMicrosatelliteRepeatRegion (SeqFeatPtr sfp)
{
  GBQualPtr qual;
  Boolean   rval = FALSE;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_repeat_region) {
    return FALSE;
  }
  for (qual = sfp->qual; qual != NULL && !rval; qual = qual->next) {
    if (StringICmp (qual->qual, "satellite") == 0 && StringNICmp (qual->val, "microsatellite", 14) == 0) {
      rval = TRUE;
    }
  }
  return rval;
}


static void FindUnwantedSetWrappersCallback(BioseqPtr bsp, Pointer data)
{
  IsUnwantedPtr up;
  SeqFeatPtr sfp;
  SeqMgrFeatContext context;
  BioSourcePtr biop;
  SubSourcePtr ssp;

  if (bsp == NULL || ISA_aa(bsp->mol) || (up = (IsUnwantedPtr) data) == NULL) {
    return;
  }

  biop = GetBiopForBsp(bsp);
  if (biop != NULL) {
    for (ssp = biop->subtype; ssp != NULL && !up->has_rearranged; ssp = ssp->next) {
      if (ssp->subtype == SUBSRC_rearranged) {
        up->has_rearranged = TRUE;
      }
    }
  }

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
       sfp != NULL && (!up->has_sat_feat || !up->has_non_sat_feat);
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context)) {
    if (IsMicrosatelliteRepeatRegion(sfp)) {
      up->has_sat_feat = TRUE;
    } else {
      up->has_non_sat_feat = TRUE;
    }
  }
}


static void FindUnwantedSetWrappersInSep(SeqEntryPtr sep, ValNodePtr PNTR pList)
{
  BioseqSetPtr bssp;
  IsUnwantedData ud;

  if (sep == NULL || !IS_Bioseq_set(sep) || (bssp = (BioseqSetPtr) sep->data.ptrvalue) == NULL || pList == NULL) {
    return;
  }
 
  if (bssp->_class == BioseqseqSet_class_eco_set
      || bssp->_class == BioseqseqSet_class_mut_set
      || bssp->_class == BioseqseqSet_class_phy_set
      || bssp->_class == BioseqseqSet_class_pop_set) {
    MemSet (&ud, 0, sizeof (IsUnwantedData));
    VisitBioseqsInSep (sep, &ud, FindUnwantedSetWrappersCallback);

    if (ud.has_rearranged || (ud.has_sat_feat && !ud.has_non_sat_feat)) {
      ValNodeAddPointer (pList, OBJ_BIOSEQSET, bssp);
    }
  } else {
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
      FindUnwantedSetWrappersInSep (sep, pList);
    }
  }
}


static void FindUnwantedSetWrappers (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    FindUnwantedSetWrappersInSep (vnp->data.ptrvalue, &item_list);
  }

  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_UNWANTED_SET_WRAPPER, "%d unwanted set wrappers", item_list));
  }
}


static Boolean IsMissingPrimerValue (BioSourcePtr biop)
{
  PCRReactionSetPtr set;
  PCRPrimerPtr fwd, rev;
  Boolean rval = FALSE;

  if (biop == NULL) {
    return FALSE;
  }
  for (set = biop->pcr_primers; set != NULL && !rval; set = set->next) {
    for (fwd = set->forward, rev = set->reverse;
         fwd != NULL && rev != NULL && !rval;
         fwd = fwd->next, rev = rev->next) {
      if ((StringHasNoText(fwd->name) && !StringHasNoText(rev->name))
          || (!StringHasNoText (fwd->name) && StringHasNoText (rev->name))
          || (StringHasNoText(fwd->seq) && !StringHasNoText(rev->seq))
          || (!StringHasNoText (fwd->seq) && StringHasNoText (rev->seq))) {
        rval = TRUE;
      }
    }
    if (fwd != NULL || rev != NULL) {
      rval = TRUE;
    }
  }
  return rval;
}


static void FindMissingPrimerValues (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    ValNodeLink (&item_list, RunBioSourceTest (vnp->data.ptrvalue, IsMissingPrimerValue));
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_MISSING_PRIMER, "%d biosources have primer sets with missing values", item_list));
  }
}


static void FindUnexpectedMiscRNABioseq (BioseqPtr bsp, Pointer data)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext context;
  CharPtr   product;

  if (bsp == NULL || data == NULL) {
    return;
  }

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_otherRNA, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_otherRNA, &context)) {
     product = GetRNARefProductString(sfp->data.value.ptrvalue, NULL);
     if (StringSearch (product, "ITS") == NULL && StringSearch (product, "internal transcribed spacer") == NULL) {
       ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
     }
     product = MemFree (product);
  }
}


static void FindUnexpectedMiscRNA (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, FindUnexpectedMiscRNABioseq);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_UNUSUAL_MISC_RNA, "%d unexpected misc_RNA features found.  misc_RNAs are unusual in a genome, consider using ncRNA, misc_binding, or misc_feature as appropriate.", item_list));
  }
}


static Boolean AmpPrimersNoEnvSample (BioSourcePtr biop)
{
  OrgModPtr mod;
  SubSourcePtr ssp;
  Boolean has_note = FALSE;

  if (biop == NULL) {
    return FALSE;
  }

  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_environmental_sample) {
      return FALSE;
    } else if (ssp->subtype == SUBSRC_other 
               && StringISearch (ssp->name, "amplified with species-specific primers") != NULL) {
      has_note = TRUE;
    }
  }

  if (!has_note && biop->org != NULL && biop->org->orgname != NULL) {
    for (mod = biop->org->orgname->mod; mod != NULL && !has_note; mod = mod->next) {
      if (mod->subtype == SUBSRC_other
          && StringISearch (mod->subname, "amplified with species-specific primers") != NULL) {
        has_note = TRUE;
      }
    }
  }

  return has_note;
}


static void FindAmpPrimersNoEnvSample (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    ValNodeLink (&item_list, RunBioSourceTest (vnp->data.ptrvalue, AmpPrimersNoEnvSample));
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_AMPLIFIED_PRIMERS_NO_ENVIRONMENTAL_SAMPLE, "%d biosources have 'amplified with species-specific primers' note but no environmental-sample qualifier.", item_list));
  }
}


static void FindDuplicateGenesOnOppositeStrandsCallback (BioseqPtr bsp, Pointer data)
{
  SeqFeatPtr sfp, sfp_prev = NULL;
  SeqMgrFeatContext context;
  Boolean sfp_prev_listed = FALSE;

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL) {
    return;
  }
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_GENE, 0, &context)) {
    if (sfp_prev != NULL) {
      if (SeqLocCompare (sfp_prev->location, sfp->location) == SLC_A_EQ_B
          && SeqLocStrand (sfp_prev->location) != SeqLocStrand (sfp->location)) {
        if (!sfp_prev_listed) {
          ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp_prev);
        }
        ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
        sfp_prev_listed = TRUE;
      } else {
        sfp_prev_listed = FALSE;
      }
    }
    sfp_prev = sfp;
  }
}


static void FindDuplicateGenesOnOppositeStrands (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, item_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &item_list, FindDuplicateGenesOnOppositeStrandsCallback);
  }
  if (item_list != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_DUP_GENES_OPPOSITE_STRANDS, "%d genes match other genes in the same location, but on the opposite strand", item_list));
  }
}


static void FindSmallGenomeSetCallback (BioseqSetPtr bssp, Pointer data)
{
  if (bssp != NULL && bssp->_class == BioseqseqSet_class_small_genome_set && data != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQSET, bssp);
  }
}


static void ListBioSources(SeqDescrPtr sdp, Pointer data)
{
  if (sdp != NULL && sdp->choice == Seq_descr_source) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}


static void FindSmallGenomeSetProblems (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr vnp, src_list = NULL, s;
  CharPtr    taxname = NULL, strain = NULL, isolate = NULL;
  CharPtr    tmp;
  BioSourcePtr biop;
  ValNodePtr   tax_qual, strain_qual, isolate_qual, segment_qual, div_qual;
  ValNodePtr   missing_segment = NULL;
  Boolean      all_taxnames_same = TRUE;
  Boolean      all_isolates_same = TRUE;
  Boolean      all_strains_same = TRUE;
  ValNodePtr   set_list = NULL, vnp_s;
  BioseqSetPtr bssp;
  SeqEntryPtr  sep;

  tax_qual = ValNodeNew (NULL);
  tax_qual->choice = SourceQualChoice_textqual;
  tax_qual->data.intvalue = Source_qual_taxname;
  strain_qual = ValNodeNew (NULL);
  strain_qual->choice = SourceQualChoice_textqual;
  strain_qual->data.intvalue = Source_qual_strain;
  isolate_qual = ValNodeNew (NULL);
  isolate_qual->choice = SourceQualChoice_textqual;
  isolate_qual->data.intvalue = Source_qual_isolate;
  segment_qual = ValNodeNew (NULL);
  segment_qual->choice = SourceQualChoice_textqual;
  segment_qual->data.intvalue = Source_qual_segment;
  div_qual = ValNodeNew (NULL);
  div_qual->choice = SourceQualChoice_textqual;
  div_qual->data.intvalue = Source_qual_division;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitSetsInSep (vnp->data.ptrvalue, &set_list, FindSmallGenomeSetCallback);
    for (vnp_s = set_list; vnp_s != NULL; vnp_s = vnp_s->next) {
      bssp = (BioseqSetPtr) vnp_s->data.ptrvalue;
      sep = SeqMgrGetSeqEntryForData (bssp);
      VisitDescriptorsInSep (sep, &src_list, ListBioSources);
      all_taxnames_same = TRUE;
      all_isolates_same = TRUE;
      all_strains_same = TRUE;
      for (s = src_list; s != NULL; s = s->next) {
        biop = GetBioSourceFromObject(s->choice, s->data.ptrvalue);
        if (biop != NULL) {
          /* look for segment when required */
          if (IsViralBioSource(biop)) {
            tmp = GetSourceQualFromBioSource(biop, segment_qual, NULL);
            if (tmp == NULL) {
              ValNodeAddPointer (&missing_segment, OBJ_SEQDESC, s->data.ptrvalue);
            }
            tmp = MemFree (tmp);
          }
          /* are taxnames all the same */
          if (all_taxnames_same) {
            tmp = GetSourceQualFromBioSource(biop, tax_qual, NULL);
            if (tmp != NULL) {
              if (s == src_list) {
                taxname = tmp;
                tmp = NULL;
              } else if (StringCmp (taxname, tmp) != 0) {
                all_taxnames_same = FALSE;
              }
              tmp = MemFree (tmp);
            }
          }
          /* are isolates all the same */
          if (all_isolates_same) {
            tmp = GetSourceQualFromBioSource(biop, isolate_qual, NULL);
            if (tmp != NULL) {
              if (s == src_list) {
                isolate = tmp;
                tmp = NULL;
              } else if (StringCmp (isolate, tmp) != 0) {
                all_isolates_same = FALSE;
              }
              tmp = MemFree (tmp);
            }
          }
          /* are strains all the same */
          if (all_strains_same) {
            tmp = GetSourceQualFromBioSource(biop, strain_qual, NULL);
            if (tmp != NULL) {
              if (s == src_list) {
                strain = tmp;
                tmp = NULL;
              } else if (StringCmp (strain, tmp) != 0) {
                all_strains_same = FALSE;
              }
              tmp = MemFree (tmp);
            }
          }
        }
      }
      src_list = FreeObjectList (src_list);
      taxname = MemFree (taxname);
      isolate = MemFree (isolate);
      strain = MemFree (strain);

      if (!all_taxnames_same) {
        ValNodeAddPointer (discrepancy_list, 0, NewClickableItemNoList (TEST_SMALL_GENOME_SET_PROBLEM, "Not all biosources have same taxname"));
      }
      if (!all_isolates_same) {
        ValNodeAddPointer (discrepancy_list, 0, NewClickableItemNoList (TEST_SMALL_GENOME_SET_PROBLEM, "Not all biosources have same isolate"));
      }
      if (!all_strains_same) {
        ValNodeAddPointer (discrepancy_list, 0, NewClickableItemNoList (TEST_SMALL_GENOME_SET_PROBLEM, "Not all biosources have same strain"));
      }
    }
    set_list = ValNodeFree (set_list);
  }
  if (missing_segment != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_SMALL_GENOME_SET_PROBLEM, "%d biosources should have segment qualifier but do not", missing_segment));
  }


}


static void FindOverlappingrRNAs (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr         sfp, sfp_compare;
  SeqMgrFeatContext  context;
  ValNodePtr PNTR    overlapping_rrnas = NULL, non_overlap;
  ValNodePtr         rrna_list = NULL, vnp, vnp_next;
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  overlapping_rrnas = (ValNodePtr PNTR) userdata;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_rRNA, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_rRNA, &context))
  {
    ValNodeAddPointer (&rrna_list, 0, sfp);
  }
  
  for (vnp = rrna_list; vnp != NULL && vnp->next != NULL; vnp = vnp->next)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    for (vnp_next = vnp->next; vnp_next != NULL; vnp_next = vnp_next->next)
    {
      sfp_compare = (SeqFeatPtr) vnp_next->data.ptrvalue;
            
      if (SeqLocCompare (sfp->location, sfp_compare->location) != SLC_NO_MATCH)
      {
        vnp->choice = OBJ_SEQFEAT;
        vnp_next->choice = OBJ_SEQFEAT;
      }
    }
  }
  
  non_overlap = ValNodeExtractList (&rrna_list, 0);
  non_overlap = ValNodeFree (non_overlap);
  ValNodeLink (overlapping_rrnas, rrna_list);
  
}


extern void AddOverlappingrRNADiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  CharPtr            bad_fmt = "%d rRNA features overlap another rRNA feature.";
  ValNodePtr         overlapping_rrnas = NULL, vnp;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &overlapping_rrnas, FindOverlappingrRNAs);
  }
  
  if (overlapping_rrnas != NULL)
  {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_OVERLAPPING_RRNAS, bad_fmt, overlapping_rrnas));
  }
}


static void FindMrnaSequencesWithMinusStrandFeaturesCallback (BioseqPtr bsp, Pointer data)
{
  SeqMgrFeatContext context;
  SeqFeatPtr sfp;
  Boolean found = FALSE;

  if (bsp == NULL || !IsMrnaSequence(bsp) || data == NULL) {
    return;
  }

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
       sfp != NULL && !found;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context)) {
    if (context.strand == Seq_strand_minus) {
      found = TRUE;
    }
  }
  if (found) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_BIOSEQ, bsp);
  }
}


static void FindMrnaSequencesWithMinusStrandFeatures (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  CharPtr            bad_fmt = "%d mRNA sequences have features on the complement strand.";
  ValNodePtr         seqs = NULL, vnp;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &seqs, FindMrnaSequencesWithMinusStrandFeaturesCallback);
  }
  
  if (seqs != NULL)
  {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_MRNA_SEQUENCE_MINUS_STRAND_FEATURES, bad_fmt, seqs));
  }
}


static void FindTaxnameMissingFromDeflineCallback (BioseqPtr bsp, Pointer data)
{
  SeqMgrDescContext context;
  SeqDescPtr sdp;
  BioSourcePtr biop;
  CharPtr cp;
  Int4    len;
  CharPtr lookfor;

  if (bsp == NULL || ISA_aa(bsp->mol) || data == NULL) {
    return;
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp == NULL || (biop = (BioSourcePtr) sdp->data.ptrvalue) == NULL
      || biop->org == NULL
      || StringHasNoText (biop->org->taxname)) {
    return;
  }

  lookfor = biop->org->taxname;
  if (StringICmp (lookfor, "Human immunodeficiency virus 1") == 0) {
    lookfor = "HIV-1";
  } else if (StringICmp (lookfor, "Human immunodeficiency virus 2") == 0) {
    lookfor = "HIV-2";
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &context);
  if (sdp != NULL) {
    cp = StringISearch (sdp->data.ptrvalue, lookfor);
    if (cp == NULL) {
      /* taxname not in defline at all */
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
    } else {
      /* capitalization must match for all but the first letter */
      len = StringLen (lookfor);
      if (StringNCmp (cp + 1, lookfor + 1, len - 1) != 0) {
        ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
      }
    }
  }
}


static void FindTaxnameMissingFromDefline (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  CharPtr            bad_fmt = "%d deflines do not contain the complete taxname.";
  ValNodePtr         seqs = NULL, vnp;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &seqs, FindTaxnameMissingFromDeflineCallback);
  }
  
  if (seqs != NULL)
  {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_TAXNAME_NOT_IN_DEFLINE, bad_fmt, seqs));
  }
}


static Boolean IsUnverified (BioseqPtr bsp)
{
  if (bsp != NULL && !ISA_aa (bsp->mol)) {
    return BioseqHasKeyword (bsp, "UNVERIFIED");
  } else {
    return FALSE;
  }
}


static void CountUnverifiedSequences  (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  CharPtr            bad_fmt = "%d sequences are unverified.";
  ValNodePtr         seqs = NULL, vnp;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    ValNodeLink (&seqs, RunBioseqTest (vnp->data.ptrvalue, IsUnverified));
  }
  
  if (seqs != NULL)
  {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (TEST_COUNT_UNVERIFIED, bad_fmt, seqs));
  }
}


static void 
RemoveUnwantedDiscrepancyItems 
(ValNodePtr PNTR      discrepancy_list,
 DiscrepancyConfigPtr dcp)
{
  ValNodePtr         vnp, prev = NULL, vnp_next;
  ClickableItemPtr dip;
  
  if (dcp == NULL || discrepancy_list == NULL || *discrepancy_list == NULL)
  {
    return;
  }
  
  for (vnp = *discrepancy_list; vnp != NULL; vnp = vnp_next)
  {
    vnp_next = vnp->next;
    dip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (dip == NULL || ! dcp->conf_list[dip->clickable_item_type])
    {
      if (prev == NULL)
      {
        *discrepancy_list = vnp_next;
      }
      else
      {
        prev->next = vnp_next;
      }
      vnp->next = NULL;
      vnp = FreeClickableList (vnp);
    }
    else
    {
      prev = vnp;
    }
  }
  
}


extern void SetDiscrepancyLevels (ValNodePtr discrepancy_list, Int4 level)
{
  ClickableItemPtr dip;
  
  while (discrepancy_list != NULL)
  {
    dip = (ClickableItemPtr) discrepancy_list->data.ptrvalue;
    if (dip != NULL)
    {
      dip->level = level;
      SetDiscrepancyLevels (dip->subcategories, level + 1);
    }
    discrepancy_list = discrepancy_list->next;
  }
}


typedef struct discrepancyinfo 
{
  CharPtr                conf_name;
  CharPtr                setting_name;
  PerformDiscrepancyTest test_func;
  AutofixCallback        autofix_func;
} DiscrepancyInfoData, PNTR DiscrepancyInfoPtr;


/*  "Runs of 20 or more Ns" -> "Runs of 10 or more Ns": per Larissa's request, by J. Chen */
/*  "Show translation exception": test if code-break exists, by J. Chen */
/*  "Show hypothetic protein having a gene name":  J. Chen */
/*  "Test defline existence": J. Chen */
/*  "Remove mRNA overlapping a pseudogene": J. Chen */
static DiscrepancyInfoData discrepancy_info_list[] = 
{
  { "Missing Genes", "MISSING_GENES", AddMissingAndSuperfluousGeneDiscrepancies, NULL },
  { "Extra Genes", "EXTRA_GENES", AddMissingAndSuperfluousGeneDiscrepancies, NULL },
  { "Missing Locus Tags", "MISSING_LOCUS_TAGS", AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags, NULL },
  { "Duplicate Locus Tags", "DUPLICATE_LOCUS_TAGS", AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags, NULL },
  { "Bad Locus Tag Format", "BAD_LOCUS_TAG_FORMAT", AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags, NULL },
  { "Inconsistent Locus Tag Prefix", "INCONSISTENT_LOCUS_TAG_PREFIX", AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags, NULL },
  { "Nongene Locus Tag", "NON_GENE_LOCUS_TAG", AddDiscrepanciesForNonGeneLocusTags, NULL },
  { "Count nucleotide sequences", "DISC_COUNT_NUCLEOTIDES", CountNucSeqs, NULL},
  { "Missing Protein ID", "MISSING_PROTEIN_ID", FindMissingProteinIDs, NULL },
  { "Inconsistent Protein ID", "INCONSISTENT_PROTEIN_ID", FindMissingProteinIDs, NULL },
  { "Feature Location Conflict", "FEATURE_LOCATION_CONFLICT", FindCDSmRNAGeneLocationDiscrepancies, NULL },
  { "Gene Product Conflict", "GENE_PRODUCT_CONFLICT", FindCDSGeneProductConflicts, NULL },
  { "Duplicate Gene Locus", "DUPLICATE_GENE_LOCUS", FindDuplicateGeneLocus, NULL },
  { "EC Number Note", "EC_NUMBER_NOTE", AddECNumberNoteDiscrepancies, NULL },
  { "Pseudo Mismatch", "PSEUDO_MISMATCH", FindPseudoDiscrepancies, OncallerToolPseudoDiscrepanciesFix },
  { "Joined Features", "JOINED_FEATURES", AddJoinedFeatureDiscrepancies, NULL },
  { "Overlapping Genes", "OVERLAPPING_GENES", AddOverlappingGeneDiscrepancies, NULL },
  { "Overlapping CDS", "OVERLAPPING_CDS", AddOverlappingCodingRegionDiscrepancies, MarkOverlappingCDSs },
  { "Contained CDS", "CONTAINED_CDS", AddContainedCodingRegionDiscrepancies, NULL },
  { "CDS RNA Overlap", "RNA_CDS_OVERLAP", AddRNACDSOverlapDiscrepancies, NULL },
  { "Short Contig", "SHORT_CONTIG", FindShortContigs, RemoveShortContigsWithoutAnnotation },
  { "Inconsistent BioSource", "INCONSISTENT_BIOSOURCE", FindNonmatchingContigSources, NULL },
  { "Suspect Product Name", "SUSPECT_PRODUCT_NAMES", FindSuspectProductNames, NULL },
  { "Suspect Product Name Typo", "DISC_PRODUCT_NAME_TYPO", FindSuspectProductNames, FixSuspectProductNameTypos },
  { "Suspect Product Name QuickFix", "DISC_PRODUCT_NAME_QUICKFIX", FindSuspectProductNames, FixSuspectProductNameQuickFixes },
  { "Inconsistent Source And Definition Line", "INCONSISTENT_SOURCE_DEFLINE", FindInconsistentSourceAndDefline, NULL },
  { "Partial CDSs in Complete Sequences", "PARTIAL_CDS_COMPLETE_SEQUENCE", FindParticalCDSsInCompleteSequences, NULL },
  { "Hypothetical or Unknown Protein with EC Number", "EC_NUMBER_ON_UNKNOWN_PROTEIN", FindUnknownProteinsWithECNumbers, NULL },
  { "Find Missing Tax Lookups", "TAX_LOOKUP_MISSING", NULL, NULL } ,
  { "Find Tax Lookup Mismatches", "TAX_LOOKUP_MISMATCH", NULL, NULL },
  { "Find Short Sequences", "SHORT_SEQUENCES", FindShortSequences, NULL },
  { "Suspect Phrases", "SUSPECT_PHRASES", FindSuspectPhrases, NULL },
  { "Find Suspicious Phrases in Note Text", "DISC_SUSPICIOUS_NOTE_TEXT", FindSuspiciousPhraseInNoteText, NULL},
  { "Count tRNAs", "COUNT_TRNAS", tRNACountFeaturesAndFindDups, NULL },
  { "Find Duplicate tRNAs", "FIND_DUP_TRNAS", tRNACountFeaturesAndFindDups, NULL },
  { "Find short and long tRNAs", "FIND_BADLEN_TRNAS", tRNAFindBadLength, NULL},
  { "Find tRNAs on the same strand", "FIND_STRAND_TRNAS", FindtRNAsOnSameStrand, NULL},
  { "Count rRNAs", "COUNT_RRNAS", rRNACountFeaturesAndFindDups, NULL },
  { "Find Duplicate rRNAs", "FIND_DUP_RRNAS", rRNACountFeaturesAndFindDups, NULL },
  { "Find RNAs without Products", "RNA_NO_PRODUCT", FindRNAsWithoutProducts, NULL },
  { "Transl_except without Note", "TRANSL_NO_NOTE", FindTranslExceptNotes, NULL },
  { "Note without Transl_except", "NOTE_NO_TRANSL", FindTranslExceptNotes, NULL },
  { "Transl_except longer than 3", "TRANSL_TOO_LONG", FindTranslExceptNotes, NULL },
  { "CDS tRNA overlaps", "CDS_TRNA_OVERLAP", FindCDSOverlappingtRNAs, NULL },
  { "Count Proteins", "COUNT_PROTEINS", CountProteins, NULL },
  { "Features Intersecting Source Features", "DISC_FEAT_OVERLAP_SRCFEAT", FindFeaturesOverlappingSrcFeatures, NULL },
  { "CDS on GenProdSet without protein", "MISSING_GENPRODSET_PROTEIN", CheckListForGenProdSets, NULL},
  { "Multiple CDS on GenProdSet, same protein", "DUP_GENPRODSET_PROTEIN", CheckListForGenProdSets, NULL},
  { "mRNA on GenProdSet without transcript ID", "MISSING_GENPRODSET_TRANSCRIPT_ID", CheckListForGenProdSets, NULL},
  { "mRNA on GenProdSet with duplicate ID", "DISC_DUP_GENPRODSET_TRANSCRIPT_ID", CheckListForGenProdSets, NULL},
  { "Greater than 5 percent Ns", "DISC_PERCENT_N", PercentNDiscrepanciesForSeqEntry, NULL},
  { "Runs of 10 or more Ns", "N_RUNS", BaseCountAndNRunDiscrepancies, NULL},
  { "Zero Base Counts", "ZERO_BASECOUNT", BaseCountAndNRunDiscrepancies, NULL},
  { "Adjacent PseudoGenes with Identical Text", "ADJACENT_PSEUDOGENES", FindAdjacentPseudoGenes, NULL},
  { "Bioseqs without Annotations", "NO_ANNOTATION", FindBioseqsWithoutAnnotation, NULL},
  { "Influenza Strain/Collection Date Mismatch", "DISC_INFLUENZA_DATE_MISMATCH", FindInfluenzaStrainCollectionDateMismatches, NULL},
  { "Introns shorter than 10 nt", "DISC_SHORT_INTRON", FindShortIntrons, AddExceptionsToShortIntrons},
  { "Viruses should specify collection-date, country, and specific-host", "DISC_MISSING_VIRAL_QUALS", FindMissingViralQuals, NULL},
  { "Source Qualifier Report", "DISC_SRC_QUAL_PROBLEM", CheckBioSourceQuals, NULL},
  { "All sources in a record should have the same qualifier set", "DISC_MISSING_SRC_QUAL", CheckBioSourceQuals, NULL},
  { "Each source in a record should have unique values for qualifiers", "DISC_DUP_SRC_QUAL", CheckBioSourceQuals, NULL},
  { "Each qualifier on a source should have different values", "DISC_DUP_SRC_QUAL_DATA", CheckBioSourceQuals, NULL},
  { "Sequences with the same haplotype should match", "DISC_HAPLOTYPE_MISMATCH", ReportHaplotypeSequenceMismatch, NULL},
  { "Sequences with rRNA or misc_RNA features should be genomic DNA", "DISC_FEATURE_MOLTYPE_MISMATCH", ReportFeatureMoltypeMismatch, ChangeMoltypeToGenomicDNA},
  { "Coding regions on eukaryotic genomic DNA should have mRNAs with matching products", "DISC_CDS_WITHOUT_MRNA", ReportCDSWithoutmRNA, AddMissingmRNA},
  { "Exon and intron locations should abut (unless gene is trans-spliced)", "DISC_EXON_INTRON_CONFLICT", CheckIntronAndExonLocations, NULL},
  { "Count features present or missing from sequences", "DISC_FEATURE_COUNT", CountFeaturesOnSequences, NULL},
  { "BioSources with the same specimen voucher should have the same taxname", "DISC_SPECVOUCHER_TAXNAME_MISMATCH", CollectSpecVoucherTaxnameDiscrepancies, NULL},
  { "Feature partialness should agree with gene partialness if endpoints match", "DISC_GENE_PARTIAL_CONFLICT", ReportPartialConflicts, NULL},
  { "Flatfile representation of object contains suspect text", "DISC_FLATFILE_FIND_ONCALLER", FindTextInFlatfileOncaller, OncallerToolSpellFix},
  { "Coding region product contains suspect text", "DISC_CDS_PRODUCT_FIND", FindTextInCDSProduct, NULL},
  { "Definition lines should be unique", "DISC_DUP_DEFLINE", FindDupDeflines, NULL},
  { "ATCC strain should also appear in culture collection", "DUP_DISC_ATCC_CULTURE_CONFLICT", CheckATCCStrainCultureCollConflict, AddATCCStrainToCultureColl},
  { "For country USA, state should be present and abbreviated", "DISC_USA_STATE", CheckUSAStates, FixUSAStates},
  { "All non-protein sequences in a set should have the same moltype", "DISC_INCONSISTENT_MOLTYPES", CheckMoltypes, NULL},
  { "Records should have identical submit-blocks", "DISC_SUBMITBLOCK_CONFLICT", CheckSubmitBlockConflicts, NULL},
  { "Possible linker sequence after poly-A tail", "DISC_POSSIBLE_LINKER", CheckForLinkerSequence, NULL},
  { "Publications with the same titles should have the same authors", "DISC_TITLE_AUTHOR_CONFLICT", CheckForTitleAuthorConflicts, NULL},
  { "Genes and features that share endpoints should be on the same strand", "DISC_BAD_GENE_STRAND", CheckGeneFeatureStrandConflicts, NULL},
  { "Eukaryotic sequences with a map source qualifier should also have a chromosome source qualifier", "DISC_MAP_CHROMOSOME_CONFLICT", CheckForMapChromosomeConflicts, NULL},
  { "RBS features should have an overlapping gene", "DISC_RBS_WITHOUT_GENE", CheckForRBSWithoutGene, NULL},
  { "All Cit-subs should have identical affiliations", "DISC_CITSUBAFFIL_CONFLICT", FindMismatchedCitSubAffiliations, NULL},
  { "Uncultured or environmental sources should have clone", "DISC_REQUIRED_CLONE", FindRequiredClones, NULL},
  { "Source Qualifier test for Asndisc", "DISC_SOURCE_QUALS_ASNDISC", CheckBioSourceQualsAsnDisc, NULL},
  { "Eukaryotic sequences that are not genomic or macronuclear should not have mRNA features", "DISC_mRNA_ON_WRONG_SEQUENCE_TYPE", ReportmRNAOnNonGenomicEukaryoticSequences, NULL},
  { "When the organism lineage contains 'Retroviridae' and the molecule type is 'DNA', the location should be set as 'proviral'", "DISC_RETROVIRIDAE_DNA", CheckRetroviridaeDNA, MakeLocationProviral},
  { "Check for correct capitalization in author names", "DISC_CHECK_AUTH_CAPS", CheckAuthCaps, FixAuthCaps},
  { "Check for gene or genes in rRNA and tRNA products and comments", "DISC_CHECK_RNA_PRODUCTS_AND_COMMENTS", CheckRNAProductsAndComments, NULL},
  { "Microsatellites must have repeat type of tandem", "DISC_MICROSATELLITE_REPEAT_TYPE", CheckMicrosatelliteRepeatType, AddRepeatTypeTandem},
  { "If D-loop or control region misc_feat is present, source must be mitochondrial", "DISC_MITOCHONDRION_REQUIRED", CheckMitochondrionRequired, MakeLocationMitochondrial},
  { "Unpublished pubs should have titles", "DISC_UNPUB_PUB_WITHOUT_TITLE", FindUnpubPubsWithoutTitles, NULL},
  { "Check for Quality Scores", "DISC_QUALITY_SCORES", CheckForQualityScores, NULL},
  { "rRNA product names should not contain 'internal', 'transcribed', or 'spacer'", "DISC_INTERNAL_TRANSCRIBED_SPACER_RRNA", InternalTranscribedSpacerrRNA, NULL},
  { "Find partial feature ends on sequences that could be extended", "DISC_PARTIAL_PROBLEMS", FindExtendablePartials, FixExtendablePartials},
  { "Find partial feature ends on bacterial sequences that cannot be extended", "DISC_BACTERIAL_PARTIAL_NONEXTENDABLE_PROBLEMS", FindBacterialNonExtendablePartials, FixBacterialNonExtendablePartials},
  { "Find partial feature ends on bacterial sequences that cannot be extended but have exceptions", "DISC_BACTERIAL_PARTIAL_NONEXTENDABLE_EXCEPTION", FindBacterialNonExtendablePartialsWithExceptions, NULL},
  { "rRNA product names should not contain 'partial' or 'domain'", "DISC_SUSPECT_RRNA_PRODUCTS", FindSuspectrRNAProducts, NULL},
  { "suspect misc_feature comments", "DISC_SUSPECT_MISC_FEATURES", FindBadMiscFeatures, NULL},
  { "Missing strain on bacterial 'Genus sp. strain'", "DISC_BACTERIA_MISSING_STRAIN", FindMissingBacteriaStrain, NULL},
  { "Missing definition lines", "DISC_MISSING_DEFLINES", FindMissingDefinitionLines, NULL},
  { "Missing affiliation", "DISC_MISSING_AFFIL", FindMissingAffiliations, NULL},
  { "Bacterial sources should not have isolate", "DISC_BACTERIA_SHOULD_NOT_HAVE_ISOLATE", FindBacteriaIsolate, NULL},
  { "Bacterial sequences should not have mRNA features", "DISC_BACTERIA_SHOULD_NOT_HAVE_MRNA", FindBacteriamRNA, NULL},
  { "Coding region has new exception", "DISC_CDS_HAS_NEW_EXCEPTION", FindCDSNewException, NULL},
  { "Trinomial sources should have corresponding qualifier", "DISC_TRINOMIAL_SHOULD_HAVE_QUALIFIER", FindTrinomialWithoutQualifier, NULL},
  { "Source has metagenomic qualifier", "DISC_METAGENOMIC", FindMetagenomic, NULL},
  { "Source has metagenome_source qualifier", "DISC_METAGENOME_SOURCE", FindMetagenomeSource, NULL},
  { "Missing genes", "ONCALLER_GENE_MISSING", OnCallerMissingAndSuperfluousGenes, NULL},
  { "Superfluous genes", "ONCALLER_SUPERFLUOUS_GENE", OnCallerMissingAndSuperfluousGenes, NULL},
  { "Short rRNA Features", "DISC_SHORT_RRNA", FindShortrRNAs, NULL},
  { "Authority and Taxname should match first two words", "ONCALLER_CHECK_AUTHORITY", CheckAuthorityTaxnameConflict, NULL},
  { "Submitter blocks and publications have consortiums", "ONCALLER_CONSORTIUM", FindConsortiums, RemoveConsortiums},
  { "Strain and culture-collection values conflict", "ONCALLER_STRAIN_CULTURE_COLLECTION_MISMATCH", FindStrainCultureCollectionMismatch, NULL},
  { "Comma or semicolon appears in strain or isolate", "ONCALLER_MULTISRC", FindMultiSrc, NULL} ,
  { "Multiple culture-collection quals", "ONCALLER_MULTIPLE_CULTURE_COLLECTION", FindMultipleCultureCollection, NULL},
  { "Segsets present", "DISC_SEGSETS_PRESENT", FindSegSets, NULL},
  { "Eco, mut, phy or pop sets present", "DISC_NONWGS_SETS_PRESENT", FindNonWGSSets, FixNonWGSSets},
  { "Feature List", "DISC_FEATURE_LIST", GetFeatureList, NULL},
  { "Category Header", "DISC_CATEGORY_HEADER", NULL, NULL},
  { "Mismatched Comments", "DISC_MISMATCHED_COMMENTS", FindMismatchedComments, FixMismatchedComments},
  { "BioSources with the same strain should have the same taxname", "DISC_STRAIN_TAXNAME_MISMATCH", CollectStrainTaxnameDiscrepancies, NULL},
  { "'Human' in host should be 'Homo sapiens'", "DISC_HUMAN_HOST", FindHumanHosts, FixHumanHosts},
  { "Genes on bacterial sequences should start with lowercase letters", "DISC_BAD_BACTERIAL_GENE_NAME", FindBadGeneNames, MoveBadGeneNames},
  { "Bad gene names", "TEST_BAD_GENE_NAME", FindBadGeneNames, MoveBadGeneNames },
  { "Location is ordered (intervals interspersed with gaps)", "ONCALLER_ORDERED_LOCATION", FindOrderedLocations, FixOrderedLocations},
  { "Comment descriptor present", "ONCALLER_COMMENT_PRESENT", FindCommentDescriptors, NULL },
  { "Titles on sets", "ONCALLER_DEFLINE_ON_SET", FindTitlesOnSets, NULL },
  { "HIV RNA location or molecule type inconsistent", "ONCALLER_HIV_RNA_INCONSISTENT", FindInconsistentHIVRNA, NULL },
  { "Protein sequences should be at least 50 aa, unless they are partial", "SHORT_PROT_SEQUENCES", FindShortProtSequences, NULL },
  { "mRNA sequences should not have exons", "TEST_EXON_ON_MRNA", FindExonsOnMrna, RemoveExonsOnMrna },
  { "Sequences with project IDs", "TEST_HAS_PROJECT_ID", FindProjectIdSequences, NULL },
  { "Feature has standard_name qualifier", "ONCALLER_HAS_STANDARD_NAME", FindStandardName, NULL },
  { "Missing structured comments", "ONCALLER_MISSING_STRUCTURED_COMMENTS", FindMissingStructuredComments, NULL },
  { "Bacteria should have strain", "DISC_REQUIRED_STRAIN", FindRequiredStrains, NULL},
  { "Bioseqs should have GenomeAssembly structured comments", "MISSING_GENOMEASSEMBLY_COMMENTS", FindMissingGenomeAssemblyStructuredComments, NULL },
  { "Bacterial taxnames should end with strain", "DISC_BACTERIAL_TAX_STRAIN_MISMATCH", FindBacterialTaxStrainMismatch, NULL },
  { "CDS has CDD Xref", "TEST_CDS_HAS_CDD_XREF", FindCDSWithCDDXref, NULL },
  { "Sequence contains unusual nucleotides", "TEST_UNUSUAL_NT", FindUnusualNT, NULL },
  { "Sequence contains regions of low quality", "TEST_LOW_QUALITY_REGION", FindLowQualityRegions, NULL },
  { "Organelle location should have genomic moltype", "TEST_ORGANELLE_NOT_GENOMIC", FindOrganelleNotGenomic, NULL },
  { "Intergenic spacer without plastid location", "TEST_UNWANTED_SPACER", FindUnwantedSpacers, NULL },
  { "Organelle products on non-organelle sequence", "TEST_ORGANELLE_PRODUCTS", FindOrganelleProducts, NULL },
  { "Organism ending in sp. needs tax consult", "TEST_SP_NOT_UNCULTURED", FindSpNotUncultured, NULL },
  { "mRNA sequence contains rearranged or germline", "TEST_BAD_MRNA_QUAL", FindBadMrnaQual, NULL },
  { "Unnecessary environmental qualifier present", "TEST_UNNECESSARY_ENVIRONMENTAL", FindUnnecessaryEnvironmental, NULL },
  { "Unnecessary gene features on virus", "TEST_UNNECESSARY_VIRUS_GENE", FindUnnecessaryVirusGene, NULL },
  { "Set wrapper on microsatellites or rearranged genes", "TEST_UNWANTED_SET_WRAPPER", FindUnwantedSetWrappers, NULL},
  { "Missing values in primer set", "TEST_MISSING_PRIMER", FindMissingPrimerValues, NULL},
  { "Unexpected misc_RNA features", "TEST_UNUSUAL_MISC_RNA", FindUnexpectedMiscRNA, NULL},
  { "Species-specific primers, no environmental sample", "TEST_AMPLIFIED_PRIMERS_NO_ENVIRONMENTAL_SAMPLE", FindAmpPrimersNoEnvSample, NULL},
  { "Duplicate genes on opposite strands", "TEST_DUP_GENES_OPPOSITE_STRANDS", FindDuplicateGenesOnOppositeStrands, NULL},
  { "Problems with small genome sets", "TEST_SMALL_GENOME_SET_PROBLEM", FindSmallGenomeSetProblems, NULL},
  { "Overlapping rRNA features", "TEST_OVERLAPPING_RRNAS", AddOverlappingrRNADiscrepancies, NULL},
  { "mRNA sequences have CDS/gene on the complement strand", "TEST_MRNA_SEQUENCE_MINUS_STRAND_FEATURES", FindMrnaSequencesWithMinusStrandFeatures, NULL},
  { "Complete taxname should be present in definition line", "TEST_TAXNAME_NOT_IN_DEFLINE", FindTaxnameMissingFromDefline, NULL},
  { "Count number of unverified sequences", "TEST_COUNT_UNVERIFIED", CountUnverifiedSequences, NULL},
  { "Show translation exception", "SHOW_TRANSL_EXCEPT", ShowTranslExcept, NULL},
  { "Show hypothetic protein having a gene name", "SHOW_HYPOTHETICAL_CDS_HAVING_GENE_NAME", ShowCDsHavingGene, NULL},
  { "Test defline existence", "TEST_DEFLINE_PRESENT", TestDeflineExistence, NULL},
  { "Remove mRNA overlapping a pseudogene", "TEST_MRNA_OVERLAPPING_PSEUDO_GENE", TestMrnaOverlappingPseudoGene, RmvMrnaOverlappingPseudoGene}
};


extern Boolean IsTestTypeAppropriateForReportType (Int4 test_type, EDiscrepancyReportType report_type)
{
  Boolean rval = FALSE;

  switch (report_type) {
    case eReportTypeDiscrepancy:
      if (test_type == DISC_SOURCE_QUALS_ASNDISC
          || test_type == TEST_MRNA_OVERLAPPING_PSEUDO_GENE
          || test_type == DISC_MISSING_VIRAL_QUALS
          || test_type == DISC_MISSING_SRC_QUAL
          || test_type == DISC_DUP_SRC_QUAL
          || test_type == DISC_DUP_SRC_QUAL_DATA
          || test_type == DISC_HAPLOTYPE_MISMATCH
          || test_type == DISC_FEATURE_MOLTYPE_MISMATCH
          || test_type == DISC_CDS_WITHOUT_MRNA
          || test_type == DISC_EXON_INTRON_CONFLICT
          || test_type == DISC_FEATURE_COUNT
          || test_type == DISC_SPECVOUCHER_TAXNAME_MISMATCH
          || test_type == DISC_GENE_PARTIAL_CONFLICT
          || test_type == DISC_FLATFILE_FIND_ONCALLER
          || test_type == DISC_CDS_PRODUCT_FIND
          || test_type == DISC_DUP_DEFLINE
          || test_type == DISC_COUNT_NUCLEOTIDES
          || test_type == DUP_DISC_ATCC_CULTURE_CONFLICT
          || test_type == DISC_USA_STATE
          || test_type == DISC_INCONSISTENT_MOLTYPES
          || test_type == DISC_SRC_QUAL_PROBLEM
          || test_type == DISC_SUBMITBLOCK_CONFLICT
          || test_type == DISC_POSSIBLE_LINKER
          || test_type == DISC_TITLE_AUTHOR_CONFLICT
          || test_type == DISC_BAD_GENE_STRAND
          || test_type == DISC_MAP_CHROMOSOME_CONFLICT
          || test_type == DISC_RBS_WITHOUT_GENE
          || test_type == DISC_CITSUBAFFIL_CONFLICT
          || test_type == DISC_REQUIRED_CLONE
          || test_type == DISC_SUSPICIOUS_NOTE_TEXT
          || test_type == DISC_mRNA_ON_WRONG_SEQUENCE_TYPE
          || test_type == DISC_RETROVIRIDAE_DNA
          || test_type == DISC_CHECK_AUTH_CAPS
          || test_type == DISC_CHECK_RNA_PRODUCTS_AND_COMMENTS
          || test_type == DISC_MICROSATELLITE_REPEAT_TYPE
          || test_type == DISC_MITOCHONDRION_REQUIRED
          || test_type == DISC_UNPUB_PUB_WITHOUT_TITLE
          || test_type == DISC_INTERNAL_TRANSCRIBED_SPACER_RRNA
          || test_type == DISC_BACTERIA_MISSING_STRAIN
          || test_type == DISC_MISSING_DEFLINES
          || test_type == DISC_MISSING_AFFIL
          || test_type == DISC_BACTERIA_SHOULD_NOT_HAVE_ISOLATE
          || test_type == DISC_BACTERIA_SHOULD_NOT_HAVE_MRNA
          || test_type == DISC_CDS_HAS_NEW_EXCEPTION
          || test_type == DISC_TRINOMIAL_SHOULD_HAVE_QUALIFIER
          || test_type == DISC_METAGENOMIC
          || test_type == DISC_METAGENOME_SOURCE
          || test_type == ONCALLER_GENE_MISSING
          || test_type == ONCALLER_SUPERFLUOUS_GENE
          || test_type == ONCALLER_CHECK_AUTHORITY
          || test_type == ONCALLER_CONSORTIUM
          || test_type == ONCALLER_STRAIN_CULTURE_COLLECTION_MISMATCH
          || test_type == ONCALLER_MULTISRC
          || test_type == ONCALLER_MULTIPLE_CULTURE_COLLECTION
          || test_type == DISC_STRAIN_TAXNAME_MISMATCH
          || test_type == DISC_HUMAN_HOST
          || test_type == ONCALLER_ORDERED_LOCATION
          || test_type == ONCALLER_COMMENT_PRESENT
          || test_type == ONCALLER_DEFLINE_ON_SET
          || test_type == ONCALLER_HIV_RNA_INCONSISTENT
          || test_type == TEST_EXON_ON_MRNA
          || test_type == TEST_HAS_PROJECT_ID
          || test_type == ONCALLER_HAS_STANDARD_NAME
          || test_type == ONCALLER_MISSING_STRUCTURED_COMMENTS
          || test_type == TEST_ORGANELLE_PRODUCTS
          || test_type == TEST_SP_NOT_UNCULTURED
          || test_type == TEST_BAD_MRNA_QUAL
          || test_type == TEST_UNNECESSARY_ENVIRONMENTAL
          || test_type == TEST_UNNECESSARY_VIRUS_GENE
          || test_type == TEST_UNWANTED_SET_WRAPPER
          || test_type == TEST_MISSING_PRIMER
          || test_type == TEST_AMPLIFIED_PRIMERS_NO_ENVIRONMENTAL_SAMPLE
          || test_type == TEST_SMALL_GENOME_SET_PROBLEM
          || test_type == TEST_MRNA_SEQUENCE_MINUS_STRAND_FEATURES
          || test_type == TEST_TAXNAME_NOT_IN_DEFLINE
          || test_type == TEST_COUNT_UNVERIFIED) {
        rval = FALSE;
      } else {
        rval = TRUE;
      }
      break;
    case eReportTypeOnCaller:
      if (test_type == DISC_RNA_NO_PRODUCT 
          || test_type == TEST_MRNA_OVERLAPPING_PSEUDO_GENE
          || test_type == DISC_BADLEN_TRNA
          || test_type == DISC_MISSING_VIRAL_QUALS
          || test_type == DISC_MISSING_SRC_QUAL
          || test_type == DISC_DUP_SRC_QUAL
          || test_type == DISC_DUP_SRC_QUAL_DATA
          || test_type == DISC_NON_GENE_LOCUS_TAG
          || test_type == DISC_PSEUDO_MISMATCH
          || test_type == DISC_SHORT_INTRON
          || test_type == DISC_INFLUENZA_DATE_MISMATCH
          || test_type == DISC_HAPLOTYPE_MISMATCH
          || test_type == DISC_FEATURE_MOLTYPE_MISMATCH
          || test_type == DISC_CDS_WITHOUT_MRNA
          || test_type == DISC_EXON_INTRON_CONFLICT
          || test_type == DISC_FEATURE_COUNT
          || test_type == DISC_SPECVOUCHER_TAXNAME_MISMATCH
          || test_type == DISC_GENE_PARTIAL_CONFLICT
          || test_type == DISC_FLATFILE_FIND_ONCALLER
          || test_type == DISC_CDS_PRODUCT_FIND
          || test_type == DISC_DUP_DEFLINE
          || test_type == DISC_COUNT_NUCLEOTIDES
          || test_type == DUP_DISC_ATCC_CULTURE_CONFLICT
          || test_type == DISC_USA_STATE
          || test_type == DISC_INCONSISTENT_MOLTYPES
          || test_type == DISC_SRC_QUAL_PROBLEM
          || test_type == DISC_SUBMITBLOCK_CONFLICT
          || test_type == DISC_POSSIBLE_LINKER
          || test_type == DISC_TITLE_AUTHOR_CONFLICT
          || test_type == DISC_BAD_GENE_STRAND
          || test_type == DISC_MAP_CHROMOSOME_CONFLICT
          || test_type == DISC_RBS_WITHOUT_GENE
          || test_type == DISC_CITSUBAFFIL_CONFLICT
          || test_type == DISC_REQUIRED_CLONE
          || test_type == DISC_SUSPICIOUS_NOTE_TEXT
          || test_type == DISC_mRNA_ON_WRONG_SEQUENCE_TYPE
          || test_type == DISC_RETROVIRIDAE_DNA
          || test_type == DISC_CHECK_AUTH_CAPS
          || test_type == DISC_CHECK_RNA_PRODUCTS_AND_COMMENTS
          || test_type == DISC_MICROSATELLITE_REPEAT_TYPE
          || test_type == DISC_MITOCHONDRION_REQUIRED
          || test_type == DISC_UNPUB_PUB_WITHOUT_TITLE
          || test_type == DISC_INTERNAL_TRANSCRIBED_SPACER_RRNA
          || test_type == DISC_BACTERIA_MISSING_STRAIN
          || test_type == DISC_MISSING_DEFLINES
          || test_type == DISC_MISSING_AFFIL
          || test_type == DISC_BACTERIA_SHOULD_NOT_HAVE_ISOLATE
          || test_type == DISC_BACTERIA_SHOULD_NOT_HAVE_MRNA
          || test_type == DISC_CDS_HAS_NEW_EXCEPTION
          || test_type == DISC_TRINOMIAL_SHOULD_HAVE_QUALIFIER
          || test_type == DISC_METAGENOMIC
          || test_type == DISC_METAGENOME_SOURCE
          || test_type == ONCALLER_GENE_MISSING
          || test_type == ONCALLER_SUPERFLUOUS_GENE
          || test_type == ONCALLER_CHECK_AUTHORITY
          || test_type == ONCALLER_CONSORTIUM
          || test_type == ONCALLER_STRAIN_CULTURE_COLLECTION_MISMATCH
          || test_type == ONCALLER_MULTISRC
          || test_type == ONCALLER_MULTIPLE_CULTURE_COLLECTION
          || test_type == DISC_STRAIN_TAXNAME_MISMATCH
          || test_type == DISC_HUMAN_HOST
          || test_type == ONCALLER_ORDERED_LOCATION
          || test_type == ONCALLER_COMMENT_PRESENT
          || test_type == ONCALLER_DEFLINE_ON_SET
          || test_type == ONCALLER_HIV_RNA_INCONSISTENT
          || test_type == TEST_EXON_ON_MRNA
          || test_type == TEST_HAS_PROJECT_ID
          || test_type == ONCALLER_HAS_STANDARD_NAME
          || test_type == ONCALLER_MISSING_STRUCTURED_COMMENTS
          || test_type == TEST_ORGANELLE_NOT_GENOMIC
          || test_type == TEST_UNWANTED_SPACER
          || test_type == TEST_ORGANELLE_PRODUCTS
          || test_type == TEST_SP_NOT_UNCULTURED
          || test_type == TEST_BAD_MRNA_QUAL
          || test_type == TEST_UNNECESSARY_ENVIRONMENTAL
          || test_type == TEST_UNNECESSARY_VIRUS_GENE
          || test_type == TEST_UNWANTED_SET_WRAPPER
          || test_type == TEST_MISSING_PRIMER
          || test_type == TEST_AMPLIFIED_PRIMERS_NO_ENVIRONMENTAL_SAMPLE
          || test_type == TEST_SMALL_GENOME_SET_PROBLEM
          || test_type == TEST_MRNA_SEQUENCE_MINUS_STRAND_FEATURES
          || test_type == TEST_TAXNAME_NOT_IN_DEFLINE
          || test_type == TEST_COUNT_UNVERIFIED
          || test_type == DISC_SHORT_RRNA) {
        rval = TRUE;
      }
      break;
    case eReportTypeMegaReport:
      rval = TRUE;
      break;
  }
  return rval;
}


extern void PrintDiscrepancyTestList (FILE *fp)
{
  Int4 i;
  CharPtr tmp;

  /* discrepancy report */
  fprintf (fp, "Discrepancy Report Tests\n");
  for (i = 0; i < MAX_DISC_TYPE; i++) {
    if (IsTestTypeAppropriateForReportType (i, eReportTypeDiscrepancy)) {
      fprintf (fp, "%s  %s  %s\n", discrepancy_info_list[i].setting_name, 
                                   discrepancy_info_list[i].conf_name,
                                   discrepancy_info_list[i].autofix_func == NULL ? "" : "Has Autofix");
    }
  }
  fprintf (fp, "\n");

  /* on-caller tool */
  fprintf (fp, "On-Caller Tool Tests\n");
  for (i = 0; i < MAX_DISC_TYPE; i++) {
    if (IsTestTypeAppropriateForReportType (i, eReportTypeOnCaller)) {
      fprintf (fp, "%s  %s  %s\n", discrepancy_info_list[i].setting_name, 
                                   discrepancy_info_list[i].conf_name,
                                   discrepancy_info_list[i].autofix_func == NULL ? "" : "Has Autofix");
    }
  }
  fprintf (fp, "\n");

  fprintf (fp, "Terms searched for by SUSPECT_PRODUCT_NAMES:\n");
  for (i = 0; i < num_suspect_product_terms; i++) {
    fprintf (fp, "'%s':%s (Category: %s)\n", 
             suspect_product_terms[i].pattern, 
             SummarizeSuspectPhraseFunc(suspect_product_terms[i].search_func),
             suspect_name_category_names[suspect_product_terms[i].fix_type]);
  }
  fprintf (fp, "\n");

  fprintf (fp, "Replacements for SUSPECT_PRODUCT_NAMES:\n");
  fprintf (fp, "Typos:\n");
  for (i = 0; i < num_suspect_product_terms; i++) {
    if (suspect_product_terms[i].replace_func != NULL && suspect_product_terms[i].fix_type == eSuspectNameType_Typo) {
      tmp = SummarizeSuspectReplacementPhrase (suspect_product_terms[i].replace_func, suspect_product_terms[i].replace_phrase);
      fprintf (fp, "'%s':%s (%s)\n", 
               suspect_product_terms[i].pattern, 
               SummarizeSuspectPhraseFunc(suspect_product_terms[i].search_func),
               tmp);
      tmp = MemFree (tmp);
    }
  }
  fprintf (fp, "QuickFixes:\n");
  for (i = 0; i < num_suspect_product_terms; i++) {
    if (suspect_product_terms[i].replace_func != NULL && suspect_product_terms[i].fix_type == eSuspectNameType_QuickFix) {
      tmp = SummarizeSuspectReplacementPhrase (suspect_product_terms[i].replace_func, suspect_product_terms[i].replace_phrase);
      fprintf (fp, "'%s':%s (%s)\n", 
               suspect_product_terms[i].pattern, 
               SummarizeSuspectPhraseFunc(suspect_product_terms[i].search_func),
               tmp);
      tmp = MemFree (tmp);
    }
  }
  fprintf (fp, "\n");

  fprintf (fp, "Terms searched for by SUSPECT_PHRASES:\n");
  for (i = 0; i < num_suspect_phrases; i++) {
    fprintf (fp, "%s\n", suspect_phrases[i]);
  }
  fprintf (fp, "\n");

  fprintf (fp, "Terms searched for by DISC_SUSPICIOUS_NOTE_TEXT:\n");
  for (i = 0; i < num_suspicious_note_phrases; i++) {
    fprintf (fp, "%s\n", suspicious_note_phrases[i]);
  }
  fprintf (fp, "\n");

  fprintf (fp, "Terms searched for by DISC_FLATFILE_FIND_ONCALLER:\n");
  for (i = 0; flatfile_find_list_oncaller[i] != NULL; i++) {
    fprintf (fp, "%s\n", flatfile_find_list_oncaller[i]);
  }
  fprintf (fp, "\n");

  fprintf (fp, "Terms searched for by DISC_CDS_PRODUCT_FIND:\n");
  for (i = 0; i < num_cds_product_find; i++) {
    fprintf (fp, "'%s':%s\n",
             cds_product_find[i].pattern,
             SummarizeSuspectPhraseFunc(cds_product_find[i].search_func));
  }
  fprintf (fp, "\n");
  
}


extern CharPtr GetDiscrepancyTestConfName (DiscrepancyType dtype) 
{
  return discrepancy_info_list[dtype].conf_name;
}

extern CharPtr GetDiscrepancyTestSettingName (DiscrepancyType dtype) 
{
  return discrepancy_info_list[dtype].setting_name;
}

extern DiscrepancyType GetDiscrepancyTypeFromSettingName (CharPtr setting_name)
{
  Int4 i;

  if (StringHasNoText (setting_name)) {
    return MAX_DISC_TYPE;
  }
  for (i = 0; i < MAX_DISC_TYPE; i++) {
    if (StringICmp (setting_name, discrepancy_info_list[i].setting_name) == 0) {
      return (DiscrepancyType) i;
    }
  }
  return MAX_DISC_TYPE;
}


extern void ConfigureForBigSequence (DiscrepancyConfigPtr dcp)
{
  Int4 i;

  if (dcp == NULL) {
    return;
  }

  for (i = 0; i < MAX_DISC_TYPE; i++) {
    dcp->conf_list[i] = FALSE;
  }
  dcp->conf_list[DISC_SHORT_CONTIG] = TRUE;
  dcp->conf_list[DISC_INCONSISTENT_BIOSRC] = TRUE;
  dcp->conf_list[DISC_SHORT_SEQUENCE] = TRUE;
  dcp->conf_list[DISC_PERCENTN] = TRUE;
  dcp->conf_list[DISC_N_RUNS] = TRUE;
  dcp->conf_list[DISC_ZERO_BASECOUNT] = TRUE;
  dcp->conf_list[DISC_NO_ANNOTATION] = TRUE;
  dcp->conf_list[DISC_COUNT_NUCLEOTIDES] = TRUE;
  dcp->conf_list[DISC_QUALITY_SCORES] = TRUE;
  dcp->conf_list[TEST_DEFLINE_PRESENT] = TRUE;
  dcp->conf_list[DISC_INCONSISTENT_BIOSRC_DEFLINE] = TRUE;
  dcp->conf_list[DISC_SHORT_SEQUENCE] = TRUE;
  dcp->conf_list[DISC_PERCENTN] = TRUE;
  dcp->conf_list[DISC_SRC_QUAL_PROBLEM] = TRUE;
  dcp->conf_list[DISC_MISSING_SRC_QUAL] = TRUE;
  dcp->conf_list[DISC_DUP_SRC_QUAL] = TRUE;
  dcp->conf_list[DISC_DUP_SRC_QUAL_DATA] = TRUE;
  dcp->conf_list[DISC_HAPLOTYPE_MISMATCH] = TRUE;
  dcp->conf_list[DISC_SPECVOUCHER_TAXNAME_MISMATCH] = TRUE;
  dcp->conf_list[DUP_DISC_ATCC_CULTURE_CONFLICT] = TRUE;
  dcp->conf_list[DISC_USA_STATE] = TRUE;
  dcp->conf_list[DISC_INCONSISTENT_MOLTYPES] = TRUE;
  dcp->conf_list[DISC_SUBMITBLOCK_CONFLICT] = TRUE;
  dcp->conf_list[DISC_TITLE_AUTHOR_CONFLICT] = TRUE;
  dcp->conf_list[DISC_MAP_CHROMOSOME_CONFLICT] = TRUE;
  dcp->conf_list[DISC_CITSUBAFFIL_CONFLICT] = TRUE;
  dcp->conf_list[DISC_REQUIRED_CLONE] = TRUE;
  dcp->conf_list[DISC_SOURCE_QUALS_ASNDISC] = TRUE;
  dcp->conf_list[DISC_CHECK_AUTH_CAPS] = TRUE;
  dcp->conf_list[DISC_UNPUB_PUB_WITHOUT_TITLE] = TRUE;
  dcp->conf_list[DISC_QUALITY_SCORES] = TRUE;
  dcp->conf_list[DISC_MISSING_DEFLINES] = TRUE;
  dcp->conf_list[DISC_MISSING_AFFIL] = TRUE;
  dcp->conf_list[DISC_BACTERIA_SHOULD_NOT_HAVE_ISOLATE] = TRUE;
  dcp->conf_list[DISC_TRINOMIAL_SHOULD_HAVE_QUALIFIER] = TRUE;
  dcp->conf_list[ONCALLER_CHECK_AUTHORITY] = TRUE;
  dcp->conf_list[ONCALLER_CONSORTIUM] = TRUE;
  dcp->conf_list[ONCALLER_STRAIN_CULTURE_COLLECTION_MISMATCH] = TRUE;
  dcp->conf_list[ONCALLER_MULTISRC] = TRUE;
  dcp->conf_list[ONCALLER_MULTIPLE_CULTURE_COLLECTION] = TRUE;
  dcp->conf_list[DISC_SEGSETS_PRESENT] = TRUE;
  dcp->conf_list[DISC_NONWGS_SETS_PRESENT] = TRUE;
  dcp->conf_list[DISC_MISMATCHED_COMMENTS] = TRUE;
  dcp->conf_list[DISC_STRAIN_TAXNAME_MISMATCH] = TRUE;
  dcp->conf_list[DISC_HUMAN_HOST] = TRUE;
  dcp->conf_list[ONCALLER_COMMENT_PRESENT] = TRUE;
  dcp->conf_list[ONCALLER_DEFLINE_ON_SET] = TRUE;
  dcp->conf_list[TEST_HAS_PROJECT_ID] = TRUE;
  dcp->conf_list[ONCALLER_MISSING_STRUCTURED_COMMENTS] = TRUE;
  dcp->conf_list[DISC_REQUIRED_STRAIN] = TRUE;
  dcp->conf_list[MISSING_GENOMEASSEMBLY_COMMENTS] = TRUE;
  dcp->conf_list[DISC_BACTERIAL_TAX_STRAIN_MISMATCH] = TRUE;
  dcp->conf_list[TEST_UNUSUAL_NT] = TRUE;
  dcp->conf_list[TEST_LOW_QUALITY_REGION] = TRUE;
  dcp->conf_list[TEST_SP_NOT_UNCULTURED] = TRUE;
  dcp->conf_list[TEST_UNNECESSARY_ENVIRONMENTAL] = TRUE;
  dcp->conf_list[TEST_UNWANTED_SET_WRAPPER] = TRUE;
  dcp->conf_list[TEST_AMPLIFIED_PRIMERS_NO_ENVIRONMENTAL_SAMPLE] = TRUE;
}


extern void ConfigureForGenomes (DiscrepancyConfigPtr dcp)
{
  Int4 i;

  if (dcp == NULL) {
    return;
  }

  for (i = 0; i < MAX_DISC_TYPE; i++) {
    dcp->conf_list[i] = TRUE;
  }
  dcp->conf_list[DISC_STRAIN_TAXNAME_MISMATCH] = FALSE;
  dcp->conf_list[DISC_CITSUBAFFIL_CONFLICT] = FALSE;
  dcp->conf_list[DISC_OVERLAPPING_GENES] = FALSE;
  dcp->conf_list[DISC_INCONSISTENT_BIOSRC_DEFLINE] = FALSE;
  dcp->conf_list[DISC_NO_TAXLOOKUP] = FALSE;
  dcp->conf_list[DISC_BAD_TAXLOOKUP] = FALSE;
  dcp->conf_list[DISC_COUNT_TRNA] = FALSE;
  dcp->conf_list[DISC_BADLEN_TRNA] = FALSE;
  dcp->conf_list[DISC_STRAND_TRNA] = FALSE;
  dcp->conf_list[DISC_COUNT_RRNA] = FALSE;
  dcp->conf_list[DISC_CDS_OVERLAP_TRNA] = FALSE;
  dcp->conf_list[DISC_FEAT_OVERLAP_SRCFEAT] = FALSE;
  dcp->conf_list[DISC_INFLUENZA_DATE_MISMATCH] = FALSE;   
  dcp->conf_list[DISC_MISSING_VIRAL_QUALS] = FALSE;          
  dcp->conf_list[DISC_SRC_QUAL_PROBLEM] = FALSE;
  dcp->conf_list[DISC_MISSING_SRC_QUAL] = FALSE;
  dcp->conf_list[DISC_DUP_SRC_QUAL] = FALSE;
  dcp->conf_list[DISC_HAPLOTYPE_MISMATCH] = FALSE;
  dcp->conf_list[DISC_FEATURE_MOLTYPE_MISMATCH] = FALSE;
  dcp->conf_list[DISC_SPECVOUCHER_TAXNAME_MISMATCH] = FALSE;
  dcp->conf_list[DISC_FLATFILE_FIND_ONCALLER] = FALSE;
  dcp->conf_list[DISC_CDS_PRODUCT_FIND] = FALSE;
  dcp->conf_list[DISC_DUP_DEFLINE] = FALSE;
  dcp->conf_list[DISC_INCONSISTENT_MOLTYPES] = FALSE;
  dcp->conf_list[DISC_SUBMITBLOCK_CONFLICT] = FALSE;
  dcp->conf_list[DISC_POSSIBLE_LINKER] = FALSE;
  dcp->conf_list[DISC_TITLE_AUTHOR_CONFLICT] = FALSE;
  dcp->conf_list[DISC_MAP_CHROMOSOME_CONFLICT] = FALSE;
  dcp->conf_list[DISC_REQUIRED_CLONE] = FALSE;
  dcp->conf_list[DISC_mRNA_ON_WRONG_SEQUENCE_TYPE] = FALSE;
  dcp->conf_list[DISC_RETROVIRIDAE_DNA] = FALSE;
  dcp->conf_list[DISC_MISSING_DEFLINES] = FALSE;
  dcp->conf_list[ONCALLER_GENE_MISSING] = FALSE;
  dcp->conf_list[ONCALLER_SUPERFLUOUS_GENE] = FALSE;
  dcp->conf_list[ONCALLER_CONSORTIUM] = FALSE;
  dcp->conf_list[DISC_FEATURE_LIST] = FALSE;
  dcp->conf_list[TEST_ORGANELLE_PRODUCTS] = FALSE;

  /* mitochondrial tests */
  dcp->conf_list[DISC_DUP_TRNA] = FALSE;
  dcp->conf_list[DISC_DUP_RRNA] = FALSE;
  dcp->conf_list[DISC_TRANSL_NO_NOTE] = FALSE;
  dcp->conf_list[DISC_NOTE_NO_TRANSL] = FALSE;
  dcp->conf_list[DISC_TRANSL_TOO_LONG] = FALSE;
  dcp->conf_list[DISC_COUNT_PROTEINS] = FALSE;

  /* on-caller specific tests */
  dcp->conf_list[DISC_SRC_QUAL_PROBLEM] = FALSE;
  dcp->conf_list[DISC_CATEGORY_HEADER] = FALSE;
  dcp->conf_list[TEST_TAXNAME_NOT_IN_DEFLINE] = FALSE;
}


extern void ConfigureForReportType (DiscrepancyConfigPtr dcp, EDiscrepancyReportType report_type)
{
  Int4 i;

  if (dcp == NULL) {
    return;
  }

  for (i = 0; i < MAX_DISC_TYPE; i++) {
    dcp->conf_list[i] = IsTestTypeAppropriateForReportType (i, report_type);
  }
}


/* Note that this function contains a hack - it assumes that all of the
 * test types that use the same collection function are listed together.
 */
extern ValNodePtr CollectDiscrepancies (DiscrepancyConfigPtr dcp, ValNodePtr sep_list, PerformDiscrepancyTest taxlookup)
{
  ValNodePtr             discrepancy_list = NULL;
  Int4                   i;
  PerformDiscrepancyTest last_test_func = NULL;

  discrepancy_info_list[DISC_NO_TAXLOOKUP].test_func = taxlookup;
  discrepancy_info_list[DISC_BAD_TAXLOOKUP].test_func = taxlookup;

  for (i = 0; i < MAX_DISC_TYPE; i++)
  {
    if ((dcp == NULL || dcp->conf_list[i])
        && discrepancy_info_list[i].test_func != NULL
        && discrepancy_info_list[i].test_func != last_test_func)
    {
      discrepancy_info_list[i].test_func (&discrepancy_list, sep_list);
      last_test_func = discrepancy_info_list[i].test_func;
    }
  }
  
  /* because some tests are run together, need to remove unwanted results */
  RemoveUnwantedDiscrepancyItems (&discrepancy_list, dcp);

  /* normalize the discrepancy levels so that they will be correctly displayed */
  SetDiscrepancyLevels (discrepancy_list, 0);
  return discrepancy_list;  
}


extern void AutofixDiscrepancies (ValNodePtr vnp, Boolean fix_all, LogInfoPtr lip)
{
  ClickableItemPtr cip;

  while (vnp != NULL) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL) {
      if (cip->chosen || fix_all) {
        if (discrepancy_info_list[cip->clickable_item_type].autofix_func != NULL) {
          (discrepancy_info_list[cip->clickable_item_type].autofix_func) (cip->item_list, NULL, lip);
        }
        if (cip->autofix_func != NULL) {
          (cip->autofix_func)(cip->item_list, cip->autofix_data, lip);
        }
      }
      AutofixDiscrepancies (cip->subcategories, fix_all || cip->chosen, lip);
    }
    vnp = vnp->next;
  }
}

extern void ChooseFixableDiscrepancies (ValNodePtr vnp)
{
  ClickableItemPtr cip;

  while (vnp != NULL) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL && !cip->chosen) {
      if (discrepancy_info_list[cip->clickable_item_type].autofix_func != NULL
          || cip->autofix_func != NULL) {
        cip->chosen = TRUE;
      } else {
        ChooseFixableDiscrepancies (cip->subcategories);
      }
    }
    vnp = vnp->next;
  }
}


static CharPtr GetLocusTagForFeature (SeqFeatPtr sfp)
{
  GeneRefPtr grp = NULL;
  SeqFeatPtr gene;
  
  if (sfp == NULL) {
    return NULL;
  }
  if (sfp->data.choice == SEQFEAT_GENE) {
    grp = sfp->data.value.ptrvalue;
  } else {
    grp = SeqMgrGetGeneXref (sfp);
    if (grp == NULL) {
      gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
      if (gene != NULL) {
        grp = (GeneRefPtr) gene->data.value.ptrvalue;
      }
    }
  }

  if (grp == NULL) {
    return NULL;
  } else {
    return grp->locus_tag;
  }
}


extern CharPtr GetBioseqLabel (BioseqPtr bsp)
{
  Char        id_str[45];

  if (bsp == NULL) {
    return NULL;
  }
  
  SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_str, PRINTID_REPORT, 39);
  return StringSave (id_str);
}


extern CharPtr GetBioseqSetLabel (BioseqSetPtr bssp)
{
  Char        id_str[45];
  CharPtr     tmp, set_fmt = "Set containing %s", id_label;
  BioseqPtr   bsp;

  if (bssp == NULL) {
    return NULL;
  }
  if (bssp->_class == BioseqseqSet_class_segset) {
    sprintf (id_str, "ss|");
  } else if (bssp->_class == BioseqseqSet_class_nuc_prot) {
    sprintf (id_str, "np|");
  } else if (bssp->seq_set != NULL && bssp->seq_set->data.ptrvalue != NULL && IS_Bioseq (bssp->seq_set)) {
    bsp = bssp->seq_set->data.ptrvalue;
    SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_str, PRINTID_REPORT, 39);
    tmp = MemNew (sizeof (Char) * (StringLen (set_fmt) + StringLen (id_str)));
    sprintf (tmp, set_fmt, id_str);
    return tmp;
  } else if (bssp->seq_set != NULL && bssp->seq_set->data.ptrvalue != NULL && IS_Bioseq_set (bssp->seq_set)) {
    id_label = GetBioseqSetLabel (bssp->seq_set->data.ptrvalue);
    tmp = MemNew (sizeof (Char) * (StringLen (set_fmt) + StringLen (id_label)));
    sprintf (tmp, set_fmt, id_label);
    id_label = MemFree (id_label);
    return tmp;
  } else {
    return StringSave ("BioseqSet");
  }
  bsp = GetRepresentativeBioseqFromBioseqSet (bssp);
  if (bsp != NULL) {
    SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_str + 3, PRINTID_REPORT, 39);
    return StringSave (id_str);
  }
  return NULL;
}


static void LIBCALLBACK CountNonATGCNTProc (CharPtr sequence, Pointer userdata)
{
  Int4Ptr p_i;
  CharPtr cp;

  if (sequence == NULL || userdata == NULL) return;
  p_i = (Int4Ptr) userdata;

  for (cp = sequence; *cp != 0; cp++)
  {
    if (*cp != 'A' && *cp != 'T' && *cp != 'G' && *cp != 'C')
    {
      (*p_i) ++;
    }
  }
}


extern CharPtr GetDiscrepancyItemTextEx (ValNodePtr vnp, CharPtr filename)
{
  CharPtr           row_text = NULL, tmp, fmt = "%s:%s";
  SeqFeatPtr        sfp, cds, sfp_index = NULL;
  BioseqPtr         bsp;
  SeqMgrFeatContext context;
  CharPtr           location;
  CharPtr           label;
  SeqDescrPtr       sdp;
  CharPtr           locus_tag = "";
  CharPtr           bsp_fmt = "%s (length %d)\n";
  CharPtr           bsp_unusual_fmt = "%s (length %d, %d other)\n";
  ObjValNodePtr     ovn;
  SeqEntryPtr       sep;
  SeqSubmitPtr      ssp;
  Boolean           special_flag = FALSE;
  Uint1             data_choice;
  ValNodePtr        extra_fields = NULL, field, field_strings = NULL, field_values, val_vnp;
  Int4              field_len = 0, label_len, num_bad;
  
  if (vnp == NULL)
  {
    return NULL;
  }

  if (vnp->extended > 0) {
    ovn = (ObjValNodePtr) vnp;
    extra_fields = ovn->idx.scratch;
  }

  data_choice = vnp->choice;
  if (data_choice > OBJ_MAX) {
    special_flag = TRUE;
    data_choice -= OBJ_MAX;
  }

  if (data_choice == OBJ_SEQFEAT)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    if (sfp != NULL)
    {
      if (SeqMgrFeaturesAreIndexed(sfp->idx.entityID) == 0) {
        SeqMgrIndexFeatures (sfp->idx.entityID, NULL);
      }

      sfp_index = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &context);
      if (sfp_index != NULL && sfp_index->idx.subtype == FEATDEF_PROT) {
        bsp = BioseqFindFromSeqLoc (sfp_index->location);
        if (bsp != NULL) {
          cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
          if (cds != NULL) {
            sfp = cds;
          }
        }
      }
      if (sfp != NULL)
      {
        location = SeqLocPrintUseBestID (sfp->location);
        if (location == NULL) {
          location = StringSave ("Unknown location");
        }
        label = (CharPtr) FeatDefTypeLabel(sfp);
        if (label == NULL) {
          label = "Unknown label";
        }
        locus_tag = GetLocusTagForFeature (sfp);
        if (sfp_index == NULL) {
          context.label = "Unknown context label";
        } else if (context.label == NULL) {
          context.label = "Unknown context label";
        }

        row_text = (CharPtr) MemNew (sizeof (Char) * 
                                     (StringLen (label) 
                                      + StringLen (context.label) 
                                      + StringLen (location) 
                                      + StringLen (locus_tag)
                                      + 6));
        sprintf (row_text, "%s\t%s\t%s\t%s\n", label,
                                               context.label,
                                               location,
                                               locus_tag == NULL ? "" : locus_tag);
        location = MemFree (location);
      }
    }
  }
  else if (data_choice == OBJ_BIOSEQ)
  {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (bsp != NULL)
    {
      tmp = GetBioseqLabel (vnp->data.ptrvalue);
      num_bad = 0;
      SeqPortStream (bsp, 0, (Pointer) &num_bad, CountNonATGCNTProc);
      if (num_bad > 0) {
        row_text = (CharPtr) MemNew (sizeof(Char) * (StringLen (bsp_unusual_fmt) + StringLen (tmp) + 47));
        sprintf (row_text, bsp_unusual_fmt, tmp, bsp->length, num_bad);
      } else {
        row_text = (CharPtr) MemNew (sizeof(Char) * (StringLen (bsp_fmt) + StringLen (tmp) + 32));
        sprintf (row_text, bsp_fmt, tmp, bsp->length);
      }
      tmp = MemFree (tmp);
    }
  }
  else if (data_choice == OBJ_BIOSEQSET) 
  {
    tmp = GetBioseqSetLabel (vnp->data.ptrvalue);
    row_text = (CharPtr) MemNew (sizeof(Char) * (StringLen (tmp) + 2));
    sprintf (row_text, "%s\n", tmp);
    tmp = MemFree (tmp);
  }
  else if (data_choice == OBJ_SEQENTRY)
  {
    sep = (SeqEntryPtr) vnp->data.ptrvalue;
    if (sep != NULL && sep->data.ptrvalue != NULL) {
      tmp = NULL;
      if (IS_Bioseq(sep)) {
        tmp = GetBioseqLabel (sep->data.ptrvalue);
      } else if (IS_Bioseq_set (sep)) {
        tmp = GetBioseqSetLabel (sep->data.ptrvalue);
      }
      if (tmp != NULL) {
        row_text = (CharPtr) MemNew (sizeof(Char) * (StringLen (tmp) + 2));
        sprintf (row_text, "%s\n", tmp);
        tmp = MemFree (tmp);
      }
    }
  }
  else if (data_choice == OBJ_SEQDESC)
  {
    sdp = (SeqDescrPtr) vnp->data.ptrvalue;
    if (sdp != NULL)
    {
      bsp = NULL;
      if (sdp->extended != 0) {
        ovn = (ObjValNodePtr) sdp;
        if (ovn->idx.parenttype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) ovn->idx.parentptr;
        } else if (ovn->idx.parenttype == OBJ_BIOSEQSET && ovn->idx.parentptr != NULL) {
          bsp = GetRepresentativeBioseqFromBioseqSet (ovn->idx.parentptr);
        }
      }
      if (bsp == NULL) {
        if (sdp->choice == Seq_descr_title || sdp->choice == Seq_descr_comment) {
          row_text = (CharPtr) MemNew (sizeof (Char) * (StringLen ((CharPtr)(sdp->data.ptrvalue)) + 2));
          StringCpy (row_text, (CharPtr)(sdp->data.ptrvalue));
        } else {
          row_text = (CharPtr) MemNew (sizeof (Char) * 61);
          SeqDescLabel (sdp, row_text, 59, TRUE);
        }
      } else {
        label_len = 61;
        if (sdp->choice == Seq_descr_title || sdp->choice == Seq_descr_comment) {
          label_len = StringLen (sdp->data.ptrvalue) + 3;
        }
        row_text = (CharPtr) MemNew (sizeof (Char) * (label_len + 41));
        SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), row_text, PRINTID_REPORT, 39);
        row_text[39] = 0;
        StringCat (row_text, ":");
        if (sdp->choice == Seq_descr_title || sdp->choice == Seq_descr_comment) {
          StringCat (row_text, (CharPtr)(sdp->data.ptrvalue));
        } else {
          SeqDescLabel (sdp, row_text + StringLen (row_text), 59, TRUE);
        }
      }
      StringCat (row_text, "\n");
    }
  } else if (data_choice == OBJ_SEQSUB) {
    ssp = (SeqSubmitPtr) vnp->data.ptrvalue;
    if (ssp != NULL && ssp->datatype == 1 && (sep = ssp->data) != NULL && sep->data.ptrvalue != NULL) {
      tmp = NULL;
      if (IS_Bioseq(sep)) {
        tmp = GetBioseqLabel (sep->data.ptrvalue);
      } else if (IS_Bioseq_set (sep)) {
        tmp = GetBioseqSetLabel (sep->data.ptrvalue);
      }
      if (tmp != NULL) {
        row_text = (CharPtr) MemNew (sizeof(Char) * (StringLen (tmp) + 14));
        sprintf (row_text, "Cit-sub for %s\n", tmp);
        tmp = MemFree (tmp);
      }
    }
  }

  if (extra_fields != NULL) {
    for (field = extra_fields; field != NULL; field = field->next) {
      field_values = GetMultipleFieldValuesForObject (vnp->choice, vnp->data.ptrvalue, field, NULL, NULL);
      if (field_values != NULL) {
        label = SummarizeFieldType (field);
        for (val_vnp = field_values; val_vnp != NULL; val_vnp = val_vnp->next) {
          if (!StringHasNoText (label)) {
            tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (label) + 6));
            sprintf (tmp, "    %s:", label);
            ValNodeAddPointer (&field_strings, 0, tmp);
            field_len += StringLen (tmp) + 1;
          }
          if (label == NULL) {
            tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (val_vnp->data.ptrvalue) + 6));
            sprintf (tmp, "    %s\n", (CharPtr) val_vnp->data.ptrvalue);
          } else {
            tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (val_vnp->data.ptrvalue) + 2));
            sprintf (tmp, "%s\n", (CharPtr) val_vnp->data.ptrvalue);
          }
          ValNodeAddPointer (&field_strings, 0, tmp);
          field_len += StringLen (tmp);
        }
        label = MemFree (label);
        field_values = ValNodeFreeData (field_values);
      }
    }
    if (field_strings != NULL) {
      tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (row_text) + field_len + 3));
      StringCpy (tmp, row_text);
      /* replace trailing carriage return with space */
      for (field = field_strings; field != NULL; field = field->next) {
        StringCat (tmp, field->data.ptrvalue);
      }
      row_text = MemFree (row_text);
      row_text = tmp;
      field_strings = ValNodeFreeData (field_strings);
    }
  }

  if (!StringHasNoText (row_text) && !StringHasNoText (filename)) {
    tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (filename) + StringLen (row_text)));
    sprintf (tmp, fmt, filename, row_text);
    row_text = MemFree (row_text);
    row_text = tmp;
  }
  
  return row_text;
}

extern CharPtr GetDiscrepancyItemText (ValNodePtr vnp)
{
  return GetDiscrepancyItemTextEx (vnp, NULL);
}

extern CharPtr GetParentLabelForDiscrepancyItem (ValNodePtr vnp)
{
  CharPtr label = NULL;
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;
  ObjValNodePtr ovn;
  BioseqPtr  bsp;

  if (vnp == NULL || vnp->data.ptrvalue == NULL) return NULL;

  switch (vnp->choice)
  {
    case OBJ_SEQFEAT:
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      bsp = BioseqFindFromSeqLoc (sfp->location);
      label = GetBioseqLabel (bsp);
      break;
    case OBJ_SEQDESC:
      sdp = (SeqDescrPtr) vnp->data.ptrvalue;
      if (sdp != NULL)
      {
        if (sdp->extended != 0) {
          ovn = (ObjValNodePtr) sdp;
          if (ovn->idx.parenttype == OBJ_BIOSEQ) {
            label = GetBioseqLabel ((BioseqPtr) ovn->idx.parentptr);
          } else if (ovn->idx.parenttype == OBJ_BIOSEQSET) {
            label = GetBioseqSetLabel ((BioseqSetPtr) ovn->idx.parentptr);
          }
        }
      }
      break;
    case OBJ_BIOSEQ:
      label = GetBioseqLabel (vnp->data.ptrvalue);
      break;
    case OBJ_BIOSEQSET:
      label = GetBioseqSetLabel (vnp->data.ptrvalue);
      break;
  }
  return label;
}


static int StringCompareWithNumbers (CharPtr str1, CharPtr str2)
{
  int rval = 0;
  CharPtr cp1, cp2;
  int val1, val2;

  if (str1 == NULL && str2 == NULL) {
    rval = 0;
  } else if (str1 == NULL) {
    rval = -1;
  } else if (str2 == NULL) {
    rval = 1;
  } else {
    cp1 = str1;
    cp2 = str2;
    while (*cp1 != 0 && *cp2 != 0 && rval == 0) {
      if (isdigit (*cp1) && isdigit (*cp2)) {
        val1 = atoi (cp1);
        val2 = atoi (cp2);
        if (val1 < val2) {
          rval = -1;
        } else if (val1 > val2) {
          rval = 1;
        }
        while (isdigit (*cp1)) {
          cp1++;
        }
        while (isdigit (*cp2)) {
          cp2++;
        }
      } else if (*cp1 < *cp2) {
        rval = -1;
      } else if (*cp1 > *cp2) {
        rval = 1;
      } else {
        cp1++;
        cp2++;
      }
    }
    if (*cp1 == 0 && *cp2 != 0) {
      rval = -1;
    } else if (*cp1 != 0 && *cp2 == 0) {
      rval = 1;
    }
  }
  return rval;
}


extern int LIBCALLBACK SortVnpByDiscrepancyDescription (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  ClickableItemPtr cip1, cip2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);

    if (vnp1->data.ptrvalue != NULL && vnp2->data.ptrvalue != NULL) {
      cip1 = vnp1->data.ptrvalue;
      cip2 = vnp2->data.ptrvalue;
      rval = StringCompareWithNumbers (cip1->description, cip2->description);
    }
  }

  return rval;
}


extern int LIBCALLBACK SortVnpByDiscrepancyItemText (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  CharPtr str1, str2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);

    str1 = GetDiscrepancyItemText (vnp1);
    str2 = GetDiscrepancyItemText (vnp2);
    rval = StringCompareWithNumbers (str1, str2);
    str1 = MemFree (str1);
    str2 = MemFree (str2);
  }

  return rval;
}



extern void ValNodeReverse (ValNodePtr PNTR list)
{
  ValNodePtr vnp_next, vnp, vnp_start = NULL;

  if (list == NULL) {
    return;
  }

  vnp = *list;
  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = vnp_start;
    vnp_start = vnp;
    vnp = vnp_next;
  }
  *list = vnp_start;
}


static ValNodePtr ValNodePointerDup (ValNodePtr vnp)
{
  ValNodePtr vnp_new = NULL;
  
  if (vnp != NULL)
  {
    vnp_new = ValNodeNew (NULL);
    vnp_new->choice = vnp->choice;
    vnp_new->data.ptrvalue = vnp->data.ptrvalue;
    vnp_new->next = ValNodePointerDup (vnp->next);
  }
  return vnp_new;
}


typedef struct ftstrings {
  CharPtr header;
  CharPtr desc;
} FTStringsData, PNTR FTStringsPtr;


static FTStringsPtr FTStringsNew (CharPtr header, CharPtr desc)
{
  FTStringsPtr f;

  f = (FTStringsPtr) MemNew (sizeof (FTStringsData));
  f->header = header;
  f->desc = desc;
  return f;
}


static FTStringsPtr FTStringsFree (FTStringsPtr f)
{
  if (f != NULL) {
    f->header = MemFree (f->header);
    f->desc = MemFree (f->desc);
    f = MemFree (f);
  }
  return f;
}


extern ValNodePtr ReplaceDiscrepancyItemWithFeatureTableStrings (ValNodePtr feat_list)
{
  BioseqPtr       bsp, prot_bsp;
  CstType         custom_flags = 0;
  Asn2gbJobPtr    ajp;
  BaseBlockPtr    bbp;
  XtraBlock       extra;
  Int4            index;
  SeqFeatPtr      sfp, cds;
  ValNodePtr      vnp, list_copy = NULL, list_vnp;
  CharPtr         feature_table_header = NULL, feat_desc;
  FTStringsPtr    fts;
  
  if (feat_list == NULL) return NULL;
  
  list_copy = ValNodePointerDup (feat_list);
  for (vnp = list_copy; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != OBJ_SEQFEAT || vnp->data.ptrvalue == NULL) continue;

    sfp = (SeqFeatPtr) vnp->data.ptrvalue;

    if (sfp->idx.subtype == FEATDEF_PROT) {
      prot_bsp = BioseqFindFromSeqLoc (sfp->location);
      if (prot_bsp != NULL) {
        cds = SeqMgrGetCDSgivenProduct (prot_bsp, NULL);
        if (cds != NULL) {
          sfp = cds;
        }
      }      
    }      
    bsp = BioseqFindFromSeqLoc (sfp->location);
    feature_table_header = NULL;
    MemSet ((Pointer) &extra, 0, sizeof (XtraBlock));
    ajp = asn2gnbk_setup (bsp, NULL, NULL, FTABLE_FMT, DUMP_MODE, NORMAL_STYLE,
                          0, 0, custom_flags, &extra);
    if (ajp == NULL) {
      continue;
    }

    for (index = 0; index < ajp->numParagraphs; index++) 
    {
      bbp = ajp->paragraphArray [index];
      if (bbp->blocktype == FEATHEADER_BLOCK) {
        feature_table_header = asn2gnbk_format (ajp, (Int4) index);
      } else if (bbp->blocktype == FEATURE_BLOCK) {
        for (list_vnp = vnp; list_vnp != NULL; list_vnp = list_vnp->next)
        {
          if (list_vnp->choice != OBJ_SEQFEAT || list_vnp->data.ptrvalue == NULL) continue;
          sfp = (SeqFeatPtr) list_vnp->data.ptrvalue;
          if (sfp != NULL && sfp->idx.subtype == FEATDEF_PROT) {
            prot_bsp = BioseqFindFromSeqLoc (sfp->location);
            if (prot_bsp != NULL) {
              cds = SeqMgrGetCDSgivenProduct (prot_bsp, NULL);
              if (cds != NULL) {
                sfp = cds;
              }
            }
          }      
           
          if (sfp != NULL 
              && bbp->entityID == sfp->idx.entityID
              && bbp->itemtype == sfp->idx.itemtype
              && bbp->itemID == sfp->idx.itemID)
          {
            /* replace list feature with description, change choice */
            list_vnp->choice = 0;
            feat_desc = asn2gnbk_format (ajp, (Int4) index);
            list_vnp->data.ptrvalue = FTStringsNew (StringSave (feature_table_header), feat_desc);
          }
        }
      }
    }
    asn2gnbk_cleanup (ajp);
    feature_table_header = MemFree (feature_table_header);
  }

  /* now remove redundant headers */
  for (list_vnp = list_copy; list_vnp != NULL; list_vnp = list_vnp->next) {
    if (list_vnp->choice != 0) continue;
    fts = (FTStringsPtr) list_vnp->data.ptrvalue;
    if (feature_table_header == NULL
        || StringCmp (feature_table_header, fts->header) != 0) {
      feature_table_header = MemFree (feature_table_header);
      feature_table_header = fts->header;
      fts->header = NULL;
      list_vnp->data.ptrvalue = (CharPtr) MemNew (sizeof (Char) * (StringLen (feature_table_header) + StringLen (fts->desc) + 2));
      StringCpy (list_vnp->data.ptrvalue, feature_table_header);
      StringCat (list_vnp->data.ptrvalue, fts->desc);
    } else {
      list_vnp->data.ptrvalue = fts->desc;
      fts->desc = NULL;
    }
    fts = FTStringsFree (fts);
  }
  feature_table_header = MemFree (feature_table_header);
  return list_copy;
}

static void StandardWriteDiscrepancy (FILE *fp, ClickableItemPtr dip, Boolean use_feature_table_fmt, CharPtr descr_prefix, Boolean list_features_if_subcat)
{
  ValNodePtr vnp, list_copy = NULL;
  CharPtr    row_text;
  
  if (fp == NULL || dip == NULL)
  {
    return;
  }
  
  if (!StringHasNoText (descr_prefix)) {
    fprintf (fp, "%s:", descr_prefix);
  }
  fprintf (fp, "%s\n", dip->description);

  if (dip->subcategories == NULL || list_features_if_subcat) {
    vnp = dip->item_list;
    
    if (use_feature_table_fmt)
    {
      list_copy = ReplaceDiscrepancyItemWithFeatureTableStrings (vnp);
      vnp = list_copy;
    }

    while (vnp != NULL)
    {
      if (vnp->choice == 0)
      {
        row_text = StringSave (vnp->data.ptrvalue);
      }
      else
      {
        row_text = GetDiscrepancyItemText (vnp);
      }
      if (row_text != NULL)
      {
        fprintf (fp, "%s", row_text);
        row_text = MemFree (row_text);
      }
      vnp = vnp->next;
    }
    
    fprintf (fp, "\n");
  }
}


static Boolean SuppressItemListForFeatureTypeForOutputFiles (Uint4 test_type)
{
  if (test_type == DISC_FEATURE_COUNT 
    || test_type == DISC_MISSING_SRC_QUAL
    || test_type == DISC_DUP_SRC_QUAL
    || test_type == DISC_DUP_SRC_QUAL_DATA
    || test_type == DISC_SOURCE_QUALS_ASNDISC) {
    return TRUE;
  } else {
    return FALSE;
  }
}

extern void WriteDiscrepancyEx (FILE *fp, ClickableItemPtr dip, Boolean use_feature_table_fmt, Boolean cmdline, CharPtr descr_prefix, Boolean list_features_if_subcat)
{
  ValNodePtr vnp;

  if (fp == NULL || dip == NULL) {
    return;
  }

  if (cmdline && SuppressItemListForFeatureTypeForOutputFiles (dip->clickable_item_type)) {
    if (!StringHasNoText (descr_prefix)) {
      fprintf (fp, "%s:", descr_prefix);
    }
    fprintf (fp, "%s\n", dip->description);
    for (vnp = dip->subcategories; vnp != NULL; vnp = vnp->next) {
      dip = vnp->data.ptrvalue;
      if (dip != NULL) {
        if (!StringHasNoText (descr_prefix)) {
          fprintf (fp, "%s:", descr_prefix);
        }
        fprintf (fp, "%s\n", dip->description);
      }
    }
  } else {
    StandardWriteDiscrepancy (fp, dip, use_feature_table_fmt, descr_prefix, list_features_if_subcat);
  }
}


extern void WriteDiscrepancy (FILE *fp, ClickableItemPtr dip, Boolean use_feature_table_fmt)
{
  WriteDiscrepancyEx (fp, dip, use_feature_table_fmt, FALSE, NULL, TRUE);
}









/* DiscrepancyConfig functions */
extern DiscrepancyConfigPtr DiscrepancyConfigFree (DiscrepancyConfigPtr dcp)
{
  return MemFree (dcp);  
}

extern void DisableTRNATests (DiscrepancyConfigPtr dcp)
{
  if (dcp != NULL) {
    dcp->conf_list[DISC_COUNT_TRNA] = FALSE;
    dcp->conf_list[DISC_DUP_TRNA] = FALSE;
    dcp->conf_list[DISC_BADLEN_TRNA] = FALSE;
    dcp->conf_list[DISC_COUNT_RRNA] = FALSE;
    dcp->conf_list[DISC_DUP_RRNA] = FALSE;
    dcp->conf_list[DISC_TRANSL_NO_NOTE] = FALSE;
    dcp->conf_list[DISC_NOTE_NO_TRANSL] = FALSE;
    dcp->conf_list[DISC_TRANSL_TOO_LONG] = FALSE;
    dcp->conf_list[DISC_CDS_OVERLAP_TRNA] = FALSE;
    dcp->conf_list[DISC_COUNT_PROTEINS] = FALSE;
  }
}

extern DiscrepancyConfigPtr DiscrepancyConfigNew (void)
{
  DiscrepancyConfigPtr dcp;
  Int4                 i;
  
  dcp = (DiscrepancyConfigPtr) MemNew (sizeof (DiscrepancyConfigData));
  for (i = 0; i < MAX_DISC_TYPE; i++)
  {
    dcp->conf_list[i] = TRUE;
  }

  dcp->use_feature_table_format = FALSE;
  return dcp;
}


extern DiscrepancyConfigPtr DiscrepancyConfigCopy (DiscrepancyConfigPtr dcp)
{
  DiscrepancyConfigPtr cpy = NULL;

  if (dcp != NULL) {
    cpy = (DiscrepancyConfigPtr) MemNew (sizeof (DiscrepancyConfigData));
    MemCpy (cpy, dcp, sizeof (DiscrepancyConfigData));
  }
  return cpy;
}


extern DiscrepancyConfigPtr ReadDiscrepancyConfigEx (CharPtr report_config_name)
{
  DiscrepancyConfigPtr dcp;
  Int4                 i;
  Char                 str[20];
  
  dcp = DiscrepancyConfigNew();
  if (StringCmp (report_config_name, "DISCREPANCY_REPORT") == 0) {
    DisableTRNATests (dcp);
  }
  if (dcp != NULL)
  {
    for (i = 0; i < MAX_DISC_TYPE; i++)
    {
      if (GetAppParam ("SEQUINCUSTOM", report_config_name, discrepancy_info_list[i].setting_name, NULL, str, sizeof (str))) {
        if (StringICmp (str, "FALSE") == 0) {
          dcp->conf_list[i] = FALSE;
        } else if (StringICmp (str, "TRUE") == 0) {
          dcp->conf_list[i] = TRUE;
        }
      }
    }
    if (GetAppParam ("SEQUINCUSTOM", report_config_name, "USE_FEATURE_TABLE_FORMAT", NULL, str, sizeof (str))) {
      if (StringICmp (str, "TRUE") == 0) {
        dcp->use_feature_table_format = TRUE;
      }
    }
  }
  return dcp;
}

extern DiscrepancyConfigPtr ReadDiscrepancyConfig (void)
{
  return ReadDiscrepancyConfigEx ("DISCREPANCY_REPORT");
}

extern void SaveDiscrepancyConfigEx (DiscrepancyConfigPtr dcp, CharPtr report_name)
{
  Int4 i;
  
  if (dcp == NULL)
  {
    return;
  }

  if (report_name == NULL) {
    report_name = "DISCREPANCY_REPORT";
  }
  
  for (i = 0; i < MAX_DISC_TYPE; i++)
  {
    if (dcp->conf_list[i])
    {
      SetAppParam ("SEQUINCUSTOM", report_name, discrepancy_info_list[i].setting_name, "TRUE");
    }
    else
    {
      SetAppParam ("SEQUINCUSTOM", report_name, discrepancy_info_list[i].setting_name, "FALSE");
    }
  }
  if (dcp->use_feature_table_format)
  {
    SetAppParam ("SEQUINCUSTOM", report_name, "USE_FEATURE_TABLE_FORMAT", "TRUE");
  }
  else
  {
    SetAppParam ("SEQUINCUSTOM", report_name, "USE_FEATURE_TABLE_FORMAT", "FALSE");
  }
}


extern void SaveDiscrepancyConfig (DiscrepancyConfigPtr dcp)
{
  SaveDiscrepancyConfigEx (dcp, "DISCREPANCY_REPORT");
}


extern CharPtr SetDiscrepancyReportTestsFromString (CharPtr list, Boolean enable, DiscrepancyConfigPtr dcp)
{
  CharPtr         ptr, tmp, name_start, err_msg;
  DiscrepancyType test_type;
  CharPtr         err_fmt = "%s is an unrecognized test name";
  Int4            i;
  
  if (dcp == NULL) return StringSave ("Unable to configure");

  if (!StringDoesHaveText (list)) {
      return StringSave ("No tests specified!");
  }

  tmp = StringSave (list);
  name_start = tmp;
  if (StringICmp (name_start, "ALL") == 0) {
    for (i = 0; i < MAX_DISC_TYPE; i++) {
      dcp->conf_list[i] = enable;
    }
  } else {
    while (name_start != NULL && StringDoesHaveText (name_start)) {
      ptr = StringChr (name_start, ',');
      if (ptr != NULL) {
        *ptr = 0;
      }
      TrimSpacesAroundString (name_start);
      test_type = GetDiscrepancyTypeFromSettingName (name_start);
      if (test_type == MAX_DISC_TYPE) {
        err_msg = (CharPtr) MemNew (StringLen (err_fmt) + StringLen (name_start));
        sprintf (err_msg, err_fmt, name_start);
        tmp = MemFree (tmp);
        return err_msg;
      }
      dcp->conf_list[test_type] = enable;
      if (ptr == NULL) {
        name_start = NULL;
      } else {
        name_start = ptr + 1;
      }
    }
  }
  tmp = MemFree (tmp);
  return NULL;  
}


static Boolean OkToExpand (ClickableItemPtr cip, DiscReportOutputConfigPtr oc)
{

  if (cip == NULL || oc == NULL) {
    return FALSE;
  } else if (cip->clickable_item_type == DISC_FEATURE_COUNT) {
    return FALSE;
  } else if ((cip->item_list == NULL || oc->expand_report_categories[cip->clickable_item_type]) 
             && cip->subcategories != NULL) {
    return TRUE;
  } else {
    return FALSE;
  }
}


/* functions for writing discrepancy report to file */
static void WriteAsnDiscReportEx (ValNodePtr discrepancy_list, FILE *ofp, DiscReportOutputConfigPtr oc, Boolean use_flag, Boolean subcategory)
{
  ValNodePtr       vnp;
  ClickableItemPtr cip;
  CharPtr          setting_name, prefix;
  CharPtr          prefix_fmt = "DiscRep%s:%s:";

  if (ofp == NULL || oc == NULL) return;

  for (vnp = discrepancy_list; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL) {
      prefix = NULL;
      if (use_flag) {
        setting_name = GetDiscrepancyTestSettingName ((DiscrepancyType) cip->clickable_item_type);
        if (StringHasNoText (setting_name)) {
          if (subcategory) {
            prefix = StringSave ("DiscRep_SUB:");
          } else {
            prefix = StringSave ("DiscRep_ALL:");
          }
        } else {
          prefix = (CharPtr) MemNew (sizeof (Char) * (StringLen (prefix_fmt) + StringLen (setting_name) + 4));
          sprintf (prefix, prefix_fmt, subcategory ? "_SUB" : "_ALL", setting_name);
        }
      }

      if (oc->summary_report) {
        fprintf (ofp, "%s%s\n", prefix == NULL ? "" : prefix, cip->description);           
      } else {
        WriteDiscrepancyEx (ofp, cip, oc->use_feature_table_format, use_flag, prefix,
                            !oc->expand_report_categories[cip->clickable_item_type]);
      }
      prefix = MemFree (prefix);
      if (OkToExpand (cip, oc)) {
        if (use_flag && cip->clickable_item_type == DISC_INCONSISTENT_BIOSRC_DEFLINE) {
          WriteAsnDiscReport (cip->subcategories, ofp, oc, FALSE);
        } else {
          WriteAsnDiscReportEx (cip->subcategories, ofp, oc, use_flag, TRUE);
        }
      }
    }
  }

}

extern void WriteAsnDiscReport (ValNodePtr discrepancy_list, FILE *ofp, DiscReportOutputConfigPtr oc, Boolean use_flag)
{
  WriteAsnDiscReportEx (discrepancy_list, ofp, oc, use_flag, FALSE);
}


static int LIBCALLBACK SortVnpByDiscrepancyType (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  ClickableItemPtr c1, c2;
  CharPtr          cp1, cp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      c1 = (ClickableItemPtr) vnp1->data.ptrvalue;
      c2 = (ClickableItemPtr) vnp2->data.ptrvalue;
      if (c1 != NULL && c2 != NULL) {
        if (c1->clickable_item_type < c2->clickable_item_type) {
          return -1;
        } else if (c1->clickable_item_type > c2->clickable_item_type) {
          return 1;
        } else {
          if (c1->description == NULL && c2->description == NULL) {
            return 0;
          } else if (c1->description == NULL) {
            return -1;
          } else if (c2->description == NULL) {
            return 1;
          } else {
            cp1 = c1->description;
            while (isdigit (*cp1)) {
              cp1++;
            }
            cp2 = c2->description;
            while (isdigit (*cp2)) {
              cp2++;
            }
            return StringCmp (cp1, cp2);
          }
        }
      }
    }
  }
  return 0;
}


static ClickableItemPtr CombineDiscrepancyReports (ClickableItemPtr cip1, ClickableItemPtr cip2)
{
  CharPtr cp1, cp2, num_start1, num_start2, num_buf;
  Char    fixed_buf[15];
  Int4    common_start_len = 0;
  Int4    num_len1, num_len2, num_items1, num_items2;
  ClickableItemPtr combined = NULL;
  

  if (cip1 == NULL || cip2 == NULL || cip1->clickable_item_type != cip2->clickable_item_type
      || StringHasNoText (cip1->description) || StringHasNoText (cip2->description)) {
    return NULL;
  }

  if (cip1->clickable_item_type == DISC_QUALITY_SCORES) {
    /* special case for quality scores */
    combined = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    combined->clickable_item_type = cip1->clickable_item_type;
    if (StringCmp (cip1->description, cip2->description) == 0) {
      combined->description = StringSave (cip1->description);
    } else {
      combined->description= StringSave ("Quality scores are missing on some sequences.");
    }
    combined->item_list = cip1->item_list;
    cip1->item_list = NULL;
    combined->subcategories = cip1->subcategories;
    cip1->subcategories = NULL;
    ValNodeLink (&(combined->item_list), cip2->item_list);
    cip2->item_list = NULL;
    ValNodeLink (&(combined->subcategories), cip2->subcategories);
    cip2->subcategories = NULL;
  } else {
    /* all other tests */
    cp1 = cip1->description;
    cp2 = cip2->description;
      
    while (*cp1 == *cp2 && *cp1 != 0 && *cp2 != 0 && !isdigit (*cp1)) {
      cp1++;
      cp2++;
      common_start_len++;
    }
    if (*cp1 == 0 && *cp2 == 0) {
      /* entire description matches */
      combined = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      combined->clickable_item_type = cip1->clickable_item_type;
      combined->description = StringSave (cip1->description);
      combined->item_list = cip1->item_list;
      cip1->item_list = NULL;
      combined->subcategories = cip1->subcategories;
      cip1->subcategories = NULL;
      ValNodeLink (&(combined->item_list), cip2->item_list);
      cip2->item_list = NULL;
      ValNodeLink (&(combined->subcategories), cip2->subcategories);
      cip2->subcategories = NULL;
    } else if (isdigit (*cp1) && isdigit (*cp2) && (cp1 == cip1->description || isspace (*(cp1 - 1)))) {
      num_start1 = cp1;
      num_len1 = 0;
      while (isdigit (*cp1)) {
        cp1++;
        num_len1++;
      }
      num_start2 = cp2;
      num_len2 = 0;
      while (isdigit (*cp2)) {
        cp2++;
        num_len2++;
      }
      if ((*cp1 == 0 || isspace (*cp1)) && StringCmp (cp1, cp2) == 0) {
        /* matches on the other side of the number */
        /* build combined description */
        if (num_len1 < sizeof (fixed_buf)) {
          StringNCpy (fixed_buf, num_start1, num_len1);
          fixed_buf[num_len1] = 0;
          num_items1 = atoi(fixed_buf);
        } else {
          num_buf = (CharPtr) MemNew (sizeof (Char) * (num_len1 + 1));
          StringNCpy (num_buf, num_start1, num_len1);
          num_buf[num_len1] = 0;
          num_items1 = atoi (num_buf);
          num_buf = MemFree (num_buf);
        }
        if (num_len2 < sizeof (fixed_buf) - 1) {
          StringNCpy (fixed_buf, num_start2, num_len2);
          fixed_buf[num_len2] = 0;
          num_items2 = atoi(fixed_buf);
        } else {
          num_buf = (CharPtr) MemNew (sizeof (Char) * (num_len2 + 1));
          StringNCpy (num_buf, num_start2, num_len2);
          num_buf[num_len2] = 0;
          num_items2 = atoi (num_buf);
          num_buf = MemFree (num_buf);
        }
        
        combined = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));

        combined->description = (CharPtr) MemNew (sizeof (Char) * (common_start_len + sizeof (fixed_buf) + StringLen (cp1) + 1));
        StringNCpy (combined->description, cip1->description, common_start_len);
        sprintf (fixed_buf, "%d", num_items1 + num_items2);
        StringCat (combined->description, fixed_buf);
        StringCat (combined->description, cp1);

        combined->clickable_item_type = cip1->clickable_item_type;  
        combined->item_list = cip1->item_list;
        cip1->item_list = NULL;
        combined->subcategories = cip1->subcategories;
        cip1->subcategories = NULL;
        ValNodeLink (&(combined->item_list), cip2->item_list);
        cip2->item_list = NULL;
        ValNodeLink (&(combined->subcategories), cip2->subcategories);
        cip2->subcategories = NULL;
      } else {
        combined = NULL;
      }
    }
  }
  if (combined != NULL && combined->subcategories != NULL) {
    CollateDiscrepancyReports (&(combined->subcategories));
  }
  return combined;
}


extern void CollateDiscrepancyReports (ValNodePtr PNTR discrepancy_reports)
{
  ValNodePtr vnp, tmp;
  ClickableItemPtr combined;

  *discrepancy_reports = ValNodeSort (*discrepancy_reports, SortVnpByDiscrepancyType);

  vnp = *discrepancy_reports;
  while (vnp != NULL && vnp->next != NULL) {
    combined = CombineDiscrepancyReports (vnp->data.ptrvalue, vnp->next->data.ptrvalue);
    if (combined != NULL) {
      vnp->data.ptrvalue = ClickableItemFree (vnp->data.ptrvalue);
      vnp->next->data.ptrvalue = ClickableItemFree (vnp->next->data.ptrvalue);
      tmp = vnp->next;
      vnp->next = vnp->next->next;
      tmp->next = NULL;
      tmp = ValNodeFree (tmp);
      vnp->data.ptrvalue = combined;
    } else {
      vnp = vnp->next;
    }
  }
}


extern CharPtr ExpandDiscrepancyReportTestsFromString (CharPtr list, Boolean expand, DiscReportOutputConfigPtr dcp)
{
  CharPtr         ptr, tmp, name_start, err_msg;
  Int4            i;
  DiscrepancyType test_type;
  CharPtr         err_fmt = "%s is an unrecognized test name";
  
  if (dcp == NULL) return StringSave ("Unable to configure");

  if (!StringDoesHaveText (list)) {
    return NULL;
  } else if (StringICmp (list, "all") == 0) {
    for (i = 0; i < MAX_DISC_TYPE; i++) {
      dcp->expand_report_categories[i] = expand;
    }
  } else {
    tmp = StringSave (list);
    name_start = tmp;
    while (name_start != NULL && StringDoesHaveText (name_start)) {
      ptr = StringChr (name_start, ',');
      if (ptr != NULL) {
        *ptr = 0;
      }
      TrimSpacesAroundString (name_start);
      test_type = GetDiscrepancyTypeFromSettingName (name_start);
      if (test_type == MAX_DISC_TYPE) {
        err_msg = (CharPtr) MemNew (StringLen (err_fmt) + StringLen (name_start));
        sprintf (err_msg, err_fmt, name_start);
        tmp = MemFree (tmp);
        return err_msg;
      }
      dcp->expand_report_categories[test_type] = expand;
      if (ptr == NULL) {
        name_start = NULL;
      } else {
        name_start = ptr + 1;
      }
    }
    tmp = MemFree (tmp);
  }
  return NULL;  
}


NLM_EXTERN DiscReportOutputConfigPtr DiscReportOutputConfigNew ()
{
  DiscReportOutputConfigPtr c;

  c = (DiscReportOutputConfigPtr) MemNew (sizeof (DiscReportOutputConfigData));
  MemSet (c, 0, sizeof (DiscReportOutputConfigData));

  return c;
}


NLM_EXTERN DiscReportOutputConfigPtr DiscReportOutputConfigFree (DiscReportOutputConfigPtr c)
{
  if (c != NULL) {
    c = MemFree (c);
  }
  return c;
}



/* The following section is for creating discrepancy reports for a large number of seq-entries,
 * which will not be available after each seq-entry has been added to the report.  Therefore
 * all item lists must be represented as strings.
 */

typedef struct globalsrcval {
  CharPtr src_id_txt;
  CharPtr val;
  ValNodePtr qual;
} GlobalSrcValData, PNTR GlobalSrcValPtr;

static GlobalSrcValPtr GlobalSrcValNew ()
{
  GlobalSrcValPtr g;

  g = (GlobalSrcValPtr) MemNew (sizeof (GlobalSrcValData));
  g->src_id_txt = NULL;
  g->val = NULL;
  g->qual = NULL;
  return g;
}


static GlobalSrcValPtr GlobalSrcValFree (GlobalSrcValPtr g)
{
  if (g != NULL) {
    g->src_id_txt = MemFree (g->src_id_txt);
    g->val = MemFree (g->val);
    g->qual = FieldTypeFree (g->qual);
    g = MemFree (g);
  }
  return g;
}


static ValNodePtr GlobalSrcValListFree (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = GlobalSrcValFree (list->data.ptrvalue);
    list = ValNodeFree (list);
    list = list_next;
  }
  return list;
}


static ValNodePtr GlobalSrcValListFromObject (ValNodePtr obj, ValNodePtr quals, CharPtr filename)
{
  GlobalSrcValPtr g;
  ValNodePtr vnp, list = NULL;
  CharPtr    str;

  if (obj == NULL || quals == NULL) {
    return NULL;
  }

  for (vnp = quals; vnp != NULL; vnp = vnp->next) {
    str = GetFieldValueForObject (obj->choice, obj->data.ptrvalue, vnp, NULL);
    if (StringHasNoText (str)) {
      str = MemFree (str);
    } else {
      g = GlobalSrcValNew ();
      g->src_id_txt = GetDiscrepancyItemTextEx (obj, filename);
      g->val = str;
      g->qual = (AsnIoMemCopy) (vnp, (AsnReadFunc) FieldTypeAsnRead, (AsnWriteFunc) FieldTypeAsnWrite);
      ValNodeAddPointer (&list, 0, g);
    }
  }
  return list;
}


static ValNodePtr SrcListFromGlobalSrcValList (ValNodePtr start, ValNodePtr stop)
{
  ValNodePtr src_list = NULL, vnp;
  GlobalSrcValPtr g;

  for (vnp = start; vnp != NULL; vnp = vnp->next) {
    g = (GlobalSrcValPtr) vnp->data.ptrvalue;
    ValNodeAddPointer (&src_list, 0, StringSave (g->src_id_txt));
    if (vnp == stop) {
      break;
    }
  }
  return src_list;
}


static int CompareGlobalSrcVal (GlobalSrcValPtr dq1, GlobalSrcValPtr dq2)
{
  int         rval = 0;

  if (dq1 != NULL && dq2 != NULL) {
    rval = CompareFieldTypes (dq1->qual, dq2->qual);
    if (rval == 0) {
      rval = StringCmp (dq1->val, dq2->val);
    }
    if (rval == 0) {
      rval = StringCmp (dq1->src_id_txt, dq2->src_id_txt);
    }
  }
  return rval;
}


static int LIBCALLBACK SortVnpByGlobalSrcVal (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  int         rval = 0;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);

    if (vnp1->data.ptrvalue != NULL && vnp2->data.ptrvalue != NULL) {
      rval = CompareGlobalSrcVal (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    }
  }

  return rval;
}


static ClickableItemPtr AnalyzeGlobalSrcVals (ValNodePtr src_list, ValNodePtr start, ValNodePtr stop, ClickableItemPtr cip_multi)
{
  ValNodePtr vnp_s, vnp, missing = NULL, present_src;
  ValNodePtr repeated = NULL, unique = NULL, dup_list = NULL;
  ClickableItemPtr missing_cip = NULL, cip, cip_dup;
  CharPtr          qual, fmt, missing_fmt = "%%d sources are missing %s", dup_fmt = "%%d sources have '%s' for %s";
  CharPtr          some_missing_some_dup = "%s (some missing, some duplicate%s)";
  CharPtr          some_missing = "%s (some missing, all unique%s)";
  CharPtr          some_dup = "%s (all present, some duplicate%s)";
  CharPtr          good = "%s (all present, all unique%s)";
  CharPtr          some_missing_all_same = "%s (some missing, all same%s)";
  CharPtr          all_present_all_same = "%s (all present, all same%s)";
  CharPtr          unique_fmt = "%%d sources have unique values for %s";
  CharPtr          some_multi = ", some multi";
  GlobalSrcValPtr  g1, g2;

  if (src_list == NULL || start == NULL || stop == NULL) {
    return NULL;
  }

  g1 = start->data.ptrvalue;
  qual = SummarizeFieldType (g1->qual);

  /* first, find missing quals */
  present_src = SrcListFromGlobalSrcValList (start, stop);
  present_src = ValNodeSort (present_src, SortVnpByString);

  vnp_s = src_list;
  vnp = present_src;
  while (vnp_s != NULL) {
    if (vnp == NULL) {
      ValNodeAddPointer (&missing, 0, StringSave (vnp_s->data.ptrvalue));
    } else if (StringCmp (vnp_s->data.ptrvalue, vnp->data.ptrvalue) != 0) {
      ValNodeAddPointer (&missing, 0, StringSave (vnp_s->data.ptrvalue));
    } else {
      vnp = vnp->next;
    }
    vnp_s = vnp_s->next;
  }
  present_src = ValNodeFreeData (present_src);

  if (missing != NULL) {
    fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (missing_fmt) + StringLen (qual)));
    sprintf (fmt, missing_fmt, qual);
    missing_cip = NewClickableItem (DISC_SOURCE_QUALS_ASNDISC, fmt, missing);
    fmt = MemFree (fmt);
  }

  /* now look for duplicates and unique values */
  g1 = start->data.ptrvalue;
  ValNodeAddPointer (&repeated, 0, g1->src_id_txt);
  if (start != stop) {
    for (vnp = start->next; vnp != NULL; vnp = vnp->next) {
      g2 = vnp->data.ptrvalue;
      if (StringCmp (g1->val, g2->val) != 0) {
        if (repeated->next == NULL) {
          ValNodeLink (&unique, repeated);
        } else {
          fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (dup_fmt) + StringLen (qual) + StringLen (g1->val)));
          sprintf (fmt, dup_fmt, g1->val, qual);
          ValNodeAddPointer (&dup_list, 0, NewClickableItem (DISC_SOURCE_QUALS_ASNDISC, fmt, repeated));
          fmt = MemFree (fmt);
        }
        repeated = NULL;
      }
      ValNodeAddPointer (&repeated, 0, g2->src_id_txt);
      g1 = g2;
      if (vnp == stop) {
        break;
      }
    }
  }

  if (repeated != NULL) {
    if (repeated->next == NULL) {
      ValNodeLink (&unique, repeated);
    } else {
      fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (dup_fmt) + StringLen (qual) + StringLen (g1->val)));
      sprintf (fmt, dup_fmt, g1->val, qual);
      ValNodeAddPointer (&dup_list, 0, NewClickableItem (DISC_SOURCE_QUALS_ASNDISC, fmt, repeated));
      fmt = MemFree (fmt);
    }
    repeated = NULL;
  }

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  cip->clickable_item_type = DISC_SOURCE_QUALS_ASNDISC;
  cip->item_list = NULL;
  cip->callback_func = NULL;
  cip->datafree_func = NULL;
  cip->callback_data = NULL;
  cip->chosen = 0;
  cip->expanded = FALSE;
  cip->level = 0;
  cip->subcategories = NULL;

  if (dup_list == NULL && missing == NULL) {
    fmt = good;
    cip->item_list = ValNodeDupStringList (src_list);
  } else if (dup_list != NULL && missing != NULL) {
    if (dup_list->next == NULL
        && (cip_dup = dup_list->data.ptrvalue) != NULL
        && ValNodeLen (cip_dup->item_list) == ValNodeLen (src_list) - ValNodeLen (missing)) {
      fmt = some_missing_all_same;
    } else {
      fmt = some_missing_some_dup;
    }
    ValNodeAddPointer (&(cip->subcategories), 0, missing_cip);
    ValNodeLink (&(cip->subcategories), dup_list);
  } else if (dup_list != NULL) {
    if (dup_list->next == NULL
        && (cip_dup = dup_list->data.ptrvalue) != NULL
        && ValNodeLen (cip_dup->item_list) == ValNodeLen (src_list)) {
      fmt = all_present_all_same;
    } else {      
      fmt = some_dup;
    }
    ValNodeLink (&(cip->subcategories), dup_list);
  } else if (missing != NULL) {
    fmt = some_missing;
    ValNodeAddPointer (&(cip->subcategories), 0, missing_cip);
  }

  if (cip_multi) {
    ValNodeAddPointer (&(cip->subcategories), 0, cip_multi);
  }

  if (fmt != NULL) {
    cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (qual) + (cip_multi == NULL ? 0 : StringLen (some_multi))));
    sprintf (cip->description, fmt, qual, cip_multi == NULL ? "" : some_multi);
  }

  if (unique != NULL) {
    fmt = (CharPtr) MemNew (sizeof (Char) * (StringLen (unique_fmt) + StringLen (qual)));
    sprintf (fmt, unique_fmt, qual);
    ValNodeAddPointer (&(cip->subcategories), 0, NewClickableItem (DISC_SOURCE_QUALS_ASNDISC, fmt, unique));
    fmt = MemFree (fmt);
  }

  qual = MemFree (qual);

  return cip;
}


static ValNodePtr ExtractMultiForQualFromList (ValNodePtr PNTR list, CharPtr qual)
{
  ValNodePtr vnp, vnp_prev = NULL;
  ClickableItemPtr cip;
  Char             found_qual[100];
  Int4             n;

  if (list == NULL || *list == NULL || StringHasNoText (qual)) {
    return NULL;
  }

  for (vnp = *list; vnp != NULL; vnp = vnp->next) {
    cip = vnp->data.ptrvalue;
    if (cip != NULL 
        && sscanf (cip->description, 
                     "%d sources have multiple %s qualifiers", &n, found_qual) == 2
                     && StringCmp (qual, found_qual) == 0) {
      if (vnp_prev == NULL) {
        *list = vnp->next;
      } else {
        vnp_prev->next = vnp->next;
      }
      vnp->next = NULL;
      return vnp;
    } else {
      vnp_prev = vnp;
    }
  }
  return NULL;
}
  

/* NOTE - I don't think we're actually going to need the qual list, but we will need the src_list */
static ValNodePtr 
GetMissingAndInconsistentDiscrepanciesFromGlobalSrcValList (ValNodePtr PNTR val_list, ValNodePtr PNTR src_list, ValNodePtr PNTR multi_list)
{
  ValNodePtr disc_list = NULL, start, last;
  ValNodePtr vnp, vnp_multi;
  GlobalSrcValPtr g1, g2;
  ClickableItemPtr cip_multi;
  CharPtr          qual;

  if (val_list == NULL || *val_list == NULL || src_list == NULL || *src_list == NULL) {
    return NULL;
  }
  
  *val_list = ValNodeSort (*val_list, SortVnpByGlobalSrcVal);
  *src_list = ValNodeSort (*src_list, SortVnpByString);

  g1 = (*val_list)->data.ptrvalue;
  start = *val_list;
  last = *val_list;
  for (vnp = *val_list; vnp != NULL; vnp = vnp->next) {
    g2 = vnp->data.ptrvalue;
    if (CompareFieldTypes (g1->qual, g2->qual) == 0) {
      last = vnp;
    } else {
      /* analyze from start to last */
      qual = SummarizeFieldType (g1->qual);
      vnp_multi = ExtractMultiForQualFromList (multi_list, qual);
      qual = MemFree (qual);
      if (vnp_multi == NULL) {
        cip_multi = NULL;
      } else {
        cip_multi = vnp_multi->data.ptrvalue;
      }
      vnp_multi = ValNodeFree (vnp_multi);
      ValNodeAddPointer (&disc_list, 0, AnalyzeGlobalSrcVals (*src_list, start, last, cip_multi));
      start = vnp;
      last = vnp;
      g1 = vnp->data.ptrvalue;
    }
  }

  /* analyze from start to last for last field*/
  qual = SummarizeFieldType (g1->qual);
  vnp_multi = ExtractMultiForQualFromList (multi_list, qual);
  qual = MemFree (qual);
  if (vnp_multi == NULL) {
    cip_multi = NULL;
  } else {
    cip_multi = vnp_multi->data.ptrvalue;
  }
  vnp_multi = ValNodeFree (vnp_multi);
  ValNodeAddPointer (&disc_list, 0, AnalyzeGlobalSrcVals (*src_list, start, last, cip_multi));
  
  return disc_list;
}


typedef struct globaldiscrepancylists {
  ValNodePtr locus_tag_list;
  ValNodePtr missing_locus_tag;
} GlobalDiscrepancyListsData, PNTR GlobalDiscrepancyListPtr;

static void CollectGlobalDiscrepancyData (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  GeneRefPtr         grp;
  GlobalDiscrepancyListPtr tbl;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_GENE) return;
  tbl = (GlobalDiscrepancyListPtr) userdata;
  if (tbl == NULL) return;

  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp != NULL) {
    if (grp->pseudo) return;
    if (StringDoesHaveText (grp->locus_tag)) {
      ValNodeAddPointer (&(tbl->locus_tag_list), 0, 
                          GlobalDiscrepancyNew (grp->locus_tag, OBJ_SEQFEAT, sfp));
    } else {
      ValNodeAddPointer (&(tbl->missing_locus_tag), 0,
                          GlobalDiscrepancyNew (NULL, OBJ_SEQFEAT, sfp));
    }
  }
}


static void SaveStringsForDiscrepancyItemList (ValNodePtr list, Boolean use_feature_fmt, CharPtr filename);

static void SaveStringsForDiscrepancyItems (ClickableItemPtr cip, Boolean use_feature_fmt, CharPtr filename)
{
  ValNodePtr vnp, list_copy;
  CharPtr    str = NULL;

  if (cip == NULL) return;
  if (cip->clickable_item_type == DISC_GENE_CDS_mRNA_LOCATION_CONFLICT)
  {
    str = str;
  }
  if (use_feature_fmt) {
    list_copy = ReplaceDiscrepancyItemWithFeatureTableStrings (cip->item_list);
    cip->item_list = ValNodeFree (cip->item_list);
    cip->item_list = list_copy;
  } else {
    for (vnp = cip->item_list; vnp != NULL; vnp = vnp->next) {
      str = GetDiscrepancyItemTextEx (vnp, filename);
      vnp->choice = 0;      
      vnp->data.ptrvalue = str;
    }
  }
  SaveStringsForDiscrepancyItemList (cip->subcategories, use_feature_fmt, filename);
}


static void SaveStringsForDiscrepancyItemList (ValNodePtr list, Boolean use_feature_fmt, CharPtr filename)
{
  while (list != NULL) {
    SaveStringsForDiscrepancyItems (list->data.ptrvalue, use_feature_fmt, filename);
    list = list->next;
  }
}


NLM_EXTERN GlobalDiscrepReportPtr GlobalDiscrepReportNew ()
{
  GlobalDiscrepReportPtr g;

  g = (GlobalDiscrepReportPtr) MemNew (sizeof (GlobalDiscrepReportData));
  MemSet (g, 0, sizeof (GlobalDiscrepReportData));
  g->output_config = DiscReportOutputConfigNew ();
  return g;
}


NLM_EXTERN GlobalDiscrepReportPtr GlobalDiscrepReportFree (GlobalDiscrepReportPtr g)
{
  if (g != NULL) {
    g->locus_tag_list = FreeGlobalDiscrepancyList (g->locus_tag_list);
    g->missing_locus_tag = FreeGlobalDiscrepancyList (g->missing_locus_tag);
    g->cds_product_list = FreeGlobalDiscrepancyList (g->cds_product_list);
    g->missing_cds_product = FreeGlobalDiscrepancyList (g->missing_cds_product);
    g->mrna_product_list = FreeGlobalDiscrepancyList (g->mrna_product_list);
    g->missing_mrna_product = FreeGlobalDiscrepancyList (g->missing_mrna_product);
    g->missing_gnl_list = FreeGlobalDiscrepancyList (g->missing_gnl_list);
    g->gnl_list = FreeGlobalDiscrepancyList (g->gnl_list);

    g->global_srcs = ValNodeFreeData (g->global_srcs);
    g->global_src_qual_vals = GlobalSrcValListFree (g->global_src_qual_vals);
    g->feature_count_list = FeatureCountListFree (g->feature_count_list); 

    g->src_qual_repeated_list = FreeClickableList (g->src_qual_repeated_list);
    g->src_qual_multi_list = FreeClickableList (g->src_qual_multi_list);
    g->discrepancy_list = FreeClickableList (g->discrepancy_list);
    g->output_config = DiscReportOutputConfigFree (g->output_config);
    g->test_config = DiscrepancyConfigFree (g->test_config);
    g = MemFree (g);
  }
  return g;
}


static void GetLocalSourceQualReportItems (ValNodePtr src_list, ValNodePtr qual_list, CharPtr filename, ValNodePtr PNTR repeated_list, ValNodePtr PNTR multi_list)
{
  ValNodePtr combo_list = NULL, disc_list = NULL, vnp;
  ValNodePtr vnp_q, vnp_s;
  DuplicateQualPtr dq;
  ClickableItemPtr cip_multi;

  /* get all values for all organisms */
  for (vnp_q = qual_list; vnp_q != NULL; vnp_q = vnp_q->next) {
    for (vnp_s = src_list; vnp_s != NULL; vnp_s = vnp_s->next) {
      dq = DuplicateQualNew (vnp_s->choice, vnp_s->data.ptrvalue, vnp_q);
      if (StringHasNoText (dq->val)) {
        dq = DuplicateQualFree (dq);
      } else {
        ValNodeAddPointer (&combo_list, 0, dq);
      }
    }
  }  
  /* now look for repeated field values in individual organisms */
  FindRepeatedFieldValues (&disc_list, &combo_list, DISC_SOURCE_QUALS_ASNDISC);
  combo_list = DuplicateQualListFree (combo_list);
  SaveStringsForDiscrepancyItemList (disc_list, FALSE, filename);
  ValNodeLink (repeated_list, disc_list);
  disc_list = NULL;

  /* also look for multiple quals */
  for (vnp = qual_list; vnp != NULL; vnp = vnp->next) {
    cip_multi = FindMultipleSourceQuals (vnp, src_list);
    if (cip_multi != NULL) {
      ValNodeAddPointer (&disc_list, 0, cip_multi);
    }
  }
  SaveStringsForDiscrepancyItemList (disc_list, FALSE, filename);
  ValNodeLink (multi_list, disc_list);
}


static void AddSourceQualReportInfoToGlobalDiscrepReport (SeqEntryPtr sep, GlobalDiscrepReportPtr g, CharPtr filename)
{
  ValNodePtr src_list, qual_list, vnp;

  src_list = GetObjectListForFieldType (FieldType_source_qual, sep);
  qual_list = GetSourceQualSampleFieldList (sep);
  AdjustSourceQualSampleFieldListForOnCallerTest (&qual_list, src_list);

  /* add local items */
  GetLocalSourceQualReportItems (src_list, qual_list, filename, &(g->src_qual_repeated_list), &(g->src_qual_multi_list));

  /* add to global src qual value list */
  for (vnp = src_list; vnp != NULL; vnp = vnp->next) {
    ValNodeLink (&(g->global_src_qual_vals), GlobalSrcValListFromObject (vnp, qual_list, filename));
    ValNodeAddPointer (&(g->global_srcs), 0, GetDiscrepancyItemTextEx (vnp, filename));
  }

  src_list = ValNodeFree (src_list);
}


NLM_EXTERN void AddSeqEntryToGlobalDiscrepReport (SeqEntryPtr sep, GlobalDiscrepReportPtr g, CharPtr filename)
{
  ClickableItemPtr adjacent_cip = NULL;
  ValNode          sep_list;
  ValNodePtr       local_discrepancy_list = NULL, local_counts = NULL;
  Uint2            entityID;
  DiscrepancyConfigPtr dcp;
  GlobalDiscrepancyListsData lists;
  GenProdSetDiscrepancyListsData gps_lists;
  ProtIdListsData                prot_lists;

  if (g == NULL || sep == NULL) return;

  entityID = SeqMgrGetEntityIDForSeqEntry (sep);
  if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
    SeqMgrIndexFeatures (entityID, NULL);
  }

  /* todo - some global tests might not be performed if they have been disabled */

  MemSet (&lists, 0, sizeof (GlobalDiscrepancyListsData));
  VisitGenProdSetFeatures (sep, &lists, CollectGlobalDiscrepancyData);
  MemSet (&gps_lists, 0, sizeof (GenProdSetDiscrepancyListsData));
  CheckGenProdSetsInSeqEntry (sep, &gps_lists);
  MemSet (&prot_lists, 0, sizeof (ProtIdListsData));
  VisitBioseqsInSep (sep, &prot_lists, FindProteinIDCallback);

  if (lists.locus_tag_list != NULL) {
    /* collect adjacent genes */
    lists.locus_tag_list = ValNodeSort (lists.locus_tag_list, SortVnpByGlobalDiscrepancyString);
    adjacent_cip = FindAdjacentDuplicateLocusTagGenes (lists.locus_tag_list);
    if (adjacent_cip != NULL) {
      SaveStringsForDiscrepancyItems (adjacent_cip, g->output_config->use_feature_table_format, filename);
      ValNodeAddPointer (&(g->adjacent_locus_tag_disc_list), 0, adjacent_cip);
    }
  }

  /* convert lists to strings and add to global lists */
  ConvertGlobalDiscrepancyListToText (lists.locus_tag_list, g->output_config->use_feature_table_format, filename);
  ValNodeLink (&(g->locus_tag_list), lists.locus_tag_list);
  ConvertGlobalDiscrepancyListToText (lists.missing_locus_tag, g->output_config->use_feature_table_format, filename);
  ValNodeLink (&(g->missing_locus_tag), lists.missing_locus_tag);
  ConvertGlobalDiscrepancyListToText (gps_lists.cds_product_list, g->output_config->use_feature_table_format, filename);
  ValNodeLink (&(g->cds_product_list), gps_lists.cds_product_list);
  ConvertGlobalDiscrepancyListToText (gps_lists.missing_protein_id, g->output_config->use_feature_table_format, filename);
  ValNodeLink (&(g->missing_cds_product), gps_lists.missing_protein_id);
  ConvertGlobalDiscrepancyListToText (gps_lists.mrna_product_list, g->output_config->use_feature_table_format, filename);
  ValNodeLink (&(g->mrna_product_list), gps_lists.mrna_product_list);
  ConvertGlobalDiscrepancyListToText (gps_lists.missing_mrna_product, g->output_config->use_feature_table_format, filename);
  ValNodeLink (&(g->missing_mrna_product), gps_lists.missing_mrna_product);
  ConvertGlobalDiscrepancyListToText (prot_lists.gnl_list, g->output_config->use_feature_table_format, filename);
  ValNodeLink (&g->gnl_list, prot_lists.gnl_list);
  ConvertGlobalDiscrepancyListToText (prot_lists.missing_gnl_list, g->output_config->use_feature_table_format, filename);
  ValNodeLink (&g->missing_gnl_list, prot_lists.missing_gnl_list);

  if (g->test_config->conf_list[DISC_SOURCE_QUALS_ASNDISC]) {
    AddSourceQualReportInfoToGlobalDiscrepReport (sep, g, filename);
  }

  if (g->test_config->conf_list[DISC_FEATURE_COUNT]) {
    VisitBioseqsInSep (sep, &local_counts, CountFeaturesOnSequenceCallback);
    SaveFeatureCountSequenceIds (local_counts, filename);
    ValNodeLink (&(g->feature_count_list), local_counts);
    local_counts = NULL;
  }

  dcp = DiscrepancyConfigCopy (g->test_config);

  /* disable tests that are global */
  dcp->conf_list[DISC_GENE_MISSING_LOCUS_TAG] = FALSE;
  dcp->conf_list[DISC_GENE_DUPLICATE_LOCUS_TAG] = FALSE;
  dcp->conf_list[DISC_GENE_LOCUS_TAG_BAD_FORMAT] = FALSE;
  dcp->conf_list[DISC_GENE_LOCUS_TAG_INCONSISTENT_PREFIX] = FALSE;
  dcp->conf_list[DISC_MISSING_GENPRODSET_PROTEIN] = FALSE;
  dcp->conf_list[DISC_DUP_GENPRODSET_PROTEIN] = FALSE;
  dcp->conf_list[DISC_MISSING_GENPRODSET_TRANSCRIPT_ID] = FALSE;
  dcp->conf_list[DISC_DUP_GENPRODSET_TRANSCRIPT_ID] = FALSE;
  dcp->conf_list[DISC_MISSING_PROTEIN_ID] = FALSE;
  dcp->conf_list[DISC_INCONSISTENT_PROTEIN_ID_PREFIX] = FALSE;
  dcp->conf_list[DISC_SOURCE_QUALS_ASNDISC] = FALSE;
  dcp->conf_list[DISC_FEATURE_COUNT] = FALSE;

  sep_list.data.ptrvalue = sep;
  sep_list.next = NULL;
  local_discrepancy_list = CollectDiscrepancies (dcp, &sep_list, g->taxlookup);

  dcp = DiscrepancyConfigFree (dcp);

  SaveStringsForDiscrepancyItemList (local_discrepancy_list, g->output_config->use_feature_table_format, filename);
  ValNodeLink (&(g->discrepancy_list), local_discrepancy_list);
}


static void PrintDiscrepancyReportSubcategories (ValNodePtr discrepancy_list, FILE *fp, Int4 indent)
{
  ValNodePtr vnp;
  ClickableItemPtr cip_sub;
  Int4 i;

  for (vnp = discrepancy_list; vnp != NULL; vnp = vnp->next) {
    cip_sub = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip_sub != NULL) {
      for (i = 0; i < indent; i++) {
        fprintf (fp, "\t");
      }
      fprintf (fp, "%s\n", cip_sub->description);
      PrintDiscrepancyReportSubcategories (cip_sub->subcategories, fp, indent + 1);
    }
  }
}


static void WriteDiscrepancyReportSummary (ValNodePtr discrepancy_list, FILE *fp)
{
  ClickableItemPtr cip;
  CharPtr          setting_name;

  while (discrepancy_list != NULL) {
    cip = discrepancy_list->data.ptrvalue;
    if (cip != NULL) {
      setting_name = GetDiscrepancyTestSettingName ((DiscrepancyType) cip->clickable_item_type);
      fprintf (fp, "%s:%s\n", setting_name, cip->description);
      if (cip->clickable_item_type == DISC_SUSPECT_PRODUCT_NAME) {
        PrintDiscrepancyReportSubcategories (cip->subcategories, fp, 1);
      }
    }
    discrepancy_list = discrepancy_list->next;
  }
}


static ClickableItemPtr GetReportForFeatdefList (ValNodePtr PNTR list, Int4 featdef)
{
  ClickableItemPtr cip = NULL;

  if (list == NULL || *list == NULL) {
    return NULL;
  }

  cip = AddFeatureTypeSummary (featdef, NULL, GetNumFeaturesInList (*list));
  return cip;
}


static ValNodePtr CreateGlobalFeatureCountReports (ValNodePtr PNTR feature_count_list)
{
  ValNodePtr feat_list, vnp, tmp_list, orig_list = NULL, disc_list = NULL;
  Int4       featdef;
  ClickableItemPtr cip;

  if (feature_count_list == NULL || *feature_count_list == NULL) {
    return NULL;
  }

  InsertMissingFeatureCountsWithSeqIdTxt (feature_count_list);
  feat_list = GetFeatureTypesFromFeatureCounts (*feature_count_list);
  for (vnp = feat_list; vnp != NULL; vnp = vnp->next) {
    featdef = vnp->data.intvalue;
    tmp_list = ValNodeExtractListByFunction (feature_count_list, FeatureCountHasFeatdef, &featdef);
    cip = GetReportForFeatdefList (&tmp_list, featdef);
    if (cip != NULL) {
      ValNodeAddPointer (&disc_list, 0, cip);
    }
    ValNodeLink (&orig_list, tmp_list);
    tmp_list = NULL;
  }

  feat_list = ValNodeFree (feat_list);
  *feature_count_list = orig_list;
  return disc_list;
}


NLM_EXTERN void WriteGlobalDiscrepancyReport (GlobalDiscrepReportPtr g, FILE *fp)
{
  ValNodePtr  local_list = NULL;
  ClickableItemPtr cip;

  if (g == NULL || fp == NULL) return;

  g->locus_tag_list = ValNodeSort (g->locus_tag_list, SortVnpByGlobalDiscrepancyString);
  g->missing_locus_tag = ValNodeSort (g->missing_locus_tag, SortVnpByGlobalDiscrepancyString);
  g->cds_product_list = ValNodeSort (g->cds_product_list, SortVnpByGlobalDiscrepancyString);
  g->missing_cds_product = ValNodeSort (g->missing_cds_product, SortVnpByGlobalDiscrepancyString);
  g->mrna_product_list = ValNodeSort (g->mrna_product_list, SortVnpByGlobalDiscrepancyString);
  g->missing_mrna_product = ValNodeSort (g->missing_mrna_product, SortVnpByGlobalDiscrepancyString);

  if (g->locus_tag_list != NULL) {    
    if (g->missing_locus_tag != NULL) {
      cip = ReportMissingFields (g->missing_locus_tag, discReportMissingLocusTags, DISC_GENE_MISSING_LOCUS_TAG);
      if (cip != NULL) {
        ValNodeAddPointer (&local_list, 0, cip);
      }
    }
    CollateDiscrepancyReports (&(g->adjacent_locus_tag_disc_list));
    cip = ReportNonUniqueGlobalDiscrepancy (g->locus_tag_list, 
                                            discReportDuplicateLocusTagFmt, 
                                            discReportOneDuplicateLocusTagFmt,
                                            DISC_GENE_DUPLICATE_LOCUS_TAG,
                                            TRUE);
    if (cip != NULL) {
      ValNodeAddPointer (&local_list, 0, cip);
      if (g->adjacent_locus_tag_disc_list != NULL) {
        ValNodeLink (&(cip->subcategories), g->adjacent_locus_tag_disc_list);
      }
    } else if (g->adjacent_locus_tag_disc_list != NULL) {
      ValNodeLink (&local_list, g->adjacent_locus_tag_disc_list);
    }
    g->adjacent_locus_tag_disc_list = NULL;

    /* inconsistent locus tags */
    ValNodeLink (&local_list,
                 ReportInconsistentGlobalDiscrepancyPrefixes (g->locus_tag_list,
                                                              discReportInconsistentLocusTagPrefixFmt,
                                                              DISC_GENE_LOCUS_TAG_INCONSISTENT_PREFIX));
    /* bad formats */
    cip = ReportBadLocusTagFormat (g->locus_tag_list);
    if (cip != NULL) {
      ValNodeAddPointer (&local_list, 0, cip);
    }
  }

  if (g->cds_product_list != NULL) {
    /* report duplicates */
    cip = ReportNonUniqueGlobalDiscrepancy (g->cds_product_list, 
                                            discReportDuplicateProteinIDFmt,
                                            discReportOneDuplicateProteinIDFmt,
                                            DISC_DUP_GENPRODSET_PROTEIN,
                                            TRUE);
    if (cip != NULL) {
      ValNodeAddPointer (&local_list, 0, cip);
    }

    /* report inconsistent IDs */
    ValNodeLink (&local_list,
                 ReportInconsistentGlobalDiscrepancyPrefixes (g->cds_product_list,
                                                              discReportInconsistentProteinIDPrefixFmt,
                                                              DISC_INCONSISTENT_PROTEIN_ID_PREFIX));
  }

  if (g->mrna_product_list != NULL) {
    if (g->missing_locus_tag != NULL) {
      cip = ReportMissingFields (g->mrna_product_list, discReportMissingTranscriptIDFmt, DISC_MISSING_GENPRODSET_TRANSCRIPT_ID);
      if (cip != NULL) {
        ValNodeAddPointer (&local_list, 0, cip);
      }
    }

    cip = ReportNonUniqueGlobalDiscrepancy (g->mrna_product_list,
                                            discReportDuplicateTranscriptIdFmt,
                                            discReportOneDuplicateTranscriptIdFmt,
                                            DISC_DUP_GENPRODSET_TRANSCRIPT_ID,
                                            TRUE);
    if (cip != NULL) {
      ValNodeAddPointer (&local_list, 0, cip);
    }
  }

  /* missing gnl protein IDs */
  cip = ReportMissingFields (g->missing_gnl_list, discReportBadProteinIdFmt, DISC_MISSING_PROTEIN_ID);
  if (cip != NULL) {
    ValNodeAddPointer (&local_list, 0, cip);
  }
  g->gnl_list = ValNodeSort (g->gnl_list, SortVnpByGlobalDiscrepancyString);
  ValNodeLink (&local_list, 
                ReportInconsistentGlobalDiscrepancyStrings (g->gnl_list,
                                                            discReportInconsistentProteinIDPrefixFmt,
                                                            DISC_INCONSISTENT_PROTEIN_ID_PREFIX));


  g->locus_tag_list = FreeGlobalDiscrepancyList (g->locus_tag_list);
  g->missing_locus_tag = FreeGlobalDiscrepancyList (g->missing_locus_tag);
  g->cds_product_list = FreeGlobalDiscrepancyList (g->cds_product_list);
  g->missing_cds_product = FreeGlobalDiscrepancyList (g->missing_cds_product);
  g->mrna_product_list = FreeGlobalDiscrepancyList (g->mrna_product_list);
  g->missing_mrna_product = FreeGlobalDiscrepancyList (g->missing_mrna_product);
  g->missing_gnl_list = FreeGlobalDiscrepancyList (g->missing_gnl_list);
  g->gnl_list = FreeGlobalDiscrepancyList (g->gnl_list);

  /* create discrepancies for inconsistent and missing values from global lists */
  ValNodeLink (&local_list, GetMissingAndInconsistentDiscrepanciesFromGlobalSrcValList (&(g->global_src_qual_vals), &(g->global_srcs), &(g->src_qual_multi_list)));
  /* note - be sure to include local discrepancy reports */
  CollateDiscrepancyReports (&(g->src_qual_repeated_list));
  ValNodeLink (&local_list, g->src_qual_repeated_list);
  g->src_qual_repeated_list = NULL;

  /* create report for feature counts */
  ValNodeLink (&local_list, CreateGlobalFeatureCountReports (&(g->feature_count_list)));

  /* data collected for some tests with global components should not be displayed */
  RemoveUnwantedDiscrepancyItems (&local_list, g->test_config);

  /* group discrepany reports from separate files */
  CollateDiscrepancyReports (&(g->discrepancy_list));

  fprintf (fp, "Discrepancy Report Results\n\n");
  fprintf (fp, "Summary\n");
  WriteDiscrepancyReportSummary (local_list, fp);
  WriteDiscrepancyReportSummary (g->discrepancy_list, fp);

  fprintf (fp, "\n\nDetailed Report\n\n");
  WriteAsnDiscReport (local_list, fp, g->output_config, TRUE);
  local_list = FreeClickableList (local_list);

  WriteAsnDiscReport (g->discrepancy_list, fp, g->output_config, TRUE);
}


/* Barcode Discrepancy Function */

/* 
 * list of names for the individual tests.
 * Note - this array should have eBarcodeTest_LAST elements (see sqnutils.h for value of eBarcodeTest_LAST).
 */
static CharPtr BarcodeTestNames[] = 
{ "Too Short",
  "Missing Primers",
  "Missing Country",
  "Missing Voucher",
  "Too Many Ns",
  "Bad Collection Date",
  "Missing Order Assignment",
  "Low Trace",
  "Frame Shift"
};


extern CharPtr GetBarcodeTestName (Int4 i)
{
  if (i < 0 || i >= sizeof (BarcodeTestNames) / sizeof (CharPtr)) 
  {
    return NULL;
  } 
  else 
  {
    return BarcodeTestNames[i];
  }
}


extern Int4 GetBarcodeTestNumFromBarcodeTestName (CharPtr test_name)
{
  Int4 i;

  if (StringHasNoText (test_name)) {
    return eBarcodeTest_LAST;
  }
  for (i = 0; i < eBarcodeTest_LAST; i++) {
    if (StringICmp (test_name, BarcodeTestNames[i]) == 0) {
      return i;
    }
  }
  return eBarcodeTest_LAST;
}


/* Functions for creating and freeing configurations for the Barcode Tests. */

extern BarcodeTestConfigPtr BarcodeTestConfigNew()
{
  BarcodeTestConfigPtr cfg;
  Int4                 i;

  cfg = (BarcodeTestConfigPtr) MemNew (sizeof (BarcodeTestConfigData));
  for (i = 0; i < eBarcodeTest_LAST; i++)
  {
    cfg->conf_list[i] = TRUE;
  }
  cfg->min_length = 500;
  cfg->min_n_percent = 1.0;
  cfg->require_keyword = TRUE;
  return cfg;
}


extern BarcodeTestConfigPtr BarcodeTestConfigFree (BarcodeTestConfigPtr cfg)
{
  if (cfg != NULL)
  {
    cfg = MemFree (cfg);
  }
  return cfg;
}


/* A BarcodeTestResults lists the Bioseq that the test was performed on,
 * indicates whether each test passed, and give the percentage of Ns
 * (if the value is above the minimum in the configuration).
 */
extern BarcodeTestResultsPtr BarcodeTestResultsNew ()
{
  BarcodeTestResultsPtr res;

  res = (BarcodeTestResultsPtr) MemNew (sizeof (BarcodeTestResultsData));
  MemSet (res, 0, sizeof (BarcodeTestResultsData));
  return res;
}


extern BarcodeTestResultsPtr BarcodeTestResultsFree (BarcodeTestResultsPtr res)
{
  if (res != NULL)
  {
    res = MemFree (res);
  }
  return res;
}


extern BarcodeTestResultsPtr BarcodeTestResultsCopy (BarcodeTestResultsPtr res)
{
  BarcodeTestResultsPtr res_new = NULL;

  if (res != NULL)
  {
    res_new = BarcodeTestResultsNew();
    MemCopy (res_new, res, sizeof (BarcodeTestResultsData));
  }
  return res_new;
}


extern ValNodePtr BarcodeTestResultsListFree (ValNodePtr res_list)
{
  ValNodePtr vnp;

  while (res_list != NULL) 
  {
    vnp = res_list->next;
    res_list->next = NULL;
    res_list->data.ptrvalue = BarcodeTestResultsFree (res_list->data.ptrvalue);
    res_list = ValNodeFree (res_list);
    res_list = vnp;
  }
  return res_list;
}


extern ValNodePtr BarcodeTestResultsExtractPass (ValNodePtr PNTR res_list)
{
  ValNodePtr   vnp, pass_list = NULL;
  BarcodeTestResultsPtr res;

  if (res_list == NULL || *res_list == NULL) {
    return NULL;
  }
  for (vnp = *res_list; vnp != NULL; vnp = vnp->next) {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    if (PassBarcodeTests(res)) {
      vnp->choice = 1;
    }
  }
  pass_list = ValNodeExtractList (res_list, 1);
  return pass_list;
}


/* determines whether barcode tests should be performed on a sequence -
 * no barcode keyword, no barcode tests needed.
 */
extern Boolean HasBARCODETech (BioseqPtr bsp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  Boolean           found = FALSE;
  MolInfoPtr        mip;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
       sdp != NULL && !found;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_molinfo, &dcontext))
  {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL && mip->tech == MI_TECH_barcode)
    {
      found = TRUE;
    }
  }
  return found;
}

/*
 * Finds the MolInfo descriptor for the Bioseq and removes the BARCODE technique.
 * Returns true if the BARCODE technique was present before it was removed.
 */
NLM_EXTERN Boolean RemoveBarcodeTechFromBioseq (BioseqPtr bsp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  Boolean           found = FALSE;
  MolInfoPtr        mip;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
       sdp != NULL && !found;
       sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext))
  {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL && mip->tech == MI_TECH_barcode)
    {
      mip->tech = MI_TECH_unknown;
      found = TRUE;
    }
  }
  return found;
}


NLM_EXTERN Boolean RemoveBarcodeKeywordFromBioseq (BioseqPtr bsp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  Boolean           found = FALSE;
  GBBlockPtr        gb;
  StringConstraint  sc;
  ObjValNodePtr     ovn;

  MemSet (&sc, 0, sizeof (StringConstraint));
  sc.case_sensitive = FALSE;
  sc.match_location = String_location_equals;
  sc.match_text = "BARCODE";
  sc.not_present = FALSE;
  sc.whole_word = FALSE;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &dcontext);
       sdp != NULL && !found;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_genbank, &dcontext))
  {
    gb = (GBBlockPtr) sdp->data.ptrvalue;
    if (gb != NULL)
    {
      found |= RemoveValNodeStringMatch (&(gb->keywords), &sc);
      if (gb->extra_accessions == NULL
          && gb->keywords == NULL
          && gb->source == NULL
          && gb->origin == NULL
          && gb->date == NULL
          && gb->div == NULL
          && gb->taxonomy == NULL
          && gb->entry_date == NULL
          && sdp->extended) {
        ovn = (ObjValNodePtr) sdp;
        ovn->idx.deleteme = TRUE;
      }
    }
  }
  return found;
}


/* 
 * Adds the BARCODE technique to the MolInfo descriptor for the Bioseq.
 * Will create a new MolInfo descriptor for the Bioseq if it doesn't 
 * find one already there.
 */
static void ApplyBarcodeTechToBioseq (BioseqPtr bsp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  Boolean           found = FALSE;
  MolInfoPtr        mip;
  SeqEntryPtr       sep;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
       sdp != NULL && !found;
       sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext))
  {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip == NULL)
    {
      mip = MolInfoNew ();
      sdp->data.ptrvalue = mip;
    }
    mip->tech = MI_TECH_barcode;
    found = TRUE;
  }

  if (!found) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    sdp = CreateNewDescriptor (sep, Seq_descr_molinfo);
    mip = MolInfoNew();
    mip->tech = MI_TECH_barcode;
    sdp->data.ptrvalue = mip;
  }
}


NLM_EXTERN Boolean BioseqHasKeyword (BioseqPtr bsp, CharPtr keyword)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  Boolean           found = FALSE;
  GBBlockPtr        gb;
  ValNodePtr        vnp;
	UserObjectPtr     uop;

  if (StringICmp (keyword, "UNVERIFIED") == 0) 
  {
    /* special case for unverified */
    for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
         sdp != NULL && !found;
         sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext))
    {
      if ((uop = (UserObjectPtr) sdp->data.ptrvalue) != NULL
          && uop->type != NULL
          && StringICmp (uop->type->str, "Unverified") == 0) 
      {
        found = TRUE;
      }
    }
  }

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &dcontext);
       sdp != NULL && !found;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_genbank, &dcontext))
  {
    gb = (GBBlockPtr) sdp->data.ptrvalue;
    if (gb != NULL)
    {
      for (vnp = gb->keywords; vnp != NULL && !found; vnp = vnp->next) 
      {
        if (StringICmp (vnp->data.ptrvalue, keyword) == 0) 
        {
          found = TRUE;
        }
      }
    }
  }
  return found;
}


NLM_EXTERN void ApplyBarcodeKeywordToBioseq (BioseqPtr bsp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  Boolean           found = FALSE;
  GBBlockPtr        gb;
  SeqEntryPtr       sep;

  if (BioseqHasKeyword (bsp, "UNVERIFIED"))
  {
    return;
  }

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &dcontext);
       sdp != NULL && !found;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_genbank, &dcontext))
  {
    gb = (GBBlockPtr) sdp->data.ptrvalue;
    if (gb == NULL)
    {
      gb = GBBlockNew ();
      sdp->data.ptrvalue = gb;
    }
    SetStringsInValNodeStringList (&(gb->keywords), NULL, "BARCODE", ExistingTextOption_add_qual);
    found = TRUE;
  }

  if (!found) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    sdp = CreateNewDescriptor (sep, Seq_descr_genbank);
    gb = GBBlockNew ();
    SetStringsInValNodeStringList (&(gb->keywords), NULL, "BARCODE", ExistingTextOption_add_qual);
    sdp->data.ptrvalue = gb;
  }
}


NLM_EXTERN Boolean BioseqHasBarcodeKeyword (BioseqPtr bsp)
{
  return BioseqHasKeyword (bsp, "BARCODE");
}


NLM_EXTERN void RemoveBarcodeKeywordsFromObjectList (FILE *fp, ValNodePtr object_list)
{
  BioseqPtr  bsp;
  ValNodePtr vnp;
  Char       id_txt[100];

  for (vnp = object_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == OBJ_BIOSEQ && (bsp = (BioseqPtr) vnp->data.ptrvalue) != NULL)
    {
      if (RemoveBarcodeKeywordFromBioseq (bsp))
      {
        if (fp != NULL)
        {
          SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
          fprintf (fp, "%s\n", id_txt);
        }
      }
    }
  }
}  


/* Used for generating Discrepancy Report style data,
 * where Bioseqs are listed separately for each test they fail.
 */
typedef struct barcodesearch {
  ValNodePtr           bioseq_list;
  BarcodeTestConfigPtr cfg;
} BarcodeSearchData, PNTR BarcodeSearchPtr;

NLM_EXTERN Boolean IsIBOL (BioseqPtr bsp)
{
  Boolean           is_ibol = FALSE;
  SeqMgrDescContext context;
  SeqDescPtr        sdp;
  UserObjectPtr     uop;
  UserFieldPtr      curr;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
       sdp != NULL && !is_ibol;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &context)) 
  {
    if ((uop = (UserObjectPtr) sdp->data.ptrvalue) != NULL
        && uop->type != NULL
        && StringICmp (uop->type->str, "StructuredComment") == 0)
    {
      for (curr = uop->data; curr != NULL && !is_ibol; curr = curr->next) 
      {
        if (curr->label != NULL
            && curr->choice == 1
            && StringICmp (curr->label->str, "StructuredCommentPrefix") == 0
                   && StringICmp (curr->data.ptrvalue, "##International Barcode of Life (iBOL)Data-START##") == 0) 
        {
          is_ibol = TRUE;
        }
      }
    }
  }
  return is_ibol;
}


static Boolean HasOrderAssignment (BioseqPtr bsp)
{
  Boolean           has_order = FALSE, is_ibol = FALSE;
  SeqMgrDescContext context;
  SeqDescPtr        sdp;
  UserObjectPtr     uop;
  UserFieldPtr      curr;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
       sdp != NULL && !has_order;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &context)) 
  {
    if ((uop = (UserObjectPtr) sdp->data.ptrvalue) != NULL
        && uop->type != NULL
        && StringICmp (uop->type->str, "StructuredComment") == 0)
    {
      is_ibol = FALSE;
      for (curr = uop->data; curr != NULL && (!has_order || !is_ibol); curr = curr->next) 
      {
        if (curr->label != NULL
            && curr->choice == 1) 
        {
          if (StringICmp (curr->label->str, "Order Assignment") == 0
              && !StringHasNoText (curr->data.ptrvalue))
          {
            has_order = TRUE;
          }
          else if (StringICmp (curr->label->str, "StructuredCommentPrefix") == 0
                   && StringICmp (curr->data.ptrvalue, "##International Barcode of Life (iBOL)Data-START##") == 0) 
          {
            is_ibol = TRUE;
          }
        }
      }
    }
  }
  if (is_ibol && !has_order) {
    return FALSE;
  } else {
    return TRUE;
  }
}


static Boolean HasFrameShift (BioseqPtr bsp)
{
  SeqDescrPtr sdp;
  SeqMgrDescContext context;
  Boolean rval = FALSE;
  UserObjectPtr uop;
  UserFieldPtr  ufp;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
       sdp != NULL && !rval;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &context)) 
  {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL && uop->type != NULL && StringICmp (uop->type->str, "multalin") == 0) 
    {
      ufp = uop->data;
      while (ufp != NULL && !rval) {
        if (ufp->label != NULL 
            && StringICmp (ufp->label->str, "frameshift-nuc") == 0
            && ufp->choice == 1
            && StringICmp (ufp->data.ptrvalue, "fail") == 0) {
          rval = TRUE;
        }
        ufp = ufp->next;
      }
    }
  }
  return rval;
}


typedef Boolean (*BarcodeBioSourceTestFunc) PROTO ((BioSourcePtr));

static Boolean HasForwardAndReversePrimers (BioSourcePtr biop)
{
  Boolean has_forward = FALSE, has_reverse = FALSE;
  PCRReactionPtr primers; 

  if (biop == NULL) return FALSE;

  for (primers = biop->pcr_primers; primers != NULL && (!has_forward || !has_reverse); primers = primers->next)
  {
    if (primers->forward != NULL) 
    {
      has_forward = TRUE;
    }
    if (primers->reverse != NULL)
    {
      has_reverse = TRUE;
    }
  }
  
  return has_forward && has_reverse;

}


static Boolean HasCountry (BioSourcePtr biop)
{
  SubSourcePtr      ssp;
  Boolean           found = FALSE;

  if (biop == NULL || biop->subtype == NULL) return FALSE;

  for (ssp = biop->subtype; ssp != NULL && !found; ssp = ssp->next)
  {
    if (ssp->subtype == SUBSRC_country)
    {
      found = TRUE;
    }
  }

  return found;
}

static Boolean HasVoucher (BioSourcePtr biop)
{
  OrgModPtr mod;
  Boolean   rval = FALSE;

  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) return FALSE;

  for (mod = biop->org->orgname->mod; 
       mod != NULL && !rval; 
       mod = mod->next)
  {
     if (mod->subtype == ORGMOD_specimen_voucher
         || mod->subtype == ORGMOD_bio_material
         || mod->subtype == ORGMOD_culture_collection)
     {
       rval = TRUE;
     }
  }
  return rval;
}


static CharPtr GetDash (CharPtr str)

{
  Char  ch;

  if (str == NULL) return NULL;
  ch = *str;
  while (ch != '\0') {
    if (ch == '-') return str;
    str++;
    ch = *str;
  }

  return NULL;
}

static CharPtr legalMonths [] = {
  "Jan",
  "Feb",
  "Mar",
  "Apr",
  "May",
  "Jun",
  "Jul",
  "Aug",
  "Sep",
  "Oct",
  "Nov",
  "Dec",
  NULL
};

static Int2 daysPerMonth [] = {
  31,
  28,
  31,
  30,
  31,
  30,
  31,
  31,
  30,
  31,
  30,
  31
};

NLM_EXTERN Boolean CollectionDateIsValid (CharPtr name)

{
  Char      ch;
  Int2      dy = 0, dpm = 0, mn = 0;
  Int2      i;
  CharPtr   ptr1, ptr2, month = NULL, day = NULL, year = NULL;
  Char      str [256];
  long int  val;
  Int4      yr = 0;

  if (StringHasNoText (name)) return FALSE;

  StringNCpy_0 (str, name, sizeof (str));
  ptr1 = GetDash (str);
  if (ptr1 != NULL) {
    *ptr1 = '\0';
    ptr1++;
    ptr2 = GetDash (ptr1);
    if (ptr2 != NULL) {
      *ptr2 = '\0';
      ptr2++;
      day = str;
      month = ptr1;
      year = ptr2;
    } else {
      month = str;
      year = ptr1;
    }
  } else {
    year = str;
  }

  if (day != NULL) {
    if (sscanf (day, "%ld", &val) != 1 || val < 1 || val > 31) return FALSE;
    if (StringLen (day) != 2 || !isdigit(day[0]) || !isdigit(day[1])) return FALSE;
    dy = (Int4) val;
  }

  if (month != NULL) {
    for (i = 0; legalMonths [i] != NULL; i++) {
      if (StringCmp (month, legalMonths [i]) == 0) {
        mn = i + 1;
        break;
      }
    }
    if (legalMonths [i] == NULL) return FALSE;
    dpm = daysPerMonth [i];
  }

  if (year != NULL) {
    ptr1 = year;
    ch = *ptr1;
    while (ch != '\0') {
      if (! (IS_DIGIT (ch))) return FALSE;
      ptr1++;
      ch = *ptr1;
    }
    if (sscanf (year, "%ld", &val) == 1) {
      yr = (Int4) val;
      if (val >= 1700 && val < 2100) {
        if (dy > 0 && dpm > 0 && dy > dpm) {
          if (mn != 2 || dy != 29 || (yr % 4) != 0) return FALSE;
        }
        return TRUE;
      }
    }
  }

  return FALSE;
}


/* This mimics a portion of the DatePtr structure,
 * but allows dates with years before 1900 because
 * the year value is Int4 instead of Uint1.
 *   data [0] : Set to 1
 *        [1] - year (- 1900)
 *        [2] - month (1-12)  optional
 *        [3] - day (1-31)     optional
 * Not bothering with time.
 */

typedef struct betterdate {
    Int4 data[8];      /* see box above */
} BetterDateData, PNTR BetterDatePtr;

static BetterDatePtr BetterDateNew()
{
  BetterDatePtr dp;

  dp = (BetterDatePtr) MemNew (sizeof (BetterDateData));
  return dp;
}

static BetterDatePtr BetterDateFree (BetterDatePtr dp)
{
  if (dp != NULL) {
    dp = MemFree (dp);
  }
  return dp;
}


static BetterDatePtr CollectionDateFromString (CharPtr name)
{
  Char      ch;
  Int2      i;
  CharPtr   ptr1, ptr2, month = NULL, day = NULL, year = NULL;
  Char      str [256];
  long int  day_val = 0;
  Int2      month_num = 0;
  long int  val, year_val = 0;
  BetterDatePtr   dp;

  if (StringHasNoText (name)) return NULL;

  StringNCpy_0 (str, name, sizeof (str));
  ptr1 = GetDash (str);
  if (ptr1 != NULL) {
    *ptr1 = '\0';
    ptr1++;
    ptr2 = GetDash (ptr1);
    if (ptr2 != NULL) {
      *ptr2 = '\0';
      ptr2++;
      day = str;
      month = ptr1;
      year = ptr2;
    } else {
      month = str;
      year = ptr1;
    }
  } else {
    year = str;
  }

  if (day != NULL) {
    if (sscanf (day, "%ld", &day_val) != 1 || day_val < 1 || day_val > 31) return NULL;
  }

  if (month != NULL) {
    for (i = 0; legalMonths [i] != NULL; i++) {
      if (StringCmp (month, legalMonths [i]) == 0) {
        month_num = i + 1;
        break;
      }
    }
    if (legalMonths [i] == NULL) return NULL;
  }

  if (year != NULL) {
    ptr1 = year;
    ch = *ptr1;
    while (ch != '\0') {
      if (! (IS_DIGIT (ch))) return NULL;
      ptr1++;
      ch = *ptr1;
    }
    if (sscanf (year, "%ld", &val) == 1) {
      if (val < 1700 || val > 2100) return NULL;
      year_val = val - 1900;
    }
    else
    {
      return NULL;
    }
  }

  dp = BetterDateNew();
  dp->data[0] = 1;
  dp->data[1] = year_val;
  dp->data[2] = month_num;
  dp->data[3] = day_val;
  return dp;
}


NLM_EXTERN Boolean CollectionDateIsInTheFuture (CharPtr name)

{
  DatePtr   dp_now;
  BetterDatePtr dp_coll_date;
  Boolean   rval = FALSE;

  dp_coll_date = CollectionDateFromString (name);
  if (dp_coll_date == NULL) return FALSE;

  if (dp_coll_date->data[1] < 0)
  {
    /* year before 1900 */
    dp_coll_date = BetterDateFree (dp_coll_date);
    return FALSE;
  }

  dp_now = DateCurr();

  /* compare years */
  if (dp_now->data[1] < dp_coll_date->data[1])
  {
    rval = TRUE;
  }
  else if (dp_now->data[1] > dp_coll_date->data[1])
  {
    rval = FALSE;
  }
  /* years are equal - compare months */
  else if (dp_now->data[2] < dp_coll_date->data[2])
  {
    rval = TRUE;
  }
  else if (dp_now->data[2] > dp_coll_date->data[2])
  {
    rval = FALSE;
  }
  /* years and months are equal - compare days */
  else if (dp_now->data[3] < dp_coll_date->data[3])
  {
    rval = TRUE;
  }
  else
  {
    rval = FALSE;
  }

  dp_now = DateFree (dp_now);
  dp_coll_date = BetterDateFree (dp_coll_date);
  return rval;
}


/* collection date is not required, but if present must be valid and in the past */
static Boolean HasCollectionDate (BioSourcePtr biop)
{
  SubSourcePtr ssp;
  Boolean      rval = TRUE;

  if (biop == NULL) {
    return FALSE;
  }
  ssp = biop->subtype;
  while (ssp != NULL && rval) {
    if (ssp->subtype == SUBSRC_collection_date) {
      if (!CollectionDateIsValid(ssp->name) || CollectionDateIsInTheFuture(ssp->name)) {
        rval = FALSE;
      }
    }
    ssp = ssp->next;
  }
  return rval;
}


static Boolean BarcodeGPSOkay (BioSourcePtr biop)
{
  SubSourcePtr ssp;
  Boolean      rval = TRUE;
  CharPtr      country = NULL;
  Boolean      format_ok, lat_in_range, lon_in_range, precision_ok;
  FloatHi      lat, lon;
  Char         buf [256];
  CharPtr      ptr;
  CharPtr      guess;

  if (biop == NULL) {
    return FALSE;
  }
  ssp = biop->subtype;
  /* first find country */
  for (ssp = biop->subtype; ssp != NULL && country == NULL; ssp = ssp->next) 
  {
    if (ssp->subtype == SUBSRC_country) 
    {
      country = ssp->name;
    }
  }
  if (StringContainsBodyOfWater (country)) 
  {
    return TRUE;
  }
  if (country != NULL) 
  {
    StringNCpy_0 (buf, country, sizeof (buf));
    ptr = StringChr (buf, ':');
    if (ptr != NULL) {
      *ptr = '\0';
    }
    country = buf;
  }

  for (ssp = biop->subtype; ssp != NULL && rval; ssp = ssp->next) 
  {
    if (ssp->subtype == SUBSRC_lat_lon) 
    {
      IsCorrectLatLonFormat (ssp->name, &format_ok, &precision_ok, &lat_in_range, &lon_in_range);
      if (!format_ok || !precision_ok || !lat_in_range || !lon_in_range) 
      {
        rval = FALSE;
      } 
      else if (!ParseLatLon (ssp->name, &lat, &lon))
      {
        rval = FALSE;
      }
      else if (country != NULL && IsCountryInLatLonList (country))
      {
        if (!(TestLatLonForCountry (country, lat, lon))) 
        {
          guess = GuessCountryForLatLon (lat, lon);
          if (StringHasNoText (guess) || !CountryBoxesOverlap (country, guess)) {
            rval = FALSE;
          }
        }
      }
    }
  }

  return rval;
}


static Boolean BarcodeBioSourceTest (BioseqPtr bsp, BarcodeBioSourceTestFunc test_func, Boolean require_keyword)
{
  SeqDescrPtr       sdp;
  BioSourcePtr      biop;
  SeqMgrDescContext context;
  Boolean           found = FALSE;
  
  if (bsp == NULL || ISA_aa (bsp->mol) || (require_keyword && !HasBARCODETech (bsp)) || test_func == NULL)
  {
    return FALSE;
  }
    
  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
       sdp != NULL && !found;
       sdp = SeqMgrGetNextDescriptor (bsp,sdp, Seq_descr_source, &context))
  {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    found = test_func(biop);
  }

  return !found;
}


static void BarcodeBioSourceTestCallback (BioseqPtr bsp, Pointer userdata, BarcodeBioSourceTestFunc test_func)
{  
  BarcodeSearchPtr bsd;

  if (bsp == NULL || ISA_aa (bsp->mol)
      || (bsd = (BarcodeSearchPtr) userdata) == NULL
      || bsd->cfg == NULL
      || test_func == NULL)
  {
    return;
  }
  
  if (BarcodeBioSourceTest (bsp, test_func, bsd->cfg->require_keyword))
  {
    ValNodeAddPointer (&(bsd->bioseq_list), OBJ_BIOSEQ, bsp);
  }
}


static void FindBadGPS (BioseqPtr bsp, Pointer userdata)
{
  BarcodeBioSourceTestCallback (bsp, userdata, BarcodeGPSOkay);
}



static void 
BarcodeTestForSeqEntry 
(SeqEntryPtr          sep,
 ValNodePtr PNTR      discrepancy_list, 
 VisitBioseqsFunc     callback, 
 CharPtr              fmt,
 BarcodeTestConfigPtr cfg)
{
  BarcodeSearchData bsd;

  bsd.bioseq_list = NULL;
  bsd.cfg = cfg;
  if (bsd.cfg == NULL)
  {
    bsd.cfg = BarcodeTestConfigNew ();
  }
  VisitBioseqsInSep (sep, &bsd, callback);

  if (bsd.bioseq_list != NULL)
  {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (0, fmt, bsd.bioseq_list));
  }
  if (cfg != bsd.cfg)
  {
    bsd.cfg = BarcodeTestConfigFree (bsd.cfg);
  }
}

 
static void BarcodePercentNDiscrepanciesForSeqEntry (ValNodePtr results, ValNodePtr PNTR discrepancy_list, FloatLo min_n_percent)
{
  BarcodeTestResultsPtr res;
  ValNodePtr subcategories = NULL, bioseq_list = NULL, vnp;
  ClickableItemPtr cip;
  CharPtr fmt = "Sequence has %.1f%% percent Ns";
  CharPtr top_fmt = "%d sequences have > %.1f%% Ns";

  for (vnp = results; vnp != NULL; vnp = vnp->next)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    if (res->n_percent < min_n_percent) {
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      MemSet (cip, 0, sizeof (ClickableItemData));
      cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + 5));
      sprintf (cip->description, fmt, res->n_percent);
      ValNodeAddPointer (&bioseq_list, OBJ_BIOSEQ, res->bsp);
      ValNodeAddPointer (&(cip->item_list), OBJ_BIOSEQ, res->bsp);
      ValNodeAddPointer (&subcategories, 0, cip);
    }
  }

  if (bioseq_list != NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (top_fmt) + 10));
    sprintf (cip->description, fmt, ValNodeLen (bioseq_list), min_n_percent);
    cip->item_list = bioseq_list;
    cip->subcategories = subcategories;
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


static void GetBarcodeDiscrepanciesForSeqEntry (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list, BarcodeTestConfigPtr cfg)
{
  ValNodePtr results, vnp;
  ValNodePtr PNTR lists;
  BarcodeTestResultsPtr res;
  Int4 i;
  CharPtr fmts[] = {"%d sequences are shorter than 500 nucleotides",
                       "%d sequences are missing forward and/or reverse primers",
                       "%d sequences are missing country", 
                       "%d sequences are missing specimen voucher",
                       NULL,
                       "%d sequences have invalid collection date",
                       "%d sequences are missing order assignment",
                       "%d sequences have low trace",
                       "%d sequences have frameshift" };



  if (cfg == NULL) return;

  results = GetBarcodePassFail(sep, cfg);

  lists = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * eBarcodeTest_LAST);
  MemSet (lists, 0, sizeof (ValNodePtr) * eBarcodeTest_LAST);
  for (vnp = results; vnp != NULL; vnp = vnp->next) {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    for (i = 0; i < eBarcodeTest_LAST; i++) {
      if (cfg->conf_list[i] && res->failed_tests[i] && fmts[i] != NULL) {
        ValNodeAddPointer (&(lists[i]), OBJ_BIOSEQ, res->bsp);
      }
    }
  }
  for (i = 0; i < eBarcodeTest_LAST; i++) {
    if (cfg->conf_list[i] && lists[i] != NULL) {
      if (fmts[i] != NULL) {
        ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (0, fmts[i], lists[i]));
      }
    }
  }
  lists = MemFree(lists);

  if (cfg->conf_list[eBarcodeTest_PercentN])
  {
    BarcodePercentNDiscrepanciesForSeqEntry (sep, discrepancy_list, cfg->min_n_percent);
  }

  results = BarcodeTestResultsListFree(results);
}


extern ValNodePtr GetBarcodeDiscrepancies (ValNodePtr sep_list, BarcodeTestConfigPtr cfg)
{
  ValNodePtr    vnp, discrepancy_list = NULL;
  SeqEntryPtr   sep;
  BarcodeTestConfigPtr local_cfg;

  if (cfg == NULL) 
  {
    local_cfg = BarcodeTestConfigNew();
  }
  else
  {
    local_cfg = cfg;
  }

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next)
  {
    sep = vnp->data.ptrvalue;
    GetBarcodeDiscrepanciesForSeqEntry (sep, &discrepancy_list, local_cfg);
  }

  if (local_cfg != cfg)
  {
    local_cfg = BarcodeTestConfigFree (local_cfg);
  }

  /* normalize the discrepancy levels so that they will be correctly displayed */
  SetDiscrepancyLevels (discrepancy_list, 0);
  return discrepancy_list;  
}


/* This section is used for generating the "Failure Report" and "Compliance Report". */

typedef struct barcodebioseqsearch {
  ValNodePtr           results_list;
  BarcodeTestConfigPtr cfg;
  Boolean              collect_positives;
} BarcodeBioseqSearchData, PNTR BarcodeBioseqSearchPtr;


extern Boolean IsBarcodeID (SeqIdPtr sip)
{
  DbtagPtr dbt;

  if (sip == NULL) return FALSE;
  if (sip->choice != SEQID_GENERAL) return FALSE;
  dbt = (DbtagPtr) sip->data.ptrvalue;
  if (dbt == NULL) return FALSE;
  if (StringICmp (dbt->db, "uoguelph") == 0) return TRUE;
  return FALSE;
}

#define cMaxBarcodeIDStringLen 200
#define cMaxGenbankIDStringLen 200

extern CharPtr BarcodeTestBarcodeIdString (BioseqPtr bsp)
{
  SeqIdPtr barcode_id;
  Char     barcode_id_str[cMaxBarcodeIDStringLen];

  if (bsp == NULL) return NULL;

  barcode_id = bsp->id;
  while (barcode_id != NULL && !IsBarcodeID (barcode_id))
  {
    barcode_id = barcode_id->next;
  }

  if (barcode_id == NULL) 
  {
    barcode_id = bsp->id;
    while (barcode_id != NULL && barcode_id->choice != SEQID_LOCAL)
    {
      barcode_id = barcode_id->next;
    }
  }

  if (barcode_id == NULL) 
  {
    sprintf (barcode_id_str, "NO");
  }
  else
  {
    SeqIdWrite (barcode_id, barcode_id_str, PRINTID_FASTA_SHORT, sizeof (barcode_id_str) - 1);
  }
  return StringSave (barcode_id_str);
}

extern CharPtr BarcodeTestGenbankIdString (BioseqPtr bsp)
{
  SeqIdPtr genbank_id;
  Char     genbank_id_str[cMaxGenbankIDStringLen];
  CharPtr  src, dst;

  genbank_id = bsp->id;
  while (genbank_id != NULL && genbank_id->choice != SEQID_GENBANK)
  {
    genbank_id = genbank_id->next;
  }
  if (genbank_id == NULL) 
  {
    sprintf (genbank_id_str, "NO");
  }
  else
  {
    SeqIdWrite (genbank_id, genbank_id_str, PRINTID_FASTA_SHORT, sizeof (genbank_id_str) - 1);
    if (StringNICmp (genbank_id_str, "gb|", 3) == 0) {
      src = genbank_id_str + 3;
      dst = genbank_id_str;
      while (*src != 0) {
        *dst = *src;
        dst++;
        src++;
      }
      dst[0] = 0;
    }
    if (genbank_id_str[StringLen (genbank_id_str) - 1] == '|') {
      genbank_id_str[StringLen (genbank_id_str) - 1] = 0;
    }
  }
  return StringSave (genbank_id_str);
}



NLM_EXTERN CharPtr GetBarcodeTestFailureReasons (BarcodeTestResultsPtr res)
{
  Int4             i, msg_len = 0;
  Boolean          any_failed = FALSE;
  CharPtr          msg;
  Char             pct[10];

  if (res == NULL || res->bsp == NULL) return NULL;

  for (i = 0; i < eBarcodeTest_LAST; i++)
  {
    if (res->failed_tests[i]) 
    {
      msg_len += StringLen (GetBarcodeTestName (i)) + 2;
      if (i == eBarcodeTest_PercentN)
      {
        msg_len += 5;
      }
      any_failed = TRUE;
    }
  }
  if (!any_failed) return NULL;
   
  msg = (CharPtr) MemNew (sizeof (Char) * msg_len);
  for (i = 0; i < eBarcodeTest_LAST; i++)
  {
    if (res->failed_tests[i])
    {
      StringCat (msg, GetBarcodeTestName(i));
      if (i == eBarcodeTest_PercentN) 
      {
        sprintf (pct, ":%.1f%%", res->n_percent);
        StringCat (msg, pct);
      }
      StringCat (msg, ",");
    }
  }
  /* remove trailing comma */
  msg[StringLen(msg) - 1] = 0;
          
  return msg;

}


static CharPtr SummaryTextFromBarcodeTestResults (BarcodeTestResultsPtr res)
{
  Int4             i, msg_len = 0;
  Boolean          any_failed = FALSE;
  CharPtr          msg, genbank_id, barcode_id;
  Char             pct[10];

  if (res == NULL || res->bsp == NULL) return NULL;

  for (i = 0; i < eBarcodeTest_LAST; i++)
  {
    if (res->failed_tests[i]) 
    {
      msg_len += StringLen (GetBarcodeTestName (i)) + 2;
      if (i == eBarcodeTest_PercentN)
      {
        msg_len += 6;
      }
      any_failed = TRUE;
    }
  }
  if (!any_failed) return NULL;

  genbank_id = BarcodeTestGenbankIdString (res->bsp);
  barcode_id = BarcodeTestBarcodeIdString (res->bsp);
  
  msg_len += StringLen (genbank_id) + StringLen (barcode_id) + 2;
 
  msg = (CharPtr) MemNew (sizeof (Char) * msg_len);
  sprintf (msg, "%s\t%s\t", barcode_id, genbank_id);
  for (i = 0; i < eBarcodeTest_LAST; i++)
  {
    if (res->failed_tests[i])
    {
      StringCat (msg, GetBarcodeTestName(i));
      if (i == eBarcodeTest_PercentN) 
      {
        sprintf (pct, ":%.1f%%", res->n_percent);
        StringCat (msg, pct);
      }
      StringCat (msg, ",");
    }
  }
  /* remove trailing comma */
  msg[StringLen(msg) - 1] = 0;
          
  return msg;
}


extern Boolean PassBarcodeTests (BarcodeTestResultsPtr res)
{
  if (res == NULL
      || res->failed_tests[eBarcodeTest_Length]
      || res->failed_tests[eBarcodeTest_Primers]
      || res->failed_tests[eBarcodeTest_Country]
      || res->failed_tests[eBarcodeTest_SpecimenVoucher]
      || res->failed_tests[eBarcodeTest_PercentN]
      || res->failed_tests[eBarcodeTest_CollectionDate]
      || res->failed_tests[eBarcodeTest_OrderAssignment]
      || res->failed_tests[eBarcodeTest_LowTrace]
      || res->failed_tests[eBarcodeTest_FrameShift])
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


/* NOTE - this no longer performs the low trace test - that test needs to be done for the seq-entry as a whole */
static BarcodeTestResultsPtr BarcodeTestResultsForBioseq (BioseqPtr bsp, BarcodeTestConfigPtr cfg)
{
  BarcodeTestResultsPtr res = NULL;

  if (bsp == NULL || ISA_aa (bsp->mol) || cfg == NULL || (cfg->require_keyword && !HasBARCODETech (bsp)))
  {
    return NULL;
  }

  res = BarcodeTestResultsNew ();

  res->bsp = bsp;

  if (bsp->length < cfg->min_length && cfg->conf_list[eBarcodeTest_Length])
  {
    res->failed_tests[eBarcodeTest_Length] = TRUE;
  } 

  if (cfg->conf_list[eBarcodeTest_Primers])
  {
    res->failed_tests[eBarcodeTest_Primers] = BarcodeBioSourceTest(bsp, HasForwardAndReversePrimers, cfg->require_keyword);
  }
  if (cfg->conf_list[eBarcodeTest_Country])
  {
    res->failed_tests[eBarcodeTest_Country] = BarcodeBioSourceTest(bsp, HasCountry, cfg->require_keyword);
  }
  if (cfg->conf_list[eBarcodeTest_SpecimenVoucher])
  {
    res->failed_tests[eBarcodeTest_SpecimenVoucher] = BarcodeBioSourceTest(bsp, HasVoucher, cfg->require_keyword);
  }
  if (cfg->conf_list[eBarcodeTest_CollectionDate]) 
  {
    res->failed_tests[eBarcodeTest_CollectionDate] = BarcodeBioSourceTest(bsp, HasCollectionDate, cfg->require_keyword);
  }
  if (cfg->conf_list[eBarcodeTest_OrderAssignment]) 
  {
    res->failed_tests[eBarcodeTest_OrderAssignment] = !HasOrderAssignment (bsp);
  }
  if (cfg->conf_list[eBarcodeTest_FrameShift])
  {
    res->failed_tests[eBarcodeTest_FrameShift] = IsIBOL(bsp) && HasFrameShift (bsp);
  }

  if (cfg->conf_list[eBarcodeTest_PercentN])
  {
    res->n_percent = PercentNInBioseq (bsp, TRUE);
    res->failed_tests[eBarcodeTest_PercentN] = (Boolean)(res->n_percent > cfg->min_n_percent);
  }
 
  return res;
}


static void DoBarcodeTestsExceptLowTrace (BioseqPtr bsp, Pointer userdata)
{
  BarcodeBioseqSearchPtr sp;
  BarcodeTestResultsPtr  res = NULL;

  if (bsp == NULL || ISA_aa (bsp->mol) 
      || (sp = (BarcodeBioseqSearchPtr) userdata) == NULL
      || sp->cfg == NULL
      || (sp->cfg->require_keyword && !HasBARCODETech (bsp)))
  {
    return;
  }

  res = BarcodeTestResultsForBioseq (bsp, sp->cfg);
  if (res == NULL) return;

  ValNodeAddPointer (&(sp->results_list), 0, res);
}


#ifdef OS_MSWIN
#include <undefwin.h>
#include <windows.h>

NLM_EXTERN Int4 RunSilent(const char *cmdline) {
    int status = -1;

    STARTUPINFO         StartupInfo;
    PROCESS_INFORMATION ProcessInfo;

    DWORD dwCreateFlags;

#ifndef COMP_METRO
    /* code warrior headers do not have this, so comment out to allow compilation */
    _flushall();
#endif

    /* Set startup info */
    memset(&StartupInfo, 0, sizeof(StartupInfo));
    StartupInfo.cb          = sizeof(STARTUPINFO);
    StartupInfo.dwFlags     = STARTF_USESHOWWINDOW;
    StartupInfo.wShowWindow = SW_HIDE;
    dwCreateFlags           = CREATE_NEW_CONSOLE;

    /* Run program */
    if (CreateProcess(NULL, (LPSTR)cmdline, NULL, NULL, FALSE,
                      dwCreateFlags, NULL, NULL, &StartupInfo, &ProcessInfo))
    {
        /* wait running process */
        DWORD exitcode = -1;
        WaitForSingleObject(ProcessInfo.hProcess, INFINITE);
        GetExitCodeProcess(ProcessInfo.hProcess, &exitcode);
        status = exitcode;
        CloseHandle(ProcessInfo.hProcess);
        CloseHandle(ProcessInfo.hThread);
    }
    else
    {
	DWORD dw = GetLastError();
	/* check for common errors first */
	if(dw == ERROR_FILE_NOT_FOUND)
	    Message(MSG_ERROR, "CreateProcess() failed: file not found.");
	else
	    /* generic error message */
	    Message(MSG_ERROR, "CreateProcess() failed, error code %d.",
		    (int)dw);
    }

    return status;
}
#endif

static CharPtr tracefetchcmd = NULL;

static void FillInMissingTraces (ValNodePtr trace_check_list)
{
  Char     path_in [PATH_MAX];
  Char     path_out [PATH_MAX];
  FILE     *fp;
  Char     id_txt[255];
  Char     cmmd [256];
  ValNodePtr vnp;
  BarcodeTestResultsPtr res;
  ReadBufferData        rbd;
  CharPtr               line, cp;

  if (tracefetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "TRACECOUNT", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	tracefetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (tracefetchcmd == NULL) return;

  TmpNam (path_in);
  fp = FileOpen (path_in, "w");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to open temporary file %s, unable to get trace results", path_in);
  } else {
    /* make list of accessions to check */
    for (vnp = trace_check_list; vnp != NULL; vnp = vnp->next) {
      res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
      if (res != NULL) {
        SeqIdWrite (SeqIdFindBest (res->bsp->id, SEQID_GENBANK), id_txt, PRINTID_TEXTID_ACC_ONLY, sizeof (id_txt) - 1);
        fprintf (fp, "%s\n", id_txt);
      }
    }
    FileClose (fp);
    TmpNam (path_out);
    /* launch script */
#ifdef OS_UNIX
    sprintf (cmmd, "%s -i %s -o %s", tracefetchcmd, path_in, path_out);
    system (cmmd);
#endif
#ifdef OS_MSWIN
    sprintf (cmmd, "%s -i %s -o %s", tracefetchcmd, path_in, path_out);
    RunSilent (cmmd);
#endif
    /* read results */
    fp = FileOpen (path_out, "r");
    if (fp == NULL) {
      Message (MSG_ERROR, "Unable to open temporary file %s for results", path_out);
    } else {
      rbd.current_data = NULL;
      rbd.fp = fp;

      line = AbstractReadFunction (&rbd); 
      vnp = trace_check_list;
      res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
      SeqIdWrite (SeqIdFindBest (res->bsp->id, SEQID_GENBANK), id_txt, PRINTID_TEXTID_ACC_ONLY, sizeof (id_txt) - 1);

      while (line != NULL && line[0] != EOF && vnp != NULL) {
        if (!StringHasNoText (line)) {
          cp = StringChr (line, '\t');
          if (cp != NULL) {
            *cp = 0;
            while (StringCmp (id_txt, line) != 0 && vnp != NULL) {
              if (res->num_trace < 2) {
                res->failed_tests[eBarcodeTest_LowTrace] = TRUE;
              }
              vnp = vnp->next;
              res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
              SeqIdWrite (SeqIdFindBest (res->bsp->id, SEQID_GENBANK), id_txt, PRINTID_TEXTID_ACC_ONLY, sizeof (id_txt) - 1);
            }
            if (vnp != NULL) {
              res->num_trace++;
            }
          }
        }
        line = MemFree (line);
        line = AbstractReadFunction (&rbd);
      }
      while (vnp != NULL) {
        res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
        if (res->num_trace < 2) {
          res->failed_tests[eBarcodeTest_LowTrace] = TRUE;
        }
        vnp = vnp->next;
      }
          
      FileClose (fp);
      FileRemove (path_out);
    }
    FileRemove (path_in);
  }
}


extern ValNodePtr GetBarcodePassFail (SeqEntryPtr sep, BarcodeTestConfigPtr cfg)
{
  BarcodeBioseqSearchData sd;
  ValNodePtr              vnp;
  BarcodeTestResultsPtr   res;
  SeqDescPtr              sdp;
  SeqMgrDescContext       context;
  UserObjectPtr           uop;
  UserFieldPtr            ufp;
  ObjectIdPtr             oip;
  Boolean                 has_low_trace, has_object;
  int                     num_trace = 0;
  ValNodeBlock            trace_check_list;

  if (cfg == NULL)
  {
    sd.cfg = BarcodeTestConfigNew();
  }
  else
  {
    sd.cfg = cfg;
  }

  sd.results_list = NULL;

  VisitBioseqsInSep (sep, &sd, DoBarcodeTestsExceptLowTrace);
  InitValNodeBlock (&trace_check_list, NULL);

  /* now do low trace test */
  /* first, loop through list - if bioseq has submission object with trace statement,
   * get result from that.  otherwise add to list. */
  for (vnp = sd.results_list; vnp != NULL; vnp = vnp->next) {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    if (res != NULL) {
      /* look for user object */
      has_low_trace = FALSE;
      has_object = FALSE;
      for (sdp = SeqMgrGetNextDescriptor (res->bsp, NULL, Seq_descr_user, &context);
            sdp != NULL && !has_low_trace;
            sdp = SeqMgrGetNextDescriptor (res->bsp, sdp, Seq_descr_user, &context)) {
        uop = (UserObjectPtr) sdp->data.ptrvalue;
        if (uop != NULL && uop->type != NULL && StringICmp (uop->type->str, "Submission") == 0) {
          for (ufp = uop->data; ufp != NULL && !has_low_trace; ufp = ufp->next) {
            oip = ufp->label;
            if (oip != NULL && StringCmp (oip->str, "AdditionalComment") == 0) {
              if ( sscanf (ufp->data.ptrvalue, "Traces: %d", &num_trace) == 1) {
                res->num_trace = num_trace;
                if (num_trace < 2) {
                  has_low_trace = TRUE;
                }
                has_object = TRUE;
              }
            }
          }
        }
      }
      if (has_low_trace) {
        res->failed_tests[eBarcodeTest_LowTrace] = TRUE;
      } else if (!has_object) {
        ValNodeAddPointerToEnd (&trace_check_list, 0, res);
      }
    }
  }
  
  /* then put IDs in list, use script to collect from trace, add to results. */
  if (trace_check_list.head != NULL) {
    FillInMissingTraces (trace_check_list.head);
    /* NOTE - do NOT free barcode result data, since this list points to data in sd.results list */
    trace_check_list.head = ValNodeFree (trace_check_list.head);
  }

  if (sd.cfg != cfg)
  {
    sd.cfg = BarcodeTestConfigFree (sd.cfg);
  }
  return sd.results_list;
}


/* Report lists each Bioseq and whether the Bioseq passed all tests
 * or failed at least one.
 */
extern void WriteBarcodeTestComplianceEx (FILE *fp, ValNodePtr results_list, Boolean low_trace_fail)
{
  BarcodeTestResultsPtr res;
  ValNodePtr            vnp;
  CharPtr               barcode_id, genbank_id;
  Boolean               pass;

  if (fp == NULL) return;

  for (vnp = results_list; vnp != NULL; vnp = vnp->next)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    barcode_id = BarcodeTestBarcodeIdString (res->bsp);
    genbank_id = BarcodeTestGenbankIdString (res->bsp);
    pass = PassBarcodeTests (res);
    fprintf (fp, "%s\t%s\t%s\n", barcode_id, genbank_id, pass ? "PASS" : "FAIL");
    barcode_id = MemFree (barcode_id);
    genbank_id = MemFree (genbank_id);
  }
}


extern void WriteBarcodeTestCompliance (FILE *fp, ValNodePtr results_list)
{
  WriteBarcodeTestComplianceEx (fp, results_list, FALSE);
}


/* Report lists each Bioseq and whether the Bioseq passed all tests
 * or failed at least one.
 */
extern void WriteBarcodeTestComprehensive (FILE *fp, ValNodePtr results_list)
{
  BarcodeTestResultsPtr res;
  ValNodePtr            vnp;
  CharPtr               barcode_id, genbank_id, reason;

  if (fp == NULL) return;

  for (vnp = results_list; vnp != NULL; vnp = vnp->next)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    barcode_id = BarcodeTestBarcodeIdString (res->bsp);
    genbank_id = BarcodeTestGenbankIdString (res->bsp);
    reason = GetBarcodeTestFailureReasons (res);
    fprintf (fp, "%s\t%s\t%s\t%s\n", barcode_id, genbank_id, 
                                 PassBarcodeTests (res) ? "PASS" : "FAIL",
                                 reason == NULL ? "" : reason);
    barcode_id = MemFree (barcode_id);
    genbank_id = MemFree (genbank_id);
    reason = MemFree (reason);
  }
}


/* Create a tag table for updates */
extern void WriteBarcodeTagTable (FILE *fp, ValNodePtr results_list)
{
  BarcodeTestResultsPtr res;
  ValNodePtr            vnp;
  CharPtr               barcode_id, genbank_id;

  if (fp == NULL) return;

  for (vnp = results_list; vnp != NULL; vnp = vnp->next)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    barcode_id = BarcodeTestBarcodeIdString (res->bsp);
    genbank_id = BarcodeTestGenbankIdString (res->bsp);
    fprintf (fp, "%s\t%s\t\n", genbank_id, barcode_id);
    barcode_id = MemFree (barcode_id);
    genbank_id = MemFree (genbank_id);
  }
}


/* Report lists the individual tests that each Bioseq failed. */
extern void WriteBarcodeDiscrepancies (FILE *fp, ValNodePtr results_list)
{
  BarcodeTestResultsPtr res;
  ValNodePtr            vnp;
  CharPtr               msg;

  if (fp == NULL) return;

  for (vnp = results_list; vnp != NULL; vnp = vnp->next)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    
    msg = SummaryTextFromBarcodeTestResults (res);
    fprintf (fp, "%s\n", msg);
    msg = MemFree (msg);
  }
}


static CharPtr FailureTextFromBarcodeTestResults (BarcodeTestResultsPtr res)
{
  Int4             i, msg_len = 0;
  Boolean          any_failed = FALSE;
  CharPtr          msg, genbank_id, barcode_id;

  if (res == NULL || res->bsp == NULL) return NULL;

  for (i = 0; i < eBarcodeTest_LAST; i++)
  {
    if (res->failed_tests[i]) 
    {
      msg_len += StringLen (GetBarcodeTestName (i)) + 2;
      any_failed = TRUE;
    }
  }
  if (!any_failed) return NULL;

  genbank_id = BarcodeTestGenbankIdString (res->bsp);
  barcode_id = BarcodeTestBarcodeIdString (res->bsp);
  
  msg_len += StringLen (genbank_id) + StringLen (barcode_id) + 2;
 
  msg = (CharPtr) MemNew (sizeof (Char) * msg_len);
  sprintf (msg, "%s\t%s\t", barcode_id, genbank_id);
  for (i = 0; i < eBarcodeTest_LAST; i++)
  {
    if (res->failed_tests[i])
    {
      StringCat (msg, GetBarcodeTestName(i));
      StringCat (msg, ",");
    }
  }
  /* remove trailing comma */
  msg[StringLen(msg) - 1] = 0;
          
  return msg;
}


extern void WriteBarcodeFailureReport (FILE *fp, ValNodePtr results_list)
{
  ValNodePtr vnp;
  BarcodeTestResultsPtr res;
  CharPtr msg;

  if (fp == NULL) return;

  for (vnp = results_list; vnp != NULL; vnp = vnp->next)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    msg = FailureTextFromBarcodeTestResults (res);
    if (msg != NULL) {
      fprintf (fp, "%s\n", msg);
      msg = MemFree (msg);
    }
  }
}


static void BarcodeValPrintStr (FILE *fp, CharPtr fmt, CharPtr str)
{
  if (fp == NULL) {
    if (fmt == NULL) {
      ErrPost (0, 0, str);
    } else {
      ErrPost (0, 0, fmt, str);
    }
  } else {
    if (fmt == NULL) {
      fprintf (fp, "%s", str);
    } else {
      fprintf (fp, fmt, str);
    }
  }
}


NLM_EXTERN Boolean 
BarcodeValidateOneSeqEntry 
(FILE *ofp,
 SeqEntryPtr sep,
 Boolean show_all,
 Boolean use_xml,
 Boolean show_header,
 CharPtr xml_header_text)

{
  BarcodeTestConfigPtr  cfg;
  ValNodePtr            pass_fail_list = NULL, vnp;
  BarcodeTestResultsPtr res;
  Char                  id_buf[255];
  Char                  num_buf[255];
  CharPtr               reason;
  Boolean               any_failures = FALSE;
  Int4                  i;

  if (sep == NULL) return FALSE;

  cfg = BarcodeTestConfigNew();
  cfg->require_keyword = FALSE;
  pass_fail_list = GetBarcodePassFail (sep, cfg);

  for (vnp = pass_fail_list; vnp != NULL && !any_failures; vnp = vnp->next) {
    if (!PassBarcodeTests (vnp->data.ptrvalue)) {
      any_failures = TRUE;
    }
  }

  if (pass_fail_list != NULL && (show_all || any_failures)) {
    if (use_xml) {
      if (show_header) {
        BarcodeValPrintStr (ofp, "<%s>\n", xml_header_text);
      }
      for (vnp = pass_fail_list; vnp != NULL; vnp = vnp->next) {
        res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
        SeqIdWrite (SeqIdFindBest (res->bsp->id, SEQID_GENBANK), id_buf, PRINTID_REPORT, sizeof (id_buf) - 1);
        for (i = 0; i < eBarcodeTest_LAST; i++) {
          if (res->failed_tests[i]) {
            BarcodeValPrintStr (ofp, "  <message severity=\"WARNING\" seq-id=\"%s\">", id_buf);
            BarcodeValPrintStr (ofp, " %s", GetBarcodeTestName (i));
            if (i == eBarcodeTest_PercentN) {
              sprintf (num_buf, ":%.1f%%", res->n_percent);
              BarcodeValPrintStr (ofp, NULL, num_buf);
            } else if (i == eBarcodeTest_Length) {
              sprintf (num_buf, ":Length should be at least %d", cfg->min_length);
              BarcodeValPrintStr (ofp, NULL, num_buf);
            }
            BarcodeValPrintStr (ofp, NULL, "</message>\n");
          }
        }
      }
      if (show_all) {
        for (vnp = pass_fail_list; vnp != NULL; vnp = vnp->next) {
          res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
          SeqIdWrite (SeqIdFindBest (res->bsp->id, SEQID_GENBANK), id_buf, PRINTID_REPORT, sizeof (id_buf) - 1);
          reason = GetBarcodeTestFailureReasons (res);
          BarcodeValPrintStr (ofp, "  <message severity=\"INFO\" seq-id=\"%s\">", id_buf);
          if (PassBarcodeTests(res)) {
            BarcodeValPrintStr (ofp, NULL, "PASS");
          } else {
            BarcodeValPrintStr (ofp, "FAIL (%s)", reason == NULL ? "" : reason);
          }
          BarcodeValPrintStr (ofp, NULL, "</message>\n");
          reason = MemFree (reason);
        }
      }
    } else {
      if (show_header) { 
        if (ofp == NULL) {
          ErrPost (0, 0, "\n\nBarcode Validation Test Results\n");
        } else {
          fprintf (ofp, "\n\nBarcode Validation Test Results\n");
        }
        if (show_all) {
          if (ofp == NULL) {
            ErrPost (0, 0, "ID\tPassed?\tReason\n");
          } else {
            fprintf (ofp, "ID\tPassed?\tReason\n");
          }
        } else {
          if (ofp == NULL) {
            ErrPost (0, 0, "ID\tReason\n");
          } else {
            fprintf (ofp, "ID\tReason\n");
          }
        }
      }
      for (vnp = pass_fail_list; vnp != NULL; vnp = vnp->next) {
        res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
        SeqIdWrite (SeqIdFindBest (res->bsp->id, SEQID_GENBANK), id_buf, PRINTID_REPORT, sizeof (id_buf) - 1);
        reason = GetBarcodeTestFailureReasons (res);
        if (show_all) {
          if (ofp == NULL) {
            ErrPost (0, 0, "%s\t%s\t%s\n", id_buf, 
                                      PassBarcodeTests (res) ? "PASS" : "FAIL",
                                      reason == NULL ? "" : reason);
          } else {
            fprintf (ofp, "%s\t%s\t%s\n", id_buf, 
                                      PassBarcodeTests (res) ? "PASS" : "FAIL",
                                      reason == NULL ? "" : reason);
          }
        } else {
          if (!PassBarcodeTests (res)) {
            if (ofp == NULL) {
              ErrPost (0, 0, "%s\t%s\n", id_buf, 
                                        reason == NULL ? "" : reason);
            } else {
              fprintf (ofp, "%s\t%s\n", id_buf, 
                                        reason == NULL ? "" : reason);
            }
          }
        }      
        reason = MemFree (reason);
      }
    }
    pass_fail_list = BarcodeTestResultsListFree (pass_fail_list);      
  }
  cfg = BarcodeTestConfigFree (cfg);

  return !any_failures;
}


static void LIBCALLBACK CountPolymorphismProc (CharPtr sequence, Pointer userdata)
{
  Int4Ptr p_i;
  CharPtr cp;

  if (sequence == NULL || userdata == NULL) return;
  p_i = (Int4Ptr) userdata;

  for (cp = sequence; *cp != 0; cp++)
  {
    if (*cp != 'N' && *cp != 'A' && *cp != 'T' && *cp != 'G' && *cp != 'C')
    {
      (*p_i) ++;
    }
  }
}


extern Int4 CountPolymorphismsInBioseq (BioseqPtr bsp)
{
  Int4 num_p = 0;
  
  if (bsp->length == 0 || IsDeltaSeqWithFarpointers (bsp)) return 0;

  /* if delta sequence, ignore Ns from gaps */

  SeqPortStream (bsp, 0, (Pointer) &num_p, CountPolymorphismProc);

  return num_p;
}


/* Removes Barcode tech from all Bioseqs in supplied list */
extern void RemoveBarcodeTech (FILE *fp, ValNodePtr results_list)
{
  BarcodeTestResultsPtr res;
  ValNodePtr            vnp;
  CharPtr               msg;

  for (vnp = results_list; vnp != NULL; vnp = vnp->next)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    if (res != NULL && res->bsp != NULL)
    {
      if (RemoveBarcodeTechFromBioseq (res->bsp))
      {
        if (fp != NULL)
        {
          msg = SummaryTextFromBarcodeTestResults (res);
          fprintf (fp, "%s\n", msg);
          msg = MemFree (msg);
        }
      }
    }
  }
}  


extern void RemoveBarcodeKeywords (FILE *fp, ValNodePtr results_list)
{
  BarcodeTestResultsPtr res;
  ValNodePtr            vnp;
  CharPtr               msg;

  for (vnp = results_list; vnp != NULL; vnp = vnp->next)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    if (res != NULL && res->bsp != NULL)
    {
      if (RemoveBarcodeKeywordFromBioseq (res->bsp))
      {
        if (fp != NULL)
        {
          msg = SummaryTextFromBarcodeTestResults (res);
          fprintf (fp, "%s\n", msg);
          msg = MemFree (msg);
        }
      }
    }
  }
}  



/* Applies Barcode technique to all Bioseqs in supplied list */
/* Used by Barcode Discrepancy Tool for the UNDO button.     */
extern void ApplyBarcodeKeywords (FILE *fp, ValNodePtr results_list)
{
  BarcodeTestResultsPtr res;
  ValNodePtr            vnp;
  CharPtr               msg;

  for (vnp = results_list; vnp != NULL; vnp = vnp->next)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    if (res != NULL && res->bsp != NULL)
    {
      ApplyBarcodeKeywordToBioseq (res->bsp);
      if (fp != NULL)
      {
        msg = SummaryTextFromBarcodeTestResults (res);
        fprintf (fp, "%s\n", msg);
        msg = MemFree (msg);
      }
    }
  }
}  


extern void ApplyBarcodeTech (FILE *fp, ValNodePtr results_list)
{
  BarcodeTestResultsPtr res;
  ValNodePtr            vnp;
  CharPtr               msg;

  for (vnp = results_list; vnp != NULL; vnp = vnp->next)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    if (res != NULL && res->bsp != NULL)
    {
      ApplyBarcodeTechToBioseq (res->bsp);
      if (fp != NULL)
      {
        msg = SummaryTextFromBarcodeTestResults (res);
        fprintf (fp, "%s\n", msg);
        msg = MemFree (msg);
      }
    }
  }
}  


#if defined (WIN32)
extern char * __stdcall AbstractReadFunction (Pointer userdata)
#else
extern char * AbstractReadFunction (Pointer userdata)
#endif
{
  ReadBufferPtr rbp;

  if (userdata == NULL) return NULL;

  rbp = (ReadBufferPtr) userdata;

  return MyFGetLine (rbp->fp, &(rbp->current_data));
}

#if defined (WIN32)
extern void __stdcall AbstractReportError (
#else
extern void AbstractReportError (
#endif
  TErrorInfoPtr err_ptr,
  Pointer      userdata
)
{
  TErrorInfoPtr PNTR list;
  TErrorInfoPtr last;

  if (err_ptr == NULL || userdata == NULL) return;

  list = (TErrorInfoPtr PNTR) userdata;

  if (*list == NULL)
  {
    *list = err_ptr;
  }
  else
  {
    for (last = *list; last != NULL && last->next != NULL; last = last->next)
    {}
    if (last != NULL) {
      last->next = err_ptr;
    }
  }

}



extern Boolean ParseLatLon (
  CharPtr lat_lon,
  FloatHi PNTR latP,
  FloatHi PNTR lonP
)

{
  char    ew;
  double  lat;
  double  lon;
  char    ns;

  if (latP != NULL) {
    *latP = 0.0;
  }
  if (lonP != NULL) {
    *lonP = 0.0;
  }

  if (StringHasNoText (lat_lon)) return FALSE;

  if (sscanf (lat_lon, "%lf %c %lf %c", &lat, &ns, &lon, &ew) == 4) {
    if (lon < -180.0) {
      lon = -180.0;
    }
    if (lat < -90.0) {
      lat = -90.0;
    }
    if (lon > 180.0) {
      lon = 180.0;
    }
    if (lat > 90.0) {
      lat = 90.0;
    }
    if (ew == 'W') {
      lon = -lon;
    }
    if (ns == 'S') {
      lat = -lat;
    }

    if (latP != NULL) {
      *latP = (FloatHi) lat;
    }
    if (lonP != NULL) {
      *lonP = (FloatHi) lon;
    }

    return TRUE;
  }

  return FALSE;
}

static CharPtr print_lat_lon_fmt = "%.*lf %c %.*lf %c";

static CharPtr MakeLatLonFromParts (FloatHi lat, Char ns, Int4 prec1, FloatHi lon, Char ew, Int4 prec2)
{
  Char buf [256];

  /* choose default directions when none supplied */
  if (ns == 0 && ew == 0)
  {
    ns = 'N';
    ew = 'E';
  }
  else if (ns == 0)
  {
    if (ew == 'E' || ew == 'W') 
    {
      ns = 'N';
    }
    else
    {
      ns = 'E';
    }
  }
  else if (ew == 0)
  {
    if (ns == 'N' || ns == 'S') 
    {
      ew = 'E';
    }
    else 
    {
      ew = 'N';
    }
  }

  /* correct -E to +W, -W to +W, -N to +S, -S to +S */
  if (lat < 0.0) 
  {
    if (ns == 'E') 
    {
      ns = 'W';
    }
    else if (ns == 'N')
    {
      ns = 'S';
    }
    lat = 0.0 - lat;
  }

  if (lon < 0.0) 
  {
    if (ew == 'E') 
    {
      ew = 'W';
    }
    else if (ew == 'N')
    {
      ew = 'S';
    }
    lon = 0.0 - lon;
  }

  if (ns == 'E' || ns == 'W')
  {
    sprintf (buf, print_lat_lon_fmt, prec2, lon, ew, prec1, lat, ns);
  }
  else
  {
    sprintf (buf, print_lat_lon_fmt, prec1, lat, ns, prec2, lon, ew);
  }
  return StringSave (buf);
}


static Int4 GetPrecisionFromNumberString (CharPtr str)
{
  CharPtr cp;
  Int4    prec = 0;

  if (StringHasNoText (str)) {
    return 0;
  }
  cp = str;

  while (isdigit (*cp)) {
    cp++;
  }
  if (*cp != '.') {
    return 0;
  }
  
  cp++;
  while (isdigit (*cp)) {
    prec++;
    cp++;
  }
  return prec;
}


extern void IsCorrectLatLonFormat (CharPtr lat_lon, BoolPtr format_correct, BoolPtr precision_correct, BoolPtr lat_in_range, BoolPtr lon_in_range)
{
  FloatHi  ns, ew;
  Char     lon, lat;
  Boolean  format_ok = FALSE, lat_ok = FALSE, lon_ok = FALSE, precision_okay = FALSE;
  Int4     processed, len, orig_len, ns_prec, ew_prec;
  CharPtr  buf, cp;

  if (StringHasNoText (lat_lon))
  {
    format_ok = FALSE;
  }
  else if (sscanf (lat_lon, "%lf %c %lf %c%n", &ns, &lat, &ew, &lon, &processed) != 4
           || processed != StringLen (lat_lon))
  {
    format_ok = FALSE;
  }
  else if ((lat != 'N' && lat != 'S') || (lon != 'E' && lon != 'W'))
  {
    format_ok = FALSE;
  }
  else 
  {
    cp = StringChr (lat_lon, ' ');
    if (cp != NULL) {
      cp = StringChr (cp + 1, ' ');
      if (cp != NULL) {
        cp++; 
      }
    }
    if (cp == NULL) {
      format_ok = FALSE;
    } else {
      ns_prec = GetPrecisionFromNumberString (lat_lon);
      ew_prec = GetPrecisionFromNumberString (cp);
      buf = MakeLatLonFromParts (ns, lat, ns_prec, ew, lon, ew_prec);
      len = StringLen (buf);
      orig_len = StringLen (lat_lon);
      if (StringNCmp (buf, lat_lon, len) == 0 &&
          (orig_len == len || (len < orig_len && lat_lon[len] == ';')))
      {
        format_ok = TRUE;
        if (ns <= 90 && ns >= 0)
        {
          lat_ok = TRUE;
        }
        if (ew <= 180 && ew >= 0)
        {
          lon_ok = TRUE;
        }
        if (ns_prec < 3 && ew_prec < 3) {
          precision_okay = TRUE;
        }
      }
      buf = MemFree (buf);
    }
  }

  if (format_correct != NULL)
  {
    *format_correct = format_ok;
  }
  if (precision_correct != NULL)
  {
    *precision_correct = precision_okay;
  }
  if (lat_in_range != NULL)
  {
    *lat_in_range = lat_ok;
  }
  if (lon_in_range != NULL)
  {
    *lon_in_range = lon_ok;
  }
}


static Boolean IsDirectionChar (Char dir)
{
  if (dir == 'E' || dir == 'W' || dir == 'N' || dir == 'S')
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static Boolean ParseNumericFromDToken (CharPtr dtoken, FloatHiPtr val, Int4Ptr prec)
{
  FloatLo  a, b, c;
  FloatHi  f = 0.0;
  Int4     i, j, k;
  Boolean  rval = FALSE;
  Int4     processed, len, dec_size;
  CharPtr  cp;

  if (StringHasNoText (dtoken) || val == NULL || prec == NULL)
  {
    return FALSE;
  }

  len = StringLen (dtoken);
  if ((sscanf (dtoken, "%d.%d.%d%n", &i, &j, &k, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%d.%d.%d'%n", &i, &j, &k, &processed) == 3 && processed == len))
  {
    f = (FloatHi) i + (FloatHi)j / (FloatHi)60.0 + (FloatHi)k / (FloatHi)3600.0;
    *prec = 4;
    rval = TRUE;
  }
  else if ((sscanf (dtoken, "%f:%f:%f%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f %f:%f%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f %f %f%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f %f'%f''%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f %f'%f\"%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f %f'%f%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f %f'%f'%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f'%f'%f%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f'%f'%f'%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f-%f-%f%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f %f-%f%n", &a, &b, &c, &processed) == 3 && processed == len))
  {
    f = a +  b / (FloatHi)60.0 + c / (FloatHi)3600.0;
    *prec = 4;
    rval = TRUE;
  }
  else if ((sscanf (dtoken, "%f %f%n", &a, &b, &processed) == 2 && processed == len)
           || (sscanf (dtoken, "%f:%f%n", &a, &b, &processed) == 2 && processed == len)
           || (sscanf (dtoken, "%f %f'%n", &a, &b, &processed) == 2 && processed == len)
           || (sscanf (dtoken, "%f'%f'%n", &a, &b, &processed) == 2 && processed == len)
           || (sscanf (dtoken, "%f'%f%n", &a, &b, &processed) == 2 && processed == len))
  {
    if (a < 0)
    {
      f = (FloatHi) a - b / (FloatHi) 60.0;
    }
    else
    {
      f = (FloatHi) a + b / (FloatHi) 60.0;
    }
    cp = StringRChr (dtoken, '.');
    if (cp == NULL)
    {
      dec_size = 0;
    }
    else
    {
      dec_size = StringLen (StringChr (dtoken, '.') + 1);
      if (dtoken[StringLen(dtoken) - 1] == '\'')
      {
        dec_size--;
      }
    }

    *prec = 2 + dec_size;
    rval = TRUE;
  }
  else if ((sscanf (dtoken, "%f%n", &a, &processed) == 1 && processed == len)
           || (sscanf (dtoken, "%f'%n", &a, &processed) == 1 && processed == len))
  {
    rval = TRUE;
    cp = StringChr (dtoken, '.');
    if (cp == NULL)
    {
      *prec = 2;
    }
    else
    {
      *prec = MAX (2, StringLen (cp + 1));
    }
    f = (FloatHi) a;
  }

  if (rval)
  {
    *val = f;
  } 
  return rval;
}

static Boolean ParseFromDToken (CharPtr dtoken, FloatHiPtr val, CharPtr d, Int4Ptr prec)
{
  FloatHi f;
  Char    dir = 0;
  Boolean rval = FALSE;
  Int4    token_len;

  if (StringHasNoText (dtoken) || val == NULL || d == NULL) 
  {
    return FALSE;
  }

  token_len = StringLen (dtoken);

  if (IsDirectionChar (dtoken[0])) 
  {
    dir = dtoken[0];
    rval = ParseNumericFromDToken (dtoken + 1, &f, prec);
  }
  else if (IsDirectionChar (dtoken[token_len - 1]))
  {
    dir = dtoken[token_len - 1];
    dtoken[token_len - 1] = 0;
    token_len--;
    while (token_len > 0 && isspace (dtoken[token_len - 1]))
    {
      dtoken[token_len - 1] = 0;
      token_len --; 
    }
    rval = ParseNumericFromDToken (dtoken, &f, prec);
  }
  else
  {
    rval = ParseNumericFromDToken (dtoken, &f, prec);
  }
  if (rval)
  {
    *val = f;
    *d = dir;
  }
  return rval;
}


static Boolean ParseFromLToken (CharPtr ltoken, Boolean first, FloatHiPtr val, CharPtr d, Int4Ptr prec)
{
  CharPtr  dtoken;
  Boolean  rval = FALSE;
  FloatHi  f;
  Char     dir;
  Char     plus_dir, minus_dir;
  Int4     len;

  if (StringHasNoText (ltoken) || val == NULL || d == NULL)
  {
    return rval;
  }
  len = StringLen (ltoken);
  if (StringNCmp (ltoken, "LAT", 3) == 0)
  {
    dtoken = ltoken + 3;
    plus_dir = 'N';
    minus_dir = 'S';
  } 
  else if (StringNCmp (ltoken, "LONG", 4) == 0)
  {
    dtoken = ltoken + 4;
    plus_dir = 'E';
    minus_dir = 'W';
  }
  else if (len > 3 && StringCmp (ltoken + len - 3, "LAT") == 0)
  {
    ltoken[len - 3] = 0;
    dtoken = ltoken;
    plus_dir = 'N';
    minus_dir = 'S';
  }
  else if (len > 4 && StringCmp (ltoken + len - 4, "LONG") == 0)
  {
    ltoken[len - 4] = 0;
    dtoken = ltoken;
    plus_dir = 'E';
    minus_dir = 'W';
  }
  else if (first)
  {
    dtoken = ltoken;
    plus_dir = 'N';
    minus_dir = 'S';
  }
  else
  {
    dtoken = ltoken;
    plus_dir = 'E';
    minus_dir = 'W';
  }
  /* trim space and punctuation from beginning */
  while (isspace (*dtoken) || (*dtoken != '-' && ispunct(*dtoken))) 
  {
    dtoken++;
  }
  /* trim space from end */
  len = StringLen (dtoken);
  while (len > 0 && isspace (dtoken[len - 1])) {
    dtoken[len - 1] = 0;
    len--;
  }
  if (ParseFromDToken (dtoken, &f, &dir, prec))
  {
    if (dir == 0)
    {
       if (f < 0)
       {
         dir = minus_dir;
         f = 0 - f;
       }
       else
       {
         dir = plus_dir;
       }
       rval = TRUE;
    }
    else if (dir == plus_dir || dir == minus_dir)
    {
      rval = TRUE;
    }
  }

  if (rval)
  {
    *val = f;
    *d = dir;
  }
  return rval;

}


static CharPtr MakeToken(CharPtr token1, CharPtr token2)
{
  Int4    token_len;
  CharPtr token;

  if (StringHasNoText (token1)) return NULL;
  while (isspace (*token1) || (ispunct (*token1) && *token1 != '-'))
  {
    token1++;
  }
  if (*token1 == 0)
  {
    return NULL;
  }
  if (token2 == NULL)
  {
    token_len = StringLen (token1) + 1;
  }
  else
  {
    token_len = token2 - token1 + 1;
  }
  token = (CharPtr) MemNew (sizeof (Char) * token_len);
  strncpy (token, token1, token_len - 1);
  token[token_len - 1] = 0;
  while ((isspace (token[token_len - 2]) || ispunct (token[token_len - 2])) && token_len > 2)
  {
    token[token_len - 2] = 0;
    token_len--;
  }
  return token;
}


static ReplacePairData latlon_replace_list[] = {
 { "LONGITUDE", "LONG" },
 { "LONG.",     "LONG" },
 { "LON.",      "LONG" },
 { "LATITUDE",  "LAT"  },
 { "LAT.",      "LAT"  },
 { "DEGREES",   " " },
 { "DEGREE",    " " },
 { "DEG.",      " " },
 { "DEG",       " " },
 { "MIN.",      "'" },
 { "MINUTES",    "'" },
 { "MINUTE",    "'" },
 { "MIN",       "'" },
 { "SEC.",      "''" },
 { "SEC",       "''" },
 { "NORTH",     "N"    },
 { "SOUTH",     "S"    },
 { "EAST",      "E"    },
 { "WEST",      "W"    },
};


static Int4 num_latlon_replace = sizeof (latlon_replace_list) / sizeof (ReplacePairData);


static Boolean CommaShouldBePeriod (CharPtr pStr, CharPtr pComma)
{
  CharPtr cp;
  Boolean rval = FALSE;

  if (StringHasNoText (pStr) || pComma == NULL || pComma <= pStr) return FALSE;

  cp = pComma - 1;
  if (!isdigit (*cp)) return FALSE;

  while (cp > pStr && isdigit (*cp)) 
  {
    cp--;
  }
  if (*cp != '.' && isdigit (*(pComma + 1)))
  {
    rval = TRUE;
  }
  return rval;
}

extern CharPtr FixLatLonFormat (CharPtr orig_lat_lon)
{
  FloatHi lon, lat;
  Char    ns, ew;
  CharPtr cpy, cp, dst;
  CharPtr first_dash = NULL, second_dash = NULL, first_space = NULL, second_space = NULL;
  CharPtr rval = NULL;
  CharPtr word1 = NULL, word2 = NULL;
  Boolean bad_letter_found = FALSE, replace_found, comma_sep = FALSE;
  Int4    i;
  CharPtr ltoken1 = NULL, ltoken2 = NULL;
  CharPtr dtoken1 = NULL, dtoken2 = NULL;
  Int4    find_len, replace_len, jump_len;
  Int4    prec1, prec2;
  CharPtr extra_text = NULL;
  CharPtr deg1 = NULL, deg2 = NULL;

  if (StringHasNoText (orig_lat_lon))
  {
    return NULL;
  }

  cpy = StringSave (orig_lat_lon);

  cp = cpy;
  while (*cp != 0)
  {
    if (*cp == 'O' 
        && (cp == cpy || !isalpha (*(cp - 1))))
    {
      *cp = '0';
    }
    else if (*cp == 'o' && cp != cpy && !isalpha (*(cp - 1)) && !isalpha (*(cp + 1)))
    {
      *cp = ' ';
    }
    else if (*cp == '#') /* # is sometimes used for degree, sometimes separator */
    {
      *cp = ' ';
    }
    else if (isalpha (*cp)) 
    {
      *cp = toupper (*cp);
    }
    else if (*cp == ',')
    {
      if (CommaShouldBePeriod (cpy, cp))
      {
        *cp = '.';
      }
    }
    cp++;
  }

  /* fix words */
  cp = cpy;
  dst = cpy;
  while (*cp != 0 && !bad_letter_found && extra_text == NULL)
  {
    if (isalpha (*cp)) 
    {
      replace_found = FALSE;
      for (i = 0; i < num_latlon_replace && !replace_found; i++) 
      {
        find_len = StringLen (latlon_replace_list[i].find);
        replace_len = StringLen (latlon_replace_list[i].replace);
        jump_len = 0;
        if (StringNICmp (cp, latlon_replace_list[i].find, find_len) == 0)
        {
          jump_len = find_len;
        }
        else if (StringNCmp (cp, latlon_replace_list[i].replace, replace_len) == 0)
        {
          jump_len = replace_len;
        }
        else
        {
          continue;
        }

        if (i < 5)
        {
          if (ltoken1 == NULL)
          {
            ltoken1 = dst;
          }
          else if (ltoken2 == NULL)
          {
            ltoken2 = dst;
          }
          else
          {
            bad_letter_found = TRUE;
          }
        }
        else if (i > 4 && i < 9)
        {
          if (deg1 == NULL) 
          {
            deg1 = dst;
          } 
          else if (deg2 == NULL) 
          {
            deg2 = dst;
          } 
          else 
          {
            bad_letter_found = TRUE;
          }
        }
        else if (i >= 15)
        {
          if (dtoken1 == NULL)
          {
            dtoken1 = dst;
          }
          else if (dtoken2 == NULL)
          {
            dtoken2 = dst;
          }
          else if (!comma_sep)
          {
            bad_letter_found = TRUE;
          }
        }

        if ((latlon_replace_list[i].replace[0] == '\'' || latlon_replace_list[i].replace[0] == ' ')
            && dst > cpy && isspace (*(dst - 1)))
        {
          /* no double spaces, no spaces before tick marks */ 
          dst--;
        }
        if (replace_len == 1 && latlon_replace_list[i].replace[0] == ' ')
        {
          if (isspace (*(cp + jump_len)) || *(cp + jump_len) == 0 || *(cp + jump_len) == '\'')
          {
            /* no double spaces */
          }
          else 
          {
            *dst = ' ';
            dst++;
          }
        }          
        else 
        {
          StringNCpy (dst, latlon_replace_list[i].replace, replace_len);
          dst += replace_len;
        }
        cp += jump_len;
        replace_found = TRUE;
      } 
      if (!replace_found)
      {
        bad_letter_found = 1;
      }
    }
    else if (isspace (*cp))
    {
      if (isspace (*(cp + 1)) || *(cp + 1) == 0 || *(cp + 1) == '\''
          || (dst > cpy && isspace (*(dst - 1))))
      {
        cp++;
      }
      else
      {
        *dst = ' ';
        if (first_space == NULL)
        {
          first_space = dst;
        }
        else if (second_space == NULL)
        {
          second_space = dst;
        }        
        dst++;
        cp++;
      }
    }
    else if (*cp == '-')
    {
      *dst = '-';
      if (first_dash == NULL)
      {
        first_dash = dst;
      }
      else if (second_dash == NULL)
      {
        second_dash = dst;
      }
      dst++;
      cp++;
    }
    else if (*cp == ',')
    {
      if (!comma_sep && ((ltoken1 != NULL && ltoken2 != NULL) || (dtoken1 != NULL && dtoken2 != NULL)))
      {
        extra_text = orig_lat_lon + (cp - cpy);
      } else {  
        *dst = ' ';
        dst++;
        cp++;
        if (dtoken1 != NULL && dtoken2 == NULL)
        {
          dtoken2 = dst;
          comma_sep = TRUE;
        }
        else if (ltoken1 == NULL && ltoken2 == NULL)
        {
          ltoken1 = cpy;
          ltoken2 = dst;
          comma_sep = TRUE;
        }
      }
    }
    else if (*cp == ';' && cp > cpy && isdigit(*(cp - 1)) && isdigit(*(cp + 1)))
    {
      /* replace typo semicolon with colon */
      *dst = ':';
      dst++;
      cp++;
    }
    else
    {
      *dst = *cp;
      dst++;
      cp++;
    }
  }

  *dst = 0;

  /* have to have both ltokens or none */
  if (ltoken1 != NULL && ltoken2 == NULL)
  {
    bad_letter_found = 1;
  }
  /* if no ltokens, must have both dtokens */
  else if (ltoken1 == NULL && (dtoken1 == NULL || dtoken2 == NULL))
  { 
    if (deg1 != NULL && deg2 != NULL) 
    {
      if (deg1 == cpy) {
        dtoken1 = deg1;
        dtoken2 = deg2;
      } else {
        cp = deg1;
        while (cp > cpy && isspace (*cp)) {
          cp--;
        }
        while (cp > cpy && !isspace (*cp)) {
          cp--;
        }
        if (isspace (*cp)) {
          cp++;
        }
        dtoken1 = cp;
        cp = deg2;
        while (cp > deg1 && isspace (*cp)) {
          cp--;
        }
        while (cp > deg1 && !isspace (*cp)) {
          cp--;
        }
        if (isspace (*cp)) {
          cp++;
        }
        dtoken2 = cp;
      }
    }
    /* use space to separate the two tokens */
    else if (first_space != NULL && second_space == NULL)
    {
      ltoken1 = cpy;
      ltoken2 = first_space + 1;
    }
    /* allow a dash to separate the two tokens if no spaces and only one dash */
    else if (first_space == NULL && second_space == NULL && first_dash != NULL && second_dash == NULL)
    {
      ltoken1 = cpy;
      *first_dash = ' ';
      ltoken2 = first_dash + 1;
    }
    else if (dtoken1 != NULL && dtoken2 == NULL && dtoken1 > cpy && dtoken1 < cpy + StringLen (cpy) - 1)
    {
      word1 = MakeToken (cpy, dtoken1 + 1);
      if (ParseFromDToken (word1, &lat, &ns, &prec1))
      {
        /* first portion parses ok, assume user just left off direction for second token */
        /* letters end tokens */
        dtoken2 = dtoken1 + 1;
        dtoken1 = cpy;
      }
      else
      {
        bad_letter_found = 1;
      }
      word1 = MemFree (word1);
    }
    else
    {
      bad_letter_found = 1;
    }
  }
  if (first_space == NULL && first_dash != NULL && second_dash == NULL && !comma_sep)
  {
    /* don't let the dash dividing the tokens be used as minus sign */
    *first_dash = ' ';
  }

  if (bad_letter_found)
  {
  }
  else if (ltoken1 != NULL)
  {
    /* if latitude and longitude are at end of token, change start */
    if (ltoken1 != cpy)
    {
      ltoken2 = ltoken1 + 3;
      if (*ltoken2 == 'G') 
      {
        ltoken2++;
      }
      ltoken1 = cpy;
    }
    word1 = MakeToken(ltoken1, ltoken2);
    word2 = MakeToken(ltoken2, NULL);
    if (ParseFromLToken (word1, TRUE, &lat, &ns, &prec1)
        && ParseFromLToken (word2, FALSE, &lon, &ew, &prec2))
    {
      if (prec1 > 2) {
        prec1 = 2;
      }
      if (prec2 > 2) {
        prec2 = 2;
      }
      rval = MakeLatLonFromParts (lat, ns, prec1, lon, ew, prec2);
    }
  }
  else
  {
    if (dtoken1 != cpy) 
    {
      /* letters end tokens */
      dtoken2 = dtoken1 + 1;
      dtoken1 = cpy;
    }
    word1 = MakeToken (dtoken1, dtoken2);
    word2 = MakeToken (dtoken2, NULL);
    if (ParseFromDToken (word1, &lat, &ns, &prec1)
        && ParseFromDToken (word2, &lon, &ew, &prec2))
    {
      if (prec1 > 2) {
        prec1 = 2;
      }
      if (prec2 > 2) {
        prec2 = 2;
      }
      rval = MakeLatLonFromParts (lat, ns, prec1, lon, ew, prec2);
    }
  }
      
  word1 = MemFree (word1);
  word2 = MemFree (word2);
  cpy = MemFree (cpy);
  
  if (rval != NULL && extra_text != NULL)
  {
    cpy = (CharPtr) MemNew (sizeof (Char) * (StringLen (rval) + StringLen (extra_text) + 1));
    sprintf (cpy, "%s%s", rval, extra_text);
    rval = MemFree (rval);
    rval = cpy;
  }
  return rval;
}


static void TestLatLonFormatting (FILE *fp)
{
  CharPtr tests[]  = 
  { "100.12 N 200.12 E",     /* already correct */
    "100 N 200 E",           /* correctable */
    "100.1 N 200.2 E",       /* correctable */
    "1OO.1 N 200.2 E",       /* correctable (replace capital o with zero) */
    "100.1 N, 200.2 E",      /* correctable (remove comma) */
    "E 100, S 120",          /* correctable (remove comma, reverse order, letters before numbers */
    "latitude: 200 N longitude: 100 E",
    "latitude: 200 E longitude: 100 N", /* NOT correctable */
    "N 37 45.403', 119 1.456' W",
    "38 52 56 N 84 44 53 W",
    "49 29 50 N 80 25 52 W",
    "39N 93W",
    "42:43:13N 01:0015W",
    "02deg 33min 00.7sec S 45deg 01min 38.8sec W",
    "42:24:37.9 N 85:22:11.7 W",
    "10 N 124 E",
    "41deg30'' S 145deg37' E",
    "59.30deg N 22.40deg E",
    "35 N 134 E",
    "2 S 114 E",
    "24deg 24.377' N 101deg 23.073' W'",
    "26deg 57.9' N 102deg 08.3 W'",
    "38 11 44.66 North 0 35 01.93 West",
    "62.08 N 129.682",
    "64.444 N -164.973",
    "62.033 N -146.533",
    "67 N -51",
    "69.107 N 124.195",
    "2:46:00-59:41:00",
    "64 degree 55 N 25 degree 05 E",
    "64.907 N -166.18",
    "2:46:00-59:41:00",
    "66 degree 21 N 29 degree 21 E",
    "37deg27N 121deg52'W",
    "01deg31'25''N 66''33'31''W",
    "07deg33'30''N 69deg20'W",
    "10.8439,-85.6138",
    "11.03,-85.527",
    "8 deg 45 min S, 63 deg 26 min W",
    "29deg 49' 23.7' N; 106deg 23' 15.8'W",
    "7:46S, 12:30E",
    "35deg48'50'' N; 82deg5658'' W",
    "45deg34.18''N, 122deg12.00 'W",
    "37deg27N, 121deg52'W",
    "41:00;00N 20:45:00E",
    "02 deg 28' 29# S, 56 deg 6' 31# W"
};
  Int4 test_num, num_tests = sizeof (tests) / sizeof (char *);
  CharPtr fix;
  Int4 num_pass = 0, num_formatted = 0;
  Boolean format_ok, lat_in_range, lon_in_range, precision_ok;

  if (fp == NULL) return;

  for (test_num = 0; test_num < num_tests; test_num++)
  {
    fprintf (fp, "Test %d: %s\n", test_num, tests[test_num]);
    fix = FixLatLonFormat (tests[test_num]);
    if (fix == NULL) 
    {
      fprintf (fp, "Unable to correct format\n");
    }
    else
    {
      IsCorrectLatLonFormat (fix, &format_ok, &precision_ok, &lat_in_range, &lon_in_range);
      if (format_ok && precision_ok)
      {
        num_formatted ++;
        fprintf (fp, "Correction succeeded:%s\n", fix);
        num_pass++;
      }
      else
      {
        num_formatted ++;
        fprintf (fp, "Correction failed:%s\n", fix);
      }
    }
  }
  fprintf (fp, "Formats %d out of %d, %d succeed\n", num_formatted, num_tests, num_pass);
}


static CharPtr StringFromObjectID (ObjectIdPtr oip)
{
  CharPtr    str;
  if (oip == NULL) return NULL;

  if (oip->id > 0)
  {
    str = (CharPtr) MemNew (sizeof (Char) * 20);
    sprintf (str, "%d", oip->id);
  }
  else
  {
    str = StringSave (oip->str);
  }
  return str;
}

static Boolean ApplyBarcodeDbxrefToBioSource (BioSourcePtr biop, ObjectIdPtr oip)
{
  ValNodePtr vnp;
  DbtagPtr   dbt;
  CharPtr    str, cmp;
  Boolean    found = FALSE;
  Boolean    rval = FALSE;

  if (biop == NULL || oip == NULL) return FALSE;

  if (biop->org == NULL)
  {
    biop->org = OrgRefNew();
  }

  str = StringFromObjectID (oip);

  for (vnp = biop->org->db; vnp != NULL && !found; vnp = vnp->next)
  {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt == NULL || dbt->tag == NULL) continue;
    if (StringCmp (dbt->db, "BOLD") != 0) continue;
    cmp = StringFromObjectID (dbt->tag);
    if (StringCmp (str, cmp) == 0) found = TRUE;
    cmp = MemFree (cmp);
  }
  if (found) 
  {
    str = MemFree (str);
  }
  else
  {
    dbt = DbtagNew ();
    dbt->db = StringSave ("BOLD");
    dbt->tag = ObjectIdNew();
    dbt->tag->str = str;
    ValNodeAddPointer (&(biop->org->db), 0, dbt);
    rval = TRUE;
  }
  return rval;
}


extern void ApplyBarcodeDbxrefsToBioseq (BioseqPtr bsp, Pointer data)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  SeqIdPtr          sip;
  DbtagPtr          dbt;
  Int4Ptr           p_num;

  if (bsp == NULL) return;
  for (sip = bsp->id; sip != NULL; sip = sip->next)
  {
    if (IsBarcodeID (sip) && sip->choice == SEQID_GENERAL && sip->data.ptrvalue != NULL) 
    {
      dbt = (DbtagPtr) sip->data.ptrvalue;
      
      sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
      if (sdp != NULL)
      {
        if (ApplyBarcodeDbxrefToBioSource ((BioSourcePtr) sdp->data.ptrvalue, dbt->tag)) {
          if ((p_num = (Int4Ptr) data) != NULL) {
            (*p_num)++;
          }
        }
      }
    }
  }
}


/* Code for Country Fixup */

static Boolean IsSubstringOfStringInList (CharPtr whole_str, CharPtr match_p, CharPtr match_str, CharPtr PNTR list)
{
  CharPtr cp;
  Int4    context_len, find_len;
  Boolean rval = FALSE;

  if (list == NULL || StringHasNoText (whole_str) || match_p == NULL || match_p < whole_str) {
    return FALSE;
  }
  find_len = StringLen (match_str);
  while (*list != NULL && !rval) {
    context_len = StringLen (*list);
    if (find_len < context_len) {
      cp = StringSearch (whole_str, *list);
      while (cp != NULL && !rval) {
        if (match_p < cp) {
          cp = NULL;
        } else if (cp + context_len > match_p) {
          rval = TRUE;
        } else {
          cp = StringSearch (cp + 1, *list);
        }
      }
    }
    list++;
  }
  return rval;
}


static ReplacePairData country_name_fixes[] = {
 {"Vietnam", "Viet Nam"},
 {"Ivory Coast", "Cote d'Ivoire"},
 {"United States of America", "USA"},
 {"U.S.A.", "USA"},
 {"The Netherlands", "Netherlands"},
 {"People's Republic of China", "China"},
 {"Pr China", "China" },
 {"Prchina", "China" },
 {"P.R.China", "China" },
 {"P.R. China", "China" },
 {"P, R, China", "China" },
 {"Canary Islands", "Spain: Canary Islands"},
 {"Tenerife", "Spain: Tenerife"},
 {"Gran Canaria", "Spain: Gran Canaria"},
 {"Fuerteventura", "Spain: Fuerteventura"},
 {"Lanzarote", "Spain: Lanzarote"},
 {"La Palma", "Spain: La Palma"},
 {"La Gomera", "Spain: La Gomera"},
 {"El Hierro", "Spain: El Hierro"},
 {"La Graciosa", "Spain: La Graciosa"},
 {"Madeira", "Portugal: Madeira"},
 {"Azores", "Portugal: Azores"},
 {"Autonomous Region of the Azores", "Portugal: Azores"},
 {"St. Lucia", "Saint Lucia"},
 {"St Lucia", "Saint Lucia"},
 {"St. Thomas", "USA: Saint Thomas"},
 {"St Thomas", "USA: Saint Thomas"},
 {"Saint Kitts & Nevis", "Saint Kitts and Nevis"},
 {"Saint Kitts", "Saint Kitts and Nevis: Saint Kitts"},
 {"St. Kitts", "Saint Kitts and Nevis: Saint Kitts"},
 {"St Kitts", "Saint Kitts and Nevis: Saint Kitts"},
 {"Nevis", "Saint Kitts and Nevis: Nevis"},
 {"St. Helena", "Saint Helena"},
 {"St Helena", "Saint Helena"},
 {"Saint Pierre & Miquelon", "Saint Pierre and Miquelon"},
 {"St. Pierre", "Saint Pierre and Miquelon: Saint Pierre"},
 {"St Pierre", "Saint Pierre and Miquelon: Saint Pierre"},
 {"Saint Pierre", "Saint Pierre and Miquelon: Saint Pierre"},
 {"Miquelon", "Saint Pierre and Miquelon: Miquelon"},
 {"St. Pierre and Miquelon", "Saint Pierre and Miquelon"},
 {"St Pierre and Miquelon", "Saint Pierre and Miquelon"},
 {"Saint Vincent & the Grenadines", "Saint Vincent and the Grenadines"},
 {"Saint Vincent & Grenadines", "Saint Vincent and the Grenadines"},
 {"Saint Vincent and Grenadines", "Saint Vincent and the Grenadines"},
 {"St. Vincent and the Grenadines", "Saint Vincent and the Grenadines"},
 {"St Vincent and the Grenadines", "Saint Vincent and the Grenadines"},
 {"Grenadines", "Saint Vincent and the Grenadines: Grenadines"},
 {"St Vincent", "Saint Vincent and the Grenadines: Saint Vincent"},
 {"St. Vincent", "Saint Vincent and the Grenadines: Saint Vincent"},
 {"Saint Vincent", "Saint Vincent and the Grenadines: Saint Vincent"},
 {"Cape Verde Islands", "Cape Verde"},
 {"Trinidad & Tobago", "Trinidad and Tobago"},
 {"Trinidad", "Trinidad and Tobago: Trinidad"},
 {"Tobago", "Trinidad and Tobago: Tobago"},
 {"Ashmore & Cartier Islands", "Ashmore and Cartier Islands"},
 {"Ashmore Island", "Ashmore and Cartier Islands: Ashmore Island"},
 {"Cartier Island", "Ashmore and Cartier Islands: Cartier Island"},
 {"Heard Island & McDonald Islands", "Heard Island and McDonald Islands"},
 {"Heard Island", "Heard Island and McDonald Islands: Heard Island"},
 {"McDonald Islands", "Heard Island and McDonald Islands: McDonald Islands"},
 {"McDonald Island", "Heard Island and McDonald Islands: McDonald Island"},
 {"Sao Tome & Principe", "Sao Tome and Principe"},
 {"Principe", "Sao Tome and Principe: Principe"},
 {"Sao Tome", "Sao Tome and Principe: Sao Tome"},
 {"South Sandwich Islands", "South Georgia and the South Sandwich Islands: South Sandwich Islands"},
 {"Turks & Caicos", "Turks and Caicos Islands"},
 {"Turks & Caicos Islands", "Turks and Caicos Islands"},
 {"Turks and Caicos", "Turks and Caicos Islands"},
 {"Turks Islands", "Turks and Caicos Islands: Turks Islands"},
 {"Caicos Islands", "Turks and Caicos Islands: Caicos Islands"},
 {"Antigua & Barbuda", "Antigua and Barbuda"},
 {"Antigua", "Antigua and Barbuda: Antigua"},
 {"Barbuda", "Antigua and Barbuda: Barbuda"},
 {"Falkland Islands", "Falkland Islands (Islas Malvinas)"},
 {"French Southern & Antarctic Lands", "French Southern and Antarctic Lands"},
 {"Ile Amsterdam", "French Southern and Antarctic Lands: Ile Amsterdam"},
 {"Ile Saint-Paul", "French Southern and Antarctic Lands: Ile Saint-Paul"},
 {"Iles Crozet", "French Southern and Antarctic Lands: Iles Crozet"},
 {"Iles Kerguelen", "French Southern and Antarctic Lands: Iles Kerguelen"},
 {"Bassas da India", "French Southern and Antarctic Lands: Bassas da India"},
 {"Europa Island", "French Southern and Antarctic Lands: Europa Island"},
 {"Glorioso Islands", "French Southern and Antarctic Lands: Glorioso Islands"},
 {"Juan de Nova Island", "French Southern and Antarctic Lands: Juan de Nova Island"},
 {"Tromelin Island", "French Southern and Antarctic Lands: Tromelin Island"},
 {"South Georgia & the South Sandwich Islands", "South Georgia and the South Sandwich Islands"},
 {"South Georgia & South Sandwich Islands", "South Georgia and the South Sandwich Islands"},
{"ABW", "Aruba"},
{"AFG", "Afghanistan"},
{"AGO", "Angola"},
{"AIA", "Anguilla"},
{"ALA", "Aland Islands"},
{"ALB", "Albania"},
{"AND", "Andorra"},
{"ARE", "United Arab Emirates"},
{"ARG", "Argentina"},
{"ARM", "Armenia"},
{"ASM", "American Samoa"},
{"ATA", "Antarctica"},
{"ATF", "French Southern Territories"},
{"ATG", "Antigua and Barbuda"},
{"AUS", "Australia"},
{"AUT", "Austria"},
{"AZE", "Azerbaijan"},
{"BDI", "Burundi"},
{"BEL", "Belgium"},
{"BEN", "Benin"},
{"BES", "Bonaire, Sint Eustatius and Saba"},
{"BFA", "Burkina Faso"},
{"BGD", "Bangladesh"},
{"BGR", "Bulgaria"},
{"BHR", "Bahrain"},
{"BHS", "Bahamas"},
{"BIH", "Bosnia and Herzegovina"},
{"BLM", "Saint Barthelemy"},
{"BLR", "Belarus"},
{"BLZ", "Belize"},
{"BMU", "Bermuda"},
{"BOL", "Bolivia"},
{"BRA", "Brazil"},
{"BRB", "Barbados"},
{"BRN", "Brunei"},
{"BTN", "Bhutan"},
{"BVT", "Bouvet Island"},
{"BWA", "Botswana"},
{"CAF", "Central African Republic"},
{"CAN", "Canada"},
{"CCK", "Cocos Islands"},
{"CHE", "Switzerland"},
{"CHL", "Chile"},
{"CHN", "China"},
{"CIV", "Cote d'Ivoire"},
{"CMR", "Cameroon"},
{"COD", "Democratic Republic of the Congo"},
{"COG", "Republic of the Congo"},
{"COK", "Cook Islands"},
{"COL", "Colombia"},
{"COM", "Comoros"},
{"CPV", "Cape Verde"},
{"CRI", "Costa Rica"},
{"CUB", "Cuba"},
{"CUW", "Curacao"},
{"CXR", "Christmas Island"},
{"CYM", "Cayman Islands"},
{"CYP", "Cyprus"},
{"CZE", "Czech Republic"},
{"DEU", "Germany"},
{"DJI", "Djibouti"},
{"DMA", "Dominica"},
{"DNK", "Denmark"},
{"DOM", "Dominican Republic"},
{"DZA", "Algeria"},
{"ECU", "Ecuador"},
{"EGY", "Egypt"},
{"ERI", "Eritrea"},
{"ESH", "Western Sahara"},
{"ESP", "Spain"},
{"EST", "Estonia"},
{"ETH", "Ethiopia"},
{"FIN", "Finland"},
{"FJI", "Fiji"},
{"FLK", "Falkland Islands (Islas Malvinas)"},
{"FRA", "France"},
{"FRO", "Faroe Islands"},
{"FSM", "Micronesia"},
{"GAB", "Gabon"},
{"GBR", "United Kingdom"},
{"GEO", "Georgia"},
{"GGY", "Guernsey"},
{"GHA", "Ghana"},
{"GIB", "Gibraltar"},
{"GIN", "Guinea"},
{"GLP", "Guadeloupe"},
{"GMB", "Gambia"},
{"GNB", "Guinea-Bissau"},
{"GNQ", "Equatorial Guinea"},
{"GRC", "Greece"},
{"GRD", "Grenada"},
{"GRL", "Greenland"},
{"GTM", "Guatemala"},
{"GUF", "French Guiana"},
{"GUM", "Guam"},
{"GUY", "Guyana"},
{"HKG", "Hong Kong"},
{"HMD", "Heard Island and McDonald Islands"},
{"HND", "Honduras"},
{"HRV", "Croatia"},
{"HTI", "Haiti"},
{"HUN", "Hungary"},
{"IDN", "Indonesia"},
{"IMN", "Isle of Man"},
{"IND", "India"},
{"IOT", "British Indian Ocean Territory"},
{"IRL", "Ireland"},
{"IRN", "Iran"},
{"IRQ", "Iraq"},
{"ISL", "Iceland"},
{"ISR", "Israel"},
{"ITA", "Italy"},
{"JAM", "Jamaica"},
{"JEY", "Jersey"},
{"JOR", "Jordan"},
{"JPN", "Japan"},
{"KAZ", "Kazakhstan"},
{"KEN", "Kenya"},
{"KGZ", "Kyrgyzstan"},
{"KHM", "Cambodia"},
{"KIR", "Kiribati"},
{"KNA", "Saint Kitts and Nevis"},
{"KOR", "South Korea"},
{"KWT", "Kuwait"},
{"LAO", "Lao People's Democratic Republic"},
{"LBN", "Lebanon"},
{"LBR", "Liberia"},
{"LBY", "Libyan Arab Jamahiriya"},
{"LCA", "Saint Lucia"},
{"LIE", "Liechtenstein"},
{"LKA", "Sri Lanka"},
{"LSO", "Lesotho"},
{"LTU", "Lithuania"},
{"LUX", "Luxembourg"},
{"LVA", "Latvia"},
{"MAC", "Macao"},
{"MAF", "Saint Martin (French part)"},
{"MAR", "Morocco"},
{"MCO", "Monaco"},
{"MDA", "Moldova"},
{"MDG", "Madagascar"},
{"MDV", "Maldives"},
{"MEX", "Mexico"},
{"MHL", "Marshall Islands"},
{"MKD", "Macedonia"},
{"MLI", "Mali"},
{"MLT", "Malta"},
{"MMR", "Myanmar"},
{"MNE", "Montenegro"},
{"MNG", "Mongolia"},
{"MNP", "Northern Mariana Islands"},
{"MOZ", "Mozambique"},
{"MRT", "Mauritania"},
{"MSR", "Montserrat"},
{"MTQ", "Martinique"},
{"MUS", "Mauritius"},
{"MWI", "Malawi"},
{"MYS", "Malaysia"},
{"MYT", "Mayotte"},
{"NAM", "Namibia"},
{"NCL", "New Caledonia"},
{"NER", "Niger"},
{"NFK", "Norfolk Island"},
{"NGA", "Nigeria"},
{"NIC", "Nicaragua"},
{"NIU", "Niue"},
{"NLD", "Netherlands"},
{"NOR", "Norway"},
{"NPL", "Nepal"},
{"NRU", "Nauru"},
{"NZL", "New Zealand"},
{"OMN", "Oman"},
{"PAK", "Pakistan"},
{"PAN", "Panama"},
{"PCN", "Pitcairn"},
{"PER", "Peru"},
{"PHL", "Philippines"},
{"PLW", "Palau"},
{"PNG", "Papua New Guinea"},
{"POL", "Poland"},
{"PRI", "Puerto Rico"},
{"PRK", "North Korea"},
{"PRT", "Portugal"},
{"PRY", "Paraguay"},
{"PSE", "Palestinian Territory"},
{"PYF", "French Polynesia"},
{"QAT", "Qatar"},
{"REU", "Reunion"},
{"ROU", "Romania"},
{"RUS", "Russia"},
{"RWA", "Rwanda"},
{"SAU", "Saudi Arabia"},
{"SDN", "Sudan"},
{"SEN", "Senegal"},
{"SGP", "Singapore"},
{"SGS", "South Georgia and the South Sandwich Islands"},
{"SHN", "Saint Helena"},
{"SJM", "Svalbard and Jan Mayen"},
{"SLB", "Solomon Islands"},
{"SLE", "Sierra Leone"},
{"SLV", "El Salvador"},
{"SMR", "San Marino"},
{"SOM", "Somalia"},
{"SPM", "Saint Pierre and Miquelon"},
{"SRB", "Serbia"},
{"SSD", "South Sudan"},
{"STP", "Sao Tome and Principe"},
{"SUR", "Suriname"},
{"SVK", "Slovakia"},
{"SVN", "Slovenia"},
{"SWE", "Sweden"},
{"SWZ", "Swaziland"},
{"SXM", "Sint Maarten (Dutch part)"},
{"SYC", "Seychelles"},
{"SYR", "Syrian Arab Republic"},
{"TCA", "Turks and Caicos Islands"},
{"TCD", "Chad"},
{"TGO", "Togo"},
{"THA", "Thailand"},
{"TJK", "Tajikistan"},
{"TKL", "Tokelau"},
{"TKM", "Turkmenistan"},
{"TLS", "Timor-Leste"},
{"TON", "Tonga"},
{"TTO", "Trinidad and Tobago"},
{"TUN", "Tunisia"},
{"TUR", "Turkey"},
{"TUV", "Tuvalu"},
{"TWN", "Taiwan"},
{"TZA", "Tanzania"},
{"UGA", "Uganda"},
{"UKR", "Ukraine"},
{"UMI", "United States Minor Outlying Islands"},
{"URY", "Uruguay"},
{"USA", "United States"},
{"UZB", "Uzbekistan"},
{"VAT", "Holy See (Vatican City State)"},
{"VCT", "Saint Vincent and the Grenadines"},
{"VEN", "Venezuela"},
{"VGB", "British Virgin Islands"},
{"VIR", "Virgin Islands"},
{"VNM", "Viet Nam"},
{"VUT", "Vanuatu"},
{"WLF", "Wallis and Futuna"},
{"WSM", "Samoa"},
{"YEM", "Yemen"},
{"ZAF", "South Africa"},
{"ZMB", "Zambia"},
{"ZWE", "Zimbabwe"},
 {NULL, NULL}
};

NLM_EXTERN CharPtr GetStateAbbreviation (CharPtr state)
{
  ReplacePairPtr fix;
  CharPtr        abbrev = NULL;

  fix = us_state_abbrev_fixes;
  while (fix->find != NULL && abbrev == NULL) {
    if (StringICmp (fix->replace, state) == 0) {
      abbrev = fix->find;
    } 
    fix++;
  }
  return abbrev;
}


static Boolean ContainsMultipleCountryNames (CharPtr PNTR list, CharPtr search_str)
{
  CharPtr PNTR  ptr;
  Int4          len_match;
  CharPtr       cp;
  Boolean       found_one = FALSE;
  
  if (list == NULL || search_str == NULL) return FALSE;
  
  for (ptr = list; ptr != NULL && *ptr != NULL; ptr++)
  {
    cp = StringISearch (search_str, *ptr);
    len_match = StringLen (*ptr);
    while (cp != NULL) {
      /* if character after match is alpha, continue */
      if (isalpha ((Int4)(cp [len_match]))
          /* if character before match is alpha, continue */
          || (cp > search_str && isalpha ((Int4)(*(cp - 1))))
        /* if is shorter match for other item, continue */
        || IsSubstringOfStringInList (search_str, cp, *ptr, list)) {
        cp = StringSearch (cp + len_match, *ptr);
      } else if (found_one) {
        return TRUE;
      } else {
        found_one = TRUE;
        cp = StringSearch (cp + len_match, *ptr);
      }
    }
  }
  return FALSE;
}


static CharPtr NewFixCountry (CharPtr country, CharPtr PNTR country_list)
{
  CharPtr cp, next_sep, start_after;
  CharPtr valid_country = NULL, new_country = NULL, tmp;
  Char    ch;
  CharPtr separator_list = ",:";
  Boolean too_many_countries = FALSE, bad_cap = FALSE;
  Int4    len_country, len_before, len_after, len_diff;
  ReplacePairPtr fix;
  Boolean fix_found;

  country = StringSave (country);
  cp = country;
  while (*cp != 0 && !too_many_countries) {
    next_sep = cp + StringCSpn (cp, separator_list);
    ch = *next_sep;
    *next_sep = 0;
    
    if (CountryIsValid (cp, NULL, &bad_cap)) {
      if (valid_country == NULL) {
        valid_country = cp;
      } else {
        too_many_countries = TRUE;
      }
    } else {
      /* see if this is a fixable country */
      fix = country_name_fixes;
      fix_found = FALSE;
      while (fix->find != NULL && !fix_found) {
        if (StringCmp (fix->find, cp) == 0) {
          fix_found = TRUE;
          if (valid_country == NULL) {
            len_before = cp - country;
            if (ch == 0) {
              len_after = 0;
            } else {
              len_after = StringLen (next_sep + 1) + 1;
            }
            len_diff = StringLen (fix->replace) - StringLen (fix->find);
            len_country = StringLen (country) + len_diff + len_after + 1;       
            tmp = (CharPtr) MemNew (sizeof (Char) * len_country);
            if (len_before > 0) {
              StringNCpy (tmp, country, len_before);
            }
            StringCpy (tmp + len_before, fix->replace);
            if (len_after > 0) {
              StringCpy (tmp + len_before + StringLen (fix->replace) + 1, next_sep + 1);
            }
            cp = tmp + len_before;
            valid_country = cp;
            next_sep = tmp + (next_sep - country) + len_diff;
            country = MemFree (country);
            country = tmp;
          } else {
            too_many_countries = TRUE;
          }
        }
        fix++;
      }
    }

    *next_sep = ch;
    if (*next_sep == 0) {
      cp = next_sep;
    } else {
      cp = next_sep + 1;
      while (isspace (*cp)) {
        cp++;
      }
    }
  }
  if (valid_country != NULL && !too_many_countries) {
    too_many_countries = ContainsMultipleCountryNames (country_list, country);
  }

  if (valid_country != NULL && too_many_countries && valid_country == country) {
    len_country = StringCSpn (valid_country, separator_list);
    if (country[len_country] == ':' && !isspace (country[len_country + 1])) {
      new_country = MemNew (sizeof (Char) * (StringLen (country) + 2));
      StringNCpy (new_country, country, len_country + 1);
      StringCat (new_country, " ");
      StringCat (new_country, country + len_country + 1);
    }   
  } else if (valid_country != NULL && !too_many_countries) {
    len_country = StringCSpn (valid_country, separator_list);
    len_before = valid_country - country;

    while (len_before > 0 
           && (isspace (country [len_before - 1]) 
               || StringChr (separator_list, country [len_before - 1]) != NULL)) {
      len_before--;
    }
    start_after = valid_country + len_country;
    while (*start_after != 0 
           && (isspace (*start_after)
               || StringChr (separator_list, *start_after) != NULL)) {
      start_after++;
    }

    len_after = StringLen (start_after);

    new_country = MemNew (sizeof (Char) * (len_country + len_before + len_after + 5));
    
    if (bad_cap && valid_country != NULL) {
      next_sep = valid_country + StringCSpn (valid_country, separator_list);
      ch = *next_sep;
      *next_sep = 0;
      tmp = GetCorrectedCountryCapitalization(valid_country);
      *next_sep = ch;
      if (tmp == NULL) {
        StringNCpy (new_country, valid_country, len_country);
      } else {
        StringNCpy (new_country, tmp, len_country);
      }
    } else {
      StringNCpy (new_country, valid_country, len_country);
    }
    if (len_before > 0 || len_after > 0) {
      StringCat (new_country, ": ");
      if (len_before > 0) {
        StringNCat (new_country, country, len_before);
        if (len_after > 0) {
          StringCat (new_country, ", ");
        }
      }
      if (len_after > 0) {
        StringCat (new_country, start_after);
      }
    }
  }
  country = MemFree (country);
  return new_country;
}


extern CharPtr GetCountryFix (CharPtr country, CharPtr PNTR country_list)
{
  CharPtr new_country;

  if (StringHasNoText (country)) return NULL;
  new_country = NewFixCountry (country, country_list);
  return new_country;
}


typedef struct countryfixup {
  CharPtr PNTR country_list;
  ValNodePtr warning_list;
  Boolean capitalize_after_colon;
  Boolean any_changed;
  FILE *log_fp;
} CountryFixupData, PNTR CountryFixupPtr;


static void CapitalizeFirstLetterOfEveryWord (CharPtr pString)
{
  CharPtr pCh;

  pCh = pString;
  if (pCh == NULL) return;
  if (*pCh == '\0') return;
  
  while (*pCh != 0)
  {
    /* skip over spaces */
    while (isspace(*pCh))
    {
      pCh++;
    }
  
    /* capitalize first letter after white space */
    if (isalpha (*pCh))
    {
      *pCh = toupper (*pCh);
      pCh++;
    }
    /* skip over rest of word */
    while (*pCh != 0 && !isspace (*pCh))
    {
      if (isalpha (*pCh)) {
        *pCh = tolower (*pCh);
      }
      pCh++;
    }
  }
}


static void CountryFixupItem (Uint1 choice, Pointer data, CountryFixupPtr c)
{
  BioSourcePtr biop;
  SubSourcePtr ssp;
  CharPtr      new_country;
  CharPtr      cp;
  CharPtr      tmp;
  Int4         country_len;

  if (data == NULL || c == NULL) return;

  biop = GetBioSourceFromObject (choice, data);
  if (biop == NULL) return;

  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) 
  {
  	if (ssp->subtype == SUBSRC_country && !StringHasNoText (ssp->name))
    {
      new_country = GetCountryFix (ssp->name, c->country_list);
      if (new_country == NULL) {
        ValNodeAddPointer (&c->warning_list, choice, data);
      } else {
        cp = StringChr (new_country, ':');
	      if (cp != NULL) {
          country_len = cp - new_country;
	        /* skip colon */
 	        cp++;
	        /* skip over space after colon */
	        cp += StringSpn (cp, " \t");
          if (c->capitalize_after_colon) {       	  
  	        /* reset capitalization */
  	        CapitalizeFirstLetterOfEveryWord (cp);
          }
          if (*(new_country + country_len + 1) != 0 && !isspace (*(new_country + country_len + 1))) {
            tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (new_country) + 2));
            StringNCpy (tmp, new_country, country_len + 1);
            StringCat (tmp, " ");
            StringCat (tmp, cp + 1);
            new_country = MemFree (new_country);
            new_country = tmp;
          }
        }
        if (StringCmp (ssp->name, new_country) == 0) {
          new_country = MemFree (new_country);
        } else {
          c->any_changed = TRUE;
          if (c->log_fp != NULL) {
            fprintf (c->log_fp, "Changed '%s' to '%s'\n", ssp->name, new_country);
          }
          ssp->name = MemFree (ssp->name);
          ssp->name = new_country;
        }
      }
    }
  }
}


static void CountryFixupDesc (SeqDescrPtr sdp, Pointer userdata)
{
  if (sdp != NULL && userdata != NULL && sdp->choice == Seq_descr_source) {
    CountryFixupItem (OBJ_SEQDESC, sdp, (CountryFixupPtr) userdata);
  }
}


static void CountryFixupFeat (SeqFeatPtr sfp, Pointer userdata)
{
  if (sfp != NULL && userdata != NULL && sfp->data.choice == SEQFEAT_BIOSRC) {
    CountryFixupItem (OBJ_SEQFEAT, sfp, (CountryFixupPtr) userdata);
  }
}


NLM_EXTERN ValNodePtr FixupCountryQuals (SeqEntryPtr sep, Boolean fix_after_colon)
{
  CountryFixupData c;

  MemSet (&c, 0, sizeof (CountryFixupData));
  c.country_list = GetValidCountryList ();
  if (c.country_list == NULL) return NULL;
  c.capitalize_after_colon = fix_after_colon;
  c.warning_list = NULL;
  VisitDescriptorsInSep (sep, &c, CountryFixupDesc);
  VisitFeaturesInSep (sep, &c, CountryFixupFeat);
  return c.warning_list;
}


NLM_EXTERN Boolean FixupCountryQualsWithLog (SeqEntryPtr sep, Boolean fix_after_colon, FILE *log_fp)
{
  CountryFixupData c;

  MemSet (&c, 0, sizeof (CountryFixupData));
  c.log_fp = log_fp;
  c.country_list = GetValidCountryList ();
  if (c.country_list == NULL) return FALSE;
  c.capitalize_after_colon = fix_after_colon;
  c.warning_list = NULL;
  VisitDescriptorsInSep (sep, &c, CountryFixupDesc);
  VisitFeaturesInSep (sep, &c, CountryFixupFeat);
  c.warning_list = ValNodeFree (c.warning_list);
  return c.any_changed;
}


typedef struct qualfixup {
  SourceConstraintPtr scp;
  ReplacePairPtr fix_list;
  Boolean case_counts;
  Boolean whole_word;
  Boolean is_orgmod;
  Uint1   subtype;
  Boolean any_changed;
  FILE *log_fp;
} QualFixupData, PNTR QualFixupPtr;

static void FixupBioSourceQuals (BioSourcePtr biop, Pointer data)
{
  QualFixupPtr qf;
  OrgModPtr    mod;
  SubSourcePtr ssp;
  ReplacePairPtr fix;
  CharPtr        orig = NULL;

  if (biop == NULL || (qf = (QualFixupPtr) data) == NULL 
      || qf->fix_list == NULL
      || !DoesBiosourceMatchConstraint(biop, qf->scp)) {
    return;
  }

  if (qf->is_orgmod) {
    if (biop->org == NULL || biop->org->orgname == NULL) {
      return;
    }
    for (mod = biop->org->orgname->mod; mod != NULL; mod = mod->next) {
      if (mod->subtype == qf->subtype) {
        for (fix = qf->fix_list; fix->find != NULL; fix++) {
          orig = StringSave (mod->subname);
          FindReplaceString (&(mod->subname), fix->find, fix->replace, qf->case_counts, qf->whole_word);
          if (StringCmp (orig, mod->subname) != 0) {
            qf->any_changed = TRUE;
            if (qf->log_fp != NULL) {
              fprintf (qf->log_fp, "Changed '%s' to '%s'\n", orig, mod->subname);
            }
          }
          orig = MemFree (orig);
        }
      }
    }
  } else {
    for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
      if (ssp->subtype == qf->subtype) {
        for (fix = qf->fix_list; fix->find != NULL; fix++) {
          orig = StringSave (ssp->name);
          FindReplaceString (&(ssp->name), fix->find, fix->replace, qf->case_counts, qf->whole_word);
          if (StringCmp (orig, ssp->name) != 0) {
            qf->any_changed = TRUE;
            if (qf->log_fp != NULL) {
              fprintf (qf->log_fp, "Changed '%s' to '%s'\n", orig, ssp->name);
            }
          }
          orig = MemFree (orig);
        }
      }
    }
  }
}


static ReplacePairData mouse_strain_fixes[] = {
  {"129/Sv",   "129/Sv"} ,         
  {"129/SvJ",  "129/SvJ"} ,                     
  {"BALB/c",   "BALB/c"} ,                   
  {"C57BL/6",  "C57BL/6"} ,                  
  {"C57BL/6J", "C57BL/6J"} ,              
  {"CD-1",     "CD-1"} ,              
  {"CZECHII",  "CZECHII"} ,
  {"FVB/N",    "FVB/N"} ,                     
  {"FVB/N-3",  "FVB/N-3"} ,                 
  {"ICR",      "ICR"} ,             
  {"NMRI",     "NMRI"} ,                    
  {"NOD",      "NOD"} , 
  {"C3H",      "C3H"} ,                       
  {"C57BL",    "C57BL"} ,                     
  {"C57BL/6",  "C57BL/6"} ,                    
  {"C57BL/6J", "C57BL/6J" } ,                
  {"DBA/2",    "DBA/2"} ,      
  {NULL, NULL}};

NLM_EXTERN Boolean FixupMouseStrains (SeqEntryPtr sep, FILE *log_fp)
{
  QualFixupData qd;

  MemSet (&qd, 0, sizeof (QualFixupData));

  qd.case_counts = FALSE;
  qd.whole_word = TRUE;
  qd.is_orgmod = TRUE;
  qd.subtype = ORGMOD_strain;
  qd.scp = SourceConstraintNew ();
  qd.scp->constraint = StringConstraintNew ();
  qd.scp->constraint->match_text = StringSave ("Mus musculus");
  qd.scp->constraint->match_location = String_location_starts;
  qd.scp->field1 = ValNodeNew (NULL);
  qd.scp->field1->choice = SourceQualValChoice_textqual;
  qd.scp->field1->data.intvalue = Source_qual_taxname;
  qd.log_fp = log_fp;
  qd.fix_list = mouse_strain_fixes;

  VisitBioSourcesInSep (sep, &qd, FixupBioSourceQuals);
  qd.scp = SourceConstraintFree (qd.scp);
  return qd.any_changed;
}


typedef struct srcqualfixlist {
  Int4 src_qual;
  CharPtr PNTR fix_list;
} SrcQualFixListData, PNTR SrcQualFixListPtr;


static CharPtr src_qual_sex_words[] = {
  "asexual",
  "bisexual",
  "diecious",
  "dioecious",
  "female",
  "hermaphrodite",
  "male",
  "monecious",
  "monoecious",
  "unisexual",
  NULL };

static CharPtr src_qual_host_words[] = {
  "alfalfa",
  "almond",
  "apple",
  "bitter melon",
  "blackberry",
  "blueberry",
  "bovine",
  "brinjal",
  "broad bean",
  "cabbage",
  "caprine",
  "cattle",
  "canine",
  "cantalope",
  "cassava",
  "cauliflower",
  "chicken",
  "chimpanzee",
  "clover",
  "corn",
  "cotton",
  "cow",
  "cowpea",
  "cucumber",
  "dairy cow",
  "dog",
  "duck",
  "equine",
  "feline",
  "fish",
  "fox",
  "goat",
  "goldfish",
  "goose",
  "honeydew",
  "horse",
  "juniper",
  "lily",
  "maize",
  "mango",
  "mangrove",
  "marine sponge",
  "mulberry",
  "mungbean",
  "nematode",
  "ovine",
  "peacock",
  "pear",
  "pepper",
  "pig",
  "pomegranate",
  "porcine",
  "potato",
  "rice",
  "salmon",
  "sesame",
  "sheep",
  "soybean",
  "sponge",
  "squash",
  "sunflower",
  "swine",
  "tomato",
  "turkey",
  "turtle",
  "watermelon",
  "wheat",
  "white clover",
  "wolf",
  "yak",
  NULL };

  static CharPtr src_qual_lab_host_words[] = {
  "porcine",
  "caprine",
  "ovine",
  "cattle",
  "canine",
  "feline",
  "bovine",
  "tomato",
  "pepper",
  "yak",
  "horse",
  "pig",
  "cow",
  "rice",
  "turkey",
  "chicken",
  "sheep",
  "yak",
  "salmon",
  "wolf",
  "nematode",
  "fox",
  "swine",
  "fish",
  "maize",
  "soybean",
  "wheat",
  NULL };

static CharPtr src_qual_isolation_source_words[] = {
  "abdomen",
  "abdominal fluid",
  "acne",
  "activated sludge",
  "agricultural soil",
  "alfalfa",
  "almond",
  "amniotic fluid",
  "apple",
  "biofilm",
  "bitter melon",
  "blackberry",
  "blood",
  "blueberry",
  "bovine",
  "bovine milk",
  "brain",
  "brain abscess",
  "brain tissue",
  "brinjal",
  "bronchial mucosa",
  "bronchoalveolar lavage",
  "buccal mucosa",
  "bursa of fabricus",
  "cabbage",
  "callus",
  "canine",
  "cantalope",
  "caprine",
  "cassava",
  "cattle",
  "cauliflower",
  "cave sediment",
  "cave sediments",
  "cerbrospinal fluid",
  "chicken",
  "chimpanzee",
  "clinical",
  "clinical isolate",
  "clinical isolates",
  "clinical sample",
  "clinical samples",
  "compost",
  "corn",
  "corn rhizosphere",
  "cornea",
  "cotton",
  "cotton rhizosphere",
  "cow",
  "cowpea",
  "cucumber",
  "dairy cow",
  "dairy cow rumen",
  "dog",
  "duck",
  "egg",
  "embryogenic callus",
  "epithelium",
  "equine",
  "esophageal mucosa",
  "fecal sample",
  "fecal samples",
  "feces",
  "feline",
  "fermented food",
  "fish",
  "flooded rice soil",
  "flower",
  "food",
  "food sample",
  "food samples",
  "forest",
  "fox",
  "freshwater stream",
  "fruit",
  "gastric mucosa",
  "gill",
  "gills",
  "goat",
  "goat milk",
  "goldfish",
  "head",
  "head kidney",
  "heart",
  "hepatocyte",
  "honeydew",
  "horse",
  "hot spring",
  "hot springs",
  "intestinal mucosa",
  "intestines",
  "juniper",
  "kidney",
  "lake sediment",
  "lake soil",
  "lake water",
  "leaf",
  "leaves",
  "lily",
  "liver",
  "liver abscess",
  "lung",
  "lymph node",
  "lymphocyte",
  "mammary gland",
  "mango",
  "mangrove sediment",
  "mangrove sediments",
  "manure",
  "marine sediment",
  "marine sediments",
  "marine sponge",
  "meat",
  "milk",
  "mitral valve",
  "mucosa",
  "mucus",
  "mulberry",
  "mungbean",
  "muscle",
  "muscle tissue",
  "nasal mucosa",
  "nasal sample",
  "nasal samples",
  "nasal swab",
  "nasopharyngeal aspirate",
  "nasopharynx",
  "nematode",
  "nodule",
  "nodules",
  "nose swab",
  "olfactory mucosa",
  "oral lexion",
  "oral mucosa",
  "ovary",
  "ovary",
  "oviduct",
  "ovine",
  "paddy soil",
  "patient",
  "pear",
  "pepper",
  "pharnyx",
  "pig",
  "placenta",
  "plasma",
  "pleura",
  "pomegranate",
  "porcine",
  "potato",
  "rhizosphere",
  "rhizosphere soil",
  "rice",
  "rice rhizosphere",
  "rice soil",
  "river sediment",
  "river sediments",
  "river water",
  "root",
  "root nodule",
  "root nodules",
  "root tip",
  "rumen",
  "salivary gland",
  "salmon",
  "seafood",
  "seawater",
  "sediment",
  "sediments",
  "seedling",
  "sera",
  "serum",
  "sesame",
  "sheep",
  "shrimp pond",
  "skeletal muscle",
  "skin",
  "skin lesion",
  "sludge",
  "soil",
  "soybean",
  "spleen",
  "sponge",
  "squash",
  "stem",
  "stomach",
  "stool",
  "stool sample",
  "stool samples",
  "sunflower",
  "swab",
  "swine",
  "testes",
  "testis",
  "throat",
  "throat swab",
  "thymus",
  "tomato",
  "turkey",
  "turtle",
  "urine",
  "uterine mucosa",
  "wastewater",
  "water",
  "watermelon",
  "wheat",
  "whole blood",
  "whole cell/tissue lysate",
  "wolf",
  "wound",
  "wound",
  "yak",
  "yogurt",
  NULL };

static CharPtr src_qual_tissue_type_words[] = {
  "blood",
  "brain tissue",
  "brain",
  "bursa of fabricus",
  "callus",
  "cornea",
  "dairy cow rumen",
  "embryogenic callus",
  "epithelium",
  "flower",
  "fruit",
  "gill",
  "gills",
  "head kidney",
  "heart",
  "intestines",
  "kidney",
  "leaf",
  "leaves",
  "liver",
  "lung",
  "lymph node",
  "mammary gland",
  "mammary gland",
  "muscle tissue",
  "muscle",
  "nodule",
  "nodules",
  "ovary",
  "oviduct",
  "paddy soil",
  "pharnyx",
  "placenta",
  "plasma",
  "root nodule",
  "root nodules",
  "root tip",
  "root",
  "rumen",
  "salivary gland",
  "skeletal muscle",
  "skin",
  "skin",
  "spleen",
  "stem",
  "stomach",
  "testes",
  "testis",
  "thymus",
  "whole blood",
  NULL };

static CharPtr src_qual_dev_stage_words[] = {
  "adult",
  NULL };

static CharPtr src_qual_cell_type_words[] = {
  "hepatocyte",
  "lymphocyte",
  NULL };

static SrcQualFixListData src_qual_fixes[] = {
  {Source_qual_sex, src_qual_sex_words} ,
  {Source_qual_nat_host, src_qual_host_words},
  {Source_qual_isolation_source, src_qual_isolation_source_words},
  {Source_qual_lab_host, src_qual_lab_host_words},
  {Source_qual_tissue_type, src_qual_tissue_type_words},
  {Source_qual_dev_stage, src_qual_dev_stage_words},
  {Source_qual_cell_type, src_qual_cell_type_words},
  {0, NULL}
};

typedef struct srcqualfix {
  Boolean any_change;
  FILE *log_fp;
  CharPtr PNTR fix_list;
  ValNode vn;
} SrcQualFixData, PNTR SrcQualFixPtr;


static void FixSourceQualCaps (BioSourcePtr biop, Pointer data)
{
  CharPtr val, orig;
  SrcQualFixPtr sq;
  Int4 i;
  StringConstraint sd;

  if (biop == NULL || (sq = (SrcQualFixPtr) data) == NULL || sq->fix_list == NULL) {
    return;
  }
  val = GetSourceQualFromBioSource (biop, &(sq->vn), NULL);
  if (val == NULL) {
    return;
  }
  orig = StringSave (val);
  for (i = 0; sq->fix_list[i] != NULL; i++) {
    if (StringICmp (val, sq->fix_list[i]) == 0) {
      val = MemFree (val);
      val = StringSave (sq->fix_list[i]);
    }
  }
  if (StringCmp (orig, val) != 0) {
    MemSet (&sd, 0, sizeof (StringConstraint));
    sd.match_text = orig;
    sd.match_location = String_location_equals;
    if (SetSourceQualInBioSource (biop, &(sq->vn), &sd, val, ExistingTextOption_replace_old)) {
      sq->any_change = TRUE;
      if (sq->log_fp != NULL) {
        fprintf (sq->log_fp, "Changed '%s' to '%s'\n", orig, val);
      }
    }
  }
  orig = MemFree (orig);
  val = MemFree (val);
}


NLM_EXTERN Boolean FixSrcQualCaps (SeqEntryPtr sep, Int4 src_qual, FILE *log_fp)
{
  Int4 i;
  SrcQualFixData sd;

  MemSet (&sd, 0, sizeof (SrcQualFixData));
  sd.log_fp = log_fp;
  sd.any_change = FALSE;
  MemSet (&sd.vn, 0, sizeof (ValNode));
  sd.vn.choice = SourceQualChoice_textqual;

  /* find fix function */
  for (i = 0; src_qual_fixes[i].fix_list != NULL; i++) {
    if (src_qual_fixes[i].src_qual == src_qual) {
      sd.fix_list = src_qual_fixes[i].fix_list;
      sd.vn.data.intvalue = src_qual;     
      VisitBioSourcesInSep (sep, &sd, FixSourceQualCaps);
    }
  }

  return sd.any_change;
}


extern ValNodePtr ListFeaturesInLocation (BioseqPtr bsp, SeqLocPtr slp, Uint1 seqfeatChoice, Uint1 featdefChoice)
{
  ValNodePtr        feat_list = NULL;
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;
  Int4              loc_left, loc_right, tmp;

  if (bsp == NULL || slp == NULL) return NULL;

  loc_left = SeqLocStart (slp);
  loc_right = SeqLocStop (slp);
  if (loc_left > loc_right) {
    tmp = loc_left;
    loc_left = loc_right;
    loc_right = tmp;
  }
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, seqfeatChoice, featdefChoice, &fcontext);
       sfp != NULL && fcontext.left <= loc_right;
       sfp = SeqMgrGetNextFeature (bsp, sfp, seqfeatChoice, featdefChoice, &fcontext))
  {
    if (fcontext.right < loc_left) continue;
    if (SeqLocCompare (sfp->location, slp) == SLC_A_IN_B)
    {
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    }
  }
  return feat_list;
}


extern ValNodePtr ListFeaturesOverlappingLocationEx (BioseqPtr bsp, SeqLocPtr slp, Uint1 seqfeatChoice, Uint1 featdefChoice, ValNodePtr constraint)
{
  ValNodePtr        feat_list = NULL;
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;
  Int4              loc_left, loc_right, tmp;
  Int4              cmp;

  if (bsp == NULL || slp == NULL) return NULL;

  loc_left = SeqLocStart (slp);
  loc_right = SeqLocStop (slp);
  if (loc_left > loc_right) {
    tmp = loc_left;
    loc_left = loc_right;
    loc_right = tmp;
  }
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, seqfeatChoice, featdefChoice, &fcontext);
       sfp != NULL && fcontext.left <= loc_right;
       sfp = SeqMgrGetNextFeature (bsp, sfp, seqfeatChoice, featdefChoice, &fcontext))
  {
    if (!DoesObjectMatchConstraintChoiceSet(OBJ_SEQFEAT, sfp, constraint)) {
      continue;
    }
    cmp = SeqLocCompare (sfp->location, slp);
    if (cmp != SLC_NO_MATCH)
    {
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    }
  }
  return feat_list;
}


extern ValNodePtr ListFeaturesOverlappingLocation (BioseqPtr bsp, SeqLocPtr slp, Uint1 seqfeatChoice, Uint1 featdefChoice)
{
  return ListFeaturesOverlappingLocationEx (bsp, slp, seqfeatChoice, featdefChoice, NULL);
}


static void CDSInSrcFeatCallback (BioseqPtr bsp, Pointer data)
{
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_BIOSRC, &fcontext);
  while (sfp != NULL)
  {
    ValNodeLink ((ValNodePtr PNTR) data, ListFeaturesInLocation (bsp, sfp->location, 0, FEATDEF_CDS));
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_BIOSRC, &fcontext);
  }
}

extern ValNodePtr ListCodingRegionsContainedInSourceFeatures (SeqEntryPtr sep)
{
  ValNodePtr feat_list = NULL;

  VisitBioseqsInSep (sep, &feat_list, CDSInSrcFeatCallback);
  return feat_list;
}


extern void CountNsInSequence (BioseqPtr bsp, Int4Ptr p_total, Int4Ptr p_max_stretch, Boolean expand_gaps)
{
  Int4      ctr, pos, i;
  Char      buf1[51];
  Int4      len = 50, total = 0, max_stretch = 0, this_stretch = 0;
  StreamFlgType flags = STREAM_CORRECT_INVAL;

  if (p_total != NULL) {
    *p_total = 0;
  }
  if (p_max_stretch != NULL) {
    *p_max_stretch = 0;
  }
  if (bsp == NULL) {
    return;
  }

  if (expand_gaps) {
    flags |= STREAM_EXPAND_GAPS;
  }
  pos = 0;
  while (pos < bsp->length) {
    ctr = SeqPortStreamInt (bsp, pos, MIN(pos + len - 1, bsp->length - 1), Seq_strand_plus,
                            flags, (Pointer) buf1, NULL);
    for (i = 0; i < ctr; i++) {
      if (buf1[i] == 'N') {
        total++;
        this_stretch++;
      } else {
        if (this_stretch > max_stretch) {
          max_stretch = this_stretch;
        }
        this_stretch = 0;
      }
    }
    pos += len;
  }
  if (p_total != NULL) {
    *p_total = total;
  }
  if (p_max_stretch != NULL) {
    *p_max_stretch = max_stretch;
  }
}


extern LogInfoPtr OpenLog (CharPtr display_title)
{
  LogInfoPtr lip;
  
  lip = (LogInfoPtr) MemNew (sizeof (LogInfoData));
  if (lip == NULL)
  {
    return NULL;
  }
  TmpNam (lip->path);
  lip->fp = FileOpen (lip->path, "w");
  lip->data_in_log = FALSE;
  lip->display_title = StringSave (display_title);
  return lip;
}

extern LogInfoPtr FreeLog (LogInfoPtr lip)
{
  if (lip != NULL)
  {
    lip->display_title = MemFree (lip->display_title);
    if (lip->fp != NULL)
    {
      FileClose (lip->fp);
      lip->fp = NULL;
      FileRemove (lip->path);
    }
    lip = MemFree (lip);
  }
  return lip;
}


static void ParseOneFromTaxnameToQuals (OrgRefPtr org, CharPtr qual_name, CharPtr start, Int4 val_len)
{
  Int4        q_type, s_type;
  Boolean     found = FALSE;
  OrgModPtr   mod;

  if (val_len > 0) {
    q_type = GetSourceQualTypeByName (qual_name);
    if (q_type > -1) {
      s_type = GetOrgModQualFromSrcQual (q_type, NULL);
      if (s_type > -1) {
        /* look for existing value */
        if (org->orgname != NULL) {
          for (mod = org->orgname->mod; mod != NULL && !found; mod = mod->next) {
            if (mod->subtype == s_type) {
              found = TRUE;
            }
          }
        }
        if (!found) {
          if (org->orgname == NULL) {
            org->orgname = OrgNameNew ();
          }
          mod = OrgModNew ();
          mod->subtype = s_type;
          mod->subname = (CharPtr) MemNew (sizeof (Char) * (val_len));
          StringNCpy (mod->subname, start + 1, val_len - 1);
          mod->subname[val_len - 1] = 0;
          mod->next = org->orgname->mod;
          org->orgname->mod = mod;
        }
      }
    }
  }
}


NLM_EXTERN void ParseTaxNameToQuals (OrgRefPtr org, TextFsaPtr tags)
{
  Char        ch;
  CharPtr     ptr;
  Int4        state;
  ValNodePtr  matches;
  CharPtr     last_hit = NULL, last_pos = NULL;
  Int4        val_len, match_len;

  if (tags == NULL || org == NULL || StringHasNoText (org->taxname)) return;

  if (StringSearch (org->taxname, " x ") != NULL) {
    /* ignore cross, applies only to one parent, do not parse */
    return;
  }
  state = 0;
  ptr = org->taxname;
  ch = *ptr;
  while (ch != '\0') {
    matches = NULL;
    state = TextFsaNext (tags, state, ch, &matches);
    if (matches != NULL && isspace (*(ptr + 1)) && (match_len = StringLen (matches->data.ptrvalue)) > 0
        && (isspace (*(ptr - match_len)) || ispunct (*(ptr - match_len)))) {
      if (last_pos != NULL) {
        val_len = ptr - last_pos - 1 - match_len;
        ParseOneFromTaxnameToQuals (org, last_hit, last_pos + 1, val_len);
      }
      last_pos = ptr;
      last_hit = (CharPtr) matches->data.ptrvalue;
    }
    ptr++;
    ch = *ptr;
  }
  if (last_pos != NULL) {
    val_len = ptr - last_pos;
    ParseOneFromTaxnameToQuals (org, last_hit, last_pos, val_len);
  }
}


static void GetLocusTagPrefixListCallback (SeqFeatPtr sfp, Pointer data)
{
  GeneRefPtr grp;
  CharPtr    cp, prefix;
  Int4       len;

  if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE 
      && (grp = (GeneRefPtr) sfp->data.value.ptrvalue) != NULL
      && (cp = StringChr (grp->locus_tag, '_')) != NULL
      && (len = cp - grp->locus_tag) > 0) {
    prefix = (CharPtr) MemNew (sizeof (Char) * (len + 1));
    StringNCpy (prefix, grp->locus_tag, len);
    prefix[len] = 0;
    ValNodeAddPointer ((ValNodePtr PNTR) data, 0, prefix);
  }
}


NLM_EXTERN ValNodePtr GetLocusTagPrefixList (SeqEntryPtr sep)
{
    ValNodePtr list = NULL;

    VisitFeaturesInSep (sep, &list, GetLocusTagPrefixListCallback);
    list = ValNodeSort (list, SortVnpByString);
    ValNodeUnique (&list, SortVnpByString, ValNodeFreeData);
    return list;
}


static CharPtr RemovableCultureNotes[] = {
 "[uncultured (using universal primers)]",
 "[uncultured (using universal primers) bacterial source]",
 "[cultured bacterial source]",
 "[enrichment culture bacterial source]",
 "[mixed bacterial source (cultured and uncultured)]",
 "[uncultured]; [universal primers]",
 "[mixed bacterial source]",
 "[virus wizard]",
 "[cDNA derived from mRNA, purified viral particles]",
 "[cDNA derived from mRNA, whole cell/tissue lysate]",
 "[cDNA derived from genomic RNA, whole cell/tissue lysate]",
 "[cDNA derived from genomic RNA, purified viral particles]",
 "[universal primers]",
 "[uncultured; wizard]",
 "[uncultured; wizard; spans unknown]",
 "[cultured; wizard]",
 "[cultured; wizard; spans unknown]",
 NULL
};

static CharPtr ReplaceableCultureNotes[] = {
 "[uncultured (with species-specific primers)]",
 "[uncultured]; [amplified with species-specific primers]",
 "[uncultured (using species-specific primers) bacterial source]",
 "[amplified with species-specific primers]",
 NULL
};


static Boolean RemoveCultureNotesFromText (CharPtr PNTR p_txt)
{
  CharPtr txt, cp, src, dst;
  Int4    i, len, extra_len;
  Boolean any_removed = FALSE;

  if (p_txt == NULL || (txt = *p_txt) == NULL) {
    return FALSE;
  }
  for (i = 0; RemovableCultureNotes[i] != NULL; i++) {
    len = StringLen (RemovableCultureNotes[i]);
    cp = StringISearch (txt, RemovableCultureNotes[i]);
    while (cp != NULL) {
      extra_len = StringSpn (cp + len, " ;");
      src = cp + len + extra_len;
      dst = cp;
      while (*src != 0) {
        *dst = *src;
        ++dst;
        ++src;
      }
      *dst = 0;
      any_removed = TRUE;
      cp = StringISearch (txt, RemovableCultureNotes[i]);
    }
  }
  /* remove leading/trailing semicolons */
  TrimSpacesAndSemicolons (txt);

  for (i = 0; ReplaceableCultureNotes[i] != NULL; i++) {
    if (StringICmp (txt, ReplaceableCultureNotes[i]) == 0) {
      *p_txt = MemFree (*p_txt);
      *p_txt = StringSave ("amplified with species-specific primers");
      txt = *p_txt;
      any_removed = TRUE;
      break;
    }
  }
  if (StringHasNoText (txt)) {
    *p_txt = MemFree (*p_txt);
    any_removed = TRUE;
  }
  return any_removed;
}


static void RemoveCultureNotesBioSourceCallback (BioSourcePtr biop, Pointer data)
{
  BoolPtr p_rval;
  Boolean      rval = FALSE;
  SubSourcePtr ssp, ssp_prev = NULL, ssp_next;

  if (biop == NULL) {
    return;
  }
  p_rval = (BoolPtr) data;

  for (ssp = biop->subtype; ssp != NULL; ssp = ssp_next) {
    ssp_next = ssp->next;
    if (ssp->subtype == 255) {
      rval |= RemoveCultureNotesFromText(&(ssp->name));
      if (StringHasNoText (ssp->name)) {
        ssp->next = NULL;
        ssp = SubSourceFree (ssp);
        if (ssp_prev == NULL) {
          biop->subtype = ssp_next;
        } else {
          ssp_prev->next = ssp_next;
        }
      } else {
        ssp_prev = ssp;
      }
    } else {
      ssp_prev = ssp;
    }
  }

  if (p_rval != NULL) {
    *p_rval |= rval;
  }
}


NLM_EXTERN Boolean RemoveCultureNotes (SeqEntryPtr sep)
{
  Boolean rval = FALSE;

  VisitBioSourcesInSep (sep, &rval, RemoveCultureNotesBioSourceCallback);
  return rval;
}


static CharPtr s_CorrectProductCaps[] = {
  "ABC",
  "AAA",
  "ATP",
  "ATPase",
  "A/G",
  "AMP",
  "CDP",
  "coproporphyrinogen III",
  "cytochrome BD",
  "cytochrome C",
  "cytochrome C2",
  "cytochrome C550",
  "cytochrome D",
  "cytochrome O",
  "cytochrome P450",
  "cytochrome P460",
  "D-alanine",
  "D-alanyl",
  "D-amino",
  "D-beta",
  "D-cysteine",
  "D-lactate",
  "D-ribulose",
  "D-xylulose",
  "endonuclease I",
  "endonuclease II",
  "endonuclease III",
  "endonuclease V",
  "EPS I",
  "Fe-S",
  "ferredoxin I",
  "ferredoxin II",
  "GTP",
  "GTPase",
  "H+",
  "hemolysin I",
  "hemolysin II",
  "hemolysin III",
  "L-allo",
  "L-arabinose",
  "L-asparaginase",
  "L-aspartate",
  "L-carnitine",
  "L-fuculose",
  "L-glutamine",
  "L-histidinol",
  "L-isoaspartate",
  "L-serine",
  "MFS",
  "FAD/NAD(P)",
  "MCP",
  "Mg+",
  "Mg chelatase",
  "Mg-protoporphyrin IX",
  "N(5)",
  "N,N-",
  "N-(",
  "N-acetyl",
  "N-acyl",
  "N-carb",
  "N-form",
  "N-iso",
  "N-succ",
  "NADP",
  "Na+/H+",
  "NAD",
  "NAD(P)",
  "NADPH",
  "O-sial",
  "O-succ",
  "pH",
  "ribonuclease BN",
  "ribonuclease D",
  "ribonuclease E",
  "ribonuclease G",
  "ribonuclease H",
  "ribonuclease I",
  "ribonuclease II",
  "ribonuclease III",
  "ribonuclease P",
  "ribonuclease PH",
  "ribonuclease R",
  "RNAse",
  "S-adeno",
  "type I",
  "type II",
  "type III",
  "type IV",
  "type V",
  "type VI",
  "UDP",
  "UDP-N",
  "Zn",
  NULL};

NLM_EXTERN void FixProductWordCapitalization (CharPtr PNTR pProduct)
{
  Int4 i;

  if (pProduct == NULL || *pProduct == NULL) {
    return;
  }

  for (i = 0; s_CorrectProductCaps[i] != NULL; i++) {
    FindReplaceString (pProduct, s_CorrectProductCaps[i], s_CorrectProductCaps[i], FALSE, TRUE);
  }
}


NLM_EXTERN Boolean IsNCBIFileID (SeqIdPtr sip)
{
  DbtagPtr dbt;

  if (sip == NULL || sip->choice != SEQID_GENERAL) return FALSE;
  dbt = (DbtagPtr) sip->data.ptrvalue;
  if (dbt == NULL) return FALSE;
  if (StringCmp (dbt->db, "NCBIFILE") == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


