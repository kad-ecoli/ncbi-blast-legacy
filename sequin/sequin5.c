/*   sequin5.c
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
* File Name:  sequin5.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   8/26/97
*
* $Revision: 6.767 $
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

#define USE_BLAST3

#include "sequin.h"
#include <document.h>
#include <biosrc.h>
#include <jzcoll.h>
#include <objsub.h>
#ifdef USE_BLAST3
#include <netblap3.h> 
#else
#include <netblap2.h>
#include <blast2.h>
#endif
#include <dust.h>
#include <pobutil.h>
#include <saledit.h>
#include <salstruc.h>
#include <salfiles.h>
#include <salign.h>
#include <salsap.h>
#include <gbfeat.h>
#include <gbftdef.h>
#include <satutil.h>
#include <rpsutil.h>
#include <subutil.h>
#include <explore.h>
#include <import.h>
#include <salutil.h>

#include <asn2gnbp.h> /* added for  parse to flatfile */
#include <sqnutils.h> /* added to include SeqEdTranslateOneCDS */
#include <seqpanel.h>
#include <salpanel.h>
#include <tax3api.h> /* added for specific-host corrections */
#include <findrepl.h>
#include <valid.h> /* added for latloncountry conflict checking */
#include <vsmutil.h>

#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>


static void CommonLaunchBioseqViewer (SeqEntryPtr sep, CharPtr path, Boolean directToEditor)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Pointer       dataptr;
  Uint2         datatype;
  Uint2         entityID;
  Int2          handled;

  if (sep != NULL) {
    if (IS_Bioseq (sep)) {
      datatype = OBJ_BIOSEQ;
    } else if (IS_Bioseq_set (sep)) {
      datatype = OBJ_BIOSEQSET;
    } else {
      sep = SeqEntryFree (sep);
      return;
    }
    dataptr = (Pointer) sep->data.ptrvalue;
    entityID = ObjMgrRegister (datatype, dataptr);
    if (dataptr != NULL && entityID > 0) {
      if (datatype == OBJ_SEQSUB || datatype == OBJ_SEQENTRY ||
          datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET) {
        WatchCursor ();
        sep = GetTopSeqEntryForEntityID (entityID);
        if (sep == NULL) {
          sep = SeqEntryNew ();
          if (sep != NULL) {
            if (datatype == OBJ_BIOSEQ) {
              bsp = (BioseqPtr) dataptr;
              sep->choice = 1;
              sep->data.ptrvalue = bsp;
              SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
            } else if (datatype == OBJ_BIOSEQSET) {
              bssp = (BioseqSetPtr) dataptr;
              sep->choice = 2;
              sep->data.ptrvalue = bssp;
              SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) bssp, sep);
            } else {
              sep = SeqEntryFree (sep);
            }
          }
          sep = GetTopSeqEntryForEntityID (entityID);
        }
        if (sep != NULL) {
          SeqEntryPack (sep);
        }
        if (directToEditor) {
          handled = GatherProcLaunch (OMPROC_EDIT, FALSE, entityID, 1,
                                      OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
        } else {
          if (sep != NULL) {
            if (! leaveAsOldAsn) {
              MySeqEntryToAsn3 (sep, TRUE, FALSE, FALSE);
            }
          }
          seqviewprocs.filepath = path;
          seqviewprocs.forceSeparateViewer = TRUE;
          handled = GatherProcLaunch (OMPROC_VIEW, FALSE, entityID, 1,
                                      OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
          seqviewprocs.filepath = NULL;
        }
        ArrowCursor ();
        if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
          Message (MSG_FATAL, "Unable to launch viewer.");
          SeqEntryFree (sep);
          return;
        } else {
          SendHelpScrollMessage (helpForm, "Editing the Record", NULL);
        }
        ObjMgrSetOptions (OM_OPT_FREE_IF_NO_VIEW, entityID);
        ObjMgrSetDirtyFlag (entityID, TRUE);
      } else {
        Message (MSG_ERROR, "Unable to process object type %d.", (int) datatype);
        ObjMgrDelete (datatype, dataptr);
      }
    }
  }
}

extern void FastaNucDirectToSeqEdProc (IteM i)

{
  FILE         *fp;
  Char         path [PATH_MAX];
  SeqEntryPtr  sep;

  if (! GetInputFileName (path, sizeof (path), "", "TEXT")) return;
  fp = FileOpen (path, "r");
  if (fp == NULL) return;
  sep = SequinFastaToSeqEntryEx (fp, TRUE, NULL, TRUE, NULL);
  FileClose (fp);
  CommonLaunchBioseqViewer (sep, path, TRUE);
}

#define BUFSIZE 128

static CharPtr GetSequenceString (BioseqPtr bsp)

{
  Char          buf [BUFSIZE + 5];
  Int2          j;
  ByteStorePtr  raw;
  CharPtr       sequence = NULL;
  SeqPortPtr    spp;
  Uint1         u1Residue;

  if (bsp == NULL) return NULL;
  if (! ISA_na (bsp->mol)) return NULL;
  spp = SeqPortNew (bsp, 0, -1, 0, Seq_code_iupacna);
  if (spp == NULL) return NULL;
  raw = BSNew (bsp->length + 5);
  if (raw == NULL) {
    SeqPortFree (spp);
    return NULL;
  }

  j = 0;
  buf [0] = '\0';
  while ((u1Residue = SeqPortGetResidue (spp)) != SEQPORT_EOF) {
    if (IS_residue (u1Residue)) {
      u1Residue = TO_UPPER(u1Residue);
      if (u1Residue == 'U') {
        u1Residue = 'T';
      }
      buf [j] = (Char) u1Residue;
      j++;
      if (j >= BUFSIZE) {
        BSWrite (raw, buf, j * sizeof (Char));
        j = 0;
      }
    }
  }
  buf [j] = '\0';
  if (j >= 0) {
    BSWrite (raw, buf, j * sizeof (Char));
  }
  sequence = BSMerge (raw, NULL);
  BSFree (raw);
  SeqPortFree (spp);
  return sequence;
}

static Boolean LookForSearchString (CharPtr title, CharPtr str, CharPtr tmp, size_t maxsize)

{
  CharPtr  ptr;

  ptr = StringISearch (title, str);
  if (ptr != NULL) {
    StringNCpy_0 (tmp, ptr + StringLen (str), maxsize);
     ptr = StringChr (tmp, ']');
     if (ptr != NULL) {
       *ptr = '\0';
     }
    return TRUE;
  }
  return FALSE;
}

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
          if (/* SeqLocCompare (gep->slp, sfp->location) != SLC_NO_MATCH */
              SeqLocAinB (gep->slp, sfp->location) <= 0) {
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

static Boolean OligoExtendGene (GeneRefPtr grp, SeqEntryPtr nsep, SeqLocPtr slp)

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

static void AddGBQual (SeqFeatPtr sfp, CharPtr qual, CharPtr val)

{
  GBQualPtr  gbq;
  GBQualPtr  last;

  if (sfp == NULL || StringHasNoText (qual) || StringHasNoText (val)) return;
  gbq = GBQualNew ();
  if (gbq == NULL) return;
  gbq->qual = StringSave (qual);
  gbq->val = StringSave (val);
  if (sfp->qual == NULL) {
    sfp->qual = gbq;
  } else {
    last = sfp->qual;
    while (last->next != NULL) {
      last = last->next;
    }
    last->next = gbq;
  }
}

static void ProcessOligoDefline (SeqFeatPtr sfp, CharPtr title, SeqEntryPtr nsep)

{
  BioseqPtr   bsp;
  GeneRefPtr  grp;
  SeqFeatPtr  gsfp;
  SeqLocPtr   gslp;
  Boolean     hasNulls;
  SeqIdPtr    sip;
  Char        tmp [256];

  if (sfp == NULL || StringHasNoText (title) || nsep == NULL) return;
  if (LookForSearchString (title, "[comment=", tmp, sizeof (tmp) - 1)) {
    sfp->comment = StringSave (tmp);
  }
  if (LookForSearchString (title, "[stdname=", tmp, sizeof (tmp) - 1)) {
    AddGBQual (sfp, "standard_name", tmp);
  }
  if (LookForSearchString (title, "[conditions=", tmp, sizeof (tmp) - 1)) {
    AddGBQual (sfp, "PCR_conditions", tmp);
  }
  if (LookForSearchString (title, "[gene=", tmp, sizeof (tmp) - 1)) {
    grp = CreateNewGeneRef (tmp, NULL, NULL, FALSE);
    if (grp != NULL) {
      /*
      slp = AsnIoMemCopy ((Pointer) sfp->location,
                          (AsnReadFunc) SeqLocAsnRead,
                          (AsnWriteFunc) SeqLocAsnWrite);
      if (slp != NULL || slp->choice == SEQLOC_INT) {
        sintp = (SeqIntPtr) slp->data.ptrvalue;
        if (sintp != NULL) {
          if (sintp->strand == Seq_strand_minus) {
            sintp->strand = Seq_strand_plus;
          } else if (sintp->strand == Seq_strand_plus) {
            sintp->strand = Seq_strand_minus;
          }
        }
      }
      */
      if (OligoExtendGene (grp, nsep, sfp->location) /* || OligoExtendGene (grp, nsep, slp) */) {
        grp = GeneRefFree (grp);
      } else {
        gsfp = CreateNewFeature (nsep, NULL, SEQFEAT_GENE, NULL);
        if (gsfp != NULL) {
          gsfp->data.value.ptrvalue = (Pointer) grp;
          gsfp->location = SeqLocFree (gsfp->location);
          gsfp->location = AsnIoMemCopy ((Pointer) sfp->location,
                                        (AsnReadFunc) SeqLocAsnRead,
                                        (AsnWriteFunc) SeqLocAsnWrite);
          sip = SeqLocId (gsfp->location);
          if (sip != NULL) {
            bsp = BioseqFind (sip);
          } else {
            bsp = (BioseqPtr) nsep->data.ptrvalue;
          }
          if (bsp != NULL) {
            gslp = SeqLocMerge (bsp, gsfp->location, NULL, TRUE, FALSE, FALSE);
            if (gslp != NULL) {
              gsfp->location = SeqLocFree (gsfp->location);
              gsfp->location = gslp;
              if (bsp->repr == Seq_repr_seg) {
                gslp = SegLocToPartsEx (bsp, gsfp->location, TRUE);
                gsfp->location = SeqLocFree (gsfp->location);
                gsfp->location = gslp;
                hasNulls = LocationHasNullsBetween (gsfp->location);
                gsfp->partial = (gsfp->partial || hasNulls);
              }
              FreeAllFuzz (gslp);
            }
          }
        }
      }
      /* SeqLocFree (slp); */
    }
  }
}

static void ProcessOligo (SeqEntryPtr sep, CharPtr sequence, CharPtr title)

{
  BioseqSetPtr  bssp;
  ImpFeatPtr    ifp;
  BioseqPtr     nbsp;
  SeqEntryPtr   nsep;
  SeqFeatPtr    sfp;
  SeqLocPtr     slp;

  if (sep == NULL || sequence == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)))) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        ProcessOligo (sep, sequence, title);
      }
      return;
    }
  }
  nsep = FindNucSeqEntry (sep);
  if (nsep == NULL) return;
  if (! IS_Bioseq (nsep)) return;
  nbsp = (BioseqPtr) nsep->data.ptrvalue;
  if (nbsp == NULL) return;
  slp = StringSearchInBioseq (nbsp->id, sequence);
  if (slp == NULL) return;
  ifp = ImpFeatNew ();
  if (ifp != NULL) {
    ifp->key = StringSave ("primer_bind");
    sfp = CreateNewFeature (nsep, NULL, SEQFEAT_IMP, NULL);
    if (sfp != NULL) {
      sfp->data.value.ptrvalue = (Pointer) ifp;
      sfp->location = SeqLocFree (sfp->location);
      sfp->location = slp;
      ProcessOligoDefline (sfp, title, nsep);
    }
  }
}

extern void ParseInOligoPrimers (IteM i)

{
  BaseFormPtr  bfp;
  BioseqPtr    bsp;
  CharPtr      extension;
  FILE         *fp;
  SeqEntryPtr  nextsep;
  Char         path [PATH_MAX];
  SeqEntryPtr  sep;
  CharPtr      sequence;
  CharPtr      title;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  extension = GetAppProperty ("FastaNucExtension");
  if (! GetInputFileName (path, sizeof (path), extension, "TEXT")) return;
  fp = FileOpen (path, "r");
  if (fp == NULL) return;
  WatchCursor ();
  Update ();

  nextsep = SequinFastaToSeqEntryEx (fp, TRUE, NULL, FALSE, NULL);
  while (nextsep != NULL) {
    if (IS_Bioseq (nextsep)) {
      bsp = (BioseqPtr) nextsep->data.ptrvalue;
      sequence = GetSequenceString (bsp);
      title = SeqEntryGetTitle (nextsep);
      ProcessOligo (sep, sequence, title);
      MemFree (sequence);
    }
    SeqEntryFree (nextsep);
    nextsep = SequinFastaToSeqEntryEx (fp, TRUE, NULL, FALSE, NULL);
  }

  FileClose (fp);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
}

typedef struct edlocdata {
  FORM_MESSAGE_BLOCK
  SeqIdPtr      sip;
  BioseqPtr     bsp;
  TexT          locusname;
  PrompT        accnnumber;
  Boolean       refgene;
} EdLocData, PNTR EdLocPtr;

static CharPtr AllToUpper (CharPtr str)

{
  Char          ch;
  CharPtr       ptr;

  if (str != NULL) {
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      *ptr = TO_UPPER (ch);
      ptr++;
      ch = *ptr;
    }
  }
  return str;
}

static void EditLocusActnProc (ForM f)

{
  BioseqPtr     bsp;
  Uint1         choice = SEQID_GENBANK;
  Int2          count;
  SeqIdPtr      id;
  EdLocPtr      elp;
  Int2          len;
  Char          num [16];
  BioseqPtr     part;
  CharPtr       ptr;
  Char          str [64];
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  Char          tmp [76];
  TextSeqIdPtr  tsip;
  ValNode       vn;

  elp = (EdLocPtr) GetObjectExtra (f);
  if (elp != NULL) {
    sip = elp->sip;
    bsp = elp->bsp;
    if (elp->refgene) {
      choice = SEQID_OTHER;
    }
    if (indexerVersion && sip == NULL && bsp != NULL && bsp->id != NULL) {
      id = bsp->id;
      while (id->next != NULL) {
        id = sip->next;
      }
      tsip = TextSeqIdNew ();
      if (tsip != NULL) {
        sip = ValNodeNew (NULL);
        if (sip != NULL) {
          sip->choice = choice;
          sip->data.ptrvalue = (Pointer) tsip;
        } else {
          TextSeqIdFree (tsip);
        }
        id->next = sip;
        SeqMgrReplaceInBioseqIndex (bsp);
      }
    }
    if (sip == NULL) return;
    tsip = (TextSeqIdPtr) sip->data.ptrvalue;
    if (tsip == NULL) return;
    tsip->accession = MemFree (tsip->accession);
    GetTitle (elp->accnnumber, str, sizeof (str) - 1);
    TrimSpacesAroundString (str);
    if (! StringHasNoText (str)) {
      tsip->accession = StringSave (AllToUpper (str));
    }
    tsip->name = MemFree (tsip->name);
    GetTitle (elp->locusname, str, sizeof (str) - 1);
    if (str [0] != '\0') {
      if (elp->refgene) {
        tsip->name = StringSave (str);
      } else {
        tsip->name = StringSave (AllToUpper (str));
      }
    }
    if (bsp != NULL && bsp->repr == Seq_repr_seg && bsp->seq_ext != NULL) {
      if (StringNCmp (str, "SEG_", 4) != 0) {
        sprintf (tmp, "SEG_%s", str);
        tsip->name = MemFree (tsip->name);
        tsip->name = StringSave (tmp);
        SeqMgrReplaceInBioseqIndex (bsp);
      }
      vn.choice = SEQLOC_MIX;
      vn.next = NULL;
      vn.data.ptrvalue = bsp->seq_ext;
      count = 0;
      slp = SeqLocFindNext (&vn, NULL);
      while (slp != NULL) {
        if (slp->choice != SEQLOC_NULL) {
          count++;
        }
        slp = SeqLocFindNext (&vn, slp);
      }
      if (count < 10) {
        len = 1;
      } else if (count < 100) {
        len = 2;
      } else {
        len = 3;
      }
      count = 0;
      slp = SeqLocFindNext (&vn, NULL);
      while (slp != NULL) {
        if (slp->choice != SEQLOC_NULL) {
          count++;
          sprintf (num, "%*d", (int) len, (int) count);
          for (ptr = num; *ptr != '\0'; ptr++) {
            if (*ptr == ' ') {
              *ptr = '0';
            }
          }
          sprintf (tmp, "%s%s", str, num);
          sip = SeqLocId (slp);
          if (sip != NULL) {
            part = BioseqFind (sip);
            if (part != NULL && part->id != NULL) {
              id = part->id;
              while (id->choice != choice && id->next != NULL) {
                id = id->next;
              }
              if (id->choice != choice) {
                tsip = TextSeqIdNew ();
                if (tsip != NULL) {
                  sip = ValNodeNew (NULL);
                  if (sip != NULL) {
                    sip->choice = choice;
                    sip->data.ptrvalue = (Pointer) tsip;
                  } else {
                    TextSeqIdFree (tsip);
                  }
                  id->next = sip;
                  id = sip;
                }
              }
              if (id->choice == choice) {
                tsip = (TextSeqIdPtr) id->data.ptrvalue;
                if (tsip != NULL) {
                  tsip->name = MemFree (tsip->name);
                  if (tmp [0] != '\0') {
                    if (elp->refgene) {
                      tsip->name = StringSave (tmp);
                    } else {
                      tsip->name = StringSave (AllToUpper (tmp));
                    }
                  }
                }
              }
              SeqMgrReplaceInBioseqIndex (part);
            }
          }
        }
        slp = SeqLocFindNext (&vn, slp);
      }
    }
    ObjMgrSetDirtyFlag (elp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, elp->input_entityID, 0, 0);
  }
}

static void EditLocusMessage (ForM f, Int2 mssg)

{
  EdLocPtr  elp;

  elp = (EdLocPtr) GetObjectExtra (f);
  if (elp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (elp->appmessage != NULL) {
          elp->appmessage (f, mssg);
        }
        break;
    }
  }
}

#define NUM_ORDER 16

static SeqIdPtr SeqIdFindGenBankOrOther (SeqIdPtr sip)

{
  Uint1  order [NUM_ORDER];

  SeqIdBestRank (order, NUM_ORDER);
  order [SEQID_LOCAL] = 255;
  order [SEQID_GENBANK] = 1;
  order [SEQID_EMBL] = 2;
  order [SEQID_PIR] = 255;
  order [SEQID_SWISSPROT] = 255;
  order [SEQID_DDBJ] = 3;
  order [SEQID_PRF] = 255;
  order [SEQID_PDB] = 255;
  order [SEQID_PATENT] = 255;
  order [SEQID_OTHER] = 2;
  order [SEQID_GENERAL] = 255;
  order [SEQID_GIBBSQ] = 255;
  order [SEQID_GIBBMT] = 255;
  order [SEQID_GIIM] = 255;
  order [SEQID_GI] = 255;
  return SeqIdSelect (sip, order, NUM_ORDER);
}

void EditLocusProc (IteM i)

{
  ButtoN             b;
  BaseFormPtr        bfp;
  BioseqPtr          bsp;
  GrouP              c;
  EdLocPtr           elp;
  GrouP              g;
  StdEditorProcsPtr  sepp;
  SeqIdPtr           sip;
  TextSeqIdPtr       tsip;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  if (bfp->input_itemtype != OBJ_BIOSEQ) return;
  bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp == NULL) return;
  sip = SeqIdFindGenBankOrOther (bsp->id);
  /*if (sip == NULL) return;*/
  elp = (EdLocPtr) MemNew (sizeof (EdLocData));
  if (elp == NULL) return;
  w = MovableModalWindow (-50, -33, -10, -10,
                          "Locus and Accession Number",
                          StdCloseWindowProc);
  SetObjectExtra (w, elp, StdCleanupFormProc);
  elp->form = (ForM) w;
  elp->actproc = EditLocusActnProc;
  elp->formmessage = EditLocusMessage;

  elp->input_entityID = bfp->input_entityID;
  elp->input_itemID = bfp->input_itemID;
  elp->input_itemtype = bfp->input_itemtype;

#ifndef WIN_MAC
  CreateStdEditorFormMenus (w);
#endif

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    elp->appmessage = sepp->handleMessages;
  }

  g = HiddenGroup (w, 2, 0, NULL);
  elp->sip = sip;
  elp->bsp = bsp;
  StaticPrompt (g, "Locus: ", 0, dialogTextHeight, systemFont, 'l');
  elp->locusname = DialogText (g, "", 20, NULL);
  StaticPrompt (g, "Accession: ", 0, stdLineHeight, systemFont, 'l');
  elp->accnnumber = StaticPrompt (g, "", 20 * stdCharWidth, stdLineHeight, systemFont, 'l');
  c = HiddenGroup (w, 2, 0, NULL);
  b = DefaultButton (c, "Accept", StdAcceptFormButtonProc);
  SetObjectExtra (b, elp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  if (sip != NULL) {
    tsip = (TextSeqIdPtr) sip->data.ptrvalue;
    if (tsip != NULL) {
      SetTitle (elp->locusname, tsip->name);
      SetTitle (elp->accnnumber, tsip->accession);
    }
    if (sip->choice == SEQID_OTHER) {
      elp->refgene = TRUE;
    }
  }
  RealizeWindow (w);
  Select (elp->locusname);
  Show (w);
  Select (w);
  Update ();
}

typedef struct accntolocal 
{
  Boolean make_secondary;
  Boolean convert_prots;
  Boolean convert_nucs;
} AccnToLocalData, PNTR AccnToLocalPtr;

static SeqIdPtr LocalSeqIdFromAccession (SeqIdPtr sip, CharPtr str, Int4 str_size)
{
  TextSeqIdPtr  tsip;
  CharPtr       text;
  SeqIdPtr      new_sip;
  ObjectIdPtr   oip;

  if (sip == NULL || str == NULL || str_size < 1)
  {
    return NULL;
  }
  if (sip->choice != SEQID_GENBANK 
      && sip->choice != SEQID_EMBL 
      && sip->choice != SEQID_DDBJ)
  {
    return NULL;
  }
  if (sip->data.ptrvalue == NULL)
  {
    return NULL;
  }
  new_sip = ValNodeNew (NULL);
  if (new_sip == NULL)
  {
    return NULL;
  }
  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
  text = NULL;
  str [0] = '\0';
  if (tsip->accession != NULL) {
    text = tsip->accession;
  } else if (tsip->name != NULL) {
    text = tsip->name;
  }
  if (text != NULL) {
    StringNCpy_0 (str, text, str_size);
    oip = ObjectIdNew ();
    if (oip != NULL) {
      oip->str = StringSave (text);
      new_sip->choice = SEQID_LOCAL;
      new_sip->data.ptrvalue = (Pointer) oip;
    }
  }
  return new_sip;
}

static void SeqIdSwap (SeqIdPtr sip1, SeqIdPtr sip2)
{
  ValNode       tmp_sip;

  if (sip1 == NULL || sip2 == NULL)
  {
    return;
  }
  
  tmp_sip.choice = sip1->choice;
  tmp_sip.data.ptrvalue = sip1->data.ptrvalue;
  sip1->choice = sip2->choice;
  sip1->data.ptrvalue = sip2->data.ptrvalue;
  sip2->choice = tmp_sip.choice;
  sip2->data.ptrvalue = tmp_sip.data.ptrvalue;
}

static void ConvertSeqIdListToLocalID (SeqIdPtr sip, SeqEntryPtr sep, AccnToLocalPtr atlp)

{
  GBBlockPtr    gbp;
  Char          str [32];
  ValNodePtr    vnp;
  BioseqPtr     bsp;
  SeqIdPtr      new_sip = NULL;
  
  if (atlp == NULL) return;

  while (sip != NULL) 
  {
    bsp = BioseqFind (sip);
    if (bsp == NULL)
    {
      new_sip = LocalSeqIdFromAccession (sip, str, sizeof (str));
      if (new_sip != NULL)
      {
        bsp = BioseqFind (new_sip);
      }
    }
    if (bsp == NULL)
    {
      /* have to skip */
    }
    else if (ISA_na (bsp->mol) && !atlp->convert_nucs)
    {
      /* don't convert the nucleotide ID */
    }
    else if (ISA_aa (bsp->mol) && ! atlp->convert_prots)
    {
      /* don't convert the protein ID */
    }
    else if ((sip->choice == SEQID_GENBANK 
          || sip->choice == SEQID_EMBL 
          || sip->choice == SEQID_DDBJ)
          && sip->data.ptrvalue != NULL)
    {
      new_sip = LocalSeqIdFromAccession (sip, str, sizeof (str));
      if (new_sip != NULL)
      {
        /* replace old sip values with new  - new_sip will now be ready for freeing */
        SeqIdSwap (sip, new_sip);
                
        /* make secondary if needed */
        if (atlp->make_secondary && sep != NULL)
        {
          vnp = SeqEntryGetSeqDescr (sep, Seq_descr_genbank, NULL);
          if (vnp == NULL) {
            vnp = CreateNewDescriptor (sep, Seq_descr_genbank);
            if (vnp != NULL) {
              gbp = GBBlockNew ();
              vnp->data.ptrvalue = (Pointer) gbp;
            }
          }
          if (vnp != NULL) {
            gbp = (GBBlockPtr) vnp->data.ptrvalue;
            if (gbp != NULL) {
              for (vnp = gbp->extra_accessions;
                   vnp != NULL && StringICmp (vnp->data.ptrvalue, str) != 0;
                   vnp = vnp->next) continue;
              if (vnp == NULL) {
                vnp = ValNodeCopyStr (&(gbp->extra_accessions), 0, str);
              }
            }
          }
        }
      }
    }
    new_sip = SeqIdFree (new_sip);
    sip = sip->next;
  }
}

static void ConvertSeqLocListToLocalID (SeqLocPtr slp, AccnToLocalPtr atlp)

{
  SeqLocPtr      loc;
  PackSeqPntPtr  psp;
  SeqBondPtr     sbp;
  SeqIntPtr      sinp;
  SeqIdPtr       sip;
  SeqPntPtr      spp;

  while (slp != NULL) {
    switch (slp->choice) {
      case SEQLOC_NULL :
        break;
      case SEQLOC_EMPTY :
      case SEQLOC_WHOLE :
        sip = (SeqIdPtr) slp->data.ptrvalue;
        ConvertSeqIdListToLocalID (sip, NULL, atlp);
        break;
      case SEQLOC_INT :
        sinp = (SeqIntPtr) slp->data.ptrvalue;
        if (sinp != NULL) {
          sip = sinp->id;
          ConvertSeqIdListToLocalID (sip, NULL, atlp);
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL) {
          sip = spp->id;
          ConvertSeqIdListToLocalID (sip, NULL, atlp);
        }
        break;
      case SEQLOC_PACKED_PNT :
        psp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (psp != NULL) {
          sip = psp->id;
          ConvertSeqIdListToLocalID (sip, NULL, atlp);
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        loc = (SeqLocPtr) slp->data.ptrvalue;
        while (loc != NULL) {
          ConvertSeqLocListToLocalID (loc, atlp);
          loc = loc->next;
        }
        break;
      case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          spp = (SeqPntPtr) sbp->a;
          if (spp != NULL) {
            sip = spp->id;
            ConvertSeqIdListToLocalID (sip, NULL, atlp);
          }
          spp = (SeqPntPtr) sbp->b;
          if (spp != NULL) {
            sip = spp->id;
            ConvertSeqIdListToLocalID (sip, NULL, atlp);
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
}

static Boolean ChangeAccessionToLocalID (GatherContextPtr gcp)

{
  BioseqPtr      bsp;
  AccnToLocalPtr atlp;
  SeqEntryPtr    sep;
  SeqFeatPtr     sfp;
  SeqIdPtr       sip;
  SeqLocPtr      slp;

  if (gcp == NULL || gcp->userdata == NULL) return TRUE;
  atlp = (AccnToLocalPtr) gcp->userdata;
  if (gcp->thisitem == NULL) return TRUE;
  switch (gcp->thistype) {
    case OBJ_BIOSEQ :
      sep = NULL;
      bsp = (BioseqPtr) gcp->thisitem;
      if (atlp->make_secondary)
      {
        sep = SeqMgrGetSeqEntryForData (bsp);
      }
      sip = bsp->id;
      ConvertSeqIdListToLocalID (sip, sep, atlp);
      SeqMgrReplaceInBioseqIndex (bsp);
      break;
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) gcp->thisitem;
      slp = sfp->location;
      ConvertSeqLocListToLocalID (slp, atlp);
      slp = sfp->product;
      ConvertSeqLocListToLocalID (slp, atlp);
      break;
    default :
      break;
  }
  return TRUE;
}

static void ConvertToLocalProc (IteM i, Boolean convert_nucs, Boolean convert_prots)

{
  MsgAnswer       ans;
  BaseFormPtr     bfp;
  GatherScope     gs;
  AccnToLocalData atld;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  if (bfp->input_itemtype != OBJ_BIOSEQ) return;
  ans = Message (MSG_OKC, "Are you sure you want to convert accessions to local IDs?");
  if (ans == ANS_CANCEL) return;
  ans = Message (MSG_OKC, "Are you REALLY sure you want to convert accessions to local IDs?");
  if (ans == ANS_CANCEL) return;
  ans = Message (MSG_YN, "Do you want to make the original accessions secondary?");
  atld.make_secondary = (Boolean) (ans == ANS_YES);
  atld.convert_nucs = convert_nucs;
  atld.convert_prots = convert_prots;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  MemSet((Pointer) (gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherEntity (bfp->input_entityID, (Pointer) &atld, ChangeAccessionToLocalID, &gs);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

extern void ConvertToLocalProcOnlyNucs (IteM i)
{
  ConvertToLocalProc (i, TRUE, FALSE);
}

extern void ConvertToLocalProcOnlyProts (IteM i)
{
  ConvertToLocalProc (i, FALSE, TRUE);
}

extern void ConvertToLocalProcAll (IteM i)
{
  ConvertToLocalProc (i, TRUE, TRUE);
}
 
static SeqIdPtr SeqIdFindBestForPromotion (SeqIdPtr sip)

{
  return SeqIdFindBest (sip, 0);
}

static SeqIdPtr SeqIdFindWorstForPromotion (SeqIdPtr sip)

{
  return SeqIdFindWorst (sip);
}

static void PromoteSeqId (SeqIdPtr sip, Boolean alsoCheckLocalAccn, Boolean findWorst)

{
  SeqIdPtr     bestid, newid, oldid;
  BioseqPtr    bsp;
  ObjectIdPtr  oip;
  TextSeqId    tsi;
  SeqId        vn;

  bsp = BioseqFind (sip);
  if (bsp == NULL && alsoCheckLocalAccn && sip->choice == SEQID_LOCAL) {
    oip = (ObjectIdPtr) sip->data.ptrvalue;
    if (oip != NULL && (! StringHasNoText (oip->str))) {
      MemSet ((Pointer) &vn, 0, sizeof (SeqId));
      MemSet ((Pointer) &tsi, 0, sizeof (TextSeqId));
      tsi.accession = oip->str;
      vn.choice = SEQID_GENBANK;
      vn.data.ptrvalue = (Pointer) &tsi;
      bsp = BioseqFind (&vn);
    }
  }
  if (bsp == NULL) return;

  if (findWorst) {
    bestid = SeqIdFindWorstForPromotion (bsp->id);
  } else {
    bestid = SeqIdFindBestForPromotion (bsp->id);
  }
  if (bestid == NULL) return;
  newid = SeqIdDup (bestid);
  if (newid == NULL) return;

  oldid = ValNodeNew (NULL);
  if (oldid == NULL) return;

  MemCopy (oldid, sip, sizeof (ValNode));
  oldid->next = NULL;

  sip->choice = newid->choice;
  sip->data.ptrvalue = newid->data.ptrvalue;

  SeqIdFree (oldid);
  ValNodeFree (newid);

  SeqIdStripLocus (sip);
}

static void PromoteSeqIdList (SeqIdPtr sip, Boolean alsoCheckLocalAccn, Boolean findWorst)

{
  while (sip != NULL) {
    PromoteSeqId (sip, alsoCheckLocalAccn, findWorst);
    sip = sip->next;
  }
}

static void PromoteSeqLocList (SeqLocPtr slp, Boolean alsoCheckLocalAccn, Boolean findWorst)

{
  SeqLocPtr      loc;
  PackSeqPntPtr  psp;
  SeqBondPtr     sbp;
  SeqIntPtr      sinp;
  SeqIdPtr       sip;
  SeqPntPtr      spp;

  while (slp != NULL) {
    switch (slp->choice) {
      case SEQLOC_NULL :
        break;
      case SEQLOC_EMPTY :
      case SEQLOC_WHOLE :
        sip = (SeqIdPtr) slp->data.ptrvalue;
        PromoteSeqIdList (sip, alsoCheckLocalAccn, findWorst);
        break;
      case SEQLOC_INT :
        sinp = (SeqIntPtr) slp->data.ptrvalue;
        if (sinp != NULL) {
          sip = sinp->id;
          PromoteSeqIdList (sip, alsoCheckLocalAccn, findWorst);
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL) {
          sip = spp->id;
          PromoteSeqIdList (sip, alsoCheckLocalAccn, findWorst);
        }
        break;
      case SEQLOC_PACKED_PNT :
        psp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (psp != NULL) {
          sip = psp->id;
          PromoteSeqIdList (sip, alsoCheckLocalAccn, findWorst);
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        loc = (SeqLocPtr) slp->data.ptrvalue;
        while (loc != NULL) {
          PromoteSeqLocList (loc, alsoCheckLocalAccn, findWorst);
          loc = loc->next;
        }
        break;
      case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          spp = (SeqPntPtr) sbp->a;
          if (spp != NULL) {
            sip = spp->id;
            PromoteSeqIdList (sip, alsoCheckLocalAccn, findWorst);
          }
          spp = (SeqPntPtr) sbp->b;
          if (spp != NULL) {
            sip = spp->id;
            PromoteSeqIdList (sip, alsoCheckLocalAccn, findWorst);
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
}

static Boolean PromoteIDsProc (GatherObjectPtr gop, Boolean findWorst)

{
  CodeBreakPtr  cbp;
  CdRegionPtr   crp;
  RnaRefPtr     rrp;
  SeqFeatPtr    sfp;
  tRNAPtr       trp;

  if (gop->itemtype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gop->dataptr;
  if (sfp == NULL) return TRUE;

  PromoteSeqLocList (sfp->location, FALSE, findWorst);

  PromoteSeqLocList (sfp->product, FALSE, findWorst);

  switch (sfp->data.choice) {
    case SEQFEAT_CDREGION :
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL && crp->code_break != NULL) {
        for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
          PromoteSeqLocList (cbp->loc, FALSE, findWorst);
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->type == 3 && rrp->ext.choice == 2) {
        trp = rrp->ext.value.ptrvalue;
        if (trp != NULL && trp->anticodon != NULL) {
          PromoteSeqLocList (trp->anticodon, FALSE, findWorst);
        }
      }
      break;
    default :
      break;
  }

  return TRUE;
}

static Boolean PromoteBestIDsProc (GatherObjectPtr gop)

{
  return PromoteIDsProc (gop, FALSE);
}

void PromoteToBestIDProc (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  oldscope, sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  oldscope = SeqEntrySetScope (sep);

  GatherObjectsInEntity (bfp->input_entityID, 0, NULL, PromoteBestIDsProc, NULL, NULL);

  SeqEntrySetScope (oldscope);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static Boolean PromoteWorstIDsProc (GatherObjectPtr gop)

{
  return PromoteIDsProc (gop, TRUE);
}

void PromoteToWorstIDProc (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  oldscope, sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  oldscope = SeqEntrySetScope (sep);

  GatherObjectsInEntity (bfp->input_entityID, 0, NULL, PromoteWorstIDsProc, NULL, NULL);

  SeqEntrySetScope (oldscope);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void ChangeGBNameIDs (SeqIdPtr sip, Pointer userdata)

{
  Char          buf [80];
  ObjectIdPtr   oip;
  TextSeqIdPtr  tsip;

  if (sip == NULL || sip->choice != SEQID_GENBANK) return;
  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
  if (tsip == NULL) return;
  if (tsip->accession != NULL || StringHasNoText (tsip->name)) return;
  StringNCpy (buf, tsip->name, sizeof (buf));
  TrimSpacesAroundString (buf);
  oip = ObjectIdNew ();
  if (oip == NULL) return;
  oip->str = StringSave (buf);
  sip->choice = SEQID_LOCAL;
  sip->data.ptrvalue = (Pointer) oip;
  TextSeqIdFree (tsip);
}

static void ChangeGBNameBioseq (BioseqPtr bsp, Pointer userdata)

{
  VisitSeqIdsInBioseq (bsp, userdata, ChangeGBNameIDs);
  SeqMgrReplaceInBioseqIndex (bsp);
}

static void ChangeGBNameFeature (SeqFeatPtr sfp, Pointer userdata)

{
  VisitSeqIdsInSeqFeat (sfp, userdata, ChangeGBNameIDs);
}

void ChangeGenBankNameToLocal (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  oldscope, sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  oldscope = SeqEntrySetScope (sep);

  VisitBioseqsInSep (sep, NULL, ChangeGBNameBioseq);
  VisitFeaturesInSep (sep, NULL, ChangeGBNameFeature);

  SeqEntrySetScope (oldscope);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


static SeqEntryPtr GetParentForSeqAnnot (SeqAnnotPtr sap)
{
  if (sap == NULL || sap->idx.parentptr == NULL) {
    return NULL;
  }
  if (sap->idx.parenttype == OBJ_BIOSEQSET) {
    return SeqMgrGetSeqEntryForData (sap->idx.parentptr);
  } else if (sap->idx.parenttype == OBJ_BIOSEQ) {
    return SeqMgrGetSeqEntryForData (sap->idx.parentptr);
  } else {
    return NULL;
  }
}


extern void PromoteAlignIDsProc (SeqAnnotPtr sp, Pointer data)

{
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  PackSegPtr    psp;
  SeqAlignPtr   sap;
  SeqLocPtr     slp;
  StdSegPtr     ssp;
  SeqEntryPtr   orig_scope, new_scope = NULL;

  if (sp == NULL || sp->type != 2) {
    return;
  }
  sap = (SeqAlignPtr) sp->data;
  if (sap == NULL) return;

  new_scope = GetParentForSeqAnnot(sp);
 
  orig_scope = SeqEntrySetScope (new_scope);

  switch (sap->segtype) {
    case SAS_DENDIAG :
      for (ddp = sap->segs; ddp != NULL; ddp = ddp->next) {
        PromoteSeqIdList (ddp->id, TRUE, FALSE);
      }
      break;
    case SAS_DENSEG :
      dsp = sap->segs;
      if (dsp != NULL) {
        PromoteSeqIdList (dsp->ids, TRUE, FALSE);
      }
      break;
    case SAS_STD :
      for (ssp = sap->segs; ssp != NULL; ssp = ssp->next) {
        PromoteSeqIdList (ssp->ids, TRUE, FALSE);
        for (slp = ssp->loc; slp != NULL; slp = slp->next) {
          PromoteSeqLocList (slp, TRUE, FALSE);
        }
      }
      break;
    case SAS_PACKED :
      psp = (PackSegPtr) sap->segs;
      if (psp != NULL) {
        PromoteSeqIdList (psp->ids, TRUE, FALSE);
      }
      break;
    default :
      break;
  }

  SeqEntrySetScope (orig_scope);
}


extern void PromoteAlignsToBestIDProc (IteM i);
void PromoteAlignsToBestIDProc (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  oldscope, sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  oldscope = SeqEntrySetScope (sep);
  VisitAnnotsInSep (sep, NULL, PromoteAlignIDsProc);

  SeqEntrySetScope (oldscope);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static Boolean RemoveGBIDsProc (GatherObjectPtr gop)

{
  BioseqPtr      bsp;
  SeqIdPtr       nextid, sip;
  SeqIdPtr PNTR  previd;

  if (gop->itemtype != OBJ_BIOSEQ) return TRUE;
  bsp = (BioseqPtr) gop->dataptr;
  if (bsp == NULL) return TRUE;

  previd = (SeqIdPtr PNTR) &(bsp->id);
  sip = bsp->id;
  while (sip != NULL) {
    nextid = sip->next;
    if (sip->choice == SEQID_GENBANK) {
      *previd = sip->next;
      sip->next = NULL;
      SeqIdFree (sip);
    } else {
      previd = (SeqIdPtr PNTR) &(sip->next);
    }
    sip = sip->next;
  }

  return TRUE;
}

void RemoveGBIDsFromBioseqs (IteM i)

{
  MsgAnswer    ans;
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  ans = Message (MSG_OKC, "Currently deletes all GenBank SeqIDs");
  if (ans == ANS_CANCEL) return;

  GatherObjectsInEntity (bfp->input_entityID, 0, NULL, RemoveGBIDsProc, NULL, NULL);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static Boolean RemoveGBProtIDsProc (GatherObjectPtr gop)

{
  BioseqPtr      bsp;
  SeqIdPtr       nextid, sip;
  SeqIdPtr PNTR  previd;

  if (gop->itemtype != OBJ_BIOSEQ) return TRUE;
  bsp = (BioseqPtr) gop->dataptr;
  if (bsp == NULL) return TRUE;
  if (! ISA_aa (bsp->mol)) return TRUE;

  previd = (SeqIdPtr PNTR) &(bsp->id);
  sip = bsp->id;
  while (sip != NULL) {
    nextid = sip->next;
    if (sip->choice == SEQID_GENBANK) {
      *previd = sip->next;
      sip->next = NULL;
      SeqIdFree (sip);
    } else {
      previd = (SeqIdPtr PNTR) &(sip->next);
    }
    sip = sip->next;
  }

  return TRUE;
}

void RemoveGBIDsFromProteins (IteM i)

{
  MsgAnswer    ans;
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  ans = Message (MSG_OKC, "Currently deletes all GenBank SeqIDs");
  if (ans == ANS_CANCEL) return;

  GatherObjectsInEntity (bfp->input_entityID, 0, NULL, RemoveGBProtIDsProc, NULL, NULL);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static Boolean RemoveGIProc (GatherObjectPtr gop)

{
  BioseqPtr      bsp;
  SeqIdPtr       nextid, sip;
  SeqIdPtr PNTR  previd;

  if (gop->itemtype != OBJ_BIOSEQ) return TRUE;
  bsp = (BioseqPtr) gop->dataptr;
  if (bsp == NULL) return TRUE;

  previd = (SeqIdPtr PNTR) &(bsp->id);
  sip = bsp->id;
  while (sip != NULL) {
    nextid = sip->next;
    if (sip->choice == SEQID_GI) {
      *previd = sip->next;
      sip->next = NULL;
      SeqIdFree (sip);
    } else {
      previd = (SeqIdPtr PNTR) &(sip->next);
    }
    sip = sip->next;
  }

  return TRUE;
}

void RemoveGIsFromBioseqs (IteM i)

{
  MsgAnswer    ans;
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  ans = Message (MSG_OKC, "Currently deletes all GI SeqIDs");
  if (ans == ANS_CANCEL) return;

  GatherObjectsInEntity (bfp->input_entityID, 0, NULL, RemoveGIProc, NULL, NULL);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


typedef struct blastfields {
  BlastNet3Hptr        bl3hp;
  BLAST_OptionsBlkPtr  options;
} BlastFields, PNTR BlastFieldsPtr;

#define EXPECT_VALUE 0.01

static void BlastCDD (BioseqPtr bsp, Pointer userdata)

{
  BlastFieldsPtr       bfp;
  BlastNet3Hptr        bl3hp;
  BLAST_OptionsBlkPtr  options;
  ValNodePtr           error_returns = NULL;
  SeqAlignPtr          salp;

  if (! ISA_aa (bsp->mol)) return;
  bfp = (BlastFieldsPtr) userdata;
  bl3hp = bfp->bl3hp;
  options = bfp->options;
  if (bl3hp == NULL) return;

  /* do blast search */

  salp = BlastBioseqNet (bl3hp, bsp, "blastp", "cdd", options,
                         NULL, &error_returns, NULL);

  /* BlastErrorPrintExtra (error_returns, TRUE, stdout); */

  /* annotation function now moved to rpsutil in toolkit for common use */

  AnnotateRegionsFromCDD (bsp, salp, EXPECT_VALUE);

  /* clean up */

  SeqAlignSetFree (salp);
}

extern void SimpleCDDBlastProc (IteM i)

{
  BaseFormPtr          bfp;
  SeqEntryPtr          sep;
  BlastFields          bf;
  BlastNet3Hptr        bl3hp = NULL;
  BLAST_OptionsBlkPtr  options = NULL;
  BlastResponsePtr     response = NULL;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  TransientSetAppParam ("NCBI", "NETBLAST", "SERVICE_NAME", "rpsblast");

  if (! BlastInit ("Sequin", &bl3hp, &response)) return;

  response = BlastResponseFree (response); 
  options = BLASTOptionNew ("blastp", TRUE);
  if (options == NULL) return;

  options->is_rps_blast = TRUE;
  options->filter = FILTER_SEG;
  options->expect_value  = (Nlm_FloatHi) EXPECT_VALUE;
  options->hitlist_size = 5000;

  /* blast fetch enable needed to retrieve by general SeqID */

  BlastNetBioseqFetchEnable (bl3hp, "cdd", FALSE, TRUE);

  bf.bl3hp = bl3hp;
  bf.options = options;

  WatchCursor ();
  Update ();

  FreeCDDRegions (sep);

  VisitBioseqsInSep (sep, (Pointer) &bf, BlastCDD);

  RemoveDuplicateCDDs (sep);

  BlastFini (bl3hp);
  options = BLASTOptionDelete (options);
  BlastNetBioseqFetchDisable (bl3hp, "cdd", FALSE);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
}
/*#endif*/


extern Boolean DoBuildContig (void)

{
  MsgAnswer    ans;
  BioseqPtr    bsp;
  Pointer      dataptr;
  Uint2        datatype;
  Uint2        entityID;
  FILE         *fp;
  Int2         handled;
  Char         path [PATH_MAX];
  SeqEntryPtr  sep;

  if (! GetInputFileName (path, sizeof (path), NULL, "TEXT")) return FALSE;
  fp = FileOpen (path, "r");
  if (fp == NULL) return FALSE;
  ans = Message (MSG_YN, "Are coordinates on the master (as opposed to the individual accession)?");
  sep = ReadContigList (fp, (Boolean) (ans == ANS_YES));
  FileClose (fp);
  if (sep == NULL || sep->choice != 1) return FALSE;

  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return FALSE;
  datatype = OBJ_BIOSEQ;
  dataptr = (Pointer) bsp;
  entityID = ObjMgrRegister (datatype, dataptr);
  if (dataptr == NULL || entityID == 0) return FALSE;

  WatchCursor ();
  Update ();
  handled = GatherProcLaunch (OMPROC_VIEW, FALSE, entityID, 1,
                              OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
  ArrowCursor ();
  Update ();

  if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
    Message (MSG_FATAL, "Unable to launch viewer.");
    SeqEntryFree (sep);
    return FALSE;
  }

  ObjMgrSetOptions (OM_OPT_FREE_IF_NO_VIEW, entityID);
  ObjMgrSetDirtyFlag (entityID, TRUE);
  return TRUE;
}

static void DoParseTrinomial (BioSourcePtr biop, Pointer userdata)

{
  BinomialOrgNamePtr  bonp;
  OrgModPtr           omp;
  OrgNamePtr          onp;
  OrgRefPtr           orp;
  CharPtr             str;

  if (biop == NULL) return;
  orp = biop->org;
  if (orp == NULL) return;
  onp = orp->orgname;
  if (onp == NULL) return;
  if (onp->choice != 1) return;
  bonp = (BinomialOrgNamePtr) onp->data;
  if (bonp == NULL) return;
  if (StringHasNoText (bonp->subspecies)) return;
  for (omp = onp->mod; omp != NULL; omp = omp->next) {
    if (omp->subtype == ORGMOD_sub_species) return;
  }
  str = bonp->subspecies;
  if (StringNICmp (str, "subsp. ", 7) == 0) {
    str += 7;
    if (StringHasNoText (str)) return;
  }
  omp = OrgModNew ();
  if (omp == NULL) return;
  omp->subtype = ORGMOD_sub_species;
  omp->subname = StringSave (str);
  omp->next = onp->mod;
  onp->mod = omp;
}

extern void ParseTrinomial (IteM i);
extern void ParseTrinomial (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioSourcesInSep (sep, NULL, DoParseTrinomial);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}



#define NUM_SITES 27
static CharPtr siteString [NUM_SITES] = {
  "", "active site", "binding site", "cleavage site", "inhibition site", "modified site",
  "glycosylation site", "myristoylation site", "mutagenized site", "metal binding site",
  "phosphorylation site", "acetylation site", "amidation site", "methylation site",
  "hydroxylation site", "sulfatation site", "oxidative deamination site",
  "pyrrolidone carboxylic acid site", "gamma carboxyglutamic acid site",
  "blocked site", "lipid binding site", "np binding site", "DNA binding site",
  "signal peptide", "transit peptide", "transmembrane region", "nitrosylation"
};

/*=========================================================================*/
/*                                                                         */
/* SeqLocAdjustByOffset ()                                                 */
/*                                                                         */
/*=========================================================================*/

extern void SeqLocAdjustByOffset (SeqLocPtr slp,
                  Int4 offset)
{
  SeqIntPtr sinp;

  switch (slp->choice) {
  case SEQLOC_INT :
    sinp = (SeqIntPtr) slp->data.ptrvalue;
    if (NULL == sinp)
      break;
    sinp->from += offset;
    sinp->to   += offset;
    break;
  case SEQLOC_EMPTY :
  case SEQLOC_NULL :
  case SEQLOC_WHOLE :
  case SEQLOC_PNT :
  case SEQLOC_PACKED_PNT :
  case SEQLOC_PACKED_INT :
  case SEQLOC_MIX :
  case SEQLOC_EQUIV :
  case SEQLOC_BOND :
  case SEQLOC_FEAT :
  default :
    break;
  }
}

      
static Boolean ConvertImpToSpecialRNA 
(SeqFeatPtr sfp,
 Uint2      featdef_to,
 Pointer    extradata)
{
  RnaRefPtr          rrp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_IMP)
  {
    return FALSE;
  }
  rrp = RnaRefNew ();
  if (rrp != NULL) {
    sfp->data.value.ptrvalue = ImpFeatFree ((ImpFeatPtr) sfp->data.value.ptrvalue);
    sfp->data.choice = SEQFEAT_RNA;
    sfp->data.value.ptrvalue = (Pointer) rrp;
    if (featdef_to == FEATDEF_precursor_RNA) {
      rrp->type = 1;
    } else {
      rrp->type = 255;
    }
  }
  return TRUE;
}

static Boolean ConvertRegionToRNA 
(SeqFeatPtr sfp,
 Uint2      featdef_to,
 Pointer    extradata)
{
  RnaRefPtr  rrp;
  CharPtr    str;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_REGION)
  {
    return FALSE;
  }

  rrp = RnaRefNew ();
  if (NULL == rrp)
  {
    return FALSE;
  }

  str = (CharPtr) sfp->data.value.ptrvalue;
  sfp->data.choice = SEQFEAT_RNA;
  sfp->data.value.ptrvalue = (Pointer) rrp;

  if (featdef_to == FEATDEF_precursor_RNA) {
    rrp->type = 1;
  } else {
    rrp->type = 255;
  }

  if (! StringHasNoText (str)) {
    rrp->ext.choice = 1;
    rrp->ext.value.ptrvalue = str;
  }
  return TRUE;
}


static Boolean CommentShouldBeProduct (CharPtr comment)
{
  CharPtr search1 = "internal transcribed spacer ";
  CharPtr search2 = "ITS";
  CharPtr found;
  Int4 len;

  if (StringHasNoText (comment)) {
    return FALSE;
  }
  if ((found = StringISearch (comment, search1)) != NULL
      && (len = StringLen (search1)) < StringLen (found)
      && (*(found + len) == '1' || *(found + len) == '2' || *(found + len) == '3')) {
    return TRUE;
  }

  if ((found = StringISearch (comment, search2)) != NULL
      && (len = StringLen (search2)) < StringLen (found)
      && (*(found + len) == '1' || *(found + len) == '2' || *(found + len) == '3')) {
    return TRUE;
  }
  return FALSE;
}


static Boolean ConvertImpToRNA 
(SeqFeatPtr sfp,
 Uint2      featdef_to,
 Pointer    extradata)
{
  RnaRefPtr  rrp;
  RNAGenPtr  rgp;
  Boolean    was_misc_feat = FALSE;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_IMP)
  {
    return FALSE;
  }
  if (sfp->idx.subtype == FEATDEF_misc_feature) {
    was_misc_feat = TRUE;
  }

  rrp = RnaRefNew ();
  if (NULL == rrp)
    return FALSE;

  sfp->data.value.ptrvalue =
    ImpFeatFree ((ImpFeatPtr) sfp->data.value.ptrvalue);
  sfp->data.choice = SEQFEAT_RNA;
  sfp->data.value.ptrvalue = (Pointer) rrp;

  switch (featdef_to) {
  case FEATDEF_preRNA :
    rrp->type = 1;
    break;
  case FEATDEF_mRNA :
    rrp->type = 2;
    break;
  case FEATDEF_tRNA :
    rrp->type = 3;
    break;
  case FEATDEF_rRNA :
    rrp->type = 4;
    break;
  case FEATDEF_snRNA :
    rrp->type = 5;
    break;
  case FEATDEF_scRNA :
    rrp->type = 6;
    break;
  case FEATDEF_snoRNA :
    rrp->type = 7;
    break;
  case FEATDEF_ncRNA :
    rrp->type = 8;
    rrp->ext.choice = 3;
    rrp->ext.value.ptrvalue = RNAGenNew ();
    break;
  case FEATDEF_tmRNA :
    rrp->type = 9;
    rrp->ext.choice = 3;
    rrp->ext.value.ptrvalue = RNAGenNew ();
    break;
  case FEATDEF_misc_RNA:
    rrp->type = 10;
    rrp->ext.choice = 3;
    rgp = RNAGenNew ();
    rrp->ext.value.ptrvalue = rgp;
    if (was_misc_feat && CommentShouldBeProduct(sfp->comment)) {
      rgp->product = sfp->comment;
      sfp->comment = NULL;
    }
    break;
  case FEATDEF_otherRNA :
    if (was_misc_feat && CommentShouldBeProduct(sfp->comment)) {
      rrp->type = 10;
      rrp->ext.choice = 3;
      rgp = RNAGenNew ();
      rrp->ext.value.ptrvalue = rgp;
      rgp->product = sfp->comment;
      sfp->comment = NULL;
    } else {      
      rrp->type = 255;
      rrp->ext.choice = 1;
      rrp->ext.value.ptrvalue = StringSave ("misc_RNA");
    }
    break;
  default :
    break;
  }
  return TRUE;
}

static Boolean ConvertCommentToMiscFeat (SeqFeatPtr sfp, Uint2 featdef_to, Pointer extradata)
{
  ImpFeatPtr ifp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_COMMENT || sfp->data.value.ptrvalue != NULL)
  {
    return FALSE;
  }
  
  ifp = ImpFeatNew ();
  if (ifp != NULL) {
    ifp->key = StringSave ("misc_feature");
    sfp->data.choice = SEQFEAT_IMP;
    sfp->data.value.ptrvalue = (Pointer) ifp;
    return TRUE;
  }
  return FALSE;
}

static Boolean ConvertGeneToMiscFeat 
(SeqFeatPtr sfp,
 Uint2      featdef_to,
 Pointer    extradata)
{
  return ConvertGeneToMiscFeatFunc (sfp, featdef_to);
}

static Boolean ConvertRNAToMiscFeat
(SeqFeatPtr sfp,
 Uint2      featdef_to,
 Pointer    extradata)
{
  return ConvertRNAToImpFeat (sfp, featdef_to);
}

static Boolean ConvertSiteToMiscFeat
(SeqFeatPtr sfp,
 Uint2      featdef_to,
 Pointer    extradata)
{
  GBQualPtr  gbqual;
  ImpFeatPtr ifp;
  Int2       sitetype;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_SITE)
  {
    return FALSE;
  }

  ifp = ImpFeatNew ();
  if (NULL == ifp)
  {
    return FALSE;
  }

  sitetype = (Int2) sfp->data.value.intvalue;
  sfp->data.choice = SEQFEAT_IMP;
  sfp->data.value.ptrvalue = (Pointer) ifp;
  ifp->key = StringSave ("misc_feature");
  if (sitetype > 0 && sitetype < NUM_SITES) {
    gbqual = GBQualNew ();
    if (gbqual != NULL) {
      gbqual->qual = StringSave ("note");
      gbqual->val = StringSave (siteString [sitetype]);
      gbqual->next = sfp->qual;
      sfp->qual = gbqual;
    }
  }
  return TRUE;
}

static Boolean ConvertProtToRegion 
(SeqFeatPtr sfp,
 Uint2      featdef_to,
 Pointer    extradata)
{
  ProtRefPtr prp;
  ValNodePtr vnp;
  CharPtr    str;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT)
  {
    return FALSE;
  }
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (NULL == prp)
  {
    return FALSE;
  }

  vnp = prp->name;
  if (vnp != NULL && vnp->next == NULL) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (! StringHasNoText (str)) {
      vnp->data.ptrvalue = NULL;
      sfp->data.value.ptrvalue = ProtRefFree (prp);
      sfp->data.choice = SEQFEAT_REGION;
      sfp->data.value.ptrvalue = (Pointer) str;
    }
  }
  return TRUE;
}


static Boolean ConvertRegionToProt
(SeqFeatPtr sfp,
 Uint2      featdef_to,
 Pointer    extradata)
{
  return ConvertRegionToProtFunc (sfp, featdef_to);
}


static Boolean ConvertRNAToNcRNA
(SeqFeatPtr sfp)
{
  RnaRefPtr  rrp;
  tRNAPtr    trp;
  RNAGenPtr  rgp = NULL;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return FALSE;
  
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (NULL == rrp) {
    rrp = RnaRefNew ();
    sfp->data.value.ptrvalue = rrp;
  }

  if (rrp->ext.choice == 2) {
    trp = (tRNAPtr) rrp->ext.value.ptrvalue;
    if (trp != NULL) {
      trp->anticodon = SeqLocFree (trp->anticodon);
      trp = MemFree (trp);
      rrp->ext.value.ptrvalue = trp;
    }
    rrp->ext.choice = 0;
  } 

  /* might be converting from other ncRNA, tmRNA, or misc_RNA */
  if (rrp->ext.choice == 1) {
    rgp = RNAGenNew ();
    rgp->product = rrp->ext.value.ptrvalue;
    rrp->ext.choice = 3;
    rrp->ext.value.ptrvalue = rgp;
  } else if (rrp->ext.choice == 0) {
    rgp = RNAGenNew ();
    rrp->ext.choice = 3;
    rrp->ext.value.ptrvalue = rgp;
  } else if (rrp->ext.choice == 3) {
    rgp = rrp->ext.value.ptrvalue;
    if (rgp == NULL) {
      rgp = RNAGenNew ();
      rrp->ext.value.ptrvalue = rgp;
    }
  }

  if (rgp != NULL) {
    if (rrp->type == 5) {
      rgp->_class = StringSave ("snRNA");
    } else if (rrp->type == 6) {
      rgp->_class = StringSave ("scRNA");
    } else if (rrp->type == 7) {
      rgp->_class = StringSave ("snoRNA");
    }
  }

  rrp->type = 8;
  return TRUE;
}


static void ConvertFromncRNAtoMiscRNA (SeqFeatPtr sfp)
{
  RnaRefPtr rrp;
  RNAGenPtr rgp;
  CharPtr tmp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return;
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL || rrp->ext.choice != 3) {
    return;
  }

  rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
  if (rgp == NULL) {
    return;
  }

  /* if no product but has comment, promote comment to product */
  if (rgp->product == NULL && !StringHasNoText (sfp->comment)) {
    rgp->product = sfp->comment;
    sfp->comment = NULL;
  }
  
  /* move ncRNA_class elsewhere */
  if (!StringHasNoText (rgp->_class)) {
    if (StringCmp (rgp->_class, "other") == 0) {
      rgp->_class = MemFree (rgp->_class);
    } else if (rgp->product == NULL) {
      rgp->product = rgp->_class;
      rgp->_class = NULL;
    } else {
      if (StringHasNoText (sfp->comment)) {
        sfp->comment = MemFree (sfp->comment);
        sfp->comment = StringSave (rgp->_class);
      } else {
        tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (sfp->comment) + StringLen (rgp->_class) + 3));
        sprintf (tmp, "%s; %s", sfp->comment, rgp->_class);
        sfp->comment = MemFree (sfp->comment);
        sfp->comment = tmp;
      }
      rgp->_class = MemFree (rgp->_class);
    }
  }    
}


static Boolean ConvertRNAToRNA 
(SeqFeatPtr sfp,
 Uint2      featdef_to,
 Pointer    extradata)
{
  RnaRefPtr  rrp;
  Boolean    from_misc = FALSE, to_misc = FALSE;
  GBQualPtr  gbq, gbq_p = NULL;

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (NULL == rrp)
    return FALSE;

  if (rrp->type == 255) {
    from_misc = TRUE;
  }

  switch (featdef_to) {
    case FEATDEF_preRNA :
      rrp->type = 1;
      break;
    case FEATDEF_mRNA :
      rrp->type = 2;
      break;
    case FEATDEF_tRNA :
      rrp->type = 3;
      break;
    case FEATDEF_rRNA :
      rrp->type = 4;
      break;
    case FEATDEF_snRNA :
      rrp->type = 5;
      to_misc = TRUE;
      break;
    case FEATDEF_scRNA :
      rrp->type = 6;
      to_misc = TRUE;
      break;
    case FEATDEF_snoRNA :
      rrp->type = 7;
      to_misc = TRUE;
      break;
    case FEATDEF_ncRNA :
      return ConvertRNAToNcRNA (sfp);
      break;
    case FEATDEF_otherRNA :
      rrp->type = 255;
      ConvertFromncRNAtoMiscRNA (sfp);
      break;
    default:
      break;
  }
  if (from_misc && !to_misc && rrp->ext.choice == 1 && StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") == 0) {
    /* move product qual to string */
    gbq = sfp->qual;
    while (gbq != NULL && StringCmp (gbq->qual, "product") != 0) {
      gbq_p = gbq;
      gbq = gbq->next;
    }
    if (gbq == NULL) {
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
      rrp->ext.choice = 0;
    } else {
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
      rrp->ext.value.ptrvalue = gbq->val;
      gbq->val = NULL;
      if (gbq_p == NULL) {
        sfp->qual = gbq->next;
      } else {
        gbq_p->next = gbq->next;
      }
      gbq->next = NULL;
      gbq = GBQualFree (gbq);
    }
  }
  return TRUE;
}

static Boolean ConvertncRNAToMiscBinding
(SeqFeatPtr sfp,
 Uint2      featdef_to,
 Pointer    extradata)
{
  RnaRefPtr  rrp;
  RNAGenPtr  rgp;
  Boolean    asked = TRUE, replace = FALSE, use_semicolon = TRUE;
  ImpFeatPtr ifp;

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (NULL == rrp)
    return FALSE;

  if (rrp->ext.choice == 1) {
    /* move product to note */
    AppendOrReplaceString (&(sfp->comment), rrp->ext.value.ptrvalue, &asked, &replace, &use_semicolon);
  } else if (rrp->ext.choice == 3 && (rgp = (RNAGenPtr) rrp->ext.value.ptrvalue) != NULL
             && !StringHasNoText (rgp->product)) {
    AppendOrReplaceString (&(sfp->comment), rgp->product, &asked, &replace, &use_semicolon);
  }
  rrp = RnaRefFree (rrp);
  sfp->data.choice = SEQFEAT_IMP;
  ifp = ImpFeatNew ();
  ifp->key = StringSave ("misc_binding");
  sfp->data.value.ptrvalue = ifp;
  sfp->idx.subtype = 0;

  return TRUE;
}

static Boolean ConvertProtToProt 
(SeqFeatPtr sfp,
 Uint2      featdef_to,
 Pointer    extradata)
{
  return ConvertProtToProtFunc (sfp, featdef_to);
}


static Boolean ConvertImpToProt (SeqFeatPtr sfp, Uint2 featdef_to, Pointer extradata)
{
  return ConvertImpToProtFunc (sfp, featdef_to);
}


static Boolean ConvertProtToImp (SeqFeatPtr sfp, Uint2 featdef_to, Pointer extradata)
{
  return ConvertProtToImpFunc (sfp, featdef_to);
}


static Boolean ConvertBioSrcToRepeatRegionEx (SeqFeatPtr sfp, Uint2 featdef_to, Pointer extradata)
{
  return ConvertBioSrcToRepeatRegion (sfp, featdef_to);
}


static Boolean ConvertCDSToMatPeptide (SeqFeatPtr sfp, Uint2 featdef_to, Pointer extradata)
{
  return AutoConvertCDSToMiscFeat (sfp, extradata == NULL ? TRUE : !(*((BoolPtr) extradata)));
}


static Boolean ConvertMiscFeatureToCDSFunction (SeqFeatPtr sfp, Uint2 featdef_to, Pointer extradata)
{
  return ConvertMiscFeatToCodingRegion (sfp);
}


static Boolean ConvertMiscFeatureToGeneFunction (SeqFeatPtr sfp, Uint2 featdef_to, Pointer extradata)
{
  return ConvertMiscFeatToGene (sfp);
}


extern EnumFieldAssoc  enum_bond_alist [];
extern EnumFieldAssoc  enum_site_alist [];


static GbFeatName EditQualifierList[] = {
 {"allele", Class_text}, 
 {"anticodon", Class_pos_aa},
 {"bound_moiety", Class_text},
 {"chromosome", Class_text},
 {"citation", Class_bracket_int},
 {"codon", Class_seq_aa},
 {"codon_start", Class_int_or}, 
 {"compare", Class_text },
 {"cons_splice", Class_site},
 {"db_xref", Class_text},
 {"direction", Class_L_R_B}, 
 {"EC_number", Class_ecnum},
 {"evidence", Class_exper}, 
 {"exception", Class_text},
 {"experiment", Class_text},
 {"frequency", Class_text}, 
 {"function", Class_text},
 {"gene", Class_text}, 
 {"gdb_xref", Class_text},
 {"label", Class_token},
 {"map", Class_text},
 {"mobile_element", Class_text}, 
 {"mod_base", Class_token}, {"note", Class_note},
 {"number", Class_number}, 
#ifdef INTERNAL_NCBI_SEQUIN
 {"old_locus_tag", Class_text},
#endif
 {"operon", Class_text},
 {"organism", Class_text},
 {"partial", Class_none}, {"PCR_conditions", Class_text},
 {"phenotype", Class_text},
 {"plasmid", Class_text}, {"product", Class_text},
 {"pseudo", Class_none},
 {"rearranged", Class_none}, { "replace", Class_text},
 {"rpt_family", Class_text}, {"rpt_type", Class_rpt},
 {"rpt_unit", Class_token}, {"rpt_unit_seq", Class_token}, {"rpt_unit_range", Class_token},
 {"satellite", Class_text},
 {"sequenced_mol", Class_text},
 {"standard_name", Class_text},
 {"translation", Class_text}, {"transl_except", Class_pos_aa},
 {"transl_table", Class_int},
 {"usedin", Class_token},
 {"focus", Class_none},
 {"protein_id", Class_text},
 {"organelle", Class_text}, {"transcript_id", Class_text},
 {"transgenic", Class_none}, {"environmental_sample", Class_none},
 {"locus_tag", Class_text}, {"mol_type", Class_text},
 {"segment", Class_text}
};

const Int4 NumEditQualifiers = sizeof (EditQualifierList) / sizeof (GbFeatName);
#define QUAL_EVIDENCE 13
#define QUAL_EXPERIMENT 15


static ENUM_ALIST(biosource_genome_alistX)
  {" ",                    0},
  {"Apicoplast",          16},
  {"Chloroplast",          2},
  {"Chromatophore",       22},
  {"Chromoplast",          3},
  {"Chromosome",          21},
  {"Cyanelle",            12},
  {"Endogenous-virus",    19},
  {"Extrachromosomal",     8},
  {"Genomic",              1},
  {"Hydrogenosome",       20},
  {"Insertion Sequence",  11},
  {"Kinetoplast",          4},
  {"Leucoplast",          17},
  {"Macronuclear",         7},
  {"Mitochondrion",        5},
  {"Nucleomorph",         15},
  {"Plasmid",              9},
  {"Plastid",              6},
  {"Proplastid",          18},
  {"Proviral",            13},
  {"Transposon",          10},
END_ENUM_ALIST

static ENUM_ALIST(orgref_textfield_alist)
  {" ",                    0},
  {"Scientific Name",      1},
  {"Common Name",          2},
  {"Lineage",              3},
  {"Division",             4},
END_ENUM_ALIST

typedef struct sourceformdata {
  FEATURE_FORM_BLOCK

  Int2           type;
  ButtoN         applyToParts;
  TexT           onlyThisPart;
  GrouP          sourceGroup;
  GrouP          modGrp;
  GrouP          genGrp;
  GrouP          refGrp;
  GrouP          txtGrp;
  GrouP          originGrp;
  PopuP          frommod;
  PopuP          tomod;
  PopuP          fromgen;
  PopuP          togen;
  PopuP          fromref;
  PopuP          toref;
  PopuP          toorigin;
  PopuP          fromorigin;
  TexT           findthis;
  TexT           replacewith;
  GrouP          applyChoiceGrp;  /* This is the group for the radio buttons for the
                                   * constraint on which biosources to operate on.
                                   */
  PopuP          qual_to_have;    /* This is the control that indicates which qualifier
                                   * should be present to operate on a biosource.
                                   */
  TexT           text_to_have;    /* This is the control that indicates the text that
                                   * should be present in the biosource to be operated on.
                                   */
  Int4           applyChoiceVal;  /* This is the value of the applyChoiceGrp choice. */
  CharPtr        text_to_have_str;/* This is the contents of the text_to_have box. */
  Int4           qual_to_have_val;/* This is the subtype for the qual to have. */

  Int2           choice;
  Int2           fromval;
  Int2           toval;
  Int2           onlythis;
  CharPtr        findStr;
  CharPtr        replaceStr;

  Boolean        replaceOldAsked;
  Boolean        doReplaceAll;
  Boolean        use_semicolon;
  ButtoN         leaveDlgUp;
} SourceFormData, PNTR SourceFormPtr;

Uint2 mod_widths [] = {
  0, 25
};

Uint2 mod_types [] = {
  TAGLIST_POPUP, TAGLIST_TEXT
};


extern Uint1 FindTypeForModNameText (CharPtr cp)
{
  Uint1 subtype;

  subtype = EquivalentSubSourceEx (cp, TRUE);
  if (subtype == 0) {
    subtype = EquivalentOrgMod (cp);
    if (subtype > 200) {
      /* function that calls this expects less than 100 vals for orgmod */
      /* need to adjust for old-lineage and old-name */
      subtype -= 200;
    }
  } else {
    /* function that calls this uses >100 to decide if subsource */
    subtype += 100;
  }
  if (subtype == 0) {
    subtype = SUBSRC_other;
  }
  return subtype;
}


extern void AppendOrReplaceString (
  CharPtr PNTR string_loc,
  CharPtr new_value,
  Boolean PNTR asked_question,
  Boolean PNTR do_replace,
  Boolean PNTR use_semicolon
)
{
  MsgAnswer ans;
  CharPtr   tmp_value, tmp_new_value;

  if (string_loc == NULL
    || new_value == NULL
    || asked_question == NULL
    || do_replace == NULL
    || use_semicolon == NULL)
  {
    return;
  }

  if (! *asked_question && !StringHasNoText (*string_loc))
  {
    *asked_question = TRUE;
    ArrowCursor ();
    ans = Message (MSG_YN, "Do you wish to overwrite existing content?");
    *do_replace = (Boolean) (ans == ANS_YES);
    if (! *use_semicolon)
    {
      if (! *do_replace)
      {
        ans = Message (MSG_YN, "Separate items with semicolon?");
        *use_semicolon = (Boolean) (ans == ANS_YES);
      }
    }
    WatchCursor ();
    Update ();
  }
  if (*do_replace || StringHasNoText (*string_loc))
  {
    MemFree (*string_loc);
    *string_loc = StringSave (new_value);
  }
  else
  {
    tmp_value = MemNew (StringLen (*string_loc) + StringLen ( new_value) + 3);
    if (tmp_value == NULL) return;
    StringCpy (tmp_value, *string_loc);
    TrimSpacesAroundString (tmp_value);
    if (*use_semicolon)
    {
      StringCat (tmp_value, "; ");
    }
    else
    {
      StringCat (tmp_value, " ");
    }
    tmp_new_value = StringSave (new_value);
    TrimSpacesAroundString (tmp_new_value);
    StringCat (tmp_value, tmp_new_value);
    MemFree (tmp_new_value);
    MemFree (*string_loc);
    *string_loc = tmp_value;
  }
}

static ENUM_ALIST(origin_alist)
{" ",               ORG_DEFAULT    },
{"Natural",         ORG_NATURAL    },
{"Natural Mutant",  ORG_NATMUT     },
{"Mutant",          ORG_MUT        },
{"Artificial",      ORG_ARTIFICIAL },
{"Synthetic",       ORG_SYNTHETIC  },
{"Other",           ORG_OTHER      },
END_ENUM_ALIST


typedef struct setsample
{
  GetFeatureFieldString    fieldstring_func;
  GetDescriptorFieldString descrstring_func;
  Uint2                    entityID;
  ValNodePtr               field_list;
  FreeValNodeProc          free_vn_proc;
  CopyValNodeDataProc      copy_vn_proc;
  MatchValNodeProc         match_vn_proc;
  NameFromValNodeProc      label_vn_proc;  
  
  FilterSetPtr             fsp;
} SetSampleData, PNTR SetSamplePtr;

static ValNodePtr IntValNodeCopy (ValNodePtr vnp)
{
  ValNodePtr vnp_new;
  
  if (vnp == NULL)
  {
    return NULL;
  }
  vnp_new = ValNodeNew (NULL);
  if (vnp_new != NULL)
  {
    vnp_new->choice = vnp->choice;
    vnp_new->data.intvalue = vnp->data.intvalue;
  }
  return vnp_new;
}

static Boolean IntValNodeMatch (ValNodePtr vnp1, ValNodePtr vnp2)
{
  if (vnp1 == NULL && vnp2 == NULL)
  {
    return TRUE;
  }
  else if (vnp1 == NULL || vnp2 == NULL)
  {
    return FALSE;
  }
  else if (vnp1->data.intvalue == vnp2->data.intvalue)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

static void LocationConstraintXClearText (LocationConstraintXPtr lcp)
{
  if (lcp != NULL)
  {
    lcp->left = -1;
    lcp->right = -1;
  }
}

static LocationConstraintXPtr LocationConstraintXFree (LocationConstraintXPtr lcp)
{
  lcp = MemFree (lcp);
  return lcp;
}

static LocationConstraintXPtr LocationConstraintXCopy (LocationConstraintXPtr lcp)
{
  LocationConstraintXPtr lcp_new;
  if (lcp == NULL)
  {
    return NULL;
  }
  lcp_new = (LocationConstraintXPtr) MemNew (sizeof (LocationConstraintXData));
  if (lcp_new != NULL)
  {
    lcp_new->left = lcp->left;
    lcp_new->right = lcp->right;
    lcp_new->interval_end_choice = lcp->interval_end_choice;
    lcp_new->match_choice = lcp->match_choice;
    lcp_new->strand = lcp->strand;
    lcp_new->sequence_type = lcp->sequence_type;
  }
  return lcp_new;
}

extern StringConstraintXPtr StringConstraintXFree (StringConstraintXPtr scp)
{
  if (scp == NULL) return NULL;
  scp->match_text = MemFree (scp->match_text);
  scp = MemFree (scp);
  return scp;  
}

static void StringConstraintClearText (StringConstraintXPtr scp)
{
  if (scp != NULL)
  {
    scp->match_text = MemFree (scp->match_text);
  }
}

static StringConstraintXPtr StringConstraintCopy (StringConstraintXPtr scp)
{
  StringConstraintXPtr scp_new;
  if (scp == NULL)
  {
    return NULL;
  }
  
  scp_new = (StringConstraintXPtr) MemNew (sizeof (StringConstraintData));
  if (scp_new != NULL)
  {
    scp_new->match_text = StringSave (scp->match_text);
    scp_new->match_location = scp->match_location;
    scp_new->insensitive = scp->insensitive;
    scp_new->whole_word = scp->whole_word;
    scp_new->not_present = scp->not_present;   
  }
  return scp_new;
}

extern ValNodePtr ValNodeFuncFree (ValNodePtr vnp, FreeValNodeProc free_vn_proc)
{
  if (vnp == NULL)
  {
    return NULL;
  }
  vnp->next = ValNodeFuncFree (vnp->next, free_vn_proc);
  if (free_vn_proc != NULL)
  {
    free_vn_proc (vnp);
  }
  ValNodeFree (vnp);
  return NULL;
}

extern ChoiceConstraintPtr ChoiceConstraintFree (ChoiceConstraintPtr scp)
{
  if (scp == NULL) return NULL;
  scp->qual_choice = ValNodeFuncFree (scp->qual_choice, scp->free_vn_proc);
  scp->qual_choice_match = ValNodeFuncFree (scp->qual_choice, scp->free_vn_proc);
  scp->string_constraint = StringConstraintXFree (scp->string_constraint);
  scp->pseudo_constraint = MemFree (scp->pseudo_constraint);
  scp = MemFree (scp);
  return scp;
}

static ValNodePtr QualChoiceCopy (ValNodePtr vnp, CopyValNodeDataProc copy_vn_proc)
{
  ValNodePtr new_list = NULL;
  
  if (vnp == NULL || copy_vn_proc == NULL)
  {
    return NULL;
  }
  
  new_list = copy_vn_proc (vnp);
  new_list->next = QualChoiceCopy (vnp->next, copy_vn_proc);
  return new_list;
}

static ChoiceConstraintPtr ChoiceConstraintCopy (ChoiceConstraintPtr ccp)
{
  ChoiceConstraintPtr ccp_new;
  if (ccp == NULL || ccp->copy_vn_proc == NULL)
  {
    return NULL;
  }
  ccp_new = (ChoiceConstraintPtr) MemNew (sizeof (ChoiceConstraintData));
  if (ccp_new != NULL)
  {
    ccp_new->constraint_type = ccp->constraint_type;
    ccp_new->qual_choice = QualChoiceCopy (ccp->qual_choice, ccp->copy_vn_proc);
    ccp_new->qual_choice_match = QualChoiceCopy (ccp->qual_choice_match, ccp->copy_vn_proc);
    ccp_new->string_constraint = StringConstraintCopy (ccp->string_constraint);
    ccp_new->copy_vn_proc = ccp->copy_vn_proc;
    ccp_new->free_vn_proc = ccp->free_vn_proc;
  }
  return ccp_new;
}

static void ChoiceConstraintClearText (ChoiceConstraintPtr scp)
{
  if (scp != NULL)
  {
    StringConstraintClearText (scp->string_constraint);
  }
}


extern void FilterSetClearText (FilterSetPtr fsp)
{
  if (fsp != NULL)
  {
    StringConstraintClearText (fsp->scp);
    ChoiceConstraintClearText (fsp->ccp);
    LocationConstraintXClearText (fsp->lcp);
    StringConstraintClearText (fsp->id_list);
  }
}

extern FilterSetPtr FilterSetNew (void)
{
  FilterSetPtr fsp_new;
  fsp_new = (FilterSetPtr) MemNew (sizeof (FilterSetData));
  if (fsp_new != NULL)
  {
    fsp_new->scp = NULL;
    fsp_new->ccp = NULL;
    fsp_new->lcp = NULL;
    fsp_new->cgp = NULL;
    fsp_new->id_list = NULL;
  }
  return fsp_new;
}

extern FilterSetPtr FilterSetFree (FilterSetPtr fsp)
{
  if (fsp != NULL)
  {
    fsp->scp = StringConstraintXFree (fsp->scp);
    fsp->ccp = ChoiceConstraintFree (fsp->ccp);
    fsp->lcp = LocationConstraintXFree (fsp->lcp);
    fsp->cgp = ChoiceConstraintFree (fsp->cgp);
    fsp->id_list = StringConstraintXFree (fsp->id_list);
    fsp = MemFree (fsp);
  }
  return fsp;
}

static FilterSetPtr FilterSetCopy (FilterSetPtr fsp)
{
  FilterSetPtr fsp_new;
  
  if (fsp == NULL)
  {
    return NULL;
  }
  fsp_new = (FilterSetPtr) MemNew (sizeof (FilterSetData));
  if (fsp_new != NULL)
  {
    fsp_new->scp = StringConstraintCopy (fsp->scp);
    fsp_new->ccp = ChoiceConstraintCopy (fsp->ccp);
    fsp_new->lcp = LocationConstraintXCopy (fsp->lcp);
    fsp_new->cgp = ChoiceConstraintCopy (fsp->cgp);
    fsp_new->id_list = StringConstraintCopy (fsp->id_list);
  }
  return fsp_new;
}


static SetSamplePtr SetSampleFree (SetSamplePtr ssp)
{
  ValNodePtr vnp;
  
  if (ssp == NULL)
  {
    return NULL;
  }
  if (ssp->free_vn_proc != NULL)
  {
    for (vnp = ssp->field_list; vnp != NULL; vnp = vnp->next)
    {
      (ssp->free_vn_proc)(vnp);
    }
  }
  ssp->field_list = ValNodeFree (ssp->field_list);
   
  ssp = MemFree (ssp);
  return ssp;
}

static SetSamplePtr SetSampleCopy (SetSamplePtr ssp)
{
  SetSamplePtr ssp_new;
  ValNodePtr   vnp_old, vnp_new, vnp_last = NULL;
  
  if (ssp == NULL)
  {
    return NULL;
  }
  ssp_new = (SetSamplePtr) MemNew (sizeof (SetSampleData));
  if (ssp_new != NULL)
  {
    ssp_new->fieldstring_func = ssp->fieldstring_func;
    ssp_new->descrstring_func = ssp->descrstring_func;
    ssp_new->free_vn_proc = ssp->free_vn_proc;
    ssp_new->copy_vn_proc = ssp->copy_vn_proc;
    ssp_new->match_vn_proc = ssp->match_vn_proc;
    ssp_new->label_vn_proc = ssp->label_vn_proc;  
    
    ssp_new->entityID = ssp->entityID;
    ssp_new->field_list = NULL;
    for (vnp_old = ssp->field_list; vnp_old != NULL; vnp_old = vnp_old->next)
    {
      vnp_new = (ssp_new->copy_vn_proc)(vnp_old);
      if (vnp_last == NULL)
      {
        ssp_new->field_list = vnp_new;
      }
      else
      {
        vnp_last->next = vnp_new;
      }
      vnp_last = vnp_new;
    }
    ssp_new->fsp = FilterSetCopy (ssp->fsp);
  }
  return ssp_new;
}

#define PARSE_FIELD_DEFLINE       1
#define PARSE_FIELD_BIOSRC_STRING 2
#define PARSE_FIELD_SOURCE_QUAL   3
#define PARSE_FIELD_GENE_FIELD    4
#define PARSE_FIELD_RNA_FIELD    5
#define PARSE_FIELD_CDS_COMMENT   6
#define PARSE_FIELD_PROTEIN_FIELD 7
#define PARSE_FIELD_IMPORT_QUAL   8
#define PARSE_FIELD_FEATURE_NOTE  9
#define PARSE_FIELD_COMMENT_DESC  10
#define PARSE_FIELD_DBXREF        11

#define MAX_PARSE_FIELD_TYPE   11
#define SEARCH_FIELD_PUBLICATION  12

#define PARSE_FIELD_FIRST_FEATURE 4
#define PARSE_FIELD_LAST_FEATURE  9

extern Boolean DoesLocationMatchConstraint (SeqLocPtr slp, LocationConstraintXPtr lcp);

typedef struct acceptpolicy 
{
  Boolean leave_dlg_up;
} AcceptPolicyData, PNTR AcceptPolicyPtr;

typedef struct acceptcanceldialog 
{
  DIALOG_MESSAGE_BLOCK
  Nlm_AcceptActnProc    accept_actn;
  Nlm_CancelActnProc    cancel_actn;
  Nlm_ClearActnProc     clear_actn;
  Nlm_ClearTextActnProc clear_text_actn;
  WindoW                w;
  
  ButtoN                accept_btn;
  
  ButtoN                leave_dlg_up_btn;
} AcceptCancelDialogData, PNTR AcceptCancelDialogPtr;

static void AcceptCancelDialogAccept (ButtoN b)
{
  AcceptCancelDialogPtr acdp;

  acdp = (AcceptCancelDialogPtr) GetObjectExtra (b);
  if (acdp == NULL) return;
  
  if (acdp->accept_actn != NULL)
  {
    if (!((acdp->accept_actn) (acdp->userdata)))
    {
      return;
    }
  }
  if (! GetStatus (acdp->leave_dlg_up_btn))
  {
    Remove (acdp->w);
  }
}

static void AcceptCancelDialogCancel (ButtoN b)
{
  AcceptCancelDialogPtr acdp;

  acdp = (AcceptCancelDialogPtr) GetObjectExtra (b);
  if (acdp == NULL) return;
  
  if (acdp->cancel_actn != NULL)
  {
    (acdp->cancel_actn) (acdp->userdata);    
  }
  Remove (acdp->w);
}

static void AcceptCancelDialogClear (ButtoN b)
{
  AcceptCancelDialogPtr acdp;

  acdp = (AcceptCancelDialogPtr) GetObjectExtra (b);
  if (acdp == NULL || acdp->clear_actn == NULL) return;
  
  (acdp->clear_actn) (acdp->userdata);
}

static void AcceptCancelDialogClearText (ButtoN b)
{
  AcceptCancelDialogPtr acdp;

  acdp = (AcceptCancelDialogPtr) GetObjectExtra (b);
  if (acdp == NULL || acdp->clear_text_actn == NULL) return;
  
  (acdp->clear_text_actn) (acdp->userdata);
}

static void AcceptPolicyToAcceptCancelDialog (DialoG d, Pointer data)

{
  AcceptCancelDialogPtr   acdp;
  AcceptPolicyPtr         app;

  acdp = (AcceptCancelDialogPtr) GetObjectExtra (d);
  app = (AcceptPolicyPtr) data;
  if (acdp == NULL || app == NULL)
  {
    return;
  }

  SetStatus (acdp->leave_dlg_up_btn, app->leave_dlg_up);  
}

static Pointer AcceptCancelDialogToAcceptPolicy (DialoG d)
{
  AcceptCancelDialogPtr acdp;
  AcceptPolicyPtr       app;
  
  acdp = (AcceptCancelDialogPtr) GetObjectExtra (d);
  if (acdp == NULL) 
  {
    return NULL;
  }
  
  app = (AcceptPolicyPtr) MemNew (sizeof (AcceptPolicyData));
  if (app != NULL)
  {
    app->leave_dlg_up = GetStatus (acdp->leave_dlg_up_btn);
  }
  return app;
}

static void AcceptCancelMessage (DialoG d, Int2 mssg)

{
  AcceptCancelDialogPtr acdp;

  acdp = (AcceptCancelDialogPtr) GetObjectExtra (d);
  if (acdp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        /* reset accept policy */
        SetStatus (acdp->leave_dlg_up_btn, FALSE);
        break;
      case VIB_MSG_ENTER :
        Select (acdp->accept_btn);
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestAcceptCancelDialog (DialoG d)

{
  return NULL;
}

extern DialoG AcceptCancelDialog 
(GrouP                 parent,
 Nlm_AcceptActnProc    accept_actn,
 Nlm_CancelActnProc    cancel_actn,
 Nlm_ClearActnProc     clear_actn,
 Nlm_ClearTextActnProc clear_text_actn,
 Pointer               userdata,
 WindoW                w)
{
  AcceptCancelDialogPtr acdp;
  GrouP                 grp, policy_grp = NULL, other_grp;
  ButtoN                b;
  
  acdp = (AcceptCancelDialogPtr) MemNew (sizeof (AcceptCancelDialogData));
  grp = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (grp, acdp, StdCleanupExtraProc);
  SetGroupSpacing (grp, 10, 10);

  acdp->dialog = (DialoG) grp;
  acdp->todialog = AcceptPolicyToAcceptCancelDialog;
  acdp->fromdialog = AcceptCancelDialogToAcceptPolicy;
  acdp->dialogmessage = AcceptCancelMessage;
  acdp->testdialog = TestAcceptCancelDialog;

  acdp->accept_actn             = accept_actn;
  acdp->cancel_actn             = cancel_actn;
  acdp->clear_actn              = clear_actn;
  acdp->clear_text_actn         = clear_text_actn;
  acdp->userdata                = userdata;
  acdp->w                       = w;

  policy_grp = HiddenGroup (grp, 3, 0, NULL);
  SetGroupSpacing (policy_grp, 10, 10);
  
  b = PushButton (policy_grp, "Clear", AcceptCancelDialogClear);
  SetObjectExtra (b, acdp, NULL);
 
  b = PushButton (policy_grp, "Clear Text", AcceptCancelDialogClearText);
  SetObjectExtra (b, acdp, NULL);
  
  acdp->leave_dlg_up_btn = CheckBox (policy_grp, "Leave Dialog Up", NULL);
  
  other_grp = HiddenGroup (grp, 5, 0, NULL);
  SetGroupSpacing (other_grp, 10, 10);
  acdp->accept_btn = PushButton (other_grp, "Accept", AcceptCancelDialogAccept);
  SetObjectExtra (acdp->accept_btn, acdp, NULL);
  
  b = PushButton (other_grp, "Cancel", AcceptCancelDialogCancel);
  SetObjectExtra (b, acdp, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) policy_grp, (HANDLE) other_grp, NULL);
  
  return (DialoG) grp;
}


extern void EnableAcceptCancelDialogAccept (DialoG d)
{
  AcceptCancelDialogPtr acdp;

  acdp = (AcceptCancelDialogPtr) GetObjectExtra (d);

  if (acdp == NULL) return;
  
  InvalObject (acdp->accept_btn);
  SafeEnable (acdp->accept_btn);
  Update ();
}

extern void DisableAcceptCancelDialogAccept (DialoG d)
{
  AcceptCancelDialogPtr acdp;

  acdp = (AcceptCancelDialogPtr) GetObjectExtra (d);

  if (acdp == NULL) return;
  
  InvalObject (acdp->accept_btn);
  SafeDisable (acdp->accept_btn);
  Update ();
}





static ValNodePtr ValNodeAppend (ValNodePtr list, ValNodePtr second_list)
{
  ValNodePtr last_vnp;
  
  if (list == NULL)
  {
    list = second_list;
  }
  else if (second_list != NULL)
  {
    last_vnp = list;
    while (last_vnp != NULL && last_vnp->next != NULL)
    {
      last_vnp = last_vnp->next;
    }
    last_vnp->next = second_list;
  }
  return list;
}

static void FindFeatureTypeCallback (SeqFeatPtr sfp, Pointer userdata)
{
  BoolPtr type_list;
  
  if (sfp == NULL || userdata == NULL) return;
  
  type_list = (BoolPtr) userdata;
  type_list[sfp->idx.subtype] = TRUE;
}


static void AddFeatureToList (ValNodePtr PNTR list, FeatDefPtr curr)
{
  Char        str [256];

  if (list == NULL || curr == NULL) return;
  
  if (curr->featdef_key == FEATDEF_PUB)
  {
    StringNCpy_0 (str, curr->typelabel, sizeof (str) - 15);
    StringCat (str, " (Publication)");
    ValNodeAddPointer (list, curr->featdef_key, StringSave (str));
  }
  else if (curr->featdef_key == FEATDEF_BIOSRC)
  {
    ValNodeAddPointer (list, FEATDEF_BIOSRC, StringSave ("BioSource"));
  }
  else if (curr->featdef_key == FEATDEF_otherRNA)
  {
    ValNodeAddPointer (list, FEATDEF_otherRNA, StringSave ("misc_RNA"));
  }
  else if (curr->featdef_key == FEATDEF_misc_RNA ||
           curr->featdef_key == FEATDEF_precursor_RNA ||
           curr->featdef_key == FEATDEF_mat_peptide ||
           curr->featdef_key == FEATDEF_sig_peptide ||
           curr->featdef_key == FEATDEF_transit_peptide ||
           curr->featdef_key == FEATDEF_Imp_CDS)
  {
    StringNCpy_0 (str, curr->typelabel, sizeof (str) - 10);
    StringCat (str, "_imp");
    ValNodeAddPointer (list, curr->featdef_key, StringSave (str));
  }
  else
  {
    ValNodeAddPointer (list, curr->featdef_key, StringSave (curr->typelabel));
  }
}

extern ValNodePtr BuildFeatureDialogList (Boolean list_most_used_first, SeqEntryPtr sep)
{
  ValNodePtr  feature_choice_list = NULL;
  ValNodePtr  import_feature_list = NULL;
  ValNodePtr  least_liked_feature_list = NULL;
  ValNodePtr  most_used_feature_list = NULL;
  ValNodePtr  unused_list = NULL;
  FeatDefPtr  curr;
  Uint1       key;
  CharPtr     label = NULL;
  Boolean     feature_type_list[256];
  Boolean     only_most_used = TRUE;
  Int4        i;

  if (sep == NULL) {
      for (i = 0; i < 255; i++) {
          feature_type_list[i] = TRUE;
      }
  } else {
      for (i = 0; i < 255; i++) {
          feature_type_list[i] = FALSE;
      }
      VisitFeaturesInSep (sep, feature_type_list, FindFeatureTypeCallback);
  }

  curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
  while (curr != NULL) {
    if (key == FEATDEF_BAD
        || key == FEATDEF_Imp_CDS
        || key == FEATDEF_source
        || key == FEATDEF_ORG
        || IsUnwantedFeatureType (key))
    {
      curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
      continue;
    }
    else if (!feature_type_list[key]) 
    {
      AddFeatureToList (&unused_list, curr);
      curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
      continue;
    }
    
    if (curr->featdef_key == FEATDEF_otherRNA) {
      AddFeatureToList (&feature_choice_list, curr);
      if (list_most_used_first) {
        AddFeatureToList (&most_used_feature_list, curr);
      }
    }
    else if (curr->featdef_key == FEATDEF_misc_RNA ||
            curr->featdef_key == FEATDEF_precursor_RNA ||
            curr->featdef_key == FEATDEF_mat_peptide ||
            curr->featdef_key == FEATDEF_sig_peptide ||
            curr->featdef_key == FEATDEF_transit_peptide ||
            curr->featdef_key == FEATDEF_Imp_CDS)
    {
      AddFeatureToList (&import_feature_list, curr);
    }
    else if (curr->featdef_key == FEATDEF_TXINIT)
    {
      AddFeatureToList (&least_liked_feature_list, curr);
    }
    else
    {
      ValNodeAddPointer (&feature_choice_list, curr->featdef_key, StringSave (curr->typelabel));
    }
    
    if (list_most_used_first)
    {
      /* also build most used list */
      if (curr->featdef_key == FEATDEF_CDS
          || curr->featdef_key == FEATDEF_exon
          || curr->featdef_key == FEATDEF_GENE
          || curr->featdef_key == FEATDEF_intron
          || curr->featdef_key == FEATDEF_mRNA
          || curr->featdef_key == FEATDEF_rRNA
          || curr->featdef_key == FEATDEF_PROT)
      {
        ValNodeAddPointer (&most_used_feature_list, curr->featdef_key, 
                           StringSave (curr->typelabel));
      }
      else
      {
        only_most_used = FALSE;
      }
    }
    
    curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
  }

  most_used_feature_list = SortValNode (most_used_feature_list,
                                        CompareFeatureValNodeStrings);
  feature_choice_list = SortValNode (feature_choice_list,
                                     CompareFeatureValNodeStrings);
  import_feature_list = SortValNode (import_feature_list,
                                     CompareFeatureValNodeStrings);
  least_liked_feature_list = SortValNode (least_liked_feature_list,
                                          CompareFeatureValNodeStrings);
  unused_list = SortValNode (unused_list, CompareFeatureValNodeStrings);       

  if (most_used_feature_list != NULL && only_most_used) {
    feature_choice_list = ValNodeFreeData (feature_choice_list);
    feature_choice_list = most_used_feature_list;
    least_liked_feature_list = ValNodeFreeData(least_liked_feature_list);
    import_feature_list = ValNodeFreeData(least_liked_feature_list);
  } else {
    feature_choice_list = ValNodeAppend (most_used_feature_list, feature_choice_list);
    feature_choice_list = ValNodeAppend (feature_choice_list, least_liked_feature_list);
    feature_choice_list = ValNodeAppend (feature_choice_list, import_feature_list);
  }
  feature_choice_list = ValNodeAppend (feature_choice_list, unused_list);
  
  return feature_choice_list;
}

static ValNodePtr BuildImportList (void)
{
  ValNodePtr  import_feature_list = NULL;
  FeatDefPtr  curr;
  Uint1       key;
  CharPtr     label = NULL;
  Char        str [256];

  curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
  while (curr != NULL) {
    if (FindFeatFromFeatDefType (key) == SEQFEAT_IMP && !IsUnwantedFeatureType(key))
    {
      ValNodeAddPointer (&import_feature_list, curr->featdef_key, StringSave (str));
    }
    
    curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
  } 
  return import_feature_list; 
}

static void RemoveFeatureSelectionDuplicates (ValNodePtr orig_list)
{
  ValNodePtr vnp1, vnp2, prev_vnp, next_vnp;
  
  if (orig_list == NULL || orig_list->next == NULL)
  {
    return;
  }
  
  vnp1 = orig_list;
  while (vnp1 != NULL && vnp1->next != NULL)
  {
    vnp2 = vnp1->next;
    prev_vnp = vnp1;
    while (vnp2 != NULL)
    {
      next_vnp = vnp2->next;
      if (vnp2->choice == vnp1->choice)
      {
        prev_vnp->next = next_vnp;
        vnp2->next = NULL;
        vnp2 = ValNodeFreeData (vnp2);
      }
      else
      {
        prev_vnp = vnp2;
      }
      vnp2 = next_vnp;
    }
    vnp1 = vnp1->next;
  }
}

static ValNodePtr FeatureSelectionRemap (ValNodePtr orig_list)
{
  ValNodePtr vnp, vnp_next, prev_vnp, repl_vnp = NULL;
  
  if (orig_list == NULL)
  {
    return NULL;
  }
  
  vnp = orig_list;
  prev_vnp = NULL;
  while (vnp != NULL && repl_vnp == NULL)
  {
    vnp_next = vnp->next;
    /* replace FEATDEF_IMP with list of import features */
    if (vnp->choice == FEATDEF_IMP)
    {
      /* build replacement list */
      repl_vnp = BuildImportList ();
      
      ValNodeLink (&repl_vnp, vnp->next);
      
      if (prev_vnp == NULL)
      {
        orig_list = repl_vnp;
      }
      else
      {
        prev_vnp->next = repl_vnp;
      }

      vnp->next = NULL;
      ValNodeFreeData (vnp);
    }
    else
    {
      prev_vnp = vnp;
    }
    vnp = vnp_next;
  }
  
  RemoveFeatureSelectionDuplicates (orig_list);
  
  return orig_list; 
}

extern DialoG 
FeatureSelectionDialogEx 
(GrouP                    h,
 Boolean                  allow_multi,
 SeqEntryPtr              sep,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  DialoG      dlg;
  ValNodePtr  feature_choice_list, vnp, vnp_next, vnp_prev;
  Boolean     found_import = FALSE;

  feature_choice_list = BuildFeatureDialogList (TRUE, sep);
  
  if (!allow_multi)
  {
    /* remove Import option */
    vnp_prev = NULL;
    for (vnp = feature_choice_list; vnp != NULL && !found_import; vnp = vnp_next)
    {
      vnp_next = vnp->next;
      if (vnp->choice == FEATDEF_IMP)
      {
        if (vnp_prev == NULL)
        {
          feature_choice_list = vnp_next;
        }
        else
        {
          vnp_prev->next = vnp_next;
        }
        vnp->next = NULL;
        vnp = ValNodeFreeData (vnp);
        found_import = TRUE;
      }
      else
      {
        vnp_prev = vnp;
      }
    }
  }
  
  /* note - the ValNodeSelectionDialog will free the feature_choice_list when done */
  
  dlg = ValNodeSelectionDialogEx (h, feature_choice_list, TALL_SELECTION_LIST, 
                                ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "feature list", 
                                change_notify, change_userdata, allow_multi, FALSE,
                                allow_multi ? FeatureSelectionRemap : NULL);
  return dlg;
}

extern DialoG 
FeatureSelectionDialog 
(GrouP                    h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  return FeatureSelectionDialogEx(h, allow_multi, NULL, change_notify, change_userdata);
}

extern DialoG
DescriptorSelectionDialog
(GrouP                    h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  ValNodePtr descriptor_list = BuildDescriptorValNodeList ();
  return ValNodeSelectionDialogEx (h, descriptor_list, TALL_SELECTION_LIST,
                                ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "descriptor list", 
                                change_notify, change_userdata, allow_multi, FALSE,
                                NULL);
}


extern CharPtr SourceQualValNodeName (ValNodePtr vnp)
{
  SourceQualDescPtr           sqdp;
  
  if (vnp == NULL || vnp->data.ptrvalue == NULL)
  {
    return NULL;
  }
  else if (vnp->choice == 0)
  {
    sqdp = (SourceQualDescPtr) vnp->data.ptrvalue;
    return StringSave (sqdp->name);
  }
  else
  {
    return StringSave (vnp->data.ptrvalue);
  }
}

extern ValNodePtr SourceQualValNodeDataCopy (ValNodePtr vnp)
{
  SourceQualDescPtr           sqdp, new_sqdp = NULL;
  ValNodePtr                  vnp_copy = NULL;
  
  if (vnp != NULL && vnp->data.ptrvalue != NULL)
  {
    if (vnp->choice == 0)
    {
      sqdp = (SourceQualDescPtr) vnp->data.ptrvalue;
      new_sqdp = (SourceQualDescPtr) MemNew (sizeof (SourceQualDescData));
      if (new_sqdp != NULL)
      {
        MemCpy (new_sqdp, sqdp, sizeof (SourceQualDescData));
      }
      ValNodeAddPointer (&vnp_copy, vnp->choice, new_sqdp);
    }
    else
    {
      ValNodeAddPointer (&vnp_copy, vnp->choice, StringSave (vnp->data.ptrvalue));
    }
  }
  return vnp_copy;
}


NLM_EXTERN ValNodePtr LIBCALL SourceQualValNodeCopy (ValNodePtr vnp)
{
  return SourceQualValNodeDataCopy(vnp);
}


extern Boolean SourceQualValNodeMatch (ValNodePtr vnp1, ValNodePtr vnp2)
{
  SourceQualDescPtr  sqdp1, sqdp2;
  Boolean            rval = FALSE;
  
  if (vnp1 == NULL || vnp2 == NULL
      || vnp1->data.ptrvalue == NULL || vnp2->data.ptrvalue == NULL)
  {
    rval = FALSE;
  }
  else if (vnp1->choice != vnp2->choice)
  {
    rval = FALSE;
  }
  else if (vnp1->choice == 0)
  {
    sqdp1 = (SourceQualDescPtr) vnp1->data.ptrvalue;
    sqdp2 = (SourceQualDescPtr) vnp2->data.ptrvalue;
  
    if (sqdp1->subtype == sqdp2->subtype
        && ((sqdp1->isOrgMod && sqdp2->isOrgMod)
            || (!sqdp1->isOrgMod && ! sqdp2->isOrgMod)))
    {
      rval = TRUE;
    }
    else
    {
      rval = FALSE;
    }
  }
  else
  {
    if (StringCmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue) == 0)
    {
      rval = TRUE;
    }
    else
    {
      rval = FALSE;
    }
  }
  return rval;
}

static DialoG SourceQualTypeSelectionDialogEx
(GrouP h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Int2                     list_height,
 Boolean                  show_discouraged,
 Boolean                  show_discontinued)

{
  DialoG     dlg;
  ValNodePtr qual_choice_list;

  if (allow_multi)
  {
    qual_choice_list = ValNodeNew (NULL);
    qual_choice_list->choice = 1;
    qual_choice_list->data.ptrvalue = StringSave ("All Notes");
    qual_choice_list->next = GetSourceQualDescList (TRUE, TRUE, show_discouraged, show_discontinued);
  }
  else
  {
    qual_choice_list = GetSourceQualDescList (TRUE, TRUE, show_discouraged, show_discontinued);
  }
  /* note - the ValNodeSelectionDialog will free the qual_choice_list when done */                                            
  dlg = ValNodeSelectionDialog (h, qual_choice_list, list_height, SourceQualValNodeName,
                                ValNodeSimpleDataFree, SourceQualValNodeDataCopy,
                                SourceQualValNodeMatch, "feature list", 
                                change_notify, change_userdata, allow_multi);

  return dlg;
}

extern DialoG SourceQualTypeSelectionDialog 
(GrouP h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
 
{
  return SourceQualTypeSelectionDialogEx (h, allow_multi, change_notify,
                                          change_userdata, 6, FALSE, FALSE);
}

static DialoG SourceQualTypeDiscSelectionDialog 
(GrouP h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
 
{
  return SourceQualTypeSelectionDialogEx (h, allow_multi, change_notify,
                                          change_userdata, 6, TRUE, TRUE);
}

static DialoG SourceQualTypeConstraintSelectionDialog
(GrouP h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)

{
  return SourceQualTypeSelectionDialogEx (h, allow_multi, change_notify,
                                          change_userdata, 4, TRUE, FALSE);
}

static DialoG SourceStringConstraintSelectionDialog
(GrouP h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)

{
  DialoG     dlg;
  ValNodePtr qual_choice_list = NULL;

  /* ignore allow_multi */
  
  ValNodeAddPointer (&qual_choice_list, 1, StringSave ("Organism or Any Qual"));
  ValNodeAddPointer (&qual_choice_list, 1, StringSave ("Organism"));
  qual_choice_list->next->next = GetSourceQualDescList (TRUE, TRUE, TRUE, FALSE);
  ValNodeAddPointer (&qual_choice_list, 1, StringSave ("Lineage"));
  ValNodeAddPointer (&qual_choice_list, 1, StringSave ("Location"));

  /* note - the ValNodeSelectionDialog will free the qual_choice_list when done */                                            
  dlg = ValNodeSelectionDialog (h, qual_choice_list, 4, SourceQualValNodeName,
                                ValNodeSimpleDataFree, SourceQualValNodeDataCopy,
                                SourceQualValNodeMatch, "feature list", 
                                change_notify, change_userdata, FALSE);

  return dlg;
}

static Boolean OrgModSpecialMatch (OrgModPtr mod, FilterSetPtr fsp)
{
  SourceQualDescPtr sqdp;
  
  if (mod == NULL)
  {
    return FALSE;
  }
  if (fsp == NULL
      || fsp->ccp == NULL 
      || fsp->ccp->constraint_type != CHOICE_CONSTRAINT_STRING
      || fsp->ccp->string_constraint == NULL
      || fsp->ccp->qual_choice == NULL
      || fsp->ccp->qual_choice->data.ptrvalue == NULL)
  {
    return TRUE;
  }
  
  sqdp = (SourceQualDescPtr) fsp->ccp->qual_choice->data.ptrvalue;
  if (sqdp->subtype != mod->subtype)
  {
    return TRUE;
  }
  
  if (DoesStringMatchConstraintX (mod->subname, fsp->ccp->string_constraint))
  {
    if (fsp->ccp->string_constraint->not_present)
    {
      return FALSE;
    }
    else
    {
      return TRUE;
    }
  }
  else
  {
    if (fsp->ccp->string_constraint->not_present)
    {
      return TRUE;
    }
    else
    {
      return FALSE;
    }
  }
}

static Boolean SubSrcSpecialMatch (SubSourcePtr ssp, FilterSetPtr fsp)
{
  SourceQualDescPtr sqdp;
  
  if (ssp == NULL)
  {
    return FALSE;
  }
  if (fsp == NULL
      || fsp->ccp == NULL 
      || fsp->ccp->constraint_type != CHOICE_CONSTRAINT_STRING
      || fsp->ccp->string_constraint == NULL
      || fsp->ccp->qual_choice == NULL
      || fsp->ccp->qual_choice->data.ptrvalue == NULL)
  {
    return TRUE;
  }
  
  sqdp = (SourceQualDescPtr) fsp->ccp->qual_choice->data.ptrvalue;
  if (sqdp->subtype != ssp->subtype)
  {
    return TRUE;
  }
  if (DoesStringMatchConstraintX (ssp->name, fsp->ccp->string_constraint))
  {
    if (fsp->ccp->string_constraint->not_present)
    {
      return FALSE;
    }
    else
    {
      return TRUE;
    }
  }
  else
  {
    if (fsp->ccp->string_constraint->not_present)
    {
      return TRUE;
    }
    else
    {
      return FALSE;
    }
  }
}


static CharPtr GetSubfield (CharPtr str, Uint1 subfield)
{
  CharPtr cp, cp2;
  Int4    num_colons = 0;
  CharPtr new_val = NULL;

  if (StringHasNoText (str)) {
    return NULL;
  }

  cp = StringChr (str, ':');
  while (cp != NULL) {
    num_colons ++;
    cp = StringChr (cp + 1, ':');
  }

  if (subfield == 0) {
    new_val = StringSave (str);
  } else if (subfield == 1) {
    if (num_colons == 0) {
      return NULL;
    } else {
      cp = StringChr (str, ':');
      new_val = (CharPtr) MemNew (sizeof (Char) * (cp - str + 1));
      StringNCpy (new_val, str, cp - str);
      new_val[cp - str] = 0;
    }
  } else if (subfield == 2) {
    if (num_colons == 0 || num_colons == 1) {
      return NULL;
    } else {
      cp = StringChr (str, ':');
      cp2 = StringChr (cp + 1, ':');
      new_val = (CharPtr) MemNew (sizeof (Char) * (cp2 - cp));
      StringNCpy (new_val, cp + 1, cp2 - cp - 1);
      new_val[cp2 - cp - 1] = 0;
    }
  } else {
    if (num_colons == 0) {
      new_val = StringSave (str);
    } else {
      cp = StringRChr (str, ':');
      new_val = StringSave (cp + 1);
    }
  }
  return new_val;
}


static CharPtr 
GetSourceQualQualValueFromBiop 
(BioSourcePtr      biop,
 SourceQualDescPtr sqdp,
 FilterSetPtr      fsp)
{
  CharPtr            str = NULL, field, tmp_str;
  OrgModPtr          mod;
  SubSourcePtr       ssp;
  
  if (biop == NULL || sqdp == NULL)
  {
    return NULL;
  }

  if (sqdp->isOrgMod)
  {
    if (biop->org != NULL && biop->org->orgname != NULL)
    {
      mod = biop->org->orgname->mod;
      while (mod != NULL)
      {
        if (mod->subtype == sqdp->subtype
            && !StringHasNoText (mod->subname)
            && OrgModSpecialMatch (mod, fsp))
        {
          field = GetSubfield (mod->subname, sqdp->subfield);
          if (!StringHasNoText (field)) {
            if (str == NULL) {
              str = field;
            } else {
              tmp_str = (CharPtr) MemNew ((StringLen (str) + StringLen (field) + 3) * sizeof (Char));
              StringCpy (tmp_str, str);
              StringCat (tmp_str, "; ");
              StringCat (tmp_str, field);
              field = MemFree (field);
              str = MemFree (str);
              str = tmp_str;
            }
          } else {
            field = MemFree (field);
          }
        }
        mod = mod->next;
      }
    }
  }
  else
  {
    ssp = biop->subtype;
    while (ssp != NULL)
    {
      if (ssp->subtype == sqdp->subtype
          && !StringHasNoText (ssp->name)
          && SubSrcSpecialMatch (ssp, fsp))
      {
        field = GetSubfield (ssp->name, sqdp->subfield);
        if (!StringHasNoText (field)) {
          if (str == NULL) {
            str = field;
          } else {
            tmp_str = (CharPtr) MemNew ((StringLen (str) + StringLen (field) + 3) * sizeof (Char));
            StringCpy (tmp_str, str);
            StringCat (tmp_str, "; ");
            StringCat (tmp_str, field);
            field = MemFree (field);
            str = MemFree (str);
            str = tmp_str;
          }
        } else {
          field = MemFree (field);
        }
      }
      ssp = ssp->next;
    }
  }
  return str;  
}

static CharPtr GetLocationFromSource (BioSourcePtr biop, ValNodePtr vnp)
{
  Nlm_EnumFieldAssocPtr eap;
  
  if (biop == NULL)
  {
    return NULL;
  }
  
  for (eap = biosource_genome_alistX;
       eap->name != NULL && eap->value != biop->genome;
       eap++)
  {
  }
  
  if (eap != NULL && eap->name != NULL && !StringHasNoText (eap->name))
  {
    return StringSave (eap->name);
  }
  else
  {
    return NULL;
  }  
}

static CharPtr GetLocationFromSourceFeature (SeqFeatPtr sfp, ValNodePtr vnp, FilterSetPtr fsp)
{
  if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC)
  {
    return GetLocationFromSource (sfp->data.value.ptrvalue, vnp);
  }
  else
  {
    return NULL;
  }
}

static CharPtr 
GetLocationFromSourceDescriptor 
(SeqDescrPtr  sdp,
 ValNodePtr   vnp, 
 FilterSetPtr fsp)
{
  if (sdp != NULL && sdp->choice == Seq_descr_source)
  {
    return GetLocationFromSource (sdp->data.ptrvalue, vnp);
  }
  else
  {
    return NULL;
  }
}

static CharPtr GetOriginFromSource (BioSourcePtr biop, ValNodePtr vnp)
{
  Nlm_EnumFieldAssocPtr eap;
  
  if (biop == NULL)
  {
    return NULL;
  }
  
  for (eap = origin_alist;
       eap->name != NULL && eap->value != biop->origin;
       eap++)
  {
  }
  
  if (eap != NULL && eap->name != NULL && !StringHasNoText (eap->name))
  {
    return StringSave (eap->name);
  }
  else
  {
    return NULL;
  }  
}

static CharPtr 
GetOriginFromSourceFeature 
(SeqFeatPtr   sfp,
 ValNodePtr   vnp,
 FilterSetPtr fsp)
{
  if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC)
  {
    return GetOriginFromSource (sfp->data.value.ptrvalue, vnp);
  }
  else
  {
    return NULL;
  }
}

static CharPtr 
GetOriginFromSourceDescriptor 
(SeqDescrPtr  sdp,
 ValNodePtr   vnp,
 FilterSetPtr fsp)
{
  if (sdp != NULL && sdp->choice == Seq_descr_source)
  {
    return GetOriginFromSource (sdp->data.ptrvalue, vnp);
  }
  else
  {
    return NULL;
  }
}

static CharPtr GetStringFromSource (BioSourcePtr biop, ValNodePtr vnp)
{
  Nlm_EnumFieldAssocPtr eap;
  
  if (biop == NULL)
  {
    return NULL;
  }
  
  for (eap = orgref_textfield_alist;
       eap->name != NULL && StringCmp (vnp->data.ptrvalue, eap->name) != 0;
       eap++)
  {
  }
  
  if (eap != NULL && eap->name != NULL)
  {
    switch (eap->value)
    {
      case 1:
        if (biop->org != NULL && !StringHasNoText (biop->org->taxname))
        {
          return StringSave (biop->org->taxname);
        }
        break;
      case 2:
        if (biop->org != NULL && !StringHasNoText (biop->org->common))
        {
          return StringSave (biop->org->common);
        }
        break;
      case 3:
        if (biop->org != NULL && biop->org->orgname != NULL
            && !StringHasNoText (biop->org->orgname->lineage))
        {
          return StringSave (biop->org->orgname->lineage);
        }
        break;
      case 4:
        if (biop->org != NULL && biop->org->orgname != NULL
            && !StringHasNoText (biop->org->orgname->div))
        {
          return StringSave (biop->org->orgname->div);
        }
        break;
    }
  }
  return NULL;  
}

static CharPtr GetStringFromSourceFeature (SeqFeatPtr sfp, ValNodePtr vnp, FilterSetPtr fsp)
{
  if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC)
  {
    return GetStringFromSource (sfp->data.value.ptrvalue, vnp);
  }
  else
  {
    return NULL;
  }
}

static CharPtr StringAppend (CharPtr str1, CharPtr str2)
{
  CharPtr new_str;
  
  if (str1 == NULL && str2 == NULL)
  {
    return NULL;
  }
  else if (str1 == NULL)
  {
    return StringSave (str2);
  }
  else if (str2 == NULL)
  {
    return StringSave (str1);
  }
  else
  {
    new_str = (CharPtr) MemNew ((StringLen (str1) + StringLen (str2) + 3) * sizeof (Char));
    if (new_str != NULL)
    {
      sprintf (new_str, "%s; %s", str1, str2);  
    }
    return new_str;
  }
}

static CharPtr GetStringFromSourceDescriptor (SeqDescrPtr sdp, ValNodePtr vnp, FilterSetPtr fsp)
{
  if (sdp != NULL && sdp->choice == Seq_descr_source)
  {
    return GetStringFromSource (sdp->data.ptrvalue, vnp);
  }
  else
  {
    return NULL;
  }
}

static CharPtr 
GetSourceQualStringFromBioSource 
(BioSourcePtr biop,
 ValNodePtr   vnp,
 FilterSetPtr fsp)
{
  CharPtr            str = NULL, str1 = NULL, str2 = NULL;
  SourceQualDescData sqdd;
  
  if (biop == NULL || vnp == NULL)
  {
    return NULL;
  }
  
  while (vnp != NULL && str == NULL)
  {
    if (vnp->choice == 0)
    {
      str = GetSourceQualQualValueFromBiop (biop, vnp->data.ptrvalue, fsp);   
    }
    else if (vnp->choice == 1)
    {
      if (StringCmp (vnp->data.ptrvalue, "Location") == 0)
      {
        str = GetLocationFromSource (biop, vnp);
      }
      else if (StringCmp (vnp->data.ptrvalue, "Origin") == 0)
      {
        str = GetOriginFromSource (biop, vnp);
      }
      else if (StringCmp (vnp->data.ptrvalue, "All Notes") == 0)
      {
        sqdd.name = "note";
        sqdd.isOrgMod = TRUE;
        sqdd.subtype = 255;
        str1 = GetSourceQualQualValueFromBiop (biop, &sqdd, fsp);
        sqdd.isOrgMod = FALSE;
        str2 = GetSourceQualQualValueFromBiop (biop, &sqdd, fsp);
        str = StringAppend (str1, str2);
        str1 = MemFree (str1);
        str2 = MemFree (str2);
      }
      else
      {
        str = GetStringFromSource (biop, vnp);
      }
    }
    vnp = vnp->next;
  }
  return str;
}

static CharPtr 
GetSourceQualDescrString 
(SeqDescrPtr  sdp,
 ValNodePtr   vnp,
 FilterSetPtr fsp)
{
  if (sdp == NULL || sdp->choice != Seq_descr_source || vnp == NULL)
  {
    return NULL;
  }
  return GetSourceQualStringFromBioSource (sdp->data.ptrvalue, vnp, fsp);
}

static CharPtr 
GetSourceQualFeatureString 
(SeqFeatPtr   sfp,
 ValNodePtr   vnp,
 FilterSetPtr fsp)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || vnp == NULL)
  {
    return NULL;
  }
  
  return GetSourceQualStringFromBioSource (sfp->data.value.ptrvalue, vnp, fsp);
}

static void ApplyLocationToBiop (BioSourcePtr biop, ApplyValuePtr avp)
{
  if (biop == NULL || avp == NULL || avp->field_list == NULL)
  {
    return;
  }
  
  if (avp->etp == NULL 
      || avp->etp->existing_text_choice == eExistingTextChoiceReplaceOld
      || biop->genome == 0)
  {
    biop->genome = avp->field_list->choice;
  }
}

static void ApplyLocationToSourceFeature (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC)
  {
    ApplyLocationToBiop (sfp->data.value.ptrvalue, userdata);
  }
}

static void ApplyLocationToSourceDescriptor (SeqDescrPtr sdp, Pointer userdata, FilterSetPtr fsp)
{
  if (sdp != NULL && sdp->choice == Seq_descr_source)
  {
    ApplyLocationToBiop (sdp->data.ptrvalue, userdata);
  }
}

static void RemoveLocationFromBiop (BioSourcePtr biop, ApplyValuePtr avp)
{
  if (biop == NULL || avp == NULL || avp->field_list == NULL)
  {
    return;
  }
  
  if (biop->genome == avp->field_list->choice)
  {
    biop->genome = 0;
  }
}

static void RemoveLocationFromSourceFeature (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC)
  {
    RemoveLocationFromBiop (sfp->data.value.ptrvalue, userdata);
  }
}

static void RemoveLocationFromSourceDescriptor (SeqDescrPtr sdp, Pointer userdata, FilterSetPtr fsp)
{
  if (sdp != NULL && sdp->choice == Seq_descr_source)
  {
    RemoveLocationFromBiop (sdp->data.ptrvalue, userdata);
  }
}

static void ApplyOriginToBiop (BioSourcePtr biop, ApplyValuePtr avp)
{
  if (biop == NULL || avp == NULL || avp->field_list == NULL)
  {
    return;
  }
  
  if (avp->etp == NULL 
      || avp->etp->existing_text_choice == eExistingTextChoiceReplaceOld
      || biop->origin == 0)
  {
    biop->origin = avp->field_list->choice;
  }
}

static void ApplyOriginToSourceFeature (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC)
  {
    ApplyOriginToBiop (sfp->data.value.ptrvalue, userdata);
  }
}

static void ApplyOriginToSourceDescriptor (SeqDescrPtr sdp, Pointer userdata, FilterSetPtr fsp)
{
  if (sdp != NULL && sdp->choice == Seq_descr_source)
  {
    ApplyOriginToBiop (sdp->data.ptrvalue, userdata);
  }
}

static void RemoveOriginFromBiop (BioSourcePtr biop, ApplyValuePtr avp)
{
  if (biop == NULL || avp == NULL || avp->field_list == NULL)
  {
    return;
  }
  
  if (biop->origin == avp->field_list->choice)
  {
    biop->origin = 0;
  }
}

static void RemoveOriginFromSourceFeature (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC)
  {
    RemoveOriginFromBiop (sfp->data.value.ptrvalue, userdata);
  }
}

static void RemoveOriginFromSourceDescriptor (SeqDescrPtr sdp, Pointer userdata, FilterSetPtr fsp)
{
  if (sdp != NULL && sdp->choice == Seq_descr_source)
  {
    RemoveOriginFromBiop (sdp->data.ptrvalue, userdata);
  }
}

static void ApplyStringToSource (BioSourcePtr biop, Pointer userdata)
{
  ApplyValuePtr         avp;

  if (biop == NULL || userdata == NULL)
  {
    return;
  }
  
  avp = (ApplyValuePtr) userdata;
  if (avp == NULL || avp->field_list == NULL)
  {
    return;
  }

  if (StringHasNoText (avp->text_to_replace))
  {
    if (biop->org == NULL)
    {
      biop->org = OrgRefNew();
    }
    if (biop->org->orgname == NULL 
        && (avp->field_list->choice == 3
            || avp->field_list->choice == 4))
    {
      biop->org->orgname = OrgNameNew ();
    }
  }
  
  switch (avp->field_list->choice)
  {
    case 1:
      biop->org->taxname = HandleApplyValue (biop->org->taxname, avp);
      break;
    case 2:
      biop->org->common = HandleApplyValue (biop->org->common, avp);
      break;
    case 3:
      biop->org->orgname->lineage = HandleApplyValue (biop->org->orgname->lineage, avp);
      break;
    case 4:
      biop->org->orgname->div = HandleApplyValue (biop->org->orgname->div, avp);
      break;
  }  
}

static void ApplyStringToSourceFeature (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC)
  {
    ApplyStringToSource (sfp->data.value.ptrvalue, userdata);
  }
}

static void ApplyStringToSourceDescriptor (SeqDescrPtr sdp, Pointer userdata, FilterSetPtr fsp)
{
  if (sdp != NULL && sdp->choice == Seq_descr_source)
  {
    ApplyStringToSource (sdp->data.ptrvalue, userdata);
  }
}

static void RemoveStringFromSource (BioSourcePtr biop, Pointer userdata)
{
  ApplyValuePtr avp;
  
  if (biop == NULL || biop->org == NULL || userdata == NULL)
  {
    return;
  }
  
  avp = (ApplyValuePtr) userdata;
  if (avp == NULL || avp->field_list == NULL)
  {
    return;
  }
  
  switch (avp->field_list->choice)
  {
    case 1:
      biop->org->taxname = MemFree (biop->org->taxname);
      break;
    case 2:
      biop->org->common = MemFree (biop->org->common);
      break;
    case 3:
      if (biop->org->orgname != NULL)
      {
        biop->org->orgname->lineage = MemFree (biop->org->orgname->lineage);
      }
      break;
    case 4:
      if (biop->org->orgname != NULL)
      {
        biop->org->orgname->div = MemFree (biop->org->orgname->div);
      }
      break;
  }   
}

static void RemoveStringFromSourceFeature (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC)
  {
    RemoveStringFromSource (sfp->data.value.ptrvalue, userdata);
  }
}

static void RemoveStringFromSourceDescriptor (SeqDescrPtr sdp, Pointer userdata, FilterSetPtr fsp)
{
  if (sdp != NULL && sdp->choice == Seq_descr_source)
  {
    RemoveStringFromSource (sdp->data.ptrvalue, userdata);
  }
}

static void SqnRemoveSourceQualFromBioSource (BioSourcePtr biop, Pointer userdata, FilterSetPtr fsp);

static CharPtr RemoveSubfield (CharPtr orig, Uint1 subfield)
{
  CharPtr cp, cp2;
  Int4    num_colons = 0;
  CharPtr new_val;

  cp = StringChr (orig, ':');
  while (cp != NULL) {
    num_colons ++;
    cp = StringChr (cp + 1, ':');
  }
  new_val = orig;

  if (subfield == 0) {
    new_val = MemFree (new_val);
  } else if (subfield == 1) {
    if (num_colons == 0) {
      /* no INST, just specID, nothing to remove */
    } else if (num_colons == 1) {
      cp = StringChr (orig, ':');
      new_val = StringSave (cp + 1);
      orig = MemFree (orig);
    } else {
      cp = StringChr (orig, ':');
      new_val = StringSave (cp);
      orig = MemFree (orig);
    }
  } else if (subfield == 2) {
    if (num_colons == 0 || num_colons == 1) {
      /* no COLL, nothing to remove */
    } else {
      cp = StringChr (orig, ':');
      cp2 = StringChr (cp + 1, ':');
      new_val = (CharPtr) MemNew (sizeof (Char) * (cp - orig + StringLen (cp2) + 1));
      StringNCpy (new_val, orig, cp - orig);
      StringCat (new_val, cp2);
      orig = MemFree (orig);
      if (StringCmp (new_val, ":") == 0) {
        new_val = MemFree (new_val);
      }
    }
  } else {
    /* specid */
    if (num_colons == 0) {
      new_val = MemFree (new_val);
    } else {
      cp = StringRChr (new_val, ':');
      *(cp + 1) = 0;
    }
  }
  return new_val;
}


static CharPtr ApplyToSubfield (CharPtr orig, ApplyValuePtr avp, Uint1 subfield)
{
  CharPtr cp, cp2;
  Int4    num_colons = 0, len;
  CharPtr new_field, new_val;

  cp = StringChr (orig, ':');
  while (cp != NULL) {
    num_colons ++;
    cp = StringChr (cp + 1, ':');
  }

  new_val = orig;
  if (subfield == 0) {
    /* whole field */
    new_val = HandleApplyValue (orig, avp);
  } else if (subfield == 1) {
    /* INST */
    if (num_colons == 0) {
      /* all we have is the specID, we are adding INST */
      new_field = HandleApplyValue (NULL, avp);
      if (!StringHasNoText (new_field)) {
        new_val = (CharPtr) MemNew (sizeof (Char) * (StringLen (new_field) + StringLen (orig) + 2));
        sprintf (new_val, "%s:%s", new_field, orig == NULL ? "" : orig);
        orig = MemFree (orig);
      } else {
        new_val = orig;
      }
      new_field = MemFree (new_field);
    } else {
      cp = StringChr (orig, ':');
      len = cp - orig + 1;
      new_field = (CharPtr) MemNew (sizeof (Char) * len);
      StringNCpy (new_field, orig, len - 1);
      new_field[len - 1] = 0;
      new_field = HandleApplyValue (new_field, avp);
      new_val = (CharPtr) MemNew (sizeof (Char) * (StringLen (new_field) + StringLen (cp) + 1));
      sprintf (new_val, "%s%s", new_field == NULL ? "" : new_field, cp);
      orig = MemFree (orig);
      new_field = MemFree (new_field);
    }
  } else if (subfield == 2) {
    /* COLL */
    if (num_colons == 0) {
      /* add blank INST, then COLL, then existing specID */
      new_field = HandleApplyValue (NULL, avp);
      if (!StringHasNoText (new_field)) {
        new_val = (CharPtr) MemNew (sizeof (Char) * (3 + StringLen (new_field) +  StringLen (orig)));
        sprintf (new_val, ":%s:%s", new_field, orig == NULL ? "" : orig);
        orig = MemFree (orig);
      } else {
        new_val = orig;
      }
      new_field = MemFree (new_field);
    } else if (num_colons == 1) {
      /* have INST and specID */
      new_field = HandleApplyValue (NULL, avp);
      if (!StringHasNoText (new_field)) {
        new_val = (CharPtr) MemNew (sizeof (Char) * (3 + StringLen (new_field) +  StringLen (orig)));
        cp = StringChr (orig, ':');
        len = cp - orig;
        StringNCpy (new_val, orig, len);
        StringCat (new_val, ":");
        StringCat (new_val, new_field);
        StringCat (new_val, cp);
        orig = MemFree (orig);
      } else {
        new_val = orig;
      }
      new_field = MemFree (new_field);
    } else {
      cp = StringChr (orig, ':');
      cp2 = StringChr (cp + 1, ':');
      len = cp2 - cp;
      new_field = (CharPtr) MemNew (sizeof (Char) * len);
      StringNCpy (new_field, cp + 1, len - 1);
      new_field[len - 1] = 0;
      new_field = HandleApplyValue (new_field, avp);
      new_val = (CharPtr) MemNew (sizeof (Char) * (cp - orig + 1 + StringLen (new_field) + StringLen (cp2) + 1));
      StringNCpy (new_val, orig, cp - orig + 1);
      StringCat (new_val, new_field);
      StringCat (new_val, cp2);
      orig = MemFree (orig);
      new_field = MemFree (new_field);
    }
  } else {
    /* specid */
    if (num_colons == 0) {
      new_val = HandleApplyValue (orig, avp);
    } else {
      cp = StringRChr (orig, ':');
      new_field = StringSave (cp + 1);
      new_field = HandleApplyValue (new_field, avp);
      new_val = (CharPtr) MemNew (sizeof (Char) + (cp - orig + 1 + StringLen (new_field) + 1));
      StringNCpy (new_val, orig, cp - orig + 1);
      StringCat (new_val, new_field);
      orig = MemFree (orig);
      new_field = MemFree (new_field);
    }
  }
  return new_val;
}


static void ApplySourceQualBioSourceCallback (BioSourcePtr biop, Pointer userdata, FilterSetPtr fsp)
{
  ApplyValuePtr         avp;
  SourceQualDescPtr     sqdp;
  OrgModPtr             mod = NULL, last_mod = NULL;
  SubSourcePtr          ssp = NULL, last_ssp = NULL;
  Boolean               found = FALSE;
  Boolean               is_nontext;
  
  if (biop == NULL || userdata == NULL)
  {
    return;
  }
  
  avp = (ApplyValuePtr) userdata;
  if (avp->field_list == NULL || avp->field_list->data.ptrvalue == NULL)
  {
    return;
  }
  
  /* don't apply to All Notes */
  if (avp->field_list->choice != 0)
  {
    return;
  }
  
  sqdp = (SourceQualDescPtr) avp->field_list->data.ptrvalue;
  
  is_nontext = IsNonTextModifier (sqdp->name);
  if (is_nontext && StringHasNoText (avp->new_text))
  {
    SqnRemoveSourceQualFromBioSource (biop, avp, fsp);
    return;
  }
  
  if (!StringHasNoText (avp->text_to_replace))
  {
    found = TRUE;
  }
  
  if (sqdp->isOrgMod)
  {
    if (biop->org == NULL)
    {
      biop->org = OrgRefNew ();
    }
    if (biop->org != NULL && biop->org->orgname == NULL)
    {
      biop->org->orgname = OrgNameNew ();
    }
    if (biop->org != NULL && biop->org->orgname != NULL)
    {
      mod = biop->org->orgname->mod;
      while (mod != NULL)
      {
        if (mod->subtype == sqdp->subtype)
        {
          if (!is_nontext)
          {
            mod->subname = ApplyToSubfield (mod->subname, avp, sqdp->subfield);
          }
          found = TRUE;
        }
        last_mod = mod;
        mod = mod->next;
      }
      if (!found)
      {
        mod = OrgModNew ();
        mod->subtype = sqdp->subtype;
        if (is_nontext)
        {
          mod->subname = StringSave ("");
        }
        else
        {
          mod->subname = ApplyToSubfield (mod->subname, avp, sqdp->subfield);
        }
        if (last_mod == NULL)
        {
          biop->org->orgname->mod = mod;
        }
        else
        {
          last_mod->next = mod;
        }
      }
    }
  }
  else
  {
    ssp = biop->subtype;
    while (ssp != NULL)
    {
      if (ssp->subtype == sqdp->subtype)
      {
        if (! is_nontext)
        {
          ssp->name = ApplyToSubfield (ssp->name, avp, sqdp->subfield);
        }
        found = TRUE;
      }
      last_ssp = ssp;
      ssp = ssp->next;
    }
    if (!found)
    {
      ssp = SubSourceNew ();
      ssp->subtype = sqdp->subtype;
      if (is_nontext)
      {
        ssp->name = StringSave ("");
      }
      else
      {
        ssp->name = ApplyToSubfield (ssp->name, avp, sqdp->subfield);
      }
      if (last_ssp == NULL)
      {
        biop->subtype = ssp;
      }
      else
      {
        last_ssp->next = ssp;
      }
    }
  }
}

static void ApplySourceQualFeatureCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || userdata == NULL)
  {
    return;
  }
  ApplySourceQualBioSourceCallback (sfp->data.value.ptrvalue, userdata, fsp);
}

static void ApplySourceQualDescriptorCallback (SeqDescrPtr sdp, Pointer userdata, FilterSetPtr fsp)
{
  if (sdp == NULL || sdp->choice != Seq_descr_source || userdata == NULL)
  {
    return;
  }
  ApplySourceQualBioSourceCallback (sdp->data.ptrvalue, userdata, fsp);
}

static void RemoveOrgModFromBioSource (BioSourcePtr biop, Uint1 subtype, Uint1 subfield, FilterSetPtr fsp)
{
  OrgModPtr mod = NULL, prev_mod = NULL, next_mod;
  
  if (biop != NULL && biop->org != NULL && biop->org->orgname != NULL)
  {
    mod = biop->org->orgname->mod;
    while (mod != NULL)
    {
      next_mod = mod->next;
      if (mod->subtype == subtype
          && OrgModSpecialMatch (mod, fsp))
      {
        mod->subname = RemoveSubfield (mod->subname, subfield);
        if (StringHasNoText (mod->subname)) {
          if (prev_mod == NULL)
          {
            biop->org->orgname->mod = mod->next;
          }
          else
          {
            prev_mod->next = mod->next;
          }
          mod->next = NULL;
          OrgModFree (mod);
        } else {
          prev_mod = mod;
        }
      }
      else
      {
        prev_mod = mod;
      }
      mod = next_mod;
    }
  }
}

static void RemoveSubSourceFromBioSource (BioSourcePtr biop, Uint2 subtype, Uint1 subfield, FilterSetPtr fsp)
{
  SubSourcePtr ssp = NULL, prev_ssp = NULL, next_ssp;

  if (biop != NULL)
  {
    ssp = biop->subtype;
    while (ssp != NULL)
    { 
      next_ssp = ssp->next;
      if (ssp->subtype == subtype
          && SubSrcSpecialMatch (ssp, fsp))
      {
        ssp->name = RemoveSubfield (ssp->name, subfield);
        if (StringHasNoText (ssp->name)) {
          if (prev_ssp == NULL)
          {
            biop->subtype = ssp->next;
          }
          else
          {
            prev_ssp->next = ssp->next;
          }
          ssp->next = NULL;
          SubSourceFree (ssp);
        } else {
          prev_ssp = ssp;
        }
      }
      else
      {
        prev_ssp = ssp;
      }
      ssp = next_ssp;
    }
  }
}

static void SqnRemoveSourceQualFromBioSource (BioSourcePtr biop, Pointer userdata, FilterSetPtr fsp)
{
  ApplyValuePtr         avp;
  SourceQualDescPtr     sqdp;
  
  if (biop == NULL || userdata == NULL)
  {
    return;
  }
  
  avp = (ApplyValuePtr) userdata;
  if (avp->field_list == NULL || avp->field_list->data.ptrvalue == NULL)
  {
    return;
  }
  if (avp->field_list->choice == 0)
  {
    sqdp = (SourceQualDescPtr) avp->field_list->data.ptrvalue;
    if (sqdp->isOrgMod)
    {
      RemoveOrgModFromBioSource (biop, sqdp->subtype, sqdp->subfield, fsp);
    }
    else
    {
      RemoveSubSourceFromBioSource (biop, sqdp->subtype, sqdp->subfield, fsp);
    }
  }
  else if (avp->field_list->choice == 1 
           && StringCmp (avp->field_list->data.ptrvalue, "All Notes") == 0)
  {
    RemoveOrgModFromBioSource (biop, 255, 0, fsp);
    RemoveSubSourceFromBioSource (biop, 255, 0, fsp);
  }
}

static void RemoveSourceQualFromSourceFeature (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || userdata == NULL)
  {
    return;
  }
  SqnRemoveSourceQualFromBioSource (sfp->data.value.ptrvalue, userdata, fsp);
}

static void RemoveSourceQualFromSourceDescriptor (SeqDescrPtr sdp, Pointer userdata, FilterSetPtr fsp)
{
  if (sdp == NULL || sdp->choice != Seq_descr_source || userdata == NULL)
  {
    return;
  }
  SqnRemoveSourceQualFromBioSource (sdp->data.ptrvalue, userdata, fsp);
}

static DialoG SourceLocationSelectionDialogEx 
(GrouP h,
 Boolean                  allow_multi,
 Boolean                  allow_discontinued,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)

{
  DialoG     dlg;
  ValNodePtr choice_list = NULL;
  Nlm_EnumFieldAssocPtr eap = biosource_genome_alistX;
  
  if (eap == NULL)
  {
    return NULL;
  }

  while (eap->name != NULL)
  {
    if (!StringHasNoText (eap->name) 
        && (allow_discontinued ||
            (StringCmp (eap->name, "Insertion Sequence") != 0
             && StringCmp (eap->name, "Transposon") != 0)))
    {
      ValNodeAddPointer (&choice_list, eap->value, StringSave (eap->name));
    }
    eap++;
  }
  
  /* note - the ValNodeSelectionDialog will free the qual_choice_list
   * when done */                                            
  dlg = ValNodeSelectionDialog (h, choice_list, TALL_SELECTION_LIST,
                                ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "Location",
                                change_notify, change_userdata, allow_multi);

  return dlg;
}


static DialoG SourceLocationSelectionDialog
(GrouP h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  return SourceLocationSelectionDialogEx (h, allow_multi, FALSE,
                                          change_notify, change_userdata);
}


static DialoG SourceLocationDiscSelectionDialog
(GrouP h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  return SourceLocationSelectionDialogEx (h, allow_multi, TRUE,
                                          change_notify, change_userdata);
}


static DialoG SourceOriginSelectionDialog
(GrouP h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)

{
  return EnumAssocSelectionDialog (h, origin_alist, "Origin",
                                   allow_multi, change_notify, change_userdata);
}

static DialoG SourceStringSelectionDialog
(GrouP h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)

{
  return EnumAssocSelectionDialog (h, orgref_textfield_alist, "string",
                                   allow_multi, change_notify, change_userdata);
}


static CharPtr 
GetStringFromStringDescriptor 
(SeqDescrPtr  sdp,
 ValNodePtr   vnp,
 FilterSetPtr fsp)
{
  if (sdp == NULL || vnp == NULL || sdp->choice != vnp->data.intvalue)
  {
    return NULL;
  }
  return StringSave (sdp->data.ptrvalue);
}

static DialoG GBQualSelectionDialog 
(GrouP                    h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)

{
  DialoG                    dlg;
  ValNodePtr            choice_name_list = NULL;
  Int4                      i;
     
  for (i = 0; i < NumEditQualifiers; i++)
  {
    ValNodeAddPointer (&choice_name_list, 0, EditQualifierList[i].name);
  }

  dlg = SelectionDialog (h, change_notify, change_userdata, 
                         allow_multi, "genbank qualifier",
                         choice_name_list, TALL_SELECTION_LIST);
  ValNodeFree (choice_name_list); 
  return dlg;
}

static Int4 GetDbtagStringLen (DbtagPtr db_tag)
{
  Int4 len;
  
  if (db_tag == NULL)
  {
    return 0;
  }
  
  len = StringLen (db_tag->db) + 2;
  if (db_tag->tag != NULL)
  {
    if (db_tag->tag->str != NULL)
    {
      len += StringLen (db_tag->tag->str);
    }
    else
    {
      len += 10;
    }
  }
  return len;
}

static CharPtr GetDbxrefString (SeqFeatPtr sfp)
{
  ValNodePtr vnp;
  Int4       len = 0;
  CharPtr    str = NULL, cp;
  
  if (sfp == NULL || sfp->dbxref == NULL)
  {
    return NULL;
  }
  
  for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next)
  {
    len += GetDbtagStringLen (vnp->data.ptrvalue) + 1;
  }
  
  if (len == 0)
  {
    return NULL;
  }
  
  str = (CharPtr) MemNew ((len + 1) * sizeof (Char));
  if (str != NULL)
  {
    for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next)
    {
      cp = GetDbtagString (vnp->data.ptrvalue);
      if (cp != NULL)
      {
        StringCat (str, cp);
        StringCat (str, ";");
      }
    }
  }
  if (StringLen (str) >1)
  {
    /* remove final semicolon */
    str [StringLen (str) - 2] = 0;
  }
  return str;
}

static CharPtr GetGBQualString (SeqFeatPtr sfp, ValNodePtr vnp, FilterSetPtr fsp)
{
  GBQualPtr  gbqual;
  CharPtr    str = NULL;
  Int4       field_choice;
  RnaRefPtr  rrp;
  GeneRefPtr grp;
  
  if (sfp == NULL)
  {
    return NULL;
  }
  while (vnp != NULL && str == NULL)
  {
    field_choice = vnp->data.intvalue;
    if (StringICmp (EditQualifierList[field_choice - 1].name, "note") == 0)
    {
      if (!StringHasNoText (sfp->comment))
      {
        str = StringSave (sfp->comment);
      }
    }
    else if (StringICmp (EditQualifierList[field_choice - 1].name, "exception") == 0)
    {
      if (!StringHasNoText (sfp->except_text))
      {
        str = StringSave (sfp->except_text);
      }
    }
    else if (StringICmp (EditQualifierList[field_choice - 1].name, "product") == 0
             && sfp->data.choice == SEQFEAT_RNA)
    {
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->ext.choice == 1
          && !StringHasNoText (rrp->ext.value.ptrvalue))
      {
        str = StringSave (rrp->ext.value.ptrvalue);        
      }
    }
    else if (StringICmp (EditQualifierList [field_choice - 1].name, "map") == 0
             && sfp->data.choice == SEQFEAT_GENE)
    {
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (grp != NULL && !StringHasNoText (grp->maploc))
      {
        str = StringSave (grp->maploc);
      }
    }
    else if (StringICmp (EditQualifierList [field_choice - 1].name, "db_xref") == 0)
    {
        str = GetDbxrefString (sfp);
    }
    else if (StringICmp (EditQualifierList [field_choice - 1].name, "evidence") == 0)
    {
      if (sfp->exp_ev == 1)
      {
        str = StringSave ("experimental");
      }
      else if (sfp->exp_ev == 2)
      {
        str = StringSave ("non-experimental");
      }
    }
    else
    {
      gbqual = sfp->qual;
      while (gbqual != NULL && str == NULL)
      {
        if (field_choice > 0 && field_choice < NumEditQualifiers + 1
            && StringCmp (gbqual->qual, EditQualifierList[field_choice - 1].name) == 0
            && !StringHasNoText (gbqual->val))
        {
          str = StringSave (gbqual->val);
        }
        gbqual = gbqual->next;
      }
    }
    vnp = vnp->next;
  }
  return str;
}

static DbtagPtr DbtagFromString (CharPtr str)
{
  DbtagPtr db_tag;
  CharPtr  cp, idp;
  Boolean  use_id = TRUE;
  
  if (StringHasNoText (str))
  {
    return NULL;
  }
  
  cp = StringChr (str, ':');
  if (cp == NULL)
  {
    return NULL;
  }
  
  db_tag = DbtagNew ();
  if (db_tag != NULL)
  {
    db_tag->db = (CharPtr) MemNew ((cp - str + 1) * sizeof (Char));
    StringNCpy (db_tag->db, str, cp - str);
    db_tag->db [cp - str] = 0;
    
    idp = cp + 1;
    while (*idp != 0 && use_id)
    {
      if (!isdigit (*idp))
      {
        use_id = FALSE;
      }
      idp++;
    }
    db_tag->tag = ObjectIdNew ();
    if (use_id)
    {
      db_tag->tag->id = atoi (cp + 1);
      db_tag->tag->str = NULL;
    }
    else
    {
      db_tag->tag->id = 0;
      db_tag->tag->str = StringSave (cp + 1);
    }
  }
  return db_tag;
}

static void HandleApplyValueForObjectID (ObjectIdPtr oip, ApplyValuePtr avp)
{
  CharPtr tmp_str = NULL;
  CharPtr cp;
  Boolean is_num = TRUE;
  
  if (oip == NULL || avp == NULL)
  {
    return;
  }
  
  if (oip->str == NULL)
  {
    tmp_str = (CharPtr) MemNew (10 * sizeof (Char));
    sprintf (tmp_str, "%d", oip->id);
    
  }
  else
  {
    tmp_str = oip->str;
    oip->str = NULL;
  }
  tmp_str = HandleApplyValue (tmp_str, avp);
  
  for (cp = tmp_str; cp != NULL && *cp != 0 && is_num; cp++)
  {
    if (!isdigit (*cp))
    {
      is_num = FALSE;
    }
  }
  if (is_num)
  {
    oip->id = atoi (tmp_str);
    tmp_str = MemFree (tmp_str);
  }
  else
  {
    oip->str = tmp_str;
  }
}

static void EditDbxref (SeqFeatPtr sfp, ApplyValuePtr avp, FilterSetPtr fsp)
{
  ValNodePtr vnp;
  DbtagPtr   db_tag_new, db_tag_to_replace, db_tag;
  CharPtr    cp;
  
  if (sfp == NULL || avp == NULL)
  {
    return;
  }
  
  db_tag_new = DbtagFromString (avp->new_text);
  db_tag_to_replace = DbtagFromString (avp->text_to_replace);
  
  if (db_tag_new == NULL && db_tag_to_replace == NULL)
  {
    for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next)
    {
      db_tag = (DbtagPtr) vnp->data.ptrvalue;
      if (db_tag != NULL)
      {
        db_tag->db = HandleApplyValue (db_tag->db, avp);
        HandleApplyValueForObjectID (db_tag->tag, avp);
      }
    }
  }
  else if (db_tag_new != NULL && db_tag_to_replace != NULL)
  {
    for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next)
    {
      db_tag = (DbtagPtr) vnp->data.ptrvalue;
      if (db_tag != NULL)
      {
        cp = GetDbtagString (db_tag);
        cp = HandleApplyValue (cp, avp);
        db_tag = DbtagFree (db_tag);
        db_tag = DbtagFromString (cp);
        vnp->data.ptrvalue = db_tag;
        cp = MemFree (cp);
      }
    }
  }  
}

static ValNodePtr RemoveDbxrefList (ValNodePtr vnp)

{
  ValNodePtr  next;

  while (vnp != NULL) {
    next = vnp->next;
    DbtagFree ((DbtagPtr) vnp->data.ptrvalue);
    MemFree (vnp);
    vnp = next;
  }
  return NULL;
}

static void SetDbxrefFromString (SeqFeatPtr sfp, ApplyValuePtr avp, FilterSetPtr fsp)
{
  DbtagPtr new_tag;
  
  if (sfp == NULL || avp == NULL)
  {
    return;
  }
  
  if (StringHasNoText (avp->text_to_replace))
  {
    /* this is an apply */
    if (avp->etp != NULL)
    {
      if (avp->etp->existing_text_choice == eExistingTextChoiceLeaveOld
          && sfp->dbxref != NULL)
      {
        return;
      }
      else if (avp->etp->existing_text_choice == eExistingTextChoiceReplaceOld)
      {
        sfp->dbxref = RemoveDbxrefList (sfp->dbxref);
      }
    }
    
    new_tag = DbtagFromString (avp->new_text);
    if (new_tag != NULL)
    {
      ValNodeAddPointer (&(sfp->dbxref), 0, new_tag);
    }
  }
  else
  {
    EditDbxref (sfp, avp, fsp);
  }
  
}

static void SetRNAProduct (SeqFeatPtr sfp, ApplyValuePtr avp);

static void SetGBQualString (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  GBQualPtr     gbqual, gbqual_last = NULL;
  ApplyValuePtr avp;
  Boolean       found = FALSE;
  GeneRefPtr    grp;
  
  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  avp = (ApplyValuePtr) userdata;
  if (avp->field_list == NULL 
      || avp->field_list->data.intvalue < 1 
      || avp->field_list->data.intvalue >= NumEditQualifiers + 1)
  {
    return;
  }

  if (StringICmp (EditQualifierList[avp->field_list->data.intvalue - 1].name, "note") == 0)
  {
    sfp->comment = HandleApplyValue (sfp->comment, avp);
    return;
  }
  else if (StringICmp (EditQualifierList[avp->field_list->data.intvalue - 1].name, "exception") == 0)
  {
    sfp->except_text = HandleApplyValue (sfp->except_text, avp);
    return;
  }
  else if (StringICmp (EditQualifierList[avp->field_list->data.intvalue - 1].name, "product") == 0
           && sfp->data.choice == SEQFEAT_RNA)
  {
    SetRNAProduct (sfp, avp);
    return;
  }
  else if (StringICmp (EditQualifierList [avp->field_list->data.intvalue - 1].name, "map") == 0
           && sfp->data.choice == SEQFEAT_GENE)
  {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (grp != NULL)
    {
      grp->maploc = HandleApplyValue (grp->maploc, avp);
      return;
    }
  }
  else if (StringICmp (EditQualifierList [avp->field_list->data.intvalue - 1].name, "db_xref") == 0)
  {
    SetDbxrefFromString (sfp, avp, fsp);
    return;
  }


  if (!StringHasNoText (avp->text_to_replace))
  {
    found = TRUE;
  }
    
  gbqual = sfp->qual;
  while (gbqual != NULL)
  {
    gbqual_last = gbqual;
    if (StringCmp (gbqual->qual, EditQualifierList[avp->field_list->data.intvalue - 1].name) == 0)
    {
      gbqual->val = HandleApplyValue (gbqual->val, avp);
      found = TRUE;
    }
    gbqual = gbqual->next;
  }
  if (!found)
  {
    gbqual = GBQualNew ();
    if (gbqual != NULL)
    {
      gbqual->qual = StringSave (EditQualifierList[avp->field_list->data.intvalue - 1].name);
      gbqual->val = StringSave (avp->new_text);
      if (gbqual_last == NULL)
      {
        sfp->qual = gbqual;
      }
      else
      {
        gbqual_last->next = gbqual;
      }
    }
  }
}

static void RemoveGBQualField (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  GBQualPtr     gbqual, gbqual_last = NULL, gbqual_next;
  ApplyValuePtr avp;
  RnaRefPtr     rrp;
  GeneRefPtr    grp;
  
  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  avp = (ApplyValuePtr) userdata;
  if (avp->field_list == NULL 
      || avp->field_list->data.intvalue < 1 
      || avp->field_list->data.intvalue >= NumEditQualifiers + 1)
  {
    return;
  }
  
  if (StringICmp (EditQualifierList[avp->field_list->data.intvalue - 1].name, "note") == 0)
  {
    sfp->comment = MemFree (sfp->comment);
    return;
  }
  else if (StringICmp (EditQualifierList[avp->field_list->data.intvalue - 1].name, "exception") == 0)
  {
    sfp->except_text = MemFree (sfp->except_text);
    return;
  }
  else if (StringICmp (EditQualifierList[avp->field_list->data.intvalue - 1].name, "product") == 0
           && sfp->data.choice == SEQFEAT_RNA)
  {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rrp != NULL && rrp->ext.choice == 1)
    {
      rrp->ext.choice = 0;
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
    }
    return;
  }
  else if (StringICmp (EditQualifierList [avp->field_list->data.intvalue - 1].name, "map") == 0
           && sfp->data.choice == SEQFEAT_GENE)
  {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (grp != NULL)
    {
      grp->maploc = MemFree (grp->maploc);
    }
  }
  else if (StringICmp (EditQualifierList [avp->field_list->data.intvalue - 1].name, "db_xref") == 0)
  {
    sfp->dbxref = RemoveDbxrefList (sfp->dbxref);
  }


  gbqual = sfp->qual;
  while (gbqual != NULL)
  {
    gbqual_next = gbqual->next;
    if (StringCmp (gbqual->qual, EditQualifierList[avp->field_list->data.intvalue - 1].name) == 0)
    {
      if (gbqual_last == NULL)
      {
        sfp->qual = gbqual->next;
      }
      else
      {
        gbqual_last->next = gbqual->next;
      }
      gbqual->next = NULL;
      GBQualFree (gbqual);
    }
    else
    {
      gbqual_last = gbqual;
    }
    gbqual = gbqual_next;
  }
}


static void RemoveGoTerm (SeqFeatPtr sfp, CharPtr term_name)
{
  UserObjectPtr       uop;
  ObjectIdPtr         oip;
  UserFieldPtr  curr, curr_next, prev = NULL;
  
  if (sfp == NULL || sfp->ext == NULL || StringHasNoText (term_name))
  {
    return;
  }

  uop = (UserObjectPtr) sfp->ext;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "GeneOntology") != 0) return;
  
  for (curr = uop->data; curr != NULL; curr = curr_next)
  {
    curr_next = curr->next;
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, term_name) == 0)
    {
      if (prev == NULL)
      {
        uop->data = curr->next;
      }
      else
      {
        prev->next = curr->next;
      }
      curr->next = NULL;
      UserFieldFree (curr);
    }
    else
    {
      prev = curr;
    }
  }
 
  if (uop->data == NULL)
  {
    sfp->ext = UserObjectFree (sfp->ext);
  }
 
}


static void RemoveGBQualByNameConstraint (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  GBQualPtr           gbqual, gbqual_last = NULL, gbqual_next;
  RnaRefPtr           rrp;
  GeneRefPtr          grp;
  StringConstraintXPtr scp;
  
  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  scp = (StringConstraintXPtr) userdata;
  if (scp == NULL)
  {
    return;
  }
  
  if (DoesStringMatchConstraintX ("note", scp))
  {
    sfp->comment = MemFree (sfp->comment);
  }
  
  if (DoesStringMatchConstraintX ("exception", scp))
  {
    sfp->except_text = MemFree (sfp->except_text);
  }

  if (DoesStringMatchConstraintX ("product", scp) && sfp->data.choice == SEQFEAT_RNA)
  {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rrp != NULL && rrp->ext.choice == 1)
    {
      rrp->ext.choice = 0;
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
    }
  }
  
  if (DoesStringMatchConstraintX ("map", scp) && sfp->data.choice == SEQFEAT_GENE)
  {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (grp != NULL)
    {
      grp->maploc = MemFree (grp->maploc);
    }
  }
  
  if (DoesStringMatchConstraintX ("db_xref", scp))
  {
    sfp->dbxref = RemoveDbxrefList (sfp->dbxref);
  }
  
  if (DoesStringMatchConstraintX ("go_process", scp))
  {
    RemoveGoTerm (sfp, "process");
  }
  
  if (DoesStringMatchConstraintX ("go_function", scp))
  {
    RemoveGoTerm (sfp, "function");
  }

  if (DoesStringMatchConstraintX ("go_component", scp))
  {
    RemoveGoTerm (sfp, "component");
  }

  gbqual = sfp->qual;
  while (gbqual != NULL)
  {
    gbqual_next = gbqual->next;
    if (DoesStringMatchConstraintX (gbqual->qual, scp))
    {
      if (gbqual_last == NULL)
      {
        sfp->qual = gbqual->next;
      }
      else
      {
        gbqual_last->next = gbqual->next;
      }
      gbqual->next = NULL;
      GBQualFree (gbqual);
    }
    else
    {
      gbqual_last = gbqual;
    }
    gbqual = gbqual_next;
  }
}


static CharPtr GetFeatureNote (SeqFeatPtr sfp, ValNodePtr vnp, FilterSetPtr fsp)
{
  if (sfp == NULL || StringHasNoText (sfp->comment))
  {
    return NULL;
  }
  else
  {
    return StringSave (sfp->comment);
  }
}


static void StripFieldNameFromText(CharPtr text,
                                   NameFromValNodeProc field_name_func, 
                                   ValNodePtr          field_vnp)
{
  CharPtr src, dst;
  CharPtr field_name = NULL;
  Uint4    field_name_len;
  
  if (StringHasNoText (text) || field_name_func == NULL)
  {
    return;
  }
  
  field_name = field_name_func(field_vnp);
  field_name_len = StringLen (field_name);
    
  if (!StringHasNoText (field_name) && StringNICmp(text, field_name, field_name_len) == 0
        && StringLen (text) > field_name_len 
        && text[field_name_len] == ' ')
  {
    src = text + field_name_len + 1;
    while (*src == ' ')
    {
      src++;
    }
    dst = text;
    while (*src != 0)
    {
      *dst = *src;
      dst++;
      src++;
    }
    *dst = 0;
  }
  field_name = MemFree (field_name);
}


static CharPtr ReplaceStringForParse(CharPtr src_text, TextPortionXPtr text_portion)
{
  CharPtr         found_loc;
  Int4            found_len;
  CharPtr         dst_txt = NULL, cp;
  
  if (src_text == NULL || text_portion == NULL) {
    return NULL;
  }
  /* for parse, replace the source text with the portion we want */
  found_loc = NULL;
  found_len = 0;
  FindTextPortionXInString (src_text, text_portion, &found_loc, &found_len);
  if (found_loc != NULL) 
  {
    dst_txt = (CharPtr)MemNew (found_len + 1);
    StringNCpy (dst_txt, found_loc, found_len);
    dst_txt[found_len] = 0;
    cp = found_loc + found_len;
    while (*cp != 0) {
      *found_loc = *cp;
      found_loc++;
      cp++;
    }
    *found_loc = 0;
  }
  return dst_txt;
}


static void ConvertFeatureFieldCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  ConvertFieldPtr cfp;
  CharPtr         src_text;
  CharPtr         dst_text = NULL;
  ApplyValueData  avd;
  
  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  
  cfp = (ConvertFieldPtr) userdata;
  
  if (cfp->src_field_list == NULL || cfp->dst_field_list == NULL
      || cfp->get_str_func == NULL || cfp->set_str_func == NULL
      || cfp->remove_str_func == NULL)
  {
    return;
  }
  
  avd.text_to_replace = NULL;
  avd.where_to_replace = EditApplyFindLocation_anywhere;
  src_text = (cfp->get_str_func)(sfp, cfp->src_field_list, cfp->fsp);
  
  if (cfp->strip_name_from_text)
  {
    StripFieldNameFromText(src_text, cfp->name_field_func, cfp->dst_field_list);
  }
  
  if (StringHasNoText (src_text))
  {
    src_text = MemFree (src_text);
    return;
  }
  
  if (cfp->convert_type == CONVERT_TYPE_SWAP)
  {
    /* for swap, we already know what we're doing with existing text - 
     * we're putting it where the new text used to be. */
    avd.etp = NULL;
    /* and we need to get the destination text before it's replaced */
    dst_text = (cfp->get_str_func)(sfp, cfp->dst_field_list, cfp->fsp); 
    if (cfp->strip_name_from_text)
    {
      StripFieldNameFromText(dst_text, cfp->name_field_func, cfp->src_field_list);
    }
  }
  else if (cfp->convert_type == CONVERT_TYPE_PARSE)
  {
    /* for parse, replace the source text with the portion we want */
    dst_text = src_text;
    src_text = ReplaceStringForParse(src_text, cfp->text_portion);
    if (src_text == NULL) 
    {
      dst_text = MemFree (dst_text);
      return;
    }
    if (!cfp->remove_parsed) {
      dst_text = MemFree (dst_text);
    }
    avd.etp = cfp->etp;
  }
  else
  {
    avd.etp = cfp->etp;
  }
  
    
  /* put src_text in the dst_field */
  avd.field_list = cfp->dst_field_list;
  avd.new_text = src_text;
  (cfp->set_str_func) (sfp, &avd, cfp->fsp);
  
  if (cfp->convert_type == CONVERT_TYPE_SWAP
      || (cfp->convert_type == CONVERT_TYPE_PARSE && cfp->remove_parsed))
  {
    /* put dst_text in the src_field */
    avd.field_list = cfp->src_field_list;
    avd.new_text = dst_text;
    (cfp->set_str_func) (sfp, &avd, cfp->fsp);
  }
  else if (cfp->convert_type == CONVERT_TYPE_MOVE)
  {
    /* remove old src_field */
    avd.etp = NULL;
    avd.field_list = cfp->src_field_list;
    avd.new_text = NULL;
    (cfp->remove_str_func) (sfp, &avd, cfp->fsp);
  }
}

static void ConvertDescriptorFieldCallback (SeqDescrPtr sdp, Pointer userdata, FilterSetPtr fsp)
{
  ConvertFieldPtr cfp;
  CharPtr         src_text;
  CharPtr         dst_text;
  ApplyValueData  avd;
  
  if (sdp == NULL || userdata == NULL)
  {
    return;
  }
  
  cfp = (ConvertFieldPtr) userdata;
  
  if (cfp->src_field_list == NULL || cfp->dst_field_list == NULL
      || cfp->get_d_str_func == NULL || cfp->set_d_str_func == NULL
      || cfp->remove_d_str_func == NULL)
  {
    return;
  }
  
  avd.text_to_replace = NULL;
  avd.where_to_replace = EditApplyFindLocation_anywhere;
  src_text = (cfp->get_d_str_func)(sdp, cfp->src_field_list, cfp->fsp);
  if (cfp->strip_name_from_text)
  {
    StripFieldNameFromText(src_text, cfp->name_field_func, cfp->dst_field_list);
  }

  if (StringHasNoText (src_text))
  {
    src_text = MemFree (src_text);
    return;
  }
  
  if (cfp->convert_type == CONVERT_TYPE_SWAP)
  {
    /* for swap, we already know what we're doing with existing text - 
     * we're putting it where the new text used to be. */
    avd.etp = NULL;
    /* and we need to get the destination text before it's replaced */
    dst_text = (cfp->get_d_str_func)(sdp, cfp->dst_field_list, cfp->fsp);    
    if (cfp->strip_name_from_text)
    {
      StripFieldNameFromText(dst_text, cfp->name_field_func, cfp->src_field_list);
    }
    
  }
  else if (cfp->convert_type == CONVERT_TYPE_PARSE)
  {
    /* for parse, replace the source text with the portion we want */
    dst_text = src_text;
    src_text = ReplaceStringForParse(src_text, cfp->text_portion);
    if (src_text == NULL) 
    {
      dst_text = MemFree (dst_text);
      return;
    }
    if (!cfp->remove_parsed) {
      dst_text = MemFree (dst_text);
    }
    avd.etp = cfp->etp;
  }
  else
  {
    avd.etp = cfp->etp;
  }
    
  /* put src_text in the dst_field */
  avd.field_list = cfp->dst_field_list;
  avd.new_text = src_text;
  (cfp->set_d_str_func) (sdp, &avd, cfp->fsp);
  
  if (cfp->convert_type == CONVERT_TYPE_SWAP
      || (cfp->convert_type == CONVERT_TYPE_PARSE && cfp->remove_parsed))
  {
    /* put dst_text in the src_field */
    avd.etp = NULL;
    avd.field_list = cfp->src_field_list;
    avd.new_text = dst_text;
    (cfp->set_d_str_func) (sdp, &avd, cfp->fsp);
  }
  else if (cfp->convert_type == CONVERT_TYPE_MOVE)
  {
    /* remove old src_field */
    avd.etp = NULL;
    avd.field_list = cfp->src_field_list;
    avd.new_text = NULL;
    (cfp->remove_d_str_func) (sdp, &avd, cfp->fsp);
  }
}

static void ConvertNonTextFeatureFieldCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  ConvertFieldPtr cfp;
  CharPtr         curr_val = NULL;
  CharPtr         src_val;
  ApplyValueData  avd;
  
  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  
  cfp = (ConvertFieldPtr) userdata;
  
  if (cfp->src_field_list == NULL || cfp->dst_field_list == NULL
      || cfp->get_str_func == NULL || cfp->set_str_func == NULL
      || cfp->name_field_func == NULL)
  {
    return;
  }
  
  curr_val = (cfp->get_str_func)(sfp, NULL, NULL);
  src_val = (cfp->name_field_func)(cfp->src_field_list);
  if (StringCmp (curr_val, src_val) == 0)
  {
    avd.where_to_replace = EditApplyFindLocation_anywhere;
    avd.field_list = cfp->dst_field_list;
    avd.etp = NULL;
    (cfp->set_str_func) (sfp, &avd, cfp->fsp);
  }
  curr_val = MemFree (curr_val);
  src_val = MemFree (src_val);
}

static void ConvertNonTextDescriptorFieldCallback (SeqDescrPtr sdp, Pointer userdata, FilterSetPtr fsp)
{
  ConvertFieldPtr cfp;
  CharPtr         curr_val = NULL;
  CharPtr         src_val;
  ApplyValueData  avd;
  
  if (sdp == NULL || userdata == NULL)
  {
    return;
  }
  
  cfp = (ConvertFieldPtr) userdata;
  
  if (cfp->src_field_list == NULL || cfp->dst_field_list == NULL
      || cfp->get_d_str_func == NULL || cfp->set_d_str_func == NULL
      || cfp->name_field_func == NULL)
  {
    return;
  }
  
  curr_val = (cfp->get_d_str_func)(sdp, NULL, NULL);
  src_val = (cfp->name_field_func)(cfp->src_field_list);
  if (StringCmp (curr_val, src_val) == 0)
  {
    avd.where_to_replace = EditApplyFindLocation_anywhere;
    avd.field_list = cfp->dst_field_list;
    avd.etp = NULL;
    (cfp->set_d_str_func) (sdp, &avd, cfp->fsp);
  }
  curr_val = MemFree (curr_val);
  src_val = MemFree (src_val);
}

#define FEATUREFIELD_NONE        0

typedef struct featurefieldselection 
{
  DIALOG_MESSAGE_BLOCK
  PopuP                    feature_field;

  Boolean                  allow_none;
  Int4                     num_fields;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} FeatureFieldSelectionData, PNTR FeatureFieldSelectionPtr;


static void FeatureFieldSelectionChange (PopuP p)
{
  FeatureFieldSelectionPtr  dlg;

  dlg = (FeatureFieldSelectionPtr) GetObjectExtra (p);
  if (dlg == NULL)
  {
    return;
  }
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}

static void ResetFeatureFieldSelection (FeatureFieldSelectionPtr dlg)
{
  if (dlg == NULL)
  {
    return;
  }
  SetValue (dlg->feature_field, 1);
  FeatureFieldSelectionChange (dlg->feature_field);
}

static void FeatureFieldToDialog (DialoG d, Pointer data)
{
  FeatureFieldSelectionPtr  dlg;
  ValNodePtr             vnp;
  Int4                   feature_field;

  dlg = (FeatureFieldSelectionPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  vnp = (ValNodePtr) data;
  if (vnp == NULL)
  {
    ResetFeatureFieldSelection (dlg);
  }
  else
  {
    feature_field = vnp->data.intvalue;
    if (feature_field < FEATUREFIELD_NONE || feature_field > dlg->num_fields + 1)
    {
      feature_field = FEATUREFIELD_NONE;
    }
    if (dlg->allow_none)
    {
      feature_field ++;
    }
    SetValue (dlg->feature_field, feature_field);
  }  
}

static Pointer DialogToFeatureField (DialoG d)
{
  FeatureFieldSelectionPtr  dlg;
  ValNodePtr             vnp = NULL;
  Int4                   feature_field;

  dlg = (FeatureFieldSelectionPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  feature_field = GetValue (dlg->feature_field);
  if (dlg->allow_none)
  {
    feature_field --;
  }
  if (feature_field > FEATUREFIELD_NONE)
  {
    vnp = ValNodeNew (NULL);
    if (vnp != NULL)
    {
      vnp->data.intvalue = feature_field;
    }
  }
  return vnp;
}

static void FeatureFieldMessage (DialoG d, Int2 mssg)

{
  FeatureFieldSelectionPtr  dlg;

  dlg = (FeatureFieldSelectionPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        /* reset list */
        ResetFeatureFieldSelection (dlg);
        break;
      case VIB_MSG_ENTER :
        Select (dlg->feature_field);
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestFeatureFieldDialog (DialoG d)

{
  return NULL;
}

static DialoG FeatureFieldSelectionDialogEx 
(GrouP                    h,
 Boolean                  allow_none,
 CharPtr                  none_txt,
 Int4                     num_fields,
 CharPtr PNTR             field_names,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  FeatureFieldSelectionPtr  dlg;
  GrouP                     p;
  Int4                      k;
  
  dlg = (FeatureFieldSelectionPtr) MemNew (sizeof (FeatureFieldSelectionData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, 1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = FeatureFieldToDialog;
  dlg->fromdialog = DialogToFeatureField;
  dlg->dialogmessage = FeatureFieldMessage;
  dlg->testdialog = TestFeatureFieldDialog;
  dlg->allow_none = allow_none;
  dlg->num_fields = num_fields;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->feature_field = PopupList (p, TRUE, FeatureFieldSelectionChange);
  if (dlg->allow_none)
  {
    PopupItem (dlg->feature_field, none_txt);
  }
  for (k = 0; k < dlg->num_fields; k++)
  {
    PopupItem (dlg->feature_field, field_names[k]);
  }
  SetValue (dlg->feature_field, 1);
  SetObjectExtra (dlg->feature_field, dlg, NULL);
  return (DialoG) p;  
}

extern DialoG 
FeatureFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_none,
 Int4                     num_fields,
 CharPtr PNTR             field_names,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  return FeatureFieldSelectionDialogEx (h, allow_none, "None", 
                                        num_fields, field_names, 
                                        change_notify, change_userdata);
}

static DialoG FeatureFieldSelectionDialogAny
(GrouP                    h,
 Boolean                  allow_none,
 Int4                     num_fields,
 CharPtr PNTR             field_names,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  return FeatureFieldSelectionDialogEx (h, allow_none, "Any Field", 
                                        num_fields, field_names, 
                                        change_notify, change_userdata);
}


static CharPtr gene_field_list [] = 
{
  "locus", "description", "comment", "allele", "maploc", "locus_tag", "synonym", "old_locus_tag"  
};

#define GENEFIELD_LOCUS         1
#define GENEFIELD_DESCRIPTION   2
#define GENEFIELD_COMMENT       3
#define GENEFIELD_ALLELE        4
#define GENEFIELD_MAPLOC        5
#define GENEFIELD_LOCUS_TAG     6
#define GENEFIELD_SYNONYM       7
#define GENEFIELD_OLD_LOCUS_TAG 8

static int num_gene_fields = sizeof (gene_field_list) / sizeof (CharPtr);

extern DialoG 
GeneFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_none,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  return FeatureFieldSelectionDialog (h, allow_none,
                                      num_gene_fields, gene_field_list, 
                                      change_notify, change_userdata);
}

extern CharPtr 
GetGeneFieldString 
(SeqFeatPtr   sfp,
 ValNodePtr   gene_field,
 FilterSetPtr fsp)
{
  ValNodePtr vnp;
  CharPtr    str = NULL;
  GeneRefPtr grp;
  GBQualPtr  qual;
  
  if (sfp == NULL || gene_field == NULL)
  {
    return NULL;
  }
  
  if (sfp->data.choice == SEQFEAT_GENE)
  {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  }
  else
  {
    grp = SeqMgrGetGeneXref (sfp);
  }
  
  if (sfp->data.choice != SEQFEAT_GENE && grp != NULL)
  {
    vnp = NULL;
  }
  
  vnp = gene_field;
  while (vnp != NULL && vnp->data.intvalue != FEATUREFIELD_NONE && StringHasNoText (str))
  {
    str = NULL;
    switch (vnp->data.intvalue)
    {
      case GENEFIELD_LOCUS:
        if (grp != NULL)
        {
          str = grp->locus;
        }
        break;
      case GENEFIELD_DESCRIPTION:
        if (grp != NULL)
        {
          str = grp->desc;
        }
        break;
      case GENEFIELD_COMMENT:
        if (sfp->data.choice == SEQFEAT_GENE)
        {
          str = sfp->comment;
        }
        break;
      case GENEFIELD_ALLELE:
        if (grp != NULL)
        {
          str = grp->allele;
        }
        break;
      case GENEFIELD_MAPLOC:
        if (grp != NULL)
        {
          str = grp->maploc;
        }
        break;
      case GENEFIELD_SYNONYM:
        if (grp != NULL && grp->syn != NULL)
        {
          str = grp->syn->data.ptrvalue;
        }
        break;
      case GENEFIELD_LOCUS_TAG:
        if (grp != NULL)
        {
          str = grp->locus_tag;
        }
        break;
      case GENEFIELD_OLD_LOCUS_TAG:
        qual = sfp->qual;
        while (qual != NULL && StringICmp (qual->qual, "old_locus_tag") != 0)
        {
          qual = qual->next;
        }
        if (qual != NULL)
        {
          str = qual->val;
        }
        break;
    }
    vnp = vnp->next;
  }
  if (StringHasNoText (str))
  {
    str = NULL;
  }
  else
  {
    str = StringSave (str);
  }
  return str;
}

extern void RemoveGeneFieldString (SeqFeatPtr sfp, ValNodePtr gene_field)
{
  ValNodePtr vnp;
  Boolean    found_nonempty = FALSE;
  GeneRefPtr grp;
  ValNodePtr syn_remove;
  GBQualPtr  qual, prevqual;
  
  if (sfp == NULL || gene_field == NULL)
  {
    return;
  }
  
  if (sfp->data.choice == SEQFEAT_GENE)
  {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  }
  else
  {
    grp = SeqMgrGetGeneXref (sfp);
  }
  
  if (grp == NULL)
  {
    return;
  }
  
  vnp = gene_field;
  while (vnp != NULL && vnp->data.intvalue != FEATUREFIELD_NONE && !found_nonempty)
  {
    switch (vnp->data.intvalue)
    {
      case GENEFIELD_LOCUS:
        if (grp != NULL && !StringHasNoText (grp->locus))
        {
          found_nonempty = TRUE;
          if (vnp->choice != 0)
          {
            grp->locus = MemFree (grp->locus);
          }
        }
        break;
      case GENEFIELD_DESCRIPTION:
        if (grp != NULL && !StringHasNoText (grp->desc))
        {
          found_nonempty = TRUE;
          if (vnp->choice != 0)
          {
            grp->desc = MemFree (grp->desc);
          }
        }
        break;
      case GENEFIELD_COMMENT:
        if (!StringHasNoText (sfp->comment))
        {
          found_nonempty = TRUE;
          if (vnp->choice != 0)
          {
            sfp->comment = MemFree (sfp->comment);
          }
        }
        break;
      case GENEFIELD_ALLELE:
        if (grp != NULL && !StringHasNoText (grp->allele))
        {
          found_nonempty = TRUE;
          if (vnp->choice != 0)
          {
            grp->allele = MemFree (grp->allele);  
          }
        }
        break;
      case GENEFIELD_MAPLOC:
        if (grp != NULL && !StringHasNoText (grp->maploc))
        {
          found_nonempty = TRUE;
          if (vnp->choice != 0)
          {
            grp->maploc = MemFree (grp->maploc);            
          }
        }
        break;
      case GENEFIELD_SYNONYM:
        if (grp != NULL && grp->syn != NULL && !StringHasNoText (grp->syn->data.ptrvalue))
        {
          found_nonempty = TRUE;
          if (vnp->choice != 0)
          {
            syn_remove = grp->syn;
            grp->syn = grp->syn->next;
            syn_remove->next = NULL;
            ValNodeFreeData (syn_remove);
          }
        }
        break;
      case GENEFIELD_LOCUS_TAG:
        if (grp != NULL && !StringHasNoText (grp->locus_tag))
        {
          found_nonempty = TRUE;
          if (vnp->choice != 0)
          {
            grp->locus_tag = MemFree (grp->locus_tag);
          }
        }
        break;
      case GENEFIELD_OLD_LOCUS_TAG:
        prevqual = NULL;
        qual = sfp->qual;
        while (qual != NULL)
        {
          if (StringICmp (qual->qual, "old_locus_tag") == 0)
          {
            if (prevqual == NULL)
            {
              sfp->qual = qual->next;
              qual->next = NULL;
              qual = GBQualFree (qual);
              qual = sfp->qual;
            }
            else
            {
              prevqual->next = qual->next;
              qual->next = NULL;
              qual = GBQualFree (qual);
              qual = prevqual->next;
            }
          }
          else
          {
            prevqual = qual;
            qual = qual->next;
          }
        }
        break;
    }
    vnp = vnp->next;
  }
}

static void SetGeneFieldString (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  GeneRefPtr    grp;
  ApplyValuePtr avp;
  Boolean       found;
  GBQualPtr     qual, qual_last;
  
  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  avp = (ApplyValuePtr) userdata;
  if (avp->field_list == NULL)
  {
    return;
  }
  
  if (sfp->data.choice == SEQFEAT_GENE)
  {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  }
  else
  {
    grp = SeqMgrGetGeneXref (sfp);
  }
  
  if (grp == NULL && sfp->data.choice == SEQFEAT_GENE)
  {
    grp = GeneRefNew ();
    sfp->data.value.ptrvalue = grp;
  }
  
  if (grp == NULL)
  {
    return;
  }
  
  switch (avp->field_list->data.intvalue)
  {
    case GENEFIELD_LOCUS:
      grp->locus = HandleApplyValue (grp->locus, avp);
      break;
    case GENEFIELD_DESCRIPTION:
      grp->desc = HandleApplyValue (grp->desc, avp);
      break;
    case GENEFIELD_COMMENT:
      sfp->comment = HandleApplyValue (sfp->comment, avp);
      break;
    case GENEFIELD_ALLELE:
      grp->allele = HandleApplyValue (grp->allele, avp);  
      break;
    case GENEFIELD_MAPLOC:
      grp->maploc = HandleApplyValue (grp->maploc, avp);            
      break;
    case GENEFIELD_SYNONYM:
      if (grp->syn == NULL || !StringHasNoText (avp->text_to_replace)) {
        grp->syn = ApplyValueToValNodeStringList (grp->syn, 0, avp);
      } else {
        grp->syn->data.ptrvalue = HandleApplyValue (grp->syn->data.ptrvalue, avp);
      }
      break;
    case GENEFIELD_LOCUS_TAG:
      grp->locus_tag = HandleApplyValue (grp->locus_tag, avp);
      break;
    case GENEFIELD_OLD_LOCUS_TAG:
      qual = sfp->qual;
      qual_last = NULL;
      found = FALSE;
      while (qual != NULL && !found)
      {
        qual_last = qual;
        if (StringCmp (qual->qual, "old_locus_tag") == 0)
        {
          qual->val = HandleApplyValue (qual->val, avp);
          found = TRUE;
        }
        qual = qual->next;
      }
      if (!found)
      {
        qual = GBQualNew ();
        if (qual != NULL)
        {
          qual->qual = StringSave ("old_locus_tag");
          qual->val = StringSave (avp->new_text);
          if (qual_last == NULL)
          {
            sfp->qual = qual;
          }
          else
          {
            qual_last->next = qual;
          }
        }
      }    
      break;
  }
}

static CharPtr mrna_field_list [] = 
{
  "product", "comment"  
};

#define MRNAFIELD_PRODUCT     1
#define MRNAFIELD_COMMENT     2

static int num_mrna_fields = sizeof (mrna_field_list) / sizeof (CharPtr);

extern DialoG 
MRNAFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_none,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  return FeatureFieldSelectionDialog (h, allow_none,
                                      num_mrna_fields, mrna_field_list, 
                                      change_notify, change_userdata);
}

extern CharPtr 
GetmRNAFieldString 
(SeqFeatPtr   sfp,
 ValNodePtr   mrna_field,
 FilterSetPtr fsp)
{
  ValNodePtr vnp;
  CharPtr    str = NULL;
  RnaRefPtr  rrp;
  
  if (sfp == NULL || mrna_field == NULL || sfp->idx.subtype != FEATDEF_mRNA)
  {
    return NULL;
  }
  
  vnp = mrna_field;
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;

  while (vnp != NULL && vnp->data.intvalue != FEATUREFIELD_NONE && StringHasNoText (str))
  {
    str = NULL;
    switch (vnp->data.intvalue)
    {
      case MRNAFIELD_PRODUCT:
        if (rrp != NULL)
        {
          str = rrp->ext.value.ptrvalue;
        }
        break;
      case MRNAFIELD_COMMENT:
        str = sfp->comment;
        break;
    }
    vnp = vnp->next;
  }
  if (StringHasNoText (str))
  {
    str = NULL;
  }
  else
  {
    str = StringSave (str);
  }
  return str;
}

extern void RemovemRNAFieldString (SeqFeatPtr sfp, ValNodePtr mrna_field)
{
  ValNodePtr vnp;
  Boolean    found_nonempty = FALSE;
  RnaRefPtr  rrp;
  
  if (sfp == NULL || mrna_field == NULL || sfp->idx.subtype != FEATDEF_mRNA)
  {
    return;
  }
  
  vnp = mrna_field;
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  while (vnp != NULL && vnp->data.intvalue != FEATUREFIELD_NONE && !found_nonempty)
  {
    switch (vnp->data.intvalue)
    {
      case MRNAFIELD_PRODUCT:
        if (rrp != NULL && !StringHasNoText (rrp->ext.value.ptrvalue))
        {
          found_nonempty = TRUE;
          if (vnp->choice != 0)
          {
            rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
          }
        }
        break;
      case MRNAFIELD_COMMENT:
        if (!StringHasNoText (sfp->comment))
        {
          found_nonempty = TRUE;
          if (vnp->choice != 0)
          {
            sfp->comment = MemFree (sfp->comment);
          }
        }
        break;
    }
    vnp = vnp->next;
  }
}

static void SetmRNAFieldString (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  RnaRefPtr     rrp;
  ApplyValuePtr avp;
  Int4          field_choice;
  
  if (sfp == NULL || userdata == NULL || sfp->idx.subtype != FEATDEF_mRNA)
  {
    return;
  }
  
  avp = (ApplyValuePtr) userdata;
  if (avp->field_list == NULL)
  {
    return;
  }
  
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL)
  {
    rrp = RnaRefNew ();
    if (rrp == NULL)
    {
      return;
    }
    sfp->data.value.ptrvalue = rrp;
  }
  field_choice = avp->field_list->data.intvalue;
  switch (avp->field_list->data.intvalue)
  {
    case MRNAFIELD_PRODUCT:
      if (rrp->ext.choice == 0)
      {
        rrp->ext.choice = 1;
        rrp->ext.value.ptrvalue = NULL;
      }    
      rrp->ext.value.ptrvalue = HandleApplyValue (rrp->ext.value.ptrvalue, avp);
      break;
    case MRNAFIELD_COMMENT:
      sfp->comment = HandleApplyValue (sfp->comment,avp);
      break;
  }
}

static CharPtr GetCDSComment (SeqFeatPtr sfp, ValNodePtr vnp, FilterSetPtr fsp)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || StringHasNoText (sfp->comment))
  {
    return NULL;
  }
  else
  {
    return StringSave (sfp->comment);
  }
}

static void SetCDSComment (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  ApplyValuePtr avp;

  if (sfp == NULL || userdata == NULL || sfp->data.choice != SEQFEAT_CDREGION)
  {
    return;
  }
  
  avp = (ApplyValuePtr) userdata;

  sfp->comment = HandleApplyValue (sfp->comment, avp);  
}

static void RemoveCDSComment (SeqFeatPtr sfp, ValNodePtr vnp)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION)
  {
    return;
  }
  sfp->comment = MemFree (sfp->comment);
}

extern CharPtr 
GetCDSFieldString 
(SeqFeatPtr   sfp,
 ValNodePtr   cds_field,
 FilterSetPtr fsp)
{
  ValNodePtr vnp;
  CharPtr    str = NULL;
  BioseqPtr  prot_bsp;
  SeqFeatPtr prot_sfp;
  SeqMgrFeatContext prot_context;
  ProtRefPtr        prp;
  
  if (sfp == NULL || cds_field == NULL || sfp->idx.subtype != FEATDEF_CDS)
  {
    return NULL;
  }
  
  vnp = cds_field;

  while (vnp != NULL && vnp->data.intvalue != FEATUREFIELD_NONE && StringHasNoText (str))
  {
    str = NULL;
    switch (vnp->data.intvalue)
    {
      case MRNAFIELD_PRODUCT:
        prot_bsp = BioseqFind (SeqLocId (sfp->product));
        prot_sfp = SeqMgrGetNextFeature (prot_bsp, NULL, SEQFEAT_PROT, 0, &prot_context);
        if (prot_sfp != NULL && prot_sfp->data.value.ptrvalue != NULL) {
          prp = (ProtRefPtr) prot_sfp->data.value.ptrvalue;
          if (prp->name != NULL && !StringHasNoText (prp->name->data.ptrvalue)) {
            str = prp->name->data.ptrvalue;
          }
        }
        break;
      case MRNAFIELD_COMMENT:
        str = sfp->comment;
        break;
    }
    vnp = vnp->next;
  }
  if (StringHasNoText (str))
  {
    str = NULL;
  }
  else
  {
    str = StringSave (str);
  }
  return str;
}

extern void RemoveCDSFieldString (SeqFeatPtr sfp, ValNodePtr cds_field)
{
  ValNodePtr vnp;
  Boolean    found_nonempty = FALSE;
  BioseqPtr  prot_bsp;
  SeqFeatPtr prot_sfp;
  SeqMgrFeatContext prot_context;
  ProtRefPtr        prp;
  
  if (sfp == NULL || cds_field == NULL || sfp->idx.subtype != FEATDEF_CDS)
  {
    return;
  }
  
  vnp = cds_field;
  while (vnp != NULL && vnp->data.intvalue != FEATUREFIELD_NONE && !found_nonempty)
  {
    switch (vnp->data.intvalue)
    {
      case MRNAFIELD_PRODUCT:
        prot_bsp = BioseqFind (SeqLocId (sfp->product));
        prot_sfp = SeqMgrGetNextFeature (prot_bsp, NULL, SEQFEAT_PROT, 0, &prot_context);
        if (prot_sfp != NULL && prot_sfp->data.value.ptrvalue != NULL) {
          prp = (ProtRefPtr) prot_sfp->data.value.ptrvalue;
          if (prp->name != NULL && !StringHasNoText (prp->name->data.ptrvalue)) {
            found_nonempty = TRUE;
            if (vnp->choice != 0) {
              prp->name = ValNodeFreeData (prp->name);
            }
          }
        }
        break;
      case MRNAFIELD_COMMENT:
        if (!StringHasNoText (sfp->comment))
        {
          found_nonempty = TRUE;
          if (vnp->choice != 0)
          {
            sfp->comment = MemFree (sfp->comment);
          }
        }
        break;
    }
    vnp = vnp->next;
  }
}

#define PROTEINFIELD_NAME               1
#define PROTEINFIELD_DESC               2
#define PROTEINFIELD_EC_NUM             3
#define PROTEINFIELD_ACTIVITY           4
#define PROTEINFIELD_COMMENT            5
#define PROTEINFIELD_MATPEPTIDE_NAME    6
#define PROTEINFIELD_MATPEPTIDE_DESC    7
#define PROTEINFIELD_MATPEPTIDE_COMMENT 8

static CharPtr protein_field_list [] =
{
  "name", "description", "E.C. number", "activity", "comment", "mat_peptide name", "mat_peptide description", "mat_peptide comment"
};

static int num_protein_fields = sizeof (protein_field_list) / sizeof (CharPtr);

extern DialoG 
ProteinFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_none,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  return FeatureFieldSelectionDialog (h, allow_none,
                                      num_protein_fields, protein_field_list, 
                                      change_notify, change_userdata);
}

extern CharPtr GetProteinFieldString (SeqFeatPtr sfp, ValNodePtr protein_field, FilterSetPtr fsp)
{
  ValNodePtr     vnp, val_vnp;
  CharPtr        str = NULL;
  ProtRefPtr     prp = NULL;
  Int4           field_choice;
  SeqFeatXrefPtr xref;
  
  if (sfp == NULL || protein_field == NULL || (sfp->data.choice != SEQFEAT_PROT && sfp->data.choice != SEQFEAT_CDREGION))
  {
    return NULL;
  }
  
  vnp = protein_field;
  if (sfp->data.choice == SEQFEAT_PROT) {
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  } else if (sfp->data.choice == SEQFEAT_CDREGION) {
      xref = sfp->xref;
      while (xref != NULL && xref->data.choice != SEQFEAT_PROT) {
          xref = xref->next;
      }
      if (xref != NULL) {
          prp = xref->data.value.ptrvalue;
      }
  }


  while (vnp != NULL && vnp->data.intvalue != FEATUREFIELD_NONE && StringHasNoText (str))
  {
    str = NULL;
    field_choice = vnp->data.intvalue;
    switch (field_choice)
    {
      case PROTEINFIELD_NAME:
        if (prp != NULL && prp->name != NULL && (sfp->idx.subtype == FEATDEF_PROT || sfp->idx.subtype == FEATDEF_CDS))
        {
          val_vnp = prp->name;
          str = val_vnp->data.ptrvalue;
        }
        break;
      case PROTEINFIELD_DESC:
        if (prp != NULL && sfp->idx.subtype == FEATDEF_PROT)
        {
          str = prp->desc;
        }
        break;
      case PROTEINFIELD_EC_NUM:
        if (prp != NULL && prp->ec != NULL && sfp->idx.subtype == FEATDEF_PROT)
        {
          val_vnp = prp->ec;
          str = val_vnp->data.ptrvalue;
        }
        break;
      case PROTEINFIELD_ACTIVITY:
        if (prp != NULL && prp->activity != NULL && sfp->idx.subtype == FEATDEF_PROT)
        {
          val_vnp = prp->activity;
          str = val_vnp->data.ptrvalue;
        }
        break;
      case PROTEINFIELD_COMMENT:
        if (sfp->idx.subtype == FEATDEF_PROT)
        {
          str = sfp->comment;
        }
        break;
      case PROTEINFIELD_MATPEPTIDE_NAME:
        if (prp != NULL && prp->name != NULL && sfp->idx.subtype == FEATDEF_mat_peptide_aa)
        {
          val_vnp = prp->name;
          str = val_vnp->data.ptrvalue;
        }
        break;
      case PROTEINFIELD_MATPEPTIDE_DESC:
        if (sfp->idx.subtype == FEATDEF_mat_peptide_aa && prp != NULL)
        {
           str = prp->desc;
        }
        break;
      case PROTEINFIELD_MATPEPTIDE_COMMENT:
        if (sfp->idx.subtype == FEATDEF_mat_peptide_aa)
        {
          str = sfp->comment;
        }
        break;
    }
    vnp = vnp->next;
  }
  if (StringHasNoText (str))
  {
    str = NULL;
  }
  else
  {
    str = StringSave (str);
  }
  return str;
}

static void RemoveProteinFieldString (SeqFeatPtr sfp, ValNodePtr protein_field)
{
  ValNodePtr vnp;
  ProtRefPtr prp = NULL;
  Int4       field_choice;
  SeqFeatXrefPtr xref = NULL, prev;
  
  if (sfp == NULL || protein_field == NULL)
  {
    return;
  }
  
  vnp = protein_field;
  if (sfp->idx.subtype == FEATDEF_PROT
      || sfp->idx.subtype == FEATDEF_mat_peptide_aa) {
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  } else if (sfp->idx.subtype == FEATDEF_CDS) {
      xref = sfp->xref;
      while (xref != NULL && xref->data.choice != SEQFEAT_PROT) {
          xref = xref->next;
      }
      if (xref != NULL) {
          prp = xref->data.value.ptrvalue;
      }
  }      

  while (vnp != NULL && vnp->data.intvalue != FEATUREFIELD_NONE)
  {
    field_choice = vnp->data.intvalue;
    switch (field_choice)
    {
      case PROTEINFIELD_NAME:
        if (prp != NULL && prp->name != NULL && sfp->idx.subtype == FEATDEF_PROT)
        {
          prp->name = ValNodeFreeData (prp->name);
        }
        break;
      case PROTEINFIELD_DESC:
        if (prp != NULL && sfp->idx.subtype == FEATDEF_PROT)
        {
          prp->desc = MemFree (prp->desc);
        }
        break;
      case PROTEINFIELD_EC_NUM:
        if (prp != NULL && prp->ec != NULL && sfp->idx.subtype == FEATDEF_PROT)
        {
          prp->ec = ValNodeFreeData (prp->ec);
        }
        break;
      case PROTEINFIELD_ACTIVITY:
        if (prp != NULL && prp->activity != NULL && sfp->idx.subtype == FEATDEF_PROT)
        {
          prp->activity = ValNodeFreeData (prp->activity);
        }
        break;
      case PROTEINFIELD_COMMENT:
        if (sfp->idx.subtype == FEATDEF_PROT)
        {
          sfp->comment = MemFree (sfp->comment);
        }
        break;
      case PROTEINFIELD_MATPEPTIDE_NAME:
        if (prp != NULL && prp->name != NULL && sfp->idx.subtype == FEATDEF_mat_peptide_aa)
        {
          prp->name = ValNodeFreeData (prp->name);
        }
        break;
      case PROTEINFIELD_MATPEPTIDE_DESC:
        if (prp != NULL && sfp->idx.subtype == FEATDEF_mat_peptide_aa)
        {
          prp->desc = MemFree (prp->desc);
        }
        break;
      case PROTEINFIELD_MATPEPTIDE_COMMENT:
        if (sfp->idx.subtype == FEATDEF_mat_peptide_aa)
        {
          sfp->comment = MemFree (sfp->comment);
        }
        break;
    }
    vnp = vnp->next;
  }
  
  /* remove protein xref if empty */
  if (prp != NULL && xref != NULL
      && prp->name == NULL
      && StringHasNoText (prp->desc)
      && prp->ec == NULL
      && prp->activity == NULL
      && prp->db == NULL) {
    prp = ProtRefFree(prp);
    xref->data.value.ptrvalue = NULL;
    if (sfp->xref == xref) {
        sfp->xref = xref->next;
    } else {
        prev = sfp->xref;
        while (prev != NULL && prev->next != xref) {
            prev = prev->next;
        }
        if (prev != NULL) {
            prev->next = xref->next;
        }
    }
    xref = SeqFeatXrefFree(xref);
  }      
}

extern ValNodePtr ApplyValueToValNodeStringListAsText (ValNodePtr list, Int2 choice, ApplyValuePtr avp);

static void SetProteinFieldString (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  ProtRefPtr    prp;
  ApplyValuePtr avp;
  
  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  if (sfp->idx.subtype != FEATDEF_PROT && sfp->idx.subtype != FEATDEF_mat_peptide_aa)
  {
    return;
  }
  
  avp = (ApplyValuePtr) userdata;
  if (avp->field_list == NULL)
  {
    return;
  }

  if (sfp->idx.subtype == FEATDEF_PROT)
  {
    if (avp->field_list->data.intvalue == PROTEINFIELD_MATPEPTIDE_NAME
        || avp->field_list->data.intvalue == PROTEINFIELD_MATPEPTIDE_DESC
        || avp->field_list->data.intvalue == PROTEINFIELD_MATPEPTIDE_COMMENT)
    {
      return;
    }
  }
  else if (sfp->idx.subtype == FEATDEF_mat_peptide_aa)
  {
    if (avp->field_list->data.intvalue != PROTEINFIELD_MATPEPTIDE_NAME
        && avp->field_list->data.intvalue != PROTEINFIELD_MATPEPTIDE_DESC
        && avp->field_list->data.intvalue != PROTEINFIELD_MATPEPTIDE_COMMENT)
    {
      return;
    }
  }
  else    
  {
    return; 
  }
  
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL)
  {
    prp = ProtRefNew ();
    if (prp == NULL)
    {
      return;
    }
    sfp->data.value.ptrvalue = prp;
  }
  
  switch (avp->field_list->data.intvalue)
  {
    case PROTEINFIELD_NAME:
    case PROTEINFIELD_MATPEPTIDE_NAME:
      if (prp->name == NULL || !StringHasNoText(avp->text_to_replace))
      {
        prp->name = ApplyValueToValNodeStringList (prp->name, 0, avp);
      }
      else
      {
        prp->name->data.ptrvalue = HandleApplyValue (prp->name->data.ptrvalue, avp);
      }
      break;
    case PROTEINFIELD_DESC:
    case PROTEINFIELD_MATPEPTIDE_DESC:
      prp->desc = HandleApplyValue (prp->desc, avp);
      break;
    case PROTEINFIELD_EC_NUM:
      prp->ec = ApplyValueToValNodeStringListAsText (prp->ec, 0, avp);
      break;
    case PROTEINFIELD_ACTIVITY:
      if (prp->activity == NULL || !StringHasNoText(avp->text_to_replace))
      {
        prp->activity = ApplyValueToValNodeStringList (prp->activity, 0, avp);
      }
      else
      {
        prp->activity->data.ptrvalue = HandleApplyValue (prp->activity->data.ptrvalue, avp);
      }
      break;
    case PROTEINFIELD_COMMENT:
    case PROTEINFIELD_MATPEPTIDE_COMMENT:
      sfp->comment = HandleApplyValue (sfp->comment, avp);
      break;
  }
}

static CharPtr PNTR BuildCDSGeneFieldList (Int4Ptr num_fields)
{
  CharPtr PNTR field_name_list;
  Int4         i, k;
  Char         tmp[100];

  if (num_fields == NULL) 
  {
    return NULL;
  }
  *num_fields = 1 + num_gene_fields + num_mrna_fields + num_protein_fields;
  field_name_list = (CharPtr PNTR) MemNew (*num_fields * sizeof (CharPtr));
  if (field_name_list == NULL)
  {
    return NULL;
  }
  
  field_name_list [0] = StringSave ("CDS comment");
  k = 1;
  for (i = 0; i < num_gene_fields; i++)
  {
    sprintf (tmp, "Gene %s", gene_field_list [i]);
    field_name_list [k++] = StringSave (tmp);
  }
  for (i = 0; i < num_mrna_fields; i++)
  {
    sprintf (tmp, "mRNA %s", mrna_field_list [i]);
    field_name_list [k++] = StringSave (tmp);
  }
  for (i = 0; i < num_protein_fields; i++)
  {
    if (StringNCmp (protein_field_list [i], "mat_peptide", 11) == 0)
    {
      sprintf (tmp, "%s", protein_field_list [i]);
      tmp [0] = toupper (tmp [0]);
    }
    else
    {
      sprintf (tmp, "Protein %s", protein_field_list [i]);
    }
    field_name_list [k++] = StringSave (tmp);
  }
  return field_name_list;  
}

static void FreeCDSGeneFieldList (CharPtr PNTR field_name_list, Int4 num_fields)
{
  Int4 i;
  
  if (field_name_list == NULL || num_fields == 0)
  {
    return;
  }
  
  for (i = 0; i < num_fields; i++)
  {
    field_name_list [i] = MemFree (field_name_list [i]);
  }
  field_name_list = MemFree (field_name_list);
}

extern DialoG 
CDSGeneProtFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_none,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  CharPtr PNTR field_name_list;
  Int4         num_fields = 0;
  DialoG       d;

  field_name_list = BuildCDSGeneFieldList (&num_fields);
  
  d = FeatureFieldSelectionDialog (h, allow_none,
                                      num_fields, field_name_list, 
                                      change_notify, change_userdata);
  FreeCDSGeneFieldList (field_name_list, num_fields);
  return d;
}

static DialoG
CDSGeneProtFieldConstraintSelectionDialog
(GrouP                    h,
 Boolean                  allow_none,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  CharPtr PNTR field_name_list;
  Int4         num_fields = 0;
  DialoG       d;

  field_name_list = BuildCDSGeneFieldList (&num_fields);
  
  d = FeatureFieldSelectionDialogAny (h, allow_none,
                                      num_fields, field_name_list, 
                                      change_notify, change_userdata);
  FreeCDSGeneFieldList (field_name_list, num_fields);
  return d; 
}


extern Boolean IsCDSetProteinProductChoice (ValNodePtr vnp)
{
  if (vnp != NULL && vnp->data.intvalue == 2 + num_gene_fields + num_mrna_fields) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean IsCDSetProteinQualChoice (ValNodePtr vnp)
{
  if (vnp != NULL 
      && vnp->data.intvalue > 1 + num_gene_fields + num_mrna_fields
      && vnp->data.intvalue <= 1 + num_gene_fields + num_mrna_fields + num_protein_fields) 
  {
    return TRUE;
  } 
  else 
  {
    return FALSE;
  }
}

static Boolean IsCDSetMatPeptideQualChoice (ValNodePtr vnp)
{
  if (IsCDSetProteinQualChoice (vnp)
      && StringNICmp (protein_field_list[vnp->data.intvalue - num_gene_fields - num_mrna_fields - 2],
                      "mat_peptide", 11) == 0) 
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  } 
}

static Boolean IsCDSetMRNAQualChoice (ValNodePtr vnp)
{
  if (vnp != NULL 
      && vnp->data.intvalue > 1 + num_gene_fields
      && vnp->data.intvalue <= 1 + num_gene_fields + num_mrna_fields) 
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

static Boolean IsCDSetGeneQualChoice (ValNodePtr vnp)
{
  if (vnp != NULL 
      && vnp->data.intvalue > 1
      && vnp->data.intvalue <= 1 + num_gene_fields)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

static Boolean IsCDSetCDSQualChoice (ValNodePtr vnp)
{
  if (vnp != NULL && vnp->data.intvalue == 1) 
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

extern CharPtr GetCDSGeneProtField (SeqFeatPtr sfp, ValNodePtr vnp, FilterSetPtr fsp)
{
  CharPtr str = NULL;
  ValNode vn;
  Int4    field_choice;
  
  if (sfp == NULL || vnp == NULL)
  {
    return NULL;
  }
  
  while (vnp != NULL && str == NULL)
  {
    vn.choice = vnp->choice;
    vn.next = NULL;
    
    field_choice = vnp->data.intvalue;
  
    if (IsCDSetCDSQualChoice(vnp))
    {
      /* CDS Comment */
      str = GetCDSComment (sfp, NULL, NULL);
    }
    else if (IsCDSetGeneQualChoice (vnp))
    {
      vn.data.intvalue = vnp->data.intvalue - 1;
      str = GetGeneFieldString (sfp, &vn, NULL);
    }
    else if (IsCDSetMRNAQualChoice (vnp))
    {
      vn.data.intvalue = vnp->data.intvalue - num_gene_fields - 1;
      str = GetmRNAFieldString (sfp, &vn, NULL);
    }
    else if (IsCDSetProteinQualChoice (vnp))
    {
      vn.data.intvalue = vnp->data.intvalue - num_gene_fields - num_mrna_fields - 1;
      str = GetProteinFieldString (sfp, &vn, NULL);
    }
    vnp = vnp->next;
  }
  return str;
}

/* we could have multiple mat_peptides in a CDSset, some that match the constraint and some that do not */
static Boolean DoesConstraintDisqualifyFeature (SeqFeatPtr sfp, ChoiceConstraintPtr ccp)
{
  Boolean is_disqualified = FALSE;
  Boolean does_match = FALSE;
  CharPtr str;

  if (sfp == NULL || ccp == NULL 
      || ccp->constraint_type == CHOICE_CONSTRAINT_ANY
      || sfp->idx.subtype != FEATDEF_mat_peptide_aa 
      || !IsCDSetMatPeptideQualChoice(ccp->qual_choice))
  {
    is_disqualified = FALSE;
  }
  else if (ccp->constraint_type == CHOICE_CONSTRAINT_QUAL_PRESENT && ccp->qual_choice != NULL)
  {
    str = GetCDSGeneProtField (sfp, ccp->qual_choice, NULL);
    if (str == NULL)
    {
      is_disqualified = TRUE;
    }
    MemFree (str);
  }
  else if (ccp->constraint_type == CHOICE_CONSTRAINT_STRING)
  {
    str = GetCDSGeneProtField (sfp, ccp->qual_choice, NULL);
    does_match = DoesStringMatchConstraintX (str, ccp->string_constraint);
    MemFree (str);
    if (ccp->string_constraint != NULL && ccp->string_constraint->not_present)
    {
      does_match = ! does_match;
    }
    if (!does_match)
    {
      is_disqualified = TRUE;
    }
  }
  return is_disqualified;

}

extern void RemoveCDSGeneProtField (SeqFeatPtr sfp, ValNodePtr vnp, FilterSetPtr fsp)
{
  ValNode vn;
  
  if (sfp == NULL || vnp == NULL)
  {
    return;
  }

  if (fsp != NULL && fsp->cgp != NULL && DoesConstraintDisqualifyFeature (sfp, fsp->cgp)) 
  {
    return;
  }
  
  while (vnp != NULL)
  {
    vn.choice = vnp->choice;
    vn.next = NULL;
  
    if (vnp->data.intvalue == 1)
    {
      /* CDS Comment */
      RemoveCDSComment (sfp, NULL);
    }
    else if (vnp->data.intvalue <= 1 + num_gene_fields)
    {
      vn.data.intvalue = vnp->data.intvalue - 1;
      RemoveGeneFieldString (sfp, &vn);
    }
    else if (vnp->data.intvalue <= 1 + num_gene_fields + num_mrna_fields)
    {
      vn.data.intvalue = vnp->data.intvalue - num_gene_fields - 1;
      RemovemRNAFieldString (sfp, &vn);
    }
    else if (vnp->data.intvalue <= 1 + num_gene_fields + num_mrna_fields + num_protein_fields)
    {
      vn.data.intvalue = vnp->data.intvalue - num_gene_fields - num_mrna_fields - 1;
      RemoveProteinFieldString (sfp, &vn);
    }
    vnp = vnp->next;
  }
}


extern Uint2 FeatDefTypeFromFieldList (ValNodePtr vnp)
{
  if (vnp == NULL) 
  {
    return 0;
  }
  else if (IsCDSetCDSQualChoice(vnp))
  {
    return FEATDEF_CDS;
  }
  else if (IsCDSetGeneQualChoice (vnp))
  {
    return FEATDEF_GENE;
  }
  else if (IsCDSetMRNAQualChoice(vnp))
  {
    return FEATDEF_mRNA;
  }
  else if (IsCDSetProteinQualChoice (vnp))
  {
    if (IsCDSetMatPeptideQualChoice (vnp))
    {
      return FEATDEF_mat_peptide_aa;
    }
    else
    {
      return FEATDEF_PROT;
    }
  }
  else
  {
    return 0;
  }
}


/* return TRUE if actually set field, FALSE otherwise */
extern Boolean 
SetCDSGeneProtField 
(SeqFeatPtr      sfp,
 ValNodePtr      vnp, 
 ApplyValuePtr   avp,
 FilterSetPtr    fsp)
{
  ValNode        vn;
  ApplyValueData local_avd;
  Uint2          choice;
  
  if (sfp == NULL || vnp == NULL || avp == NULL)
  {
    return FALSE;
  }

  choice = FeatDefTypeFromFieldList (vnp);
  if (FindFeatFromFeatDefType(choice) != sfp->data.choice)
  {
    return FALSE;
  }

  if (fsp != NULL && fsp->cgp != NULL && DoesConstraintDisqualifyFeature (sfp, fsp->cgp)) 
  {
    return FALSE;
  }
  
  vn.choice = vnp->choice;
  vn.next = NULL;
  local_avd.field_list = &vn;
  local_avd.new_text = avp->new_text;
  local_avd.etp = avp->etp;
  local_avd.text_to_replace = avp->text_to_replace;
  local_avd.where_to_replace = avp->where_to_replace;

  if (vnp->data.intvalue == 1)
  {
    /* CDS Comment */
    SetCDSComment (sfp, &local_avd, NULL);
  }
  else if (vnp->data.intvalue <= 1 + num_gene_fields)
  {
    vn.data.intvalue = vnp->data.intvalue - 1;
    SetGeneFieldString (sfp, &local_avd, NULL);
  }
  else if (vnp->data.intvalue <= 1 + num_gene_fields + num_mrna_fields)
  {
    vn.data.intvalue = vnp->data.intvalue - num_gene_fields - 1;
    SetmRNAFieldString (sfp, &local_avd, NULL);
  }
  else if (vnp->data.intvalue <= 1 + num_gene_fields + num_mrna_fields + num_protein_fields)
  {
    vn.data.intvalue = vnp->data.intvalue - num_gene_fields - num_mrna_fields - 1;
    SetProteinFieldString (sfp, &local_avd, NULL);
  }
  
  /* no need to free any part of local_avd */

  return TRUE;
}

#define RNAFIELD_PRODUCT           1
#define RNAFIELD_COMMENT           2
#define RNAFIELD_CODONS_RECOGNIZED 3
#define RNAFIELD_NCRNA_CLASS       4
#define RNAFIELD_ANTICODON         5
#define RNAFIELD_TRANSCRIPT_ID     6

static CharPtr rna_field_list [] = 
{
  "Product", "Comment", "Codons Recognized", "ncRNA class", "Anticodon", "Transcript ID"
};

static Int4 num_rna_fields = sizeof (rna_field_list) / sizeof (CharPtr);

static CharPtr GetRNAFieldName (ValNodePtr vnp)
{
  CharPtr label = NULL;
  
  if (vnp != NULL && vnp->data.intvalue >= 1)
  {
    if (vnp->data.intvalue <= num_rna_fields)
    {
      label = StringSave (rna_field_list [vnp->data.intvalue - 1]);
    }
    else if (vnp->data.intvalue <= num_rna_fields + num_gene_fields)
    {
      label = (CharPtr) MemNew ((StringLen (gene_field_list [vnp->data.intvalue - 1 - num_rna_fields]) + 6) * sizeof (Char));
      if (label != NULL)
      {
        sprintf (label, "Gene %s", gene_field_list [vnp->data.intvalue - 1 - num_rna_fields]);
      }
    }
  }
  return label;
}

extern DialoG
RNAAddFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  ValNodePtr choice_list = NULL;
  DialoG     dlg;
  Int4       i;
  
  for (i = 0; i < 4; i++)
  {
    ValNodeAddInt (&choice_list, 0, i + 1);
  }
  for (i = 0; i < num_gene_fields; i++)
  {
    ValNodeAddInt (&choice_list, 0, i + num_rna_fields + 1);
  }

  dlg = ValNodeSelectionDialog (h, choice_list, TALL_SELECTION_LIST,
                                GetRNAFieldName, NULL, 
                                IntValNodeCopy, IntValNodeMatch,
                                "RNA field", change_notify, change_userdata, 
                                FALSE);
                                
  return dlg;
}

extern DialoG
RNARemoveFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  ValNodePtr choice_list = NULL;
  DialoG     dlg;
  Int4       i;
  
  for (i = 0; i < num_rna_fields; i++)
  {
    ValNodeAddInt (&choice_list, 0, i + 1);
  }
  for (i = 0; i < num_gene_fields; i++)
  {
    ValNodeAddInt (&choice_list, 0, i + num_rna_fields + 1);
  }

  dlg = ValNodeSelectionDialog (h, choice_list, TALL_SELECTION_LIST,
                                GetRNAFieldName, NULL, 
                                IntValNodeCopy, IntValNodeMatch,
                                "RNA field", change_notify, change_userdata, 
                                FALSE);

  return dlg;
}

extern DialoG
RNAFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  ValNodePtr choice_list = NULL;
  DialoG     dlg;
  Int4       i;
  
  for (i = 0; i < 4; i++)
  {
    ValNodeAddInt (&choice_list, 0, i + 1);
  }
  for (i = 0; i < num_gene_fields; i++)
  {
    ValNodeAddInt (&choice_list, 0, i + num_rna_fields + 1);
  }

  dlg = ValNodeSelectionDialog (h, choice_list, TALL_SELECTION_LIST,
                                GetRNAFieldName, NULL, 
                                IntValNodeCopy, IntValNodeMatch,
                                "RNA field", change_notify, change_userdata, 
                                FALSE);
                                
  return dlg;
}

static CharPtr rna_subtype_list [] = { "misc_RNA", "preRna", "mRNA", "tRNA",
                                        "rRNA", "ncRNA", "tmRNA"};
static Int2    num_rna_subtypes = sizeof (rna_subtype_list) / sizeof (CharPtr); 

static CharPtr GetRNASubtypeName (ValNodePtr vnp)
{
  CharPtr label = NULL;
  
  if (vnp != NULL && vnp->choice < num_rna_subtypes)
  {
    label = StringSave (rna_subtype_list[vnp->choice]);
  }
  return label;
}

static DialoG RNASubtypeSelectionDialog
(GrouP                    h,
 Boolean                  allow_multi,
 SeqEntryPtr              sep,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  ValNodePtr choice_list = NULL;
  DialoG     dlg;
  
  ValNodeAddInt (&choice_list, 0, FEATDEF_otherRNA);
  ValNodeAddInt (&choice_list, 1, FEATDEF_preRNA);
  ValNodeAddInt (&choice_list, 2, FEATDEF_mRNA);
  ValNodeAddInt (&choice_list, 3, FEATDEF_tRNA);
  ValNodeAddInt (&choice_list, 4, FEATDEF_rRNA);
  ValNodeAddInt (&choice_list, 5, FEATDEF_ncRNA);
  ValNodeAddInt (&choice_list, 6, FEATDEF_tmRNA);

  dlg = ValNodeSelectionDialog (h, choice_list, TALL_SELECTION_LIST,
                                GetRNASubtypeName, NULL, 
                                IntValNodeCopy, IntValNodeMatch,
                                "RNA subtype", change_notify, change_userdata, 
                                allow_multi);
  return dlg;
}


extern CharPtr GetRNAFieldString (SeqFeatPtr sfp, ValNodePtr vnp, FilterSetPtr fsp)
{
  RnaRefPtr  rrp;
  CharPtr    str = NULL;
  ValNode    vn;
  GeneRefPtr grp;
  SeqFeatPtr gene;
  SeqMgrFeatContext context;
  Int4              field_choice;
  GBQualPtr         qual;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || vnp == NULL)
  {
    return NULL;
  }
  
  rrp = sfp->data.value.ptrvalue;
  
  while (vnp != NULL && str == NULL)
  {
    field_choice = vnp->data.intvalue;
    switch (vnp->data.intvalue)
    {
      case RNAFIELD_PRODUCT :        
        str = GetRNAProductString (sfp, NULL);
        break;
      case RNAFIELD_COMMENT :
        if (!StringHasNoText (sfp->comment))
        {
          str = StringSave (sfp->comment);
        }
        break;
      case RNAFIELD_NCRNA_CLASS :
        for (qual = sfp->qual; qual != NULL && str == NULL; qual = qual->next)
        {
          if (StringCmp (qual->qual, "ncRNA_class") == 0
              && !StringHasNoText (qual->val))
          {
            str = StringSave (qual->val);
          }
        }
        break;
      default:
        if (vnp->data.intvalue >= num_rna_fields + 1)
        {
          vn.choice = 0;
          vn.next = NULL;
          vn.data.intvalue = vnp->data.intvalue - num_rna_fields;
          grp = SeqMgrGetGeneXref (sfp);
          if (grp == NULL)
          {
            gene = SeqMgrGetOverlappingGene (sfp->location, &context);
            str = GetGeneFieldString (gene, &vn, NULL);
          }
          else
          {
            str = GetGeneFieldString (sfp, &vn, NULL);
          }
        }
        break;
          
    }
    vnp = vnp->next;
  }
  return str;
}


static Boolean IsParseabletRNAName (CharPtr name_string)
{
  if (StringHasNoText(name_string)) 
  {
    return TRUE;
  }
  else if (StringNICmp (name_string, "trna-", 5) != 0)
  {
    return FALSE;
  }
  else if (StringLen (name_string) != 8)
  {
    return FALSE;
  }
  else if (ParseTRnaString (name_string, NULL, NULL, TRUE) == 0)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


static void SettRNAName (SeqFeatPtr sfp, ApplyValuePtr avp)
{
  RnaRefPtr     rrp;
  tRNAPtr           trp;
  Uint1             new_aa;
  CharPtr           name_string;
  Boolean           justTrnaText = FALSE;
  Uint1             codon [6];
  GBQualPtr         gbqual, gbqual_last = NULL;

  if (sfp == NULL || sfp->data.value.ptrvalue == NULL 
      || sfp->idx.subtype != FEATDEF_tRNA
      || avp == NULL)
  {
    return;
  }
  
  rrp = sfp->data.value.ptrvalue;
  if (rrp->ext.choice == 0)
  {
    trp = MemNew (sizeof (tRNA));
    trp->aatype = 0;
    MemSet (trp->codon, 255, sizeof (trp->codon));
    trp->anticodon = NULL;
    rrp->ext.value.ptrvalue = trp;
    rrp->ext.choice = 2;
  }
  if (rrp->ext.choice != 2)
  {
    return;
  }
  trp = (tRNAPtr) rrp->ext.value.ptrvalue;
  
  name_string = GetRNAProductString (sfp, NULL);
  name_string = HandleApplyValue (name_string, avp);
  if (!IsParseabletRNAName(name_string))
  {
    if (trp->anticodon == NULL
        && trp->codon[0] == 255
        && trp->codon[1] == 255
        && trp->codon[2] == 255
        && trp->codon[3] == 255
        && trp->codon[4] == 255
        && trp->codon[5] == 255)
    {
      trp = MemFree (trp);
      rrp->ext.choice = 1;
      rrp->ext.value.ptrvalue = name_string;
    }
    else
    {
      trp->aa = 0;
      gbqual = sfp->qual;
      while (gbqual != NULL)
      {
        gbqual_last = gbqual;
        gbqual = gbqual->next;
      }
      gbqual = GBQualNew ();
      gbqual->qual = StringSave ("product");
      gbqual->val = name_string;
      if (gbqual_last == NULL)
      {
        sfp->qual = gbqual;
      }
      else
      {
        gbqual->next = gbqual_last->next;
        gbqual_last->next = gbqual;
      }
    }
  }
  else
  {
    new_aa = ParseTRnaString (name_string, &justTrnaText, codon, TRUE);
    trp->aa = new_aa;
    trp->aatype = 2;
    name_string = MemFree (name_string);
  }
}


static void SetRNAProduct (SeqFeatPtr sfp, ApplyValuePtr avp)
{
  RnaRefPtr rrp;
  GBQualPtr gbq, prev_gbq = NULL;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL || avp == NULL)
  {
    return;
  }

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;

  gbq = sfp->qual;
  while (gbq != NULL && StringCmp (gbq->qual, "product") != 0)
  {
    prev_gbq = gbq;
    gbq = gbq->next;
  }
  if (gbq != NULL)
  {
    gbq->val = HandleApplyValue (gbq->val, avp);
    if (StringHasNoText (gbq->val))
    {
      if (prev_gbq == NULL)
      {
        sfp->qual = gbq->next;
      }
      else
      {
        prev_gbq->next = gbq->next;
      }
      gbq->next = NULL;
      gbq = GBQualFree (gbq);
    }
  }
  else if (sfp->idx.subtype == FEATDEF_tRNA || rrp->ext.choice == 2)
  {
    SettRNAName (sfp, avp);
  }
  else if (rrp->ext.choice == 1
           && (StringCmp (rrp->ext.value.ptrvalue, "ncRNA") == 0
               || StringCmp (rrp->ext.value.ptrvalue, "tmRNA") == 0
               || StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") == 0))
  {
    gbq = GBQualNew();
    gbq->qual = StringSave ("product");
    gbq->val = HandleApplyValue (gbq->val, avp);
    gbq->next = sfp->qual;
    sfp->qual = gbq;
  }
  else if (rrp->ext.choice == 1 || rrp->ext.choice == 0)
  {
    rrp->ext.value.ptrvalue = HandleApplyValue (rrp->ext.value.ptrvalue, avp);
    if (StringHasNoText (rrp->ext.value.ptrvalue))
    {
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
      rrp->ext.choice = 0;
    }
    else
    {
      rrp->ext.choice = 1;
    }
  }
}


static void SetRNAFieldString (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  RnaRefPtr     rrp;
  ApplyValuePtr avp;
  SeqFeatPtr    gene;
  GeneRefPtr    grp;
  SeqMgrFeatContext context;
  GBQualPtr         qual;
  Boolean           found = FALSE;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || userdata == NULL)
  {
    return;
  }
  
  avp = (ApplyValuePtr) userdata;
  if (avp->field_list == NULL)
  {
    return;
  }

  rrp = sfp->data.value.ptrvalue;
    
  switch (avp->field_list->data.intvalue)
  {
    case RNAFIELD_PRODUCT :
      SetRNAProduct (sfp, avp);
      break;
    case RNAFIELD_COMMENT :
      sfp->comment = HandleApplyValue (sfp->comment, avp);
      break;
    case RNAFIELD_NCRNA_CLASS :
      for (qual = sfp->qual; qual != NULL; qual = qual->next)
      {
        if (StringCmp (qual->qual, "ncRNA_class") == 0)
        {
          qual->val = HandleApplyValue (qual->val, avp);
          found = TRUE;
        }
      }
      if (!found)
      {
        qual = GBQualNew ();
        qual->qual = StringSave ("ncRNA_class");
        qual->val = StringSave (avp->new_text);
        qual->next = sfp->qual;
        sfp->qual = qual;
      }
      break;
    default:
      if (avp->field_list->data.intvalue >= num_rna_fields + 1)
      {
        avp->field_list->data.intvalue -= num_rna_fields;
        grp = SeqMgrGetGeneXref (sfp);
        if (grp == NULL)
        {
          gene = SeqMgrGetOverlappingGene (sfp->location, &context);
          SetGeneFieldString (gene, avp, fsp);
        }
        else
        {
          SetGeneFieldString (sfp, avp, fsp);
        }
        avp->field_list->data.intvalue += num_rna_fields;
      }
      break;
  }
}


static void RemoveRNAProduct (SeqFeatPtr sfp)
{
  RnaRefPtr rrp;
  GBQualPtr gbq, prev_gbq = NULL;
  tRNAPtr     trp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL)
  {
    return;
  }

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;

  gbq = sfp->qual;
  while (gbq != NULL && StringCmp (gbq->qual, "product") != 0)
  {
    prev_gbq = gbq;
    gbq = gbq->next;
  }
  if (gbq != NULL)
  {
    if (prev_gbq == NULL)
    {
      sfp->qual = gbq->next;
    }
    else
    {
      prev_gbq->next = gbq->next;
    }
    gbq->next = NULL;
    gbq = GBQualFree (gbq);
  }
  else if (rrp->ext.choice == 1
           && (StringCmp (rrp->ext.value.ptrvalue, "ncRNA") == 0
               || StringCmp (rrp->ext.value.ptrvalue, "tmRNA") == 0
               || StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") == 0))
  {
    /* do nothing - already have no product */
  }
  else if (rrp->ext.choice == 0)
  {
    /* do nothing - already have no product */
  }
  else if (rrp->ext.choice == 1)
  {
    rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
    rrp->ext.choice = 0;
  }
  else if (rrp->ext.choice == 2)
  {
    trp = (tRNAPtr) rrp->ext.value.ptrvalue;
    if (trp != NULL)
    {
      trp->aatype = 0;
    }
  }
}


static void RemoveRNAField (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  RnaRefPtr     rrp;
  tRNAPtr       trp;
  ApplyValuePtr avp;
  SeqFeatPtr    gene;
  GeneRefPtr    grp;
  ValNode       vn;
  SeqMgrFeatContext context;
  GBQualPtr         qual, qual_prev = NULL, qual_next;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || userdata == NULL)
  {
    return;
  }
  
  avp = (ApplyValuePtr) userdata;
  if (avp->field_list == NULL)
  {
    return;
  }
  
  rrp = sfp->data.value.ptrvalue;

  switch (avp->field_list->data.intvalue)
  {
    case RNAFIELD_PRODUCT :
      RemoveRNAProduct (sfp);
      break;
    case RNAFIELD_COMMENT :
      sfp->comment = MemFree (sfp->comment);
      break;
    case RNAFIELD_NCRNA_CLASS :
      qual = sfp->qual;
      while (qual != NULL)
      {
        qual_next = qual->next;
        if (StringCmp (qual->qual, "ncRNA_class") == 0)
        {
          if (qual_prev == NULL)
          {
            sfp->qual = qual_next;
          }
          else
          {
            qual_prev->next = qual_next;
          }
          qual->next = NULL;
          qual = GBQualFree (qual);
        }
        else
        {
          qual_prev = qual;
        }
        qual = qual_next;
      }
      break;
    case RNAFIELD_ANTICODON :
      if (rrp->ext.choice == 2) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        if (trp != NULL && trp->anticodon != NULL) {
          trp->anticodon = SeqLocFree (trp->anticodon);
        }
      }
      break;
    case RNAFIELD_CODONS_RECOGNIZED :
      if (rrp->ext.choice == 2) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        if (trp != NULL) {
          trp->codon [0] = 255;
          trp->codon [1] = 255;
          trp->codon [2] = 255;
          trp->codon [3] = 255;
          trp->codon [4] = 255;
          trp->codon [5] = 255;
        }
      }
      break;
    case RNAFIELD_TRANSCRIPT_ID :
      sfp->product = SeqLocFree (sfp->product);
      break;
    default:
      if (avp->field_list->data.intvalue >= num_rna_fields + 1)
      {
        vn.choice = 1;
        vn.next = NULL;
        vn.data.intvalue = avp->field_list->data.intvalue - num_rna_fields;
        grp = SeqMgrGetGeneXref (sfp);
        if (grp == NULL)
        {
          gene = SeqMgrGetOverlappingGene (sfp->location, &context);
          RemoveGeneFieldString (gene, &vn);
        }
        else
        {
          RemoveGeneFieldString (sfp, &vn);
        }
      }
      break;
  }
}

typedef struct exonfieldselection 
{
  DIALOG_MESSAGE_BLOCK
  PopuP                    exon_field;

  Boolean                  allow_none;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;  
} ExonFieldSelectionData, PNTR ExonFieldSelectionPtr;

#define EXONFIELD_ALLELE            1
#define EXONFIELD_COMMENT           2
#define EXONFIELD_EC_NUMBER         3
#define EXONFIELD_FUNCTION          4
#define EXONFIELD_OLD_LOCUS_TAG     5
#define EXONFIELD_NUMBER            6
#define EXONFIELD_PRODUCT           7

static CharPtr exon_field_list [] = 
{
  "allele", "comment", "EC_number", "function", "old_locus_tag", "number", "product"  
};
#define NUM_EXON_FIELDS 7

extern DialoG
ExonFieldSelectionDialog
(GrouP                    h,
 Boolean                  allow_none,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  return FeatureFieldSelectionDialog (h, allow_none,
                                      NUM_EXON_FIELDS, exon_field_list, 
                                      change_notify, change_userdata);

}

extern CharPtr GetExonFieldString (SeqFeatPtr sfp, ValNodePtr exon_field)
{
  ValNodePtr vnp;
  CharPtr    str = NULL;
  GBQualPtr  gbqual;
  
  if (sfp == NULL || exon_field == NULL || sfp->idx.subtype != FEATDEF_exon)
  {
    return NULL;
  }
  
  vnp = exon_field;
  while (vnp != NULL && StringHasNoText (str))
  {
    str = NULL;
    if (vnp->data.intvalue == EXONFIELD_COMMENT)
    {
      str = sfp->comment;
    }
    else
    {
      gbqual = sfp->qual;
      while (gbqual != NULL && StringHasNoText (str))
      {
        if (vnp->data.intvalue < NUM_EXON_FIELDS + 1
            && vnp->data.intvalue >= 1
            && StringCmp (gbqual->qual, exon_field_list[vnp->data.intvalue - 1]) == 0)
        {
          str = gbqual->val;
        }
        gbqual = gbqual->next;
      }    
    }
    vnp = vnp->next;
  }
  if (StringHasNoText (str))
  {
    str = NULL;
  }
  else
  {
    str = StringSave (str);
  }
  return str;
}

extern void RemoveExonFieldString (SeqFeatPtr sfp, ValNodePtr exon_field)
{
  ValNodePtr vnp;
  Boolean    found_nonempty = FALSE;
  GBQualPtr  gbqual, prev_qual;
  
  if (sfp == NULL || exon_field == NULL || sfp->idx.subtype != FEATDEF_exon)
  {
    return;
  }
  
  vnp = exon_field;
  while (vnp != NULL && !found_nonempty)
  {
    if (vnp->data.intvalue == EXONFIELD_COMMENT)
    {
      if (!StringHasNoText (sfp->comment))
      {
        found_nonempty = TRUE;
        if (vnp->choice != 0)
        {
          sfp->comment = MemFree (sfp->comment);
        }
      }
    }
    else
    {     
      gbqual = sfp->qual;
      prev_qual = NULL;
      while (gbqual != NULL)
      {
        if (vnp->data.intvalue <= NUM_EXON_FIELDS + 1
            && vnp->data.intvalue >= 1
            && StringCmp (gbqual->qual, exon_field_list[vnp->data.intvalue - 1]) == 0)
        {
           if (!StringHasNoText (gbqual->val))
           {
            found_nonempty = TRUE;
            if (vnp->choice != 0)
            {
              if (prev_qual == NULL)
              {
                sfp->qual = gbqual->next;
              } else {
                prev_qual->next = gbqual->next;
              }
              gbqual->next = NULL;
              GBQualFree (gbqual);
              return;
            }
          }
        }
        prev_qual = gbqual;
        gbqual = gbqual->next;
      }
    }
    vnp = vnp->next;
  }
}


#define PARSE_FIELD_BIOSRC_TAXNAME  1
#define PARSE_FIELD_BIOSRC_LINEAGE  2
#define PARSE_FIELD_BIOSRC_DIVISION 3

static CharPtr biosrc_string_list [] =
{
  "Organism Name", "Lineage", "Division"
};

static int num_biosrc_strings = sizeof (biosrc_string_list) / sizeof (CharPtr);

extern DialoG BioSourceStringDialog 
(GrouP                    h,
 Boolean                  allow_multi,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  return FeatureFieldSelectionDialog (h, allow_multi,
                                      num_biosrc_strings, biosrc_string_list, 
                                      change_notify, change_userdata);
}

static CharPtr GetSourceStringFromBioSource (BioSourcePtr biop, ValNodePtr vnp)
{
  if (biop == NULL || vnp == NULL || biop->org == NULL)
  {
    return NULL;
  }
  
  if (vnp->data.intvalue == PARSE_FIELD_BIOSRC_TAXNAME
      && !StringHasNoText (biop->org->taxname))
  {
    return StringSave (biop->org->taxname);
  }
  else if (vnp->data.intvalue == PARSE_FIELD_BIOSRC_LINEAGE
           && biop->org->orgname != NULL
           && !StringHasNoText (biop->org->orgname->lineage))
  {
    return StringSave (biop->org->orgname->lineage);
  }
  else if (vnp->data.intvalue == PARSE_FIELD_BIOSRC_DIVISION
           && biop->org->orgname != NULL
           && !StringHasNoText (biop->org->orgname->div))
  {
    return StringSave (biop->org->orgname->div);
  }
  else
  {
    return NULL;
  }
}

static CharPtr 
GetSourceFeatureString 
(SeqFeatPtr   sfp,
 ValNodePtr   vnp,
 FilterSetPtr fsp)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || vnp == NULL)
  {
    return NULL;
  }
  
  return GetSourceStringFromBioSource (sfp->data.value.ptrvalue, vnp);
}

static CharPtr 
GetSourceDescriptorString 
(SeqDescrPtr  sdp,
 ValNodePtr   vnp,
 FilterSetPtr fsp)
{
  if (sdp == NULL || sdp->choice != Seq_descr_source || vnp == NULL)
  {
    return NULL;
  }
  return GetSourceStringFromBioSource (sdp->data.ptrvalue, vnp);
}


static CharPtr GetDbxrefFromBioSource (BioSourcePtr biop, ValNodePtr vnp)
{
  ValNodePtr dbx;
  DbtagPtr   dbtag;
  Char       buf[20];

  if (biop == NULL || vnp == NULL || vnp->data.ptrvalue == NULL || biop->org == NULL)
  {
    return NULL;
  }
  
  dbx = biop->org->db;
  while (dbx != NULL)
  {
    dbtag = (DbtagPtr) dbx->data.ptrvalue;
    if (dbtag != NULL && dbtag->tag != NULL
        && StringCmp (dbtag->db, vnp->data.ptrvalue) == 0)
    {
      if (dbtag->tag->str == NULL)
      {
        sprintf (buf, "%d", dbtag->tag->id);
        return StringSave (buf);
      }
      else
      {
        return StringSave (dbtag->tag->str);
      }
    }
    dbx = dbx->next;
  }
  return NULL;
}

static CharPtr 
GetBioSourceFeatureDbxrefString 
(SeqFeatPtr   sfp,
 ValNodePtr   vnp,
 FilterSetPtr fsp)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || vnp == NULL)
  {
    return NULL;
  }
  
  return GetDbxrefFromBioSource (sfp->data.value.ptrvalue, vnp);
}

static CharPtr 
GetBioSourceDescriptorDbxrefString 
(SeqDescrPtr  sdp,
 ValNodePtr   vnp,
 FilterSetPtr fsp)
{
  if (sdp == NULL || sdp->choice != Seq_descr_source || vnp == NULL)
  {
    return NULL;
  }
  return GetDbxrefFromBioSource (sdp->data.ptrvalue, vnp);
}


#define TEXT_PORTION_START_AFTER 1
#define TEXT_PORTION_START_WITH  2
#define TEXT_PORTION_END_AFTER   1
#define TEXT_PORTION_END_WITH    2

typedef struct textportiondialog
{
  DIALOG_MESSAGE_BLOCK

  GrouP  start_choice;
  TexT   start_text;
  GrouP  end_choice;
  TexT   end_text;
  ButtoN rem_before;
  ButtoN also_rem_before;
  ButtoN rem_after;
  ButtoN also_rem_after;
  ButtoN insensitive;
  ButtoN whole_word;

  Boolean inside;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} TextPortionXDialogData, PNTR TextPortionXDialogPtr;


static void OutsideEnableDisable (TextPortionXDialogPtr tp)
{
  if (tp == NULL || tp->inside)
  {
    return;
  }
  if (GetStatus (tp->rem_before)) 
  {
    Enable (tp->start_text);
    Enable (tp->also_rem_before);
  }
  else
  {
    Disable (tp->start_text);
    Disable (tp->also_rem_before);
  }
  if (GetStatus (tp->rem_after)) 
  {
    Enable (tp->end_text);
    Enable (tp->also_rem_after);
  }
  else
  {
    Disable (tp->end_text);
    Disable (tp->also_rem_after);
  }
}

static void ResetTextPortionXDialog (TextPortionXDialogPtr tp)
{
  if (tp == NULL)
  {
    return;
  }
  if (tp->inside)
  {
    SetValue (tp->start_choice, TEXT_PORTION_START_AFTER);
    SetValue (tp->end_choice, TEXT_PORTION_END_AFTER);
  }
  else
  {
    SetStatus (tp->rem_before, FALSE);
    SetStatus (tp->also_rem_before, FALSE);
    SetStatus (tp->rem_after, FALSE);
    SetStatus (tp->also_rem_after, FALSE);
  }
  SetTitle (tp->start_text, "");
  SetTitle (tp->end_text, "");
  SetStatus (tp->insensitive, FALSE);
  SetStatus (tp->whole_word, FALSE);
  OutsideEnableDisable (tp);
}

static void TextPortionXToDialog (DialoG d, Pointer data)
{
  TextPortionXDialogPtr tdlg;
  TextPortionXPtr       tdata;
  Int4                 start_choice, end_choice;
  
  tdlg = (TextPortionXDialogPtr) GetObjectExtra (d);
  if (tdlg == NULL)
  {
    return;
  }
  tdata = (TextPortionXPtr) data;
  ResetTextPortionXDialog (tdlg);  
  if (tdata != NULL)
  {
    start_choice = tdata->end_choice;
    end_choice = tdata->end_choice;
    if (start_choice < TEXT_PORTION_START_AFTER || start_choice > TEXT_PORTION_START_WITH)
    {
      start_choice = TEXT_PORTION_START_AFTER;
    }
    if (end_choice < TEXT_PORTION_END_AFTER || end_choice > TEXT_PORTION_END_WITH)
    {
      end_choice = TEXT_PORTION_END_AFTER;
    }
    if (tdlg->inside)
    {
      SetValue (tdlg->start_choice, start_choice);
      SetValue (tdlg->end_choice, end_choice);
    }
    else
    {
      SetStatus (tdlg->rem_before, StringHasNoText (tdata->start_text));
      SetStatus (tdlg->rem_after, StringHasNoText (tdata->end_text));
      if (start_choice == TEXT_PORTION_START_AFTER)
      {
        SetStatus (tdlg->also_rem_before, TRUE);
      }
      else
      {
        SetStatus (tdlg->also_rem_before, FALSE);
      }
      if (end_choice == TEXT_PORTION_END_AFTER)
      {
        SetStatus (tdlg->also_rem_after, TRUE);
      }
      else
      {
        SetStatus (tdlg->also_rem_after, FALSE);
      }
    }
    if (tdata->start_text != NULL)
    {
      SetTitle (tdlg->start_text, tdata->start_text);
    }
    if (tdata->end_text != NULL)
    {
      SetTitle (tdlg->end_text, tdata->end_text);
    }
    SetStatus (tdlg->insensitive, tdata->insensitive);
    SetStatus (tdlg->whole_word, tdata->whole_word);
  }
  OutsideEnableDisable (tdlg);
}

static Pointer DialogToTextPortionX (DialoG d)
{
  TextPortionXDialogPtr tdlg;
  TextPortionXPtr       tdata;

  tdlg = (TextPortionXDialogPtr) GetObjectExtra (d);
  if (tdlg == NULL)
  {
    return NULL;
  }

  tdata = (TextPortionXPtr) MemNew (sizeof (TextPortionXData));
  if (tdata != NULL)
  {
    if (tdlg->inside)
    {
      tdata->start_choice = GetValue (tdlg->start_choice);
      tdata->end_choice = GetValue (tdlg->end_choice);
      tdata->start_text = JustSaveStringFromText (tdlg->start_text);
      tdata->end_text = JustSaveStringFromText (tdlg->end_text);
    }
    else
    {
      if (GetStatus (tdlg->rem_before))
      {
        tdata->start_text = JustSaveStringFromText (tdlg->start_text);
      }
      else
      {
        tdata->start_text = NULL;
      }
      if (GetStatus (tdlg->rem_after))
      {
        tdata->end_text = JustSaveStringFromText (tdlg->end_text);
      }
      else
      {
        tdata->end_text = NULL;
      }

      if (GetStatus (tdlg->also_rem_before))
      {
        tdata->start_choice = TEXT_PORTION_START_AFTER;
      }
      else
      {
        tdata->start_choice = TEXT_PORTION_START_WITH;
      }
      if (GetStatus (tdlg->also_rem_after))
      {
        tdata->end_choice = TEXT_PORTION_END_AFTER;
      }
      else
      {
        tdata->end_choice = TEXT_PORTION_END_WITH;
      }

    }
        
    tdata->insensitive = GetStatus (tdlg->insensitive);
    tdata->whole_word = GetStatus (tdlg->whole_word);
  }
  return tdata;
}

static void TextPortionXMessage (DialoG d, Int2 mssg)

{
  TextPortionXDialogPtr tp;

  tp = (TextPortionXDialogPtr) GetObjectExtra (d);
  if (tp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        ResetTextPortionXDialog (tp);        
        break;
      case VIB_MSG_ENTER :
        Select (tp->start_text);
        break;
      case NUM_VIB_MSG + 1 :
        SetTitle (tp->start_text, "");
        SetTitle (tp->end_text, "");
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestTextPortionXDialog (DialoG d)

{
  return NULL;
}

static void ChangeTextPortionXGroup (GrouP g)
{
  TextPortionXDialogPtr tp;

  tp = (TextPortionXDialogPtr) GetObjectExtra (g);
  if (tp == NULL) return;

  if (tp->change_notify != NULL)
  {
    (tp->change_notify) (tp->change_userdata);
  }
}

static void ChangeTextPortionXText (TexT t)
{
  TextPortionXDialogPtr tp;

  tp = (TextPortionXDialogPtr) GetObjectExtra (t);
  if (tp == NULL) return;

  if (tp->change_notify != NULL)
  {
    (tp->change_notify) (tp->change_userdata);
  }
}


static void ChangeTextPortionXBtn (ButtoN b)
{
  TextPortionXDialogPtr tp;

  tp = (TextPortionXDialogPtr) GetObjectExtra (b);
  if (tp == NULL) return;
  OutsideEnableDisable (tp);

  if (tp->change_notify != NULL)
  {
    (tp->change_notify) (tp->change_userdata);
  }
}



extern DialoG TextPortionXDialogEx (GrouP h, Boolean inside, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  TextPortionXDialogPtr tp;
  GrouP                p, g1, g2;
  
  tp = (TextPortionXDialogPtr) MemNew (sizeof (TextPortionXDialogData));
  if (tp == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, tp, StdCleanupExtraProc);

  tp->dialog = (DialoG) p;
  tp->todialog = TextPortionXToDialog;
  tp->fromdialog = DialogToTextPortionX;
  tp->dialogmessage = TextPortionXMessage;
  tp->testdialog = TestTextPortionXDialog;

  tp->change_notify = change_notify;
  tp->change_userdata = change_userdata;
  tp->inside = inside;

  g1 = HiddenGroup (p, 3, 0, NULL);
  SetGroupSpacing (g1, 10, 10);

  if (inside) 
  {
    StaticPrompt (g1, "Between", 0, popupMenuHeight, programFont, 'r');
    tp->start_choice = HiddenGroup (g1, 2, 0, ChangeTextPortionXGroup);
    RadioButton (tp->start_choice, "just after");
    RadioButton (tp->start_choice, "starting at");
    SetValue (tp->start_choice, TEXT_PORTION_START_AFTER);
    SetObjectExtra (tp->start_choice, tp, NULL);

    tp->start_text = DialogText (g1, "", 10, ChangeTextPortionXText);
    SetObjectExtra (tp->start_text, tp, NULL);
    
    StaticPrompt (g1, "And", 0, popupMenuHeight, programFont, 'r');
    tp->end_choice = HiddenGroup (g1, 2, 0, ChangeTextPortionXGroup);
    RadioButton (tp->end_choice, "up to");
    RadioButton (tp->end_choice, "including");
    SetValue (tp->end_choice, TEXT_PORTION_END_AFTER);
    SetObjectExtra (tp->end_choice, tp, NULL);
      
    tp->end_text = DialogText (g1, "", 10, ChangeTextPortionXText);
    SetObjectExtra (tp->end_text, tp, NULL);
  }
  else
  {
    tp->rem_before = CheckBox (g1, "Before", ChangeTextPortionXBtn);
    SetObjectExtra (tp->rem_before, tp, NULL);
    tp->start_text = DialogText (g1, "", 10, ChangeTextPortionXText);
    SetObjectExtra (tp->start_text, tp, NULL);
    tp->also_rem_before = CheckBox (g1, "Also Remove Entered Text", ChangeTextPortionXBtn);

    tp->rem_after = CheckBox (g1, "After", ChangeTextPortionXBtn);
    SetObjectExtra (tp->rem_after, tp, NULL);
    tp->end_text = DialogText (g1, "", 10, ChangeTextPortionXText);
    SetObjectExtra (tp->end_text, tp, NULL);
    tp->also_rem_after = CheckBox (g1, "Also Remove Entered Text", ChangeTextPortionXBtn);    
  }
  
  g2 = HiddenGroup (p, 2, 0, NULL);
  tp->insensitive = CheckBox (g2, "Case insensitive", ChangeTextPortionXBtn);
  SetObjectExtra (tp->insensitive, tp, NULL);
  tp->whole_word = CheckBox (g2, "Whole word", ChangeTextPortionXBtn);
  SetObjectExtra (tp->whole_word, tp, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) g2, NULL);
  
  ResetTextPortionXDialog (tp);
  return (DialoG) p;
}

extern DialoG TextPortionXDialog (GrouP h)
{
  return TextPortionXDialogEx (h, TRUE, NULL, NULL);
}

extern TextPortionXPtr TextPortionXFree (TextPortionXPtr tp)
{
  if (tp == NULL) return NULL;
  tp->start_text = MemFree (tp->start_text);
  tp->end_text = MemFree (tp->end_text);
  tp = MemFree (tp->end_text);
  return tp;
}


static TextPortionXPtr TextPortionXNew (void)
{
  TextPortionXPtr tp;

  tp = (TextPortionXPtr) MemNew (sizeof (TextPortionXData));
  return tp;
}


static Boolean IsWholeWordMatch (CharPtr start, CharPtr found, Int4 match_len)
{
  Boolean rval = TRUE;
  Char    char_after;
  Char    char_before;
  
  if (match_len == 0)
  {
    rval = TRUE;
  }
  else if (start == NULL || found == NULL)
  {
    rval = FALSE;
  }
  else
  {
      char_after = *(found + match_len);
    if (found != start)
      {
        char_before = *(found - 1);
        if (isalpha ((Int4) char_before) || isdigit ((Int4) char_before))
        {
          rval = FALSE;
        }
      }
      if (char_after != 0 && (isalpha ((Int4) char_after) || isdigit ((Int4)char_after)))
      {
        rval = FALSE;
      }   
  }
  return rval;
}

extern void 
FindTextPortionXInString 
(CharPtr        str, 
 TextPortionXPtr tp, 
 CharPtr PNTR   ploc, 
 Int4Ptr        plen)
{
  CharPtr found_start, found_end;
  Int4    found_len;
  
  if (ploc == NULL || plen == NULL || tp == NULL)
  {
    return;
  }
  *ploc = NULL;
  *plen = 0;

  if (str == NULL)
  {
    return;
  }
  
  if (tp->start_text == NULL || tp->start_text [0] == 0)
  {
    found_start = str;
  }
  else
  {
    if (tp->insensitive)
    {
      found_start = StringISearch (str, tp->start_text);
    }
    else
    {
      found_start = StringSearch (str, tp->start_text);
    }
    
    if (tp->whole_word && ! IsWholeWordMatch (str, found_start, StringLen (tp->start_text)))
    {
      found_start = NULL;
    }
  }
  
  if (found_start == NULL)
  {
    return;
  }
  

  
  if (tp->start_choice == TEXT_PORTION_START_AFTER)
  {
    found_start += StringLen (tp->start_text);
  }
  
  if (tp->end_text == NULL || tp->end_text [0] == 0)
  {
    found_len = StringLen (found_start);
  }
  else
  {
    if (tp->insensitive)
    {
      found_end = StringISearch (found_start, tp->end_text);
    }
    else
    {
      found_end = StringSearch (found_start, tp->end_text);
    }
    if (tp->whole_word && ! IsWholeWordMatch (str, found_end, StringLen (tp->end_text)))
    {
      found_end = NULL;
    }    
    
    if (found_end == NULL)
    {
      found_len = 0;
    }
    else if (tp->end_choice == TEXT_PORTION_END_WITH)
    {
      found_len = (Int4)(found_end - found_start) + StringLen (tp->end_text);
    }
    else
    {
      found_len = found_end - found_start;
    }
  }

  if (found_len > 0)
  {
    *ploc = found_start;
    *plen = found_len;
  }
}

static void ClearDialogBtn (ButtoN b)
{
  DialoG d;
  
  d = (DialoG) GetObjectExtra (b);
  
  PointerToDialog (d, NULL);
}

typedef struct changecasedlg {
  DIALOG_MESSAGE_BLOCK
  GrouP change_case;
} ChangeCaseDlgData, PNTR ChangeCaseDlgPtr;

typedef enum {
  eChangeCaseNone = 1,
  eChangeCaseAllLower,
  eChangeCaseCapFirst
} EChangeCaseRule;

typedef struct changecase {
  EChangeCaseRule change;
} ChangeCaseData, PNTR ChangeCasePtr;

static void DataToChangeCaseDialog (DialoG d, Pointer data)
{
  ChangeCasePtr    ccp = (ChangeCasePtr) data;
  ChangeCaseDlgPtr dlg = (ChangeCaseDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  if (ccp == NULL || ccp->change == eChangeCaseNone) {
    SetValue (dlg->change_case, eChangeCaseNone);
  } else {
    SetValue (dlg->change_case, ccp->change);
  }
}

static Pointer ChangeCaseDialogToData (DialoG d)
{
  ChangeCaseDlgPtr dlg = (ChangeCaseDlgPtr) GetObjectExtra (d);
  ChangeCasePtr    ccp = (ChangeCasePtr) MemNew (sizeof (ChangeCaseData));

  if (dlg == NULL) {
    ccp->change = eChangeCaseNone;
  } else {
    ccp->change = GetValue (dlg->change_case);
  }
  return ccp;
}

static DialoG ChangeCaseDialog (GrouP h)
{
  ChangeCaseDlgPtr dlg;
  GrouP            p;
  
  dlg = (ChangeCaseDlgPtr) MemNew (sizeof (ChangeCaseDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = DataToChangeCaseDialog;
  dlg->fromdialog = ChangeCaseDialogToData;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;

  dlg->change_case = NormalGroup (p, 3, 0, "Capitalization", programFont, NULL);
  SetGroupSpacing (dlg->change_case, 10, 10);

  RadioButton (dlg->change_case, "No change");
  RadioButton (dlg->change_case, "To lower");
  RadioButton (dlg->change_case, "First cap, rest lower");
  SetValue (dlg->change_case, eChangeCaseNone);
  
  return (DialoG) p;
}


static void ChangeCase (CharPtr PNTR pText, ChangeCasePtr ccp, ValNodePtr orgnames)
{
  if (ccp != NULL) {
    switch (ccp->change) {
      case eChangeCaseNone:
        break;
      case eChangeCaseAllLower:
        FixCapitalizationInTitle (pText, FALSE, orgnames);
        break;
      case eChangeCaseCapFirst:
        FixCapitalizationInTitle (pText, TRUE, orgnames);
        break;
    }
  }
}

typedef struct stringconstraintdialog 
{
  DIALOG_MESSAGE_BLOCK
  PopuP  match_choice;
  TexT   match_text;
  ButtoN insensitive;
  ButtoN whole_word;   
} StringConstraintDialogXData, PNTR StringConstraintDialogXPtr;

static void ResetStringConstraintDialogX (StringConstraintDialogXPtr scdp)
{
  if (scdp == NULL) return;
  
  SetValue (scdp->match_choice, 1);
  SetTitle (scdp->match_text, "");
  SetStatus (scdp->insensitive, FALSE);
  SetStatus (scdp->whole_word, FALSE);
}


static Int4 GetPopupPosForStringConstraint (StringConstraintXPtr scp)
{
  Int4 rval = 1;

  if (scp == NULL) return 1;

  switch (scp->match_location)
  {
    case eStringConstraintContains:
      if (scp->not_present) 
      {
        rval = 2;
      } else {
        rval = 1;
      }
      break;
    case eStringConstraintEquals:
      if (scp->not_present) {
        rval = 4;
      } else {
        rval = 3;
      }
      break;
    case eStringConstraintStarts:
      if (scp->not_present) {
        rval = 9;
      } else {
        rval = 5;
      }
      break;
    case eStringConstraintEnds:
      if (scp->not_present) {
        rval = 10;
      } else {
        rval = 6;
      }
      break;
    case eStringConstraintInList:
      if (scp->not_present) {
        rval = 8;
      } else {
        rval = 7;
      }
      break;
  }
  return rval;
}


static void StringConstraintToDialog (DialoG d, Pointer data)

{
  StringConstraintDialogXPtr scdp;
  StringConstraintXPtr       scp;

  scdp = (StringConstraintDialogXPtr) GetObjectExtra (d);
  scp = (StringConstraintXPtr) data;
  if (scdp == NULL)
  {
    return;
  }

  if (scp == NULL)
  {
    ResetStringConstraintDialogX (scdp);
  }
  else
  {
    SetValue (scdp->match_choice, GetPopupPosForStringConstraint (scp));
    
    if (StringHasNoText (scp->match_text))
    {
      SetTitle (scdp->match_text, " ");
    }
    else
    {
      SetTitle (scdp->match_text, scp->match_text);
    }
    
    SetStatus (scdp->insensitive, scp->insensitive);
    SetStatus (scdp->whole_word, scp->whole_word);
  }
}

static Pointer DialogToStringConstraint (DialoG d)

{
  StringConstraintDialogXPtr scdp;
  StringConstraintXPtr       scp;
  Int4                      match_choice;

  scdp = (StringConstraintDialogXPtr) GetObjectExtra (d);
  if (scdp == NULL)
  {
    return NULL;
  }
  scp = (StringConstraintXPtr) MemNew (sizeof (StringConstraintData));
  if (scp != NULL)
  {
    scp->match_text = SaveStringFromText (scdp->match_text);
    scp->insensitive = GetStatus (scdp->insensitive);
    scp->whole_word = GetStatus (scdp->whole_word);
    match_choice = GetValue (scdp->match_choice);
    switch (match_choice)
    {
      case 1:
        scp->match_location = eStringConstraintContains;
        scp->not_present = FALSE;
        break;
      case 2:
        scp->match_location = eStringConstraintContains;
        scp->not_present = TRUE;
        break;
      case 3:
        scp->match_location = eStringConstraintEquals;
        scp->not_present = FALSE;
        break;
      case 4:
        scp->match_location = eStringConstraintEquals;
        scp->not_present = TRUE;
        break;
      case 5:
        scp->match_location = eStringConstraintStarts;
        scp->not_present = FALSE;
        break;
      case 6:
        scp->match_location = eStringConstraintEnds;
        scp->not_present = FALSE;
        break;
      case 7:
        scp->match_location = eStringConstraintInList;
        scp->not_present = FALSE;
        break;
      case 8:
        scp->match_location = eStringConstraintInList;
        scp->not_present = TRUE;
        break;
      case 9:
        scp->match_location = eStringConstraintStarts;
        scp->not_present = TRUE;
        break;
      case 10:
        scp->match_location = eStringConstraintEnds;
        scp->not_present = TRUE;
        break;
      default:
        scp->match_location = eStringConstraintContains;
        scp->not_present = FALSE;
        break;
    }
  }
  return scp;
}

static void StringConstraintMessage (DialoG d, Int2 mssg)

{
  StringConstraintDialogXPtr scdp;

  scdp = (StringConstraintDialogXPtr) GetObjectExtra (d);
  if (scdp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        ResetStringConstraintDialogX (scdp);        
        break;
      case VIB_MSG_ENTER :
        Select (scdp->match_text);
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestStringConstraintDialogX (DialoG d)

{
  return NULL;
}

extern DialoG StringConstraintDialogX (GrouP h, CharPtr label, Boolean clear_btn)

{
  StringConstraintDialogXPtr scdp;
  GrouP                     p, g, k;
  ButtoN                    b = NULL;
  
  scdp = (StringConstraintDialogXPtr) MemNew (sizeof (StringConstraintDialogXData));
  if (scdp == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, scdp, StdCleanupExtraProc);

  scdp->dialog = (DialoG) p;
  scdp->todialog = StringConstraintToDialog;
  scdp->fromdialog = DialogToStringConstraint;
  scdp->dialogmessage = StringConstraintMessage;
  scdp->testdialog = TestStringConstraintDialogX;

  g = HiddenGroup (p, 3, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  
  if (!StringHasNoText (label))
  {
    StaticPrompt (g, label, 0, dialogTextHeight, systemFont, 'l');
  }
  
  scdp->match_choice = PopupList (g, TRUE, NULL);
  PopupItem (scdp->match_choice, "Contains");
  PopupItem (scdp->match_choice, "Does not contain");
  PopupItem (scdp->match_choice, "Equals");
  PopupItem (scdp->match_choice, "Does not equal");
  PopupItem (scdp->match_choice, "Starts with");
  PopupItem (scdp->match_choice, "Ends with");
  PopupItem (scdp->match_choice, "Is one of");
  PopupItem (scdp->match_choice, "Is not one of");
  PopupItem (scdp->match_choice, "Does not start with");
  PopupItem (scdp->match_choice, "Does not end with");
  SetValue (scdp->match_choice, 1);
  scdp->match_text = DialogText (g, "", 15, NULL);
  
  k = HiddenGroup (p, 3, 0, NULL);
  SetGroupSpacing (k, 10, 10);
  scdp->insensitive = CheckBox (k, "Case Insensitive", NULL);
  scdp->whole_word = CheckBox (k, "Whole Word", NULL);
  
  if (clear_btn)
  {
    b = PushButton (p, "Clear Constraint", ClearDialogBtn);
    SetObjectExtra (b, p, NULL);
  }

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) k, (HANDLE) b, NULL);
    
  return (DialoG) p;
}


typedef struct pseudoconstraintdialog 
{
  DIALOG_MESSAGE_BLOCK
  GrouP  pseudo_choice;  
  DialoG feature_choice;
} PseudoConstraintDialogData, PNTR PseudoConstraintDialogPtr;

static ENUM_ALIST(pseudoconstraint_alist)
{"Any Feature",   FEATDEF_ANY  },
{"CDS",   FEATDEF_CDS  },
{"Gene",  FEATDEF_GENE },
{"mRNA",  FEATDEF_mRNA },
END_ENUM_ALIST

static void ResetPseudoConstraintDialog (PseudoConstraintDialogPtr pcdp)
{
  if (pcdp != NULL)
  {
    SetValue (pcdp->pseudo_choice, 1);
    SendMessageToDialog (pcdp->feature_choice, VIB_MSG_INIT);
  }
}

static void PseudoConstraintToDialog (DialoG d, Pointer data)
{
  PseudoConstraintDialogPtr pcdp;
  PseudoConstraintPtr       pcp;
  ValNode                   vn;

  pcdp = (PseudoConstraintDialogPtr) GetObjectExtra (d);
  
  if (pcdp != NULL)
  {
    ResetPseudoConstraintDialog (pcdp);
    pcp = (PseudoConstraintPtr) data;
    if (pcp != NULL)
    {
      if (pcp->is_pseudo)
      {
        SetValue (pcdp->pseudo_choice, 1);
      }
      else
      {
        SetValue (pcdp->pseudo_choice, 2);
      }
      vn.choice = pcp->featdef_type;
      vn.data.ptrvalue = NULL;
      vn.next = NULL;
      PointerToDialog (pcdp->feature_choice, &vn);
    }
  }
}

static Pointer DialogToPseudoConstraint (DialoG d)
{
  PseudoConstraintDialogPtr pcdp;
  PseudoConstraintPtr       pcp;  
  ValNodePtr                vnp;

  pcdp = (PseudoConstraintDialogPtr) GetObjectExtra (d);
  if (pcdp == NULL) return NULL;
  pcp = (PseudoConstraintPtr) MemNew (sizeof (PseudoConstraintData));
  if (pcp != NULL)
  {
    if (GetValue (pcdp->pseudo_choice) == 1)
    {
      pcp->is_pseudo = TRUE;
    }
    else
    {
      pcp->is_pseudo = FALSE;
    }
    vnp = DialogToPointer (pcdp->feature_choice);
    if (vnp == NULL)
    {
      pcp->featdef_type = FEATDEF_ANY;
    }
    else
    {
      pcp->featdef_type = vnp->choice;
    }
    ValNodeFreeData (vnp);
  }
  return pcp;
}

static void PseudoConstraintMessage (DialoG d, Int2 mssg)

{
  PseudoConstraintDialogPtr pcdp;

  pcdp = (PseudoConstraintDialogPtr) GetObjectExtra (d);
  if (pcdp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        ResetPseudoConstraintDialog (pcdp);        
        break;
      case VIB_MSG_ENTER :
        Select (pcdp->feature_choice);
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestPseudoConstraintDialog (DialoG d)

{
  return NULL;
}

static DialoG PseudoConstraintDialog (GrouP h)
{
  PseudoConstraintDialogPtr pcdp;
  GrouP                     p;
  
  pcdp = (PseudoConstraintDialogPtr) MemNew (sizeof (PseudoConstraintDialogData));
  if (pcdp == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, 2, 0, NULL);
  SetObjectExtra (p, pcdp, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  pcdp->dialog = (DialoG) p;
  pcdp->todialog = PseudoConstraintToDialog;
  pcdp->fromdialog = DialogToPseudoConstraint;
  pcdp->dialogmessage = PseudoConstraintMessage;
  pcdp->testdialog = TestPseudoConstraintDialog;

  pcdp->feature_choice = EnumAssocSelectionDialog (p, pseudoconstraint_alist,
                                                   "feature", FALSE, NULL, NULL);

  pcdp->pseudo_choice = HiddenGroup (p, 2, 0, NULL);
  RadioButton (pcdp->pseudo_choice, "Is Pseudo");
  RadioButton (pcdp->pseudo_choice, "Is Not Pseudo");
  SetValue (pcdp->pseudo_choice, 1);
  
  return (DialoG) p;  
}

typedef struct choiceconstraintdialog 
{
  DIALOG_MESSAGE_BLOCK
  GrouP               constraint_type;
  DialoG              present_qual_choice;
  DialoG              string_qual_choice;
  DialoG              string_constraint;
  DialoG              pseudo_constraint;
  DialoG              match_1;
  DialoG              match_2;
  FreeValNodeProc     free_vn_proc;
  CopyValNodeDataProc copy_vn_proc;
} ChoiceConstraintDialogData, PNTR ChoiceConstraintDialogPtr;

static void ChangeChoiceConstraintType (GrouP p)
{
  ChoiceConstraintDialogPtr scdp;
  Int4                      constraint_type;

  scdp = (ChoiceConstraintDialogPtr) GetObjectExtra (p);
  if (scdp == NULL)
  {
    return;
  }
  
  constraint_type = GetValue (scdp->constraint_type);
  if (constraint_type == CHOICE_CONSTRAINT_QUAL_PRESENT)
  {
    Enable (scdp->present_qual_choice);
    Disable (scdp->string_qual_choice);
    Disable (scdp->string_constraint);
    Disable (scdp->match_1);
    Disable (scdp->match_2);
    SafeDisable (scdp->pseudo_constraint);
  }
  else if (constraint_type == CHOICE_CONSTRAINT_STRING)
  {
    Enable (scdp->string_qual_choice);
    Enable (scdp->string_constraint);
    Disable (scdp->present_qual_choice);
    Disable (scdp->match_1);
    Disable (scdp->match_2);
    SafeDisable (scdp->pseudo_constraint);
  }
  else if (constraint_type == CHOICE_CONSTRAINT_MATCH)
  {
    Enable (scdp->match_1);
    Enable (scdp->match_2);
    Disable (scdp->string_qual_choice);
    Disable (scdp->string_constraint);
    Disable (scdp->present_qual_choice);
    SafeDisable (scdp->pseudo_constraint);
  }
  else if (constraint_type == CHOICE_CONSTRAINT_PSEUDO)
  {
    Disable (scdp->present_qual_choice);
    Disable (scdp->string_qual_choice);
    Disable (scdp->string_constraint);
    Disable (scdp->match_1);
    Disable (scdp->match_2);
    Enable (scdp->pseudo_constraint);
  }
  else
  {
    Disable (scdp->present_qual_choice);
    Disable (scdp->string_qual_choice);
    Disable (scdp->string_constraint);
    Disable (scdp->match_1);
    Disable (scdp->match_2);
    SafeDisable (scdp->pseudo_constraint);
  }
}

static void ResetChoiceConstraintDialog (ChoiceConstraintDialogPtr scdp)
{
  if (scdp != NULL)  
  {
    SetValue (scdp->constraint_type, CHOICE_CONSTRAINT_ANY);
    PointerToDialog (scdp->present_qual_choice, NULL);
    PointerToDialog (scdp->string_qual_choice, NULL);
    PointerToDialog (scdp->string_constraint, NULL);
    PointerToDialog (scdp->pseudo_constraint, NULL);
    ChangeChoiceConstraintType (scdp->constraint_type);
  }
}

static void ChoiceConstraintToDialog (DialoG d, Pointer data)
{
  ChoiceConstraintDialogPtr src_dlg;
  ChoiceConstraintPtr       src_data;
  Int4 constraint_type;
  
  src_dlg = (ChoiceConstraintDialogPtr) GetObjectExtra (d);
  if (src_dlg == NULL)
  {
    return;
  }
  src_data = (ChoiceConstraintPtr) data;
  if (src_data == NULL)
  {
    ResetChoiceConstraintDialog (src_dlg);
  }
  else
  {
    constraint_type = src_data->constraint_type;
    if (constraint_type < CHOICE_CONSTRAINT_ANY || constraint_type > CHOICE_CONSTRAINT_PSEUDO)
    {
      constraint_type = CHOICE_CONSTRAINT_ANY;
    }
    SetValue (src_dlg->constraint_type, constraint_type);
    if (constraint_type == CHOICE_CONSTRAINT_QUAL_PRESENT)
    {
      PointerToDialog (src_dlg->present_qual_choice, src_data->qual_choice);
    }
    else if (constraint_type == CHOICE_CONSTRAINT_STRING)
    {
      PointerToDialog (src_dlg->string_qual_choice, src_data->qual_choice);
      PointerToDialog (src_dlg->string_constraint, src_data->string_constraint);
    }
    else if (constraint_type == CHOICE_CONSTRAINT_MATCH)
    {
      PointerToDialog (src_dlg->match_1, src_data->qual_choice);
      PointerToDialog (src_dlg->match_2, src_data->qual_choice_match);
    }
    else if (constraint_type == CHOICE_CONSTRAINT_PSEUDO)
    {
      PointerToDialog (src_dlg->pseudo_constraint, src_data->pseudo_constraint);
    }
    ChangeChoiceConstraintType (src_dlg->constraint_type);
  }
}

static Pointer DialogToChoiceConstraint (DialoG d)
{
  ChoiceConstraintDialogPtr src_dlg;
  ChoiceConstraintPtr       src_data;

  src_dlg = (ChoiceConstraintDialogPtr) GetObjectExtra (d);
  if (src_dlg == NULL) 
  {
    return NULL;
  }
  
  src_data = (ChoiceConstraintPtr) MemNew (sizeof (ChoiceConstraintData));
  if (src_data != NULL)
  {
    src_data->constraint_type = GetValue (src_dlg->constraint_type);
    if (src_data->constraint_type == CHOICE_CONSTRAINT_QUAL_PRESENT)
    {
      src_data->qual_choice = DialogToPointer (src_dlg->present_qual_choice);
    }
    else if (src_data->constraint_type == CHOICE_CONSTRAINT_STRING)
    {
      src_data->qual_choice = DialogToPointer (src_dlg->string_qual_choice);
      src_data->string_constraint = DialogToPointer (src_dlg->string_constraint);
    }
    else if (src_data->constraint_type == CHOICE_CONSTRAINT_MATCH)
    {
      src_data->qual_choice = DialogToPointer (src_dlg->match_1);
      src_data->qual_choice_match = DialogToPointer (src_dlg->match_2);
    }
    else if (src_data->constraint_type == CHOICE_CONSTRAINT_PSEUDO)
    {
      src_data->pseudo_constraint = DialogToPointer (src_dlg->pseudo_constraint);
    }
    src_data->free_vn_proc = src_dlg->free_vn_proc;
    src_data->copy_vn_proc = src_dlg->copy_vn_proc;
  }
  return src_data;
}

static void ChoiceConstraintMessage (DialoG d, Int2 mssg)

{
  ChoiceConstraintDialogPtr scdp;

  scdp = (ChoiceConstraintDialogPtr) GetObjectExtra (d);
  if (scdp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        ResetChoiceConstraintDialog (scdp);        
        break;
      case VIB_MSG_ENTER :
        Select (scdp->constraint_type);
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestChoiceConstraintDialog (DialoG d)

{
  return NULL;
}

extern DialoG 
ConstraintChoiceDialog 
(GrouP                     h, 
 FeatureFieldSelectionProc func_present,
 FeatureFieldSelectionProc func_string,
 FreeValNodeProc           free_vn_proc,
 CopyValNodeDataProc       copy_vn_proc,
 CharPtr                   any_name,
 CharPtr                   present_name,
 Boolean                   clear_btn,
 Boolean                   use_pseudo)
{
  ChoiceConstraintDialogPtr scdp;
  GrouP                     p, q;
  ButtoN                    b = NULL;
  
  scdp = (ChoiceConstraintDialogPtr) MemNew (sizeof (ChoiceConstraintDialogData));
  if (scdp == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, scdp, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  scdp->dialog = (DialoG) p;
  scdp->todialog = ChoiceConstraintToDialog;
  scdp->fromdialog = DialogToChoiceConstraint;
  scdp->dialogmessage = ChoiceConstraintMessage;
  scdp->testdialog = TestChoiceConstraintDialog;

  scdp->constraint_type = HiddenGroup (p, 2, 0, ChangeChoiceConstraintType);
  RadioButton (scdp->constraint_type, any_name);
  StaticPrompt (scdp->constraint_type, "", 0, dialogTextHeight, systemFont, 'l');  
  
  RadioButton (scdp->constraint_type, present_name);
  scdp->present_qual_choice = func_present (scdp->constraint_type, FALSE, NULL, NULL);
    
  RadioButton (scdp->constraint_type, "When");
  q = HiddenGroup (scdp->constraint_type, 2, 0, NULL);
  scdp->string_qual_choice = func_string (q, TRUE, NULL, NULL);
  scdp->string_constraint = StringConstraintDialogX (q, NULL, FALSE);
  
  RadioButton (scdp->constraint_type, "When");
  q = HiddenGroup (scdp->constraint_type, 3, 0, NULL);
  scdp->match_1 = func_string (q, TRUE, NULL, NULL);
  StaticPrompt (q, "Equals", 0, dialogTextHeight, systemFont, 'l');
  scdp->match_2 = func_string (q, TRUE, NULL, NULL);

  if (use_pseudo)
  {
    RadioButton (scdp->constraint_type, "When");
    scdp->pseudo_constraint = PseudoConstraintDialog (scdp->constraint_type);
  }
  else
  {
    scdp->pseudo_constraint = NULL;
  }
  
  SetValue (scdp->constraint_type, CHOICE_CONSTRAINT_ANY);
  SetObjectExtra (scdp->constraint_type, scdp, NULL);

  Disable (scdp->present_qual_choice);
  Disable (scdp->string_qual_choice);
  Disable (scdp->string_constraint);
  Disable (scdp->match_1);
  Disable (scdp->match_2);
  
  scdp->free_vn_proc = free_vn_proc;
  scdp->copy_vn_proc = copy_vn_proc;
  
  if (clear_btn)
  {
    b = PushButton (p, "Clear Constraint", ClearDialogBtn);
    SetObjectExtra (b, p, NULL);    
  }
  
  AlignObjects (ALIGN_CENTER, (HANDLE) scdp->constraint_type, (HANDLE) b, NULL);
   
  return (DialoG) p;
}

extern DialoG SourceConstraintDialogX (GrouP h, Boolean clear_btn)
{
  return ConstraintChoiceDialog (h, SourceQualTypeConstraintSelectionDialog, 
                                 SourceStringConstraintSelectionDialog,
                                 ValNodeSimpleDataFree,
                                 SourceQualValNodeDataCopy,
                                 "For any source", "When qualifier present:",
                                 clear_btn,
                                 FALSE);
}


extern DialoG CDSGeneProtConstraintDialog (GrouP h, Boolean clear_btn)
{
  return ConstraintChoiceDialog (h, CDSGeneProtFieldSelectionDialog,
                                 CDSGeneProtFieldConstraintSelectionDialog,
                                 NULL, IntValNodeCopy,
                                 "For any CDS-Gene-Prot-mRNA set", 
                                 "When field present:", 
                                 clear_btn, TRUE);
}

typedef struct LocationConstraintXdialog 
{
  DIALOG_MESSAGE_BLOCK
  PopuP  interval_end_choice;
  PopuP  interval_match_choice;
  PopuP  endpoint_match_choice;
  TexT   only_val;
  TexT   first_val;
  TexT   second_val;
  PrompT second_val_prompt;
  PopuP  strand;
  PopuP  sequence_type;
  Boolean show_interval_controls;
} LocationConstraintXDialogData, PNTR LocationConstraintXDialogPtr;

static void ShowLocationChoiceControls (PopuP p)

{
  LocationConstraintXDialogPtr lcdp;
  Int4                        match_choice;

  lcdp = (LocationConstraintXDialogPtr) GetObjectExtra (p);
  if (lcdp == NULL || lcdp->interval_end_choice == NULL) return;
  
  if (GetValue (lcdp->interval_end_choice) == LOCATION_CONSTRAINT_WHOLE_INTERVAL)
  {
    Show (lcdp->interval_match_choice);
    Hide (lcdp->endpoint_match_choice);
    if (p == lcdp->interval_end_choice)
    {
      SetValue (lcdp->interval_match_choice, GetValue (lcdp->endpoint_match_choice)); 
    }
    match_choice = GetValue (lcdp->interval_match_choice);
  }
  else
  {
    Hide (lcdp->interval_match_choice);
    Show (lcdp->endpoint_match_choice);
    if (p == lcdp->interval_end_choice)
    {
      SetValue (lcdp->endpoint_match_choice, GetValue (lcdp->interval_match_choice)); 
    }
    match_choice = GetValue (lcdp->endpoint_match_choice);
  }
  
  switch (match_choice)
  {
    case LOCATION_CONSTRAINT_ANY :
      Hide (lcdp->only_val);
      Hide (lcdp->first_val);
      Hide (lcdp->second_val_prompt);
      Hide (lcdp->second_val);
      break;
    case LOCATION_CONSTRAINT_UPSTREAM :
    case LOCATION_CONSTRAINT_DOWNSTREAM :
      Show (lcdp->only_val);
      Hide (lcdp->first_val);
      Hide (lcdp->second_val_prompt);
      Hide (lcdp->second_val);
      break;
    case LOCATION_CONSTRAINT_CONTAINED :
    case LOCATION_CONSTRAINT_NOT_IN :
    case LOCATION_CONSTRAINT_OVERLAP :
    case LOCATION_CONSTRAINT_EQUAL :
      Hide (lcdp->only_val);
      Show (lcdp->first_val);
      Show (lcdp->second_val_prompt);
      Show (lcdp->second_val);
      break;
  }
}

static void ResetLocationConstraintXDialog (LocationConstraintXDialogPtr lcdp)
{
  if (lcdp == NULL) return;
  
  SafeSetValue (lcdp->interval_end_choice, LOCATION_CONSTRAINT_WHOLE_INTERVAL);
  SafeSetValue (lcdp->interval_match_choice, LOCATION_CONSTRAINT_ANY);
  SafeSetTitle (lcdp->only_val, "");
  SafeSetTitle (lcdp->first_val, "");
  SafeSetTitle (lcdp->second_val, "");
  SetValue (lcdp->strand, LOCATION_CONSTRAINT_ANY_STRAND);
  SetValue (lcdp->sequence_type, LOCATION_CONSTRAINT_ANY_SEQ);
  ShowLocationChoiceControls (lcdp->interval_end_choice);
}

static void LocationConstraintXToDialog (DialoG d, Pointer data)

{
  LocationConstraintXDialogPtr lcdp;
  LocationConstraintXPtr       lcp;
  Int4                        interval_end_choice;
  Int4                        match_choice;
  Char                        tmp [15];
  Int4                        strand;
  Int4                        sequence_type;

  lcdp = (LocationConstraintXDialogPtr) GetObjectExtra (d);
  if (lcdp == NULL)
  {
    return;
  }
  lcp = (LocationConstraintXPtr) data;
  
  ResetLocationConstraintXDialog (lcdp);
  if (lcp != NULL)
  {
    interval_end_choice = lcp->interval_end_choice;
    if (interval_end_choice < LOCATION_CONSTRAINT_WHOLE_INTERVAL
        || interval_end_choice > LOCATION_CONSTRAINT_STOP_ENDPOINT)
    {
      interval_end_choice = LOCATION_CONSTRAINT_WHOLE_INTERVAL;
    }
    match_choice = lcp->match_choice;
    if (match_choice < LOCATION_CONSTRAINT_ANY)
    {
      match_choice = LOCATION_CONSTRAINT_ANY;
    }
    else if (match_choice > LOCATION_CONSTRAINT_EQUAL)
    {
      match_choice = LOCATION_CONSTRAINT_EQUAL;
    }
    if (match_choice > LOCATION_CONSTRAINT_NOT_IN && interval_end_choice != LOCATION_CONSTRAINT_WHOLE_INTERVAL)
    {
      interval_end_choice = LOCATION_CONSTRAINT_WHOLE_INTERVAL;
    }
    
    SetValue (lcdp->interval_end_choice, interval_end_choice);
    if (interval_end_choice == LOCATION_CONSTRAINT_WHOLE_INTERVAL)
    {
      SetValue (lcdp->interval_match_choice, match_choice);
    }
    else
    {
      SetValue (lcdp->endpoint_match_choice, match_choice);
    }
    switch (match_choice)
    {
      case LOCATION_CONSTRAINT_ANY :
        break;
      case LOCATION_CONSTRAINT_UPSTREAM :
      case LOCATION_CONSTRAINT_DOWNSTREAM :
        if (lcp->left >= 0)
        {
          sprintf (tmp, "%d", lcp->left + 1);
          SetTitle (lcdp->only_val, tmp);
        }
        break;
      case LOCATION_CONSTRAINT_CONTAINED :
      case LOCATION_CONSTRAINT_NOT_IN :
      case LOCATION_CONSTRAINT_OVERLAP :
      case LOCATION_CONSTRAINT_EQUAL :
        if (lcp->left >= 0)
        {
          sprintf (tmp, "%d", lcp->left + 1);
          SetTitle (lcdp->first_val, tmp);
        }
        if (lcp->right >= 0)
        {
          sprintf (tmp, "%d", lcp->right + 1);
          SetTitle (lcdp->second_val, tmp);
        }
        break;
    }
    strand = lcp->strand;
    if (strand < LOCATION_CONSTRAINT_ANY_STRAND || strand > LOCATION_CONSTRAINT_MINUS_STRAND)
    {
      strand = LOCATION_CONSTRAINT_ANY_STRAND;
    }
    SetValue (lcdp->strand, strand);
    sequence_type = lcp->sequence_type;
    if (sequence_type < LOCATION_CONSTRAINT_ANY_SEQ || sequence_type > LOCATION_CONSTRAINT_PROT_SEQ)
    {
      sequence_type = LOCATION_CONSTRAINT_ANY_SEQ;
    }
    SetValue (lcdp->sequence_type, LOCATION_CONSTRAINT_ANY_SEQ);
  }
  ShowLocationChoiceControls (lcdp->interval_match_choice);
}

static Pointer DialogToLocationConstraintX (DialoG d)

{
  LocationConstraintXDialogPtr lcdp;
  LocationConstraintXPtr       lcp;
  CharPtr                     tmp;

  lcdp = (LocationConstraintXDialogPtr) GetObjectExtra (d);
  if (lcdp == NULL)
  {
    return NULL;
  }
  lcp = (LocationConstraintXPtr) MemNew (sizeof (LocationConstraintXData));
  if (lcp != NULL)
  {
    if (lcdp->show_interval_controls)
    {
      lcp->interval_end_choice = GetValue (lcdp->interval_end_choice);
      if (lcp->interval_end_choice == LOCATION_CONSTRAINT_WHOLE_INTERVAL)
      {
        lcp->match_choice = GetValue (lcdp->interval_match_choice);
      }
      else
      {
        lcp->match_choice = GetValue (lcdp->endpoint_match_choice);
      }
      switch (lcp->match_choice)
      {
        case LOCATION_CONSTRAINT_ANY :
          break;
        case LOCATION_CONSTRAINT_UPSTREAM :
        case LOCATION_CONSTRAINT_DOWNSTREAM :
          tmp = SaveStringFromText (lcdp->only_val);
          if (StringHasNoText (tmp))
          {
            lcp->left = 0;
          }
          else
          {
            lcp->left = atoi (tmp) - 1;
          }
          tmp = MemFree (tmp);
          break;
        case LOCATION_CONSTRAINT_CONTAINED :
        case LOCATION_CONSTRAINT_NOT_IN :
        case LOCATION_CONSTRAINT_OVERLAP :
        case LOCATION_CONSTRAINT_EQUAL :
          tmp = SaveStringFromText (lcdp->first_val);
          if (StringHasNoText (tmp))
          {
            lcp->left = 0;
          }
          else
          {
            lcp->left = atoi (tmp) - 1;
          }
          tmp = MemFree (tmp);
          tmp = SaveStringFromText (lcdp->second_val);
          if (StringHasNoText (tmp))
          {
            lcp->right = 0;
          }
          else
          {
            lcp->right = atoi (tmp) - 1;
          }
          tmp = MemFree (tmp);
          break;
      }
    }
    else
    {
      lcp->match_choice = LOCATION_CONSTRAINT_ANY;
      lcp->left = 0;
      lcp->right = 0;
    }
    lcp->strand = GetValue (lcdp->strand);
    lcp->sequence_type = GetValue (lcdp->sequence_type);
  }
  return lcp;
}

static void LocationConstraintXMessage (DialoG d, Int2 mssg)

{
  LocationConstraintXDialogPtr lcdp;

  lcdp = (LocationConstraintXDialogPtr) GetObjectExtra (d);
  if (lcdp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        ResetLocationConstraintXDialog (lcdp);        
        break;
      case VIB_MSG_ENTER :
        if (lcdp->show_interval_controls)
        {
          Select (lcdp->interval_end_choice);
        }
        else
        {
          Select (lcdp->strand);
        }
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestLocationConstraintXDialog (DialoG d)

{
  return NULL;
}


extern DialoG 
LocationConstraintXDialog 
(GrouP   h,
 Boolean show_interval_controls, 
 Boolean clear_btn)

{
  LocationConstraintXDialogPtr lcdp;
  GrouP                       p, g, k, val_grp, strand_grp, g2, g3, g4;
  ButtoN                      b = NULL;
  
  lcdp = (LocationConstraintXDialogPtr) MemNew (sizeof (LocationConstraintXDialogData));
  if (lcdp == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, lcdp, StdCleanupExtraProc);

  lcdp->dialog = (DialoG) p;
  lcdp->todialog = LocationConstraintXToDialog;
  lcdp->fromdialog = DialogToLocationConstraintX;
  lcdp->dialogmessage = LocationConstraintXMessage;
  lcdp->testdialog = TestLocationConstraintXDialog;
  
  lcdp->show_interval_controls = show_interval_controls;

  if (lcdp->show_interval_controls)
  {
    g = HiddenGroup (p, 6, 0, NULL);
    SetGroupSpacing (g, 10, 10);
    lcdp->interval_end_choice = PopupList (g, TRUE, ShowLocationChoiceControls);
    PopupItem (lcdp->interval_end_choice, "Entire location");
    PopupItem (lcdp->interval_end_choice, "Start");
    PopupItem (lcdp->interval_end_choice, "Stop");
    SetValue (lcdp->interval_end_choice, 1);
    SetObjectExtra (lcdp->interval_end_choice, lcdp, NULL);
  
    StaticPrompt (g, "is", 0, dialogTextHeight, systemFont, 'l');
  
    k = HiddenGroup (g, 0, 0, NULL);
    lcdp->interval_match_choice = PopupList (k, TRUE, ShowLocationChoiceControls);
    PopupItem (lcdp->interval_match_choice, "Any location");
    PopupItem (lcdp->interval_match_choice, "Upstream from");
    PopupItem (lcdp->interval_match_choice, "Downstream from");
    PopupItem (lcdp->interval_match_choice, "Contained in");
    PopupItem (lcdp->interval_match_choice, "Not in");
    PopupItem (lcdp->interval_match_choice, "Overlaps");
    PopupItem (lcdp->interval_match_choice, "Equal to");
    SetValue (lcdp->interval_match_choice, 1);
    SetObjectExtra (lcdp->interval_match_choice, lcdp, NULL);
  
    lcdp->endpoint_match_choice = PopupList (k, TRUE, ShowLocationChoiceControls);
    PopupItem (lcdp->endpoint_match_choice, "Any location");
    PopupItem (lcdp->endpoint_match_choice, "Upstream from");
    PopupItem (lcdp->endpoint_match_choice, "Downstream from");
    PopupItem (lcdp->endpoint_match_choice, "Contained in");
    PopupItem (lcdp->endpoint_match_choice, "Not in");
    SetValue (lcdp->endpoint_match_choice, 1);
    SetObjectExtra (lcdp->endpoint_match_choice, lcdp, NULL);

    g2 = HiddenGroup (g, 2, 0, NULL);
    SetGroupSpacing (g2, 10, 10);
    val_grp = HiddenGroup (g2, 0, 0, NULL);
    lcdp->only_val = DialogText (val_grp, "", 5, NULL);
    g3 = HiddenGroup (val_grp, 3, 0, NULL);
    lcdp->first_val = DialogText (g3, "", 5, NULL);
    lcdp->second_val_prompt = StaticPrompt (g3, "to", 0, dialogTextHeight, systemFont, 'l');
    lcdp->second_val = DialogText (g3, "", 5, NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) lcdp->only_val, (HANDLE) g3, NULL);
  }
  else
  {
    g = HiddenGroup (p, 2, 0, NULL);
  }
  
  strand_grp = HiddenGroup (g, 2, 0, NULL);
  StaticPrompt (strand_grp, "on", 0, dialogTextHeight, systemFont, 'l');
  lcdp->strand = PopupList (strand_grp, TRUE, NULL);
  PopupItem (lcdp->strand, "Any strand");
  PopupItem (lcdp->strand, "Plus strand");
  PopupItem (lcdp->strand, "Minus strand");
  SetValue (lcdp->strand, LOCATION_CONSTRAINT_ANY_STRAND);
  
  g4 = HiddenGroup (g, 3, 0, NULL);
  SetGroupSpacing (g4, 10, 10);
  StaticPrompt (g4, "on", 0, dialogTextHeight, systemFont, 'l');
  lcdp->sequence_type = PopupList (g4, TRUE, NULL);
  PopupItem (lcdp->sequence_type, "nucleotide and protein sequences");
  PopupItem (lcdp->sequence_type, "nucleotide sequences only");
  PopupItem (lcdp->sequence_type, "protein sequences only");
  SetValue (lcdp->sequence_type, LOCATION_CONSTRAINT_ANY_SEQ);

  if (clear_btn)
  {
    b = PushButton (p, "Clear Constraint", ClearDialogBtn);
    SetObjectExtra (b, p, NULL);
  }

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) b, NULL);
  
  ShowLocationChoiceControls (lcdp->interval_match_choice);
  return (DialoG) p;
}

static Boolean DoesStrandMatchConstraint (SeqLocPtr slp, LocationConstraintXPtr lcp)
{
  Uint2 strand;
  Boolean rval = FALSE;
  
  if (slp == NULL)
  {
    rval = FALSE;
  }
  else if (lcp == NULL || lcp->strand == LOCATION_CONSTRAINT_ANY_STRAND)
  {
    rval = TRUE;
  }
  else
  {
    strand = SeqLocStrand (slp);
    if (strand == Seq_strand_minus)
    {
      if (lcp->strand == LOCATION_CONSTRAINT_MINUS_STRAND)
      {
        rval = TRUE;
      }
      else
      {
        rval = FALSE;
      }
    }
    else
    {
      if (lcp->strand == LOCATION_CONSTRAINT_PLUS_STRAND)
      {
        rval = TRUE;
      }
      else
      {
        rval = FALSE;
      }
    }
  }
  return rval;
}

static Boolean DoesSequenceTypeMatchContraint (SeqLocPtr slp, LocationConstraintXPtr lcp)
{
  Boolean   rval = FALSE;
  BioseqPtr bsp;
  
  if (slp == NULL)
  {
    rval = FALSE;
  }
  else if (lcp == NULL || lcp->sequence_type == LOCATION_CONSTRAINT_ANY_SEQ)
  {
    rval = TRUE;
  }
  else
  {
    bsp = BioseqFindFromSeqLoc (slp);
    if (bsp != NULL)
    {
      if (ISA_na (bsp->mol) && lcp->sequence_type == LOCATION_CONSTRAINT_NUC_SEQ)
      {
        rval = TRUE;
      }
      else if (ISA_aa (bsp->mol) && lcp->sequence_type == LOCATION_CONSTRAINT_PROT_SEQ)
      {
        rval = TRUE;
      }
    }
  }
  return rval;
}

extern Boolean DoesLocationMatchConstraint (SeqLocPtr slp, LocationConstraintXPtr lcp)

{
  Boolean rval = FALSE;
  Int4    loc_start, loc_stop, endpoint = 0;
  
  if (slp == NULL)
  {
    return FALSE;
  }
  
  if (lcp == NULL)
  {
    rval = TRUE;
  }
  else if (! DoesStrandMatchConstraint (slp, lcp))
  {
    rval = FALSE;
  }
  else if (! DoesSequenceTypeMatchContraint (slp, lcp))
  {
    rval = FALSE;
  }
  else
  {
    switch (lcp->match_choice)
    {
      case LOCATION_CONSTRAINT_ANY :
        rval = TRUE;
        break;
      case LOCATION_CONSTRAINT_UPSTREAM :
        if (lcp->interval_end_choice == LOCATION_CONSTRAINT_WHOLE_INTERVAL
            || lcp->interval_end_choice == LOCATION_CONSTRAINT_STOP_ENDPOINT)
        {
          endpoint = SeqLocStop (slp);
        }
        else
        {
          endpoint = SeqLocStart (slp);
        }
        if (endpoint < lcp->left)
        {
          rval = TRUE;
        }
        break;
      case LOCATION_CONSTRAINT_DOWNSTREAM :
        if (lcp->interval_end_choice == LOCATION_CONSTRAINT_WHOLE_INTERVAL
            || lcp->interval_end_choice == LOCATION_CONSTRAINT_START_ENDPOINT)
        {
          endpoint = SeqLocStart (slp);
        }
        else
        {
          endpoint = SeqLocStop (slp);
        }
        if (endpoint > lcp->left)
        {
          rval = TRUE;
        }
        break;
      case LOCATION_CONSTRAINT_CONTAINED :
        if (lcp->interval_end_choice == LOCATION_CONSTRAINT_WHOLE_INTERVAL)
        {
          loc_start = SeqLocStart (slp);
          loc_stop = SeqLocStop (slp);
        }
        else
        {
          if (lcp->interval_end_choice == LOCATION_CONSTRAINT_START_ENDPOINT)
          {
            endpoint = SeqLocStart (slp);
          }
          else
          {
            endpoint = SeqLocStop (slp);
          }
          loc_start = endpoint;
          loc_stop = endpoint;
        }
        if (loc_start >= lcp->left && loc_stop <= lcp->right)
        {
          rval = TRUE;
        }
        break;
      case LOCATION_CONSTRAINT_NOT_IN :
        if (lcp->interval_end_choice == LOCATION_CONSTRAINT_WHOLE_INTERVAL)
        {
          loc_start = SeqLocStart (slp);
          loc_stop = SeqLocStop (slp);
        }
        else
        {
          if (lcp->interval_end_choice == LOCATION_CONSTRAINT_START_ENDPOINT)
          {
            endpoint = SeqLocStart (slp);
          }
          else
          {
            endpoint = SeqLocStop (slp);
          }
          loc_start = endpoint;
          loc_stop = endpoint;
        }
        if (loc_stop <= lcp->left || loc_start >= lcp->right)
        {
          rval = TRUE;
        }
        break;
      case LOCATION_CONSTRAINT_OVERLAP :
        loc_start = SeqLocStart (slp);
        loc_stop = SeqLocStop (slp);
        if (loc_start <= lcp->left && loc_stop >= lcp->left)
        {
          rval = TRUE;
        }
        else if (loc_start <= lcp->right && loc_stop >= lcp->right)
        {
          rval = TRUE;
        }
        else if (loc_start >= lcp->left && loc_stop <= lcp->right)
        {
          rval = TRUE;
        }
        break;
      case LOCATION_CONSTRAINT_EQUAL :
        loc_start = SeqLocStart (slp);
        loc_stop = SeqLocStop (slp);
        if (loc_start == lcp->left && loc_stop == lcp->right)
        {
          rval = TRUE;
        }
        break;
    }
  }
  return rval; 
}

typedef struct filterform
{
  DIALOG_MESSAGE_BLOCK
  
  DialoG string_constraint;
  DialoG source_constraint;
  DialoG location_constraint;
  DialoG cds_gene_prot_constraint;
  DialoG id_list_constraint;
  
  DialoG tbs;
  DialoG pages[6];
  Int4   current_page;
  
} FilterFormData, PNTR FilterFormPtr;

static void FilterToDialog (DialoG d, Pointer userdata)
{
  FilterFormPtr dlg;
  FilterSetPtr dlg_data;
  
  dlg = (FilterFormPtr) GetObjectExtra (d);
  
  if (dlg == NULL)
  {
    return;
  }
  dlg_data = (FilterSetPtr) userdata;
  PointerToDialog (dlg->string_constraint, NULL);
  PointerToDialog (dlg->source_constraint, NULL);
  PointerToDialog (dlg->location_constraint, NULL);
  PointerToDialog (dlg->cds_gene_prot_constraint, NULL);
  PointerToDialog (dlg->id_list_constraint, NULL);
  
  if (dlg_data != NULL)
  {
    PointerToDialog (dlg->string_constraint, dlg_data->scp);
    PointerToDialog (dlg->source_constraint, dlg_data->ccp);
    PointerToDialog (dlg->location_constraint, dlg_data->lcp);
    PointerToDialog (dlg->cds_gene_prot_constraint, dlg_data->cgp);
    PointerToDialog (dlg->id_list_constraint, dlg_data->id_list);
  }
}

static void ApplyNthPage (FilterFormPtr dlg, FilterSetPtr dlg_data, Int4 page)
{
  if (page < 0) return;

  /* adjust page value to reflect missing pages */
  if (dlg->source_constraint == NULL)
  {
    page++;
  }
  if (dlg->string_constraint == NULL && page >= 1)
  {
    page++;
  }
  if (dlg->location_constraint == NULL && page >= 2)
  {
    page++;
  }
  if (dlg->cds_gene_prot_constraint == NULL && page >= 3)
  {
    page++;
  }
  if (dlg->id_list_constraint == NULL && page >= 4)
  {
    page++;
  }
  
  switch (page)
  {
    case 0:
      dlg_data->ccp = DialogToPointer (dlg->source_constraint);
      break;
    case 1:
      dlg_data->scp = DialogToPointer (dlg->string_constraint);
      break;
    case 2:
      dlg_data->lcp = DialogToPointer (dlg->location_constraint);
      break;
    case 3:
      dlg_data->cgp = DialogToPointer (dlg->cds_gene_prot_constraint);
      break;
    case 4:
      dlg_data->id_list = DialogToPointer (dlg->id_list_constraint);
      break;
  }
  
}


static Pointer DialogToFilter (DialoG d)
{
  FilterFormPtr dlg;
  FilterSetPtr dlg_data;
  
  dlg = (FilterFormPtr) GetObjectExtra (d);
  
  if (dlg == NULL)
  {
    return NULL;
  }
  dlg_data = (FilterSetPtr) MemNew (sizeof (FilterSetData));
  
  if (dlg_data != NULL)
  {
    if (dlg->tbs == NULL) 
    {
      dlg_data->scp = DialogToPointer (dlg->string_constraint);
      dlg_data->ccp = DialogToPointer (dlg->source_constraint);
      dlg_data->lcp = DialogToPointer (dlg->location_constraint);
      dlg_data->cgp = DialogToPointer (dlg->cds_gene_prot_constraint);
      dlg_data->id_list = DialogToPointer (dlg->id_list_constraint);
    }
    else
    {
      /* only take one value */
      ApplyNthPage (dlg, dlg_data, dlg->current_page);
    }
  }
  return dlg_data;
}

static void FilterFormMessage (DialoG d, Int2 mssg)

{
}

static ValNodePtr TestFilterFormDialog (DialoG d)

{
  return NULL;
}

static void ChangeFilterPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  FilterFormPtr ffp;

  ffp = (FilterFormPtr) data;
  if (ffp != NULL) {
    ffp->current_page = newval;
    SafeHide (ffp->pages [oldval]);
    Update ();
    SafeShow (ffp->pages [newval]);
    Update ();
  }
}

static void ClearAllConstraints (ButtoN b)
{
  FilterFormPtr ffp;

  ffp = (FilterFormPtr) GetObjectExtra (b);
  if (ffp != NULL)
  {
    PointerToDialog (ffp->dialog, NULL);
  }
}

extern DialoG 
FilterGroup 
(GrouP h,
 Boolean has_string_constraint,
 Boolean has_source_constraint,
 Boolean has_location_constraint,
 Boolean has_cds_gene_prot_constraint,
 Boolean has_id_list_constraint,
 CharPtr string_constraint_label)
{
  GrouP         g, k;
  FilterFormPtr ffp;
  Int4          num_pages = 0;
  CharPtr       filterTabs[5];
  Int4          i;
  ButtoN        clear_constraints = NULL;
  
  if (! has_string_constraint
      && ! has_source_constraint
      && ! has_location_constraint
      && ! has_cds_gene_prot_constraint
      && ! has_id_list_constraint)
  {
    return NULL;
  }
  
  ffp = (FilterFormPtr) MemNew (sizeof (FilterFormData));
  if (ffp == NULL)
  {
    return NULL;
  }

  g = NormalGroup (h, -1, 0, NULL, programFont, NULL);
  SetObjectExtra (g, ffp, StdCleanupExtraProc);

  ffp->dialog = (DialoG) g;
  ffp->todialog = FilterToDialog;
  ffp->fromdialog = DialogToFilter;
  ffp->dialogmessage = FilterFormMessage;
  ffp->testdialog = TestFilterFormDialog;
  
  ffp->current_page = 0;

  if (has_source_constraint)
  {
    filterTabs [num_pages++] = "Source Constraint";
  }
  if (has_string_constraint)
  {
    filterTabs [num_pages++] = "String Constraint";
  }
  if (has_location_constraint)
  {
    filterTabs [num_pages++] = "Location Constraint";
  }
  if (has_cds_gene_prot_constraint)
  {
    filterTabs [num_pages++] = "CDS-Gene-Prot-mRNA Set Constraint";
  }
  if (has_id_list_constraint)
  {
    filterTabs [num_pages++] = "Sequence ID Constraint";
  }
  filterTabs [num_pages] = NULL;
  
  if (num_pages > 1)
  {
  
    ffp->tbs = CreateFolderTabs (g, filterTabs, 0,
                                 0, 0, PROGRAM_FOLDER_TAB,
                                 ChangeFilterPage, (Pointer) ffp);
    k = HiddenGroup (g, 0, 0, NULL);
    num_pages = 0;
    if (has_source_constraint)
    {
      ffp->source_constraint = SourceConstraintDialogX (k, FALSE);
      ffp->pages [num_pages++] = ffp->source_constraint;
    }
    if (has_string_constraint)
    {
      ffp->string_constraint = StringConstraintDialogX (k, string_constraint_label, FALSE);
      ffp->pages [num_pages++] = ffp->string_constraint;
    }
    if (has_location_constraint)
    {
      ffp->location_constraint = LocationConstraintXDialog (k, FALSE, FALSE);
      ffp->pages [num_pages++] = ffp->location_constraint;
    }
    if (has_cds_gene_prot_constraint)
    {
      ffp->cds_gene_prot_constraint = CDSGeneProtConstraintDialog (k, FALSE);
      ffp->pages [num_pages++] = ffp->cds_gene_prot_constraint;
    }
    if (has_id_list_constraint)
    {
      ffp->id_list_constraint = StringConstraintDialogX (k, "Where sequence ID", FALSE);
      ffp->pages[num_pages++] = ffp->id_list_constraint;
    }
    for (i = 1; i < num_pages; i++)
    {
      Hide (ffp->pages [i]);
    }
    AlignObjects (ALIGN_CENTER, (HANDLE) ffp->pages [0],
                                (HANDLE) ffp->pages [1],
                                (HANDLE) ffp->pages [2],
                                (HANDLE) ffp->pages [3],
                                NULL);
    clear_constraints = PushButton (g, "Clear Constraints", ClearAllConstraints);
    SetObjectExtra (clear_constraints, ffp, NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) ffp->tbs, (HANDLE) k, (HANDLE) clear_constraints, NULL);
  }
  else
  {
    if (has_source_constraint)
    {
      ffp->source_constraint = SourceConstraintDialogX (g, TRUE);
    }
    else if (has_string_constraint)
    {
      ffp->string_constraint = StringConstraintDialogX (g, string_constraint_label, TRUE);
    }
    else if (has_location_constraint)
    {
      ffp->location_constraint = LocationConstraintXDialog (g, FALSE, TRUE);
    }
    else if (has_cds_gene_prot_constraint)
    {
      ffp->cds_gene_prot_constraint = CDSGeneProtConstraintDialog (g, TRUE);
    }
    else if (has_id_list_constraint)
    {
      ffp->id_list_constraint = StringConstraintDialogX (g, "Where sequence ID", TRUE);
    }
  }
  return (DialoG) g;
}

/* Operations on constrained features and descriptors */
typedef struct objecthasstring
{
  StringConstraintXPtr scp;
  Boolean             found;
} ObjectHasStringData, PNTR ObjectHasStringPtr;

typedef struct constraintop 
{
  FilterSetPtr           fsp;
  ObjMgrPtr              omp;
  ObjMgrTypePtr          omtp;
  AsnIoPtr               aip;
  Uint2                  entityID;
  ObjectHasStringData    ohsd;
  Pointer                userdata;
  FeatureActionProc      feature_action;
  Uint1                  seqFeatChoice;
  Uint1                  featDefChoice;
  DescriptorActionProc   descriptor_action;
  Uint1                  descriptorChoice;
} ConstraintOpData, PNTR ConstraintOpPtr;

extern Boolean DoesStringMatchConstraintX (CharPtr pchSource, StringConstraintXPtr scp)
{
  CharPtr pFound;
  Boolean rval = FALSE;
  Char    char_after = 0;
  
  if (pchSource == NULL) return FALSE;
  
  if (scp == NULL || StringHasNoText (scp->match_text)) return TRUE;

  switch (scp->match_location) 
  {
    case eStringConstraintContains:
        if (scp->insensitive)
        {
          pFound = StringISearch (pchSource, scp->match_text);
        }
        else
        {
          pFound = StringSearch (pchSource, scp->match_text);
        }
      if (pFound == NULL) 
      {
        rval = FALSE;
      }
      else if (scp->whole_word) 
      {
        rval = IsWholeWordMatch (pchSource, pFound, StringLen (scp->match_text));
        while (!rval && pFound != NULL) 
        {
            if (scp->insensitive)
            {
              pFound = StringISearch (pFound + 1, scp->match_text);
            }
            else
            {
              pFound = StringSearch (pFound + 1, scp->match_text);
            }
          if (pFound != NULL)
          {
            rval = IsWholeWordMatch (pchSource, pFound, StringLen (scp->match_text));
          }
        }
      }
      else
      {
        rval = TRUE;
      }
      break;
    case eStringConstraintStarts:
        if (scp->insensitive)
        {
          pFound = StringISearch (pchSource, scp->match_text);
        }
        else
        {
          pFound = StringSearch (pchSource, scp->match_text);
        }
      if (pFound == pchSource)
      {
        if (scp->whole_word) 
        {
          rval = IsWholeWordMatch (pchSource, pFound, StringLen (scp->match_text));
        }
        else
        {
          rval = TRUE;
        }
      }
      break;
    case eStringConstraintEnds:
        if (scp->insensitive)
        {
          pFound = StringISearch (pchSource, scp->match_text);
        }
        else
        {
          pFound = StringSearch (pchSource, scp->match_text);
        }
      while (pFound != NULL && !rval) {
          char_after = *(pFound + StringLen (scp->match_text));
        if (char_after == 0)
        {
          if (scp->whole_word) 
          {
            rval = IsWholeWordMatch (pchSource, pFound, StringLen (scp->match_text));
          }
          else
          {
            rval = TRUE;
          }
          /* stop the search, we're at the end of the string */
          pFound = NULL;
        }
        else
        {
            if (scp->insensitive)
            {
              pFound = StringISearch (pFound + 1, scp->match_text);
            }
            else
            {
              pFound = StringSearch (pFound + 1, scp->match_text);
            }
        }
      }
      break;
    case eStringConstraintEquals:
      if (scp->insensitive) 
      {
        if (StringICmp (pchSource, scp->match_text) == 0) 
        {
          rval = TRUE;
        }
      }
      else
      {
        if (StringCmp (pchSource, scp->match_text) == 0) 
        {
          rval = TRUE;
        }
      }
      break;
    case eStringConstraintInList:
        if (scp->insensitive)
        {
          pFound = StringISearch (scp->match_text, pchSource);
        }
        else
        {
          pFound = StringSearch (scp->match_text, pchSource);
        }
      if (pFound == NULL) 
      {
        rval = FALSE;
      }
      else
      {
        rval = IsWholeWordMatch (scp->match_text, pFound, StringLen (pchSource));
        while (!rval && pFound != NULL) 
        {
            if (scp->insensitive)
            {
              pFound = StringISearch (pFound + 1, pchSource);
            }
            else
            {
              pFound = StringSearch (pFound + 1, pchSource);
            }
          if (pFound != NULL)
          {
            rval = IsWholeWordMatch (scp->match_text, pFound, StringLen (pchSource));
          }
        }
      }
      if (!rval) {
        /* look for spans */
        rval = IsStringInSpanInList (pchSource, scp->match_text);
      }
      break;
    }
    return rval;
}

static void LIBCALLBACK AsnWriteConstraintCallBack (AsnExpOptStructPtr pAEOS)

{
  CharPtr            pchSource;
  ObjectHasStringPtr ohsp;

  ohsp = (ObjectHasStringPtr) pAEOS->data;
  if (ISA_STRINGTYPE (AsnFindBaseIsa (pAEOS->atp))) 
  {
      pchSource = (CharPtr) pAEOS->dvp->ptrvalue;
      ohsp->found |= DoesStringMatchConstraintX (pchSource, ohsp->scp);
  }
}

static Boolean 
DoesObjectMatchStringConstraint 
(ObjMgrTypePtr      omtp,
 AsnIoPtr           aip, 
 Pointer            ptr, 
 ObjectHasStringPtr ohsp)

{
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  CharPtr           search_txt = NULL;
  
  if (omtp == NULL || ohsp == NULL)
  {
    return FALSE;
  }
  ohsp->found = FALSE;
  (omtp->asnwrite) (ptr, aip, NULL);
  
  if (!ohsp->found && omtp->datatype == OBJ_SEQFEAT)
  {
    sfp = (SeqFeatPtr) ptr;
    if (SeqMgrFeaturesAreIndexed(sfp->idx.entityID) == 0) {
      SeqMgrIndexFeatures (sfp->idx.entityID, NULL);
    }

    sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &fcontext);
    ohsp->found = DoesStringMatchConstraintX (fcontext.label, ohsp->scp);
    if (!ohsp->found && sfp != NULL && sfp->idx.subtype == FEATDEF_tRNA)
    {
      search_txt = (CharPtr) MemNew ((StringLen (fcontext.label) + 6) * sizeof (Char));
      if (search_txt != NULL)
      {
        sprintf (search_txt, "tRNA-%s", fcontext.label);
        ohsp->found = DoesStringMatchConstraintX (search_txt, ohsp->scp);
        search_txt = MemFree (search_txt);
      }
    }
  }
  return ohsp->found;
}


static Boolean DoesSourceHaveOneQualPresent (BioSourcePtr biop, SourceQualDescPtr sqdp)
{
  OrgModPtr    mod;
  SubSourcePtr ssp;
  Boolean      is_present = FALSE;
  
  if (biop == NULL || sqdp == NULL)
  {
    return FALSE;
  }
  
  if (sqdp->isOrgMod)
  {
    if (biop->org != NULL && biop->org->orgname != NULL)
    {
      mod = biop->org->orgname->mod;
      while (mod != NULL && mod->subtype != sqdp->subtype)
      {
        mod = mod->next;
      }
      if (mod != NULL && mod->subtype == sqdp->subtype)
      {
        is_present = TRUE;
      }
    }
  }
  else
  {
    ssp = biop->subtype;
    while (ssp != NULL && ssp->subtype != sqdp->subtype)
    {
      ssp = ssp->next;
    }
    if (ssp != NULL && ssp->subtype == sqdp->subtype)
    {
      is_present = TRUE;
    }
  }
  return is_present;
}

static Boolean DoesSourceHaveQualPresent (BioSourcePtr biop, ValNodePtr qual_list)
{
  ValNodePtr        vnp;
  Boolean           qual_found = FALSE;
  
  if (biop == NULL)
  {
    return FALSE;
  }
  if (qual_list == NULL)
  {
    return TRUE;
  }
  
  for (vnp = qual_list; vnp != NULL && !qual_found; vnp = vnp->next)
  {
    qual_found = DoesSourceHaveOneQualPresent (biop, (SourceQualDescPtr) vnp->data.ptrvalue);
  }
  
  return qual_found;
}

static Boolean IsMatchAny (ValNodePtr vnp)
{
  if (vnp == NULL
      || vnp->data.ptrvalue == NULL
      || (vnp->choice > 0 && StringICmp (vnp->data.ptrvalue, "Organism or Any Qual") == 0))
  {
    return TRUE;
  } else {
    return FALSE;
  }
}

static Boolean DoFieldsMatch (BioSourcePtr biop, ChoiceConstraintPtr ccp)
{
  ValNodePtr swap, swap2;
  ValNodePtr qual_choice_list = NULL, vnp1, vnp2;
  Boolean    rval = FALSE;
  ChoiceConstraintData tmp_const;
  StringConstraintData scd;
  SourceQualDescPtr sqdp;
  CharPtr    location;
  OrgModPtr              mod = NULL;
  SubSourcePtr           ssp;

  if (biop == NULL) return FALSE;
  if (ccp == NULL) return TRUE;

  /* Init temporary choice constraint */
  MemSet (&tmp_const, 0, sizeof (ChoiceConstraintData));
  MemSet (&scd, 0, sizeof (StringConstraintData));
  scd.match_location = eStringConstraintEquals;
  tmp_const.string_constraint = &scd;
  tmp_const.constraint_type = CHOICE_CONSTRAINT_STRING;

  if (IsMatchAny (ccp->qual_choice) && IsMatchAny (ccp->qual_choice_match)) {
    /* both are wildcards */
    ValNodeAddPointer (&qual_choice_list, 1, StringSave ("Organism"));
    qual_choice_list->next = GetSourceQualDescList (TRUE, TRUE, TRUE, FALSE);
    ValNodeAddPointer (&qual_choice_list, 1, StringSave ("Lineage"));
    ValNodeAddPointer (&qual_choice_list, 1, StringSave ("Location"));

    /* save original wildcards */
    swap = ccp->qual_choice;
    swap2 = ccp->qual_choice_match;
    for (vnp1 = qual_choice_list; vnp1 != NULL && vnp1->next != NULL && !rval; vnp1 = vnp1->next) {
      ccp->qual_choice = vnp1;
      for (vnp2 = vnp1->next; vnp2 != NULL && !rval; vnp2 = vnp2->next) {
        ccp->qual_choice_match = vnp2;
        if (DoFieldsMatch (biop, ccp)) {
          rval = TRUE;
        }
      }
    }
    /* put original wildcards back */
    ccp->qual_choice = swap;
    ccp->qual_choice_match = swap2;
    qual_choice_list = ValNodeFreeData (qual_choice_list);
  } else if (IsMatchAny (ccp->qual_choice) || IsMatchAny (ccp->qual_choice_match)) {
    /* one is a wild card */
    ValNodeAddPointer (&qual_choice_list, 1, StringSave ("Organism"));
    qual_choice_list->next = GetSourceQualDescList (TRUE, TRUE, TRUE, FALSE);
    ValNodeAddPointer (&qual_choice_list, 1, StringSave ("Lineage"));
    ValNodeAddPointer (&qual_choice_list, 1, StringSave ("Location"));

    /* if one wild card, make it the first one */    
    if (IsMatchAny (ccp->qual_choice_match)) {
      swap = ccp->qual_choice;
      ccp->qual_choice = ccp->qual_choice_match;
      ccp->qual_choice_match = swap;
    }
    /* save original wildcard */
    swap = ccp->qual_choice;
    for (vnp1 = qual_choice_list; vnp1 != NULL && !rval; vnp1 = vnp1->next) {
      if (SourceQualValNodeMatch (vnp1, ccp->qual_choice_match)) {
        /* don't compare something to itself */
        continue;
      } else {
        ccp->qual_choice = vnp1;
        rval = DoFieldsMatch (biop, ccp);
      }
    }
    /* put original wildcard back */
    ccp->qual_choice = swap;
    qual_choice_list = ValNodeFreeData (qual_choice_list);
  } 
  else if (ccp->qual_choice->choice > 0) 
  {  
    tmp_const.qual_choice = ccp->qual_choice_match;
    if (StringICmp (ccp->qual_choice->data.ptrvalue, "Organism") == 0)
    {
      if (biop->org != NULL && !StringHasNoText (biop->org->taxname))
      {
        tmp_const.string_constraint->match_text = biop->org->taxname;
        rval = DoesOneSourceMatchConstraint (biop, &tmp_const);
      }
    }
    else if (StringICmp (ccp->qual_choice->data.ptrvalue, "Lineage") == 0)
    {
      if (biop->org != NULL && biop->org->orgname != NULL && !StringHasNoText (biop->org->orgname->lineage))
      {
        tmp_const.string_constraint->match_text = biop->org->orgname->lineage;
        rval = DoesOneSourceMatchConstraint (biop, &tmp_const);
      }
    }
    else if (StringICmp (ccp->qual_choice->data.ptrvalue, "Location") == 0)
    {
      location = GetLocationFromSource (biop, NULL);
      tmp_const.string_constraint->match_text = location;
      if (!StringHasNoText (location)) {
        rval = DoesOneSourceMatchConstraint (biop, &tmp_const);
      }
      location = MemFree (location);
    }
  }
  else
  {
    tmp_const.qual_choice = ccp->qual_choice_match;
    sqdp = (SourceQualDescPtr) ccp->qual_choice->data.ptrvalue;
    if (sqdp->isOrgMod)
    {
      if (biop->org != NULL && biop->org->orgname != NULL)
      {
        mod = biop->org->orgname->mod;
      }
      while (! rval && mod != NULL)
      {
        if (mod->subtype == sqdp->subtype && !StringHasNoText (mod->subname))
        {
          tmp_const.string_constraint->match_text = mod->subname;
          rval = DoesOneSourceMatchConstraint (biop, &tmp_const);
        }
        mod = mod->next;
      }
    }
    else
    {
      ssp = biop->subtype;
      while (!rval && ssp != NULL)
      {
        if (ssp->subtype == sqdp->subtype && !StringHasNoText (ssp->name))
        {
          tmp_const.string_constraint->match_text = ssp->name;
          rval = DoesOneSourceMatchConstraint (biop, &tmp_const);
        }
        ssp = ssp->next;
      }
    }
  }
  return rval;
}

extern Boolean 
DoesOneSourceMatchConstraint 
(BioSourcePtr biop, ChoiceConstraintPtr scp)
{
  Boolean                does_match = FALSE;
  OrgModPtr              mod = NULL;
  SubSourcePtr           ssp;
  SourceQualDescPtr      sqdp;
  CharPtr                location;
  
  if (scp == NULL || scp->constraint_type == CHOICE_CONSTRAINT_ANY)
  {
    return TRUE;
  }
  if (biop == NULL)
  {
    return FALSE;
  }
  
  if (scp->constraint_type == CHOICE_CONSTRAINT_STRING)
  {
    if (IsMatchAny (scp->qual_choice))
    {
      if (biop->org != NULL)
      {
        does_match = DoesStringMatchConstraintX (biop->org->taxname, scp->string_constraint);
        if (biop->org->orgname != NULL)
        {
          mod = biop->org->orgname->mod;
        }
        while (! does_match && mod != NULL)
        {
          does_match = DoesStringMatchConstraintX (mod->subname, scp->string_constraint);
          mod = mod->next;
        }
      }
      ssp = biop->subtype;
      while (!does_match && ssp != NULL)
      {
        does_match = DoesStringMatchConstraintX (ssp->name, scp->string_constraint);
        ssp = ssp->next;
      }
    }
    else if (scp->qual_choice->choice > 0 && StringICmp (scp->qual_choice->data.ptrvalue, "Organism") == 0)
    {
      if (biop->org != NULL)
      {
        does_match = DoesStringMatchConstraintX (biop->org->taxname, scp->string_constraint);
      }
    }
    else if (scp->qual_choice->choice > 0 && StringICmp (scp->qual_choice->data.ptrvalue, "Lineage") == 0)
    {
      if (biop->org != NULL && biop->org->orgname != NULL)
      {
        does_match = DoesStringMatchConstraintX (biop->org->orgname->lineage,
                                                scp->string_constraint);
      }
    }
    else if (scp->qual_choice->choice > 0 && StringICmp (scp->qual_choice->data.ptrvalue, "Location") == 0)
    {
      location = GetLocationFromSource (biop, NULL);
      does_match = DoesStringMatchConstraintX (location, scp->string_constraint);
      location = MemFree (location);
    }
    else
    {
      sqdp = (SourceQualDescPtr) scp->qual_choice->data.ptrvalue;
      if (sqdp->isOrgMod)
      {
        if (biop->org != NULL && biop->org->orgname != NULL)
        {
          mod = biop->org->orgname->mod;
        }
        while (! does_match && mod != NULL)
        {
          if (mod->subtype == sqdp->subtype)
          {
            does_match = DoesStringMatchConstraintX (mod->subname, scp->string_constraint);
          }
          mod = mod->next;
        }
      }
      else
      {
        ssp = biop->subtype;
        while (!does_match && ssp != NULL)
        {
          if (ssp->subtype == sqdp->subtype)
          {
            does_match = DoesStringMatchConstraintX (ssp->name, scp->string_constraint);
          }
          ssp = ssp->next;
        }
      }
    }
    if (scp->string_constraint->not_present)
    {
      does_match = ! does_match;
    }
  }
  else if (scp->constraint_type == CHOICE_CONSTRAINT_MATCH)
  {
    does_match = DoFieldsMatch (biop, scp);
  }
  else
  {
    does_match = DoesSourceHaveQualPresent (biop, scp->qual_choice);
  }
  return does_match;
}

static void SeqEntryConstrainedFeaturesCallback (SeqFeatPtr sfp, Pointer userdata)
{
  SeqMgrFeatContext fcontext;
  ConstraintOpPtr   cop;
  Boolean           feature_matches = TRUE;
  
  if (sfp == NULL || userdata == NULL) return;
  cop = (ConstraintOpPtr) userdata;
  sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &fcontext);
  if (sfp == NULL) return;
  
  if (cop->seqFeatChoice != 0 && cop->seqFeatChoice != sfp->data.choice) return;
  if (cop->featDefChoice != 0 && cop->featDefChoice != sfp->idx.subtype) return;
  
  if (cop->fsp != NULL && cop->fsp->scp != NULL)
  {
    if (cop->fsp->scp->not_present)
    {
      if (DoesObjectMatchStringConstraint (cop->omtp, cop->aip, (Pointer) sfp, &(cop->ohsd))
          || DoesStringMatchConstraintX (fcontext.label, cop->fsp->scp))
      {
        feature_matches = FALSE; 
      }
    }
    else
    {
      if (! DoesObjectMatchStringConstraint (cop->omtp, cop->aip, (Pointer) sfp, &(cop->ohsd))
          && ! DoesStringMatchConstraintX (fcontext.label, cop->fsp->scp))
      {
        feature_matches = FALSE;
      }
    }
  }
  if (cop->fsp != NULL && cop->fsp->lcp != NULL)
  {
    if (! DoesLocationMatchConstraint (sfp->location, cop->fsp->lcp))
    {
      feature_matches = FALSE;
    }
  }
  
  if (feature_matches)
  {
    (cop->feature_action) (sfp, cop->userdata, cop->fsp);
  }
}


static void BioseqConstrainedFeaturesCallback (BioseqPtr bsp, Pointer userdata)
{
  ConstraintOpPtr   cop;
  SeqFeatPtr        sfp;
  ObjMgrDataPtr     omdp;
  BioseqExtraPtr    bspextra;
  Int4              index;
  
  if (bsp == NULL || userdata == NULL) return;
  cop = (ConstraintOpPtr) userdata;

  omdp = SeqMgrGetOmdpForBioseq (bsp);
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return;

  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) {
    SeqMgrIndexFeatures (bsp->idx.entityID, NULL);
    omdp = SeqMgrGetOmdpForBioseq (bsp);
    if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return;

    bspextra = (BioseqExtraPtr) omdp->extradata;
  }
  if (bspextra == NULL) return;

  for (index = 0; index < bspextra->numfeats; index++)
  {
    sfp = bspextra->featsByPos[index]->sfp;
    SeqEntryConstrainedFeaturesCallback (sfp, userdata);
  }
}


static void SeqEntryConstrainedDescriptorsCallback (SeqDescrPtr sdp, Pointer userdata)
{
  ConstraintOpPtr cop;
  Boolean         descriptor_matches = TRUE;
  
  if (sdp == NULL || userdata == NULL) return;
  cop = (ConstraintOpPtr) userdata;
  
  if (cop != NULL && cop->descriptorChoice != 0 && cop->descriptorChoice != sdp->choice)
  {
    return;
  }
  
  if (cop != NULL && cop->fsp != NULL && cop->fsp->scp != NULL)
  {
    if (cop->fsp->scp->not_present)
    {
      if (DoesObjectMatchStringConstraint (cop->omtp, cop->aip, (Pointer) sdp, &(cop->ohsd)))
      {
        descriptor_matches = FALSE;
      }
    }
    else
    {
      if (! DoesObjectMatchStringConstraint (cop->omtp, cop->aip, (Pointer) sdp, &(cop->ohsd)))
      {
        descriptor_matches = FALSE;
      }
    }
  }
  if (descriptor_matches)
  {
    (cop->descriptor_action) (sdp, cop->userdata, cop->fsp);
  }
}


static void ConstrainedSourceCallback (SeqEntryPtr sep, ConstraintOpPtr cop)
{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp = NULL;
  SeqAnnotPtr     sap = NULL;
  SeqDescrPtr     sdp = NULL;
  SeqFeatPtr      sfp;
  Boolean         found = FALSE;
  
  if (sep == NULL || sep->data.ptrvalue == NULL || cop == NULL || cop->fsp == NULL)
  {
    return;
  }
   
  if (IS_Bioseq (sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    sdp = bsp->descr;
  }
  else if (IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    sdp = bssp->descr;
  }
  while (sap != NULL && !found)
  {
    if (sap->type == 1)
    {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL && ! found)
      {
        if (sfp->data.choice == SEQFEAT_BIOSRC 
            && DoesOneSourceMatchConstraint (sfp->data.value.ptrvalue, cop->fsp->ccp))
        {
          found = TRUE;
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
  
  while (sdp != NULL && !found)
  {
    if (sdp->choice == Seq_descr_source
        && DoesOneSourceMatchConstraint (sdp->data.ptrvalue, cop->fsp->ccp))
    {
      found = TRUE;
    }
    sdp = sdp->next;
  }
  if (found)
  {
    if (cop->feature_action != NULL)
    {
      cop->omtp = ObjMgrTypeFind (cop->omp, OBJ_SEQFEAT, NULL, NULL);
      if (cop->omtp != NULL)
      {
        VisitBioseqsInSep (sep, cop, BioseqConstrainedFeaturesCallback);
      }
    }
  
    if (cop->descriptor_action != NULL)
    {
      cop->omtp = ObjMgrTypeFind (cop->omp, OBJ_SEQDESC, NULL, NULL);
      if (cop->omtp != NULL)
      {
        VisitDescriptorsInSep (sep, cop, SeqEntryConstrainedDescriptorsCallback);
      }
    }
  }
  else if (bssp != NULL)
  {
    ConstrainedSourceCallback (bssp->seq_set, cop);
  }
  ConstrainedSourceCallback (sep->next, cop);
}


static void ConstrainedBioseqCallback (SeqEntryPtr sep, ConstraintOpPtr cop)
{
  BioseqPtr       bsp = NULL;
  BioseqSetPtr    bssp = NULL;
  SeqFeatPtr      sfp;
  SeqDescPtr      sdp;
  SeqMgrFeatContext fcontext;
  SeqMgrDescContext dcontext;
  
  if (sep == NULL || sep->data.ptrvalue == NULL || cop == NULL || cop->fsp == NULL)
  {
    return;
  }
   
  if (IS_Bioseq (sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
  }
  else if (IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
  }

  if (bsp != NULL && DoesIDListMeetStringConstraint (bsp->id, cop->fsp->id_list))
  {
    if (cop->feature_action != NULL)
    {
      cop->omtp = ObjMgrTypeFind (cop->omp, OBJ_SEQFEAT, NULL, NULL);
      if (cop->omtp != NULL)
      {
        for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
             sfp != NULL;
             sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext))
        {
          SeqEntryConstrainedFeaturesCallback (sfp, cop);
        }
      }
    }
  
    if (cop->descriptor_action != NULL)
    {
      cop->omtp = ObjMgrTypeFind (cop->omp, OBJ_SEQDESC, NULL, NULL);
      if (cop->omtp != NULL)
      {
        for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, 0, &dcontext);
             sdp != NULL;
             sdp = SeqMgrGetNextDescriptor (bsp, sdp, 0, &dcontext))
        {
          SeqEntryConstrainedDescriptorsCallback (sdp, cop);
        }
      }
    }
  }
  else if (bssp != NULL)
  {
    ConstrainedBioseqCallback (bssp->seq_set, cop);
  }
  ConstrainedBioseqCallback (sep->next, cop);
}


extern void 
OperateOnSeqEntryConstrainedObjects 
(SeqEntryPtr           sep,
 FilterSetPtr          fsp,
 FeatureActionProc     feature_action,
 DescriptorActionProc  descriptor_action,
 Uint1                 seqFeatChoice,
 Uint1                 featDefChoice,
 Uint1                 descriptorChoice,
 Pointer               userdata)
{
  ConstraintOpData  cod;
  AsnExpOptPtr      aeop;
  SeqEntryPtr       old_scope;
  BioseqPtr         bsp = NULL;
  SeqFeatPtr        sfp;
  SeqDescrPtr       sdp;
  SeqMgrFeatContext fcontext;
  SeqMgrDescContext dcontext;
  
  if (sep == NULL) return;
  if (feature_action == NULL && descriptor_action == NULL) return;
  
  cod.omp = ObjMgrGet ();
  if (cod.omp == NULL) return;
  
  cod.fsp               = fsp;
  cod.feature_action    = feature_action;
  cod.descriptor_action = descriptor_action;
  cod.userdata          = userdata;
  cod.entityID          = SeqMgrGetEntityIDForSeqEntry (sep);
  cod.seqFeatChoice     = seqFeatChoice;
  cod.featDefChoice     = featDefChoice;
  cod.descriptorChoice  = descriptorChoice;
  
  old_scope = SeqEntrySetScope (sep);

  cod.aip = AsnIoNullOpen ();
  if (fsp == NULL)
  {
    cod.ohsd.scp = NULL;
  }
  else
  {
    cod.ohsd.scp = fsp->scp;
  }
  
  aeop = AsnExpOptNew (cod.aip, NULL, NULL, AsnWriteConstraintCallBack);
  if (aeop != NULL) {
    aeop->user_data = (Pointer) &(cod.ohsd);
  }
  
  if (IS_Bioseq (sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
  }
  
  if (fsp == NULL || (fsp->ccp == NULL && fsp->id_list == NULL)
      || (fsp->ccp != NULL && fsp->ccp->constraint_type == CHOICE_CONSTRAINT_ANY)
      || (fsp->id_list != NULL && bsp != NULL && DoesIDListMeetStringConstraint (bsp->id, fsp->id_list)))
  {
    if (feature_action != NULL)
    {
      cod.omtp = ObjMgrTypeFind (cod.omp, OBJ_SEQFEAT, NULL, NULL);
      if (cod.omtp != NULL)
      {
        if (bsp == NULL)
        {
          VisitBioseqsInSep (sep, &cod, BioseqConstrainedFeaturesCallback);
        }
        else
        {
          sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
          while (sfp != NULL)
          {
            SeqEntryConstrainedFeaturesCallback (sfp, &cod);
            sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
          }
        }
      }
    }
  
    if (descriptor_action != NULL)
    {
      cod.omtp = ObjMgrTypeFind (cod.omp, OBJ_SEQDESC, NULL, NULL);
      if (cod.omtp != NULL)
      {
        if (bsp == NULL)
        {
          VisitDescriptorsInSep (sep, &cod, SeqEntryConstrainedDescriptorsCallback);
        }
        else
        {
          sdp = SeqMgrGetNextDescriptor (bsp, NULL, 0, &dcontext);
          while (sdp != NULL)
          {
            SeqEntryConstrainedDescriptorsCallback (sdp, &cod);
            sdp = SeqMgrGetNextDescriptor (bsp, sdp, 0, &dcontext);
          }
        }
      }
    }
  }
  else if (fsp->id_list != NULL)
  {
    ConstrainedBioseqCallback (sep, &cod);
  }
  else
  {
    ConstrainedSourceCallback (sep, &cod);
  }

  AsnIoClose (cod.aip);
  SeqEntrySetScope (old_scope);
}

extern GetSamplePtr GetSampleNew (void)
{
  GetSamplePtr gsp;
  
  gsp = (GetSamplePtr) MemNew (sizeof (GetSampleData));
  if (gsp != NULL)
  {
    gsp->all_same = TRUE;
    gsp->sample_text = NULL;
    gsp->num_found = 0;
    gsp->feat_dest_list = NULL;
    gsp->descr_dest_list = NULL;
    
    gsp->fieldstring_func = NULL;
    gsp->descrstring_func = NULL;
    gsp->requested_field = NULL;
    gsp->free_vn_proc = NULL;
    gsp->copy_vn_proc = NULL;
  }
  return gsp;
}

extern GetSamplePtr GetSampleFree (GetSamplePtr gsp)
{
  if (gsp != NULL)
  {
    gsp->sample_text = MemFree (gsp->sample_text);
    if (gsp->free_vn_proc != NULL)
    {
      (gsp->free_vn_proc)(gsp->requested_field);
    }
    gsp->requested_field = ValNodeFree (gsp->requested_field);
    gsp->feat_dest_list = ValNodeFree (gsp->feat_dest_list);
    gsp->descr_dest_list = ValNodeFree (gsp->descr_dest_list);
    gsp = MemFree (gsp);
  }
  return gsp;
}

static GetSamplePtr GetSampleCopy (GetSamplePtr gsp)
{
  GetSamplePtr copy_gsp;
  ValNodePtr   copy_vnp;
  
  if (gsp == NULL)
  {
    return NULL;
  }
  
  copy_gsp = (GetSamplePtr) MemNew (sizeof (GetSampleData));
  if (copy_gsp != NULL)
  {
    copy_gsp->fieldstring_func = gsp->fieldstring_func;
    copy_gsp->descrstring_func = gsp->descrstring_func;
    copy_gsp->free_vn_proc = gsp->free_vn_proc;
    copy_gsp->copy_vn_proc = gsp->copy_vn_proc;
    if (copy_gsp->copy_vn_proc != NULL)
    {
      copy_gsp->requested_field = (gsp->copy_vn_proc) (gsp->requested_field);
    }
    copy_gsp->sample_text = StringSave (gsp->sample_text);
    copy_gsp->num_found = gsp->num_found;
    copy_gsp->all_same = gsp->all_same;
    
    copy_vnp = gsp->feat_dest_list;
    copy_gsp->feat_dest_list = NULL;
    while (copy_vnp != NULL)
    {
      ValNodeAddPointer (&(copy_gsp->feat_dest_list), copy_vnp->choice, copy_vnp->data.ptrvalue);
      copy_vnp = copy_vnp->next;
    }
    
    copy_vnp = gsp->descr_dest_list;
    copy_gsp->descr_dest_list = NULL;
    while (copy_vnp != NULL)
    {
      ValNodeAddPointer (&(copy_gsp->descr_dest_list), copy_vnp->choice, copy_vnp->data.ptrvalue);
      copy_vnp = copy_vnp->next;
    }
    
  }
  return copy_gsp;
}

static ValNodePtr AddValNodeLists (ValNodePtr vnp1, ValNodePtr vnp2)
{
  ValNodePtr new_list = NULL, copy_vnp, check_vnp;
  
  /* copy first list */
  copy_vnp = vnp1;
  while (copy_vnp != NULL)
  {
    ValNodeAddPointer (&(new_list), copy_vnp->choice, copy_vnp->data.ptrvalue);
    copy_vnp = copy_vnp->next;
  }
      
  copy_vnp = vnp2;
  while (copy_vnp != NULL)
  {
    /* check to see if destination is already in the list */
      check_vnp = new_list;
       while (check_vnp != NULL && check_vnp->data.ptrvalue != copy_vnp->data.ptrvalue)
       {
      check_vnp = check_vnp->next;
    }
      
    if (check_vnp == NULL)
    {
      /* add to list if not present */
      ValNodeAddPointer (&new_list, copy_vnp->choice, copy_vnp->data.ptrvalue);
    }
    else if (check_vnp->choice < 250)
    {
      /* do not increase choice beyond 250 */
      check_vnp->choice ++;
    }
          
    copy_vnp = copy_vnp->next;
  }
  return new_list;
}

static GetSamplePtr GetSampleAdd (GetSamplePtr gsp1, GetSamplePtr gsp2)
{
  GetSamplePtr sum_gsp;
  
  if (gsp1 == NULL && gsp2 == NULL)
  {
    sum_gsp =  NULL;
  }
  else if (gsp1 == NULL)
  {
    sum_gsp =  GetSampleCopy (gsp2);
  }
  else if (gsp2 == NULL)
  {
    sum_gsp =  GetSampleCopy (gsp1);
  }
  else
  {
    sum_gsp = (GetSamplePtr) MemNew (sizeof (GetSampleData));
    if (sum_gsp != NULL)
    {
      sum_gsp->fieldstring_func = gsp1->fieldstring_func;
      sum_gsp->descrstring_func = gsp1->descrstring_func;
      sum_gsp->requested_field = gsp1->requested_field;
      gsp1->requested_field = NULL;
      sum_gsp->sample_text = StringSave (gsp1->sample_text);
      sum_gsp->num_found = gsp1->num_found + gsp2->num_found;
      if (gsp1->num_found == 0)
      {
        sum_gsp->all_same = gsp2->all_same;
      }
      else if (gsp2->num_found == 0)
      {
        sum_gsp->all_same = gsp1->all_same;
      }
      else if (!gsp1->all_same || ! gsp2->all_same)
      {
        sum_gsp->all_same = FALSE;
      }
      else if (StringCmp (gsp1->sample_text, gsp2->sample_text) == 0)
      {
        sum_gsp->all_same = TRUE;
      }
      else
      {
        sum_gsp->all_same = FALSE;
      }
      
      /* combine destination lists */
      sum_gsp->feat_dest_list = AddValNodeLists (gsp1->feat_dest_list, gsp2->feat_dest_list);
            
    }
  }
  return sum_gsp;
}


static void GetSampleFeatureCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  GetSamplePtr cp;
  CharPtr      str;
  ValNodePtr   check_vnp;
  
  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  
  cp = (GetSamplePtr) userdata;
  if (cp->fieldstring_func == NULL)
  {
    return;
  }
  
  check_vnp = cp->feat_dest_list;
  while (check_vnp != NULL && check_vnp->data.ptrvalue != sfp)
  {
    check_vnp = check_vnp->next;
  }
  if (check_vnp == NULL)
  {
    ValNodeAddPointer (&(cp->feat_dest_list), 1, sfp);
  }
  else
  {
    if (check_vnp->choice < 250)
    {
      check_vnp->choice ++;
    }
    return;
  }
  
  str = cp->fieldstring_func (sfp, cp->requested_field, NULL);
  if (!StringHasNoText (str) || (fsp != NULL && fsp->scp != NULL && fsp->scp->not_present))
  {
    cp->num_found ++;
    if (cp->sample_text == NULL)
    {
      cp->sample_text = str;
    }
    else
    {
      if (StringCmp (str, cp->sample_text) != 0)
      {
        cp->all_same = FALSE;
      }
      str = MemFree (str);
    }
  }
}

static void GetSampleDescriptorCallback (SeqDescrPtr sdp, Pointer userdata, FilterSetPtr fsp)
{
  GetSamplePtr cp;
  CharPtr      str;
  ValNodePtr   check_vnp;
  
  if (sdp == NULL || userdata == NULL)
  {
    return;
  }
  
  cp = (GetSamplePtr) userdata;
  if (cp->descrstring_func == NULL)
  {
    return;
  }
  
  check_vnp = cp->descr_dest_list;
  while (check_vnp != NULL && check_vnp->data.ptrvalue != sdp)
  {
    check_vnp = check_vnp->next;
  }
  if (check_vnp == NULL)
  {
    ValNodeAddPointer (&(cp->descr_dest_list), 1, sdp);
  }
  else
  {
    if (check_vnp->choice < 250)
    {
      check_vnp->choice ++;
    }
    return;
  }

  str = cp->descrstring_func (sdp, cp->requested_field, NULL);
  if (!StringHasNoText (str))
  {
    cp->num_found ++;
    if (cp->sample_text == NULL)
    {
      cp->sample_text = str;
    }
    else
    {
      if (StringCmp (str, cp->sample_text) != 0)
      {
        cp->all_same = FALSE;
      }
      str = MemFree (str);
    }
  }
}

static GetSamplePtr CheckForExistingTextInSeqEntry
(SeqEntryPtr              sep,
 ValNodePtr               requested_field,
 GetFeatureFieldString    fieldstring_func,
 GetDescriptorFieldString descrstring_func,
 FreeValNodeProc          free_vn_proc,
 CopyValNodeDataProc      copy_vn_proc,
 FilterSetPtr             fsp,
 Uint1                    seqfeat_choice,
 Uint1                    featdef_choice,
 Uint1                    descr_choice)
{
  GetSamplePtr gsp;
  
  gsp = (GetSamplePtr) MemNew (sizeof (GetSampleData));
  if (gsp != NULL)
  {
    gsp->sample_text = NULL;
    gsp->fieldstring_func = fieldstring_func;
    gsp->descrstring_func = descrstring_func;
    gsp->free_vn_proc = free_vn_proc;
    gsp->copy_vn_proc = copy_vn_proc;
    gsp->requested_field = (gsp->copy_vn_proc) (requested_field);
    gsp->num_found = 0;
    gsp->all_same = TRUE;
    OperateOnSeqEntryConstrainedObjects (sep, fsp, GetSampleFeatureCallback, 
                                         GetSampleDescriptorCallback,
                                         seqfeat_choice, featdef_choice,
                                         descr_choice, gsp);
  }
  return gsp;  
}

static GetSamplePtr CheckForExistingText
(Uint2                    entityID,
 ValNodePtr               requested_field,
 GetFeatureFieldString    fieldstring_func,
 GetDescriptorFieldString descrstring_func,
 FreeValNodeProc          free_vn_proc,
 CopyValNodeDataProc      copy_vn_proc,
 FilterSetPtr             fsp,
 Uint1                    seqfeat_choice,
 Uint1                    featdef_choice,
 Uint1                    descr_choice)
{
  SeqEntryPtr  sep;

  sep = GetTopSeqEntryForEntityID (entityID);
  return CheckForExistingTextInSeqEntry (sep, requested_field, 
                                         fieldstring_func, descrstring_func,
                                         free_vn_proc, copy_vn_proc, fsp,
                                         seqfeat_choice, featdef_choice, descr_choice);
}

extern CharPtr HandleApplyValue (CharPtr orig_text, ApplyValuePtr avp)
{
  CharPtr new_str, cp_found;
  Int4    found_len, replace_len, new_len;
  
  if (avp == NULL)
  {
    return orig_text;
  }
  else if (StringHasNoText (orig_text))
  {
    if (StringHasNoText (avp->text_to_replace))
    {
      MemFree (orig_text);
      return StringSave (avp->new_text);
    }
    else
    {
      return orig_text;
    }
  }
  else if (StringHasNoText (avp->text_to_replace))
  {
    return HandleExistingText (orig_text, StringSave (avp->new_text), avp->etp);
  }
  else
  {
    found_len = StringLen (avp->text_to_replace);
    replace_len = StringLen (avp->new_text);
    cp_found = StringISearch (orig_text, avp->text_to_replace);
    if (avp->where_to_replace == EditApplyFindLocation_beginning
        && cp_found != orig_text) {
      cp_found = NULL;
    } 
    while (cp_found != NULL)
    {
      if (avp->where_to_replace == EditApplyFindLocation_end
          && cp_found != orig_text + StringLen (orig_text) - found_len) {
        cp_found = StringISearch (cp_found + found_len, avp->text_to_replace);
      } else {
        new_len = StringLen (orig_text) + 1 - found_len + replace_len;
        new_str = (CharPtr) MemNew (new_len * sizeof (Char));
        if (new_str != NULL)
        {
          if (cp_found != orig_text)
          {
            StringNCpy (new_str, orig_text, cp_found - orig_text);
          }
          StringCat (new_str, avp->new_text);
          StringCat (new_str, cp_found + found_len);
          cp_found = new_str + (cp_found - orig_text) + replace_len;
          orig_text = MemFree (orig_text);
          orig_text = new_str;
        }
        cp_found = StringISearch (cp_found, avp->text_to_replace);
      }
    }
    return orig_text;
  }
}

extern ValNodePtr ApplyValueToValNodeStringList (ValNodePtr list, Int2 choice, ApplyValuePtr avp)
{
  ValNodePtr vnp;
  
  if (avp == NULL)
  {
    return NULL;
  }
  
  if (!StringHasNoText (avp->text_to_replace))
  {
    for (vnp = list; vnp != NULL; vnp = vnp->next)
    {
      vnp->data.ptrvalue = HandleApplyValue (vnp->data.ptrvalue, avp);
    }
  }
  else
  {
    ValNodeAddPointer (&list, choice, StringSave (avp->new_text));
  }
  return list;
}



static void MakeSemicolonStringsIntoSeparateItems (ValNodePtr PNTR list)
{
  ValNodePtr prev = NULL, vnp_next, vnp, vnp_new;
  CharPtr    cp;

  if (list == NULL || *list == NULL) return;

  vnp = *list;
  while (vnp != NULL)
  {
    vnp_next = vnp->next;
    if (StringHasNoText (vnp->data.ptrvalue))
    {
      if (prev == NULL)
      {
        *list = vnp_next;
      }
      else
      {
        prev->next = vnp_next;
      }
      vnp->next = NULL;
      vnp = ValNodeFreeData (vnp);
    }
    else if ((cp = StringChr (vnp->data.ptrvalue, ';')) != NULL)
    {
      vnp_new = ValNodeNew (NULL);
      vnp_new->choice = vnp->choice;
      vnp_new->data.ptrvalue = StringSave (cp + 1);
      vnp_new->next = vnp_next;
      *cp = 0;
      vnp->next = vnp_new;
      vnp_next = vnp_new;
      prev = vnp;
    }
    else
    {
      prev = vnp;
    }
    vnp = vnp_next;
  }
}


extern ValNodePtr ApplyValueToValNodeStringListAsText (ValNodePtr list, Int2 choice, ApplyValuePtr avp)
{
  ValNodePtr vnp;
  
  if (avp == NULL)
  {
    return NULL;
  }
  
  if (!StringHasNoText (avp->text_to_replace))
  {
    for (vnp = list; vnp != NULL; vnp = vnp->next)
    {
      vnp->data.ptrvalue = HandleApplyValue (vnp->data.ptrvalue, avp);
    }
  }
  else
  {
    if (avp->etp == NULL
        || avp->etp->existing_text_choice == eExistingTextChoiceReplaceOld)
    {
      list = ValNodeFree (list);
      ValNodeAddPointer (&list, choice, StringSave (avp->new_text));
    }
    else if (avp->etp->existing_text_choice == eExistingTextChoiceAppendSemi)
    {
      ValNodeAddPointer (&list, choice, StringSave (avp->new_text));
    }
    else if (avp->etp->existing_text_choice == eExistingTextChoicePrefixSemi)
    {
      vnp = ValNodeNew(NULL);
      vnp->choice = (Uint1) choice;
      vnp->data.ptrvalue = StringSave (avp->new_text);
      vnp->next = list;
      list = vnp;
    }
    else if (avp->etp->existing_text_choice == eExistingTextChoiceLeaveOld)
    {
      if (list == NULL)
      {
        ValNodeAddPointer (&list, choice, StringSave (avp->new_text));
      }
    } else {
      for (vnp = list; vnp != NULL; vnp = vnp->next)
      {
        vnp->data.ptrvalue = HandleApplyValue (vnp->data.ptrvalue, avp);
      }
    }
  }

  if (avp->etp == NULL || avp->etp->existing_text_choice != eExistingTextChoiceLeaveOld) 
  {
    MakeSemicolonStringsIntoSeparateItems (&list);
  }
  return list;
}


typedef struct parsefielddialog
{
  DIALOG_MESSAGE_BLOCK
  DialoG                   parse_field_type;
  DialoG                   biosrc_string_choice;
  DialoG                   source_qual_choice;
  DialoG                   gene_field;
  DialoG                   rna_subtype;
  DialoG                   rna_field;
  DialoG                   protein_field;
  DialoG                   import_feature;
  DialoG                   import_qual;
  DialoG                   feature;
  GrouP                    dbxref_grp;
  TexT                     dbxref_db;
  GrouP                    feat_or_desc;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
  Boolean                  is_search_field;
} ParseFieldDialogData, PNTR ParseFieldDialogPtr;

static void ChangeParseFieldType (Pointer data)
{
  ParseFieldDialogPtr dlg;
  ValNodePtr          vnp;

  dlg = (ParseFieldDialogPtr) data;
  if (dlg == NULL)
  {
    return;
  }
  vnp = DialogToPointer (dlg->parse_field_type);
  Hide (dlg->biosrc_string_choice);
  Hide (dlg->feat_or_desc);
  Hide (dlg->source_qual_choice);
  Hide (dlg->gene_field);
  Hide (dlg->rna_subtype);
  Hide (dlg->rna_field);
  Hide (dlg->protein_field);
  Hide (dlg->import_feature);
  Hide (dlg->import_qual);
  Hide (dlg->feature);
  Hide (dlg->dbxref_grp);
  if (vnp != NULL)
  {
    switch (vnp->choice)
    {
      case PARSE_FIELD_BIOSRC_STRING:
        Show (dlg->biosrc_string_choice);
        Show (dlg->feat_or_desc);
        break;
      case PARSE_FIELD_SOURCE_QUAL:
        Show (dlg->source_qual_choice);
        Show (dlg->feat_or_desc);
        break;
      case PARSE_FIELD_GENE_FIELD:
        Show (dlg->gene_field);
        break;
      case PARSE_FIELD_RNA_FIELD:
        Show (dlg->rna_subtype);
        Show (dlg->rna_field);
        break;
      case PARSE_FIELD_PROTEIN_FIELD:
        Show (dlg->protein_field);
        break;
      case PARSE_FIELD_IMPORT_QUAL:
        Show (dlg->import_feature);
        Show (dlg->import_qual);
        break;
      case PARSE_FIELD_FEATURE_NOTE:
        Show (dlg->feature);
        break;
      case PARSE_FIELD_DBXREF:
        Show (dlg->dbxref_grp);
        break;
    }
  }
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify)(dlg->change_userdata);
  }
}

static void ResetParseFieldDialog (ParseFieldDialogPtr dlg)
{
  ValNode vn;

  if (dlg == NULL)
  {
    return;
  }
  vn.choice = PARSE_FIELD_DEFLINE;
  vn.data.ptrvalue = NULL;
  vn.next = NULL;
  PointerToDialog (dlg->parse_field_type, &vn);
  ChangeParseFieldType (dlg);
}

static void ParseFieldToDialog (DialoG d, Pointer data)
{
  ParseFieldDialogPtr dlg;
  ParseFieldPtr       parse_field;
  Int4                parse_field_type;
  ValNode             vn;

  dlg = (ParseFieldDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  ResetParseFieldDialog (dlg);
  parse_field = (ParseFieldPtr) data;
  if (parse_field != NULL)
  {
    parse_field_type = parse_field->parse_field_type;
    if (parse_field_type < PARSE_FIELD_DEFLINE 
        || (dlg->is_search_field && parse_field_type > SEARCH_FIELD_PUBLICATION)
        || (!dlg->is_search_field && parse_field_type > MAX_PARSE_FIELD_TYPE))
    {
      parse_field_type = PARSE_FIELD_DEFLINE;
    }
    vn.choice = parse_field_type;
    vn.data.ptrvalue = NULL;
    vn.next = NULL;
    PointerToDialog (dlg->parse_field_type, &vn);
    switch (parse_field_type)
    {
      case PARSE_FIELD_BIOSRC_STRING:
        PointerToDialog (dlg->biosrc_string_choice, parse_field->feature_field);
        if (parse_field->do_desc && parse_field->do_feat)
        {
          SetValue (dlg->feat_or_desc, 1);
        }
        else if (parse_field->do_desc)
        {
          SetValue (dlg->feat_or_desc, 2);
        }
        else if (parse_field->do_feat)
        {
          SetValue (dlg->feat_or_desc, 3);
        }
        else
        {
          SetValue (dlg->feat_or_desc, 1);
        }
        break;
      case PARSE_FIELD_SOURCE_QUAL:
        PointerToDialog (dlg->source_qual_choice, parse_field->feature_field);
        if (parse_field->do_desc && parse_field->do_feat)
        {
          SetValue (dlg->feat_or_desc, 1);
        }
        else if (parse_field->do_desc)
        {
          SetValue (dlg->feat_or_desc, 2);
        }
        else if (parse_field->do_feat)
        {
          SetValue (dlg->feat_or_desc, 3);
        }
        else
        {
          SetValue (dlg->feat_or_desc, 1);
        }
        break;
      case PARSE_FIELD_DBXREF:
        if (parse_field->feature_field == NULL || parse_field->feature_field->data.ptrvalue == NULL)
        {
          SetTitle (dlg->dbxref_db, "");
        }
        else
        {
          SetTitle (dlg->dbxref_db, parse_field->feature_field->data.ptrvalue);
        }
        if (parse_field->do_desc && parse_field->do_feat)
        {
          SetValue (dlg->feat_or_desc, 1);
        }
        else if (parse_field->do_desc)
        {
          SetValue (dlg->feat_or_desc, 2);
        }
        else if (parse_field->do_feat)
        {
          SetValue (dlg->feat_or_desc, 3);
        }
        else
        {
          SetValue (dlg->feat_or_desc, 1);
        }
        break;
      case PARSE_FIELD_GENE_FIELD:
        PointerToDialog (dlg->gene_field, parse_field->feature_field);
        break;
      case PARSE_FIELD_RNA_FIELD:
        PointerToDialog (dlg->rna_subtype, parse_field->feature_subtype);
        PointerToDialog (dlg->rna_field, parse_field->feature_field);
        break;
      case PARSE_FIELD_PROTEIN_FIELD:
        PointerToDialog (dlg->protein_field, parse_field->feature_field);
        break;
      case PARSE_FIELD_IMPORT_QUAL:
        PointerToDialog (dlg->import_feature, parse_field->feature_subtype);
        PointerToDialog (dlg->import_qual, parse_field->feature_field);
        break;
      case PARSE_FIELD_FEATURE_NOTE:
        PointerToDialog (dlg->feature, parse_field->feature_field);
        break;
      case SEARCH_FIELD_PUBLICATION:
        /* nothing to do here */
        break;
    }
  }
  ChangeParseFieldType (dlg);
}

static Pointer DialogToParseField (DialoG d)
{
  ParseFieldDialogPtr dlg;
  ParseFieldPtr       parse_field;
  Int4                val;
  ValNodePtr          vnp;

  dlg = (ParseFieldDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  vnp = DialogToPointer (dlg->parse_field_type);
  if (vnp == NULL)
  {
    return NULL;
  }

  parse_field = (ParseFieldPtr) MemNew (sizeof (ParseFieldData));
  if (parse_field != NULL)
  {
    parse_field->parse_field_type = vnp->choice;
    parse_field->feature_field = NULL;
    switch (parse_field->parse_field_type)
    {
      case PARSE_FIELD_BIOSRC_STRING:
        parse_field->feature_field = DialogToPointer (dlg->biosrc_string_choice);
        val = GetValue (dlg->feat_or_desc);
        switch (val)
        {
          case 2:
            parse_field->do_desc = TRUE;
            parse_field->do_feat = FALSE;
            break;
          case 3:
            parse_field->do_desc = FALSE;
            parse_field->do_feat = TRUE;
            break;
          case 1:
          default:
            parse_field->do_desc = TRUE;
            parse_field->do_feat = TRUE;
            break;
        }
        break;
      case PARSE_FIELD_SOURCE_QUAL:
        parse_field->feature_field = DialogToPointer (dlg->source_qual_choice);
        val = GetValue (dlg->feat_or_desc);
        switch (val)
        {
          case 2:
            parse_field->do_desc = TRUE;
            parse_field->do_feat = FALSE;
            break;
          case 3:
            parse_field->do_desc = FALSE;
            parse_field->do_feat = TRUE;
            break;
          case 1:
          default:
            parse_field->do_desc = TRUE;
            parse_field->do_feat = TRUE;
            break;
        }
        break;
      case PARSE_FIELD_DBXREF:        
        if (!TextHasNoText (dlg->dbxref_db))
        {
          ValNodeAddPointer (&(parse_field->feature_field), 0, SaveStringFromText (dlg->dbxref_db));
        }
        val = GetValue (dlg->feat_or_desc);
        switch (val)
        {
          case 2:
            parse_field->do_desc = TRUE;
            parse_field->do_feat = FALSE;
            break;
          case 3:
            parse_field->do_desc = FALSE;
            parse_field->do_feat = TRUE;
            break;
          case 1:
          default:
            parse_field->do_desc = TRUE;
            parse_field->do_feat = TRUE;
            break;
        }
        break;
      case PARSE_FIELD_GENE_FIELD:
        parse_field->feature_field = DialogToPointer (dlg->gene_field);
        break;
      case PARSE_FIELD_RNA_FIELD:
        parse_field->feature_subtype = DialogToPointer (dlg->rna_subtype);
        parse_field->feature_field = DialogToPointer (dlg->rna_field);
        break;
      case PARSE_FIELD_PROTEIN_FIELD:
        parse_field->feature_field = DialogToPointer (dlg->protein_field);
        break;
      case PARSE_FIELD_IMPORT_QUAL:
        parse_field->feature_subtype = DialogToPointer (dlg->import_feature);
        parse_field->feature_field = DialogToPointer (dlg->import_qual);
        break;
      case PARSE_FIELD_FEATURE_NOTE:
        parse_field->feature_field = DialogToPointer (dlg->feature);
        break;
      case SEARCH_FIELD_PUBLICATION:
        /* nothing to do here */
        break;
    }
  }
  vnp = ValNodeFreeData (vnp);
  return parse_field;
}

static void ParseFieldMessage (DialoG d, Int2 mssg)

{
  ParseFieldDialogPtr dlg;

  dlg = (ParseFieldDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        ResetParseFieldDialog (dlg);        
        break;
      case VIB_MSG_ENTER :
        Select (dlg->parse_field_type);
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestParseFieldDialog (DialoG d)

{
  ValNodePtr          head = NULL;
  ParseFieldDialogPtr dlg;
  Int4                parse_field_type;
  ValNodePtr          vnp;

  dlg = (ParseFieldDialogPtr) GetObjectExtra (d);
  
  if (dlg != NULL)
  {
    vnp = DialogToPointer (dlg->parse_field_type);
    if (vnp == NULL) 
    {
      head = AddStringToValNodeChain (head, "destination type", 1);
    }
    else
    {
      parse_field_type = vnp->choice;
      vnp = ValNodeFreeData (vnp);
      switch (parse_field_type)
      {
        case PARSE_FIELD_DEFLINE:
        case PARSE_FIELD_CDS_COMMENT:
        case PARSE_FIELD_COMMENT_DESC:
          /* don't need any extra information */
          break;
        case PARSE_FIELD_BIOSRC_STRING:
          vnp = DialogToPointer (dlg->biosrc_string_choice);
          if (vnp == NULL)
          {
            head = AddStringToValNodeChain (head, "biosrc string type", 1);
          }
          ValNodeFree (vnp);
          break;
        case PARSE_FIELD_SOURCE_QUAL:
          vnp = DialogToPointer (dlg->source_qual_choice);
          if (vnp == NULL)
          {
            head = AddStringToValNodeChain (head, "source qual type", 1);
          }
          ValNodeFree (vnp);
          break;
        case PARSE_FIELD_GENE_FIELD:
          vnp = DialogToPointer (dlg->gene_field);
          if (vnp == NULL || vnp->data.intvalue == FEATUREFIELD_NONE)
          {
            head = AddStringToValNodeChain (head, "gene field", 1);
          }
          ValNodeFree (vnp);
          break;
        case PARSE_FIELD_RNA_FIELD:
          vnp = DialogToPointer (dlg->rna_subtype);
          if (vnp == NULL)
          {
            head = AddStringToValNodeChain (head, "RNA subtype", 1);
          }
          ValNodeFree (vnp);
          vnp = DialogToPointer (dlg->rna_field);
          if (vnp == NULL || vnp->data.intvalue == FEATUREFIELD_NONE)
          {
            head = AddStringToValNodeChain (head, "RNA field", 1);
          }
          ValNodeFree (vnp);          
          break;
        case PARSE_FIELD_PROTEIN_FIELD:
          vnp = DialogToPointer (dlg->protein_field);
          if (vnp == NULL || vnp->data.intvalue == FEATUREFIELD_NONE)
          {
            head = AddStringToValNodeChain (head, "protein field", 1);
          }
          ValNodeFree (vnp);
          break;
        case PARSE_FIELD_IMPORT_QUAL:
          vnp = DialogToPointer (dlg->import_qual);
          if (vnp == NULL)
          {
            ValNodeAddPointer (&head, 1, "import qualifier");
          }
          ValNodeFree (vnp);
          vnp = DialogToPointer (dlg->import_feature);
          if (vnp == NULL)
          {
            ValNodeAddPointer (&head, 1, "import feature");
          }
          ValNodeFree (vnp);
          break;
        case PARSE_FIELD_FEATURE_NOTE:
          vnp = DialogToPointer (dlg->feature);
          if (vnp == NULL)
          {
            ValNodeAddPointer (&head, 1, "feature");
            ValNodeFree (vnp);
          }
          break;
        case PARSE_FIELD_DBXREF:
          if (TextHasNoText (dlg->dbxref_db))
          {
            ValNodeAddPointer (&head, 1, "specify database");
          }
          break;
        case SEARCH_FIELD_PUBLICATION:
          if (!dlg->is_search_field) 
          {
            head = AddStringToValNodeChain (head, "field type", 1);
          }
          break;
        default:
          head = AddStringToValNodeChain (head, "field type", 1);
          break;
      }
    }
  }      

  return head;
}


static ValNodePtr GetParseFieldDestList (Boolean include_pub, Boolean include_dbxref)
{
  ValNodePtr list = NULL;

  ValNodeAddPointer (&list, PARSE_FIELD_DEFLINE, StringSave ("Definition Line"));
  ValNodeAddPointer (&list, PARSE_FIELD_BIOSRC_STRING, StringSave ("Biosource"));
  ValNodeAddPointer (&list, PARSE_FIELD_SOURCE_QUAL, StringSave ("Source Qualifier"));
  ValNodeAddPointer (&list, PARSE_FIELD_GENE_FIELD, StringSave ("Gene Field"));
  ValNodeAddPointer (&list, PARSE_FIELD_RNA_FIELD, StringSave ("RNA Field"));
  ValNodeAddPointer (&list, PARSE_FIELD_CDS_COMMENT, StringSave ("CDS Comment"));
  ValNodeAddPointer (&list, PARSE_FIELD_PROTEIN_FIELD, StringSave ("Protein Field"));
  ValNodeAddPointer (&list, PARSE_FIELD_IMPORT_QUAL, StringSave ("Import Feature"));
  ValNodeAddPointer (&list, PARSE_FIELD_FEATURE_NOTE, StringSave ("Feature Note"));
  ValNodeAddPointer (&list, PARSE_FIELD_COMMENT_DESC, StringSave ("Comment Descriptor"));
  if (include_dbxref)
  {
    ValNodeAddPointer (&list, PARSE_FIELD_DBXREF, StringSave ("Dbxref"));
  }
  if (include_pub)
  {
    ValNodeAddPointer (&list, SEARCH_FIELD_PUBLICATION, StringSave ("Publication"));
  }
  return list;  
}


static void ChangeParseFieldDestText (TexT t)
{
  ParseFieldDialogPtr dlg;

  dlg = (ParseFieldDialogPtr) GetObjectExtra (t);
  if (dlg != NULL && dlg->change_notify != NULL)
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


extern DialoG ParseFieldDestDialogEx 
(GrouP                    h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  is_search_field,
 Boolean                  include_dbxref)
{
  ParseFieldDialogPtr dlg;
  GrouP               p, m, g;
  GrouP               subgrp;
  ValNodePtr          choice_list;
  ValNode             vn;
  
  dlg = (ParseFieldDialogPtr) MemNew (sizeof (ParseFieldDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = ParseFieldToDialog;
  dlg->fromdialog = DialogToParseField;
  dlg->dialogmessage = ParseFieldMessage;
  dlg->testdialog = TestParseFieldDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->is_search_field = is_search_field;
  
  m = HiddenGroup (p, 2, 0, NULL);

  choice_list = GetParseFieldDestList (is_search_field, include_dbxref);
  dlg->parse_field_type = ValNodeSelectionDialog (m, choice_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "parse destination", 
                                ChangeParseFieldType, dlg, FALSE);

  subgrp = HiddenGroup (m, 0, 0, NULL);
  dlg->biosrc_string_choice = BioSourceStringDialog (subgrp, FALSE, change_notify, change_userdata);
  Hide (dlg->biosrc_string_choice);
  dlg->source_qual_choice = SourceQualTypeSelectionDialog (subgrp, FALSE, change_notify, change_userdata);
  Hide (dlg->source_qual_choice);
  dlg->gene_field = GeneFieldSelectionDialog (subgrp, FALSE, change_notify, change_userdata);
  Hide (dlg->gene_field);
  g = HiddenGroup (subgrp, 2, 0, NULL);
  dlg->rna_subtype = RNASubtypeSelectionDialog (g, FALSE, NULL, change_notify, change_userdata);
  Hide (dlg->rna_subtype);
  dlg->rna_field = RNAFieldSelectionDialog (g, FALSE, 
                                            change_notify, 
                                            change_userdata);
  Hide (dlg->rna_field);
  dlg->protein_field = ProteinFieldSelectionDialog (subgrp, FALSE, change_notify, change_userdata);
  Hide (dlg->protein_field);
  g = HiddenGroup (subgrp, 2, 0, NULL);
  dlg->import_feature = FeatureSelectionDialog (g, TRUE, change_notify, change_userdata);
  Hide (dlg->import_feature);
  dlg->import_qual = GBQualSelectionDialog (g, FALSE, change_notify, change_userdata);
  Hide (dlg->import_qual);
  dlg->feature = FeatureSelectionDialog (subgrp, TRUE, change_notify, change_userdata);
  Hide (dlg->feature);
  if (include_dbxref)
  {
    dlg->dbxref_grp = HiddenGroup (subgrp, 2, 0, NULL);
    StaticPrompt (dlg->dbxref_grp, "Database:", 0, 0, programFont, 'r');
    dlg->dbxref_db = DialogText (dlg->dbxref_grp, "", 5, ChangeParseFieldDestText);
    SetObjectExtra (dlg->dbxref_db, dlg, NULL);
    Hide (dlg->dbxref_grp);
  }
  else
  {
    dlg->dbxref_db = NULL;
  }

  vn.choice = PARSE_FIELD_DEFLINE;
  vn.data.ptrvalue = NULL;
  vn.next = NULL;
  PointerToDialog (dlg->parse_field_type, &vn);

  dlg->feat_or_desc = HiddenGroup (p, 3, 0, NULL);
  RadioButton (dlg->feat_or_desc, "Descriptors and Features");
  RadioButton (dlg->feat_or_desc, "Descriptors Only");
  RadioButton (dlg->feat_or_desc, "Features Only");
  SetValue (dlg->feat_or_desc, 1);
  Hide (dlg->feat_or_desc);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) m, (HANDLE) dlg->feat_or_desc, NULL);

  return (DialoG) p;  
}

extern DialoG ParseFieldDestDialog 
(GrouP                    h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  return ParseFieldDestDialogEx (h, change_notify, change_userdata, FALSE, FALSE);
}
/* Search Field dialog is identical to ParseFieldDestDialog except that it allows
 * an additional option, Publication.
 */
extern DialoG SearchFieldDialog 
(GrouP                    h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  return ParseFieldDestDialogEx (h, change_notify, change_userdata, TRUE, FALSE);
}


typedef enum
{
 eParseFieldSrcDefLine = 1,
 eParseFieldSrcGenBankFlatFile,
 eParseFieldSrcLocalID,
 eParseFieldSrcTaxName,
 eParseFieldSrcSourceQual,
 eParseFieldSrcComment,
 eParseFieldSrcBankitComment,
 eParseFieldSrcStructuredComment,
 eParseFieldSrcTaxNameAfterNomial,
 eParseFieldSrcFileID
} EParseFieldSrc;

const Int4 MAX_PARSE_FIELD = eParseFieldSrcFileID;

static CharPtr parse_src_field_list [] = 
{
  "Definition Line", "GenBank Flatfile", "Local ID", "Organism Name", "Source Qualifier", "Comment", "Bankit Comment", 
  "Structured Comment", "Taxname after Binomial/Trinomial", "File ID"
};

static int num_parse_src_fields = sizeof (parse_src_field_list) / sizeof (CharPtr);

typedef struct parsefieldsrc 
{
  Int4              parse_field;
  SourceQualDescPtr sqdp;
  CharPtr           comment_field;
} ParseFieldSrcData, PNTR ParseFieldSrcPtr;

static ParseFieldSrcPtr ParseFieldSrcFree (ParseFieldSrcPtr pfsp)
{
  if (pfsp != NULL) {
    pfsp->comment_field = MemFree (pfsp->comment_field);
    pfsp->sqdp = MemFree (pfsp->sqdp);
    pfsp = MemFree (pfsp);
  }
  return pfsp;
}

typedef struct parsefieldsrcdlg
{
  DIALOG_MESSAGE_BLOCK
  DialoG     parse_field;
  PopuP      comment_field;
  ValNodePtr comment_field_list;
  DialoG     src_qual;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;      
} ParseFieldSrcDlgData, PNTR ParseFieldSrcDlgPtr;

static void CleanupParseFieldSrcDialog (GraphiC g, VoidPtr data)
{
  ParseFieldSrcDlgPtr dlg;

  dlg = (ParseFieldSrcDlgPtr) data;
  if (dlg != NULL) {
    dlg->comment_field_list = ValNodeFreeData(dlg->comment_field_list);
  }
  StdCleanupExtraProc (g, data);
}

static void ResetParseFieldSourceDialog (ParseFieldSrcDlgPtr dlg)
{
  if (dlg == NULL) {
      return;
  }
  SendMessageToDialog (dlg->parse_field, VIB_MSG_INIT);
  SetValue (dlg->comment_field, 0);
  Hide (dlg->comment_field);
}

static void ParseFieldSourceDialogMessage (DialoG d, Int2 mssg)

{
  ParseFieldSrcDlgPtr dlg;

  dlg = (ParseFieldSrcDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        /* reset list */
        ResetParseFieldSourceDialog (dlg);
        break;
      default :
        break;
    }
  }
}

static void ParseFieldChange (Pointer userdata) 
{
  ParseFieldSrcDlgPtr dlg;
  ValNodePtr          vnp;
  
  dlg = (ParseFieldSrcDlgPtr) userdata;
  if (dlg == NULL) return;
   
  vnp = DialogToPointer (dlg->parse_field);
  if (vnp == NULL || vnp->data.intvalue != eParseFieldSrcStructuredComment) {
    Hide (dlg->comment_field);
  } else {
    Show (dlg->comment_field);
  }
  if (vnp == NULL || vnp->data.intvalue != eParseFieldSrcSourceQual) {
    Hide (dlg->src_qual);
  } else {
    Show (dlg->src_qual);
  }
  
  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }  
}

static Int4 FindCommentFieldPosition (CharPtr comment_field, ValNodePtr comment_list)
{
  Int4       pos = 1;
  ValNodePtr vnp;
  
  if (StringHasNoText (comment_field) || comment_list == NULL) {
    return 0;
  }
  
  vnp = comment_list;
  while (vnp != NULL && StringCmp (comment_field, vnp->data.ptrvalue) != 0) {
    vnp = vnp->next;
    pos++;
  }
  if (vnp == NULL) {
    pos = 0;
  }
  return pos;
}

static void DataToParseFieldSrcDialog(DialoG d, Pointer data)
{
  ParseFieldSrcDlgPtr dlg;
  ParseFieldSrcPtr    pfsp;
  ValNode             vn;

  dlg = (ParseFieldSrcDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    pfsp = (ParseFieldSrcPtr) data;
    if (pfsp == NULL) {
      ResetParseFieldSourceDialog (dlg);
    } else {
      vn.data.intvalue = pfsp->parse_field;
      vn.next = NULL;
      vn.choice = 0;
      PointerToDialog (dlg->parse_field, &vn);
      SetValue (dlg->comment_field, FindCommentFieldPosition (pfsp->comment_field, dlg->comment_field_list));
      vn.data.ptrvalue = pfsp->sqdp;
      PointerToDialog (dlg->src_qual, &vn);
    }
  }
  ParseFieldChange (dlg);     
}

static Pointer ParseFieldSrcDialogToData(DialoG d)
{
  ParseFieldSrcDlgPtr dlg;
  ParseFieldSrcPtr    pfsp;
  ValNodePtr          vnp;
  Int4                field_pos;

  dlg = (ParseFieldSrcDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }
  
  pfsp = (ParseFieldSrcPtr) MemNew (sizeof (ParseFieldSrcData));
  vnp = DialogToPointer (dlg->parse_field);
  if (vnp != NULL) {
    pfsp->parse_field = vnp->data.intvalue;
    vnp = ValNodeFree (vnp);
  }
  pfsp->comment_field = NULL;
  pfsp->sqdp = NULL;
  if (pfsp->parse_field == eParseFieldSrcStructuredComment) {
    field_pos = GetValue (dlg->comment_field);
    if (field_pos > 0) {
      vnp = dlg->comment_field_list;
      field_pos --;
      while (field_pos > 0 && vnp != NULL) {
        vnp = vnp->next;
        field_pos --;
      }
      if (vnp != NULL) {
        pfsp->comment_field = StringSave (vnp->data.ptrvalue);
      }
    }
  } else if (pfsp->parse_field == eParseFieldSrcSourceQual) {
    vnp = DialogToPointer (dlg->src_qual);
    if (vnp != NULL) {
      pfsp->sqdp = vnp->data.ptrvalue;
      vnp->data.ptrvalue = NULL;
      vnp = ValNodeFreeData (vnp);
    }
  }
  
  return pfsp;
}

static void AddCommentFieldName (SeqDescrPtr sdp, Pointer userdata)
{
  UserObjectPtr uop;
  ObjectIdPtr   oip;
  UserFieldPtr  ufp;
  ValNodePtr    vnp;
  Boolean       found;
  ValNodePtr PNTR comment_field_list = (ValNodePtr PNTR) userdata;
  
  if (comment_field_list == NULL 
      || sdp->choice != Seq_descr_user 
      || sdp->extended == 0
      || sdp->data.ptrvalue == NULL) {
    return;
  }

  uop = (UserObjectPtr) sdp->data.ptrvalue;
  oip = uop->type;
  if (oip != NULL && StringCmp (oip->str, "StructuredComment") == 0)
  {
    for (ufp = uop->data; ufp != NULL; ufp = ufp->next)
    {
      oip = ufp->label;
      if (oip != NULL && !StringHasNoText (oip->str))
      {
        found = FALSE;
        vnp = *comment_field_list;
        while (vnp != NULL && !found) {
          if (StringCmp (vnp->data.ptrvalue, oip->str) == 0) {
            found = TRUE;
          }
          vnp = vnp->next;
        }
        if (!found) {
          ValNodeAddPointer (comment_field_list, 0, StringSave (oip->str));
        }
      }
    }
  }
}


static ValNodePtr TestParseFieldSource (DialoG d)

{
  ParseFieldSrcDlgPtr dlg;
  ValNodePtr          head = NULL, vnp, vnp2;

  dlg = (ParseFieldSrcDlgPtr) GetObjectExtra (d);
  
  vnp = DialogToPointer (dlg->parse_field);
  if (vnp == NULL)
  {
    head = AddStringToValNodeChain (head, "field choice", 1);
  }
  else if (vnp->data.intvalue == eParseFieldSrcSourceQual)
  {
    vnp2 = DialogToPointer (dlg->src_qual);
    if (vnp2 == NULL)
    {
      head = AddStringToValNodeChain (head, "source qual", 1);
    }
    vnp2 = ValNodeFree (vnp2);
  }

  vnp = ValNodeFree (vnp);
  
  return head;  
}


extern DialoG ParseFieldSourceDialog
(GrouP                    h,
 SeqEntryPtr              sep,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  ParseFieldSrcDlgPtr dlg;
  GrouP               p, k;
  ValNodePtr          vnp;
  
  dlg = (ParseFieldSrcDlgPtr) MemNew (sizeof (ParseFieldSrcDlgData));
  if (dlg == NULL) 
  {
    return NULL;
  }
  
  p = HiddenGroup (h, 2, 0, NULL);

  SetObjectExtra (p, dlg, CleanupParseFieldSrcDialog);

  dlg->dialog = (DialoG) p;
  dlg->todialog = DataToParseFieldSrcDialog;
  dlg->fromdialog = ParseFieldSrcDialogToData;
  dlg->dialogmessage = ParseFieldSourceDialogMessage;
  dlg->testdialog = TestParseFieldSource;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->parse_field = FeatureFieldSelectionDialog (p, FALSE,
                                                  num_parse_src_fields, parse_src_field_list, 
                                                  ParseFieldChange, dlg);

  k = HiddenGroup (p, 0, 0, NULL);
  dlg->comment_field = PopupList (k, TRUE, NULL);
  dlg->src_qual = SourceQualTypeSelectionDialogEx (k, FALSE, change_notify,
                                                   change_userdata, 6, FALSE, FALSE);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->comment_field, (HANDLE) dlg->src_qual, NULL);
  
  /* Add to comment_field list all user fields found in SeqEntry sep */
  VisitDescriptorsInSep (sep, &(dlg->comment_field_list), AddCommentFieldName);
  for (vnp = dlg->comment_field_list; vnp != NULL; vnp = vnp->next) {
    PopupItem (dlg->comment_field, vnp->data.ptrvalue);
  }
    
  return (DialoG) p;  
}

typedef struct sampledialog
{
  DIALOG_MESSAGE_BLOCK
  DoC                doc;
  Int2               start_item;
  Int2               start_row;
  Int2               start_col;
  Int2               end_item;
  Int2               end_row;
  Int2               end_col;
  SetSamplePtr       ssp;
} SampleDialogData, PNTR SampleDialogPtr;

static void TruncateStringAfterTwoLines (CharPtr str, Uint4 char_per_line)
{
  CharPtr cp;
  Uint4    string_len;
  
  if (StringHasNoText (str))
  {
    return;
  }
  
  string_len = StringLen (str);
  if (string_len < char_per_line)
  {
    return;
  }
  
  cp = str + char_per_line;
  while (! isspace ((Int4)(*cp)) && cp > str)
  {
    cp--;
  }
  if (cp == str)
  {
    if (string_len > 2 * char_per_line)
    {
      str [2 * char_per_line] = 0;
    }
  }
  else
  {
    if (StringLen (cp) > char_per_line)
    {
      cp[char_per_line] = 0;
    }
  }
  
}

typedef struct exportsample
{
  ValNodePtr   row_list;
  SetSamplePtr ssp;
} ExportSampleData, PNTR ExportSamplePtr;

static void ExportSampleDialogCallback (BioseqPtr bsp, Pointer userdata)
{
  ExportSamplePtr   esp;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  Char              tmp[128];
  CharPtr           val_str = NULL, new_val_str, tmp_val_str;
  ValNodePtr        sample_field, sample_field_next;
  ValNodePtr        new_row = NULL;
  SeqIdPtr          sip;
  
  if (bsp == NULL || userdata == NULL || ! ISA_na (bsp->mol))
  {
    return;
  }
  
  esp = (ExportSamplePtr) userdata;
  
  sip = SeqIdFindBest (bsp->id, 0);
  SeqIdWrite (sip, tmp, PRINTID_REPORT, sizeof (tmp));
  
  ValNodeAddPointer (&new_row, 0, StringSave (tmp));
  
  for (sample_field = esp->ssp->field_list;
       sample_field != NULL; 
       sample_field = sample_field->next)
  {
    sample_field_next = sample_field->next;
    sample_field->next = NULL;
    val_str = NULL;
  
    if (esp->ssp->descrstring_func != NULL)
    {
      sdp = SeqMgrGetNextDescriptor (bsp, NULL, 0, &dcontext);
      while (sdp != NULL)
      {
        new_val_str = (esp->ssp->descrstring_func) (sdp, sample_field, NULL);
        tmp_val_str = StringAppend (val_str, new_val_str);
        val_str = MemFree (val_str);
        new_val_str = MemFree (new_val_str);
        val_str = tmp_val_str;
        tmp_val_str = NULL;
        sdp = SeqMgrGetNextDescriptor (bsp, sdp, 0, &dcontext);
      }
    }
    if (esp->ssp->fieldstring_func != NULL)
    {
      sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
      while (sfp != NULL)
      {
        new_val_str = (esp->ssp->fieldstring_func)(sfp, sample_field, NULL);
        tmp_val_str = StringAppend (val_str, new_val_str);
        val_str = MemFree (val_str);
        new_val_str = MemFree (new_val_str);
        val_str = tmp_val_str;
        tmp_val_str = NULL;
        sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
      }
    }
    ValNodeAddPointer (&new_row, 0, val_str);
    sample_field->next = sample_field_next;
  }
  ValNodeAddPointer (&(esp->row_list), 0, new_row);
}

static void GetTableDataColumnWidths (ValNodePtr row_list)
{
  ValNodePtr row_vnp, col_vnp, header_vnp;

  if (row_list == NULL || row_list->next == NULL)
  {
    return;
  }
  
  /* only collect widths for data columns, not header, so we can detect empty rows */
  for (row_vnp = row_list->next; row_vnp != NULL; row_vnp = row_vnp->next)
  {
    for (col_vnp = row_vnp->data.ptrvalue, header_vnp = row_list->data.ptrvalue;
         col_vnp != NULL && header_vnp != NULL;
         col_vnp = col_vnp->next, header_vnp = header_vnp->next)
    {
      if (!StringHasNoText (col_vnp->data.ptrvalue))
      {
        header_vnp->choice = MAX (header_vnp->choice, StringLen (col_vnp->data.ptrvalue));
      }
    }
  }
}

static void RemoveTableEmptyColumns (ValNodePtr row_list)
{
  ValNodePtr row_vnp, col_vnp, header_vnp, col_vnp_next, col_vnp_prev;
  ValNodePtr column_list;

  if (row_list == NULL)
  {
    return;
  }
  
  /* remove data columns */
  for (row_vnp = row_list->next; row_vnp != NULL; row_vnp = row_vnp->next)
  {
    column_list = row_vnp->data.ptrvalue;
    col_vnp_prev = NULL;
    for (col_vnp = column_list, header_vnp = row_list->data.ptrvalue;
         col_vnp != NULL && header_vnp != NULL;
         col_vnp = col_vnp_next, header_vnp = header_vnp->next)
    {
      col_vnp_next = col_vnp->next;
      if (header_vnp->choice == 0)
      {
        if (col_vnp_prev == NULL)
        {
          row_vnp->data.ptrvalue = col_vnp_next;
          column_list = row_vnp->data.ptrvalue;
        }
        else
        {
          col_vnp_prev->next = col_vnp_next;
        }
        col_vnp->next = NULL;
        ValNodeFreeData (col_vnp);
      }
      else
      {
        col_vnp_prev = col_vnp;
      }
    }
  }
  
  /* remove headers for empty columns */
  col_vnp_prev = NULL;
  for (col_vnp = row_list->data.ptrvalue;
       col_vnp != NULL;
       col_vnp = col_vnp_next)
  {
    col_vnp_next = col_vnp->next;
    if (col_vnp->choice == 0)
    {
      if (col_vnp_prev == NULL)
      {
        row_list->data.ptrvalue = col_vnp_next;
      }
      else
      {
        col_vnp_prev->next = col_vnp_next;
      }
      col_vnp->next = NULL;
      ValNodeFreeData (col_vnp);
    }
    else
    {
      col_vnp_prev = col_vnp;
    }
  }
}

static void ExportSampleDialogData (SampleDialogPtr dlg)
{
  Char             path [PATH_MAX];
  ExportSampleData esd;
  SeqEntryPtr      sep;
  CharPtr          sample_label;
  ValNodePtr       sample_field;
  FILE             *fp;
  ValNodePtr       new_row = NULL;
  
  if (dlg == NULL || dlg->ssp == NULL || dlg->ssp->label_vn_proc == NULL)
  {
    return;
  }
  
  sep = GetTopSeqEntryForEntityID (dlg->ssp->entityID);
  if (sep == NULL)
  {
    return;
  }
  if (! GetOutputFileName (path, sizeof (path), NULL)) return;
  fp = FileOpen (path, "w");
  if (fp == NULL)
  {
    Message (MSG_ERROR, "Unable to open %s", path);
    return;
  }
  WatchCursor ();
  Update ();
  esd.ssp = dlg->ssp;
  ValNodeAddPointer (&new_row, 0, StringSave ("Sequence"));
    
  for (sample_field = dlg->ssp->field_list;
       sample_field != NULL; 
       sample_field = sample_field->next)
  {
    sample_label = (dlg->ssp->label_vn_proc) (sample_field);
    ValNodeAddPointer (&new_row, 0, StringSave (sample_label));
  }
  esd.row_list = NULL;
  ValNodeAddPointer (&(esd.row_list), 0, new_row);
  
  VisitBioseqsInSep (sep, &esd, ExportSampleDialogCallback);

  GetTableDataColumnWidths (esd.row_list);
  RemoveTableEmptyColumns (esd.row_list);

  PrintTableDisplayRowListToFile (esd.row_list, fp);
  esd.row_list = FreeTableDisplayRowList (esd.row_list);
  
  FileClose (fp);
  ArrowCursor();
  Update();
}

static void RefreshSampleDialog (SampleDialogPtr dlg)
{
  SeqEntryPtr     sep;
  RecT            r;
  ParData         ParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
  ColData         ColFmt[] = 
  {
    {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE}
  };
  GetSampleData   gsd;
  CharPtr         new_line;
  ValNodePtr      sample_field, sample_field_next;
  CharPtr         sample_label = NULL;
  Int4            max_char_per_line;

  if (dlg == NULL)
  {
    return;
  }
  Reset (dlg->doc);
  if (dlg->ssp == NULL || dlg->ssp->field_list == NULL
      || (dlg->ssp->fieldstring_func == NULL && dlg->ssp->descrstring_func == NULL))
  {
    return;
  }
  
  sep = GetTopSeqEntryForEntityID (dlg->ssp->entityID);
  if (sep == NULL)
  {
    return;
  }
  gsd.sample_text = NULL;
  gsd.fieldstring_func = dlg->ssp->fieldstring_func;
  gsd.descrstring_func = dlg->ssp->descrstring_func;
    
  ObjectRect (dlg->doc, &r);
  InsetRect (&r, 4, 4);
  ColFmt[0].pixWidth = (r.right - r.left) / 3;
  ColFmt[1].pixWidth = (r.right - r.left) / 4;
  ColFmt[2].pixWidth = (r.right - r.left)  - ColFmt[0].pixWidth - ColFmt[1].pixWidth;

  SelectFont (programFont);
  max_char_per_line = ColFmt[2].pixWidth / CharWidth ('0');
  
  
  for (sample_field = dlg->ssp->field_list; sample_field != NULL; sample_field = sample_field->next)
  {
    sample_field_next = sample_field->next;
    sample_field->next = NULL;
    gsd.requested_field = sample_field;
    gsd.num_found = 0;
    gsd.all_same = TRUE;
    gsd.feat_dest_list = NULL;
    gsd.descr_dest_list = NULL;
    
    OperateOnSeqEntryConstrainedObjects (sep, dlg->ssp->fsp, GetSampleFeatureCallback, 
                                         GetSampleDescriptorCallback,
                                         0, 0, 0, &gsd);
    if (gsd.num_found > 0)
    {
      sample_label = (dlg->ssp->label_vn_proc) (sample_field);
      new_line = (CharPtr) MemNew (
           (StringLen (sample_label) 
            + StringLen (gsd.sample_text) + 20)
            * sizeof (Char));
      if (new_line != NULL)
      {
        if (gsd.all_same)
        {
          sprintf (new_line, "%s\t(%d same)\t", sample_label, gsd.num_found);
        }
        else
        {
          sprintf (new_line, "%s\t(%d mixed)\t", sample_label, gsd.num_found);
        }
        if (!StringHasNoText (gsd.sample_text))
        {
          TruncateStringAfterTwoLines (gsd.sample_text, max_char_per_line);
          StringCat (new_line, gsd.sample_text);
        }
        StringCat (new_line, "\n");
          AppendText (dlg->doc, new_line, &ParFmt, ColFmt, programFont);
          MemFree (new_line);
      }
      sample_label = MemFree (sample_label);
    }
    gsd.sample_text = MemFree (gsd.sample_text);
    sample_field->next = sample_field_next;
  }
  InvalDocument (dlg->doc);
}

static void PopulateSampleDialog (DialoG d, Pointer userdata)
{
  SampleDialogPtr dlg;
  SetSamplePtr    ssp;
  
  dlg = (SampleDialogPtr) GetObjectExtra (d);
  if (dlg == NULL || userdata == NULL)
  {
    return;
  }
  ssp = (SetSamplePtr) userdata;
  dlg->ssp = SetSampleFree (dlg->ssp);
  dlg->ssp = SetSampleCopy (ssp);
}

static void 
GetSampleDialogStartAndEnd 
(SampleDialogPtr dlg,
 Int2Ptr         start_item,
 Int2Ptr         start_row,
 Int2Ptr         start_col,
 Int2Ptr         end_item,
 Int2Ptr         end_row,
 Int2Ptr         end_col)
{
  if (dlg == NULL || start_item == NULL || start_row == NULL
      || start_col == NULL || end_item == NULL
      || end_row == NULL || end_col == NULL)
  {
    return;
  }
  
  if (dlg->start_item < 0)
  {
    *start_item = -1;
    *start_row = -1;
    *start_col = -1;
    *end_item = -1;
    *end_row = -1;
    *end_col = -1;
  }
  else if (dlg->end_item < 0)
  {
    *start_item = dlg->start_item;
    *start_row = dlg->start_row;
    *start_col = dlg->start_col;
    *end_item = *start_item;
    *end_row = *start_row;
    *end_col = *start_col;
  }
  else if (dlg->start_item <= dlg->end_item)
  {
    *start_item = dlg->start_item;
    *end_item = dlg->end_item;
    if (dlg->start_row <= dlg->end_row)
    {
      *start_row = dlg->start_row;
      *end_row = dlg->end_row;
      if (dlg->start_col <= dlg->end_col)
      {
        *start_col = dlg->start_col;
        *end_col = dlg->end_col;
      }
      else
      {
        *start_col = dlg->end_col;
        *end_col = dlg->start_col;
      }
    }
    else
    {
      *start_row = dlg->end_row;
      *start_col = dlg->end_col;
      *end_row = dlg->start_row;
      *end_col = dlg->start_col;
    }
  }
  else
  {
    *start_item = dlg->end_item;
    *start_row = dlg->end_row;
    *start_col = dlg->end_col;
    *end_item = dlg->start_item;
    *end_row = dlg->start_row;
    *end_col = dlg->start_col;
  }
  
}

static void SampleDialogCopy (SampleDialogPtr dlg)
{
  Int2       start_row = 0, end_row = 0, tmp_row, first_row, last_row;
  Int2       start_col = 0, end_col = 0, first_col, last_col;
  Int2       start_item = 0, end_item = 0, tmp_item;
  CharPtr    str;
  ValNodePtr strings = NULL;
  Int2       num_rows, num_cols;
  
  if (dlg->start_row < 0)
  {
    return;
  }
  
  GetSampleDialogStartAndEnd (dlg, &start_item, &start_row, &start_col,
                              &end_item, &end_row, &end_col);

  first_row = start_row;
  first_col = start_col;
  for (tmp_item = start_item; tmp_item <= end_item; tmp_item++)
  {
    GetItemParams (dlg->doc, tmp_item, NULL, &num_rows, &num_cols, NULL, NULL);
    if (tmp_item == end_item)
    {
      last_row = end_row;
    }
    else
    {
      last_row = num_rows;
    }
    for (tmp_row = first_row; tmp_row <= last_row; tmp_row++)
    {
      if (tmp_row == last_row && tmp_item == end_item)
      {
        last_col = end_col;
      }
      else
      {
        last_col = num_cols;
      }
      str = GetDocText (dlg->doc, tmp_item, tmp_row, 3);
      ValNodeAddPointer (&strings, 0, str);
    }
    first_row = 1;
  }
  str = MergeValNodeStrings (strings, FALSE);
  
  StringToClipboard (str);
  MemFree (str);
  strings = ValNodeFreeData (strings);
}

static void InvalidateSampleDialogRows (DoC d, Int2 start_item, Int2 end_item)
{
  Int2 num_rows;
  if (d == NULL)
  {
    return;
  }
  if (start_item < 1)
  {
    start_item = 1;
  }
  if (end_item < 1)
  {
    end_item = start_item;
  }
  while (start_item <= end_item)
  {
    GetItemParams (d, start_item, NULL, &num_rows, NULL, NULL, NULL);
    InvalDocRows (d, start_item, 1, num_rows); 
    start_item++;       
  }
}

static void SampleDialogOnClick (DoC d, PoinT pt)
{
  SampleDialogPtr   dlg;
  Int2              pos_item, pos_row, pos_col;
  Int2              start_item = 0, start_row, start_col;
  Int2              end_item = -1, end_row, end_col;
  
  dlg = (SampleDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  
  MapDocPoint (d, pt, &(pos_item), &(pos_row), &(pos_col), NULL);
  if (dlg->start_item == pos_item
      && dlg->start_row == pos_row 
      && dlg->start_col == pos_col)
  {
    dlg->start_row = -1;
    dlg->start_col = -1;
    
  }
  else
  {
    if (dlg->start_item > -1)
    {
      GetSampleDialogStartAndEnd (dlg, &start_item, &start_row, &start_col,
                                  &end_item, &end_row, &end_col);
    }
    dlg->start_item = pos_item;
    dlg->start_row = pos_row;
    dlg->start_col = pos_col;
  }
  dlg->end_item = -1;
  dlg->end_row = -1;
  dlg->end_col = -1;
  InvalidateSampleDialogRows (d, start_item, end_item);
  InvalidateSampleDialogRows (d, pos_item, pos_item);
}

static void SampleDialogOnDrag (DoC d, PoinT pt)
{
  SampleDialogPtr   dlg;
  
  dlg = (SampleDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  
  MapDocPoint (d, pt, &(dlg->end_item), &(dlg->end_row), &(dlg->end_col), NULL);
  InvalDocument (d);
}

static void SampleDialogOnRelease (DoC d, PoinT pt)
{
  SampleDialogPtr   dlg;
  
  dlg = (SampleDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  
  MapDocPoint (d, pt, &(dlg->end_item), &(dlg->end_row), &(dlg->end_col), NULL);
  InvalDocument (d);
}

static Boolean SampleDialogInvert (DoC d, Int2 item, Int2 row, Int2 col)
{
  SampleDialogPtr   dlg;
  Int2              start_item = 0, start_row = 0, start_col = 0;
  Int2              end_item = 0, end_row = 0, end_col = 0;
  Boolean           rval = FALSE;
  
  dlg = (SampleDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return FALSE;
  
  if (dlg->start_row < 0)
  {
    return FALSE;
  }
  
  if (dlg->end_item == -1)
  {
    if(dlg->start_item == item
       && dlg->start_row == row
       && dlg->start_col == col)
    {
      return TRUE;
    }
    else
    {
      return FALSE;
    }
  }
  GetSampleDialogStartAndEnd (dlg, &start_item, &start_row, &start_col,
                              &end_item, &end_row, &end_col);
  if (start_item <= item && end_item >= item)
  {
    rval = TRUE;
    if (start_item == item)
    {
      if (start_row == row)
      {
        if (start_col > col)
        {
          rval = FALSE;
        }
      }
      else if (start_row > row)
      {
        rval = FALSE;
      }
    }
    
    if (end_item == item)
    {
      if (end_row == row)
      {
        if (end_col < col)
        {
          rval = FALSE;
        }
      }
      else if (end_row < row)
      {
        rval = FALSE;
      }
    }
  }
  if (col == 3)
  {
    return rval;
  }
  else
  {
    return FALSE;
  }
}

static void SampleDialogOnKey (SlatE s, Char ch)
{
  DoC             d;
  SampleDialogPtr dlg;

  d = (DoC) s;
  Select (d);
  dlg = (SampleDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  CaptureSlateFocus ((SlatE) s);
  
#ifdef WIN_MSWIN
  if (ch == 3)
  {
    SampleDialogCopy (dlg);
  }
#else
  if (ctrlKey)
  {
    if (ch == 'c')
    {
      SampleDialogCopy (dlg);
    }
  }
#endif
}

static void SampleDialogMessage (DialoG d, Int2 mssg)

{
  SampleDialogPtr dlg;

  dlg = (SampleDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        /* reset list */
        RefreshSampleDialog (dlg);
        break;
      default :
        break;
    }
  }
}

static void CleanupSampleDialog (GraphiC g, VoidPtr data)
{
  SampleDialogPtr dlg;

  dlg = (SampleDialogPtr) data;
  if (dlg != NULL) {
    dlg->ssp = SetSampleFree (dlg->ssp);
  }
  StdCleanupExtraProc (g, data);
}

static void RefreshSampleBtn (ButtoN b)
{
  SampleDialogPtr dlg;
  
  dlg = (SampleDialogPtr) GetObjectExtra (b);
  RefreshSampleDialog (dlg);  
}

static void ExportSampleBtn (ButtoN b)
{
  SampleDialogPtr dlg;
  
  dlg = (SampleDialogPtr) GetObjectExtra (b);
  ExportSampleDialogData (dlg);  
}

extern DialoG SampleDialog (GrouP h)
{
  SampleDialogPtr dlg;
  GrouP           p, g;
  ButtoN          b;
  
  dlg = (SampleDialogPtr) MemNew (sizeof (SampleDialogData));
  if (dlg == NULL) 
  {
    return NULL;
  }
  
  p = HiddenGroup (h, -1, 0, NULL);

  SetObjectExtra (p, dlg, CleanupSampleDialog);

  dlg->dialog = (DialoG) p;
  dlg->todialog = PopulateSampleDialog;
  dlg->fromdialog = NULL;
  dlg->dialogmessage = SampleDialogMessage;
  dlg->testdialog = NULL;

  dlg->doc = DocumentPanel (p, stdCharWidth * 27, stdLineHeight * 8);
  SetObjectExtra (dlg->doc, dlg, NULL);
  SetDocProcs (dlg->doc, SampleDialogOnClick, SampleDialogOnDrag, SampleDialogOnRelease, NULL);  
  SetDocShade (dlg->doc, NULL, NULL, SampleDialogInvert, NULL);
  SetDocAutoAdjust (dlg->doc, TRUE);
  SetSlateChar ((SlatE) dlg->doc, SampleDialogOnKey);
  
  g = HiddenGroup (p, 2, 0, NULL);
  b = PushButton (g, "Refresh Qualifiers", RefreshSampleBtn);
  SetObjectExtra (b, dlg, NULL);
  b = PushButton (g, "Export Qualifiers", ExportSampleBtn);
  SetObjectExtra (b, dlg, NULL);  
  
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->doc, (HANDLE) g, NULL);
  
  return (DialoG) p;  
}

typedef struct featurefieldchoicedialog 
{
  DIALOG_MESSAGE_BLOCK
  DialoG first_choice;
  DialoG second_choice;
  DialoG third_choice;
  ButtoN remove_first;
  ButtoN remove_second;
  ButtoN remove_third;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;      
} FeatureFieldChoiceDialogData, PNTR FeatureFieldChoiceDialogPtr;

static void FeatureFieldChoiceChange (Pointer data)
{
  FeatureFieldChoiceDialogPtr dlg;
  ValNodePtr                  vnp;

  dlg = (FeatureFieldChoiceDialogPtr) data;
  if (dlg == NULL)
  {
    return;
  }
  
  Disable (dlg->remove_first);
  Disable (dlg->remove_second);
  Disable (dlg->remove_third);
  
  Disable (dlg->second_choice);
  Disable (dlg->third_choice);
  
  vnp = DialogToPointer (dlg->first_choice);
  if (vnp != NULL)
  {
    ValNodeFree (vnp);
    Enable (dlg->second_choice);
    Enable (dlg->remove_first);
    vnp = DialogToPointer (dlg->second_choice);
    if (vnp != NULL)
    {
      ValNodeFree (vnp);
      Enable (dlg->third_choice);
      Enable (dlg->remove_second);
      vnp = DialogToPointer (dlg->third_choice);
      if (vnp != NULL)
      {
        ValNodeFree (vnp);
        Enable (dlg->remove_third);
      }
    }
  }
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify)(dlg->change_userdata);
  }
}

static void ResetFeatureFieldChoiceDialog (FeatureFieldChoiceDialogPtr dlg)
{
  if (dlg == NULL) 
  {
    return;
  }
  PointerToDialog (dlg->first_choice, NULL);
  PointerToDialog (dlg->second_choice, NULL);
  PointerToDialog (dlg->third_choice, NULL);
  if (dlg->remove_first != NULL)
  {
    SetStatus (dlg->remove_first, FALSE);
  }
  if (dlg->remove_second != NULL)
  {
    SetStatus (dlg->remove_second, FALSE);
  }
  if (dlg->remove_third != NULL)
  {
    SetStatus (dlg->remove_third, FALSE);
  }
  FeatureFieldChoiceChange (dlg);
}

static void FeatureFieldChoiceToDialog (DialoG d, Pointer data)
{
  FeatureFieldChoiceDialogPtr dlg;
  ValNodePtr                  feature_field;
  
  dlg = (FeatureFieldChoiceDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  ResetFeatureFieldChoiceDialog (dlg);
  feature_field = (ValNodePtr) data;
  if (feature_field != NULL)
  {
    PointerToDialog (dlg->first_choice, feature_field);
    if (feature_field->choice != 0 && dlg->remove_first != NULL)
    {
      SetStatus (dlg->remove_first, TRUE);
    }
    feature_field = feature_field->next;
    if (feature_field != NULL)
    {
      PointerToDialog (dlg->second_choice, feature_field);
      if (feature_field->choice != 0 && dlg->remove_second != NULL)
      {
        SetStatus (dlg->remove_second, TRUE);
      }
      feature_field = feature_field->next;
      if (feature_field != NULL)
      {
        if (feature_field->choice != 0 && dlg->remove_third != NULL)
        {
          SetStatus (dlg->remove_third, TRUE);
        }
        PointerToDialog (dlg->third_choice, feature_field);
      }
    }
  }
  FeatureFieldChoiceChange (dlg);
}

static Pointer DialogToFeatureFieldChoice (DialoG d)
{
  FeatureFieldChoiceDialogPtr dlg;
  ValNodePtr                  choice_list = NULL;
  ValNodePtr                  feature_field, vnp_prev;

  dlg = (FeatureFieldChoiceDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  feature_field = DialogToPointer (dlg->first_choice);
  if (feature_field != NULL)
  {
    if (dlg->remove_first != NULL && GetStatus (dlg->remove_first))
    {
      feature_field->choice = 1;
    }
    choice_list = feature_field;
    
    feature_field = DialogToPointer (dlg->second_choice);
    if (feature_field != NULL)
    {
      if (dlg->remove_second != NULL && GetStatus (dlg->remove_second))
      {
        feature_field->choice = 1;
      }
      choice_list->next = feature_field;
      vnp_prev = feature_field;
      
      feature_field = DialogToPointer (dlg->third_choice);
      if (feature_field != NULL)
      {
        if (dlg->remove_third != NULL && GetStatus (dlg->remove_third))
        {
          feature_field->choice = 1;
        }
        vnp_prev->next = feature_field;
      }
      else
      {
        ValNodeFree (feature_field);
      }
    }
    else
    {
      ValNodeFree (feature_field);
    }
  }
  else
  {
    ValNodeFree (feature_field);
  }
  return choice_list;
}

static void FeatureFieldChoiceMessage (DialoG d, Int2 mssg)

{
  FeatureFieldChoiceDialogPtr  dlg;

  dlg = (FeatureFieldChoiceDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        /* reset list */
        ResetFeatureFieldChoiceDialog (dlg);
        break;
      case VIB_MSG_ENTER :
        Select (dlg->first_choice);
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestFeatureFieldChoice (DialoG d)

{
  FeatureFieldChoiceDialogPtr dlg;
  ValNodePtr                  head = NULL;
  ValNodePtr                  vnp;

  dlg = (FeatureFieldChoiceDialogPtr) GetObjectExtra (d);
  vnp = DialogToPointer (d);
  if (vnp == NULL)
  {
    head = AddStringToValNodeChain (head, "feature field choice", 1);
  }
  else
  {
    ValNodeFree (vnp);
  }
  return head;  
}

extern DialoG FeatureFieldChoiceDialog 
(GrouP h,
 FeatureFieldSelectionProc make_fieldlist_dlg,
 Boolean                   offer_to_remove,
 Nlm_ChangeNotifyProc      change_notify,
 Pointer                   change_userdata)
{
  FeatureFieldChoiceDialogPtr  dlg;
  GrouP                     p, g;

  if (make_fieldlist_dlg == NULL)
  {
    return NULL;
  }
  
  dlg = (FeatureFieldChoiceDialogPtr) MemNew (sizeof (FeatureFieldChoiceDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, 0, 6, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = FeatureFieldChoiceToDialog;
  dlg->fromdialog = DialogToFeatureFieldChoice;
  dlg->dialogmessage = FeatureFieldChoiceMessage;
  dlg->testdialog = TestFeatureFieldChoice;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  g = HiddenGroup (p, 3, 0, NULL);
  StaticPrompt (g, "1st Choice", 0, 0, programFont, 'c');
  dlg->first_choice = (make_fieldlist_dlg) (g, TRUE, FeatureFieldChoiceChange, dlg);
  if (offer_to_remove)
  {
    dlg->remove_first = CheckBox (g, "Remove if used", NULL);
    Disable (dlg->remove_first);
  }
  else
  {
    dlg->remove_first = NULL;
  }
  
  g = HiddenGroup (p, 3, 0, NULL);
  StaticPrompt (g, "2nd Choice", 0, 0, programFont, 'c');
  dlg->second_choice = (make_fieldlist_dlg) (g, TRUE, FeatureFieldChoiceChange, dlg);
  if (offer_to_remove)
  {
    dlg->remove_second = CheckBox (g, "Remove if used", NULL);
    Disable (dlg->remove_second);
  }
  else
  {
    dlg->remove_second = NULL;
  }
  Disable (dlg->second_choice);
  
  g = HiddenGroup (p, 3, 0, NULL);
  StaticPrompt (g, "3rd Choice", 0, 0, programFont, 'c');
  dlg->third_choice = (make_fieldlist_dlg) (g, TRUE, FeatureFieldChoiceChange, dlg);
  if (offer_to_remove)
  {
    dlg->remove_third = CheckBox (g, "Remove if used", NULL);
    Disable (dlg->remove_third);
  }
  else
  {
    dlg->remove_third = NULL;
  }
  Disable (dlg->third_choice);

  return (DialoG) p;
}


typedef ValNodePtr (*AdjustConvertFeatureListFunc) PROTO ((ValNodePtr, Pointer, BaseFormPtr, BoolPtr cancel));
typedef Pointer (*CollectConvertFeatureOptionsFunc) PROTO ((ValNodePtr, BaseFormPtr, BoolPtr cancel));
typedef Pointer (*FreeConvertFeatureOptionsFunc) PROTO ((Pointer));
typedef Boolean (*ConvertFeatureTypeFunc) PROTO ((SeqFeatPtr, Uint2, Pointer));
typedef void (*ConversionCleanupFunc) PROTO ((Uint2, Pointer));

typedef struct convertfeatureprocs {
  Uint2 seqfeat_from;
  Uint2 featdef_from;
  Uint2 seqfeat_to;
  Uint2 featdef_to;
  AdjustConvertFeatureListFunc adjust_features;
  CollectConvertFeatureOptionsFunc collect_options;
  FreeConvertFeatureOptionsFunc    free_options;
  ConvertFeatureTypeFunc           convert_func;
  ConversionCleanupFunc            cleanup_func;
  CharPtr                          help_text;
} ConvertFeatureProcsData, PNTR ConvertFeatureProcsPtr;

static Pointer SimpleOptsFree (Pointer);
static Boolean ConvertCDSToRNAEx (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertGeneToRNAEx (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static ValNodePtr AdjustConvertFeaturesForBioSrcToRepeatRegion (ValNodePtr feat_list, Pointer extradata, BaseFormPtr bfp, BoolPtr cancel);
static Boolean ConvertBioSrcToRepeatRegionEx (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static ValNodePtr AdjustCDSToMiscFeatList (ValNodePtr feat_list, Pointer extradata, BaseFormPtr bfp, BoolPtr cancel);
static Pointer GetCDSConversionOpts (ValNodePtr feat_list, BaseFormPtr bfp, BoolPtr cancel);
static Boolean CDSToMiscFeatConvertFunc (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static void CDSToMiscFeatConversionCleanup (Uint2 entityID, Pointer extradata);
static Boolean ConvertImpToProt (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertProtToImp (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertImpToSpecialRNA (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertRegionToImp (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertRegionToRNA (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertImpToRNA (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertCommentToMiscFeat (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertGeneToMiscFeat (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertRNAToMiscFeat (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertSiteToMiscFeat (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertProtToRegion (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertRegionToProt (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Pointer GetBondOpts (ValNodePtr feat_list, BaseFormPtr bfp, BoolPtr cancel);
static Pointer GetSiteOpts (ValNodePtr feat_list, BaseFormPtr bfp, BoolPtr cancel);
static Boolean ConvertToSiteOrBond (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Pointer GetRegionOpts (ValNodePtr feat_list, BaseFormPtr bfp, BoolPtr cancel);
static ValNodePtr AdjustRegionConversionFeatures (ValNodePtr feat_list, Pointer extradata, BaseFormPtr bfp, BoolPtr cancel);
static Boolean ConvertAnyToRegion (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertImpToImp (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertRNAToRNA (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertProtToProt (SeqFeatPtr, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertCDSToMatPeptide (SeqFeatPtr sfp, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertMiscFeatureToCDSFunction (SeqFeatPtr sfp, Uint2 featdef_to, Pointer extradata);
static Boolean ConvertMiscFeatureToGeneFunction (SeqFeatPtr sfp, Uint2 featdef_to, Pointer extradata);

static ConvertFeatureProcsData ConvertFeaturesTable[] = {
  { SEQFEAT_CDREGION, FEATDEF_CDS,                SEQFEAT_RNA,    FEATDEF_ANY,
        NULL,  GetCDSConversionOpts, SimpleOptsFree, ConvertCDSToRNAEx, NULL, 
        "Delete protein product sequence.\nClear product field if transcript ID removal was requested.\nIf converting to tRNA and anticodon value can be parsed from label, set aa value, and add any text that could not be parsed into an anticodon value to the feature note.\nIf converting to other RNA, put label in RNA product." },
  { SEQFEAT_GENE,     FEATDEF_GENE,               SEQFEAT_RNA,    FEATDEF_ANY,
        NULL,  NULL,                   NULL,           ConvertGeneToRNAEx, NULL,
        "If converting to tRNA and anticodon value can be parsed from label, set aa value, and add any text that could not be parsed into an anticodon value to the feature note.  If converting to other RNA, put label in RNA product.  Also append gene locus, allele, description, map location, and locus tag to comment (as long as these values are not already in the label and therefore in the RNA product)." },
  { SEQFEAT_BIOSRC,   FEATDEF_BIOSRC,             SEQFEAT_IMP,    FEATDEF_repeat_region, 
       AdjustConvertFeaturesForBioSrcToRepeatRegion, NULL, NULL, ConvertBioSrcToRepeatRegionEx, NULL,
        "Creates a repeat_region with mobile_element qualifiers for the transposon and/or insertion sequence qualifiers on the BioSource.  All other BioSource information is discarded." },
  { SEQFEAT_CDREGION, FEATDEF_CDS,                SEQFEAT_IMP,    FEATDEF_misc_feature, 
       AdjustCDSToMiscFeatList, GetCDSConversionOpts, SimpleOptsFree, CDSToMiscFeatConvertFunc, CDSToMiscFeatConversionCleanup,
       "Copy comment from coding region to new misc_feature and remove product field.  If not pseudo coding region, add product name from protein feature to new misc_feature comment and delete product sequence." },
  { SEQFEAT_IMP,      FEATDEF_ANY,                SEQFEAT_PROT,   FEATDEF_ANY, 
       NULL, NULL, NULL, ConvertImpToProt, NULL,
       "Original feature must be on nucleotide sequence and be contained in coding region location.  Coding region must have product protein sequence.  New feature is created on product protein sequence so that the translated location will be as close as possible to the original nucleotide location (may not be exact because of codon boundaries)." },
  { SEQFEAT_PROT,     FEATDEF_mat_peptide_aa,     SEQFEAT_IMP,    FEATDEF_ANY,
       NULL, NULL, NULL, ConvertProtToImp, NULL,
       "Original feature must be on a protein sequence that is a product of a coding region.\nNew feature will be created on same sequence as coding region.\n"
       "If protein feature has name, this will be saved as /product qualifier on new feature.\nIf protein feature does not have name but does have description, this will be saved as /product qualifier on new feature.\n"
       "EC_number values from the protein feature will be saved as /EC_number qualifiers on the new feature.\nActivity values will be saved as /function qualifiers on the new feature.\n"
       "Db_xref values from the protein feature will be saved as /db_xref qualifers on the new feature." },
  { SEQFEAT_PROT,     FEATDEF_sig_peptide_aa,     SEQFEAT_IMP,    FEATDEF_ANY,
       NULL, NULL, NULL, ConvertProtToImp, NULL,
       "Original feature must be on a protein sequence that is a product of a coding region.\nNew feature will be created on same sequence as coding region.\n"
       "If protein feature has name, this will be saved as /product qualifier on new feature.\nIf protein feature does not have name but does have description, this will be saved as /product qualifier on new feature.\n"
       "EC_number values from the protein feature will be saved as /EC_number qualifiers on the new feature.\nActivity values will be saved as /function qualifiers on the new feature.\n"
       "Db_xref values from the protein feature will be saved as /db_xref qualifers on the new feature." },
  { SEQFEAT_PROT,     FEATDEF_transit_peptide_aa, SEQFEAT_IMP,    FEATDEF_ANY,
       NULL, NULL, NULL, ConvertProtToImp, NULL,
       "Original feature must be on a protein sequence that is a product of a coding region.\nNew feature will be created on same sequence as coding region.\n"
       "If protein feature has name, this will be saved as /product qualifier on new feature.\nIf protein feature does not have name but does have description, this will be saved as /product qualifier on new feature.\n"
       "EC_number values from the protein feature will be saved as /EC_number qualifiers on the new feature.\nActivity values will be saved as /function qualifiers on the new feature.\n"
       "Db_xref values from the protein feature will be saved as /db_xref qualifers on the new feature." },
  { SEQFEAT_IMP,      FEATDEF_misc_feature,                SEQFEAT_CDREGION,    FEATDEF_CDS,
       NULL, NULL, NULL, ConvertMiscFeatureToCDSFunction, NULL,
       "Use misc_feature comment for coding region product name." },
  { SEQFEAT_IMP,      FEATDEF_misc_feature,                SEQFEAT_GENE,    FEATDEF_GENE,
       NULL, NULL, NULL, ConvertMiscFeatureToGeneFunction, NULL,
       "Use misc_feature comment for gene locus." },
  { SEQFEAT_IMP,      FEATDEF_ANY,                SEQFEAT_RNA,    FEATDEF_misc_RNA,
       NULL, NULL, NULL, ConvertImpToSpecialRNA, NULL,
       "Creates a misc_RNA.  Import feature key is discarded." },
  { SEQFEAT_IMP,      FEATDEF_ANY,                SEQFEAT_RNA,    FEATDEF_precursor_RNA,
       NULL, NULL, NULL, ConvertImpToSpecialRNA, NULL,
       "Creates a precursor RNA.  Import feature key is discarded." },
  { SEQFEAT_REGION,   FEATDEF_REGION,             SEQFEAT_IMP,    FEATDEF_ANY,
       NULL, NULL, NULL, ConvertRegionToImp, NULL,
       "Creates a misc_feature with the region name saved as a /note qualifier." },
  { SEQFEAT_REGION,   FEATDEF_REGION,             SEQFEAT_RNA,    FEATDEF_misc_RNA,
       NULL, NULL, NULL, ConvertRegionToRNA, NULL,
       "Creates a misc_RNA feature with the region name as the product name." },
  { SEQFEAT_REGION,   FEATDEF_REGION,             SEQFEAT_RNA,    FEATDEF_precursor_RNA,
       NULL, NULL, NULL, ConvertRegionToRNA, NULL,
       "Creates a precursor RNA feature with the region name as the product name." },
  { SEQFEAT_IMP,      FEATDEF_ANY,                SEQFEAT_RNA,    FEATDEF_ANY,
       NULL, NULL, NULL, ConvertImpToRNA, NULL,
       "Creates an RNA feature of the specified type." },
  { SEQFEAT_COMMENT,  FEATDEF_ANY,                SEQFEAT_IMP,    FEATDEF_misc_feature,
       NULL, NULL, NULL, ConvertCommentToMiscFeat, NULL,
       "Creates a misc_feature with the same note as the original.  Note the flatfile display for the feature is the same." },
  { SEQFEAT_GENE,     FEATDEF_GENE,               SEQFEAT_IMP,    FEATDEF_misc_feature,
       NULL, NULL, NULL, ConvertGeneToMiscFeat, NULL,
       "Creates a misc_feature with the gene description and locus prepended to the original comment, separated by semicolons." },
  { SEQFEAT_RNA,      FEATDEF_ANY,                SEQFEAT_IMP,    FEATDEF_misc_feature,
       NULL, NULL, NULL, ConvertRNAToMiscFeat, NULL,
       "Creates a misc_feature and adds the RNA product name to the comment." } ,
  { SEQFEAT_SITE,     FEATDEF_ANY,                SEQFEAT_IMP,    FEATDEF_misc_feature,
       NULL, NULL, NULL, ConvertSiteToMiscFeat, NULL,
       "Creates a misc_feature with the site type name as a /note qualifier." } ,
  { SEQFEAT_PROT,     FEATDEF_mat_peptide_aa,     SEQFEAT_REGION, FEATDEF_REGION,
       NULL, NULL, NULL, ConvertProtToRegion, NULL,
       "Creates a Region feature with the protein name as the region name." },
  { SEQFEAT_PROT,     FEATDEF_sig_peptide_aa,     SEQFEAT_REGION, FEATDEF_REGION,
       NULL, NULL, NULL, ConvertProtToRegion, NULL,
       "Creates a Region feature with the protein name as the region name." },
  { SEQFEAT_PROT,     FEATDEF_transit_peptide_aa, SEQFEAT_REGION, FEATDEF_REGION,
       NULL, NULL, NULL, ConvertProtToRegion, NULL,
       "Creates a Region feature with the protein name as the region name." },
  { SEQFEAT_REGION,   FEATDEF_REGION,             SEQFEAT_PROT,   FEATDEF_mat_peptide_aa,
       NULL, NULL, NULL, ConvertRegionToProt, NULL,
       "If feature is on nucleotide sequence, will create feature on protein product sequence for overlapping coding region.  Protein name will be region name." },
  { SEQFEAT_REGION,   FEATDEF_REGION,             SEQFEAT_PROT,   FEATDEF_sig_peptide_aa,
       NULL, NULL, NULL, ConvertRegionToProt, NULL,
       "If feature is on nucleotide sequence, will create feature on protein product sequence for overlapping coding region.  Protein name will be region name." },
  { SEQFEAT_REGION,   FEATDEF_REGION,             SEQFEAT_PROT,   FEATDEF_transit_peptide_aa,
       NULL, NULL, NULL, ConvertRegionToProt, NULL,
       "If feature is on nucleotide sequence, will create feature on protein product sequence for overlapping coding region.  Protein name will be region name." },
  { SEQFEAT_REGION,   FEATDEF_REGION,             SEQFEAT_PROT,   FEATDEF_preprotein,
       NULL, NULL, NULL, ConvertRegionToProt, NULL,
       "If feature is on nucleotide sequence, will create feature on protein product sequence for overlapping coding region.  Protein name will be region name." },
  { 0,                FEATDEF_ANY,                SEQFEAT_BOND,    FEATDEF_BOND,
       NULL, GetBondOpts, SimpleOptsFree, ConvertToSiteOrBond, NULL,
       "Create Bond feature with specified bond type.  Location is a SeqLocBond with a point at the start of the original location and a point at the end of the original location.  All feature ID, partialness, except, comment, product, location, genbank qualifiers, title, citation, experimental evidence, gene xrefs, db xrefs, and pseudo-ness information is discarded." },
  { 0,                FEATDEF_ANY,                SEQFEAT_SITE,    FEATDEF_SITE,
       NULL, GetSiteOpts, SimpleOptsFree, ConvertToSiteOrBond, NULL,
       "Create Site feature with specified site type.  All feature ID, partialness, except, comment, product, location, genbank qualifiers, title, citation, experimental evidence, gene xrefs, db xrefs, and pseudo-ness information is discarded." },
  { 0,                FEATDEF_ANY,                SEQFEAT_REGION,    FEATDEF_REGION,
       AdjustRegionConversionFeatures, GetRegionOpts, SimpleOptsFree, ConvertAnyToRegion, NULL,
       "Create Region feature on nucleotide sequence or protein product sequence of overlapping coding region as specified.  Use comment on feature for region name.\n"
       "All feature ID, partialness, except, comment, product, location, genbank qualifiers, title, citation, experimental evidence, gene xrefs, db xrefs, and pseudo-ness information is discarded." },
  { SEQFEAT_IMP,      FEATDEF_ANY,                SEQFEAT_IMP,    FEATDEF_ANY,
       NULL, NULL, NULL, ConvertImpToImp, NULL,
       "Changes type of import feature." },
  { SEQFEAT_RNA,      FEATDEF_ANY,                SEQFEAT_RNA,    FEATDEF_ANY,
       NULL, NULL, NULL, ConvertRNAToRNA, NULL,
       "Changes type of RNA feature." },
  { SEQFEAT_RNA,      FEATDEF_ncRNA,              SEQFEAT_IMP,    FEATDEF_misc_binding,
       NULL, NULL, NULL, ConvertncRNAToMiscBinding, NULL,
       "Changes ncRNA to misc_binding." },
  { SEQFEAT_PROT,     FEATDEF_ANY,                SEQFEAT_PROT,   FEATDEF_ANY,
       NULL, NULL, NULL, ConvertProtToProt, NULL,
       "Changes type of protein feature." },
  { SEQFEAT_CDREGION, FEATDEF_CDS,                SEQFEAT_PROT,   FEATDEF_mat_peptide_aa,
       NULL, NULL, NULL, ConvertCDSToMatPeptide, NULL,
       "If coding region is overlapped by another coding region, will convert the coding region to a mat-peptide on the overlapping coding region's protein sequence, otherwise if you have checked \"Leave Original Feature\" it will create a mat-peptide with the same protein names and description on the protein sequence for the coding region." }
};

static Int4 num_convert_feature_table_lines = sizeof (ConvertFeaturesTable) / sizeof (ConvertFeatureProcsData);

static Int4 GetConvertFeatureTableLine (Uint2 seqfeat_from, Uint2 featdef_from, Uint2 seqfeat_to, Uint2 featdef_to)
{
  Int4 i, table_line_num = -1;

  for (i = 0; i < num_convert_feature_table_lines && table_line_num == -1; i++)
  {
    if ((ConvertFeaturesTable[i].seqfeat_from == 0 || ConvertFeaturesTable[i].seqfeat_from == seqfeat_from)
        && (ConvertFeaturesTable[i].featdef_from == FEATDEF_ANY || ConvertFeaturesTable[i].featdef_from == featdef_from)
        && (ConvertFeaturesTable[i].seqfeat_to == 0 || ConvertFeaturesTable[i].seqfeat_to == seqfeat_to)
        && (ConvertFeaturesTable[i].featdef_to == FEATDEF_ANY || ConvertFeaturesTable[i].featdef_to == featdef_to))
    {
      table_line_num = i;
    }
  }
  return table_line_num;
}


typedef struct featureselremconform
{
  FORM_MESSAGE_BLOCK
  DialoG  feature_select;
  DialoG  constraints;
  DialoG  accept_cancel;
  PopuP   action_choice;
  DialoG  feature_select_to;
  PrompT  from_prompt;
  PrompT  to_prompt;
  DialoG  feature_select_from;
  GrouP   remove_grp;
  GrouP   convert_grp;
  ButtoN  clear_constraints_on_action_change;
  ButtoN  leave_original_feature_btn;  /* checkbox for conversion only */

  DoC     help_text; /* for conversion only - describes conversion */  

  /* used for convert callback */
  Uint1      featdef_to;
  Boolean    renormalize_nucprot;
  Boolean    leave_original_feature;
  ValNodePtr feat_list;
  BaseFormPtr bfp;
} FeatureSelRemConvFormData, PNTR FeatureSelRemConvFormPtr;

static void CleanupFeatureSelRemConvForm (GraphiC g, VoidPtr data)

{
  FeatureSelRemConvFormPtr  mrfp;

  mrfp = (FeatureSelRemConvFormPtr) data;
  if (mrfp != NULL) {
  }
  StdCleanupFormProc (g, data);
}

static void FeatureRemoveClearText (Pointer data)
{
  FeatureSelRemConvFormPtr   mrfp;
  FilterSetPtr          fsp;

  mrfp = (FeatureSelRemConvFormPtr) data;
  if (mrfp == NULL) return;
 
  fsp = DialogToPointer (mrfp->constraints);
  FilterSetClearText (fsp);
  PointerToDialog (mrfp->constraints, fsp);
  FilterSetFree (fsp);
}

static void FeatureRemoveClear (Pointer data)
{
  FeatureSelRemConvFormPtr mrfp;

  mrfp = (FeatureSelRemConvFormPtr) data;
  if (mrfp == NULL) return;
 
  PointerToDialog (mrfp->feature_select, NULL);
  PointerToDialog (mrfp->constraints, NULL);
  PointerToDialog (mrfp->feature_select_to, NULL);
  PointerToDialog (mrfp->feature_select_from, NULL);
  SetStatus (mrfp->leave_original_feature_btn, FALSE);
}

#define FEATURE_REMOVE   1
#define FEATURE_CONVERT  2
#define FEATURE_SELECT   3
#define FEATURE_DESELECT 4

static void FeatureRemoveChangeNotify (Pointer userdata)
{
  FeatureSelRemConvFormPtr mrfp;
  ValNodePtr          err_list = NULL, vnp_from, vnp_to;
  Int4                table_line;
  RecT                r;
  Boolean             need_original = FALSE;

  mrfp = (FeatureSelRemConvFormPtr) userdata;
  if (mrfp == NULL) return;

  if (mrfp->action_choice == NULL)
  {
    /* just the remove dialog */
    err_list = TestDialog (mrfp->feature_select);
  }
  else if (GetValue (mrfp->action_choice) == FEATURE_CONVERT)
  {
    /* convert dialog in choice */
    if (mrfp->feature_select_from != NULL)
    {
      err_list = TestDialog (mrfp->feature_select_from);
    }
    else
    {
      err_list = TestDialog (mrfp->feature_select);
    }
    if (err_list == NULL && mrfp->feature_select_to != NULL)
    {
      err_list = TestDialog (mrfp->feature_select_to); 
    }
    /* set help text */
    Reset (mrfp->help_text);
    vnp_from = (ValNodePtr) DialogToPointer (mrfp->feature_select_from);
    vnp_to = (ValNodePtr) DialogToPointer (mrfp->feature_select_to);
    if (vnp_from == NULL || vnp_to == NULL) 
    {
      AppendText (mrfp->help_text, "No conversion selected", NULL, NULL, systemFont);
    }
    else if (vnp_from->choice == vnp_to->choice)
    {
      AppendText (mrfp->help_text, "Can't convert from and to the same type", NULL, NULL, systemFont);
      ValNodeAddPointer (&err_list, 0, "Can't convert from and to the same type");
    }
    else
    {
      /* set help text */
      table_line = GetConvertFeatureTableLine (FindFeatFromFeatDefType (vnp_from->choice),
                                               vnp_from->choice,
                                               FindFeatFromFeatDefType (vnp_to->choice),
                                               vnp_to->choice);
      if (table_line < 0)
      {
        AppendText (mrfp->help_text, "Conversion not supported", NULL, NULL, systemFont);
      }
      else if (ConvertFeaturesTable[table_line].help_text == NULL)
      {
        AppendText (mrfp->help_text, "No description available", NULL, NULL, systemFont);
      }
      else
      {
        AppendText (mrfp->help_text, ConvertFeaturesTable[table_line].help_text, NULL, NULL, systemFont);
      }

      /* if converting from CDS to site, bond, or region, "leave original" must be checked */
      if ((vnp_from->choice == FEATDEF_CDS
           || vnp_from->choice == FEATDEF_ANY)
          && (vnp_to->choice == FEATDEF_SITE 
              || vnp_to->choice == FEATDEF_BOND
              || vnp_to->choice == FEATDEF_REGION))
      {
        need_original = TRUE;
        AppendText (mrfp->help_text, "If converting from CDS to site or bond, or region, you must check Leave Original Feature.", NULL, NULL, systemFont);
      }
    }
    vnp_from = ValNodeFreeData (vnp_from);
    vnp_to = ValNodeFreeData (vnp_to);
    ObjectRect (mrfp->help_text, &r);
    InvalRect (&r);  

  }
  else
  {
    /* single feature dialog in choice */
    if (mrfp->feature_select != NULL)
    {
      err_list = TestDialog (mrfp->feature_select);
    }
  }
  
  if (err_list == NULL && (!need_original || GetStatus (mrfp->leave_original_feature_btn)))
  {
    EnableAcceptCancelDialogAccept (mrfp->accept_cancel);
  }
  else
  {
    DisableAcceptCancelDialogAccept (mrfp->accept_cancel);
  }
  ValNodeFree (err_list);
}


static void FeatureRemoveChangeButton (ButtoN b)
{
  FeatureSelRemConvFormPtr mrfp;
  
  mrfp = (FeatureSelRemConvFormPtr) GetObjectExtra (b);
  if (mrfp == NULL) return;
  FeatureRemoveChangeNotify (mrfp);
}


static void FeatureRemoveOrConvertCenterAction (PopuP p)
{
  FeatureSelRemConvFormPtr mrfp;
  
  mrfp = (FeatureSelRemConvFormPtr) GetObjectExtra (p);
  if (mrfp == NULL) return;
  
  if (GetValue (p) == 2)
  {
    Hide (mrfp->remove_grp);
    Show (mrfp->convert_grp);
  }
  else
  {
    Show (mrfp->remove_grp);
    Hide (mrfp->convert_grp);
  }
  
  if (GetStatus (mrfp->clear_constraints_on_action_change))
  {
    PointerToDialog (mrfp->feature_select, NULL);
    PointerToDialog (mrfp->feature_select_from, NULL);
    PointerToDialog (mrfp->feature_select_to, NULL);
    PointerToDialog (mrfp->constraints, NULL);
  }
  FeatureRemoveChangeNotify ((Pointer) mrfp);
}



static void SelectFeatureCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)

{
  if (sfp == NULL) return;
  ObjMgrAlsoSelect (sfp->idx.entityID, sfp->idx.itemID, OBJ_SEQFEAT, 0, NULL);
}

static void DeselectFeatureCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)

{
  if (sfp == NULL) return;
  ObjMgrDeSelect (sfp->idx.entityID, sfp->idx.itemID, OBJ_SEQFEAT, 0, NULL);
}

typedef struct removefeature
{
  ValNodePtr bsplist;
  ValNodePtr bssplist;  
} RemoveFeatureData, PNTR RemoveFeaturePtr;

static void RemoveFeatureCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  RemoveFeaturePtr rfp;
  SeqIdPtr         sip;
  BioseqPtr        productbsp, productcdna;
  BioseqSetPtr     productnps;
  
  if (sfp == NULL || userdata == NULL) return;
  rfp = (RemoveFeaturePtr) userdata;
  
  sfp->idx.deleteme = TRUE;
  if (sfp->data.choice == SEQFEAT_CDREGION)
  {
    if (sfp->product != NULL)
    {
      sip = SeqLocId (sfp->product);
      if (sip != NULL)
      {
        productbsp = BioseqFind (sip);
        if (productbsp != NULL)
        {
          ValNodeAddPointer (&(rfp->bsplist), 0, (Pointer) productbsp);
        }
      }
    }
  }
  else if (sfp->data.choice == SEQFEAT_RNA)
  {
    if (sfp->product != NULL)
    {
      sip = SeqLocId (sfp->product);
      if (sip != NULL)
      {
        productcdna = BioseqFind (sip);
        if (productcdna != NULL && productcdna->idx.parenttype == OBJ_BIOSEQSET)
        {
          productnps = (BioseqSetPtr) productcdna->idx.parentptr;
          if (productnps != NULL && productnps->_class == BioseqseqSet_class_nuc_prot)
          {
            ValNodeAddPointer (&(rfp->bssplist), 0, (Pointer) productnps);
          }
        }
      }
    }
  }
}


static Boolean ConvertRegionToImp (SeqFeatPtr sfp, Uint2 featdef_to, Pointer extradata)
{
  return ConvertRegionToImpFunc (sfp, featdef_to);
}


static Boolean 
ConvertToBondSiteOrRegion 
(SeqFeatPtr sfp, 
 Int4       site_or_bond_type,
 Uint2      featdef_to,
 Boolean    create_prot_feats)
{
  BioseqPtr          bsp;
  SeqLocPtr          slp;
  SeqEntryPtr        sep;
  SeqMgrFeatContext  context;
  SeqFeatPtr         newsfp;
  SeqIdPtr           sip;
  CharPtr            comment;
  SeqBondPtr         sbp;
  SeqPntPtr          spp;
  Boolean            no_cds = FALSE;
  
  if (sfp == NULL)
  {
    return FALSE;
  }
  
  if (site_or_bond_type == -1 && featdef_to != FEATDEF_REGION)
  {
    return FALSE;
  }
  
  bsp = BioseqFindFromSeqLoc (sfp->location);
  sfp = SeqMgrGetDesiredFeature (0, bsp, 0, 0, sfp, &context);
  if (sfp == NULL || bsp == NULL) return FALSE;

  if (ISA_aa (bsp->mol))
  {
    if (create_prot_feats)
    {
      slp = (SeqLocPtr) AsnIoMemCopy (sfp->location, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
    }
    else
    {
      slp = FindNucleotideLocationForProteinFeatureConversion (sfp->location);
    }
  }
  else if (create_prot_feats) 
  {
    slp = GetProteinLocationForNucleotideFeatureConversion (sfp->location, &no_cds);
    if (no_cds) {
      Message (MSG_ERROR, "Unable to find CDS covering %s", context.label);
      return FALSE;
    }
  } else {
    slp = (SeqLocPtr) AsnIoMemCopy (sfp->location, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
  }

  bsp = GetBioseqGivenSeqLoc (slp, context.entityID);
  if (bsp == NULL) {
    Message (MSG_ERROR, "Unable to find target %s bioseq for %s", create_prot_feats ? "protein" : "nucleotide", context.label);
    slp = SeqLocFree (slp);
    return FALSE;
  } 

  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) {
    Message (MSG_ERROR, "Unable to find target %s seq-entry for %s", create_prot_feats ? "protein" : "nucleotide", context.label);
    return FALSE;
  }
  
  comment = sfp->comment;

  newsfp = NULL;
  if (featdef_to == FEATDEF_BOND || featdef_to == FEATDEF_SITE) {
    newsfp = SeqFeatNew ();
    if (newsfp != NULL) 
    {
      newsfp->data.choice = FindFeatFromFeatDefType(featdef_to);
      newsfp->data.value.intvalue = site_or_bond_type;
        sfp->comment = NULL;
    }
  } else if (featdef_to == FEATDEF_REGION) {
    newsfp = SeqFeatNew ();
    if (newsfp != NULL) 
    {
      newsfp->data.choice = FindFeatFromFeatDefType(featdef_to);
      newsfp->data.value.ptrvalue = sfp->comment;
      sfp->comment = NULL;
        comment = NULL;
    }
  }
  if (newsfp != NULL) {
    sfp->idx.deleteme = TRUE;
    CreateNewFeature (sep, NULL, FindFeatFromFeatDefType(featdef_to), newsfp);
    newsfp->location = slp;
    newsfp->comment = comment;
    if (featdef_to == FEATDEF_BOND) {
      if (slp->choice != SEQLOC_BOND) {
        sip = SeqLocId (slp);
        if (sip != NULL) {
          sbp = SeqBondNew ();
          if (sbp != NULL) {
            slp = ValNodeNew (NULL);
            if (slp != NULL) {
              slp->choice = SEQLOC_BOND;
              slp->data.ptrvalue = (Pointer) sbp;
              spp = SeqPntNew ();
              if (spp != NULL) {
                spp->strand = SeqLocStrand (newsfp->location);
                spp->id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (sip, 0)));
                spp->point = SeqLocStart (newsfp->location);
                sbp->a = spp;
              }
              spp = SeqPntNew ();
              if (spp != NULL) {
                spp->strand = SeqLocStrand (newsfp->location);
                spp->id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (sip, 0)));
                spp->point = SeqLocStop (newsfp->location);
                sbp->b = spp;
              }
              newsfp->location = SeqLocFree (newsfp->location);
              newsfp->location = slp;
            }
          }
        }
      }
    }
  }
  if (newsfp == NULL)
  {
    return FALSE;
  } 
  else
  {
    return TRUE;
  }
}


typedef struct origfeat 
{
  SeqEntryPtr sep;
  SeqFeatPtr  sfp;
} OrigFeatData, PNTR OrigFeatPtr;


static ValNodePtr AddDuplicateFeature (SeqFeatPtr sfp_orig, ValNodePtr PNTR feat_list)
{
  SeqEntryPtr sep;
  BioseqPtr   bsp;
  OrigFeatPtr ofp;
  
  if (sfp_orig == NULL || feat_list == NULL) return NULL;
  bsp = BioseqFindFromSeqLoc (sfp_orig->location);
  if (bsp == NULL) return NULL;
  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) return NULL;
  
  ofp = (OrigFeatPtr) MemNew (sizeof (OrigFeatData));
  if (ofp == NULL) return NULL;
  
  ofp->sep = sep;
  ofp->sfp = (SeqFeatPtr) AsnIoMemCopy (sfp_orig, (AsnReadFunc) SeqFeatAsnRead,
                                                 (AsnWriteFunc) SeqFeatAsnWrite);
  
  return ValNodeAddPointer (feat_list, 0, ofp);
}

static Pointer SimpleOptsFree (Pointer data)
{
  return MemFree (data);
}

static Boolean 
ConvertCDSToRNAEx 
(SeqFeatPtr  sfp,
 Uint2       toFeatSubType,
 Pointer     extradata)
{
  CDSConversionOptsPtr   opts;
  Boolean     remove_mRNA = FALSE;
  Boolean     remove_gene = FALSE;
  Boolean     remove_transcript_id = FALSE;
  Boolean     remove_only_pseudo = FALSE;


  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return FALSE;

  opts = (CDSConversionOptsPtr) extradata;
  if (opts != NULL) {
    remove_mRNA = opts->remove_mrna;
    remove_gene = opts->remove_gene;
    remove_transcript_id = opts->remove_transcript_id;
    remove_only_pseudo = opts->only_pseudo;
  }
  if (remove_only_pseudo && !sfp->pseudo) {
    return FALSE;
  }
  ApplyCDSOptionsToFeature (sfp, remove_mRNA, remove_gene, remove_transcript_id, FALSE);
  return ConvertCDSToRNA (sfp, toFeatSubType);
}


static Boolean ConvertGeneToRNAEx (SeqFeatPtr sfp, Uint2 featdef_to, Pointer extradata)
{
  return ConvertGeneToRNA (sfp, featdef_to);
}


static ValNodePtr CreateClickableListForSourceFeats (ValNodePtr feat_list)
{
  ValNodePtr transposon_only_list = NULL;
  ValNodePtr insertion_seq_only_list = NULL;
  ValNodePtr transposon_and_insertion_seq_list = NULL;
  ValNodePtr note_only_list = NULL;
  ValNodePtr no_data_list = NULL;
  ValNodePtr vnp;
  CharPtr transposon_txt, insertion_seq_txt, note_txt;
  SeqFeatPtr sfp;
  BioSourcePtr biop;
  Boolean is_transposon;
  Boolean is_insertion_seq;
  ClickableItemPtr cip;
  ValNodePtr clickable_list = NULL;
  
  for (vnp = feat_list; vnp != NULL; vnp = vnp->next) {   
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    is_transposon = FALSE;
    is_insertion_seq = FALSE;
    transposon_txt = SubSourceText (biop, SUBSRC_transposon_name, &is_transposon);
    insertion_seq_txt = SubSourceText (biop, SUBSRC_insertion_seq_name, &is_insertion_seq);
    note_txt = NoteText (biop, sfp->comment);
    
    if (is_transposon && is_insertion_seq) {
      ValNodeAddPointer (&transposon_and_insertion_seq_list, OBJ_SEQFEAT, vnp->data.ptrvalue);
    } else if (is_transposon) {
      ValNodeAddPointer (&transposon_only_list, OBJ_SEQFEAT, vnp->data.ptrvalue);
    } else if (is_insertion_seq) {
      ValNodeAddPointer (&insertion_seq_only_list, OBJ_SEQFEAT, vnp->data.ptrvalue);
    } else if (!StringHasNoText (note_txt)) {
      ValNodeAddPointer (&note_only_list, OBJ_SEQFEAT, vnp->data.ptrvalue);
    } else {
      ValNodeAddPointer (&no_data_list, OBJ_SEQFEAT, vnp->data.ptrvalue);
    }
  }

  if (transposon_and_insertion_seq_list != NULL) {
    cip = NewClickableItem (0, "%d source features have both transposon-name and insertion-seq-name.", transposon_and_insertion_seq_list);
    ValNodeAddPointer (&clickable_list, 0, cip);
  }
  if (transposon_only_list != NULL) {
    cip = NewClickableItem (0, "%d source features have transposon-name.", transposon_only_list);
    cip->chosen = TRUE;
    ValNodeAddPointer (&clickable_list, 0, cip);
  }
  if (insertion_seq_only_list != NULL) {
    cip = NewClickableItem (0, "%d source features have insertion-seq-name.", insertion_seq_only_list);
    cip->chosen = TRUE;
    ValNodeAddPointer (&clickable_list, 0, cip);
  }
  if (note_only_list != NULL) {
    cip = NewClickableItem (0, "%d source features have notes, but neither transposon-name nor insertion-seq-name.", note_only_list);
    ValNodeAddPointer (&clickable_list, 0, cip);
  }
  if (no_data_list != NULL) {
    cip = NewClickableItem (0, "%d source features have no notes, no transposon-name, and no insertion-seq-name.", no_data_list);
    ValNodeAddPointer (&clickable_list, 0, cip);
  }

  return clickable_list;
}


static ValNodePtr AdjustConvertFeaturesForBioSrcToRepeatRegion (ValNodePtr feat_list, Pointer extradata, BaseFormPtr bfp, BoolPtr cancel)
{
  ValNodePtr clickable_list;

  *cancel = FALSE;
  clickable_list = CreateClickableListForSourceFeats (feat_list);
  feat_list = ValNodeFree (feat_list);
  feat_list = ChooseFeaturesForConversion (clickable_list, bfp, "Available Data", "Features");
  /* NOTE - the clickable_list will be freed by the dialog  - do not free it again here */
  if (feat_list == NULL)
  {
    *cancel = TRUE;
  }
  return feat_list;
}


static void RemoveNonPseudoFeatures (ValNodePtr PNTR feat_list)
{
  ValNodePtr prev = NULL, next_vnp, vnp;
  SeqFeatPtr sfp;

  vnp = *feat_list;
  while (vnp != NULL) {
    next_vnp = vnp->next;
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    if (sfp != NULL && !sfp->pseudo) {
      if (prev == NULL) {
        *feat_list = vnp->next;
      } else {
        prev->next = vnp->next;
      }
      vnp->next = NULL;
      vnp = ValNodeFree (vnp);
    } else {
      prev = vnp;
    }
    vnp = next_vnp;
  }    
}


static Pointer GetCDSConversionOpts (ValNodePtr feat_list, BaseFormPtr bfp, BoolPtr cancel)
{
  Boolean all_are_pseudo = TRUE, any_are_pseudo = FALSE, any_gps = FALSE;
  CDSConversionOptsPtr opts;
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  
  for (vnp = feat_list; vnp != NULL; vnp = vnp->next)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    if (sfp != NULL)
    {
      if (sfp->pseudo)
      {
        any_are_pseudo = TRUE;
      }
      else
      {
        all_are_pseudo = FALSE;
      }
      if (!any_gps)
      {
        any_gps = IsFeatInGPS (sfp);
      }
    }
  }

  opts = GetCDSConversionOptions (all_are_pseudo, any_are_pseudo, any_gps, cancel);

  return opts;
}


static ValNodePtr AdjustCDSToMiscFeatList (ValNodePtr feat_list, Pointer extradata, BaseFormPtr bfp, BoolPtr cancel)
{
  CDSConversionOptsPtr opts;

  if (feat_list == NULL || extradata == NULL || cancel == NULL) return feat_list;
  opts = (CDSConversionOptsPtr) extradata;
  *cancel = FALSE;
  if (opts->only_pseudo) {
    /* remove features that are not pseudo */
    RemoveNonPseudoFeatures (&feat_list);
  }
  return feat_list;
}


static void CDSToMiscFeatConversionCleanup (Uint2 entityID, Pointer extradata)
{
  SeqEntryPtr sep;

  DeleteMarkedObjects (entityID, 0, NULL);
  sep = GetTopSeqEntryForEntityID (entityID);
  RenormalizeNucProtSets (sep, TRUE);  
}


static Boolean CDSToMiscFeatConvertFunc (SeqFeatPtr sfp, Uint2 featdef_to, Pointer extradata)
{
  CDSConversionOptsPtr opts;
  Boolean              rval = FALSE;

  opts = (CDSConversionOptsPtr) extradata;
  if (opts == NULL || sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION)
  {
    return FALSE;
  }
  else if (sfp->pseudo) 
  {
    rval = ConvertOnePseudoCDSToMiscFeat (sfp);
  }
  else
  {
    /* do other here */
  
    rval = ConvertOneCDSToMiscFeat (sfp, FALSE, FALSE, opts);
  }
  return TRUE;
}


typedef struct siteorbond {
  Int4 site_or_bond_type;
} SiteOrBondData, PNTR SiteOrBondPtr;

static Pointer GetSiteOrBondOpts (ValNodePtr feat_list, BaseFormPtr bfp, BoolPtr cancel, Uint2 toFeat)
{
  EnumFieldAssocPtr  al;
  CharPtr            mssg;
  UIEnum             val;
  SiteOrBondPtr      opts = NULL;

  *cancel = FALSE;
  al = NULL;
  mssg = NULL;
  if (toFeat == FEATDEF_BOND) {
    al = enum_bond_alist;
    mssg = "Select type of bond";
  } else if (toFeat == FEATDEF_SITE) {
    al = enum_site_alist;
    mssg = "Select type of site";
  } else {
    *cancel = TRUE;
    return NULL;
  }

  if (al == NULL || ! AlistMessage (al, &val, 1, mssg)) 
  {
    *cancel = TRUE;
  }
  else
  {
    opts = (SiteOrBondPtr) MemNew (sizeof (SiteOrBondData));
    opts->site_or_bond_type = (Int4) val; 
  }
  return opts;
}


static Pointer GetSiteOpts (ValNodePtr feat_list, BaseFormPtr bfp, BoolPtr cancel)
{
  return GetSiteOrBondOpts (feat_list, bfp, cancel, FEATDEF_SITE);
}


static Pointer GetBondOpts (ValNodePtr feat_list, BaseFormPtr bfp, BoolPtr cancel)
{
  return GetSiteOrBondOpts (feat_list, bfp, cancel, FEATDEF_BOND);
}


static Boolean ConvertToSiteOrBond (SeqFeatPtr sfp, Uint2 featdef_to, Pointer extradata)
{
  SiteOrBondPtr      opts = NULL;

  opts = (SiteOrBondPtr) extradata;
  if (sfp == NULL || opts == NULL) return FALSE;

  return ConvertToBondSiteOrRegion (sfp, opts->site_or_bond_type, featdef_to, TRUE);
}

typedef struct regionopts {
  Boolean create_prot_feats;
} RegionOptsData, PNTR RegionOptsPtr;


static ValNodePtr GetBadProteinConversionFeatures (ValNodePtr feat_list)
{
  ValNodePtr no_cds_list = NULL, bad_coordinates_list = NULL, err_list = NULL, vnp;
  Boolean    no_cds;
  SeqLocPtr  slp;
  SeqFeatPtr sfp;
  BioseqPtr  bsp;

  for (vnp = feat_list; vnp != NULL; vnp = vnp->next) {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL && ISA_aa (bsp->mol))
    {
      slp = FindNucleotideLocationForProteinFeatureConversion (sfp->location);
      if (slp == NULL)
      {
        no_cds = TRUE;
      }
    }
    else
    {
      slp = GetProteinLocationForNucleotideFeatureConversion (sfp->location, &no_cds);
    }
    if (slp == NULL) {
      if (no_cds) {
        ValNodeAddPointer (&no_cds_list, OBJ_SEQFEAT, sfp);
      } else {
        ValNodeAddPointer (&bad_coordinates_list, OBJ_SEQFEAT, sfp);
      }
    } else {
      slp = SeqLocFree (slp);
    }
  }

  if (no_cds_list != NULL) {
    ValNodeAddPointer (&err_list, 0, NewClickableItem (0, "%d features have no overlapping coding region", no_cds_list));
  }
  if (bad_coordinates_list != NULL) {
    ValNodeAddPointer (&err_list, 0, NewClickableItem (0, "%d features do not stop or start on codon boundaries", bad_coordinates_list));
  }
  return err_list;
}


static Pointer GetRegionOpts (ValNodePtr feat_list, BaseFormPtr bfp, BoolPtr cancel)
{
  ValNodePtr            err_list;
  WindoW                w;
  GrouP                 h, nuc_prot_choice, c;
  DialoG                bad_prot_feats = NULL;
  ButtoN                b;
  ModalAcceptCancelData acd;
  Boolean               rval = FALSE;
  RegionOptsPtr         opts = NULL;
  
  *cancel = FALSE;

  if (feat_list == NULL) return NULL;
  
  err_list = GetBadProteinConversionFeatures (feat_list);
    
  w = MovableModalWindow(-20, -13, -10, -10, "Create Nucleotide or Protein Features", NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  if (err_list != NULL) {
    bad_prot_feats = CreateClickableListDialog (h, "Reasons", "Original Features",
                                                ScrollToDiscrepancyItem, EditDiscrepancyItem, NULL,
                                                GetDiscrepancyItemText);
    PointerToDialog (bad_prot_feats, err_list);
  }

  nuc_prot_choice = HiddenGroup (h, 0, 2, NULL);
  RadioButton (nuc_prot_choice, "Create Nucleotide Features");
  RadioButton (nuc_prot_choice, err_list == NULL ? "Create Protein Features" : "Create Protein Features (features listed above will be skipped)");
  SetValue (nuc_prot_choice, err_list == NULL ? 2 : 1);
  
  c = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) nuc_prot_choice,
                              (HANDLE) c, 
                              (HANDLE) bad_prot_feats, 
                              NULL);
  Show (w);
  Select (w);
  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  if (acd.cancelled)
  {
    *cancel = TRUE;
  }
  else
  {
    opts = (RegionOptsPtr) MemNew (sizeof (RegionOptsData));
    if (GetValue (nuc_prot_choice) == 2)
    {
      opts->create_prot_feats = TRUE;
    }
    else
    {
      opts->create_prot_feats = FALSE;
    }
  }
  Remove (w);
  
  return opts;
}


static void RemoveClickableItemFeaturesFromFeatureList (ClickableItemPtr cip, ValNodePtr PNTR feat_list)
{
  ValNodePtr vnp_i, vnp_p, vnp_r;

  if (cip == NULL || feat_list == NULL || *feat_list == NULL) return;
  for (vnp_i = cip->item_list; vnp_i != NULL; vnp_i = vnp_i->next) {
    if (vnp_i->choice != OBJ_SEQFEAT || vnp_i->data.ptrvalue == NULL) continue;
    vnp_p = NULL;
    vnp_r = *feat_list;
    while (vnp_r != NULL && vnp_r->data.ptrvalue != vnp_i->data.ptrvalue) {
      vnp_p = vnp_r;
      vnp_r = vnp_r->next;
    }
    if (vnp_r != NULL) {
      if (vnp_p == NULL) {
        *feat_list = vnp_r->next;
      } else {
        vnp_p->next = vnp_r->next;
      }
      vnp_r->next = NULL;
      vnp_r = ValNodeFree (vnp_r);
    }
  }
  for (vnp_i = cip->subcategories; vnp_i != NULL; vnp_i = vnp_i->next) {
    RemoveClickableItemFeaturesFromFeatureList (vnp_i->data.ptrvalue, feat_list);
  }
}


static ValNodePtr AdjustRegionConversionFeatures (ValNodePtr feat_list, Pointer extradata, BaseFormPtr bfp, BoolPtr cancel)
{
  RegionOptsPtr opts = (RegionOptsPtr) extradata;
  ValNodePtr    err_list, vnp_c;

  *cancel = FALSE;
  if (opts == NULL || !opts->create_prot_feats) return feat_list;
  
  err_list = GetBadProteinConversionFeatures (feat_list);
  
  for (vnp_c = err_list; vnp_c != NULL; vnp_c = vnp_c->next) {
    RemoveClickableItemFeaturesFromFeatureList (vnp_c->data.ptrvalue, &feat_list);
  }
  return feat_list;  
}


static Boolean ConvertAnyToRegion (SeqFeatPtr sfp, Uint2 featdef_to, Pointer extradata)
{
  RegionOptsPtr opts = (RegionOptsPtr) extradata;

  if (opts == NULL || sfp == NULL) return FALSE;
  return ConvertToBondSiteOrRegion (sfp, 0, FEATDEF_REGION, opts->create_prot_feats);
}


static Boolean ConvertImpToImp (SeqFeatPtr sfp, Uint2 featdef_to, Pointer extradata)
{
  return ConvertImpToImpFunc (sfp, featdef_to);
}


static void CollectFeaturesToBeConverted (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  ValNodePtr PNTR feat_list = (ValNodePtr PNTR) userdata;
  
  ValNodeAddPointer (feat_list, OBJ_SEQFEAT, sfp);
}


static void MoveInferenceToOverlappingGene (SeqFeatPtr sfp)
{
  GBQualPtr qual, qual_prev= NULL, qual_next;
  SeqFeatPtr gene;

  if (sfp == NULL) return;

  gene = GetGeneForFeature (sfp);
  if (gene == NULL || gene->idx.deleteme || gene == sfp) return;

  qual = sfp->qual;
  while (qual != NULL) {
    qual_next = qual->next;
    if (StringCmp (qual->qual, "inference") == 0) {
      if (qual_prev == NULL) {
        sfp->qual = qual->next;
      } else {
        qual_prev->next = qual->next;
      }
      qual->next = NULL;
      qual->next = gene->qual;
      gene->qual = qual;
    } else {
      qual_prev = qual;
    }
    qual = qual_next;
  }
}


static Boolean NewConvertFeaturesByList (SeqEntryPtr sep, FilterSetPtr fsp, Uint1 featdef_from, FeatureSelRemConvFormPtr mrfp)
{
  ValNodePtr feat_list = NULL, sel_feat_list = NULL, vnp;
  Boolean not_match = FALSE;
  Boolean    only_pseudo = FALSE;
  SelStructPtr sel;
  SeqFeatPtr   sfp;
  SeqMgrFeatContext fcontext;
  Int4              table_line_num = -1;
  Uint2             seqfeat_from, seqfeat_to;
  Pointer           extradata = NULL;
  Boolean           cancel = FALSE;
  ValNodePtr        dup_feat_vnp = NULL;
  OrigFeatPtr       ofp;
  Boolean           removed_cds = FALSE;
  BioseqPtr         protBsp;

  OperateOnSeqEntryConstrainedObjects (sep, fsp, 
                                       CollectFeaturesToBeConverted,
                                       NULL, 0, featdef_from, 0, &feat_list);
  sel = ObjMgrGetSelected ();
  while (sel != NULL) {
    if (sel->entityID == mrfp->input_entityID && sel->itemtype == OBJ_SEQFEAT) {
      sfp = SeqMgrGetDesiredFeature (sel->entityID, NULL, sel->itemID, 0, NULL, &fcontext);
      if (sfp != NULL) {
        vnp = feat_list; 
        while (vnp != NULL && vnp->data.ptrvalue != sfp) {
          vnp = vnp->next;
        }
        if (vnp == NULL) {
          not_match = TRUE;
        }
        ValNodeAddPointer (&sel_feat_list, OBJ_SEQFEAT, sfp);
      }
    }
    sel = sel->next;
  }
  if (sel_feat_list != NULL) {
    feat_list = ValNodeFree (feat_list);
    if (not_match 
        && Message (MSG_YN, 
                    "Selected feature%s do%s not match options - apply to selected feature anyway?",
                    sel_feat_list->next == NULL ? "" : "s",
                    sel_feat_list->next == NULL ? "es" : "") == ANS_NO) {
      sel_feat_list = ValNodeFree (sel_feat_list);
      return FALSE;
    }
    feat_list = sel_feat_list;
  }

  if (feat_list == NULL) {
    Message (MSG_ERROR, "No features of the specified type exist to be converted!");
    return FALSE;
  }

  /* find line in table that applies */
  seqfeat_from = FindFeatFromFeatDefType (featdef_from);
  seqfeat_to =   FindFeatFromFeatDefType (mrfp->featdef_to);

  table_line_num = GetConvertFeatureTableLine (seqfeat_from, featdef_from, seqfeat_to, mrfp->featdef_to);
  if (table_line_num < 0 || ConvertFeaturesTable[table_line_num].convert_func == NULL)
  {
    Message (MSG_ERROR, "This conversion not supported - contact"
         " sequindev for instructions.");
    return FALSE;
  }

  /* collect extra data if necessary */
  if (ConvertFeaturesTable[table_line_num].collect_options != NULL)
  {
    extradata = (ConvertFeaturesTable[table_line_num].collect_options) (feat_list, mrfp->bfp, &cancel);
    if (cancel)
    {
      feat_list = ValNodeFree (feat_list);
      /* free extradata */
      if (ConvertFeaturesTable[table_line_num].free_options != NULL)
      {
        (ConvertFeaturesTable[table_line_num].free_options) (extradata);
      }
      return FALSE;
    }
  }
  else if (featdef_from == FEATDEF_CDS && mrfp->featdef_to == FEATDEF_mat_peptide_aa) 
  {
    extradata = &(mrfp->leave_original_feature);
  }

  /* adjust feature list if necessary */  
  if (ConvertFeaturesTable[table_line_num].adjust_features != NULL) 
  {
    feat_list = (ConvertFeaturesTable[table_line_num].adjust_features) (feat_list, extradata, mrfp->bfp, &cancel);
    if (cancel || feat_list == NULL)
    {
      feat_list = ValNodeFree (feat_list);
      /* free extradata */
      if (ConvertFeaturesTable[table_line_num].free_options != NULL)
      {
        (ConvertFeaturesTable[table_line_num].free_options) (extradata);
      }
      return FALSE;
    }
  }

  /* process list */
  for (vnp = feat_list; vnp != NULL; vnp = vnp->next)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    if (sfp == NULL) continue;
    protBsp = NULL;
    if (sfp->data.choice == SEQFEAT_CDREGION && !mrfp->leave_original_feature)
    {
      protBsp = BioseqFindFromSeqLoc (sfp->product);
      if (protBsp != NULL)
      {
        protBsp->idx.deleteme = TRUE;
      }
      removed_cds = TRUE;
    }
      
    dup_feat_vnp = NULL;
    if (mrfp->leave_original_feature)
    {
      /* need a list of duplicate features that we can add to the appropriate SeqEntrys
       * when we are done with the conversions.
       */
      dup_feat_vnp = AddDuplicateFeature (sfp, &(mrfp->feat_list)); 
    }

    if (ConvertFeaturesTable[table_line_num].convert_func (sfp, mrfp->featdef_to, extradata))
    {
      /* set the subtype to zero so that it will be reindexed */
      sfp->idx.subtype = 0;
      if (mrfp->featdef_to == FEATDEF_misc_feature)
      {
        MoveInferenceToOverlappingGene (sfp);
      }
    }
    else if (dup_feat_vnp != NULL && dup_feat_vnp->data.ptrvalue != NULL)
    {
      /* feature wasn't removed, for whatever reason */
      /* remove the duplicate feature from the list to be added later */
      ofp = (OrigFeatPtr) dup_feat_vnp->data.ptrvalue;
      ofp->sfp = SeqFeatFree (ofp->sfp);
      ofp->sep = NULL;

      /* don't remove protein product if feature wasn't removed */
      if (protBsp != NULL)
      {
        protBsp->idx.deleteme = FALSE;
      }
    }
    else if (sfp->data.choice == SEQFEAT_CDREGION) 
    {
      /* failed to convert feature, don't remove product */
      if (protBsp != NULL)
      {
        protBsp->idx.deleteme = FALSE;
      }
    }
  }

  feat_list = ValNodeFree (feat_list);

  /* do cleanup tasks */
  if (ConvertFeaturesTable[table_line_num].cleanup_func != NULL)
  {
    (ConvertFeaturesTable[table_line_num].cleanup_func) (mrfp->input_entityID, extradata);
  }

  /* free extradata */
  if (ConvertFeaturesTable[table_line_num].free_options != NULL)
  {
    (ConvertFeaturesTable[table_line_num].free_options) (extradata);
  }

  if (removed_cds)
  {
    DeleteMarkedObjects (mrfp->input_entityID, 0, NULL);
    RenormalizeNucProtSets (sep, TRUE);       
  }

  return TRUE;
}


static void RemoveProteinProducts (ValNodePtr bsplist, Uint2 entityID, SeqEntryPtr sep)
{
  MsgAnswer  ans;
  ValNodePtr tmp;
  BioseqPtr  bsp;
  Uint4      itemID;
  OMProcControl  ompc;
  
  if (bsplist == NULL) 
  {
    return;
  }
  ans = Message (MSG_YN, "Remove protein products?");
  if (ans == ANS_YES) {
    for (tmp = bsplist; tmp != NULL; tmp = tmp->next) {
      bsp = (BioseqPtr) tmp->data.ptrvalue;
      itemID = GetItemIDGivenPointer (entityID, OBJ_BIOSEQ, (Pointer) bsp);
      if (itemID > 0) {
        MemSet ((Pointer) (&ompc), 0, sizeof (OMProcControl));
        ompc.do_not_reload_from_cache = TRUE;
        ompc.input_entityID = entityID;
        ompc.input_itemID = itemID;
        ompc.input_itemtype = OBJ_BIOSEQ;
        if (! DetachDataForProc (&ompc, FALSE)) {
          Message (MSG_POSTERR, "DetachDataForProc failed");
        }
        SeqMgrDeleteFromBioseqIndex (bsp);
      }
    }
    ans = Message (MSG_YN, "Renormalize Nuc-Prot sets?");
    if (ans == ANS_YES)
    {
      RemoveOrphanProteins (entityID, sep);
      RenormalizeNucProtSets (sep, TRUE);
    }
  }
}

static void RemoveCDNANucProtProducts (ValNodePtr bssplist, Uint2 entityID)
{
  MsgAnswer     ans;
  ValNodePtr    tmp;
  BioseqSetPtr  bssp;
  Uint4         itemID;
  OMProcControl ompc;

  if (bssplist == NULL)
  {
    return;
  }
  
  ans = Message (MSG_YN, "Remove cDNA nuc-prot products?");
  if (ans != ANS_YES)
  {
    return;
  }

  for (tmp = bssplist; tmp != NULL; tmp = tmp->next) 
  {
    bssp = (BioseqSetPtr) tmp->data.ptrvalue;
    itemID = GetItemIDGivenPointer (entityID, OBJ_BIOSEQSET, (Pointer) bssp);
    if (itemID > 0) {
      MemSet ((Pointer) (&ompc), 0, sizeof (OMProcControl));
      ompc.do_not_reload_from_cache = TRUE;
      ompc.input_entityID = entityID;
      ompc.input_itemID = itemID;
      ompc.input_itemtype = OBJ_BIOSEQSET;
      if (! DetachDataForProc (&ompc, FALSE)) {
        Message (MSG_POSTERR, "DetachDataForProc failed");
      }
    }
  }
}


static Boolean FeatureRemoveOrConvertAction (Pointer userdata)
{
  FeatureSelRemConvFormPtr   mrfp;
  FilterSetPtr          fsp;
  SeqEntryPtr           sep;
  ValNodePtr            feature_type_list, vnp;
  Int4                  window_action;
  Uint1                 feat_def_choice;
  RemoveFeatureData     rfd;
  OrigFeatPtr           ofp;
  SeqFeatPtr            sfp;
  Boolean               rval = TRUE;
  SeqEntryPtr           create_sep;
  
  if (userdata == NULL) return FALSE;
  
  mrfp = (FeatureSelRemConvFormPtr) userdata;
  window_action = GetValue (mrfp->action_choice);
  
  fsp = (FilterSetPtr) DialogToPointer (mrfp->constraints);
  
  sep = GetTopSeqEntryForEntityID (mrfp->input_entityID);
  if (sep == NULL) return FALSE;
  
  if (window_action == FEATURE_CONVERT)
  {
    vnp = (ValNodePtr) DialogToPointer (mrfp->feature_select_to);
    if (vnp == NULL) 
    {
      FilterSetFree (fsp);
      return FALSE;
    }
    mrfp->featdef_to = vnp->choice;
    /* NOTE - I do not need to use ValNodeFreeData because I'm using
     * the data.ptrvalue in mrfp->featname_to and will free it when I free
     * mrfp.
     */
    vnp = ValNodeFree (vnp);
    feature_type_list = (ValNodePtr) DialogToPointer (mrfp->feature_select_from);
    mrfp->leave_original_feature = GetStatus (mrfp->leave_original_feature_btn);
  }
  else
  {
    feature_type_list = (ValNodePtr) DialogToPointer (mrfp->feature_select);
  }
  if (feature_type_list == NULL)
  {
    FilterSetFree (fsp);
    return FALSE;
  }
    
  for (vnp = feature_type_list; vnp != NULL; vnp = vnp->next)
  {
    feat_def_choice = vnp->choice;
    if (feat_def_choice == 255)
    {
      feat_def_choice = 0;
    }
    switch (window_action)
    {
      case FEATURE_CONVERT:
        rval = NewConvertFeaturesByList (sep, fsp, feat_def_choice, mrfp);
        break;
      case FEATURE_REMOVE:
        rfd.bsplist = NULL;
        rfd.bssplist = NULL;
        OperateOnSeqEntryConstrainedObjects (sep, fsp, 
                                             RemoveFeatureCallback,
                                             NULL, 0, feat_def_choice, 0, &rfd);
        RemoveProteinProducts (rfd.bsplist, mrfp->input_entityID, sep);
        RemoveCDNANucProtProducts (rfd.bssplist, mrfp->input_entityID);
        rfd.bsplist = ValNodeFree (rfd.bsplist);
        rfd.bssplist = ValNodeFree (rfd.bssplist);
        break;
      case FEATURE_SELECT:
        OperateOnSeqEntryConstrainedObjects (sep, fsp, 
                                             SelectFeatureCallback,
                                             NULL, 0, feat_def_choice, 0, NULL);
        break;
      case FEATURE_DESELECT:
        OperateOnSeqEntryConstrainedObjects (sep, fsp, 
                                             DeselectFeatureCallback,
                                             NULL, 0, feat_def_choice, 0, NULL);
        break;
    }
  }
  
  ValNodeFree (feature_type_list);
  FilterSetFree (fsp);
  
  if (window_action == FEATURE_REMOVE || window_action == FEATURE_CONVERT)
  {
    DeleteMarkedObjects (mrfp->input_entityID, 0, NULL);
  }
  
  if (window_action == FEATURE_CONVERT && mrfp->leave_original_feature)
  {
    for (vnp = mrfp->feat_list; vnp != NULL; vnp = vnp->next)
    {
      ofp = (OrigFeatPtr) vnp->data.ptrvalue;
      if (ofp == NULL || ofp->sep == NULL || ofp->sfp == NULL)
      {
        continue;
      }
      if (IS_Bioseq_set (ofp->sep)) {
        create_sep = FindNucSeqEntry (ofp->sep);
      } else {
        create_sep = ofp->sep;
      }
      sfp = CreateNewFeature (create_sep, NULL, ofp->sfp->data.choice, ofp->sfp);
    }
  }
  mrfp->feat_list = ValNodeFreeData (mrfp->feat_list);

  if (mrfp->renormalize_nucprot) {
    RenormalizeNucProtSets (sep, TRUE);       
  }
  
  ObjMgrSetDirtyFlag (mrfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, mrfp->input_entityID, 0, 0);  
  Update ();
  return rval;
}

static void FeatureRemoveOrConvert (IteM i, Int4 first_action)
{
  BaseFormPtr         bfp;
  FeatureSelRemConvFormPtr mrfp;
  WindoW              w;
  GrouP               h, k, g, n;
  SeqEntryPtr         sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  mrfp = (FeatureSelRemConvFormPtr) MemNew (sizeof (FeatureSelRemConvFormData));
  if (mrfp == NULL) return;
  
  mrfp->bfp = bfp;
  
  w = FixedWindow (-50, -33, -10, -10, "Remove or Convert Features", StdCloseWindowProc);
  SetObjectExtra (w, mrfp, CleanupFeatureSelRemConvForm);
  mrfp->form = (ForM) w;
  mrfp->input_entityID = bfp->input_entityID;
  
  sep = GetTopSeqEntryForEntityID(bfp->input_entityID);
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  if (first_action < FEATURE_REMOVE || first_action > FEATURE_DESELECT)
  {
    first_action = FEATURE_REMOVE;
  }
  mrfp->action_choice = PopupList (h, TRUE, FeatureRemoveOrConvertCenterAction);
  PopupItem (mrfp->action_choice, "Remove Features");
  PopupItem (mrfp->action_choice, "Convert Features");
  SetValue (mrfp->action_choice, first_action);
  SetObjectExtra (mrfp->action_choice, mrfp, NULL);
  
  mrfp->clear_constraints_on_action_change = CheckBox (h, "Clear when changing actions", NULL);
  SetStatus (mrfp->clear_constraints_on_action_change, TRUE);

  n = HiddenGroup (h, 0, 0, NULL);
  mrfp->remove_grp = HiddenGroup (n, 0, 1, NULL);
  mrfp->feature_select =  FeatureSelectionDialogEx (mrfp->remove_grp, TRUE, sep,
                                                  FeatureRemoveChangeNotify, 
                                                  mrfp);
  
  mrfp->convert_grp = HiddenGroup (n, 3, 0, NULL);
  k = HiddenGroup (mrfp->convert_grp, 0, 2, NULL);
  mrfp->from_prompt = StaticPrompt (k, "From", 0, dialogTextHeight, systemFont, 'l');
  mrfp->feature_select_from =  FeatureSelectionDialogEx (k, FALSE, sep,
                                                  FeatureRemoveChangeNotify, 
                                                  mrfp);
  AlignObjects (ALIGN_CENTER, (HANDLE) mrfp->from_prompt, (HANDLE) mrfp->feature_select_from, NULL);
  k = HiddenGroup (mrfp->convert_grp, 0, 2, NULL);
  mrfp->to_prompt = StaticPrompt (k, "To", 0, dialogTextHeight, systemFont, 'l');
  mrfp->feature_select_to =  FeatureSelectionDialog (k, FALSE,
                                                  FeatureRemoveChangeNotify, 
                                                  mrfp);
  AlignObjects (ALIGN_CENTER, (HANDLE) mrfp->to_prompt, (HANDLE) mrfp->feature_select_to, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) mrfp->remove_grp, (HANDLE) mrfp->convert_grp, NULL);

  g = HiddenGroup (mrfp->convert_grp, 0, 2, NULL);
  StaticPrompt (g, "Conversion Function", 0, dialogTextHeight, systemFont, 'c');
  SelectFont (systemFont);
  mrfp->help_text = DocumentPanel (g, stdCharWidth * 12, LineHeight () * (TALL_SELECTION_LIST + 1));

  mrfp->leave_original_feature_btn = CheckBox (mrfp->convert_grp, "Leave Original Feature", FeatureRemoveChangeButton);
  SetObjectExtra (mrfp->leave_original_feature_btn, mrfp, NULL);
  SetStatus (mrfp->leave_original_feature_btn, FALSE);
  AlignObjects (ALIGN_CENTER, (HANDLE) k, (HANDLE) mrfp->leave_original_feature_btn, NULL);
  
  mrfp->feat_list = NULL;
  mrfp->constraints = FilterGroup (h, TRUE, FALSE, TRUE, FALSE, TRUE, "Where feature text");
  mrfp->accept_cancel = AcceptCancelDialog (h, FeatureRemoveOrConvertAction, NULL, FeatureRemoveClear, FeatureRemoveClearText, (Pointer)mrfp, w);
  AlignObjects (ALIGN_CENTER, (HANDLE) mrfp->action_choice,
                              (HANDLE) mrfp->clear_constraints_on_action_change,
                              (HANDLE) n,
                              (HANDLE) mrfp->constraints,
                              (HANDLE) mrfp->accept_cancel, NULL);
                              
  FeatureRemoveOrConvertCenterAction (mrfp->action_choice);
  
  Show (w);
}

extern void FeatureRemove (IteM i)
{
  FeatureRemoveOrConvert (i, FEATURE_REMOVE);
}

extern void ConvertFeatures (IteM i)
{
  FeatureRemoveOrConvert (i, FEATURE_CONVERT);
}

static Boolean SelectFeatureAction (Pointer userdata)
{
  FeatureSelRemConvFormPtr   mrfp;
  FilterSetPtr          fsp;
  SeqEntryPtr           sep;
  ValNodePtr            feature_type_list, vnp;
  Uint1                 feat_def_choice;
  
  if (userdata == NULL) return FALSE;
  
  mrfp = (FeatureSelRemConvFormPtr) userdata;
  fsp = (FilterSetPtr) DialogToPointer (mrfp->constraints);
  
  sep = GetTopSeqEntryForEntityID (mrfp->input_entityID);
  if (sep == NULL) return FALSE;
  
  feature_type_list = (ValNodePtr) DialogToPointer (mrfp->feature_select);
  if (feature_type_list == NULL)
  {
    FilterSetFree (fsp);
    return FALSE;
  }
    
  for (vnp = feature_type_list; vnp != NULL; vnp = vnp->next)
  {
    feat_def_choice = vnp->choice;
    if (feat_def_choice == 255)
    {
      feat_def_choice = 0;
    }
    OperateOnSeqEntryConstrainedObjects (sep, fsp, 
                                         SelectFeatureCallback,
                                         NULL, 0, feat_def_choice, 0, NULL);
  }
  
  ValNodeFree (feature_type_list);
  FilterSetFree (fsp);
  
  ObjMgrSetDirtyFlag (mrfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, mrfp->input_entityID, 0, 0);  
  Update ();
  return TRUE;
}

extern void SelectFeatures (IteM i)
{
  BaseFormPtr              bfp;
  FeatureSelRemConvFormPtr mrfp;
  WindoW                   w;
  GrouP                    h;
  SeqEntryPtr              sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  mrfp = (FeatureSelRemConvFormPtr) MemNew (sizeof (FeatureSelRemConvFormData));
  if (mrfp == NULL) return;
  
  w = FixedWindow (-50, -33, -10, -10, "Select Features", StdCloseWindowProc);
  SetObjectExtra (w, mrfp, CleanupFeatureSelRemConvForm);
  mrfp->form = (ForM) w;
  mrfp->input_entityID = bfp->input_entityID;
  sep = GetTopSeqEntryForEntityID(bfp->input_entityID);
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  mrfp->feature_select =  FeatureSelectionDialogEx (h, TRUE, sep,
                                                  FeatureRemoveChangeNotify, 
                                                  mrfp);
    
  mrfp->constraints = FilterGroup (h, TRUE, FALSE, TRUE, FALSE, FALSE, "Where feature text");
  mrfp->accept_cancel = AcceptCancelDialog (h, SelectFeatureAction, NULL, FeatureRemoveClear, FeatureRemoveClearText, (Pointer)mrfp, w);
  AlignObjects (ALIGN_CENTER, (HANDLE) mrfp->feature_select,
                              (HANDLE) mrfp->constraints,
                              (HANDLE) mrfp->accept_cancel, NULL);
                              
  Show (w);
}


typedef struct reverseintervalsform
{
  FORM_MESSAGE_BLOCK
  DialoG  feature_select;
  DialoG  constraints;
  DialoG  accept_cancel;

  ValNodePtr feat_list;
  BaseFormPtr bfp;
} ReverseIntervalsFormData, PNTR ReverseIntervalsFormPtr;

static void ReverseIntervalsChangeNotify (Pointer userdata)
{
  ReverseIntervalsFormPtr fp;
  ValNodePtr          err_list;

  fp = (ReverseIntervalsFormPtr) userdata;
  if (fp == NULL) return;

  err_list = TestDialog (fp->feature_select);
  
  if (err_list == NULL)
  {
    EnableAcceptCancelDialogAccept (fp->accept_cancel);
  }
  else
  {
    DisableAcceptCancelDialogAccept (fp->accept_cancel);
  }
  ValNodeFree (err_list);
}

static void ReverseIntervalsClearText (Pointer data)
{
  ReverseIntervalsFormPtr   mrfp;
  FilterSetPtr          fsp;

  mrfp = (ReverseIntervalsFormPtr) data;
  if (mrfp == NULL) return;
 
  fsp = DialogToPointer (mrfp->constraints);
  FilterSetClearText (fsp);
  PointerToDialog (mrfp->constraints, fsp);
  FilterSetFree (fsp);
}

static void ReverseIntervalsClear (Pointer data)
{
  ReverseIntervalsFormPtr mrfp;

  mrfp = (ReverseIntervalsFormPtr) data;
  if (mrfp == NULL) return;
 
  PointerToDialog (mrfp->feature_select, NULL);
  PointerToDialog (mrfp->constraints, NULL);
}

static void ReverseLocationIntervalOrder (SeqLocPtr slp)
{
  SeqLocPtr     subslp_list = NULL, subslp, subslp_next;
  PackSeqPntPtr pspp;
  Int4 pnt_num, swap;

  if (slp == NULL 
      || slp->choice == SEQLOC_NULL
      || slp->choice == SEQLOC_EMPTY
      || slp->choice == SEQLOC_WHOLE
      || slp->choice == SEQLOC_INT
      || slp->choice == SEQLOC_PNT
      || slp->choice == SEQLOC_BOND
      || slp->choice == SEQLOC_FEAT) {
    return;
  }
  if (slp->choice == SEQLOC_MIX
      || slp->choice == SEQLOC_PACKED_INT
      || slp->choice == SEQLOC_EQUIV) {
    for (subslp = slp->data.ptrvalue;
         subslp != NULL;
         subslp = subslp_next) {
      subslp_next = subslp->next;
      subslp->next = subslp_list;
      subslp_list = subslp;
      ReverseLocationIntervalOrder (subslp);
    }
    slp->data.ptrvalue = subslp_list;
  } else if (slp->choice == SEQLOC_PACKED_PNT) {
    pspp = (PackSeqPntPtr) slp->data.ptrvalue;
    if (pspp != NULL) {
      for (pnt_num = 0; pnt_num < pspp->used / 2; pnt_num++) {
        swap = pspp->pnts[pnt_num];
        pspp->pnts[pnt_num] = pspp->pnts[pspp->used - 1 - pnt_num];
        pspp->pnts[pspp->used - 1 - pnt_num] = swap;
      }
    }
  }
}

static void ReverseIntervalsCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)

{
  if (sfp == NULL || sfp->location == NULL) return;

  ReverseLocationIntervalOrder (sfp->location);
}


static Boolean ReverseIntervalsAction (Pointer userdata)
{
  ReverseIntervalsFormPtr   mrfp;
  FilterSetPtr          fsp;
  SeqEntryPtr           sep;
  ValNodePtr            feature_type_list, vnp;
  Uint1                 feat_def_choice;
  
  if (userdata == NULL) return FALSE;
  
  mrfp = (ReverseIntervalsFormPtr) userdata;
  fsp = (FilterSetPtr) DialogToPointer (mrfp->constraints);
  
  sep = GetTopSeqEntryForEntityID (mrfp->input_entityID);
  if (sep == NULL) return FALSE;
  
  feature_type_list = (ValNodePtr) DialogToPointer (mrfp->feature_select);
  if (feature_type_list == NULL)
  {
    FilterSetFree (fsp);
    return FALSE;
  }
    
  for (vnp = feature_type_list; vnp != NULL; vnp = vnp->next)
  {
    feat_def_choice = vnp->choice;
    if (feat_def_choice == 255)
    {
      feat_def_choice = 0;
    }
    OperateOnSeqEntryConstrainedObjects (sep, fsp, 
                                         ReverseIntervalsCallback,
                                         NULL, 0, feat_def_choice, 0, NULL);
  }
  
  ValNodeFree (feature_type_list);
  FilterSetFree (fsp);
  
  ObjMgrSetDirtyFlag (mrfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, mrfp->input_entityID, 0, 0);  
  Update ();
  return TRUE;
}

extern void ReverseFeatureIntervals (IteM i)
{
  BaseFormPtr              bfp;
  ReverseIntervalsFormPtr  fp;
  WindoW                   w;
  GrouP                    h;
  SeqEntryPtr              sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  fp = (ReverseIntervalsFormPtr) MemNew (sizeof (ReverseIntervalsFormData));
  if (fp == NULL) return;
  
  w = FixedWindow (-50, -33, -10, -10, "Reverse Feature Intervals", StdCloseWindowProc);
  SetObjectExtra (w, fp, StdCleanupFormProc);
  fp->form = (ForM) w;
  fp->input_entityID = bfp->input_entityID;
  sep = GetTopSeqEntryForEntityID(bfp->input_entityID);
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  fp->feature_select =  FeatureSelectionDialogEx (h, TRUE, sep,
                                                  ReverseIntervalsChangeNotify, 
                                                  fp);
    
  fp->constraints = FilterGroup (h, TRUE, FALSE, TRUE, FALSE, FALSE, "Where feature text");
  fp->accept_cancel = AcceptCancelDialog (h, ReverseIntervalsAction, NULL, ReverseIntervalsClear, ReverseIntervalsClearText, (Pointer)fp, w);
  AlignObjects (ALIGN_CENTER, (HANDLE) fp->feature_select,
                              (HANDLE) fp->constraints,
                              (HANDLE) fp->accept_cancel, NULL);
                              
  Show (w);
}



extern ParseFieldPtr ParseFieldFree (ParseFieldPtr pfp)
{
  if (pfp != NULL)
  {
    switch (pfp->parse_field_type) {
      case PARSE_FIELD_SOURCE_QUAL :
      case PARSE_FIELD_FEATURE_NOTE:
      case PARSE_FIELD_DBXREF:
        pfp->feature_field = ValNodeFreeData (pfp->feature_field);
        break;
      case PARSE_FIELD_DEFLINE:
      case PARSE_FIELD_BIOSRC_STRING:
      case PARSE_FIELD_GENE_FIELD:
      case PARSE_FIELD_RNA_FIELD:
      case PARSE_FIELD_CDS_COMMENT:
      case PARSE_FIELD_COMMENT_DESC:
      case PARSE_FIELD_IMPORT_QUAL:
      case PARSE_FIELD_PROTEIN_FIELD:
        pfp->feature_field = ValNodeFree (pfp->feature_field);
        break;
      default:
        pfp->feature_field = ValNodeFree (pfp->feature_field);
    }

    if (pfp->parse_field_type == PARSE_FIELD_RNA_FIELD)
    {
      pfp->feature_subtype = ValNodeFree (pfp->feature_subtype);
    }
    else
    {
      pfp->feature_subtype = ValNodeFreeData (pfp->feature_subtype);
    }
    pfp = MemFree (pfp);
  }
  return pfp;
}


#define AECR_APPLY   1
#define AECR_EDIT    2
#define AECR_CONVERT 3
#define AECR_SWAP    4
#define AECR_PARSE   5
#define AECR_REMOVE  6
#define NUM_AECR     6


static void FeatureEditor (IteM i, Int4 first_action)
{
  BaseFormPtr     bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  FeatureEditorBaseForm (bfp, first_action);
}


extern void FeatureEvidenceEditor (IteM i)
{
  FeatureEditor(i, FEAT_ED_EVIDENCE);
}

extern void FeatureExceptionEditor (IteM i)
{
  FeatureEditor(i, FEAT_ED_EXCEPTION);
}

extern void FeaturePartialEditor (IteM i)
{
  FeatureEditor(i, FEAT_ED_PARTIAL);
}

extern void FeatureStrandEditor (IteM i)
{
  FeatureEditor(i, FEAT_ED_STRAND);
}

extern void FeatureCitationEditor (IteM i)
{
  FeatureEditor(i, FEAT_ED_CITATION);
}

extern void FeatureExperimentEditor (IteM i)
{
  FeatureEditor(i, FEAT_ED_EXPERIMENT);
}

extern void FeatureInferenceEditor (IteM i)
{
  FeatureEditor(i, FEAT_ED_INFERENCE);
}

extern void FeaturePseudoEditor (IteM i)
{
  FeatureEditor (i, FEAT_ED_PSEUDO);
}

typedef struct cdset 
{
  ValNodePtr cds_list;
  ValNodePtr gene_list;
  ValNodePtr prot_list;
  ValNodePtr mrna_list;
} CDSetData, PNTR CDSetPtr;

static ValNodePtr FreeCDSetList (ValNodePtr vnp)
{
  CDSetPtr cdsp;
  
  if (vnp != NULL)
  {
    FreeCDSetList (vnp->next);
    vnp->next = NULL;
    cdsp = (CDSetPtr) vnp->data.ptrvalue;
    if (cdsp != NULL)
    {
      cdsp->cds_list = ValNodeFree (cdsp->cds_list);
      cdsp->gene_list = ValNodeFree (cdsp->gene_list);
      cdsp->prot_list = ValNodeFree (cdsp->prot_list);
      cdsp->mrna_list = ValNodeFree (cdsp->mrna_list);
    }
    ValNodeFreeData (vnp);
  }
  return NULL;
}

extern CharPtr GetCDSGeneProtFieldName (ValNodePtr vnp)
{
  CharPtr PNTR      field_name_list = NULL;
  Int4              num_fields = 0;
  CharPtr           label = NULL;

  field_name_list = BuildCDSGeneFieldList (&num_fields);
  if (field_name_list != NULL
      && vnp->data.intvalue > 0 && vnp->data.intvalue <= num_fields)
  {
    label = StringSave (field_name_list [vnp->data.intvalue - 1]);
  }
  FreeCDSGeneFieldList (field_name_list, num_fields);
  return label;
}

static DialoG GetCDSGeneProtSample (GrouP h, Uint2 entityID)
{
  CharPtr PNTR      field_name_list = NULL;
  Int4              num_fields = 0, i;
  DialoG            d = NULL;
  SetSampleData     ssd;
  
  field_name_list = BuildCDSGeneFieldList (&num_fields);
  if (field_name_list != NULL)
  {    
    ssd.entityID = entityID;
    ssd.field_list = NULL;
    for (i = 0; i < num_fields; i++)
    {
      ValNodeAddInt (&ssd.field_list, 0, i + 1);
    }
    ssd.fieldstring_func = GetCDSGeneProtField;
    ssd.descrstring_func = NULL;
    ssd.fsp = NULL;
    ssd.free_vn_proc = NULL;
    ssd.copy_vn_proc = IntValNodeCopy;
    ssd.match_vn_proc = IntValNodeMatch;
    ssd.label_vn_proc = GetCDSGeneProtFieldName;  
    
    
    d = SampleDialog (h);
    PointerToDialog (d, &ssd);
    FreeCDSGeneFieldList (field_name_list, num_fields);
  }
  return d;
}


typedef struct buildcdset2
{
  ValNodePtr cds_list;
  ValNodePtr mrna_list;
  ValNodePtr gene_list;
} BuildCDSet2Data, PNTR BuildCDSet2Ptr;

static void BuildCDSet2Callback (SeqFeatPtr sfp, Pointer userdata)
{
  BuildCDSet2Ptr b;

  if (sfp == NULL || sfp->idx.deleteme || userdata == NULL) return;
  b = (BuildCDSet2Ptr) userdata;
  if (sfp->data.choice == SEQFEAT_CDREGION)
  {
    ValNodeAddPointer (&(b->cds_list), OBJ_SEQFEAT, sfp);
  }
  else if (sfp->data.choice == SEQFEAT_GENE)
  {
    ValNodeAddPointer (&(b->gene_list), OBJ_SEQFEAT, sfp);
  }
  else if (sfp->idx.subtype == FEATDEF_mRNA)
  {
    ValNodeAddPointer (&(b->mrna_list), OBJ_SEQFEAT, sfp);
  }
  else if (SeqMgrGetGeneXref (sfp) != NULL)
  {
    ValNodeAddPointer (&(b->gene_list), OBJ_SEQFEAT, sfp);
  }
}


static void UnmarkFeatureList (ValNodePtr list)
{
  SeqFeatPtr sfp;

  while (list != NULL)
  {
    sfp = list->data.ptrvalue;
    if (sfp != NULL)
    {
      sfp->idx.deleteme = FALSE;
    }
    list = list->next;
  }
}


static SeqFeatPtr GetGeneForCDSetFeature (SeqFeatPtr sfp)
{
  GeneRefPtr        grp;
  SeqFeatPtr        gene = NULL;
  SeqMgrFeatContext fcontext;

  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL)
  {
    gene = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
  }
  else if (!SeqMgrGeneIsSuppressed (grp))
  {
    if (!StringHasNoText (grp->locus_tag))
    {
      gene = SeqMgrGetGeneByLocusTag (BioseqFindFromSeqLoc (sfp->location), grp->locus_tag, &fcontext);
    }
  }
  return gene;
}


static CDSetPtr BuildCDSetFromCodingRegion (SeqFeatPtr cds, BoolPtr indexing_needed)
{
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        gene = NULL, mrna, prot;
  BioseqPtr         protbsp;
  CDSetPtr          cdsp;
  ProtRefPtr        prp;

  if (cds == NULL || cds->data.choice != SEQFEAT_CDREGION) return NULL;

  cdsp = (CDSetPtr) MemNew (sizeof (CDSetData));
  ValNodeAddPointer (&(cdsp->cds_list), 0, cds);

  gene = GetGeneForCDSetFeature (cds);
  if (gene != NULL)
  {
    ValNodeAddPointer (&(cdsp->gene_list), 0, gene);
    /* mark gene, so that we'll know it isn't lonely */
    gene->idx.deleteme = TRUE;
  }

  mrna = SeqMgrGetOverlappingmRNA (cds->location, &fcontext);
  if (mrna != NULL)
  {
    ValNodeAddPointer (&(cdsp->mrna_list), 0, mrna);
    /* mark mrna, so that we'll know it's already in a set */
    mrna->idx.deleteme = TRUE;
  }

  if (cds->product != NULL)
  {
    protbsp = BioseqFindFromSeqLoc (cds->product);
    if (protbsp != NULL)
    {
      prot = SeqMgrGetNextFeature (protbsp, NULL, SEQFEAT_PROT, FEATDEF_PROT, &fcontext);
      /* if there is no full-length protein feature, make one */
      if (prot == NULL)
      {
        prp = ProtRefNew ();
        prot = CreateNewFeatureOnBioseq (protbsp, SEQFEAT_PROT, NULL);
        if (prot != NULL)
        {
          prot->data.value.ptrvalue = prp;
          if (indexing_needed != NULL)
          {
            *indexing_needed = TRUE;
          }
        }
      }
      if (prot != NULL)
      {
        ValNodeAddPointer (&(cdsp->prot_list), 0, prot);
      }
      
      /* also add in mat_peptides from protein feature */
      prot = SeqMgrGetNextFeature (protbsp, NULL, SEQFEAT_PROT, FEATDEF_mat_peptide_aa, &fcontext);
      while (prot != NULL)
      {
        ValNodeAddPointer (&(cdsp->prot_list), 0, prot);
        prot = SeqMgrGetNextFeature (protbsp, prot, SEQFEAT_PROT, FEATDEF_mat_peptide_aa, &fcontext);
      }
    }
  }  
  return cdsp;
}


static CDSetPtr BuildCDSetFrommRNA (SeqFeatPtr mrna)
{
  SeqFeatPtr        gene;
  CDSetPtr          cdsp;

  if (mrna == NULL || mrna->idx.deleteme || mrna->idx.subtype != FEATDEF_mRNA) return NULL;

  cdsp = (CDSetPtr) MemNew (sizeof (CDSetData));
  ValNodeAddPointer (&(cdsp->mrna_list), 0, mrna);

  gene = GetGeneForCDSetFeature (mrna);
  if (gene != NULL)
  {
    ValNodeAddPointer (&(cdsp->gene_list), 0, gene);
    /* mark gene, so that we'll know it isn't lonely */
    gene->idx.deleteme = TRUE;
  }

  return cdsp;
}


static Boolean DoesCDSetMatchConstraint (CDSetPtr cdsp, ChoiceConstraintPtr ccp);

static ValNodePtr BuildCDSet2List (Uint2 entityID, ChoiceConstraintPtr  ccp)
{
  SeqEntryPtr    sep;
  BuildCDSet2Data b;
  CDSetPtr       cdsp;
  ValNodePtr     vnp, vnp_next, vnp_prev;
  ValNodePtr     cdset_list = NULL;
  SeqFeatPtr     cds, gene, mrna;
  Boolean        need_indexing = FALSE;
  
  sep = GetTopSeqEntryForEntityID (entityID);

  b.cds_list = NULL;
  b.gene_list = NULL;
  b.mrna_list = NULL;
  
  VisitFeaturesInSep (sep, &b, BuildCDSet2Callback);

  /* build cdsets that have coding regions */
  for (vnp = b.cds_list; vnp != NULL; vnp = vnp->next)
  {
    cds = (SeqFeatPtr) vnp->data.ptrvalue;
    if (cds == NULL) continue;
    cdsp = BuildCDSetFromCodingRegion (cds, &need_indexing);
    if (cdsp != NULL)
    {
      ValNodeAddPointer (&cdset_list, 0, cdsp);
    }
  }
  if (need_indexing)
  {
    /* indexing because we have created full-length protein features */
    SeqMgrIndexFeatures (entityID, NULL);
  }

  /* build cdsets for mrna features that don't have coding regions */
  for (vnp = b.mrna_list; vnp != NULL; vnp = vnp->next)
  {
    mrna = (SeqFeatPtr) vnp->data.ptrvalue;
    if (mrna == NULL || mrna->idx.deleteme) continue;
    cdsp = BuildCDSetFrommRNA (mrna);
    if (cdsp != NULL)
    {
      ValNodeAddPointer (&cdset_list, 0, cdsp);
    }
  }

  /* build cdsets for lonely genes / features with gene xrefs that are not coding regions or mrnas */
  for (vnp = b.gene_list; vnp != NULL; vnp = vnp->next)
  {
    gene = (SeqFeatPtr) vnp->data.ptrvalue;
    if (gene == NULL || gene->idx.deleteme) continue;
    cdsp = (CDSetPtr) MemNew (sizeof (CDSetData));
    ValNodeAddPointer (&(cdsp->gene_list), 0, gene);
    ValNodeAddPointer (&cdset_list, 0, cdsp);
  }

  /* now unmark features */
  UnmarkFeatureList (b.cds_list);
  UnmarkFeatureList (b.mrna_list);
  UnmarkFeatureList (b.gene_list);

  b.cds_list = ValNodeFree (b.cds_list);
  b.mrna_list = ValNodeFree (b.mrna_list);
  b.gene_list = ValNodeFree (b.gene_list);

  /* now remove sets that don't match our choice constraint */
  vnp_prev = NULL;
  for (vnp = cdset_list; vnp != NULL; vnp = vnp_next)
  {
    vnp_next = vnp->next;
    if (!DoesCDSetMatchConstraint (vnp->data.ptrvalue, ccp))
    {
      if (vnp_prev == NULL)
      {
        cdset_list = vnp->next;
      }
      else
      {
        vnp_prev->next = vnp->next;
      }
      vnp->next = NULL;
      FreeCDSetList (vnp);     
    }
    else
    {
      vnp_prev = vnp;
    }
  }
  
  return cdset_list;
}

typedef struct buildcdset
{
  ValNodePtr cdset_list;
  ValNodePtr lonely_gene_list;
} BuildCDSetData, PNTR BuildCDSetPtr;


static CharPtr GetCDSetField (CDSetPtr cdsp, ValNodePtr vnp)
{
  ValNodePtr sfp_vnp;
  CharPtr    str = NULL;
  
  if (cdsp == NULL || vnp == NULL)
  {
    return NULL;
  }
  
  for (sfp_vnp = cdsp->cds_list; sfp_vnp != NULL && str == NULL; sfp_vnp = sfp_vnp->next)
  {
    str = GetCDSGeneProtField (sfp_vnp->data.ptrvalue, vnp, NULL);
  }
  for (sfp_vnp = cdsp->gene_list; sfp_vnp != NULL && str == NULL; sfp_vnp = sfp_vnp->next)
  {
    str = GetCDSGeneProtField (sfp_vnp->data.ptrvalue, vnp, NULL);
  }
  for (sfp_vnp = cdsp->prot_list; sfp_vnp != NULL && str == NULL; sfp_vnp = sfp_vnp->next)
  {
    str = GetCDSGeneProtField (sfp_vnp->data.ptrvalue, vnp, NULL);
  }
  for (sfp_vnp = cdsp->mrna_list; sfp_vnp != NULL && str == NULL; sfp_vnp = sfp_vnp->next)
  {
    str = GetCDSGeneProtField (sfp_vnp->data.ptrvalue, vnp, NULL);
  }
  return str;
}

static Boolean SetCDSetField (CDSetPtr cdsp, ValNodePtr vnp, ApplyValuePtr avp, FilterSetPtr fsp)
{
  ValNodePtr sfp_vnp;
  Boolean    rval = FALSE;
  
  if (cdsp == NULL || vnp == NULL || avp == NULL)
  {
    return FALSE;
  }
  
  for (sfp_vnp = cdsp->cds_list; sfp_vnp != NULL ; sfp_vnp = sfp_vnp->next)
  {
    rval |= SetCDSGeneProtField (sfp_vnp->data.ptrvalue, vnp, avp, fsp);
  }
  for (sfp_vnp = cdsp->gene_list; sfp_vnp != NULL; sfp_vnp = sfp_vnp->next)
  {
    if (sfp_vnp->choice == 0) /* don't set again if already set */
    {
      rval |= SetCDSGeneProtField (sfp_vnp->data.ptrvalue, vnp, avp, fsp);
    }
  }
  for (sfp_vnp = cdsp->prot_list; sfp_vnp != NULL; sfp_vnp = sfp_vnp->next)
  {
    rval |= SetCDSGeneProtField (sfp_vnp->data.ptrvalue, vnp, avp, fsp);
  }
  for (sfp_vnp = cdsp->mrna_list; sfp_vnp != NULL; sfp_vnp = sfp_vnp->next)
  {
    rval |= SetCDSGeneProtField (sfp_vnp->data.ptrvalue, vnp, avp, fsp);
  }
  return rval;
}

static void RemoveCDSetField (CDSetPtr cdsp, ValNodePtr vnp, FilterSetPtr fsp)
{
  ValNodePtr sfp_vnp;
  CharPtr    str = NULL;
  
  if (cdsp == NULL || vnp == NULL)
  {
    return;
  }
  
  for (sfp_vnp = cdsp->cds_list; sfp_vnp != NULL && str == NULL; sfp_vnp = sfp_vnp->next)
  {
    RemoveCDSGeneProtField (sfp_vnp->data.ptrvalue, vnp, fsp);
  }
  for (sfp_vnp = cdsp->gene_list; sfp_vnp != NULL && str == NULL; sfp_vnp = sfp_vnp->next)
  {
    RemoveCDSGeneProtField (sfp_vnp->data.ptrvalue, vnp, fsp);
  }
  for (sfp_vnp = cdsp->prot_list; sfp_vnp != NULL && str == NULL; sfp_vnp = sfp_vnp->next)
  {
    RemoveCDSGeneProtField (sfp_vnp->data.ptrvalue, vnp, fsp);
  }
  for (sfp_vnp = cdsp->mrna_list; sfp_vnp != NULL && str == NULL; sfp_vnp = sfp_vnp->next)
  {
    RemoveCDSGeneProtField (sfp_vnp->data.ptrvalue, vnp, fsp);
  }
}

static Boolean GetCDSetPseudoMatch (CDSetPtr cdsp, Int4 featdef_type, Boolean is_pseudo)
{
  ValNodePtr sfp_vnp;
  CharPtr    str = NULL;
  SeqFeatPtr sfp;
  GeneRefPtr grp;
  Boolean    gene_is_pseudo;
  
  if (cdsp == NULL)
  {
    return FALSE;
  }
  
  if (featdef_type == FEATDEF_CDS || featdef_type == FEATDEF_ANY)
  {
    for (sfp_vnp = cdsp->cds_list; sfp_vnp != NULL && str == NULL; sfp_vnp = sfp_vnp->next)
    {
      sfp = sfp_vnp->data.ptrvalue;
      if (sfp != NULL && ((is_pseudo && sfp->pseudo) || (!is_pseudo && ! sfp->pseudo)))
      {
        return TRUE;
      }
    }
  }
  if (featdef_type == FEATDEF_GENE || featdef_type == FEATDEF_ANY)
  {
    for (sfp_vnp = cdsp->gene_list; sfp_vnp != NULL && str == NULL; sfp_vnp = sfp_vnp->next)
    {
      sfp = sfp_vnp->data.ptrvalue;
      if (sfp != NULL)
      {
        gene_is_pseudo = FALSE;
        if (sfp->pseudo)
        {
          gene_is_pseudo = TRUE;
        }
        if (sfp->data.choice == SEQFEAT_GENE)
        {
          grp = (GeneRefPtr) sfp->data.value.ptrvalue;
          if (grp != NULL && grp->pseudo)
          {
            gene_is_pseudo = TRUE;
          }
        }
        if ((is_pseudo && gene_is_pseudo) || (!is_pseudo && ! gene_is_pseudo))
        {
          return TRUE;
        }
      }
    }
  }
  if (featdef_type == FEATDEF_mRNA || featdef_type == FEATDEF_ANY)
  {
    for (sfp_vnp = cdsp->mrna_list; sfp_vnp != NULL && str == NULL; sfp_vnp = sfp_vnp->next)
    {
      sfp = sfp_vnp->data.ptrvalue;
      if (sfp != NULL && ((is_pseudo && sfp->pseudo) || (!is_pseudo && ! sfp->pseudo)))
      {
        return TRUE;
      }
    }
  }
  return FALSE;
}

static Boolean DoCDSFieldsMatch (CDSetPtr cdsp, ValNodePtr qual_choice, ValNodePtr qual_choice_match)
{
  CharPtr str1, str2;
  Int4    field_choice, field_choice2;
  ValNode field_vn, field_vn2;
  Boolean does_match = FALSE;

  if (qual_choice == NULL && qual_choice_match == NULL) 
  {
    /* both are wildcards */
    field_vn.next = NULL;
    field_vn.choice = 0;
    field_vn2.next = NULL;
    field_vn2.choice = 0;

    for (field_choice = 1;
         ! does_match && field_choice < num_gene_fields + num_mrna_fields + num_protein_fields;
         field_choice++)
    {
      field_vn.data.intvalue = field_choice;
      str1 = GetCDSetField (cdsp, &field_vn);
      if (!StringHasNoText (str1)) {
        for (field_choice2 = field_choice + 1;
             ! does_match && field_choice2 < 1 + num_gene_fields + num_mrna_fields + num_protein_fields;
             field_choice2++)
        {
          field_vn2.data.intvalue = field_choice2;
          str2 = GetCDSetField (cdsp, &field_vn2);
          if (StringCmp (str1, str2) == 0) 
          {
            does_match = TRUE;
          }
          str2 = MemFree (str2);
        }
      }
      str1 = MemFree (str1);
    }
  }
  else if (qual_choice == NULL || qual_choice_match == NULL) 
  {
    /* one is a wildcard */
    /* make sure it is the first one */
    if (qual_choice_match == NULL) 
    {
      qual_choice_match = qual_choice;
      qual_choice = NULL;
    }
    field_vn.next = NULL;
    field_vn.choice = 0;
    str2 = GetCDSetField (cdsp, qual_choice_match);
    if (!StringHasNoText (str2)) {
      for (field_choice = 1;
           ! does_match && field_choice < num_gene_fields + num_mrna_fields + num_protein_fields;
           field_choice++)
      {
        if (field_choice != qual_choice_match->data.intvalue) {
          /* don't compare something to itself */
          field_vn.data.intvalue = field_choice;
          str1 = GetCDSetField (cdsp, &field_vn);
          if (StringCmp (str1, str2) == 0) 
          {
            does_match = TRUE;
          }
          str1 = MemFree (str1);
        }
      }
    }
    str2 = MemFree (str2);
  } else {
    str1 = GetCDSetField (cdsp, qual_choice);
    str2 = GetCDSetField (cdsp, qual_choice_match);
    if (!StringHasNoText (str1) && StringCmp (str1, str2) == 0)
    {
      does_match = TRUE;
    }
    str1 = MemFree (str1);
    str2 = MemFree (str2);
  }
  return does_match;
}


static Boolean DoesCDSetMatchConstraint (CDSetPtr cdsp, ChoiceConstraintPtr ccp)
{
  Boolean                does_match = FALSE;
  CharPtr                str;
  Int4                   field_choice;
  ValNode                field_vn;
  ValNodePtr             sfp_vnp;
  
  if (cdsp == NULL)
  {
    return FALSE;
  }
  if (ccp == NULL || ccp->constraint_type == CHOICE_CONSTRAINT_ANY)
  {
    does_match = TRUE;
  }
  else if (ccp->constraint_type == CHOICE_CONSTRAINT_QUAL_PRESENT)
  {
    if (ccp->qual_choice == NULL)
    {
      does_match = TRUE;
    }
    else
    {
      str = GetCDSetField (cdsp, ccp->qual_choice);
      if (str != NULL)
      {
        MemFree (str);
        does_match = TRUE;
      }
    }
  }
  else if (ccp->constraint_type == CHOICE_CONSTRAINT_STRING)
  {
    if (ccp->qual_choice == NULL) 
    {
      field_vn.next = NULL;
      field_vn.choice = 0;
      for (field_choice = 1;
           ! does_match && field_choice < 1 + num_gene_fields + num_mrna_fields + num_protein_fields;
           field_choice++)
      {
        field_vn.data.intvalue = field_choice;
        str = GetCDSetField (cdsp, &field_vn);
        does_match = DoesStringMatchConstraintX (str, ccp->string_constraint);
        MemFree (str);
      }
      if (ccp->string_constraint != NULL && ccp->string_constraint->not_present)
      {
        does_match = ! does_match;
      }
    }
    else if (IsCDSetMatPeptideQualChoice(ccp->qual_choice))
    {
      /* if we're matching protein fields, we need to look at all of them */
      for (sfp_vnp = cdsp->prot_list; sfp_vnp != NULL && !does_match; sfp_vnp = sfp_vnp->next)
      {
        str = GetCDSGeneProtField (sfp_vnp->data.ptrvalue, ccp->qual_choice, NULL);
        does_match = DoesStringMatchConstraintX (str, ccp->string_constraint);
        if (ccp->string_constraint != NULL && ccp->string_constraint->not_present)
        {
          does_match = ! does_match;
        }
        MemFree (str);
      }
    } 
    else
    {
      str = GetCDSetField (cdsp, ccp->qual_choice);
      does_match = DoesStringMatchConstraintX (str, ccp->string_constraint);
      MemFree (str);
      if (ccp->string_constraint != NULL && ccp->string_constraint->not_present)
      {
        does_match = ! does_match;
      }
    }
  }
  else if (ccp->constraint_type == CHOICE_CONSTRAINT_MATCH)
  {
    does_match = DoCDSFieldsMatch (cdsp, ccp->qual_choice, ccp->qual_choice_match);
  }
  else if (ccp->constraint_type == CHOICE_CONSTRAINT_PSEUDO)
  {
    if (ccp->pseudo_constraint == NULL)
    {
      does_match = TRUE;
    }
    else
    {
      does_match = GetCDSetPseudoMatch (cdsp, 
                                        ccp->pseudo_constraint->featdef_type, 
                                        ccp->pseudo_constraint->is_pseudo);
    }
  }
  return does_match;
}

static SeqLocPtr GetNewFeatureLocation(ValNodePtr feat_list, Uint2 choice, Uint2 choice_from)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  BioseqPtr  bsp = NULL;
  SeqLocPtr  slp = NULL;
  
  for (vnp = feat_list; vnp != NULL && slp == NULL; vnp = vnp->next)
  {
    sfp = vnp->data.ptrvalue;
    if (sfp != NULL)
    {
      if (choice == FEATDEF_PROT || choice == FEATDEF_mat_peptide_aa)
      {
        bsp = BioseqFindFromSeqLoc (sfp->product);
        slp = SeqLocWholeNew (bsp);
      }
      else
      {
        if (choice_from == sfp->idx.subtype || (sfp->idx.subtype == FEATDEF_CDS && (choice_from == FEATDEF_PROT || choice_from == FEATDEF_mat_peptide_aa)))
        {
          slp = (SeqLocPtr) AsnIoMemCopy (sfp->location, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
        }
        else
        { 
          bsp = BioseqFindFromSeqLoc (sfp->location);
          slp = SeqLocWholeNew (bsp);
        }
      }
    }
  }  
  return slp;  
}


static SeqFeatPtr BuildOneNewCDSGeneProtFeature (CDSetPtr cdsp, Uint2 choice, Uint2 choice_from)
{
  BioseqPtr bsp = NULL;
  SeqLocPtr slp = NULL, slp_tmp, slp_sum;
  SeqFeatPtr sfp = NULL;
  RnaRefPtr  rrp;
  ProtRefPtr prp;
  
  if (cdsp == NULL)
  {
    return NULL;
  }

  if (choice == FEATDEF_PROT || choice == FEATDEF_mat_peptide_aa)
  {
    slp = GetNewFeatureLocation(cdsp->cds_list, choice, choice_from);
  }
  else
  {
    /* get location from CDS features */
    if (choice_from == 0 || choice_from == FEATDEF_CDS || choice_from == FEATDEF_PROT || choice_from == FEATDEF_mat_peptide_aa)
    {
      slp = GetNewFeatureLocation(cdsp->cds_list, choice, choice_from);
    }

    /* get location from mRNA features */
    if ((choice_from == 0 && slp == NULL) || choice_from == FEATDEF_mRNA)
    {    
      slp_tmp = GetNewFeatureLocation(cdsp->mrna_list, choice, choice_from);
      if (slp == NULL) 
      {
        slp = slp_tmp;
      }
      else
      {
        bsp = BioseqFindFromSeqLoc (sfp->location);
        if (bsp != NULL)
        {
          slp_sum = SeqLocMerge (bsp, slp_tmp, slp, TRUE, FALSE, FALSE);
          slp = SeqLocFree (slp);
          slp_tmp = SeqLocFree (slp_tmp);
          slp = slp_sum;
          slp_sum = NULL;
        }
      }
    }

    /* get location from gene features */
    if ((choice_from == 0 && slp == NULL) || choice_from == FEATDEF_GENE)
    {    
      slp_tmp = GetNewFeatureLocation(cdsp->gene_list, choice, choice_from);
      if (slp == NULL) 
      {
        slp = slp_tmp;
      }
      else
      {
        bsp = BioseqFindFromSeqLoc (sfp->location);
        if (bsp != NULL)
        {
          slp_sum = SeqLocMerge (bsp, slp_tmp, slp, TRUE, FALSE, FALSE);
          slp = SeqLocFree (slp);
          slp_tmp = SeqLocFree (slp_tmp);
          slp = slp_sum;
          slp_sum = NULL;
        }
      }
    }
  }  
  
  bsp = BioseqFindFromSeqLoc (slp);
  if (bsp != NULL && slp != NULL)
  {
    sfp = CreateNewFeatureOnBioseq (bsp, FindFeatFromFeatDefType(choice), slp);
    switch (choice)
    {
      case FEATDEF_mRNA:
        rrp = RnaRefNew ();
        rrp->type = 2;
        sfp->data.value.ptrvalue = rrp;
        break;
      case FEATDEF_GENE:
        sfp->data.value.ptrvalue = GeneRefNew();
        break;
      case FEATDEF_CDS:
        sfp->data.value.ptrvalue = CdRegionNew();
        break;
      case FEATDEF_PROT:
      case FEATDEF_mat_peptide_aa:
        prp = ProtRefNew();        
        if (choice == FEATDEF_mat_peptide_aa)
        {
          prp->processed = 2;
        }
        sfp->data.value.ptrvalue = prp;
        break;
    }
  }
  return sfp;
}

static SeqFeatPtr BuildNewCDSGeneProtFeature (CDSetPtr cdsp, ValNodePtr vnp, ValNodePtr field_from)
{
  SeqFeatPtr new_sfp = NULL, sfp;
  Uint2      choice_to, choice_from;
  ValNodePtr v;
  Boolean    found;
  
  if (cdsp == NULL || vnp == NULL)
  {
    return NULL;
  }

  choice_to = FeatDefTypeFromFieldList (vnp);
  choice_from = FeatDefTypeFromFieldList (field_from);

  switch (choice_to)
  {
    case FEATDEF_CDS:
      if (choice_from != 0 && cdsp->cds_list == NULL)
      {
        /* add new CDS */
        new_sfp = BuildOneNewCDSGeneProtFeature (cdsp, choice_to, choice_from);
        if (new_sfp != NULL)
        {
          ValNodeAddPointer (&(cdsp->cds_list), 0, new_sfp);
        }
      }
      break;
    case FEATDEF_GENE:
      if ((choice_from != 0 || vnp->data.intvalue - 1 == GENEFIELD_LOCUS) && cdsp->gene_list == NULL)
      {
        /* add new gene */
        new_sfp = BuildOneNewCDSGeneProtFeature (cdsp, choice_to, choice_from);
        if (new_sfp != NULL)
        {
          ValNodeAddPointer (&(cdsp->gene_list), 0, new_sfp);
        }
      }
      break;
    case FEATDEF_mRNA:
      if (choice_from != 0 && cdsp->mrna_list == NULL)
      {
        /* add new mrna */
        new_sfp = BuildOneNewCDSGeneProtFeature (cdsp, choice_to, choice_from);
        if (new_sfp != NULL)
        {
          ValNodeAddPointer (&(cdsp->mrna_list), 0, new_sfp);
        }
      }
      break;
    case FEATDEF_PROT:
    case FEATDEF_mat_peptide_aa:
      if (choice_from != 0 && cdsp->cds_list != NULL)
      {
        found = FALSE;
        for (v = cdsp->prot_list; v != NULL && !found; v = v->next)
        {
          sfp = (SeqFeatPtr) v->data.ptrvalue;
          if (sfp != NULL && sfp->idx.subtype == choice_to)
          {
            found = FALSE;
          }
        }
        if (!found)
        {
          new_sfp = BuildOneNewCDSGeneProtFeature (cdsp, choice_to, choice_from);
          if (new_sfp != NULL)
          {
            /* set subtype so other functions will work */
            new_sfp->idx.subtype = choice_to;
            ValNodeAddPointer (&(cdsp->mrna_list), 0, new_sfp);
          }
        }
      }
      break;
  }
  return new_sfp;        
}

static Boolean CDSetOnBioseq (CDSetPtr cdsp, BioseqPtr bsp)
{
  SeqFeatPtr sfp;
  Boolean    is_on_bioseq = FALSE;
  ValNodePtr sfp_vnp;

  if (cdsp == NULL || bsp == NULL) return FALSE;

  for (sfp_vnp = cdsp->cds_list; sfp_vnp != NULL && !is_on_bioseq; sfp_vnp = sfp_vnp->next)
  {
    sfp = sfp_vnp->data.ptrvalue;
    if (sfp != NULL && BioseqFindFromSeqLoc (sfp->location) == bsp) {
      is_on_bioseq = TRUE;
    }
  }
  for (sfp_vnp = cdsp->gene_list; sfp_vnp != NULL && !is_on_bioseq; sfp_vnp = sfp_vnp->next)
  {
    sfp = sfp_vnp->data.ptrvalue;
    if (sfp != NULL && BioseqFindFromSeqLoc (sfp->location) == bsp) {
      is_on_bioseq = TRUE;
    }
  }
  for (sfp_vnp = cdsp->mrna_list; sfp_vnp != NULL && !is_on_bioseq; sfp_vnp = sfp_vnp->next)
  {
    sfp = sfp_vnp->data.ptrvalue;
    if (sfp != NULL && BioseqFindFromSeqLoc (sfp->location) == bsp) {
      is_on_bioseq = TRUE;
    }
  }
  return is_on_bioseq;
}

static Boolean ValNodeDataAlreadyInValNodeList (ValNodePtr list, ValNodePtr findme)
{
  if (findme == NULL) return FALSE;
  while (list != NULL) 
  {
    if (list->data.ptrvalue == findme->data.ptrvalue)
    {
      return TRUE;
    }
    list = list->next;
  }
  return FALSE;
}

static void BuildGeneForSegSet (BioseqSetPtr segset, ValNodePtr PNTR cdset_list, Boolean lonely_ok)
{
  SeqEntryPtr       this_sep, master_sep = NULL;
  BioseqSetPtr      part_set = NULL;
  SeqLocPtr         slp = NULL, slp_new, slp_last = NULL;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext context;
  ValNodePtr        vnp, add_gene_to = NULL;
  Boolean           on_part;
  CDSetPtr          cdsp;
  
  if (segset == NULL || segset->seq_set == NULL)
  {
    return;
  }
  
  this_sep = segset->seq_set;
  while (this_sep != NULL && (part_set == NULL || master_sep == NULL))
  {
    if (IS_Bioseq_set (this_sep))
    {
      part_set = (BioseqSetPtr) this_sep->data.ptrvalue;
      if (part_set != NULL && part_set->_class != BioseqseqSet_class_parts)
      {
        part_set = NULL;
      }
    }
    else if (IS_Bioseq (this_sep))
    {
      master_sep = this_sep;
    }
    this_sep = this_sep->next;
  }
  if (part_set != NULL && master_sep != NULL)
  {
    this_sep = part_set->seq_set;
    sfp = SeqMgrGetNextFeature (master_sep->data.ptrvalue, NULL, 
                                SEQFEAT_GENE, FEATDEF_GENE, &context);
    
    while (this_sep != NULL && sfp == NULL)
    {
      if (IS_Bioseq (this_sep))
      {
        sfp = SeqMgrGetNextFeature (this_sep->data.ptrvalue, NULL, 
                                    SEQFEAT_GENE, FEATDEF_GENE, &context);
      }
      this_sep = this_sep->next;
    }
    if (sfp != NULL)
    {
      return;
    }

    /* no other genes anywhere in segmented set.  Any cdsets? */
    for (vnp = *cdset_list; vnp != NULL; vnp = vnp->next)
    {
      cdsp = (CDSetPtr) vnp->data.ptrvalue;
      if (CDSetOnBioseq(cdsp, master_sep->data.ptrvalue)
          && !ValNodeDataAlreadyInValNodeList(add_gene_to, vnp)) 
      {
        ValNodeAddPointer (&add_gene_to, 0, cdsp);
      }
      else
      {
        on_part = FALSE;
        for (this_sep = part_set->seq_set; this_sep != NULL && !on_part; this_sep = this_sep->next) 
        {
          if (CDSetOnBioseq(cdsp, this_sep->data.ptrvalue)
              && !ValNodeDataAlreadyInValNodeList(add_gene_to, vnp)) 
          {
            ValNodeAddPointer (&add_gene_to, 0, cdsp);
            on_part = TRUE;
          }
        }
      }
    }

    if (add_gene_to != NULL || lonely_ok) {    
      /* create gene */
      this_sep = part_set->seq_set;
      while (this_sep != NULL)
      {
        if (IS_Bioseq (this_sep))
        {
          slp_new = SeqLocWholeNew (this_sep->data.ptrvalue);
          if (slp_new != NULL)
          {
            if (slp_last == NULL)
            {
              slp = slp_new;
            }
            else
            {
              slp_last = ValNodeNew (slp);
              slp_last->choice = SEQLOC_NULL;
              slp_last->next = slp_new;
            }
            slp_last = slp_new;
          }
        }
        this_sep = this_sep->next;
      }
      slp_new = ValNodeNew (NULL);
      slp_new->choice = SEQLOC_MIX;
      slp_new->data.ptrvalue = slp;
      sfp = CreateNewFeature (master_sep, NULL, SEQFEAT_GENE, NULL);
      sfp->location = SeqLocFree (sfp->location);
      sfp->location = slp_new;
      sfp->data.value.ptrvalue = GeneRefNew ();

      /* add to sets */
      if (add_gene_to == NULL)
      {
        cdsp = (CDSetPtr) MemNew (sizeof (CDSetData));
        MemSet (cdsp, 0, sizeof (CDSetData));
        ValNodeAddPointer (&(cdsp->gene_list), 0, sfp);
        ValNodeAddPointer (cdset_list, 0, cdsp);
      }
      else
      {
        for (vnp = add_gene_to; vnp != NULL; vnp = vnp->next)
        {
          cdsp = (CDSetPtr) vnp->data.ptrvalue;
          ValNodeAddPointer (&(cdsp->gene_list), 0, sfp);
        }
      }
    }
    add_gene_to = ValNodeFree (add_gene_to);
  }
}

static Boolean LoneGeneOkWithFilter (FilterSetPtr fsp)
{
  if (fsp == NULL || fsp->cgp == NULL
      || fsp->cgp->constraint_type == CHOICE_CONSTRAINT_ANY
      || (fsp->cgp->constraint_type == CHOICE_CONSTRAINT_QUAL_PRESENT
          && fsp->cgp->qual_choice == NULL))
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

static void BuildNewGenes (SeqEntryPtr sep, ValNodePtr PNTR cdset_list, Boolean lonely_ok)
{
  BioseqSetPtr      bssp = NULL;
  BioseqPtr         bsp = NULL;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext context;
  SeqLocPtr         slp = NULL;
  ValNodePtr        vnp;
  CDSetPtr          cdsp;
  ValNodePtr        add_gene_to;
  
  if (sep == NULL) 
  {
    return;
  }
  
  if (IS_Bioseq (sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (ISA_na (bsp->mol))
    {
      sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, FEATDEF_GENE, &context);
      if (sfp == NULL)
      {
        /* are there any sets already on this Bioseq? */
        add_gene_to = NULL;
        for (vnp = *cdset_list; vnp != NULL; vnp = vnp->next) 
        {
          cdsp = (CDSetPtr) vnp->data.ptrvalue;
          if (CDSetOnBioseq(cdsp, bsp))
          {
            if (cdsp->gene_list == NULL)
            {
              ValNodeAddPointer (&add_gene_to, 0, cdsp);
            }
          }
        }
        if (add_gene_to != NULL || lonely_ok) 
        {
          slp = SeqLocWholeNew (bsp);
          sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, slp);    
          sfp->data.value.ptrvalue = GeneRefNew();  
          if (add_gene_to == NULL)
          {
            cdsp = (CDSetPtr) MemNew (sizeof (CDSetData));
            MemSet (cdsp, 0, sizeof (CDSetData));
            ValNodeAddPointer (&(cdsp->gene_list), 0, sfp);
            ValNodeAddPointer (cdset_list, 0, cdsp);
          }
          else
          {
            for (vnp = add_gene_to; vnp != NULL; vnp = vnp->next)
            {
              cdsp = (CDSetPtr) vnp->data.ptrvalue;
              ValNodeAddPointer (&(cdsp->gene_list), 0, sfp);
            }
          }
        }
        add_gene_to = ValNodeFree (add_gene_to);
      }
    }
  }
  else if (IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp->_class == BioseqseqSet_class_segset)
    {
      BuildGeneForSegSet (bssp, cdset_list, lonely_ok);
    }
    else
    {
      BuildNewGenes (bssp->seq_set, cdset_list, lonely_ok);
    }
  }
  BuildNewGenes (sep->next, cdset_list, lonely_ok);
}

typedef Boolean (*Nlm_TestAECRChoiceProc) PROTO((Pointer));
typedef DialoG  (*Nlm_AECRMakeDlgProc) PROTO ((GrouP, Int4, Nlm_ChangeNotifyProc, Pointer, Uint2));
typedef DialoG  (*Nlm_SampleDialogProc) PROTO ((GrouP, Uint2));
typedef DialoG  (*Nlm_FieldListDlgProc) PROTO ((GrouP, Boolean, Nlm_ChangeNotifyProc, Pointer));
typedef DialoG  (*Nlm_SubtypeListDlgProc) PROTO ((GrouP, Boolean, SeqEntryPtr, Nlm_ChangeNotifyProc, Pointer));

#define AECR_VIB_MSG_SET_DEFAULT (NUM_VIB_MSG + 1)
#define AECR_VIB_MSG_CLEAR_TEXT  (NUM_VIB_MSG + 2)
#define AECR_VIB_MSG_AUTOPOPULATE (NUM_VIB_MSG + 3)

typedef void (*Nlm_ChangeAECRActionProc) PROTO ((DialoG, DialoG));

typedef struct applyeditconvertremove
{
  FORM_MESSAGE_BLOCK
  DialoG                   aecr_pages [NUM_AECR];
  DialoG                   accept_cancel;
  PopuP                    action_popup;
  DialoG                   constraints;
  ButtoN                   clear_constraints_on_action_change;
  Nlm_ChangeAECRActionProc change_action;
  Int4                     prev_page;
  Boolean                  crippled;
  Int4                     crippled_action;
} ApplyEditConvertRemoveData, PNTR ApplyEditConvertRemovePtr;

static void ApplyEditConvertRemoveChangeNotify (Pointer userdata)
{
  ApplyEditConvertRemovePtr mp;
  Int4                      action_choice;
  ValNodePtr                err_list;
  
  mp = (ApplyEditConvertRemovePtr) userdata;
  if (mp == NULL)
  {
    DisableAcceptCancelDialogAccept (mp->accept_cancel);
    return;
  }
  
  if (mp->crippled)
  {
    action_choice = mp->crippled_action;
  }
  else
  {
    action_choice = GetValue (mp->action_popup);
  }
  
  if (action_choice > 0 && action_choice <= NUM_AECR)
  {
    err_list = TestDialog (mp->aecr_pages [action_choice - 1]);
    if (err_list != NULL)
    {
      err_list = ValNodeFree (err_list);
      DisableAcceptCancelDialogAccept (mp->accept_cancel);
    }
    else
    {
      EnableAcceptCancelDialogAccept (mp->accept_cancel);
    }
  }
  else
  {
    DisableAcceptCancelDialogAccept (mp->accept_cancel);
    return;
  }  
}


static void ApplyEditConvertRemoveClearText (Pointer userdata)
{
  ApplyEditConvertRemovePtr mp;
  FilterSetPtr              fsp;
  Int4                      i;

  mp = (ApplyEditConvertRemovePtr) userdata;
  if (mp != NULL)
  {
    for (i = 0; i < NUM_AECR; i++)
    {
      SendMessageToDialog (mp->aecr_pages[i], AECR_VIB_MSG_CLEAR_TEXT);
    }
    fsp = DialogToPointer (mp->constraints);
    FilterSetClearText (fsp);
    PointerToDialog (mp->constraints, fsp);
    FilterSetFree (fsp);
    ApplyEditConvertRemoveChangeNotify (mp);
  }
}

static void ApplyEditConvertRemoveClear (Pointer userdata)
{
  ApplyEditConvertRemovePtr mp;
  Int4                      i;

  mp = (ApplyEditConvertRemovePtr) userdata;
  if (mp != NULL)
  {  
    for (i = 0; i < NUM_AECR; i++)
    {
      PointerToDialog (mp->aecr_pages[i], NULL);
    }
    PointerToDialog (mp->constraints, NULL);    
    ApplyEditConvertRemoveChangeNotify (mp);
  }
}

static void 
ChangeApplyEditConvertRemoveAction 
(ApplyEditConvertRemovePtr mp,
 Int4                      action_choice)
{
  Int4 i;

  if (mp == NULL || action_choice < AECR_APPLY || action_choice > NUM_AECR)
  {
    return;
  }
    
  for (i = 0; i < NUM_AECR; i++)
  {
    if (i + 1 == action_choice)
    {
      Show (mp->aecr_pages[i]);
    }
    else
    {
      Hide (mp->aecr_pages[i]);
    }
  }
  
  if (GetStatus (mp->clear_constraints_on_action_change))
  {
    PointerToDialog (mp->constraints, NULL);
    PointerToDialog (mp->aecr_pages [action_choice - 1], NULL);
  }
  
  if (mp->change_action != NULL)
  {
    if (mp->prev_page < AECR_APPLY || mp->prev_page > NUM_AECR)
    {
      (mp->change_action) (NULL, mp->aecr_pages [action_choice - 1]);
    }
    else
    {
      (mp->change_action) (mp->aecr_pages [mp->prev_page - 1],
                           mp->aecr_pages [action_choice - 1]);
    }
  }
  
  mp->prev_page = action_choice;
  
  if (action_choice == AECR_APPLY)
  {
    SendMessageToDialog (mp->aecr_pages[AECR_APPLY - 1], AECR_VIB_MSG_AUTOPOPULATE);
  }
  else if (action_choice == AECR_EDIT)
  {
    SendMessageToDialog (mp->aecr_pages[AECR_EDIT - 1], AECR_VIB_MSG_AUTOPOPULATE);
  }

  ApplyEditConvertRemoveChangeNotify (mp);
}

static void ChangeApplyEditConvertRemoveActionPopup (PopuP p)
{
  ApplyEditConvertRemovePtr mp;
  Int4                      action_choice;

  mp = (ApplyEditConvertRemovePtr) GetObjectExtra (p);
  if (mp == NULL)
  {
    return;
  }
  
  action_choice = GetValue (mp->action_popup);
  ChangeApplyEditConvertRemoveAction (mp, action_choice);
}
typedef struct samplebtn
{
  Uint2                entityID;
  Nlm_SampleDialogProc make_sample_dlg;
  CharPtr              title;
} SampleBtnData, PNTR SampleBtnPtr;

static void AECRSampleButton (ButtoN b)
{
  SampleBtnPtr sbp;
  WindoW       w;
  
  sbp = GetObjectExtra (b);
  if (sbp != NULL && sbp->make_sample_dlg != NULL)
  {
    w = FixedWindow (-50, -33, -10, -10, sbp->title, StdCloseWindowProc); 
    (sbp->make_sample_dlg) (w, sbp->entityID);
    Show (w);
  }
}

static void 
ApplyEditConvertRemoveCombo 
(IteM                     i,
 Int4                     first_action_choice,
 Boolean                  crippled,
 CharPtr                  title,
 Nlm_AECRMakeDlgProc      make_dlg,
 Boolean                  show_sample,
 Nlm_SampleDialogProc     make_sample_dlg,
 Boolean                  string_constraint,
 Boolean                  source_constraint,
 Boolean                  location_constraint,
 Boolean                  cds_gene_prot_constraint,
 Boolean                  id_list_constraint,
 CharPtr                  string_constraint_label,
 Nlm_AcceptActnProc       accept_actn,
 Nlm_CancelActnProc       cancel_actn,
 Nlm_ChangeAECRActionProc change_action,
 OkToPreSample            presample_check_proc)
{
  BaseFormPtr               bfp;
  ApplyEditConvertRemovePtr mp;
  WindoW                    w;
  GrouP                     h, g1, k, n, m;
  Int4                      page_num;
  SampleBtnPtr              sbp;
  ButtoN                    b;
  DialoG                    sample_dlg = NULL;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  mp = (ApplyEditConvertRemovePtr) MemNew (sizeof (ApplyEditConvertRemoveData));
  if (mp == NULL) return;
  
  if (first_action_choice < AECR_APPLY || first_action_choice > AECR_REMOVE)
  {
    first_action_choice = AECR_APPLY;
  }
  mp->crippled = crippled;
  mp->crippled_action = first_action_choice;
  
  w = FixedWindow (-50, -33, -10, -10, title, StdCloseWindowProc); 
  SetObjectExtra (w, mp, StdCleanupExtraProc); 
  SetObjectExtra (w, mp, NULL);
  mp->form = (ForM) w;
  mp->input_entityID = bfp->input_entityID;

#ifndef WIN_MAC
  CreateStandardEditMenu (w);
#endif
  
  n = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (n, 10, 10);
  
  k = HiddenGroup (n, 2, 0, NULL);
  SetGroupSpacing (k, 10, 10);
  if (show_sample && make_sample_dlg != NULL && ! mp->crippled)
  {
    sample_dlg = make_sample_dlg (k, mp->input_entityID);
  }
  
  h = HiddenGroup (k, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);  
    
  if (!mp->crippled)
  {
    m = HiddenGroup (h, 2, 0, NULL);
    SetGroupSpacing (m, 10, 10);
  
    if (!show_sample && make_sample_dlg != NULL)
    {
      sbp = (SampleBtnPtr) MemNew (sizeof (SampleBtnData));
      sbp->entityID = mp->input_entityID;
      sbp->title = title;
      sbp->make_sample_dlg = make_sample_dlg;
      b = PushButton (m, "Show Qualifiers", AECRSampleButton);
      SetObjectExtra (b, sbp, StdCleanupExtraProc);
    }
 
    mp->action_popup = PopupList (m, TRUE, ChangeApplyEditConvertRemoveActionPopup);
    SetObjectExtra (mp->action_popup, mp, NULL);
    PopupItem (mp->action_popup, "Apply");
    PopupItem (mp->action_popup, "Edit");
    PopupItem (mp->action_popup, "Convert");
    PopupItem (mp->action_popup, "Swap");
    PopupItem (mp->action_popup, "Parse");
    PopupItem (mp->action_popup, "Remove");
    SetValue (mp->action_popup, first_action_choice);
  
    mp->change_action = change_action;
  
    mp->clear_constraints_on_action_change = CheckBox (h, "Clear when changing actions", NULL);
    SetStatus (mp->clear_constraints_on_action_change, TRUE);
  }
  
  g1 = HiddenGroup (h, 0, 0, NULL);
  for (page_num = 0; page_num < NUM_AECR; page_num++)
  {
    mp->aecr_pages [page_num] = make_dlg (g1, page_num + 1, 
                                          ApplyEditConvertRemoveChangeNotify,
                                          mp, mp->input_entityID);
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) mp->aecr_pages[0],
                              (HANDLE) mp->aecr_pages[1],
                              (HANDLE) mp->aecr_pages[2],
                              (HANDLE) mp->aecr_pages[3],
                              (HANDLE) mp->aecr_pages[4],
                              (HANDLE) mp->aecr_pages[5],
                               NULL);

  if (! mp->crippled)
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) m,
                                (HANDLE) mp->clear_constraints_on_action_change,
                                (HANDLE) g1,
                                (HANDLE) NULL);
    mp->constraints = FilterGroup (n, string_constraint, source_constraint, 
                                   location_constraint, cds_gene_prot_constraint,
                                   id_list_constraint,
                                   string_constraint_label);
  
  }
  
  mp->accept_cancel = AcceptCancelDialog (n, accept_actn, cancel_actn, 
                                          ApplyEditConvertRemoveClear, 
                                          ApplyEditConvertRemoveClearText, 
                                          (Pointer)mp, w);
                                          
  if (mp->crippled)
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) k, (HANDLE) mp->accept_cancel, NULL);
  }
  else
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) k,
                                (HANDLE) mp->constraints,
                                (HANDLE) mp->accept_cancel,
                                NULL);
  }
                              
  mp->prev_page = -1;                              
  ChangeApplyEditConvertRemoveAction (mp, first_action_choice);
  /* initialize sample values if conditions are met */
  if (sample_dlg != NULL && (presample_check_proc == NULL || presample_check_proc(mp->input_entityID))) {
    SendMessageToDialog (sample_dlg, VIB_MSG_INIT);
  }
  Show (w);   
}


static void CountFeaturesCallback(SeqFeatPtr sfp, Pointer userdata)
{
  Int4Ptr pTotal;
  
  if (userdata == NULL || sfp == NULL) {
      return;
  }
  
  pTotal = (Int4Ptr) userdata;
  
  (*pTotal) ++;
}


static Boolean CheckFeaturesForPresample (Uint2 entityID)
{
  SeqEntryPtr sep;
  Int4        total = 0;
  
  sep = GetTopSeqEntryForEntityID(entityID);
  if (sep == NULL) {
    return FALSE;
  }
  
  VisitFeaturesInSep(sep, &total, CountFeaturesCallback);
  
  if (total > 1500) {
    return FALSE;
  } else {
    return TRUE;
  }
}


typedef struct conversionconflict 
{
  ConvertFieldPtr cfp;
  GetSamplePtr    gsp;
} ConversionConflictData, PNTR ConversionConflictPtr;

static void GetConversionConflictsFeatureCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  ConversionConflictPtr ccp;
  CharPtr               src_str, dst_str;
  
  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  
  ccp = (ConversionConflictPtr) userdata;
  
  if (ccp->cfp == NULL
      || ccp->cfp->src_field_list == NULL 
      || ccp->cfp->dst_field_list == NULL
      || ccp->cfp->get_str_func == NULL
      || ccp->gsp == NULL)
  {
    return;
  }
  
  src_str = (ccp->cfp->get_str_func) (sfp, ccp->cfp->src_field_list, ccp->cfp->fsp);
  
  if (ccp->cfp->convert_type == CONVERT_TYPE_PARSE) 
  {
    dst_str = src_str;
    src_str = ReplaceStringForParse(src_str, ccp->cfp->text_portion);
    dst_str = MemFree (dst_str);
    if (src_str == NULL) 
    {
      return;
    }
  }
  
  if (!StringHasNoText (src_str))
  {
    dst_str = (ccp->cfp->get_str_func) (sfp, ccp->cfp->dst_field_list, ccp->cfp->fsp);
    if (!StringHasNoText (dst_str))
    {
      ccp->gsp->num_found ++;
      if (ccp->gsp->sample_text == NULL)
      {
        ccp->gsp->sample_text = dst_str;
      }
      else
      {
        if (StringCmp (dst_str, ccp->gsp->sample_text) != 0)
        {
          ccp->gsp->all_same = FALSE;
        }
        dst_str = MemFree (dst_str);
      }
    }
    src_str = MemFree (src_str);   
  }
}

static void GetConversionConflictsDescriptorCallback (SeqDescPtr sdp, Pointer userdata, FilterSetPtr fsp)
{
  ConversionConflictPtr ccp;
  CharPtr               src_str, dst_str;
  
  if (sdp == NULL || userdata == NULL)
  {
    return;
  }
  
  ccp = (ConversionConflictPtr) userdata;
  
  if (ccp->cfp == NULL
      || ccp->cfp->src_field_list == NULL 
      || ccp->cfp->dst_field_list == NULL
      || ccp->cfp->get_d_str_func == NULL
      || ccp->gsp == NULL)
  {
    return;
  }
  
  src_str = (ccp->cfp->get_d_str_func) (sdp, ccp->cfp->src_field_list, ccp->cfp->fsp);
  
  if (ccp->cfp->convert_type == CONVERT_TYPE_PARSE) 
  {
    dst_str = src_str;
    src_str = ReplaceStringForParse(src_str, ccp->cfp->text_portion);
    dst_str = MemFree (dst_str);
    if (src_str == NULL) 
    {
      return;
    }  
  }
  
  if (!StringHasNoText (src_str))
  {
    dst_str = (ccp->cfp->get_d_str_func) (sdp, ccp->cfp->dst_field_list, ccp->cfp->fsp);
    if (!StringHasNoText (dst_str))
    {
      ccp->gsp->num_found ++;
      if (ccp->gsp->sample_text == NULL)
      {
        ccp->gsp->sample_text = dst_str;
      }
      else
      {
        if (StringCmp (dst_str, ccp->gsp->sample_text) != 0)
        {
          ccp->gsp->all_same = FALSE;
        }
        dst_str = MemFree (dst_str);
      }
    }
    src_str = MemFree (src_str);   
  }
}

static GetSamplePtr CheckForConversionExistingTextInSeqEntry
(SeqEntryPtr              sep,
 ConvertFieldPtr          cfp,
 FilterSetPtr             fsp,
 Uint1                    seqfeat_choice,
 Uint1                    featdef_choice,
 Uint1                    descr_choice)
{
  ConversionConflictData ccd;
  
  ccd.cfp = cfp;
  ccd.gsp = (GetSamplePtr) MemNew (sizeof (GetSampleData));
  if (ccd.gsp != NULL)
  {
    ccd.gsp->sample_text = NULL;
    ccd.gsp->fieldstring_func = NULL;
    ccd.gsp->descrstring_func = NULL;
    ccd.gsp->free_vn_proc = NULL;
    ccd.gsp->copy_vn_proc = NULL;
    ccd.gsp->requested_field = NULL;
    ccd.gsp->num_found = 0;
    ccd.gsp->all_same = TRUE;
    OperateOnSeqEntryConstrainedObjects (sep, fsp, GetConversionConflictsFeatureCallback, 
                                         GetConversionConflictsDescriptorCallback,
                                         seqfeat_choice, featdef_choice, descr_choice, &ccd);
  }
  return ccd.gsp;  
}

NLM_EXTERN void 
AddEditApplyDataToApplyValue 
(Int4 action_choice,
 EditApplyPtr edit_apply, 
 ApplyValuePtr avp)
{
  if (avp == NULL)
  {
    return;
  }
  
  if (edit_apply == NULL)
  {
    avp->text_to_replace = NULL;
    avp->new_text = NULL;
    avp->where_to_replace = EditApplyFindLocation_anywhere;

    return;
  }

  if (action_choice == AECR_EDIT)
  {
    avp->text_to_replace = StringSave (edit_apply->find_txt);
    avp->new_text = StringSave (edit_apply->repl_txt);
    avp->where_to_replace = edit_apply->find_location;
  }
  else if (action_choice == AECR_APPLY)
  {
    avp->new_text = StringSave (edit_apply->apply_txt);
    avp->text_to_replace = NULL; 
    avp->where_to_replace = EditApplyFindLocation_anywhere;
  }
  else
  {
    avp->text_to_replace = NULL;
    avp->new_text = NULL;
    avp->where_to_replace = EditApplyFindLocation_anywhere;
  } 
}

#define AECR_DATA_BLOCK \
  ValNodePtr      field_list;           \
  ValNodePtr      field_list_to;        \
  ValNodePtr      subtype_list;         \
  Boolean         leave_on_original;    \
  Boolean         strip_name_from_text; \
  EditApplyPtr    edit_apply;           \
  TextPortionXPtr  text_portion;         \
  Boolean         remove_parsed;        \
  FreeValNodeProc free_field_vn_proc;   \
  FreeValNodeProc free_subtype_vn_proc;

typedef struct simpleaecr
{
  AECR_DATA_BLOCK
} SimpleAECRData, PNTR SimpleAECRPtr;

static void AECRDataBlockFreeContents (SimpleAECRPtr sp)
{
  if (sp != NULL)
  {
    sp->field_list = ValNodeFuncFree (sp->field_list, sp->free_field_vn_proc);
    sp->field_list_to = ValNodeFuncFree (sp->field_list_to, sp->free_field_vn_proc);
    sp->subtype_list = ValNodeFuncFree (sp->subtype_list, sp->free_subtype_vn_proc);
    sp->edit_apply = EditApplyFree (sp->edit_apply);
    sp->text_portion = TextPortionXFree(sp->text_portion);
  }
}

static SimpleAECRPtr SimpleAECRFree (SimpleAECRPtr sp)
{
  if (sp != NULL)
  {
    AECRDataBlockFreeContents (sp);
    sp = MemFree (sp);
  }
  return sp;
}

/* This function is used to autopopulate the new_text and find fields */
typedef CharPtr (*GetAutoPopulateTextFunc) PROTO ((ValNodePtr, ValNodePtr, Uint2));

#define AECR_BLOCK \
  DIALOG_MESSAGE_BLOCK                                \
  DialoG          field_list;                         \
  DialoG          field_list_to;                      \
  DialoG          subtype_list;                       \
  ButtoN          leave_on_original;                  \
  ButtoN          strip_name_from_text;               \
  DialoG          edit_apply;                         \
  DialoG          text_portion;                       \
  ButtoN          remove_parsed;                      \
  Int4            action_choice;                      \
  /* These are necessary to produce the SimpleAECR structure for output */ \
  FreeValNodeProc free_field_vn_proc;                 \
  FreeValNodeProc free_subtype_vn_proc;               \
  /* These are used for autopopulating text fields */ \
  GetAutoPopulateTextFunc  get_autopopulate_text;     \
  Nlm_ChangeNotifyProc     change_notify;             \
  Pointer                  change_userdata;           \
  Uint2                    entityID;

typedef struct simpleaecrdlg
{
  AECR_BLOCK
} SimpleAECRDlgData, PNTR SimpleAECRDlgPtr;

static void ClearTextSimpleAECRDlg (SimpleAECRDlgPtr dlg)
{
  if (dlg != NULL)
  {
    PointerToDialog (dlg->edit_apply, NULL);
    PointerToDialog (dlg->text_portion, NULL);
  }
}

static void ResetSimpleAECRDlg (SimpleAECRDlgPtr dlg)
{
  if (dlg != NULL)
  {
    PointerToDialog (dlg->field_list, NULL);
    PointerToDialog (dlg->field_list_to, NULL);
    PointerToDialog (dlg->subtype_list, NULL);
    PointerToDialog (dlg->text_portion, NULL);
    SafeSetStatus (dlg->leave_on_original, FALSE);
    SafeSetStatus (dlg->strip_name_from_text, FALSE);
    SafeSetStatus (dlg->remove_parsed, FALSE);
    ClearTextSimpleAECRDlg (dlg);
  }
}

static void SimpleAECRToDialog (DialoG d, Pointer userdata)
{
  SimpleAECRDlgPtr dlg;
  SimpleAECRPtr    data;
  
  dlg = (SimpleAECRDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  ResetSimpleAECRDlg (dlg);
  data = (SimpleAECRPtr) userdata;
  if (data != NULL)
  {
    PointerToDialog (dlg->field_list, data->field_list);
    PointerToDialog (dlg->field_list_to, data->field_list_to);
    SendMessageToDialog (dlg->field_list_to, NUM_VIB_MSG + 1);
    PointerToDialog (dlg->subtype_list, data->subtype_list);
    PointerToDialog (dlg->edit_apply, data->edit_apply);
    PointerToDialog (dlg->text_portion, data->text_portion);
    SafeSetValue (dlg->leave_on_original, data->leave_on_original);
    SafeSetValue (dlg->strip_name_from_text, data->strip_name_from_text);
    SafeSetStatus (dlg->remove_parsed, data->remove_parsed);
    dlg->free_field_vn_proc = data->free_field_vn_proc;
    dlg->free_subtype_vn_proc = data->free_subtype_vn_proc;
  }
  else
  {
    SendMessageToDialog (dlg->field_list, NUM_VIB_MSG + 1);
    SendMessageToDialog (dlg->field_list_to, NUM_VIB_MSG + 1);
  }
}

static Pointer DialogToSimpleAECR (DialoG d)
{
  SimpleAECRDlgPtr dlg;
  SimpleAECRPtr    data;
  
  dlg = (SimpleAECRDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  data = (SimpleAECRPtr) MemNew (sizeof (SimpleAECRData));
  if (data != NULL)
  {
    data->field_list = DialogToPointer (dlg->field_list);
    data->field_list_to = DialogToPointer (dlg->field_list_to);
    data->subtype_list = DialogToPointer (dlg->subtype_list);
    data->edit_apply = DialogToPointer (dlg->edit_apply);
    data->text_portion = DialogToPointer (dlg->text_portion);
    if (dlg->remove_parsed == NULL) {
      data->remove_parsed = FALSE;
    } 
    else
    {
      data->remove_parsed = GetStatus (dlg->remove_parsed);
    }
    if (dlg->leave_on_original != NULL)
    {
      data->leave_on_original = GetStatus (dlg->leave_on_original);
    }
    else
    {
      data->leave_on_original = FALSE;
    }
    if (dlg->strip_name_from_text != NULL)
    {
      data->strip_name_from_text = GetStatus (dlg->strip_name_from_text);
    }
    else
    {
      data->strip_name_from_text = FALSE;
    }
    data->free_field_vn_proc = dlg->free_field_vn_proc;
    data->free_subtype_vn_proc = dlg->free_subtype_vn_proc;
  }
  return data;
}

static void SimpleAECRDialogChangeNotify (Pointer userdata)
{
  SimpleAECRDlgPtr dlg;
  CharPtr          autopopulate_text;
  EditApplyPtr     eap;
  ValNodePtr       subtype_list, field_list;

  dlg = (SimpleAECRDlgPtr) userdata;
  if (dlg == NULL)
  { 
    return;
  }
  
  if (dlg->get_autopopulate_text != NULL
      && (dlg->action_choice == AECR_APPLY || dlg->action_choice == AECR_EDIT))
  {
    eap = DialogToPointer (dlg->edit_apply);
    subtype_list = DialogToPointer (dlg->subtype_list);
    field_list = DialogToPointer (dlg->field_list);
    autopopulate_text = (dlg->get_autopopulate_text) (subtype_list, field_list, dlg->entityID);
    subtype_list = ValNodeFuncFree (subtype_list, dlg->free_subtype_vn_proc);
    field_list = ValNodeFuncFree (field_list, dlg->free_field_vn_proc);
    if (dlg->action_choice == AECR_APPLY)
    {
      eap->apply_txt = MemFree (eap->find_txt);
      eap->apply_txt = autopopulate_text;
    }
    else
    {
      eap->find_txt = MemFree (eap->find_txt);
      eap->find_txt = autopopulate_text;
    }
    PointerToDialog (dlg->edit_apply, eap);
    EditApplyFree (eap);
  }
  
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void SimpleAECRMessage (DialoG d, Int2 mssg)

{
  SimpleAECRDlgPtr  dlg;

  dlg = (SimpleAECRDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) 
    {
      case VIB_MSG_INIT :
        /* reset list */
        ResetSimpleAECRDlg (dlg);
        break;
      case VIB_MSG_ENTER :
        if (dlg->subtype_list != NULL)
        {
          Select (dlg->subtype_list);
        }
        else
        {
          Select (dlg->field_list);
        }
        break;
      case AECR_VIB_MSG_SET_DEFAULT :
        SendMessageToDialog (dlg->field_list, AECR_VIB_MSG_SET_DEFAULT);
        SendMessageToDialog (dlg->field_list_to, AECR_VIB_MSG_SET_DEFAULT);
        SendMessageToDialog (dlg->subtype_list, AECR_VIB_MSG_SET_DEFAULT);
        ClearTextSimpleAECRDlg (dlg);
        break;
      case AECR_VIB_MSG_AUTOPOPULATE:
        SimpleAECRDialogChangeNotify (dlg);
        break;
      case AECR_VIB_MSG_CLEAR_TEXT :
        ClearTextSimpleAECRDlg (dlg);
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestSimpleAECR (DialoG d)
{
  SimpleAECRDlgPtr dlg;
  ValNodePtr       err_list = NULL;
  ValNodePtr       total_err_list = NULL;
  
  dlg = (SimpleAECRDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return FALSE;
  }

  total_err_list = TestDialog (dlg->field_list);
  
  if (dlg->subtype_list != NULL)
  {
    err_list = TestDialog (dlg->subtype_list);
    total_err_list = ValNodeAppend (total_err_list, err_list);
  }
    
  err_list = NULL;  
  switch (dlg->action_choice)
  {
    case AECR_CONVERT:
    case AECR_SWAP:
      err_list = TestDialog (dlg->field_list_to);
      break;
    case AECR_APPLY:
    case AECR_EDIT:
      err_list = TestDialog (dlg->edit_apply);
      break;
  }

  total_err_list = ValNodeAppend (total_err_list, err_list);
    
  return total_err_list;
}

static void SimpleAECRButtonChange (ButtoN b)
{
  SimpleAECRDlgPtr dlg;

  dlg = (SimpleAECRDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) return;

  if (dlg->change_notify != NULL) 
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}

static DialoG SimpleAECRDialogEx
(GrouP                    h,
 Int4                     action_choice, 
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 FreeValNodeProc          free_field_vn_proc,
 FreeValNodeProc          free_subtype_vn_proc,
 Nlm_FieldListDlgProc     fieldlist_dlg_proc,
 Nlm_SubtypeListDlgProc   subtypelist_dlg_proc,
 GetAutoPopulateTextFunc  get_autopopulate_text,
 CharPtr                  field_label,
 CharPtr                  subtype_label,
 Boolean                  is_text,
 Boolean                  strip_name,
 Uint2                    entityID)
{
  SimpleAECRDlgPtr dlg;
  GrouP            p, g1;
  SeqEntryPtr      sep;
  
  dlg = (SimpleAECRDlgPtr) MemNew (sizeof (SimpleAECRDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = SimpleAECRToDialog;
  dlg->fromdialog = DialogToSimpleAECR;
  dlg->dialogmessage = SimpleAECRMessage;
  dlg->testdialog = TestSimpleAECR;
  dlg->action_choice = action_choice;
  dlg->free_field_vn_proc = free_field_vn_proc;
  dlg->free_subtype_vn_proc = free_subtype_vn_proc;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->get_autopopulate_text = get_autopopulate_text;
  dlg->entityID = entityID;
  dlg->strip_name_from_text = NULL;
  
  dlg->subtype_list = NULL;
  
  sep = GetTopSeqEntryForEntityID(entityID);
  
  if (action_choice == AECR_CONVERT || action_choice == AECR_SWAP || action_choice == AECR_PARSE)
  {
    if (subtypelist_dlg_proc == NULL)
    {
      g1 = HiddenGroup (p, 2, 0, NULL);
    }
    else
    {
      g1 = HiddenGroup (p, 3, 0, NULL);
      StaticPrompt (g1, subtype_label, 0, dialogTextHeight, systemFont, 'l');
    }
    StaticPrompt (g1, "From", 0, dialogTextHeight, systemFont, 'l');
    StaticPrompt (g1, "To", 0, dialogTextHeight, systemFont, 'l');
    
    if (subtypelist_dlg_proc != NULL)
    {
      dlg->subtype_list = (subtypelist_dlg_proc) (g1, TRUE, sep, change_notify, change_userdata);
    }
    dlg->field_list = fieldlist_dlg_proc (g1, FALSE, 
                                          change_notify, 
                                          change_userdata);
    dlg->field_list_to = fieldlist_dlg_proc (g1, FALSE, 
                                             change_notify, 
                                             change_userdata);
    if (action_choice == AECR_CONVERT && is_text)
    {
      dlg->leave_on_original = CheckBox (p, "Leave on original", SimpleAECRButtonChange);
      SetObjectExtra (dlg->leave_on_original, dlg, NULL);
    }
    if (strip_name)
    {
      dlg->strip_name_from_text = CheckBox (p, "Strip name from text", SimpleAECRButtonChange);
      SetObjectExtra (dlg->strip_name_from_text, dlg, NULL);
    }
  }
  else
  {
    if (subtypelist_dlg_proc == NULL)
    {
      g1 = HiddenGroup (p, 1, 0, NULL);
      
    }
    else
    {
      g1 = HiddenGroup (p, 2, 0, NULL);
      StaticPrompt (g1, subtype_label, 0, dialogTextHeight, systemFont, 'l');
    }
    StaticPrompt (g1, field_label, 0, dialogTextHeight, systemFont, 'l');
    if (action_choice == AECR_REMOVE)
    {
      if (subtypelist_dlg_proc != NULL)
      {
        dlg->subtype_list = (subtypelist_dlg_proc) (g1, TRUE, sep, change_notify, change_userdata);
      }
      dlg->field_list = fieldlist_dlg_proc (g1, TRUE, 
                                            change_notify, 
                                            change_userdata);
    }
    else
    {
      if (subtypelist_dlg_proc != NULL)
      {
        dlg->subtype_list = (subtypelist_dlg_proc) (g1, TRUE, sep,
                                                    SimpleAECRDialogChangeNotify,
                                                    dlg);
      }
      dlg->field_list = fieldlist_dlg_proc (g1, FALSE, 
                                            SimpleAECRDialogChangeNotify, 
                                            dlg);
    }
  }
  
  if (is_text && (action_choice == AECR_APPLY || action_choice == AECR_EDIT))
  {
    dlg->edit_apply = EditApplyDialog (p, 
                                       action_choice == AECR_APPLY ? eEditApplyChoice_Apply : eEditApplyChoice_Edit,
                                       "New Text", NULL,
                                       change_notify, change_userdata);
    dlg->text_portion = NULL;
    dlg->remove_parsed = NULL;
  }
  else if (action_choice == AECR_PARSE)
  {
    dlg->edit_apply = NULL;
    if (is_text) {
        dlg->text_portion = TextPortionXDialog(p);
    }
    dlg->remove_parsed = CheckBox (p, "Remove from parsed field", SimpleAECRButtonChange);
    SetObjectExtra (dlg->remove_parsed, dlg, NULL);
  }
  else
  {
    dlg->edit_apply = NULL;
    dlg->text_portion = NULL;
    dlg->remove_parsed = NULL;
  }

  if (action_choice == AECR_CONVERT)
  {
    if (is_text)
    {
      AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) dlg->leave_on_original, (HANDLE) dlg->strip_name_from_text, NULL);
    }
    else
    {
      AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) dlg->strip_name_from_text, NULL);
    }
  }
  else if (action_choice == AECR_SWAP)
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) dlg->strip_name_from_text, NULL);
  }
  else if (action_choice == AECR_PARSE)
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) dlg->text_portion, (HANDLE) dlg->remove_parsed, NULL);
  }
  else
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) dlg->edit_apply, NULL);
  }
  return (DialoG) p;
}

static DialoG SimpleAECRDialog
(GrouP                    h,
 Int4                     action_choice, 
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 FreeValNodeProc          free_field_vn_proc,
 FreeValNodeProc          free_subtype_vn_proc,
 Nlm_FieldListDlgProc     fieldlist_dlg_proc,
 Nlm_SubtypeListDlgProc   subtypelist_dlg_proc,
 GetAutoPopulateTextFunc  get_autopopulate_text,
 CharPtr                  field_label,
 CharPtr                  subtype_label,
 Boolean                  is_text,
 Uint2                    entityID)
{
  return SimpleAECRDialogEx(h, action_choice, change_notify, change_userdata,
                            free_field_vn_proc, free_subtype_vn_proc, 
                            fieldlist_dlg_proc, subtypelist_dlg_proc,
                            get_autopopulate_text,
                            field_label, subtype_label,
                            is_text, FALSE, entityID);
}

static CharPtr 
CDSGeneProtAutopopulate 
(ValNodePtr subtype_list,
 ValNodePtr field_list, 
 Uint2 entityID)
{
  GetSamplePtr gsp;
  CharPtr      rval;
  
  gsp = CheckForExistingText (entityID, field_list, 
                              GetCDSGeneProtField,
                              NULL,
                              NULL,
                              IntValNodeCopy,
                              NULL,
                              0, 0, 0);
  rval = gsp->sample_text;
  gsp->sample_text = NULL;
  GetSampleFree (gsp);
  return rval;
}

typedef struct cdsgeneprotaecrdlg
{
  AECR_BLOCK
  ButtoN also_mrna;
} CDSGeneProtAECRDlgData, PNTR CDSGeneProtAECRDlgPtr;

typedef struct cdsgeneprotaecr 
{
  AECR_DATA_BLOCK
  Boolean also_mrna;
} CDSGeneProtAECRData, PNTR CDSGeneProtAECRPtr;


static CDSGeneProtAECRPtr CDSGeneProtAECRFree (CDSGeneProtAECRPtr sp)
{
  if (sp != NULL)
  {
    AECRDataBlockFreeContents ((SimpleAECRPtr) sp);
    sp = MemFree (sp);
  }
  return sp;
}


static void PointerToCDSGeneProtAECRDialog (DialoG d, Pointer userdata)
{
  CDSGeneProtAECRDlgPtr dlg;
  CDSGeneProtAECRPtr    data;

  dlg = (CDSGeneProtAECRDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  data = (CDSGeneProtAECRPtr) userdata;
  SimpleAECRToDialog (d, userdata);
  if (data != NULL)
  {
    SetStatus (dlg->also_mrna, data->also_mrna);
  }
  
}


static Boolean ShowingProteinName (CDSGeneProtAECRDlgPtr dlg)
{
  ValNodePtr vnp_from, vnp_to;
  Int4       protein_name_field_num = 2 + num_gene_fields + num_mrna_fields;
  Boolean    rval = FALSE;

  if (dlg == NULL) return FALSE;

  vnp_from = DialogToPointer (dlg->field_list);
  vnp_to = DialogToPointer (dlg->field_list_to);

  if ((vnp_from != NULL && vnp_from->data.intvalue == protein_name_field_num)
      || (vnp_to != NULL && vnp_to->data.intvalue == protein_name_field_num))
  {
    rval = TRUE;
  }
  vnp_from = ValNodeFree (vnp_from);
  vnp_to = ValNodeFree (vnp_to);

  return rval;
}


static Pointer CDSGeneProtAECRDialogToPointer (DialoG d)
{
  CDSGeneProtAECRDlgPtr dlg;
  SimpleAECRPtr         subset;
  CDSGeneProtAECRPtr    data = NULL;

  dlg = (CDSGeneProtAECRDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  subset = DialogToSimpleAECR (d);
  if (subset != NULL)
  {
    data = (CDSGeneProtAECRPtr) MemNew (sizeof (CDSGeneProtAECRData));
    MemCpy (data, subset, sizeof (SimpleAECRData));
    MemSet (subset, 0, sizeof (SimpleAECRData));
    subset = MemFree (subset);
    if (ShowingProteinName (dlg))
    {
      data->also_mrna = GetStatus (dlg->also_mrna);
    }
    else
    {
      data->also_mrna = FALSE;
    }
  }
  return data;
}


static void CDSGeneProtAECRDialogChange (Pointer data)
{
  CDSGeneProtAECRDlgPtr dlg;

  dlg = (CDSGeneProtAECRDlgPtr) data;
  if (dlg == NULL) return;
  
  if (ShowingProteinName (dlg))
  {
    Show (dlg->also_mrna);
  }
  else
  {
    Hide (dlg->also_mrna);
  }

  SimpleAECRDialogChangeNotify (dlg);
}


static DialoG CDSGeneProtAECRDialog
(GrouP                    h,
 Int4                     action_choice, 
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Uint2                    entityID)
{
  CDSGeneProtAECRDlgPtr dlg;
  GrouP                 p, g1;
  SeqEntryPtr           sep;
  
  dlg = (CDSGeneProtAECRDlgPtr) MemNew (sizeof (CDSGeneProtAECRDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = PointerToCDSGeneProtAECRDialog;
  dlg->fromdialog = CDSGeneProtAECRDialogToPointer;
  dlg->dialogmessage = SimpleAECRMessage;
  dlg->testdialog = TestSimpleAECR;
  dlg->action_choice = action_choice;
  dlg->free_field_vn_proc = NULL;
  dlg->free_subtype_vn_proc = NULL;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->get_autopopulate_text = CDSGeneProtAutopopulate;
  dlg->entityID = entityID;
  dlg->strip_name_from_text = NULL;
  
  dlg->subtype_list = NULL;
  
  sep = GetTopSeqEntryForEntityID(entityID);
  
  if (action_choice == AECR_CONVERT || action_choice == AECR_SWAP || action_choice == AECR_PARSE)
  {
    g1 = HiddenGroup (p, 2, 0, NULL);
    StaticPrompt (g1, "From", 0, dialogTextHeight, systemFont, 'l');
    StaticPrompt (g1, "To", 0, dialogTextHeight, systemFont, 'l');
    
    dlg->field_list = CDSGeneProtFieldSelectionDialog (g1, FALSE, 
                                          CDSGeneProtAECRDialogChange, 
                                          dlg);
    dlg->field_list_to = CDSGeneProtFieldSelectionDialog (g1, FALSE, 
                                             CDSGeneProtAECRDialogChange, 
                                             dlg);
    if (action_choice == AECR_CONVERT)
    {
      dlg->leave_on_original = CheckBox (p, "Leave on original", NULL);
    }
  }
  else
  {
    g1 = HiddenGroup (p, 1, 0, NULL);
    StaticPrompt (g1, "Field", 0, dialogTextHeight, systemFont, 'l');
    if (action_choice == AECR_REMOVE)
    {
      dlg->field_list = CDSGeneProtFieldSelectionDialog (g1, TRUE, 
                                            CDSGeneProtAECRDialogChange, 
                                            dlg);
    }
    else
    {
      dlg->field_list = CDSGeneProtFieldSelectionDialog (g1, FALSE, 
                                            CDSGeneProtAECRDialogChange, 
                                            dlg);
    }
  }
  
  if (action_choice == AECR_APPLY || action_choice == AECR_EDIT)
  {
    dlg->edit_apply = EditApplyDialog (p, 
                                       action_choice == AECR_APPLY ? eEditApplyChoice_Apply : eEditApplyChoice_Edit,
                                       "New Text", NULL,
                                       change_notify, change_userdata);
    dlg->text_portion = NULL;
    dlg->remove_parsed = NULL;
  }
  else if (action_choice == AECR_PARSE)
  {
    dlg->edit_apply = NULL;
    dlg->text_portion = TextPortionXDialog(p);
    dlg->remove_parsed = CheckBox (p, "Remove from parsed field", NULL);
  }
  else
  {
    dlg->edit_apply = NULL;
    dlg->text_portion = NULL;
    dlg->remove_parsed = NULL;
  }

  dlg->also_mrna = CheckBox (p, "Make mRNA product match CDS protein name", NULL);

  if (action_choice == AECR_CONVERT)
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) dlg->leave_on_original, (HANDLE) dlg->also_mrna, NULL);
  }
  else if (action_choice == AECR_SWAP)
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) dlg->strip_name_from_text, (HANDLE) dlg->also_mrna, NULL);
  }
  else if (action_choice == AECR_PARSE)
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) dlg->text_portion, (HANDLE) dlg->remove_parsed, (HANDLE) dlg->also_mrna, NULL);
  }
  else
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) dlg->edit_apply, (HANDLE) dlg->also_mrna, NULL);
  }
  return (DialoG) p;
}


static void AddStringToSample (CharPtr str, GetSamplePtr gsp)
{
  if (gsp == NULL) return;

  if (str != NULL)
  {
    gsp->num_found ++;
  }
  if (gsp->sample_text == NULL)
  {
    gsp->sample_text = StringSave (str);
  }
  else
  {
    if (StringCmp (str, gsp->sample_text) != 0)
    {
      gsp->all_same = FALSE;
    }
  }
}

static GetSamplePtr 
GetCDSGeneProtAECRSample 
(ValNodePtr cdset_list,
 ValNodePtr field_list, 
 ValNodePtr field_list_to,
 FilterSetPtr fsp)
{
  GetSamplePtr gsp;
  ValNodePtr   cdset_vnp, sfp_vnp;
  CharPtr      str;
  CDSetPtr       cdsp;
 
  gsp = (GetSamplePtr) MemNew (sizeof (GetSampleData));
  if (gsp == NULL)
  {
    return NULL;
  }
  gsp->free_vn_proc = NULL;
  gsp->copy_vn_proc = IntValNodeCopy;
  for (cdset_vnp = cdset_list; cdset_vnp != NULL; cdset_vnp = cdset_vnp->next)
  {
    cdsp = (CDSetPtr) cdset_vnp->data.ptrvalue;
    if (IsCDSetMatPeptideQualChoice(field_list) || IsCDSetMatPeptideQualChoice(field_list_to))
    {
      /* need to look at all mat_peptide features */
      for (sfp_vnp = cdsp->prot_list; sfp_vnp != NULL; sfp_vnp = sfp_vnp->next)
      {
        if (fsp != NULL && fsp->cgp != NULL && DoesConstraintDisqualifyFeature (sfp_vnp->data.ptrvalue, fsp->cgp))
        {
          /* this feature will be skipped in the action, so it should not be included in the sample */
        }
        else
        {
          str = GetCDSGeneProtField (sfp_vnp->data.ptrvalue, field_list, fsp);
          if (str != NULL && field_list_to != NULL)
          {
            str = GetCDSGeneProtField (sfp_vnp->data.ptrvalue, field_list_to, fsp);
          }
          AddStringToSample (str, gsp);     
          str = MemFree (str);
        }
      }
    } else {
      str = GetCDSetField (cdsp, field_list);
      if (str != NULL && field_list_to != NULL)
      {
        str = MemFree (str);
        str = GetCDSetField (cdset_vnp->data.ptrvalue, field_list_to);
      }
      AddStringToSample (str, gsp);     
      str = MemFree (str);
    }
  }
  return gsp;
}

static Boolean FeatureTypeMissingFromCDSet (CDSetPtr cdsp, Uint2 choice)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  Boolean    is_missing = TRUE;

  if (cdsp == NULL) return TRUE;
  switch (choice)
  {
    case FEATDEF_CDS:
      if (cdsp->cds_list != NULL)
      {
        is_missing = FALSE;
      }
      break;
    case FEATDEF_GENE:
      if (cdsp->gene_list != NULL)
      {
        is_missing = FALSE;
      }
      break;
    case FEATDEF_mRNA:
      if (cdsp->mrna_list != NULL)
      {
        is_missing = FALSE;
      }
      break;
    case FEATDEF_PROT:
    case FEATDEF_mat_peptide_aa:
      for (vnp = cdsp->prot_list; vnp != NULL; vnp = vnp->next)
      {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->idx.subtype == choice)
        {
          is_missing = FALSE;
        }
      }
      break;
  }
  return is_missing;
}

static Int4 
CountMissingFeatures 
(ValNodePtr cdset_list,
 ValNodePtr field_list, 
 ValNodePtr field_list_to,
 FilterSetPtr fsp)
{
  ValNodePtr   cdset_vnp, sfp_vnp;
  CharPtr      str;
  CDSetPtr     cdsp;
  Uint2        choice_from, choice_to;
  Int4         num_missing = 0;

  choice_from = FeatDefTypeFromFieldList (field_list);
  choice_to = FeatDefTypeFromFieldList (field_list_to);

  if (choice_from == 0 || choice_to == 0 || choice_from == choice_to)
  {
    return 0;
  }

  for (cdset_vnp = cdset_list; cdset_vnp != NULL; cdset_vnp = cdset_vnp->next)
  {
    cdsp = (CDSetPtr) cdset_vnp->data.ptrvalue;

    if (choice_from == FEATDEF_mat_peptide_aa)
    {
      /* need to look at all mat_peptide features */
      for (sfp_vnp = cdsp->prot_list; sfp_vnp != NULL; sfp_vnp = sfp_vnp->next)
      {
        if (fsp != NULL && fsp->cgp != NULL && DoesConstraintDisqualifyFeature (sfp_vnp->data.ptrvalue, fsp->cgp))
        {
          /* this feature will be skipped in the action, so it doesn't count towards missing sets */
        }
        else
        {
          str = GetCDSGeneProtField (sfp_vnp->data.ptrvalue, field_list, fsp);
          if (str != NULL)
          {
            if (FeatureTypeMissingFromCDSet (cdsp, choice_to)) 
            {
              num_missing++;
            }
          }
          str = MemFree (str);
        }
      }
    } else {
      str = GetCDSetField (cdsp, field_list);
      if (str != NULL)
      {
        if (FeatureTypeMissingFromCDSet (cdsp, choice_to)) 
        {
          num_missing++;
        }
      }
      str = MemFree (str);
    }
  }
  return num_missing;
}


/* things get tricksy when we're dealing with mat_peptides, because there
 * can be more than one of them in the CDS-mRNA-Gene-Prot group.
 */

static Boolean CDSGeneProtApplyToAllThatMatchFilter 
(CharPtr         str_src,
 CDSetPtr        cdsp,
 ValNodePtr      field_list_to,
 FilterSetPtr    fsp,
 ExistingTextPtr etp)
{
  ApplyValueData avd;

  avd.where_to_replace = EditApplyFindLocation_anywhere;

  /* build gene if we're applying to gene but gene doesn't already exist */
  BuildNewCDSGeneProtFeature (cdsp, field_list_to, NULL);
  avd.field_list = NULL;
  avd.new_text = str_src;
  avd.text_to_replace = NULL;
  avd.etp = etp;
          
  return SetCDSetField (cdsp, field_list_to, &avd, fsp);
}

static Boolean CDSGeneProtApplyToOneFeature
(CharPtr         str_src,
 SeqFeatPtr      sfp,
 ValNodePtr      field_list_to,
 FilterSetPtr    fsp,
 ExistingTextPtr etp)
{
  ApplyValueData avd;

  avd.where_to_replace = EditApplyFindLocation_anywhere;

  avd.field_list = NULL;
  avd.new_text = str_src;
  avd.text_to_replace = NULL;
  avd.etp = etp;
          
  return SetCDSGeneProtField (sfp, field_list_to, &avd, fsp);
}


static void CopyProteinNameTomRNAName (CDSetPtr cdsp)
{
  Int4    protein_name_field_num = 2 + num_gene_fields + num_mrna_fields;
  Int4    mrna_product_field_num = 2 + num_gene_fields;
  ValNode vn_src, vn_dst;
  ApplyValueData avd;

  vn_src.data.intvalue = protein_name_field_num;
  vn_src.next = NULL;
  vn_src.choice = 0;

  vn_dst.data.intvalue = mrna_product_field_num;
  vn_dst.next = NULL;
  vn_dst.choice = 0;

  avd.new_text = GetCDSetField (cdsp, &vn_src);
  avd.etp = NULL;
  avd.field_list = NULL;
  avd.text_to_replace = NULL;
  avd.where_to_replace = EditApplyFindLocation_anywhere;
  SetCDSetField (cdsp, &vn_dst, &avd, NULL);
}

static void CDSGeneProtConvert 
(CDSetPtr        cdsp,
 CDSGeneProtAECRPtr   sp,
 FilterSetPtr    fsp,
 ExistingTextPtr etp,
 Boolean         build_missing)
{
  CharPtr        str_src = NULL;
  CharPtr        str_dst = NULL;
  ValNodePtr     sfp_vnp;
  SeqFeatPtr     sfp;
  Boolean        able_to_apply = FALSE;
  Boolean        change_made = FALSE;

  if (IsCDSetMatPeptideQualChoice(sp->field_list)) 
  {
    for (sfp_vnp = cdsp->prot_list; sfp_vnp != NULL; sfp_vnp = sfp_vnp->next)
    {
      sfp = sfp_vnp->data.ptrvalue;
      able_to_apply = FALSE;
      if (fsp == NULL || !DoesConstraintDisqualifyFeature(sfp, fsp->cgp)) {
        str_src = GetCDSGeneProtField (sfp_vnp->data.ptrvalue, sp->field_list, NULL);
        if (str_src != NULL)
        {
          if (IsCDSetMatPeptideQualChoice (sp->field_list_to)) 
          {
            CDSGeneProtApplyToOneFeature (str_src, sfp, sp->field_list_to, fsp, etp);
            change_made = TRUE;
          } 
          else 
          {
            if (build_missing)
            {
              BuildNewCDSGeneProtFeature (cdsp, sp->field_list_to, sp->field_list);
            }
            able_to_apply = CDSGeneProtApplyToAllThatMatchFilter (str_src, cdsp, sp->field_list_to, fsp, etp);
            if (able_to_apply)
            {
              change_made = TRUE;
            }
          }
          if (!sp->leave_on_original && able_to_apply)
          {
            sp->field_list->choice = 1;
            RemoveCDSGeneProtField (sfp, sp->field_list, fsp);
            change_made = TRUE;
          }
        }
      }
    }    
  }
  else 
  {
    str_src = GetCDSetField (cdsp, sp->field_list);
    if (str_src != NULL)
    {
      if (build_missing)
      {
        BuildNewCDSGeneProtFeature (cdsp, sp->field_list_to, sp->field_list);
      }
      able_to_apply = CDSGeneProtApplyToAllThatMatchFilter (str_src, cdsp, sp->field_list_to, fsp, etp);
      if (!sp->leave_on_original && able_to_apply)
      {
        sp->field_list->choice = 1;
        RemoveCDSetField (cdsp, sp->field_list, fsp);
      }
      if (able_to_apply)
      {
        change_made = TRUE;
      }
    }   
  }
  if (sp->also_mrna && change_made)
  {
    CopyProteinNameTomRNAName (cdsp);
  }

}

static void CDSGeneProtSwap
(CDSetPtr        cdsp,
 CDSGeneProtAECRPtr   sp,
 FilterSetPtr    fsp,
 ExistingTextPtr etp,
 Boolean         build_missing)
{
  CharPtr        str_src = NULL;
  CharPtr        str_dst = NULL;
  ValNodePtr     sfp_vnp_src, sfp_vnp_dst, vnp_swap;
  SeqFeatPtr     sfp_src, sfp_dst;
  ApplyValueData avd;
  Uint2          choice_to, choice_from;
  Boolean        change_made = FALSE;

  avd.where_to_replace = EditApplyFindLocation_anywhere;

  choice_from = FeatDefTypeFromFieldList (sp->field_list);
  choice_to = FeatDefTypeFromFieldList (sp->field_list_to);

  if (choice_from == FEATDEF_mat_peptide_aa || choice_to == FEATDEF_mat_peptide_aa)
  {
    /* make sure mat_peptide qual if in field_list */
    if (!IsCDSetMatPeptideQualChoice(sp->field_list)) 
    {
      vnp_swap = sp->field_list;
      sp->field_list = sp->field_list_to;
      sp->field_list_to = vnp_swap;
    }
    for (sfp_vnp_src = cdsp->prot_list; sfp_vnp_src != NULL; sfp_vnp_src = sfp_vnp_src->next)
    {
      sfp_src = sfp_vnp_src->data.ptrvalue;
      if (fsp == NULL || !DoesConstraintDisqualifyFeature(sfp_src, fsp->cgp)) 
      {
        str_src = GetCDSGeneProtField (sfp_src, sp->field_list, NULL);
        if (IsCDSetMatPeptideQualChoice(sp->field_list_to)) 
        {
          for (sfp_vnp_dst = cdsp->prot_list; sfp_vnp_dst != NULL; sfp_vnp_dst = sfp_vnp_dst->next)
          {
            sfp_dst = sfp_vnp_dst->data.ptrvalue;
            if (fsp == NULL || !DoesConstraintDisqualifyFeature(sfp_dst, fsp->cgp))
            {
              str_dst = GetCDSGeneProtField (sfp_dst, sp->field_list_to, NULL);
              if (!StringHasNoText (str_src) || !StringHasNoText (str_dst))
              {
                if (build_missing)
                {
                  BuildNewCDSGeneProtFeature (cdsp, sp->field_list, sp->field_list_to);
                  BuildNewCDSGeneProtFeature (cdsp, sp->field_list_to, sp->field_list);
                }
                if (!FeatureTypeMissingFromCDSet(cdsp, choice_from) 
                    && !FeatureTypeMissingFromCDSet(cdsp, choice_to))
                {
                  avd.etp = NULL;
                  avd.text_to_replace = NULL;
                  avd.new_text = StringSave (str_src);
                  avd.field_list = NULL;
                  SetCDSGeneProtField (sfp_dst, sp->field_list_to, &avd, fsp);
                  avd.new_text = StringSave (str_dst);
                  SetCDSGeneProtField (sfp_src, sp->field_list, &avd, fsp);
                  change_made = TRUE;
                }
              }
              str_dst = MemFree (str_dst);
            }
          }
        } 
        else
        {
          str_dst = GetCDSetField (cdsp, sp->field_list_to);
          if (!StringHasNoText (str_src) || !StringHasNoText (str_dst))
          {
            if (build_missing) 
            {
              BuildNewCDSGeneProtFeature (cdsp, sp->field_list, sp->field_list_to);
              BuildNewCDSGeneProtFeature (cdsp, sp->field_list_to, sp->field_list);
            }
            if (!FeatureTypeMissingFromCDSet(cdsp, choice_from)
                && !FeatureTypeMissingFromCDSet(cdsp, choice_to))
            {
              avd.etp = NULL;
              avd.text_to_replace = NULL;
              avd.new_text = StringSave (str_src);
              avd.field_list = NULL;
              SetCDSetField (cdsp, sp->field_list_to, &avd, fsp);
              avd.new_text = StringSave (str_dst);
              SetCDSGeneProtField (sfp_src, sp->field_list, &avd, fsp);
              change_made = TRUE;
            }
          }
          str_dst = MemFree (str_dst);
        }
        str_src = MemFree (str_src);
      }
    } 
  }
  else
  {
    str_src = GetCDSetField (cdsp, sp->field_list);
    str_dst = GetCDSetField (cdsp, sp->field_list_to);
    if (!StringHasNoText (str_src) || !StringHasNoText (str_dst))
    {
      if (build_missing)
      {
        BuildNewCDSGeneProtFeature (cdsp, sp->field_list, sp->field_list_to);
        BuildNewCDSGeneProtFeature (cdsp, sp->field_list_to, sp->field_list);
      }
      if (!FeatureTypeMissingFromCDSet(cdsp, choice_from)
          && !FeatureTypeMissingFromCDSet(cdsp, choice_to))
      {
        avd.etp = NULL;
        avd.text_to_replace = NULL;
        avd.new_text = str_src;
        avd.field_list = NULL;
        SetCDSetField (cdsp, sp->field_list_to, &avd, fsp);
        avd.new_text = str_dst;
        SetCDSetField (cdsp, sp->field_list, &avd, fsp);
        change_made = TRUE;
      }
    }
  }
  if (change_made && sp->also_mrna)
  {
    CopyProteinNameTomRNAName (cdsp);
  }
}


static void CDSGeneProtParse
(CDSetPtr        cdsp,
 CDSGeneProtAECRPtr   sp,
 FilterSetPtr    fsp,
 ExistingTextPtr etp,
 Boolean         build_missing)
{
  CharPtr        str_src = NULL;
  CharPtr        str_dst = NULL;
  ValNodePtr     sfp_vnp;
  SeqFeatPtr     sfp;
  ApplyValueData avd;
  Boolean        able_to_apply;
  Boolean        change_made = FALSE;

  avd.where_to_replace = EditApplyFindLocation_anywhere;

  if (IsCDSetMatPeptideQualChoice(sp->field_list)) 
  {
    for (sfp_vnp = cdsp->prot_list; sfp_vnp != NULL; sfp_vnp = sfp_vnp->next)
    {
      sfp = sfp_vnp->data.ptrvalue;
      if (fsp == NULL || !DoesConstraintDisqualifyFeature(sfp, fsp->cgp)) 
      {
        str_src = GetCDSGeneProtField (sfp, sp->field_list, NULL);
        str_dst = str_src;
        str_src = ReplaceStringForParse(str_src, sp->text_portion);
        if (!sp->remove_parsed) {
          str_dst = MemFree (str_dst);
        }

        if (str_src != NULL) {
          if (build_missing) 
          {
            BuildNewCDSGeneProtFeature (cdsp, sp->field_list_to, sp->field_list);
          }
          able_to_apply = FALSE;
          if (IsCDSetMatPeptideQualChoice (sp->field_list_to)) 
          {
            able_to_apply = CDSGeneProtApplyToOneFeature (str_src, sfp, sp->field_list_to, fsp, etp);
          } 
          else 
          {
            able_to_apply = CDSGeneProtApplyToAllThatMatchFilter (str_src, cdsp, sp->field_list_to, fsp, etp);
          }
          if (sp->remove_parsed && able_to_apply) {
            CDSGeneProtApplyToOneFeature (str_dst, sfp, sp->field_list, fsp, NULL);
          }
          if (able_to_apply)
          {
            change_made = TRUE;
          }
          str_src = MemFree (str_src);
        }
        str_dst = MemFree (str_dst);
      }
    }
  }
  else
  {
    str_src = GetCDSetField (cdsp, sp->field_list);
    str_dst = str_src;
    str_src = ReplaceStringForParse(str_src, sp->text_portion);
    if (!sp->remove_parsed) {
      str_dst = MemFree (str_dst);
    }
    if (str_src != NULL) {
      if (build_missing) 
      {
        BuildNewCDSGeneProtFeature (cdsp, sp->field_list_to, sp->field_list);
      }
      able_to_apply = CDSGeneProtApplyToAllThatMatchFilter (str_src, cdsp, sp->field_list_to, fsp, etp);
      if (sp->remove_parsed && able_to_apply) {
        CDSGeneProtApplyToAllThatMatchFilter (str_dst, cdsp, sp->field_list, fsp, NULL);
      }
      if (able_to_apply)
      {
        change_made = TRUE;
      }
      str_src = MemFree (str_src);
    }   
    str_dst = MemFree (str_dst);
  }

  if (sp->also_mrna && change_made)
  {
    CopyProteinNameTomRNAName (cdsp);
  }
}     

static void CDSGeneProtEdit
(CDSetPtr        cdsp,
 CDSGeneProtAECRPtr   sp,
 FilterSetPtr    fsp,
 ExistingTextPtr etp)
{
  ApplyValueData avd;

  avd.where_to_replace = EditApplyFindLocation_anywhere;

  avd.field_list = sp->field_list;
  avd.etp = NULL;
  AddEditApplyDataToApplyValue (AECR_EDIT, sp->edit_apply, &avd);
  SetCDSetField (cdsp, sp->field_list, &avd, fsp);
  avd.text_to_replace = MemFree (avd.text_to_replace);
  avd.new_text = MemFree (avd.new_text);
  if (sp->also_mrna)
  {
    CopyProteinNameTomRNAName (cdsp);
  }
}

static ValNodePtr CollectCDSetList (Uint2 entityID, FilterSetPtr fsp)
{
  ValNodePtr cdset_list;

  /* need to create a list of sets */
  cdset_list = BuildCDSet2List (entityID, fsp == NULL || fsp->cgp == NULL ? NULL : fsp->cgp);

  return cdset_list;
}


static void MarkAlreadyEditedGenes (CDSetPtr cdsp, ValNodePtr cdset_list)
{
  ValNodePtr gene_vnp, list_vnp, genecheck_vnp;
  CDSetPtr   check_cdsp;

  if (cdsp == NULL || cdset_list == NULL || cdsp->gene_list == NULL) return;

  for (gene_vnp = cdsp->gene_list; gene_vnp != NULL; gene_vnp = gene_vnp->next)
  {
    if (gene_vnp->choice > 0) continue;
    for (list_vnp = cdset_list; list_vnp != NULL; list_vnp = list_vnp->next)
    {
      check_cdsp = (CDSetPtr) list_vnp->data.ptrvalue;
      if (check_cdsp != NULL) 
      {
        for (genecheck_vnp = check_cdsp->gene_list; genecheck_vnp != NULL; genecheck_vnp = genecheck_vnp->next)
        {
          if (genecheck_vnp->data.ptrvalue == gene_vnp->data.ptrvalue)
          {
            genecheck_vnp->choice = 1;
          }
        }
      }
    }
  }  
}


static Boolean CDSGeneProtAECRAction (Pointer userdata)
{
  ApplyEditConvertRemovePtr ap;
  CDSGeneProtAECRPtr        sp;
  Int4                      action_choice;
  SeqEntryPtr               sep;
  FilterSetPtr              fsp;
  GetSamplePtr              gsp = NULL;
  Boolean                   rval = TRUE;
  ValNodePtr                cdset_list, cdset_vnp;
  ExistingTextPtr           etp = NULL;
  Int4                      num_missing;
  Boolean                   build_missing = FALSE;
  MsgAnswer                 ans;

  ap = (ApplyEditConvertRemovePtr) userdata;
  if (ap == NULL)
  {
    return FALSE;
  }
  
  WatchCursor ();
  Update ();
  
  if (ap->crippled)
  {
    action_choice = ap->crippled_action;
  }
  else
  {
    action_choice = GetValue (ap->action_popup);
  }
  if (action_choice < 1 || action_choice > NUM_AECR)
  {
    return FALSE;
  }
  sp = (CDSGeneProtAECRPtr) DialogToPointer (ap->aecr_pages[action_choice - 1]);

  sep = GetTopSeqEntryForEntityID (ap->input_entityID);
  if (sp->field_list == NULL
      || ((action_choice == AECR_CONVERT || action_choice == AECR_SWAP) 
           && sp->field_list == NULL))
  {
    sp = CDSGeneProtAECRFree (sp);
    ArrowCursor ();
    Update ();
    return FALSE;
  }

  fsp = (FilterSetPtr) DialogToPointer (ap->constraints);
  
  /* need to create a list of sets */
  cdset_list = CollectCDSetList (ap->input_entityID, fsp);
  
  /* find existing text and how to handle it */
  if (action_choice == AECR_CONVERT || action_choice == AECR_APPLY || action_choice == AECR_PARSE)
  {
    gsp = GetCDSGeneProtAECRSample (cdset_list, sp->field_list, sp->field_list_to, fsp);
    etp = GetExistingTextHandlerInfo (gsp == NULL ? 0 : gsp->num_found, FALSE);
    gsp = GetSampleFree (gsp);
    if (etp != NULL && etp->existing_text_choice == eExistingTextChoiceCancel)
    {
      rval = FALSE;
    }
  }
  
  num_missing = CountMissingFeatures (cdset_list, sp->field_list, sp->field_list_to, fsp);
  if (action_choice == AECR_SWAP) 
  {
    num_missing += CountMissingFeatures (cdset_list, sp->field_list_to, sp->field_list, fsp);
  }
  if (num_missing > 0) 
  {
    ans = Message (MSG_YNC, "There are %d fields where the feature to which the field should be applied does not exist.  Do you want to create these features?", num_missing);
    if (ans == ANS_CANCEL)
    {
      rval = FALSE;
    } else if (ans == ANS_YES) {
      build_missing = TRUE;
    }
  }
  
  if (rval)
  {
    /* if we are applying gene locus, we need to make sure every sequence has a gene
     * and we need to recollect the cdset_list
     */
    if ((action_choice == AECR_APPLY || action_choice == AECR_PARSE)
        && sp->field_list != NULL
        && sp->field_list->data.intvalue == 1 + GENEFIELD_LOCUS)
    {
      BuildNewGenes (sep, &cdset_list, LoneGeneOkWithFilter(fsp));
    }
    
    for (cdset_vnp = cdset_list; cdset_vnp != NULL; cdset_vnp = cdset_vnp->next)
    {
      if (action_choice == AECR_CONVERT)
      {
        CDSGeneProtConvert (cdset_vnp->data.ptrvalue, sp, fsp, etp, build_missing);
      }
      else if (action_choice == AECR_SWAP)
      {
        CDSGeneProtSwap (cdset_vnp->data.ptrvalue, sp, fsp, etp, build_missing);
      }
      else if (action_choice == AECR_PARSE)
      {
        CDSGeneProtParse (cdset_vnp->data.ptrvalue, sp, fsp, etp, build_missing);
      }      
      else if (action_choice == AECR_APPLY)
      {
        CDSGeneProtApplyToAllThatMatchFilter (sp->edit_apply->apply_txt, cdset_vnp->data.ptrvalue, sp->field_list, fsp, etp);
      }
      else if (action_choice == AECR_REMOVE)
      {
        sp->field_list->choice = 1;
        RemoveCDSetField (cdset_vnp->data.ptrvalue, sp->field_list, fsp);
        if (sp->also_mrna)
        {
          CopyProteinNameTomRNAName (cdset_vnp->data.ptrvalue);
        }
      }
      else if (action_choice == AECR_EDIT)
      {
        CDSGeneProtEdit (cdset_vnp->data.ptrvalue, sp, fsp, etp);
      }
      if (((action_choice == AECR_APPLY || action_choice == AECR_EDIT) && IsCDSetGeneQualChoice (sp->field_list))
          || ((action_choice == AECR_PARSE || action_choice == AECR_CONVERT) && IsCDSetGeneQualChoice (sp->field_list_to))
          || (action_choice == AECR_SWAP && (IsCDSetGeneQualChoice (sp->field_list) || IsCDSetGeneQualChoice (sp->field_list_to))))
      {
        MarkAlreadyEditedGenes (cdset_vnp->data.ptrvalue, cdset_vnp->next);
      }
    }
  }
  FilterSetFree (fsp);
  sp = CDSGeneProtAECRFree (sp);
  FreeCDSetList (cdset_list);
  
  if (rval)
  {
    ObjMgrSetDirtyFlag (ap->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, ap->input_entityID, 0, 0);
  }
  ArrowCursor ();
  Update (); 
  
  return rval;

}

static void 
CDSGeneProtApplyEditConvertRemove 
(IteM    i, 
 Int4    first_action_choice, 
 Boolean crippled)
{
  ApplyEditConvertRemoveCombo (i, first_action_choice, crippled,
                               "CDS-Gene-Prot-mRNA",
                               CDSGeneProtAECRDialog, 
                               TRUE, GetCDSGeneProtSample,
                               FALSE, FALSE, FALSE, TRUE, FALSE, NULL,
                               CDSGeneProtAECRAction, NULL, NULL,
                               CheckFeaturesForPresample);
}

extern void ApplyCDSGeneProt (IteM i)
{
  CDSGeneProtApplyEditConvertRemove (i, AECR_APPLY, FALSE);
}

extern void PublicApplyCDSGeneProt (IteM i)
{
  CDSGeneProtApplyEditConvertRemove (i, AECR_APPLY, TRUE);
}

extern void EditCDSGeneProt (IteM i)
{
  CDSGeneProtApplyEditConvertRemove (i, AECR_EDIT, FALSE);
}

extern void PublicEditCDSGeneProt (IteM i)
{
  CDSGeneProtApplyEditConvertRemove (i, AECR_EDIT, TRUE);
}

extern void ConvertCDSGeneProt (IteM i)
{
  CDSGeneProtApplyEditConvertRemove (i, AECR_CONVERT, FALSE);
}

extern void SwapCDSGeneProt (IteM i)
{
  CDSGeneProtApplyEditConvertRemove (i, AECR_SWAP, FALSE);
}

extern void RemoveCDSGeneProt (IteM i)
{
  CDSGeneProtApplyEditConvertRemove (i, AECR_REMOVE, FALSE);
}

#define QUALKIND_SOURCE_QUAL 1
#define QUALKIND_STRINGS     2
#define QUALKIND_LOCATION    3
#define QUALKIND_ORIGIN      4

typedef struct sourcequalaecrdata
{
  Int4          qual_kind;
  SimpleAECRPtr source_qual;
  SimpleAECRPtr location;
  SimpleAECRPtr origin;
  SimpleAECRPtr strings;
  
  TaxnameOptionsPtr taxname_options;
} SourceQualAECRData, PNTR SourceQualAECRPtr;

static SourceQualAECRPtr SourceQualAECRFree (SourceQualAECRPtr sp)
{
  if (sp != NULL)
  {
    sp->source_qual = SimpleAECRFree (sp->source_qual);
    sp->taxname_options = MemFree (sp->taxname_options);
    sp = MemFree (sp);
  }
  return sp;
}

typedef struct sourcequalaecrdlg
{
  DIALOG_MESSAGE_BLOCK
  GrouP  qual_kind;
  
  DialoG source_qual;
  DialoG location;
  DialoG origin;
  DialoG strings;

  DialoG taxname_options;
  
  Int4   action_choice;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
  
} SourceQualAECRDlgData, PNTR SourceQualAECRDlgPtr;


static Boolean AreTaxNameOptionsRelevant (Int4 action_choice, SimpleAECRPtr sp)
{
  Boolean relevant = FALSE;
  if (sp == NULL) return FALSE;

  switch (action_choice)
  {
    case AECR_APPLY:
    case AECR_EDIT:
    case AECR_REMOVE:
      if (sp->field_list != NULL && sp->field_list->choice == 1)
      {
        relevant = TRUE;
      }
      break;
    case AECR_SWAP:
      if ((sp->field_list != NULL && sp->field_list->choice == 1) 
          || (sp->field_list_to != NULL && sp->field_list_to->choice == 1))
      {
        relevant = TRUE;
      }
      break;
    case AECR_CONVERT:
      if ((sp->field_list != NULL && sp->field_list->choice == 1 && !sp->leave_on_original) 
          || (sp->field_list_to != NULL && sp->field_list_to->choice == 1))
      {
        relevant = TRUE;
      }
      break;
    case AECR_PARSE:
      if ((sp->field_list != NULL && sp->field_list->choice == 1 && sp->remove_parsed)
          || (sp->field_list_to != NULL && sp->field_list_to->choice == 1))
      {
        relevant = TRUE;
      }
      break;
    default:
      relevant = FALSE;
      break;
  }
  return relevant;
}

static void ChangeQualKind (GrouP g)
{
  SourceQualAECRDlgPtr dlg;
  Int4                 qual_kind;
  SimpleAECRPtr        sp;

  dlg = (SourceQualAECRDlgPtr) GetObjectExtra (g);
  if (dlg == NULL)
  {
    return;
  }
  
  Hide (dlg->source_qual);
  Hide (dlg->location);
  Hide (dlg->origin);
  Hide (dlg->strings);
  Hide (dlg->taxname_options);
  qual_kind = GetValue (dlg->qual_kind);
  switch (qual_kind)
  {
    case QUALKIND_SOURCE_QUAL :
      Show (dlg->source_qual);
      break;
    case QUALKIND_LOCATION :
      Show (dlg->location);
      break;
    case QUALKIND_ORIGIN :
      Show (dlg->origin);
      break;
    case QUALKIND_STRINGS :
      Show (dlg->strings);
      Show (dlg->taxname_options);
      sp = DialogToPointer (dlg->strings);
      if (AreTaxNameOptionsRelevant (dlg->action_choice, sp))
      {
        Enable (dlg->taxname_options);
      }
      else
      {
        Disable (dlg->taxname_options);
      }
      sp = SimpleAECRFree (sp);
      break;
  }
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify)(dlg->change_userdata);
  }
}

static void ResetSourceQualAECRDlg (SourceQualAECRDlgPtr dlg)
{
  if (dlg == NULL)
  {
    return;
  }
  
  PointerToDialog (dlg->source_qual, NULL);
  PointerToDialog (dlg->location, NULL);
  PointerToDialog (dlg->origin, NULL);
  PointerToDialog (dlg->strings, NULL);
  ChangeQualKind (dlg->qual_kind);
  PointerToDialog (dlg->taxname_options, NULL);
}

static void ClearTextSourceQualAECRDlg (SourceQualAECRDlgPtr dlg)
{
  if (dlg == NULL)
  {
    return;
  }
  SendMessageToDialog (dlg->source_qual, AECR_VIB_MSG_CLEAR_TEXT);
  SendMessageToDialog (dlg->location, AECR_VIB_MSG_CLEAR_TEXT);
  SendMessageToDialog (dlg->origin, AECR_VIB_MSG_CLEAR_TEXT);
  SendMessageToDialog (dlg->strings, AECR_VIB_MSG_CLEAR_TEXT);
}

static void SourceQualAECRToDialog (DialoG d, Pointer userdata)
{
  SourceQualAECRDlgPtr dlg;
  SourceQualAECRPtr    data;
  
  dlg = (SourceQualAECRDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  ResetSourceQualAECRDlg (dlg);
  
  data = (SourceQualAECRPtr) userdata;
  if (data == NULL)
  {
    return;
  }
  PointerToDialog (dlg->source_qual, data->source_qual);
  PointerToDialog (dlg->location, data->location);
  PointerToDialog (dlg->origin, data->origin);
  PointerToDialog (dlg->strings, data->strings);
  
  if (data->qual_kind < QUALKIND_SOURCE_QUAL)
  {
    SetValue (dlg->qual_kind, QUALKIND_SOURCE_QUAL);
  }
  else if (data->qual_kind > QUALKIND_STRINGS 
           && (dlg->action_choice == AECR_EDIT || dlg->action_choice == AECR_SWAP))
  {
    SetValue (dlg->qual_kind, QUALKIND_STRINGS);
    PointerToDialog (dlg->taxname_options, data->taxname_options);
  }
  else if (data->qual_kind > QUALKIND_ORIGIN)
  {
    SetValue (dlg->qual_kind, QUALKIND_ORIGIN);
  }
  else
  {
    SetValue (dlg->qual_kind, data->qual_kind);
  }
  ChangeQualKind (dlg->qual_kind);  
}

static Pointer DialogToSourceQualAECR (DialoG d)
{
  SourceQualAECRDlgPtr dlg;
  SourceQualAECRPtr    data;
  
  dlg = (SourceQualAECRDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  data = (SourceQualAECRPtr) MemNew (sizeof (SourceQualAECRData));
  if (data != NULL)
  {
    data->taxname_options = NULL;
    data->qual_kind = GetValue (dlg->qual_kind);
    switch (data->qual_kind)
    {
      case QUALKIND_SOURCE_QUAL :
        data->source_qual = DialogToPointer (dlg->source_qual);
        break;
      case QUALKIND_LOCATION :
        data->location = DialogToPointer (dlg->location);
        break;
      case QUALKIND_ORIGIN :
        data->origin = DialogToPointer (dlg->origin);
        break;
      case QUALKIND_STRINGS :
        data->strings = DialogToPointer (dlg->strings);
        if (AreTaxNameOptionsRelevant (dlg->action_choice, data->strings))
        {
          data->taxname_options = DialogToPointer (dlg->taxname_options);
        }
        break;
    }
    
  }
  return data;
}

static void SourceQualAECRMessage (DialoG d, Int2 mssg)

{
  SourceQualAECRDlgPtr  dlg;

  dlg = (SourceQualAECRDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) 
    {
      case VIB_MSG_INIT :
        /* reset list */
        ResetSourceQualAECRDlg (dlg);
        break;
      case VIB_MSG_ENTER :
        Select (dlg->qual_kind);
        break;
      case AECR_VIB_MSG_SET_DEFAULT :
        SetValue (dlg->qual_kind, QUALKIND_SOURCE_QUAL);
        SendMessageToDialog (dlg->source_qual, AECR_VIB_MSG_SET_DEFAULT);
        SendMessageToDialog (dlg->location, AECR_VIB_MSG_SET_DEFAULT);
        SendMessageToDialog (dlg->origin, AECR_VIB_MSG_SET_DEFAULT);
        SendMessageToDialog (dlg->strings, AECR_VIB_MSG_SET_DEFAULT);
        PointerToDialog (dlg->taxname_options, NULL);
        ChangeQualKind (dlg->qual_kind);
        break;
      case AECR_VIB_MSG_CLEAR_TEXT :
        ClearTextSourceQualAECRDlg (dlg);
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestSourceQualAECR (DialoG d)
{
  SourceQualAECRDlgPtr dlg;
  ValNodePtr           total_err_list = NULL;
  Int4                 qual_kind;
  
  dlg = (SourceQualAECRDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return FALSE;
  }

  qual_kind = GetValue (dlg->qual_kind);
  switch (qual_kind)
  {
    case QUALKIND_SOURCE_QUAL :
      total_err_list = TestDialog (dlg->source_qual);
      break;
    case QUALKIND_LOCATION :
      total_err_list = TestDialog (dlg->location);
      break;
    case QUALKIND_ORIGIN :
      total_err_list = TestDialog (dlg->origin);
      break;
    case QUALKIND_STRINGS :
      total_err_list = TestDialog (dlg->strings);
      break;
    default:
      ValNodeAddPointer (&total_err_list, 0, "qual kind");
      break;
  }
  return total_err_list;
}

typedef struct applysourcequaldlg
{
  DIALOG_MESSAGE_BLOCK
  DialoG               field_list;
  DialoG               edit_apply;
  PopuP                true_false;
  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
  Uint2                entityID;
  
} ApplySourceQualDlgData, PNTR ApplySourceQualDlgPtr;

static void SourceQualChangeNotify (Pointer userdata)
{
  ApplySourceQualDlgPtr dlg;
  ValNodePtr            vnp;
  SourceQualDescPtr     sqdp;
  
  dlg = (ApplySourceQualDlgPtr) userdata;
  if (dlg == NULL)
  {
    return;
  }
  
  vnp = (ValNodePtr) DialogToPointer (dlg->field_list);
  if (vnp == NULL || vnp->data.ptrvalue == NULL)
  {
    Hide (dlg->edit_apply);
    Hide (dlg->true_false);
  }
  else 
  {
    sqdp = (SourceQualDescPtr) vnp->data.ptrvalue;
    if (IsNonTextModifier (sqdp->name))
    {
      Hide (dlg->true_false);
      Hide (dlg->edit_apply);
    }
    else
    {
      Show (dlg->edit_apply);
      Hide (dlg->true_false);
    }
  }
  vnp = ValNodeFreeData (vnp);
  
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}

static void ClearTextApplySourceQualDialog (ApplySourceQualDlgPtr dlg)
{
  if (dlg != NULL)
  {
    PointerToDialog (dlg->edit_apply, NULL);
  }
}

static void ResetApplySourceQualDialog (ApplySourceQualDlgPtr dlg)
{
  if (dlg == NULL)
  {
    return;
  }
  PointerToDialog (dlg->field_list, NULL);
  PointerToDialog (dlg->edit_apply, NULL);
  SetValue (dlg->true_false, 1);
  SourceQualChangeNotify (dlg);
}

static void ApplySourceQualToDialog (DialoG d, Pointer data)
{
  ApplySourceQualDlgPtr dlg;
  SimpleAECRPtr         sap;
  SourceQualDescPtr     sqdp;
  
  dlg = (ApplySourceQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  ResetApplySourceQualDialog (dlg);
  
  sap = (SimpleAECRPtr) data;
  if (sap != NULL)
  {
    PointerToDialog (dlg->field_list, sap->field_list);

    if (sap->field_list != NULL && sap->field_list->data.ptrvalue != NULL)
    {
      sqdp = (SourceQualDescPtr) sap->field_list->data.ptrvalue;
      if (IsNonTextModifier (sqdp->name))
      {
        if (sap->edit_apply == NULL 
            || StringHasNoText (sap->edit_apply->apply_txt))
        {
          SetValue (dlg->true_false, 2);
        }
        else
        {
          SetValue (dlg->true_false, 1);
        }
      }
      else
      {
        PointerToDialog (dlg->edit_apply, sap->edit_apply);
      }
    }
  }
  SourceQualChangeNotify (dlg);
}

static Pointer DialogToApplySourceQual (DialoG d)
{
  ApplySourceQualDlgPtr dlg;
  SimpleAECRPtr         sap;
  SourceQualDescPtr     sqdp;
  
  dlg = (ApplySourceQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  sap = (SimpleAECRPtr) MemNew (sizeof (SimpleAECRData));
  if (sap == NULL)
  {
    return NULL;
  }
  sap->field_list = DialogToPointer (dlg->field_list);
  sap->edit_apply = DialogToPointer (dlg->edit_apply);
  if (sap->field_list != NULL && sap->field_list->data.ptrvalue != NULL)
  {
    sqdp = (SourceQualDescPtr) sap->field_list->data.ptrvalue;
    if (IsNonTextModifier (sqdp->name))
    {
      sap->edit_apply->apply_txt = MemFree (sap->edit_apply->apply_txt);
      if (GetValue (dlg->true_false) == 1)
      {
        sap->edit_apply->apply_txt = StringSave ("TRUE");
      }
    }
  }
  
  sap->free_field_vn_proc = ValNodeSimpleDataFree;
  sap->free_subtype_vn_proc = NULL;
  return sap;  
}

static void ApplySourceQualMessage (DialoG d, Int2 mssg)

{
  ApplySourceQualDlgPtr  dlg;

  dlg = (ApplySourceQualDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) 
    {
      case VIB_MSG_INIT :
        /* reset list */
        ResetApplySourceQualDialog (dlg);
        break;
      case VIB_MSG_ENTER :
        Select (dlg->field_list);
        break;
      case AECR_VIB_MSG_SET_DEFAULT :
        SendMessageToDialog (dlg->field_list, AECR_VIB_MSG_SET_DEFAULT);
        ClearTextApplySourceQualDialog (dlg);
        break;
      case AECR_VIB_MSG_CLEAR_TEXT :
        ClearTextApplySourceQualDialog (dlg);
        break;
      case AECR_VIB_MSG_AUTOPOPULATE:
        SourceQualChangeNotify (dlg);
        break; 
      default :
        break;
    }
  }
}

static ValNodePtr TestApplySourceQual (DialoG d)
{
  ApplySourceQualDlgPtr dlg;
  ValNodePtr            err_list = NULL;
  ValNodePtr            total_err_list = NULL;
  ValNodePtr            field_list;
  SourceQualDescPtr     sqdp;
   
  dlg = (ApplySourceQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }

  err_list = TestDialog (dlg->field_list);
  total_err_list = ValNodeAppend (total_err_list, err_list);
      
  field_list = DialogToPointer (dlg->field_list);
  if (field_list != NULL && field_list->data.ptrvalue != NULL)
  {
    sqdp = (SourceQualDescPtr) field_list->data.ptrvalue;
    if (IsNonTextModifier (sqdp->name))
    {
      /* any true or false value is valid */
    }
    else
    {
      err_list = TestDialog (dlg->edit_apply);
      total_err_list = ValNodeAppend (total_err_list, err_list);
    }
  }
  field_list = ValNodeFreeData (field_list);
      
  return total_err_list;
}

static DialoG ApplySourceQualDialog
(GrouP                    h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Uint2                    entityID)
{
  ApplySourceQualDlgPtr dlg;  
  GrouP                 p, g1, g2;
  
  dlg = (ApplySourceQualDlgPtr) MemNew (sizeof (ApplySourceQualDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = ApplySourceQualToDialog;
  dlg->fromdialog = DialogToApplySourceQual;
  dlg->dialogmessage = ApplySourceQualMessage;
  dlg->testdialog = TestApplySourceQual;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->entityID = entityID;
  
  g1 = HiddenGroup (p, 1, 0, NULL);
  StaticPrompt (g1, "Qualifier", 0, dialogTextHeight, systemFont, 'l');

  dlg->field_list = SourceQualTypeSelectionDialog (g1, FALSE, 
                                                SourceQualChangeNotify,
                                                dlg);
  
  g2 = HiddenGroup (p, 0, 0, NULL);
  dlg->edit_apply = EditApplyDialog (g2, eEditApplyChoice_Apply, "New Text", NULL,
                                     change_notify, change_userdata);
  
  dlg->true_false = PopupList (g2, TRUE, NULL);
  PopupItem (dlg->true_false, "TRUE");
  PopupItem (dlg->true_false, "FALSE");
  SetValue (dlg->true_false, 1);
  
  Hide (dlg->edit_apply);
  Hide (dlg->true_false);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->edit_apply, (HANDLE) dlg->true_false, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) g2, NULL);
  
  SourceQualChangeNotify (dlg);
  
  return (DialoG) p;
}

static DialoG ConvertSourceQualDialog
(GrouP                    h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Uint2                    entityID)
{
  SimpleAECRDlgPtr dlg;
  GrouP            p, g1;
  
  dlg = (SimpleAECRDlgPtr) MemNew (sizeof (SimpleAECRDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = SimpleAECRToDialog;
  dlg->fromdialog = DialogToSimpleAECR;
  dlg->dialogmessage = SimpleAECRMessage;
  dlg->testdialog = TestSimpleAECR;
  dlg->action_choice = AECR_CONVERT;
  dlg->free_field_vn_proc = ValNodeSimpleDataFree;
  dlg->free_subtype_vn_proc = NULL;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->get_autopopulate_text = NULL;
  dlg->entityID = entityID;
  
  dlg->subtype_list = NULL;
  
  g1 = HiddenGroup (p, 2, 0, NULL);
  StaticPrompt (g1, "From", 0, dialogTextHeight, systemFont, 'l');
  StaticPrompt (g1, "To", 0, dialogTextHeight, systemFont, 'l');
    
  dlg->field_list = SourceQualTypeDiscSelectionDialog (g1, FALSE, 
                                          change_notify, 
                                          change_userdata);
  dlg->field_list_to = SourceQualTypeSelectionDialog (g1, FALSE, 
                                             change_notify, 
                                             change_userdata);
  dlg->leave_on_original = CheckBox (p, "Leave on original", NULL);
  dlg->strip_name_from_text = CheckBox (p, "Strip name from text", NULL);  
  
  dlg->edit_apply = NULL;

  AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) dlg->leave_on_original, (HANDLE) dlg->strip_name_from_text, NULL);
  return (DialoG) p;
}


static DialoG ConvertSourceLocationDialog
(GrouP                    h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Uint2                    entityID)
{
  SimpleAECRDlgPtr dlg;
  GrouP            p, g1;
  
  dlg = (SimpleAECRDlgPtr) MemNew (sizeof (SimpleAECRDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = SimpleAECRToDialog;
  dlg->fromdialog = DialogToSimpleAECR;
  dlg->dialogmessage = SimpleAECRMessage;
  dlg->testdialog = TestSimpleAECR;
  dlg->action_choice = AECR_CONVERT;
  dlg->free_field_vn_proc = ValNodeSimpleDataFree;
  dlg->free_subtype_vn_proc = NULL;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->get_autopopulate_text = NULL;
  dlg->entityID = entityID;
  
  dlg->subtype_list = NULL;
  
  g1 = HiddenGroup (p, 2, 0, NULL);
  StaticPrompt (g1, "From", 0, dialogTextHeight, systemFont, 'l');
  StaticPrompt (g1, "To", 0, dialogTextHeight, systemFont, 'l');
    
  dlg->field_list = SourceLocationDiscSelectionDialog (g1, FALSE, 
                                          change_notify, 
                                          change_userdata);
  dlg->field_list_to = SourceLocationSelectionDialog (g1, FALSE, 
                                             change_notify, 
                                             change_userdata);
  dlg->leave_on_original = CheckBox (p, "Leave on original", NULL);
  
  dlg->edit_apply = NULL;

  AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) dlg->leave_on_original, NULL);
  return (DialoG) p;
}


static void ChangeTaxonomyStrings (Pointer data)
{
  SourceQualAECRDlgPtr dlg;
  SimpleAECRPtr        sp;

  dlg = (SourceQualAECRDlgPtr) data;
  if (dlg == NULL) return;

  if (GetValue (dlg->qual_kind) == QUALKIND_STRINGS)
  {
    sp = DialogToPointer (dlg->strings);
    if (AreTaxNameOptionsRelevant (dlg->action_choice, sp))
    {
      Enable (dlg->taxname_options);
    }
    else
    {
      Disable (dlg->taxname_options);
    }
    sp = SimpleAECRFree (sp);
  }
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify)(dlg->change_userdata);
  }

}

static DialoG SourceQualAECRDialog 
(GrouP                    h,
 Int4                     action_choice, 
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Uint2                    entityID)
{

  SourceQualAECRDlgPtr dlg;
  GrouP                p, g1;
  
  dlg = (SourceQualAECRDlgPtr) MemNew (sizeof (SourceQualAECRDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = SourceQualAECRToDialog;
  dlg->fromdialog = DialogToSourceQualAECR;
  dlg->dialogmessage = SourceQualAECRMessage;
  dlg->testdialog = TestSourceQualAECR;
  dlg->action_choice = action_choice;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->qual_kind = HiddenGroup (p, 4, 0, ChangeQualKind);
  SetObjectExtra (dlg->qual_kind, dlg, NULL);
  RadioButton (dlg->qual_kind, "Qualifiers");
  RadioButton (dlg->qual_kind, "Taxonomy");
  if (action_choice != AECR_EDIT && action_choice != AECR_SWAP)
  {
    RadioButton (dlg->qual_kind, "Location");
    RadioButton (dlg->qual_kind, "Origin");
  }
  SetValue (dlg->qual_kind, QUALKIND_SOURCE_QUAL);
    
  g1 = HiddenGroup (p, 0, 0, NULL);
  /* Source Qualifiers dialog */
  if (action_choice == AECR_APPLY)
  {
    dlg->source_qual = ApplySourceQualDialog (g1, 
                                              change_notify, 
                                              change_userdata,
                                              entityID);
  }
  else if (action_choice == AECR_REMOVE)
  {
    dlg->source_qual = SimpleAECRDialog (g1, action_choice, 
                           change_notify, change_userdata,
                           ValNodeSimpleDataFree, NULL,
                           SourceQualTypeDiscSelectionDialog, NULL, NULL,
                           "Qualifier", NULL, TRUE, entityID);    
  }
  else if (action_choice == AECR_CONVERT)
  {
    dlg->source_qual = ConvertSourceQualDialog (g1,
                                                change_notify, 
                                                change_userdata,
                                                entityID);
  }
  else if (action_choice == AECR_SWAP)
  {
    dlg->source_qual = SimpleAECRDialogEx (g1, action_choice, 
                           change_notify, change_userdata,
                           ValNodeSimpleDataFree, NULL,
                           SourceQualTypeSelectionDialog, NULL, NULL,
                           "Qualifier", NULL, TRUE, TRUE, entityID);    
  }
  else
  {
    dlg->source_qual = SimpleAECRDialog (g1, action_choice, 
                           change_notify, change_userdata,
                           ValNodeSimpleDataFree, NULL,
                           SourceQualTypeSelectionDialog, NULL, NULL,
                           "Qualifier", NULL, TRUE, entityID);    
  }
  
  /* Taxonomy Strings dialog */
  dlg->strings = SimpleAECRDialog (g1, action_choice, ChangeTaxonomyStrings, dlg,
                           ValNodeSimpleDataFree, NULL,
                           SourceStringSelectionDialog, NULL, NULL,
                           "Taxonomy", NULL, TRUE, entityID);
                          
  if (action_choice != AECR_EDIT && action_choice != AECR_SWAP)
  {
    /* Location Dialog */
    if (action_choice == AECR_APPLY)
    {
      dlg->location = SimpleAECRDialog (g1, action_choice,
                                        change_notify, change_userdata,
                                        ValNodeSimpleDataFree, NULL,
                                        SourceLocationSelectionDialog, NULL,
                                        NULL, "Location", NULL, FALSE, entityID);
    }
    else if (action_choice == AECR_CONVERT)
    {
      dlg->location = ConvertSourceLocationDialog (g1,
                                                    change_notify,
                                                    change_userdata, 
                                                    entityID);
    }
    else
    {
      dlg->location = SimpleAECRDialog (g1, action_choice,
                                        change_notify, change_userdata,
                                        ValNodeSimpleDataFree, NULL,
                                        SourceLocationDiscSelectionDialog, 
                                        NULL, NULL, "Location", NULL, FALSE,
                                        entityID);
    }

    /* Origin Dialog */
    dlg->origin = SimpleAECRDialog (g1, action_choice, change_notify, change_userdata,
                           ValNodeSimpleDataFree, NULL,
                           SourceOriginSelectionDialog, NULL, NULL,
                           "Origin", NULL, FALSE, entityID);
  }

  dlg->taxname_options = TaxnameOptionsDialog (p, change_notify, change_userdata);

  Disable (dlg->taxname_options);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->source_qual,
                              (HANDLE) dlg->location,
                              (HANDLE) dlg->origin,
                              (HANDLE) dlg->strings,
                              NULL);
                              
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->qual_kind,
                              (HANDLE) g1,
                              (HANDLE) dlg->taxname_options,
                              NULL);
                           
  ChangeQualKind (dlg->qual_kind);                         
  return (DialoG) p;
}


static Boolean BioSourceHasMultipleQuals (BioSourcePtr biop, SourceQualDescPtr sqdp, FilterSetPtr fsp)
{
  OrgModPtr    mod;
  SubSourcePtr ssp;
  Boolean      found_one = FALSE;

  if (biop == NULL || sqdp == NULL) return FALSE;

  if (sqdp->isOrgMod)
  {
    if (biop->org != NULL && biop->org->orgname != NULL)
    {
      mod = biop->org->orgname->mod;
      while (mod != NULL)
      {
        if (mod->subtype == sqdp->subtype
            && !StringHasNoText (mod->subname)
            && OrgModSpecialMatch (mod, fsp))
        {
          if (found_one) {
            return TRUE;
          } else {
            found_one = TRUE;
          }
        }
        mod = mod->next;
      }
    }
  }
  else
  {
    ssp = biop->subtype;
    while (ssp != NULL)
    {
      if (ssp->subtype == sqdp->subtype
          && !StringHasNoText (ssp->name)
          && SubSrcSpecialMatch (ssp, fsp))
      {
        if (found_one) {
          return TRUE;
        } else {
          found_one = TRUE;
        }
      }
      ssp = ssp->next;
    }
  }
  return FALSE;
}

static Boolean CheckBioSrcListForMultipleQuals (ValNodePtr biosrc_list, SourceQualDescPtr sqdp, FilterSetPtr fsp)
{
  ValNodePtr   vnp;
  SeqFeatPtr   sfp;
  SeqDescrPtr  sdp;
  Boolean      has_multiple = FALSE;

  if (sqdp == NULL) return FALSE;

  for (vnp = biosrc_list; vnp != NULL && !has_multiple; vnp = vnp->next)
  {
    if (vnp->data.ptrvalue == NULL) continue;
    if (vnp->choice == OBJ_SEQFEAT)
    {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      if (sfp->data.choice == SEQFEAT_BIOSRC)
      {
        has_multiple = BioSourceHasMultipleQuals ((BioSourcePtr)sfp->data.value.ptrvalue, sqdp, fsp);
      }
    }
    else if (vnp->choice == OBJ_SEQDESC)
    {
      sdp = (SeqDescrPtr) vnp->data.ptrvalue;
      if (sdp->choice == Seq_descr_source)
      {
        has_multiple = BioSourceHasMultipleQuals ((BioSourcePtr)sdp->data.ptrvalue, sqdp, fsp);
      }
    }
  }
  return has_multiple;
}


typedef struct taxnameoptionsdlg
{
  DIALOG_MESSAGE_BLOCK
  ButtoN               remove_taxref;
  ButtoN               remove_common;
  ButtoN               remove_old_name;
  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
} TaxnameOptionsDlgData, PNTR TaxnameOptionsDlgPtr;

static void TaxnameOptionsToDialog (DialoG d, Pointer data)
{
  TaxnameOptionsDlgPtr dlg;
  TaxnameOptionsPtr    top;

  dlg = (TaxnameOptionsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  top = (TaxnameOptionsPtr) data;
  if (top == NULL)
  {
    SetStatus (dlg->remove_taxref, TRUE);
    SetStatus (dlg->remove_common, TRUE);
    SetStatus (dlg->remove_old_name, TRUE);
  }
  else
  {
    SetStatus (dlg->remove_taxref, top->remove_taxref);
    SetStatus (dlg->remove_common, top->remove_common);
    SetStatus (dlg->remove_old_name, top->remove_old_name);
  }
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}

static Pointer TaxnameOptionsFromDialog (DialoG d)
{
  TaxnameOptionsDlgPtr dlg;
  TaxnameOptionsPtr    top;

  dlg = (TaxnameOptionsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  top = (TaxnameOptionsPtr) MemNew (sizeof (TaxnameOptionsData));
  top->remove_taxref = GetStatus (dlg->remove_taxref);
  top->remove_common = GetStatus (dlg->remove_common);
  top->remove_old_name = GetStatus (dlg->remove_old_name);
  return top;
}


static void ChangeTaxnameOptions (ButtoN b)
{
  TaxnameOptionsDlgPtr dlg;

  dlg = (TaxnameOptionsDlgPtr) GetObjectExtra (b);
  if (dlg != NULL && dlg->change_notify != NULL)
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


extern DialoG TaxnameOptionsDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  TaxnameOptionsDlgPtr dlg;
  GrouP                p;
  
  dlg = (TaxnameOptionsDlgPtr) MemNew (sizeof (TaxnameOptionsDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, 3, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = TaxnameOptionsToDialog;
  dlg->fromdialog = TaxnameOptionsFromDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->remove_taxref = CheckBox (p, "Remove taxonomy xref", ChangeTaxnameOptions);
  SetObjectExtra (dlg->remove_taxref, dlg, NULL);
  SetStatus (dlg->remove_taxref, TRUE);
  dlg->remove_common = CheckBox (p, "Remove GenBank common name", ChangeTaxnameOptions);
  SetObjectExtra (dlg->remove_common, dlg, NULL);
  SetStatus (dlg->remove_common, TRUE);
  dlg->remove_old_name = CheckBox (p, "Remove old_name qualifier", ChangeTaxnameOptions);
  SetObjectExtra (dlg->remove_old_name, dlg, NULL);
  SetStatus (dlg->remove_old_name, TRUE);
  return (DialoG) p;
}

extern void ApplyTaxnameOptionsToBioSource (BioSourcePtr biop, TaxnameOptionsPtr top)
{
  if (biop == NULL || top == NULL) return;

  if (top->remove_taxref)
  {
    RemoveTaxRef (biop->org);
  }
  if (top->remove_old_name)
  {
    RemoveOldName (biop->org);
  }
  if (top->remove_common) 
  {
    if (biop->org != NULL)
    {
      biop->org->common = MemFree (biop->org->common);
    }
  }
}

static void HandleTaxnameOptionsFeature (SeqFeatPtr sfp, TaxnameOptionsPtr top)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || top == NULL)
  {
    return;
  }
  ApplyTaxnameOptionsToBioSource (sfp->data.value.ptrvalue, top);
}

static void HandleTaxnameOptionsDescriptorCallback (SeqDescrPtr sdp, TaxnameOptionsPtr top)
{
  if (sdp == NULL || sdp->choice != Seq_descr_source || top == NULL)
  {
    return;
  }
  ApplyTaxnameOptionsToBioSource (sdp->data.ptrvalue, top);
}


static void HandleTaxnameOptionsForBioSourceList (ValNodePtr biosrc_list, TaxnameOptionsPtr top)
{
  ValNodePtr vnp;

  if (top == NULL) return;
  for (vnp = biosrc_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == OBJ_SEQFEAT) 
    {
      HandleTaxnameOptionsFeature (vnp->data.ptrvalue, top);
    }
    else if (vnp->choice == OBJ_SEQDESC)
    {
      HandleTaxnameOptionsDescriptorCallback (vnp->data.ptrvalue, top);
    }
  }  
}

static void CollectBioSourceFeatureFieldCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_SEQFEAT, sfp);  
}

static void CollectBioSourceDescriptorFieldCallback (SeqDescrPtr sdp, Pointer userdata, FilterSetPtr fsp)
{  
  if (sdp == NULL || userdata == NULL)
  {
    return;
  }
  
  ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_SEQDESC, sdp);  
}


static void RemoveSourceQualEditNonMatches (ValNodePtr PNTR biosrc_list, ApplyValuePtr avp)
{
  ValNodePtr vnp_prev = NULL, vnp_next, vnp_this;
  CharPtr    val1, val2;
  Boolean    do_remove;

  if (biosrc_list == NULL || avp == NULL || StringHasNoText (avp->text_to_replace))
  {
    return;
  }

  vnp_this = *biosrc_list;
  while (vnp_this != NULL)
  {
    vnp_next = vnp_this->next;
    do_remove = FALSE;
    if (vnp_this->choice == OBJ_SEQFEAT)
    {
      val1 = GetSourceQualFeatureString (vnp_this->data.ptrvalue, avp->field_list, NULL);
    }
    else if (vnp_this->choice == OBJ_SEQDESC)
    {
      val1 = GetSourceQualDescrString (vnp_this->data.ptrvalue, avp->field_list, NULL);
    } 
    else
    {
      val1 = NULL;
    }
    if (val1 == NULL)
    {
      do_remove = TRUE;
    }
    else
    {
      val2 = HandleApplyValue (StringSave (val1), avp);
      if (StringCmp (val1, val2) == 0)
      {
        do_remove = TRUE;
      }
      val1 = MemFree (val1);
      val2 = MemFree (val2);
    }
    if (do_remove)
    {
      if (vnp_prev == NULL)
      {
        *biosrc_list = vnp_this->next;
      }
      else
      {
        vnp_prev->next = vnp_this->next;
      }
      vnp_this->next = NULL;
      vnp_this = ValNodeFree (vnp_this);
    }
    else
    {
      vnp_prev = vnp_this;
    }
    vnp_this = vnp_next;
  }
}


static Boolean 
SourceQualQualApplyEditConvertRemoveAction
(ApplyEditConvertRemovePtr   ap,
 SimpleAECRPtr               sp,
 Int4                        action_choice,
 SetFeatureFieldString       feature_apply_action,
 SetDescriptorFieldString    descriptor_apply_action,
 RemoveFeatureFieldString    feature_remove_action,
 RemoveDescriptorFieldString descriptor_remove_action,
 CopyValNodeDataProc         copy_vn_proc,
 GetFeatureFieldString       fieldstring_func,
 GetDescriptorFieldString    descrstring_func,
 Boolean                     non_text,
 NameFromValNodeProc         name_field_func,
 TaxnameOptionsPtr           top)

{
  SeqEntryPtr               sep;
  ApplyValueData            avd;
  ConvertFieldData          cfd;
  FilterSetPtr              fsp;
  GetSamplePtr              gsp;
  Boolean                   rval = TRUE;
  ValNodePtr                vnp;
  CharPtr                   field_name = NULL;
  MsgAnswer                 ans;
  ValNodePtr                biosrc_list = NULL;
  SourceQualDescPtr         multi_sqdp = NULL;
  
  sep = GetTopSeqEntryForEntityID (ap->input_entityID);
  if (sp == NULL
      || sp->field_list == NULL
      || feature_apply_action == NULL
      || descriptor_apply_action == NULL
      || feature_remove_action == NULL
      || descriptor_remove_action == NULL
      || ((action_choice == AECR_CONVERT || action_choice == AECR_SWAP || action_choice == AECR_PARSE) 
           && sp->field_list == NULL))
  {
    return FALSE;
  }
  
  WatchCursor ();
  Update ();
  
  avd.where_to_replace = EditApplyFindLocation_anywhere;

  fsp = (FilterSetPtr) DialogToPointer (ap->constraints);

  /* collect list of BioSources that will be acted upon */
  OperateOnSeqEntryConstrainedObjects (sep, fsp, 
                                       CollectBioSourceFeatureFieldCallback,
                                       CollectBioSourceDescriptorFieldCallback,
                                       SEQFEAT_BIOSRC, 0, Seq_descr_source,
                                       &biosrc_list);    
  
  if (action_choice == AECR_CONVERT || action_choice == AECR_SWAP || action_choice == AECR_PARSE)
  {
    /* check for multiple sources */
    if ((sp->field_list->choice == 0
         && (multi_sqdp = sp->field_list->data.ptrvalue) != NULL
         && CheckBioSrcListForMultipleQuals(biosrc_list, multi_sqdp, fsp))
        || (action_choice == AECR_SWAP && sp->field_list_to->choice == 0
            && (multi_sqdp = sp->field_list_to->data.ptrvalue) != NULL
            && CheckBioSrcListForMultipleQuals(biosrc_list, multi_sqdp, fsp))) {
      if (name_field_func != NULL) {  
        if (action_choice == AECR_SWAP && multi_sqdp == sp->field_list_to->data.ptrvalue) {
          field_name = name_field_func(sp->field_list_to);
        } else {
          field_name = name_field_func(sp->field_list);
        }
      }
      if (field_name != NULL) {
        ans = Message (MSG_OKC, "One or more BioSources has multiple source quals.  These will be combined.  Do you want to continue?");
      } else {
        ans = Message (MSG_OKC, "One or more BioSources has multiple %s source quals.  These will be combined.  Do you want to continue?", field_name);
      }
      field_name = MemFree (field_name);
      if (ans == ANS_CANCEL) {
        FilterSetFree (fsp);
        ArrowCursor();
        Update();
        return FALSE;
      }
    }
    cfd.src_field_list = sp->field_list;
    cfd.dst_field_list = sp->field_list_to;
    cfd.etp = NULL;
    cfd.get_str_func = fieldstring_func;
    cfd.set_str_func = feature_apply_action;
    cfd.remove_str_func = feature_remove_action;
    cfd.get_d_str_func = descrstring_func;
    cfd.set_d_str_func = descriptor_apply_action;
    cfd.remove_d_str_func = descriptor_remove_action;
    cfd.fsp = fsp;
    cfd.strip_name_from_text = sp->strip_name_from_text;
    cfd.remove_parsed = sp->remove_parsed;
    cfd.name_field_func = name_field_func;
    cfd.text_portion = sp->text_portion;
    if (action_choice == AECR_CONVERT)
    {
      if (sp->leave_on_original)
      {
        cfd.convert_type = CONVERT_TYPE_COPY;
      }
      else
      {
        cfd.convert_type = CONVERT_TYPE_MOVE;
      }
    }
    else if (action_choice == AECR_SWAP)
    {
      cfd.convert_type = CONVERT_TYPE_SWAP;
    }
    else if (action_choice == AECR_PARSE)
    {
      cfd.convert_type = CONVERT_TYPE_PARSE;
    }
    if ((cfd.convert_type == CONVERT_TYPE_MOVE || cfd.convert_type == CONVERT_TYPE_COPY || cfd.convert_type == CONVERT_TYPE_PARSE)
        && ! non_text)
    {
      gsp = CheckForConversionExistingTextInSeqEntry (sep, &cfd, fsp, 
                                                      SEQFEAT_BIOSRC, 0,
                                                      Seq_descr_source);
      cfd.etp = GetExistingTextHandlerInfo (gsp == NULL ? 0 : gsp->num_found, non_text);
      gsp = GetSampleFree (gsp);
      if (cfd.etp != NULL && cfd.etp->existing_text_choice == eExistingTextChoiceCancel)
      {
        rval = FALSE;
      }
    }
    else
    {
      cfd.etp = NULL;
    }  

    if (rval)
    {
      if (non_text)
      {
        OperateOnSeqEntryConstrainedObjects (sep, fsp, 
                                             ConvertNonTextFeatureFieldCallback,
                                             ConvertNonTextDescriptorFieldCallback,
                                             SEQFEAT_BIOSRC, 0, Seq_descr_source,
                                             &cfd);    
      }
      else
      {
        OperateOnSeqEntryConstrainedObjects (sep, fsp, 
                                             ConvertFeatureFieldCallback,
                                             ConvertDescriptorFieldCallback,
                                             SEQFEAT_BIOSRC, 0, Seq_descr_source,
                                             &cfd);    
      }
      
    }
    cfd.etp = MemFree (cfd.etp);
  }
  else
  {
    avd.field_list = sp->field_list;
  
    /* get handling for existing text */
    if (action_choice == AECR_APPLY)
    {
      gsp = CheckForExistingText (ap->input_entityID,
                                  avd.field_list,
                                  fieldstring_func,
                                  descrstring_func,
                                  ValNodeSimpleDataFree,
                                  copy_vn_proc,
                                  fsp, SEQFEAT_BIOSRC, 0, Seq_descr_source);
      avd.etp = GetExistingTextHandlerInfo (gsp == NULL ? 0 : gsp->num_found, non_text);
      gsp = GetSampleFree (gsp);
    
      if (avd.etp != NULL 
          && avd.etp->existing_text_choice == eExistingTextChoiceCancel)
      {
        rval = FALSE;
      }
    }
    else
    {
      avd.etp = NULL;
    }
  
    if (rval)
    {
      AddEditApplyDataToApplyValue (action_choice, sp->edit_apply, &avd);

      if (action_choice == AECR_EDIT)
      {
        /* remove items that don't contain text that will be replaced */
        RemoveSourceQualEditNonMatches (&biosrc_list, &avd);
      }

   
      if (action_choice == AECR_EDIT || action_choice == AECR_APPLY)
      {
        OperateOnSeqEntryConstrainedObjects (sep, fsp, 
                                             feature_apply_action,
                                             descriptor_apply_action,
                                             SEQFEAT_BIOSRC, 0, Seq_descr_source, &avd);  
      }
      else if (action_choice == AECR_REMOVE)
      {
        for (vnp = sp->field_list; vnp != NULL; vnp = vnp->next)
        {
          avd.field_list = vnp;
          OperateOnSeqEntryConstrainedObjects (sep, fsp, 
                                               feature_remove_action,
                                               descriptor_remove_action,
                                               SEQFEAT_BIOSRC, 0, Seq_descr_source, &avd);
        }
      }

      avd.text_to_replace = MemFree (avd.text_to_replace);
      avd.new_text = MemFree (avd.new_text);
    }
    avd.etp = MemFree (avd.etp);
  }

  FilterSetFree (fsp);
  
  if (rval)
  {
    /* Handle taxname options */
    HandleTaxnameOptionsForBioSourceList (biosrc_list, top);

    ObjMgrSetDirtyFlag (ap->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, ap->input_entityID, 0, 0);
  }
  biosrc_list = ValNodeFree (biosrc_list);
  ArrowCursor ();
  Update (); 
  
  return rval;
}

static Boolean SourceQualApplyEditConvertRemoveAction (Pointer userdata)
{
  ApplyEditConvertRemovePtr ap;  
  SourceQualAECRPtr         sp;
  Int4                      action_choice;
  Boolean                   rval = FALSE;
  
  ap = (ApplyEditConvertRemovePtr) userdata;
  if (ap == NULL)
  {
    return FALSE;
  }
  
  if (ap->crippled)
  {
    action_choice = ap->crippled_action;
  }
  else
  {
    action_choice = GetValue (ap->action_popup);
  }
  if (action_choice < 1 || action_choice > NUM_AECR)
  {
    return FALSE;
  }
  sp = (SourceQualAECRPtr) DialogToPointer (ap->aecr_pages[action_choice - 1]);
  
  switch (sp->qual_kind)
  {
    case QUALKIND_SOURCE_QUAL :
      rval = SourceQualQualApplyEditConvertRemoveAction (ap, sp->source_qual, 
                                                         action_choice,
                                       ApplySourceQualFeatureCallback,
                                       ApplySourceQualDescriptorCallback,
                                       RemoveSourceQualFromSourceFeature,
                                       RemoveSourceQualFromSourceDescriptor,
                                       SourceQualValNodeDataCopy,
                                       GetSourceQualFeatureString,
                                       GetSourceQualDescrString,
                                       FALSE,
                                       SourceQualValNodeName,
                                       NULL);
      break;
    case QUALKIND_STRINGS :
      rval = SourceQualQualApplyEditConvertRemoveAction (ap, sp->strings, 
                                                         action_choice,
                                       ApplyStringToSourceFeature,
                                       ApplyStringToSourceDescriptor,
                                       RemoveStringFromSourceFeature,
                                       RemoveStringFromSourceDescriptor,
                                       ValNodeStringCopy,
                                       GetStringFromSourceFeature,
                                       GetStringFromSourceDescriptor,
                                       FALSE,
                                       NULL,
                                       sp->taxname_options);
      break;
    case QUALKIND_LOCATION :
      rval = SourceQualQualApplyEditConvertRemoveAction (ap, sp->location, 
                                                         action_choice,
                                       ApplyLocationToSourceFeature,
                                       ApplyLocationToSourceDescriptor,
                                       RemoveLocationFromSourceFeature,
                                       RemoveLocationFromSourceDescriptor,
                                       ValNodeStringCopy,
                                       GetLocationFromSourceFeature,
                                       GetLocationFromSourceDescriptor,
                                       TRUE,
                                       ValNodeStringName,
                                       NULL);
      break;      
    case QUALKIND_ORIGIN :
      rval = SourceQualQualApplyEditConvertRemoveAction (ap, sp->origin, 
                                                         action_choice,
                                       ApplyOriginToSourceFeature,
                                       ApplyOriginToSourceDescriptor,
                                       RemoveOriginFromSourceFeature,
                                       RemoveOriginFromSourceDescriptor,
                                       ValNodeStringCopy,
                                       GetOriginFromSourceFeature,
                                       GetOriginFromSourceDescriptor,
                                       TRUE,
                                       ValNodeStringName,
                                       NULL);
      break;      
  }
  sp = SourceQualAECRFree (sp);
  return rval;
}

static void ChangeSourceQualAction (DialoG prev_page, DialoG new_page)
{
  SourceQualAECRPtr         prev_sp, new_sp;
  
  if (prev_page == NULL || new_page == NULL)
  {
    return;
  }
  
  prev_sp = (SourceQualAECRPtr) DialogToPointer (prev_page);
  new_sp = (SourceQualAECRPtr) DialogToPointer (new_page);
  if (prev_sp != NULL && new_sp != NULL)
  {
    new_sp->qual_kind = prev_sp->qual_kind;
    PointerToDialog (new_page, new_sp);
  }
  new_sp = SourceQualAECRFree (new_sp);
  prev_sp = SourceQualAECRFree (prev_sp);
}

static DialoG GetSourceQualSample (GrouP h, Uint2 entityID)
{
  DialoG                d;
  SetSampleData         ssd;

  d = SampleDialog (h);
  ssd.fieldstring_func = GetSourceQualFeatureString;
  ssd.descrstring_func = GetSourceQualDescrString;
  ssd.entityID = entityID;
  ssd.field_list = GetSourceQualDescList (TRUE, TRUE, TRUE, FALSE);
  ValNodeAddPointer (&(ssd.field_list), 1, StringSave ("Scientific Name"));
  ValNodeAddPointer (&(ssd.field_list), 1, StringSave ("Location"));
  
  ssd.fsp = NULL;
  ssd.free_vn_proc = ValNodeSimpleDataFree;
  ssd.copy_vn_proc = SourceQualValNodeDataCopy;
  ssd.match_vn_proc = SourceQualValNodeMatch;
  ssd.label_vn_proc = SourceQualValNodeName;  
  PointerToDialog (d, &ssd);
  return d;
}

static void 
SourceQualApplyEditConvertRemove 
(IteM    i,
 Int4    first_action_choice,
 Boolean crippled)
{
  ApplyEditConvertRemoveCombo (i, first_action_choice, crippled,
                               "Source Qual",
                               SourceQualAECRDialog, 
                               FALSE, GetSourceQualSample,
                               FALSE, TRUE, FALSE, FALSE, FALSE, NULL,
                               SourceQualApplyEditConvertRemoveAction, NULL,
                               ChangeSourceQualAction, NULL);
}

extern void ApplySourceQual (IteM i)
{
  SourceQualApplyEditConvertRemove (i, AECR_APPLY, FALSE);
}

extern void PublicApplySourceQual (IteM i)
{
  SourceQualApplyEditConvertRemove (i, AECR_APPLY, TRUE);
}

extern void EditSourceQual (IteM i)
{
  SourceQualApplyEditConvertRemove (i, AECR_EDIT, FALSE);
}

extern void PublicEditSourceQual (IteM i)
{
  SourceQualApplyEditConvertRemove (i, AECR_EDIT, TRUE);
}

extern void ConvertSourceQual (IteM i)
{
  SourceQualApplyEditConvertRemove (i, AECR_CONVERT, FALSE);
}

extern void SwapSourceQual (IteM i)
{
  SourceQualApplyEditConvertRemove (i, AECR_SWAP, FALSE);
}

extern void RemoveSourceQual (IteM i)
{
  SourceQualApplyEditConvertRemove (i, AECR_REMOVE, FALSE);
}

static CharPtr rna_sample_field_list [] = { "Product", 
                                            "Comment" };
static Int2    num_rna_sample_fields = sizeof (rna_sample_field_list) / sizeof (CharPtr);


static CharPtr GetRNASampleFieldString (SeqFeatPtr sfp, ValNodePtr vnp, FilterSetPtr fsp)
{
  RnaRefPtr  rrp;
  CharPtr    str = NULL;
  Int2       rna_subtype, rna_field_type;
  GeneRefPtr grp;
  SeqFeatPtr gene;
  ValNode    vn;
  SeqMgrFeatContext    fcontext;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || vnp == NULL)
  {
    return NULL;
  }
  
  rrp = sfp->data.value.ptrvalue;
  
  while (vnp != NULL && str == NULL)
  {
    rna_subtype = (vnp->data.intvalue - 1) / (num_rna_sample_fields + num_gene_fields);
    rna_field_type = (vnp->data.intvalue - 1) % (num_rna_sample_fields + num_gene_fields);
    switch (rna_subtype)
    {
      case 0:
        if (rrp->type != 255 
            || (rrp->ext.choice == 1
                && (StringCmp (rrp->ext.value.ptrvalue, "ncRNA") == 0
                    || StringCmp (rrp->ext.value.ptrvalue, "tmRNA") == 0)))
        {
          vnp = vnp->next;
          continue;
        }
        break;
      case 1:
      case 2:
      case 3:
      case 4:
        if (rna_subtype != rrp->type)
        {
          vnp = vnp->next;
          continue;
        }
        break;
      case 5:
        if (rrp->type != 255 
            || rrp->ext.choice != 1
            || StringCmp (rrp->ext.value.ptrvalue, "ncRNA") != 0)
        {
          vnp = vnp->next;
          continue;
        }
        break;
      case 6:
        if (rrp->type != 255 
            || rrp->ext.choice != 1
            || StringCmp (rrp->ext.value.ptrvalue, "tmRNA") != 0)
        {
          vnp = vnp->next;
          continue;
        }
        break;
    }

    switch (rna_field_type)
    {
      case 0 :
        str = GetRNAProductString (sfp, NULL);
        break;
      case 1 :
        if (!StringHasNoText (sfp->comment))
        {
          str = StringSave (sfp->comment);
        }
        break;
      default :
        vn.next = NULL;
        vn.choice = 0;
        vn.data.intvalue = rna_field_type - num_rna_sample_fields + 1;
        grp = SeqMgrGetGeneXref (sfp);
        if (grp == NULL)
        {
          gene = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
          str = GetGeneFieldString (gene, &vn, NULL);
        }
        else
        {
          str = GetGeneFieldString (sfp, &vn, NULL);     
        }
        break;
    }
    vnp = vnp->next;
  }
  return str;
}

static CharPtr GetRNAQualSampleName (ValNodePtr vnp)
{
  CharPtr label = NULL;
  Int2    rna_subtype, rna_field_type;
  
  if (vnp == NULL)
  {
    return NULL;
  }
  rna_subtype = (vnp->data.intvalue - 1) / (num_rna_sample_fields + num_gene_fields);
  rna_field_type = (vnp->data.intvalue - 1) % (num_rna_sample_fields + num_gene_fields);
  
  if (rna_field_type < num_rna_sample_fields)
  {
    label = (CharPtr) MemNew ((StringLen (rna_subtype_list [rna_subtype])
                              + StringLen (rna_sample_field_list[rna_field_type])
                              + 2) * sizeof (Char));
    if (label != NULL)
    {
      sprintf (label, "%s %s", rna_subtype_list [rna_subtype],
               rna_sample_field_list[rna_field_type]);
    }
  }
  else if (rna_field_type < num_rna_sample_fields + num_gene_fields)
  {
    label = (CharPtr) MemNew ((StringLen (rna_subtype_list [rna_subtype])
                              + StringLen (gene_field_list[rna_field_type - num_rna_sample_fields])
                              + 7) * sizeof (Char));
    if (label != NULL)
    {
      sprintf (label, "%s Gene %s", rna_subtype_list [rna_subtype],
               gene_field_list[rna_field_type - num_rna_sample_fields]);
    }
  }
  return label;
}

static DialoG GetRNAQualSample (GrouP h, Uint2 entityID)
{
  DialoG        d;
  SetSampleData ssd;
  Int2          rna_subtype, field_type;
  Int4          list_index;      

  d = SampleDialog (h);
  ssd.fieldstring_func = GetRNASampleFieldString;
  ssd.descrstring_func = NULL;
  ssd.entityID = entityID;
  
  /* construct RNA field list for separate types of RNAs */
  ssd.free_vn_proc = NULL;
  ssd.copy_vn_proc = IntValNodeCopy;
  ssd.match_vn_proc = IntValNodeMatch;
  ssd.label_vn_proc = GetRNAQualSampleName;
  
  ssd.field_list = NULL;
  for (rna_subtype = 0; rna_subtype < num_rna_subtypes; rna_subtype++)
  {
    for (field_type = 0; field_type < num_rna_sample_fields; field_type++)
    {
      list_index = (rna_subtype * (num_rna_sample_fields + num_gene_fields)) + field_type;
      ValNodeAddInt (&ssd.field_list, 0, list_index + 1);
    }
    for (field_type = 0; field_type < num_gene_fields; field_type ++)
    {
      list_index = (rna_subtype * (num_rna_sample_fields + num_gene_fields))
                    + field_type + num_rna_sample_fields;
      ValNodeAddInt (&ssd.field_list, 0, list_index + 1);
    }
  }
  ssd.fsp = NULL;
  PointerToDialog (d, &ssd);
  
  /* now free up field list (SampleDialog maintains its own copy */
  ssd.field_list = ValNodeFree (ssd.field_list);
  
  return d;
}

typedef struct aecrrnaqual
{
  ValNodePtr   field_list;
  ValNodePtr   field_list_to;
  RnaTypePtr   rna_type;
  EditApplyPtr edit_apply;
  ValNodePtr   codons;
  SeqLocPtr    anticodon;
  TextPortionXPtr text_portion;
  Boolean        leave_on_original;
  Boolean        remove_parsed;
} AECRRNAQualData, PNTR AECRRNAQualPtr;

static AECRRNAQualPtr AECRRNAQualFree (AECRRNAQualPtr ap)
{
  if (ap != NULL)
  {
    ap->field_list = ValNodeFree (ap->field_list);
    ap->field_list_to = ValNodeFree (ap->field_list_to);
    ap->rna_type = RnaTypeFree (ap->rna_type);
    ap->edit_apply = EditApplyFree (ap->edit_apply);
    ap->codons = ValNodeFreeData (ap->codons);
    ap = MemFree (ap);
  }
  return ap;
}

typedef struct aecrrnaqualdlg
{
  DIALOG_MESSAGE_BLOCK
  DialoG               field_list;
  DialoG               field_list_to;
  DialoG               rna_type;
  DialoG               edit_apply;
  DialoG               codons;
  DialoG               ncrna_class;
  DialoG               text_portion;
  ButtoN               leave_on_original;
  ButtoN               remove_parsed;

  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
  Uint2                entityID;
  
  /* The following is used to create the anticodon location dialog */
  FeatureForm          ff;
  DialoG               anticodon;
} AECRRNAQualDlgData, PNTR AECRRNAQualDlgPtr;

static void ClearTextAECRRNAQualDlg (AECRRNAQualDlgPtr dlg)
{
  if (dlg != NULL)
  {
    PointerToDialog (dlg->edit_apply, NULL);
    PointerToDialog (dlg->codons, NULL);
    PointerToDialog (dlg->ncrna_class, NULL);
  }
}

static void ResetAECRRNAQualDlg (AECRRNAQualDlgPtr dlg)
{
  if (dlg != NULL)
  {
    PointerToDialog (dlg->field_list, NULL);
    PointerToDialog (dlg->field_list_to, NULL);
    SendMessageToDialog (dlg->field_list, NUM_VIB_MSG + 1);
    PointerToDialog (dlg->rna_type, NULL);
    PointerToDialog (dlg->text_portion, NULL);
    SafeSetStatus (dlg->leave_on_original, FALSE);
    SafeSetStatus (dlg->remove_parsed, FALSE);
    ClearTextAECRRNAQualDlg (dlg);
  }
}

static void AECRRNAQualToDialog (DialoG d, Pointer userdata)
{
  AECRRNAQualDlgPtr dlg;
  AECRRNAQualPtr    data;
  
  dlg = (AECRRNAQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  ResetAECRRNAQualDlg (dlg);
  data = (AECRRNAQualPtr) userdata;
  if (data != NULL)
  {
    PointerToDialog (dlg->field_list, data->field_list);
    PointerToDialog (dlg->field_list_to, data->field_list_to);
    PointerToDialog (dlg->rna_type, data->rna_type);
    PointerToDialog (dlg->edit_apply, data->edit_apply);
    PointerToDialog (dlg->codons, data->codons);
    if (data->field_list != NULL && data->field_list->data.intvalue == RNAFIELD_NCRNA_CLASS) {
      if (data->edit_apply == NULL) {
        PointerToDialog (dlg->ncrna_class, NULL);
      } else {
        PointerToDialog (dlg->ncrna_class, data->edit_apply->apply_txt);
      }
    }
    PointerToDialog (dlg->anticodon, data->anticodon);
    PointerToDialog (dlg->text_portion, data->text_portion);
    SafeSetStatus (dlg->leave_on_original, data->leave_on_original);
    SafeSetStatus (dlg->remove_parsed, data->remove_parsed);
  }
}

static Pointer DialogToAECRRNAQual (DialoG d)
{
  AECRRNAQualDlgPtr dlg;
  AECRRNAQualPtr    data;
  
  dlg = (AECRRNAQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  data = (AECRRNAQualPtr) MemNew (sizeof (AECRRNAQualData));
  if (data != NULL)
  {
    data->field_list = DialogToPointer (dlg->field_list);
    data->field_list_to = DialogToPointer (dlg->field_list_to);
    data->rna_type = DialogToPointer (dlg->rna_type);
    if (data->field_list != NULL && data->field_list->data.intvalue == RNAFIELD_NCRNA_CLASS
        && dlg->ncrna_class != NULL) {
      data->edit_apply = EditApplyNew ();
      data->edit_apply->apply_txt = DialogToPointer (dlg->ncrna_class);
    } else {
      data->edit_apply = DialogToPointer (dlg->edit_apply);
    }
    data->codons = DialogToPointer (dlg->codons);
    data->anticodon = DialogToPointer (dlg->anticodon);
    data->text_portion = DialogToPointer (dlg->text_portion);
    data->leave_on_original = (dlg->leave_on_original == NULL) ? FALSE : GetStatus (dlg->leave_on_original);
    data->remove_parsed = (dlg->remove_parsed == NULL) ? FALSE : GetStatus (dlg->remove_parsed);
  }
  return data;
}

static ValNodePtr GetRNAFeatureList (SeqEntryPtr sep, FilterSetPtr fsp, RnaTypePtr rtp);
static GetSamplePtr GetRnaExistingTextForList (ValNodePtr rna_list, ValNodePtr requested_field);

static CharPtr 
RNAQualAutopopulate 
(RnaTypePtr rna_type,
 ValNodePtr field_list,
 Uint2      entityID)
{
  SeqEntryPtr  sep;
  GetSamplePtr gsp;
  ValNodePtr   rna_list;
  CharPtr      rval = NULL;

  sep = GetTopSeqEntryForEntityID (entityID);
  rna_list = GetRNAFeatureList (sep, NULL, rna_type);
  gsp = GetRnaExistingTextForList (rna_list, field_list);
  if (gsp != NULL && !StringHasNoText (gsp->sample_text))
  {
    rval = gsp->sample_text;
    gsp->sample_text = NULL;
  }
  GetSampleFree (gsp);
  rna_list = ValNodeFree (rna_list);
  return rval;
}

static void RNAQualChangeNotify (Pointer userdata)
{
  AECRRNAQualDlgPtr  dlg;
  ValNodePtr         field_list = NULL;
  RnaTypePtr         rna_type;
  EditApplyPtr       eap;
  CharPtr            ncrna_class;

  dlg = (AECRRNAQualDlgPtr) userdata;
  if (dlg == NULL)
  {
    return;
  }
  
  field_list = DialogToPointer (dlg->field_list);
  if (field_list == NULL)
  {
    Hide (dlg->edit_apply);
    Hide (dlg->anticodon);
    Hide (dlg->codons);
    Hide (dlg->ncrna_class);
  }  
  else if (field_list->data.intvalue == RNAFIELD_CODONS_RECOGNIZED)
  {
    Show (dlg->codons);
    Hide (dlg->edit_apply);
    Hide (dlg->anticodon);
    Hide (dlg->ncrna_class);
  }
  else if (field_list->data.intvalue == RNAFIELD_ANTICODON)
  {
    Show (dlg->anticodon);
    Hide (dlg->edit_apply);
    Hide (dlg->codons);
    Hide (dlg->ncrna_class);
  }
  else if (field_list->data.intvalue == RNAFIELD_NCRNA_CLASS && dlg->ncrna_class != NULL)
  {
    Show (dlg->ncrna_class);
    Hide (dlg->codons);
    Hide (dlg->edit_apply);
    Hide (dlg->anticodon);
    /* autopopulate */
    rna_type = DialogToPointer (dlg->rna_type);
    ncrna_class = RNAQualAutopopulate (rna_type, field_list, dlg->entityID);
    PointerToDialog (dlg->ncrna_class, ncrna_class);
    ncrna_class = MemFree (ncrna_class);
  }
  else
  {
    Show (dlg->edit_apply);
    Hide (dlg->codons);
    Hide (dlg->anticodon);
    Hide (dlg->ncrna_class);

    /* autopopulate text fields */
    if (dlg->edit_apply != NULL)
    {    
      rna_type = DialogToPointer (dlg->rna_type);
      eap = DialogToPointer (dlg->edit_apply);
      eap->apply_txt = MemFree (eap->apply_txt);
      eap->apply_txt = RNAQualAutopopulate (rna_type, field_list, dlg->entityID);
      PointerToDialog (dlg->edit_apply, eap);
      EditApplyFree (eap);    
      rna_type = RnaTypeFree (rna_type);
    }
  }
  
  field_list = ValNodeFree (field_list);
  
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify)(dlg->change_userdata);
  }  
}

static void AECRRNAQualMessage (DialoG d, Int2 mssg)

{
  AECRRNAQualDlgPtr  dlg;

  dlg = (AECRRNAQualDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) 
    {
      case VIB_MSG_INIT :
        /* reset list */
        ResetAECRRNAQualDlg (dlg);
        break;
      case VIB_MSG_ENTER :
        Select (dlg->rna_type);
        break;
      case AECR_VIB_MSG_SET_DEFAULT :
        SendMessageToDialog (dlg->field_list, AECR_VIB_MSG_SET_DEFAULT);
        SendMessageToDialog (dlg->rna_type, AECR_VIB_MSG_SET_DEFAULT);
        ClearTextAECRRNAQualDlg (dlg);
        break;
      case AECR_VIB_MSG_CLEAR_TEXT :
        ClearTextAECRRNAQualDlg (dlg);
        break;
      case AECR_VIB_MSG_AUTOPOPULATE:
        RNAQualChangeNotify (dlg);
        break; 
      default :
        break;
    }
  }
}

static ValNodePtr TestAECRRNAQual (DialoG d)
{
  AECRRNAQualDlgPtr dlg;
  ValNodePtr         total_err_list = NULL;
  ValNodePtr         field_list;
  
  dlg = (AECRRNAQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }

  total_err_list = TestDialog (dlg->rna_type);
  
  ValNodeLink (&total_err_list, TestDialog (dlg->field_list));
      
  field_list = DialogToPointer (dlg->field_list);
  if (field_list != NULL) {
    if (field_list->data.intvalue == RNAFIELD_NCRNA_CLASS && dlg->ncrna_class != NULL) {
      ValNodeLink (&total_err_list, TestDialog (dlg->ncrna_class));
    } else if (field_list->data.intvalue != RNAFIELD_CODONS_RECOGNIZED) {
      ValNodeLink (&total_err_list, TestDialog (dlg->edit_apply));
    }
  }
      
  return total_err_list;
}


static void AECRRNAQualButtonChange (ButtoN b)
{
  AECRRNAQualDlgPtr dlg;

  dlg = (AECRRNAQualDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) return;

  if (dlg->change_notify != NULL) 
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static DialoG RNAQualListDialog
(GrouP                    h,
 Int4                     action_choice, 
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Uint2                    entityID)
{
  AECRRNAQualDlgPtr dlg;
  GrouP              p, g1, g2 = NULL;
  SeqEntryPtr        sep;
  
  dlg = (AECRRNAQualDlgPtr) MemNew (sizeof (AECRRNAQualDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = AECRRNAQualToDialog;
  dlg->fromdialog = DialogToAECRRNAQual;
  dlg->dialogmessage = AECRRNAQualMessage;
  dlg->testdialog = TestAECRRNAQual;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->entityID = entityID;


  if (action_choice == AECR_CONVERT || action_choice == AECR_SWAP || action_choice == AECR_PARSE)
  {
    g1 = HiddenGroup (p, 3, 0, NULL);
    StaticPrompt (g1, "", 0, dialogTextHeight, systemFont, 'l');
    StaticPrompt (g1, "From", 0, dialogTextHeight, systemFont, 'l');
    StaticPrompt (g1, "To", 0, dialogTextHeight, systemFont, 'l');
    dlg->rna_type = RnaTypeDialog (g1, TRUE, RNAQualChangeNotify, dlg);

    dlg->field_list = RNAFieldSelectionDialog (g1, FALSE, 
                                          change_notify, 
                                          change_userdata);
    dlg->field_list_to = RNAFieldSelectionDialog (g1, FALSE, 
                                             change_notify, 
                                             change_userdata);
    if (action_choice == AECR_CONVERT)
    {
      dlg->leave_on_original = CheckBox (p, "Leave on original", AECRRNAQualButtonChange);
      SetObjectExtra (dlg->leave_on_original, dlg, NULL);
    }
  }
  else
  {
    g1 = HiddenGroup (p, 2, 0, NULL);
    StaticPrompt (g1, "", 0, dialogTextHeight, systemFont, 'l');
    StaticPrompt (g1, "RNA Field", 0, dialogTextHeight, systemFont, 'l');
    dlg->rna_type = RnaTypeDialog (g1, TRUE, RNAQualChangeNotify, dlg);
    if (action_choice == AECR_REMOVE)
    {
      dlg->field_list = RNARemoveFieldSelectionDialog (g1, TRUE, 
                                            RNAQualChangeNotify, dlg);
    }
    else
    {
      dlg->field_list = RNAFieldSelectionDialog (g1, FALSE, 
                                            RNAQualChangeNotify, 
                                            dlg);
    }
  }

  if (action_choice == AECR_APPLY)
  {
    g2 = HiddenGroup (p, 0, 0, NULL);
    dlg->edit_apply = EditApplyDialog (g2, eEditApplyChoice_Apply, "New Text", NULL,
                                       change_notify, change_userdata);
  
    dlg->codons = CreateVisibleStringDialog (g2, 3, -1, 4);

    dlg->ncrna_class = CreatencRNAClassDialog (g2, FALSE, change_notify, change_userdata);
  
    sep = GetTopSeqEntryForEntityID (entityID);
  
    dlg->anticodon = NULL;
  /*
    dlg->anticodon = CreateIntervalEditorDialogEx (g2, NULL, 2, 2, sep, TRUE, FALSE,
                                                    TRUE, TRUE, FALSE,
                                                    &(dlg->ff),
                                                    NULL); */

    AlignObjects (ALIGN_CENTER, (HANDLE) dlg->edit_apply,
                                (HANDLE) dlg->codons,
                                (HANDLE) dlg->ncrna_class,
                                (HANDLE) dlg->anticodon,
                                NULL);
  }
  else if (action_choice == AECR_EDIT)
  {
    g2 = HiddenGroup (p, 1, 0, NULL);
    dlg->edit_apply = EditApplyDialog (g2, eEditApplyChoice_Edit,
                                       "New Text", NULL,
                                       change_notify, change_userdata);
    dlg->text_portion = NULL;
    dlg->remove_parsed = NULL;
  }
  else if (action_choice == AECR_PARSE)
  {
    g2 = HiddenGroup (p, -1, 0, NULL);
    dlg->edit_apply = NULL;
    dlg->text_portion = TextPortionXDialog(p);
    dlg->remove_parsed = CheckBox (p, "Remove from parsed field", SimpleAECRButtonChange);
    SetObjectExtra (dlg->remove_parsed, dlg, NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) dlg->edit_apply, (HANDLE) dlg->remove_parsed, NULL);
  }
  else
  {
    dlg->edit_apply = NULL;
    dlg->text_portion = NULL;
    dlg->remove_parsed = NULL;
  }
  

  AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) g2, NULL);
  
  RNAQualChangeNotify (dlg);
  
  return (DialoG) p;
}


static void ApplyRNACodonsCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  ValNodePtr codons;
  RnaRefPtr  rrp;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA)
  {
    return;
  }
  
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  codons = (ValNodePtr) userdata;
  if (rrp != NULL && rrp->type == 3 && rrp->ext.choice == 2) 
  {
    AddCodonListTotRNA (rrp->ext.value.ptrvalue, codons);
  }
}

static void ApplyRNAAnticodonCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  RnaRefPtr  rrp;
  tRNAPtr    trp;
  SeqLocPtr  anticodon_loc, slp;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA 
      || sfp->data.value.ptrvalue == NULL)
  {
    return;
  }
  
  anticodon_loc = (SeqLocPtr) userdata;

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;

  if (rrp->ext.choice == 2) 
  {
    trp = (tRNAPtr) rrp->ext.value.ptrvalue;
    if (trp != NULL) 
    {
      trp->anticodon = SeqLocFree (trp->anticodon);
      if (anticodon_loc != NULL)
      {
        slp = (SeqLocPtr) AsnIoMemCopy (anticodon_loc,
                                        (AsnReadFunc) SeqLocAsnRead,
                                        (AsnWriteFunc) SeqLocAsnWrite);

        slp = SeqLocReplaceID (slp, SeqLocId (sfp->location));
        trp->anticodon = slp;
      }
    }
  }
}


static void CollectRNAFeaturesCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{  
  if (sfp != NULL && sfp->data.choice == SEQFEAT_RNA && userdata != NULL)
  {
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_SEQFEAT, sfp);
  }
  
}


static ValNodePtr GetRNAFeatureList (SeqEntryPtr sep, FilterSetPtr fsp, RnaTypePtr rtp)
{
  ValNodePtr rna_list = NULL, vnp_rna;
  ValNodePtr not_match;
  SeqFeatPtr sfp;
  Boolean    matches;

  if (sep == NULL) return NULL;

  OperateOnSeqEntryConstrainedObjects (sep, fsp, CollectRNAFeaturesCallback,
                                       NULL, SEQFEAT_RNA, 0, 0, &rna_list);  

  if (rtp != NULL)
  {
    for (vnp_rna = rna_list; vnp_rna != NULL; vnp_rna = vnp_rna->next)
    {  
      matches = FALSE;
      sfp = vnp_rna->data.ptrvalue;
      if (!MatchesRnaType (sfp, rtp))
      {
        vnp_rna->choice = 0;
      }
    }

    not_match = ValNodeExtractList (&rna_list, 0);
    not_match = ValNodeFree (not_match);
  }
  return rna_list;
}


static GetSamplePtr GetRnaExistingTextForList (ValNodePtr rna_list, ValNodePtr requested_field)
{
  GetSamplePtr gsp;
  ValNodePtr   vnp;
  SeqFeatPtr   sfp;

  if (rna_list == NULL || requested_field == NULL) return NULL;

  gsp = GetSampleNew ();
  gsp->sample_text = NULL;
  gsp->fieldstring_func = GetRNAFieldString;
  gsp->descrstring_func = NULL;
  gsp->free_vn_proc = NULL;
  gsp->copy_vn_proc = IntValNodeCopy;
  gsp->requested_field = (gsp->copy_vn_proc) (requested_field);
  gsp->num_found = 0;
  gsp->all_same = TRUE;

  for (vnp = rna_list; vnp != NULL; vnp = vnp->next)
  {  
    sfp = vnp->data.ptrvalue;
    GetSampleFeatureCallback (sfp, gsp, NULL);
  }
  return gsp;
}


static Boolean ApplyToRNAList (ValNodePtr rna_list, AECRRNAQualPtr rp, FilterSetPtr fsp)
{
  ValNodePtr      vnp_rna;
  SeqFeatPtr      sfp;
  ApplyValueData  avd;
  Boolean         rval = TRUE;
  GetSamplePtr    gsp;

  if (rna_list == NULL || rp == NULL || rp->field_list == NULL)
  {
    return FALSE;
  }

  if (rp->field_list->data.intvalue == RNAFIELD_CODONS_RECOGNIZED)
  {
    for (vnp_rna = rna_list; vnp_rna != NULL; vnp_rna = vnp_rna->next)
    {
      sfp = vnp_rna->data.ptrvalue;
      if (sfp->idx.subtype == FEATDEF_tRNA)
      {
        ApplyRNACodonsCallback (sfp, rp->codons, fsp);
      }
    }
  }
  else if (rp->field_list->data.intvalue == RNAFIELD_ANTICODON)
  {
    for (vnp_rna = rna_list; vnp_rna != NULL; vnp_rna = vnp_rna->next)
    {
      sfp = vnp_rna->data.ptrvalue;
      if (sfp->idx.subtype == FEATDEF_tRNA)
      {
        ApplyRNAAnticodonCallback (sfp, rp->anticodon, fsp);
      }
    }
  }
  else
  {
    avd.field_list = rp->field_list;
    avd.where_to_replace = EditApplyFindLocation_anywhere;
    
    /* get handling for existing text */
    gsp = GetRnaExistingTextForList (rna_list, rp->field_list);
    avd.etp = GetExistingTextHandlerInfo (gsp == NULL ? 0 : gsp->num_found, FALSE);
    gsp = GetSampleFree (gsp);
  
    if (avd.etp != NULL 
        && avd.etp->existing_text_choice == eExistingTextChoiceCancel)
    {
      rval = FALSE;
    }

    if (rval)
    {
      AddEditApplyDataToApplyValue (AECR_APPLY, rp->edit_apply, &avd);
     
      for (vnp_rna = rna_list; vnp_rna != NULL; vnp_rna = vnp_rna->next)
      {
        sfp = vnp_rna->data.ptrvalue;
        SetRNAFieldString (sfp, &avd, fsp);
      }
      avd.text_to_replace = MemFree (avd.text_to_replace);
      avd.new_text = MemFree (avd.new_text);
    }
    avd.etp = MemFree (avd.etp);
  }
  return rval;
}


static Boolean 
RNAQualConvertSwapParseList
(ValNodePtr      rna_list,
 AECRRNAQualPtr  rp,
 Int4            action_choice,
 FilterSetPtr    fsp)
{
  ConvertFieldData       cfd;
  ConversionConflictData ccd;
  ValNodePtr             vnp_rna;
  ApplyValueData         avd;
  Boolean                rval = TRUE;
  SeqFeatPtr             sfp;
  
  if (rna_list == NULL
      || rp == NULL
      || action_choice < AECR_EDIT 
      || action_choice > NUM_AECR
      || rp->field_list == NULL
      || rp->field_list_to == NULL)
  {
    return FALSE;
  }

  avd.where_to_replace = EditApplyFindLocation_anywhere;
  
  cfd.src_field_list = rp->field_list;
  cfd.dst_field_list = rp->field_list_to;
  cfd.etp = NULL;
  cfd.get_str_func = GetRNAFieldString;
  cfd.set_str_func = SetRNAFieldString;
  cfd.remove_str_func = RemoveRNAField;
  cfd.get_d_str_func = NULL;
  cfd.set_d_str_func = NULL;
  cfd.remove_d_str_func = NULL;
  cfd.fsp = fsp;
  cfd.strip_name_from_text = FALSE;
  cfd.name_field_func = NULL;
  cfd.text_portion = rp->text_portion;
  cfd.remove_parsed = rp->remove_parsed;

  if (action_choice == AECR_CONVERT)
  {
    if (rp->leave_on_original)
    {
      cfd.convert_type = CONVERT_TYPE_COPY;
    }
    else
    {
      cfd.convert_type = CONVERT_TYPE_MOVE;
    }
  }
  else if (action_choice == AECR_SWAP)
  {
    cfd.convert_type = CONVERT_TYPE_SWAP;
  }
  else if (action_choice == AECR_PARSE)
  {
    cfd.convert_type = CONVERT_TYPE_PARSE;
  }
  
  if (cfd.convert_type == CONVERT_TYPE_MOVE || cfd.convert_type == CONVERT_TYPE_COPY || cfd.convert_type == CONVERT_TYPE_PARSE)
  {
    /* get existing text sample */
    ccd.cfp = &cfd;
    ccd.gsp = GetSampleNew();

    for (vnp_rna = rna_list; vnp_rna != NULL; vnp_rna = vnp_rna->next)
    {
      sfp = vnp_rna->data.ptrvalue;
      GetConversionConflictsFeatureCallback (sfp, &ccd, fsp);
    }
    cfd.etp = GetExistingTextHandlerInfo (ccd.gsp == NULL ? 0 : ccd.gsp->num_found, FALSE);
    ccd.gsp = GetSampleFree (ccd.gsp);
    if (cfd.etp != NULL && cfd.etp->existing_text_choice == eExistingTextChoiceCancel)
    {
      rval = FALSE;
    }
  }
  else
  {
    cfd.etp = NULL;
  }  

  if (rval)
  {
    for (vnp_rna = rna_list; vnp_rna != NULL; vnp_rna = vnp_rna->next)
    {
      sfp = vnp_rna->data.ptrvalue;
      ConvertFeatureFieldCallback (sfp, &cfd, fsp);
    }
  }

  return rval;
}


static Boolean RNAQualEditRemoveList 
(ValNodePtr    rna_list,
 AECRRNAQualPtr  rp,
 Int4          action_choice,
 FilterSetPtr  fsp)
{
  SeqFeatPtr       sfp;
  ValNodePtr       vnp_rna;
  ApplyValueData   avd;
  Boolean          rval = TRUE;
  
  if (rna_list == NULL
      || rp == NULL
      || (action_choice != AECR_EDIT && action_choice != AECR_REMOVE)
      || rp->field_list == NULL)
  {
    return FALSE;
  }

  avd.where_to_replace = EditApplyFindLocation_anywhere;
  avd.field_list = rp->field_list;
  avd.etp = NULL;
  AddEditApplyDataToApplyValue (action_choice, rp->edit_apply, &avd);

  for (vnp_rna = rna_list; vnp_rna != NULL; vnp_rna = vnp_rna->next)
  {
    sfp = vnp_rna->data.ptrvalue;
    if (action_choice == AECR_EDIT)
    {
      SetRNAFieldString (sfp, &avd, fsp);
    }
    else if (action_choice == AECR_REMOVE)
    {
      RemoveRNAField (sfp, &avd, fsp);
    }
  }
    
  avd.text_to_replace = MemFree (avd.text_to_replace);
  avd.new_text = MemFree (avd.new_text);
  return TRUE;
}




static Boolean AcceptRNAQualApplyEditRemoveConvert (Pointer userdata)
{
  ApplyEditConvertRemovePtr ap;
  AECRRNAQualPtr            rp;
  Int4                      action_choice;
  Boolean                   rval = TRUE;
  FilterSetPtr              fsp;
  ValNodePtr                rna_list;
  SeqEntryPtr               sep;

  ap = (ApplyEditConvertRemovePtr) userdata;
  if (ap == NULL)
  {
    return FALSE;
  }
  
  WatchCursor ();
  Update ();
  
  if (ap->crippled)
  {
    action_choice = ap->crippled_action;
  }
  else
  {
    action_choice = GetValue (ap->action_popup);
  }

  rp = (AECRRNAQualPtr) DialogToPointer (ap->aecr_pages[action_choice - 1]);
  if (rp == NULL) return FALSE;

  fsp = (FilterSetPtr) DialogToPointer (ap->constraints);
  sep = GetTopSeqEntryForEntityID (ap->input_entityID);
  rna_list = GetRNAFeatureList (sep, fsp, rp->rna_type);
  
  if (action_choice == AECR_APPLY)
  {
    rval = ApplyToRNAList (rna_list, rp, fsp);
  }
  else if (action_choice == AECR_CONVERT || action_choice == AECR_SWAP
           || action_choice == AECR_PARSE)
  {
    rval = RNAQualConvertSwapParseList (rna_list, rp, action_choice, fsp);
  }
  else if (action_choice > AECR_APPLY && action_choice <= NUM_AECR)
  {
    rval = RNAQualEditRemoveList (rna_list, rp, action_choice, fsp);
  }
    
  if (rval)
  {
    ObjMgrSetDirtyFlag (ap->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, ap->input_entityID, 0, 0);
  }
  ArrowCursor ();
  Update (); 
  
  return rval;
}

static void 
RNAQualApplyEditRemoveConvert 
(IteM    i, 
 Int4    first_action_choice,
 Boolean crippled)
{
  ApplyEditConvertRemoveCombo (i, first_action_choice, crippled,
                               "RNA Qualifiers",
                               RNAQualListDialog, TRUE, GetRNAQualSample,
                               TRUE, FALSE, TRUE, FALSE, FALSE, "Where feature text",
                               AcceptRNAQualApplyEditRemoveConvert, NULL, NULL,
                               CheckFeaturesForPresample);
}

extern void ApplyRNAQual (IteM i)
{
  RNAQualApplyEditRemoveConvert (i, AECR_APPLY, FALSE);
}

extern void PublicApplyRNAQual (IteM i)
{
  RNAQualApplyEditRemoveConvert (i, AECR_APPLY, TRUE);
}

extern void EditRNAQual (IteM i)
{
  RNAQualApplyEditRemoveConvert (i, AECR_EDIT, FALSE);
}

extern void PublicEditRNAQual (IteM i)
{
  RNAQualApplyEditRemoveConvert (i, AECR_EDIT, TRUE);
}

extern void ConvertRNAQual (IteM i)
{
  RNAQualApplyEditRemoveConvert (i, AECR_CONVERT, FALSE);
}

extern void SwapRNAQual (IteM i)
{
  RNAQualApplyEditRemoveConvert (i, AECR_SWAP, FALSE);
}

extern void RemoveRNAQual (IteM i)
{
  RNAQualApplyEditRemoveConvert (i, AECR_REMOVE, FALSE);
}


typedef struct removegbqualdlg 
{
  DIALOG_MESSAGE_BLOCK
  DialoG          field_list;
  DialoG          subtype_list;
  GrouP           list_or_txt_grp;
  DialoG          qual_name_constraint;
  
  /* These are necessary to produce the SimpleAECR structure for output */
  FreeValNodeProc free_field_vn_proc;
  FreeValNodeProc free_subtype_vn_proc;
  
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
  Uint2                    entityID;
  
} RemoveGBQualDlgData, PNTR RemoveGBQualDlgPtr;


typedef struct removegbqual
{
  ValNodePtr          field_list;
  ValNodePtr          subtype_list;
  StringConstraintXPtr scp;
  Boolean             use_field_list;
} RemoveGBQualData, PNTR RemoveGBQualPtr;

static RemoveGBQualPtr RemoveGBQualFree (RemoveGBQualPtr ap)
{
  if (ap != NULL)
  {
    ap->field_list = ValNodeFree (ap->field_list);
    ap->subtype_list = ValNodeFree (ap->subtype_list);
    ap->scp = StringConstraintXFree (ap->scp);
    ap = MemFree (ap);
  }
  return ap;
}


static void ClearTextRemoveGBQualDlg (RemoveGBQualDlgPtr dlg)
{
  if (dlg != NULL)
  {
    PointerToDialog (dlg->qual_name_constraint, NULL);    
  }
}

static void ResetRemoveGBQualDlg (RemoveGBQualDlgPtr dlg)
{
  if (dlg != NULL)
  {
    PointerToDialog (dlg->field_list, NULL);
    SendMessageToDialog (dlg->field_list, NUM_VIB_MSG + 1);
    PointerToDialog (dlg->subtype_list, NULL);
    ClearTextRemoveGBQualDlg (dlg);
    SetValue (dlg->list_or_txt_grp, 1);
    Enable (dlg->subtype_list);
    Disable (dlg->qual_name_constraint);
  }
}


static void RemoveGBQualToDialog (DialoG d, Pointer userdata)
{
  RemoveGBQualDlgPtr dlg;
  RemoveGBQualPtr    data;
  
  dlg = (RemoveGBQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
    
  data = (RemoveGBQualPtr) userdata;
  if (data == NULL)
  {
    PointerToDialog (dlg->field_list, NULL);
    PointerToDialog (dlg->subtype_list, NULL);
    ClearTextRemoveGBQualDlg (dlg);
    SetValue (dlg->list_or_txt_grp, 1);
    Enable (dlg->subtype_list);
    Disable (dlg->qual_name_constraint);
  }
  else
  {
    PointerToDialog (dlg->field_list, data->field_list);
    PointerToDialog (dlg->subtype_list, data->subtype_list);
    PointerToDialog (dlg->qual_name_constraint, data->scp);
    if (data->use_field_list)
    {
      SetValue (dlg->list_or_txt_grp, 1);
      Enable (dlg->subtype_list);
      Disable (dlg->qual_name_constraint);
    }
    else
    {
      SetValue (dlg->list_or_txt_grp, 2);
      Disable (dlg->subtype_list);
      Enable (dlg->qual_name_constraint);
    }
  }
}


static Pointer DialogToRemoveGBQual (DialoG d)
{
  RemoveGBQualDlgPtr dlg;
  RemoveGBQualPtr    data;
  
  dlg = (RemoveGBQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  data = (RemoveGBQualPtr) MemNew (sizeof (RemoveGBQualData));
  if (data != NULL)
  {
    data->field_list = DialogToPointer (dlg->field_list);
    data->subtype_list = DialogToPointer (dlg->subtype_list);
    data->scp = DialogToPointer (dlg->qual_name_constraint);
    if (GetValue (dlg->list_or_txt_grp) == 1)
    {
      data->use_field_list = TRUE;
    }
    else
    {
      data->use_field_list = FALSE;
    }
  }
  return data;
}


static void RemoveGBQualMessage (DialoG d, Int2 mssg)

{
  RemoveGBQualDlgPtr  dlg;

  dlg = (RemoveGBQualDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) 
    {
      case VIB_MSG_INIT :
        /* reset list */
        ResetRemoveGBQualDlg (dlg);
        break;
      case VIB_MSG_ENTER :
        if (dlg->subtype_list != NULL)
        {
          Select (dlg->subtype_list);
        }
        else
        {
          Select (dlg->field_list);
        }
        break;
      case AECR_VIB_MSG_SET_DEFAULT :
        SendMessageToDialog (dlg->field_list, AECR_VIB_MSG_SET_DEFAULT);
        SendMessageToDialog (dlg->subtype_list, AECR_VIB_MSG_SET_DEFAULT);
        ClearTextRemoveGBQualDlg (dlg);
        break;
      case AECR_VIB_MSG_CLEAR_TEXT :
        ClearTextRemoveGBQualDlg (dlg);
        break;
      case AECR_VIB_MSG_AUTOPOPULATE:
        break; 
      default :
        break;
    }
  }
}

static ValNodePtr TestRemoveGBQual (DialoG d)
{
  RemoveGBQualDlgPtr dlg;
  ValNodePtr         err_list = NULL;
  ValNodePtr         total_err_list = NULL;
  
  dlg = (RemoveGBQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }

  total_err_list = TestDialog (dlg->subtype_list);

  if (GetValue (dlg->list_or_txt_grp) == 1)
  {
    err_list = TestDialog (dlg->field_list);
    total_err_list = ValNodeAppend (total_err_list, err_list);
  }      
      
  return total_err_list;
}


static void GBQualChangeNotify (GrouP b)
{
  RemoveGBQualDlgPtr dlg;

  dlg = (RemoveGBQualDlgPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }
  
  if (GetValue (dlg->list_or_txt_grp) == 1)
  {
    Enable (dlg->field_list);
    Disable (dlg->qual_name_constraint);
  }
  else
  {
    Disable (dlg->field_list);
    Enable (dlg->qual_name_constraint);
  }
    
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify)(dlg->change_userdata);
  }  
}


static DialoG GBQualRemoveDialog
(GrouP                    h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Uint2                    entityID)
{
  RemoveGBQualDlgPtr dlg;
  GrouP              p, g1;
  SeqEntryPtr        sep;
  
  dlg = (RemoveGBQualDlgPtr) MemNew (sizeof (RemoveGBQualDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = RemoveGBQualToDialog;
  dlg->fromdialog = DialogToRemoveGBQual;
  dlg->dialogmessage = RemoveGBQualMessage;
  dlg->testdialog = TestRemoveGBQual;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->entityID = entityID;
  sep = GetTopSeqEntryForEntityID(entityID);
  
  g1 = HiddenGroup (p, 2, 0, NULL);
  StaticPrompt (g1, "Feature", 0, dialogTextHeight, systemFont, 'l');
  StaticPrompt (g1, "Qualifier", 0, dialogTextHeight, systemFont, 'l');
  dlg->subtype_list = FeatureSelectionDialogEx (g1, TRUE, sep,
                                                 change_notify,
                                                 change_userdata);
  dlg->list_or_txt_grp = HiddenGroup (g1, 2, 0, GBQualChangeNotify);
  SetObjectExtra (dlg->list_or_txt_grp, dlg, NULL);
  RadioButton (dlg->list_or_txt_grp, "Choose from Standard List");
  RadioButton (dlg->list_or_txt_grp, "Select Qual Name Constraint");
  dlg->field_list = GBQualSelectionDialog (dlg->list_or_txt_grp, FALSE, 
                                                change_notify,
                                                change_userdata);
  dlg->qual_name_constraint = StringConstraintDialogX (dlg->list_or_txt_grp, "Where qualifier name contains", FALSE);
  SetValue (dlg->list_or_txt_grp, 1);
  Disable (dlg->qual_name_constraint);
    
  return (DialoG) p;
}


static DialoG 
GBQualListDialog
(GrouP                    h,
 Int4                     action_choice,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Uint2                    entityID)
{
  if (action_choice == AECR_REMOVE)
  {
    return GBQualRemoveDialog (h, change_notify, change_userdata, entityID);
  }
  else
  {
    return SimpleAECRDialog (h, action_choice, change_notify, change_userdata,
                             NULL, ValNodeSimpleDataFree,
                             GBQualSelectionDialog, FeatureSelectionDialogEx, NULL,
                             "Qualifier", "Feature", TRUE, entityID);
  }
}

static GetSamplePtr 
GetGBQualExistingText 
(SeqEntryPtr  sep,
 ValNodePtr   feature_type_list,
 ValNodePtr   requested_field,
 FilterSetPtr fsp)
{
  ValNodePtr   vnp;
  GetSamplePtr gsp;
  
  if (sep == NULL || feature_type_list == NULL)
  {
    return NULL;
  }
  
  gsp = (GetSamplePtr) MemNew (sizeof (GetSampleData));
  if (gsp == NULL)
  {
    return NULL;
  }
  gsp->sample_text = NULL;
  gsp->fieldstring_func = GetGBQualString;
  gsp->descrstring_func = NULL;
  gsp->free_vn_proc = NULL;
  gsp->copy_vn_proc = IntValNodeCopy;
    
  gsp->requested_field = (gsp->copy_vn_proc) (requested_field);
  gsp->num_found = 0;
  gsp->all_same = TRUE;
  for (vnp = feature_type_list; vnp != NULL; vnp = vnp->next)
  {
    OperateOnSeqEntryConstrainedObjects (sep, fsp, GetSampleFeatureCallback, 
                                         NULL,
                                         0, vnp->choice, 0, gsp);
  }
  return gsp;  
}


static Boolean AcceptGBQualApplyEditRemoveConvert (Pointer userdata)
{
  ApplyEditConvertRemovePtr ap;
  SimpleAECRPtr             gp = NULL;
  RemoveGBQualPtr           rgp = NULL;
  Int4                      action_choice;
  SeqEntryPtr               sep;
  ApplyValueData            avd;
  ConvertFieldData          cfd;
  FilterSetPtr              fsp;
  GetSamplePtr              gsp = NULL, gsp_new, gsp_sum;
  ValNodePtr                vnp, qual_vnp;
  Boolean                   rval = TRUE;

  ValNodePtr                field_list = NULL;
  ValNodePtr                subtype_list = NULL;
  EditApplyPtr              eap = NULL;


  ap = (ApplyEditConvertRemovePtr) userdata;
  if (ap == NULL)
  {
    return FALSE;
  }
  
  if (ap->crippled)
  {
    action_choice = ap->crippled_action;
  }
  else
  {
    action_choice = GetValue (ap->action_popup);
  }
  if (action_choice < 1 || action_choice > NUM_AECR)
  {
    return FALSE;
  }

  avd.where_to_replace = EditApplyFindLocation_anywhere;
  
  if (action_choice == AECR_REMOVE)
  {
    rgp = (RemoveGBQualPtr) DialogToPointer (ap->aecr_pages[action_choice - 1]);
    if (rgp == NULL)
    {
      return FALSE;
    }
    field_list = rgp->field_list;
    subtype_list = rgp->subtype_list;
  }
  else
  {
    gp = (SimpleAECRPtr) DialogToPointer (ap->aecr_pages[action_choice - 1]);
    if (gp == NULL)
    {
      return FALSE;
    }
    field_list = gp->field_list;
    subtype_list = gp->subtype_list;
    eap = gp->edit_apply;
  }

  WatchCursor ();
  Update ();
  
  sep = GetTopSeqEntryForEntityID (ap->input_entityID);
  if (subtype_list == NULL || (action_choice != AECR_REMOVE && field_list == NULL)
      || (action_choice == AECR_REMOVE && rgp->use_field_list && field_list == NULL)
      || ((action_choice == AECR_CONVERT || action_choice == AECR_SWAP || action_choice == AECR_PARSE) 
           && gp->field_list_to == NULL))
  {
    gp = SimpleAECRFree (gp);
    rgp = RemoveGBQualFree (rgp);
    ArrowCursor ();
    Update ();
    return FALSE;
  }

  fsp = (FilterSetPtr) DialogToPointer (ap->constraints);
  
  if (action_choice == AECR_CONVERT || action_choice == AECR_SWAP || action_choice == AECR_PARSE)
  {
    cfd.src_field_list = gp->field_list;
    cfd.dst_field_list = gp->field_list_to;
    cfd.etp = NULL;
    cfd.get_str_func = GetGBQualString;
    cfd.set_str_func = SetGBQualString;
    cfd.remove_str_func = RemoveGBQualField;
    cfd.get_d_str_func = NULL;
    cfd.set_d_str_func = NULL;
    cfd.remove_d_str_func = NULL;
    cfd.fsp = fsp;
    cfd.strip_name_from_text = gp->strip_name_from_text;
    cfd.remove_parsed = gp->remove_parsed;
    cfd.name_field_func = NULL;
    cfd.text_portion = gp->text_portion;
    
    if (action_choice == AECR_CONVERT)
    {
      if (gp->leave_on_original)
      {
        cfd.convert_type = CONVERT_TYPE_COPY;
      }
      else
      {
        cfd.convert_type = CONVERT_TYPE_MOVE;
      }
    }
    else if (action_choice == AECR_SWAP)
    {
      cfd.convert_type = CONVERT_TYPE_SWAP;
    }
    else if (action_choice == AECR_PARSE)
    {
      cfd.convert_type = CONVERT_TYPE_PARSE;
    }
    
    if (cfd.convert_type == CONVERT_TYPE_MOVE || cfd.convert_type == CONVERT_TYPE_COPY || cfd.convert_type == CONVERT_TYPE_PARSE)
    {
      /* get existing text data */
      gsp_sum = NULL;
      for (vnp = gp->subtype_list; vnp != NULL; vnp = vnp->next)
      {
        gsp = CheckForConversionExistingTextInSeqEntry (sep, &cfd, fsp,
                                                        0,
                                                        vnp->choice,
                                                        0);
        gsp_new = GetSampleAdd (gsp_sum, gsp);
        gsp = GetSampleFree (gsp);
        gsp_sum = GetSampleFree (gsp_sum);
        gsp_sum = gsp_new;
      }
      cfd.etp = GetExistingTextHandlerInfo (gsp_sum == NULL ? 0 : gsp_sum->num_found, FALSE);
      gsp_sum = GetSampleFree (gsp_sum);
      if (cfd.etp != NULL && cfd.etp->existing_text_choice == eExistingTextChoiceCancel)
      {
        rval = FALSE;
      }
    }
    else
    {
      cfd.etp = NULL;
    }        
    
    if (rval)
    {
      for (vnp = gp->subtype_list; vnp != NULL; vnp = vnp->next)
      {
        OperateOnSeqEntryConstrainedObjects (sep, fsp, ConvertFeatureFieldCallback,
                                             NULL, 0, vnp->choice, 0, &cfd);  
      }
    }
    cfd.etp = MemFree (cfd.etp);
  }
  else if (action_choice == AECR_REMOVE && rgp != NULL && ! rgp->use_field_list)
  {
    for (vnp = subtype_list; vnp != NULL; vnp = vnp->next)
    {
      OperateOnSeqEntryConstrainedObjects (sep, fsp, RemoveGBQualByNameConstraint,
                                           NULL, 0, vnp->choice, 0, rgp->scp);  
    }
  }
  else
  {
    avd.field_list = field_list;
  
    /* get handling for existing text */
    if (action_choice == AECR_APPLY)
    {
      gsp = GetGBQualExistingText (sep, subtype_list,
                                   avd.field_list,
                                   fsp);
      avd.etp = GetExistingTextHandlerInfo (gsp == NULL ? 0 : gsp->num_found, FALSE);
      gsp = GetSampleFree (gsp);
    
      if (avd.etp != NULL 
        && avd.etp->existing_text_choice == eExistingTextChoiceCancel)
      {
        rval = FALSE;
      }
    }
    else
    {
      avd.etp = NULL;
    }
  
    AddEditApplyDataToApplyValue (action_choice, eap, &avd);
        
    for (vnp = subtype_list; vnp != NULL; vnp = vnp->next)
    {
      if (action_choice == AECR_EDIT || action_choice == AECR_APPLY)
      {
        OperateOnSeqEntryConstrainedObjects (sep, fsp, SetGBQualString,
                                             NULL, 0, vnp->choice, 0, &avd);  
      }
      else if (action_choice == AECR_REMOVE)
      {
        for (qual_vnp = field_list; qual_vnp != NULL; qual_vnp = qual_vnp->next)
        {
          avd.field_list = qual_vnp;
          OperateOnSeqEntryConstrainedObjects (sep, fsp, RemoveGBQualField,
                                               NULL, 0, vnp->choice, 0, &avd);  
        }
      }
    }

    avd.text_to_replace = MemFree (avd.text_to_replace);
    avd.new_text = MemFree (avd.new_text);
  }
  FilterSetFree (fsp);
  gp = SimpleAECRFree (gp);
  rgp = RemoveGBQualFree (rgp);
  
  if (rval)
  {
    ObjMgrSetDirtyFlag (ap->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, ap->input_entityID, 0, 0);
  }
  ArrowCursor ();
  Update (); 
  
  return rval;
}


static CharPtr GetGBQualSampleName (ValNodePtr vnp)
{
  CharPtr    label = NULL;
  Int2       feature_subtype, gb_field_choice;
  ValNodePtr feature_list, f_vnp;
  Int4       num_feature_choices;
  CharPtr    feature_label = "";
  
  if (vnp == NULL)
  {
    return NULL;
  }
  
  feature_list = BuildFeatureDialogList (FALSE, NULL);
  num_feature_choices = ValNodeLen (feature_list);
  
  feature_subtype = (vnp->data.intvalue - 1) / NumEditQualifiers;
  gb_field_choice = (vnp->data.intvalue - 1) % NumEditQualifiers;
  
  for (f_vnp = feature_list; f_vnp != NULL; f_vnp = f_vnp->next)
  {
    if (f_vnp->choice == feature_subtype)
    {
      feature_label = f_vnp->data.ptrvalue;
    }
  }
  
  
  label = (CharPtr) MemNew ((StringLen (feature_label)
                            + StringLen (EditQualifierList[gb_field_choice].name)
                            + 2) * sizeof (Char));
  if (label != NULL)
  {
    sprintf (label, "%s %s", feature_label,
             EditQualifierList[gb_field_choice].name);
  }
  ValNodeFreeData (feature_list);
  return label;
}

static CharPtr GetGBQualSampleFieldString (SeqFeatPtr sfp, ValNodePtr vnp, FilterSetPtr fsp)
{
  CharPtr   str = NULL;
  Int2      feature_subtype, gb_field_choice;
  ValNode   vn;
  
  if (sfp == NULL || vnp == NULL)
  {
    return NULL;
  }

  feature_subtype = (vnp->data.intvalue - 1) / NumEditQualifiers;
  gb_field_choice = (vnp->data.intvalue - 1) % NumEditQualifiers;
    
  if (feature_subtype != sfp->idx.subtype)
  {
    return NULL;
  }

  vn.next = NULL;
  vn.choice = 0;
  vn.data.intvalue = gb_field_choice + 1;
  str = GetGBQualString (sfp, &vn, NULL);
  
  return str;
}

static DialoG GetGBQualSample (GrouP h, Uint2 entityID)
{
  DialoG        d;
  SetSampleData ssd;
  Int4          qual_choice, list_index;  
  ValNodePtr    feature_list, f_vnp;    

  d = SampleDialog (h);
  ssd.fieldstring_func = GetGBQualSampleFieldString;
  ssd.descrstring_func = NULL;
  ssd.entityID = entityID;
  ssd.free_vn_proc = NULL;
  ssd.copy_vn_proc = IntValNodeCopy;
  ssd.match_vn_proc = IntValNodeMatch;
  ssd.label_vn_proc = GetGBQualSampleName;
  ssd.fsp = NULL;
 
  /* construct list of Features and GBQualifiers */
  ssd.field_list = NULL;

  feature_list = BuildFeatureDialogList (FALSE, NULL);
  
  for (f_vnp = feature_list; f_vnp != NULL; f_vnp = f_vnp->next)
  {
    for (qual_choice = 0; qual_choice < NumEditQualifiers; qual_choice++)
    {
      list_index = f_vnp->choice * NumEditQualifiers
                   + qual_choice;
      ValNodeAddInt (&ssd.field_list, 0, list_index + 1);    
    }
  }

  PointerToDialog (d, &ssd);
  
  /* now free up field list (SampleDialog maintains its own copy */
  ssd.field_list = ValNodeFree (ssd.field_list);
  feature_list = ValNodeFreeData (feature_list);
  return d;
}

static void 
GBQualApplyEditRemoveConvert 
(IteM    i,
 Int4    first_action_choice,
 Boolean crippled)
{
  ApplyEditConvertRemoveCombo (i, first_action_choice, crippled,
                               "Import Qualifiers",
                               GBQualListDialog, FALSE, GetGBQualSample,
                               TRUE, FALSE, TRUE, FALSE, FALSE, "Where feature text",
                               AcceptGBQualApplyEditRemoveConvert, NULL, NULL,
                               CheckFeaturesForPresample);
}

extern void ApplyGBQual (IteM i)
{
  GBQualApplyEditRemoveConvert (i, AECR_APPLY, FALSE);
}

extern void PublicApplyGBQual (IteM i)
{
  GBQualApplyEditRemoveConvert (i, AECR_APPLY, TRUE);
}

extern void EditGBQual (IteM i)
{
  GBQualApplyEditRemoveConvert (i, AECR_EDIT, FALSE);
}

extern void PublicEditGBQual (IteM i)
{
  GBQualApplyEditRemoveConvert (i, AECR_EDIT, TRUE);
}

extern void ConvertGBQual (IteM i)
{
  GBQualApplyEditRemoveConvert (i, AECR_CONVERT, FALSE);
}

extern void SwapGBQual (IteM i)
{
  GBQualApplyEditRemoveConvert (i, AECR_SWAP, FALSE);
}

extern void RemoveGBQual (IteM i)
{
  GBQualApplyEditRemoveConvert (i, AECR_REMOVE, FALSE);
}


typedef enum 
{
  eAllTargetsAllActions_BioSource = 0,
  eAllTargetsAllActions_GBQual,
  eAllTargetsAllActions_RNAQual,
  eAllTargetsAllActions_CDSGeneProtmRNA,
  eNumAllTargetsAllActions
} EAllTargetsAllActions;

typedef struct allactionsalltargetsdlg
{
  DIALOG_MESSAGE_BLOCK
  DialoG                   tbs; /* folder tab for switching between target types */
  Int2                     currentPage;
  GrouP                    target_pages [eNumAllTargetsAllActions];
  DialoG                   aecr_pages [NUM_AECR * eNumAllTargetsAllActions];
  PopuP                    action_popup;

  DialoG                   accept_cancel;
  DialoG                   constraints;
  ButtoN                   clear_constraints_on_action_change;
  Nlm_ChangeAECRActionProc change_action;
  Int4                     prev_page;
  Boolean                  crippled;
  Int4                     crippled_action;
  Uint2                    entityID;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} AllActionsAllTargetsDlgData, PNTR AllActionsAllTargetsDlgPtr;


extern void AllActionsAllTargets (IteM i)
{

}

typedef struct locustagtooldlocustag 
{
  FORM_MESSAGE_BLOCK
  DialoG         constraints;
  DialoG         accept_cancel;
  
  SeqEntryPtr    sep;
} LocusTagToOldLocusTagData, PNTR LocusTagToOldLocusTagPtr;

static void LocusTagToOldLocusTagActionClear (Pointer userdata)
{
  LocusTagToOldLocusTagPtr mp;
  
  mp = (LocusTagToOldLocusTagPtr) userdata;
  if (mp == NULL)
  {
    return;
  }
  PointerToDialog (mp->constraints, NULL);
}

static void LocusTagToOldLocusTagActionClearText (Pointer userdata)
{
  LocusTagToOldLocusTagPtr mp;
  FilterSetPtr             fsp;
  
  mp = (LocusTagToOldLocusTagPtr) userdata;
  if (mp == NULL)
  {
    return;
  }
  fsp = (FilterSetPtr) DialogToPointer (mp->constraints);
  FilterSetClearText (fsp);
  PointerToDialog (mp->constraints, fsp);
  FilterSetFree (fsp);
}


static void 
ConvertGeneLocusTagToOldLocusTagCallback 
(SeqFeatPtr   gene,
 Pointer      userdata, 
 FilterSetPtr fsp)
{
  GeneRefPtr               grp;
  CharPtr                  locus_tag;
  GBQualPtr                new_qual, prev_qual;
  LocusTagToOldLocusTagPtr mp;
  
  if (gene == NULL || gene->data.choice != SEQFEAT_GENE || userdata == NULL) return;
  
  mp = (LocusTagToOldLocusTagPtr)userdata;
  
  grp = (GeneRefPtr) gene->data.value.ptrvalue;
  if (grp == NULL || StringHasNoText (grp->locus_tag))
  {
      return;
  }
  locus_tag = grp->locus_tag;
  grp->locus_tag = NULL;
  
  new_qual = GBQualNew ();
  if (new_qual == NULL) return;
  new_qual->qual = StringSave ("old_locus_tag");
  new_qual->val = locus_tag;
  new_qual->next = NULL;
  
  if (gene->qual == NULL)
  {
      gene->qual = new_qual;
  }
  else
  {
      for (prev_qual = gene->qual; prev_qual->next != NULL; prev_qual = prev_qual->next)
      {
      }
      prev_qual->next = new_qual;
  }
}

static Boolean LocusTagToOldLocusTagAction (Pointer userdata)
{
  LocusTagToOldLocusTagPtr mp;
  FilterSetPtr             fsp;
  
  mp = (LocusTagToOldLocusTagPtr) userdata;
  if (mp == NULL)
  {
    return FALSE;
  }
  
  WatchCursor ();
  Update ();
  
  mp->sep = GetTopSeqEntryForEntityID (mp->input_entityID);
  fsp = (FilterSetPtr) DialogToPointer (mp->constraints);
  OperateOnSeqEntryConstrainedObjects (mp->sep, fsp, ConvertGeneLocusTagToOldLocusTagCallback, 
                                       NULL, SEQFEAT_GENE, 0, 0, mp);
  ObjMgrSetDirtyFlag (mp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, mp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();  
  return TRUE;
}

extern void ConvertLocusTagToOldLocusTag (IteM i)
{
  BaseFormPtr         bfp;
  LocusTagToOldLocusTagPtr     mp;
  WindoW              w;
  GrouP               h;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  mp = (LocusTagToOldLocusTagPtr) MemNew (sizeof (LocusTagToOldLocusTagData));
  if (mp == NULL) return;
  
  w = FixedWindow (-50, -33, -10, -10, "Convert Locus Tag to Old Locus Tag", StdCloseWindowProc); 
  SetObjectExtra (w, mp, StdCleanupExtraProc); 
  SetObjectExtra (w, mp, NULL);
  mp->form = (ForM) w;
  mp->input_entityID = bfp->input_entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  mp->constraints = FilterGroup (h, TRUE, FALSE, TRUE, FALSE, FALSE, "Where feature text"); 
  mp->accept_cancel = AcceptCancelDialog (h, LocusTagToOldLocusTagAction, NULL, 
                                          LocusTagToOldLocusTagActionClear,
                                          LocusTagToOldLocusTagActionClearText,
                                          (Pointer)mp, w);  
  AlignObjects (ALIGN_CENTER, (HANDLE) mp->constraints,
                              (HANDLE) mp->accept_cancel,
                              NULL);
  Show (w);   
}

static CharPtr GetLastLineageFromString (CharPtr lineage)
{
  CharPtr last_semicolon, lineage_start, penultimate_semicolon;
  
  
  if (StringHasNoText (lineage))
  {
    return lineage;
  }
  
  last_semicolon = StringRChr (lineage, ';');
  if (last_semicolon == NULL)
  {
    return lineage;
  }
  
  /* skip past semicolon */
  lineage_start = last_semicolon + 1;

  /* skip over whitespace */
  lineage_start += StringSpn (lineage_start, " \t");
  
  if (StringICmp (lineage_start, "environmental samples") == 0 
      && last_semicolon != lineage)
  {
    *last_semicolon = 0;
    penultimate_semicolon = StringRChr (lineage, ';');
    if (penultimate_semicolon == NULL)
    {
      lineage_start = lineage;
    }
    else
    {
      lineage_start = penultimate_semicolon + 1;
    }
    lineage_start += StringSpn (lineage_start, " \t");
  }
  return lineage_start;
}

static void ExportLastLineageCallback (BioseqPtr bsp, Pointer userdata)
{
  LogInfoPtr        lip;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  SeqIdPtr          sip;
  Char              tmp[128];
  CharPtr           val_str = NULL;
  ValNode           lineage_field;
  CharPtr           last_lineage;
  Boolean           found_text = FALSE;
  
  if (bsp == NULL || userdata == NULL || ISA_aa (bsp->mol))
  {
    return;
  }
  
  lip = (LogInfoPtr) userdata;
  
  sip = SeqIdFindBest (bsp->id, 0);
  SeqIdWrite (sip, tmp, PRINTID_REPORT, sizeof (tmp));
  
  fprintf (lip->fp, "%s", tmp);
  
  lineage_field.next = NULL;
  lineage_field.choice = 1;
  lineage_field.data.ptrvalue = "Lineage";

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  while (sdp != NULL)
  {
    val_str = GetSourceQualDescrString (sdp, &lineage_field, NULL);
    last_lineage = GetLastLineageFromString (val_str);
    if (last_lineage != NULL)
    {
      if (found_text)
      {
        fprintf (lip->fp, ";%s", last_lineage);
      }
      else
      {
        fprintf (lip->fp, "\t%s", last_lineage);
      }
      val_str = MemFree (val_str);
      found_text = TRUE;
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, 0, &dcontext);
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
  while (sfp != NULL)
  {
    val_str = GetSourceQualFeatureString (sfp, &lineage_field, NULL);
    last_lineage = GetLastLineageFromString (val_str);
    if (last_lineage != NULL)
    {
      if (found_text)
      {
        fprintf (lip->fp, ";%s", last_lineage);
      }
      else
      {
        fprintf (lip->fp, "\t%s", last_lineage);
      }
      val_str = MemFree (val_str);
      found_text = TRUE;
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_BIOSRC, 0, &fcontext);
  }
  fprintf (lip->fp, "\n");
  lip->data_in_log = TRUE;
}

extern void ExportLastLineage (IteM i)
{
  BaseFormPtr         bfp;
  SeqEntryPtr         sep;
  LogInfoPtr          lip;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  
  lip = OpenLog ("Lineage");
  WatchCursor ();
  Update ();
  VisitBioseqsInSep (sep, lip, ExportLastLineageCallback);
  CloseLog (lip);
  lip = FreeLog (lip);    
  ArrowCursor ();
  Update ();
}

typedef struct combinecds
{
  FORM_MESSAGE_BLOCK
  GrouP               action_choice_grp;
  GrouP               new_cds_grp;
  ButtoN              cover_sequence_btn;
  DialoG              string_constraint_dlg;
  
  Int2                action_choice;
  Boolean             cover_sequence;
  StringConstraintXPtr string_constraint;
  LogInfoPtr          lip;
  
  /* for string constraint */
  ObjectHasStringData ohsd;
  ObjMgrTypePtr       omtp;
  AsnIoPtr            aip;
  
  /* for finding location for master CDS */
  Int4                left_end;
  Int4                right_end;
  Boolean             partial5;
  Boolean             partial3;  
  Uint1               strand;
  SeqFeatPtr          first_cds;
  
  /* for product name for master CDS */
  GrouP               prod_name_choice_grp;
  TexT                prod_name_txt;
} CombineCDSData, PNTR CombineCDSPtr;


static void GetCombinedCDSLocationCallback (SeqFeatPtr sfp, Pointer userdata)
{
  CombineCDSPtr     ccp;
  Int4              start, stop;
  Boolean           partial5, partial3;
  Uint1             strand;
  SeqMgrFeatContext fcontext;

  if (sfp == NULL 
      || sfp->data.choice != SEQFEAT_CDREGION 
      || sfp->location == NULL
      || userdata == NULL)
  {
    return;
  }
  
  ccp = (CombineCDSPtr) userdata;
  
  if (ccp->ohsd.scp != NULL
      && !DoesObjectMatchStringConstraint (ccp->omtp, ccp->aip, sfp, &(ccp->ohsd)))
  {
    sfp = SeqMgrGetDesiredFeature (ccp->input_entityID, NULL, 0, 0, sfp, &fcontext);
    if (sfp == NULL || !DoesStringMatchConstraintX (fcontext.label, ccp->ohsd.scp))
    {
      return;
    }
  }
  
  start = SeqLocStart (sfp->location);
  stop = SeqLocStop (sfp->location);
  strand = SeqLocStrand (sfp->location);
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  if (start < ccp->left_end)
  {
    ccp->first_cds = sfp;
    ccp->left_end = start;
    if (strand == Seq_strand_minus)
    {
      ccp->partial3 = partial3;
    }
    else
    {
      ccp->partial5 = partial5;
    }
    ccp->strand = strand;
  }
  
  if (stop > ccp->right_end)
  {
    ccp->right_end = stop;
    if (strand == Seq_strand_minus)
    {
      ccp->partial5 = partial5;
    }
    else
    {
      ccp->partial3 = partial3;
    }
    ccp->strand = strand;
  }
}


static void ApplyProductName (CombineCDSPtr ccp, SeqFeatPtr new_cds)
{
  BioseqPtr         first_prot_bsp, new_prot_bsp;
  SeqFeatPtr        first_prot = NULL;
  ProtRefPtr        prp;
  CharPtr           product_name = NULL;
  
  if (ccp == NULL || new_cds == NULL)
  {
    return;
  }

  if (GetValue (ccp->prod_name_choice_grp) == 1)
  {
    first_prot_bsp = BioseqFindFromSeqLoc (ccp->first_cds->product);
    if (first_prot_bsp != NULL)
    {
      first_prot = GetProtFeature (first_prot_bsp);
      if (first_prot != NULL)
      {
        prp = (ProtRefPtr) first_prot->data.value.ptrvalue;
        if (prp != NULL && prp->name != NULL 
            && ! StringHasNoText (prp->name->data.ptrvalue))
        {
          product_name = StringSave (prp->name->data.ptrvalue);
        }
      }
    }
  }
  else
  {
    product_name = SaveStringFromText (ccp->prod_name_txt);
    if (StringHasNoText (product_name))
    {
      product_name = MemFree (product_name);
    }
  }
  
  if (product_name != NULL)
  {
    new_prot_bsp = BioseqFindFromSeqLoc (new_cds->product);
    if (new_prot_bsp != NULL)
    {
      first_prot = GetProtFeature (new_prot_bsp);
      if (first_prot != NULL)
      {
        prp = (ProtRefPtr) first_prot->data.value.ptrvalue;
        if (prp != NULL)
        {
          ValNodeAddPointer (&(prp->name), 0, product_name);
          product_name = NULL;
        }
      }
    }
  }
  product_name = MemFree (product_name);
}

static SeqLocPtr GetMasterCDSLocation (BioseqPtr bsp, CombineCDSPtr ccp)
{
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  SeqLocPtr         cover_loc;
  
  if (bsp == NULL || ! ISA_na (bsp->mol) || ccp == NULL)
  {
    return NULL;
  }
  
  /* create location for master CDS */
  ccp->left_end = bsp->length - 1;
  ccp->right_end = 0;
  ccp->partial5 = FALSE;
  ccp->partial3 = FALSE;
  
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
  if (sfp == NULL)
  {
    return NULL;
  }
  if (ccp->cover_sequence)
  {
    ccp->left_end = 0;
    ccp->right_end = bsp->length - 1;
    ccp->partial5 = FALSE;
    ccp->partial3 = FALSE;
    ccp->first_cds = sfp;
  }
  else
  {
    while (sfp != NULL)
    {
      GetCombinedCDSLocationCallback (sfp, ccp);
      sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
    }
    if (ccp->left_end > ccp->right_end || ccp->first_cds == NULL)
    {
      return NULL;
    }
  }
   
  cover_loc = SeqLocIntNew (ccp->left_end, ccp->right_end, ccp->strand, bsp->id);
  SetSeqLocPartial (cover_loc, ccp->partial5, ccp->partial3);

  return cover_loc;  
}


static void CreateMasterCodingRegions (BioseqPtr bsp, Pointer userdata)
{
  CombineCDSPtr ccp;
  SeqLocPtr     cover_loc;
  SeqFeatPtr    cds;
  
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  ccp = (CombineCDSPtr) userdata;
  
  cover_loc = GetMasterCDSLocation (bsp,  ccp);
  if (cover_loc == NULL)
  {
    return;
  }
  
  cds = CreateNewFeatureOnBioseq (bsp, SEQFEAT_CDREGION, cover_loc);
  cds->partial = (ccp->partial5 | ccp->partial3);
  
  cds->data.value.ptrvalue = AsnIoMemCopy (ccp->first_cds->data.value.ptrvalue,
                                            (AsnReadFunc) CdRegionAsnRead,
                                            (AsnWriteFunc) CdRegionAsnWrite);
  RetranslateOneCDS (cds, ccp->first_cds->idx.entityID, TRUE, FALSE);

  ApplyProductName (ccp, cds);    
}


static Boolean ConvertCDSToMatPeptideNoOverlap (SeqFeatPtr sfp, LogInfoPtr lip)
{
  SeqMgrFeatContext fcontext;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) {
    return FALSE;
  }
  if (!CreateMatPeptideFromCDS (sfp)) {
    SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &fcontext);
    if (lip != NULL && lip->fp != NULL) {
      fprintf (lip->fp, "Coding region %s has no product sequence\n", fcontext.label);
      lip->data_in_log = TRUE;
    } else {
      Message (MSG_ERROR, "Coding region %s has no product sequence\n", fcontext.label);
    }
    return FALSE;
  } else {
    return TRUE;
  }
}


static void ConvertCDSToMatPeptideNoOverlapCallback (BioseqPtr bsp, Pointer userdata)
{
  LogInfoPtr lip;
  SeqMgrFeatContext fcontext;
  SeqFeatPtr sfp;

  if (bsp == NULL || ISA_aa (bsp->mol)) {
    return;
  }

  lip = (LogInfoPtr) userdata;
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext)) {
    ConvertCDSToMatPeptideNoOverlap (sfp, lip);
  }
}


static void ConvertCDSToMatPeptideOverlap (SeqFeatPtr sfp, SeqFeatPtr top_cds, LogInfoPtr lip)
{
  SeqMgrFeatContext gene_context;
  SeqFeatPtr gene;

  if (sfp == NULL || top_cds == NULL || sfp->data.choice != SEQFEAT_CDREGION || top_cds->data.choice != SEQFEAT_CDREGION) {
    return;
  }

  if (!ConvertCDSToMatPeptideForOverlappingCDS (sfp, top_cds, TRUE)) {
    if (lip != NULL && lip->fp != NULL)
    {
      fprintf (lip->fp, "Invalid coordinates for mat_peptide conversion\n");
      lip->data_in_log = TRUE;
    }
    else
    {
      Message (MSG_ERROR, "Invalid coordinates for mat_peptide conversion");
    }
  } else {
    /* Mark gene with same location for deletion */
    gene = SeqMgrGetOverlappingGene (sfp->location, &gene_context);
    if (gene != NULL && SeqLocCompare (gene->location, sfp->location) == SLC_A_EQ_B)
    {
      gene->idx.deleteme = 1;
    }
  }
}


static void ConvertInnerCDSsToMatPeptidesCallback (BioseqPtr bsp, Pointer userdata)
{
  LogInfoPtr        lip;
  SeqFeatPtr        sfp = NULL;
  ValNodePtr        top_level_cds_list = NULL;
  SeqMgrFeatContext context;
  ValNodePtr        vnp;
  SeqFeatPtr        top_cds;
  
  if (bsp == NULL || ! ISA_na (bsp->mol)) return;
  lip = (LogInfoPtr) userdata;
  
  sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &context);
  while (sfp != NULL)
  {
    top_cds = NULL;
    for (vnp = top_level_cds_list;
         vnp != NULL && top_cds == NULL; 
         vnp = vnp->next)
    {
      top_cds = (SeqFeatPtr) vnp->data.ptrvalue;
      if (top_cds != NULL)
      {
        if (SeqLocCompare (top_cds->location, sfp->location) != SLC_B_IN_A)
        {
          top_cds = NULL;
        }
      }
    }

    if (top_cds == NULL)
    {
      /* add to list of top level CDSs */
      ValNodeAddPointer (&top_level_cds_list, 0, sfp);
    }
    else
    {
      ConvertCDSToMatPeptideOverlap (sfp, top_cds, lip);
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &context);
  }
  ValNodeFree (top_level_cds_list);
}


static void DoCombineCDS (ButtoN b)
{
  CombineCDSPtr      ccp;
  SeqEntryPtr        sep, old_scope;
  ObjMgrPtr          omp;
  AsnExpOptPtr       aeop;

  ccp = (CombineCDSPtr) GetObjectExtra (b);
  if (ccp == NULL)
  {
    return;
  }
  
  sep = GetTopSeqEntryForEntityID (ccp->input_entityID);
  if (NULL == sep)
    return;
  
  old_scope = SeqEntrySetScope (sep);
  
  omp = ObjMgrGet ();
  if (omp == NULL) return;

  ccp->action_choice = GetValue (ccp->action_choice_grp);
  ccp->cover_sequence = GetStatus (ccp->cover_sequence_btn);
  ccp->string_constraint = DialogToPointer (ccp->string_constraint_dlg);

  WatchCursor ();
  Update ();
  
  ccp->lip = OpenLog ("Combine CDS Features");
  
  if (ccp->string_constraint != NULL)
  {
    ccp->omtp = ObjMgrTypeFind (omp, OBJ_SEQFEAT, NULL, NULL);
    ccp->aip = AsnIoNullOpen ();
    ccp->ohsd.scp = ccp->string_constraint;
    aeop = AsnExpOptNew (ccp->aip, NULL, NULL, AsnWriteConstraintCallBack);
    if (aeop != NULL) {
      aeop->user_data = (Pointer) &(ccp->ohsd);
    }
  }
  else
  {
    ccp->omtp = NULL;
    ccp->aip = NULL;
    ccp->ohsd.scp = NULL;
    aeop = NULL;
  }

  if (ccp->action_choice == 2) {
    VisitBioseqsInSep (sep, ccp, CreateMasterCodingRegions);
    SeqMgrIndexFeatures (ccp->input_entityID, NULL);
  }
  
  if (ccp->action_choice == 1 || ccp->action_choice == 2) {
    VisitBioseqsInSep (sep, ccp->lip, ConvertInnerCDSsToMatPeptidesCallback);
  } else {
    VisitBioseqsInSep (sep, ccp->lip, ConvertCDSToMatPeptideNoOverlapCallback);
  }

  AsnIoClose (ccp->aip);
  SeqEntrySetScope (old_scope);
  
  CloseLog (ccp->lip);
  FreeLog (ccp->lip);  
  ccp->string_constraint = StringConstraintXFree (ccp->string_constraint);
  Remove (ccp->form);
  DeleteMarkedObjects (ccp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (ccp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ccp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
}

static void EnableNewCDSControls (GrouP g)
{
  CombineCDSPtr ccp;

  ccp = (CombineCDSPtr) GetObjectExtra (g);
  if (ccp == NULL) return;
  
  if (GetValue (ccp->action_choice_grp) == 2)
  {
    Enable (ccp->new_cds_grp);
  }
  else
  {
    Disable (ccp->new_cds_grp);
  }
}

static void EnableProdName (GrouP g)
{
  CombineCDSPtr ccp;

  ccp = (CombineCDSPtr) GetObjectExtra (g);
  if (ccp == NULL) return;
  
  if (GetValue (ccp->prod_name_choice_grp) == 1)
  {
    Disable (ccp->prod_name_txt);
  }
  else
  {
    Enable (ccp->prod_name_txt);
  }
}

extern void CombineMultipleCDS (IteM i)
{
  BaseFormPtr   bfp;
  CombineCDSPtr ccp;
  WindoW        w;
  GrouP         h, c;
  ButtoN        b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  ccp = (CombineCDSPtr) MemNew (sizeof (CombineCDSData));
  if (ccp == NULL) return;
  
  w = FixedWindow (-50, -33, -10, -10, "Convert CDS to Mat-peptide", StdCloseWindowProc);
  SetObjectExtra (w, ccp, StdCleanupExtraProc);
  ccp->form = (ForM) w;
  ccp->input_entityID = bfp->input_entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  ccp->action_choice_grp = HiddenGroup (h, 0, 3, EnableNewCDSControls);
  RadioButton (ccp->action_choice_grp, "Convert inner CDSs to mat_peptides");
  RadioButton (ccp->action_choice_grp, "Merge multiple CDSs into one CDS and convert inner CDSs to mat_peptides");
  RadioButton (ccp->action_choice_grp, "Create mat-peptide from protein on each CDS");
  SetValue (ccp->action_choice_grp, 1);
  SetObjectExtra (ccp->action_choice_grp, ccp, NULL);
  
  ccp->new_cds_grp = HiddenGroup (h, -1, 0, NULL);
  ccp->cover_sequence_btn = CheckBox (ccp->new_cds_grp, "New CDS should cover entire sequence", NULL);
  SetStatus (ccp->cover_sequence_btn, FALSE);
  
  ccp->prod_name_choice_grp = HiddenGroup (ccp->new_cds_grp, 0, 2, EnableProdName);
  SetObjectExtra (ccp->prod_name_choice_grp, ccp, NULL);
  RadioButton (ccp->prod_name_choice_grp, "Use product name from first CDS");
  RadioButton (ccp->prod_name_choice_grp, "Use this product name:");
  StaticPrompt (ccp->prod_name_choice_grp, "", 0, 0, programFont, 'l');
  ccp->prod_name_txt = DialogText (ccp->prod_name_choice_grp, " ", 10, NULL);
  SetValue (ccp->prod_name_choice_grp, 1);
  Disable (ccp->prod_name_txt);
  
  ccp->string_constraint_dlg = StringConstraintDialogX (ccp->new_cds_grp, "Where CDS field", TRUE);
  Disable (ccp->new_cds_grp);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) ccp->cover_sequence_btn,
                              (HANDLE) ccp->prod_name_choice_grp,
                              (HANDLE) ccp->string_constraint_dlg,
                              NULL);
  
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Accept", DoCombineCDS);
  SetObjectExtra (b, ccp, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) ccp->action_choice_grp,
                              (HANDLE) ccp->new_cds_grp,
                              (HANDLE) c,
                              NULL);
  RealizeWindow (w);
  Show (w);
  Select (w);
  Update ();
}


static void MoveToProteinSequenceCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  SeqLocPtr         prot_loc;
  SeqFeatPtr        overlapping_cds, new_sfp;
  SeqMgrFeatContext cds_context, fcontext;
  Int4              offset;
  CdRegionPtr       crp;
  BioseqPtr         prot_bsp;
  LogInfoPtr        lip;
  
  if (sfp == NULL) {
    return;
  }
  
  lip = (LogInfoPtr) userdata;
  
  overlapping_cds = SeqMgrGetOverlappingCDS(sfp->location, &cds_context);
  if (overlapping_cds == NULL) {
    SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &fcontext);
    if (lip != NULL && lip->fp != NULL)
    {
      fprintf (lip->fp, "No overlapping coding region for %s\n", fcontext.label);
      lip->data_in_log = TRUE;
    }
    else
    {
      Message (MSG_ERROR, "No overlapping coding region for %s\n", fcontext.label);
    }
    return;
  }
  
  prot_bsp = BioseqFindFromSeqLoc(overlapping_cds->product);
  if (prot_bsp == NULL) {
    if (lip != NULL && lip->fp != NULL)
    {
      fprintf (lip->fp, "No protein product for %s\n", cds_context.label);
      lip->data_in_log = TRUE;
    }
    else
    {
      Message (MSG_ERROR, "No protein product for %s\n", cds_context.label);
    }
    return;
  }
  
  offset = cds_context.left;
  crp = (CdRegionPtr) overlapping_cds->data.value.ptrvalue;
  if (crp != NULL) {
    if (crp->frame == 2)
    {
      offset += 1;
    }
    else if (crp->frame == 3)
    {
      offset += 2;
    }
  }
  
  prot_loc = BuildProtLoc (overlapping_cds, sfp->location, NULL);
  if (prot_loc == NULL)
  {
    SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &fcontext);
    if (lip != NULL && lip->fp != NULL)
    {
      fprintf (lip->fp, "Unable to translate coordinates for %s\n", fcontext.label);
      lip->data_in_log = TRUE;
    }
    else
    {
      Message (MSG_ERROR, "Unable to translate coordinates for %s\n", fcontext.label);
    }
  } 
  else 
  {
    new_sfp = CreateNewFeatureOnBioseq (prot_bsp, sfp->data.choice, prot_loc);
    if (new_sfp != NULL)
    {
      new_sfp->data.value.ptrvalue = sfp->data.value.ptrvalue;
      sfp->data.value.ptrvalue = NULL;
      sfp->idx.deleteme = TRUE;
    }
  }
  
}


typedef struct maptoprotein {
  FORM_MESSAGE_BLOCK
  DialoG feature_selection;
  DialoG constraints;
} MapToProteinData, PNTR MapToProteinPtr;


static void MapToProteinSequence(ButtoN b)
{
  SeqEntryPtr     sep;
  MapToProteinPtr mtpp;
  FilterSetPtr    fsp;
  LogInfoPtr      lip;
  ValNodePtr      feature_type_list, vnp;
  Uint1           feat_def_choice;
  
  mtpp = (MapToProteinPtr) GetObjectExtra (b);
  if (mtpp == NULL) {
    return;
  }
  
  feature_type_list = (ValNodePtr) DialogToPointer (mtpp->feature_selection);
  if (feature_type_list != NULL)
  {
    sep = GetTopSeqEntryForEntityID (mtpp->input_entityID); 
    fsp = DialogToPointer (mtpp->constraints);
    lip = OpenLog ("Map to Protein Sequence");
    
    for (vnp = feature_type_list; vnp != NULL; vnp = vnp->next)
    {
      feat_def_choice = vnp->choice;
      if (feat_def_choice == 255)
      {
        feat_def_choice = 0;
      }
      OperateOnSeqEntryConstrainedObjects (sep, fsp, MoveToProteinSequenceCallback, 
                                           NULL, 0, feat_def_choice, 0, lip);
    }
    fsp = FilterSetFree (fsp);
    CloseLog (lip);
    lip = FreeLog (lip);
    DeleteMarkedObjects (mtpp->input_entityID, 0, NULL);
    ObjMgrSetDirtyFlag (mtpp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, mtpp->input_entityID, 0, 0);
    Update ();
    Remove (mtpp->form);
    ValNodeFree (feature_type_list);
  }
}


extern void MapFeaturesToProteinSequence(IteM i)
{
  ButtoN             b;
  BaseFormPtr        bfp;
  GrouP              c;
  GrouP              g;
  WindoW             w;
  MapToProteinPtr    mtpp;
  PrompT             p;
  SeqEntryPtr        sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  
  mtpp = (MapToProteinPtr) MemNew (sizeof (MapToProteinData));
  if (mtpp == NULL) return;
  mtpp->input_entityID = bfp->input_entityID;
  sep = GetTopSeqEntryForEntityID(bfp->input_entityID);
  w = MovableModalWindow (-50, -33, -10, -10,
                          "Map to Protein Sequence",
                          StdCloseWindowProc);
  SetObjectExtra (w, mtpp, StdCleanupFormProc);
  mtpp->form = (ForM) w;

  g = HiddenGroup (w, -1, 0, NULL);
  p = StaticPrompt (g, "Move features to protein of type", 0, dialogTextHeight, systemFont, 'l');

  mtpp->feature_selection = FeatureSelectionDialogEx (g, TRUE, sep, NULL, NULL);
  mtpp->constraints = FilterGroup (g, TRUE, FALSE, FALSE, FALSE, FALSE, "Where feature text");
  AlignObjects (ALIGN_CENTER, (HANDLE) p, (HANDLE) mtpp->feature_selection, (HANDLE) mtpp->constraints, NULL);

  c = HiddenGroup (w, 2, 0, NULL);
  b = DefaultButton (c, "Accept", MapToProteinSequence);
  SetObjectExtra (b, mtpp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Select (w);
  Update ();
}


typedef struct uniquevaluecollection {
  ValNodePtr               requested_field;
  GetFeatureFieldString    fieldstring_func;
  GetDescriptorFieldString descrstring_func;
  FreeValNodeProc          free_vn_proc;
  CopyValNodeDataProc      copy_vn_proc;
  ValNodePtr               value_lists;
} UniqueValueCollectionData, PNTR UniqueValueCollectionPtr;


static ClickableItemPtr FindExistingCategory (CharPtr description, ValNodePtr category_list)
{
  ClickableItemPtr cip;

  while (category_list != NULL) {
    cip = (ClickableItemPtr) category_list->data.ptrvalue;
    if (cip != NULL && 
        (StringCmp (description, cip->description) == 0
         || (StringHasNoText (description) && StringHasNoText (cip->description)))) {
      return cip;
    }
    category_list = category_list->next;
  }
  return NULL;
}


static ValNodePtr DataInValNodeList (ValNodePtr list, Uint1 choice, Pointer ptrvalue)
{
  while (list != NULL) {
    if (choice == list->choice && ptrvalue == list->data.ptrvalue) {
      return list;
    }
    list = list->next;
  }
  return NULL;
}


static void CombineClickableItemLists (ValNodePtr PNTR list1, ValNodePtr list2)
{
  ValNodePtr add_vnp, prev_vnp = NULL, next_vnp;

  if (*list1 == NULL) {
    *list1 = list2;
    return;
  } else if (list2 == NULL) {
    return;
  }

  add_vnp = list2;
  while (add_vnp != NULL) {
    next_vnp = add_vnp->next;
    if (DataInValNodeList (*list1, add_vnp->choice, add_vnp->data.ptrvalue)) {
      if (prev_vnp == NULL) {
        list2 = next_vnp;
      } else {
        prev_vnp->next = next_vnp;
      }
      add_vnp->next = NULL;
      ValNodeFree (add_vnp);
    } else {
      prev_vnp = add_vnp;
    }
    add_vnp = next_vnp;
  }
  ValNodeLink (list1, list2); 
}


static void CombineClickableLists (ValNodePtr PNTR list1, ValNodePtr list2)
{
  ClickableItemPtr cip1, cip2;
  ValNodePtr vnp;
  for (vnp = list2; vnp != NULL; vnp = vnp->next) {
    cip2 = vnp->data.ptrvalue;
    if (cip2 != NULL) {
      cip1 = FindExistingCategory (cip2->description, *list1);
      if (cip1 == NULL) {
        ValNodeAddPointer (list1, 0, cip2);
        vnp->data.ptrvalue = NULL;
      } else {
        CombineClickableItemLists (&(cip1->item_list), cip2->item_list);
        cip2->item_list = NULL;
      }
    }
  }
  list2 = FreeClickableList (list2);
}


static void GetUniqueValuesFeatureCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  UniqueValueCollectionPtr uvcp;
  CharPtr      str;
  ValNodePtr   check_vnp;
  ClickableItemPtr cip;
  
  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  
  uvcp = (UniqueValueCollectionPtr) userdata;
  if (uvcp->fieldstring_func == NULL)
  {
    return;
  }
  
  str = uvcp->fieldstring_func (sfp, uvcp->requested_field, NULL);
  cip = FindExistingCategory (str, uvcp->value_lists);
  if (cip == NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->description = str;
    ValNodeAddPointer (&cip->item_list, OBJ_SEQFEAT, sfp);
    ValNodeAddPointer (&(uvcp->value_lists), 0, cip);
  } else {
    str = MemFree (str);
    /* make sure feature isn't already in the list */
    for (check_vnp = cip->item_list; check_vnp != NULL; check_vnp = check_vnp->next) {
      if (check_vnp->choice == OBJ_SEQFEAT && sfp == check_vnp->data.ptrvalue) break;
    }
    if (check_vnp == NULL) {
      ValNodeAddPointer (&cip->item_list, OBJ_SEQFEAT, sfp);
    }
  }
}

static void GetUniqueValuesDescriptorCallback (SeqDescrPtr sdp, Pointer userdata, FilterSetPtr fsp)
{
  UniqueValueCollectionPtr uvcp;
  ClickableItemPtr         cip;
  CharPtr      str;
  ValNodePtr   check_vnp;
  
  if (sdp == NULL || userdata == NULL)
  {
    return;
  }
  
  uvcp = (UniqueValueCollectionPtr) userdata;
  if (uvcp->descrstring_func == NULL)
  {
    return;
  }

  str = uvcp->descrstring_func (sdp, uvcp->requested_field, NULL);

  cip = FindExistingCategory (str, uvcp->value_lists);

  if (cip == NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->description = str;
    ValNodeAddPointer (&cip->item_list, OBJ_SEQDESC, sdp);
    ValNodeAddPointer (&(uvcp->value_lists), 0, cip);
  } else {
    str = MemFree (str);
    /* make sure feature isn't already in the list */
    for (check_vnp = cip->item_list; check_vnp != NULL; check_vnp = check_vnp->next) {
      if (check_vnp->choice == OBJ_SEQDESC && sdp == check_vnp->data.ptrvalue) break;
    }
    if (check_vnp == NULL) {
      ValNodeAddPointer (&cip->item_list, OBJ_SEQDESC, sdp);
    }
  }
}


static ValNodePtr CheckForUniqueValuesInSeqEntry
(SeqEntryPtr              sep,
 ValNodePtr               requested_field,
 GetFeatureFieldString    fieldstring_func,
 GetDescriptorFieldString descrstring_func,
 FreeValNodeProc          free_vn_proc,
 CopyValNodeDataProc      copy_vn_proc,
 FilterSetPtr             fsp,
 Uint1                    seqfeat_choice,
 Uint1                    featdef_choice,
 Uint1                    descr_choice)
{
  UniqueValueCollectionData uvcd;
  
  uvcd.fieldstring_func = fieldstring_func;
  uvcd.descrstring_func = descrstring_func;
  uvcd.free_vn_proc = free_vn_proc;
  uvcd.copy_vn_proc = copy_vn_proc;
  uvcd.requested_field = (copy_vn_proc) (requested_field);
  uvcd.value_lists = NULL;
  OperateOnSeqEntryConstrainedObjects (sep, fsp, GetUniqueValuesFeatureCallback, 
                                         GetUniqueValuesDescriptorCallback,
                                         seqfeat_choice, featdef_choice,
                                         descr_choice, &uvcd);
  return uvcd.value_lists;  
}


static ValNodePtr
GetUniqueValueListForSeqEntry
(SeqEntryPtr   sep,
 Uint2         entityID, 
 ParseFieldPtr dst_field_data,
 FilterSetPtr  fsp)
{
  ValNodePtr   unique_values = NULL;
  ValNodePtr   requested_field = NULL, vnp;
  SeqEntryPtr  orig_sep;

  if (sep == NULL || dst_field_data == NULL)
  {
    return NULL;
  }
  
  switch (dst_field_data->parse_field_type)
  {
    case PARSE_FIELD_SOURCE_QUAL :
      unique_values = CheckForUniqueValuesInSeqEntry (sep, 
                                            dst_field_data->feature_field, 
                                            GetSourceQualFeatureString,
                                            GetSourceQualDescrString,
                                            ValNodeSimpleDataFree,
                                            SourceQualValNodeDataCopy,
                                            fsp,
                                            SEQFEAT_BIOSRC, 0,
                                            Seq_descr_source);
      break;
    case PARSE_FIELD_DEFLINE:
      requested_field = ValNodeNew (NULL);
      requested_field->data.intvalue = Seq_descr_title;
      unique_values = CheckForUniqueValuesInSeqEntry (sep, requested_field, 
                                            NULL,
                                            GetStringFromStringDescriptor,
                                            NULL, IntValNodeCopy,
                                            fsp, 0, 0, Seq_descr_title);
      requested_field = ValNodeFree (requested_field);
      break;
    case PARSE_FIELD_BIOSRC_STRING:
      unique_values = CheckForUniqueValuesInSeqEntry (sep, 
                                            dst_field_data->feature_field, 
                                            GetSourceFeatureString,
                                            GetSourceDescriptorString,
                                            NULL, IntValNodeCopy,
                                            fsp,
                                            SEQFEAT_BIOSRC, 0,
                                            Seq_descr_source);
      break;
    case PARSE_FIELD_DBXREF:
      unique_values = CheckForUniqueValuesInSeqEntry (sep, dst_field_data->feature_field,
                                            GetBioSourceFeatureDbxrefString,
                                            GetBioSourceDescriptorDbxrefString,
                                            ValNodeSimpleDataFree,
                                            ValNodeStringCopy,
                                            fsp,
                                            SEQFEAT_BIOSRC, 0, Seq_descr_source);
      break;
    case PARSE_FIELD_GENE_FIELD:
      unique_values = CheckForUniqueValuesInSeqEntry (sep, 
                                            dst_field_data->feature_field, 
                                            GetGeneFieldString,
                                            NULL,
                                            NULL, IntValNodeCopy,
                                            fsp,
                                            SEQFEAT_GENE, 0, 0);
      break;
    case PARSE_FIELD_RNA_FIELD:
      unique_values = CheckForUniqueValuesInSeqEntry (sep, 
                                            dst_field_data->feature_field, 
                                            GetRNAFieldString,
                                            NULL,
                                            NULL, IntValNodeCopy,
                                            fsp,
                                            SEQFEAT_RNA, 
                                            dst_field_data->feature_subtype == NULL ? 0 : dst_field_data->feature_subtype->data.intvalue,
                                            0);
      break;
    case PARSE_FIELD_CDS_COMMENT:
      unique_values = CheckForUniqueValuesInSeqEntry (sep, 
                                            NULL, 
                                            GetCDSComment,
                                            NULL,
                                            NULL, IntValNodeCopy,
                                            fsp,
                                            SEQFEAT_CDREGION, FEATDEF_CDS, 0);
      break;
    case PARSE_FIELD_COMMENT_DESC:
      requested_field = ValNodeNew (NULL);
      requested_field->data.intvalue = Seq_descr_comment;

      unique_values = CheckForUniqueValuesInSeqEntry (sep, requested_field, 
                                            NULL,
                                            GetStringFromStringDescriptor,
                                            NULL, IntValNodeCopy,
                                            fsp, 0, 0, Seq_descr_comment);
      requested_field = ValNodeFree (requested_field);
      break;
    case PARSE_FIELD_PROTEIN_FIELD:
      unique_values = CheckForUniqueValuesInSeqEntry (sep, 
                                            dst_field_data->feature_field, 
                                            GetProteinFieldString,
                                            NULL,
                                            NULL, IntValNodeCopy,
                                            fsp,
                                            SEQFEAT_PROT, 0, 0);
      break;
    case PARSE_FIELD_IMPORT_QUAL:
      orig_sep = sep;
      for (vnp = dst_field_data->feature_subtype; vnp != NULL; vnp = vnp->next)
      {
        CombineClickableLists (&unique_values, 
                               CheckForUniqueValuesInSeqEntry (sep, 
                                                               dst_field_data->feature_field, 
                                                               GetGBQualString,
                                                               NULL,
                                                               NULL, IntValNodeCopy,
                                                               fsp, 0, vnp->choice, 0));
        /* if we are also looking at features other than mat_peptides, use
         * original SeqEntry */
        sep = orig_sep;
      }
      break;
    case PARSE_FIELD_FEATURE_NOTE:
      orig_sep = sep;      
      for (vnp = dst_field_data->feature_field; vnp != NULL; vnp = vnp->next)
      {
        CombineClickableLists (&unique_values, 
                               CheckForUniqueValuesInSeqEntry (sep, vnp, 
                                                               GetFeatureNote,
                                                               NULL,
                                                               ValNodeSimpleDataFree,
                                                               ValNodeStringCopy,
                                                               fsp, 0, vnp->choice, 0));
        /* if we are also looking at features other than mat_peptides, use
         * original SeqEntry */
        sep = orig_sep;
      }
      break;
  }
  return unique_values;  
}


static ValNodePtr GetSeqEntryListForItem (ValNodePtr vnp)
{
  SeqFeatPtr       sfp;
  SeqDescrPtr      sdp;
  ObjValNodePtr    ovp;
  SeqEntryPtr      sep = NULL;
  BioseqPtr        bsp;
  BioseqSetPtr     bssp;
  ValNodePtr       sep_list = NULL;

  if (vnp == NULL || vnp->data.ptrvalue == NULL) return NULL;

  if (vnp->choice == OBJ_SEQFEAT) {
    sfp = vnp->data.ptrvalue;
    sep = GetBestTopParentForData (sfp->idx.entityID, BioseqFindFromSeqLoc (sfp->location));
    ValNodeAddPointer (&sep_list, OBJ_SEQENTRY, sep);
  } else if (vnp->choice == OBJ_SEQDESC) {
    sdp = vnp->data.ptrvalue;
    if (sdp->extended != 0) {
      ovp = (ObjValNodePtr) sdp;
      if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
        bssp = (BioseqSetPtr) ovp->idx.parentptr;
        if (bssp != NULL) {
          /* nuc_prot sets, segsets, and members of parts should travel together */
          if (bssp->_class == BioseqseqSet_class_nuc_prot) {
            sep = SeqMgrGetSeqEntryForData (bssp);
            if (sep != NULL) {
              ValNodeAddPointer (&sep_list, OBJ_SEQENTRY, sep);
            }
          } else if (bssp->_class == BioseqseqSet_class_segset
                     || bssp->_class == BioseqseqSet_class_parts) {
            sep = bssp->seq_set;
            if (sep != NULL && IS_Bioseq (sep)) {
              bsp = sep->data.ptrvalue;
              sep = GetBestTopParentForData (bsp->idx.entityID, bsp);
              if (sep != NULL) {
                ValNodeAddPointer (&sep_list, OBJ_SEQENTRY, sep);
              }
            }
          } else {            
            sep = bssp->seq_set;
            while (sep != NULL) {
              ValNodeAddPointer (&sep_list, OBJ_SEQENTRY, sep);
              sep = sep->next;
            }
          }
        }
      } else if (ovp->idx.parenttype == OBJ_BIOSEQ) {
        sep = GetBestTopParentForData (ovp->idx.entityID, ovp->idx.parentptr);
        if (sep != NULL) {
          ValNodeAddPointer (&sep_list, OBJ_SEQENTRY, sep);
        }
      }
    }
  } else if (vnp->choice == OBJ_BIOSEQ) {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    sep = GetBestTopParentForData (bsp->idx.entityID, bsp);
    ValNodeAddPointer (&sep_list, OBJ_SEQENTRY, sep);
  } else if (vnp->choice == OBJ_SEQENTRY) {
    sep = vnp->data.ptrvalue;
    ValNodeAddPointer (&sep_list, OBJ_SEQENTRY, sep);
  }
  return sep_list;  
}


static Boolean RemoveItemAndCategory (ClickableItemPtr cip, Uint1 choice, Pointer ptrvalue)
{
  ValNodePtr vnp, prev = NULL;
  Boolean    found = FALSE;
  CharPtr    description = NULL;
  ClickableItemPtr subcat;

  vnp = cip->item_list;
  while (vnp != NULL) {
    if (vnp->choice == choice && vnp->data.ptrvalue == ptrvalue) {
      found = TRUE;
      if (prev == NULL) {
        cip->item_list = vnp->next;
      } else {
        prev->next = vnp->next;
      }
      vnp->next = NULL;
      description = GetDiscrepancyItemText (vnp);
      vnp = ValNodeFree (vnp);
    } else {
      prev = vnp;
      vnp = vnp->next;
    }
  }
  if (found) {
    /* also remove subcategory */
    prev = NULL;
    vnp = cip->subcategories;
    while (vnp != NULL) {
      subcat = (ClickableItemPtr) vnp->data.ptrvalue;
      if (subcat != NULL && StringCmp (subcat->description, description) == 0) {
        if (prev == NULL) {
          cip->subcategories = vnp->next;
        } else {
          prev->next = vnp->next;
        }
        vnp->next = NULL;
        vnp = FreeClickableList (vnp);
      } else {
        prev = vnp;
        vnp = vnp->next;
      }
    }
  }
  return found;
}

static ValNodePtr RearrangeForSeqEntriesInMultipleCategories (ValNodePtr value_lists)
{
  ClickableItemPtr cip1, cip2, existing_cat;
  ValNodePtr       vnp, check_vnp, last_vnp, vnp_item1, prev, vnp_next;
  Boolean          found;
  CharPtr          new_description;

  /* now check for SeqEntries that appear in multiple lists */
  for (vnp = value_lists; vnp != NULL && vnp->next != NULL; vnp = vnp->next) {
    if (vnp->data.ptrvalue == NULL) continue;
    last_vnp = vnp->next;
    while (last_vnp->next != NULL) {
      last_vnp = last_vnp->next;
    }
    cip1 = (ClickableItemPtr) vnp->data.ptrvalue;
    for (vnp_item1 = cip1->item_list; 
         vnp_item1 != NULL;
         vnp_item1 = vnp_next) {
      vnp_next = vnp_item1->next;
      found = FALSE;
      if (vnp_item1->data.ptrvalue != NULL && vnp_item1->choice == OBJ_SEQENTRY) {
        check_vnp = vnp->next;
        while (check_vnp != NULL && !found) {
          if (check_vnp->data.ptrvalue != NULL) {
            cip2 = check_vnp->data.ptrvalue;
            found = RemoveItemAndCategory (cip2, vnp_item1->choice, vnp_item1->data.ptrvalue);
          } 
          if (!found) {
            if (check_vnp == last_vnp) {
              check_vnp = NULL;
            } else {
              check_vnp = check_vnp->next;
            }
          }
        }
      }
      if (found) {
        /* need to take this SeqEntry out of both lists, put in new category */
        new_description = (CharPtr) MemNew (sizeof (Char) * (StringLen (cip1->description)
                                                             + StringLen (cip2->description) + 6));
        if (StringCmp (cip1->description, cip2->description) < 0) {
          sprintf (new_description, "%s AND %s", cip1->description, cip2->description);
        } else {
          sprintf (new_description, "%s AND %s", cip2->description, cip1->description);
        }
        existing_cat = FindExistingCategory (new_description, check_vnp->next);
        if (existing_cat == NULL) {
          /* make new category, add to end of list */
          existing_cat = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
          MemSet (existing_cat, 0, sizeof (ClickableItemData));
          existing_cat->description = new_description;
          ValNodeAddPointer (&value_lists, 0, existing_cat);
        } else {
          new_description = MemFree (new_description);
        }
        /* don't add item twice */
        if (!DataInValNodeList (existing_cat->item_list, vnp_item1->choice, vnp_item1->data.ptrvalue)) {
          ValNodeAddPointer (&(existing_cat->item_list), vnp_item1->choice, vnp_item1->data.ptrvalue);
        }
        /* remove from cip1 */
        RemoveItemAndCategory (cip1, vnp_item1->choice, vnp_item1->data.ptrvalue);
      }  
    }   
  }

  /* remove empty categories */
  for (vnp = value_lists, prev = NULL; vnp != NULL; vnp = vnp_next) {
    vnp_next = vnp->next;
    cip1 = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip1 == NULL || cip1->item_list == NULL) {
      if (prev == NULL) {
        value_lists = vnp->next;
      } else {
        prev->next = vnp->next;
      }
      vnp->next = NULL;
      vnp = FreeClickableList (vnp);
    } else {
      prev = vnp;
    }
  }
  return value_lists;
}


static ValNodePtr ChangeUniqueValueListsToSeqEntryLists (ValNodePtr value_lists)
{
  ValNodePtr vnp, vnp_item, new_item_list, new_vnp, sep_list, sep_vnp;
  ClickableItemPtr cip, subcat;
  SeqEntryPtr      sep;
  CharPtr          description;

  /* First, add SeqEntry for each feature or descriptor to the item list */
  for (vnp = value_lists; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    new_item_list = NULL;
    if (cip != NULL) {
      for (vnp_item = cip->item_list; vnp_item != NULL; vnp_item = vnp_item->next) {
        if (vnp_item->data.ptrvalue == NULL) continue;
        sep_list = GetSeqEntryListForItem (vnp_item);
        for (sep_vnp = sep_list; sep_vnp != NULL; sep_vnp = sep_vnp->next) {
          sep = sep_vnp->data.ptrvalue;
          new_vnp = DataInValNodeList (new_item_list, OBJ_SEQENTRY, sep);
          /* add to list for this category */
          if (new_vnp == NULL) {
            new_vnp = ValNodeAddPointer (&new_item_list, OBJ_SEQENTRY, sep);
          }
          description = GetDiscrepancyItemText(new_vnp);
          subcat = FindExistingCategory (description, cip->subcategories);
          if (subcat == NULL) {
            subcat = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
            MemSet (subcat, 0, sizeof (ClickableItemData));
            subcat->description = description;
            ValNodeAddPointer (&(cip->subcategories), 0, subcat);
          } else {
            description = MemFree (description);
          }
          /* don't add item twice */
          if (!DataInValNodeList (subcat->item_list, vnp_item->choice, vnp_item->data.ptrvalue)) {
            ValNodeAddPointer (&(subcat->item_list), vnp_item->choice, vnp_item->data.ptrvalue);
          }
        }
        sep_list = ValNodeFree (sep_list);
      }
      
      cip->item_list = ValNodeFree (cip->item_list);
      cip->item_list = new_item_list;
    }
  }
  value_lists = RearrangeForSeqEntriesInMultipleCategories (value_lists);
  return value_lists;
}


static void GetMolTypes (SeqDescrPtr sdp, Pointer userdata)
{
  ValNodePtr PNTR   type_list;
  ClickableItemPtr  cip;
  CharPtr           str;
  MolInfoPtr        mip;
  
  if (sdp == NULL || sdp->choice != Seq_descr_molinfo || sdp->data.ptrvalue == NULL || userdata == NULL)
  {
    return;
  }
  
  type_list = (ValNodePtr PNTR) userdata;

  mip = (MolInfoPtr) sdp->data.ptrvalue;
  /* don't look for peptides */
  if (mip->biomol == 8) return;

  str = GetMoleculeTypeName (mip->biomol);

  cip = FindExistingCategory (str, *type_list);

  if (cip == NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->description = StringSave (str);
    ValNodeAddPointer (&cip->item_list, OBJ_SEQDESC, sdp);
    ValNodeAddPointer (type_list, 0, cip);
  } else {
    /* make sure feature isn't already in the list */
    if (DataInValNodeList (cip->item_list, OBJ_SEQDESC, sdp) == NULL) {
      ValNodeAddPointer (&cip->item_list, OBJ_SEQDESC, sdp);
    }
  }
}

static void GetMolClasses (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR   type_list;
  ClickableItemPtr  cip;
  CharPtr           str;
  
  if (bsp == NULL || ISA_aa(bsp->mol) || userdata == NULL)
  {
    return;
  }
  
  type_list = (ValNodePtr PNTR) userdata;

  str = GetMoleculeClassName (bsp->mol);

  cip = FindExistingCategory (str, *type_list);

  if (cip == NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->description = StringSave (str);
    ValNodeAddPointer (&cip->item_list, OBJ_BIOSEQ, bsp);
    ValNodeAddPointer (type_list, 0, cip);
  } else {
    /* make sure feature isn't already in the list */
    if (DataInValNodeList (cip->item_list, OBJ_BIOSEQ, bsp) == NULL) {
      ValNodeAddPointer (&cip->item_list, OBJ_BIOSEQ, bsp);
    }
  }
}

static void GetMolTopologies (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR   type_list;
  ClickableItemPtr  cip;
  CharPtr           str;
  
  if (bsp == NULL || ISA_aa(bsp->mol) || userdata == NULL)
  {
    return;
  }
  
  type_list = (ValNodePtr PNTR) userdata;

  str = GetMoleculeTopologyName (bsp->topology);

  cip = FindExistingCategory (str, *type_list);

  if (cip == NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->description = StringSave (str);
    ValNodeAddPointer (&cip->item_list, OBJ_BIOSEQ, bsp);
    ValNodeAddPointer (type_list, 0, cip);
  } else {
    /* make sure feature isn't already in the list */
    if (DataInValNodeList (cip->item_list, OBJ_BIOSEQ, bsp) == NULL) {
      ValNodeAddPointer (&cip->item_list, OBJ_BIOSEQ, bsp);
    }
  }
}


typedef struct segfeatures {
  Uint2      featdef;
  ValNodePtr feature_list;
} SegFeaturesData, PNTR SegFeaturesPtr;

static void FindFeaturesForSeg (SeqFeatPtr sfp, Pointer userdata)
{
  SegFeaturesPtr fp;

  if (sfp == NULL || userdata == NULL) return;

  fp = (SegFeaturesPtr) userdata;

  if (sfp->idx.subtype == fp->featdef) {
    ValNodeAddPointer (&(fp->feature_list), OBJ_SEQFEAT, sfp);
  }
}

typedef struct segdescriptors {
  Uint2      desc_type;
  ValNodePtr descriptor_list;
} SegDescriptorsData, PNTR SegDescriptorsPtr;

static void FindDescriptorsForSeg (SeqDescrPtr sdp, Pointer userdata)
{
  SegDescriptorsPtr fp;

  if (sdp == NULL || userdata == NULL) return;

  fp = (SegDescriptorsPtr) userdata;

  if (sdp->choice == fp->desc_type) {
    ValNodeAddPointer (&(fp->descriptor_list), OBJ_SEQDESC, sdp);
  }
}


extern ValNodePtr 
CreateSeqEntryListsForUniqueValues
(SeqEntryPtr   sep,
 Uint2         entityID, 
 ParseFieldPtr dst_field_data,
 FilterSetPtr  fsp)
{
  ValNodePtr value_lists;

  value_lists = GetUniqueValueListForSeqEntry (sep, entityID, dst_field_data, fsp);

  value_lists = ChangeUniqueValueListsToSeqEntryLists (value_lists);
  return value_lists;
}

static void AddPubFeatureWithText (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR pub_list;

  if (sfp == NULL || userdata == NULL) return;
  pub_list = (ValNodePtr PNTR) userdata;

  ValNodeAddPointer (pub_list, OBJ_SEQFEAT, sfp);
}

static void AddPubDescriptorWithText (SeqDescrPtr sdp, Pointer userdata)
{
  ValNodePtr PNTR pub_list;

  if (sdp == NULL || userdata == NULL || sdp->choice != Seq_descr_pub) return;
  pub_list = (ValNodePtr PNTR) userdata;

  ValNodeAddPointer (pub_list, OBJ_SEQDESC, sdp);
}

static ValNodePtr GetTextPubListsForSeqEntry (SeqEntryPtr sep, CharPtr search_text)
{
  FeaturesWithTextData    ftd;
  DescriptorsWithTextData dtd;
  ValNodePtr              obj_list = NULL;

  ftd.seqFeatChoice = SEQFEAT_PUB;
  ftd.featDefChoice = FEATDEF_PUB;
  ftd.search_text = search_text;
  ftd.act_when_string_not_present = FALSE;
  ftd.case_insensitive = TRUE;
  ftd.whole_word = FALSE;
  ftd.no_text = FALSE;
  ftd.callback = AddPubFeatureWithText;
  ftd.userdata = &obj_list;

  dtd.search_text = search_text;
  dtd.act_when_string_not_present = FALSE;
  dtd.case_insensitive = TRUE;
  dtd.whole_word = FALSE;
  dtd.no_text = FALSE;
  dtd.callback = AddPubDescriptorWithText;
  dtd.userdata = &obj_list;

  OperateOnSeqEntryFeaturesWithText (sep, &ftd);
  OperateOnSeqEntryDescriptorsWithText (sep, &dtd);

  return obj_list;
}

static ValNodePtr GetTextSeqEntryListForSeqEntry (SeqEntryPtr sep, Uint2 entityID, CharPtr search_text, ParseFieldPtr dst_field_data)
{
  ValNodePtr value_lists = NULL;
  ValNodePtr obj_list = NULL, tmp_list, vnp;
  ClickableItemPtr cip;
  CharPtr          desc_fmt = "%d sequences contain %s";
  StringConstraintData scd;
  ChoiceConstraintData ccd;
  FilterSetData        fsd;

  if (dst_field_data == NULL) {
    /* do nothing */
  } else if (dst_field_data->parse_field_type == SEARCH_FIELD_PUBLICATION) {
    obj_list = GetTextPubListsForSeqEntry(sep, search_text);
  } else {
    scd.match_text = search_text;
    scd.match_location = 1;
    scd.insensitive = TRUE;
    scd.whole_word = FALSE;
    scd.not_present = FALSE;

    ccd.constraint_type = CHOICE_CONSTRAINT_STRING;
    ccd.qual_choice = 
    ccd.qual_choice_match = NULL;
    ccd.string_constraint = &scd;
    ccd.pseudo_constraint = NULL;
    ccd.free_vn_proc = NULL;
    ccd.copy_vn_proc = NULL;

    MemSet (&fsd, 0, sizeof (FilterSetData));
    fsd.ccp = &ccd;

    tmp_list = GetUniqueValueListForSeqEntry (sep, entityID, dst_field_data, NULL);

    for (vnp = tmp_list; vnp != NULL; vnp = vnp->next) {
      cip = (ClickableItemPtr) vnp->data.ptrvalue;
      if (cip != NULL && StringISearch (cip->description, search_text) != NULL) {
        ValNodeLink (&obj_list, cip->item_list);
        cip->item_list = NULL;
      }
    }
    tmp_list = FreeClickableList (tmp_list);
  }
  if (obj_list != NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->description = StringSave (search_text);
    cip->item_list = obj_list;
    cip->chosen = TRUE;
    ValNodeAddPointer (&value_lists, 0, cip);
    value_lists = ChangeUniqueValueListsToSeqEntryLists (value_lists);
    /* rename description */
    if (value_lists != NULL) {
      cip = value_lists->data.ptrvalue;
      if (cip != NULL) {
        cip->description = MemFree (cip->description);
        cip->description = (CharPtr) MemNew (sizeof(Char) * (StringLen (desc_fmt) + StringLen (search_text) + 15));
        sprintf (cip->description, desc_fmt, ValNodeLen (cip->item_list), search_text);
      }
    }
  }

  return value_lists;
}

typedef enum vectorhitlevel {
  eVectorHitLevelSuspect = 1,
  eVectorHitLevelWeak,
  eVectorHitLevelModerate,
  eVectorHitLevelStrong } EVectorHitLevel;

typedef struct vectorcontaminationoptions {
  Int4    internal_level;
  Int4    end5_level;
  Int4    end3_level;
  ValNodePtr bioseq_list;
} VectorContaminationOptionsData, PNTR VectorContaminationOptionsPtr;

typedef struct vectorcontaminationoptionsdlg {
  DIALOG_MESSAGE_BLOCK

  ButtoN internal;
  ButtoN end5;
  ButtoN end3;
  GrouP  internal_level;
  GrouP  end5_level;
  GrouP  end3_level;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} VectorContaminationOptionsDlgData, PNTR VectorContaminationOptionsDlgPtr;

static void DataToVectorContaminationOptionsDlg (DialoG d, Pointer data)
{
  VectorContaminationOptionsDlgPtr dlg;
  VectorContaminationOptionsPtr    options;

  dlg = (VectorContaminationOptionsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  options = (VectorContaminationOptionsPtr) data;
  if (options == NULL) {
    SetStatus (dlg->internal, FALSE);
    SetStatus (dlg->end5, FALSE);
    SetStatus (dlg->end3, FALSE);
    SetValue (dlg->internal_level, eVectorHitLevelStrong);
    SetValue (dlg->end5_level, eVectorHitLevelStrong);
    SetValue (dlg->end3_level, eVectorHitLevelStrong);
  } else {
    if (options->internal_level == 0) {
      SetStatus (dlg->internal, FALSE);
    } else {
      SetStatus (dlg->internal, TRUE);
      SetValue (dlg->internal_level, options->internal_level);
    }
      
    if (options->end5_level == 0) {
      SetStatus (dlg->end5, FALSE);
    } else {
      SetStatus (dlg->end5, TRUE);
      SetValue (dlg->end5_level, options->end5_level);
    }
    if (options->end3_level == 0) {
      SetStatus (dlg->end3, FALSE);
    } else {
      SetStatus (dlg->end3, TRUE);
      SetValue (dlg->end3_level, options->end3_level);
    }
  }
}


static Pointer VectorContaminationOptionsDlgToData (DialoG d)
{
  VectorContaminationOptionsDlgPtr dlg;
  VectorContaminationOptionsPtr    options = NULL;

  dlg = (VectorContaminationOptionsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  options = (VectorContaminationOptionsPtr) MemNew (sizeof (VectorContaminationOptionsData));
  if (GetStatus (dlg->internal)) {
    options->internal_level = GetValue (dlg->internal_level);
  } else {
    options->internal_level = 0;
  }
  if (GetStatus (dlg->end5)) {
    options->end5_level = GetValue (dlg->end5_level);
  } else {
    options->end5_level = 0;
  }
  if (GetStatus (dlg->end3)) {
    options->end3_level = GetValue (dlg->end3_level);
  } else {
    options->end3_level = 0;
  }

  options->bioseq_list = NULL;
  return (Pointer) options;
}


static void ChangeVectorContaminationButton (ButtoN b)
{
  VectorContaminationOptionsDlgPtr dlg;

  dlg = (VectorContaminationOptionsDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  if (GetStatus (dlg->internal)) {
    Enable (dlg->internal_level);
  } else {
    Disable (dlg->internal_level);
  }
  if (GetStatus (dlg->end5)) {
    Enable (dlg->end5_level);
  } else {
    Disable (dlg->end5_level);
  }
  if (GetStatus (dlg->end3)) {
    Enable (dlg->end3_level);
  } else {
    Disable (dlg->end3_level);
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static void ChangeVectorContaminationGroup (GrouP g)
{
  VectorContaminationOptionsDlgPtr dlg;

  dlg = (VectorContaminationOptionsDlgPtr) GetObjectExtra (g);
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static DialoG VectorContaminationOptionsDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  VectorContaminationOptionsDlgPtr dlg;
  GrouP                     p;
  
  dlg = (VectorContaminationOptionsDlgPtr) MemNew (sizeof (VectorContaminationOptionsDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, 2, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = DataToVectorContaminationOptionsDlg;
  dlg->fromdialog = VectorContaminationOptionsDlgToData;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->internal = CheckBox (p, "Internal", ChangeVectorContaminationButton);
  SetObjectExtra (dlg->internal, dlg, NULL);
  dlg->internal_level = HiddenGroup (p, 4, 0, ChangeVectorContaminationGroup);
  SetObjectExtra (dlg->internal_level, dlg, NULL);
  RadioButton (dlg->internal_level, "Suspect");
  RadioButton (dlg->internal_level, "Weak");
  RadioButton (dlg->internal_level, "Moderate");
  RadioButton (dlg->internal_level, "Strong");
  SetValue (dlg->internal_level, 1);
  Disable (dlg->internal_level);

  dlg->end5 = CheckBox (p, "5' End", ChangeVectorContaminationButton);
  SetObjectExtra (dlg->end5, dlg, NULL);
  dlg->end5_level = HiddenGroup (p, 4, 0, ChangeVectorContaminationGroup);
  SetObjectExtra (dlg->end5_level, dlg, NULL);
  RadioButton (dlg->end5_level, "Suspect");
  RadioButton (dlg->end5_level, "Weak");
  RadioButton (dlg->end5_level, "Moderate");
  RadioButton (dlg->end5_level, "Strong");
  SetValue (dlg->end5_level, 1);
  Disable (dlg->end5_level);

  dlg->end3 = CheckBox (p, "3' End", ChangeVectorContaminationButton);
  SetObjectExtra (dlg->end3, dlg, NULL);
  dlg->end3_level = HiddenGroup (p, 4, 0, ChangeVectorContaminationGroup);
  SetObjectExtra (dlg->end3_level, dlg, NULL);
  RadioButton (dlg->end3_level, "Suspect");
  RadioButton (dlg->end3_level, "Weak");
  RadioButton (dlg->end3_level, "Moderate");
  RadioButton (dlg->end3_level, "Strong");
  SetValue (dlg->end3_level, 1);
  Disable (dlg->end3_level);

  return (DialoG) p;
}


typedef struct chooseseqsbylen {
  ValNodePtr bsp_list;
  Int4 min;
  Int4 max;
} ChooseSeqsByLenData, PNTR ChooseSeqsByLenPtr;


typedef struct chooseseqsbylendlg {
  DIALOG_MESSAGE_BLOCK
  TexT min;
  TexT max;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} ChooseSeqsByLenDlgData, PNTR ChooseSeqsByLenDlgPtr;


static void DataToChooseSeqsByLenDlg (DialoG d, Pointer data)
{
  ChooseSeqsByLenDlgPtr dlg;
  ChooseSeqsByLenPtr    options;
  Char                  buf[255];

  dlg = (ChooseSeqsByLenDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  options = (ChooseSeqsByLenPtr) data;
  if (options == NULL) {
    SetTitle (dlg->min, "0");
    SetTitle (dlg->max, "0");
  } else {
    sprintf (buf, "%d", options->min);
    SetTitle (dlg->min, buf);
    sprintf (buf, "%d", options->max);
    SetTitle (dlg->max, buf);
  }
}


static Pointer ChooseSeqsByLenDlgToData (DialoG d)
{
  ChooseSeqsByLenDlgPtr dlg;
  ChooseSeqsByLenPtr    options = NULL;
  CharPtr               buf;

  dlg = (ChooseSeqsByLenDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  options = (ChooseSeqsByLenPtr) MemNew (sizeof (ChooseSeqsByLenData));
  if (TextHasNoText (dlg->min)) {
    options->min = 0;
  } else {
    buf = SaveStringFromText (dlg->min);
    options->min = atoi (buf);
    buf = MemFree (buf);
  }
  if (TextHasNoText (dlg->max)) {
    options->max = 0;
  } else {
    buf = SaveStringFromText (dlg->max);
    options->max = atoi (buf);
    buf = MemFree (buf);
  }
  return (Pointer) options;
}


static void ChangeChooseSeqsByLenText (TexT t)
{
  ChooseSeqsByLenDlgPtr dlg;

  dlg = (ChooseSeqsByLenDlgPtr) GetObjectExtra (t);
  if (dlg == NULL) {
    return;
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static DialoG ChooseSeqsByLenDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  ChooseSeqsByLenDlgPtr dlg;
  GrouP                     p;
  
  dlg = (ChooseSeqsByLenDlgPtr) MemNew (sizeof (ChooseSeqsByLenDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, 2, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = DataToChooseSeqsByLenDlg;
  dlg->fromdialog = ChooseSeqsByLenDlgToData;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  StaticPrompt (p, "Minimum length", 0, dialogTextHeight, programFont, 'c');
  dlg->min = DialogText (p, "0", 10, ChangeChooseSeqsByLenText);
  SetObjectExtra (dlg->min, dlg, NULL);
  StaticPrompt (p, "Maximum length", 0, dialogTextHeight, programFont, 'c');
  dlg->max = DialogText (p, "50", 10, ChangeChooseSeqsByLenText);
  SetObjectExtra (dlg->max, dlg, NULL);

  return (DialoG) p;
}


typedef enum {
  eSegPageId = 1,
  eSegPageText,
  eSegPageField,
  eSegPageMolInfo,
  eSegPageFeatureType,
  eSegPageDescriptorType,
  eSegPageLength,
  eSegPageVectorContamination,
  eSegPageNumFinalSets,
  eSegPageFile,
  eSegPageConstraint,
  eSegPageStructuredCommentField,
  eNumSegPages
} ESegPage;

typedef struct segbyfield {
  FORM_MESSAGE_BLOCK
  GrouP        seg_choice_grp;

  GrouP        pages[eNumSegPages];
  /* for segregating by field */
  DialoG       field_dlg;
  DialoG       group_list_dlg;

  /* for segregating by molinfo */
  ButtoN       mol_type_btn;
  ButtoN       mol_class_btn;
  ButtoN       mol_topology_btn;

  /* for segregating by feature */
  DialoG       feature_select;

  /* for segregating by descriptor */
  DialoG       descriptor_select;

  /* for segregating by text */
  DialoG       search_field;
  TexT         search_text;

  /* for segregating by ID */
  DialoG       id_constraint_dlg;

  /* for segregating by length */
  DialoG       len_options;

  /* for segregating by vector contamination */
  DialoG       vector_options;

  /* for segregating into a specific number of sets */
  TexT         num_sets;

  /* for segregating by any constraint */
  DialoG       seq_constraint_dlg;

  /* for segregating by structured comment field */
  DialoG       struc_comm_field;

  ButtoN       accept;
  ButtoN       leave_dlg_up;
  BioseqSetPtr target_set;
  ValNodePtr   value_lists;
} SegByFieldData, PNTR SegByFieldPtr;

static void ChangeSegChoice (GrouP g);


static void PullChosenIntoOneGroup (ValNodePtr PNTR value_list, ClickableItemPtr chosen)
{
  ValNodePtr vnp_prev = NULL, vnp_next, vnp_this;
  ClickableItemPtr cip;

  if (value_list == NULL || chosen == NULL) return;

  vnp_this = *value_list;
  while (vnp_this != NULL)
  {
    vnp_next = vnp_this->next;
    cip = (ClickableItemPtr) vnp_this->data.ptrvalue;
    if (cip == NULL)
    {
      vnp_this->next = NULL;
      if (vnp_prev == NULL) 
      {
        *value_list = vnp_next;
      }
      else
      {
        vnp_prev->next = vnp_next;
      }
      vnp_this = ValNodeFree (vnp_this);
    }
    else if (cip->chosen)
    {
      vnp_this->next = NULL;
      if (vnp_prev == NULL) 
      {
        *value_list = vnp_next;
      }
      else
      {
        vnp_prev->next = vnp_next;
      }
      ValNodeLink (&(chosen->item_list), cip->item_list);
      cip->item_list = NULL;
      vnp_this = FreeClickableList (vnp_this);
    }
    else
    {
      PullChosenIntoOneGroup (&(cip->subcategories), chosen);
      vnp_prev = vnp_this;
    }
    vnp_this = vnp_next;
  }
}


static void SegregateByField_Callback (ButtoN b)
{
  SegByFieldPtr sfp;
  Int4          num_chosen = 0;
  ClickableItemPtr cip;
  Int2             val;
  
  sfp = (SegByFieldPtr) GetObjectExtra (b);
  if (sfp == NULL) return;

  num_chosen = CountChosenDiscrepancies (sfp->value_lists, FALSE);
  if (num_chosen < 1) {
    if (ANS_OK == Message (MSG_OKC, "You have not chosen any groups.  Create all?")) {
      ChooseAllDiscrepancies (sfp->value_lists);
    } else {
      return;
    }
  }

  RemovePopsetTitles (GetTopSeqEntryForEntityID (sfp->input_entityID));

  val = GetValue (sfp->seg_choice_grp);
  if (val == eSegPageId || val == eSegPageVectorContamination || val == eSegPageLength || val == eSegPageConstraint) {
    /* put all selected sequences into one group */
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    cip->chosen = TRUE;
    PullChosenIntoOneGroup (&(sfp->value_lists), cip);
    ValNodeAddPointer (&(sfp->value_lists), 0, cip);
  }

  sfp->target_set = MakeGroupsForUniqueValues (sfp->target_set, sfp->value_lists);

  if (GetStatus (sfp->leave_dlg_up)) {
    ChangeSegChoice (sfp->seg_choice_grp);
  } else {
    Remove (sfp->form);
  }
}

static void ChangeSegField (Pointer userdata)
{
  SegByFieldPtr sfp = (SegByFieldPtr) userdata;
  SeqEntryPtr   sep;
  ParseFieldPtr dst_field_data;

  if (sfp == NULL || sfp->target_set == NULL) return;

  WatchCursor();
  Update();
  sep = SeqMgrGetSeqEntryForData (sfp->target_set);

  PointerToDialog (sfp->group_list_dlg, NULL);
  sfp->value_lists = FreeClickableList (sfp->value_lists);

  dst_field_data = DialogToPointer (sfp->field_dlg);
  if (dst_field_data == NULL) {
    Disable (sfp->accept);
    return;
  }

  sfp->value_lists = CreateSeqEntryListsForUniqueValues (sep, sfp->target_set->idx.entityID, dst_field_data, NULL);
  dst_field_data = ParseFieldFree (dst_field_data);
  PointerToDialog (sfp->group_list_dlg, sfp->value_lists);
  Enable (sfp->accept);
  ArrowCursor();
  Update();
}


static void ChangeSegStructuredCommentField (Pointer userdata)
{
  SegByFieldPtr sfp = (SegByFieldPtr) userdata;
  SeqEntryPtr   sep;
  ValNodePtr    vnp;

  if (sfp == NULL || sfp->target_set == NULL) return;

  sep = SeqMgrGetSeqEntryForData (sfp->target_set);

  PointerToDialog (sfp->group_list_dlg, NULL);
  sfp->value_lists = FreeClickableList (sfp->value_lists);

  vnp = DialogToPointer (sfp->struc_comm_field);
  if (vnp == NULL) {
    Disable (sfp->accept);
    return;
  }

  WatchCursor();
  Update();

  sfp->value_lists = GetUniqueValuesForStructuredCommentField (sep, vnp->data.ptrvalue);
  sfp->value_lists = ChangeUniqueValueListsToSeqEntryLists (sfp->value_lists);
  vnp = ValNodeFreeData (vnp);

  PointerToDialog (sfp->group_list_dlg, sfp->value_lists);
  Enable (sfp->accept);
  ArrowCursor();
  Update();
}


static void SegTextChangeNotify (Pointer userdata)
{
  SegByFieldPtr sfp = (SegByFieldPtr) userdata;
  SeqEntryPtr   sep;
  ParseFieldPtr dst_field_data;
  CharPtr       search_text;

  if (sfp == NULL || sfp->target_set == NULL) return;

  WatchCursor();
  Update();
  sep = SeqMgrGetSeqEntryForData (sfp->target_set);

  PointerToDialog (sfp->group_list_dlg, NULL);
  sfp->value_lists = FreeClickableList (sfp->value_lists);

  if (TextHasNoText (sfp->search_text)) {
    Disable (sfp->accept);
    ArrowCursor();
    Update();
    return;
  }
  dst_field_data = DialogToPointer (sfp->search_field);
  if (dst_field_data == NULL) {
    Disable (sfp->accept);
    ArrowCursor();
    Update();
    return;
  }
  search_text = SaveStringFromText (sfp->search_text);

  sfp->value_lists = GetTextSeqEntryListForSeqEntry (sep, sfp->target_set->idx.entityID, search_text, dst_field_data);

  search_text = MemFree (search_text);

  dst_field_data = ParseFieldFree (dst_field_data);
  PointerToDialog (sfp->group_list_dlg, sfp->value_lists);
  if (sfp->value_lists == NULL) {
    Disable(sfp->accept);
  } else {
    Enable (sfp->accept);
  }
  ArrowCursor();
  Update();
}

static void SegTextChange (TexT t)
{
  SegByFieldPtr sfp;

  sfp = (SegByFieldPtr) GetObjectExtra (t);

  SegTextChangeNotify (sfp);

}


static Boolean IsBankitId (SeqIdPtr sip)
{
  DbtagPtr dbtag;

  if (sip == NULL || sip->choice != SEQID_GENERAL) {
    return FALSE;
  }
  if ((dbtag = (DbtagPtr) sip->data.ptrvalue) == NULL) {
    return FALSE;
  }
  if (StringCmp (dbtag->db, "BankIt") == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


NLM_EXTERN ClickableItemPtr ClickableItemForBioseq (BioseqPtr bsp)
{
  ClickableItemPtr cip = NULL;
  SeqIdPtr         sip, sip2;
  Char             id_str[3][100];
  Int4             num_id = 1, len = 0;
  Int4             i;

  if (bsp != NULL && ! ISA_aa (bsp->mol))
  {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    sip = SeqIdFindBest (bsp->id, SEQID_GENBANK);
    SeqIdWrite (sip, id_str[0], PRINTID_REPORT, sizeof (id_str[0]) - 1);
    len = StringLen (id_str[0]) + 1;
    if (!IsBankitId(sip)) {
      sip2 = bsp->id;
      while (sip2 != NULL && !IsBankitId(sip2)) {
        sip2 = sip2->next;
      }
      if (sip2 != NULL) {
        SeqIdWrite (sip2, id_str[1], PRINTID_REPORT, sizeof (id_str[1]) - 1);
        len += StringLen (id_str[1]) + 1;
        num_id++;
      }
    }
    if (!IsNCBIFileID(sip)) {
      sip2 = bsp->id;
      while (sip2 != NULL && !IsNCBIFileID(sip2)) {
        sip2 = sip2->next;
      }
      if (sip2 != NULL) {
        SeqIdWrite (sip2, id_str[1], PRINTID_REPORT, sizeof (id_str[1]) - 1);
        len += StringLen (id_str[1]) + 1;
        num_id++;
      }
    }

    /* space for length */
    len += 20;

    cip->description = (CharPtr) MemNew (sizeof (Char) * len);      
    sprintf (cip->description, "%s", id_str[0]);
    for (i = 1; i < num_id; i++) {
      StringCat (cip->description, "|");
      StringCat (cip->description, id_str[i]);
    }
    sprintf (cip->description + StringLen (cip->description), ":%d", bsp->length);

    ValNodeAddPointer (&(cip->item_list), OBJ_BIOSEQ, bsp);
  }
  return cip;
}


extern void ListAllSequences (BioseqPtr bsp, Pointer userdata)
{
  ClickableItemPtr cip;

  cip = ClickableItemForBioseq (bsp);
  if (cip != NULL) 
  {
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, 0, cip);
  }
}


static void ChooseSegId (Pointer userdata)
{
  SegByFieldPtr sfp = (SegByFieldPtr) userdata;
  SeqEntryPtr   sep;

  if (sfp == NULL || sfp->target_set == NULL) return;

  WatchCursor();
  Update();
  sep = SeqMgrGetSeqEntryForData (sfp->target_set);

  PointerToDialog (sfp->group_list_dlg, NULL);
  sfp->value_lists = FreeClickableList (sfp->value_lists);

  VisitBioseqsInSep (sep, &sfp->value_lists, ListAllSequences);

  PointerToDialog (sfp->group_list_dlg, sfp->value_lists);
  if (sfp->value_lists == NULL) {
    Disable(sfp->accept);
  } else {
    Enable (sfp->accept);
  }
  ArrowCursor();
  Update();
}


static Boolean DoesLevelStringMatchLevel (CharPtr str, Int4 level)
{
  Boolean rval = FALSE;

  switch (level) {
    case eVectorHitLevelSuspect:
      if (StringISearch (str, "suspect") != NULL) {
        rval = TRUE;
      }
      /* note - absence of break is deliberate */
    case eVectorHitLevelWeak:
      if (StringISearch (str, "weak") != NULL) {
        rval = TRUE;
      }
      /* note - absence of break is deliberate */
    case eVectorHitLevelModerate:
      if (StringISearch (str, "moderate") != NULL) {
        rval = TRUE;
      }
      /* note - absence of break is deliberate */
    case eVectorHitLevelStrong:
      if (StringISearch (str, "strong") != NULL) {
        rval = TRUE;
      }
      break;
  }
  return rval;
}


static Boolean DoesVectorMatchOptions (ClickableItemPtr cip, VectorContaminationOptionsPtr options)
{
  Boolean rval = FALSE;

  if (cip == NULL || options == NULL || StringHasNoText (cip->description)) {
    return FALSE;
  }

  if (options->internal_level > 0 
      && StringNICmp (cip->description, "Internal", 8) == 0
      && DoesLevelStringMatchLevel (cip->description, options->internal_level)) {
    rval = TRUE;
  } else if (options->end5_level > 0
      && StringNICmp (cip->description, "5'", 2) == 0
      && DoesLevelStringMatchLevel (cip->description, options->end5_level)) {
    rval = TRUE;
  } else if (options->end3_level > 0
      && StringNICmp (cip->description, "3'", 2) == 0
      && DoesLevelStringMatchLevel (cip->description, options->end3_level)) {
    rval = TRUE;
  }
  return rval;
}


static void GetVectorBioseqList (BioseqPtr bsp, Pointer data)
{
  VectorContaminationOptionsPtr options;
  SeqFeatPtr                    sfp;
  SeqMgrFeatContext             fcontext;
  Boolean                       isvector;
  GBQualPtr                     gbq;
  Int4                          last_right;
  ClickableItemPtr              cip = NULL, cip_match;
  ValNodePtr                    vector_list = NULL, vnp;
  Int4                          desc_len = 0;
  Boolean                       is_match = FALSE;

  if (bsp == NULL || data == NULL) {
    return;
  }

  options = (VectorContaminationOptionsPtr) data;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_IMP, FEATDEF_misc_feature, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_IMP, FEATDEF_misc_feature, &fcontext)) {         
    for (isvector = FALSE, gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
      if (StringCmp (gbq->qual, "standard_name") == 0) {
        if (StringCmp (gbq->val, "Vector Contamination") == 0) {
          isvector = TRUE;
        }
      }
    }
    if (isvector) {
      if (cip != NULL && last_right >= fcontext.left - 1) {
        ValNodeAddPointer (&(cip->item_list), OBJ_SEQFEAT, sfp);
        last_right = fcontext.right;
      } else {        
        /* calculate description for previous list */
        CalculateVectorDescription (cip);
        cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        if (cip != NULL) {
          cip->clickable_item_type = 0;
          cip->callback_func = NULL;
          cip->datafree_func = NULL;
          cip->callback_data = NULL;
          cip->item_list = NULL;
          cip->subcategories = NULL;
          cip->expanded = FALSE;
          cip->level = 0;
          cip->description = NULL;
  
          ValNodeAddPointer (&(cip->item_list), OBJ_SEQFEAT, sfp);

          ValNodeAddPointer (&vector_list, 0, cip);
          last_right = fcontext.right;
        }
      }
    }
  }
  /* calculate description for last item */
  CalculateVectorDescription (cip);

  /* now determine whether we have a hit that matches our requirements */
  for (vnp = vector_list; vnp != NULL; vnp = vnp->next) {
    cip = vnp->data.ptrvalue;
    if (DoesVectorMatchOptions (cip, options)) {
      is_match = TRUE;
    }
    desc_len += StringLen (cip->description) + 2;
  }

  if (is_match) {
    cip_match = MemNew (sizeof (ClickableItemData));
    cip_match->description = (CharPtr) MemNew (sizeof (Char) * desc_len);
    cip_match->description[0] = 0;
    for (vnp = vector_list; vnp != NULL; vnp = vnp->next) {
      cip = vnp->data.ptrvalue;
      StringCat (cip_match->description, cip->description);
      if (vnp->next != NULL) {
        StringCat (cip_match->description, ";");
      }
    }
    ValNodeAddPointer (&(cip_match->item_list), OBJ_BIOSEQ, bsp);
    cip_match->chosen = TRUE;
    ValNodeAddPointer (&(options->bioseq_list), 0, cip_match);
  }
  vector_list = FreeClickableList (vector_list);
}


static void ChooseVectorSeqs (Pointer userdata)
{
  SegByFieldPtr sfp = (SegByFieldPtr) userdata;
  SeqEntryPtr   sep;
  VectorContaminationOptionsPtr options;

  if (sfp == NULL || sfp->target_set == NULL) return;

  WatchCursor();
  Update();
  sep = SeqMgrGetSeqEntryForData (sfp->target_set);

  PointerToDialog (sfp->group_list_dlg, NULL);
  sfp->value_lists = FreeClickableList (sfp->value_lists);

  options = DialogToPointer (sfp->vector_options);

  VisitBioseqsInSep (sep, options, GetVectorBioseqList);

  sfp->value_lists = options->bioseq_list;
  options->bioseq_list = NULL;

  options = MemFree (options);

  PointerToDialog (sfp->group_list_dlg, sfp->value_lists);
  if (sfp->value_lists == NULL) {
    Disable(sfp->accept);
  } else {
    Enable (sfp->accept);
  }
  ArrowCursor();
  Update();
}


static void SeqsByLenCallback (BioseqPtr bsp, Pointer userdata)
{
  ChooseSeqsByLenPtr p = (ChooseSeqsByLenPtr) userdata;
  ClickableItemPtr   cip_match;
  CharPtr            fmt = "%s %d";
  Char               id_str[45];

  if (bsp == NULL || ISA_aa (bsp->mol) || p == NULL) {
    return;
  }

  if (bsp->length >= p->min && bsp->length <= p->max) {
    cip_match = MemNew (sizeof (ClickableItemData));

    /* create description (use ID and length ) */
    SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_str, PRINTID_REPORT, sizeof (id_str) - 1);
    cip_match->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (id_str) + 15));
    sprintf (cip_match->description, fmt, id_str, bsp->length);

    ValNodeAddPointer (&(cip_match->item_list), OBJ_BIOSEQ, bsp);
    cip_match->chosen = TRUE;
    ValNodeAddPointer (&(p->bsp_list), 0, cip_match);
  }
}


static void ChooseSeqsByLen (Pointer userdata)
{
  SegByFieldPtr sfp = (SegByFieldPtr) userdata;
  SeqEntryPtr   sep;
  ChooseSeqsByLenPtr options;

  if (sfp == NULL || sfp->target_set == NULL) return;

  WatchCursor();
  Update();
  sep = SeqMgrGetSeqEntryForData (sfp->target_set);

  PointerToDialog (sfp->group_list_dlg, NULL);
  sfp->value_lists = FreeClickableList (sfp->value_lists);

  options = DialogToPointer (sfp->len_options);

  VisitBioseqsInSep (sep, options, SeqsByLenCallback);

  sfp->value_lists = options->bsp_list;
  options->bsp_list = NULL;

  options = MemFree (options);
  PointerToDialog (sfp->group_list_dlg, sfp->value_lists);
  if (sfp->value_lists == NULL) {
    Disable(sfp->accept);
  } else {
    Enable (sfp->accept);
  }
  ArrowCursor();
  Update();
}


static void ChangeMol (ButtoN b)
{
  SegByFieldPtr sfp;
  SeqEntryPtr   sep;

  sfp = (SegByFieldPtr) GetObjectExtra (b);

  if (sfp == NULL) return;

  WatchCursor ();
  Update();
  sep = SeqMgrGetSeqEntryForData (sfp->target_set);

  PointerToDialog (sfp->group_list_dlg, NULL);
  sfp->value_lists = FreeClickableList (sfp->value_lists);

  if (GetStatus (sfp->mol_class_btn)) {
    VisitBioseqsInSep (sep, &(sfp->value_lists), GetMolClasses);
  }
  if (GetStatus (sfp->mol_type_btn)) {
    VisitDescriptorsInSep (sep, &(sfp->value_lists), GetMolTypes);
  }
  if (GetStatus (sfp->mol_topology_btn)) {
   VisitBioseqsInSep (sep, &(sfp->value_lists), GetMolTopologies);
  }

  sfp->value_lists = ChangeUniqueValueListsToSeqEntryLists (sfp->value_lists);

  PointerToDialog (sfp->group_list_dlg, sfp->value_lists);
  if (ValNodeLen (sfp->value_lists) < 2) {
    Disable(sfp->accept);
  } else {
    Enable (sfp->accept);
  }
  ArrowCursor();
  Update();
}


static void SegFeatureChangeNotify (Pointer userdata)
{
  SegByFieldPtr sfp = (SegByFieldPtr) userdata;
  SeqEntryPtr   sep;
  SegFeaturesData sfd;
  ClickableItemPtr cip;
  ValNodePtr       vnp;
  CharPtr          desc_fmt = "%d %ss";

  if (sfp == NULL) return;

  vnp = DialogToPointer (sfp->feature_select);
  if (vnp == NULL) {
    PointerToDialog (sfp->group_list_dlg, NULL);
    sfp->value_lists = FreeClickableList (sfp->value_lists);
    Disable (sfp->accept);
    return;
  }

  WatchCursor ();
  Update();
  sep = SeqMgrGetSeqEntryForData (sfp->target_set);

  PointerToDialog (sfp->group_list_dlg, NULL);
  sfp->value_lists = FreeClickableList (sfp->value_lists);

  sfd.featdef = vnp->choice;
  sfd.feature_list = NULL;
  VisitFeaturesInSep (sep, &sfd, FindFeaturesForSeg);

  if (sfd.feature_list != NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->description = (CharPtr) MemNew (sizeof(Char) * (StringLen (desc_fmt) + StringLen (vnp->data.ptrvalue) + 15));
    sprintf (cip->description, desc_fmt, ValNodeLen (sfd.feature_list), vnp->data.ptrvalue);
    cip->item_list = sfd.feature_list;
    cip->chosen = TRUE;
    ValNodeAddPointer (&(sfp->value_lists), 0, cip);
    sfp->value_lists = ChangeUniqueValueListsToSeqEntryLists (sfp->value_lists);
  }
  vnp = ValNodeFree (vnp);

  PointerToDialog (sfp->group_list_dlg, sfp->value_lists);

  if (sfp->value_lists == NULL) {
    Disable (sfp->accept);
  } else {
    Enable (sfp->accept);
  }
  ArrowCursor();
  Update();
}

static void SegDescriptorChangeNotify (Pointer data)
{
  SegByFieldPtr      sfp;
  ValNodePtr         vnp;
  SegDescriptorsData sdd;
  ClickableItemPtr   cip;
  CharPtr            desc_fmt = "%d %ss";
  SeqEntryPtr        sep;

  sfp = (SegByFieldPtr) data;
  if (sfp == NULL) return;

  vnp = DialogToPointer (sfp->descriptor_select);
  if (vnp == NULL) {
    PointerToDialog (sfp->group_list_dlg, NULL);
    sfp->value_lists = FreeClickableList (sfp->value_lists);
    Disable (sfp->accept);
    return;
  }

  WatchCursor ();
  Update();
  sep = SeqMgrGetSeqEntryForData (sfp->target_set);

  PointerToDialog (sfp->group_list_dlg, NULL);
  sfp->value_lists = FreeClickableList (sfp->value_lists);

  sdd.desc_type = vnp->choice;
  sdd.descriptor_list = NULL;
  VisitDescriptorsInSep (sep, &sdd, FindDescriptorsForSeg);

  if (sdd.descriptor_list != NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->description = (CharPtr) MemNew (sizeof(Char) * (StringLen (desc_fmt) + StringLen (vnp->data.ptrvalue) + 15));
    sprintf (cip->description, desc_fmt, ValNodeLen (sdd.descriptor_list), vnp->data.ptrvalue);
    cip->item_list = sdd.descriptor_list;
    cip->chosen = TRUE;
    ValNodeAddPointer (&(sfp->value_lists), 0, cip);
    sfp->value_lists = ChangeUniqueValueListsToSeqEntryLists (sfp->value_lists);
  }
  vnp = ValNodeFree (vnp);

  PointerToDialog (sfp->group_list_dlg, sfp->value_lists);

  if (sfp->value_lists == NULL) {
    Disable (sfp->accept);
  } else {
    Enable (sfp->accept);
  }
  ArrowCursor();
  Update();

}

static void ChooseSeqsForNumSets (Pointer data)
{
  SegByFieldPtr      sfp;
  SeqEntryPtr        sep;
  Int4               num_sets;
  CharPtr            str;

  sfp = (SegByFieldPtr) data;
  if (sfp == NULL) return;

  str = SaveStringFromText (sfp->num_sets);
  if (str == NULL) {
    Disable (sfp->accept);
    return;
  }

  num_sets = atoi (str);
  str = MemFree (str);
  if (num_sets < 1) {
    Disable (sfp->accept);
    return;
  }

  WatchCursor ();
  Update();
  sep = SeqMgrGetSeqEntryForData (sfp->target_set);

  PointerToDialog (sfp->group_list_dlg, NULL);
  sfp->value_lists = FreeClickableList (sfp->value_lists);

  sfp->value_lists = PrepareSequenceListForSegregateByNumberOfSets (num_sets, sep);

  PointerToDialog (sfp->group_list_dlg, sfp->value_lists);

  if (sfp->value_lists == NULL) {
    Disable (sfp->accept);
  } else {
    Enable (sfp->accept);
  }
  ArrowCursor();
  Update();
}

static void SegNumSetsChange (TexT t)
{
  SegByFieldPtr sfp;

  sfp = (SegByFieldPtr) GetObjectExtra (t);

  ChooseSeqsForNumSets (sfp);
}


static void GroupSeqByFileCallback (BioseqPtr bsp, Pointer data)
{
  ValNodePtr PNTR list;
  ValNodePtr vnp;
  SeqIdPtr   sip;
  DbtagPtr   dbtag;
  Char       id_buf[15];
  CharPtr    id = NULL, last_slash;
  ClickableItemPtr cip;
  Boolean          found = FALSE;

  if (bsp == NULL || ISA_aa (bsp->mol) || (list = (ValNodePtr PNTR) data) == NULL) {
    return;
  }

  sip = bsp->id; 
  while (sip != NULL && !IsNCBIFileID(sip)) {
    sip = sip->next;
  }

  if (sip != NULL && (dbtag = (DbtagPtr) sip->data.ptrvalue) != NULL) {
    if (dbtag->tag->id > 0) {
      sprintf (id_buf, "%d", dbtag->tag->id);
      id = id_buf;
    } else {
      id = dbtag->tag->str;
    }
    id = StringSave (id);
    if ((last_slash = StringRChr (id, '/')) != NULL) {
      *last_slash = 0;
    }
    for (vnp = *list; vnp != NULL && !found; vnp = vnp->next) {
      cip = (ClickableItemPtr) vnp->data.ptrvalue;
      if (cip != NULL && StringCmp (id, cip->description) == 0) {
        ValNodeAddPointer (&cip->item_list, OBJ_BIOSEQ, bsp);
        found = TRUE;
        id = MemFree (id);
      }
    }
    if (!found) {
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      cip->description = id;
      ValNodeAddPointer (&cip->item_list, OBJ_BIOSEQ, bsp);
      ValNodeAddPointer (list, 0, cip);
    }
  }
}


static void ChooseFile (SegByFieldPtr sfp)
{
  SeqEntryPtr        sep;

  if (sfp == NULL) return;

  WatchCursor ();
  Update();
  sep = SeqMgrGetSeqEntryForData (sfp->target_set);

  PointerToDialog (sfp->group_list_dlg, NULL);
  sfp->value_lists = FreeClickableList (sfp->value_lists);

  VisitBioseqsInSep (sep, &(sfp->value_lists), GroupSeqByFileCallback);

  PointerToDialog (sfp->group_list_dlg, sfp->value_lists);

  if (sfp->value_lists == NULL) {
    Disable (sfp->accept);
  } else {
    Enable (sfp->accept);
  }
  ArrowCursor();
  Update();
}


static void ChangeSegChoice (GrouP g)
{
  SegByFieldPtr sfp;
  Int4          seg_choice, i;

  sfp = (SegByFieldPtr) GetObjectExtra (g);
  if (sfp == NULL) return;

  seg_choice = GetValue (sfp->seg_choice_grp);

    for (i = 0; i < eNumSegPages; i++) {
    if (i != seg_choice - 1) {
        SafeHide (sfp->pages[i]);
      }
    }
    Update ();

    for (i = 0; i < eNumSegPages; i++) {
    if (i == seg_choice - 1) {
        SafeShow (sfp->pages[i]);
      }
    }
    Update ();

    switch (seg_choice) {
      case eSegPageField:
        ChangeSegField (sfp);
          break;
      case eSegPageMolInfo:
        ChangeMol(sfp->mol_class_btn);
        break;
      case eSegPageFeatureType:
        SegFeatureChangeNotify (sfp);
        break;
    case eSegPageDescriptorType:
        SegDescriptorChangeNotify (sfp);
        break;
    case eSegPageText:
        SegTextChangeNotify (sfp);
        break;
    case eSegPageId:
    case eSegPageConstraint:
        ChooseSegId (sfp);
        break;
    case eSegPageVectorContamination:
        ChooseVectorSeqs (sfp);
        break;
    case eSegPageLength:
        ChooseSeqsByLen (sfp);
        break;
    case eSegPageNumFinalSets:
        ChooseSeqsForNumSets (sfp);
        break;
    case eSegPageFile:
        ChooseFile(sfp);
        break;
    case eSegPageStructuredCommentField:
        ChangeSegStructuredCommentField (sfp);
        break;
  }
  if (seg_choice == eSegPageId || seg_choice == eSegPageConstraint 
      || seg_choice == eSegPageVectorContamination || seg_choice == eSegPageLength) {
    SetClickableListDialogTitles (sfp->group_list_dlg, "Sequences for New Set", "",
                                                       "Use checkbox to mark sequences for new set",
                                                       "Single click to navigate to sequence in record");
    SetClickableListDisplayChosen (sfp->group_list_dlg, TRUE);
  } else {
    SetClickableListDialogTitles (sfp->group_list_dlg, "New Sets", "Sequences in Highlighted Set",
                                                       "Use checkbox to mark sets to create",
                                                       "");
    SetClickableListDisplayChosen (sfp->group_list_dlg, FALSE);
  }
}


extern void ChooseCategories (ValNodePtr value_list, Boolean do_choose)
{
  ClickableItemPtr cip;

  while (value_list != NULL) {
    cip = (ClickableItemPtr) value_list->data.ptrvalue;
      if (cip != NULL) {
        cip->chosen = do_choose;
        ChooseCategories (cip->subcategories, FALSE);
      }
      value_list = value_list->next;
  }
}


static void SelectSegCategories (ButtoN b)
{
  SegByFieldPtr sfp;

  sfp = (SegByFieldPtr) GetObjectExtra (b);
  if (sfp != NULL) {
    ChooseCategories (sfp->value_lists, TRUE);
    PointerToDialog (sfp->group_list_dlg, sfp->value_lists);
  }
}


static void UnselectSegCategories (ButtoN b)
{
  SegByFieldPtr sfp;

  sfp = (SegByFieldPtr) GetObjectExtra (b);
  if (sfp != NULL) {
    ChooseCategories (sfp->value_lists, FALSE);
    PointerToDialog (sfp->group_list_dlg, sfp->value_lists);
  }
}


static void ResortMarkedSegCategories (ButtoN b)
{
  SegByFieldPtr sfp;

  sfp = (SegByFieldPtr) GetObjectExtra (b);
  if (sfp != NULL) {
    sfp->value_lists = ValNodeSort (sfp->value_lists, SortVnpByClickableItemChosen);
    PointerToDialog (sfp->group_list_dlg, sfp->value_lists);
  }
}


static void SelectSequenceIDsForSegregate (ButtoN b)
{
  SegByFieldPtr       sfp;

  sfp = (SegByFieldPtr) GetObjectExtra (b);
  if (sfp == NULL) return;

  ChooseCategoriesForSegregateIdDialog(sfp->id_constraint_dlg, sfp->value_lists);
  PointerToDialog (sfp->group_list_dlg, sfp->value_lists);
}


static void SelectSequenceIDsByConstraintForSegregate (ButtoN b)
{
  SegByFieldPtr       sfp;

  sfp = (SegByFieldPtr) GetObjectExtra (b);
  if (sfp == NULL) return;

  ChooseCategoriesForSegregateConstraintDialog(sfp->seq_constraint_dlg, sfp->value_lists);
  PointerToDialog (sfp->group_list_dlg, sfp->value_lists);
}


static void UnselectSequenceIdsByConstraintForSegregate (ButtoN b)
{
  SegByFieldPtr       sfp;

  sfp = (SegByFieldPtr) GetObjectExtra (b);
  if (sfp == NULL) return;

  UnchooseCategoriesForSegregateConstraintDialog(sfp->seq_constraint_dlg, sfp->value_lists);
  PointerToDialog (sfp->group_list_dlg, sfp->value_lists);
}


static void CleanupSegregateByFieldForm (GraphiC g, VoidPtr data)

{
  SegByFieldPtr      sfp;

  sfp = (SegByFieldPtr) data;
  if (sfp != NULL) {
    sfp->value_lists = FreeClickableList (sfp->value_lists);
  }
  StdCleanupFormProc (g, data);
}



static Int2 NewSegregateBioseqSet (BioseqSetPtr bssp, Nlm_BtnActnProc actn, ESegPage default_type)
{
  GrouP              c;
  SegByFieldPtr      cfp;
  GrouP              h, k, page_grp, g1, g2, g3;
  PrompT             p, p2;
  StdEditorProcsPtr  sepp;
  WindoW             w;
  ButtoN             b;


  /* Create a new window, and a struct */
  /* to pass around the data in.       */

  cfp = (SegByFieldPtr) MemNew (sizeof (SegByFieldData));
  if (cfp == NULL)
    return OM_MSG_RET_ERROR;
  cfp->target_set = bssp;

  w = FixedWindow (-50, -33, -10, -10, "Segregate By Field",
           StdCloseWindowProc);
  SetObjectExtra (w, cfp, CleanupSegregateByFieldForm);
  cfp->form = (ForM) w;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    cfp->appmessage = sepp->handleMessages;
  }

  cfp->input_entityID = bssp->idx.entityID;
  cfp->input_itemID = bssp->idx.itemID;
  cfp->input_itemtype = OBJ_BIOSEQSET;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    cfp->appmessage = sepp->handleMessages;
  }

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  cfp->group_list_dlg = CreateClickableListDialogEx (h, "New Sets", "Sequences for New Set",
                                                      "Use checkbox to mark sets to create",
                                                      "",
                                                      ScrollToDiscrepancyItem, EditDiscrepancyItem, NULL,
                                                      GetDiscrepancyItemText,
                                                      stdCharWidth * 30,
                                                      stdCharWidth * 30 + 5,
                                                      TRUE, FALSE);

  g3 = NormalGroup (h, 2, 0, "", programFont, NULL);
  SetGroupSpacing (g3, 10, 10);
  StaticPrompt (g3, "Segregate Based on", 0, dialogTextHeight, programFont, 'c');
  cfp->seg_choice_grp = HiddenGroup (g3, eNumSegPages, 0, ChangeSegChoice);
  SetObjectExtra (cfp->seg_choice_grp, cfp, NULL);
  SetGroupSpacing (cfp->seg_choice_grp, 10, 10);

  RadioButton (cfp->seg_choice_grp, "ID");
  RadioButton (cfp->seg_choice_grp, "Text");
  RadioButton (cfp->seg_choice_grp, "Field");
  RadioButton (cfp->seg_choice_grp, "MolInfo");
  RadioButton (cfp->seg_choice_grp, "Feature Type");
  RadioButton (cfp->seg_choice_grp, "Descriptor Type");
  RadioButton (cfp->seg_choice_grp, "Nucleotide Sequence Length");
  RadioButton (cfp->seg_choice_grp, "Vector Contamination");
  RadioButton (cfp->seg_choice_grp, "Number of Sets");
  RadioButton (cfp->seg_choice_grp, "File Name");
  RadioButton (cfp->seg_choice_grp, "Constraint");
  RadioButton (cfp->seg_choice_grp, "Structured Comment Field");

  page_grp = HiddenGroup (h, 0, 0, NULL);
  /* identical fields */
  cfp->pages[eSegPageField - 1] = HiddenGroup (page_grp, 2, 0, NULL);
  SetGroupSpacing (cfp->pages[eSegPageField - 1], 10, 10);
  p = StaticPrompt (cfp->pages[eSegPageField - 1], "Segregate sequences with identical", 0, dialogTextHeight,
                    programFont, 'c');  
  cfp->field_dlg = ParseFieldDestDialog (cfp->pages[eSegPageField - 1], ChangeSegField, cfp);

  /* molinfo */
  cfp->pages[eSegPageMolInfo - 1] = HiddenGroup (page_grp, 0, 4, NULL);
  SetGroupSpacing (cfp->pages[eSegPageMolInfo - 1], 10, 10);
  cfp->mol_class_btn = CheckBox (cfp->pages[eSegPageMolInfo - 1], "Molecule Class", ChangeMol);
  SetObjectExtra (cfp->mol_class_btn, cfp, NULL);
  SetStatus (cfp->mol_class_btn, TRUE);
  cfp->mol_type_btn = CheckBox (cfp->pages[eSegPageMolInfo - 1], "Molecule Type", ChangeMol);
  SetObjectExtra (cfp->mol_type_btn, cfp, NULL);
  SetStatus (cfp->mol_type_btn, TRUE);
  cfp->mol_topology_btn = CheckBox (cfp->pages[eSegPageMolInfo - 1], "Molecule Topology", ChangeMol);
  SetObjectExtra (cfp->mol_topology_btn, cfp, NULL);
  SetStatus (cfp->mol_topology_btn, TRUE);

  /* features */
  cfp->pages[eSegPageFeatureType - 1] = HiddenGroup (page_grp, 2, 0, NULL);
  SetGroupSpacing (cfp->pages[eSegPageFeatureType - 1], 10, 10);
  StaticPrompt (cfp->pages[eSegPageFeatureType - 1], "Segregate by Feature", 0, dialogTextHeight,
                programFont, 'c');
  cfp->feature_select =  FeatureSelectionDialogEx (cfp->pages[eSegPageFeatureType - 1], FALSE, NULL,
                                                   SegFeatureChangeNotify, 
                                                   cfp);

  /* descriptors */
  cfp->pages[eSegPageDescriptorType - 1] = HiddenGroup (page_grp, 2, 0, NULL);
  SetGroupSpacing (cfp->pages[eSegPageDescriptorType - 1], 10, 10);
  StaticPrompt (cfp->pages[eSegPageDescriptorType - 1], "Segregate by Descriptor", 0, dialogTextHeight,
                programFont, 'c');
  cfp->descriptor_select = DescriptorSelectionDialog (cfp->pages[eSegPageDescriptorType - 1], FALSE, 
                                                      SegDescriptorChangeNotify, cfp);

  /* text */
  cfp->pages[eSegPageText - 1] = HiddenGroup (page_grp, 4, 0, NULL);
  SetGroupSpacing (cfp->pages[eSegPageText - 1], 10, 10);
  StaticPrompt (cfp->pages[eSegPageText - 1], "Create set with sequences where", 0, dialogTextHeight, programFont, 'l');
  cfp->search_text = DialogText (cfp->pages[eSegPageText - 1], "", 10, SegTextChange);
  SetObjectExtra (cfp->search_text, cfp, NULL);
  StaticPrompt (cfp->pages[eSegPageText - 1], "appears in", 0, dialogTextHeight, programFont, 'l');
  cfp->search_field = SearchFieldDialog (cfp->pages[eSegPageText - 1], SegTextChangeNotify, cfp);

  /* id */
  cfp->pages[eSegPageId - 1] = HiddenGroup (page_grp, -1, 0, NULL);
  SetGroupSpacing (cfp->pages[eSegPageId - 1], 10, 10);
  p2 = StaticPrompt (cfp->pages[eSegPageId - 1], "Create set with marked sequences", 0, dialogTextHeight, programFont, 'l');
  g1 = HiddenGroup (cfp->pages[eSegPageId - 1], 3, 0, NULL);
  g2 = HiddenGroup (cfp->pages[eSegPageId - 1], 2, 0, NULL);
  cfp->id_constraint_dlg = SegregateIdDialog (g2);
  b = PushButton (g2, "Mark", SelectSequenceIDsForSegregate);
  SetObjectExtra (b, cfp, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) p2, (HANDLE) g1, (HANDLE) g2, NULL);

  /* length */
  cfp->pages[eSegPageLength - 1] = HiddenGroup (page_grp, -1, 0, NULL);
  SetGroupSpacing (cfp->pages[eSegPageLength - 1], 10, 10);
  cfp->len_options = ChooseSeqsByLenDialog (cfp->pages[eSegPageLength - 1], ChooseSeqsByLen, cfp);

  /* vector contamination */
  cfp->pages[eSegPageVectorContamination - 1] = HiddenGroup (page_grp, -1, 0, NULL);
  SetGroupSpacing (cfp->pages[eSegPageVectorContamination - 1], 10, 10);
  cfp->vector_options = VectorContaminationOptionsDialog (cfp->pages[eSegPageVectorContamination - 1], ChooseVectorSeqs, cfp);

  /* number of sets */
  cfp->pages[eSegPageNumFinalSets - 1] = HiddenGroup (page_grp, 2, 0, NULL);
  SetGroupSpacing (cfp->pages[eSegPageNumFinalSets - 1], 10, 10);
  StaticPrompt (cfp->pages[eSegPageNumFinalSets - 1], "Number of sets to create", 0, dialogTextHeight, programFont, 'l');
  cfp->num_sets = DialogText (cfp->pages[eSegPageNumFinalSets - 1], "", 10, SegNumSetsChange);
  SetObjectExtra (cfp->num_sets, cfp, NULL);

  /* don't need controls for File */
  cfp->pages[eSegPageFile - 1] = NULL;

  /* constraint */
  cfp->pages[eSegPageConstraint- 1] = HiddenGroup (page_grp, -1, 0, NULL);
  SetGroupSpacing (cfp->pages[eSegPageConstraint - 1], 10, 10);
  p2 = StaticPrompt (cfp->pages[eSegPageConstraint - 1], "Create set with marked sequences", 0, dialogTextHeight, programFont, 'l');
  g1 = HiddenGroup (cfp->pages[eSegPageConstraint - 1], 2, 0, NULL);
  cfp->seq_constraint_dlg = SegregateConstraintDialog (g1);
  g2 = HiddenGroup (g1, 0, 2, NULL);
  SetGroupSpacing (g2, 10, 10);
  b = PushButton (g2, "Mark", SelectSequenceIDsByConstraintForSegregate);
  SetObjectExtra (b, cfp, NULL);
  b = PushButton (g2, "Unmark", UnselectSequenceIdsByConstraintForSegregate);
  SetObjectExtra (b, cfp, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) p2, (HANDLE) g1, NULL);

  cfp->pages[eSegPageStructuredCommentField - 1] = HiddenGroup (page_grp, 2, 0, NULL);
  SetGroupSpacing (cfp->pages[eSegPageStructuredCommentField - 1], 10, 10);
  StaticPrompt (cfp->pages[eSegPageStructuredCommentField - 1], "Choose Structured Comment Field Name", 0, dialogTextHeight, programFont, 'l');
  cfp->struc_comm_field = ValNodeSelectionDialog (cfp->pages[eSegPageStructuredCommentField - 1],
                                                  GetStructuredCommentFieldNames(SeqMgrGetSeqEntryForData(bssp)),
                                                  8, ValNodeStringName, ValNodeSimpleDataFree,
                                                  ValNodeStringCopy, ValNodeStringMatch, "fieldname", 
                                                  ChangeSegStructuredCommentField, cfp, FALSE);

  AlignObjects (ALIGN_CENTER, (HANDLE) cfp->pages[0],
                                (HANDLE) cfp->pages[1], 
                                (HANDLE) cfp->pages[2], 
                              (HANDLE) cfp->pages[3],
                              (HANDLE) cfp->pages[4],
                              (HANDLE) cfp->pages[5],
                              (HANDLE) cfp->pages[6],
                              (HANDLE) cfp->pages[7],
                              (HANDLE) cfp->pages[8],
                              (HANDLE) cfp->pages[9],
                              (HANDLE) cfp->pages[10],
                              (HANDLE) cfp->pages[11],
                              (HANDLE) cfp->pages[12],
                              NULL);

  SetValue (cfp->seg_choice_grp, default_type);

  /* add select all/unselect all buttons */
  k = HiddenGroup (h, 3, 0, NULL);
  b = PushButton (k, "Mark All", SelectSegCategories);
  SetObjectExtra (b, cfp, NULL);
  b = PushButton (k, "Unmark All", UnselectSegCategories);
  SetObjectExtra (b, cfp, NULL);  
  b = PushButton (k, "Resort Marked", ResortMarkedSegCategories);
  SetObjectExtra (b, cfp, NULL);

  /* Add Accept and Cancel buttons */

  c = HiddenGroup (h, 3, 0, NULL);
  cfp->accept = DefaultButton (c, "Accept", actn);
  SetObjectExtra (cfp->accept, cfp, NULL);
  Disable (cfp->accept);
  PushButton (c, "Cancel", StdCancelButtonProc);
  cfp->leave_dlg_up = CheckBox (c, "Leave Dialog Up", NULL);

  /* Line things up nicely */

  AlignObjects (ALIGN_CENTER, (HANDLE) g3,
                              (HANDLE) cfp->group_list_dlg,
                              (HANDLE) page_grp,
                              (HANDLE) k,
                              (HANDLE) c, NULL);

  /* initialize the display */
  ChangeSegChoice (cfp->seg_choice_grp);

  /* Display the window now */

  RealizeWindow (w);
  Show (w);
  Select (w);
  Select (cfp->accept);
  Update ();
  return OM_MSG_RET_OK;
}

static Int2 LIBCALLBACK SegregateSetsByFieldEx (Pointer data, Nlm_BtnActnProc actn, ESegPage default_type)
{
  OMProcControlPtr   ompcp;

  /* Check parameters and get a pointer to the current data */

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL
      || ompcp->input_itemtype != OBJ_BIOSEQSET 
      || ompcp->input_data == NULL) {
    Message (MSG_ERROR, "You must select a set to segregate!");
    return OM_MSG_RET_ERROR;
  } 
  return NewSegregateBioseqSet (ompcp->input_data, actn, default_type);
}


extern Int2 LIBCALLBACK SegregateSetsByField (Pointer data)
{
  return SegregateSetsByFieldEx (data, SegregateByField_Callback, eSegPageText);
}


NLM_EXTERN void SegregateBioseqSetEntityID (Uint2 entityID)
{
  SeqEntryPtr sep;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL || !IS_Bioseq_set (sep)) {
    Message (MSG_ERROR, "This record does not have a top-levelset!");
  } else {
    NewSegregateBioseqSet (FindTopLevelSetForDesktopFunction((BioseqSetPtr) sep->data.ptrvalue), SegregateByField_Callback, eSegPageText);
  }
}


extern void NewSegregateBioseqSetMenuItem (IteM i)
{
  BaseFormPtr bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  SegregateBioseqSetEntityID (bfp->input_entityID);
}


typedef struct unsequester{
  SeqEntryPtr  sep;
  BioseqSetPtr orig_parent;
  Int4         set_pos;
} UnsequesterData, PNTR UnsequesterPtr;


static UnsequesterPtr UnsequesterNew (SeqEntryPtr sep, BioseqSetPtr parent, Int4 set_pos)
{
  UnsequesterPtr u;

  u = (UnsequesterPtr) MemNew (sizeof (UnsequesterData));
  u->sep = sep; 
  u->orig_parent = parent;
  u->set_pos = set_pos;
  return u;
}


static UnsequesterPtr GetUnsequesterForSep (SeqEntryPtr sep, ValNodePtr list)
{
  UnsequesterPtr u;

  while (list != NULL) {
    if ((u = list->data.ptrvalue) != NULL && u->sep == sep) {
      return u;
    } else {
      list = list->next;
    }
  }
  return NULL;
}


static ValNodePtr unsequester_list = NULL;

static BioseqSetPtr sequester_set = NULL;
static const CharPtr SEQUESTER_BACKUP_FILE = "sequestr.bkp";
static Boolean doReportAfterUnsequester = 0;

NLM_EXTERN void SetDoReportAfterUnsequester (Int4 report_type)
{
  doReportAfterUnsequester = report_type;
}


static const CharPtr placeholder_id_string = "ABCDEFGHIJKLMNOPQRSTUVWXYZMagicCookiePlaceHolderABCDEFGHIJKLMNOPQRSTUVWXYZ";

static BioseqPtr PlaceholderBioseq (void)
{
  BioseqPtr bsp;
  ObjectIdPtr oip;

  oip = ObjectIdNew ();
  oip->str = StringSave (placeholder_id_string);
  bsp = BioseqNew ();
  bsp->id = ValNodeNew (NULL);
  bsp->id->choice = SEQID_LOCAL;
  bsp->id->data.ptrvalue = oip;
  return bsp;
}


static Boolean IsPlaceholderBioseq (BioseqPtr bsp)
{
  ObjectIdPtr oip;

  if (bsp == NULL || bsp->id == NULL 
      || bsp->id->choice != SEQID_LOCAL 
      || (oip = (ObjectIdPtr) bsp->id->data.ptrvalue) == NULL
      || StringCmp (oip->str, placeholder_id_string) != 0) {
    return FALSE;
  } else {
    return TRUE;
  }
}


/* sequesters is different from segregate - we aren't creating a permanent new set with
 * the combined descriptors of the previous set; we're creating a temporary holding set.
 */
static void SequesterCategorySeqEntries (BioseqSetPtr newset, ClickableItemPtr category)
{
  ValNodePtr vnp_item;
  SeqEntryPtr sep, last_sep, prev_sep, remove_sep, set_sep, tmp_sep;
  BioseqSetPtr bssp, orig_parent, prev_new_parent = NULL, prev_orig_parent = NULL;
  BioseqPtr bsp;
  Int4      set_pos;

  if (newset == NULL || category == NULL || (category->item_list == NULL && category->subcategories == NULL)) return;

  if (category->chosen) {
    last_sep = newset->seq_set;
    while (last_sep != NULL && last_sep->next != NULL) {
      last_sep = last_sep->next;
    }

    for (vnp_item = category->item_list; vnp_item != NULL; vnp_item = vnp_item->next) {
      sep = GetBestSeqEntryForItem (vnp_item);
      if (sep == NULL || sep->data.ptrvalue == NULL) continue;
      orig_parent = NULL;
      if (IS_Bioseq (sep)) {
        bsp = sep->data.ptrvalue;
        if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
          orig_parent = bsp->idx.parentptr;
          bsp->idx.parentptr = NULL;
        }
      } else if (IS_Bioseq_set (sep)) {
        bssp = sep->data.ptrvalue;
        if (bssp->idx.parenttype == OBJ_BIOSEQSET) {
          orig_parent = bssp->idx.parentptr;
          bssp->idx.parentptr = NULL;
        }
      } else {
        continue;
      }

      if (orig_parent == NULL) {
        continue;
      }
         
      /* remove this seq-entry from the original parent */
      prev_sep = NULL;
      set_pos = 0;
      for (remove_sep = orig_parent->seq_set;
            remove_sep != NULL && remove_sep != sep;
            remove_sep = remove_sep->next) {
        prev_sep = remove_sep;
        set_pos++;
      }
      if (remove_sep == sep) {
        if (prev_sep == NULL) {
          orig_parent->seq_set = orig_parent->seq_set->next;
          if (orig_parent->seq_set == NULL) {
            orig_parent->seq_set = SeqEntryNew ();
            orig_parent->seq_set->choice = 1;
            orig_parent->seq_set->data.ptrvalue = PlaceholderBioseq();     
            SeqMgrLinkSeqEntry (orig_parent->seq_set, OBJ_BIOSEQSET, orig_parent);
          }
        } else {
          prev_sep->next = sep->next;
        }
        sep->next = NULL;
      }

      /* if this seq-entry is from the same parent as the previous one, just add this to
       * our current placeholder.
       */
      if (prev_orig_parent != orig_parent) {
        prev_new_parent = BioseqSetNew ();
        prev_new_parent->_class = orig_parent->_class;
        /* add descriptors from the orig_parent to the new parent */
        AddNewUniqueDescriptors (&(prev_new_parent->descr), orig_parent->descr);

        /* add annotations from the orig_parent to the new parent */
        AddNewUniqueAnnotations (&(prev_new_parent->annot), orig_parent->annot);

        /* add new parent to set */
        set_sep = SeqEntryNew ();
        set_sep->choice = 2;
        set_sep->data.ptrvalue = prev_new_parent;

        if (newset->seq_set == NULL) {
          newset->seq_set = set_sep;
        } else {
          for (tmp_sep = newset->seq_set; tmp_sep->next != NULL; tmp_sep = tmp_sep->next) {
          }
          tmp_sep->next = set_sep;
        }
        SeqMgrLinkSeqEntry (set_sep, OBJ_BIOSEQSET, newset);
        last_sep = NULL;
      }
      /* add seqentry to new parent */
      if (last_sep == NULL) {
        prev_new_parent->seq_set = sep;
      } else {
        last_sep->next = sep;
      }
      last_sep = sep;

      SeqMgrLinkSeqEntry (sep, OBJ_BIOSEQSET, prev_new_parent);

      /* create unsequester entry */
      ValNodeAddPointer (&unsequester_list, 0, UnsequesterNew (sep, orig_parent, set_pos));

      prev_orig_parent = orig_parent;
    }
  } else {
    for (vnp_item = category->subcategories; vnp_item != NULL; vnp_item = vnp_item->next) {
      SequesterCategorySeqEntries (newset, vnp_item->data.ptrvalue);
    }
  }
}


static void RestoreSequesteredSets (ForM f)
{
  Uint2 new_entityID;
  BaseFormPtr bfp;
  BioseqViewFormPtr vfp;

  bfp = GetBaseFormForEntityID (sequester_set->idx.entityID);
  if (bfp == NULL) {
    return;
  }
  SequinSeqViewFormMessage (f, VIB_MSG_CLOSE);

  new_entityID = RestoreEntityIDFromFile (SEQUESTER_BACKUP_FILE, bfp->input_entityID);
  if (new_entityID != 0) {
    bfp->input_entityID = new_entityID;
    ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  }
  Show (bfp->form);
  vfp = (BioseqViewFormPtr) bfp;
  Show (vfp->toolForm);

  sequester_set = NULL;
}


static Boolean IsSeqDescInSet (SeqDescrPtr sdp, SeqDescrPtr set)
{
  while (set != NULL) {
    if (AsnIoMemComp (sdp, set, (AsnWriteFunc) SeqDescAsnWrite)) {
      return TRUE;
    } else {
      set = set->next;
    }
  }
  return FALSE;
}

static SeqDescrPtr PropagateToSeqEntryIfNotDup (SeqEntryPtr sep, SeqDescrPtr sdp)
{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  SeqDescrPtr  new_sdp = NULL;
  
  if (sep != NULL && sep->data.ptrvalue != NULL) {
    if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (!IsSeqDescInSet (sdp, bsp->descr)) {
        new_sdp = AsnIoMemCopy ((Pointer) sdp,
                                (AsnReadFunc) SeqDescrAsnRead,
                                (AsnWriteFunc) SeqDescrAsnWrite);
        ValNodeLink (&(bsp->descr), new_sdp);
      }
    } else if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (!IsSeqDescInSet (sdp, bssp->descr)) {
        new_sdp = AsnIoMemCopy ((Pointer) sdp,
                                (AsnReadFunc) SeqDescrAsnRead,
                                (AsnWriteFunc) SeqDescrAsnWrite);
        ValNodeLink (&(bssp->descr), new_sdp);
      }
    }
  }
  return new_sdp;
}


static void AddSepBackToSet (SeqEntryPtr add_back_sep, BioseqSetPtr bssp, SeqDescrPtr parent_descriptors, Int4 set_pos)
{
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;
  SeqEntryPtr    sep, prev, tmp;
  SeqDescrPtr    sdp, new_sdp, sdp_tmp;
  ObjValNodePtr  ovp, new_ovp;
  Boolean        found;
  Int4           pos;
  BioseqPtr      bsp;

  /* any descriptors currently on bssp that are not in parent_descriptors must be propagated */
  /* mark parent descriptor as found, so we don't add it again */
  for (sdp = bssp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->extended == 0) {
      continue;
    }
    /* don't copy set titles */
    if (sdp->choice == Seq_descr_title) {
      continue;
    }
    ovp = (ObjValNodePtr) sdp;
    found = FALSE;
    for (new_sdp = parent_descriptors; new_sdp != NULL && !found; new_sdp = new_sdp->next) {
      if (new_sdp->extended == 0) {
        continue;
      }
      new_ovp = (ObjValNodePtr) new_sdp;
      if (new_ovp->idx.deleteme) {
        continue;
      }
      if (AsnIoMemComp (sdp, new_sdp, (AsnWriteFunc) SeqDescAsnWrite)) {
        found = TRUE;
        new_ovp->idx.deleteme = TRUE;
      }
    }
    if (!found) {
      /* propagate descriptor to members of bssp */
      sdp_tmp = sdp->next;
      sdp->next = NULL;
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        PropagateToSeqEntryIfNotDup (sep, sdp);
      }
      sdp->next = sdp_tmp;
      ovp->idx.deleteme = TRUE;
    }
  }
  /* now add any parent descriptors that weren't marked to the contents of add_back_sep
   * (and unmark those that were marked) */
  for (new_sdp = parent_descriptors; new_sdp != NULL; new_sdp = new_sdp->next) {
    if (new_sdp->extended == 0) {
      continue;
    }
    new_ovp = (ObjValNodePtr) new_sdp;
    if (new_ovp->idx.deleteme) {
      new_ovp->idx.deleteme = FALSE;
    } else {
      /* add to added_back_sep */
      PropagateToSeqEntryIfNotDup (add_back_sep, new_sdp);
    }
  }

  /* mark placeholders in destination set for removal */
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    if (IS_Bioseq (tmp) && (bsp = tmp->data.ptrvalue) != NULL && IsPlaceholderBioseq (bsp)) {
      bsp->idx.deleteme = TRUE;
    }
  }

  sep = SeqMgrGetSeqEntryForData (bssp);
  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);
  prev = NULL;
  pos = 0;
  for (tmp = bssp->seq_set; tmp != NULL && pos < set_pos; tmp = tmp->next) {
    prev = tmp;
    pos++;
  }
  if (prev == NULL) {
    add_back_sep->next = bssp->seq_set;
    bssp->seq_set = add_back_sep;
  } else {
    add_back_sep->next = tmp;
    prev->next = add_back_sep;
  }
  SeqMgrLinkSeqEntry (add_back_sep, OBJ_BIOSEQSET, bssp);

  SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
}


static void ReintegrateSequesteredSets (Uint2 entityID)
{
  BaseFormPtr bfp;
  BioseqViewFormPtr vfp;
  SeqEntryPtr new_parent_sep, sep, add_back_sep, next_sep;
  BioseqSetPtr   top_bssp, bssp;
  UnsequesterPtr u;

  add_back_sep = GetTopSeqEntryForEntityID (entityID);
  if (!IS_Bioseq_set (add_back_sep)) {
    return;
  }
  bfp = GetBaseFormForEntityID (entityID);
  RemoveSeqEntryViewer (bfp->form);

  /* remove any popset titles created */
  RemovePopsetTitles (add_back_sep);

  /* we created a genbank set to hold the actual parent sets */
  top_bssp = (BioseqSetPtr) add_back_sep->data.ptrvalue;

  /* need to process entries in seq_set in reverse, so that the set_pos will be 
   * correct relative to when the seq-entry was removed */
  ValNodeReverse (&(top_bssp->seq_set));
  new_parent_sep = top_bssp->seq_set;

  while (new_parent_sep != NULL) {
    bssp = new_parent_sep->data.ptrvalue;
    ValNodeReverse (&(bssp->seq_set));
    sep = bssp->seq_set;

    while (sep != NULL) {
      next_sep = sep->next;
      sep->next = NULL;
      u = GetUnsequesterForSep (sep, unsequester_list);
      if (u == NULL || u->orig_parent == NULL) {
        AddSepBackToSet (sep, sequester_set, bssp->descr, 0);
      } else {
        AddSepBackToSet (sep, u->orig_parent, bssp->descr, u->set_pos);
      }
      sep = next_sep;
    }
    bssp->seq_set = NULL;
    new_parent_sep = new_parent_sep->next;
  }

  ObjMgrDelete (OBJ_SEQENTRY, add_back_sep);
  top_bssp->seq_set = NULL;
  top_bssp->idx.deleteme = TRUE;
  DeleteMarkedObjects (entityID, 0, NULL);

  sep = GetTopSeqEntryForEntityID (sequester_set->idx.entityID);
  DeleteMarkedObjects (sequester_set->idx.entityID, 0, NULL);
  ObjMgrSetDirtyFlag (sequester_set->idx.entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, sequester_set->idx.entityID, 0, 0);

  bfp = GetBaseFormForEntityID (sequester_set->idx.entityID);
  if (bfp != NULL) {
    Show (bfp->form);
    vfp = (BioseqViewFormPtr) bfp;
    Show (vfp->toolForm);
  }

}


static void FinishedSequester (void)
{
  sequester_set = NULL;
  unsequester_list = ValNodeFreeData (unsequester_list);
  FileRemove(SEQUESTER_BACKUP_FILE);
  if (doReportAfterUnsequester > 0) {
    CreateReportWindow (doReportAfterUnsequester);
    doReportAfterUnsequester = 0;
  }
}


static void SequesterShutdown (ForM f) 
{
  MsgAnswer     ans;
  BaseFormPtr   bfp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    ans = Message (MSG_YNC, "Accept changes and add sequences back to original set?  Hit No to continue working, Cancel to restore original set");
    if (ans == ANS_YES) {
      ReintegrateSequesteredSets (bfp->input_entityID);
      FinishedSequester();
    } else if (ans == ANS_CANCEL) {
      RestoreSequesteredSets (f);
      FinishedSequester();
    }
  }
}


static void SequesterFormMessage (ForM f, Int2 mssg)

{
  BaseFormPtr   bfp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
      case VIB_MSG_QUIT :
      case VIB_MSG_ACCEPT :
        SequesterShutdown (f);
        break;
      default:
        SequinSeqViewFormMessage (f, mssg);
        break;
    }
  }
}


static void SequesterSets (ButtoN b)
{
  SegByFieldPtr sfp;
  Int4          num_chosen = 0;
  ValNodePtr       vnp;
  ClickableItemPtr cip;
  Int2             val;
  SeqEntryPtr      top_sep, pulled_sep;
  BioseqSetPtr     new_set;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;
  Uint2          entityID;
  FormMessageFunc old_message_handler;
  BaseFormPtr     bfp;
  BioseqViewFormPtr vfp;

  sfp = (SegByFieldPtr) GetObjectExtra (b);
  if (sfp == NULL) return;

  if (unsequester_list != NULL || sequester_set != NULL) {
    Message (MSG_ERROR, "You are already sequestering a set.  You cannot sequester until you have restored these sequences to the original set.");
    return;
  }

  WriteTheEntityID (sfp->input_entityID, SEQUESTER_BACKUP_FILE, FALSE);

  num_chosen = CountChosenDiscrepancies (sfp->value_lists, FALSE);
  if (num_chosen < 1) {
    if (ANS_OK == Message (MSG_OKC, "You have not chosen any groups.  Create all?")) {
      ChooseAllDiscrepancies (sfp->value_lists);
    } else {
      return;
    }
  }

  val = GetValue (sfp->seg_choice_grp);
  if (val == eSegPageId || val == eSegPageVectorContamination || val == eSegPageLength) {
    /* put all selected sequences into one group */
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    cip->chosen = TRUE;
    PullChosenIntoOneGroup (&(sfp->value_lists), cip);
    ValNodeAddPointer (&(sfp->value_lists), 0, cip);
  }

  sequester_set = sfp->target_set;

  /* create a new set with just the chosen sequences */
  top_sep = SeqMgrGetSeqEntryForData (sfp->target_set);
  SaveSeqEntryObjMgrData (top_sep, &omdptop, &omdata);
  GetSeqEntryParent (top_sep, &parentptr, &parenttype);

  new_set = BioseqSetNew ();
  new_set->_class = BioseqseqSet_class_genbank;

  for (vnp = sfp->value_lists; vnp != NULL; vnp = vnp->next) {      
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip == NULL || (!cip->chosen && ! AnyDiscrepanciesChosen (cip->subcategories))) {
      continue;
    }

    /* add SeqEntries for this category here */
    SequesterCategorySeqEntries (new_set, cip);

  }

  RestoreSeqEntryObjMgrData (top_sep, omdptop, &omdata); 
  DeleteMarkedObjects (sfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (sfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, sfp->input_entityID, 0, 0);

  pulled_sep = SeqEntryNew ();
  pulled_sep->choice = 2;
  pulled_sep->data.ptrvalue = new_set;
  SeqMgrLinkSeqEntry (pulled_sep, 0, NULL);
  
  entityID = ObjMgrRegister (OBJ_SEQENTRY, pulled_sep);
  seqviewprocs.filepath = NULL;
  seqviewprocs.forceSeparateViewer = TRUE;
  old_message_handler = seqviewprocs.handleMessages;
  seqviewprocs.handleMessages = SequesterFormMessage;


  SeqEntrySetScope (NULL);
  GatherProcLaunch (OMPROC_VIEW, FALSE, entityID, 1,
                              OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);

  seqviewprocs.handleMessages = old_message_handler;

  /* hide viewer with the remaining sequences */
  bfp = GetBaseFormForEntityID (sfp->input_entityID);
  if (bfp != NULL) {
    Hide (bfp->form);
    vfp = (BioseqViewFormPtr) bfp;
    Hide (vfp->toolForm);
  }

  if (GetStatus (sfp->leave_dlg_up)) {
    ChangeSegChoice (sfp->seg_choice_grp);
  } else {
    Remove (sfp->form);
  }
}

NLM_EXTERN Boolean OkToSequester (void)
{
  if (unsequester_list != NULL || sequester_set != NULL) {
    return FALSE;
  } else {
    return TRUE;
  }
}


NLM_EXTERN Uint2 SequesterClickableItem (Uint2 entityID, ClickableItemPtr cip)
{
  SeqEntryPtr sep;
  BioseqSetPtr bssp, new_set;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;
  Uint2          new_entityID;
  FormMessageFunc old_message_handler;
  BaseFormPtr     bfp;
  BioseqViewFormPtr vfp;
  SeqEntryPtr       pulled_sep;

  if (unsequester_list != NULL || sequester_set != NULL) {
    Message (MSG_ERROR, "You are already sequestering a set.  You cannot sequester until you have restored these sequences to the original set.");
    return 0;
  }

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL || !IS_Bioseq_set (sep)) {
    Message (MSG_ERROR, "This record does not have a top-levelset!");
    return 0;
  }

  bssp = FindTopLevelSetForDesktopFunction((BioseqSetPtr) sep->data.ptrvalue);
  if (bssp == NULL) {
    Message (MSG_ERROR, "Unable to find top level set");
    return 0;
  }

  WriteTheEntityID (entityID, SEQUESTER_BACKUP_FILE, FALSE);
  sequester_set = bssp;

  /* create a new set with just the chosen sequences */
  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);

  new_set = BioseqSetNew ();
  new_set->_class = BioseqseqSet_class_genbank;

  SequesterCategorySeqEntries (new_set, cip);

  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata); 
  DeleteMarkedObjects (entityID, 0, NULL);
  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);

  pulled_sep = SeqEntryNew ();
  pulled_sep->choice = 2;
  pulled_sep->data.ptrvalue = new_set;
  SeqMgrLinkSeqEntry (pulled_sep, 0, NULL);
  
  new_entityID = ObjMgrRegister (OBJ_SEQENTRY, pulled_sep);
  seqviewprocs.filepath = NULL;
  seqviewprocs.forceSeparateViewer = TRUE;
  old_message_handler = seqviewprocs.handleMessages;
  seqviewprocs.handleMessages = SequesterFormMessage;


  SeqEntrySetScope (NULL);
  GatherProcLaunch (OMPROC_VIEW, FALSE, new_entityID, 1,
                              OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);

  seqviewprocs.handleMessages = old_message_handler;

  /* hide viewer with the remaining sequences */
  bfp = GetBaseFormForEntityID (entityID);
  if (bfp != NULL) {
    Hide (bfp->form);
    vfp = (BioseqViewFormPtr) bfp;
    Hide (vfp->toolForm);
  }  
  return new_entityID;
}
  

NLM_EXTERN void LIBCALLBACK SequesterSequenceList (Uint2 entityID, ValNodePtr bsp_list)
{
  ClickableItemData cid;
  Uint2          new_entityID;
  BaseFormPtr     bfp;

  if (bsp_list == NULL) {
    return;
  }

  MemSet (&cid, 0, sizeof (ClickableItemData));
  cid.chosen = TRUE;
  cid.item_list = bsp_list;
  if ((new_entityID = SequesterClickableItem (entityID, &cid)) > 0) {
    /* get rid of existing validator window, replace with validator for new sequences */
    FreeValidateWindow ();
    bfp = GetBaseFormForEntityID (new_entityID);
    ValSeqEntryForm (bfp->form);
  }
}


extern Int2 LIBCALLBACK SequesterSequences (Pointer data)
{
  /* note - in future, will need to make sequester dialog modal. */
  return SegregateSetsByFieldEx (data, SequesterSets, eSegPageId);
}


extern void SequesterSequencesMenuItem (IteM i)
{
  BaseFormPtr bfp;
  SeqEntryPtr sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL || !IS_Bioseq_set (sep)) {
    Message (MSG_ERROR, "This record does not have a top-level set!");
  } else {
    NewSegregateBioseqSet (FindTopLevelSetForDesktopFunction((BioseqSetPtr) sep->data.ptrvalue), SequesterSets, eSegPageId);
  }
}


NLM_EXTERN void LIBCALLBACK SegregateSequenceList (Uint2 entityID, ValNodePtr bsp_list)
{
  SeqEntryPtr sep;
  BioseqSetPtr bssp;
  ClickableItemData cid;
  ValNodePtr list = NULL;

  if (bsp_list == NULL) {
    return;
  }

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL || !IS_Bioseq_set (sep)) {
    Message (MSG_ERROR, "This record does not have a top-level set!");
    return;
  }

  bssp = FindTopLevelSetForDesktopFunction((BioseqSetPtr) sep->data.ptrvalue);
  if (bssp == NULL) {
    Message (MSG_ERROR, "Unable to find top level set");
    return;
  }

  MemSet (&cid, 0, sizeof (ClickableItemData));
  cid.item_list = bsp_list;
  cid.chosen = TRUE;
  ValNodeAddPointer (&list, 0, &cid);

  MakeGroupsForUniqueValues (bssp, list);
  list = ValNodeFree (list);
}


static void FormatBoxColumn (ColPtr col)
{
  if (col == NULL) return;
  col->pixWidth = stdCharWidth;
  col->pixInset = 0;
  col->charWidth = 0;
  col->charInset = 0;
  col->font = NULL;
  col->just = 'l';
  col->wrap = 1;
  col->bar = 0;
  col->underline = 0;
  col->left = 0;
  col->last = FALSE;
}

static ColPtr AllocateDocColumnsForFieldList (BulkEdFieldPtr field_list, Int4Ptr p_len, Int4 last_sorted_col, Int4 display_offset, Int4 columns_shown)
{
  Int4   i, num_columns = 4, len = 0, col;
  ColPtr col_list;

  if (field_list == NULL) {
    return NULL;
  }

  /* set the font so that the pixwidth will be correct */
  SelectFont (programFont);

  /* start with 3 columns, one for expand/collapse and one for checkbox for row selection, one for right scroll, and one for final blank */
  for (i = 0; field_list[i].name != NULL; i++) {
    num_columns ++;
  }
  col_list = (ColPtr) MemNew (num_columns * sizeof (ColData));
  col = 0;
  FormatBoxColumn (col_list + col);
  len += col_list[col].pixWidth;
  col++;

  for (i = 0; field_list[i].name != NULL; i++) {
    if (i < display_offset) {
      continue;
    } else if (i >= display_offset + columns_shown) {
      break;
    }
    if (i == last_sorted_col) {
      FormatBoxColumn (col_list + col);
      len += col_list[col].pixWidth;
      col++;
    }
    len += (field_list[i].format_col_func) (col_list + col, field_list[i].name);
    col_list[col].last = FALSE;
    col++;
  }
  
  /* right scroll */
  FormatBoxColumn (col_list + col);
  len += col_list[col].pixWidth;
  col++;

  /* final blank column */
  col_list[col].pixWidth = 0;
  col_list[col].pixInset = 0;
  col_list[col].charWidth = 0;
  col_list[col].charInset = 0;
  col_list[col].font = NULL;
  col_list[col].just = 'l';
  col_list[col].wrap = 1;
  col_list[col].bar = 0;
  col_list[col].underline = 0;
  col_list[col].left = 0;
  col_list[col].last = TRUE;
  if (p_len != NULL) {
    *p_len = len;
  }
  return col_list;
}

typedef struct bulkeditorrow {
  ValNodePtr values_list;
  ValNodePtr object_list;
  struct bulkeditorrow PNTR subrows;
  Int4       copied_to_meta;
  Int4       copied_to_sub;
  Boolean    expanded;
  Boolean    selected;
} BulkEditorRowData, PNTR BulkEditorRowPtr;

static ValNodePtr FreeBulkEditorValues (ValNodePtr values, BulkEdFieldPtr field_list)
{
  ValNodePtr       vnp;
  Int4             col;

  for (col = 0, vnp = values;
       vnp != NULL;
       col++, vnp = vnp->next) {
    if (field_list[col].free_func == NULL) {
      vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
    } else {
      (field_list[col].free_func) (vnp->data.ptrvalue);
      vnp->data.ptrvalue = NULL;
    }
  }
  values = ValNodeFree (values);
  return values;
}

static void FreeBulkEditorValuesList (ValNodePtr PNTR values_list, Int4 num_rows, BulkEdFieldPtr field_list)
{
  Int4             i;

  if (values_list == NULL || field_list == NULL) return;
  for (i = 0; i < num_rows; i++) {
    values_list[i] = FreeBulkEditorValues (values_list[i], field_list);
  }
  values_list = MemFree (values_list);
}

static BulkEditorRowPtr FreeBulkEditorRows (BulkEditorRowPtr rows, BulkEdFieldPtr field_list)
{
  Int4 row_num;

  if (rows == NULL) return NULL;
  for (row_num = 0; rows[row_num].object_list != NULL; row_num++) {
    rows[row_num].values_list = FreeBulkEditorValues (rows[row_num].values_list, field_list);
    rows[row_num].object_list = ValNodeFree (rows[row_num].object_list);
    FreeBulkEditorRows (rows[row_num].subrows, field_list);
    rows[row_num].subrows = NULL;
  }
  rows = MemFree (rows);
  return rows;
}


static void CopyToSingleBulkEditRow (BulkEditorRowPtr dst, BulkEditorRowPtr src, BulkEdFieldPtr field_list)
{
  ValNodePtr vnp;
  Int4       col;

  if (dst == NULL || src == NULL) return;
  
  /* copy values */
  dst->values_list = NULL;
  for (vnp = src->values_list, col = 0; vnp != NULL; vnp = vnp->next, col++) {
    ValNodeAddPointer (&(dst->values_list), vnp->choice, field_list[col].copy_func(vnp->data.ptrvalue));
  }

  /* copy object */
  dst->object_list = NULL;
  ValNodeAddPointer (&(dst->object_list), src->object_list->choice, src->object_list->data.ptrvalue);

  /* copy selection */
  dst->selected = src->selected;
}


static void ExpandAllBulkEditRows (BulkEditorRowPtr row_list) 
{
  Int4 row_num;

  if (row_list == NULL) return;
  for (row_num = 0; row_list[row_num].object_list != NULL; row_num++) {
    if (row_list[row_num].subrows != NULL) {
      row_list[row_num].expanded = TRUE;
    }
  }
}

static void CollapseAllBulkEditRows (BulkEditorRowPtr row_list) 
{
  Int4 row_num;

  if (row_list == NULL) return;
  for (row_num = 0; row_list[row_num].object_list != NULL; row_num++) {
    if (row_list[row_num].subrows != NULL) {
      row_list[row_num].expanded = FALSE;
    }
  }
}


/* sort_col is zero-based data column, already corrected for selection and expansion columns */
static CharPtr GetTextForBulkEditorRow (BulkEditorRowPtr berp, BulkEdFieldPtr field_list, Int4 sort_col, Int4 display_offset, Int4 columns_shown)
{
  ValNodePtr column_values = NULL, vnp;
  Int4 text_len = 6, col;
  CharPtr str, line_text = NULL, tmp_str;
 
  if (berp == NULL) return NULL;
  for (vnp = berp->values_list, col = 0; vnp != NULL; vnp = vnp->next, col++) {
    if (col < display_offset) {
      continue;
    } else if (col >= display_offset + columns_shown) {
      break;
    }
    if (field_list[col].display_func == NULL) {
      str = NULL;
    } else {
      str = (field_list[col].display_func)(vnp->data.ptrvalue);
    }
    if (str == NULL) {
      str = StringSave ("");
    }
    if (berp->subrows != NULL && col == sort_col) {
      tmp_str = (CharPtr) MemNew (sizeof (Char) * (18 + StringLen (str)));
      sprintf (tmp_str, "%s%s(%d)", field_list[col].draw_col_func == NULL ? "" : " ",               
                                    str, ValNodeLen (berp->object_list));
      str = MemFree (str);
      str = tmp_str;
    }
    ValNodeAddPointer (&column_values, 0, str);
    text_len += StringLen (str) + 1;
  }
  line_text = (CharPtr) MemNew (text_len * sizeof (Char));
  line_text[0] = 0;
  /* add blank column for row selector */
  StringCat (line_text, "\t");
  for (vnp = column_values, col = display_offset; vnp != NULL; vnp = vnp->next, col++) {
    if (col == sort_col) {
      /* add blank column for expand/collapse */
      StringCat (line_text, "\t");
    }
    StringCat (line_text, vnp->data.ptrvalue);
    StringCat (line_text, "\t");
  }
  column_values = ValNodeFreeData (column_values);
  /* add blank column at end, to allow last real column to be highlighted when blank */
  StringCat (line_text, "\t\n");
  return line_text;
}


static Boolean DataPosFromBulkEdDlgRow (Int2 row, BulkEditorRowPtr row_list, Int4Ptr pRowNum, Int4Ptr pSubNum)
{
  Int4 row_num = 0, sub_num;
  Int2 display_row = 1;

  if (row < 1 || row_list == NULL || pRowNum == NULL || pSubNum == NULL) return FALSE;

  while (row_list[row_num].object_list != NULL) {
    if (display_row == row) {
      *pRowNum = row_num;
      *pSubNum = -1;
      return TRUE;
    }
    display_row++;
    if (row_list[row_num].subrows != NULL && row_list[row_num].expanded) {
      sub_num = 0;
      while (display_row != row && row_list[row_num].subrows[sub_num].object_list != NULL) {
        display_row++;
        sub_num++;
      }
      if (row_list[row_num].subrows[sub_num].object_list == NULL) {
        sub_num = -1;
      } else {
        if (display_row == row) {
          *pRowNum = row_num;
          *pSubNum = sub_num;
          return TRUE;
        }
      }
    }
    row_num++;
  }
  *pRowNum = -1;
  *pSubNum = -1;
  return FALSE;
}


static Int4 BulkEdDlgRowFromDataPos (Int4 row_num, Int4 sub_num, BulkEditorRowPtr row_list)
{
  Int4 i, j, display_row = 1;
  
  if (row_list == NULL) return 0;
  /* count subrows also passed */
  for (i = 0; row_list[i].object_list != NULL && i < row_num; i++) {
    display_row++;
    if (row_list[i].subrows != NULL && row_list[i].expanded) {
      for (j = 0; row_list[i].subrows[j].object_list != NULL; j++) {
        display_row++;
      }
    }
  }
  if (sub_num > -1) {
    display_row += sub_num + 1;
  }
  return display_row; 
}


static Int4 BulkEdDataColumnFromDocColumn (Int4 doc_col, Int4 sort_col, Int4 display_offset)
{
  /* subtract offset (document maps first column to 1) */
  doc_col--;
  /* subtract selection column */
  doc_col--;

  /* skip over columns not shown */
  doc_col += display_offset;

  /* subtract expand/collapse column */
  if (doc_col > sort_col && sort_col >= display_offset) {
    doc_col--;
  }
  return doc_col;
}


static ValNodePtr GetBulkEditorRowValueByColumn (BulkEditorRowPtr berp, Int4 col)
{
  Int4       i;
  ValNodePtr vnp;

  if (col < 0) return NULL;

  for (i = 0, vnp = berp->values_list;
       i < col && vnp != NULL;
       i++, vnp = vnp->next) {}
  return vnp;
}


static void SelectBulkEditorRows (BulkEditorRowPtr berp, Boolean val)
{
  if (berp == NULL) return;

  while (berp->object_list != NULL) {
    berp->selected = val;
    SelectBulkEditorRows (berp->subrows, val);
    berp++;
  }
}


static Int4 CountBulkEditorRowsSelected (BulkEditorRowPtr berp)
{
  Int4 num = 0;
  if (berp == NULL) return 0;

  while (berp->object_list != NULL) {
    if (berp->subrows == NULL) {
      if (berp->selected) {
        num++;
      }
    } else {
      num+= CountBulkEditorRowsSelected (berp->subrows);
    }
    berp++;
  }
  return num;
}


static void GetBulkEditorSelectedObjects (BulkEditorRowPtr berp, ValNodePtr PNTR object_list)
{
  ValNodePtr vnp;
  if (berp == NULL || object_list == NULL) return;

  while (berp->object_list != NULL) {
    if (berp->subrows == NULL) {
      if (berp->selected) {
        for (vnp = berp->object_list; vnp != NULL; vnp = vnp->next) {
          ValNodeAddPointer (object_list, vnp->choice, vnp->data.ptrvalue);
        }
      }
    } else {
      GetBulkEditorSelectedObjects (berp->subrows, object_list);
    }
    berp++;
  }
}


static Boolean AllBulkEditorRowsSelected (BulkEditorRowPtr berp)
{
  if (berp == NULL) return FALSE;

  while (berp->object_list != NULL) {
    /* note - don't need to check subrows, because metarow would
     * have been unselected if any subrows were unselected */
    if (!berp->selected) {
      return FALSE;
    }
    berp++;
  }
  return TRUE;
}


static void SelectBulkEditorRowsByStringConstraint (BulkEditorRowPtr berp, BulkEdFieldPtr field_list, Int4 col, StringConstraintXPtr scp, Boolean val)
{
  ValNodePtr vnp;
  CharPtr    str = NULL;
  Boolean    match;

  if (berp == NULL || field_list == NULL || col < 0 || field_list[col].display_func == NULL) return;

  while (berp->object_list != NULL) {
    if (berp->subrows == NULL) {
      vnp = GetBulkEditorRowValueByColumn (berp, col);
      if (vnp != NULL) {
        str = field_list[col].display_func(vnp->data.ptrvalue);
      }
      match = DoesStringMatchConstraintX (str, scp);
      str = MemFree (str);
      if (scp != NULL && scp->not_present) {
        match = !match;
      }
      if (match) {
        berp->selected = val;
      }
    } else {
      SelectBulkEditorRowsByStringConstraint (berp->subrows, field_list, col, scp, val);
      berp->selected = AllBulkEditorRowsSelected (berp->subrows);
    }
    berp++;
  }
}


/* call nulls and null values bigger */
/* return 0 if a and b are equal, -1 if a is smaller, and 1 if b is smaller */
static Int4 CompareBulkRowValues (Pointer val_a, Pointer val_b, BulkDisplayFieldFunc display_func)
{
  CharPtr str_a, str_b;
  Int4    rval = 0;

  if (val_a == NULL && val_b == NULL) {
    rval = 0;
  } else if (val_a == NULL) {
    rval = 1;
  } else if (val_b == NULL) {
    rval = -1;
  } else if (display_func == NULL) {
    rval = 0;
  } else {
    str_a = display_func(val_a);
    str_b = display_func(val_b);
    rval = StringCmp (str_a, str_b);
  }
  return rval;
}

/* call nulls and null values bigger */
/* return 0 if a and b are equal, -1 if a is smaller, and 1 if b is smaller */
static Int4 CompareBulkRows (BulkEditorRowPtr row_a, BulkEditorRowPtr row_b, Int4 sort_column, BulkEdFieldPtr field_list)
{
  Int4 col;
  ValNodePtr vnp_a, vnp_b;
  Int4 rval = 0;

  if (row_a == NULL && row_b == NULL) {
    rval = 0;
  } else if (row_a == NULL) {
    rval = 1;
  } else if (row_b == NULL) {
    rval = -1;
  } else if (row_a->values_list == NULL && row_b->values_list == NULL) {
    rval = 0;
  } else if (row_a->values_list == NULL) {
    rval = 1;
  } else if (row_b->values_list == NULL) {
    rval = -1;
  } else {
    col = 0;
    vnp_a = GetBulkEditorRowValueByColumn (row_a, sort_column);
    vnp_b = GetBulkEditorRowValueByColumn (row_b, sort_column);

    if (vnp_a == NULL && vnp_b == NULL) {
      rval = 0;
    } else if (vnp_a == NULL) {
      rval = 1;
    } else if (vnp_b == NULL) {
      rval = -1;
    } else {
      rval = CompareBulkRowValues (vnp_a->data.ptrvalue, vnp_b->data.ptrvalue, field_list[sort_column].display_func);
    }
  }
  return rval;
}

static Int4 CountIndividualRows (BulkEditorRowPtr row_list)
{
  Int4 num_rows = 0, row_num;

  if (row_list == NULL) return 0;

  for (row_num = 0; row_list[row_num].object_list != NULL; row_num++) {
    if (row_list[row_num].subrows == NULL) {
      num_rows++;
    } else {
      num_rows += CountIndividualRows(row_list[row_num].subrows);
    }
  }
  return num_rows;
}


static void SetDefaultCopiedToValues (BulkEditorRowPtr row_list)
{
  Int4 row_num;

  if (row_list == NULL) return;

  for (row_num = 0; row_list[row_num].object_list != NULL; row_num++) {
    row_list[row_num].copied_to_meta = -1;
    row_list[row_num].copied_to_sub = -1;
    if (row_list[row_num].subrows != NULL) {
      SetDefaultCopiedToValues(row_list[row_num].subrows);
    }
  }
}


static void SummarizeOneSubrowSet (BulkEditorRowPtr row, Int4 last_sorted_col)
{
  Int4       i, num_columns, row_num;
  ValNodePtr vnp;
  BoolPtr    all_present, any_present, all_same;
  CharPtr PNTR values;

  if (row == NULL || row->subrows == NULL) 
  {
    return;
  }
  /* count columns */
  for (num_columns = 0, vnp = row->values_list; vnp != NULL; num_columns++, vnp = vnp->next)
  {}

  all_present = (BoolPtr) MemNew (sizeof (Boolean) * num_columns);
  any_present = (BoolPtr) MemNew (sizeof (Boolean) * num_columns);
  all_same = (BoolPtr) MemNew (sizeof (Boolean) * num_columns);
  values = (CharPtr PNTR) MemNew (sizeof (CharPtr) * num_columns);

  for (i = 0; i < num_columns; i++) 
  {
    all_present[i] = TRUE;
    any_present[i] = FALSE;
    all_same[i] = TRUE;
    values[i] = NULL;
  }

  for (row_num = 0; row->subrows[row_num].object_list != NULL; row_num++) 
  {
    for (i = 0, vnp = row->subrows[row_num].values_list; vnp != NULL && i < num_columns; i++, vnp = vnp->next)
    {
      if (StringHasNoText (vnp->data.ptrvalue)) 
      {
        all_present[i] = FALSE;
      } 
      else 
      {
        any_present[i] = TRUE;
        if (values[i] == NULL) 
        {
          values[i] = vnp->data.ptrvalue;
        } 
        else if (StringCmp (values[i], vnp->data.ptrvalue) != 0) 
        {
          all_same[i] = FALSE;
        }
      }
    }
    while (i < num_columns) 
    {
      all_present[i] = FALSE;
      i++;
    }
  }

  for (i = 0, vnp = row->values_list; vnp != NULL; i++, vnp = vnp->next) 
  {
    if (i == last_sorted_col) {
      continue;
    }
    vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
    if (any_present[i]) {
      if (all_present[i]) {
        if (all_same[i]) {
          vnp->data.ptrvalue = StringSave (values[i]);
        } else {
          vnp->data.ptrvalue = StringSave ("All present, mixed");
        }
      } else if (all_same[i]) {
        vnp->data.ptrvalue = StringSave ("Some missing, one unique");
      } else {
        vnp->data.ptrvalue = StringSave ("Some missing, mixed");
      }
    }
  }
  all_present = MemFree (all_present);
  any_present = MemFree (any_present);
  all_same = MemFree (all_same);
  values = MemFree (values);
}


static void AddSummaryInfoToMetaRows (BulkEditorRowPtr row_list, Int4 last_sorted_col)
{
  Int4 row_num;

  if (row_list == NULL) return;

  for (row_num = 0; row_list[row_num].object_list != NULL; row_num++) {
    SummarizeOneSubrowSet(row_list + row_num, last_sorted_col);
    AddSummaryInfoToMetaRows(row_list[row_num].subrows, last_sorted_col);
  }
}


static BulkEditorRowPtr ResortBulkEditorRows (BulkEditorRowPtr orig_rows, BulkEdFieldPtr field_list, Int4 new_sort_column)
{
  Int4 num_individual_rows = 0;
  Int4 row_num, sub_num, last_inserted = 0, comp, col;
  BulkEditorRowPtr new_rows = NULL, smallest, berp;
  ValNodePtr       matches, vnp;
  Boolean          any_left, all_selected;

  /* count individual rows (will need this many or fewer) */
  num_individual_rows = CountIndividualRows(orig_rows);

  /* set default copied_to values */
  SetDefaultCopiedToValues (orig_rows);
  
  new_rows = (BulkEditorRowPtr) MemNew (sizeof (BulkEditorRowData) * (num_individual_rows + 1));
  MemSet (new_rows, 0, sizeof (BulkEditorRowData) * (num_individual_rows + 1));
  any_left = TRUE;
  while (any_left) {
    smallest = NULL;
    matches = NULL;
    for (row_num = 0; orig_rows[row_num].object_list != NULL; row_num++) {
      if (orig_rows[row_num].subrows != NULL) {
        /* look in subrows */
        for (sub_num = 0; orig_rows[row_num].subrows[sub_num].object_list != NULL; sub_num++) {
          if (orig_rows[row_num].subrows[sub_num].copied_to_meta > -1) continue; /* already copied to other list */
          if (smallest == NULL) {
            smallest = orig_rows[row_num].subrows + sub_num;
            ValNodeAddPointer (&matches, 0, orig_rows[row_num].subrows + sub_num);
          } else {
            comp = CompareBulkRows (orig_rows[row_num].subrows + sub_num, smallest, new_sort_column, field_list);
            if (comp == 0) {
              /* equal - add to list of matches */
              ValNodeAddPointer (&matches, 0, orig_rows[row_num].subrows + sub_num);
            } else if (comp < 0) {
              /* smaller - make new smallest and list of matches */
              smallest = orig_rows[row_num].subrows + sub_num;
              matches = ValNodeFree (matches);
              ValNodeAddPointer (&matches, 0, orig_rows[row_num].subrows + sub_num);
            }
            /* if it's bigger, leave it alone for now */
          }
        }
      } else {
        /* look at just this row */
        if (orig_rows[row_num].copied_to_meta > -1) continue; /* already copied to other list */
        if (smallest == NULL) {
          smallest = orig_rows + row_num;
          ValNodeAddPointer (&matches, 0, orig_rows + row_num);
        } else {
          comp = CompareBulkRows (orig_rows + row_num, smallest, new_sort_column, field_list);
          if (comp == 0) {
            /* equal - add to list of matches */
            ValNodeAddPointer (&matches, 0, orig_rows + row_num);
          } else if (comp < 0) {
            /* smaller - make new smallest and list of matches */
            smallest = orig_rows + row_num;
            matches = ValNodeFree (matches);
            ValNodeAddPointer (&matches, 0, orig_rows + row_num);
          }
          /* if it's bigger, leave it alone for now */
        }
      }
    }
    if (smallest == NULL) {
      any_left = FALSE;
    } else {
      if (ValNodeLen (matches) > 1) {
        /* make metarow */
        /* add value for only sort column */
        for (col = 0, vnp = smallest->values_list; vnp != NULL; col++, vnp = vnp->next) {
          if (col == new_sort_column) {
            ValNodeAddPointer (&(new_rows[last_inserted].values_list), vnp->choice, field_list[col].copy_func(vnp->data.ptrvalue));
          } else {
            ValNodeAddPointer (&(new_rows[last_inserted].values_list), 0, NULL);
          }
        }
        /* add features and subrows */
        new_rows[last_inserted].subrows = (BulkEditorRowPtr) MemNew (sizeof(BulkEditorRowData) * (ValNodeLen(matches) + 1));
        MemSet (new_rows[last_inserted].subrows, 0, sizeof(BulkEditorRowData) * (ValNodeLen(matches) + 1));
        all_selected = TRUE;
        for (vnp = matches, sub_num = 0; vnp != NULL; vnp = vnp->next, sub_num++) {
          berp = (vnp->data.ptrvalue);
          ValNodeAddPointer (&(new_rows[last_inserted].object_list), berp->object_list->choice, berp->object_list->data.ptrvalue);
          CopyToSingleBulkEditRow (new_rows[last_inserted].subrows + sub_num, berp, field_list);
          berp->copied_to_meta = last_inserted;
          berp->copied_to_sub = sub_num;
          if (!berp->selected) {
            all_selected = FALSE;
          }
        }
        new_rows[last_inserted].selected = all_selected;
      } else {
        CopyToSingleBulkEditRow (new_rows + last_inserted, smallest, field_list);
        smallest->copied_to_meta = last_inserted;
        smallest->copied_to_sub = 0;
      }
      last_inserted++;
      matches = ValNodeFree (matches);
    }
  }
  AddSummaryInfoToMetaRows (new_rows, new_sort_column);
  return new_rows;     
}


/* sort_column is zero-based data column, already corrected for selection and expansion columns */
static Boolean AllSubrowsSameValue (BulkEditorRowPtr berp, BulkEdFieldPtr field_list, Int4 sort_column)
{
  Int4 i;
  Boolean all_same = TRUE;

  if (berp == NULL || berp->subrows == NULL || berp->subrows[0].object_list == NULL) {
    return TRUE;
  }
  
  for (i = 1; berp->subrows[i].object_list != NULL && all_same; i++) {
    if (CompareBulkRows (berp->subrows + i - 1, berp->subrows + i, sort_column, field_list) != 0) {
      all_same = FALSE;
    }
  }

  return all_same;
}


/* sort_column is zero-based data column, already corrected for selection and expansion columns */
/* look to see if data in columns is different from data in editor */
static Boolean AnyChangeInBulkFields (BulkEditorRowPtr berp, BulkEdFieldPtr field_list, DialoG editor, Int4 sort_column)
{
  Pointer    m_data;
  ValNodePtr vnp;
  Int4       comp;

  if (berp == NULL || berp->subrows == NULL || berp->subrows[0].object_list == NULL) return FALSE;
  if (field_list == NULL || editor == NULL) return FALSE;

  if (!AllSubrowsSameValue (berp, field_list, sort_column)) {
    return TRUE;
  }

  /* find value for sort_column from first subrow */
  vnp = GetBulkEditorRowValueByColumn (berp, sort_column);
  
  if (vnp == NULL) return TRUE;
  
  /* compare with data from editor */
  m_data = DialogToPointer (editor);

  comp = CompareBulkRowValues (m_data, vnp->data.ptrvalue, field_list[sort_column].display_func);
  
  /* free data from editor */
  if (field_list[sort_column].free_func) {
    (field_list[sort_column].free_func)(m_data);
  }
  
  if (comp == 0) {
    return FALSE;
  } else {
    return TRUE;
  }
}


/* col is zero-based data column, already corrected for selection and expansion columns */
static void ApplyValueToOneBulkEditorRow (BulkEditorRowPtr berp, Int2 col, BulkEdFieldPtr field_list, DialoG editor)
{
  ValNodePtr vnp;

  if (berp == NULL || field_list == NULL ||  editor == NULL) return;

  vnp = GetBulkEditorRowValueByColumn (berp, col);

  if (vnp != NULL) {
    if (field_list[col].free_func == NULL) {
      vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
    } else {
      (field_list[col].free_func) (vnp->data.ptrvalue);
      vnp->data.ptrvalue = NULL;
    }
    vnp->data.ptrvalue = DialogToPointer (editor);
  }
}


static void PopulateSingleBulkEditRowFromData (BulkEditorRowPtr berp, Uint1 data_choice, Pointer data, BulkEdFieldPtr field_list, Pointer metadata)
{
  Int4 col;

  /* add values for all columns */
  for (col = 0; field_list[col].name != NULL; col++) {
    ValNodeAddPointer (&(berp->values_list), 0, (field_list[col].get_func (data_choice, data, metadata)));
  }
  /* add feature */
  ValNodeAddPointer (&(berp->object_list), data_choice, data);
}

static BulkEditorRowPtr MakeBulkEditorRowsFromObjectList (ValNodePtr object_list, BulkEdFieldPtr field_list, Int4 sort_column, Pointer metadata)
{
  ValNodePtr       vnp;
  BulkEditorRowPtr bulk_ed_rows = NULL, sorted_rows = NULL;
  Int4             num_rows, row_num;

  if (object_list == NULL || field_list == NULL) return NULL;

  /* first, create list with one row per object */
  num_rows = ValNodeLen (object_list);
  
  bulk_ed_rows = (BulkEditorRowPtr) MemNew (sizeof (BulkEditorRowData) * (num_rows + 1));
  MemSet (bulk_ed_rows, 0, sizeof (BulkEditorRowData) * (num_rows + 1));

  for (vnp = object_list, row_num = 0; vnp != NULL; vnp = vnp->next, row_num++)
  {
    PopulateSingleBulkEditRowFromData (bulk_ed_rows + row_num, vnp->choice, vnp->data.ptrvalue, field_list, metadata); 
  }

  if (sort_column >= 0)
  {
    /* Now resort rows for column */
    sorted_rows = ResortBulkEditorRows (bulk_ed_rows, field_list, sort_column);
    bulk_ed_rows = FreeBulkEditorRows (bulk_ed_rows, field_list);
    bulk_ed_rows = sorted_rows;
  }
  return bulk_ed_rows;
}


static Boolean IsBulkEditorColumnEmpty (BulkEditorRowPtr row_list, Int4 column);

static Boolean IsBulkEditorRowColumnEmpty (BulkEditorRowPtr row, Int4 column)
{
  ValNodePtr vnp;
  Int4 pos;

  if (row == NULL) {
    return TRUE;
  }

  for (pos = 0, vnp = row->values_list;
       pos < column && vnp != NULL;
       pos++, vnp = vnp->next)
  {
  }
  if (pos == column && vnp != NULL && !StringHasNoText (vnp->data.ptrvalue)) 
  {
    return FALSE;
  } 
  else 
  {
    return IsBulkEditorColumnEmpty(row->subrows, column);
  }
}


static Boolean IsBulkEditorColumnEmpty (BulkEditorRowPtr row_list, Int4 column)
{
  Boolean rval = TRUE;
  Int4    i;

  if (row_list == NULL) 
  {
    return TRUE;
  }
  for (i = 0; row_list[i].object_list != NULL && rval; i++) 
  {
    rval = IsBulkEditorRowColumnEmpty (row_list + i, column);
  }
  return rval;
}


static Int4 
AddRowsToListBySelection 
(BulkEditorRowPtr dest,
 BulkEditorRowPtr src,
 BulkEdFieldPtr   field_list,
 Int4             last_inserted,
 Boolean          selected)
{
  Int4 row_num;

  if (dest == NULL || src == NULL) return last_inserted;

  for (row_num = 0; src[row_num].object_list != NULL; row_num++) {
    if (src[row_num].subrows == NULL) {
      if ((src[row_num].selected && selected)
          || (!src[row_num].selected && !selected)) {
        CopyToSingleBulkEditRow (dest + last_inserted, src + row_num, field_list);
        src[row_num].copied_to_meta = last_inserted;
        src[row_num].copied_to_sub = 0;
        last_inserted++;
      }
    } else {
      last_inserted = AddRowsToListBySelection(dest, src[row_num].subrows, field_list, last_inserted, selected);
    }
  }
  return last_inserted;
}

static BulkEditorRowPtr ResortBulkEditorRowsBySelected (BulkEditorRowPtr orig_rows, BulkEdFieldPtr field_list)
{
  Int4 num_individual_rows;
  Int4 last_inserted = 0;
  BulkEditorRowPtr new_rows = NULL;

  /* count individual rows (will need this many or fewer) */
  num_individual_rows = CountIndividualRows(orig_rows);

  /* set default copied_to values */
  SetDefaultCopiedToValues (orig_rows);
  
  new_rows = (BulkEditorRowPtr) MemNew (sizeof (BulkEditorRowData) * (num_individual_rows + 1));
  MemSet (new_rows, 0, sizeof (BulkEditorRowData) * (num_individual_rows + 1));

  /* add selected rows first */
  last_inserted = AddRowsToListBySelection(new_rows, orig_rows, field_list, last_inserted, TRUE);

  /* add unselected rows */
  last_inserted = AddRowsToListBySelection(new_rows, orig_rows, field_list, last_inserted, FALSE);

  return new_rows;     
}





typedef enum {
  eBulkApplyField = 1,
  eBulkEditField,
  eBulkConvertField,
  eBulkParseField,
  eBulkSwapField,
  eBulkRemoveField,
  eBulkCapitalizeField
} EBulkAction;

static Int4 GetFieldNumFromPopupValue (Int4 popup_val, BulkEdFieldPtr field_list, Boolean edit_only)
{
  Int4 i, pval = 1;

  for (i = 0; field_list[i].name != NULL; i++) {
    if (field_list[i].create_dlg_func != NULL || (!edit_only && field_list[i].display_func != NULL)) {
      if (pval == popup_val) {
        return i;
      } else {
        pval++;
      }
    }
  }
  return -1;
}


static PopuP MakeBulkEditorFieldListPopup (GrouP g, BulkEdFieldPtr field_list, Boolean edit_only)
{
  PopuP p;
  Int4  i;

  if (field_list == NULL) return NULL;

  p = PopupList (g, TRUE, NULL);
  for (i = 0; field_list[i].name != NULL; i++) {
    if (field_list[i].create_dlg_func != NULL || (!edit_only && field_list[i].display_func != NULL)) {
      PopupItem (p, field_list[i].name);
    }
  }
  return p;
}


static DialoG MakeBulkEditorFieldListDialog (GrouP g, BulkEdFieldPtr field_list, Boolean edit_only)
{
  ValNodeBlock block;
  Int4  i;

  if (field_list == NULL) return NULL;
  InitValNodeBlock (&block, NULL);

  for (i = 0; field_list[i].name != NULL; i++) {
    if (field_list[i].create_dlg_func != NULL || (!edit_only && field_list[i].display_func != NULL)) {
      ValNodeAddPointerToEnd (&block, i, StringSave (field_list[i].name));
    }
  }
  return ValNodeSelectionDialog (g, block.head, TALL_SELECTION_LIST, ValNodeStringName, ValNodeSimpleDataFree, ValNodeStringCopy, ValNodeChoiceMatch, "field", NULL, NULL, FALSE);
}


static void SetBulkEditorFieldListDialogValue (DialoG d, Int4 val)
{
  ValNode vn;

  MemSet (&vn, 0, sizeof (ValNode));
  vn.choice = val;
  PointerToDialog (d, &vn);
}


static void ResetBulkEditorFieldListDialogValue (DialoG d, BulkEdFieldPtr field_list)
{
  SetBulkEditorFieldListDialogValue(d, GetFieldNumFromPopupValue(1, field_list, TRUE));
}


static Int4 GetBulkEditorFieldListDialogValue (DialoG d)
{
  Int4 val = 0;
  ValNodePtr vnp;

  vnp = DialogToPointer (d);
  if (vnp != NULL) {
    val = vnp->choice;
    vnp = ValNodeFree (vnp);
  }
  return val;
}



typedef struct onefieldbulkedit
{
  Int4          field;
  EditApplyPtr  edit_apply;
  ChangeCasePtr capitalization;
} OneFieldBulkEditData, PNTR OneFieldBulkEditPtr;


static OneFieldBulkEditPtr OneFieldBulkEditNew (void)
{
  OneFieldBulkEditPtr ofp;

  ofp = (OneFieldBulkEditPtr) MemNew (sizeof (OneFieldBulkEditData));
  return ofp;
}


static OneFieldBulkEditPtr OneFieldBulkEditFree (OneFieldBulkEditPtr ofp)
{
  if (ofp != NULL) {
    ofp->edit_apply = EditApplyFree (ofp->edit_apply);
    ofp->capitalization = MemFree (ofp->capitalization);
    ofp = MemFree (ofp);
  }
  return ofp;
}


static Int4 GetDlgFieldNumFromMasterFieldNum (BulkEdFieldPtr master_field_list, BulkEdFieldPtr field_list, Int4 master_field)
{
  Int4 m = 0, f = 0;

  while (master_field_list[m].name != NULL && m < master_field) {
    if (StringCmp (master_field_list[m].name, field_list[f].name) == 0) {
      f++;
    }
    m++;
  }

  return f;
}


static Int4 GetMasterFieldNumFromDlgFieldNum (BulkEdFieldPtr master_field_list, BulkEdFieldPtr field_list, Int4 dlg_field)
{
  Int4 m = 0;

  while (master_field_list[m].name != NULL && StringCmp (master_field_list[m].name, field_list[dlg_field].name) != 0) {
    m++;
  }

  return m;

}


static void AdjustOneFieldForMaster (OneFieldBulkEditPtr ofp, BulkEdFieldPtr master_field_list, BulkEdFieldPtr field_list)
{
  if (ofp == NULL || master_field_list == NULL || field_list == NULL) {
    return;
  }
  
  ofp->field = GetDlgFieldNumFromMasterFieldNum (master_field_list, field_list, ofp->field);
}


typedef CharPtr  (*Nlm_GetAutopopulateTextProc) PROTO ((Pointer, Int4));

typedef struct onefieldbulkeditdlg
{
  DIALOG_MESSAGE_BLOCK
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
  BulkEdFieldPtr           field_list;
  DialoG                   edit_apply;
  DialoG                   capitalization;
  DialoG                   field_dlg;

  Nlm_GetAutopopulateTextProc autopopulate_func;
  Pointer                     autopopulate_data;
} OneFieldBulkEditDialogData, PNTR OneFieldBulkEditDialogPtr;


static void OneFieldBulkEditToDialog (DialoG d, Pointer data)
{
  OneFieldBulkEditDialogPtr dlg;
  OneFieldBulkEditPtr       dlg_data;

  dlg = (OneFieldBulkEditDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  dlg_data = (OneFieldBulkEditPtr) data;
  if (dlg_data == NULL) {
    ResetBulkEditorFieldListDialogValue (dlg->field_dlg, dlg->field_list);
    PointerToDialog (dlg->edit_apply, NULL);
    PointerToDialog (dlg->capitalization, NULL);
  } else {
    SetBulkEditorFieldListDialogValue (dlg->field_dlg, dlg_data->field);
    PointerToDialog (dlg->edit_apply, dlg_data->edit_apply);
    PointerToDialog (dlg->capitalization, dlg_data->capitalization);
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static Pointer OneFieldBulkEditDialogToPointer (DialoG d)
{
  OneFieldBulkEditDialogPtr dlg;
  OneFieldBulkEditPtr       dlg_data;

  dlg = (OneFieldBulkEditDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  dlg_data = (OneFieldBulkEditPtr) MemNew (sizeof (OneFieldBulkEditData));
  dlg_data->field = GetBulkEditorFieldListDialogValue (dlg->field_dlg);
  dlg_data->edit_apply = DialogToPointer (dlg->edit_apply);
  dlg_data->capitalization = DialogToPointer (dlg->capitalization);

  return dlg_data;
}


static void OneFieldBulkEditMsg (DialoG d, Int2 mssg)

{
  OneFieldBulkEditDialogPtr dlg;

  dlg = (OneFieldBulkEditDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    if (mssg == AECR_VIB_MSG_CLEAR_TEXT) {
      SendMessageToDialog (dlg->edit_apply, NUM_VIB_MSG + 1);
    }
  }
}


static void OneFieldAutopopulate(ButtoN b)
{
  OneFieldBulkEditDialogPtr dlg;
  CharPtr                   txt;
  EditApplyPtr  edit_apply;
  Int4          field;

  dlg = (OneFieldBulkEditDialogPtr) GetObjectExtra (b);
  if (dlg != NULL && dlg->autopopulate_func != NULL) {
    field = GetBulkEditorFieldListDialogValue (dlg->field_dlg);
    txt = (dlg->autopopulate_func) (dlg->autopopulate_data, field);
    edit_apply = DialogToPointer (dlg->edit_apply);
    if (edit_apply == NULL) {
      edit_apply = EditApplyNew();
    }
    edit_apply->find_txt = MemFree (edit_apply->find_txt);
    edit_apply->find_txt = StringSave (txt);
    edit_apply->apply_txt = MemFree (edit_apply->apply_txt);
    edit_apply->apply_txt = StringSave (txt);
    PointerToDialog (dlg->edit_apply, edit_apply);
    edit_apply = EditApplyFree (edit_apply);
    txt = MemFree (txt);
  }
}


static void OneFieldClearText(ButtoN b)
{
  OneFieldBulkEditDialogPtr dlg;
  EditApplyPtr  edit_apply;

  dlg = (OneFieldBulkEditDialogPtr) GetObjectExtra (b);
  if (dlg != NULL) {
    edit_apply = DialogToPointer (dlg->edit_apply);
    if (edit_apply != NULL) {
      edit_apply->find_txt = MemFree (edit_apply->find_txt);
      edit_apply->repl_txt = MemFree (edit_apply->repl_txt);
      edit_apply->apply_txt = MemFree (edit_apply->apply_txt);
      PointerToDialog (dlg->edit_apply, edit_apply);
      edit_apply = EditApplyFree (edit_apply);
    }
  }
}


static DialoG OneFieldBulkEditDialog
(GrouP                    h,
 EBulkAction              action_choice,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 BulkEdFieldPtr           field_list,
 Nlm_GetAutopopulateTextProc autopopulate_func,
 Pointer                     autopopulate_data
)

{
  OneFieldBulkEditDialogPtr dlg;
  GrouP                     p, g1, g2 = NULL;
  ButtoN                    b;
  
  dlg = (OneFieldBulkEditDialogPtr) MemNew (sizeof (OneFieldBulkEditDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = OneFieldBulkEditToDialog;
  dlg->fromdialog = OneFieldBulkEditDialogToPointer;
  dlg->dialogmessage = OneFieldBulkEditMsg;
  dlg->testdialog = NULL;

  dlg->field_list = field_list;
  dlg->autopopulate_func = autopopulate_func;
  dlg->autopopulate_data = autopopulate_data;
  
  g1 = HiddenGroup (p, 2, 0, NULL);
  SetGroupSpacing (g1, 10, 10);
  dlg->field_dlg = MakeBulkEditorFieldListDialog (g1, field_list, TRUE);
  ResetBulkEditorFieldListDialogValue (dlg->field_dlg, dlg->field_list);

  if (action_choice == eBulkApplyField) 
  {
    dlg->edit_apply = EditApplyDialog (g1, eEditApplyChoice_Apply, NULL, NULL, NULL, NULL);
  } 
  else if (action_choice == eBulkEditField)
  {
    dlg->edit_apply = EditApplyDialog (g1, eEditApplyChoice_Edit, NULL, NULL, NULL, NULL);
  }
  else if (action_choice == eBulkCapitalizeField)
  {
    dlg->capitalization = ChangeCaseDialog (g1);
  }

  if (autopopulate_func != NULL) {
    g2 = HiddenGroup (p, 2, 0, NULL);
    SetGroupSpacing (g2, 10, 10);
    b = PushButton (g2, "Autopopulate", OneFieldAutopopulate);
    SetObjectExtra (b, dlg, NULL);
    b = PushButton (g2, "Clear Text", OneFieldClearText);
    SetObjectExtra (b, dlg, NULL);
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) g2, NULL);

  return (DialoG) p;
}


typedef struct twofieldbulkedit
{
  Int4 field_from;
  Int4 field_to;
  Boolean strip_name;
  Boolean leave_on_original;
  Boolean remove_parsed;
  TextPortionXPtr text_portion;
} TwoFieldBulkEditData, PNTR TwoFieldBulkEditPtr;


static TwoFieldBulkEditPtr TwoFieldBulkEditNew (void)
{
  TwoFieldBulkEditPtr tfp;

  tfp = (TwoFieldBulkEditPtr) MemNew (sizeof (TwoFieldBulkEditData));
  return tfp;
}


static TwoFieldBulkEditPtr TwoFieldBulkEditFree (TwoFieldBulkEditPtr tfp)
{
  if (tfp != NULL) {
    tfp->text_portion = TextPortionXFree (tfp->text_portion);
    tfp = MemFree (tfp);
  }
  return tfp;
}


static void AdjustTwoFieldsForMaster (TwoFieldBulkEditPtr tfp, BulkEdFieldPtr master_field_list, BulkEdFieldPtr field_list)
{
  if (tfp == NULL || master_field_list == NULL || field_list == NULL) {
    return;
  }
  
  tfp->field_to = GetDlgFieldNumFromMasterFieldNum (master_field_list, field_list, tfp->field_to);
  tfp->field_from = GetDlgFieldNumFromMasterFieldNum (master_field_list, field_list, tfp->field_from);
}


typedef struct twofieldbulkeditdlg
{
  DIALOG_MESSAGE_BLOCK
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
  BulkEdFieldPtr           field_list;
  ButtoN                   strip_name;
  ButtoN                   leave_on_original;
  ButtoN                   remove_parsed;
  DialoG                   text_portion;
  DialoG                   field_from;
  DialoG                   field_to;

  Nlm_GetAutopopulateTextProc autopopulate_func;
  Pointer                     autopopulate_data;
} TwoFieldBulkEditDialogData, PNTR TwoFieldBulkEditDialogPtr;


static void TwoFieldBulkEditToDialog (DialoG d, Pointer data)
{
  TwoFieldBulkEditDialogPtr dlg;
  TwoFieldBulkEditPtr       dlg_data;

  dlg = (TwoFieldBulkEditDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  dlg_data = (TwoFieldBulkEditPtr) data;
  if (dlg_data == NULL) {
    ResetBulkEditorFieldListDialogValue (dlg->field_from, dlg->field_list);
    ResetBulkEditorFieldListDialogValue (dlg->field_to, dlg->field_list);
    if (dlg->leave_on_original != NULL) {
      SetStatus (dlg->leave_on_original, FALSE);
    }
    if (dlg->strip_name != NULL) {
      SetStatus (dlg->strip_name, FALSE);
    }
    if (dlg->remove_parsed != NULL) {
      SetStatus (dlg->remove_parsed, FALSE);
    }
    PointerToDialog (dlg->text_portion, NULL);
  } else {
    SetBulkEditorFieldListDialogValue (dlg->field_from, dlg_data->field_from);
    SetBulkEditorFieldListDialogValue (dlg->field_to, dlg_data->field_from);
    if (dlg->leave_on_original != NULL) {
      SetStatus (dlg->leave_on_original, dlg_data->leave_on_original);
    }
    if (dlg->strip_name != NULL) {
      SetStatus (dlg->strip_name, dlg_data->strip_name);
    }
    if (dlg->remove_parsed != NULL) {
      SetStatus (dlg->remove_parsed, dlg_data->remove_parsed);
    }
    PointerToDialog (dlg->text_portion, dlg_data->text_portion);
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static Pointer TwoFieldBulkEditDialogToPointer (DialoG d)
{
  TwoFieldBulkEditDialogPtr dlg;
  TwoFieldBulkEditPtr       dlg_data;

  dlg = (TwoFieldBulkEditDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  dlg_data = (TwoFieldBulkEditPtr) MemNew (sizeof (TwoFieldBulkEditData));
  dlg_data->field_from = GetBulkEditorFieldListDialogValue (dlg->field_from);
  dlg_data->field_to = GetBulkEditorFieldListDialogValue (dlg->field_to);
  dlg_data->leave_on_original = dlg->leave_on_original == NULL ? FALSE : GetStatus (dlg->leave_on_original);
  dlg_data->strip_name = dlg->strip_name == NULL ? FALSE : GetStatus (dlg->strip_name);
  dlg_data->remove_parsed = dlg->remove_parsed == NULL ? FALSE : GetStatus (dlg->remove_parsed);
  dlg_data->text_portion = DialogToPointer (dlg->text_portion);

  return dlg_data;
}


static void TwoFieldBulkEditMsg (DialoG d, Int2 mssg)

{
  TwoFieldBulkEditDialogPtr dlg;

  dlg = (TwoFieldBulkEditDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    if (mssg == AECR_VIB_MSG_CLEAR_TEXT) {
      SendMessageToDialog (dlg->text_portion, NUM_VIB_MSG + 1);
    }
  }
}


static void TwoFieldAutopopulate(ButtoN b)
{
  TwoFieldBulkEditDialogPtr dlg;
  CharPtr                   txt;
  TextPortionXPtr  text_portion;
  Int4             field;

  dlg = (TwoFieldBulkEditDialogPtr) GetObjectExtra (b);
  if (dlg != NULL && dlg->autopopulate_func != NULL) {
    field = GetBulkEditorFieldListDialogValue (dlg->field_from);
    txt = (dlg->autopopulate_func) (dlg->autopopulate_data, field);
    text_portion = DialogToPointer (dlg->text_portion);
    if (text_portion == NULL) {
      text_portion = TextPortionXNew();
    }
    text_portion->start_text = MemFree (text_portion->start_text);
    text_portion->start_text = StringSave (txt);
    text_portion->end_text = MemFree (text_portion->end_text);
    text_portion->end_text = StringSave (txt);
    PointerToDialog (dlg->text_portion, text_portion);
    text_portion = TextPortionXFree (text_portion);
    txt = MemFree (txt);
  }
}


static void TwoFieldClearText(ButtoN b)
{
  TwoFieldBulkEditDialogPtr dlg;
  TextPortionXPtr  text_portion;

  dlg = (TwoFieldBulkEditDialogPtr) GetObjectExtra (b);
  if (dlg != NULL) {
    text_portion = DialogToPointer (dlg->text_portion);
    if (text_portion != NULL) {
      text_portion->start_text = MemFree (text_portion->start_text);
      text_portion->end_text = MemFree (text_portion->end_text);
      PointerToDialog (dlg->text_portion, text_portion);
      text_portion = TextPortionXFree (text_portion);
    }
  }
}


static DialoG TwoFieldBulkEditDialog
(GrouP                    h,
 EBulkAction              action_choice,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 BulkEdFieldPtr           field_list,
 Boolean                  strip_name,
 Nlm_GetAutopopulateTextProc autopopulate_func,
 Pointer                     autopopulate_data
)

{
  TwoFieldBulkEditDialogPtr dlg;
  GrouP            p, g1, g2, g3 = NULL;
  ButtoN           b;
  
  dlg = (TwoFieldBulkEditDialogPtr) MemNew (sizeof (TwoFieldBulkEditDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = TwoFieldBulkEditToDialog;
  dlg->fromdialog = TwoFieldBulkEditDialogToPointer;
  dlg->dialogmessage = TwoFieldBulkEditMsg;
  dlg->testdialog = NULL;

  dlg->field_list = field_list;
  dlg->autopopulate_func = autopopulate_func;
  dlg->autopopulate_data = autopopulate_data;
  
  g1 = HiddenGroup (p, 2, 0, NULL);
  StaticPrompt (g1, "From", 0, dialogTextHeight, systemFont, 'l');
  StaticPrompt (g1, "To", 0, dialogTextHeight, systemFont, 'l');
  
  dlg->field_from = MakeBulkEditorFieldListDialog (g1, field_list, TRUE);
  ResetBulkEditorFieldListDialogValue (dlg->field_from, field_list);
  dlg->field_to = MakeBulkEditorFieldListDialog (g1, field_list, TRUE);
  ResetBulkEditorFieldListDialogValue (dlg->field_to, field_list);

  if (action_choice == eBulkParseField) 
  {
    dlg->text_portion = TextPortionXDialog(p);
  }

  g2 = HiddenGroup (p, 3, 0, NULL);

  if (action_choice == eBulkParseField)
  {
    dlg->remove_parsed = CheckBox (g2, "Remove from parsed field", NULL);

    if (autopopulate_func != NULL) {
      g3 = HiddenGroup (p, 2, 0, NULL);
      SetGroupSpacing (g3, 10, 10);
      b = PushButton (g3, "Autopopulate", TwoFieldAutopopulate);
      SetObjectExtra (b, dlg, NULL);
      b = PushButton (g3, "Clear Text", TwoFieldClearText);
      SetObjectExtra (b, dlg, NULL);
    }
  }

  if (action_choice == eBulkConvertField)
  {
    dlg->leave_on_original = CheckBox (g2, "Leave on original", NULL);
  }
  if (strip_name)
  {
    dlg->strip_name = CheckBox (g2, "Strip name from text", NULL);
  }

  AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE)g2, (HANDLE) dlg->text_portion, (HANDLE) g3, NULL);
  return (DialoG) p;
}


static void BulkDrawCollapseExpand (Boolean val, RectPtr r)
{
  RecT rct;

  if (r == NULL) return;

  rct = *r;
  rct.bottom = rct.top + stdLineHeight - 4;
  rct.right = rct.left + stdLineHeight - 4;
  InsetRect (&rct, 2, 2);
  FrameRect (&rct);
  
  Gray();
  if (val) {
    /* draw minus */
    MoveTo (rct.left, (rct.top + rct.bottom) / 2);
    LineTo (rct.right - 1, (rct.top + rct.bottom) / 2);
  } else {
    /* draw plus */
    MoveTo (rct.left, (rct.top + rct.bottom) / 2);
    LineTo (rct.right - 1, (rct.top + rct.bottom) / 2);
    MoveTo ((rct.left + rct.right) / 2, rct.top);
    LineTo ((rct.left + rct.right) / 2, rct.bottom);
  }
  Black();
}


static void BulkDrawSelect (Boolean val, RectPtr r)
{
  RecT rct;

  if (r == NULL) return;

  rct = *r;
  rct.bottom = rct.top + stdLineHeight - 4;
  rct.right = rct.left + stdLineHeight - 4;
  InsetRect (&rct, 2, 2);
  FrameRect (&rct);
  
  if (val) {
    /* draw X */
    MoveTo (rct.left, rct.top);
    LineTo (rct.right - 1, rct.bottom - 1);
    MoveTo (rct.right, rct.top);
    LineTo (rct.left, rct.bottom - 1);
  }
}


typedef struct bulkeditordlg {
  DIALOG_MESSAGE_BLOCK
  DoC               title_doc;
  DoC               doc;
  BulkEdFieldPtr    master_field_list;
  BulkEdFieldPtr    field_list;
  ColPtr            col_list;
  ParData           par_fmt;
  BulkEditorRowPtr  row_list;
  DialoG PNTR       editor_list;
  GrouP             ed_btn_grp;
  Int4              num_columns;
  Int4              num_master_columns;
  Int2              last_viewed_row;
  Int2              last_viewed_col;  /* this will always refer to the data column */
  Int2              last_sorted_col;  /* this will always refer to the data column */
  /* for columns shown in display when there are extra */
  Int4              display_offset;  
  Int4              columns_shown;

  /* These are for the bulk edit actions */
  PopuP             bulk_action;
  GrouP             action_pages[8];

  DialoG            bulk_apply;
  DialoG            bulk_edit;
  DialoG            bulk_convert;
  DialoG            bulk_parse;
  DialoG            bulk_swap;
  DialoG            bulk_remove;
  DialoG            bulk_capitalize;

  GrouP             sel_unsel_grp;
  PopuP             check_constraint_field;
  DialoG            check_constraint;
  PrompT            check_status;

  Boolean           any_editing;
  Boolean           collapse_by_default;
  ClickableCallback single_click_func;
  ClickableCallback double_click_func;
  Pointer           metadata;
  SeqEntryPtr       sep;
} BulkEditorDlgData, PNTR BulkEditorDlgPtr;


static void InitBulkEdTitleDoc (BulkEditorDlgPtr dlg, DoC doc)
{
  Int4       i, text_len = 9;
  CharPtr    line_text;

  if (dlg == NULL || doc == NULL) {
    return;
  }

  Reset (doc);
  for (i = 0; i < dlg->display_offset + dlg->columns_shown; i++) {
    if (i < dlg->display_offset) {
      continue;
    }
    text_len += StringLen (dlg->field_list[i].name) + 1;
  }

  line_text = (CharPtr) MemNew (text_len * sizeof (Char));
  line_text[0] = 0;
  /* add blank column for row selection/scrolling */
  if (dlg->display_offset > 0) {
    StringCat (line_text, "<");
  }
  StringCat (line_text, "\t");

  for (i = 0; i < dlg->display_offset + dlg->columns_shown - 1; i++) {
    if (i < dlg->display_offset) {
      continue;
    }
    if (i == dlg->last_sorted_col) {
      /* add blank column for expand/collapse */
      StringCat (line_text, "\t");
    }
    StringCat (line_text, dlg->field_list[i].name);
    StringCat (line_text, "\t");
  }
  if (i == dlg->last_sorted_col) {
    /* add blank column for expand/collapse */
    StringCat (line_text, "\t");
  }
  StringCat (line_text, dlg->field_list[i].name);
  /* add blank column at end */
  if (dlg->display_offset + dlg->columns_shown < dlg->num_columns) {
    StringCat (line_text, "\t>\n");
  } else {
    StringCat (line_text, "\t\n");
  }

  AppendText (doc, line_text, &(dlg->par_fmt), dlg->col_list, programFont);
  line_text = MemFree (line_text);
}


/* redraw, set scroll position */
static void UpdateBulkEdDisplay (BulkEditorDlgPtr dlg)
{
  Int4       i, text_len = 0, sub_num;
  CharPtr    line_text;
  BaR        sb;
  Int4       offset = 0;
  Int2       first_shown = 0;

  if (dlg == NULL) {
    return;
  }

  /* Get original scroll position */
  sb = GetSlateVScrollBar ((SlatE)dlg->doc);
  offset = 0;
  GetScrlParams4 (dlg->doc, &offset, &first_shown, NULL);

  Reset (dlg->doc);

  if (dlg->row_list == NULL) {
    return;
  }

  for (i = 0; dlg->row_list[i].object_list != NULL; i++) {
    line_text = GetTextForBulkEditorRow(dlg->row_list + i, dlg->field_list, dlg->last_sorted_col, dlg->display_offset, dlg->columns_shown);
    AppendText (dlg->doc, line_text, &(dlg->par_fmt), dlg->col_list, programFont);
    line_text = MemFree (line_text);
    if (dlg->row_list[i].subrows != NULL && dlg->row_list[i].expanded) {
      for (sub_num = 0; dlg->row_list[i].subrows[sub_num].object_list != NULL; sub_num++) {
        line_text = GetTextForBulkEditorRow(dlg->row_list[i].subrows + sub_num, dlg->field_list, dlg->last_sorted_col, dlg->display_offset, dlg->columns_shown);
          AppendText (dlg->doc, line_text, &(dlg->par_fmt), dlg->col_list, programFont);
        line_text = MemFree (line_text);
      }
    }
  }
  
  /* restore original scroll position */
  GetItemParams4 (dlg->doc, first_shown, &offset, NULL, NULL, NULL, NULL);
  offset = MIN (offset, GetBarMax (sb));
  SetScrlParams4 (dlg->doc, offset);
}


static void RedrawAfterColumnShift (BulkEditorDlgPtr dlg)
{
  RecT r;
  Int4 width = 0;

  ResetClip();
  /* reset column widths */
  dlg->col_list = MemFree (dlg->col_list);
  SelectFont (systemFont);
  dlg->col_list = AllocateDocColumnsForFieldList (dlg->field_list, &width, dlg->last_sorted_col, dlg->display_offset, dlg->columns_shown);

  /* set data */
  InitBulkEdTitleDoc (dlg, dlg->title_doc);
  UpdateBulkEdDisplay (dlg);
  /* redraw highlight etc. */
  ObjectRect (dlg->title_doc, &r);
  InvalRect (&r);  
  ObjectRect (dlg->doc, &r);
  InvalRect (&r);  
  Update ();
}


/* col is zero-based data column, already corrected for selection and expansion columns */
static void BulkEdDlgToField (Int2 row, Int2 col, BulkEditorDlgPtr dlg)
{
  Int4       i;
  Int4       bulk_row, sub_num;
  Int2       edit_col, last_edit_col;

  if (row < 1 || col < 0 || dlg == NULL) return;

  if (!DataPosFromBulkEdDlgRow (row, dlg->row_list, &bulk_row, &sub_num)) return;

  /* don't apply values if editing master row and column other than sort */
  if (sub_num < 0 && col != dlg->last_sorted_col && dlg->row_list[bulk_row].subrows != NULL) return;

  edit_col = GetMasterFieldNumFromDlgFieldNum (dlg->master_field_list, dlg->field_list, col);

  if (dlg->editor_list[edit_col] == NULL) return;

  if (sub_num < 0) {
    /* if applying to all, get confirmation */
    if (dlg->row_list[bulk_row].subrows != NULL) {
      /* if nothing changed, no need to ask, nothing to do */
      last_edit_col = GetMasterFieldNumFromDlgFieldNum (dlg->master_field_list, dlg->field_list, dlg->last_sorted_col);
      if (!AnyChangeInBulkFields(dlg->row_list + bulk_row, dlg->field_list, dlg->editor_list[last_edit_col], dlg->last_sorted_col)) {
        return;
      } else if (ANS_CANCEL == Message (MSG_OKC, "Apply value to all features in collapsible group?")) {
        return;
      }
    }
    ApplyValueToOneBulkEditorRow (dlg->row_list + bulk_row, col, dlg->field_list, dlg->editor_list[edit_col]);
    if (dlg->row_list[bulk_row].subrows != NULL) {
      for (i = 0; dlg->row_list[bulk_row].subrows[i].object_list != NULL; i++) {
        ApplyValueToOneBulkEditorRow (dlg->row_list[bulk_row].subrows + i, col, dlg->field_list, dlg->editor_list[edit_col]);
      }
    }
  } else {
    ApplyValueToOneBulkEditorRow (dlg->row_list[bulk_row].subrows + sub_num, col, dlg->field_list, dlg->editor_list[edit_col]);
  }
  UpdateBulkEdDisplay (dlg);
}


/* col is zero-based data column, already corrected for selection and expansion columns */
/* if col is -1, hide everything (sorted by selection) */
static void ShowBulkFieldEditor (BulkEditorDlgPtr dlg, Int4 bulk_row, Int4 sub_num, Int2 col)
{
  BulkEditorRowPtr berp;
  ValNodePtr       vnp;
  Int2             i;
  Int4             edit_col;

  for (i = 0; i < dlg->num_master_columns; i++) {
    SafeHide (dlg->editor_list[i]);
  }
  SafeHide (dlg->ed_btn_grp);

  if (col < 0) return;

  if (sub_num < 0) {
    berp = dlg->row_list + bulk_row;
  } else {
    berp = dlg->row_list[bulk_row].subrows + sub_num;
  }

  /* only show editor dialog if sorted column or not master row */
  if (col == dlg->last_sorted_col || berp->subrows == NULL) {
    vnp = GetBulkEditorRowValueByColumn (berp, col);

    edit_col = GetMasterFieldNumFromDlgFieldNum (dlg->master_field_list, dlg->field_list, col);

    if (vnp == NULL) {
      PointerToDialog (dlg->editor_list[edit_col], NULL);
    } else {
      PointerToDialog (dlg->editor_list[edit_col], vnp->data.ptrvalue);
    }
    if (dlg->editor_list[edit_col] != NULL) {
      Show (dlg->editor_list[edit_col]);
      Show (dlg->ed_btn_grp);
    }
  }
}


NLM_EXTERN void SetBulkEditorSortColumn (DialoG dialog, Int4 col)
{
  BulkEditorDlgPtr dlg;
  Int4             bulk_row_num = -1, bulk_sub_num = -1;
  RecT             r;
  Int4             offset = 0, i;
  BulkEditorRowPtr sorted_rows;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (dialog);
  if (dlg == NULL) return;

  if (col < -1 || col >= dlg->num_columns) return;
  ResetClip();

  /* save value of column listed */
  dlg->last_sorted_col = col;

  if (col < 0) {
    sorted_rows = ResortBulkEditorRowsBySelected (dlg->row_list, dlg->field_list);
  } else {
    sorted_rows = ResortBulkEditorRows (dlg->row_list, dlg->field_list, dlg->last_sorted_col);
    if (dlg->collapse_by_default) {
      CollapseAllBulkEditRows (sorted_rows);
    } else {
      ExpandAllBulkEditRows (sorted_rows);
    }
  }
  /* get new position from old row */
  if (DataPosFromBulkEdDlgRow (dlg->last_viewed_row, dlg->row_list, &bulk_row_num, &bulk_sub_num)) {
    if (bulk_sub_num > -1) {
      dlg->last_viewed_row = BulkEdDlgRowFromDataPos (dlg->row_list[bulk_row_num].subrows[bulk_sub_num].copied_to_meta,
                                                      dlg->row_list[bulk_row_num].subrows[bulk_sub_num].copied_to_sub,
                                                      sorted_rows);
    } else {
      dlg->last_viewed_row = BulkEdDlgRowFromDataPos (dlg->row_list[bulk_row_num].copied_to_meta,
                                                      dlg->row_list[bulk_row_num].copied_to_sub,
                                                      sorted_rows);
    }
    dlg->last_viewed_col = dlg->last_sorted_col;
  }

  dlg->row_list = FreeBulkEditorRows (dlg->row_list, dlg->field_list);
  dlg->row_list = sorted_rows;


  /* reformat the columns */
  dlg->col_list = AllocateDocColumnsForFieldList (dlg->field_list, NULL, dlg->last_sorted_col, dlg->display_offset, dlg->columns_shown);
  ResetClip();
  InitBulkEdTitleDoc (dlg, dlg->title_doc);
  ObjectRect (dlg->title_doc, &r);
  InvalRect (&r);  

  UpdateBulkEdDisplay (dlg);

  if (dlg->last_viewed_row > 0) {
    /* show appropriate editor for highlighted item */
    if (DataPosFromBulkEdDlgRow (dlg->last_viewed_row, dlg->row_list, &bulk_row_num, &bulk_sub_num)) {
      ShowBulkFieldEditor (dlg, bulk_row_num, bulk_sub_num, dlg->last_sorted_col);
    } else {
      dlg->last_viewed_row = 0;
      dlg->last_viewed_col = 0;
      for (i = 0; i < dlg->num_master_columns; i++) {
        SafeHide (dlg->editor_list[i]);
      }
      Hide (dlg->ed_btn_grp);
    }

    /* scroll to new position of last_viewed */
    GetItemParams4 (dlg->doc, dlg->last_viewed_row, &offset, NULL, NULL, NULL, NULL);
    SetScrlParams4 (dlg->doc, offset);
  }

  /* redraw highlight etc. */
  ObjectRect (dlg->title_doc, &r);
  InvalRect (&r);  
  ObjectRect (dlg->doc, &r);
  InvalRect (&r);  
  Update ();
}


/* clicking on a column in the title indicates that we should sort by this column. */
static void ClickBulkTitle (DoC d, PoinT pt)
{
  BulkEditorDlgPtr dlg;
  Int2             item;
  Int2             row;
  Int2             col;

  MapDocPoint (d, pt, &item, &row, &col, NULL);
  if (item < 1) {
    return;
  }
  dlg = (BulkEditorDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  if (col == 1) {
    if (dlg->display_offset > 0) {
      dlg->display_offset --;
      RedrawAfterColumnShift (dlg);
      return;
    }
  } else if (dlg->col_list[col].last) {
    if (dlg->display_offset < dlg->num_columns - dlg->columns_shown) {
        dlg->display_offset ++;
        RedrawAfterColumnShift(dlg);
    }
    return;
  }


  if (col - 2 == dlg->last_sorted_col) {
    /* clicked on checkbox column */
    col = -1;
  } else {
    col = BulkEdDataColumnFromDocColumn(col, dlg->last_sorted_col, dlg->display_offset);
  }

  SetBulkEditorSortColumn (dlg->dialog, col);
}


static void UpdateCheckStatus (BulkEditorDlgPtr dlg)
{
  Int4 num_selected;
  CharPtr status_fmt = "%d feature%s currently checked";
  Char status_str[50];

  if (dlg == NULL) return;

  num_selected = CountBulkEditorRowsSelected (dlg->row_list);
  sprintf (status_str, status_fmt, num_selected, num_selected == 1 ? "" : "s");
  SetTitle (dlg->check_status, status_str);
}


/* when clicking on field, show editor (if available), change selection, collapse/expand,
 * set last viewed row, populate data from previously displayed editor 
 */
static void ClickBulkEdField (DoC d, PoinT pt)
{
  BulkEditorDlgPtr dlg;
  BulkEditorRowPtr berp;
  Int2            item;
  Int2            row;
  Int2            col;
  RecT            r;
  ValNodePtr      vnp;
  Int4            bulk_row = -1, sub_num = -1;
  Boolean         dblclick = dblClick;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  MapDocPoint (d, pt, &item, &row, &col, NULL);

  if (!DataPosFromBulkEdDlgRow (item, dlg->row_list, &bulk_row, &sub_num)) {
    return;
  }

  if (sub_num < 0) {
    berp = dlg->row_list + bulk_row;
  } else {
    berp = dlg->row_list[bulk_row].subrows + sub_num;
  }

  if (col == 1) {
    /* do selection */
    berp->selected = !berp->selected;
    if (berp->subrows != NULL) {
      SelectBulkEditorRows (berp->subrows, berp->selected);
    } else if (sub_num > -1) {
      /* this is a subrow - metarow can no longer be selected
       * if this is not selected.
       */
      if (dlg->row_list[bulk_row].selected && !berp->selected) {
        dlg->row_list[bulk_row].selected = FALSE;
      } else if (AllBulkEditorRowsSelected(dlg->row_list[bulk_row].subrows)) {
        /* all child rows are now selected, so select metarow */
        dlg->row_list[bulk_row].selected = TRUE;
      }
    }
    UpdateCheckStatus (dlg);
    ObjectRect (dlg->doc, &r);
    InvalRect (&r);  
    Update ();
    return;
  } else if (col == dlg->last_sorted_col + 2 - dlg->display_offset) {
    /* do expand/collapse */
    if (berp->subrows != NULL) {
      berp->expanded = !berp->expanded;
      UpdateBulkEdDisplay (dlg);
      ObjectRect (dlg->doc, &r);
      InvalRect (&r);  
      Update ();
    }
    return;
  } else  {
    /* adjust to account for checkbox column */
    col = BulkEdDataColumnFromDocColumn(col, dlg->last_sorted_col, dlg->display_offset);
  }
  if (item < 1 || col < 0 || col >= dlg->num_columns) return;
  /* took out this line - use Apply button instead, less confusing */
  /*BulkEdDlgToField (dlg->last_viewed_row, dlg->last_viewed_col, dlg); */
  ResetClip();

  ShowBulkFieldEditor (dlg, bulk_row, sub_num, col);

  if (dlg->field_list[col].release_cell_func != NULL) {
    vnp = GetBulkEditorRowValueByColumn (berp, col);
    if (vnp != NULL) {
      vnp->data.ptrvalue = (dlg->field_list[col].release_cell_func) (vnp->data.ptrvalue);
    }
  }
  dlg->last_viewed_row = item;
  dlg->last_viewed_col = col;

  /* redraw highlight etc. */
  ObjectRect (dlg->doc, &r);
  InvalRect (&r);  
  Update ();


  if (dblclick) {
    if (dlg->double_click_func != NULL) {
      (dlg->double_click_func) (berp->object_list, NULL);
    }
  } else {
    if (dlg->single_click_func != NULL) {
      (dlg->single_click_func) (berp->object_list, NULL);
    } 
  }
}


static void ApplyBulkEditingChange (ButtoN b)
{
  BulkEditorDlgPtr dlg;
  RecT r;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) return;

  BulkEdDlgToField (dlg->last_viewed_row, dlg->last_viewed_col, dlg);

  /* redraw highlight etc. */
  Select (dlg->dialog);
  ObjectRect (dlg->doc, &r);
  InvalRect (&r);  
  Update ();

}


static Boolean HighlightBulkEdField (DoC d, Int2 item, Int2 row, Int2 col)

{
  BulkEditorDlgPtr dlg;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return FALSE;

  if (col == 1 || col == dlg->last_sorted_col + 1) {
    /* don't highlight selection column or expansion column */
    return FALSE;
  } 

  col = BulkEdDataColumnFromDocColumn (col, dlg->last_sorted_col, dlg->display_offset);
  
  if (item == dlg->last_viewed_row && col == dlg->last_viewed_col) {
    return TRUE;
  } else {
    return FALSE;
  }
}

static void DrawBulkEdCell (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  Int2  lineHeight;
  RecT  rct;
  Int2  format_col, data_col;
  BulkEditorDlgPtr dlg;
  BulkEditorRowPtr berp;
  ValNodePtr       vnp;
  Int4             bulk_row_num = -1, bulk_sub_num = -1;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (d);
  if (dlg == NULL || r == NULL || item <= 0 || firstLine != 0) {
    return;
  }

  GetItemParams (d, item, NULL, NULL, NULL, &lineHeight, NULL);
  rct = *r;

  if (!DataPosFromBulkEdDlgRow (item, dlg->row_list, &bulk_row_num, &bulk_sub_num)) {
    return;
  }

  /* draw lines at top of each row after the first that does not have subrows */
  if (bulk_row_num > 0 && bulk_sub_num < 0) {
    MoveTo (rct.left, rct.top);
    LineTo (rct.right, rct.top);
  }
  
  if (bulk_sub_num < 0) {
    berp = dlg->row_list + bulk_row_num;
  } else {
    berp = dlg->row_list[bulk_row_num].subrows + bulk_sub_num;
  }
  
  /* for column 0, draw selection box */
  format_col = 0;
  rct.right = rct.left + dlg->col_list[format_col].pixWidth - 1;
  BulkDrawSelect (berp->selected, &rct);
  rct.left = rct.right + 1;
  format_col++; 
 
  for (data_col = 0, vnp = berp->values_list;
       data_col < dlg->num_columns && vnp != NULL; 
       data_col++, vnp = vnp->next) {

    if (data_col < dlg->display_offset) {
      continue;
    } else if (data_col >= dlg->display_offset + dlg->columns_shown) {
      break;
    }

    rct.right = rct.left + dlg->col_list[format_col].pixWidth - 1;
    if (data_col == dlg->last_sorted_col) {
      /* draw box if expandable */
      if (bulk_sub_num == -1
          && dlg->row_list[bulk_row_num].subrows != NULL) {
        BulkDrawCollapseExpand (dlg->row_list[bulk_row_num].expanded, &rct);
      }
      rct.left = rct.right + 1;
      format_col++;
      rct.right = rct.left + dlg->col_list[format_col].pixWidth - 1;
    }  
    if (dlg->field_list[data_col].draw_col_func != NULL) {
      (dlg->field_list[data_col].draw_col_func)(vnp->data.ptrvalue, &rct);
    }
    rct.left = rct.right + 1;
    format_col++;
  }
}


static CharPtr GetFirstBulkEditValue (Pointer data, Int4 col)
{
  BulkEditorDlgPtr dlg;
  CharPtr txt = NULL;
  Int4    i;

  if ((dlg = (BulkEditorDlgPtr) data) == NULL 
      || col < 0 || col >= dlg->num_master_columns 
      || dlg->master_field_list[col].get_func == NULL)
  {
    return NULL;
  }

  for (i = 0; dlg->row_list[i].object_list != NULL && txt == NULL; i++) 
  {
    txt = (dlg->master_field_list[col].get_func)(dlg->row_list[i].object_list->choice, dlg->row_list[i].object_list->data.ptrvalue, dlg->metadata);
    if (StringHasNoText (txt)) {
      txt = MemFree (txt);
    }
  }
  return txt;
}


static void CleanupBulkEditorDialog (GraphiC g, VoidPtr data)
{
  BulkEditorDlgPtr dlg;

  dlg = (BulkEditorDlgPtr) data;
  if (dlg != NULL) {
    dlg->row_list = FreeBulkEditorRows(dlg->row_list, dlg->field_list);
    /* field list has pointers to non-freeable data, just free container */
    dlg->field_list = MemFree (dlg->field_list);
    /* col_list is structures */
    dlg->col_list = MemFree (dlg->col_list);
  }
  StdCleanupExtraProc (g, data);
}

static void ApplyBulkEditorRowToObject (
  BulkEditorRowPtr berp, 
  ValNodePtr        target, 
  BulkEdFieldPtr   field_list, 
  Pointer          metadata, 
  LogInfoPtr       lip
)
{
  ValNodePtr vnp_val;
  Int4       i;
  Pointer    orig_val;
  CharPtr    orig_txt, new_txt, txt;

  if (berp == NULL || target == NULL || field_list == NULL) return;

  for (i = 0, vnp_val = berp->values_list;
       vnp_val != NULL;
       i++, vnp_val = vnp_val->next) {
    if (field_list[i].set_func != NULL) {
      if (lip != NULL) {
        orig_val = (field_list[i].get_func) (target->choice, target->data.ptrvalue, metadata);
        orig_txt = (field_list[i].display_func) (orig_val);
        new_txt = (field_list[i].display_func)(vnp_val->data.ptrvalue);
        if (StringCmp (orig_txt, new_txt) != 0) {
          txt = GetDiscrepancyItemText (target);
          fprintf (lip->fp, "Changed '%s' to '%s' for %s on %s\n",
                   orig_txt, new_txt, field_list[i].name, txt);
          txt = MemFree (txt);
          lip->data_in_log = TRUE;
        }
        orig_txt = MemFree (orig_txt);
        new_txt = MemFree (new_txt);
        field_list[i].free_func(orig_val);
      }
      (field_list[i].set_func) (target->data.ptrvalue, vnp_val->data.ptrvalue);
    }
  }
}

NLM_EXTERN void ApplyBulkEditorToObjectList (DialoG d, Boolean do_log)
{
  BulkEditorDlgPtr dlg;
  Int4             row, i;
  LogInfoPtr       lip = NULL;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  if (do_log) {
    lip = OpenLog ("Bulk Edit Changes");
  }

  for (row = 0; dlg->row_list[row].object_list != NULL; row++) {
    if (dlg->row_list[row].subrows == NULL) {
      ApplyBulkEditorRowToObject (dlg->row_list + row, 
                                  dlg->row_list[row].object_list, 
                                  dlg->field_list, 
                                  dlg->metadata,
                                  lip);
    } else {
      for (i = 0; dlg->row_list[row].subrows[i].object_list != NULL; i++) {
        ApplyBulkEditorRowToObject (dlg->row_list[row].subrows + i,
                                    dlg->row_list[row].subrows[i].object_list,
                                    dlg->field_list, 
                                    dlg->metadata,
                                    lip);
      }
    }
  }
  if (do_log) {
    CloseLog (lip);
    lip = FreeLog (lip);
  }
}


static void ExportBulkEditorRow (
  BulkEditorRowPtr berp, 
  ValNodePtr       target, 
  BulkEdFieldPtr   field_list, 
  Pointer          metadata,
  FILE *fp
)
{
  ValNodePtr vnp_val;
  Int4       i;
  Pointer    orig_val;
  CharPtr    orig_txt, new_txt;

  if (berp == NULL || target == NULL || field_list == NULL) return;

  /* first, print original values */
  for (i = 0, vnp_val = berp->values_list;
       vnp_val != NULL;
       i++, vnp_val = vnp_val->next) {
    orig_val = (field_list[i].get_func) (target->choice, target->data.ptrvalue, metadata);
    orig_txt = (field_list[i].display_func) (orig_val);
    fprintf (fp, "%s\t", orig_txt == NULL ? "" : orig_txt);
    orig_txt = MemFree (orig_txt);
    field_list[i].free_func(orig_val);
  }

  /* now print values that can be set */
  for (i = 0, vnp_val = berp->values_list;
       vnp_val != NULL;
       i++, vnp_val = vnp_val->next) {  
    if (field_list[i].set_str_func != NULL && field_list[i].set_func != NULL) {
      new_txt = (field_list[i].display_func)(vnp_val->data.ptrvalue);
      fprintf (fp, "%s\t", new_txt == NULL ? "" : new_txt);
      new_txt = MemFree (new_txt);
    }
  }
  fprintf (fp, "\n");
}


NLM_EXTERN void ExportBulkEditorTable (DialoG d)
{
  BulkEditorDlgPtr dlg;
  Int4             row, i;
  Char             path [PATH_MAX];
  FILE *fp;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  if (! GetOutputFileName (path, sizeof (path), NULL)) return;
  fp = FileOpen (path, "w");
  if (fp == NULL)
  {
    Message (MSG_ERROR, "Unable to open %s", path);
    return;
  }

  /* print header */
  /* first, print original values */
  for (i = 0; dlg->field_list[i].name != NULL; i++) {
    fprintf (fp, "Original %s\t", dlg->field_list[i].name);
  }
  /* now print values that can be changed */
  for (i = 0; dlg->field_list[i].name != NULL; i++) {
    if (dlg->field_list[i].set_str_func != NULL && dlg->field_list[i].set_func != NULL) {
      fprintf (fp, "New %s\t", dlg->field_list[i].name);
    }
  }
  fprintf (fp, "\n");

  /* now print table values */
  for (row = 0; dlg->row_list[row].object_list != NULL; row++) {
    if (dlg->row_list[row].subrows == NULL) {
      ExportBulkEditorRow (dlg->row_list + row, 
                                  dlg->row_list[row].object_list, 
                                  dlg->field_list, 
                                  dlg->metadata,
                                  fp);
    } else {
      for (i = 0; dlg->row_list[row].subrows[i].object_list != NULL; i++) {
        ExportBulkEditorRow (dlg->row_list[row].subrows + i,
                                    dlg->row_list[row].subrows[i].object_list,
                                    dlg->field_list, 
                                    dlg->metadata,
                                    fp);
      }
    }
  }
  FileClose (fp);
}


static void AddColumnValuesToRowList (BulkEditorRowPtr row_list, BulkEdFieldPtr master_field_list, Int4 master_pos, Int4 insert_pos, Pointer metadata)
{
  Int4 i, j;
  ValNodePtr vnp_prev, vnp;

  if (row_list == NULL || master_field_list == NULL) {
    return;
  }

  for (i = 0; row_list[i].object_list != NULL; i++) {
    vnp_prev = NULL;
    for (j = 0, vnp = row_list[i].values_list; j < insert_pos; j++, vnp = vnp->next) {
      vnp_prev = vnp;
    }
    vnp = ValNodeNew (NULL);
    vnp->choice = 0;
    if (row_list[i].subrows == NULL) {
      vnp->data.ptrvalue = master_field_list[master_pos].get_func (row_list[i].object_list->choice, row_list[i].object_list->data.ptrvalue, metadata);
    } 
    if (vnp_prev == NULL) {
      vnp->next = row_list[i].values_list;
      row_list[i].values_list = vnp;
    } else {
      vnp->next = vnp_prev->next;
      vnp_prev->next = vnp;
    }
    AddColumnValuesToRowList (row_list[i].subrows, master_field_list, master_pos, insert_pos, metadata);
  }
}


static Int4 AddColumnToRowList (BulkEditorRowPtr row_list, BulkEdFieldPtr master_field_list, BulkEdFieldPtr PNTR p_field_list, Int4 field, Int4 last_sort_col, Pointer metadata)
{
  Int4 i, m = 0, f = 0, num_columns = 1;
  BulkEdFieldPtr field_list, new_field_list;

  if (row_list == NULL || master_field_list == NULL || p_field_list == NULL || (field_list = *p_field_list) == NULL ||  field < 0) {
    return -1;
  }

  /* first, find insertion location */
  f = GetDlgFieldNumFromMasterFieldNum (master_field_list, *p_field_list, field);
  if (StringCmp (master_field_list[field].name, field_list[f].name) == 0) {
    /* already have column */
    return -1;
  }

  AddColumnValuesToRowList (row_list, master_field_list, field, f, metadata);

  /* insert column in field_list */
  for (i = 0; field_list[i].name != NULL; i++) {
    num_columns++;
  }

  num_columns++;
  new_field_list = (BulkEdFieldPtr) MemNew (sizeof (BulkEdFieldData) * num_columns);
  for (i = 0; i < f; i++) {
    new_field_list[i] = field_list[i];
  }
  new_field_list[i] = master_field_list[field];
  while (field_list[i].name != NULL) {
    new_field_list[i + 1] = field_list[i];
    i++;
  }
  new_field_list[i + 1] = field_list[i];
  field_list = MemFree (field_list);
  *p_field_list = new_field_list;

  AddSummaryInfoToMetaRows (row_list, last_sort_col);

  return f;
}


static void ChangeBulkEditorAction (PopuP p)
{
  BulkEditorDlgPtr dlg;
  Int4             page_num;
  Int4             num_pages;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (p);
  if (dlg == NULL || !dlg->any_editing) return;

  num_pages = sizeof (dlg->action_pages) / sizeof (GrouP);
  for (page_num = 0; page_num < num_pages; page_num++) {
    Hide (dlg->action_pages[page_num]);
  }
  page_num = GetValue (dlg->bulk_action) - 1;
  if (page_num > -1 && page_num < num_pages) {
    Show (dlg->action_pages[page_num]);
  }
}

static void RemoveColumnFromSelectedBulkEditRows (BulkEditorRowPtr row_list, Int4 col, BulkEdFieldPtr field_list, Boolean select_all)
{
  ValNodePtr vnp;

  if (row_list == NULL || col < 0 || field_list == NULL) return;
  while (row_list->object_list != NULL) {
    if (row_list->selected || select_all) {
      if (row_list->subrows == NULL || AllBulkEditorRowsSelected(row_list->subrows) || select_all) {
        vnp = GetBulkEditorRowValueByColumn (row_list, col);
        if (vnp != NULL) {
          field_list[col].free_func(vnp->data.ptrvalue);
          vnp->data.ptrvalue = NULL;
        }
      }
    }
    RemoveColumnFromSelectedBulkEditRows (row_list->subrows, col, field_list, select_all);
    row_list++;
  }
}

static void 
ApplyEditColumnFromSelectedBulkEditRows 
(BulkEditorRowPtr berp, 
 BulkEdFieldPtr   field_list,
 Int4             col,
 ApplyValuePtr    avp,
 Boolean select_all)
{
  ValNodePtr     vnp;

  if (berp == NULL) return;

  while (berp->object_list != NULL) {
    if (berp->selected || select_all) {
      if (berp->subrows == NULL || AllBulkEditorRowsSelected(berp->subrows) || select_all) {
        vnp = GetBulkEditorRowValueByColumn (berp, col);
        if (vnp != NULL) {
          vnp->data.ptrvalue = (field_list[col].set_str_func) (vnp->data.ptrvalue, avp);
        }
      }
    }
    ApplyEditColumnFromSelectedBulkEditRows (berp->subrows, field_list, col, avp, select_all);
    berp++;
  }
}


static void 
CapitalizeColumnsForSelectedBulkEditRows
(BulkEditorRowPtr berp, 
 BulkEdFieldPtr   field_list,
 Int4             col,
 ChangeCasePtr    capitalization,
 Boolean          select_all)
{
  ValNodePtr     vnp;
  CharPtr        val_str;
  ApplyValueData avd;
  SeqFeatPtr     sfp;
  SeqEntryPtr    sep, last_sep = NULL;
  ValNodePtr     orgname_list = NULL;

  if (berp == NULL) return;

  avd.field_list = NULL;
  avd.text_to_replace = NULL;
  avd.where_to_replace = EditApplyFindLocation_anywhere;
  avd.etp = NULL;

  while (berp->object_list != NULL) {
    if (berp->selected || select_all) {
      if (berp->subrows == NULL || AllBulkEditorRowsSelected(berp->subrows) || select_all) {
        vnp = GetBulkEditorRowValueByColumn (berp, col);
        if (vnp != NULL) {
          val_str = (field_list[col].display_func)(vnp->data.ptrvalue);
          if (!StringHasNoText (val_str)) {
            if (berp->object_list->choice == OBJ_SEQFEAT && berp->object_list->data.ptrvalue != NULL) {
              sfp = (SeqFeatPtr) berp->object_list->data.ptrvalue;
              sep = GetTopSeqEntryForEntityID (sfp->idx.entityID);
              if (sep != last_sep) {
                orgname_list = ValNodeFree (orgname_list);
                /* get org names to use in fixes */
                VisitBioSourcesInSep (sep, &orgname_list, GetOrgNamesInRecordCallback);
                last_sep = sep;
              }
            }
            ChangeCase (&val_str, capitalization, orgname_list);
            avd.new_text = val_str;
            vnp->data.ptrvalue = (field_list[col].set_str_func) (vnp->data.ptrvalue, &avd);
          }
          val_str = MemFree (val_str);
        }
      }
    }
    CapitalizeColumnsForSelectedBulkEditRows (berp->subrows, field_list, col, capitalization, select_all);
    berp++;
  }
  orgname_list = ValNodeFree (orgname_list);
}

static void RemoveNameFromText (CharPtr orig, CharPtr remove)
{
  CharPtr cp, cp_dst;
  Int4    len;

  if (StringHasNoText (orig) || StringHasNoText (remove)) {
    return;
  }
  cp = StringISearch (orig, remove);
  if (cp != NULL) {
    len = StringLen (remove);
    cp_dst = cp;
    cp = cp + len;
    while (*cp != 0 && len > 0) {
      *cp_dst = *cp;
      cp_dst++;
      cp++;
      len--;
    }
    *cp_dst = 0;
  }
}

static void ConvertColumnsForSelectedBulkEditRows 
(BulkEditorRowPtr    berp,
 BulkEdFieldPtr      field_list,
 TwoFieldBulkEditPtr tfp,
 ExistingTextPtr     etp,
 Boolean             select_all)
{
  ValNodePtr vnp_from, vnp_to;
  CharPtr    val_from;
  ApplyValueData avd;
  if (berp == NULL || tfp == NULL) return;

  avd.field_list = NULL;
  avd.text_to_replace = NULL;
  avd.where_to_replace = EditApplyFindLocation_anywhere;
  avd.etp = etp;

  while (berp->object_list != NULL) {
    if ((berp->selected || select_all) && berp->subrows == NULL) {
      vnp_from = GetBulkEditorRowValueByColumn (berp, tfp->field_from);
      vnp_to   = GetBulkEditorRowValueByColumn (berp, tfp->field_to);
      if (vnp_from != NULL && vnp_to != NULL) {
        val_from = (field_list[tfp->field_from].display_func)(vnp_from->data.ptrvalue);
        if (StringDoesHaveText (val_from)) {
          if (tfp->strip_name) {
            RemoveNameFromText (val_from, field_list[tfp->field_to].name);
          }

          avd.new_text = val_from;
          vnp_to->data.ptrvalue = (field_list[tfp->field_to].set_str_func) (vnp_to->data.ptrvalue, &avd);
          
          if (!tfp->leave_on_original) {
            (field_list[tfp->field_from].free_func)(vnp_from->data.ptrvalue);
            vnp_from->data.ptrvalue = NULL;
          }
      
        }
        val_from = MemFree (val_from);
      }
    }
    ConvertColumnsForSelectedBulkEditRows (berp->subrows, field_list, tfp, etp, select_all);
    berp++;
  }
}


static void ParseColumnsForSelectedBulkEditRows 
(BulkEditorRowPtr    berp,
 BulkEdFieldPtr      field_list,
 TwoFieldBulkEditPtr tfp,
 ExistingTextPtr     etp,
 Boolean             select_all)
{
  ValNodePtr vnp_from, vnp_to;
  CharPtr    val_from, val_to;
  ApplyValueData avd;
  if (berp == NULL || tfp == NULL) return;

  avd.field_list = NULL;
  avd.text_to_replace = NULL;
  avd.where_to_replace = EditApplyFindLocation_anywhere;
  avd.etp = etp;

  while (berp->object_list != NULL) {
    if ((berp->selected || select_all) && berp->subrows == NULL) {
      vnp_from = GetBulkEditorRowValueByColumn (berp, tfp->field_from);
      vnp_to   = GetBulkEditorRowValueByColumn (berp, tfp->field_to);
      if (vnp_from != NULL && vnp_to != NULL) {
        val_from = (field_list[tfp->field_from].display_func)(vnp_from->data.ptrvalue);
        val_to = NULL;
        if (StringDoesHaveText (val_from)) {
          val_to = ReplaceStringForParse(val_from, tfp->text_portion);
          if (StringDoesHaveText (val_to)) {
            avd.new_text = val_to;
            vnp_to->data.ptrvalue = (field_list[tfp->field_to].set_str_func) (vnp_to->data.ptrvalue, &avd);
          }
          if (tfp->remove_parsed) {
            avd.new_text = val_from;
            avd.etp = NULL;
            vnp_from->data.ptrvalue = (field_list[tfp->field_from].set_str_func) (vnp_from->data.ptrvalue, &avd);
            avd.etp = etp;
          }      
        }
        val_from = MemFree (val_from);
        val_to = MemFree (val_to);
      }
    }
    ParseColumnsForSelectedBulkEditRows (berp->subrows, field_list, tfp, etp, select_all);
    berp++;
  }
}


static void SwapColumnsForSelectedBulkEditRows
(BulkEditorRowPtr berp,
 BulkEdFieldPtr   field_list,
 Int4             col_from,
 Int4             col_to,
 Boolean          select_all)
{
  ValNodePtr vnp_from, vnp_to;
  CharPtr    val_from, val_to;
  ApplyValueData avd;
  if (berp == NULL) return;

  avd.etp = NULL;
  avd.field_list = NULL;
  avd.text_to_replace = NULL;
  avd.where_to_replace = EditApplyFindLocation_anywhere;

  while (berp->object_list != NULL) {
    if ((berp->selected || select_all) && berp->subrows == NULL) {
      vnp_from = GetBulkEditorRowValueByColumn (berp, col_from);
      vnp_to   = GetBulkEditorRowValueByColumn (berp, col_to);
      if (vnp_from != NULL && vnp_to != NULL) {
        val_from = (field_list[col_from].display_func)(vnp_from->data.ptrvalue);
        val_to = (field_list[col_to].display_func)(vnp_to->data.ptrvalue);
        if (StringDoesHaveText (val_from) || StringDoesHaveText (val_to)) {
          avd.new_text = val_to;
          vnp_from->data.ptrvalue = (field_list[col_from].set_str_func) (vnp_from->data.ptrvalue, &avd);
      
          avd.new_text = val_from;
          vnp_to->data.ptrvalue = (field_list[col_to].set_str_func) (vnp_to->data.ptrvalue, &avd);
        }
        val_from = MemFree (val_from);
        val_to = MemFree (val_to);
      }
    }
    SwapColumnsForSelectedBulkEditRows (berp->subrows, field_list, col_from, col_to, select_all);
    berp++;
  }
}



static void 
AddSampleDataForBulkEditorSelectedFeatures 
(BulkEditorRowPtr berp,
 BulkEdFieldPtr   field_list,
 Int4             col,
 Int4             col_src,
 GetSamplePtr     gsp,
 TextPortionXPtr   text_portion,
 Boolean           select_all)
{
  CharPtr    str, str_src, str_parse;
  ValNodePtr vnp, vnp_src;
  Boolean    skip;

  if (berp == NULL || field_list == NULL || col < 0 || field_list[col].display_func == NULL || gsp == NULL) return;
  
  while (berp->object_list != NULL) {
    if (berp->subrows == NULL && (berp->selected || select_all)) {
      vnp = GetBulkEditorRowValueByColumn (berp, col);
      if (vnp != NULL) {
        str = (field_list[col]).display_func (vnp->data.ptrvalue);
        if (!StringHasNoText (str)) {
          skip = FALSE;
          /* only add values if col_src < 0 or data in col_src is not blank */
          if (col_src >= 0) {
            skip = TRUE;
            vnp_src = GetBulkEditorRowValueByColumn (berp, col_src);
            if (vnp_src != NULL) {
              str_src = (field_list[col_src]).display_func(vnp_src->data.ptrvalue);
              if (StringDoesHaveText (str_src)) {
                if (text_portion == NULL) {
                  skip = FALSE;
                } else {
                  str_parse = ReplaceStringForParse(str_src, text_portion);
                  if (StringDoesHaveText (str_parse)) {
                    skip = FALSE;
                  }
                  str_parse = MemFree (str_parse);
                }
              }
              str_src = MemFree (str_src);
            }
          }
          if (!skip) {
            gsp->num_found++;
            if (gsp->sample_text == NULL) {
              gsp->sample_text = str;
              str = NULL;
            } else if (StringCmp (str, gsp->sample_text) != 0) {
              gsp->all_same = FALSE;
            }
          }
        }
        str = MemFree (str);
      }
    }
    AddSampleDataForBulkEditorSelectedFeatures (berp->subrows, field_list, col, col_src, gsp, text_portion, select_all);
    berp++;
  }
}


static GetSamplePtr 
GetSampleDataForBulkEditorSelectedFeatures 
(BulkEditorRowPtr berp,
 BulkEdFieldPtr   field_list,
 Int4             col,
 Int4             col_src,
 TextPortionXPtr   text_portion,
 Boolean           select_all)
{
  GetSamplePtr gsp;

  if (berp == NULL || field_list == NULL || col < 0) return NULL;

  gsp = GetSampleNew ();
 
  AddSampleDataForBulkEditorSelectedFeatures(berp, field_list, col, col_src, gsp, text_portion, select_all);

  return gsp;
}


static void AdjustBulkEditPositions (BulkEditorDlgPtr dlg, Int4 added_col)
{
  if (added_col > -1) {
    if (added_col <= dlg->last_viewed_col) {
      dlg->last_viewed_col++;
    }
    if (added_col <= dlg->last_sorted_col) {
      dlg->last_sorted_col++;
    }
    dlg->num_columns++;
  }
}


static void ApplyBulkEditCommon (BulkEditorDlgPtr dlg, Boolean select_all)
{
  RecT             r;
  GetSamplePtr     gsp = NULL;
  ExistingTextPtr  etp = NULL;
  ApplyValueData   avd;
  OneFieldBulkEditPtr ofp = NULL;
  TwoFieldBulkEditPtr tfp = NULL;
  Int4                bulk_row_num, bulk_sub_num, i;
  Int4                added_col1 = -1, added_col2 = -1;
  Boolean             need_redraw = FALSE;

  if (dlg == NULL) return;

  switch (GetValue (dlg->bulk_action)) {
    case eBulkRemoveField:
      ofp = (OneFieldBulkEditPtr) DialogToPointer (dlg->bulk_remove);
      if (ofp != NULL) {
        added_col1 = AddColumnToRowList (dlg->row_list, dlg->master_field_list, &(dlg->field_list), ofp->field, dlg->last_sorted_col, dlg->metadata);
        AdjustOneFieldForMaster (ofp, dlg->master_field_list, dlg->field_list);
        RemoveColumnFromSelectedBulkEditRows (dlg->row_list, ofp->field, dlg->field_list, select_all);
      }
      break;
    case eBulkApplyField:
      ofp = (OneFieldBulkEditPtr) DialogToPointer (dlg->bulk_apply);
      if (ofp != NULL) {
        /* add field to table, if it's not in the table yet */
        added_col1 = AddColumnToRowList (dlg->row_list, dlg->master_field_list, &(dlg->field_list), ofp->field, dlg->last_sorted_col, dlg->metadata);
        AdjustOneFieldForMaster (ofp, dlg->master_field_list, dlg->field_list);

        /* warn about existing values */
        gsp = GetSampleDataForBulkEditorSelectedFeatures (dlg->row_list, dlg->field_list, ofp->field, -1, NULL, select_all);
        etp = GetExistingTextHandlerInfo (gsp == NULL ? 0 : gsp->num_found, FALSE);
        if (etp == NULL || etp->existing_text_choice != eExistingTextChoiceCancel) {
          avd.etp = etp;
          avd.field_list = NULL;
          AddEditApplyDataToApplyValue (AECR_APPLY, ofp->edit_apply, &avd);
          ApplyEditColumnFromSelectedBulkEditRows (dlg->row_list, dlg->field_list, ofp->field, &avd, select_all);
        }
      }
      break;
    case eBulkEditField:
      ofp = (OneFieldBulkEditPtr) DialogToPointer (dlg->bulk_edit);
      if (ofp != NULL) {
        added_col1 = AddColumnToRowList (dlg->row_list, dlg->master_field_list, &(dlg->field_list), ofp->field, dlg->last_sorted_col, dlg->metadata);
        AdjustOneFieldForMaster (ofp, dlg->master_field_list, dlg->field_list);
        avd.etp = NULL;
        avd.field_list = NULL;
        AddEditApplyDataToApplyValue (AECR_EDIT, ofp->edit_apply, &avd);
        ApplyEditColumnFromSelectedBulkEditRows (dlg->row_list, dlg->field_list, ofp->field, &avd, select_all);
      }
      break;     
    case eBulkCapitalizeField:
      ofp = (OneFieldBulkEditPtr) DialogToPointer (dlg->bulk_capitalize);
      if (ofp != NULL && ofp->capitalization != NULL && ofp->capitalization->change != eChangeCaseNone) {
        added_col1 = AddColumnToRowList (dlg->row_list, dlg->master_field_list, &(dlg->field_list), ofp->field, dlg->last_sorted_col, dlg->metadata);
        AdjustOneFieldForMaster (ofp, dlg->master_field_list, dlg->field_list);
        CapitalizeColumnsForSelectedBulkEditRows (dlg->row_list, dlg->field_list, ofp->field, ofp->capitalization, select_all);
      }
      break;
    case eBulkConvertField:
      tfp = (TwoFieldBulkEditPtr) DialogToPointer (dlg->bulk_convert);
      if (tfp != NULL) {
        added_col1 = AddColumnToRowList (dlg->row_list, dlg->master_field_list, &(dlg->field_list), tfp->field_to, dlg->last_sorted_col, dlg->metadata);
        added_col2 = AddColumnToRowList (dlg->row_list, dlg->master_field_list, &(dlg->field_list), tfp->field_from, dlg->last_sorted_col, dlg->metadata);
        AdjustTwoFieldsForMaster (tfp, dlg->master_field_list, dlg->field_list);

        gsp = GetSampleDataForBulkEditorSelectedFeatures (dlg->row_list, dlg->field_list, tfp->field_to, tfp->field_from, NULL, select_all);
        etp = GetExistingTextHandlerInfo (gsp == NULL ? 0 : gsp->num_found, FALSE);
        if (etp == NULL || etp->existing_text_choice != eExistingTextChoiceCancel) {
          ConvertColumnsForSelectedBulkEditRows (dlg->row_list, dlg->field_list, tfp, etp, select_all);
        }
      }
      break;
    case eBulkParseField:
      tfp = (TwoFieldBulkEditPtr) DialogToPointer (dlg->bulk_parse);
      if (tfp != NULL) {
        added_col1 = AddColumnToRowList (dlg->row_list, dlg->master_field_list, &(dlg->field_list), tfp->field_to, dlg->last_sorted_col, dlg->metadata);
        added_col2 = AddColumnToRowList (dlg->row_list, dlg->master_field_list, &(dlg->field_list), tfp->field_from, dlg->last_sorted_col, dlg->metadata);
        AdjustTwoFieldsForMaster (tfp, dlg->master_field_list, dlg->field_list);
        gsp = GetSampleDataForBulkEditorSelectedFeatures (dlg->row_list, dlg->field_list, tfp->field_to, tfp->field_from, tfp->text_portion, select_all);
        etp = GetExistingTextHandlerInfo (gsp == NULL ? 0 : gsp->num_found, FALSE);        
        if (etp == NULL || etp->existing_text_choice != eExistingTextChoiceCancel) {
          ParseColumnsForSelectedBulkEditRows (dlg->row_list, dlg->field_list, tfp, etp, select_all);
        }        
      }
      break;
    case eBulkSwapField:
      tfp = (TwoFieldBulkEditPtr) DialogToPointer (dlg->bulk_swap);
      if (tfp != NULL) {
        added_col1 = AddColumnToRowList (dlg->row_list, dlg->master_field_list, &(dlg->field_list), tfp->field_to, dlg->last_sorted_col, dlg->metadata);
        added_col2 = AddColumnToRowList (dlg->row_list, dlg->master_field_list, &(dlg->field_list), tfp->field_from, dlg->last_sorted_col, dlg->metadata);
        AdjustTwoFieldsForMaster (tfp, dlg->master_field_list, dlg->field_list);
        SwapColumnsForSelectedBulkEditRows (dlg->row_list, dlg->field_list, tfp->field_from, tfp->field_to, select_all);
      }
      break;
  }

  ofp = OneFieldBulkEditFree (ofp);
  tfp = TwoFieldBulkEditFree (tfp);
  gsp = GetSampleFree (gsp);
  etp = MemFree (etp);

  AdjustBulkEditPositions (dlg, added_col1);
  AdjustBulkEditPositions (dlg, added_col2);

  AddSummaryInfoToMetaRows (dlg->row_list, dlg->last_sorted_col);

  /* update individual field editor */
  if (dlg->last_viewed_row > 0) {
    /* show appropriate editor for highlighted item */
    if (DataPosFromBulkEdDlgRow (dlg->last_viewed_row, dlg->row_list, &bulk_row_num, &bulk_sub_num)) {
      ShowBulkFieldEditor (dlg, bulk_row_num, bulk_sub_num, dlg->last_viewed_col);
    } else {
      dlg->last_viewed_row = 0;
      dlg->last_viewed_col = 0;
      for (i = 0; i < dlg->num_master_columns; i++) {
        SafeHide (dlg->editor_list[i]);
      }
      Hide (dlg->ed_btn_grp);
    }
  }

  if (added_col1 > -1 || added_col2 > -1) {
    RedrawAfterColumnShift (dlg);
  } else {
    UpdateBulkEdDisplay (dlg);
    ObjectRect (dlg->doc, &r);
    InvalRect (&r);  
    Update ();
  }
}


static void ApplyBulkEditAll (ButtoN b)
{
  BulkEditorDlgPtr dlg;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) return;

  ApplyBulkEditCommon (dlg, TRUE);
}

static void ApplyBulkEdit (ButtoN b)
{
  BulkEditorDlgPtr dlg;
  Boolean             select_all = FALSE;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) return;

  if (CountBulkEditorRowsSelected(dlg->row_list) == 0) {
    if (ANS_CANCEL == Message (MSG_OKC, "No items selected!  Apply to all?")) {
      return;
    } else {
      select_all = TRUE;
    }
  }

  ApplyBulkEditCommon (dlg, select_all);
}


extern void BulkEditorCheckAllDialog (DialoG d)
{
  BulkEditorDlgPtr dlg;
  RecT             r;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  SelectBulkEditorRows (dlg->row_list, TRUE);
  UpdateCheckStatus (dlg);
  ObjectRect (dlg->doc, &r);
  InvalRect (&r);  
  Update ();
}

static void BulkEditorCheckAll (ButtoN b)
{
  BulkEditorDlgPtr dlg;
  RecT             r;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) return;

  SelectBulkEditorRows (dlg->row_list, TRUE);
  UpdateCheckStatus (dlg);
  ObjectRect (dlg->doc, &r);
  InvalRect (&r);  
  Update ();
}


static void BulkEditorUncheckAll (ButtoN b)
{
  BulkEditorDlgPtr dlg;
  RecT             r;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) return;

  SelectBulkEditorRows (dlg->row_list, FALSE);
  UpdateCheckStatus (dlg);
  ObjectRect (dlg->doc, &r);
  InvalRect (&r);  
  Update ();
}

static void BulkEditorStringConstraintCheck (ButtoN b)
{
  BulkEditorDlgPtr    dlg;
  RecT                r;
  Int4                field_num;
  StringConstraintXPtr scp;
  Boolean             val;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) return;

  field_num = GetFieldNumFromPopupValue(GetValue (dlg->check_constraint_field), dlg->field_list, FALSE);
  scp = DialogToPointer (dlg->check_constraint);

  if (GetValue (dlg->sel_unsel_grp) == 1) {
    val = TRUE;
  } else {
    val = FALSE;
  }

  SelectBulkEditorRowsByStringConstraint (dlg->row_list, dlg->field_list, field_num, scp, val);
  
  scp = StringConstraintXFree (scp);
  UpdateCheckStatus (dlg);

  ObjectRect (dlg->doc, &r);
  InvalRect (&r);  
  Update ();
}


static void ObjectListToBulkEditorDialog (DialoG d, Pointer data)
{
  BulkEditorDlgPtr dlg;
  ValNodePtr       object_list;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  object_list = (ValNodePtr) data;
    
  dlg->row_list = FreeBulkEditorRows(dlg->row_list, dlg->field_list);
  dlg->row_list = MakeBulkEditorRowsFromObjectList (object_list, dlg->field_list, dlg->last_sorted_col, dlg->metadata);
  if (dlg->collapse_by_default) {
    CollapseAllBulkEditRows (dlg->row_list);
  } else {
    ExpandAllBulkEditRows (dlg->row_list);
  }
  ResetClip();
  UpdateBulkEdDisplay (dlg);
  UpdateCheckStatus (dlg);
}


static Pointer GetSelectedObjectsFromBulkEditorDialog (DialoG d)
{
  BulkEditorDlgPtr dlg;
  ValNodePtr       object_list = NULL;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  GetBulkEditorSelectedObjects (dlg->row_list, &object_list);
  return (Pointer) object_list;
}


static void CopyEditToClipboard (ButtoN b)
{
  BulkEditorDlgPtr dlg;
  CharPtr str;
  Int2    edit_col;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }
  if (dlg->last_viewed_col < 0 
    || dlg->last_viewed_col >= dlg->num_columns) {
    return;
  }
  edit_col = GetMasterFieldNumFromDlgFieldNum (dlg->master_field_list, dlg->field_list, dlg->last_viewed_col);
  if (dlg->editor_list[edit_col] == NULL) {
    return;
  }
  str = DialogToPointer (dlg->editor_list[edit_col]);
  if (!StringHasNoText (str)) {
    StringToClipboard (str);
  }
  str = MemFree (str);
}


static void CopyEditFromClipboard (ButtoN b)
{
  BulkEditorDlgPtr dlg;
  CharPtr str;
  Int2    edit_col;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }
  if (dlg->last_viewed_col < 0 
    || dlg->last_viewed_col >= dlg->num_columns) {
    return;
  }
  edit_col = GetMasterFieldNumFromDlgFieldNum (dlg->master_field_list, dlg->field_list, dlg->last_viewed_col);
  if (dlg->editor_list[edit_col] == NULL) {
    return;
  }
  str = ClipboardToString ();
  if (!StringHasNoText (str)) {
    PointerToDialog (dlg->editor_list[edit_col], str);
  }
  str = MemFree (str);
}


static BulkEdFieldPtr GetPopulatedFieldsFromMasterList (BulkEdFieldPtr master_field_list, BulkEditorRowPtr row_list)
{
  Int4 i, num_fields = 0;
  BoolPtr empties;
  Int4 num_full = 0;
  BulkEdFieldPtr pop_fields;

  for (i = 0; master_field_list[i].name != NULL; i++) {
    num_fields++;
  }

  empties = (BoolPtr) MemNew (sizeof (Boolean) * num_fields);
  for (i = 0; i < num_fields; i++) {
    if (IsBulkEditorColumnEmpty (row_list, i)) {
      empties[i] = TRUE;
    } else {
      empties[i] = FALSE;
      num_full++;
    }
  }

  pop_fields = (BulkEdFieldPtr) MemNew (sizeof (BulkEdFieldData) * (num_full + 1));
  num_full = 0;
  for (i = 0; i < num_fields; i++) {
    if (!empties[i]) {
      pop_fields[num_full++] = master_field_list[i];
    }
  }
  pop_fields[num_full] = master_field_list[num_fields];
  empties = MemFree (empties);
  return pop_fields;
}


static BulkEdFieldPtr CopyAllFieldsFromList (BulkEdFieldPtr master_field_list)
{
  Int4 i, num_fields = 0;
  Int4 num_full = 0;
  BulkEdFieldPtr pop_fields;

  for (i = 0; master_field_list[i].name != NULL; i++) {
    num_fields++;
  }


  pop_fields = (BulkEdFieldPtr) MemNew (sizeof (BulkEdFieldData) * (num_fields + 1));
  for (i = 0; i < num_fields; i++) {
    pop_fields[i] = master_field_list[i];
  }
  pop_fields[num_fields] = master_field_list[num_fields];
  return pop_fields;
}

static void FinishBulkEditorDialogSetup (GrouP p, BulkEditorDlgPtr dlg)
{
  Int4             width, num_columns = 0;
  Int4             i, page_num;
  GrouP            ed_grp, ed2_grp, check_grp, action_grp, page_grp, box_grp;
  GrouP            g1, g2, a, g_tmp, g_tmp2;
  ButtoN           b; 

  if (p == NULL || dlg == NULL) {
    return;
  }

  dlg->dialog = (DialoG) p;
  dlg->todialog = ObjectListToBulkEditorDialog; 
  dlg->fromdialog = GetSelectedObjectsFromBulkEditorDialog;
  dlg->testdialog = NULL;

  /* look for editing functions */
  dlg->num_master_columns = 0;
  dlg->any_editing = FALSE;
  for (i = 0; dlg->master_field_list[i].name != NULL; i++) {
    if (dlg->master_field_list[i].create_dlg_func != NULL) {
      dlg->any_editing = TRUE;
    }
    dlg->num_master_columns++;
  }
  /* count columns */
  for (i = 0; dlg->field_list[i].name != NULL; i++) {
    num_columns++;
  }

  dlg->num_columns = num_columns;
  dlg->columns_shown = num_columns;
  dlg->display_offset = 0;
  /* populate values list */

  SelectFont (systemFont);
  dlg->col_list = AllocateDocColumnsForFieldList (dlg->field_list, &width, dlg->last_sorted_col, dlg->display_offset, dlg->columns_shown);
  dlg->par_fmt.minLines = 1;
  dlg->par_fmt.minHeight = 1;

  dlg->title_doc = DocumentPanel (p, width, stdLineHeight);
  InitBulkEdTitleDoc (dlg, dlg->title_doc);
  SetObjectExtra (dlg->title_doc, dlg, NULL);
  SetDocProcs (dlg->title_doc, ClickBulkTitle, NULL, NULL, NULL);

  dlg->doc = DocumentPanel (p, width, 10 * stdLineHeight);
  SetObjectExtra (dlg->doc, dlg, NULL);
  SetDocProcs (dlg->doc, ClickBulkEdField, NULL, NULL, NULL);
  SetDocShade (dlg->doc, DrawBulkEdCell, NULL, HighlightBulkEdField, NULL);

  ed_grp = HiddenGroup (p, 4, 0, NULL);
  SetGroupSpacing (ed_grp, 10, 10);
  ed2_grp = HiddenGroup (ed_grp, 0, 0, NULL);
  dlg->editor_list = (DialoG PNTR) MemNew (sizeof(DialoG) * dlg->num_master_columns);
  for (i = 0; dlg->master_field_list[i].name != NULL; i++) {
    if (dlg->master_field_list[i].create_dlg_func == NULL) {
      dlg->editor_list[i] = NULL;
    } else {
      dlg->editor_list[i] = (dlg->master_field_list[i].create_dlg_func) (ed2_grp, dlg->master_field_list[i].name, dlg->sep);
      Hide (dlg->editor_list[i]);
    }
  }
  
  dlg->ed_btn_grp = HiddenGroup (ed_grp, 1, 0, NULL);
  SetGroupSpacing (dlg->ed_btn_grp, 10, 10);
  
  b = PushButton (dlg->ed_btn_grp, "Apply Change", ApplyBulkEditingChange);
  SetObjectExtra (b, dlg, NULL);
  g_tmp = HiddenGroup (dlg->ed_btn_grp, 2, 0, NULL);
  SetGroupSpacing (g_tmp, 10, 10);
  b = PushButton (g_tmp, "Copy to Clipboard", CopyEditToClipboard);
  SetObjectExtra (b, dlg, NULL);
  b = PushButton (g_tmp, "Copy from Clipboard", CopyEditFromClipboard);
  SetObjectExtra (b, dlg, NULL);
  Hide (dlg->ed_btn_grp);

  box_grp = HiddenGroup (p, 0, 2, NULL);

  /* group that contains controls for checking and unchecking features */
  check_grp = NormalGroup (box_grp, -1, 0, "Check Features", programFont, NULL);
  SetGroupSpacing (check_grp, 10, 10);
  g1 = HiddenGroup (check_grp, 3, 0, NULL);
  SetGroupSpacing (g1, 10, 10);
  b = PushButton (g1, "Check All Features", BulkEditorCheckAll);
  SetObjectExtra (b, dlg, NULL);
  b = PushButton (g1, "Uncheck All Features", BulkEditorUncheckAll);
  SetObjectExtra (b, dlg, NULL);

  g2 = HiddenGroup (check_grp, 4, 0, NULL);
  dlg->sel_unsel_grp = HiddenGroup (g2, 0, 2, NULL);
  RadioButton (dlg->sel_unsel_grp, "Check features where");
  RadioButton (dlg->sel_unsel_grp, "Uncheck features where");
  SetValue (dlg->sel_unsel_grp, 1);
  dlg->check_constraint_field = MakeBulkEditorFieldListPopup (g2, dlg->field_list, FALSE);
  dlg->check_constraint = StringConstraintDialogX (g2, NULL, FALSE);
  b = PushButton (g2, "Apply", BulkEditorStringConstraintCheck);
  SetObjectExtra (b, dlg, NULL);

  if (dlg->any_editing) {
    /* group that contains controls for applying actions to check features */
    action_grp = NormalGroup (box_grp, -1, 0, "Action for Checked Features", programFont, NULL);
    a = HiddenGroup (action_grp, 2, 0, NULL);
    SetGroupSpacing (a, 10, 10);
    dlg->bulk_action = PopupList (a, TRUE, ChangeBulkEditorAction);
    SetObjectExtra (dlg->bulk_action, dlg, NULL);
    PopupItem (dlg->bulk_action, "Apply");
    PopupItem (dlg->bulk_action, "Edit");
    PopupItem (dlg->bulk_action, "Convert");
    PopupItem (dlg->bulk_action, "Parse");
    PopupItem (dlg->bulk_action, "Swap");
    PopupItem (dlg->bulk_action, "Remove");
    PopupItem (dlg->bulk_action, "Change Case");
    SetValue (dlg->bulk_action, 2);
    page_grp = HiddenGroup (a, 0, 0, NULL);
    page_num = 0;
    dlg->action_pages[page_num] = HiddenGroup (page_grp, 2, 0, NULL);
    dlg->bulk_apply = OneFieldBulkEditDialog (dlg->action_pages[page_num], eBulkApplyField, NULL, NULL, dlg->master_field_list, GetFirstBulkEditValue, dlg);

    page_num++;
    dlg->action_pages[page_num] = HiddenGroup (page_grp, 2, 0, NULL);
    dlg->bulk_edit = OneFieldBulkEditDialog (dlg->action_pages[page_num], eBulkEditField, NULL, NULL, dlg->master_field_list, GetFirstBulkEditValue, dlg);

    page_num++;
    dlg->action_pages[page_num] = HiddenGroup (page_grp, 2, 0, NULL);
    dlg->bulk_convert = TwoFieldBulkEditDialog (dlg->action_pages[page_num], eBulkConvertField, NULL, NULL, dlg->master_field_list, FALSE, NULL, NULL);

    page_num++;
    dlg->action_pages[page_num] = HiddenGroup (page_grp, 2, 0, NULL);
    dlg->bulk_parse = TwoFieldBulkEditDialog (dlg->action_pages[page_num], eBulkParseField, NULL, NULL, dlg->master_field_list, FALSE, GetFirstBulkEditValue, dlg);

    page_num++;
    dlg->action_pages[page_num] = HiddenGroup (page_grp, 2, 0, NULL);
    dlg->bulk_swap = TwoFieldBulkEditDialog (dlg->action_pages[page_num], eBulkSwapField, NULL, NULL, dlg->master_field_list, FALSE, NULL, NULL);

    page_num++;
    dlg->action_pages[page_num] = HiddenGroup (page_grp, 2, 0, NULL);  
    dlg->bulk_remove = OneFieldBulkEditDialog (dlg->action_pages[page_num], eBulkRemoveField, NULL, NULL, dlg->master_field_list, NULL, NULL);

    page_num++;
    dlg->action_pages[page_num] = HiddenGroup (page_grp, 2, 0, NULL);  
    dlg->bulk_capitalize = OneFieldBulkEditDialog (dlg->action_pages[page_num], eBulkCapitalizeField, NULL, NULL, dlg->master_field_list, NULL, NULL);

    AlignObjects (ALIGN_CENTER, (HANDLE) dlg->action_pages[0], 
                                (HANDLE) dlg->action_pages[1], 
                                (HANDLE) dlg->action_pages[2], 
                                (HANDLE) dlg->action_pages[3], 
                                (HANDLE) dlg->action_pages[4], 
                                (HANDLE) dlg->action_pages[5], 
                                (HANDLE) dlg->action_pages[6], 
                                NULL);

    dlg->check_status = StaticPrompt (action_grp, "0 features are checked", 20 * stdCharWidth, dialogTextHeight, programFont, 'l');

    g_tmp2 = HiddenGroup (action_grp, 3, 0, NULL);
    SetGroupSpacing (g_tmp2, 10, 10);
    b = PushButton (g_tmp2, "Apply to All Features", ApplyBulkEditAll);
    SetObjectExtra (b, dlg, NULL);
    b = PushButton (g_tmp2, "Apply to Checked Features", ApplyBulkEdit);
    SetObjectExtra (b, dlg, NULL);

    AlignObjects (ALIGN_CENTER, (HANDLE) a, (HANDLE) dlg->check_status, (HANDLE) g_tmp2, NULL);
  }


  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->title_doc, (HANDLE) dlg->doc, (HANDLE) ed_grp, (HANDLE) box_grp, NULL);

  UpdateBulkEdDisplay (dlg);
  UpdateCheckStatus (dlg);
  ChangeBulkEditorAction (dlg->bulk_action);
}


NLM_EXTERN DialoG 
CreateBulkEditorDialogEx 
(GrouP               h,
 BulkEdFieldPtr      field_list,
 ValNodePtr          feat_list,
 SeqEntryPtr         sep, 
 Boolean             collapse_by_default,
 ClickableCallback   single_click_func,
 ClickableCallback   double_click_func,
 Pointer             metadata,
 Boolean             show_only_present)
{
  GrouP            p;
  BulkEditorDlgPtr dlg;

  p = HiddenGroup (h, -1, 0, NULL);
  dlg = (BulkEditorDlgPtr) MemNew (sizeof (BulkEditorDlgData));

  SetObjectExtra (p, dlg, CleanupBulkEditorDialog); /* need to cleanup values_list */

  dlg->collapse_by_default = collapse_by_default;
  dlg->single_click_func = single_click_func;
  dlg->double_click_func = double_click_func;
  dlg->master_field_list = field_list;
  dlg->metadata = metadata;
  dlg->sep = sep;

  dlg->last_viewed_row = 0;
  dlg->last_viewed_col = 0;
  dlg->last_sorted_col = 0;
  dlg->row_list = MakeBulkEditorRowsFromObjectList(feat_list, dlg->master_field_list, dlg->last_sorted_col, dlg->metadata);

  /* need master list and actual list */
  if (show_only_present) {
    dlg->field_list = GetPopulatedFieldsFromMasterList (dlg->master_field_list, dlg->row_list);
    dlg->row_list = FreeBulkEditorRows(dlg->row_list, dlg->master_field_list);
    dlg->row_list = MakeBulkEditorRowsFromObjectList(feat_list, dlg->field_list, dlg->last_sorted_col, dlg->metadata);
  } else {
    dlg->field_list = CopyAllFieldsFromList (dlg->master_field_list);
  }

  if (dlg->collapse_by_default) {
    CollapseAllBulkEditRows (dlg->row_list);
  } else {
    ExpandAllBulkEditRows (dlg->row_list);
  }

  FinishBulkEditorDialogSetup (p, dlg);

  return (DialoG) p;
}


NLM_EXTERN DialoG 
CreateBulkEditorDialog 
(GrouP               h,
 BulkEdFieldPtr      field_list,
 ValNodePtr          feat_list,
 SeqEntryPtr         sep, 
 Boolean             collapse_by_default,
 ClickableCallback   single_click_func,
 ClickableCallback   double_click_func)
{
  return CreateBulkEditorDialogEx (h, field_list, feat_list, sep, collapse_by_default,
                                   single_click_func, double_click_func, NULL, FALSE);
}


static Pointer BulkGetLocation (Uint1 data_choice, Pointer data, Pointer metadata)
{
  SeqLocPtr slp;
  SeqFeatPtr sfp = (SeqFeatPtr) data;

  if (sfp == NULL) return NULL;
  slp = (SeqLocPtr) AsnIoMemCopy (sfp->location, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
  return slp;
}

static void BulkSetLocation (Pointer target, Pointer data)
{
  SeqFeatPtr sfp = (SeqFeatPtr) target;
  SeqLocPtr slp = (SeqLocPtr) data;
  if (sfp == NULL) return;
  sfp->location = SeqLocFree (sfp->location);
  sfp->location = (SeqLocPtr) AsnIoMemCopy (slp, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
}

static Int4 BulkFormatLocation (ColPtr col, CharPtr name)
{
  if (col == NULL) return 0;

  col->pixWidth = MAX (15, StringLen (name)) * stdCharWidth;
  col->pixInset = 1;
  col->charWidth = 0;
  col->charInset = 0;
  col->font = NULL;
  col->just = 'l';
  col->wrap = 1;
  col->bar = 0;
  col->underline = 0;
  col->left = 0;
  return col->pixWidth;
}

static void BulkFreeLocation (Pointer data)
{
  SeqLocPtr slp = (SeqLocPtr) data;

  slp = SeqLocFree (slp);
}

static CharPtr BulkDisplayLocation (Pointer data)
{
  SeqLocPtr slp = (SeqLocPtr) data;

  return SeqLocPrint (slp);
}

static DialoG BulkNucLocationDialog (GrouP g, CharPtr name, SeqEntryPtr sep)
{
  return CreateIntervalEditorDialog (g, NULL, 4, 2, sep, TRUE, FALSE);
}

static Pointer BulkCopyLocation (Pointer data)
{
  SeqLocPtr orig = (SeqLocPtr) data;
  SeqLocPtr new_loc;

  if (orig == NULL) return NULL;
  new_loc = (SeqLocPtr) AsnIoMemCopy ((Pointer) orig,
                          (AsnReadFunc) SeqLocAsnRead,
                          (AsnWriteFunc) SeqLocAsnWrite);

  return new_loc;
}

NLM_EXTERN Int4 BulkFormatSimpleText (ColPtr col, CharPtr name)
{
  if (col == NULL) return 0;

  col->pixWidth = (MAX (15, StringLen (name)) + 3) * MaxCharWidth();
  col->pixInset = 0;
  col->charWidth = 0;
  col->charInset = 0;
  col->font = NULL;
  col->just = 'l';
  col->wrap = 1;
  col->bar = 0;
  col->underline = 0;
  col->left = 0;
  return col->pixWidth;
}

NLM_EXTERN void BulkFreeSimpleText (Pointer data)
{
  CharPtr str = (CharPtr) data;

  str = MemFree (str);
}

NLM_EXTERN CharPtr BulkDisplaySimpleText (Pointer data)
{
  CharPtr str = (CharPtr) data;
  if (StringHasNoText (str)) {
    return NULL;
  } else {
    return StringSave ((CharPtr) data);
  }
}

typedef struct bulksimpletext {
  DIALOG_MESSAGE_BLOCK
  TexT txt;
} BulkSimpleTextData, PNTR BulkSimpleTextPtr;

static void StringToBulkSimpleTextDialog (DialoG d, Pointer data)
{
  BulkSimpleTextPtr dlg;
  CharPtr           str;

  dlg = (BulkSimpleTextPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  str = (CharPtr) data;
  if (str == NULL) {
    SetTitle (dlg->txt, "");
  } else {
    SetTitle (dlg->txt, str);
  }
}

static Pointer BulkSimpleTextDialogToString (DialoG d)
{
  BulkSimpleTextPtr dlg;

  dlg = (BulkSimpleTextPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  return (Pointer) SaveStringFromText (dlg->txt);
}


static DialoG BulkSimpleTextDialogEx (GrouP g, CharPtr name, SeqEntryPtr sep, Boolean use_scroll)
{
  BulkSimpleTextPtr dlg;
  GrouP           p;

  p = HiddenGroup (g, 4, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  dlg = (BulkSimpleTextPtr) MemNew (sizeof(BulkSimpleTextData));

  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = StringToBulkSimpleTextDialog;
  dlg->fromdialog = BulkSimpleTextDialogToString;
  dlg->testdialog = NULL;

  StaticPrompt (p, name, 0, dialogTextHeight, programFont, 'l');
  if (use_scroll) {
    dlg->txt = ScrollText (p, 20, 4, programFont, TRUE, NULL);
  } else {
    dlg->txt = DialogText (p, "", 20, NULL);
  }

  return (DialoG) p;
}

NLM_EXTERN DialoG BulkSimpleTextDialog (GrouP g, CharPtr name, SeqEntryPtr sep)
{
  return BulkSimpleTextDialogEx (g, name, sep, FALSE);
}

static DialoG BulkSimpleTextScrollDialog (GrouP g, CharPtr name, SeqEntryPtr sep)
{
  return BulkSimpleTextDialogEx (g, name, sep, TRUE);
}

NLM_EXTERN Pointer BulkSimpleTextCopy (Pointer data)
{
  CharPtr orig = (CharPtr) data;
  CharPtr new_str = StringSave (orig);
  return new_str;
}

static Int4 BulkFormatTrueFalse (ColPtr col, CharPtr name)
{
  size_t box_width = stdLineHeight - 3;
  if (col == NULL) return 0;

  col->pixWidth = MAX (box_width, (StringLen (name) + 2) * MaxCharWidth());
  col->pixInset = 0;
  col->charWidth = 0;
  col->charInset = 0;
  col->font = NULL;
  col->just = 'l';
  col->wrap = 0;
  col->bar = 0;
  col->underline = 0;
  col->left = 0;
  return col->pixWidth;
}

static void BulkDrawTrueFalse (Pointer p, RectPtr r)
{
  RecT rct;

  if (r == NULL) return;

  rct = *r;
  rct.bottom = rct.top + stdLineHeight - 4;
  rct.right = rct.left + stdLineHeight - 4;
  Gray();
  InsetRect (&rct, 2, 2);
  FrameRect (&rct);
  
  if (p != NULL) {
    MoveTo (rct.left, rct.top);
    LineTo (rct.right - 1, rct.bottom - 1);
    MoveTo (rct.left, rct.bottom - 1);
    LineTo (rct.right - 1, rct.top);
  }
  Black();
}

static Pointer BulkReleaseTrueFalse (Pointer p)
{
  if (p == NULL) {
    p = (Pointer) StringSave ("TRUE");
  } else {
    p = MemFree (p);
  }
  return p;
}


static void BulkSetPseudo (Pointer target, Pointer val)
{
  SeqFeatPtr sfp = (SeqFeatPtr) target;
  if (sfp == NULL) return;
  if (StringHasNoText ((CharPtr)val)) {
    sfp->pseudo = FALSE;
  } else {
    sfp->pseudo = TRUE;
  }
}

static Pointer BulkGetPseudo (Pointer data)
{
  SeqFeatPtr sfp = (SeqFeatPtr) data;
  if (sfp == NULL || ! sfp->pseudo) {
    return NULL;
  } else {
    return StringSave ("TRUE");
  }
}

static void BulkSetComment (Pointer target, Pointer data)
{
  CharPtr str = (CharPtr) data;
  SeqFeatPtr sfp = (SeqFeatPtr) target;
  if (sfp == NULL) return;

  sfp->comment = MemFree (sfp->comment);
  if (str != NULL) {
    sfp->comment = StringSave (str);
  }
}

static Pointer BulkGetComment (Uint1 data_choice, Pointer data, Pointer metadata)
{
  SeqFeatPtr sfp = (SeqFeatPtr) data;
  if (sfp == NULL || StringHasNoText (sfp->comment)) {
    return NULL;
  } else {
    return (Pointer) StringSave (sfp->comment);
  }
}


static Pointer BulkCDSGetProtein (Uint1 data_choice, Pointer data, Pointer metadata)
{
  SeqFeatXrefPtr xref;
  ProtRefPtr     prp = NULL;
  BioseqPtr      prot_bsp;
  SeqFeatPtr     prot_feat;
  SeqMgrFeatContext fcontext;
  ValNodePtr        prot_name_list = NULL, vnp;
  SeqFeatPtr        sfp = (SeqFeatPtr) data;

  if (sfp == NULL) return NULL;

  xref = sfp->xref;
  while (xref != NULL && xref->data.choice != SEQFEAT_PROT) {
    xref = xref->next;
  }
  if (xref != NULL) {
    prp = xref->data.value.ptrvalue;
  }

  if (prp == NULL) {
    if (sfp->data.choice == SEQFEAT_PROT) {
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    } else if (sfp->data.choice == SEQFEAT_CDREGION) {
      prot_bsp = BioseqFindFromSeqLoc (sfp->product);
      prot_feat = SeqMgrGetNextFeature (prot_bsp, NULL, SEQFEAT_PROT, FEATDEF_PROT, &fcontext);
      if (prot_feat != NULL) {
        prp = (ProtRefPtr) prot_feat->data.value.ptrvalue;
      }
    }
  }

  if (prp != NULL) {
    for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
      ValNodeAddPointer (&prot_name_list, vnp->choice, StringSave (vnp->data.ptrvalue));
    }
  }
    
  return prot_name_list;
}

static Pointer BulkCDSCopyProtein (Pointer data)
{
  ValNodePtr        new_prot_name_list = NULL, prot_name_list, vnp;
 
  prot_name_list = (ValNodePtr) data;
  
  for (vnp = prot_name_list; vnp != NULL; vnp = vnp->next) {
    ValNodeAddPointer (&new_prot_name_list, vnp->choice, StringSave (vnp->data.ptrvalue));
  }
  return new_prot_name_list;
}  

static Int4 BulkCDSFormatProtein (ColPtr col, CharPtr name)
{
  if (col == NULL) return 0;

  col->pixWidth = MAX (15, StringLen (name)) * MaxCharWidth();
  col->pixInset = 0;
  col->charWidth = 0;
  col->charInset = 0;
  col->font = NULL;
  col->just = 'l';
  col->wrap = 1;
  col->bar = 0;
  col->underline = 0;
  col->left = 0;
  return col->pixWidth;
}

static CharPtr BulkCDSDisplayProtein (Pointer data)
{
  ValNodePtr prot_name_list, vnp;
  CharPtr    prot_text;
  Int4       text_len = 1;

  if (data == NULL) return NULL;
  prot_name_list = (ValNodePtr) data;
  for (vnp = prot_name_list; vnp != NULL; vnp = vnp->next) {
    text_len += StringLen (vnp->data.ptrvalue) + 1;
  }
  prot_text = (CharPtr) MemNew (sizeof (Char) * text_len);
  prot_text[0] = 0;
  for (vnp = prot_name_list; vnp != NULL; vnp = vnp->next) {
    StringCat (prot_text, vnp->data.ptrvalue);
    if (vnp->next != NULL) {
      StringCat (prot_text, ";");
    }
  }
  return prot_text;
}


static Pointer BulkCDSSetProteinString (Pointer curr_val, ApplyValuePtr avp)
{
  CharPtr        curr_val_str;
  CharPtr        new_str, name_str, cp;
  ValNodePtr     name_list;

  name_list = (ValNodePtr) curr_val;
  curr_val_str = BulkCDSDisplayProtein (name_list);
  name_list = ValNodeFreeData (name_list);

  new_str = HandleApplyValue (curr_val_str, avp);

  name_str = new_str;
  cp = StringChr (name_str, ';');
  while (cp != NULL) {
    *cp = 0;
    ValNodeAddPointer (&name_list, 0, StringSave (name_str));
    *cp = ';';
    name_str = cp + 1;
    cp = StringChr (name_str, ';');
  }
  ValNodeAddPointer (&name_list, 0, StringSave (name_str));
  new_str = MemFree (new_str);

  return name_list;
}


static void BulkCDSFreeProtein (Pointer data)
{
  ValNodePtr prot_name_list = (ValNodePtr) data;
 
  ValNodeFreeData (prot_name_list);
}

static DialoG BulkCDSProteinDialog (GrouP g, CharPtr name, SeqEntryPtr sep)
{
  return CreateVisibleStringDialog (g, 3, -1, 25);
}

static void BulkCDSSetProteinDesc (Pointer target, Pointer data)
{
  SeqFeatXrefPtr xref;
  ProtRefPtr     prp = NULL;
  BioseqPtr      prot_bsp;
  SeqFeatPtr     prot_feat;
  SeqMgrFeatContext fcontext;
  CharPtr           prot_desc;
  SeqFeatPtr sfp = (SeqFeatPtr) target;
  if (sfp == NULL) return;

  prot_desc = (CharPtr) data;

  xref = sfp->xref;
  while (xref != NULL && xref->data.choice != SEQFEAT_PROT) {
    xref = xref->next;
  }
  if (xref != NULL) {
    prp = xref->data.value.ptrvalue;
  }

  if (prp == NULL) {
    if (sfp->data.choice == SEQFEAT_PROT) {
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    } else if (sfp->data.choice == SEQFEAT_CDREGION) {
      prot_bsp = BioseqFindFromSeqLoc (sfp->product);
      prot_feat = SeqMgrGetNextFeature (prot_bsp, NULL, SEQFEAT_PROT, FEATDEF_PROT, &fcontext);
      if (prot_feat == NULL && !StringHasNoText (prot_desc)) {
        prot_feat = CreateNewFeatureOnBioseq (prot_bsp, SEQFEAT_PROT, NULL);
        prp = ProtRefNew ();
        prot_feat->data.value.ptrvalue = prp;
          SeqMgrIndexFeatures (prot_bsp->idx.entityID, NULL);
      }
      if (prot_feat != NULL) {
        prp = (ProtRefPtr) prot_feat->data.value.ptrvalue;
      }
    }
  }

  if (prp == NULL && !StringHasNoText (prot_desc)) {
    xref = SeqFeatXrefNew ();
    prp = ProtRefNew ();
    xref->data.choice = SEQFEAT_PROT;
    xref->data.value.ptrvalue = prp;
    xref->next = sfp->xref;
    sfp->xref = xref;
  }
  if (prp != NULL) {
    prp->desc = StringSave (prot_desc);
  }

}


static Pointer BulkCDSGetProteinDesc (Uint1 data_choice, Pointer data, Pointer metadata)
{
  SeqFeatXrefPtr xref;
  ProtRefPtr     prp = NULL;
  BioseqPtr      prot_bsp;
  SeqFeatPtr     prot_feat;
  SeqMgrFeatContext fcontext;
  CharPtr           prot_desc = NULL;
  SeqFeatPtr sfp = (SeqFeatPtr) data;

  if (sfp == NULL) return NULL;

  xref = sfp->xref;
  while (xref != NULL && xref->data.choice != SEQFEAT_PROT) {
    xref = xref->next;
  }
  if (xref != NULL) {
    prp = xref->data.value.ptrvalue;
  }

  if (prp == NULL) {
    if (sfp->data.choice == SEQFEAT_PROT) {
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    } else if (sfp->data.choice == SEQFEAT_CDREGION) {
      prot_bsp = BioseqFindFromSeqLoc (sfp->product);
      prot_feat = SeqMgrGetNextFeature (prot_bsp, NULL, SEQFEAT_PROT, FEATDEF_PROT, &fcontext);
      if (prot_feat != NULL) {
        prp = (ProtRefPtr) prot_feat->data.value.ptrvalue;
      }
    }
  }

  if (prp != NULL) {
    prot_desc = StringSave (prp->desc);
  }
    
  return prot_desc;
}


static void BulkGeneSetLocus (Pointer target, Pointer data)
{
  CharPtr    str = (CharPtr) data;
  GeneRefPtr grp;
  SeqFeatPtr sfp = (SeqFeatPtr) target;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_GENE) return;

  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) {
    grp = GeneRefNew();
    sfp->data.value.ptrvalue = grp;
  }
  grp->locus = MemFree (grp->locus);
  if (str != NULL) {
    grp->locus = StringSave (str);
  }
}

static Pointer BulkGeneGetLocus (Uint1 data_choice, Pointer data, Pointer metadata)
{
  GeneRefPtr grp; 
  SeqFeatPtr sfp = (SeqFeatPtr) data;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_GENE || sfp->data.value.ptrvalue == NULL) {
    return NULL;
  }
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (StringHasNoText (grp->locus)) {
    return NULL;
  } else {
    return (Pointer) StringSave (grp->locus);
  }
}

static void BulkGeneSetLocusTag (Pointer target, Pointer data)
{
  CharPtr    str = (CharPtr) data;
  GeneRefPtr grp;
  SeqFeatPtr sfp = (SeqFeatPtr) target;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_GENE) return;

  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) {
    grp = GeneRefNew();
    sfp->data.value.ptrvalue = grp;
  }
  grp->locus_tag = MemFree (grp->locus_tag);
  if (str != NULL) {
    grp->locus_tag = StringSave (str);
  }
}

static Pointer BulkGeneGetLocusTag (Uint1 data_choice, Pointer data, Pointer metadata)
{
  GeneRefPtr grp;
  SeqFeatPtr sfp = (SeqFeatPtr) data;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_GENE || sfp->data.value.ptrvalue == NULL) {
    return NULL;
  }
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (StringHasNoText (grp->locus_tag)) {
    return NULL;
  } else {
    return (Pointer) StringSave (grp->locus_tag);
  }
}

static void BulkGeneSetDescription (Pointer target, Pointer data)
{
  CharPtr    str = (CharPtr) data;
  GeneRefPtr grp;
  SeqFeatPtr sfp = (SeqFeatPtr) target;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_GENE) return;

  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) {
    grp = GeneRefNew();
    sfp->data.value.ptrvalue = grp;
  }
  grp->desc = MemFree (grp->desc);
  if (str != NULL) {
    grp->desc = StringSave (str);
  }
}

static Pointer BulkGeneGetDescription (Uint1 data_choice, Pointer data, Pointer metadata)
{
  GeneRefPtr grp;
  SeqFeatPtr sfp = (SeqFeatPtr) data;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_GENE || sfp->data.value.ptrvalue == NULL) {
    return NULL;
  }
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (StringHasNoText (grp->desc)) {
    return NULL;
  } else {
    return (Pointer) StringSave (grp->desc);
  }
}


static void BulkRNASetProduct (Pointer target, Pointer data)
{
  RnaRefPtr  rrp;
  CharPtr    rna_product;
  SeqFeatPtr sfp = (SeqFeatPtr) target;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return;

  rna_product = (CharPtr) data;

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL && !StringHasNoText (rna_product)) {
    rrp = RnaRefNew ();
    rrp->ext.choice = 1;
    sfp->data.value.ptrvalue = rrp;
  }
  if (rrp->ext.choice == 0 && !StringHasNoText (rna_product)) {
    rrp->ext.choice = 1;
  } else if (rrp->ext.choice == 1 && StringHasNoText (rna_product)) {
    rrp->ext.choice = 0;
    rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
  }
 
  if (rrp != NULL && rrp->ext.choice == 1) {    
    rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
    rrp->ext.value.ptrvalue = StringSave (rna_product);
  }
}

static Pointer BulkRNAGetProduct (Uint1 data_choice, Pointer data, Pointer metadata)
{
  RnaRefPtr  rrp;
  CharPtr    rna_product = NULL;
  SeqFeatPtr sfp = (SeqFeatPtr) data;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return NULL;

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp != NULL) {
    if (rrp->ext.choice == 1) {
      rna_product = StringSave( rrp->ext.value.ptrvalue);
    }
  }

    
  return rna_product;
}

static Boolean IsBulkEditableRNA (SeqFeatPtr sfp)
{
  RnaRefPtr rrp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) {
    return FALSE;
  }
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp != NULL && rrp->ext.choice != 1 && rrp->ext.choice != 0) {
    return FALSE;
  } else {
    return TRUE;
  }
}


/* template for setting simple text strings */
NLM_EXTERN Pointer BulkSetSimpleTextString (Pointer curr_val, ApplyValuePtr avp)
{
  CharPtr        curr_val_str;

  curr_val_str = (CharPtr)curr_val;
  return HandleApplyValue (curr_val_str, avp);
}

/* for now, disallow pseudo editing and display.  To add back, put in this line:
  { "pseudo", BulkSetPseudo, BulkGetPseudo, NULL, NULL, NULL, BulkFormatTrueFalse, BulkDrawTrueFalse, BulkReleaseTrueFalse },
*/
/* for now, disallow location editing (but continue display).  To add back, use BulkNucLocationDialog in dialog. */
static BulkEdFieldData bulk_cds_fields[] = {
  { "protein name                     ", BulkCDSSetProtein, BulkCDSSetProteinString, BulkCDSGetProtein, BulkCDSDisplayProtein, BulkCDSFreeProtein, BulkCDSProteinDialog, BulkCDSFormatProtein, NULL, NULL, BulkCDSCopyProtein },
  { "protein description", BulkCDSSetProteinDesc, BulkSetSimpleTextString, BulkCDSGetProteinDesc, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "location", BulkSetLocation, NULL, BulkGetLocation, BulkDisplayLocation, BulkFreeLocation, NULL, BulkFormatLocation, NULL, NULL, BulkCopyLocation },
  { "comment", BulkSetComment, BulkSetSimpleTextString, BulkGetComment, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextScrollDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

static BulkEdFieldData bulk_gene_fields[] = {
  { "locus", BulkGeneSetLocus, BulkSetSimpleTextString, BulkGeneGetLocus, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "locus_tag", BulkGeneSetLocusTag, BulkSetSimpleTextString, BulkGeneGetLocusTag, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "description", BulkGeneSetDescription, BulkSetSimpleTextString, BulkGeneGetDescription, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "location", BulkSetLocation, NULL, BulkGetLocation, BulkDisplayLocation, BulkFreeLocation, NULL, BulkFormatLocation, NULL, NULL, BulkCopyLocation },
  { "comment", BulkSetComment, BulkSetSimpleTextString, BulkGetComment, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextScrollDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

static BulkEdFieldData bulk_rna_fields[] = {
  { "product", BulkRNASetProduct, BulkSetSimpleTextString, BulkRNAGetProduct, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "location", BulkSetLocation, NULL, BulkGetLocation, BulkDisplayLocation, BulkFreeLocation, NULL, BulkFormatLocation, NULL, NULL, BulkCopyLocation },
  { "comment", BulkSetComment, BulkSetSimpleTextString, BulkGetComment, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextScrollDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};


static BulkEdFieldPtr GetBulkEditorFieldListForFeature (SeqFeatPtr sfp)
{
  if (sfp == NULL) {
    return NULL;
  } else if (sfp->data.choice == SEQFEAT_CDREGION) {
    return bulk_cds_fields;
  } else if (sfp->data.choice == SEQFEAT_GENE) {
    return bulk_gene_fields;
  } else if (sfp->data.choice == SEQFEAT_RNA && IsBulkEditableRNA(sfp)) {
    /* Note - using FEATDEF_rRNA to represent all editable rRNA features */
    return bulk_rna_fields;
  } else {
    return NULL;
  }
}


static void BulkSrcSetTaxNameDesc (Pointer target, Pointer data)
{
  CharPtr      taxname;
  SeqDescrPtr  sdp = (SeqDescrPtr) target;
  BioSourcePtr biop;

  if (sdp == NULL || sdp->choice != Seq_descr_source || sdp->data.ptrvalue == NULL) return;

  biop = (BioSourcePtr) sdp->data.ptrvalue;

  taxname = (CharPtr) data;

  if (biop->org == NULL) {
    biop->org = OrgRefNew();
  }
  SetTaxNameAndRemoveTaxRef (biop->org, StringSave (taxname)); 
}


NLM_EXTERN BioSourcePtr GetBioSourceFromObject (Uint1 choice, Pointer data);

static Pointer BulkSrcGetTaxNameDesc (Uint1 data_choice, Pointer data, Pointer metadata)
{
  CharPtr      taxname = NULL;
  BioSourcePtr biop = NULL;

  biop = GetBioSourceFromObject (data_choice, data);
  if (biop == NULL) return NULL;

  if (biop->org != NULL) {
    taxname = StringSave (biop->org->taxname);
  }

  return taxname;
}


static void BulkSrcSetTaxNameFeat (Pointer target, Pointer data)
{
  CharPtr      taxname;
  SeqFeatPtr   sfp = (SeqFeatPtr) target;
  BioSourcePtr biop;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || sfp->data.value.ptrvalue == NULL) return;

  biop = (BioSourcePtr) sfp->data.value.ptrvalue;

  taxname = (CharPtr) data;

  if (biop->org == NULL) {
    biop->org = OrgRefNew();
  }
  SetTaxNameAndRemoveTaxRef (biop->org, StringSave (taxname)); 
}


static Pointer BulkSrcGetTaxNameFeat (Uint1 data_choice, Pointer data, Pointer metadata)
{
  CharPtr      taxname = NULL;
  SeqFeatPtr   sfp = (SeqFeatPtr) data;
  BioSourcePtr biop;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || sfp->data.value.ptrvalue == NULL) return NULL;

  biop = (BioSourcePtr) sfp->data.value.ptrvalue;

  if (biop->org != NULL) {
    taxname = StringSave (biop->org->taxname);
  }

  return taxname;
}


static Pointer BulkSrcGetSubSrcModDesc (Uint1 data_choice, Pointer data, Uint2 mod_type)
{
  CharPtr      val = NULL;
  BioSourcePtr biop = NULL;
  SubSourcePtr ssp;

  biop = GetBioSourceFromObject (data_choice, data);
  if (biop == NULL) return NULL;

  ssp = biop->subtype;
  while (ssp != NULL && ssp->subtype != mod_type) {
    ssp = ssp->next;
  }
  if (ssp != NULL) {
    val = StringSave (ssp->name);
  }

  return val;
}


static void BulkSrcSetSubSrcModDesc (Pointer target, Pointer data, Uint2 mod_type)
{
  CharPtr      val;
  SeqDescrPtr  sdp = (SeqDescrPtr) target;
  BioSourcePtr biop;
  SubSourcePtr ssp, ssp_last = NULL;

  if (sdp == NULL || sdp->choice != Seq_descr_source || sdp->data.ptrvalue == NULL) return;

  biop = (BioSourcePtr) sdp->data.ptrvalue;

  val = (CharPtr) data;

  ssp = biop->subtype;
  while (ssp != NULL && ssp->subtype != mod_type) {
    ssp_last = ssp;
    ssp = ssp->next;
  }
  if (ssp == NULL) {
    ssp = SubSourceNew();
    ssp->subtype = mod_type;
    if (ssp_last == NULL) {
      biop->subtype = ssp;
    } else {
      ssp_last->next = ssp;
    }
  }
  ssp->name = MemFree (ssp->name);
  ssp->name = StringSave (val);
}


#define BULK_SRC_GET_AND_SET_SUBSRC(Name,Type) \
static Pointer BulkSrcGet##Name (Uint1 data_choice, Pointer data, Pointer metadata) \
{ \
  return BulkSrcGetSubSrcModDesc(data_choice, data, SUBSRC_##Type); \
} \
static void BulkSrcSet##Name (Pointer target, Pointer data) \
{ \
  BulkSrcSetSubSrcModDesc(target, data, SUBSRC_##Type); \
}


static Pointer BulkSrcGetSubSrcModDescTF (Uint1 data_choice, Pointer data, Uint2 mod_type)
{
  CharPtr      val = NULL;
  BioSourcePtr biop = NULL;
  SubSourcePtr ssp;

  biop = GetBioSourceFromObject (data_choice, data);
  if (biop == NULL) return NULL;

  ssp = biop->subtype;
  while (ssp != NULL && ssp->subtype != mod_type) {
    ssp = ssp->next;
  }
  if (ssp != NULL) {
    val = StringSave ("TRUE");
  }

  return val;
}


static void BulkSrcSetSubSrcModDescTF (Pointer target, Pointer data, Uint2 mod_type)
{
  CharPtr      val;
  SeqDescrPtr  sdp = (SeqDescrPtr) target;
  BioSourcePtr biop;
  SubSourcePtr ssp, ssp_last = NULL;

  if (sdp == NULL || sdp->choice != Seq_descr_source || sdp->data.ptrvalue == NULL) return;

  biop = (BioSourcePtr) sdp->data.ptrvalue;

  val = (CharPtr) data;

  ssp = biop->subtype;
  while (ssp != NULL && ssp->subtype != mod_type) {
    ssp_last = ssp;
    ssp = ssp->next;
  }
  if (StringHasNoText ((CharPtr)val)) {
    if (ssp != NULL) {
      if (ssp_last == NULL) {
        biop->subtype = ssp->next;
      } else {
        ssp_last->next = ssp->next;
      }
      ssp->next = NULL;
      ssp = SubSourceFree (ssp);
    }
  } else {
    if (ssp == NULL) {
      ssp = SubSourceNew();
      ssp->subtype = mod_type;
      if (ssp_last == NULL) {
        biop->subtype = ssp;
      } else {
        ssp_last->next = ssp;
      }
    }
    ssp->name = MemFree (ssp->name);
    ssp->name = StringSave (" ");
  }
}


#define BULK_SRC_GET_AND_SET_SUBSRC_TF(Name,Type) \
static Pointer BulkSrcGet##Name (Uint1 data_choice, Pointer data, Pointer metadata) \
{ \
  return BulkSrcGetSubSrcModDescTF(data_choice, data, SUBSRC_##Type); \
} \
static void BulkSrcSet##Name (Pointer target, Pointer data) \
{ \
  BulkSrcSetSubSrcModDescTF(target, data, SUBSRC_##Type); \
}


static Pointer BulkSrcGetOrgModDesc (Uint1 data_choice, Pointer data, Uint2 mod_type)
{
  CharPtr      val = NULL;
  BioSourcePtr biop = NULL;
  OrgModPtr    mod;

  biop = GetBioSourceFromObject (data_choice, data);
  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) return NULL;

  mod = biop->org->orgname->mod;
  while (mod != NULL && mod->subtype != mod_type) {
    mod = mod->next;
  }
  if (mod != NULL) {
    val = StringSave (mod->subname);
  }

  return val;
}


static void BulkSrcSetOrgModDesc (Pointer target, Pointer data, Uint2 mod_type)
{
  CharPtr      val;
  SeqDescrPtr  sdp = (SeqDescrPtr) target;
  BioSourcePtr biop;
  OrgModPtr    mod, mod_last = NULL;

  if (sdp == NULL || sdp->choice != Seq_descr_source || sdp->data.ptrvalue == NULL) return;

  biop = (BioSourcePtr) sdp->data.ptrvalue;

  val = (CharPtr) data;

  if (biop->org == NULL) {
    biop->org = OrgRefNew ();
  }
  if (biop->org->orgname == NULL) {
    biop->org->orgname = OrgNameNew ();
  }

  mod = biop->org->orgname->mod;
  while (mod != NULL && mod->subtype != mod_type) {
    mod_last = mod;
    mod = mod->next;
  }
  if (mod == NULL) {
    mod = OrgModNew();
    mod->subtype = mod_type;
    if (mod_last == NULL) {
      biop->org->orgname->mod = mod;
    } else {
      mod_last->next = mod;
    }
  }
  mod->subname = MemFree (mod->subname);
  mod->subname = StringSave (val);
}


#define BULK_SRC_GET_AND_SET_ORGMOD(Name,Type) \
static Pointer BulkSrcGet##Name (Uint1 data_choice, Pointer data, Pointer metadata) \
{ \
  return BulkSrcGetOrgModDesc(data_choice, data, ORGMOD_##Type); \
} \
static void BulkSrcSet##Name (Pointer target, Pointer data) \
{ \
  BulkSrcSetOrgModDesc(target, data, ORGMOD_##Type); \
}


static Pointer BulkSrcGetOrgModDescTF (Uint1 data_choice, Pointer data, Uint2 mod_type)
{
  CharPtr      val = NULL;
  BioSourcePtr biop = NULL;
  OrgModPtr    mod;

  biop = GetBioSourceFromObject (data_choice, data);
  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) return NULL;

  mod = biop->org->orgname->mod;
  while (mod != NULL && mod->subtype != mod_type) {
    mod = mod->next;
  }
  if (mod != NULL) {
    val = StringSave ("TRUE");
  }

  return val;
}


static void BulkSrcSetOrgModDescTF (Pointer target, Pointer data, Uint2 mod_type)
{
  CharPtr      val;
  SeqDescrPtr  sdp = (SeqDescrPtr) target;
  BioSourcePtr biop;
  OrgModPtr    mod, mod_last = NULL;

  if (sdp == NULL || sdp->choice != Seq_descr_source || sdp->data.ptrvalue == NULL) return;

  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;

  val = (CharPtr) data;

  if (biop->org == NULL) {
    biop->org = OrgRefNew ();
  }
  if (biop->org->orgname == NULL) {
    biop->org->orgname = OrgNameNew ();
  }

  mod = biop->org->orgname->mod;
  while (mod != NULL && mod->subtype != mod_type) {
    mod = mod->next;
  }

  if (StringHasNoText ((CharPtr)val)) {
    /* remove */
    if (mod != NULL) {
      if (mod_last == NULL) {
        biop->org->orgname->mod = mod->next;
      } else {
        mod_last->next = mod->next;
      }
      mod->next = NULL;
      mod = OrgModFree (mod);
    }
  } else {
    /* set */
    if (mod == NULL) {
      mod = OrgModNew();
      mod->subtype = mod_type;
      if (mod_last == NULL) {
        biop->org->orgname->mod = mod;
      } else {
        mod_last->next = mod;
      }
    }
    mod->subname = MemFree (mod->subname);
    mod->subname = StringSave (" ");
  }
}


#define BULK_SRC_GET_AND_SET_ORGMOD_TF(Name,Type) \
static Pointer BulkSrcGet##Name (Uint1 data_choice, Pointer data, Pointer metadata) \
{ \
  return BulkSrcGetOrgModDescTF(data_choice, data, SUBSRC_##Type); \
} \
static void BulkSrcSet##Name (Pointer target, Pointer data) \
{ \
  BulkSrcSetOrgModDescTF(target, data, SUBSRC_##Type); \
}


BULK_SRC_GET_AND_SET_ORGMOD(Acronym,acronym)
BULK_SRC_GET_AND_SET_ORGMOD(Anamorph,anamorph)
BULK_SRC_GET_AND_SET_ORGMOD(Authority,authority)
BULK_SRC_GET_AND_SET_ORGMOD(Bio_material,bio_material)
BULK_SRC_GET_AND_SET_ORGMOD(Biotype,biotype)
BULK_SRC_GET_AND_SET_ORGMOD(Biovar,biovar)
BULK_SRC_GET_AND_SET_ORGMOD(Breed,breed)
BULK_SRC_GET_AND_SET_ORGMOD(Chemovar,chemovar)
BULK_SRC_GET_AND_SET_ORGMOD(Common,common)
BULK_SRC_GET_AND_SET_ORGMOD(Cultivar,cultivar)
BULK_SRC_GET_AND_SET_ORGMOD(Culture_collection,culture_collection)
BULK_SRC_GET_AND_SET_ORGMOD(Ecotype,ecotype)
BULK_SRC_GET_AND_SET_ORGMOD(Forma,forma)
BULK_SRC_GET_AND_SET_ORGMOD(Forma_specialis,forma_specialis)
BULK_SRC_GET_AND_SET_ORGMOD(Group,group)
BULK_SRC_GET_AND_SET_ORGMOD(Host,nat_host)
BULK_SRC_GET_AND_SET_ORGMOD(Isolate,isolate)
BULK_SRC_GET_AND_SET_ORGMOD(Metagenome_source,metagenome_source)
BULK_SRC_GET_AND_SET_ORGMOD(Pathovar,pathovar)
BULK_SRC_GET_AND_SET_ORGMOD(Serogroup,serogroup)
BULK_SRC_GET_AND_SET_ORGMOD(Serotype,serotype)
BULK_SRC_GET_AND_SET_ORGMOD(Serovar,serovar)
BULK_SRC_GET_AND_SET_ORGMOD(Specimen_voucher,specimen_voucher)
BULK_SRC_GET_AND_SET_ORGMOD(Strain,strain)
BULK_SRC_GET_AND_SET_ORGMOD(Subgroup,subgroup)
BULK_SRC_GET_AND_SET_ORGMOD(Sub_species,sub_species)
BULK_SRC_GET_AND_SET_ORGMOD(Substrain,substrain)
BULK_SRC_GET_AND_SET_ORGMOD(Subtype,subtype)
BULK_SRC_GET_AND_SET_ORGMOD(Synonym,synonym)
BULK_SRC_GET_AND_SET_ORGMOD(Teleomorph,teleomorph)
BULK_SRC_GET_AND_SET_ORGMOD(Type,type)
BULK_SRC_GET_AND_SET_ORGMOD(Variety,variety)
BULK_SRC_GET_AND_SET_ORGMOD(OrgModNote,other)
BULK_SRC_GET_AND_SET_SUBSRC(Cell_line,cell_line)
BULK_SRC_GET_AND_SET_SUBSRC(Cell_type,cell_type)
BULK_SRC_GET_AND_SET_SUBSRC(Chromosome,chromosome)
BULK_SRC_GET_AND_SET_SUBSRC(Clone,clone)
BULK_SRC_GET_AND_SET_SUBSRC(Clone_lib,clone_lib)
BULK_SRC_GET_AND_SET_SUBSRC(Collected_by,collected_by)
BULK_SRC_GET_AND_SET_SUBSRC(Collection_date,collection_date)
BULK_SRC_GET_AND_SET_SUBSRC(Country,country)
BULK_SRC_GET_AND_SET_SUBSRC(Dev_stage,dev_stage)
BULK_SRC_GET_AND_SET_SUBSRC(Endogenous_virus_name,endogenous_virus_name)
BULK_SRC_GET_AND_SET_SUBSRC_TF(Environmental_sample,environmental_sample)
BULK_SRC_GET_AND_SET_SUBSRC(Frequency,frequency)
BULK_SRC_GET_AND_SET_SUBSRC(Genotype,genotype)
BULK_SRC_GET_AND_SET_SUBSRC_TF(Germline,germline)
BULK_SRC_GET_AND_SET_SUBSRC(Haplogroup,haplogroup)
BULK_SRC_GET_AND_SET_SUBSRC(Haplotype,haplotype)
BULK_SRC_GET_AND_SET_SUBSRC(Identified_by,identified_by)
BULK_SRC_GET_AND_SET_SUBSRC(Isolation_source,isolation_source)
BULK_SRC_GET_AND_SET_SUBSRC(Lab_host,lab_host)
BULK_SRC_GET_AND_SET_SUBSRC(Lat_Lon,lat_lon)
BULK_SRC_GET_AND_SET_SUBSRC(Linkage_group,linkage_group)
BULK_SRC_GET_AND_SET_SUBSRC(Map,map)
BULK_SRC_GET_AND_SET_SUBSRC(Mating_type,mating_type)
BULK_SRC_GET_AND_SET_SUBSRC_TF(Metagenomic,metagenomic)
BULK_SRC_GET_AND_SET_SUBSRC(Plasmid_name,plasmid_name)
BULK_SRC_GET_AND_SET_SUBSRC(Pop_variant,pop_variant)
BULK_SRC_GET_AND_SET_SUBSRC_TF(Rearranged,rearranged)
BULK_SRC_GET_AND_SET_SUBSRC(Segment,segment)
BULK_SRC_GET_AND_SET_SUBSRC(Sex,sex)
BULK_SRC_GET_AND_SET_SUBSRC(Subclone,subclone)
BULK_SRC_GET_AND_SET_SUBSRC(Tissue_lib,tissue_lib)
BULK_SRC_GET_AND_SET_SUBSRC(Tissue_type,tissue_type)
BULK_SRC_GET_AND_SET_SUBSRC_TF(Transgenic,transgenic)
BULK_SRC_GET_AND_SET_SUBSRC(SubSrcNote,other)


#define BULK_SRC_SUBSRC_FIELD_TF(Name) \
 BulkSrcSet##Name, BulkSetSimpleTextString, BulkSrcGet##Name, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy


static ValNodePtr StringsToDbxrefList (ValNodePtr strings)
{
  ValNodePtr dbxref_list = NULL, vnp;
  DbtagPtr dbtag;

  for (vnp = strings; vnp != NULL; vnp = vnp->next)
  {
    dbtag = DbtagNew();
    SetDbtagString (dbtag, vnp->data.ptrvalue, ExistingTextOption_replace_old);
    ValNodeAddPointer (&dbxref_list, vnp->choice, dbtag);
  }
  return dbxref_list;
}


static ValNodePtr DbxrefListToStrings (ValNodePtr dbxref_list)
{
  ValNodePtr strings = NULL, vnp;
  DbtagPtr dbtag;

  for (vnp = dbxref_list; vnp != NULL; vnp = vnp->next)
  {
    dbtag = (DbtagPtr) vnp->data.ptrvalue;
    ValNodeAddPointer (&strings, vnp->choice, GetDbtagString (dbtag));
  }
  return strings;
}


static ValNodePtr FreeDbxrefList (ValNodePtr list)
{
  ValNodePtr vnp;
 
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    vnp->data.ptrvalue = DbtagFree (vnp->data.ptrvalue);
  }
  list = ValNodeFree (list);
  return list;
}




static Pointer BulkSrcGetDbxref (Uint1 data_choice, Pointer data, Pointer metadata)
{
  SeqDescrPtr  sdp = (SeqDescrPtr) data;
  BioSourcePtr biop;
  ValNodePtr   dbxref_list = NULL;

  if (sdp == NULL || sdp->choice != Seq_descr_source 
      || (biop = (BioSourcePtr)sdp->data.ptrvalue) == NULL
      || biop->org == NULL) 
  {
    return NULL;
  }

  dbxref_list = DbxrefListToStrings (biop->org->db);
    
  return dbxref_list;
}


static void BulkSrcSetDbxref (Pointer target, Pointer data)
{
  ValNodePtr   dbxref_list;
  SeqDescrPtr  sdp = (SeqDescrPtr) target;
  BioSourcePtr biop;

  if (sdp == NULL || sdp->choice != Seq_descr_source || sdp->data.ptrvalue == NULL) return;

  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;
  dbxref_list = (ValNodePtr) data;
  if (biop->org == NULL && dbxref_list == NULL) {
    return;
  }
  if (biop->org == NULL) {
    biop->org = OrgRefNew();
  }

  /* erase previous dbxrefs */
  biop->org->db = FreeDbxrefList (biop->org->db);
  biop->org->db = StringsToDbxrefList (dbxref_list);
}


static Int4 BulkSrcFormatDbxref (ColPtr col, CharPtr name)
{
  if (col == NULL) return 0;

  col->pixWidth = MAX (15, StringLen (name)) * MaxCharWidth();
  col->pixInset = 0;
  col->charWidth = 0;
  col->charInset = 0;
  col->font = NULL;
  col->just = 'l';
  col->wrap = 1;
  col->bar = 0;
  col->underline = 0;
  col->left = 0;
  return col->pixWidth;
}


static CharPtr BulkSrcDisplayDbxref (Pointer data)
{
  ValNodePtr dbxref_list, vnp;
  CharPtr    dbxref_text;
  Int4       text_len = 1;

  if (data == NULL) return NULL;
  dbxref_list = (ValNodePtr) data;
  for (vnp = dbxref_list; vnp != NULL; vnp = vnp->next) {
    text_len += StringLen (vnp->data.ptrvalue) + 1;
  }
  dbxref_text = (CharPtr) MemNew (sizeof (Char) * text_len);
  dbxref_text[0] = 0;
  for (vnp = dbxref_list; vnp != NULL; vnp = vnp->next) {
    StringCat (dbxref_text, vnp->data.ptrvalue);
    if (vnp->next != NULL) {
      StringCat (dbxref_text, ";");
    }
  }
  return dbxref_text;
}


static Pointer BulkSrcSetDbxrefString (Pointer curr_val, ApplyValuePtr avp)
{
  CharPtr        curr_val_str;
  CharPtr        new_str, name_str, cp;
  ValNodePtr     name_list;

  name_list = (ValNodePtr) curr_val;
  curr_val_str = BulkSrcDisplayDbxref (name_list);
  name_list = ValNodeFreeData (name_list);

  new_str = HandleApplyValue (curr_val_str, avp);

  name_str = new_str;
  cp = StringChr (name_str, ';');
  while (cp != NULL) {
    *cp = 0;
    ValNodeAddPointer (&name_list, 0, StringSave (name_str));
    *cp = ';';
    name_str = cp + 1;
    cp = StringChr (name_str, ';');
  }
  ValNodeAddPointer (&name_list, 0, StringSave (name_str));
  new_str = MemFree (new_str);

  return name_list;
}


static void BulkSrcFreeDbxref (Pointer data)
{
  ValNodePtr prot_name_list = (ValNodePtr) data;
 
  ValNodeFreeData (prot_name_list);
}


static Pointer BulkSrcCopyDbxref (Pointer data)
{
  ValNodePtr new_list = NULL, list, vnp;
 
  list = (ValNodePtr) data;
  
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    ValNodeAddPointer (&new_list, vnp->choice, StringSave (vnp->data.ptrvalue));
  }
  return new_list;
}  


static DialoG BulkSrcDbxrefDialog (GrouP g, CharPtr name, SeqEntryPtr sep)
{
  return CreateVisibleStringDialog (g, 3, -1, 25);
}


static BulkEdFieldData bulk_src_desc_fields[] = {
  { "tax name", BulkSrcSetTaxNameDesc, BulkSetSimpleTextString, BulkSrcGetTaxNameDesc, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "acronym", BulkSrcSetAcronym, BulkSetSimpleTextString, BulkSrcGetAcronym, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "anamorph", BulkSrcSetAnamorph, BulkSetSimpleTextString, BulkSrcGetAnamorph, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "authority", BulkSrcSetAuthority, BulkSetSimpleTextString, BulkSrcGetAuthority, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "biotype", BulkSrcSetBiotype, BulkSetSimpleTextString, BulkSrcGetBiotype, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "biovar", BulkSrcSetBiovar, BulkSetSimpleTextString, BulkSrcGetBiovar, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "bio_material", BulkSrcSetBio_material, BulkSetSimpleTextString, BulkSrcGetBio_material, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "breed", BulkSrcSetBreed, BulkSetSimpleTextString, BulkSrcGetBreed, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "cell_line", BulkSrcSetCell_line, BulkSetSimpleTextString, BulkSrcGetCell_line, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "cell_type", BulkSrcSetCell_type, BulkSetSimpleTextString, BulkSrcGetCell_type, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "chemovar", BulkSrcSetChemovar, BulkSetSimpleTextString, BulkSrcGetChemovar, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "chromosome", BulkSrcSetChromosome, BulkSetSimpleTextString, BulkSrcGetChromosome, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "clone", BulkSrcSetClone, BulkSetSimpleTextString, BulkSrcGetClone, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "clone_lib", BulkSrcSetClone_lib, BulkSetSimpleTextString, BulkSrcGetClone_lib, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "collected_by", BulkSrcSetCollected_by, BulkSetSimpleTextString, BulkSrcGetCollected_by, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "collection_date", BulkSrcSetCollection_date, BulkSetSimpleTextString, BulkSrcGetCollection_date, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "common", BulkSrcSetCommon, BulkSetSimpleTextString, BulkSrcGetCommon, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "country", BulkSrcSetCountry, BulkSetSimpleTextString, BulkSrcGetCountry, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "cultivar", BulkSrcSetCultivar, BulkSetSimpleTextString, BulkSrcGetCultivar, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "culture_collection", BulkSrcSetCulture_collection, BulkSetSimpleTextString, BulkSrcGetCulture_collection, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "dbxref", BulkSrcSetDbxref, BulkSrcSetDbxrefString, BulkSrcGetDbxref, BulkSrcDisplayDbxref, BulkSrcFreeDbxref, BulkSrcDbxrefDialog, BulkSrcFormatDbxref, NULL, NULL, BulkSrcCopyDbxref },
  { "dev_stage", BulkSrcSetDev_stage, BulkSetSimpleTextString, BulkSrcGetDev_stage, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "ecotype", BulkSrcSetEcotype, BulkSetSimpleTextString, BulkSrcGetEcotype, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "endogenous_virus_name", BulkSrcSetEndogenous_virus_name, BulkSetSimpleTextString, BulkSrcGetEndogenous_virus_name, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "environmental-sample", BULK_SRC_SUBSRC_FIELD_TF(Environmental_sample) },
  { "forma", BulkSrcSetForma, BulkSetSimpleTextString, BulkSrcGetForma, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "forma_specialis", BulkSrcSetForma_specialis, BulkSetSimpleTextString, BulkSrcGetForma_specialis, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "frequency", BulkSrcSetFrequency, BulkSetSimpleTextString, BulkSrcGetFrequency, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "genotype", BulkSrcSetGenotype, BulkSetSimpleTextString, BulkSrcGetGenotype, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "germline", BULK_SRC_SUBSRC_FIELD_TF(Germline) },
  { "group", BulkSrcSetGroup, BulkSetSimpleTextString, BulkSrcGetGroup, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "haplogroup", BulkSrcSetHaplogroup, BulkSetSimpleTextString, BulkSrcGetHaplogroup, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "haplotype", BulkSrcSetHaplotype, BulkSetSimpleTextString, BulkSrcGetHaplotype, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "host", BulkSrcSetHost, BulkSetSimpleTextString, BulkSrcGetHost, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "identified_by", BulkSrcSetIdentified_by, BulkSetSimpleTextString, BulkSrcGetIdentified_by, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "isolate", BulkSrcSetIsolate, BulkSetSimpleTextString, BulkSrcGetIsolate, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "isolation_source", BulkSrcSetIsolation_source, BulkSetSimpleTextString, BulkSrcGetIsolation_source, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "lab_host", BulkSrcSetLab_host, BulkSetSimpleTextString, BulkSrcGetLab_host, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "lat_Lon", BulkSrcSetLat_Lon, BulkSetSimpleTextString, BulkSrcGetLat_Lon, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "linkage_group", BulkSrcSetLinkage_group, BulkSetSimpleTextString, BulkSrcGetLinkage_group, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "map", BulkSrcSetMap, BulkSetSimpleTextString, BulkSrcGetMap, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "mating_type", BulkSrcSetMating_type, BulkSetSimpleTextString, BulkSrcGetMating_type, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "metagenome_source", BulkSrcSetMetagenome_source, BulkSetSimpleTextString, BulkSrcGetMetagenome_source, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "metagenomic", BULK_SRC_SUBSRC_FIELD_TF(Metagenomic) },
  { "note (orgmod)", BulkSrcSetOrgModNote, BulkSetSimpleTextString, BulkSrcGetOrgModNote, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "note (subsrc)", BulkSrcSetSubSrcNote, BulkSetSimpleTextString, BulkSrcGetSubSrcNote, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "pathovar", BulkSrcSetPathovar, BulkSetSimpleTextString, BulkSrcGetPathovar, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "plasmid_name", BulkSrcSetPlasmid_name, BulkSetSimpleTextString, BulkSrcGetPlasmid_name, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "pop_variant", BulkSrcSetPop_variant, BulkSetSimpleTextString, BulkSrcGetPop_variant, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "rearranged", BULK_SRC_SUBSRC_FIELD_TF(Rearranged) } ,
  { "segment", BulkSrcSetSegment, BulkSetSimpleTextString, BulkSrcGetSegment, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "serogroup", BulkSrcSetSerogroup, BulkSetSimpleTextString, BulkSrcGetSerogroup, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "serotype", BulkSrcSetSerotype, BulkSetSimpleTextString, BulkSrcGetSerotype, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "serovar", BulkSrcSetSerovar, BulkSetSimpleTextString, BulkSrcGetSerovar, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "sex", BulkSrcSetSex, BulkSetSimpleTextString, BulkSrcGetSex, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "specimen_voucher", BulkSrcSetSpecimen_voucher, BulkSetSimpleTextString, BulkSrcGetSpecimen_voucher, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "strain", BulkSrcSetStrain, BulkSetSimpleTextString, BulkSrcGetStrain, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "subclone", BulkSrcSetSubclone, BulkSetSimpleTextString, BulkSrcGetSubclone, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "subgroup", BulkSrcSetSubgroup, BulkSetSimpleTextString, BulkSrcGetSubgroup, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "substrain", BulkSrcSetSubstrain, BulkSetSimpleTextString, BulkSrcGetSubstrain, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "subtype", BulkSrcSetSubtype, BulkSetSimpleTextString, BulkSrcGetSubtype, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "sub_species", BulkSrcSetSub_species, BulkSetSimpleTextString, BulkSrcGetSub_species, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "synonym", BulkSrcSetSynonym, BulkSetSimpleTextString, BulkSrcGetSynonym, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "teleomorph", BulkSrcSetTeleomorph, BulkSetSimpleTextString, BulkSrcGetTeleomorph, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "tissue_lib", BulkSrcSetTissue_lib, BulkSetSimpleTextString, BulkSrcGetTissue_lib, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "tissue_type", BulkSrcSetTissue_type, BulkSetSimpleTextString, BulkSrcGetTissue_type, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "transgenic", BULK_SRC_SUBSRC_FIELD_TF(Transgenic) },
  { "type", BulkSrcSetType, BulkSetSimpleTextString, BulkSrcGetType, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "variety", BulkSrcSetVariety, BulkSetSimpleTextString, BulkSrcGetVariety, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};


static BulkEdFieldData bulk_src_feat_fields[] = {
  { "tax name", BulkSrcSetTaxNameFeat, BulkSetSimpleTextString, BulkSrcGetTaxNameFeat, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

typedef struct bulkeditor {
  FORM_MESSAGE_BLOCK
  DialoG         bulk_ed;
  ButtoN         do_log;
  BulkEdFieldPtr field_list;
} BulkEditorData, PNTR BulkEditorPtr;

typedef struct collectfeat {
  Uint1 featdef;
  Uint1 featchoice;
  ValNodePtr feat_list;
} CollectFeatData, PNTR CollectFeatPtr;

static void CollectFeaturesForBulkEditor (BioseqPtr bsp, Pointer data)
{
  CollectFeatPtr cfp;
  SeqFeatPtr     sfp;
  SeqMgrFeatContext fcontext;

  if (bsp == NULL || data == NULL) return;
  cfp = (CollectFeatPtr) data;

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, cfp->featchoice, cfp->featdef, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, cfp->featchoice, cfp->featdef, &fcontext)) {
    ValNodeAddPointer (&(cfp->feat_list), OBJ_SEQFEAT, sfp);
  }
}

static void AcceptBulkEditor (ButtoN b)
{
  BulkEditorPtr      bep;

  bep = (BulkEditorPtr) GetObjectExtra (b);
  if (bep == NULL) return;

  WatchCursor();
  Update();

  ApplyBulkEditorToObjectList (bep->bulk_ed, bep->do_log == NULL ? FALSE : GetStatus (bep->do_log));

  ObjMgrSetDirtyFlag (bep->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bep->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
  Remove (bep->form);
}

extern Uint1 GetSubtypeForBulkEdit (ValNodePtr feat_list)
{
  ValNodePtr vnp_check;
  SeqFeatPtr sfp;
  Uint1      subtype = FEATDEF_BAD;

  for (vnp_check = feat_list; vnp_check != NULL; vnp_check = vnp_check->next) {
    if (vnp_check->choice != OBJ_SEQFEAT) continue;
    sfp = vnp_check->data.ptrvalue;
    if (sfp->data.choice == SEQFEAT_RNA) {
      /* note - using FEATDEF_rRNA to represent all editable RNA features */
      if (!IsBulkEditableRNA(sfp)) {
        return FEATDEF_BAD;
      } else if (subtype == FEATDEF_BAD) {
        subtype = FEATDEF_rRNA;
      } else if (subtype != FEATDEF_rRNA) {
        return FEATDEF_BAD;
      }
    } else if (subtype == FEATDEF_BAD) {
      subtype = sfp->idx.subtype;
    } else if (subtype != sfp->idx.subtype) {
      return FEATDEF_BAD;
    }
  }
  return subtype;
}

static void ExportBulkEditorTableBtn (ButtoN b)
{
  BulkEditorPtr bep;

  bep = (BulkEditorPtr) GetObjectExtra (b);
  if (bep == NULL) {
    return;
  }

  ExportBulkEditorTable (bep->bulk_ed);

}


extern void BulkEditorFeatList (Uint2 entityID, ValNodePtr feat_list)
{
  GrouP              c, g;
  ButtoN             b;
  BulkEditorPtr      bep;
  GrouP              h;
  WindoW             w;
  SeqEntryPtr        sep;
  BulkEdFieldPtr     field_list = NULL;
  Uint1              subtype = 0;

  if (feat_list == NULL) {
    Message (MSG_ERROR, "No features found!");
    return;
  }
  /* Create a new window, and a struct */
  /* to pass around the data in.       */

  sep = GetTopSeqEntryForEntityID (entityID);

  subtype = GetSubtypeForBulkEdit (feat_list);
  if (subtype == FEATDEF_BAD) {
    Message (MSG_ERROR, "Can't bulk edit a list of mixed feature types!");
    return;
  }
  if (subtype == FEATDEF_CDS) {
    field_list = bulk_cds_fields;
  } else if (subtype == FEATDEF_GENE) {
    field_list = bulk_gene_fields;
  } else if (subtype == FEATDEF_rRNA) {
    /* Note - using FEATDEF_rRNA to represent all editable rRNA features */
    field_list = bulk_rna_fields;
  } else if (subtype == FEATDEF_BIOSRC) {
    field_list = bulk_src_feat_fields;
  } else {
    Message (MSG_ERROR, "No bulk editor for this feature type!");
    return;
  }

  bep = (BulkEditorPtr) MemNew (sizeof (BulkEditorData));
  if (bep == NULL)
    return;

  bep->field_list = field_list;

  w = FixedWindow (-50, -33, -10, -10, "Bulk Editor",
           StdCloseWindowProc);
  SetObjectExtra (w, bep, StdCleanupFormProc); 
  bep->form = (ForM) w;
  bep->input_entityID = entityID;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  /* want: Document with columns for pseudo (when appropriate), location, fields */
  /* when user clicks on field, should be able to edit contents */
  
  bep->bulk_ed = CreateBulkEditorDialog (h, bep->field_list, feat_list, sep, TRUE, NULL, NULL);

  if (subtype == FEATDEF_CDS) {
    bep->do_log = CheckBox (h, "Log Changes", NULL);
  } else {
    bep->do_log = NULL;
  }

  g = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (g, "Export Table", ExportBulkEditorTableBtn);
  SetObjectExtra (b, bep, NULL);

  /* Add Accept and Cancel buttons */

  c = HiddenGroup (h, 3, 0, NULL);
  b = PushButton (c, "Accept", AcceptBulkEditor);
  SetObjectExtra (b, bep, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  /* Line things up nicely */

  AlignObjects (ALIGN_CENTER, (HANDLE) bep->bulk_ed,
                              (HANDLE) c, (HANDLE) g, 
                              (HANDLE) bep->do_log, NULL);


  /* Display the window now */

  RealizeWindow (w);
  Show (w);
  Select (w);
  Update ();
}


extern Uint1 GetSubtypeForBulkDescrEdit (ValNodePtr descr_list)
{
  ValNodePtr vnp_check;
  SeqDescrPtr sdp;
  Uint1      subtype = 0;

  for (vnp_check = descr_list; vnp_check != NULL; vnp_check = vnp_check->next) {
    if (vnp_check->choice != OBJ_SEQDESC) continue;
    sdp = vnp_check->data.ptrvalue;
    if (subtype == 0) {
      subtype = sdp->choice;
    } else if (sdp->choice != subtype) {
      return 0;
    }
  }
  return subtype;
}


extern void BulkEditorDescrList (Uint2 entityID, ValNodePtr descr_list)
{
  GrouP              c;
  ButtoN             b;
  BulkEditorPtr      bep;
  GrouP              h;
  WindoW             w;
  SeqEntryPtr        sep;
  BulkEdFieldPtr     field_list = NULL;
  Uint1              subtype = 0;

  if (descr_list == NULL) {
    Message (MSG_ERROR, "No descriptors found!");
    return;
  }
  /* Create a new window, and a struct */
  /* to pass around the data in.       */

  sep = GetTopSeqEntryForEntityID (entityID);

  subtype = GetSubtypeForBulkDescrEdit (descr_list);
  if (subtype == 0) {
    Message (MSG_ERROR, "Can't bulk edit a list of mixed descriptor types!");
    return;
  }
  if (subtype == Seq_descr_source) {
    field_list = bulk_src_desc_fields;
  } else {
    Message (MSG_ERROR, "No bulk editor for this descriptor type!");
    return;
  }

  bep = (BulkEditorPtr) MemNew (sizeof (BulkEditorData));
  if (bep == NULL)
    return;

  bep->field_list = field_list;

  w = FixedWindow (-50, -33, -10, -10, "Bulk Editor",
           StdCloseWindowProc);
  SetObjectExtra (w, bep, StdCleanupFormProc); 
  bep->form = (ForM) w;
  bep->input_entityID = entityID;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  /* want: Document with columns for pseudo (when appropriate), location, fields */
  /* when user clicks on field, should be able to edit contents */
  
  bep->bulk_ed = CreateBulkEditorDialogEx (h, bep->field_list, descr_list, sep, TRUE, NULL, NULL, NULL, subtype == Seq_descr_source);

  /* Add Accept and Cancel buttons */

  c = HiddenGroup (h, 3, 0, NULL);
  b = PushButton (c, "Accept", AcceptBulkEditor);
  SetObjectExtra (b, bep, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  /* Line things up nicely */

  AlignObjects (ALIGN_CENTER, (HANDLE) bep->bulk_ed,
                              (HANDLE) c, NULL);


  /* Display the window now */

  RealizeWindow (w);
  Show (w);
  Select (w);
  Update ();
}


NLM_EXTERN void BulkEditorObjectList (Uint2 entityID, CharPtr title, ValNodePtr feat_list, BulkEdFieldPtr field_list)
{
  GrouP              c;
  ButtoN             b;
  BulkEditorPtr      bep;
  GrouP              h;
  WindoW             w;
  SeqEntryPtr        sep;
  Uint1              subtype = 0;

  if (feat_list == NULL) {
    Message (MSG_ERROR, "No features found!");
    return;
  }
  /* Create a new window, and a struct */
  /* to pass around the data in.       */

  sep = GetTopSeqEntryForEntityID (entityID);

  bep = (BulkEditorPtr) MemNew (sizeof (BulkEditorData));
  if (bep == NULL)
    return;

  bep->field_list = field_list;

  w = FixedWindow (-50, -33, -10, -10, title,
           StdCloseWindowProc);
  SetObjectExtra (w, bep, StdCleanupFormProc); 
  bep->form = (ForM) w;
  bep->input_entityID = entityID;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  /* want: Document with columns for pseudo (when appropriate), location, fields */
  /* when user clicks on field, should be able to edit contents */
  
  bep->bulk_ed = CreateBulkEditorDialog (h, bep->field_list, feat_list, sep, TRUE, NULL, NULL);

  /* Add Accept and Cancel buttons */

  c = HiddenGroup (h, 3, 0, NULL);
  b = PushButton (c, "Accept", AcceptBulkEditor);
  SetObjectExtra (b, bep, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  /* Line things up nicely */

  AlignObjects (ALIGN_CENTER, (HANDLE) bep->bulk_ed,
                              (HANDLE) c, NULL);


  /* Display the window now */

  RealizeWindow (w);
  Show (w);
  Select (w);
  Update ();
}


static void BulkEditor (IteM i, Uint1 featchoice, Uint1 featdef)
{
  BaseFormPtr        bfp;
  SeqEntryPtr        sep;
  CollectFeatData    cfd;
  BulkEdFieldPtr     field_list = NULL;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  /* Create a new window, and a struct */
  /* to pass around the data in.       */

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  if (featchoice == SEQFEAT_CDREGION || featdef == FEATDEF_CDS) {
    field_list = bulk_cds_fields;
  } else if (featchoice == SEQFEAT_GENE || featdef == FEATDEF_GENE) {
    field_list = bulk_gene_fields;
  } else if (featchoice == SEQFEAT_RNA) {
    field_list = bulk_rna_fields;
  } else {
    Message (MSG_ERROR, "No bulk editor for this feature type!");
    return;
  }

  cfd.featchoice = featchoice;
  cfd.featdef = featdef;
  cfd.feat_list = NULL;
  VisitBioseqsInSep (sep, &cfd, CollectFeaturesForBulkEditor);

  if (cfd.feat_list == NULL) {
    Message (MSG_ERROR, "No features found!");
    return;
  }

  BulkEditorFeatList (bfp->input_entityID, cfd.feat_list);
}

extern void BulkEditCDS (IteM i)
{
  BulkEditor (i, SEQFEAT_CDREGION, FEATDEF_CDS);
}

extern void BulkEditGene (IteM i)
{
  BulkEditor (i, SEQFEAT_GENE, FEATDEF_GENE);
}

extern void BulkEditRNA (IteM i)
{
  BulkEditor (i, SEQFEAT_RNA, 0);
}


static void GetBulkEditSourceCallback (SeqDescrPtr sdp, Pointer data)
{
  if (sdp != NULL && sdp->choice == Seq_descr_source && data != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQDESC, sdp);
  }
}


NLM_EXTERN void BulkEditSourceBaseForm (BaseFormPtr bfp)
{
  SeqEntryPtr  sep;
  ValNodePtr   list = NULL;

  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitDescriptorsInSep (sep, &list, GetBulkEditSourceCallback);

  BulkEditorDescrList (bfp->input_entityID, list);
}


extern void BulkEditSource (IteM i)
{
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  BulkEditSourceBaseForm (bfp);
}


static Pointer GetBarcodeTestBarcodeID (Uint1 data_choice, Pointer data, Pointer metadata)
{
  BarcodeTestResultsPtr res = (BarcodeTestResultsPtr) data;
  if (res == NULL || res->bsp == NULL) {
    return NULL;
  } else {
    return BarcodeTestBarcodeIdString (res->bsp);
  }
}

  
static Pointer GetBarcodeTestGenbankID (Uint1 data_choice, Pointer data, Pointer metadata)
{
  BarcodeTestResultsPtr res = (BarcodeTestResultsPtr) data;
  if (res == NULL || res->bsp == NULL) {
    return NULL;
  } else {
    return BarcodeTestGenbankIdString (res->bsp);
  }
}


static Pointer GetBarcodeTestLengthResult (Uint1 data_choice, Pointer data, Pointer metadata)
{
  BarcodeTestResultsPtr res = (BarcodeTestResultsPtr) data;
  if (res == NULL || !res->failed_tests[eBarcodeTest_Length]) {
    return NULL;
  } else {
    return StringSave ("TRUE");
  }
}


static Pointer GetBarcodeTestPrimersResult (Uint1 data_choice, Pointer data, Pointer metadata)
{
  BarcodeTestResultsPtr res = (BarcodeTestResultsPtr) data;
  if (res == NULL || !res->failed_tests[eBarcodeTest_Primers]) {
    return NULL;
  } else {
    return StringSave ("TRUE");
  }
}


static Pointer GetBarcodeTestCountryResult (Uint1 data_choice, Pointer data, Pointer metadata)
{
  BarcodeTestResultsPtr res = (BarcodeTestResultsPtr) data;
  if (res == NULL || !res->failed_tests[eBarcodeTest_Country]) {
    return NULL;
  } else {
    return StringSave ("TRUE");
  }
}


static Pointer GetBarcodeTestSpecimenVoucherResult (Uint1 data_choice, Pointer data, Pointer metadata)
{
  BarcodeTestResultsPtr res = (BarcodeTestResultsPtr) data;
  if (res == NULL || !res->failed_tests[eBarcodeTest_SpecimenVoucher]) {
    return NULL;
  } else {
    return StringSave ("TRUE");
  }
}


static Pointer GetBarcodeTestCollectionDateResult (Uint1 data_choice, Pointer data, Pointer metadata)
{
  BarcodeTestResultsPtr res = (BarcodeTestResultsPtr) data;
  if (res == NULL || !res->failed_tests[eBarcodeTest_CollectionDate]) {
    return NULL;
  } else {
    return StringSave ("TRUE");
  }
}


static Pointer GetBarcodeTestOrderAssignmentResult (Uint1 data_choice, Pointer data, Pointer metadata)
{
  BarcodeTestResultsPtr res = (BarcodeTestResultsPtr) data;
  if (res == NULL || !res->failed_tests[eBarcodeTest_OrderAssignment]) {
    return NULL;
  } else {
    return StringSave ("TRUE");
  }
}


static Pointer GetBarcodeTestLowTraceResult (Uint1 data_choice, Pointer data, Pointer metadata)
{
  BarcodeTestResultsPtr res = (BarcodeTestResultsPtr) data;
  if (res == NULL || !res->failed_tests[eBarcodeTest_LowTrace]) {
    return NULL;
  } else {
    return StringSave ("TRUE");
  }
}


static Pointer GetBarcodeTestFrameShiftResult (Uint1 data_choice, Pointer data, Pointer metadata)
{
  BarcodeTestResultsPtr res = (BarcodeTestResultsPtr) data;
  if (res == NULL || !res->failed_tests[eBarcodeTest_FrameShift]) {
    return NULL;
  } else {
    return StringSave ("TRUE");
  }
}


static Pointer GetBarcodeTestPercentNsResult (Uint1 data_choice, Pointer data, Pointer metadata)
{
  BarcodeTestResultsPtr res = (BarcodeTestResultsPtr) data;
  Char txt[5];
  if (res == NULL || !res->failed_tests[eBarcodeTest_PercentN]) {
    return NULL;
  } else {
    sprintf (txt, "%.1f", res->n_percent);
    return StringSave (txt);
  }
}


static Pointer HasBarcodeKeyword (Uint1 data_choice, Pointer data, Pointer metadata)
{
  Boolean rval = FALSE;

  BarcodeTestResultsPtr res = (BarcodeTestResultsPtr) data;
  if (res != NULL && res->bsp != NULL) {
    rval = BioseqHasBarcodeKeyword(res->bsp);
  }
  if (rval) {
    return StringSave ("TRUE");
  } else {
    return NULL;
  }
}



static Int4 BarcodeFormatBarcodeID (ColPtr col, CharPtr name)
{
  if (col == NULL) return 0;

  col->pixWidth = MAX (21, StringLen (name)) * MaxCharWidth();
  col->pixInset = 0;
  col->charWidth = 0;
  col->charInset = 0;
  col->font = NULL;
  col->just = 'l';
  col->wrap = 1;
  col->bar = 0;
  col->underline = 0;
  col->left = 0;
  return col->pixWidth;
}


static Int4 BarcodeFormatGenbankID (ColPtr col, CharPtr name)
{
  if (col == NULL) return 0;

  col->pixWidth = MAX (10, StringLen (name)) * MaxCharWidth();
  col->pixInset = 0;
  col->charWidth = 0;
  col->charInset = 0;
  col->font = NULL;
  col->just = 'l';
  col->wrap = 1;
  col->bar = 0;
  col->underline = 0;
  col->left = 0;
  return col->pixWidth;
}


static Int4 BarcodeFormatPercentN (ColPtr col, CharPtr name)
{
  if (col == NULL) return 0;

  col->pixWidth = MAX (5, StringLen (name)) * MaxCharWidth();
  col->pixInset = 0;
  col->charWidth = 0;
  col->charInset = 0;
  col->font = NULL;
  col->just = 'l';
  col->wrap = 1;
  col->bar = 0;
  col->underline = 0;
  col->left = 0;
  return col->pixWidth;
}


static BulkEdFieldData barcode_test_fields[] = {
  { "Barcode ID", NULL, NULL, GetBarcodeTestBarcodeID, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BarcodeFormatBarcodeID, NULL, NULL, BulkSimpleTextCopy },
  { "Genbank Accession", NULL, NULL, GetBarcodeTestGenbankID, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BarcodeFormatGenbankID, NULL, NULL, BulkSimpleTextCopy },
  { "Length", NULL, NULL, GetBarcodeTestLengthResult, NULL, BulkFreeSimpleText, NULL, BulkFormatTrueFalse, BulkDrawTrueFalse, NULL, BulkSimpleTextCopy }, 
  { "Primers", NULL, NULL, GetBarcodeTestPrimersResult, NULL, BulkFreeSimpleText, NULL, BulkFormatTrueFalse, BulkDrawTrueFalse, NULL, BulkSimpleTextCopy }, 
  { "Country", NULL, NULL, GetBarcodeTestCountryResult, NULL, BulkFreeSimpleText, NULL, BulkFormatTrueFalse, BulkDrawTrueFalse, NULL, BulkSimpleTextCopy }, 
  { "Voucher", NULL, NULL, GetBarcodeTestSpecimenVoucherResult, NULL, BulkFreeSimpleText, NULL, BulkFormatTrueFalse, BulkDrawTrueFalse, NULL, BulkSimpleTextCopy }, 
  { "Percent Ns", NULL, NULL, GetBarcodeTestPercentNsResult, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BarcodeFormatPercentN, NULL, NULL, BulkSimpleTextCopy }, 
  { "Collection Date", NULL, NULL, GetBarcodeTestCollectionDateResult, NULL, BulkFreeSimpleText, NULL, BulkFormatTrueFalse, BulkDrawTrueFalse, NULL, BulkSimpleTextCopy }, 
  { "Order Assignment", NULL, NULL, GetBarcodeTestOrderAssignmentResult, NULL, BulkFreeSimpleText, NULL, BulkFormatTrueFalse, BulkDrawTrueFalse, NULL, BulkSimpleTextCopy }, 
  { "Low Trace", NULL, NULL, GetBarcodeTestLowTraceResult, NULL, BulkFreeSimpleText, NULL, BulkFormatTrueFalse, BulkDrawTrueFalse, NULL, BulkSimpleTextCopy }, 
  { "Frame Shift", NULL, NULL, GetBarcodeTestFrameShiftResult, NULL, BulkFreeSimpleText, NULL, BulkFormatTrueFalse, BulkDrawTrueFalse, NULL, BulkSimpleTextCopy }, 
  { "Has BARCODE Keyword", NULL, NULL, HasBarcodeKeyword, NULL, BulkFreeSimpleText, NULL, BulkFormatTrueFalse, BulkDrawTrueFalse, NULL, BulkSimpleTextCopy }, 
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

static void NavigateToBarcodeTestResultBioseq (ValNodePtr object, Pointer userdata)
{
  BaseFormPtr           bfp;
  BarcodeTestResultsPtr res;

  if (object == NULL) return;

  res = (BarcodeTestResultsPtr) object->data.ptrvalue;

  if (res == NULL || res->bsp == NULL) return;
  bfp = GetBaseFormForEntityID (res->bsp->idx.entityID);
  if (bfp != NULL) {
    Select (bfp->form);
    SetBioseqViewTargetByBioseq (bfp, res->bsp);
  }
}


extern DialoG BarcodeTestResultsDisplay (GrouP h, BarcodeTestConfigPtr cfg)
{
  return CreateBulkEditorDialog (h, barcode_test_fields, NULL, NULL, FALSE, NavigateToBarcodeTestResultBioseq, NULL);
}

static Pointer GetCurrentLatLon (Uint1 data_choice, Pointer data, Pointer metadata)
{
  SubSourcePtr bad_ssp;
  BioSourcePtr biop;

  biop = GetBioSourceFromObject (data_choice, data);
  if (biop == NULL) return NULL;

  bad_ssp = FindBadLatLon (biop);
  if (bad_ssp == NULL) 
  {
    return NULL;
  }
  else
  {
    return StringSave (bad_ssp->name);
  }
}


static Pointer GetCorrectedLatLon (Uint1 data_choice, Pointer data, Pointer metadata)
{
  SubSourcePtr bad_ssp;
  CharPtr      fix;
  BioSourcePtr biop;

  biop = GetBioSourceFromObject (data_choice, data);
  if (biop == NULL) return NULL;

  bad_ssp = FindBadLatLon (biop);
  if (bad_ssp == NULL) 
  {
    return NULL;
  }
  else
  {
    fix = FixLatLonFormat (bad_ssp->name);
    if (fix == NULL)
    {
      return StringSave ("Unable to autocorrect");
    }
    else
    {
      return fix;
    }
  }
}


static BulkEdFieldData latlon_fields[] = {
  { "Current Lat-lon", NULL, NULL, GetCurrentLatLon, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "Suggested Correction", NULL, NULL, GetCorrectedLatLon, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

extern DialoG LatLonTestResultsDisplay (GrouP h)
{
  return CreateBulkEditorDialog (h, latlon_fields, NULL, NULL, FALSE, ScrollToDiscrepancyItem, EditDiscrepancyItem);
}


static BulkEdFieldData trim_sequence_fields[] = {
  { "Current PercentN", NULL, NULL, GetCurrentPercentN, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "Bases to Remove", NULL, NULL, GetTrimRemove, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "Remaining Bases", NULL, NULL, GetTrimRemaining, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

extern DialoG TrimSequenceResultsDisplay (GrouP h)
{
  return CreateBulkEditorDialog (h, trim_sequence_fields, NULL, NULL, FALSE, ScrollToTrimSequenceItem, NULL);
}


static void DefaultCheckTrimSequenceResults (BulkEditorRowPtr berp)
{
  ValNodePtr vnp;

  if (berp == NULL) {
    return;
  }
  while (berp->object_list != NULL) {
    if (berp->subrows == NULL) {
      vnp = berp->object_list;
      if (vnp != NULL) {
        if (TrimShouldDelete (vnp->choice, vnp->data.ptrvalue, NULL)) {
          berp->selected = FALSE;
        } else {
          berp->selected = TRUE;
        }
      }
    } else {
      DefaultCheckTrimSequenceResults (berp->subrows);
      berp->selected = AllBulkEditorRowsSelected (berp->subrows);
    }

    berp++;
  }
}


extern void DefaultCheckTrimSequenceResultsDisplay (DialoG d)
{
  BulkEditorDlgPtr    dlg;
  RecT                r;

  dlg = (BulkEditorDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  DefaultCheckTrimSequenceResults (dlg->row_list);
    
  UpdateCheckStatus (dlg);

  ObjectRect (dlg->doc, &r);
  InvalRect (&r);  
  Update ();
}


static BulkEdFieldData trim_sequence_end_fields[] = {
  { "Accession", NULL, NULL, GetTrimAccession, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "Remove from 5'", NULL, NULL, GetTrim5, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "Remove from 3'", NULL, NULL, GetTrim3, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

extern DialoG TrimSequenceEndDisplay (GrouP h)
{
  return CreateBulkEditorDialog (h, trim_sequence_end_fields, NULL, NULL, FALSE, ScrollToTrimSequenceItem, NULL);
}



static Pointer GetTaxName (Uint1 data_choice, Pointer data, Pointer metadata)
{
  TaxFixItemPtr t;

  if ((t = (TaxFixItemPtr) data) == NULL || t->taxname == NULL) {
    return NULL;
  } else {
    return StringSave (t->taxname);
  }
}


extern Boolean IsTaxNameBad (OrgRefPtr org)
{
  Boolean rval = TRUE;
  T3DataPtr        tdp;
  Taxon3RequestPtr t3rq;
  Taxon3ReplyPtr   t3ry;
  T3ReplyPtr       trp;
  ValNodePtr       val;
  T3StatusFlagsPtr tfp;

  if (org == NULL) {
    return TRUE;
  }
  t3rq = CreateTaxon3Request (0, NULL, org);
  if (t3rq != NULL) {
    t3ry = Tax3SynchronousQuery (t3rq);
    Taxon3RequestFree (t3rq);
    if (t3ry != NULL) {
      for (trp = t3ry->reply; trp != NULL; trp = trp->next) {
        if (trp->choice == T3Reply_data) {
          tdp = (T3DataPtr) trp->data.ptrvalue;
          if (tdp != NULL) {
            for (tfp = tdp->status; tfp != NULL; tfp = tfp->next) {
              if (StringICmp (tfp->property, "is_species_level") == 0) {
                val = tfp->Value_value;
                if (val != NULL && val->choice == Value_value_bool) {
                  if (val->data.intvalue != 0) {
                    rval = FALSE;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return rval;
}


static Pointer GetTaxNameCorrection (Uint1 data_choice, Pointer data, Pointer metadata)
{
  TaxFixItemPtr t;

  if ((t = (TaxFixItemPtr) data) == NULL || t->suggested_fix == NULL) {
    return NULL;
  } else {
    return StringSave (t->suggested_fix);
  }
}


static void TaxFixSetSuggestedFix (Pointer target, Pointer data)
{
  TaxFixItemPtr  t;

  t = (TaxFixItemPtr) target;
  t->suggested_fix = MemFree (t->suggested_fix);
  t->suggested_fix = StringSave((CharPtr) data);
}


static void ScrollToTaxFixItem (ValNodePtr vnp, Pointer userdata)
{
  TaxFixItemPtr t;
  ValNode vn;

  if (vnp != NULL && (t = (TaxFixItemPtr) vnp->data.ptrvalue) != NULL) {
    vn.choice = t->data_choice;
    vn.data.ptrvalue = t->data;
    vn.next = NULL;
    ScrollToDiscrepancyItem (&vn, NULL);
  }
}


static void EditTaxFixItem(ValNodePtr vnp, Pointer userdata)
{
  TaxFixItemPtr t;
  ValNode vn;

  if (vnp != NULL && (t = (TaxFixItemPtr) vnp->data.ptrvalue) != NULL) {
    vn.choice = t->data_choice;
    vn.data.ptrvalue = t->data;
    vn.next = NULL;
    EditDiscrepancyItem (&vn, NULL);
  }
}


static BulkEdFieldData taxfix_fields[] = {
  { "Tax Name", NULL, NULL, GetTaxName, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "Suggested Correction", TaxFixSetSuggestedFix, BulkSetSimpleTextString, GetTaxNameCorrection, BulkDisplaySimpleText, BulkFreeSimpleText, BulkSimpleTextDialog, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

extern DialoG TaxFixDisplay (GrouP h)
{
  return CreateBulkEditorDialog (h, taxfix_fields, NULL, NULL, TRUE, ScrollToTaxFixItem, EditTaxFixItem);
}


static CharPtr GetSubSource (Uint1 data_choice, Pointer data, Uint1 subtype)
{
  SubSourcePtr ssp;
  BioSourcePtr biop;

  biop = GetBioSourceFromObject (data_choice, data);
  if (biop == NULL) return NULL;

  for (ssp = biop->subtype; ssp != NULL && ssp->subtype != subtype; ssp = ssp->next)
  {}

  if (ssp == NULL) 
  {
    return NULL;
  }
  else
  {
    return StringSave (ssp->name);
  }
}


static Pointer GetCountry (Uint1 data_choice, Pointer data, Pointer metadata)
{
  return GetSubSource (data_choice, data, SUBSRC_country);
}

static Pointer GetLatLon (Uint1 data_choice, Pointer data, Pointer metadata)
{
  return GetSubSource (data_choice, data, SUBSRC_lat_lon);
}


extern Pointer GetLatLonCountryCorrection (Uint1 data_choice, Pointer data, Pointer metadata)
{
  SubSourcePtr ssp;
  BioSourcePtr biop;
  CharPtr      country = NULL, cp;
  FloatHi      lat, lon;
  Boolean      found_lat_lon = FALSE;
  CharPtr      msg = NULL;
  CharPtr      guess_fmt = "Lat_lon does not map to '%s', but may be in '%s'";
  CharPtr      guess;

  biop = GetBioSourceFromObject (data_choice, data);
  if (biop == NULL) return NULL;

  for (ssp = biop->subtype; ssp != NULL && (country == NULL || !found_lat_lon); ssp = ssp->next)
  {
    if (ssp->subtype == SUBSRC_country && !StringHasNoText (ssp->name))
    {
      country = StringSave (ssp->name);
    }
    else if (ssp->subtype == SUBSRC_lat_lon)
    {
      if (ParseLatLon (ssp->name, &lat, &lon))
      {
        found_lat_lon = TRUE;
      }
    }
  }

  cp = StringChr (country, ':');
  if (cp != NULL) 
  {
    *cp = 0;
  }

  if (country == NULL && !found_lat_lon)
  {
    msg = StringSave ("Country and lat-lon not specified");
  }
  else if (country == NULL) 
  {
    msg = StringSave ("Country not specified");
  }
  else if (!found_lat_lon)
  {
    msg = StringSave ("Lat-lon not specified.");
  }
  else if (!IsCountryInLatLonList (country)) 
  {
    msg = StringSave ("Country not in lat-lon list");
  }
  else if (TestLatLonForCountry (country, lat, lon))
  {
    msg = StringSave ("No conflict");
  }
  else if (TestLatLonForCountry (country, -lat, lon)) 
  {
    if (lat < 0.0) {
      msg = StringSave ("Latitude should be set to N (northern hemisphere)");
    } else {
      msg = StringSave ("Latitude should be set to S (southern hemisphere)");
    }
  } else if (TestLatLonForCountry (country, lat, -lon)) {
    if (lon < 0.0) {
      msg = StringSave ("Longitude should be set to E (eastern hemisphere)");
    } else {
      msg = StringSave ("Longitude should be set to W (western hemisphere)");
    }
  } else if (TestLatLonForCountry (country, lon, lat)) {
    msg = StringSave ("Latitude and longitude values appear to be exchanged");
  } else {
    guess = GuessCountryForLatLon (lat, lon);
    if (guess == NULL)
    {
      msg = StringSave ("Lat-lon does not map to country");
    }
    else
    {
      msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (guess_fmt) + StringLen (country) + StringLen (guess)));
      sprintf (msg, guess_fmt, country, guess);
    }
  }
  return msg;
}

static BulkEdFieldData latloncountry_fields[] = {
  { "Lat-lon", NULL, NULL, GetLatLon, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "Country", NULL, NULL, GetCountry, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "Suggested Correction", NULL, NULL, GetLatLonCountryCorrection, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

extern DialoG LatLonCountryResultsDisplay (GrouP h)
{
  return CreateBulkEditorDialog (h, latloncountry_fields, NULL, NULL, FALSE, ScrollToDiscrepancyItem, EditDiscrepancyItem);
}


static Pointer GetCurrentSpecificHost (Uint1 data_choice, Pointer data, Pointer metadata)
{
  SpecificHostFixPtr s = (SpecificHostFixPtr) data;

  if (s == NULL || StringHasNoText (s->bad_specific_host))
  {
    return NULL;
  } 
  else
  {
    return StringSave (s->bad_specific_host);
  }
}


static Pointer GetCorrectedSpecificHost (Uint1 data_choice, Pointer data, Pointer metadata)
{
  SpecificHostFixPtr s = (SpecificHostFixPtr) data;
  CharPtr            fmt = "Unable to suggest correction for '%s'";
  CharPtr            rval = NULL;

  if (s == NULL)
  {
    return NULL;
  } 
  else if (StringHasNoText (s->new_taxname))
  {
    if (StringHasNoText (s->old_taxname)) {
      rval = StringSave ("Unable to suggest correction");
    } else {
      rval = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (s->old_taxname)));
      sprintf (rval, fmt, s->old_taxname);
    }
  }
  else
  {
    rval = StringSave (s->bad_specific_host);
    FindReplaceString (&rval, s->old_taxname, s->new_taxname, TRUE, TRUE);
  }
  return rval;
}


static Pointer GetSpecificHostCorrectionCategory  (Uint1 data_choice, Pointer data, Pointer metadata)
{
  SpecificHostFixPtr s = (SpecificHostFixPtr) data;
  CharPtr            fmt = "Unable to suggest correction for '%s'";
  CharPtr            rval = NULL;

  if (s == NULL)
  {
    return NULL;
  } 
  else 
  {
    switch (s->fix_type) {
      case eSpecificHostFix_unrecognized:
        rval = StringSave ("Unrecognized");
        break;
      case eSpecificHostFix_spelling:
        rval = StringSave ("Spelling");
        break;
      case eSpecificHostFix_capitalization:
        rval = StringSave ("Capitalization");
        break;
      case eSpecificHostFix_truncation:
        rval = StringSave ("Truncation");
        break;
      case eSpecificHostFix_replacement:
        rval = StringSave ("Replacement");
        break;
      case eSpecificHostFix_ambiguous:
        rval = StringSave ("Ambiguous");
        break;
    }
  }
  return rval;
}


static void ScrollToSpecificHostFix (ValNodePtr vnp, Pointer userdata)
{
  SpecificHostFixPtr s;
  
  if (vnp == NULL || vnp->data.ptrvalue == NULL) return;

  s = (SpecificHostFixPtr) vnp->data.ptrvalue;
  ScrollToDiscrepancyItem (s->feat_or_desc, NULL);
}


static void EditSpecificHostFix (ValNodePtr vnp, Pointer userdata)
{
  SpecificHostFixPtr s;
  
  if (vnp == NULL || vnp->data.ptrvalue == NULL) return;

  s = (SpecificHostFixPtr) vnp->data.ptrvalue;
  EditDiscrepancyItem (s->feat_or_desc, NULL);
}


static BulkEdFieldData specifichost_fields[] = {
  { "Current Specific-host", NULL, NULL, GetCurrentSpecificHost, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "Suggested Correction", NULL, NULL, GetCorrectedSpecificHost, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "Category", NULL, NULL, GetSpecificHostCorrectionCategory, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

extern DialoG SpecificHostResultsDisplay (GrouP h)
{
  return CreateBulkEditorDialog (h, specifichost_fields, NULL, NULL, FALSE, ScrollToSpecificHostFix, EditSpecificHostFix);
}


