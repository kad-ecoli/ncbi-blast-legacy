/*   asn2all.c
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
* File Name:  asn2all.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   7/26/04
*
* $Revision: 1.91 $
*
* File Description:
*
* Modifications:
* --------------------------------------------------------------------------
* ==========================================================================
*/

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <objgbseq.h>
#include <objtseq.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <explore.h>
#include <asn2gnbi.h>
#include <tofasta.h>
#include <pmfapi.h>
#include <lsqfetch.h>

#define ASN2ALL_APP_VER "7.9"

CharPtr ASN2ALL_APPLICATION = ASN2ALL_APP_VER;

static ValNodePtr DoLockFarComponents (
  SeqEntryPtr sep,
  Boolean useThreads
)

{
  ValNodePtr  rsult;
  time_t      start_time, stop_time;

  start_time = GetSecs ();

  if (NlmThreadsAvailable () && useThreads) {
    rsult = AdvcLockFarComponents (sep, TRUE, FALSE, FALSE, NULL, TRUE);
  } else if (useThreads) {
    Message (MSG_POST, "Threads not available in this executable");
    rsult = AdvcLockFarComponents (sep, TRUE, FALSE, FALSE, NULL, FALSE);
  } else {
    rsult = AdvcLockFarComponents (sep, TRUE, FALSE, FALSE, NULL, FALSE);
  }

  stop_time = GetSecs ();

  return rsult;
}

typedef enum {
  FLATFILE_FORMAT = 1,
  FASTA_FORMAT,
  CDS_FORMAT,
  GENE_FORMAT,
  TABLE_FORMAT,
  TINY_FORMAT,
  INSDSEQ_FORMAT,
  ASN_FORMAT,
  XML_FORMAT,
  CACHE_COMPONENTS
} AppFormat;

typedef struct appflags {
  AppFormat     format;
  Boolean       automatic;
  Boolean       catenated;
  Boolean       batch;
  Boolean       binary;
  Boolean       compressed;
  Boolean       lock;
  Boolean       useThreads;
  Int2          type;
  Int2          linelen;
  Int2          nearpolicy;
  ModType       mode;
  Boolean       extended;
  Boolean       failed;
  Uint4         cdsID;
  Uint4         geneID;
  Int4          filterLen;
  ValNodePtr    filterList;
  CharPtr PNTR  filterArray;
  Boolean       go_on;
  Int4          from;
  Int4          to;
  Uint1         strand;
  FILE          *nt;
  FILE          *aa;
  AsnIoPtr      an;
  AsnIoPtr      ap;
  AsnModulePtr  amp;
  AsnTypePtr    atp_bss;
  AsnTypePtr    atp_bsss;
  AsnTypePtr    atp_se;
  AsnTypePtr    atp_bsc;
  AsnTypePtr    bssp_atp;
  AsnTypePtr    atp_inst;
  AsnTypePtr    atp_insd;
  AsnTypePtr    atp_insde;
  AsnTypePtr    atp_tss;
  AsnTypePtr    atp_tsse;
  BioseqSet     bss;
  GBSeq         gbsq;
  GBSet         gbst;
  XtraBlock     xtran;
  XtraBlock     xtrap;
  TSeqSet       tss;
} AppFlagData, PNTR AppFlagPtr;

NLM_EXTERN void AsnPrintNewLine PROTO((AsnIoPtr aip));

static void DoProtFtables (
  BioseqPtr bsp,
  Pointer userdata
)

{
  AppFlagPtr  afp;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;
  afp = (AppFlagPtr) userdata;
  BioseqToGnbk (bsp, NULL, FTABLE_FMT, afp->mode, NORMAL_STYLE, 0, 0, SHOW_PROT_FTABLE, NULL, afp->aa);
}

static void SaveTinyNucStreams (
  BioseqPtr bsp,
  Pointer userdata
)

{
  AppFlagPtr  afp;

  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;
  afp = (AppFlagPtr) userdata;

  BioseqAsnWriteAsTSeq (bsp, afp->an, afp->atp_tsse);
  /*
  AsnPrintNewLine (afp->an);
  AsnIoFlush (afp->an);
  */
}

static void SaveTinyPrtStreams (
  BioseqPtr bsp,
  Pointer userdata
)

{
  AppFlagPtr  afp;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;
  afp = (AppFlagPtr) userdata;

  BioseqAsnWriteAsTSeq (bsp, afp->ap, afp->atp_tsse);
  /*
  AsnPrintNewLine (afp->ap);
  AsnIoFlush (afp->ap);
  */
}

static Boolean A2ADeltaLitOnly (
  BioseqPtr bsp
)

{
  ValNodePtr  vnp;

  if (bsp == NULL || bsp->repr != Seq_repr_delta) return FALSE;
  for (vnp = (ValNodePtr)(bsp->seq_ext); vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 1) return FALSE;
  }
  return TRUE;
}

static Boolean A2ASegHasParts (
  BioseqPtr bsp
)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   sep;

  if (bsp == NULL || bsp->repr != Seq_repr_seg) return FALSE;
  sep = bsp->seqentry;
  if (sep == NULL) return FALSE;
  sep = sep->next;
  if (sep == NULL || (! IS_Bioseq_set (sep))) return FALSE;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp != NULL && bssp->_class == BioseqseqSet_class_parts) return TRUE;
  return FALSE;
}

static void IsItFar (
  BioseqPtr bsp,
  Pointer userdata
)

{
  BoolPtr  bp;

  if (bsp == NULL || userdata == NULL) return;
  bp = (BoolPtr) userdata;

  if (bsp->repr == Seq_repr_seg && (! A2ASegHasParts (bsp))) {
    *bp = TRUE;
  } else if (bsp->repr == Seq_repr_delta && (! A2ADeltaLitOnly (bsp))) {
    *bp = TRUE;
  }
}

static void DoCDSFasta (
  BioseqPtr bsp,
  Pointer userdata
)

{
  AppFlagPtr         afp;
  Char               buf [32];
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         sfp;

  if (bsp == NULL || ! ISA_na (bsp->mol)) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
  while (sfp != NULL) {
    afp->cdsID++;
    sprintf (buf, "_cds_%ld", (long) afp->cdsID);
    CdRegionFastaStream (sfp, afp->nt,
                         STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL,
                         afp->linelen, 0, 0, TRUE, buf);
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
  }
}

static void DoTransFasta (
  BioseqPtr bsp,
  Pointer userdata
)

{
  AppFlagPtr         afp;
  Char               buf [32];
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         sfp;

  if (bsp == NULL || ! ISA_na (bsp->mol)) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
  while (sfp != NULL) {
    afp->cdsID++;
    sprintf (buf, "_prt_%ld", (long) afp->cdsID);
    TranslationFastaStream (sfp, afp->aa,
                            STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL,
                            afp->linelen, 0, 0, TRUE, buf);
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
  }
}

static void DoGeneFasta (
  BioseqPtr bsp,
  Pointer userdata
)

{
  AppFlagPtr         afp;
  Char               buf [32];
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         sfp;

  if (bsp == NULL || ! ISA_na (bsp->mol)) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &fcontext);
  while (sfp != NULL) {
    afp->geneID++;
    sprintf (buf, "_gene_%ld", (long) afp->geneID);
    GeneFastaStream (sfp, afp->nt,
                     STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL,
                     afp->linelen, 0, 0, TRUE, buf);
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_GENE, 0, &fcontext);
  }
}

  Int4          filterLen;
  ValNodePtr    filterList;
  CharPtr PNTR  filterArray;

static Boolean IdInFilter (
  CharPtr id,
  AppFlagPtr afp
)

{
  CharPtr PNTR  array;
  Int2          L, R, mid;

  if (StringHasNoText (id) || afp == NULL) return FALSE;

  array = afp->filterArray;
  if (array == NULL) return FALSE;

  L = 0;
  R = afp->filterLen - 1;

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (array [mid], id) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (array [R], id) == 0) return TRUE;

  return FALSE;
}

static void CheckFilter (
  BioseqPtr bsp,
  Pointer userdata
)

{
  AppFlagPtr  afp;
  Char        id [64];
  CharPtr     ptr;
  SeqIdPtr    sip;

  if (bsp == NULL || userdata == NULL) return;
  afp = (AppFlagPtr) userdata;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    SeqIdWrite (sip, id, PRINTID_REPORT, sizeof (id));
    ptr = StringChr (id, '.');
    if (ptr != NULL) {
      *ptr = '\0';
    }
    if (IdInFilter (id, afp)) {
      afp->go_on = TRUE;
      return;
    }
  }
}

static void GetFirstGoodBioseq (
  BioseqPtr bsp,
  Pointer userdata
)

{
  BioseqPtr PNTR bspp;

  bspp = (BioseqPtr PNTR) userdata;
  if (*bspp != NULL) return;
  *bspp = bsp;
}

static SeqLocPtr AfpToSeqLoc (
  SeqEntryPtr sep,
  AppFlagPtr afp
)

{
  BioseqPtr   bsp = NULL;
  Int4        from;
  SeqIntPtr   sintp;
  SeqLocPtr   slp = NULL;
  Uint1       strand;
  Int4        to;

  if (sep == NULL || afp == NULL) return NULL;

  if ((afp->from < 1 || afp->to < 1) && afp->strand == Seq_strand_plus) return NULL;

  if (afp->nt != NULL && afp->aa == NULL) {
    VisitSequencesInSep (sep, (Pointer) &bsp, VISIT_NUCS, GetFirstGoodBioseq);
  } else if (afp->aa != NULL && afp->nt == NULL) {
    VisitSequencesInSep (sep, (Pointer) &bsp, VISIT_PROTS, GetFirstGoodBioseq);
  }
  if (bsp == NULL) return NULL;

  from = afp->from;
  to = afp->to;
  strand = afp->strand;

  if (strand == Seq_strand_minus && from == 0 && to == 0) {
	from = 1;
	to = bsp->length;
  }
  if (from < 0) {
	from = 1;
  } else if (from > bsp->length) {
	from = bsp->length;
  }
  if (to < 0) {
	to = 1;
  } else if (to > bsp->length) {
	to = bsp->length;
  }

  sintp = SeqIntNew ();
  if (sintp == NULL) return NULL;

  sintp->from = from - 1;
  sintp->to = to - 1;
  sintp->strand = strand;
  sintp->id = SeqIdFindBest (bsp->id, 0);

  slp = ValNodeNew (NULL);
  if (slp == NULL) return NULL;

  slp->choice = SEQLOC_INT;
  slp->data.ptrvalue = (Pointer) sintp;

  return slp;
}

static void FormatRecord (
  SeqEntryPtr sep,
  AppFlagPtr afp,
  ValNodePtr bsplist
)

{
  BioseqPtr    bsp;
  CstType      custom = 0;
  Uint2        entityID;
  FlgType      flags = 0;
  Boolean      is_far = FALSE;
  LckType      locks = 0;
  SeqLocPtr    slp = NULL;
  SeqEntryPtr  top;
  ValNodePtr   vnp;

  if (sep == NULL || afp == NULL) return;

  if (afp->filterArray != NULL) {
    afp->go_on = FALSE;
    VisitBioseqsInSep (sep, (Pointer) afp, CheckFilter);
    if (! afp->go_on) return;
  }

  VisitBioseqsInSep (sep, (Pointer) &is_far, IsItFar);

  if (afp->nearpolicy == 2 && is_far) {
    flags = SHOW_CONTIG_FEATURES | ONLY_NEAR_FEATURES;
  } else {
    flags = SHOW_CONTIG_FEATURES;
  }
  if (is_far && (! afp->lock)) {
    locks = LOOKUP_FAR_COMPONENTS;
  }
  if (afp->extended) {
    flags |= REFSEQ_CONVENTIONS | SHOW_TRANCRIPTION | SHOW_PEPTIDE;
  }

  slp = AfpToSeqLoc (sep, afp);

  switch (afp->format) {
    case FLATFILE_FORMAT :
      if (afp->nt != NULL) {
        SeqEntryToGnbk (sep, slp, GENBANK_FMT, afp->mode, NORMAL_STYLE,
                        flags, locks, custom, NULL, afp->nt);
      }
      if (afp->aa != NULL) {
        SeqEntryToGnbk (sep, slp, GENPEPT_FMT, afp->mode, NORMAL_STYLE,
                        flags, 0, custom, NULL, afp->aa);
      }
      break;
    case FASTA_FORMAT :
      if (afp->nt != NULL) {
        if (afp->nearpolicy == 1 ||
            (afp->nearpolicy == 2 && (! is_far)) ||
            (afp->nearpolicy == 3 && is_far)) {
          if (slp != NULL) {
            SeqLocFastaStream (slp, afp->nt, STREAM_EXPAND_GAPS, afp->linelen, 0, 0);
          } else {
            SeqEntryFastaStream (sep, afp->nt, STREAM_EXPAND_GAPS, afp->linelen,
                                 0, 0, TRUE, FALSE, FALSE);
          }
        }
      }
      if (afp->aa != NULL) {
        if (slp != NULL) {
          SeqLocFastaStream (slp, afp->aa, STREAM_EXPAND_GAPS, afp->linelen, 0, 0);
        } else {
          SeqEntryFastaStream (sep, afp->aa, STREAM_EXPAND_GAPS, afp->linelen,
                               0, 0, FALSE, TRUE, FALSE);
        }
      }
      break;
    case CDS_FORMAT :
      if (afp->nt != NULL) {
        entityID = ObjMgrGetEntityIDForChoice (sep);
        top = GetTopSeqEntryForEntityID (entityID);
        if (top != NULL) {
          SeqMgrIndexFeatures (0, top->data.ptrvalue);
          afp->cdsID = 0;
          VisitBioseqsInSep (top, (Pointer) afp, DoCDSFasta);
        }
      }
      if (afp->aa != NULL) {
        entityID = ObjMgrGetEntityIDForChoice (sep);
        top = GetTopSeqEntryForEntityID (entityID);
        if (top != NULL) {
          SeqMgrIndexFeatures (0, top->data.ptrvalue);
          afp->cdsID = 0;
          VisitBioseqsInSep (top, (Pointer) afp, DoTransFasta);
        }
      }
      break;
    case GENE_FORMAT :
      if (afp->nt != NULL) {
        entityID = ObjMgrGetEntityIDForChoice (sep);
        top = GetTopSeqEntryForEntityID (entityID);
        if (top != NULL) {
          SeqMgrIndexFeatures (0, top->data.ptrvalue);
          afp->geneID = 0;
          VisitBioseqsInSep (top, (Pointer) afp, DoGeneFasta);
        }
      }
      break;
    case TABLE_FORMAT :
      if (afp->nt != NULL) {
        SeqEntryToGnbk (sep, slp, FTABLE_FMT, afp->mode, NORMAL_STYLE,
                        flags, locks, 0, NULL, afp->nt);
      }
      if (afp->aa != NULL) {
        VisitBioseqsInSep (sep, (Pointer) afp, DoProtFtables);
      }
      break;
    case TINY_FORMAT :
      if (afp->an != NULL) {
        VisitBioseqsInSep (sep, (Pointer) afp, SaveTinyNucStreams);
      }
      if (afp->ap != NULL) {
        VisitBioseqsInSep (sep, (Pointer) afp, SaveTinyPrtStreams);
      }
      break;
    case INSDSEQ_FORMAT :
      if (afp->an != NULL) {
        SeqEntryToGnbk (sep, slp, GENBANK_FMT, afp->mode, NORMAL_STYLE,
                        flags, locks, custom, &(afp->xtran), NULL);
      }
      if (afp->ap != NULL) {
        SeqEntryToGnbk (sep, slp, GENPEPT_FMT, afp->mode, NORMAL_STYLE,
                        flags, 0, custom, &(afp->xtrap), NULL);
      }
      break;
    case ASN_FORMAT :
    case XML_FORMAT :
      SeqEntryAsnWrite (sep, afp->an, NULL);
      break;
    case CACHE_COMPONENTS :
      if (afp->an != NULL) {
        for (vnp = bsplist; vnp != NULL; vnp = vnp->next) {
          bsp = (BioseqPtr) vnp->data.ptrvalue;
          if (bsp == NULL) continue;
          entityID = ObjMgrGetEntityIDForPointer (bsp);
          if (entityID < 1) continue;
          top = GetTopSeqEntryForEntityID (entityID);
          if (top == NULL) continue;
          SeqEntryAsnWrite (top, afp->an, afp->atp_se);
        }
      }
      break;
    default :
      break;
  }

  SeqLocFree (slp);
}

static void ProcessSingleRecord (
  CharPtr filename,
  AppFlagPtr afp
)

{
  AsnIoPtr      aip;
  BioseqPtr     bsp;
  ValNodePtr    bsplist;
  BioseqSetPtr  bssp;
  Pointer       dataptr = NULL;
  Uint2         datatype = 0, entityID = 0;
  FILE          *fp;
  ObjMgrPtr     omp;
  SeqEntryPtr   sep;

  if (afp == NULL) return;

  if (StringHasNoText (filename)) return;

  if (afp->type == 1) {
    fp = FileOpen (filename, "r");
    if (fp == NULL) {
      Message (MSG_POSTERR, "Failed to open '%s'", filename);
      return;
    }

    dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, TRUE, FALSE);

    FileClose (fp);

    entityID = ObjMgrRegister (datatype, dataptr);

  } else if (afp->type >= 2 && afp->type <= 5) {
    aip = AsnIoOpen (filename, afp->binary? "rb" : "r");
    if (aip == NULL) {
      Message (MSG_POSTERR, "AsnIoOpen failed for input file '%s'", filename);
      return;
    }

    SeqMgrHoldIndexing (TRUE);
    switch (afp->type) {
      case 2 :
        dataptr = (Pointer) SeqEntryAsnRead (aip, NULL);
        datatype = OBJ_SEQENTRY;
        break;
      case 3 :
        dataptr = (Pointer) BioseqAsnRead (aip, NULL);
        datatype = OBJ_BIOSEQ;
        break;
      case 4 :
        dataptr = (Pointer) BioseqSetAsnRead (aip, NULL);
        datatype = OBJ_BIOSEQSET;
        break;
      case 5 :
        dataptr = (Pointer) SeqSubmitAsnRead (aip, NULL);
        datatype = OBJ_SEQSUB;
        break;
      default :
        break;
    }
    SeqMgrHoldIndexing (FALSE);

    AsnIoClose (aip);

    entityID = ObjMgrRegister (datatype, dataptr);

  } else {
    Message (MSG_POSTERR, "Input format type '%d' unrecognized", (int) afp->type);
    return;
  }

  if (entityID < 1 || dataptr == NULL) {
    Message (MSG_POSTERR, "Data read failed for input file '%s'", filename);
    return;
  }

  if (datatype == OBJ_SEQSUB || datatype == OBJ_SEQENTRY ||
        datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET) {

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
      bsplist = NULL;
      if (afp->lock) {
        bsplist = DoLockFarComponents (sep, afp->useThreads);
      }

      FormatRecord (sep, afp, bsplist);

      bsplist = UnlockFarComponents (bsplist);
    }

  } else {

    Message (MSG_POSTERR, "Datatype %d not recognized", (int) datatype);
  }

  ObjMgrFree (datatype, dataptr);

  omp = ObjMgrGet ();
  ObjMgrReapOne (omp);
  SeqMgrClearBioseqIndex ();
  ObjMgrFreeCache (0);
  FreeSeqIdGiCache ();

  SeqEntrySetScope (NULL);
}

static void ProcessMultipleRecord (
  CharPtr filename,
  AppFlagPtr afp
)

{
  AsnIoPtr       aip, aop = NULL;
  AsnTypePtr     atp;
  BioseqPtr      bsp;
  ValNodePtr     bsplist;
  DataVal        dv;
  FILE           *fp;
  Boolean        io_failure = FALSE;
  ObjMgrPtr      omp;
  SeqEntryPtr    sep;
#ifdef OS_UNIX
  Char           cmmd [256];
  CharPtr        gzcatprog;
  int            ret;
  Boolean        usedPopen = FALSE;
#endif

  if (afp == NULL) return;

  if (StringHasNoText (filename)) return;

#ifndef OS_UNIX
  if (afp->compressed) {
    Message (MSG_POSTERR, "Can only decompress on-the-fly on UNIX machines");
    return;
  }
#endif

#ifdef OS_UNIX
  if (afp->compressed) {
    gzcatprog = getenv ("NCBI_UNCOMPRESS_BINARY");
    if (gzcatprog != NULL) {
      sprintf (cmmd, "%s %s", gzcatprog, filename);
    } else {
      ret = system ("gzcat -h >/dev/null 2>&1");
      if (ret == 0) {
        sprintf (cmmd, "gzcat %s", filename);
      } else if (ret == -1) {
        Message (MSG_POSTERR, "Unable to fork or exec gzcat in ScanBioseqSetRelease");
        return;
      } else {
        ret = system ("zcat -h >/dev/null 2>&1");
        if (ret == 0) {
          sprintf (cmmd, "zcat %s", filename);
        } else if (ret == -1) {
          Message (MSG_POSTERR, "Unable to fork or exec zcat in ScanBioseqSetRelease");
          return;
        } else {
          Message (MSG_POSTERR, "Unable to find zcat or gzcat in ScanBioseqSetRelease - please edit your PATH environment variable");
          return;
        }
      }
    }
    fp = popen (cmmd, /* afp->binary? "rb" : */ "r");
    usedPopen = TRUE;
  } else {
    fp = FileOpen (filename, afp->binary? "rb" : "r");
  }
#else
  fp = FileOpen (filename, afp->binary? "rb" : "r");
#endif
  if (fp == NULL) {
    Message (MSG_POSTERR, "FileOpen failed for input file '%s'", filename);
    return;
  }

  aip = AsnIoNew (afp->binary? ASNIO_BIN_IN : ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
  if (aip == NULL) {
    Message (MSG_POSTERR, "AsnIoNew failed for input file '%s'", filename);
    return;
  }

  switch (afp->format) {
    case ASN_FORMAT :
      aop = afp->an;
      break;
    case XML_FORMAT :
      aop = afp->an;
      break;
    default :
      break;
  }

  atp = afp->atp_bss;

  if (aop != NULL) {

    if (afp->format == XML_FORMAT) {
      while ((! io_failure) && (atp = AsnReadId (aip, afp->amp, atp)) != NULL) {
        if (aip->io_failure) {
          io_failure = TRUE;
          aip->io_failure = FALSE;
        }
        if (atp == afp->atp_inst) {
          /* converts compressed sequences to iupac like asn2xml */
          bsp = BioseqNew ();
          BioseqInstAsnRead (bsp, aip, atp);
          BioseqInstAsnWrite (bsp, aop, atp);
          bsp = BioseqFree (bsp);
        } else {
          AsnReadVal (aip, atp, &dv);
          AsnWrite (aop, atp, &dv);
          AsnKillValue (atp, &dv);
        }
        if (aip->io_failure) {
          io_failure = TRUE;
          aip->io_failure = FALSE;
        }
      }
    } else {
      while ((! io_failure) && (atp = AsnReadId (aip, afp->amp, atp)) != NULL) {
        if (aip->io_failure) {
          io_failure = TRUE;
          aip->io_failure = FALSE;
        }
        AsnReadVal (aip, atp, &dv);
        AsnWrite (aop, atp, &dv);
        AsnKillValue (atp, &dv);
        if (aip->io_failure) {
          io_failure = TRUE;
          aip->io_failure = FALSE;
        }
      }
    }

  } else {

    while ((! io_failure) && (atp = AsnReadId (aip, afp->amp, atp)) != NULL) {
      if (aip->io_failure) {
        io_failure = TRUE;
        aip->io_failure = FALSE;
      }
      if (atp == afp->atp_se) {

        SeqMgrHoldIndexing (TRUE);
        sep = SeqEntryAsnRead (aip, atp);
        SeqMgrHoldIndexing (FALSE);

        if (sep != NULL) {
          bsplist = NULL;
          if (afp->lock) {
            bsplist = DoLockFarComponents (sep, afp->useThreads);
          }

          FormatRecord (sep, afp, bsplist);

          bsplist = UnlockFarComponents (bsplist);
        }

        SeqEntryFree (sep);
        omp = ObjMgrGet ();
        ObjMgrReapOne (omp);
        SeqMgrClearBioseqIndex ();
        ObjMgrFreeCache (0);
        FreeSeqIdGiCache ();

        SeqEntrySetScope (NULL);

      } else {

        AsnReadVal (aip, atp, NULL);
      }

      if (aip->io_failure) {
        io_failure = TRUE;
        aip->io_failure = FALSE;
      }
    }
  }

  if (aip->io_failure) {
    io_failure = TRUE;
  }

  if (io_failure) {
    Message (MSG_POSTERR, "Asn io_failure for input file '%s'", filename);
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
}

static void FormatWrapper (
  SeqEntryPtr sep,
  Pointer userdata
)

{
  AppFlagPtr  afp;
  ValNodePtr  bsplist;

  if (sep == NULL) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  bsplist = NULL;
  if (afp->lock) {
    bsplist = DoLockFarComponents (sep, afp->useThreads);
  }

  FormatRecord (sep, afp, bsplist);

  bsplist = UnlockFarComponents (bsplist);
}

static void ProcessOneRecord (
  CharPtr filename,
  Pointer userdata
)

{
  AppFlagPtr   afp;
  Pointer      dataptr;
  Uint2        datatype;
  Uint2        entityID;
  FILE         *fp;
  SeqEntryPtr  sep;

  if (StringHasNoText (filename)) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  if (afp->automatic) {
    ReadSequenceAsnFile (filename, afp->binary, afp->compressed, (Pointer) afp, FormatWrapper);
  } else if (afp->catenated) {
    fp = FileOpen (filename, "r");
    if (fp != NULL) {
      while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE)) != NULL) {
        sep = GetTopSeqEntryForEntityID (entityID);
        FormatWrapper (sep, afp);
      }
      FileClose (fp);
    }
  } else if (afp->batch) {
    ProcessMultipleRecord (filename, afp);
  } else {
    ProcessSingleRecord (filename, afp);
  }
}

static SeqEntryPtr SeqEntryFromAccnOrGi (
  CharPtr str,
  AppFlagPtr afp
)

{
  CharPtr      accn;
  Boolean      alldigits;
  BioseqPtr    bsp;
  Char         buf [64];
  Char         ch;
  Int4         flags = 0;
  CharPtr      ptr;
  Int2         retcode = 0;
  SeqEntryPtr  sep = NULL;
  SeqIdPtr     sip;
  CharPtr      tmp1 = NULL;
  CharPtr      tmp2 = NULL;
  Int4         uid = 0;
  long int     val;
  ValNode      vn;

  if (StringHasNoText (str)) return NULL;
  StringNCpy_0 (buf, str, sizeof (buf));
  TrimSpacesAroundString (buf);

  accn = buf;
  tmp1 = StringChr (accn, ',');
  if (tmp1 != NULL) {
    *tmp1 = '\0';
    tmp1++;
    tmp2 = StringChr (tmp1, ',');
    if (tmp2 != NULL) {
      *tmp2 = '\0';
      tmp2++;
      if (StringDoesHaveText (tmp2) && sscanf (tmp2, "%ld", &val) == 1) {
        flags = (Int4) val;
      }
    }
    if (StringDoesHaveText (tmp1) && sscanf (tmp1, "%ld", &val) == 1) {
      retcode = (Int2) val;
    }
  }

  alldigits = TRUE;
  ptr = accn;
  ch = *ptr;
  while (ch != '\0') {
    if (! IS_DIGIT (ch)) {
      alldigits = FALSE;
    }
    ptr++;
    ch = *ptr;
  }

  if (alldigits) {
    if (sscanf (accn, "%ld", &val) == 1) {
      uid = (Int4) val;
    }
  } else {
    sip = SeqIdFromAccessionDotVersion (accn);
    if (sip != NULL) {
      uid = GetGIForSeqId (sip);
      SeqIdFree (sip);
    }
  }

  if (uid > 0) {
    sep = PubSeqSynchronousQuery (uid, retcode, flags);
    if (sep != NULL) {
      MemSet ((Pointer) &vn, 0, sizeof (ValNode));
      vn.choice = SEQID_GI;
      vn.data.intvalue = uid;
      bsp = BioseqFind (&vn);
      if (bsp != NULL) {
        if (afp != NULL) {
          if (afp->format == ASN_FORMAT || afp->format == XML_FORMAT) return sep;
        }
        sep = SeqMgrGetSeqEntryForData ((Pointer) bsp);
      }
    }
  }

  return sep;
}

static Boolean IsAllDigits (CharPtr str)

{
  CharPtr cp;

  if (StringHasNoText (str)) return FALSE;

  cp = str;
  while (*cp != 0 && isdigit (*cp)) {
    cp++;
  }
  if (*cp == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}

static void ReadFilterFile (
  CharPtr filterfile,
  AppFlagPtr afp
)

{
  CharPtr PNTR  array;
  FileCache     fc;
  FILE          *fp;
  ValNodePtr    head = NULL;
  Int4          i;
  ValNodePtr    last = NULL;
  Int4          len;
  Char          line [1023];
  CharPtr       ptr;
  CharPtr       str;
  Char          tmp [64];
  ValNodePtr    vnp;

  if (StringHasNoText (filterfile) || afp == NULL) return;


  fp = FileOpen (filterfile, "r");
  if (fp == NULL) return;

  if (FileCacheSetup (&fc, fp)) {
    for (str = FileCacheGetString (&fc, line, sizeof (line));
         str != NULL;
         str = FileCacheGetString (&fc, line, sizeof (line))) {
      TrimSpacesAroundString (str);
      if (StringHasNoText (str)) continue;

      if (IsAllDigits (str)) {
        sprintf (tmp, "%s", str);
        str = tmp;
      }
      ptr = StringChr (str, '.');
      if (ptr != NULL) {
        *ptr = '\0';
      }

      vnp = ValNodeCopyStr (&last, 0, str);
      if (head == NULL) {
        head = vnp;
      }
      last = vnp;
    }
  }

  FileClose (fp);

  if (head == NULL) return;

  if (! ValNodeIsSorted (head, SortVnpByString)) {
    head = ValNodeSort (head, SortVnpByString);
  }
  ValNodeUnique (&head, SortVnpByString, ValNodeFreeData);

  len = ValNodeLen (head);
  array = (CharPtr PNTR) MemNew (sizeof (CharPtr) * (len + 1));
  if (array == NULL) return;

  for (i = 0, vnp = head; i < len && vnp != NULL; i++, vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    array [i] = str;
  }
  afp->filterLen = len;
  afp->filterList = head;
  afp->filterArray = array;
}

static CharPtr helpLines [] = {
  "asn2all is primarily intended for generating reports from the binary",
  "ASN.1 Bioseq-set release files downloaded from the NCBI ftp site",
  "(ncbi-asn1 directory). It can produce GenBank and GenPept flatfiles,",
  "FASTA sequence files, INSDSet structured XML, TinySeq XML, and 5-column",
  "feature table format.",
  "",
  "The release files (which have .aso.gz suffix), should be uncompressed",
  "with gunzip, resulting in files with suffix .aso. For example,",
  "gbpri1.aso is the first file in the primate division, so the command",
  "",
  "  gunzip gbpri1.aso.gz",
  "",
  "will result in gbpri1.aso being created. The original gbpri1.aso.gz",
  "file is removed after successful decompression.",
  "",
  "In asn2all, the name of the file to be processed is specified by the -i",
  "command line argument. Use -a t to indicate that it is a release file",
  "and -b to indicate that it is binary ASN.1. A text ASN.1 file obtained",
  "from Entrez can be processed by using -a a instead of -a t -b.",
  "",
  "Nucleotide and protein records can be processed simultaneously. Use the",
  "-o argument to indicate the nucleotide output file, and the -v argument",
  "for the protein output file.",
  "",
  "The -f argument determines the format to be generated. Legal values of",
  "-f and the resulting formats are:",
  "",
  "  g  GenBank (nucleotide) or GenPept (protein)",
  "  f  FASTA",
  "  d  CDS FASTA (nucleotide) or Translated FASTA (protein)",
  "  t  5-column feature table",
  "  y  TinySet XML",
  "  s  INSDSet XML",
  "  a  ASN.1 of entire record",
  "  x  XML version of entire record",
  "",
  "The command",
  "",
  "  asn2all -i gbpri1.aso -a t -b -f g -o gbpri1.nuc -v gbpri1.prt",
  "",
  "will generate GenBank and GenPept reports from gbpri1.aso.",
  NULL
};

static void DisplayHelpText (
  void
)

{
  Int2  i;

  for (i = 0; helpLines [i] != NULL; i++) {
    printf ("%s\n", helpLines [i]);
  }
  printf ("\n");
}

/* Args structure contains command-line arguments */

typedef enum {
  p_argInputPath = 0,
  i_argInputFile,
  o_argNtOutFile,
  v_argAaOutFile,
  x_argSuffix,
  f_argFormat,
  a_argType,
  b_argBinary,
  c_argCompressed,
  r_argRemote,
  k_argLocal,
  d_argAsnIdx,
  l_argLockFar,
  T_argThreads,
  n_argNear,
  X_argExtended,
  A_argAccession,
  F_argFilterFile,
  h_argHelp,
  J_argFrom,
  K_argTo,
  M_argStrand,
} Arguments;


Args myargs [] = {
  {"Path to Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Input File Name", "stdin", NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Nucleotide Output File Name", NULL, NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Protein Output File Name", NULL, NULL, NULL,
    TRUE, 'v', ARG_FILE_OUT, 0.0, 0, NULL},
  {"File Selection Suffix", ".aso", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Format\n"
   "      g GenBank/GenPept\n"
   "      f FASTA\n"
   "      d CDS FASTA\n"
   "      e Gene FASTA\n"
   "      t Feature Table\n"
   "      y TinySet XML\n"
   "      s INSDSet XML\n"
   "      a ASN.1\n"
   "      x XML\n"
   "      c Cache Components\n", NULL, NULL, NULL,
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"ASN.1 Type\n"
   "      a Automatic\n"
   "      c Catenated\n"
   "      z Any\n"
   "      e Seq-entry\n"
   "      b Bioseq\n"
   "      s Bioseq-set\n"
   "      m Seq-submit\n"
   "      t Batch Processing\n", "a", NULL, NULL,
    TRUE, 'a', ARG_STRING, 0.0, 0, NULL},
  {"Bioseq-set is Binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Bioseq-set is Compressed", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Remote Fetching", "F", NULL, NULL,
    TRUE, 'r', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Local Fetching", "F", NULL, NULL,
    TRUE, 'k', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Path to Indexed Binary ASN.1 Data", NULL, NULL, NULL,
    TRUE, 'd', ARG_STRING, 0.0, 0, NULL},
  {"Lock Components in Advance", "F", NULL, NULL,
    TRUE, 'l', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Use Threads", "F", NULL, NULL,
    TRUE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Near Fasta Policy\n"
   "      a All\n"
   "      n Near Only\n"
   "      f Far Only\n", "n", NULL, NULL,
    TRUE, 'n', ARG_STRING, 0.0, 0, NULL},
  {"Extended Qualifier Output", "F", NULL, NULL,
    TRUE, 'X', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Accession to Fetch (or Accession,retcode,flags where flags -1 fetches external features)", NULL, NULL, NULL,
    TRUE, 'A', ARG_STRING, 0.0, 0, NULL},
  {"Accession Filter File", NULL, NULL, NULL,
    TRUE, 'F', ARG_FILE_IN, 0.0, 0, NULL},
  {"Display Help Message", "F", NULL, NULL,
    TRUE, 'h', ARG_BOOLEAN, 0.0, 0, NULL},
  {"SeqLoc From", "0", NULL, NULL,
    TRUE, 'J', ARG_INT, 0.0, 0, NULL},
  {"SeqLoc To", "0", NULL, NULL,
    TRUE, 'K', ARG_INT, 0.0, 0, NULL},
  {"SeqLoc Minus Strand", "F", NULL, NULL,
    TRUE, 'M', ARG_BOOLEAN, 0.0, 0, NULL},
};

Int2 Main (void)

{
  CharPtr      asnin, aaout, directory, suffix, ntout, accn, filterfile, asnidx, str;
  AppFlagData  afd;
  Char         app [64], format, nearpolicy, type, xmlbuf [128];
  DataVal      av;
  ValNodePtr   bsplist;
  Boolean      help, indexed, local, remote;
  SeqEntryPtr  sep;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  ErrSetLogfile ("stderr", ELOG_APPEND);
  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

  if (! AllObjLoad ()) {
    Message (MSG_FATAL, "AllObjLoad failed");
    return 1;
  }
  if (! SubmitAsnLoad ()) {
    Message (MSG_FATAL, "SubmitAsnLoad failed");
    return 1;
  }
  if (! FeatDefSetLoad ()) {
    Message (MSG_FATAL, "FeatDefSetLoad failed");
    return 1;
  }
  if (! SeqCodeSetLoad ()) {
    Message (MSG_FATAL, "SeqCodeSetLoad failed");
    return 1;
  }
  if (! GeneticCodeTableLoad ()) {
    Message (MSG_FATAL, "GeneticCodeTableLoad failed");
    return 1;
  }

  /* process command line arguments */

  sprintf (app, "asn2all %s", ASN2ALL_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  /* additional setup modifications */

  help = (Boolean) myargs [h_argHelp].intvalue;
  if (help) {
    DisplayHelpText ();
    return 0;
  }

  if (! objgbseqAsnLoad ()) {
    Message (MSG_POSTERR, "objgbseqAsnLoad failed");
    return 1;
  }
  if (! objinsdseqAsnLoad ()) {
    Message (MSG_POSTERR, "objinsdseqAsnLoad failed");
    return 1;
  }

  if (GetAppParam ("NCBI", "SETTINGS", "XMLPREFIX", NULL, xmlbuf, sizeof (xmlbuf))) {
    AsnSetXMLmodulePrefix (StringSave (xmlbuf));
  }

  MemSet ((Pointer) &afd, 0, sizeof (AppFlagData));

  remote = (Boolean ) myargs [r_argRemote].intvalue;
  local = (Boolean) myargs [k_argLocal].intvalue;
  asnidx = (CharPtr) myargs [d_argAsnIdx].strvalue;
  indexed = (Boolean) StringDoesHaveText (asnidx);
  accn = (CharPtr) myargs [A_argAccession].strvalue;
  filterfile = (CharPtr) myargs [F_argFilterFile].strvalue;

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  asnin = (CharPtr) myargs [i_argInputFile].strvalue;
  ntout = (CharPtr) myargs [o_argNtOutFile].strvalue;
  aaout = (CharPtr) myargs [v_argAaOutFile].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;

  /* default to stdout for nucleotide output if nothing specified */

  if (StringHasNoText (ntout) &&
      StringHasNoText (aaout)) {
    ntout = "stdout";
  }

  /* populate parameter structure */

  afd.automatic = FALSE;
  afd.catenated = FALSE;
  afd.batch = FALSE;
  afd.binary = (Boolean) myargs [b_argBinary].intvalue;
  afd.compressed = (Boolean) myargs [c_argCompressed].intvalue;
  afd.lock = (Boolean) myargs [l_argLockFar].intvalue;
  afd.useThreads = (Boolean) myargs [T_argThreads].intvalue;
  afd.type = 1;
  afd.linelen = 70;
  afd.nearpolicy = 1;
  afd.mode = ENTREZ_MODE;
  afd.extended = (Boolean) myargs [X_argExtended].intvalue;
  afd.failed = FALSE;

  str = myargs [f_argFormat].strvalue;
  TrimSpacesAroundString (str);
  if (StringDoesHaveText (str)) {
    format = str [0];
  } else {
    Message (MSG_POSTERR, "You must indicate a format with the -f parameter");
    return 1;
  }

  format = TO_LOWER (format);
  switch (format) {
    case 'g' :
      afd.format = FLATFILE_FORMAT;
      break;
    case 'f' :
      afd.format = FASTA_FORMAT;
      break;
    case 'd' :
      afd.format = CDS_FORMAT;
      break;
    case 'e' :
      afd.format = GENE_FORMAT;
      break;
    case 't' :
      afd.format = TABLE_FORMAT;
      break;
    case 'y' :
      afd.format = TINY_FORMAT;
      break;
    case 's' :
      afd.format = INSDSEQ_FORMAT;
      break;
    case 'a' :
      afd.format = ASN_FORMAT;
      break;
    case 'x' :
      afd.format = XML_FORMAT;
      break;
    case 'c' :
      afd.format = CACHE_COMPONENTS;
      break;
    default :
      afd.format = FLATFILE_FORMAT;
      break;
  }

  str = myargs [a_argType].strvalue;
  TrimSpacesAroundString (str);
  if (StringDoesHaveText (str)) {
    type = str [0];
  } else {
    type = 'a';
  }

  type = TO_LOWER (type);
  switch (type) {
    case 'a' :
      afd.type = 1;
      afd.automatic = TRUE;
      break;
    case 'c' :
      afd.type = 1;
      afd.catenated = TRUE;
      break;
    case 'z' :
      afd.type = 1;
      break;
    case 'e' :
      afd.type = 2;
      break;
    case 'b' :
      afd.type = 3;
      break;
    case 's' :
      afd.type = 4;
      break;
    case 'm' :
      afd.type = 5;
      break;
    case 't' :
      afd.type = 1;
      afd.batch = TRUE;
      afd.mode = RELEASE_MODE;
      break;
    default :
      afd.type = 1;
      break;
  }

  afd.from = myargs [J_argFrom].intvalue;
  afd.to = myargs [K_argTo].intvalue;
  if (myargs [M_argStrand].intvalue) {
    afd.strand = Seq_strand_minus;
  } else {
    afd.strand = Seq_strand_plus;
  }

  str = myargs [n_argNear].strvalue;
  TrimSpacesAroundString (str);
  if (StringDoesHaveText (str)) {
    nearpolicy = str [0];
  } else {
    nearpolicy = 'a';
  }

  nearpolicy = TO_LOWER (nearpolicy);
  switch (nearpolicy) {
    case 'a' :
      afd.nearpolicy = 1;
      break;
    case 'n' :
      afd.nearpolicy = 2;
      break;
    case 'f' :
      afd.nearpolicy = 3;
      break;
    default :
      afd.nearpolicy = 1;
      break;
  }

  afd.nt = NULL;
  afd.aa = NULL;
  afd.an = NULL;
  afd.ap = NULL;

  afd.amp = AsnAllModPtr ();
  afd.atp_bss = AsnFind ("Bioseq-set");
  afd.atp_bsss = AsnFind ("Bioseq-set.seq-set");
  afd.atp_se = AsnFind ("Bioseq-set.seq-set.E");
  afd.atp_inst = AsnFind ("Bioseq.inst");
  afd.atp_bsc = AsnFind ("Bioseq-set.class");
  afd.bssp_atp = AsnLinkType (NULL, afd.atp_bss);
  afd.atp_insd = AsnLinkType (NULL, AsnFind ("INSDSet"));
  afd.atp_insde = AsnLinkType (NULL, AsnFind ("INSDSet.E"));
  afd.atp_tss = AsnLinkType (NULL, AsnFind ("TSeqSet"));
  afd.atp_tsse = AsnLinkType (NULL, AsnFind ("TSeqSet.E"));

  /* open output files */

  switch (afd.format) {
    case FLATFILE_FORMAT :
    case FASTA_FORMAT :
    case CDS_FORMAT :
    case GENE_FORMAT :
    case TABLE_FORMAT :
      if (! StringHasNoText (ntout)) {
        afd.nt = FileOpen (ntout, "w");
        if (afd.nt == NULL) {
          Message (MSG_FATAL, "Unable to open nucleotide output file");
          return 1;
        }
      }
      if (! StringHasNoText (aaout)) {
        afd.aa = FileOpen (aaout, "w");
        if (afd.aa == NULL) {
          Message (MSG_FATAL, "Unable to open protein output file");
          return 1;
        }
      }
      break;
    case TINY_FORMAT :
    case INSDSEQ_FORMAT :
      if (! StringHasNoText (ntout)) {
        afd.an = AsnIoOpen (ntout, "wx");
        if (afd.an == NULL) {
          Message (MSG_FATAL, "Unable to open nucleotide output file");
          return 1;
        }
      }
      if (! StringHasNoText (aaout)) {
        afd.ap = AsnIoOpen (aaout, "wx");
        if (afd.ap == NULL) {
          Message (MSG_FATAL, "Unable to open protein output file");
          return 1;
        }
      }
      break;
    case ASN_FORMAT :
      if (! StringHasNoText (ntout)) {
        afd.an = AsnIoOpen (ntout, "w");
        if (afd.an == NULL) {
          Message (MSG_FATAL, "Unable to open output file");
          return 1;
        }
      }
      break;
    case XML_FORMAT :
      if (! StringHasNoText (ntout)) {
        afd.an = AsnIoOpen (ntout, "wx");
        if (afd.an == NULL) {
          Message (MSG_FATAL, "Unable to open output file");
          return 1;
        }
      }
      break;
    case CACHE_COMPONENTS :
      if (! StringHasNoText (ntout)) {
        afd.an = AsnIoOpen (ntout, "wb");
        if (afd.an == NULL) {
          Message (MSG_FATAL, "Unable to open output file");
          return 1;
        }
      }
      break;
    default :
      break;
  }

  /* register fetch functions */

  if (remote) {
    PubSeqFetchEnable ();
    PubMedFetchEnable ();
  }

  if (local) {
    LocalSeqFetchInit (FALSE);
  }

  if (indexed) {
    AsnIndexedLibFetchEnable (asnidx, TRUE);
  }

  /* open output structures */

  switch (afd.format) {
    case TINY_FORMAT :
      if (afd.an != NULL) {
        AsnOpenStruct (afd.an, afd.atp_tss, (Pointer) &(afd.tss));
      }
      if (afd.ap != NULL) {
        AsnOpenStruct (afd.ap, afd.atp_tss, (Pointer) &(afd.tss));
      }
      break;
    case INSDSEQ_FORMAT :
      if (afd.an != NULL) {
        afd.xtran.gbseq  = &(afd.gbsq);
        afd.xtran.aip = afd.an;
        afd.xtran.atp = afd.atp_insde;
        AsnOpenStruct (afd.an, afd.atp_insd, (Pointer) &(afd.gbst));
      }
      if (afd.ap != NULL) {
        afd.xtrap.gbseq  = &(afd.gbsq);
        afd.xtrap.aip = afd.ap;
        afd.xtrap.atp = afd.atp_insde;
        AsnOpenStruct (afd.ap, afd.atp_insd, (Pointer) &(afd.gbst));
      }
      break;
    case CACHE_COMPONENTS :
      if (afd.an != NULL) {
        AsnOpenStruct (afd.an, afd.bssp_atp, (Pointer) &(afd.bss));
        av.intvalue = 7;
        AsnWrite (afd.an, afd.atp_bsc, &av);
        AsnOpenStruct (afd.an, afd.atp_bsss, (Pointer) &(afd.bss.seq_set));
      }
      break;
    default :
      break;
  }

  if (StringDoesHaveText (filterfile)) {
    ReadFilterFile (filterfile, &afd);
  }

  /* process input file or download accession */

  if (StringDoesHaveText (accn)) {

    if (remote) {
      sep = SeqEntryFromAccnOrGi (accn, &afd);
      if (sep != NULL) {
        bsplist = NULL;
        if (afd.lock) {
          bsplist = DoLockFarComponents (sep, afd.useThreads);
        }

        FormatRecord (sep, &afd, bsplist);

        bsplist = UnlockFarComponents (bsplist);

        SeqEntryFree (sep);
      }
    }

  } else if (StringDoesHaveText (directory)) {

    DirExplore (directory, NULL, suffix, TRUE, ProcessOneRecord, (Pointer) &afd);

  } else {

    ProcessOneRecord (asnin, &afd);
  }

  /* close output structures */

  switch (afd.format) {
    case TINY_FORMAT :
      if (afd.an != NULL) {
        AsnCloseStruct (afd.an, afd.atp_tss, NULL);
        AsnPrintNewLine (afd.an);
      }
      if (afd.ap != NULL) {
        AsnCloseStruct (afd.ap, afd.atp_tss, NULL);
        AsnPrintNewLine (afd.ap);
      }
      break;
    case INSDSEQ_FORMAT :
      if (afd.an != NULL) {
        AsnCloseStruct (afd.an, afd.atp_insd, NULL);
        AsnPrintNewLine (afd.an);
      }
      if (afd.ap != NULL) {
        AsnCloseStruct (afd.ap, afd.atp_insd, NULL);
        AsnPrintNewLine (afd.ap);
      }
      break;
    case CACHE_COMPONENTS :
      if (afd.an != NULL) {
        AsnCloseStruct (afd.an, afd.atp_bsss, (Pointer) &(afd.bss.seq_set));
        AsnCloseStruct (afd.an, afd.bssp_atp, (Pointer) &(afd.bss));
        AsnPrintNewLine (afd.an);
      }
      break;
    default :
      break;
  }

  /* close output files */

  if (afd.nt != NULL) {
    FileClose (afd.nt);
  }
  if (afd.aa != NULL) {
    FileClose (afd.aa);
  }

  if (afd.an != NULL) {
    AsnIoClose (afd.an);
  }

  if (afd.ap != NULL) {
    AsnIoClose (afd.ap);
  }

  if (afd.format == CACHE_COMPONENTS) {
    CreateAsnIndex (ntout, NULL, TRUE);
  }

  if (afd.filterList != NULL) {
    ValNodeFreeData (afd.filterList);
  }

  /* close fetch functions */

  if (indexed) {
    AsnIndexedLibFetchDisable ();
  }

  if (local) {
    LocalSeqFetchDisable ();
  }

  if (remote) {
    PubMedFetchDisable ();
    PubSeqFetchDisable ();
  }

  if (afd.failed) {
    return 1;
  }

  return 0;
}

