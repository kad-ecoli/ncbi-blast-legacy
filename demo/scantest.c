/*   scantest.c
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
* File Name:  scantest.c
*
* Author:  Kans
*
* Version Creation Date:   1/20/95
*
* $Revision: 6.61 $
*
* File Description: 
*       template for custom scans of ASN.1 release files
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
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <explore.h>
#include <gather.h>
#include <toasn3.h>
#include <tofasta.h>
#include <asn2gnbi.h>

/*
#define DEBUG_SCANTEST
*/

/* UTILITIES */

static ByteStorePtr Se2Bs (
  SeqEntryPtr sep
)

{
  AsnIoBSPtr    aibp;
  ByteStorePtr  bs;

  if (sep == NULL) return NULL;

  bs = BSNew (1000);
  if (bs == NULL) return NULL;
  aibp = AsnIoBSOpen ("w", bs);
  if (aibp == NULL) return NULL;

  SeqEntryAsnWrite (sep, aibp->aip, NULL);

  AsnIoFlush (aibp->aip);
  AsnIoBSClose (aibp);

  return bs;
}

/*
static CharPtr Se2Str (
  SeqEntryPtr sep
)

{
  AsnIoBSPtr    aibp;
  ByteStorePtr  bs;
  CharPtr       str;

  if (sep == NULL) return NULL;

  bs = BSNew (1000);
  if (bs == NULL) return NULL;
  aibp = AsnIoBSOpen ("w", bs);
  if (aibp == NULL) return NULL;

  SeqEntryAsnWrite (sep, aibp->aip, NULL);

  AsnIoFlush (aibp->aip);
  AsnIoBSClose (aibp);

  str = BSMerge (bs, NULL);
  BSFree (bs);

  return str;
}
*/

static TNlmMutex  print_report_mutex = NULL;

static void TSPrintLine (
  FILE *fp,
  CharPtr prefix,
  CharPtr str,
  CharPtr suffix,
  CharPtr final,
  CharPtr delim
)

{
  if (fp == NULL || str == NULL) return;

  if (NlmMutexLockEx (&print_report_mutex)) {
    ErrPostEx (SEV_FATAL, 0, 0, "TSPrintLine failed");
    return;
  }

  if (StringDoesHaveText (prefix)) {
    fprintf (fp, "%s", prefix);
    if (delim != NULL) {
      fprintf (fp, "%s", delim);
    }
  }
  if (StringDoesHaveText (str)) {
    fprintf (fp, "%s", str);
  }
  if (StringDoesHaveText (suffix)) {
    if (delim != NULL) {
      fprintf (fp, "%s", delim);
    }
    fprintf (fp, "%s", suffix);
  }
  if (StringDoesHaveText (final)) {
    if (delim != NULL) {
      fprintf (fp, "%s", delim);
    }
    fprintf (fp, "%s", final);
  }
  fprintf (fp, "\n");
  fflush (fp);

  NlmMutexUnlock (print_report_mutex);
}

static void PrintFirstSeqId (
  SeqEntryPtr sep,
  CharPtr buf,
  size_t buflen,
  Uint2Ptr eIdP,
  BoolPtr isRefSeqP
)

{
  BioseqPtr    fbsp;
  SeqEntryPtr  fsep;
  SeqIdPtr     sip, siphead;

  if (eIdP != NULL) {
    *eIdP = 0;
  }
  if (isRefSeqP != NULL) {
    *isRefSeqP = FALSE;
  }
  if (buf == NULL || buflen < 1) return;
  *buf = '\0';
  if (sep == NULL) return;

  fsep = FindNthBioseq (sep, 1);
  if (fsep == NULL) return;
  fbsp = (BioseqPtr) fsep->data.ptrvalue;
  if (fbsp == NULL) return;

  siphead = SeqIdSetDup (fbsp->id);
  for (sip = siphead; sip != NULL; sip = sip->next) {
    SeqIdStripLocus (sip);
    if (sip->choice == SEQID_OTHER) {
      if (isRefSeqP != NULL) {
        *isRefSeqP = TRUE;
      }
    }
  }
  SeqIdWrite (siphead, buf, PRINTID_FASTA_LONG, buflen);
  SeqIdSetFree (siphead);

  if (eIdP != NULL) {
    *eIdP = fbsp->idx.entityID;
  }
}

/* REPORT SECTION */

typedef struct thrddata {
  Boolean      log;
  Boolean      verbose;
  Boolean      normal;
  Boolean      extended;
  FILE         *fp;
  Char         id [64];
  Boolean      is_refseq;
  SeqEntryPtr  top;
} ThrdData, PNTR ThrdDataPtr;

typedef struct chgdata {
  Boolean      rubisco;
  Boolean      rubiscoL;
  Boolean      rubiscoS;
  Boolean      rbc;
  Boolean      rbcL;
  Boolean      rbcS;
  Boolean      its;
  Boolean      sgml;
  Boolean      miscrna;
  Boolean      msrna;
  Boolean      ncrna;
  Boolean      tmrna;
  Boolean      rnaother;
  Boolean      trnanote;
  Boolean      oldbiomol;
  Boolean      badname;
  Boolean      hasLarge;
  Boolean      hasSmall;
  Boolean      strucSpace;
  Boolean      badDbxref;
  Boolean      refDbxref;
  Boolean      srcDbxref;
  Boolean      capDbxref;
  Boolean      oldDbxref;
  Boolean      privDbxref;
  Boolean      multDbxref;
  Boolean      rareDbxref;
  Boolean      badOrg;
  Int4         protdesc;
  Int4         sfpnote;
  Int4         gbsource;
  Int4         cdsconf;
  Int4         cdscodon;
  ThrdDataPtr  tdp;
} ChangeData, PNTR ChangeDataPtr;



static Boolean IsITS (
  CharPtr name
)

{
  return (StringICmp (name, "its1") == 0 ||
          StringICmp (name, "its 1") == 0 ||
          StringICmp (name, "its2") == 0 ||
          StringICmp (name, "its 2") == 0 ||
          StringICmp (name, "its3") == 0 ||
          StringICmp (name, "its 3") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 1") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 2") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 3") == 0 ||
          StringICmp (name, "internal transcribed spacer 1 (ITS1)") == 0 ||
          StringICmp (name, "internal transcribed spacer 2 (ITS2)") == 0 ||
          StringICmp (name, "internal transcribed spacer 3 (ITS3)") == 0);
}

static Boolean HasSgml (
  CharPtr str,
  ThrdDataPtr tdp
)

{
  Int2  ascii_len;
  Char  buf [1024];

  if (StringHasNoText (str) || tdp == NULL) return FALSE;

  ascii_len = Sgml2AsciiLen (str);
  if (ascii_len + 2 > sizeof (buf)) return FALSE;

  Sgml2Ascii (str, buf, ascii_len + 1);
  if (StringCmp (str, buf) != 0) {
    if (tdp->verbose) {
      TSPrintLine (tdp->fp, "GML", tdp->id, str, NULL, "\t");
    }
    return TRUE;
  }

  return FALSE;
}

static void ScoreFeature (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  CharPtr        comment;
  CdRegionPtr    crp;
  CharPtr        desc;
  GBQualPtr      gbq;
  GeneRefPtr     grp;
  CharPtr        name;
  ProtRefPtr     prp;
  Uint1          residue;
  RnaRefPtr      rrp;
  CharPtr        str;
  ThrdDataPtr    tdp;
  ValNodePtr     vnp;

  if (sfp == NULL) return;
   cdp = (ChangeDataPtr) userdata;
   if (cdp == NULL) return;
  tdp = cdp->tdp;
  if (tdp == NULL || tdp->fp == NULL) return;

  comment = sfp->comment;
  if (StringDoesHaveText (comment)) {
    (cdp->sfpnote)++;
  }

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
    case SEQFEAT_GENE:
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (HasSgml (grp->locus, tdp)) {
        cdp->sgml = TRUE;
      }
      if (HasSgml (grp->desc, tdp)) {
        cdp->sgml = TRUE;
      }
      for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (HasSgml (str, tdp)) {
          cdp->sgml = TRUE;
        }
      }
      break;
    case SEQFEAT_CDREGION:
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp->conflict) {
        (cdp->cdsconf)++;
      }
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "codon") != 0) continue;
        if (StringHasNoText (gbq->val)) continue;
        (cdp->cdscodon)++;
        if (tdp->verbose) {
          TSPrintLine (tdp->fp, "CDN", tdp->id, gbq->val, NULL, "\t");
        }
      }
      break;
    case SEQFEAT_PROT:
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      desc = prp->desc;
      if (StringDoesHaveText (desc)) {
        (cdp->protdesc)++;
      }
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (StringICmp (str, "rubisco large subunit") == 0) {
          cdp->rubiscoL = TRUE;
        }
        if (StringICmp (str, "rubisco small subunit") == 0) {
          cdp->rubiscoS = TRUE;
        }
        if (StringICmp (str, "rubisco") == 0) {
          cdp->rubisco = TRUE;
        }
        if (StringICmp (str, "RbcL") == 0) {
          cdp->rbcL = TRUE;
        }
        if (StringICmp (str, "RbcS") == 0) {
          cdp->rbcS = TRUE;
        }
        if (StringICmp (str, "Rbc") == 0) {
          cdp->rbc = TRUE;
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp->type == 255 && rrp->ext.choice == 1) {
        name = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringCmp (name, "misc_RNA") == 0) {
          cdp->miscrna = TRUE;
          for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
            if (StringCmp (gbq->qual, "product") != 0) continue;
            name = gbq->val;
            if (StringHasNoText (name)) continue;
            if (IsITS (name)) {
              cdp->its = TRUE;
            }
          }
        } else if (StringCmp (name, "miscRNA") == 0) {
          cdp->msrna = TRUE;
        } else if (StringCmp (name, "ncRNA") == 0) {
          cdp->ncrna = TRUE;
        } else if (StringCmp (name, "tmRNA") == 0) {
          cdp->tmrna = TRUE;
        } else {
          cdp->rnaother = TRUE;
          if (IsITS (name)) {
            cdp->its = TRUE;
          }
        }
      } else if (rrp->type == 3 && rrp->ext.choice == 2) {
        if (StringDoesHaveText (comment)) {
          if (StringNCmp (comment, "aa: ", 4) == 0) {
            comment += 4;
          }
          residue = FindTrnaAA3 (comment);
          if (residue > 0 && residue != 255) {
            cdp->trnanote = TRUE;
            if (tdp->verbose) {
              TSPrintLine (tdp->fp, "TR3", tdp->id, comment, NULL, "\t");
            }
          }
          residue = FindTrnaAA (comment);
          if (residue > 0 && residue != 255) {
            cdp->trnanote = TRUE;
            if (tdp->verbose) {
              TSPrintLine (tdp->fp, "TR1", tdp->id, comment, NULL, "\t");
            }
          }
        }
      }
      break;
    default:
      break;
  }
}

static void DoKeywords (
  ChangeDataPtr cdp,
  ValNodePtr keywords
)

{
  CharPtr     str;
  ValNodePtr  vnp;

  if (cdp == NULL || keywords == NULL) return;

  for (vnp = keywords; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    if (StringStr (str, "rbcL gene") != NULL || StringStr (str, "large subunit") != NULL) {
      cdp->hasLarge = TRUE;
    }
    if (StringStr (str, "rbcS gene") != NULL || StringStr (str, "small subunit") != NULL) {
      cdp->hasSmall = TRUE;
    }
  }
}

static void ScoreDescriptor (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  EMBLBlockPtr   ebp;
  GBBlockPtr     gbp;
  MolInfoPtr     mip;

  if (sdp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;

  switch (sdp->choice) {
    case Seq_descr_genbank :
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      if (gbp != NULL) {
        if (StringDoesHaveText (gbp->source)) {
          (cdp->gbsource)++;
        }
        DoKeywords (cdp, gbp->keywords);
      }
      break;
    case Seq_descr_embl :
      ebp = (EMBLBlockPtr) sdp->data.ptrvalue;
      if (ebp != NULL) {
        DoKeywords (cdp, ebp->keywords);
      }
      break;
    case Seq_descr_molinfo :
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL) {
        switch (mip->biomol) {
          case MOLECULE_TYPE_SNRNA:
          case MOLECULE_TYPE_SCRNA:
          case MOLECULE_TYPE_SNORNA:
            cdp->oldbiomol = TRUE;
            break;
          default :
            break;
        }
      }
      break;
    default :
      break;
  }
}

static void CheckForChanges (
  SeqEntryPtr sep,
  ChangeDataPtr cdp
)

{
  if (sep == NULL || cdp == NULL) return;

  VisitFeaturesInSep (sep, (Pointer) cdp, ScoreFeature);
  VisitDescriptorsInSep (sep, (Pointer) cdp, ScoreDescriptor);
}

static void ModGenes (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ModernizeGeneFields (sfp);
}

static void ModRNAs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ModernizeRNAFields (sfp);
}

static void ModPCRs (
  BioSourcePtr biop,
  Pointer userdata
)

{
  BoolPtr         namP;
  PCRPrimerPtr    ppp;
  PCRReactionPtr  prp;

  if (biop == NULL) return;

  ModernizePCRPrimers (biop);

  namP = (BoolPtr) userdata;
  if (namP == NULL) return;

  for (prp = biop->pcr_primers; prp != NULL; prp = prp->next) {
    if (prp->forward == NULL || prp->reverse == NULL) {
      *namP = TRUE;
      return;
    }
    for (ppp = prp->forward; ppp != NULL; ppp = ppp->next) {
      if (StringHasNoText (ppp->seq) && StringDoesHaveText (ppp->name)) {
        *namP = TRUE;
        return;
      }
    }
    for (ppp = prp->reverse; ppp != NULL; ppp = ppp->next) {
      if (StringHasNoText (ppp->seq) && StringDoesHaveText (ppp->name)) {
        *namP = TRUE;
        return;
      }
    }
  }
}

typedef struct bsp2cds {
  BioseqPtr   bsp;
  SeqFeatPtr  cds;
} Bsp2Cds, PNTR Bsp2CdsPtr;

static void ParentCDSProc (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  Bsp2CdsPtr  bcp;
  BioseqPtr   bsp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;
  if (sfp->product == NULL) return;
  bcp = (Bsp2CdsPtr) userdata;
  if (bcp == NULL) return;

  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp != NULL) {
    if (bsp == bcp->bsp) {
      bcp->cds = sfp;
    }
  }
}

static SeqFeatPtr FindParentCDS (
  BioseqPtr bsp,
  ThrdDataPtr tdp
)

{
  Bsp2Cds  b2c;

  if (bsp == NULL || tdp == NULL || tdp->top == NULL) return NULL;

  b2c.bsp = bsp;
  b2c.cds = NULL;

  VisitFeaturesInSep (tdp->top, (Pointer) &b2c, ParentCDSProc);

  return b2c.cds;
}

static void TestForRubisco (
  SeqFeatPtr sfp,
  CharPtr str,
  ChangeDataPtr cdp
)

{
  BioseqPtr    bsp;
  SeqFeatPtr   cds;
  ThrdDataPtr  tdp;

  if (StringHasNoText (str)) return;
  if (cdp == NULL) return;
  tdp = cdp->tdp;
  if (tdp == NULL || tdp->fp == NULL) return;

  if (StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit") == 0) return;
  if (StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase small subunit") == 0) return;
  if (StringStr (str, "ribulose") == NULL || StringStr (str, "bisphosphate") == NULL) return;

  if (StringStr (str, "methyltransferase") != NULL) {
    if (tdp->verbose) {
      TSPrintLine (tdp->fp, "METH", tdp->id, str, NULL, "\t");
    } else {
      TSPrintLine (tdp->fp, "METH", tdp->id, NULL, NULL, " ");
    }
    return;
  }

  if (StringStr (str, "activase") != NULL) {
    if (tdp->verbose) {
      TSPrintLine (tdp->fp, "ACTV", tdp->id, str, NULL, "\t");
    } else {
      TSPrintLine (tdp->fp, "ACTV", tdp->id, NULL, NULL, " ");
    }
    /*
    return;
    */
  }

  if (StringStr (str, "methyltransferase") == NULL) {
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
      if (tdp->verbose) {
        TSPrintLine (tdp->fp, "RIBBIS", tdp->id, str, NULL, "\t");
      } else {
        TSPrintLine (tdp->fp, "RIBBIS", tdp->id, NULL, NULL, " ");
      }
      return;
    }
  }

  if (StringStr (str, "large") != NULL && StringStr (str, "small") == NULL) {
    if (tdp->verbose) {
      TSPrintLine (tdp->fp, "RIBLRG", tdp->id, str, NULL, "\t");
    } else {
      TSPrintLine (tdp->fp, "RIBLRG", tdp->id, NULL, NULL, " ");
    }
    return;
  }

  if (StringStr (str, "small") != NULL && StringStr (str, "large") == NULL) {
    if (tdp->verbose) {
      TSPrintLine (tdp->fp, "RIBSML", tdp->id, str, NULL, "\t");
    } else {
      TSPrintLine (tdp->fp, "RIBSML", tdp->id, NULL, NULL, " ");
    }
    return;
  }

  if (sfp != NULL && sfp->location != NULL) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      cds = FindParentCDS (bsp, tdp);
      if (cds != NULL && StringDoesHaveText (cds->comment)) {
        if (StringStr (cds->comment, "large") != NULL && StringStr (cds->comment, "small") == NULL) {
          if (tdp->verbose) {
            TSPrintLine (tdp->fp, "CDSLRG", tdp->id, str, NULL, "\t");
          } else {
            TSPrintLine (tdp->fp, "CDSLRG", tdp->id, NULL, NULL, " ");
          }
          return;
        }
        if (StringStr (cds->comment, "small") != NULL && StringStr (cds->comment, "large") == NULL) {
          if (tdp->verbose) {
            TSPrintLine (tdp->fp, "CDSSML", tdp->id, str, NULL, "\t");
          } else {
            TSPrintLine (tdp->fp, "CDSSML", tdp->id, NULL, NULL, " ");
          }
          return;
        }
      }
    }
  }

  if (cdp->hasSmall && cdp->hasLarge) {
    if (tdp->verbose) {
      TSPrintLine (tdp->fp, "RIBAMB", tdp->id, str, NULL, "\t");
    } else {
      TSPrintLine (tdp->fp, "RIBAMB", tdp->id, NULL, NULL, " ");
    }
    return;
  }

  if (cdp->hasLarge && (! cdp->hasSmall)) {
    if (tdp->verbose) {
      TSPrintLine (tdp->fp, "KEYLRG", tdp->id, str, NULL, "\t");
    } else {
      TSPrintLine (tdp->fp, "KEYLRG", tdp->id, NULL, NULL, " ");
    }
    return;
  }

  if (cdp->hasSmall && (! cdp->hasLarge)) {
    if (tdp->verbose) {
      TSPrintLine (tdp->fp, "KEYSML", tdp->id, str, NULL, "\t");
    } else {
      TSPrintLine (tdp->fp, "KEYSML", tdp->id, NULL, NULL, " ");
    }
    return;
  }

  if (tdp->verbose) {
    TSPrintLine (tdp->fp, "RIBREM", tdp->id, str, NULL, "\t");
  } else {
    TSPrintLine (tdp->fp, "RIBREM", tdp->id, NULL, NULL, " ");
  }
}

static void TrailingCommaFix (
  CharPtr str,
  ThrdDataPtr tdp,
  CharPtr prefix
)

{
  Char    ch;
  size_t  len;

  if (StringHasNoText (str)) return;
  len = StringLen (str);
  if (len < 1) return;
  ch = str [len - 1];
  while (ch == ' ' && len > 2) {
    len--;
    ch = str [len - 1];
  }
  if (ch == ',') {
    if (tdp != NULL && tdp->verbose && tdp->fp != NULL) {
      str [len] = '\0';
      if (StringHasNoText (prefix)) {
        prefix = "?";
      }
      TSPrintLine (tdp->fp, prefix, tdp->id, str, NULL, "\t");
    }
    str [len - 1] = '_';
    str [len] = '\0';
  }
}

static void FindCommaInGene (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  FILE           *fp;
  GeneRefPtr     grp;
  CharPtr        str;
  ThrdDataPtr    tdp;
  ValNodePtr     vnp;

  if (sfp == NULL) return;
  if (sfp->data.choice != SEQFEAT_GENE) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  tdp = cdp->tdp;
  if (tdp == NULL) return;
  fp = tdp->fp;
  if (fp == NULL) return;
  if (! tdp->verbose) return;

  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return;
  str = grp->locus;
  if (StringDoesHaveText (str)) {
    if (StringChr (str, ',') != NULL) {
      TSPrintLine (tdp->fp, "LOCCOM", tdp->id, str, NULL, "\t");
    }
    if (StringChr (str, ';') != NULL) {
      TSPrintLine (tdp->fp, "LOCSEM", tdp->id, str, NULL, "\t");
    }
  }
  for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    if (StringChr (str, ',') != NULL) {
      TSPrintLine (tdp->fp, "SYNCOM", tdp->id, str, NULL, "\t");
    }
    if (StringChr (str, ';') != NULL) {
      TSPrintLine (tdp->fp, "SYNSEM", tdp->id, str, NULL, "\t");
    }
  }
}

static void FindMuidCitations (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  Boolean        has_muid_cit = FALSE;
  ValNodePtr     ppr;
  ValNodePtr     psp;
  ThrdDataPtr    tdp;
  ValNodePtr     vnp;

  if (sfp == NULL) return;
  psp = sfp->cit;
  if (psp == NULL || psp->data.ptrvalue == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  tdp = cdp->tdp;
  if (tdp == NULL) return;
  if (tdp->fp == NULL) return;
  if (! tdp->verbose) return;

  for (ppr = (ValNodePtr) psp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
    if (ppr->choice == PUB_Muid) {
      has_muid_cit = TRUE;
    } else if (ppr->choice == PUB_Equiv) {
      for (vnp = (ValNodePtr) ppr->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == PUB_Equiv) {
          has_muid_cit = TRUE;
        }
      }
    }
  }

  if (has_muid_cit) {
    TSPrintLine (tdp->fp, "CITMUID", tdp->id, NULL, NULL, "\t");
  }
}

static void FindWholeGraphLocs (
  SeqGraphPtr sgp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  SeqLocPtr      slp;
  ThrdDataPtr    tdp;

  if (sgp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  tdp = cdp->tdp;
  if (tdp == NULL) return;
  if (tdp->fp == NULL) return;

  slp = sgp->loc;
  if (slp == NULL) {
    if (tdp->verbose) {
      TSPrintLine (tdp->fp, "GPHLOC", tdp->id, NULL, NULL, "\t");
    } else {
      TSPrintLine (tdp->fp, "GPHLOC", tdp->id, NULL, NULL, " ");
    }
  } else if (slp->choice == SEQLOC_WHOLE) {
    if (tdp->verbose) {
      TSPrintLine (tdp->fp, "GPHWHL", tdp->id, NULL, NULL, "\t");
    } else {
      TSPrintLine (tdp->fp, "GPHWHL", tdp->id, NULL, NULL, " ");
    }
  }
}

static void RnaProtCmntTrailingCommaFix (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  ProtRefPtr     prp;
  RnaRefPtr      rrp;
  CharPtr        str;
  ThrdDataPtr    tdp;
  ValNodePtr     vnp;

  if (sfp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  tdp = cdp->tdp;
  if (tdp == NULL) return;

  str = sfp->comment;
  if (StringDoesHaveText (str)) {
    TrailingCommaFix (str, tdp, "SFPCOMM");
  }

  if (sfp->data.choice == SEQFEAT_PROT) {
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    /* turn trailing space into trailing underscore for validator */
    for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      TrailingCommaFix (str, tdp, "PRTCOMM");
      TestForRubisco (sfp, str, cdp);
    }
  } else if (sfp->data.choice == SEQFEAT_RNA) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    /* turn trailing space into trailing underscore for validator */
    if (rrp->ext.choice == 1) {
      str = rrp->ext.value.ptrvalue;
      if (StringDoesHaveText (str)) {
        TrailingCommaFix (str, tdp, "RNACOMM");
      }
    }
  }
}

static void ReportObsoleteDbxref (
  ValNodePtr list,
  ChangeDataPtr cdp,
  ThrdDataPtr tdp
)

{
  DbtagPtr     dp;
  ObjectIdPtr  oip;
  CharPtr      str;
  ValNodePtr   vnp;

  if (list == NULL || cdp == NULL || tdp == NULL) return;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    dp = (DbtagPtr) vnp->data.ptrvalue;
    if (dp != NULL && StringDoesHaveText (dp->db)) {
      str = dp->db;
      if (StringICmp (str, "PID") == 0 ||
          StringICmp (str, "PIDg") == 0 ||
          StringICmp (str, "PIDd") == 0 ||
          StringICmp (str, "PIDe") == 0 ||
          StringICmp (str, "NID") == 0 ||
          StringICmp (str, "GI") == 0) {
        cdp->privDbxref = TRUE;
        if (tdp->verbose) {
          TSPrintLine (tdp->fp, "PRVDBX", tdp->id, str, NULL, "\t");
        }
      }
      if (StringICmp (str, "SWISS-PROT") == 0 ||
          StringICmp (str, "SWISSPROT") == 0 ||
          StringICmp (str, "SPTREMBL") == 0 ||
          StringICmp (str, "SUBTILIS") == 0 ||
          StringICmp (str, "MGD") == 0 ||
          StringCmp (str, "cdd") == 0 ||
          StringICmp (str, "TrEMBL") == 0 ||
          StringICmp (str, "LocusID") == 0 ||
          StringICmp (str, "MaizeDB") == 0 ||
          StringICmp (str, "UniProt/Swiss-Prot") == 0 ||
          StringICmp (str, "UniProt/TrEMBL") == 0 ||
          StringICmp (str, "Genew") == 0 ||
          StringICmp (str, "GENEDB") == 0 ||
          StringICmp (str, "GreengenesID") == 0 ||
          StringICmp (str, "HMPID") == 0 ||
          StringICmp (str, "IFO") == 0 ||
          StringICmp (str, "BHB") == 0 ||
          StringICmp (str, "BioHealthBase") == 0) {
        cdp->oldDbxref = TRUE;
        if (tdp->verbose) {
          TSPrintLine (tdp->fp, "OLDDBX", tdp->id, str, NULL, "\t");
        }
      }
      if (StringICmp (str, "ATCC(dna)") == 0 ||
          StringICmp (str, "ATCC(in host)") == 0 ||
          StringICmp (str, "BDGP_EST") == 0 ||
          StringICmp (str, "BDGP_INS") == 0 ||
          StringICmp (str, "CGNC") == 0 ||
          StringICmp (str, "CloneID") == 0 ||
          StringICmp (str, "ENSEMBL") == 0 ||
          StringICmp (str, "ESTLIB") == 0 ||
          StringICmp (str, "GDB") == 0 ||
          StringICmp (str, "GOA") == 0 ||
          StringICmp (str, "IMGT/HLA") == 0 ||
          StringICmp (str, "PIR") == 0 ||
          StringICmp (str, "PSEUDO") == 0 ||
          StringICmp (str, "RZPD") == 0 ||
          StringICmp (str, "SoyBase") == 0 ||
          StringICmp (str, "UNILIB") == 0) {
        cdp->rareDbxref = TRUE;
        if (tdp->verbose) {
          TSPrintLine (tdp->fp, "RARDBX", tdp->id, str, NULL, "\t");
        }
      }
      if (StringICmp (str, "MGD") == 0 || StringICmp (str, "MGI") == 0) {
        oip = dp->tag;
        if (oip != NULL && StringDoesHaveText (oip->str)) {
          str = oip->str;
          if (StringNICmp (str, "MGI:", 4) == 0 || StringNICmp (str, "MGD:", 4) == 0) {
            cdp->oldDbxref = TRUE;
            if (tdp->verbose) {
              TSPrintLine (tdp->fp, "OLDDBX", tdp->id, str, NULL, "\t");
            }
          }
        }
      } else if (StringICmp (str, "HPRD") == 0) {
        oip = dp->tag;
        if (oip != NULL && StringDoesHaveText (oip->str)) {
          str = oip->str;
          if (StringNICmp (str, "HPRD_", 5) == 0) {
            cdp->oldDbxref = TRUE;
            if (tdp->verbose) {
              TSPrintLine (tdp->fp, "OLDDBX", tdp->id, str, NULL, "\t");
            }
          }
        }
      }
      oip = dp->tag;
      if (oip != NULL && StringDoesHaveText (oip->str)) {
        if (StringChr (oip->str, ':') != NULL) {
          cdp->multDbxref = TRUE;
          if (tdp->verbose) {
            TSPrintLine (tdp->fp, "MLTDBX", tdp->id, dp->db, oip->str, "\t");
          }
        }
      }
    }
  }
}

static void ReportInvalidDbxref (
  ValNodePtr list,
  ChangeDataPtr cdp,
  ThrdDataPtr tdp,
  Boolean is_source
)

{
  Boolean     cap;
  DbtagPtr    dp;
  CharPtr     good;
  Boolean     ref;
  Boolean     src;
  CharPtr     str;
  ValNodePtr  vnp;

  if (list == NULL || cdp == NULL || tdp == NULL) return;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    dp = (DbtagPtr) vnp->data.ptrvalue;
    if (dp != NULL && StringDoesHaveText (dp->db)) {
      str = dp->db;

      if (is_source && StringCmp (str, "taxon") == 0) continue;

      if (DbxrefIsValid (str, &ref, &src, &cap, &good)) {
        if (ref && (! tdp->is_refseq)) {
          cdp->refDbxref = TRUE;
          if (tdp->verbose) {
            TSPrintLine (tdp->fp, "REFDBX", tdp->id, str, NULL, "\t");
          }
        }
        if (is_source && (! src)) {
          cdp->srcDbxref = TRUE;
          if (tdp->verbose) {
            TSPrintLine (tdp->fp, "SRCDBX", tdp->id, str, NULL, "\t");
          }
        }
        if (cap) {
          cdp->capDbxref = TRUE;
          if (tdp->verbose) {
            TSPrintLine (tdp->fp, "CAPDBX", tdp->id, str, NULL, "\t");
          }
        }
      } else {
        cdp->badDbxref = TRUE;
        if (tdp->verbose) {
          TSPrintLine (tdp->fp, "BADDBX", tdp->id, str, NULL, "\t");
        }
      }
    }
  }
}

static void LookForObsoleteFeatDbxref (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioSourcePtr   biop;
  ChangeDataPtr  cdp;
  GeneRefPtr     grp;
  OrgRefPtr      orp = NULL;
  ProtRefPtr     prp;
  ThrdDataPtr    tdp;

  if (sfp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  tdp = cdp->tdp;
  if (tdp == NULL) return;

  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (grp != NULL) {
        ReportObsoleteDbxref (grp->db, cdp, tdp);
        ReportInvalidDbxref (grp->db, cdp, tdp, FALSE);
      }
      break;
    case SEQFEAT_PROT :
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      if (prp != NULL) {
        ReportObsoleteDbxref (prp->db, cdp, tdp);
        ReportInvalidDbxref (prp->db, cdp, tdp, FALSE);
      }
      break;
    case SEQFEAT_ORG :
      orp = (OrgRefPtr) sfp->data.value.ptrvalue;
      cdp->badOrg = TRUE;
      break;
    case SEQFEAT_BIOSRC :
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
      }
    default :
      break;
  }

  if (orp != NULL) {
    ReportObsoleteDbxref (orp->db, cdp, tdp);
    ReportInvalidDbxref (orp->db, cdp, tdp, TRUE);
  }

  ReportObsoleteDbxref (sfp->dbxref, cdp, tdp);
  ReportInvalidDbxref (sfp->dbxref, cdp, tdp, FALSE);
}

static void LookForObsoleteDescDbxref (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  BioSourcePtr   biop;
  ChangeDataPtr  cdp;
  OrgRefPtr      orp = NULL;
  ThrdDataPtr    tdp;

  if (sdp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  tdp = cdp->tdp;
  if (tdp == NULL) return;

  switch (sdp->choice) {
    case Seq_descr_org :
      orp = (OrgRefPtr) sdp->data.ptrvalue;
      cdp->badOrg = TRUE;
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

  if (orp != NULL) {
    ReportObsoleteDbxref (orp->db, cdp, tdp);
    ReportInvalidDbxref (orp->db, cdp, tdp, TRUE);
  }
}

static void LookForSemicolonedStrains (
  BioSourcePtr biop,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  OrgModPtr      omp;
  OrgNamePtr     onp;
  OrgRefPtr      orp;
  CharPtr        str;
  ThrdDataPtr    tdp;

  if (biop == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  tdp = cdp->tdp;
  if (tdp == NULL) return;

  orp = biop->org;
  if (orp == NULL) return;
  onp = orp->orgname;
  if (onp == NULL) return;

  for (omp = onp->mod; omp != NULL; omp = omp->next) {
    if (omp->subtype != ORGMOD_strain) continue;
    str = omp->subname;
    if (StringHasNoText (str)) continue;
    if (StringChr (str, ';') == NULL) continue;
    if (tdp->verbose) {
      TSPrintLine (tdp->fp, "STR", tdp->id, str, NULL, "\t");
    } else {
      TSPrintLine (tdp->fp, "STR", tdp->id, NULL, NULL, " ");
    }
  }
}

static void LookForSemicolonedVouchers (
  BioSourcePtr biop,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  OrgModPtr      omp;
  OrgNamePtr     onp;
  OrgRefPtr      orp;
  CharPtr        str;
  ThrdDataPtr    tdp;

  if (biop == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  tdp = cdp->tdp;
  if (tdp == NULL) return;


  orp = biop->org;
  if (orp == NULL) return;
  onp = orp->orgname;
  if (onp == NULL) return;

  for (omp = onp->mod; omp != NULL; omp = omp->next) {
    if (omp->subtype != ORGMOD_specimen_voucher) continue;
    str = omp->subname;
    if (StringHasNoText (str)) continue;
    if (StringChr (str, ';') == NULL) continue;
    if (tdp->verbose) {
      TSPrintLine (tdp->fp, "VCH", tdp->id, str, NULL, "\t");
    } else {
      TSPrintLine (tdp->fp, "VCH", tdp->id, NULL, NULL, " ");
    }
  }
}

static void LookForBadAuth (
  NameStdPtr nsp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  Char           ch;
  Int2           i;
  Boolean        is_bad = FALSE;
  CharPtr        str;
  ThrdDataPtr    tdp;

  if (nsp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  tdp = cdp->tdp;
  if (tdp == NULL) return;

  for (i = 0; i < 6; i++) {
    str = nsp->names [i];
    if (StringHasNoText (str)) continue;
    ch = *str;
    while (ch != '\0') {
      if (IS_DIGIT (ch)) {
        cdp->badname = TRUE;
        is_bad = TRUE;
      }
      str++;
      ch = *str;
    }
  }

  if (is_bad && tdp->fp != NULL && tdp->verbose) {
    TSPrintLine (tdp->fp, "AUTHOR", tdp->id, NULL, NULL, "\t");
    for (i = 0; i < 6; i++) {
      str = nsp->names [i];
      if (StringHasNoText (str)) continue;
      TSPrintLine (tdp->fp, "PERSON", tdp->id, str, NULL, "\t");
    }
  }
}

static void LookForBadPub (
  PubdescPtr pdp,
  Pointer userdata
)

{
  VisitAuthorsInPub (pdp, userdata, LookForBadAuth);
}

static void LookForStrucComment (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  ObjectIdPtr    oip;
  CharPtr        str;
  ThrdDataPtr    tdp;
  UserFieldPtr   ufp;
  UserObjectPtr  uop;

  if (sdp == NULL || sdp->choice != Seq_descr_user) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  tdp = cdp->tdp;
  if (tdp == NULL) return;

  if (tdp->fp == NULL) return;

  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "StructuredComment") != 0) return;

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    if (ufp->choice != 1) continue;
    oip = ufp->label;
    if (oip == NULL) continue;
    str = oip->str;
    if (StringHasNoText (str)) continue;
    if (StringStr (str, " ") != NULL) {
      cdp->strucSpace = TRUE;
      if (tdp->verbose) {
        TSPrintLine (tdp->fp, "STCCOM", tdp->id, str, NULL, "\t");
      }
    }
  }
}

static void CommentDescrTrailingCommaFix (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  CharPtr        str;
  ThrdDataPtr    tdp;

  if (sdp == NULL || sdp->choice != Seq_descr_comment) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  tdp = cdp->tdp;
  if (tdp == NULL) return;

  str = (CharPtr) sdp->data.ptrvalue;
  if (StringDoesHaveText (str)) {
    TrailingCommaFix (str, tdp, "DSCCOMM");
  }
}

static void LookForGeneOverlapChange (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  SeqFeatPtr     gene1, gene2;
  GeneRefPtr     grp;
  ThrdDataPtr    tdp;

  if (sfp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  tdp = cdp->tdp;
  if (tdp == NULL) return;

  if (sfp->data.choice != SEQFEAT_CDREGION) return;

  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL) return;

  gene1 = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_GENE, NULL, 0, NULL, CONTAINED_WITHIN, NULL);
  gene2 = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_GENE, NULL, 0, NULL, LOCATION_SUBSET, NULL);

  if (gene1 != gene2) {
    if (tdp->verbose) {
      TSPrintLine (tdp->fp, "OVL", tdp->id, NULL, NULL, "\t");
    } else {
      TSPrintLine (tdp->fp, "OVL", tdp->id, NULL, NULL, " ");
    }
  }
}

static void StripBadProtTitles (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CharPtr            buf;
  size_t             buflen = 1001;
  ObjValNodePtr      ovp;
  SeqIdPtr           sip;
  CharPtr            title;
  ValNodePtr         vnp;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) return;
  }

  vnp = BioseqGetSeqDescr (bsp, Seq_descr_title, NULL);
  if (vnp == NULL) return;
  title = (CharPtr) vnp->data.ptrvalue;
  if (StringHasNoText (title)) return;

  buf = MemNew (sizeof (Char) * (buflen + 1));
  if (buf == NULL) return;

  if (NewCreateDefLineBuf (NULL, bsp, buf, buflen, TRUE, FALSE)) {
    if (StringICmp (buf, title) != 0) {
      if (vnp->extended != 0) {
        ovp = (ObjValNodePtr) vnp;
        ovp->idx.deleteme = TRUE;
      }
    }
  }

  MemFree (buf);
}

static void BadProtTitleProc (
  SeqEntryPtr sep,
  Pointer mydata,
  Int4 index,
  Int2 indent
)

{
  BioseqSetPtr  bssp;

  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp->_class != BioseqseqSet_class_nuc_prot) return;
  VisitBioseqsInSep (sep, NULL, StripBadProtTitles);
}

static void DoReport (
  SeqEntryPtr sep,
  ThrdDataPtr tdp
)

{
  ByteStorePtr  bs = NULL, tmp = NULL;
  Boolean       bsec = FALSE, cma = FALSE, norm = FALSE, ttl = FALSE, ssec = FALSE;
  Boolean       gen = FALSE, ncr = FALSE, pcr = FALSE, nam = FALSE, gbsrc = FALSE;
  ChangeData    cdbefore, cdafter;
  Uint2         entityID;

  if (sep == NULL || tdp == NULL) return;

  RemoveAllNcbiCleanupUserObjects (sep);

  MemSet ((Pointer) &cdbefore, 0, sizeof (ChangeData));
  MemSet ((Pointer) &cdafter, 0, sizeof (ChangeData));

  cdbefore.tdp = tdp;
  cdafter.tdp = tdp;

  CheckForChanges (sep, &cdbefore);

  if (tdp->extended) {
    bs = Se2Bs (sep);
    NormalizeDescriptorOrder (sep);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      norm = TRUE;
    }
    BSFree (bs);
    bs = tmp;
  }

  VisitFeaturesInSep (sep, (Pointer) &cdbefore, RnaProtCmntTrailingCommaFix);
  VisitDescriptorsInSep (sep, (Pointer) &cdbefore, CommentDescrTrailingCommaFix);
  VisitDescriptorsInSep (sep, (Pointer) &cdbefore, LookForStrucComment);
  VisitFeaturesInSep (sep, (Pointer) &cdbefore, LookForObsoleteFeatDbxref);
  VisitDescriptorsInSep (sep, (Pointer) &cdbefore, LookForObsoleteDescDbxref);
  VisitBioSourcesInSep (sep, (Pointer) &cdbefore, LookForSemicolonedStrains);
  VisitBioSourcesInSep (sep, (Pointer) &cdbefore, LookForSemicolonedVouchers);
  VisitFeaturesInSep (sep, (Pointer) &cdbefore, FindCommaInGene);
  VisitFeaturesInSep (sep, (Pointer) &cdbefore, FindMuidCitations);
  VisitGraphsInSep (sep, (Pointer) &cdbefore, FindWholeGraphLocs);
  VisitFeaturesInSep (sep, (Pointer) &cdbefore, LookForGeneOverlapChange);

  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    cma = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  BasicSeqEntryCleanup (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    bsec = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  EntryChangeGBSource (sep);
  EntryCheckGBBlock (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    gbsrc = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitPubdescsInSep (sep, (Pointer) &cdbefore, LookForBadPub);

  if (tdp->extended) {
    VisitFeaturesInSep (sep, NULL, ModGenes);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      gen = TRUE;
    }
    BSFree (bs);
    bs = tmp;

    VisitFeaturesInSep (sep, NULL, ModRNAs);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      ncr = TRUE;
    }
    BSFree (bs);
    bs = tmp;

    VisitBioSourcesInSep (sep, (Pointer) &nam, ModPCRs);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      pcr = TRUE;
    }
    BSFree (bs);
    bs = tmp;

    entityID = ObjMgrGetEntityIDForChoice (sep);
    SeqMgrIndexFeatures (entityID, NULL);
    SeqEntryExplore (sep, NULL, BadProtTitleProc);
    DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
    SeqMgrIndexFeatures (entityID, NULL);
    InstantiateProteinTitles (entityID, NULL);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    BasicSeqEntryCleanup (sep);
    NormalizeDescriptorOrder (sep);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      ttl = TRUE;
    }
    BSFree (bs);
    bs = tmp;

    SeriousSeqEntryCleanup (sep, NULL, NULL);

    RemoveAllNcbiCleanupUserObjects (sep);

    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      ssec = TRUE;
    }
    BSFree (bs);
    bs = tmp;
  }

  CheckForChanges (sep, &cdafter);

  BSFree (bs);

  if (tdp->extended) {
    if (ssec) {
      TSPrintLine (tdp->fp, "SSEC", tdp->id, NULL, NULL, " ");
    } else if (ttl) {
      TSPrintLine (tdp->fp, "TITL", tdp->id, NULL, NULL, " ");
    } else if (gbsrc) {
      TSPrintLine (tdp->fp, "GBSR", tdp->id, NULL, NULL, " ");
    } else if (bsec) {
      TSPrintLine (tdp->fp, "BSEC", tdp->id, NULL, NULL, " ");
    } else if (norm) {
      TSPrintLine (tdp->fp, "NORM", tdp->id, NULL, NULL, " ");
    } else {
      /*
      TSPrintLine (tdp->fp, "OKAY", tdp->id, NULL, NULL, " ");
      */
    }
  }

  if (cma) {
    TSPrintLine (tdp->fp, "CMA", tdp->id, NULL, NULL, " ");
  }

  if (tdp->extended) {
    if (gen) {
      TSPrintLine (tdp->fp, "GEN", tdp->id, NULL, NULL, " ");
    }
    if (ncr) {
      TSPrintLine (tdp->fp, "MDR", tdp->id, NULL, NULL, " ");
    }
    if (pcr) {
      TSPrintLine (tdp->fp, "PCR", tdp->id, NULL, NULL, " ");
    }
    if (nam) {
      TSPrintLine (tdp->fp, "NAM", tdp->id, NULL, NULL, " ");
    }
  }

  if (cdbefore.rubiscoL) {
    TSPrintLine (tdp->fp, "RUL", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.rubiscoS) {
    TSPrintLine (tdp->fp, "RUS", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.rubisco) {
    TSPrintLine (tdp->fp, "RUB", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.rbcL) {
    TSPrintLine (tdp->fp, "RBL", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.rbcS) {
    TSPrintLine (tdp->fp, "RBS", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.rbc) {
    TSPrintLine (tdp->fp, "RBC", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.its) {
    TSPrintLine (tdp->fp, "ITS", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.sgml) {
    TSPrintLine (tdp->fp, "SGM", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.rnaother) {
    TSPrintLine (tdp->fp, "RNA", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.miscrna) {
    TSPrintLine (tdp->fp, "MCR", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.msrna) {
    TSPrintLine (tdp->fp, "MSR", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.ncrna) {
    TSPrintLine (tdp->fp, "NCR", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.tmrna) {
    TSPrintLine (tdp->fp, "TMR", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.trnanote) {
    TSPrintLine (tdp->fp, "TRN", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.oldbiomol) {
    TSPrintLine (tdp->fp, "MOL", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.badname) {
    TSPrintLine (tdp->fp, "AUT", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.strucSpace) {
    TSPrintLine (tdp->fp, "SPC", tdp->id, NULL, NULL, " ");
  }
  if (! tdp->verbose) {
    if (cdbefore.badDbxref) {
      TSPrintLine (tdp->fp, "BDBX", tdp->id, NULL, NULL, " ");
    }
    if (cdbefore.refDbxref) {
      TSPrintLine (tdp->fp, "RDBX", tdp->id, NULL, NULL, " ");
    }
    if (cdbefore.srcDbxref) {
      TSPrintLine (tdp->fp, "SDBX", tdp->id, NULL, NULL, " ");
    }
    if (cdbefore.capDbxref) {
      TSPrintLine (tdp->fp, "CDBX", tdp->id, NULL, NULL, " ");
    }
    if (cdbefore.privDbxref) {
      TSPrintLine (tdp->fp, "PDBX", tdp->id, NULL, NULL, " ");
    }
    if (cdbefore.oldDbxref) {
      TSPrintLine (tdp->fp, "ODBX", tdp->id, NULL, NULL, " ");
    }
    if (cdbefore.multDbxref) {
      TSPrintLine (tdp->fp, "MDBX", tdp->id, NULL, NULL, " ");
    }
    if (cdbefore.rareDbxref) {
      TSPrintLine (tdp->fp, "RDBX", tdp->id, NULL, NULL, " ");
    }
    if (cdbefore.cdscodon > 0) {
    TSPrintLine (tdp->fp, "CDN", tdp->id, NULL, NULL, " ");
    }
  }
  if (cdbefore.badOrg) {
    TSPrintLine (tdp->fp, "ORG", tdp->id, NULL, NULL, " ");
  }

  if (cdbefore.protdesc != cdafter.protdesc) {
    TSPrintLine (tdp->fp, "PRT", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.sfpnote != cdafter.sfpnote) {
    TSPrintLine (tdp->fp, "COM", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.gbsource != cdafter.gbsource) {
    TSPrintLine (tdp->fp, "SRC", tdp->id, NULL, NULL, " ");
  }
  if (cdbefore.cdsconf != cdafter.cdsconf) {
    TSPrintLine (tdp->fp, "CNF", tdp->id, NULL, NULL, " ");
  }
}

static void DoFlatfile (
  SeqEntryPtr sep,
  ThrdDataPtr tdp
)

{
#ifdef OS_UNIX
  Char    buf [256];
  Char    cmmd [256];
  size_t  ct;
  FILE    *fp;
  Char    path1 [PATH_MAX];
  Char    path2 [PATH_MAX];
  Char    path3 [PATH_MAX];

  if (sep == NULL || tdp == NULL || tdp->fp == NULL) return;

  if (NlmMutexLockEx (&print_report_mutex)) {
    ErrPostEx (SEV_FATAL, 0, 0, "DoFlatfile mutex lock failed");
    return;
  }

  fprintf (tdp->fp, "%s\n", tdp->id);
  fflush (tdp->fp);

  TmpNam (path1);
  TmpNam (path2);
  TmpNam (path3);

  fp = FileOpen (path1, "w");
  if (fp != NULL) {
    SeqEntryToGnbk (sep, NULL, GENBANK_FMT, ENTREZ_MODE, NORMAL_STYLE, 1, 0, 0, NULL, fp);
  }
  FileClose (fp);
  /*
  SeriousSeqEntryCleanupBulk (sep);
  */
  fp = FileOpen (path2, "w");
  if (fp != NULL) {
    SeqEntryToGnbk (sep, NULL, GENBANK_FMT, ENTREZ_MODE, NORMAL_STYLE, 262145, 0, 0, NULL, fp);
  }
  FileClose (fp);

  sprintf (cmmd, "diff %s %s > %s", path1, path2, path3);
  system (cmmd);

  sprintf (cmmd, "cat %s", path3);
  fp = popen (cmmd, "r");
  if (fp != NULL) {
    while ((ct = fread (buf, 1, sizeof (buf), fp)) > 0) {
      fwrite (buf, 1, ct, tdp->fp);
      fflush (tdp->fp);
    }
    pclose (fp);
  }

  sprintf (cmmd, "rm %s; rm %s; rm %s", path1, path2, path3);
  system (cmmd);

  NlmMutexUnlock (print_report_mutex);
#endif
}

static void DoTest (
  SeqEntryPtr sep,
  ThrdDataPtr tdp
)

{
  ByteStorePtr  bs, tmp;
  Uint2         entityID;
#ifdef DEBUG_SCANTEST
  CharPtr       final = NULL;
  Char          pfx [32];
  Char          sfx [64];
#endif

  if (sep == NULL || tdp == NULL) return;

  entityID = ObjMgrGetEntityIDForChoice (sep);

#ifdef DEBUG_SCANTEST
  sprintf (pfx, "     %5d", (int) entityID);
  sprintf (sfx, "%ld - %ld", (long) sep, (long) sep->data.ptrvalue);
#endif

  bs = Se2Bs (sep);
  NormalizeDescriptorOrder (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
#ifdef DEBUG_SCANTEST
    final = "*";
#endif
  }
  BSFree (bs);
  BSFree (tmp);

#ifdef DEBUG_SCANTEST
  TSPrintLine (tdp->fp, pfx, tdp->id, sfx, final, " ");
#endif

  SeriousSeqEntryCleanupBulk (sep);
}

/* CONSUMER CALLBACK */

typedef struct userflags {
  Boolean  log;
  Boolean  verbose;
  Boolean  debug;
  Boolean  normal;
  Boolean  extended;
  Boolean  flatfile;
  FILE     *fp;
} UserFlagData, PNTR UserFlagPtr;

static void ProcessOneRecord (
  SeqEntryPtr sep,
  Pointer userdata
)

{
  ThrdData     thd;
  UserFlagPtr  ufp;

  if (sep == NULL) return;
  ufp = (UserFlagPtr) userdata;
  if (ufp == NULL) return;

  if (ufp->fp == NULL) return;

  /* ufp structure is common to all consumer threads */

  /* thd structure is specific to each consumer thread */

  MemSet ((Pointer) &thd, 0, sizeof (ThrdData));

  thd.log = ufp->log;
  thd.verbose = ufp->verbose;
  thd.normal = ufp->normal;
  thd.extended = ufp->extended;
  thd.fp = ufp->fp;

  /* assign IDs needed for GetNextDescriptorUnindexed */

  AssignIDsInEntityEx (0, OBJ_SEQENTRY, (Pointer) sep, NULL);

  /* print first SeqID and set is_refseq flag */

  PrintFirstSeqId (sep, thd.id, sizeof (thd.id), NULL, &(thd.is_refseq));

  if (ufp->log) {
    TSPrintLine (ufp->fp, "LOG", thd.id, NULL, NULL, " ");
  }

  thd.top = sep;

  /* process the record */

  if (ufp->debug) {
    DoTest (sep, &thd);
    return;
  }

  if (ufp->flatfile) {
    DoFlatfile (sep, &thd);
    return;
  }

  if (ufp->normal || ufp->extended) {
    DoReport (sep, &thd);
    return;
  }
}

/* PRODUCER SECTION */

typedef void (*ConsumerFunc) (SeqEntryPtr sep, Pointer userdata);

typedef struct appflags {
  CharPtr            directory;
  CharPtr            infile;
  CharPtr            filter;
  CharPtr            suffix;
  Boolean            dorecurse;
  Boolean            binary;
  Boolean            compressed;
  Boolean            useThreads;
  Int2               numConsumers;
  FILE               *fp;
  ScanBioseqSetFunc  scanproc;
  ConsumerFunc       consumerproc;
  VoidPtr            consumerdata;
} AppFlagData, PNTR AppFlagPtr;

/* PRODUCER-CONSUMER DATA BUFFER */

#define DATABUFSIZE 16

static VoidPtr  dataBuf [DATABUFSIZE];

static int  producer_position = 0;
static int  consumer_position = 0;

/* SEMAPHORE AND MUTEX PROTECTIONS */

static TNlmSemaphore  producer_semaphore = NULL;
static TNlmSemaphore  consumer_semaphore = NULL;
static TNlmMutex      buffer_mutex = NULL;

/* THREAD-SAFE ACCESS FUNCTIONS */

static void QueueItem (
  VoidPtr data
)

{
  NlmSemaWait (producer_semaphore);
  NlmMutexLock (buffer_mutex);

  dataBuf [producer_position] = data;
  producer_position = (producer_position + 1) % DATABUFSIZE;

  NlmMutexUnlock (buffer_mutex);
  NlmSemaPost (consumer_semaphore);
}

static VoidPtr FetchItem (
  void
)

{
  VoidPtr  data = NULL;

  NlmSemaWait (consumer_semaphore);
  NlmMutexLock (buffer_mutex);

  data = dataBuf [consumer_position];
  consumer_position = (consumer_position + 1) % DATABUFSIZE;

  NlmMutexUnlock (buffer_mutex);
  NlmSemaPost (producer_semaphore);

  return data;
}

/* CONSUMER CALLING FUNCTIONS */

static VoidPtr DoConsume (
  VoidPtr arg
)

{
  AppFlagPtr    afp;
  ConsumerFunc  consumerproc;
  SeqEntryPtr   sep;
#ifdef DEBUG_SCANTEST
  Uint2         entityID;
  Char          id [64];
  Char          pfx [32];
  Char          sfx [64];
#endif

  afp = (AppFlagPtr) arg;
  if (afp == NULL) return NULL;
  consumerproc = afp->consumerproc;
  if (consumerproc == NULL) return NULL;

  /* threads consume data until poison pill is encountered */

  sep = (SeqEntryPtr) FetchItem ();

  while (sep != NULL) {

    consumerproc (sep, afp->consumerdata);

#ifdef DEBUG_SCANTEST
    PrintFirstSeqId (sep, id, sizeof (id), &entityID, NULL);

    sprintf (pfx, "%5d", (int) entityID);
    sprintf (sfx, "%ld - %ld", (long) sep, (long) sep->data.ptrvalue);
    TSPrintLine (afp->fp, "FREE", pfx, id, sfx, " ");
#endif

    FreeScanSeqEntryMT (sep);

#ifdef DEBUG_SCANTEST
    TSPrintLine (afp->fp, "DONE", pfx, id, sfx, " ");
#endif

    /*
    omp = ObjMgrGet ();
    ObjMgrReapOne (omp);
    SeqMgrClearBioseqIndex ();
    ObjMgrFreeCache (0);
    FreeSeqIdGiCache ();
    */

    sep = (SeqEntryPtr) FetchItem ();
  }

  return NULL;
}

static void QueueForLater (
  SeqEntryPtr sep,
  Pointer userdata
)

{
  if (sep == NULL) return;

  /* threaded version, queue for processing by consumer */

  QueueItem ((Pointer) sep);
}

static void ConsumeImmediately (
  SeqEntryPtr sep,
  Pointer userdata
)

{
  AppFlagPtr    afp;
  ConsumerFunc  consumerproc;

  if (sep == NULL) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;
  consumerproc = afp->consumerproc;
  if (consumerproc == NULL) return;

  /* non-threaded version, process immediately */

  consumerproc (sep, afp->consumerdata);
  SeqEntryFree (sep);
}

/* RELEASE FILE PROCESSING */

static void DoOneReleaseFile (
  CharPtr filename,
  Pointer userdata
)

{
  AppFlagPtr  afp;

  if (StringHasNoText (filename)) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  if (StringStr (filename, "gbcon") != NULL ||
      StringStr (filename, "gbest") != NULL ||
      StringStr (filename, "gbgss") != NULL ||
      StringStr (filename, "gbhtg") != NULL ||
      StringStr (filename, "gbsts") != NULL) {
    if (afp->fp != NULL) {
      TSPrintLine (afp->fp, "SKIP", filename, NULL, NULL, " ");
    }
#ifdef DEBUG_SCANTEST
    printf ("SKIP %s\n", filename);
    fflush (stdout);
#endif
    return;
  }

  if (afp->fp != NULL) {
    TSPrintLine (afp->fp, "FILE", filename, NULL, NULL, " ");
  }

#ifdef DEBUG_SCANTEST
  printf ("FILE %s\n", filename);
  fflush (stdout);
#endif

  ScanBioseqSetReleaseMT (filename, afp->binary, afp->compressed,
                          (Pointer) afp, afp->scanproc);
}

/* DIRECTORY EXPLORATION */

static VoidPtr DoProduce (
  VoidPtr arg
)

{
  AppFlagPtr  afp;

  afp = (AppFlagPtr) arg;
  if (afp == NULL) return NULL;

  if (StringDoesHaveText (afp->directory)) {

    DirExplore (afp->directory, afp->filter, afp->suffix, afp->dorecurse,
                DoOneReleaseFile, (Pointer) afp);

  } else if (StringDoesHaveText (afp->infile)) {

    DoOneReleaseFile (afp->infile, (Pointer) afp);
  }

  return NULL;
}

/* JOB PROCESSING */

#define MAX_PROCESSING_THREADS 8

static void DoJob (
  AppFlagPtr afp,
  ConsumerFunc consumerproc,
  VoidPtr consumerdata
)

{
  TNlmThread  cthds [MAX_PROCESSING_THREADS];
  Int2        i;
  Int4        numThreadsToUse;
  TNlmThread  pthd;
  VoidPtr     status;

  if (afp == NULL) return;

  afp->consumerproc = consumerproc;
  afp->consumerdata = consumerdata;

  if (NlmThreadsAvailable () && afp->useThreads) {

    numThreadsToUse = NlmCPUNumber ();
    /*
    numThreadsToUse = 2 * NlmCPUNumber () - 1;
    */
    if (afp->numConsumers > 0) {
      numThreadsToUse = afp->numConsumers;
    }
    if (numThreadsToUse > MAX_PROCESSING_THREADS) {
      numThreadsToUse = MAX_PROCESSING_THREADS;
    }

    afp->scanproc = QueueForLater;

    /* initialize semaphores and mutex */

    producer_semaphore = NlmSemaInit (DATABUFSIZE);
    consumer_semaphore = NlmSemaInit (0);
    NlmMutexInit (&buffer_mutex);

    pthd = NlmThreadCreate (DoProduce, (Pointer) afp);

    for (i = 0; i < numThreadsToUse; i++) {
      cthds [i] = NlmThreadCreate (DoConsume, (Pointer) afp);
    }

    /* wait for producer to finish */

    NlmThreadJoin (pthd, &status);

    /* once producer is done send one poison pill per consumer thread */

    for (i = 0; i < numThreadsToUse; i++) {
      QueueItem (NULL);
    }

    /* wait for all consumers to finish */

    for (i = 0; i < numThreadsToUse; i++) {
      NlmThreadJoin (cthds [i], &status);
    }

  } else if (afp->useThreads) {

    Message (MSG_POSTERR, "Threads unavailable, continuing...");

    afp->scanproc = ConsumeImmediately;

    DoProduce ((Pointer) afp);

  } else {

    afp->scanproc = ConsumeImmediately;

    DoProduce ((Pointer) afp);
  }
}

/* COMMAND-LINE ARGUMENTS */

#define p_argInputPath    0
#define i_argInputFile    1
#define o_argOutputFile   2
#define f_argFilter       3
#define x_argSuffix       4
#define u_argRecurse      5
#define b_argBinary       6
#define c_argCompressed   7
#define T_argThreads      8
#define C_argConsumers    9
#define l_argLog         10
#define v_argVerbose     11
#define d_argDebug       12
#define n_argNormal      13
#define e_argExtended    14
#define q_argFlatfile    15

Args myargs [] = {
  {"Path to Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Input File Name", NULL, NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Output File Name", NULL, NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Substring Filter", NULL, NULL, NULL,
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"File Selection Suffix", ".aso", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Recurse", "F", NULL, NULL,
    TRUE, 'u', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Bioseq-set is Binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Bioseq-set is Compressed", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Use Threads", "F", NULL, NULL,
    TRUE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Consumer Thread Count", "0", NULL, NULL,
    FALSE, 'C', ARG_INT, 0.0, 0, NULL},
  {"Log SeqID", "F", NULL, NULL,
    TRUE, 'l', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Verbose Output", "F", NULL, NULL,
    TRUE, 'v', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Debugging Tests", "F", NULL, NULL,
    TRUE, 'd', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Normal Tests", "F", NULL, NULL,
    TRUE, 'n', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Extended Tests", "F", NULL, NULL,
    TRUE, 'e', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Flatfile Report", "F", NULL, NULL,
    TRUE, 'q', ARG_BOOLEAN, 0.0, 0, NULL},
};

/* MAIN FUNCTION */

extern Int2 Main (void)

{
  AppFlagData   afd;
  FILE          *fp;
  CharPtr       outfile;
  time_t        runtime;
  time_t        starttime;
  time_t        stoptime;
  UserFlagData  ufd;

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

  if (! GetArgs ("scantest", sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  fp = FileOpen (outfile, "w");
  if (fp == NULL) {
    return 0;
  }

  MemSet ((Pointer) &afd, 0, sizeof (AppFlagData));

  afd.directory = (CharPtr) myargs [p_argInputPath].strvalue;
  afd.infile = (CharPtr) myargs [i_argInputFile].strvalue;
  afd.filter = (CharPtr) myargs [f_argFilter].strvalue;
  afd.suffix = (CharPtr) myargs [x_argSuffix].strvalue;
  afd.dorecurse = (Boolean) myargs [u_argRecurse].intvalue;
  afd.binary = (Boolean) myargs [b_argBinary].intvalue;
  afd.compressed = (Boolean) myargs [c_argCompressed].intvalue;
  afd.useThreads = (Boolean) myargs [T_argThreads].intvalue;
  afd.numConsumers = (Int2) myargs [C_argConsumers].intvalue;
  afd.fp = fp;

  MemSet ((Pointer) &ufd, 0, sizeof (UserFlagData));

  ufd.log = (Boolean) myargs [l_argLog].intvalue;
  ufd.verbose = (Boolean) myargs [v_argVerbose].intvalue;
  ufd.debug = (Boolean) myargs [d_argDebug].intvalue;
  ufd.normal = (Boolean) myargs [n_argNormal].intvalue;
  ufd.extended = (Boolean) myargs [e_argExtended].intvalue;
  ufd.flatfile = (Boolean) myargs [q_argFlatfile].intvalue;
  ufd.fp = fp;

  starttime = GetSecs ();

  /* main processing function */

  DoJob (&afd, ProcessOneRecord, (VoidPtr) &ufd);

  /* clean up */

  stoptime = GetSecs ();
  runtime = stoptime - starttime;

  fprintf (fp, "Finished in %ld seconds\n", (long) runtime);

  FileClose (fp);

  return 0;
}


