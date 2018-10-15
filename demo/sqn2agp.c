/*   sqn2agp.c
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
* File Name:  sqn2agp.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   3/14/11
*
* $Revision: 1.12 $
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

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objfdef.h>
#include <objsub.h>
#include <sequtil.h>
#include <gather.h>
#include <sqnutils.h>
#include <explore.h>
#include <gather.h>
#include <seqport.h>
#include <tofasta.h>

#define SQN2AGP_APP_VER "1.5"

CharPtr SQN2AGP_APPLICATION = SQN2AGP_APP_VER;

typedef struct s2aflags {
  CharPtr      results;
  Int2         known;
  Boolean      unknown;
  Boolean      justreport;
  Boolean      contigsqn;
  Boolean      alteredsqn;
  FILE         *afp;
  FILE         *ffp;
  FILE         *rfp;
  ValNodePtr   conlenhead, conlentail;
  ValNodePtr   gaplenhead, gaplentail;
  ValNodePtr   scaffhead, scafftail;
  SeqEntryPtr  sephead, septail;
} S2AFlagData, PNTR S2AFlagPtr;

static BioseqPtr SeqLitToBsp (
  SeqLitPtr lit,
  CharPtr id,
  Int4 seg,
  Uint1 mol
)

{
  BioseqPtr  bsp;
  Char       buf [64];
  SeqLitPtr  slp;

  if (lit == NULL || StringHasNoText (id)) return NULL;

  slp = AsnIoMemCopy (lit, (AsnReadFunc) SeqLitAsnRead, (AsnWriteFunc) SeqLitAsnWrite);
  if (slp == NULL) return NULL;

  bsp = BioseqNew ();
  if (bsp == NULL) return NULL;

  sprintf (buf, "lcl|%s_%ld", id, (long) seg);
  bsp->id = SeqIdParse (buf);
  bsp->repr = Seq_repr_raw;
  bsp->mol = mol;

  bsp->length = slp->length;
  bsp->fuzz = slp->fuzz;
  bsp->seq_data_type = slp->seq_data_type;
  bsp->seq_data = slp->seq_data;

  MemFree (slp); /* fields moved to bioseq, do not free SeqLit components */

  return bsp;
}

static ValNodePtr LIBCALL ValNodeAddIntExAgp (
  ValNodePtr PNTR head,
  ValNodePtr PNTR tail,
  Int2 choice,
  Int4 value
)

{
	ValNodePtr newnode = NULL;
    ValNodePtr vnp;

    newnode = ValNodeNew (NULL);
    if (newnode == NULL) return NULL;

    if (head != NULL) {
        if (*head == NULL) {
            *head = newnode;
        }
    }

    if (tail != NULL) {
        if (*tail != NULL) {
            vnp = *tail;
            while (vnp->next != NULL) {
              vnp = vnp->next;
            }
            vnp->next = newnode;
        }
        *tail = newnode;
    }

	if (newnode != NULL)
	{
		newnode->choice = (Nlm_Uint1)choice;
		newnode->data.intvalue = value;
	}

	return newnode;
}

static Boolean IsStructuredComment (
  UserObjectPtr uop
)

{
  ObjectIdPtr  oip;

  if (uop == NULL) return FALSE;
  oip = uop->type;
  if (oip == NULL) return FALSE;

  if (StringICmp (oip->str, "StructuredComment") == 0) return TRUE;

  return FALSE;
}

static void AddContigToList (
  BioseqPtr bsp,
  SeqLitPtr lit,
  S2AFlagPtr sfp,
  CharPtr id,
  Int4 seg
)

{
  BioSourcePtr   biop;
  CharPtr        comm;
  MolInfoPtr     mip;
  BioseqPtr      newbsp;
  PubdescPtr     pdp;
  SeqEntryPtr    sep;
  UserObjectPtr  uop;
  ValNodePtr     vnp;

  if (bsp == NULL || lit == NULL || sfp == NULL) return;

  newbsp = SeqLitToBsp (lit, id, seg, bsp->mol);
  if (newbsp == NULL) return;
  sep = SeqEntryNew ();
  if (sep == NULL) return;
  sep->choice = 1;
  sep->data.ptrvalue = (Pointer) newbsp;
  if (sfp->sephead == NULL) {
    sfp->sephead = sep;
  }
  if (sfp->septail != NULL) {
    sfp->septail->next = sep;
  }
  sfp->septail = sep;
  vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_source, NULL);
  if (vnp != NULL) {
    if (vnp->data.ptrvalue != NULL) {
      biop = (BioSourcePtr) AsnIoMemCopy (vnp->data.ptrvalue, (AsnReadFunc) BioSourceAsnRead, (AsnWriteFunc) BioSourceAsnWrite);
      if (biop != NULL) {
        SeqDescrAddPointer (&(newbsp->descr), Seq_descr_source, biop);
      }
    }
  }
  vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
  if (vnp != NULL) {
    if (vnp->data.ptrvalue != NULL) {
      mip = (MolInfoPtr) AsnIoMemCopy (vnp->data.ptrvalue, (AsnReadFunc) MolInfoAsnRead, (AsnWriteFunc) MolInfoAsnWrite);
      if (mip != NULL) {
        SeqDescrAddPointer (&(newbsp->descr), Seq_descr_molinfo, mip);
      }
    }
  }
  vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_pub, NULL);
  while (vnp != NULL) {
    if (vnp->data.ptrvalue != NULL) {
      pdp = (PubdescPtr) AsnIoMemCopy (vnp->data.ptrvalue, (AsnReadFunc) PubdescAsnRead, (AsnWriteFunc) PubdescAsnWrite);
      if (pdp != NULL) {
        SeqDescrAddPointer (&(newbsp->descr), Seq_descr_pub, pdp);
      }
    }
    vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_pub, vnp);
  }
  vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_comment, NULL);
  while (vnp != NULL) {
    if (vnp->data.ptrvalue != NULL) {
      comm = StringSave ((CharPtr) vnp->data.ptrvalue);
      if (comm != NULL) {
        SeqDescrAddPointer (&(newbsp->descr), Seq_descr_comment, comm);
      }
    }
    vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_comment, vnp);
  }
  vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_user, NULL);
  while (vnp != NULL) {
    if (vnp->data.ptrvalue != NULL && IsStructuredComment ((UserObjectPtr) vnp->data.ptrvalue)) {
      uop = (UserObjectPtr) AsnIoMemCopy (vnp->data.ptrvalue, (AsnReadFunc) UserObjectAsnRead, (AsnWriteFunc) UserObjectAsnWrite);
      if (uop != NULL) {
        SeqDescrAddPointer (&(newbsp->descr), Seq_descr_user, uop);
      }
    }
    vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_user, vnp);
  }
}

static void DoOneBioseq (
  BioseqPtr bsp,
  Pointer userdata
)

{
  Char         buf [512];
  Int4         cumulative = 0;
  DeltaSeqPtr  dsp;
  Int4         gap_sizes [2];
  Char         id [64];
  Int4         len;
  SeqLitPtr    lit;
  Int4         part = 0;
  Int4         seg = 0;
  S2AFlagPtr   sfp;
  SeqLit       sl;
  Char         space_or_star = ' ';

  if (bsp == NULL) return;
  sfp = (S2AFlagPtr) userdata;
  if (sfp == NULL) return;

  if (sfp->afp == NULL || sfp->ffp == NULL) return;

  if (! ISA_na (bsp->mol)) return;

  if (sfp->unknown) {
    gap_sizes [0] = 100;
  } else {
    gap_sizes [0] = 0;
  }
  gap_sizes [1] = -(sfp->known);

  ConvertNsToGaps (bsp, (Pointer) gap_sizes);

  SeqIdWrite (bsp->id, id, PRINTID_REPORT, sizeof (id) - 1);

  if (bsp->repr == Seq_repr_raw) {

    ValNodeAddIntExAgp (&(sfp->conlenhead), &(sfp->conlentail), 0, bsp->length);

    seg++;
    if (! sfp->justreport) {
      fprintf (sfp->ffp, ">%s_%ld\n", id, (long) seg);
    }

    if (bsp->length < 200) {
      buf [0] = '\0';
      sprintf (buf, "*%s\t%s_%ld\t%ld\t%ld\t%ld", id, id, (long) seg, (long) 1, (long) bsp->length - 1, (long) bsp->length);
      ValNodeCopyStrEx (&(sfp->scaffhead), &(sfp->scafftail), 0, buf);
    }

    if (! sfp->justreport) {
      BioseqFastaStream (bsp, sfp->ffp, 0, 60, 0, 0, FALSE);

      fprintf (sfp->afp, "%s\t1\t%ld\t1\tW\t%s_%ld\t1\t%ld\t+\n", id,
               (long) bsp->length,
               id, (long) seg, (long) bsp->length);
    }

    if (sfp->contigsqn) {
      MemSet ((Pointer) &sl, 0, sizeof (SeqLit));
      sl.length = bsp->length;
      sl.fuzz = bsp->fuzz;
      sl.seq_data_type = bsp->seq_data_type;
      sl.seq_data = bsp->seq_data;
      AddContigToList (bsp, &sl, sfp, id, seg);
    }

    return;
  }

  if (bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4) return;

  for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp != NULL; dsp = dsp->next, cumulative += len) {

    space_or_star = ' ';
    if (dsp == (DeltaSeqPtr) bsp->seq_ext || dsp->next == NULL) {
      space_or_star = '*';
    }

    len = 0;

    if (dsp->choice != 2) continue;
    lit = (SeqLitPtr) dsp->data.ptrvalue;
    if (lit == NULL) continue;

    if (lit->length < 1) {
      /* unknown length */
      continue;
    }

    len = lit->length;

    if (sfp->unknown && len == 100) {

      ValNodeAddIntExAgp (&(sfp->gaplenhead), &(sfp->gaplentail), 0, len);

      /* designated unknown length */
      part++;
      if (! sfp->justreport) {
        fprintf (sfp->afp, "%s\t%ld\t%ld\t%ld\tU\t%ld\tfragment\tyes\t\n", id,
                 (long) cumulative + 1, (long) cumulative + len,
                 (long) part, (long) len);
      }
      continue;
    }

    if (lit->seq_data == NULL || lit->seq_data_type == Seq_code_gap) {

      ValNodeAddIntExAgp (&(sfp->gaplenhead), &(sfp->gaplentail), 0, len);

      /* known length */
      part++;
      if (! sfp->justreport) {
        fprintf (sfp->afp, "%s\t%ld\t%ld\t%ld\tN\t%ld\tfragment\tyes\t\n", id,
                 (long) cumulative + 1, (long) cumulative + len,
                 (long) part, (long) len);
      }
      continue;
    }

    ValNodeAddIntExAgp (&(sfp->conlenhead), &(sfp->conlentail), 0, len);

    seg++;
    if (! sfp->justreport) {
      fprintf (sfp->ffp, ">%s_%ld\n", id, (long) seg);
    }

    if (len < 200) {
      buf [0] = '\0';
      sprintf (buf, "%c%s\t%s_%ld\t%ld\t%ld\t%ld", space_or_star, id, id, (long) seg,
               (long) cumulative + 1, (long) cumulative + len, (long) bsp->length);
      ValNodeCopyStrEx (&(sfp->scaffhead), &(sfp->scafftail), 0, buf);
    }

    part++;
    if (! sfp->justreport) {
      SeqLitFastaStream (lit, sfp->ffp, 0, 60, 0, 0);

      fprintf (sfp->afp, "%s\t%ld\t%ld\t%ld\tW\t%s_%ld\t1\t%ld\t+\n", id,
               (long) cumulative + 1, (long) cumulative + len,
               (long) part, id, (long) seg, (long) len);
    }

    if (sfp->contigsqn) {
      AddContigToList (bsp, lit, sfp, id, seg);
    }
  }
}

static Int4 DoUniqCountInt (
  ValNodePtr list,
  FILE *fp
)

{
  Int4        count = 0;
  ValNodePtr  curr, vnp;
  Int4        shorts = 0;

  if (list == NULL || fp == NULL) return 0;

  curr = list;
  while (curr != NULL) {
    count = 1;
    for (vnp = curr->next; vnp != NULL && curr->data.intvalue == vnp->data.intvalue; vnp = vnp->next) {
      count++;
    }
    fprintf (fp, "%ld\t%ld\n", (long) curr->data.intvalue, (long) count);
    /*
    fprintf (fp, "%6ld\t%6ld\n", (long) curr->data.intvalue, (long) count);
    */
    if (curr->data.intvalue < 200) {
      shorts += count;
    }
    curr = vnp;
    /*
    count++;
    fprintf (fp, "%6ld         %6ld\n", (long) curr->data.intvalue, (long) count);
    curr = curr->next;
    */
  }

  return shorts;
}

static void DoWriteAgpAndFsa (
  SeqEntryPtr sep,
  SeqSubmitPtr ssp,
  Uint2 entityID,
  S2AFlagPtr sfp,
  CharPtr agppath,
  CharPtr altpath,
  CharPtr conpath,
  CharPtr fsapath,
  CharPtr reppath
)

{
  AsnIoPtr      aip;
  BioseqSetPtr  bssp;
  Int4          shorter;
  SeqSubmit     ss;
  CharPtr       str;
  SeqEntryPtr   tmp;
  ValNodePtr    vnp;

  if (sep == NULL || sfp == NULL) return;

  sfp->afp = FileOpen (agppath, "w");
  sfp->ffp = FileOpen (fsapath, "w");
  sfp->rfp = FileOpen (reppath, "w");

  sfp->conlenhead = NULL;
  sfp->conlentail = NULL;
  sfp->gaplenhead = NULL;
  sfp->gaplentail = NULL;
  sfp->scaffhead = NULL;
  sfp->scafftail = NULL;

  sfp->sephead = NULL;
  sfp->septail = NULL;

  VisitBioseqsInSep (sep, (Pointer) sfp, DoOneBioseq);

  if (sfp->contigsqn) {
    tmp = SeqEntryNew ();
    if (tmp != NULL) {
      bssp = BioseqSetNew ();
      if (bssp != NULL) {
        bssp->_class = BioseqseqSet_class_genbank;
        bssp->seq_set = sfp->sephead;
        tmp->choice = 2;
        tmp->data.ptrvalue = (Pointer) bssp;
        aip = AsnIoOpen (conpath, "w");
        if (aip != NULL) {
          if (ssp != NULL) {
            MemSet ((Pointer) &ss, 0, sizeof (SeqSubmit));
            ss.sub = ssp->sub;
            ss.datatype = 1;
            ss.data = (Pointer) tmp;
            SeqSubmitAsnWrite (&ss, aip, NULL);
          } else {
            SeqEntryAsnWrite (tmp, aip, NULL);
          }
          AsnIoClose (aip);
        }
      }
    }
    SeqEntryFree (tmp);
  }

  if (sfp->alteredsqn) {
    aip = AsnIoOpen (altpath, "w");
    if (aip != NULL) {
      if (ssp != NULL) {
        SeqSubmitAsnWrite (ssp, aip, NULL);
      } else {
        SeqEntryAsnWrite (sep, aip, NULL);
      }
      AsnIoClose (aip);
    }
  }

  sfp->conlenhead = ValNodeSort (sfp->conlenhead, SortByIntvalue);
  sfp->gaplenhead = ValNodeSort (sfp->gaplenhead, SortByIntvalue);

  if (sfp->gaplenhead != NULL) {
    fprintf (sfp->rfp, "Ns run length    # of occurrences\n");
    fprintf (sfp->rfp, "-------------    ----------------\n");

    DoUniqCountInt (sfp->gaplenhead, sfp->rfp);

    fprintf (sfp->rfp, "\n\n");
  }

  if (sfp->conlenhead != NULL) {
    fprintf (sfp->rfp, "Contig length    # of occurrences\n");
    fprintf (sfp->rfp, "-------------    ----------------\n");

    shorter = DoUniqCountInt (sfp->conlenhead, sfp->rfp);

    fprintf (sfp->rfp, "\n%ld contigs shorter than 200 nt\n", (long) shorter);

    fprintf (sfp->rfp, "\n\n");
  }

  if (sfp->scaffhead != NULL) {
    fprintf (sfp->rfp, "Scaffold coordinates for contigs shorter than 200 nt\n\n");
    fprintf (sfp->rfp, "Scaffold name   contig name    component start   component end   scaffold length\n");
    fprintf (sfp->rfp, "--------------------------------------------------------------------------------\n");

    for (vnp = sfp->scaffhead; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) return;
      fprintf (sfp->rfp, "%s\n", str);
    }

    fprintf (sfp->rfp, "\n\n");
  }

  sfp->conlenhead = ValNodeFree (sfp->conlenhead);
  sfp->conlentail = NULL;
  sfp->gaplenhead = ValNodeFree (sfp->gaplenhead);
  sfp->gaplentail = NULL;
  sfp->scaffhead = ValNodeFreeData (sfp->scaffhead);
  sfp->scafftail = NULL;

  sfp->sephead = NULL;
  sfp->septail = NULL;

  FileClose (sfp->afp);
  FileClose (sfp->ffp);
  FileClose (sfp->rfp);

  if (sfp->justreport) {
    FileRemove (agppath);
    FileRemove (fsapath);
  }
}

static void ProcessOneRecord (
  CharPtr filename,
  Pointer userdata
)

{
  Char          agppath [PATH_MAX], altpath [PATH_MAX], conpath [PATH_MAX], fsapath [PATH_MAX], reppath [PATH_MAX];
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Pointer       dataptr = NULL;
  Uint2         datatype, entityID = 0;
  FILE          *fp;
  CharPtr       ptr;
  SeqEntryPtr   sep;
  S2AFlagPtr    sfp;
  SeqSubmitPtr  ssp = NULL;

  if (StringHasNoText (filename)) return;
  sfp = (S2AFlagPtr) userdata;
  if (sfp == NULL) return;

  fp = FileOpen (filename, "r");
  if (fp == NULL) {
    Message (MSG_POSTERR, "Failed to open '%s'", filename);
    return;
  }

  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, FALSE, FALSE);

  FileClose (fp);

  entityID = ObjMgrRegister (datatype, dataptr);

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

    if (datatype == OBJ_SEQSUB) {
      ssp = (SeqSubmitPtr) dataptr;
    }

    if (sep != NULL) {

      ptr = StringRChr (filename, '.');
      if (ptr != NULL) {
        *ptr = '\0';
      }

      agppath [0] = '\0';
      altpath [0] = '\0';
      conpath [0] = '\0';
      fsapath [0] = '\0';
      reppath [0] = '\0';
      if (StringDoesHaveText (sfp->results)) {

        StringNCpy_0 (agppath, sfp->results, sizeof (agppath));
        StringNCpy_0 (altpath, sfp->results, sizeof (altpath));
        StringNCpy_0 (conpath, sfp->results, sizeof (conpath));
        StringNCpy_0 (fsapath, sfp->results, sizeof (fsapath));
        StringNCpy_0 (reppath, sfp->results, sizeof (reppath));

        ptr = StringRChr (filename, DIRDELIMCHR);
        if (ptr != NULL) {
          ptr++;
          filename = ptr;
        }
      }
      FileBuildPath (agppath, NULL, filename);
      FileBuildPath (altpath, NULL, filename);
      FileBuildPath (conpath, NULL, filename);
      FileBuildPath (fsapath, NULL, filename);
      FileBuildPath (reppath, NULL, filename);
      StringCat (agppath, ".agp");
      StringCat (altpath, ".alt.asn");
      StringCat (conpath, ".contig.asn");
      StringCat (fsapath, ".contigs.fsa");
      StringCat (reppath, ".rep");

      sep = GetTopSeqEntryForEntityID (entityID);
      if (sep != NULL) {

        AssignIDsInEntity (entityID, 0, NULL);
        DoWriteAgpAndFsa (sep, ssp, entityID, sfp, agppath, altpath, conpath, fsapath, reppath);

      }

      ObjMgrFreeByEntityID (entityID);
    }

  } else {

    Message (MSG_POSTERR, "Datatype %d not recognized", (int) datatype);
  }
}

/* Args structure contains command-line arguments */

typedef enum {
  p_argInputPath = 0,
  r_argOutputPath,
  i_argInputFile,
  f_argFilter,
  x_argSuffix,
  k_argKnown,
  u_argUnknown,
  j_argJustReport,
  c_argSaveContig,
  s_argSaveAltered
} Arguments;

Args myargs [] = {
  {"Path to Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Path for Results", NULL, NULL, NULL,
    TRUE, 'r', ARG_STRING, 0.0, 0, NULL},
  {"Single Input File", "stdin", NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Substring Filter", NULL, NULL, NULL,
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"File Selection Suffix", ".sqn", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Known gap length",  "20", NULL, NULL,
    TRUE, 'k', ARG_INT, 0.0, 0, NULL},
  {"Unknown gap length 100",  "F", NULL, NULL,
    TRUE, 'u', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Just Report",  "F", NULL, NULL,
    TRUE, 'j', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Save Contig Sqn File",  "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Save Altered Sqn File",  "F", NULL, NULL,
    TRUE, 's', ARG_BOOLEAN, 0.0, 0, NULL},
};

Int2 Main (void)

{
  Char         app [64];
  CharPtr      directory, filter, infile, results, suffix;
  S2AFlagData  sfd;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

  /* finish resolving internal connections in ASN.1 parse tables */

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

  sprintf (app, "sqn2agp %s", SQN2AGP_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  MemSet ((Pointer) &sfd, 0, sizeof (S2AFlagData));

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  results = (CharPtr) myargs [r_argOutputPath].strvalue;
  if (StringHasNoText (results)) {
    results = directory;
  }
  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  filter = (CharPtr) myargs [f_argFilter].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;

  sfd.known = (Int2) myargs [k_argKnown].intvalue;
  sfd.unknown = (Boolean) myargs [u_argUnknown].intvalue;
  sfd.justreport = (Boolean) myargs [j_argJustReport].intvalue;
  sfd.contigsqn = (Boolean) myargs [c_argSaveContig].intvalue;
  sfd.alteredsqn = (Boolean) myargs [s_argSaveAltered].intvalue;
  if (sfd.contigsqn) {
    sfd.alteredsqn = TRUE;
  }

  if (StringDoesHaveText (directory)) {

    sfd.results = results;
    DirExplore (directory, filter, suffix, FALSE, ProcessOneRecord, (Pointer) &sfd);

  } else if (StringDoesHaveText (infile)) {

    sfd.results = results;
    ProcessOneRecord (infile, (Pointer) &sfd);
  }

  return 0;
}

