/*  tofasta.c
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
* File Name:  tofasta.c
*
* Author:  James Ostell
*
* Version Creation Date: 7/12/91
*
* $Revision: 6.235 $
*
* File Description:  various sequence objects to fasta output
*
* Modifications:
* --------------------------------------------------------------------------
* Date       Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/
#include <tofasta.h>
#include <gather.h>
#include <sqnutils.h>  /* MakeSeqID */
#include <subutil.h>   /* MOLECULE_TYPE_GENOMIC */
#include <explore.h>
#include <objloc.h>
#include <objfdef.h>
#include <asn2gnbi.h>

#ifdef OS_UNIX_DARWIN
#define NLM_GETC fgetc
#else
#define NLM_GETC getc
#endif
#define SeqLocNew(_a) ValNodeNew((_a))

static Uint1 na_order[NUM_SEQID] = {   /* order of nucleic acid deflines */
    255, /* 0 = not set */
    230, /* 1 = local Object-id */
    30,  /* 2 = gibbsq */
    30,  /* 3 = gibbmt */
    255, /* 4 = giim Giimport-id */
    20, /* 5 = genbank */
    20, /* 6 = embl */
    255, /* 7 = pir */
    255, /* 8 = swissprot */
    40,  /* 9 = patent */
    15, /* 10 = other TextSeqId (RefGene) */
    50, /* 11 = general Dbtag */
    120,  /* 12 = gi */
    20, /* 13 = ddbj */
    255, /* 14 = prf */
    30, /* 15 = pdb */
    20,  /* 16 = tpg */
    20,  /* 17 = tpe */
    20,  /* 18 = tpd */
    20,  /* 19 = gpp */
    20   /* 30 = nat */
    };
static Uint1 aa_order[NUM_SEQID] = {   /* order of nucleic acid deflines */
    255, /* 0 = not set */
    230, /* 1 = local Object-id */
    40,  /* 2 = gibbsq */
    40,  /* 3 = gibbmt */
    255, /* 4 = giim Giimport-id */
    60, /* 5 = genbank */
    60, /* 6 = embl */
    30, /* 7 = pir */
    20, /* 8 = swissprot */
    80,  /* 9 = patent */
    15, /* 10 = other TextSeqId (RefGene) */
    90, /* 11 = general Dbtag */
    120,  /* 12 = gi */
    60, /* 13 = ddbj */
    70, /* 14 = prf */
    50, /* 15 = pdb */
    60,  /* 16 = tpg */
    60,  /* 17 = tpe */
    60,  /* 18 = tpd */
    60,  /* 19 = gpp */
    60   /* 20 = nat */
    };
#define FASTA_BUFFER_LEN 524288
#define PATENT_ORDER 110         /* order for any patent */
/*****************************************************************************
*
*   The above sets the ordering to be, lowest to highest
*
Nucleic Acids:
    GenBank/EMBL/DDBJ
    PDB
    Patents
    Anything else
Proteins:
    SWISSPROT
    PIR
    NCBI BackBone (but not in GenBank)
    PDB
    GenBank/EMBL/DDBJ translations
    PRF
    Patents
    Anything else
*
*****************************************************************************/
Int4 GetOrderBySeqId(Int4 choice, Boolean is_prot)
{
    if(choice > NUM_SEQID)
        return -1;
    if(is_prot)
        return aa_order[choice];
    else
        return na_order[choice];
}
/*****************************************************************************
*
*   Traversal routine for SeqEntryToFasta
*
*****************************************************************************/
void SeqEntryFasta (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  FastaPtr tfa;
  BioseqPtr bsp = NULL;
  BioseqSetPtr bssp = NULL;
  MyFsaPtr mfp;
  Boolean is_na;
  SeqIdPtr sip;
  TextSeqIdPtr tsip;
  ValNodePtr vnp;
  OrgRefPtr orp;
  MolInfoPtr mip;
  tfa = (FastaPtr) data;
  mfp = tfa->mfp;
  if (tfa->group_segs == 2)       /* put out only segments */
    {
      if (tfa->parts != -1)      /* in parts set */
        {
          if (indent <= tfa->parts)   /* out of parts set */
            {
              tfa->parts = -1;
              tfa->seg = -1;
            }
        }
    }
  if (IS_Bioseq(sep))
    {
      bsp = (BioseqPtr)(sep->data.ptrvalue);
      vnp = bsp->descr;
    }
  else
    {
      bssp = (BioseqSetPtr)(sep->data.ptrvalue);
      vnp = bssp->descr;
    }
  orp = NULL;
  mip = NULL;
  while (vnp != NULL)   /* check for organism info */
    {
      switch (vnp->choice)
        {
        case Seq_descr_source:
          orp = ((BioSourcePtr)(vnp->data.ptrvalue))->org;
          break;
        case Seq_descr_org:
          orp = (OrgRefPtr)(vnp->data.ptrvalue);
          break;
        case Seq_descr_molinfo:
          mip = (MolInfoPtr)(vnp->data.ptrvalue);
          break;
        default:
          break;
        }
      vnp = vnp->next;
    }
  if (orp != NULL)
    {
      if (orp->taxname != NULL)
        mfp->organism = orp->taxname;
      else if (orp->common != NULL)
        mfp->organism = orp->common;
    }
  if (mip != NULL)
    mfp->tech = mip->tech;
  else
    mfp->tech = 0 ;
  if (! IS_Bioseq(sep))    /* check for taking only parts of seg seqs */
    {
      if (tfa->group_segs == 2)    /* put out only segments */
        {
          if (bssp->_class == 2)   /* segset */
            tfa->seg = indent;
          else if (bssp->_class == 4)   /* parts */
            {
              if ((tfa->seg >= 0) && (tfa->seg < indent))
                {
                  tfa->parts = indent;   /* in parts set */
                }
            }
        }
      return;
    }
  is_na = tfa->is_na;
  if ((! is_na) && (! ISA_aa(bsp->mol))) /* check for translations */
    {
      for (sip = bsp->id; sip != NULL; sip = sip->next)
        {
          switch (sip->choice)
            {
            case SEQID_GENBANK:
            case SEQID_EMBL:
            case SEQID_DDBJ:
            case SEQID_OTHER:
            case SEQID_TPG:
            case SEQID_TPE:
            case SEQID_TPD:
            case SEQID_GPIPE:
              tsip = (TextSeqIdPtr)(sip->data.ptrvalue);
              if (tsip->accession != NULL)
                mfp->accession = tsip->accession;
              break;
            default:
              break;
            }
        }
    }
  if (tfa->last_indent != -1)   /* putting out segments together */
    {
      if (indent > tfa->last_indent)
        return;
      tfa->last_indent = -1;
    }
  /* do raw bioseqs only */
  if (! tfa->group_segs)
    {
      if (BioseqRawToFastaX(bsp, mfp, is_na))
        tfa->got_one = TRUE;
    }
  else if (tfa->group_segs == 1)    /* do segmented sets */
    {
      if (BioseqToFastaX(bsp, mfp, is_na))
        {
          tfa->got_one = TRUE;
          if (bsp->repr == Seq_repr_seg)
            tfa->last_indent = indent;
        }
    }
  else if (tfa->group_segs == 2)    /* take only the parts */
    {
      if (tfa->parts >= 0)    /* in segmented parts set */
        {
          if (BioseqRawToFastaX(bsp, mfp, is_na))
            tfa->got_one = TRUE;
        }
    }
  return;
}
/*****************************************************************************
*
*   SeqEntryToFasta(sep, fp, is_na)
*
*****************************************************************************/
NLM_EXTERN Boolean SeqEntryToFasta (SeqEntryPtr sep, FILE *fp, Boolean is_na)
{
    if (IS_Bioseq(sep))
        return SeqEntrysToFasta(sep, fp, is_na, 3);
    else
        return SeqEntrysToFasta(sep, fp, is_na, 0);
}
static Boolean SeqEntrysToFastaXX (SeqEntryPtr sep, FILE *fp, Boolean is_na, Uint1 group_segs, Boolean printid_general);
NLM_EXTERN Boolean SeqEntryToFastaEx (SeqEntryPtr sep, FILE *fp, Boolean is_na, Boolean printid_general)
{
    if (IS_Bioseq(sep))
        return SeqEntrysToFastaXX(sep, fp, is_na, 3, printid_general);
    else
        return SeqEntrysToFastaXX(sep, fp, is_na, 0, printid_general);
}
/*****************************************************************************
*
*   FastaFileFunc(key, buf, data)
*       standard "write to file" callback
*
*****************************************************************************/
NLM_EXTERN Boolean FastaFileFunc (BioseqPtr bsp, Int2 key, CharPtr buf,
                                  Uint4 buflen, Pointer data)
{
    FILE * fp;
    fp = (FILE *)data;
    switch (key)
    {
        case FASTA_ID:
            fprintf(fp, ">%s ", buf);
            break;
        case FASTA_DEFLINE:
            fprintf(fp, "%s\n", buf);
            break;
        case FASTA_SEQLINE:
            fprintf(fp, "%s\n", buf);
            break;
        case FASTA_EOS:   /* end of sequence */
            break;
        default:
            break;
    }
    return TRUE;
}
/*****************************************************************************
*
*   FastaFileFunc(key, buf, data)
*       standard "write to file" callback
*
*    Used for BLAST (FASTA) databases.  If the defline is
*    longer than buflen, then check that an ID is not
*    truncated in the middle.
*
*****************************************************************************/
NLM_EXTERN Boolean FastaDumpFileFunc (BioseqPtr bsp, Int2 key, CharPtr buf,
                                  Uint4 buflen, Pointer data)
{
    FILE * fp;
    fp = (FILE *)data;
    switch (key)
    {
        case FASTA_ID:
            fprintf(fp, ">%s ", buf);
            break;
        case FASTA_DEFLINE:
            if (buflen >= FASTA_BUFFER_LEN-1)
            {
                Uint4 index=buflen;
                while (index > 0 && buf[index] != ' ')
                {
                    if (buf[index] == '\001')
                    {
                        buf[index] = NULLB;
                        break;
                    }
                    index--;
                }
            }
            fprintf(fp, "%s\n", buf);
            break;
        case FASTA_SEQLINE:
            fprintf(fp, "%s\n", buf);
            break;
        case FASTA_EOS:   /* end of sequence */
            break;
        default:
            break;
    }
    return TRUE;
}
/*****************************************************************************
*
*   SeqEntrysToFasta(sep, fp, is_na, group_segs)
*
*       group_segs = 0 ... take only raw Bioseqs
*       group_segs = 1 ... group segmented seqs into single entry.. no parts
*       group_segs = 2 ... show only parts of segmented seqs
*       group_segs = 3 ... like 1, but instantiate virtual Bioseqs
*
*****************************************************************************/
static Boolean SeqEntrysToFastaXX (SeqEntryPtr sep, FILE *fp, Boolean is_na, Uint1 group_segs, Boolean printid_general)
{
    FastaDat tfa;
    MyFsa mfa;
    Char buf[FASTA_BUFFER_LEN+1];
    if ((sep == NULL) || (fp == NULL))
        return FALSE;
    MemSet ((Pointer) (&mfa), 0, sizeof (MyFsa));
    mfa.buf = buf;
    mfa.buflen = FASTA_BUFFER_LEN;
    mfa.seqlen = 70;
    mfa.mydata = (Pointer)fp;
    mfa.myfunc = FastaFileFunc;
    mfa.bad_asn1 = FALSE;
    mfa.order = 0;
    mfa.accession = NULL;
    mfa.organism = NULL;
    mfa.do_virtual = FALSE;
    mfa.tech = 0;
    mfa.no_sequence = FALSE;
    mfa.formatdb    = FALSE;
    mfa.printid_general = printid_general;
    mfa.seqloc = NULL;
    tfa.mfp = &mfa;
    tfa.is_na = is_na;
    if (is_na)
        mfa.code = Seq_code_iupacna;
    else
        mfa.code = Seq_code_ncbieaa;
    if (group_segs == 3)  /* do 2 things */
    {
        mfa.do_virtual = TRUE;
        group_segs = 1;
    }
    tfa.group_segs = group_segs;
    tfa.last_indent = -1;
    tfa.parts = -1;
    tfa.seg = -1;
    tfa.got_one = FALSE;
    SeqEntryExplore(sep, (Pointer)&tfa, SeqEntryFasta);
    return tfa.got_one;
}
NLM_EXTERN Boolean SeqEntrysToFasta (SeqEntryPtr sep, FILE *fp, Boolean is_na, Uint1 group_segs)
{
    return SeqEntrysToFastaXX (sep, fp, is_na, group_segs, FALSE);
}
/*****************************************************************************
*
*   SeqEntrysToFastaX(sep, mfa, is_na, group_segs)
*
*****************************************************************************/
NLM_EXTERN Boolean SeqEntrysToFastaX (SeqEntryPtr sep, MyFsaPtr mfp, Boolean is_na, Uint1 group_segs)
{
    FastaDat tfa;
    if ((sep == NULL) || (mfp == NULL))
        return FALSE;
    tfa.mfp = mfp;
    tfa.is_na = is_na;
    if (group_segs == 3)  /* do 2 things */
    {
        mfp->do_virtual = TRUE;
        group_segs = 1;
    }
    tfa.group_segs = group_segs;
    tfa.last_indent = -1;
    tfa.parts = -1;
    tfa.seg = -1;
    tfa.got_one = FALSE;
        SeqEntryExplore(sep, (Pointer)&tfa, SeqEntryFasta);
    return tfa.got_one;
}
/*****************************************************************************
*
*   SeqEntrysToDefline(sep, mfa, is_na, group_segs)
*
*****************************************************************************/
#define DEFLINE_MAX_LEN FASTA_BUFFER_LEN
NLM_EXTERN Boolean SeqEntrysToDefline(SeqEntryPtr sep,
                           FILE *fp, Boolean is_na, Uint1 group_segs)
{
  FastaDat tfa;
  MyFsa mfa;
  if ((sep == NULL) || (fp == NULL))
    return FALSE;
  MemSet ((Pointer) (&mfa), 0, sizeof (MyFsa));
  mfa.buf = (CharPtr) MemNew(DEFLINE_MAX_LEN);
  mfa.buflen = DEFLINE_MAX_LEN-1;
  mfa.seqlen = DEFLINE_MAX_LEN;
  mfa.mydata = (Pointer)fp;
  mfa.myfunc = FastaFileFunc;
  mfa.no_sequence = TRUE;
  mfa.bad_asn1 = FALSE;
  mfa.order = 0;
  mfa.accession = NULL;
  mfa.organism = NULL;
  mfa.do_virtual = FALSE;
  mfa.formatdb = FALSE;
  mfa.tech = 0;
  mfa.printid_general = FALSE;
  mfa.seqloc = NULL;
  tfa.mfp = &mfa;
  tfa.is_na = is_na;
  if (group_segs == 3)  /* do 2 things */
    {
      mfa.do_virtual = TRUE;
      group_segs = 1;
    }
  tfa.group_segs = group_segs;
  tfa.last_indent = -1;
  tfa.parts = -1;
  tfa.seg = -1;
  tfa.got_one = FALSE;
  SeqEntryExplore(sep, (Pointer)&tfa, SeqEntryFasta);
  MemFree(mfa.buf);
  return tfa.got_one;
}
/*****************************************************************************
*
*   Boolean BioseqRawToFasta(bsp, fp, is_na)
*
*****************************************************************************/
NLM_EXTERN Boolean BioseqRawToFasta (BioseqPtr bsp, FILE *fp, Boolean is_na)
{
    return BioseqRawToFastaExtra(bsp, fp, 80);
}
NLM_EXTERN Boolean BioseqRawToFastaExtra (BioseqPtr bsp, FILE *fp, Int2 line_length)
{
     return BioseqRawToFastaExtraEx (bsp, fp, line_length, NULL);
}
NLM_EXTERN Boolean BioseqRawToFastaExtraEx(BioseqPtr bsp, FILE *fp, Int2 line_length, SeqLocPtr slp)
{
    MyFsa mfa;
    Char buf[FASTA_BUFFER_LEN+1];
    if ((bsp == NULL) || (fp == NULL))
        return FALSE;
    MemSet ((Pointer) (&mfa), 0, sizeof (MyFsa));
    mfa.buf = buf;
    mfa.buflen = FASTA_BUFFER_LEN;
    mfa.seqlen = line_length;
    mfa.mydata = (Pointer)fp;
    mfa.myfunc = FastaFileFunc;
    mfa.bad_asn1 = FALSE;
    mfa.order = 0;
    mfa.accession = NULL;
    mfa.organism = NULL;
    mfa.do_virtual = FALSE;
    mfa.tech = 0;
    mfa.no_sequence = FALSE;
    mfa.formatdb = FALSE;
    mfa.printid_general = FALSE;
    mfa.seqloc = slp;
     return BioseqRawToFastaX(bsp, &mfa, ISA_na(bsp->mol));
}
/*****************************************************************************
*
*   Boolean BioseqRawToFastaX(bsp, mfp, is_na)
*
*****************************************************************************/
NLM_EXTERN Boolean BioseqRawToFastaX (BioseqPtr bsp, MyFsaPtr mfp, Boolean is_na)
{
    Uint1 repr;
    if ((bsp == NULL) || (mfp == NULL))
        return FALSE;
    repr = Bioseq_repr(bsp);
    if (! ((repr == Seq_repr_raw) || (repr == Seq_repr_const)))
        return FALSE;
    return BioseqToFastaX(bsp, mfp, is_na);
}
/*****************************************************************************
*
*   Boolean BioseqToFasta(bsp, fp, is_na)
*
*****************************************************************************/
NLM_EXTERN Boolean BioseqToFasta (BioseqPtr bsp, FILE *fp, Boolean is_na)
{
    MyFsa mfa;
    Char buf[FASTA_BUFFER_LEN+1];
    if ((bsp == NULL) || (fp == NULL))
        return FALSE;
    MemSet ((Pointer) (&mfa), 0, sizeof (MyFsa));
    mfa.buf = buf;
    mfa.buflen = FASTA_BUFFER_LEN;
    mfa.seqlen = 80;
    mfa.mydata = (Pointer)fp;
    mfa.myfunc = FastaFileFunc;
    mfa.bad_asn1 = FALSE;
    mfa.order = 0;
    mfa.accession = NULL;
    mfa.organism = NULL;
    mfa.do_virtual = FALSE;
    mfa.tech = 0;
    mfa.no_sequence = FALSE;
    mfa.formatdb = FALSE;
    mfa.printid_general = FALSE;
    mfa.seqloc = NULL;
    return BioseqToFastaX(bsp, &mfa, is_na);
}
/*****************************************************************************
*
*   Boolean BioseqToFastaDump(bsp, fp, is_na)
*
*****************************************************************************/
NLM_EXTERN Boolean BioseqToFastaDump (BioseqPtr bsp, FILE *fp, Boolean is_na)
{
    MyFsa mfa;
    Char buf[FASTA_BUFFER_LEN+1];
    if ((bsp == NULL) || (fp == NULL))
        return FALSE;
    MemSet ((Pointer) (&mfa), 0, sizeof (MyFsa));
    mfa.buf = buf;
    mfa.buflen = FASTA_BUFFER_LEN;
    mfa.seqlen = 80;
    mfa.mydata = (Pointer)fp;
    mfa.myfunc = FastaDumpFileFunc;
    mfa.bad_asn1 = FALSE;
    mfa.order = 0;
    mfa.accession = NULL;
    mfa.organism = NULL;
    mfa.do_virtual = FALSE;
    mfa.tech = 0;
    mfa.no_sequence = FALSE;
    mfa.formatdb = FALSE;
    mfa.printid_general = FALSE;
    mfa.seqloc = NULL;
    return BioseqToFastaX(bsp, &mfa, is_na);
}
/*****************************************************************************
*
*   Boolean BioseqToFastaX(bsp, mfp, is_na)
*
*****************************************************************************/
static Boolean FastaIdX(BioseqPtr bsp, CharPtr buf, Uint4 buflen, Boolean printid_general, SeqLocPtr seqloc);
NLM_EXTERN Boolean BioseqToFastaX (BioseqPtr bsp, MyFsaPtr mfp, Boolean is_na)
{
    SeqPortPtr spp;
    Uint1 repr, code;
    Char buf[41];
    SeqIdPtr sip;
    Uint1 order = 255;
    Boolean is_patent = FALSE, is_genbank = FALSE;
    Uint1Ptr order_array;
    int i;
    CharPtr organism = NULL;
    if ((bsp == NULL) || (mfp == NULL))
        return FALSE;
    repr = Bioseq_repr(bsp);
    if (ISA_na(bsp->mol))
    {
        if (! is_na)
            return FALSE;
        order_array = na_order;
    }
    else if (ISA_aa(bsp->mol))
    {
        if (is_na)
            return FALSE;
        order_array = aa_order;
        if (mfp->accession != NULL)           /* translated genbank */
        {
            order = order_array[SEQID_GENBANK];
            is_genbank = TRUE;
            organism = mfp->organism;
        }
    }
    else
    {
        buf[0] = '\0';
          SeqIdWrite(SeqIdFindBest(bsp->id, 0), buf, PRINTID_FASTA_LONG, 40);
          ErrPostEx(SEV_ERROR,0,0,"ToFasta: [%s] Unrecognized bsp->mol = %d",
          buf, (int)(bsp->mol));
    mfp->bad_asn1 = TRUE;
    return FALSE;
    }
    mfp->bsp = bsp;
    for (sip = bsp->id; sip != NULL; sip = sip->next)
    {
        i=(int)(sip->choice);
        if (! is_genbank)    /* don't change order for translated genbank */
        {
            if (order_array[i] < order)
                order = order_array[i];
        }
        if (i == (int)SEQID_PATENT)
            is_patent = TRUE;
        else if (i == (int)SEQID_PRF)
            organism = mfp->organism;
    }
    if (is_patent)
        order = PATENT_ORDER;
    mfp->order = order;
    switch (mfp->tech)
    {
        case MI_TECH_est:
        case MI_TECH_sts:
        case MI_TECH_survey:
        case MI_TECH_htgs_1:
        case MI_TECH_htgs_2:
        case MI_TECH_htgs_3:
            organism = mfp->organism;
            break;
        default:
            break;
    }
    if (! FastaIdX(bsp, mfp->buf, mfp->buflen, mfp->printid_general, mfp->seqloc))
        return FALSE;
    (*(mfp->myfunc))(bsp, FASTA_ID, mfp->buf, StringLen(mfp->buf), mfp->mydata);
       if (! CreateDefLine(NULL, bsp, mfp->buf, mfp->buflen, mfp->tech, mfp->accession, organism))
           return FALSE;
    (*(mfp->myfunc))(bsp, FASTA_DEFLINE, mfp->buf, StringLen(mfp->buf), mfp->mydata);
        if (mfp->formatdb && is_na) {
            (*(mfp->myfunc))(bsp, FASTA_FORMATDB_AMB, mfp->buf, StringLen(mfp->buf), mfp->mydata);
        }
        else if(!mfp->no_sequence) {
        if (!mfp->formatdb) {
        if (is_na)
            code = Seq_code_iupacna;
        else
            code = Seq_code_ncbieaa;
        } else {
        code = mfp->code;
        }
        if (repr == Seq_repr_virtual && (! mfp->do_virtual)) {
            StringCpy (mfp->buf, "-");
                (*(mfp->myfunc))(bsp, FASTA_SEQLINE, mfp->buf, StringLen(mfp->buf),
                                 mfp->mydata);
            (*(mfp->myfunc))(bsp, FASTA_EOS, mfp->buf, StringLen(mfp->buf),
                             mfp->mydata);
            return TRUE;
        }
            spp = FastaSeqPortEx (bsp, is_na, mfp->do_virtual, code, mfp->seqloc);
            if (spp == NULL) return FALSE;
            while (FastaSeqLineEx(spp, mfp->buf, mfp->seqlen, is_na, mfp->do_virtual))
                (*(mfp->myfunc))(bsp, FASTA_SEQLINE, mfp->buf, StringLen(mfp->buf),
                                 mfp->mydata);
            SeqPortFree(spp);
            (*(mfp->myfunc))(bsp, FASTA_EOS, mfp->buf, StringLen(mfp->buf),
                             mfp->mydata);
        }
        return TRUE;
}

/*****************************************************************************
*
*   BioseqFastaStream (bsp, fp, flags, linelen, blocklen, grouplen, do_defline)
*
*       Rapid FASTA generator using SeqPortStream
*
*****************************************************************************/

typedef struct streamfsa {
  FILE          *fp;
  ByteStorePtr  bs;
  Char          buf [512];
  Int2          idx;
  Int2          lin;
  Int2          blk;
  Int2          grp;
  Int2          linelen;
  Int2          blocklen;
  Int2          grouplen;
  Int2          skip;
  Int4          gi;
  Int4          start;
  Int4          seqpos;
  Boolean       seqspans;
} StreamFsa, PNTR StreamFsaPtr;

static void LIBCALLBACK FsaStreamProc (
  CharPtr sequence,
  Pointer userdata
)

{
  Char          ch;
  StreamFsaPtr  sfp;
  Char          spn [64];

  if (StringHasNoText (sequence) || userdata == NULL) return;
  sfp = (StreamFsaPtr) userdata;
  ch = *sequence;
  while (ch != '\0' && sfp->skip > 0) {
    (sfp->skip)--;
    (sfp->seqpos)++;
    sequence++;
    ch = *sequence;
  }
  while (ch != '\0') {
    /* optionally separate blocks with space */
    if (sfp->blk >= sfp->blocklen && sfp->blocklen > 0) {
      sfp->buf [sfp->idx] = ' ';
      (sfp->idx)++;
      sfp->blk = 0;
    }
    /* save sequence character to buffer */
    sfp->buf [sfp->idx] = ch;
    (sfp->idx)++;
    (sfp->lin)++;
    (sfp->blk)++;
    /* write sequence as soon as we have line of characters */
    if (sfp->lin >= sfp->linelen) {
      sfp->buf [sfp->idx] = '\0';
      /* optionally separate groups with blank line */
      if (sfp->grp >= sfp->grouplen && sfp->grouplen > 0) {
        if (sfp->fp != NULL) {
          fprintf (sfp->fp, "\n");
        } else if (sfp->bs != NULL) {
          BSWrite (sfp->bs, "\n", sizeof ("\n"));
        }
        sfp->grp = 0;
      }
      /* print actual sequence line here */
      if (sfp->fp != NULL) {
        if (sfp->seqspans) {
          fprintf (sfp->fp, "<span class=\"ff_line\" id=\"gi_%ld_%ld\">", (long) sfp->gi, (long) (sfp->start + 1));
        }
        fprintf (sfp->fp, "%s", sfp->buf);
        if (sfp->seqspans) {
          fprintf (sfp->fp, "</span>");
        }
        fprintf (sfp->fp, "\n");
      } else if (sfp->bs != NULL) {
        if (sfp->seqspans) {
          sprintf (spn, "<span class=\"ff_line\" id=\"gi_%ld_%ld\">", (long) sfp->gi, (long) (sfp->start + 1));
          BSWrite (sfp->bs, spn, StringLen (spn));
        }
        BSWrite (sfp->bs, sfp->buf, StringLen (sfp->buf));
        if (sfp->seqspans) {
          BSWrite (sfp->bs, "</span>", sizeof ("</span>"));
        }
        BSWrite (sfp->bs, "\n", sizeof ("\n"));
      }
      sfp->start = sfp->seqpos + 1;
      sfp->idx = 0;
      sfp->lin = 0;
      sfp->blk = 0;
      (sfp->grp)++;
    }
    (sfp->seqpos)++;
    sequence++;
    ch = *sequence;
  }
}

/* If Bioseq is a protein, does not have an accession,
 * and is the only protein in the nuc-prot set or is part of
 * a sorted protein file, use the nuc bioseq ID in the FASTA
 * defline.
 */

static SeqIdPtr ChooseFastaID (BioseqPtr bsp, Boolean allow_mult)

{
  BioseqSetPtr bssp;
  BioseqPtr    nuc_bsp = NULL;
  SeqIdPtr     sip;
  if (bsp == NULL) return NULL;
  if (!ISA_aa(bsp->mol) || bsp->idx.parenttype != OBJ_BIOSEQSET || bsp->idx.parentptr == NULL) {
    return bsp->id;
  }
  /* if protein sequence has an accession, do not use nucleotide ID */
  sip = bsp->id;
  while (sip != NULL) {
    if (sip->choice == SEQID_GENBANK) {
      return sip;
    } else {
      sip = sip->next;
    }
  }
  bssp = (BioseqSetPtr) bsp->idx.parentptr;
  if (bssp->_class != BioseqseqSet_class_nuc_prot /* not in nuc-prot set */
      || bssp->seq_set == NULL /* no sequences in set - bad indexing */
      || bssp->seq_set->next == NULL /* only one sequence in nuc-prot set, degenerate */
      || (!allow_mult && bssp->seq_set->next->next != NULL) /* more than one protein in nuc-prot set */) {
    return bsp->id;
  }
  if (IS_Bioseq (bssp->seq_set)) {
    nuc_bsp = bssp->seq_set->data.ptrvalue;
  } else if (IS_Bioseq_set (bssp->seq_set)) {
    bssp = bssp->seq_set->data.ptrvalue;
    if (bssp->_class == BioseqseqSet_class_segset
        && bssp->seq_set != NULL
        && IS_Bioseq (bssp->seq_set)) {
      nuc_bsp = bssp->seq_set->data.ptrvalue;
    }
  }
  if (nuc_bsp == NULL) {
    return bsp->id;
  } else {
    return nuc_bsp->id;
  }
}

static Int4 BioseqFastaStreamInternal (
  BioseqPtr bsp,
  SeqLocPtr slp,
  SeqLitPtr lit,
  CharPtr str,
  FILE *fp,
  ByteStorePtr bs,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen,
  Boolean do_defline,
  Boolean substitute_ids,
  Boolean sorted_prot,
  Int2 skip
)

{
  Char         acc [41];
  SeqIdPtr     accn = NULL;
  Char         buf [4096];
  Char         ch, ch1, ch2, ch3;
  Int4         count = 0;
  Int4         gi = -1;
  SeqIdPtr     gpp = NULL;
  Char         id [128];
  Uint1        id_format = PRINTID_FASTA_LONG;
  CharPtr      ptr;
  StreamFsa    sf;
  SeqIdPtr     sip = NULL;
  Char         spn [64];
  CharPtr      tmp;

  if (bsp == NULL && slp == NULL && lit == NULL && str == NULL) return 0;
  if (fp == NULL && bs == NULL) return 0;
  if (bsp != NULL && bsp->repr == Seq_repr_virtual) return 0;
  if (linelen > 128) {
    linelen = 128;
  }
  if (linelen < 1) {
    linelen = 60;
  }
  if (blocklen > 100) {
    blocklen = 100;
  }
  if (blocklen < 1) {
    blocklen = 0;
  }
  if (grouplen > 100) {
    grouplen = 100;
  }
  if (grouplen < 1) {
    grouplen = 0;
  }
  acc [0] = '\0';
  MemSet ((Pointer) &sf, 0, sizeof (StreamFsa));
  sf.fp = fp;
  sf.bs = bs;
  sf.idx = 0;
  sf.lin = 0;
  sf.blk = 0;
  sf.grp = 0;
  sf.linelen = linelen;
  sf.blocklen = blocklen;
  sf.grouplen = grouplen;
  sf.skip = skip;
  sf.gi = 0;
  sf.start = 0;
  sf.seqpos = 0;
  sf.seqspans = (Boolean) ((flags & STREAM_HTML_SPANS) != 0);
  if (sf.seqspans) {
    if (bsp != NULL) {
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        switch (sip->choice) {
          case SEQID_GI :
            gi = sip->data.intvalue;
            break;
          case SEQID_GENBANK :
          case SEQID_EMBL :
          case SEQID_DDBJ :
          case SEQID_OTHER :
            accn = sip;
            break;
          case SEQID_PIR :
          case SEQID_SWISSPROT :
          case SEQID_PRF :
          case SEQID_PDB :
            accn = sip;
            break;
          case SEQID_TPG :
          case SEQID_TPE :
          case SEQID_TPD :
            accn = sip;
            break;
          case SEQID_GPIPE :
            /* should not override better accession */
            gpp = sip;
            break;
          default :
            break;
        }
      }
    } else if (slp != NULL) {
      /* PUBSEQ_OS will send a SeqInt with a chain of Seq-ids */
      for (sip = SeqLocId (slp); sip != NULL; sip = sip->next) {
        switch (sip->choice) {
          case SEQID_GI :
            gi = sip->data.intvalue;
            break;
          case SEQID_GENBANK :
          case SEQID_EMBL :
          case SEQID_DDBJ :
          case SEQID_OTHER :
            accn = sip;
            break;
          case SEQID_PIR :
          case SEQID_SWISSPROT :
          case SEQID_PRF :
          case SEQID_PDB :
            accn = sip;
            break;
          case SEQID_TPG :
          case SEQID_TPE :
          case SEQID_TPD :
            accn = sip;
            break;
          case SEQID_GPIPE :
            /* should not override better accession */
            gpp = sip;
            break;
          default :
            break;
        }
      }
      if (sip != NULL && sip->choice == SEQID_GI) {
        sf.gi = sip->data.intvalue;
      }
    }
    if (gi > 0) {
      sf.gi = gi;
    }
    if (accn == NULL) {
      accn = gpp;
    }
    if (accn != NULL) {
      SeqIdWrite (accn, acc, PRINTID_TEXTID_ACC_ONLY, sizeof (acc) - 1);

      if (accn->choice == SEQID_PDB) {
        ptr = StringChr (acc, '_');
        if (ptr != NULL) {
          ch1 = ptr [1];
          if (ch1 != '\0') {
            ch2 = ptr [2];
            if (ch2 != '\0') {
              ch3 = ptr [3];
              if (ch3 == '\0') {
                if (ch1 == ch2) {
                  if (IS_UPPER (ch1)) {
                    ptr [1] = TO_LOWER (ch1);
                    ptr [2] = '\0';
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (do_defline) {
    id [0] = '\0';
    if (substitute_ids) {
      sip = ChooseFastaID (bsp, sorted_prot);
    } else if (bsp != NULL) {
      sip = bsp->id;
    }
    if ((flags & STREAM_ALL_FASTA_IDS) != 0) {
      id_format = PRINTID_FASTA_ALL;
    }
    SeqIdWrite (sip, id, id_format, sizeof (id) - 1);
    /* no longer need to do feature indexing if title not present to speed up creation */
    /*
    sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_title, NULL);
    if (sdp == NULL) {
      entityID = ObjMgrGetEntityIDForPointer (bsp);
      if (! SeqMgrFeaturesAreIndexed (entityID)) {
        SeqMgrIndexFeatures (entityID, NULL);
      }
    }
    */
    buf [0] = '\0';
    NewCreateDefLineBuf (NULL, bsp, buf, sizeof (buf), FALSE, FALSE);
    tmp = buf;
    ch = *tmp;
    while (ch != '\0') {
      if (ch == '>') {
        *tmp = '_';
      }
      tmp++;
      ch = *tmp;
    }
    if (sf.fp != NULL) {
      fprintf (fp, ">%s %s\n", id, buf);
    } else if (sf.bs != NULL) {
      BSWrite (sf.bs, ">", sizeof (">"));
      BSWrite (sf.bs, id, StringLen (id));
      BSWrite (sf.bs, " ", sizeof (" "));
      BSWrite (sf.bs, buf, StringLen (buf));
      BSWrite (sf.bs, "\n", sizeof ("\n"));
    }
  }
  if (bsp != NULL) {
    count = SeqPortStream (bsp, flags, (Pointer) &sf, FsaStreamProc);
  } else if (slp != NULL) {
    count = SeqPortStreamLoc (slp, flags, (Pointer) &sf, FsaStreamProc);
  } else if (lit != NULL) {
    count = SeqPortStreamLit (lit, flags, (Pointer) &sf, FsaStreamProc);
  } else if (str != NULL) {
    count = StringLen (str);
    FsaStreamProc (str, (Pointer) &sf);
  }
  /* print any remaining sequence */
  if (sf.lin > 0) {
    sf.buf [sf.idx] = '\0';
    if (sf.grp >= sf.grouplen && sf.grouplen > 0) {
      if (sf.fp != NULL) {
        fprintf (fp, "\n");
      } else if (sf.bs != NULL) {
        BSWrite (sf.bs, "\n", sizeof ("\n"));
      }
    }
    if (sf.fp != NULL) {
      if (sf.seqspans) {
        fprintf (sf.fp, "<span class=\"ff_line\" id=\"gi_%ld_%ld\">", (long) sf.gi, (long) (sf.start + 1));
      }
      fprintf (sf.fp, "%s", sf.buf);
      if (sf.seqspans) {
        fprintf (sf.fp, "</span>");
      }
      fprintf (sf.fp, "\n");
      if (sf.seqspans) {
        fprintf (sf.fp, "<script type=\"text/javascript\">");
        fprintf (sf.fp, "if (typeof(oData) == \"undefined\") oData = []; ");
        fprintf (sf.fp, "oData.push({gi:%ld,acc:\"%s\"})", (long) sf.gi, acc);
        fprintf (sf.fp, "</script>\n");
      }
    } else if (sf.bs != NULL) {
      if (sf.seqspans) {
        sprintf (spn, "<span class=\"ff_line\" id=\"gi_%ld_%ld\">", (long) sf.gi, (long) (sf.start + 1));
        BSWrite (sf.bs, spn, StringLen (spn));
      }
      BSWrite (sf.bs, sf.buf, StringLen (sf.buf));
      if (sf.seqspans) {
        BSWrite (sf.bs, "</span>", sizeof ("</span>"));
      }
      BSWrite (sf.bs, "\n", sizeof ("\n"));
      if (sf.seqspans) {
        sprintf (spn, "<script type=\"text/javascript\">");
        BSWrite (sf.bs, spn, StringLen (spn));
        sprintf (spn, "if (typeof(oData) == \"undefined\") oData = []; ");
        BSWrite (sf.bs, spn, StringLen (spn));
        sprintf (spn, "oData.push({gi:%ld,acc:\"%s\"})", (long) sf.gi, acc);
        BSWrite (sf.bs, spn, StringLen (spn));
        sprintf (spn, "</script>\n");
        BSWrite (sf.bs, spn, StringLen (spn));
      }
    }
  }
  return count;
}

NLM_EXTERN Int4 BioseqFastaStream (
  BioseqPtr bsp,
  FILE *fp,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen,
  Boolean do_defline
)

{
  return BioseqFastaStreamInternal (bsp, NULL, NULL, NULL, fp, NULL, flags,
                                    linelen, blocklen, grouplen,
                                    do_defline, FALSE, FALSE, 0);
}

NLM_EXTERN Int4 BioseqFastaStreamEx (
  BioseqPtr bsp,
  FILE *fp,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen,
  Boolean do_defline,
  Boolean substitute_ids,
  Boolean sorted_protein
)

{
  return BioseqFastaStreamInternal (bsp, NULL, NULL, NULL, fp, NULL, flags,
                                    linelen, blocklen, grouplen,
                                    do_defline, substitute_ids, sorted_protein, 0);
}

NLM_EXTERN Int4 BioseqFastaMemStream (
  BioseqPtr bsp,
  ByteStorePtr bs,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen,
  Boolean do_defline
)

{
  return BioseqFastaStreamInternal (bsp, NULL, NULL, NULL, NULL, bs, flags,
                                    linelen, blocklen, grouplen,
                                    do_defline, FALSE, FALSE, 0);
}

NLM_EXTERN Int4 SeqLocFastaStream (
  SeqLocPtr slp,
  FILE *fp,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen
)

{
  if (slp == NULL || fp == NULL) return 0;

  return BioseqFastaStreamInternal (NULL, slp, NULL, NULL, fp, NULL, flags, 
                                    linelen, blocklen, grouplen,
                                    FALSE, FALSE, FALSE, 0);
}

NLM_EXTERN Int4 SeqLitFastaStream (
  SeqLitPtr lit,
  FILE *fp,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen
)

{
  if (lit == NULL || fp == NULL) return 0;

  return BioseqFastaStreamInternal (NULL, NULL, lit, NULL, fp, NULL, flags, 
                                    linelen, blocklen, grouplen,
                                    FALSE, FALSE, FALSE, 0);
}

static void DoSpecialDefline (
  SeqFeatPtr sfp,
  FILE *fp,
  CdRegionPtr crp,
  CharPtr idSuffix
)

{
  BioseqPtr          bsp = NULL;
  Char               buf [512];
  SeqMgrFeatContext  cdscontext;
  Boolean            do_defline = TRUE;
  Uint2              entityID;
  SeqFeatPtr         gene = NULL;
  SeqMgrFeatContext  genecontext;
  Int4               gi;
  GeneRefPtr         grp;
  IntAsn2gbJob       iaj;
  Boolean            partial5;
  Boolean            partial3;
  BioseqPtr          prod;
  SeqIdPtr           sip;
  CharPtr            str;
  Char               tmp [64];

  if (sfp == NULL || fp == NULL || crp == NULL) return;

  MemSet ((Pointer) &genecontext, 0, sizeof (SeqMgrFeatContext));
  MemSet ((Pointer) &cdscontext, 0, sizeof (SeqMgrFeatContext));

  if (do_defline) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp == NULL) {
      do_defline = FALSE;
      StringCpy (buf, "lcl|");
      sip = SeqLocId (sfp->location);
      if (sip != NULL) {
        SeqIdWrite (sip, tmp, PRINTID_TEXTID_ACC_VER, sizeof (tmp) - 1);
        StringCat (buf, tmp);
      }
      if (StringDoesHaveText (idSuffix) && StringLen (idSuffix) < 200) {
        StringCat (buf, idSuffix);
      }
      FastaFileFunc (bsp, FASTA_ID, buf, sizeof (buf), (Pointer) fp);
      StringCpy (buf, "?");
      FastaFileFunc (bsp, FASTA_DEFLINE, buf, sizeof (buf), (Pointer) fp);
      fflush (fp);
    }
  }

  if (do_defline && bsp != NULL) {
    if (sfp != SeqMgrGetDesiredFeature (0, bsp, 0, 0, sfp, &cdscontext)) {
      do_defline = FALSE;
      StringCpy (buf, "lcl|");
      sip = SeqIdFindWorst (bsp->id);
      if (sip != NULL) {
        SeqIdWrite (sip, tmp, PRINTID_TEXTID_ACC_VER, sizeof (tmp) - 1);
        StringCat (buf, tmp);
      }
      if (StringDoesHaveText (idSuffix) && StringLen (idSuffix) < 200) {
        StringCat (buf, idSuffix);
      }
      FastaFileFunc (bsp, FASTA_ID, buf, sizeof (buf), (Pointer) fp);
      StringCpy (buf, "?");
      FastaFileFunc (bsp, FASTA_DEFLINE, buf, sizeof (buf), (Pointer) fp);
      fflush (fp);
    }
  }

  if (do_defline) {
    entityID = ObjMgrGetEntityIDForPointer (bsp);
    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, NULL);
    }

    CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  
    grp = SeqMgrGetGeneXref (sfp);
    if (grp == NULL || (! SeqMgrGeneIsSuppressed (grp))) {
      gene = SeqMgrGetOverlappingGene (sfp->location, &genecontext);
    }
  
    MemSet ((Pointer) &iaj, 0, sizeof (IntAsn2gbJob));
    iaj.flags.iupacaaOnly = FALSE;
    iaj.relModeError = FALSE;
  
    StringCpy (buf, "lcl|");
    sip = SeqIdFindWorst (bsp->id);
    if (sip != NULL) {
      SeqIdWrite (sip, tmp, PRINTID_TEXTID_ACC_VER, sizeof (tmp) - 1);
      StringCat (buf, tmp);
    }
    if (StringDoesHaveText (idSuffix) && StringLen (idSuffix) < 200) {
      StringCat (buf, idSuffix);
    }
  
    FastaFileFunc (bsp, FASTA_ID, buf, sizeof (buf), (Pointer) fp);
  
    buf [0] = '\0';
    if (StringDoesHaveText (genecontext.label)) {
      StringCat (buf, "[gene=");
      StringCat (buf, genecontext.label);
      StringCat (buf, "] ");
    }
    if (StringDoesHaveText (cdscontext.label)) {
      StringCat (buf, "[protein=");
      StringCat (buf, cdscontext.label);
      StringCat (buf, "] ");
    }
    if (crp->frame == 2) {
      StringCat (buf, "[frame=2] ");
    } else if (crp->frame == 3) {
      StringCat (buf, "[frame=3] ");
    }
    if (partial5 && partial3) {
      StringCat (buf, "[partial=5',3'] ");
    } else if (partial5) {
      StringCat (buf, "[partial=5'] ");
    } else if (partial3) {
      StringCat (buf, "[partial=3'] ");
    }
    if (sfp->product != NULL) {
      tmp [0] = '\0';
      sip = SeqLocId (sfp->product);
      if (sip != NULL && sip->choice == SEQID_GI) {
        prod = BioseqFind (sip);
        if (prod != NULL) {
          sip = SeqIdFindWorst (prod->id);
          SeqIdWrite (sip, tmp, PRINTID_TEXTID_ACC_VER, sizeof (tmp) - 1);
        } else {
          gi = sip->data.intvalue;
          sip = GetSeqIdForGI (gi);
          SeqIdWrite (sip, tmp, PRINTID_TEXTID_ACC_VER, sizeof (tmp));
          SeqIdFree (sip);
        }
      } else if (sip != NULL) {
        SeqIdWrite (sip, tmp, PRINTID_TEXTID_ACC_VER, sizeof (tmp));
      }
      if (StringDoesHaveText (tmp)) {
        StringCat (buf, "[protein_id=");
        StringCat (buf, tmp);
        StringCat (buf, "] ");
      }
    }
    str = FFFlatLoc (&iaj, bsp, sfp->location, FALSE, FALSE);
    if (str != NULL && StringLen (str) + StringLen (buf) < sizeof (buf) - 10) {
      StringCat (buf, "[location=");
      StringCat (buf, str);
      StringCat (buf, "] ");
      MemFree (str);
    }
    TrimSpacesAroundString (buf);
  
    FastaFileFunc (bsp, FASTA_DEFLINE, buf, sizeof (buf), (Pointer) fp);
  
    fflush (fp);
  }
}

NLM_EXTERN Int4 CdRegionFastaStream (
  SeqFeatPtr sfp,
  FILE *fp,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen,
  Boolean do_defline,
  CharPtr idSuffix
)

{
  CdRegionPtr  crp;
  Int2         skip = 0;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return 0;
  if (fp == NULL) return 0;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL) return 0;

  if (do_defline) {
    DoSpecialDefline (sfp, fp, crp, idSuffix);
  }

  if (crp->frame == 2) {
    skip = 1;
  } else if (crp->frame == 3) {
    skip = 2;
  }

  return BioseqFastaStreamInternal (NULL, sfp->location, NULL, NULL, fp, NULL, flags,
                                    linelen, blocklen, grouplen,
                                    FALSE, FALSE, FALSE, skip);
}

NLM_EXTERN Int4 TranslationFastaStream (
  SeqFeatPtr sfp,
  FILE *fp,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen,
  Boolean do_defline,
  CharPtr idSuffix
)

{
  ByteStorePtr  bs;
  Char          ch;
  Int4          count = 0;
  CdRegionPtr   crp;
  size_t        prtlen;
  CharPtr       ptr;
  CharPtr       str;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return 0;
  if (fp == NULL) return 0;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL) return 0;

  if (do_defline) {
    DoSpecialDefline (sfp, fp, crp, idSuffix);
  }

  str = NULL;
  bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
  str = BSMerge (bs, NULL);
  bs = BSFree (bs);

  if (str != NULL) {
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      *ptr = TO_UPPER (ch);
      ptr++;
      ch = *ptr;
    }
    prtlen = StringLen (str);
    if (prtlen > 1) {
       if (str [prtlen - 1] == '*') {
         str [prtlen - 1] = '\0';
       }
    }
  }

  count = BioseqFastaStreamInternal (NULL, NULL, NULL, str, fp, NULL, flags,
                                     linelen, blocklen, grouplen,
                                     FALSE, FALSE, FALSE, 0);

  MemFree (str);
 
  return count;
}

static void DoGeneDefline (
  SeqFeatPtr sfp,
  FILE *fp,
  GeneRefPtr grp,
  CharPtr idSuffix
)

{
  BioseqPtr          bsp = NULL;
  Char               buf [512];
  Boolean            do_defline = TRUE;
  Uint2              entityID;
  SeqMgrFeatContext  genecontext;
  IntAsn2gbJob       iaj;
  Boolean            partial5;
  Boolean            partial3;
  SeqIdPtr           sip;
  CharPtr            str;
  Char               tmp [64];

  if (sfp == NULL || fp == NULL || grp == NULL) return;
  if (sfp == NULL || fp == NULL || sfp->data.choice != SEQFEAT_GENE) return;
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return;

  if (do_defline) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp == NULL) {
      do_defline = FALSE;
      StringCpy (buf, "lcl|");
      sip = SeqLocId (sfp->location);
      if (sip != NULL) {
        SeqIdWrite (sip, tmp, PRINTID_TEXTID_ACC_VER, sizeof (tmp) - 1);
        StringCat (buf, tmp);
      }
      if (StringDoesHaveText (idSuffix) && StringLen (idSuffix) < 200) {
        StringCat (buf, idSuffix);
      }
      FastaFileFunc (bsp, FASTA_ID, buf, sizeof (buf), (Pointer) fp);
      StringCpy (buf, "?");
      FastaFileFunc (bsp, FASTA_DEFLINE, buf, sizeof (buf), (Pointer) fp);
      fflush (fp);
    }
  }

  if (do_defline && bsp != NULL) {
    if (sfp != SeqMgrGetDesiredFeature (0, bsp, 0, 0, sfp, &genecontext)) {
      do_defline = FALSE;
      StringCpy (buf, "lcl|");
      sip = SeqIdFindWorst (bsp->id);
      if (sip != NULL) {
        SeqIdWrite (sip, tmp, PRINTID_TEXTID_ACC_VER, sizeof (tmp) - 1);
        StringCat (buf, tmp);
      }
      if (StringDoesHaveText (idSuffix) && StringLen (idSuffix) < 200) {
        StringCat (buf, idSuffix);
      }
      FastaFileFunc (bsp, FASTA_ID, buf, sizeof (buf), (Pointer) fp);
      StringCpy (buf, "?");
      FastaFileFunc (bsp, FASTA_DEFLINE, buf, sizeof (buf), (Pointer) fp);
      fflush (fp);
    }
  }

  if (do_defline) {
    entityID = ObjMgrGetEntityIDForPointer (bsp);
    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, NULL);
    }

    CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  
    MemSet ((Pointer) &iaj, 0, sizeof (IntAsn2gbJob));
    iaj.flags.iupacaaOnly = FALSE;
    iaj.relModeError = FALSE;
  
    StringCpy (buf, "lcl|");
    sip = SeqIdFindWorst (bsp->id);
    if (sip != NULL) {
      SeqIdWrite (sip, tmp, PRINTID_TEXTID_ACC_VER, sizeof (tmp) - 1);
      StringCat (buf, tmp);
    }
    if (StringDoesHaveText (idSuffix) && StringLen (idSuffix) < 200) {
      StringCat (buf, idSuffix);
    }
  
    FastaFileFunc (bsp, FASTA_ID, buf, sizeof (buf), (Pointer) fp);
  
    buf [0] = '\0';
    if (StringDoesHaveText (grp->locus)) {
      StringCat (buf, "[gene=");
      StringCat (buf, grp->locus);
      StringCat (buf, "] ");
    }
    if (StringDoesHaveText (grp->locus_tag)) {
      StringCat (buf, "[locus_tag=");
      StringCat (buf, grp->locus_tag);
      StringCat (buf, "] ");
    }
    if (StringLen (buf) == 0 && StringDoesHaveText (genecontext.label)) {
      StringCat (buf, "[gene=");
      StringCat (buf, genecontext.label);
      StringCat (buf, "] ");
    }
    str = FFFlatLoc (&iaj, bsp, sfp->location, FALSE, FALSE);
    if (str != NULL && StringLen (str) + StringLen (buf) < sizeof (buf) - 10) {
      StringCat (buf, "[location=");
      StringCat (buf, str);
      StringCat (buf, "] ");
      MemFree (str);
    }
    TrimSpacesAroundString (buf);
  
    FastaFileFunc (bsp, FASTA_DEFLINE, buf, sizeof (buf), (Pointer) fp);
  
    fflush (fp);
  }
}

NLM_EXTERN Int4 GeneFastaStream (
  SeqFeatPtr sfp,
  FILE *fp,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen,
  Boolean do_defline,
  CharPtr idSuffix
)

{
  GeneRefPtr  grp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE) return 0;
  if (fp == NULL) return 0;
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return 0;

  if (do_defline) {
    DoGeneDefline (sfp, fp, grp, idSuffix);
  }

  return BioseqFastaStreamInternal (NULL, sfp->location, NULL, NULL, fp, NULL, flags,
                                    linelen, blocklen, grouplen,
                                    FALSE, FALSE, FALSE, 0);
}

/*****************************************************************************
*
*   SeqEntryFastaStream (bsp, fp, flags, linelen, blocklen, grouplen,
*                        do_na, do_aa, master_style)
*
*       Rapid FASTA generator on ASN.1 record including GenBank release set
*
*****************************************************************************/

typedef struct fastastreamdata {
  FILE           *fp;
  StreamFlgType  flags;
  Int2           linelen;
  Int2           blocklen;
  Int2           grouplen;
  Boolean        do_na;
  Boolean        do_aa;
  Boolean        master_style;
  Boolean        failed;
  Int4           count;
  Boolean        substitute_ids;
  Boolean        sorted_prot;
} FastaStreamData, PNTR FastaStreamPtr;

static BioseqSetPtr GetSegParts (
  BioseqPtr bsp
)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   sep;
  if (bsp == NULL || bsp->repr != Seq_repr_seg) return NULL;
  sep = bsp->seqentry;
  if (sep == NULL) return NULL;
  sep = sep->next;
  if (sep == NULL || (! IS_Bioseq_set (sep))) return NULL;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp != NULL && bssp->_class == BioseqseqSet_class_parts) return bssp;
  return NULL;
}

static void FastaOneBioseq (
  BioseqPtr bsp,
  Pointer userdata
)

{  Int4            count;
  FastaStreamPtr  fsp;
  BioseqSetPtr    parts;
  if (bsp == NULL) return;
  fsp = (FastaStreamPtr) userdata;
  if (fsp == NULL) return;
  /* return if molecule not right for format */
  if (ISA_na (bsp->mol)) {
    if (! fsp->do_na) return;
  } else if (ISA_aa (bsp->mol)) {
    if (! fsp->do_aa) return;
  }
  if (bsp->repr == Seq_repr_seg && (! fsp->master_style)) {
    /* if bsp followed by parts set, recurse to make FASTA from individual parts */
    parts = GetSegParts (bsp);
    if (parts != NULL) {
      VisitBioseqsInSet (parts, (Pointer) fsp, FastaOneBioseq);
      return;
    }
  }
  if (bsp->repr == Seq_repr_raw ||
      bsp->repr == Seq_repr_seg ||
      bsp->repr == Seq_repr_const ||
      bsp->repr == Seq_repr_delta ||
      bsp->repr == Seq_repr_virtual) {
    count = BioseqFastaStreamEx (bsp, fsp->fp, fsp->flags, fsp->linelen, fsp->blocklen, fsp->grouplen,
                                 TRUE, fsp->substitute_ids, fsp->sorted_prot);
    if (count < 0) {
      fsp->failed = TRUE;
      fsp->count -= count;
    } else {
      fsp->count += count;
    }
  }
}

NLM_EXTERN Int4 SeqEntryFastaStreamEx (
  SeqEntryPtr sep,
  FILE *fp,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen,
  Boolean do_na,
  Boolean do_aa,
  Boolean master_style,
  Boolean substitute_ids,
  Boolean sorted_prot
)

{  BioseqPtr        bsp = NULL;
  BioseqSetPtr     bssp = NULL;
  Uint2            entityID = 0;
  FastaStreamData  fsd;
  SeqEntryPtr      oldscope;
  if (sep == NULL || fp == NULL) return 0;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    entityID = ObjMgrGetEntityIDForPointer (bsp);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    entityID = ObjMgrGetEntityIDForPointer (bssp);
  }
  if (entityID == 0) return 0;
  /* AssignIDs sets bsp->seqentry so GetSegParts can work */
  AssignIDsInEntity (entityID, 0, NULL);
  fsd.fp = fp;
  fsd.flags = flags;
  fsd.linelen = linelen;
  fsd.blocklen = blocklen;
  fsd.grouplen = grouplen;
  fsd.do_na = do_na;
  fsd.do_aa = do_aa;
  fsd.master_style = master_style;
  fsd.failed = FALSE;
  fsd.count = 0;
  fsd.substitute_ids = substitute_ids;
  fsd.sorted_prot = sorted_prot;
  oldscope = SeqEntrySetScope (sep);
  if (bssp != NULL) {
    /* handle all components of a pop/phy/mut/eco set */
    sep = SeqMgrGetSeqEntryForData (bssp);
    VisitSequencesInSep (sep, (Pointer) &fsd, VISIT_MAINS, FastaOneBioseq);
  } else {
    /* handle single bioseq, which may be segmented or a local part */
    FastaOneBioseq (bsp, (Pointer) &fsd);
  }
  SeqEntrySetScope (oldscope);
  if (fsd.failed) {
    return -fsd.count;
  }
  return fsd.count;
}

NLM_EXTERN Int4 SeqEntryFastaStream (
  SeqEntryPtr sep,
  FILE *fp,
  StreamFlgType flags,
  Int2 linelen,
  Int2 blocklen,
  Int2 grouplen,
  Boolean do_na,
  Boolean do_aa,
  Boolean master_style
)

{  return SeqEntryFastaStreamEx (sep, fp, flags, linelen, blocklen, grouplen, do_na, do_aa, master_style, FALSE, FALSE);
}

/*****************************************************************************
*
*   Here are functions that convert FASTA format from file or from memory
*
*****************************************************************************/
/********* DEFINES *********/
#define FTSE_BUFF_CHUNK 4096
#define BIOSEQ 1
/********* INTERNAL FUNCTIONS *********/
NLM_EXTERN SeqEntryPtr FastaToSeqEntryInternalEx
(
 VoidPtr input,          /* input pointer (file or memory) */
 Int4 type,              /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR next_char, /* returned pointer to next FASTA sequence */
 Boolean is_na,          /* type of sequence */
 CharPtr PNTR errormsg,  /* error messge for debugging */
 Boolean parseSeqId,     /* Parse SeqID from def line */
 CharPtr special_symbol, /* Returns special symbol if no SeqEntry */
 CharPtr prefix,         /* prefix for localID if not parsable */
 Int2Ptr ctrptr,         /* starting point for constructing unique ID */
 SeqLocPtr PNTR mask_ptr /* Pointer to a SeqLoc to Fill with Masking information */
 );
NLM_EXTERN SeqEntryPtr FastaToSeqEntryInternal
(
 VoidPtr input,          /* input pointer (file or memory) */
 Int4 type,              /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR last_char, /* returned pointer to next FASTA sequence */
 Boolean is_na,          /* type of sequence */
 CharPtr PNTR errormsg,  /* error messge for debugging */
 Boolean parseSeqId,     /* Parse SeqID from def line */
 CharPtr special_symbol  /* Returns special symbol if no SeqEntry */
 );
static Boolean FastaReadSequenceInternal
(
 VoidPtr input,          /* input pointer (file or memory) */
 Int4 type,              /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR last_char, /* returned pointer to next FASTA sequence */
 Boolean is_na,          /* type of sequence */
 Int4Ptr seq_length,     /* Returned length of sequence in residues */
 ByteStorePtr PNTR,      /* Returned pointer to sequence ByteStore */
 CharPtr PNTR errormsg,  /* error messge for debugging */
 CharPtr special_symbol  /* Returns special symbol if no SeqEntry */
 );
static Boolean FastaReadSequenceInternalEx
(
 VoidPtr input,          /* input pointer (file or memory) */
 Int4 type,              /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR last_char, /* returned pointer to next FASTA sequence */
 Boolean is_na,          /* type of sequence */
 Int4Ptr seq_length,     /* Returned length of sequence in residues */
 ByteStorePtr PNTR,      /* Returned pointer to sequence ByteStore */
 CharPtr PNTR errormsg,  /* error messge for debugging */
 CharPtr special_symbol, /* Returns special symbol if no SeqEntry */
 SeqLocPtr PNTR mask_ptr,/* Pointer to a SeqLoc to Fill with Masking information */
 SeqIdPtr sip            /* SeqId of current sequence used for Masking Info */
 );
static Int4 FastaReadSequenceChunk
(
 VoidPtr input,          /* input pointer (file or memory) */
 Int4    type,           /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR next_char, /* returned pointer to next FASTA sequence */
 Uint1Ptr sequence,      /* buffer to read sequence to */
 Int4     length,        /* size of buffer */
 CharPtr special_symbol  /* Returns special symbol if no SeqEntry */
 );
static SeqEntryPtr FastaToSeqEntryInternalExEx
(
 VoidPtr input,           /* input pointer (file or memory) */
 Int4 type,               /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR next_char,  /* returned pointer to next FASTA sequence */
 Boolean is_na,           /* type of sequence */
 CharPtr PNTR errormsg,   /* error messge for debugging */
 Boolean parseSeqId,      /* Parse SeqID from def line */
 CharPtr special_symbol,  /* Returns special symbol if no SeqEntry */
 CharPtr prefix,          /* prefix for localID if not parsable */
 Int2Ptr ctrptr,          /* starting point for constructing unique ID */
 SeqLocPtr PNTR mask_ptr, /* Pointer to a SeqLoc to Fill with Masking information */
 Boolean trustID
 );
/********* FINCTIONS *********/
/*****************************************************************************
*
*   SeqEntryPFtr FastaToSeqBuffEx() - function to return SeqEntryPtr from
*                                    buffer with error handling
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr FastaToSeqBuffEx
  (
    CharPtr buffer,         /* buffer in memory with FASTA sequence */
    CharPtr PNTR last_char, /* here returned pointer to next FASTA if any */
    Boolean is_na,          /* type of sequence */
    CharPtr PNTR errormsg,  /* error message for debugging */
    Boolean parseSeqId      /* Parse SeqID from def line */
  )
{
  return FastaToSeqEntryInternal((void *)buffer, FASTA_MEM_IO ,
                                 last_char, is_na, errormsg, parseSeqId, NULL);
}
NLM_EXTERN SeqEntryPtr FastaToSeqBuffForDb
  (
    CharPtr buffer,         /* buffer in memory with FASTA sequence */
    CharPtr PNTR last_char, /* here returned pointer to next FASTA if any */
    Boolean is_na,          /* type of sequence */
    CharPtr PNTR errormsg,  /* error message for debugging */
    Boolean parseSeqId,     /* Parse SeqID from def line */
    CharPtr prefix,         /* prefix for localID if not parsable */
    Int2Ptr ctrptr,         /* starting point for constructing unique ID */
    SeqLocPtr PNTR mask_ptr /* Pointer to a SeqLoc to Fill with Masking information from lowercased letters */
  )
{
  return FastaToSeqEntryInternalExEx((void *)buffer, FASTA_MEM_IO ,
                     last_char, is_na, errormsg, parseSeqId,
                     NULL, prefix, ctrptr, mask_ptr, TRUE);
}
/*****************************************************************************
*
*   SeqEntryPtr FastaToSeqEntryEx() - function to return SeqEntryPtr from
*                                     file with error handling
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr FastaToSeqEntryEx
  (
    FILE *fp,               /* file to get sequence from */
    Boolean is_na,          /* type of sequence */
    CharPtr PNTR errormsg,  /* error message for debugginq */
    Boolean parseSeqId      /* Parse SeqID from def line */
  )
{
  return FastaToSeqEntryInternal((void *)fp, FASTA_FILE_IO,
                                 NULL,is_na, errormsg, parseSeqId, NULL);
}
/*****************************************************************************
*
*   SeqEntryPtr FastaToSeqEntryForDb() - function to return SeqEntryPtr from
*                                     file with error handling and with control
*                                     over generation of unique SeqIDs
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr FastaToSeqEntryForDb
  (
    FILE *fp,               /* file to get sequence from */
    Boolean is_na,          /* type of sequence */
    CharPtr PNTR errormsg,  /* error message for debugginq */
    Boolean parseSeqId,     /* Parse SeqID from def line */
    CharPtr prefix,         /* prefix for localID if not parsable */
    Int2Ptr ctrptr,         /* starting point for constructing unique ID */
    SeqLocPtr PNTR mask_ptr /* Pointer to a SeqLoc to Fill with Masking information from lowercased letters */
  )
{
  return FastaToSeqEntryInternalExEx ((void *) fp, FASTA_FILE_IO,
                                 NULL, is_na, errormsg, parseSeqId,
                                 NULL, prefix, ctrptr, mask_ptr, TRUE);
}
/*****************************************************************************
*
*   SeqEntryPtr FastaToSeqEntry() - function to return SeqEntryPtr from
*                                   file without error handling
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr FastaToSeqEntry (FILE *fp, Boolean is_na)
{
  return FastaToSeqEntryEx (fp, is_na, NULL, TRUE);
}
/*****************************************************************************
*
*   SeqEntryPtr FastaToSeqBuff() - function to return SeqEntryPtr from
*                                   buffer without error handling
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr FastaToSeqBuff (CharPtr buffer, CharPtr PNTR last_char,
                                   Boolean is_na)
{
  return FastaToSeqBuffEx (buffer, last_char, is_na, NULL, TRUE);
}
/*****************************************************************************
*
*   Boolean FastaReadSequenceChunk() - read sequence chunkfrom
*                                      file or buffer for use in
*                                      FastaReadsequenceInternal()
*****************************************************************************/
static Int4 FastaReadSequenceChunk
(
 VoidPtr input,          /* input pointer (file or memory) */
 Int4    type,           /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR next_char, /* returned pointer to next FASTA sequence */
 Uint1Ptr sequence,      /* buffer to read sequence to */
 Int4     length,        /* size of buffer */
 CharPtr special_symbol  /* Returns special symbol if no SeqEntry */
 )
{
    const Char PNTR firstchar;
    FILE *fd;
    register Int4 i;
    Int2 ch = 0;
    /* Type of input depends upon calling function */
    if(type == FASTA_FILE_IO) {
    fd = (FILE *) input;
    /* Skip empty lines and lines starting with a comment symbol. */
    while (1) {
        ch = NLM_GETC(fd);
        /* Ignore lines starting with a comment symbol. */
        if (ch == '!' || ch == '#') {
            do {
                ch = NLM_GETC(fd);
            } while (ch != '\n' && ch != '\r' && ch != '\0' && ch != EOF);
        }
        /* If end of file reached, return 0. */
        if (ch == EOF)
            return 0;
        /* If line not empty, break out of this loop. */
        if (ch != '\n' && ch != '\r')
            break;
    }
    if(ch == '>' || ch == '&' || ch == '{' || ch == '}' || ch == '[' || ch == ']')
    {
        ungetc(ch, fd);
        if (special_symbol != NULL) {
            *special_symbol = (Char) ch;
        }
        return 0;
    }
    sequence[0] = (Uint1) ch;
    if ((fgets((CharPtr) sequence+1, length-1, fd)) == NULL) {
        sequence [1] = '\0';
    }
    } else {   /* type == FASTA_MEM_IO */
        if((firstchar = (const Char PNTR) input) == NULL)
            return 0;
    }
    if(type == FASTA_FILE_IO) {
        for(i=0; i < length; i++) {
        if (sequence[i] == '\n' || sequence[i] == '\r' || sequence[i] == '\0')
        break;
        }
    } else { /* type = FASTA_MEM_IO */
        for(i =0; i < length && (ch = *firstchar) != NULLB; firstchar++, i++) {
            if((sequence[i] = (Char) ch) == '>' || (Char) ch == '&' || (Char) ch == '{' ||
                (Char) ch == '}' || (Char) ch == '[' || (Char) ch == ']') {
                if((i == 0) ||
                   (i > 0 && (sequence[i-1] == '\n' ||
                              sequence[i-1] == '\r'))) {
                    if (special_symbol != NULL) {
                        *special_symbol = (Char) ch;
                    }
                    break;
                }
            }
        }
        if(ch == NULLB) /* the end of buffer */
            *next_char = NULL;
        else
            *next_char = (CharPtr) firstchar;
    }
    return i;
}
/*****************************************************************************
*
*   Boolean FastaReadSequence() - read sequence from file
*
*****************************************************************************/
Boolean FastaReadSequence
(
 FILE *fd,            /* input pointer (file or memory) */
 Boolean is_na,            /* type of sequence */
 Int4Ptr seq_length,       /* Returned length of sequence in residues */
 ByteStorePtr PNTR bs_out, /* Returned pointer to sequence ByteStore */
 CharPtr PNTR errormsg     /* error message for debugging */
 )
{
    return  FastaReadSequenceInternal((VoidPtr) fd, FASTA_FILE_IO, NULL,
                                      is_na, seq_length, bs_out, errormsg, NULL);
}
/*****************************************************************************
*
*   Boolean FastaReadSequenceMem() - read sequence from buffer
*
*****************************************************************************/
Boolean FastaReadSequenceMem
(
 CharPtr buffer,           /* input buffer with sequence */
 CharPtr PNTR next_char,   /* returned pointer to next FASTA sequence */
 Boolean is_na,            /* type of sequence */
 Int4Ptr seq_length,       /* Returned length of sequence in residues */
 ByteStorePtr PNTR bs_out, /* Returned pointer to sequence ByteStore */
 CharPtr PNTR errormsg     /* error message for debugging */
)
{
    return  FastaReadSequenceInternal((VoidPtr) buffer, FASTA_MEM_IO,
                                      next_char, is_na, seq_length, bs_out,
                                      errormsg, NULL);
}
/*****************************************************************************
*
*   Boolean FastaReadSequenceInternal() - read sequence from
*                                         file or buffer for internal use
*
*****************************************************************************/
static Boolean FastaReadSequenceInternal
(
 VoidPtr input,            /* input pointer (file or memory) */
 Int4 type,                /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR next_char,   /* returned pointer to next FASTA sequence */
 Boolean is_na,            /* type of sequence */
 Int4Ptr seq_length,       /* Returned length of sequence in residues */
 ByteStorePtr PNTR bs_out, /* Returned pointer to sequence ByteStore */
 CharPtr PNTR errormsg,    /* error message for debugging */
 CharPtr special_symbol  /* Returns special symbol if no SeqEntry */
)
{
    return FastaReadSequenceInternalEx(input,type,next_char,is_na,seq_length,bs_out,errormsg,special_symbol,NULL, NULL);
}
/*****************************************************************************
*
*   Boolean FastaReadSequenceInternalEx() - read sequence from
*                                         file or buffer for internal use
*                                         and Create Masked SeqLoc of Lowercase sequences.
*
*****************************************************************************/
static Boolean FastaReadSequenceInternalEx
(
 VoidPtr input,            /* input pointer (file or memory) */
 Int4 type,                /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR next_char,   /* returned pointer to next FASTA sequence */
 Boolean is_na,            /* type of sequence */
 Int4Ptr seq_length,       /* Returned length of sequence in residues */
 ByteStorePtr PNTR bs_out, /* Returned pointer to sequence ByteStore */
 CharPtr PNTR errormsg,    /* error message for debugging */
 CharPtr special_symbol,   /* Returns special symbol if no SeqEntry */
 SeqLocPtr PNTR mask_ptr,  /* Pointer to a SeqLoc to Fill with Masking information */
 SeqIdPtr sip            /* SeqId of current sequence used for Masking Info */
)
{
    SeqMapTablePtr smtp;
    Uint1Ptr       in_buff, out_buff;
    CharPtr        ptr, chptr;
    Int2           ch;
    Uint1          byte_from, uch;
    register Int4  i;
    CharPtr        badchar = NULL;
    Int4           in_index, out_index, total_read, badchars = 0;
    Int4           total_length = 0;
    Int4           mask_to;
    Char           tmp[32];
    ValNodePtr     mask_head,mask,mask_new;
    SeqIntPtr      mask_sint;
    Boolean        Second, skip_to_eol, last_was_star;
    Boolean        this_char_masked;
    if (input == NULL)     /* empty input */
        return FALSE;
    /* Initializing conversion tables */
    if(is_na) {
        if((smtp = SeqMapTableFind(Seq_code_ncbi4na,
                                   Seq_code_iupacna)) == NULL) {
            return FALSE;
        }
    } else {
        if((smtp = SeqMapTableFind(Seq_code_ncbistdaa,
                                   Seq_code_ncbieaa)) == NULL) {
            return FALSE;
        }
    }
    /* Allocationg error message buffers if required */
    if (errormsg != NULL) {
        *errormsg = NULL;
        if((badchar = (CharPtr) MemNew(256)) == NULL)
            return FALSE;
    }
    if((in_buff = (Uint1Ptr) MemNew(FTSE_BUFF_CHUNK)) == NULL)
        return FALSE;
    if((out_buff = (Uint1Ptr) MemNew(FTSE_BUFF_CHUNK)) == NULL)
        return FALSE;
    if((*bs_out = BSNew(FTSE_BUFF_CHUNK)) == NULL)
        return FALSE;
    Second = FALSE;
    skip_to_eol = FALSE;
    last_was_star = FALSE;
    in_index = out_index = total_read = 0;
    if(mask_ptr) {
        mask_head=mask=NULL;
        mask_sint=NULL;
        this_char_masked=FALSE;
    }
    while(TRUE) {
        if (in_index == total_read) {
            if((total_read = FastaReadSequenceChunk(input, type,
                                                    next_char, in_buff,
                                                    FTSE_BUFF_CHUNK, special_symbol)) == 0)
                break; /* Here is exit from the loop */
            if(type == FASTA_MEM_IO)
                input = (VoidPtr) *next_char;
            in_index = 0;
        }
          byte_from = in_buff[in_index];
          in_index++;
        if ((! is_na) && (! last_was_star) && byte_from == '*') {
            last_was_star = TRUE;
        } else if(byte_from != ';' && !skip_to_eol) {
            if(mask_ptr) {
                if(IS_LOWER(byte_from)) {
                    if(this_char_masked) {
                        mask_to++;
                    } else { /* First lowercase character in this segment */
                        this_char_masked = TRUE;
                        /* save previous segment if any */
                        mask_new = ValNodeNew(NULL);
                        mask_new->choice = SEQLOC_INT;
                        if(mask_sint) {
                            mask_sint->to = mask_to;
                            mask->next = mask_new;
                        } else {
                            mask_head = mask_new;
                        }
                        mask = mask_new;
                        mask_sint = SeqIntNew();
                        mask_sint->from = total_length;
                        mask_sint->to = total_length;
                        mask_to = total_length;
                        mask_sint->strand = Seq_strand_both;
                        mask_sint->id = SeqIdDup(sip);
                        mask_new->data.ptrvalue = mask_sint;
                    }
                } else {
                    this_char_masked = FALSE;
                }
            }
            byte_from = TO_UPPER (byte_from);
            if (is_na && byte_from == 'U') byte_from = 'T';
            if (is_na && byte_from == 'X') byte_from = 'N';
            if((uch = SeqMapTableConvert(smtp, byte_from)) !=
               INVALID_RESIDUE && byte_from != '-') {
                if (last_was_star) {
                    total_length++;
                    out_buff[out_index] = SeqMapTableConvert(smtp, '*');
                    out_index++;
                    if(out_index == FTSE_BUFF_CHUNK) {
                        if(BSWrite(*bs_out, out_buff, out_index) != out_index) {
                            MemFree(badchar);
                            MemFree(in_buff);
                            MemFree(out_buff);
                            return FALSE;
                        }
                        out_index = 0;
                    }
                    last_was_star = FALSE;
                }
                total_length++;
                if(is_na) {
                    if(!Second) {
                        uch <<= 4;
                        out_buff[out_index] = uch;
                    } else {
                        out_buff[out_index] += uch;
                        out_index++;
                    }
                    Second = !Second;
                } else {
                    out_buff[out_index] = uch;
                    out_index++;
                }
            } else if (errormsg != NULL){
                if(IS_ALPHA(byte_from) || byte_from == '-' || byte_from == '?') {
                    (badchar [(int) (byte_from)])++;
                    badchars++;
                }
            }
        } else {    /* ch == ';' */
            /* We have to ignore rest of the line */
            skip_to_eol = TRUE;
            while(in_index < total_read  &&
                  (byte_from = in_buff[in_index]) != '\n' &&
                  byte_from != '\r')
                in_index++;
    /* Do not skip other passes if a line-return has
    been encountered as shown by examining less than the total
    (for FASTA_MEM_IO) or finding a line-return (FASTA_FILE_IO). */
            if(in_index < total_read ||
        (in_index < FTSE_BUFF_CHUNK &&
        (in_buff[in_index] == '\n' || in_buff[in_index] == '\r')))
                skip_to_eol = FALSE;
        }
        if(out_index == FTSE_BUFF_CHUNK) {
            if(BSWrite(*bs_out, out_buff, out_index) != out_index) {
                MemFree (badchar);
                MemFree(in_buff);
                MemFree(out_buff);
                return FALSE;
            }
            out_index = 0;
        }
    }  /* while (TRUE) */
    /* We have to write remaining stuff in out_buff */
    if(is_na && Second) out_index++; /* Partial byte for DNA */
    if(BSWrite(*bs_out, out_buff, out_index) != out_index) {
        MemFree (badchar);
        MemFree(in_buff);
        MemFree(out_buff);
        return FALSE;
    }
    *seq_length = total_length;
    /* If required bad characters statistics */
    if (errormsg != NULL && badchars > 0) {
        if((ptr = (CharPtr) MemNew (sizeof(Char)*512)) == NULL)
            return FALSE;
        chptr = "";
        sprintf (ptr, "%ld illegal %s %s removed:\n", (long) badchars,
                 badchars == 1 ? "character" : "characters",
                 badchars == 1 ? "was" : "were"
                 );
        for (ch = 'A', i =0; ch <= 'Z'; ch++, i++) {
            if ((badchar[ch]) > 0) {
                sprintf (tmp, "%s%d %c%s",
                         chptr, (int) badchar[ch], ch,
                         badchar[ch] == 1 ? "" : "s");
                StringCat (ptr, tmp);
                chptr = ", ";
            }
        }
        ch = '-';
        if ((badchar[ch]) > 0) {
            sprintf (tmp, "%s%d %c%s",
                     chptr, badchar[ch], ch,
                     badchar[ch] == 1 ? "" : "s");
            StringCat (ptr, tmp);
            chptr = ", ";
        }
        ch = '?';
        if ((badchar[ch]) > 0) {
            sprintf (tmp, "%s%d %c%s",
                     chptr, badchar[ch], ch,
                     badchar[ch] == 1 ? "" : "s");
            StringCat (ptr, tmp);
            chptr = ", ";
        }
        *errormsg = StringSave (ptr);
        MemFree (ptr);
    }
    MemFree (badchar);
    MemFree(in_buff);
    MemFree(out_buff);
    if(mask_ptr && mask_head) {
        SeqLocPtr slp;
        if(mask_sint) {
            mask_sint->to = mask_to;
            mask->next = NULL;
        }
        slp = SeqLocNew(NULL);
        slp->choice = SEQLOC_PACKED_INT;
        slp->data.ptrvalue = mask_head;
        *mask_ptr = slp;
    }
    return TRUE;
}
static SeqIdPtr MakeTrustedID (CharPtr prefix, Int2Ptr ctrptr)
{
  Char buf[40];
  ValNodePtr newid;
  ObjectIdPtr oid;
  Int2 start = 1;
    if (ctrptr != NULL) {
        start = *ctrptr;
    }
    if (start < 1) {
        start = 1;
    }
        if (prefix)
           sprintf (buf, "%d_%.32s", (int) start, prefix);
        else
           sprintf(buf, "%d", (int) start);
    newid = ValNodeNew (NULL);
    oid = ObjectIdNew ();
    newid->choice = SEQID_LOCAL;
    newid->data.ptrvalue = oid;
    oid->str = StringSave (buf);
    if (ctrptr != NULL) {
        *ctrptr = start + 1;
    }
    return newid;
}
static SeqEntryPtr FastaToSeqEntryInternalExEx
(
 VoidPtr input,           /* input pointer (file or memory) */
 Int4 type,               /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR next_char,  /* returned pointer to next FASTA sequence */
 Boolean is_na,           /* type of sequence */
 CharPtr PNTR errormsg,   /* error messge for debugging */
 Boolean parseSeqId,      /* Parse SeqID from def line */
 CharPtr special_symbol,  /* Returns special symbol if no SeqEntry */
 CharPtr prefix,          /* prefix for localID if not parsable */
 Int2Ptr ctrptr,          /* starting point for constructing unique ID */
 SeqLocPtr PNTR mask_ptr, /* Pointer to a SeqLoc to Fill with Masking information */
 Boolean trustID
 )
{
    SeqEntryPtr    sep = NULL;
    BioseqPtr      bsp = NULL;
    ValNodePtr     vnp = NULL;
    Int2           ch;
    CharPtr        chptr = NULL, ptr = NULL;
    register Int4  i;
    CharPtr        defline, buffer= NULL;   /* Working buffers */
    Int4           BuffSize = FTSE_BUFF_CHUNK;
    long           len = 0;
    FILE           *fd;
    const Char     PNTR firstchar;
    Boolean        is_gap = FALSE;
    if (special_symbol != NULL) {
        *special_symbol = '\0';
    }
    if (input == NULL)     /* empty input */
        return NULL;
    /* Type of input depends upon calling function */
    if(type == FASTA_FILE_IO)
        fd = (FILE *) input;
    else    /* type == FASTA_MEM_IO */
        firstchar = (const Char PNTR) input;
    /* Rolling spaces to check first non-space character */
    if(type == FASTA_FILE_IO) {
        do {
            ch = NLM_GETC(fd);
            if (ch == '!' || ch == '#') { /* comment symbol - ignore rest of line */
                do {
                    ch = NLM_GETC(fd);
                } while (ch != '\n' && ch != '\r' && ch != '\0' && ch != EOF);
            }
        } while (IS_WHITESP(ch));
    } else {   /* if(type == FASTA_MEM_IO*/
        while (IS_WHITESP(ch = *firstchar)) /* Rolling spaces */
            firstchar++;
    }
    if(ch == EOF || ch == NULLB || ch == '&' || ch == '{' ||
       ch == '}' || ch == '[' || ch == ']') {
        /* This is empty FILE or buffer or special symbol detected */
        if (special_symbol != NULL) {
            *special_symbol = ch;
        }
        return NULL;
    }
    /* First character is valid: initializing main structures */
    /* Initializing Seq-entry structure */
    if((sep = SeqEntryNew()) == NULL) {
        MemFree(buffer);
        return NULL;
    }
    sep->choice = BIOSEQ;  /* == 1 */
    if((bsp = BioseqNew()) == NULL) {
        MemFree(buffer);
        return NULL;
    }
    sep->data.ptrvalue = bsp;
    SeqMgrSeqEntry (SM_BIOSEQ, (Pointer)bsp, sep);
    if (is_na) {
        bsp->mol = Seq_mol_na;
        bsp->seq_data_type = Seq_code_ncbi4na;
    } else {
        bsp->mol = Seq_mol_aa;
        bsp->seq_data_type = Seq_code_ncbistdaa;
    }
    bsp->repr = Seq_repr_raw;
    /*  ------------- */
    /* Now reading defline into memory */
  /* DEFLINE PROCCESSING*/
    if(ch == '>') {     /* Defline is present - processing */
        if((buffer = (CharPtr) MemNew(BuffSize+1)) == NULL)
            return NULL;
        if(type == FASTA_FILE_IO) {    /* File */
            buffer[0] = (Char) ch;
        i = 0;
        fgets(buffer+1, BuffSize, fd);
        while (1)
        {
        while (i<BuffSize-1)
        {
                    if(buffer[i] == '\n' || buffer[i] == '\r' || buffer[i] == NULLB)
            {
                buffer[i] = NULLB;
                            break;
            }
            i++;
        }
        if (i == BuffSize-1 && (buffer[i] == '\n' || buffer[i] == '\r'))
        {
            buffer[i] = NULLB;
        }
        if (buffer[i] == NULLB)
            break;
                BuffSize = i + FTSE_BUFF_CHUNK;
                if((buffer = (CharPtr)Realloc(buffer, BuffSize+1)) == NULL)
        {
                          ErrLogPrintf("Error re-allocating memory in FastaToSeqEntry");
                        MemFree(buffer);
                        return NULL;
        }
            fgets(buffer+i+1, FTSE_BUFF_CHUNK, fd);
        }
        } else {  /* type = FASTA_MEM_IO */
            for(i =0; (ch = *firstchar) != NULLB; firstchar++, i++) {
                if (i >= BuffSize) {
                    BuffSize = i + FTSE_BUFF_CHUNK;
                    buffer = (CharPtr) Realloc(buffer, BuffSize);
                }
                if((buffer[i] = (Char) ch) == '\n' || ch == '\r') {
                    break;
                }
            }
            buffer[i] = NULLB;
            if(ch == NULLB) {/* the end of buffer */
                *next_char = NULL;
                input =  (VoidPtr) "\0";
            } else {
                *next_char = (CharPtr) firstchar;
                input = (VoidPtr) firstchar;
            }
        }
        defline = buffer+1;   /* Character after '>' */
        if(defline[0] != '?') {
            /* Creating standard Seq-id */
            ptr = defline;
            while (IS_WHITESP(*ptr))
                ptr++;
            if (parseSeqId) {
                if (*ptr == '"') {
                    ptr++;
                    chptr = StringChr (ptr, '"');
                } else {
                   for (chptr = ptr; *chptr != NULLB && !IS_WHITESP(*chptr);
                        chptr++) continue;
                   if (*chptr == NULLB)
                      chptr = NULL;
                }
            }
            if (!parseSeqId) {
                chptr = ptr;
            } else if (chptr != NULL) {
                *chptr = NULLB;
                chptr++;
                bsp->id = MakeSeqID (ptr);
            } else if (*ptr != NULLB) {
                bsp->id = MakeSeqID (ptr);
            }
            if (bsp->id == NULL) {
                if (trustID) {
                  bsp->id = MakeTrustedID (prefix, ctrptr);
                } else {
                  bsp->id = MakeNewProteinSeqIdExMT (NULL, NULL, prefix, ctrptr, TRUE);
                }
            }
            if (chptr != NULL) {
                if((vnp = SeqDescrNew(NULL)) != NULL) {
                    vnp->choice = Seq_descr_title;
                    while (IS_WHITESP(*chptr))
                        chptr++;
                    vnp->data.ptrvalue = StringSave (chptr);
                }
                bsp->descr = vnp;
            }
        } else {
            /* Unknown Seq-id */
            ptr = defline + 1;
            while (IS_WHITESP(*ptr))
                ptr++;
            if (StringNCmp (ptr, "unk100", 6) == 0) {
              bsp->id = MakeSeqID ("lcl|unk100");
              ptr += 3;
            } else {
              bsp->id = MakeSeqID ("lcl|gap");
            }
            bsp->repr = Seq_repr_virtual;
            if(*ptr != '\0' && sscanf(ptr, "%ld", &len) == 1 && len > 0) {
                bsp->length =  (Int4) len;
            } else {
                bsp->length = -1;
            }
            is_gap = TRUE;
        }
        MemFree(buffer);
    } else {  /* if ch == '>' EMPTY DEFLINE */
        /* Defline is upsent - creating default defline */
        if (trustID) {
          bsp->id = MakeTrustedID (prefix, ctrptr);
        } else {
          bsp->id = MakeNewProteinSeqIdExMT (NULL, NULL, prefix, ctrptr, TRUE);
        }
        if(type == FASTA_FILE_IO)
            ungetc(ch, fd);
    }
    SeqMgrAddToBioseqIndex (bsp);
    /* OK, now processing sequence */
    if (! is_gap) {
        if(!FastaReadSequenceInternalEx(input, type, next_char, is_na,
                                        &bsp->length,
                                        (ByteStorePtr PNTR) &(bsp->seq_data),
                                        errormsg, special_symbol,
                                        mask_ptr, bsp->id)) {
            ErrPostEx(SEV_FATAL, 0, 0, "Failure to read sequence. "
                      "FastaToSeqEntry() failed.\n");
            return NULL;
        }
    }
    BioseqPack(bsp);     /* Trying to pack Bioseq more */
    return sep;
}
NLM_EXTERN SeqEntryPtr FastaToSeqEntryInternalEx
(
 VoidPtr input,          /* input pointer (file or memory) */
 Int4 type,              /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR next_char, /* returned pointer to next FASTA sequence */
 Boolean is_na,          /* type of sequence */
 CharPtr PNTR errormsg,  /* error messge for debugging */
 Boolean parseSeqId,     /* Parse SeqID from def line */
 CharPtr special_symbol, /* Returns special symbol if no SeqEntry */
 CharPtr prefix,         /* prefix for localID if not parsable */
 Int2Ptr ctrptr,         /* starting point for constructing unique ID */
 SeqLocPtr PNTR mask_ptr /* Pointer to a SeqLoc to Fill with Masking information */
 )
{
  return FastaToSeqEntryInternalExEx (input, type, next_char, is_na, errormsg,
                                      parseSeqId, special_symbol, prefix, ctrptr,
                                      mask_ptr, FALSE);
}
NLM_EXTERN SeqEntryPtr FastaToSeqEntryInternal
(
 VoidPtr input,          /* input pointer (file or memory) */
 Int4 type,              /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR next_char, /* returned pointer to next FASTA sequence */
 Boolean is_na,          /* type of sequence */
 CharPtr PNTR errormsg,  /* error messge for debugging */
 Boolean parseSeqId,     /* Parse SeqID from def line */
 CharPtr special_symbol  /* Returns special symbol if no SeqEntry */
 )
{
    return FastaToSeqEntryInternalEx (input, type, next_char, is_na, errormsg,
                                          parseSeqId, special_symbol, NULL, NULL,NULL);
}
/*****************************************************************************
*
*   FastaId(bsp, buf, buflen)
*      Makes the string for the id part of fasta format.
*      buf should be at least 40 bytes
*
*****************************************************************************/
NLM_EXTERN Boolean FastaId(BioseqPtr bsp, CharPtr buf, Uint4 buflen)
{
    if ((bsp == NULL) || (buf == NULL)) return FALSE;
    SeqIdWrite(bsp->id, buf, PRINTID_FASTA_LONG, buflen);
    return TRUE;
}
static Boolean FastaIdX(BioseqPtr bsp, CharPtr buf, Uint4 buflen, Boolean printid_general, SeqLocPtr seqloc)
{
    Int4 length;
    if ((bsp == NULL) || (buf == NULL)) return FALSE;
    if (seqloc == NULL || SeqLocLen(seqloc) == bsp->length)
    { /* Full sequence is being dumped. */
        if (printid_general) {
            SeqIdWrite(bsp->id, buf, PRINTID_FASTA_GENERAL, buflen);
        } else {
            SeqIdWrite(bsp->id, buf, PRINTID_FASTA_LONG, buflen);
        }
    }
    else
    {
        SeqIdWrite(bsp->id, buf, PRINTID_FASTA_SHORT, buflen);
        length = StringLen(buf);
        sprintf(buf+length, ":%ld-%ld", (long) (SeqLocStart(seqloc)+1), (long) (SeqLocStop(seqloc)+1));
    }
    return TRUE;
}
/*****************************************************************************
*
*   FastaDefLine(bsp, buf, buflen, accession, organism)
*       Finds or makes a FASTA format defline (just locates the string)
*       buf should be very long if possible
*       function truncates if buf not long enough
*       a few deflines are longer than 255
*
*****************************************************************************/
NLM_EXTERN Boolean FastaDefLine (BioseqPtr bsp, CharPtr buf, Uint4 buflen,
                                          CharPtr accession, CharPtr organism, Uint1 tech)
{
    BioseqContextPtr bcp;
    ValNodePtr vnp;
    CharPtr tmp;
    PdbBlockPtr pbp;
    PatentSeqIdPtr psip;
    Uint4 diff, phase;
    Int4 num_segs, num_gaps;
    Char tbuf[80];
    static CharPtr htgs[2] = {
        "unordered", "ordered" };
    if ((bsp == NULL) || (buf == NULL)) return FALSE;
    buflen--;
    buf[buflen] = '\0';
    if (accession != NULL)
    {
        diff = LabelCopyExtra(buf, accession, buflen, "(", ") ");
        buflen -= diff;
        buf += diff;
    }
    bcp = BioseqContextNew(bsp);
    diff = 0;
    if ((tmp = BioseqContextGetTitle(bcp)) != NULL) {
        diff = LabelCopy(buf, tmp, buflen);
                                /* remove trailing blanks and periods */
        tmp = buf + diff - 1;   /* point at last character */
        while (((*tmp <= ' ') || (*tmp == '.')) && (diff))
        {
            *tmp = '\0';
            tmp--; diff--;
        }
    }
    else
        if ((vnp = BioseqContextGetSeqDescr(bcp, Seq_descr_pdb, NULL, NULL)) != NULL)
    {
        pbp = (PdbBlockPtr)(vnp->data.ptrvalue);
        diff = LabelCopy(buf, (CharPtr)(pbp->compound->data.ptrvalue), buflen);
    }
    else
    {
        for (vnp = bsp->id; vnp != NULL; vnp = vnp->next)
        {
            if (vnp->choice == SEQID_PATENT)
            {
                psip = (PatentSeqIdPtr)(vnp->data.ptrvalue);
                sprintf(tbuf, "Sequence %d from Patent %s %s",
                    (int)psip->seqid, psip->cit->country, psip->cit->number);
                diff = LabelCopy(buf, tbuf, buflen);
                break;
            }
        }
        if (vnp == NULL)
            diff = LabelCopy(buf, "No definition line found", buflen);
    }
    buflen -= diff;
    buf += diff;
    BioseqContextFree(bcp);
    if (((tech >= MI_TECH_htgs_1) && (tech <= MI_TECH_htgs_3)) ||
        (tech == MI_TECH_htgs_0))
    {
      if (tech == MI_TECH_htgs_0) {
        phase = 0;
        StringMove(tbuf, ", LOW-PASS SEQUENCE SAMPLING.");
      }
      else {
        phase = (Int2)(tech - MI_TECH_htgs_1 + 1);
        if (phase != 3)
          StringMove(tbuf, ", WORKING DRAFT SEQUENCE");
      }
      if (phase != 3) {
        diff = LabelCopy(buf, tbuf, buflen);
        buflen -= diff;
        buf += diff;
      }
        if (phase == 3)
        {
            if (tmp && StringStr(tmp, "complete sequence") == NULL) {
                diff = LabelCopy(buf, ", complete sequence", buflen);
                buflen -= diff;
                buf += diff;
            }
        }
        else if ((bsp->repr == Seq_repr_delta) && (phase != 0))
        {
            if (CountGapsInDeltaSeq(bsp, &num_segs, &num_gaps, NULL, NULL, NULL, 0))
            {
                if (num_gaps > 0) {
                    sprintf(tbuf, ", %ld %s pieces", (long)(num_gaps + 1), htgs[phase - 1]);
                } else {
                    sprintf(tbuf, ", %ld %s piece", (long)(num_gaps + 1), htgs[phase - 1]);
                }
                diff = LabelCopy(buf, tbuf, buflen);
                buflen -= diff;
                buf += diff;
            }
        }
    }
    if (organism != NULL)
    {
        LabelCopyExtra(buf, organism, buflen, " [", "]");
    }
    return TRUE;
}
static Boolean is_pdb(BioseqPtr bsp)
{
    SeqIdPtr id;
    if (bsp ==NULL)
        return FALSE;
    for (id = bsp->id; id; id=id->next)
    {
        if (id->choice == SEQID_PDB)
            return TRUE;
    }
    return FALSE;
}
static ValNodePtr tie_next(ValNodePtr head, ValNodePtr next)
{
   ValNodePtr v;
   if (head == NULL) {
      return next;
   }
   for (v = head; v->next != NULL; v = v->next)
           continue;
   v->next = next;
   return head;
}
static Boolean get_descr_on_top (GatherContextPtr gcp)
{
    ValNodePtr    tmp;
    DescrInfoPtr    PNTR dspp;
    DescrInfoPtr    dsp;
    ItemInfoPtr     iip;
    dspp = (DescrInfoPtr PNTR) gcp->userdata;
    dsp = *dspp;
    switch (gcp->thistype) {
    case OBJ_SEQDESC:
        tmp = (ValNodePtr) (gcp->thisitem);
        if (tmp->choice == dsp->choice) {
            if (tmp->data.ptrvalue != NULL) {
                dsp->vnp = tmp;
                iip = (ItemInfoPtr) MemNew(sizeof(ItemInfo));
                if(dsp->iip != NULL)
                    MemFree(dsp->iip);
                dsp->iip = iip;
                iip->entityID = gcp->entityID;
                iip->itemID = gcp->itemID;
                iip->itemtype = gcp->thistype;
            }
        }
        break;
    default:
        break;
    }
    return TRUE;
}
static Boolean get_descr (GatherContextPtr gcp)
{
    ValNodePtr    tmp;
    DescrInfoPtr    PNTR dspp;
    DescrInfoPtr    dsp;
    ItemInfoPtr     iip;
    BioseqPtr         bsp;
    dspp = (DescrInfoPtr PNTR) gcp->userdata;
    dsp = *dspp;
    switch (gcp->thistype)
    {
        case OBJ_SEQDESC:
            tmp = (ValNodePtr) (gcp->thisitem);
            if (tmp->choice == dsp->choice) {
                bsp = (BioseqPtr) (gcp->parentitem);
                if (dsp->bsp != bsp) {
                    break;
                }
                if (tmp->data.ptrvalue != NULL) {
                    dsp->vnp = tmp;
                    iip = (ItemInfoPtr) MemNew(sizeof(ItemInfo));
                    dsp->iip = iip;
                    iip->entityID = gcp->entityID;
                    iip->itemID = gcp->itemID;
                    iip->itemtype = gcp->thistype;
                }
            }
            break;
        default:
            break;
    }
    return TRUE;
}
static Boolean GetFeatProt (GatherContextPtr gcp)
{
    ValNodePtr    PNTR vnpp;
    ValNodePtr tmp;
    SeqFeatPtr    sfp;
    vnpp = (ValNodePtr PNTR) gcp->userdata;
    switch (gcp->thistype)
    {
        case OBJ_SEQFEAT:
            sfp = (SeqFeatPtr) (gcp->thisitem);
            if (sfp->data.choice == SEQFEAT_PROT) {
                tmp = ValNodeNew(NULL);
                tmp->data.ptrvalue = sfp;
                *vnpp = tie_next(*vnpp, tmp);
            }
            break;
        default:
            break;
    }
    return TRUE;
}
static Boolean GetFeatCDS (GatherContextPtr gcp)
{
    SeqFeatPtr    PNTR sfpp;
    SeqFeatPtr    sfp;
    sfpp = (SeqFeatPtr PNTR) gcp->userdata;
    switch (gcp->thistype)
    {
        case OBJ_SEQFEAT:
            sfp = (SeqFeatPtr) (gcp->thisitem);
            if (sfp->data.choice == SEQFEAT_CDREGION) {
                *sfpp = sfp;
                return FALSE;
            }
            break;
        default:
            break;
    }
    *sfpp = NULL;
    return TRUE;
}
static Boolean GetFeatGenes (GatherContextPtr gcp)
{
    ValNodePtr    PNTR vnpp;
    ValNodePtr tmp;
    SeqFeatPtr    sfp;
    vnpp = (ValNodePtr PNTR) gcp->userdata;
    switch (gcp->thistype)
    {
        case OBJ_SEQFEAT:
            sfp = (SeqFeatPtr) (gcp->thisitem);
            if (sfp->data.choice == SEQFEAT_GENE) {
                tmp = ValNodeNew(NULL);
                tmp->data.ptrvalue = sfp;
                *vnpp = tie_next(*vnpp, tmp);
            }
            break;
        default:
            break;
    }
    return TRUE;
}
static ValNodePtr IndexedGatherDescrOnBioseq (ItemInfoPtr iip, BioseqPtr bsp, Uint1 choice)
{
    SeqMgrDescContext  dcontext;
    SeqDescrPtr        sdp;
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, choice, &dcontext);
    if (sdp == NULL) return NULL;
    if (ISA_aa(bsp->mol) && !is_pdb(bsp)) {
        if (dcontext.level != 0) return NULL;
    }
    if (iip != NULL) {
        iip->entityID = dcontext.entityID;
        iip->itemID = dcontext.itemID;
        iip->itemtype = OBJ_SEQDESC;
    }
    return sdp;
}
static ValNodePtr GatherDescrOnBioseq(ItemInfoPtr iip, BioseqPtr bsp, Uint1 choice, Boolean get_first)
{
    ValNodePtr    vnp = NULL;
    /*
    GatherScope   gsc;
    SeqLocPtr     slp;
    Uint2         bspID;
    DescrInfoPtr  dsp;
    Uint2         entityID;
    */
    ObjValNodePtr ovp;
  if (ISA_aa(bsp->mol) && !is_pdb(bsp)) {
    vnp = BioseqGetSeqDescr (bsp, choice, NULL);
  } else {
    vnp = GetNextDescriptorUnindexed (bsp, choice, NULL);
  }
  if (vnp != NULL) {
    if (iip != NULL) {
      if (vnp->extended != 0) {
        ovp = (ObjValNodePtr) vnp;
        iip->entityID = ovp->idx.entityID;
        iip->itemtype = ovp->idx.itemtype;
        iip->itemID = ovp->idx.itemID;
      }
    }
  }
  return vnp;
#if 0
    entityID = ObjMgrGetEntityIDForPointer (bsp);
    if (SeqMgrFeaturesAreIndexed (entityID)) {
        return IndexedGatherDescrOnBioseq (iip, bsp, choice);
    }
    /*
    if (iip==NULL && (get_first || (ISA_aa(bsp->mol) && !is_pdb(bsp))) ) {
        for(vnp=bsp->descr;vnp && vnp->choice != choice; vnp=vnp->next){}
        return vnp;
    }
    */
    if (iip==NULL && get_first)
          {
            for(vnp=bsp->descr;vnp; vnp=vnp->next)
              if(vnp->choice == choice)
                return vnp;
          }
    dsp = (DescrInfoPtr) MemNew(sizeof(DescrInfo));
    dsp->choice = choice;
    dsp->bsp = bsp;
      MemSet ((Pointer) (&gsc), 0, sizeof (GatherScope));
    MemSet ((Pointer) (gsc.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
    gsc.ignore[OBJ_SEQDESC] = FALSE;
    bspID = ObjMgrGetEntityIDForPointer(bsp);
    slp = ValNodeNew(NULL);
    slp->choice = SEQLOC_WHOLE;
    slp->data.ptrvalue = (SeqIdPtr) SeqIdDup (SeqIdFindBest (bsp->id, 0));
    gsc.target = slp;
    if (ISA_aa(bsp->mol) && !is_pdb(bsp)) {
        GatherEntity(bspID, &dsp, get_descr, &gsc);
    } else {
        GatherEntity(bspID, &dsp, get_descr_on_top, &gsc);
    }
    SeqLocFree(slp);
    vnp = dsp->vnp;
    if (vnp && vnp->data.ptrvalue) {
        if (iip != NULL) {
            iip->entityID = dsp->iip->entityID;
            iip->itemID = dsp->iip->itemID;
            iip->itemtype = dsp->iip->itemtype;
        }
        MemFree(dsp->iip);
        MemFree(dsp);
        return vnp;
    }
    MemFree(dsp->iip);
    MemFree(dsp);
    return NULL;
#endif
}
/* more efficient versions of feature gather functions for protein defline */
typedef struct unidxfeatdata {
  SeqIdPtr    bspid;
  SeqLocPtr   loc;
  Int4        longest;
  Int4        shortest;
  SeqFeatPtr  sfp;
} UndxFeatData, PNTR UndxFeatPtr;
static void GetLongestProtFeat (
  SeqFeatPtr sfp,
  Pointer userdata
)
{
  Int4         len;
  SeqIdPtr     sip;
  UndxFeatPtr  ufp;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT) return;
  ufp = (UndxFeatPtr) userdata;
  if (ufp == NULL) return;
  sip = SeqLocId (sfp->location);
  if (sip == NULL) return;
  if (! SeqIdIn (sip, ufp->bspid)) return;
  len = SeqLocLen (sfp->location);
  if (len == -1) return;
  if (len > ufp->longest) {
    ufp->sfp = sfp;
    ufp->longest = len;
  }
}
static SeqFeatPtr GetLongestProteinUnindexed (
  BioseqPtr bsp
)
{
  BioseqSetPtr  bssp = NULL;
  UndxFeatData  ufd;
  if (bsp == NULL) return NULL;
  MemSet ((Pointer) &ufd, 0, sizeof (UndxFeatData));
  ufd.bspid = bsp->id;
  ufd.longest = 0;
  ufd.sfp = NULL;
  VisitFeaturesOnBsp (bsp, (Pointer) &ufd, GetLongestProtFeat);
  if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) bsp->idx.parentptr;
  }
  if (bssp != NULL && bssp->_class == BioseqseqSet_class_parts) {
    VisitFeaturesOnSet (bssp, (Pointer) &ufd, GetLongestProtFeat);
    if (bssp->idx.parenttype == OBJ_BIOSEQSET) {
      bssp = (BioseqSetPtr) bssp->idx.parentptr;
    }
  }
  if (bssp != NULL && bssp->_class == BioseqseqSet_class_segset) {
    VisitFeaturesOnSet (bssp, (Pointer) &ufd, GetLongestProtFeat);
  }
  return ufd.sfp;
}
static void GetCDSProtFeat (
  SeqFeatPtr sfp,
  Pointer userdata
)
{
  SeqIdPtr     sip;
  UndxFeatPtr  ufp;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;
  ufp = (UndxFeatPtr) userdata;
  if (ufp == NULL) return;
  sip = SeqLocId (sfp->product);
  if (sip == NULL) return;
  if (! SeqIdIn (sip, ufp->bspid)) return;
  ufp->sfp = sfp;
}
static SeqFeatPtr GetCDSProtUnindexed (
  BioseqPtr bsp
)
{
  Uint2         entityID;
  SeqEntryPtr   sep;
  UndxFeatData  ufd;
  if (bsp == NULL) return NULL;
  entityID = ObjMgrGetEntityIDForPointer (bsp);
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return NULL;
  MemSet ((Pointer) &ufd, 0, sizeof (UndxFeatData));
  ufd.bspid = bsp->id;
  ufd.sfp = NULL;
  VisitFeaturesInSep (sep, (Pointer) &ufd, GetCDSProtFeat);
  return ufd.sfp;
}
static void GetBestGeneFeat (
  SeqFeatPtr sfp,
  Pointer userdata
)
{
  Int4         diff;
  SeqIdPtr     sip;
  UndxFeatPtr  ufp;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE) return;
  ufp = (UndxFeatPtr) userdata;
  if (ufp == NULL) return;
  sip = SeqLocId (sfp->location);
  if (sip == NULL) return;
  if (! SeqIdIn (sip, ufp->bspid)) return;
  diff = SeqLocAinB (ufp->loc, sfp->location);
  if (diff >= 0) {
    if (diff < ufp->shortest) {
      ufp->sfp = sfp;
      ufp->shortest = diff;
    }
  }
}
static SeqFeatPtr GetBestGeneUnindexed (
  SeqLocPtr slp,
  Uint2 entityID
)
{
  BioseqPtr     bsp;
  SeqEntryPtr   sep;
  SeqIdPtr      sip;
  UndxFeatData  ufd;
  if (slp == NULL) return NULL;
  sip = SeqLocId (slp);
  if (sip == NULL) return NULL;
  bsp = BioseqFindCore (sip);
  if (bsp == NULL) return NULL;
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return NULL;
  MemSet ((Pointer) &ufd, 0, sizeof (UndxFeatData));
  ufd.bspid = bsp->id;
  ufd.loc = slp;
  ufd.shortest = INT4_MAX;
  ufd.sfp = NULL;
  VisitFeaturesInSep (sep, (Pointer) &ufd, GetBestGeneFeat);
  return ufd.sfp;
}
/* GatherProtCDS is still faster than GetCDSProtUnindexed for some reason */
static SeqFeatPtr GatherProtCDS(BioseqPtr bsp)
{
    GatherScope gsc;
    SeqLocPtr slp = NULL;
    Uint2 bspID;
    SeqFeatPtr sfp;
      MemSet ((Pointer) (&gsc), 0, sizeof (GatherScope));
    MemSet ((Pointer) (gsc.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
    gsc.ignore[OBJ_SEQFEAT] = FALSE;
    gsc.ignore[OBJ_SEQANNOT] = FALSE;
    gsc.get_feats_product = TRUE;
    bspID = ObjMgrGetEntityIDForPointer(bsp);
    slp = ValNodeNew(NULL);
    slp->choice = SEQLOC_WHOLE;
    slp->data.ptrvalue = (SeqIdPtr) SeqIdDup (SeqIdFindBest (bsp->id, 0));
    gsc.target = slp;
    sfp = NULL;
    GatherEntity(bspID, &sfp, GetFeatCDS, &gsc);
    SeqLocFree(slp);
    return sfp;
}
/* obsolete functions, replaced by Unindexed versions */
static SeqFeatPtr GatherSeqFeatProt(BioseqPtr bsp)
{
    GatherScope gsc;
    SeqLocPtr slp = NULL;
    Uint2 bspID;
    SeqFeatPtr sfp = NULL;
    SeqFeatPtr f;
    ValNodePtr prot, v;
    Int4 length, longest_length=0;
      MemSet ((Pointer) (&gsc), 0, sizeof (GatherScope));
    MemSet ((Pointer) (gsc.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
    gsc.ignore[OBJ_SEQFEAT] = FALSE;
    gsc.ignore[OBJ_SEQANNOT] = FALSE;
    gsc.get_feats_location = TRUE;
    bspID = ObjMgrGetEntityIDForPointer(bsp);
    slp = ValNodeNew(NULL);
    slp->choice = SEQLOC_WHOLE;
    slp->data.ptrvalue = (SeqIdPtr) SeqIdDup (SeqIdFindBest (bsp->id, 0));
    gsc.target = slp;
    prot = NULL;
    GatherEntity(bspID, &prot, GetFeatProt, &gsc);
    for (v=prot; v; v=v->next) {
        f = (SeqFeatPtr) v->data.ptrvalue;
        if ((length=SeqLocLen(f->location)) == -1)
            continue;
        if (length > longest_length) {
            sfp = f;
            longest_length = length;
        }
    }
    ValNodeFree(prot);
    SeqLocFree(slp);
    return sfp;
}
static ValNodePtr GatherGenesForCDS(SeqLocPtr slp)
{
    GatherScope gsc;
    Uint2 bspID;
    ValNodePtr vnp;
    BioseqPtr bsp;
    bsp = BioseqFindCore(SeqLocId(slp));
    if (bsp == NULL)
        return NULL;
    bspID = ObjMgrGetEntityIDForPointer(bsp);
      MemSet ((Pointer) (&gsc), 0, sizeof (GatherScope));
    MemSet ((Pointer) (gsc.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
    gsc.ignore[OBJ_SEQFEAT] = FALSE;
    gsc.ignore[OBJ_SEQANNOT] = FALSE;
    gsc.get_feats_location = TRUE;
    gsc.target = slp;
    vnp = NULL;
    GatherEntity(bspID, &vnp, GetFeatGenes, &gsc);
    return vnp;
}
typedef struct nmdef {
  SeqFeatPtr  gene;
  SeqFeatPtr  cds;
  SeqFeatPtr  prot;
  Int4        protlen;
  Int2        numgenes;
  Int2        numcds;
  Int2        numprots;
} NMDef, PNTR NMDefPtr;
static void FindNMFeats (SeqFeatPtr sfp, Pointer userdata)
{
  Int4      len;
  NMDefPtr  ndp;
  if (sfp == NULL) return;
  ndp = (NMDefPtr) userdata;
  if (ndp == NULL) return;
  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      ndp->gene = sfp;
      (ndp->numgenes)++;
      break;
    case SEQFEAT_CDREGION :
      ndp->cds = sfp;
      (ndp->numcds++);
      break;
    case SEQFEAT_PROT :
      len = SeqLocLen (sfp->location);
      if (len > ndp->protlen) {
        ndp->prot = sfp;
        ndp->protlen = len;
        (ndp->numprots)++;
      }
      break;
    default :
      break;
  }
}
static Boolean IsFlyCG (CharPtr str)
{
  Char  ch;
  if (StringHasNoText (str)) return FALSE;
  ch = *str;
  if (ch != 'C') return FALSE;
  str++;
  ch = *str;
  if (ch != 'G') return FALSE;
  str++;
  ch = *str;
  while (IS_DIGIT (ch)) {
    str++;
    ch = *str;
  }
  if (ch != '-') return FALSE;
  str++;
  ch = *str;
  if (ch != 'P') return FALSE;
  str++;
  ch = *str;
  if (IS_ALPHA (ch)) {
    str++;
    ch = *str;
    if (ch == '\0' || ch == ' ' || ch == ',' || ch == ';') return TRUE;
  }
  return FALSE;
}
static void ReplaceFlyDashPwithDashR (CharPtr str)
{
  Char     ch;
  CharPtr  ptr;
  while (StringDoesHaveText (str)) {
    ch = *str;
    while (IS_WHITESP (ch)) {
      str++;
      ch = *str;
    }
    if (IsFlyCG (str)) {
      ptr = StringStr (str, "-P");
      if (ptr != NULL) {
        ptr [1] = 'R';
        return;
      }
    }
    while (ch != '\0' && (! IS_WHITESP (ch))) {
      str++;
      ch = *str;
    }
  }
}
static CharPtr FindNMDefLine (BioseqPtr bsp)
{
  BioSourcePtr  biop;
  Char          buf [512], buf2 [600];
  CharPtr       cds = NULL;
  Uint2         entityID;
  CharPtr       gene;
  Boolean       is_refseq = FALSE;
  size_t        len;
  NMDef         nd;
  OrgRefPtr     orp;
  CharPtr       ptr;
  SeqEntryPtr   sep;
  SeqIdPtr      sip;
  CharPtr       str;
  ValNodePtr    vnp;
  MemSet ((Pointer) &nd, 0, sizeof (NMDef));
  entityID = ObjMgrGetEntityIDForPointer (bsp);
  sep = GetBestTopParentForDataEx (entityID, bsp, TRUE);
  VisitFeaturesInSep (sep, (Pointer) &nd, FindNMFeats);
  if (nd.numgenes != 1 || nd.numcds != 1 || nd.numprots < 1) return NULL;
  vnp = GatherDescrOnBioseq (NULL, bsp, Seq_descr_source, FALSE);
  if (vnp == NULL) return NULL;
  biop = (BioSourcePtr) vnp->data.ptrvalue;
  orp = biop->org;
  if (orp == NULL || StringHasNoText (orp->taxname)) return NULL;
  FeatDefLabel (nd.gene, buf, sizeof (buf) - 1, OM_LABEL_CONTENT);
  gene = StringSaveNoNull (buf);
  FeatDefLabel (nd.cds, buf, sizeof (buf) - 1, OM_LABEL_CONTENT);
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      is_refseq = TRUE;
    }
  }
  if (is_refseq) {
    /* special case Drosophila RefSeq NM titles */
    if (StringICmp (orp->taxname, "Drosophila melanogaster") == 0) {
      ReplaceFlyDashPwithDashR (buf);
    }
    ptr = StringStr (buf, "isoform ");
    if (ptr != NULL) {
      *ptr = '\0';
      ptr += 8;
      StringCpy (buf2, buf);
      StringCat (buf2, "transcript variant ");
      StringCat (buf2, ptr);
      cds = StringSaveNoNull (buf2);
    } else {
      cds = StringSaveNoNull (buf);
    }
  } else {
    cds = StringSaveNoNull (buf);
  }
  len = StringLen (orp->taxname) + StringLen (cds) +
        StringLen (gene) + StringLen ("  (), mRNA") + 10;
  str = (CharPtr) MemNew (len);
  if (str != NULL) {
    sprintf (str, "%s %s (%s), mRNA", orp->taxname, cds, gene);
  }
  MemFree (gene);
  MemFree (cds);
  return str;
}
static CharPtr FindNRDefLine (BioseqPtr bsp)
{
  BioSourcePtr  biop;
  Char          buf [512];
  Uint2         entityID;
  CharPtr       gene;
  size_t        len;
  MolInfoPtr    mip;
  NMDef         nd;
  OrgRefPtr     orp;
  CharPtr       rna = "miscRNA";
  SeqEntryPtr   sep;
  CharPtr       str;
  ValNodePtr    vnp;
  MemSet ((Pointer) &nd, 0, sizeof (NMDef));
  entityID = ObjMgrGetEntityIDForPointer (bsp);
  sep = GetBestTopParentForDataEx (entityID, bsp, TRUE);
  VisitFeaturesInSep (sep, (Pointer) &nd, FindNMFeats);
  if (nd.numgenes < 1) return NULL;
  vnp = GatherDescrOnBioseq (NULL, bsp, Seq_descr_source, FALSE);
  if (vnp == NULL) return NULL;
  biop = (BioSourcePtr) vnp->data.ptrvalue;
  orp = biop->org;
  if (orp == NULL || StringHasNoText (orp->taxname)) return NULL;
  FeatDefLabel (nd.gene, buf, sizeof (buf) - 1, OM_LABEL_CONTENT);
  gene = StringSaveNoNull (buf);
  vnp = GatherDescrOnBioseq (NULL, bsp, Seq_descr_molinfo,TRUE);
  if (vnp != NULL) {
    mip = (MolInfoPtr) vnp->data.ptrvalue;
    if (mip != NULL) {
      switch (mip->biomol) {
        case MOLECULE_TYPE_PRE_MRNA :
          rna = "precursorRNA";
          break;
        case MOLECULE_TYPE_MRNA :
          rna = "mRNA";
          break;
        case MOLECULE_TYPE_RRNA :
          rna = "rRNA";
          break;
        case MOLECULE_TYPE_TRNA :
          rna = "tRNA";
          break;
        case MOLECULE_TYPE_SNRNA :
          rna = "snRNA";
          break;
        case MOLECULE_TYPE_SCRNA :
          rna = "scRNA";
          break;
        case MOLECULE_TYPE_CRNA :
          rna = "cRNA";
          break;
        case MOLECULE_TYPE_SNORNA :
          rna = "snoRNA";
          break;
        case MOLECULE_TYPE_TRANSCRIBED_RNA :
          rna = "miscRNA";
          break;
        case MOLECULE_TYPE_NCRNA :
          rna = "ncRNA";
          break;
        case MOLECULE_TYPE_TMRNA :
          rna = "tmRNA";
          break;
        default :
          break;
      }
    }
  }
  len = StringLen (orp->taxname) + StringLen (gene) +
        StringLen (", ") + 30;
  str = (CharPtr) MemNew (len);
  if (str != NULL) {
    sprintf (str, "%s %s, %s", orp->taxname, gene, rna);
  }
  MemFree (gene);
  return str;
}

static CharPtr TrimPunctuationFromEnd (CharPtr str)

{
  Uchar    ch;      /* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == ' ' || ch == ';' || ch == ',' || ch == '~' || ch == '.') {
        if (dst == NULL) {
          dst = ptr;
        }
      } else  {
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

static CharPtr TrimNonPeriodPunctuationFromEnd (CharPtr str)

{
  Uchar    ch;      /* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == ' ' || ch == ';' || ch == ',' || ch == '~') {
        if (dst == NULL) {
          dst = ptr;
        }
      } else  {
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

static CharPtr FindProtDefLine(BioseqPtr bsp, Boolean extProtTitle)
{
    SeqFeatPtr sfp = NULL /* , f */;
    ProtRefPtr prp;
    SeqFeatXrefPtr xref;
    GeneRefPtr grp=NULL;
    ValNodePtr vnp, /* v, */ syn;
    SeqLocPtr loc;
    CharPtr title = NULL, s, geneprod;
    /*
    Int4 diff_lowest = INT4_MAX, diff_current;
    */
    Int2 length = 0;
    SeqFeatPtr best_gene = NULL;
    Uint2 entityID;
    Boolean indexed;
    if (bsp == NULL) {
        return NULL;
    }
    entityID = ObjMgrGetEntityIDForPointer (bsp);
    indexed = (Boolean)SeqMgrFeaturesAreIndexed (entityID);
    sfp = NULL;
    if (indexed) {
        sfp = SeqMgrGetBestProteinFeature (bsp, NULL);
    } else {
        sfp = GetLongestProteinUnindexed (bsp);
        /*
        if (sfp == NULL) {
            sfp = GatherSeqFeatProt(bsp);
        }
        */
    }
    if (sfp != NULL) {
        prp = (ProtRefPtr) sfp->data.value.ptrvalue;
        if (prp && prp->name) {
            for (vnp=prp->name; vnp; vnp=vnp->next) {
                length += StringLen((CharPtr)vnp->data.ptrvalue) + 2;
            }
            s = title = (CharPtr) MemNew(length + 1);
            if (prp->name->data.ptrvalue) {
                sprintf(title, "%s",
                                        (CharPtr) prp->name->data.ptrvalue);
            }
            s += StringLen(title);
            if (extProtTitle) {
                for (vnp=prp->name->next; vnp; vnp=vnp->next) {
                    sprintf(s, "; %s",
                                            (CharPtr) vnp->data.ptrvalue);
                    s += StringLen((CharPtr)vnp->data.ptrvalue) + 2;
                }
            }
            TrimPunctuationFromEnd (title);
            /* if hypothetical protein, append locus_tag */
            if (StringICmp (title, "hypothetical protein") == 0) {
                sfp = NULL;
                if (indexed) {
                    sfp = SeqMgrGetCDSgivenProduct (bsp, NULL);
                } else {
                    /*
                    sfp = GetCDSProtUnindexed (bsp);
                    */
                    sfp = GatherProtCDS(bsp);
                }
                if (sfp != NULL) {
                    grp = SeqMgrGetGeneXref (sfp);
                    if (grp == NULL) {
                        loc = sfp->location;
                        best_gene = NULL;
                        if (indexed) {
                            best_gene = SeqMgrGetOverlappingGene (loc, NULL);
                        } else {
                            best_gene = GetBestGeneUnindexed (loc, entityID);
                            /*
                            vnp = GatherGenesForCDS(loc);
                            for (v=vnp; v; v=v->next) {
                                f = (SeqFeatPtr) v->data.ptrvalue;
                                diff_current = SeqLocAinB(loc, f->location);
                                if (! diff_current) {
                                    best_gene = f;
                                    break;
                                } else if (diff_current > 0) {
                                    if ((diff_lowest == -1) || (diff_current<diff_lowest)) {
                                        diff_lowest = diff_current;
                                        best_gene = f;
                                    }
                                }
                            }
                            ValNodeFree(vnp);
                            */
                        }
                        if (best_gene != NULL) {
                            grp = (GeneRefPtr) best_gene->data.value.ptrvalue;
                        }
                    }
                }
                if (grp != NULL) {
                    geneprod = NULL;
                    if (grp->locus_tag != NULL) {
                        geneprod = grp->locus_tag;
                    }
                    if (geneprod != NULL) {
                        s = (CharPtr) MemNew (StringLen (geneprod) + StringLen (title) + 20);
                        if (s != NULL) {
                            sprintf (s, "%s %s", title, geneprod);
                            MemFree (title);
                            title = s;
                        }
                    }
                }
            }
        } else if (prp && prp->desc) {
            title = StringSave(prp->desc);
        } else if (prp && prp->activity) {
            if (prp->activity->data.ptrvalue) {
                title = StringSave (prp->activity->data.ptrvalue);
            }
        }
    }
    if (title == NULL) {
        sfp = NULL;
        if (indexed) {
            sfp = SeqMgrGetCDSgivenProduct (bsp, NULL);
        } else {
            /*
            sfp = GetCDSProtUnindexed (bsp);
            */
            sfp = GatherProtCDS(bsp);
        }
        if (sfp != NULL) {
            loc = sfp->location;
            for (xref = sfp->xref; xref; xref=xref->next) {
                if (xref->data.choice == SEQFEAT_GENE) {
                    grp = (GeneRefPtr) xref->data.value.ptrvalue;
                }
            }
            if (grp) {
                geneprod = NULL;
                if (grp->locus != NULL) {
                    geneprod = grp->locus;
                } else if (grp->syn != NULL) {
                    syn = grp->syn;
                    geneprod = (CharPtr) syn->data.ptrvalue;
                } else if (grp->desc != NULL) {
                    geneprod = (CharPtr) grp->desc;
                }
                if (geneprod != NULL) {
                    s = (CharPtr) MemNew(StringLen(geneprod) + 15);
                    sprintf(s, "%s gene product", geneprod);
                    title = s;
                }
            }
            if (title == NULL) {
                best_gene = NULL;
                if (indexed) {
                    best_gene = SeqMgrGetOverlappingGene (loc, NULL);
                } else {
                    best_gene = GetBestGeneUnindexed (loc, entityID);
                    /*
                    vnp = GatherGenesForCDS(loc);
                    for (v=vnp; v; v=v->next) {
                        f = (SeqFeatPtr) v->data.ptrvalue;
                        diff_current = SeqLocAinB(loc, f->location);
                        if (! diff_current) {
                            best_gene = f;
                            break;
                        } else if (diff_current > 0) {
                            if ((diff_lowest == -1) || (diff_current<diff_lowest)) {
                                diff_lowest = diff_current;
                                best_gene = f;
                            }
                        }
                    }
                    ValNodeFree(vnp);
                    */
                }
                if (best_gene != NULL) {
                    grp = (GeneRefPtr) best_gene->data.value.ptrvalue;
                    if (grp) {
                        geneprod = NULL;
                        if (grp->locus != NULL) {
                            geneprod = grp->locus;
                        } else if (grp->syn != NULL) {
                            syn = grp->syn;
                            geneprod = (CharPtr) syn->data.ptrvalue;
                        } else if (grp->desc != NULL) {
                            geneprod = (CharPtr) grp->desc;
                        }
                        if (geneprod != NULL) {
                            s = (CharPtr) MemNew(StringLen(geneprod) + 15);
                            sprintf(s, "%s gene product", geneprod);
                            title = s;
                        }
                    }
                }
            }
        }
    }
    if (title != NULL) {
      TrimPunctuationFromEnd (title);
    }
    if (title == NULL) {
      title = StringSave ("unnamed protein product");
    }
    return title;
}
static Boolean StrainNotAtEndOfTaxname (CharPtr name, CharPtr strain)
{
  size_t   len;
  CharPtr  ptr;
  char ch;

  if (StringHasNoText (name) || StringHasNoText (strain)) return TRUE;
  ptr = StringChr (name, ' ');
  if (ptr == NULL) return TRUE;
  ptr++;
  ptr = StringChr (ptr, ' ');
  if (ptr == NULL) return TRUE;
  ptr++;
  ptr = StringISearch (ptr, strain);
  if (ptr == NULL) return TRUE;
  len = StringLen (strain);
  ptr += len;
  if (! StringHasNoText (ptr)) {
    if (StringCmp (ptr, "'") == 0) {
      ptr -= len + 1;
      if (*ptr == '\'') return FALSE;
    }
    return TRUE;
  }
  ptr -= len + 1;
  ch = *ptr;
  if (ch != ' ' && ch != '-' && ch != ':' && ch != ';' && ch != '.') return TRUE;
  return FALSE;
}
static Int2 GetNumClones (CharPtr str)
{
  Char  ch;
  Int2  count;
  if (StringHasNoText (str)) return 0;
  count = 1;
  ch = *str;
  while (ch != '\0') {
    if (ch == ';') {
      count++;
    }
    str++;
    ch = *str;
  }
  return count;
}
static CharPtr SimpleSegSeqTitle (BioseqPtr bsp)
{
  BioSourcePtr       biop;
  SeqMgrFeatContext  ccontext;
  SeqFeatPtr         cds;
  CharPtr            clone = NULL;
  CharPtr            complete = "gene, complete cds";
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  gcontext;
  SeqFeatPtr         gene;
  GeneRefPtr         grp;
  CharPtr            isolate = NULL;
  CharPtr            label = NULL;
  size_t             len;
  CharPtr            locus = NULL;
  OrgModPtr          mod;
  CharPtr            modifier = NULL;
  Int2               numclones;
  ObjMgrDataPtr      omdp;
  ObjMgrPtr          omp;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  CharPtr            organism = NULL;
  CharPtr            product = NULL;
  SeqDescrPtr        sdp;
  SubSourcePtr       ssp;
  CharPtr            str;
  CharPtr            strain = NULL;
  CharPtr            title;
  ValNodePtr         vnp;
  Uint2              entityID;

  if (bsp == NULL) return NULL;
  /* check to see if feature indexing has been called */
  omdp = (ObjMgrDataPtr) bsp->omdp;
  if (omdp == NULL) return NULL;
  omp = ObjMgrReadLock ();
  omdp = ObjMgrFindTop (omp, omdp);
  ObjMgrUnlock ();
  if (omdp == NULL) return NULL;
  /*
  if (omdp->indexed == 0) return NULL;
  */

  entityID = ObjMgrGetEntityIDForPointer (bsp);
  if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
    SeqMgrIndexFeatures (entityID, NULL);
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL) return NULL;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return NULL;
  orp = biop->org;
  if (orp != NULL && (! StringHasNoText (orp->taxname))) {
    organism = orp->taxname;
    onp = orp->orgname;
    if (onp != NULL) {
      mod = onp->mod;
      if (mod != NULL) {
        if (mod->subtype == ORGMOD_strain) {
          if (mod->subname != NULL && StrainNotAtEndOfTaxname (organism, mod->subname)) {
            strain = (CharPtr) mod->subname;
          }
        } else if (mod->subtype == ORGMOD_isolate) {
          isolate = (CharPtr) mod->subname;
        }
      }
    }
  } else {
    organism = "Unknown";
  }
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_clone) {
      if (ssp->name != NULL) {
        numclones = GetNumClones (ssp->name);
        if (numclones < 4) {
          clone = (CharPtr) ssp->name;
        }
      }
    }
  }
  cds = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &ccontext);
  if (cds != NULL) {
    if (cds->partial) {
      complete = "gene, partial cds";
    }
    product = ccontext.label;
    grp = SeqMgrGetGeneXref (cds);
    if (grp != NULL) {
      if (! StringHasNoText (grp->locus)) {
        locus = grp->locus;
      } else {
        vnp = grp->syn;
        if (vnp != NULL) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (! StringHasNoText (str)) {
            locus = str;
          }
        }
      }
    }
    if (locus == NULL) {
      gene = SeqMgrGetOverlappingGene (cds->location, &gcontext);
      if (gene != NULL) {
        locus = gcontext.label;
      }
    }
  } else {
    if (StringDoesHaveText (strain)) {
      modifier = strain;
      label = "strain";
    } else if (StringDoesHaveText (clone)) {
      modifier = clone;
      label = "clone";
    } else if (StringDoesHaveText (isolate)) {
      modifier = isolate;
      label = "isolate";
    }
  }
  len = StringLen (organism) + StringLen (label) + StringLen (modifier) +
        StringLen (product) + StringLen (locus) + StringLen (complete);
  title = (CharPtr) MemNew (len + 10);
  if (organism != NULL) {
    StringCat (title, organism);
  }
  if (modifier != NULL) {
    StringCat (title, " ");
    StringCat (title, label);
    StringCat (title, " ");
    StringCat (title, modifier);
  }
  if (product != NULL) {
    StringCat (title, " ");
    StringCat (title, product);
  }
  if (locus != NULL) {
    StringCat (title, " (");
    StringCat (title, locus);
    StringCat (title, ")");
  }
  if (product != NULL || locus != NULL) {
    StringCat (title, " ");
    StringCat (title, complete);
  }
  TrimSpacesAroundString (title);
  return title;
}

static CharPtr UseOrgMods(BioseqPtr bsp, CharPtr suffix, Uint1 tech, Boolean htgs_pooled_multiclone)
{
    ItemInfoPtr         iip = NULL;
    ValNodePtr             vnp;
    BioSourcePtr         biop;
    OrgModPtr       mod;
    OrgNamePtr      onp;
    OrgRefPtr       orp;
    SubSourcePtr       ssp;
    Char            ch;
    CharPtr                  name = NULL, chr = NULL, str = NULL,
                              cln = NULL, map = NULL, pls = NULL, def = NULL, ptr;
    Int2                     deflen = 0;
    Int2            numclones;
    if (bsp == NULL) {
        return NULL;
    }
    if ((vnp=GatherDescrOnBioseq(iip, bsp, Seq_descr_source,FALSE)) == NULL) {
        return NULL;
    }
    biop = (BioSourcePtr) vnp->data.ptrvalue;
  orp = biop->org;
  if (orp && orp->taxname) {
    name = StringSave(orp->taxname);
    deflen += StringLen(orp->taxname);
  }
    for (ssp = biop->subtype; ssp; ssp=ssp->next) {
        if (ssp->subtype == SUBSRC_chromosome) { /* chromosome */
            if (ssp->name != NULL) {
                chr = (CharPtr) MemNew(StringLen(ssp->name) + 13);
                deflen += StringLen(ssp->name) + 13;
                sprintf(chr, " chromosome %s", ssp->name);
            }
        } else if (ssp->subtype == SUBSRC_clone) { /* clone */
            if (ssp->name != NULL) {
                numclones = GetNumClones (ssp->name);
                if (htgs_pooled_multiclone) {
                    cln = (CharPtr) MemNew (30);
                    sprintf (cln, ", pooled multiple clones");
                    deflen += StringLen (cln) + 2;
                } else if (numclones > 3) {
                    cln = (CharPtr) MemNew (20);
                    sprintf (cln, ", %d clones", (int) numclones);
                    deflen += StringLen (cln) + 2;
                } else {
                    cln = (CharPtr) MemNew(StringLen(ssp->name) + 8);
                    deflen += StringLen(ssp->name) + 8;
                    sprintf(cln, " clone %s", ssp->name);
                }
            }
        } else if (ssp->subtype == SUBSRC_map) { /* map */
            if (ssp->name != NULL) {
                map = (CharPtr) MemNew(StringLen(ssp->name) + 7);
                deflen += StringLen(ssp->name) + 7;
                sprintf(map, " map %s", ssp->name);
            }
        } else if (ssp->subtype == SUBSRC_plasmid_name) { /* plasmid name */
            if (ssp->name != NULL) {
                pls = (CharPtr) MemNew(StringLen(ssp->name) + 10);
                deflen += StringLen(ssp->name) + 10;
                sprintf(pls, " plasmid %s", ssp->name);
            }
        }
    }
  if (orp != NULL) {
        onp = orp->orgname;
        if (onp != NULL) {
            for (mod = onp->mod; mod != NULL; mod = mod->next) {
                if (mod->subtype != ORGMOD_strain) continue; /* strain */
                if (StringDoesHaveText (str)) continue;
                if (mod->subname != NULL && StrainNotAtEndOfTaxname (name, mod->subname)) {
                    str = (CharPtr) MemNew(StringLen(mod->subname) + 9);
                    deflen += StringLen(mod->subname) + 9;
                    sprintf(str, " strain %s", mod->subname);
                    ptr = StringChr (str, ';');
                    if (ptr != NULL) {
                      *ptr = '\0';
                    }
                    TrimNonPeriodPunctuationFromEnd (str);
                }
            }
        }
    }
    deflen += StringLen (suffix) + 2;
    def = (CharPtr) MemNew(deflen+1);
    if (def == NULL) return NULL;
    if (name) {
        def = StringCat(def, name);
        MemFree(name);
    }
    if (str) {
        def = StringCat(def, str);
        MemFree(str);
    }
    if (chr) {
        def = StringCat(def, chr);
        MemFree(chr);
    }
    if (cln) {
        def = StringCat(def, cln);
        MemFree(cln);
    }
    if (map) {
        def = StringCat(def, map);
        MemFree(map);
    }
  if (pls) {
    if (tech == MI_TECH_wgs) {
          def = StringCat(def, pls);
    }
        MemFree(pls);
  }
    if (suffix) {
        def = StringCat(def, " ");
        def = StringCat(def, suffix);
    }
    TrimSpacesAroundString (def);
    ch = def [0];
    def [0] = TO_UPPER (ch);
    return def;
}

/*
   The following lists need endogenous virus, hydrogenosome, chromosome, and chromatophore
*/

static CharPtr organelleByItself [] = {
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
  "provirus",
  "virus",
  "nucleomorph",
  "apicoplast",
  "leucoplast",
  "protoplast",
  NULL,
  NULL,
  NULL,
  NULL
};
static CharPtr organelleWithPlasmid [] = {
  NULL,
  NULL,
  "chloroplast",
  "chromoplast",
  "kinetoplast",
  "mitochondrial",
  "plastid",
  "macronuclear",
  "extrachrom",
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
  NULL,
  NULL,
  NULL,
  NULL
};
static CharPtr organelleForWGS [] = {
  NULL,
  NULL,
  "chloroplast",
  "chromoplast",
  "kinetoplast",
  "mitochondrial",
  "plastid",
  "",
  "",
  "",
  "",
  "",
  "cyanelle",
  "proviral",
  "virus",
  "",
  "apicoplast",
  "leucoplast",
  "proplastid",
  "endogenous virus",
  "hydrogenosome",
  "chromosome",
  "chromatophore"
};

const Int4 kNumWGSOrganelles = sizeof (organelleForWGS) / sizeof (CharPtr);


static void LowercasePlasmidOrElement (CharPtr def)
{
  CharPtr  ptr;
  if (StringHasNoText (def)) return;
  def++;
  ptr = StringISearch (def, "plasmid");
  while (ptr != NULL) {
    if (*ptr == 'P') {
      *ptr = 'p';
    }
    ptr = StringISearch (ptr + 7, "plasmid");
  }
  ptr = StringISearch (def, "element");
  while (ptr != NULL) {
    if (*ptr == 'E') {
      *ptr = 'e';
    }
    ptr = StringISearch (ptr + 7, "element");
  }
}


NLM_EXTERN CharPtr MakeCompleteChromTitle (BioseqPtr bsp, Uint1 biomol, Uint1 completeness)
{
    CharPtr       completeseq = ", complete sequence";
    CharPtr       completegen = ", complete genome";
    ItemInfoPtr   iip = NULL;
    ValNodePtr    vnp;
    BioSourcePtr  biop;
    OrgRefPtr     orp;
    SubSourcePtr  ssp;
    CharPtr       name = NULL, chr = NULL, orgnl = NULL,
                  seg = NULL, pls = NULL, def = NULL;
    Int2          deflen = 80; /* starts with space for all fixed text */
    Char          ch;
    Boolean       plasmid;
    Uint1         genome;
    if (bsp == NULL) {
        return NULL;
    }
    if ((vnp=GatherDescrOnBioseq(iip, bsp, Seq_descr_source,TRUE)) == NULL) {
        return NULL;
    }
    biop = (BioSourcePtr) vnp->data.ptrvalue;
    if (biop == NULL) {
        return NULL;
    }
        orp = biop->org;
        if (orp == NULL || orp->taxname == NULL) {
          return NULL;
        }
    name = orp->taxname;
    deflen += StringLen(orp->taxname);
    genome = biop->genome;
    plasmid = (Boolean) (biop->genome == GENOME_plasmid);
    for (ssp = biop->subtype; ssp; ssp=ssp->next) {
        if (ssp->subtype == SUBSRC_chromosome) {
            if (ssp->name != NULL) {
                chr = ssp->name;
                deflen += StringLen(ssp->name);
            }
        } else if (ssp->subtype == SUBSRC_segment) {
            if (ssp->name != NULL) {
                seg = ssp->name;
                deflen += StringLen(ssp->name);
            }
        } else if (ssp->subtype == SUBSRC_plasmid_name) {
            if (ssp->name != NULL) {
                pls = ssp->name;
                deflen += StringLen(ssp->name);
            }
        }
    }
    if (genome < kNumWGSOrganelles) {
        if (pls != NULL) {
            orgnl = organelleWithPlasmid [genome];
        } else {
            orgnl = organelleByItself [genome];
        }
        if (StringISearch (name, "virus") != NULL || StringISearch (name, "phage") != NULL) {
            if (genome == GENOME_proviral || genome == GENOME_virion) {
                orgnl = NULL;
            }
        }
    }
    if (completeness == 2 ||
        completeness == 3 ||
        completeness == 4 ||
        completeness == 5) {
        /* remove "complete" component */
        completeseq = ", partial sequence";
        completegen = ", genome";
    }
    def = (CharPtr) MemNew(deflen+1);
    if (StringISearch (name, "plasmid") != NULL) {
        StringCat(def, name);
        StringCat (def, completeseq);
        ch = *def;
        *def = TO_UPPER (ch);
        LowercasePlasmidOrElement (def);
        return def;
    } else if (plasmid) {
        if (name && (! pls)) {
            StringCat (def, name);
            StringCat (def, " unnamed plasmid");
            StringCat (def, completeseq);
            ch = *def;
            *def = TO_UPPER (ch);
            return def;
        }
        if (pls) {
            if (name) {
                StringCat (def, name);
                StringCat (def, " ");
            }
            if (StringISearch (pls, "plasmid") == NULL && StringISearch (pls, "element") == NULL) {
                StringCat(def, "plasmid ");
            }
            StringCat (def, pls);
            StringCat (def, completeseq);
            ch = *def;
            *def = TO_UPPER (ch);
            LowercasePlasmidOrElement (def);
            return def;
        }
    } else if (pls) {
        if (name) {
             StringCat (def, name);
            StringCat (def, " ");
        }
        if (orgnl != NULL) {
            StringCat (def, orgnl);
            StringCat (def, " ");
        }
        if (StringISearch (pls, "plasmid") == NULL && StringISearch (pls, "element") == NULL) {
            StringCat (def, "plasmid ");
        }
        StringCat (def, pls);
        StringCat (def, completeseq);
        ch = *def;
        *def = TO_UPPER (ch);
        LowercasePlasmidOrElement (def);
        return def;
    } else if (name) {
        StringCat (def, name);
    }
    if (orgnl != NULL) {
        if (chr != NULL) {
            StringCat (def, " ");
            StringCat (def, orgnl);
            StringCat (def, " chromosome ");
            StringCat(def, chr);
            StringCat (def, completeseq);
            ch = *def;
            *def = TO_UPPER (ch);
            return def;
        }
        StringCat (def, " ");
        StringCat (def, orgnl);
        StringCat (def, completegen);
        ch = *def;
        *def = TO_UPPER (ch);
        return def;
    }
    if (seg != NULL) {
        StringCat (def, " ");
        if (StringStr (seg, "DNA") == NULL &&
            StringStr (seg, "RNA") == NULL &&
            StringStr (seg, "segment") == NULL &&
            StringStr (seg, "Segment") == NULL) {
          StringCat (def, "segment ");
        }
        StringCat(def, seg);
        StringCat (def, completeseq);
        ch = *def;
        *def = TO_UPPER (ch);
        return def;
    }
    if (chr != NULL) {
        StringCat (def, " chromosome ");
        StringCat(def, chr);
        StringCat (def, completeseq);
        ch = *def;
        *def = TO_UPPER (ch);
        return def;
    }
    StringCat (def, completegen);
    ch = *def;
    *def = TO_UPPER (ch);
    return def;
}
static Boolean NotSpecialTaxName (CharPtr taxname)
{
  if (StringHasNoText (taxname)) return TRUE;
  if (StringICmp (taxname, "synthetic construct") == 0) return FALSE;
  if (StringICmp (taxname, "artificial sequence") == 0) return FALSE;
  if (StringStr (taxname, "vector") != NULL) return FALSE;
  if (StringStr (taxname, "Vector") != NULL) return FALSE;
  return TRUE;
}
static Boolean DoTpaPrefix (
  CharPtr title,
  CharPtr PNTR ttl,
  CharPtr PNTR pfx,
  Boolean is_tpa,
  Boolean tpa_exp,
  Boolean tpa_inf,
  Boolean is_tsa
)
{
  /* must be called with ttl and pfx pointing to stack variables */
  /* string literals declared here will persist and can be passed to calling function */
  *ttl = title;
  *pfx = NULL;
  if (title == NULL || *title == '\0') return FALSE;
  if (is_tsa) {
    if (StringNICmp (title, "TSA: ", 5) == 0) return FALSE;
    *pfx = "TSA: ";
    return TRUE;
  } else if (is_tpa) {
    if (tpa_exp) {
      if (StringNICmp (title, "TPA_exp: ", 9) == 0) return FALSE;
      *pfx = "TPA_exp: ";
      if (StringNICmp (title, "TPA: ", 5) == 0) {
        *ttl = title +  5;
      }
      return TRUE;
    } else if (tpa_inf) {
      if (StringNICmp (title, "TPA_inf: ", 9) == 0) return FALSE;
      *pfx = "TPA_inf: ";
      if (StringNICmp (title, "TPA: ", 5) == 0) {
        *ttl = title +  5;
      }
      return TRUE;
    } else {
      if (StringNICmp (title, "TPA: ", 5) == 0) return FALSE;
      *pfx = "TPA: ";
      return TRUE;
    }
  }
  return FALSE;
}

/*****************************************************************************
*
*   CreateDefLine(iip, bsp, buf, buflen, tech)
*       Finds or makes a FASTA format defline using Gather functions
*       buf should be very long if possible
*       function truncates if buf not long enough
*       a few deflines are longer than 255
*
*        ItemInfoPtr iip is used in flat file generator to keep entityId, itemId
*        and itemtype
*****************************************************************************/
NLM_EXTERN Boolean CreateDefLineExEx (ItemInfoPtr iip, BioseqPtr bsp, CharPtr buf, Uint4 buflen, Uint1 tech,
                                      CharPtr accession, CharPtr organism, Boolean ignoreTitle, Boolean extProtTitle)
{
    ValNodePtr vnp = NULL;
    CharPtr tmp = NULL, title = NULL, ttl = NULL, pfx = NULL;
    PdbBlockPtr pbp;
    PatentSeqIdPtr psip;
    PDBSeqIdPtr    pdbip;
    Uint4 diff, phase, i;
    Boolean doit;
    Int4 num_segs, num_gaps;
    static Char tbuf[80];
    static CharPtr htgs[2] = {
        "unordered", "ordered" };
    static CharPtr htg_phrase[3] = {
        "LOW-PASS SEQUENCE SAMPLING",
        "WORKING DRAFT SEQUENCE",
        "*** SEQUENCING IN PROGRESS ***" };
    Boolean htg_tech = FALSE, htgs_draft = FALSE, htgs_cancelled = FALSE,
            htgs_pooled_multiclone = FALSE, is_nc = FALSE, is_nm = FALSE,
            is_nr = FALSE, is_tpa = FALSE, tpa_exp = FALSE, tpa_inf = FALSE,
            is_tsa = FALSE;
    MolInfoPtr mip;
    GBBlockPtr gbp = NULL;
    EMBLBlockPtr ebp = NULL;
    ValNodePtr keywords = NULL;
    Boolean wgsmaster = FALSE;
    CharPtr suffix = NULL;
    SeqIdPtr sip;
    TextSeqIdPtr tsip;
    DbtagPtr general = NULL, dbt;
    ObjectIdPtr oip;
    ItemInfo ii;
    BioSourcePtr biop = NULL;
    OrgRefPtr orp;
    CharPtr taxname = NULL;
    SeqMgrDescContext dcontext;
    SeqMgrFeatContext fcontext;
    SeqFeatPtr sfp, src;
    Uint2 entityID;
    Uint1 genome;
    CharPtr orgnl = NULL;

    if ((bsp == NULL) || (buf == NULL) || buflen == 0) return FALSE;
    /* now using GetNextDescriptorUnindexed, so need to have called AssignIDsInEntityEx */
    if (bsp->idx.entityID == 0) {
      entityID = ObjMgrGetEntityIDForPointer (bsp);
      if (entityID != 0) {
        AssignIDsInEntityEx (entityID, 0, NULL, NULL);
      }
    }
    entityID = bsp->idx.entityID;
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
        switch (sip->choice) {
            case SEQID_OTHER :
                tsip = (TextSeqIdPtr) sip->data.ptrvalue;
                if (tsip != NULL && tsip->accession != NULL) {
                    if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
                        is_nc = TRUE;
                    } else if (StringNICmp (tsip->accession, "NM_", 3) == 0) {
                        is_nm = TRUE;
                    } else if (StringNICmp (tsip->accession, "NR_", 3) == 0) {
                        is_nr = TRUE;
                    }
                }
                break;
            case SEQID_TPG :
            case SEQID_TPE :
            case SEQID_TPD :
                is_tpa = TRUE;
                break;
            case SEQID_GENERAL :
                dbt = (DbtagPtr) sip->data.ptrvalue;
                if (dbt != NULL && (! IsSkippableDbtag (dbt))) {
                  general = dbt;
                }
                break;
            case SEQID_GENBANK :
            case SEQID_EMBL :
            case SEQID_DDBJ :
                tsip = (TextSeqIdPtr) sip->data.ptrvalue;
                if (tsip != NULL && tsip->accession != NULL) {
                    if (StringLen (tsip->accession) == 12) {
                        if (StringCmp (tsip->accession + 6, "000000") == 0) {
                            wgsmaster = TRUE;
                        }
                    }
                }
                break;
            case SEQID_GPIPE :
                tsip = (TextSeqIdPtr) sip->data.ptrvalue;
                break;
            default :
                break;
        }
    }
    buflen--;
    buf[buflen] = '\0';
    tbuf[0] = '\0';

    if (tech == 0) {
      vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
      if (vnp != NULL) {
        mip = (MolInfoPtr) vnp->data.ptrvalue;
        if (mip != NULL) {
          tech = mip->tech;
        }
      }
    }

    if (((tech >= MI_TECH_htgs_1) && (tech <= MI_TECH_htgs_3)) ||
        (tech == MI_TECH_htgs_0)) {
        htg_tech = TRUE;
    } else if (tech == MI_TECH_tsa) {
      is_tsa = TRUE;
    }
    if (iip == NULL && accession != NULL) {
        diff = LabelCopyExtra(buf, accession, buflen, "(", ") ");
        buflen -= diff;
        buf += diff;
    }
    diff = 0;
    if (htg_tech || is_tpa) {
        vnp=GatherDescrOnBioseq(iip, bsp, Seq_descr_genbank,TRUE);
        if (vnp != NULL) {
            gbp = (GBBlockPtr) vnp->data.ptrvalue;
            if (gbp != NULL) {
              keywords = gbp->keywords;
            }
        }
        vnp=GatherDescrOnBioseq(iip, bsp, Seq_descr_embl,TRUE);
        if (vnp != NULL) {
            ebp = (EMBLBlockPtr) vnp->data.ptrvalue;
            if (ebp != NULL) {
              keywords = ebp->keywords;
            }
        }
    }
    if (keywords != NULL) {
        for (vnp = keywords; vnp != NULL; vnp = vnp->next) {
            if (StringICmp ((CharPtr) vnp->data.ptrvalue, "HTGS_DRAFT") == 0) {
                htgs_draft = TRUE;
            } else if (StringICmp ((CharPtr) vnp->data.ptrvalue, "HTGS_CANCELLED") == 0) {
                htgs_cancelled = TRUE;
            } else if (StringICmp ((CharPtr) vnp->data.ptrvalue, "HTGS_POOLED_MULTICLONE") == 0 && htg_tech) {
                htgs_pooled_multiclone = TRUE;
            } else if (StringICmp ((CharPtr) vnp->data.ptrvalue, "TPA:experimental") == 0) {
                tpa_exp = TRUE;
            } else if (StringICmp ((CharPtr) vnp->data.ptrvalue, "TPA:inferential") == 0) {
                tpa_inf = TRUE;
            }
        }
    }
    if (! ignoreTitle)
          {
            vnp=GatherDescrOnBioseq(iip, bsp, Seq_descr_title,TRUE);
            if (vnp != NULL)
              title = StringSaveNoNull((CharPtr)vnp->data.ptrvalue);
              if (title != NULL) {
                TrimSpacesAroundString (title);
                TrimPunctuationFromEnd (title);
              }
          }
    if (tech == MI_TECH_htgs_0 || tech == MI_TECH_htgs_1 || tech == MI_TECH_htgs_2) {
        MemFree(title);  /* manufacture all HTG titles */
        title = NULL;
        if (iip != NULL) {
          iip->entityID = 0;
          iip->itemID = 0;
          iip->itemtype = 0;
        }
        if (title == NULL || *title == '\0') {
            title = UseOrgMods(bsp, NULL, tech, htgs_pooled_multiclone);
            organism = NULL;
        }
    } else if (tech == MI_TECH_est || tech == MI_TECH_sts || tech == MI_TECH_survey) {
        if (title == NULL || *title == '\0') {
            title = UseOrgMods(bsp, NULL, tech, FALSE);
            organism = NULL;
        }
    } else if (tech == MI_TECH_wgs) {
        if (title == NULL || *title == '\0') {
            if (! wgsmaster) {
                if (general != NULL) {
                    oip = general->tag;
                    if (oip != NULL) {
                        if (! StringHasNoText (oip->str)) {
                            suffix = oip->str;
                        }
                    }
                }
            }
            title = UseOrgMods(bsp, suffix, tech, FALSE);
            organism = NULL;
        }
    } else if (tech == MI_TECH_tsa) {
        if (title == NULL || *title == '\0') {
            if (general != NULL) {
                oip = general->tag;
                if (oip != NULL) {
                    if (! StringHasNoText (oip->str)) {
                        suffix = oip->str;
                    }
                }
            }
            title = UseOrgMods(bsp, suffix, tech, FALSE);
            organism = NULL;
        }
    } else if (is_nc && title == NULL) {
        /* manufacture complete chromosome titles if not already present */
        vnp = GatherDescrOnBioseq (&ii, bsp, Seq_descr_molinfo,TRUE);
        if (vnp != NULL) {
            mip = (MolInfoPtr) vnp->data.ptrvalue;
            if (mip != NULL &&
                (mip->biomol == MOLECULE_TYPE_GENOMIC || mip->biomol == MOLECULE_TYPE_OTHER_GENETIC_MATERIAL) /* && mip->completeness == 1 */) {
                title = MakeCompleteChromTitle (bsp, mip->biomol, mip->completeness);
                organism = NULL;
                if (iip != NULL) {
                    iip->entityID = ii.entityID;
                    iip->itemID = ii.itemID;
                    iip->itemtype = ii.itemtype;
                }
            }
        }
    } else if (is_nm && title == NULL) {
      title = FindNMDefLine (bsp);
      if (title != NULL && iip != NULL) {
        iip->entityID = 0;
        iip->itemID = 0;
        iip->itemtype = 0;
      }
    } else if (is_nr && title == NULL) {
      title = FindNRDefLine (bsp);
      if (title != NULL && iip != NULL) {
        iip->entityID = 0;
        iip->itemID = 0;
        iip->itemtype = 0;
      }
    }
/* some titles may have zero length */
    if (title != NULL && *title != '\0') {
        ttl = title;
        pfx = NULL;
        if (DoTpaPrefix (title, &ttl, &pfx, is_tpa, tpa_exp, tpa_inf, is_tsa)) {
            diff = LabelCopy (buf, pfx, buflen);
            buflen -= diff;
            buf += diff;
        }
        diff = LabelCopy (buf, ttl, buflen);
                                /* remove trailing blanks BUT NOT periods */
        tmp = buf + diff - 1;   /* point at last character */
        while (tmp >= buf && ((*tmp <= ' ') /* || (*tmp == '.') */)) {
            *tmp = '\0';
            tmp--;
            diff--;
        }
    } else if ((vnp = GatherDescrOnBioseq(iip, bsp, Seq_descr_pdb,TRUE)) != NULL) {
        pbp = (PdbBlockPtr)(vnp->data.ptrvalue);
        for (vnp = bsp->id; vnp != NULL; vnp = vnp->next) {
            if (vnp->choice == SEQID_PDB) {
                pdbip = (PDBSeqIdPtr)(vnp->data.ptrvalue);
                if (pdbip && pdbip->chain > 32) {
                    sprintf(tbuf, "Chain %c, ", pdbip->chain);
                    diff = LabelCopy(buf, tbuf, buflen);
                    buflen -= diff;
                    buf += diff;
                    break;
                }
            }
        }
        if (pbp && pbp->compound) {
            tmp = StringSave ((CharPtr)(pbp->compound->data.ptrvalue));
            TrimNonPeriodPunctuationFromEnd (tmp);
            diff = LabelCopy(buf, tmp, buflen);
            MemFree (tmp);
        }
    } else {
        for (vnp = bsp->id; vnp != NULL; vnp = vnp->next) {
            if (vnp->choice == SEQID_PATENT)
            {
                psip = (PatentSeqIdPtr)(vnp->data.ptrvalue);
                if (psip) {
                    sprintf(tbuf, "Sequence %d from Patent %s %s",
                    (int)psip->seqid, psip->cit->country, psip->cit->number);
                    diff = LabelCopy(buf, tbuf, buflen);
                    break;
                }
            }
        }
        if (vnp == NULL) {
            if (ISA_aa(bsp->mol)) {
                title = FindProtDefLine(bsp, extProtTitle);
                vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
                if (vnp != NULL && organism == NULL) {
                    biop = (BioSourcePtr) vnp->data.ptrvalue;
                    if (biop != NULL) {
                        orp = biop->org;
                        if (orp != NULL) {
                            taxname = orp->taxname;
                        }
                    }
                    if (taxname == NULL || NotSpecialTaxName (taxname)) {
                        if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
                            SeqMgrIndexFeatures (entityID, NULL);
                        }
                        sfp = SeqMgrGetCDSgivenProduct (bsp, NULL);
                        if (sfp != NULL) {
                            src = SeqMgrGetOverlappingSource (sfp->location, &fcontext);
                            if (src != NULL) {
                                biop = (BioSourcePtr) src->data.value.ptrvalue;
                                if (biop != NULL) {
                                    orp = biop->org;
                                    if (orp != NULL) {
                                        taxname = orp->taxname;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (title != NULL) {
                /*
                if (! StringHasNoText (taxname)) {
                    diff = LabelCopy(buf, taxname, buflen);
                    buflen -= diff;
                    buf += diff;
                    diff = LabelCopy(buf, " ", buflen);
                    buflen -= diff;
                    buf += diff;
                    diff = LabelCopy(buf, title, buflen);
                } else {
                    diff = LabelCopy(buf, title, buflen);
                }
                */
              ttl = title;
              pfx = NULL;
              if (DoTpaPrefix (title, &ttl, &pfx, is_tpa, tpa_exp, tpa_inf, is_tsa)) {
                    diff = LabelCopy (buf, pfx, buflen);
                    buflen -= diff;
                    buf += diff;
                }
                diff = LabelCopy (buf, ttl, buflen);
                if (organism == NULL && taxname != NULL) {
                    organism = taxname;
                    iip = NULL;
                }
            } else if (!htg_tech) {
                if (bsp->repr == Seq_repr_seg) {
                    title = SimpleSegSeqTitle (bsp);
                }
                if (title == NULL) {
                    title = UseOrgMods(bsp, NULL, tech, FALSE);
                }
                ttl = title;
                pfx = NULL;
                if (DoTpaPrefix (title, &ttl, &pfx, is_tpa, tpa_exp, tpa_inf, is_tsa)) {
                    diff = LabelCopy (buf, pfx, buflen);
                    buflen -= diff;
                    buf += diff;
                }
                if (ttl != NULL) {
                    diff = LabelCopy (buf, ttl, buflen);
                } else {
                    diff = LabelCopy (buf, "No definition line found", buflen);
                }
            }
        }
    }
    if (title != NULL) {
      TrimNonPeriodPunctuationFromEnd (title);
    }
    buflen -= diff;
    buf += diff;
    if (htg_tech) {
        if (tech == MI_TECH_htgs_0)
            phase = 0;
        else
            phase = (Uint4)(tech - MI_TECH_htgs_1 + 1);
        if (title == NULL|| *title == '\0') {
            title = UseOrgMods(bsp, NULL, tech, htgs_pooled_multiclone);
            organism = NULL;
            if (title != NULL) {
                ttl = title;
                pfx = NULL;
                if (DoTpaPrefix (title, &ttl, &pfx, is_tpa, tpa_exp, tpa_inf, is_tsa)) {
                    diff = LabelCopy (buf, pfx, buflen);
                    buflen -= diff;
                    buf += diff;
                }
                diff = LabelCopy (buf, ttl, buflen);
                buflen -= diff;
                buf += diff;
            }
        }
        if (phase == 3)
        {
            if (title) {
                if (title && StringStr(title, "complete sequence") == NULL) {
                    diff = LabelCopy(buf, ", complete sequence", buflen);
                    buflen -= diff;
                    buf += diff;
                }
            }
        } else {
            doit = FALSE;
            if (phase == 0) {
                if (StringStr(title, "LOW-PASS") == NULL) {
                    doit = TRUE;
                    i = 0;
                }
            } else {
                if (htgs_draft) {
                    if (StringStr(title, "WORKING DRAFT") == NULL) {
                        doit = TRUE;
                        i = 1;
                    }
                } else if (! htgs_cancelled) {
                    if (StringStr(title, "SEQUENCING IN") == NULL) {
                        doit = TRUE;
                        i = 2;
                    }
                }
            }
            if (doit)
            {
                if (diff != 0) {
                    diff = LabelCopy(buf, ", ", buflen);
                    buflen -= diff;
                    buf += diff;
                }
                diff = LabelCopy(buf, htg_phrase[i], buflen);
                buflen -= diff;
                buf += diff;
            }
            if ((phase != 0) && (bsp->repr == Seq_repr_delta)) {
                if (CountGapsInDeltaSeq(bsp,
                        &num_segs, &num_gaps, NULL, NULL, NULL, 0))
                {
                    if (num_gaps > 0) {
                        sprintf(tbuf, ", %ld %s pieces", (long)(num_gaps + 1), htgs[phase - 1]);
                    } else {
                        /*
                        sprintf(tbuf, ", %ld %s piece", (long)(num_gaps + 1), htgs[phase - 1]);
                        */
                    }
                    diff = LabelCopy(buf, tbuf, buflen);
                    buflen -= diff;
                    buf += diff;
                }
            }
            else if (phase != 0) {
                /*
                sprintf(tbuf, ", in %s pieces", htgs[phase-1]);
                diff = LabelCopy(buf, tbuf, buflen);
                buflen -= diff;
                buf += diff;
                */
            }
        }
    } else if (tech == MI_TECH_est || tech == MI_TECH_sts || tech == MI_TECH_survey || tech == MI_TECH_wgs) {
        if (title == NULL|| *title == '\0') {
            title = UseOrgMods(bsp, NULL, tech, FALSE);
            organism = NULL;
            if (title != NULL) {
                ttl = title;
                pfx = NULL;
                if (DoTpaPrefix (title, &ttl, &pfx, is_tpa, tpa_exp, tpa_inf, is_tsa)) {
                    diff = LabelCopy (buf, pfx, buflen);
                    buflen -= diff;
                    buf += diff;
                }
                diff = LabelCopy (buf, ttl, buflen);
                buflen -= diff;
                buf += diff;
            }
        }
        if (tech == MI_TECH_est) {
            if (title) {
                if (title && StringStr(title, "mRNA sequence") == NULL) {
                    diff = LabelCopy(buf, ", mRNA sequence", buflen);
                    buflen -= diff;
                    buf += diff;
                }
            }
        } else if (tech == MI_TECH_sts) {
            if (title) {
                if (title && StringStr(title, "sequence tagged site") == NULL) {
                    diff = LabelCopy(buf, ", sequence tagged site", buflen);
                    buflen -= diff;
                    buf += diff;
                }
            }
        } else if (tech == MI_TECH_survey) {
            if (title) {
                if (title && StringStr(title, "genomic survey sequence") == NULL) {
                    diff = LabelCopy(buf, ", genomic survey sequence", buflen);
                    buflen -= diff;
                    buf += diff;
                }
            }
        } else if (tech == MI_TECH_wgs) {
            if (title) {
                vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
                if (vnp != NULL) {
                    biop = (BioSourcePtr) vnp->data.ptrvalue;
                    if (biop != NULL) {
                        genome = biop->genome;
                        if (genome < kNumWGSOrganelles) {
                            orgnl = organelleForWGS [genome];
                        }
                    }
                }
                if (wgsmaster) {
                    if (title && StringStr (title, "whole genome shotgun sequencing project") == NULL) {
                        diff = LabelCopy(buf, " whole genome shotgun sequencing project", buflen);
                        buflen -= diff;
                        buf += diff;
                    }
                } else if (title && StringStr (title, "whole genome shotgun sequence") == NULL) {
                    if (orgnl != NULL && StringStr (title, orgnl) == NULL) {
                        diff = LabelCopy(buf, " ", buflen);
                        buflen -= diff;
                        buf += diff;
                        diff = LabelCopy(buf, orgnl, buflen);
                        buflen -= diff;
                        buf += diff;
                    }
                    diff = LabelCopy(buf, ", whole genome shotgun sequence", buflen);
                    buflen -= diff;
                    buf += diff;
                }
            }
        }
    }
    if (iip == NULL && organism != NULL) {
        doit = TRUE;
        if (title) {
            if (StringStr(title, organism) != NULL)
                doit = FALSE;
        }
        if (doit)
            LabelCopyExtra(buf, organism, buflen, " [", "]");
    }
        MemFree(title);
    return TRUE;
}
NLM_EXTERN Boolean CreateDefLineEx (ItemInfoPtr iip, BioseqPtr bsp, CharPtr buf, Uint4 buflen, Uint1 tech,
                                    CharPtr accession, CharPtr organism, Boolean ignoreTitle)
{
  return CreateDefLineExEx (iip, bsp, buf, buflen, tech, accession, organism, ignoreTitle, FALSE);
}
NLM_EXTERN Boolean CreateDefLine (ItemInfoPtr iip, BioseqPtr bsp, CharPtr buf, Uint4 buflen,
                                  Uint1 tech, CharPtr accession, CharPtr organism)
{
  return CreateDefLineExEx (iip, bsp, buf, buflen, tech, accession, organism, FALSE, FALSE);
}
/*****************************************************************************
*
*   FastaSeqPort(bsp, is_na, do_virtual)
*       opens a SeqPort for a fasta output of bsp
*
*****************************************************************************/
NLM_EXTERN SeqPortPtr FastaSeqPort(BioseqPtr bsp, Boolean is_na, Boolean do_virtual,
                                   Uint1 code)
{
    SeqPortPtr spp = NULL;
    if (bsp == NULL) return spp;
    spp = SeqPortNew(bsp, 0, -1, 0, code);
    if (do_virtual)
        SeqPortSet_do_virtual(spp, TRUE);
    SeqPortSeek(spp, 0, SEEK_SET);
    return spp;
}
/*****************************************************************************
*
*   FastaSeqPortEx(bsp, is_na, do_virtual, slp)
*       opens a SeqPort for a fasta output of bsp constrained to slp
*
*****************************************************************************/
NLM_EXTERN SeqPortPtr FastaSeqPortEx(BioseqPtr bsp, Boolean is_na, Boolean do_virtual,
                                     Uint1 code, SeqLocPtr slp)
{
    SeqPortPtr spp = NULL;
    if (bsp == NULL) return spp;
    if (slp == NULL) return FastaSeqPort (bsp, is_na, do_virtual, code);
    spp = SeqPortNew(bsp, SeqLocStart(slp), SeqLocStop(slp),
            SeqLocStrand(slp), code);
    if (do_virtual)
        SeqPortSet_do_virtual(spp, TRUE);
    SeqPortSeek(spp, 0, SEEK_SET);
    return spp;
}
/*****************************************************************************
*
*   FastaSeqLine(spp, buf, linelen)
*     an open seqport is passed in.
*     fills buf with linelen bases
*     assumes buf[linelen] = '\0'
*     returns FALSE when no more residues to print
*
*****************************************************************************/
NLM_EXTERN Boolean FastaSeqLine(SeqPortPtr spp, CharPtr buf, Int2 linelen, Boolean is_na)
{
    return FastaSeqLineEx(spp, buf, linelen, is_na, FALSE);
}
NLM_EXTERN Boolean FastaSeqLineEx(SeqPortPtr spp, CharPtr buf, Int2 linelen, Boolean is_na, Boolean
do_virtual)
{
    Int2 ctr = 0;
    Uint1 residue;
    Int4 pos;
    Char idbuf[40];
    if ((spp == NULL) || (buf == NULL)) return FALSE;
    while ((residue = SeqPortGetResidue(spp)) != SEQPORT_EOF)
    {
        if (! IS_residue(residue))
        {
            if (residue == INVALID_RESIDUE)
            {
                if (is_na)
                    residue = 'N';
                else
                    residue = 'X';
                FastaId(spp->bsp, idbuf, 39);
                pos = SeqPortTell(spp);
                ErrPostEx(SEV_ERROR,0,0, "ToFastA: Invalid residue at position %ld in %s",
                    (long) pos, idbuf);
            }
            else
            {
                if (residue == SEQPORT_VIRT)  /* gap */
                {
                    if (ctr)  /* got some residues already */
                    {
                        buf[ctr] = '\0';
                        SeqPortSeek(spp, -1, SEEK_CUR);  /* back up one */
                                     /* can only seek to a real residue, so go past it */
                                                residue = SeqPortGetResidue(spp);
                                                if (residue == SEQPORT_VIRT)
                                                   SeqPortSeek(spp, -1, SEEK_CUR);
                        return TRUE;
                    }
                    else if (! do_virtual)       /* first one */
                    {
                        buf[ctr] = '-';
                        buf[ctr + 1] = '\0';
                        return TRUE;
                    }
                }
                residue = '\0';
            }
        }
        if (residue != '\0')
        {
            buf[ctr] = residue;
            ctr++;
            if (ctr == linelen)
            {
                buf[ctr] = '\0';
                return TRUE;
            }
        }
    }
    buf[ctr] = '\0';
    if (ctr)
        return TRUE;
    else
        return FALSE;
}
/*****************************************************************************
*
*   NC_Cleanup (entityID, ptr)
*     internal function for genome RefSeq processing
*
*****************************************************************************/
static Boolean RemoveAllTitles (GatherObjectPtr gop)
{
  ObjValNodePtr  ovp;
  SeqDescrPtr    sdp;
  if (gop == NULL ||
      gop->itemtype != OBJ_SEQDESC ||
      gop->subtype != Seq_descr_title) return TRUE;
  sdp = (SeqDescrPtr) gop->dataptr;
  if (sdp == NULL || sdp->extended == 0) return TRUE;
  ovp = (ObjValNodePtr) sdp;
  ovp->idx.deleteme = TRUE;
  return TRUE;
}
static Boolean AddNcTitles (GatherObjectPtr gop)
{
  BioseqPtr     bsp;
  Char          buf [512];
  Boolean       is_nc;
  /*
  MolInfoPtr    mip;
  SeqDescrPtr   sdp;
  */
  SeqIdPtr      sip;
  CharPtr       str;
  TextSeqIdPtr  tsip;
  if (gop == NULL ||
      gop->itemtype != OBJ_BIOSEQ) return TRUE;
  bsp = (BioseqPtr) gop->dataptr;
  if (bsp == NULL) return TRUE;
  is_nc = FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
          is_nc = TRUE;
        }
      }
    }
  }
  if (! is_nc) return TRUE;
  if (NewCreateDefLineBuf (NULL, bsp, buf, sizeof (buf), FALSE, FALSE)) {
    if (! StringHasNoText (buf)) {
      str = StringSaveNoNull (buf);
      if (str != NULL) {
        SeqDescrAddPointer (&(bsp->descr), Seq_descr_title, (Pointer) str);
      }
    }
  }
  /*
  for (sdp = bsp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_molinfo) {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL &&
          mip->biomol == MOLECULE_TYPE_GENOMIC &&
          mip->completeness == 1) {
        mip->completeness = 0;
      }
    }
  }
  */
  return TRUE;
}
static void ClearKeywordsProc (SeqDescrPtr sdp, Pointer userdata)
{
  GBBlockPtr     gbp;
  ObjValNodePtr  ovn;
  if (sdp == NULL || sdp->choice != Seq_descr_genbank) return;
  gbp = (GBBlockPtr) sdp->data.ptrvalue;
  if (gbp == NULL) return;
  gbp->keywords = ValNodeFreeData (gbp->keywords);
  if (gbp->extra_accessions == NULL && gbp->source == NULL &&
      gbp->keywords == NULL && gbp->origin == NULL &&
      gbp->date == NULL && gbp->entry_date == NULL &&
      gbp->div == NULL && gbp->taxonomy == NULL) {
  }
  if (sdp->extended == 0) return;
  ovn = (ObjValNodePtr) sdp;
  ovn->idx.deleteme = TRUE;
}
NLM_EXTERN void ClearGenBankKeywords (Uint2 entityID, Pointer ptr)
{
  SeqEntryPtr  sep;
  if (entityID == 0) {
    entityID = ObjMgrGetEntityIDForPointer (ptr);
  }
  if (entityID == 0) return;
  sep = GetTopSeqEntryForEntityID (entityID);
  VisitDescriptorsInSep (sep, NULL, ClearKeywordsProc);
  DeleteMarkedObjects (entityID, 0, NULL);
}
NLM_EXTERN void NC_Cleanup (Uint2 entityID, Pointer ptr)
{
  Boolean      objMgrFilt [OBJ_MAX];
  SeqEntryPtr  sep;
  if (entityID == 0) {
    entityID = ObjMgrGetEntityIDForPointer (ptr);
  }
  if (entityID == 0) return;
  AssignIDsInEntity (entityID, 0, NULL);
  MemSet ((Pointer) objMgrFilt, FALSE, sizeof (objMgrFilt));
  objMgrFilt [OBJ_SEQDESC] = TRUE;
  GatherObjectsInEntity (entityID, 0, NULL, RemoveAllTitles, NULL, objMgrFilt);
  sep = GetTopSeqEntryForEntityID (entityID);
  VisitDescriptorsInSep (sep, NULL, ClearKeywordsProc);
  DeleteMarkedObjects (entityID, 0, NULL);
  MemSet ((Pointer) objMgrFilt, FALSE, sizeof (objMgrFilt));
  objMgrFilt [OBJ_BIOSEQ] = TRUE;
  GatherObjectsInEntity (entityID, 0, NULL, AddNcTitles, NULL, objMgrFilt);
}

NLM_EXTERN void InstantiateNCTitle (Uint2 entityID, Pointer ptr)
{
  Boolean      objMgrFilt [OBJ_MAX];

  if (entityID == 0) {
    entityID = ObjMgrGetEntityIDForPointer (ptr);
  }
  if (entityID == 0) return;

  AssignIDsInEntity (entityID, 0, NULL);
  MemSet ((Pointer) objMgrFilt, FALSE, sizeof (objMgrFilt));
  objMgrFilt [OBJ_BIOSEQ] = TRUE;
  GatherObjectsInEntity (entityID, 0, NULL, AddNcTitles, NULL, objMgrFilt);
}

static Boolean AddNmTitles (GatherObjectPtr gop)
{
  BioseqPtr     bsp;
  Char          buf [512];
  Boolean       is_nm;
  SeqIdPtr      sip;
  CharPtr       str;
  TextSeqIdPtr  tsip;
  if (gop == NULL ||
      gop->itemtype != OBJ_BIOSEQ) return TRUE;
  bsp = (BioseqPtr) gop->dataptr;
  if (bsp == NULL) return TRUE;
  is_nm = FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        if (StringNICmp (tsip->accession, "NM_", 3) == 0) {
          is_nm = TRUE;
        } else if (StringNICmp (tsip->accession, "XM_", 3) == 0) {
          is_nm = TRUE;
        }
      }
    }
  }
  if (! is_nm) return TRUE;
  if (NewCreateDefLineBuf (NULL, bsp, buf, sizeof (buf), FALSE, FALSE)) {
    if (! StringHasNoText (buf)) {
      str = StringSaveNoNull (buf);
      if (str != NULL) {
        SeqDescrAddPointer (&(bsp->descr), Seq_descr_title, (Pointer) str);
      }
    }
  }
  return TRUE;
}

NLM_EXTERN void InstantiateNMTitles (Uint2 entityID, Pointer ptr)
{
  Boolean      objMgrFilt [OBJ_MAX];

  if (entityID == 0) {
    entityID = ObjMgrGetEntityIDForPointer (ptr);
  }
  if (entityID == 0) return;

  AssignIDsInEntity (entityID, 0, NULL);
  MemSet ((Pointer) objMgrFilt, FALSE, sizeof (objMgrFilt));
  objMgrFilt [OBJ_BIOSEQ] = TRUE;
  GatherObjectsInEntity (entityID, 0, NULL, AddNmTitles, NULL, objMgrFilt);
}

static void ClearProtTitlesProc (BioseqPtr bsp, Pointer userdata)
{
  ObjValNodePtr  ovp;
  SeqDescrPtr    sdp;
  SeqIdPtr       sip;
  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) return;
  }
  for (sdp = bsp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_title) {
      if (sdp->extended != 0) {
        ovp = (ObjValNodePtr) sdp;
        ovp->idx.deleteme = TRUE;
      }
    }
  }
}
static void ClearProtTitlesNPS (BioseqSetPtr bssp, Pointer userdata)
{
  if (bssp->_class != BioseqseqSet_class_nuc_prot) return;
  VisitBioseqsInSet (bssp, NULL, ClearProtTitlesProc);
}
NLM_EXTERN void ClearProteinTitlesInNucProts (Uint2 entityID, Pointer ptr)
{
  SeqEntryPtr  sep;
  if (entityID == 0) {
    entityID = ObjMgrGetEntityIDForPointer (ptr);
  }
  if (entityID == 0) return;
  sep = GetTopSeqEntryForEntityID (entityID);
  VisitSetsInSep (sep, NULL, ClearProtTitlesNPS);
  DeleteMarkedObjects (entityID, 0, NULL);
}


static void AddProtTitles (BioseqPtr bsp, Pointer userdata)
{
  Char         buf [512];
  SeqDescrPtr  sdp;
  SeqIdPtr     sip;
  CharPtr      str;
  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_PIR ||
        sip->choice == SEQID_SWISSPROT ||
        sip->choice == SEQID_PATENT ||
        sip->choice == SEQID_PRF ||
        sip->choice == SEQID_PDB) return;
  }
  for (sdp = bsp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_title) return;
  }
  if (NewCreateDefLineBuf (NULL, bsp, buf, sizeof (buf), FALSE, FALSE)) {
    if (! StringHasNoText (buf)) {
      str = StringSaveNoNull (buf);
      if (str != NULL) {
        SeqDescrAddPointer (&(bsp->descr), Seq_descr_title, (Pointer) str);
      }
    }
  }
}
NLM_EXTERN void InstantiateProteinTitles (Uint2 entityID, Pointer ptr)
{
  SeqEntryPtr  sep;
  if (entityID == 0) {
    entityID = ObjMgrGetEntityIDForPointer (ptr);
  }
  if (entityID == 0) return;
  AssignIDsInEntity (entityID, 0, NULL);
  sep = GetTopSeqEntryForEntityID (entityID);
  VisitBioseqsInSep (sep, NULL, AddProtTitles);
}


NLM_EXTERN void UpdateProteinTitle (BioseqPtr bsp)
{
  Char         buf [512];
  SeqDescrPtr  sdp;
  ObjValNodePtr ovp;
  SeqIdPtr     sip;
  CharPtr      str;

  if (bsp == NULL || !ISA_aa (bsp->mol)) {
    return;
  }

  /* we don't create protein titles for these IDs */
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_PIR ||
        sip->choice == SEQID_SWISSPROT ||
        sip->choice == SEQID_PATENT ||
        sip->choice == SEQID_PRF ||
        sip->choice == SEQID_PDB) return;
  }

  sdp = BioseqGetSeqDescr (bsp, Seq_descr_title, NULL);
  if (sdp == NULL) {
    /* we only update a title if it already exists */
    return;
  }
  if (sdp->extended) {
    ovp = (ObjValNodePtr) sdp;
    ovp->idx.deleteme = TRUE;
    DeleteMarkedObjects (bsp->idx.entityID, OBJ_BIOSEQ, bsp);
  }

  if (NewCreateDefLineBuf (NULL, bsp, buf, sizeof (buf), FALSE, FALSE)) {
    if (! StringHasNoText (buf)) {
      str = StringSaveNoNull (buf);
      if (str != NULL) {
        SeqDescrAddPointer (&(bsp->descr), Seq_descr_title, (Pointer) str);
      }
    }
  }
}


/* NEW DEFLINE GENERATOR */

typedef struct deflinestruct {
  /* instance variables */
  ItemInfoPtr  m_iip;
  BioseqPtr    m_bioseq;

  /* ignore existing title is forced for certain types */
  Boolean  m_reconstruct;
  Boolean  m_allprotnames;

  /* seq-inst fields */
  Boolean m_is_na;
  Boolean m_is_aa;

  Boolean m_is_seg;
  Boolean m_is_delta;

  /* seq-id fields */
  Boolean m_is_nc;
  Boolean m_is_nm;
  Boolean m_is_nr;
  Boolean m_is_patent;
  Boolean m_is_pdb;
  Boolean m_third_party;
  Boolean m_wgs_master;

  CharPtr m_general_str;
  CharPtr m_patent_country;
  CharPtr m_patent_number;
  int     m_patent_sequence;

  int     m_pdb_chain;

  /* molinfo fields */
  Uint1   m_mi_biomol;
  Uint1   m_mi_tech;
  Uint1   m_mi_completeness;

  Boolean m_htg_tech;
  Boolean m_htgs_unfinished;
  Boolean m_is_tsa;
  Boolean m_is_wgs;
  Boolean m_is_est_sts_gss;

  Boolean m_use_biosrc;

  /* genbank or embl block keyword fields */
  Boolean m_htgs_cancelled;
  Boolean m_htgs_draft;
  Boolean m_htgs_pooled;
  Boolean m_tpa_exp;
  Boolean m_tpa_inf;
  Boolean m_tpa_reasm;

  /* pdb block fields */
  CharPtr m_pdb_compound;

  /* biosource fields */
  CharPtr m_taxname;
  int     m_genome;

  /* subsource fields */
  CharPtr m_chromosome;
  CharPtr m_clone;
  Boolean m_has_clone;
  CharPtr m_map;
  CharPtr m_plasmid;
  CharPtr m_segment;

  /* orgmod fields */
  CharPtr m_isolate;
  CharPtr m_strain;

  /* user object fields */
  Boolean m_is_unverified;

  /* exception fields */
  TextFsaPtr m_low_quality_fsa;
} DefLineData, PNTR DefLinePtr;

static Boolean x_CDShasLowQualityException (
  DefLinePtr dlp,
  SeqFeatPtr sfp
)

{
  Char        ch;
  TextFsaPtr  fsa;
  ValNodePtr  matches;
  CharPtr     ptr;
  Int4        state;

  if (dlp == NULL || sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return FALSE;

  if (! sfp->excpt) return FALSE;
  if (StringHasNoText (sfp->except_text)) return FALSE;

  fsa = dlp->m_low_quality_fsa;
  if (fsa == NULL) return FALSE;

  state = 0;
  matches = NULL;
  for (ptr = sfp->except_text, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    state = TextFsaNext (fsa, state, ch, &matches);
    if (matches != NULL) {
      return TRUE;
    }
  }

  return FALSE;
}

/* set instance variables from Seq-inst, Seq-ids, MolInfo, etc., but not BioSource */
static void x_SetFlags (
  DefLinePtr dlp
)

{
  BioseqPtr       bsp;
  IdPatPtr        cit;
  ValNodePtr      compound;
  DbtagPtr        dbt;
  EMBLBlockPtr    ebp;
  GBBlockPtr      gbp;
  DbtagPtr        general;
  ValNodePtr      keywords;
  MolInfoPtr      mip;
  ObjectIdPtr     oip;
  PdbBlockPtr     pbp;
  PDBSeqIdPtr     pdbip;
  PatentSeqIdPtr  psip;
  SeqDescrPtr     sdp;
  SeqIdPtr        sip;
  CharPtr         str;
  TextSeqIdPtr    tsip;
  UserObjectPtr   uop;
  ValNodePtr      vnp;

  if (dlp == NULL) return;

  bsp = dlp->m_bioseq;
  if (bsp == NULL) return;

  dlp->m_is_na = (Boolean) ISA_na (bsp->mol);
  dlp->m_is_aa = (Boolean) ISA_aa (bsp->mol);

  dlp->m_is_seg = (Boolean) (bsp->repr == Seq_repr_seg);
  dlp->m_is_delta = (Boolean) (bsp->repr == Seq_repr_delta);

  /* process Seq-ids */
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
        case SEQID_OTHER :
          tsip = (TextSeqIdPtr) sip->data.ptrvalue;
          if (tsip != NULL && tsip->accession != NULL) {
            if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
              dlp->m_is_nc = TRUE;
            } else if (StringNICmp (tsip->accession, "NM_", 3) == 0) {
              dlp->m_is_nm = TRUE;
            } else if (StringNICmp (tsip->accession, "XM_", 3) == 0) {
              dlp->m_is_nm = TRUE;
            } else if (StringNICmp (tsip->accession, "NR_", 3) == 0) {
              dlp->m_is_nr = TRUE;
            }
          }
          break;
        case SEQID_GENBANK :
        case SEQID_EMBL :
        case SEQID_DDBJ :
          tsip = (TextSeqIdPtr) sip->data.ptrvalue;
          if (tsip != NULL && tsip->accession != NULL) {
            if (StringLen (tsip->accession) == 12) {
              if (StringCmp (tsip->accession + 6, "000000") == 0) {
                dlp->m_wgs_master = TRUE;
              }
            }
          }
          break;
        case SEQID_GENERAL :
          dbt = (DbtagPtr) sip->data.ptrvalue;
          if (dbt != NULL && (! IsSkippableDbtag (dbt))) {
            general = dbt;
            if (general != NULL) {
              oip = general->tag;
              if (oip != NULL) {
                if (! StringHasNoText (oip->str)) {
                  dlp->m_general_str = oip->str;
                }
              }
            }
          }
          break;
        case SEQID_TPG :
        case SEQID_TPE :
        case SEQID_TPD :
          dlp->m_third_party = TRUE;
          break;
        case SEQID_PDB :
          dlp->m_is_pdb = TRUE;
          pdbip = (PDBSeqIdPtr) sip->data.ptrvalue;
          if (pdbip && pdbip->chain > 32) {
            dlp->m_pdb_chain = pdbip->chain;
          }
          break;
        case SEQID_PATENT :
          dlp->m_is_patent = TRUE;
          psip = (PatentSeqIdPtr) sip->data.ptrvalue;
          if (psip != NULL) {
            dlp->m_patent_sequence = (int) psip->seqid;
            cit = psip->cit;
            if (cit != NULL) {
              dlp->m_patent_country = cit->country;
              if (StringDoesHaveText (cit->number)) {
                dlp->m_patent_number = cit->number;
              } else if (StringDoesHaveText (cit->app_number)) {
                dlp->m_patent_number = cit->app_number;
              }
            }
          }
          break;
        case SEQID_GPIPE :
          break;
        default :
          break;
    }
  }

  /* process MolInfo tech */
  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
  if (sdp != NULL && sdp->choice == Seq_descr_molinfo) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      dlp->m_mi_biomol = mip->biomol;
      dlp->m_mi_tech = mip->tech;
      dlp->m_mi_completeness = mip->completeness;
      switch (dlp->m_mi_tech) {
          case MI_TECH_htgs_0 :
          case MI_TECH_htgs_1 :
          case MI_TECH_htgs_2 :
            dlp->m_htgs_unfinished = TRUE;
            /* manufacture all titles for unfinished HTG sequences */
            dlp->m_reconstruct = TRUE;
            /* fall through */
          case MI_TECH_htgs_3 :
            dlp->m_htg_tech = TRUE;
            dlp->m_use_biosrc = TRUE;
            break;
          case MI_TECH_est :
          case MI_TECH_sts :
          case MI_TECH_survey :
            dlp->m_is_est_sts_gss = TRUE;
            dlp->m_use_biosrc = TRUE;
            break;
          case MI_TECH_wgs :
            dlp->m_is_wgs = TRUE;
            dlp->m_use_biosrc = TRUE;
            break;
          case MI_TECH_tsa :
            dlp->m_is_tsa = TRUE;
            dlp->m_use_biosrc = TRUE;
            break;
          default :
            break;
      }
    }
  }

  /* process Unverified user object */
  for (sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_user, NULL);
       sdp != NULL;
       sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_user, sdp)) {
    if (sdp->choice != Seq_descr_user) continue;
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop == NULL) continue;
    oip = uop->type;
    if (oip == NULL) continue;
    if (StringICmp (oip->str, "Unverified") != 0) continue;
    dlp->m_is_unverified = TRUE;
  }

  if (dlp->m_htg_tech || dlp->m_third_party) {
    /* process keywords */
    keywords = NULL;

    sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_genbank, NULL);
    if (sdp != NULL && sdp->choice == Seq_descr_genbank) {
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      if (gbp != NULL) {
        keywords = gbp->keywords;
      }
    }
    if (keywords == NULL) {
      sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_embl, NULL);
      if (sdp != NULL && sdp->choice == Seq_descr_embl) {
        ebp = (EMBLBlockPtr) sdp->data.ptrvalue;
        if (ebp != NULL) {
          keywords = ebp->keywords;
        }
      }
    }
    if (keywords != NULL) {
      for (vnp = keywords; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (StringICmp (str, "HTGS_DRAFT") == 0) {
          dlp->m_htgs_draft = TRUE;
        } else if (StringICmp (str, "HTGS_CANCELLED") == 0) {
          dlp->m_htgs_cancelled = TRUE;
        } else if (StringICmp (str, "HTGS_POOLED_MULTICLONE") == 0) {
          dlp->m_htgs_pooled = TRUE;
        } else if (StringICmp (str, "TPA:experimental") == 0) {
          dlp->m_tpa_exp = TRUE;
        } else if (StringICmp (str, "TPA:inferential") == 0) {
          dlp->m_tpa_inf = TRUE;
        } else if (StringICmp (str, "TPA:reassembly") == 0) {
          dlp->m_tpa_reasm = TRUE;
        }
      }
    }
  }

  if (dlp->m_is_pdb) {

    /* process PDB block */
    sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_pdb, NULL);
    if (sdp != NULL && sdp->choice == Seq_descr_pdb) {
      pbp = (PdbBlockPtr) sdp->data.ptrvalue;
      if (pbp != NULL) {
        compound = pbp->compound;
        if (compound != NULL) {
          dlp->m_pdb_compound = (CharPtr) compound->data.ptrvalue;
        }
      }
    }
  }
}

/* set instance variables from BioSource */
static void x_SetSrcClone (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioSourcePtr  biop;
  DefLinePtr    dlp;
  SubSourcePtr  ssp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC) return;
  dlp = (DefLinePtr) userdata;
  if (dlp == NULL) return;

  biop = (BioSourcePtr) sfp->data.value.ptrvalue;
  if (biop == NULL) return;

  /* look for clones on source features */
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (StringHasNoText (ssp->name)) continue;
    if (ssp->subtype != SUBSRC_clone) continue;
    dlp->m_has_clone = TRUE;
  }
}

static void x_SetBioSrc (
  DefLinePtr dlp
)

{
  BioSourcePtr  biop;
  BioseqPtr     bsp;
  OrgModPtr     omp;
  OrgNamePtr    onp;
  OrgRefPtr     orp;
  SeqDescrPtr   sdp;
  SubSourcePtr  ssp;

  if (dlp == NULL) return;

  bsp = dlp->m_bioseq;
  if (bsp == NULL) return;

  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_source, NULL);
  if (sdp != NULL && sdp->choice == Seq_descr_source) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL) {
      orp = biop->org;
      if (orp != NULL) {
        if (StringDoesHaveText (orp->taxname)) {
          dlp->m_taxname = orp->taxname;
        }
      }
      dlp->m_genome = biop->genome;

      /* process SubSource */
      for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
        if (StringHasNoText (ssp->name)) continue;
        switch (ssp->subtype) {
            case SUBSRC_chromosome :
              dlp->m_chromosome = ssp->name;
              break;
            case SUBSRC_clone :
              dlp->m_clone = ssp->name;
              dlp->m_has_clone = TRUE;
              break;
            case SUBSRC_map :
              dlp->m_map = ssp->name;
              break;
            case SUBSRC_plasmid_name :
              dlp->m_plasmid = ssp->name;
              break;
            case SUBSRC_segment :
              dlp->m_segment = ssp->name;
              break;
            default :
              break;
        }
      }

      /* process OrgMod */
      if (orp != NULL) {
        onp = orp->orgname;
        if (onp != NULL) {
          for (omp = onp->mod; omp != NULL; omp = omp->next) {
            if (StringHasNoText (omp->subname)) continue;
            switch (omp->subtype) {
                case ORGMOD_strain :
                  if (StringHasNoText (dlp->m_strain)) {
                    dlp->m_strain = omp->subname;
                  }
                  break;
                case ORGMOD_isolate :
                  if (StringHasNoText (dlp->m_isolate)) {
                    dlp->m_isolate = omp->subname;
                  }
                  break;
                default :
                  break;
            }
          }
        }
      }
    }
  }

  if (dlp->m_has_clone) return;

  VisitFeaturesOnBsp (bsp, (Pointer) dlp, x_SetSrcClone);
}

static CharPtr x_TrimFirstNCharacters (
  CharPtr str,
  Int2 count
)

{
  Uchar    ch;      /* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && count > 0) {
      count--;
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
  return str;
}

static CharPtr x_TrimPunctuationFromEnd (
  CharPtr str
)

{
  Uchar    ch;      /* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == ' ' || ch == ';' || ch == ',' || ch == '~' || ch == '.') {
        if (dst == NULL) {
          dst = ptr;
        }
      } else  {
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

static CharPtr x_TrimMostPunctFromEnd (
  CharPtr str
)

{
  Uchar    ch;      /* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == ' ' || ch == ';' || ch == ',' || ch == '~') {
        if (dst == NULL) {
          dst = ptr;
        }
      } else  {
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

static CharPtr x_CatenateValNodeStrings (
  ValNodePtr list
)

{
  size_t      len;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     tmp;
  ValNodePtr  vnp;


  ptr = NULL;
  if (list != NULL) {
    vnp = list;
    len = 0;
    while (vnp != NULL) {
      if (vnp->data.ptrvalue != NULL) {
        len += StringLen ((CharPtr) vnp->data.ptrvalue) + 1;
      }
      vnp = vnp->next;
    }
    if (len > 0) {
      ptr = MemNew (sizeof (Char) * (len + 2));
      if (ptr != NULL) {
        vnp = list;
        tmp = ptr;
        while (vnp != NULL) {
          str = (CharPtr) vnp->data.ptrvalue;
          /* do not use StringDoesHaveText because generalID must be prefixed by space */
          if (str != NULL) {
            tmp = StringMove (tmp, str);
          }
          vnp = vnp->next;
        }
      }
    }
  }
  return ptr;
}

static CharPtr x_DescribeClones (
  DefLinePtr dlp
)

{
  Char     buf [40];
  Char     ch;
  Int4     count;
  size_t   len;
  CharPtr  result = NULL;
  CharPtr  str;

  if (dlp == NULL) return NULL;

  if (dlp->m_htgs_unfinished && dlp->m_htgs_pooled && dlp->m_has_clone) {
    result = StringSave (", pooled multiple clones");
    return result;
  }

  str = dlp->m_clone;
  if (StringHasNoText (str)) return NULL;

  count = 1;
  ch = *str;
  while (ch != '\0') {
    if (ch == ';') {
      count++;
    }
    str++;
    ch = *str;
  }

  if (count > 3) {
    sprintf (buf, ", %d clones", (int) count);
    result = StringSave (buf);
  } else {
    len = StringLen (dlp->m_clone) + 20;
    result = (CharPtr) MemNew (sizeof (Char) * len);
    if (result != NULL) {
      StringCat (result, " clone ");
      StringCat (result, dlp->m_clone);
    }
  }

  return result;
}

static Boolean x_EndsWithStrain (
  DefLinePtr dlp
)

{
  Char     ch;
  size_t   len;
  CharPtr  nxt;
  CharPtr  ptr;

  if (dlp == NULL) return FALSE;

  len = StringLen (dlp->m_strain);
  if (len >= StringLen (dlp->m_taxname)) return FALSE;

  ptr = StringChr (dlp->m_taxname, ' ');
  if (ptr == NULL) return FALSE;
  ptr++;
  ptr = StringChr (ptr, ' ');
  if (ptr == NULL) return FALSE;
  ptr++;

  ptr = StringISearch (dlp->m_taxname, dlp->m_strain);
  if (ptr == NULL) return FALSE;

  nxt = StringISearch (ptr + 1, dlp->m_strain);
  while (nxt != NULL) {
    ptr = nxt;
    nxt = StringISearch (ptr + 1, dlp->m_strain);
  }

  ptr += len;
  if (! StringHasNoText (ptr)) {
    if (StringCmp (ptr, "'") == 0) {
      ptr -= len + 1;
      if (*ptr == '\'') return TRUE;
    }
    return FALSE;
  }
  ptr -= len + 1;
  ch = *ptr;
  /*
  if (ch == ' ' || ch == '-' || ch == '_' || ch == ':' ||
      ch == ';' || ch == '.' || ch == '/') {
     return TRUE;
  }
  */
  if (ispunct (ch) || isspace (ch)) {
    return TRUE;
  }

  return FALSE;
}

static CharPtr x_TitleFromBioSrc (
  DefLinePtr dlp
)

{
  Char        ch;
  CharPtr     result = NULL, cln, stn, ptr;
  ValNodePtr  strings = NULL;

  if (dlp == NULL) return NULL;

  ValNodeCopyStr (&strings, 0, dlp->m_taxname);

  if (StringDoesHaveText (dlp->m_strain)) {
    if (! x_EndsWithStrain (dlp)) {
      ValNodeCopyStr (&strings, 0, " strain ");
      stn = StringSave (dlp->m_strain);
      ptr = StringChr (stn, ';');
      if (ptr != NULL) {
        *ptr = '\0';
      }
      ValNodeCopyStr (&strings, 0, stn);
      MemFree (stn);
    }
  }

  if (StringDoesHaveText (dlp->m_chromosome)) {
    ValNodeCopyStr (&strings, 0, " chromosome ");
    ValNodeCopyStr (&strings, 0, dlp->m_chromosome);
  }

  cln = x_DescribeClones (dlp);
  if (StringDoesHaveText (cln)) {
    ValNodeCopyStr (&strings, 0, cln);
  }
  MemFree (cln); 

  if (StringDoesHaveText (dlp->m_map)) {
    ValNodeCopyStr (&strings, 0, " map ");
    ValNodeCopyStr (&strings, 0, dlp->m_map);
  }

  if (StringDoesHaveText (dlp->m_plasmid)) {
    if (dlp->m_is_wgs) {
      ValNodeCopyStr (&strings, 0, " plasmid ");
      ValNodeCopyStr (&strings, 0, dlp->m_plasmid);
    }
  }

  result = x_CatenateValNodeStrings (strings);
  ValNodeFreeData (strings);
  if (result == NULL) return NULL;

  ch = result [0];
  if (IS_LOWER (ch)) {
    result [0] = TO_UPPER (ch);
  }

  return result;
}

static CharPtr x_OrganelleName (
  DefLinePtr dlp,
  Boolean has_plasmid,
  Boolean virus_or_phage,
  Boolean wgs_suffix
)

{
  CharPtr  result = NULL;

  if (dlp == NULL) return NULL;

  switch (dlp->m_genome) {
      case GENOME_chloroplast :
        result = "chloroplast";
        break;
      case GENOME_chromoplast :
        result = "chromoplast";
        break;
      case GENOME_kinetoplast :
        result = "kinetoplast";
        break;
      case GENOME_mitochondrion :
        {
          if (has_plasmid || wgs_suffix) {
            result = "mitochondrial";
          } else {
            result = "mitochondrion";
          }
          break;
        }
      case GENOME_plastid :
        result = "plastid";
        break;
      case GENOME_macronuclear :
        {
          if (! wgs_suffix) {
            result = "macronuclear";
          }
          break;
        }
      case GENOME_extrachrom :
        {
          if (! wgs_suffix) {
            result = "extrachromosomal";
          }
          break;
        }
      case GENOME_plasmid :
        {
          if (! wgs_suffix) {
            result = "plasmid";
          }
          break;
        }
        /* transposon and insertion-seq are obsolete */
      case GENOME_cyanelle :
        result = "cyanelle";
        break;
      case GENOME_proviral :
        {
          if (! virus_or_phage) {
            if (has_plasmid || wgs_suffix) {
              result = "proviral";
            } else {
              result = "provirus";
            }
          }
          break;
        }
      case GENOME_virion :
        {
          if (! virus_or_phage) {
            result = "virus";
          }
          break;
        }
      case GENOME_nucleomorph :
        {
          if (! wgs_suffix) {
            result = "nucleomorph";
          }
          break;
        }
      case GENOME_apicoplast :
        result = "apicoplast";
        break;
      case GENOME_leucoplast :
        result = "leucoplast";
        break;
      case GENOME_proplastid :
        result = "proplastid";
        break;
      case GENOME_endogenous_virus :
        result = "endogenous virus";
        break;
      case GENOME_hydrogenosome :
        result = "hydrogenosome";
        break;
      case GENOME_chromosome :
        result = "chromosome";
        break;
      case GENOME_chromatophore :
        result = "chromatophore";
        break;
  }

  return result;
}

static void x_LowercasePlasmidOrElement (
  CharPtr def
)

{
  CharPtr  ptr;

  if (StringHasNoText (def)) return;

  def++;

  ptr = StringISearch (def, "plasmid");
  while (ptr != NULL) {
    if (*ptr == 'P') {
      *ptr = 'p';
    }
    ptr = StringISearch (ptr + 7, "plasmid");
  }

  ptr = StringISearch (def, "element");
  while (ptr != NULL) {
    if (*ptr == 'E') {
      *ptr = 'e';
    }
    ptr = StringISearch (ptr + 7, "element");
  }
}

static CharPtr x_TitleFromNC (
  DefLinePtr dlp
)

{
  Char        ch;
  CharPtr     completeseq = ", complete sequence";
  CharPtr     completegen = ", complete genome";
  Boolean     is_plasmid, has_plasmid = FALSE, virus_or_phage = FALSE;
  CharPtr     result = NULL, orgnl = NULL, pls_pfx = "";
  ValNodePtr  strings = NULL;

  if (dlp == NULL) return NULL;

  if (dlp->m_mi_biomol != MOLECULE_TYPE_GENOMIC &&
      dlp->m_mi_biomol != MOLECULE_TYPE_OTHER_GENETIC_MATERIAL) return NULL;

  if (StringHasNoText (dlp->m_taxname)) return NULL;

  if (StringISearch (dlp->m_taxname, "virus") != NULL ||
      StringISearch (dlp->m_taxname, "phage") != NULL) {
    virus_or_phage = TRUE;
  }

  if (StringDoesHaveText (dlp->m_plasmid)) {
    has_plasmid = TRUE;
  }

  orgnl = x_OrganelleName (dlp, has_plasmid, virus_or_phage, FALSE);

  is_plasmid = (Boolean) (dlp->m_genome == GENOME_plasmid);

  if (dlp->m_mi_completeness == 2 ||
      dlp->m_mi_completeness == 3 ||
      dlp->m_mi_completeness == 4 ||
      dlp->m_mi_completeness == 5) {
    /* remove "complete" component */
    completeseq = ", partial sequence";
    completegen = ", genome";
  }

  if (StringDoesHaveText (dlp->m_plasmid)) {
    if (StringISearch (dlp->m_plasmid, "plasmid") == NULL &&
        StringISearch (dlp->m_plasmid, "element") == NULL) {
      pls_pfx = "plasmid ";
    }
  }

  if (StringISearch (dlp->m_taxname, "plasmid") != NULL) {

    ValNodeCopyStr (&strings, 0, dlp->m_taxname);
    ValNodeCopyStr (&strings, 0, completeseq);

  } else if (is_plasmid) {

    if (StringDoesHaveText (dlp->m_plasmid)) {
      ValNodeCopyStr (&strings, 0, dlp->m_taxname);
      ValNodeCopyStr (&strings, 0, " ");
      ValNodeCopyStr (&strings, 0, pls_pfx);
      ValNodeCopyStr (&strings, 0, dlp->m_plasmid);
      ValNodeCopyStr (&strings, 0, completeseq);
    } else {
      ValNodeCopyStr (&strings, 0, dlp->m_taxname);
      ValNodeCopyStr (&strings, 0, " unnamed plasmid");
      ValNodeCopyStr (&strings, 0, completeseq);
    }

  } else if (StringDoesHaveText (dlp->m_plasmid)) {

    if (StringDoesHaveText (orgnl)) {
      ValNodeCopyStr (&strings, 0, dlp->m_taxname);
      ValNodeCopyStr (&strings, 0, " ");
      ValNodeCopyStr (&strings, 0, orgnl);
      ValNodeCopyStr (&strings, 0, " ");
      ValNodeCopyStr (&strings, 0, pls_pfx);
      ValNodeCopyStr (&strings, 0, dlp->m_plasmid);
      ValNodeCopyStr (&strings, 0, completeseq);
    } else {
      ValNodeCopyStr (&strings, 0, dlp->m_taxname);
      ValNodeCopyStr (&strings, 0, " ");
      ValNodeCopyStr (&strings, 0, pls_pfx);
      ValNodeCopyStr (&strings, 0, dlp->m_plasmid);
      ValNodeCopyStr (&strings, 0, completeseq);
    }

  } else if (StringDoesHaveText (orgnl)) {

    if (StringDoesHaveText (dlp->m_chromosome)) {
      ValNodeCopyStr (&strings, 0, dlp->m_taxname);
      ValNodeCopyStr (&strings, 0, " ");
      ValNodeCopyStr (&strings, 0, orgnl);
      ValNodeCopyStr (&strings, 0, " chromosome ");
      ValNodeCopyStr (&strings, 0, dlp->m_chromosome);
      ValNodeCopyStr (&strings, 0, completeseq);
    } else {
      ValNodeCopyStr (&strings, 0, dlp->m_taxname);
      ValNodeCopyStr (&strings, 0, " ");
      ValNodeCopyStr (&strings, 0, orgnl);
      ValNodeCopyStr (&strings, 0, completegen);
    }

  } else if (StringDoesHaveText (dlp->m_segment)) {

    if (StringStr (dlp->m_segment, "DNA") == NULL &&
        StringStr (dlp->m_segment, "RNA") == NULL &&
        StringStr (dlp->m_segment, "segment") == NULL &&
        StringStr (dlp->m_segment, "Segment") == NULL) {
      ValNodeCopyStr (&strings, 0, dlp->m_taxname);
      ValNodeCopyStr (&strings, 0, " segment ");
      ValNodeCopyStr (&strings, 0, dlp->m_segment);
      ValNodeCopyStr (&strings, 0, completegen);
    } else {
      ValNodeCopyStr (&strings, 0, dlp->m_taxname);
      ValNodeCopyStr (&strings, 0, " ");
      ValNodeCopyStr (&strings, 0, dlp->m_segment);
      ValNodeCopyStr (&strings, 0, completegen);
    }

  } else if (StringDoesHaveText (dlp->m_chromosome)) {

    ValNodeCopyStr (&strings, 0, dlp->m_taxname);
    ValNodeCopyStr (&strings, 0, " chromosome ");
    ValNodeCopyStr (&strings, 0, dlp->m_chromosome);
    ValNodeCopyStr (&strings, 0, completegen);

  } else {

    ValNodeCopyStr (&strings, 0, dlp->m_taxname);
    ValNodeCopyStr (&strings, 0, completegen);
  }

  result = x_CatenateValNodeStrings (strings);
  ValNodeFreeData (strings);
  if (result == NULL) return NULL;

  x_LowercasePlasmidOrElement (result);

  ch = result [0];
  if (IS_LOWER (ch)) {
    result [0] = TO_UPPER (ch);
  }

  return result;
}

typedef struct nmfeatdata {
  SeqFeatPtr  gene;
  SeqFeatPtr  cds;
  Int2        numgenes;
  Int2        numcds;
  Int2        numprots;
} NmFeatData, PNTR NmFeatPtr;

static void x_FindNMFeats (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  NmFeatPtr  nfp;

  if (sfp == NULL) return;
  nfp = (NmFeatPtr) userdata;
  if (nfp == NULL) return;

  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      nfp->gene = sfp;
      (nfp->numgenes)++;
      break;
    case SEQFEAT_CDREGION :
      nfp->cds = sfp;
      (nfp->numcds++);
      break;
    case SEQFEAT_PROT :
      (nfp->numprots)++;
      break;
    default :
      break;
  }
}

static Boolean x_IsFlyCG (
  CharPtr str
)
{
  Char  ch;

  if (StringHasNoText (str)) return FALSE;

  ch = *str;
  if (ch != 'C') return FALSE;

  str++;
  ch = *str;
  if (ch != 'G') return FALSE;

  str++;
  ch = *str;
  while (IS_DIGIT (ch)) {
    str++;
    ch = *str;
  }
  if (ch != '-') return FALSE;

  str++;
  ch = *str;
  if (ch != 'P') return FALSE;

  str++;
  ch = *str;
  if (IS_ALPHA (ch)) {
    str++;
    ch = *str;
    if (ch == '\0' || ch == ' ' || ch == ',' || ch == ';') return TRUE;
  }

  return FALSE;
}

static void x_FlyCG_PtoR (
  CharPtr str
)

{
  Char     ch;
  CharPtr  ptr;

  while (StringDoesHaveText (str)) {
    ch = *str;
    while (IS_WHITESP (ch)) {
      str++;
      ch = *str;
    }
    if (x_IsFlyCG (str)) {
      ptr = StringStr (str, "-P");
      if (ptr != NULL) {
        ptr [1] = 'R';
        return;
      }
    }
    while (ch != '\0' && (! IS_WHITESP (ch))) {
      str++;
      ch = *str;
    }
  }
}

static CharPtr x_TitleFromNM (
  DefLinePtr dlp
)

{
  Char         buf [512], buf2 [600];
  CharPtr      cds = NULL, gene = NULL, ptr, result = NULL;
  Uint2        entityID;
  size_t       len;
  NmFeatData   nfd;
  SeqEntryPtr  sep;

  if (dlp == NULL) return NULL;

  if (StringHasNoText (dlp->m_taxname)) return NULL;

  MemSet ((Pointer) &nfd, 0, sizeof (NmFeatData));

  entityID = ObjMgrGetEntityIDForPointer (dlp->m_bioseq);
  sep = GetBestTopParentForDataEx (entityID, dlp->m_bioseq, TRUE);

  VisitFeaturesInSep (sep, (Pointer) &nfd, x_FindNMFeats);
  if (nfd.numgenes != 1 || nfd.numcds != 1 || nfd.numprots < 1) return NULL;

  FeatDefLabel (nfd.gene, buf, sizeof (buf) - 1, OM_LABEL_CONTENT);
  gene = StringSaveNoNull (buf);

  FeatDefLabel (nfd.cds, buf, sizeof (buf) - 1, OM_LABEL_CONTENT);

  /* special case Drosophila RefSeq NM titles */
  if (StringICmp (dlp->m_taxname, "Drosophila melanogaster") == 0) {
    x_FlyCG_PtoR (buf);
  }
  ptr = StringStr (buf, "isoform ");
  if (ptr != NULL) {
    *ptr = '\0';
    ptr += 8;
    StringCpy (buf2, buf);
    StringCat (buf2, "transcript variant ");
    StringCat (buf2, ptr);
    cds = StringSaveNoNull (buf2);
  } else {
    cds = StringSaveNoNull (buf);
  }

  len = StringLen (dlp->m_taxname) + StringLen (cds) +
        StringLen (gene) + StringLen ("  (), mRNA") + 10;

  result = (CharPtr) MemNew (sizeof (Char) * len);

  if (result != NULL) {
    sprintf (result, "%s %s (%s), mRNA", dlp->m_taxname, cds, gene);
  }

  MemFree (gene);
  MemFree (cds);

  return result;
}

static CharPtr x_TitleFromNR (
  DefLinePtr dlp
)

{
  Char         buf [512];
  Uint2        entityID;
  CharPtr      gene = NULL,  rna = "miscRNA", result = NULL;
  size_t       len;
  NmFeatData   nfd;
  SeqEntryPtr  sep;

  if (dlp == NULL) return NULL;

  if (StringHasNoText (dlp->m_taxname)) return NULL;

  MemSet ((Pointer) &nfd, 0, sizeof (NmFeatData));

  entityID = ObjMgrGetEntityIDForPointer (dlp->m_bioseq);
  sep = GetBestTopParentForDataEx (entityID, dlp->m_bioseq, TRUE);

  VisitFeaturesInSep (sep, (Pointer) &nfd, x_FindNMFeats);
  if (nfd.numgenes < 1) return NULL;

  FeatDefLabel (nfd.gene, buf, sizeof (buf) - 1, OM_LABEL_CONTENT);
  gene = StringSaveNoNull (buf);

  switch (dlp->m_mi_biomol) {
      case MOLECULE_TYPE_PRE_MRNA :
        rna = "precursorRNA";
        break;
      case MOLECULE_TYPE_MRNA :
        rna = "mRNA";
        break;
      case MOLECULE_TYPE_RRNA :
        rna = "rRNA";
        break;
      case MOLECULE_TYPE_TRNA :
         rna = "tRNA";
        break;
      case MOLECULE_TYPE_SNRNA :
        rna = "snRNA";
        break;
      case MOLECULE_TYPE_SCRNA :
        rna = "scRNA";
        break;
      case MOLECULE_TYPE_CRNA :
        rna = "cRNA";
        break;
      case MOLECULE_TYPE_SNORNA :
        rna = "snoRNA";
        break;
      case MOLECULE_TYPE_TRANSCRIBED_RNA :
        rna = "miscRNA";
        break;
      case MOLECULE_TYPE_NCRNA :
        rna = "ncRNA";
        break;
      case MOLECULE_TYPE_TMRNA :
        rna = "tmRNA";
        break;
      default :
        break;
  }

  len = StringLen (dlp->m_taxname) + StringLen (gene) +
        StringLen (", ") + 30;

  result = (CharPtr) MemNew (sizeof (Char) * len);
  if (result != NULL) {
    sprintf (result, "%s %s, %s", dlp->m_taxname, gene, rna);
  }

  MemFree (gene);

  return result;
}

static CharPtr x_TitleFromPatent (
  DefLinePtr dlp
)

{
  Char  buf [80];

  if (dlp == NULL) return NULL;

  sprintf (buf, "Sequence %d from Patent %s %s",
           (int) dlp->m_patent_sequence,
           dlp->m_patent_country,
           dlp->m_patent_number);

  return StringSave (buf);
}

static CharPtr x_TitleFromPDB (
  DefLinePtr dlp
)

{
  Char        buf [40];
  Char        ch;
  CharPtr     result = NULL;
  ValNodePtr  strings = NULL;

  if (dlp == NULL) return NULL;

  ch = dlp->m_pdb_chain;
  if (IS_PRINT (ch)) {
    sprintf (buf, "Chain %c, ", ch);
    ValNodeCopyStr (&strings, 0, buf);
  }
  ValNodeCopyStr (&strings, 0, dlp->m_pdb_compound);

  result = x_CatenateValNodeStrings (strings);
  ValNodeFreeData (strings);

  return result;
}

typedef struct udxfeatdata {
  SeqIdPtr    bspid;
  Int4        longest;
  Uint1       processed;
  SeqFeatPtr  sfp;
} UdxFeatData, PNTR UdxFeatPtr;

static void x_GetLongestProtFeat (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  Int4        len;
  ProtRefPtr  prp;
  SeqIdPtr    sip;
  UdxFeatPtr  ufp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT) return;
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL) return;

  ufp = (UdxFeatPtr) userdata;
  if (ufp == NULL) return;

  sip = SeqLocId (sfp->location);
  if (sip == NULL) return;

  if (! SeqIdIn (sip, ufp->bspid)) return;
  len = SeqLocLen (sfp->location);
  if (len == -1) return;

  if (len > ufp->longest) {
    ufp->sfp = sfp;
    ufp->longest = len;
    ufp->processed = prp->processed;
  } else if (len == ufp->longest) {
    /* unprocessed 0 preferred over preprotein 1 preferred over mat peptide 2 */
    if (prp->processed < ufp->processed) {
      ufp->sfp = sfp;
      ufp->longest = len;
      ufp->processed = prp->processed;
    }
  }
}

static SeqFeatPtr x_GetLongestProteinUnindexed (
  BioseqPtr bsp
)

{
  BioseqSetPtr  bssp = NULL;
  UdxFeatData   ufd;

  if (bsp == NULL) return NULL;

  MemSet ((Pointer) &ufd, 0, sizeof (UdxFeatData));
  ufd.bspid = bsp->id;
  ufd.longest = 0;
  ufd.sfp = NULL;

  VisitFeaturesOnBsp (bsp, (Pointer) &ufd, x_GetLongestProtFeat);

  if (ufd.sfp != NULL && ufd.longest == bsp->length) return ufd.sfp;

  if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) bsp->idx.parentptr;
  }

  if (bssp != NULL && bssp->_class == BioseqseqSet_class_parts) {
    VisitFeaturesOnSet (bssp, (Pointer) &ufd, x_GetLongestProtFeat);

    if (bssp->idx.parenttype == OBJ_BIOSEQSET) {
      bssp = (BioseqSetPtr) bssp->idx.parentptr;
    }
  }

  if (bssp != NULL && bssp->_class == BioseqseqSet_class_segset) {
    VisitFeaturesOnSet (bssp, (Pointer) &ufd, x_GetLongestProtFeat);
  }

  return ufd.sfp;
}

static Boolean x_NotSpecialTaxName (
  CharPtr taxname
)

{
  if (StringHasNoText (taxname)) return TRUE;

  if (StringICmp (taxname, "synthetic construct") == 0) return FALSE;
  if (StringICmp (taxname, "artificial sequence") == 0) return FALSE;
  if (StringStr (taxname, "vector") != NULL) return FALSE;
  if (StringStr (taxname, "Vector") != NULL) return FALSE;

  return TRUE;
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

static CharPtr x_TitleFromProtein (
  DefLinePtr dlp
)

{
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  SeqFeatPtr         cds = NULL;
  Uint2              entityID;
  SeqMgrFeatContext  fcontext;
  GeneRefPtr         grp;
  Boolean            indexed;
  size_t             len;
  CharPtr            low_qual = "LOW QUALITY PROTEIN: ";
  CharPtr            organelle = NULL;
  OrgRefPtr          orp;
  Boolean            partial = FALSE;
  CharPtr            prefix = "";
  ProtRefPtr         prp;
  CharPtr            result = NULL;
  SeqFeatPtr         sfp = NULL;
  CharPtr            str;
  ValNodePtr         strings = NULL;
  CharPtr            taxname = NULL;
  CharPtr            title = NULL;
  CharPtr            tmp;
  ValNodePtr         vnp;

  if (dlp == NULL) return NULL;

  bsp = dlp->m_bioseq;
  if (bsp == NULL) return NULL;

  entityID = ObjMgrGetEntityIDForPointer (bsp);
  indexed = (Boolean) (SeqMgrFeaturesAreIndexed (entityID) != 0);

  if (indexed) {
    sfp = SeqMgrGetBestProteinFeature (bsp, NULL);
  } else {
    if (dlp->m_is_seg) {
      SeqMgrIndexFeatures (entityID, NULL);
      indexed = TRUE;
      sfp = SeqMgrGetBestProteinFeature (bsp, NULL);
    } else {
      sfp = x_GetLongestProteinUnindexed (bsp);
    }
  }

  if (dlp->m_mi_completeness > 1 && dlp->m_mi_completeness < 6) {
    partial = TRUE;
  }

  if (sfp != NULL) {
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    if (prp != NULL) {
      if (prp->name != NULL) {
        if (dlp->m_allprotnames) {
          for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
            str = (CharPtr) vnp->data.ptrvalue;
            ValNodeCopyStr (&strings, 0, prefix);
            ValNodeCopyStr (&strings, 0, str);
            prefix = "; ";
          }
          title = x_CatenateValNodeStrings (strings);
          strings = ValNodeFreeData (strings);
        } else {
          vnp = prp->name;
          /* although vnp should not be NULL, a compiler/optimizer bug might let it, so check again */
          if (vnp != NULL && vnp->data.ptrvalue != NULL) {
            str = (CharPtr) vnp->data.ptrvalue;
            title = StringSave (str);
          }
        }
        x_TrimPunctuationFromEnd (title);

        /* if hypothetical protein, append locus_tag */
        if (StringICmp (title, "hypothetical protein") == 0) {
          if (! indexed) {
            SeqMgrIndexFeatures (entityID, NULL);
            indexed = TRUE;
          }
          if (cds == NULL) {
            cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
          }
          if (cds != NULL) {
            grp = SeqMgrGetGeneXref (cds);
            if (grp == NULL) {
              sfp = SeqMgrGetOverlappingFeature (cds->location, FEATDEF_GENE, NULL, 0, NULL, LOCATION_SUBSET, NULL);
              if (sfp != NULL) {
                grp = (GeneRefPtr) sfp->data.value.ptrvalue;
              }
            }
            if (grp != NULL) {
              if (grp->locus_tag != NULL) {
                len = StringLen (title) + StringLen (grp->locus_tag) + 20;
                str = (CharPtr) MemNew (sizeof (Char) * len);
                if (str != NULL) {
                  StringCat (str, title);
                  StringCat (str, " ");
                  StringCat (str, grp->locus_tag);
                  MemFree (title);
                  title = str;
                }
              }
            }
          }
        }
      }

      if ( title == NULL && prp->desc != NULL) {
        title = StringSave (prp->desc);
      }

      if ( title == NULL && prp->activity != NULL) {
        vnp = prp->activity;
        str = (CharPtr) vnp->data.ptrvalue;
        title = StringSave (str);
      }
    }
  }

  if (title == NULL) {
    if (! indexed) {
      SeqMgrIndexFeatures (entityID, NULL);
      indexed = TRUE;
    }
    if (cds == NULL) {
      cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
    }
    if (cds != NULL) {
      grp = SeqMgrGetGeneXref (cds);
      if (grp == NULL) {
        sfp = SeqMgrGetOverlappingFeature (cds->location, FEATDEF_GENE, NULL, 0, NULL, LOCATION_SUBSET, NULL);
        if (sfp != NULL) {
          grp = (GeneRefPtr) sfp->data.value.ptrvalue;
        }
      }
      if (grp != NULL) {
        str = NULL;
        if (grp->locus != NULL) {
          str = grp->locus;
        } else if (grp->syn != NULL) {
          vnp = grp->syn;
          str = (CharPtr) vnp->data.ptrvalue;
        } else if (grp->desc != NULL) {
          str = grp->desc;
        }
        if (StringDoesHaveText (str)) {
          ValNodeCopyStr (&strings, 0, str);
          ValNodeCopyStr (&strings, 0, " gene product");
          title = x_CatenateValNodeStrings (strings);
          strings = ValNodeFreeData (strings);
        }
      }
    }
  }

  if (title == NULL) {
    title = StringSave ("unnamed protein product");
  }

  if (title != NULL) {
    x_TrimPunctuationFromEnd (title);
  }

  taxname = dlp->m_taxname;
  if (StringHasNoText (taxname) || x_NotSpecialTaxName (taxname)) {
    if (! indexed) {
      SeqMgrIndexFeatures (entityID, NULL);
      indexed = TRUE;
    }
    if (cds == NULL) {
      cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
    }
    if (cds != NULL) {
      sfp = SeqMgrGetOverlappingSource (cds->location, &fcontext);
      if (sfp != NULL) {
        biop = (BioSourcePtr) sfp->data.value.ptrvalue;
        if (biop != NULL) {
          orp = biop->org;
          if (orp != NULL) {
            taxname = orp->taxname;
          }
        }
      }
    }
  }

  if (dlp->m_genome >= GENOME_chloroplast && dlp->m_genome <= GENOME_chromatophore) {
    organelle = proteinOrganellePrefix [dlp->m_genome];
    if (StringNICmp (organelle, taxname, StringLen (organelle)) == 0) {
      organelle = NULL;
    }
  }

  if (cds == NULL) {
    if (! indexed) {
      SeqMgrIndexFeatures (entityID, NULL);
      indexed = TRUE;
    }
    cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
  }
  if (cds != NULL) {
    if (x_CDShasLowQualityException (dlp, cds)) {
      if (StringStr (title, low_qual) == NULL) {
        len = StringLen (title) + StringLen (low_qual) + 6;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          StringCat (tmp, low_qual);
          StringCat (tmp, title);
          MemFree (title);
          title = tmp;
        }
      }
    }
  }

  if (partial) {
    if (StringStr (title, ", partial") == NULL) {
      len = StringLen (title) + 12;
      tmp = (CharPtr) MemNew (sizeof (Char) * len);
      if (tmp != NULL) {
        StringCat (tmp, title);
        StringCat (tmp, ", partial");
        MemFree (title);
        title = tmp;
      }
    }
  }
  if (StringDoesHaveText (organelle)) {
    if (StringStr (title, organelle) == NULL) {
      len = StringLen (title) + StringLen (organelle) + 6;
      tmp = (CharPtr) MemNew (sizeof (Char) * len);
      if (tmp != NULL) {
        StringCat (tmp, title);
        StringCat (tmp, " (");
        StringCat (tmp, organelle);
        StringCat (tmp, ")");
        MemFree (title);
        title = tmp;
      }
    }
  }

  if (StringDoesHaveText (taxname)) {
    if (StringStr (title, taxname) == NULL) {
      len = StringLen (title) + StringLen (taxname) + 6;
      tmp = (CharPtr) MemNew (sizeof (Char) * len);
      if (tmp != NULL) {
        StringCat (tmp, title);
        StringCat (tmp, " [");
        StringCat (tmp, taxname);
        StringCat (tmp, "]");
        MemFree (title);
        title = tmp;
      }
    }
  }

  if (result == NULL) {
    result = StringSave (title);
  }

  MemFree (title);

  return result;
}

static CharPtr x_TitleFromSegSeq (
  DefLinePtr dlp
)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  ccontext;
  SeqFeatPtr         cds;
  CharPtr            cln = NULL;
  CharPtr            complete = "gene, complete cds";
  Uint2              entityID;
  SeqMgrFeatContext  gcontext;
  SeqFeatPtr         gene;
  GeneRefPtr         grp;
  CharPtr            label = NULL;
  size_t             len;
  CharPtr            locus = NULL;
  CharPtr            modifier = NULL;
  CharPtr            product = NULL;
  CharPtr            result = NULL;
  CharPtr            str;
  CharPtr            taxname = NULL;
  ValNodePtr         vnp;

  if (dlp == NULL) return NULL;

  bsp = dlp->m_bioseq;
  if (bsp == NULL) return NULL;

  entityID = ObjMgrGetEntityIDForPointer (bsp);
  if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
    SeqMgrIndexFeatures (entityID, NULL);
  }

  cds = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &ccontext);

  if (cds != NULL) {
    if (cds->partial) {
      complete = "gene, partial cds";
    }
    product = ccontext.label;
    grp = SeqMgrGetGeneXref (cds);
    if (grp != NULL) {
      if (StringDoesHaveText (grp->locus)) {
        locus = grp->locus;
      } else {
        vnp = grp->syn;
        if (vnp != NULL) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringDoesHaveText (str)) {
            locus = str;
          }
        }
      }
    }
    if (locus == NULL) {
      gene = SeqMgrGetOverlappingGene (cds->location, &gcontext);
      if (gene != NULL) {
        locus = gcontext.label;
      }
    }
  } else {
    if (StringDoesHaveText (dlp->m_strain) && (! x_EndsWithStrain (dlp))) {
      modifier = dlp->m_strain;
      label = " strain ";
    } else if (StringDoesHaveText (dlp->m_clone)) {
      cln = x_DescribeClones (dlp);
      modifier = cln;
    } else if (StringDoesHaveText (dlp->m_isolate)) {
      modifier = dlp->m_isolate;
      label = " isolate ";
    }
  }

  taxname = dlp->m_taxname;
  if (StringHasNoText (taxname)) {
    taxname = "Unknown";
  }

  len = StringLen (taxname) + StringLen (label) + StringLen (modifier) +
        StringLen (product) + StringLen (locus) + StringLen (complete) + 10;

  result = (CharPtr) MemNew (sizeof (Char) * len);
  if (result == NULL) {
    MemFree (cln); 
    return NULL;
  }

  if (taxname != NULL) {
    StringCat (result, taxname);
  }

  if (modifier != NULL) {
    if (label != NULL) {
      StringCat (result, label);
    }
    StringCat (result, modifier);
  }

  if (product != NULL) {
    StringCat (result, " ");
    StringCat (result, product);
  }
  if (locus != NULL) {
    StringCat (result, " (");
    StringCat (result, locus);
    StringCat (result, ")");
  }
  if (product != NULL || locus != NULL) {
    StringCat (result, " ");
    StringCat (result, complete);
  }
  TrimSpacesAroundString (result);

  MemFree (cln); 

  return result;
}

static CharPtr x_TitleFromWGS (
  DefLinePtr dlp
)

{
  Char        ch;
  CharPtr     result = NULL, cln, stn, ptr;
  ValNodePtr  strings = NULL;

  if (dlp == NULL) return NULL;

  ValNodeCopyStr (&strings, 0, dlp->m_taxname);

  if (StringDoesHaveText (dlp->m_strain)) {
    if (! x_EndsWithStrain (dlp)) {
      ValNodeCopyStr (&strings, 0, " strain ");
      stn = StringSave (dlp->m_strain);
      ptr = StringChr (stn, ';');
      if (ptr != NULL) {
        *ptr = '\0';
      }
      ValNodeCopyStr (&strings, 0, stn);
      MemFree (stn);
    }
  }

  if (StringDoesHaveText (dlp->m_chromosome)) {
    ValNodeCopyStr (&strings, 0, " chromosome ");
    ValNodeCopyStr (&strings, 0, dlp->m_chromosome);
  }

  cln = x_DescribeClones (dlp);
  if (StringDoesHaveText (cln)) {
    ValNodeCopyStr (&strings, 0, cln);
  }
  MemFree (cln); 

  if (StringDoesHaveText (dlp->m_map)) {
    ValNodeCopyStr (&strings, 0, " map ");
    ValNodeCopyStr (&strings, 0, dlp->m_map);
  }

  if (StringDoesHaveText (dlp->m_plasmid)) {
    if (dlp->m_is_wgs) {
      ValNodeCopyStr (&strings, 0, " plasmid ");
      ValNodeCopyStr (&strings, 0, dlp->m_plasmid);
    }
  }

  if (StringDoesHaveText (dlp->m_general_str)) {
    ValNodeCopyStr (&strings, 0, " ");
    ValNodeCopyStr (&strings, 0, dlp->m_general_str);
  }

  result = x_CatenateValNodeStrings (strings);
  ValNodeFreeData (strings);
  if (result == NULL) return NULL;

  ch = result [0];
  if (IS_LOWER (ch)) {
    result [0] = TO_UPPER (ch);
  }

  return result;
}

static CharPtr x_SetPrefix (
  DefLinePtr dlp,
  CharPtr title
)

{
  CharPtr  prefix = NULL;

  if (dlp == NULL) return NULL;

  if (dlp->m_is_unverified) {
    if (StringStr (title, "UNVERIFIED") == NULL) {
      prefix = "UNVERIFIED: ";
    }
  } else if (dlp->m_is_tsa) {
    prefix = "TSA: ";
  } else if (dlp->m_third_party) {
    if (dlp->m_tpa_exp) {
      prefix = "TPA_exp: ";
    } else if (dlp->m_tpa_inf) {
      prefix = "TPA_inf: ";
    } else if (dlp->m_tpa_reasm) {
      prefix = "TPA_reasm: ";
    } else {
      prefix = "TPA: ";
    }
  }

  return StringSave (prefix);
}

static CharPtr x_SetSuffix (
  DefLinePtr dlp,
  CharPtr title
)

{
  Char     buf [80];
  size_t   len;
  Int4     num_segs, num_gaps;
  CharPtr  orgnl = NULL, str, suffix = "", un = "";

  if (dlp == NULL) return NULL;

  switch (dlp->m_mi_tech) {
      case MI_TECH_htgs_0 :
        if (StringStr (title, "LOW-PASS") == NULL) {
          suffix = ", LOW-PASS SEQUENCE SAMPLING";
        }
        break;
      case MI_TECH_htgs_1 :
        un = "un";
        /* fall through */
      case MI_TECH_htgs_2 :
        if (dlp->m_htgs_draft) {
          if (StringStr (title, "WORKING DRAFT") == NULL) {
            suffix = ", WORKING DRAFT SEQUENCE";
          }
        } else if (! dlp->m_htgs_cancelled) {
          if (StringStr (title, "SEQUENCING IN") == NULL) {
            suffix = ", *** SEQUENCING IN PROGRESS ***";
          }
        }
        if (dlp->m_is_delta) {
          if (CountGapsInDeltaSeq (dlp->m_bioseq, &num_segs, &num_gaps, NULL, NULL, NULL, 0)) {
            if (num_gaps > 0) {
              sprintf (buf, "%s, %ld %sordered pieces", suffix, (long) (num_gaps + 1), un);
              suffix = StringSave (buf);
              return suffix;
            }
          }
        }
        break;
      case MI_TECH_htgs_3 :
        if (StringStr (title, "complete sequence") == NULL) {
          suffix = ", complete sequence";
        }
        break;
      case MI_TECH_est :
        if (StringStr (title, "mRNA sequence") == NULL) {
          suffix = ", mRNA sequence";
        }
        break;
      case MI_TECH_sts :
        if (StringStr (title, "sequence tagged site") == NULL) {
          suffix = ", sequence tagged site";
        }
        break;
      case MI_TECH_survey :
        if (StringStr (title, "genomic survey sequence") == NULL) {
          suffix = ", genomic survey sequence";
        }
        break;
      case MI_TECH_wgs :
        if (dlp->m_wgs_master) {
          if (StringStr (title, "whole genome shotgun sequencing project") == NULL) {
            suffix = ", whole genome shotgun sequencing project";
          }
        } else {
          if (StringStr (title, "whole genome shotgun sequence") == NULL) {
            orgnl = x_OrganelleName (dlp, FALSE, FALSE, TRUE);
            len = StringLen (", whole genome shotgun sequence") + StringLen (orgnl) + 5;
            str = (CharPtr) MemNew (sizeof (Char) * len);
            if (str != NULL) {
              if (StringDoesHaveText (orgnl)) {
                StringCat (str, " ");
                StringCat (str, orgnl);
              }
              StringCat (str, ", whole genome shotgun sequence");
              return str;
            }
          }
        }
        break;
      default :
        break;
  }

  return StringSave (suffix);
}

NLM_EXTERN CharPtr NewCreateDefLine (
  ItemInfoPtr iip,
  BioseqPtr bsp,
  Boolean ignoreTitle,
  Boolean extProtTitle
)

{
  DefLinePtr     dlp;
  Uint2          entityID;
  size_t         len;
  ObjValNodePtr  ovp;
  CharPtr        result = NULL, prefix = NULL, suffix = NULL, title = NULL;
  SeqDescrPtr    sdp = NULL;
  CharPtr        str = NULL;

  if (bsp == NULL) return NULL;

  /* now using GetNextDescriptorUnindexed, so need to have called AssignIDsInEntityEx */
  if (bsp->idx.entityID == 0) {
    entityID = ObjMgrGetEntityIDForPointer (bsp);
    if (entityID != 0) {
      AssignIDsInEntityEx (entityID, 0, NULL, NULL);
    }
  }

  dlp = (DefLinePtr) MemNew (sizeof (DefLineData));
  if (dlp == NULL) return NULL;

  dlp->m_low_quality_fsa = TextFsaNew ();
  TextFsaAdd (dlp->m_low_quality_fsa, "heterogeneous population sequenced");
  TextFsaAdd (dlp->m_low_quality_fsa, "low-quality sequence region");
  TextFsaAdd (dlp->m_low_quality_fsa, "unextendable partial coding region");

  /* set flags from record components */
  dlp->m_iip = iip;
  dlp->m_bioseq = bsp;

  dlp->m_reconstruct = ignoreTitle;
  dlp->m_allprotnames = extProtTitle;

  /* clear ItemInfo fields */
  if (iip != NULL) {
    iip->entityID = 0;
    iip->itemID = 0;
    iip->itemtype = 0;
  }

  /* set flags from record components */
  x_SetFlags (dlp);

  if (! dlp->m_reconstruct) {
    /* look for existing instantiated title */
    if (dlp->m_is_aa && (! dlp->m_is_pdb)) {
      sdp = BioseqGetSeqDescr (bsp, Seq_descr_title, NULL);
      if (sdp != NULL && sdp->choice == Seq_descr_title) {
        str = (CharPtr) sdp->data.ptrvalue;
      }
    } else {
      sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_title, NULL);
      if (sdp != NULL && sdp->choice == Seq_descr_title) {
        str = (CharPtr) sdp->data.ptrvalue;
      }
    }
    if (StringDoesHaveText (str)) {
      title = StringSave (str);
      /* strip trailing periods, commas, semicolons, etc. */
      x_TrimPunctuationFromEnd (title);

      /* set ItemInfo fields for selection */
      if (iip != NULL && sdp != NULL && sdp->extended != 0) {
        ovp = (ObjValNodePtr) sdp;
        iip->entityID = ovp->idx.entityID;
        iip->itemtype = ovp->idx.itemtype;
        iip->itemID = ovp->idx.itemID;
      }
    }
  }

  /* use appropriate algorithm if title needs to be generated */
  if (StringHasNoText (title)) {
    /* PDB and patent records do not normally need source data */
    if (dlp->m_is_pdb) {
      title = x_TitleFromPDB (dlp);
    } else if (dlp->m_is_patent) {
      title = x_TitleFromPatent (dlp);
    }

    if (StringHasNoText (title)) {
      /* set fields from source information */
      x_SetBioSrc (dlp);

      /* several record types have specific methods */
      if (dlp->m_is_nc) {
        title = x_TitleFromNC (dlp);
      } else if (dlp->m_is_nm) {
        title = x_TitleFromNM (dlp);
      } else if (dlp->m_is_nr) {
        title = x_TitleFromNR (dlp);
      } else if (dlp->m_is_aa) {
        title = x_TitleFromProtein (dlp);
      } else if (dlp->m_is_seg && (! dlp->m_is_est_sts_gss)) {
        title = x_TitleFromSegSeq (dlp);
      } else if (dlp->m_is_tsa || (dlp->m_is_wgs && (! dlp->m_wgs_master))) {
        title = x_TitleFromWGS (dlp);
      }
    }

    if (StringHasNoText (title)) {
      /* default title using source fields */
      title = x_TitleFromBioSrc (dlp);
    }

    if (StringHasNoText (title)) {
      /* last resort title created here */
      title = StringSave ("No definition line found");
    }
  }

  /* remove TPA or TSA prefix, will rely on other data in record to set */
  if (StringNICmp (title, "TPA:", 4) == 0) {
    x_TrimFirstNCharacters (title, 4);
  } else if (StringNICmp (title, "TPA_exp:", 8) == 0) {
    x_TrimFirstNCharacters (title, 8);
  } else if (StringNICmp (title, "TPA_inf:", 8) == 0) {
    x_TrimFirstNCharacters (title, 8);
  } else if (StringNICmp (title, "TPA_reasm:", 10) == 0) {
    x_TrimFirstNCharacters (title, 10);
  } else if (StringNICmp (title, "TSA:", 4) == 0) {
    x_TrimFirstNCharacters (title, 4);
  } else if (StringNICmp (title, "UNVERIFIED:", 11) == 0) {
    x_TrimFirstNCharacters (title, 11);
  }

  /* strip leading spaces remaining after removal of old TPA or TSA prefixes */
  TrimSpacesAroundString (title);

  /* strip trailing commas, semicolons, and spaces (period may be an sp. species) */
  x_TrimMostPunctFromEnd (title);

  /* calcualte prefix */
  prefix = x_SetPrefix (dlp, title);

  /* calculate suffix */
  suffix = x_SetSuffix (dlp, title);

  len = StringLen (prefix) + StringLen (title) + StringLen (suffix) + 4;
  result = (CharPtr) MemNew (sizeof (Char) * len);

  if (result != NULL) {
    StringCat (result, prefix);
    StringCat (result, title);
    StringCat (result, suffix);
  }

  MemFree (prefix);
  MemFree (title);
  MemFree (suffix);

  TextFsaFree (dlp->m_low_quality_fsa);

  dlp = MemFree (dlp);

  Asn2gnbkCompressSpaces (result);

  return result;
}

NLM_EXTERN Boolean NewCreateDefLineBuf (
  ItemInfoPtr iip,
  BioseqPtr bsp,
  CharPtr buf,
  Uint4 buflen,
  Boolean ignoreTitle,
  Boolean extProtTitle
)

{
  CharPtr  title = NULL;

  if (bsp == NULL || buf == NULL|| buflen == 0) return FALSE;

  title = NewCreateDefLine (iip, bsp, ignoreTitle, extProtTitle);
  StringNCpy_0 (buf, title, buflen);
  MemFree (title);

  return TRUE;
}

